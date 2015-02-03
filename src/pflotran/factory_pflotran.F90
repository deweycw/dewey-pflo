module Factory_PFLOTRAN_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: PFLOTRANInitializePrePetsc, &
            PFLOTRANInitializePostPetsc, &
            PFLOTRANFinalize

contains

! ************************************************************************** !

subroutine PFLOTRANInitializePrePetsc(multisimulation,option)
!
! Sets up PFLOTRAN subsurface simulation framework prior to PETSc 
!   initialization
! Author: Glenn Hammond
! Date: 06/07/13
!
  use Option_module
  use Input_Aux_module
  use Multi_Simulation_module
  
  implicit none
  
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: bool_flag
  PetscBool :: option_found
  
  ! NOTE: Cannot add anything that requires PETSc in this routine as PETSc 
  !       has not yet been initialized.
  
  call PFLOTRANInitCommandLineSettings(option)
  ! initialize stochastic realizations here
  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    multisimulation => MultiSimulationCreate()
    call MultiSimulationInitialize(multisimulation,option)
  endif
  
end subroutine PFLOTRANInitializePrePetsc

! ************************************************************************** !

subroutine PFLOTRANInitializePostPetsc(simulation,multisimulation,option)
!
! Sets up PFLOTRAN subsurface simulation framework after PETSc initialization
! Author: Glenn Hammond
! Date: 06/17/13
!
  use Option_module
  use Multi_Simulation_module
  use Simulation_Base_class
  use Logging_module
  use EOS_module
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type), pointer :: option
  
  character(len=MAXSTRINGLENGTH) :: filename
  PetscErrorCode :: ierr

  ! must come after logging is created
  call LoggingSetupComplete()
  call MultiSimulationIncrement(multisimulation,option)
  call OptionBeginTiming(option)

  ! popped in SimulationBaseInitializeRun()
  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_init,ierr);CHKERRQ(ierr)
  
  call EOSInit()
  filename = trim(option%global_prefix) // trim(option%group_prefix) // &
             '.out'
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename, action="write", status="unknown")
  endif
  
  call PFLOTRANReadSimulation(simulation,option)
  
end subroutine PFLOTRANInitializePostPetsc

! ************************************************************************** !

subroutine PFLOTRANReadSimulation(simulation,option)
!
! Sets up PFLOTRAN subsurface simulation framework after PETSc initialization
! Author: Glenn Hammond
! Date: 06/17/13
!
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Simulation_Base_class
  use PM_Base_class
  use PM_Surface_TH_class
  use PMC_Base_class
  
  use Factory_Subsurface_module
  use Factory_Hydrogeophysics_module
  use Factory_Surf_Subsurf_module
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation
  type(option_type), pointer :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: simulation_type
  
  class(pm_base_type), pointer :: pm_master
  class(pm_base_type), pointer :: cur_pm
  class(pm_base_type), pointer :: new_pm

  class(pmc_base_type), pointer :: pmc_master
  
  nullify(pm_master)
  nullify(cur_pm)
  nullify(new_pm)
  
  nullify(pmc_master)
  
  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = 'SIMULATION'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  word = ''
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'PROCESS_MODEL','SIMULATION')
    
    call StringToUpper(word)
    select case(trim(word))
      case('SIMULATION_TYPE')
          call InputReadWord(input,option,simulation_type,PETSC_TRUE)
          call InputErrorMsg(input,option,'simulation_type', &
                             'SIMULATION')
      case('PROCESS_MODELS')
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'process_model', &
                             'SIMULATION,PROCESS_MODELS')
          call InputReadWord(input,option,name,PETSC_TRUE)
          call InputErrorMsg(input,option,'name','SIMULATION,PROCESS_MODEL')
          call StringToUpper(word)
          select case(trim(word))
            case('SUBSURFACE_FLOW')
              call SubsurfaceReadFlowPM(input, option, new_pm)
            case('SUBSURFACE_TRANSPORT')
              call SubsurfaceReadRTPM(input, option, new_pm)
            case('WASTE_FORM')
              call SubsurfaceReadWasteFormPM(input, option,new_pm)
            case('HYDROGEOPHYSICS')
            case('SURFACE_SUBSURFACE')
              option%surf_flow_on = PETSC_TRUE
              new_pm => PMSurfaceTHCreate()
            case('GEOMECHANICS')
            case default
              option%io_buffer =  'PROCESS_MODEL " ' // trim(word) // &
                '" not recognized.'
              call printErrMsg(option)
          end select
          new_pm%name = name
          new_pm%option => option
          if (associated(cur_pm)) then
            cur_pm%next => new_pm
          else
            cur_pm => new_pm
          endif
          if (.not.associated(pm_master)) then
            pm_master => new_pm
          endif
          cur_pm => new_pm
          nullify(new_pm)
        enddo
      case('MASTER')
        call PFLOTRANSetupPMCHierarchy(input,option,pmc_master)
    end select
  enddo

  select case(simulation_type)
    case('CUSTOM')
      ! link process models with their respective couplers
      cur_pm => pm_master
      do
        if (.not.associated(cur_pm)) exit
        call PFLOTRANLinkPMToPMC(input,option,pmc_master,cur_pm)
        cur_pm => cur_pm%next
      enddo
      call pmc_master%CheckNullPM(option)
    case('SUBSURFACE')
      call SubsurfaceInitialize(simulation,pm_master,option)  
    case('HYDROGEOPHYSICS')
      call HydrogeophysicsInitialize(simulation,pm_master,option)
    case('SURFACE_SUBSURFACE')
      call SurfSubsurfaceInitialize(simulation,pm_master,option)
    case default
      option%io_buffer =  'SIMULATION_TYPE " ' // trim(word) // &
        '" not recognized.'
      call printErrMsg(option)
  end select
  
  call InputDestroy(input)
  
end subroutine PFLOTRANReadSimulation

! ************************************************************************** !

recursive subroutine PFLOTRANSetupPMCHierarchy(input,option,pmc)
!
! Forms a linked list of named dummy pmcs as placeholders
! Author: Glenn Hammond
! Date: 12/10/14
!
  use Option_module
  use Input_Aux_module
  use PMC_Base_class
  use String_module
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  class(pmc_base_type), pointer :: pmc
  
  character(len=MAXWORDLENGTH) :: word
  
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'PMC name','SIMULATION')
    ! at this point, we are creating a 
  pmc => PMCBaseCreate()
  pmc%name = word

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'CHILD or PEER','SIMULATION')
    call StringToUpper(word)
    select case(trim(word))
      case('PEER')
        call PFLOTRANSetupPMCHierarchy(input,option,pmc%peer)
      case('CHILD')
        call PFLOTRANSetupPMCHierarchy(input,option,pmc%child)
      case default
        call InputKeywordUnrecognized(word,'PFLOTRANSetupPMCHierarchy',option)
    end select    
  enddo
  
end subroutine PFLOTRANSetupPMCHierarchy

! ************************************************************************** !

recursive subroutine PFLOTRANLinkPMToPMC(input,option,pmc,pm)
!
! Forms a linked list of named dummy pmcs as placeholders
! Author: Glenn Hammond
! Date: 12/10/14
!
  use Option_module
  use Input_Aux_module
  use String_module
  use PM_Base_class
  use PMC_Base_class
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  class(pmc_base_type), pointer :: pmc
  class(pm_base_type), pointer :: pm

  if (.not.associated(pmc)) return
  
  print *, pmc%name, pm%name
  if (StringCompareIgnoreCase(pmc%name,pm%name)) then
    pmc%pms => pm
    return
  endif
  
  call PFLOTRANLinkPMToPMC(input,option,pmc%peer,pm)
  call PFLOTRANLinkPMToPMC(input,option,pmc%child,pm)
  
end subroutine PFLOTRANLinkPMToPMC

! ************************************************************************** !

subroutine PFLOTRANFinalize(option)
!
! Destroys PFLOTRAN subsurface simulation framework
! Author: Glenn Hammond
! Date: 06/07/13
!
  use Option_module
  use Logging_module
  
  implicit none
  
  type(option_type) :: option
  PetscErrorCode :: ierr
  
  ! pushed in FinalizeRun()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)
  call OptionEndTiming(option)
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif

end subroutine PFLOTRANFinalize

! ************************************************************************** !

subroutine PFLOTRANInitCommandLineSettings(option)
  ! 
  ! Initializes PFLOTRAN output filenames, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscBool :: option_found
  PetscBool :: bool_flag
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  PetscInt :: i
  PetscErrorCode :: ierr
  
  ! check for non-default input filename
  option%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename, &
                                 pflotranin_option_found,option)
  string = '-input_prefix'
  call InputGetCommandLineString(string,option%input_prefix, &
                                 input_prefix_option_found,option)
  
  if (pflotranin_option_found .and. input_prefix_option_found) then
    option%io_buffer = 'Cannot specify both "-pflotranin" and ' // &
      '"-input_prefix" on the command lines.'
    call printErrMsg(option)
  else if (pflotranin_option_found) then
    strings => StringSplit(option%input_filename,'.')
    option%input_prefix = strings(1)
    deallocate(strings)
    nullify(strings)
  else if (input_prefix_option_found) then
    option%input_filename = trim(option%input_prefix) // '.in'
  endif
  
  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)
  if (.not.option_found) option%global_prefix = option%input_prefix  
  
  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-v'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%verbosity = i
 
  ! this will get overwritten later if stochastic
  string = '-realization_id'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) then
    if (i < 1) then
      option%io_buffer = 'realization_id must be greater than zero.'
      call printErrMsg(option)
    endif
    option%id = i
  endif
  
  ! this will get overwritten later if stochastic
  string = '-simulation_mode'
  call InputGetCommandLineString(string,string2, &
                                 option_found,option)
  if (option_found) then
    call StringToUpper(string2)
    option%simulation_mode = string2
  endif

end subroutine PFLOTRANInitCommandLineSettings

end module Factory_PFLOTRAN_module