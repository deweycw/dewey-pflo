module Reaction_Sandbox_JinBethke_O2aq_DOC_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_jinbethke_o2aq_doc_type
    PetscInt :: auxiliary_offset
    PetscInt :: bicarbonate_id
    PetscInt :: o2aq_id
    PetscInt :: doc_id
    PetscInt :: h_ion_id
    PetscInt :: xim_id
    PetscReal :: rmax
    PetscReal :: Kdonor
    PetscReal :: Kacceptor
    PetscReal :: Y
    PetscReal :: m
    PetscReal :: chi
    PetscReal :: o2_threshold
    PetscReal :: dg0
    

  contains
    procedure, public :: ReadInput => JinBethkeO2aqDOCReadInput
    procedure, public :: Setup => JinBethkeO2aqDOCSetup
    procedure, public :: AuxiliaryPlotVariables => O2aqAuxiliaryPlotVariables
    procedure, public :: Evaluate => JinBethkeO2aqDOCEvaluate
  end type reaction_sandbox_jinbethke_o2aq_doc_type

  public :: JinBethkeO2aqDOCCreate, &
            JinBethkeO2aqDOCSetup
contains
! ************************************************************************** !
function JinBethkeO2aqDOCCreate()
  !
  ! Allocates Ferrihydrite reaction object.
  !
  implicit none
  class(reaction_sandbox_jinbethke_o2aq_doc_type), pointer :: JinBethkeO2aqDOCCreate
  allocate(JinBethkeO2aqDOCCreate)
  JinBethkeO2aqDOCCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  JinBethkeO2aqDOCCreate%doc_id = UNINITIALIZED_INTEGER
  JinBethkeO2aqDOCCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  JinBethkeO2aqDOCCreate%o2aq_id = UNINITIALIZED_INTEGER
  JinBethkeO2aqDOCCreate%h_ion_id = UNINITIALIZED_INTEGER
  JinBethkeO2aqDOCCreate%xim_id = UNINITIALIZED_INTEGER

  JinBethkeO2aqDOCCreate%rmax = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%Kdonor = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%Kacceptor = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%Y = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%m = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%chi = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%o2_threshold = UNINITIALIZED_DOUBLE
  JinBethkeO2aqDOCCreate%dg0 = UNINITIALIZED_DOUBLE

  nullify(JinBethkeO2aqDOCCreate%next)
end function JinBethkeO2aqDOCCreate
! ************************************************************************** !
subroutine JinBethkeO2aqDOCReadInput(this,input,option)
  !
  ! Reads calcite reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_jinbethke_o2aq_doc_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,JINBETHKE_O2AQ_DOC'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('RMAX')
        call InputReadDouble(input,option,this%rmax)
        call InputErrorMsg(input,option,word,error_string)
      case('DG0')
        call InputReadDouble(input,option,this%dg0)
        call InputErrorMsg(input,option,word,error_string)
      case('K_DONOR')
        call InputReadDouble(input,option,this%Kdonor)
        call InputErrorMsg(input,option,word,error_string)
      case('K_ACCEPTOR')
        call InputReadDouble(input,option,this%Kacceptor)
        call InputErrorMsg(input,option,word,error_string)
      case('Y')
        call InputReadDouble(input,option,this%Y)
        call InputErrorMsg(input,option,word,error_string)
      case('M')
        call InputReadDouble(input,option,this%m)
        call InputErrorMsg(input,option,word,error_string)
      case('CHI')
        call InputReadDouble(input,option,this%chi)
        call InputErrorMsg(input,option,word,error_string)
      case('O2_THRESHOLD')
        call InputReadDouble(input,option,this%o2_threshold)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rmax) .or. &
      Uninitialized(this%Kdonor) .or. &
      Uninitialized(this%Kacceptor) .or. &
      Uninitialized(this%Y) .or. &
      Uninitialized(this%m) .or. &
      Uninitialized(this%o2_threshold) .or. &
      Uninitialized(this%dg0) .or. &
      Uninitialized(this%chi)) then
    option%io_buffer = 'RMAX, K_DONOR, K_ACCEPTOR, Y, M, CHI, DG0, and O2_THRESHOLD must be set for &
      JINBETHKE_O2AQ_DOC.'
    call PrintErrMsg(option)
  endif
end subroutine JinBethkeO2aqDOCReadInput
! ************************************************************************** !
subroutine JinBethkeO2aqDOCSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Reaction_Immobile_Aux_module, only: GetImmobileSpeciesIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_jinbethke_o2aq_doc_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 3
  ! Aqueous species
  word = 'DOC-'
  this%doc_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H+'
  this%h_ion_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%bicarbonate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Xim'
  this%xim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

end subroutine JinBethkeO2aqDOCSetup
! ************************************************************************** !
subroutine O2aqAuxiliaryPlotVariables(this,list,reaction,option)
  !
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_jinbethke_o2aq_doc_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'JB O2(aq) resp. Sandbox Rate'
  units = 'mol-C/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  word = 'dG-rxn_O2aq-resp Sandbox Generic'
  units = 'kJ/mol-C'
  call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+2)

  word = 'Ft_O2aq-resp Sandbox Generic'
  units = 'unitless'
  call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+3)
end subroutine O2aqAuxiliaryPlotVariables
! ************************************************************************** !
subroutine JinBethkeO2aqDOCEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !Jacobian,compute_derivative,
  ! Evaluates aerobic respiration reaction storing residual but no Jacobian
  !
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/10/19

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  implicit none
  class(reaction_sandbox_jinbethke_o2aq_doc_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

  type(material_auxvar_type) :: material_auxvar
  PetscInt, parameter :: iphase = 1
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: porosity             ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: L_water              ! L water

  PetscReal :: Bicarbonate, Xim, yield, O2aq, DOC, Proton
  PetscReal :: Rate, Rate_DOC, Rate_Bicarbonate, Rate_O2aq, Rate_Proton
  PetscReal :: stoi_doc, stoi_bicarbonate, stoi_o2aq, stoi_proton
  PetscReal :: m, chi, k_rmax
  PetscReal :: temp_K, RT
  PetscReal :: Ft, Ftr, Fdonor, Facceptor
  PetscReal :: reaction_Q
  PetscReal :: dG0, dGr, dG_ATP

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: imnrl
  PetscInt :: iauxiliary

  PetscBool :: calculate_rate

  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  ! Rxn:    1.00 CHO- + 1.00 O2(aq) = 1.00 HCO3-  + 1.00 H+ 

  DOC = rt_auxvar%pri_molal(this%doc_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%doc_id) 
  Proton = rt_auxvar%pri_molal(this%h_ion_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  Bicarbonate = rt_auxvar%pri_molal(this%bicarbonate_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%bicarbonate_id) 
  O2aq = rt_auxvar%pri_molal(this%o2aq_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%o2aq_id) 
  Xim = rt_auxvar%immobile(this%xim_id)

  m = this%m
  chi = this%chi
  k_rmax = this%rmax
  
  stoi_bicarbonate = 1.d0
  stoi_doc = 1.d0
  stoi_o2aq = 1.d0 
  stoi_proton = 1.d0

  RT = (8.314e-3) * (global_auxvar%temp + 273.15d0)
  dG0 = this%dg0 ! kJ / mol C; set by user
  dG_ATP = 50.d0 ! kJ / mol ATP

  reaction_Q = ( (Proton**stoi_proton) * (Bicarbonate**stoi_bicarbonate)) / &
    ((DOC**stoi_doc) * (O2aq**stoi_o2aq))

  dGr = dG0 + RT*log(reaction_Q)
  
  ! Monod expressions for acetate
  Fdonor = DOC / (DOC + this%Kdonor)

  ! Monod expressions for oxygen
  Facceptor = O2aq / (O2aq + this%Kacceptor)

  ! Thermodynamic factor 
  Ft = 1.d0 - exp((dGr + m*dG_ATP) / (chi * RT))

  if (Ft < 0) then 
    Ftr = 0.d0
  else
    Ftr = Ft
  endif

  calculate_rate = (O2aq > (this%o2_threshold))

  Rate = 0.d0
  
  if (calculate_rate) then
    ! base rate, mol/sec/m^3 bulk
    ! units on k: mol/sec/mol-bio

    Rate = -k_rmax * Ftr * Facceptor * Fdonor * Xim * L_water  ! mol/sec
    ! species-specifc 
    Rate_DOC = Rate * stoi_doc
    Rate_Proton = Rate * stoi_proton 
    Rate_Bicarbonate = Rate * stoi_bicarbonate 
    Rate_O2aq = Rate * stoi_o2aq
    !Rate_xim = Rate * yield

    rt_auxvar%auxiliary_data(iauxiliary) = Rate_DOC
    rt_auxvar%auxiliary_data(iauxiliary+1) = dGr
    rt_auxvar%auxiliary_data(iauxiliary+2) = Ftr
    
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate_Proton
    Residual(this%doc_id) = Residual(this%doc_id) - Rate_DOC
    Residual(this%o2aq_id) = Residual(this%o2aq_id) - Rate_O2aq
    Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + Rate_Bicarbonate
    
  endif

end subroutine JinBethkeO2aqDOCEvaluate

end module Reaction_Sandbox_JinBethke_O2aq_DOC_class