module PMC_Subsurface_OSRT_class

#include "petsc/finclude/petscts.h"
  use petscts

  use PMC_Subsurface_class
  use PM_OSRT_class

  use PFLOTRAN_Constants_module

  implicit none

  
  private

  type, public, extends(pmc_subsurface_type) :: pmc_subsurface_osrt_type
  contains
    procedure, public :: Init => PMCSubsurfaceOSRTInit
    procedure, public :: SetupSolvers => PMCSubsurfaceOSRTSetupSolvers
    procedure, public :: StepDT => PMCSubsurfaceOSRTStepDT
    procedure, public :: Destroy => PMCSubsurfaceOSRTDestroy
  end type pmc_subsurface_osrt_type
  
  public :: PMCSubsurfaceOSRTCreate, &
            PMCSubsurfaceOSRTInit
  
contains

! ************************************************************************** !

function PMCSubsurfaceOSRTCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  implicit none
  
  class(pmc_subsurface_osrt_type), pointer :: PMCSubsurfaceOSRTCreate
  
  class(pmc_subsurface_osrt_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceOSRTCreate => pmc  
  
end function PMCSubsurfaceOSRTCreate

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  
  implicit none
  
  class(pmc_subsurface_osrt_type) :: this
  
  call PMCSubsurfaceInit(this)
  this%name = 'PMCSubsurfaceOSRT'

end subroutine PMCSubsurfaceOSRTInit

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTSetupSolvers(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module
  use Solver_module
  use Discretization_module
  use PM_RT_class
  use Timestepper_KSP_class

  implicit none

  class(pmc_subsurface_osrt_type) :: this

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver
  character(len=MAXSTRINGLENGTH) :: string
  class(pm_rt_type), pointer :: pm_rt
  PetscErrorCode :: ierr

  option => this%option
  solver => this%timestepper%solver

  select type(ts=>this%timestepper) 
    class is(timestepper_KSP_type)
    class default
      option%io_buffer = 'A KSP timestepper must be used for operator-split &
        &reactive transport'
      call PrintErrMsg(option)
  end select

  select type(pm=>this%pm_ptr%pm)
    class is(pm_rt_type)
      pm_rt => pm
  end select

  call SolverCreateKSP(solver,option%mycomm)

  call PrintMsg(option,"  Beginning setup of TRAN KSP")
  call KSPSetOptionsPrefix(solver%ksp, "tran_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  solver%M_mat_type = MATAIJ
  solver%Mpre_mat_type = MATAIJ
  call DiscretizationCreateMatrix(pm_rt%realization%discretization, &
                                  ONEDOF, &
                                  solver%Mpre_mat_type, &
                                  solver%Mpre,option)
  call MatSetOptionsPrefix(solver%Mpre,"tran_",ierr);CHKERRQ(ierr)
  solver%M = solver%Mpre

  ! Have PETSc do a KSP_View() at the end of each solve if 
  ! verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-tran_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
  endif

  call SolverSetKSPOptions(solver,option)

end subroutine PMCSubsurfaceOSRTSetupSolvers

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTStepDT(this,stop_flag)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Solver_module
  use Field_module
  use Reactive_Transport_Aux_module
  use Reactive_Transport_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Aux_module
  use Reaction_module
  use PM_RT_class
  use Debug_module
  use Timestepper_Base_class
  use Timestepper_KSP_class
  use Output_module
  use String_module

  implicit none

  class(pmc_subsurface_osrt_type) :: this
  PetscInt :: stop_flag

  type(option_type), pointer :: option
  class(pm_osrt_type), pointer :: process_model
  class(realization_subsurface_type), pointer :: realization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(solver_type), pointer :: solver
  class(timestepper_KSP_type), pointer :: timestepper
  class(reaction_rt_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt, parameter :: iphase = 1

  PetscViewer :: viewer

  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: tran_xx_p(:)
  PetscInt :: idof
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset_global
  PetscInt :: istart, iend
  PetscInt :: nimmobile
  PetscInt :: num_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: sum_linear_iterations
  PetscInt :: sum_linear_iterations_temp
  PetscInt :: sum_wasted_linear_iterations
  PetscInt :: icut
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  PetscInt :: rreact_error
  PetscInt :: tempint
  PetscReal :: tconv
  PetscReal :: tempreal

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: tunit

  KSPConvergedReason :: ksp_reason
  PetscLogDouble :: log_outer_start_time
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_ksp_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr

  call PetscTime(log_outer_start_time,ierr);CHKERRQ(ierr)
  call this%PrintHeader()

  option => this%option
  realization => this%realization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  reaction => patch%reaction
  select type(ts=>this%timestepper)
    class is(timestepper_KSP_type)
      timestepper => ts
  end select
  solver => timestepper%solver
  select type(pm_osrt=>this%pm_list)
    class is(pm_osrt_type)
      process_model => pm_osrt
  end select

  material_auxvars => patch%aux%Material%auxvars
  global_auxvars => patch%aux%Global%auxvars
  rt_auxvars => patch%aux%RT%auxvars

  nimmobile = reaction%immobile%nimmobile

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_linear_iterations_temp = 0
  sum_wasted_linear_iterations = 0
  icut = 0

  option%dt = timestepper%dt
  option%time = timestepper%target_time - timestepper%dt

  call process_model%InitializeTimestep()

  ! from RTUpdateRHSCoefs
  ! update time derivative on RHS
#if 0
  call VecGetArrayF90(process_model%rhs_coef,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    vec_ptr(local_id) = material_auxvars(ghosted_id)%porosity* &
                        global_auxvars(ghosted_id)%sat(iphase)* &
                        1000.d0* &
                        material_auxvars(ghosted_id)%volume
  enddo
  call VecRestoreArrayF90(process_model%rhs_coef,vec_ptr,ierr);CHKERRQ(ierr)
#else
  call VecGetArrayF90(process_model%fixed_accum,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset_global = (local_id-1)*reaction%ncomp
    istart = offset_global + 1
    iend = offset_global + reaction%naqcomp
    vec_ptr(istart:iend) = material_auxvars(ghosted_id)%porosity* &
                           global_auxvars(ghosted_id)%sat(iphase)* &
                           1000.d0* &
                           material_auxvars(ghosted_id)%volume* &
                           rt_auxvars(ghosted_id)%total(:,iphase)
    !TODO(geh): do immobile need to be stored; they are in tran_xx
    if (nimmobile > 0) then
      ! need to store immobile concentrations for the purpose of reverting
      ! when reaction fails
      istart = istart + reaction%offset_immobile
      iend = istart + nimmobile - 1
      vec_ptr(istart:iend) = rt_auxvars(ghosted_id)%immobile
    endif
  enddo
  call VecRestoreArrayF90(process_model%fixed_accum,vec_ptr,ierr);CHKERRQ(ierr)
#endif

  do

    call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

    rreact_error = 0

!    call RTCalculateRHS_t0(realization)

    call PMRTWeightFlowParameters(process_model,TIME_TpDT)
    ! update diffusion/dispersion coefficients
    call RTUpdateTransportCoefs(realization)
    ! RTCalculateRHS_t1() updates aux vars (cell and boundary) to k+1 
    ! and calculates RHS fluxes and src/sinks
    call VecCopy(process_model%fixed_accum,process_model%rhs,ierr);CHKERRQ(ierr)
    tempreal = 1.d0/option%tran_dt
    call VecScale(process_model%rhs,tempreal,ierr);CHKERRQ(ierr)
    call RTCalculateRHS_t1(realization,process_model%rhs)

    ! RTCalculateTransportMatrix() calculates flux coefficients and the
    ! t^(k+1) coefficient in accumulation term
    if (rt_parameter%species_dependent_diffusion) then
      option%io_buffer = 'Operator-split reactive transport is not &
        &currently configured to handle species-dependent diffusion.'
      call PrintErrMsg(option)
    endif
    call RTCalculateTransportMatrix(realization,solver%M)
    call KSPSetOperators(solver%ksp,solver%M,solver%Mpre,ierr);CHKERRQ(ierr)

    ! loop over chemical component and transport
    do idof = 1, rt_parameter%naqcomp

      tempint = idof-1
      call VecStrideGather(process_model%rhs,tempint,field%work, &
                           INSERT_VALUES,ierr);CHKERRQ(ierr)

!debug      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
!debug      print *, 'Trhs: ', trim(StringWrite(vec_ptr))
!debug      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

      if (realization%debug%vecview_residual) then
        string = 'Trhs'
        call DebugCreateViewer(realization%debug,string,option,viewer)
        call VecView(field%work,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
      endif

      call PetscTime(log_ksp_start_time,ierr);CHKERRQ(ierr)      
      call KSPSolve(solver%ksp,field%work,field%work,ierr);CHKERRQ(ierr)
      call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
      timestepper%cumulative_solver_time = &
        timestepper%cumulative_solver_time + (log_end_time - log_ksp_start_time)
      
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
!debug      print *, 'Tsol: ', trim(StringWrite('(es17.8)',vec_ptr))
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle
        rt_auxvars(ghosted_id)%total(idof,iphase) = vec_ptr(local_id)
      enddo
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call KSPGetIterationNumber(solver%ksp,num_linear_iterations,ierr)
      call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
      sum_linear_iterations_temp = sum_linear_iterations_temp + &
        num_linear_iterations
    enddo

    if (option%compute_mass_balance_new) then
      call RTZeroMassBalanceDelta(realization)
      call RTComputeBCMassBalanceOS(realization)
    endif

    sum_linear_iterations = sum_linear_iterations + sum_linear_iterations_temp
    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
    process_model%cumulative_transport_time = &
      process_model%cumulative_transport_time + log_end_time - log_start_time
    log_start_time = log_end_time

    ! react all chemical components
    call VecGetArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset_global = (local_id-1)*reaction%ncomp
      istart = offset_global + 1
      iend = offset_global + reaction%ncomp
      if (nimmobile > 0) then
        tempint = istart+reaction%offset_immobile
        rt_auxvars(ghosted_id)%immobile = tran_xx_p(tempint:tempint+nimmobile-1)
      endif
      call RReact(tran_xx_p(istart:iend),rt_auxvars(ghosted_id), &
                  global_auxvars(ghosted_id),material_auxvars(ghosted_id), &
                  num_iterations,reaction,grid%nG2A(ghosted_id),option, &
                  PETSC_TRUE,PETSC_TRUE,rreact_error)
      if (rreact_error /= 0) exit
      ! set primary dependent var back to free-ion molality
      iend = offset_global + reaction%naqcomp
      tran_xx_p(istart:iend) = rt_auxvars(ghosted_id)%pri_molal(:)
      if (nimmobile > 0) then
        tempint = istart+reaction%offset_immobile
        tran_xx_p(tempint:tempint+nimmobile-1) = rt_auxvars(ghosted_id)%immobile
      endif
    enddo
    call VecRestoreArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)

    call MPI_Allreduce(MPI_IN_PLACE,rreact_error,ONE_INTEGER_MPI, &
                       MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
    call MPI_Barrier(option%mycomm,ierr)
    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
    process_model%cumulative_reaction_time = &
      process_model%cumulative_reaction_time + log_end_time - log_start_time

    if (rreact_error /= 0) then
      !TODO(geh): move to timestepper base and call from daughters.
      sum_wasted_linear_iterations = sum_wasted_linear_iterations + &
        sum_linear_iterations_temp
      icut = icut + 1
      timestepper%time_step_cut_flag = PETSC_TRUE
      ! if a cut occurs on the last time step, the stop_flag will have been
      ! set to TS_STOP_END_SIMULATION.  Set back to TS_CONTINUE to prevent
      ! premature ending of simulation.
      if (stop_flag /= TS_STOP_MAX_TIME_STEP) stop_flag = TS_CONTINUE

      if (icut > timestepper%max_time_step_cuts .or. &
          timestepper%dt < timestepper%dt_min) then

        if (icut > timestepper%max_time_step_cuts) then
          option%io_buffer = ' Stopping: Time step cut criteria exceeded.'
          call PrintMsg(option)
          write(option%io_buffer, &
                '("    icut =",i3,", max_time_step_cuts=",i3)') &
                icut,timestepper%max_time_step_cuts
          call PrintMsg(option)
        endif
        if (timestepper%dt < timestepper%dt_min) then
          option%io_buffer = ' Stopping: Time step size is less than the &
                             &minimum allowable time step.'
          call PrintMsg(option)
          write(option%io_buffer, &
                '("    dt   =",es15.7,", dt_min=",es15.7," [",a,"]")') &
               timestepper%dt/tconv,timestepper%dt_min/tconv,trim(tunit)
          call PrintMsg(option)
       endif

        ! print snapshot only 
        process_model%output_option%plot_name = trim(process_model%name) // &
          '_cut_to_failure'
        snapshot_plot_flag = PETSC_TRUE
        observation_plot_flag = PETSC_FALSE
        massbal_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,snapshot_plot_flag, &
                    observation_plot_flag,massbal_plot_flag)
        stop_flag = TS_STOP_FAILURE
        return
      endif

      timestepper%target_time = this%timestepper%target_time - timestepper%dt

      timestepper%dt = timestepper%time_step_reduction_factor * timestepper%dt

      write(option%io_buffer,'(''-> Cut time step: icut='',i2, &
           &   ''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
           &   1pe12.5)')  icut, &
           timestepper%cumulative_time_step_cuts+icut, &
           option%time/tconv, &
           timestepper%dt/tconv
      call PrintMsg(option)
      timestepper%target_time = timestepper%target_time + timestepper%dt
      option%dt = timestepper%dt
      call process_model%TimeCut()
    else
      exit
    endif

  enddo

  timestepper%steps = timestepper%steps + 1
  timestepper%cumulative_linear_iterations = &
    timestepper%cumulative_linear_iterations + sum_linear_iterations
  timestepper%cumulative_wasted_linear_iterations = &
    timestepper%cumulative_wasted_linear_iterations + &
    sum_wasted_linear_iterations

  if (option%print_screen_flag) then
      write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
           & " [",a,"]", " cuts = ",i2, " [",i4,"]")') &
           this%timestepper%steps, &
           this%timestepper%target_time/tconv, &
           this%timestepper%dt/tconv,trim(tunit),icut, &
           this%timestepper%cumulative_time_step_cuts
  endif

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
  this%cumulative_time = this%cumulative_time + &
    log_end_time - log_outer_start_time

end subroutine PMCSubsurfaceOSRTStepDT

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  
  implicit none
  
  class(pmc_subsurface_osrt_type) :: this

  call PMCSubsurfaceStrip(this)
  nullify(this%realization)

end subroutine PMCSubsurfaceOSRTStrip

! ************************************************************************** !

recursive subroutine PMCSubsurfaceOSRTDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  use Option_module

  implicit none
  
  class(pmc_subsurface_osrt_type) :: this
  
  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif
  
  !TODO(geh): place this routine in PMC_Base_class and redirect Strip() to
  !           avoid creating all these Destroy routines
  call PMCSubsurfaceOSRTStrip(this)
  
end subroutine PMCSubsurfaceOSRTDestroy

! ************************************************************************** !
  
end module PMC_Subsurface_OSRT_class
