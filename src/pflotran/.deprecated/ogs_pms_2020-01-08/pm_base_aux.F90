module PM_Base_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use AuxVars_Base_module
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  type, public :: pm_base_aux_type 
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(matrix_zeroing_type), pointer :: matrix_zeroing

    !if required to operate with auxvar_base objects use pointers below
    !they must point to doughter classes auxvars (e.g. see: well_flow, well_flow_eenergy) 
    !class(auxvar_base_type), pointer :: auxvars_base(:,:)
    !class(auxvar_base_type), pointer :: auxvars_bc_base(:)
    !class(auxvar_base_type), pointer :: auxvars_ss_base(:)
    !gfrotran 4.8.4 bugs. class arrays cannot be passed corretly as dummy arguments
    ! workaround: covert type ouside calling fucntion using a class pointers
    ! and use compatible "type()" to declair the dummy arguments in the called function. 
  contains
    !add here type-bound-procedure
  end type pm_base_aux_type

 
  public :: PMBaseAuxInit, PMBaseAuxSetup, &
            PMBaseAuxStrip 

contains

! ************************************************************************** !

subroutine PMBaseAuxInit(this)
  ! 
  ! Initialize pm auvars 
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/27/16
  ! 

  implicit none

  class(pm_base_aux_type) :: this

  this%auxvars_up_to_date = PETSC_FALSE
  this%inactive_cells_exist = PETSC_FALSE
  this%num_aux = 0
  this%num_aux_bc = 0
  this%num_aux_ss = 0
  !nullify(this%auxvars_base)
  !nullify(this%auxvars_bc_base)
  !nullify(this%auxvars_ss_base)
  nullify(this%matrix_zeroing)

end subroutine PMBaseAuxInit

! ************************************************************************** !

subroutine PMBaseAuxSetup(this,grid,option)
  ! 
  ! Initialize pm_aux  
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/27/16
  ! 

  use Option_module
  use Grid_module 

  implicit none

  class(pm_base_aux_type) :: this
  type(grid_type) :: grid
  type(option_type) :: option

  call MatrixZeroingInitRowZeroing(this%matrix_zeroing,grid%nlmax)

  !note
  !inactive_rows_local and inactive_rows_local_ghosted
  ! are allocated in InitSubsurfaceCreateZeroArray 

end subroutine PMBaseAuxSetup

! ************************************************************************** !

subroutine PMBaseAuxStrip(this)
  ! 
  ! destroy pm_base_aux  
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/30/16
  ! 

  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_base_aux_type) :: this

  call MatrixZeroingDestroy(this%matrix_zeroing)
 
end subroutine PMBaseAuxStrip

! ************************************************************************** !

end module PM_Base_Aux_module
