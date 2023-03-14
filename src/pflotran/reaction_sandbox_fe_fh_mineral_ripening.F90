module Reaction_Sandbox_Fe_Fh_Mineral_Ripening_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class 
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_fe_fh_mineral_ripening_type
    PetscInt :: auxiliary_offset
    PetscInt :: fe2_id
    PetscInt :: fh_mineral_id
    PetscInt :: feoh_s_mineral_id
    PetscReal :: rate_constant_feoh_fh

  contains
    procedure, public :: ReadInput => FeFhMineralRipeningReadInput
    procedure, public :: Setup => FeFhMineralRipeningSetup
    procedure, public :: AuxiliaryPlotVariables => FeFhMineralRipeningAuxiliaryPlotVariables
    procedure, public :: Evaluate => FeFhMineralRipeningEvaluate
    procedure, public :: UpdateKineticState => FeFhMineralRipeningUpdateKineticState
  end type reaction_sandbox_fe_fh_mineral_ripening_type

  public :: FeFhMineralRipeningCreate, &
            FeFhMineralRipeningSetup
contains
! ************************************************************************** !
function FeFhMineralRipeningCreate()
  !
  implicit none
  class(reaction_sandbox_fe_fh_mineral_ripening_type), pointer :: FeFhMineralRipeningCreate
  allocate(FeFhMineralRipeningCreate)
  FeFhMineralRipeningCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  FeFhMineralRipeningCreate%fe2_id = UNINITIALIZED_INTEGER
  FeFhMineralRipeningCreate%feoh_s_mineral_id = UNINITIALIZED_INTEGER
  FeFhMineralRipeningCreate%fh_mineral_id = UNINITIALIZED_INTEGER
  FeFhMineralRipeningCreate%rate_constant_feoh_fh = UNINITIALIZED_DOUBLE
  nullify(FeFhMineralRipeningCreate%next)
end function FeFhMineralRipeningCreate
! ************************************************************************** !
subroutine FeFhMineralRipeningReadInput(this,input,option)
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_fe_fh_mineral_ripening_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FE_FH_MINERAL_RIPENING'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('RATE_CONSTANT_FE_TO_FH')
        call InputReadDouble(input,option,this%rate_constant_feoh_fh)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rate_constant_feoh_fh)) then
    option%io_buffer = 'RATE_CONSTANT_FE_TO_FH must be set for &
      FE_FH_MINERAL_RIPENING.'
    call PrintErrMsg(option)
  endif
end subroutine FeFhMineralRipeningReadInput
! ************************************************************************** !
subroutine FeFhMineralRipeningSetup(this,reaction,option)
  !
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_fe_fh_mineral_ripening_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1

  word = 'Fe++'
  this%fe2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Fe(OH)3(s)'
  this%feoh_s_mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)
  word = 'Ferrihydrite'
  this%fh_mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine FeFhMineralRipeningSetup
! ************************************************************************** !
subroutine FeFhMineralRipeningAuxiliaryPlotVariables(this,list,reaction,option)
  !
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_fe_fh_mineral_ripening_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'Fe(OH)3(s) to Fh Ripening Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine FeFhMineralRipeningAuxiliaryPlotVariables
! ************************************************************************** !
subroutine FeFhMineralRipeningEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)

  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/3/13


  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Mineral_Aux_module
  implicit none
  class(reaction_sandbox_fe_fh_mineral_ripening_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt, parameter :: iphase = 1
  type(mineral_type), pointer :: mineral
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: porosity             ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: L_water              ! L water

  PetscReal :: Fe2
  PetscReal :: Rate
  PetscReal :: k_ripen

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: imnrl, jmnrl
  PetscInt :: iauxiliary

  mineral => reaction%mineral

  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  imnrl = this%feoh_s_mineral_id
  jmnrl = this%fh_mineral_id

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  k_ripen = this%rate_constant_feoh_fh

  Fe2 = rt_auxvar%pri_molal(this%fe2_id) * &
    rt_auxvar%pri_act_coef(this%fe2_id) 
    
  Rate = k_ripen * Fe2
  
  !Rate = Rate * material_auxvar%volume ! mol/sec

  rt_auxvar%auxiliary_data(iauxiliary) = Rate


end subroutine FeFhMineralRipeningEvaluate
! ************************************************************************** !
subroutine FeFhMineralRipeningUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
  !
  ! Updates mineral volume fraction at end converged timestep based on latest
  ! rate
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  implicit none
  class(reaction_sandbox_fe_fh_mineral_ripening_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  PetscInt :: imnrl, jmnrl 
  PetscReal :: delta_volfrac_imnrl, delta_volfrac_jmnrl
  imnrl = this%feoh_s_mineral_id
  jmnrl = this%fh_mineral_id
  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
  !                                mol_vol (m^3 mnrl/mol mnrl)
  delta_volfrac_imnrl = (-1) * rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                  option%tran_dt

  delta_volfrac_jmnrl = rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(jmnrl)* &
                  option%tran_dt
  ! m^3 mnrl/m^3 bulk
  rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                  delta_volfrac_imnrl

  rt_auxvar%mnrl_volfrac(jmnrl) = rt_auxvar%mnrl_volfrac(jmnrl) + &
                                  delta_volfrac_jmnrl
  ! zero to avoid negative volume fractions
  if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
    rt_auxvar%mnrl_volfrac(imnrl) = 0.d0
end subroutine FeFhMineralRipeningUpdateKineticState
end module Reaction_Sandbox_Fe_Fh_Mineral_Ripening_class