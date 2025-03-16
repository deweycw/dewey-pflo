module Reaction_Sandbox_Fe_Mineral_Precipitation_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class 
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_fe_mineral_precipitation_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: fe3_id
    PetscInt :: mineral_id
    PetscReal :: rate_constant

  contains
    procedure, public :: ReadInput => FeMineralPrecipitationReadInput
    procedure, public :: Setup => FeMineralPrecipitationSetup
    procedure, public :: AuxiliaryPlotVariables => FeMineralPrecipitationAuxiliaryPlotVariables
    procedure, public :: Evaluate => FeMineralPrecipitationEvaluate
    procedure, public :: UpdateKineticState => FeMineralPrecipitationUpdateKineticState
  end type reaction_sandbox_fe_mineral_precipitation_type

  public :: FeMineralPrecipitationCreate, &
            FeMineralPrecipitationSetup
contains
! ************************************************************************** !
function FeMineralPrecipitationCreate()
  !
  !
  implicit none
  class(reaction_sandbox_fe_mineral_precipitation_type), pointer :: FeMineralPrecipitationCreate
  allocate(FeMineralPrecipitationCreate)
  FeMineralPrecipitationCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  FeMineralPrecipitationCreate%h_ion_id = UNINITIALIZED_INTEGER
  FeMineralPrecipitationCreate%fe3_id = UNINITIALIZED_INTEGER
  FeMineralPrecipitationCreate%mineral_id = UNINITIALIZED_INTEGER

  FeMineralPrecipitationCreate%rate_constant = UNINITIALIZED_DOUBLE

  nullify(FeMineralPrecipitationCreate%next)
end function FeMineralPrecipitationCreate
! ************************************************************************** !
subroutine FeMineralPrecipitationReadInput(this,input,option)
  !
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_fe_mineral_precipitation_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FE_MINERAL_PRECIPITATION'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be set for &
      FE_MINERAL_PRECIPITATION.'
    call PrintErrMsg(option)
  endif
end subroutine FeMineralPrecipitationReadInput
! ************************************************************************** !
subroutine FeMineralPrecipitationSetup(this,reaction,option)
  !
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_fe_mineral_precipitation_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1
  ! Aqueous species
  word = 'H+'
  this%h_ion_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Fe+++'
  this%fe3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Fe(OH)3(s)'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine FeMineralPrecipitationSetup
! ************************************************************************** !
subroutine FeMineralPrecipitationAuxiliaryPlotVariables(this,list,reaction,option)
  !
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_fe_mineral_precipitation_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'Fe(OH)3(s) Precip. Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine FeMineralPrecipitationAuxiliaryPlotVariables
! ************************************************************************** !
subroutine FeMineralPrecipitationEvaluate(this, Residual,Jacobian,compute_derivative, &
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
  use Material_Aux_module
  use Reaction_Mineral_Aux_module
  implicit none
  class(reaction_sandbox_fe_mineral_precipitation_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt, parameter :: iphase = 1
  type(mineral_type), pointer :: mineral
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: porosity             ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: L_water              ! L water

  PetscReal :: Proton, Fe3
  PetscReal :: Rate, Rate_Proton, Rate_Fe3
  PetscReal :: stoi_proton, stoi_fe3
  PetscReal :: k_precip
  PetscReal :: lnQK, sign_, affinity_factor, QK
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: imnrl
  PetscInt :: iauxiliary

  PetscBool :: calculate_precip

  mineral => reaction%mineral

  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
  imnrl = this%mineral_id

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  k_precip = this%rate_constant

    ! Rxn:    1.00 Fe+++ + 3.00 H2O = 1.00 Fe(OH)3(s) + 3.00 H+ 

  Proton = rt_auxvar%pri_molal(this%h_ion_id) * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  Fe3 = rt_auxvar%pri_molal(this%fe3_id) * &
    rt_auxvar%pri_act_coef(this%fe3_id) 


  stoi_fe3 = 1.d0
  stoi_proton = 3.d0 

  
  ! TST for precipitation; uses pkeq from database
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN

  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)* &
                  rt_auxvar%ln_act_h2o
  endif

  ! activity of other species
  ncomp = mineral%kinmnrlspecid(0,imnrl)
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo

  QK = exp(lnQK)
  affinity_factor = 1.d0-QK
  sign_ = sign(1.d0,affinity_factor)

  calculate_precip = (sign_<0)

  Rate = 0.d0

  if (calculate_precip) then
    
    Rate = this%rate_constant * Fe3

    rt_auxvar%auxiliary_data(iauxiliary) = Rate

    Rate = Rate * material_auxvar%volume ! mol/sec
      
    ! species-specifc 
    Rate_Fe3 = Rate 
    Rate_Proton = Rate * (3.d0) 
    
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate_Proton
    Residual(this%fe3_id) = Residual(this%fe3_id) - Rate_Fe3
  endif

end subroutine FeMineralPrecipitationEvaluate
! ************************************************************************** !
subroutine FeMineralPrecipitationUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
  !
  ! Updates mineral volume fraction at end converged timestep based on latest
  ! rate
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  implicit none
  class(reaction_sandbox_fe_mineral_precipitation_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  PetscInt :: imnrl
  PetscReal :: delta_volfrac
  imnrl = this%mineral_id
  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
  !                                mol_vol (m^3 mnrl/mol mnrl)
  delta_volfrac = rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                  option%tran_dt
  ! m^3 mnrl/m^3 bulk
  rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                  delta_volfrac
  ! zero to avoid negative volume fractions
  if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
    rt_auxvar%mnrl_volfrac(imnrl) = 0.d0
end subroutine FeMineralPrecipitationUpdateKineticState
end module Reaction_Sandbox_Fe_Mineral_Precipitation_class