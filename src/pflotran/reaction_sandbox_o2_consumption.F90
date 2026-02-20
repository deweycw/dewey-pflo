module Reaction_Sandbox_O2_Consumption_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private

! ************************************************************************** !
! Combined O2-consuming reactions sandbox
!
! Handles aerobic respiration, Fe2+ oxidation, and HS- oxidation in a
! single module with coordinated O2 consumption to avoid Newton solver
! stiffness.
!
! Reactions:
!   SOC(aq) + O2(aq) + H2O -> HCO3- + H+         (aerobic respiration)
!   Fe2+ + 0.25 O2(aq) + H+ -> Fe3+ + 0.5 H2O    (Fe2+ oxidation)
!   HS-  + 2.0  O2(aq)       -> SO4-- + H+         (HS- oxidation)
!
! Design:
!   - Uses smooth sigmoid O2 activation (not Monod on O2) to avoid
!     stiff Jacobian entries near the half-saturation concentration
!   - Substrate limitation via Monod for SOC, Fe2+, and HS- only
!   - All three reactions evaluated together with coordinated O2 cap
!     to prevent over-consumption and Newton oscillation
!   - Hard cutoff on f_o2 at 1e-20 to avoid denormalized float issues
!
! Input block:
!   REACTION_SANDBOX
!     O2_CONSUMPTION
!       AERO_VMAX         5.d-10   ! mol/L-water/s  aerobic respiration
!       AERO_HALF_SAT     1.d-5    ! mol/L SOC(aq)
!       FE_VMAX           5.d-9    ! mol/L-water/s
!       FE_HALF_SAT       1.d-5    ! mol/L Fe2+
!       HS_VMAX           2.d-8    ! mol/L-water/s
!       HS_HALF_SAT       1.d-5    ! mol/L HS-
!       O2_THRESHOLD      1.d-5    ! mol/L - sigmoid midpoint
!       O2_SAFETY_FACTOR  0.8d0    ! max fraction of O2 consumed per dt
!     /
!   /
!
! Author: Christian Dewey
! Date: 2026/02/19
! ************************************************************************** !

  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_o2_consumption_type
    ! Species IDs
    PetscInt :: auxiliary_offset
    PetscInt :: fe2_id
    PetscInt :: fe3_id
    PetscInt :: o2_id
    PetscInt :: h_id
    PetscInt :: hs_id
    PetscInt :: so4_id
    PetscInt :: soc_id
    PetscInt :: hco3_id
    ! Aerobic respiration parameters
    PetscReal :: aero_vmax
    PetscReal :: aero_half_sat
    ! Fe2+ oxidation parameters
    PetscReal :: fe_vmax
    PetscReal :: fe_half_sat
    ! HS- oxidation parameters
    PetscReal :: hs_vmax
    PetscReal :: hs_half_sat
    ! Shared O2 parameters
    PetscReal :: o2_threshold
    PetscReal :: o2_safety_factor

  contains
    procedure, public :: ReadInput => O2ConsumptionReadInput
    procedure, public :: Setup => O2ConsumptionSetup
    procedure, public :: AuxiliaryPlotVariables => O2ConsumptionAuxiliaryPlotVariables
    procedure, public :: Evaluate => O2ConsumptionEvaluate
  end type reaction_sandbox_o2_consumption_type

  public :: O2ConsumptionCreate, &
            O2ConsumptionSetup
contains

! ************************************************************************** !
function O2ConsumptionCreate()
  !
  ! Allocates O2 consumption reaction object.
  !
  implicit none
  class(reaction_sandbox_o2_consumption_type), pointer :: O2ConsumptionCreate
  allocate(O2ConsumptionCreate)
  O2ConsumptionCreate%auxiliary_offset = UNINITIALIZED_INTEGER

  O2ConsumptionCreate%fe2_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%fe3_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%o2_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%h_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%hs_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%so4_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%soc_id = UNINITIALIZED_INTEGER
  O2ConsumptionCreate%hco3_id = UNINITIALIZED_INTEGER

  O2ConsumptionCreate%aero_vmax = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%aero_half_sat = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%fe_vmax = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%fe_half_sat = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%hs_vmax = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%hs_half_sat = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%o2_threshold = UNINITIALIZED_DOUBLE
  O2ConsumptionCreate%o2_safety_factor = 0.8d0  ! default

  nullify(O2ConsumptionCreate%next)
end function O2ConsumptionCreate

! ************************************************************************** !
subroutine O2ConsumptionReadInput(this,input,option)
  !
  ! Reads O2 consumption reaction parameters from input file
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_o2_consumption_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,O2_CONSUMPTION'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('AERO_VMAX')
        call InputReadDouble(input,option,this%aero_vmax)
        call InputErrorMsg(input,option,word,error_string)
      case('AERO_HALF_SAT')
        call InputReadDouble(input,option,this%aero_half_sat)
        call InputErrorMsg(input,option,word,error_string)
      case('FE_VMAX')
        call InputReadDouble(input,option,this%fe_vmax)
        call InputErrorMsg(input,option,word,error_string)
      case('FE_HALF_SAT')
        call InputReadDouble(input,option,this%fe_half_sat)
        call InputErrorMsg(input,option,word,error_string)
      case('HS_VMAX')
        call InputReadDouble(input,option,this%hs_vmax)
        call InputErrorMsg(input,option,word,error_string)
      case('HS_HALF_SAT')
        call InputReadDouble(input,option,this%hs_half_sat)
        call InputErrorMsg(input,option,word,error_string)
      case('O2_THRESHOLD')
        call InputReadDouble(input,option,this%o2_threshold)
        call InputErrorMsg(input,option,word,error_string)
      case('O2_SAFETY_FACTOR')
        call InputReadDouble(input,option,this%o2_safety_factor)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%aero_vmax) .or. &
      Uninitialized(this%aero_half_sat) .or. &
      Uninitialized(this%fe_vmax) .or. &
      Uninitialized(this%fe_half_sat) .or. &
      Uninitialized(this%hs_vmax) .or. &
      Uninitialized(this%hs_half_sat) .or. &
      Uninitialized(this%o2_threshold)) then
    option%io_buffer = 'AERO_VMAX, AERO_HALF_SAT, FE_VMAX, FE_HALF_SAT, ' // &
      'HS_VMAX, HS_HALF_SAT, and O2_THRESHOLD must be set for O2_CONSUMPTION.'
    call PrintErrMsg(option)
  endif
end subroutine O2ConsumptionReadInput

! ************************************************************************** !
subroutine O2ConsumptionSetup(this,reaction,option)
  !
  ! Maps species names to internal IDs
  !
  use Reaction_Aux_module, only : reaction_rt_type, &
                                  ReactionAuxGetPriSpecIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_o2_consumption_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word

  ! Reserve auxiliary data slots: aero_rate, fe_rate, hs_rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 3

  word = 'Fe++'
  this%fe2_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Fe+++'
  this%fe3_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'H+'
  this%h_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'HS-'
  this%hs_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'SO4--'
  this%so4_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'SOC(aq)'
  this%soc_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%hco3_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

end subroutine O2ConsumptionSetup

! ************************************************************************** !
subroutine O2ConsumptionAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_o2_consumption_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units

  word = 'Aerobic Respiration Rate'
  units = 'mol/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)

  word = 'Abiotic Fe2+ Oxidation Rate'
  units = 'mol/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+2)

  word = 'Abiotic HS- Oxidation Rate'
  units = 'mol/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+3)

end subroutine O2ConsumptionAuxiliaryPlotVariables

! ************************************************************************** !
subroutine O2ConsumptionEvaluate(this,Residual,Jacobian,compute_derivative, &
                             rt_auxvar,global_auxvar,material_auxvar, &
                             reaction,option)
  !
  ! Evaluates combined O2-consuming reactions
  !
  ! SOC(aq) + O2(aq) + H2O -> HCO3- + H+          (aerobic respiration)
  ! Fe2+ + 0.25 O2(aq) + H+  -> Fe3+ + 0.5 H2O    (Fe2+ oxidation)
  ! HS-  + 2.0  O2(aq)        -> SO4-- + H+         (HS- oxidation)
  !
  ! Key stability features:
  !   1. Smooth sigmoid for O2 activation with hard cutoff at 1e-20
  !   2. Coordinated O2 consumption cap across all three reactions
  !   3. No Monod on O2 (avoids stiff dRate/d[O2] near K)
  !
  ! Author: Christian Dewey
  ! Date: 2026/02/19

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  implicit none
  class(reaction_sandbox_o2_consumption_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume
  PetscReal :: porosity
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity
  PetscReal :: L_water

  PetscReal :: Fe2, O2aq, HS, SOC
  PetscReal :: f_o2, f_fe, f_hs, f_soc
  PetscReal :: rate_aero, rate_fe, rate_hs
  PetscReal :: total_o2_demand, o2_available, scale
  PetscReal :: threshold_sq

  PetscReal :: Rate_aero_mol, Rate_fe_mol, Rate_hs_mol
  PetscReal :: stoi_aero_o2, stoi_fe_o2, stoi_hs_o2

  ! Hard cutoff for sigmoid activation — below this, treat as zero
  ! to avoid denormalized float arithmetic
  PetscReal, parameter :: F_O2_CUTOFF = 1.d-20

  ! Stoichiometric coefficients for O2 consumption
  stoi_aero_o2 = 1.0d0   ! mol O2 per mol SOC(aq) respired
  stoi_fe_o2 = 0.25d0    ! mol O2 per mol Fe2+ oxidized
  stoi_hs_o2 = 2.0d0     ! mol O2 per mol HS- oxidized

  volume = material_auxvar%volume
  molality_to_molarity = global_auxvar%den_kg(iphase) * 1.d-3
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  L_water = porosity * liquid_saturation * volume * 1.d3

  ! Get species concentrations (molarity)
  Fe2 = rt_auxvar%pri_molal(this%fe2_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%fe2_id)
  O2aq = rt_auxvar%pri_molal(this%o2_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%o2_id)
  HS = rt_auxvar%pri_molal(this%hs_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%hs_id)
  SOC = rt_auxvar%pri_molal(this%soc_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%soc_id)

  ! Only compute in saturated zone
  if (liquid_saturation < 0.95d0) then
    rt_auxvar%auxiliary_data(this%auxiliary_offset+1) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+2) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+3) = 0.d0
    return
  endif

  ! ---- Smooth sigmoid O2 activation ----
  ! f_o2 = O2^2 / (O2^2 + threshold^2)
  ! Gives: ~0 when O2 << threshold, ~1 when O2 >> threshold
  ! Continuous first derivative (unlike hard if/else)
  ! Avoids Monod-style dRate/d[O2] stiffness
  threshold_sq = this%o2_threshold * this%o2_threshold
  if (O2aq > 0.d0) then
    f_o2 = (O2aq * O2aq) / (O2aq * O2aq + threshold_sq)
  else
    f_o2 = 0.d0
  endif

  ! Hard cutoff: skip all reactions if f_o2 is negligibly small
  if (f_o2 < F_O2_CUTOFF) then
    rt_auxvar%auxiliary_data(this%auxiliary_offset+1) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+2) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+3) = 0.d0
    return
  endif

  ! ---- Substrate Monod terms (SOC, Fe2+, HS- only — NOT O2) ----
  if (SOC > 0.d0) then
    f_soc = SOC / (SOC + this%aero_half_sat)
  else
    f_soc = 0.d0
  endif

  if (Fe2 > 0.d0) then
    f_fe = Fe2 / (Fe2 + this%fe_half_sat)
  else
    f_fe = 0.d0
  endif

  if (HS > 0.d0) then
    f_hs = HS / (HS + this%hs_half_sat)
  else
    f_hs = 0.d0
  endif

  ! ---- Unconstrained rates (mol/L-water/s) ----
  rate_aero = this%aero_vmax * f_o2 * f_soc
  rate_fe = this%fe_vmax * f_o2 * f_fe
  rate_hs = this%hs_vmax * f_o2 * f_hs

  ! ---- Coordinated O2 consumption cap ----
  ! Total O2 demand from all three reactions
  ! Prevent combined demand from exceeding a safe fraction of available O2
  total_o2_demand = rate_aero * stoi_aero_o2 + &
                    rate_fe * stoi_fe_o2 + &
                    rate_hs * stoi_hs_o2

  if (total_o2_demand > 0.d0 .and. option%tran_dt > 0.d0) then
    ! mol O2 / L available over this timestep, converted to rate
    o2_available = O2aq / option%tran_dt
    if (total_o2_demand > this%o2_safety_factor * o2_available) then
      scale = this%o2_safety_factor * o2_available / total_o2_demand
      rate_aero = rate_aero * scale
      rate_fe = rate_fe * scale
      rate_hs = rate_hs * scale
    endif
  endif

  ! ---- Convert to mol/sec (multiply by L_water) ----
  Rate_aero_mol = rate_aero * L_water
  Rate_fe_mol = rate_fe * L_water
  Rate_hs_mol = rate_hs * L_water

  ! Store rates for output
  rt_auxvar%auxiliary_data(this%auxiliary_offset+1) = Rate_aero_mol
  rt_auxvar%auxiliary_data(this%auxiliary_offset+2) = Rate_fe_mol
  rt_auxvar%auxiliary_data(this%auxiliary_offset+3) = Rate_hs_mol

  ! ---- Aerobic respiration residuals ----
  ! SOC(aq) + O2(aq) + H2O -> HCO3- + H+
  ! Convention: += for consumed species, -= for produced species (positive rates)
  if (Rate_aero_mol > 0.d0) then
    Residual(this%soc_id) = Residual(this%soc_id) + Rate_aero_mol
    Residual(this%o2_id) = Residual(this%o2_id) + Rate_aero_mol * stoi_aero_o2
    Residual(this%hco3_id) = Residual(this%hco3_id) - Rate_aero_mol
    Residual(this%h_id) = Residual(this%h_id) - Rate_aero_mol
  endif

  ! ---- Fe2+ oxidation residuals ----
  ! Fe2+ + 0.25 O2(aq) + H+ -> Fe3+ + 0.5 H2O
  if (Rate_fe_mol > 0.d0) then
    Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_fe_mol
    Residual(this%fe3_id) = Residual(this%fe3_id) - Rate_fe_mol
    Residual(this%o2_id) = Residual(this%o2_id) + Rate_fe_mol * stoi_fe_o2
    Residual(this%h_id) = Residual(this%h_id) + Rate_fe_mol
  endif

  ! ---- HS- oxidation residuals ----
  ! HS- + 2.0 O2(aq) -> SO4-- + H+
  if (Rate_hs_mol > 0.d0) then
    Residual(this%hs_id) = Residual(this%hs_id) + Rate_hs_mol
    Residual(this%so4_id) = Residual(this%so4_id) - Rate_hs_mol
    Residual(this%o2_id) = Residual(this%o2_id) + Rate_hs_mol * stoi_hs_o2
    Residual(this%h_id) = Residual(this%h_id) - Rate_hs_mol
  endif

end subroutine O2ConsumptionEvaluate

end module Reaction_Sandbox_O2_Consumption_class
