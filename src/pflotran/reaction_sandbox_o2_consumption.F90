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
!   - Smooth sigmoid O2 activation: f_o2 = O2^2/(O2^2 + K^2)
!   - Substrate limitation via Monod for SOC, Fe2+, and HS- only
!   - Analytical Jacobian provided for Newton convergence
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
  ! Evaluates combined O2-consuming reactions with analytical Jacobian
  !
  ! SOC(aq) + O2(aq) + H2O -> HCO3- + H+          (aerobic respiration)
  ! Fe2+ + 0.25 O2(aq) + H+  -> Fe3+ + 0.5 H2O    (Fe2+ oxidation)
  ! HS-  + 2.0  O2(aq)        -> SO4-- + H+         (HS- oxidation)
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
  PetscReal :: threshold_sq, o2_sq, denom_o2

  PetscReal :: Rate_aero_mol, Rate_fe_mol, Rate_hs_mol
  PetscReal :: stoi_aero_o2, stoi_fe_o2, stoi_hs_o2

  ! Jacobian variables
  PetscReal :: df_o2_dO2, df_soc_dSOC, df_fe_dFe2, df_hs_dHS
  PetscReal :: dr_aero_dm_o2, dr_aero_dm_soc
  PetscReal :: dr_fe_dm_o2, dr_fe_dm_fe2
  PetscReal :: dr_hs_dm_o2, dr_hs_dm_hs
  PetscReal :: cf_o2, cf_fe2, cf_soc, cf_hs

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

  ! Get species concentrations (activity in molarity)
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
  ! f_o2 = O2^2 / (O2^2 + K^2)
  threshold_sq = this%o2_threshold * this%o2_threshold
  if (O2aq > 0.d0) then
    o2_sq = O2aq * O2aq
    denom_o2 = o2_sq + threshold_sq
    f_o2 = o2_sq / denom_o2
  else
    f_o2 = 0.d0
  endif

  if (f_o2 < F_O2_CUTOFF) then
    rt_auxvar%auxiliary_data(this%auxiliary_offset+1) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+2) = 0.d0
    rt_auxvar%auxiliary_data(this%auxiliary_offset+3) = 0.d0
    return
  endif

  ! ---- Substrate Monod terms ----
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

  ! ---- Rates (mol/L-water/s) ----
  rate_aero = this%aero_vmax * f_o2 * f_soc
  rate_fe = this%fe_vmax * f_o2 * f_fe
  rate_hs = this%hs_vmax * f_o2 * f_hs

  ! ---- Convert to mol/sec ----
  Rate_aero_mol = rate_aero * L_water
  Rate_fe_mol = rate_fe * L_water
  Rate_hs_mol = rate_hs * L_water

  ! Store rates for output
  rt_auxvar%auxiliary_data(this%auxiliary_offset+1) = Rate_aero_mol
  rt_auxvar%auxiliary_data(this%auxiliary_offset+2) = Rate_fe_mol
  rt_auxvar%auxiliary_data(this%auxiliary_offset+3) = Rate_hs_mol

  ! ---- Residuals ----
  ! Convention: += consumed, -= produced (positive rates)

  ! Aerobic: SOC(aq) + O2(aq) + H2O -> HCO3- + H+
  Residual(this%soc_id) = Residual(this%soc_id) + Rate_aero_mol
  Residual(this%o2_id) = Residual(this%o2_id) + Rate_aero_mol * stoi_aero_o2
  Residual(this%hco3_id) = Residual(this%hco3_id) - Rate_aero_mol
  Residual(this%h_id) = Residual(this%h_id) - Rate_aero_mol

  ! Fe2+ oxidation: Fe2+ + 0.25 O2 + H+ -> Fe3+ + 0.5 H2O
  Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_fe_mol
  Residual(this%fe3_id) = Residual(this%fe3_id) - Rate_fe_mol
  Residual(this%o2_id) = Residual(this%o2_id) + Rate_fe_mol * stoi_fe_o2
  Residual(this%h_id) = Residual(this%h_id) + Rate_fe_mol

  ! HS- oxidation: HS- + 2 O2 -> SO4-- + H+
  Residual(this%hs_id) = Residual(this%hs_id) + Rate_hs_mol
  Residual(this%so4_id) = Residual(this%so4_id) - Rate_hs_mol
  Residual(this%o2_id) = Residual(this%o2_id) + Rate_hs_mol * stoi_hs_o2
  Residual(this%h_id) = Residual(this%h_id) - Rate_hs_mol

  ! ---- Analytical Jacobian ----
  if (compute_derivative) then

    ! d(conc_j)/d(pri_molal_j) conversion factors
    cf_o2 = molality_to_molarity * rt_auxvar%pri_act_coef(this%o2_id)
    cf_soc = molality_to_molarity * rt_auxvar%pri_act_coef(this%soc_id)
    cf_fe2 = molality_to_molarity * rt_auxvar%pri_act_coef(this%fe2_id)
    cf_hs = molality_to_molarity * rt_auxvar%pri_act_coef(this%hs_id)

    ! df_o2/dO2 = 2*O2*K^2 / (O2^2 + K^2)^2
    if (O2aq > 0.d0) then
      df_o2_dO2 = 2.d0 * O2aq * threshold_sq / (denom_o2 * denom_o2)
    else
      df_o2_dO2 = 0.d0
    endif

    ! df_soc/dSOC = K_soc / (SOC + K_soc)^2
    if (SOC > 0.d0) then
      df_soc_dSOC = this%aero_half_sat / &
        ((SOC + this%aero_half_sat) * (SOC + this%aero_half_sat))
    else
      df_soc_dSOC = 0.d0
    endif

    ! df_fe/dFe2 = K_fe / (Fe2 + K_fe)^2
    if (Fe2 > 0.d0) then
      df_fe_dFe2 = this%fe_half_sat / &
        ((Fe2 + this%fe_half_sat) * (Fe2 + this%fe_half_sat))
    else
      df_fe_dFe2 = 0.d0
    endif

    ! df_hs/dHS = K_hs / (HS + K_hs)^2
    if (HS > 0.d0) then
      df_hs_dHS = this%hs_half_sat / &
        ((HS + this%hs_half_sat) * (HS + this%hs_half_sat))
    else
      df_hs_dHS = 0.d0
    endif

    ! Rate derivatives w.r.t. pri_molal (mol/L/s per molal)
    ! d(rate)/d(molal_j) = d(rate)/d(conc_j) * cf_j

    ! Aerobic: rate_aero = aero_vmax * f_o2 * f_soc
    dr_aero_dm_o2 = this%aero_vmax * df_o2_dO2 * f_soc * cf_o2
    dr_aero_dm_soc = this%aero_vmax * f_o2 * df_soc_dSOC * cf_soc

    ! Fe2+: rate_fe = fe_vmax * f_o2 * f_fe
    dr_fe_dm_o2 = this%fe_vmax * df_o2_dO2 * f_fe * cf_o2
    dr_fe_dm_fe2 = this%fe_vmax * f_o2 * df_fe_dFe2 * cf_fe2

    ! HS-: rate_hs = hs_vmax * f_o2 * f_hs
    dr_hs_dm_o2 = this%hs_vmax * df_o2_dO2 * f_hs * cf_o2
    dr_hs_dm_hs = this%hs_vmax * f_o2 * df_hs_dHS * cf_hs

    ! ---- Aerobic respiration Jacobian ----
    ! R(soc) += rate_aero * L_water
    Jacobian(this%soc_id,this%o2_id) = &
      Jacobian(this%soc_id,this%o2_id) + dr_aero_dm_o2 * L_water
    Jacobian(this%soc_id,this%soc_id) = &
      Jacobian(this%soc_id,this%soc_id) + dr_aero_dm_soc * L_water
    ! R(o2) += rate_aero * stoi_aero_o2 * L_water
    Jacobian(this%o2_id,this%o2_id) = &
      Jacobian(this%o2_id,this%o2_id) + dr_aero_dm_o2 * stoi_aero_o2 * L_water
    Jacobian(this%o2_id,this%soc_id) = &
      Jacobian(this%o2_id,this%soc_id) + dr_aero_dm_soc * stoi_aero_o2 * L_water
    ! R(hco3) -= rate_aero * L_water
    Jacobian(this%hco3_id,this%o2_id) = &
      Jacobian(this%hco3_id,this%o2_id) - dr_aero_dm_o2 * L_water
    Jacobian(this%hco3_id,this%soc_id) = &
      Jacobian(this%hco3_id,this%soc_id) - dr_aero_dm_soc * L_water
    ! R(h) -= rate_aero * L_water
    Jacobian(this%h_id,this%o2_id) = &
      Jacobian(this%h_id,this%o2_id) - dr_aero_dm_o2 * L_water
    Jacobian(this%h_id,this%soc_id) = &
      Jacobian(this%h_id,this%soc_id) - dr_aero_dm_soc * L_water

    ! ---- Fe2+ oxidation Jacobian ----
    ! R(fe2) += rate_fe * L_water
    Jacobian(this%fe2_id,this%o2_id) = &
      Jacobian(this%fe2_id,this%o2_id) + dr_fe_dm_o2 * L_water
    Jacobian(this%fe2_id,this%fe2_id) = &
      Jacobian(this%fe2_id,this%fe2_id) + dr_fe_dm_fe2 * L_water
    ! R(fe3) -= rate_fe * L_water
    Jacobian(this%fe3_id,this%o2_id) = &
      Jacobian(this%fe3_id,this%o2_id) - dr_fe_dm_o2 * L_water
    Jacobian(this%fe3_id,this%fe2_id) = &
      Jacobian(this%fe3_id,this%fe2_id) - dr_fe_dm_fe2 * L_water
    ! R(o2) += rate_fe * stoi_fe_o2 * L_water
    Jacobian(this%o2_id,this%o2_id) = &
      Jacobian(this%o2_id,this%o2_id) + dr_fe_dm_o2 * stoi_fe_o2 * L_water
    Jacobian(this%o2_id,this%fe2_id) = &
      Jacobian(this%o2_id,this%fe2_id) + dr_fe_dm_fe2 * stoi_fe_o2 * L_water
    ! R(h) += rate_fe * L_water
    Jacobian(this%h_id,this%o2_id) = &
      Jacobian(this%h_id,this%o2_id) + dr_fe_dm_o2 * L_water
    Jacobian(this%h_id,this%fe2_id) = &
      Jacobian(this%h_id,this%fe2_id) + dr_fe_dm_fe2 * L_water

    ! ---- HS- oxidation Jacobian ----
    ! R(hs) += rate_hs * L_water
    Jacobian(this%hs_id,this%o2_id) = &
      Jacobian(this%hs_id,this%o2_id) + dr_hs_dm_o2 * L_water
    Jacobian(this%hs_id,this%hs_id) = &
      Jacobian(this%hs_id,this%hs_id) + dr_hs_dm_hs * L_water
    ! R(so4) -= rate_hs * L_water
    Jacobian(this%so4_id,this%o2_id) = &
      Jacobian(this%so4_id,this%o2_id) - dr_hs_dm_o2 * L_water
    Jacobian(this%so4_id,this%hs_id) = &
      Jacobian(this%so4_id,this%hs_id) - dr_hs_dm_hs * L_water
    ! R(o2) += rate_hs * stoi_hs_o2 * L_water
    Jacobian(this%o2_id,this%o2_id) = &
      Jacobian(this%o2_id,this%o2_id) + dr_hs_dm_o2 * stoi_hs_o2 * L_water
    Jacobian(this%o2_id,this%hs_id) = &
      Jacobian(this%o2_id,this%hs_id) + dr_hs_dm_hs * stoi_hs_o2 * L_water
    ! R(h) -= rate_hs * L_water
    Jacobian(this%h_id,this%o2_id) = &
      Jacobian(this%h_id,this%o2_id) - dr_hs_dm_o2 * L_water
    Jacobian(this%h_id,this%hs_id) = &
      Jacobian(this%h_id,this%hs_id) - dr_hs_dm_hs * L_water

  endif ! compute_derivative

end subroutine O2ConsumptionEvaluate

end module Reaction_Sandbox_O2_Consumption_class
