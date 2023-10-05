module Reaction_Sandbox_Mackinawite_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_mackinawite_type
    PetscInt :: auxiliary_offset
    PetscInt :: fe2_id
    PetscInt :: o2_id
    PetscInt :: sulfate_id
    PetscInt :: proton_id
    PetscInt :: hs_id
    PetscInt :: mineral_id
    PetscReal :: precip_ksp
    PetscReal :: precip_rate
    PetscReal :: diss_ksp
    PetscReal :: diss_rate
    PetscReal :: k_o2aq
    PetscReal :: o2_threshold

  contains
    procedure, public :: ReadInput => MackinawiteReadInput
    procedure, public :: Setup => MackinawiteSetup
    procedure, public :: AuxiliaryPlotVariables => MackinawiteAuxiliaryPlotVariables
    procedure, public :: Evaluate => MackinawiteEvaluate
    procedure, public :: UpdateKineticState => MackinawiteUpdateKineticState
  end type reaction_sandbox_mackinawite_type

  public :: MackinawiteCreate, &
            MackinawiteSetup
contains
! ************************************************************************** !
function MackinawiteCreate()
  !
  !
  implicit none
  class(reaction_sandbox_mackinawite_type), pointer :: MackinawiteCreate
  allocate(MackinawiteCreate)
  
  MackinawiteCreate%auxiliary_offset = UNINITIALIZED_INTEGER

  MackinawiteCreate%fe2_id = UNINITIALIZED_INTEGER
  MackinawiteCreate%sulfate_id = UNINITIALIZED_INTEGER
  MackinawiteCreate%o2_id = UNINITIALIZED_INTEGER
  MackinawiteCreate%proton_id = UNINITIALIZED_INTEGER
  MackinawiteCreate%hs_id = UNINITIALIZED_INTEGER
  
  MackinawiteCreate%mineral_id = UNINITIALIZED_INTEGER

  MackinawiteCreate%precip_rate = UNINITIALIZED_DOUBLE
  MackinawiteCreate%precip_ksp = UNINITIALIZED_DOUBLE
  MackinawiteCreate%diss_rate = UNINITIALIZED_DOUBLE
  MackinawiteCreate%diss_ksp = UNINITIALIZED_DOUBLE

  MackinawiteCreate%k_o2aq = UNINITIALIZED_DOUBLE
  MackinawiteCreate%o2_threshold = UNINITIALIZED_DOUBLE
  nullify(MackinawiteCreate%next)
end function MackinawiteCreate
! ************************************************************************** !
subroutine MackinawiteReadInput(this,input,option)
  !
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_mackinawite_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,MACKINAWITE'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('PRECIP_KSP')
        call InputReadDouble(input,option,this%precip_ksp)
        call InputErrorMsg(input,option,word,error_string)
      case('PRECIP_RATE')
        call InputReadDouble(input,option,this%precip_rate)
        call InputErrorMsg(input,option,word,error_string)
      case('DISS_KSP')
        call InputReadDouble(input,option,this%diss_ksp)
        call InputErrorMsg(input,option,word,error_string)
      case('DISS_RATE')
        call InputReadDouble(input,option,this%diss_rate)
        call InputErrorMsg(input,option,word,error_string)
      case('O2_THRESHOLD_M')
        call InputReadDouble(input,option,this%o2_threshold)
        call InputErrorMsg(input,option,word,error_string)
      case('HALF_SAT_O2')
        call InputReadDouble(input,option,this%k_o2aq)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%precip_rate) .or. &
      Uninitialized(this%precip_ksp)) .or. &
      Uninitialized(this%diss_rate) .or. &
      Uninitialized(this%diss_ksp)) .or. &
      Uninitialized(this%o2_threshold) .or. &
      Uninitialized(this%k_o2aq)) then
    option%io_buffer = 'RATES, KSPS, THRESHOLDS, O2_HALFSAT must be set for &
      MACKINAWITE DISSOLUTION/PRECIP.'
    call PrintErrMsg(option)
  endif
end subroutine MackinawiteReadInput
! ************************************************************************** !
subroutine MackinawiteSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_mackinawite_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1
  ! Aqueous species
  word = 'Fe++'
  this%fe2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'SO4--'
  this%sulfate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'HS-'
  this%hs_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H+'
  this%proton_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Mackinawite'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine MackinawiteSetup
! ************************************************************************** !
subroutine MackinawiteAuxiliaryPlotVariables(this,list,reaction,option)
  !
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_mackinawite_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'Mackinawite Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine MackinawiteAuxiliaryPlotVariables
! ************************************************************************** !
subroutine MackinawiteEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !Jacobian,compute_derivative,
  !
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/10/5

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Mineral_Aux_module
  implicit none
  class(reaction_sandbox_mackinawite_type) :: this
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
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: L_water              ! L water

  PetscReal :: Fe2, O2aq, sulfate, HS, Proton, K_O2aq
  PetscReal :: Rate, Rate_Fe, Rate_O2aq, Rate_Sulfate
  PetscReal :: Rate_HS, Rate_Proton
  PetscReal :: stoi_o2, stoi_fe2, stoi_sulfate
  PetscReal :: stoi_proton, stoi_hs
  PetscReal :: diss_rate_from_user, precip_rate_from_user
  PetscReal :: diss_ksp, precip_ksp
  PetscReal :: diss_Q, precip_Q
  PetscReal :: sat_index_diss, sat_index_precip

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: imnrl
  PetscInt :: iauxiliary
  mineral => reaction%mineral

  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  imnrl = this%mineral_id

  if (dabs(rt_auxvar%mnrl_rate(imnrl)) > 1.d-40) then
    option%io_buffer = 'For MACKINAWITE DISS/PRECIP to function correctly, &
      &the MACKINAWITE RATE_CONSTANT in the default MINERAL_KINETICS block must be set &
      &to zero.'
    call PrintErrMsg(option)
  endif

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  precip_rate_from_user = this%precip_rate
  diss_rate_from_user = this%diss_rate
  precip_ksp = this%precip_ksp
  diss_ksp = this%diss_ksp
  K_O2aq = this%k_o2aq
  threshold = this%o2_threshold

  Fe2 = rt_auxvar%pri_molal(this%fe2_id) * &
    rt_auxvar%pri_act_coef(this%fe2_id) 
  O2aq = rt_auxvar%pri_molal(this%o2_id) * &
    rt_auxvar%pri_act_coef(this%o2_id) 
  sulfate = rt_auxvar%pri_molal(this%sulfate_id) * &
    rt_auxvar%pri_act_coef(this%sulfate_id) 
  HS = rt_auxvar%pri_molal(this%hs_id) * &
    rt_auxvar%pri_act_coef(this%hs_id) 
  Proton = rt_auxvar%pri_molal(this%proton_id) * &
    rt_auxvar%pri_act_coef(this%proton_id) 

  stoi_fe2 = 1.d0
  stoi_sulfate = 1.d0
  stoi_o2 = 2.d0
  stoi_proton = 1.d0
  stoi_hs = 1.d0

  diss_Q = ((sulfate**stoi_sulfate) * (Fe2**stoi_fe2) ) / &
    (O2aq**stoi_o2)
  sat_index_diss = log(diss_Q / diss_ksp)

  precip_Q = ((HS ** stoi_hs) * (Fe2 ** stoi_fe2))/(Proton ** stoi_proton)
  sat_index_precip = log(precip_Q / precip_ksp)

  Rate = 0.d0 

  ! base rate, mol/sec/m^3 bulk
  ! units on k: mol/sec/mol-bio
  
  if ((O2aq <= threshold) .and. (sat_index_diss < 0)) then
    Rate = (-1.d0) * diss_rate_from_user * (O2aq  / (K_O2aq + O2aq)) 
    ! negative for dissolution 
    rt_auxvar%auxiliary_data(iauxiliary) = Rate 

    Rate = Rate * material_auxvar%volume ! mol/sec
      
    Rate_Fe = Rate * stoi_fe2  
    Rate_Sulfate = Rate * stoi_sulfate
    Rate_O2aq = Rate * stoi_o2

    Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_Fe
    Residual(this%o2_id) = Residual(this%o2_id) - Rate_O2aq
    Residual(this%sulfate_id) = Residual(this%sulfate_id) + Rate_Sulfate

  else if ((O2aq > threshold) .and. (sat_index_precip > 0))  then
    Rate = precip_rate_from_user * sat_index_precip
    ! positive for precip 
    rt_auxvar%auxiliary_data(iauxiliary) = Rate 

    Rate = Rate * material_auxvar%volume ! mol/sec
      
    Rate_Fe = Rate * stoi_fe2  
    Rate_HS = Rate * stoi_hs
    Rate_Proton = Rate * stoi_proton

    Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_Fe
    Residual(this%proton_id) = Residual(this%proton_id) - Rate_Proton
    Residual(this%hs_id) = Residual(this%hs_id) + Rate_HS

  endif
end subroutine MackinawiteEvaluate
! ************************************************************************** !
subroutine MackinawiteUpdateKineticState(this,rt_auxvar,global_auxvar, &
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
  class(reaction_sandbox_mackinawite_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
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
end subroutine MackinawiteUpdateKineticState
end module Reaction_Sandbox_Mackinawite_class