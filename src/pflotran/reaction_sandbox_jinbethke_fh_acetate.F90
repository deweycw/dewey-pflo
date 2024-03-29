module Reaction_Sandbox_JinBethke_Ferrihydrite_Acetate_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_jinbethke_ferrihydrite_acetate_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: fe2_id
    PetscInt :: acetate_id
    PetscInt :: bicarbonate_id
    PetscInt :: o2aq_id
    PetscInt :: mineral_id
    PetscInt :: fim_id
    PetscReal :: rmax
    PetscReal :: rate_precip
    PetscReal :: Kdonor
    PetscReal :: Y
    PetscReal :: m
    PetscReal :: chi
    PetscReal :: o2_threshold

  contains
    procedure, public :: ReadInput => JinBethkeFerrihydriteAcetateReadInput
    procedure, public :: Setup => JinBethkeFerrihydriteAcetateSetup
    procedure, public :: AuxiliaryPlotVariables => FerrihydriteAuxiliaryPlotVariables
    procedure, public :: Evaluate => JinBethkeFerrihydriteAcetateEvaluate
    procedure, public :: UpdateKineticState => JinBethkeFerrihydriteAcetateUpdateKineticState
  end type reaction_sandbox_jinbethke_ferrihydrite_acetate_type

  public :: JinBethkeFerrihydriteAcetateCreate, &
            JinBethkeFerrihydriteAcetateSetup
contains
! ************************************************************************** !
function JinBethkeFerrihydriteAcetateCreate()
  !
  ! Allocates Ferrihydrite reaction object.
  !
  implicit none
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type), pointer :: JinBethkeFerrihydriteAcetateCreate
  allocate(JinBethkeFerrihydriteAcetateCreate)
  JinBethkeFerrihydriteAcetateCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%h_ion_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%fe2_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%acetate_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%o2aq_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%mineral_id = UNINITIALIZED_INTEGER
  JinBethkeFerrihydriteAcetateCreate%fim_id = UNINITIALIZED_INTEGER

  JinBethkeFerrihydriteAcetateCreate%rmax = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%rate_precip = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%Kdonor = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%Y = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%m = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%chi = UNINITIALIZED_DOUBLE
  JinBethkeFerrihydriteAcetateCreate%o2_threshold = UNINITIALIZED_DOUBLE

  nullify(JinBethkeFerrihydriteAcetateCreate%next)
end function JinBethkeFerrihydriteAcetateCreate
! ************************************************************************** !
subroutine JinBethkeFerrihydriteAcetateReadInput(this,input,option)
  !
  ! Reads calcite reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,JINBETHKE_FERRIHYDRITE_ACETATE'
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
      case('K_DONOR')
        call InputReadDouble(input,option,this%Kdonor)
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
      case('K_PRECIPITATION')
        call InputReadDouble(input,option,this%rate_precip)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rmax) .or. &
      Uninitialized(this%rate_precip) .or. &
      Uninitialized(this%Kdonor) .or. &
      Uninitialized(this%Y) .or. &
      Uninitialized(this%m) .or. &
      Uninitialized(this%o2_threshold) .or. &
      Uninitialized(this%chi)) then
    option%io_buffer = 'RMAX, K_PRECIPITATION, K_DONOR, Y, M, CHI, and O2_THRESHOLD must be set for &
      JINBETHKE_FERRIHYDRITE_ACETATE.'
    call PrintErrMsg(option)
  endif
end subroutine JinBethkeFerrihydriteAcetateReadInput
! ************************************************************************** !
subroutine JinBethkeFerrihydriteAcetateSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Reaction_Immobile_Aux_module, only: GetImmobileSpeciesIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type) :: this
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
  word = 'Fe++'
  this%fe2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Ac-'
  this%acetate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%bicarbonate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Ferrihydrite'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)
  word = 'Fim'
  this%fim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

end subroutine JinBethkeFerrihydriteAcetateSetup
! ************************************************************************** !
subroutine FerrihydriteAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds ferrihydrite auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'JB Ferrihydrite Reduction / Acetate Oxidation Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine FerrihydriteAuxiliaryPlotVariables
! ************************************************************************** !
subroutine JinBethkeFerrihydriteAcetateEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !Jacobian,compute_derivative,
  ! Evaluates ferrihydrite reaction storing residual but no Jacobian
  !
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2022/10/4
  ! Modified 2023/2/7

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Mineral_Aux_module
  implicit none
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type) :: this
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

  PetscReal :: Ac, Proton, Fe2, Bicarbonate
  PetscReal :: fim, yield, O2aq
  PetscReal :: Rate, Rate_Ac, Rate_Proton, Rate_fh
  PetscReal :: Rate_Fe2, Rate_Bicarbonate, Rate_O2aq
  PetscReal :: stoi_ac, stoi_proton
  PetscReal :: stoi_fe2, stoi_bicarbonate
  PetscReal :: k_diss, k_precip, m, chi
  PetscReal :: temp_K, RT
  PetscReal :: Ft, Ftr, Fa
  PetscReal :: reaction_Q
  PetscReal :: dG0, dGr, dG_ATP

  PetscReal :: lnQK, sign_, affinity_factor, QK

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: imnrl
  PetscInt :: iauxiliary
  PetscBool :: calculate_dissolution
  PetscBool :: calculate_precip
  mineral => reaction%mineral

  iauxiliary = this%auxiliary_offset + 1

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  imnrl = this%mineral_id

  if (dabs(rt_auxvar%mnrl_rate(imnrl)) > 1.d-40) then
    option%io_buffer = 'For JINBETHKE_FERRIHYDRITE_ACETATE to function correctly, &
      &the RATE_CONSTANT in the default MINERAL_KINETICS block must be set &
      &to zero.'
    call PrintErrMsg(option)
  endif

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  k_diss = this%rmax
  k_precip = this%rate_precip


  ! Rxn:    1.00 Ac- + 8.00 FHY + 15.00 H+ = 8.00 Fe++ + 2.00 HCO3- + 20.00 H2O 
  ! dG0 =   -612.0 kJ per mol Ac-, Kocar & Fendorf 2009

  Ac = rt_auxvar%pri_molal(this%acetate_id) * &
    rt_auxvar%pri_act_coef(this%acetate_id) 
  Proton = rt_auxvar%pri_molal(this%h_ion_id) * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  Fe2 = rt_auxvar%pri_molal(this%fe2_id) * &
    rt_auxvar%pri_act_coef(this%fe2_id) 
  Bicarbonate = rt_auxvar%pri_molal(this%bicarbonate_id) * &
    rt_auxvar%pri_act_coef(this%bicarbonate_id) 
  O2aq = rt_auxvar%pri_molal(this%o2aq_id) * &
    rt_auxvar%pri_act_coef(this%o2aq_id) 

  fim = rt_auxvar%immobile(this%fim_id)

  m = this%m
  chi = this%chi

  stoi_fe2 = 8.d0
  stoi_bicarbonate = 2.d0
  
  stoi_ac = 1.d0
  stoi_proton = 15.d0 !+ 2.d0  ! +2.d0 to account for H+ consumed in database formulation 

  RT = (8.314e-3) * (global_auxvar%temp + 273.15d0)
  dG0 = (-612.0d0) ! kJ / mol acetate; dG0 for FeIII in ferrihydrite as electron acceptor
  dG_ATP = 50.d0 ! kJ / mol ATP

  reaction_Q = ( (Fe2**stoi_fe2) * (Bicarbonate**stoi_bicarbonate)) / &
    ((Ac**stoi_ac) * (Proton**stoi_proton))

  dGr = dG0 + RT*log(reaction_Q)
  
  ! Monod expressions for acetate
  Fa = Ac / (Ac + this%Kdonor)

  ! Thermodynamic factor 
  Ft = 1.d0 - exp((dGr + m*dG_ATP) / (chi * RT))

  if (Ft < 0) then 
    Ftr = 0.d0
  else
    Ftr = Ft
  endif

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

  ! only calculate diss / iron reduciton rate if mineral is present and O2(aq) below threhsold
  calculate_dissolution = (rt_auxvar%mnrl_volfrac(imnrl) > 0 .and. &
    O2aq < (this%o2_threshold))

  Rate = 0.d0
  
  if (calculate_dissolution) then
    ! base rate, mol/sec/m^3 bulk
    ! units on k: mol/sec/mol-bio

    Rate = -k_diss *  Fa * Ftr * fim  

    Rate_fh = Rate
    
    rt_auxvar%auxiliary_data(iauxiliary) = Rate_fh

    Rate = Rate * material_auxvar%volume ! mol/sec
      
    ! species-specifc 
    Rate_Ac = Rate * stoi_ac  
    Rate_Proton = Rate * stoi_proton 
    Rate_Fe2 = Rate * stoi_fe2 
    Rate_Bicarbonate = Rate * stoi_bicarbonate 
    !Rate_fim = Rate * yield
    
    Residual(this%h_ion_id) = Residual(this%h_ion_id) - Rate_Proton
    Residual(this%acetate_id) = Residual(this%acetate_id) - Rate_Ac
    Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_Fe2
    Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + Rate_Bicarbonate

  else

    if (calculate_precip) then
      
      Rate = (-1.d0) * sign_ * abs(affinity_factor) * this%rate_precip

      !multiple Rate by 8 for Fe stoichiometry?
      rt_auxvar%auxiliary_data(iauxiliary) = Rate

      Rate = Rate * material_auxvar%volume ! mol/sec
        
      ! species-specifc 
      Rate_O2aq = Rate * (0.0625d0) 
      Rate_Fe2 = Rate 
      Rate_Proton = Rate * (2.d0) 
      !Rate_fim = Rate * yield
      
      Residual(this%h_ion_id) = Residual(this%h_ion_id) - Rate_Proton
      Residual(this%o2aq_id) = Residual(this%o2aq_id) + Rate_O2aq
      Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_Fe2
    endif

  endif

end subroutine JinBethkeFerrihydriteAcetateEvaluate
! ************************************************************************** !
subroutine JinBethkeFerrihydriteAcetateUpdateKineticState(this,rt_auxvar,global_auxvar, &
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
  class(reaction_sandbox_jinbethke_ferrihydrite_acetate_type) :: this
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
end subroutine JinBethkeFerrihydriteAcetateUpdateKineticState
end module Reaction_Sandbox_JinBethke_Ferrihydrite_Acetate_class