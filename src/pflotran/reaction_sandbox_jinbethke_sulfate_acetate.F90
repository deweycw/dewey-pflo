module Reaction_Sandbox_JinBethke_Sulfate_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_jinbethke_sulfate_type
    PetscInt :: auxiliary_offset
    PetscInt :: so4_id
    PetscInt :: acetate_id
    PetscInt :: bicarbonate_id
    PetscInt :: hs_id
    PetscInt :: sim_id
    PetscInt :: o2aq_id
    PetscReal :: rmax
    PetscReal :: Kdonor
    PetscReal :: Kacceptor
    PetscReal :: Y
    PetscReal :: m
    PetscReal :: chi
    PetscReal :: o2_threshold

  contains
    procedure, public :: ReadInput => JinBethkeSulfateReadInput
    procedure, public :: Setup => JinBethkeSulfateSetup
    procedure, public :: AuxiliaryPlotVariables => JinBethkeSulfateAuxiliaryPlotVariables
    procedure, public :: Evaluate => JinBethkeSulfateEvaluate
  end type reaction_sandbox_jinbethke_sulfate_type

  public :: JinBethkeSulfateCreate, &
            JinBethkeSulfateSetup
contains
! ************************************************************************** !
function JinBethkeSulfateCreate()
  !
  ! Allocates JinBethkeSulfate reaction object.
  !
  implicit none
  class(reaction_sandbox_jinbethke_sulfate_type), pointer :: JinBethkeSulfateCreate
  allocate(JinBethkeSulfateCreate)
  JinBethkeSulfateCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  ! reactants
  JinBethkeSulfateCreate%so4_id = UNINITIALIZED_INTEGER
  JinBethkeSulfateCreate%acetate_id = UNINITIALIZED_INTEGER
  ! products
  JinBethkeSulfateCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  JinBethkeSulfateCreate%hs_id = UNINITIALIZED_INTEGER

  JinBethkeSulfateCreate%o2aq_id = UNINITIALIZED_INTEGER
  JinBethkeSulfateCreate%sim_id = UNINITIALIZED_INTEGER

  JinBethkeSulfateCreate%rmax = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%Kdonor = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%Kacceptor = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%Y = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%m = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%chi = UNINITIALIZED_DOUBLE
  JinBethkeSulfateCreate%o2_threshold = UNINITIALIZED_DOUBLE

  nullify(JinBethkeSulfateCreate%next)
end function JinBethkeSulfateCreate
! ************************************************************************** !
subroutine JinBethkeSulfateReadInput(this,input,option)
  !
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_jinbethke_sulfate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,JINBETHKE_SULFATE_ACETATE'
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
      Uninitialized(this%m) .or. &
      Uninitialized(this%chi) .or. &
      Uninitialized(this%Y)) then
    option%io_buffer = 'RMAX, Kdonor, Kacceptor, m, chi, and Y must be set for &
      JINBETHKE_SULFATE_ACETATE'
    call PrintErrMsg(option)
  endif
end subroutine JinBethkeSulfateReadInput
! ************************************************************************** !
subroutine JinBethkeSulfateSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only: GetImmobileSpeciesIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_jinbethke_sulfate_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1
  ! Aqueous species
  word = 'HS-'
  this%hs_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'SO4--'
  this%so4_id = &
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
  word = 'Sim'
  this%sim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

end subroutine JinBethkeSulfateSetup
! ************************************************************************** !
subroutine JinBethkeSulfateAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds sulfate reduction auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_jinbethke_sulfate_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'JB Sulfafte Acetate Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine JinBethkeSulfateAuxiliaryPlotVariables
! ************************************************************************** !
subroutine JinBethkeSulfateEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  ! 
  ! Evaluates sulfate reduction reaction storing residual but no Jacobian
  ! Based on Jin & Bethke (2003, 2005)
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/2/7
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  implicit none
  class(reaction_sandbox_jinbethke_sulfate_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt, parameter :: iphase = 1
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: porosity             ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: L_water              ! L water

  PetscReal :: Ac, hs, so4, Bicarbonate
  PetscReal :: Sim, yield, O2aq
  PetscReal :: Rate, Rate_Ac, Rate_hs, Rate_so4, Rate_sulf
  PetscReal :: Rate_Bicarbonate
  PetscReal :: stoi_ac, stoi_so4
  PetscReal :: stoi_hs, stoi_bicarbonate
  PetscReal :: k_rmax, m, chi
  PetscReal :: temp_K, RT
  PetscReal :: Ft, Ftr, Fdonor, Facceptor
  PetscReal :: reaction_Q
  PetscReal :: dG0, dGr, dG_ATP

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
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

  k_rmax = this%rmax

  ! Rxn:    1.00 Ac- + 1.00 SO4-- = 1.00 HS- + 2.00 HCO3- 
  ! dG0 =   -48.1 kJ per mol Ac-, Kocar & Fendorf 2009

  Ac = rt_auxvar%pri_molal(this%acetate_id) * &
    rt_auxvar%pri_act_coef(this%acetate_id) 
  so4 = rt_auxvar%pri_molal(this%so4_id) * &
    rt_auxvar%pri_act_coef(this%so4_id) 
  
  hs = rt_auxvar%pri_molal(this%hs_id) * &
    rt_auxvar%pri_act_coef(this%hs_id) 
  Bicarbonate = rt_auxvar%pri_molal(this%bicarbonate_id) * &
    rt_auxvar%pri_act_coef(this%bicarbonate_id) 

  O2aq = rt_auxvar%pri_molal(this%o2aq_id) * &
    rt_auxvar%pri_act_coef(this%o2aq_id) 

  Sim = rt_auxvar%immobile(this%sim_id)

  m = this%m
  chi = this%chi

  stoi_ac = 1.d0
  stoi_so4 = 1.d0

  stoi_bicarbonate = 2.d0
  stoi_hs = 1.d0 

  RT = (8.314e-3) * (global_auxvar%temp + 273.15d0)
  dG0 = (-47.64d0) ! kJ / mol acetate; dG0 for SO4-- as electron acceptor
  dG_ATP = 50.d0 ! kJ / mol ATP

  reaction_Q = ( (hs**stoi_hs) * (Bicarbonate**stoi_bicarbonate)) / &
    ((Ac**stoi_ac) * (so4**stoi_so4))

  dGr = dG0 + RT*log(reaction_Q)
  
  ! Monod expressions for acetate
  Fdonor = Ac / (Ac + this%Kdonor)

  ! Monod expressions for sulfate
  Facceptor = so4 / (so4 + this%Kacceptor)

  ! Thermodynamic factor 
  Ft = 1.d0 - exp((dGr + m*dG_ATP) / (chi * RT))

  if (Ft <= 0.d0) then 
    Ftr = 0.d0
  else
    Ftr = Ft
  endif

  ! only sulfate reduction if O2(aq) below threhsold
  calculate_rate = O2aq < (this%o2_threshold)

  Rate = 0.d0
  
  if (calculate_rate) then
    ! base rate, mol/sec/m^3 bulk
    ! units on k: mol/sec/mol-bio

    Rate = -k_rmax *  Facceptor * Fdonor * Ftr * Sim  
 
    Rate_sulf = Rate
    
    rt_auxvar%auxiliary_data(iauxiliary) = Rate_sulf

    Rate = Rate * material_auxvar%volume ! mol/sec
      
    ! species-specifc 
    Rate_Ac = Rate * stoi_ac  
    Rate_so4 = Rate * stoi_so4 
    Rate_hs = Rate * stoi_hs 
    Rate_Bicarbonate = Rate * stoi_bicarbonate 
    !Rate_Nim = Rate * yield
    
    Residual(this%so4_id) = Residual(this%so4_id) - Rate_so4
    Residual(this%acetate_id) = Residual(this%acetate_id) - Rate_Ac
    Residual(this%hs_id) = Residual(this%hs_id) + Rate_hs
    Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + Rate_Bicarbonate

  endif
 
end subroutine JinBethkeSulfateEvaluate

end module Reaction_Sandbox_JinBethke_Sulfate_class