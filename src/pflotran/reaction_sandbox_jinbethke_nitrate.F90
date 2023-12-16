module Reaction_Sandbox_JinBethke_Nitrate_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_jinbethke_nitrate_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: no3_id
    PetscInt :: acetate_id
    PetscInt :: bicarbonate_id
    PetscInt :: n2aq_id
    PetscInt :: nim_id
    PetscInt :: o2aq_id
    PetscReal :: rmax
    PetscReal :: Kdonor
    PetscReal :: Kacceptor
    PetscReal :: Y
    PetscReal :: m
    PetscReal :: chi
    PetscReal :: o2_threshold

  contains
    procedure, public :: ReadInput => JinBethkeNitrateReadInput
    procedure, public :: Setup => JinBethkeNitrateSetup
    procedure, public :: AuxiliaryPlotVariables => JinBethkeNitrateAuxiliaryPlotVariables
    procedure, public :: Evaluate => JinBethkeNitrateEvaluate
  end type reaction_sandbox_jinbethke_nitrate_type

  public :: JinBethkeNitrateCreate, &
            JinBethkeNitrateSetup
contains
! ************************************************************************** !
function JinBethkeNitrateCreate()
  !
  ! Allocates JinBethkeNitrate reaction object.
  !
  implicit none
  class(reaction_sandbox_jinbethke_nitrate_type), pointer :: JinBethkeNitrateCreate
  allocate(JinBethkeNitrateCreate)
  JinBethkeNitrateCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%h_ion_id = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%no3_id = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%acetate_id = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%o2aq_id = UNINITIALIZED_INTEGER
  JinBethkeNitrateCreate%nim_id = UNINITIALIZED_INTEGER

  JinBethkeNitrateCreate%rmax = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%Kdonor = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%Kacceptor = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%Y = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%m = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%chi = UNINITIALIZED_DOUBLE
  JinBethkeNitrateCreate%o2_threshold = UNINITIALIZED_DOUBLE

  nullify(JinBethkeNitrateCreate%next)
end function JinBethkeNitrateCreate
! ************************************************************************** !
subroutine JinBethkeNitrateReadInput(this,input,option)
  !
  ! Reads calcite reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_jinbethke_nitrate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,JINBETHKE_NITRATE_ACETATE'
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
      JINBETHKE_NITRATE_ACETATE'
    call PrintErrMsg(option)
  endif
end subroutine JinBethkeNitrateReadInput
! ************************************************************************** !
subroutine JinBethkeNitrateSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only: GetImmobileSpeciesIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_jinbethke_nitrate_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entry for storing the rate
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 3
  ! Aqueous species
  word = 'H+'
  this%h_ion_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO3-'
  this%no3_id = &
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
  word = 'N2(aq)'
  this%n2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Nim'
  this%nim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

end subroutine JinBethkeNitrateSetup
! ************************************************************************** !
subroutine JinBethkeNitrateAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds denit auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_jinbethke_nitrate_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'JB Nitrate Acetate Sandbox Rate'
  units = 'mol-Ac/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1) 
  word = 'dG-rxn_Nitrate_Acetate Sandbox Generic'
  units = 'kJ/mol-Ac'
  call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+2)

  word = 'Ft_Nitrate_Acetate Sandbox Generic'
  units = 'unitless'
  call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+3)
end subroutine JinBethkeNitrateAuxiliaryPlotVariables
! ************************************************************************** !
subroutine JinBethkeNitrateEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  ! 
  ! Evaluates denitriciation reaction storing residual but no Jacobian
  ! Based on Jin & Bethke (2003, 2005)
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/2/6
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  implicit none
  class(reaction_sandbox_jinbethke_nitrate_type) :: this
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

  PetscReal :: Ac, Proton, no3, Bicarbonate
  PetscReal :: Nim, yield, O2aq, n2aq, Rate_b
  PetscReal :: Rate, Rate_Ac, Rate_Proton, Rate_n2aq
  PetscReal :: Rate_no3, Rate_Bicarbonate
  PetscReal :: stoi_ac, stoi_proton, stoi_n2aq
  PetscReal :: stoi_no3, stoi_bicarbonate
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


  ! Rxn:    1.00 Ac- + 1.60 NO3- + 0.60 H+ = 0.80 N2(aq) + 2.00 HCO3- + 0.80 H2O 
  ! dG0 =   -816.034 kJ per mol Ac-, from Gf0 values, Stumm & Morgan; Kocar & Fendorf 2009 

  Ac = rt_auxvar%pri_molal(this%acetate_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%acetate_id) 
  Proton = rt_auxvar%pri_molal(this%h_ion_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  no3 = rt_auxvar%pri_molal(this%no3_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%no3_id) 
  Bicarbonate = rt_auxvar%pri_molal(this%bicarbonate_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%bicarbonate_id) 
  O2aq = rt_auxvar%pri_molal(this%o2aq_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%o2aq_id) 
  n2aq = rt_auxvar%pri_molal(this%n2aq_id) * molality_to_molarity * &
    rt_auxvar%pri_act_coef(this%n2aq_id) 

  Nim = rt_auxvar%immobile(this%nim_id)

  m = this%m
  chi = this%chi

  stoi_n2aq = 0.8d0
  stoi_bicarbonate = 2.d0

  stoi_no3 = 1.6d0
  stoi_ac = 1.0d0
  stoi_proton = 0.6d0 

  RT = (8.314e-3) * (global_auxvar%temp + 273.15d0)
  dG0 = (-816.034d0) ! kJ / mol acetate; dG0 for NO3- as electron acceptor
  dG_ATP = 50.d0 ! kJ / mol ATP

  reaction_Q = ( (n2aq**stoi_n2aq) * (Bicarbonate**stoi_bicarbonate)) / &
    ((Ac**stoi_ac) * (Proton**stoi_proton) * (no3**stoi_no3))

  dGr = dG0 + RT*log(reaction_Q)
  
  ! Monod expressions for acetate
  Fdonor = Ac / (Ac + this%Kdonor)

  ! Monod expressions for nitrate
  Facceptor = no3 / (no3 + this%Kacceptor)

  ! Thermodynamic factor 
  Ft = 1.d0 - exp((dGr + m*dG_ATP) / (chi * RT))

  if (Ft < 0) then 
    Ftr = 0.d0
  else
    Ftr = Ft
  endif

  ! only denitrification if O2(aq) below threhsold
  calculate_rate = O2aq < (this%o2_threshold)


  Rate = 0.d0
  
  if (calculate_rate) then
    ! base rate, mol/sec/m^3 bulk
    ! units on k: mol/sec/mol-bio

    Rate_b = -k_rmax *  Facceptor * Fdonor * Ftr * Nim * L_water

    Rate = Rate_b  ! mol/sec
      
    ! species-specifc 
    Rate_Ac = Rate * stoi_ac  
    Rate_Proton = Rate * stoi_proton 
    Rate_no3 = Rate * stoi_no3 
    Rate_n2aq = Rate * stoi_n2aq
    Rate_Bicarbonate = Rate * stoi_bicarbonate 
    !Rate_Nim = Rate * yield
       
    rt_auxvar%auxiliary_data(iauxiliary) = Rate_Ac
    rt_auxvar%auxiliary_data(iauxiliary+1) = dGr
    rt_auxvar%auxiliary_data(iauxiliary+2) = Ftr
    
    Residual(this%h_ion_id) = Residual(this%h_ion_id) - Rate_Proton
    Residual(this%acetate_id) = Residual(this%acetate_id) - Rate_Ac
    Residual(this%no3_id) = Residual(this%no3_id) - Rate_no3
    Residual(this%n2aq_id) = Residual(this%n2aq_id) + Rate_n2aq
    Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + Rate_Bicarbonate

  endif
 
end subroutine JinBethkeNitrateEvaluate

end module Reaction_Sandbox_JinBethke_Nitrate_class