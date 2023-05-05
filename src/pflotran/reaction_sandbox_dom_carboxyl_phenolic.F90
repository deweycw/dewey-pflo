module Reaction_Sandbox_Functionalized_DOM_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_functionalized_dom_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: dom_id
    PetscInt :: hdom_id
    PetscInt :: h2dom_id
    PetscReal :: n_site1
    PetscReal :: n_site2
    PetscReal :: ka1
    PetscReal :: ka2
    PetscReal :: kf
    PetscReal :: kr


  contains
    procedure, public :: ReadInput => FcnDOMReadInput
    procedure, public :: Setup => FcnDOMSetup
    procedure, public :: AuxiliaryPlotVariables => FcnDOMAuxiliaryPlotVariables
    procedure, public :: Evaluate => FcnDOMEvaluate
  end type reaction_sandbox_functionalized_dom_type

  public :: FcnDOMCreate, &
            FcnDOMSetup
contains
! ************************************************************************** !
function FcnDOMCreate()
  !
  ! Allocates Functionalized DOM reaction object.
  !
  implicit none
  class(reaction_sandbox_functionalized_dom_type), pointer :: FcnDOMCreate
  allocate(FcnDOMCreate)
  FcnDOMCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  FcnDOMCreate%h_ion_id = UNINITIALIZED_INTEGER
  FcnDOMCreate%dom_id = UNINITIALIZED_INTEGER
  FcnDOMCreate%hdom_id = UNINITIALIZED_INTEGER
  FcnDOMCreate%h2dom_id = UNINITIALIZED_INTEGER

  FcnDOMCreate%ka1 = UNINITIALIZED_DOUBLE
  FcnDOMCreate%ka2 = UNINITIALIZED_DOUBLE
  FcnDOMCreate%kf = UNINITIALIZED_DOUBLE
  FcnDOMCreate%kr = UNINITIALIZED_DOUBLE
  FcnDOMCreate%n_site1 = UNINITIALIZED_DOUBLE
  FcnDOMCreate%n_site2 = UNINITIALIZED_DOUBLE

  nullify(FcnDOMCreate%next)
end function FcnDOMCreate
! ************************************************************************** !
subroutine FcnDOMReadInput(this,input,option)
  !
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_functionalized_dom_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FUNCTIONALIZED_DOM'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('KA_SITE_1')
        call InputReadDouble(input,option,this%ka1)
        call InputErrorMsg(input,option,word,error_string)
      case('KA_SITE_2')
        call InputReadDouble(input,option,this%ka2)
        call InputErrorMsg(input,option,word,error_string)
      case('FRACTION_SITE_1')
        call InputReadDouble(input,option,this%n_site1)
        call InputErrorMsg(input,option,word,error_string)
      case('FRACTION_SITE_2')
        call InputReadDouble(input,option,this%n_site2)
        call InputErrorMsg(input,option,word,error_string)
    case('K_FORWARD')
        call InputReadDouble(input,option,this%kf)
        call InputErrorMsg(input,option,word,error_string)
      case('K_REVERSE')
        call InputReadDouble(input,option,this%kr)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%ka1) .or. &
      Uninitialized(this%ka2) .or. &
      Uninitialized(this%n_site1) .or. &
      Uninitialized(this%n_site2) .or. &
      Uninitialized(this%kf) .or. &
      Uninitialized(this%kr)) then
    option%io_buffer = 'KA1&2, Kf&r, site fractions all must be set for &
      FUNCTIONALIZED_DOM'
    call PrintErrMsg(option)
  endif
end subroutine FcnDOMReadInput
! ************************************************************************** !
subroutine FcnDOMSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_functionalized_dom_type) :: this
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
  word = 'HDOM-'
  this%hdom_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'DOM--'
  this%dom_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H2DOM(aq)'
  this%h2dom_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

end subroutine FcnDOMSetup
! ************************************************************************** !
subroutine FcnDOMAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds dissociation auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_functionalized_dom_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'Ka1 DOM Dissociation Sandbox Rate'
  units = 'mol/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  word = 'Ka2 DOM Dissociation Sandbox Rate'
  units = 'mol/sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+2) 
end subroutine FcnDOMAuxiliaryPlotVariables
! ************************************************************************** !
subroutine FcnDOMEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !Jacobian,compute_derivative,
  ! Evaluates dissociation reaction storing residual but no Jacobian
  !
  ! 
  !
  !
  ! Author: Christian Dewey
  ! Date: 2023/3/18

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  implicit none
  class(reaction_sandbox_functionalized_dom_type) :: this
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
  PetscReal :: L_water              ! L water

  PetscReal :: DOM, HDOM, H2DOM, Proton
  PetscReal :: kf, kr
  PetscReal :: ka_1, ka_2
  PetscReal :: reaction1_Q, reaction2_Q

  PetscReal :: Rate1, Rate2
  PetscReal :: lnQK_1, sign1_, affinity_factor_1, QK_1
  PetscReal :: lnQK_2, sign2_, affinity_factor_2, QK_2
  PetscReal :: fraction1, fraction2
  PetscBool :: calculate_forward_1, calculate_reverse_1
  PetscBool :: calculate_forward_2, calculate_reverse_2

  PetscInt :: jcomp, icomp
  PetscInt :: ncomp, i 
  PetscInt :: iauxiliary, jauxiliary

  iauxiliary = this%auxiliary_offset + 1
  iauxiliary = this%auxiliary_offset + 2

  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  Proton = rt_auxvar%pri_molal(this%h_ion_id) * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  DOM = rt_auxvar%pri_molal(this%dom_id) * &
    rt_auxvar%pri_act_coef(this%dom_id) 
  HDOM = rt_auxvar%pri_molal(this%hdom_id) * &
    rt_auxvar%pri_act_coef(this%hdom_id) 
  H2DOM = rt_auxvar%pri_molal(this%h2dom_id) * &
    rt_auxvar%pri_act_coef(this%h2dom_id) 

  ! Rxn1:    1.00 H2DOM(aq) = n_1 HDOM- + n_1 H+ 
  ! ka_1 ~= 4

  fraction1 = this%n_site1 / (this%n_site1 + this%n_site2)

  reaction1_Q = ( (Proton**fraction1) * (HDOM**fraction1)) / &
    (H2DOM) 
  
  lnQK_1 = this%ka1*LOG_TO_LN  ! uses ka_1 defined by user

  lnQK_1 = lnQK_1 + log(reaction1_Q)

  QK_1 = exp(lnQK_1)
  affinity_factor_1 = 1.d0-QK_1
  sign1_ = sign(1.d0,affinity_factor_1)

  calculate_forward_1 = (sign1_>0)
  calculate_reverse_1 = (sign1_<0)

  ! Rxn2:    1.00 HDOM-(aq) = n_2 HDOM- + n_2 H+ 
  ! ka_2 ~= 8
  
  fraction2 = this%n_site2 / (this%n_site1 + this%n_site2)

  reaction2_Q = ( (Proton**fraction2) * (DOM**fraction2)) / &
    (HDOM) 
  
  lnQK_2 = this%ka2*LOG_TO_LN  ! uses ka_2 defined by user

  lnQK_2 = lnQK_2 + log(reaction2_Q)

  QK_2 = exp(lnQK_2)
  affinity_factor_2 = 1.d0-QK_2
  sign2_ = sign(1.d0,affinity_factor_2)

  calculate_forward_2 = (sign2_>0)
  calculate_reverse_2 = (sign2_<0)


  Rate1 = 0.d0
  
  if (calculate_forward_1) then

    Rate1 = kf * H2DOM * sign1_ * abs(affinity_factor_1) 
   
    rt_auxvar%auxiliary_data(iauxiliary) = Rate1

    Rate1 = Rate1 * material_auxvar%volume ! mol/sec
         
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate1
    Residual(this%hdom_id) = Residual(this%hdom_id) + Rate1
    Residual(this%h2dom_id) = Residual(this%h2dom_id) - Rate1

  end if 

  if (calculate_reverse_1) then

    Rate1 = kr * HDOM * Proton * sign1_ * abs(affinity_factor_1) 
   
    rt_auxvar%auxiliary_data(iauxiliary) = Rate1

    Rate1 = Rate1 * material_auxvar%volume ! mol/sec
         
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate1
    Residual(this%hdom_id) = Residual(this%hdom_id) + Rate1
    Residual(this%h2dom_id) = Residual(this%h2dom_id) - Rate1

  end if 


  Rate2 = 0.d0

  if (calculate_forward_2) then

    Rate2 = kf * HDOM * sign2_ * abs(affinity_factor_2) 
   
    rt_auxvar%auxiliary_data(jauxiliary) = Rate2

    Rate2 = Rate2 * material_auxvar%volume ! mol/sec
         
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate2
    Residual(this%dom_id) = Residual(this%dom_id) + Rate2
    Residual(this%hdom_id) = Residual(this%hdom_id) - Rate2

  end if 

  if (calculate_reverse_2) then

    Rate2 = kr * DOM * Proton * sign2_ * abs(affinity_factor_2) 
   
    rt_auxvar%auxiliary_data(jauxiliary) = Rate2

    Rate2 = Rate2 * material_auxvar%volume ! mol/sec
         
    Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate2
    Residual(this%dom_id) = Residual(this%dom_id) + Rate2
    Residual(this%hdom_id) = Residual(this%hdom_id) - Rate2

  end if 

end subroutine FcnDOMEvaluate
end module Reaction_Sandbox_Functionalized_DOM_class