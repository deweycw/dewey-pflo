module Reaction_Sandbox_SOM_Acetate_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  use String_module

  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_som_acetate_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: ac_id
    PetscInt :: acetate_id
    PetscInt :: mineral_id
    PetscReal :: rate
    PetscReal :: Km
    PetscReal :: Ct
   

  contains
    procedure, public :: ReadInput => SOMAcetateReadInput
    procedure, public :: Setup => SOMAcetateSetup
    procedure, public :: AuxiliaryPlotVariables => SOMAcetateAuxiliaryPlotVariables
    procedure, public :: Evaluate => SOMAcetateEvaluate
    procedure, public :: UpdateKineticState => SOMAcetateUpdateKineticState
  end type reaction_sandbox_som_acetate_type

  public :: SOMAcetateCreate, &
            SOMAcetateSetup
contains
! ************************************************************************** !
function SOMAcetateCreate()
  !
  ! Allocates SOM_Acetate reaction object.
  !
  implicit none
  class(reaction_sandbox_som_acetate_type), pointer :: SOMAcetateCreate
  allocate(SOMAcetateCreate)
  SOMAcetateCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  SOMAcetateCreate%h_ion_id = UNINITIALIZED_INTEGER
  SOMAcetateCreate%acetate_id = UNINITIALIZED_INTEGER
  SOMAcetateCreate%mineral_id = UNINITIALIZED_INTEGER

  SOMAcetateCreate%rate = UNINITIALIZED_DOUBLE
  SOMAcetateCreate%Km = UNINITIALIZED_DOUBLE
  SOMAcetateCreate%Ct = UNINITIALIZED_DOUBLE

  nullify(SOMAcetateCreate%next)
end function SOMAcetateCreate
! ************************************************************************** !
subroutine SOMAcetateReadInput(this,input,option)
  !
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_som_acetate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,SOM_ACETATE'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(word)
      case('RATE')
        call InputReadDouble(input,option,this%rate)
        call InputErrorMsg(input,option,word,error_string)
      case('INVERSE_KM')
        call InputReadDouble(input,option,this%Km)
        call InputErrorMsg(input,option,word,error_string)
      case('THRESHOLD')
        call InputReadDouble(input,option,this%Ct)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rate) .or. &
      Uninitialized(this%Km) .or. &
      Uninitialized(this%Ct)) then
    option%io_buffer = 'RATE, KM_INVERSE, and THRESHOLD must be set for &
      SOM_ACETATE.'
    call PrintErrMsg(option)
  endif
end subroutine SOMAcetateReadInput
! ************************************************************************** !
subroutine SOMAcetateSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_som_acetate_type) :: this
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
  word = 'Ac-'
  this%acetate_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'SOM'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine SOMAcetateSetup
! ************************************************************************** !
subroutine SOMAcetateAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds ferrihydrite auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_som_acetate_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'SOM Acetate Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
end subroutine SOMAcetateAuxiliaryPlotVariables
! ************************************************************************** !
subroutine SOMAcetateEvaluate(this, Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !Jacobian,compute_derivative,
  ! Evaluates SOM fermentation reaction storing residual but no Jacobian
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
  class(reaction_sandbox_som_acetate_type) :: this
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

  PetscReal :: Ac, Proton
  PetscReal :: Rate, Rate_Ac, Rate_Proton
  PetscReal :: stoi_ac, stoi_proton
  PetscReal :: threshold, rate_from_user, km

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
    option%io_buffer = 'For SOM_ACETATE to function correctly, &
      &the SOM RATE_CONSTANT in the default MINERAL_KINETICS block must be set &
      &to zero.'
    call PrintErrMsg(option)
  endif

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  rate_from_user = this%rate
  threshold = this%Ct
  km = this%Km

  Ac = rt_auxvar%pri_molal(this%acetate_id) * &
    rt_auxvar%pri_act_coef(this%acetate_id) 
  Proton = rt_auxvar%pri_molal(this%h_ion_id) * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  
  stoi_ac = 1.d0
  stoi_proton = 1.d0

  Rate = 0.d0 
  
  ! calculate rate if acetate concentration below threshold
  ! negative for dissolution 
  if (Ac < threshold) then
    Rate = (-1.d0) * rate_from_user * ( km - Ac) / km) 
  endif 

  ! base rate, mol/sec/m^3 bulk
  ! units on k: mol/sec/mol-bio

  rt_auxvar%auxiliary_data(iauxiliary) = Rate 

  Rate = Rate * material_auxvar%volume ! mol/sec
    
  Rate_Ac = Rate * stoi_ac  
  Rate_Proton = Rate * stoi_proton 

  Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate_Proton
  Residual(this%acetate_id) = Residual(this%acetate_id) + Rate_Ac  

end subroutine SOMAcetateEvaluate
! ************************************************************************** !
subroutine SOMAcetateUpdateKineticState(this,rt_auxvar,global_auxvar, &
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
  class(reaction_sandbox_som_acetate_type) :: this
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
end subroutine SOMAcetateUpdateKineticState
end module Reaction_Sandbox_SOM_Acetate_class