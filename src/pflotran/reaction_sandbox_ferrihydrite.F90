module Reaction_Sandbox_Ferrihydrite_class
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  implicit none
  private
  type, public, &
    extends(reaction_sandbox_base_type) :: &
      reaction_sandbox_ferrihydrite_type
    PetscInt :: auxiliary_offset
    PetscInt :: h_ion_id
    PetscInt :: fe2_id
    PetscInt :: acetate_id
    PetscInt :: bicarbonate_id
    PetscInt :: mineral_id
    PetscInt :: xim_id
    PetscReal :: rate_constant1
    PetscReal :: Kd
    PetscReal :: Y
    PetscReal :: m
    PetscReal :: chi
  contains
    procedure, public :: ReadInput => FerrihydriteReadInput
    procedure, public :: Setup => FerrihydriteSetup
    procedure, public :: AuxiliaryPlotVariables => FerrihydriteAuxiliaryPlotVariables
    procedure, public :: Evaluate => FerrihydriteEvaluate
    procedure, public :: UpdateKineticState => FerrihydriteUpdateKineticState
  end type reaction_sandbox_ferrihydrite_type

  public :: FerrihydriteCreate, &
            FerrihydriteSetup
contains
! ************************************************************************** !
function FerrihydriteCreate()
  !
  ! Allocates Ferrihydrite reaction object.
  !
  implicit none
  class(reaction_sandbox_ferrihydrite_type), pointer :: FerrihydriteCreate
  allocate(FerrihydriteCreate)
  FerrihydriteCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  FerrihydriteCreate%h_ion_id = UNINITIALIZED_INTEGER
  FerrihydriteCreate%fe2_id = UNINITIALIZED_INTEGER
  FerrihydriteCreate%acetate_id = UNINITIALIZED_INTEGER
  FerrihydriteCreate%bicarbonate_id = UNINITIALIZED_INTEGER
  FerrihydriteCreate%mineral_id = UNINITIALIZED_INTEGER
  FerrihydriteCreate%xim_id = UNINITIALIZED_INTEGER

  FerrihydriteCreate%rate_constant1 = UNINITIALIZED_DOUBLE
  FerrihydriteCreate%Kd = UNINITIALIZED_DOUBLE
  FerrihydriteCreate%Y = UNINITIALIZED_DOUBLE
  FerrihydriteCreate%m = UNINITIALIZED_DOUBLE
  FerrihydriteCreate%chi = UNINITIALIZED_DOUBLE

  nullify(FerrihydriteCreate%next)
end function FerrihydriteCreate
! ************************************************************************** !
subroutine FerrihydriteReadInput(this,input,option)
  !
  ! Reads calcite reaction parameters
  !
  use Option_module
  use Input_Aux_module
  use String_module
  implicit none
  class(reaction_sandbox_ferrihydrite_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FERRIHYDRITE'
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
        call InputReadDouble(input,option,this%rate_constant1)
        call InputErrorMsg(input,option,word,error_string)
      case('KD')
        call InputReadDouble(input,option,this%Kd)
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
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (Uninitialized(this%rate_constant1) .or. &
      Uninitialized(this%Kd) .or. &
      Uninitialized(this%Y) .or. &
      Uninitialized(this%m) .or. &
      Uninitialized(this%chi)) then
    option%io_buffer = 'RATE_CONSTANT, KD, Y, M, and CHI must be set for &
      REACTION_SANDBOX_FERRIHYDRITE.'
    call PrintErrMsg(option)
  endif
end subroutine FerrihydriteReadInput
! ************************************************************************** !
subroutine FerrihydriteSetup(this,reaction,option)
  !
  ! Sets up the calcite reaction with hardwired parameters
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module
  implicit none
  class(reaction_sandbox_ferrihydrite_type) :: this
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
  word = 'Ferrihydrite'
  this%mineral_id = &
    GetMineralIDFromName(word,reaction%mineral,option)

end subroutine FerrihydriteSetup
! ************************************************************************** !
subroutine FerrihydriteAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds ferrihydrite auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY
  class(reaction_sandbox_ferrihydrite_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  word = 'Ferrihydrite Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)
  
  !word = 'Thermodynamic Factor (FT)'
  !units = 'mol/m^3-sec'
  !call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
   !                             REACTION_AUXILIARY, &
    !                            this%auxiliary_offset+1)
  
end subroutine FerrihydriteAuxiliaryPlotVariables
! ************************************************************************** !
subroutine FerrihydriteEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates ferrihydrite reaction storing residual but no Jacobian
  !
  ! 
  !
  ! Rxn:    1.00 Ac- + 8.00 FHY + 15.00 H+ = 8.00 Fe++ + 2.00 HCO3- + 20.00 H2O 
  ! dG0 =   -586.85 kJ per mol Ac- 
  !
  ! Author: Christian Dewey
  ! Date: 2022/10/4
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Mineral_Aux_module
  implicit none
  class(reaction_sandbox_ferrihydrite_type) :: this
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
  !PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: L_water              ! L water
 !PetscReal :: ln_act(reaction%ncomp)
 !PetscReal :: drate_xim
  PetscReal :: drate, drate_ac
  PetscReal :: drate_bicarbonate
  PetscReal :: drate_fe2, drate_h

  PetscReal :: Ac, Proton, Fe2, Bicarbonate
  PetscReal :: Xim, yield
  PetscReal :: Rate, Rate_Ac, Rate_Proton !Rate_Xim,
  PetscReal :: Rate_Fe2, Rate_Bicarbonate
  PetscReal :: Rate_Fh_dissoluton
  PetscReal :: stoi_ac, stoi_proton
  PetscReal :: stoi_fe2, stoi_bicarbonate
  PetscReal :: k, m, chi
  PetscReal :: temp_K, RT
  PetscReal :: Ft, Ftr, Fd
  PetscReal :: reaction_Q
  PetscReal :: dG0, dGr, dG_ATP

  PetscInt :: jcomp
  PetscInt :: imnrl
  PetscInt :: iauxiliary
  PetscBool :: calculate_rate

  Ac = 0.d0
  Proton = 0.d0
  Fe2 = 0.d0
  Bicarbonate = 0.d0
  Xim = 0.d0
  yield = 0.d0
  Rate = 0.d0
  !Rate_Xim = 0.d0
  Rate_Ac = 0.d0
  Rate_Proton = 0.d0
  Rate_Fe2 = 0.d0
  Rate_Bicarbonate = 0.d0
  Rate_Fh_dissoluton = 0.d0
  drate = 0.d0
  drate_ac = 0.d0
  drate_bicarbonate = 0.d0
  drate_h = 0.d0
  drate_fe2 = 0.d0
  !drate_xim = 0.d0

  k = 0.d0
  Ft = 0.d0
  Ftr = 0.d0
  Fd = 0.d0
  reaction_Q = 0.d0
  dGr = 0.d0

  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  volume = material_auxvar%volume
  L_water = porosity*liquid_saturation*volume*1.d3

  Ac = rt_auxvar%pri_molal(this%acetate_id) * &
    rt_auxvar%pri_act_coef(this%acetate_id) 
  Proton = rt_auxvar%pri_molal(this%h_ion_id) * &
    rt_auxvar%pri_act_coef(this%h_ion_id) 
  Fe2 = rt_auxvar%pri_molal(this%fe2_id) * &
    rt_auxvar%pri_act_coef(this%fe2_id) 
  Bicarbonate = rt_auxvar%pri_molal(this%bicarbonate_id) * &
    rt_auxvar%pri_act_coef(this%bicarbonate_id) 

  Xim = rt_auxvar%immobile(this%xim_id)

  m = this%m
  chi = this%chi

  stoi_fe2 = 8.d0
  stoi_bicarbonate = 2.d0
  
  stoi_ac = 1.d0
  stoi_proton = 15.d0

  RT = (8.314e-3) * (global_auxvar%temp + 273.15d0)
  dG0 = -586.85d0 ! kJ / mol acetate; dG0 for FeIII in ferrihydrite as electron acceptor
  dG_ATP = -50.d0 ! kJ / mol ATP

  reaction_Q = ( Fe2**stoi_fe2 * Bicarbonate**stoi_bicarbonate) / &
    (Ac**stoi_ac * Proton**stoi_proton)

  dGr = dG0 + RT*log(reaction_Q)
  
  ! Monod expressions
  Fd = Ac / (Ac + this%Kd)

  ! Thermodynamic factor 
  Ft = 1.d0 - exp((dGr + m*dG_ATP) / &
    (chi * RT))

  if (Ft < 0) then 
    Ftr = 0.d0
  else
    Ftr = Ft
  endif

  k = this%rate_constant1

  mineral => reaction%mineral
  iauxiliary = this%auxiliary_offset + 1
  volume = material_auxvar%volume        ! den_kg [kg fluid / m^3 fluid]
  !molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3  ! kg water/L water
 
  ! only calculate rate if mineral is present 
  calculate_rate = rt_auxvar%mnrl_volfrac(imnrl) > 0 
  if (calculate_rate) then
    ! mol/sec/m^3 bulk
    !rate = -rt_auxvar%mnrl_area(imnrl) * &
          ! sign_ * abs(affinity_factor) * this%rate_constant2

    ! overall rate 
    Rate = k * Xim * Fd * Ftr * L_water

    ! species-specifc 
    !Rate_Xim = Rate * yield
    Rate_Ac = Rate * (-1.d0)*stoi_ac
    Rate_Proton = Rate * (-1.d0)*stoi_proton
    Rate_Fe2 = Rate * stoi_fe2
    Rate_Bicarbonate = Rate * stoi_bicarbonate

    ! ferrihydrite dissolution rate 
    Rate_Fh_dissoluton = Rate_Fe2 * (-1.d0)

  endif
  rt_auxvar%auxiliary_data(iauxiliary) = &
    rt_auxvar%auxiliary_data(iauxiliary) + Rate_Fh_dissoluton
  ! mol/sec
  Rate_Fh_dissoluton = Rate_Fh_dissoluton * material_auxvar%volume
  
  Residual(this%h_ion_id) = Residual(this%h_ion_id) + Rate_Proton
  Residual(this%acetate_id) = Residual(this%acetate_id) + Rate_Ac
  Residual(this%fe2_id) = Residual(this%fe2_id) + Rate_Fe2
  Residual(this%bicarbonate_id) = Residual(this%bicarbonate_id) + Rate_Bicarbonate
  if (compute_derivative .and. calculate_rate) then
    ! derivative of rate wrt affinity factor (1-QK)
    ! mol/sec   ! m^2 mnrl/m^3 bulk          ! mol/m^2 mnrl/sec
    drate = rt_auxvar%mnrl_area(imnrl) * this%rate_constant1 * &
                ! m^3 bulk
                material_auxvar%volume

    drate_h = drate * stoi_proton
    drate_ac = drate * stoi_ac
    drate_fe2 = drate * stoi_fe2
    drate_bicarbonate = drate * stoi_bicarbonate
    !drate_xim = drate * yield

    ! derivative wrt H+
    jcomp = this%h_ion_id
                                    ! subtract due to H+ stoichiometry
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_h
                                    ! subtract due to Ac- stoichiometry
    Jacobian(this%acetate_id,jcomp) = &
      Jacobian(this%acetate_id,jcomp) - drate_ac
                                    ! add due to HCO3- stoichiometry
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_bicarbonate
                                    ! add due to Fe2+ stoichiometry
    Jacobian(this%fe2_id,jcomp) = &
      Jacobian(this%fe2_id,jcomp) + drate_fe2
                                    ! add due to Xim stoichiometry
    !Jacobian(this%xim_id,jcomp) = &
    !  Jacobian(this%xim_id,jcomp) + drate_xim

    ! derivative wrt Ac-
    jcomp = this%acetate_id
                                    ! subtract due to H+ stoichiometry
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_h
                                    ! subtract due to Ac- stoichiometry
    Jacobian(this%acetate_id,jcomp) = &
      Jacobian(this%acetate_id,jcomp) - drate_ac
                                    ! add due to HCO3- stoichiometry
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_bicarbonate
                                    ! add due to Fe2+ stoichiometry
    Jacobian(this%fe2_id,jcomp) = &
      Jacobian(this%fe2_id,jcomp) + drate_fe2
                                    ! add due to Xim stoichiometry
    !Jacobian(this%xim_id,jcomp) = &
    !  Jacobian(this%xim_id,jcomp) + drate_xim

    ! derivative wrt HCO3-
    jcomp = this%bicarbonate_id
                                    ! subtract due to H+ stoichiometry
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_h
                                    ! subtract due to Ac- stoichiometry
    Jacobian(this%acetate_id,jcomp) = &
      Jacobian(this%acetate_id,jcomp) - drate_ac
                                    ! add due to HCO3- stoichiometry
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_bicarbonate
                                    ! add due to Fe2+ stoichiometry
    Jacobian(this%fe2_id,jcomp) = &
      Jacobian(this%fe2_id,jcomp) + drate_fe2
                                    ! add due to Xim stoichiometry
    !Jacobian(this%xim_id,jcomp) = &
    !  Jacobian(this%xim_id,jcomp) + drate_xim

    ! derivative wrt Fe2+
    jcomp = this%fe2_id
                                    ! subtract due to H+ stoichiometry
    Jacobian(this%h_ion_id,jcomp) = &
      Jacobian(this%h_ion_id,jcomp) - drate_h
                                    ! subtract due to Ac- stoichiometry
    Jacobian(this%acetate_id,jcomp) = &
      Jacobian(this%acetate_id,jcomp) - drate_ac
                                    ! add due to HCO3- stoichiometry
    Jacobian(this%bicarbonate_id,jcomp) = &
      Jacobian(this%bicarbonate_id,jcomp) + drate_bicarbonate
                                    ! add due to Fe2+ stoichiometry
    Jacobian(this%fe2_id,jcomp) = &
      Jacobian(this%fe2_id,jcomp) + drate_fe2
                                    ! add due to Xim stoichiometry
    !Jacobian(this%xim_id,jcomp) = &
    !  Jacobian(this%xim_id,jcomp) + drate_xim

    ! derivative wrt Xim
    !jcomp = this%xim_id
                                    ! subtract due to H+ stoichiometry
    !Jacobian(this%h_ion_id,jcomp) = &
    !  Jacobian(this%h_ion_id,jcomp) - drate_h
                                    ! subtract due to Ac- stoichiometry
    !Jacobian(this%acetate_id,jcomp) = &
    !  Jacobian(this%acetate_id,jcomp) - drate_ac
                                    ! add due to HCO3- stoichiometry
    !Jacobian(this%bicarbonate_id,jcomp) = &
    !  Jacobian(this%bicarbonate_id,jcomp) + drate_bicarbonate
                                    ! add due to Fe2+ stoichiometry
    !Jacobian(this%fe2_id,jcomp) = &
    !  Jacobian(this%fe2_id,jcomp) + drate_fe2
                                    ! add due to Xim stoichiometry
    !Jacobian(this%xim_id,jcomp) = &
    !  Jacobian(this%xim_id,jcomp) + drate_xim

  endif
end subroutine FerrihydriteEvaluate
! ************************************************************************** !
subroutine FerrihydriteUpdateKineticState(this,rt_auxvar,global_auxvar, &
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
  class(reaction_sandbox_ferrihydrite_type) :: this
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
end subroutine FerrihydriteUpdateKineticState
end module Reaction_Sandbox_Ferrihydrite_class