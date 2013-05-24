#ifdef GEOMECH

module Geomechanics_Material_module

  implicit none
  
  private
  
#include "definitions.h"  

  type, public :: geomech_material_property_type
    character(len=MAXWORDLENGTH)   :: name
    PetscInt                       :: id
    PetscReal                      :: youngs_modulus
    PetscReal                      :: poissons_ratio
  
    type(geomech_material_property_type), pointer :: next
  end type geomech_material_property_type

  type, public :: geomech_material_property_ptr_type
    type(geomech_material_property_type), pointer :: ptr
  end type geomech_material_property_ptr_type
  
  public :: GeomechanicsMaterialPropertyCreate, &
            GeomechanicsMaterialPropertyDestroy, &
            GeomechanicsMaterialPropertyAddToList, &
            GeomechanicsMaterialPropertyRead, &
            GeomechanicsMaterialPropConvertListToArray, &
            GeomechanicsMaterialPropGetPtrFromArray

contains

! ************************************************************************** !
!
! GeomechanicsMaterialPropertyCreate: Creates a geomechanics material property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
function GeomechanicsMaterialPropertyCreate()

  implicit none
  
  type(geomech_material_property_type), &
    pointer :: GeomechanicsMaterialPropertyCreate
  type(geomech_material_property_type), pointer :: geomech_material_property
  
  allocate(geomech_material_property)
  
  geomech_material_property%name = ''
  geomech_material_property%id = 0
  geomech_material_property%youngs_modulus = 0.d0
  geomech_material_property%poissons_ratio = 0.d0
  
  nullify(geomech_material_property%next)
  
  GeomechanicsMaterialPropertyCreate => geomech_material_property

end function GeomechanicsMaterialPropertyCreate

! ************************************************************************** !
!
! GeomechanicsMaterialPropertyRead: Reads geomechanics material properties
! property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsMaterialPropertyRead(geomech_material_property, &
                                            input,option)

  use Option_module
  use Input_module
  use String_module
  
  implicit none
  
  type(geomech_material_property_type) :: geomech_material_property
  type(input_type)                     :: input
  type(option_type)                    :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string
  
  do
    call InputReadFlotranString(input,option)
    
    if(InputCheckExit(input,option)) exit
  
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS_MATERIAL_PROPERTY')
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('ID')
        call InputReadInt(input,option,geomech_material_property%id)
        call InputErrorMsg(input,option,'id','GEOMECHANICS_MATERIAL_PROPERTY')
      case('YOUNGS_MODULUS')
        call InputReadDouble(input,option,geomech_material_property% &
                             youngs_modulus)
        call InputErrorMsg(input,option,'YOUNGS_MODULUS', &
                           'GEOMECHANICS_MATERIAL_PROPERTY')
      case('POISSONS_RATIO')
        call InputReadDouble(input,option,geomech_material_property% &
                             poissons_ratio)
        call InputErrorMsg(input,option,'POISSONS_RATIO', &
                           'GEOMECHANICS_MATERIAL_PROPERTY')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
          ' not recognized in geomechanics_material_property'
        call printErrMsg(option)
      end select
  enddo
  
end subroutine GeomechanicsMaterialPropertyRead

! ************************************************************************** !
!
! GeomechanicsMaterialPropertyAddToList: Destroys a geomechanics material
! property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsMaterialPropertyAddToList(geomech_material_property,list)

  implicit none
  
  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: list
  type(geomech_material_property_type), pointer :: cur_geomech_material_property
  
  if (associated(list)) then
    cur_geomech_material_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_geomech_material_property%next)) exit
      cur_geomech_material_property => cur_geomech_material_property%next
    enddo
    cur_geomech_material_property%next => geomech_material_property
  else
    list => geomech_material_property
  endif
  
end subroutine GeomechanicsMaterialPropertyAddToList

! ************************************************************************** !
!
! GeomechanicsMaterialPropConvertListToArray: Destroys a geomechanics material
! property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsMaterialPropConvertListToArray(list,array,option)

  use Option_module
  use String_module

  implicit none

  type(geomech_material_property_type), pointer :: list
  type(geomech_material_property_ptr_type), pointer :: array(:)
  type(option_type) :: option

  type(geomech_material_property_type), pointer :: cur_material_property
  type(geomech_material_property_type), pointer :: prev_material_property
  type(geomech_material_property_type), pointer :: next_material_property
  PetscInt :: i, j, length1,length2, max_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

  max_id = 0
  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    max_id = max(max_id,cur_material_property%id)
    cur_material_property => cur_material_property%next
  enddo

  allocate(array(max_id))
  do i = 1, max_id
    nullify(array(i)%ptr)
  enddo

  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_id))
  id_count = 0

  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    id_count(cur_material_property%id) = &
      id_count(cur_material_property%id) + 1
    array(cur_material_property%id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo

  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Material ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call printMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Material IDs.'
    call printErrMsg(option)
  endif

  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Material name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call printMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Material names.'
    call printErrMsg(option)
  endif

end subroutine GeomechanicsMaterialPropConvertListToArray


! ************************************************************************** !
!
! GeomechanicsMaterialPropGetPtrFromArray: Destroys a geomechanics material
! property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
function GeomechanicsMaterialPropGetPtrFromArray( &
                                            geomech_material_property_name, &
                                            geomech_material_property_array)

  use String_module

  implicit none

  type(geomech_material_property_type), &
    pointer :: GeomechanicsMaterialPropGetPtrFromArray
  type(geomech_material_property_ptr_type), &
    pointer :: geomech_material_property_array(:)
  character(len=MAXWORDLENGTH) :: geomech_material_property_name
  PetscInt :: length
  PetscInt :: igeomech_material_property

  nullify(GeomechanicsMaterialPropGetPtrFromArray)

  do igeomech_material_property = 1, size(geomech_material_property_array)
    length = len_trim(geomech_material_property_name)
    if (.not.associated(geomech_material_property_array &
      (igeomech_material_property)%ptr)) cycle
    if (length == &
        len_trim(geomech_material_property_array &
          (igeomech_material_property)%ptr%name) .and. &
        StringCompare(geomech_material_property_array &
          (igeomech_material_property)%ptr%name, &
                        geomech_material_property_name,length)) then
      GeomechanicsMaterialPropGetPtrFromArray => &
        geomech_material_property_array(igeomech_material_property)%ptr
      return
    endif
  enddo

end function GeomechanicsMaterialPropGetPtrFromArray

! ************************************************************************** !
!
! GeomechanicsMaterialPropertyDestroy: Destroys a geomechanics material
! property
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
recursive subroutine GeomechanicsMaterialPropertyDestroy(geomech_material_property)

  implicit none
  
  type(geomech_material_property_type), pointer :: geomech_material_property
  
  if(.not.associated(geomech_material_property)) return
  
  call GeomechanicsMaterialPropertyDestroy(geomech_material_property%next)
  
  deallocate(geomech_material_property)
  nullify(geomech_material_property)
  
end subroutine GeomechanicsMaterialPropertyDestroy

end module Geomechanics_Material_module
#endif 
!GEOMECH