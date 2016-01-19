module array_swap_mod

implicit none

interface swap
! This interface allows us to call swap(array, index1, index2) with an array
! of default kind integers, reals or doubles.
   module procedure swap_int, swap_real, swap_double, &
                    swap_list_int, swap_list_real, swap_list_double
end interface swap

contains

subroutine swap_int(array, index1, index2)
! Rearranges two array elements of an array.
!
! Example call:
!
! array = (/ 1, 2, 3, 5 /)
! call swap_int(array, 2, 4)
! array: (/ 1, 5, 3, 2 /) 
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element
! index2(in): index of the second element
!
   implicit none
! Arguments
   integer, dimension(:), intent(in out)  :: array
   integer, intent(in)                    :: index1 
   integer, intent(in)                    :: index2 
! Variables
   integer                                :: temp
! Do the swap
   temp          = array(index1)
   array(index1) = array(index2)
   array(index2) = temp
end subroutine swap_int

subroutine swap_real(array, index1, index2)
! Rearranges two array elements of an array.
!
! Example call:
!
! array = (/ 1.0, 2.0, 3.0, 5.0 /)
! call swap_real(array, 2, 4)
! array: (/ 1.0, 5.0, 3.0, 2.0 /) 
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element
! index2(in): index of the second element
!
   implicit none
! Arguments
   real, dimension(:), intent(in out)  :: array
   integer, intent(in)                 :: index1 
   integer, intent(in)                 :: index2 
! Variables
   real                                :: temp
! Do the swap
   temp          = array(index1)
   array(index1) = array(index2)
   array(index2) = temp
end subroutine swap_real

subroutine swap_double(array, index1, index2)
! Rearranges two array elements of an array.
!
! Example call:
!
! array = (/ 1.0d0, 2.0d0, 3.0d0, 5.0d0 /)
! call swap_double(array, 2, 4)
! array: (/ 1.0d0, 5.0d0, 3.0d0, 2.0d0 /) 
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element
! index2(in): index of the second element
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in out)  :: array
   integer, intent(in)                             :: index1 
   integer, intent(in)                             :: index2 
! Variables
   double precision                                :: temp
! Do the swap
   temp          = array(index1)
   array(index1) = array(index2)
   array(index2) = temp
end subroutine swap_double

subroutine swap_list_int(array, index1, index2)
! Rearranges two array elements of an array of arrays.
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element, into the first dimension of the array
! index2(in): index of the second element, into the second dimension of the
!             array
!
   implicit none
! Arguments
   integer, dimension(:,:), intent(in out)   :: array
   integer, intent(in)                       :: index1 
   integer, intent(in)                       :: index2 
! Variables
   integer, dimension(:), allocatable        :: temp

! Allocate the temporary array.
   allocate(temp(size(array(index1,:))))

! Now swap the arrays in the array of arrays.
   temp(:)         = array(index1,:)
   array(index1,:) = array(index2,:)
   array(index2,:) = temp(:)
end subroutine swap_list_int

subroutine swap_list_real(array, index1, index2)
! Rearranges two array elements of an array of arrays.
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element, into the first dimension of the array
! index2(in): index of the second element, into the second dimension of the
!             array
!
   implicit none
! Arguments
   real, dimension(:,:), intent(in out)   :: array
   integer, intent(in)                    :: index1 
   integer, intent(in)                    :: index2 
! Variables
   real, dimension(:), allocatable        :: temp

! Allocate the temporary array.
   allocate(temp(size(array(index1,:))))

! Now swap the arrays in the array of arrays.
   temp(:)         = array(index1,:)
   array(index1,:) = array(index2,:)
   array(index2,:) = temp(:)
end subroutine swap_list_real

subroutine swap_list_double(array, index1, index2)
! Rearranges two array elements of an array of arrays.
!
! Arguments
! ---------
! array(in/out): array whose elements are to be swapped
! index1(in): index of the first element, into the first dimension of the array
! index2(in): index of the second element, into the second dimension of the
!             array
!
   implicit none
! Arguments
   double precision, dimension(:,:), intent(in out)   :: array
   integer, intent(in)                                :: index1 
   integer, intent(in)                                :: index2 
! Variables
   double precision, dimension(:), allocatable        :: temp

! Allocate the temporary array.
   allocate(temp(size(array(index1,:))))

! Now swap the arrays in the array of arrays.
   temp(:)         = array(index1,:)
   array(index1,:) = array(index2,:)
   array(index2,:) = temp(:)
end subroutine swap_list_double

end module array_swap_mod
