module sorts_mod

use array_swap_mod

implicit none

interface quicksort
! This interface allows us to call quicksort(data_array, order) with an array
! of default kind integers, reals or doubles.
   module procedure quicksort_int, quicksort_real, quicksort_double, &
                    quicksort_int_list, quicksort_real_list, &
                    quicksort_double_list
end interface quicksort

contains

recursive subroutine quicksort_int(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   integer, dimension(:), intent(in out)           :: data_array
   integer, dimension(:), intent(in out)           :: order
! Variables
   real                                            :: random
   integer                                         :: pivot
   integer                                         :: i
   integer                                         :: store_index
   integer                                         :: right

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array) /= size(order)) then
      stop 'Error in quicksort_int, the data_array and order arrays&are not &
           &the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         if (data_array(i) <= data_array(right)) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1), order(:store_index-1))
      call quicksort(data_array(store_index+1:), order(store_index+1:))
   end if
end subroutine quicksort_int

recursive subroutine quicksort_real(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   real, dimension(:), intent(in out)              :: data_array
   integer, dimension(:), intent(in out)           :: order
! Variables
   real                                            :: random
   integer                                         :: pivot
   integer                                         :: i
   integer                                         :: store_index
   integer                                         :: right

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array) /= size(order)) then
      stop 'Error in quicksort_real, the data_array and order arrays are not &
           &the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         if (data_array(i) <= data_array(right)) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1), order(:store_index-1))
      call quicksort(data_array(store_index+1:), order(store_index+1:))
   end if
end subroutine quicksort_real

recursive subroutine quicksort_double(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in out)  :: data_array
   integer, dimension(:), intent(in out)           :: order
! Variables
   real                                            :: random
   integer                                         :: pivot
   integer                                         :: i
   integer                                         :: store_index
   integer                                         :: right

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array) /= size(order)) then
      stop 'Error in quicksort_double, the data_array and order arrays are not &
           &the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         if (data_array(i) <= data_array(right)) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1), order(:store_index-1))
      call quicksort(data_array(store_index+1:), order(store_index+1:))
   end if
end subroutine quicksort_double

recursive subroutine quicksort_int_list(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   integer, dimension(:,:), intent(in out)            :: data_array
   integer, dimension(:), intent(in out)              :: order
! Variables
   real                                               :: random
   integer                                            :: pivot
   integer                                            :: i
   integer                                            :: j
   integer                                            :: store_index
   integer                                            :: right
   logical, dimension(:), allocatable                 :: array_lt
   logical, dimension(:), allocatable                 :: array_gt

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array, 1) /= size(order)) then
      stop 'Error in quicksort_int_list, the data_array and order arrays &
           &are not the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array, 1) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array, 1) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array, 1)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Allocate two arrays to store the list of less than and greater than
   ! results for comparing each subarray.
      allocate(array_lt(size(data_array, 2)))
      allocate(array_gt(size(data_array, 2)))
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         array_lt = data_array(i,:) < data_array(right,:) 
         array_gt = data_array(i,:) > data_array(right,:) 
      ! First do the simple check whether or not all the values of array_lt and
      ! array_gt are false (this means the arrays are equal).
         if (all(.not.(array_lt .or. array_gt))) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         else
         ! We order lists with the most significant value in the smallest index
            do j = 1, size(array_lt)
               if (array_lt(j) .and. .not.(array_gt(j))) then
                  call swap(data_array, i, store_index)
                  call swap(order, i, store_index)
                  store_index = store_index + 1
                  exit
               else if (array_gt(j) .and. .not.(array_lt(j))) then
                  exit
               end if
            end do
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1,:), order(:store_index-1))
      call quicksort(data_array(store_index+1:,:), order(store_index+1:))
   end if
end subroutine quicksort_int_list

recursive subroutine quicksort_real_list(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   real, dimension(:,:), intent(in out)               :: data_array
   integer, dimension(:), intent(in out)              :: order
! Variables
   real                                               :: random
   integer                                            :: pivot
   integer                                            :: i
   integer                                            :: j
   integer                                            :: store_index
   integer                                            :: right
   logical, dimension(:), allocatable                 :: array_lt
   logical, dimension(:), allocatable                 :: array_gt

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array, 1) /= size(order)) then
      stop 'Error in quicksort_real_list, the data_array and order arrays &
           &are not the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array, 1) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array, 1) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array, 1)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Allocate two arrays to store the list of less than and greater than
   ! results for comparing each subarray.
      allocate(array_lt(size(data_array, 2)))
      allocate(array_gt(size(data_array, 2)))
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         array_lt = data_array(i,:) < data_array(right,:) 
         array_gt = data_array(i,:) > data_array(right,:) 
      ! First do the simple check whether or not all the values of array_lt and
      ! array_gt are false (this means the arrays are equal).
         if (all(.not.(array_lt .or. array_gt))) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         else
         ! We order lists with the most significant value in the smallest index
            do j = 1, size(array_lt)
               if (array_lt(j) .and. .not.(array_gt(j))) then
                  call swap(data_array, i, store_index)
                  call swap(order, i, store_index)
                  store_index = store_index + 1
                  exit
               else if (array_gt(j) .and. .not.(array_lt(j))) then
                  exit
               end if
            end do
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1,:), order(:store_index-1))
      call quicksort(data_array(store_index+1:,:), order(store_index+1:))
   end if
end subroutine quicksort_real_list

recursive subroutine quicksort_double_list(data_array, order)
! Apply the in-place quicksort algorithm to data_array, returning the sorted
! array and an array of integers mapping the new order of the elements of 
! data_array in terms of the input order.
!
! We use a random pivot point, to prevent the worst-case performance seen when
! trying to sort already sorted arrays.
!
! Arguments
! ---------
! data_array(in/out): array of data to be sorted
! order(in/out): array of integers showing the new order of data_array in terms
! of the original order.
!
   implicit none
! Arguments
   double precision, dimension(:,:), intent(in out)   :: data_array
   integer, dimension(:), intent(in out)              :: order
! Variables
   real                                               :: random
   integer                                            :: pivot
   integer                                            :: i
   integer                                            :: j
   integer                                            :: store_index
   integer                                            :: right
   logical, dimension(:), allocatable                 :: array_lt
   logical, dimension(:), allocatable                 :: array_gt

! Check that the arrays are of the same size, otherwise there's a problem with
! the sort routine.
   if (size(data_array, 1) /= size(order)) then
      stop 'Error in quicksort_double_list, the data_array and order arrays &
           &are not the same size.'
   end if

! If the array is of size smaller than 1 (i.e. 0 or 1), then the array is
! already sorted
   if (size(data_array, 1) > 1) then
   ! Choose a random pivot.
      call random_number(random)
      pivot = int(random * real(size(data_array, 1) - 1)) + 1
   ! Put the pivot on the right
      right = size(data_array, 1)
      call swap(data_array, right, pivot)
      call swap(order, right, pivot)
   ! Allocate two arrays to store the list of less than and greater than
   ! results for comparing each subarray.
      allocate(array_lt(size(data_array, 2)))
      allocate(array_gt(size(data_array, 2)))
   ! Set store_index to 1
      store_index = 1
   ! Now do some swaps to ensure that all the elements to the left of the pivot
   ! are smaller and all those on the right are larger
   ! N.B. data_array(right) now contains the value of the pivot
      do i = 1, right - 1
         array_lt = data_array(i,:) < data_array(right,:) 
         array_gt = data_array(i,:) > data_array(right,:) 
      ! First do the simple check whether or not all the values of array_lt and
      ! array_gt are false (this means the arrays are equal).
         if (all(.not.(array_lt .or. array_gt))) then
            call swap(data_array, i, store_index)
            call swap(order, i, store_index)
            store_index = store_index + 1
         else
         ! We order lists with the most significant value in the smallest index
            do j = 1, size(array_lt)
               if (array_lt(j) .and. .not.(array_gt(j))) then
                  call swap(data_array, i, store_index)
                  call swap(order, i, store_index)
                  store_index = store_index + 1
                  exit
               else if (array_gt(j) .and. .not.(array_lt(j))) then
                  exit
               end if
            end do
         end if
      end do
   ! Put the pivot back in its original place
      call swap(data_array, store_index, right)
      call swap(order, store_index, right)
   ! Sort the arrays to the left and right of the pivot
      call quicksort(data_array(:store_index-1,:), order(:store_index-1))
      call quicksort(data_array(store_index+1:,:), order(store_index+1:))
   end if
end subroutine quicksort_double_list

end module sorts_mod
