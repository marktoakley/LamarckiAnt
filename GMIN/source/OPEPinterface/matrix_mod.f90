module matrix_mod

contains

! This subroutine computes the elements of matrix mtr_c, here mtr_c is a 3*3 matrix. 
! mtr_c(1,1) is the sum of xa(i)*xb(i), mtr_c(1,2) is the sum of ya(i)*xb(i), mtr_c(1,3) 
! is the sum of za(i)*xb(i) and so on 
subroutine matrix_c(posa,posb,mtr_c)
  use defs
  implicit none

  double precision, dimension(vecsize), intent(in), target :: posa, posb

  double precision, dimension(3,3), intent(out) :: mtr_c
  double precision, dimension(:), pointer :: xa, ya, za, xb, yb, zb
  integer :: i,j,k

  ! We first set-up pointers for the x, y, z components for posref and pos
  xa => posa(1:3*NATOMS:3)
  ya => posa(2:3*NATOMS:3)
  za => posa(3:3*NATOMS:3)

  xb => posb(1:3*NATOMS:3)
  yb => posb(2:3*NATOMS:3)
  zb => posb(3:3*NATOMS:3)


  do j=1,3
     do k=1,3
        mtr_c(j,k)=0.d0
     end do
  end do
    
! we get the matrix C(mtr_c), the matrix element mtr_c(k,j)=sum(i)((x_ji*y_ki))   
  
  do i=1, NATOMS
     mtr_c(1,1) = mtr_c(1,1)+xa(i)*xb(i)
     mtr_c(1,2) = mtr_c(1,2)+ya(i)*xb(i)
     mtr_c(1,3) = mtr_c(1,3)+za(i)*xb(i)

     mtr_c(2,1) = mtr_c(2,1)+xa(i)*yb(i)
     mtr_c(2,2) = mtr_c(2,2)+ya(i)*yb(i)
     mtr_c(2,3) = mtr_c(2,3)+za(i)*yb(i)  

     mtr_c(3,1) = mtr_c(3,1)+xa(i)*zb(i)
     mtr_c(3,2) = mtr_c(3,2)+ya(i)*zb(i)
     mtr_c(3,3) = mtr_c(3,3)+za(i)*zb(i)
  end do

end subroutine

end module






