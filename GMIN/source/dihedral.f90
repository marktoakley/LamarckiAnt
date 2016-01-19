module dihedral_mod

implicit none

contains

function cross_product(vec1, vec2) result(output)
! Calculate the cross product of vec1 and vec2
   use prec, only: real64
   implicit none
! Arguments
   real(real64), intent(in)   :: vec1(3), vec2(3)
   real(real64)               :: output(3)

   output(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
   output(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
   output(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

end function cross_product

subroutine dihedral_with_grad(coords, angle, grad)
! Calculate the dihedral angle and the gradient of the angle wrt each x, y and z
! coordinate. The gradient is only calculated if an appropriate argument is
! passed.
   use prec, only: int32, real64
   implicit none
! Arguments
   real(real64), intent(in)            :: coords(12)
   real(real64), intent(out)           :: angle
   real(real64), intent(out), optional :: grad(12)
! Variables
   real(real64)                        :: b1(3), b2(3), b3(3)
   real(real64)                        :: b1xb2(3), b2xb3(3)
   real(real64)                        :: ddot, dnrm2
   integer(int32)                      :: i

! Calculate the values of variables as in function dihedral
   b1(1:3) = coords(4:6)   - coords(1:3)
   b2(1:3) = coords(7:9)   - coords(4:6)
   b3(1:3) = coords(10:12) - coords(7:9)

! Calculate cross products of b1 with b2 and b2 with b3
   b1xb2 = cross_product(b1(1:3), b2(1:3))
   b2xb3 = cross_product(b2(1:3), b3(1:3))

! Dihedral is given by:
!
! phi = atan2( ([b1 x b2] x [b2 x b3]) . (b2/|b2|),
!               [b1 x b2] . [b2 x b3] )
!
   angle = atan2( ddot(3, cross_product(b1xb2, b2xb3), 1, b2/dnrm2(3, b2, 1), 1), &
                  ddot(3, b1xb2, 1, b2xb3, 1))

! If the grad argument is provided, calculate the gradient.
! We follow the form given by:
!
! van Schaik, R.C.; Berendsen, H.J.C.; Torda, A.E.; van Gunsteren W.F.
! J. Mol. Biol. 1993 234(3) pp751-762
!
! doi: http://dx.doi.org/10.1006/jmbi.1993.1624
!
   if (present(grad)) then
! dphi/dr_i
      grad(1:3)   = -1.0 * (dnrm2(3, b2, 1) / dnrm2(3, b1xb2, 1)**2) * b1xb2(1:3)
! dphi/dr_l
      grad(10:12) =        (dnrm2(3, b2, 1) / dnrm2(3, b2xb3, 1)**2) * b2xb3(1:3)
! dphi/dr_j
      grad(4:6)   = grad(10:12) * ddot(3, b3, 1, b2, 1) / dnrm2(3, b2, 1)**2 &
                  - grad(1:3)   * ddot(3, b1, 1, b2, 1) / dnrm2(3, b2, 1)**2 &
                  - grad(1:3)   
! dphi/dr_k
      grad(7:9)   = grad(1:3)   * ddot(3, b1, 1, b2, 1) / dnrm2(3, b2, 1)**2 &
                  - grad(10:12) * ddot(3, b3, 1, b2, 1) / dnrm2(3, b2, 1)**2 &
                  - grad(10:12)   
   end if

end subroutine dihedral_with_grad

subroutine test_dihedral_grad(coords, delta, tol, passed)
! Test for the dihedral gradient which uses the numerical gradient.
   use prec, only: int32, real64
   implicit none
! Arguments
   real(real64), intent(in)   :: coords(12)
   real(real64), intent(in)   :: delta
   real(real64), intent(in)   :: tol
   logical, intent(out)       :: passed
! Variables
   real(real64)               :: dih_plus, dih_minus, dih_zero
   real(real64)               :: tmp_coords(12)
   real(real64)               :: analytic_grad(12), numeric_grad(12)
   integer(int32)             :: i

! Copy coordinates to a temporary array to save them from being changed.
   tmp_coords(1:12) = coords(1:12)

   do i = 1, 12
! Calculate the angle at +delta
      tmp_coords(i) = tmp_coords(i) + delta
      call dihedral_with_grad(tmp_coords, dih_plus)
! Calculate the angle at -delta
      tmp_coords(i) = tmp_coords(i) - 2.0 * delta
      call dihedral_with_grad(tmp_coords, dih_minus)
! Reset tmp_coords
      tmp_coords(i) = coords(i)
      numeric_grad(i) = (dih_plus - dih_minus) / (2.0 * delta)
   end do

! Calculate analytic gradient
   call dihedral_with_grad(coords, dih_zero, analytic_grad)

! Compare the gradients
   passed = all(abs(analytic_grad - numeric_grad) < tol)
   if (.not. passed) then
      do i = 1, 12
         write(*, '(3E10.3)') analytic_grad(i), numeric_grad(i), analytic_grad(i) - numeric_grad(i)
      end do
   else
      write(*, *) "Dihedral test passed!"
   end if

end subroutine test_dihedral_grad

end module dihedral_mod
