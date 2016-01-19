!
! XY model with Periodic Boundary conditions. ch558
!

SUBROUTINE Energy_1d_PBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE COMMONS, ONLY : NONEDAPBC , NATOMS , XYPHI
  IMPLICIT NONE
  INTEGER N
  DOUBLE PRECISION , DIMENSION(3*NATOMS) :: THETA, GRAD
  DOUBLE PRECISION :: ENERGY
  LOGICAL GTEST, SECT

  N = NONEDAPBC

  theta(n)=0
  energy = sum( cos(xyphi(1:n-1) + theta(2:n) - theta(1:n-1)))
  energy = energy + cos( xyphi(n) + theta(1) )
  energy = 1 - (energy/n)

  IF (.NOT.GTEST) RETURN

  grad(1) = -sin(xyphi(1) + theta(2) - theta(1)) + sin( xyphi(n) + theta(1) - theta(n))
  
  grad(2:(n-1)) = sin( xyphi(1:(n-2)) + theta(2:(n-1)) - theta(1:(n-2))) - sin( xyphi(2:(n-1)) + theta(3:n) - theta(2:(n-1)))

  grad(n)=0

END SUBROUTINE Energy_1d_PBC
