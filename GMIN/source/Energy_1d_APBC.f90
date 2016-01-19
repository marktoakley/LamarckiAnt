! 
! One-dimensional periodic XY model. ch558
!
SUBROUTINE Energy_1d_APBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE COMMONS, ONLY : NONEDAPBC, NATOMS, XYPHI
  IMPLICIT NONE
  INTEGER N
  DOUBLE PRECISION, dimension(3*NATOMS) :: theta, GRAD
  DOUBLE PRECISION :: Energy
  LOGICAL GTEST,SECT

  n=NONEDAPBC

 Energy = sum(cos(xyphi(1:n-1)+theta(2:n)-theta(1:n-1)))  
 Energy = Energy + cos(xyphi(n) - theta(1)-theta(n))
 Energy = 1 - (Energy/n)

 IF (.NOT.GTEST) RETURN

grad(1)=-sin(xyphi(1) + theta(2) - theta(1)) - sin( xyphi(n) - theta(1) - theta(n))

grad(n)= sin(xyphi(n-1) + theta(n) - theta(n-1)) - sin(xyphi(n) - theta(1) -theta(n))

grad(2:(n-1))= sin(xyphi(1:(n-2)) + theta(2:(n-1)) - theta(1:(n-2))) - sin( xyphi(2:(n-1)) + theta(3:n) - theta(2:n-1))

END SUBROUTINE ENERGY_1d_APBC

