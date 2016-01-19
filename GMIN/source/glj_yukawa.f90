!   GMIN: A program for finding global minima                           
!   Copyright (C) 1999-2006 David J. Wales                              
!   This file is part of GMIN.                                          
!                                                                       
!   GMIN is free software; you can redistribute it and/or modify        
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or   
!   (at your option) any later version.                                 
!                                                                       
!   GMIN is distributed in the hope that it will be useful,             
!   but WITHOUT ANY WARRANTY; without even the implied warranty of      
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
!   GNU General Public License for more details.                        
!                                                                       
!   You should have received a copy of the GNU General Public License   
!   along with this program; if not, write to the Free Software         
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!                                                                       
!-----------------------------------------------------------------------
! 
! Energy and gradient for generalised Lennard-Jones and Yukawa potentials.
! [See Mossa et al. Langmuir 20, 10756 (2004).]
!
! 28/5/2013 ds656
!
MODULE GLJYMOD
  !
  INTEGER :: GLJ_EXP
  DOUBLE PRECISION :: YUK_A, YUK_XI
  !
END MODULE GLJYMOD
!
SUBROUTINE GLJYPOT(X,GRAD,POT,GTEST)
  !
  USE COMMONS, ONLY : NATOMS, VT, DEBUG
  USE GLJYMOD
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (IN) :: GTEST
  DOUBLE PRECISION, INTENT (IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT (OUT) :: POT, GRAD(3*NATOMS)
  !
  INTEGER :: J1, J2, J13, J23, I
  DOUBLE PRECISION, PARAMETER :: EPS=1.0D0, SIG=1.0D0
  DOUBLE PRECISION :: DIST, DX(3), GLJ_ATT, GLJ_REP, YUK_RAT, YUK_EXP
  DOUBLE PRECISION :: POTIJ, EPS4, FACT, DERIV
  !
  EPS4 = 4.0D0*EPS
  FACT = EPS4*GLJ_EXP
  !
  POT=0.0D0
  VT(1:NATOMS)=0.0D0
  !
  IF (GTEST) THEN ! With gradient
     !
     GRAD(:) = 0.0D0
     !
     DO J1=1,NATOMS-1
        J13=3*(J1-1)
        DO J2=J1+1,NATOMS
           J23=3*(J2-1)
           !
           DIST = 0.0D0
           DO I = 1,3
              DX(I) = X(J13 + I) - X(J23 + I)
              DIST = DIST + DX(I)*DX(I)
           ENDDO
           DIST = DSQRT(DIST)
           !
           GLJ_ATT = (SIG/DIST)**GLJ_EXP
           GLJ_REP = GLJ_ATT*GLJ_ATT
           YUK_RAT = DIST /  YUK_XI
           YUK_EXP = YUK_A*DEXP(-YUK_RAT)
           !
           POTIJ = EPS4*(GLJ_REP - GLJ_ATT) + YUK_EXP/YUK_RAT
           !
           VT(J1) = VT(J1) + POTIJ
           VT(J2) = VT(J2) + POTIJ
           POT = POT + POTIJ
           !
           DERIV = - ( FACT*(2.0D0*GLJ_REP - GLJ_ATT) + &
                YUK_EXP*(1.0D0 + 1.0D0/YUK_RAT) ) / DIST 
           !
           IF(DEBUG) WRITE(*,*)'gljy> DVDR = ', DERIV
           DO I = 1,3
              DX(I) = DX(I) / DIST
              GRAD(J13+I) = GRAD(J13+I) + DERIV*DX(I)
              GRAD(J23+I) = GRAD(J23+I) - DERIV*DX(I)
           ENDDO
           !
        ENDDO
     ENDDO
     !
  ELSE ! Without gradient
     !
     DO J1=1,NATOMS-1
        J13=3*(J1-1)
        DO J2=J1+1,NATOMS
           J23=3*(J2-1)
           !
           DIST = 0.0D0
           DO I = 1,3
              DX(I) = X(J13 + I) - X(J23 + I)
              DIST = DIST + DX(I)*DX(I)
           ENDDO
           DIST = DSQRT(DIST)
           !
           GLJ_ATT = (SIG/DIST)**GLJ_EXP
           GLJ_REP = GLJ_ATT*GLJ_ATT
           YUK_RAT = DIST /  YUK_XI
           YUK_EXP = YUK_A*DEXP(-YUK_RAT)
           !
           POTIJ = EPS4*(GLJ_REP - GLJ_ATT) + YUK_EXP/YUK_RAT
           !
           VT(J1) = VT(J1) + POTIJ
           VT(J2) = VT(J2) + POTIJ
           POT = POT + POTIJ
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE GLJYPOT
