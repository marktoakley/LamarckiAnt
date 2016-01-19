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
!  ===============================================================
!  Multicomponent Lennard-Jones without a cutoff.
!  Assumed reduced units with SIGMA_AA = EPSILON_AA = 1.
!  Per-atom energy is stored in array VT(NATOMS).
!
!  ds656> 6/12/2014
!  ===============================================================
!
SUBROUTINE MLJ(X,GRAD,POT,GTEST)
  !
  USE COMMONS, ONLY : NATOMS,VT,ATOMLISTS, NSPECIES
  USE POT_PARAMS, ONLY : MLJ_EPS, MLJ_SIG
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (IN) :: GTEST
  DOUBLE PRECISION, INTENT (IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT (INOUT) :: POT, GRAD(3*NATOMS)
  !
  LOGICAL :: SELF
  INTEGER :: I,J1,J13,J2,J23,T1,T2,G1,G2,G2START,&
       K1,K2,LI1,LI2,UI1,UI2
  DOUBLE PRECISION :: DIST2, DX(3), IR6, IR12, DUMMY, SIG2, EPS4, EPS24
  !
  ! Zero the potential and the gradient
  VT(1:NATOMS) = 0.0D0
  POT = 0.0D0
  IF(GTEST) GRAD(:) = 0.0D0
  !
  ! Double loop over atom types
  LT1: DO T1=1,NSPECIES(0)
     LT2: DO T2=1,T1 
        !
        ! This can be precomputed elsewhere...
        SIG2 = MLJ_SIG(T1,T2)*MLJ_SIG(T1,T2)
        EPS4 = 4.0D0*MLJ_EPS(T1,T2)
        IF(GTEST) EPS24 = 24.0D0*MLJ_EPS(T1,T2)
        !
        ! Double loop over mobile and frozen atom groups
        LG1: DO G1 = 1,2 ! 1 -> mobile; 2 -> frozen
           !
           IF(T1 == T2) THEN
              G2START = G1
           ELSE
              G2START = 1
           ENDIF
           !
           LG2: DO G2 = G2START,2 ! include frozen-frozen
              !
              LI1=1
              UI1=ATOMLISTS(T1,G1,0)
              !
              IF(T1==T2.AND.G1==G2) THEN
                 SELF = .TRUE.
                 UI1 = UI1 - 1 ! exclude self-interaction
              ELSE
                 SELF = .FALSE.
              ENDIF
              !
              LA1: DO K1=LI1,UI1
                 !
                 J1 = ATOMLISTS(T1,G1,K1) ! actual atom 1 index
                 J13 = 3*(J1-1)
                 !
                 ! Atom loop bounds depend on type and group.
                 ! Can this IF block be moved outside of loop LA1?
                 IF(SELF) THEN
                    LI2=K1+1
                    UI2=UI1+1
                 ELSE
                    LI2=1
                    UI2=ATOMLISTS(T2,G2,0)
                 ENDIF
                 !
                 LA2: DO K2=LI2,UI2
                    !
                    J2 = ATOMLISTS(T2,G2,K2) ! actual atom 2 index
                    J23 = 3*(J2-1)
                    !
                    DIST2 = 0.0D0
                    DO I = 1,3
                       DX(I) = X(J13 + I) - X(J23 + I)
                       DIST2 = DIST2 + DX(I)*DX(I)
                    ENDDO
                    !
                    ! ---- Potential-specific bit -------------
                    IR6 = (SIG2/DIST2)**3
                    IR12 = IR6*IR6
                    DUMMY = EPS4*(IR12 - IR6)
                    !WRITE(*,*) "DUMMY FOR POT", DUMMY
                    ! -----------------------------------------
                    !
                    VT(J1) = VT(J1) + DUMMY
                    VT(J2) = VT(J2) + DUMMY
                    POT = POT + DUMMY
                    !
                    IF(GTEST) THEN ! Calculate gradient for groups 1 and 2
                       DUMMY = EPS24*(IR6 - 2.0D0*IR12)/DIST2
                       !WRITE(*,*) "DUMMY FOR GRAD", DUMMY
                       DO I = 1,3
                          GRAD(J13+I) = GRAD(J13+I) + DUMMY*DX(I)
                          GRAD(J23+I) = GRAD(J23+I) - DUMMY*DX(I)
                       ENDDO
                    ENDIF
                    !
                 ENDDO LA2 ! Double loop over atoms
              ENDDO LA1
              !
           ENDDO LG2 ! Double loop over groups
        ENDDO LG1
        !
     ENDDO LT2 ! Double loop over types
  ENDDO LT1
  !
  IF(GTEST) THEN
     DO T1=1,NSPECIES(0) ! Atom types
        DO J2=1,ATOMLISTS(T1,2,0) ! span group 2 only
           J1=ATOMLISTS(T1,2,J2)
           J13=3*(J1-1)
           GRAD(J13+1) = 0.0D0 ! Reset to zero
           GRAD(J13+2) = 0.0D0
           GRAD(J13+3) = 0.0D0
        ENDDO
     ENDDO
     !
     ! Test per-atom energies and gradients:
     ! DO T1=1,NSPECIES(0)
     !    DO G1=1,2
     !       WRITE(*,'(A,3(1X,I5))') "BLJ_CLUST> T,G,NAT(T,G)=", &
     !            T1, G1, ATOMLISTS(T1,G1,0)
     !       DO J2=1,ATOMLISTS(T1,G1,0)
     !          J1=ATOMLISTS(T1,G1,J2)
     !          J13=3*(J1-1)
     !          WRITE(*,'(A,1X, I4, 4(1X,F12.6))') "BLJ_CLUST>", &
     !               J1, VT(J1), GRAD(J13+1), GRAD(J13+2), GRAD(J13+3)
     !       ENDDO
     !    ENDDO
     ! ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MLJ
