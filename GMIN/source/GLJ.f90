!   GMIN: A program for finding global minima
!   Copyright (C) 1999- David J. Wales
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
!  ===================================================================
!  Energy and gradient for Binary Lennard-Jones without a cutoff.
!  Cleaner and faster alternative to LJPSHIFTBINC, which requires a 
!  cutoff and is called when the keyword BLJCLUSTER is used. 
!  Assumed reduced units with SIGMA_AA = EPSILON_AA = 1.
!  Per-atom energy is stored in array VT(NATOMS).
!
!  ds656> 20/06/2013
!  ===================================================================
!
SUBROUTINE GLJ(X,GRAD,POT,GTEST)
  !
  ! USE COMMONS, ONLY : NATOMS,GLJEPS,GLJSIG,VT,NSPECIES,MYUNIT,ZSTAR
  USE COMMONS, ONLY : NATOMS,GLJEPS,GLJSIG,VT,NSPECIES,ZSTAR
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (IN) :: GTEST
  DOUBLE PRECISION, INTENT (IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT (INOUT) :: POT, GRAD(3*NATOMS)
  !
  INTEGER :: I,J13,J23,T1,T2,NT1,NT2,ND1,ND2,J2,J1,J,IJ,J3,J4
  DOUBLE PRECISION :: DIST2, DX(3), IR6, IR12, DUMMY, SIG2, EPS4, EPS24
  DOUBLE PRECISION DIST(NATOMS,NATOMS), SDIST(NATOMS,NATOMS), DOT(NATOMS,NATOMS)

  DOUBLE PRECISION ABBC, ABAC, ACBC, ABI, ACI, BCI, R2(NATOMS,NATOMS), RX, RY, RZ, &
     &                 VEC(NATOMS,NATOMS,3), TEMP, RR2(NATOMS,NATOMS), TDOT, RAB, RRAB, RAC, RRAC, RBC, RRBC, P3

  !
  ! Zero the potential and the gradient
  VT(1:NATOMS) = 0.0D0
  POT = 0.0D0
  IF (GTEST) GRAD(:) = 0.0D0
  NT1=1
  ND1=1
  LT1: DO T1=1,NATOMS
     J13 = 3*(T1-1)
     NT2=1
     ND2=1
542 CONTINUE
     IF (ND1.GT.NSPECIES(NT1)) THEN
        NT1=NT1+1
        ND1=1
        GOTO 542 ! there could be zero atoms of more than one species, though this is silly input!
     ENDIF
     LT2: DO T2=1,T1-1
543     CONTINUE
        IF (ND2.GT.NSPECIES(NT2)) THEN
           NT2=NT2+1
           ND2=1
           GOTO 543 ! there could be zero atoms of more than one species, though this is silly input!
        ENDIF
        SIG2 = GLJSIG(NT2,NT1)*GLJSIG(NT2,NT1)
        EPS4 = 4.0D0*GLJEPS(NT2,NT1)
        J23 = 3*(T2-1)
        DIST2 = 0.0D0
        DO I = 1,3
           DX(I) = X(J13 + I) - X(J23 + I)
           DIST2 = DIST2 + DX(I)*DX(I)
        ENDDO
        IR6 = (SIG2/DIST2)**3
        IR12 = IR6*IR6
        DUMMY = EPS4*(IR12 - IR6)
        VT(T1) = VT(T1) + DUMMY
        VT(T2) = VT(T2) + DUMMY
        POT = POT + DUMMY
!       WRITE(MYUNIT,'(A,2I6,3F15.4)') 'T1,T2,eps,sig,pot=',T1,T2,GLJEPS(NT2,NT1),GLJSIG(NT2,NT1),POT
        IF (GTEST) THEN ! Calculate gradient 
           EPS24 = 24.0D0*GLJEPS(NT2,NT1)
           DUMMY = EPS24*(IR6 - 2.0D0*IR12)/DIST2
           DO I = 1,3
              GRAD(J13+I) = GRAD(J13+I) + DUMMY*DX(I)
              GRAD(J23+I) = GRAD(J23+I) - DUMMY*DX(I)
           ENDDO
        ENDIF
        ND2=ND2+1
     ENDDO LT2
     ND1=ND1+1
  ENDDO LT1

IF (ZSTAR.EQ.0.0D0) RETURN
!  Add an Axilrod-Teller three-body term.
      DO J1=1,NATOMS
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO J2=J1+1,NATOMS
            SDIST(J1,J2)=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 + &
     &                   ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 + &
     &                   ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
            DIST(J1,J2)=DSQRT(SDIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
         ENDDO
      ENDDO

      DO J1=1,NATOMS
         DO J2=J1,NATOMS
            DOT(J2,J1)=X(3*(J2-1)+1)*X(3*(J1-1)+1) &
     &                +X(3*(J2-1)+2)*X(3*(J1-1)+2) &
     &                +X(3*(J2-1)+3)*X(3*(J1-1)+3)
            DOT(J1,J2)=DOT(J2,J1)
         ENDDO
      ENDDO

      P3=0.0D0
      DO I=1,NATOMS
         DO J=1,NATOMS
            IF (I.NE.J) THEN
               DO IJ=1,NATOMS
                  IF ((IJ.NE.I).AND.(IJ.NE.J)) THEN
                     P3=P3+(1.0-3.0* &
     &                   (DOT(J,I)+DOT(IJ,J)-DOT(IJ,I)-DOT(J,J)) &
     &                  *(DOT(J,I)+DOT(IJ,IJ)-DOT(IJ,J)-DOT(IJ,I)) &
     &                  *(DOT(I,I)+DOT(IJ,J)-DOT(J,I)-DOT(IJ,I)) &
     &                  / (DIST(I,J)*DIST(I,IJ)*DIST(J,IJ))**2 )  &
     &                  / (DIST(I,J)*DIST(I,IJ)*DIST(J,IJ))**3
                   ENDIF
                ENDDO
             ENDIF
         ENDDO
      ENDDO
      P3=P3*ZSTAR/6.0
      POT=POT+P3
      IF (GTEST) THEN
      DO J1=1,NATOMS
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO J2=J1+1,NATOMS
            R2(J2,J1)= &
     &                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 + &
     &                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 + &
     &                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
            RR2(J2,J1)=1.0D0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            R2(J1,J2)=R2(J2,J1)
            RR2(J1,J2)=RR2(J2,J1)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
         ENDDO
      ENDDO

      DO J1=1,NATOMS
         DO J2=1,3
            TEMP=0.0D0
            DO J3=1,NATOMS
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABI=VEC(J3,J1,J2)
               RX=VEC(J3,J1,1)
               RY=VEC(J3,J1,2)
               RZ=VEC(J3,J1,3)
               DO J4=J3+1,NATOMS
                  ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                  ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3)
                  ACBC=VEC(J4,J1,1)*VEC(J4,J3,1) &
     &                +VEC(J4,J1,2)*VEC(J4,J3,2) &
     &                +VEC(J4,J1,3)*VEC(J4,J3,3)
                  TDOT=ABAC*ACBC*ABBC
                  BCI=VEC(J4,J3,J2)
                  ACI=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)

                  TEMP=TEMP+ DSQRT(RRAB*RRAC*RRBC)**5 *          ( &
     &    3*(ABAC*ACBC*BCI + ABBC*(ACBC*(ABI + ACI) + ABAC*BCI) +  &
     &                            (ACI*RAB + ABI*RAC)*RBC) -  &
     &                         15*(ABI*RRAB + ACI*RRAC)*TDOT   )

               ENDDO
            ENDDO
            GRAD(3*(J1-1)+J2)=GRAD(3*(J1-1)+J2)+ZSTAR*TEMP
         ENDDO
      ENDDO
      ENDIF



  RETURN
END SUBROUTINE GLJ
