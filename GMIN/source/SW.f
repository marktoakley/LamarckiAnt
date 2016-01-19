C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C*************************************************************************
C
C  Subroutine SW2 calculates the energy and Cartesian gradient
C  for the SW Si potential.
C
C*************************************************************************
C
      SUBROUTINE SWTWO(X, V, ENERGY, GTEST)
      USE commons
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST
      DOUBLE PRECISION X(3*NATOMS), ENERGY, 
     1                 V(3*NATOMS), R(NATOMS,NATOMS), R4T, LVEC(NATOMS,NATOMS,3),
     2                 DUMMY,
     3                 ER(NATOMS,NATOMS), R2(NATOMS,NATOMS)
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)

      N=NATOMS
      DO J1=1,NATOMS
         VT(J1)=0.0D0
      ENDDO
C
C  Deal with any atoms that have left the box.
C
      IF (PERIODIC.AND.(.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C 
C  Store distance matrices.
C
      ENERGY=0.0D0
         DO J1=1, N
            LVEC(J1,J1,1)=0.0D0
            LVEC(J1,J1,2)=0.0D0
            LVEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               LVEC(J2,J1,1)=X(J3+1)-X(J4+1)
               LVEC(J2,J1,2)=X(J3+2)-X(J4+2)
               LVEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (PERIODIC.AND.(.NOT.FIXIMAGE)) THEN
                  ANV(J2,J1,1)=NINT(LVEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(LVEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(LVEC(J2,J1,3)/BOXLZ)
               ENDIF
               LVEC(J2,J1,1)=LVEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               LVEC(J2,J1,2)=LVEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               LVEC(J2,J1,3)=LVEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               LVEC(J1,J2,1)=-LVEC(J2,J1,1)
               LVEC(J1,J2,2)=-LVEC(J2,J1,2)
               LVEC(J1,J2,3)=-LVEC(J2,J1,3)
            ENDDO
         ENDDO 
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R(J1,J1)=0.0D0
            ER(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=LVEC(J2,J1,1)**2+LVEC(J2,J1,2)**2+LVEC(J2,J1,3)**2
               R(J2,J1)=SQRT(R2(J2,J1)) 
               IF (R(J2,J1).LT.SMALLA) THEN
                  R2(J2,J1)=1.0D0/R2(J2,J1)
                  R4T=R2(J2,J1)**2 
                  ER(J2,J1)=EXP(1.0D0/(-SMALLA + R(J2,J1)))
                  DUMMY=ER(J2,J1)*(-1.0D0 + BIGB*R4T)
                  ENERGY=ENERGY+DUMMY
                  VT(J1)=VT(J1)+DUMMY*BIGA
                  VT(J2)=VT(J2)+DUMMY*BIGA
               ENDIF
               R(J1,J2)=R(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               ER(J1,J2)=ER(J2,J1)
            ENDDO 
         ENDDO

      ENERGY=BIGA*ENERGY

      IF (.NOT.GTEST) RETURN
C
      DO J1=1,3*N
         V(J1)=0.0D0
      ENDDO
      DO J1=1,N
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               DUMMY=(BIGA*ER(J2,J1)*(-4*SMALLA**2*BIGB + (-1.0D0 + 8.0D0*SMALLA)*BIGB*R(J2,J1) - 
     1                4*BIGB/R2(J2,J1) + R(J2,J1)**5))*R2(J2,J1)**3/(SMALLA - R(J2,J1))**2
               V(3*(J1-1)+1)=V(3*(J1-1)+1)+DUMMY*LVEC(J2,J1,1)
               V(3*(J1-1)+2)=V(3*(J1-1)+2)+DUMMY*LVEC(J2,J1,2)
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+DUMMY*LVEC(J2,J1,3)
               V(3*(J2-1)+1)=V(3*(J2-1)+1)-DUMMY*LVEC(J2,J1,1)
               V(3*(J2-1)+2)=V(3*(J2-1)+2)-DUMMY*LVEC(J2,J1,2)
               V(3*(J2-1)+3)=V(3*(J2-1)+3)-DUMMY*LVEC(J2,J1,3)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine SW3 calculates the energy and Cartesian gradient 
C  derivative matrix analytically for the SW Si potential.
C
C*************************************************************************
C
      SUBROUTINE SWTHREE(V, ENERGY, GTEST)
      USE commons
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, K1, K2, K3
      LOGICAL GTEST
      DOUBLE PRECISION ENERGY, ER2131, E3,
     1                 V(3*NATOMS), R(NATOMS,NATOMS), R122, R132, R12S, R13S,
     2                 r12DOTr13PT, r12DOTr13, ERLDOT,
     3                 ER(NATOMS,NATOMS), R2(NATOMS,NATOMS), THIRD,
     4                 ER21, ER31, R12, R13, DOT2, DUMMY1, DUMMY2, DUMMY3,
     5                 LVEC(NATOMS,NATOMS,3), V211, V212, V213, V311, V312, V313, V211b, V212b, V213b
      DOUBLE PRECISION XLAMBDA, LGAMMA, SMALLA
      PARAMETER (LGAMMA=1.2D0, SMALLA=1.8D0, THIRD=0.33333333333333333333D0, XLAMBDA=21.0D0)
C
C  In the Mousseau version XLAMBDA=21.0D0*1.5D0
C
C 
C  Store things
C
      N=NATOMS
      DO J1=1,N
         ER(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               ER(J2,J1)=EXP(LGAMMA/(-SMALLA + R(J2,J1)))
            ELSE
               ER(J2,J1)=0.0D0
            ENDIF
            ER(J1,J2)=ER(J2,J1)
         ENDDO
      ENDDO
      E3=0.0D0

      DO J1=1,N
         K1=3*(J1-1)
         DO J2=1,N
            IF ((R(J2,J1).LT.SMALLA).AND.(J2.NE.J1)) THEN
               K2=3*(J2-1)
               ER21=ER(J2,J1)
               R12=R(J2,J1)
               r12s=(r12 - smalla)**2
               V211=LVEC(J2,J1,1)
               V211b=V211/(r12*r12s)
               V212=LVEC(J2,J1,2)
               V212b=V212/(r12*r12s)
               V213=LVEC(J2,J1,3)
               V213b=V213/(r12*r12s)
               
               DO J3=1,N
                  IF (((R(J3,J1).LT.SMALLA).AND.(J3.NE.J1)).AND.(J3.NE.J2)) THEN
                     K3=3*(J3-1)
                     ER31=ER(J3,J1)
                     ER2131=ER21*ER31
                     R13=R(J3,J1)
                     V311=LVEC(J3,J1,1)
                     V312=LVEC(J3,J1,2)
                     V313=LVEC(J3,J1,3)
                     r12DOTr13=(V211*V311+V212*V312+V213*V313)/(r12*r13)
                     r12DOTr13PT=(r12DOTr13+third)
                     DOT2=r12DOTr13PT**2

                     DUMMY1=ER2131*DOT2
                     E3=E3+DUMMY1
                     VT(J1)=VT(J1)+DUMMY1*XLAMBDA/2.0D0
                     VT(J2)=VT(J2)+DUMMY1*XLAMBDA/2.0D0
                     VT(J3)=VT(J3)+DUMMY1*XLAMBDA/2.0D0

                     IF (GTEST) THEN
                        R122=R12**2
                        R132=R13**2
                        r13s=(r13 - smalla)**2
                        ERLDOT=ER2131*XLAMBDA*r12DOTr13PT
                        DUMMY1=2*r12DOTr13+(lgamma*r12*r12DOTr13PT)/r12s
                        DUMMY2=2*r12/r13
                        DUMMY3=2*r12DOTr13+(lgamma*r12DOTr13PT*r13)/r13s

                        V(K1+1)=V(K1+1)+ (ERLDOT*((2*r13*(r12 - r12DOTr13*r13)*V211 + 
     &                           2*r12*(-(r12*r12DOTr13) + r13)*V311)/(r122*r132) + lgamma*r12DOTr13PT*(-V211b - 
     &                           V311/(r13*r13s))))/2.0D0
                        V(K1+2)=V(K1+2)+ (ERLDOT*((2*r13*(r12 - r12DOTr13*r13)*V212 + 
     &                           2*r12*(-(r12*r12DOTr13) + r13)*V312)/(r122*r132) + lgamma*r12DOTr13PT*(-V212b - 
     &                           V312/(r13*r13s))))/2.0D0
                        V(K1+3)=V(K1+3)+ (ERLDOT*((2*r13*(r12 - r12DOTr13*r13)*V213 + 
     &                           2*r12*(-(r12*r12DOTr13) + r13)*V313)/(r122*r132) + lgamma*r12DOTr13PT*(-V213b - 
     &                           V313/(r13*r13s))))/2.0D0
                        V(K2+1)=V(K2+1)+ ERLDOT*(DUMMY1*V211-DUMMY2*V311)/(2.*r122)
                        V(K2+2)=V(K2+2)+ ERLDOT*(DUMMY1*V212-DUMMY2*V312)/(2.*r122)
                        V(K2+3)=V(K2+3)+ ERLDOT*(DUMMY1*V213-DUMMY2*V313)/(2.*r122)
                        V(K3+1)=V(K3+1)+ ERLDOT*((-2*r13*V211)/r12+DUMMY3*V311)/(2.*r132)
                        V(K3+2)=V(K3+2)+ ERLDOT*((-2*r13*V212)/r12+DUMMY3*V312)/(2.*r132)
                        V(K3+3)=V(K3+3)+ ERLDOT*((-2*r13*V213)/r12+DUMMY3*V313)/(2.*r132)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ENERGY=ENERGY+E3*XLAMBDA/2.0D0

      RETURN
      END
