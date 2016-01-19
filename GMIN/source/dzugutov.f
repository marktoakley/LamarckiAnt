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
C-----------------------------------------------------------------------*
C
C  Energy and gradient for the Dzugutov potential.
C
C  16.02.2000
C  Sergei Simdyankin <ssim@nada.kth.se>
C
C correspondence with the parameters in the file "potential.f"
C V     = GRAD
C EDZ   = EREAL
C GTEST = GRADT
C STEST 
C
C meaning of the variables
C DIST - distance
C EDZ  - total interaction energy
C VT   - (3*N) vector of interaction energies of each particle
C G    - (N,N) storage for gradient intermediate bits
C V    - (3*N) vector - gradient (force acting on a particle)
! 
C-----------------------------------------------------------------------*
      SUBROUTINE DZPOT(X,V,EDZ,GTEST,STEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST, STEST, DZT
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 DDUMMY, DDDUMMY, NEARD(NATOMS), DZP72
      INTEGER NEAREST(NATOMS)
      COMMON /DZ/ DZT

      DZP72=DZP7**2
      DZT=.FALSE.
      EDZ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
         NEARD(J1)=1.0D100
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2
               ENDIF
               DDUMMY=0.0D0
               IF (DIST.LT.DZP72) THEN
                  call derphi(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,DZP1,DZP2,DZP3,DZP4,DZP5,DZP6,DZP7)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
               G(J2,J1)=DDUMMY
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2
               ENDIF
               IF (DIST.LT.DZP72) THEN
                  call derphi(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,DZP1,DZP2,DZP3,DZP4,DZP5,DZP6,DZP7)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) GOTO 10

      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
            V(J3-2)=0.0D0
            V(J3-1)=0.0D0
            V(J3)=0.0D0
         ELSE
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

10    DO J1=1,NATOMS 
C        IF (VT(J1).EQ.0.0D0) THEN
         IF (NEARD(J1).GT.2.00D0) THEN
            DZT=.TRUE.
C           DIST=SQRT(X(3*(J1-1)+1)**2+X(3*(J1-1)+2)**2+X(3*(J1-1)+3)**2)
C           X(3*(J1-1)+1)=X(3*(J1-1)+1)*(DIST-1.0D0)/DIST
C           X(3*(J1-1)+2)=X(3*(J1-1)+2)*(DIST-1.0D0)/DIST
C           X(3*(J1-1)+3)=X(3*(J1-1)+3)*(DIST-1.0D0)/DIST
C           PRINT*,'J1,VT,DIST=',J1,VT(J1),DIST
C           PRINT*,'moving atom ',J1,' nearer to atom ',NEAREST(J1)
            X(3*(J1-1)+1)=X(3*(NEAREST(J1)-1)+1)+(X(3*(J1-1)+1)-X(3*(NEAREST(J1)-1)+1))/NEARD(J1)
            X(3*(J1-1)+2)=X(3*(NEAREST(J1)-1)+2)+(X(3*(J1-1)+2)-X(3*(NEAREST(J1)-1)+2))/NEARD(J1)
            X(3*(J1-1)+3)=X(3*(NEAREST(J1)-1)+3)+(X(3*(J1-1)+3)-X(3*(NEAREST(J1)-1)+3))/NEARD(J1)
            NEARD(NEAREST(J1))=1.0D0
         ENDIF
      ENDDO

      RETURN
      END

C%%%%%%%
C values of the potential, first, and second derivatives
C if secder = .false. ddphi = 0 on output
C%%%%%%%
C-----------------------------------------------------------------------*
      subroutine derphi(r2,phi,dphi,ddphi,secder,m,A,c,aa,B,d,bb)
      implicit none
      DOUBLE PRECISION r2, phi, dphi, ddphi, m
      logical secder
      integer m2

      DOUBLE PRECISION A, aa, B, bb, c, d 
C     parameter(m=16, A=5.82d0, c=1.10d0, aa=1.87d0, B=1.28d0,  d=0.27d0, bb=1.94d0, m2=8)
C     parameter(m=4,  A=3.00d0, c=0.52d0, aa=1.65d0, B=2.109d0, d=0.55d0, bb=1.94d0, m2=2)
      DOUBLE PRECISION r, V1, dV1, ddV1, V2, dV2, ddV2, DUMMY
      DOUBLE PRECISION expo, exparg

      m2=m/2
      r = sqrt(r2)
C
C  Initialization needed for Sun!
C
      V1=0.0D0
      V2=0.0D0
      dV1=0.0D0
      ddV1=0.0D0
      dV2=0.0D0
      ddV2=0.0D0

      IF (R .LT. BB) THEN
         EXPARG = D/(R-BB)
C        IF (DABS(EXPARG).GT.708.0D0) THEN 
C           EXPO = 0.0D0
C        ELSE
            EXPO = DEXP(EXPARG)
C        ENDIF
         V2 =  B*EXPO
         DV2 = -V2*D/(R-BB)**2
         IF (SECDER) DDV2 = -DV2*D/(R-BB)**2 + 2*V2*D/(R-BB)**3
         IF (R.LT.AA) THEN 
            EXPARG = C/(R-AA)
C           IF (DABS(EXPARG).GT.708.0D0) THEN 
C              EXPO=0.0D0
C           ELSE
               EXPO = DEXP(EXPARG)
C           ENDIF
            DUMMY=1.0D0/R2**M2
            V1=A*(DUMMY-B)*EXPO
            DV1=-A*M*(DUMMY/R)*EXPO - V1*C/(R-AA)**2
            IF (SECDER) THEN
               DDV1 = A*M*(M+1)*(R**(-M-2))*EXPO 
     1             +  A*M*(R**(-M-1))*EXPO*C/((R-AA)**2) 
     2             -  DV1*C/((R-AA)**2) + 2*V1*C/((R-AA)**3) 
            ENDIF
         ENDIF
      ENDIF
  
      PHI=V1+V2
      DPHI=(DV1+DV2)/R
      IF (SECDER) DDPHI = DDV1 + DDV2

      RETURN
      end                       
