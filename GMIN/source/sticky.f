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
C  Energy and gradient for sticky patches model.
C
      SUBROUTINE STICKY(X,V,ESTICKY,GTEST,SECT)
      USE commons

      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, NAT2, J3, J4, J3MIN, J4MIN
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, 
     1                 ESTICKY, DUMMY, RALPHA12, RALPHA22, SA1, SA2,
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, C3A1DOT1, C3A2DOT2,
     4                 GX1, GY1, GZ1, C2A2, C2A1, DOT1, DOT2, DIST, 
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, RSQ, R123, R1212, R126, 
     7           SIGMA22, R128, R1214, R12,
     8           OMG1,OMG2,CSOMG1,CSOMG2,PX,PY,PZ,CSOMGM1,CSOMGM2,
     9           FOMG1,FOMG2,
     +           OMG1B,OMG2B,CSOMG1B,CSOMG2B,DUMMYB


      SIGMA22=STICKYSIG**2
      NAT2=NATOMS/2
      ESTICKY=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO

C     DO J1=1,NRBSITES
C        PRINT '(A,I5,3G20.10)','patch ',J1, site(J1,1),site(J1,2),site(J1,3)
C     ENDDO

C  Potential energy first.
C
      DO J1=1,NAT2-1
         X1=X(3*(J1-1)+1)
         Y1=X(3*(J1-1)+2)
         Z1=X(3*(J1-1)+3)
         L1=X(3*(NAT2+J1-1)+1)
         M1=X(3*(NAT2+J1-1)+2)
         N1=X(3*(NAT2+J1-1)+3)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
!        RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-20)
         CA1=COS(ALPHA1)
         SA1=SIN(ALPHA1)
!        C2A1=CA1
!	another change here! there was a mistake in the limit case of C3A2
         IF (ALPHA1.LT.0.0010D0) THEN
C           C3A1=-1/2+ALPHA1**2/24  ! bug !! DJW
            C3A1=-0.5D0+ALPHA1**2/24
            S1=1.0D0-ALPHA1**2/6
c            C3A1=-1/2+ALPHA1**2/24-ALPHA1**4/720+ALPHA1**6/40320 ! the 1/2 is wrong here as well!
c            S1=1.0D0-ALPHA1**2/6+ALPHA1**4/120-ALPHA1**6/5040
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF
C        WRITE(*,'(A,6F15.5)') 'ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1=',ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1
         DO J2=J1+1,NAT2
            X2=X(3*(J2-1)+1)
            Y2=X(3*(J2-1)+2)
            Z2=X(3*(J2-1)+3)
            L2=X(3*(NAT2+J2-1)+1)
            M2=X(3*(NAT2+J2-1)+2)
            N2=X(3*(NAT2+J2-1)+3)
            L22=L2**2
            M22=M2**2
            N22=N2**2
            ALPHA2=SQRT(L22+M22+N22)
!           RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-20)
            CA2=COS(ALPHA2)
            SA2=SIN(ALPHA2)
!           C2A2=CA2
!	another change here! there was a mistake in the limit case of C3A2
            IF (ALPHA2.LT.0.00001D0) THEN
C              C3A2=-1/2+ALPHA2**2/24  ! bug !!
               C3A2=-0.5D0+ALPHA2**2/24
               S2=1.0D0-ALPHA2**2/6
c               C3A2=-1/2+ALPHA2**2/24-ALPHA2**4/720+ALPHA2**6/40320 ! 1/2 is a bug here as well
c               S2=1.0D0-ALPHA2**2/6+ALPHA2**4/120-ALPHA2**6/5040
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
            RSQ=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
            R12=SQRT(RSQ)
            R123=R12*RSQ
            R126=RSQ**3
            R1212=R126*R126
            R128=R126*RSQ
            R1214=R1212*RSQ
C
C  Identify closest patches and load into P1X, P1Y, P1Z, P2X, P2Y, P2Z
C
cC
c            DMIN=1.0D100
c            DO J3=1,NRBSITES
c               P1X=SITE(J3,1)
c               P1Y=SITE(J3,2)
c               P1Z=SITE(J3,3)
c               DO J4=1,NRBSITES
c                  P2X=SITE(J4,1)
c                  P2Y=SITE(J4,2)
c                  P2Z=SITE(J4,3)
c                  DIST=
c     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+
c     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 +
c     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+
c     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 +
c     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+
c     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
c                  IF (DIST.LT.DMIN) THEN
c                     J3MIN=J3
c                     J4MIN=J4
c                     DMIN=DIST
c                  ENDIF
cC                 PRINT '(A,2I5,2G20.10)','J3,J4,DIST,DMIN=',J3,J4,DIST,DMIN
c               ENDDO
c            ENDDO
            CSOMGM1=-1.0D100
            CSOMGM2=-1.0D100
            DO J3=1,NRBSITES
                  PX=SITE(J3,1)
                  PY=SITE(J3,2)
                  PZ=SITE(J3,3)
                  DOT1=l1*PX + m1*PY + n1*PZ
                  C3A1DOT1=C3A1*DOT1
                  DOT2=l2*PX + m2*PY + n2*PZ
                  C3A2DOT2=C3A2*DOT2
                  CSOMG1=((C3A1DOT1*l1 - CA1*PX - n1*PY*S1 +
     @                 m1*PZ*S1)*(x1 - x2) +
     @              (C3A1DOT1*m1 - CA1*PY + n1*PX*S1 - l1*PZ*S1)*
     @               (y1 - y2) + (C3A1DOT1*n1 - CA1*PZ - m1*PX*S1 +
     @                 l1*PY*S1)*(z1 - z2))/r12
                  CSOMG2=((-(C3A2DOT2*l2) + CA2*PX + n2*PY*S2 - m2*PZ*S2)*
     @                (x1 - x2) + (-(C3A2DOT2*m2) + CA2*PY - n2*PX*S2 +
     @                  l2*PZ*S2)*(y1 - y2) +
     @               (-(C3A2DOT2*n2) + CA2*PZ + m2*PX*S2 - l2*PY*S2)*
     @                (z1 - z2))/r12
                  IF (CSOMG1.GT.CSOMGM1) THEN
                     J3MIN=J3
                     CSOMGM1=CSOMG1
                  ENDIF
                  IF (CSOMG2.GT.CSOMGM2) THEN
                     J4MIN=J3
                     CSOMGM2=CSOMG2
                  ENDIF
            ENDDO
c           PRINT*,'J3MIN,J4MIN=',J3MIN,J4MIN
           P1X=SITE(J3MIN,1)
           P1Y=SITE(J3MIN,2)
           P1Z=SITE(J3MIN,3)
           P2X=SITE(J4MIN,1)
           P2Y=SITE(J4MIN,2)
            P2Z=SITE(J4MIN,3)

            DOT1=l1*P1X + m1*P1Y + n1*P1Z 
            C3A1DOT1=C3A1*DOT1
            DOT2=l2*P2X + m2*P2Y + n2*P2Z 
            C3A2DOT2=C3A2*DOT2
c            ojo, esto lo he comentado yo
c            P1X=SITE(1,1)
c            P1Y=SITE(1,2)
c            P1Z=SITE(1,3)
c            P2X=SITE(1,1)
c            P2Y=SITE(1,2)
c            P2Z=SITE(1,3)
C           PRINT*,'coordinates of molecule pair:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
            CSOMG1=((C3A1DOT1*l1 - CA1*P1X - n1*P1Y*S1 + 
     @           m1*P1Z*S1)*(x1 - x2) +
     @        (C3A1DOT1*m1 - CA1*P1Y + n1*P1X*S1 - l1*P1Z*S1)*
     @         (y1 - y2) + (C3A1DOT1*n1 - CA1*P1Z - m1*P1X*S1 +
     @           l1*P1Y*S1)*(z1 - z2))/r12
            CSOMG1B= (-P1X*(x1-x2)-P1Y*(y1-y2)-P1z*(z1-z2))/r12
            OMG1B=ACOS(CSOMG1B)
            OMG1=ACOS(CSOMG1)
            IF(OMG1.LT.0.001D00) THEN
                FOMG1=1+OMG1**2/6
            ELSE
                FOMG1=OMG1/SIN(OMG1)
            ENDIF
            CSOMG2=((-C3A2DOT2*l2 + CA2*P2X + n2*P2Y*S2 - m2*P2Z*S2)*
     @         (x1 - x2) + (-(C3A2DOT2*m2) + CA2*P2Y - n2*P2X*S2 +
     @           l2*P2Z*S2)*(y1 - y2) +
     @        (-(C3A2DOT2*n2) + CA2*P2Z + m2*P2X*S2 - l2*P2Y*S2)*
     @         (z1 - z2))/r12
            CSOMG2B= (P2X*(x1-x2)+P2Y*(y1-y2)+P2z*(z1-z2))/r12
            OMG2B=ACOS(CSOMG2B)
            OMG2=ACOS(CSOMG2)
            IF(OMG2.LT.0.001D00) THEN
                FOMG2=1+OMG2**2/6
            ELSE
                FOMG2=OMG2/SIN(OMG2)
            ENDIF
            DUMMY= (-4*(r1212 - r126)*exp(-(OMG1**2+OMG2**2)/(2.*sigma22)))/(r1212*r126)
            DUMMYB= (-4*(r1212 - r126)*exp(-(OMG1B**2+OMG2B**2)/(2.*sigma22)))/(r1212*r126)
c	print*,'dummy=',dummy,'dummyb=',dummyb
            IF (GTEST.OR.SECT) THEN
               GX1=
     @ (4*(-12/r1214 + 6/r128)*(x1 - x2) - 
     @ (4*(1/r1212 - 1/r126)*( (-((-C3A2DOT2*l2 + CA2*P2X + 
     @              (n2*P2Y - m2*P2Z)*S2)/r12) + 
     @         ((x1 - x2)*CSOMG2)/RSQ)*FOMG2 +
     @    ( -((C3A1DOT1*l1 - CA1*P1X - (n1*P1Y - m1*P1Z)*S1)/r12) + 
     @         ((x1 - x2)*CSOMG1)/RSQ)*FOMG1))/sigma22)*
     @    exp((-OMG1**2-OMG2**2)/(2.*sigma22))

c	GX1=(4*(-12/r1214 + 6/r128)*(x1 - x2) - 
c     @ (4*(1/r1212 - 1/r126)*(( -P2X/r12+ ((x1 - x2)*CSOMG2)/RSQ)*OMG2B/SIN(OMG2B) +
c     @   (P1X/R12+((X1-X2)*CSOMG1B)/RSQ)*OMG1B/SIN(OMG1B)))/sigma22)*
c     @    exp((-OMG1B**2-OMG2B**2)/(2.*sigma22))


               V(3*(J1-1)+1)=V(3*(J1-1)+1)+GX1

               GY1=
     @  (4*(-12/r1214 + 6/r128)*(y1 - y2) - 
     @  (4*(1/r1212 - 1/r126)*( (-((-C3A2DOT2*m2 + CA2*P2Y + 
     @              (-(n2*P2X) + l2*P2Z)*S2)/r12) + 
     @         ((y1 - y2)*CSOMG2)/RSQ)*FOMG2+
     @    ( -((C3A1DOT1*m1 - CA1*P1Y - (-(n1*P1X) + l1*P1Z)*S1)/r12) + 
     @         ((y1 - y2)*CSOMG1)/RSQ)*FOMG1))/sigma22)*
     @   exp((-OMG1**2-OMG2**2)/(2.*sigma22))

               V(3*(J1-1)+2)=V(3*(J1-1)+2)+GY1

         GZ1=
     @   (4*(-12/r1214 + 6/r128)*(z1 - z2) - 
     @   (4*(1/r1212 - 1/r126)*( (-((-C3A2DOT2*n2 + CA2*P2Z + 
     @              (m2*P2X - l2*P2Y)*S2)/r12) + 
     @         ((z1-z2)*CSOMG2)/RSQ)*FOMG2+
     @    ( -((C3A1DOT1*n1 - CA1*P1Z - (m1*P1X - l1*P1Y)*S1)/r12) + 
     @         ((z1 - z2)*CSOMG1)/RSQ)*FOMG1))/sigma22)*
     @  exp((-OMG1**2-OMG2**2)/(2.*sigma22))
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+GZ1

               V(3*(NAT2+J1-1)+1)=V(3*(NAT2+J1-1)+1)+
     @   (4*(1/r1212 - 1/r126)*
     @   ((-C3A1DOT1 - (2*(1 - CA1)*DOT1*l1**2)/ALPHA1**4 - C3A1*l1*P1X + 
     @    (CA1*l1*(n1*P1Y - m1*P1Z))/ALPHA1**2 - l1*P1X*S1 + 
     @    (DOT1*l1**2*SA1)/ALPHA1**3 - 
     @    (l1*(n1*P1Y - m1*P1Z)*SA1)/(ALPHA1**3))*(-x1 + x2) + 
     @ (P1Z*S1 + l1*(-(P1Y*S1) + 
     @       (-(n1*P1X) + l1*P1Z)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) 
     @  + m1*(-(C3A1*P1X) + DOT1*l1*((-2*(1 - CA1))/ALPHA1**4 + 
     @          SA1/(ALPHA1**3))))*(-y1 + y2) + 
     @ (-(P1Y*S1) + l1*((CA1*(m1*P1X - l1*P1Y))/ALPHA1**2 - P1Z*S1) - 
     @    (l1*(m1*P1X - l1*P1Y)*SA1)/(ALPHA1**3) + 
     @    n1*((-2*(1 - CA1)*DOT1*l1)/ALPHA1**4 - C3A1*P1X + 
     @       (DOT1*l1*SA1)/(ALPHA1**3)))*(-z1 + z2))*FOMG1*
     @   exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/
     @   (r12*sigma22)

               V(3*(NAT2+J1-1)+2)=V(3*(NAT2+J1-1)+2)+
     @   (4*(1/r1212 - 1/r126)*
     @   ((-P1Z*S1 + m1*((CA1*(n1*P1Y - m1*P1Z))/ALPHA1**2 - P1X*S1) - 
     @    (m1*(n1*P1Y - m1*P1Z)*SA1)/ALPHA1**3 + 
     @    l1*((-2*(1 - CA1)*DOT1*m1)/ALPHA1**4 - C3A1*P1Y + 
     @       (DOT1*m1*SA1)/ALPHA1**3))*(-x1 + x2) + 
     @ (-C3A1DOT1 - (2*(1 - CA1)*DOT1*m1**2)/ALPHA1**4 - C3A1*m1*P1Y + 
     @    (CA1*m1*(-n1*P1X + l1*P1Z))/ALPHA1**2 - m1*P1Y*S1 + 
     @    (DOT1*m1**2*SA1)/ALPHA1**3 - 
     @    (m1*(-n1*P1X + l1*P1Z)*SA1)/ALPHA1**3)*(-y1 + y2) + 
     @ (P1X*S1 + m1*(-(P1Z*S1) + 
     @       (m1*P1X - l1*P1Y)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) + 
     @    n1*(-(C3A1*P1Y) + DOT1*m1*
     @        ((-2*(1 - CA1))/ALPHA1**4 + SA1/ALPHA1**3)))*(-z1 + z2)
     @  )*FOMG1*
     @   exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/
     @  (r12*sigma22)

               V(3*(NAT2+J1-1)+3)=V(3*(NAT2+J1-1)+3)+
     @ (4*(1/r1212 - 1/r126)* ((P1Y*S1 + n1*(-(P1X*S1) + 
     @       (n1*P1Y - m1*P1Z)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) + 
     @    l1*(-(C3A1*P1Z) + DOT1*n1*
     @        ((-2*(1 - CA1))/ALPHA1**4 + SA1/ALPHA1**3)))*
     @  (-x1 + x2) + (-(P1X*S1) + 
     @    n1*((CA1*(-(n1*P1X) + l1*P1Z))/ALPHA1**2 - P1Y*S1) - 
     @    (n1*(-(n1*P1X) + l1*P1Z)*SA1)/ALPHA1**3 + 
     @    m1*((-2*(1 - CA1)*DOT1*n1)/ALPHA1**4 - C3A1*P1Z + 
     @       (DOT1*n1*SA1)/ALPHA1**3))*(-y1 + y2) + 
     @ (-C3A1DOT1 - (2*(1 - CA1)*DOT1*n1**2)/ALPHA1**4 + 
     @    (CA1*n1*(m1*P1X - l1*P1Y))/ALPHA1**2 - C3A1*n1*P1Z - n1*P1Z*S1 + 
     @    (DOT1*n1**2*SA1)/ALPHA1**3 - 
     @    (n1*(m1*P1X - l1*P1Y)*SA1)/ALPHA1**3)*(-z1 + z2))*
     @  FOMG1*exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/
     @  (r12*sigma22)

               V(3*(J2-1)+1)=V(3*(J2-1)+1)-GX1

               V(3*(J2-1)+2)=V(3*(J2-1)+2)-GY1

               V(3*(J2-1)+3)=V(3*(J2-1)+3)-GZ1

               V(3*(NAT2+J2-1)+1)=V(3*(NAT2+J2-1)+1)+
     @  (4*(1/r1212 - 1/r126)*
     @  ((-C3A2DOT2 - (2*(1 - CA2)*DOT2*l2**2)/ALPHA2**4 - C3A2*l2*P2X + 
     @    (CA2*l2*(n2*P2Y - m2*P2Z))/ALPHA2**2 - l2*P2X*S2 + (DOT2*l2**2*SA2)/ALPHA2**3 - 
     @    (l2*(n2*P2Y - m2*P2Z)*SA2)/ALPHA2**3)*(x1 - x2) + (P2Z*S2 + l2*(-(P2Y*S2) + 
     @       (-(n2*P2X) + l2*P2Z)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) 
     @    + m2*(-(C3A2*P2X) + DOT2*l2*((-2*(1 - CA2))/ALPHA2**4 + 
     @          SA2/ALPHA2**3)))*(y1 - y2) + 
     @ (-(P2Y*S2) + l2*((CA2*(m2*P2X - l2*P2Y))/ALPHA2**2 - P2Z*S2) - 
     @    (l2*(m2*P2X - l2*P2Y)*SA2)/ALPHA2**3 + 
     @    n2*((-2*(1 - CA2)*DOT2*l2)/ALPHA2**4 - C3A2*P2X + 
     @       (DOT2*l2*SA2)/ALPHA2**3))*(z1 - z2))*FOMG2*
     @    exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/(r12*sigma22)

               V(3*(NAT2+J2-1)+2)=V(3*(NAT2+J2-1)+2)+
     @     (4*(1/r1212 - 1/r126)* ((-(P2Z*S2) + m2*((CA2*(n2*P2Y - m2*P2Z))/ALPHA2**2 - P2X*S2) - 
     @    (m2*(n2*P2Y - m2*P2Z)*SA2)/ALPHA2**3 + 
     @    l2*((-2*(1 - CA2)*DOT2*m2)/ALPHA2**4 - C3A2*P2Y + 
     @       (DOT2*m2*SA2)/ALPHA2**3))*(x1 - x2) + 
     @ (-C3A2DOT2 - (2*(1 - CA2)*DOT2*m2**2)/ALPHA2**4 - C3A2*m2*P2Y + 
     @    (CA2*m2*(-(n2*P2X) + l2*P2Z))/ALPHA2**2 - m2*P2Y*S2 + 
     @    (DOT2*m2**2*SA2)/ALPHA2**3 - 
     @    (m2*(-(n2*P2X) + l2*P2Z)*SA2)/ALPHA2**3)*(y1 - y2) + 
     @ (P2X*S2 + m2*(-(P2Z*S2) + 
     @       (m2*P2X - l2*P2Y)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) + 
     @    n2*(-(C3A2*P2Y) + DOT2*m2*
     @        ((-2*(1 - CA2))/ALPHA2**4 + SA2/ALPHA2**3)))*(z1 - z2))*
     @  FOMG2*exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/ (r12*sigma22)

               V(3*(NAT2+J2-1)+3)=V(3*(NAT2+J2-1)+3)+
     @  (4*(1/r1212 - 1/r126)* ((P2Y*S2 + n2*(-(P2X*S2) + 
     @       (n2*P2Y - m2*P2Z)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) + l2*(-(C3A2*P2Z) + DOT2*n2*
     @        ((-2*(1 - CA2))/ALPHA2**4 + SA2/ALPHA2**3)))*(x1 - x2) 
     @   + (-(P2X*S2) + n2*((CA2*(-(n2*P2X) + l2*P2Z))/ALPHA2**2 - P2Y*S2) - 
     @    (n2*(-(n2*P2X) + l2*P2Z)*SA2)/ALPHA2**3 + 
     @    m2*((-2*(1 - CA2)*DOT2*n2)/ALPHA2**4 - C3A2*P2Z + 
     @       (DOT2*n2*SA2)/ALPHA2**3))*(y1 - y2) + 
     @ (-C3A2DOT2 - (2*(1 - CA2)*DOT2*n2**2)/ALPHA2**4 + 
     @    (CA2*n2*(m2*P2X - l2*P2Y))/ALPHA2**2 - C3A2*n2*P2Z - n2*P2Z*S2 + 
     @    (DOT2*n2**2*SA2)/ALPHA2**3 - 
     @    (n2*(m2*P2X - l2*P2Y)*SA2)/ALPHA2**3)*(z1 - z2))*FOMG2*
     @   exp((-OMG1**2-OMG2**2)/(2.*sigma22)))/ (r12*sigma22)

            ENDIF
   
            ESTICKY=ESTICKY+DUMMY
            VT(J1)=VT(J1)+DUMMY
            VT(J2)=VT(J2)+DUMMY

         ENDDO
      ENDDO

c      WRITE(*,'(A,G20.10)') 'energy=',ESTICKY
c      PRINT*,'coords:'
c      WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
c      PRINT*,'gradient:'
c      WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

      RETURN
      END
