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
C  Energy and gradient for the Stockmayer potential using two polar
C  coordinates for the dipole direction
C
      SUBROUTINE STOCK(X,V,ESTOCK,GTEST)
      USE commons

      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2, NAT2, J3, J4, REALNATOMS
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, ESTOCK, T1, P1, T2, P2,
     &                 R12, R122, R127, R129, R1212, R1214, MU2, R125, R126, N1DOTR12, N2DOTR12,
     &                 N1DOTN2, DUMMY, CT1, CT2, ST1, ST2, CP1, CP2, SP1, SP2

      REALNATOMS=NATOMS/2
      MU2=STOCKMU**2
      ESTOCK=0.0D0
      VT(1:REALNATOMS)=0.0D0
      V(1:6*REALNATOMS)=0.0D0
      
      DO J1=1,REALNATOMS
         J3=3*J1
         X1=X(J3-2)
         Y1=X(J3-1)
         Z1=X(J3)
         T1=X(3*REALNATOMS+J3-2)
         CT1=cos(T1)
         ST1=sin(T1)
         P1=X(3*REALNATOMS+J3-1)
         CP1=cos(P1)
         SP1=sin(P1)
         DO J2=J1+1,REALNATOMS
            J4=3*J2
            X2=X(J4-2)
            Y2=X(J4-1)
            Z2=X(J4)
            T2=X(3*REALNATOMS+J4-2)
            CT2=cos(T2)
            ST2=sin(T2)
            P2=X(3*REALNATOMS+J4-1)
            CP2=cos(P2)
            SP2=sin(P2)
            n1dotr12=st1*cp1*(x1-x2)+st1*sp1*(y1-y2)+ct1*(z1-z2)
            n2dotr12=st2*cp2*(x1-x2)+st2*sp2*(y1-y2)+ct2*(z1-z2)
            n1dotn2=st1*cp1*st2*cp2 + st1*sp1*st2*sp2 + ct1*ct2
            R122=(X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2
            R12=SQRT(R122)
            R126=R122**3
            R127=R126*R12
            R129=R127*R122
            R1212=R126**2
            R125=R126/R12
            R1214=R1212*R122
           
            DUMMY= (MU2*R127*(-3*n1dotr12*n2dotr12 + n1dotn2*R122) +
     &             (4 - 4*stocklambda*R126))/R1212
         
            VT(J1)=VT(J1)+DUMMY
            VT(J2)=VT(J2)+DUMMY
            ESTOCK=ESTOCK+DUMMY

! derivatives for positions

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(x1 - x2))
     &   + mu2*n2dotr12*R129*CP1*ST1 + mu2*n1dotr12*R129*CP2*ST2))/R1214
      V(J3-2)=V(J3-2)+DUMMY
      V(J4-2)=V(J4-2)-DUMMY

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(y1 - y2))
     &   + mu2*n2dotr12*R129*SP1*ST1 + mu2*n1dotr12*R129*SP2*ST2))/R1214
      V(J3-1)=V(J3-1)+DUMMY
      V(J4-1)=V(J4-1)-DUMMY

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(z1 - z2))
     &   + mu2*n2dotr12*R129*CT1 + mu2*n1dotr12*R129*CT2))/R1214
      V(J3)=V(J3)+DUMMY
      V(J4)=V(J4)-DUMMY

! derivatives for angular variables of atom J1

      V(3*REALNATOMS+J3-2) = V(3*REALNATOMS+J3-2) +
     &   (mu2*((3*n2dotr12*z1 - 3*n2dotr12*z2 - R122*CT2)*ST1 + 
     &   CP1*CT1*(3*n2dotr12*(-x1 + x2) + R122*CP2*ST2) + 
     &   CT1*SP1*(-3*n2dotr12*y1 + 3*n2dotr12*y2 + R122*SP2*ST2)))/R125

      V(3*REALNATOMS+J3-1) = V(3*REALNATOMS+J3-1) +
     &   (mu2*ST1*(SP1*(3*n2dotr12*(x1 - x2) - R122*CP2*ST2) + 
     &   CP1*(3*n2dotr12*(-y1 + y2) + R122*SP2*ST2)))/R125

! derivatives for angular variables of atom J2

      V(3*REALNATOMS+J4-2) = V(3*REALNATOMS+J4-2) +
     &   (mu2*(CP2*CT2*(3*n1dotr12*(-x1 + x2) + R122*CP1*ST1) + 
     &   CT2*SP2*(-3*n1dotr12*y1 + 3*n1dotr12*y2 + R122*SP1*ST1) + 
     &   (3*n1dotr12*z1 - 3*n1dotr12*z2 - R122*CT1)*ST2))/R125

      V(3*REALNATOMS+J4-1) = V(3*REALNATOMS+J4-1) +
     &   (mu2*(SP2*(3*n1dotr12*(x1 - x2) - R122*CP1*ST1) + 
     &   CP2*(3*n1dotr12*(-y1 + y2) + R122*SP1*ST1))*ST2)/R125

         ENDDO
      ENDDO

      RETURN
      END
