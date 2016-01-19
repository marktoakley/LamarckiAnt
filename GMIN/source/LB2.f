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
C  Energy and gradient for Lynden-Bell^2 potential LB2
C
      SUBROUTINE LB2(X,V,ELB2,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4, I
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), ELB2, XMASS, YMASS, ZMASS

      ELB2=0.0D0
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO I=1,NATOMS
         XMASS=XMASS+X(3*(I-1)+1)
         YMASS=YMASS+X(3*(I-1)+2)
         ZMASS=ZMASS+X(3*(I-1)+3)
      ENDDO
      XMASS=XMASS/NATOMS
      YMASS=YMASS/NATOMS
      ZMASS=ZMASS/NATOMS
      VT(1:NATOMS)=0.0D0
      IF (GTEST) THEN
         V(1:3*NATOMS)=0.0D0
         DO J1=1,NATOMS
            J3=3*J1
            DIST=(X(J3-2)-XMASS)**2+(X(J3-1)-YMASS)**2+(X(J3)-ZMASS)**2
            ELB2=ELB2+DIST
            VT(J1)=VT(J1)+DIST    ! ignore the factor of 1/2 in the pair energies
            V(J3-2)=V(J3-2)+X(J3-2)-XMASS
            V(J3-1)=V(J3-1)+X(J3-1)-YMASS
            V(J3)=V(J3)+X(J3)-ZMASS
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=1.0D0/((X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2)
               ELB2=ELB2+DIST
               VT(J1)=VT(J1)+DIST
               VT(J2)=VT(J2)+DIST
               DIST=DIST**2
               V(J3-2)=V(J3-2)-(X(J3-2)-X(J4-2))*DIST
               V(J3-1)=V(J3-1)-(X(J3-1)-X(J4-1))*DIST
               V(J3)=V(J3)-(X(J3)-X(J4))*DIST
               V(J4-2)=V(J4-2)+(X(J3-2)-X(J4-2))*DIST
               V(J4-1)=V(J4-1)+(X(J3-1)-X(J4-1))*DIST
               V(J4)=V(J4)+(X(J3)-X(J4))*DIST
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=(X(J3-2)-XMASS)**2+(X(J3-1)-YMASS)**2+(X(J3)-ZMASS)**2
            ELB2=ELB2+DIST
            VT(J1)=VT(J1)+DIST
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               ELB2=ELB2+1.0D0/DIST
               VT(J1)=VT(J1)+1.0D0/DIST
               VT(J2)=VT(J2)+1.0D0/DIST
            ENDDO
         ENDDO
      ENDIF
      ELB2=ELB2/2.0D0

      RETURN
      END
