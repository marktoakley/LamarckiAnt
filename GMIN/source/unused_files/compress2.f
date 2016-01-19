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
C  Add energy and gradient correction terms for compression
C
      SUBROUTINE COMPRESS(X,V,ENERGY,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*MXATMS), DIST, V(3*MXATMS), ENERGY, XMASS, YMASS, ZMASS

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+X(3*(J1-1)+1)
         YMASS=YMASS+X(3*(J1-1)+2)
         ZMASS=ZMASS+X(3*(J1-1)+3)
      ENDDO
      XMASS=XMASS/NATOMS
      YMASS=YMASS/NATOMS
      ZMASS=ZMASS/NATOMS

      IF (PERIODIC) RETURN

      DO J1=1,NATOMS
         J3=3*J1
C        DIST=(X(J3-2)-XMASS)**2+(X(J3-1)-YMASS)**2+(X(J3)-ZMASS)**2
         ENERGY=ENERGY+COMP*(ABS(X(J3-2)-XMASS)+ABS(X(J3-1)-YMASS)+ABS(X(J3)-ZMASS))
         IF (GTEST) THEN
	    IF (X(J3-2)-XMASS.NE.0.0D0) V(J3-2)=V(J3-2)+COMP*(X(J3-2)-XMASS)/ABS(X(J3-2)-XMASS)
	    IF (X(J3-1)-YMASS.NE.0.0D0) V(J3-1)=V(J3-1)+COMP*(X(J3-1)-YMASS)/ABS(X(J3-1)-YMASS)
	    IF (X(J3)  -ZMASS.NE.0.0D0) V(J3)=  V(J3)  +COMP*(X(J3)  -ZMASS)/ABS(X(J3)  -ZMASS)
         ENDIF
      ENDDO

      RETURN
      END
