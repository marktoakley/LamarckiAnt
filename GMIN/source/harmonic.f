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
C  Energy and gradient for harmonic field
C
      SUBROUTINE HARMONICFIELD(X,X0,V,E)
      USE COMMONS, ONLY : BOXLX, BOXLY, BOXLZ, HARMONICSTR
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION, INTENT(IN) :: X(3), X0(3)
      DOUBLE PRECISION, INTENT(OUT) :: V(3), E
      DOUBLE PRECISION DX,DY,DZ,DR

      !WRITE(*,'(1x,6F25.16)') X(1), X(2), X(3), X0(1), X0(2), X0(3)

      DX = X(1)-X0(1) 
      DY = X(2)-X0(2) 
      DZ = X(3)-X0(3) 
      !apply periodic boundary conditions
      DX = DX - ANINT(DX/BOXLX)*BOXLX
      DY = DY - ANINT(DY/BOXLY)*BOXLY
      DZ = DZ - ANINT(DZ/BOXLZ)*BOXLZ
      DR = SQRT(DX**2 + DY**2 + DZ**2)
      

      E = E + HARMONICSTR*DR**2
      V(1) = V(1) + HARMONICSTR*2.0D0*DX
      V(2) = V(2) + HARMONICSTR*2.0D0*DY
      V(3) = V(3) + HARMONICSTR*2.0D0*DZ

      !WRITE(*,'(1x,3F25.16)') DX, V(1), HARMONICSTR*2.0D0*DX

      END SUBROUTINE HARMONICFIELD
