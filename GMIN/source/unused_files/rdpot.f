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
C  Energy and gradient for model staircase surface.
C
      SUBROUTINE RDPOT(X,ENERGY)
      USE commons
      IMPLICIT NONE
      INTEGER J1, NDIV
      DOUBLE PRECISION X(3*MXATMS), ENERGY, DUM1, DUM2, DUM3, DUM4, RT
      PARAMETER (RT=17.320508D0)
C     PARAMETER (RT=DSQRT(300.0D0))
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      EVAP=.FALSE.
      ENERGY=0.0D0
      DUM1=0.0D0
      DUM2=0.0D0
      DUM3=0.0D0
      DUM4=0.0D0
      DO J1=1,3*NATOMS
         IF (X(J1)**2.GT.RADIUS) THEN
            ENERGY=ENERGY-(X(J1)**2-RADIUS)**2
         ENDIF
         DUM1=DUM1-(X(J1)+2.0D0)**2
         DUM2=DUM2-(X(J1)-1.0D0)**2
         DUM3=DUM3-(X(J1)-3.5D0)**2
         DUM4=DUM4+DSIN(RT*X(J1))
      ENDDO
      NDIV=3*NATOMS
      ENERGY=ENERGY+DEXP(5.0D0*DUM1/NDIV)+2.0D0*DEXP(DUM2/(2.0D0*NDIV))+2.2D0*DEXP(30.0D0*DUM3/NDIV)
     1       -(DUM4-NDIV)/(4.0D0*NDIV)
      ENERGY=ANINT(-5.0D0*ENERGY)

      RETURN
      END
