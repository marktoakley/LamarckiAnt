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
C  Energy and gradient for harmonic oscillator of dimension 3*NATOMS
C
      SUBROUTINE QDTEST(X,V,EQD,GTEST)
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), EQD
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT
  
      EVAP=.FALSE.
      EVAPREJECT=.FALSE.

      EQD=0.0D0
!
! Shift one minimum from the origin.
!
!     EQD=EQD+0.5D0*(X(1)-1.0D0)**2
!     V(1)=X(1)-1.0D0

      DO J1=1,3*NATOMS
         EQD=EQD+0.5D0*X(J1)**2
         V(J1)=X(J1)
!
! Try it with force constant = J1
!
!        EQD=EQD+0.5D0*J1*X(J1)**2
!        V(J1)=J1*X(J1)
!
! All degrees of freedom shifted from the origin
!
!        EQD=EQD+0.5D0*J1*(X(J1)-1.0D0/J1)**2
!        V(J1)=J1*(X(J1)-1.0D0/J1)
      ENDDO

      RETURN
      END
