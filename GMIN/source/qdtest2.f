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
C  Energy and gradient for model system of dimension 3*NATOMS
C
      SUBROUTINE QDTEST2(X,V,EQD,GTEST)
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

!     EQD=0.25D0*X(1)**4-0.25D0*X(1)**2+(1.0D0/20.0D0)*X(1)
!     V(1)=X(1)**3-0.5D0*X(1)+(1.0D0/20.0D0)
!     EQD=0.25D0*X(1)**4-0.25D0*X(1)**2+(1.0D0/5.0D0)*X(1)
!     V(1)=X(1)**3-0.5D0*X(1)+(1.0D0/5.0D0)
!     EQD=0.25D0*X(1)**4
!     V(1)=X(1)**3
!     EQD=0.5D0*X(1)**2
!     V(1)=X(1)
!
! Try it with force constant = J1
!
      DO J1=2,3*NATOMS
         EQD=EQD+0.5D0*J1*X(J1)**2
         V(J1)=J1*X(J1)
      ENDDO

!     DO J1=1,3*NATOMS
!        EQD=EQD+0.25D0*J1*(X(J1)-1.0D0/J1)**4
!        V(J1)=J1*(X(J1)-1.0D0/J1)**3
!     ENDDO

      DO J1=1,1
!        EQD=EQD+0.25D0*X(J1)**4-0.2D0*X(J1)**2+0.03D0*X(J1) ! double well
!        V(J1)=X(J1)**3-0.4D0*X(J1)+0.03D0
         EQD=EQD+0.25D0*X(J1)**4-0.2D0*X(J1)**2+0.2D0*X(J1)  ! single well
         V(J1)=X(J1)**3-0.4D0*X(J1)+0.2D0
!        EQD=EQD+0.25D0*(X(J1)-0.80522D0)**4-0.2D0*(X(J1)-0.80522D0)**2+0.2D0*(X(J1)-0.80522D0)  ! single well shifted to origin
!        V(J1)=(X(J1)-0.80522D0)**3-0.4D0*(X(J1)-0.80522D0)+0.2D0
      ENDDO


      RETURN
      END
