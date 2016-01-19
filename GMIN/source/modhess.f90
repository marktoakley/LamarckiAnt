!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE MODHESS
      IMPLICIT NONE
      SAVE

      DOUBLE PRECISION, ALLOCATABLE :: HESS(:, :)  !  3*MXATMS,3*MXATMS

      ! hk286 - toggle between moving and stationary frame
      LOGICAL :: RBAANORMALMODET
      LOGICAL :: MASS_WEIGHTED = .FALSE.

CONTAINS

SUBROUTINE MASSWT()
      USE COMMONS, ONLY: NATOMS, ATMASS

      IMPLICIT NONE

      INTEGER :: J1, J2
      INTEGER :: X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION :: AMASS, BMASS, FMASS

! Only apply mass weighting if it hasn't already been applied.
! Currently disabled.
!      IF (.NOT. MASS_WEIGHTED) THEN
          DO J1 = 1, NATOMS
             X1 = 3 * J1 - 2
             Y1 = 3 * J1 - 1
             Z1 = 3 * J1
             AMASS = 1.D0/SQRT(ATMASS(J1))
             DO J2 = J1, NATOMS
                X2 = 3 * J2 - 2
                Y2 = 3 * J2 - 1
                Z2 = 3 * J2
                BMASS = 1.D0/SQRT(ATMASS(J2))
                FMASS = AMASS*BMASS
                HESS(X1,X2) = FMASS*HESS(X1,X2)
                HESS(X1,Y2) = FMASS*HESS(X1,Y2)
                HESS(X1,Z2) = FMASS*HESS(X1,Z2)
                HESS(Y1,X2) = FMASS*HESS(Y1,X2)
                HESS(Y1,Y2) = FMASS*HESS(Y1,Y2)
                HESS(Y1,Z2) = FMASS*HESS(Y1,Z2)
                HESS(Z1,X2) = FMASS*HESS(Z1,X2)
                HESS(Z1,Y2) = FMASS*HESS(Z1,Y2)
                HESS(Z1,Z2) = FMASS*HESS(Z1,Z2)
                IF (J1 .NE. J2) THEN
                   HESS(X2,X1) = HESS(X1,X2)
                   HESS(Y2,X1) = HESS(X1,Y2)
                   HESS(Z2,X1) = HESS(X1,Z2)
                   HESS(X2,Y1) = HESS(Y1,X2)
                   HESS(Y2,Y1) = HESS(Y1,Y2)
                   HESS(Z2,Y1) = HESS(Y1,Z2)
                   HESS(X2,Z1) = HESS(Z1,X2)
                   HESS(Y2,Z1) = HESS(Z1,Y2)
                   HESS(Z2,Z1) = HESS(Z1,Z2)
                ENDIF
             ENDDO
          ENDDO
!          MASS_WEIGHTED = .TRUE.
!      END IF    
END SUBROUTINE MASSWT

END MODULE MODHESS
