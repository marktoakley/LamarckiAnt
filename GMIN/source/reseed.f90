!   GMIN: A program for finding global minima
!   Copyright (C) 1999 - 2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111 - 1307  USA
!
SUBROUTINE RESEED(NATOMS, P, RADIUS)
    USE PREC
    IMPLICIT NONE
! Arguments
    INTEGER(INT32), INTENT(IN)  :: NATOMS
    REAL(REAL64), INTENT(IN)    :: RADIUS
    REAL(REAL64), INTENT(OUT)   :: P(3 * NATOMS)
! Variables
    INTEGER(INT32)  :: J1
    REAL(REAL64)    :: DPRAND, SR_RADIUS_3
    EXTERNAL DPRAND
      
    SR_RADIUS_3 = DSQRT(RADIUS / 3.0D0)
    DO J1 = 1, NATOMS
! x
        P(3 * J1 - 2) = (DPRAND() - 0.5D0) * 2.0D0 * SR_RADIUS_3
! y
        P(3 * J1 - 1) = (DPRAND() - 0.5D0) * 2.0D0 * SR_RADIUS_3
! z
        P(3 * J1)     = (DPRAND() - 0.5D0) * 2.0D0 * SR_RADIUS_3
    ENDDO

END SUBROUTINE RESEED

