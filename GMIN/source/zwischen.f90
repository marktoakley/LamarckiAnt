!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
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
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
! the function DF1DIM has to accompany LINMIN
!

SUBROUTINE ZWISCHEN(XX, DF1DIM, POTEL1)
    USE PREC
    USE COMMONS, ONLY: NATOMS
    USE F1COM, ONLY: NCOM, XICOM, PCOM
    IMPLICIT NONE
! Arguments
    REAL(REAL64), INTENT(IN)    :: XX
    REAL(REAL64), INTENT(OUT)   :: DF1DIM
    REAL(REAL64), INTENT(OUT)   :: POTEL1
! Variables
    REAL(REAL64)                :: XT(3 * NATOMS), DF(3 * NATOMS)
    REAL(REAL64)                :: POTEL
    INTEGER(INT32)              :: J

    DO J = 1, NCOM
        XT(J) = PCOM(J) + XX * XICOM(J)
    END DO
    CALL POTENTIAL(XT, DF, POTEL, .TRUE., .FALSE.)
    POTEL1 = POTEL
    DF1DIM = 0.0D0
    DO J = 1, NCOM
        DF1DIM = DF1DIM + DF(J) * XICOM(J)
    END DO
!    WRITE(*, *) 'energy in zwischen=', POTEL
END SUBROUTINE ZWISCHEN
