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
MODULE GAUSS_MOD

USE PREC, ONLY: REAL64
IMPLICIT NONE

REAL(REAL64)   :: PI = ATAN(1.0_REAL64) * 4.0_REAL64

CONTAINS

SUBROUTINE KEGEN()
   USE COMMONS, ONLY: GAUSSKK, GAUSSEE, GMODES, NATOMS, GKSMALL
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Variables
   INTEGER(INT32) :: I, J
   REAL(REAL64)   :: DPRAND

   DO J = 1, GMODES
      GAUSSEE(J) = 2.0_REAL64 * PI * DPRAND()  
      DO I = 1, 3*NATOMS
         GAUSSKK(I, J) = GASDEV()
      END DO
   END DO
   DO I = 1, 3*NATOMS
! GKSMALL(I) sets the length scale for the corresponding variable.
! The effective range is 2*Pi/GKSMALL(I)
      GKSMALL(I) = MINVAL(ABS(GAUSSKK(I, 1:GMODES)))
      WRITE(*, '(A, I8, G20.10)') 'I,smallest kk=', I, GKSMALL(I)
   END DO

END SUBROUTINE KEGEN

FUNCTION GASDEV()
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Arguments
   REAL(REAL64)               :: GASDEV
! Variables
   REAL(REAL64)               :: DPRAND
   INTEGER(INT32)             :: ISET = 0
   REAL(REAL64), SAVE         :: GSET
   REAL(REAL64)               :: FAC, RSQ, V1, V2

   IF (ISET == 0) THEN
      DO WHILE ( RSQ >= 1.0 .OR. RSQ == 0.0)
         V1 = 2.*DPRAND()-1.
         V2 = 2.*DPRAND()-1.
         RSQ = V1**2+V2**2
      END DO
      FAC = SQRT(-2.0_REAL64 * LOG(RSQ) / RSQ)
      GSET = V1 * FAC
      GASDEV = V2 * FAC
      ISET = 1
   ELSE
      GASDEV = GSET
      ISET = 0
   END IF

END FUNCTION GASDEV

SUBROUTINE GFIELD(X, V, FUNCVALUE, GTEST)
   USE COMMONS, ONLY: NATOMS, GMODES, GAUSSKK, GAUSSEE
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Arguments
   REAL(REAL64), INTENT(IN)   :: X(3*NATOMS)
   REAL(REAL64), INTENT(OUT)  :: V(3*NATOMS), FUNCVALUE
   LOGICAL, INTENT(IN)        :: GTEST
! Variables 
   INTEGER(INT32)             :: I
   REAL(REAL64)               :: NORM
   REAL(REAL64)               :: SPROD(GMODES)

   NORM = SQRT(2.0_REAL64/GMODES)

   DO I = 1, GMODES
      SPROD(I) = SUM(GAUSSKK(1:3*NATOMS, I) * X(1:3*NATOMS))
   END DO

   FUNCVALUE = NORM * SUM(COS(SPROD(1:GMODES) + GAUSSEE(1:GMODES)))

   DO I = 1, 3*NATOMS
      V(I) = -NORM * SUM(GAUSSKK(I, 1:GMODES)*SIN(SPROD(1:GMODES) + GAUSSEE(1:GMODES)))
   END DO

END SUBROUTINE GFIELD

END MODULE GAUSS_MOD
