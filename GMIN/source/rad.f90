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
MODULE RAD_MOD

CONTAINS

SUBROUTINE RAD(X, V, ENERGY, GTEST)
   USE COMMONS, ONLY: NATOMS, BSPT, PERIODIC, PERCOLATET, RADIUS, DEBUG, MYUNIT, TWOD
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Arguments
   REAL(REAL64), INTENT(INOUT)   :: X(3*NATOMS)
   REAL(REAL64), INTENT(IN)      :: V(3*NATOMS)
   REAL(REAL64), INTENT(IN)      :: ENERGY
   LOGICAL, INTENT(IN)           :: GTEST
! Variables
   INTEGER(INT32)                :: I
   REAL(REAL64)                  :: DIST, RESCALE_FACTOR
   LOGICAL                       :: EVAP, EVAPREJECT
! TODO: Delete common block
   COMMON /EV/ EVAP, EVAPREJECT


   !WRITE(MYUNIT,'(A)') 'rad> A here' !ds656 commented this out!
   IF (BSPT .OR. PERIODIC .OR. PERCOLATET) RETURN ! container is accounted for in bspt by recounting previous configuration
                                                  ! jwrm2> we don't do radius checks if PERCOLATE is used
   EVAP = .FALSE.
   EVAPREJECT = .FALSE.
   DO I = 1, NATOMS
      DIST = SUM(X(3*I-2:3*I) ** 2)
      IF (DIST > RADIUS) THEN
!         WRITE(MYUNIT,'(A,I5,5G20.10)') 'J1,DIST,RADIUS in rad = ',J1,DIST,RADIUS,X(3*I-2),X(3*I-1),X(3*I)
         EVAP = .TRUE.
!         IF (EVAP.AND.(BSWL.OR.BSPT)) then
!            EVAPREJECT=.TRUE.
!            IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'EVAP: atom, radius=',J1,SQRT(DIST)
!            RETURN
!         END IF

         IF (DEBUG)  WRITE(MYUNIT,'(A,2G20.10,L10)') 'rad> EVAP: atom, radius, EVAP=', I, SQRT(DIST), EVAP
!         WRITE(*, *) 'EVAP: atom, radius=',J1,SQRT(DIST)
!         ENERGY=ENERGY+1.0D5*(DIST-RADIUS)**2
!         IF (GTEST.AND.(.NOT.(SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE))) THEN
!            DUMMYX=1.0D5*4.0D0*(DIST-RADIUS)*X(3*I-2)
!            DUMMYY=1.0D5*4.0D0*(DIST-RADIUS)*X(3*I-1)
!            DUMMYZ=1.0D5*4.0D0*(DIST-RADIUS)*X(3*I)
!            V(3*I-2)=V(3*I-2)+DUMMYX
!            V(3*I-1)=V(3*I-1)+DUMMYY
!            V(3*I)=V(3*I)+DUMMYZ
!         END IF

         RESCALE_FACTOR = (SQRT(RADIUS)-0.5D0) / SQRT(DIST)
! Move x, y and z coordinates
         X(3*I-2) = X(3*I-2) * RESCALE_FACTOR
         X(3*I-1) = X(3*I-1) * RESCALE_FACTOR
         IF (.NOT. TWOD) X(3*I) = X(3*I) * RESCALE_FACTOR
!         WRITE(MYUNIT,'(A,3G20.10)') 'rad> reset coords: ',X(3*I-2),X(3*I-1),X(3*I)

!  Put it back in at the opposite end of a diameter
!         X(3*I-2)=-X(3*I-2)*0.8D0
!         X(3*I-1)=-X(3*I-1)*0.8D0
!         X(3*I)=-X(3*I)*0.8D0
      END IF
   END DO

END SUBROUTINE RAD

!
!  For rigid-body angle-axis coordinates, just move the fixed site
!
SUBROUTINE RADR(X, V, ENERGY, GTEST)
   USE COMMONS, ONLY: NATOMS, PERIODIC, RADIUS, MYUNIT
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Arguments
   REAL(REAL64), INTENT(INOUT)   :: X(3*NATOMS)
   REAL(REAL64), INTENT(IN)      :: V(3*NATOMS)
   REAL(REAL64), INTENT(IN)      :: ENERGY
   LOGICAL, INTENT(IN)           :: GTEST
! Variables
   INTEGER(INT32)                :: I
   REAL(REAL64)                  :: DIST
   LOGICAL                       :: EVAP, EVAPREJECT
! TODO: Delete common block
   COMMON /EV/ EVAP, EVAPREJECT

   IF (PERIODIC) RETURN
   EVAP = .FALSE.
   DO I = 1, NATOMS/2
      DIST = SUM(X(3*I-2:3*I) ** 2)
      IF (DIST > RADIUS) THEN
         EVAP = .TRUE.
         WRITE(MYUNIT,'(A,I5,5G17.8)') 'EVAP: molecule, coords, dist, radius=', I, X(3*I-2:3*I), SQRT(DIST), SQRT(RADIUS)
!         RESCALE_FACTOR=SQRT(RADIUS*0.9D0/DIST)
!         X(3*I-2)=X(3*I-2)*RESCALE_FACTOR
!         X(3*I-1)=X(3*I-1)*RESCALE_FACTOR
!         X(3*I)=X(3*I)*RESCALE_FACTOR

!  Put it back in at the opposite end of a diameter
         X(3*I-2:3*I) = X(3*I-2:3*I) * (-0.8D0)
      END IF
   END DO

END SUBROUTINE RADR

SUBROUTINE RADCOM(X, MOVET)
   USE COMMONS, ONLY : NATOMS, DEBUG, MYUNIT, BSPT, PERIODIC, PERCOLATET, RADIUS
   USE PREC, ONLY: INT32, REAL64
   IMPLICIT NONE
! Arguments
   REAL(REAL64), INTENT(INOUT)   :: X(3*NATOMS)
   LOGICAL, INTENT(IN)           :: MOVET
! Variables
   INTEGER                       :: I
   REAL(REAL64)                  :: COM(3)
   REAL(REAL64)                  :: DIST
   LOGICAL                       :: EVAP, EVAPREJECT
!TODO: Remove common block
   COMMON /EV/ EVAP, EVAPREJECT

   IF (BSPT .OR. PERIODIC .OR. PERCOLATET) RETURN ! container is accounted for in bspt by recounting previous configuration
                                                  ! jwrm2> we don't do radius checks if PERCOLATE is used
   DO I = 1, 3
      COM(I) = SUM(X(I:3*NATOMS:3)) / NATOMS
   END DO

   EVAP = .FALSE.
   DO I = 1, NATOMS
      DIST = SUM((X(3*I-2:3*I)-COM(1:3)) ** 2)
      IF (DIST > RADIUS) THEN
         EVAP = .TRUE.
         IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10,L10)') 'radcom> EVAP: atom, radius, EVAP=', I, SQRT(DIST), EVAP
         IF (MOVET) THEN
            X(3*I-2:3*I) = X(3*I-2:3*I) * (SQRT(RADIUS) - 0.5D0) / SQRT(DIST)
         END IF
      END IF
   END DO

END SUBROUTINE RADCOM

END MODULE RAD_MOD
