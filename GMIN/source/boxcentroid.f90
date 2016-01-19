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

SUBROUTINE BOXCENTROID(X)
  !
  !ds656> Box centroid to specified rectangular region.
  !
  USE COMMONS, ONLY : NATOMS, BOXCENTROID_X, BOXCENTROID_DX, &
       BOXCENTROID_DISCRETE, MYUNIT
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS)
  LOGICAL :: SHIFT
  INTEGER :: I,J,K
  DOUBLE PRECISION :: CENTROID(3), OFFSET(3), DUMMY(3*NATOMS), E0, E1
  !
  ! First find centre of geometry
  CENTROID(:)=0.0D0; K=0
  DO I=1,NATOMS
     DO J=1,3
        K=K+1
        CENTROID(J) = CENTROID(J) + X(K)
     ENDDO
  ENDDO
  CENTROID(1:3) = CENTROID(1:3)/DBLE(NATOMS)
  !
  ! Now determine the offset by which we will be shifting
  SHIFT=.FALSE.
  DO I=1,3
     IF( ABS(CENTROID(I)-BOXCENTROID_X(I)) > BOXCENTROID_DX(I) ) THEN
        SHIFT=.TRUE.
        OFFSET(I) = CENTROID(I) - BOXCENTROID_X(I)
        IF(BOXCENTROID_DISCRETE(I)) THEN
           ! This rounding is for substrates with periodicity...
           ! in which case the DX ought to be equal to the period!
           OFFSET(I) = BOXCENTROID_DX(I)*DBLE(NINT(OFFSET(I)/BOXCENTROID_DX(I)))
        ENDIF
     ELSE
        OFFSET(I)=0.0D0
     ENDIF
  ENDDO
  !
  IF(SHIFT) THEN
     !
     !CALL POTENTIAL(X,DUMMY,E0,.false.,.false.)
     !WRITE(MYUNIT,'(A,3(1X,F15.8))') &
     !     'boxcog> Old CoG: ', (CENTROID(J), J=1,3)
     !WRITE(MYUNIT,'(A,3(1X,F15.8))') &
     !     'boxcog> Move by: ', (OFFSET(J), J=1,3)
     !
     ! Apply shift and compute new CoG
     CENTROID(:)=0.0D0; K=0
     DO I=1,NATOMS
        DO J=1,3
           K=K+1
           X(K) = X(K) - OFFSET(J)
           CENTROID(J) = CENTROID(J) + X(K)
        ENDDO
     ENDDO
     CENTROID(1:3) = CENTROID(1:3)/DBLE(NATOMS)
     WRITE(MYUNIT,'(A,3(1X,F15.8))') &
          'boxcentroid> Shifted to ', (CENTROID(J), J=1,3)
     !
     !CALL POTENTIAL(X,DUMMY,E1,.false.,.false.)
     !WRITE(MYUNIT,'(A,2(1X,F20.10))') &
     !     'boxcog> E0 and E1:', E0, E1
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE BOXCENTROID
