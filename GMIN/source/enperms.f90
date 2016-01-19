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
!----------------------------------------------------------------------- 
!
! ds656> Use Lehmer's algorithm to step through distinct permutations 
! of N_A repetitions of A and N_B(=N-N_A) repetitions of B. 
!
! Details: "Generating Multiset Permutations in Constant Time"
! by Korsh and Lipschutz in J. Alg. 25, 321-335 (1997).
!
SUBROUTINE ENPERMS_TAKESTEP(NP,FINISH)
  !
  USE COMMONS, ONLY : NATOMS, LEHMER_LIST, LEHMER_ILASTB, LEHMER_LAST
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  LOGICAL, INTENT(INOUT) :: FINISH
  LOGICAL :: CURSOR
  INTEGER :: I, ICUR
  !
  ! Find the moving B item (i.e. cursor)
  CURSOR = .FALSE.
  ICUR = LEHMER_ILASTB(NP)
  DO WHILE(.NOT. CURSOR .AND. ICUR > 1)
     IF(LEHMER_LIST(ICUR-1,NP) == 'A') THEN
        CURSOR=.TRUE.
     ELSE
        ICUR = ICUR - 1
     ENDIF
  ENDDO
  !
  ! Terminate if no cursor found
  IF( .NOT. CURSOR) THEN
     FINISH = .TRUE.
     RETURN
  ENDIF
  !
  ! Update LEHMER_LIST while LEHMER_COORDS remains fixed
  !write(*,*) ICUR, FINISH, LEHMER_ILASTB(NP)
  IF(ICUR == LEHMER_ILASTB(NP)) THEN 
     CALL MOVE_CURSOR(NP,ICUR)
     LEHMER_ILASTB(NP) = LEHMER_ILASTB(NP) - 1
  ELSEIF( LEHMER_LIST(ICUR+1,NP)=='B' .AND. LEHMER_LAST(NP)=='A') THEN 
     CALL INVERT_TAIL(NP,ICUR)
     LEHMER_ILASTB(NP) = NATOMS
     CALL MOVE_CURSOR(NP,ICUR)
  ELSE
     CALL MOVE_CURSOR(NP,ICUR)
  ENDIF
  !
  LEHMER_LAST(NP) = LEHMER_LIST(NATOMS,NP) ! update the last entry
  !
  CALL UPDATE_COORDS(NP)
  !
  RETURN
  !
END SUBROUTINE ENPERMS_TAKESTEP
!
SUBROUTINE INVERT_TAIL(NP,ICUR)
  !
  USE COMMONS, ONLY : LEHMER_LIST, NATOMS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP,ICUR
  !
  INTEGER :: I,I1,I2,HTAIL
  CHARACTER(LEN=1) :: C
  !
  HTAIL = (NATOMS-ICUR)/2 ! NOTE: integer division rounds down
  !
  DO I = 1,HTAIL 
     !
     I1 = ICUR+I
     I2 = NATOMS-I+1
     C = LEHMER_LIST(I1,NP)
     LEHMER_LIST(I1,NP) = LEHMER_LIST(I2,NP)
     LEHMER_LIST(I2,NP) = C
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE INVERT_TAIL
!
SUBROUTINE MOVE_CURSOR(NP,ICUR)
  !
  USE COMMONS, ONLY : LEHMER_LIST
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP,ICUR
  !
  INTEGER :: ICM
  CHARACTER(LEN=1) :: C
  DOUBLE PRECISION :: Y
  !
  ICM = ICUR - 1
  C = LEHMER_LIST(ICUR,NP)
  LEHMER_LIST(ICUR,NP) = LEHMER_LIST(ICM,NP)
  LEHMER_LIST(ICM,NP) = C
  !
  RETURN
  !
END SUBROUTINE MOVE_CURSOR
!
SUBROUTINE UPDATE_COORDS(NP)
  !
  USE COMMONS, ONLY : LEHMER_COORDS,LEHMER_LIST,COORDS,NATOMS,NTYPEA
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  !
  INTEGER :: I,I3,J,J3,NA,NB
  !
  NA=0
  NB=0
  !
  DO I=1,NATOMS
     !
     I3 = 3*(I-1)
     !
     IF(LEHMER_LIST(I,NP) == 'A') THEN
        NA=NA+1
        J3=3*(NA-1)
     ELSEIF(LEHMER_LIST(I,NP) == 'B') THEN
        NB=NB+1
        J3=3*(NTYPEA+NB-1)
     ELSE
        WRITE(*,*) "enperms> Bad LIST entry:", LEHMER_LIST(I,NP)
        STOP
     ENDIF
     !
     DO J=1,3
        COORDS(J3+J,NP) = LEHMER_COORDS(I3+J,NP)
     ENDDO
     !
  ENDDO
  !
  !do I=1,NATOMS
  !   J3=3*(I-1)
  !   write(*,*) LEHMER_LIST(I,NP), (LEHMER_COORDS(J3+J,NP), J=1,3)
  !end do
  !
  ! Sanity checks:
  IF(NA /= NTYPEA) THEN
     WRITE(*,*) "enperms> Inconsistency: NA /= NTYPEA"
     STOP
  END IF
  !
  RETURN
  !
END SUBROUTINE UPDATE_COORDS
