!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!=============================================================
!   All routines in this file were implemented by
!   Dmitri Schebarchov (ds656).
!=============================================================    
!
SUBROUTINE QALCS_SYM(NP,ITER,TIME,BRUN,QDONE,SCREENC)
  !
  ! Symmetrise current structure and keep the lowest-energy one.
  ! Many subroutines used here are contained in QALCS_SYM_MOD.
  !
  USE COMMONS, ONLY : NATOMS, COORDS, NQ, ECONV, MYUNIT, NSPECIES, QALCSV
  USE QALCS_SYM_MOD
  !
  IMPLICIT NONE
  !
  ! Parse passed variables (mainly for QUENCH)
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE 
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
  ! 
  ! Declare variables for internal use
  LOGICAL :: PROGRESS, COMMENSURATE
  INTEGER :: NQTOT, J
  DOUBLE PRECISION :: POTEL, X0(3*NATOMS), E0
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  IPROC = NP ! to avoid passing it back and forth for COORDS.
  !
  CALL ALLOCATE_ALL()
  CALL FIND_SYM_CORE(PROGRESS)
  !
  J=0
  DO WHILE(PROGRESS)
     !
     ! Store initial energy and coordinates
     X0(:) = COORDS(:,NP)
     E0 = POTEL
     !
     J=J+1
     !
     COMMENSURATE = .FALSE.
     CALL BUILD_LATTICE(COMMENSURATE)
     !
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '  ')
     !
     WRITE(MYUNIT, '(A,I3)', ADVANCE='NO') &
          'QALCS_sym> Iteration', J
     IF(POTEL < E0 - ECONV) THEN
        WRITE(MYUNIT,'(A)') ' led to lower energy!'
     ELSE
        WRITE(MYUNIT,'(A)') ' led to higher energy.'
     ENDIF
     !
     CALL FIND_SYM_CORE(PROGRESS)
     !
  ENDDO
  !
  !
  CALL DEALLOCATE_ALL()
  !
  RETURN
  !
END SUBROUTINE QALCS_SYM
