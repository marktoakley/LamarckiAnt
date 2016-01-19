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

MODULE NEIGHBOR_LIST_MOD
   !a module to create and maintain a neighbor list
   !
   !atoms are listed as neighbors if they are closer then rlist = rcut + rskin
   !
   !if any atom moves more than REDO_DISPLACEMENT then the list is recreated
   !
   IMPLICIT NONE
   PRIVATE

   INTEGER, ALLOCATABLE :: NL_LIST(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: COORDSOLD(:)
   INTEGER NL_NLIST
   INTEGER NATOMS
   DOUBLE PRECISION RSKIN, RCUT, REDO_DISPLACEMENT, RLIST, RLIST2
   DOUBLE PRECISION BOXL(3), IBOXL(3)

   INTEGER NCHECK, NREDO
   !DOUBLE PRECISION TREDO 

   logical use_cell_lists ! use cell lists to construct the list

   PUBLIC :: NL_SETUP, NL_UPDATE, NL_LIST, NL_NLIST

   CONTAINS

   SUBROUTINE NL_SETUP(NATOMS_I, COORDS, RCUT_I, BOXLX, BOXLY, BOXLZ )
      use cell_lists_mod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), RCUT_I, BOXLX, BOXLY, BOXLZ
      INTEGER, INTENT(IN) :: NATOMS_I
      !LOGICAL, INTENT(IN) FREEZE, FROZEN(NATOMS)
      LOGICAL :: FIRST = .TRUE.

      IF (.NOT. FIRST) RETURN
      FIRST = .FALSE.
      !WRITE( *,*) "IN NL_SETUP"
      
      NATOMS = NATOMS_I
      !FREEZE = FREEZE_I

      BOXL(1) = BOXLX
      BOXL(2) = BOXLY
      BOXL(3) = BOXLZ
      IBOXL(:) = 1.D0/BOXL(:)

      RCUT = RCUT_I !THE CUTOFF FOR THE INTERACTION
      RSKIN = 0.5 !THIS COULD BE AN INPUT PARAMETER.
      RLIST = RSKIN + RCUT ! INCLUDE PAIRS IN THE INTERACTION LIST WITH SEPARATION LESS THAN RLIST
      RLIST2 = RLIST*RLIST
      !REDO_DISPLACEMENT = RSKIN / SQRT(3.D0) / 2.D0 !THIS SEEMS EXTREME
      REDO_DISPLACEMENT = RSKIN / 2.D0

      ALLOCATE( COORDSOLD(NATOMS*3) )

      !this list is way larger than it needs to be.  for very large systems you may have memory problems
      ALLOCATE( NL_LIST(2,NATOMS * (NATOMS-1)/2) ) 

      !TREDO = 0.D0

      use_cell_lists = .true.
      if (use_cell_lists) call cell_lists_setup(natoms, boxl(1), rlist)

      write(*,*) "neighbor_list_mod> potential cutoff ", RCUT
      write(*,*) "neighbor_list_mod> skin width       ", RSKIN
      if (use_cell_lists) write(*,*) "neighbor_list_mod> using cell lists to construct neighbor list"

      CALL CREATE_LIST(COORDS)
      !CALL CREATE_LIST_FROM_CELL_LISTS(COORDS)
   END SUBROUTINE NL_SETUP

   FUNCTION GET_SEP2(V1, V2) RESULT( R2 )
      DOUBLE PRECISION, INTENT(IN) :: V1(3), V2(3)
      DOUBLE PRECISION :: x
      DOUBLE PRECISION R2
      INTEGER K
      R2 = 0.d0
      DO K=1,3
         x = V1(k) - V2(k)
         x = x- BOXL(K) * NINT( x * IBOXL(K) )
         R2 = R2 + x*x
      ENDDO

   
   END FUNCTION GET_SEP2

   FUNCTION NEED_UPDATE( COORDS ) RESULT( REDO )
      !IF ANY OF THE ATOMS HAVE MOVED MORE THAN RMOVE
      !THEN RETURN TRUE
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      LOGICAL REDO
      DOUBLE PRECISION MAXDIST, D2
      INTEGER J1, K
      NCHECK = NCHECK + 1
      MAXDIST = 0.D0
      DO J1 = 1,NATOMS
         K = 3*(J1-1)
         D2 = GET_SEP2( COORDS(K+1 : K+3) , COORDSOLD(K+1 : K+3) )
         IF (D2 > MAXDIST) MAXDIST = D2
      ENDDO
      MAXDIST = SQRT(MAXDIST)
      REDO = (MAXDIST .GT. REDO_DISPLACEMENT)
   END FUNCTION NEED_UPDATE

   SUBROUTINE CREATE_LIST_from_cell_lists( COORDS )
      use cell_lists_mod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      DOUBLE PRECISION R2
      INTEGER J1, J2, K1, K2, n
      !DOUBLE PRECISION TIME0, TIME1
      NREDO = NREDO + 1

      !CALL MYCPU_TIME(TIME0)

      call make_neighbor_list(coords)

      NL_NLIST = 0

      DO n=1,cell_lists_neib_nlist
         j1 = cell_lists_neib_list(1, n)
         j2 = cell_lists_neib_list(2, n)
         K1 = 3*(j1-1)
         K2 = 3*(j2-1)
         R2 = GET_SEP2( COORDS(K1+1:K1+3) , COORDS(K2+1:K2+3) )
         IF (R2 .LE. RLIST2) THEN
            NL_NLIST = NL_NLIST + 1
            NL_LIST(1,NL_NLIST) = j1
            NL_LIST(2,NL_NLIST) = j2
         ENDIF
      ENDDO
      !CALL MYCPU_TIME(TIME1)
      !TREDO = TREDO + TIME1 - TIME0

      !write(*,*) "neighbor list", (natoms*natoms-natoms)/2
      !write(*,*) "              ", NL_NLIST
      !write(*,*) "time redo", TREDO, nredo, ncheck, nl_nlist
      !write(*,*) "         ", (natoms*(natoms-1))/2
      !write(*,*) "         ", cell_lists_neib_nlist
      !write(*,*) "         ", nl_nlist

   end SUBROUTINE CREATE_LIST_from_cell_lists

   SUBROUTINE CREATE_LIST( COORDS )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      DOUBLE PRECISION R2
      INTEGER J1, J2, K1, K2
      !DOUBLE PRECISION time0, time1
      if (use_cell_lists) then
         call create_list_from_cell_lists(coords)
         return
      endif
      NREDO = NREDO + 1
      !CALL MYCPU_TIME(TIME0)

      NL_NLIST = 0

      DO J1=1,NATOMS
         K1 = 3*(J1-1)
         DO J2=1,J1-1
            K2 = 3*(J2-1)
            R2 = GET_SEP2( COORDS(K1+1:K1+3) , COORDS(K2+1:K2+3) )
            IF (R2 .LE. RLIST2) THEN
               NL_NLIST = NL_NLIST + 1
               NL_LIST(1,NL_NLIST) = J1
               NL_LIST(2,NL_NLIST) = J2
            ENDIF
         ENDDO
      ENDDO
      !CALL MYCPU_TIME(TIME1)
      !TREDO = TREDO + TIME1 - TIME0

      !write(*,*) "neighbor list", (natoms*natoms-natoms)/2
      !write(*,*) "              ", NL_NLIST
      !write(*,*) "time redo", TREDO, nredo, ncheck, nl_nlist

   END SUBROUTINE CREATE_LIST

   SUBROUTINE NL_UPDATE(COORDS, CHANGED )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      LOGICAL, INTENT(OUT) :: CHANGED

      CHANGED = NEED_UPDATE(COORDS)
      !WRITE(*,*) "CHANGED", CHANGED
      IF (.NOT. CHANGED) RETURN

      !CREATE THE INTERACTION LIST
      CALL CREATE_LIST(COORDS)
      !CALL CREATE_LIST_FROM_CELL_LISTS(COORDS)
      !COORDSOLD(1:NATOMS) = COORDS(1:NATOMS)
      COORDSOLD(:) = COORDS(:)

   END SUBROUTINE NL_UPDATE

END MODULE NEIGHBOR_LIST_MOD
