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

module NL_BIN_MOVEONE
   !neighbor lists for a binary system using one atom moves.  Maintain 3
   !neighbor lists for each atom.  One for each type of interaction.  When any
   !atom moves more than a given amount remake all the lists.
   !
   !Along with the neighbor lists, maintain lists of periodic correcitons for
   !each pair.  i.e. the boxl shift, so we can avoid NINT() calls 
   !

   !
   !NL_BIN_MOVEONE_LISTSAA( j, k ) is the jth type A neighbor of the kth atom
   !
   IMPLICIT NONE
   PRIVATE

   INTEGER, ALLOCATABLE :: NL_BIN_MOVEONE_LISTSAA(:,:)
   INTEGER, ALLOCATABLE :: NL_BIN_MOVEONE_LISTSBB(:,:)
   INTEGER, ALLOCATABLE :: NL_BIN_MOVEONE_LISTSAB(:,:)
   !
   !NLISTAA( k ) is the number of type A neighbros of atom k
   !
   INTEGER, ALLOCATABLE :: NLISTSAA(:)
   INTEGER, ALLOCATABLE :: NLISTSBB(:)
   INTEGER, ALLOCATABLE :: NLISTSAB(:)

   !
   ! NL_BIN_MOVEONE_XOFFSET_LISTSAA( 1:3, j, k ) is the (x,y,z) correction from
   !        periodic boundaries between atom k and the jth type A neighbor of atom k,
   !        stored in NL_BIN_MOVEONE_LISTSAA( j, k ) 
   !
   double precision, ALLOCATABLE :: NL_BIN_MOVEONE_XOFFSET_LISTSAA(:,:,:)
   double precision, ALLOCATABLE :: NL_BIN_MOVEONE_XOFFSET_LISTSBB(:,:,:)
   double precision, ALLOCATABLE :: NL_BIN_MOVEONE_XOFFSET_LISTSAB(:,:,:)

   integer NL_BIN_MOVEONE_nmoved
   integer, allocatable :: NL_BIN_MOVEONE_moved_atoms(:)

   double precision boxl, boxlvec(3), iboxlvec(3)
   integer natoms, ntypea
   DOUBLE PRECISION RSKIN, RCUT, REDO_DISPLACEMENT, REDO_DISPLACEMENT2, RLIST, RLIST2

   double precision, allocatable :: nl_bin_moveone_oldcoords(:) !the coords before the current move
   double precision, allocatable :: buildcoords(:) !the coords the last time the lists were rebuilt

   integer build_count, update_count
   integer nl_myunit

   PUBLIC :: NL_BIN_MOVEONE_SETUP, NL_BIN_MOVEONE_NMOVED, NL_BIN_MOVEONE_MOVED_ATOMS, &
      nl_bin_moveone_oldcoords, &
      NL_BIN_MOVEONE_XOFFSET_LISTSAA, &
      NL_BIN_MOVEONE_XOFFSET_LISTSBB, &
      NL_BIN_MOVEONE_XOFFSET_LISTSAB, nl_myunit, &
      NL_BIN_MOVEONE_UPDATE_LISTS, &
      NL_BIN_MOVEONE_LISTSAA, &
      NL_BIN_MOVEONE_LISTSBB, &
      NL_BIN_MOVEONE_LISTSAB, &
      NLISTSAA, NLISTSBB, NLISTSAB

   contains

   SUBROUTINE NL_BIN_MOVEONE_SETUP(COORDS, NATOMS_I, BOXL_I, RCUT_I, NTYPEA_i)
      implicit none
      integer, intent(in) :: natoms_i, ntypea_i
      double precision, intent(in) :: boxl_i, rcut_i
      double precision, intent(in) :: coords(3*natoms_i)

      natoms = natoms_i
      ntypea = ntypea_i
      build_count = 0
      update_count = 0

      !i'm allocating much too much memory here, but it makes indexing simpler
      allocate(NL_BIN_MOVEONE_LISTSAA( natoms, natoms) )
      allocate(NL_BIN_MOVEONE_LISTSBB( natoms, natoms) )
      allocate(NL_BIN_MOVEONE_LISTSAB( natoms, natoms) )
      allocate(NLISTSAA( natoms) )
      allocate(NLISTSBB( natoms) )
      allocate(NLISTSAB( natoms) )

      allocate(NL_BIN_MOVEONE_XOFFSET_LISTSAA( 3, natoms, natoms) )
      allocate(NL_BIN_MOVEONE_XOFFSET_LISTSBB( 3, natoms, natoms) )
      allocate(NL_BIN_MOVEONE_XOFFSET_LISTSAB( 3, natoms, natoms) )

      allocate(NL_BIN_MOVEONE_moved_atoms(natoms))

      allocate(nl_bin_moveone_oldcoords( 3*natoms) )
      allocate(buildcoords( 3*natoms) )
      nl_bin_moveone_oldcoords(:) = coords(:)

      BOXL = BOXL_i
      BOXLVEC(:) = BOXL
      IBOXLVEC(:) = 1.D0/BOXL

      RCUT = RCUT_I !THE CUTOFF FOR THE INTERACTION
      RSKIN = 0.8 !THIS COULD BE AN INPUT PARAMETER.
      RLIST = RSKIN + RCUT ! INCLUDE PAIRS IN THE INTERACTION LIST WITH SEPARATION LESS THAN RLIST
      RLIST2 = RLIST*RLIST
      !REDO_DISPLACEMENT = RSKIN / SQRT(3.D0) / 2.D0 !THIS SEEMS EXTREME
      REDO_DISPLACEMENT = RSKIN / 2.D0
      REDO_DISPLACEMENT2 = REDO_DISPLACEMENT**2

      call build_lists(coords)
   end subroutine NL_BIN_MOVEONE_SETUP


   subroutine build_lists(coords)
      !build all lists
      implicit none
      DOUBLE PRECISION, intent(in) :: coords(3*natoms)
      integer j1, j2
      DOUBLE PRECISION R2, xvec(3), dx(3)
      build_count = build_count + 1
      if (mod(build_count,1000) .eq. 0) then
         write(nl_myunit,*) "building lists", build_count, update_count
      endif
      !write(*,*) "NL_BIN_MOVEONE> building lists"
      !AA
      NLISTSAA(:) = 0
      DO J1=1,NTYPEA
         DO J2=1,j1-1
            XVEC(:) = COORDS(3*(J1-1)+1:3*(J1-1)+3) - COORDS(3*(J2-1)+1:3*(J2-1)+3)
            dx(:) = - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
            R2 = SUM((XVEC+dx)**2)
            IF (R2 .LE. RLIST2) THEN
               NLISTSAA(J1) = NLISTSAA(J1) + 1
               NLISTSAA(J2) = NLISTSAA(J2) + 1
               NL_BIN_MOVEONE_LISTSAA( NLISTSAA(J1) , J1 ) = J2
               NL_BIN_MOVEONE_LISTSAA( NLISTSAA(J2) , J2 ) = J1
               NL_BIN_MOVEONE_XOFFSET_LISTSAA(1:3, NLISTSAA(J1), J1) = DX
               NL_BIN_MOVEONE_XOFFSET_LISTSAA(1:3, NLISTSAA(J2), J2) = -DX
               !if ((j1.eq.879 .and. j2.eq.1042) .or.(j2.eq.879 .and. j1.eq.1042) ) then
                  !write(*,*) "setting up", j1, j2
                  !write(*,*) "      r2  ", sqrt(r2)
                  !write(*,*) COORDS(3*(J1-1)+1:3*(J1-1)+3)  
                  !write(*,*) COORDS(3*(J2-1)+1:3*(J2-1)+3)
                  !write(*,*) xvec
                  !write(*,*) dx
               !endif
            ENDIF
         ENDDO
      ENDDO

      !BB
      NLISTSBB(:) = 0
      DO J1=NTYPEA+1, NATOMS
         DO J2=J1+1, NATOMS
            XVEC(:) = COORDS(3*(J1-1)+1:3*(J1-1)+3) - COORDS(3*(J2-1)+1:3*(J2-1)+3)
            dx(:) = - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
            R2 = SUM((XVEC+dx)**2)
            IF (R2 .LE. RLIST2) THEN
               NLISTSBB(J1) = NLISTSBB(J1) + 1
               NLISTSBB(J2) = NLISTSBB(J2) + 1
               NL_BIN_MOVEONE_LISTSBB( NLISTSBB(J1) , J1 ) = J2
               NL_BIN_MOVEONE_LISTSBB( NLISTSBB(J2) , J2 ) = J1
               NL_BIN_MOVEONE_XOFFSET_LISTSBB(1:3, NLISTSBB(J1), J1) = DX
               NL_BIN_MOVEONE_XOFFSET_LISTSBB(1:3, NLISTSBB(J2), J2) = -DX
            ENDIF
         ENDDO
      ENDDO

      !AB
      NLISTSAB(:) = 0
      DO J1=1,NTYPEA
         DO J2=NTYPEA+1,NATOMS
            XVEC(:) = COORDS(3*(J1-1)+1:3*(J1-1)+3) - COORDS(3*(J2-1)+1:3*(J2-1)+3)
            dx(:) = - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
            R2 = SUM((XVEC+dx)**2)
            IF (R2 .LE. RLIST2) THEN
               NLISTSAB(J1) = NLISTSAB(J1) + 1
               NLISTSAB(J2) = NLISTSAB(J2) + 1
               NL_BIN_MOVEONE_LISTSAB( NLISTSAB(J1) , J1 ) = J2
               NL_BIN_MOVEONE_LISTSAB( NLISTSAB(J2) , J2 ) = J1
               NL_BIN_MOVEONE_XOFFSET_LISTSAB(1:3, NLISTSAB(J1), J1) = DX
               NL_BIN_MOVEONE_XOFFSET_LISTSAB(1:3, NLISTSAB(J2), J2) = -DX
            ENDIF
         ENDDO
      ENDDO
      buildcoords(:) = coords(:)
   END SUBROUTINE BUILD_LISTS

   SUBROUTINE NL_BIN_MOVEONE_UPDATE_LISTS(COORDS, NMOVED, MOVED_ATOMS)
      !update the lists if they need to be updated.
      !also update oldcoords
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      INTEGER, INTENT(IN) :: NMOVED, MOVED_ATOMS(NMOVED)
      INTEGER J1, J2
      DOUBLE PRECISION R2
      LOGICAL REBUILD
      update_count = update_count + 1
      rebuild = .false.
      if (nmoved .eq. natoms ) then
         rebuild = .true.
      else
         do j1=1,nmoved
            j2 = moved_atoms(j1)
            r2 = sum((coords(3*(j2-1)+1:3*(j2-1)+3) - buildcoords(3*(j2-1)+1:3*(j2-1)+3) )**2)
            if (r2 .gt. REDO_DISPLACEMENT2) then
               rebuild = .true.
               !write(*,*) "need to rebuild lists", r2, REDO_DISPLACEMENT2, j2
               !write(*,*) coords(3*(j2-1)+1:3*(j2-1)+3)
               !write(*,*) buildcoords(3*(j2-1)+1:3*(j2-1)+3)
               exit !break loop
            endif
            !write(*,*) "moved", sqrt(r2), REDO_DISPLACEMENT
         enddo
      endif
      if (rebuild) then
         call build_lists(coords)
      endif
      !update oldcoords
      if (nmoved .eq. natoms) then
         !if many atoms have moved, then just update all of them
         nl_bin_moveone_oldcoords(:) = coords(:)
      else
         do j1=1,nmoved
            j2 = moved_atoms(j1)
            nl_bin_moveone_oldcoords(3*(j2-1)+1:3*(j2-1)+3) = coords(3*(j2-1)+1:3*(j2-1)+3)
         enddo
      endif
   end SUBROUTINE NL_BIN_MOVEONE_UPDATE_LISTS
      
end module NL_BIN_MOVEONE
