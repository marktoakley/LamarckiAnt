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

module cell_lists_binary_mod
   !a module for creating cell lists using linked-lists
   ! http://en.wikipedia.org/wiki/Cell_lists

   ! The module currently only works for cubic periodic systems, but could be
   ! easily modified
   !
   ! This version is specialized for binary systems.  Two copies of each list
   ! are kept.  
   !
   ! This version is further specialize for one atom moves.  I.e. many
   ! subroutines deal with updating only one atom at a time
   implicit none
   private

   double precision rcut  !potential cutoff
   double precision rcell !lenth of a cell
   double precision boxl  !box size (assume cubic)
   integer ncellx !number of cells in the x direction (assume cubic box)
   integer ncellx2
   integer ncells !total number of cells

   integer, allocatable :: cell_neib(:,:)  ! a list containing neighboring cells, no repeats
   integer :: ncell_neib                   ! the length of cell_neib
   integer, allocatable :: cell_neib_all(:,:)  ! list all neighbors for each cell
   integer, allocatable :: ncell_neib_all(:)  ! the number of neighbors of each cell

   !here we gather information specific to constructing a cell lists
   integer natoms , ntypeA, ntypeB
   integer, allocatable :: hocA(:) !head of cell:  hoc(icell) is the first atom in cell icell
   integer, allocatable :: llA(:)  !linked list:    ll(i) is the next atom after i in the cell
                                  !             if ll(i)==0 then there are no more atoms in the cell
   integer, allocatable :: hocB(:)
   integer, allocatable :: llB(:)

   integer, allocatable :: CL_listAA(:,:)  ! a list containing neighboring atoms
   integer, allocatable :: CL_listBB(:,:)  ! a list containing neighboring atoms
   integer, allocatable :: CL_listAB(:,:)  ! a list containing neighboring atoms
   integer CL_nlistAA !length of the above lists
   integer CL_nlistBB
   integer CL_nlistAB

   integer, allocatable :: typeAlist(:)  !A list of type A atoms
   integer, allocatable :: typeBlist(:)  !A list of type B atoms


   !logical debug
   DOUBLE PRECISION, ALLOCATABLE :: cl_OLDCOORDS(:) !the old coordinates.  
      !Used to calculate old atom-atom interaction energies in order to subract
      !out
      ! If atom i has moved
      ! new_energy = old_energy + sum_j ( Enew_{ij} - Eold_{ij} )

   double precision :: iboxl, ircell !inverse box length and cell length

   integer, allocatable :: in_cell(:) !save which cell each atom is in

   integer, allocatable :: CL_MOVED_ATOMS(:)  !the list of which atoms have moved
   INTEGER CL_NMOVED   !the number of moved atoms

   public :: cell_lists_binary_setup, CL_NLISTAA, CL_LISTAA, &
      CL_NLISTBB, CL_LISTBB, CL_NLISTAB, CL_LISTAB, &
      CL_NMOVED, CL_MOVED_ATOMS, &
      cl_rebuild_cell_lists, cl_oldcoords, cl_get_neighbors_several, &
      cl_update_cell_lists

   contains

   subroutine cell_lists_binary_setup(coords, natoms_i, boxl_i, rcut_i, ntypeA_i)
      implicit none
      integer, intent(in) :: natoms_i, ntypeA_i
      double precision, intent(in) :: boxl_i, rcut_i, coords(3*natoms)
      integer i, j
      logical :: first = .true.
      if (.not. first) return
      first = .false.
      write(*,*) "rcell before assignment", rcell


      rcut = rcut_i
      natoms = natoms_i
      boxl = boxl_i
      iboxl = 1.d0/boxl

      !determine the number of cells in the linear direction.  This is arbitrary,
      !and should be chosen based on what is fastest.  rcell should not be larger
      !than rcut, but could be smaller.
      ncellx = int( boxl / rcut ) * 2
      if (ncellx .lt. 8) ncellx = 8
      ncellx2 = ncellx**2

      rcell = boxl / ncellx !rcell must be greater than rcut
      ncells = ncellx**3
      ircell = 1.d0/rcell

      write(*,*) "cell_lists_binary_mod> rcut   ", rcut
      write(*,*) "cell_lists_binary_mod> rcell  ", rcell
      write(*,*) "cell_lists_binary_mod> ncellx ", ncellx
      write(*,*) "cell_lists_binary_mod> ncells ", ncells
      write(*,*) "cell_lists_binary_mod> atoms per cell ", float(natoms) / ncells

      allocate(hocA( ncells ) )
      allocate(llA( natoms ) )
      allocate(hocB( ncells ) )
      allocate(llB( natoms ) )
      allocate(cell_neib( 2, ncells*(ncells+1)/2 ) ) !this allocates way too much memory
      allocate(cell_neib_all( ncells, ncells ) ) !this allocates way too much memory
      allocate(ncell_neib_all( ncells ) )
      ntypeA = ntypeA_i
      ntypeB = natoms - ntypeA
      allocate(CL_listAA( 2,ntypeA*(ntypeA-1)/2) )
      allocate(CL_listBB( 2,ntypeB*(ntypeB-1)/2) )
      allocate(CL_listAB( 2,ntypeA*ntypeB) )
      allocate(typeAlist(ntypeA))
      allocate(typeBlist(ntypeB))
      do i=1,ntypeA
         typeAlist(i) = i
      enddo
      do i=1,ntypeB
         j=ntypeA+i
         typeBlist(i) = j
      enddo
      allocate(in_cell(natoms))


      CL_NMOVED = 0
      allocate(CL_MOVED_ATOMS(natoms))

      call setup_cell_neibs()
      call cl_rebuild_cell_lists(coords)

      !debug = .true.
      allocate(cl_oldcoords(3*natoms))
      cl_oldcoords(:) = coords(:)


   end subroutine cell_lists_binary_setup

   subroutine xyz2icell( xyz_i, icell )
      !determine the cell index from the xyz coordinates
      implicit none
      double precision, intent(in) :: xyz_i(3)
      integer, intent(out) :: icell
      double precision :: xyz(3)
      integer ijk(3), i
      xyz(:) = xyz_i(:)
      !put x, y, z in box [0,boxl)
      do i=1,3
         xyz(i) = xyz(i) - boxl * nint( xyz(i) * iboxl )
         if (xyz(i) .lt. 0) xyz(i) = xyz(i) + boxl
         ijk(i) = floor(xyz(i) * ircell)
      enddo
      icell = ijk(1) + ijk(2)*ncellx + ijk(3) * ncellx2 +1 !+1 for fortran indexing
   end subroutine xyz2icell

   subroutine icell2xyz( icell_i, xyz )
      !return the xyz coordinates of the corner of the cell
      implicit none
      double precision, intent(out) :: xyz(3)
      integer, intent(in) :: icell_i
      integer icell, ijk(3)
      icell = icell_i - 1

      ijk(3) = int( (icell) / ncellx2 ) 
      ijk(2) = int( (icell - ijk(3)*ncellx2 ) / ncellx ) 
      ijk(1) = int( (icell - ijk(3)*ncellx2 - ijk(2)*ncellx ) ) 

      xyz(3) = rcell * ijk(3)
      xyz(2) = rcell * ijk(2)
      xyz(1) = rcell * ijk(1)
      !write(*,*) "icell2xyz", icell_i, xyz, ijk, ncellx2, icell/ncellx2
   end subroutine icell2xyz

   function are_cell_neibs(icell, jcell) result( b )
      !determine if icell and jcell are close enough together
      !find the minimum distance between the two cells
      implicit none
      double precision :: xyz1(3), xyz2(3), r2, x, dxmin
      logical b
      integer i,k, icell, jcell
      call icell2xyz( icell, xyz1 )
      call icell2xyz( jcell, xyz2 )
      xyz2 = xyz2(:) - xyz1(:)
      !write(*,*) "xyz2", xyz2(:)
      do i=1,3 !loop over x,y,z
         dxmin = 1000.d0
         do k=-1,1 !determine minimum distance in this direction
            x = xyz2(i) + dble(k) * rcell
            x = x - boxl*nint(x/boxl)
            if (abs(x) .lt. abs(dxmin)) then
               dxmin = x
               !write(*,*) "dxmin k", dxmin, k
            endif
         enddo
         xyz2(i) = dxmin
         !write(*,*) "xyz2(i), i", xyz2(i), i
      enddo
      !write(*,*) "xyz2", xyz2(:)
      r2 = sum( xyz2(:)**2 )
      !write(*,*) r2, icell, jcell
      b = (r2 .lt. rcut**2)
      !if ( b) write(*,*) sqrt(r2), rcut, icell, jcell
   end function are_cell_neibs


   subroutine setup_cell_neibs()
      !determine which cells are neighbors.
      !a cell is it's own neighbor
      implicit none
      integer n, icell, jcell

      ncell_neib_all(:) = 1
      do n=1,ncells !every cell is it's own neighbor
         cell_neib_all(1, n) = n
      enddo


      n = 0
      do icell = 1,ncells
         n = n + 1
         !write(*,*) "cell neibs", icell, icell, n
         cell_neib(1,n) = icell
         cell_neib(2,n) = icell
      enddo
      do icell = 1,ncells
         do jcell = icell+1,ncells
            if ( are_cell_neibs(icell, jcell) ) then
               n = n + 1
               !write(*,*) "cell neibs", icell, jcell, n
               cell_neib(1,n) = icell
               cell_neib(2,n) = jcell
               !update cell_neib_all
               ncell_neib_all(icell) = ncell_neib_all(icell) + 1
               ncell_neib_all(jcell) = ncell_neib_all(jcell) + 1
               cell_neib_all( ncell_neib_all(icell), icell ) = jcell
               cell_neib_all( ncell_neib_all(jcell), jcell ) = icell
            endif
         enddo
      enddo
      ncell_neib = n
      write(*,*) "number of cell-cell neighbors", ncell_neib, ncells*14
   end subroutine setup_cell_neibs

   subroutine generate_cell_lists(coords, hoc, ll, natomlist, atomlist)
      !determine which cells all atoms in atomlist are in.
      !Reset and populate hoc and ll.  
      !hoc and ll are passed so that this function can be used for hocA, and
      !hocB, etc.
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: natomlist, atomlist(natomlist)
      integer, intent(out) :: hoc(ncells), ll(natoms)
      integer i, j, icell
      hoc(:) = 0
      ll(:) = 0
      !write(*,*) "starting generate_cell_lists", natomlist, 4, atomlist(4)
      do i=1,natomlist
         j = atomlist(i)
         !if (j .eq. 516) write(*,*) "generate_cell_lists atomlist", i, j
         !write(*,*) i, j
         call xyz2icell( coords(3*(j-1) + 1 : 3*(j-1) + 3), icell )
         in_cell(j) = icell
         call add_to_cell(coords, hoc, ll, j, icell)
         !ll(j) = hoc(icell)
         !hoc(icell) = j
      enddo
      !write(*,*) "ending generate_cell_lists atomlist", 4, (atomlist(4))
   end subroutine generate_cell_lists

   subroutine CL_rebuild_cell_lists(coords)
      !rebuild all cell lists for all atoms
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      call generate_cell_lists(coords, hocA, llA, ntypeA, typeAlist)
      call generate_cell_lists(coords, hocB, llB, ntypeB, typeBlist)
   end subroutine CL_rebuild_cell_lists

   subroutine remove_from_cell( coords, hoc, ll, imove, icell)
      !update hoc and ll to reflect that atom imove is no longer in cell icell
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: imove, icell
      integer, intent(inout) :: hoc(ncells), ll(natoms)
      integer :: j, prev
      !write(*,*) "removing", imove, "from cell", icell, hoc(icell), hoc(1)
      j = hoc(icell)
      if (j .eq. imove) then
         !write(*,*) imove, "was head of cell", icell, hoc(icell)
         hoc(icell) = ll(j)
         !write(*,*) ll(j), "is now head of cell", icell, hoc(icell)
         return
      endif
      if ( j .ne. 0 ) then
         prev = j
         j = ll(j)
         do while( j .ne. 0 )
            !write(*,*) "remove_from_cell>", j, imove
            if (j .eq. imove) then
               !make ll skip over imove
               !currently ll(prev) == imove
               !write(*,*) "ll(prev) = ll(imove)", prev, imove, ll(prev), ll(imove)
               !write(*,*) prev, "->", imove, ll(prev) == imove
               !write(*,*) imove, "->", ll(imove)
               ll(prev) = ll(imove)
               !write(*,*) prev, "->", ll(prev)
               return
            endif
            prev = j
            j = ll(j)
         enddo
      endif
      !should not get here
      write(*,*) "remove_from_cell> warning: atom", imove, "not in cell", icell, in_cell(imove)

      j = hoc(icell)
      write(*,*) "        hoc( ", icell,") = ", j
      do while (j .ne. 0)
         write(*,*) "     ", j, "->", ll(j)
         j = ll(j)
      enddo
   end subroutine remove_from_cell

   subroutine add_to_cell( coords, hoc, ll, imove, icell)
      !update hoc and ll to reflect that atom imove is now in cell icell
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: imove, icell
      integer, intent(inout) :: hoc(ncells), ll(natoms)
      integer :: j, prev
      !add it to the front of the cell
      ll(imove) = hoc(icell)
      hoc(icell) = imove
      return
   end subroutine add_to_cell

   subroutine update_cell_lists(coords, hoc, ll, imove)
      !update hoc and ll for the current position of atom imove.
      !If it's in a new cell, remove it from it's old cell and add it to the new
      !one
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: imove
      integer, intent(inout) :: hoc(ncells), ll(natoms)
      integer i, j, icell, imovecell

      !is imove in the same cell?
      call xyz2icell( coords(3*(imove-1) + 1 : 3*(imove-1) + 3), imovecell )
      if ( imovecell .eq. in_cell(imove) ) then
         !no need to update anything
         !write(*,*) "update_cell_lists> nothing to do", imove
         return
      endif
  
      !write(*,*) "new cell", imove, imovecell, in_cell(imove)
      call remove_from_cell( coords, hoc, ll, imove, in_cell(imove) )
      call add_to_cell( coords, hoc, ll, imove, imovecell )
      in_cell(imove) = imovecell
   end subroutine update_cell_lists

   subroutine cl_update_cell_lists(coords, nmove, movelist)
      !update all cell lists to reflect the new positions of atoms in movelist
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: nmove, movelist(nmove)
      integer n, imove
      do n=1,nmove
         imove = movelist(n)
         if (imove .le. ntypeA) then
            call update_cell_lists(coords, hocA, llA, imove)
         else
            call update_cell_lists(coords, hocB, llB, imove)
         endif
      enddo
   end subroutine cl_update_cell_lists

!   subroutine add_to_atom_neighbor_list(coords, i, hoc, ll, nlist, list)
!      !Find all neighbors of i in hoc and ll and add them to the list nlist
!      implicit none
!      double precision, intent(in) :: coords(3*natoms)
!      integer, intent(in) :: i
!      integer, intent(in) :: hoc(ncells), ll(natoms)
!      integer, intent(inout) :: nlist, list(2,natoms*natoms) !this length is wrong, but I think it doesn't matter
!      integer :: icell, jcell, j, n
!
!      call xyz2icell( coords(3*(i-1) + 1 : 3*(i-1) + 3), icell )
!      do n=1, ncell_neib_all( icell )
!         jcell = cell_neib_all( n, icell )
!         !loop through atoms in cell jcell
!         j = hoc( jcell )
!         do while (j .ne. 0)
!            if (i .ne. j) then
!               !add pair j to list of neighbors of i
!               nlist = nlist + 1
!               list(1, nlist ) = i
!               list(2, nlist ) = j
!            endif
!            !go to next atom in cell
!            j = ll(j)
!         enddo
!      enddo
!   end subroutine add_to_atom_neighbor_list

   function is_in_array(i, n, list) result ( b)
      !return true if i is in the list, else return false
      !assume the list is sorted from low to high (to make it run faster)
      implicit none
      integer, intent(in) :: i, n, list(n)
      integer j
      logical b
      b = .false.
      do j=1,n
         if ( i .le. list(j) ) then
            if ( i .eq. list(j) ) then
               b = .true.
               !write(*,*) "is_in_array>", i, j
            endif
            return
         endif
      enddo
   end function is_in_array

   subroutine make_atom_neighbor_list2(coords, i, hoc, ll, nlist, list, navoid, avoidlist)
      !add the neighbors of i to the list of interactions
      !to avoid duplicates, don't add if the neighbor is in avoidlist
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: i
      integer, intent(in) :: hoc(ncells), ll(natoms)
      integer, intent(inout) :: nlist, list(2,natoms*natoms) !this is incorrect, but I hope it wont matter
      integer, intent(in) :: navoid, avoidlist(navoid)
      integer :: icell, jcell, j, n
      !logical is_in_array

      !nlist = 0

      call xyz2icell( coords(3*(i-1) + 1 : 3*(i-1) + 3), icell )
      do n=1, ncell_neib_all( icell )
         jcell = cell_neib_all( n, icell )
         !loop through atoms in cell jcell
         j = hoc( jcell )
         !write(*,*) "hoc(",jcell,") = ", j
         do while (j .ne. 0)
            if (i .ne. j .and. (.not. is_in_array(j, navoid, avoidlist))) then
               !add pair j to list of neighbors of i
               nlist = nlist + 1
               list(1, nlist ) = i
               list(2, nlist ) = j
            endif
            j = ll(j)
         enddo
      enddo
   end subroutine make_atom_neighbor_list2

!   subroutine make_neighbor_list_all(coords, hoc, ll, nlist, list)
!      !add the neighbors of i to the list of interactions
!      !to avoid duplicates, don't add if the neighbor is in avoidlist
!      implicit none
!      double precision, intent(in) :: coords(3*natoms)
!      integer, intent(in) :: hoc(ncells), ll(natoms)
!      integer, intent(inout) :: nlist, list(2,natoms*natoms) !this is incorrect, but I hope it wont matter
!      integer :: icell, jcell, j, n, i
!      nlist = 0
!      do n = 1,ncell_neib
!         icell = cell_neib(1,n)
!         jcell = cell_neib(2,n)
!         if (icell .ne. jcell) then
!            !loop over all atoms in each cell
!            i = hoc( icell )
!            do while (i .ne. 0)
!               j = hoc( jcell )
!               do while (j .ne. 0)
!                  nlist = nlist + 1
!                  list(1,nlist) = i
!                  list(2,nlist) = j
!                  j = ll(j)
!               enddo
!               i = ll(i)
!            enddo
!         else
!            !loop over all atoms in each cell, avoiding duplicates
!            i = hoc( icell )
!            do while (i .ne. 0)
!               j = hoc( jcell )
!               do while (j .ne. i)
!                  nlist = nlist + 1
!                  list(1,nlist) = i
!                  list(2,nlist) = j
!                  j = ll(j)
!               enddo
!               i = ll(i)
!            enddo
!         endif
!      enddo
!   end subroutine make_neighbor_list_all


!   subroutine check_coords(coords, i)
!      implicit none
!      double precision, intent(in) :: coords(3*natoms)
!      integer, intent(in) :: i
!      integer j, k
!      double precision r
!      write(*,*) "check_coords", i
!      if (cl_oldcoords(1) .eq. 0.d0) then
!         cl_oldcoords(:) = coords(:)
!      else
!         do j=1,natoms
!            k=3*(j-1)
!            r = sum( (coords(k+1:k+3) - cl_oldcoords(k+1:k+3))**2 )
!            if (j .ne. i) then
!               if ( r .ne. 0.d0 ) then
!                  write(*,*) "ERROR: coords changed", i, j, sqrt(r)
!                  write(*,'(A,3G16.9)') "               ", cl_oldcoords(k+1:k+3)
!                  write(*,'(A,3G16.9)') "               ", coords(k+1:k+3)
!                  write(*,'(A,3G16.9)') "               ", coords(k+1:k+3)- cl_oldcoords(k+1:k+3)
!                  k=3*(i-1)
!                  write(*,'(A,3G25.12)') "               ", cl_oldcoords(k+1:k+3)
!                  write(*,'(A,3G16.9)') "               ", coords(k+1:k+3)
!               endif
!            endif
!         enddo
!      endif
!   end subroutine check_coords

!   subroutine which_atom_moved(coords )
!      implicit none
!      double precision, intent(in) :: coords(3*natoms)
!      integer j, k, i
!      double precision r
!      !write(*,*) "coords(891) begin", coords((891-1)*3+1), cl_oldcoords((891-1)*3+1)
!      !write(*,*) "which_atom_moved> ", coords(1), cl_oldcoords(1)
!      if (cl_oldcoords(1) .eq. 0.d0) then
!         write(*,*) "which_atom_moved> warning: cl_oldcoords was not initialized"
!      endif
!      CL_NMOVED = 0
!      do j=1,natoms
!         k=3*(j-1)
!         do i=1,3
!            if ( coords(k+i) .ne. cl_oldcoords(k+i) ) then
!            !r = sum( (coords(k+1:k+3) - cl_oldcoords(k+1:k+3))**2 )
!            !if ( r .ne. 0.d0 ) then
!               CL_NMOVED = CL_NMOVED + 1
!               CL_MOVED_ATOMS(CL_NMOVED) = j
!               !write(*,*) j, "moved", sqrt(r)
!               exit !exit loop
!            endif
!         enddo
!      enddo
!      !write(*,*) "coords(891) mid  ", coords((891-1)*3+1), cl_oldcoords((891-1)*3+1)
!      cl_oldcoords(1:3*natoms) = coords(1:3*natoms)
!      !write(*,*) "coords(891) end  ", coords((891-1)*3+1), cl_oldcoords((891-1)*3+1)
!   end subroutine which_atom_moved

   subroutine cl_get_neighbors_several(coords, nmove, movelist_in)
      !find all neighbors of all atoms in movelist_in.
      !
      !save the lists in arrays
      !CL_LISTAA, CL_LISTBB, CL_LISTAB
      !
      !We must take care we have only one of each pair.  This is quite easy to
      !add the pair i,j twice if both i and j are in movelist_in
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer, intent(in) :: nmove, movelist_in(nmove)
      integer             :: movelist(nmove)
      integer :: i, j, n, n2

      movelist(:) = movelist_in(:)
      call sort4_int( nmove, movelist ) !sort movelist
      
      CL_NLISTAA = 0
      CL_NLISTBB = 0
      CL_NLISTAB = 0
      do n=1,nmove
         i = movelist(n)
         !find neighbors of i
         !write(*,*) "find neighbors of", i
         if (i .le. ntypeA) then
            !find A-A neighbors
            call make_atom_neighbor_list2(coords, i, hocA, llA, cl_nlistAA, cl_listAA, n, movelist)
            !write(*,*) "AA pairs", cl_nlistaa
            !do n2=1,cl_nlistaa
               !write(*,*) n2, cl_listAA(:,n2)
            !enddo
            !find A-B neighbors
            call make_atom_neighbor_list2(coords, i, hocB, llB, cl_nlistAB, cl_listAB, n, movelist)
         else 
            !find B-B neighbors
            call make_atom_neighbor_list2(coords, i, hocB, llB, cl_nlistBB, cl_listBB, n, movelist)
            !find B-A neighbors
            call make_atom_neighbor_list2(coords, i, hocA, llA, cl_nlistAB, cl_listAB, n, movelist)
         endif
         !write(*,"(A,4I8)") "neighbors", i, cl_nlistaa, cl_nlistbb, cl_nlistab
      enddo

   end subroutine cl_get_neighbors_several

!   subroutine get_neighbors_all(coords)
!      implicit none
!      double precision, intent(in) :: coords(3*natoms)
!      integer i
!      
!      CL_NLISTAA = 0
!      CL_NLISTBB = 0
!      CL_NLISTAB = 0
!      call make_neighbor_list_all(coords, hocA, llA, cl_nlistAA, cl_listAA)
!      call make_neighbor_list_all(coords, hocB, llB, cl_nlistBB, cl_listBB)
!      do i=ntypeA+1,natoms
!         !find the A neighbors of all b atoms
!         call add_to_atom_neighbor_list(coords, i, hocA, llA, cl_nlistAB, cl_listAB)
!      enddo
!
!   end subroutine get_neighbors_all

end module cell_lists_binary_mod

