module chirality
   use amber12_interface_mod, only: amber12_atom
   use porfuncs
   use commons, only: myunit

   implicit none

   logical                                :: chirality_debug = .false.
   double precision                       :: cis_trans_tol = 30.0d0
! These variables store the S/R and cis/trans states of all relevant atoms
! For chirality: R = .true., S = .false.
! For cis/trans: cis = 'C', trans = 'T', neither = 'N'
   logical, dimension(:), allocatable     :: sr_states
   logical, dimension(:), allocatable     :: sr_states_initial
   character, dimension(:), allocatable   :: cis_trans_states
   character, dimension(:), allocatable   :: cis_trans_states_initial
! These variables store the indices of the atoms to be checked for the
! appropriate type of isomerism. These lists can be generated automatically
! for proteins and nucleic acids using the relevant subroutines below.
! The central atom is first in the sr_atoms list
   integer, dimension(:, :), allocatable  :: sr_atoms
   integer, dimension(:, :), allocatable  :: cis_trans_atoms
! These variables contains the location of the chirality and cis/trans Python
! scripts, if they're not provided by a source build.
   character, dimension(:), allocatable   :: chirality_script
   character, dimension(:), allocatable   :: cis_trans_script

   private :: cross_product, dihedral
   public

contains

function cross_product(vector1, vector2) result(output)
   double precision, dimension(3), intent(in)   :: vector1, vector2
   double precision, dimension(3)               :: output

   output(1) = vector1(2) * vector2(3) - vector1(3) * vector2(2)
   output(2) = vector1(3) * vector2(1) - vector1(1) * vector2(3)
   output(3) = vector1(1) * vector2(2) - vector1(2) * vector2(1)

end function cross_product

function dihedral(coords) result(angle)
   double precision, dimension(12), intent(in)  :: coords
   double precision                             :: angle
   double precision, dimension(9)               :: vectors
   integer                                      :: i
   double precision, dimension(3)               :: b1xb2, b2xb3, b1xb2_x_b2xb3, b2_norm
   double precision                             :: ddot, dnrm2

! i counts through x, y and z
   do i = 1, 3
   ! vectors(1:3) = b1, vectors(4:6) = b2, vectors(7:9) = b3
      vectors(i)     = coords(i+3) - coords(i)
      vectors(i+3)   = coords(i+6) - coords(i+3)
      vectors(i+6)   = coords(i+9) - coords(i+6)
   end do

! calculate cross products of b1 with b2 and b2 with b3
   b1xb2 = cross_product(vectors(1:3), vectors(4:6)) 
   b2xb3 = cross_product(vectors(4:6), vectors(7:9)) 
   b1xb2_x_b2xb3 = cross_product(b1xb2, b2xb3)

! normalise b2
   b2_norm = vectors(4:6) / dnrm2(3, vectors(4:6), 1)

! calculate angle according to formula:
!
! phi = atan2( ([b1 x b2] x [b2 x b3]) . (b2/|b2|), [b1 x b2] . [b2 x b3] )
!
! @article {BlondelKarplus1996,
! author = {Blondel, Arnaud and Karplus, Martin},
! title = {New formulation for derivatives of torsion angles and improper torsion angles in molecular mechanics: Elimination of singularities},
! journal = {Journal of Computational Chemistry},
! volume = {17},
! number = {9},
! publisher = {John Wiley & Sons, Inc.},
! issn = {1096-987X},
! doi = {10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T},
! pages = {1132--1141},
! year = {1996},
! }
   angle = atan2(ddot(3, b1xb2_x_b2xb3, 1, b2_norm, 1), ddot(3, b1xb2, 1, b2xb3, 1))

end function dihedral

function cis_trans(coords, tolerance_in) result(state)
! Works out whether a molecule is cis, trans or neither, to within a certain tolerance
! (in degrees). If tolerance isn't specified, uses the default of 5.0 degrees.
   double precision, dimension(12), intent(in)  :: coords
   double precision, intent(in), optional       :: tolerance_in
   double precision                             :: tolerance
   character                                    :: state
   double precision                             :: angle_degrees
   double precision                             :: pi

! Set tolerance appropriately to 5 degrees if not specified
   if (present(tolerance_in)) then
      tolerance = tolerance_in
   else
      tolerance = cis_trans_tol
   end if


! Calculate the dihedral angle in degrees
   pi = 4.0d0 * atan(1.0d0) 
   angle_degrees = 180.0d0 * (dihedral(coords) / pi)

! If we're debugging, print the angle.
!   if (chirality_debug) then
!      print *, "Angle calculated by cis_trans function: ", angle_degrees, "degrees"
!   end if

! Assign state to be "C" (cis), "T" (trans) or "N" (neither)
   if (abs(angle_degrees) < tolerance) then
      state = "C"
   else if ((180.0d0 - abs(angle_degrees)) < tolerance) then
      state = "T"
   else
      state = "N"
   end if

end function cis_trans

function chirality_sr(coords, centre) result(right_handed)
! Works out whether a molecule is left- or right-handed (i.e. S or R)
! This is calculated by applying the CIP rules to the system and working out the
! dihedral angle between the atoms formed by:
!
! 1 - centre - 4 - 2
!
! It assumes that the coordinates have been given in CIP priority order. Thus, you
! need to ensure that the atoms have been put in the correct order. In general, for
! amino acids this is:
!
! N, carbonyl C, C-beta, H
!
! However, for cysteine, the C-alpha has higher priority groups attached to it 
! (i.e. sulphur) so the priority order changes to:
!
! N, C-beta, carbonyl C, H
!
! Accordingly, most amino acids are S at the C-alpha, but cysteine is R and glycine
! is non-chiral.
!
   double precision, dimension(12), intent(in)  :: coords
   double precision, dimension(3), intent(in)   :: centre
   logical                                      :: right_handed
   double precision                             :: angle
   double precision, dimension(12)              :: dihedral_coords

! Select the set of coordinates we wish to calculate the dihedral angle for
   dihedral_coords(1:3)    = coords(1:3)
   dihedral_coords(4:6)    = centre(:)
   dihedral_coords(7:9)    = coords(10:12)
   dihedral_coords(10:12)  = coords(4:6)

! Calculate the angle
   angle = dihedral(dihedral_coords)

! If we're debugging, print the angle.
   if (chirality_debug) then
      print *, "Angle calculated by chirality_sr function: ", 180 * angle / (4.0d0 * atan(1.0d0)), "degrees"
   end if

! If the angle is positive, this should indicate a right-handed (R) system, if it's
! negative, this indicates a left-handed (S) system.
   right_handed = (angle > 0.0d0)

end function chirality_sr

subroutine init_chiral(coords)
!
! Runs a Python script ($SVN/SCRIPTS/AMBER/chirality/chirality.py) which 
! identifies the chiral centres and their atoms, ranked in order.
!
! Then calculates the initial SR states and stores them in sr_states_initial 
!
! Arguments
! ---------
! coords(in): starting coordinates
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in)   :: coords
! Variables
   integer                                      :: num_chiral_centres
   integer                                      :: i
   integer                                      :: j
   integer                                      :: atom_number
   double precision, dimension(3)               :: centre_coords
   double precision, dimension(12)              :: neighbour_coords
   integer                                      :: file_unit
! Functions
   integer                                      :: file_length

#ifdef _SVN_ROOT_
   call system('python ' // _SVN_ROOT_ // '/SCRIPTS/AMBER/chirality/chirality.py' // ' coords.prmtop')
#else
   call system('python ' // chirality_script // ' coords.prmtop')
#endif

! Work out the number of chiral centres by reading the .chirality_list file 
   num_chiral_centres = file_length('.chirality_list')
   if (.not. allocated(sr_atoms)) allocate(sr_atoms(num_chiral_centres, 5))

! Now read the chiral centres into sr_atoms
   call file_open('.chirality_list', file_unit, .false.)
   do i = 1, num_chiral_centres
      read(file_unit, '(5i8)') sr_atoms(i, :)
   end do
   close(file_unit)

! Print a test copy
!   call file_open('chirality_list_copy', file_unit, .false.)
!   do i = 1, num_chiral_centres
!      write(file_unit, '(5i10)') sr_atoms(i, :)
!   end do
!   close(file_unit)

! Now calculate the chirality of the centres and save it in sr_states_initial
   if (.not. allocated(sr_states_initial)) allocate(sr_states_initial(num_chiral_centres))
   do i = 1, num_chiral_centres
      atom_number = sr_atoms(i, 1)
      centre_coords(1) = coords(3 * atom_number - 2)
      centre_coords(2) = coords(3 * atom_number - 1)
      centre_coords(3) = coords(3 * atom_number    )
      do j = 1, 4
         atom_number = sr_atoms(i, j + 1) 
         neighbour_coords(3 * j - 2) = coords(3 * atom_number - 2)
         neighbour_coords(3 * j - 1) = coords(3 * atom_number - 1)
         neighbour_coords(3 * j    ) = coords(3 * atom_number    )
      end do
      sr_states_initial(i) = chirality_sr(neighbour_coords, centre_coords)
   end do 

! Print the initial chiral states
   call file_open('initial_chiral_states', file_unit, .false.)
   do i = 1, num_chiral_centres
      write(file_unit, '(5I10)') sr_atoms(i, :)
      write(file_unit, '(L2)') sr_states_initial(i)
   end do
   close(file_unit)

end subroutine init_chiral

subroutine init_cis_trans(coords)
!
! Runs a Python script ($SVN/SCRIPTS/AMBER/chirality/cistrans.py) which 
! identifies peptide bonds.
!
! Then calculates the initial cis/trans states and stores them in cis_trans_states_initial.
!
! Arguments
! ---------
! coords(in): starting coordinates
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in)   :: coords
! Variables
   integer                                      :: num_peptide_bonds
   integer                                      :: i
   integer                                      :: j
   integer                                      :: atom_number
   double precision, dimension(12)              :: peptide_coords
   integer                                      :: file_unit
! Functions
   integer                                      :: file_length

#ifdef _SVN_ROOT_
   call system('python ' // _SVN_ROOT_ // '/SCRIPTS/AMBER/chirality/cistrans.py' // ' coords.prmtop')
#else
   call system('python ' // cis_trans_script // ' coords.prmtop')
#endif

! Work out the number of peptide bonds by reading the .cis_trans_list file 
   num_peptide_bonds = file_length('.cis_trans_list')
   if (.not. allocated(cis_trans_atoms)) allocate(cis_trans_atoms(num_peptide_bonds, 4))

! Now read the chiral centres into cis_trans_atoms
   call file_open('.cis_trans_list', file_unit, .false.)
   do i = 1, num_peptide_bonds
      read(file_unit, '(4i8)') cis_trans_atoms(i, :)
   end do
   close(file_unit)

! Print a test copy
!   call file_open('cis_trans_list_copy', file_unit, .false.)
!   do i = 1, num_peptide_bonds
!      write(file_unit, '(4i8)') cis_trans_atoms(i, :)
!   end do
!   close(file_unit)

! Now calculate the isomerism of the peptide bonds and save it in cis_trans_states_initial
   if (.not. allocated(cis_trans_states_initial)) allocate(cis_trans_states_initial(num_peptide_bonds))
   do i = 1, num_peptide_bonds
      do j = 1, 4
         atom_number = cis_trans_atoms(i, j) 
         peptide_coords(3 * j - 2) = coords(3 * atom_number - 2)
         peptide_coords(3 * j - 1) = coords(3 * atom_number - 1)
         peptide_coords(3 * j    ) = coords(3 * atom_number    )
      end do
      cis_trans_states_initial(i) = cis_trans(peptide_coords)
   end do 

! Print the initial cis trans states
   call file_open('initial_cis_trans_states', file_unit, .false.)
   do i = 1, num_peptide_bonds
      write(file_unit, '(5I10)') cis_trans_atoms(i, :)
      write(file_unit, '(A2)') cis_trans_states_initial(i)
   end do
   close(file_unit)

end subroutine init_cis_trans

subroutine chirality_check(coords, same_chirality)
!
! Compares the chirality of the atoms in sr_atoms with the states recorded in
! sr_states_initial. Returns .false. if the chirality has changed and .true.
! if it is the same.
!
! Arguments
! ---------
! coords(in): starting coordinates
! same_chirality(out): .true. if the chiral states are the same
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in)   :: coords
   logical, intent(out)                         :: same_chirality
! Variables
   integer                                      :: num_chiral_centres
   integer                                      :: i
   integer                                      :: j
   integer                                      :: atom_number
   double precision, dimension(3)               :: centre_coords
   double precision, dimension(12)              :: neighbour_coords
   integer                                      :: file_unit
   logical, dimension(:), allocatable           :: same_chirality_array

! The number of chiral centres should be equal to the number of entries in
! sr_states_initial and the first dimension of sr_atoms.
   num_chiral_centres = size(sr_states_initial)
! Do the sanity check.
   if (num_chiral_centres /= size(sr_atoms, dim = 1)) then
      stop 'Error in chirality_check: num_chiral_centres /= size(sr_atoms, dim = 1)'
   end if

! Calculate the chirality of the centres and save it in sr_states
   if (.not. allocated(sr_states)) allocate(sr_states(num_chiral_centres))
   do i = 1, num_chiral_centres
      atom_number = sr_atoms(i, 1)
      centre_coords(1) = coords(3 * atom_number - 2)
      centre_coords(2) = coords(3 * atom_number - 1)
      centre_coords(3) = coords(3 * atom_number    )
      do j = 1, 4
         atom_number = sr_atoms(i, j + 1) 
         neighbour_coords(3 * j - 2) = coords(3 * atom_number - 2)
         neighbour_coords(3 * j - 1) = coords(3 * atom_number - 1)
         neighbour_coords(3 * j    ) = coords(3 * atom_number    )
      end do
      sr_states(i) = chirality_sr(neighbour_coords, centre_coords)
   end do 

! Print the chiral states
!   call file_open('chiral_states', file_unit, .true.)
!   write(file_unit, '(10l2)') sr_states(:)
!   close(file_unit)

! Compare the values in sr_states to those in sr_states_initial. If all the
! values are the same, set same_chirality to .true.
   if (.not. allocated(same_chirality_array)) allocate(same_chirality_array(num_chiral_centres))
   same_chirality_array = (sr_states .eqv. sr_states_initial)
   same_chirality = all(same_chirality_array) 

! For debugging, print out chiral states that invert.
   do i = 1, num_chiral_centres
      if (.not. same_chirality_array(i)) then
         write(myunit, '(a23, i6, a20)') 'chirality_check> Atom', sr_atoms(i, 1), 'inverted.'
         if (sr_states_initial(i) .eqv. .true.) then
            write(myunit, '(a42)') 'chirality_check> Originally R chirality.'
         else
            write(myunit, '(a42)') 'chirality_check> Originally S chirality.'
         end if            
      end if
   end do

end subroutine chirality_check

subroutine cis_trans_check(coords, same_cis_trans)
!
! Runs a Python script ($SVN/SCRIPTS/AMBER/chirality/cistrans.py) which 
! identifies peptide bonds.
!
! Then calculates the initial cis/trans states and stores them in cis_trans_states_initial.
!
! Arguments
! ---------
! coords(in): starting coordinates
!
   implicit none
! Arguments
   double precision, dimension(:), intent(in)   :: coords
   logical, intent(out)                         :: same_cis_trans
! Variables
   integer                                      :: num_peptide_bonds
   integer                                      :: i
   integer                                      :: j
   integer                                      :: atom_number
   double precision, dimension(12)              :: peptide_coords
   integer                                      :: file_unit
   logical, dimension(:), allocatable           :: same_cis_trans_array

! The number of peptide bonds should be equal to the number of entries in
! cis_trans_states_initial and the first dimension of cis_trans_atoms.
   num_peptide_bonds = size(cis_trans_states_initial)
! Do the sanity check.
   if (num_peptide_bonds /= size(cis_trans_atoms, dim = 1)) then
      stop 'Error in chirality_check: num_peptide_bonds /= size(cis_trans_atoms, dim = 1)'
   end if

! Now calculate the isomerism of the peptide bonds and save it in cis_trans_states_initial
   if (.not. allocated(cis_trans_states)) allocate(cis_trans_states(num_peptide_bonds))
   do i = 1, num_peptide_bonds
      do j = 1, 4
         atom_number = cis_trans_atoms(i, j) 
         peptide_coords(3 * j - 2) = coords(3 * atom_number - 2)
         peptide_coords(3 * j - 1) = coords(3 * atom_number - 1)
         peptide_coords(3 * j    ) = coords(3 * atom_number    )
      end do
      cis_trans_states(i) = cis_trans(peptide_coords)
   end do 

! Print the cis/trans states
!   call file_open('cis_trans_states', file_unit, .true.)
!   write(file_unit, '(10A2)') cis_trans_states(:)
!   close(file_unit)

! Compare the values in cis_trans_states to those in cis_trans_states_initial. If all the
! values are the same, set same_cis_trans to .true.
   if (.not. allocated(same_cis_trans_array)) allocate(same_cis_trans_array(num_peptide_bonds))
   same_cis_trans_array = (cis_trans_states == cis_trans_states_initial)
   same_cis_trans = all(same_cis_trans_array) 

! For debugging, print out cis/trans states that change.
   do i = 1, num_peptide_bonds
      if (.not. same_cis_trans_array(i)) then
         if (cis_trans_states(i) == 'C') then
            write(myunit, '(a31, 4i6, a11, a5)') 'cis_trans_check> Peptide bond', cis_trans_atoms(i, :), 'changed to', 'cis.'
         else if (cis_trans_states(i) == 'T') then
            write(myunit, '(a31, 4i6, a11, a7)') 'cis_trans_check> Peptide bond', cis_trans_atoms(i, :), 'changed to', 'trans.'
         else
            write(myunit, '(a31, 4i6, a11, a9)') 'cis_trans_check> Peptide bond', cis_trans_atoms(i, :), 'changed to', 'neither.'
         end if
         if (cis_trans_states_initial(i) == 'C') then
            write(myunit, '(a34)') 'cis_trans_check> Originally cis.'
         else if (cis_trans_states_initial(i) == 'T') then
            write(myunit, '(a36)') 'cis_trans_check> Originally trans.'
         else
            write(myunit, '(a38)') 'cis_trans_check> Originally neither.'
         end if
      end if
   end do

end subroutine cis_trans_check

! =================================================================================================
! =================================================================================================
!
! END OF FUNCTIONALITY, START OF TESTS
!
! =================================================================================================
! =================================================================================================

subroutine test_topology_read()
   use amber12_interface_mod, only: amber12_atom, amber12_residue
   use graph_mod
   implicit none

   type(amber12_atom), dimension(:), allocatable      :: atoms
   type(amber12_residue), dimension(:), allocatable   :: residues

   integer, dimension(:,:), allocatable               :: adjacency
   integer                                            :: root_atom
   integer, dimension(:,:), allocatable               :: paths
   integer, dimension(:,:), allocatable               :: chiral
   integer, dimension(:,:), allocatable               :: cis_trans

   integer                                            :: i
   integer                                            :: j

   print *, 'I currently do nothing.'

end subroutine test_topology_read

subroutine test_chirality()
! series of tests for this module which should pass appropriately
   double precision                 :: pi
   double precision, dimension(12)  :: plane_coords_1, plane_coords_2
   double precision, dimension(12)  :: right_angle_coords_1, right_angle_coords_2
   double precision, dimension(12)  :: alanine_coords
   double precision, dimension(3)   :: alanine_centre
   double precision, dimension(12)  :: alanine_cis_trans
   double precision                 :: output
   logical                          :: right_handed

   plane_coords_1 = (/ 2.0, 1.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 1.5, 5.5 /) 
   plane_coords_2 = (/ 2.0, 1.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, -1.5, 5.5 /) 
   right_angle_coords_1 = (/ 2.0, 1.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0, -1.0, 0.0, 2.0 /) 
   right_angle_coords_2 = (/ 2.0, 1.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0, 4.0, 0.0, 2.0 /) 

! Alanine S/R coords from the following pdb:
!================================================================================
! ATOM      7  N   ALA     2       3.555   3.970  -0.000  1.00  0.00           N
! ATOM      9  CA  ALA     2       2.404   4.850  -0.000  1.00  0.00           C
! ATOM     10  HA  ALA     2       1.803   4.663  -0.890  1.00  0.00           H
! ATOM     11  CB  ALA     2       1.536   4.618   1.232  1.00  0.00           C
! ATOM     15  C   ALA     2       2.831   6.311  -0.000  1.00  0.00           C
!================================================================================
   alanine_coords = (/ 3.555, 3.970, 0.000, 2.831, 6.311, 0.000, 1.536, 4.618, 1.232, 1.803, 4.663, -0.890 /)
   alanine_centre = (/ 2.404, 4.850, 0.000 /)

   right_handed = chirality_sr(alanine_coords, alanine_centre)
   if (right_handed) then
      print *, "Alanine: R"
   else
      print *, "Alanine: S"
   end if

! Alanine cis/trans coords from the following pdb:
!================================================================================
! ATOM      9  CA  ALA     2       2.404   4.850  -0.000  1.00  0.00           C
! ATOM     15  C   ALA     2       2.831   6.311  -0.000  1.00  0.00           C
! ATOM     17  N   ALA     3       1.854   7.219  -0.000  1.00  0.00           N
! ATOM     19  CA  ALA     3       2.130   8.642  -0.000  1.00  0.00           C
!================================================================================

   alanine_cis_trans = (/ 2.404, 4.850, 0.000, 2.831, 6.311, 0.000, 1.854, 7.219, 0.000, 2.130, 8.642, 0.000 /)
   print *, "Alanine cis/trans: ", cis_trans(alanine_cis_trans)

! assign pi
   pi = 4.0d0 * atan(1.0d0)

   print *, "Angle:" , 180.0 * dihedral(plane_coords_1) / pi
   print *, "Angle:" , 180.0 * dihedral(plane_coords_2) / pi
   print *, "Angle:", 180.0 * dihedral(right_angle_coords_1) / pi
   print *, "Angle:", 180.0 * dihedral(right_angle_coords_2) / pi

   print *, "Right angle cis/trans: ", cis_trans(right_angle_coords_1)

end subroutine test_chirality

end module chirality
