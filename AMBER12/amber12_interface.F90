

module amber12_interface_mod
  use iso_c_binding

  implicit none

!******************************************************************************
! Types and parameters for representing AMBER12 atoms.
  type amber12_atom
    character (len = 4)    :: name
    double precision       :: charge
    double precision       :: mass
    integer                :: res_index
    character (len = 4)    :: res_name
    integer                :: type_index
!    character (len = 4)    :: atom_type
!    character (len = 4)    :: tree_class
    double precision       :: gb_radii
    double precision       :: gb_screen
    integer, dimension(4)  :: bonded_atoms
    integer                :: num_bonds
  end type amber12_atom

  integer, parameter                   :: amber12_atom_size = 10
!******************************************************************************

!******************************************************************************
! Types and parameters for representing AMBER12 residues.
  type amber12_residue
    character (len = 4)                :: name
    integer                            :: num_atoms
    type(amber12_atom), dimension(100) :: atoms
    integer                            :: start_index
    integer                            :: end_index
  end type amber12_residue

  integer, parameter                   :: amber12_residue_size = 5
!******************************************************************************

!******************************************************************************
! Types and parameters for writing pdb ATOM coordinate records to pdb files.
  type pdb_atom_data
    character (len = 6)    :: record_name
    integer                :: atom_number
    character (len = 4)    :: atom_name
    character (len = 1)    :: alt_loc_indicator
    character (len = 3)    :: residue_name
    character (len = 1)    :: chain_id
    integer                :: residue_number
    character (len = 1)    :: insertion_code
    double precision       :: x_coord
    double precision       :: y_coord
    double precision       :: z_coord
    double precision       :: occupancy
    double precision       :: b_factor
    character (len = 2)    :: element
    character (len = 2)    :: charge
  end type pdb_atom_data

! This defines the format to be used when writing ATOM lines for PDB files.    
  integer, parameter                   :: pdb_atom_data_size = 15
  type(pdb_atom_data), parameter       :: null_pdb_atom_data = &
    pdb_atom_data('ATOM  ',0,'    ',' ','   ',' ',0,' ',0.d0,0.d0,0.d0,1.d0,&
                 &0.d0,'  ','  ')
  character (len=*), parameter         :: atom_string_format = &
  &'(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2,10x,a2,a2)'
!******************************************************************************

!******************************************************************************
! Subset of gb_pot_ene_rec and pme_pot_ene_rec to interoperate with C and GMIN.
  type, bind(c) :: pot_ene_rec_c
    real(c_double)    :: total
    real(c_double)    :: vdw_tot
    real(c_double)    :: elec_tot
    real(c_double)    :: gb       ! Only GB
    real(c_double)    :: surf     ! Only GB
    real(c_double)    :: hbond    ! Only PME
    real(c_double)    :: bond
    real(c_double)    :: angle
    real(c_double)    :: dihedral
    real(c_double)    :: vdw_14
    real(c_double)    :: elec_14
    real(c_double)    :: restraint
    real(c_double)    :: angle_ub !CHARMM Urey-Bradley Energy
    real(c_double)    :: imp      !CHARMM Improper Energy
    real(c_double)    :: cmap     !CHARMM CMAP
  end type pot_ene_rec_c

!******************************************************************************

  double precision, allocatable        :: amber12_masses(:)
  type(amber12_atom), allocatable      :: amber12_atoms(:)
  type(amber12_residue), allocatable   :: amber12_residues(:)

contains

subroutine amber12_setup(natoms) bind(C, name='amber12_setup')
!
! Initialises all the stuff needed for AMBER to run... 
!
! Arguments
! ---------
!
! natoms(out): number of atoms from topology file
!
  use iso_c_binding, only: c_int
  implicit none
! Arguments
  integer(c_int), intent(out) :: natoms
! Variables  
! If we haven't compiled AMBER12, set natoms to 0.
  natoms = 0

end subroutine amber12_setup

subroutine amber12_get_coords(natoms, coords) bind(C, name='amber12_get_coords')
!
! Places the coordinates detected by AMBER12 into the coords array provided.
!
! Arguments
! ---------
!
! natoms(in): number of atoms
! coords(out): array into which to write coords
!
  use iso_c_binding, only : c_int, c_double

  implicit none
! Arguments
  integer(c_int), intent(in)              :: natoms
  real(c_double), intent(out)             :: coords(3 * natoms)

end subroutine amber12_get_coords

subroutine amber12_energy_and_gradient(natoms, coords, energy, gradients, decomposed_energy) &
           bind(C, name='amber12_energy_and_gradient')
!
! Calculates the energy and gradients using the GB potential.
!
! Arguments
! ---------
!
! natoms (in): number of atoms
! coords (in): atomic coordinates
! energy (out): the total energy of the system
! gradients (out): the gradient of the energy wrt each atomic coordinate
! decomposed_energy (out): the energy of the system decomposed into its 
! components (see the gb_force_mod module for information on the available 
! components)
!

  use iso_c_binding, only: c_int, c_double

  implicit none

! Arguments 
  integer(c_int), intent(in)              :: natoms
  real(c_double), intent(in)              :: coords(3 * natoms)
  real(c_double), intent(out)             :: energy
  real(c_double), intent(out)             :: gradients(3 * natoms)
  type(pot_ene_rec_c), intent(out)        :: decomposed_energy
end subroutine amber12_energy_and_gradient

subroutine amber12_finish() bind(C, name='amber12_finish')
end subroutine amber12_finish

subroutine amber12_write_restart(coords, rst_name, rst_name_length) bind(C, name='amber12_write_restart')
!
! Writes a restart file containing the coordinates provided in the coords array.
!
! Arguments
! ---------
!
! coords(in): coordinates of the atoms
! rst_name(in): file name for the output restart file.
! rst_name_length(in): length of the file name
!
! The restart file format contains atomic coordinates and velocities. In this 
! case, we use the built-in AMBER routine for writing the restart file, as this
! takes care of everything. We provide a file name, coordinates and a velocity
! array full of zeros.
!
! The velocity array hopefully shouldn't matter, because we have imin != 0,
! for which velocities will not be written. However, they are included to
! ensure that we don't fail any checks which write_restart may perform.
!
! This used to use Fortran strings (i.e. character (len=*)), rather than
! character arrays. However, to interoperate with C for Pele, we want to use
! character arrays. This means that we must also specify a maximum length for
! the input string and take caution when dealing with C versus Fortran arrays
! (e.g. use of the C null character \0, which terminates C strings). 
!
  use iso_c_binding, only : c_int, c_double, c_char 

  implicit none
! Arguments
  real(c_double), intent(in)              :: coords(1)
  integer(c_int), intent(in)              :: rst_name_length
  character(kind=c_char, len=1), &
            dimension(rst_name_length), &
            intent(in)                    :: rst_name
end subroutine amber12_write_restart

subroutine amber12_write_pdb(coords, pdb_name, pdb_name_length) bind(C, name='amber12_write_pdb')
!
! Writes a PDB file containing the coordinates provided in the coords array.
!
! Arguments
! ---------
!
! coords(in): coordinates of the atoms
! pdb_name(in): file name for the output PDB file
! pdb_name_length(in): length of the file name
!
! The PDB file format for atomic coordinates is as follows, relevant types are 
! defined at the top of the amber12_interface_mod module:
!
! COLUMNS  DATA  TYPE    FIELD       DEFINITION
! ------------------------------------------------------------------------------
!  1 -  6  Record name   "ATOM  "
!  7 - 11  Integer       serial      Atom  serial number.
! 12       Blank         ------
! 13 - 16  Atom          name        Atom name.
! 17       Character     altLoc      Alternate location indicator.
! 18 - 20  Residue name  resName     Residue name.
! 21       Blank         ------
! 22       Character     chainID     Chain identifier.
! 23 - 26  Integer       resSeq      Residue sequence identifier.
! 27       AChar         iCode       Code for insertion of residues.
! 28 - 30  Blank         ------
! 31 - 38  Real(8.3)     x           Orthogonal coordinates for X in Angstroms.
! 39 - 46  Real(8.3)     y           Orthogonal coordinates for Y in Angstroms.
! 47 - 54  Real(8.3)     z           Orthogonal coordinates for Z in Angstroms.
! 55 - 60  Real(6.2)     occupancy   Occupancy.
! 61 - 66  Real(6.2)     tempFactor  Temperature  factor.
! 67 - 76  Blank         ------
! 77 - 78  LString(2)    element     Element symbol, right-justified.
! 79 - 80  LString(2)    charge      Charge  on the atom.
!
! Also, take note of the notes on string/character arrays as mentioned in the
! amber12_write_restart subroutine above.
!

  use iso_c_binding, only : c_int, c_double, c_char

  implicit none
! Arguments
  real(c_double), intent(in)              :: coords(1)
  integer(c_int), intent(in)              :: pdb_name_length
  character(kind=c_char, len=1), &
            dimension(pdb_name_length), &
            intent(in)                    :: pdb_name
end subroutine amber12_write_pdb

subroutine amber12_write_xyz(coords, xyz_name, xyz_name_length, header) bind(C, name='amber12_write_xyz')
!
! Writes a xyz file containing the coordinates provided in the coords array.
!
! Arguments
! ---------
!
! coords(in): coordinates of the atoms
! xyz_name(in): file name for the output xyz file
! xyz_name_length(in): length of the file name
! header(in): whether or not to write a header for the file 
!
! The xyz file format contains a single line with the number of atoms, then a 
! line containing a comment, followed by a line per atom containing atomic 
! names and coordinates in the following format:
!
! COLUMNS  DATA  TYPE    FIELD       DEFINITION
! ------------------------------------------------------------------------------
!  1 -  4  String        name        Atom name
!  5 - 24  Real (20.10)  x           Atom x-coordinate 
! 25 - 44  Real (20.10)  y           Atom y-coordinate 
! 45 - 64  Real (20.10)  z           Atom z-coordinate 
! 
! The xyz file format is used primarily for OPTIM jobs.
!
! For more information, see Wikipedia:
! http://en.wikipedia.org/wiki/XYZ_file_format
!
! Also, take note of the notes on string/character arrays as mentioned in the
! amber12_write_restart subroutine above.
!

  use iso_c_binding, only : c_int, c_double, c_char, c_bool

  implicit none
! Arguments
  real(c_double), intent(in)              :: coords(1)
  integer(c_int), intent(in)              :: xyz_name_length
  character(kind=c_char, len=1), &
            dimension(xyz_name_length), &
            intent(in)                    :: xyz_name
  logical(c_bool), intent(in)             :: header
end subroutine amber12_write_xyz

subroutine populate_atom_data(atoms, residues)
!
! Populates arrays containing the atoms and residues in the amber12_atom and 
! amber12_residue formats defined above.
!
! Arguments
! ---------
!
! atoms(out): an allocatable array containing the atoms in amber12_atom format
! residues(out): an allocatable array containing the residues in amber12_residue
! format
!


  implicit none
! Arguments
  type(amber12_atom), allocatable, intent(out)     :: atoms(:)
  type(amber12_residue), allocatable, intent(out)  :: residues(:)
end subroutine populate_atom_data
      
subroutine amber12_num_hess(natoms, coords, delta, hessian)
   use iso_c_binding, only : c_int, c_double
   implicit none
! Arguments
   real(c_double), intent(in)    :: coords(1)
   real(c_double), intent(out)   :: hessian(1, 1)
   integer(c_int), intent(in)    :: natoms
   real(c_double), intent(in)    :: delta
   
   write(*, *) "This binary is not compiled to use AMBER12."
   stop

end subroutine amber12_num_hess

end module amber12_interface_mod
