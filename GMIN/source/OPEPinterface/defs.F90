module defs
  ! This module defines all variables used accross the program ART01
  !
  !  Copyright N. Mousseau, May 2001
  implicit none

  save
  ! ! Numero de version
  double precision, parameter :: version = 2.0 
  ! Lists of parameters
  double precision :: target_temperature
  integer :: NFRAG 
  integer :: NATOMS 
  integer :: VECSIZE
  integer :: VECSIZE1                     ! Length of the force and position vectors of one fragment 

  double precision :: degree_freedom               ! degree of freedom, for nonlinear molecules, degree_freedom = n*3 - 6 - num_constraints

  integer :: NUMBER_EVENTS                ! Total number of events in this run
  character(len=3)  :: SIMULATION_METHOD  ! Type of simulation ART or MD 
  character(len=10)  :: SIMULATION_TYPE    ! serial or replica (exchange) or tempering
  character(len=10) :: REPLICA_TYPE       ! T_Exchange or E_Scale


  logical :: constrained_fragments, restrained_fragments
  double precision :: k_spring

  character(len=5), dimension(:), allocatable :: atomic_type ! Atomic type
  double precision, dimension(:), allocatable, target :: force       ! Working forces on the atoms
  double precision, dimension(:), allocatable, target :: pos         ! Working positions of the atoms
  double precision, dimension(:), allocatable, target :: posref      ! Reference position
  double precision, dimension(:), allocatable, target :: mass        ! masses
  double precision :: force_scaling_factor  ! Factor for rescaling the potential

  double precision, dimension(:), pointer :: x, y, z
  double precision, dimension(:), pointer :: xref, yref, zref
  double precision, dimension(:), pointer :: fx, fy, fz

  integer, dimension(:,:), allocatable :: list_fragments
  
  integer :: mincounter                     ! Counter for output files
  integer :: refcounter                     ! Id of reference file
  integer :: niter                         
  integer :: nitt
  integer :: evalf_number  = 0              ! Number of force evalutions
  integer :: ndigits_filename

  double precision :: total_energy, ref_energy       ! Energies
  character(len=20) :: restart

  logical :: usextc
  logical :: singlefile 
 
  integer :: T_id
  integer :: E_scale
  logical :: init_single_file
  logical :: PBC  ! periodic boundary condition
  logical :: C_M  ! center of mass for writing non-broken chains in pdb file
  double precision :: BL   ! box length 

  logical :: RNA_simulation

  ! Units for writing/reading
  integer, parameter :: FCONF        = 1         
  integer, parameter :: FCOUNTER     = 2        
  integer, parameter :: FLIST        = 3         
  integer, parameter :: FLOG         = 4         
  integer, parameter :: FREFCONFIG   = 11
  integer, parameter :: FENER        = 12
  integer, parameter :: FRESTART     = 13
  integer, parameter :: FREP         = 17

  character(len=20) :: conf_saddle, conf_final, conf_initial
  character(len=20) :: debug_status
  logical :: use_qbug ! If true, calculates the forces the first time and stop

  character(len=20) :: CHAIN_FILE
  character(len=20) :: MASTER_LOGFILE
  character(len=20) :: EVENTSLIST
  character(len=20) :: REFCONFIG
  character(len=3)  :: FINAL
  character(len=3)  :: SADDLE
  character(len=11) :: COUNTER
  character(len=11) :: RESTARTFILE
  character(len=13) :: REPLICAFILE

  character(len=4)  :: PDB_EXT   = '.pdb'
  
  integer :: N_REPLICA      ! number of temperature replica (T exchange)
  integer :: N_E_REPLICA    ! number of energy replica (hamiltonian exchange)
  integer :: N_TEMP         ! number of temperature tempering)

  integer :: n_step_exchange
  double precision, dimension(:), allocatable :: T_replica
  double precision, dimension(:), allocatable :: T_tempering
  double precision, dimension(:), allocatable :: F_tempering
  double precision, dimension(:), allocatable :: E_tempering
  double precision, dimension(:), allocatable :: N_tempering

  integer :: n_step_E_exchange
  integer :: last_id
  double precision, dimension(:), allocatable :: E_replica

  type t_conformations
    character(len=20) :: path
    character(len=20) :: logfile
    integer ::  id, counter
    double precision, dimension(:), allocatable :: pos, posref, vel
    double precision, dimension(:), allocatable :: temperatures,scales
    double precision :: energy, temperature, energyscale,free_energy
  end type t_conformations
  
  integer :: ntasks, taskid, Error
end module defs
