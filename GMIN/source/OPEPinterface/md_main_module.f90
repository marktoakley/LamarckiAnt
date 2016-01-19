! This set of routines serves for the initialisation 
!
! Copyright Normand Mousseau January 2006

module md_defs
  use defs
  use geometric_corrections

  implicit none
  save

  double precision :: timeunit, dt, dt1, dt2, timestep
  double precision, dimension(:), allocatable, target :: vel
  double precision, dimension(:), allocatable :: invmass
  double precision, dimension(:), allocatable :: dt1_mass, dt2_mass
  double precision :: friction_coef
  double precision :: langevin_scalev
  double precision, dimension(:), allocatable :: temperatures
  double precision, dimension(:), pointer :: vx, vy, vz
  double precision :: initial_simulation_time
  integer :: restart_simulation_time

  logical :: confining_sphere
  double precision :: radius_confining_sphere,radius_confining_sphere2, k_confine

  integer :: n_production, n_correct_com, n_correct_rotation
  integer :: n_stats, n_stats_average, berendsen_time
  integer :: n_save_configs, n_equilibration, n_rescale_v_equil
  integer :: n_steps_thermalization, n_save_restarts
  double precision :: simulation_time
  character(len=20) :: thermostat
  character(len=20) :: LOGFILE

  double precision :: avstat, avpot, avkin, avetot
  double precision :: avtemp, avpot2, avkin2, avetot2, avtemp2
  
  logical :: save_unrotated_structures
  
  ! Gaussian number
  integer :: gaussian_flag
  double precision :: gaussian_number2

  ! for rattle
  logical :: rattle
  integer :: nrattle
  integer :: nbondh, nbonda, nbondt, nbonds

  ! For contrained motion of hydrogen
  logical :: control_hydrogen

  double precision, allocatable, dimension(:) :: beq, beq2, ibeq2, rij2
  double precision, allocatable, dimension(:,:) :: rij1
  double precision, allocatable, dimension(:) :: redu1_mass, redu2_mass, reduA, reduB, reduAA, reduBB
  integer, allocatable, dimension(:) :: ia, ib
  double precision :: epspos, epsvel
 
  logical :: thermo
  logical :: lj_test
  logical :: force_calc_test
  logical :: chk_ene
  double precision, allocatable, dimension(:) :: tot_ene

  contains

  ! Routine that interfaces with OPEP
  subroutine calcforce(scale,xa,fa,etot)
    use calcforces
    double precision, intent(in) :: scale
    double precision, intent(out) :: etot
    double precision, dimension(vecsize), intent(in) :: xa
    double precision, dimension(vecsize), intent(out) :: fa

    ! Call the force field
    if (RNA_simulation) then
      call calcforce_RNA(scale,xa,fa,etot,.false., 0.0d0)
    else
      call calcforce_protein(scale,xa,fa,etot)
    endif
 
  end subroutine calcforce
  

  ! Routine that interfaces with OPEP
  ! this one is modified to allow for logging of the energy terms
  subroutine calcforce_log(scale,xa,fa,etot,logener, simuPercent)
    use calcforces
    double precision, intent(in) :: scale
    double precision, intent(out) :: etot
    double precision, dimension(vecsize), intent(in) :: xa
    double precision, dimension(vecsize), intent(out) :: fa
    logical :: logener
    double precision simuPercent

    ! Call the force field
    if (RNA_simulation) then
        call calcforce_RNA(scale,xa,fa,etot,logener, simuPercent)
    else
      call calcforce_protein(scale,xa,fa,etot)
    endif
 
  end subroutine calcforce_log

end module md_defs
