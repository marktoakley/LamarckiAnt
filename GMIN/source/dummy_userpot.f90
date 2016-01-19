subroutine userpot_init
    print *,'ERROR: you are using USERPOT with a standard GMIN binary'
    stop
end subroutine

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
!    integer, intent(out) :: num_atoms
    integer :: num_atoms
    print *,'ERROR: you are using USERPOT with a standard GMIN binary'
    stop
end subroutine

! copy all coordinates to gmin coords(:,1) array (in commons)
! TODO: what about MPI?
! for very simple standard cases this is just initial configuration,
! but can be also more complicated, e.g. initalizing generalized rigid
! body framework
subroutine userpot_initialize_gmin(dof, x)
!    integer, intent(in) :: dof
!    double precision, intent(out) :: x(dof)
    integer :: dof
    double precision :: x(dof)
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

! called by gmin to calculate potential
subroutine userpot_potential(dof,X,GRAD,EREAL,GRADT)
!    integer, intent(in) :: dof                 ! number of degrees of freedom
!    double precision, intent(in) :: X(dof)     ! current coordinates
!    double precision, intent(out) :: GRAD(dof) ! gradient
!    double precision, intent(out) :: EREAL     ! energy
!    logical, intent(in) :: gradt              ! is the gradient needed?
    integer :: dof                 ! number of degrees of freedom
    double precision :: X(dof)     ! current coordinates
    double precision :: GRAD(dof) ! gradient
    double precision :: EREAL     ! energy
    logical :: gradt              ! is the gradient needed?
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine userpot_dump_configuration(filename, coords) 
DOUBLE PRECISION COORDS(*)
CHARACTER(*) :: filename   
end subroutine

subroutine userpot_dump_lowest
end subroutine
