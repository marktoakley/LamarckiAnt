subroutine userpot_init
	use main_mod
	use list_mod
	use init_mod
	implicit none
	
	call init_lattice(.false.,'quad')
	return
end subroutine userpot_init	

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
	use main_mod
	implicit none
	integer, intent(out) :: num_atoms
	
	num_atoms = nats
	return
end subroutine

! copy all coordinates to gmin coords(:,1) array (in commons)
! TODO: what about MPI?
! for very simple standard cases this is just initial configuration,
! but can be also more complicated, e.g. initalizing generalized rigid
! body framework
subroutine userpot_initialize_gmin(dof,x)
	use init_mod
	implicit none
	integer, intent(in) :: dof
	double precision, intent(out) :: x(dof)
	
	call init_shiftspins(x)
	return
end subroutine

! called by gmin to calculate potential
subroutine userpot_potential(dof,X,GRAD,EREAL,GRADT)
	use main_mod
	use force_mod
	integer, intent(in) :: dof                 ! number of degrees of freedom
	double precision, intent(in) :: X(dof)     ! current coordinates
	double precision, intent(out) :: GRAD(dof) ! gradient
	double precision, intent(out) :: EREAL     ! energy
	logical, intent(in) :: gradt              ! is the gradient needed?
	
	if (dof .ne. n) then
		print *,'Failure in initializing userpot_potential!'
		stop
	endif
	call get_energy(x,ereal)
	print *,x
	print *,"energy", ereal
	if (gradt) then
		call get_gradient(x,grad)
	endif
	return
end subroutine

subroutine userpot_dump_configuration(filename, coords)
end subroutine

subroutine userpot_dump_lowest
end subroutine

