! Called from keywords.f to initialize DMACRYS when specified
! in input file
subroutine dmacrys_setup
    use dmacrys_interface
    print *,'DMACRYS INTERFACE> Using dmacrys, setting up interface'
    CALL dmacrys_initialize
end subroutine

subroutine dmacrys_get_natoms(num_atoms_out)
    use dmacrys_interface
    integer, intent(out) :: num_atoms_out
    print *,"dmacrys_get_natoms",num_atoms
    num_atoms_out = num_atoms
end subroutine

subroutine dmacrys_initialize_gmin
    use dmacrys_interface
    print *,'DMACRYS INTERFACE> Setting up GMIN from DMACRYS structures'
    CALL dmacrys_init_gmin
end subroutine

! calculate the potential using dmacrys
subroutine dmacrys_potential(X,GRAD,EREAL,GRADT)
    use dmacrys_interface
    use commons
    implicit none
    double precision x(1:3*NATOMS)
    double precision grad(1:3*NATOMS)
    double precision ereal
    logical GRADT

    grad(1:3*NATOMS) = 0
    ereal = 0
!    call dmacrys_set_coords(X)
!    EREAL=dmagrys_calc_energy(X, GRAD, GRADT)
    call dmacrys_calc_potential(X,GRAD,EREAL,GRADT)
!    print *,"Energy",ereal
!    if(dmadidstep) then
!        call dmacrys_check_derivatives(3*NATOMS,x,grad)
!    endif
    dmadidstep=.false.
end subroutine

!     vr274> dump lowest configurations in cif file.
subroutine dmacrys_dump_lowest()
    USE QMODULE
    USE DMACRYS_INTERFACE, only : dmacrys_get_lattice,DMACRYS_DUMP_CIF_NO_TRANSFORM
    USE COMMONS, ONLY : NSAVE,NATOMS
    use genrigid, only : get_lattice_matrix

    IMPLICIT NONE
    INTEGER I
    CHARACTER(LEN=10) :: CNUM
    CHARACTER(LEN=16) :: cifname
    DOUBLE PRECISION MLATTICE(3,3)
    print *,"DUMP_LOWEST_DMACRYS"
    DO I=1,NSAVE
        WRITE (CNUM,'(I3.3)') I
        print *,I,"lowest"//trim(cnum)//".cif"
        call dmacrys_reduce_cell(NATOMS, QMINP(I,:))

        call get_lattice_matrix(QMINP(I,3*NATOMS-5:3*NATOMS), MLATTICE)
        WRITE (cifname,'(E16.6)') QMIN(I)
        CALL DMACRYS_DUMP_CIF_NO_TRANSFORM(&
            "lowest"//trim(cnum)//".cif",&
            NATOMS-2,QMINP(I,1:3*NATOMS-6), MLATTICE, "E_"//trim(adjustl(cifname)))
    END DO
end subroutine

subroutine dmacrys_dump(filename, coords)
    USE DMACRYS_INTERFACE, only : dmacrys_get_lattice,DMACRYS_DUMP_CIF_NO_TRANSFORM
    USE COMMONS, ONLY : NATOMS
    use genrigid, only : get_lattice_matrix
    IMPLICIT NONE
    INTEGER I
    DOUBLE PRECISION coords(3*NATOMS)
    CHARACTER*(*) :: filename
    DOUBLE PRECISION MLATTICE(3,3)
    call dmacrys_reduce_cell(NATOMS, coords)
    call get_lattice_matrix(coords(3*NATOMS-5:3*NATOMS), MLATTICE)
    CALL DMACRYS_DUMP_CIF_NO_TRANSFORM(filename,&
            NATOMS-2,coords(1:3*NATOMS-6), MLATTICE, "noname")
end subroutine

subroutine dmacrys_takestep(np)
    use commons
    use genrigid
    use cell_reduction
    use dmacrys_interface, only: dmadidstep, strain_to_matrix, initial_volume
    use rotations
    use vec3
    implicit none

    integer          :: np
    integer          :: j1, j3, j5, i
    double precision :: dprand, random
    double precision :: translate(3), E(6), LTRANSFORM(3,3), lattice(3,3), vmax, volume

    print *,"DMACRYS_TAKESTEP"

!    call dmacrys_genrandom(coords(:,np))
!    return
    np = 1
    call transformctorigid(coords(:,np), rigidcoords)

    do j1 = 1,3*natoms
        coordso(j1,np) = coords(j1,np)
    end do

    ! displace and rotate rigid bodies
    do j1 = 1, nrigidbody
        j3 = 3*j1
        j5 = 3*nrigidbody + j3

        ! translational steps
        if (tmove(np)) then
            do i=1,3
                translate(i) = step(np)*(2.0d0 * dprand() - 1.0d0)
            end do
            rigidcoords(j3-2:j3) = rigidcoords(j3-2:j3) + translate
        endif

        ! rotational steps
        if (omove(np)) then
            call rot_takestep_aa(rigidcoords(j5-2:j5), ostep(np))
        endif
    enddo

    ! expand the lattice
    rigidcoords(6*nrigidbody+1:6*nrigidbody+6) =  DMACRYS_EXPAND* rigidcoords(6*nrigidbody+1:6*nrigidbody+6)

    ! random displacement in lattice
    do j1=1,6
        E(j1) = dprand()
    end do
    E(1:3) = 2d0*E(1:3)-1d0
    E = E * DMACRYS_LATTICE_STEP

    ! construct strain matrix
    call strain_to_matrix(E, LTRANSFORM)

    ! get the lattice vector (matrix)
    call get_lattice_matrix(rigidcoords(6*nrigidbody+1:6*nrigidbody+6), lattice)

    ! apply strain matrix
    do i=1,3
       lattice(:,i) = matmul(LTRANSFORM, lattice(:,i))
    end do

    ! make lattice to upper (lower?) triangular
    call cell_set_lattice(lattice)
    call cell_get_lattice(lattice)

    ! scale the volume to twice the old volume
    volume = dot_product(lattice(:,1), vec_cross(lattice(:,2),lattice(:,3)))
    vmax =  1.5d0*DMACRYS_EXPAND*initial_volume
    if(volume > vmax) then
        lattice = lattice * (vmax / volume)**(1d0/3d0)
    endif

    ! set the lattice coordinates
    call set_lattice_matrix(rigidcoords(6*nrigidbody+1:6*nrigidbody+6), lattice)

    ! transform back to rigid coordinates
    call transformrigidtoc (1,nrigidbody, coords(:,np), rigidcoords)

    dmadidstep=.true.
end subroutine

subroutine dmacrys_reduce_cell_rigid(rcoords,cell_was_changed)
    use genrigid
    use cell_reduction
    use vec3
    use commons, only: natoms
    implicit none
    integer i
    double precision rcoords(DEGFREEDOMS)
    double precision coords(1:3*NATOMS)
    double precision lattice(3,3), tmp(3,3),transform(3,3),lattice_inv(3,3)
    logical cell_was_changed

    cell_was_changed = .false.
    call get_lattice_matrix(rcoords(DEGFREEDOMS-5:DEGFREEDOMS), lattice)

    call cell_set_lattice(lattice)
    call cell_minimum_reduction

    if(.not.cell_has_changed()) return
    cell_was_changed = .true.

    call invert3x3(lattice, lattice_inv)

    tmp = cell_get_transformation()
    call cell_get_lattice(lattice)

    transform = matmul(lattice,matmul(tmp, lattice_inv))
    call TRANSFORMRIGIDTOC (1, NRIGIDBODY, coords, rcoords(1:DEGFREEDOMS))
    do i=1,natoms
        coords(3*i-2:3*i) = matmul(transform, coords(3*i-2:3*i))
    end do
    coords(3*NATOMS-5) = lattice(1,1)
    coords(3*NATOMS-4) = lattice(2,1)
    coords(3*NATOMS-3) = lattice(3,1)
    coords(3*NATOMS-2) = lattice(2,2)
    coords(3*NATOMS-1) = lattice(3,2)
    coords(3*NATOMS-0) = lattice(3,3)

    call TRANSFORMCTORIGID (coords, rcoords(1:DEGFREEDOMS))

    do i=1,3*NRIGIDBODY
        rcoords(i) = rcoords(i) - entier(rcoords(i))
    end do
end subroutine

subroutine dmacrys_reduce_cell(natoms, atomcoords)
    use cell_reduction
    use genrigid
    use vec3
    implicit none
    integer natoms, i
    double precision atomcoords(1:3*NATOMS)
    double precision lattice(3,3)
    double precision transform(3,3)
    double precision lattice_inv(3,3)
    double precision tmp(3,3)
    double precision rcoords(DEGFREEDOMS)

    call get_lattice_matrix(atomcoords(3*natoms-5:3*natoms), lattice)

    call cell_set_lattice(lattice)
    call cell_minimum_reduction

    if(.not.cell_has_changed()) return

    call invert3x3(lattice, lattice_inv)

    tmp = cell_get_transformation()
    call cell_get_lattice(lattice)

    transform = matmul(lattice,matmul(tmp, lattice_inv))
    do i=1,natoms
        atomcoords(3*i-2:3*i) = matmul(transform, atomcoords(3*i-2:3*i))
    end do
    atomcoords(3*natoms-6+1) = lattice(1,1)
    atomcoords(3*natoms-6+2) = lattice(2,1)
    atomcoords(3*natoms-6+3) = lattice(3,1)
    atomcoords(3*natoms-6+4) = lattice(2,2)
    atomcoords(3*natoms-6+5) = lattice(3,2)
    atomcoords(3*natoms-6+6) = lattice(3,3)

    call TRANSFORMCTORIGID (atomcoords, rcoords(1:DEGFREEDOMS))
    do i=1,3*NRIGIDBODY
        rcoords(i) = rcoords(i) - entier(rcoords(i))
    end do
    call TRANSFORMRIGIDTOC (1, NRIGIDBODY, atomcoords, rcoords(1:DEGFREEDOMS))
end subroutine

! generate a random configuration
subroutine dmacrys_genrandom(coords)
    use commons, only : natoms
    use genrigid
    use rotations
    use vec3
    use cell_reduction
    use dmacrys_interface, only: initial_volume, remove_overlap

    implicit none
    double precision, intent(inout) :: coords(3*NATOMS)
    double precision :: rcoords(degfreedoms)
    integer i, offset,j
    double precision :: dprand
    double precision :: lattice(3,3), u(6)
    double precision :: volume1, volume2, volume_new

    do i = 1, nrigidbody
        rcoords(3*i-2) = dprand()
        rcoords(3*i-1) = dprand()
        rcoords(3*i-0) = dprand()

        offset = 3*nrigidbody + 3*i
        rcoords(offset-2:offset) = rot_random_aa()
    end do

    call get_lattice_matrix(coords(3*natoms-5:3*natoms), lattice)
    volume1 = dot_product(lattice(:,1), vec_cross(lattice(:,2),lattice(:,3)))


    ! generate a random lattice
    ! first do a random rectangular box
    do i=1,3
      u(i) =  1.7d0 * dprand() ! the maximum ration for axis length is 4:1
    end do
    lattice(:,1) = (/1d0, u(2), u(3)/)
    lattice(:,2) = (/0d0, 1d0, u(1)/)
    lattice(:,3) = (/0d0, 0d0, 1d0/)

    ! execute cell reduction
    call cell_set_lattice(lattice)
    !call cell_reduce_cell()
    call cell_get_lattice(lattice)

    volume_new = min(volume1 * 1.2d0**3,2d0**3 * initial_volume)

    ! scale the volume to twice the old volume
    volume2 = dot_product(lattice(:,1), vec_cross(lattice(:,2),lattice(:,3)))

    lattice = lattice * (volume_new / volume2)**(1d0/3d0)

    rcoords(degfreedoms-6+1) = lattice(1,1)
    rcoords(degfreedoms-6+2) = lattice(2,1)
    rcoords(degfreedoms-6+3) = lattice(3,1)
    rcoords(degfreedoms-6+4) = lattice(2,2)
    rcoords(degfreedoms-6+5) = lattice(3,2)
    rcoords(degfreedoms-6+6) = lattice(3,3)
!    CALL MYMYLBFGS(N,M,rcoords,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)

    CALL transformrigidtoc(1,nrigidbody,coords,rcoords)

    !call remove_overlap(coords)

end subroutine

subroutine dmacrys_minimize(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
    use commons, only : natoms, dmacryst
    use genrigid
    use dma_minimize, only: minimize
    implicit none

    integer :: n,m,itmax,itdone,np
    double precision :: xcoords(3*natoms), eps,energy
    logical :: diagco, reset,  mflag

    call minimize(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
end subroutine
