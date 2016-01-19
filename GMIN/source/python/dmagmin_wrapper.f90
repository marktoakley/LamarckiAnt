function py_get_dof() result(n)
    use genrigid, only: degfreedoms
    integer n
    n = degfreedoms
end function

function py_get_natoms() result(n)
    use commons, only: natoms
    integer n
    n = natoms
end function

function py_get_energy(x) result(energy)
    use commons, only: natoms
    double precision :: x(3*natoms)
    double precision :: grad(3*natoms)
    double precision energy
    call dmacrys_potential(x, grad, energy, .false.,.false.)
end function

function py_get_energy_gradient(x, grad) result(energy)
    use commons, only: natoms
    double precision :: x(3*natoms)
    double precision :: grad(3*natoms)
    double precision energy
    call dmacrys_potential(x, grad, energy, .true.,.false.)
end function

subroutine py_dmagmin_get_coords(x)
    use genrigid
    use commons, only: coords, natoms
    double precision :: x(degfreedoms)

    call transformctorigid(coords(:,1), x)
end subroutine

subroutine py_dmagmin_takestep(x)
    use genrigid
!    use dmacrys_interface
    use commons, only: coords, natoms, tmove, omove
    double precision :: x(degfreedoms)
    integer np
    np=1
    allocate(tmove(1),omove(1))
    call transformrigidtoc(1, nrigidbody, coords(:,1), x)
    call dmacrys_takestep(np)
    print *,"ts done"
    call transformctorigid(coords(:,1), x)
    deallocate(tmove)
    deallocate(omove)
end subroutine
