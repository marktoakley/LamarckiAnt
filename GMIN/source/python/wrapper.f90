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
    call potential(x, grad, energy, .FALSE.,.FALSE.)
end function

function py_get_energy_gradient(x, grad) result(energy)
    use commons, only: natoms
    double precision :: x(3*natoms)
    double precision :: grad(3*natoms)
    double precision energy
    call potential(x, grad, energy, .TRUE.,.FALSE.)
end function

subroutine py_gmin_get_coords(x)
    use commons, only: coords, natoms
    double precision :: x(3*natoms)

    x(1:3*natoms) = coords(:,1)
end subroutine