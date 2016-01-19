subroutine dmacrys_setup
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_get_natoms(num_atoms_out)
    integer :: num_atoms_out
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_initialize_gmin
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

! calculate the potential using dmacrys
subroutine dmacrys_potential(X,GRAD,EREAL,GRADT)
    DOUBLE PRECISION X(*), GRAD(*), EREAL
    LOGICAL GRADT
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_dump_lowest()
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_takestep(np)
    integer :: np
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_dump(filename, coords)
    DOUBLE PRECISION :: COORDS(:)
    CHARACTER(*) :: filename
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_reduce_cell_rigid(coords,cell_has_changed)
    double precision :: coords(:)
    integer :: cell_has_changed
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_genrandom(coords)
    DOUBLE PRECISION COORDS(*)
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine dmacrys_minimize(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
    integer :: n,m,itmax,itdone,np
    double precision :: xcoords(:), eps,energy
    logical :: diagco, reset,  mflag

    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine
