! this is a routine to minimize dmacrys input. it basically is a wrapper for mylbfgs, which first tries
! to minimize the system, and if it gets cold fusion, first tries to turn of multipoles, ... to be continued when running

module dma_minimize
contains
    subroutine minimize(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
        use commons, only : natoms
        use genrigid
        implicit none

        integer :: n,m,itmax,itdone,np
        double precision :: xcoords(3*natoms), eps,energy
        logical :: diagco, reset,  mflag

        double precision :: xrigidcoords(degfreedoms)
        ! first transform to rigid bodies
        call transformctorigid(xcoords, xrigidcoords)

        ! now go for minimize_rigid
        call minimize_rigid(n,m,xrigidcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)

        ! transform back
        !xrigidcoords(1:degfreedoms) = xcoords(1:degfreedoms)
        call transformrigidtoc(1,nrigidbody, xcoords, xrigidcoords)
        !print *,"after transform",xcoords(3*natoms-5:)
    end subroutine

    subroutine minimize_rigid(n,m,xrigidcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
        use commons, only : natoms, coldfusion
        use genrigid
        use dmacrys_interface
        implicit none

        logical tmp

        integer :: n,m,itmax,itdone,np
        double precision :: xcoords(3*natoms), gtmp(3*natoms),eps,energy
        logical :: diagco, reset,  mflag
        double precision :: xrigidcoords(degfreedoms)
        integer cycle

        tmp = ATOMRIGIDCOORDT
        ATOMRIGIDCOORDT = .FALSE.
        cycle=0
        xcoords(1:degfreedoms) = xrigidcoords(1:degfreedoms)
        xcoords(degfreedoms+1:3*natoms) = 0.0d0
        do cycle = 1,5
!            print *,"Cycle",cycle
            ! always start with initial coordinates for each cycle
            COLDFUSION=.FALSE.

            ! the first minimization has failed, try without multipoles
            if(cycle == 2) then
                !interaction_flags=IBCLR(interaction_flags, flag_shortrange)
                !interaction_flags=2
                use_remove_overlap_potential=.true.
                !call dmacrys_dump_cif("1.cif", xcoords(1:degfreedoms))
                call run_lbfgs(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset,np)
                use_remove_overlap_potential=.false.
!                print *, "overlap", energy
                !interaction_flags=IBSET(interaction_flags, flag_shortrange)
                !interaction_flags=63
                !remove_overlap(xcoords)
!               call dmacrys_potential(xcoords, gtmp, energy,.false.)
!                print *,energy
!                call dmacrys_dump_cif("bla.cif", xcoords(1:degfreedoms))
!                stop
!                if(energy < -2000) stop
            endif
            ! that last thing is a normal quench
            call run_lbfgs(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset,np)
!            print *, "real", itdone, energy

            if(COLDFUSION) then
                !print *,"Coldfusion"
            else
                exit
            endif
        end do
        if(COLDFUSION) then
            print *,"Coldfusion could not be removed"
        endif

        ATOMRIGIDCOORDT = tmp
        xrigidcoords(1:degfreedoms) = xcoords(1:degfreedoms)
    end subroutine

    subroutine run_lbfgs(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset, np)
        use commons, only : natoms
        use genrigid
        implicit none

        integer :: n,m,itmax,itdone,np
        double precision :: xcoords(3*natoms), eps,energy
        logical :: diagco, reset,  mflag
        call mymylbfgs(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset,np)
    end subroutine

end module
