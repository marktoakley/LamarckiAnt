module opepdata
    use defs
    use md_initialise
    use calcforces
    type (t_conformations), save :: conf
end module opepdata

subroutine userpot_init
    use opepdata

    call definitions()
    call initialise(conf)
end subroutine

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
    use md_initialise
    use calcforces
    integer, intent(out) :: num_atoms
    num_atoms = natoms
end subroutine

! copy all coordinates to gmin coords(:,1) array (in commons)
! TODO: what about MPI?
! for very simple standard cases this is just initial configuration,
! but can be also more complicated, e.g. initalizing generalized rigid
! body framework
subroutine userpot_initialize_gmin(dof, xvec)
    use opepdata
    integer, intent(in) :: dof
    double precision, intent(out) :: xvec(dof)
    xvec = pos
end subroutine

! called by gmin to calculate potential
subroutine userpot_potential(dof,xvec,GRAD,EREAL,GRADT)
    use md_defs
    integer, intent(in) :: dof                 ! number of degrees of freedom
    double precision, intent(in) :: xvec(dof)     ! current coordinates
    double precision, intent(out) :: GRAD(dof) ! gradient
    double precision, intent(out) :: EREAL     ! energy
    logical, intent(in) :: gradt              ! is the gradient needed?
    call calcforce(1.0d0, xvec, GRAD, EREAL)
    GRAD = -GRAD 
end subroutine

subroutine userpot_dump_configuration(filename, coords)
    use md_initialise
    use calcforces
    
end subroutine


subroutine userpot_dump_lowest() ! could be done more elegantly...
     USE QMODULE , ONLY : QMINP
     USE COMMONS , ONLY : NATOMS,COORDS,NSAVE
     IMPLICIT NONE
     character(len=100) :: line,name
     CHARACTER(LEN=10) :: fmt
     character(len=7), DIMENSION(3000) :: text2
     character(len=5), DIMENSION(3000) :: text3
     character(len=7), DIMENSION(3000) :: text4
     integer, DIMENSION(3000)          :: numres,Id_atom
     double precision, DIMENSION(3000) :: X
     integer l,i,J1
      
     OPEN(UNIT=14,FILE="conf_initiale.pdb",status="unknown")
     l = 1
     do while (l < natoms + 1)
       read(14,'(a)') line
       if(line(1:4) .eq. 'ATOM') then
         read(line,5000) text2(l),Id_atom(l),text3(l),text4(l),numres(l)
         l = l + 1
       endif
     enddo
 5000 format(a,I4,a,a,I3) 
     close(14)

    DO J1=1,NSAVE
       WRITE (fmt,'(I3.3)') J1
       name="lowest"//trim(fmt)//".pdb"
       open(unit=41,file=name,status='replace',position='append')
       write(41,'(a)') "MODEL"
     DO I=1,NATOMS
       write(41,'(a,I4,a,a,I3,f12.3,2F8.3)') text2(i),Id_atom(i), text3(i),text4(i),numres(i),QMINP(J1,3*I-2),QMINP(J1,3*I-1),QMINP(J1,3*I)
     END DO
    write(41,'(a)') "ENDMDL"
    END DO
end subroutine
