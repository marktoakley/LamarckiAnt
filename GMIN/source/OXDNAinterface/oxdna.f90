! this wrapper routine is here since I'm not sure how fortran logical translates to C
! called by gmin to calculate potential
subroutine userpot_potential(dof,X,GRAD,EREAL,GRADT)
    integer, intent(in) :: dof                 ! number of degrees of freedom
    double precision, intent(in) :: X(dof)     ! current coordinates
    double precision, intent(out) :: GRAD(dof) ! gradient
    double precision, intent(out) :: EREAL     ! energy
    logical, intent(in) :: gradt              ! is the gradient needed?

    call calc_potential(x, grad, ereal)
end subroutine

subroutine userpot_dump_lowest()
    USE QMODULE
    USE COMMONS, ONLY : NSAVE,NATOMS

    IMPLICIT NONE
    INTEGER I
    CHARACTER(LEN=10) :: CNUM

    DO I=1,NSAVE
        WRITE (CNUM,'(I3.3)') I
        print *,I,"lowest"//trim(cnum)//".dat"

        CALL USERPOT_DUMP_CONFIGURATION(&
            "lowest"//trim(cnum)//".dat",&
            QMINP(I,1:3*NATOMS))
    END DO
end subroutine
