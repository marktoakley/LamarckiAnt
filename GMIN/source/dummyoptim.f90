! dummy routines that GMIN is not calling, but the compilation for AMBER interface needs them

!SUBROUTINE CHIRALH_ALIGN(COORDSB, COORDSA, i1, DEBUG)
!use commons, only : MYUNIT, NATOMS
!implicit none
!DOUBLE PRECISION        :: COORDSB(3*NATOMS), COORDSA(3*NATOMS)
!LOGICAL                 :: DEBUG
!INTEGER                 :: i1
!
!WRITE(MYUNIT,*) ' Should not be in this dummy routine!!! (dummyoptim.f90)'
!STOP
!END SUBROUTINE CHIRALH_ALIGN

SUBROUTINE INTH_ALIGN(COORDSB,COORDSA,centre)
implicit none
DOUBLE PRECISION        :: COORDSB(:), COORDSA(:)
INTEGER                 :: centre 

STOP 'Should not be in this dummy routine!!! (dummyoptim.f90)'
END SUBROUTINE INTH_ALIGN
