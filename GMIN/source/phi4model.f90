SUBROUTINE Phi4Model(XC, DX, ENERGY, GTEST)

    ! *****************************************************************************************


    !USE MODHESS
    USE COMMONS, ONLY: NATOMS
    !USE KEY, ONLY: FROZEN, NFREEZE

    IMPLICIT NONE
    DOUBLE PRECISION :: XC(3*NATOMS), DX(3*NATOMS) !Cartesian coords, derivatives
    DOUBLE PRECISION :: ENERGY, ENERGYNN, ENERGYDUMMY1, ENERGYDUMMY2, DUMMYDX
    LOGICAL          :: GTEST !, STEST
    INTEGER          :: i,j,N
    REAL, PARAMETER  :: LAMBDA = 0.6 ! LAMBDA = 3/5
    REAL, PARAMETER  :: MUSQUARE = 1.00
    REAL, PARAMETER  :: JPARAM = 0.001


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D Nearest-Neighbour Phi^4 model

!    ENERGY = 0.00    
!    DO i = 1, NATOMS
!      ENERGY = ENERGY + (LAMBDA/24.0)* XC(i)**4 - (MUSQUARE/2.0)*XC(i)**2

!       ENERGYNN = 0.00
!         IF (i==1) THEN
!            ENERGYNN = ENERGYNN + (JPARAM/4.0)*(XC(i) - XC(NATOMS))**2 + (JPARAM/4.0)*(XC(i) - XC(i+1))**2
!             ELSE IF (i==NATOMS) THEN
!                ENERGYNN = ENERGYNN + (JPARAM/4.0)*(XC(i) - XC(i-1))**2 + (JPARAM/4.0)*(XC(i) - XC(1))**2
!             ELSE
!          ENERGYNN = ENERGYNN + (JPARAM/4.0)*(XC(i) - XC(i-1))**2 + (JPARAM/4.0)*(XC(i) - XC(i+1))**2
!          END IF
!       ENERGY = ENERGY + ENERGYNN
!    ENDDO
       


!   IF (.NOT.GTEST) RETURN
!       DO i = 1, NATOMS
!          DX(i) = (LAMBDA/6.0)* XC(i)**3 + (4.0 * JPARAM - MUSQUARE)*XC(i)
!          DUMMYDX = 0.00
!          IF (i==1) THEN
!             DUMMYDX = DUMMYDX + XC(NATOMS) + XC(i+1)
!             ELSE IF (i==NATOMS) THEN
!                DUMMYDX = DUMMYDX + XC(i-1) + XC(1)
!                ELSE
!                   DUMMYDX = DUMMYDX + XC(i-1) + XC(i+1)
!                END IF
!                DUMMYDX = JPARAM*DUMMYDX

!                DX(i) = DX(i) - DUMMYDX
!             ENDDO
    
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mean-field phi^4 model

N=3*NATOMS

ENERGY=0.0

DO i = 1, N
   ENERGY = ENERGY + (-0.5*XC(i)**2 + 0.25*XC(i)**4)
ENDDO
ENERGYDUMMY1 = 0.0
DO i = 1, N
   ENERGYDUMMY1 = ENERGYDUMMY1 + XC(i)
ENDDO
ENERGY = ENERGY - (JPARAM/(2.0*N)) *(ENERGYDUMMY1**2)

IF (.NOT.GTEST) RETURN

DUMMYDX = 0.0
DO i = 1, N
DUMMYDX = DUMMYDX + (JPARAM/N)*XC(i)
ENDDO


DO i = 1, N
   DX(i) = - XC(i) + XC(i)**3 - DUMMYDX
ENDDO



!    IF (STEST) THEN
!        PRINT*, 'Warning: There is no analytical hessian implemented for the phi4 model yet.'
!        CALL MAKENUMHESS(XC,NATOMS)
!    ENDIF

!call writespec_xyz(17,XC)

END SUBROUTINE Phi4Model

subroutine writespec_xyz(out_unit,COORDS)
    USE COMMONS, ONLY: NATOMS, ZSYM
    implicit none
    DOUBLE PRECISION, intent(in):: coords(3*NATOMS)
    INTEGER :: IAT
    INTEGER, INTENT(IN):: out_unit

    ! **********************************************************************
    !  open (out_unit,file="Nimetpath.xyz",ACTION="WRITE")

    write(out_unit,*) NATOMS
    write(out_unit,*)
    DO  IAT = 1,NATOMS
        write(out_unit,'(a2,3f12.7)') ZSYM(IAT),coords(3*IAT-2),&
        coords(3*IAT-1),coords(3*IAT)
    ENDDO
!  call flush(out_unit)
!  close(out_unit)
end subroutine writespec_xyz
