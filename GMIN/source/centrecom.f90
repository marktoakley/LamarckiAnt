SUBROUTINE CENTRECOM(X)
    USE COMMONS, ONLY: NATOMS, ATMASS, DEBUG, MYUNIT
    USE PREC
    IMPLICIT NONE
! Arguments
    REAL(REAL64), INTENT(INOUT) :: X(3 * NATOMS)
! Variables
    INTEGER(INT32)              :: I
    REAL(REAL64)                :: COM(3), TOTMASS

    TOTMASS = SUM(ATMASS)

    DO I = 1, 3
        COM(I) = SUM(ATMASS(:) * X(I::3)) / TOTMASS
        X(I::3) = X(I::3) - COM(I)
    END DO

    IF (DEBUG) THEN 
        WRITE(MYUNIT,'(A,3F15.10)') 'Centre of mass reset to the origin from ', &
                                    COM(:)
    ENDIF
!    PRINT*,'final coordinates in centrecom:'
!    WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)

END SUBROUTINE CENTRECOM
