MODULE SHIFT_HESS
    USE PREC
    IMPLICIT NONE

CONTAINS

SUBROUTINE SHIFT_HESS_ZEROS(INPUT_COORDS, SHIFTS)
    USE COMMONS, ONLY: ATMASS
    USE INERTIA_MOD, ONLY: MOM_INERTIA, CENTRE_OF_MASS
    USE MODHESS, ONLY: HESS, MASS_WEIGHTED, MASSWT
    IMPLICIT NONE
! Arguments
    REAL(REAL64), DIMENSION(:), INTENT(IN)  :: INPUT_COORDS
    REAL(REAL64), DIMENSION(6), INTENT(IN)  :: SHIFTS
! Variables
    INTEGER                                 :: NUM_COORDS
    REAL(REAL64), DIMENSION(:, :), &
    ALLOCATABLE                             :: ZERO_VECTORS
    INTEGER                                 :: N, I
    INTEGER                                 :: X, Y, Z
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: MASSES
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: COORDS
    REAL(REAL64), DIMENSION(3)              :: COM
    REAL(REAL64), DIMENSION(3, 3)           :: INERTIA_TENS
    REAL(REAL64), DIMENSION(3)              :: I_123
    REAL(REAL64), DIMENSION(3, 3)           :: PRINCIPAL_AXES
    REAL(REAL64)                            :: M
! DSYEV variables
    REAL(REAL64), DIMENSION(1000)           :: WORK
    INTEGER, PARAMETER                      :: LWORK = 1000
    INTEGER                                 :: INFO
! Parameters
    REAL(REAL64), PARAMETER                 :: MAX_INERTIA_ZEROS = 1.0D-3

! Check that the number of coordinates and size of hessian correspond
    NUM_COORDS = SIZE(INPUT_COORDS)
    IF (NUM_COORDS .NE. SIZE(HESS, 1)) THEN
        STOP 'Coords/hessian size mismatch in subroutine SHIFT_HESS_ZEROS.'
    END IF

! Set masses, if MASS_WEIGHTED, otherwise set to 1.0D0.
    IF (.NOT. ALLOCATED(MASSES)) ALLOCATE(MASSES(NUM_COORDS/3))
    MASSES(:) = ATMASS(:)

! Transform INPUT_COORDS into the centre of mass frame, storing the result in COORDS
    IF (.NOT. ALLOCATED(COORDS)) ALLOCATE(COORDS(NUM_COORDS))
    CALL CENTRE_OF_MASS(INPUT_COORDS, COM)
    DO I = 1, 3
        COORDS(I:NUM_COORDS:3) = INPUT_COORDS(I:NUM_COORDS:3) - COM(I)
    END DO

! Calculate the moment of inertia tensor
    CALL MOM_INERTIA(INERTIA_TENS, COORDS)
! Diagonalise the moment of inertia tensor and find the principal axes
    PRINCIPAL_AXES(:, :) = INERTIA_TENS(:, :)
! DSYEV changes the value of the input matrix on output, hence copying in advance.
! This gives the eigenvalues are in ascending order in I_123 and their corresponding principal axes
! are in the columns of PRINCIPAL_AXES. i.e. PRINCIPAL_AXES(:, N) is the eigenvector corresponding
! to I_123(N).
    CALL DSYEV('V', &
               'L', &
               3, &
               PRINCIPAL_AXES, &
               3, &
               I_123, &
               WORK, &
               LWORK, &
               INFO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For debugging
! 
!    PRINT *, 'Calculated principal axes'
!    PRINT '(3F15.7)', PRINCIPAL_AXES(1, :)
!    PRINT '(3F15.7)', I_123(1)
!    PRINT '(3F15.7)', PRINCIPAL_AXES(2, :)
!    PRINT '(3F15.7)', I_123(2)
!    PRINT '(3F15.7)', PRINCIPAL_AXES(3, :)
!    PRINT '(3F15.7)', I_123(3)
!    PRINT *, 'End of calculated principal axes'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create zero eigenvectors
    IF (.NOT. ALLOCATED(ZERO_VECTORS)) ALLOCATE(ZERO_VECTORS(6, NUM_COORDS))
    ZERO_VECTORS(:, :) = 0.0D0
! However, we must ensure that the system is oriented correctly, prior to calculating the hessian, otherwise
! these are not the rotational eigenvectors.
    DO N = 1, NUM_COORDS / 3
        X = 3 * N - 2
        Y = 3 * N - 1
        Z = 3 * N
        M = MASSES(N)
        DO I = 1, 3
        ! translations along principal axes
            ZERO_VECTORS(I, X)     = SQRT(MASSES(N)) * PRINCIPAL_AXES(1, I)
            ZERO_VECTORS(I, Y)     = SQRT(MASSES(N)) * PRINCIPAL_AXES(2, I)
            ZERO_VECTORS(I, Z)     = SQRT(MASSES(N)) * PRINCIPAL_AXES(3, I)
        ! rotation about principal axes
            ZERO_VECTORS(I + 3, X) = SQRT(M) * (PRINCIPAL_AXES(3, I) * COORDS(Y) - PRINCIPAL_AXES(2, I) * COORDS(Z))
            ZERO_VECTORS(I + 3, Y) = SQRT(M) * (PRINCIPAL_AXES(1, I) * COORDS(Z) - PRINCIPAL_AXES(3, I) * COORDS(X))
            ZERO_VECTORS(I + 3, Z) = SQRT(M) * (PRINCIPAL_AXES(2, I) * COORDS(X) - PRINCIPAL_AXES(1, I) * COORDS(Y))
        END DO
    END DO

! Normalise the zero vectors.
    DO I = 1, 6
        IF(ALL(ZERO_VECTORS(I, :) == 0.0D0)) CYCLE
        ZERO_VECTORS(I, :) = ZERO_VECTORS(I, :) / SQRT(SUM(ZERO_VECTORS(I, :)**2))
    END DO

! Now apply the shifting to the hessian.
    DO I = 1, 6
        CALL SHIFT_EIGENVALUE(HESS, SHIFTS(I), ZERO_VECTORS(I, :))
    END DO 

END SUBROUTINE SHIFT_HESS_ZEROS        

SUBROUTINE SHIFT_EIGENVALUE(HESSIAN, SHIFT_VALUE, EIGENVECTOR)
! Shifts the eigenvalue of HESSIAN corresponding to EIGENVECTOR by
! +SHIFT_VALUE.
    IMPLICIT NONE
! Arguments
    REAL(REAL64), DIMENSION(:, :), INTENT(INOUT)    :: HESSIAN
    REAL(REAL64), INTENT(IN)                        :: SHIFT_VALUE
    REAL(REAL64), DIMENSION(:), INTENT(IN)          :: EIGENVECTOR
! Variables
    INTEGER                                         :: A, B
    INTEGER                                         :: NUM_COORDS

! Check the number of coordinates.
    NUM_COORDS = SIZE(HESSIAN, 1)
    IF (NUM_COORDS .NE. SIZE(EIGENVECTOR)) THEN
        STOP 'Invalid eigenvector/hessian combination in subroutine SHIFT_EIGENVALUE.'
    END IF

! Apply formula Hout_ab = H_ab + L e_a e_b
! See Energy Landscapes: Eq 6.3
    DO A = 1, NUM_COORDS
        DO B = A, NUM_COORDS
            HESSIAN(A, B) = HESSIAN(A, B) + SHIFT_VALUE * EIGENVECTOR(A) * EIGENVECTOR(B)
            HESSIAN(B, A) = HESSIAN(A, B)
        END DO
    END DO

END SUBROUTINE SHIFT_EIGENVALUE

END MODULE SHIFT_HESS
