MODULE INERTIA_MOD
    USE PREC
    IMPLICIT NONE

CONTAINS

SUBROUTINE MOM_INERTIA(INER_TENS, COORDS)
    USE COMMONS, ONLY: NATOMS, ATMASS, DEBUG
    IMPLICIT NONE
! Arguments
    REAL(REAL64), INTENT(OUT)   :: INER_TENS(3, 3)
    REAL(REAL64), INTENT(IN)    :: COORDS(:)
! Variables
    INTEGER                     :: I
    REAL(REAL64)                :: M
    REAL(REAL64), ALLOCATABLE   :: MASSES(:)
    REAL(REAL64)                :: X, Y, Z
    REAL(REAL64)                :: COM(3)
    REAL(REAL64)                :: FIRST_MOM(3)

! Get the centre of mass
    CALL CENTRE_OF_MASS(COORDS, COM)

! Now assign masses
    IF (.NOT. ALLOCATED(MASSES)) ALLOCATE(MASSES(NATOMS))
    IF (ALLOCATED(ATMASS)) THEN
        MASSES(:) = ATMASS(:)
    ELSE
        MASSES(:) = 1.0D0
    END IF

! Loop through the atoms, calculating contributions to the moment of inertia tensor.
    INER_TENS = 0.0D0
    DO I = 1, NATOMS
        M = MASSES(I)
        X = COORDS(3*I - 2) - COM(1)
        Y = COORDS(3*I - 1) - COM(2)
        Z = COORDS(3*I)     - COM(3)
        ! Diagonal terms
        INER_TENS(1, 1) = M*(Y*Y + Z*Z) + INER_TENS(1, 1)
        INER_TENS(2, 2) = M*(Z*Z + X*X) + INER_TENS(2, 2)
        INER_TENS(3, 3) = M*(X*X + Y*Y) + INER_TENS(3, 3)
        ! Off-diagonal terms
        INER_TENS(1, 2) = -M*X*Y + INER_TENS(1, 2)
        INER_TENS(2, 3) = -M*Y*Z + INER_TENS(2, 3)
        INER_TENS(3, 1) = -M*Z*X + INER_TENS(3, 1)
    END DO

! Symmetrise the tensor.
    INER_TENS(2, 1) = INER_TENS(1, 2)
    INER_TENS(3, 2) = INER_TENS(2, 3)
    INER_TENS(1, 3) = INER_TENS(3, 1)
 
END SUBROUTINE MOM_INERTIA

SUBROUTINE CENTRE_OF_MASS(COORDS, COM)
    USE COMMONS, ONLY: NATOMS, ATMASS, DEBUG
    IMPLICIT NONE
! Arguments
    REAL(REAL64), INTENT(IN)    :: COORDS(:)
    REAL(REAL64), INTENT(OUT)   :: COM(3)
! Variables
    INTEGER                     :: I
    REAL(REAL64)                :: M                 
    REAL(REAL64), ALLOCATABLE   :: MASSES(:)
    REAL(REAL64)                :: X, Y, Z

! Assign masses
    IF (.NOT. ALLOCATED(MASSES)) ALLOCATE(MASSES(NATOMS))
    IF (ALLOCATED(ATMASS)) THEN
        MASSES(:) = ATMASS(:)
    ELSE
        MASSES(:) = 1.0D0
    END IF

! Initialise centre of mass to zero and use as an accumulator.
    COM(:) = 0.0D0
    DO I = 1, NATOMS
        M = MASSES(I)
        X = COORDS(3*I - 2)
        Y = COORDS(3*I - 1)
        Z = COORDS(3*I)
        COM(1) = COM(1) + M*X
        COM(2) = COM(2) + M*Y
        COM(3) = COM(3) + M*Z
    END DO

! Divide by total mass, since we have used COM as an accumulator.
    COM(:) = COM(:) / SUM(MASSES)

END SUBROUTINE CENTRE_OF_MASS

SUBROUTINE ORIENT(COORDS_IN, COORDS_OUT)
    USE COMMONS, ONLY: NATOMS
    IMPLICIT NONE
! Arguments
    REAL(REAL64), INTENT(IN)    :: COORDS_IN(:)
    REAL(REAL64), INTENT(OUT)   :: COORDS_OUT(:)
! Variables
    REAL(REAL64), ALLOCATABLE   :: COORDS(:)
    REAL(REAL64)                :: COM_COORDS(3)
    REAL(REAL64)                :: INERTIA_TENS(3, 3)
    REAL(REAL64)                :: PRINCIPAL_AXES(3, 3)
    REAL(REAL64)                :: I_XYZ(3)
    REAL(REAL64)                :: MAG
    REAL(REAL64)                :: ROT_MAT(3, 3)
    REAL(REAL64)                :: X, Y, Z
    INTEGER                     :: I
! LAPACK variables
    REAL(REAL64)                :: WORK(1000)
    INTEGER                     :: LWORK = 1000
    INTEGER                     :: INFO
    INTEGER                     :: PIV_IDXS(3)

! Transform the coordinates to the centre of mass frame
    CALL CENTRE_OF_MASS(COORDS_IN, COM_COORDS)
    IF (.NOT. ALLOCATED(COORDS)) THEN
        ALLOCATE(COORDS(SIZE(COORDS_IN)))
    END IF
    DO I = 1, NATOMS
        COORDS(3*I - 2) = COORDS_IN(3*I - 2) - COM_COORDS(1)
        COORDS(3*I - 1) = COORDS_IN(3*I - 1) - COM_COORDS(2)
        COORDS(3*I)     = COORDS_IN(3*I)     - COM_COORDS(3)
    END DO

! Calculate the moment of inertia tensor
    CALL MOM_INERTIA(INERTIA_TENS, COORDS)
! Diagonalise the moment of inertia tensor and find the principal axes
    PRINCIPAL_AXES(:, :) = INERTIA_TENS(:, :)
! DSYEV changes the value of the input matrix on output, hence copying in advance.
    CALL DSYEV('V', &
               'L', &
               3, &
               PRINCIPAL_AXES, &
               3, &
               I_XYZ, &
               WORK, &
               LWORK, &
               INFO)

! Normalise the principal axes (just in case)
    DO I = 1, 3
        MAG = PRINCIPAL_AXES(1, I)**2 + &
              PRINCIPAL_AXES(2, I)**2 + &
              PRINCIPAL_AXES(3, I)**2
        PRINCIPAL_AXES(1, I) = PRINCIPAL_AXES(1, I) / MAG
        PRINCIPAL_AXES(2, I) = PRINCIPAL_AXES(2, I) / MAG
        PRINCIPAL_AXES(3, I) = PRINCIPAL_AXES(3, I) / MAG
    END DO

! We wish to find a transformation that orients the coordinates such that:
! 
! I = R P
!
! where P is the matrix containing the principal axes and I is the identity matrix.
!
! R is thus P^-1
    ROT_MAT = PRINCIPAL_AXES

! Calculate the LU factorisation using DGETRF
    CALL DGETRF(3, 3, ROT_MAT, 3, PIV_IDXS, INFO)       

! Use this to calculate the inverse.
    CALL DGETRI(3, ROT_MAT, 3, PIV_IDXS, WORK, 3, INFO)

! Now apply our rotation matrix to the coordinates
    DO I = 1, NATOMS
        X = COORDS(3*I - 2) 
        Y = COORDS(3*I - 1) 
        Z = COORDS(3*I) 
        COORDS(3*I - 2) = ROT_MAT(1, 1)*X + ROT_MAT(1, 2)*Y + ROT_MAT(1, 3)*Z
        COORDS(3*I - 1) = ROT_MAT(2, 1)*X + ROT_MAT(2, 2)*Y + ROT_MAT(2, 3)*Z
        COORDS(3*I)     = ROT_MAT(3, 1)*X + ROT_MAT(3, 2)*Y + ROT_MAT(3, 3)*Z
    END DO

! Assign to COORDS_OUT
    COORDS_OUT = COORDS

END SUBROUTINE ORIENT

END MODULE INERTIA_MOD
