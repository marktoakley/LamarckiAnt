MODULE MOVES
! This module is for the magical new way of doing custom movesets and also fixing our
! horrible old way of doing moves with 15 different implementations of Cartesian steps.
! Please code nicely in here...GOTOs will be removed, as will unclear variable names.
! Consistent indentation is mandatory!

CONTAINS

SUBROUTINE CARTESIAN_SPHERE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random spherically symmetric displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved
 
! The VEC3 module (vec3.f90) contains helper functions for handling vectors and matricies
   USE VEC3
! The SANITY module contains sanity check functions
   USE SANITY
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   DOUBLE PRECISION                              :: DPRAND
   DOUBLE PRECISION, INTENT(IN)                  :: MAX_STEP
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SPHERE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
      XYZ(3*I-2:3*I)=XYZ(3*I-2:3*I)+VEC_RANDOM()*(DPRAND()**(1.0D0/3.0D0))*MAX_STEP
   ENDDO

END SUBROUTINE CARTESIAN_SPHERE

SUBROUTINE CARTESIAN_SIMPLE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved

! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   DOUBLE PRECISION                              :: DPRAND
   DOUBLE PRECISION                              :: RANDOMX
   DOUBLE PRECISION                              :: RANDOMY
   DOUBLE PRECISION                              :: RANDOMZ
   DOUBLE PRECISION, INTENT(IN)                  :: MAX_STEP
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SIMPLE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
! Draw a random number between -1 and 1 for each coordinate
      RANDOMX=(DPRAND()-0.5D0)*2.0D0
      RANDOMY=(DPRAND()-0.5D0)*2.0D0
      RANDOMZ=(DPRAND()-0.5D0)*2.0D0
! Displace each coordinate
      XYZ(3*I-2)=XYZ(3*I-2)+MAX_STEP*RANDOMX
      XYZ(3*I-1)=XYZ(3*I-1)+MAX_STEP*RANDOMY
      XYZ(3*I  )=XYZ(3*I  )+MAX_STEP*RANDOMZ
   ENDDO

END SUBROUTINE CARTESIAN_SIMPLE

SUBROUTINE ROTATION_ABOUT_AXIS(XYZ, VECTOR_START_XYZ, &
                               VECTOR_END_XYZ, ANGLE_DEGS, ATOM_LIST)
!
! Rotate the coordinates of the atoms in ATOM_LIST about the line from VECTOR_START_XYZ
! to VECTOR_END_XYZ through ANGLE degrees.
!
! Arguments
! ---------
!
! Required:
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! VECTOR_START_XYZ(in): start of the line about which to rotate, in Cartesian coords
! VECTOR_END_XYZ(in): end of the line about which to rotate, in Cartesian coords
! ANGLE_DEGS(in): angle through which to rotate, in degrees
!
! Optional:
! ATOM_LIST(in): list of atoms to be rotated
!
! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: XYZ
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)     :: ATOM_LIST
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)      :: VECTOR_START_XYZ
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)      :: VECTOR_END_XYZ
   DOUBLE PRECISION, INTENT(IN)                    :: ANGLE_DEGS
! Constants
! Variables
   DOUBLE PRECISION                                :: PI
   DOUBLE PRECISION                                :: DEGS_OVER_RADS
   DOUBLE PRECISION                                :: COS_THETA, SIN_THETA
   DOUBLE PRECISION                                :: A, B, C
   DOUBLE PRECISION                                :: D, E, F
   DOUBLE PRECISION                                :: U, V, W
   DOUBLE PRECISION                                :: X, Y, Z
   DOUBLE PRECISION                                :: VECTOR_MAG
   INTEGER                                         :: NUM_ATOMS
   INTEGER                                         :: I
   LOGICAL, ALLOCATABLE, DIMENSION(:)              :: ATOM_MASK
   LOGICAL                                         :: TEST
! Function declarations
   DOUBLE PRECISION                                :: DNRM2

! Define PI using ATAN and the conversion factor between degrees and radians.
! x degrees = x * DEGS_OVER_RADS radians
   PI = 4.0D0 * ATAN(1.0D0)
   DEGS_OVER_RADS = PI / 180.0D0

! Calculate cos and sin of the angle.
   SIN_THETA = SIN(ANGLE_DEGS * DEGS_OVER_RADS)
   COS_THETA = COS(ANGLE_DEGS * DEGS_OVER_RADS)

! Assign variables according to the functional form described on:
! http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html

   A = VECTOR_START_XYZ(1)
   B = VECTOR_START_XYZ(2)
   C = VECTOR_START_XYZ(3)
   
   D = VECTOR_END_XYZ(1)
   E = VECTOR_END_XYZ(2)
   F = VECTOR_END_XYZ(3)

! DNRM2(INTEGER N, REAL X(N), INTEGER INCX)
   VECTOR_MAG = DNRM2(3, VECTOR_END_XYZ - VECTOR_START_XYZ, 1)

   U = (D - A) / VECTOR_MAG
   V = (E - B) / VECTOR_MAG
   W = (F - C) / VECTOR_MAG

! Work out how many atoms the coordinates passed describe
   NUM_ATOMS = SIZE(XYZ) / 3 
! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to ROTATION_ABOUT_AXIS'
   ENDIF

! Convert ATOM_LIST (a list of atoms to be rotated) into ATOM_MASK, which can
! applied in a WHERE loop.
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Loop through and apply the formula if the atom described is in ATOM_LIST
! N.B. I've rearranged the formula so that the variables cycle, since I think it
! makes checking a bit easier.
   DO I = 1, NUM_ATOMS
      IF (.NOT. ATOM_MASK(I)) CYCLE
      X = XYZ(3 * I - 2)
      Y = XYZ(3 * I - 1)
      Z = XYZ(3 * I    )
      XYZ(3 * I - 2) = (A * (V**2 + W**2) - U * (B*V + C*W - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       X * COS_THETA + &
                       (-C*V + B*W + V*Z - W*Y) * SIN_THETA
      XYZ(3 * I - 1) = (B * (W**2 + U**2) - V * (C*W + A*U - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Y * COS_THETA + &
                       (-A*W + C*U + W*X - U*Z) * SIN_THETA
      XYZ(3 * I    ) = (C * (U**2 + V**2) - W * (A*U + B*V - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Z * COS_THETA + &
                       (-B*U + A*V + U*Y - V*X) * SIN_THETA
   END DO

END SUBROUTINE ROTATION_ABOUT_AXIS

END MODULE MOVES
