MODULE TWIST_MOD

USE PREC, ONLY: INT32, REAL64

! TWIST_K is force constants, TWIST_THETA0 is equilibrium dihedral angles.
REAL(REAL64), ALLOCATABLE     :: TWIST_K(:), TWIST_THETA0(:)
! TWIST_ATOMS contains indices of atoms involved in the dihedral.
INTEGER(INT32), ALLOCATABLE   :: TWIST_ATOMS(:, :)
LOGICAL                       :: TWISTT
INTEGER(INT32)                :: NTWISTGROUPS
REAL(REAL64)                  :: PI = 4.0_REAL64 * ATAN(1.0_REAL64)

CONTAINS

SUBROUTINE TWIST(COORDS, NATOMS, GRAD, ENERGY, GTEST)
! Update the energy and gradient with contributions from a twisting potential
! which adds torsional terms. These have the functional form:
!
!
!
   USE PREC, ONLY: INT32, REAL64
   USE COMMONS, ONLY: MYUNIT
   USE DIHEDRAL_MOD, ONLY: DIHEDRAL_WITH_GRAD
   IMPLICIT NONE
! Arguments
   INTEGER(INT32), INTENT(IN)    :: NATOMS
   REAL(REAL64), INTENT(IN)      :: COORDS(3*NATOMS)
   LOGICAL, INTENT(IN)           :: GTEST
   REAL(REAL64), INTENT(INOUT)   :: GRAD(:), ENERGY
! Variables
   REAL(REAL64), ALLOCATABLE     :: ANGLES(:)
   REAL(REAL64), ALLOCATABLE     :: DIHEDRAL_COORDS(:, :)
   INTEGER(INT32)                :: NUM_DIHEDRALS
   INTEGER(INT32)                :: I, J, IDX
   REAL(REAL64), ALLOCATABLE     :: TWIST_ENERGIES(:)
   REAL(REAL64), ALLOCATABLE     :: DE_DFS(:)
   REAL(REAL64), ALLOCATABLE     :: DF_DRS(:, :)
   REAL(REAL64), ALLOCATABLE     :: TWIST_GRADS(:, :)

! First check that there are the same number of force constants, equilibrium angles
! and defined dihedral atom sets.
   NUM_DIHEDRALS = SIZE(TWIST_K)
   IF (.NOT. (       (NUM_DIHEDRALS == SIZE(TWIST_THETA0)) &
               .AND. (NUM_DIHEDRALS == SIZE(TWIST_ATOMS, 2)))) THEN
      WRITE(MYUNIT, *) 'Check the twist input file. The numbers of force constants, equilibrium angles' // &
                       ' and defined dihedrals do not match.'
      STOP
   END IF

! Calculate dihedral angles for each set of atoms.
   ALLOCATE(ANGLES(NUM_DIHEDRALS))
   ALLOCATE(DIHEDRAL_COORDS(12, NUM_DIHEDRALS))
   ALLOCATE(DF_DRS(12, NUM_DIHEDRALS))
   DO I = 1, NUM_DIHEDRALS
      DO J = 1, 4
         IDX = TWIST_ATOMS(J, I)
         DIHEDRAL_COORDS(3*J-2:3*J, I) = COORDS(3*IDX-2:3*IDX)
      END DO
! Calculates the dihedral angles *and* the gradients of the dihedral wrt coordinates.
      CALL DIHEDRAL_WITH_GRAD(DIHEDRAL_COORDS(1:12, I), ANGLES(I), DF_DRS(1:12, I))
! Restrict each angle to be between pi and -pi
      CALL ENFORCE_ANGLE_RANGE(ANGLES(I))
   END DO

! Calculate the contribution to the energy from each dihedral
! E_i = k_i * (theta_i - theta0_i)^2
   ALLOCATE(TWIST_ENERGIES(NUM_DIHEDRALS))
   DO I = 1, NUM_DIHEDRALS
      TWIST_ENERGIES(I) = TWIST_K(I) * (ANGLES(I) - TWIST_THETA0(I))**2
   END DO

! Calculate gradient terms using the chain rule
!
! dE   dE   df
! -- = -- . --
! dr   df   dr

! dE/df is simply k_i * cos(theta_i - theta0_i)
! df/dr was calculated earlier
   ALLOCATE(DE_DFS(NUM_DIHEDRALS))
   DO I = 1, NUM_DIHEDRALS
      DE_DFS(I) = TWIST_K(I) * COS(ANGLES(I) - TWIST_THETA0(I))
   END DO

! Each 12-coordinate term of TWIST_GRADS is the product of each DE_DFS term with
! the corresponding 12-coordinate DF_DRS term. 
   ALLOCATE(TWIST_GRADS(12, NUM_DIHEDRALS))
   DO I = 1, NUM_DIHEDRALS
      TWIST_GRADS(1:12, I) = DE_DFS(I) * DF_DRS(1:12, I)
   END DO

! Map each of these back to the corresponding parts of the gradient.
   DO I = 1, NUM_DIHEDRALS
      DO J = 1, 4
! IDX refers to the index of the atom whose gradient is currently being mapped.
         IDX = TWIST_ATOMS(J, I)
         GRAD(3*IDX-2:3*IDX) = GRAD(3*IDX-2:3*IDX) + TWIST_GRADS(3*J-2:3*J, I)
      END DO
   END DO

! Free dynamic memory
   IF (ALLOCATED(ANGLES)) DEALLOCATE(ANGLES)
   IF (ALLOCATED(DIHEDRAL_COORDS)) DEALLOCATE(DIHEDRAL_COORDS)
   IF (ALLOCATED(DF_DRS)) DEALLOCATE(DF_DRS)
   IF (ALLOCATED(TWIST_ENERGIES)) DEALLOCATE(TWIST_ENERGIES)
   IF (ALLOCATED(DE_DFS)) DEALLOCATE(DE_DFS)
   IF (ALLOCATED(TWIST_GRADS)) DEALLOCATE(TWIST_GRADS)
END SUBROUTINE TWIST

SUBROUTINE ENFORCE_ANGLE_RANGE(ANGLE)
! Forces an angle to be between pi and -pi
   IMPLICIT NONE
! Arguments
   REAL(REAL64), INTENT(INOUT)   :: ANGLE
   
   ANGLE = MOD(ANGLE, 2*PI)
   IF (ANGLE < -1.0 * PI) THEN
      ANGLE = ANGLE + 2*PI
   ELSE IF (ANGLE > PI) THEN
      ANGLE = ANGLE - 2*PI
   END IF

END SUBROUTINE ENFORCE_ANGLE_RANGE

END MODULE TWIST_MOD
