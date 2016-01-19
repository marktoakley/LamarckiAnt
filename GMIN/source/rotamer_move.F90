MODULE ROTAMER
   IMPLICIT NONE
   LOGICAL                 :: ROTAMER_MOVET
   CHARACTER(LEN=100)      :: ROTAMER_SCRIPT

CONTAINS

SUBROUTINE ROTAMER_INIT()
   USE PORFUNCS
   IMPLICIT NONE

! If the Makefile/CMake have defined _SVN_ROOT_, we can use that to show where the
! rotamer python script should be. If not, rely on the user-defined location.
#ifdef _SVN_ROOT_
   CALL SYSTEM(_SVN_ROOT_ // '/SCRIPTS/AMBER/rotamer/rotamer.py' // ' init coords.prmtop')
#else
   CALL SYSTEM(TRIM(ADJUSTL(ROTAMER_SCRIPT)) // ' init coords.prmtop')
#endif 

END SUBROUTINE ROTAMER_INIT

SUBROUTINE ROTAMER_STEP(COORDS)
! Main driver routine for rotamer moves.
!
! COORDS is the coordinates to be changed by the step. Most of the work is done by the rotamer.py
! script, which is also used for initialisation. Settings are controlled by the file
! rotamer_settings in the working directory.
   USE PORFUNCS
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, INTENT(INOUT)  :: COORDS(:)

   CALL DUMP_COORDS(COORDS)
#ifdef _SVN_ROOT_
   CALL SYSTEM(_SVN_ROOT_ // '/SCRIPTS/AMBER/rotamer/rotamer.py' // ' step coords.prmtop')
#else
   CALL SYSTEM(TRIM(ADJUSTL(ROTAMER_SCRIPT)) // ' step coords.prmtop')
#endif 
   CALL READ_NEW_COORDS(COORDS)

END SUBROUTINE ROTAMER_STEP

SUBROUTINE DUMP_COORDS(COORDS)
! Dump coordinates in Amber restart format to .coords_before_rotamer.rst for rotamer moves.
   USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_WRITE_RESTART
   USE COMMONS, ONLY: AMBERT, AMBER12T, NATOMS
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, INTENT(IN)  :: COORDS(:)

   IF (AMBERT) THEN
   ! If not using AMBER 12, then use DUMPRST from hbondmatrix.f90
      CALL DUMPRST(COORDS(:), NATOMS, '.coords_before_rotamer.rst')
   ELSE IF (AMBER12T) THEN
      CALL AMBER12_WRITE_RESTART(COORDS(:), &
                                 '.coords_before_rotamer.rst', &
                                 LEN('.coords_before_rotamer.rst'))
   ELSE
      WRITE(*, *) 'Rotamer moves can only currently be used with the Amber potentials.'
      STOP
   END IF

END SUBROUTINE DUMP_COORDS

SUBROUTINE READ_NEW_COORDS(COORDS)
! Read the new coordinates (after rotamer moves).
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, INTENT(OUT) :: COORDS(:)
! Variables
   INTEGER                       :: UNIT_NUM
   CHARACTER(LEN=100)            :: FILE_NAME
   INTEGER                       :: N_LINES, I
   UNIT_NUM = 5827
   FILE_NAME = '.coords_after_rotamer.xyz'
   OPEN(UNIT=UNIT_NUM, FILE=FILE_NAME, STATUS='UNKNOWN')
   READ(UNIT_NUM, '(I10)') N_LINES
   DO I = 1, N_LINES
      READ(UNIT_NUM, '(3F20.10)') COORDS(3*I-2:3*I)
   END DO
   CLOSE(UNIT_NUM)

END SUBROUTINE READ_NEW_COORDS

END MODULE ROTAMER
