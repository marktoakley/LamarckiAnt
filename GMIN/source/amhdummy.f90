! Dummies for compiling GMIN without AMH

! Used by keywords.o
SUBROUTINE WALESAMH_INITIAL()
    IMPLICIT NONE

    STOP 'This is an AMH dummy routine. You should not be here.'
END SUBROUTINE WALESAMH_INITIAL

! Used by finalio.o
SUBROUTINE NUM_TO_CHAR(COUNT_INT, COUNT_CHAR)
    IMPLICIT NONE
    INTEGER             :: COUNT_INT
    CHARACTER (LEN=1)   :: COUNT_CHAR

    STOP 'This is an AMH dummy routine. You should not be here.'
END SUBROUTINE NUM_TO_CHAR

! Used by potential.o
SUBROUTINE WALESAMH_INTERFACE(DUMMY_COORDS, DUMMY_GRAD, DUMMY_ENERGY)
    IMPLICIT NONE
    DOUBLE PRECISION    :: DUMMY_COORDS(:)
    DOUBLE PRECISION    :: DUMMY_GRAD(:)
    DOUBLE PRECISION    :: DUMMY_ENERGY

    STOP 'This is an AMH dummy routine. You should not be here.'
END SUBROUTINE WALESAMH_INTERFACE

MODULE AMH_INTERFACES
CONTAINS
! Used by mc.o
    SUBROUTINE E_WRITE(AVEP, T, NUMPRO, I_TEMP)
        IMPLICIT NONE
        DOUBLE PRECISION :: AVEP(:, :, :)
        DOUBLE PRECISION :: T
        INTEGER          :: NUMPRO
        INTEGER          :: I_TEMP

        STOP 'This is an AMH dummy routine. You should not be here.'
    END SUBROUTINE E_WRITE
END MODULE AMH_INTERFACES

MODULE AMHGLOBALS
! Used by finalio.o and mc.o
      INTEGER           :: NMRES
      INTEGER           :: IRES(1)
! Used by finalio.o
      INTEGER           :: OMOVI
      INTEGER           :: NUMPRO
      DOUBLE PRECISION  :: AVEP(1,1,1)
END MODULE AMHGLOBALS
