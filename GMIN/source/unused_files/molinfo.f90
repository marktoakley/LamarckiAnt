MODULE MOLINFO
    IMPLICIT NONE

    TYPE ATOM
        INTEGER            :: ATOM_ID    ! Atom number
        CHARACTER (LEN=10) :: ATOM_NAME  ! Atom name
        INTEGER            :: RES_ID     ! Residue number
        CHARACTER (LEN=10) :: RES_NAME   ! Residue name
        INTEGER            :: CHAIN_ID   ! Chain id
        CHARACTER (LEN=10) :: SEG_NAME   ! Segment name
        CHARACTER (LEN=5)  :: ELEMENT    ! Element name
        DOUBLE PRECISION   :: OCCUPANCY  ! Site occupancy
        DOUBLE PRECISION   :: B_FACTOR   ! Crystallographic temp.
        DOUBLE PRECISION   :: CHARGE     ! Atomic charge
        DOUBLE PRECISION   :: WEIGHT     ! CHARMM weight
    END TYPE ATOM

    TYPE(ATOM), DIMENSION(:), PRIVATE, SAVE :: MOL_INFO

    CONTAINS

        TYPE(ATOM) FUNCTION GET_ATOMS(START, FINISH) RESULT(ATOMS)
            INTEGER, INTENT(IN) :: START, FINISH

            ATOMS(:) = MOL_INFO(START:FINISH)

            RETURN
        END FUNCTION GET_ATOMS

        LOGICAL FUNCTION SET_ATOMS(START, FINISH, ATOMS)
            INTEGER, INTENT(IN)                  :: START, FINISH
            TYPE(ATOM), DIMENSION(:), INTENT(IN) :: ATOMS

            MOL_INFO(START:FINISH) = ATOMS(:)
            IF (!!!ARRAYS ARE NOT EQUAL!!! ) THEN
                SET_ATOMS = .TRUE.
            ELSE
                SET_ATOMS = .FALSE.
            END IF

            RETURN
        END FUNCTION SET_ATOMS

END MODULE MOLINFO
