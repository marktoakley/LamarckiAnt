MODULE MODSPARSEHESS
USE, INTRINSIC :: ISO_C_BINDING
USE MODHESS

INTERFACE
    SUBROUTINE PRINT_HESSIAN(DIMENSIONS) BIND(C, NAME='print_hess')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
        INTEGER(C_INT) :: DIMENSIONS
    END SUBROUTINE PRINT_HESSIAN

    SUBROUTINE GET_DETERMINANT(DIMENSIONS, DETERMINANT, QUENCH_NUM) BIND(C, NAME='get_determinant')
        USE, INTRINSIC     :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
        INTEGER(C_INT)     :: DIMENSIONS
        REAL(C_DOUBLE)     :: DETERMINANT
        INTEGER(C_INT)     :: QUENCH_NUM
    END SUBROUTINE GET_DETERMINANT
END INTERFACE

CONTAINS

    SUBROUTINE GET_POINTER_IN_C(ARRAY_PTR) BIND(C, NAME='get_pointer')
    !  C signature: void get_pointer(double **);
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE, C_LOC
        IMPLICIT NONE
    ! Arguments
        TYPE(C_PTR), INTENT(OUT)    :: ARRAY_PTR
    ! Local variables
        REAL(C_DOUBLE), POINTER     :: HESS_COPY(:, :)

    ! Copy the hessian into HESS_COPY
        ALLOCATE(HESS_COPY(SIZE(HESS, 1), SIZE(HESS, 2)))
        HESS_COPY(:, :) = HESS(:, :)

    ! Pass to C
        ARRAY_PTR = C_LOC(HESS_COPY(1, 1))

    END SUBROUTINE GET_POINTER_IN_C

    SUBROUTINE FREE_POINTER_IN_C(ARRAY_PTR) BIND(C, NAME='free_pointer')
    !  C signature: void free_pointer(double **);
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_DOUBLE
        IMPLICIT NONE
    ! Arguments
        TYPE(C_PTR), INTENT(OUT)    :: ARRAY_PTR
    ! Local variables
        REAL(C_DOUBLE), POINTER     :: HESS_COPY(:, :)

    ! Get the memory allocated to ARRAY_PTR and deallocate it.
        CALL C_F_POINTER(ARRAY_PTR, HESS_COPY, SHAPE(HESS))
        DEALLOCATE(HESS_COPY)
    END SUBROUTINE FREE_POINTER_IN_C

    SUBROUTINE FILTER_ZEROS(ARRAY, MAGNITUDE)    
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
        USE COMMONS, ONLY: DEBUG
        IMPLICIT NONE
    ! Arguments
        REAL(C_DOUBLE), INTENT(INOUT)   :: ARRAY(:, :)
        REAL(C_DOUBLE), INTENT(IN)      :: MAGNITUDE
    ! Local variables
        INTEGER                         :: I, J, TOT
        INTEGER, PARAMETER              :: HESS_OUT = 4738
   
        TOT = 0 
    ! Set all values in ARRAY which are less than MAGNITUDE to 0.
        DO I = 1, SIZE(ARRAY, 1)
            DO J = I, SIZE(ARRAY, 2)
                IF (ABS(ARRAY(J, I)) < MAGNITUDE) THEN
                    ARRAY(J, I) = 0.0D0
                    ARRAY(I, J) = 0.0D0
                    TOT = TOT + 2
                ELSE
                END IF
            END DO
        END DO

    ! Printing, if we use debug
        IF (DEBUG) THEN
            OPEN(HESS_OUT, FILE='hess_sparsity')
            DO I = 1, SIZE(ARRAY, 1)
                DO J = 1, SIZE(ARRAY, 2)
                    IF (ABS(ARRAY(J, I)) > 0.0D0) THEN
                        WRITE(HESS_OUT, '(A1)', ADVANCE='NO') '@'
                    ELSE
                        WRITE(HESS_OUT, '(A1)', ADVANCE='NO') '.'
                    END IF
                END DO
                WRITE(HESS_OUT, '(A1)') ' '
            END DO
            CLOSE(HESS_OUT)
        END IF
        
    END SUBROUTINE

END MODULE MODSPARSEHESS
