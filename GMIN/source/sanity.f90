MODULE SANITY
! This module is for the sanity check functions we should be using as much as possible.
! The module name may or may not have been choosen to allow us to USE SANITY. 
! Please code nicely in here...GOTOs will be removed, as will unclear variable names.
! Consistent indentation is mandatory!

CONTAINS

FUNCTION CHECK_DIMENSION(ARRAYSIZE, EXPECTED)
! Function to check the dimension of the ARRAY passed in.
! Returns .TRUE. if the dimension is as EXPECTED.
!
! Arguments:
! 
! ARRAYSIZE: The size of the array to be checked
! EXPECTED: The expected dimension of the ARRAY
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: EXPECTED
   INTEGER, INTENT(IN) :: ARRAYSIZE
   LOGICAL             :: CHECK_DIMENSION
   
   CHECK_DIMENSION=.FALSE.

   IF (MODULO(ARRAYSIZE, EXPECTED) == 0) THEN
      CHECK_DIMENSION=.TRUE.
   ENDIF 
END FUNCTION 

END MODULE SANITY
