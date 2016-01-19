! operations.f90 is intended to contain simple transforms and operations that
! are often done elsewhere in the code. It is to be gradually expanded to reduce
! the number of times we are repeating outselves.

! CONTENTS

! FUNCTION PAIRDISTANCE(ATOM1,ATOM2) - returns the distance between atoms given
! their coordinates


!
!> \brief Function to determine the distance between two points in 3D space 
!> \author Chris Whittleston, csw34@cam.ac.uk
!
FUNCTION PAIRDISTANCE(ATOM1,ATOM2)
IMPLICIT NONE
! PAIRDISTANCE is defined as it is the function name and hence the returned
! value
DOUBLE PRECISION :: PAIRDISTANCE 
DOUBLE PRECISION, INTENT(IN) :: ATOM1(3),ATOM2(3)
   
PAIRDISTANCE=DSQRT((ATOM1(1)-ATOM2(1))**2+(ATOM1(2)-ATOM2(2))**2+(ATOM1(3)-ATOM2(3))**2)
RETURN 
END FUNCTION PAIRDISTANCE
