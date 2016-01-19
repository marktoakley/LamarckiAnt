!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE TRANSFORMS
    
CONTAINS

!!! General transformation

    FUNCTION MATRIX_TRANSFORM(INPUT_COORDS, TRANSFORMATION_MATRIX) RESULT(OUTPUT_COORDS)
    ! General transform, requires that the transformation matrix be square and the same size as the input coordinates.
        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(:) :: OUTPUT_COORDS
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: INPUT_COORDS, TRANSFORMATION_MATRIX

!        IF (!the sizes don't match...

!       MATMUL

    END FUNCTION MATRIX_TRANSFORM

!!! Translations start here

    FUNCTION TRANSLATE(INPUT_COORDS, TRANSLATION) RESULT(OUTPUT_COORDS)
    ! Translates the input coordinates by TRANSLATION (X1,Y1,Z1,X2,Y2,Z2...Xn,Yn,Zn)
    ! (x, y, z) -> (x+X, y+Y, z+Z)
        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(:) :: OUTPUT_COORDS
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: INPUT_COORDS, TRANSLATION

        IF (SIZE(INPUT_COORDS) .NE. SIZE(TRANSLATION)) THEN
            STOP 'Error: Sizes of the input coordinates and translation arrays are different in function TRANSLATE of module TRANSFORMS.'
        END IF

        OUTPUT_COORDS = INPUT_COORDS + TRANSLATION
    END FUNCTION TRANSLATE

!!! (Proper) rotations start here

    FUNCTION ROTATE_ABOUT_AXIS_ORIGIN(INPUT_COORDS, ROTATION_AXIS, ROTATION_ANGLE) RESULT(OUTPUT_COORDS)
    ! Rotates the input coordinates through the given angle about the axis provided.
    ! CAN THE INPUT_COORDS VECTOR BE LARGER THAN 3?  NEED TO THINK ABOUT COORDS, SHOULD IT BE AN ARRAY - COORDS(INDEX, ATOM) 

        
    END FUNCTION ROTATE_ABOUT_AXIS_ORIGIN

    FUNCTION ROTATE_ABOUT_XYZ(INPUT_COORDS, ROTATION_ANGLES) RESULT(OUTPUT_COORDS)
    ! Rotates the input coordinates about the x, y and z axes by the angles specified.

    END FUNCTION ROTATE_ABOUT_XYZ

!!! Reflections start here

!!! Inversions start here

    FUNCTION INVERSION_THROUGH_ORIGIN(INPUT_COORDS) RESULT(OUTPUT_COORDS)
    ! Inverts the input coordinates through (0,0,0)
    ! (x, y, z) -> (-x, -y, -z)
        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(:) :: OUTPUT_COORDS
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: INPUT_COORDS

        OUTPUT_COORDS = -INPUT_COORDS
    END FUNCTION INVERSION_THROUGH_ORIGIN

    FUNCTION INVERSION_THROUGH_POINT(INPUT_COORDS, POINT_COORDS) RESULT(OUTPUT_COORDS)
    ! Inverts the input coordinates through POINT_COORDS (X1,Y1,Z1,X2,Y2,Z2...Xn,Yn,Zn)
    ! (x, y, z) -> (2X-x, 2Y-y, 2Z-z)
        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(:) :: OUTPUT_COORDS
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: INPUT_COORDS, POINT_COORDS

        IF (SIZE(INPUT_COORDS) .NE. SIZE(POINT_COORDS)) THEN
            STOP 'Error: Sizes of the input coordinates and the point of inversion coordinates &
                  & are different in function INVERSION_THROUGH_POINT of module TRANSFORMS.'
        END IF

        OUTPUT_COORDS = (2 * POINT_COORDS) - INPUT_COORDS
    END FUNCTION INVERSION_THROUGH_POINT

!!! Improper rotations start here

END MODULE TRANSFORMS
