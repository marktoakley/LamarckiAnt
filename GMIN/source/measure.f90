MODULE MEASURE
! khs26: Created the module to measure distances, angles and dihedral angles using a generic set of interfaces.

   IMPLICIT NONE

   ! This set of functions calculates the distance between a set of two points.
   INTERFACE DISTANCE
      MODULE PROCEDURE DOUBLE_DISTANCE
      MODULE PROCEDURE ATOMS_DISTANCE
   END INTERFACE DISTANCE

   ! This set of functions calculates the angle between a set of three points or
   ! between two vectors.
   INTERFACE ANGLE
      MODULE PROCEDURE DOUBLE_CART_ANGLE
      MODULE PROCEDURE DOUBLE_VECTOR_ANGLE
      MODULE PROCEDURE ATOMS_ANGLE
   END INTERFACE ANGLE

   ! This set of functions calculates the dihedral angle between a set of four points or
   ! between three vectors.
   INTERFACE DIHEDRAL
      MODULE PROCEDURE DOUBLE_CART_DIHEDRAL
      MODULE PROCEDURE DOUBLE_VECTOR_DIHEDRAL
      MODULE PROCEDURE ATOMS_DIHEDRAL
   END INTERFACE DIHEDRAL

CONTAINS

   ! Cross product
   FUNCTION CROSS_PRODUCT(A, B) RESULT(C)
      ! This function is just the cross product of two vectors.  It gives A x B = C.
      DOUBLE PRECISION, DIMENSION(1:3), INTENT(IN)    :: A, B
      DOUBLE PRECISION, DIMENSION(1:3)                :: C
      
      C(1) = (A(2) * B(3)) - (A(3) * B(2))
      C(2) = (A(3) * B(1)) - (A(1) * B(3))
      C(3) = (A(1) * B(2)) - (A(2) * B(1))

      RETURN
   END FUNCTION CROSS_PRODUCT

   ! Distance functions
   FUNCTION DOUBLE_DISTANCE(X1, Y1, Z1, X2, Y2, Z2) RESULT(DISTANCE)
      ! This one takes 6 doubles (2 sets of 3 Cartesian coordinates) and returns the magnitude of the distance as a double.
      DOUBLE PRECISION, INTENT(IN)  :: X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION              :: DISTANCE

      DISTANCE = SQRT((X1 - X2)**2 + (Y1 - Y2)**2 + (Z1 - Z2)**2)
      RETURN
   END FUNCTION DOUBLE_DISTANCE

   FUNCTION ATOMS_DISTANCE(ATOM_I, ATOM_J, COORDS, NATOMS) RESULT(DISTANCE)
      ! This one takes calculates the distance between the atoms with index ATOM_I and ATOM_J from COORDS and returns the magnitude as a double.
      INTEGER, INTENT(IN)                          :: ATOM_I, ATOM_J, NATOMS
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)   :: COORDS
      DOUBLE PRECISION                             :: DISTANCE
     
      DISTANCE = SQRT((COORDS(3 * ATOM_I - 2) - COORDS(3 * ATOM_J - 2))**2 + &
                    & (COORDS(3 * ATOM_I - 1) - COORDS(3 * ATOM_J - 1))**2 + &
                    & (COORDS(3 * ATOM_I)     - COORDS(3 * ATOM_J))**2)
      RETURN 
   END FUNCTION ATOMS_DISTANCE

   ! Angle functions
   FUNCTION DOUBLE_CART_ANGLE(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3) RESULT(ANGLE)
   ! This function takes 9 doubles (3 sets of 3 Cartesian coordinates) and returns the angle 1-2-3 as a double.
      DOUBLE PRECISION, INTENT(IN)     :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
      DOUBLE PRECISION                 :: ANGLE
      ! This is a set of vectors
      DOUBLE PRECISION, DIMENSION(1:3) :: D21, D23, CROSS

      D21(1) = X1 - X2
      D21(2) = Y1 - Y2
      D21(3) = Z1 - Z2
      D23(1) = X3 - X2
      D23(2) = Y3 - Y2
      D23(3) = Z3 - Z2

      CROSS = CROSS_PRODUCT(D21, D23)

      ANGLE = ATAN2(SQRT(DOT_PRODUCT(CROSS, CROSS)), DOT_PRODUCT(D21, D23))
      RETURN      
   END FUNCTION DOUBLE_CART_ANGLE

   FUNCTION DOUBLE_VECTOR_ANGLE(X1, Y1, Z1, X2, Y2, Z2) RESULT(ANGLE)
   ! This function takes 6 doubles (2 sets of 3 Cartesian vector directions) and returns the angle between the two vectors as a double.
      DOUBLE PRECISION, INTENT(IN)     :: X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION                 :: ANGLE
      ! This is a set of vectors
      DOUBLE PRECISION, DIMENSION(1:3) :: V1, V2, CROSS
   
      V1(1) = X1
      V1(2) = Y1
      V1(3) = Z1
      V2(1) = X2
      V2(2) = Y2
      V2(3) = Z2

      CROSS = CROSS_PRODUCT(V1, V2)

      ANGLE = ATAN2(SQRT(DOT_PRODUCT(CROSS, CROSS)), DOT_PRODUCT(V1, V2))
      RETURN      
   END FUNCTION DOUBLE_VECTOR_ANGLE

   FUNCTION ATOMS_ANGLE(ATOM_I, ATOM_J, ATOM_K, COORDS, NATOMS) RESULT(ANGLE)
   ! This function calculates the angle I-J-K between the atoms with index ATOM_I, ATOM_J and ATOM_K and returns the angle as a double.
      INTEGER, INTENT(IN)                          :: ATOM_I, ATOM_J, ATOM_K, NATOMS
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)   :: COORDS
      DOUBLE PRECISION                             :: ANGLE
      ! This is a set of vectors
      DOUBLE PRECISION, DIMENSION(1:3)             :: D_JI, D_JK, CROSS

      D_JI(1) = COORDS(3 * ATOM_I - 2) - COORDS(3 * ATOM_J - 2)
      D_JI(2) = COORDS(3 * ATOM_I - 1) - COORDS(3 * ATOM_J - 1)
      D_JI(3) = COORDS(3 * ATOM_I)     - COORDS(3 * ATOM_J)
      D_JK(1) = COORDS(3 * ATOM_K - 2) - COORDS(3 * ATOM_J - 2)
      D_JK(2) = COORDS(3 * ATOM_K - 1) - COORDS(3 * ATOM_J - 1)
      D_JK(3) = COORDS(3 * ATOM_K)     - COORDS(3 * ATOM_J)

      CROSS = CROSS_PRODUCT(D_JI, D_JK)

      ANGLE = ATAN2(SQRT(DOT_PRODUCT(CROSS, CROSS)), DOT_PRODUCT(D_JI, D_JK))
      RETURN
   END FUNCTION ATOMS_ANGLE

   FUNCTION DOUBLE_CART_DIHEDRAL(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4) RESULT(DIHEDRAL)
   ! This function takes 12 doubles (4 sets of 3 Cartesian coordinates) and returns the dihedral angle 1-2=3-4 as a double.
      DOUBLE PRECISION, INTENT(IN)  :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4
      DOUBLE PRECISION              :: DIHEDRAL
      ! Vectors and cross products
      DOUBLE PRECISION, DIMENSION(1:3) :: D12, D23, D34, CROSS_12_23, CROSS_23_34
      ! Args for atan2
      DOUBLE PRECISION                 :: FIRST, SECOND

      D12(1) = X2 - X1
      D12(2) = Y2 - Y1
      D12(3) = Z2 - Z1
      D23(1) = X3 - X2
      D23(2) = Y3 - Y2
      D23(3) = Z3 - Z2
      D34(1) = X4 - X3
      D34(2) = Y4 - Y3
      D34(3) = Z4 - Z3
      
      CROSS_12_23 = CROSS_PRODUCT(D12, D23)
      CROSS_23_34 = CROSS_PRODUCT(D23, D34)
      
      ! Equation:
      ! See FUNCTION DOUBLE_VECTOR_DIHEDRAL below.
      FIRST = SQRT(DOT_PRODUCT(D23, D23)) * DOT_PRODUCT(D12, CROSS_23_34)
      SECOND = DOT_PRODUCT(CROSS_12_23, CROSS_23_34)

      DIHEDRAL = ATAN2(FIRST, SECOND) 
      RETURN
   END FUNCTION DOUBLE_CART_DIHEDRAL

   FUNCTION DOUBLE_VECTOR_DIHEDRAL(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3) RESULT(DIHEDRAL)
   ! This function takes 9 doubles (3 sets of 3 Cartesian vector directions) and returns the dihedral angle between 1 and 3 about 2 as a double.
      DOUBLE PRECISION, INTENT(IN)     :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
      DOUBLE PRECISION                 :: DIHEDRAL
      ! Vectors and cross products
      DOUBLE PRECISION, DIMENSION(1:3) :: V1, V2, V3, CROSS_12, CROSS_23
      ! Args for atan2
      DOUBLE PRECISION                 :: FIRST, SECOND

      V1(1) = X1
      V1(2) = Y1
      V1(3) = Z1
      V2(1) = X2
      V2(2) = Y2
      V2(3) = Z2
      V3(1) = X3
      V3(2) = Y3
      V3(3) = Z3

      CROSS_12 = CROSS_PRODUCT(V1, V2)
      CROSS_23 = CROSS_PRODUCT(V2, V3)

      ! Equation:
      ! dihedral = atan2(|V2|(V1 . [V2 x V3]), [V1 x V2] . [V2 x V3])
      FIRST = SQRT(DOT_PRODUCT(V2, V2)) * DOT_PRODUCT(V1, CROSS_23)
      SECOND = DOT_PRODUCT(CROSS_12, CROSS_23)

      DIHEDRAL = ATAN2(FIRST, SECOND) 
      RETURN
   END FUNCTION DOUBLE_VECTOR_DIHEDRAL

   FUNCTION ATOMS_DIHEDRAL(ATOM_I, ATOM_J, ATOM_K, ATOM_L, COORDS, NATOMS) RESULT(DIHEDRAL)
   ! This function calculates the dihedral angle I-J=K-L between the atoms with index ATOM_I, ATOM_J, ATOM_K and ATOM_L and returns the angle as a double.
      INTEGER, INTENT(IN)                          :: ATOM_I, ATOM_J, ATOM_K, ATOM_L, NATOMS
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)   :: COORDS
      DOUBLE PRECISION                             :: DIHEDRAL
      ! Vectors and cross products
      DOUBLE PRECISION, DIMENSION(1:3)             :: D_IJ, D_JK, D_KL, CROSS_IJ_JK, CROSS_JK_KL
      ! Args for atan2
      DOUBLE PRECISION                             :: FIRST, SECOND
      
      D_IJ(1) = COORDS(3 * ATOM_J - 2) - COORDS(3 * ATOM_I - 2)
      D_IJ(2) = COORDS(3 * ATOM_J - 1) - COORDS(3 * ATOM_I - 1)
      D_IJ(3) = COORDS(3 * ATOM_J)     - COORDS(3 * ATOM_I)
      D_JK(1) = COORDS(3 * ATOM_K - 2) - COORDS(3 * ATOM_J - 2)
      D_JK(2) = COORDS(3 * ATOM_K - 1) - COORDS(3 * ATOM_J - 1)
      D_JK(3) = COORDS(3 * ATOM_K)     - COORDS(3 * ATOM_J)
      D_KL(1) = COORDS(3 * ATOM_L - 2) - COORDS(3 * ATOM_K - 2)
      D_KL(2) = COORDS(3 * ATOM_L - 1) - COORDS(3 * ATOM_K - 1)
      D_KL(3) = COORDS(3 * ATOM_L)     - COORDS(3 * ATOM_K)
      
      CROSS_IJ_JK = CROSS_PRODUCT(D_IJ, D_JK)
      CROSS_JK_KL = CROSS_PRODUCT(D_JK, D_KL)

      ! Equation:
      ! See FUNCTION DOUBLE_VECTOR_DIHEDRAL above.
      FIRST = SQRT(DOT_PRODUCT(D_JK, D_JK)) * DOT_PRODUCT(D_IJ, CROSS_JK_KL)
      SECOND = DOT_PRODUCT(CROSS_IJ_JK, CROSS_JK_KL)

      DIHEDRAL = ATAN2(FIRST, SECOND) 
      RETURN
   END FUNCTION ATOMS_DIHEDRAL

END MODULE MEASURE
