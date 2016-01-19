!
! Apply the projection operator for the totally symmetric IR of point group I
! to vector X, dimension 3*NATOMS
!
SUBROUTINE PROJI(X,NATOMS)
USE COMMONS,ONLY : FROZEN,PERMUTE,DEBUG
IMPLICIT NONE
INTEGER NATOMS, J1, J2, J3, J4, JMATCH
INTEGER NMATCH(NATOMS)
DOUBLE PRECISION X(3*NATOMS), RMAT(3,3,60), XDUMMY(3*NATOMS), DIFF, MINDIST, DIST, DX, DY, DZ
DATA RMAT / 1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0/

XDUMMY(1:3*NATOMS)=0.0D0
DO J1=1,60
!  PRINT *,'RMAT ',J1
!  PRINT '(9G12.5)',RMAT(1:3,1:3,J1)
   DO J2=1,NATOMS
      IF (FROZEN(J2)) CYCLE
      J3=3*(J2-1)
      J4=3*(PERMUTE(J1,J2)-1)
      XDUMMY(J4+1)=XDUMMY(J4+1)+RMAT(1,1,J1)*X(J3+1)+RMAT(1,2,J1)*X(J3+2)+RMAT(1,3,J1)*X(J3+3)
      XDUMMY(J4+2)=XDUMMY(J4+2)+RMAT(2,1,J1)*X(J3+1)+RMAT(2,2,J1)*X(J3+2)+RMAT(2,3,J1)*X(J3+3)
      XDUMMY(J4+3)=XDUMMY(J4+3)+RMAT(3,1,J1)*X(J3+1)+RMAT(3,2,J1)*X(J3+2)+RMAT(3,3,J1)*X(J3+3)
   ENDDO
ENDDO
XDUMMY(1:3*NATOMS)=XDUMMY(1:3*NATOMS)/60.0D0
DIFF=0.0D0
DO J1=1,3*NATOMS
   DIFF=DIFF+(X(J1)-XDUMMY(J1))**2
ENDDO
DO J1=1,13
   J2=3*(J1-1)
!  PRINT '(I8,6G20.10))',J1,X(J2+1),X(J2+2),X(J2+3),XDUMMY(J2+1),XDUMMY(J2+2),XDUMMY(J2+3)
ENDDO
IF (DEBUG) PRINT '(A,G20.10)','projI> RMS difference=',SQRT(DIFF/(3*NATOMS))
X(1:3*NATOMS)=XDUMMY(1:3*NATOMS)

END SUBROUTINE PROJI
!
! Apply the projection operator for the totally symmetric IR of point group I
! to vector X, dimension 3*NATOMS
!
SUBROUTINE PROJIINIT(X,NATOMS)
USE COMMONS, ONLY : PERMUTE
IMPLICIT NONE
INTEGER NATOMS, J1, J2, J3, J4, JMATCH
INTEGER NMATCH(NATOMS)
DOUBLE PRECISION X(3*NATOMS), RMAT(3,3,60), XDUMMY(3*NATOMS), DIFF, MINDIST, DIST, DX, DY, DZ
DATA RMAT / 1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.000000000000D0,   0.000000000000D0,   1.000000000000D0,  &
           -1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,   0.000000000000D0,  -1.000000000000D0,  &
            1.000000000000D0,   0.000000000000D0,   0.000000000000D0,  &
            0.000000000000D0,  -1.000000000000D0,   0.000000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,   0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.309016994375D0,   0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
           -0.500000000000D0,  -0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,  -0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,   0.309016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,   0.809016994375D0,  &
           -0.809016994375D0,   0.309016994375D0,  -0.500000000000D0,  &
           -0.309016994375D0,   0.500000000000D0,   0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0,  &
            0.809016994375D0,  -0.309016994375D0,   0.500000000000D0,  &
            0.309016994375D0,  -0.500000000000D0,  -0.809016994375D0,  &
            0.500000000000D0,   0.809016994375D0,  -0.309016994375D0/

ALLOCATE(PERMUTE(60,NATOMS))
DO J1=1,60
   NMATCH(1:NATOMS)=0
   DO J2=1,NATOMS
      J3=3*(J2-1)
      MINDIST=1.0D100
!
!  Where does this atom transform to under operation J1?
!
      DX=RMAT(1,1,J1)*X(J3+1)+RMAT(1,2,J1)*X(J3+2)+RMAT(1,3,J1)*X(J3+3)
      DY=RMAT(2,1,J1)*X(J3+1)+RMAT(2,2,J1)*X(J3+2)+RMAT(2,3,J1)*X(J3+3)
      DZ=RMAT(3,1,J1)*X(J3+1)+RMAT(3,2,J1)*X(J3+2)+RMAT(3,3,J1)*X(J3+3)
!
!  Which atom does this match in the original structure?
!
      DO J3=1,NATOMS
         DIST=(DX-X(3*(J3-1)+1))**2+(DY-X(3*(J3-1)+2))**2+(DZ-X(3*(J3-1)+3))**2
         IF (DIST.LT.MINDIST) THEN
            JMATCH=J3
            MINDIST=DIST
         ENDIF
      ENDDO
!     PRINT '(A,I8,A,I8,A,I8)','projI> operation ',J1,' atom ',J2,' maps to ',JMATCH
      NMATCH(JMATCH)=NMATCH(JMATCH)+1
      PERMUTE(J1,J2)=JMATCH
   ENDDO
   DO J2=1,NATOMS-1
      IF (NMATCH(J2).NE.NMATCH(J2+1)) THEN
         PRINT '(A)','projI> ERROR - NMATCH:'
         PRINT '(10I6)',NMATCH(1:NATOMS)
         STOP
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE PROJIINIT
