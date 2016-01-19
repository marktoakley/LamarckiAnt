      SUBROUTINE INERTIA2(IT,Q)
!     Builds inertia tensor IT from the coordinates Q
      USE COMMONS, ONLY : NATOMS, ATMASS

      IMPLICIT NONE

      INTEGER          :: I, J, K
      DOUBLE PRECISION :: Q(3*NATOMS), IT(3,3), R(3) 

      IT(:,:) = 0.D0
      DO K = 1, NATOMS
         R(:) = Q(3*K-2:3*K)
         DO I = 1, 3
            IT(I,I) = IT(I,I) + ATMASS(K)*(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))
            DO J = 1, 3       
               IT(I,J) = IT(I,J) - ATMASS(K)*R(I)*R(J)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE INERTIA2

!     ==============================================================================================

      SUBROUTINE ROTCON(IT,IPRNT,IERR,PTEST)
!
!  COMPUTES ROTATIONAL CONSTANTS (IN CM-1), OPTIONALLY
!  PRINTS THEM AND RETURNS AN ERROR MESSAGE IF INERTIA
!  TENSOR PASSED IS NOT DIAGONAL.
!
      IMPLICIT NONE
      DOUBLE PRECISION IT(3,3),RC(3),PI,FACTOR,Z
      INTEGER NORD(3),IPRNT,IERR,IBOT,I,J
      DATA FACTOR /189.123013D0/
      DATA PI /3.141592654D0/
      LOGICAL PTEST

      Z=0.D0
      IBOT=0
      IERR=0
      DO I=1,3
         DO J=I+1,3
            Z=DABS(IT(I,J))+DABS(IT(J,I))+Z
         ENDDO
      ENDDO
      IF(Z.GT.1.D-9)THEN
       WRITE(6,100)
 100   FORMAT(' ***PROGRAM ERROR***, Inertia tensor not diagonal in ROTCON.')
       PRINT*,'Z=',Z
       IERR=1
       IBOT=1
      ELSE
       IBOT=1
      ENDIF
      DO 20 I=1,3
       IF(IT(I,I).GT.0)THEN
        RC(I)=FACTOR/IT(I,I)
        RC(I)=RC(I)/PI
       ELSE
        IBOT=2
       ENDIF
 20   CONTINUE
      CALL PIKSR2(3,RC,NORD)
      IF(IPRNT.NE.0)THEN
       IF (PTEST) WRITE(6,200)
       IF (PTEST) WRITE(6,300) (RC(J),J=IBOT,3),IT(1,1),IT(2,2),IT(3,3)
 200   FORMAT(' rotcon> Rotational constants (in cm-1) and principal moments of inertia: ')
 300   FORMAT((6(F15.5,1X)))
      ENDIF
      RETURN
      END SUBROUTINE ROTCON

!     ==============================================================================================

!
!      SUBROUTINE CALCVEC(A,B,V,IX)
!C
!C CALCULATES THE VECTOR V BETWEEN CARTESIAN POINTS A AND B
!C
!      IMPLICIT NONE
!      DOUBLE PRECISION A(3),B(3),V(3)
!      INTEGER IX,I
!
!      DO I=1,3
!         V(I)=B(I)-A(I)
!      ENDDO
!      IF(IX.EQ.1)CALL NORMAL(V,3)
!
!      RETURN
!      END SUBROUTINE CALCVEC

      SUBROUTINE VSTAT(V,ZQ,LENGTH,N)
!
! RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
!     ZQ(1)  Largest absolute magnitude
!     ZQ(2)  Smallest absolute magnitude
!     ZQ(3)  Largest value
!     ZQ(4)  Smallest value
!     ZQ(5)  2-norm
!     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
!   

      IMPLICIT NONE
      INTEGER N,I,LENGTH
      DOUBLE PRECISION V(N),ZQ(6),U

      U=0.D0
      ZQ(1:6)=0.0D0
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
!     PRINT*,'Statistics reported for vector:'
!     WRITE(*,10) (J,V(J),J=1,LENGTH)
!10    FORMAT(I4,F20.10)
      DO I=1,LENGTH
         ZQ(1)=MAX(ZQ(1),DABS(V(I)))
         ZQ(2)=MIN(ZQ(2),DABS(V(I)))
         ZQ(3)=MAX(ZQ(3),V(I))
         ZQ(4)=MIN(ZQ(4),V(I))
         U=U+V(I)*V(I)
      ENDDO
      If (Length .ne. 0.0) ZQ(5)=DSQRT(U/LENGTH)
      ZQ(6)=ZQ(2)-ZQ(1)
      RETURN
      END SUBROUTINE VSTAT

       INTEGER FUNCTION ATOI (STRING)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Purpose:      Convert STRING to an integer value
!
! Arguments:    STRING   character string (input only)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! 
! Revision 4.0  89/03/14  01:15:20  bernhold
! Baseline for Sun & VAX prior to porting everywhere
! 
! Revision 3.0  89/01/29  23:09:23  bernhold
! First working release for VAX
! 
! Revision 2.1  89/01/02  20:35:02  bernhold
! To keep consistent with .u file just checked in.
! 
!
! System:       Standard FORTRAN 77
!
! Copyright 1988 David E. Bernholdt
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       CHARACTER(LEN=80) STRING
       CHARACTER(LEN=1) C
       INTEGER I
       LOGICAL NEG
!
       ATOI = 0
       NEG = .FALSE.
       I = 1
       C = STRING(I:I)
!
!      Pass over any leading spaces
!
 100   IF (C .EQ. ' ') THEN
          I = I + 1
          C = STRING(I:I)
          GOTO 100
       ENDIF
!
!      See if first character makes the number negative
!      Accept + as valid character before the digits start
!
       IF (C .EQ. '-') THEN
          NEG = .TRUE.
          I = I + 1
          C = STRING(I:I)
       ELSEIF (C .EQ. '+') THEN
          I = I + 1
          C = STRING(I:I)
       ENDIF
!
!      Continue as long as its a digit...
!
 200   IF (LGE(C, '0') .AND. LLE(C,'9')) THEN
!            Shift number over & add new digit
          ATOI = 10*ATOI + ICHAR(C)-48
          I = I + 1
          C = STRING(I:I)
          GOTO 200
       ENDIF
!
!      Negate the result if necessary
!
       IF (NEG) ATOI = -ATOI
       RETURN
       END FUNCTION ATOI
!*****************************************************
!
! multiply matrices: A^T=B^T * C
!
      FUNCTION DOTOPT(A,B,N)
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION A(N),B(N), DOTOPT

      DOTOPT=0.D0
      DO I=1,N
         DOTOPT=DOTOPT+A(I)*B(I)
      ENDDO
      RETURN
      END FUNCTION DOTOPT

!     ==============================================================================================

      SUBROUTINE CROSSOPT(A,B,C,IX)
!
! CALCULATES THE (OPTIONALLY) NORMALIZED VECTOR CROSS PRODUCT C=A x B
!
      IMPLICIT NONE
      INTEGER IX
      DOUBLE PRECISION A(3),B(3),C(3)

      C(3)=A(1)*B(2)-B(1)*A(2)
      C(2)=-A(1)*B(3)+A(3)*B(1)
      C(1)=A(2)*B(3)-A(3)*B(2)
      IF(IX.EQ.1)CALL NORMAL(C,3)
      RETURN
      END SUBROUTINE CROSSOPT
!
! CALCULATES THE DISTANCE BETWEEN TWO POINTS IN CARTESIAN SPACE
!
      FUNCTION DIST(A,B)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION A(3),B(3),Z,DIST

      Z=0.D0
      DO I=1,3
         Z=Z+(A(I)-B(I))**2
      ENDDO
      DIST=DSQRT(Z)
      RETURN
      END FUNCTION DIST

!
!     ROBUST EQUIVALENCE CHECK - DO WELL DEFINED SORT ON COORDINATE
!     MATRIX AND COMPARE ELEMENT BY ELEMENT.  SHOULD BE FOOLPROOF.
!
!     VEC      coordinate vector to be checked (modified)
!     VECR     sorted reference coordinate vector (input only)
!     NORD     ???
!     ICOMP    number of coordinates outside of TOL (output only)
!     TOL      tolerance for comparison of coords (input only)
!
      SUBROUTINE COMPARE3(VEC,VECR,NORD,ICOMP,TOL)
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE
!     Maximum number of atoms currently allowed
      INTEGER NORD(2*NATOMS+1),I,JAP,ICOMP
      DOUBLE PRECISION VEC(3*NATOMS),VECR(3*NATOMS),Z,TOL
!
      ICOMP=0
!      IF (IPRNT.GT.130) THEN
!         WRITE(6,*)'Input vector before sorting'
!         WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
! 80   FORMAT(3(1X,F10.5))
!      ENDIF
      CALL SORTXYZ2(VEC,VEC,NORD(NATOMS+1),TOL)
!      IF (IPRNT.GT.130) THEN
!         WRITE(6,*)'Comparison vector'
!         WRITE(6,80)(VECR(JAP),JAP=1,3*NATOMS)
!         WRITE(6,*)'Sorted input vector'
!         WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
!      ENDIF
      DO 30 I=1,NATOMS*3
         Z = DABS( VECR(I)-VEC(I) )
!        PRINT*,'I,Z,VECR,VEC=',I,Z,VECR(I),VEC(I)
         IF (Z.GT.TOL) THEN
            ICOMP = ICOMP + 1
            RETURN
         ENDIF
30    CONTINUE
      RETURN
      END SUBROUTINE COMPARE3

      SUBROUTINE MTRANSOPT(A,AT,NR,NC,NTR,NTC)
      IMPLICIT NONE
      INTEGER I,J,NTC,NTR,NR,NC
      DOUBLE PRECISION A(NR,NC),AT(NC,NR)

      DO I=1,NTR
         DO J=1,NTC
            AT(J,I)=A(I,J)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MTRANSOPT
!
! REFLECTS POINTS IN PLANE
!
!     ==============================================================================================

      SUBROUTINE SORTXYZ2(XX,Y,NORD,TOL)
! Doxygen {{{
!> \name SORTXYZ
!> \brief Sort vector of nuclear coordinates 
!> \param[in] XX DP(3*NATOMS)  - input vector to be sorted
!> \param Y  DP(3*NATOMS)      - output sorted vector
!> \param NORD INT(2*NATOMS+1) - vector containing permutations 
!> \param TOL DP - tolerance factor
! }}}
! Declarations {{{ 
      USE COMMONS, ONLY: NATOMS, ATMASS
      IMPLICIT NONE
! subroutine parameters  
      DOUBLE PRECISION XX(3*NATOMS),Y(3*NATOMS),TOL
      INTEGER NORD(2*NATOMS+1)
! local parameters 
!      DOUBLE PRECISION X(3*NATOMS),XX(*),Y(*),TOL
!      INTEGER I, J, JK, NORD(NATOMS), ILINE
      DOUBLE PRECISION X(3*NATOMS)
      INTEGER I, J, JK, ILINE
! }}}
! Subroutine body {{{
!
      ILINE(J)=1+J/3
!
!     Sort on the X - if two X's are equivalent, sort on Y and so on.
!     If the coordinate < the tolerance we should ignore it! However,
!     if the tolerance is sloppy that can lead to the sorting ignoring
!     genuine small differences between coordinates. Sigh. 
!
      DO I=1,3*NATOMS
         X(I)=XX(I)
      ENDDO
!
!     FIRST GIVE DUMMY ATOMS RIDICULOUS SCRATCH COORDINATES - ENSURES
!     THAT THEY WILL WIND UP AT THE BOTTOM OF THE LIST
!
      DO I=1,3*NATOMS-2,3
         IF (DABS(ATMASS(ILINE(I))).LT.1D-3) THEN
            DO J=0,2
               X(J+I) = -99995.0D0
            ENDDO
         ENDIF
      ENDDO
      JK=1
40    J=1
      DO I=1,3*NATOMS-2,3
         IF (X(I)-X(J).GT.TOL) J=I
         IF (DABS(X(I)-X(J)).LT.TOL) THEN
            IF (X(I+1)-X(J+1).GT.TOL) J=I
            IF (DABS(X(I+1)-X(J+1)).LT.TOL) THEN
               IF (X(I+2)-X(J+2).GT.TOL) J=I
            ENDIF
         ENDIF
      ENDDO
!
!     MASS-WEIGHT SORTED VECTOR - WILL ZERO ELEMENTS CORRESPONDING
!     TO DUMMY ATOMS SO THEY DONT MUCK UP THE SYMMETRY CHECKING.
!
      DO I=0,2
         Y(3*JK-2+I)=X(J+I)*NINT(ATMASS(ILINE(J)))
         X(J)=-99999.D0
      ENDDO
      NORD(JK)=(J+2)/3
      JK=JK+1
      IF (JK.EQ.NATOMS+1) GOTO 70
      GOTO 40
70    CONTINUE

! }}}
      RETURN
      END SUBROUTINE SORTXYZ2

