!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!js850> OVERLAP: calculate the overlap between two structures.
!The overlap is defined by dividing space into cubic cells of size DR, and
!counting the number of occupied cells in common between the two structures.
!overlap of 1 means the strucures are the same, or very similar.  overlap of 0
!means they are significantly different structurally.  This is not a robust
!comparison if there is translational or rotational symmetry, so it is really
!only useful if some particles are frozen, or the symmetry is broken in some
!other manner.  There are other definitions of overlap which can fix this
!problem
!
! usage: initialize the module with the first configuration to compare using
! OVERLAP_SETUP.  Then compare against this configuration using
! OVERLAP_GET_OVERLAP.
!

MODULE CLASS_OVERLAP
      USE COMMONS , ONLY : NATOMS, NFREEZE, FREEZE, BOXLX, BOXLY, BOXLZ, FROZENLIST, NFREEZETYPEA
      IMPLICIT NONE
      !SAVE
      PRIVATE
      INTEGER :: OVERLAP_NBINSX
      INTEGER :: OVERLAP_NBINSY
      INTEGER :: OVERLAP_NBINSZ
      INTEGER :: OVERLAP_NBINSZ2
      INTEGER :: OVERLAP_NBINSTOT
      DOUBLE PRECISION :: OVERLAP_DX
      DOUBLE PRECISION :: OVERLAP_DY
      DOUBLE PRECISION :: OVERLAP_DZ
      DOUBLE PRECISION, ALLOCATABLE :: OVERLAP_COORDS1(:)
      !declare a type to hold multiple instances of the same data structures
      TYPE OVERLAP_DATA
        INTEGER N
        INTEGER, ALLOCATABLE :: LIST(:)
        LOGICAL, ALLOCATABLE :: OCCUPIED_CELLS(:)
      END TYPE OVERLAP_DATA
      TYPE(OVERLAP_DATA) :: OVERLAP_A
      TYPE(OVERLAP_DATA) :: OVERLAP_B
      TYPE(OVERLAP_DATA) :: OVERLAP_RA
      TYPE(OVERLAP_DATA) :: OVERLAP_RB


      PUBLIC :: OVERLAP_SETUP, OVERLAP_GET_OVERLAP, OVERLAP_GET_OVERLAP2, OVERLAP_GET_OVERLAP2_R

      CONTAINS

      !put x,y,z in range [0,boxL)
      SUBROUTINE OVERLAP_APPLY_PERIODIC(X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: X,Y,Z
      X=MOD(MOD(X,BOXLX)+BOXLX,BOXLX)
      Y=MOD(MOD(Y,BOXLY)+BOXLY,BOXLY)
      Z=MOD(MOD(Z,BOXLZ)+BOXLZ,BOXLY)
      END SUBROUTINE OVERLAP_APPLY_PERIODIC

      FUNCTION OVERLAP_XYZ2INDEX(X,Y,Z) RESULT(IND)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: X,Y,Z !INOUT so that x, y, z are modifyable
      INTEGER IND, IX, IY, IZ
      CALL OVERLAP_APPLY_PERIODIC(X,Y,Z)
      IX = INT(X/OVERLAP_DX)
      IY = INT(Y/OVERLAP_DY)
      IZ = INT(Z/OVERLAP_DZ)
      IND = IX + IY*OVERLAP_NBINSY + IZ*OVERLAP_NBINSZ2
      END FUNCTION OVERLAP_XYZ2INDEX

      SUBROUTINE OVERLAP_POPULATE(COORDS, A)
      !set OCCUPIED_CELLS for A
      IMPLICIT NONE
      TYPE(OVERLAP_DATA), INTENT(INOUT) :: A
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      INTEGER J1,J2, IND
      DOUBLE PRECISION X,Y,Z

        !WRITE(*,*) "OPN", A%N
      DO J2=1,A%N
        J1 = A%LIST(J2)
        !WRITE(*,*) "j1", J1, j2
        X = COORDS(3*(J1-1)+1)
        Y = COORDS(3*(J1-1)+2)
        Z = COORDS(3*(J1-1)+3)

        IND = OVERLAP_XYZ2INDEX(X,Y,Z)
        !WRITE(*,"(A,I10,3F20.10)"), "IND", IND, X, Y, Z
        A%OCCUPIED_CELLS(IND) = .TRUE.
      ENDDO
      END SUBROUTINE OVERLAP_POPULATE

      FUNCTION OVERLAP_TEST_OVERLAP( X, Y, Z, A ) RESULT(RES)
      IMPLICIT NONE
      TYPE(OVERLAP_DATA), INTENT(IN) :: A
      DOUBLE PRECISION, INTENT(INOUT) :: X,Y,Z !INOUT so that x, y, z are modifyable
      INTEGER IND
      LOGICAL RES
      IND = OVERLAP_XYZ2INDEX( X, Y, Z );
      RES = A%OCCUPIED_CELLS( IND )
      END FUNCTION OVERLAP_TEST_OVERLAP

      FUNCTION OVERLAP_GET( COORDS, A ) RESULT(Q)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      TYPE(OVERLAP_DATA), INTENT(IN) :: A
      DOUBLE PRECISION Q
      INTEGER NCOMPARE
      INTEGER J1,J2
      DOUBLE PRECISION :: X,Y,Z
      NCOMPARE=0
      DO J2=1,A%N
        J1 = A%LIST(J2)
        X = COORDS(3*(J1-1)+1)
        Y = COORDS(3*(J1-1)+2)
        Z = COORDS(3*(J1-1)+3)
        IF ( OVERLAP_TEST_OVERLAP(X,Y,Z,A) ) THEN
          NCOMPARE = NCOMPARE + 1
        ENDIF
      ENDDO
      Q = DBLE(NCOMPARE)/A%N
      END FUNCTION OVERLAP_GET

      !a second definition of overlap between structures R1 and R2
      ! Q = sum_i sum_j exp(-(R1(i)-R2(j))^2 / a^2 ) / NATOMS
      ! If R1 and R2 are identical then Q will be 1 plus some very small
      ! contributions of order NEIBS*exp(-1/L2) from the nearest neighbors
      FUNCTION OVERLAP_CALC_OVERLAP2( COORDS, A, B, L ) RESULT(Q)
      !test the atoms in A agains the atoms in B.  They don't have to be the
      !same length
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), L
      TYPE(OVERLAP_DATA), INTENT(IN) :: A, B
      DOUBLE PRECISION Q, L2
      INTEGER J1,J2, J3, J4, NORM
      DOUBLE PRECISION :: X,Y,Z, R2
      L2 = L*L
      !0.35 seems to be a good value for L2 for a system of 64 free particles.
      !L2 could in principle scale with SIGBB, or SIGAA
      !L2 = (0.35D0)**2 
      Q = 0.D0
      DO J3=1,A%N
        J1 = A%LIST(J3)
        DO J4=1,B%N
          J2 = B%LIST(J4)
          X = COORDS(3*(J2-1)+1) - OVERLAP_COORDS1(3*(J1-1)+1)
          Y = COORDS(3*(J2-1)+2) - OVERLAP_COORDS1(3*(J1-1)+2)
          Z = COORDS(3*(J2-1)+3) - OVERLAP_COORDS1(3*(J1-1)+3)
          R2 = X**2 + Y**2 + Z**2
          !write(*,*) "R ", SQRT(R2)
          Q = Q + EXP( - R2/L2)
        ENDDO
      ENDDO
      NORM = MIN(A%N, B%N)
      Q = Q/NORM
      END FUNCTION OVERLAP_CALC_OVERLAP2

      SUBROUTINE OVERLAP_COPY_RESTRICT_R(COORDS, A, B, Rmax)
      !Using B as a template, set set LIST and N for A, but only for atoms not
      !further than Rmax from the origin  
      IMPLICIT NONE
      TYPE(OVERLAP_DATA), INTENT(INOUT) :: A
      TYPE(OVERLAP_DATA), INTENT(IN) :: B
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), RMAX
      INTEGER J1,J2, IND
      DOUBLE PRECISION X,Y,Z, R2, R2MAX
      R2MAX = RMAX*RMAX
      IF ( .NOT. ALLOCATED(A%LIST ) ) THEN
        ALLOCATE(A%LIST( B%N )) !USE B%N AS AN UPPER BOUND ON THE NUMBER OF ATOMS
      ENDIF
      A%N = 0
      DO J2=1,B%N
        J1 = B%LIST(J2)
        !WRITE(*,*) "j1", J1, j2
        X = COORDS(3*(J1-1)+1)
        Y = COORDS(3*(J1-1)+2)
        Z = COORDS(3*(J1-1)+3)
        R2 = X*X + Y*Y + Z*Z
        IF ( R2 .LE. R2MAX ) THEN
          A%N = A%N + 1
          A%LIST( A%N ) = J1
          !write(*,*) "add ", j1, x,y,z, R2, R2max
        ENDIF
      ENDDO
      END SUBROUTINE OVERLAP_COPY_RESTRICT_R

      !
      ! the subroutines above this opperate on general intances of, or are
      ! independent of TYPE(OVERLAP_DATA) objects.  The subroutines below this 
      ! depend on specifics: OVERLAP_A, OVERLAP_B, etc
      !

      SUBROUTINE OVERLAP_SETUP(COORDS, DR_I )
      USE COMMONS, ONLY : NTYPEA
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: DR_I, COORDS(3*NATOMS)
      INTEGER J1, J2, IND
      DOUBLE PRECISION X,Y,Z
      !these declarations are only because I need NTYPEA
      DOUBLE PRECISION RMAX

      !put COORDS into OVERLAP_COORDS1
      IF (.NOT. ALLOCATED(OVERLAP_COORDS1) ) THEN
        ALLOCATE(OVERLAP_COORDS1(3*NATOMS))
      ENDIF
      OVERLAP_COORDS1(1:3*NATOMS) = COORDS(1:3*NATOMS)

      !try to make drx dry drz as close as possible to DR_I while making sure
      !that boxlx is a multiple of drx
      OVERLAP_NBINSX = BOXLX/DR_I + 1
      OVERLAP_NBINSY = BOXLY/DR_I + 1
      OVERLAP_NBINSZ = BOXLZ/DR_I + 1
      OVERLAP_NBINSZ2 = OVERLAP_NBINSZ**2
      OVERLAP_NBINSTOT = OVERLAP_NBINSX* OVERLAP_NBINSY* OVERLAP_NBINSZ

      OVERLAP_DX = BOXLX/OVERLAP_NBINSX 
      OVERLAP_DY = BOXLY/OVERLAP_NBINSY 
      OVERLAP_DZ = BOXLZ/OVERLAP_NBINSZ 

      !allocate and initialize OVERLAP_OCCUPIED_CELLSA
      IF (.NOT. ALLOCATED(OVERLAP_A%OCCUPIED_CELLS) ) THEN
        !might fail if previously allocated with different values of nbins
        ALLOCATE(OVERLAP_A%OCCUPIED_CELLS(OVERLAP_NBINSTOT))
      ENDIF
      IF (.NOT. ALLOCATED(OVERLAP_B%OCCUPIED_CELLS) ) THEN
        !might fail if previously allocated with different values of nbins
        ALLOCATE(OVERLAP_B%OCCUPIED_CELLS(OVERLAP_NBINSTOT))
      ENDIF
      OVERLAP_A%OCCUPIED_CELLS(:) = .FALSE.
      OVERLAP_B%OCCUPIED_CELLS(:) = .FALSE.

      IF (FREEZE .AND. .NOT. ALLOCATED(FROZENLIST) ) THEN
        WRITE(*,*) "overlap> error: frozenlist not allocated"
        STOP 
      ENDIF


      !get the number of mobile type A and mobile type B atoms
      IF ( FREEZE ) THEN
        OVERLAP_A%N = NTYPEA - NFREEZETYPEA
        !WRITE(*,*) "NA", OVERLAP_A%N,NTYPEA - NFREEZETYPEA
        OVERLAP_B%N = (NATOMS - NFREEZE) - OVERLAP_A%N
      ELSE
        OVERLAP_A%N = NTYPEA
        OVERLAP_B%N = NATOMS - OVERLAP_A%N
      ENDIF
        !WRITE(*,*) "NA1", OVERLAP_A%N,NTYPEA - NFREEZETYPEA
      !allocate spae for LIST
      IF (.NOT. ALLOCATED(OVERLAP_A%LIST) ) THEN
        ALLOCATE(OVERLAP_A%LIST(OVERLAP_A%N))
      ENDIF
      IF (.NOT. ALLOCATED(OVERLAP_B%LIST) ) THEN
        ALLOCATE(OVERLAP_B%LIST(OVERLAP_B%N))
      ENDIF
      !fill OVERLAP_A%LIST OVERLAP_B%LIST witht a list of mobile A and B atoms
      IF ( FREEZE ) THEN
        OVERLAP_A%LIST(:) = FROZENLIST(NFREEZE+1:NFREEZE+OVERLAP_A%N)
        OVERLAP_B%LIST(:) = FROZENLIST(NFREEZE+OVERLAP_A%N+1:NATOMS)
      ELSE
        DO J1=1,OVERLAP_B%N
          OVERLAP_A%LIST(J1) = J1
        ENDDO
        DO J1=1,OVERLAP_B%N
          OVERLAP_B%LIST(J1) = J1+NTYPEA
        ENDDO
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !set up OVERLAP_RA and OVERLAP_RB, with the A and B atoms in the central
      !core.  i.e. with those atoms less than Rcore from the center.
      !Work from the lists in OVERLAP_A and OVERLAP_B
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RMAX = 1.5D0
      CALL OVERLAP_COPY_RESTRICT_R(COORDS, OVERLAP_RA, OVERLAP_A, Rmax)
      CALL OVERLAP_COPY_RESTRICT_R(COORDS, OVERLAP_RB, OVERLAP_B, Rmax)
      ALLOCATE(OVERLAP_RA%OCCUPIED_CELLS(OVERLAP_NBINSTOT))
      ALLOCATE(OVERLAP_RB%OCCUPIED_CELLS(OVERLAP_NBINSTOT))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !OVERLAP_NINSERTA=0
      !OVERLAP_NINSERTB=0

        !WRITE(*,*) "NA2", OVERLAP_A%N,NTYPEA - NFREEZETYPEA
      CALL OVERLAP_POPULATE(COORDS, OVERLAP_A)
      CALL OVERLAP_POPULATE(COORDS, OVERLAP_B)
      CALL OVERLAP_POPULATE(COORDS, OVERLAP_RA)
      CALL OVERLAP_POPULATE(COORDS, OVERLAP_RB)

      WRITE(*,*) "NA", OVERLAP_A%N
      WRITE(*,*) "NB", OVERLAP_B%N
      WRITE(*,*) "DX", OVERLAP_DX, OVERLAP_NBINSX
      WRITE(*,*) "DY", OVERLAP_DY, OVERLAP_NBINSY
      WRITE(*,*) "DZ", OVERLAP_DZ, OVERLAP_NBINSZ
      !WRITE(*,*) "NB", OVERLAP_B%N
      WRITE(*,*) "NA_R", OVERLAP_RA%N
      WRITE(*,*) "NB_R", OVERLAP_RB%N
      !WRITE(*,*) "OVERLAP_RA%LIST"
      !DO J2=1,OVERLAP_RA%N
        !write(*,*)  "    ", OVERLAP_RA%LIST(J2)
      !enddo
      !WRITE(*,*) "OVERLAP_RB%LIST"
      !DO J2=1,OVERLAP_RB%N
        !write(*,*)  "    ", OVERLAP_RB%LIST(J2)
      !enddo

      END SUBROUTINE


      SUBROUTINE OVERLAP_GET_OVERLAP(COORDS, QA, QB, QAB )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: QA, QB, QAB
      QA = OVERLAP_GET(COORDS, OVERLAP_A)
      QB = OVERLAP_GET(COORDS, OVERLAP_B)
      !QAB = weighted average
      QAB = ( QA*OVERLAP_A%N + QB*OVERLAP_B%N) / (OVERLAP_A%N+ OVERLAP_B%N)
      END SUBROUTINE OVERLAP_GET_OVERLAP

      SUBROUTINE OVERLAP_GET_OVERLAP2(COORDS, QA, QB, QAB, L )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), L
      DOUBLE PRECISION, INTENT(OUT) :: QA, QB, QAB
      QA = OVERLAP_CALC_OVERLAP2(COORDS, OVERLAP_A, OVERLAP_A, L)
      QB = OVERLAP_CALC_OVERLAP2(COORDS, OVERLAP_B, OVERLAP_B, L)
      !QAB = weighted average
      QAB = ( QA*OVERLAP_A%N + QB*OVERLAP_B%N) / (OVERLAP_A%N+ OVERLAP_B%N)
      END SUBROUTINE OVERLAP_GET_OVERLAP2

      SUBROUTINE OVERLAP_GET_OVERLAP2_R(COORDS, QA, QB, QAB, L )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), L
      DOUBLE PRECISION, INTENT(OUT) :: QA, QB, QAB
      QA = OVERLAP_CALC_OVERLAP2(COORDS, OVERLAP_RA, OVERLAP_A, L)
      QB = OVERLAP_CALC_OVERLAP2(COORDS, OVERLAP_RB, OVERLAP_B, L)
      !QAB = weighted average
      QAB = ( QA*OVERLAP_RA%N + QB*OVERLAP_RB%N) / (OVERLAP_RA%N+ OVERLAP_RB%N)
      END SUBROUTINE OVERLAP_GET_OVERLAP2_R

END MODULE CLASS_OVERLAP

SUBROUTINE PTMC_OVERLAP_DUMP(XYZ, OVERLAP_UNIT, MYUNIT, OVERLAP_COUNT, &
   total_TIME, quench_time, VNEW, IMCSTEP)
   USE COMMONS, ONLY : NATOMS, COORDS, MYNODE, NQ, DEBUG
   USE CLASS_OVERLAP
   USE PORFUNCS
   implicit none
   integer, intent(IN) :: OVERLAP_UNIT, myunit
   integer, intent(INOUT) :: OVERLAP_COUNT
   double precision, intent(IN) :: xyz(3*natoms), vnew, imcstep
   double precision, intent(OUT) :: total_time, quench_time
   DOUBLE PRECISION OVERLAP_VALA, OVERLAP_VALB, OVERLAP_VALAB
   DOUBLE PRECISION OVERLAP_VAL2A, OVERLAP_VAL2B, OVERLAP_VAL2AB
   DOUBLE PRECISION OVERLAP_VAL2RA, OVERLAP_VAL2RB, OVERLAP_VAL2RAB
   DOUBLE PRECISION OVERLAP_VALQA, OVERLAP_VALQB, OVERLAP_VALQAB
   DOUBLE PRECISION OVERLAP_VALQ2A, OVERLAP_VALQ2B, OVERLAP_VALQ2AB
   DOUBLE PRECISION OVERLAP_VALQ2RA, OVERLAP_VALQ2RB, OVERLAP_VALQ2RAB
   DOUBLE PRECISION timestart, dummy, time1, time2
   integer LBFGS_ITERATIONS, converged, ndummy

   DOUBLE PRECISION POTEL
   COMMON /MYPOT/ POTEL
   CALL MYCPU_TIME(time1)
   quench_time = 0.d0
   total_time = 0.d0

   if (DEBUG) write(myunit,*) "checking overlap"
   OVERLAP_COUNT = OVERLAP_COUNT +1
   CALL OVERLAP_GET_OVERLAP( XYZ, OVERLAP_VALA, OVERLAP_VALB, OVERLAP_VALAB)
   CALL OVERLAP_GET_OVERLAP2( XYZ, OVERLAP_VAL2A, OVERLAP_VAL2B, OVERLAP_VAL2AB, 0.35D0)
   CALL OVERLAP_GET_OVERLAP2_R( XYZ, OVERLAP_VAL2RA, OVERLAP_VAL2RB, OVERLAP_VAL2RAB, 0.35D0)
   IF ( .true. ) THEN
   !
   ! Print the quenched overlap also.  If a normal PTMC run then
   ! quench, using COORDS as scratch space. The quenched energy
   ! will be in POTEL because QUENCH has access to POTEL through a
   ! common block
   !
   ! copy XYZ into coords to use coords as scratch space for doing the quench
   !
      CALL MYCPU_TIME(TIMESTART)
      COORDS(:,MYNODE+1) = XYZ(:)
      CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1))
      NQ(MYNODE+1)=NQ(MYNODE+1)+1
      quench_time = (DUMMY-TIMESTART)
      !IF (.TRUE.) THEN !print some info about the quench
         !IF (CONVERGED.NE.1) WRITE(MYUNIT, '(A)') 'bspt> WARNING - quench did not converge' 
         !WRITE(MYUNIT,'(A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'bspt> overlap> E=', &
            !POTEL,' steps=',LBFGS_ITERATIONS,' RMS=',RMS,' unquenched E=',VNEW,' ttot=',OVERLAP_TIMETOT
      !ENDIF
      CALL OVERLAP_GET_OVERLAP( COORDS(1:3*NATOMS, MYNODE+1), OVERLAP_VALQA, OVERLAP_VALQB, OVERLAP_VALQAB )
      CALL OVERLAP_GET_OVERLAP2( COORDS(1:3*NATOMS, MYNODE+1), OVERLAP_VALQ2A, &
                OVERLAP_VALQ2B, OVERLAP_VALQ2AB, 0.2D0 )
      CALL OVERLAP_GET_OVERLAP2_R( COORDS(1:3*NATOMS, MYNODE+1), OVERLAP_VALQ2RA, &
                OVERLAP_VALQ2RB, OVERLAP_VALQ2RAB, 0.2D0 )
      WRITE(OVERLAP_UNIT,"(F15.1,4F20.10,14F9.4)") IMCSTEP, OVERLAP_VALAB, VNEW, OVERLAP_VALQAB, POTEL, &
           OVERLAP_VALA, OVERLAP_VALB, &
           OVERLAP_VAL2A, OVERLAP_VAL2B, OVERLAP_VAL2AB, &
           OVERLAP_VALQ2A, OVERLAP_VALQ2B, OVERLAP_VALQ2AB, &
           OVERLAP_VAL2RA, OVERLAP_VAL2RB, OVERLAP_VAL2RAB, &
           OVERLAP_VALQ2RA, OVERLAP_VALQ2RB, OVERLAP_VALQ2RAB
   ELSE
      WRITE(OVERLAP_UNIT,"(F15.1,2F20.10)") IMCSTEP, OVERLAP_VALAB, VNEW
   ENDIF
   CALL MYCPU_TIME(time2)
   total_time = time2 - time1
   CALL FLUSH(OVERLAP_UNIT)
end subroutine ptmc_overlap_dump
