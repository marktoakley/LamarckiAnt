      SUBROUTINE CHECKD(X)

      USE COMMONS, ONLY: NATOMS, COMPRESST, PERCOLATET, CHECKDID, GTHOMSONT

      USE MODHESS
      IMPLICIT NONE

      INTEGER          :: IVRNO, IVRNO1, IVRNO2
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY, FM, FP, DFA, DFN, TMPCOORDS(3*NATOMS)
      LOGICAL          :: GTEST, STEST, COMPON
      DOUBLE PRECISION, PARAMETER :: ERRLIM = 1.D-05, DELX = 1.D-06
      COMMON /CO/ COMPON

! jwrm2> Turning compression on, if required
      IF (COMPRESST .OR. PERCOLATET) COMPON = .TRUE.

! jwrm2> Converting GTHOMSON coordinates to polars
      IF (GTHOMSONT) THEN
        CALL GTHOMSONCTOANG(X(1:3*NATOMS), TMPCOORDS(1:3*NATOMS), NATOMS)
        X(1:3*NATOMS) = TMPCOORDS(1:3*NATOMS)
      END IF

      STEST = .FALSE.

      IF (CHECKDID == 0) THEN
         GTEST = .FALSE.
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         WRITE(*, *) 'Energy  = ', ENERGY

      ELSEIF (CHECKDID == 1) THEN

!     Checks gradients

      DO IVRNO = 1, 3*NATOMS

         WRITE(*, *) IVRNO

         GTEST    = .FALSE.
         X(IVRNO) = X(IVRNO) - DELX
         CALL POTENTIAL (X, G, FM, GTEST, STEST)
         WRITE(*, *) 'Energy minus = ', FM

         X(IVRNO) = X(IVRNO) + 2.D0*DELX
         CALL POTENTIAL (X, G,  FP, GTEST, STEST)
         WRITE(*, *) 'Energy plus  = ', FP

         GTEST = .TRUE.
         X(IVRNO) = X(IVRNO) - DELX
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         DFN = (FP - FM) / (2.D0*DELX)
         IF (ABS(DFN) .LT. 1.0D-10) DFN = 0.D0
         DFA = G(IVRNO)

         WRITE(*, *) 'Gradient numerical  = ', DFN
         WRITE(*, *) 'Gradient analytical = ', DFA

         IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) IVRNO, DFN, DFA, ABS(DFN-DFA)

      ENDDO

      ELSE IF (CHECKDID == 2) THEN

         IF (.NOT. ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))

         DO IVRNO1 = 1, 3*NATOMS
            DO IVRNO2 = 1, 3*NATOMS
               WRITE(*,*) IVRNO1, IVRNO2
               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE.)
               FM   = G(IVRNO2)
!              WRITE(*, *) 'Gradient minus = ', FM

               X(IVRNO1) = X(IVRNO1) + 2.D0*DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE.)
               FP   = G(IVRNO2)
!              WRITE(*, *) 'Gradient plus = ', FP

               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.TRUE.)
               DFN  = (FP - FM) / (2.D0*DELX)
               DFA  = HESS(IVRNO1,IVRNO2)

               WRITE(*, *) 'Hessian numerical  = ', DFN
               WRITE(*, *) 'Hessian analytical = ', DFA

               IF (ABS(DFN - DFA) > ERRLIM) WRITE(*,*) 'Error:', IVRNO1, IVRNO2, DFN, DFA, ABS(DFN-DFA)
            ENDDO
         ENDDO
      ENDIF

      STOP

      END SUBROUTINE CHECKD
