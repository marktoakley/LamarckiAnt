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
!
!  Conjugate gradient driver. 
!  CFLAG convergence test
!  CTEST checks for changes in chirality for AMBER runs
!
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
      USE MODHESS  
#ifdef __SPARSE
      USE MODSPARSEHESS
      USE SHIFT_HESS
      USE INERTIA_MOD
#endif /* __SPARSE */
      use COMMONS
      USE MODAMBER9, ONLY : cisarray1, cisarray2, chiarray1, chiarray2, dihedralsave, atomindex, exclude,
     &                      setchiral, setchiralgeneric, nocistransdna, nocistransrna
      USE QMODULE
      use porfuncs
      USE CHIRALITY, ONLY: CIS_TRANS_CHECK, CHIRALITY_CHECK
      IMPLICIT NONE

      INTEGER I, J1, NSQSTEPS, NP, IFLAG, ITER, NOPT, J2, NDUMMY, CSMIT, J5, NUM_ZERO_EVS, HORDER, NMOL
      DOUBLE PRECISION P(3*NATOMS),POTEL,TIME,EREAL,RBCOORDS(18),TMPCOORDS(3*NATOMS), DIST, QE, QX, AVVAL, CSMRMS
      LOGICAL QTEST, CFLAG, RES, COMPON, EVAPREJECT, EVAP, PASS, FAIL
      DOUBLE PRECISION  GRAD(3*NATOMS), DUMMY, DUM(3*NATOMS), DISTMIN, SSAVE, DIST2, RMAT(3,3), EVALUES(3*NATOMS), IT(3,3)
      DOUBLE PRECISION LOG_PROD_EV, LOG_PROD_FREQ
      DOUBLE PRECISION, ALLOCATABLE :: QS(:)
!     DOUBLE PRECISION  WORK(60*NATOMS)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0
      DOUBLE PRECISION RRX,RRY,RRZ,RRR,ITDET
! khs26> Added these to calculate the separation of zero frequencies for free energy BH
      DOUBLE PRECISION LARGEST_ZERO, SMALLEST_NONZERO
      INTEGER ATTEMPTS
      LOGICAL TS_FOUND
      DOUBLE PRECISION INERTIA_TENS(3, 3)
      DOUBLE PRECISION SHIFT_ARRAY(6)

      CHARACTER (LEN=20) QUENCHNUM
      CHARACTER(LEN=80) DSTRING
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET, DUMMYL
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      ! ds656> EVAP was missing from here before 6/9/13
      COMMON /EV/ EVAP, EVAPREJECT
      DOUBLE PRECISION QSTART, QFINISH
      COMMON /Q4C/ QSTART, QFINISH
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT

!
!   sf344> gradually changing parameters to prevent dissociation of PY ellipsoids with repulsive sites 
!
      DOUBLE PRECISION epssave(3)

!
!  Data for the screen saver.
!
      INTEGER BRUN, QDONE,ii
!     the following required to call the LAPACK routine DSYEV
      INTEGER          :: INFO
      INTEGER, PARAMETER :: LWORK=10000 ! the dimension is set arbitrarily
      DOUBLE PRECISION :: WORK(LWORK)

! Timing variables
      DOUBLE PRECISION  :: T_LBFGS_START, T_LBFGS_END, T_DIAG_START, T_DIAG_END

!
!  Turn on guiding potentials. These get turned off in potential.F when
!  the RMS force is small enough.
!

      SSAVE=STEP(NP)
!
! csw34 Reset the NFIX counter
!
      NFIX=0

11    IF (WELCH) TOSI=.TRUE.
      IF (PACHECO) AXTELL=.FALSE.
      IF (CPMD) SCT=.TRUE.
      IF (ZETT1.OR.ZETT2) THEN
         MORSET=.TRUE.
         RHO=6.0D0
      ENDIF
      IF (PERCOLATET) THEN
        COMPON=.TRUE.
      ENDIF
      IF (NATBT.AND.GUPTAT) GUIDET=.TRUE.
      IF (NATBT.AND.GUIDET) GUPTAT=.TRUE.
      IF (DFTBCT.AND.LJATT) GUIDET=.TRUE.
      IF (DFTBCT.AND.GUIDET) THEN
         LJATT=.TRUE.
         IF (DEBUG) WRITE(MYUNIT,'(A)') 'quench> Turning on LJAT guiding potential and rescaling coordinates'
         COORDS(1:3*NATOMS,NP)=COORDS(1:3*NATOMS,NP)/LJATTOC
      ENDIF
      IF (CSMGUIDET) CSMDOGUIDET=.TRUE.
      NOPT=3*NATOMS
      IF (WENZEL) NOPT=2
      IF (MULLERBROWNT) NOPT=2
!
!  QTEST is set for the final quenches with tighter convergence criteria.
!
      IF (QTEST) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      QDONE=0
      DO I=1,3*NATOMS
         P(I)=COORDS(I,NP)
      ENDDO
      IF (QCIPOTT) THEN ! initial setup
         CALL CHECKREPQCIPOT(P,3*NATOMS,0,1)
         CALL QCIPOT(EREAL,P,GRAD)
      ENDIF

!
!     IF (TIP) THEN
!        WRITE(DUMPXYZUNIT(NP),'(I6)') (NATOMS/2)*3
!        WRITE(DUMPXYZUNIT(NP),70) NP,NQ(NP), EREAL, RMS
!        DO J2=1,NATOMS/2
!           CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
!    1           P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
!           WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
!           WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
!           WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
!        ENDDO
!     ENDIF


      IF (COMPRESST.AND.(.NOT.QTEST)) THEN
         COMPON=.TRUE.
         IF (PATCHY) THEN
            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,1.D1*GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ELSE
            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         END IF
         POTEL=EREAL
         IF (.NOT.CFLAG) WRITE(MYUNIT,'(A,I7,A)') ' WARNING - compressed quench ',NQ(NP),'  did not converge'
         WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,F15.7,A,I4,A,F12.2)') 'Comp Q ',NQ(NP),' energy=',
     1              POTEL,' steps=',ITER,' RMS force=',RMS
      ENDIF

      IF (.NOT.PERCOLATET) COMPON=.FALSE.


10    IF (PERMOPT.OR.PERMINVOPT.OR.DISTOPT) THEN ! lb415
         !IF ( NQ(NP) .eq. 1) THEN
         IF (DUMPT) THEN
            IF (NP.EQ.1) WRITE(MYUNIT,'(A,4I6)') 'quench> initial NP,DUMPXYZUNIT=',NP,DUMPXYZUNIT(NP)
            WRITE(DUMPXYZUNIT(NP),'(I6)') NDUMMY
            WRITE(DUMPXYZUNIT(NP),'(A,I6)') 'quench> initial points before quench ',NQ(NP)
            WRITE(DUMPXYZUNIT(NP),'(A,3G20.10)') ('LA ',P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),J2=1,NATOMS)
         ENDIF
         CALL POTENTIAL(P,GRAD,EREAL,.FALSE.,.FALSE.)
         CFLAG=.TRUE.
!        ITER=1
!        RMS=0.0D0
         RMS=CSMRMS
         ITER=CSMIT
         !ELSE
         !   CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP) ! minimize structure
         !   write(*,*) 'permdist mylbfgs', EREAL, ITER, RMS
         !   POTEL=EREAL
         !   IF (.NOT.CFLAG) WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
         !   DO II=1,NSAVE
         !      IF ( II .GE. NQ(NP) ) EXIT ! There's no need to check further, there's nothing
         !      CALL MINPERMDIST(P,QMINP(II,:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID,RMAT)
         !      write(*,*) DUMMY, 'dummy',ii
         !      IF (DUMMY .LT. 0.5D0) THEN
         !         !DO NOT ACCEPT THIS QUENCH
         !         WRITE(MYUNIT,*) 'This quench ended in a known minimum. It won`t be counted.'
         !         RETURN
         !      ENDIF
         !   ENDDO
         !ENDIF 
      ELSEIF (MODEL1T) THEN
         CALL MODEL1(P,GRAD,EREAL,QE,QX)
         EREAL=QE
         CFLAG=.TRUE.
         ITER=1
         RMS=0.0D0
         P(1)=QX
      ELSE IF (DL_POLY) THEN
!
!  Need to make DL_POLY input file for current coordinates.
!
         OPEN (UNIT=91,FILE='CONFIG',STATUS='OLD')
         OPEN (UNIT=92,FILE='config',STATUS='UNKNOWN')
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         DO J1=1,NATOMS
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(3G20.10)') P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3)
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
         ENDDO
         CLOSE(91)
         CLOSE(92)
         CALL SYSTEM('cp CONFIG CONFIG.old; cp config CONFIG')
         CALL SYSTEM('DLPOLY.X > output.DL_POLY ; tail -9 STATIS > energy')
         OPEN (UNIT=91,FILE='energy',STATUS='OLD')
         READ(91,*) EREAL
         WRITE(MYUNIT,'(A,G20.10)') 'energy=',EREAL
         CLOSE(91)
         OPEN(UNIT=91,FILE='REVCON',STATUS='OLD')
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
         NATOMS=0
13       READ(91,'(A1)',END=14) DUMMY
         NATOMS=NATOMS+1
         READ(91,*) P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
!        WRITE(MYUNIT,'(3G20.10)') P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         GOTO 13
14       CONTINUE
         CLOSE(91)
         CFLAG=.TRUE.
!
!  Read the coordinates of the minimised geometry into vector P.
!
!     ELSE IF (BFGS .AND.(.NOT.QTEST)) THEN
      ELSE IF (BFGS) THEN
!        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
         CALL DFPMIN(MAXIT,P,3*NATOMS,GMAX,ITER,EREAL,CFLAG)
      ELSEIF (TNT) THEN
!        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
          WRITE(MYUNIT, '(A)') 'subroutine tn does not compile with NAG/PG'
         STOP
!        CALL TN(IFLAG,3*NATOMS,P,EREAL,GRAD,WORK,60*NATOMS,GMAX,ITER,MAXIT,CFLAG,DEBUG)
      ELSEIF (CONJG) THEN
         CALL CGMIN(MAXIT,P,CFLAG,ITER,EREAL,NP)
    ! 
! Compute quantum energy with Variation Gaussian Wavepacket.
! Coords are scaled by VGW LJ sigma (LJSIGMA) inputed with VGW params.
! Coords are then scaled back to unit sigma.
! 
      ELSEIF (VGW) THEN    
        IF(QTEST) THEN              
          CALL VGWQUENCH(P,EREAL,CFLAG)
          ELSE
            CALL VGWQUENCHSP(P,EREAL,CFLAG)
        ENDIF 

      ELSEIF (MYSDT) THEN
         CALL MYSD(MAXIT,P,CFLAG,ITER,EREAL)
      ELSEIF (RKMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSEIF (BSMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSE
!        CALL CGMIN(5,P,CFLAG,ITER,EREAL,NP)
         IF (CHRMMT.AND.INTMINT) THEN
            CALL MYLBFGS(NINTS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ELSE IF (THOMSONT .AND. (.NOT. GTHOMSONT)) THEN
            TMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)
            CALL THOMSONCTOANG(TMPCOORDS,P,NATOMS)
            CALL MYLBFGS(2*NATOMS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
            CALL THOMSONANGTOC(P,NATOMS)
         ELSE IF (THOMSONT .AND. GTHOMSONT ) THEN
! jwrm2> For the initial quench, the coordinates come in as Cartesians,
!        during the run as spherical polars. The conversion routines should
!        check before converting the coordinates.
            CALL GTHOMSONCTOANG(COORDS(1:3*NATOMS,1),P(1:3*NATOMS),NATOMS)
            CALL MYLBFGS(3*NATOMS,MUPDATE,P(1:3*NATOMS),.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!            CALL GTHOMSONANGTOC(TMPCOORDS(1:3*NATOMS), P(1:3*NATOMS), NATOMS)
!            P(1:3*NATOMS) = TMPCOORDS(1:3*NATOMS)
!         ELSE IF(PYBINARYT) THEN
!! sf344> trying out some sort of systematic parameter change to prevent particles from dissociating:
!! first decrease repulsive epsilon values, converge, then gradually increase them
!           epssave(:)=PEPSILON1(:)
!          IF(.NOT.QTEST) THEN
!           WRITE(MYUNIT,*) 'first iteration: decreasing epsilon_rep values by a factor of 10000' 
!           PEPSILON1(:)=PEPSILON1(:)/10000.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!           WRITE(MYUNIT,*) 'second iteration: increasing epsilon_rep values by a factor of 100' 
!           PEPSILON1(:)=PEPSILON1(:)*100.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!            WRITE(MYUNIT,*) 'third iteration: increasing epsilon_rep values by a factor of 100' 
!           PEPSILON1(:)=PEPSILON1(:)*100.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!          END IF
         ELSEIF ((.NOT.QCIPOTT).AND.INTCONSTRAINTT) THEN
            WRITE(MYUNIT,'(A,2L5)') 'QCIPOTT,INTCONSTRAINTT=',QCIPOTT,INTCONSTRAINTT
            CALL CPU_TIME(T_LBFGS_START)
            CALL INTLBFGS(P,FINISH)
            CALL CPU_TIME(T_LBFGS_END)
            STOP
         ELSE
! Added timing for LBFGS call
            CALL CPU_TIME(T_LBFGS_START)
            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
            CALL CPU_TIME(T_LBFGS_END)
            IF (FEBHT.AND.DEBUG) WRITE(MYUNIT, '(A, F10.3)') 'Time to minimise (s):', T_LBFGS_END - T_LBFGS_START
         ENDIF
         IF (EVAPREJECT) RETURN
      ENDIF
      !
      IF (FEBHT .AND. CFLAG) THEN
      ! Calculate the free energy
         NUM_ZERO_EVS=6
         IF (NATOMS.EQ.2) NUM_ZERO_EVS=5
         IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
         ALLOCATE(HESS(3*NATOMS,3*NATOMS))
         CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
         CALL MASSWT()
#ifdef __SPARSE
         IF (SPARSET) THEN
             IF (DEBUG) THEN
                WRITE(MYUNIT, *) 'Using sparse cholesky decomposition to calculate determinant.'
             END IF
             SHIFT_ARRAY(:) = 1.0D0
             CALL CPU_TIME(T_DIAG_START)
             CALL SHIFT_HESS_ZEROS(P, SHIFT_ARRAY)
             CALL FILTER_ZEROS(HESS, ZERO_THRESH)
             CALL GET_DETERMINANT(3*NATOMS, LOG_PROD_EV, NQ(NP))
             CALL CPU_TIME(T_DIAG_END)
             IF (DEBUG) WRITE(MYUNIT, '(A, F10.3)') 'Time to diagonalise with sparse routines (s):', T_DIAG_END - T_DIAG_START
             IF (LOG_PROD_EV < SMALL_DOUBLE) THEN
                WRITE(MYUNIT, '(A)') 'Log product from SuiteSparse is -infinity. The hessian ' //
     &                               'is probably not positive definite. If using a threshold ' //
     &                               'on zero values, try reducing it.'
!
! Changed to convergence failure so that the run can continue, as for ts location. DJW
!
                CFLAG=.FALSE.
                RETURN
             END IF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For debugging only, this *should* be very slow for cases where we use the sparse approach: 
!             PRINT *, "Sparse:", LOG_PROD_EV
!             CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
!             CALL MASSWT()
!             CALL SHIFT_HESS_ZEROS(P, (/ 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0 /))
!             CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)             
!             PRINT *, "DSYEV:", SUM(LOG(EVALUES))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ELSE
#endif /* __SPARSE */
            CALL CPU_TIME(T_DIAG_START)
            CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)
            CALL CPU_TIME(T_DIAG_END)
            IF (DEBUG) WRITE(MYUNIT, '(A, F10.3)') 'Time to diagonalise with DSYEV (s):', T_DIAG_END - T_DIAG_START
            IF (DEBUG) THEN
               WRITE(MYUNIT, *) 'Using DSYEV to calculate determinant.'
            END IF
            IF (DEBUG) THEN
                WRITE(MYUNIT, '(A)') "Eigenvalues:"
                WRITE(MYUNIT, '(3F12.3)') EVALUES
            END IF
! Check that we don't have any proper negative eigenvalues (the magnitude should drop).
! We do it this way, since some of the zeros come out negative as well, so you can't 
! just check the sign.
            TS_FOUND = (ABS(EVALUES(1)) > ABS(EVALUES(NUM_ZERO_EVS + 1)))
            IF (TS_FOUND) THEN
                WRITE(MYUNIT, '(A, I10, A)') 'Quench ', NQ(NP), ' converged to a transition state.'
                IF (DEBUG) THEN
                   WRITE(MYUNIT, '(A)') 'Eigenvalues'
                   WRITE(MYUNIT, '(A)') '======================='
                   WRITE(MYUNIT, '(6F20.12)') EVALUES(1:(3+NUM_ZERO_EVS))
                   WRITE(MYUNIT, '(A)') '======================='
                END IF
                CFLAG=.FALSE.
                RETURN
            END IF
            LOG_PROD_EV = SUM(DLOG(EVALUES(NUM_ZERO_EVS + 1:)))
            IF ((MIN_ZERO_SEP .GT. 0.0D0) .AND. (MAX_ATTEMPTS .GT. 0)) THEN
! Square the minimum zero separation, since this refers to separation of frequencies.
                DO ATTEMPTS = 1, MAX_ATTEMPTS            
! Identify the largest zero and smallest non-zero
                    EVALUES(:) = ABS(EVALUES(:))
                    LARGEST_ZERO = 0.D0
                    LARGEST_ZERO = MAX(LARGEST_ZERO, MAXVAL(EVALUES(1:NUM_ZERO_EVS)))
                    SMALLEST_NONZERO = 1.0D100
                    SMALLEST_NONZERO = MIN(SMALLEST_NONZERO, MINVAL(EVALUES(NUM_ZERO_EVS+1:)))
! If the separation of zeros and non-zeros is too small, reduce the convergence
! threshold by an order of magnitude, along with the corresponding sloppy or tight
! convergence threshold (for future quenches).
                    IF ((SMALLEST_NONZERO / LARGEST_ZERO) < MIN_ZERO_SEP) THEN
                        WRITE(MYUNIT, '(A,I8)') 'Attempt: ', ATTEMPTS
                        WRITE(MYUNIT, '(A,E15.7)') 'Current separation is ', (SMALLEST_NONZERO / LARGEST_ZERO)
                        WRITE(MYUNIT, '(A,E15.7)') 'Target separation is ', MIN_ZERO_SEP
                        GMAX = GMAX * 1.0D-1
                        IF (QTEST) THEN
                            CQMAX = GMAX
                            WRITE(MYUNIT, '(A,F12.8)') 'Lowering tight convergence to ', CQMAX
                        ELSE
                            BQMAX = GMAX
                            WRITE(MYUNIT, '(A,F12.8)') 'Lowering sloppy convergence to ', BQMAX
                        END IF
                        CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
                        CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
                        CALL MASSWT()
                        CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)
                        EVALUES(:) = ABS(EVALUES(:))
! Recalculate separation.
                        LARGEST_ZERO = 0.D0
                        LARGEST_ZERO = MAX(LARGEST_ZERO, MAXVAL(EVALUES(1:NUM_ZERO_EVS)))
                        SMALLEST_NONZERO = 1.0D100
                        SMALLEST_NONZERO = MIN(SMALLEST_NONZERO, MINVAL(EVALUES(NUM_ZERO_EVS+1:)))
                    END IF
                    IF ((SMALLEST_NONZERO / LARGEST_ZERO) >= MIN_ZERO_SEP) THEN
                        ! If we have converged properly, then we're ok.
                        IF (DEBUG) THEN
                           WRITE(MYUNIT, '(A,E15.7)') 'Converged. Separation of zeros: ', (SMALLEST_NONZERO / LARGEST_ZERO)
                        END IF
                        EXIT
                    END IF
                    IF (ATTEMPTS == MAX_ATTEMPTS) THEN
                        WRITE(MYUNIT, '(A)') 'Failed to achieve desired separation of zeros and non-zeros.'
                        WRITE(MYUNIT, '(A)') 'Magnitudes of lowest eigenvalues:'
                        WRITE(MYUNIT, '(F20.12)') EVALUES(1:NUM_ZERO_EVS + 6)
                        STOP 'Cannot achieve desired separation of zeros and non-zeros. Please check your input.'
                    END IF
                END DO
                LOG_PROD_EV = SUM(DLOG(EVALUES(NUM_ZERO_EVS+1:)))
            END IF
#ifdef __SPARSE
         END IF
#endif /* __SPARSE */
! Log product of frequencies is log product of sqrt of evalues
         LOG_PROD_FREQ = 0.5D0 * LOG_PROD_EV 

! Continue with FEBH calculation
         ITDET=0.0D0
         IF (NTIPT) THEN
            ALLOCATE(QS(9*NATOMS/2))
            J1=NRBSITES
            NRBSITES=3
            CALL SITEPOS(P,QS)
            NRBSITES=J1
            NMOL =NATOMS/2
            NATOMS=(NATOMS/2)*3
            CALL DETSYMMETRY(QS,HORDER,IT,.FALSE.)
            NATOMS=2*NMOL
            DEALLOCATE(QS)
         ELSE 
            CALL DETSYMMETRY(P,HORDER,IT,.FALSE.)
            IF (USEROT) CALL MYINERTIA(P,ITDET)
         ENDIF
         IF (DEBUG) THEN
            IF (QTEST) WRITE(MYUNIT,'(A,I4)') 'Order of the point group of minimum=', HORDER
         ENDIF
         IF (ITDET.NE.0.0D0) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'quench> Inertia determinant=',ITDET
            ITDET=LOG(ITDET)/2.0D0 
         ENDIF
         FEBH_POT_ENE=EREAL
! khs26> At zero temperature, free energy is potential energy. However, DLOG(0)=NaN
         
         IF (FETEMP .GT. 0.D0) THEN
            IF (DEBUG) WRITE(MYUNIT, '(A,F20.12)') 'Potential energy=', FEBH_POT_ENE
            EREAL=-FETEMP*(SYMFCTR+LOG(1.0D0/HORDER) - EREAL/FETEMP + (3*NATOMS-NUM_ZERO_EVS)*DLOG(FETEMP) - LOG_PROD_EV + ITDET)
            IF (DEBUG) WRITE(MYUNIT, '(A,F20.12)') 'Harmonic superposition contribution=', EREAL - FEBH_POT_ENE
         END IF
      ELSE IF (FEBHT .AND. (.NOT. CFLAG)) THEN
         WRITE(MYUNIT, '(A)') 'Quench did not converge, not calculating free energy and adding 1E10 to energy.'
         EREAL=EREAL + 1.0D10
      ENDIF

! ds656> If (relative) chemical potentials are present, then we tack
!        them on to the (free) energy for semi-grand canonical (FE)BH.
      IF (SEMIGRAND_MUT) THEN
         DUMMY=0.0D0
         DO J1=2, NSPECIES(0)
            DUMMY = DUMMY + NSPECIES(J1)*SEMIGRAND_MU(J1)
         ENDDO
         EREAL = EREAL - DUMMY
         IF(FEBHT) FEBH_POT_ENE = FEBH_POT_ENE - DUMMY 
      ENDIF

      POTEL=EREAL
!fh301>{{{ 
! for more effective restraints just compare the CamShift part of the energies
!      IF ((CHEMSHIFT2) .AND. (DABS(ENERGYCAMSHIFT).GE.0.0001)) THEN
!        FULLENERGY=POTEL
!        POTEL=ENERGYCAMSHIFT
!      ENDIF
!fh301>}}}

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ELSE
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ENDIF
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

! khs26> Write free energies to an output file after the quench.
! Quench=NQ(NP) 
! Potential energy=FEBH_POT_ENE
! Free energy=EREAL
! Markov energy=EPREV(NP)
! Harmonic superposition contribution=EREAL - FEBH_POT_ENE
! Time=TIME
      IF (FEBHT.AND.(MOD(NQ(NP)-1,PRTFRQ).EQ.0)) THEN
         IF (.NOT. QTEST) THEN
            WRITE(FE_FILE_UNIT, '(I12,8X,4(F18.12,2X),F18.1,2X)')         NQ(NP), FEBH_POT_ENE, EREAL - FEBH_POT_ENE,
     &                                                               EREAL, EPREV(NP), TIME
         ELSE
            WRITE(FE_FILE_UNIT, '(A1,I11,8X,4(F18.12,2X),F18.1,2X)') 'F', NQ(NP), FEBH_POT_ENE, EREAL - FEBH_POT_ENE,
     &                                                               EREAL, EPREV(NP), TIME
         END IF
      END IF

      RES=.FALSE.
      IF (TABOOT.AND.(.NOT.QTEST).AND.(.NOT.RENORM)) THEN
         CALL TABOO(EREAL,POTEL,P,NP,RES)
         IF (RES) GOTO 10
      ENDIF

!     PRINT*,'Taboo lists:'
!     DO J1=1,NPAR
!        PRINT*,'Parallel run ',J1
!        WRITE(*,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
!     ENDDO
!     PRINT*,'Inertia lists:'
!     DO J1=1,NPAR
!        PRINT*,'Parallel run ',J1
!        WRITE(MYUNIT,'(6F15.7)') (XINSAVE(J2,J1),J2=1,NT(J1))
!     ENDDO

! csw34> CHIRALITY AND PEPTIDE BOND CHECKS - reports GOODSTRUCTURE
!        If the checks pass (or are not done!), GSAVEIT is called to 
!        add the quenches structure to QMIN and QMINP if low enough E
      GOODSTRUCTURE=.TRUE.
! csw34> SAVEQ is .TRUE. for quenches, and so the checks will be
!        applied, and structures possibly saved. For final quenches, it 
!        is set to .FALSE. (in finalq.f) and so the checks are skipped.
      IF (SAVEQ) THEN
!
! csw34> CHARMM TESTS
!
         IF(CHRMMT) THEN
!
! csw34> CHIRALITY
!
            FAIL=.FALSE.
            IF (CHECKCHIRALITY) CALL CHECKCHIRAL(P,FAIL)
            IF (FAIL) THEN
               GOODSTRUCTURE=.FALSE.
               WRITE(MYUNIT,*) ' quench> CHIRALITY CHECK FAILED - discarding structure'
            ENDIF
!
! csw34> PEPTIDE BOND
!
            FAIL=.FALSE.
            IF (NOCISTRANS) CALL CHECKOMEGA(P,FAIL)
            IF (FAIL) THEN
               GOODSTRUCTURE=.FALSE.
               WRITE(MYUNIT,*) ' quench> PEPTIDE BOND CHECK FAILED - discarding structure'
            ENDIF
         ENDIF
!
! csw34> AMBER tests
!        The AMBER tests are a bit fancier, they are designed to be
!        check if the chirality and peptide bond geometry has been
!        maintained from the starting structure, not just look at an
!        absolute value. This is done automatically for the cis/trans
!        checks and with by specifying SETCHIRAL for chirality checks
         IF(AMBERT) THEN
!
! csw34> CHIRALITY
!
            IF (CHECKCHIRALITY) THEN
               PASS=.TRUE.
               FAIL=.FALSE.
               IF (SETCHIRALGENERIC) THEN
                  IF (NQ(NP)==0) THEN
!                    WRITE(*,*) 'calling check_chirality_generic, initial call' 
                     ALLOCATE(dihedralsave(NATOMS),atomindex(NATOMS,4),exclude(NATOMS))
                     CALL check_chirality_generic(P,NATOMS,GOODSTRUCTURE,.TRUE.)
                  ELSE
                     WRITE(*,*) 'calling check_chirality_generic' 
                     CALL check_chirality_generic(P,NATOMS,GOODSTRUCTURE,.FALSE.)
                     IF (.NOT.GOODSTRUCTURE) THEN
                        FAIL=.TRUE.
                        PASS=.FALSE.
                     END IF
                  END IF
               END IF
               IF (SETCHIRAL) THEN 
! csw34> First - if SETCHIRAL is specified
                  CALL set_check_chiral(P,NATOMS,DUMMYL,chiarray2)
                  DO J5=1,NATOMS
! csw34> Compare the chiarrays, a change in sign will result in a non
!        zero output and so signal a change in chirality from the
!        initial structure. 
                     IF ((CHIARRAY1(J5)-CHIARRAY2(J5))/=0) FAIL=.TRUE.
                  ENDDO
                  IF (FAIL) THEN
                     WRITE(MYUNIT,*) ' quench> WARNING: chirality differs from initial structure!'
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
! csw34> If SETCHIRAL is NOT specified
               ELSE
                  CALL check_chirality(P,NATOMS,PASS)
                  IF (.NOT.PASS) THEN
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
               ENDIF
! csw34> If either test fails, print a warning
               IF (FAIL.OR.(.NOT.PASS)) THEN 
                  WRITE(MYUNIT,*) ' quench> CHIRALITY CHECK FAILED - discarding structure'
               ENDIF
! csw34> End of CHECKCHIRALITY block
            ENDIF
!
! csw34> PEPTIDE BOND
!        There are seperate checks for DNA, RNA and proteins!
!
            IF (NOCISTRANS) THEN
               PASS=.TRUE.
               FAIL=.FALSE.
! csw34> DNA
               IF (NOCISTRANSDNA) THEN
                  CALL CHECK_CISTRANS_DNA(P,NATOMS,ZSYM,PASS)
                  IF (.NOT.PASS) THEN
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
! csw34> RNA
               ELSEIF (NOCISTRANSRNA) THEN
                  CALL CHECK_CISTRANS_RNA(P,NATOMS,ZSYM,PASS)
                  IF (.NOT.PASS) THEN
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
! csw34> PROTEIN
               ELSE
! csw34> As before, AMBER stores the initial cis/trans setup and
!        compares it to the current structure, flagging up any changes 
                  CALL check_cistrans_protein(P,NATOMS,DUMMYL,MINOMEGA,cisarray2)
                  DO J5=1,NATOMS
! csw34> Compare the cisarrays, a change in sign will result in a non
!        zero output and so signal a change in isomer from the
!        initial structure. 
                     IF ((CISARRAY1(J5)-CISARRAY2(J5))/=0) FAIL=.TRUE.
                  ENDDO
                  IF (FAIL) THEN
                     WRITE(MYUNIT,*) ' quench> WARNING: cis/trans differs from initial structure!'
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
               ENDIF
! csw34> If either test fails, print a warning
               IF (FAIL.OR.(.NOT.PASS)) THEN 
                  WRITE(MYUNIT,*) ' quench> CIS/TRANS CHECK FAILED - discarding structure'
               ENDIF
            ENDIF
         ENDIF

!================================ AMBER 12 =====================================
!
! khs26> AMBER 12 tests
! The AMBER 12 tests compare the chirality and peptide bond geometry of the
! molecule with those at the start, provided that CHECKCHIRALITY, SETCHIRAL
! and NOCISTRANS are .TRUE..
!
! The subroutines for these are contained in chirality.F90. Scripts used to
! identify peptide bonds and chiral centres are contained in 
! $SVN/SCRIPTS/AMBER/chirality or distributed separately, depending on how you
! obtained GMIN.

         IF(AMBER12T) THEN
! 
! khs26> Start of chirality checks
!
            IF (CHECKCHIRALITY) THEN
               PASS=.TRUE.
               IF (SETCHIRAL) THEN 
                  CALL CHIRALITY_CHECK(P, PASS)
                  IF (.NOT. PASS) THEN
                     WRITE(MYUNIT,*) ' quench> WARNING: chirality differs from initial structure!'
                     WRITE(MYUNIT,*) ' quench> CHIRALITY CHECK FAILED - discarding structure'
                     GOODSTRUCTURE=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
!
! khs26> End of chirality checks
!

! 
! khs26> Start of peptide bond checks.
!
            IF (NOCISTRANS) THEN
               PASS=.TRUE.
               FAIL=.FALSE.
               CALL CIS_TRANS_CHECK(P, PASS)
               IF (.NOT. PASS) THEN
                  WRITE(MYUNIT,*) ' quench> WARNING: cis/trans differs from initial structure!'
                  WRITE(MYUNIT,*) ' quench> CIS/TRANS CHECK FAILED - discarding structure'
                  GOODSTRUCTURE=.FALSE.
               END IF
            ENDIF
         ENDIF
!============================ END AMBER 12 =====================================

! js850> Check to see if RESTRICTREGION constraints have been violated - if they have, reject the step
         IF ( RESTRICTREGION ) THEN
           RESTRICTREGIONTEST=.FALSE.
!          J1=13
!          WRITE(*,'(1x,1I5,3F25.16,2I5)') J1,
!    &       ( P(3*(J1-1)+1) ),
!    &       ( P(3*(J1-1)+2) ),
!    &       ( P(3*(J1-1)+3) ), 1, NP
           DO J1=1,NATOMS
             IF ( .NOT. FROZEN(J1) .AND. .NOT. HARMONICFLIST(J1) 
     &       .AND. .NOT. DONTMOVE(J1) ) THEN
               RRX=( P(3*(J1-1)+1)-RESTRICTREGIONX0 )
               RRY=( P(3*(J1-1)+2)-RESTRICTREGIONY0 )
               RRZ=( P(3*(J1-1)+3)-RESTRICTREGIONZ0 )
               RRX=RRX - ANINT(RRX/BOXLX)*BOXLX
               RRY=RRY - ANINT(RRY/BOXLY)*BOXLY
               IF ( RESTRICTCYL ) THEN
                 RRZ=0
                 RRR=DSQRT(RRX**2+RRY**2 )
               ELSE
                 RRZ=RRZ - ANINT(RRZ/BOXLZ)*BOXLZ
                 RRR=DSQRT(RRX**2+RRY**2+RRZ**2 )
               ENDIF
               IF ( RRR > RESTRICTREGIONRADIUS ) THEN
!                WRITE(*,'(1x,1I5,7F25.16)') J1,
!    &             ( P(3*(J1-1)+1) ),
!    &             ( P(3*(J1-1)+2) ),
!    &             ( P(3*(J1-1)+3) ), RRX, RRY, RRZ, RRR
!                WRITE(*,*) "GOODSTRUCTURE", GOODSTRUCTURE
!                 write(*,*) "not saving structure"
                 GOODSTRUCTURE=.FALSE.
                 RESTRICTREGIONTEST=.TRUE.
                 EXIT
               ENDIF
             ENDIF
           ENDDO
         ENDIF

! jwrm2> Check percolation. If the structure is disconnected, don't save it.
         PERCT=.TRUE.
         IF (PERCOLATET) THEN
           IF (GTHOMSONT) THEN
! jwrm2> Need to use Cartesian rather than polar coordinates for this
             CALL GTHOMSONANGTOC(TMPCOORDS(1:3*NATOMS), P(1:3*NATOMS), NATOMS)
             CALL PERC(TMPCOORDS(1:3*NATOMS), NATOMS, PERCCUT, PERCT, DEBUG, MYUNIT, .FALSE.)
           ELSE
             CALL PERC(P,NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)
           END IF
         ENDIF
         
! csw34> If all tests have been passed, save the structure!        
         IF (GOODSTRUCTURE .AND. PERCT .AND. (.NOT. SAVEMULTIMINONLY) .AND. (.NOT. GRADPROBLEMT)) THEN
!fh301>{{{
!           IF (CHEMSHIFT2) THEN
!             IF (DABS(ENERGYCAMSHIFT).GE.0.0001) THEN
!               CALL GSAVEIT(POTEL,P,NP)
!               IF (MONITORT) CALL MSAVEIT(POTEL,P,NP)
!             ENDIF
!           ELSE
!fh301>}}}
             IF(NSPECIES(0) > 1) THEN
                CALL GSAVEIT_MC(EREAL,P,INVATOMLISTS(1:NATOMS,1), NP)
             ELSE
                CALL GSAVEIT(EREAL,P,NP)
             ENDIF
             
             IF (NSAVE==0) QMIN(1)=min(QMIN(1),EREAL)
             IF (MONITORT) CALL MSAVEIT(EREAL,P,NP)
!fh301>{{{
!           ENDIF
!fh301>}}}
         ENDIF

      ENDIF
! csw34> END OF CHIRALITY AND PEPTIDE BOND CHECKS


!     IF (QDONE.EQ.0) THEN
!        PRINT '(A)','WARNING quench did not converge from starting coodinates:'
!        WRITE(MYUNIT,'(3G20.10)') (COORDS(J1,NP),J1=1,3*NATOMS)
!     ENDIF
!
!  If EPSSPHERE is non-zero we are presumably doing a calculation of the 
!  energy density of local minima. We need to know the minimum distance
!  between the starting point and the quenched minima.
!
      IF ((EPSSPHERE.NE.0.0D0).OR.BSWL) THEN
         DO J1=1,3*NATOMS
            DUM(J1)=COORDS(J1,NP)
         ENDDO
!
!  DUM is returned in the closest orientation to P; P should not change.
!  This is nearly the same mind as OPTIM! To execute a random walk we must take 
!  another step and minimise until the distance between the starting point
!  and the quench minimum is less than EPSSPHERE.
!
!        CALL MINDGMIN(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD)
         CALL NEWMINDIST(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
      ENDIF
!
!  Deal with EPSSPHERE sampling.
!
      IF (EPSSPHERE.NE.0.0D0) THEN
         IF ((DISTMIN.GT.EPSSPHERE).OR.(ABS(EREAL-EPREV(NP)).LE.ECONV)) THEN
            WRITE(MYUNIT,'(A,F12.5,A,4F14.5)') 'step ',STEP(NP),' EREAL, EPREV, DISTMIN, EPSSPHERE=',
     1                                     EREAL, EPREV(NP), DISTMIN, EPSSPHERE
            DO J1=1,3*NATOMS
               COORDS(J1,NP)=COORDSO(J1,NP)
            ENDDO
            CALL TAKESTEP(NP)
             WRITE(MYUNIT,'(A,G20.10)' ) 'reseeding step, maximum displacement reset to ',STEP(NP)
            GOTO 11
         ELSE
            WRITE(MYUNIT,'(A,2F20.10)') 'valid step, DISTMIN, EPSSPHERE=',DISTMIN, EPSSPHERE
         ENDIF
      ENDIF
!
!  If we are provided with target minimum coordinates in file coords.target then
!  calculate the minimum distances. May be useful for algorithm development.
!  If we get close, we don;t want to escape without a hit!
!
!     IF (ALLOCATED(TCOORDS)) THEN
!        DO J1=1,NTARGETS
!           TMPCOORDS(1:3*NATOMS)=TCOORDS(J1,1:3*NATOMS)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID)
!           WRITE(MYUNIT, '(A,I5,A,F15.3,A,F15.3,A,F20.10)') 'for target structure ',J1,' dist=',DUMMY,' dist2=',DIST2,' V=',POTEL
!        ENDDO
!        DO J1=1,MIN(NMSBSAVE,MAXSAVE)
!           TMPCOORDS(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J1)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID)
!           PRINT '(A,I5,A,F15.3,A,F15.3,A,F20.10)','for taboo  structure ',J1,' dist=',DUMMY,' dist2=',DIST2,' V=',POTEL
!        ENDDO
!     ENDIF
!
!  NORESET true does not set the configuration point to the quench geometry
!  A relaxed frozen core does not get saved, but the lowest minima are saved
!  by GSAVEIT.
!
      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      IF (CALCQT) CALL ORDERQ4(NATOMS,P,QFINISH)
!
!  Calling CENTRE here without an evaporation check can put particles
!  outside the container, and make a valid step in takestep impossible.
!
!     PRINT*,'Calling centre from quench'
!     IF ((.NOT.FIELDT).AND.(.NOT.SEEDT).AND.CENT) CALL CENTRE2(COORDS(1:3*NATOMS,NP))

      IF (DUMPT) THEN
         IF (ARNO) THEN
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS+2
            WRITE(DUMPXYZUNIT(NP),70) NP,NQ(NP),EREAL,RMS
            WRITE(DUMPXYZUNIT(NP),'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(DUMPXYZUNIT(NP),'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            IF (NS.NE.0) WRITE(DUMPXYZUNIT(NP),65) (P(I),I=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TIP) THEN
            WRITE(DUMPVUNIT(NP),'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT(NP),'(I6)') (NATOMS/2)*3
            WRITE(DUMPXYZUNIT(NP),70) NP,NQ(NP), EREAL, RMS
            DO J2=1,NATOMS/2
               CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
     1              P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ELSE IF (CHRMMT) THEN
            CALL CHARMMDUMP3(P)
            CALL CHARMMDUMP2(P,DUMPXYZUNIT(NP)) ! xyz
         ELSEIF (NCORE(NP).GT.0) THEN
            WRITE(DUMPVUNIT(NP),'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
            WRITE(DUMPXYZUNIT(NP),70) NQ(NP), EREAL, RMS
!           WRITE(DUMPXYZUNIT(NP),80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NCORE(NP))
!           WRITE(DUMPXYZUNIT(NP),80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NCORE(NP)+1,NATOMS)
            WRITE(DUMPXYZUNIT(NP),80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NCORE(NP))
            IF (NCORE(NP).GT.0) WRITE(DUMPXYZUNIT(NP),80) 
     &                     ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NCORE(NP)+1,NATOMS)
         ELSE IF (DMACRYST) THEN
            WRITE(QUENCHNUM,*) NQ(NP)
            CALL DMACRYS_DUMP('dump.'//TRIM(ADJUSTL(QUENCHNUM))//'.cif', P)
         ELSE
            WRITE(DUMPVUNIT(NP),'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
            WRITE(DUMPXYZUNIT(NP),70) NQ(NP), EREAL, RMS
            WRITE(DUMPXYZUNIT(NP),80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.NE.0) WRITE(DUMPXYZUNIT(NP),80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
70          FORMAT(1X,'QUENCH NUMBER ',I6,' final energy=',F20.10,' RMS force=',E20.10)
80          FORMAT(A2,3F20.10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG DJW
         ENDIF
      ENDIF

      IF (SQUEEZET) THEN
         IF ((EREAL.GT.0.0D0).AND.(SQUEEZED.LT.1.0D0)) THEN
            SQUEEZED=2.0D0-SQUEEZED
            NSQSTEPS=NQ(NP)
         ELSE
            NSQSTEPS=100000
         ENDIF
         DO J1=1,3*NVEC
            VEC(J1)=VEC(J1)*SQUEEZED
         ENDDO
         IF (NQ(NP).GT.2*NSQSTEPS) SQUEEZET=.FALSE.
      ENDIF
    
      IF ((NQ(NP).GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(MYUNIT,'(I6,A,G20.10)') NSSTOP,' quenches completed, setting coordinates to the lowest minimum, E=',QMIN(1)
         DO J1=1,3*QMINNATOMS(1)
            COORDS(J1,NP)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      END SUBROUTINE QUENCH
