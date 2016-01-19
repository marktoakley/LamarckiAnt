!  GMIN: A program for finding global minima
!
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
! 
!---======================================---
      SUBROUTINE PTBASINSAMPLING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements 
!  
!  1. PTMC or MC Replica Exchange 
!  2. PTMC with Reservoir of Minima 
!  3. BSPT 
!
! MC trial move is all atom random cartesian perturbation. 
!
! Should be compiled with MPI. 
!
! Change history: 
!   04/09/201, ss2029:
!                 - MC step variable name changed from I to IMCSTEP 
!                 - implemented system time seeding for random number
!                 generator        
!                 - implemented keyword PTMCDUMPENER       
!                 - collecting and printing time taken per step 
!                 - several changes to output messages (e.g. end of
!                 equi, step size adjustment,.. )       
!                 - several explanatory comments and todo tasks!
!
! IMPORTANT VARIABLES
!
! XO, YO, ZO, VOLD = the old coordinates and energy. i.e. the saved markov state
! at the end of the previous step
!
! X, Y, Z, VNEW    = the new perterbed coordinates and energy.  If the step is
! rejected, these are set back to XO, YO, ZO, VOLD, so at the end of the
! iteration these are the state of the markov chain
!
! RECOUNT  =  the boolean which records whether the step is accepted or
! rejected.  If the step is rejected, RECOUNT==.TRUE., and the previous markov
! state is recounted.
!
! MINCOORDS, VMINNEW   =  the coordinates and energy of the perterbed coordinates (X,Y,Z) after a quenched
!
! VMINOLD    =  (?) the energy of the previous quench that wasn't rejecected. only used for MINDENSITY (I think)
!
! COORDS     = scratch space.  could be anything.  Used to do quenches
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE MODCHARMM
      USE COMMONS
      USE TETHERFUNC
      USE CLASS_OVERLAP
      IMPLICIT NONE

#ifdef MPI
      INCLUDE 'mpif.h'
      !
      ! NHIST=100 and NHISTE=1000 were hardcoded. They are used to allocate huge arrays later on. Changing to small numbers for now. 
      !     todo - read from data file
      !
      INTEGER, PARAMETER :: NHIST=2, NHISTE=2
      INTEGER :: IACCEPT(0:NPAR-1), MPIERR, J,K, TEMPUNIT
      INTEGER NHISTQ4(NHIST,0:NPAR-1), NDUMMY, NDUMMY2, NTOT, NH, IQE, IQ4, IQ6, J1, J2, J3, J4, JLOW,
     1        NHISTQ6(NHIST,0:NPAR-1), NHISTQE(NHISTE, 0:NPAR-1), N1,
     2        ITRAJ, ITRAJO,NEACCEPT, RNDSEED, NUPDATE,
     3        CONVERGED,LBFGS_ITERATIONS, BININDEX, MINIMANUMBER(HBINS,0:NPAR-1),
     4        NHISTALLQ(NHIST, NHIST, 0:NPAR-1), IBININDEX, IBININDEX2,
     5        NHISTGRAND(NHISTE, NHIST, NHIST, 0:NPAR-1),LBFGS_ITERATIONSO,PEVISITS2(NENRPER, HBINS, 0:NPAR-1),
     7        NOUTQBIN, NOUTPEBIN, LBFGS_ITERATIONSSAVE, HBINMAX, HBINMIN, GETUNIT, LUNIT
      DOUBLE PRECISION IMCSTEP, XOUT
      DOUBLE PRECISION QVISITS(HBINS, 0:NPAR-1), NACCEPTPT(0:NPAR-1)
      DOUBLE PRECISION PEVISITS(NENRPER,0:NPAR-1) ! instantaneous energy histogram from different replicas  

      DOUBLE PRECISION TEMPTRAJ(0:NPAR-1), BETA(0:NPAR-1), 
     1        EAV(0:NPAR-1), EAV2(0:NPAR-1), Q(3,NATOMS), Q4AV(0:NPAR-1), Q4AV2(0:NPAR-1), 
     2        Q6AV(0:NPAR-1), Q6AV2(0:NPAR-1),
     3        X(NATOMS), Y(NATOMS), Z(NATOMS), 
     4        CTE, T, VOLD, VNEW, POTEL, GRAD(3*NATOMS), Q4, Q6, RANDOM, DPRAND, Q4MAX, Q6MAX, 
     5        DQ4, DQ6, DHISTE, XO(NATOMS), YO(NATOMS), ZO(NATOMS), DDX, DDY, DDZ, DE, 
     6        W, WCOMP, WAC, E, CV, FQ4, FQ6, RMAX, DDXN, DDYN, DDZN, R2, 
     7        DUMMY, BINLABEL(HBINS), VNEWSAVE, DIHEORDERPARAM, SASAORDERPARAM, PEINT, HISTINT, 
     8        DIHEORDERPARAM_AV(0:NPAR-1), SASAORDERPARAM_AV(0:NPAR-1), MINCOORDS(3*NATOMS,NPAR),
     9        DIST,
     A        SR3, OPOTEL,
     B        DPRAND_UNIVERSAL

       DOUBLE PRECISION DIHEORDERPARAM_MIN, DIHEORDERPARAM_MAX, SASAORDERPARAM_MIN, SASAORDERPARAM_MAX, 
     1          DDIHE, DSASA, RGYR_MIN, RGYR_MAX, EINT_MIN, EINT_MAX, DRGYR, DEINT,
     2          RGYR_AV(0:NPAR-1), EINT_AV(0:NPAR-1), RGYR_AV2(0:NPAR-1), EINT_AV2(0:NPAR-1),
     3          ORDERPARAM1, ORDERPARAM2,
     5          VMINOLD, VMINNEW, CX, CY, CZ,
     6          DIHEORDERPARAM_AV2(0:NPAR-1), SASAORDERPARAM_AV2(0:NPAR-1), RGYR, EINT, VMINNEWSAVE

      CHARACTER (LEN=256)  FILENAME, FILENAME2,FILENAME3,FILENAME4,FILENAME5,FILENAME6,
     1                      FILENAME7,FILENAME8,FILENAME9, FILENAME10, FILENAME11, FILENAME12, 
     2                      FILENAME100,FILENAME101,FILENAME102, FILENAME103, FILENAME104, FILENAME105,
     3                      FILENAME106, FILENAMEDUMMY
      CHARACTER (LEN=80)  ISTR, SDUMMY
      CHARACTER (LEN=80)  CFNAME
      LOGICAL FITS, NEWENERGY, RECOUNT, OUTSIDE, JUMPT, EXCHANGEACCEPT
      LOGICAL BAD_CONFIG
      DOUBLE PRECISION RRX,RRY,RRZ,RRR

      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),ITIME

! <--ss2029 variables for writing markov energies from each replica 
      CHARACTER (LEN=256) FILENAME_ENER 
      INTEGER             FILENAME_ENER_LUNIT 
! ss2029--> 

! ab2111>:
      DOUBLE PRECISION VNEWab,POTELtest,deltalnJ

      DOUBLE PRECISION PREVSTEPS, PFNORM1, PFNORM2, NORM1, NORM2
      INTEGER OVERLAP_UNIT, OVERLAP_COUNT, KAPPA
      CHARACTER (LEN=200)  OVERLAP_FNAME
      DOUBLE PRECISION TIMESTART, OVERLAP_TIMETOT
      DOUBLE PRECISION TIME2, MCSTEPSTART, time1
      DOUBLE PRECISION timequench, time_tryexchange, time_potential, npot_call
      LOGICAL NORESET_TEMP
      INTEGER EXAB_COUNT, EXAB_ACC

!     DOUBLE PRECISION PTCAND(PTSTEPS/10,3*NATOMS), PTCANDE(PTSTEPS/10)
      DOUBLE PRECISION PEQMIN(NRESMIN), QV(NENRPER), ZT, VZERO, THISV, PEQV(NENRPER), XDUMMY
      INTEGER NGMIN, NCHOSEN
      DOUBLE PRECISION :: LASTEXUP(0:NPAR-1), LASTEXDOWN(0:NPAR-1)
      LOGICAL MY_ORDER_PARAMETERT

      COMMON /MYPOT/ POTEL

!
! Initialisation
!

      ! js850> MY_ORDER_PARAMETERT is a flag designed to protect
      ! the blocks of code that seem to deal with order parameters.
      ! They are system depended and shouldn't be run for all systems.
      ! Someone who knows more about this should protect them in a more
      ! inteligent way.
      IF (CALCQT .OR. ODIHET .OR. CHRMMT .OR. OSASAT .OR. ORGYT.or.OEINTT) THEN
         MY_ORDER_PARAMETERT = .TRUE.
      ELSE
         MY_ORDER_PARAMETERT = .FALSE.
      ENDIF

      EXAB_COUNT = 0
      EXAB_ACC = 0
      timequench = 0.
      time_tryexchange = 0.
      TIME_POTENTIAL = 0.
      npot_call = 0.

      !ab2111> jacobian factor set to zero 
      deltalnJ = 0.0

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NDUMMY,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      XOUT=0.0D0
      NOUTQBIN=0
      NOUTPEBIN=0
      LASTEXUP(0:NPAR-1)=-1.0D0
      LASTEXDOWN(0:NPAR-1)=-1.0D0

      WRITE(MYUNIT, '(A,I10,A,I10)') "bspt> This is processor", MYNODE+1, " of", NPAR
      WRITE(MYUNIT, '(A,I10)') 'bspt> Number of atoms', natoms
      IF (FIXSTEP(1)) FIXSTEP(2:NPAR)=.TRUE.
      IF (FIXTEMP(1)) FIXTEMP(2:NPAR)=.TRUE.
      IF (PERIODIC) THEN
         WRITE(MYUNIT, '(A,6G20.10)') 'bspt> Binary data', ntypea, epsab, epsbb, sigab, sigbb, cutoff
         WRITE(MYUNIT, '(A,3G20.10)') 'bspt> Box data', boxlx, boxly, boxlz
      ELSEIF(CHRMMT) THEN
         WRITE(MYUNIT, '(A)') 'bspt> CHARMM job'
         IF (FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1)) THEN
            WRITE(MYUNIT, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),' (fixed)'
         ELSE
            WRITE(MYUNIT, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),
     &                            ' (will be dynamically adjusted during equilibration)'
         ENDIF
      ELSEIF (.NOT.MODEL1T) THEN
         IF (PERCOLATET) THEN
            WRITE(MYUNIT, '(A,2G20.10)') 'bspt> Checking for a percolating system with threshold ',SQRT(PERCCUT)
         ELSE
            WRITE(MYUNIT, '(A,2G20.10)') 'bspt> Radius and Radius**2: ',SQRT(RADIUS),radius
         ENDIF
      ENDIF

      ITRAJ=MYNODE
      NEACCEPT=0 ! number of accepted PT exchanges
!     MINIMANUMBER=0 ! Number of distinct minima for quench bins in each replica. Not currently used.

      ! PEINT is the PE bin width.
      PEINT=(PTEMAX-PTEMIN)/NENRPER
      MAXEFALL=-1.0D100
      !HISTINT is the quenched energy bin width
      HISTINT=(HISTMAX-HISTMIN)/HBINS
      IF (MINDENSITYT.AND.(PTSTEPS*1.0D0.GT.0.0D0)) THEN
         WRITE(MYUNIT, '(A,2I8)') 'bspt> WARNING **** PTSTEPS reset to zero for density of minima run' 
         PTSTEPS=0
      ENDIF
!
! Limit for quench energy above which steps will be rejected.
!
!     HBINMAX=(MIN(BSPTQMAX,HISTMAX)-HISTMIN)/HISTINT+1
!     HBINMIN=(MAX(BSPTQMIN,HISTMIN)-HISTMIN)/HISTINT
      HBINMIN=1
      HBINMAX=HBINS
      WRITE(MYUNIT, '(A,2I8)') 'bspt> Minimum and maximum quench bin limits: ',HBINMIN,HBINMAX
      DO J1=1, HBINS
!        BINLABEL(J1)=HISTMIN + HISTINT*(J1-0.5D0) ! these energies point to the middle of the quench bin
         BINLABEL(J1)=HISTMIN + HISTINT*(J1-1.0D0) ! these energies point to the bottom of the quench bin
      ENDDO

!js850> initialize overlap stuff
      IF (OVERLAPK) THEN
        OVERLAP_COUNT = 0
        OVERLAP_UNIT=GETUNIT()
        WRITE(OVERLAP_FNAME,*) 'overlap.'
        WRITE(SDUMMY,'(I3)') MYNODE+1 
        OVERLAP_FNAME=TRIM(ADJUSTL(OVERLAP_FNAME)) // TRIM(ADJUSTL(SDUMMY))
        OPEN(UNIT=OVERLAP_unit,FILE=OVERLAP_FNAME, STATUS="unknown", form="formatted")
        WRITE(MYUNIT,'(A)') "bspt> overlap> IMCSTEP Qsum E Qsumquench Equench QA QB QexpA QexpB QexpAB" 
        WRITE(OVERLAP_UNIT,'(2A)') "#IMCSTEP, QAB, E, quenchQAB, Equench, QA, QB, QexpA, QexpB, QexpAB, quenchQexpA, ", 
     &         "quenchQexpB, quenchQexpAB Qexp2RA Qexp2RB Qexp2RAB quenchQexp2RA quenchQexp2RB quenchQexp2RAB"
      ENDIF

      KAPPA=3*NATOMS-6 !warning: only true for cluster systems

      ! set up the temperatures
      ! ab2111> RESERVOIR hard-coded temperature progression: USERES=0, with T_0 = T_1
      !IF (RESERVOIRT) THEN
      IF (.FALSE.) THEN
         CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-2)
      ELSE 
      ! default 
         CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
      ENDIF

      CTE=EXP(CTE)
      !ab2111> debug
!      CTE=EXP(LOG(0.6/(0.0125))/(NPAR-1))

      DO J1=0, NPAR-1
         !IF (RESERVOIRT) THEN
         IF (.FALSE.) THEN
            IF ((J1.EQ.0).OR.(J1.EQ.1)) THEN 
               TEMPTRAJ(J1) = PTTMIN
            ELSE
               TEMPTRAJ(J1)=PTTMIN*CTE**(J1-1)
            ENDIF
		   !default (non RESERVOIR) temperature progression
         ELSE
            TEMPTRAJ(J1)=PTTMIN*CTE**J1
         ENDIF
         BETA(J1)=1.0D0/TEMPTRAJ(J1)
      ENDDO
      IF (BETA_RES.EQ.(0.0D0)) BETA_RES = BETA(USERES)
     
      ! RESERVOIR initialisation 
      IF (RESERVOIRT) THEN
      ! ab2111> debug: initialize coordinates at well 1 / 2
         DO J1=1,3*NATOMS
            IF (MYNODE.EQ.1) THEN 
               COORDS(J1,MYNODE+1)=RESPOINTS(J1,1) 
            ELSEIF (MYNODE.EQ.2) THEN
               COORDS(J1,MYNODE+1)=RESPOINTS(J1,2) 
            ENDIF
         ENDDO
      ENDIF

      WRITE(MYUNIT, '(A, I4)') 'bspt> Coordinates and temperatures initialized for replica = ', MYNODE+1

      ! each node must create it's own directory because they may be run on
      ! local scratch on different nodes on the cluster
      WRITE(SDUMMY,'(A,I6,A)') 'mkdir ',mynode+1,' &> /dev/null'
      CALL SYSTEM(SDUMMY)

!
! MPI synchronise to make sure that the directories have been created.
! Needed so that dumpenergy files can be opened in the directories.
! ss2029. 
!
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)

!
! ss2029> Open dumpenergy file for each replica 
!
      IF ( PTMCDUMPENERT ) THEN 
          WRITE (ISTR, '(I2,A1)') MYNODE+1
          ISTR=TRIM(ADJUSTL(ISTR)) // '/dumpenergy'
          WRITE(FILENAME_ENER, *) ISTR 
          ! print *,'MYNODE = ', MYNODE+1, ' file name = ', FILENAME_ENER 
          ! open file 
          FILENAME_ENER_LUNIT=GETUNIT()
          OPEN(UNIT=FILENAME_ENER_LUNIT,FILE=FILENAME_ENER, STATUS="unknown", form="formatted")
          WRITE (MYUNIT,'(A)') 'bspt> Instantaneous potential energies file = '//TRIM(ADJUSTL(FILENAME_ENER))
          WRITE (MYUNIT,'(A,I10)') 'bspt> Instantaneous potential energy output frequency = ', PTMCDUMPENERFRQ 
      ENDIF 

! For restart we need to get the current configuration, its pe, the pe of the minimum it quenched to,
! if applicable, the number of steps already done, the maximum step size, and the Visits and Visits2
! histograms. If we dump using BSPTDUMPFRQ then we can restore from the last such file. We can work
! out what the last dump was once we know how many steps have been done!
!
      IF (BSPTRESTART) THEN
         !
         !read file [mynode+1]/bsptrestart
         !
         CALL BSPT_READ_BSPTRESTART( COORDS(1:3*NATOMS,MYNODE+1),
     &      PREVSTEPS, VOLD, VMINOLD, NACCEPTPT, NEACCEPT, NTOT, NOUTQBIN,
     &      NOUTPEBIN )

         !
         !read file Visits.his file
         !
         IF (PREVSTEPS.GT.NEQUIL) THEN
            CALL BSPT_READ_VISTS_HIS( PREVSTEPS, PEVISITS, QVISITS, T )
            TEMPTRAJ(MYNODE) = T
            BETA(MYNODE) = 1.0D0/T
         ENDIF
         IF (BSPT.AND.(PREVSTEPS.GT.NEQUIL+PTSTEPS).AND.(.NOT.MINDENSITYT)) THEN
            call bspt_read_vists_his2( pevisits2 )
         ENDIF
         CONVERGED=1
      ELSE
!
! This is the block executed if we are not restarting.
!
12       CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
         VOLD=POTEL
         WRITE(MYUNIT,'(2(A,G20.10))') 'bspt> Initial configuration energy             = ',VOLD
!        WRITE(MYUNIT,'(2(A,G20.10))') 'bspt> NOT QUENCHING initial config'
         CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1))
         WRITE(MYUNIT,'(A,G20.10)')    'bspt> Initial configuration quenched to energy = ',POTEL
!
!  If we start with an invalid configuration we may never recover!
!
         IF (.NOT.(CHRMMT.OR.MODEL1T.OR.PERIODIC)) THEN
             IF (PERCOLATET) THEN
               CALL PERC(COORDS(1:3*NATOMS,MYNODE+1),NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)
               IF (.NOT.PERCT) THEN
                  WRITE(MYUNIT,'(A,I6,A)') 'bspt> ERROR *** Initial coordinates do not define a percolating system'
                  STOP
               ENDIF
            ELSE
               DO K=1,NATOMS
                  DIST=COORDS(3*(K-1)+1,MYNODE+1)**2+COORDS(3*(K-1)+2,MYNODE+1)**2+COORDS(3*(K-1)+3,MYNODE+1)**2
                  IF (DIST.GT.RADIUS) THEN
                     IF (MOD(IMCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,I6,A)') 'bspt> Atom ',K,
     &                       ' outside container, reseed and try again'
                     SR3=DSQRT(3.0D0)
                     DO J1=1,NATOMS
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                        COORDS(3*(K-1)+1,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                        COORDS(3*(K-1)+2,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                        COORDS(3*(K-1)+3,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                     ENDDO
                     GOTO 12
                  ENDIF
               ENDDO 
            ENDIF
         ENDIF
      ENDIF

      ! write temperatures to file
      ! Bookkeeping is performed by node 0. 
      IF (MYNODE.EQ.0) THEN
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT, FILE='temperatures',STATUS='UNKNOWN')
         WRITE(LUNIT,'(G20.10)') TEMPTRAJ(0:NPAR-1)
         CLOSE(LUNIT)
      ENDIF

      LBFGS_ITERATIONS=0
      LBFGS_ITERATIONSO=0 ! so that we don't use quench DoS statistics from arbitrary high energy
      !VOLD=POTEL  !js850> POTEL is not defined yet for bsptrestart
      VNEW=VOLD
      DO J1=1,NATOMS
         X(J1)=COORDS(3*(J1-1)+1,MYNODE+1)
         Y(J1)=COORDS(3*(J1-1)+2,MYNODE+1)
         Z(J1)=COORDS(3*(J1-1)+3,MYNODE+1)
      ENDDO

      IF (MY_ORDER_PARAMETERT) THEN
         DO J1=1,NATOMS
            Q(1,J1)=X(J1)
            Q(2,J1)=Y(J1)
            Q(3,J1)=Z(J1)
         ENDDO

         IF (PERIODIC) THEN
            CALL QORDER_BLJ(Q,Q4,Q6)
         ELSE IF (CHRMMT) THEN
              IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
              IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
              IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.FALSE.) 
              IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
         ELSE
            CALL QORDER_LJ(Q,Q4,Q6)
         ENDIF
      ENDIF

      WRITE(MYUNIT, '(A,2G20.10)') 'bspt> Temperature range             = ', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A,G20.10)')  'bspt> This temperature trajectory   = ', TEMPTRAJ(MYNODE)
      IF (.NOT.BSPTRESTART) VMINOLD=VOLD
      IF ((.NOT.BSPTRESTART).AND.(NEQUIL+PTSTEPS.GT.0)) VMINOLD=0.0D0 
      WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of equilibration steps = ',NEQUIL
      WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of PT steps            = ',PTSTEPS
      WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of BSPT steps          = ',NQUENCH
      WRITE (MYUNIT,'(A,G20.10)')   'bspt> Replica Exchange Interval    = ',EXCHINT 
      WRITE (MYUNIT,'(A,G20.10)')   'bspt> Replica Exchange Prob        = ',EXCHPROB

      IF (MINDENSITYT) WRITE (MYUNIT,'(A,F20.1)') 'bspt> Accept/reject based on quench energy'
      IF (PERIODIC) WRITE(MYUNIT, '(A, 2G20.10)') 'bspt> Starting Q4, Q6=', Q4, Q6
      IF (CHRMMT.and.MY_ORDER_PARAMETERT) THEN
        IF (OSASAT)  WRITE(MYUNIT, '(A, 2G20.10)') 'bspt> Starting Met-enk order params=', 
     &                                         DIHEORDERPARAM, SASAORDERPARAM 
        IF (ORGYT)  WRITE(MYUNIT, '(A, G20.10)') 'bspt> Starting Rgyr order params=', 
     &                                          RGYR 
        IF (OEINTT) WRITE(MYUNIT, '(A, G20.10)') 'bspt> Starting Eint order params=',
     &                                          EINT
      ENDIF
!
! Initialisation complete
!
!   ss2029> implementing system time seed as done in mc.F  
!
      IF (RANSEEDT) THEN
        CONTINUE
      ELSE IF (RANDOMSEEDT) THEN
        CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
        ITIME= values(6)*60 +  values(7)
        CALL SDPRND(ITIME)
        CALL SDPRND_UNIVERSAL(ITIME+NPAR)
        WRITE(MYUNIT, '(A)') 'bspt> Using system time to seed random number generator' 
        WRITE(MYUNIT, '(A)') 'bspt> and system time + NPAR to seed universal random number generator'
      ELSE
        WRITE(MYUNIT, '(A)') 'bspt> Using hardcoded value (=2002) as random number seed' 
        WRITE(MYUNIT, '(A)') 'bspt> and 2002 + NPAR as universal random number seed' 
        RNDSEED=2002+MYNODE
        CALL SDPRND(RNDSEED)
        CALL SDPRND_UNIVERSAL(2002+NPAR)
      ENDIF
      RANDOM=DPRAND()
      WRITE(MYUNIT, '(A, G20.10)') 'bspt> Starting local random number     = ', RANDOM
      RANDOM=DPRAND_UNIVERSAL()
      WRITE(MYUNIT, '(A, G20.10)') 'bspt> Starting universal random number = ', RANDOM


      IF (FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1)) THEN
         WRITE(MYUNIT, '(A,G20.10,A)') 'bspt> Starting stepsize            = ', STEP(MYNODE+1),' (fixed)'
      ELSE
         WRITE(MYUNIT, '(A,G20.10)') 'bspt> Starting stepsize (adjusted) = ', STEP(MYNODE+1)
      ENDIF
      CALL FLUSH(MYUNIT)

      IF (.NOT.BSPTRESTART) then
         NTOT=0
         NACCEPTPT(MYNODE)=0.0D0
         QVISITS(:,MYNODE)=0.0D0
         PEVISITS(:,MYNODE)=0.0D0
         PEVISITS2(:,:,MYNODE)=0
      ENDIF
      IACCEPT(MYNODE)=0
      !js850> these are not used anymore
      !EAV(MYNODE)=0.
      !EAV2(MYNODE)=0.

      IF (MY_ORDER_PARAMETERT) THEN
         !js850> todo: this should go to a separate subroutine, or subroutines
         Q4AV(MYNODE)=0.
         Q6AV(MYNODE)=0.
         Q4AV2(MYNODE)=0.
         Q6AV2(MYNODE)=0.
         DIHEORDERPARAM_AV(MYNODE)=0.
         SASAORDERPARAM_AV(MYNODE)=0.
         DIHEORDERPARAM_AV2(MYNODE)=0.
         SASAORDERPARAM_AV2(MYNODE)=0.
         RGYR_AV(MYNODE)=0.
         EINT_AV(MYNODE)=0.
         RGYR_AV2(MYNODE)=0.
         EINT_AV2(MYNODE)=0.
         DO J1=1,NHIST
            NHISTQ4(J1,MYNODE)=0
            NHISTQ6(J1,MYNODE)=0
         ENDDO
         DO J1=1,NHISTE
            NHISTQE(J1,MYNODE)=0
         ENDDO
         NHISTALLQ(:,:,MYNODE)=0
         NHISTGRAND(:,:,:,MYNODE)=0
         !DISTANCE(:,MYNODE)=0.0D0
         !DISTANCE_AV(:,MYNODE)=0.0D0
         !DELTA_EAV(MYNODE)=0.0D0

         Q4MAX=0.1
         Q6MAX=0.5 ! Will have to change as is rather system specific
         DIHEORDERPARAM_MIN=0.4D0
         DIHEORDERPARAM_MAX=1.0D0
         SASAORDERPARAM_MIN=300.0D0
         SASAORDERPARAM_MAX=800.0D0
         RGYR_MIN=4.0D0
         RGYR_MAX=12.0D0
         EINT_MIN=-35.0
         EINT_MAX=10.0
         DDIHE=(DIHEORDERPARAM_MAX-DIHEORDERPARAM_MIN)/(NHIST-1)
         DSASA=(SASAORDERPARAM_MAX-SASAORDERPARAM_MIN)/(NHIST-1)
         DRGYR=(RGYR_MAX-RGYR_MIN)/(NHIST-1)
         DEINT=(EINT_MAX-EINT_MIN)/(NHIST-1)
         DQ4=Q4MAX/(NHIST-1)
         DQ6=Q6MAX/(NHIST-1)
         DHISTE=(PTEMAX-PTEMIN)/(NHISTE-1)
      ENDIF !closes if (MY_ORDER_PARAMETERT)

!
! NUPDATE specifies the interval for dynamically altering the maximum step size.
! Only used if step size isn't fixed.
!
! ss2029> NUPDATE is hardcoded. read from a keyword. todo.  
!
      NUPDATE=100
      IF (.NOT.(FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1)))  THEN 
          WRITE(MYUNIT, '(A,I6,A)') 'bspt> Step size will be adjusted every ',NUPDATE ,' MC steps'
      ENDIF

      IF (DEBUGss2029) THEN
          print *,'bspt> ss2029> MYNODE = ', MYNODE+1, ' just before the main MC loop.' 
      ENDIF

      WRITE (MYUNIT,'(A)') 'bspt> Entering main MC loop ... ' 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop over steps. We have to use double precision for the number of steps, because it
! may exceed integer capacity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMCSTEP=0.0D0    ! = MC step number 
      IF (BSPTRESTART) IMCSTEP=PREVSTEPS
      MCSTEPSTART = IMCSTEP
      CALL MYCPU_TIME(TIMESTART)
      DO 
         IMCSTEP=IMCSTEP+1.0D0
         IF (IMCSTEP.GT.NEQUIL+PTSTEPS+NQUENCH) EXIT
         RECOUNT=.FALSE.
         DO K=1, NATOMS
            XO(K)=X(K)
            YO(K)=Y(K)
            ZO(K)=Z(K)
            COORDSO(3*(K-1)+1,MYNODE+1)=XO(K)
            COORDSO(3*(K-1)+2,MYNODE+1)=YO(K)
            COORDSO(3*(K-1)+3,MYNODE+1)=ZO(K)
         ENDDO
!
! Should probably worry about centre of coordinates versus centre of mass
! for heteroatomic systems.
!
         IF (CENT) THEN
            CX=0.0D0; CY=0.0D0; CZ=0.0D0
            DO K=1,NATOMS
               CX=CX+X(K); CY=CY+Y(K); CZ=CZ+Z(K)
            ENDDO
            CX=CX/NATOMS; CY=CY/NATOMS; CZ=CZ/NATOMS
            DO K=1,NATOMS
               X(K)=X(K)-CX; Y(K)=Y(K)-CY; Z(K)=Z(K)-CZ
            ENDDO
         ENDIF  
         IF ((DEBUG.OR.CHECKMARKOVT).AND.(IMCSTEP.GT.1.0D0)  
     &        ) THEN
            CALL POTENTIAL(COORDSO(:,MYNODE+1),GRAD,OPOTEL,.FALSE.,.FALSE.)
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,VNEWab,.FALSE.,.FALSE.)
            IF (ABS(OPOTEL-VNEW).GT.ECONV) THEN
               IF(RESERVOIRT)  THEN
                  WRITE(MYUNIT,'(4(A,G20.10))') 'bspt> possible ERROR - energy for coordinates in COORDSO=',
     &                     OPOTEL,' but Markov energy=',VNEW,'IMCSTEP=',IMCSTEP, 'VNEWab',VNEWab
               ELSE 
                  WRITE(MYUNIT,'(3(A,G20.10))') 'bspt> possible ERROR - energy for coordinates in COORDSO=',
     &                     OPOTEL,' but Markov energy=',VNEW,'IMCSTEP=',IMCSTEP
               STOP
               ENDIF
            ENDIF
            IF (PERCOLATET) THEN
               CALL PERC(COORDSO(1:3*NATOMS,MYNODE+1),NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)
               IF (.NOT.PERCT) THEN
                  IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) WRITE(MYUNIT,'(A)') 
     &              'bspt> Markov structure is not a percolating network'
               ENDIF
            ENDIF
         ENDIF

         JUMPT=.FALSE.
!
! If a move is rejected, the result for the current configuration
! should be recounted. This allows for the unsymmetrical forward
! and backward move probabilities for an atom near the surface of
! the constraining sphere.
!             

         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, VOLD, .TRUE., .FALSE.)

         !
         ! take a monte trial carlo step
         !
         !ab2111>
         IF (AMBERT.EQV..TRUE.) THEN
            CALL  BSPT_TAKESTEP_amber(X, Y, Z, IMCSTEP, 
     &               JUMPT, 
     &               PEINT, PEQV, VNEW, 
     &               EXAB_COUNT,BETA,deltalnJ)
         ELSE 
            CALL  BSPT_TAKESTEP(X, Y, Z, IMCSTEP, 
     &            JUMPT, 
     &            PEINT, PEQV, VNEW, 
     &            EXAB_COUNT)

         ENDIF
         ! make sure no frozen atoms have moved
         IF (FREEZE) THEN
            DO J1=1,NFREEZE
               J2=FROZENLIST(J1)
               X(J2)=XO(J2)
               Y(J2)=YO(J2)
               Z(J2)=ZO(J2)
            ENDDO
         ENDIF

         ! copying post-move coordinates 
         DO K=1,NATOMS
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO

         !
         ! do the configuration checks on the new coordinates.
         ! the result will be saved in BAD_CONFIG
         ! if any of the checks fail then BAD_CONFIG will be TRUE
         ! if they all pass, then BAD_CONFIG will be FALSE
         !
         CALL BSPT_CHECK_CONFIGURATION(X, Y, Z, BAD_CONFIG,
     &               IMCSTEP, XOUT, OUTSIDE)
         RECOUNT = BAD_CONFIG

!
! At this point all we have done is take a step, energies have not been computed yet and 
! also Metropolis check has not been done. 
!
! The perturbed coordinates are in both COORDS and X, Y, Z. 
! The old coordinates are in XO, YO, ZO. 
!
! For node 0 and reservoir, we have instead chosen a PE bin.
!
! New and old quench energies will be in        VMINNEW VMINOLD
! New and old instantaneous energies will be in    VNEW    VOLD
!
! ss2029> Compute energy of perturbed coordinates and apply Metropolis check.   
!         If the step is to be rejected, RECOUNT is set to TRUE 
!
         IF (RECOUNT) THEN
            VNEW=0.0D0
         ELSE 
            CALL MYCPU_TIME(TIME1)
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
            CALL MYCPU_TIME(TIME2)
            npot_call = npot_call + 1
            TIME_POTENTIAL = TIME_POTENTIAL + TIME2 - TIME1
            VNEW=POTEL
            IF (.NOT.MINDENSITYT) THEN
            ! ss2029> this is where Metropolis check on instantaneous energy is done. 
               WCOMP=(VNEW-VOLD)*BETA(MYNODE) ! use difference in instantaneous energies
               ! ab2111 > include jacobian factor
               WCOMP=WCOMP - deltalnJ
               W=MIN(1.0D0,EXP(-WCOMP))
               RANDOM=DPRAND()
               IF (RANDOM.GT.W) RECOUNT=.TRUE. ! RECOUNT is initialised to .FALSE. at the top of the loop
               IF (JUMPT)  WRITE(MYUNIT, '(A,4G20.10,L5)') 'bspt> VOLD,VNEW,W,RANDOM,RECOUNT=',VOLD,VNEW,W,RANDOM,RECOUNT
               IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) THEN
                  WRITE(MYUNIT, '(A,4G20.10,L5,F10.1)') 'bspt> VOLD,VNEW,W,RANDOM,RECOUNT=',VOLD,VNEW,W,RANDOM,RECOUNT,IMCSTEP
               ENDIF
            ENDIF
         ENDIF

! histogram bounds check on VNEW.  Reject the step if the energy is out of
! bounds of the histogram
         IF (.NOT.RECOUNT) THEN
            IBININDEX=INT((VNEW-PTEMIN)/PEINT)+1
            IF ((IBININDEX.GT.NENRPER).OR.(IBININDEX.LT.1)) THEN
               IF (.NOT.MINDENSITYT) RECOUNT=.TRUE.
!              RECOUNT=.TRUE.
               IF (.NOT.OUTSIDE) NOUTPEBIN=NOUTPEBIN+1
            ENDIF
         ENDIF

!
! Quenching part if required.
!
         IF (BSPT) THEN
            CALL MYCPU_TIME(TIME1)
            CALL BSPT_DO_QUENCHING(X, Y, Z, XO, YO, ZO, MINCOORDS, IMCSTEP,
     &         VMINNEW, VMINOLD, LBFGS_ITERATIONS, RECOUNT, HISTINT, NOUTQBIN,
     &         BETA, OUTSIDE )
            CALL MYCPU_TIME(TIME2)
            timequench = timequench + (time2 - time1)
         ENDIF

         VNEWSAVE=VNEW  ! this value is saved so it can be printed if the step is rejected
         IF (BSPT) THEN
            VMINNEWSAVE=VMINNEW
            LBFGS_ITERATIONSSAVE=LBFGS_ITERATIONS
         ENDIF

!
! ss2029> finally! this is where the move is accepted/rejected 
!
         IF (RECOUNT) THEN ! reject move
            ! 
            ! ss2029> overwriting post-move coordinates X,Y,Z why over writing post-move coordinates X,Y,Z?
            !         needed because XO, YO and ZO are initialized to X,Y,Z at the top of the loop. 
            !         Post-move coordinates are lost at this point. If needed, add another variable for post-move coords.
            !          
            DO K=1, NATOMS
               X(K)=XO(K)
               Y(K)=YO(K)
               Z(K)=ZO(K)
            ENDDO
            DO K=1,NATOMS
              COORDS(3*(K-1)+1,MYNODE+1)=X(K)
              COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
              COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
            ENDDO
            VNEW=VOLD
            VMINNEW=VMINOLD
            LBFGS_ITERATIONS=LBFGS_ITERATIONSO
         ELSE ! accept move
            ! since move is accepted, no change to COORDS(..,MYNODE+1)
            NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1.0D0
            IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
            IF ( BINARY_EXAB .AND. (MOD(IMCSTEP-1.0D0,BINARY_EXAB_FRQ*1.0D0).EQ.0.0D0) ) EXAB_ACC = EXAB_ACC + 1
         ENDIF ! closes IF (RECOUNT)
   
         
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
         

!
! ss2029> Local MC move has been accepted/rejected. Stats can be collected now.
! No replica exchange  so far. 
!
         SDUMMY='ACC'
         IF (RECOUNT) SDUMMY='REJ'
         IF (MOD(IMCSTEP-1.0D0*1,1.0D0*PRTFRQ).EQ.1.0D0*0) THEN
            IF ( PTMC .AND. .NOT. BSPT ) THEN !js850> no need to print energies of minima
               WRITE(MYUNIT, '(F15.1,A,F20.10,A,F20.10)') 
     &            IMCSTEP,' Vn=', VNEWSAVE,' Vo=',VOLD
            ELSE
               WRITE(MYUNIT, '(F15.1,A,G16.6,A,G16.6,A,G16.6,A,G16.6,A,I6,1X,A)') 
     &            IMCSTEP,' En= ', VMINNEWSAVE, ' Eo= ',VMINOLD, ' Vn=', VNEWSAVE,' Vo=',VOLD,' iter ',
     &            LBFGS_ITERATIONSSAVE,TRIM(ADJUSTL(SDUMMY))
            ENDIF
            !write the time per step
            CALL MYCPU_TIME(TIME2)
            WRITE(MYUNIT, '(4(A,G15.5))') "bspt> time (sec) ", time2-timestart, 
     &              " per step ", (time2-timestart)/(IMCSTEP-mcstepstart), 
     &              " in tryexchange ", time_tryexchange, " computing energy ", time_potential
         ENDIF
         IF (DEBUG) THEN
!           IF (.NOT.((MYNODE.LE.USERES).AND.RESERVOIRT)) THEN
!              CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
!              WRITE(MYUNIT,'(A,G20.10)') 'bspt> potential energy after acc/rej=',POTEL
!           ELSE
!              WRITE(MYUNIT,'(A,G20.10)') 'bspt> potential energy after acc/rej=',VNEW
!           ENDIF
            IF (MOD(IMCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,2F15.1)') 
     &         'bspt> Last exchanges to replicas above and below were at steps: ',LASTEXUP(MYNODE),LASTEXDOWN(MYNODE)
         ENDIF
         IF (DEBUG) CALL FLUSH(MYUNIT)
!
!  At this point the quench and instantaneous energies for the current
!  configuration in the Markov chain are VMINNEW and VNEW.
!
!  QUENCHFRQ must be initialised to one for PTMC to avoid division by zero!
!
         IF (BSPT.and.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.1.0D0*0)) BININDEX=INT((VMINNEW-HISTMIN)/HISTINT)+1
         IBININDEX=INT((VNEW-PTEMIN)/PEINT)+1
!
!  Must not accumulate statistics until we have equilibrated for NEQUIL steps.
!
         IF (IMCSTEP.GT.NEQUIL) THEN
            ! BSPT stuff 
            IF (BSPT.AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(IMCSTEP.GT.NEQUIL+PTSTEPS)) THEN 
               IF (BINSTRUCTURES.AND.(MOD(IMCSTEP, 1.0D0*SAVENTH).EQ.1.0D0*0).AND.(.NOT.RECOUNT)) THEN 
                  CALL SAVEBINSTRUCTURESMPI(POTEL,MINCOORDS(:,MYNODE+1),BININDEX,.TRUE.,MYNODE,NEWENERGY,MYUNIT)
!                 IF (NEWENERGY) THEN 
!                     MINIMANUMBER(BININDEX, MYNODE)=MINIMANUMBER(BININDEX, MYNODE)+1
!                 ENDIF
               ENDIF
            ENDIF 
            ! closes IF (BSPT.AND...) 

!           IF (BSPT.AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(IMCSTEP.GT.NEQUIL+PTSTEPS)) THEN
!              IF ((BININDEX.GT.HBINS).OR.(BININDEX.LT.1)) GOTO 888
!           ENDIF
!           IF ((IBININDEX.GT.NENRPER).OR.(IBININDEX.LT.1)) GOTO 888

            ! Must protect histograms from out-of-bounds error if we are not recounting.
            ! 
            IF ((IBININDEX.LE.NENRPER).AND.(IBININDEX.GE.1)) THEN
               PEVISITS(IBININDEX,MYNODE)=PEVISITS(IBININDEX,MYNODE)+1.0D0
            ENDIF

            ! BSPT stuff 
            IF (BSPT.AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(IMCSTEP.GT.NEQUIL+PTSTEPS)) THEN ! otherwise BININDEX could be out of range
               IF ((BININDEX.LE.HBINS).AND.(BININDEX.GE.1)) THEN
                  QVISITS(BININDEX,MYNODE)=QVISITS(BININDEX,MYNODE)+1.0D0
                  IF ((IBININDEX.LE.NENRPER).AND.(IBININDEX.GE.1)) THEN
                     PEVISITS2(IBININDEX,BININDEX,MYNODE)=PEVISITS2(IBININDEX,BININDEX,MYNODE)+1
                  ENDIF
               ENDIF
            ENDIF


! at this point the state points are saved in the following variables (from what
! I can tell -- js850).  COORDS is used as scratch space
!  ---The saved markov state unquenched:
!              XO, YO, ZO      energy=VOLD
!  ---The saved markov state quenched: (or unquenched if PTMC)
!              ?????           energy=VMINOLD
!  ---After the step, but unquenched:
!              not saved        energy=VNEWSAVE?
!  ---After the quench:
!              MINCOORDS       energy=VMINNEWSAVE?
!  ---The current markov state unquenched: i.e. after the step was accepted or rejected.
!              X, Y, Z          energy=VNEW
!  ---The current markov state quenched (or unquenched if PTMC): i.e. after the step was accepted or rejected.
!              ?????????        energy=VMINNEW
!

! 
! ss2029> another summary, reconciled with js850's summary above. 
!
! Just considering PTMC (no quenching), infact, this is just
! MC, PT has not been done yet. That is later around the TRYEXCHANGE call
! So just for PT, important variables are: 
!  
!     --- pre-move 
!            coord = XO,YO,ZO   energy = VOLD  
!     --- post-move and post-Metropolis  
!            coord = X, Y, Z    energy = VNEW (if accepted) 
!            coord and energy are lost if rejected  
!     --- Outcome of Metropolis and other checks 
!            RECOUNT = FALSE if accepted 
!                    = TRUE  if rejected
!     --- Markov state post Metropolis 
!            coord = COORDS      energy = VNEW              
!

!
! Should take out the distance if we aren't going to use it DJW
!
            !IF (BSPT.AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.1.0D0*0).AND.(IMCSTEP.GT.NEQUIL+PTSTEPS)) THEN
               !DISTANCEOLD=CALCULATEDDISTANCE(COORDS(:,MYNODE+1), COORDSO(:,MYNODE+1))
               !DISTANCE(BININDEX,MYNODE)=DISTANCE(BININDEX,MYNODE)+DISTANCEOLD
            !ENDIF

!
! ss2029> write Markov energy to repNum/dumpenergy 
!
           IF ( PTMCDUMPENERT ) THEN
             IF ( (IMCSTEP.GE.NEQUIL) .AND. (MOD(IMCSTEP-1.0D0,PTMCDUMPENERFRQ*1.0D0).EQ.0.0D0) ) THEN
             ! print *, 'MCstep = ', IMCSTEP , '   Markov Ener = ', VNEW 
               WRITE(FILENAME_ENER_LUNIT,'(G20.10,A,G20.06)') IMCSTEP,'      ',  VNEW 
             ENDIF 
           ENDIF 

            IF ( PTMCDUMPSTRUCT .AND. (MOD(IMCSTEP-1.0D0,PTMCDS_FRQ*1.0D0).EQ.0.0D0) ) THEN
               ! js850> PTMCDUMPSTRUCT
               CALL PTMC_DUMPSTRUCT_DUMP(COORDS(1:3*NATOMS, MYNODE+1), VNEW, IMCSTEP)
            ENDIF
            IF ( OVERLAPK .AND. (MOD(IMCSTEP-1.0D0,OVERLAP_FRQ*1.0D0).EQ.0.0D0) ) THEN
               ! js850> OVERLAP keyword
               CALL PTMC_OVERLAP_DUMP(COORDS(1:3*NATOMS, MYNODE+1), OVERLAP_UNIT, MYUNIT, OVERLAP_COUNT,
     &            time1, time2, VNEW, IMCSTEP)
               overlap_timetot = overlap_timetot + time1 
               timequench = timequench + time2 
            ENDIF

         ENDIF
!
! If we are not calculating distances then the following lines may be unnecessary.
!
888      DO K=1,NATOMS ! COORDS have to be resaved because we have to use noreset to be able to calc distance 
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO

!
! ss2029> Found it! This is where the step size is adjusted. All adjustments
! are make during equilibration (IMCSTEP<NEQUIL). Updates are
! made every NUPDATE steps based on acceptance ratio during previous
! NUPDATE steps.
!
         IF ((IMCSTEP.LE.NEQUIL).AND.(MOD(IMCSTEP,1.0D0*NUPDATE).EQ.1.0D0*0)) THEN ! update MC step size if not fixed
            WAC=1.0*IACCEPT(MYNODE)/NUPDATE
            IF (.NOT.(FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1))) THEN
               IF (WAC.LT.ACCRAT(MYNODE+1)-0.1D0) THEN
                  STEP(MYNODE+1)=STEP(MYNODE+1)*0.9D0
               ENDIF
               IF (WAC.GT.ACCRAT(MYNODE+1)+0.1D0) THEN
                  STEP(MYNODE+1)=STEP(MYNODE+1)*1.1D0
               ENDIF
            ENDIF
            IACCEPT(MYNODE)=0
            ! ss2029> commenting next line, due to mod(..,PRTFRQ) it was
            !           being printed only if PRTFRQ was a multiple of
            !           NUPDATE 
            ! IF (MOD(IMCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,G20.10)') 'bspt> maximum step size is now ',STEP(MYNODE+1)
            
            ! 
            ! ss2029> Now printing everytime the step size is changed.
            !   The output gives an idea of the length of equilibration
            !   required to converge the step size
            ! 
            WRITE(MYUNIT, '(A,G20.10,A,G20.10,A,I6,A)') 'bspt> adjusting step-size> current step size = ', STEP(MYNODE+1) ,
     &                        ' acceptance ratio = ', WAC ,
     &                        ' over ', NUPDATE, ' steps'
         ENDIF ! step size update 

         !
         ! ss2029> print end of equilibration 
         !
         IF (IMCSTEP.EQ.NEQUIL) THEN 
               WRITE(MYUNIT, '(A)') 'bspt> ---------- Equilibration done ' 
               WRITE(MYUNIT, '(A,I6,A,F20.10,A,G20.10,A,G20.10)') 'bspt> Replica = ', MYNODE+1 , 
     &                           ' Temperature = ', TEMPTRAJ(MYNODE) ,
     &                           ' MCStep = ', IMCSTEP ,
     &                           ' MarkovEner = ', VNEW 
               WRITE(MYUNIT, '(A,G20.10,A,G20.10,A,I6,A)') 'bspt> Final step size = ', STEP(MYNODE+1) ,
     &                        '  corresponding to acceptance ratio = ', WAC , ! WAC calculated in step-adjustment block above 
     &                        ' over previous ', NUPDATE, ' steps'
               WRITE(MYUNIT, '(A,G20.10)') 'bspt>   compare with target acceptance ratio = ', ACCRAT(MYNODE+1)
               WRITE(MYUNIT, '(A)') 'bspt> ---------- Starting production run ' 
         ENDIF 

         !!!!!!!!!!!!!!!!!!! REPLICA EXCHANGE begin !!!!!!!!!!!!!!!!!!
         
         VOLD=VNEW         ! saving current Markov state 
         VMINOLD=VMINNEW   !

         ! IF ( IMCSTEP.GT.NEQUIL) THEN 
         LBFGS_ITERATIONSO=LBFGS_ITERATIONS
         E=VNEW
!
!  Replica exchange part is done in TRYEXCHANGE.
!
         IF (EXCHPROB.GT.0.0D0) THEN
            CALL MYCPU_TIME(TIME1)
            CALL TRYEXCHANGE(E,X,Y,Z,XO,YO,ZO,VOLD,EXCHANGEACCEPT,JLOW,
     &                  VNEW,GRAD,VNEWSAVE,VMINOLD,VMINNEW,BETA,ITRAJ,ITRAJO,NTOT, 
     &                  LBFGS_ITERATIONS,NEACCEPT,LBFGS_ITERATIONSO,QV,XDUMMY,PEINT,NCHOSEN,LASTEXDOWN,IMCSTEP)
            CALL MYCPU_TIME(TIME2)
            time_tryexchange = time_tryexchange + (time2 - time1)
            IF (EXCHANGEACCEPT) then
               IF (MYNODE.EQ.JLOW)    LASTEXUP(MYNODE)=IMCSTEP
               IF (MYNODE.EQ.JLOW+1)  LASTEXDOWN(MYNODE)=IMCSTEP
            ENDIF
         ENDIF
         !  js850>At this point VOLD, VNEW  are the markov energy (perhaps changed)
         !                      X           is the markov coords
         !                      E, VNEWSAVE, XO are not determined.  They should not be used after this
         E = VNEW !js850> fix a bug here
         !ENDIF

         !!!!!!!!!!!!!!!!!!! REPLICA EXCHANGE end !!!!!!!!!!!!!!!!!!

         ! 
         ! Stats (EAV, etc) for each replica are collected within the next IF block. 
         !  Logical flow (MC->stats->PT) might be better if this block is moved up to before replica
         !  exchange is done. 
         !  todo - clean up the Q-stuff by moving it behind CALCQT flag 
         !

         IF (IMCSTEP.GT.NEQUIL) THEN 
            !js850> these are not used anymore
            !EAV(MYNODE)=EAV(MYNODE)+E
            !EAV2(MYNODE)=EAV2(MYNODE)+E**2

            IF (MY_ORDER_PARAMETERT) THEN
               !js850> todo: this should go to a separate subroutine, or subroutines
               IQE=INT((E-PTEMIN)/DHISTE+1)
               IF (IQE.GT.0.AND.IQE.LT.NHISTE) THEN
                  NHISTQE(IQE,MYNODE)=NHISTQE(IQE,MYNODE)+1
               ENDIF
               DO K=1,NATOMS
                  Q(1,K)=X(K)
                  Q(2,K)=Y(K)
                  Q(3,K)=Z(K)
               ENDDO
               IF (PERIODIC.AND.CALCQT) THEN
                  CALL QORDER_BLJ(Q,Q4,Q6)
               ELSE IF (CHRMMT) THEN
                    IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
                    IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                      COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
                    IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                      COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.FALSE.) 
                    IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     &                      COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
               ELSEIF (CALCQT) THEN
                  CALL QORDER_LJ(Q,Q4,Q6)
               ENDIF
               Q4AV(MYNODE)=Q4AV(MYNODE)+Q4
               Q4AV2(MYNODE)=Q4AV2(MYNODE)+Q4**2
               Q6AV(MYNODE)=Q6AV(MYNODE)+Q6
               Q6AV2(MYNODE)=Q6AV2(MYNODE)+Q6**2

               !
               !todo move behind approriate flags - ODIHET, CHRMMT, OSASAT, ORGYT, etc 
               !
               DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)+DIHEORDERPARAM
               DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)+DIHEORDERPARAM**2
               SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)+SASAORDERPARAM
               SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)+SASAORDERPARAM**2
               RGYR_AV(MYNODE)=RGYR_AV(MYNODE)+RGYR
               RGYR_AV2(MYNODE)=RGYR_AV2(MYNODE)+RGYR**2
               EINT_AV(MYNODE)=EINT_AV(MYNODE)+EINT
               EINT_AV2(MYNODE)=EINT_AV2(MYNODE)+EINT**2
               IF (CHRMMT) THEN
                  IF (ODIHET) THEN
                     IQ4=INT((DIHEORDERPARAM-DIHEORDERPARAM_MIN)/DDIHE+1)
                     IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                        NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                     ENDIF
                     IQ6=INT((SASAORDERPARAM-SASAORDERPARAM_MIN)/DSASA+1)             
                     IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                        NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                     ENDIF
                  ELSE IF (ORGYT) THEN
                     IQ4=INT((RGYR-RGYR_MIN)/DRGYR+1)
                     IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                        NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                     ENDIF
                     IQ6=INT((EINT-EINT_MIN)/DEINT+1)             
                     IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                        NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                     ENDIF
                  ENDIF 
               ELSE
                  IQ4=INT(Q4/DQ4+1)
                  IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                     NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                  ENDIF
                  IQ6=INT(Q6/DQ6+1)             
                  IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                     NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                  ENDIF
               ENDIF

               ! Free energy statistic
               IF (CALCQT) THEN 
                 IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                    NHISTALLQ(IQ4,IQ6,MYNODE)=NHISTALLQ(IQ4,IQ6,MYNODE)+1
                 ENDIF
                 IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST.AND.
     &              IQE.GT.0.AND.IQE.LT.NHISTE) THEN
                    NHISTGRAND(IQE,IQ4,IQ6,MYNODE)=NHISTGRAND(IQE,IQ4,IQ6,MYNODE)+1
                 ENDIF
               ENDIF 
            ENDIF !closes if (MY_ORDER_PARAMETERT)
              
            !
            ! Dump visits histograms and restart information
            ! 
            IF (BSPTDUMPFRQ.GT.0) THEN
               IF (MOD(IMCSTEP,1.0D0*BSPTDUMPFRQ).EQ.0.0D0) THEN
                  ! dump visits histograms
                  CALL PTMC_DUMP_HISTOGRAM(PTEMIN, PEINT, PEVISITS, BINLABEL, QVISITS,
     &              TEMPTRAJ, IMCSTEP, .false.)

                  IF (BSPT.AND.(IMCSTEP.GT.NEQUIL+PTSTEPS).AND.(.NOT.MINDENSITYT)) THEN
                     call bspt_dump_histogram2(pevisits2, imcstep, .false.)
                  ENDIF

                  !  dump Restart information to bsptrestart
                  CALL PTMC_DUMP_RESTART_INFO(X, Y, Z, IMCSTEP, VOLD, VMINOLD, 
     &               NACCEPTPT, NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN)

                  ! js850> use this time to print a bit more information
                  IF (BINARY_EXAB) WRITE(MYUNIT, '(A,G15.1,A,2G15.1)') "bspt> binary_exab> tried ", 
     &               EXAB_COUNT, " accepted ", EXAB_ACC, IMCSTEP
               ENDIF
            ENDIF ! closes IF (BSPTDUMPFRQ.GT.0) 
         ENDIF ! closes IF (IMCSTEP.GT.NEQUIL), i think .. this IF block spans 2-3 pages!  
      ENDDO 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! End of main loop over BSPT or PT steps.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      WRITE (MYUNIT,'(A)') 'bspt> Exited main MC loop. ' 

! compute the averages over any plain PT steps

!     IF (PTSTEPS.GT.0) THEN
!        EAV(MYNODE)=EAV(MYNODE)/PTSTEPS
!        EAV2(MYNODE)=EAV2(MYNODE)/PTSTEPS
!        CV=(EAV2(MYNODE)-EAV(MYNODE)**2)*BETA(MYNODE)**2

!        WRITE (ISTR, '(i10)') MYNODE+1
!        FILENAME3=trim(adjustl(istr)) // "/T.Ev.Cv.Ev2.Steps"
!        LUNIT=GETUNIT()
!        OPEN(UNIT=LUNIT,FILE=FILENAME3, STATUS="unknown", form="formatted")
!        IF (CHRMMT) THEN
!           WRITE(LUNIT,'(6G20.10)') TEMPTRAJ(MYNODE),TEMPTRAJ(MYNODE)/0.001987, 
!    &                                EAV(MYNODE), CV, EAV2(MYNODE),PTSTEPS
!        ELSE
!           WRITE(LUNIT,'(5G20.10)') TEMPTRAJ(MYNODE),EAV(MYNODE),CV, EAV2(MYNODE),PTSTEPS
!        ENDIF
!        CALL FLUSH(LUNIT)
!        CLOSE(LUNIT)
!     ENDIF

      IF (MY_ORDER_PARAMETERT) THEN
         !js850> todo: this should go to a separate subroutine, or subroutines
         Q4AV(MYNODE)=Q4AV(MYNODE)/PTSTEPS
         Q6AV(MYNODE)=Q6AV(MYNODE)/PTSTEPS
         Q4AV2(MYNODE)=Q4AV2(MYNODE)/PTSTEPS
         Q6AV2(MYNODE)=Q6AV2(MYNODE)/PTSTEPS
         DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)/PTSTEPS
         SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)/PTSTEPS
         DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)/PTSTEPS
         SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)/PTSTEPS
         RGYR_AV(MYNODE)=RGYR_AV(MYNODE)/PTSTEPS
         EINT_AV(MYNODE)=EINT_AV(MYNODE)/PTSTEPS
         RGYR_AV2(MYNODE)=RGYR_AV2(MYNODE)/PTSTEPS
         EINT_AV2(MYNODE)=EINT_AV2(MYNODE)/PTSTEPS
         FQ4=SQRT(Q4AV2(MYNODE)-Q4AV(MYNODE)**2)
         FQ6=SQRT(Q6AV2(MYNODE)-Q6AV(MYNODE)**2)
         IF (BINARY .OR. SOFT_SPHERE) THEN  
            FILENAME4="T.Q4Av.Q6Av.Q4Av2.Q6Av2.Steps."//trim(adjustl(istr))
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE=FILENAME4, STATUS="unknown", form="formatted")
            WRITE(LUNIT,'(6G20.10)') TEMPTRAJ(MYNODE), Q4AV(MYNODE), Q6AV(MYNODE), Q4AV2(MYNODE),Q6AV2(MYNODE), PTSTEPS
            CLOSE(LUNIT)
         ENDIF
         IF (CHRMMT) THEN 
            IF (ODIHET) THEN 
               FILENAME4="T.DiheAv.SasaAv.DiheAv2.SasaAv2.Steps."//trim(adjustl(istr))
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE=FILENAME4, STATUS="unknown", form="formatted")
               WRITE(LUNIT,'(6G20.10)') TEMPTRAJ(MYNODE), DIHEORDERPARAM_AV(MYNODE), 
     &                 SASAORDERPARAM_AV(MYNODE), DIHEORDERPARAM_AV2(MYNODE),SASAORDERPARAM_AV2(MYNODE), PTSTEPS
               CLOSE(LUNIT)
               ELSE IF (ORGYT) THEN
               FILENAME4="T.RgyrAv.EintAv.RgyrAv2.Eint2.Steps."//trim(adjustl(istr))
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE=FILENAME4, STATUS="unknown", form="formatted")
               WRITE(LUNIT,'(6G20.10)') TEMPTRAJ(MYNODE), RGYR_AV(MYNODE),EINT_AV(MYNODE), RGYR_AV2(MYNODE),  
     &                 EINT_AV2(MYNODE), PTSTEPS
               CLOSE(LUNIT)
            ENDIF
         ENDIF

!     FILENAME5=TRIM(ADJUSTL(ISTR)) // '/profile_E'
!     LUNIT=GETUNIT()
!     OPEN(UNIT=LUNIT,FILE=FILENAME5, STATUS="unknown", form="formatted")
!     DO K=1,NHISTE
!        WRITE(LUNIT,'(2G20.10)') PTEMIN+(K-1)*DHISTE,NHISTQE(K,MYNODE)
!     ENDDO
!     CLOSE(LUNIT)

         IF (BINARY .OR. SOFT_SPHERE) THEN  
            FILENAME6="profile_Q4."//trim(adjustl(istr))
            FILENAME7="profile_Q6."//trim(adjustl(istr))
            LUNIT=GETUNIT()
            TEMPUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE=FILENAME6, STATUS="unknown", form="formatted")
            OPEN(UNIT=TEMPUNIT,FILE=FILENAME7, STATUS="unknown", form="formatted")
            DO K=1,NHIST
               WRITE(LUNIT,'(2G20.10)') (K-1)*DQ4,NHISTQ4(K,MYNODE)
               WRITE(TEMPUNIT,'(2G20.10)') (K-1)*DQ6,NHISTQ6(K,MYNODE)
            ENDDO
            CLOSE(LUNIT)
            CLOSE(TEMPUNIT)
         ELSE IF (CHRMMT) THEN
            IF (ODIHET) THEN
               FILENAME6="profile_dihe."//trim(adjustl(istr))
               FILENAME7="profile_sasa."//trim(adjustl(istr))
               LUNIT=GETUNIT()
               TEMPUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE=FILENAME6, STATUS="unknown", form="formatted")
               OPEN(UNIT=TEMPUNIT,FILE=FILENAME7, STATUS="unknown", form="formatted")
               DO K=1,NHIST
                  WRITE(LUNIT,'(2G20.10)') diheorderparam_min+(K-1)*Ddihe,NHISTQ4(K,MYNODE)
                  WRITE(TEMPUNIT,'(2G20.10)') SASAorderparam_min+(K-1)*Dsasa,NHISTQ6(K,MYNODE)
               ENDDO
               CLOSE(LUNIT)
               CLOSE(TEMPUNIT)
            ELSE IF (ORGYT) THEN
               FILENAME6="profile_rgyr."//trim(adjustl(istr))
               FILENAME7="profile_eint."//trim(adjustl(istr))
               LUNIT=GETUNIT()
               TEMPUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE=FILENAME6, STATUS="unknown", form="formatted")
               OPEN(UNIT=TEMPUNIT,FILE=FILENAME7, STATUS="unknown", form="formatted")
               DO K=1,NHIST
                  WRITE(LUNIT,'(2G20.10)') rgyr_min+(K-1)*Drgyr,NHISTQ4(K,MYNODE)
                  WRITE(TEMPUNIT,'(2G20.10)') eint_min+(K-1)*Deint,NHISTQ6(K,MYNODE)
               ENDDO
               CLOSE(LUNIT)
               CLOSE(TEMPUNIT)
            ENDIF
         ENDIF
      ENDIF !closes if (MY_ORDER_PARAMETERT)
!
! ss2029> close dumpenergy file
!
      IF ( PTMCDUMPENERT ) CLOSE(FILENAME_ENER_LUNIT) 
!
! Printing summary
!
      IF (OVERLAPK) WRITE(MYUNIT, '(A,3G15.3)') "bspt> overlap> summary ", OVERLAP_TIMETOT, 
     &          OVERLAP_COUNT, OVERLAP_TIMETOT/OVERLAP_COUNT
      IF (BINARY_EXAB) WRITE(MYUNIT, '(A,G15.1,A,G15.1)') "bspt> binary_exab> tried ", EXAB_COUNT, " accepted ", EXAB_ACC
      WRITE(MYUNIT, '(A,F15.1,A,G15.1,A,F15.5,A)') 'bspt> ',NACCEPTPT(MYNODE), ' steps accepted out of ', 
     &            PTSTEPS+NEQUIL+NQUENCH, ' i.e. ',NACCEPTPT(MYNODE)*100.0D0/(PTSTEPS+NEQUIL+NQUENCH),'%'
      WRITE(MYUNIT, '(A,G20.10)') 'bspt> Final stepsize ', STEP(MYNODE+1)
      WRITE(MYUNIT, '(A,G20.10,A,2G20.10,A)') 'bspt> ',NEACCEPT, ' PT exchanges accepted out of ', 
     &                                                 NTOT,NEACCEPT*100.0D0/(1.0D0*MAX(1,NTOT)),' %'
      IF (BSPT) WRITE(MYUNIT, '(A,G20.10,A,G20.10,A)') 'bspt> ',NOUTQBIN, ' quenches outside quench bin range: ',
     &                                      NOUTQBIN*100.0D0/(NQUENCH),' %'
      WRITE(MYUNIT, '(A,I8,A,G20.10,A)') 'bspt> ',NOUTPEBIN, ' potential energies outside bin range: ',
     &                                      NOUTPEBIN*100.0D0/(PTSTEPS+NEQUIL+NQUENCH),' %'
      IF (BSPT) THEN
         IF (PERCOLATET) THEN
            IF (MINDENSITYT) THEN
               WRITE(MYUNIT, '(A,G20.10,A,G20.10,A)') 'bspt> ',XOUT, ' steps failed to produce a percolating system: ',
     &                                         XOUT*100.0D0/(NQUENCH+NEQUIL),' %'
            ELSE
               WRITE(MYUNIT, '(A,G20.10,A,G20.10,A)') 'bspt> ',XOUT, ' steps failed to produce a percolating system: ',
     &                                            XOUT*100.0D0/(NQUENCH),' %'
            ENDIF
         ELSE
            IF (MINDENSITYT) THEN
               WRITE(MYUNIT, '(A,G20.10,A,G20.10,A)') 'bspt> ',XOUT, ' steps had an atom outside the container: ',
     &                                            XOUT*100.0D0/(NQUENCH+NEQUIL),' %'
            ELSE
               WRITE(MYUNIT, '(A,G20.10,A,G20.10,A)') 'bspt> ',XOUT, ' steps had an atom outside the container: ',
     &                                            XOUT*100.0D0/(NQUENCH),' %'
            ENDIF
         ENDIF
      ENDIF
      WRITE(MYUNIT, '(A,I8)') 'bspt> total number of quenches actually required=',NQ(MYNODE+1)
      CALL MYCPU_TIME(TIME2)
      write(myunit, '(A,G20.10)') "bspt> total time elapsed ", time2 - timestart
      write(myunit, '(A,G20.10,A,G20.10)') "bspt> time spent quenching ", 
     &   timequench, " % ", 100. * timequench / (time2 - timestart)
      write(myunit, '(A,G20.10,A,G20.10)') "bspt> time spent in tryexchange ", 
     &   time_tryexchange, " % ", 100. * time_tryexchange / (time2 - timestart)
      write(myunit, '(A,G20.10,A,G20.10,A,G20.10)') "bspt> time spent computing potential for metropolis ", 
     &   time_potential, " % ", 100. * time_potential / (time2 - timestart), " per call ", time_potential / npot_call

      IF (MY_ORDER_PARAMETERT) THEN
         DO K=1,NATOMS
            Q(1,K)=X(K)
            Q(2,K)=Y(K)
            Q(3,K)=Z(K)
         ENDDO
         IF (PERIODIC) THEN
            CALL QORDER_BLJ(Q,Q4,Q6)
         ELSE
            CALL QORDER_LJ(Q,Q4,Q6)
         ENDIF
      ENDIF
!
! Dump visits histograms. Energies are written for the middle of the bins.
! This is the really important part!
!
      CALL PTMC_DUMP_HISTOGRAM(PTEMIN, PEINT, PEVISITS, BINLABEL, QVISITS,
     &  TEMPTRAJ, IMCSTEP, .true.)

      IF (BSPT.AND.(.NOT.MINDENSITYT)) THEN
         call bspt_dump_histogram2(pevisits2, imcstep, .true.)
      ENDIF
!
! Dump final restart information.
!
      CALL PTMC_DUMP_RESTART_INFO(X, Y, Z, IMCSTEP, VOLD, VMINOLD, 
     &   NACCEPTPT, NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN)

      IF (OVERLAPK) CLOSE(OVERLAP_UNIT)

      RETURN
#else
      RETURN
#endif
      END SUBROUTINE PTBASINSAMPLING
