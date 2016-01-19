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
! 
!---======================================---
      SUBROUTINE PTMC_ONE_ATOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements parallel tempering with a one atom takestep
! routine.  The benefit of a one atom takestep routine is that we can take much
! larger steps and still get an energy that is accepted with reasonable
! probability. This leads to much more efficient sampling.  See Daan Frenkel's
! book.  The downside is that we need NATOMS monte carlo moves to do one monte
! carlo sweep.  Thus we need to take many more monte carlo moves, so each monte
! carlo move must be very fast.  A typical energy calculation is of order
! O(NATOMS**2), but this will destroy all the efficiency gained by the one atom
! takestep.  So the potential calculation needs to be rewritten to update the
! energy incrementally for the moved atoms, a calculation of order O(NATOMS).
! This requires a lot of book keeping to keep track of which atoms have moved,
! since the last potential evaluation.
!
! Finally, we must ensure that the everything else in the monte carlo step is as
! quick as possible.  Copying arrays is of order O(NATOMS), which, when compared
! with the standard potential call is trivial, but with one atoms moves can
! become a significan bottleneck.  So all array copies must be done
! incrementally as well --- only updating those atoms that have moved.  Any
! other O(NATOMS) operations should also be removed.
!
! Currently the only potential that works with one atom takestep is the binary
! lennard jones with smooth cutoff LJPSHIFT().  I use an auxiliary function
! LJPSHIFT_SET_MOVED2() which communicates with a neighbor list module to help
! keep track of which atoms have moved.  This would need to be modified to get
! it to work with another potential.
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
      INTEGER :: IACCEPT(0:NPAR-1), MPIERR, J,K, TEMPUNIT
      INTEGER NDUMMY, NDUMMY2, NTOT, NH, J1, J2, J3, J4, JLOW,
     1        N1,
     2        ITRAJ, ITRAJO,NEACCEPT, RNDSEED, NUPDATE,
     3        CONVERGED, LBFGS_ITERATIONS, BININDEX,
     4        IBININDEX, IBININDEX2,
     5        LBFGS_ITERATIONSO, PEVISITS2(NENRPER, HBINS, 0:NPAR-1),
     7        NOUTQBIN, NOUTPEBIN, LBFGS_ITERATIONSSAVE, HBINMAX, HBINMIN, GETUNIT, LUNIT
      DOUBLE PRECISION IMCSTEP
      DOUBLE PRECISION QVISITS(HBINS, 0:NPAR-1), NACCEPTPT(0:NPAR-1)
      DOUBLE PRECISION PEVISITS(NENRPER,0:NPAR-1) ! instantaneous energy histogram from different replicas  

      DOUBLE PRECISION TEMPTRAJ(0:NPAR-1), BETA(0:NPAR-1), 
     3        X(NATOMS), Y(NATOMS), Z(NATOMS), 
     4        CTE, T, VOLD,VNEW, POTEL, GRAD(3*NATOMS), RANDOM, DPRAND,
     5        DQ4, DQ6, DHISTE, XO(NATOMS), YO(NATOMS), ZO(NATOMS), DDX, DDY, DDZ, DE, 
     6        W, WCOMP, WAC, E, CV, FQ4, FQ6, RMAX, DDXN, DDYN, DDZN, R2, 
     7        DUMMY, BINLABEL(HBINS), VNEWSAVE, PEINT, HISTINT, 
     8        MINCOORDS(3*NATOMS,NPAR),
     9        DIST,
     A        DCOORDS(3*NATOMS), SR3, OPOTEL,
     B        DPRAND_UNIVERSAL

       DOUBLE PRECISION 
     5          VMINOLD, VMINNEW, CX, CY, CZ,
     6          VMINNEWSAVE
      double precision xyz(3*natoms) 

      CHARACTER (LEN=256)  FILENAME, FILENAME2,FILENAME3,FILENAME4,FILENAME5,FILENAME6,
     1                      FILENAME7,FILENAME8,FILENAME9, FILENAME10, FILENAME11, FILENAME12, 
     2                      FILENAME100,FILENAME101,FILENAME102, FILENAME103, FILENAME104, FILENAME105,
     3                      FILENAME106, FILENAMEDUMMY
      CHARACTER (LEN=80)  ISTR, SDUMMY
      CHARACTER (LEN=80)  CFNAME
      LOGICAL FITS, NEWENERGY, RECOUNT, OUTSIDE
      DOUBLE PRECISION RRX,RRY,RRZ,RRR

      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),ITIME
      DOUBLE PRECISION PREVSTEPS
      INTEGER OVERLAP_UNIT, OVERLAP_COUNT
      CHARACTER (LEN=200)  OVERLAP_FNAME
      DOUBLE PRECISION TIMESTART, OVERLAP_TIMETOT
      DOUBLE PRECISION TIME2, MCSTEPSTART, time1
      DOUBLE PRECISION timequench, time_tryexchange, time_potential, npot_call
      DOUBLE PRECISION RESTRICT_REGION_COUNT
      LOGICAL NORESET_TEMP
      INTEGER EXAB_COUNT, EXAB_ACC
      DOUBLE PRECISION IBOXLX, IBOXLY, IBOXLZ  !inverse box lengths
      INTEGER :: MOVED_ATOMS(NATOMS), NMOVED
      INTEGER :: MOVED_ATOMS_RESET(NATOMS), NMOVED_RESET
      INTEGER :: MOVED_ATOMS_TOT(NATOMS), NMOVED_TOT
      LOGICAL MOVEALL
      LOGICAL :: LDUM, CHECK_EQUAL, check_equal13 !a function
      LOGICAL EXCHANGEACCEPT
      LOGICAL SKIPPED_POTENTIAL
      logical :: debuglong = .false.

!     DOUBLE PRECISION PTCAND(PTSTEPS/10,3*NATOMS), PTCANDE(PTSTEPS/10)
      DOUBLE PRECISION QV(NENRPER), XDUMMY
      INTEGER NGMIN, NCHOSEN
      DOUBLE PRECISION :: LASTEXDOWN(0:NPAR-1)
      DOUBLE PRECISION :: BINARY_EXAB_RAND
      LOGICAL :: BINARY_EXAB_TRIED

      COMMON /MYPOT/ POTEL

!
! Initialisation
!

      RESTRICT_REGION_COUNT = 0.D0
      IBOXLX = 1.d0/BOXLX
      IBOXLY = 1.d0/BOXLY
      IBOXLZ = 1.d0/BOXLZ
      NMOVED = 0
      NMOVED_RESET = 0
      NMOVED_TOT = 0
      MOVEALL = .FALSE.

      EXAB_COUNT = 0
      EXAB_ACC = 0
      timequench = 0.
      time_tryexchange = 0.
      TIME_POTENTIAL = 0.
      npot_call = 0.

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NDUMMY,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      NOUTQBIN=0
      NOUTPEBIN=0
      LASTEXDOWN(0:NPAR-1)=-1.0D0

      WRITE(MYUNIT, '(A,I10,A,I10)') "bspt> This is processor", MYNODE+1, " of", NPAR
      WRITE(MYUNIT, '(A,I10)') 'bspt> Number of atoms', natoms
      IF (FIXSTEP(1)) FIXSTEP(2:NPAR)=.TRUE.
      IF (FIXTEMP(1)) FIXTEMP(2:NPAR)=.TRUE.
      IF (PERIODIC) THEN
         WRITE(MYUNIT, '(A,6G20.10)') 'bspt> Binary data', ntypea, epsab, epsbb, sigab, sigbb, cutoff
         WRITE(MYUNIT, '(A,3G20.10)') 'bspt> Box data', boxlx, boxly, boxlz
      ENDIF

      ITRAJ=MYNODE
      NEACCEPT=0 ! number of accepted PT exchanges

      ! PEINT is the PE bin width.
      PEINT=(PTEMAX-PTEMIN)/NENRPER
      MAXEFALL=-1.0D100
      IF (BSPT) THEN
         !HISTINT is the quenched energy bin width
         HISTINT=(HISTMAX-HISTMIN)/HBINS
!
! Limit for quench energy above which steps will be rejected.
!
!     HBINMAX=(MIN(BSPTQMAX,HISTMAX)-HISTMIN)/HISTINT+1
!     HBINMIN=(MAX(BSPTQMIN,HISTMIN)-HISTMIN)/HISTINT
         HBINMIN=1
         HBINMAX=HBINS
         WRITE(MYUNIT, '(A,2I8)') 'bspt> Minimum and maximum quench bin limits: ',HBINMIN,HBINMAX
         WRITE(MYUNIT, '(A,2G20.10)') 'bspt> Minimum and maximum quench energies: ',histmin, histmax
         DO J1=1, HBINS
            !BINLABEL(J1) is the lowest energy in a quench bin J1
!           BINLABEL(J1)=HISTMIN + HISTINT*(J1-0.5D0) ! these energies point to the middle of the quench bin
            BINLABEL(J1)=HISTMIN + HISTINT*(J1-1.0D0) ! these energies point to the bottom of the quench bin
         ENDDO
      ENDIF

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


      ! set up the temperatures
      IF (PTREADTEMPS) THEN
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT, FILE='temperatures.init',STATUS='old')
         DO J1=0, NPAR-1
            READ(LUNIT,*) T
            TEMPTRAJ(J1) = T
            BETA(J1)=1.0D0/T
         ENDDO
         CLOSE(LUNIT)
      ELSE
         CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
         CTE=EXP(CTE)
         DO J1=0, NPAR-1
            TEMPTRAJ(J1)=PTTMIN*CTE**J1
            T=TEMPTRAJ(J1)
            BETA(J1)=1.0D0/T
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
      XYZ(:) = COORDS(:,MYNODE+1)


      WRITE(MYUNIT, '(A,2G20.10)') 'bspt> Temperature range             = ', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A,G20.10)')  'bspt> This temperature trajectory   = ', TEMPTRAJ(MYNODE)
      IF (.NOT.BSPTRESTART) VMINOLD=VOLD
      IF ((.NOT.BSPTRESTART).AND.(NEQUIL+PTSTEPS.GT.0)) VMINOLD=0.0D0 
      WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of equilibration steps = ', NEQUIL
      WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of PT steps            = ', PTSTEPS
      IF (BSPT) THEN
         WRITE (MYUNIT,'(A,F20.1)')   'bspt> Number of BSPT steps          = ',NQUENCH
         WRITE (MYUNIT,'(A,G20.10)')   'bspt> Replica Exchange Interval    = ',EXCHINT 
         WRITE (MYUNIT,'(A,G20.10)')   'bspt> Replica Exchange Prob        = ',EXCHPROB
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
         IF (BSPT) PEVISITS2(:,:,MYNODE)=0
      ENDIF
      IACCEPT(MYNODE)=0

!
! NUPDATE specifies the interval for dynamically altering the maximum step size.
! Only used if step size isn't fixed.
!
! ss2029> NUPDATE is hardcoded. read from a keyword. todo.  
!
      NUPDATE=100
      IF (.NOT.(FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1).OR.((MYNODE.LE.USERES).AND.RESERVOIRT)))  THEN 
          WRITE(MYUNIT, '(A,I6,A)') 'bspt> Step size will be adjusted every ',NUPDATE ,' MC steps'
      ENDIF
      WRITE (MYUNIT,'(A)') 'bspt> Entering main MC loop ... ' 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop over steps. We have to use double precision for the number of steps, because it
! may exceed integer capacity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMCSTEP=0.0D0    ! = MC step number 
      IF (BSPTRESTART) IMCSTEP=PREVSTEPS
      MCSTEPSTART = IMCSTEP
      CALL MYCPU_TIME(timestart)
      DO 
         IMCSTEP=IMCSTEP+1.0D0
         IF (IMCSTEP.GT.NEQUIL+PTSTEPS+NQUENCH) EXIT
         RECOUNT=.FALSE.
         IF (ONE_ATOM_TAKESTEP .AND. NMOVED .GT. 0 .AND. .NOT. MOVEALL) THEN
            DO J1=1, NMOVED
               K = MOVED_ATOMS(J1)
               XO(K)=X(K)
               YO(K)=Y(K)
               ZO(K)=Z(K)
            ENDDO
         ELSE
            DO K=1, NATOMS
               XO(K)=X(K)
               YO(K)=Y(K)
               ZO(K)=Z(K)
            ENDDO
         ENDIF
         IF (DEBUG) then
            LDUM = .true.
            if (.not. check_equal(XO, X, natoms)) LDUM = .false.
            if (.not. check_equal(YO, Y, natoms)) LDUM = .false.
            if (.not. check_equal(ZO, Z, natoms)) LDUM = .false.
            if (.not. LDUM) then
               write(MYUNIT, *) "bspt> warning: problem coppying array XO <- X"
            endif
         endif
!
! Should probably worry about centre of coordinates versus centre of mass
! for heteroatomic systems.
!
         !IF ((DEBUG.OR.CHECKMARKOVT).AND.(IMCSTEP.GT.1.0D0).AND.(.NOT.CANDIDATELISTT)) THEN
            !CALL POTENTIAL(COORDSO(:,MYNODE+1),GRAD,OPOTEL,.TRUE.,.FALSE.)
            !IF (ABS(OPOTEL-VNEW).GT.ECONV) THEN
               !WRITE(MYUNIT,'(3(A,G20.10))') 'bspt> possible ERROR - energy for coordinates in COORDSO=',
!     &                        OPOTEL,' but Markov energy=',VNEW,' IMCSTEP=',IMCSTEP
               !!STOP
            !ENDIF
         !ENDIF
!
! If a move is rejected, the result for the current configuration
! should be recounted. This allows for the unsymmetrical forward
! and backward move probabilities for an atom near the surface of
! the constraining sphere.
!             

!
! take the step
!
         IF (BINARY_EXAB) THEN
            BINARY_EXAB_RAND = DPRAND()
            BINARY_EXAB_TRIED = .FALSE.
         ENDIF

         IF ( BINARY_EXAB .AND. (BINARY_EXAB_RAND .le. 1.D0/BINARY_EXAB_FRQ) ) THEN
            BINARY_EXAB_TRIED = .TRUE.
            !js850> try to exchange a type A and type B particle as the MC step
            EXAB_COUNT = EXAB_COUNT + 1
            !choose an A and a B particle
            IF ( FREEZE ) THEN
              !get a mobile type B particle
              !NMOBILETYPEB = NTYPEB - NFREEZETYPEB = (NATOMS - NTYPEA) - (NFREEZE - NFREEZETYPEA)
              J3 = (NATOMS - NTYPEA) - (NFREEZE - NFREEZETYPEA) !NMOBILETYPEB
        !      write(*,*) NATOMS , NTYPEA,  NFREEZE , NFREEZETYPEA
              !J4=       *( NMOBILETYPEB ) +1 + NFREEZE + (NMOBILETYPEA)
              RANDOM=DPRAND()
              J4 = RANDOM*( J3           ) +1 + NFREEZE + (NTYPEA-NFREEZETYPEA)
              J2 = FROZENLIST(J4)
              !get a mobile type A particle
              !NMOBILETYPEA = NTYPEA - NFREEZETYPEA
              !J1=       *( NMOBILETYPEA ) +1 + NFREEZE
              RANDOM=DPRAND()
              J4 = RANDOM*( NTYPEA-NFREEZETYPEA ) + 1 + NFREEZE
              J1 = FROZENLIST(J4)
            ELSE
              RANDOM=DPRAND()
              J1 = RANDOM*NTYPEA + 1
              RANDOM=DPRAND()
              J2 = RANDOM*(NATOMS-NTYPEA) + 1 + NTYPEA
            ENDIF
            !WRITE(MYUNIT,'(A,2I8,F20.1)') "BINARY_EXAB> trying exchange ", J1, J2, IMCSTEP
            !IF ( FROZEN(J1) .OR. FROZEN(J2) .OR. J1 .GT. NTYPEA .OR. J2 .LE.  NTYPEA ) THEN
              !WRITE(*,*) "problem with binary_exab", j1, j2, frozen(j1), frozen(j2), nfreeze
              !CALL EXIT()
            !endif
            !exchange the xyz coordinates of J1 and J2
            NMOVED=2
            MOVED_ATOMS(1) = J1
            MOVED_ATOMS(2) = J2
            CALL LJPSHIFT_SET_MOVED2( NMOVED, MOVED_ATOMS, .FALSE., myunit)
            DUMMY = X(J1)
            X(J1) = X(J2)
            X(J2) = DUMMY
            DUMMY = Y(J1)
            Y(J1) = Y(J2)
            Y(J2) = DUMMY
            DUMMY = Z(J1)
            Z(J1) = Z(J2)
            Z(J2) = DUMMY
            IF (DEBUG .and. debuglong) then
               write(myunit,*) "binary_exab> exchanging atoms", j1,j2
            endif
         ELSE
            !js850> move only one atom               
            IF (FREEZE) THEN
               !select from among the unfrozen atoms
               RANDOM=DPRAND()
               J2 = INT(RANDOM*(NATOMS - NFREEZE)) + 1
               J2 = J2 + NFREEZE
               K = FROZENLIST(J2)
               if (FROZEN(K) .or. K .gt. natoms) then
                  write(myunit,*) "error: ONE_ATOM_TAKESTEP> problem with atom selection", K
               endif
            ELSE
               RANDOM=DPRAND()
               K = INT(RANDOM*NATOMS) + 1
            ENDIF
            NMOVED = 1
            MOVED_ATOMS(NMOVED) = K
            CALL LJPSHIFT_SET_MOVED2( NMOVED, MOVED_ATOMS, .FALSE., myunit)
            !write(*,*) "moving atom ", K, MYNODE+1
            RANDOM=DPRAND()
            X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            RANDOM=DPRAND()
            Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            RANDOM=DPRAND()
            Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            IF (DEBUG .and. debuglong) then
               write(myunit,*) "one_atoms_takestep> moving atom", K
            endif
         ENDIF

         OUTSIDE=.FALSE.


         ! copy X Y Z into XYZ for use in POTENTIAL
         IF (ONE_ATOM_TAKESTEP .AND. .NOT. MOVEALL) THEN
            !THIS PARTIAL COPY IS DANGEROUS BECAUSE XYZ IS USED AS SCRATCH
            DO J1=1, NMOVED
               K = MOVED_ATOMS(J1)
               XYZ(3*(K-1)+1) = X(K)
               XYZ(3*(K-1)+2) = Y(K)
               XYZ(3*(K-1)+3) = Z(K)
            ENDDO
         ELSE
            !write(MYUNIT,*) "doing full copying 2 XYZ <- X"
            DO K=1,NATOMS
               XYZ(3*(K-1)+1) = X(K)
               XYZ(3*(K-1)+2) = Y(K)
               XYZ(3*(K-1)+3) = Z(K)
            ENDDO
         ENDIF
         IF (DEBUG) THEN
            if (.not. check_equal13(XYZ, x, y, z, natoms)) then
               write(MYUNIT, *) "bspt> warning problem copying array XYZ <- X, Y, Z"
            endif
         endif

! js850> Check to see if RESTRICTREGION constraints have been violated - if they have, reject the step
         IF ( RESTRICTREGION ) THEN
           !RESTRICTREGIONTEST=.FALSE.
           DO K=1,NMOVED
              J1 = MOVED_ATOMS(K)
              RRX = ( X(J1)-RESTRICTREGIONX0 )
              RRY = ( Y(J1)-RESTRICTREGIONY0 )
              RRZ = ( Z(J1)-RESTRICTREGIONZ0 )
              RRX = RRX - ANINT(RRX*IBOXLX)*BOXLX
              RRY = RRY - ANINT(RRY*IBOXLY)*BOXLY
              IF ( RESTRICTCYL ) THEN
                 RRZ = 0.D0
                 RRR = DSQRT(RRX**2+RRY**2 )
              ELSE
                 RRZ = RRZ - ANINT(RRZ*IBOXLZ)*BOXLZ
                 RRR = DSQRT(RRX**2+RRY**2+RRZ**2 )
              ENDIF
              IF ( RRR > RESTRICTREGIONRADIUS ) THEN
                 RECOUNT=.TRUE.
                 RESTRICT_REGION_COUNT = RESTRICT_REGION_COUNT + 1.D0
                 RRX = sqrt(XO(j1)**2+YO(J1)**2+zo(j1)**2) ! for debugging
                 IF (DEBUG) WRITE(MYUNIT, '(A,F20.1,I8,F20.12,5F14.8)') 'bspt> restrictregion> rejecting step ', IMCSTEP
     &              , j1, Vold, X(j1), y(j1), z(j1), RRR, RRX
                 !GOODSTRUCTURE=.FALSE.
                 !RESTRICTREGIONTEST=.TRUE.
                 !EXIT !exit loop
              ENDIF
           ENDDO
         ENDIF
!
! At this point all we have done is take a step, energies have not been computed yet and 
! also Metropolis check has not been done. 
!
! The perturbed coordinates are in both XYZ and X, Y, Z. 
! The old coordinates are in XO, YO, ZO. 
!
!
! New and old quench energies will be in        VMINNEW VMINOLD
! New and old instantaneous energies will be in    VNEW    VOLD
!
         SKIPPED_POTENTIAL = .FALSE.
         IF (RECOUNT) THEN !no need to call POTENTIAL if step is already rejected
            VNEW=0.0D0
            SKIPPED_POTENTIAL = .TRUE.
         ELSE
            IF (DEBUG) THEN
               if (.not. check_equal13(XYZ, x, y, z, natoms)) then
                  write(MYUNIT, *) "bspt> warning problem about to call potential: XYZ .NE. X,Y,Z"
               endif
            endif
            CALL MYCPU_TIME(TIME1)
            CALL POTENTIAL(XYZ, GRAD, POTEL, .FALSE., .FALSE.)
            CALL MYCPU_TIME(TIME2)
            npot_call = npot_call + 1
            TIME_POTENTIAL = TIME_POTENTIAL + TIME2 - TIME1
            VNEW=POTEL
            IF (DEBUG .OR. INT(MOD(IMCSTEP,1000000.D0)).EQ.0) THEN
               CALL POTENTIAL(XYZ, GRAD, DUMMY, .TRUE., .FALSE.)
               IF (ABS(POTEL-DUMMY).GT.1.D-7) THEN
                  WRITE(MYUNIT,'(A,G20.10,A,3G20.10)') 'bspt> ERROR - one atom potential gives energy',
     &                           POTEL,' but standard potential gives ',DUMMY, POTEL-DUMMY,IMCSTEP
                  !STOP
               else
                  WRITE(MYUNIT,'(A,G20.10,A,3G20.10)') 'bspt> ALL OK- one atom potential gives energy',
     &                           POTEL,' and standard potential gives ',DUMMY, POTEL-DUMMY,IMCSTEP
               ENDIF
            ENDIF
            ! ss2029> this is where Metropolis check on instantaneous energy is done. 
            WCOMP=(VNEW-VOLD)*BETA(MYNODE) ! use difference in instantaneous energies
            W=MIN(1.0D0,EXP(-WCOMP))
            RANDOM=DPRAND()
            IF (RANDOM.GT.W) RECOUNT=.TRUE. ! RECOUNT is initialised to .FALSE. at the top of the loop
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
            !only need to move a few atoms back
            DO J1=1,NMOVED
               K = MOVED_ATOMS(J1)
               X(K)=XO(K)
               Y(K)=YO(K)
               Z(K)=ZO(K)
               XYZ(3*(K-1)+1)=X(K)
               XYZ(3*(K-1)+2)=Y(K)
               XYZ(3*(K-1)+3)=Z(K)
            ENDDO

            IF (.NOT. SKIPPED_POTENTIAL) THEN
               !record which atoms have been put back in their old positions
               !but only if the potential was told about the moved atoms
               NMOVED_RESET = NMOVED
               MOVED_ATOMS_RESET(1:NMOVED) = MOVED_ATOMS(1:NMOVED)
               CALL LJPSHIFT_SET_MOVED2( NMOVED_RESET, MOVED_ATOMS_RESET, .FALSE., myunit)
            ENDIF

            VNEW=VOLD
            IF (BSPT) THEN
               VMINNEW=VMINOLD
               LBFGS_ITERATIONS=LBFGS_ITERATIONSO
            ENDIF
         ELSE ! accept move
            NMOVED_RESET = 0
            NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1.0D0
            IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
            IF ( BINARY_EXAB .AND. BINARY_EXAB_TRIED ) THEN
               EXAB_ACC = EXAB_ACC + 1
               write(MYUNIT,*) "BINARY_EXAB> exchange accepted"
            ENDIF
         ENDIF
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
            IF (.NOT. CHECK_EQUAL13(XYZ, X, Y, Z, NATOMS)) THEN
               write(MYUNIT, *) "bspt> warning problem after acc/rej: XYZ .NE. X,Y,Z"
            ENDIF
            !CALL POTENTIAL(XYZ,GRAD, POTEL, .TRUE., .FALSE.)
            !WRITE(MYUNIT,'(A,3G20.10)') 'bspt> potential energy after acc/rej=',POTEL, VOLD, VNEW
         ENDIF
         IF (DEBUG) CALL FLUSH(MYUNIT)
!
!  At this point the quench and instantaneous energies for the current
!  configuration in the Markov chain are VMINNEW and VNEW.
!
!
         IF (BSPT.and.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.1.0D0*0)) BININDEX=INT((VMINNEW-HISTMIN)/HISTINT)+1
         IBININDEX=INT((VNEW-PTEMIN)/PEINT)+1
!
!  Must not accumulate statistics until we have equilibrated for NEQUIL steps.
!
         IF (IMCSTEP.GT.NEQUIL) THEN
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
!              X, Y, Z, XYZ          energy=VNEW
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


            IF ( PTMCDUMPSTRUCT .AND. (MOD(IMCSTEP-1.0D0,PTMCDS_FRQ*1.0D0).EQ.0.0D0) ) THEN
               ! js850> PTMCDUMPSTRUCT
               CALL PTMC_DUMPSTRUCT_DUMP(XYZ, VNEW, IMCSTEP)
            ENDIF
            IF ( OVERLAPK .AND. (MOD(IMCSTEP-1.0D0,OVERLAP_FRQ*1.0D0).EQ.0.0D0) ) THEN
               ! js850> OVERLAP keyword
               CALL PTMC_OVERLAP_DUMP(XYZ, OVERLAP_UNIT, MYUNIT, OVERLAP_COUNT,
     &            time1, time2, VNEW, IMCSTEP)
               overlap_timetot = overlap_timetot + time1 
               timequench = timequench + time2 
            ENDIF

         ENDIF

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
            WRITE(MYUNIT, '(A,G20.10,A,G20.10,A,I6,A)') 'bspt> adjusting step-size> current step size = ', STEP(MYNODE+1) ,
     &                        ' acceptance ratio = ', WAC ,
     &                        ' over ', NUPDATE, ' steps'
         ENDIF
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

         LBFGS_ITERATIONSO=LBFGS_ITERATIONS
         E=VNEW
!
!  Replica exchange part is done in TRYEXCHANGE.
!
         MOVEALL = .FALSE.
         IF (EXCHPROB.GT.0.0D0) THEN
            CALL MYCPU_TIME(TIME1)
            CALL TRYEXCHANGE(E,X,Y,Z,XO,YO,ZO,VOLD,EXCHANGEACCEPT,JLOW,
     &                  VNEW,GRAD,VNEWSAVE,VMINOLD,VMINNEW,BETA,ITRAJ,ITRAJO,NTOT, 
     &                  LBFGS_ITERATIONS,NEACCEPT,LBFGS_ITERATIONSO,QV,XDUMMY,PEINT,NCHOSEN,LASTEXDOWN,IMCSTEP)
            CALL MYCPU_TIME(TIME2)
            time_tryexchange = time_tryexchange + (time2 - time1)
            IF (EXCHANGEACCEPT) then
               !write(myunit,*) "moveall set to true: exchange accepted"
               MOVEALL = .TRUE.
               CALL LJPSHIFT_SET_MOVED2( 0, MOVED_ATOMS_RESET, MOVEALL, myunit)
               IF (DEBUG .OR. MOD(INT(NTOT),1000).EQ.0) THEN
                  !check to make sure the markov energy agrees with the energy
                  !in coords
                  DO K=1,NATOMS
                     DCOORDS(3*(K-1)+1) = X(K)
                     DCOORDS(3*(K-1)+2) = Z(K)
                     DCOORDS(3*(K-1)+3) = Y(K)
                  ENDDO
                  CALL POTENTIAL(DCOORDS(:),GRAD,DUMMY,.TRUE.,.FALSE.)
                  IF (ABS(DUMMY-VOLD).GT.1.D-7) THEN
                     WRITE(MYUNIT,'(3(A,G20.10))') 'bspt> ERROR - after PT EXCHANGE markov energy is',
     &                              VNEW, 'but energy of X,Y,Z is', DUMMY,' IMCSTEP=',IMCSTEP
                     !STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         !  js850>At this point VOLD, VNEW  are the markov energy (perhaps changed)
         !                      X           is the markov coords
         !                      E, VNEWSAVE, XO are not determined.  They should not be used after this
         E = VNEW !js850> fix a bug here
         !ENDIF

         !!!!!!!!!!!!!!!!!!! REPLICA EXCHANGE end !!!!!!!!!!!!!!!!!!

         IF (IMCSTEP.GT.NEQUIL) THEN 
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
                  IF (RESTRICTREGION) WRITE(MYUNIT, '(A,G20.1,A,F20.1)') 
     &              "bspt> restrict_region> rejected ", RESTRICT_REGION_COUNT, " step ", IMCSTEP
               ENDIF
            ENDIF ! closes IF (BSPTDUMPFRQ.GT.0) 
         ENDIF ! closes IF (IMCSTEP.GT.NEQUIL), i think .. this IF block spans 2-3 pages!  
      ENDDO 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! End of main loop over BSPT or PT steps.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      WRITE (MYUNIT,'(A)') 'bspt> Exited main MC loop. ' 

! compute the averages over any plain PT steps

!
! Printing summary
!

      IF (RESTRICTREGION) WRITE(MYUNIT, '(A,F20.1,A,F20.1)') 
     &              "bspt> restrict_region> rejected ", RESTRICT_REGION_COUNT, " step ", IMCSTEP
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
      WRITE(MYUNIT, '(A,I8)') 'bspt> total number of quenches actually required=',NQ(MYNODE+1)
      CALL MYCPU_TIME(TIME2)
      write(myunit, '(A,G20.10)') "bspt> total time elapsed ", time2 - timestart
      write(myunit, '(A,G20.10,A,G20.10)') "bspt> time spent quenching ", 
     &   timequench, " % ", 100. * timequench / (time2 - timestart)
      write(myunit, '(A,G20.10,A,G20.10)') "bspt> time spent in tryexchange ", 
     &   time_tryexchange, " % ", 100. * time_tryexchange / (time2 - timestart)
      write(myunit, '(A,G20.10,A,G20.10,A,G20.10)') "bspt> time spent computing potential for metropolis ", 
     &   time_potential, " % ", 100. * time_potential / (time2 - timestart), " per call ", time_potential / npot_call


!
! Dump visits histograms. Energies are written for the middle of the bins.
! This is the really important part!
!
      CALL PTMC_DUMP_HISTOGRAM(PTEMIN, PEINT, PEVISITS, BINLABEL, QVISITS,
     &  TEMPTRAJ, IMCSTEP, .true.)

      IF (BSPT) THEN
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
      END SUBROUTINE PTMC_ONE_ATOM

      function check_equal(array1, array2, n) result(b)
         implicit none
         integer n
         double precision array1(n), array2(n)
         integer i
         logical b
         b = .true.
         do i=1,n
            if (array1(i) .ne. array2(i))  then
               b = .false.
               return
            endif
         enddo
      end function check_equal

      function check_equal13(coords, x, y, z, n) result(b)
         implicit none
         integer n
         double precision coords(3*n), x(n), y(n), z(n)
         integer i
         logical b
         b = .true.
         do i=1,n
            if (coords(3*(i-1)+1) .ne. x(i))  then
               b = .false.
               return
            endif
            if (coords(3*(i-1)+2) .ne. y(i))  then
               b = .false.
               return
            endif
            if (coords(3*(i-1)+3) .ne. z(i))  then
               b = .false.
               return
            endif
         enddo
      end function check_equal13

