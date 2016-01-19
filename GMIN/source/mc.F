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
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
      USE COMMONS
      use genrigid
      use moves
      use grouprotmod
      use mbpolmod, only: mbpolstep

      USE QMODULE , ONLY : QMIN, QMINP, INTEQMIN
      USE modcharmm
      USE MODAMBER9, ONLY : MDSTEPT,CISARRAY1,CHIARRAY1,NOCISTRANSDNA,NOCISTRANSRNA,
     &                      SETCHIRAL,AMCHPMAX,DOLIGMOVE,LIGMOVEFREQ, AMCHNMAX, ROTAMERT,
     &                      E_IGB, E_BOND, E_ANGLE, E_DIHEDRAL,
     &                      E_VDW, E_14_VDW, E_ELEC, E_14_ELEC, IGB
      USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_WRITE_RESTART, AMBER12_WRITE_PDB, POT_ENE_REC_C
      USE CHIRALITY, ONLY : INIT_CHIRAL, INIT_CIS_TRANS
      USE porfuncs
      USE AMHGLOBALS, ONLY: NMRES,OMOVI,AVEP,NUMPRO,IRES
      USE AMH_INTERFACES, ONLY:E_WRITE
      USE ROTAMER

      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
      INTEGER MPIERR
      LOGICAL HITANY
#endif

      INTEGER J1, NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), J2, NSTEPS, JP, 
     1        UNT, ITERATIONS, NSUPERCOUNT, NQTOT, JACCPREV, NREN, NLAST, NSTEPREN, BRUN,QDONE,JBEST(NPAR),GCJBEST(NPAR),
     2        NRMS, NDONE, I, RNDSEED, J, NTOT, ITRAJ, NEACCEPT, J3, J4, ISTAT, LOCALCOUNT, QNEWRES, QGCBH, GCNATOMSBEST(NPAR),
     &        GCJESTP(NPAR), GCNATOMSBESTP(NPAR), GCJBESTP(NPAR)
      INTEGER :: NSYMCALL=0, NULLMOVES(NPAR), NULLMOVEST(NPAR), NSWAPS=0
      DOUBLE PRECISION POTEL, SCALEFAC, RANDOM, DPRAND, DPRAND_UNIVERSAL, SAVECOORDS(3*NATOMSALLOC), TEMPCOORDS(3*NATOMSALLOC),
     1                 TIME, SPOTEL(NSUPER), SCOORDS(3*NATOMSALLOC,NSUPER), SCREENC(3*NATOMSALLOC),
     2                 EPPREV(NPAR), QSTART, QFINISH, RANNJ, RMIN, RMINO, RCOORDS(3*NATOMSALLOC),ELASTSYM(NPAR),
     3                 RCOORDSO(3*NATOMSALLOC), RVAT(NATOMSALLOC), RVATO(NATOMSALLOC), EPSSAVE, EBEST(NPAR), GCEBEST(NPAR), 
     4                 GCEBESTP(NPAR),
     4                 BESTCOORDS(3*NATOMSALLOC,NPAR), RMSD, VINIT, CTE, TEMPTRAJ(0:NPAR-1), GCBESTCOORDS(3*NATOMSALLOC,NPAR),
     5                 T, BETA(0:NPAR-1), GRAD(3*NATOMSALLOC), E, W, GCBESTCOORDSP(3*NATOMSALLOC,NPAR),
     &                 DUMMY1, DUMMY2, DUMMY3, INTE, OPOTEL, GCBESTVAT(NATOMSALLOC,NPAR), GCBESTVATP(NATOMSALLOC,NPAR), 
     &                 XMASS, YMASS, ZMASS, CMMAX
      LOGICAL CHANGEDE, RES1, RES2
      LOGICAL CHIRALFAIL,AMIDEFAIL, LOGDUMMY, DISTOK
      CHARACTER FNAME*9
      CHARACTER (LEN= 3)  ISTR
      CHARACTER (LEN=20) QUENCHNUM, QUNAME
      CHARACTER (LEN=20) BESTNAME
      LOGICAL  CHANGED
c  AMH 
      INTEGER :: GLY_COUNT,III,I2,I500,SNAPCOUNT
      DOUBLE PRECISION prcord(NATOMSALLOC,3,3,3)
      DOUBLE PRECISION :: mctemp
!  csw34> PAIRDIST variables
      INTEGER :: PAIRCOUNTER
      DOUBLE PRECISION, EXTERNAL :: PAIRDISTANCE
      DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      LOGICAL EVAP, ATEST, STAY, evapreject, LOPEN, COMPLETE
      COMMON /EV/ EVAP, evapreject
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      COMMON /Q4C/ QSTART, QFINISH

      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),itime1
      DOUBLE PRECISION :: DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE,DUMMY
      integer :: J6, L, M, SUMSQUAREDIFF, SUMSQUAREDIFF1D, HBONDDIFF
      INTEGER GETUNIT
      INTEGER DUMPUNIQUEUNIT

! khs26> Energy decomposition for AMBER 12
      TYPE(POT_ENE_REC_C) :: AMBER12_ENERGY_DECOMP

! ds656 > energy before and after BGUPTA swap (without relaxation)
      DOUBLE PRECISION :: EBSWAP=0.0D0, EASWAP=0.0D0, EASWAPQ=0.0D0, DE1=0.0D0, DE2=0.0D0

      QNEWRES=0
      QGCBH=0

! csw34> initialise RANDOM to 0.0D0 to avoid printing uninitialised variable if not set in
!        TRANSITION (e.g. if Qu 1 causes cold fusion)
      RANDOM=0.0D0

      IF (FREEZEGROUPT) THEN
         OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
!
! Work out the distance from GROUPCENTRE to the current atom J1
! 
            DISTGROUPX2=(COORDS(3*GROUPCENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*GROUPCENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*GROUPCENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
            IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
! If working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
            ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
            END IF
         END DO
         CLOSE(4431)
! Prevent it doing this again
         FREEZEGROUPT=.FALSE.     
      ENDIF

! Write a list of DONTMOVE atoms for use in an (o)data file
      IF (DONTMOVEGROUPT) THEN
              OPEN(UNIT=4431,FILE='dontmove.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
!
! Work out the distance from DONTMOVECENTRE to the current atom J1
! 
            DISTGROUPX2=(COORDS(3*DONTMOVECENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*DONTMOVECENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*DONTMOVECENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
! If working in GT mode (default), DONTMOVE all atoms >GROUPRADIUS from the DONTMOVECENTRE atom
            IF((DONTMOVEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
! IF working in LT mode, DONTMOVE all atoms <GROUPRADIUS from the DONTMOVECENTRE atom
       ELSE IF((DONTMOVEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
            END IF
         END DO
         CLOSE(4431)
! Prevent it doing this again
         DONTMOVEGROUPT=.FALSE.     
      ENDIF
      
! js850> initialize DUMPUNIQUE variables
      IF ( DUMPUNIQUE ) THEN
        DUMPUNIQUEEPREV = 0.D0
        DUMPUNIQUEUNIT = GETUNIT()
        OPEN(UNIT=DUMPUNIQUEUNIT,FILE="dumpunique",STATUS='UNKNOWN')  
      END IF
      
!     csw34> Added defaults to prevent accidentaly discarding
!     structures for AMH      

      CHIRALFAIL=.FALSE.
      AMIDEFAIL=.FALSE.
      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mc> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      ALLOCATE(TMOVE(NPAR), OMOVE(NPAR))
      snapcount=0
      NSTEPREN=0
      EVAPREJECT=.FALSE.
      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mc> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

#ifdef MPI
      IF (MPIT) THEN 
         IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'MPIERR=',MPIERR
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
         IF (DEBUG) WRITE(MYUNIT, '(A,2I6)') 'NPAR,MPIERR=',NPAR,MPIERR
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
         JP=MYNODE+1
         ITRAJ=MYNODE
         IF (DEBUG) WRITE(MYUNIT, '(A,3I6)') 'In mc after MPI_MPI_COMM_RANK MPIERR,MYNODE,JP=',MPIERR,MYNODE,JP
      ENDIF
      NEACCEPT=0
      HITANY=.FALSE.
#endif

      NDONE=0
      IF (RESTORET) THEN
#ifdef MPI
         CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
#else
         DO JP=1,NPAR
            CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
         ENDDO
#endif
         WRITE(MYUNIT, '(A,I10)') 'MC> restore NDONE=',NDONE
!     csw34> Sets the quench counter so that the GMIN_out file makes sense after using RESTORE!
!         NQ(:)=NDONE
      ENDIF
      NQ(:)=NDONE

! tvb requesting a basin-sampling MC run: 
      
      IF (BSWL.AND.(.NOT.TETHER)) then
         CALL BASINSAMPLING()
         RETURN
      ELSEIF (TETHER) THEN
         CALL TETHEREDWL()
         RETURN
      ENDIF

#ifdef MPI
      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Processor", mynode+1, " of", NPAR, " speaking:"
      WRITE(MYUNIT, '(A,I10)') 'Number of atoms', natoms
#endif

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      NRMS=0
      NLAST=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM
#ifdef MPI
#else
      DO JP=1,NPAR 
#endif
         TMOVE(JP)=.TRUE.
         OMOVE(JP)=.TRUE.
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
         IF (JDUMP(JP).AND.(.NOT.NEWJUMP)) THEN
            WRITE(FNAME,'(A,I1)') 'ebuffer.',JP
            UNT=70+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
            WRITE(FNAME,'(A,I1)') 'cbuffer.',JP
            UNT=70+NPAR+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
         ENDIF
#ifdef MPI
#else
      ENDDO
#endif

      IF (AMHT) THEN
         write(omovi,1334)nmres,3,1,INT(real(NSTEPS)/real(NINT_AMH))
1334     format(4(i8,1x),' nmres nmcrd numpro nmsnap')
      ENDIF
    
      IF (.NOT.RESTORET) THEN
!       csw34> Set the centre of mass to be at the specified location
!       contained in the SETCENTRE X Y Z keyword
         IF (SETCENT) CALL SETCENTRE(COORDS)
!
! For MAKEOLIGOT and MAKEOLIGOSTART=TRUE: generate oligomers by placing new segments.
         IF (CHRMMT.AND.MAKEOLIGOT.AND.MAKEOLIGOSTART) THEN
#ifdef MPI
!        seed the random number generator with system time  + MYNODE (for MPI runs)
             IF (RANDOMSEEDT) THEN
                CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                itime1= values(6)*60 + values(7)
                CALL SDPRND(itime1+MYNODE)
             ENDIF
             CALL CHMAKEOLIGOMER(JP)
             call flush(6)
#else
             DO JP=1,NPAR
                CALL CHMAKEOLIGOMER(JP)
             ENDDO
#endif
         ENDIF
      ENDIF

!  Calculate the initial energy and save in EPREV
#ifdef MPI
      WRITE(MYUNIT, '(A)')  'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
      CALL FLUSH(MYUNIT)

 
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
      EPREV(JP)=POTEL
      EPPREV(JP)=0.0D0
      ELASTSYM(JP)=0.0D0
      IF (.NOT.RESTORET) EBEST(JP)=POTEL
      BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      JBEST(JP)=0
      RMIN=POTEL
      RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
      COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
      EPSSPHERE=EPSSAVE

! Initialisation 

      IF (PTTMIN < 1.0D-6 ) PTTMIN = 1.0D-6 ! to avoid devision by zero
      CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
      CTE=EXP(CTE)

      DO I=0, NPAR-1
         TEMPTRAJ(I)=PTTMIN*CTE**I
         T=TEMPTRAJ(I)
         BETA(I)=1.0D0/T
      ENDDO
      DO I=1, NPAR
         TEMP(I)=TEMPTRAJ(I-1)
      ENDDO
      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
      VINIT=POTEL

      WRITE(MYUNIT, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A, G20.10)') 'For this replica T=', TEMPTRAJ(MYNODE)
      WRITE(MYUNIT, '(A, G20.10)') 'Starting potential energy=', VINIT

      IF (.NOT.RANSEEDT) THEN
         RNDSEED=2002+MYNODE
         CALL SDPRND(RNDSEED)
         CALL SDPRND_UNIVERSAL(2002+NPAR)
      ENDIF
      RANDOM=DPRAND()
      WRITE(MYUNIT, '(A, G20.10)') 'Starting random number=', RANDOM
      RANDOM=DPRAND_UNIVERSAL()
      WRITE(MYUNIT, '(A, G20.10)') 'Starting universal random number=', RANDOM
      WRITE(ISTR,'(I3)') JP
      TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
      IF ((CHRMMT).AND.(.NOT.RESTORET)) CALL CHARMMDUMP(TEMPCOORDS,'initialmin.'//TRIM(ADJUSTL(ISTR)))

#else
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
!fh301>{{{
         COUT=.FALSE.
!fh301>}}}
         IF (PERMOPT) THEN
            PERMOPT=.FALSE.
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.) 
            PERMOPT=.TRUE.
         ENDIF
! abc> Added initial call to set_check_chiral for L/D info for initial structure
         IF (AMBERT.AND.SETCHIRAL) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing chiral information for initial structure'
            CALL set_check_chiral(COORDS(:,1),NATOMS,LOGDUMMY,chiarray1)
         END IF
! khs26> Generating lists of chiral atoms, and checking chirality for the initial structure.
         IF (AMBER12T .AND. SETCHIRAL) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing chiral information for initial (quench) structure'
            CALL INIT_CHIRAL(COORDS(:,1))
         END IF
         IF (AMBER12T.AND.NOCISTRANS) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing cis/trans information for initial (quench) structure'
            CALL INIT_CIS_TRANS(COORDS(:,1))
         ENDIF
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
! swap coordinates of atom 2 and 3 to test the chirality check
!                  TEMPCOORDS=COORDS(:,1)
!                  COORDS(4:6,1)=TEMPCOORDS(7:9)
!                  COORDS(7:9,1)=TEMPCOORDS(4:6)
        NQTOT=NQTOT+1

!cl457 Qu 0 output
         IF (CHEMSHIFT2) THEN
           IF(NPAR.GT.1) THEN
              WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1           POTEL,' ECAM=', ENERGYCAMSHIFT,' steps=',ITERATIONS,' RMS=',RMS,
     1           ' Markov E=',POTEL,' t=',TIME-TSTART
           ELSE
              WRITE(MYUNIT,'(A,I10,A,G20.10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' ECAM=', ENERGYCAMSHIFT,' steps=',ITERATIONS,' RMS=',RMS,
     1           ' Markov E=',POTEL,' t=',TIME-TSTART
           END IF
         ELSE
!cl457
           IF (NPAR.GT.1) THEN
              WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
           ELSE
              WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
           ENDIF
         ENDIF

! js850> initial dump
      IF (DUMPUNIQUE) THEN
        CALL DUMPUNIQUE_DUMP(DUMPUNIQUEUNIT, JP, POTEL, DUMPUNIQUEEPREV )
      ENDIF

! Added dump of the initial structure, this is very useful for flu! csw34
      IF (CHRMMT) CALL CHARMMDUMP(COORDS,'initialmin')

!     csw34> Added initial call to check_cistrans_protein to store cis/trans info for initial structure
         IF (AMBERT.AND.NOCISTRANS.AND.(.NOT.NOCISTRANSDNA).AND.(.NOT.NOCISTRANSRNA)) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing cis/trans information for initial structure'
            CALL check_cistrans_protein(COORDS(:,1),NATOMS,LOGDUMMY,MINOMEGA,cisarray1)
         ENDIF
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
         EPREV(JP)=POTEL
         EPPREV(JP)=0.0D0
         ELASTSYM(JP)=0.0D0
         IF (.NOT.RESTORET) EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=0
         RMIN=POTEL
         RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
         COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
!hk286 - otherwise NEWRESTART is always called when the keyword RESTORE is also used
         IF (.NOT.RESTORET) THEN
            EBEST(JP)=POTEL
            JBEST(JP)=0
            BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         ENDIF
!
! Set up data for relxation blocks with grand canonical BH
!
         IF (GCBHT) THEN
            GCEBEST(JP)=POTEL
            GCJBEST(JP)=0
            GCBESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
            GCBESTVAT(1:NATOMS,JP)=VAT(1:NATOMS,JP)
            GCNATOMSBEST(JP)=NATOMS
            GCEBESTP(JP)=POTEL
            GCJBESTP(JP)=0
            GCBESTCOORDSP(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
            GCBESTVATP(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               GCNATOMSBESTP(JP)=NATOMS
         ENDIF
! hk286
      ENDDO
      EPSSPHERE=EPSSAVE
#endif

      IF (THOMSONT .AND. (.NOT. GTHOMSONT)) THEN
!
! Scale maximum step size for the Thomson problem according to the mean nearest-neighbour 
! distance after the first quench.
!
         DUMMY1=0.0D0
         DO J1=1,NATOMS
            DUMMY2=1.0D100
            DO J2=1,NATOMS
               IF (J2.EQ.J1) CYCLE
               DUMMY3=(COORDS(3*(J1-1)+1,1)-COORDS(3*(J2-1)+1,1))**2 +  
     &                (COORDS(3*(J1-1)+2,1)-COORDS(3*(J2-1)+2,1))**2 + 
     &                (COORDS(3*(J1-1)+3,1)-COORDS(3*(J2-1)+3,1))**2 
               IF (DUMMY3.LT.DUMMY2) DUMMY2=DUMMY3
            ENDDO
            DUMMY1=DUMMY1+DUMMY2
         ENDDO
         DUMMY1=SQRT(DUMMY1/NATOMS)
         DO J1=1,NPAR
            STEP(J1)=STEP(J1)*DUMMY1
         ENDDO
         WRITE(MYUNIT, '(2(A,G20.10))') 'Maximum step size scaled by mean nearest neighbour distance of ',DUMMY1,' to ',STEP(1)
      ENDIF

      IF (HIT) GOTO 37 

      IF (NPAR.EQ.1) THEN
         WRITE(MYUNIT,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      ELSE
         WRITE(MYUNIT,'(A,I3,A,I10,A)') 'Starting ',NPAR,' parallel MC runs of ',NSTEPS,' steps'
      ENDIF
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'
!
! csw34> Before we start BH steps, 
! check the structure satisfies the conditions in LOCALSAMPLE if specified
!        This is important to prevent an infinite loop occuring!
!        Also, make sure that the ligand has been specified in movableatoms
      IF (LOCALSAMPLET.AND.AMBERT) THEN
         IF(NMOVABLEATOMS==0) THEN
            WRITE(MYUNIT,*) 'must have MOVABLEATOMS specified when using LOCALSAMPLE' 
            STOP
         ENDIF
         DISTOK=.FALSE.
         CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
         IF (.NOT.DISTOK) THEN
            WRITE(MYUNIT,*) 'initial structure violates LOCALSAMPLE conditions - exiting!'
            STOP
         ENDIF
         DISTOK=.FALSE.
      ENDIF   
      NSUPERCOUNT=NSUPER
! HBONDMATRIX
! csw34> Use the initial quenched structure to set up the first hydrogen-bond group if using HBONDMATRIX
      IF ((AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))).AND.HBONDMATRIX) THEN
! Dump the coordinates in .rst format to hbond.rst
         CALL DUMPRST(COORDS(:,1),NATOMS,'hbond')
! Calculate hydrogen-bond matrix using ptraj defaults (specified in initialization.f90)
         CALL GETHBONDMATRIX('coords.prmtop','hbond.rst',TRIM(ADJUSTL(HBONDDONORSACCEPTORS)), 
     &                       TRIM(ADJUSTL(HBONDRESIDUES)),TRIM(ADJUSTL(HBONDDCUT)),
     &                       TRIM(ADJUSTL(HBONDACUT)),HBONDMAT)
! Set the number of groups to 1
         NHBONDGROUPS=1
! Save the group matrix
         HBONDGROUPS(1,:,:)=HBONDMAT
! Initialize the population of the group to 1
         HBONDGROUPPOP(1)=1
! Save the initial energy for final output
         HBONDQUZEROE=POTEL
! Set the group max energy
         HBONDMAXE(1)=POTEL
! Set the group best and markov energy and save the coordinates
         HBONDBEST(1)=POTEL
         HBONDBESTCOORDS(1,:)=COORDS(:,1)
         HBONDMARKOV(1)=POTEL
! Write out the matrix for reference
         CALL WRITEHBONDMATRIX(HBONDMAT(:,:),HBONDNRES,'group1')
! Print some output
         WRITE(MYUNIT,'(A)') " HBONDMATRIX> Group 1 created using initial quench structure"
       ENDIF
! END OF HBONDMATRIX

! csw34> Dump the initial quench structure to .rst and .pdb files if DUMPQU is specified
       IF(AMBERT.AND.DUMPQUT) THEN
          CALL DUMPRST(COORDS(:,1),NATOMS,'initialmin') 
          CALL A9DUMPPDB(COORDS(:,1),'initialmin') 
       ENDIF

! khs26> AMBER12 dumping restart and PDB files
       IF(AMBER12T.AND.DUMPQUT) THEN
          CALL AMBER12_WRITE_RESTART(COORDS(:,1),'initialmin.rst',LEN('initialmin.rst')) 
          CALL AMBER12_WRITE_PDB(COORDS(:,1),'initialmin.pdb',LEN('initialmin.pdb')) 
       ENDIF

! jdf43>
      IF (TBP) THEN
         WRITE(MYUNIT,'(A,I0.2,A)')"[",JP,"]TBP>Initialising basins..."
         TBPBASINS(:) = 0
      ENDIF
      IF (RATIOT) THEN
         NULLMOVES=0
         NULLMOVEST=0
      ENDIF

! mo361> Use SYMMETRIZE in GA searches if we're not doing any basin-hopping
      IF (GENALT.AND.(MCSTEPS(1).EQ.0).AND.SYMMETRIZE) THEN
         JP=1
         J1=1
         NSYMCALL=0
         ELASTSYM(JP)=EPREV(JP)
         IF (SYMMETRIZECSM) THEN
            CALL SYMMETRYCSM(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,
     &                       J1,EBEST,BESTCOORDS,JBEST)
         ELSE
            CALL SYMMETRY(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,
     &                    J1,NSUCCESS,NFAIL,EBEST,BESTCOORDS,JBEST,EPPREV)

         ENDIF
      ENDIF

      J1=0
      COMPLETE=.FALSE.

!  Main basin-hopping loop 
      DO WHILE (.NOT.COMPLETE)
         J1=J1+1
!sf344> the line below was commented out, and caused infinite step taking for MPI runs. COMPLETE is not set
!       to .TRUE. in this routine at all, only on this line! 18.10.2013 
         IF ((J1.GE.NSTEPS).AND.(.NOT.REPMATCHT)) COMPLETE=.TRUE.
         IF (MAX_NPCALLT .AND. (NPCALL .GE. MAX_NPCALL)) THEN
            write(myunit,"(A,I18,A,I18,A)") "mc> MAX_NPCALL: maximum number of potential calls reached ( ", 
     &           npcall, " > ", max_npcall, " ).  Ending MC loop"
            GOTO 37
         ENDIF
!      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1
         QNEWRES=QNEWRES+1
         QGCBH=QGCBH+1

! hk286 - corrected 9/11/11
         IF (RELAXRIGIDT) THEN
            IF( (MOD(J1,NRELAXRIGIDR + NRELAXRIGIDA) ==  1) ) THEN
               RIGIDINIT = .TRUE.
               WRITE (MYUNIT, *) 'calling rigid quench'
               CALL CHECKSITES(.TRUE.,COORDS(:,1))
            ELSE IF ( (MOD(J1,NRELAXRIGIDR + NRELAXRIGIDA) == NRELAXRIGIDR + 1) ) THEN
               RIGIDINIT = .FALSE.
               WRITE (MYUNIT, *) 'calling atom quench'
               !CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            ENDIF
         ENDIF
! hk286   

         CALL FLUSH(MYUNIT)
         IF (NEWJUMP) RANNJ=DPRAND()
!
!  ********************************* Loop over NPAR parallel runs ******************************
!
#ifdef MPI
#else
         DO JP=1,NPAR
#endif

!       csw34> If QUCENTRE is specified, move the centre of coordinates
!       to (0,0,0) before taking the next step (improve this so that you
!       can specify where to move the centre like SETCENTRE?
            IF (QUCENTRE) THEN 
!       csw34> David mentioned a possible compiler bug causing problems
!       when you just pass a part of the COORDS array so reading the
!       right bit into TEMPCOORDS first and then back out.
               TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
               CALL CENTRE2(TEMPCOORDS)
               COORDS(1:3*NATOMS,JP)=TEMPCOORDS(1:3*NATOMS)
            ENDIF
            IF (RIGID.AND.(BLOCK(JP).GT.0)) THEN
               IF (MOD(J1-1,BLOCK(JP)).EQ.0) THEN
                  IF (MOD((J1-1)/BLOCK(JP),2).EQ.0) THEN
                     IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body translational moves'
                     TMOVE(JP)=.TRUE.
                     OMOVE(JP)=.FALSE.
                  ELSE 
                     IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body angular moves'
                     OMOVE(JP)=.TRUE.
                     TMOVE(JP)=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
!
!  MAM1000> The default temperature used for the MC acceptance criterion is the one derived in
!           the initialisation section above.  MCTEMP is passed to the subroutine TRANSITION, where
!           the acceptance/rejection decision is made.  However, individual MC moves can override
!           this temperature by setting MCTEMP to a different value for the current step.
            MCTEMP = TEMP(JP)

! jdf43> 07/01/12 TEMPERATURE BASIN PAVING.
!        Adjusts MCTEMP based on frequency of visitation of energy bins
!        TBPBASINS.
            IF (TBP) THEN
               IF (MOD(J1-TBPCI,TBPCI).EQ.1) THEN
                  IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I0.2,A)')"[",JP,"]TBP>Cooling basins"
                  TBPBASINS = INT(TBPBASINS-TBPCF*TBPBASINS)
               ENDIF
               TBPBASIN = INT((EPREV(JP)-TBPMIN)/TBPSTEP)+1
               IF (TBPBASIN.LE.TBPSTEPS) THEN
                  TBPINCR=TBPHF*TBPBASINS(TBPBASIN)
                  MCTEMP = TEMP(JP)*(1+TBPINCR)
                  TBPBASINS(TBPBASIN) = TBPBASINS(TBPBASIN)+1
               ELSEIF (EPREV(JP).LT.TBPMIN) THEN
                  IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I0.2,A)')"[",JP,"]TBP>TBPMIN is too damn high"
               ELSE
                  IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I0.2,A)')"[",JP,
     &            "]TBP>EPREV doesn't fit into a basin - consider increasing TBPSTEP or TBPSTEPS"
               ENDIF
            ENDIF
! 
!  Jump-moves.
! 
            IF (JUMPMOVE(JP).AND.(MOD(J1,JUMPINT(JP)).EQ.0)) CALL JUMPM(RANNJ,J1,JP,EPPREV)
! 
!  Ordinary steps.
! 
23          CONTINUE
!
! Don;t call symmetry unless the minimum in the Markov chain has changed.
! We should really check if the minimum has changed since the last call to SYMMETRY,
! which can be done with ABS(ELASTSYM(JP)-EPREV(JP)) if NSYMINTERVAL=1.
!
!           IF (SYMMETRIZE.AND.(MOD(J1-1,NSYMINTERVAL).EQ.0).AND.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) THEN
!           WRITE(MYUNIT,'(A,3G20.10)') 'ELASTSYM,EPREV,diff=',ELASTSYM(JP),EPREV(JP),ABS(ELASTSYM(JP)-EPREV(JP))
            IF (SYMMETRIZE.AND.(MOD(NQ(JP),NSYMINTERVAL).EQ.0).AND.
     &               ((SYMMETRIZECSM).OR.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV))) THEN
               IF ((ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) NSYMCALL=0
               ELASTSYM(JP)=EPREV(JP)
!              IF ((.NOT.MOVESHELLT).OR.(NSURFMOVES(JP).LT.0)) THEN
                  IF (SYMMETRIZECSM) THEN
                     CALL SYMMETRYCSM(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,  
     &                                J1,EBEST,BESTCOORDS,JBEST)
                  ELSE
                     CALL SYMMETRY(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,  
     &                                J1,NSUCCESS,NFAIL,EBEST,BESTCOORDS,JBEST,EPPREV)
                  ENDIF
                  IF (HIT) GOTO 37 
!              ELSE
!                 CALL SYMMETRY2(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL)
!              ENDIF
!              WRITE(MYUNIT,'(A,I2,A,2I6)') '[',JP,']mc> NCORE: ',NCORE(1:NPAR)
!              IF (HIT) GOTO 37 ! hit cannot change in symmetry2 
!
!  Check for reseeding.
! 
               POTEL=EPREV(JP) ! NEWRES assumes POTEL is the energy of the current structure in COORDS
               IF (CHANGEDE.AND.NEWRESTART) THEN
                  CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST,RES1,RES2)
                  IF(RES1.OR.RES2) QNEWRES=0
               ENDIF
            ELSEIF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
               NSYMREM=0                                       ! Should therefore reset NSYMREM.
            ENDIF
!
!  Check for potential improvement for grand canonical relaxation block.
!
            IF (GCBHT) THEN
               IF (POTEL.LT.GCEBEST(JP)) THEN
                  IF (GCEBEST(JP)-POTEL.GT.ECONV) GCJBEST(JP)=J1
                  GCEBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
                  GCBESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                  GCBESTVAT(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               ENDIF
               IF (MOD(J1-1,PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,2G20.10)') 'mc> POTEL,GCEBEST(JP)=',POTEL,GCEBEST(JP)
            ENDIF
! csw34> Temporary calls to new step taking routines for testing
            IF (NEWMOVEST) THEN
! GROUPROTATION moves are in here for testing CUDAGMIN only
               IF (GROUPROTT) CALL GROUPROTSTEP(JP) 
! Each new move should take a final, optional argument - ATOM_LIST
! This is the list of atoms in a integer array to which you want to
! apply the move. 
               CALL CARTESIAN_SPHERE(COORDS(:,JP),STEP(JP))
!               CALL CARTESIAN_SPHERE(COORDS(:,JP),STEP(JP),ATOM_LIST)
!               CALL CARTESIAN_SIMPLE(COORDS(:,JP),STEP(JP))
!               CALL CARTESIAN_SIMPLE(COORDS(:,JP),STEP(JP),ATOM_LIST)
            ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF STEP TAKING CALLS!                                                                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set switch variables for moves that are not always performed
               IF(GROUPROTT.AND.MOD(J1,GROUPROTFREQ).EQ.0) THEN
                       DOGROUPROT=.TRUE.
               ENDIF
               IF(DIHEDRALROTT.AND.MOD(J1,DIHEDRALROTFREQ).EQ.0) THEN
                       DODIHEDRALROT=.TRUE.
               ENDIF
! hk286 - damped group moves
               IF(DAMPEDGMOVET.AND.MOD(J1,DMOVEFREQ).EQ.0) THEN
                       DODGROUPMOVET=.TRUE.
               ENDIF
! csw34> rigid body expansion move switch 
               IF(EXPANDRIGIDT.AND.MOD(J1,EXPANDRIGIDFREQ).EQ.0) THEN
                       DOEXPANDRIGID=.TRUE.
               ENDIF
! csw34> rigid body rotation move switch
               IF(ROTATERIGIDT.AND.MOD(J1,ROTATERIGIDFREQ).EQ.0) THEN
                       DOROTATERIGID=.TRUE.
               ENDIF
! mo361> rigid body rotation move switch
               IF(TRANSLATERIGIDT.AND.MOD(J1,TRANSLATERIGIDFREQ).EQ.0) THEN
                       DOTRANSLATERIGID=.TRUE.
               ENDIF
!
! csw34> Coordinates are saved so that moves can be undone
!
               SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
! csw34> If you want to look at the effect of moves, you can dump out
! the structure BEFORE the move here.
!                 CALL A9DUMPPDB(COORDS(:,JP),"beforemove")
!                 CALL CHARMMDUMP(COORDS(:,JP),'beforemove')

!
! UNIVERSAL STEP TAKING
!

! csw34> Rigid body EXPAND moves
                     IF (RIGIDINIT.AND.EXPANDRIGIDT.AND.DOEXPANDRIGID) THEN
                        CALL GENRIGID_EXPAND(COORDS(:,JP), EXPANDFACTOR) 
                     ENDIF
! csw34> Rigid body rotation moves
                     IF (RIGIDINIT.AND.ROTATERIGIDT.AND.DOROTATERIGID) THEN
                        CALL GENRIGID_ROTATE(COORDS(:,JP), ROTATEFACTOR)
                     ENDIF
! mo361> Rigid body translation moves
                     IF (RIGIDINIT.AND.TRANSLATERIGIDT.AND.DOTRANSLATERIGID) THEN
                        CALL GENRIGID_TRANSLATE(COORDS(:,JP),TRANSLATEFACTOR)
                     ENDIF
!
! CHARMM STEP TAKING
!
               IF (CHRMMT) THEN
                  IF (CHMDT.AND.MOD(J1,CHMDFREQ).EQ.0) THEN
                     CALL CHMD(JP)
                  ELSE
                     IF (CHRIGIDTRANST.AND.MOD(J1,FTRANS).EQ.0) CALL MKRIGIDTRANS(JP)
                     IF (CHRIGIDROTT.AND.MOD(J1,FROT).EQ.0) CALL MKRIGIDROT(JP)
                     IF (MOD(J1,CHFREQ).EQ.0) THEN
! hk286 - Do damped group moves
                        IF (DAMPEDGMOVET.AND.DODGROUPMOVET) THEN
                           CALL DAMPEDMOVESTEP(JP)
                        ELSEIF (CHCARTMOVET) THEN
! Cartesian moves with GROUPROTATION; settings of CHBBT and CHSCT do not matter 
                           IF (GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP) 
! Cartesian moves without GROUPROTATION; settings of CHBBT and CHSCT do not matter 
                           CALL TAKESTEPCH(JP)
                        ELSE
! dihedral moves
                        IF (CHBBT.AND.(.NOT.CHSCT)) THEN
! backbone moves only
                           CALL TAKESTEPCH(JP)
                        ELSEIF ((.NOT.CHBBT).AND.CHSCT) THEN
! sidechain moves only
                           IF (GROUPROTT.AND.DOGROUPROT) THEN
                              CALL GROUPROTSTEP(JP)
                           ELSEIF (ROTAMERT) THEN
!bs360: rotamer library is not working yet for CHARMM
                              WRITE(*,*)'Rotamer library is not working yet for CHARMM'
                              STOP
!                              CALL CHROTMOVE(JP)
                           ELSE
                              CALL TAKESTEPCH(JP)
                           ENDIF
                        ELSEIF (CHBBT.AND.CHSCT) THEN
                           IF (CHFREQBB.EQ.1 .AND. CHFREQSC.EQ.1) THEN
! BB and SC moves at the same time; no GROUPROTATIONS in this case
                              CALL TAKESTEPCH(JP)
                           ELSEIF (CHFREQBB.EQ.2 .AND. CHFREQSC.EQ.2) THEN
! alternating BB and SC moves; start with backbone moves
                              IF (MOD(J1,2).EQ.0) THEN
                                 CHBBT=.FALSE.
                                 IF (GROUPROTT.AND.DOGROUPROT) THEN
                                    CALL GROUPROTSTEP(JP)
                                 ELSE
                                    CALL TAKESTEPCH(JP)
                                 ENDIF
                                 CHBBT=.TRUE.
                              ELSE
                                 CHSCT=.FALSE.
                                 CALL TAKESTEPCH(JP)
                                 CHSCT=.TRUE.
                              ENDIF
                           ELSEIF (CHFREQBB.EQ.1 .AND. CHFREQSC.GE.3) THEN
! every CHFREQSCth step SC move, otherwise BB move
                              IF (MOD(J1,CHFREQSC).EQ.0) THEN
                                 CHBBT=.FALSE.
                                 IF (GROUPROTT.AND.DOGROUPROT) THEN
                                    CALL GROUPROTSTEP(JP)
                                 ELSE
                                    CALL TAKESTEPCH(JP)
                                 ENDIF
                                 CHBBT=.TRUE.
                              ELSE
                                 CHSCT=.FALSE.
                                 CALL TAKESTEPCH(JP)
                                 CHSCT=.TRUE.
                              ENDIF
                           ELSEIF (CHFREQBB.GE.3 .AND. CHFREQSC.EQ.1) THEN
! every CHFREQBBth step BB move, otherwise SC move
                              IF (MOD(J1,CHFREQBB).EQ.0) THEN
                                 CHSCT=.FALSE.
                                 CALL TAKESTEPCH(JP)
                                 CHSCT=.TRUE.
                              ELSE
                                 CHBBT=.FALSE.
                                 IF (GROUPROTT.AND.DOGROUPROT) THEN
                                    CALL GROUPROTSTEP(JP)
                                 ELSE
                                    CALL TAKESTEPCH(JP)
                                 ENDIF
                                 CHBBT=.TRUE.
                              ENDIF
                           ELSE
                              WRITE(MYUNIT,*) 'Your choices for CHFREQBB and CHFREQSC are not supported.'
                              STOP
                           ENDIF
                        ENDIF
                        ENDIF ! if (dampedgmovet.and.dodgroupmovet)
                     ENDIF    ! if (mod(j1,chfreq).eq.0)
                  ENDIF       ! if (chmdt.and.mod(j1,chmdfreq).eq.0)
! bs360: check the perturbed structure prior minimization (for debugging reasons)
#ifdef MPI
                     IF (DEBUG) CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'aftermove.'//TRIM(ADJUSTL(ISTR)))
#else
                     IF (DEBUG) CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'aftermove')
#endif
!               seed the random number generator with system time  + MYNODE (for MPI runs)
                  IF (RANDOMSEEDT) THEN
                     CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                     itime1= values(6)*60 + values(7)
                     CALL SDPRND(itime1+MYNODE)
                  END IF
!
! ANISOTROPIC STEP TAKING
!
               ELSE IF (CHIROT) THEN
                CALL CHIRO_TAKESTEP(JP) 
        
               ELSE IF (ELLIPSOIDT.AND..NOT.PYGPERIODICT.AND..NOT.PYBINARYT) THEN
                CALL TAKESTEPELLIPSOIDS(JP) 

!-------------------------------------------
! hk286 - generalised THOMSON
!-----------------------------------------------
               ELSE IF (GTHOMSONT) THEN
                CALL TAKESTEPGTHOMSON ()

               ELSE IF (GAYBERNEDCT) THEN
                CALL TAKESTEPGB(JP)

               ELSE IF (PYT) THEN
!               seed the random number generator with system time + MYNODE (for MPI runs)
                  IF(RANDOMSEEDT) THEN
                     CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                     itime1= values(6)*60 + values(7)
                     CALL SDPRND(itime1+MYNODE)
                  END IF
                call py_takestep(jp)

               ELSE IF (PYGPERIODICT.OR.PYBINARYT.OR.LJCAPSIDT) THEN
!               seed the random number generator with system time + MYNODE (for MPI runs)
                  IF(RANDOMSEEDT) THEN
                     CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                     itime1= values(6)*60 + values(7)
                     CALL SDPRND(itime1+MYNODE)
                  END IF
                  IF(SWAPMOVEST) THEN
                     CALL TAKESTEPSWAPMOVES(JP)
                     CALL TAKESTEPELPSD(JP)
                  ELSE
                     CALL TAKESTEPELPSD(JP)
                  END IF

!
! MARK'S MIXED CLUSTER STEP TAKING
!
!     MAM1000 >  Charged/neutral particle swaps for LJ+Coulomb mixed clusters
               ELSE IF (LJCOULT) THEN
                  RANDOM=DPRAND()
                  IF (RANDOM < COULSWAP) THEN
                     CALL TAKESTEPLJC(JP)
                     MCTEMP = COULTEMP
                  ELSE
                     CALL TAKESTEP(JP)
                  END IF
                  
!
! ds656> Identity swaps for binary systems. If BASWAP is FALSE or
!        J1 <= NWAIT then the usual TAKESTEP will be called further down.

               ELSE IF (BASWAP .AND. (J1 > BASWAP_NWAIT)) THEN
                  !
                  RANDOM=DPRAND()
                  !IF(BASWAP_FRAC > 0.0D0) THEN ! no longer needed!!!
                  !   RANDOM=RANDOM/DEXP( (QMIN(1)-EPREV(JP))/BASWAP_TEMP )
                  !ENDIF 
                  !
                  IF(BASWAP_NMAX > 0 .AND. RANDOM < BASWAP_FRAC ) THEN
                     !
                     NSWAPS = 1 + INT(BASWAP_NMAX*DPRAND())
                     !     
                     DO J4=1,NSWAPS
                        CALL ABSWAP(JP, DE1, DE2)
                     ENDDO
                     !
                     IF(BASWAPTEST) THEN
                        !
                        IF(BGUPTAT) THEN
                           CALL BGUPTA(COORDS,GRAD,EASWAP,.FALSE.)
                        ELSEIF(BLJCLUSTER_NOCUT) THEN
                           CALL BLJ_CLUST(COORDS,GRAD,EASWAP,.FALSE.)
                        ENDIF
                        !
                     ENDIF
                     !
                  ELSEIF(BASWAP_NMAX==0 .AND. RANDOM < BASWAP_FRAC) THEN
                     !
                     NSWAPS = 1
                     !
                     !CALL BGUPTA(COORDS,GRAD,EBSWAP,.FALSE.)
                     CALL ABSWAP2(JP, DE1, DE2)
                     !
                     IF(BASWAPTEST .AND. BGUPTAT) THEN
                        CALL BGUPTA(COORDS,GRAD,EASWAP,.FALSE.)
                     ELSEIF(BASWAPTEST .AND. BLJCLUSTER_NOCUT) THEN
                        CALL BLJ_CLUST(COORDS,GRAD,EASWAP,.FALSE.)
                     ENDIF
                     !     
                  ELSE
                     CALL TAKESTEP(JP)
                  END IF
! <ds656
!               ELSE IF (LWOTPT) THEN
!                  CALL TAKESTEPLWOTP(JP)

               ELSE IF (GBT .OR. GBDT .OR. GBDPT) THEN
                  CALL TKSTDCELPSD(JP)
      
               ELSE IF (MSGBT) THEN
                  CALL TAKESTEPMSGB(JP) 

               ELSE IF (SANDBOXT) THEN
                  IF(RANDOMSEEDT) THEN
                    CALL DATE_AND_TIME(DATECHAR,TIMECHAR,ZONECHAR,VALUES)
                    ITIME1 = VALUES(6) * 60 + VALUES(7)
                    CALL SDPRND(ITIME1 + MYNODE)
                  END IF
                  CALL SANDBOX_TAKESTEP(JP)

!
! AMBER STEP TAKING
!               
               ELSE IF (AMBERT) THEN
                  IF (.NOT.ALLOCATED(MOVABLEATOMLIST)) THEN 
                     ALLOCATE(MOVABLEATOMLIST(NATOMS))
                     NMOVABLEATOMS=NATOMS
                  ENDIF
                  DISTOK=.FALSE.
                  LOCALCOUNT=0
                  DO WHILE (.NOT.DISTOK)
                     LOCALCOUNT=LOCALCOUNT+1
                     IF (LIGMOVET.AND.MOD(J1,ligmovefreq).EQ.0) THEN
                        doligmove=.TRUE.
                     ENDIF
                     
                     IF (DEBUG) THEN
                        WRITE(MYUNIT, '(A)') '=== PRE-TAKESTEP COORDS ENERGY ==='
                        CALL POTENTIAL(COORDS(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                        WRITE(MYUNIT, '(A)') '=== PRE-TAKESTEP COORDSO ENERGY ==='
                        CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL

                        OPEN(1473, FILE="prestepcoordsfile", POSITION="APPEND", STATUS="unknown", form="formatted")

                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A, I6)') 'Step ', NQ(JP)
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A)') 'coordso'
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(3F20.10)') COORDSO(:,JP)
                        CLOSE(1473)
                     END IF
!
!  seed the random number generator with system time + MYNODE (for MPI runs)
!  sf344> moved this routine before calling TAKESTEPAMBER
                     IF (RANDOMSEEDT) THEN
                        CALL DATE_AND_TIME(DATECHAR,TIMECHAR,ZONECHAR,VALUES)
                        ITIME1= VALUES(6)*60 + VALUES(7)
                        CALL SDPRND(ITIME1+MYNODE)
                     ENDIF

                     ! dc550
                     IF (SELECTMOVET) THEN
                        CALL TAKESTEPSELECTMOVE(JP)
                     ENDIF

                     CALL TAKESTEPAMBER(JP,COORDS(:,JP),MOVABLEATOMLIST,NMOVABLEATOMS,LIGMOVET,MDSTEPT,RANDOMSEEDT, 
     &                                  BLOCKMOVET,NBLOCKS,ATOMSINBLOCK)


                     IF (DEBUG) THEN
                        OPEN(1483, FILE="poststepcoordsfile", POSITION="APPEND", STATUS="unknown", form="formatted")
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A, I6)') 'Step ', NQ(JP)
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A)') 'coordso'
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(3F20.10)') COORDSO(:,JP)
                        CLOSE(1483)                    
 
                     
                        WRITE(MYUNIT, '(A)') '=== POST-TAKESTEP COORDS ENERGY ==='
                        CALL POTENTIAL(COORDS(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                        WRITE(MYUNIT, '(A)') '=== POST-TAKESTEP COORDSO ENERGY ==='
                        CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                     END IF
                     
                     doligmove=.FALSE.

                     IF (AMCHPMAX.EQ.0) THEN
                        CALL TAKESTEP(JP)
                     ELSE
!     msb50> New AMBER dihedral move routine
                        CALL TAKESTEPAMM(COORDS(:,JP), DEBUG, STEP(JP))
                     ENDIF
! Do group rotation moves
                     IF (GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP)
                     IF (DIHEDRALROTT.AND.DODIHEDRALROT) CALL DIHEDRALROTSTEP(DUMMY1,JP,DUMMY2,DUMMY3)
! hk286 - Do damped group moves
                     IF (DAMPEDGMOVET.AND.DODGROUPMOVET) CALL DAMPEDMOVESTEP(JP) 
! Rotamer moves
                     IF (ROTAMER_MOVET) CALL ROTAMER_STEP(COORDS(:, JP))

!     csw34> Check distances for groups defined in movableatoms file
!            If A->B > ABTHRESH or A->C > ACTHRESH, the step is discarded and we try again 
                     IF (LOCALSAMPLET) THEN
                        CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
                        IF (.NOT.DISTOK) COORDS(:,JP)=SAVECOORDS(:)
                     ELSE
                        DISTOK=.TRUE.
                     ENDIF
                  ENDDO
! AMBER12 step taking
               ELSE IF (AMBER12T) THEN
! Do group rotation moves
                  CALL CARTESIAN_SPHERE(COORDS(:,JP), STEP(JP))
                  IF (GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP)
                  IF (DIHEDRALROTT.AND.DODIHEDRALROT) CALL DIHEDRALROTSTEP(DUMMY1,JP,DUMMY2,DUMMY3)
                  IF (ROTAMER_MOVET) CALL ROTAMER_STEP(COORDS(:, JP))
! jdf43> MWFILM
               ELSE IF ((MWFILMT).AND.(PERIODIC)) THEN
                  CALL MWSTEP(COORDS,JP,NPAR,NATOMS,STEP(JP),BOXLX,BOXLY,BOXLZ,LAT)
! vr274> DMACRYS steps taking
               ELSE IF (DMACRYST) THEN
                  call DMACRYS_TAKESTEP(JP)
!ds656> Enumerate permutations
               ELSE IF (ENPERMST) THEN
                  CALL ENPERMS_TAKESTEP(JP,COMPLETE)
                  !write(*,*) COMPLETE
               ELSE IF (MBPOLT) THEN
                  CALL MBPOLSTEP(COORDS,STEP(JP),OSTEP(JP),TMOVE(JP),OMOVE(JP))
!
! Grand canonical size change.
!
               ELSE IF (GCBHT.AND.(MOD(J1,GCINT).EQ.0)) THEN
                  RANDOM=DPRAND()
                  CHANGED=.FALSE.
                  IF (RANDOM.LT.GCPLUS) THEN
                     IF (NATOMS+1.LE.GCNATOMS) THEN
                        CMMAX=-1.0D0
                        IF (RIGID) NATOMS=NATOMS/2
!
! Find centre of coordinates 
!
                        XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
                        DO J2=1,NATOMS
                           XMASS=XMASS+COORDS(3*(J2-1)+1,JP)
                           YMASS=YMASS+COORDS(3*(J2-1)+2,JP)
                           ZMASS=ZMASS+COORDS(3*(J2-1)+3,JP)
                        ENDDO
                        XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS

                        DO J3=1,NATOMS-NCORE(JP)
                           J2=3*J3
                           DUMMY1=SQRT((COORDS(J2-2,JP)-XMASS)**2+(COORDS(J2-1,JP)-YMASS)**2+(COORDS(J2,JP)-ZMASS)**2)
                           IF (DUMMY1.GT.CMMAX) CMMAX=DUMMY1
                        ENDDO
                        CMMAX=CMMAX+SQRT(PERCCUT) ! add the percolation distance for random atom addition radial distance
                        CHANGED=.TRUE.
                        NATOMS=NATOMS+1
                        SYMFCTR=SYMFCTR+LOG(1.0D0*NATOMS)
                        IF (MOD(J1-1,PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,G20.10)') 
     &  'mc> Attempting to add an extra atom at radial distance from centre of coordinates ',CMMAX
!
! To generate uniformly distributed random points on the (n âˆ’ 1)-sphere (i.e., 
! the surface of the n-ball), Marsaglia (1972) gives the following algorithm.
! Generate an n-dimensional vector of normal deviates (it suffices to use N(0, 1), although in fact the choice of the variance is arbitrary),
! Now calculate the "radius" of this point. ï¿¼
! The vector is uniformly distributed over the surface of the unit n-ball.
!
! Marsaglia polar method is a modification of the BoxâMuller method algorithm, 
! which does not require computation of functions sin() and cos(). 
! In this method U and V are drawn from the uniform (â1,1) distribution, 
! and then S = U2 + V2 is computed. If S is greater or equal to one then the method starts over, otherwise two quantities
! are returned. Again, X and Y will be independent and standard normally distributed.
!
! The BoxâMuller method uses two independent random numbers U and V distributed uniformly on (0,1). Then the two random variables X and Y
! will both have the standard normal distribution, and will be independent. This formulation arises because for a bivariate normal random 
! vector (X Y) the squared norm X^2 + Y^2 will have the chi-squared distribution with two degrees of freedom, which is an easily generated 
! exponential random variable corresponding to the quantity â2ln(U) in these equations; and the angle is distributed uniformly around the circle, 
! chosen by the random variable V.
!
567                     DUMMY1=DPRAND()
                        DUMMY2=DPRAND()
                        COORDS(3*NATOMS-2,JP)=SQRT(-2.0D0*LOG(DUMMY1))*COS(6.2831853071796D0*DUMMY2)
                        COORDS(3*NATOMS-1,JP)=SQRT(-2.0D0*LOG(DUMMY1))*SIN(6.2831853071796D0*DUMMY2)
                        DUMMY1=DPRAND()
                        DUMMY2=DPRAND()
                        COORDS(3*NATOMS,JP)=SQRT(-2.0D0*LOG(DUMMY1))*COS(6.2831853071796D0*DUMMY2)
                        DUMMY1=SQRT(COORDS(3*NATOMS-2,JP)**2+COORDS(3*NATOMS-1,JP)**2+COORDS(3*NATOMS,JP)**2)
                        COORDS(3*NATOMS-2,JP)=COORDS(3*NATOMS-2,JP)*CMMAX/DUMMY1
                        COORDS(3*NATOMS-1,JP)=COORDS(3*NATOMS-1,JP)*CMMAX/DUMMY1
                        COORDS(3*NATOMS,JP)=COORDS(3*NATOMS,JP)*CMMAX/DUMMY1
                        IF (DEBUG) WRITE(MYUNIT,'(A,3F20.10)') 'mc> Adding atom at ',COORDS(3*NATOMS-2:3*NATOMS,JP)
                     ELSE
                        CALL TAKESTEP(JP) ! normal step
                     ENDIF
                  ELSE
!
! For a sorted list this should do what we want - remove the most weakly bound atom.
!
                     IF (MOD(J1-1,PRTFRQ).EQ.0) WRITE(MYUNIT,'(A)') 'mc> Attempting to remove the most weakly bound atom'
                     NATOMS=NATOMS-1
                     SYMFCTR=SYMFCTR-LOG(1.0D0*NATOMS)
                     IF (NATOMS.GT.1) THEN
                        CHANGED=.TRUE.
                     ELSE
                        NATOMS=NATOMS+1
                        SYMFCTR=SYMFCTR+LOG(1.0D0*NATOMS)
                        CALL TAKESTEP(JP) ! normal step
                     ENDIF
                  ENDIF
                  IF (MOD(J1,PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,I6,A)') 'mc> Grand canonical size is now',NATOMS,' atoms'

                  IF (CHANGED) THEN
                     CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
                     IF (RMS.GT.BQMAX) THEN
                        WRITE(MYUNIT,'(A)') 'mc> Quench from new size geometry failed - reject GC step attempt'
                        WRITE(MYUNIT,'(A,2G20.10)') 'mc> RMS and BQMAX=',RMS,BQMAX
                        CHANGED=.FALSE.
                        IF (NATOMS.GT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR-LOG(1.0D0*NATOMS)
                        IF (NATOMS.LT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR+LOG(1.0D0*GCNATOMSBEST(JP))
                        NATOMS=GCNATOMSBEST(JP)
                        COORDS(1:3*NATOMS,JP)=COORDSO(1:3*NATOMS,JP)
                        VAT(1:NATOMS,JP)=VATO(1:NATOMS,JP)
                        POTEL=EPREV(JP)
                        WRITE(MYUNIT,'(A,I6)') 'mc> Number of atoms reset to ',NATOMS
                     ENDIF

                     IF (PERCOLATET) THEN
                        IF (.NOT. PERCT) THEN
                           WRITE(MYUNIT,'(A)') 'mc> Disconnected atoms, rejecting step. '
                           CHANGED=.FALSE.
                           IF (NATOMS.GT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR-LOG(1.0D0*NATOMS)
                           IF (NATOMS.LT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR+LOG(1.0D0*GCNATOMSBEST(JP))
                           NATOMS=GCNATOMSBEST(JP)
                           COORDS(1:3*NATOMS,JP)=COORDSO(1:3*NATOMS,JP)
                           VAT(1:NATOMS,JP)=VATO(1:NATOMS,JP)
                           POTEL=EPREV(JP)
                           WRITE(MYUNIT,'(A,I6)') 'mc> Number of atoms reset to ',NATOMS
                        ENDIF
                     ENDIF
                     IF (MOD(J1-1,PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,G20.10)') 'mc> New quench energy=',POTEL
                  ENDIF

                  IF (CHANGED) THEN
!
! Need to save best values for previous block for later block accept/reject test
! We may need to restore them.
!
                     GCEBESTP(JP)=GCEBEST(JP)
                     GCJBESTP(JP)=GCJBEST(JP)
                     GCBESTCOORDSP(1:3*GCNATOMSBEST(JP),JP)=GCBESTCOORDS(1:3*GCNATOMSBEST(JP),JP)
                     GCBESTVATP(1:GCNATOMSBEST(JP),JP)=GCBESTVAT(1:GCNATOMSBEST(JP),JP)
                     GCNATOMSBESTP(JP)=GCNATOMSBEST(JP)
!
! Need to minimise new configuration as in newrestart and save new energy, coords, etc.
!
                     QGCBH=0 ! reset relaxation block start
!
! need to ensure step is accepted
!
                     NSUCCESS(JP)=0
                     NFAIL(JP)=0
                     EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
                     BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                     JBEST(JP)=J1
                     EPREV(JP)=POTEL
                     EPPREV(JP)=0.0D0
                     NSYMREM=0
                     IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
                     COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
                     VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

                     GCEBEST(JP)=POTEL
                     GCJBEST(JP)=J1
                     GCBESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                     GCBESTVAT(1:NATOMS,JP)=VAT(1:NATOMS,JP)
                     GCNATOMSBEST(JP)=NATOMS
                  ENDIF
!
! ALL OTHER STEP TAKING
!
               ELSE
 
!  These coordinates are overwritten if we try to call MPI_IPROBE
!
!                 WRITE(MYUNIT,'(I6)') NATOMS
!                 WRITE(MYUNIT,'(A,G20.10)') 'eprev=',EPREV
!                 WRITE(MYUNIT,'(A,3G20.10)') ('LA ',COORDS(3*(J3-1)+1:3*(J3-1)+3,JP),J3=1,NATOMS)
                  IF (USERPOTT.AND.GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP)
!ab2111> dihedralrotation with user-defined potential
                  IF (USERPOTT.AND.DIHEDRALROTT.AND.DODIHEDRALROT) CALL DIHEDRALROTSTEP(DUMMY1,JP,DUMMY2,DUMMY3)
                  CALL TAKESTEP(JP)
!ds656>
                  IF (RANDPERMT .AND. MOD(J1,RANDPERM_STEP)==0) CALL RANDPERM(JP)
                  IF (RANDMULTIPERMT.AND.MOD(J1,RANDMULTIPERM_STEP)==0) CALL RANDMULTIPERM(JP,1)

               ENDIF
! jdf43> 
               IF (DDMT) THEN
                  DO I2=1,3*NATOMS
                     CALL DDMCONDENSE(COORDS)
                  ENDDO
               ENDIF
! Restore atom coordinates if atom is FROZEN or DONTMOVE as long as
! we're not taking internal coordinate moves in CHARMM or AMBER
               IF ((CHNMAX.LE.0.0D0).OR.(AMCHNMAX.LE.0.0D0)) THEN
                  IF(FREEZE) THEN
                     DO J2=1,NATOMS
                        IF (FROZEN(J2)) THEN
! hk286 - need to remove this
                           !COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        ENDIF
                     ENDDO
                  ENDIF
!
! Same for DONTMOVE atoms
! 
                  IF(DONTMOVET) THEN
                     DO J2=1,NATOMS
                        IF (DONTMOVE(J2)) THEN
                           COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        ENDIF
                     ENDDO
                  ENDIF
! js850> Same for HARMONICFLIST atoms if HARMONICDONTMOVE is specified
                  IF (HARMONICDONTMOVE) THEN
                     DO J2=1,NATOMS
                        IF (HARMONICFLIST(J2)) THEN
                           COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
               CALL A9DUMPPDB(COORDS(:,1),'afterstep')
! csw34> Redefine rigid bodies if they are being changed by steps
               IF (RIGIDINIT.AND.UPDATERIGIDREFT) THEN
                  CALL GENRIGID_UPDATE_REFERENCE(COORDS(:,JP))
               ENDIF 
 
!
! Reset switch variables for steps not done every time
!
               DOGROUPROT=.FALSE.
               DODGROUPMOVET=.FALSE.
               DOEXPANDRIGID=.FALSE.
               DOROTATERIGID=.FALSE.
               DOTRANSLATERIGID=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF STEP TAKING CALLS!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! csw34> end of temporary IF block for new moves - look for NEWMOVEST
            ENDIF
            
! ds656> Translate centre of geometry (COG) if needed
            IF(BOXCENTROIDT) CALL BOXCENTROID(COORDS(:,JP))

! csw34> If you want to look at the effect of moves, you can dump out
! the structure AFTER the move here.
!              CALL A9DUMPPDB(COORDS(:,JP),"aftermove")
!              CALL CHARMMDUMP(COORDS(:,JP),"aftermove")
            
! KEYWORD <DUMPSTEPS> BLOCK            
! csw34> Dump the coordinates after every step in AMBER pdb and rst format 
                  IF (DUMPSTEPST) THEN
                      WRITE(QUENCHNUM,*) NQ(JP)
                      IF(DMACRYST) THEN
                          CALL DMACRYS_DUMP('afterstep'//TRIM(ADJUSTL(QUENCHNUM))//'.cif', COORDS(:,JP))
                      ELSE IF(USERPOTT) THEN
                          CALL USERPOT_DUMP_CONFIGURATION('afterstep'//TRIM(ADJUSTL(QUENCHNUM))//'.dat', COORDS(:,JP))
                      ELSE IF (AMBERT) THEN
                          QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                          OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')
                          WRITE(20,'(a20)') QUNAME
                          WRITE(20,'(i5)') NATOMS
                          WRITE(20,'(6f12.7)') COORDS(:,JP)
                          CLOSE(20)
! csw34> Dump to PDB using routine in amberinterface.f
                          QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))
                          CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                      ELSE IF (AMBER12T) THEN
                          QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))
                          CALL AMBER12_WRITE_RESTART(COORDS(:,JP), TRIM(ADJUSTL(QUNAME))//'.rst', 
     &                                                         LEN(TRIM(ADJUSTL(QUNAME))//'.rst'))
                          CALL AMBER12_WRITE_PDB(COORDS(:,JP), TRIM(ADJUSTL(QUNAME))//'.pdb', 
     &                                                     LEN(TRIM(ADJUSTL(QUNAME))//'.pdb'))
                      ENDIF
                  ENDIF
! END OF KEYWORD <DUMPSTEPS> BLOCK


               NQ(JP)=NQ(JP)+1
               IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
!fh301>{{{
               STRUC = NQ(JP)
!fh301>}}}
               CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
               NQTOT=NQTOT+1

!  Check for results of taboo list. SELFT is set in taboo.
 
               IF (SELFT) THEN
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,'(A)') 'Taboo list:'
                     WRITE(MYUNIT,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
                     WRITE(MYUNIT,73) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
73                   FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' TABOO')
                  ENDIF
                  GOTO 23
               ENDIF

!ds656> print some energies to test BASWAP
               IF(BASWAPTEST) THEN
                  !CALL BGUPTA(COORDS,GRAD,EASWAPQ,.FALSE.)
                  IF(NSWAPS == 1) THEN
                     WRITE(MYUNIT, '(2(A,1X,E12.6,1X))') 
     1                    "mc> 1 swap: de_A=",DE1,"and de_B=",DE2
                  ENDIF
                  WRITE(MYUNIT, '(2(A,1X,E12.6,1X),A)') 
     2                 "mc> BA swaps: dE=", EASWAP-EPREV(JP), "before and dE=",
     3                 POTEL-EPREV(JP), "after quench."
                  !WRITE(MYUNIT, '(A,3(1X,F15.8))') "mc> ds656:", EPREV(JP)-EBSWAP, EASWAP, EASWAPQ-POTEL 
                  NSWAPS=0
               ENDIF               
!<ds656        
               
!
!  Output
!
!fh301>{{{
               IF (CHEMSHIFT2) THEN
                  CHEMSHIFT=.TRUE.
!     CHEMSHIFTITER=0
               ENDIF
               IF (MOD(J1-1,PRTFRQ).EQ.0) THEN !<start of looong block of print statements
                  IF (CHEMSHIFT) THEN
                     IF (NPAR.GT.1) THEN
                        WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1                       POTEL,' ECAM=', ENERGYCAMSHIFT,' steps=',ITERATIONS,' RMS=',RMS,
     1                       ' Markov E=',EPREV(JP),' t=',TIME-TSTART
                     ELSE
                        WRITE(MYUNIT,'(A,I10,A,G20.10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1                       POTEL,' ECAM=', ENERGYCAMSHIFT,' steps=',ITERATIONS,' RMS=',RMS,
     1                       ' Markov E=',EPREV(JP),' t=',TIME-TSTART
                     ENDIF
                  ELSE
!     fh301>}}}
                     IF (NPAR.GT.1) THEN
                        IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,
     1                       ']Qu ',NQ(JP),' E=',POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
                     ELSE
                        IF (.NOT.SUPPRESST) WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1                       POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
!fh301>{{{
                     ENDIF
!fh301>}}}
                  ENDIF
                  
!                       
! khs26> Printing energy decomposition for AMBER energies if we're using AMBER.
                  IF (ENERGY_DECOMPT) THEN
                     IF (AMBERT) THEN
                        WRITE(MYUNIT,'(A10,F20.10,A10,F20.10,A10,F20.10)') 'bond:     ', E_BOND, ' angle:   ',  
     1                       E_ANGLE, ' dihed:   ', E_DIHEDRAL
                        WRITE(MYUNIT,'(A10,F20.10,A10,F20.10)') 'vdW:      ', E_VDW,  ' 14 vdW:  ', E_14_VDW
                        WRITE(MYUNIT,'(A10,F20.10,A10,F20.10)') 'elec:     ', E_ELEC, ' 14 elec: ', E_14_ELEC
!     khs26> igb = 0 is no implicit solvation
                        IF (IGB.EQ.0) THEN
                           WRITE(MYUNIT,'(A10,F20.10)') 'hbond:    ', E_IGB
!     khs26> igb = 10 is the Poisson-Boltzmann case
                        ELSE IF (IGB.EQ.10) THEN
                           WRITE(MYUNIT,'(A10,F20.10)') 'ipb:      ', E_IGB
!     khs26> other values of igb are detailed in the AMBER manual
                        ELSE
                           WRITE(MYUNIT,'(A10,F20.10)') 'igb:      ', E_IGB
                        END IF
                     ELSE IF (AMBER12T) THEN
! khs26> Add the energy decomposition terms from AMBER12 here! I need to think of how to get access to it from
! potential.f without adding it to commons (eugh).
                        WRITE(MYUNIT,'(A80)') 'Energy decomposition not yet implemented for AMB12GMIN.' 
                     END IF
                  END IF
               ENDIF !<End of looong block of print statements
               
!     mp466>  writes structure and energetic data at regular increments
!             to *plot and movie files for AMH potential
!     
               IF ((MOD(J1,NINT_AMH).EQ.0).AND.AMHT) THEN
                  
                  GLY_COUNT = 0
                  SNAPCOUNT = SNAPCOUNT + 1
                  
                  DO III = 1,NMRES
                     IF (IRES(III).EQ.8) THEN
                        PRCORD(III,1,1,1) = COORDS(9*(III-1)+1- GLY_COUNT*3,JP) !  CA X
                        PRCORD(III,2,1,1) = COORDS(9*(III-1)+2- GLY_COUNT*3,JP) !  CA Y
                        PRCORD(III,3,1,1) = COORDS(9*(III-1)+3- GLY_COUNT*3,JP) !  CA Z
!     SWAP  CA for CB
                        PRCORD(III,1,1,2) = COORDS(9*(III-1)+1- GLY_COUNT*3,JP) !  CB X
                        PRCORD(III,2,1,2) = COORDS(9*(III-1)+2- GLY_COUNT*3,JP) !  CB Y
                        PRCORD(III,3,1,2) = COORDS(9*(III-1)+3- GLY_COUNT*3,JP) !  CB Z
                        PRCORD(III,1,1,3) = COORDS(9*(III-1)+4- GLY_COUNT*3,JP) !  O X
                        PRCORD(III,2,1,3) = COORDS(9*(III-1)+5- GLY_COUNT*3,JP) !  O Y
                        PRCORD(III,3,1,3) = COORDS(9*(III-1)+6- GLY_COUNT*3,JP) !  O Z
                        GLY_COUNT = GLY_COUNT +1
                     ELSE
                        PRCORD(III,1,1,1) = COORDS(9*(III-1)+1 - GLY_COUNT*3,JP) !  CA X
                        PRCORD(III,2,1,1) = COORDS(9*(III-1)+2 - GLY_COUNT*3,JP) !  CA Y
                        PRCORD(III,3,1,1) = COORDS(9*(III-1)+3 - GLY_COUNT*3,JP) !  CA Z
                        PRCORD(III,1,1,2) = COORDS(9*(III-1)+4 - GLY_COUNT*3,JP) !  CB X
                        PRCORD(III,2,1,2) = COORDS(9*(III-1)+5 - GLY_COUNT*3,JP) !  CB Y
                        PRCORD(III,3,1,2) = COORDS(9*(III-1)+6 - GLY_COUNT*3,JP) !  CB Z
                        PRCORD(III,1,1,3) = COORDS(9*(III-1)+7 - GLY_COUNT*3,JP) !  O X
                        PRCORD(III,2,1,3) = COORDS(9*(III-1)+8 - GLY_COUNT*3,JP) !  O Y
                        PRCORD(III,3,1,3) = COORDS(9*(III-1)+9 - GLY_COUNT*3,JP) !  O Z
                     ENDIF
                  ENDDO
                  
                  WRITE(OMOVI,683)1,SNAPCOUNT,1,TEMP(JP),1
 683              FORMAT(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')
                  
                  DO  I500=1,NMRES
                     WRITE(OMOVI,332)(PRCORD(I500,I2,1,1),I2=1,3),(PRCORD(I500,I2,1,2),I2=1,3),(PRCORD(I500,I2,1,3),I2=1,3)
                  ENDDO
 332              FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1x))

                  CALL E_WRITE(AVEP(:,:,:),TEMP(JP),NUMPRO,SNAPCOUNT)
               ENDIF            !                  IF (MOD(J1,NINT_AMH).EQ.0)
               
               CALL FLUSH(MYUNIT)
               
! ds656> Homotop refinement 
! Don't bother if last quench failed to converge,
! which will happen when SYMMETRY puts two atoms
! in one symmetry site. In this case QDONE=0.
!     
               IF ( LSWAPST .AND.
     1              (MOD(J1,LSWAPS_N1)==0) ) THEN
                  CALL USWAPS(JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
               ENDIF

               IF ( LFLIPST .AND.
     1              (MOD(J1,LFLIPS_N1)==0) ) THEN
                  CALL FLIPSEQ(JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
               ENDIF

               IF ((QDONE.EQ.1).AND.(.NOT.SYMMETRIZE)) THEN
                  IF (HOMOREF_BHT) CALL HOMOREF_BH(JP, ITERATIONS, TIME, BRUN, QDONE, SCREENC)
                  IF (HOMOREFT) CALL HOMOREF(JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)                 
                  IF (QALCST.OR.QALCS_SURFT.OR.QALCS_SYMT) CALL QALCS(JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
               ENDIF
!<ds656
               
!  RMS compared to reference structure 'compare'
          IF (RMST.AND.CHRMMT) THEN
             CALL CHRMS(JP,RMSD)
             IF (DEBUG) WRITE(MYUNIT,'(A,F15.5)')'RMSD = ',RMSD
             IF (RMSD.LE.RMSLIMIT) CALL SAVERMS(JP,POTEL,RMSD)
!                NRMS=NRMS+1
!                WRITE(CNRMS,'(I6)') NRMS
!                CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'rms.'//TRIM(ADJUSTL(CNRMS)))
!                OPEN(UNIT=20,FILE='rms.'//TRIM(ADJUSTL(CNRMS)),POSITION='APPEND',STATUS='OLD')
!                WRITE(20,'(A,I6,A,F15.5,A,F15.5)') '*   Qu ',NQ(JP),' E=',POTEL,' RMSD=',RMSD
!                CLOSE(20)
!            ENDIF
          ENDIF

!
!     csw34> PAIRDIST prints the distance between specified atom pairs
!     after each quench. 
!
          IF (PAIRDISTT) THEN
!     Write end of previous line as using ADVANCE="NO"
                WRITE(MYPUNIT,*) " "
!     Write current quench numberWRITE(97,'(A,I0.3,A,I5,A)')"found in process
!     ",RANK," after ",I," attempts."
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ(JP)
!     For each pair, assign ATOM1 and ATOM2 arrays containing coordinates
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1),JP)
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2),JP)
!     Call PAIRDISTANCE with (x,y,z) for each atom
                WRITE(MYPUNIT,'(F10.4)',ADVANCE="NO") PAIRDISTANCE(ATOM1,ATOM2) 
             ENDDO
             CALL FLUSH(MYPUNIT)
          ENDIF
          
!
!     csw34> TRACKDATA keyword prints the quench energy, markov energy
!     and energy of the current lowest minimum to files for viewing during a run. 
!
          IF (TRACKDATAT) THEN
             WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV(JP)
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             CALL FLUSH(MYEUNIT)
             CALL FLUSH(MYMUNIT)
             CALL FLUSH(MYBUNIT)
             IF (A9INTET.AND.(NQ(JP).GT.2)) THEN
                WRITE(3999,'(I10,G20.10)') NQ(JP),INTEQMIN(1)
                CALL FLUSH(3999)
             ENDIF
!
!     csw34> if RMS is also specified, prints the RMSD from the
!     comparison structure to the file 'rmsd'. Which RMSD depends on the
!     final arguement of the RMS keyword - see documentation!

             IF (RMST.AND.CHRMMT) THEN
                WRITE(MYRUNIT,'(I10,F15.5)') NQ(JP),RMSD
                CALL FLUSH(MYRUNIT)
             ENDIF
          ENDIF

! DAESTAT keyword prints two files :
! stat.all which contains all quenches and their energies
! stat.acc which contains only accepted quenches and their energies
! used to analyse how many new minima are being found
 
            IF (DAESTAT) THEN
               PRINT*,'DAESTAT block in mc.f not implemented'
               STOP
!              CALL CALCMIND(JP,MIND)
!              CALL CALCDIHE(DIHE)
!              WRITE(MYUNIT,'(A,I6,3F20.10)') 'NQALL POTEL',NQ(JP),POTEL,MIND,DIHE
!              WRITE(36,'(I6,3F20.10)') NQ(JP),POTEL,MIND,DIHE
            ENDIF

! js850>
          IF (DUMPUNIQUE) THEN
            CALL DUMPUNIQUE_DUMP(DUMPUNIQUEUNIT, JP, POTEL, DUMPUNIQUEEPREV )
          ENDIF
 
!  Check for reseeding.
 
            IF (NEWRESTART.AND.(.NOT.SEEDT)) THEN
              CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,RMIN,RVAT,BRUN,SCREENC,QDONE,
     &                    JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST,RES1,RES2)
              IF (RES1.OR.RES2) QNEWRES=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
!             IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
!             COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!             WRITE(MYUNIT,'(A,2G20.10)'),'mc> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)     
!             VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
            ENDIF
            IF (STAY) THEN
            ELSE IF (EVAP .and. .not.evapreject) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL MYRESET(JP,NATOMS,NPAR,NSEED)
               IF (DEBUG) THEN
                  WRITE(MYUNIT,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
!     csw34> A series of tests start here to check if a structure should
!     be allowed into the markov chain. If it fails a test, the ATEST
!     variable will be set to .FALSE. 
               ATEST=.TRUE. 
!     csw34> Test for cold fusion in mylbfgs.f - if so, set ATEST to .FALSE.
               IF (COLDFUSION) THEN
                  ATEST=.FALSE.
               ENDIF
               COLDFUSION=.FALSE.

               IF (GRADPROBLEMT) THEN
                  ATEST=.FALSE.
               ENDIF

!     jwrm2> If using percolation, and it's not a connected network, reject the step, unless no
!            step has yet been connected.
!            If this is the first step that's fully connected, set PERCACCEPTED to true.
               IF (PERCOLATET) THEN
                 IF(PERCACCEPTED .AND. (.NOT. PERCT)) THEN
                   WRITE(MYUNIT,*) 'mc> Disconnected atoms, rejecting step. '
                   IF (DEBUG) WRITE(MYUNIT,'(3G20.10)') COORDS(1:3*NATOMS,JP)
                   ATEST=.FALSE.
                 ELSE IF((.NOT. PERCACCEPTED) .AND. PERCT) THEN
                   PERCACCEPTED=.TRUE.
                   WRITE(MYUNIT,*) 'perc> Found first connected configuration'
                 ENDIF
               ENDIF

!     csw34> If the chirality or peptide bond checks in quench.f have
!            failed, GOODSTRUCTURE should be .FALSE. and we should set ATEST to
!            .FALSE. to prevent these bad structures entering the markov chain!
               IF (.NOT.GOODSTRUCTURE) ATEST=.FALSE.

! djw34> Only for AMBER9, AMBER12 or CUDA with AMBER!
               IF (AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))) THEN
! HBONDMATRIX
                  IF (HBONDMATRIX.AND.ATEST) THEN
! Save the quench energy for later analysis
                     HBONDQUE(NQ(JP))=POTEL
! Dump the coordinates in .rst format to hbond.rst
                     CALL DUMPRST(COORDS(:,JP),NATOMS,'hbond')
! Set cutoffs to default values
                     HBONDDCUT='3.00'
                     HBONDACUT='120.0'
! Calculate hydrogen-bond matrix
                     CALL GETHBONDMATRIX('coords.prmtop','hbond.rst',TRIM(ADJUSTL(HBONDDONORSACCEPTORS)), 
     &                                   TRIM(ADJUSTL(HBONDRESIDUES)),TRIM(ADJUSTL(HBONDDCUT)),
     &                                   TRIM(ADJUSTL(HBONDACUT)),HBONDMAT)
! Check to see if this structure is part of an existing group 
                     HBONDMATCH=.FALSE.
                     HBONDMATCHN=0
                     DO L=1,NHBONDGROUPS
                        IF(HBONDLIGAND) THEN
                           HBONDDIFF=SUMSQUAREDIFF1D(HBONDMAT(HBONDLIGANDN,:),HBONDGROUPS(L,HBONDLIGANDN,:),HBONDNRES)
                        ELSE
                           HBONDDIFF=SUMSQUAREDIFF(HBONDMAT,HBONDGROUPS(L,:,:),HBONDNRES)
                        ENDIF
                        IF(HBONDDIFF.EQ.0) THEN
! If SUMSQUAREDIFF returns 0, the matricies are identical 
                           HBONDMATCH=.TRUE.
                           HBONDMATCHN=L
                           WRITE(MYUNIT,'(A,I4)') " HBONDMATRIX> Structure matches group ",HBONDMATCHN
! Increment the group population
                           HBONDGROUPPOP(HBONDMATCHN)=HBONDGROUPPOP(HBONDMATCHN)+1
                           HBONDGROUPIDS(NQ(JP))=HBONDMATCHN
                           WRITE(MYUNIT,'(A,I4,A,I5)') " HBONDMATRIX> Group ",HBONDMATCHN,
     &                                                 " now has population ",HBONDGROUPPOP(HBONDMATCHN)
                           EXIT
                        ENDIF
                     ENDDO
! If no match with ptraj defaults, use soft cutoff and try again
! Looser
                     IF(.NOT.HBONDMATCH) THEN
                        HBONDDCUT=HBONDDCUTLOOSE
                        HBONDACUT=HBONDACUTLOOSE
! Recalculate hydrogen-bond matrix using looser cutoffs
                        CALL GETHBONDMATRIX('coords.prmtop','hbond.rst',TRIM(ADJUSTL(HBONDDONORSACCEPTORS)), 
     &                                      TRIM(ADJUSTL(HBONDRESIDUES)),TRIM(ADJUSTL(HBONDDCUT)), 
     &                                      TRIM(ADJUSTL(HBONDACUT)),HBONDSOFTMAT)
! Recheck all groups with the new matrix HBONDSOFTMAT
                        HBONDMATCHN=0
                        HBONDDIFF=1
                        DO L=1,NHBONDGROUPS
                           IF(HBONDLIGAND) THEN 
                              HBONDDIFF=SUMSQUAREDIFF1D(HBONDSOFTMAT(HBONDLIGANDN,:),HBONDGROUPS(L,HBONDLIGANDN,:),HBONDNRES)
                           ELSE
                              HBONDDIFF=SUMSQUAREDIFF(HBONDSOFTMAT,HBONDGROUPS(L,:,:),HBONDNRES)
                           ENDIF
                           IF(HBONDDIFF.EQ.0) THEN
! If SUMSQUAREDIFF returns 0, the matricies are identical 
                              HBONDMATCH=.TRUE.
                              HBONDMATCHN=L
                              WRITE(MYUNIT,'(A,I4)') 
     &              " HBONDMATRIX> Using looser hydrogen-bond cutoffs, structure matches group ",HBONDMATCHN
! Increment the group population
                              HBONDGROUPPOP(HBONDMATCHN)=HBONDGROUPPOP(HBONDMATCHN)+1
                              HBONDGROUPIDS(NQ(JP))=HBONDMATCHN
                              WRITE(MYUNIT,'(A,I4,A,I5)') " HBONDMATRIX> Group ",HBONDMATCHN,
     &                                                 " now has population ",HBONDGROUPPOP(HBONDMATCHN)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
! Tighter
                     IF(.NOT.HBONDMATCH) THEN
                        HBONDDCUT=HBONDDCUTTIGHT
                        HBONDACUT=HBONDACUTTIGHT
! Recalculate hydrogen-bond matrix using tighter cutoffs
                        CALL GETHBONDMATRIX('coords.prmtop','hbond.rst',TRIM(ADJUSTL(HBONDDONORSACCEPTORS)), 
     &                                      TRIM(ADJUSTL(HBONDRESIDUES)),TRIM(ADJUSTL(HBONDDCUT)), 
     &                                      TRIM(ADJUSTL(HBONDACUT)),HBONDSOFTMAT)
! Recheck all groups with the new matrix HBONDSOFTMAT
                        HBONDMATCHN=0
                        DO L=1,NHBONDGROUPS
                           IF(HBONDLIGAND) THEN 
                              HBONDDIFF=SUMSQUAREDIFF1D(HBONDSOFTMAT(HBONDLIGANDN,:),HBONDGROUPS(L,HBONDLIGANDN,:),HBONDNRES)
                           ELSE
                              HBONDDIFF=SUMSQUAREDIFF(HBONDSOFTMAT,HBONDGROUPS(L,:,:),HBONDNRES)
                           ENDIF
                           IF(HBONDDIFF.EQ.0) THEN
! If SUMSQUAREDIFF returns 0, the matricies are identical 
                              HBONDMATCH=.TRUE.
                              HBONDMATCHN=L
                              WRITE(MYUNIT,'(A,I4)') 
     &              " HBONDMATRIX> Using tighter hydrogen-bond cutoffs, structure matches group ",HBONDMATCHN
! Increment the group population
                              HBONDGROUPPOP(HBONDMATCHN)=HBONDGROUPPOP(HBONDMATCHN)+1
                              HBONDGROUPIDS(NQ(JP))=HBONDMATCHN
                              WRITE(MYUNIT,'(A,I4,A,I5)') " HBONDMATRIX> Group ",HBONDMATCHN,
     &                                                 " now has population ",HBONDGROUPPOP(HBONDMATCHN)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
! If no match was found, even with the soft cutoff, create a new group
                     IF(.NOT.HBONDMATCH) THEN
                        NHBONDGROUPS=NHBONDGROUPS+1
! Save the group matrix
                        HBONDGROUPS(NHBONDGROUPS,:,:)=HBONDMAT
! Initialize the population of the group to 1
                        HBONDGROUPPOP(NHBONDGROUPS)=1
! Add the initial structure to the group
                        HBONDGROUPIDS(NQ(JP))=NHBONDGROUPS
! Set the group max energy
                        HBONDMAXE(NHBONDGROUPS)=POTEL
! Set the group best and markov energy
                        HBONDBEST(NHBONDGROUPS)=POTEL
                        HBONDMARKOV(NHBONDGROUPS)=POTEL
! Set the best coordinates
                        HBONDBESTCOORDS(NHBONDGROUPS,:)=COORDS(:,JP)
! Print some output
                        WRITE(MYUNIT,'(A,I4)') " HBONDMATRIX> New group found, creating group ",NHBONDGROUPS
! Write out the matrix for reference
                        WRITE(NHBONDGROUPSCHAR,*) NHBONDGROUPS
                        HBONDGROUPNAME='group'//TRIM(ADJUSTL(NHBONDGROUPSCHAR))
                        CALL WRITEHBONDMATRIX(HBONDMAT(:,:),HBONDNRES,HBONDGROUPNAME)
                     ENDIF
                  ENDIF
! END OF HBONDMATRIX
!     csw34> Calculate interaction energy to residue specified in the AMBGMINintE.sh script
                  IF (A9INTET.AND.ATEST) THEN
                     CALL A9INTE(COORDS(:,JP),INTE)
                     CALL A9INTESAVEIT(INTE,COORDS(:,JP),JP)
                     CALL A9INTEDUMPLOWEST()
                     IF (DEBUG) WRITE(MYUNIT,*) 'intE=',INTE
                     IF (TRACKDATAT) THEN
                        WRITE(3998,'(I10,G20.10)') NQ(JP),INTE
                        CALL FLUSH(3998)
                     ENDIF
                  ENDIF
!     csw34> Dumping every quench in AMBER9 format - only dumps if no chiral problems found 
!     csw34> Dump the quench coordinates in AMBER pdb and rst format if DUMPQ is present and chirality checks satisfied
                  IF (DUMPQUT.AND.ATEST) THEN
                      IF (AMBERT) THEN
                         WRITE(QUENCHNUM,*) NQ(JP)
                         QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                         OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                         WRITE(20,'(a20)') QUNAME
                         WRITE(20,'(i5)') NATOMS
                         WRITE(20,'(6f12.7)') COORDS(:,JP) 
                         CLOSE(20)
!     csw34> Dump to PDB using routine in amberinterface.f
                         QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))
                         CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                      ELSEIF (AMBER12T) THEN
                         WRITE(QUENCHNUM,*) NQ(JP)
                         QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))
                         CALL AMBER12_WRITE_RESTART(COORDS(:,JP), TRIM(ADJUSTL(QUNAME))//'.rst',
     &                                                        LEN(TRIM(ADJUSTL(QUNAME))//'.rst'))
                         CALL AMBER12_WRITE_PDB(COORDS(:,JP), TRIM(ADJUSTL(QUNAME))//'.pdb',
     &                                                    LEN(TRIM(ADJUSTL(QUNAME))//'.pdb'))
                      ENDIF
                  ELSEIF (DUMPQUT) THEN
                      WRITE(*,'(A,I6)') 'DUMPQU> chirality/cis-trans problem detected, not dumping quench ',NQ(JP)
                  ENDIF
!     khs26> Dump current best structure to PDB using routine in amberinterface.f
                  IF (DUMPBESTT) THEN
                     WRITE(QUENCHNUM,*) NQ(JP)
                     BESTNAME='best_'//TRIM(ADJUSTL(QUENCHNUM))
                     IF (AMBERT) THEN
                        CALL A9DUMPPDB(QMINP(1,:), BESTNAME)
                        CALL A9DUMPPDB(QMINP(1,:), 'current_best')
                     ELSEIF (AMBER12T) THEN
                        CALL AMBER12_WRITE_PDB(QMINP(1,:), TRIM(ADJUSTL(BESTNAME))//'.pdb',
     &                                                 LEN(TRIM(ADJUSTL(BESTNAME))//'.pdb'))
                        CALL AMBER12_WRITE_PDB(QMINP(1,:), 'current_best.pdb', LEN('current_best.pdb'))
                     ENDIF
                  ENDIF
               ENDIF
!     csw34> Check to see if LOCALSAMPLE constraints have been violated - if they have, reject the step
               IF (LOCALSAMPLET) THEN
                  DISTOK=.FALSE.
                  CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
                  IF (.NOT.DISTOK) ATEST=.FALSE.
               ENDIF

!     js850> Check to see if RESTRICTREGION constraints have been violated - if they have, reject the step
               IF ( RESTRICTREGION .AND. RESTRICTREGIONTEST ) THEN
!                write(*,*) "rejecting from markov chain"
                 ATEST=.FALSE.
                 WRITE(MYUNIT,* ) 'RESTRICTREGION> rejecting quench ', NQ(JP) 
               ENDIF
!     csw34> If structure matched a previous hydrogen-bond group, we want to compare it to the last 
!            accepted energy in that group
               IF ((AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))).AND.HBONDMATRIX.AND.HBONDMATCH.AND.ATEST) THEN
                  HBONDESAVE=EPREV(JP)
                  EPREV(JP)=HBONDMARKOV(HBONDMATCHN)
                  WRITE(MYUNIT,'(A,I3,A,G20.10)') ' HBONDMATRIX> Setting Markov energy to previous group ',
     &                                             HBONDMATCHN,' value of ',HBONDMARKOV(HBONDMATCHN) 
               ENDIF
               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
!ds656> When enumerating perumations, every single step should be accepted
               ELSEIF(ENPERMST) THEN
                  ATEST=.TRUE.
               ELSEIF (ATEST.AND..NOT.RATIOT) THEN
                  CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM,MCTEMP)
!     csw34> Auto-accept if we have changed to a new hydrogen-bond group (HBONDMATCH=.FALSE.)
                  IF ((AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))).AND.HBONDMATRIX) THEN
                     IF (.NOT.HBONDMATCH) THEN
!     csw34> In the default case, we want to auto-accept the step if it finds a new hydrogen-bond group
                        IF (HBONDACCEPT) THEN
                           ATEST=.TRUE.
                           WRITE(MYUNIT,'(A)') " HBONDMATRIX> New group so step auto-accepted"
!     csw34> However - if we wanted to restrict sampling to a single (the starting) hydrogen-bond group
!            we can do this by auto-rejecting any step which moves to a new group.
                        ELSE
                           ATEST=.FALSE.
                           WRITE(MYUNIT,'(A)') " HBONDMATRIX> Group changed so step rejected!"
                        ENDIF
                     ELSE
                        EPREV(JP)=HBONDESAVE
!                        WRITE(MYUNIT,'(A)') " HBONDMATRIX> Markov energy reset to pass sanity checks"
                     ENDIF
                  ENDIF
               ENDIF
!
!  Sanity check to make sure the Markov energy agrees with COORDSO. 
!  Stop if not true.
!
               IF ((DEBUG.OR.CHECKMARKOVT).AND.(.NOT.(HBONDMATRIX.OR.PERMOPT.OR.PERMINVOPT.OR.FEBHT.OR.QALCST))) THEN
                  CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV(JP)).GT.ABS(ECONV)) THEN
                     IF (EVAP) THEN
                        WRITE(MYUNIT,'(3(A,G20.10))') 'mc> WARNING - energy for saved coordinates ',OPOTEL,
     &                     ' differs from Markov energy ',EPREV(JP),' because an atom moved outside the container'
                     ELSE
                        WRITE(MYUNIT,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,
     &                                                 ' but Markov energy=',EPREV(JP)
                        IF (.NOT.(CHEMSHIFT2)) STOP
                     ENDIF
                  ENDIF
               ENDIF
               
               IF (RATIOT) THEN
                  CALL NULLMOVE(SCREENC,COORDSO(:,JP),ATEST,NULLMOVES,JP)
               ENDIF

! Accept or reject step. If the quench did not converge then allow a
! potenial move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
! However, for TIP5P some cold fusion events that had not actually reached the threshold for
! rejection were actually accepted. Must prevent this!
               IF (ATEST) THEN
#ifdef MPI
		  IF (DEBUG) THEN
                     WRITE(MYUNIT,334) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
334                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
#else
                  IF (DEBUG) THEN
! csw34> Make sure the group Markov energy is printed if using HBONDMATRIX
                     IF (HBONDMATRIX) THEN
                        IF (HBONDMATCH) THEN
                           WRITE(MYUNIT,34) JP,RANDOM,POTEL,HBONDMARKOV(HBONDMATCHN),NSUCCESS(JP),NFAIL(JP)
                        ENDIF 
                     ELSE
                        WRITE(MYUNIT,34) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
                     ENDIF 
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
#endif
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV(JP)).GT.ECONV) THEN
!                    NRELAX=J1-JACCPREV
!                    IF (RESTART) WRITE(MYUNIT,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
                     JACCPREV=J1
                  ENDIF
                  IF (QDONE.EQ.1) THEN
                     NSUCCESS(JP)=NSUCCESS(JP)+1
                  ELSE
                     NFAIL(JP)=NFAIL(JP)+1
                  ENDIF
                  EPPREV(JP)=EPREV(JP)
                  EPREV(JP)=POTEL
! HBONDMATRIX 
! csw34> If quench is part of existing group, update the Markov energy 
                  IF ((AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))).AND.HBONDMATRIX.AND.HBONDMATCH) THEN
                     HBONDMARKOV(HBONDMATCHN)=POTEL  
! csw34> Update the group max energy if necessary
                     IF (POTEL.GT.HBONDMAXE(HBONDMATCHN)) THEN
                        WRITE(MYUNIT,'(A,I3,A,G20.10)') ' HBONDMATRIX> New highest energy structure for group ',
     &                                                  HBONDMATCHN,' found. E=',POTEL
                        HBONDMAXE(HBONDMATCHN)=POTEL
                     ENDIF
! csw34> Update the group best energy and coordinates if necessary
                     IF (POTEL.LT.HBONDBEST(HBONDMATCHN)) THEN
                        WRITE(MYUNIT,'(A,I3,A,G20.10)') ' HBONDMATRIX> New lowest energy structure for group ',
     &                                                  HBONDMATCHN,' found. E=',POTEL
                        HBONDBEST(HBONDMATCHN)=POTEL
                        HBONDBESTCOORDS(HBONDMATCHN,:)=COORDS(:,JP)
! csw34> If we're identifying a group only by its ligand interactions,
! we need to update the saved matrix when sidechain interactions which
! result in a lower energy are found
                        IF (HBONDLIGAND) THEN
                           WRITE(MYUNIT,'(A,I3)') ' HBONDLIGAND> Updating matrix for group ',HBONDMATCHN
                           HBONDGROUPS(HBONDMATCHN,:,:)=HBONDMAT
! csw34> Dump an updates group*.mat file for the group
                           WRITE(NHBONDGROUPSCHAR,*) HBONDMATCHN
                           HBONDGROUPNAME='group'//TRIM(ADJUSTL(NHBONDGROUPSCHAR))
                           CALL WRITEHBONDMATRIX(HBONDMAT(:,:),HBONDNRES,HBONDGROUPNAME)
                        ENDIF
                     ENDIF
                  ENDIF
! END OF HBONDMATRIX
                  COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                  VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               ELSE
                  NFAIL(JP)=NFAIL(JP)+1
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
!           WRITE(MYUNIT,'(A,4F20.10)') 'Q4 values ',QS,QF,QSTART,QFINISH
!
!  Dump coordinates and energy if another run is attempting to jump to this one.
!
            IF ((NPAR.GT.1).AND.(.NOT.NEWJUMP)) CALL DUMPJ(JP,JUMPTO,NPAR,COORDS(1:3*NATOMS,1:NPAR),NATOMS,EPREV)
!
!  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
!
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
!
!  Check the acceptance ratio.
! 
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) THEN
               IF (RATIOT) THEN
                  CALL FIXRATIO(NQ(JP),NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               ELSE
                  CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
               ENDIF
            ENDIF

            TEMP(JP)=TEMP(JP)*SCALEFAC
#ifdef MPI
#else
            IF (HIT) GOTO 37
            IF (EPREV(JP)<EBEST(JP)) EBEST(JP)=EPREV(JP)
            IF ((REPMATCHT).AND.(ABS(MAXVAL(EBEST(:))-MINVAL(EBEST(:)))<ECONV)) GOTO 37
#endif
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(NQ(JP),EBEST,BESTCOORDS,JBEST,JP)
               ENDIF
            ENDIF
!
! Accept/reject step over grand canonical relxation blocks.
!
            IF (GCBHT.AND.(QGCBH.EQ.GCRELAX)) THEN
!
! did the GC potential decrease overall in the GC relax block?
! this is accept/reject for potential compared with previous best value in previous relaxation block
! (not the value in the Markov chain, the lowest value in the previous block)
!
               IF (MOD(J1-1,PRTFRQ).EQ.0) 
     &          WRITE(MYUNIT,'(A,2G20.10)') 'mc> MC block best and previous best potential values are ',GCEBEST(JP),GCEBESTP(JP)
               IF (GCEBEST(JP).LT.GCEBESTP(JP)) THEN
                  RANDOM=0.0D0
                  ATEST=.TRUE.
               ELSE
                  RANDOM=DPRAND()
                  DUMMY=GCEBEST(JP)-GCEBESTP(JP)
!
! FEBH already includes N! factor and beta^kappa. These cancel for FEBH, but not for size changes
! in GCBH. This is done in quench subroutine - don't duplicate here!
!
!                 IF (FEBHT) THEN
!                    IF (GCNATOMSBEST(JP).EQ.GCNATOMSBESTP(JP)+1) DUMMY=DUMMY-FETEMP*LOG(1.0D0*GCNATOMSBEST(JP))  ! we added an atom
!                    IF (GCNATOMSBEST(JP).EQ.GCNATOMSBESTP(JP)-1) DUMMY=DUMMY+FETEMP*LOG(1.0D0*GCNATOMSBESTP(JP)) ! we subtracted an atom
!                 ENDIF
                  DUMMY=EXP(-DUMMY/MAX(MCTEMP,1.0D-100))
                  IF (DUMMY.GT.RANDOM) THEN
                     ATEST=.TRUE.
                  ELSE
                     ATEST=.FALSE.
                  ENDIF
               ENDIF
               IF (ATEST) THEN
                  IF (DEBUG) WRITE(MYUNIT,2334) JP,RANDOM,GCEBEST(JP),GCEBESTP(JP)
2334              FORMAT('JP,RANDOM,GCEBEST(JP),GCEBESTP(JP)=',I6,F15.5,2G17.7,' ACC')
                  GCEBESTP(JP)=GCEBEST(JP)
                  GCJBESTP(JP)=GCJBEST(JP)
                  GCBESTCOORDSP(1:3*GCNATOMSBEST(JP),JP)=GCBESTCOORDS(1:3*GCNATOMSBEST(JP),JP)
                  GCBESTVATP(1:GCNATOMSBEST(JP),JP)=GCBESTVAT(1:GCNATOMSBEST(JP),JP)
                  GCNATOMSBESTP(JP)=GCNATOMSBEST(JP)
               ELSE ! reset to previous size and reset the Markov chain as well
                  IF (DEBUG) WRITE(MYUNIT,1335) JP,RANDOM,GCEBEST(JP),GCEBESTP(JP)
1335              FORMAT('JP,RANDOM,GCEBEST(JP),GCEBESTP(JP)=',I6,F15.5,2G17.7,' REJ')
                  GCEBEST(JP)=GCEBESTP(JP)
                  GCJBEST(JP)=GCJBESTP(JP)
                  GCBESTCOORDS(1:3*GCNATOMSBESTP(JP),JP)=GCBESTCOORDSP(1:3*GCNATOMSBESTP(JP),JP)
                  GCBESTVAT(1:GCNATOMSBESTP(JP),JP)=GCBESTVATP(1:GCNATOMSBESTP(JP),JP)
                  GCNATOMSBEST(JP)=GCNATOMSBESTP(JP)
                  IF (NATOMS.GT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR-LOG(1.0D0*NATOMS)
                  IF (NATOMS.LT.GCNATOMSBEST(JP)) SYMFCTR=SYMFCTR+LOG(1.0D0*GCNATOMSBEST(JP))
                  NATOMS=GCNATOMSBEST(JP)
                  COORDS(1:3*GCNATOMSBEST(JP),JP)=GCBESTCOORDS(1:3*GCNATOMSBEST(JP),JP)
                  COORDSO(1:3*GCNATOMSBEST(JP),JP)=GCBESTCOORDS(1:3*GCNATOMSBEST(JP),JP)
                  VAT(1:GCNATOMSBEST(JP),JP)=GCBESTVAT(1:GCNATOMSBEST(JP),JP)
                  VATO(1:GCNATOMSBEST(JP),JP)=GCBESTVAT(1:GCNATOMSBEST(JP),JP)
                  POTEL=GCEBEST(JP)
                  EPREV(JP)=GCEBEST(JP)
                  IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'mc> potel and eprev reset to ',GCEBEST(JP)
                  IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'mc> number of atoms reset to ',NATOMS
               ENDIF
            ENDIF
#ifdef MPI

            IF (TARGET) CALL BHPT_TARGET(HITANY)
            IF (HITANY) GOTO 37
            CALL BHPT_EXCHANGE(BETA,J1,NTOT,NEACCEPT)
!
!  When SYMMETRISE is turned on we can have multiple quenches on each pass through the main
!  basin-hopping loop. If we really want to limit the number of quenches to NSTEPS then we
!  need to test NQ(JP) as well.
!
  
            IF (NQ(1).GT.NSTEPS) EXIT
#else
            IF ((NQ(JP).GT.NSTEPS).AND.(.NOT.REPMATCHT)) GOTO 37

         ENDDO
#endif
!
!  ****************************** End of loop over NPAR parallel runs *****************************
!
!  Supersteps
!
!        IF (SUPERSTEP) CALL SUPERMC(SPOTEL,SCOORDS,NSUPERCOUNT,POTEL)
!
!  Renormalised type steps
!
!        IF (RENORM) CALL REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
!        IF (STAY) CALL MYRESET(1,NATOMS,NPAR,NSEED)
  
         CALL FLUSH(MYUNIT)
      ENDDO
#ifdef MPI
      IF (MPIT) THEN
! sf344> Apparently MPI_FINALIZE can only be called once during a run, so we call it at the end of main.F
!         CALL MPI_FINALIZE(MPIERR)
!         IF (DEBUG) WRITE(MYUNIT, '(A,I6)') 'In mc after MPI_FINALIZE MPIERR=',MPIERR
      ENDIF
#endif

37    CONTINUE

#ifdef MPI
!jdf43>Wait for all runs to get to 37 then collate statistics
      CALL BHPT_IO(TIME-TSTART,NQTOT-1,NPCALL-1,NTOT,NEACCEPT)

      WRITE(MYUNIT,10) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),STEP(JP),ASTEP(JP),TEMP(JP)
10    FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
      IF (RIGID) WRITE(MYUNIT,'(A,F12.5)') 'Rigid body orientational step=',OSTEP(JP)

#else
      DO JP=1,NPAR
         IF (NPAR.GT.1) THEN
            IF (RATIOT) THEN 
               IF (MOD(NQ(JP)-1,NACCEPT).NE.0) CALL FIXRATIO(NQ(JP),NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               CALL FINALRATIO(NULLMOVEST,NSUCCESST,JP,SUM(NQ(:)))
            ENDIF 
            WRITE(MYUNIT,20) JP,NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
20          FORMAT('[',I2,']Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ELSE
            IF (RATIOT) THEN
               IF (MOD(J1-1,NACCEPT).NE.0) CALL FIXRATIO(J1,NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               CALL FINALRATIO(NULLMOVEST,NSUCCESST,JP,SUM(NQ(:)))
            ENDIF
            WRITE(MYUNIT,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
21          FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ENDIF
         IF (RIGID) WRITE(MYUNIT,'(A,F12.5)') 'Rigid body orientational step=',OSTEP(JP)
      ENDDO
#endif
!mo361>Deallocating these arrays to cope with multiple runs of this subroutine in GA
      DEALLOCATE(TMOVE)
      DEALLOCATE(OMOVE)
      IF ( DUMPUNIQUE ) CLOSE(DUMPUNIQUEUNIT)
      
! beh35 > mlowest addition

      IF (TARGET .OR. MLOWEST) THEN

         IF (HIT) THEN
            WRITE(MYUNIT,'(I0,1X,I0,1X,I0,1X,F20.10)')1,NQTOT,NPCALL,TIME-TSTART
         ELSE
            WRITE(MYUNIT,'(I0,1X,I0,1X,I0,1X,F20.10)')0,NQTOT,NPCALL,TIME-TSTART
         ENDIF

      ENDIF

      RETURN
      END
C     
C
! csw34> Subroutine to dump the current lowest interaction energy structure 
!> \brief Subroutine to dump the current lowest interaction energy structure 
!> \author Chris Whittleston, csw34@cam.ac.uk
!
      SUBROUTINE A9INTEDUMPLOWEST()
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE
      OPEN(UNIT=20,FILE="bestintE.rst",STATUS='UNKNOWN')  
      WRITE(20,'(g20.10,i5)') INTEQMIN(1), INTEFF(1)
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTEQMINP(1,:) 
      CLOSE(20)
!    csw34> Dump to PDB using routine in amberinterface.f
      CALL A9DUMPPDB(INTEQMINP(1,:),"bestintE")
      END SUBROUTINE A9INTEDUMPLOWEST
! csw34> Subroutine to work out the interaction energy INTE between residue RESLIG for geometry INTECOORDS
      SUBROUTINE A9INTE(INTECOORDS,INTE)
      USE PORFUNCS
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE
      DOUBLE PRECISION :: INTECOORDS(3*NATOMS), INTE
      INTEGER ISTAT
! Dump current coordinates to file for SANDER to read in
      OPEN(UNIT=20,FILE='coords.intres',STATUS='UNKNOWN')
      WRITE(20,'(a20)') 'intE coordinates'
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTECOORDS(:)
      CLOSE(20)
! Run the script that does the interaction energy calculation (downloadable from the group website)
      CALL SYSTEM_SUBR('bash AMBGMINintE.sh',ISTAT)
! Read interaction energy in from the file called intE
      OPEN(UNIT=20,FILE="intE",STATUS='OLD')
      READ(20,*) INTE
      CLOSE(20)
      END SUBROUTINE A9INTE

!
!> \brief Subroutine to check the  distance 
!> \brief between groups of atoms defined in the movableatoms file.
!> Checks the A->B and A->C distances,
!> and if they are greater than the A(B/C)THRESH values defined
!> in the data file, the routine returns DISTOK=.FALSE.
!> \author Chris Whittleston, csw34@cam.ac.uk
!
      SUBROUTINE A9DISTCHECK(COORDS,DISTOK)
      USE MODAMBER9, ONLY : NATOMSINA,NATOMSINB,NATOMSINC,ATOMSINALIST,ATOMSINBLIST,ATOMSINCLIST
      USE COMMONS, ONLY : NATOMS,ABTHRESH,ACTHRESH,DEBUG
      IMPLICIT NONE
      LOGICAL :: DISTOK
      INTEGER :: I,J
      DOUBLE PRECISION :: COORDS(3*NATOMS)
      DOUBLE PRECISION :: CENTREOFA(3),CENTREOFB(3),CENTREOFC(3)  
      DOUBLE PRECISION :: ABDIST,ACDIST 
! initialise variables
      CENTREOFA(:)=0.0D0
      CENTREOFB(:)=0.0D0
      CENTREOFC(:)=0.0D0
      ABDIST=0.0D0
      ACDIST=0.0D0
! find centre of ligand (group A)
      DO I=1,NATOMSINA
         J=ATOMSINALIST(I)
         CENTREOFA(1)=CENTREOFA(1)+COORDS(3*J-2)
         CENTREOFA(2)=CENTREOFA(2)+COORDS(3*J-1)
         CENTREOFA(3)=CENTREOFA(3)+COORDS(3*J  )
      END DO
      CENTREOFA(:) = CENTREOFA(:)/NATOMSINA
! find centre of group B 
      DO I=1,NATOMSINB
         J=ATOMSINBLIST(I)
         CENTREOFB(1)=CENTREOFB(1)+COORDS(3*J-2)
         CENTREOFB(2)=CENTREOFB(2)+COORDS(3*J-1)
         CENTREOFB(3)=CENTREOFB(3)+COORDS(3*J  )
      END DO
      CENTREOFB(:) = CENTREOFB(:)/NATOMSINB
! find centre of group B 
      DO I=1,NATOMSINC
         J=ATOMSINCLIST(I)
         CENTREOFC(1)=CENTREOFC(1)+COORDS(3*J-2)
         CENTREOFC(2)=CENTREOFC(2)+COORDS(3*J-1)
         CENTREOFC(3)=CENTREOFC(3)+COORDS(3*J  )
      END DO
      CENTREOFC(:) = CENTREOFC(:)/NATOMSINC
! calculate A->B distance
      ABDIST=SQRT((CENTREOFA(1)-CENTREOFB(1))**2+(CENTREOFA(2)-CENTREOFB(2))**2+(CENTREOFA(3)-CENTREOFB(3))**2)
! calculate A->C distance
      ACDIST=SQRT((CENTREOFA(1)-CENTREOFC(1))**2+(CENTREOFA(2)-CENTREOFC(2))**2+(CENTREOFA(3)-CENTREOFC(3))**2)
! some DEBUG printing
      IF (DEBUG) THEN
         WRITE(*,*) 'AB distance=',ABDIST
         WRITE(*,*) 'AC distance=',ACDIST
         IF (ABDIST.LT.ABTHRESH) THEN
            WRITE(*,*) 'A->B condition met! :)'
         ELSE 
            WRITE(*,*) 'A->B condition broken :('
         ENDIF   
         IF (ACDIST.LT.ACTHRESH) THEN
            WRITE(*,*) 'A->C condition met! :)'
         ELSE 
            WRITE(*,*) 'A->C condition broken :('
         ENDIF   
      ENDIF
! do the check for both conditions
      IF ((ABDIST.LT.ABTHRESH).AND.(ACDIST.LT.ACTHRESH)) DISTOK=.TRUE.
! more debug printing
      IF (DEBUG) WRITE(*,*) 'DISTOK=',DISTOK
 
      END SUBROUTINE A9DISTCHECK


      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE
      DOUBLE PRECISION ENEW, EOLD, DUMMY, DPRAND, RANDOM, EREF, TEOLD, TENEW, RATIO,MCTEMP
      DOUBLE PRECISION TRANS, DISTMIN, DISTMINOLD
      LOGICAL ATEST, EVAP, EVAPREJECT
      INTEGER NP, J1
      DATA DISTMINOLD /0.0D0/
      COMMON /DMIN/ DISTMIN
!     COMMON /IG/ IGNOREBIN, FIXBIN
      common /ev/ evap, evapreject

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
      IF (TUNNELT) THEN
         TEOLD=TRANS(EOLD,QMIN(1),GAMMA)
         TENEW=TRANS(ENEW,QMIN(1),GAMMA)
!        WRITE(MYUNIT,'(A,4F20.10)') 'TEOLD,TENEW,QMIN(1),GAMMA=',TEOLD,TENEW,QMIN(1),GAMMA
      ELSE
         TEOLD=EOLD
         TENEW=ENEW
      ENDIF

      IF (TSALLIST) THEN
         EREF=QMIN(NP)*1.1D0
         DUMMY=(1.0D0-(1.0D0-QTSALLIS)*(TENEW-EREF)/MCTEMP)/(1.0D0-(1.0D0-QTSALLIS)*(TEOLD-EREF)/MCTEMP)
         DUMMY=DUMMY**(QTSALLIS/(1.0D0-QTSALLIS))
!        WRITE(MYUNIT,'(A,4F20.10)') 'TENEW,TEOLD,EREF,DUMMY=',TENEW,TEOLD,EREF,DUMMY
         IF (DUMMY.GE.1.0D0) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DUMMY.GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF
      ELSEIF(NEWTSALLIST) THEN
         ATEST=.FALSE.
         EREF=QMIN(NP) 
         IF(TENEW<TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSEIF((QTSALLIS-(1.0D0-QTSALLIS)*(TENEW-EREF)/MCTEMP)>0D0) THEN
            DUMMY=(QTSALLIS-(1.0D0-QTSALLIS)*(TENEW-EREF)/MCTEMP)/(QTSALLIS-(1.0D0-QTSALLIS)*(TEOLD-EREF)/MCTEMP)
            DUMMY=DUMMY**(QTSALLIS/(1.0D0-QTSALLIS))
            RANDOM=DPRAND()
            IF (DUMMY.GT.RANDOM) THEN
               ATEST=.TRUE.
            ENDIF
         ENDIF
      ELSE
!
!  Standard canonical sampling.
!
         IF (TENEW.LT.TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(TENEW-TEOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF
      ENDIF 

      RETURN 
      END 

      SUBROUTINE JUMPM(RANNJ,J1,JP,EPPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, J2, UNT, NDUM, ITERATIONS, BRUN, QDONE
      DOUBLE PRECISION RANNJ, RANDOM, DPRAND, EPPREV(NPAR), DUMMY, TIME, EJ, SCREENC(3*NATOMS)

      IF (NEWJUMP) THEN
         IF (RANNJ.LT.PNEWJUMP) THEN
            RANDOM=DPRAND()
            IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))/TEMP(JP)).GT.RANDOM) THEN
!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
               WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' accepted before quench ',J1
               DUMMY=EPREV(JP)
               EPREV(JP)=EPREV(JUMPTO(JP))
               EPREV(JUMPTO(JP))=DUMMY
               DUMMY=EPPREV(JP)
               EPPREV(JP)=EPPREV(JUMPTO(JP))
               EPPREV(JUMPTO(JP))=DUMMY
               DO J2=1,NATOMS
                  DUMMY=VATO(J2,JP)
                  VATO(J2,JP)=VATO(J2,JUMPTO(JP))
                  VATO(J2,JUMPTO(JP))=DUMMY
                  DUMMY=VAT(J2,JP)
                  VAT(J2,JP)=VAT(J2,JUMPTO(JP))
                  VAT(J2,JUMPTO(JP))=DUMMY
               ENDDO
               DO J2=1,3*NATOMS
                  DUMMY=COORDS(J2,JP)
                  COORDS(J2,JP)=COORDS(J2,JUMPTO(JP))
                  COORDS(J2,JUMPTO(JP))=DUMMY
                  DUMMY=COORDSO(J2,JP)
                  COORDSO(J2,JP)=COORDSO(J2,JUMPTO(JP))
                  COORDSO(J2,JUMPTO(JP))=DUMMY
               ENDDO
            ELSE
               WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                  ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' rejected before quench ',J1
            ENDIF
         ENDIF
      ELSE
         UNT=70+JUMPTO(JP)
         REWIND(UNT)           
         NDUM=INT(DPRAND()*(NQ(JUMPTO(JP))-1))
         IF (DEBUG) WRITE(MYUNIT,'(A, G20.10)') 'Should be choosing buffer energy number ',NDUM
         DO J2=1,NDUM
            READ(UNT,*) EJ
         ENDDO
         DO J2=1,NQ(JUMPTO(JP))-1-NDUM
            READ(UNT,*) DUMMY
         ENDDO
         RANDOM=DPRAND()
!
!  Coordinates are only read if the jump is successful.
!
         IF (DEXP((EPREV(JP)-EJ)*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
            WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' accepted before quench ',NQ(JP)
            EPREV(JP)=EJ
            UNT=70+NPAR+JUMPTO(JP)
            REWIND(UNT)           
            DO J2=1,(NDUM-1)*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
            READ(UNT,*) (COORDS(J2,JP),J2=1,3*NATOMS)
!
!  Coordinates should be converged already, but we need to reset VAT and VATO.
!
!  next line should be uncommented if routine is made availabe to use with CHARMM
!            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)

            DO J2=1,NATOMS
               VATO(J2,JP)=VAT(J2,JP)
            ENDDO
            DO J2=1,3*NATOMS
               COORDSO(J2,JP)=COORDS(J2,JP)
            ENDDO
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(A)') 'Jump coordinates:'
               WRITE(MYUNIT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
            ENDIF
            DO J2=1,NATOMS*(NQ(JUMPTO(JP))-1)-NDUM*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' rejected before quench ',NQ(JP)
         ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE DUMPJ(JP,JUMPTO,NPAR,COORDS,NATOMS,EPREV)
      IMPLICIT NONE
      LOGICAL TEST
      INTEGER NPAR, J2, JP, JUMPTO(NPAR), NATOMS, UNT
      DOUBLE PRECISION COORDS(3*NATOMS,NPAR), EPREV(NPAR)

      DO J2=1,NPAR
         IF (JUMPTO(J2).EQ.JP) TEST=.TRUE.
      ENDDO
!     TEST=.FALSE.
!     IF (TEST) THEN
!        UNT=70+JP
!        WRITE(UNT,'(F20.10)') EPREV(JP)
!        UNT=70+NPAR+JP
!        WRITE(UNT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
!     ENDIF

      RETURN
      END

      SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
      USE COMMONS,ONLY : COORDS,COORDSO,VAT,VATO
      IMPLICIT NONE
      INTEGER JP, NSEED, J2, NATOMS, NPAR

      DO J2=1,3*(NATOMS-NSEED)
         COORDS(J2,JP)=COORDSO(J2,JP)
      ENDDO
      DO J2=1,NATOMS
         VAT(J2,JP)=VATO(J2,JP)
      ENDDO

      RETURN
      END
 
 
 
      SUBROUTINE REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
      USE commons
      IMPLICIT NONE
      INTEGER ITERATIONS, J2, JACCPREV, J1, NQTOT, BRUN, QDONE
      DOUBLE PRECISION TIME, POTEL, RCOORDS(3*NATOMS), RMIN, RVAT(NATOMS), SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

10    CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,NHSRESTART)
!  next line should be uncommented if routine is made availabe to use with CHARMM
!      IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
      CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!
!  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
!
      IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
         DO J2=1,3*NATOMS
            COORDS(J2,1)=COORDSO(J2,1)
         ENDDO
         GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(MYUNIT,'(A,I7,A,G20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ(1),' E=',
     1              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
      DO J2=1,3*NATOMS
         COORDSO(J2,1)=COORDS(J2,1)
         RCOORDS(J2)=COORDS(J2,1)
      ENDDO
      DO J2=1,NATOMS
         VATO(J2,1)=VAT(J2,1)
         RVAT(J2)=VAT(J2,1)
      ENDDO
      EPREV(1)=POTEL
      RMIN=POTEL

      RETURN
      END
 
 
 
      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      USE commons
      USE modcharmm
      USE grouprotmod
      IMPLICIT NONE
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), JP, NFAILT(NPAR), NSUCCESST(NPAR)
      LOGICAL evap, evapreject
      DOUBLE PRECISION P0,FAC
!     COMMON /IG/ IGNOREBIN, FIXBIN
!     COMMON /MOVE/ TMOVE, OMOVE
      common /ev/ evap, evapreject
! Calculate the fraction of steps that have succeeded, the acceptance ratio
      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
! Set the scaling factor to be used to scale STEP/TEMPERATURE
! Calculate the scaling factor using the 'acceptance-ratio method' (ARM)
      IF(ARMT) THEN
        FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
      ELSE
! Otherwise just set it to 1.05 
        FAC=1.05D0
      ENDIF
! CASE 1: P0 > ACCRAT
! If the current acceptance ratio (P0) is larger than the target acceptance ratio (ACCRAT)
! we want to increase the step size or reduce the temperature to try to reduce it
      IF (P0.GT.ACCRAT(JP)) THEN
! If both the STEP and TEMPERATURE are fixed with FIXBOTH, do nothing
         IF (FIXBOTH(JP)) THEN
! Else, if only the STEP is fixed, scale down the TEMPERATURE
         ELSE IF (FIXSTEP(JP)) THEN
! Sanity check
            IF (.NOT.FIXTEMP(JP)) THEN
               TEMP(JP)=TEMP(JP)/FAC
            ENDIF
         ELSE
! FIXD is NOT DOCUMENTED :/ It appears to be to do with hard-sphere moves
            IF (FIXD) THEN
               NHSMOVE=NHSMOVE+1 
            ELSE
! Otherwise scale up the STEP size
! RIGID seems to be one of Szilard's lovely undocumented keywords. Great.
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)*FAC
! OSTEP is the maximum displacement of an angle-axis coordinate for a rigid body system  
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)*FAC
               ELSE
                  STEP(JP)=STEP(JP)*FAC
! GROUPROTATION scaling
                  IF(GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) THEN 
                     CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,FAC)
                  ENDIF
! Scale up the CHARMM RIGID translation and rotation step sizes 
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=TRANSMAX*FAC
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=ROTMAX*FAC  
               ENDIF
            ENDIF
! ASTEP is the maximum angular step size? Not 100% sure, may be a threshold for Morse/LJ...
            ASTEP(JP)=ASTEP(JP)*FAC
! jwrm2> limit step size for percolation to the cutoff distance for determining connectivity
            IF (PERCOLATET .AND. ( STEP(JP) .GT. PERCCUT))  STEP(JP) = PERCCUT
            IF (PERCOLATET .AND. (ASTEP(JP) .GT. PERCCUT)) ASTEP(JP) = PERCCUT
            IF (PERCOLATET .AND. (OSTEP(JP) .GT. PERCCUT)) OSTEP(JP) = PERCCUT
         ENDIF
! CASE 2: P0 < ACCRAT
! If the current acceptance ratio (P0) is smaller than the target acceptance ratio (ACCRAT)
! we want to decrease the step size or increase the temperature to try to increase it
      ELSE
! If both the STEP and TEMPERATURE are fixed with FIXBOTH, do nothing
         IF (FIXBOTH(JP)) THEN
! Else, if only the TEMPERATURE is fixed, scale up the TEMPERATURE
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) THEN 
               TEMP(JP)=TEMP(JP)*FAC
            ENDIF
         ELSE
! FIXD is NOT DOCUMENTED :/ It appears to be to do with hard-sphere moves
            IF (FIXD) THEN
               NHSMOVE=MAX(1,NHSMOVE-1)
            ELSE
! Otherwise scale down the STEP size
! RIGID seems to be one of Szilard's lovely undocumented keywords. Great.
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)/FAC
! OSTEP is the maximum displacement of an angle-axis coordinate for a rigid body system  
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/FAC
               ELSE
! This is the usual case for most GMIN runs
                  STEP(JP)=STEP(JP)/FAC
! GROUPROTATION scaling
                  IF(GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) THEN 
                     CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,1.0D0/FAC)
                  ENDIF
! Scale up the CHARMM RIGID translation and rotation step sizes 
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=TRANSMAX/FAC
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=ROTMAX/FAC
               ENDIF
            ENDIF
! ASTEP is the maximum angular step size? Not 100% sure, may be a threshold for Morse/LJ...
            ASTEP(JP)=ASTEP(JP)/FAC
         ENDIF
      ENDIF
!
! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!
!
      STEP(JP)=MIN(STEP(JP),1.0D3)
      OSTEP(JP)=MIN(OSTEP(JP),1.0D3)
      ASTEP(JP)=MIN(ASTEP(JP),1.0D3)
! Do some printing to give the user feedback
      IF (MOD(NQ(JP),PRTFRQ).EQ.0) THEN
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I0.2,A,I6,A,F8.4,A,F8.4)') '[',JP,']Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
         ELSE
            WRITE(MYUNIT,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I0.2,A)',ADVANCE='NO') '[',JP,']Steps are now:'
            ELSE
               WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Steps are now:'
            ENDIF
            WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP(JP)    
            IF(ASTEP(JP).GT.0.D0) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP(JP) 
            IF(CHRIGIDTRANST.AND.CHRMMT) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  TRANSMAX=',TRANSMAX
            IF(CHRIGIDROTT.AND.CHRMMT) WRITE(MYUNIT,'(A,F10.4)')'  ROTMAX=',ROTMAX
            IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F10.4)') ' Temperature is now:',TEMP(JP)
            IF (RIGID) WRITE(MYUNIT,'(A,F12.6,A,F12.6)') 'Maximum rigid body rotational move is now ',OSTEP(JP)

         ENDIF
         IF (FIXD) WRITE(MYUNIT,'(A,I4)') 'hard sphere collision moves=',NHSMOVE
      ENDIF
!
      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 
!
      RETURN
      END


      SUBROUTINE REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
      USE commons

      IMPLICIT NONE
      LOGICAL STAY, REJECT, METROPOLIS
      INTEGER J1, J2, NREN, ITERATIONS, NQTOT, NLAST, JACCPREV, NSTEPREN, J3, BRUN, QDONE
      DOUBLE PRECISION POTEL, RMIN, RCOORDS(3*NATOMS), RVAT(NATOMS), RANDOM, DPRAND, RMINO, RCOORDSO(3*NATOMS), RVATO(NATOMS),
     1                 TIME, XIP, DUMMY, SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

      STAY=.FALSE.
      IF (POTEL.LT.RMIN) THEN
         RMIN=POTEL          
         DO J2=1,3*NATOMS
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            RVAT(J2)=VAT(J2,1)
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'RMIN,POTEL=',RMIN,POTEL
!     PRINT*,'J1,JACCPREV+NRENSTUCK,NLAST+NREN=',J1,JACCPREV+NRENSTUCK,NLAST+NREN
      IF ((J1.GE.JACCPREV+NRENSTUCK).OR.(J1.GE.NLAST+NREN)) THEN
         JACCPREV=J1
         NLAST=J1
         RANDOM=DPRAND()
         METROPOLIS=DEXP(-(RMIN-RMINO)/TRENORM).GT.RANDOM
         REJECT=.FALSE.
!
!  Taboo list for renormalised energies. Skip if the step is going to be rejected by
!  the Metropolis condition.
C
         IF (TABOOT.AND.METROPOLIS) THEN
            IF (NSTEPREN.EQ.0) NT(1)=0
            CALL NEWINERTIA(RCOORDS,NATOMS,NATOMS,XIP)
            DO J1=1,NT(1)
               IF (DABS(RMIN-ESAVE(J1,1)).LT.ECONV) THEN
                  IF (DABS(XIP-XINSAVE(J1,1))/(XIP+XINSAVE(J1,1)).LT.1.0D-2) THEN
                     REJECT=.TRUE.
                     GOTO 20
                  ELSE
                     WRITE(MYUNIT,'(A, 2G20.10)') 'Energies nearly degenerate:',RMIN,ESAVE(J1,1)
                     WRITE(MYUNIT,'(A, 2G20.10)') 'But  different  structures:',XIP,XINSAVE(J1,1)
                  ENDIF
               ENDIF
               IF (RMIN.LT.ESAVE(J1,1)) THEN
                  NT(1)=MIN(NT(1)+1,NTAB)
                  DO J3=NT(1),J1+1,-1
                     ESAVE(J3,1)=ESAVE(J3-1,1)
                     XINSAVE(J3,1)=XINSAVE(J3-1,1)
                  ENDDO
                  ESAVE(J1,1)=RMIN
                  XINSAVE(J1,1)=XIP
                  GOTO 20
               ENDIF
            ENDDO
            
            NT(1)=MIN(NT(1)+1,NTAB)
            ESAVE(NT(1),1)=RMIN
            XINSAVE(NT(1),1)=XIP

20          CONTINUE

            WRITE(MYUNIT,'(A,I10)') ' Number of entries in taboo list=',NT(1)
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
!
!  Accept/reject for renormalised step
!
         IF (METROPOLIS.AND.(.NOT.REJECT)) THEN
            IF (NSTEPREN.GT.0) WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' accepted'
            NREN=MAX(NREN/1.1D0,NRENORM/2.0D0)
            RMINO=RMIN
            IF (.NOT.STAY) THEN
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=RCOORDS(J2)
                  COORDS(J2,1)=RCOORDS(J2)
               ENDDO
               DO J2=1,NATOMS
                  RVATO(J2)=RVAT(J2)
               ENDDO
            ELSE
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=COORDSO(J2,1)
                  COORDS(J2,1)=COORDSO(J2,1)
               ENDDO
            ENDIF
         ELSE
            IF (REJECT) THEN
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected by taboo criterion'
            ELSE
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected'
            ENDIF
            DO J2=1,3*NATOMS
               COORDS(J2,1)=RCOORDSO(J2)
            ENDDO
            DO J2=1,NATOMS
               RVAT(J2)=RVATO(J2)
            ENDDO
            NREN=NREN*1.1D0
         ENDIF
         NSTEPREN=NSTEPREN+1
         IF (NSTEPREN.EQ.1) WRITE(MYUNIT,'(A,G20.10)') ' first renorm energy is ',RMIN
         WRITE(MYUNIT,'(A,I6)') ' renormalisation interval is now ',NREN
!
!  Longer renorm step
!
10       IF ((XMOVERENORM.GT.3.0D0).OR.FIXD) THEN
            CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,INT(XMOVERENORM))
         ELSE
            DUMMY=STEP(1)
            STEP(1)=XMOVERENORM
            CALL TAKESTEP(1)
            STEP(1)=DUMMY
         ENDIF
!  next line should be uncommented if routine is made availabe to use with CHARMM
!         IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!
!  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
!
         IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=COORDSO(J2,1)
            ENDDO
            GOTO 10
         ENDIF
         NQTOT=NQTOT+1
         WRITE(MYUNIT,'(A,I7,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Renorm Qu ',NQ(1),' E=',
     1        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
         DO J2=1,3*NATOMS
            COORDSO(J2,1)=COORDS(J2,1)
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            VATO(J2,1)=VAT(J2,1)
            RVAT(J2)=VAT(J2,1)
         ENDDO
         EPREV(1)=POTEL
         RMIN=POTEL
      ENDIF

      RETURN
      END
!
!  Reseed if the energy has not improved by more than ECONV over the
!  last NRELAX mc steps.
!  If AVOID is true then save the energy and coordinates of the lowest
!  minimum achieved before each restart and restart if we get too close
!  to any one of them using bipartite matching and mind. Note that bipartite
!  matching and mind can give a local minimum of distance if the optimal 
!  permutation-inversion isn;t found. Using ORIENT to effect a standard
!  orientation first seems to help. It should at least ensure that permutation-inversion
!  isomers are always found. 
!  Would it be possible to just use the energy and inertia components to identify
!  permutation-inversion isomers instead of MINPERM and MIND?
!  This would be like using a small value of AVOIDDIST, which doesn;t seem to be
!  as good.
!
      SUBROUTINE NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST,RES1,RES2)
      USE commons
      use genrigid
      USE modamber9, only : mdstept 
      USE ROTATIONS
      USE VEC3
 
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J4, NEWRESATTEMPTS
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2, DUMMY2,
     1                 RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS), DPRAND, FCOORDS(3*NATOMS), 
     2                 DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS), EBESTP, BESTCOORDSP(3*NATOMS), RMAT(3,3), QENERGY,
     &                 ROTA(3,3), ROTINVA(3,3)
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), NORBIT1, NORBIT2, NPOSITION,  NSUCCESST(NPAR), NFAILT(NPAR)
      LOGICAL RES1, RES2, HIGHEST
      COMMON /MYPOT/ QENERGY
! hk286
      DOUBLE PRECISION :: GRAD(3*NATOMS)
! hk286

      SR3=DSQRT(3.0D0)

! csw34> initialise the attempts counter
      NEWRESATTEMPTS = 1
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
!
!  Reseed if the energy has not improved in the last NRELAX mc cycles,
!  or if the current minimum is too close to one of the NMSBSAVE structures
!  saved in MSBCOORDS.
!
!  Instead of using the current minimum, employ the current minimum in 
!  BESTCOORDS for the AVOID check. Then we need only do the AVOID check when we have
!  a new minimum in BESTCOORDS, i.e. J1.EQ.JBEST(JP).
!
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     WRITE(MYUNIT,'(A,I5,2G17.7,3I5,L8)') 'J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
!     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)

            IF (DEBUG) THEN
               CALL MINPERMDIST(FCOORDS,XMSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
               WRITE(MYUNIT,'(A,F10.3)')'newres> minimum distance is ',DISTANCE
            ENDIF

!
!  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
!  but the coordinates in XMSB do. DISTANCE is the distance in this case.
!
! If the energy is lower no reseeding regardless of separation ? DJW
!           PRINT '(A,2G20.10,L8)','POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV=',
!    1                              POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            XMSB(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J2)
            CALL MINPERMDIST(FCOORDS,XMSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(MYUNIT,'(A,G20.10,A,I6,A,G20.10,A,F10.3)') 'newres> Minimum energy ',POTEL,
     &                                    ' is too close to saved structure ',
     &                                 J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
         ENDDO savedloop
      ENDIF

20    CONTINUE 


      IF (RES1.OR.RES2) THEN
         COORDS(:,JP) = COORDSO(:, JP)
          !  Does not seem necessary to save MSB data for restarts.
          !  If we are reseeding because RES2 is true then:
          !  (1) If the AVOID test is based on the distance for BESTCOORDS, add
          !      these to the  MSB list if the distance is > 0.01D0, i.e. a different structure.
          !      Otherwise, add the coordinates of the previous lowest energy minimum to
          !      the MSB list; these have been saved in BESTCOORDSP and the energy in EBESTP.
          !  (2) Alternatively, if the AVOID test is based on the coordinates of the current
          !      structure, we should add these to the MSB list instead. Testing
          !      every step involves an overhead, which grows with the number of saved taboo
          !      structures. We could make the taboo list cyclic, so that new structures to be
          !      avoided replace old ones at the beginning of the list when the length of the
          !      list is exceeded. We are currently only doing the bipartite matching test if
          !      the energy of the lowest minimum since the last reseeding has changed, which
          !      saves time.
          !
          !  Change to a cyclic AVOID list 30/12/04.

         IF (RES2.AND.(.NOT.AVOIDRESEEDT)) THEN 
            POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
            WRITE(MYUNIT,'(A,I8,A)') 'newres> Resetting energy to a large value due to taboo condition'
            RETURN
         ELSEIF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.1D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)

            CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)

            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
! csw34> store the number of atoms in the new structure for dumping in DUMPSTATE
            AVOIDNATOMS(NPOSITION) = NATOMS

            WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving current best minimum to position ',NPOSITION,
     &                         ' in the AVOID list E=',EBEST(JP)
!           OPEN(UNIT=34,FILE='MSBdata',POSITION='APPEND')
!           WRITE(34,'(G20.10)') MSBE(NPOSITION) 
!           WRITE(34,'(3G20.10)') BESTCOORDS(1:3*NATOMS,JP)
!           CLOSE(34)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
!           IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
               MSBE(NPOSITION)=EBESTP
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)

               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
! csw34> store the number of atoms in the new structure for dumping in DUMPSTATE
               AVOIDNATOMS(NPOSITION) = NATOMS

               WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving previous best minimum to position ',NPOSITION,
     1                             ' in the AVOID list E=',EBESTP
!           ENDIF
         ENDIF ! end new condition

         IF (NEWRESTART_MD) THEN ! lb415
            IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
            WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding via a short high temperature MD run'
            CHANGE_TEMP = .true.

! th368: 20-10-2009 Extending MD-Reseeding to the CHARMM interface
! terminate if neither AMBER nor CHARMM was requested in the data file

             IF (AMBERT) THEN
               CALL TAKESTEPAMBER(JP,COORDS(:,JP),movableatomlist,nmovableatoms,ligmovet,mdstept,randomseedt,
     &                                  BLOCKMOVET,NBLOCKS,ATOMSINBLOCK)
             ELSEIF (CHRMMT) THEN
! hk286 - run MD from previously saved energy coordinates
! csw34> Reseed from Markov chain coordinates to minimise chance of getting stuck in a loop, unable to quench
               COORDS(1:3*NATOMS,JP)=COORDSO(1:3*NATOMS,JP)
               CALL POTENTIAL(COORDS(:,1),GRAD, POTEL, .TRUE., .FALSE.)
! hk286
               CALL CHMD(JP)
! hk286 - second revision
               CALL TRANSFORMCTORIGID(COORDS(:,1), RIGIDCOORDS)
               CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY, COORDS(:,1), RIGIDCOORDS)
! hk286
             ELSE
               WRITE(MYUNIT,'(A,I8,A)') 'newres> Molecular Dynamics reseeding is available for AMBER or CHARMM runs only.'
               STOP
             ENDIF
! end th368: 20-10-2009

            CHANGE_TEMP = .false.
         ELSE
            IF (NHSRESTART.GT.0) THEN
               IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
               IF (RES2) WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               IF (SHELLMOVES(JP)) WRITE(MYUNIT,'(A)') 'newres> Turning off shell moves'
               SHELLMOVES(JP)=.FALSE.
               CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
            ELSE
               IF (RES1) THEN
                  IF (AVOIDRESEEDT) THEN
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - reseeding'
                  ELSE
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP)
                  ENDIF
               ENDIF
               IF (RES2) WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               HIGHEST=.TRUE.
               DO J2=1,NPAR
                  IF (J2.EQ.JP) CYCLE
                  IF (EBEST(J2).GT.EBEST(JP)) THEN
                     HIGHEST=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF (HIGHEST) THEN
                  WRITE(MYUNIT,'(A,I6,A,F20.10)') 'newres> Parallel run ',JP,' has the highest energy ',EBEST(JP)
                  IF (NPAR.GT.1) WRITE(MYUNIT,'(6F20.10)') EBEST(1:NPAR)
                  IF (.NOT.AVOIDRESEEDT) THEN
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Resetting energy to a large value to reject this step (not reseeding)'
                     JBEST(JP)=J1
                     EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
                     POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
                     RETURN
                  ELSE
                     WRITE(MYUNIT,'(A)') 'newres> Full reseeding'
                     IF(DMACRYST) THEN ! vr274> special reseeding for DMACRYS
                        CALL DMACRYS_GENRANDOM(COORDS(:,JP))
                     ELSE IF (RIGID) THEN    ! jdf43> reseeding with uniform angle-axis distro & uniform spatial distro
                        DO J2=1,NATOMS/2
                           J4=J2*3
                           COORDS(J4-2:J4,JP)=VEC_RANDOM()*(DPRAND()**(1.D0/3.D0))*DSQRT(RADIUS)
                           J4=(J2+NATOMS/2)*3
                           COORDS(J4-2:J4,JP)=ROT_RANDOM_AA()
                        ENDDO
                     ELSE
                        DO J2=1,3*NATOMS
                           RANDOM=(DPRAND()-0.5D0)*2.0D0
                           COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
                        ENDDO
                     ENDIF
                     NCORE(JP)=0
                     PTGROUP(JP)='   '
                  ENDIF
               ELSEIF (NCORE(JP).GT.0) THEN
!         
!  Reseed everything if this is the lowest of a set of parallel runs, otherwise
!  just reseed the surface.
!         
                  WRITE(MYUNIT,'(A)') 'newres> Accepting an enforced surface reseeding'
                  DUMMY2=-1.0D0
                  DO J2=NATOMS-NCORE(JP)+1, NATOMS
                     DISTANCE=COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2
                     IF (DISTANCE.GT.DUMMY2) DUMMY2=DISTANCE
                  ENDDO
                  DISTANCE=SQRT(DISTANCE)
                  WRITE(MYUNIT,'(A,F15.5)') 'newres> largest core radius=',DISTANCE
                  DO J2=1,NATOMS-NCORE(JP)
                     COORDS(3*(J2-1)+1,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+2,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+3,JP)=(DPRAND()-0.5D0)*2.0D0
                     RANDOM=DISTANCE+1.0D0+DPRAND()*(SQRT(RADIUS)-DISTANCE-1.0D0)
                     DUMMY2=SQRT(COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2)
                     COORDS(3*(J2-1)+1,JP)=COORDS(3*(J2-1)+1,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+2,JP)=COORDS(3*(J2-1)+2,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+3,JP)=COORDS(3*(J2-1)+3,JP)*RANDOM/DUMMY2
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
! jdf43> 
         IF (DDMT) THEN
            DO J2=1,3*NATOMS
               CALL DDMCONDENSE(COORDS)
            ENDDO
         ENDIF

!
!  This step will be accepted because EPREV(JP)=POTEL. However, there
!  could be additional steps, so we need
!  to change COORDSO and VATO. Should reset EBEST and BESTCOORDS as well.
!
!  next line should be uncommented if routine is made available to use with CHARMM
!         IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
654      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
! csw34> Include a check for the number of restart attempts to prevent infinate loop behaviour
         IF ((RMS.GT.BQMAX).AND.(NEWRESATTEMPTS.LE.10)) THEN
! csw34> Increment the attempts counter
            NEWRESATTEMPTS=NEWRESATTEMPTS+1
            WRITE(MYUNIT,'(A)') 'newres> Quench from reseeded geometry failed - try again'

! hk286 - first revision
            IF (NEWRESTART_MD) THEN
               CHANGE_TEMP = .true.
               WRITE(MYUNIT,'(A)') 'newres> Running high temperature MD'
               IF (AMBERT) THEN
                  CALL TAKESTEPAMBER(JP,COORDS(:,JP),movableatomlist,nmovableatoms,ligmovet,mdstept,randomseedt,
     &                                  BLOCKMOVET,NBLOCKS,ATOMSINBLOCK)
               ELSEIF (CHRMMT) THEN
                  CALL CHMD(JP)
               ENDIF
               CHANGE_TEMP = .false.

            ELSEIF(DMACRYST) THEN ! vr274> special reseeding for DMACRYS
                  CALL DMACRYS_GENRANDOM(COORDS(:,JP))
            ELSEIF (RIGID) THEN    ! jdf43> reseeding with uniform angle-axis distro & uniform spatial distro
               DO J2=1,NATOMS/2
                  J4=J2*3
                  COORDS(J4-2:J4,JP)=VEC_RANDOM()*(DPRAND()**(1.D0/3.D0))*DSQRT(RADIUS)
                  J4=(J2+NATOMS/2)*3
                  COORDS(J4-2:J4,JP)=ROT_RANDOM_AA()
               ENDDO
            ELSE
! hk286
               DO J2=1,3*NATOMS
                  RANDOM=(DPRAND()-0.5D0)*2.0D0
                  COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
               ENDDO
            ENDIF
! jdf43> 
           IF (DDMT) THEN
              DO J2=1,3*NATOMS
                 CALL DDMCONDENSE(COORDS)
              ENDDO
           ENDIF
! hk286          
            GOTO 654
         ENDIF
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         POTEL=QENERGY
         EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=J1
         EPREV(JP)=POTEL !
         EPPREV(JP)=0.0D0
         NSYMREM=0
! th368: 20.10.2009 updates should only take place in case
! reseeding took place - transfered from subroutine MC lines 363 to 367

         IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
         COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!        WRITE(MYUNIT,'(A,2G20.10)') 'newres> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

! end th368 20.10.2009
      ENDIF

      RETURN
      END

! js850> if the structure just found is unique (i.e. if its energy hasn't been
! recently found), then dump it to file dumpunique.
! note: quench energy is in POTEL.  markov energy is in EPREV(JP)
      SUBROUTINE DUMPUNIQUE_DUMP(DUMPUNIQUEUNIT, JP, POTEL, DUMPUNIQUEEPREV )
        USE COMMONS, ONLY : NATOMS, ECONV, EPREV, NQ, COORDS
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: DUMPUNIQUEUNIT, JP  
        DOUBLE PRECISION, INTENT(IN) :: POTEL
        DOUBLE PRECISION, INTENT(INOUT) :: DUMPUNIQUEEPREV
        INTEGER J2
        
        WRITE(*,*) NQ(JP), "POTEL", POTEL, "DUMPUNIQUEEPREV", DUMPUNIQUEEPREV, "EPREV", EPREV(JP), JP
        IF ( (ABS(POTEL-DUMPUNIQUEEPREV).GT.ECONV).AND.
     &       (ABS(POTEL-EPREV(JP)).GT.ECONV) ) THEN
          WRITE(DUMPUNIQUEUNIT,*) NATOMS!,  POTEL, EPREV(JP), Jp
          WRITE(DUMPUNIQUEUNIT,*) "Qu ", NQ(JP), "Energy ", POTEL
          DO J2=1,NATOMS
            WRITE(DUMPUNIQUEUNIT,'(A,3F20.10)') 'LA  ', COORDS( 3*(J2-1)+1, JP), COORDS( 3*(J2-1)+2, JP), COORDS( 3*(J2-1)+3, JP)
            !WRITE(DUMPUNIQUEUNIT,'(A,3F20.10)') 'LA  ', SCREENC( 3*(J2-1)+1), SCREENC( 3*(J2-1)+2), SCREENC( 3*(J2-1)+3)
          END DO
        ENDIF
        DUMPUNIQUEEPREV=POTEL
      END SUBROUTINE DUMPUNIQUE_DUMP
