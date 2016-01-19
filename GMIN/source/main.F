C   GMIN: A program for finding global minima
C   Copyright (C) 1999- David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      
      PROGRAM GMIN
      USE MODMXATMS
      USE NOA
      USE COMMONS
      USE QMODULE
      USE PERMU
      USE F1COM
      USE MODAMBER
      USE MODAMBER9, only : AMBFINALIO_NODE,MDCRD_UNIT,MDINFO_UNIT,AMBPDB_UNIT, ATMASS1
      USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_MASSES
      USE MODCHARMM
      USE PORFUNCS
      USE TWIST_MOD
      USE HOMOREFMOD
      USE GENRIGID, only : RIGIDINIT, GENRIGID_READ_FROM_FILE, DEGFREEDOMS
      USE GAUSS_MOD, ONLY: KEGEN
      IMPLICIT NONE
      !EXTERNAL READ_CMD_ARGS
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      INTEGER J1,J2, JP, MPIERR, NDUMMY3,NPTOTAL,VERSIONTEMP,GETUNIT,LUNIT
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:)
      DOUBLE PRECISION POTEL, DUMMY, INTFREEZETOLSAVE, RMAT(3,3), DUMMY2, DIST2, LINTCONSTRAINTTOL
      DOUBLE PRECISION, ALLOCATABLE :: TMPCOORDS(:), ENDCOORDS(:,:)
      INTEGER, ALLOCATABLE :: NDUMMY(:), NDUMMY2(:,:)
      LOGICAL LOPEN, YESNO
      INTEGER J6, EPOCH

      CHARACTER(LEN=130) ISTR,JSTR
      CHARACTER(LEN=40) :: atom1,atom2,atompair

      CHARACTER(LEN=13) :: CUDAFILENAME1
      CHARACTER(LEN=19) :: CUDAFILENAME2
      CHARACTER(LEN=21) :: CUDAFILENAME3
      CHARACTER(LEN=24) :: CUDAFILENAME4
      CHARACTER(LEN=17) :: CUDAFILENAME5

      COMMON /MYPOT/ POTEL

      CALL CPU_TIME(TSTART)
      CALL READ_CMD_ARGS
#ifdef MPI
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      MYUNIT=22980+MYNODE
      MDCRD_UNIT=20000+MYNODE
      MDINFO_UNIT=21000+MYNODE
      AMBPDB_UNIT=18000+MYNODE
      AMBFINALIO_NODE=200+MYNODE
      WRITE (ISTR, '(I10)') MYNODE+1
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') INFILE='output'
      IF (NPAR.GT.1) THEN 
         MYFILENAME=TRIM(ADJUSTL(INFILE))//"."//TRIM(ADJUSTL(ISTR))
      ELSE
         MYFILENAME=TRIM(ADJUSTL(INFILE))
      ENDIF
      OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      WRITE(MYUNIT, '(A,I10,A,I10)') "Starting parallel execution: Processor", mynode+1, " of ",NPAR
#else
      NPAR=1
      MYNODE=0
      MYUNIT=22979+1
      MDCRD_UNIT=20000
      MDINFO_UNIT=21000
      AMBPDB_UNIT=18000
      AMBFINALIO_NODE=1
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') THEN
         MYFILENAME="output"
      ELSE
         MYFILENAME=TRIM(ADJUSTL(INFILE))
      ENDIF
      OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Starting serial execution"
#endif

! csw34> TIME LIMITING BINARIES
!        If you want to compile a time limited binary, uncomment the below
!        and replace the time in seconds after IF (EPOCH.LT. with your own
!        and update the error messages to make sense in your case

! find the number of seconds since 0000 on 1/1/1970 (unix epoch)
!      EPOCH=TIME()
! compare this to the number that will have passed at 1am on 5/11/15
!      IF (EPOCH.LT.1449277200) THEN
!         WRITE(MYUNIT, '(A)') "NOTE: Binary only licensed for use at AlgoSB2015"
!      ELSE
! if after the workshop has ended, stop execution with error
!         WRITE(MYUNIT, '(A)') "ERROR: license missing, contact dw34@cam.ac.uk. Stopping."
!         STOP   
!      END IF

! Add the GMIN version to the output - helps bug hunting :)
      !WRITE(MYUNIT, '(A,I5)') 'GMIN version r',VERSIONTEMP
!     CALL DISPLAY_VERSION(MYUNIT)
      CALL COUNTATOMS(MYUNIT, NPAR, GCBHT, GCMU, GCNATOMS)
      IF (CHRMMT) WRITE(MYUNIT,'(A,I8)') 'main> MAXAIM parameter for CHARMM MXATMS=',MXATMS
!
! NATOMS is set to NUMBER_OF_ATOMS in modcommonunit
!
      NATOMSALLOC=NUMBER_OF_ATOMS
      NATOMS=NUMBER_OF_ATOMS
      IF (GCBHT) THEN
         NATOMSALLOC=GCNATOMS
         IF (GCNATOMS.LT.NUMBER_OF_ATOMS) THEN
            WRITE(MYUNIT,'(A)') 'main> *** ERROR - number of atoms exceeds maximum allowed by GCBH keyword'
            STOP
         ENDIF
         ALLOCATE(QENERGIES(GCNATOMS),QCOORDINATES(GCNATOMS,3*GCNATOMS),QPE(GCNATOMS))
         QENERGIES(1:GCNATOMS)=HUGE(1.0D0)/1.0D2
         QPE(1:GCNATOMS)=HUGE(1.0D0)/1.0D2
         QCOORDINATES(1:GCNATOMS,1:3*GCNATOMS)=0.0D0
      ENDIF

! DMACRYS things
! Moved from commons.f90, since they forced commons to depend on countatoms.f90
      allocate( ANV(NATOMSALLOC,NATOMSALLOC,3))

! vr274> DMACRYS is used
      DMACRYST = .false.
      DMACRYS_RANDOMSTART=.false.
      DMACRYS_EXPAND = 1.0D0
      DMACRYS_LATTICE_STEP = 0.0D0
! END DMACRYS things
 
      ALLOCATE(FIN(3*NATOMSALLOC))
      ALLOCATE(XICOM(3*NATOMSALLOC),PCOM(3*NATOMSALLOC))
      ALLOCATE(SCREENC(3*NATOMSALLOC))
      ALLOCATE(IATNUM(NATOMSALLOC), VT(NATOMSALLOC), ZSYM(NATOMSALLOC))
      VT(1:NATOMSALLOC)=0.0D0 ! TO PREVENT READING FROM UNINITIALISED MEMORY
      IF (CALCQT) CALL SHINIT

      ALLOCATE(FINISH(3*NATOMSALLOC))

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      CALL KEYWORD
!     write out the atom indices to check:
!        IF(TWISTT) THEN
!              WRITE(MYUNIT,*)'NWISTGROUPS dihedrals to be constrained..'
!              DO J1=1, NTWISTGROUPS
!                WRITE(MYUNIT,*) 'Group J1  indices:'
!                WRITE(MYUNIT,*) TWIST_ATOMS(1:4, J1)
!                WRITE(MYUNIT,*) 'Group J1  k value, ref.dihedral:'
!                WRITE(MYUNIT,*) TWIST_K(J1), TWIST_THETA0(J1)
!              END DO
!        END IF

      IF (RIGIDINIT) THEN
          CALL GENRIGID_READ_FROM_FILE ()
! csw34> Tell the user how many degrees of freedom there are in the system
          WRITE(MYUNIT, '(A,I15)') " genrigid> rbodyconfig used to specifiy rigid bodies, degrees of freedom now ", DEGFREEDOMS
          IF (GCBHT) THEN
! csw34> Make sure we aren't running GCBH with rigid bodies - very dangerous! Could be done but only with great care...
              WRITE(MYUNIT, '(A)') " genrigid> ERROR: cannot use rigid bodies with GCBH! Stopping."
              STOP                
          END IF
      END IF

      IF (CUDAT) THEN
         IF ((CUDAPOT .EQ. 'A') .AND. (.NOT. AMBER12T)) THEN
            WRITE(MYUNIT,'(A)') " main> The AMBER12 keyword must be used with 'CUDA A'. "
            STOP
         END IF
         IF (DEBUG) THEN
            CUDAFILENAME1 ="GPU_debug_out"
            OPEN(UNIT=325, FILE=CUDAFILENAME1, STATUS="REPLACE")
            CLOSE(325)
         END IF
         IF (CUDATIMET) THEN
            CUDAFILENAME2 ="GPU_LBFGS_total.txt"
            CUDAFILENAME3 ="GPU_LBFGS_updates.txt"
            CUDAFILENAME4 ="GPU_LBFGS_linesearch.txt"
            CUDAFILENAME5 ="GPU_potential.txt"
            OPEN(UNIT=326, FILE=CUDAFILENAME2, STATUS="REPLACE")
            CLOSE(326)
            OPEN(UNIT=327, FILE=CUDAFILENAME3, STATUS="REPLACE")
            CLOSE(327)
            OPEN(UNIT=328, FILE=CUDAFILENAME4, STATUS="REPLACE")
            CLOSE(328)
            OPEN(UNIT=329, FILE=CUDAFILENAME5, STATUS="REPLACE")
            CLOSE(329)
         END IF
      END IF

      IF (GLJT.AND.SORTT) THEN
         WRITE(MYUNIT,'(A)') 'keyword> *** ERROR - you cannot specify SORT with GLJ'
         STOP
      ENDIF

!     related to FEBHT
      IF (FEBHT) THEN
         IF (.NOT. RIGID) THEN
            ALLOCATE(ATMASS(NATOMSALLOC)) 
         ELSE IF (NTIPT) THEN
            ALLOCATE(ATMASS(NATOMSALLOC/2*3))
         ELSE
            ALLOCATE(ATMASS(NATOMSALLOC/2*NRBSITES))
         ENDIF 
         ATMASS(:)=1.D0
         IF (NTIPT) THEN
            DO J1=1, NATOMSALLOC/2
               ATMASS(3*J1-2)=16.D0
               ATMASS(3*J1-1)=1.D0
               ATMASS(3*J1)  =1.D0
            ENDDO
         ENDIF
         IF (AMBERT) THEN
            ATMASS(:)=ATMASS1(:)
         ENDIF
         IF (AMBER12T .AND. ALLOCATED(AMBER12_MASSES)) THEN
            ATMASS(:)=AMBER12_MASSES(:)
         END IF
!
! The 2N! factor is accessed via commons and used in quench for febh.
!
         SYMFCTR=LOG(2.D0)
         IF (RIGID) THEN
            DO J1=1, NATOMSALLOC/2
               SYMFCTR=SYMFCTR+LOG(1.0D0*J1)
            ENDDO
         ELSE
            DO J1=1, NATOMSALLOC
               SYMFCTR=SYMFCTR+LOG(1.0D0*J1)
            ENDDO
         ENDIF
      ENDIF
      IF (READMASST) THEN
         INQUIRE(FILE='masses',EXIST=YESNO)
         IF (ALLOCATED(ATMASS)) DEALLOCATE(ATMASS)
         ALLOCATE(ATMASS(NATOMSALLOC)) 
         IF (YESNO) THEN
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='masses',STATUS='OLD')
            READ(LUNIT,*) (ATMASS(J1),J1=1,NATOMSALLOC)
            CLOSE(LUNIT)
            WRITE(MYUNIT,'(A)') ' fetchz> Masses replaced by values from masses file, Values are:'
            WRITE(MYUNIT,'(F20.10)') (ATMASS(J1),J1=1,NATOMSALLOC)
         ELSE
            WRITE(MYUNIT,'(A)') ' fetchz> ERROR *** keyword READMASS found, but no masses file'
            STOP
         ENDIF
      ENDIF
      IF (.NOT.ALLOCATED(ATMASS)) THEN
         IF (.NOT. RIGID) THEN
            ALLOCATE(ATMASS(NATOMSALLOC)) 
            ATMASS(1:NATOMSALLOC)=1.0D0
         ELSE IF (NTIPT) THEN
            ALLOCATE(ATMASS(NATOMSALLOC/2*3))
            ATMASS(1:NATOMSALLOC/2*3)=1.0D0
         ELSE
            ALLOCATE(ATMASS(NATOMSALLOC/2*NRBSITES))
            ATMASS(1:NATOMSALLOC/2*NRBSITES)=1.0D0
         ENDIF 
      ENDIF

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

!fh301>{{{
      IF (CHEMSHIFT) CALL CAMSHIFTREAD(20)
!fh301>}}}

      IF (CHRMMT) THEN
         IF (CHRIGIDTRANST.OR.CHRIGIDROTT) CALL CHSETSEG
! set charmm moves to standard (Mortensen & Wales) if no Cartesian moves are invoked via CHNMAX=0 
         IF ((.NOT.CHMOVET) .AND. (.NOT.CHCARTMOVET)) THEN
            CHSTANDARDT=.TRUE.
            CHBBT=.TRUE.
            CHSCT=.TRUE.
         ENDIF
         IF (SECPREDT) CALL CHSECPRED
         CALL CHSETPROB
      ENDIF

      IF (RMST) THEN
         ALLOCATE(RMSBEST(RMSSAVE,2),RMSCOOR(RMSSAVE,3*NATOMSALLOC))
         RMSBEST(1:RMSSAVE,1)=RMSLIMIT+RMSTOL
         RMSBEST(1:RMSSAVE,2)=0.D0
         RMSCOOR(1:RMSSAVE,1:3*NATOMSALLOC)=0.D0
         ALLOCATE(COORCOMP(1:3*NATOMSALLOC))
!
!        csw34> Need to add the read in for other file formats
!        (crd,pdb). Still not sure how GMIN reads in the input.crd file
!        for CHARMM. It appears to be reading into unit 7 in io1.f
!        though - but only from xyz format! Is there a hidden
!        conversion?
!         
         OPEN(UNIT=1,FILE='compare',STATUS='OLD')
         READ(1,*) (COORCOMP(J1),J1=1,3*NATOMS)
         CLOSE(1)
      ENDIF

      ALLOCATE(FF(NSAVE),QMIN(MAX(NSAVE,1)))
      ALLOCATE(QMINP(NSAVE,3*NATOMSALLOC))
      ALLOCATE(QMINNATOMS(NSAVE))
      ALLOCATE(AVOIDNATOMS(MAXSAVE))
      ALLOCATE(QMINT(NSAVE,NATOMSALLOC))
      ALLOCATE(NPCALL_QMIN(NSAVE))
      IF (MONITORT) THEN
         ALLOCATE(LOWESTC(3*NATOMSALLOC))
      ENDIF

!        csw34> ALLOCATE the interaction energy tracking arrays if A9INTE in data
      IF (A9INTET.AND.AMBERT) THEN 
         ALLOCATE(INTEFF(NSAVEINTE),INTEQMIN(NSAVEINTE))
         ALLOCATE(INTEQMINP(NSAVEINTE,3*NATOMS))
         INTEQMIN(:)=1.0D10 
         INTEQMINP(1:NSAVEINTE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
         INTEFF(1:NSAVEINTE)=0 ! to prevent reading from uninitialised memory
      ENDIF  

      IF (GAUSST) THEN
         ALLOCATE(GAUSSKK(3*NATOMSALLOC,GMODES),GAUSSEE(GMODES),GKSMALL(3*NATOMSALLOC))
         CALL KEGEN ! INITIAL SETUP 
         DO J1=1,GMODES
             PRINT *,J1,GAUSSEE(J1)
         ENDDO
      ENDIF

      QMINP(1:NSAVE,1:3*NATOMSALLOC)=0.0D0 ! to prevent reading from uninitialised memory
      QMINT(1:NSAVE,1:NATOMSALLOC)=1 ! to prevent reading from uninitialised memory
      QMINNATOMS(1:NSAVE)=NATOMS ! to prevent reading from uninitialised memory
      AVOIDNATOMS(1:MAXSAVE)=NATOMS ! as above
      COORDSO(1:3*NATOMSALLOC,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMSALLOC,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      ALLOCATE(ESAVE(NTAB,NPAR),XINSAVE(NTAB,NPAR))
      ALLOCATE(VEC(NVEC))

!     IF (SYMMETRIZE.AND.(.NOT.CENT)) THEN
!        PRINT '(A)','Probable input error - SYMMETRIZE true but CENT false'
!        STOP
!     ENDIF

      IF (DUMPT) THEN
         IF (CHRMMT) THEN
            OPEN(UNIT=719,FILE='dump.crd',STATUS='UNKNOWN')
            OPEN(UNIT=720,FILE='dump.pdb',STATUS='UNKNOWN')
         ENDIF
!
! dump.1.xyz is partly filled with control characters for no apparent reason.
! Suspect a compiler or mpi bug? DJW
! dump.1.xyz is fine with debug compilation, so it is a compiler bug!
!
         ALLOCATE(DUMPXYZUNIT(NPAR),DUMPVUNIT(NPAR))
         DO J1=1,NPAR
            WRITE (JSTR,'(I6)') J1

            ISTR='dump.' // TRIM(ADJUSTL(JSTR)) // '.V'
            DUMPVUNIT(J1)=GETUNIT()
            OPEN(UNIT=DUMPVUNIT(J1),FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')

            ISTR='dump.' // TRIM(ADJUSTL(JSTR)) // '.xyz'
            DUMPXYZUNIT(J1)=GETUNIT()
            OPEN(UNIT=DUMPXYZUNIT(J1),FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
         ENDDO
      ENDIF

      IF (PAIRDISTT) THEN
         MYPUNIT=3000+MYNODE
         IF (NPAR.GT.1) THEN
            OPEN(MYPUNIT,FILE="pairdists."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ELSE
            OPEN(MYPUNIT,FILE="pairdists",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ENDIF
         WRITE(MYPUNIT,'(A10)',ADVANCE="NO") "Quench  "
         DO J1=1,NPAIRS
            WRITE(atom1,*) PAIRDIST(J1,1)
            WRITE(atom2,*) PAIRDIST(J1,2)
            WRITE(atompair,*) TRIM(ADJUSTL(atom1))//"-"//TRIM(ADJUSTL(atom2))
            WRITE(MYPUNIT,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(atompair))//"  " 
         ENDDO
         WRITE(MYUNIT,'(A)') ""
      ENDIF

!
!     csw34> TRACKDATA keyword prints the energy and markov energy 
!     to files for viewing during a run. If RMS is also specified it
!     prints the rmsd from the comparison structure into a file.
!
      IF (TRACKDATAT) THEN
         MYEUNIT=4000+MYNODE
         MYMUNIT=6000+MYNODE
         MYBUNIT=8000+MYNODE
         MYRUNIT=10000+MYNODE
         IF (NPAR.GT.1) THEN
            IF (RESTORET) THEN
               OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               IF (RMST) OPEN(MYRUNIT,FILE="rmsd."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            ELSE
               OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               IF (RMST) OPEN(MYRUNIT,FILE="rmsd."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF
         ELSE
            IF (RESTORET) THEN
               OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               IF (RMST) OPEN(MYRUNIT,FILE='rmsd',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            ELSE
               OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED')
               IF (RMST) OPEN(MYRUNIT,FILE='rmsd',STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF
            IF (A9INTET) THEN
               OPEN(UNIT=3998,FILE='intE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(UNIT=3999,FILE='bestintE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF
         ENDIF
      ENDIF
      CALL FLUSH(6)
      CALL IO1
      IF (QCIPOTT) THEN
         ALLOCATE(TMPCOORDS(3*NATOMSALLOC))
         TMPCOORDS(1:3*NATOMSALLOC)=COORDS(1:3*NATOMSALLOC,1)
         IF (ALLOCATED(INTFROZEN)) DEALLOCATE(INTFROZEN)
         ALLOCATE(INTFROZEN(NATOMS))
         INTFROZEN(1:NATOMS)=.FALSE.
!        WRITE(MYUNIT, *) NATOMSALLOC
!        WRITE(MYUNIT, *) 'A tmpcoords before'
!        DO J1=1,NATOMSALLOC
!           WRITE(MYUNIT,'(A,3G20.10)') 'O ',TMPCOORDS(3*(J1-1)+1:3*(J1-1)+3)
!        ENDDO
         CALL CHECKPERC(TMPCOORDS,LINTCONSTRAINTTOL,0,1)
!        WRITE(MYUNIT, *) NATOMSALLOC
!        WRITE(MYUNIT, *) 'A tmpcoords after'
!        DO J1=1,NATOMSALLOC
!           WRITE(MYUNIT,'(A,3G20.10)') 'O ',TMPCOORDS(3*(J1-1)+1:3*(J1-1)+3)
!        ENDDO
         DEALLOCATE(TMPCOORDS)
      ENDIF
!
! Allocations for QCI
!
      IF (INTCONSTRAINTT.AND.(NCONGEOM.GE.2)) THEN
!
! Set up all the constraints and repulsions for zero frozen atoms.
!
         IF (.NOT.ALLOCATED(CONI)) THEN
            ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX))
            ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
         ENDIF

         ALLOCATE(TMPCOORDS(3*NATOMSALLOC))
         TMPCOORDS(1:3*NATOMSALLOC)=COORDS(1:3*NATOMSALLOC,1)
         DEALLOCATE(TMPCOORDS)
         INTFREEZETOLSAVE=INTFREEZETOL
         INTFREEZETOL=-1.0D0
         CALL MAKE_CONPOT(NCONGEOM,CONGEOM)
         INTFREEZETOL=INTFREEZETOLSAVE
!
! Now align the two endpoints to the first reference minimum.
!
         
         ALLOCATE(ENDCOORDS(2,3*NATOMSALLOC))
         ENDCOORDS(1,1:3*NATOMSALLOC)=COORDS(1:3*NATOMSALLOC,1)
         ENDCOORDS(2,1:3*NATOMSALLOC)=FINISH(1:3*NATOMSALLOC)
         CALL MAKE_CONPOT(2,ENDCOORDS)
         COORDS(1:3*NATOMSALLOC,1)=ENDCOORDS(1,1:3*NATOMSALLOC)
         FINISH(1:3*NATOMSALLOC)=ENDCOORDS(2,1:3*NATOMSALLOC)
         DEALLOCATE(ENDCOORDS)
      ENDIF
!
! If this is a CSM optimisation we now have to multiply the number of atoms by the number of
! group operations and replicate some coordinates and allowed permutations.
!
      IF (CSMT) THEN
         CALL CSMINIT
         IF (SYMMETRIZECSM) THEN
            IF (CSMMAXIT.EQ.0) CSMMAXIT=MAXIT
         ELSE
            CSMMAXIT=MAXIT
         ENDIF
      ENDIF

      CALL FLUSH(6)

      IF (PROJIT) THEN
         ALLOCATE(TMPCOORDS(3*NATOMSALLOC))
         TMPCOORDS(1:3*NATOMSALLOC)=COORDS(1:3*NATOMSALLOC,1)
         CALL PROJIINIT(TMPCOORDS,NATOMSALLOC)
         DEALLOCATE(TMPCOORDS)
      ENDIF


      IF(VGW) THEN
        CALL INITIALIZE_VGWSP(NATOMSALLOC,LJSIGMA,LJEPSILON,TAUMAX,CPFACTORSG,CPS,VGWTOL)
        CALL INITIALIZE_VGW(NATOMSALLOC,LJSIGMA,LJEPSILON,TAUMAXFULL,CPFACTORFG,CPF)
      ENDIF

      IF (SEEDT) THEN
         CALL GSEED
      ELSE
         IF ((.NOT.FIELDT).AND.CENT) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRE2(COORDS(1:3*NATOMS,J1))
            ENDDO
         ELSEIF ((.NOT.FIELDT).AND.FIXCOM) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRECOM(COORDS(1:3*NATOMS,J1))
            ENDDO
         ENDIF
      ENDIF

      IF (SUPERSTEP) NSUPERSTEP=0
      
      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
         NPCALL_QMIN(J1)=0
      ENDDO

      CALL INITIALIZATIONS()

! vr274>  call test routines.
!         By default this routine does nothing. In the testsuite, it is replaced by
!         a customized function which can run checks whether everything runs fine.
      CALL RUN_TESTS_AFTER_INIT()

!mo361> Run GA
      IF (GENALT) THEN
         NRUNS=0
         CALL MYGA_RUN()
      ENDIF
!ds656> Span label permutations
      IF (MULTIPERMT) THEN
         NRUNS=0
         CALL MULTIPERM()
      ENDIF

      IF ((NRUNS.GT.0).OR.PTMC.OR.BSPT) CALL MCRUNS(SCREENC)
C     CALL SYSTEM('rm ssdump ssave >& /dev/null')

      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      IF (ALLOCATED(XICOM)) DEALLOCATE(XICOM)
      IF (ALLOCATED(PCOM)) DEALLOCATE(PCOM)
      IF (ALLOCATED(GAUSSKK)) DEALLOCATE(GAUSSKK,GAUSSEE)
C     deallocate dynamic memory for AMBER
      IF (AMBERT) CALL AMBER_DEALLOCATE_STACKS
      IF (ALLOCATED(ANV)) DEALLOCATE(ANV)      
      CALL FLUSH(MYUNIT)

      IF(VGW) THEN
        CALL CLEANUP_VGWSP
        CALL CLEANUP_VGW
      ENDIF

      CLOSE(MYUNIT)
    
! khs26> Close file storing free energies
      IF (FEBHT) THEN
         CLOSE(FE_FILE_UNIT)
      END IF
! csw34> close pairdists.* files
      IF (PAIRDISTT) CLOSE(MYPUNIT)
      IF (TRACKDATAT) THEN
         CLOSE(MYEUNIT)
         CLOSE(MYMUNIT)
         IF (RMST) CLOSE(MYRUNIT)
         CLOSE(MYBUNIT)
         IF (A9INTET) THEN
            CLOSE(3998)
            CLOSE(3999)
         ENDIF
      ENDIF
#ifdef MPI
      CALL MPI_FINALIZE(MPIERR)
#endif
      DEALLOCATE(SCREENC)
      DEALLOCATE(IATNUM, VT, ZSYM)
      DEALLOCATE(FF,QMIN)
      DEALLOCATE(QMINP,QMINT,QMINNATOMS)
      DEALLOCATE(ESAVE,XINSAVE)
      DEALLOCATE(VEC)
      DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,NQ,EPREV,
     &           JUMPMOVE,JUMPINT,JDUMP,COORDS,COORDSO,VAT,VATO,
     &         JUMPTO,SHELLMOVES,PTGROUP,NSURFMOVES,NCORE)
      DEALLOCATE(FROZEN)
      DEALLOCATE(FROZENRES)
      DEALLOCATE(DONTMOVE)
      DEALLOCATE(DONTMOVERES)
      DEALLOCATE(HARMONICFLIST)
      DEALLOCATE(HARMONICR0)
      DEALLOCATE(BESTPERM)
      IF (ALLOCATED(ANV)) DEALLOCATE(ANV)
      IF (ALLOCATED(NSPECIES)) DEALLOCATE(NSPECIES)
      IF (ALLOCATED(NSPECIES_INI)) DEALLOCATE(NSPECIES_INI)
      IF (ALLOCATED(SPECMASS)) DEALLOCATE(SPECMASS)
      IF (ALLOCATED(MSBCOORDS)) DEALLOCATE(MSBCOORDS)
      IF (ALLOCATED(MSBE)) DEALLOCATE(MSBE)
      IF (ALLOCATED(BEADLETTER)) DEALLOCATE(BEADLETTER)
      IF (ALLOCATED(BLNSSTRUCT)) DEALLOCATE(BLNSSTRUCT)
      IF (ALLOCATED(NPERMSIZE)) DEALLOCATE(NPERMSIZE,PERMGROUP,NSETS,SETS)
      IF (ALLOCATED(HDIST)) DEALLOCATE(HDIST,HWEIGHT,HISTVALS,LHISTVALS,IGNOREBIN)
      IF (ALLOCATED(SITE)) DEALLOCATE(SITE)
      IF (ALLOCATED(MOVABLEATOMLIST)) DEALLOCATE(MOVABLEATOMLIST)
      IF (ALLOCATED(MOVABLEATOMLISTLOGICAL)) DEALLOCATE(MOVABLEATOMLISTLOGICAL)
      IF (ALLOCATED(ATOMSINBLOCK)) DEALLOCATE(ATOMSINBLOCK)
      IF (ALLOCATED(COORDS)) DEALLOCATE(COORDS)
      IF (ALLOCATED(DATOMGROUPNAMES)) DEALLOCATE(DATOMGROUPNAMES)
      IF (ALLOCATED(PYA1bin)) DEALLOCATE(PYA1bin)
      IF (ALLOCATED(PYA2bin)) DEALLOCATE(PYA2bin)
      IF (ALLOCATED(PYA1bin)) DEALLOCATE(PYA1bin)
      IF (ALLOCATED(PYA2bin)) DEALLOCATE(PYA2bin)
      IF (ALLOCATED(TARGETS)) DEALLOCATE(TARGETS)
      IF (ALLOCATED(TCOORDS)) DEALLOCATE(TCOORDS)
      IF (ALLOCATED(TBPBASINS)) DEALLOCATE(TBPBASINS)
      IF (ALLOCATED(SEMIGRAND_MU)) DEALLOCATE(SEMIGRAND_MU)
      IF (ALLOCATED(STRESS)) DEALLOCATE(STRESS)
      IF (ALLOCATED(ATOMLISTS)) DEALLOCATE(ATOMLISTS)
      IF (ALLOCATED(INVATOMLISTS)) DEALLOCATE(INVATOMLISTS)
      IF (ALLOCATED(NBRLISTS)) DEALLOCATE(NBRLISTS)
      IF (ALLOCATED(NNLISTS)) DEALLOCATE(NNLISTS)
      IF (ALLOCATED(VSITES)) DEALLOCATE(VSITES)
      IF (ALLOCATED(MIEF_SIG)) DEALLOCATE(MIEF_SIG)
      IF (ALLOCATED(MIEF_EPS)) DEALLOCATE(MIEF_EPS)
      IF (ALLOCATED(MIEF_SITES)) DEALLOCATE(MIEF_SITES)
      IF (ALLOCATED(MIEF_U_RCUT)) DEALLOCATE(MIEF_U_RCUT)
      IF (ALLOCATED(MIEF_DUDR_RCUT)) DEALLOCATE(MIEF_DUDR_RCUT)
      IF (ALLOCATED(NNBRS)) DEALLOCATE(NNBRS)
      IF (ALLOCATED(NNHIST)) DEALLOCATE(NNHIST)
      IF (ALLOCATED(ANNHIST)) DEALLOCATE(ANNHIST)
      IF (ALLOCATED(ANNHIST_MEAN)) DEALLOCATE(ANNHIST_MEAN)
      IF (ALLOCATED(ANNHIST_VAR)) DEALLOCATE(ANNHIST_VAR)
      IF (ALLOCATED(NNBOND_WEIGHTS)) DEALLOCATE(NNBOND_WEIGHTS)
      IF (ALLOCATED(IFLIPE)) DEALLOCATE(IFLIPE)
      IF (ALLOCATED(OBJ)) DEALLOCATE(OBJ)
      IF (ALLOCATED(TWIST_K)) DEALLOCATE(TWIST_K)
      IF (ALLOCATED(TWIST_THETA0)) DEALLOCATE(TWIST_THETA0)
      IF (ALLOCATED(TWIST_ATOMS)) DEALLOCATE(TWIST_ATOMS)

!
! Should also deallocate any of the arrays below, which may have been
! allocated in keywords. This is just to make nagfmcheck happy, though! DJW
!
!          ALLOCATE(DATOMGROUPAXIS(NDGROUPS,2))
!          ALLOCATE(DATOMGROUPPROBA(NDGROUPS))
!          ALLOCATE(DATOMGROUPSCALINGA(NDGROUPS))
!          ALLOCATE(DATOMGROUPPROBB(NDGROUPS))
!          ALLOCATE(DATOMGROUPSCALINGB(NDGROUPS))
!          ALLOCATE(DATOMGROUPPROBC(NDGROUPS))
!          ALLOCATE(DATOMGROUPSCALINGC(NDGROUPS))
!          ALLOCATE(DATOMGROUPPROBD(NDGROUPS))
!          ALLOCATE(DATOMGROUPSCALINGD(NDGROUPS))
!          ALLOCATE(DATOMGROUPS(NDGROUPS,NATOMS))
!          ALLOCATE(DATOMGROUPATT(NDGROUPS,NATOMS))
!          ALLOCATE(ATOMGROUPNAMES(NGROUPS))
!          ALLOCATE(ATOMGROUPAXIS(NGROUPS,2))
!          ALLOCATE(ATOMGROUPPSELECT(NGROUPS))
!          ALLOCATE(ATOMGROUPSCALING(NGROUPS))
!          ALLOCATE(ATOMGROUPS(NGROUPS,NATOMS))
!          ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
!          ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
!          ALLOCATE(SHELLMOVES(NPAR))
!          ALLOCATE(PTGROUP(NPAR))
!          ALLOCATE(NSURFMOVES(NPAR))
!          ALLOCATE(NCORE(NPAR))
!          ALLOCATE(PAIRDIST(NPAIRS,2))
!          ALLOCATE(RBOPS(4,NRBGROUP))
!          ALLOCATE(EMIN(NRESMIN),FVIBMIN(NRESMIN),PFMIN(NRESMIN),IXMIN(NRESMIN),IYMIN(NRESMIN),IZMIN(NRESMIN),HORDERMIN(NRESMIN))
!          ALLOCATE(RESPOINTS(3*NATOMS,NRESMIN))
!          ALLOCATE(RBUV(NRBSITES,3))
!          ALLOCATE(DPMU(NRBSITES))
!          ALLOCATE(RBUV(NRBSITES,3))
!          ALLOCATE(DPMU(NRBSITES))
!          ALLOCATE(RBUV(NRBSITES,3))
!          ALLOCATE(NCMP(TPAHA))
!          ALLOCATE(RBSTLA(NRBSITES,3))
!          ALLOCATE(STCHRG(NRBSITES))
!          ALLOCATE(RBSTLA(NRBSITES,3))
!          ALLOCATE(RBSTLA(NRBSITES,3))

      STOP

      END
