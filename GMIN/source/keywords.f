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
      SUBROUTINE KEYWORD
!      USE NEW_INPUT_MOD, ONLY : NEW_INPUT
      use COMMONS
      use VEC3
      use genrigid
      use MODMXATMS   ! NEEDED FOR charmm
      USE modcharmm
      USE TWIST_MOD
!       sf344> AMBER additions
      USE modamber9, only : coords1,amberstr,amberstr1,mdstept,inpcrd,amberenergiest, nocistransdna, nocistransrna,
     &                      uachiral, ligrotscale, setchiral, STEEREDMINT, SMINATOMA, SMINATOMB, SMINK, SMINKINC,
     &                      SMINDISTSTART, SMINDISTFINISH, natomsina, natomsinb, natomsinc, atomsinalist, atomsinblist,
     &                      atomsinclist, atomsinalistlogical, atomsinblistlogical, atomsinclistlogical, ligcartstep,
     &                      ligtransstep, ligmovefreq, amchnmax, amchnmin, amchpmax, amchpmin, rotamert, rotmaxchange,
     &                      rotcentre, rotpselect, rotoccuw, rotcutoff, setchiralgeneric, PRMTOP, IGB, RGBMAX, CUT,
     &                      SALTCON
      USE modamber
      USE PORFUNCS
      USE MYGA_PARAMS
      USE BGUPMOD
      USE GLJYMOD
      USE CHIRO_MODULE, ONLY: CHIRO_SIGMA, CHIRO_MU, CHIRO_GAMMA, CHIRO_L
      USE MBPOLMOD, ONLY: MBPOLINIT
      USE SWMOD, ONLY: SWINIT, MWINIT
      USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_SETUP, AMBER12_GET_COORDS, AMBER12_ATOMS,
     &                                  AMBER12_RESIDUES, POPULATE_ATOM_DATA
      USE CHIRALITY, ONLY : CIS_TRANS_TOL
      USE ISO_C_BINDING, ONLY: C_NULL_CHAR
      USE PARSE_POT_PARAMS, ONLY : PARSE_MGUPTA_PARAMS, PARSE_MSC_PARAMS, 
     &     PARSE_MLJ_PARAMS
      USE ROTAMER, ONLY: ROTAMER_MOVET, ROTAMER_SCRIPT, ROTAMER_INIT
      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, LAST, IX, J1, JP, NPCOUNT, NDUMMY, INDEX, J2, J3, J4
      INTEGER DATA_UNIT
      INTEGER MOVABLEATOMINDEX
      LOGICAL CAT, YESNO, PERMFILE, CONFILE
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT
       DOUBLE PRECISION XX, ROH, ROM, WTHETA 
      LOGICAL END, SKIPBL, CLEAR, ECHO
      CHARACTER WORD*16,PBC*3,WORD2*10
      DOUBLE PRECISION EAMLJA0, EAMLJBETA, EAMLJZ0, DUMMY
      COMMON /EAMLJCOMM/ EAMLJA0, EAMLJBETA, EAMLJZ0
      DOUBLE PRECISION SLENGTH, EPS
      INTEGER NOK, NBAD
      COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      DOUBLE PRECISION EPS2, RAD, HEIGHT
      COMMON /CAPS/ EPS2, RAD, HEIGHT

!     LOGICAL IGNOREBIN(HISTBINMAX), FIXBIN
!     COMMON /IG/ IGNOREBIN, FIXBIN
!      DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP
!      COMMON /AMBWORD/    PMAX,PMIN,NMAX,NMIN,SIDESTEP

      INTEGER NATOM, DMODE, NDUM, NONTYPEA
!
! These arrays should have dimension MXATMS
!
      DOUBLE PRECISION, ALLOCATABLE :: CHX(:), CHY(:), CHZ(:), CHMASS(:)
      CHARACTER(LEN=1) DUMMYCH
      CHARACTER(LEN=100) TOPFILE,PARFILE
      CHARACTER(LEN=20) UNSTRING
      DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      DOUBLE PRECISION LJREPBL, LJATTBL, LJREPBN, LJATTBN, LJREPLN, LJATTLN

!     DC430 >
      DOUBLE PRECISION :: LPL, LPR
      LOGICAL          :: RBSYMTEST     ! jdf43>
!
!       sf344> added stuff
!
      CHARACTER(LEN=10) check1
      CHARACTER(LEN=1) readswitch
      CHARACTER(LEN=4) J1CHAR
      CHARACTER(LEN=20) J2CHAR
      INTEGER iostatus, groupsize, groupatom,groupoffset,axis1,axis2,EOF
      INTEGER LUNIT, GETUNIT

      !ab2111> dihedral rotation stuff
      INTEGER dihedraloffset,dihedralgroupsize,A1,A2,A3,A4

! ab2111 > Reservoir stuff
      INTEGER VECLINE, VECNUM, RESCOUNT

! hk286 - DAMPED GROUP MOVES
      DOUBLE PRECISION GROUPATT

      CHARACTER(LEN=120), DIMENSION(:,:), ALLOCATABLE :: KEY_WORDS
      DOUBLE PRECISION DUMMY1(NATOMS)

!      OPEN(5120, FILE = 'data')
!      CALL NEW_INPUT(5120, KEY_WORDS)
!      CLOSE(5120)
!      PRINT *, KEY_WORDS

      AAA=0
      AAB=0
      ABB=0
      PAA=0
      PAB=0
      PBB=0
      QAA=0
      QAB=0
      QBB=0
      ZAA=0
      ZAB=0
      ZBB=0
      R0AA=0
      R0AB=0
      R0BB=0

      NPCOUNT=0
      NPCALL=0
      NSEED=0
      NS=0
      NSSTOP=0
      HIT=.FALSE.
      SAVEQ=.TRUE.
      NSAVE=5
      NSAVEINTE=0
      SAVEMULTIMINONLY=.FALSE.
      TFAC(:)=1.0D0
      RESIZET=.FALSE.
      STEPOUT=.FALSE.
      SUPERSTEP=.FALSE.
      NSUPER=10
      SUPSTEP=1.1D0
      SACCRAT=0.5D0
      NSACCEPT=100
      EVSTEPT=.FALSE.
      NEVS=100
      CEIG=0.1D0
      NEVL=100
      NVECTORS=2
      TEMPS=0.8
      NRBSITES=0
      CHFREQ=1
      CHFREQSC=1
      CHFREQBB=1
      FTRANS=1
      FROT=1

!
! QCI parameters
!
         CONDATT=.FALSE.
         QCIPOTT=.FALSE.
         QCIPOT2T=.FALSE.
         FREEZETOL=1.0D-3
         INTCONSTRAINTT=.FALSE.
         INTCONSTRAINTTOL=0.1D0
         INTCONSTRAINTDEL=10.0D0
         INTCONSTRAINTREP=100.0D0
         INTCONSTRAINREPCUT=1.7D0
         INTFREEZET=.FALSE.
         INTFREEZETOL=1.0D-3
         INTFREEZEMIN=10
         INTCONFRAC=0.9D0
         INTCONSEP=15
         INTREPSEP=0
         INTSTEPS1=300001
         INTCONSTEPS=100
         INTRELSTEPS=200
         MAXCONUSE=4
         MAXCONE=0.01D0
         INTRMSTOL=0.01D0
         IMSEPMIN=0.2D0
         IMSEPMAX=10.0D0
         INTIMAGE=3
         MAXINTIMAGE=75
         INTNTRIESMAX=2
         INTIMAGEINCR=6
         INTIMAGECHECK=25
         IMSEPMIN=0.0D0
         IMSEPMAX=HUGE(1.0D0)

         CHECKCONINT=.FALSE.
         CONCUTABS=0.15D0
         CONCUTABST=.TRUE.
         CONCUTFRAC=0.1D0
         CONCUTFRACT=.FALSE.
         CHECKREPINTERVAL=10
         CHECKREPCUTOFF=2.0D0
         DUMPINTXYZ=.FALSE.
         DUMPINTEOS=.FALSE.
         DUMPINTXYZFREQ=100
         DUMPINTEOSFREQ=100
         KINT=0.0D0

! ns644>    Adding the TWIST keyword:
      TWISTT=.FALSE.
      NTWISTGROUPS=0


      ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1),TEMP(1),ACCRAT(1),STEP(1),ASTEP(1),OSTEP(1),BLOCK(1),NT(1),NQ(1),EPREV(1),
     @         JUMPMOVE(1),JUMPINT(1),JDUMP(1),COORDS(3*NATOMSALLOC,1),COORDSO(3*NATOMSALLOC,1),VAT(NATOMSALLOC,1),
     @         VATO(NATOMSALLOC,1),  
     @         JUMPTO(1),SHELLMOVES(1),PTGROUP(1),NSURFMOVES(1),NCORE(1))
      DO JP=1,1
         EPREV(JP)=1.0D100
         FIXSTEP(JP)=.FALSE.
         FIXTEMP(JP)=.FALSE.
         FIXBOTH(JP)=.FALSE.
         TEMP(JP)=0.3D0
         ACCRAT(JP)=0.5D0
         STEP(JP)=0.3D0
         ASTEP(JP)=0.3D0
         OSTEP(JP)=0.3D0
         BLOCK(JP)=0
         NT(JP)=0
         JUMPMOVE(JP)=.FALSE.
         JUMPINT(JP)=100
         JDUMP(JP)=.FALSE.
         SHELLMOVES(JP)=.FALSE.
         PTGROUP(JP)='    '
         NSURFMOVES(JP)=0
         NCORE(JP)=0
      ENDDO
      NEWJUMP=.FALSE.
      PNEWJUMP=0.2D0
      ECONV=0.02D0
      TABOOT=.FALSE.
      NTAB=10
      CUTOFF=1.0D6
      PCUTOFF=1.0D6
      FINALCUTOFF=1.0D6
      MYPOWER=5
      NEON=.FALSE.
      RGCL2=.FALSE.
      AXTELL=.FALSE.
      ZSTAR=0.0D0
      GROUND=.FALSE.
      ARGON=.FALSE.
      ARNO=.FALSE.
      STAR=.FALSE.
      PLUS=.FALSE.
      TWOPLUS=.FALSE.
      DIPOLE=.FALSE.
      DUMPT=.FALSE.
      TARGET=.FALSE.
      SORTT=.FALSE.
      NTARGETS=0
      MSORIGT=.FALSE.
      MSTRANST=.FALSE.
      FRAUSIT=.FALSE.
      ANGST=.FALSE.
      MORSET=.FALSE.
      LB2T=.FALSE.
      DZTEST=.FALSE.
      ZETT1=.FALSE.
      ZETT2=.FALSE.
      P46=.FALSE.
      G46=.FALSE.
      BLNT=.FALSE.
      CHAPERONINT=.FALSE.
      DFTBT=.FALSE.
      DFTBCT=.FALSE.
      LJATT=.FALSE.
      SW=.FALSE.
      QUIPT =.FALSE.
      XMUL=1
      SCT=.FALSE.
      MSCT=.FALSE.
      MGUPTAT=.FALSE.
      SQUEEZET=.FALSE.
      NVEC=0
!     SQUEEZER=5.0D0
!     SQUEEZED=0.95D0
      DEBUG=.FALSE.
      DEBUGss2029=.FALSE.
      SEEDT=.FALSE.
      FREEZECORE=.TRUE.
      FREEZE=.FALSE.
      FREEZERES=.FALSE.
      FREEZEALL=.FALSE.
      UNFREEZERES =.FALSE.
! sf344> unfreeze structures at the final quench
      UNFREEZEFINALQ=.FALSE.
      NFREEZE=0
      ALLOCATE(FROZEN(NATOMSALLOC))
      FREEZESAVE=.TRUE.
! csw34> The FROZENRES array is bigger than needed
      ALLOCATE(FROZENRES(NATOMSALLOC))
      DO J1=1,NATOMSALLOC
         FROZEN(J1)=.FALSE.
         FROZENRES(J1)=.FALSE.
      ENDDO
      FREEZEGROUPTYPE='GT'
      FREEZEGROUPT=.FALSE.
! jdf43> FROZENRIGIDBODY
      ALLOCATE(FROZENRIGIDBODY(NATOMSALLOC))
      FROZENRIGIDBODY(:)=.FALSE.
! csw34> DONTMOVE defaults
      DONTMOVET=.FALSE.
      NDONTMOVE=0
      DONTMOVEREST=.FALSE.
      DONTMOVEALL=.FALSE.
      DOMOVEREST=.FALSE.
      DONTMOVEGROUPT=.FALSE.
      DONTMOVEGROUPTYPE='GT'
      ALLOCATE(DONTMOVE(NATOMSALLOC))
      ALLOCATE(DONTMOVERES(NATOMSALLOC))
      DO J1=1,NATOMSALLOC
         DONTMOVE(J1)=.FALSE.
         DONTMOVERES(J1)=.FALSE.
      ENDDO
! END of DONTMOVE defaults
! csw34> flag for new moves module
      NEWMOVEST=.FALSE.
      RESTRICTREGION=.FALSE.
      RESTRICTCYL=.FALSE.
      OVERLAPK=.FALSE.
      OVERLAP_IMPORT=.FALSE.
      ONE_ATOM_TAKESTEP=.FALSE.
      HARMONICF=.FALSE.
      HARMONICDONTMOVE=.FALSE.
      ALLOCATE(HARMONICFLIST(NATOMSALLOC))
      ALLOCATE(HARMONICR0(3*NATOMSALLOC))
      DO J1=1,NATOMSALLOC
         HARMONICFLIST(J1)=.FALSE.
      ENDDO
      PTREADTEMPS = .FALSE.
      FREEZEIL_USE=.FALSE.
      CHECKCHIRALITY=.TRUE.
      NOCISTRANS=.TRUE.
      NOCISTRANSRNA=.FALSE.
      NOCISTRANSDNA=.FALSE.
      MINOMEGA=150.D0
      CIS_TRANS_TOL=180.0D0 - MINOMEGA
      UACHIRAL=.FALSE.
      SETCHIRAL=.FALSE.
      SETCHIRALGENERIC=.FALSE.
      !
      !ds656> box centroid
      BOXCENTROIDT=.FALSE.
      BOXCENTROID_DISCRETE(1:3)=.FALSE.
      BOXCENTROID_X(1:3)=0.0D0
      BOXCENTROID_DX(1:3)=0.0D0
      !
      !ds656> substrate field(s)
      MIEFT=.FALSE.
      MIEF_PBCT=.FALSE.
      MIEF_CUTT=.FALSE.
      MIEF_BOX(1:3) = 1.0D9
      MIEF_RCUT= 1.0D9
      !
      FIELDT=.FALSE.
      OHT=.FALSE.
      IHT=.FALSE.
      TDT=.FALSE.
      D5HT=.FALSE.
      CENT=.FALSE.
      CENTXY=.FALSE.
      SETCENT=.FALSE.
      CENTX=0.0D0
      CENTY=0.0D0
      CENTZ=0.0D0
      QUCENTRE=.FALSE.
      FIXCOM=.FALSE.
      FIH=0.0D0
      FTD=0.0D0
      FD5H=0.0D0
      TOLD=0.0001D0
      TOLE=0.0001D0
      CUTT=.FALSE.
      PERIODIC=.FALSE.
      ! ds656> initialise box dimensions to very large
      BOXLX = 1.0D6
      BOXLY = 1.0D6
      BOXLZ = 1.0D6
      BOX3D(1:3) = (/BOXLX,BOXLY,BOXLZ/)
      PARAMONOVPBCX=.FALSE.
      PARAMONOVPBCY=.FALSE.
      PARAMONOVPBCZ=.FALSE.
      PARAMONOVCUTOFF=.FALSE.
      LJSITE=.FALSE.
      BLJSITE=.FALSE.
      LJSITECOORDST=.FALSE.
      LJSITEATTR=.FALSE.
      NRUNS=0
      PCUT=1.0D0
      RADIUS=0.0D0
      MAXIT=500
      MAXIT2=500
      EXPFAC=10.0D0
      EXPD=1.0D0
      CQMAX=1.0D-10
      BQMAX=1.0D-3
      RHO=6.0D0
      NACCEPT=50
      NORESET=.FALSE.
      TSALLIST=.FALSE.
      NEWTSALLIST=.FALSE.
      QTSALLIS=0.0D0
      PARALLELT=.FALSE.
      TOSI=.FALSE.
      WELCH=.FALSE.
      SOFT_SPHERE = .FALSE.
      SOFT_SPHERE_NTYPEA = 0
      BINARY=.FALSE.
      BINARY_EXAB=.FALSE.
      SHIFTCUT=.FALSE.
      FAL=.FALSE.
      FNI=.FALSE.
      PHI4MODELT=.FALSE.
      
      READMASST=.FALSE.
      SPECMASST=.FALSE.

!     AMBER=.FALSE.
      AMHT=.FALSE.
      NINT_AMH=1
      DPARAM=1.0D0
      FAKEWATER=.FALSE.
      AMCUT= .FALSE.
      MGBWATER=.FALSE.
      BIN=.FALSE.
      AMBERSEED= .FALSE.
      FIXT= .FALSE.
      FIX= .FALSE.
      CAP= .TRUE.
      WATERSTEP= .FALSE.
      QCUTOFF= 1.0D6
      RCUTOFF= 1.0D6
      REALQCUTOFF= 1.0D6
      REALRCUTOFF= 1.0D6
      RINGROTSCALE=0.0D0
      TRACKDATAT=.FALSE.
      PROGRESS=.FALSE.
      listupdate=20

      BLJCLUSTER=.FALSE.
      BLJCLUSTER_NOCUT=.FALSE.
      BGUPTAT=.FALSE.
      BGUPTATAB=.FALSE.
      BGUPTATBB=.FALSE.
      LSWAPST=.FALSE.
      LSWAPS_N1 = 1
      LSWAPS_N2 = 1
      LSWAPS_TEMP = 0.1
      LSWAPS_TFACTOR = 1.0
      LSWAPS_NUP = 1
      LFLIPST=.FALSE.
      LFLIPS_RESET=.FALSE.
      LFLIPS_N1 = 1
      LFLIPS_N2 = 1
      LFLIPS_TEMP = 0.1
      LFLIPS_TFACTOR = 1.0
      LFLIPS_NUP = 1
      BASWAP=.FALSE.
      BASWAP_FRAC=0.0D0
      !BASWAP_TEMP=TEMP(1) ! no longer used
      BASWAP_NWAIT=0
      BASWAP_NMAX=0
      BASWAPTEST=.FALSE.
      NTYPEA = NATOMS ! ds656> before this was initialised to 0
      BGUPTANAME1 = 'Au'
      BGUPTANAME2 = 'Ag'
      HOMOREFT = .FALSE.
      HOMOREFTEST = .FALSE.
      HOMOREFCHECK = .FALSE.
      HOMOREF_LSMODE = 0
      HOMOREF_FGMODE = 0
      HOMOREF_NCYCLES = 0
      HOMOREF_NFMAX = 3
      HOMOREF_NSMAX = 1
      HOMOREF_AUXT = .FALSE.
      HOMOREF_AUX_NSWAPS = 0
      HOMOREF_AUX_TEMP = 1.0D0
      HOMOREF_AUX_FACTOR = 1.0D0
      HOMOREF_AUX_NNCUT = 1.5D0
      HOMOREF_BHT = .FALSE.
      HOMOREF_BH_NSWAPMAX = 0
      HOMOREF_BH_NDUDMAX = 0
      HOMOREF_BH_TEMP = 1.0D0
      HOMOREF_BH_FACTOR = 1.0D0
      QALCST = .FALSE.
      QALCSV = .FALSE.
      QALCSMODE = 0
      QALCS_PARAM = 0
      QALCS_SURFT = .FALSE.
      QALCS_SYMT = .FALSE.
      QALCS_SYM_MINCORESIZE = 5
      SPECLABELST = .FALSE.
      RANDMULTIPERMT = .FALSE.
      RANDMULTIPERM_STEP=1
      RANDPERMT = .FALSE.
      RANDPERM_STEP=1
      MULTIPERMT = .FALSE.
      SPANSWAPST = .FALSE.
      ENPERMST = .FALSE.
      KEEPLISTS = .FALSE.
      NBRCUT1 = 0.0D0
      NBRCUT2 = 0.0D0

!     ds656> Stress tensor calculation
      STRESST=.FALSE.
      STRESS_MODE=1

!     ds656> Generalised LJ with Yukawa
      GLJY = .FALSE.
      GLJ_EXP = 6
      YUK_A = 0.0D0
      YUK_XI = 1.0D0
      
      CHRMMT=.FALSE.
      CHARMMTYPE=1
      CHARMMDFTBT=.FALSE.
      ACESOLV=.FALSE.
      ACEUPSTEP=50
      CHRIGIDTRANST=.FALSE.
      CHRIGIDROTT=.FALSE.
      CHNMIN=0.D0
      CHNMAX=HUGE(1.0D0)
      CHMDT=.FALSE.
      CHMDFREQ=HUGE(1)
      CURRENTIMP=0
      BOXT=.FALSE.
      SPHERET=.FALSE.
      RMST=.FALSE.
      NEWCONFT=.FALSE.
      INTMINT=.FALSE.
      INTSPRINGACTIVET=.FALSE.
      INTMINFAC=1.0D0
      DAESTAT=.FALSE.
      MAKEOLIGOT=.FALSE.
      MAKEOLIGOSTART=.FALSE.
      TRANSXYT=.FALSE.
      ROTZT=.FALSE.
      NREPEAT=0
      NFIXSEG=0
      OHCELLT=.FALSE.

! khs26> AMBER12 stuff
      AMBER12T=.FALSE.

!  sf344> AMBER stuff
      AMBERT=.FALSE.
      AMCHNMAX=0.0D0
      AMCHPMAX=0.0D0
      MDSTEPT=.FALSE.
      DUMPSTRUCTURES=.FALSE.
      ENERGY_DECOMPT=.FALSE.
      DUMPMINT=.FALSE.
! allocate ATOMSINBLOCK to avoid reading unallocated memory
      ALLOCATE(ATOMSINBLOCK(1))
! csw34> RANDOMSEED now works for CHARMM also!
      RANDOMSEEDT=.FALSE.
! csw34> Dumping structures after every quench
      DUMPQUT=.FALSE.
! csw34> Dumping structures after every step (before quenching)
      DUMPSTEPST=.FALSE.
! khs26> Dump best structures after every step
      DUMPBESTT=.FALSE.
! csw34> Local sampling within distance constraints 
      LOCALSAMPLET=.FALSE. 
      ABTHRESH=999.99
      ACTHRESH=999.99
! csw34> AMBER interaction energy logical
      A9INTET=.FALSE.
      INTERESTORE=.FALSE.
! csw34> set COLDFUSION flag to .FALSE.
      COlDFUSION=.FALSE.

      GRADPROBLEMT=.FALSE.

!
!  sf344> for specifying movable atoms and ligand rotating steptaking moves
!
      MOVABLEATOMST=.FALSE.
      LIGMOVET=.FALSE.
      LIGROTSCALE=0.0D0
      LIGCARTSTEP=0.0D0
      LIGTRANSSTEP=0.0D0
      LIGMOVEFREQ=1
!  sf344> this will allow an arbitrary number of ligands to be moved separately
!         as rigid units during steptaking moves.
!         Keyword BLOCKMOVE to be used in conjunction with LIGMOVE and MOVABLEATOMS.
!         should be useful for global optimisation of clusters of organic molecules with AMBER.
      BLOCKMOVET=.FALSE. 
      NBLOCKS=1
!
!  csw34> rotamer move stuff
!
      ROTAMERT=.FALSE.
!
!  csw34> some defaults (just in case)
!
      ROTMAXCHANGE=1
      ROTPSELECT=0.2
      ROTOCCUW=0.004
      ROTCENTRE=1
      ROTCUTOFF=999.99
!
!  csw34> atom group rotation moves
!
      GROUPROTT=.FALSE.
      NGROUPS=0
      GROUPROTFREQ=1
      GROUPOFFSET=0
      GROUPROT_SUPPRESS = .FALSE.
! logicals controlling move scaling
      GR_SCALEROT=.FALSE.
      GR_SCALEPROB=.FALSE.
        
!
!  ab2111> dihedral group rotation moves
!
      DIHEDRALROTT=.FALSE.
      BGSMOVE=.FALSE.
      NDIHEDRALGROUPS=0
      DIHEDRALROTFREQ=1
      DIHEDRALOFFSET=0
      PHI0=-0.33 ! * pi radians (rough alpha helix values)
      PSI0=-0.25 ! * pi radians 
      PHIk= 0.1 !  * pi radians  (standard deviation of phi sampling)
      PSIk= 0.1 !  * pi radians
!
!  csw34> HBONDMATRIX
!
      HBONDMATRIX=.FALSE.
      NHBONDGROUPS=0
      MAXHBONDGROUPS=1000
      HBONDTYPE='ACCEPT'
      HBONDACCEPT=.TRUE.
      HBONDLIGAND=.FALSE.
! cutoffs set to ptraj defaults
      HBONDDCUT='3.00'
      HBONDACUT='120.0'
! set defaults for the soft cutoff values
      HBONDDCUTLOOSE='3.06'
      HBONDDCUTTIGHT='2.94'
! note that the angle is (pi-cutoff) so lower is looser
      HBONDACUTLOOSE='118.0'
      HBONDACUTTIGHT='122.0'
! End of HBONDMATRIX defaults
!
! csw34> EXPANDRIGID
!
      EXPANDRIGIDT=.FALSE.
      NORMALISEEXPANDT=.FALSE.
      EXPANDFACTOR=2.0D0
      EXPANDRIGIDFREQ=1
!
! csw34> ROTATERIGID
!
      ROTATERIGIDT=.FALSE.
      ROTATEFACTOR=1.0D0
      ROTATERIGIDFREQ=1
      UPDATERIGIDREFT=.FALSE.
      HYBRIDMINT=.FALSE.
      EPSRIGID=1.0D-3
      COMPRESSRIGIDT=.FALSE.
      KCOMP_RIGID=0.0D0
      RIGIDCOMDIST = HUGE(1.0D0) 
      CHMOVET=.FALSE.
      OMEGAT=.FALSE.
      CHSTANDARDT=.FALSE.
      CHLOOPMODT=.FALSE.
      CHCARTMOVET=.FALSE.
      CHCLUSTERT=.FALSE.
      CHNEIGHBOURT=.FALSE.
      CHBBT=.FALSE.
      CHSCT=.FALSE.
      SECPREDT=.FALSE.
      SECPREDFILE='UNDEFINED'

      OSASAT=.FALSE.
      RPRO=1.4D0
      ODIHET=.FALSE.
      ORGYT=.FALSE.
      OEINTT=.FALSE.
      MON1(1:2)=1
      MON2(1:2)=1

      BSMIN=.FALSE.
      RKMIN=.FALSE.
      PERMDIST=.FALSE.
      PERMOPT=.FALSE.
      DISTOPT=.FALSE.
      PERMINVOPT=.FALSE.

      GAMMA=1.0D0
      TUNNELT=.FALSE.
      
      TWOD=.FALSE.
      COMPRESST=.FALSE.

      MUPDATE=4
      DGUESS=0.1D0
      INTDGUESS=0.001D0
      BFGS=.FALSE.
      LBFGST=.TRUE.
      CONJG=.FALSE.
      TNT=.FALSE.
      TOLB=0.1D0
      DBRENTT=.FALSE.
      GUIDECUT=0.0001D0
      CPMD=.FALSE.
      DL_POLY=.FALSE.
      EFAC=0.0D0
      EAMP=0.01D0
      FIXD=.FALSE.
      NHSMOVE=1
      T12FAC=1.1D0
      RENORM=.FALSE.
      NRENORM=10
      NRENSTUCK=20
      XMOVERENORM=6.0
      TRENORM=1.0D0
      PACHECO=.FALSE.
      EAMLJT=.FALSE.
      PBGLUET=.FALSE.
      EAMALT=.FALSE.
      ALGLUET=.FALSE.
      MGGLUET=.FALSE.
      GUPTAT=.FALSE.
      FST=.FALSE.
      WENZEL=.FALSE.
      RESTART=.FALSE.
      NEWRESTART=.FALSE.
      NRELAX=0
      NMSBSAVE=0
      AVOID=.FALSE.
      AVOIDDIST=1.0D0
      AVOIDRESEEDT=.TRUE.
      GEOMDIFFTOL=0.5D0 !jdf43>
      MAXSAVE=10
      NHSRESTART=0
      MAXBFGS=0.4D0

      CAPSID=.FALSE.
      STRANDT=.FALSE.
      PAHT=.FALSE.
      TIP=.FALSE.
      TTM3T=.FALSE.
      QUADT=.FALSE.
      STOCKT=.FALSE.
      LJCOULT=.FALSE.
      COULN=0
      COULQ=0.0D0
      COULSWAP = 0.0D0
      COULTEMP = 0.0D0
      GAYBERNET=.FALSE.
      ELLIPSOIDT=.FALSE.
      PYGPERIODICT=.FALSE.
      LJCAPSIDT=.FALSE.
      PYBINARYT=.FALSE.
      pyt = .false.
      PYOVERLAPTHRESH=1.0D0
      PYCFTHRESH=0.1D0 ! cold fusion threshold
      LJSITE=.FALSE.
      SWAPMOVEST=.FALSE.
      STICKYT=.FALSE.
      RIGID=.FALSE.
      TIPID=4
      HEIGHT=0.5D0
      OTPT=.FALSE.
      LJMFT=.FALSE.
      CALCQT=.FALSE.
      MONITORT=.FALSE.
      MONITORINT=1
      LOWESTE=1.0D100

! beh35> keyword for mlowest

      MLOWEST = .FALSE.

!     DC430 >
      ALLOCATE(BESTPERM(NATOMSALLOC))! jdf43>
      CHIROT      = .FALSE.
      DBPT        = .FALSE.
      DBPTDT      = .FALSE.
      DMBLPYT     = .FALSE.
      DMBLMT      = .FALSE.
      EFIELDT     = .FALSE.
      GAYBERNEDCT = .FALSE.
      GBDT        = .FALSE.
      GBDPT       = .FALSE.
      GEMT        = .FALSE.
      LINRODT     = .FALSE.
      LWOTPT      = .FALSE.
      MMRSDPT     = .FALSE.
      MORSEDPT    = .FALSE.
      MSGBT       = .FALSE. 
      MSTBINT     = .FALSE.
      MSSTOCKT    = .FALSE.
      MULTPAHAT   = .FALSE.
      NCAPT       = .FALSE.
      NPAHT       = .FALSE.
      NTIPT       = .FALSE.
      NZERO=0                   ! jdf43>
      PAHAT       = .FALSE.
      PAPT        = .FALSE.
      PTSTSTT     = .FALSE.
      RBSYMT      = .FALSE.     ! jdf43>
      SHIFTV=1.0D6              ! jdf43>
      SANDBOXT    = .FALSE.
      SILANET     = .FALSE.
      STOCKAAT    = .FALSE.
      TDHDT       = .FALSE.
      WATERDCT    = .FALSE.
      WATERKZT    = .FALSE.
!|gd351>
      PATCHY = .FALSE.
      ASAOOS = .FALSE.
!<gd351|


       
      THRESHOLDT=.FALSE.
      BSWL=.FALSE.
      BSPT=.FALSE.
      MINDENSITYT=.FALSE.
      BSPTQMAX=1.0D100
      BSPTQMIN=-1.0D100
      BSPTRESTART=.FALSE.
      HISTSMOOTH=.FALSE.
      NSpline=1
      EPSSPHERE=0.0D0
      FIXBIN=.FALSE.
      QUENCHFRQ=1
      NQUENCH=0

      DECAY=.FALSE.
      DECAYPARAM=0.0D0
      COOP=.FALSE.
      NCOOP=5
      COOPCUT=1.0D0

      UNSTRING='UNDEFINED'
      BOXSIZE=20.D0
      SPHERERAD=20.D0
      NCHENCALLS=0
      NATBT=.FALSE.
      MAXERISE=1.0D-10
      MAXERISE_SET=.FALSE.
      MAXEFALL=-HUGE(1.0D0)
      PRINT_PTGRP=.FALSE.
      SYMMETRIZE=.FALSE.
      SYMMETRIZECSM=.FALSE.
      NSYMINTERVAL=10
      SYMTOL1=0.1D0
      SYMTOL2=0.1D0
      SYMTOL3=0.1D0
      SYMTOL4=0.1D0
      SYMTOL5=0.1D0
      PGSYMTOLS(1) = 1.0D-3     ! (Normalised) eigenvalue tolerance
      PGSYMTOLS(2) = 2.0D-1     ! Overlap (distance) tolerance
      PGSYMTOLS(3) = 1.0D-1     ! Rotation matrix tolerance
      NSYMQMAX=20
      MATDIFF=0.1D0
      DISTFAC=0.0D0
      ARMA=0.4D0
      ARMB=0.4D0

      BINSTRUCTURES=.FALSE.
      TETHER=.FALSE.
      EQUIL=0
      PTMC=.FALSE.
      PTMCDUMPSTRUCT=.FALSE.
      PTEX_INDEP = .FALSE.
      VISITPROP=.FALSE.
      HWINDOWS=1

      FixedEndMoveT=.FALSE.
      PIVOTP=0.0D0
      SIDECHAINP=0.0D0

      DIFFRACTT=.FALSE.
      THOMSONT=.FALSE.
      GAUSST=.FALSE.
      MPIT=.FALSE.
      DUMPUNIQUE=.FALSE.
      DUMPUNIQUEEPREV=0.0D0
      DUMPUNIQUEEMARKOV=0.0D0
      DUMPINT=1000 ! default is to dump a restart file every 1000 cycles of mc.f
      RESTORET=.FALSE.
      DUMPFILE=''
      INTEDUMPFILE=''
      MOVESHELLT=.FALSE.
      SHELLMOVEMAX=0
      SHELLPROB=0.0D0
      COREFRAC=0.0D0
      MYSDT=.FALSE.
      TSTAR=-1.0D0
      PRTFRQ=1
      BSPTDUMPFRQ=100000
      QUENCHDOS=.FALSE.
      QDT=.FALSE.
      QD2T=.FALSE.
      MULLERBROWNT=.FALSE.
      QDLIMIT=-1
      CHARMMENERGIES=.FALSE.
      FIXDIHEFLAG=.FALSE.
      JMT=.FALSE.
      PROJIT=.FALSE.
      PROJIHT=.FALSE.
      COLDFUSIONLIMIT=-1.0D6
      MODEL1T=.FALSE.
      CHECKDID=1

!fh301>{{{
      CHEMSHIFT=.FALSE.
      CHEMSHIFT2=.FALSE.
!fh301>}}}

      VGW=.FALSE.              ! VGW PARAMETERS
      LJSIGMA=1.0D0
      LJEPSILON=1.0D0
      TAUMAX=5.0D0
      TAUMAXFULL=7.0D0
      CPFACTORSG=3.0D0
      CPFACTORFG=1.0D0
      CPS=1
      CPF=1
      VGWTOL=0.0001D0
      ACKLANDT=.FALSE.
      ACKLANDID=5
      ACK1=.FALSE.
      ACK2=.FALSE.

      STEEREDMINT=.FALSE.
      DF1T=.FALSE.
      PULLT=.FALSE.
      CSMT=.FALSE.
      CSMGUIDET=.FALSE.
      CSMEPS=1.0D-6
      CSMSTEPS=1
      CSMQUENCHES=1
      CSMMAXIT=0
      CHECKMARKOVT=.FALSE.
      PERCOLATET=.FALSE.
      PERCCUT=1.0D100
      GENALT=.FALSE.
!
! Reservoir list parameters for PT / BSPT.
!
      RESERVOIRT=.FALSE.
      USERES=0
      EXEQ=0
      PERTSTEP=0.01D0
      RES_PSWAP = 0.0D0
      BETA_RES = 0.0D0

!--------------------------------!
! hk286 > Generalised rigid body !
!--------------------------------!
      GENRIGIDT = .FALSE.
      ATOMRIGIDCOORDT = .TRUE.
      RIGIDINIT = .FALSE.
      RELAXFQ = .FALSE.
      RELAXRIGIDT = .FALSE.
      NRELAXRIGIDR = 1000000000
      NRELAXRIGIDA = 1000000000

      RIGIDOPTIMROTAT = .FALSE.
      OPTIMROTAVALUES(:) = 0.0D0
      FREEZERIGIDBODYT = .FALSE.      

      AACONVERGENCET = .FALSE.

!--------------------------------!
! hk286 > Generalised Thomson    !
!--------------------------------!
      GTHOMSONT = .FALSE.
      GTHOMPOT = 1

! hk286 > Damped group moves
      DAMPEDGMOVET = .FALSE.
      DMOVEFREQ = 1

! hk286 > YUKAWA DUMBBELLS
      DBYUKAWAT = .FALSE.
      LAMBDAYAA = 1.0D0
      LAMBDAYAB = 1.0D0
      LAMBDAYBB = 1.0D0
      YEPSFAC = 1.0D0

! hk286 > RESTRAINING POTENTIAL
      RESTRAINLT   = .FALSE.
      RESTRAINLK   = 0.0D0

! jdf43> random numbers
      RANSEEDT=.FALSE.
! jdf43> temperature basin-paving
      TBP=.FALSE.
      TBPMIN=0.0D0
      TBPSTEP=0.0D0
      TBPSTEPS=1
      TBPHF=0.0D0
      TBPCF=0.0D0
      TBPCI=0
      TBPBASIN=0
! jdf43> BHPT & BSPT
      PTRANDOM=.FALSE.
      PTINTERVAL=.FALSE.
      PTSINGLE=.FALSE.
      PTSETS=.FALSE.
! jdf43> RATIO
      RATIOT=.FALSE.
      SRATIO=-1
      SUMSTEP=0.D0
      SUMOSTEP=0.D0
      TRATIO=-1
      SUMTEMP=0.D0
! jdf43> DODECAMORSE
      DDMT=.FALSE.
      DDMCUT=100.D0
! jdf43> MW
      MWFILMT=.FALSE.
      MWPOTT=.FALSE.
      SWPOTT=.FALSE.
! jdf43> SUPPRESS
      SUPPRESST=.FALSE.
! jdf43> MFET
      MFETT=.FALSE.
! jdf43> POLIR
      POLIRT=.FALSE.
! jdf43> MBPOL
      MBPOLT=.FALSE.
      MOLECULART=.FALSE.
! jdf43>
      REPMATCHT=.FALSE.
      
      UNIFORMMOVE=.FALSE.
      ORBITTOL=1.0D-3
      NOINVERSION=.FALSE.
!
! General mixed LJ systems
      GLJT=.FALSE.
! ds656> Multicomponent LJ system (different implementation to GLJ!)
      MLJT=.FALSE.

! khs26> Free energy basin-hopping stuff
      FEBHT = .FALSE.
      FETEMP = 0.0D0
! khs26> This requires a minimum separation between the zero and non-zero
! eigenvalues for runs using free energy basin-hopping. This is based on 
! the minimum found in quench zero.
      MIN_ZERO_SEP = 0.0D0
      MAX_ATTEMPTS = 5
! Use sparse hessian methods
      SPARSET = .FALSE.
      ZERO_THRESH = 0.0D0
! khs26> Rotamer move stuff
      ROTAMER_MOVET = .FALSE.

        
      CUDAT=.FALSE.
      CUDAPOT=' '
      CUDATIMET=.FALSE.
      GCBHT=.FALSE.
      SEMIGRAND_MUT=.FALSE.
      USEROT=.FALSE.
      GCMU=0.0D0
      GCNATOMS=1
      GCINT=100
      GCRELAX=10*GCINT
      GCPLUS=0.5D0
      
      CALL FILE_OPEN('data', DATA_UNIT, .FALSE.)
      
!      OPEN (5,FILE='data',STATUS='OLD')

!190   CALL INPUT(END,5)
190   CALL INPUT(END, DATA_UNIT)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF
      IF (END .OR. WORD .EQ. 'STOP') THEN

         IF (NPCOUNT.LT.NPAR) THEN
            DO J1=NPCOUNT+1,NPAR
               STEP(J1)=STEP(1)
               ASTEP(J1)=ASTEP(1)
               OSTEP(J1)=OSTEP(1)
               BLOCK(J1)=BLOCK(1)
            ENDDO
         ENDIF
! th368: 20-10-2009 Read parameter file containing CHARMM DYNAmics 
! parameters if either CHARMM/MD or CHARMM/NEWRESTART_MD was
! requested terminate if file is not found

         IF(CHMDT .OR. ( CHRMMT .AND. NEWRESTART_MD)) THEN

           INQUIRE(FILE='chmd.par',EXIST=YESNO)

           IF (YESNO) THEN
              OPEN(99,file='chmd.par')
              READ(99,'(A)') CHMDPAR
              CLOSE(99)
           ELSE
              WRITE(MYUNIT,*) 'keywords> File chmd.par has to be provided.'
              STOP
           ENDIF
         ENDIF
! end th368: 20-10-2009

        RETURN
      ENDIF

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                  .OR.WORD.EQ.'\\'.OR.WORD.EQ."!"
     &                  .OR.WORD.EQ."#") THEN 
         GOTO 190

!
!  Remaining documented keywords should be in alphabetical order.
!
      ELSE IF (WORD.EQ.'2D') THEN
         TWOD=.TRUE.

      ELSE IF (WORD.EQ.'ACCEPTRATIO') THEN
         IF (NITEMS-1.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of acceptance ratios exceeds NPAR - quit'
            STOP
         ENDIF
         DO J1=1,NITEMS-1
            CALL READF(ACCRAT(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            IF (NPAR.GT.SIZE(ACCRAT)) THEN
               WRITE(MYUNIT,'(A,I10,A,I10)') 'NPAR=',NPAR,' but SIZE(ACCRAT)=',SIZE(ACCRAT)
               WRITE(MYUNIT,'(A,I10,A,I10)') 'Do you need to move the ACCRAT keyword before MPI in data file?'
               STOP
            ENDIF
            DO J1=NITEMS,NPAR
               ACCRAT(J1)=ACCRAT(1)
            ENDDO
         ENDIF
!
!  bs360: ACE is to be used together with CHARMM and the ACE solvent model,
!  it makes sure that the Born radii are regularly updated
!
      ELSE IF (WORD.EQ.'ACE') THEN
          ACESOLV=.TRUE.
          IF (NITEMS.GT.1) CALL READI(ACEUPSTEP)
!
!  Ackland embedded atom metal potentials.
!
      ELSE IF (WORD.EQ.'ACKLAND') THEN
         ACKLANDT=.TRUE.
         CALL READI(ACKLANDID) ! default is 5 = Fe
c
! Specification of the two possible Ackland potentials for iron.
!
      ELSE IF (WORD.EQ.'ACKLAND1') THEN
         ACK1=.TRUE.
      ELSE IF (WORD.EQ.'ACKLAND2') THEN
         ACK2=.TRUE.
      ELSE IF (WORD.EQ.'ALGLUE') THEN
          ALGLUET=.TRUE.
      ELSE IF (WORD.EQ.'CISTRANS') THEN
          NOCISTRANS=.FALSE.
      ELSE IF (WORD.EQ.'NOCISTRANSCHECKS') THEN
          NOCISTRANS=.FALSE.
!     ELSE IF (WORD.EQ.'AMBER') THEN
!        AMBER=.TRUE.
!        CALL APARAMS
!        CALL AREAD
!        NATOMS=ATOMS
!        DO J1=1,NATOMS
!           COORDS(3*(J1-1)+1,1)=x(J1)
!           COORDS(3*(J1-1)+2,1)=y(J1)
!           COORDS(3*(J1-1)+3,1)=z(J1)
!        ENDDO
!        t=0
!        ang=0
!        imp=0
!        count=0
!
! Dump Amber9 energy components at every step
!
      ELSE IF (WORD.EQ.'AMBERENERGIES') THEN
         AMBERENERGIEST=.TRUE.

      ELSE IF (WORD.EQ.'AMCHPMAX') THEN
         CALL READF(AMCHPMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHPMAX=  ',AMCHPMAX

      ELSE IF (WORD.EQ.'AMCHPMIN') THEN
         CALL READF(AMCHPMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHPMIN=  ',AMCHPMIN

      ELSE IF (WORD.EQ.'AMCHNMAX') THEN
         CALL READF(AMCHNMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHNMAX=  ',AMCHNMAX

      ELSE IF (WORD.EQ.'AMCHNMIN') THEN
         CALL READF(AMCHNMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHNMIN=  ',AMCHNMIN

      ELSE IF (WORD.EQ.'AMH') THEN
         AMHT=.TRUE.
         WRITE(MYUNIT,'(A)')'USING AMH ENERGIES FORCES'
         WRITE(MYUNIT,'(A)')'CALCULATE ENERGY AND FORCE TABLES  '
         CALL WALESAMH_INITIAL 
      ELSE IF (WORD.EQ.'NINT_AMH') THEN
         CALL READI(NINT_AMH)
         WRITE(MYUNIT,*)'NINT_AMH' , NINT_AMH
      ELSE IF (WORD.EQ.'HARM_AMH') THEN
         CALL READI(HARM_AMH)
         WRITE(MYUNIT,*)'HARM_AMH' , HARM_AMH
      ELSE IF (WORD.EQ.'ANGSTROM') THEN
         ANGST=.TRUE.
      ELSE IF (WORD.EQ.'ARGON') THEN
         ARGON=.TRUE.
  
      ELSE IF (WORD.EQ.'ARM') THEN
         ARMT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ARMA)
         IF (NITEMS.GT.2) CALL READF(ARMB)

      ELSE IF (WORD.EQ.'ARNO') THEN
         ARNO=.TRUE.
!
!  Specify resetting if the latest structure gets too close to minima saved
!  in MSBCOORDS. Use bipartite matching and closest approach distance AVOIDDIST.
!  Maximum number of saved structures is specified by MAXSAVE.
! 
      ELSE IF (WORD.EQ.'AVOID') THEN
         AVOID=.TRUE.
         IF (NITEMS.GT.1) CALL READF(AVOIDDIST)
         IF (NITEMS.GT.2) CALL READI(MAXSAVE)
         IF (NITEMS.GT.3) CALL READA(WORD2)
         WORD2=TRIM(ADJUSTL(WORD2))
         IF (WORD2(1:1).EQ.'F') AVOIDRESEEDT=.FALSE.
         IF (MAXSAVE.EQ.0) AVOID=.FALSE.
         IF (MAXSAVE.EQ.0) AVOIDRESEEDT=.TRUE.
         IF (MAXSAVE.EQ.0) MAXSAVE=10 ! to avoid division by zero later

         IF (.NOT.ALLOCATED(MSBCOORDS)) THEN
            ALLOCATE(MSBCOORDS(3*NATOMSALLOC,MAXSAVE))
         ELSE
            WRITE(MYUNIT,*) 'reallocating MSBCOORDS'
            DEALLOCATE(MSBCOORDS)
            ALLOCATE(MSBCOORDS(3*NATOMSALLOC,MAXSAVE))
         ENDIF
         IF (.NOT.ALLOCATED(MSBE)) THEN
            ALLOCATE(MSBE(MAXSAVE))
         ELSE
            WRITE(MYUNIT,*) 'reallocating MSBE'
            DEALLOCATE(MSBE)
            ALLOCATE(MSBE(MAXSAVE))
         ENDIF

      ELSE IF (WORD.EQ.'AXTELL') THEN
         AXTELL=.TRUE.
         CALL READF(ZSTAR)
      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
         IF (NITEMS.GT.1) CALL READF(BQMAX)

      ELSE IF (WORD.EQ.'BFGS') THEN
         BFGS=.TRUE.
!
! PT basin-hopping. This keyword is simply used to read in PTTMIN, PTTMAX, and EXCHPROB.
! It is used in conjunction with MPI to decide if this is a BHPT run in mc.F.
!
      ELSE IF (WORD.EQ.'BHPT') THEN
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         IF (EXCHPROB.GE.1.D0) THEN
            EXCHINT=INT(EXCHPROB)
            EXCHPROB=1.D0/EXCHPROB
         ELSEIF (EXCHPROB.GT.0.D0) THEN
            EXCHINT=INT(1.D0/EXCHPROB)
         ELSE
            EXCHINT=10000000
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RANDOM') PTRANDOM=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'INTERVAL') PTINTERVAL=.TRUE.
         ELSE
            PTRANDOM=.TRUE.
         ENDIF
         IF (NITEMS.GT.5) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SINGLE') PTSINGLE=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SETS') PTSETS=.TRUE.
         ELSE
            PTSINGLE=.TRUE.
         ENDIF
         IF (NITEMS.GT.6) THEN
            CALL READI(NDUMMY)
            CALL SDPRND(NDUMMY)
         ENDIF
!         IF (PTSETS) THEN
!            EXCHINT = INT(EXCHINT*(NPAR-1)/2)
!            EXCHPROB= EXCHPROB*2.D0/(1.D0*(NPAR-1))
!            EXCHPROB = 1.D0/EXCHINT
!         ENDIF

      ELSE IF (WORD.EQ.'BINARY') THEN
         BINARY=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
         
!js850> every BINARY_EXAB_FRQ steps, try to exchange an A atom with a B atom
      ELSE IF (WORD.EQ.'BINARY_EXAB') THEN
         BINARY_EXAB=.TRUE.
         CALL READI(BINARY_EXAB_FRQ)
!
! Saves every n'th structure to the file with corresponding bin label
!
      ELSE IF (WORD.EQ.'BINSTRUCTURES') THEN
         BINSTRUCTURES=.TRUE.
         CALL READI(SAVENTH)
         IF (SAVENTH.LT.1) BINSTRUCTURES=.FALSE.
! BLJCLUSTER 
      ELSE IF (WORD.EQ.'BLJCLUSTER') THEN
         BLJCLUSTER=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
         CALL READF(CUTOFF)
! ds656> BLJCLUSTER_NOCUT 
      ELSE IF (WORD.EQ.'BLJCLUSTER_NOCUT') THEN
         BLJCLUSTER_NOCUT=.TRUE.
         IF (NSPECIES(0) /= 2) THEN
            WRITE(MYUNIT,'(A)') 'keywords> Inconsistent species count for BLJCLUSTER_NOCUT!'
            STOP
         ENDIF
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
         NSPECIES(1) = NTYPEA
         NSPECIES(2) = NATOMS - NTYPEA
         IF (ALLOCATED(GLJEPS)) DEALLOCATE(GLJEPS)
         IF (ALLOCATED(GLJSIG)) DEALLOCATE(GLJSIG)
         ALLOCATE(GLJEPS(2,2),GLJSIG(2,2))
         GLJEPS(1,1) = 1.0D0
         GLJEPS(2,2) = EPSBB
         GLJEPS(1,2) = EPSAB
         GLJEPS(2,1) = EPSAB ! impose symmetry
         GLJSIG(1,1) = 1.0D0
         GLJSIG(2,2) = SIGBB
         GLJSIG(1,2) = SIGAB
         GLJSIG(2,1) = SIGAB ! impose symmetry

      ELSE IF (WORD .EQ. 'PHI4MODEL') THEN
         PHI4MODELT = .TRUE.

! ds656> Homotop refinement for binary systems
      ELSE IF (WORD .EQ. 'HOMOREF') THEN
         HOMOREFT = .TRUE.
         IF(NITEMS >= 5) THEN
            CALL READI(HOMOREF_LSMODE)
            CALL READI(HOMOREF_FGMODE)
            CALL READI(HOMOREF_NCYCLES)
            CALL READI(HOMOREF_NFMAX)
            IF(NITEMS > 5) CALL READI(HOMOREF_NSMAX)
            WRITE(MYUNIT,'(A,I2,A,I2,A,I5,A,I2,A,I2)') 
     2           'keywords> HOMOREF with LSMODE=', HOMOREF_LSMODE, 
     3           ', FGMODE=', HOMOREF_FGMODE, 
     4           ', NCYCLES=', HOMOREF_NCYCLES,
     5           ', NFMAX=', HOMOREF_NFMAX,
     6           ' and NSMAX=', HOMOREF_NSMAX
         ELSE
            WRITE(MYUNIT,'(A)') 
     2           "keywords> HOMOREF requires at least 5 arguments!"
            STOP
         ENDIF
      ELSE IF (WORD .EQ. 'HOMOREF_AUX') THEN
         HOMOREF_AUXT = .TRUE.
         IF(NITEMS .EQ. 5) THEN
            CALL READI(HOMOREF_AUX_NSWAPS)
            CALL READF(HOMOREF_AUX_TEMP)
            CALL READF(HOMOREF_AUX_FACTOR)
            CALL READF(HOMOREF_AUX_NNCUT)
            WRITE(MYUNIT,'(A,I4,A,F6.3,A,F6.3,A,F8.4)') 
     2           'keywords> HOMOREF_AUX with NSWAPS=',
     3           HOMOREF_AUX_NSWAPS,', TEMP=', 
     4           HOMOREF_AUX_TEMP,', FACTOR=',
     5           HOMOREF_AUX_FACTOR,' and NNCUT=',
     6           HOMOREF_AUX_NNCUT
         ELSE
            WRITE(MYUNIT,'(A)') "keywords> HOMOREF_AUX takes 4 arguments!"
            STOP            
         ENDIF
      ELSE IF (WORD .EQ. 'HOMOREF_BH') THEN
         HOMOREF_BHT = .TRUE.
         IF(NITEMS .EQ. 5) THEN
            CALL READI(HOMOREF_BH_NSWAPMAX)
            CALL READI(HOMOREF_BH_NDUDMAX)
            CALL READF(HOMOREF_BH_TEMP)
            CALL READF(HOMOREF_BH_FACTOR)
         ELSE
            WRITE(MYUNIT,'(A)') "keywords> HOMOREF_BH takes 4 arguments!"
            STOP            
         ENDIF
      ELSE IF (WORD .EQ. 'HOMOREFTEST') THEN
         HOMOREFTEST = .TRUE.
      ELSE IF (WORD .EQ. 'HOMOREFCHECK') THEN
         HOMOREFCHECK = .TRUE.
! ds656> Quanech-Assisted Local Combinatorial Search
      ELSE IF (WORD .EQ. 'QALCS') THEN
         QALCST = .TRUE.
         IF(NITEMS > 1) THEN
            CALL READI(QALCSMODE)
            IF(NITEMS > 2) THEN
               CALL READI(QALCS_PARAM)
            ENDIF
         ELSE
            WRITE(MYUNIT,'(A)') 'keywords> QALCS mode unspecified!'
            STOP
         ENDIF
      ELSE IF (WORD .EQ. 'QALCS_SURF') THEN
         QALCS_SURFT = .TRUE.
      ELSE IF (WORD .EQ. 'QALCS_SYM') THEN
         QALCS_SYMT = .TRUE.
         IF(NITEMS > 1) THEN
            CALL READI(QALCS_SYM_MINCORESIZE)
         ENDIF
      ELSE IF (WORD .EQ. 'QALCSVERBOSE') THEN
         QALCSV = .TRUE.
      ELSE IF(WORD .EQ. 'SPECLABELS') THEN
         SPECLABELST = .TRUE.
         ALLOCATE(SPECLABELS(NITEMS-1))
         DO J1=1,NITEMS-1
            CALL READA(SPECLABELS(J1))
         ENDDO
      ELSE IF(WORD .EQ. 'PGSYMTOLS') THEN
         IF(NITEMS .EQ. 4) THEN
            CALL READF(PGSYMTOLS(1))
            CALL READF(PGSYMTOLS(2))
            CALL READF(PGSYMTOLS(3))
         ELSE
            WRITE(MYUNIT,'(A)') 
     2           'keywords> PGSYMTOLS expects 4 arguments.'
            STOP
         ENDIF
      ELSE IF (WORD .EQ. 'KEEPLISTS') THEN
         KEEPLISTS = .TRUE.
         IF(NITEMS .EQ. 3) THEN
            CALL READF(NBRCUT1)
            CALL READF(NBRCUT2)
         ELSE
            WRITE(MYUNIT,'(A)') "keywords> KEEPLISTS requires 2 arguments!"
            STOP
         ENDIF
!  
      ELSE IF (WORD.EQ.'BGUPTANAME') THEN
         IF (NITEMS.GT.1) THEN
             CALL READA(BGUPTANAME1)
          ENDIF
          IF (NITEMS.GT.2) THEN
             CALL READA(BGUPTANAME2)
          ENDIF
       ELSE IF (WORD.EQ.'BGUPTAT') THEN
         BGUPTAT=.TRUE.
         IF (NSPECIES(0) /= 2) THEN
            WRITE(MYUNIT,'(A)') 'keywords> Inconsistent species count for BGUPTA!'
            STOP
         ENDIF
         CALL READI(NTYPEA)
         CALL READF(AAA)
         CALL READF(PAA)
         CALL READF(QAA)
         CALL READF(ZAA)
         CALL READF(R0AA)
         NSPECIES(1) = NTYPEA
         NSPECIES(2) = NATOMS - NTYPEA
       ELSE IF (WORD.EQ.'BGUPTATAB') THEN
         CALL READF(AAB)
         CALL READF(PAB)
         CALL READF(QAB)
         CALL READF(ZAB)
         CALL READF(R0AB)
       ELSE IF (WORD.EQ.'BGUPTATBB') THEN
         CALL READF(ABB)
         CALL READF(PBB)
         CALL READF(QBB)
         CALL READF(ZBB)
         CALL READF(R0BB)
       ELSE IF (WORD.EQ.'BASWAPTEST') THEN
         BASWAPTEST=.TRUE.
       ELSE IF (WORD.EQ.'LSWAPS') THEN
         LSWAPST=.TRUE.
         CALL READI(LSWAPS_N1)
         CALL READI(LSWAPS_N2)
         CALL READF(LSWAPS_TEMP)
         IF(NITEMS > 4) CALL READF(LSWAPS_TFACTOR)
         IF(NITEMS > 5) CALL READI(LSWAPS_NUP)
      ELSE IF (WORD.EQ.'LFLIPS') THEN
         LFLIPST=.TRUE.
         CALL READI(LFLIPS_N1)
         CALL READI(LFLIPS_N2)
         CALL READF(LFLIPS_TEMP)
         IF(NITEMS > 4) CALL READF(LFLIPS_TFACTOR)
         IF(NITEMS > 5) CALL READI(LFLIPS_NUP)
      ELSE IF (WORD.EQ.'LFLIPS_RESET') THEN
         LFLIPS_RESET=.TRUE.
      ELSE IF (WORD.EQ.'BASWAP') THEN
         BASWAP=.TRUE.
         CALL READI(BASWAP_NWAIT)
         !CALL READF(BASWAP_TEMP) ! no longer used!
         CALL READF(BASWAP_FRAC)
         CALL READI(BASWAP_NMAX)
      ELSE IF (WORD.EQ.'RANDPERM') THEN
         RANDPERMT=.TRUE.
         IF(NITEMS > 1) THEN
            CALL READI(RANDPERM_STEP)
         ENDIF
      ELSE IF (WORD.EQ.'RANDMULTIPERM') THEN
         RANDMULTIPERMT=.TRUE.
         IF(NITEMS > 1) THEN
            CALL READI(RANDMULTIPERM_STEP)
         ENDIF
      ELSE IF (WORD.EQ.'MULTIPERM') THEN
         MULTIPERMT=.TRUE.
         IF(NITEMS .GT. 1) THEN
            SPANSWAPST=.TRUE.
         ENDIF
      ELSE IF (WORD.EQ.'ENPERMS') THEN
         ENPERMST=.TRUE.
! ds656> Generalised LJ with Yukawa
      ELSE IF (WORD .EQ. 'GLJY') THEN
         GLJY = .TRUE.
         CALL READI(GLJ_EXP)
         CALL READF(YUK_A)
         CALL READF(YUK_XI)
! ds656> Stress tensor calculation
      ELSE IF (WORD .EQ. 'STRESS') THEN
         STRESST=.TRUE.
         IF(NITEMS.GT.1) CALL READI(STRESS_MODE)
! BLN 
      ELSE IF (WORD.EQ.'BLN') THEN
         BLNT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         READ(1,*) DUMMYCH
         READ(1,*) LJREPBB, LJATTBB
         READ(1,*) LJREPLL, LJATTLL
         READ(1,*) LJREPNN, LJATTNN
         READ(1,*) DUMMYCH
         READ(1,*) DUMMYCH
         READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(1)
         WRITE(MYUNIT,'(A,I8,A)') 'BLN sequence of ',NATOMS,' beads read:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,I8,A)') 'BLN dihedral types:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         WRITE(MYUNIT,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         WRITE(MYUNIT,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         WRITE(MYUNIT,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
!        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
!    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS) 
! End BLN 
! BLN-Go Model 
      ELSE IF (WORD.EQ.'BLNGO') THEN
         GOTYPE=.TRUE.
         BLNT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         IF (NITEMS.GT.3) THEN
            CALL READF(GOFACTOR)
         ENDIF
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         LJREP_BLN=0
         LJATT_BLN=0
         OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         READ(1,*) DUMMYCH
         READ(1,*) LJREPBB, LJATTBB
         READ(1,*) LJREPLL, LJATTLL
         READ(1,*) LJREPNN, LJATTNN
         READ(1,*) DUMMYCH
         READ(1,*) DUMMYCH
         READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(1)
         WRITE(MYUNIT,'(A,I8,A)') 'BLN sequence of ',NATOMS,' beads read:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,I8,A)') 'BLN dihedral types:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         WRITE(MYUNIT,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         WRITE(MYUNIT,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         WRITE(MYUNIT,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
!        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
!    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS)
! End BLN 
!
      ELSE IF (WORD.EQ.'CHAPERONIN') THEN
         CHAPERONINT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         CALL READF(RADIUS_CONTAINER)
         CALL READF(HYDROPHOBIC)
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     $        LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS)
     $        ,A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS)
     $        ,HYDRO_BLN(NATOMS))
         OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         READ(1,*) DUMMYCH
         READ(1,*) LJREPBB, LJATTBB
         READ(1,*) LJREPLL, LJATTLL
         READ(1,*) LJREPNN, LJATTNN
         READ(1,*) LJREPBL, LJATTBL
         READ(1,*) LJREPBN, LJATTBN
         READ(1,*) LJREPLN, LJATTLN
         READ(1,*) DUMMYCH
         READ(1,*) DUMMYCH
         READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(1)
         WRITE(MYUNIT,'(A,I8,A)') 'BLN/CHAPERONIN sequence of ',NATOMS,' beads read:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,I8,A)') 'BLN/CHAPERONIN dihedral types:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,F15.5)') 'Container radius:',RADIUS_CONTAINER
         WRITE(MYUNIT,'(A,F15.5)') 'Hydrophobicity coefficient:',HYDROPHOBIC
         WRITE(MYUNIT,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         WRITE(MYUNIT,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         WRITE(MYUNIT,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         WRITE(MYUNIT,'(A,2F15.5)') 'B-L LJ coefficients: ',LJREPBL, LJATTBN
         WRITE(MYUNIT,'(A,2F15.5)') 'B-N LJ coefficients: ',LJREPBN, LJATTBN
         WRITE(MYUNIT,'(A,2F15.5)') 'L-N LJ coefficients: ',LJREPLN, LJATTLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayCHAPERONIN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN
     $        ,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,LJREPBB, LJATTBB,
     $        LJREPLL, LJATTLL, LJREPNN, LJATTNN,LJREPBL, LJATTBL,
     $        LJREPBN, LJATTBN, LJREPLN, LJATTLN,HABLN, HBBLN, HCBLN,
     $        HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN,
     $        TDBLN, HYDROPHOBIC, HYDRO_BLN, NATOMS)
!        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
!    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS) 
! End CHAPERONIN 
      ELSE IF (WORD.EQ.'BSMIN') THEN
         BSMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
         WRITE(MYUNIT,'(A,2L5)') 'BSMIN branch RKMIN,BSMIN=',RKMIN,BSMIN
!
! Basin-sampling. 
!
      ELSE IF (WORD.EQ.'BSPT') THEN
         BSPT=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(PTEMIN)
         CALL READF(PTEMAX)
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         CALL READF(NEQUIL)  ! equilibration only
         CALL READF(PTSTEPS) ! PT following equilibration
         CALL READF(NQUENCH) ! combined PT and quench steps
         CALL READI(NENRPER)
         CALL READI(HBINS)
         CALL READI(QUENCHFRQ)
         write (*,*) "js850> quenchfrq", quenchfrq
         call flush(6)
         IF (NITEMS.GT.14) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'R') PTRANDOM=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'I') PTINTERVAL=.TRUE.
         ELSE
            PTRANDOM=.TRUE.
         ENDIF
         IF (NITEMS.GT.15) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SIN') PTSINGLE=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SET') PTSETS=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'IND') PTEX_INDEP=.TRUE.
            ! ** Note: if use use PTEX_INDEP then you must use PTINTERVAL.  No
            ! physical reason for this, just a result of the hacky way
            ! tryexchange is currently coded ***
         ELSE
            PTSINGLE=.TRUE.
         ENDIF
         EXCHINT=INT(1.D0/EXCHPROB)
! jdf43> the following would enforce the same average number of pair
! exchanges for the single & sets exchange schemes - ss2029 prefers this
! not to be the case.
!
!         IF (PTSETS) THEN
!            EXCHINT = INT(EXCHINT*(NPAR-1)/2)
!            EXCHPROB= EXCHPROB*2.D0/(1.D0*(NPAR-1))
!         ENDIF
! 
!
! Frequency of dumping the Visits.his.n files
! 
      ELSE IF (WORD.EQ.'BSPTDUMPFRQ') THEN
         CALL READI(BSPTDUMPFRQ)
!
! BSPT restriction on quench energy.
!
      ELSE IF (WORD.EQ.'BSPTQRANGE') THEN
         CALL READF(BSPTQMIN)
         CALL READF(BSPTQMAX)
!
! Restart from Visits and bsptrestart files.
!
      ELSE IF (WORD.EQ.'BSPTRESTART') THEN
         BSPTRESTART=.TRUE.
!
! WL Basin-sampling. 
!
      ELSE IF (WORD.EQ.'BSWL') THEN
         BSWL=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(HISTFAC)
         CALL READI(HBINS)
         CALL READF(HISTFACMUL)
         CALL READI(TargetWL)
         CALL READF(HPERCENT)
         ALLOCATE(HDIST(HBINS),HWEIGHT(HBINS),HISTVALS(HBINS),LHISTVALS(HBINS),IGNOREBIN(HBINS))
         DO J1=1,HBINS
            HISTVALS(J1)=0
            LHISTVALS(J1)=0
            HWEIGHT(J1)=1.0D0
            HDIST(J1)=0.0D0
!           DO J2=1,HBINS
!              HTRANS(J1,J2)=1.0D0 ! transition matrix
!           ENDDO
         ENDDO
!
!  During the BSWL run HDIST contains the sum of the distances found for minima in each bin. The
!  average is saved in hist.new.
!
         DO J1=1,HBINS
            HDIST(J1)=HDIST(J1)*HISTVALS(J1)
         ENDDO
!

!fh301>{{{
!
! CAMSHIFT calculates the NMR chemical shifts from atomic coordinates and is mostly used for
! chemical shift restrained simulations
!
      ELSE IF (WORD .EQ. 'CAMSHIFT') THEN
        CHEMSHIFT=.TRUE.
        CHEMSHIFT2=.TRUE.
!        CHEMSHIFTITER=0
        IF (NITEMS .LT. 4) THEN
          WRITE(*,*) 'CamShift version, path and shiftfile are needed for Camshift'
          STOP
        ELSE
          CALL READA(CSVERSION)
          CSVERSION=TRIM(ADJUSTL(CSVERSION))
          SELECT CASE(CSVERSION)
            CASE('MERGE','ORIGINAL','NOFF')
            CASE DEFAULT
              WRITE(*,*) TRIM(ADJUSTL(CSVERSION)),' is not a correct CamShift version'
              STOP
          END SELECT
          CALL READA(SVNROOT)
          SVNROOT=TRIM(ADJUSTL(SVNROOT))
          CSPATH=TRIM(ADJUSTL(SVNROOT))//'CAMSHIFTDATA/'
          CSPATH=TRIM(ADJUSTL(CSPATH))
          CALL READA(SHIFTFILE)
          SHIFTFILE=TRIM(ADJUSTL(SHIFTFILE))
          IF (NITEMS .GT. 4) THEN
            CALL READF(CSN)
          ENDIF
          IF (NITEMS .GT. 5) THEN
            CALL READF(CSALPHA)
          ENDIF
        ENDIF
!fh301>}}}
 
      ELSE IF (WORD.EQ.'CAPSID') THEN
         CAPSID=.TRUE.
         RIGID=.TRUE.
         CALL READF(RHO)
         CALL READF(EPS2)
         CALL READF(RAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)
!
!  The six reference site positions per capped pentagon. These need to
!  be multiplied by RAD, including the repulsive site!
!
         NRBSITES=6
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=1.0D0
         SITE(1,2)=0.0D0
         SITE(1,3)=0.0D0
         SITE(2,1)=((-1.0D0 + Sqrt(5.0D0)))/4.0D0
         SITE(2,2)=(Sqrt((5.0D0 + Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(2,3)=0.0D0
         SITE(3,1)=((-1.0D0 - Sqrt(5.0D0)))/4.0D0
         SITE(3,2)=(Sqrt((5.0D0 - Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(3,3)=0.0D0
         SITE(4,1)=((-1 - Sqrt(5.0D0)))/4.0D0
         SITE(4,2)=-(Sqrt((5.0D0 - Sqrt(5.0D0))/2.))/2.0D0
         SITE(4,3)=0.0D0
         SITE(5,1)=((-1 + Sqrt(5.0D0)))/4.0D0
         SITE(5,2)=-(Sqrt((5.0D0 + Sqrt(5.0D0))/2.))/2.0D0
         SITE(5,3)=0.0D0
         SITE(6,1)=0.0D0
         SITE(6,2)=0.0D0
         SITE(6,3)=HEIGHT
      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.
       
! csw34> When using the implicit membrane potential IMM1, we want to let
! the molecule move in and out of the membrane (z direction), but fix it
! in x and y.
       
      ELSE IF (WORD.EQ.'CENTREXY') THEN
         CENT=.TRUE.
         CENTXY=.TRUE.

! csw34> SETCENTRE moves the centre of coordinates to the specified
! location before the initial quench is done.

      ELSE IF (WORD.EQ.'SETCENTRE') THEN
         SETCENT=.TRUE.
         IF (NITEMS.EQ.2) THEN
            CALL READF(CENTX) 
         ELSE IF (NITEMS.EQ.3) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
         ELSE IF (NITEMS.EQ.4) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
            CALL READF(CENTZ)
         ENDIF
      ELSE IF (WORD.EQ.'QCIPOT') THEN
         QCIPOTT=.TRUE.
      ELSE IF (WORD.EQ.'QCIPOT2') THEN
         QCIPOTT=.TRUE.
         QCIPOT2T=.TRUE.

!       csw34> QUCENTRE moves the centre of coordinates to (0,0,0)
!       before each step is taken. 
      ELSE IF (WORD.EQ.'QUCENTRE') THEN
         QUCENTRE=.TRUE.
!
!  Conjugate gradient optimisation instead of default LBFGS
!
      ELSE IF (WORD.EQ.'CG') THEN
         LBFGST=.FALSE.
         CONJG=.TRUE.
! 
! sf344> Start of AMBER-related keywords 
!
      ELSE IF (WORD.EQ.'AMBERMDSTEPS') THEN
        MDSTEPT = .TRUE.
      ELSE IF (WORD .EQ. 'AMBER12') THEN
        AMBER12T = .TRUE.
! Read the coords from AMBER12 into COORDS(:,1)
        CALL AMBER12_GET_COORDS(NATOMS, COORDS(:,1))
      ELSE IF (WORD.EQ.'AMBER9') THEN
        AMBERT=.TRUE.
        WRITE(MYUNIT,'(A)') 'keyword> RADIUS set to 999 for AMBER9 run'
        RADIUS=999
        
!
! csw34> if residues are frozen with FREEZERES, call the amber routine
! to fill the FROZEN array correctly (in amberinterface.f) 
!
        IF (PERMDIST) THEN
          IF(NPERMSIZE(1).EQ.NATOMS) THEN
           PRINT '(A)','keyword> ERROR - PERMDIST is specified for AMBER, but there is no perm.allow file present'
           STOP
          ENDIF
        ENDIF

! sf344> file open unit used to conflict with AMBER's IO units (mdin opened with unit = 5),

!               call amberinterface(natom,1,trim(adjustl(inpcrd)),MYUNIT)

        IF(NITEMS==2) then
         inpcrd='coords.inpcrd'
         CALL READA(amberstr)
               call amberinterface(natom,2,inpcrd,MYUNIT)
         IF(MPIT) THEN
             WRITE(J1CHAR,'(I3)') MYNODE + 1
             WRITE(J2CHAR,'(A)') trim(adjustl(trim(adjustl(amberstr))//'.'//trim(adjustl(J1CHAR))))
             WRITE(MYUNIT,*) 'keywords> input coordinates for AMBER9 system will be read from '// trim(adjustl(J2CHAR))
             CALL FLUSH(MYUNIT)
             call amber_readcoords(J2CHAR)
         ELSE

             WRITE(MYUNIT,'(A)') 'keywords> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr))
             CALL amber_readcoords(amberstr)
         END IF
        ELSE IF(NITEMS==3) then
         CALL READA(amberstr)
         CALL READA(amberstr1)
         WRITE(MYUNIT,'(A)') 'keywords> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr)),
     &                              'type: ', trim(adjustl(amberstr1))
          IF(trim(adjustl(amberstr1)).EQ.'inpcrd') then
               inpcrd=amberstr
               call amberinterface(natom,2,inpcrd,MYUNIT)
           WRITE(MYUNIT,'(A)') 'keywords> reading AMBER inpcrd coordinate format'
          ELSE
           WRITE(MYUNIT,'(A)') 'keywords> ERROR - no other types defined currently than inpcrd'
           STOP
          END IF
        END IF
               IF(.NOT.ALLOCATED(COORDS1)) ALLOCATE(COORDS1(3*NATOM))
               IF(.NOT.ALLOCATED(MOVABLEATOMLIST)) ALLOCATE(MOVABLEATOMLIST(NATOMS))
               IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
               ALLOCATE(COORDS(3*NATOM,NPAR))
               NATOMS = NATOM
             DO J1=1,NPAR
               COORDS(:,J1) = COORDS1(:)
             END DO
!                natoms = natom

!                WRITE(MYUNIT,*) 'sf344> keywords.f, natoms = ', natoms
!        call amopen(9,AMBERPRMTOP,'O','F','R')
        IF (FREEZERES) THEN 
           IF (UNFREEZERES) THEN
              CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE,.TRUE.)
           ELSE
              CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE,.FALSE.)
           ENDIF
        ENDIF

        IF(DONTMOVEREST) THEN
           IF (DOMOVEREST) THEN
              CALL A9RESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.TRUE.)
           ELSE
              CALL A9RESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.FALSE.)
           ENDIF
        ENDIF

! khs26> Initialise MME2 calls for calculating analytical second derivatives with NAB.
        IF (IGB.EQ.6) THEN
            ! khs26> NAB doesn't recognise igb=6
            CALL MMEINITWRAPPER(TRIM(ADJUSTL(PRMTOP))//C_NULL_CHAR, 0, SALTCON, RGBMAX, CUT)
        ELSE
            CALL MMEINITWRAPPER(TRIM(ADJUSTL(PRMTOP))//C_NULL_CHAR, IGB, SALTCON, RGBMAX, CUT)
        END IF

!
! The maximum number of constraints to use in the constrained potential.
! The deafult is 4.
!
      ELSE IF (WORD.EQ.'MAXCON') THEN
         CALL READI(MAXCONUSE)

      ELSE IF(WORD.EQ.'MODEL1') THEN
         MODEL1T=.TRUE.
         CALL READF(ME1)
         CALL READF(ME2)
         CALL READF(ME3)
         CALL READF(MSTART)
         CALL READF(MFINISH)
         CALL READF(MBSTART1)
         CALL READF(MBFINISH1)
         CALL READF(MBSTART2)
         CALL READF(MBFINISH2)
         CALL READF(MBHEIGHT1)
         CALL READF(MBHEIGHT2)
      ELSE IF(WORD.EQ.'MONITOR') THEN
         MONITORT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(MONITORINT)

! sf344> keyword for taking constrained random steps of just one part of the molecule between quenches
!        should be useful for systems with ligands docked in a protein, this would be correspond to some
!        partially flexible docking scheme
      ELSE IF(WORD.EQ.'MOVABLEATOMS') THEN
         MOVABLEATOMST=.TRUE.
         nmovableatoms=0
! csw34> for steered minimisation, need to read in the atoms used to define the points force acts on/from
!        group A is the ligand, group B is the protein. They start at -1 as groups delimited by A and B
! csw34> also, for LOCALSAMPLEing, need to define a group of atoms in the ligand (A) and two groups of atoms (B and C) 
!        that are used to calculate distances to the ligand
         natomsina=-1
         natomsinb=-1
         natomsinc=-1
!        readswitch controls which part of movableatoms we're reading in i.e. 
!        M=whole ligand, A=group A, B=group B, C=group C
         readswitch='M'
           WRITE(MYUNIT,'(A)') ' keyword> list of movable atoms will be read from file <movableatoms>'
           OPEN(unit=222,file='movableatoms',status='old')
                 do
                   read (unit=222,iostat=iostatus,fmt='(A6)') check1
                     if (iostatus<0) then
                        close(222)
                        exit
                     else if (TRIM(ADJUSTL(check1)).EQ.'A') then
                        readswitch='A'
                     else if (TRIM(ADJUSTL(check1)).EQ.'B') then
                        readswitch='B'
                     else if (TRIM(ADJUSTL(check1)).EQ.'C') then
                        readswitch='C'
                     end if
                     if (readswitch.EQ.'M') then 
                        nmovableatoms = nmovableatoms + 1
                     else if (readswitch.EQ.'A') then
                        natomsina = natomsina + 1
                     else if (readswitch.EQ.'B') then
                        natomsinb = natomsinb + 1
                     else if (readswitch.EQ.'C') then
                        natomsinc = natomsinc + 1
                     endif
                 end do
! setup arrays for movableatoms
           if(.not.allocated(movableatomlist)) ALLOCATE(movableatomlist(nmovableatoms))
           if(.not.allocated(movableatomlistlogical)) ALLOCATE(movableatomlistlogical(natoms))
           movableatomlistlogical(:)=.false.
! setup arrays for steered minimisation groups
           if (natomsina.gt.0) then 
              if(.not.allocated(atomsinalist)) ALLOCATE(atomsinalist(natomsina))
              if(.not.allocated(atomsinalistlogical)) ALLOCATE(atomsinalistlogical(natoms))
              atomsinalistlogical(:)=.false.
           endif
           if (natomsinb.gt.0) then 
              if(.not.allocated(atomsinblist)) ALLOCATE(atomsinblist(natomsinb))
              if(.not.allocated(atomsinblistlogical)) ALLOCATE(atomsinblistlogical(natoms))
              atomsinblistlogical(:)=.false.
           endif
           if (natomsinc.gt.0) then 
              if(.not.allocated(atomsinclist)) ALLOCATE(atomsinclist(natomsinc))
              if(.not.allocated(atomsinclistlogical)) ALLOCATE(atomsinclistlogical(natoms))
              atomsinclistlogical(:)=.false.
           endif
! now open movableatoms for the second time to actually read in the atom indices
! khs26> Also check that the movable atom indices are between 1 and NATOMS, otherwise we end up
! writing past array bounds in amberinterface.f
           OPEN(unit=222,file='movableatoms',status='old')
              do i=1,nmovableatoms
                 read (unit=222,iostat=iostatus,fmt='(I6)') MOVABLEATOMINDEX
                 IF ((MOVABLEATOMINDEX .LE. 0) .OR. (MOVABLEATOMINDEX .GT. NATOMS)) THEN
                    WRITE(MYUNIT, '(A,I6,A,I6,A)') 'The index ', MOVABLEATOMINDEX, ' is out of bounds (NATOMS = ', NATOMS, ').'
                    STOP 'Error: movableatoms contains atom indices which are out of bounds.'
                 END IF
                 movableatomlist(i) = MOVABLEATOMINDEX
              end do
! if groups for steered minimisation are specified, read them in
           if (natomsina.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsina
                    read (unit=222,iostat=iostatus,fmt='(I6)') MOVABLEATOMINDEX
                    IF ((MOVABLEATOMINDEX .LE. 0) .OR. (MOVABLEATOMINDEX .GT. NATOMS)) THEN
                       WRITE(MYUNIT, '(A,I6,A,I6,A)') 'The index ', MOVABLEATOMINDEX, ' is out of bounds (NATOMS = ', NATOMS, ').'
                       STOP 'Error: movableatoms contains atom indices which are out of bounds.'
                    END IF
                    atomsinalist(i) = MOVABLEATOMINDEX
                 end do
           endif
           if (natomsinb.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsinb
                    read (unit=222,iostat=iostatus,fmt='(I6)') MOVABLEATOMINDEX
                    IF ((MOVABLEATOMINDEX .LE. 0) .OR. (MOVABLEATOMINDEX .GT. NATOMS)) THEN
                       WRITE(MYUNIT, '(A,I6,A,I6,A)') 'The index ', MOVABLEATOMINDEX, ' is out of bounds (NATOMS = ', NATOMS, ').'
                       STOP 'Error: movableatoms contains atom indices which are out of bounds.'
                    END IF
                    atomsinblist(i) = MOVABLEATOMINDEX
                 end do
           endif
           if (natomsinc.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsinc
                    read (unit=222,iostat=iostatus,fmt='(I6)') MOVABLEATOMINDEX
                    IF ((MOVABLEATOMINDEX .LE. 0) .OR. (MOVABLEATOMINDEX .GT. NATOMS)) THEN
                       WRITE(MYUNIT, '(A,I6,A,I6,A)') 'The index ', MOVABLEATOMINDEX, ' is out of bounds (NATOMS = ', NATOMS, ').'
                       STOP 'Error: movableatoms contains atom indices which are out of bounds.'
                    END IF
                    atomsinclist(i) = MOVABLEATOMINDEX
                 end do
           endif
! now need to set the logicals to .true. for all atoms in movableatom groups
              do i=1,nmovableatoms
                 movableatomlistlogical(movableatomlist(i))=.true.
              end do
           if (natomsina.gt.0) then
              do i=1,natomsina
              end do
           endif
           if (natomsinb.gt.0) then
              do i=1,natomsinb
                 atomsinblistlogical(atomsinblist(i))=.true.
              end do
           endif
           if (natomsinb.gt.0) then
              do i=1,natomsinc
                 atomsinclistlogical(atomsinclist(i))=.true.
              end do
           endif
! write some output for the user
           WRITE(MYUNIT,'(A15,I6,A46)') ' keyword> ligand defined with ',nmovableatoms,' atoms'
           IF(natomsina.GT.0) WRITE(MYUNIT,'(A10,I6,A64)') ' keyword> ',natomsina,' atoms in group A'
           IF(natomsinb.GT.0) WRITE(MYUNIT,'(A10,I6,A63)') ' keyword> ',natomsinb,' atoms in group B'
           IF(natomsinc.GT.0) WRITE(MYUNIT,'(A10,I6,A63)') ' keyword> ',natomsinc,' atoms in group C'
! sf344> keyword for taking moves between quenches in which one part of the molecule is moved
! csw34> updated on 29th September 2009 to allow local cartesian moves 
      ELSE IF(WORD.EQ.'LIGMOVE') THEN
        LIGMOVET=.TRUE.
        IF (NITEMS.GT.1) THEN
           CALL READF(ligrotscale)
        ENDIF
        IF (NITEMS.GT.2) THEN
           CALL READF(ligcartstep)
        ENDIF
        IF (NITEMS.GT.3) THEN
           CALL READF(ligtransstep)
        ENDIF
        IF (NITEMS.GT.4) THEN
           CALL READI(ligmovefreq)
! csw34> to prevent an arithmetic exception (divide by 0), we need to
! check that the frequency of ligand moves is > 0. Otherwise, disable
! the moves.
           IF(ligmovefreq.EQ.0) THEN 
              LIGMOVET=.FALSE.
              WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of LIGMOVE moves set to 0 - moves DISABLED!'
           ENDIF
        ENDIF

! csw34> some user info about ligand moves 
        IF (ligrotscale.gt.0) THEN
           WRITE(MYUNIT,'(A)') ' keyword> one part of the molecule (ligand) will be randomly rotated during MC steptaking moves'
           WRITE(MYUNIT,'(A,G10.3)') ' keyword> ligand rotations scaled will be scaled by ',ligrotscale 
        ENDIF
        IF (ligcartstep.gt.0) THEN 
           WRITE(MYUNIT,'(A,G10.3,A)') ' keyword> cartesian perturbations of up to ',ligcartstep,' will be applied to the ligand' 
        ENDIF 
! sf344> block moves. To be used for rigid body moves of multiple ligands or small clusters of molecules. Rigid body move 
!        parameters will be set using LIGMOVE. Syntax: BLOCKMOVE nblocks block1 block2 ... blockN. The file movableatomslist
!        has to contain all atom indices, and these will be separated into nblocks parts, each block containing block_i number
!        of atoms.
      ELSE IF(WORD.EQ.'BLOCKMOVE') THEN
        BLOCKMOVET=.TRUE.
        CALL READI(NBLOCKS)
        DEALLOCATE(ATOMSINBLOCK)
        ALLOCATE(ATOMSINBLOCK(NBLOCKS))
        DO J1=1,NBLOCKS
           CALL READI(ATOMSINBLOCK(J1))
        END DO
        WRITE(MYUNIT,*) ' keywords> atoms in movableatomlist will be moved independently in blocks of ', ATOMSINBLOCK(:)
!
! Explicit dump of image coordinates in xyz format for intlbfgs. Should
! be set .TRUE. if DEBUG is set.
!
      ELSE IF (WORD == 'DUMPINTXYZ') THEN
         DUMPINTXYZ=.TRUE.
         IF (NITEMS>1) CALL READI(DUMPINTXYZFREQ)
!
! Explicit dump of interpolation EofS for intlbfgs. Should be set .TRUE. if DEBUG is set.
!
      ELSE IF (WORD == 'DUMPINTEOS') THEN
         DUMPINTEOS=.TRUE.
         IF (NITEMS>1) CALL READI(DUMPINTEOSFREQ)
      ELSE IF (WORD.eq.'DUMPSTRUCTURES') THEN
        DUMPSTRUCTURES=.TRUE.
        WRITE(MYUNIT,'(A)') ' keywords> Final structures will be dumped in different formats (.rst, .xyz, .pdb)'
      ELSE IF (WORD.eq.'DUMPMIN') THEN
        DUMPMINT=.TRUE.
        WRITE(MYUNIT,'(A)') ' keywords> The SAVE lowest minima will be dumped every DUMPINT steps as dumpmin.x files'
      ELSE IF (WORD.eq.'RANDOMSEED') THEN
        WRITE(MYUNIT,'(A)') ' keywords> The random number generator for the random steptaking moves will be seeded with system time'
        RANDOMSEEDT=.TRUE.
! khs26> This writes the energy decomposition for AMBER 
      ELSE IF (WORD.eq.'ENERGY_DECOMP') THEN
        WRITE(MYUNIT,'(A)') ' keywords> After each quench, the energy decomposition will also be printed.'
        ENERGY_DECOMPT=.TRUE.
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(IX)
         NACCEPT=IX

! js850> if DUMPUNIQUE then dump all unique structures to file dumpunique as
! soon as they are found.  The test for uniqueness is by difference in energy
! using the criterion ECONV.  Only the last structure is tested against, so the
! uniqueness test is only approximate
      ELSE IF (WORD.eq.'DUMPUNIQUE') THEN
        DUMPUNIQUE=.TRUE.

! csw34> Dumping after every quench in AMBER rst and pdb format
      ELSE IF (WORD.eq.'DUMPQU') THEN
         DUMPQUT=.TRUE.
! csw34> Dumping after every step (before quenching) in AMBER rst and pdb format
      ELSE IF (WORD.eq.'DUMPSTEPS') THEN
         DUMPSTEPST=.TRUE.
! khs26> Dump best structures in pdb format after every step
      ELSE IF (WORD.eq.'DUMPBEST') THEN
         DUMPBESTT=.TRUE.
! csw34> AMBER interaction energy using script AMBGMINintE.sh
      ELSE IF (WORD.eq.'A9INTE') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='AMBGMINintE.sh',EXIST=YESNO)
         IF (YESNO) THEN
            A9INTET=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> The interaction enthalpy to the specified residue will be calculated after each quench'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: NEED AMBGMINintE.sh SCRIPT TO USE A9INTE - see www-wales.ch.cam.ac.uk/software'
            STOP
         ENDIF

! csw34> FREEZEGROUP 
! FREEZEGROUP lets you freeze all atoms within or beyond a radius
      ELSE IF (WORD.EQ.'FREEZEGROUP') THEN
         FREEZE=.TRUE.
         FREEZEGROUPT=.TRUE.
         CALL READI(GROUPCENTRE)
         CALL READF(GROUPRADIUS)
         IF(NITEMS.GT.3) CALL READA(FREEZEGROUPTYPE)

! csw34> Rotamer moves 
      ELSE IF (WORD.EQ.'ROTAMER') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='PdbRotamerSearch',EXIST=YESNO)
         IF (YESNO) THEN
            ROTAMERT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> AMBER rotamer moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: NEED PdbRotamerSearch TO USE ROTAMER - see www-wales.ch.cam.ac.uk/software'
            STOP
         ENDIF
         CALL READI(ROTMAXCHANGE)
         CALL READF(ROTPSELECT)
         CALL READF(ROTOCCUW)
         IF (NITEMS.GT.4) CALL READI(ROTCENTRE)
         IF (NITEMS.GT.5) CALL READF(ROTCUTOFF)
!
! Start of CHARMM-related keywords, including options for MC moves and order parameter specifications.
!
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHRMMT=.TRUE.
         IF (MXATMS.EQ.0) THEN
            WRITE(MYUNIT,'(A)') 'keyword> ERROR *** MXATMS is zero'
            STOP
         ENDIF
         CALL FLUSH(MYUNIT)

         IF (PERMDIST) THEN
            IF(NPERMSIZE(1).EQ.NATOMS) THEN
            WRITE(MYUNIT,'(A)') 'keyword> ERROR - PERMDIST is specfied for CHARMM, but there is no perm.allow file present'
            STOP
            ENDIF
         ENDIF

         ALLOCATE(CHX(MXATMS),CHY(MXATMS),CHZ(MXATMS),CHMASS(MXATMS))

         CHX(1)=13.13d13 ! this way we will tell CHARMM to save it's coords into CH. arrays; otherwise it will
                         ! use input.crd only which is the default now
         CALL FLUSH(MYUNIT)
         CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE,DATA_UNIT)
         CALL FLUSH(MYUNIT)
         CALL CHSETZSYMATMASS
         CALL CHALLOCATE(NATOMS)
! jmc49>         CALL CHSETDIHE
! Moved the above call to CHSETDIHE to below the "IF (FREEZERES) CALL CHRESTOATOM(FROZENRES,FROZEN)" below.
! Necessary for the case of freezing residues when the MC moves are made in internal coordinates, and 
! inconsequential otherwise.

! csw34> If we're freezing RESIDUES, fill the FROZEN array with atoms that are
! within the frozen residues. This requires info from a CHARMM routine
! and so will be done in charmmgmin.src

         IF (FREEZERES) THEN 
            IF (UNFREEZERES) THEN
               CALL CHRESTOATOM(FROZENRES,FROZEN,NFREEZE,.TRUE.)
            ELSE
               CALL CHRESTOATOM(FROZENRES,FROZEN,NFREEZE,.FALSE.)
            ENDIF
         ENDIF

         IF(DONTMOVEREST) THEN
            IF (DOMOVEREST) THEN
               CALL CHRESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.TRUE.)
            ELSE
               CALL CHRESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.FALSE.)
            ENDIF
         ENDIF
         CALL CHSETDIHE
         
         IF (NATOM /= NATOMS) THEN
            WRITE(MYUNIT,'(A)') 'No. of atoms in "input.crd" and file specified in CHARMM part of data conflict'
            WRITE(MYUNIT,'(A,2I8)') 'NATOM,NATOMS=',NATOM, NATOMS
            CALL EXIT(10)
         ENDIF
         IF (MPIT) THEN
            DO J1=1,NATOMS
               COORDS(3*(J1-1)+1,MYNODE+1)=CHX(J1)
               COORDS(3*(J1-1)+2,MYNODE+1)=CHY(J1)
               COORDS(3*(J1-1)+3,MYNODE+1)=CHZ(J1)
            ENDDO
         ELSE
            DO J1=1,NATOMS
               COORDS(3*(J1-1)+1,1)=CHX(J1)
               COORDS(3*(J1-1)+2,1)=CHY(J1)
               COORDS(3*(J1-1)+3,1)=CHZ(J1)
            ENDDO
         ENDIF
         DEALLOCATE(CHX,CHY,CHZ,CHMASS)
         IF (INTMINT) CALL GETNINT(NINTS)  ! DJW - this is OK because CHARMM is the last keyword!
         IF (ODIHET) CALL READREF ! likewise, this is OK and in fact essential because CHARMM is the last keyword!
!
! Dump charmm energy components at every step
!
      ELSE IF (WORD.EQ.'CHARMMENERGIES') THEN
         CHARMMENERGIES=.TRUE.
! csw34> If using the CHARMM SCC-DFTB potential, we assume that all
! atoms are QM. If you are using a mixed QM/MM system, you should either
! not use the CHARMMDFTB keyword, or re-code it to check for fully QM
! systems. This keyword essentially prevents unnessesary printing!
!
      ELSE IF (WORD.EQ.'CHARMMDFTB') THEN
         CHARMMDFTBT=.TRUE.
         WRITE(MYUNIT,'(A)') 'keyword> WARNING - All atoms assumed to be QM, NBONDS calls disabled'
!
      ELSE IF (WORD.EQ.'CHARMMTYPE') THEN
         IF (NITEMS.GT.1) THEN
            CALL READA(TOPFILE)
            TOPFILE=TRIM(ADJUSTL(TOPFILE))
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READA(PARFILE)
            PARFILE=TRIM(ADJUSTL(PARFILE))
         ELSE
            WRITE(MYUNIT,*) 'keyword> TOPFILE and PARFILE have to be defined for CHARMMTYPE'
            STOP
         ENDIF
         IF (TOPFILE(1:6).EQ."toph19") THEN
            CHARMMTYPE=2
         ELSEIF (TOPFILE(1:7).EQ."top_all") THEN
            CHARMMTYPE = 1
         ELSE
             WRITE(MYUNIT,*) 'keywords> TOPFILE ', TRIM(ADJUSTL(TOPFILE)),' is not recognised by GMIN'
             STOP
         ENDIF
         WRITE(MYUNIT,'(A,I2)') 'CHARMMTYPE set to ',CHARMMTYPE
!
! Sanity check for the energy in COORDSO.
!
      ELSE IF (WORD.EQ.'CHECKMARKOV') THEN
         CHECKMARKOVT=.TRUE.
!
! Parameters for recalculating the repulsions in INTCONSTRAINT
!
         ELSE IF (WORD.EQ.'CHECKREP') THEN
            IF (NITEMS.GT.1) CALL READI(CHECKREPINTERVAL)
            IF (NITEMS.GT.2) CALL READF(CHECKREPCUTOFF)
!
! MD for CHARMM for the generation of new structures.
!
      ELSE IF (WORD.EQ.'CHMD') THEN
         CHMDT=.TRUE.
         INQUIRE(FILE='chmd.par',EXIST=YESNO)
         IF (YESNO) THEN
            OPEN(99,file='chmd.par')
            READ(99,'(A)') CHMDPAR
            CLOSE(99)
         ENDIF
         CALL READI(CHMDFREQ) 
!
      ELSE IF (WORD.EQ.'CHMOVE') THEN
         CHMOVET=.TRUE.
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'STANDARD') THEN
            CHSTANDARDT=.TRUE.
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'LOOPMODEL') THEN
            CHLOOPMODT=.TRUE.
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'CART') THEN
            CHCARTMOVET=.TRUE.
            CHNMAX=0.D0
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'NEIGHBOURS') THEN
            CHNEIGHBOURT=.TRUE.
         ELSE
            WRITE(MYUNIT,'(A)') 'Unrecognized 1st parameter for CHMOVE. Must be STANDARD, LOOPMODEL or CART.'
            STOP
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BB') THEN
               CHBBT=.TRUE.
            ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'SC') THEN
               CHSCT=.TRUE.
            ENDIF
         ELSE
            CHBBT=.TRUE.
            CHSCT=.TRUE.
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'CLUSTER') THEN
               CHCLUSTERT=.TRUE.
            ELSE
               CHCLUSTERT=.FALSE.
               WRITE(MYUNIT,'(A,F14.10)') 'Unrecognized 3rd argument for CHMOVE. Must be CLUSTER or none.'
               STOP
           ENDIF
         ENDIF

      ELSE IF (WORD.EQ.'CHPMAX') THEN
         CALL READF(CHPMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'CHPMAX=  ',CHPMAX

      ELSE IF (WORD.EQ.'CHPMIN') THEN
         CALL READF(CHPMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'CHPMIN=  ',CHPMIN

      ELSE IF (WORD.EQ.'CHNMAX') THEN
         CALL READF(CHNMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'CHNMAX=  ',CHNMAX
         IF (CHNMAX.LE.0.D0) CHCARTMOVET=.TRUE.

      ELSE IF (WORD.EQ.'CHNMIN') THEN
         CALL READF(CHNMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'CHNMIN=  ',CHNMIN

      ELSE IF (WORD.EQ.'TOMEGA') THEN
         OMEGAT=.TRUE.
         WRITE(MYUNIT,'(A)') 'TOMEGA set : peptide bonds will be twisted along with all other dihedrals'

      ELSE IF (WORD.EQ.'CHFREQ') THEN
         CALL READI(CHFREQ)
         WRITE(MYUNIT,'(A,I4,A)') 'Every ',CHFREQ,' steps TAKESTEPCH is called'
         IF(NITEMS.GT.2) CALL READI(CHFREQBB)
         IF(NITEMS.GT.3) CALL READI(CHFREQSC)

      ELSE IF (WORD.EQ.'CHRIGIDTRANS') THEN
         CHRIGIDTRANST=.TRUE.
         CALL READF(PTRANS)
         CALL READF(TRANSMAX)
         CALL READI(FTRANS)
         IF(NITEMS.GT.4) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BOX') BOXT=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SPHERE') SPHERET=.TRUE.
         IF (BOXT.AND.NITEMS.GT.5) CALL READF(BOXSIZE)
         IF (SPHERET.AND.NITEMS.GT.5) CALL READF(SPHERERAD)

      ELSE IF (WORD.EQ.'CHRIGIDROT') THEN
         CHRIGIDROTT=.TRUE.
         CALL READF(PROT)
         CALL READF(ROTMAX)
         CALL READI(FROT)

!
! Check for internal minimum in constraint terms for INTCONSTRAINT
!
         ELSE IF ((WORD.EQ.'CONINT').OR.(WORD.EQ.'QCIINT')) THEN
            CHECKCONINT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(INTMINFAC)
            WRITE(MYUNIT,'(A,G20.10)') ' keyword> Internal minima terms will be scaled by a factor of ',INTMINFAC 
!
! Absolute distance to allow before turning on constraint potential.
!
         ELSE IF (WORD.EQ.'CONCUTABS') THEN
            CONCUTABST=.TRUE.
            CONCUTFRACT=.FALSE.
            IF (NITEMS.GT.1) CALL READF(CONCUTABS)
!
! Fraction of constraint distance to allow before turning on constraint potential.
!
         ELSE IF (WORD.EQ.'CONCUTFRAC') THEN
            CONCUTFRACT=.TRUE.
            CONCUTABST=.FALSE.
            IF (NITEMS.GT.1) CALL READF(CONCUTFRAC)


      ELSE IF (WORD.EQ.'FIXEDEND') THEN
         FixedEndMoveT = .TRUE.
         IF (NITEMS.GT.1) CALL READF(PIVOTP)
         IF (NITEMS.GT.2) CALL READF(SIDECHAINP)

      ELSE IF (WORD.EQ.'OSASA') THEN
         OSASAT=.TRUE.
         CALL READF(RPRO)
         WRITE(MYUNIT,'(A)') 'OSASA set: solvent accessible surface area order parameter will be calculated'
         WRITE(MYUNIT,'(A,F4.1)') 'using probe radius ',RPRO

      ELSE IF (WORD.EQ.'ODIHE') THEN
         ODIHET=.TRUE.
         WRITE(MYUNIT,'(A)') 'ODIHE set: dihedral-angle order parameter will be calculated'
         WRITE(MYUNIT,'(A)') 'using the reference structure supplied in ref.crd'
 
      ELSE IF (WORD.EQ.'OEINT') THEN
         OEINTT=.TRUE.
         CALL READI(MON1(1))
         CALL READI(MON1(2))
         CALL READI(MON2(1))
         CALL READI(MON2(2))
         WRITE(MYUNIT,'(A)') 'OEINTT set: interaction energy between 2 peptides will be used as an order parameter'

      ELSE IF (WORD.EQ.'ORGYR') THEN
         ORGYT=.TRUE.
         WRITE(MYUNIT,'(A)') 'ORGYT set: radius of gyration will be calculated as an order parameter'

!     ELSE IF (WORD.EQ.'NORANDOM') THEN
!        NORANDOM=.TRUE.
!        IF (NITEMS.GT.1) CALL READF(RANDOMCUTOFF)

!     ELSE IF (WORD.EQ.'PERMDIHE') THEN
!        PERMDIHET=.TRUE.
!        DO J1=1,NITEMS-1
!           CALL READI(NDUM)
!           PERMDIHE(J1)=NDUM
!        ENDDO
!        NPERMDIHE=NITEMS-1
!        DO J1=1,NITEMS-1
!           print *,'PERMDIHE',PERMDIHE(J1)
!        ENDDO
!
!  End of CHARMM block.  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!  Sometimes have to modify the cold fusion limit when using high electric fields
!
      ELSE IF (WORD.EQ.'COLDFUSION') THEN
         CALL READF(COLDFUSIONLIMIT)

      ELSE IF (WORD.EQ.'COMPRESS') THEN
         COMPRESST=.TRUE.
         CALL READF(K_COMP)

! csw34> Apply compression to a rigid body system if the COM-COM distances exceed a threshold
      ELSE IF (WORD.EQ.'COMPRESSRIGID') THEN
         COMPRESSRIGIDT=.TRUE.
! Read in KCOMP_RIGID 
         IF (NITEMS.GT.1) CALL READF(KCOMP_RIGID)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body system compression enabled'
         WRITE(MYUNIT,'(A,F20.10)') ' COMPRESSRIGID> Compression force constant =',KCOMP_RIGID
! Read in COM-COM distance cutoff
         IF (NITEMS.GT.2) CALL READF(RIGIDCOMDIST)
         WRITE(MYUNIT,'(A,F20.10)') ' COMPRESSRIGID> Compression will be applied if COM-COM distance exceeds =',RIGIDCOMDIST
!
!  sf344> keyword to explore the potential field around a RIGID building block
!         (currently by changing only the Cartesian coordinates of the second particle on a grid, and evaluating the energy)
!       IMPORTANT: coords file should contain only 4 lines, with the position lines (first two) containing only zeros)
!
      ELSE IF (WORD.EQ.'RIGIDCONTOUR') THEN
          RIGIDCONTOURT=.TRUE.
          CALL READI(GRIDSIZE)
          DO J1=1,3
            DO J2=1,2
             CALL READF(CONTOURBOUNDS(J1,J2)) ! e.g. x between (-1,1) would be CONTOURBOUNDS(1,J2)
            END DO
          END DO
!
!  Alternative correlated moves: NCOOP nearest neighbours of a randomly selected atom
!  all move by the same amount. NCOOP default is 5.
!
      ELSE IF (WORD.EQ.'COOPMOVE') THEN
         COOP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCOOP)
         IF (NITEMS.GT.2) CALL READF(COOPCUT)

      ELSE IF (WORD.EQ.'CPMD') THEN
         CPMD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(MYUNIT,'(A)') ' ERROR - no CPMD system specified'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 12
            ENDIF
         ENDDO
12       CONTINUE
!
!  Calculate the continuous symmetry measure from Pinsky, Dryzun, Casanova, Alemany and Avnir,
!  J. Comp. Chem., 29, 2712, 2008.
!
      ELSE IF (WORD.EQ.'CSM') THEN
         CSMT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(CSMGP)
            CALL MYUPCASE(CSMGP)
         ELSE
            PRINT '(A)','keyword> ERROR - point group must be specified for CMS keyword'
            STOP
         ENDIF
         IF (NITEMS.GT.2) CALL READF(CSMEPS)
         IF (NITEMS.GT.3) THEN
            CSMGUIDET=.TRUE.
            CALL READA(CSMGUIDEGP)
            CALL MYUPCASE(CSMGUIDEGP)
         ENDIF

         IF (.NOT.PERMDIST) THEN ! set permdist if not done already
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
!        ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSWAP(NATOMS),SWAP1(NATOMS,2),SWAP2(NATOMS,2))
         ALLOCATE(NPERMSIZE(3*NATOMSALLOC),PERMGROUP(3*NATOMSALLOC),NSETS(3*NATOMSALLOC),SETS(NATOMSALLOC,70))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
!
!  Sanity checks!
!
               IF (NSETS(J1).GT.70) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 70'
!                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMSALLOC) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check! This condition is now allowed.
!  
!           DO J1=1,NDUMMY
!              DO J2=J1+1,NDUMMY
!                 IF (PERMGROUP(J2).EQ.PERMGROUP(J1)) THEN
!                    PRINT '(2(A,I8))','keyword> ERROR - atom ',PERMGROUP(J1),' appears more than once'
!                    STOP
!                 ENDIF
!              ENDDO
!           ENDDO
            CLOSE(1)
!
!  And yet another!
!  
            IF (NFREEZE.GT.0) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
            ENDIF
         ELSE
            NSETS(1:NATOMSALLOC)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMSALLOC ! all atoms can be permuted - default
            DO J1=1,NATOMSALLOC
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            WRITE(MYUNIT) '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(MYUNIT,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1)
            IF (NSETS(J1).GT.0) THEN
               WRITE(MYUNIT,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(MYUNIT,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' / '
                  ENDDO
                  IF (J2.LT.NSETS(J1)) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' ; '
               ENDDO
            ENDIF
            WRITE(MYUNIT,*) ' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO

      ELSE IF (WORD .EQ. 'CUDA') THEN
         CUDAT=.TRUE.
         IF(NITEMS .EQ. 2) THEN
            CALL READA(CUDAPOT)
         ELSE IF (NITEMS .EQ. 1) THEN
            WRITE(MYUNIT,'(A)') "keywords> You must specify a potential with keyword CUDA. "
            STOP
         ELSE
            WRITE(MYUNIT,'(A)') "keywords> The keyword CUDA only takes one argument. "
            STOP
         ENDIF

      ELSE IF (WORD .EQ. 'CUDATIME') THEN
         CUDATIMET=.TRUE.

      ELSE IF (WORD.EQ.'CUTOFF') THEN
         CUTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(CUTOFF)
         FINALCUTOFF=CUTOFF
         IF (NITEMS.GT.2) CALL READF(FINALCUTOFF)
         
      ELSE IF (WORD.EQ.'BOXCENTROID') THEN
         BOXCENTROIDT=.TRUE.
         CALL READF(BOXCENTROID_X(1))
         CALL READF(BOXCENTROID_X(2))
         CALL READF(BOXCENTROID_X(3))
         CALL READF(BOXCENTROID_DX(1))
         CALL READF(BOXCENTROID_DX(2))
         CALL READF(BOXCENTROID_DX(3))
         IF(NITEMS.GT.7) THEN
            DO J1=1,3
               CALL READI(J2)
               IF(J2==1) BOXCENTROID_DISCRETE(J1) = .TRUE.
            ENDDO
         ENDIF
         
      ELSE IF (WORD.EQ.'MIE_FIELD') THEN
         MIEFT=.TRUE.
         CALL READA(MIEF_FILENAME)
         IF(NITEMS.GT.2) THEN
            MIEF_CUTT=.TRUE.
            CALL READF(MIEF_RCUT)
         ENDIF
         IF(NITEMS.GT.3) THEN
            MIEF_PBCT=.TRUE.
            CALL READF(MIEF_BOX(1))
            CALL READF(MIEF_BOX(2))
            CALL READF(MIEF_BOX(3))
         ENDIF

!
! NOT DOCUMENTED - INTENTIONAL
!
      ELSE IF (WORD.EQ.'D5H') THEN
         FIELDT=.TRUE.
         D5HT=.TRUE.
         CALL READF(XX)
         FD5H=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF

! hk286 - damped group moves
! the structure is very similar to group rotation moves due to csw34 with appropriate changes
      ELSE IF (WORD.EQ.'DAMPEDGROUPMOVES') THEN
! Check the input file is present
         YESNO=.FALSE.
         INQUIRE(FILE='Dampedatomgroups',EXIST=YESNO)
         IF (YESNO) THEN
            DAMPEDGMOVET=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> DAMPED group moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: damped atom groups must be defined in Dampedatomgroups file'
            STOP
         ENDIF

         IF (NITEMS.GT.1) CALL READI(DMOVEFREQ)
! if the frequency is 0, we need to disable the move
         IF(DMOVEFREQ.EQ.0) THEN 
            DAMPEDGMOVET=.FALSE.
            WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of DAMPED moves set to 0 - moves DISABLED!'
         ENDIF
         NDGROUPS=0
         OPEN(UNIT=222,FILE='Dampedatomgroups',status='old')
         DO
            READ(222,*,IOSTAT=iostatus) CHECK1
            IF (iostatus<0) THEN
            CLOSE(222)
            EXIT
            ELSE IF (TRIM(ADJUSTL(check1)).EQ.'DGROUP') then
               NDGROUPS=NDGROUPS+1
            ENDIF
         END DO        
         CLOSE(222)
! Allocate atom group info arrays appropriately
         ALLOCATE(DATOMGROUPNAMES(NDGROUPS))
         ALLOCATE(DATOMGROUPAXIS(NDGROUPS,2))
         ALLOCATE(DATOMGROUPPROBA(NDGROUPS))
         ALLOCATE(DATOMGROUPSCALINGA(NDGROUPS))
         ALLOCATE(DATOMGROUPPROBB(NDGROUPS))
         ALLOCATE(DATOMGROUPSCALINGB(NDGROUPS))
         ALLOCATE(DATOMGROUPPROBC(NDGROUPS))
         ALLOCATE(DATOMGROUPSCALINGC(NDGROUPS))
         ALLOCATE(DATOMGROUPPROBD(NDGROUPS))
         ALLOCATE(DATOMGROUPSCALINGD(NDGROUPS))
         ALLOCATE(DATOMGROUPS(NDGROUPS,NATOMSALLOC))
         ALLOCATE(DATOMGROUPATT(NDGROUPS,NATOMSALLOC))
! Set safe defaults
         DATOMGROUPNAMES(:)='EMPTY'
         DATOMGROUPAXIS(:,:)=0
         DATOMGROUPPROBA(:)=0.0D0
         DATOMGROUPSCALINGA(:)=1.0D0
         DATOMGROUPPROBB(:)=0.0D0
         DATOMGROUPSCALINGB(:)=1.0D0
         DATOMGROUPPROBC(:)=0.0D0
         DATOMGROUPSCALINGC(:)=1.0D0
         DATOMGROUPPROBD(:)=0.0D0
         DATOMGROUPSCALINGD(:)=1.0D0
         DATOMGROUPS(:,:)=.FALSE.
         DATOMGROUPATT(:,:) = 1.0D0
! Read in group info
! Here is an example entry:
! DGROUP ALPHA 6 5 4 ProbA ScaleA ProbB ScaleB ProbC ScaleC ProbD ScaleD
! 1 attenuation
! 2 attenuation
! 3 attenuation
! 4 attenuation
! This says that group ALPHA has 4 members. The group will be rotated about the bond from atom 6->5
! with probability ProbA and rotations of -pi->+pi are to be scaled by ScaleA. ProbB and ScaleB
! corresponds to rotation in a random axis perpendicular to 6->5 axis. ProbC and ScaleC are translational
! moves along the 6->5 axis, while ProbD and ScaleD translational moves perpendicular to it.
! Finally, the group members are specified one per line, with the attenuation, between 0.0-1.0.
         OPEN(UNIT=222,FILE='Dampedatomgroups',status='unknown')
         WRITE(MYUNIT,*) 'keyword> Reading in atom groups for DAMPEDGROUPMOVES'
         DO J1=1,NDGROUPS
            READ(222,*) CHECK1,DATOMGROUPNAMES(J1),DATOMGROUPAXIS(J1,1),DATOMGROUPAXIS(J1,2),GROUPSIZE,DATOMGROUPPROBA(J1),
     &           DATOMGROUPSCALINGA(J1),DATOMGROUPPROBB(J1),DATOMGROUPSCALINGB(J1),
     &           DATOMGROUPPROBC(J1),DATOMGROUPSCALINGC(J1),DATOMGROUPPROBD(J1),DATOMGROUPSCALINGD(J1)
            CALL FLUSH(MYUNIT)
            IF (TRIM(ADJUSTL(CHECK1)).EQ.'DGROUP') THEN
               DO J2=1,GROUPSIZE
                  READ(222,*) GROUPATOM, GROUPATT
                  DATOMGROUPS(J1,GROUPATOM)=.TRUE. 
                  DATOMGROUPATT(J1,GROUPATOM) = GROUPATT
               END DO 
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword: ERROR! DAMPED group file not formatted correctly!'
               STOP
            ENDIF
            WRITE(MYUNIT,'(3A)') '<DGROUP ',TRIM(ADJUSTL(DATOMGROUPNAMES(J1))),'>'
            WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
            WRITE(MYUNIT,'(A,I4)') 'Size: ',GROUPSIZE
            WRITE(MYUNIT,'(A,2I5)') 'Atoms defining axis: ',DATOMGROUPAXIS(J1,1),DATOMGROUPAXIS(J1,2)
            WRITE(MYUNIT,'(A,4F5.2)') 'A - Selection probablity and scaling: ',DATOMGROUPPROBA(J1),DATOMGROUPSCALINGA(J1)
            WRITE(MYUNIT,'(A,4F5.2)') 'B - Selection probablity and scaling: ',DATOMGROUPPROBB(J1),DATOMGROUPSCALINGB(J1)
            WRITE(MYUNIT,'(A,4F5.2)') 'C - Selection probablity and scaling: ',DATOMGROUPPROBC(J1),DATOMGROUPSCALINGC(J1)
            WRITE(MYUNIT,'(A,4F5.2)') 'D - Selection probablity and scaling: ',DATOMGROUPPROBD(J1),DATOMGROUPSCALINGD(J1)
            WRITE(MYUNIT,'(A)') 'Members:'
            DO J2=1,NATOMSALLOC
               IF(DATOMGROUPS(J1,J2)) WRITE(MYUNIT,*) J2, DATOMGROUPATT(J1,J2)
            ENDDO
         ENDDO
         CLOSE(222)

      ELSE IF (WORD.EQ.'DBRENT') THEN
         DBRENTT=.TRUE.

      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.

      ELSE IF (WORD.EQ.'DEBUGss2029') THEN
         DEBUG=.TRUE.
!
!  Correlated random moves, in the sense that the magnitude of the step
!  decays exponentially as DECAYPARAM from a randomly chosen atom.
!
      ELSE IF (WORD.EQ.'DECAY') THEN
         DECAY=.TRUE.
         CALL READF(DECAYPARAM)
!
!  2D binary potential - Daan Frenkel #1
!

      ELSE IF (WORD.eq.'DF1') THEN
         DF1T=.TRUE.
!
! DFTB potential for silicon.
!
      ELSE IF (WORD.EQ.'DFTB') THEN
         DFTBT=.TRUE.
!
! Tiffany's DFTB potential for carbon, transplanted from OPTIM.
!
      ELSE IF (WORD.EQ.'DFTBC') THEN
         DFTBCT=.TRUE.
!
!  Initial guess for diagonal matrix elements in LBFGS.
!
      ELSE IF (WORD.EQ.'DGUESS') THEN
         CALL READF(DGUESS)
         IF (NITEMS.GT.2) CALL READF(INTDGUESS)

!      ELSE IF (WORD.EQ.'DIELEC') THEN
!         CALL READF(XX)
!         DPARAM=XX
!         WRITE(MYUNIT,'(A,F9.5)') ' Dielectric constant = ',DPARAM
!
!  NOT DOCUMENTED - INTENTIONAL
!
      ELSE IF (WORD.EQ.'DIFFRACT') THEN
         DIFFRACTT=.TRUE.

      ELSE IF (WORD.EQ.'DIPOLES') THEN
         DIPOLE=.TRUE.
      ELSE IF (WORD.EQ.'DISTOPT') THEN
         DISTOPT=.TRUE.
!
!  NOT DOCUMENTED - INTENTIONAL
!
      ELSE IF (WORD.EQ.'DL_POLY') THEN
         DL_POLY=.TRUE.
         CALL READI(NATOMS)

! vr274> Uses DMACRYS for potential calculations
      ELSE IF (WORD.EQ.'DMACRYS') THEN
         DMACRYST=.true.
         CALL DMACRYS_INITIALIZE_GMIN
         RETURN

      ELSE IF (WORD.EQ.'DMA_RANDOMSTART') THEN
         DMACRYS_RANDOMSTART=.true.
      ELSE IF (WORD.EQ.'DMA_EXPAND') THEN
         CALL READF(DMACRYS_EXPAND)
      ELSE IF (WORD.EQ.'DMA_LSTEP') THEN
         IF(NITEMS>7) THEN
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: wrong number of parameters in DMA_LSTEP'
            STOP
         ENDIF

         DO J1=1,NITEMS-1
           CALL READF(DMACRYS_LATTICE_STEP(J1))
         END DO

! csw34> DONTMOVE is an extension of the FREEZE idea. Atoms are replaced
! after the STEP, but their gradients are NOT set to 0, hence they may
! move during minimisation.
!
      ELSE IF (WORD.EQ.'DONTMOVE') THEN
         DONTMOVET=.TRUE.
         DO J1=1,NITEMS-1
            NDONTMOVE=NDONTMOVE+1
            CALL READI(NDUMMY)
            DONTMOVE(NDUMMY)=.TRUE.
         ENDDO
         
! csw34> DONTMOVEGROUP is analagous to FREEZEGROUP 
      ELSE IF (WORD.EQ.'DONTMOVEGROUP') THEN
         DONTMOVET=.TRUE.
         DONTMOVEGROUPT=.TRUE.
         CALL READI(DONTMOVECENTRE)
         CALL READF(DONTMOVERADIUS)
         IF(NITEMS.GT.3) CALL READA(DONTMOVEGROUPTYPE)
         
      ELSE IF (WORD.EQ.'DONTMOVEALL') THEN
         DONTMOVET=.TRUE.
         DONTMOVEALL=.TRUE.
         NDONTMOVE=NATOMSALLOC
         DO J1=1,NATOMSALLOC
            DONTMOVE(J1)=.TRUE.
            DONTMOVERES(J1)=.TRUE.
         ENDDO
! csw34
! Things are then moved using the DOMOVE and DOMOVERES keywords
! This is only a valid keyword if DONTMOVEALL is also specified
      ELSE IF ((WORD.EQ.'DOMOVE').AND.DONTMOVEALL) THEN
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            DONTMOVE(NDUMMY)=.FALSE.
            NDONTMOVE=NDONTMOVE-1
         ENDDO

      ELSE IF (WORD.EQ.'DONTMOVERES') THEN
         DONTMOVET=.TRUE.
         DONTMOVEREST=.TRUE.
! The FROZENRES array is then filled with the residue number from the
! data file
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            DONTMOVERES(NDUMMY)=.TRUE.
         ENDDO
         
         ELSEIF ((WORD.EQ.'DOMOVERES').AND.DONTMOVEALL) THEN
         DOMOVEREST=.TRUE.
! Set the right parts of the DONTMOVERES array to FALSE
         DO J1=1,NITEMS-1 
            CALL READI(NDUMMY)
            DONTMOVERES(NDUMMY)=.FALSE.
         ENDDO
         DONTMOVEREST=.TRUE.
         
      ELSE IF (WORD.EQ.'DUMP') THEN
         DUMPT=.TRUE.

      ELSE IF (WORD.EQ.'DUMPINT') THEN
         CALL READI(DUMPINT)

      ELSE IF (WORD.EQ.'DZUGUTOV') THEN
         DZTEST=.TRUE.
         CALL READF(DZP1)
         CALL READF(DZP2)
         CALL READF(DZP3)
         CALL READF(DZP4)
         CALL READF(DZP5)
         CALL READF(DZP6)
         CALL READF(DZP7)

      ELSE IF (WORD.EQ.'EAMAL') THEN
         EAMALT=.TRUE.

      ELSE IF (WORD.EQ.'EAMLJ') THEN
         EAMLJT=.TRUE.
         CALL READF(EAMLJA0)
         CALL READF(EAMLJBETA)
         CALL READF(EAMLJZ0)

      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(ECONV)

!
! Accumulation of thermodynamic statistics starting after Equil steps, 
! calculated thermodynamic properties is dumped every DumpEveryNthQuench quench.

!
      ELSE IF (WORD.EQ.'EQUILIBRATION') THEN
         CALL READI(EQUIL)
         CALL READI(DumpEveryNthQuench)
!
!  Steps using transition state search-type moves. Obsolete.
!  NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'EVSTEP') THEN
         EVSTEPT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NEVS)
         IF (NITEMS.GT.2) CALL READF(CEIG)
         IF (NITEMS.GT.3) CALL READI(NEVL)
         IF (NITEMS.GT.4) CALL READI(NVECTORS)
!
! Number of steps of equilibration for the replica above a reservoir
! after a successful exchange.
!
      ELSE IF (WORD.EQ.'EXEQ') THEN
         CALL READI(EXEQ)
!
!  NOT DOCUMENTED - INTENTIONAL
!
!
! csw34> Rigid body expansion moves. Scales the distance of the COM of each rigid body
!        from the system COM by EXPANDFACTOR every EXPANDRIGIDFREQ steps. The idea here 
!        is to provide more space for subsequent large moves in dense systems. If the third 
!        arguement 'NORMALISE' is given, the expansion vector for each rigid body is normalised 
!        before it is scaled, essentially translating all bodies by EXPANDFACTOR from the system COM. 
!
      ELSE IF (WORD.EQ.'EXPANDRIGID') THEN
         EXPANDRIGIDT=.TRUE.
! Read EXPANDRIGIDFREQ 
         IF (NITEMS.GT.1) CALL READI(EXPANDRIGIDFREQ)
! Read EXPANDFACTOR
         IF (NITEMS.GT.2) CALL READF(EXPANDFACTOR)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body expansion moves enabled'
         WRITE(MYUNIT,'(A,I2,A)') ' keyword> Rigid bodies will be expanded every ',EXPANDRIGIDFREQ,' steps' 
         WRITE(MYUNIT,'(A,F20.10)') ' EXPANDRIGID> System expansion factor =',EXPANDFACTOR
         IF (NITEMS.GT.3) THEN
! Read second arguement 
            CALL READA(WORD2)
            WORD2=TRIM(ADJUSTL(WORD2))
            IF (WORD2.EQ.'NORMALISE') THEN
               WRITE(MYUNIT,'(A)') ' EXPANDRIGID> Expandsion vextor will be normalised before scaling'
               NORMALISEEXPANDT=.TRUE.
            ENDIF
         ENDIF

      ELSE IF (WORD.EQ.'EXPFAC') THEN
         CALL READF(EFAC)
         IF (NITEMS.GT.2) CALL READF(EAMP)

! Commenting out this AMBER keyword that should be used only with PNM's hand-coded AMBER
!      ELSE IF (WORD.EQ.'FAKEWATER') THEN
!         FAKEWATER=.TRUE.
!         WRITE (MYUNIT,'(A)') '**********************************************************'
!         WRITE (MYUNIT,'(A)') '* DISTANCE DEPENDENT DIELECTRIC BEING USED - FAKE WATER! *'
!         WRITE (MYUNIT,'(A)') '**********************************************************'

      ELSE IF (WORD.EQ.'FAL') THEN
         FAL=.TRUE.
     
      ELSE IF (WORD == 'FEBH') THEN
         CALL READF(FETEMP)
         FEBHT = .TRUE.
         FE_FILE_UNIT = GETUNIT()
         OPEN(UNIT = FE_FILE_UNIT, FILE = 'free_energy', STATUS = 'REPLACE')
         WRITE(FE_FILE_UNIT, '(6A20)') '       Quench       ', '  Potential energy  ',
     &   '   Harmonic term    ', '    Free energy     ', '   Markov energy    ', '        Time        '
         WRITE(FE_FILE_UNIT, '(6A20)') ' ------------------ ', ' ------------------ ',
     &   ' ------------------ ', ' ------------------ ', ' ------------------ ', ' ------------------ '
         IF (NITEMS .GT. 2) THEN
             CALL READA(WORD2)
             WORD2 = TRIM(ADJUSTL(WORD2))
             IF (WORD2 .EQ. 'SPARSE') SPARSET = .TRUE.
         END IF
         IF (NITEMS .GT. 3) THEN
             CALL READF(ZERO_THRESH)
         END IF
 
      ELSE IF (WORD.EQ.'FIXBOTH') THEN
         IF (NITEMS.EQ.1) THEN
            FIXBOTH(1)=.TRUE.
            IF (NPAR.GT.1) THEN
               DO J1=2,NPAR
                  FIXBOTH(J1)=.TRUE.
               ENDDO
            ENDIF
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXBOTH(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FIXCOM') THEN
         FIXCOM=.TRUE.
!ds656> no need to allocate masses here, since it is done in main.F
!         IF (ALLOCATED(ATMASS)) DEALLOCATE(ATMASS)
!         ALLOCATE(ATMASS(NATOMSALLOC)) 
!     
!  Take hard sphere type moves.
!  T12FAC is the fraction of the first collision time to be used in HSMOVE
!
      ELSE IF (WORD.EQ.'FIXD') THEN
         FIXD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NHSMOVE)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(T12FAC)
         ENDIF

      ELSE IF (WORD.EQ.'FIXSTEP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXSTEP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXSTEP(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FIXTEMP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXTEMP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXTEMP(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FNI') THEN
         FNI=.TRUE.

      ELSE IF (WORD.EQ.'FORCERIGID') THEN
         RIGID=.TRUE.

      ELSE IF (WORD.EQ.'FRAUSI') THEN
         FRAUSIT=.TRUE.
!
!  Frozen atoms.
!
      ELSE IF (WORD.EQ.'FREEZE') THEN
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.TRUE.
         ENDDO
      ELSE IF (WORD.EQ.'FREEZENODES') THEN
         FREEZENODEST=.TRUE.
         CALL READF(FREEZETOL)

!
! sf344> unfreeze everything at the final quenches
!
      ELSE IF (WORD.EQ.'UNFREEZEFINALQ') THEN
        UNFREEZEFINALQ=.TRUE.
!
! csw34
! Frozen residues (to be converted to frozen atoms)
!
      ELSE IF (WORD.EQ.'FREEZERES') THEN
         FREEZE=.TRUE.
         FREEZERES=.TRUE.
! The FROZENRES array is then filled with the residue number from the
! data file
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.TRUE.
         ENDDO
! Finally, the frozen residue numbers are converted into frozen atom
! numbers. This is also forcefield dependant and must be done when we
! know which forcefield to use (i.e. in the CHARMM block above)

! hk286
      ELSE IF (WORD.EQ.'FREEZERIGIDBODY') THEN
         FREEZERIGIDBODYT = .TRUE.
! jdf43
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZENRIGIDBODY(NDUMMY)=.TRUE.
         ENDDO
! csw34
! Freezing EVERYTHING and then permitting small parts to move
! This is useful for large system to prevent the data file getting silly
      ELSEIF (WORD.EQ.'FREEZEALL') THEN
         FREEZE=.TRUE.
         FREEZEALL=.TRUE.
         NFREEZE=NATOMSALLOC
         DO J1=1,NATOMSALLOC
            FROZEN(J1)=.TRUE.
            FROZENRES(J1)=.TRUE.
         ENDDO
! csw34
! Things are then UNFROZEN using the UNFREEZE and UNFREEZERES keywords
! This is only a valid keyword if FREEZEALL is also specified
      ELSEIF ((WORD.EQ.'UNFREEZE').AND.FREEZEALL) THEN
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.FALSE.
            NFREEZE=NFREEZE-1
         ENDDO

      ELSEIF ((WORD.EQ.'UNFREEZERES').AND.FREEZEALL) THEN
         UNFREEZERES=.TRUE.
! Set the right parts of the FROZENRES array to FALSE
         DO J1=1,NITEMS-1 
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.FALSE.
         ENDDO
         FREEZERES=.TRUE.
! The FREEZERES routines for AMBER and CHARMM do the rest :)
!      
! Finnis-Sinclair potential coded by James Elliott
!
      ELSE IF (WORD.EQ.'FS') THEN
         FST=.TRUE.
         CALL READI(GATOM)

      ELSE IF (WORD.EQ.'G46') THEN
         G46=.TRUE.
         BLNT=.TRUE.
! mo361> Genetic algorithm keywords
!
! Genetic algorithm
!
      ELSE IF (WORD.EQ.'GA') THEN
         CALL READI(MYGA_NSTRUC)
         CALL READI(MYGA_NOFF)
         CALL READI(MYGA_GENS)
         GENALT=.TRUE.
!
! Genetic algorithm with variable-length basin-hopping optimisation
!
      ELSE IF (WORD.EQ.'GABHINCR') THEN
         CALL READF(MYGA_BH_INCR)
!
! Select 1- or 2-point crossover for GA
!
      ELSE IF (WORD.EQ.'GACROSS') THEN
         CALL READI(MYGA_CROSS)
!
! Dump genomes of whole population to file after every generation
!
      ELSE IF (WORD.EQ.'GADUMPPOP') THEN
         MYGA_DUMP_POP=.TRUE.
!
! Genetic algorithm duplicate predator
!
      ELSE IF (WORD.EQ.'GADUPPRED') THEN
         IF (NITEMS.GT.1) THEN
            CALL READF(MYGA_DUPLICATE_ETHRESH)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(MYGA_DUPLICATE_GTHRESH)
         ENDIF
!
! Genetic algorithm epoch operator
!
      ELSE IF (WORD.EQ.'GAEPOCH') THEN
         MYGA_L_EPOCH=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(MYGA_EPOCH_THRESH)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READI(MYGA_EPOCH_SAVE)
         IF (NITEMS.GT.3) THEN
            CALL READA(WORD2)
            WORD2=TRIM(ADJUSTL(WORD2))
            IF (WORD2.EQ."DUMP") THEN
               MYGA_EPOCH_DUMP=.TRUE.
            ENDIF
         ENDIF
         ENDIF
!
! Start genetic algorithm from chain structures
!
      ELSE IF (WORD.EQ.'GAINITCHAIN') THEN
         MYGA_L_CHAIN=.TRUE.
!
! Start genetic algorithm from sphere structures
!
      ELSE IF (WORD.EQ.'GAINITSPHERE') THEN
         MYGA_L_SPHERE=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(MYGA_SPHERE_RADIUS)
         ENDIF
!
! Genetic algorithm mutation rate
!
      ELSE IF (WORD.EQ.'GAMUTRATE') THEN
         CALL READF(MYGA_MUT_RATE)
!
! Genetic algorithm roulette selection
!
      ELSE IF (WORD.EQ.'GASELROUL') THEN
         MYGA_L_ROUL=.TRUE.
!
! Genetic algorithm tournament selection
!
      ELSE IF (WORD.EQ.'GASELTOURN') THEN
         CALL READI(MYGA_TOURN_SIZE)
         IF (MYGA_TOURN_SIZE.LT.2) THEN
            MYGA_TOURN_SIZE=2
            WRITE(MYUNIT,*) 'keyword> WARNING - GA tournament size must be at least 2.'
         ENDIF
!
! Allow twinning moves (crossover between two copies of the same parent)
!
      ELSE IF (WORD.EQ.'GATWIN') THEN
         MYGA_TWIN=.TRUE.
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'GAMMA') THEN
         CALL READF(GAMMA)
         TUNNELT=.TRUE.
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'GAUSS') THEN
         GAUSST=.TRUE.
         CALL READI(GMODES) ! number of nodes
!
! Grand canonical basin-hopping.
!
      ELSE IF (WORD.EQ.'GCBH') THEN
         GCBHT=.TRUE.
         CALL READF(GCMU) ! chemical potential
         CALL READI(GCNATOMS) ! maximum number of atoms
         IF (NITEMS.GT.3) CALL READI(GCINT) 
         IF (NITEMS.GT.4) CALL READI(GCRELAX) 
         IF (NITEMS.GT.5) CALL READF(GCPLUS) 
         WRITE(MYUNIT,'(A,I6)') 'keyyword> Interval for GC size change=',GCINT
         WRITE(MYUNIT,'(A,I6)') 'keyyword> Interval for GC relaxation=',GCRELAX
         IF (GCRELAX.GT.GCINT) THEN
            WRITE(MYUNIT,'(A,I6)') 'keyyword> *** ERROR - relaxation interval should be < size change interval'
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,2G20.10)') 'keyyword> GC atom addition/subtraction probabilities=',GCPLUS,1.0D0-GCPLUS
         
!ds656> Relative chemical potential(s) for semi-grand canonical BH.
      ELSE IF (WORD.EQ.'SEMIGRAND_MU') THEN
         SEMIGRAND_MUT = .TRUE.
         WRITE(MYUNIT,'(A)',ADVANCE='NO') 
     1        'keyword> Semi-grand chem. pots.:'
         ALLOCATE(SEMIGRAND_MU(NSPECIES(0)))
         SEMIGRAND_MU(1) = 0.0D0
         WRITE(MYUNIT,'(1X,F8.5)',ADVANCE='NO') SEMIGRAND_MU(1)
         DO J1=2,NSPECIES(0)
            CALL READF(SEMIGRAND_MU(J1))
            WRITE(MYUNIT,'(1X,F8.5)',ADVANCE='NO') SEMIGRAND_MU(J1)
         ENDDO
         WRITE(MYUNIT,*)
!------------------------------------------------------------
! ds656> Multicomponent Sutton-Chen with(out) PBCs and CUTOFF
!------------------------------------------------------------
      ELSE IF (WORD.EQ.'MSC') THEN
         MSCT = .TRUE.
         CALL PARSE_MSC_PARAMS(DATA_UNIT)
!------------------------------------
! ds656> Multicomponent Gupta systems
!------------------------------------
      ELSE IF (WORD.EQ.'MGUPTA') THEN
         MGUPTAT=.TRUE.
         CALL PARSE_MGUPTA_PARAMS(DATA_UNIT)
!-------------------------------------------
! ds656> Multicomponent Lennard-Jones system
!-------------------------------------------  
      ELSE IF (WORD.EQ.'MLJ') THEN
         MLJT=.TRUE.
         CALL PARSE_MLJ_PARAMS(DATA_UNIT)
!     
! General LJ for mixed systems
!
      ELSE IF (WORD.EQ.'GLJ') THEN
         GLJT=.TRUE.
         IF (NSPECIES(0) /= NITEMS-1) THEN
            WRITE(MYUNIT,'(A)') 'keywords> Inconsistent species count for GLJ!'
            STOP
         ENDIF
         !  the number of different species is given by NITEMS-1, 
         !  since we specify respective atom counts on this line
         DO J1=1,NSPECIES(0) ! loop over species
            CALL READI(NSPECIES(J1)) ! number of atoms for each species
         ENDDO
         NTYPEA = NSPECIES(1)
         IF (ALLOCATED(GLJEPS)) DEALLOCATE(GLJEPS)
         IF (ALLOCATED(GLJSIG)) DEALLOCATE(GLJSIG)
         ALLOCATE(GLJEPS(NSPECIES(0),NSPECIES(0)))
         ALLOCATE(GLJSIG(NSPECIES(0),NSPECIES(0)))
         WRITE(MYUNIT,'(A,I8,A)') 'keyword> Mixed LJ potential with ',NSPECIES(0),' particle types. Numbers of each particle:'
         NTYPEA=NSPECIES(1)
         DO J1=1,NSPECIES(0)
            WRITE(MYUNIT,'(I6)',ADVANCE='NO') NSPECIES(J1)
         ENDDO
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A)') 'keyword> Epsilon parameters:'
         DO J1=1,NSPECIES(0)
            READ(DATA_UNIT,*) (GLJEPS(J2,J1),J2=1,J1)
            DO J2=1,J1
               GLJEPS(J1,J2)=GLJEPS(J2,J1)
               WRITE(MYUNIT,'(F15.4)',ADVANCE='NO') GLJEPS(J2,J1)
            ENDDO
            WRITE(MYUNIT,'(A)') ' '
         ENDDO
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A)') 'keyword> Sigma parameters:'
         DO J1=1,NSPECIES(0)
            READ(DATA_UNIT,*) (GLJSIG(J2,J1),J2=1,J1)
            DO J2=1,J1
               GLJSIG(J1,J2)=GLJSIG(J2,J1)
               WRITE(MYUNIT,'(F15.4)',ADVANCE='NO') GLJSIG(J2,J1)
            ENDDO
            WRITE(MYUNIT,'(A)') ' '
         ENDDO
      ELSE IF (WORD.EQ.'GROUND') THEN
         GROUND=.TRUE.

!--------------------------------!
! hk286 > Generalised Thomson    !
!--------------------------------!

      ELSE IF (WORD .EQ. 'GTHOMSON') THEN
         GTHOMSONT = .TRUE.
         THOMSONT=.TRUE.
         RIGID = .TRUE.
         CALL READI(GTHOMMET)
         CALL READF(GTHOMSONZ)
         IF (NITEMS.GT.3) THEN
            CALL READF(GTHOMSONC)
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READF(GTHOMSONC2)
            IF ( (GTHOMMET .EQ. 2) .AND. (GTHOMSONC2 > 0.0D0) ) THEN
               GTHOMSONZ = LOG((GTHOMSONC2/GTHOMSONC)+ SQRT((GTHOMSONC2/GTHOMSONC)**2-1.0D0))*GTHOMSONC
!               PRINT *, GTHOMSONZ
            ENDIF
         ENDIF
         IF (NITEMS.GT.5) CALL READF(VUNDULOID)
     
         IF ((GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4)) THEN
            CALL CONVERTUNDULOIDPARAMETERS(VUNDULOID/2.0D0)
         ENDIF

!         ODDCHARGE=1.0D0
!         IF (NITEMS.GT.6) CALL READF(ODDCHARGE)
         
         IF ((GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4)) THEN
            CALL FINDNGZ()
         ENDIF

         CALL INIGTHOMSON()

      ELSE IF (WORD .EQ. 'GTHOMSONPOT') THEN
         CALL READI(GTHOMPOT)
         IF (NITEMS.GT.2) THEN
            CALL READF(GThomsonSigma)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(GThomsonRho)
         ENDIF

      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)

! csw34> Group rotation moves (now for both AMBER and CHARMM! 
      ELSE IF (WORD.EQ.'GROUPROTATION') THEN

! csw34> Check the group file is present
         YESNO=.FALSE.
         INQUIRE(FILE='atomgroups',EXIST=YESNO)
         IF (YESNO) THEN
            GROUPROTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> AMBER group rotation moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: atom groups must be defined in atomgroups file'
            STOP
         ENDIF
         IF (NITEMS.GT.1) CALL READI(GROUPROTFREQ)
! csw34> if the frequency is 0, we need to disable the moves to present
! a divide by 0!
         IF(GROUPROTFREQ.EQ.0) THEN 
            GROUPROTT=.FALSE.
            WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of GROUPROTATION moves set to 0 - moves DISABLED!'
         ENDIF
! Specify GROUPROTATION move scaling mode amd inform the user
         IF (NITEMS.GT.2) THEN 
            CALL READA(GR_SCALEMODE)
            IF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEROT') THEN
               GR_SCALEROT=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION amplitudes will be scaled'
            ELSEIF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEPROB') THEN  
               GR_SCALEPROB=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION selection probabilities will be scaled'
            ELSEIF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEBOTH') THEN
               GR_SCALEROT=.TRUE.
               GR_SCALEPROB=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION amplitudes and selection probabilities will be scaled'
            ENDIF
         ENDIF
! Read atom offset for group definitions in atomgroups 
         IF (NITEMS.GT.3) CALL READI(GROUPOFFSET)
! csw34> Figure out how many atom groups have been defined
         NGROUPS=0
         OPEN(UNIT=222,FILE='atomgroups',status='old')
         DO
            READ(222,*,IOSTAT=iostatus) CHECK1
            IF (iostatus<0) THEN
            CLOSE(222)
            EXIT
            ELSE IF (TRIM(ADJUSTL(check1)).EQ.'GROUP') then
               NGROUPS=NGROUPS+1
            ENDIF
         END DO        
         CLOSE(222)
! csw34> Allocate atom group info arrays appropriately
         ALLOCATE(ATOMGROUPNAMES(NGROUPS))
         ALLOCATE(ATOMGROUPAXIS(NGROUPS,2))
         ALLOCATE(ATOMGROUPPSELECT(NGROUPS))
         ALLOCATE(ATOMGROUPSCALING(NGROUPS))
         ALLOCATE(ATOMGROUPS(NGROUPS,NATOMSALLOC))
! csw34> Set safe defaults
         ATOMGROUPS(:,:)=.FALSE.
         ATOMGROUPNAMES(:)='EMPTY'
         ATOMGROUPAXIS(:,:)=0
         ATOMGROUPSCALING(:)=1.0D0
         ATOMGROUPPSELECT(:)=1.0D0
! csw34> Read in group info
! Here is an example entry:
! GROUP OME 6 5 4 1.0
! 1
! 2
! 3
! 4
! This says that group OME is to be rotated about the bond from atom 6->5.
! There are 4 atoms in the OME group. Rotations of -pi->+pi are to be scaled by 1.0. 
! Finally, the group members are specified one per line
         OPEN(UNIT=222,FILE='atomgroups',status='unknown')
         WRITE(MYUNIT,*) 'keyword> Reading in atom groups for GROUPROTATION'
         IF(GROUPOFFSET.NE.0) WRITE(MYUNIT,*) 'keyword> Group atom numbering offset by ',GROUPOFFSET
         DO J1=1,NGROUPS
            READ(222,*) CHECK1,ATOMGROUPNAMES(J1),AXIS1,AXIS2,GROUPSIZE,ATOMGROUPSCALING(J1),
     &                       ATOMGROUPPSELECT(J1) 
            ATOMGROUPAXIS(J1,1)=AXIS1+GROUPOFFSET
            ATOMGROUPAXIS(J1,2)=AXIS2+GROUPOFFSET
            CALL FLUSH(MYUNIT)
            IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') THEN
               DO J2=1,GROUPSIZE
                  READ(222,*) GROUPATOM
                  IF(GROUPOFFSET.GT.0) GROUPATOM=GROUPATOM+GROUPOFFSET
! hk286 - add bound checks
                  IF (GROUPATOM > NATOMSALLOC) THEN
                     WRITE(MYUNIT,'(A)') ' keyword: ERROR! GROUPATOM > NATOMSALLOC'
                  ENDIF
! hk286
                  ATOMGROUPS(J1,GROUPATOM)=.TRUE. 
               END DO 
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword: ERROR! Group file not formatted correctly!'
               STOP
            ENDIF
            WRITE(MYUNIT,'(3A)') '<GROUP ',TRIM(ADJUSTL(ATOMGROUPNAMES(J1))),'>'
            WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
            WRITE(MYUNIT,'(A,I4)') 'Size: ',GROUPSIZE
            WRITE(MYUNIT,'(A,2I6)') 'Atoms defining axis: ',ATOMGROUPAXIS(J1,1),ATOMGROUPAXIS(J1,2)
            WRITE(MYUNIT,'(A,F5.2)') 'Rotation scaling: ',ATOMGROUPSCALING(J1)
            WRITE(MYUNIT,'(A,F5.4)') 'Selection probablity: ',ATOMGROUPPSELECT(J1)
            WRITE(MYUNIT,'(A)') 'Members:'
            DO J2=1,NATOMSALLOC
               IF(ATOMGROUPS(J1,J2)) WRITE(MYUNIT,*) J2
            ENDDO
         ENDDO
         CLOSE(222)

! Keyword to suppress output of group rotation moves
      ELSE IF (WORD .EQ. 'QUIETGROUPROT') THEN
         GROUPROT_SUPPRESS = .TRUE.


! ab2111> Group rotation moves for dihedral angles 
      ELSE IF (WORD.EQ.'DIHEDRALROTATION') THEN
! Check the group file is present
         YESNO=.FALSE.
         INQUIRE(FILE='dihedralgroups',EXIST=YESNO)
         IF (YESNO) THEN
            DIHEDRALROTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> AMBER dihedral group rotation moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: atom groups must be defined in dihedralgroups file'
            STOP
         ENDIF
         IF (NITEMS.GT.1) CALL READI(DIHEDRALROTFREQ)
! if the frequency is 0, we need to disable the moves to present
! a divide by 0!
         IF(DIHEDRALROTFREQ.EQ.0) THEN 
            DIHEDRALROTT=.FALSE.
            WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of DIHEDRALROTATION moves set to 0 - moves DISABLED!'
         ENDIF
         IF (NITEMS.GT.2) CALL READI(DIHEDRALOFFSET)

! Figure out how many atom groups have been defined
         NDIHEDRALGROUPS=0
         OPEN(UNIT=223,FILE='dihedralgroups',status='old')
         DO
            READ(223,*,IOSTAT=iostatus) CHECK1
            IF (iostatus<0) THEN
            CLOSE(223)
            EXIT
            ELSE IF (TRIM(ADJUSTL(check1)).EQ.'GROUP') then
               NDIHEDRALGROUPS=NDIHEDRALGROUPS+1
            ENDIF
         END DO        
         CLOSE(223)
! Allocate atom group info arrays appropriately
         ALLOCATE(DIHEDRALGROUPNAMES(NDIHEDRALGROUPS))
         ALLOCATE(ANGLETYPE(NDIHEDRALGROUPS))
         ALLOCATE(DIHEDRALGROUPAXIS(NDIHEDRALGROUPS,4))
         ALLOCATE(DIHEDRALGROUPPSELECT(NDIHEDRALGROUPS))
         ALLOCATE(DIHEDRALGROUPSCALING(NDIHEDRALGROUPS))
         ALLOCATE(DIHEDRALGROUPS(NDIHEDRALGROUPS,NATOMSALLOC))
!  Set safe defaults
         DIHEDRALGROUPS(:,:)=.FALSE.
         DIHEDRALGROUPNAMES(:)='EMPTY'
         ANGLETYPE(:)='EMPTY'
         DIHEDRALGROUPAXIS(:,:)=0
         DIHEDRALGROUPSCALING(:)=1.0D0
         DIHEDRALGROUPPSELECT(:)=1.0D0
!  Read in group info
! Here is an example entry:
! GROUP OME <phi> 8 7 6 5 Natoms Scaling RotProbability
! 1
! 2
! 3
! 4
! This says that group OME, angle type <phi> (used for printing purposes only) 
! is to be rotated about the dihedral angle defined by atoms 8 7 6 5.
! Natoms = number of atoms to be rotated ( number of subsequent lines before next GROUP declaration
! Scaling : uniform perturbation steps within range phi-pi*Scaling and phi+pi*Scaling
! RotProbability : probability of performing this trial move per MC step

! Finally, the group members are specified one per line
         OPEN(UNIT=223,FILE='dihedralgroups',status='unknown')
         WRITE(MYUNIT,*) 'keyword> Reading in atom groups for DIHEDRALROTATION'
         IF(DIHEDRALOFFSET.NE.0) WRITE(MYUNIT,*) 'keyword> Group atom numbering offset by ',DIHEDRALOFFSET
         DO J1=1,NDIHEDRALGROUPS
            READ(223,*) CHECK1,DIHEDRALGROUPNAMES(J1),ANGLETYPE(J1),A1,A2,A3,A4,
     &       dihedralgroupsize,DIHEDRALGROUPSCALING(J1),DIHEDRALGROUPPSELECT(J1) 
!            READ(223,*) CHECK1,DIHEDRALGROUPNAMES(J1),A1,A2,A3,A4,GROUPSIZE,ATOMGROUPSCALING(J1),
!     &                       ATOMGROUPPSELECT(J1) 
            DIHEDRALGROUPAXIS(J1,1)=A1+DIHEDRALOFFSET
            DIHEDRALGROUPAXIS(J1,2)=A2+DIHEDRALOFFSET
            DIHEDRALGROUPAXIS(J1,3)=A3+DIHEDRALOFFSET
            DIHEDRALGROUPAXIS(J1,4)=A4+DIHEDRALOFFSET
            CALL FLUSH(MYUNIT)
            IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') THEN
               DO J2=1,dihedralgroupsize
                  READ(223,*) GROUPATOM
                  IF(dihedraloffset.GT.0) GROUPATOM=GROUPATOM+dihedraloffset
! - add bound checks
                  IF (GROUPATOM > NATOMSALLOC) THEN
                     WRITE(MYUNIT,'(A)') ' keyword: ERROR! GROUPATOM > NATOMSALLOC'
                  ENDIF
                  DIHEDRALGROUPS(J1,GROUPATOM)=.TRUE. 
               END DO 
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword: ERROR! Group file not formatted correctly!'
               STOP
            ENDIF
            WRITE(MYUNIT,'(3A)') '<GROUP ',TRIM(ADJUSTL(DIHEDRALGROUPNAMES(J1))),'>'
            WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
            WRITE(MYUNIT,'(A,I4)') 'Size: ',dihedralgroupsize
            WRITE(MYUNIT,'(A,4I5)') 'Atoms defining dihedral: ',DIHEDRALGROUPAXIS(J1,1),DIHEDRALGROUPAXIS(J1,2),
     &       DIHEDRALGROUPAXIS(J1,2),DIHEDRALGROUPAXIS(J1,2)
            WRITE(MYUNIT,'(A,F5.2)') 'Rotation scaling: ',DIHEDRALGROUPSCALING(J1)
            WRITE(MYUNIT,'(A,F5.2)') 'Selection probablity: ',DIHEDRALGROUPPSELECT(J1)
            WRITE(MYUNIT,'(A)') 'Members:'
            DO J2=1,NATOMSALLOC
               IF(DIHEDRALGROUPS(J1,J2)) WRITE(MYUNIT,*) J2
            ENDDO
         ENDDO
         CLOSE(223)
      
      ELSE IF (WORD.EQ.'BGSMOVE') THEN
            BGSMOVE=.TRUE. 
            ! check that DIHEDRALROTATION is turned on
            IF (DIHEDRALROTT.EQV..FALSE.) THEN
               WRITE(MYUNIT,*) "DIHEDRALROTATION must be turned on for BGSMOVE"
               STOP
            ENDIF
         IF (NITEMS.GT.1) THEN
               CALL READF(bgsb1)
               CALL READF(bgsb2)
               CALL READF(pselectBGS)
         ENDIF
         IF (NITEMS.GT.4) THEN
               CALL READI(NDIHEDRAL_BB_GROUPS)
            ELSE
               NDIHEDRAL_BB_GROUPS = NDIHEDRALGROUPS
            ENDIF
            !WRITE(MYUNIT,*) "HELIX MOVES turned on, with (Phi0, Psi0)=",PHI0,PSI0,"*pi rad; (PHIk,PSIk)=",PHIk,PSIk

      ELSE IF (WORD.EQ.'ROTAMERMOVE') THEN
         STOP 'Rotamer moves not implemented yet.'
         ROTAMER_MOVET = .TRUE.
         IF (NITEMS .GT. 1) THEN
            CALL READA(ROTAMER_SCRIPT)
         END IF
         CALL ROTAMER_INIT()
      ELSE IF (WORD.EQ.'GUIDE') THEN
         CALL READF(GUIDECUT)

      ELSE IF (WORD.EQ.'GUPTA') THEN
         GUPTAT=.TRUE.
         CALL READI(GATOM)

!
! js850> Add a harmonic field to the given particles so they dont stray very far
!        from their initial conditions.
!
      ELSE IF (WORD.EQ.'HARMONICF') THEN
         HARMONICF=.TRUE.
         HARMONICSTR=25.D0
         DO J1=1,NITEMS-1
            !NHARMONIC=NHARMONIC+1
            CALL READI(NDUMMY)
            HARMONICFLIST(NDUMMY)=.TRUE.
         ENDDO
!
! js850> IF HARMONICDONTMOVE, then have step size = 0 for the atoms with a
!        harmonic constraint
!
      ELSE IF (WORD.EQ.'HARMONICDONTMOVE') THEN
         HARMONICDONTMOVE=.TRUE.
!
! csw34> Hydrogen-bond matrix calculation
!
      ELSE IF (WORD.EQ.'HBONDMATRIX') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='hbond_matrix.sh',EXIST=YESNO)
         IF (YESNO) THEN
            HBONDMATRIX=.TRUE.
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: NEED hbond_matrix.sh SCRIPT TO USE HBONDMATRIX'
            STOP
         ENDIF
         CALL READA(HBONDDONORSACCEPTORS)
         CALL READA(HBONDRESIDUES)
         
! Count how many residues are of interest
         OPEN(UNIT=20,FILE=TRIM(ADJUSTL(HBONDRESIDUES)),STATUS='OLD')
         HBONDNRES=0
         DO
            READ(20,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               HBONDNRES = HBONDNRES + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         CLOSE(20)
! Print some info
         WRITE(MYUNIT,'(A,I5,A)') ' keyword> Hydrogen-bond matrix includes ',HBONDNRES,' residues' 
! Optional mode switch for hydrogen-bond matrix analysis
         IF (NITEMS.GT.3) THEN
            CALL READA(HBONDTYPE)
! If the third arguement is set to REJECT - all steps which move to a new group are rejected
            IF (HBONDTYPE.EQ.'REJECT') HBONDACCEPT=.FALSE.
         ENDIF

! csw34> Consider only the ligand when accepting/rejecting based on hydrogen-bonding
      ELSE IF (WORD.EQ.'HBONDLIGAND') THEN
            IF (.NOT.HBONDMATRIX) THEN
               WRITE(MYUNIT,'(A)') ' keyword> ERROR: HBONDLIGAND must appear after HBONDMATRIX in data!'
               STOP
            ENDIF
            HBONDLIGAND=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READI(HBONDLIGANDN)
            ELSE
               HBONDLIGANDN=HBONDNRES
            ENDIF 
            WRITE(MYUNIT,'(A)') ' keyword> HBONDLIGAND: Only bonds to the ligand will be considered in HBONDMATRIX grouping'


! csw34> Set the window size for the soft cutoffs to custom values
      ELSE IF (WORD.EQ.'HBONDSOFTCUT') THEN
         IF (HBONDMATRIX) THEN  
            CALL READA(HBONDDCUTLOOSE)
            CALL READA(HBONDDCUTTIGHT)
            WRITE(MYUNIT,'(4A)') ' keyword> HBONDSOFTCUT: Soft distance cutoff set to ',TRIM(ADJUSTL(HBONDDCUTLOOSE)),
     &                          ' - ',TRIM(ADJUSTL(HBONDDCUTTIGHT))
            CALL READA(HBONDACUTLOOSE)
            CALL READA(HBONDACUTTIGHT)
            WRITE(MYUNIT,'(4A)') ' keyword> HBONDSOFTCUT: Soft angle cutoff set to ',TRIM(ADJUSTL(HBONDACUTLOOSE)),
     &                          ' - ',TRIM(ADJUSTL(HBONDACUTTIGHT))
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: HBONDSOFTCUT can only be used with HBONDMATRIX'
         ENDIF

!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'HISTSMOOTH') THEN
         CALL READI(NSpline)
!
! Parameters of the temperature range on which to calculate thermodynamic properties in Basin Sampling
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'HISTTEMP') THEN
         CALL READF(MinimalTemperature)
         CALL READF(MaximalTemperature)
         CALL READI(NTempPoints)
!
! csw34> Hybrid rigid body/all-atom minimisation - only works with RIGIDINIT, the generalised rigid
!        body framework. EPSRIGID is the RMS convergence condition for the rigid body part of the
!        minimisation. Once this is complete, the RMS condition reverts to that from SLOPPYCONV for
!        the rest of the minimisation. When using this keyword, final quenches are always done
!        atomistically.
!
      ELSE IF (WORD.EQ.'HYBRIDMIN') THEN
! Sanity check
         IF (.NOT.RIGIDINIT) THEN
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: HYBRIDMIN can only be used with RIGIDINIT!'
            STOP
         ENDIF
! Read in rigid body minimisation convergence criterion
         CALL READF(EPSRIGID)
         HYBRIDMINT=.TRUE.
         WRITE(MYUNIT,'(A)') ' keyword> Using hybrid rigid body/all-atom minimisation'
         WRITE(MYUNIT,'(A,F20.10)') ' HYBRIDMIN> Rigid body RMS force target= ',EPSRIGID
         WRITE(MYUNIT,'(A)') ' HYBRIDMIN> Final quenches will be done atomistically'
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'IH') THEN
         FIELDT=.TRUE.
         IHT=.TRUE.
         CALL READF(XX)
         FIH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'INTMIN') THEN
         INTMINT=.TRUE.
!
! Images for INTCONSTRAINT
!
         ELSE IF ((WORD.EQ.'INTIMAGE').OR.(WORD.EQ.'QCIIMAGE')) THEN
            IF (NITEMS.GT.1) CALL READF(IMSEPMIN)
            IF (NITEMS.GT.2) CALL READF(IMSEPMAX)
            IF (NITEMS.GT.3) CALL READI(INTIMAGE)
            IF (NITEMS.GT.4) CALL READI(MAXINTIMAGE)
            IF (NITEMS.GT.5) CALL READI(INTNTRIESMAX)
            IF (NITEMS.GT.6) CALL READI(INTIMAGEINCR)
            IF (NITEMS.GT.7) CALL READI(INTIMAGECHECK)
!
! Use constraint potential for initial interpolation in each cycle.
!
         ELSE IF ((WORD.EQ.'INTCONSTRAINT').OR.(WORD.EQ.'QCI')) THEN
            INTCONSTRAINTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
            IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
            IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
            IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
            IF (NITEMS.GT.5) CALL READF(INTCONFRAC)
            IF (NITEMS.GT.6) CALL READI(INTCONSEP)
            IF (NITEMS.GT.7) CALL READI(INTREPSEP)
            IF (NITEMS.GT.8) CALL READI(INTSTEPS1)
            IF (NITEMS.GT.9) CALL READI(INTCONSTEPS)
            IF (NITEMS.GT.10) CALL READI(INTRELSTEPS)
            IF (NITEMS.GT.11) CALL READF(MAXCONE)
            IF (NITEMS.GT.12) CALL READF(INTRMSTOL)

            INQUIRE(FILE='finish',EXIST=YESNO)
            IF (YESNO) THEN
               LUNIT=GETUNIT()
               OPEN (UNIT=LUNIT,FILE='finish',STATUS='OLD')
               READ(LUNIT,*) FINISH(1:3*NATOMS)
               CLOSE(LUNIT)
               WRITE(MYUNIT,'(A,I6,A)') ' keyword> Read second end point from file finish'
            ELSE
               WRITE(MYUNIT,'(A,I6,A)') ' keyword> No finish file found - assuming a QCIPOT run'
            ENDIF
!
! Parse the congeom file to obtain number of reference minima for setting up
! constraints and get their geometries.
!
            INQUIRE(FILE='congeom.dat',EXIST=YESNO)

            IF (YESNO) THEN
               CONDATT=.TRUE.
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='congeom.dat',STATUS='OLD')
               READ(LUNIT,*) NCONGEOM
               !
               ! we only need the first reference if we have congeom.dat
               ! However, CONGEOM is passed to make_conpot and declared with first dimension NCPFIT there!
               !
               ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
               CONGEOM(2:NCONGEOM,1:3*NATOMS)=0.0D0
               READ(LUNIT,*) CONGEOM(1,1:3*NATOMS)
               READ(LUNIT,*) NCONSTRAINTFIX
               ALLOCATE(CONIFIX(NCONSTRAINTFIX),CONJFIX(NCONSTRAINTFIX),
     &         CONDISTREFFIX(NCONSTRAINTFIX),CONCUTFIX(NCONSTRAINTFIX))
               READ(LUNIT,*) CONIFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONJFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONDISTREFFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONCUTFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) NREPULSIVEFIX
               ALLOCATE(REPIFIX(NREPULSIVEFIX),REPJFIX(NREPULSIVEFIX),REPCUTFIX(NREPULSIVEFIX))
               READ(LUNIT,*) REPIFIX(1:NREPULSIVEFIX)
               READ(LUNIT,*) REPJFIX(1:NREPULSIVEFIX)
               READ(LUNIT,*) REPCUTFIX(1:NREPULSIVEFIX)
               CLOSE(LUNIT)
               WRITE(MYUNIT,'(A)') ' keyword> Constraint potential parameters read from file congeom.dat'
               INTCONMAX=NCONSTRAINTFIX
               NREPMAX=NREPULSIVEFIX
               ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX))
               ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
               ALLOCATE(CONACTIVE(NCONSTRAINTFIX))
            ELSE
               INQUIRE(FILE='congeom',EXIST=CONFILE)
               NCONGEOM=0
               NCONGEOM=0
               IF (.NOT.CONFILE) THEN
                  WRITE(MYUNIT,'(A)') ' keyword> WARNING *** no congeom file found. Will use end point minima only.'
               ELSE
                  LUNIT=GETUNIT()
                  OPEN(LUNIT,FILE='congeom',STATUS='OLD')
                  DO
                     READ(LUNIT,*,END=864) DUMMY1(1)
                     NCONGEOM=NCONGEOM+1
                  ENDDO
864               CONTINUE
                  NCONGEOM=NCONGEOM/NATOMS
                  ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
                  REWIND(LUNIT)
                  DO J1=1,NCONGEOM
                     READ(LUNIT,*) CONGEOM(J1,1:3*NATOMS)
                  ENDDO
                  CLOSE(LUNIT)
                  WRITE(MYUNIT,'(A,I6,A)') 'keyword> Read ',NCONGEOM,' reference geometries from congeom file'
                  IF (NCONGEOM.LT.2) WRITE(MYUNIT,'(A)') ' WARNING *** insufficient reference geometries - using end point minima'
               ENDIF
            ENDIF
!
! Use the quasi-continuous metric for connection attempts, instead of distance.
!
            INTERPCOSTFUNCTION=.TRUE.
         ELSE IF (WORD.EQ.'INTFREEZE') THEN
            INTFREEZET=.TRUE.
            IF (NITEMS.GT.1) CALL READF(INTFREEZETOL)
            IF (NITEMS.GT.2) CALL READI(INTFREEZEMIN)
!
! Use interpolation potential for LJ.
!
         ELSE IF (WORD.EQ.'INTLJ') THEN
            I NTLJT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(INTLJSTEPS)
            IF (NITEMS.GT.2) CALL READF(INTLJDEL)
            IF (NITEMS.GT.3) CALL READF(INTLJTOL)
            IF (NITEMS.GT.4) CALL READI(INTIMAGE)
            IF (NITEMS.GT.5) CALL READF(INTLJEPS)
!
! Parse the congeom file to obtain number of reference minima for setting up
! constraints and get their geometries.
!
            INQUIRE(FILE='congeom',EXIST=CONFILE)
            NCONGEOM=0
            IF (.NOT.CONFILE) THEN
               PRINT '(A)',' keyword> WARNING *** no congeom file found. Will use end point minima only.'
            ELSE
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='congeom',STATUS='OLD')
               DO
                  READ(LUNIT,*,END=863) DUMMY1(1)
                  NCONGEOM=NCONGEOM+1
               ENDDO
863            CONTINUE
               NCONGEOM=NCONGEOM/NATOMS
               ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
               REWIND(LUNIT)
               DO J1=1,NCONGEOM
                  READ(LUNIT,*) CONGEOM(J1,1:3*NATOMS)
               ENDDO
               CLOSE(LUNIT)
               PRINT '(A,I6,A)',' keyword> Read ',NCONGEOM,' reference geometries from congeom file'
               IF (NCONGEOM.LT.2) PRINT '(A)',' WARNING *** insufficient reference geometries - using end point minima'
            ENDIF
!
! Use the quasi-continuous metric for connection attempts, instead of distance.
!
            INTERPCOSTFUNCTION=.TRUE.
!
! Only include spring gradients for active atoms.
!
      ELSE IF (WORD.EQ.'INTSPRINGACTIVE') THEN
         INTSPRINGACTIVET=.TRUE.
!
! KINT: force constant for springs in INTCONSTRAINT calculations.
! Default zero.
!
         ELSE IF (WORD.EQ.'KINT') THEN
            CALL READF(KINT)
 
      ELSE IF (WORD.EQ.'JM') THEN
         JMT=.TRUE.

      ELSE IF (WORD.EQ.'JUMPMOVE') THEN
         CALL READI(IX)
         JUMPMOVE(IX)=.TRUE.
         CALL READI(JUMPTO(IX))
         JDUMP(JUMPTO(IX))=.TRUE.
         IF (NITEMS.GT.3) CALL READI(JUMPINT(IX))

      ELSE IF (WORD.EQ.'LB2') THEN
         LB2T=.TRUE.
!
! LJAT
!
      ELSE IF (WORD.EQ.'LJAT') THEN
         LJATT=.TRUE.
         CALL READF(ZSTAR)
         LJATTOC=2.423D0
         IF (NITEMS.EQ.3) CALL READF(LJATTOC)
      ELSE IF (WORD.EQ.'LOCALSAMPLE') THEN
         LOCALSAMPLET=.TRUE.
         IF (NITEMS.EQ.2) THEN
            CALL READF(ABTHRESH)
         ELSEIF (NITEMS.EQ.3) THEN
            CALL READF(ABTHRESH)
            CALL READF(ACTHRESH)
         ENDIF

!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'LJMF') THEN
         LJMFT=.TRUE.
!        CALL LJPARAMMF
         WRITE(MYUNIT,'(A)') 'LJMF not currently maintained'
         STOP

! start = an N-oligomer is constructed by relocating NREPEAT units and placing them
! at random distance R with Rmin <= R <= Rmax and angle alpha in the xy-plane
! transxy: rigid body translation only in the xy-plane
! rotz:  rigid body rotation only around the z-axis
! dihesc: only perturbation to the side chains
!
      ELSE IF (WORD.EQ.'MAKEOLIGO') THEN
         MAKEOLIGOT=.TRUE.
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'START') THEN
            MAKEOLIGOSTART=.TRUE.
            CALL READI(NFIXSEG)
            CALL READI(NREPEAT)
            ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
            REPATF(1:NREPEAT)=0
            REPATL(1:NREPEAT)=0
            REPPHIF(1:NREPEAT)=0.D0
            REPPHIL(1:NREPEAT)=0.D0
            DO J1=1,NREPEAT
               CALL READI(REPATF(J1))
               CALL READI(REPATL(J1))
               CALL READF(REPPHIF(J1))
               CALL READF(REPPHIL(J1))
            ENDDO
            CALL READF(PLACERMIN)
            CALL READF(PLACERMAX)
            TRANSXYT=.TRUE.
            ROTZT=.TRUE.
            CHBBT=.FALSE.
            CHSCT=.TRUE.
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'INITROT') THEN
            MAKEOLIGOSTART=.TRUE.
            INITROT=.TRUE.
            CALL READI(NFIXSEG)
            CALL READI(NREPEAT)
            IF (NREPEAT.GT.0) THEN
               ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
               REPATF(1:NREPEAT)=0
               REPATL(1:NREPEAT)=0
               REPPHIF(1:NREPEAT)=0.D0
               REPPHIL(1:NREPEAT)=0.D0
               DO J1=1,NREPEAT
                  CALL READI(REPATF(J1))
                  CALL READI(REPATL(J1))
                  CALL READF(REPPHIF(J1))
                  CALL READF(REPPHIL(J1))
               ENDDO
               CALL READF(PLACERMIN)
               CALL READF(PLACERMAX)
            ENDIF
            TRANSXYT=.TRUE.
            IF (NREPEAT.GT.0) THEN
               IF (NITEMS.GT.(6+NREPEAT*4)) THEN
                  CALL READA(UNSTRING)
                  IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') THEN
                     CHSCT=.TRUE.
                     CHBBT=.FALSE.
                  ENDIF
               ENDIF
            ELSE
               IF (NITEMS.GT.4) THEN
                  CALL READA(UNSTRING)
                  IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') THEN
                     CHSCT=.TRUE.
                     CHBBT=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'REFINE') THEN
            MAKEOLIGOSTART=.FALSE.
            IF (NITEMS.GT.2) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') THEN
               CHSCT=.TRUE.
               CHBBT=.FALSE.
            ENDIF
            IF (NITEMS.GT.3) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') THEN
               CHSCT=.TRUE.
               CHBBT=.FALSE.
            ENDIF
            IF (NITEMS.GT.4) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') THEN
               CHSCT=.TRUE.
               CHBBT=.FALSE.
            ENDIF
         ELSE
            WRITE(MYUNIT,'(A)') 'The first argument to MAKEOLIGO has to be START, INITROT or REFINE - quit.'
            STOP
         ENDIF

      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)

      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)
         MAXERISE_SET=.TRUE.

      ELSE IF (WORD.EQ.'MAXIT') THEN
         CALL READI(IX)
         MAXIT=IX
         IF (NITEMS.GT.2) THEN
            CALL READI(IX)
            MAXIT2=IX
         ENDIF

      ELSE IF (WORD.EQ.'MGGLUE') THEN
         MGGLUET=.TRUE.

      ELSE IF (WORD.EQ.'MINDENSITY') THEN
         MINDENSITYT=.TRUE.

! khs26> for use with free energy basin-hopping
! this keyword ensures that minimisation is repeated to a stricter convergence criterion if
! the zero and non-zero eigenvalues are not separated by MIN_ZERO_SEP
      ELSE IF (WORD.EQ.'MIN_ZERO_SEP') THEN
         IF (NITEMS.GT.1) THEN
             CALL READF(MIN_ZERO_SEP)
! Refers to frequencies, but checks eigenvalues.
             MIN_ZERO_SEP = MIN_ZERO_SEP**2
         END IF
         IF (NITEMS.GT.2) THEN
             CALL READI(MAX_ATTEMPTS)
         END IF

      ELSE IF (WORD.EQ.'MORSE') THEN
         MORSET=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(XX)
            RHO=XX
         ENDIF


!
! beh35> MLOWEST keyword
!

      ELSE IF (WORD.EQ.'MLOWEST') THEN
         MLOWEST=.TRUE.
         IF(TARGET) THEN
             PRINT *,"Cannot have both MLOWEST and TARGET"
             STOP
         END IF

         MTARGETS=NITEMS-1
         ALLOCATE(TARGETS(MTARGETS))

         !* Here's space to accommodate file input later

         DO J1=2,NITEMS
            CALL READF(XX)
            TARGETS(J1-1)=XX
         ENDDO

         ALLOCATE(OBJ(MTARGETS))

         CALL create_mtargets_array(TARGETS)

         DEALLOCATE(TARGETS)

!
!  MPI keyword
!
      ELSE IF (WORD.EQ.'MPI') THEN
         MPIT=.TRUE.
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,JUMPMOVE,JUMPINT,JDUMP,COORDS,NQ,
     @              JUMPTO,EPREV,COORDSO,VAT,VATO,SHELLMOVES,PTGROUP,NSURFMOVES,NCORE)
         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR),OSTEP(NPAR),
     @         BLOCK(NPAR),NT(NPAR),JUMPMOVE(NPAR),JUMPINT(NPAR),JDUMP(NPAR),COORDS(3*NATOMSALLOC,NPAR),NQ(NPAR),
     @         JUMPTO(NPAR),EPREV(NPAR),
     @         COORDSO(3*NATOMSALLOC,NPAR),VAT(NATOMSALLOC,NPAR),VATO(NATOMSALLOC,NPAR))
         ALLOCATE(SHELLMOVES(NPAR))
         ALLOCATE(PTGROUP(NPAR))
         ALLOCATE(NSURFMOVES(NPAR))
         ALLOCATE(NCORE(NPAR))
         DO JP=1,NPAR
            EPREV(JP)=1.0D100
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
            OSTEP(JP)=0.3D0
            BLOCK(JP)=0
            NT(JP)=0
            JUMPMOVE(JP)=.FALSE.
            JUMPINT(JP)=100
            JDUMP(JP)=.FALSE.
            SHELLMOVES(JP)=.FALSE.
            PTGROUP(JP)='    '
            NSURFMOVES(JP)=0
            NCORE(JP)=0
         ENDDO

      ELSE IF (WORD.EQ.'MSORIG') THEN
         MSORIGT=.TRUE.

      ELSE IF (WORD.EQ.'MSTRANS') THEN
         MSTRANST=.TRUE.

      ELSE IF (WORD.EQ.'MULLERBROWN') THEN
         MULLERBROWNT=.TRUE.

      ELSE IF (WORD.EQ.'MULTIPLICITY') THEN
         CALL READI(XMUL)
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'MYSD') THEN
         MYSDT=.TRUE.
         CALL READF(SDTOL)

      ELSE IF (WORD.EQ.'NATB') THEN
         NATBT=.TRUE.

      ELSE IF (WORD.EQ.'NEON') THEN
         NEON=.TRUE.
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'NEWCONF') THEN
         NEWCONFT=.TRUE.
         CALL READI(NEWCONFST)
         CALL READF(NCWALL)

      ELSE IF (WORD.EQ.'NEWJUMP') THEN
         NEWJUMP=.TRUE.
         IF (NITEMS.GT.1) CALL READF(PNEWJUMP)
!
!  Reseed runs if the energy does not decrease within NRELAX mc steps.
!  NHSRESTART defines the number of hard sphere moves used to produce the new starting
!  configuration. If NHSRESTART=0 then the geometry is changed using RESEED.
!
      ELSE IF (WORD.EQ.'NEWMOVES') THEN
         NEWMOVEST=.TRUE.
      ELSE IF (WORD.EQ.'NEWRESTART') THEN
         NEWRESTART=.TRUE.
         NEWRESTART_MD = .FALSE.                   ! lb415
         NEWRES_TEMP = 0.0D0
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
         IF (NITEMS.GT.3) CALL READA(WORD2)        ! lb415
         IF (WORD2.EQ.'MD') NEWRESTART_MD = .TRUE. ! lb415
         IF (NITEMS.GT.4) THEN                     ! lb415
            CALL READF(NEWRES_TEMP)
         ELSE
            NEWRES_TEMP = 1000.0D0
            WRITE(MYUNIT,'(A)') 'keyword> WARNING - temperature unspecified for NEWRESTART. Default for high T MD is 1000K'
         ENDIF                                     ! lb415
         IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
         IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMSALLOC,MAXSAVE))

!      ELSE IF (WORD.EQ.'NMAX') THEN
!         CALL READF(NMAX)
!         WRITE(MYUNIT,'(A,F14.10)') 'NMAX=  ',NMAX
!
!      ELSE IF (WORD.EQ.'NMIN') THEN
!         CALL READF(NMIN)
!         WRITE(MYUNIT,'(A,F14.10)') 'NMIN=  ',NMIN
 
      ELSE IF (WORD.EQ.'NOCHIRALCHECKS') THEN
         CHECKCHIRALITY=.FALSE.

      ELSE IF (WORD.EQ.'UACHIRAL') THEN
         UACHIRAL=.TRUE.

      ELSE IF (WORD.EQ.'NOCISTRANS') THEN
         IF (NITEMS.GT.1) THEN
            CALL READF(MINOMEGA)
            CIS_TRANS_TOL = 180.0D0 - MINOMEGA
         END IF

      ELSE IF (WORD.EQ.'NOCISTRANSDNA') THEN
         NOCISTRANSDNA=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(MINOMEGA)
            CIS_TRANS_TOL = 180.0D0 - MINOMEGA
         END IF

      ELSE IF (WORD.EQ.'NOCISTRANSRNA') THEN
         NOCISTRANSRNA=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(MINOMEGA)
            CIS_TRANS_TOL = 180.0D0 - MINOMEGA
         END IF
      
      ELSE IF (WORD.EQ.'NOFREEZE') THEN
         FREEZECORE=.FALSE.
!
!  Turn off inversion in distance optimisation for minpermdist.
!
      ELSE IF (WORD.EQ.'NOINVERSION') THEN
         NOINVERSION=.TRUE.

      ELSE IF (WORD.EQ.'NORESET') THEN
         NORESET=.TRUE.

! hk286 > redefine rigid body sites every NRelaxRigid step   !
      ELSE IF (WORD.EQ.'NRELAXRIGID') THEN
         RELAXRIGIDT = .TRUE.
         CALL READI(NRELAXRIGIDR)
         CALL READI(NRELAXRIGIDA)
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'OH') THEN
         FIELDT=.TRUE.
         OHT=.TRUE.
         CALL READF(XX)
         FOH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(EXPFAC)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF
!
!  Specify Oh supercell to allow box symmetries in permutational alignment.
!
      ELSE IF (WORD.EQ.'OHCELL') THEN
         OHCELLT=.TRUE.
         WRITE(MYUNIT,'(A)') 'Octahedral supercell specfied'
!
!  Specify 1D XY APBC potential
!
      ELSE IF (WORD.EQ.'ONEDAPBC') THEN
         ONEDAPBCT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
         IF (MOD(NONEDAPBC,3).NE.0) THEN
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
            STOP
         ENDIF

!
!  Specify 1D XY PBC potential
!
      ELSE IF (WORD.EQ.'ONEDPBC') THEN
         ONEDPBCT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
         IF (MOD(NONEDAPBC,3).NE.0) THEN
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
            STOP
         ENDIF

! hk286
      ELSE IF (WORD.EQ.'OPTIMISEROTATION') THEN
         RIGIDOPTIMROTAT = .TRUE.
         CALL READF(OPTIMROTAVALUES(1))
         CALL READF(OPTIMROTAVALUES(2))
         CALL READF(OPTIMROTAVALUES(3))

!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'OTP') THEN
         OTPT=.TRUE.
         RIGID=.TRUE.
!        CALL OTPPARAMMF

!js850> OVERLAP: calculate the overlap with the initial structre every
!OVERLAP_FRQ mc steps.  The overlap is defined by dividing space into cubic
!cells of size OVERLAP_DR, and counting the number of occupied cells in common
!between the two structures.  overlap of 1 means the strucures are the same.
!overlap of 0 means they are significantly different structurally.  This is not
!a robust comparison if there is translational or rotational symmetry, so it is
!really only useful if some particles are frozen, or the symmetry is broken in
!some other manner.  There are other definitions of overlap which can fix this
!problem
!js850> update: I now calculate many different definitions of overlap.  Still
!not sure which is the best order parameter for my system
      ELSE IF (WORD.EQ.'OVERLAP') THEN
         OVERLAPK=.TRUE.
         CALL READF(OVERLAP_DR)
         CALL READI(OVERLAP_FRQ)
         IF (NITEMS.GT.3) THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'I') OVERLAP_IMPORT=.TRUE.
         ENDIF

!js850> move only one atom, randomly chosen, during takestep routine.
!Only implemented in bspt.F
      ELSE IF (WORD.EQ.'ONEATOM_TAKESTEP') THEN
         ONE_ATOM_TAKESTEP=.TRUE.

      ELSE IF (WORD.EQ.'P46') THEN
         P46=.TRUE.
         BLNT=.TRUE.

      ELSE IF (WORD.EQ.'PACHECO') THEN
         PACHECO=.TRUE.

      ELSE IF (WORD.EQ.'PAH') THEN
         PAHT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=36
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'PAIRDIST') THEN
! csw34> PAIRDIST allows the monitoring of the distance between pairs of
! atoms during a BH run. The atom pairs are read from either as
! arguements to the PAIRDIST keyword, or from the pairdist file
         NPAIRS=0
         IF (NITEMS.GT.1) THEN
! If arguements are specified, assume NITEMS/2 pairs on the line 
            PAIRDISTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            NPAIRS=(NITEMS-1)/2
            ALLOCATE(PAIRDIST(NPAIRS,2))
            DO J1=1,NPAIRS
               CALL READI(PAIRDIST(J1,1))
               CALL READI(PAIRDIST(J1,2))
               IF (PAIRDIST(J1,1).GT.NATOMSALLOC) THEN
                  WRITE(MYUNIT,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  STOP
               ELSEIF (PAIRDIST(J1,2).GT.NATOMSALLOC) THEN
                  WRITE(MYUNIT,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  STOP
               ENDIF
            ENDDO
         ELSE
! If there are no atoms specified on the PAIRDIST line, assume reading them from the 'pairdist' file
! First step - check the pairdist file is present
            YESNO=.FALSE.
            INQUIRE(FILE='pairdist',EXIST=YESNO)
            IF (YESNO) THEN
               PAIRDISTT=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword> ERROR: pairdist input file missing for PAIRDIST'
               CALL FLUSH(MYUNIT)
               STOP
            ENDIF
! Determine NPAIRS to allow allocation of the PAIRDIST array 
            OPEN(UNIT=222,FILE='pairdist',status='old')
            DO
               READ(222,*,IOSTAT=iostatus) CHECK1
               IF (iostatus<0) THEN
                  CLOSE(222)
                  EXIT
               ELSE 
                  NPAIRS=NPAIRS+1
               ENDIF
            END DO        
            CLOSE(222)
! Allocate the PAIRDIST array and read the pairs in
            ALLOCATE(PAIRDIST(NPAIRS,2))
            OPEN(UNIT=222,FILE='pairdist',status='old')
            DO J1=1,NPAIRS
               READ(222,*) PAIRDIST(J1,1),PAIRDIST(J1,2)
            ENDDO
            CLOSE(222)
         ENDIF
! Print list of pairs to GMIN output for checking
         WRITE(MYUNIT,'(A)') ' keyword> Atom pair distances to monitor:'
         DO J1=1,NPAIRS
            WRITE(MYUNIT,*) PAIRDIST(J1,:)
         ENDDO

!
!  PARALLEL must come before STEP and ACCRAT
!  This keyword is for the serial parallel implementation - now obsolete.
!
      ELSE IF (WORD.EQ.'PARALLEL') THEN
         PARALLELT=.TRUE.
         CALL READI(NPAR)
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,JUMPMOVE,JUMPINT,JDUMP,COORDS,NQ,
     @              JUMPTO,EPREV,COORDSO,VAT,VATO,SHELLMOVES,PTGROUP,NSURFMOVES,NCORE) 
         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR),OSTEP(NPAR), 
     @         BLOCK(NPAR),NT(NPAR),JUMPMOVE(NPAR),JUMPINT(NPAR),JDUMP(NPAR),NQ(NPAR),JUMPTO(NPAR),
     &         COORDS(3*NATOMSALLOC,NPAR),
     @         COORDSO(3*NATOMSALLOC,NPAR),VAT(NATOMSALLOC,NPAR),VATO(NATOMSALLOC,NPAR),EPREV(NPAR),SHELLMOVES(NPAR),PTGROUP(NPAR),
     @         NSURFMOVES(NPAR),NCORE(NPAR),REPLOW(NPAR))
         NATOMS=NATOMS/NPAR
         DO JP=1,NPAR
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
            OSTEP(JP)=0.3D0
            BLOCK(JP)=0
            NT(JP)=0
            JUMPMOVE(JP)=.FALSE.
            JUMPINT(JP)=100
            JDUMP(JP)=.FALSE.
            SHELLMOVES(JP)=.FALSE.
            PTGROUP(JP)='    '
            NSURFMOVES(JP)=0
            NCORE(JP)=0
         ENDDO
         IF (NITEMS.GT.1) THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'GMIN') REPMATCHT=.TRUE.
         ENDIF

      ELSE IF (WORD.EQ.'PBGLUE') THEN
         PBGLUET=.TRUE.

      ELSE IF (WORD.EQ.'PERCOLATE') THEN
         PERCOLATET=.TRUE.
         PERCACCEPTED=.FALSE.
         IF (NITEMS.GT.1) CALL READF(PERCCUT)
         IF (NITEMS.GT.2) CALL READF(K_COMP)
         IF (NITEMS.GT.3) CALL READF(GUIDECUT)

      ELSE IF (WORD.EQ.'PERIODIC') THEN
         PERIODIC=.TRUE.
         CALL READF(XX)
         BOXLX=XX
         BOX3D(1) = XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            BOXLY=XX
            BOX3D(2) = XX
            IF (NITEMS.GT.3) THEN
               CALL READF(XX)
               BOXLZ=XX
               BOX3D(3) = XX
            ENDIF
         ELSE
            BOXLY=BOXLX
            BOX3D(2) = BOX3D(1)
            BOXLZ=BOXLX
            BOX3D(3) = BOX3D(1)
         ENDIF
         PI = 4.D0*DATAN(1.D0)
         PALPHA=PI/2.D0
         PBETA=PALPHA
         PGAMMA=PALPHA
         IF (NITEMS.GT.4) THEN
            CALL READF(XX)
            PGAMMA=PI*XX/180.D0
         ENDIF
         IF (NITEMS.GT.6) THEN
            PALPHA=PGAMMA
            CALL READF(XX)
            PBETA=PI*XX/180.D0
            CALL READF(XX)
            PGAMMA=PI*XX/180.D0
         ENDIF
! lower triangular lattice matrix
!         LAT(3,3)=BOXLZ
!         LAT(2,2)=BOXLY*SIN(PALPHA)
!         LAT(3,2)=BOXLY*COS(PALPHA)
!         LAT(1,1)=BOXLX/SIN(PALPHA)*SQRT(1.D0-COS(PALPHA)**2-COS(PBETA)**2-COS(PGAMMA)**2+2D0*COS(PALPHA)*COS(PBETA)*COS(PGAMMA))
!         LAT(2,1)=BOXLX*(COS(PGAMMA)-COS(PBETA)*SIN(PALPHA))/SIN(PALPHA)
!         LAT(3,1)=BOXLX*COS(PBETA)

! upper triangular lattice matrix

         IF (MWFILMT) BOXLZ=BOXLZ*SIN(PGAMMA)/SQRT(1.D0-COS(PALPHA)*COS(PALPHA)-COS(PBETA)*COS(PBETA)
     1     -COS(PGAMMA)*COS(PGAMMA)+2.D0*COS(PALPHA)*COS(PBETA)*COS(PGAMMA))

         LAT(1,1)=BOXLX
         LAT(1,2)=BOXLY*COS(PGAMMA)
         LAT(1,3)=BOXLZ*COS(PBETA)

         XX=ACOS((COS(PBETA)*COS(PGAMMA)-COS(PALPHA))/SIN(PBETA)/SIN(PGAMMA))

         LAT(2,1)=0.D0
         LAT(2,2)=BOXLY*SIN(PGAMMA)
         LAT(2,3)=-BOXLZ*sin(PBETA)*COS(XX)

         XX=SQRT(1.D0-COS(PALPHA)*COS(PALPHA)-COS(PBETA)*COS(PBETA)
     1     -COS(PGAMMA)*COS(PGAMMA)+2.D0*COS(PALPHA)*COS(PBETA)*COS(PGAMMA))

         XX=SIN(PGAMMA)/BOXLZ/XX

         LAT(3,1)=0.D0
         LAT(3,2)=0.D0
         LAT(3,3)=1.D0/XX
         WRITE(MYUNIT,*)LAT(:,1)
         WRITE(MYUNIT,*)LAT(:,2)
         WRITE(MYUNIT,*)LAT(:,3)

!
! If permdist is set then distance calculations are performed with minpermdist instead
! of newmindist in procedures such as AVOID and CSM. This keyword is now independent
! from PERMOPT
!
      ELSE IF (WORD.EQ.'PERMDIST'.AND.PERMOPT) THEN
         WRITE(MYUNIT,'(A)') 'keywords> PERMDIST has already been set by PERMOPT keyword'
         IF (NITEMS.GT.1) CALL READF(ORBITTOL)
         WRITE(MYUNIT,'(A,F15.5)') ' keyword> Distance tolerance for distinguising atoms in the same orbit=',ORBITTOL

      ELSE IF (WORD.EQ.'PERMDIST'.AND.(.NOT.PERMOPT)) THEN
         PERMDIST=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ORBITTOL)
         WRITE(MYUNIT,'(A,F15.5)') ' keyword> Distance tolerance for distinguising atoms in the same orbit=',ORBITTOL
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         IF (.NOT.ALLOCATED(NPERMSIZE)) THEN
            ALLOCATE(NPERMSIZE(NATOMSALLOC),PERMGROUP(NATOMSALLOC),NSETS(NATOMSALLOC),SETS(NATOMSALLOC,70))
         ENDIF
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
!
!  Sanity checks!
!
               IF (NSETS(J1).GT.70) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 70'
                  STOP
               ENDIF
!              IF (NDUMMY+NPERMSIZE(J1).GT.NATOMSALLOC) THEN
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMSALLOC) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check!
!  
!           DO J1=1,NDUMMY
!              DO J2=J1+1,NDUMMY
!                 IF (PERMGROUP(J2).EQ.PERMGROUP(J1)) THEN
!                    PRINT '(2(A,I8))','keyword> ERROR - atom ',PERMGROUP(J1),' appears more than once'
!                    STOP
!                 ENDIF
!              ENDDO
!           ENDDO
            CLOSE(1)
!
!  And yet another!
!  
            IF (NFREEZE.GT.0) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
         ELSE
            NSETS(1:NATOMSALLOC)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMSALLOC ! all atoms can be permuted - default
            DO J1=1,NATOMSALLOC
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            WRITE(MYUNIT,'(A,3(I6,A))') ' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(MYUNIT,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1)
            IF (NSETS(J1).GT.0) THEN
               WRITE(MYUNIT,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(MYUNIT,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' / '
                  ENDDO
                  IF (J2.LT.NSETS(J1)) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' ; '
               ENDDO
            ENDIF
            WRITE(MYUNIT,'(A)') ' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO
!
!  This keyword is for optimising the distance between permutational isomers.
!  PERMINVOPT allows inversions.
!
      ELSE IF ((WORD.EQ.'PERMOPT').OR.(WORD.EQ.'PERMINVOPT')) THEN
         IF (WORD.EQ.'PERMOPT') PERMOPT=.TRUE.
         IF (WORD.EQ.'PERMINVOPT') PERMINVOPT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ORBITTOL)
         WRITE(MYUNIT,'(A,F15.5)') ' keyword> Distance tolerance for distinguising atoms in the same orbit=',ORBITTOL

         STEP(1)=6.0D0
         TEMP(1)=1.0D0
         WRITE(MYUNIT,'(A,F15.5,A)') 'keywords> Setting default step size to ',STEP(1),' for PERMOPT rotations'
         WRITE(MYUNIT,'(A,F15.5,A)') 'keywords> Setting default temperature to ',TEMP(1),' for PERMOPT'
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         IF (.NOT.ALLOCATED(PERMGROUP)) ALLOCATE(PERMGROUP(NATOMSALLOC))
         IF (.NOT.ALLOCATED(NSETS)) ALLOCATE(NSETS(NATOMSALLOC))
         IF (.NOT.ALLOCATED(SETS)) ALLOCATE(SETS(NATOMSALLOC,70))
         IF (.NOT.ALLOCATED(NPERMSIZE)) ALLOCATE(NPERMSIZE(NATOMSALLOC))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
!
!  Sanity checks!
!
               IF (NSETS(J1).GT.70) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 70'
                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMSALLOC) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),
!    &                                 ((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))

               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check!
!  
!           DO J1=1,NDUMMY
!              DO J2=J1+1,NDUMMY
!                 IF (PERMGROUP(J2).EQ.PERMGROUP(J1)) THEN
!                    PRINT '(2(A,I8))','keyword> ERROR - atom ',PERMGROUP(J1),' appears more than once'
!                    STOP
!                 ENDIF
!              ENDDO
!           ENDDO
            CLOSE(1)
!
!  And yet another!
!  
            IF (NFREEZE.GT.0) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
         ELSE
            NSETS(1:NATOMSALLOC)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMSALLOC ! all atoms can be permuted - default
            DO J1=1,NATOMSALLOC
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            WRITE(MYUNIT,'(A,3(I6,A))') ' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(MYUNIT,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1)
            IF (NSETS(J1).GT.0) THEN
               WRITE(MYUNIT,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(MYUNIT,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' / '
                  ENDDO
                  IF (J2.LT.NSETS(J1)) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' ; '
               ENDDO
            ENDIF
            WRITE(MYUNIT,'(A)') ' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO
!
! Initial step size for geometry perturbations used in RESERVOIR calculations.
!
      ELSE IF (WORD.EQ.'PERTSTEP') THEN
         CALL READF(PERTSTEP)

      ELSE IF (WORD.EQ.'PLUS') THEN
         PLUS=.TRUE.

!      ELSE IF (WORD.EQ.'PMAX') THEN
!         CALL READF(PMAX)
!         WRITE(MYUNIT,'(A,F14.10)') 'PMAX=  ',PMAX
!
!      ELSE IF (WORD.EQ.'PMIN') THEN
!         CALL READF(PMIN)
!         WRITE(MYUNIT,'(A,F14.10)') 'PMIN=  ',PMIN
!
!  POWER provides a means to set the initial premultiplication factor for the
!  gradient in MYLINMIN
!
      ELSE IF (WORD.EQ.'POWER') THEN
         CALL READI(IX)
         MYPOWER=IX
!
!  Purify the geometry in mylbfgs to preserve icosahedral (I)
!  symmetry.
!
      ELSE IF (WORD.EQ.'PROJI') THEN
         PROJIT=.TRUE.
!
!  Purify the geometry in mylbfgs to preserve icosahedral (Ih)
!  symmetry.
!
      ELSE IF (WORD.EQ.'PROJIH') THEN
         PROJIHT=.TRUE.
!
!  Frequency of printing in lbfgs to reduce size of output files
!  NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'PRTFRQ') THEN
         CALL READI(PRTFRQ)
!
!  Plain parallel tempering. Same as BSPT but without quenches.
!
      ELSE IF (WORD.EQ.'PTMC') THEN
         PTMC=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(PTEMIN)
         CALL READF(PTEMAX)
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         CALL READF(NEQUIL)
         CALL READF(PTSTEPS)
         NQUENCH=0.0D0
         CALL READI(NENRPER)
         CALL READI(HBINS)
         QUENCHFRQ=1 ! must be set to avoid division by zero in bspt.F
         IF (NITEMS.GT.12) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'R') PTRANDOM=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'I') PTINTERVAL=.TRUE.
         ELSE
            PTRANDOM=.TRUE.
         ENDIF
         IF (NITEMS.GT.13) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SIN') PTSINGLE=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SET') PTSETS=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'IND') PTEX_INDEP=.TRUE.
            ! ** Note: if use use PTEX_INDEP then you must use PTINTERVAL.  No
            ! physical reason for this, just a result of the hacky way
            ! tryexchange is currently coded ***
         ELSE
            PTSINGLE=.TRUE.
         ENDIF
         IF (NITEMS.GT.14) THEN
            CALL READI(NDUMMY)
            CALL SDPRND(NDUMMY)
         ENDIF
         EXCHINT=INT(1.D0/EXCHPROB)
! jdf43> the following would enforce the same average number of pair
! exchanges for the single & sets exchange schemes - ss2029 prefers this
! not to be the case.
!
!         IF (PTSETS) THEN
!            EXCHINT = INT(EXCHINT*(NPAR-1)/2)
!            EXCHPROB= EXCHPROB*2.D0/(1.D0*(NPAR-1))
!         ENDIF
!         CALL READI(QUENCHFRQ)

! js850 PTEX_INDEP: use alternate tryexchange routine where all the replicas
! are redistributed independently according to the appropriate distribution.
! ** must use with PTINTERVAL.  No physical reason for this, just a result of
! the hacky way tryexchange is currently coded ***
      ELSE IF (WORD.EQ.'PTEX_INDEP') THEN
         PTEX_INDEP = .TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(PTEX_INDEP_NSTEPS)
         ELSE
            PTEX_INDEP_NSTEPS = -1 !use the default
         ENDIF

! js850 PTMCDUMPSTRUCT: dump the coords every PTMCDS_FRQ MC steps.  Dump
! unquenched coords
      ELSE IF (WORD.EQ.'PTMCDUMPSTRUCT') THEN
         PTMCDUMPSTRUCT=.TRUE.
         CALL READI(PTMCDS_FRQ)

! ss2029 PTMCDUMPENER: for every replica dump Markov energy PTMCDUMPENERFRQ MC steps
      ELSE IF (WORD.EQ.'PTMCDUMPENER') THEN
         PTMCDUMPENERT=.TRUE.
         CALL READI(PTMCDUMPENERFRQ)

! js850> PTREADTEMPS: read the temperatures for the parallel tempering replicas
! from file temperatures.init
      ELSE IF (WORD.EQ.'PTREADTEMPS') THEN
         PTREADTEMPS=.TRUE.


!
!  Keyword for applied static force.
!
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
!
! Request calculation of structural order parameters Q on the fly 
! NOT DOCUMENTED.
!
      ELSE IF (WORD.EQ.'CALCQ') THEN
         CALCQT=.TRUE.
!
! Distance cut-off for Coulomb interactions in AMBER (the PNM hand-coded version).
!
!      ELSE IF (WORD.EQ.'QCUTOFF') THEN
!         AMCUT=.TRUE.
!         CALL READF(REALQCUTOFF)
!         QCUTOFF=1.1D0*REALQCUTOFF
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'QDTEST2') THEN
         QD2T=.TRUE.
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'QDTEST') THEN
         QDT=.TRUE.

      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
         CALL READF(CQMAX)

      ELSE IF (WORD.EQ.'QUAD') THEN
         QUADT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=5
         ALLOCATE(SITE(NRBSITES,3))
!
! Collect data from quenches for hopeful conversion into relative density of states.
! Requires saved lowestdirect and firstfit files to be present.
! NOT DOCUMENTED.
!
      ELSE IF (WORD.EQ.'QUENCHDOS') THEN
         QUENCHDOS=.TRUE.
         IF (NITEMS.GT.1) CALL READI(QDLIMIT)

      ELSE IF (WORD.EQ.'RADIUS') THEN
         CALL READF(XX)
         RADIUS=XX
!
!  integer seed for random number generator.
!
      ELSE IF (WORD.EQ.'RANSEED') THEN
         RANSEEDT=.TRUE.
         CALL READI(NDUMMY)
         CALL SDPRND(NDUMMY+MYNODE)
         CALL SDPRND_UNIVERSAL(NDUMMY+NPAR)
         WRITE(MYUNIT,'(A,I8)') 'keywords> Random seed = ',NDUMMY+MYNODE 

!
!  RBSYM defines the internal symmetry operations for each sort of rigid body
!  coded via RBAAT.
!
      ELSE IF (WORD.EQ.'RBSYM') THEN
         RBSYMT=.TRUE.
         INQUIRE(FILE='rbsymops',EXIST=RBSYMTEST)
         IF (RBSYMTEST) THEN
            OPEN(UNIT=111,FILE='rbsymops',STATUS='OLD')
            READ(111,*) NRBGROUP
            ALLOCATE(RBOPS(4,NRBGROUP))
            READ(111,*) ((RBOPS(J1,J2),J1=1,4),J2=1,NRBGROUP)
            WRITE(MYUNIT,'(A,I6)')' keywords> number of symmetry operations for rigid body=',NRBGROUP
            DO J1=1,NRBGROUP
               WRITE(MYUNIT, '(A,I6)')' keywords> rigid-body symmetry operation', J1
               RBOPS(4,J1) = RBOPS(4,J1)*ATAN(1.D0)/45.D0
               WRITE(MYUNIT,'(3F20.10)')RBOPS(1:4,J1)
            ENDDO
            CLOSE(111)
         ELSE
            WRITE(MYUNIT,'(A)')' keywords> ERROR *** missing file rbsymops'
            STOP
         ENDIF
!
!  If READMASS is specified we read the masses from file masses.
!
      ELSE IF (WORD.EQ.'READMASS') THEN
         READMASST=.TRUE.

! If SPECMASS is specified we read the mass for each species.
      ELSE IF(WORD.EQ.'SPECMASS') THEN
         SPECMASST=.TRUE.
         IF(ALLOCATED(SPECMASS)) DEALLOCATE(SPECMASS)
         ALLOCATE(SPECMASS(1:NSPECIES(0)))
         IF(NITEMS .GT. NSPECIES(0)) THEN
            DO J1=1,NSPECIES(0)
               CALL READF(SPECMASS(J1))
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)')' keywords> ERROR *** not enough arguments for SPECMASS'
            STOP
         ENDIF
!
! Low temperature replica will use min.data and points.min entries.
!
      ELSE IF (WORD.EQ.'RESERVOIR') THEN
         RESERVOIRT=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READI(USERES)
            CALL READF(RES_PSWAP)
         ENDIF
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','keyword> ERROR *** min.data must exist for a RESERVOIR run'
            STOP
         ENDIF
         INQUIRE(FILE='points.min',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','keyword> ERROR *** points.min must exist for a RESERVOIR run'
            STOP
         ENDIF
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE='min.data',STATUS='OLD')
         NRESMIN=0
         DO
            READ(LUNIT,*,END=31) DUMMY
            NRESMIN=NRESMIN+1
         ENDDO
31       WRITE(MYUNIT,'(A,I6,A)') 'keyword> There are ',NRESMIN,' minima in the min.data file of RESERVOIR'
         ALLOCATE(EMIN(NRESMIN),FVIBMIN(NRESMIN),PFMIN(NRESMIN),IXMIN(NRESMIN),IYMIN(NRESMIN),IZMIN(NRESMIN),HORDERMIN(NRESMIN))
         ALLOCATE(RESPOINTS(3*NATOMSALLOC,NRESMIN))
         REWIND(LUNIT)
         DO J1=1,NRESMIN
            READ(LUNIT,*) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
         ENDDO
         CLOSE(LUNIT)
         LUNIT=GETUNIT()
!         OPEN(LUNIT,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMSALLOC)
         OPEN(LUNIT,FILE='points.min',STATUS='OLD')
         DO J1=1,NRESMIN
            DO J2=1,3*NATOMS,3
               READ(LUNIT,'(3G25.15)') RESPOINTS(J2,J1),RESPOINTS(J2+1,J1),RESPOINTS(J2+2,J1)
               !READ(LUNIT,'(3G25.15)') RESPOINTS(3*J2-2,J1),RESPOINTS(3*J2-1,J1),RESPOINTS(3*J2,J1)
               !WRITE(*,'(A,2I3,G10.10)') "read points.min file", J1, J2, RESPOINTS(J2,J1)
            ENDDO
         ENDDO
!         DO J2=1,NRESMIN
!            DO J1=1,NATOMS
!               WRITE(MYUNIT,'(A,3F20.6,2I)') 'setup> RESPOINTS ', RESPOINTS(3*(J1-1)+1,J2), 
!     &          RESPOINTS(3*(J1-1)+2,J2),RESPOINTS(3*(J1-1)+3,J2),3*(J1-1)+3,J2
!            ENDDO
!         ENDDO
         CLOSE(LUNIT)
         IF (DEBUG) WRITE(MYUNIT,'(A,I6,A)') 'setup> points for ',NRESMIN,' minima read from file points.min'

         ! AB > read eigenvalues / eigenvectors for each minimum from file
         ! vector.dump
         ALLOCATE(HESSEIGVC(3*NATOMSALLOC,3*NATOMSALLOC,NRESMIN))
         ALLOCATE(HESSEIGVA(3*NATOMSALLOC,NRESMIN))
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE='vector.dump',STATUS='OLD')
         VECNUM=0
         RESCOUNT=1
         DO J1=1,(3*NATOMSALLOC-6)*(NATOMSALLOC+1)*NRESMIN
         !DO J1=1,(3*NATOMSALLOC)*(NATOMSALLOC+1)*NRESMIN
            VECLINE=MOD(J1-1,NATOMSALLOC+1)
            IF (VECLINE.EQ.0) THEN
               VECNUM=VECNUM+1
               !IF (VECNUM.GT.(3*NATOMSALLOC)) THEN
               IF (VECNUM.GT.(3*NATOMSALLOC-6)) THEN
                  !comment below loop if vector.dump stores all 3N eivenvectors
                  DO J3=3*NATOMSALLOC-5,3*NATOMSALLOC
                     HESSEIGVC(1:3*NATOMSALLOC,J3,RESCOUNT) = 0.0
                     HESSEIGVA(J3,RESCOUNT) = 0.0
                  ENDDO
                  VECNUM=1
                  RESCOUNT=RESCOUNT+1
               ENDIF
            ENDIF
            IF (VECLINE.EQ.0) THEN
                  READ(LUNIT,*) HESSEIGVA(VECNUM,RESCOUNT)
                  IF ((VECNUM.GT.(3*NATOMSALLOC-6)).AND.(VECNUM.LE.(3*NATOMSALLOC))) HESSEIGVA(VECNUM,RESCOUNT) = 0.0
            ELSE
                READ(LUNIT,*) (HESSEIGVC(J2,VECNUM,RESCOUNT),J2=3*(VECLINE-1)+1,3*(VECLINE-1)+3)
!                WRITE(MYUNIT,'(A,3I,2F)') "eigvc ",3*(VECLINE-1)+1,VECNUM,RESCOUNT,
!     &           HESSEIGVC(3*(VECLINE-1)+1,VECNUM,RESCOUNT),HESSEIGVA(VECNUM,RESCOUNT)
            ENDIF
         ENDDO

         !comment below loop if vector.dump stores all 3N eivenvectors
         DO RESCOUNT=1,NRESMIN
            DO J3=3*NATOMSALLOC-5,3*NATOMSALLOC
               HESSEIGVC(1:3*NATOMSALLOC,J3,RESCOUNT) = 0.0
               HESSEIGVA(J3,RESCOUNT) = 0.0
            ENDDO
         ENDDO

         CLOSE(LUNIT)
         WRITE(MYUNIT,'(A)') 'setup> Eigenvalues and Eigenvectors read in'

         ! set reservoir well probabilities:
         DUMMY=0.0
         ALLOCATE(PW(NRESMIN))
         DO J2=1,NRESMIN
            PW(J2)= (1./HORDERMIN(J2)) * EXP(-BETA_RES*(EMIN(J2)-EMIN(1)))
            DO J1=1,3*NATOMS-6
               !PW(J2)=PW(J2)*(1./HESSEIGVA(J1,J2)
               PW(J2)=PW(J2)*DSQRT(1./HESSEIGVA(J1,J2))
            ENDDO
            DUMMY=DUMMY+PW(J2)
         ENDDO
         ! Normalize Pw
         DO J2=1,NRESMIN
            PW(J2)=PW(J2)*1.0/DUMMY
            WRITE(MYUNIT,'(A,I4,G20.10,A,G20.10)') 'setup> reservoir well probabilities: ',J2,PW(J2),"at temp",1./(BETA_RES)
         ENDDO

!      ELSE IF (WORD.EQ.'RCUTOFF') THEN
!         AMCUT=.TRUE.
!         CALL READF(REALRCUTOFF)
!         RCUTOFF=1.1D0*REALRCUTOFF
!
!  Read data for previous geometries that lead to reseeding, which
!  probably approximate MSB bottoms.
!  NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'READMSB') THEN
         INQUIRE(FILE='MSBdata',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(MYUNIT,'(A)') 'ERROR - READMSB specified, but no MSBdata data file found'
            STOP
         ELSE
            IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMSALLOC,MAXSAVE))
            IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
            OPEN(UNIT=34,FILE='MSBdata',STATUS='OLD')
57          READ(34,*,END=56) DUMMY
            NMSBSAVE=NMSBSAVE+1
            MSBE(NMSBSAVE)=DUMMY
            READ(34,*) (MSBCOORDS(J1,NMSBSAVE),J1=1,3*NATOMSALLOC)
            IF (NMSBSAVE.LT.MAXSAVE) GOTO 57
56          WRITE(MYUNIT,'(A,I6,A)') 
     1         'Energies and coordinates read for ',NMSBSAVE,' previous structures from MSBdata'
            CLOSE(34)
         ENDIF

! hk286 > use atom coords during final quench   !
      ELSE IF (WORD.EQ.'RELAXFINALQUENCH') THEN
         RELAXFQ = .TRUE.

!
!  Renormalisation attempt
!
!  TRENORM is the temperature for the Metropolis accept/reject comparison
!          of lowest energies calculated over NRENORM steps having moved
!          XMOVERENORM atoms. NRENORM is dynamically adjusted with a
!          minimum value equal to half the original NRENORM.
!  NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'RENORM') THEN
         RENORM=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRENORM)
         IF (NITEMS.GT.2) CALL READF(XMOVERENORM)
         IF (NITEMS.GT.3) CALL READF(TRENORM)
         IF (NITEMS.GT.4) CALL READI(NRENSTUCK)

      ELSE IF (WORD.EQ.'RESIZE') THEN
         RESIZET=.TRUE.
         CALL READF(XX)
         RESIZE=XX
!
!  Reseed runs if a step is not accepted in twice the relaxation time,
!  defined in terms of a number of mc steps NRELAX. NHSRESTART defines
!  the number of hard sphere moves used to produce the new starting
!  configuration. 
!
      ELSE IF (WORD.EQ.'RESTART') THEN
         RESTART=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
!
!  Restore the state of a previous GMIN run from dumpfile.
!
      ELSE IF (WORD.EQ.'RESTORE') THEN
         RESTORET=.TRUE.
         CALL READA(DUMPFILE)
         IF (NITEMS.GT.2) THEN
            CALL READA(INTEDUMPFILE)
            INTERESTORE=.TRUE.
         ENDIF

! hk286
      ELSE IF (WORD.EQ.'RESTRAINL') THEN

         RESTRAINLT   = .TRUE.
         CALL READF(RESTRAINLK)
         CALL READRESTRAINL()

!
! js850> reject moves if mobile particles are more than 
! RESTRICTREGIONRADIUS from point (RESTRICTREGIONX0, RESTRICTREGIONY0, RESTRICTREGIONZ0).
! Check mobility by looking at FROZEN, HARMONICFLIST, and DONTMOVE
!
      ELSE IF (WORD.EQ.'RESTRICTREGION') THEN
            RESTRICTREGION=.TRUE.
            CALL READF( RESTRICTREGIONRADIUS)
            CALL READF( RESTRICTREGIONX0)
            CALL READF( RESTRICTREGIONY0)
            CALL READF( RESTRICTREGIONZ0)

!
! js850> Make the region in RESTRICTREGION a cylinder rather than a sphere
!        The cylinder will have its axis aligned with the z axis
!
      ELSE IF (WORD.EQ.'RESTRICTCYL') THEN
            RESTRICTCYL=.TRUE.

! hk286 > Generalised rigid body   !
      ELSE IF (WORD.EQ.'RIGIDINIT') THEN

         RIGIDINIT = .TRUE.
         ATOMRIGIDCOORDT = .TRUE.
         AACONVERGENCET = .TRUE.
         !RIGID    = .TRUE.

      ELSE IF (WORD.EQ.'RGCL2') THEN
         RGCL2=.TRUE.

      ELSE IF (WORD.EQ.'RINGROTSCALE') THEN
         CALL READF(RINGROTSCALE)

      ELSE IF (WORD.EQ.'RKMIN') THEN
         RKMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
         WRITE(MYUNIT,'(A,2L5)') 'RKMIN branch RKMIN,BSMIN=',RKMIN,BSMIN

      ELSE IF (WORD.EQ.'RMS') THEN
         RMST=.TRUE.
         CALL READF(RMSLIMIT)
         CALL READF(RMSTOL)
         CALL READI(RMSSAVE)
         CALL READI(J1)
         CALL READI(J2)
         IF(J1.EQ.1) THEN
           SELECTT=.TRUE.
         ELSE
           SELECTT=.FALSE.
         ENDIF
         IF (J2.EQ.1) PROGRESS=.TRUE.
         WRITE(MYUNIT,'(A)') 'RMST set'
!
! csw34> Rigid body rotation moves. Each rigid body is randomly rotated about its COM every ROTATERIGIDFREQ steps.
!        ROTATEFACTOR scales the maximum rotation with 1.0 being complete freedom to rotate.
!
      ELSE IF (WORD.EQ.'ROTATERIGID') THEN
         ROTATERIGIDT=.TRUE.
! Read ROTATERIGIDFREQ 
         IF (NITEMS.GT.1) CALL READI(ROTATERIGIDFREQ)
! Read in ROTATEFACTOR
         IF (NITEMS.GT.2) CALL READF(ROTATEFACTOR)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body rotation moves enabled'
         WRITE(MYUNIT,'(A,I2,A)') ' keyword> Rigid bodies will be rotated every ',ROTATERIGIDFREQ,' steps'
         WRITE(MYUNIT,'(A,F20.10)') ' ROTATERIGID> Rigid body ROTATEFACTOR =',ROTATEFACTOR
   
!
! mo361> Rigid body translation moves
!        TRANSLATEFACTOR sets the maximum translation distance
!
      ELSE IF (WORD.EQ.'TRANSLATERIGID') THEN
         TRANSLATERIGIDT=.TRUE.
! Read TRANSLATERIGIDFREQ 
         IF (NITEMS.GT.1) CALL READI(TRANSLATERIGIDFREQ)
! Read in TRANSLATEFACTOR
         IF (NITEMS.GT.2) CALL READF(TRANSLATEFACTOR)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body translation moves enabled'
         WRITE(MYUNIT,'(A,I2,A)') ' keyword> Rigid bodies will be translated every ',TRANSLATERIGIDFREQ,' steps'
         WRITE(MYUNIT,'(A,F20.10)') ' TRANSLATERIGID> Rigid body TRANSLATEFACTOR =',TRANSLATEFACTOR
   
      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)
         IF (A9INTET.AND.(NSAVEINTE.EQ.0)) NSAVEINTE=NSAVE

      ELSE IF (WORD.EQ.'SAVEMULTIMINONLY') THEN
         SAVEMULTIMINONLY = .TRUE.

      ELSE IF (WORD.EQ.'SAVEINTE') THEN
         CALL READI(NSAVEINTE)

      ELSE IF (WORD.EQ.'SECPRED') THEN
         SECPREDT=.TRUE.
         CALL READA(UNSTRING)
         SECPREDFILE=TRIM(ADJUSTL(UNSTRING))

! dc550
      ELSE IF (WORD.EQ.'SELECTMOVE') THEN
         SELECTMOVET = .TRUE.
         CALL READI(SELMOVNO)
         CALL INISELECTMOVE()


      ELSE IF (WORD.EQ.'SC') THEN
         SCT=.TRUE.
         CALL READI(IX)
         NN=IX
         CALL READI(IX)
         MM=IX
         CALL READF(XX)
         SIG=XX
         CALL READF(XX)
         SCEPS=XX
         CALL READF(XX)
         SCC=XX

      ELSE IF (WORD.EQ.'SEED') THEN
         SEEDT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(IX)
            NSSTOP=IX
         ENDIF
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'SHELLMOVE') THEN
         MOVESHELLT=.TRUE.
         CALL READI(SHELLMOVEMAX)
         CALL READF(SHELLPROB)
         CALL READF(COREFRAC)

      ELSE IF (WORD.EQ.'SHIFTCUT') THEN
         SHIFTCUT=.TRUE.
         CALL READF(XX)
         CUTOFF=XX

!      ELSE IF (WORD.EQ.'SIDESTEP') THEN
!         CALL READF(SIDESTEP)
!         WRITE(MYUNIT,'(A,F14.10)') 'SIDESTEP=  ',SIDESTEP

!js850> soft sphere potential. note that SHIFT_CUT_NTYPEA and NTYPEA are
!equivalent and redundant.  NTYPEA is in a common block specific to potential
!BINARY, but this is also a binary potential, so some of the routines are the
!same.
      ELSE IF (WORD.EQ.'SOFT_SPHERE') THEN
         SOFT_SPHERE = .TRUE.
         CALL READI(SOFT_SPHERE_NTYPEA)
         NTYPEA = SOFT_SPHERE_NTYPEA

      ELSE IF (WORD.EQ.'SORT') THEN
         SORTT=.TRUE.

      ELSE IF (WORD.EQ.'SPARSE') THEN
         SPARSET=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ZERO_THRESH)

!     ELSE IF (WORD.EQ.'SQUEEZE') THEN
!        CALL READI(NVEC)
!        SQUEEZET=.TRUE.
!        IF (NITEMS.GT.2) THEN
!           CALL READF(XX)
!           SQUEEZER=XX
!        ENDIF
!        IF (NITEMS.GT.3) THEN
!           CALL READF(XX)
!           SQUEEZED=XX
!        ENDIF

      ELSE IF (WORD.EQ.'STAR') THEN
         STAR=.TRUE.
!
! Read in the maximum initial step size, factor for determining angular
! moves, and for rigid bodies the angular step size and the size of the
! blocks for Cartesian and angular moves.
!
! For parallel runs different values can be used for different runs by
! adding additional "STEP" lines to the data file. Otherwise the
! parameters for subsequent parallel runs are set to the values for the
! first one.
!
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))
         IF (NITEMS.GT.3) CALL READF(OSTEP(NPCOUNT))
         IF (NITEMS.GT.4) CALL READI(BLOCK(NPCOUNT))
!
!  Steered minimisation. This is for basin-hopping steps involving two well-defined
!  objects, e.g. a ligand and a protein.
!
      ELSE IF (WORD.EQ.'STEEREDMIN') THEN
         STEEREDMINT=.TRUE.
         CALL READF(SMINK)          ! final value of force constant
         CALL READF(SMINKINC)        ! increment of force constant per LBFGS step
         CALL READF(SMINDISTSTART)  ! initial distance for atoms SMINATOMA and SMINATOMB
         CALL READF(SMINDISTFINISH) ! final distance for atoms SMINATOMA and SMINATOMB (force turned off)
         CALL READI(SMINATOMA)      ! Atom A in the body to be rotated
         CALL READI(SMINATOMB)      ! Atom B in the other body (fixed for step)

      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.     
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'STEPOUT') THEN
         STEPOUT=.TRUE.

      ELSE IF (WORD.EQ.'STEPS') THEN
         NRUNS=1
         CALL READI(IX)
         MCSTEPS(1)=IX
         IF (NITEMS.GT.2) THEN
         CALL READF(XX)
         TFAC(1)=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            NRUNS=2
            CALL READI(IX)
            MCSTEPS(2)=IX
            CALL READF(XX)
            TFAC(2)=XX
         ENDIF
         IF (NITEMS.GT.5) THEN
            NRUNS=3
            CALL READI(IX)
            MCSTEPS(3)=IX
            CALL READF(XX)
            TFAC(3)=XX
         ENDIF

      ELSE IF (WORD.EQ.'STICKY') THEN
         STICKYT=.TRUE.
         RIGID=.TRUE.
         CALL READI(NRBSITES)
         CALL READF(STICKYSIG)
!        WRITE(MYUNIT,*) 'NRBSITES=',NRBSITES 
!        WRITE(MYUNIT,*) 'STICKYSIG=',STICKYSIG 
         ALLOCATE(SITE(NRBSITES,3))
         DO J1=1,NRBSITES
            READ(DATA_UNIT,*) SITE(J1,1:3)
!           CALL READF(SITE(J1,1))
!           CALL READF(SITE(J1,2))
!           CALL READF(SITE(J1,3))
!           WRITE(MYUNIT,'(A,I5,3G20.10)') 'J1,site: ',J1,SITE(J1,1:3)
         ENDDO

      ELSE IF (WORD.EQ.'LJCOUL') THEN
         LJCOULT=.TRUE.
         CALL READI(COULN)
         CALL READF(COULQ)
         CALL READF(COULSWAP)
         CALL READF(COULTEMP)
         NRBSITES=1
         ALLOCATE(SITE(NRBSITES,3))
!        Maybe the above two lines are not necessary!

      ELSE IF (WORD.EQ.'STOCK') THEN
         STOCKT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(STOCKMU)
         CALL READF(STOCKLAMBDA)
         ALLOCATE(SITE(NRBSITES,3))

!       Anisotropic potentials:

!     DC430 >

      ELSE IF (WORD.EQ.'CHECKD') THEN
         CHECKDT = .TRUE.
         IF (NITEMS .GT. 1) CALL READI(CHECKDID)

      ELSE IF (WORD.EQ.'CAPBIN') THEN

         CAPBINT = .TRUE.
         RIGID   = .TRUE.
         CALL READI(NPS)     ! Number of pentamers
         CALL READF(CAPEPS2)
         CALL READF(CAPRHO)    
         CALL READF(CAPRAD)
         CALL READF(CAPHEIGHT1)
         CALL READF(CAPHEIGHT2)        ! comment out this line for the previous model
         NRBSITES1 = 7        ! Number of sites defining a pentamer; change it to 6 for the previous model
         IF (NPS == NATOMSALLOC/2) THEN
            NRBSITES = 7      ! change it to 6 for the previous model    
         ELSE
            NRBSITES = 15     ! Sum of the numbers of sites defining a pentamer and a hexamer; change it to 13 for the previous model
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NRBSITES1*NPS + (NRBSITES-NRBSITES1)*(NATOMSALLOC/2-NPS)
         CALL DEFCAPBIN()

      ELSE IF (WORD .EQ. 'CHIRO') THEN

         CHIROT = .TRUE.
         RIGID = .TRUE.
         CALL READF(CHIRO_SIGMA)
         CALL READF(CHIRO_MU)
         CALL READF(CHIRO_GAMMA)
         IF (NITEMS > 3) THEN
             CALL READF(CHIRO_L)
         ELSE
             CHIRO_L = 0.D0
         END IF
         CALL INITIALIZE_CHIRO(MYUNIT)

      ELSE IF (WORD.EQ.'DBP') THEN

         DBPT   = .TRUE.
         RIGID  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBSIGBB)
         CALL READF(DBPMU)
         IF (NITEMS > 4) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF 
         NRBSITES = 3
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NATOMSALLOC*NRBSITES/2    !jdf43>

      ELSE IF (WORD.EQ.'DMBLPY') THEN

         DMBLPYT = .TRUE.
         RIGID   = .TRUE.
         CALL READF(YEPS) 
         CALL READF(YKAPPA)
         CALL READF(DBSIGBB)
         CALL READF(DBPMU)
         IF (NITEMS > 5) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 3
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NATOMSALLOC*NRBSITES/2  

! hk286
      ELSE IF (WORD.EQ.'DBYUKAWA') THEN
         DBYUKAWAT   = .TRUE.
         CALL READF(LAMBDAYAA)
         CALL READF(LAMBDAYAB)
         CALL READF(LAMBDAYBB)
         CALL READF(YEPSFAC)

      ELSE IF (WORD.EQ.'DBPTD') THEN

         DBPTDT = .TRUE.
         RIGID  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBSIGBB)
         CALL READF(DBPMU)
         IF (NITEMS > 4) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 3 
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = (NATOMSALLOC/2-1)*NRBSITES + 4    !jdf43>

      ELSE IF (WORD.EQ.'DMBLM') THEN

         DMBLMT   = .TRUE.
         RIGID    = .TRUE.
         CALL READF(EPS11)
         CALL READF(EPS22)
         CALL READF(MRHO11)
         CALL READF(MRHO22)
         CALL READF(REQ11)
         CALL READF(REQ22)
         CALL READF(DBPMU)
         IF (NITEMS > 8) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 2
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         ALLOCATE(DPMU(NRBSITES))

         CALL DEFDMBL()

      ELSE IF (WORD.EQ.'LINROD') THEN

         LINRODT = .TRUE.
         RIGID   = .TRUE.

         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFLINROD()

      ELSE IF (WORD.EQ.'LWOTP') THEN

         LWOTPT = .TRUE.
         RIGID  = .TRUE.
         IF (NITEMS > 1) CALL READF(LWRCUT)
         NRBSITES = 3
         ALLOCATE(SITE(NRBSITES,3)) 
         NTSITES = NATOMSALLOC*NRBSITES/2    !jdf43>
 
         CALL DEFLWOTP()

      ELSE IF (WORD.EQ.'NCAP') THEN

         NCAPT  = .TRUE.
         RIGID  = .TRUE.
         CALL READF(EPS2)
         CALL READF(RHO)
         CALL READF(RAD)
         CALL READF(HEIGHT)
         NRBSITES = 8
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NATOMSALLOC*NRBSITES/2    !jdf43>

         CALL DEFCAPSID(RAD,HEIGHT)

      ELSE IF (WORD .EQ. 'NPAH') THEN

         CALL READI(PAHID)         
         NPAHT    = .TRUE.
         RIGID    = .TRUE.
         IF (PAHID == 1) THEN
            NRBSITES = 36
         ELSEIF (PAHID == 2) THEN
            NRBSITES = 26
         ELSEIF (PAHID == 3) THEN
            NRBSITES = 12
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ. 'NTIP') THEN
 
         NTIPT =.TRUE.
         RIGID =.TRUE.
         IF (NITEMS > 1) THEN 
            CALL READI(TIPID)
         ELSE
            PRINT *, 'ERROR, TIPID is missing'
            STOP
         ENDIF 
         IF (TIPID == 1) NRBSITES = 3
         IF (TIPID == 2) NRBSITES = 4
         IF (TIPID == 3) NRBSITES = 3
         IF (TIPID == 4) NRBSITES = 4
         IF (TIPID == 5) NRBSITES = 5

         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NATOMSALLOC*NRBSITES/2    !jdf43>
         IF (TIPID == 4) THEN
             PI = 4.D0*DATAN(1.D0)
             ROH   = 0.9572D0
             ROM   = 0.15D0
             WTHETA = 104.52D0
             WTHETA = PI*WTHETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

             SITE(1,1) = 0.D0
             SITE(1,2) = 0.D0
             SITE(1,3) = 0.D0

             SITE(2,1) = 0.D0
             SITE(2,2) = SIN(0.5D0*WTHETA)*ROH
             SITE(2,3) = COS(0.5D0*WTHETA)*ROH

             SITE(3,1) = 0.D0
             SITE(3,2) = -SIN(0.5D0*WTHETA)*ROH
             SITE(3,3) = COS(0.5D0*WTHETA)*ROH

             SITE(4,1) = 0.D0
             SITE(4,2) = 0.D0
             SITE(4,3) = ROM

         ENDIF
         IF (PERMDIST) THEN ! correct all permutations allowed if perm.allow is not given explicitly
            IF (NPERMSIZE(1).EQ.NATOMSALLOC) NPERMSIZE(1)=NATOMSALLOC/2
         ENDIF


!|gd351>

      ELSE IF (WORD .EQ. 'PATCHY') THEN
 
         PATCHY =.TRUE.
         RIGID =.TRUE.
         IF (NITEMS > 1) THEN 
            CALL READI(NRBSITES)
         ELSE
            PRINT *, 'ERROR, NRBSITES is missing'
            STOP
         ENDIF 

         ALLOCATE(SITE(NRBSITES,3))

         SIGMASQ = (1.D0)**2
         RANGESQ = (1.9D0)**2
         FACTOR =  (2*3.14159265358979D0*0.05)**2

         CALL DEFINE_PATCHES(7.298D0)

      ELSE IF (WORD .EQ. 'ASAOOS') THEN
 
         ASAOOS =.TRUE.

          IF (NITEMS > 1) THEN 
            CALL READF(SIGMAP)
         ELSE
            SIGMAP=0.1D0
         ENDIF 

         CALL ASAOOSPRINT()

      ELSE IF (WORD .EQ. 'SANDBOX') THEN

         SANDBOXT = .TRUE.

      ELSE IF (WORD .EQ. 'SILANE') THEN

         SILANET = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 5
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFSILANE()

      ELSE IF (WORD .EQ. 'STOCKAA') THEN

         STOCKAAT = .TRUE.
         RIGID    = .TRUE.
         CALL READF(STOCKMU)
         IF (NITEMS > 2) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         IF (.NOT. EFIELDT) EFIELD = 0.D0
         NRBSITES = 1
         ALLOCATE(SITE(NRBSITES,3))   !jdf43>
         NTSITES = NATOMSALLOC*NRBSITES/2

      ELSE IF (WORD .EQ. 'MORSEDP') THEN

         MORSEDPT = .TRUE.
         RIGID    = .TRUE.
         CALL READF(RHO)
         CALL READF(STOCKMU)
         IF (NITEMS > 3) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         IF (.NOT. EFIELDT) EFIELD = 0.D0
         NRBSITES = 1
         ALLOCATE(SITE(NRBSITES,3))
         NTSITES = NATOMSALLOC*NRBSITES/2

      ELSE IF (WORD .EQ. 'MSSTK') THEN

         CALL READI(NRBSITES)
         MSSTOCKT = .TRUE.
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         ALLOCATE(DPMU(NRBSITES))
         NTSITES = NATOMSALLOC*NRBSITES/2    !jdf43>
         DO J1 = 1, NRBSITES
            CALL READF(DPMU(J1))
         ENDDO    
         IF (NRBSITES == 2) THEN
            CALL DEFMULT2STOCK()
         ELSE IF (NRBSITES == 3) THEN
            CALL DEFMULT3STOCK()
         ELSE IF (NRBSITES == 4) THEN
            CALL DEFMULT4STOCK()
         ELSE IF (NRBSITES == 5) THEN
            CALL DEFMULT5STOCK()
         ELSE
            WRITE(MYUNIT,*) 'NOT ALLOWED NRBSITES=',NRBSITES
            STOP 
         ENDIF
         IF (NITEMS > (2 + NRBSITES)) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         IF (.NOT. EFIELDT) EFIELD = 0.D0

      ELSE IF (WORD .EQ. 'MSBIN') THEN

!         CALL READI(NRBSITES)
!         CALL READI(NRBSITES1)
         CALL READI(NPS)
         NRBSITES  = 11
         NRBSITES1 = 5        
         MSTBINT = .TRUE.
         RIGID    = .TRUE.
         CALL READF(STOCKMU)
         IF (NITEMS > 3) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         CALL DEFMSTBIN()
         IF (.NOT. EFIELDT) EFIELD = 0.D0

      ELSE IF (WORD .EQ. 'MULTPAHA') THEN

         TPAHA     = 4
         ALLOCATE(NCMP(TPAHA))
         CALL READI (NCMP(1))
         CALL READI (NCMP(2))
         CALL READI (NCMP(3))
         CALL READI (NCMP(4))
         MULTPAHAT = .TRUE.
         RIGID     = .TRUE.

      ELSE IF (WORD .EQ. 'TDHD') THEN

         CALL READF(RHO)
         CALL READF(MREQ)
         CALL READF(EPSR)
         TDHDT    = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFTDHD()

      ELSE IF (WORD .EQ. 'DODECAMORSE') THEN

         CALL READF(RHO)
         CALL READF(MREQ)
         CALL READF(HSEFF)
         CALL READF(BEPS)
         DDMT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=20
         ALLOCATE(SITE(NRBSITES,3))
         CALL DEFDDM()
         IF (NITEMS.GT.5) CALL READF(DDMCUT)

      ELSE IF (WORD .EQ. 'WATERDC') THEN

         WATERDCT = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ. 'WATERKZ') THEN

         WATERKZT = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'GB') THEN

         GBT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         RIGID     = .TRUE.
         ESA       = 0.5D0*GBSIGNOT*(/1.D0, 1.D0, GBKAPPA/)
         GBCHI     = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM  = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

      ELSE IF (WORD.EQ.'GD') THEN

         GBDT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         RIGID     = .TRUE.
         SIGMAF     = GBSIGNOT * GBKAPPA
         ESA        = 0.5D0*(/GBSIGNOT, GBSIGNOT, SIGMAF/)
         INVKAP     = 1.D0/GBKAPPA
         GBCHI      = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM   = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

      ELSE IF (WORD.EQ.'GBDP') THEN

         GBDPT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)
         CALL READF(GBDPMU)
         CALL READF(GBDPEPS)

         RIGID     = .TRUE.
         SIGMAF     = GBSIGNOT * GBKAPPA
         ESA        = 0.5D0*(/GBSIGNOT, GBSIGNOT, SIGMAF/)
         INVKAP     = 1.D0/GBKAPPA
         GBCHI      = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         IF (GBMU == 0.D0) THEN
            GBCHIPRM = -1.D0
         ELSE
            GBCHIPRM   = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)
         ENDIF
         GBDPFCT    = 3.D0*GBDPEPS*GBDPMU*GBDPMU*GBSIGNOT**3.D0

      ELSE IF (WORD.EQ.'GEM') THEN

         GEMT   = .TRUE.
         CALL READF(GEMRC)

!     ----------------------------------------------------------------------------------------------

      ELSE IF (WORD .EQ. 'PAHA') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 12
         ELSEIF (PAHID == 2) THEN
            NRBSITES = 18
         ELSEIF (PAHID == 3) THEN
            NRBSITES = 24
         ELSEIF (PAHID == 4) THEN
            NRBSITES = 26
         ENDIF

         PAHAT    = .TRUE.
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBSTLA(NRBSITES,3))
         ALLOCATE(STCHRG(NRBSITES))
         
         NTSITES = (NATOMSALLOC/2)*NRBSITES  !jdf43>

         CALL DEFPAHA()

         IF (PAHID == 1) THEN
            NCARBON  = 6
            CALL DEFBENZENE()
         ELSEIF (PAHID == 2) THEN
            NCARBON  = 10
            CALL DEFNAPHTHALENE()
         ELSEIF (PAHID == 3) THEN
            NCARBON  = 14
            CALL DEFANTHRACENE()
         ELSEIF (PAHID == 4) THEN
            NCARBON  = 16
            CALL DEFPYRENE()
         ENDIF

!     ----------------------------------------------------------------------------------------------

      ELSE IF (WORD .EQ. 'PAHW99') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 18
         ENDIF

         PAHW99T  = .TRUE.
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ. 'PAP') THEN

         CALL READI(NPATCH)
         CALL READF(PAPALP)
         CALL READF(PAPS)
         CALL READF(PAPCD)
         CALL READF(PAPEPS)
         
         NRBSITES = 2*NPATCH

         PAPT   = .TRUE.
         RIGID  = .TRUE.

      ELSE IF (WORD .EQ. 'PAPBIN') THEN

         CALL READF(PAPEPS)
         CALL READF(PAPS)
         CALL READF(PAPANG1)
         CALL READF(PAPANG2)
         CALL READF(YKAPPA)

         NRBSITES = 2
         ALLOCATE(RBSTLA(NRBSITES,3))

         CALL DEFPAPBIN()

         PAPBINT = .TRUE.
         RIGID   = .TRUE.

      ELSE IF (WORD .EQ. 'PAPJAN') THEN

         CALL READF(PAPEPS)
         CALL READF(PAPS)
         CALL READF(PAPANG1)
         CALL READF(YKAPPA)

         NRBSITES = 2
         ALLOCATE(RBSTLA(NRBSITES,3))

         CALL DEFPAPJANUS()

         PAPJANT = .TRUE.
         RIGID   = .TRUE.

      ELSE IF (WORD .EQ. 'PTSTST') THEN

         CALL READI(NPATCH)
         CALL READF(PAPEPS)
         CALL READF(PAPCD)
         CALL READF(YKAPPA)

         NRBSITES = NPATCH
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFPTSTST()

         PTSTSTT = .TRUE.
         RIGID   = .TRUE.

! js850> stop the monte carlo loop when MAX_NPCALL is reached.
! This can be useful for benchmarking where you want to see 
! the lowest energy reached after a certain number of energy evaluations
       ELSE IF (WORD .EQ.'MAX_NPCALL') THEN
          MAX_NPCALLT = .TRUE.
          CALL READI(MAX_NPCALL)

       ELSE IF (WORD .EQ.'MSGB') THEN

         MSGBT = .TRUE.
!         RIGID = .TRUE.
         CALL READI(NGBSITE)
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(SIGNOT)
         CALL READF(EPSNOT)
!         ALLOCATE(SITE(NRBSITES,3))

         LPL    = 0.5D0 * SIGNOT * GBKAPPA
         LPR    = 0.5D0 * SIGNOT
         LPRSQ  = LPR * LPR
         LSQDFR = LPL * LPL - LPRSQ

      ELSE IF (WORD .EQ. 'PY') THEN
         ! Syntax: PY sig_0 eps_0 [cut] [XYZ boxx boxy boxz]

         PYT = .TRUE.
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         RIGID = .TRUE.
         ! Rigid body SITE, NRBSITES, NTSITES information specified in py_input routine

         ! Specify cutoff for potential in absolute units
         IF (NITEMS.GT.3) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            WRITE(MYUNIT,*) "multisitepy cutoff: ", PCUTOFF
         ENDIF
         ! Specify periodic boundary conditions (PBCs)
         IF (NITEMS.GT.4) THEN
            ! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.            
            ! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
            CALL READA(PBC)
            BOXLX=0
            BOXLY=0
            BOXLZ=0
            IF (SCAN(PBC,'Xx').NE.0) THEN
                PARAMONOVPBCX=.TRUE.
                CALL READF(BOXLX)
                BOXLX = BOXLX*PCUTOFF
                WRITE(MYUNIT,*) "PBC X:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY = BOXLY*PCUTOFF
                WRITE(MYUNIT,*) "PBC Y:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ = BOXLZ*PCUTOFF
                WRITE(MYUNIT,*) "PBC Z:",BOXLZ
            ENDIF
         ENDIF

      ELSE IF (WORD.EQ.'PYOVERLAPTHRESH') THEN
         CALL READF(PYOVERLAPTHRESH)
         WRITE(MYUNIT,'(A,F8.3)') 'keyword> ellipsoids considered to overlap for an ECF value below ', PYOVERLAPTHRESH
         IF(NITEMS.GT.2) CALL READF(PYCFTHRESH)
         WRITE(MYUNIT,'(A,F8.3)') 'keyword> cold fusion will be diagnosed for an ECF value below ', PYCFTHRESH
 
      ELSE IF (WORD.EQ.'GAYBERNE') THEN

         GAYBERNET=.TRUE.
         ELLIPSOIDT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(GBANISOTROPYR)
         CALL READF(GBWELLDEPTHR)
         CALL READF(PSIGMA0(1))
         CALL READF(PEPSILON0)
         CALL READF(GBMU)
         CALL READF(GBNU)
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ.'PYGPERIODIC') THEN

         PYGPERIODICT = .TRUE.
         ELLIPSOIDT = .TRUE.
         RIGID = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         PARAMa1=PYA1(1)
         PARAMb1=PYA1(2)
         PARAMc1=PYA1(3)

         IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMSALLOC/2,3))
         IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMSALLOC/2,3))
         DO J1=1,NATOMSALLOC/2
           PYA1bin(J1,:)=PYA1(:)
           PYA2bin(J1,:)=PYA2(:)
         END DO
         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF

         IF (NITEMS.GT.9) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*PYSIGNOT
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         ENDIF
         IF (NITEMS.GT.10) THEN
! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.
! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
            CALL READA(PBC)
            write (*,*) "PBCs are: ",PBC
            BOXLX=0
            BOXLY=0
            BOXLZ=0
            IF (SCAN(PBC,'Xx').NE.0) THEN
                PARAMONOVPBCX=.TRUE.
                CALL READF(BOXLX)       ! BOXLX is a scaling factor, not the actual box length!
                BOXLX=BOXLX*PCUTOFF     ! now BOXLX is the actual box length
                write(*,*) "Paramonov Periodic Boundary Condition X active. BOXLX:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY=BOXLY*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Y active. BOXLY:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ=BOXLZ*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Z active. BOXLZ",BOXLZ
            ENDIF
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD .EQ.'LJCAPSID') THEN
!         Three-site Lennard-Jones based capsid model. Sites at the origin are standard LJ sites, the two apex sites 
!         are repulsive LJ sites, polarised. The site in the middle only interacts with sites in the middle of other 
!         molecules.
         LJCAPSIDT = .TRUE.
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         CALL READF(PEPSILON1(1))
         CALL READF(PSCALEFAC1(1))
         CALL READF(PSCALEFAC2(1))
         MAXINTERACTIONS=4
        
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'EXTRALJSITE') THEN
         LJSITE=.TRUE.
         CALL READF(PEPSILON1(1))
         CALL READF(PSCALEFAC1(1))
          MAXINTERACTIONS=1
         IF(NITEMS.GT.3) THEN
          CALL READF(PSCALEFAC2(1))
          WRITE(MYUNIT,'(A,3F8.3)') 'keyword> primary and secondary apex sites will be used, epsilon and heights: ', 
     &                              PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)
          IF(.NOT.LJSITEATTR) THEN
                MAXINTERACTIONS=3
          ELSE
                MAXINTERACTIONS=4
          END IF
         ELSE
          WRITE(MYUNIT,'(A,2F8.3)') 'keyword> primary apex sites will be used, epsilon and height: ', PEPSILON1(1), PSCALEFAC1(1)
         END IF
         IF(NITEMS.GT.4) THEN           ! binary ellipsoidal clusters will be set up only for two apex sites, not one
           BLJSITE=.TRUE.               ! we also won't use the sigma parameter from now on, epsilon is enough for repulsive sites
           WRITE(MYUNIT,'(A,3F8.3)') 'keyword> binary system with primary and secondary apex sites, ' //  
     &  'epsilon and heights for 2nd type particle: ', PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)

           CALL READF(PEPSILON1(2))
           CALL READF(PSCALEFAC1(2))
           CALL READF(PSCALEFAC2(2))
           CALL READF(PEPSILON1(3))     ! this is epsilon for the interaction between A and B type ellipsoids
           MAXINTERACTIONS=3 ! attractive secondary apex sites not incorporated for binary systems
         END IF
      ELSE IF (WORD.EQ.'EXTRALJSITEATTR') THEN
         LJSITE=.TRUE.
         LJSITEATTR=.TRUE.
         CALL READF(PSIGMAATTR(1))
         CALL READF(PEPSILONATTR(1))
         CALL READF(PSIGMAATTR(2))
         CALL READF(PEPSILONATTR(2))

         WRITE(MYUNIT,'(A,4F13.4)') 'keyword> primary and secondary apex sites '//
     &                             'with normal LJ attraction, sigmas and epsilons: ', 
     &                             PSIGMAATTR(1), PEPSILONATTR(1), PSIGMAATTR(2), PEPSILONATTR(2)
         MAXINTERACTIONS=4
      ELSE IF (WORD.EQ.'LJSITECOORDS') THEN
           LJSITECOORDST=.TRUE.
           CALL READF(LJSITECOORDS(1))
           CALL READF(LJSITECOORDS(2))
           CALL READF(LJSITECOORDS(3))
      ELSE IF (WORD.EQ.'SWAPMOVES') THEN
         SWAPMOVEST=.TRUE.
         IF(PYBINARYT) THEN
                PYSWAP(1) = 1
                PYSWAP(2) = PYBINARYTYPE1 + 1
                PYSWAP(3) = 1
                IF(NITEMS.GT.1) CALL READI(PYSWAP(3))
                WRITE(MYUNIT,'(A,I5,A)') 'keyword> ',PYSWAP(3), ' pairs of atoms will be swapped at once'
         END IF

      ELSE IF (WORD.EQ.'PYBINARY') THEN
         PYBINARYT=.TRUE.
         ELLIPSOIDT=.TRUE.
         RADIFT=.TRUE.
         CALL READI(PYBINARYTYPE1)
         CALL READF(PYA11(1))
         CALL READF(PYA11(2))
         CALL READF(PYA11(3))
         CALL READF(PYA21(1))
         CALL READF(PYA21(2))
         CALL READF(PYA21(3))
         CALL READF(PYA12(1))
         CALL READF(PYA12(2))
         CALL READF(PYA12(3))
         CALL READF(PYA22(1))
         CALL READF(PYA22(2))
         CALL READF(PYA22(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         IF(NITEMS.GT.16) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*PYSIGNOT
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         END IF
         IF(SWAPMOVEST) THEN
                PYSWAP(1) = 1
                PYSWAP(2) = PYBINARYTYPE1 + 1
         END IF

         IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMSALLOC/2,3))
         IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMSALLOC/2,3))
         DO J1=1,NATOMSALLOC/2
          IF(J1<=PYBINARYTYPE1) THEN
           PYA1bin(J1,:)=PYA11(:)
           PYA2bin(J1,:)=PYA21(:)
          ELSE
           PYA1bin(J1,:)=PYA12(:)
           PYA2bin(J1,:)=PYA22(:)
          END IF
         END DO

      ELSE IF (WORD.EQ.'GAYBERNEDC') THEN
         GAYBERNEDCT=.TRUE.
!         ELLIPSOIDT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(SIGNOT)
         CALL READF(EPSNOT)
         ALLOCATE(SITE(NRBSITES,3))

!
!  Eigenvalue shift parameter.
!
      ELSE IF (WORD .EQ. 'SHIFT') THEN
         CALL READF(SHIFTV)

      ELSE IF (WORD.EQ.'STRAND') THEN
         STRANDT=.TRUE.
         RIGID=.TRUE.
!
!  The nine reference site positions per strand.
!
         NRBSITES=9
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=-2.7298862082
         SITE(1,2)=2.3622865625 
         SITE(1,3)=0.6475151629
         SITE(2,1)=-1.7492122114
         SITE(2,2)=2.3331194664 
         SITE(2,3)=0.5887015133
         SITE(3,1)=-1.5963638586
         SITE(3,2)=1.4304320585 
         SITE(3,3)=0.2442792479
         SITE(4,1)=-0.6166461313
         SITE(4,2)=1.4301805389 
         SITE(4,3)=0.1327546571
         SITE(5,1)=-0.4460267836
         SITE(5,2)=0.5254809645  
         SITE(5,3)=-0.2196837962
         SITE(6,1)=0.5313983749 
         SITE(6,2)=0.5210707739  
         SITE(6,3)=-0.3409645197
         SITE(7,1)=0.7065341613  
         SITE(7,2)=-0.3914277962 
         SITE(7,3)=-0.6719579835 
         SITE(8,1)=1.6776397940  
         SITE(8,2)=-0.3830053500 
         SITE(8,3)=-0.8355266604
         SITE(9,1)=1.8162689403  
         SITE(9,2)=-1.3093381947 
         SITE(9,3)=-1.1427874015
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'SUPERSTEP') THEN
         SUPERSTEP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSUPER)
         IF (NITEMS.GT.2) CALL READF(SUPSTEP)
         IF (NITEMS.GT.3) CALL READF(TEMPS)
         IF (NITEMS.GT.4) CALL READF(SACCRAT)
         IF (NITEMS.GT.5) CALL READI(NSACCEPT)

      ELSE IF (WORD.EQ.'SW') THEN
         SW=.TRUE.
      ELSE IF (WORD.EQ.'QUIP') THEN
         QUIPT=.TRUE.
         QUIPATOMTYPE='Ag '
         QUIPEQDIST=2.0D0**(1.0D0/6.0D0)
         QUIPARGSTR='IP LJ'
         IF (NITEMS.GT.1) THEN
            CALL READA(QUIPATOMTYPE)
            CALL READF(QUIPEQDIST)
            CALL READA(QUIPARGSTR)
         ENDIF
      ELSE IF (WORD.EQ.'SETCHIRAL') THEN
         SETCHIRAL=.TRUE.
      ELSE IF (WORD.EQ.'SETCHIRALGENERIC') THEN
         SETCHIRALGENERIC=.TRUE.
!
      ELSE IF (WORD .EQ. 'PRINT_PTGRP') THEN
         PRINT_PTGRP = .TRUE.
         IF (NITEMS.GT.1) CALL READF(SYMTOL1)
         IF (NITEMS.GT.2) CALL READF(SYMTOL2)
         IF (NITEMS.GT.3) CALL READF(SYMTOL3)

!  Keyword and parameters for symmetrisation.
!
      ELSE IF (WORD.EQ.'SYMMETRISE') THEN
         SYMMETRIZE=.TRUE.
         NCORE=0
         IF (NITEMS.GT.1) CALL READI(NSYMINTERVAL)
         IF (NITEMS.GT.2) CALL READF(SYMTOL1)
         IF (NITEMS.GT.3) CALL READF(SYMTOL2)
         IF (NITEMS.GT.4) CALL READF(SYMTOL3)
         IF (NITEMS.GT.5) CALL READF(SYMTOL4)
         IF (NITEMS.GT.6) CALL READF(SYMTOL5)
         IF (NITEMS.GT.7) CALL READI(NSYMQMAX)
         IF (NITEMS.GT.8) CALL READF(MATDIFF) ! appears to have little effect now
         IF (NITEMS.GT.9) CALL READF(DISTFAC)
!
!  Keyword and parameters for symmetrisation according to a continuous symmetry measure.
!
      ELSE IF (WORD.EQ.'SYMMETRISECSM') THEN
         SYMMETRIZE=.TRUE.
         SYMMETRIZECSM=.TRUE.
         CSMT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSYMINTERVAL)
         IF (NITEMS.GT.2) THEN
            CALL READA(CSMGP)
            CALL MYUPCASE(CSMGP)
         ELSE
            PRINT '(A)','keyword> ERROR - point group must be specified for SYMMETRIZECMS keyword'
            STOP
         ENDIF
         IF (NITEMS.GT.3) CALL READF(CSMEPS)
         IF (NITEMS.GT.4) CALL READI(CSMSTEPS)
         IF (NITEMS.GT.5) CALL READI(CSMQUENCHES)
         IF (NITEMS.GT.6) CALL READI(CSMMAXIT)
         IF (.NOT.PERMDIST) THEN
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
!        ALLOCATE(NPERMSIZE(NATOMSALLOC),PERMGROUP(NATOMSALLOC),NSWAP(NATOMSALLOC),SWAP1(NATOMSALLOC,2),SWAP2(NATOMSALLOC,2))
         ALLOCATE(NPERMSIZE(3*NATOMSALLOC),PERMGROUP(3*NATOMSALLOC),NSETS(3*NATOMSALLOC),SETS(NATOMSALLOC,70))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
!
!  Sanity checks!
!
               IF (NSETS(J1).GT.13) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 70'
                  STOP
               ENDIF
!              IF (NDUMMY+NPERMSIZE(J1).GT.NATOMSALLOC) THEN
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMSALLOC) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))

               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
            CLOSE(1)
!
!  And yet another!
!  
            IF (NFREEZE.GT.0) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
         ELSE
            NSETS(1:NATOMSALLOC)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMSALLOC ! all atoms can be permuted - default
            DO J1=1,NATOMSALLOC
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            WRITE(MYUNIT,'(A,3(I6,A))') ' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(MYUNIT,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1)
            IF (NSETS(J1).GT.0) THEN
               WRITE(MYUNIT,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(MYUNIT,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' / '
                  ENDDO
                  IF (J2.LT.NSETS(J1)) WRITE(MYUNIT,'(A3)',ADVANCE='NO') ' ; '
               ENDDO
            ENDIF
            WRITE(MYUNIT,'(A)') ' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO

         ENDIF
      ELSE IF (WORD.EQ.'TABOO') THEN
         TABOOT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NTAB)

      ELSE IF (WORD.EQ.'TARGET') THEN
         TARGET=.TRUE.

  ! beh35 > new handling of TARGET logical to account for mlowest

        IF(MLOWEST) THEN
           PRINT *,"Cannot have both TARGET and MLOWEST"
           STOP
        END IF

         NTARGETS=NITEMS-1
         ALLOCATE(TARGETS(NTARGETS))
         INQUIRE(FILE='coords.target',EXIST=YESNO)
         IF (YESNO) THEN
            ALLOCATE(TCOORDS(NTARGETS,3*NATOMSALLOC))
            OPEN(UNIT=1,FILE='coords.target',STATUS='OLD')
            READ(1,*) ((TCOORDS(J1,J2),J2=1,3*NATOMSALLOC),J1=1,NTARGETS)
            CLOSE(1)
         ENDIF
         DO J1=2,NITEMS
            CALL READF(XX)
            TARGETS(J1-1)=XX
         ENDDO
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'TD') THEN
         FIELDT=.TRUE.
         TDT=.TRUE.
         CALL READF(XX)
         FTD=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF
!
! Tethered WL walk to determine anharmonic vibrational density of states
!
      ELSE IF (WORD.EQ.'TETHER') THEN
         TETHER=.TRUE.
         CALL READF(hdistconstraint) 
        CALL READI(hwindows)
         lhbins=int(hbins/hwindows)
         CALL READF(ExtrapolationPercent)
         lhbins=int(hbins/hwindows)
         sampledbins=int((1.0d0-ExtrapolationPercent)*hbins/hwindows)
         CALL READF(lnHarmFreq)
      ELSE IF (WORD.EQ.'THOMSON') THEN
         THOMSONT=.TRUE.
         ODDCHARGE=1.0D0
         IF (NITEMS.GT.1) CALL READF(ODDCHARGE)
!
!  Threshold acceptance rather than Metropolis, i.e. the energy change
!  can;t increase by more than a certain amount.
!  NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'THRESHOLD') THEN
         THRESHOLDT=.TRUE.
!        WRITE(MYUNIT,*) 'keyword THRESHOLD doesnt appear to do anything at the moment'
         STOP

      ELSE IF (WORD.EQ.'TIP') THEN
         TIP=.TRUE.
         RIGID=.TRUE.
         IF (NITEMS.GT.1) CALL READI(TIPID)
         IF (TIPID.EQ.5) NRBSITES=5
         IF (TIPID.EQ.4) NRBSITES=4
         IF (TIPID.EQ.3) NRBSITES=3
         IF (TIPID.EQ.2) NRBSITES=4
         IF (TIPID.EQ.1) NRBSITES=3
         ALLOCATE(SITE(NRBSITES,3))
!     ELSE IF (WORD.EQ.'TN') THEN
!        TNT=.TRUE.
!        WRITE(MYUNIT,'(A)') 'optimisation with tn no longer supported'
!        STOP

      ELSE IF (WORD.EQ.'TOLBRENT') THEN
         CALL READF(TOLB)
!
! NOT DOCUMENTED
!
      ELSE IF (WORD .EQ. 'TOLD') THEN
        CALL READF(XX)
        TOLD=XX
!
! NOT DOCUMENTED
!
      ELSE IF (WORD .EQ. 'TOLE') THEN
        CALL READF(XX)
        TOLE=XX

      ELSE IF (WORD.EQ.'TOSI') THEN
         TOSI=.TRUE.
         CALL READF(APP)
         CALL READF(AMM)
         CALL READF(APM)
         CALL READF(RHO)
!
!  Set Tsallis statistics with some q value.
!
      ELSE IF (WORD.EQ.'TSALLIS') THEN
         TSALLIST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(QTSALLIS)
         ENDIF
      ELSE IF (WORD.EQ.'NEWTSALLIS') THEN
         NEWTSALLIST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(QTSALLIS)
         ENDIF
         
!
!  Xantheas' TTM3-F water potential
!
      ELSE IF (WORD.EQ.'TTM3') THEN
        TTM3T=.TRUE.

!! ns644> Adding dihedral angle potential via TWIST keyword: 
      ELSE IF (WORD.EQ.'TWIST') THEN
         TWISTT=.TRUE. 
         INQUIRE(FILE='twistgroups',EXIST=YESNO)
         IF (NITEMS.EQ.7) THEN
            NTWISTGROUPS = 1
            ALLOCATE(TWIST_K(1))
            ALLOCATE(TWIST_THETA0(1))
            ALLOCATE(TWIST_ATOMS(1,1))
            CALL READI(TWIST_ATOMS(1,1))
            CALL READI(TWIST_ATOMS(2,1))
            CALL READI(TWIST_ATOMS(3,1))
            CALL READI(TWIST_ATOMS(4,1))
            CALL READF(TWIST_K(1))
            CALL READF(TWIST_THETA0(1))
         ELSE IF (NITEMS.EQ.1 .AND. YESNO) THEN
            OPEN(UNIT=444,FILE='twistgroups',status='unknown')
            WRITE(MYUNIT,*) 'keyword> Reading in twistgroups..'
            READ(444,*) NTWISTGROUPS
            ALLOCATE(TWIST_K(NTWISTGROUPS))
            ALLOCATE(TWIST_THETA0(NTWISTGROUPS))
            ALLOCATE(TWIST_ATOMS(4,NTWISTGROUPS))
            DO J1=1, NTWISTGROUPS
               READ(444,*) TWIST_ATOMS(1,J1),TWIST_ATOMS(2,J1), 
     &                     TWIST_ATOMS(3,J1),TWIST_ATOMS(4,J1),
     &                     TWIST_K(J1),TWIST_THETA0(J1)          
            END DO 
         ELSE
            WRITE(MYUNIT,*) 'Provide twistgroups file, or specify
     &                          details of single dihedral in data'
            STOP
         END IF

      ELSE IF (WORD.EQ.'TWOPLUS') THEN
         TWOPLUS=.TRUE.
!
!     Specify 2D XY APBC potential
!
      ELSE IF (WORD.EQ.'TWODAPBC') THEN                                        
         TWODAPBCT=.TRUE. 
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)  
         IF (MOD(NONEDAPBC,3).NE.0) THEN        
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'    
            STOP
         ENDIF

!
!     Specify 2D XY PBC potential 
!      
      ELSE IF (WORD.EQ.'TWODPBC') THEN   
         TWODPBCT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
         IF (MOD(NONEDAPBC,3).NE.0) THEN 
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'  
            STOP 
         ENDIF 

!                                                                                                                                                                                                          
!     Specify 3D XY APBC potential                                                                                                                                                                         
!                                                                                                                                                                                                          
      ELSE IF (WORD.EQ.'THREEDAPBC') THEN                                                                                                                                       
         THREEDAPBCT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
         IF (MOD(NONEDAPBC,3).NE.0) THEN
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'                                                                                                         
            STOP                                                                                                                                                                                            
         ENDIF 


!                                                                                                                                                                                                          
!     Specify 3D XY PBC potential                                                                                                                                                                         
!                                                                                                                                                                                                          
      ELSE IF (WORD.EQ.'THREEDPBC') THEN                                                                                                                                                                 
         THREEDPBCT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
         IF (MOD(NONEDAPBC,3).NE.0) THEN
            WRITE(MYUNIT,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'                                                                                                         
            STOP                                                                                                                                                                                            
         ENDIF 
!
!  Set proposed atomic steps in takestep to be uniformly distributed within
!  a sphere of given radius.
!
      ELSE IF (WORD.EQ.'UNIFORMMOVE') THEN
         UNIFORMMOVE=.TRUE.

!
! csw34> Update the reference coordinates for the generalised rigid bodies after a step has been taken. 
!        This allows steps to be taken WITHIN the rigid bodies, although HYBRIDMIN should also be used
!        as any bad conformation intrduced by these step will otherwise be frozen into the rigid bodies.
!
      ELSE IF (WORD.EQ.'UPDATERIGIDREF') THEN
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body reference coordinates will be updated after each step'
         WRITE(MYUNIT,'(A)') ' UPDATERIGIDREF> WARNING: make sure HYBRIDMIN is enabled!'
         UPDATERIGIDREFT=.TRUE.

!
!  Number of BFGS updates before resetting, default=4
!
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
!
!  Whether to include rotational partition function for FEBH
!
      ELSE IF (WORD.EQ.'USEROT') THEN
         USEROT=.TRUE.
!      
! vr274
! This is supposed to be a potential which is linked in from a different file
! TODO: to some checking, e.g. give a string which defines the userpotential
! and checks whether the correct binary is used, etc
!
      ELSE IF (WORD.EQ.'USERPOT') THEN
         USERPOTT=.true.
         CALL USERPOT_INITIALIZE_GMIN(3*NATOMSALLOC, COORDS(:,1))

!
!  Use VGW (Variational Gaussian Wavepacket) Minimization (Quantum Quenching)
!
      ELSE IF (WORD.EQ.'VGW') THEN
         VGW=.TRUE.
         LBFGST=.FALSE.
         CALL READF(LJSIGMA)
         CALL READF(LJEPSILON)
         CALL READF(TAUMAX)
         CALL READF(TAUMAXFULL)
                 
      ELSE IF (WORD.EQ.'VGWCPS') THEN
         CALL READI(CPS)
         CALL READF(CPFACTORSG)

      ELSE IF (WORD.EQ.'VGWCPF') THEN
         CALL READI(CPF)
         CALL READF(CPFACTORFG)

      ELSE IF (WORD.EQ.'VGWTOL') THEN
         CALL READF(VGWTOL)

!
! Choice of convergence regime for Wang-Landau runs: histogram flatness (default), 
! VisitProp - minimal number of visits proportional to 1/sqrt(ln(f))
!
      ELSE IF (WORD.EQ.'VISITPROP') THEN
         VISITPROP=.TRUE.
!
! Maximum PE for an instantaneous configuration above a basin bottom in BSPT
!
      ELSE IF (WORD.EQ.'TSTAR') THEN
         CALL READF(TSTAR)
      ELSE IF (WORD.EQ.'WELCH') THEN
         WELCH=.TRUE.
         CALL READF(APP)
         CALL READF(AMM)
         CALL READF(APM)
         CALL READF(RHO)
         CALL READF(XQP)
         CALL READF(XQM)
         CALL READF(ALPHAP)
         CALL READF(ALPHAM)
!
! NOT DOCUMENTED
!
      ELSE IF (WORD.EQ.'WENZEL') THEN
         WENZEL=.TRUE.

      ELSE IF (WORD.EQ.'ZETT1') THEN
         ZETT1=.TRUE.

      ELSE IF (WORD.EQ.'ZETT2') THEN
         ZETT2=.TRUE.

! jdf43>
      ELSE IF (WORD.EQ.'TBP') THEN
         TBP=.TRUE.
         CALL READF(TBPMIN)
         CALL READF(TBPSTEP)
         CALL READI(TBPSTEPS)
         CALL READF(TBPHF)
         CALL READF(TBPCF)
         CALL READI(TBPCI)
         ALLOCATE(TBPBASINS(TBPSTEPS),TBPBASINSO(TBPSTEPS)) 

      ELSE IF (WORD.EQ.'RATIO') THEN
         RATIOT=.TRUE.
         CALL READF(SRATIO)
         CALL READF(TRATIO)
         FIXBOTH=.TRUE.

      ELSE IF (WORD.EQ.'MWFILM') THEN
         MWFILMT=.TRUE.

      ELSE IF (WORD.EQ.'SWPOT') THEN
         SWPOTT=.TRUE.
         CALL SWINIT

      ELSE IF (WORD.EQ.'MWPOT') THEN
         MWPOTT=.TRUE.
         CALL MWINIT

      ELSE IF (WORD.EQ.'SUPPRESS') THEN
         SUPPRESST=.TRUE.

      ELSE IF (WORD.EQ.'MFET') THEN
         MFETT=.TRUE.
         CALL READI(NPAR)
         CALL READF(MFETPCTL)
         CALL READF(MFETTRGT)
         ALLOCATE(TARGETCOORDS(NATOMSALLOC*3))
!         OPEN(UNIT=626,FILE='target.coords')
!         DO J1=1,NATOMSALLOC/NPAR
!            J2=3*(J1-1)
!            READ(626,*) TARGETCOORDS(J2+1), TARGETCOORDS(J2+2), TARGETCOORDS(J2+3)
!         ENDDO
!         CLOSE(626)
         NATOMSALLOC=NATOMSALLOC/NPAR
         NRUNS=1

      ELSE IF (WORD.EQ.'POLIR') THEN
         POLIRT=.TRUE.

      ELSE IF (WORD.EQ.'MBPOL') THEN
         CALL MBPOLINIT
         MBPOLT=.TRUE.
         MOLECULART=.TRUE.
         ALLOCATE(TYPECH(NATOMSALLOC))
         TYPECH(1::3)(1:1)="O"
         TYPECH(2::3)(1:1)="H"
         TYPECH(3::3)(1:1)="H"
      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF
      CALL FLUSH(MYUNIT)
      
!ds656> NTYPEA can fluctuate in homotop refinemet routine, so we need
!     a fixed reference. It is set here instead of every single
!     IF block where NTYPEA can potentially change. 
      NTYPEA_FIX = NTYPEA
      
      GOTO 190
            
      RETURN
      END
