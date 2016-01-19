      MODULE MODAMBER9 
!      use commons, only : natoms        ! COORDS deleted

      IMPLICIT NONE
      SAVE

! Stuff for dumping normal mode info
      logical, dimension(:), allocatable :: DUMPMODEN
      logical :: KTWNT
      double precision :: KTWN

      LOGICAL MDSTEPT, readcoords, NOCISTRANSRNA, CHECKCISTRANSALWAYS, GOODSTRUCTURE1, GOODSTRUCTURE2, AMBERENERGIEST
      LOGICAL NOCISTRANSDNA, UACHIRAL, SETCHIRAL, SETCHIRALGENERIC, CHECKCISTRANSALWAYSRNA, CHECKCISTRANSALWAYSDNA
      CHARACTER(len=20) :: amberstr
      CHARACTER(len=8)  :: amberstr1
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: coords1, coords, atmass1, dihedralsave
      INTEGER, DIMENSION(:), ALLOCATABLE :: cisarray1, cisarray2, chiarray1, chiarray2
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: atomindex
      LOGICAL, DIMENSION(:), ALLOCATABLE :: exclude
      CHARACTER(len=81) :: prmtop
      DOUBLE PRECISION, PARAMETER :: TEMPPARAM = 0.0000001
! LIGAND MOVE PARAMETERS
! ligrotscale = scaling factor for rigid rotations (0->1)
! ligcartstep = size of random cartesian coordinate move for ligand
! ligtransstep = size of rigid translation move for ligand
! ligmovefreq = frequency of ligand moves i.e. 1 = every step
! doligmove = logical switch to tell amberinterface.f to do the ligand move
      DOUBLE PRECISION :: ligrotscale, ligcartstep, ligtransstep
      INTEGER :: ligmovefreq
      LOGICAL :: doligmove=.FALSE.
! END OF LIGAND MOVE PARAMETERS

! ROTAMER MOVE PARAMETERS
! rotamert = logical switch indicating rotamer moves should be used
! rotmaxchange = maximum number of rotamers to change in one step
! rotpselect = the probability of choosing a particular rotamer
! rotoccuw = the amount the occupation frequency affects the rotamer selection i.e.
!            rotoccuw = 0.0 -> all rotamers possible
!            rotoccuw = 0.004 -> only rotamers with a > 0.4% occupation selected
! rotcentre = residue id used to specify centre of rotamer moves (ligand maybe?)
! rotcuttoff = selection probability decays linearly from the centre (above) to the cutoff
      LOGICAL :: ROTAMERT
      INTEGER :: ROTMAXCHANGE, ROTCENTRE
      DOUBLE PRECISION :: ROTPSELECT, ROTOCCUW, ROTCUTOFF 
! END OF ROTAMER MOVE PARAMETERS

! GROUP ROTATION MOVE PARAMETERS
!      INTEGER :: GROUPROTFREQ, NGROUPS
!      LOGICAL :: GROUPROTT, DOGROUPROT=.FALSE.
!      CHARACTER(LEN=10), ALLOCATABLE :: ATOMGROUPNAMES(:)
!      INTEGER, ALLOCATABLE :: ATOMGROUPAXIS(:,:)
!      DOUBLE PRECISION, ALLOCATABLE :: ATOMGROUPSCALING(:),ATOMGROUPPSELECT(:)
!      LOGICAL, ALLOCATABLE :: ATOMGROUPS(:,:)
      
! END OF GROUP ROTATION MOVE PARAMETERS

!      INTEGER   NATOM, irest, ntb
!      DOUBLE PRECISION    tt,crd(3*656),vel(3,656)
!      CHARACTER(20) AMBERPRMTOP 
!      DOUBLE PRECISION xbar, ybar, zbar
!   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE        :: y
   double precision,  dimension(:), allocatable :: x
   integer, dimension(:), allocatable :: ipairs
   integer, dimension(:), allocatable, target :: ix
   character(len=4), dimension(:), allocatable :: ih
   integer :: myunitnew,myunit2,mdcrd_unit,mdinfo_unit,ambpdb_unit,ambxyz_unit,ambrst_unit,ambfinalio_node
logical master
common/extra_logical/master


integer iredir(8)
common/nmrrdr/redir,iredir

integer nmropt,iprint,noeskp,iscale,ipnlty,iuse,maxsub,jar
double precision scalm,pencut,ensave,tausw,ebdev,eadev,drjar
common/nmr1amber/scalm,pencut,ensave,tausw,ebdev,eadev,drjar, &
      nmropt,iprint,noeskp,iscale,ipnlty,iuse,maxsub,jar

integer     ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave
common/hulp/ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave

!  parameters for LES:

integer maxles,maxlestyp,maxlesadj
parameter (maxles=500000)
parameter (maxlestyp=100)
parameter (maxlesadj=100000)

integer bc_lesr,bc_lesi
parameter( bc_lesi=1+maxles*3+maxlesadj*2+1)
parameter (bc_lesr=maxlestyp*maxlestyp+1)

double precision lesfac(maxlestyp*maxlestyp),lfac

! for separate LES and non-LES temperature coupling

double precision ekinles0,temp0les,rndfles,sdfacles
double precision scaltles,tempsules,ekeles,rsdles
double precision ekmhles,ekphles

common/lesr/lesfac,temp0les

logical belly, erstop, qsetup, qpsander

logical newstyle,ok
integer BC_MDI  ! size in integers of common block mdi
integer BC_MDR  ! size in Reals of common block mdr

      integer ISGSTA,ISGEND,IPS,NNBIPST,NNBIPS
      COMMON/DTSGI/ISGSTA,ISGEND,IPS,NNBIPST,NNBIPS

      LOGICAL TSGLD,TLANGV,TEIPS,TVIPS
      COMMON/DTSGL/TSGLD,TLANGV,TEIPS,TVIPS

   ! runmin/trajene var
   double precision carrms


!+ Specification and control of Amber's working precision


double precision box,cut,scnb,scee,dielc,rad,wel,radhb,welhb, &
      cutcap,xcap,ycap,zcap,fcap,rwell,xbox0,ybox0,zbox0
common/boxr/box(3),cut,scnb,scee,dielc,xbox0,ybox0,zbox0, &
      cutcap,xcap,ycap,zcap,fcap,rwell, &
      rad(100),wel(100),radhb(200),welhb(200)

! ... integers:

integer ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp
common/boxi/ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp

double precision extraboxdim
!parameter (extraboxdim=30.d0)

! ... floats:

double precision t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, & !9
      tol,taur,dx0,drms,vlimit,rbtarg(9),tmass,tmassinv,  & !25
          kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt,  & !32
          gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon,  & !38
          solvph,rgbmax,fsmax,restraint_wt, &  !42
          skmin,skmax,vfac,gbneckscale,v11,v12,v22,kevb,evbt,Arad   !52
parameter (BC_MDR=52)
common/mdr/t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, &
      tol,taur,dx0,drms,vlimit,rbtarg,tmass,tmassinv, &
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt, &
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon, &
      solvph,rgbmax,fsmax,restraint_wt,skmin,skmax,vfac,gbneckscale, &
      v11,v12,v22,kevb,evbt,Arad

! ... strings:

character(len=4) iwtnm,iowtnm,ihwtnm
character(len=256) restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask
common/mds/ restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask,  &
            iwtnm,iowtnm,ihwtnm(2)

integer nrp,nspm,ig,ntx,ntcx,           &!5
      ntxo,ntt,ntp,ntr,init,             &!10
      ntcm,nscm,isolvp,nsolut,klambda,   &!15
      ntc,ntcc,ntf,ntid,ntn,             &!20
      ntnb,nsnb,ndfmin,nstlim,nrc,       &!25
      ntrx,npscal,imin,maxcyc,ncyc,      &!30
      ntmin,irest,jfastw,                &!33
      ibgwat,ienwat,iorwat,              &!36
      iwatpr,nsolw,igb,alpb,iyammp,           &!41
      gbsa,vrand,iwrap,nrespa,irespa,nrespai,icfe,  &!48
      rbornstat,ivcap,iconstreff,        &!51
      neb,vv,tmode,ipol,iesp,ievb,nodeid,num_noshake,    &!59
      idecomp,icnstph,ntcnstph,maxdup,numexchg,repcrd,numwatkeep     !66
parameter (BC_MDI=66)

common/mdiamber/nrp,nspm,ig, &
      ntx,ntcx,ntxo,ntt,ntp,ntr,init,ntcm,nscm, &
      isolvp,nsolut,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,ndfmin, &
      nstlim,nrc,ntrx,npscal,imin,maxcyc,ncyc,ntmin, &
      irest,jfastw,ibgwat,ienwat,iorwat, &
      iwatpr,nsolw,igb,alpb,iyammp,gbsa,vrand,numexchg,repcrd,numwatkeep, &
      iwrap,nrespa,irespa,nrespai,icfe,rbornstat, &
      ivcap,iconstreff,idecomp,klambda,icnstph,ntcnstph,maxdup,neb,vv, &
          tmode,ipol,iesp,ievb,nodeid,num_noshake


character(len=4096) groupbuffer
character(len=256) mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmrf, mincor, &
      vecs, radii, freqe,redir(8),rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file

character owrite, facc
common /files/ groupbuffer, mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmrf, mincor, &
      vecs, radii, freqe, owrite, facc,rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file

integer ier
integer       natom,nres,nbonh,nbona,ntheth,ntheta,nphih, &
      nphia,nnb,ntypes,nconp,maxmem,nwdvar,nparm, &
      natc,nattgtfit,nattgtrms,ibelly,natbel,ishake,nmxrs, &
      mxsub,natyp,npdec,i02,i04,i06,i08,i10, &
      iibh,ijbh,iicbh,iiba,ijba,iicba, &
      i24,i26,i28,i30,i32,i34,i36,i38,i40, &
      i42,i44,i46,i48,i50,i52,i54,i56,i58,ibellygp, &
      icnstrgp,itgtfitgp,itgtrmsgp,i64,i65,i68, &
      i70,i72,i74,i76,i78,i79,i80,i82,i84,i86, &
      icpstinf,icpresst,icptrsct, icpptcnt, &
      l15,lwinv,lpol,lcrd,lforce,l36,lvel,lvel2,l45,l50, &
      lcrdr,l60,l65,lmass,l75,l80,l85,l90,l95,l96,l97,l98,l99,lfrctmp, &
      l105,l110,l115,l120,l125,l130,l135,l140,l145,l150, &
      l165,l170,l175,l180,l185,l186,l187,l188,l189,l190, &
      lcpcrg,lcpene, &
      m02,m04,m06,m08,m10,m12,m14,m16,m18,i01, &
      iifstwt,iifstwr,nrealb,nintb,nholb,npairb,lastr,lasti,lasth, &
      lastpr,nbper,ngper,ndper,ifpert,lpolp, ncopy, &
      imask1,imask2,numadjst,mxadjmsk,icphidx,icptpair,lfsg,lvsg,noshake

!  1        2         3         4         5      6     7     8      9      10
common/memory/ &
 natom   ,nres     ,nbonh    ,nbona   ,ntheth,ntheta,nphih ,                       & ! 7
 nphia   ,nnb      ,ntypes   ,nconp   ,maxmem,nwdvar,nparm ,                       & !14
 natc    ,nattgtfit,nattgtrms,ibelly  ,natbel,ishake,nmxrs ,                       & !21
 mxsub   ,natyp    ,npdec    ,i02     ,i04   ,i06   ,i08   ,i10   ,                & !29
 iibh    ,ijbh     ,iicbh    ,iiba    ,ijba  ,iicba ,                              & !35
 i24     ,i26      ,i28      ,i30     ,i32   ,i34   ,i36   ,i38   ,i40   ,         & !44
 i42     ,i44      ,i46      ,i48     ,i50   ,i52   ,i54   ,i56   ,i58   ,ibellygp,& !54
 icnstrgp,itgtfitgp,itgtrmsgp,i64     ,i65   ,i68   ,                              & !60
 i70     ,i72      ,i74      ,i76     ,i78   ,i79   ,i80   ,i82   ,                & !68
 i84     ,i86      ,                                                               & !70
 icpstinf,icpresst ,icptrsct ,icpptcnt,                                            & !74
 l15      ,lwinv   ,lpol     ,lcrd    ,lforce,l36   ,lvel  ,lvel2 ,                & !82
 l45     ,l50      ,                                                               & !84
 lcrdr   ,l60      ,l65      ,lmass   ,l75   ,l80   ,l85   ,l90   ,l95   ,l96   ,  & !94
 l97     ,l98      ,l99      ,lfrctmp  ,                                           & !98
 l105    ,l110     ,l115     ,l120    ,l125  ,l130  ,l135  ,l140  ,l145  ,l150  ,  & !108
 l165    ,l170     ,l175     ,l180    ,l185  ,l186  ,l187  ,l188  ,l189  ,l190  ,  & !118
 lcpcrg  ,lcpene   ,                                                               & !120
 m02     ,m04      ,m06      ,m08     ,m10   ,m12   ,m14   ,m16   ,m18   ,i01   ,  & !130
 iifstwt ,iifstwr  ,nrealb   ,nintb   ,nholb ,npairb,lastr ,lasti ,lasth ,         & !139
 lastpr  ,nbper   ,ngper     ,ndper   ,ifpert,lpolp ,ncopy ,                       & !146
 imask1  ,imask2   ,numadjst ,mxadjmsk,icphidx,icptpair,lfsg,lvsg,noshake            !155

integer verbose,netfrc,     ew_type,    vdwmeth, &
      periodic,  use_pme,    opt_infl,   ischrgd, fix_dip, &
      fix_quad,  mpoltype,   induced,    frameon, chngmask, &
      scaldip

common/ewcntrl/ &
      verbose,   netfrc,     ew_type,    vdwmeth,    &! 4
      periodic,  use_pme,    opt_infl,   ischrgd, fix_dip,    &!9
      fix_quad,  mpoltype,   induced,    frameon, chngmask,   &!14
      scaldip


!integer    mchrg,  mamass,  mrk,    mreq,   mtk, &
!      mteq,   mfk,     mfpk,   mqeq,   mpk, &
!      mpn,    mphase,  msolty, mcn1,   mcn2, &
!      masol,  mbsol,   mhbcut, mxref,  msf, &
!      momega, mgraph,  miac,   miblo,  mico, &
!      mlbres, mipres,  mibh,   mjbh,   micbh, &
!      miba,   mjba,    micba,  mith,   mjth, &
!      mkth,   micth,   mita,   mjta,   mkta, &
!      micta,  miph,    mjph,   mkph,   mlph, &
!      micph,  mipa,    mjpa,   mkpa,   mlpa, &
!      micpa,  minb,    msymbl, mitree, mgroup, &
!      migres, miar1,   mx,     mf,     mh, &
!      mnbel,  mxbel,   miar2,  mcscr,  mcval, &
!      mcvec,  mdd,     mb,     mroots, mvect, &
!      mxdir,  mwref,   migrp2, ma,     mwr, &
!      mwi,    mz,      mfv1,   miv1,   mhrad, &
!      mgam,   mwinv,   mkpvt,  mxinit, mcn114, &
!      mcn214, mjoin,   mrotat, mpol


!common /pointr/   mchrg,  mamass,  mrk,    mreq,   mtk, &
!      mteq,   mfk,     mfpk,   mqeq,   mpk, &
!      mpn,    mphase,  msolty, mcn1,   mcn2, &
!      masol,  mbsol,   mhbcut, mxref,  msf, &
!      momega, mgraph,  miac,   miblo,  mico, &
!      mlbres, mipres,  mibh,   mjbh,   micbh, &
!      miba,   mjba,    micba,  mith,   mjth, &
!      mkth,   micth,   mita,   mjta,   mkta, &
!      micta,  miph,    mjph,   mkph,   mlph, &
!      micph,  mipa,    mjpa,   mkpa,   mlpa, &
!      micpa,  minb,    msymbl, mitree, mgroup, &
!      migres, miar1,   mx,     mf,     mh, &
!      mnbel,  mxbel,   miar2,  mcscr,  mcval, &
!      mcvec,  mdd,     mb,     mroots, mvect, &
!      mxdir,  mwref,   migrp2, ma,     mwr, &
!      mwi,    mz,      mfv1,   miv1,   mhrad, &
!      mgam,   mwinv,   mkpvt,  mxinit, mcn114, &
!      mcn214, mjoin,   mrotat, mpol

   double precision ene(51)
   integer native,nr3,nr
   ! nmrcal vars
   double precision f,enmr,devdis,devang,devtor,ag,bg,cg
   integer numphi,nttyp,nhb
! Common block containing variables relating to nmr restraints.

integer       intreq,irlreq,lnmr01,inmr02,iprr,iprw
common/nmrstf/intreq,irlreq,lnmr01,inmr02,iprr,iprw


double precision tgtrmsd,tgtmdfrc
common/tmd_real/ tgtrmsd,tgtmdfrc

!logical calledonce      ! just call subroutine amberinterface once
!common/logicals/ calledonce

! sf344> extra variables added to turn on/off the continuous smoothing of non-bonded terms
!       (that is: force switching for the electrostatics, Stoddard-Ford for the van der Waals terms)

integer ifswitch, irespa2, nrespa2
double precision fswitchbeta
!double precision, dimension (:,:),allocatable :: dedrcut,ecut
common/extrasamberint/ifswitch,irespa2,nrespa2
common/extrasamberdp/fswitchbeta

!msb50
integer, dimension(:,:),allocatable:: IC_COORDS
character(len=3), dimension(:), allocatable:: IC_IMPROP
integer:: lenic
integer, dimension(:), allocatable :: NPHIPSI
integer, dimension(:), allocatable :: NOMEGAC
integer, dimension(:), allocatable :: NSIDECHAIN
integer, dimension(:), allocatable :: NCHIRAL
integer, dimension(:), allocatable :: PHIPSI
integer, dimension(:), allocatable :: OMEGAC
integer, dimension(:), allocatable :: TW_SIDECHAIN
integer, dimension(:), allocatable :: CHIRAL
logical, dimension(:), allocatable :: IS_SIDECHAIN
integer, dimension(:), allocatable :: NICTOT !no of residues per segment
integer:: nseg
logical:: tomegac
logical:: AMBERICT, AMBSTEPT, AMBPERTT,AMBIT, AMBOLDPERTT, AMBICDNEBT
double precision:: PERTHRESH !threshhold for dihedral perturbation
                    !perturb if difference (start-end) > perthresh
integer:: NTW_ANGLES !if ambpertt: how many angles can you really twist
double precision, dimension(:), allocatable :: TW_DIFFP
!need this for natinterns
integer :: NBOND
logical, dimension (:),allocatable ::AM_BACKBONE
integer, allocatable:: IBIB(:),JBJB(:) 
logical, allocatable :: PROCHIRALH(:) !see chiralhyd
integer, allocatable :: PROCHIRALCNT(:)
logical :: NOPERMPROCHIRAL = .FALSE.
logical :: INTMINPERMT = .FALSE.
!
! Variables for steered minimisation/grouping.
!
LOGICAL :: STEEREDMINT, LOCALSTEEREDMINT
INTEGER :: SMINATOMA, SMINATOMB
DOUBLE PRECISION :: SMINK, SMINKINC, SMINDISTSTART, SMINDISTFINISH, SMINKCURRENT
!   atom groups for steered minimisation
INTEGER :: NATOMSINA, NATOMSINB, NATOMSINC
INTEGER,ALLOCATABLE,DIMENSION(:) :: ATOMSINALIST, ATOMSINBLIST, ATOMSINCLIST  
LOGICAL,ALLOCATABLE,DIMENSION(:) :: ATOMSINALISTLOGICAL, ATOMSINBLISTLOGICAL, ATOMSINCLISTLOGICAL
!LOGICAL :: BHDEBUG
!DOUBLE PRECISION :: BHSTEPSIZE
DOUBLE PRECISION :: AMCHNMAX, AMCHNMIN, AMCHPMAX, AMCHPMIN
! freezing in AMBER to speed up endhess
LOGICAL, ALLOCATABLE,DIMENSION(:) :: FROZENAMBER
INTEGER :: n_amb_calls
! khs26> IGB number (0 = no igb, for others check AMBER docs)
INTEGER :: AMBER_IGB
! khs26> Components of the energy from AMBER
DOUBLE PRECISION :: E_POTENTIAL, E_ELEC, E_IGB, E_BOND, E_ANGLE, E_DIHEDRAL, E_VDW, E_14_VDW, E_14_ELEC

! khs26> Add some topology information (e.g. atoms, names etc.)
CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:)  :: AMBER_ATOM_NAMES
CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:)  :: AMBER_RES_NAMES
CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)  :: AMBER_ELEMENT
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)  :: AMBER_MASSES
INTEGER, ALLOCATABLE, DIMENSION(:)           :: RESIDUE_INDEX
INTEGER, ALLOCATABLE, DIMENSION(:, :)        :: BONDED_ATOMS

! Arrays containing AMBER 9 topology data
CHARACTER(LEN=4), ALLOCATABLE :: AMBER9_ATOM_NAMES(:)
DOUBLE PRECISION, ALLOCATABLE :: AMBER9_CHARGES(:)
DOUBLE PRECISION, ALLOCATABLE :: AMBER9_MASSES(:)
INTEGER, ALLOCATABLE          :: AMBER9_ATOM_TYPE_INDICES(:)
INTEGER, ALLOCATABLE          :: AMBER9_ATOM_TO_RES(:)
CHARACTER(LEN=4), ALLOCATABLE :: AMBER9_RESIDUE_LABELS(:)
INTEGER, ALLOCATABLE          :: AMBER9_RESIDUE_POINTERS(:)
INTEGER, ALLOCATABLE          :: AMBER9_BONDS_INC_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_BONDS_WITHOUT_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_BONDS_ALL(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_ANGLES_INC_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_ANGLES_WITHOUT_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_ANGLES_ALL(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_DIHEDRALS_INC_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_DIHEDRALS_WITHOUT_HYDROGEN(:, :)
INTEGER, ALLOCATABLE          :: AMBER9_DIHEDRALS_ALL(:, :)

END MODULE MODAMBER9 

