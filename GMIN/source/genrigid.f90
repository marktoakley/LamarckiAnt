MODULE GENRIGID

      USE COMMONS, ONLY : AMBERT, AMBER12T 

      INTEGER :: NRIGIDBODY, DEGFREEDOMS, MAXSITE, NRELAXRIGIDR, NRELAXRIGIDA
      INTEGER :: XNATOMS
      INTEGER, ALLOCATABLE :: NSITEPERBODY(:), REFVECTOR(:), RIGIDSINGLES(:)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: RIGIDGROUPS
      DOUBLE PRECISION, ALLOCATABLE :: RIGIDCOORDS(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SITESRIGIDBODY
      DOUBLE PRECISION, ALLOCATABLE :: GR_WEIGHTS(:) ! weights for com calculation, e.g. masses
      LOGICAL :: RIGIDINIT, ATOMRIGIDCOORDT, RELAXRIGIDT
      LOGICAL :: GENRIGIDT

      DOUBLE PRECISION, ALLOCATABLE :: IINVERSE(:,:,:)

      LOGICAL :: RIGIDOPTIMROTAT, FREEZERIGIDBODYT
      DOUBLE PRECISION :: OPTIMROTAVALUES(3)
      LOGICAL :: AACONVERGENCET
      INTEGER, ALLOCATABLE :: LRBNPERMGROUP(:), LRBNPERMSIZE(:,:), LRBPERMGROUP(:,:), LRBNSETS(:,:), LRBSETS(:,:,:)

!   vr274:  added lattice coordinates
!           if HAS_LATTICE_COORDS is true, the last two atoms are treated
!           as lattice coordintes and rigidcoords is in reduced lattice units
      LOGICAL HAS_LATTICE_COORDS

!   jdf43:  RIGIDISRIGID = logical array, size NATOMS, TRUE if atom is part of
!           RB
      LOGICAL, ALLOCATABLE :: RIGIDISRIGID(:)
!   jdf43:  FROZENRIGIDBODY = logical array, size NATOMS, TRUE if RB is frozen
      LOGICAL, ALLOCATABLE :: FROZENRIGIDBODY(:)

!-----------------------------------------------------------------------------------!
! NRIGIDBODY  = number of rigid bodies
! DEGFREEDOMS = number of degrees of freedom = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS
! MAXSITE     = maximum number of sites in a rigid body
! NRELAXRIGIDR = rigid body minimisation for this number of steps
! NRELAXRIGIDA = atom minimisation for this number of steps
! NSITEPERBODY= number of rigid body sites, no need to be the same for all bodies
! REFVECTOR   = reference vector for the atomistic to rigic coordinate transformation
! RIGIDSINGLES= list of atoms not in rigid bodies
! RIGIDGROUPS = list of atoms in rigid bodies, need a file called rbodyconfig
! RIGIDCOORDS = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS coordinates
! SITESRIGIDBODY = coordinates of the rigid body sites
! RIGIDINIT   = logical variable for generalised rigid body
! ATOMRIGIDCOORDT, .TRUE. = atom coords active, .FALSE. = rigid coords active, used in mylbfgs & potential
! GENRIGIDT = generalised rigid body takestep taken if .TRUE.
!-----------------------------------------------------------------------------------!


CONTAINS
! vr274> TODO: better initialization of SITESRIDIGBODY, why use maxsite?

!-------------------------------------------
! vr274> Initializes basic structures
! hast to be the first call to GENRIGID function in order to setup basic structures.
! After that, the array which defines the sites can be filled. Then GENRIGID_INITIALIZE
! completes the initialization of rigid bodies.
!-------------------------------------------
SUBROUTINE GENRIGID_ALLOCATE(NEW_NRIGIDBODY,NEW_MAXSITE)
  USE COMMONS, only: NATOMS, AMBERT, AMBER12T
! hk286
  USE MODAMBER9, ONLY : ATMASS1
  USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_MASSES

  IMPLICIT NONE
  INTEGER, intent(in) :: NEW_NRIGIDBODY, NEW_MAXSITE
  
  NRIGIDBODY = NEW_NRIGIDBODY
  MAXSITE = NEW_MAXSITE

  ! hk286 > Allocate NSITEPERBODY
  ALLOCATE (NSITEPERBODY(NRIGIDBODY))
  ALLOCATE (SITESRIGIDBODY(MAXSITE,3,NRIGIDBODY))
  ALLOCATE (RIGIDGROUPS(MAXSITE,NRIGIDBODY))
  ALLOCATE (REFVECTOR(NRIGIDBODY))
  ALLOCATE (GR_WEIGHTS(NATOMS))
  ALLOCATE (IINVERSE(NRIGIDBODY,3,3))
!jdf43>
  ALLOCATE (RIGIDISRIGID(NATOMS))
  RIGIDISRIGID=.FALSE.

   ! If ATMASS1 is allocated and using AMBER 9, use ATMASS1 to scale.
   ! If AMBER12_MASSES is allocated and using AMBER 12, use AMBER12_MASSES.
   ! Otherwise, use centre of coordinates.
   IF ( ALLOCATED(ATMASS1) .AND. AMBERT ) THEN
      ! AMBER 9 atom masses
      GR_WEIGHTS = ATMASS1
   ELSE IF ( ALLOCATED(AMBER12_MASSES) .AND. AMBER12T ) THEN
      ! AMBER 12 atom masses
      GR_WEIGHTS = AMBER12_MASSES
   ELSE
      ! Otherwise centre of coords
      GR_WEIGHTS(:) = 1.0D0
   END IF

END SUBROUTINE

!-------------------------------------------
! vr274> Setup rigid body stuff after site definitions are done
!-------------------------------------------
SUBROUTINE GENRIGID_INITIALISE(INICOORDS)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2, J3, DUMMY
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, MASS
  LOGICAL :: SATOMT, RTEST
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION INICOORDS(3*NATOMS)
  
  INTEGER          :: INFO
  INTEGER, PARAMETER :: LWORK = 1000000 ! the dimension is set arbitrarily
  INTEGER :: I, J
  DOUBLE PRECISION :: DR(3), KBLOCK(3,3), KBEGNV(3)
  DOUBLE PRECISION :: WORK(LWORK)

! hk286 - 6/2
  DOUBLE PRECISION :: DET

! vr275> initialize coordinates for rigid bodies
  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        SITESRIGIDBODY(J2,:,J1) = INICOORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO


  ! hk286 > determine number of degrees of freedom
  DEGFREEDOMS = 0
  DO J1 = 1, NRIGIDBODY
     DEGFREEDOMS = DEGFREEDOMS + NSITEPERBODY(J1)
  ENDDO
  DEGFREEDOMS = 6 * NRIGIDBODY + 3 * (NATOMS - DEGFREEDOMS)

! hk286 > Allocate further data 
  ALLOCATE (RIGIDSINGLES((DEGFREEDOMS/3 - 2 * NRIGIDBODY)))
  ALLOCATE (RIGIDCOORDS(DEGFREEDOMS))

  DUMMY = 0
  DO J1 = 1, NATOMS
     SATOMT = .TRUE.
     DO J2 = 1, NRIGIDBODY
        DO J3 = 1, NSITEPERBODY(J2)
           IF (J1 == RIGIDGROUPS(J3,J2)) SATOMT = .FALSE.
        ENDDO
     ENDDO
     IF (SATOMT) THEN
        IF (DUMMY.EQ.DEGFREEDOMS/3-2*NRIGIDBODY) THEN
           WRITE(*,*) "genrigid> Error. More free atoms than expected."
           WRITE(*,*) "Likely problem with rbodyconfig."
           STOP
        ENDIF
        DUMMY = DUMMY + 1
        RIGIDSINGLES(DUMMY) = J1
     ENDIF
  ENDDO
  
  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0d0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + SITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + SITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + SITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS  = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        SITESRIGIDBODY(J2,1,J1) = SITESRIGIDBODY(J2,1,J1) - XMASS
        SITESRIGIDBODY(J2,2,J1) = SITESRIGIDBODY(J2,2,J1) - YMASS
        SITESRIGIDBODY(J2,3,J1) = SITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO

  ENDDO

  DO J1 = 1, NRIGIDBODY
     IINVERSE(J1,:,:) = 0.0D0
  END DO

! hk286 - 4/2
  IF (AACONVERGENCET .EQV. .TRUE.) THEN
     DO J1 = 1, NRIGIDBODY
        KBLOCK(:,:) = 0.0D0
        DO J2 = 1, NSITEPERBODY(J1)        
           DR(:)  = SITESRIGIDBODY(J2,:,J1)
           DO I = 1, 3
              KBLOCK(I,I) = KBLOCK(I,I) + (DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
        CALL DSYEV('V','U',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
        CALL RBDET(KBLOCK, DET)
        IF (DET < 0.0D0) THEN
           KBLOCK(:,3) = -KBLOCK(:,3)
           CALL RBDET(KBLOCK, DET)
           IF (DET < 0.0D0) THEN
              PRINT *, "GENRIGID> BAD ALIGNMENT", J1
              STOP
           ENDIF
        ENDIF
        KBLOCK = TRANSPOSE(KBLOCK)
!        PRINT *, KBEGNV
        IINVERSE(J1,1,1) = 1.0D0/KBEGNV(1)
        IINVERSE(J1,2,2) = 1.0D0/KBEGNV(2)
        IINVERSE(J1,3,3) = 1.0D0/KBEGNV(3)
        DO J2 = 1, NSITEPERBODY(J1)
           SITESRIGIDBODY(J2,:,J1) = MATMUL(KBLOCK,SITESRIGIDBODY(J2,:,J1))
        ENDDO
     ENDDO
  ENDIF
!  PRINT *, SITESRIGIDBODY
  
! hk286 > make sure the two atoms used as reference for rigid bodies are suitable
! Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
  DO J1 = 1, NRIGIDBODY
     REFVECTOR(J1) = 1
     RTEST = .TRUE.
     DO WHILE (RTEST)
        RTEST = .FALSE.
        DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
           PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
           IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
              RTEST = .TRUE.
           ENDIF
        ENDDO
        PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
        IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
        IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
           RTEST = .TRUE.
        ENDIF
        IF (RTEST) THEN
           REFVECTOR(J1) = REFVECTOR(J1) + 1               
        ENDIF
     ENDDO
  ENDDO


! hk286 - new 7/1/12
  IF (RIGIDOPTIMROTAT .EQV. .TRUE.) THEN
     CALL ROTATEINITIALREF ()
  ENDIF

! list of useful checks
!  P(1) = 0.0001D0
!  P(2) = 0.0003D0
!  P(3) = 0.0002D0
!  CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .TRUE.)

  !PRINT *, "TEST"
  !ATOMRIGIDCOORDT = .FALSE.
  !AMBERT = .TRUE.
  !CALL TRANSFORMCTORIGID(COORDS, RIGIDCOORDS)
  !COORDS(1:DEGFREEDOMS,1) = RIGIDCOORDS(1:DEGFREEDOMS)
  !COORDS(DEGFREEDOMS+1:3*NATOMS,1) = 0.0D0
  !CALL POTENTIAL(COORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
  !PRINT *, GRAD
  !CALL TRANSFORMGRAD(GRAD,RIGIDCOORDS,GRADR)
  !PRINT*, ENERGY
 
  
!      PRINT *, REFVECTOR
!      CHECKS  
!      PRINT*, COORDS(:,1)
!      PRINT*, " "
!      CALL TRANSFORMCTORIGID(COORDS(:,1), RIGIDCOORDS)
!      PRINT*, RIGIDCOORDS
!      PRINT*, " "
!      CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY, COORDS(:,1), RIGIDCOORDS)
!      PRINT*, COORDS(:,1)
!      PRINT*, " "
!      READ *, DUMMY
!      PRINT*, NRIGIDBODY
!      PRINT*, NSITEPERBODY
!      PRINT*, NATOMS

!      DO J1 = 1, 100
!         CALL TAKESTEPGENRIGID ()
!         PRINT*, " "
!         PRINT*, COORDS(:,1)
!         CALL TRANSFORMCTORIGID()
!         CALL TRANSFORMRIGIDTOC(1,2)
!         PRINT*, COORDS(:,1)
!      ENDDO
!      NATOMS = 4
!      CALL NEWCAPSID(RIGIDCOORDS,GRADR,ENERGY, .TRUE.)
!      NATOMS = 12
!      PRINT*, ENERGY
!      PRINT*, GRADR(1), GRADR(2), GRADR(3), GRADR(4), GRADR(5), GRADR(6)
!      PRINT*, GRADR(7), GRADR(8), GRADR(9), GRADR(10), GRADR(11), GRADR(12)

!      CALL GENCAPSID(COORDS(:,1),GRAD,ENERGY,.TRUE.)
!      CALL TRANSFORMGRAD(GRAD,RIGIDCOORDS,GRADR)
!      PRINT*, ENERGY
!      PRINT*, GRADR
!      PRINT*, " "
!      PRINT*, " "
!      PRINT*, NRIGIDBODY
!      PRINT*, NSITEPERBODY
!      PRINT*, NATOMS

      !READ*, DUMMY
    
!      CALL POTENTIALLBFGS(RIGIDCOORDS,GRADR,ENERGY,.TRUE.,.FALSE.,NATOMS)
!      PRINT*, ENERGY
!      PRINT*, RIGIDCOORDS
!      PRINT*, GRADR
!      PRINT*, " "

      !CALL VIEWGENCAPSID()

      !READ *, DUMMY
END SUBROUTINE

!-----------------------------------------------------------
SUBROUTINE GENRIGID_READ_FROM_FILE ()
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  CHARACTER(LEN=10) CHECK1
  INTEGER :: J1, J2, DUMMY, iostatus
  DOUBLE PRECISION :: INICOORDS(3*NATOMS)

! hk286 > read atomistic coordinates
! hk286 > in future, no need for separate coordsinirigid
! hk286 > currently the input coords files vary for CHARMM, AMBER, and RIGID BODIES
  IF (NATOMS == 0) THEN
     PRINT *, "ERROR STOP NOW > During generalised rigid body initialisation NATOMS = 0"
     STOP
  ENDIF

! vr274 > by standard don't do lattice coordinates
  HAS_LATTICE_COORDS = .FALSE.

  OPEN(UNIT = 28, FILE = 'coordsinirigid', STATUS = 'OLD')
  DO J1 = 1, NATOMS
     READ(28, *) INICOORDS(3*J1-2), INICOORDS(3*J1-1), INICOORDS(3*J1)
  ENDDO
  CLOSE(UNIT = 28)

! hk286 > determine no of rigid bodies
  NRIGIDBODY=0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO
     READ(222,*,IOSTAT=iostatus) CHECK1
     IF (iostatus<0) THEN
        CLOSE(222)
        EXIT
     ELSE IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') then
        NRIGIDBODY=NRIGIDBODY+1
     ENDIF
  END DO
  CLOSE(222)
  
! hk286 > determine maximum no of rigid body sites
  MAXSITE = 0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1,DUMMY
     IF (MAXSITE < DUMMY) MAXSITE = DUMMY
     DO J2 = 1, DUMMY
        READ(222,*) CHECK1
     ENDDO
  ENDDO
  CLOSE(222)

!  vr274> Calling function for allocation to make more general (setup rigid bodies from code)
  CALL GENRIGID_ALLOCATE(NRIGIDBODY,MAXSITE)

! hk286 > initialise SITESRIGIDBODY, RIGIDGROUPS, RIGIDSINGLES 
  OPEN(UNIT=222,FILE='rbodyconfig',status='unknown')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1, NSITEPERBODY(J1)
     DO J2 = 1, NSITEPERBODY(J1)
        READ(222,*) RIGIDGROUPS(J2,J1)
! csw34> check to make sure the current atom is not already in a rigid body
        IF (RIGIDISRIGID(RIGIDGROUPS(J2,J1))) THEN
            PRINT *," genrigid> ERROR: atom ",RIGIDGROUPS(J2,J1)," is in multiple rigid bodies! Stopping."
            STOP
        ELSE       
! csw34> if not, flag the current atom
            RIGIDISRIGID(RIGIDGROUPS(J2,J1))=.TRUE.
        ENDIF
! vr274> Moved initialization of coordinates to GENRIGID_INITIALISE, here only read the setup
!        SITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY,1)
     ENDDO
  ENDDO
  CLOSE(222)

  CALL GENRIGID_INITIALISE(INICOORDS)
END SUBROUTINE GENRIGID_READ_FROM_FILE

!-----------------------------------------------------------

!-----------------------------------------------------------

SUBROUTINE TRANSFORMRIGIDTOC (CMIN, CMAX, XCOORDS, XRIGIDCOORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J5, J7, J9
  INTEGER :: CMIN, CMAX
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: COM(3) ! center of mass
  LOGICAL          :: GTEST !, ATOMTEST
  DOUBLE PRECISION :: MLATTICE(3,3)
  
  GTEST = .FALSE.

! vr274 > are there additional lattice coordinates? If yes, setup transformation matrix
  IF(HAS_LATTICE_COORDS) THEN
    CALL GET_LATTICE_MATRIX(XRIGIDCOORDS(DEGFREEDOMS-5:DEGFREEDOMS), MLATTICE)
  ELSE ! vr274 > otherwise identity matrix
    MLATTICE = 0D0
    MLATTICE(1,1)=1d0
    MLATTICE(2,2)=1D0
    MLATTICE(3,3)=1D0
  ENDIF


  ! hk286 > coord transformations for rigid bodies CMIN to CMAX
  DO J1 = CMIN, CMAX
     J5   = 3*J1
     J7   = 3*NRIGIDBODY + J5
     P(:) = XRIGIDCOORDS(J7-2:J7)
     CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! vr274 > MLATTICE can have lattice transformation or be identity matrix
     COM = matmul(MLATTICE, XRIGIDCOORDS(J5-2:J5))
     DO J2 = 1,  NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        XCOORDS(3*J9-2:3*J9) = COM + MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
     ENDDO
     
  ENDDO
  
! hk286 > now the single atoms
! vr274 > this copies lattice coordinates as well which is stored in last 2 atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        J9 = RIGIDSINGLES(J1)
        XCOORDS(3*J9-2:3*J9) = XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1)
     ENDDO
  ENDIF
      
END SUBROUTINE TRANSFORMRIGIDTOC

!----------------------------------------------------------

SUBROUTINE ROTATEINITIALREF ()
IMPLICIT NONE
DOUBLE PRECISION :: P(3)
INTEGER J1

! hk286 - rotate the system - new
  P(:) = OPTIMROTAVALUES(:)
!  P(1) = -1.0D0 * 8.0D0 * ATAN(1.0D0) 
!  P(2) = 0.0D0 !4.0D0 * ATAN(1.0D0) !-(8*ATAN(1.0D0) - 5.0D0)/DSQRT(2.0D0)
!  P(3) = 0.0D0 !4.0D0 * ATAN(1.0D0)
  DO J1 = 1, NRIGIDBODY
     CALL REDEFINERIGIDREF (J1,P)
  ENDDO

END SUBROUTINE ROTATEINITIALREF

!----------------------------------------------------------

SUBROUTINE REDEFINERIGIDREF (J1,P)

  IMPLICIT NONE
  
  INTEGER :: J1, J2     !No of processor
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

!  PRINT *, "REDEFINE ", J1
  CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)  
  DO J2 = 1, NSITEPERBODY(J1)
     SITESRIGIDBODY(J2,:,J1) = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
  ENDDO

END SUBROUTINE REDEFINERIGIDREF

!----------------------------------------------------------

SUBROUTINE TRANSFORMCTORIGID (XCOORDS, XRIGIDCOORDS)
  USE COMMONS, ONLY: NATOMS, PERMDIST, MYUNIT
  USE VEC3
  USE ROTATIONS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J9     !No of processor
  DOUBLE PRECISION :: P(3)
  DOUBLE PRECISION :: COM(3), PNORM, PT(3,3), PI(3,3), MASS
  DOUBLE PRECISION :: XRIGIDCOORDS (DEGFREEDOMS), XCOORDS(3*NATOMS)

! vr274 > lattice matrix and inverse
  DOUBLE PRECISION MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

! hk286 - extra variables for minpermdist
  DOUBLE PRECISION :: D, DIST2, RMAT(3,3) 
  DOUBLE PRECISION :: PP1(3*NATOMS), PP2(3*NATOMS)
  LOGICAL :: TEMPPERMDIST

! vr274 > if has lattice coordinates, setup matrices
  IF(HAS_LATTICE_COORDS) THEN
    NLATTICECOORDS=6
    CALL GET_LATTICE_MATRIX(XCOORDS(3*NATOMS-5:3*NATOMS),MLATTICE)
  ELSE
    NLATTICECOORDS=0
    MLATTICE=0
    MLATTICE(1,1)=1
    MLATTICE(2,2)=1
    MLATTICE(3,3)=1
  ENDIF
  CALL INVERT3X3(MLATTICE, MLATTICEINV)

! loop over all rigid bodies
  DO J1 = 1, NRIGIDBODY
     COM = 0.0D0
     MASS = 0.0D0
     ! calculate center of mass
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        COM = COM + XCOORDS(3*J9-2:3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     COM = COM / MASS
     XRIGIDCOORDS(3*J1-2:3*J1) = COM

     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        PP1(3*J2-2:3*J2) = XCOORDS(3*J9-2:3*J9) - COM
        PP2(3*J2-2:3*J2) = SITESRIGIDBODY(J2,:,J1)
     ENDDO
     TEMPPERMDIST = PERMDIST
     PERMDIST = .FALSE.
     CALL MINPERMDIST(PP1(1:3*NSITEPERBODY(J1)),PP2(1:3*NSITEPERBODY(J1)),NSITEPERBODY(J1),.FALSE., &
          1.0D0,1.0D0,1.0D0,.FALSE.,.FALSE.,D,DIST2,.FALSE.,RMAT)
     PERMDIST = TEMPPERMDIST
     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = rot_mx2aa(RMAT)

     IF ( D/NSITEPERBODY(J1) > 0.1D0 ) THEN
        WRITE(MYUNIT, '(A, I3)')  'Warning: Genrigid > mapping looks bad for RB no ', J1 
        WRITE(MYUNIT, '(A)')  'Warning: Genrigid >  Often it is the permutation of the RB members, e.g. Hs in NH3'
     ENDIF

  ENDDO

! vr274> now translate everything to reduced units
  DO J1 = 1, NRIGIDBODY
    XRIGIDCOORDS(3*J1-2:3*J1) = MATMUL(MLATTICEINV, XRIGIDCOORDS(3*J1-2:3*J1))
  END DO

! hk286 > now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        ! vr274 > added lattice stuff
        XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = MATMUL(MLATTICEINV, XCOORDS(3*J9-2:3*J9))
     ENDDO
  ENDIF

! vr274 > copy lattice coords
  IF(HAS_LATTICE_COORDS) THEN
    XRIGIDCOORDS(DEGFREEDOMS - 5:DEGFREEDOMS) =  XCOORDS(3*NATOMS-5:3*NATOMS)
  ENDIF
END SUBROUTINE TRANSFORMCTORIGID

!-----------------------------------------------------------

SUBROUTINE TRANSFORMCTORIGID_OLD (XCOORDS, XRIGIDCOORDS)

  USE COMMONS, ONLY: NATOMS
  USE VEC3
  IMPLICIT NONE

  INTEGER :: J1, J2, J9
  DOUBLE PRECISION :: P(3), RMI(3,3), RMITEMP(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: TEMPP(3)
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, PT(3,3), PI(3,3), MASS
  DOUBLE PRECISION :: QR1, QR2, QR3, PCROSS(3), PCROSS21(3), PCROSS22(3), P2R(3)
  DOUBLE PRECISION :: PTEMP(3), PTEMP1(3), PTEMP2(3), DUMMY, DUMMY2
  DOUBLE PRECISION :: XRIGIDCOORDS (DEGFREEDOMS), XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: THETA, THETA2, PID, FRPISQ
  LOGICAL          :: GTEST, REDEFINET

! vr274 > lattice matrix and inverse
  DOUBLE PRECISION MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

  PID        = 4.D0*ATAN(1.0D0)
  FRPISQ     = 4.D0*PID*PID
  GTEST = .FALSE.

! vr274 > if has lattice coordinates, setup matrices
  IF(HAS_LATTICE_COORDS) THEN
    NLATTICECOORDS=6
    CALL GET_LATTICE_MATRIX(XCOORDS(3*NATOMS-5:3*NATOMS),MLATTICE)
  ELSE
    NLATTICECOORDS=0
    MLATTICE=0
    MLATTICE(1,1)=1
    MLATTICE(2,2)=1
    MLATTICE(3,3)=1
  ENDIF
  CALL INVERT3X3(MLATTICE, MLATTICEINV)

  DO J1 = 1, NRIGIDBODY

     REDEFINET = .FALSE.
10   IF (REDEFINET) THEN
        PRINT *, "REDEFINE CALLED"
!        TEMPP = (/0.2D0,0.2D0,0.2D0/)
        CALL REDEFINERIGIDREF(J1,OPTIMROTAVALUES)
        REDEFINET = .FALSE.
     ENDIF

! hk286 > compute the rotated reference vectors, save in PT(:,:)
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        IF( (J2-REFVECTOR(J1) .GE. 0) .AND. (J2-REFVECTOR(J1) .LE. 2) ) THEN
           PT(J2-REFVECTOR(J1)+1,:) = XCOORDS(3*J9-2:3*J9);
        ENDIF
        XMASS = XMASS + XCOORDS(3*J9-2)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XCOORDS(3*J9-1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XCOORDS(3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     XRIGIDCOORDS(3*J1-2) = XMASS
     XRIGIDCOORDS(3*J1-1) = YMASS
     XRIGIDCOORDS(3*J1)   = ZMASS         
     DO J2 = 1, 3
        PT(J2,1) = PT(J2,1) - XMASS
        PT(J2,2) = PT(J2,2) - YMASS
        PT(J2,3) = PT(J2,3) - ZMASS
     ENDDO
     DO J2 = 1, 3
        PNORM   = DSQRT(DOT_PRODUCT(PT(J2,:),PT(J2,:)))
        PT(J2,1) = PT(J2,1) / PNORM
        PT(J2,2) = PT(J2,2) / PNORM
        PT(J2,3) = PT(J2,3) / PNORM
        
! hk286 > the original, unrotated reference vectors PI
        PNORM   = DSQRT(DOT_PRODUCT(SITESRIGIDBODY(J2+REFVECTOR(J1)-1,:,J1),SITESRIGIDBODY(J2+REFVECTOR(J1)-1,:,J1)))
        PI(J2,1) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,1,J1) / PNORM
        PI(J2,2) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,2,J1) / PNORM
        PI(J2,3) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,3,J1) / PNORM           
     ENDDO

! hk286 > compute rotation around axis perpendicular to PI(1,:) and PT(1,:)
! hk286 > QR1 rotation around PCROSS as the axis 
     QR1 = DOT_PRODUCT(PI(1,:),PT(1,:))
     IF ( ((QR1 < 0.01D0 .AND. QR1 > -0.01D0) .OR. (QR1 < 2.0D0*PID + 0.01D0 .AND. QR1 > 2.0D0 * PID -0.01D0 )) &
          & .AND. (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF
     IF (QR1 > 1.0D0) QR1 = 1.0D0
     IF (QR1 < -1.0D0) QR1 = -1.0D0
     QR1 = ACOS(QR1)
     PCROSS(1) = PI(1,2)*PT(1,3) - PT(1,2)*PI(1,3)
     PCROSS(2) = PI(1,3)*PT(1,1) - PT(1,3)*PI(1,1)
     PCROSS(3) = PI(1,1)*PT(1,2) - PT(1,1)*PI(1,2)
     PNORM = DSQRT(DOT_PRODUCT(PCROSS,PCROSS))
     IF ((QR1 < 1D-6 .AND. QR1 > -1D-6).OR.(PNORM < 1D-6)) THEN
        QR1 = QR1 + 4.0D0*ASIN(1.0D0)
        PCROSS(:) = PT(1,:) 
     ELSE
        PCROSS(1) = PCROSS(1) / PNORM
        PCROSS(2) = PCROSS(2) / PNORM
        PCROSS(3) = PCROSS(3) / PNORM
     ENDIF
     CALL RMDRVT(QR1*PCROSS, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! hk286 > any additional rotation not accounted by the above operation must be
! hk286 > rotation around PT(1,:), now compute this rotation
! hk286 > rotate the second reference vector PI(2,:) around PCROSS by QR1
     P2R = MATMUL(RMI(:,:),PI(2,:))
! hk286 > take suitable projections
! hk286 > QR2 is the amount of rotation around PT(1,:)
     PCROSS21(1) = P2R(2)*PT(1,3) - PT(1,2)*P2R(3)
     PCROSS21(2) = P2R(3)*PT(1,1) - PT(1,3)*P2R(1)
     PCROSS21(3) = P2R(1)*PT(1,2) - PT(1,1)*P2R(2)
     PCROSS22(1) = PT(2,2)*PT(1,3) - PT(1,2)*PT(2,3)
     PCROSS22(2) = PT(2,3)*PT(1,1) - PT(1,3)*PT(2,1)
     PCROSS22(3) = PT(2,1)*PT(1,2) - PT(1,1)*PT(2,2)

     PNORM = DSQRT(DOT_PRODUCT(PCROSS21,PCROSS21))
     PCROSS21 = PCROSS21 / PNORM
     PNORM = DSQRT(DOT_PRODUCT(PCROSS22,PCROSS22))
     PCROSS22 = PCROSS22 / PNORM
     QR2 = DOT_PRODUCT(PCROSS21,PCROSS22)
     IF ( ((QR2 < 0.01D0 .AND. QR2 > -0.01D0) .OR. (QR2 < 2.0D0*PID + 0.01D0 .AND. QR2 > 2.0D0 * PID -0.01D0)) &
          & .AND. (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF     
     IF (QR2 > 1.0D0) QR2 = 1.0D0
     IF (QR2 < -1.0D0) QR2 = -1.0D0
     QR2 = ACOS(QR2)
     CALL RMDRVT(QR2*PT(1,:), RMI, DRMI1, DRMI2, DRMI3, GTEST)
     PTEMP = MATMUL(RMI(:,:),P2R(:))
     DUMMY = (PTEMP(1)-PT(2,1))**2 + (PTEMP(2)-PT(2,2))**2 + (PTEMP(3)-PT(2,3))**2
     CALL RMDRVT((2.0D0*PID-QR2)*PT(1,:), RMI, DRMI1, DRMI2, DRMI3, GTEST)
     PTEMP = MATMUL(RMI(:,:),P2R(:))
     DUMMY2 = (PTEMP(1)-PT(2,1))**2 + (PTEMP(2)-PT(2,2))**2 + (PTEMP(3)-PT(2,3))**2
     IF (DUMMY2 < DUMMY) THEN
        QR2 = 2.0D0 * PID - QR2
     ENDIF

! hk286 > Construct quarternions
!         Quarternion1 = [COS(QR1/2) PCROSS(:)*SIN(QR1/2)]
!         Quarternion2 = [COS(QR2/2) PT(1,:)*SIN(QR2/2)]
     PTEMP1 = PCROSS(:)*SIN(QR1/2)
     PTEMP2 = PT(1,:)*SIN(QR2/2)
     QR3 = COS(QR1/2)*COS(QR2/2) - DOT_PRODUCT(PTEMP1,PTEMP2)
     IF (QR3 > 1.0D0) QR3 = 1.0D0
     IF (QR3 < -1.0D0) QR3 = -1.0D0    
     QR3 = 2*ACOS(QR3)
     IF ( ((QR3 < 0.01D0 .AND. QR3 > -0.01D0) .OR. (QR3 < 2.0D0*PID + 0.01D0 .AND. QR3 > 2.0D0 * PID -0.01D0)) &
          & .AND. (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF          
     P(1) = (PT(1,2)*PCROSS(3) - PCROSS(2)*PT(1,3))*SIN(QR1/2)*SIN(QR2/2)
     P(2) = (PT(1,3)*PCROSS(1) - PCROSS(3)*PT(1,1))*SIN(QR1/2)*SIN(QR2/2)
     P(3) = (PT(1,1)*PCROSS(2) - PCROSS(1)*PT(1,2))*SIN(QR1/2)*SIN(QR2/2)
     P(:) = P(:) + COS(QR1/2)*PT(1,:)*SIN(QR2/2) + COS(QR2/2)*PCROSS(:)*SIN(QR1/2)
     IF (abs(sin(QR3)) < 1D-6 ) THEN
        PNORM = SQRT(DOT_PRODUCT(P,P))
        P(:) = P(:) * 2.0D0
     ELSE
        P(:) = P(:) / SIN(QR3/2) * QR3
     ENDIF
     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = P(:)

  ENDDO
  
! vr274> now translate everything to reduced units
  DO J1 = 1, NRIGIDBODY
    XRIGIDCOORDS(3*J1-2:3*J1) = MATMUL(MLATTICEINV, XRIGIDCOORDS(3*J1-2:3*J1))
  END DO

! hk286 > now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        ! vr274 > added lattice stuff
        XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = MATMUL(MLATTICEINV, XCOORDS(3*J9-2:3*J9))
     ENDDO
  ENDIF

! vr274 > copy lattice coords
  IF(HAS_LATTICE_COORDS) THEN
    XRIGIDCOORDS(DEGFREEDOMS - 5:DEGFREEDOMS) =  XCOORDS(3*NATOMS-5:3*NATOMS)
  ENDIF


! hk286 > reduce the amount of rotation between 2*PI and 4*PI
  DO J1 = 1, NRIGIDBODY

     J2       = 3*NRIGIDBODY + 3*J1
     THETA2  = DOT_PRODUCT(XRIGIDCOORDS(J2-2:J2),XRIGIDCOORDS(J2-2:J2))
     THETA   = DSQRT(THETA2)
     IF (THETA2 > FRPISQ) THEN
        THETA2   = DSQRT(THETA2)
        THETA    = THETA2 - INT(THETA2/(2.D0*PID))*2.D0*PID
        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA2 * THETA
     ENDIF
!     IF (THETA > PID) THEN
!        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA * (-2.0D0*PID + THETA)
!        THETA = 2.0D0 * PID - THETA
!     ENDIF
     IF (RIGIDOPTIMROTAT .EQV. .TRUE.) THEN
        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA * (2.0D0*PID + THETA)
        THETA = 2.0D0*PID + THETA
     ENDIF

  ENDDO
  
END SUBROUTINE TRANSFORMCTORIGID_OLD

!-----------------------------------------------------------

SUBROUTINE TRANSFORMGRAD (G, XR, GR)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  DOUBLE PRECISION :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: DR1(3),DR2(3),DR3(3) 
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  LOGICAL :: GTEST
  INTEGER :: NLATTICECOORDS
  DOUBLE PRECISION :: MLATTICE(3,3)
  
  NLATTICECOORDS=0
  IF(HAS_LATTICE_COORDS) THEN
      NLATTICECOORDS=6
  ENDIF

  GTEST = .TRUE.
  GR(:) = 0.0D0
  
  DO J1 = 1, NRIGIDBODY
     
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)
     
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)

! hk286 > translation
        GR(3*J1-2:3*J1) = GR(3*J1-2:3*J1) + G(3*J9-2:3*J9)
        
! hk286 > rotation
        DR1(:) = MATMUL(DRMI1,SITESRIGIDBODY(J2,:,J1))
        DR2(:) = MATMUL(DRMI2,SITESRIGIDBODY(J2,:,J1))
        DR3(:) = MATMUL(DRMI3,SITESRIGIDBODY(J2,:,J1))
        GR(3*NRIGIDBODY+3*J1-2) = GR(3*NRIGIDBODY+3*J1-2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        GR(3*NRIGIDBODY+3*J1-1) = GR(3*NRIGIDBODY+3*J1-1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        GR(3*NRIGIDBODY+3*J1)   = GR(3*NRIGIDBODY+3*J1)   + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
  ENDDO

! hk286 - testing 6/6/12
  IF (FREEZERIGIDBODYT .EQV. .TRUE.) THEN
! jdf43
     DO J1=1,NRIGIDBODY
        IF (FROZENRIGIDBODY(J1)) THEN
!           GR(3*J1-2:3*J1) = 0.0D0
!           GR(DEGFREEDOMS/2+3*J1-2:DEGFREEDOMS/2+3*J1) = 0.0D0
            CONTINUE
        ENDIF
     ENDDO
!     GR(3*NRIGIDBODY-2:3*NRIGIDBODY) = 0.0D0
!     GR(6*NRIGIDBODY-2:6*NRIGIDBODY) = 0.0D0
  ENDIF

! hk286 > single atoms
! vr274 > and lattice
  IF (DEGFREEDOMS > 6 * NRIGIDBODY - NLATTICECOORDS) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = G(3*J9-2:3*J9)
     ENDDO
  ENDIF

  IF(HAS_LATTICE_COORDS) THEN
      CALL GET_LATTICE_MATRIX(XR(DEGFREEDOMS-5:DEGFREEDOMS),MLATTICE)

      ! vr274> for lattice, go to reduced coordinates
      DO J1 = 1, NRIGIDBODY
          GR(3*J1-2:3*J1) =  matmul(transpose(mlattice), GR(3*J1-2:3*J1))
      ENDDO
      ! vr274> and single atoms
      IF (DEGFREEDOMS > 6 * NRIGIDBODY + NLATTICECOORDS) THEN
          DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
              J2 = 6*NRIGIDBODY + 3*J1
              GR(J2-2:J2) = matmul(transpose(mlattice), GR(J2-2:J2))
          ENDDO
      ENDIF
      ! copy lattice gradient
      GR(DEGFREEDOMS-5:DEGFREEDOMS) = G(3*NATOMS-5:3*NATOMS)
  ENDIF

END SUBROUTINE TRANSFORMGRAD

!-----------------------------------------------------------

!SUBROUTINE BACKTRANSFORMGRAD (G, XR, GR)
  
!  USE COMMONS, ONLY: NATOMS
!  USE MEASURE
!  IMPLICIT NONE
  
!  INTEGER          :: J1, J2, J9
!  DOUBLE PRECISION :: G(3*NATOMS), G2(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
!  DOUBLE PRECISION :: PI(3), RC(3)
!  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
!  DOUBLE PRECISION :: RMI2(3,3)
!  DOUBLE PRECISION :: TORQUE(3), OMEGA(3), OMEGAMAT(3,3)

!  G2(:) = G(:)
!  G(:) = 0.0D0

!  DO J1 = 1, NRIGIDBODY
     
!     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
!     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, .TRUE.)

!     TORQUE(:) = 0.0D0
!     DO J2 = 1, NSITEPERBODY(J1)
!        RC = MATMUL(RMI, SITESRIGIDBODY(J2,:,J1))
!        J9 = RIGIDGROUPS(J2, J1)
!        TORQUE(:) = TORQUE(:) + CROSS_PRODUCT(RC,G2(3*J9-2:3*J9))
!     ENDDO

!     RMI2 = MATMUL( RMI, IINVERSE(J1,:,:))
!     RMI2 = MATMUL( RMI2, TRANSPOSE(RMI) )
!     OMEGA(:) = MATMUL ( RMI2, TORQUE ) 
!     OMEGAMAT(:,:) = 0.0D0
!     OMEGAMAT(1,2) = -OMEGA(3)
!     OMEGAMAT(1,3) =  OMEGA(2)
!     OMEGAMAT(2,1) =  OMEGA(3)
!     OMEGAMAT(2,3) = -OMEGA(1)
!     OMEGAMAT(3,1) = -OMEGA(2)
!     OMEGAMAT(3,2) =  OMEGA(1)

!     RMI2 = MATMUL (OMEGAMAT, RMI)

!     DO J2 = 1, NSITEPERBODY(J1)
!        J9 = RIGIDGROUPS(J2, J1)
!        G(3*J9-2:3*J9) = MATMUL(RMI2, SITESRIGIDBODY(J2,:,J1))
!        G(3*J9-2:3*J9) = G(3*J9-2:3*J9) + GR(3*J1-2:3*J1) / NSITEPERBODY(J1)
!     ENDDO
     
!  ENDDO

!  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
!     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
!        J9 = RIGIDSINGLES(J1)
!        G(3*J9-2:3*J9) = GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1)
!     ENDDO
!  ENDIF

!END SUBROUTINE BACKTRANSFORMGRAD


!-----------------------------------------------------------

SUBROUTINE AACONVERGENCE (G, XR, GR, RMS)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  DOUBLE PRECISION :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: TORQUE(3)
  DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  DOUBLE PRECISION :: DR1(3),DR2(3),DR3(3), RMI3(3,3), RMS 

  RMS = 0.0D0
  PI = (/0.0D0, 0.0D0, 0.0D0/)
  CALL RMDRVT(PI, RMI0, DRMI10, DRMI20, DRMI30, .TRUE.)

  DO J1 = 1, NRIGIDBODY
     
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)

     TORQUE(:) = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        DR1(:) = MATMUL(DRMI10,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR2(:) = MATMUL(DRMI20,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR3(:) = MATMUL(DRMI30,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        TORQUE(1) = TORQUE(1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        TORQUE(2) = TORQUE(2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        TORQUE(3) = TORQUE(3) + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
     TORQUE = MATMUL(TRANSPOSE(RMI), TORQUE)
     RMS = RMS + DOT_PRODUCT(TORQUE, MATMUL(TRANSPOSE(IINVERSE(J1,:,:)),TORQUE))
!     RMI3 = MATMUL( RMI, IINVERSE(J1,:,:))
!     RMI3 = MATMUL( RMI3, TRANSPOSE(RMI) )
!     RMS = RMS + DOT_PRODUCT(TORQUE, MATMUL(TRANSPOSE(RMI3),TORQUE))
     RMS = RMS + 1.0D0/NSITEPERBODY(J1) * DOT_PRODUCT(GR(3*J1-2:3*J1),GR(3*J1-2:3*J1)) 
  ENDDO

  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        RMS = RMS + DOT_PRODUCT(GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1),GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1))
     ENDDO
  ENDIF

  RMS=MAX(DSQRT(RMS/(3*NATOMS)),1.0D-100)

END SUBROUTINE AACONVERGENCE

!--------------------------------------------------------------

! hk286 > Often we want to check if the atoms grouped in a rigid body has moved or not
! hk286 > They should not if everything is done correctly
! hk286 > REDEFINESITEST = .FALSE. then it prints to standard output
! hk286 > REDEFINESITEST = .TRUE. then regroup atoms, SITESRIGIDBODY rewritten

SUBROUTINE CHECKSITES (REDEFINESITEST, COORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  INTEGER :: J1, J2, J3, DUMMY
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, MASS
  DOUBLE PRECISION :: XSITESRIGIDBODY(MAXSITE,3,NRIGIDBODY)
  DOUBLE PRECISION :: COORDS(3*NATOMS)
  LOGICAL :: RTEST, REDEFINESITEST
  

  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY = RIGIDGROUPS(J2,J1)
        XSITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + XSITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XSITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XSITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        XSITESRIGIDBODY(J2,1,J1) = XSITESRIGIDBODY(J2,1,J1) - XMASS
        XSITESRIGIDBODY(J2,2,J1) = XSITESRIGIDBODY(J2,2,J1) - YMASS
        XSITESRIGIDBODY(J2,3,J1) = XSITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO
  ENDDO
  

  IF (REDEFINESITEST) THEN
!     PRINT *, " SITES REDEFINED "
     SITESRIGIDBODY(:,:,:) = XSITESRIGIDBODY(:,:,:)

!Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
     DO J1 = 1, NRIGIDBODY
        REFVECTOR(J1) = 1
        RTEST = .TRUE.
        DO WHILE (RTEST)
           RTEST = .FALSE.
           DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
              PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
              IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
                 RTEST = .TRUE.
              ENDIF
           ENDDO
           PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
           IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
           IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
              RTEST = .TRUE.
           ENDIF
           IF (RTEST) THEN
              REFVECTOR(J1) = REFVECTOR(J1) + 1               
           ENDIF
        ENDDO
     ENDDO
  ELSE
!     PRINT *, XSITESRIGIDBODY
  ENDIF

END SUBROUTINE CHECKSITES

!--------------------------------------------------------------

! vr274 > build the lattice matrix.
!         The matrix is an triangular matrix,  c vector is always perpendicular to z
SUBROUTINE GET_LATTICE_MATRIX(LATTICE_COORDS, M)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: LATTICE_COORDS(6)
    DOUBLE PRECISION, INTENT(OUT) :: M(3,3)
    M=0
    M(1,1) = LATTICE_COORDS(1)
    M(2,1) = LATTICE_COORDS(2)
    M(3,1) = LATTICE_COORDS(3)
    M(2,2) = LATTICE_COORDS(4)
    M(3,2) = LATTICE_COORDS(5)
    M(3,3) = LATTICE_COORDS(6)
END SUBROUTINE

! vr274 > set lattice coordinates from lattice matrix.
!         The matrix is an triangular matrix,  c vector is always perpendicular to z
SUBROUTINE SET_LATTICE_MATRIX(LATTICE_COORDS, M)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: LATTICE_COORDS(6)
    DOUBLE PRECISION, INTENT(IN) :: M(3,3)

    LATTICE_COORDS(1) = M(1,1)
    LATTICE_COORDS(2) = M(2,1)
    LATTICE_COORDS(3) = M(3,1)
    LATTICE_COORDS(4) = M(2,2)
    LATTICE_COORDS(5) = M(3,2)
    LATTICE_COORDS(6) = M(3,3)
END SUBROUTINE

! csw34> expansion move for rigid bodies
SUBROUTINE GENRIGID_EXPAND(XCOORDS, EXPANDFACTOR)

USE COMMONS, ONLY: NATOMS, NORMALISEEXPANDT
IMPLICIT NONE

INTEGER :: J1, J2, J3  
DOUBLE PRECISION :: SYSTEMCOM(3), COM(3), NEWCOM(3), MASS, EXPANDVECTOR(3) 
DOUBLE PRECISION, INTENT(INOUT) :: XCOORDS(3*NATOMS) 
DOUBLE PRECISION, INTENT(IN) :: EXPANDFACTOR

SYSTEMCOM = 0.0D0

! Calculate the centre of coordinates for the system
DO J1 = 1,NATOMS
   SYSTEMCOM = SYSTEMCOM + XCOORDS(3*J1-2:3*J1)
ENDDO
SYSTEMCOM = SYSTEMCOM / NATOMS

! Loop over all rigid bodies
DO J1 = 1, NRIGIDBODY
   COM = 0.0D0
   NEWCOM = 0.0D0
   MASS = 0.0D0
   EXPANDVECTOR = 0.0D0

! For each rigid body, calculate center of mass
   DO J2 = 1, NSITEPERBODY(J1)
      J3 = RIGIDGROUPS(J2, J1)
      COM = COM + XCOORDS(3*J3-2:3*J3)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
      MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
   ENDDO
   COM = COM / MASS

! Calculate the expansion vector
   EXPANDVECTOR = COM - SYSTEMCOM
! If NORMALISE is specified as a second arguement, normalise the vector using VECNORM
   IF (NORMALISEEXPANDT) CALL VECNORM(EXPANDVECTOR,3)
! Scale the vector by EXPANDFACTOR
   EXPANDVECTOR = EXPANDVECTOR*EXPANDFACTOR  
! Apply the expansion vector to the atoms in the rigid body
   DO J2 = 1, NSITEPERBODY(J1)
      J3 = RIGIDGROUPS(J2, J1)
      XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) + EXPANDVECTOR
   ENDDO
ENDDO

END SUBROUTINE GENRIGID_EXPAND 

SUBROUTINE RBDET(A, DET)
  
  IMPLICIT NONE
  DOUBLE PRECISION :: A (3,3), DET

  DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)) 

END SUBROUTINE RBDET

SUBROUTINE INVERSEMATRIX(A, AINVERSE)
  
  IMPLICIT NONE
  DOUBLE PRECISION :: A (3,3), AINVERSE(3,3), DET

  DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)) 
  AINVERSE(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
  AINVERSE(1,2) = A(3,2)*A(1,3)-A(3,3)*A(1,2)
  AINVERSE(1,3) = A(2,3)*A(1,2)-A(2,2)*A(1,3)
  AINVERSE(2,1) = A(3,1)*A(2,3)-A(3,3)*A(2,1)
  AINVERSE(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
  AINVERSE(2,3) = A(2,1)*A(1,3)-A(2,3)*A(1,1)
  AINVERSE(3,1) = A(3,2)*A(2,1)-A(3,1)*A(2,2)
  AINVERSE(3,2) = A(3,1)*A(1,2)-A(3,2)*A(1,1)
  AINVERSE(3,3) = A(2,2)*A(1,1)-A(2,1)*A(1,2)
  AINVERSE(:,:) = AINVERSE(:,:)/DET

END SUBROUTINE INVERSEMATRIX

SUBROUTINE TRANSFORMHESSIAN (H, G, XR, HR, RBAANORMALMODET)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J3, J4, J8, J9, K, L
  DOUBLE PRECISION :: G(3*NATOMS), H(3*NATOMS,3*NATOMS)
  DOUBLE PRECISION :: XR(DEGFREEDOMS), HR(DEGFREEDOMS,DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: AD2R11(3),AD2R22(3),AD2R33(3),AD2R12(3),AD2R23(3),AD2R31(3) 
  DOUBLE PRECISION :: ADR1(3),ADR2(3),ADR3(3) 
  DOUBLE PRECISION :: ARMI(3,3), ADRMI1(3,3), ADRMI2(3,3), ADRMI3(3,3)
  DOUBLE PRECISION :: AD2RMI11(3,3), AD2RMI22(3,3), AD2RMI33(3,3)
  DOUBLE PRECISION :: AD2RMI12(3,3), AD2RMI23(3,3), AD2RMI31(3,3)
  DOUBLE PRECISION :: BDR1(3),BDR2(3),BDR3(3) 
  DOUBLE PRECISION :: BRMI(3,3), BDRMI1(3,3), BDRMI2(3,3), BDRMI3(3,3)
  DOUBLE PRECISION :: BD2RMI11(3,3), BD2RMI22(3,3), BD2RMI33(3,3)
  DOUBLE PRECISION :: BD2RMI12(3,3), BD2RMI23(3,3), BD2RMI31(3,3)
  LOGICAL :: GTEST, STEST, RBAANORMALMODET
  DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  DOUBLE PRECISION :: D2RMI10(3,3), D2RMI20(3,3), D2RMI30(3,3), D2RMI120(3,3), D2RMI230(3,3), D2RMI310(3,3)
  
  GTEST = .TRUE.
  STEST = .TRUE.
  HR(:,:) = 0.0D0

  IF ( RBAANORMALMODET ) THEN
     PI = (/0.0D0, 0.0D0, 0.0D0/)
     CALL RMDFAS(PI, RMI0, DRMI10, DRMI20, DRMI30, D2RMI10, D2RMI20, D2RMI30, D2RMI120, D2RMI230, D2RMI310, GTEST, STEST)
  ENDIF

  DO J1 = 1, NRIGIDBODY

     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)
    
     DO J2 = J1, NRIGIDBODY
        
        PI = XR(3*NRIGIDBODY+3*J2-2 : 3*NRIGIDBODY+3*J2)
        CALL RMDFAS(PI, BRMI, BDRMI1, BDRMI2, BDRMI3, BD2RMI11, BD2RMI22, BD2RMI33, BD2RMI12, BD2RMI23, BD2RMI31, GTEST, STEST)
        
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)

           DO J4 = 1, NSITEPERBODY(J2)
              J9 = RIGIDGROUPS(J4, J2)
                            
! hk286 > translation
              HR(3*J1-2:3*J1, 3*J2-2:3*J2) = HR(3*J1-2:3*J1, 3*J2-2:3*J2) + H(3*J8-2:3*J8, 3*J9-2:3*J9)

! hk286 > rotations
              IF ( RBAANORMALMODET ) THEN
                 ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 BDR1(:) = MATMUL(DRMI10,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR2(:) = MATMUL(DRMI20,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR3(:) = MATMUL(DRMI30,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
              ELSE
                 ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
                 ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
                 ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
                 BDR1(:) = MATMUL(BDRMI1,SITESRIGIDBODY(J4,:,J2))
                 BDR2(:) = MATMUL(BDRMI2,SITESRIGIDBODY(J4,:,J2))
                 BDR3(:) = MATMUL(BDRMI3,SITESRIGIDBODY(J4,:,J2))
              ENDIF

! hk286 - mixed translation rotation
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*J1  , 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*J1  , 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*J1-2, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*J1-1, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*J1  , 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR3(:))        
              
              IF (J2 > J1) THEN
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-2) = HR(3*J2  , 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-1) = HR(3*J2  , 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1  ) = HR(3*J2  , 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
              ENDIF
! hk286 - double rotation
              DO K = 1, 3
                 DO L = 1, 3
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR3(L)   
                 ENDDO
              ENDDO
           ENDDO
           IF (J1 .EQ. J2) THEN
              IF ( RBAANORMALMODET ) THEN
                 AD2R11(:) = MATMUL(D2RMI10, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R22(:) = MATMUL(D2RMI20, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R33(:) = MATMUL(D2RMI30, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R12(:) = MATMUL(D2RMI120,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R23(:) = MATMUL(D2RMI230,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R31(:) = MATMUL(D2RMI310,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ELSE
                 AD2R11(:) = MATMUL(AD2RMI11,SITESRIGIDBODY(J3,:,J1))
                 AD2R22(:) = MATMUL(AD2RMI22,SITESRIGIDBODY(J3,:,J1))
                 AD2R33(:) = MATMUL(AD2RMI33,SITESRIGIDBODY(J3,:,J1))
                 AD2R12(:) = MATMUL(AD2RMI12,SITESRIGIDBODY(J3,:,J1))
                 AD2R23(:) = MATMUL(AD2RMI23,SITESRIGIDBODY(J3,:,J1))
                 AD2R31(:) = MATMUL(AD2RMI31,SITESRIGIDBODY(J3,:,J1))
              ENDIF
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R11(:))
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R22(:))
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R33(:))
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R12(:))
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R23(:))
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R31(:))
           ENDIF
        ENDDO

        IF (J1 .EQ. J2) THEN
           HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1-1) 
           HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1  )
           HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-2)
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
        ELSE
           HR(3*J2-2:3*J2, 3*J1-2:3*J1) = TRANSPOSE(HR(3*J1-2:3*J1, 3*J2-2:3*J2))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J2-2:3*J2) = &
                TRANSPOSE(HR(3*J2-2:3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
                TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
        ENDIF

     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY

     DO J2 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
            
        J9 = RIGIDSINGLES(J2)
        PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
        CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)
        
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)
                     
! hk286 > translation
           HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) &
                + H(3*J8-2:3*J8, 3*J9-2:3*J9)

! hk286 > rotations
           IF ( RBAANORMALMODET ) THEN
              ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
           ELSE
              ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
              ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
              ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
           ENDIF
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
        ENDDO
        
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
             TRANSPOSE(HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
             TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
     ENDDO
  ENDDO

  DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J8 = RIGIDSINGLES(J1)
     DO J2 = J1, (DEGFREEDOMS - 6*NRIGIDBODY)/3            
        J9 = RIGIDSINGLES(J2)                     
        HR(6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = H(3*J8-2:3*J8, 3*J9-2:3*J9)
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1) = H(3*J9-2:3*J9, 3*J8-2:3*J8)
     ENDDO
  ENDDO

END SUBROUTINE TRANSFORMHESSIAN

! csw34> random rotation move for rigid bodies
SUBROUTINE GENRIGID_ROTATE(XCOORDS, ROTATEFACTOR)

USE COMMONS, ONLY: NATOMS
IMPLICIT NONE

INTEGER :: J1, J2, J3
DOUBLE PRECISION :: COM(3), NEWCOM(3), MASS, EXPANDVECTOR(3), PI, TWOPI, DPRAND
DOUBLE PRECISION :: ROTATIONMATRIX(3,3), TOROTATE(3), ATOMROTATED(3)
DOUBLE PRECISION :: RANDOMPHI, RANDOMTHETA, RANDOMPSI, ST, CT, SPH, CPH, SPS, CPS
DOUBLE PRECISION, INTENT(INOUT) :: XCOORDS(3*NATOMS)
DOUBLE PRECISION, INTENT(IN) :: ROTATEFACTOR

ROTATIONMATRIX(:,:) = 0.0D0
TOROTATE(:) = 0.0D0
! Define some constants
PI=ATAN(1.0D0)*4
TWOPI=2.0D0*PI

! Loop over all rigid bodies
DO J1 = 1, NRIGIDBODY
   IF (.NOT.FROZENRIGIDBODY(J1)) THEN
      COM = 0.0D0
      MASS = 0.0D0

! For each rigid body, calculate center of mass
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         COM = COM + XCOORDS(3*J3-2:3*J3)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
         MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
      ENDDO
      COM = COM / MASS

! Move the rigid body centre of mass to the origin
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) - COM
      ENDDO

! Calculate a random rotation scaled by ROTATEFACTOR
      RANDOMPHI=(DPRAND()-0.5)*TWOPI*ROTATEFACTOR
      RANDOMTHETA=(DPRAND()-0.5)*PI*ROTATEFACTOR
      RANDOMPSI=(DPRAND()-0.5)*TWOPI*ROTATEFACTOR
      ST=SIN(RANDOMTHETA)
      CT=COS(RANDOMTHETA)
      SPH=SIN(RANDOMPHI)
      CPH=COS(RANDOMPHI)
      SPS=SIN(RANDOMPSI)
      CPS=COS(RANDOMPSI)

! Assemble the rotation matrix
      ROTATIONMATRIX(1,1)=CPS*CPH-CT*SPH*SPS
      ROTATIONMATRIX(2,1)=CPS*SPH+CT*CPH*SPS
      ROTATIONMATRIX(3,1)=SPS*ST
      ROTATIONMATRIX(1,2)=-SPS*CPH-CT*SPH*CPS
      ROTATIONMATRIX(2,2)=-SPS*SPH+CT*CPH*CPS
      ROTATIONMATRIX(3,2)=CPS*ST
      ROTATIONMATRIX(1,3)=ST*SPH
      ROTATIONMATRIX(2,3)=-ST*CPH
      ROTATIONMATRIX(3,3)=CT

! Apply the rotation matrix to the atoms in the rigid body
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         TOROTATE = XCOORDS(3*J3-2:3*J3)
         ATOMROTATED=MATMUL(ROTATIONMATRIX,TOROTATE)
         XCOORDS(3*J3-2:3*J3) = ATOMROTATED
      ENDDO

! Translate the rigid body centre of mass back to its old position
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) + COM
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_ROTATE

! mo361> random rotation move for rigid bodies
SUBROUTINE GENRIGID_TRANSLATE(XCOORDS, TRANSLATEFACTOR)

USE COMMONS, ONLY: NATOMS
IMPLICIT NONE

INTEGER :: J1, J2, J3  
DOUBLE PRECISION DPRAND
DOUBLE PRECISION, INTENT(INOUT) :: XCOORDS(3*NATOMS) 
DOUBLE PRECISION, INTENT(IN) :: TRANSLATEFACTOR
DOUBLE PRECISION:: TRANSLATEVECTOR(3),LENGTH

! Loop over all rigid bodies
DO J1 = 1, NRIGIDBODY
   IF (.NOT.FROZENRIGIDBODY(J1)) THEN
      DO J2=1,3
         TRANSLATEVECTOR(J2)=2.0*(DPRAND()-0.5)*TRANSLATEFACTOR
      ENDDO
      LENGTH = DSQRT(TRANSLATEVECTOR(1)**2+TRANSLATEVECTOR(3)**2+TRANSLATEVECTOR(3)**2)
   
! Move the rigid body
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) + TRANSLATEVECTOR
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_TRANSLATE

! csw34> subroutine to update the reference coordinates for the rigid bodies using 
! the NATOMS coordinates in XCOORDS. Note that the rigid body coordinates are relative
! to the COM of each rigid body. 
SUBROUTINE GENRIGID_UPDATE_REFERENCE(XCOORDS)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2, DUMMY
  DOUBLE PRECISION, INTENT(IN) :: XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, MASS

  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0d0
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        XMASS = XMASS + XCOORDS(3*DUMMY-2)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XCOORDS(3*DUMMY-1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XCOORDS(3*DUMMY  )*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS  = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        SITESRIGIDBODY(J2,1,J1) = XCOORDS(3*DUMMY-2) - XMASS
        SITESRIGIDBODY(J2,2,J1) = XCOORDS(3*DUMMY-1) - YMASS
        SITESRIGIDBODY(J2,3,J1) = XCOORDS(3*DUMMY  ) - ZMASS
     ENDDO
  ENDDO

END SUBROUTINE GENRIGID_UPDATE_REFERENCE

! csw34> Apply a compression to the system. KCOMP_RIGID is the compression force constant
SUBROUTINE GENRIGID_COMPRESS(XRIGIDCOORDS,XRIGIDGRAD,EREAL,KCOMP_RIGID)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1
  DOUBLE PRECISION, INTENT(IN) :: XRIGIDCOORDS(3*NATOMS),KCOMP_RIGID
  DOUBLE PRECISION, INTENT(INOUT) :: XRIGIDGRAD(3*NATOMS),EREAL
  DOUBLE PRECISION :: SYSTEMCOM(3), COM(3), ATOM(3), DIST

! Calculate system centre of mass taking into account rigid bodies and single atoms
SYSTEMCOM = 0.0D0
! For the rigid bodies
DO J1 = 1,NRIGIDBODY
   SYSTEMCOM = SYSTEMCOM + XRIGIDCOORDS(3*J1-2:3*J1)
ENDDO
! For the single atoms
DO J1 = 6*NRIGIDBODY+1,DEGFREEDOMS
   SYSTEMCOM = SYSTEMCOM + XRIGIDCOORDS(3*J1-2:3*J1)
ENDDO 
SYSTEMCOM = SYSTEMCOM / DEGFREEDOMS

! Apply compression to rigid bodies
DO J1 = 1, NRIGIDBODY
   COM = XRIGIDCOORDS(3*J1-2:3*J1) 
! Work out squared distance of rigid body COM from system COM
   DIST = (COM(1)-SYSTEMCOM(1))**2+(COM(2)-SYSTEMCOM(2))**2+(COM(3)-SYSTEMCOM(3))**2
! Add energy
   EREAL = EREAL + KCOMP_RIGID*DIST/2.0D0
! Add gradient
   XRIGIDGRAD(3*J1-2)=XRIGIDGRAD(3*J1-2)+KCOMP_RIGID*(COM(1)-SYSTEMCOM(1))
   XRIGIDGRAD(3*J1-1)=XRIGIDGRAD(3*J1-1)+KCOMP_RIGID*(COM(2)-SYSTEMCOM(2))
   XRIGIDGRAD(3*J1  )=XRIGIDGRAD(3*J1  )+KCOMP_RIGID*(COM(3)-SYSTEMCOM(3))
ENDDO

! Apply compression to single atoms
DO J1 = 6*NRIGIDBODY+1,DEGFREEDOMS
   ATOM = XRIGIDCOORDS(3*J1-2:3*J1)
! Work out squared distance of atom from system COM
   DIST = (ATOM(1)-SYSTEMCOM(1))**2+(ATOM(2)-SYSTEMCOM(2))**2+(ATOM(3)-SYSTEMCOM(3))**2
! Add energy
   EREAL = EREAL + KCOMP_RIGID*DIST/2.0D0
! Add gradient
   XRIGIDGRAD(3*J1-2)=XRIGIDGRAD(3*J1-2)+KCOMP_RIGID*(ATOM(1)-SYSTEMCOM(1))
   XRIGIDGRAD(3*J1-1)=XRIGIDGRAD(3*J1-1)+KCOMP_RIGID*(ATOM(2)-SYSTEMCOM(2))
   XRIGIDGRAD(3*J1  )=XRIGIDGRAD(3*J1  )+KCOMP_RIGID*(ATOM(3)-SYSTEMCOM(3))
ENDDO

END SUBROUTINE GENRIGID_COMPRESS

! Subroutine to calculate the distance between the centre of mass of two rigid bodies
! XCOORDS must be in RIGID BODY FORMAT!
SUBROUTINE GENRIGID_COMDISTANCE(BODY1,BODY2,XCOORDS,COMDISTANCE)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2
  INTEGER, INTENT(IN) :: BODY1, BODY2
  DOUBLE PRECISION, INTENT(IN) :: XCOORDS(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: COMDISTANCE
  DOUBLE PRECISION :: BODY1COM(3), BODY2COM(3), MASS

! For each rigid body, calculate center of mass
! Body 1
MASS = 0.0D0
BODY1COM(:) = 0.0D0
BODY2COM(:) = 0.0D0

! You MUST be using rigid body coordinates here
! Assign BODY1/BODY2 COMs
BODY1COM = XCOORDS(3*BODY1-2:3*BODY1)
BODY2COM = XCOORDS(3*BODY2-2:3*BODY2)

! Calculate the COM-COM distance
COMDISTANCE = SQRT((BODY1COM(1)-BODY2COM(1))**2+(BODY1COM(2)-BODY2COM(2))**2+(BODY1COM(3)-BODY2COM(3))**2)


END SUBROUTINE GENRIGID_COMDISTANCE

! Subroutine to check that all rigid bodies are within DIST of each other. Returns TEST=.TRUE. if they are.
! We only require each rigid body to be within DIST of one other body to allow for large homogeneous systems.
! This might need changing in the future to require a user defined number of bodies to be within DIST.
SUBROUTINE GENRIGID_DISTANCECHECK(XCOORDS,DIST,TEST)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2
  DOUBLE PRECISION, INTENT(IN) :: XCOORDS(3*NATOMS), DIST
  DOUBLE PRECISION :: COMDIST, MINCOMDIST
  LOGICAL, INTENT(OUT) :: TEST

TEST = .TRUE.
! Loop over all pairs of rigid bodies
DO J1 = 1,NRIGIDBODY
! Set the initial minimum distance to a large value
   MINCOMDIST = HUGE(1.0D0)
   DO J2 = 1,NRIGIDBODY
! Skip self-self checks
      IF (J1.EQ.J2) CYCLE
! Check the COM-COM distance
      CALL GENRIGID_COMDISTANCE(J1,J2,XCOORDS,COMDIST)
! If the calculated distance is less than the current minimum, replace it
      IF (COMDIST.LT.MINCOMDIST) THEN
         MINCOMDIST = COMDIST
      ENDIF
   ENDDO
! Once body J1 has been compare to all bodies J2, check at least one is within DIST
   IF (MINCOMDIST.GT.DIST) THEN
      TEST = .FALSE.
      RETURN
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_DISTANCECHECK

END MODULE GENRIGID
