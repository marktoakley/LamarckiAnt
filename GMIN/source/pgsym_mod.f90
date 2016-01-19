! =================================================================
! !!!!!! The following module is exclusively for PGSYM !!!!!!!!!!!!
! =================================================================
!
MODULE PGSYMMOD
  !
  USE COMMONS, ONLY : MYUNIT, DEBUG
  !
  IMPLICIT NONE
  !
  INTEGER :: NAT ! number of atoms
  INTEGER :: IDEGEN ! counter for degenerate eigenvalues
  INTEGER :: NSYMOPS ! number of symmetry operations
  INTEGER, PARAMETER :: NSYMOPS_MAX=120
  INTEGER :: NROTAXES ! number of rotation axis
  INTEGER, PARAMETER :: NROTAXES_MAX=60
  INTEGER :: ROTAXESORDER(NROTAXES_MAX) ! Order of each rot. axis
  LOGICAL :: KEEPTRYING ! Flag for the TRY loop in PGSYM.
  LOGICAL :: LINEAR ! Flag linear system (one zero eigenvalue)
  CHARACTER(LEN=4) :: FPG ! Point group label 
  DOUBLE PRECISION :: TOL(3) ! Tolerances (can change!)
  DOUBLE PRECISION, ALLOCATABLE :: XREF(:,:) ! Reference coords.
  DOUBLE PRECISION, ALLOCATABLE :: XCUR(:,:) ! Current coords.
  DOUBLE PRECISION, ALLOCATABLE :: XCUR_NORM(:,:) ! Normalised XCUR
  DOUBLE PRECISION, ALLOCATABLE :: XTMP(:,:) ! Temporary coords
  INTEGER, ALLOCATABLE :: LABELS(:) ! Atomic labels
  INTEGER, ALLOCATABLE :: MAJORITY(:) ! List of atoms with majority label
  DOUBLE PRECISION :: NUIT(3,3) ! (Normalised & Uniform) Inertia
  DOUBLE PRECISION :: EVECS(3,3) ! Eigenvectors
  DOUBLE PRECISION :: SYMOPS(NSYMOPS_MAX,3,3) ! Symmetry matrices
  DOUBLE PRECISION :: ROTAXES(3,NROTAXES_MAX) ! Rotation axes
  !
  ! NOTE: ROTAXESORDER and ROTAXES are meant to be in-sync, with the
  ! last entry (at NROTAEXES) containing the axis of highest order.
  ! So complete sorting not required, but we admit a possible
  ! interchange of two tail-end entries when a new axis is added.
  !
CONTAINS
  !
  SUBROUTINE PGSYM_DIAGNOSE() !------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: I,J
    DOUBLE PRECISION :: TEMP1,TEMP2,VTEMP(3),MTEMP(3,3),RM(3,3)
    !
    ! Compute NUIT for XREF
    CALL CALC_NUIT(NAT,XREF,MAJORITY,NUIT)
    IF(DEBUG) THEN
       WRITE(MYUNIT,'(A)') &
            'pgsym_diagnose> Normalised (Uniform) Inertia Tensor:'
       DO I=1,3
          WRITE(MYUNIT,'(A,3(2X,F12.8))') &
               'pgsym_diagnose> ', ( NUIT(I,J), J=1,3 )
       ENDDO
    ENDIF
    !
    ! Diagonalise UIT and sort eigen-vals/vectors
    CALL EIG(NUIT,EVECS,3,3,0) ! in ptgrp.f (spaghetti code!)
    !
    ! Check for degeneracy of eigenvalues. If present, then check
    ! if the unique moment is along x. If so, change the unique 
    ! axis to z by rotating the eigenvector matrix about y by 90.
    IF( ABS(NUIT(2,2)-NUIT(3,3)) < TOL(1) ) THEN
       ! Build 3x3 rotation matrix RM about 2nd (i.e. y) axis
       VTEMP(:) = 0.0D0; VTEMP(2)=1.0D0
       CALL CALC_GEN_ROTMAT(VTEMP,90.0D0,0,RM)
       ! Compute MTEMP_ij = EVECS_ik*RM_kj, so that the coordinated
       ! transformation EVECS is preceeded by the rotation RM.
       CALL PRODMM(3,3,3,EVECS,RM,MTEMP) 
       EVECS(1:3,1:3) = MTEMP(1:3,1:3)
       TEMP1 = NUIT(1,1); NUIT(1,1) = NUIT(3,3); NUIT(3,3) = TEMP1
       IF(DEBUG) WRITE(MYUNIT,'(A)') &
            'pgsym_diagnose> Realigned the unique axis with Z.'
    ENDIF
    !
    IF(DEBUG) THEN
       WRITE(MYUNIT,'(A)') &
            'pgsym_diagnose> Diagonalized Inertia Tensor:'
       DO I=1,3
          WRITE(MYUNIT,'(A,3(2X,F12.8))') &
               'pgsym_diagnose> ', ( NUIT(I,J), J=1,3 )
       ENDDO
    END IF
    !                     
    ! Compute XCUR_ij = EVECS_ki * XREF_kj
    CALL PRODMTM(3,3,NAT,EVECS,XREF,XCUR)
    !
    IF(DEBUG) THEN
       ! To check for sanity, recompute the inertia tensor.
       CALL CALC_NUIT(NAT,XCUR,MAJORITY,NUIT)
       WRITE(MYUNIT,'(A)') &
            'pgsym_diagnose> Inertia Tensor for XCUR (sanity check):'
       DO I=1,3
          WRITE(MYUNIT,'(A,3(2X,F12.8))') &
               'pgsym_diagnose> ', ( NUIT(I,J), J=1,3 )
       ENDDO
    ENDIF
    !
    ! If (EV1 x EV2) . EV3 < 0 then we have left-handed inertial 
    ! axes, in which case we flip the sign of y axis to make it 
    ! right-handed. ds656> Not sure why this is necessary...
    CALL VCROSSPROD(EVECS(1:3,1), EVECS(1:3,2), VTEMP(1:3))
    TEMP1 = DOT_PRODUCT(VTEMP(1:3),EVECS(1:3,3))
    IF(TEMP1 < 0.0D0) XCUR(2,1:NAT) = -XCUR(2,1:NAT)
    !
    ! Store XCUR_NORM as it is useful for making the distance
    ! tolerance more robust.
    DO I=1,NAT
       TEMP1=DSQRT(DOT_PRODUCT(XCUR(1:3,I),XCUR(1:3,I)))
       XCUR_NORM(1:3,I) = XCUR(1:3,I)/TEMP1
    ENDDO
    !
    ! Determine the degeneracy index and check for linearity. 
    TEMP1 = -1.0D0
    IDEGEN=0
    LINEAR=.FALSE.
    DO I=1,3
       IF(NUIT(I,I) < TOL(1)) LINEAR=.TRUE.
       TEMP2 = NUIT(I,I)
       IF( ABS(TEMP2-TEMP1) < TOL(1) ) IDEGEN = IDEGEN+1
       TEMP1 = TEMP2
    ENDDO
    !
    ! Sanity check: if LINEAR==.TRUE., then should have IDEGEN=1
    !
    RETURN
    !
  END SUBROUTINE PGSYM_DIAGNOSE !----------------------------------
  !
  !---! Batch of high-level routines starts here !-----------------
  !
  SUBROUTINE PROC_ASYM_TOP() !-------------------------------------
    !
    ! Handles asymmetric top structures, which cannot have rotation
    ! axes of order higher than 2.
    !
    IMPLICIT NONE
    !
    CALL CHECK_R2_AXES_ASYM()
    !
    IF(NROTAXES==0) THEN ! No 2-fold axis.
       CALL PROC_NO_ROT_SYM()
    ELSEIF(NROTAXES==3) THEN ! Three 2-fold axes => dihedral group.
       CALL PROC_DIHEDRAL() 
    ELSE ! One or two 2-fold axes => a cyclic group.
       CALL PROC_CYCLIC() 
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_ASYM_TOP !-----------------------------------
  !
  SUBROUTINE PROC_LINEAR() !---------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL :: INVERSION
    !
    CALL CHECK_INVERSION(INVERSION)
    !
    IF(INVERSION) THEN
       FPG='D*h'
    ELSE
       FPG='C*v'
    ENDIF    
    !
    RETURN
    !
  END SUBROUTINE PROC_LINEAR !-------------------------------------
  !
  SUBROUTINE PROC_SYM_TOP() !--------------------------------------
    !
    ! Handles a symmetric top, with one unique eigenvalue whose
    ! corresponding principal axis is a unique rotational axis. 
    !
    ! More complex handling requires looking for R2 axes 
    ! perpendicular to this unique axis. 
    ! 
    ! NOTE: The unique rotational axis would have been 
    ! aligned with the z-axis in PGSYM_DIAGNOSE.
    !
    IMPLICIT NONE
    !
    CALL CHECK_ROT_SYM(3)
    IF(NROTAXES > 0) CALL CHECK_PERP_R2_AXIS(3)
    !
    IF(NROTAXES >= 2) THEN
       CALL PROC_DIHEDRAL()
    ELSEIF(NROTAXES == 1) THEN
       CALL PROC_CYCLIC()
    ELSE
       WRITE(MYUNIT, '(A)') 'proc_sym_top> Found no rotation axes.'
       IF(TOL(1) > 1.0D0-6) THEN
          WRITE(MYUNIT, '(A)') &
               'proc_sym_top> Tightening TOL(1).'
          TOL(1) = TOL(1)/10.0D0
          KEEPTRYING=.TRUE.
       ELSE
          WRITE(MYUNIT, '(A)') &
               'proc_sym_top> Accidental degeneracy.'
          CALL PROC_NO_ROT_SYM()
       ENDIF
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_SYM_TOP !-------------------------------------
  !
  SUBROUTINE PROC_SPH_TOP() !---------------------------------------
    !
    ! Handle spherical tops, which can belog to the following point
    ! groups: T, T_d, T_h, O, O_h, I, I_h.
    !
    IMPLICIT NONE
    !
    LOGICAL :: INVERSION
    CHARACTER(LEN=1) :: MIRROR_TYPE
    !
    CALL CHECK_SPH_AXES()
    !
    IF(NROTAXES > 0) THEN
       IF(ROTAXESORDER(NROTAXES) == 3) THEN
          CALL CHECK_MIRROR(ROTAXES(1:3, NROTAXES), MIRROR_TYPE)
          IF(MIRROR_TYPE /= "") THEN
             CALL CHECK_INVERSION(INVERSION)
             IF(INVERSION) THEN
                FPG='Th '
             ELSE
                FPG='Td '
             ENDIF
          ELSE
             FPG='T  '
          ENDIF
       ELSEIF(ROTAXESORDER(NROTAXES) == 4) THEN
          CALL CHECK_INVERSION(INVERSION)
          IF(INVERSION) THEN
             FPG='Oh '
          ELSE
             FPG='O  '
          ENDIF
       ELSEIF(ROTAXESORDER(NROTAXES) == 5) THEN
          CALL CHECK_INVERSION(INVERSION)
          IF(INVERSION) THEN
             FPG='Ih '
          ELSE
             FPG='I  '
          ENDIF
       ELSE
          WRITE(MYUNIT, '(A)') &
               'proc_sph_top> Found no 3, 4 or 5-fold axes.'
          IF(TOL(1) > 1.0D0-6) THEN
             WRITE(MYUNIT, '(A)') &
                  'proc_sph_top> Tightening TOL(1).'
             TOL(1) = TOL(1)/10.0D0
             KEEPTRYING=.TRUE.
          ELSE
             WRITE(MYUNIT, '(A)') &
                  'proc_sph_top> Accidental degeneracy.'
             CALL PROC_SYM_TOP()
          ENDIF
       ENDIF
    ELSE
       WRITE(MYUNIT, '(A)') &
            'proc_sph_top> Found no rotation axes.'
       IF(TOL(1) > 1.0D0-6) THEN
          WRITE(MYUNIT, '(A)') &
                  'proc_sph_top> Tightening TOL(1).'
          TOL(1) = TOL(1)/10.0D0
          KEEPTRYING=.TRUE.
       ELSE
          WRITE(MYUNIT, '(A)') &
               'proc_sph_top> Accidental degeneracy.'
          CALL PROC_ASYM_TOP()
       ENDIF
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_SPH_TOP !-------------------------------------
  !
  !---! Batch of medium-level routines starts here !----------------
  !
  SUBROUTINE PROC_NO_ROT_SYM() !------------------------------------
    !
    ! Handles moleculs with no rotational symmetry: C1, Cs and Ci
    !
    IMPLICIT NONE
    !
    LOGICAL INVERSION
    INTEGER :: I
    DOUBLE PRECISION :: AXIS(3)
    CHARACTER(LEN=1) :: MIRROR_TYPE
    !
    FPG='C1'
    CALL CHECK_INVERSION(INVERSION)
    IF(INVERSION) THEN
       FPG='Ci'
    ELSE
       DO I=1,3
          AXIS(1:3) = 0.0D0
          AXIS(I) = 1.0D0
          CALL CHECK_MIRROR(AXIS,MIRROR_TYPE)
          IF(TRIM(MIRROR_TYPE) /= "" ) THEN 
             FPG="Cs"
             RETURN
          ENDIF
       END DO
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_NO_ROT_SYM !----------------------------------
  !
  SUBROUTINE PROC_CYCLIC() !----------------------------------------
    !
    ! Handles cyclic groups, assuming the main rotation axis 
    ! (of some highest integer order n) has been identified. 
    ! Possibilities: C_n, C_nh, C_nv or S_2n
    !
    IMPLICIT NONE
    !
    LOGICAL :: OK
    CHARACTER(LEN=1) :: MIRROR_TYPE
    INTEGER :: IORDER
    DOUBLE PRECISION :: AXIS(3)
    !
    AXIS(1:3) = ROTAXES(1:3,NROTAXES)
    IORDER = ROTAXESORDER(NROTAXES)
    WRITE(FPG,'(A1,I1)') "C", IORDER
    !
    CALL CHECK_MIRROR(AXIS,MIRROR_TYPE)
    !
    IF(MIRROR_TYPE=="") THEN
       CALL CHECK_ROTOREFLECTION(AXIS, 180.0D0/DBLE(IORDER), OK)
       IF(OK) THEN
          WRITE(FPG,'(A1,I1)') "S", 2*IORDER
       ENDIF
    ELSEIF(MIRROR_TYPE=="d") THEN ! error?!"
       WRITE(MYUNIT,'(A)') 'pgsym> ERROR in PROC_CYCLIC !!!!'
       STOP
    ELSE ! MIRROR_TYPE is 'h' or 'v'
       FPG = TRIM(FPG)//MIRROR_TYPE
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_CYCLIC !--------------------------------------
  !
  SUBROUTINE PROC_DIHEDRAL() !--------------------------------------
    !
    ! Handles dihedral groups,i.e those with intersecting C2 axes
    ! and a main axis. Assume the main rotation axis 
    ! (of some highest integer order n) has been identified. 
    ! Possibilities: D_n, D_nh or D_nd.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=1) :: MIRROR_TYPE
    !
    WRITE(FPG,'(A1,I1)') "D",ROTAXESORDER(NROTAXES)
    !
    CALL CHECK_MIRROR(ROTAXES(1:3,NROTAXES),MIRROR_TYPE)
    !
    IF(MIRROR_TYPE == "h") THEN
       FPG=TRIM(FPG)//"h"
    ELSEIF(MIRROR_TYPE /= "") THEN
       FPG=TRIM(FPG)//"d"
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE PROC_DIHEDRAL !------------------------------------
  !
  !---! Batch of low-level checks and tests starts here !-----------
  !
  SUBROUTINE CHECK_INVERSION(INVERSION) !---------------------------
    !
    ! Test for inversion 
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: INVERSION
    INTEGER :: I, J
    DOUBLE PRECISION :: MAT(3,3), DUMMY
    !
    DO I=1,3
       MAT(I,I) = -1.0D0
       DO J=I+1,3
          MAT(I,J) = 0.0D0
          MAT(J,I) = 0.0D0
       ENDDO
    ENDDO
    !
    CALL TEST_SYMMAT(MAT,INVERSION,DUMMY)
    IF(INVERSION) THEN
       CALL ADD_SYMMAT(MAT)
       WRITE(MYUNIT, '(A)') &
            'check_inversion> Inversion symmetry present.'
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE CHECK_INVERSION !----------------------------------
  !
  SUBROUTINE CHECK_R2_AXES_ASYM() !--------------------------------
    !
    ! Test for 2-fold rotation around the principal axes.
    ! Used for asymmetric tops only, so there is no need
    ! to check if the detected generators already exist.
    !
    IMPLICIT NONE
    !
    LOGICAL :: OK
    INTEGER :: I
    DOUBLE PRECISION :: MAT(3,3), AXIS(3), DUMMY
    !
    DO I=1,3
       ! Get MAT that rotates 180 about principal axis I.
       AXIS(:)=0.0D0; AXIS(I)=1.0D0
       CALL CALC_GEN_ROTMAT(AXIS,180.0D0,0,MAT)
       ! Test for symmetry (and store if OK)
       CALL TEST_SYMMAT(MAT,OK,DUMMY)
       IF(OK) THEN
          CALL ADD_SYMMAT(MAT)
          CALL ADD_ROTAXIS(2,AXIS)
          WRITE(MYUNIT,'(A, I1)') &
               'check_r2_axes_asym> Found 2-fold axis.'
       ENDIF
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE CHECK_R2_AXES_ASYM !------------------------------
  !
  SUBROUTINE CHECK_MIRROR(AXIS, MIRROR_TYPE)
    !
    ! Looks for mirror symmetry of specified type about axis. 
    ! Possibilities are "h", "v" or "d". Horizontal (h) mirrors are 
    ! perpendicular to the supplied axis while vertical (v) or 
    ! diagonal (d) mirrors are parallel. A "v" mirror has atoms 
    ! lying on the mirror plane while "d" mirrors do not.
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: AXIS(3)
    CHARACTER(LEN=1), INTENT(OUT) :: MIRROR_TYPE
    !
    LOGICAL :: OK
    INTEGER :: I,J,K
    DOUBLE PRECISION :: NORMAL(3), V(3), MAT(3,3), DUMMY
    !
    MIRROR_TYPE=""
    !
    ! First check if the supplied axis is normal to a mirror plane
    ! 1) Get MATrix rep. for reflection about plane perp to AXIS.
    CALL CALC_REFLECMAT(AXIS,MAT)
    CALL TEST_SYMMAT(MAT,OK,DUMMY)
    !
    IF(OK) THEN
       CALL ADD_SYMMAT(MAT)
       MIRROR_TYPE="h"
    ELSE ! Messy business...
       ! Loop over ALL atom pairs to find candidate mirrors, whose
       ! normal is orthogonal to the supplied AXIS.
       LOUT: DO I=1,NAT-1
          DO J=I+1,NAT
             NORMAL(1:3) = XCUR(1:3,I) - XCUR(1:3,J)
             IF(ABS(DOT_PRODUCT(NORMAL,AXIS)) < TOL(2)) THEN
                ! 1) Get MAT for reflecting w.r.t. plane || to AXIS
                CALL CALC_REFLECMAT(NORMAL,MAT)
                CALL TEST_SYMMAT(MAT,OK,DUMMY)
                IF(OK) THEN
                   CALL ADD_SYMMAT(MAT)
                   IF(NROTAXES > 1) THEN
                      MIRROR_TYPE="d"
                      LIN: DO K=1,NROTAXES
                         V(1:3)=ROTAXES(1:3,K)-AXIS(1:3)
                         IF( DSQRT(DOT_PRODUCT(V, V)) >= TOL(2) ) THEN
                            IF( ABS(DOT_PRODUCT(ROTAXES(1:3,K), &
                                 NORMAL(1:3))) < TOL(2) ) THEN
                               MIRROR_TYPE="v"
                               EXIT LIN
                            ENDIF
                         ENDIF
                      ENDDO LIN
                   ELSE
                      MIRROR_TYPE="v"
                   ENDIF
                   EXIT LOUT
                ENDIF
             ENDIF
          ENDDO
       ENDDO LOUT
       !
    ENDIF
    !
    IF(MIRROR_TYPE /= "") THEN
       WRITE(MYUNIT, '(A, A, A)') &
            'check_mirror> Found ',MIRROR_TYPE,'-type mirror plane.'
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE CHECK_MIRROR !------------------------------------
  !
  SUBROUTINE CHECK_ROTOREFLECTION(AXIS, ANG, OK) !-----------------
    !
    ! Check for rotoreflection symmetry operation defined by the
    ! suplied AXIS and ANGle. 
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: AXIS(3), ANG
    LOGICAL, INTENT(OUT) :: OK
    !
    DOUBLE PRECISION :: MAT1(3,3), MAT2(3,3), MAT3(3,3), DUMMY
    !
    CALL CALC_GEN_ROTMAT(AXIS,ANG,0,MAT1)
    CALL CALC_REFLECMAT(AXIS,MAT2)
    CALL PRODMM(3,3,3,MAT2,MAT1,MAT3)
    CALL TEST_SYMMAT(MAT3,OK,DUMMY)
    !
    IF(OK) THEN
       CALL ADD_SYMMAT(MAT3)
       WRITE(MYUNIT, '(A)') &
            'check_rotoreflection> Found improper rotation axis.'
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE CHECK_ROTOREFLECTION !----------------------------
  !
  SUBROUTINE CHECK_ROT_SYM(IAXIS) !--------------------------------
    !
    ! Determines the rotational symmetry about a principal axis.  
    ! Used only for symmetric tops, which potentially have a
    ! rotational axis of order > 2.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: IAXIS
    !
    LOGICAL :: OK
    INTEGER :: LIST(0:NAT), I
    DOUBLE PRECISION :: MAT(3,3), AXIS(3), DUMMY
    !
    AXIS(:) = 0.0D0
    AXIS(IAXIS) = 1.0D0
    CALL GET_SMALLEST_SET_NOT_ON_AXIS(AXIS,LIST)
    !
    !write(myunit,*) 'check_rot_sym> LIST=', list(1:list(0))
    !
    DO I=LIST(0),2,-1
       ! Compute the rotation matrix MAT
       CALL CALC_GEN_ROTMAT(AXIS,360.0D0/DBLE(I),0,MAT)
       CALL TEST_SYMMAT(MAT,OK,DUMMY)
       IF(OK) THEN
          CALL ADD_SYMMAT(MAT)
          CALL ADD_ROTAXIS(I,AXIS)
          WRITE(MYUNIT, '(A, I1, A)') 'check_rot_sym> Found ',&
               I,'-fold axis.'
          !ds656> return now?!?
          RETURN
       ENDIF
       !
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE CHECK_ROT_SYM !-----------------------------------
  !
  SUBROUTINE CHECK_PERP_R2_AXIS(IAXIS) !---------------------------
    !
    ! Used for symmetric tops to check for R2 axes perpendicular to 
    ! the unique axis.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: IAXIS
    !
    LOGICAL :: YES
    INTEGER :: LIST(0:NAT), I,J
    DOUBLE PRECISION :: AXIS(3),TESTAXIS(3),DX(3),MAT(3,3),DUMMY
    !
    YES=.FALSE.
    AXIS(:) = 0.0D0
    AXIS(IAXIS) = 1.0D0
    CALL GET_SMALLEST_SET_NOT_ON_AXIS(AXIS,LIST)
    !
    DO I=1,LIST(0)-1
       DO J=I+1,LIST(0)
          DX(1:3) = XCUR(1:3,LIST(I)) - XCUR(1:3,LIST(J))
          ! ds656> Had I and J instead of LIST(I) and LIST(J) !!!!
          CALL NORMAL(DX,3)
          CALL VCROSSPROD(DX,AXIS,TESTAXIS)
          IF(DSQRT(DOT_PRODUCT(TESTAXIS,TESTAXIS)) > TOL(2)) THEN
             CALL CALC_GEN_ROTMAT(TESTAXIS,180.0D0,0,MAT)
             CALL TEST_SYMMAT(MAT,YES,DUMMY)
             IF(YES) THEN
                CALL ADD_SYMMAT(MAT)
                CALL ADD_ROTAXIS(2, TESTAXIS)
                WRITE(MYUNIT,'(A, I1)') &
                     'check_perp_r2_axis> Found 2-fold axis.'
                RETURN
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !ds656 > test
    !WRITE(MYUNIT,'(A, I1)') &
    !     'check_perp_r2_axis> Did not find 2-fold axis.'
    !<ds656 tested
    RETURN
    !
  END SUBROUTINE CHECK_PERP_R2_AXIS !-------------------------------
  !
  SUBROUTINE CHECK_SPH_AXES() !-------------------------------------
    !
    ! Look for R5, R4, R3 and R2 axes in speherical top molecules.  
    ! Point group T molecules have only one unique 3-fold and one 
    ! unique 2-fold axis. O molecules have one unique 4, 3 and 
    ! 2-fold axes. I molecules have a unique 5-fold axis.
    !
    USE COMMONS, ONLY : MYUNIT
    IMPLICIT NONE
    !
    LOGICAL :: ROT_PRESENT(5), OK
    INTEGER :: LISTS(2,0:NAT), LIST(0:NAT), I, IA, J, JA, K, L
    DOUBLE PRECISION :: X1(3), X2(3), X3(3), TAXIS(3), MAG, &
         DEG(1:5), MAT(3,3), DUMMY
    !
    CALL FIND_SHELLS(LISTS)
    ! Testing ...
    !WRITE(MYUNIT, '(A)') 'check_sph_axes> Shells:'
    !WRITE(MYUNIT, *) LISTS(1,1:NAT)
    !WRITE(MYUNIT, *) LISTS(2,1:LISTS(2,0))
    ! ... Tested.
    !
    ! Now find minimal shell
    LIST(0:NAT) = 0
    DO I=1,LISTS(2,0) ! Loop over shell starts 
       ! Find the next shell start and storas K
       IF(I==LISTS(2,0)) THEN  
          K=NAT+1
       ELSE
          K=LISTS(2,I+1)
       ENDIF
       J=K-LISTS(2,I) ! Atoms in shell
       IF(LIST(0) == 0 .OR. J < LIST(0)) THEN
          LIST(0) = J
          LIST(1:J) = LISTS(1,LISTS(2,I):K-1)
       ENDIF
       ! Testing ...
       ! DO J=LISTS(2,I), K-1
       !    WRITE(*,'(A,I1,3(1X,F20.10))') &
       !         'L',I, XCUR(1,LISTS(1,J)), XCUR(2,LISTS(1,J)), &
       !         XCUR(3,LISTS(1,J)) 
       ! ENDDO
       ! ... Tested
    ENDDO
    !
    !WRITE(MYUNIT, '(A)') 'check_sph_axes> Minimal shell:'
    !WRITE(MYUNIT, *) LIST(1:LIST(0))
    !
    ROT_PRESENT(1:5) = .FALSE.
    DO I=1,5
       DEG(I)=360.0D0/DBLE(I)
    ENDDO
    !
    ! Span combination 2- and 3-tuples in LIST.
    TRIP: DO I=1, LIST(0)-1 ! or LIST(0)-2?
       IA = LIST(I) ! Actual atom index
       DO J=I+1, LIST(0)  ! or LIST(0)-1?
          JA = LIST(J) ! Actual atom index
          ! Check for R2 axes from current pair. 
          ! Test axis bisecting the current pair "bond".
          IF( .NOT. ROT_PRESENT(2) ) THEN
             !WRITE(MYUNIT, *) 'check_sph_axes> R2:', IA, JA
             TAXIS(1:3) = XCUR_NORM(1:3,IA) + XCUR_NORM(1:3,JA)
             MAG=DSQRT(DOT_PRODUCT(TAXIS,TAXIS))
             !WRITE(MYUNIT, *) 'check_sph_axes> MAG and TOL:', &
             !     MAG, TOL(2)
             IF(MAG > TOL(2)) THEN
                CALL CALC_GEN_ROTMAT(TAXIS,DEG(2),0,MAT)
                !WRITE(MYUNIT, *) 'check_sph_axes> MAT:', MAT
                CALL TEST_SYMMAT(MAT,OK,DUMMY)
                IF(OK) THEN
                   CALL ADD_SYMMAT(MAT)
                   CALL ADD_ROTAXIS(2,TAXIS)
                   ROT_PRESENT(2) = .TRUE.
                   WRITE(MYUNIT, '(A)') &
                        'check_sph_axes> Found 2-fold axis.'
                ENDIF
             ENDIF
          ENDIF
          !
          X1(1:3) = XCUR_NORM(1:3,JA) - XCUR_NORM(1:3,IA)
          !
          DO K=J+1, LIST(0)
             !
             ! Check for R3, R4, and R5 axes from current triplet.
             ! Test the normal to the plane defined by the
             ! triplet, if the triplet is not (approximately) linear.
             !
             X2(1:3) = XCUR_NORM(1:3,LIST(K)) - XCUR_NORM(1:3,IA)
             !
             !WRITE(MYUNIT, *) 'check_sph_axes> R3:', IA, JA, LIST(K)
             CALL VCROSSPROD(X1,X2,TAXIS)
             MAG=DSQRT(DOT_PRODUCT(TAXIS,TAXIS))
             !WRITE(MYUNIT, *) 'check_sph_axes> MAG and TOL:', &
             !     MAG, TOL(2)
             IF( MAG > TOL(2) ) THEN
                FOLD: DO L=2,5
                   IF( .NOT. ROT_PRESENT(L) ) THEN
                      CALL CALC_GEN_ROTMAT(TAXIS,DEG(L),0,MAT)
                      CALL TEST_SYMMAT(MAT,OK,DUMMY)
                      IF(OK) THEN
                         CALL ADD_SYMMAT(MAT)
                         CALL ADD_ROTAXIS(L,TAXIS)
                         ROT_PRESENT(L) = .TRUE.
                         WRITE(MYUNIT, '(A, I1, A)') &
                              'check_sph_axes> Found ',L,'-fold axis.'
                         EXIT FOLD
                      ENDIF
                   ENDIF
                ENDDO FOLD
             ENDIF
             !
             IF( ROT_PRESENT(2) .AND. ROT_PRESENT(3) .AND. &
                  (ROT_PRESENT(4) .OR. ROT_PRESENT(5)) ) EXIT TRIP
             !
          ENDDO
       ENDDO
    ENDDO TRIP
    !
    RETURN
    !
  END SUBROUTINE CHECK_SPH_AXES !----------------------------------
  !
  SUBROUTINE CALC_REFLECMAT(NORMAL, MAT) !-------------------------
    !
    ! Returns a reflection symmetry operation (MAT) for the 
    ! specified NORMAL vector. The origin is presumed to be on
    ! the mirror plane (so no translation is applie).
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: NORMAL(3)
    DOUBLE PRECISION, INTENT(OUT) :: MAT(3,3)
    !
    INTEGER :: I,J
    DOUBLE PRECISION :: NNORM(3)
    !
    ! Normalise the normal vector.
    NNORM(:)=NORMAL(:)/DSQRT(DOT_PRODUCT(NORMAL,NORMAL))
    ! Compute the 
    DO I=1,3
       MAT(I,I) = 1.0D0 - 2.0D0*NNORM(I)**2
       DO J=I+1,3
          MAT(I,J) = -2.0D0*NNORM(I)*NNORM(J)
          MAT(J,I) = MAT(I,J)
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE CALC_REFLECMAT !----------------------------------
  !
  SUBROUTINE CALC_GEN_ROTMAT(AXIS, ANGLE, ISRAD, MAT) !------------
    !
    ! Compute a rotation matrix about a supplied AXIS and ANGLE.
    ! Use standard formula to compute the matrix components.
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: AXIS(3), ANGLE
    INTEGER, INTENT(IN) :: ISRAD
    DOUBLE PRECISION, INTENT(OUT) :: MAT(3,3)
    !
    DOUBLE PRECISION :: LANG, COSANG, SINANG, NAXIS(3), MCOSANG
    !
    IF(ISRAD==1) THEN
       LANG = ANGLE
    ELSE
       LANG = ANGLE*ACOS(-1.0D0)/180.0D0
    ENDIF
    !
    COSANG = COS(LANG)
    MCOSANG = 1.0D0-COSANG
    SINANG = SIN(LANG)
    !
    NAXIS(:) = AXIS(:)/DSQRT(DOT_PRODUCT(AXIS,AXIS))
    !
    ! Can this be written more compactly?
    MAT(1,1) = COSANG + NAXIS(1)**2*MCOSANG
    MAT(1,2) = NAXIS(1)*NAXIS(2)*MCOSANG - NAXIS(3)*SINANG
    MAT(1,3) = NAXIS(1)*NAXIS(3)*MCOSANG + NAXIS(2)*SINANG
    MAT(2,1) = NAXIS(1)*NAXIS(2)*MCOSANG + NAXIS(3)*SINANG
    MAT(2,2) = COSANG + NAXIS(2)**2*MCOSANG
    MAT(2,3) = NAXIS(2)*NAXIS(3)*MCOSANG - NAXIS(1)*SINANG
    MAT(3,1) = NAXIS(1)*NAXIS(3)*MCOSANG - NAXIS(2)*SINANG
    MAT(3,2) = NAXIS(2)*NAXIS(3)*MCOSANG + NAXIS(1)*SINANG
    MAT(3,3) = COSANG + NAXIS(3)**2*MCOSANG
    !
    RETURN
    !
  END SUBROUTINE CALC_GEN_ROTMAT !---------------------------------
  !
  SUBROUTINE GET_SMALLEST_SET_NOT_ON_AXIS(AXIS, LIST) !------------
    !
    ! Returns the smallest list of atoms with the same distance 
    ! from origin (i.e. same shell), with no atoms lying on the 
    ! specified axis. This maximal set limits the possible 
    ! rotational symmetry operations, since atoms lying on a test 
    ! axis are irrelevant in testing rotational symmetry.
    !
    IMPLICIT NONE 
    !
    DOUBLE PRECISION, INTENT(IN) :: AXIS(3)
    INTEGER, INTENT(OUT) :: LIST(0:NAT)
    !
    INTEGER :: I,J,JMAX,K,LISTS(2,0:NAT),LIST_TMP(0:NAT)
    DOUBLE PRECISION :: V1(3),V2(3),V3(3),AXIS_NORM(3),DIST2,TOL2
    !
    CALL FIND_SHELLS(LISTS)
    !
    ! Initialise list
    LIST(0:NAT) = 0
    ! Noramlise supplied axis to avoid division inside loop.
    AXIS_NORM(1:3) = AXIS(1:3)/DSQRT(DOT_PRODUCT(AXIS,AXIS))
    ! Square the relevant tolerance for the same reason.
    TOL2 = TOL(2)**2
    !
    DO I=1,LISTS(2,0) ! Loop over shells
       IF(I == LISTS(2,0)) THEN
          JMAX=NAT
       ELSE
          JMAX=LISTS(2,I+1)-1
       ENDIF
       LIST_TMP(0:NAT) = 0
       DO J=LISTS(2,I),JMAX ! Loop over atoms in shell
          K=LISTS(1,J) ! Actual atom index
          V1(1:3) = XCUR(1:3,K)
          V2(1:3) = V1(1:3) - AXIS_NORM(1:3)
          ! Now compute shortest distance from atom to axis.
          CALL VCROSSPROD(V1,V2,V3)
          DIST2 = DOT_PRODUCT(V3,V3)
          IF (DIST2 > TOL2) THEN ! append to current list.
             LIST_TMP(0) = LIST_TMP(0) + 1
             LIST_TMP(LIST_TMP(0)) = K
          ENDIF
       ENDDO
       IF(LIST(0) == 0 .OR. &
            (LIST_TMP(0) < LIST(0) .AND. LIST_TMP(0) > 0) ) THEN
          LIST(0:NAT) = LIST_TMP(0:NAT)
       ENDIF
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE GET_SMALLEST_SET_NOT_ON_AXIS !--------------------
  !
  SUBROUTINE FIND_SHELLS(LISTS) !----------------------------------
    !
    ! Partition the list of atom indices into orbits, with each
    ! orbit containing atoms that are (approximately) equidistant 
    ! from the origin. Procedure:
    ! 1) Compute the norm for each atom position vector.
    ! 2) Sort atom indices by the norm.
    ! 3) Look for "gaps" > tolearance between adjecant norms to
    !    define a partition.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT) :: LISTS(1:2,0:NAT)
    !
    INTEGER :: I, J, K
    DOUBLE PRECISION :: R(3), DIST, DISTS(0:NAT)
    !
    LISTS(:,:) = 0
    DISTS(:) = 0.0D0
    !
    ! First loop over atoms to sort them by radial distance.
    DO I=1,NAT
       R(1:3) = XCUR(1:3,I)
       ! Append atom index and sistance to the lists.
       DISTS(I) = DSQRT(DOT_PRODUCT(R,R))
       LISTS(1,I) = I
       ! Now shuffle down the list, relying on DISTS2(0)=0.0d0.
       J=I
       DO WHILE(DISTS(J) < DISTS(J-1)) 
          K=LISTS(1,J); LISTS(1,J)=LISTS(1,J-1); LISTS(1,J-1)=K
          DIST=DISTS(J); DISTS(J)=DISTS(J-1); DISTS(J-1)=DIST
          J=J-1
       ENDDO
    ENDDO
    !
    ! Loop through sorted lists to identify gaps separating shells.
    ! If the 1st atom is at the origin, it will not be included in
    ! any shell.
    DO I=1,NAT
       IF(DISTS(I) > DISTS(I-1) + TOL(2)) THEN
          ! I is immediately preceded by a gap -> new shell!
          LISTS(2,0) = LISTS(2,0) + 1
          LISTS(2,LISTS(2,0)) = I 
       ENDIF
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE FIND_SHELLS !-------------------------------------
  !
  SUBROUTINE TEST_SYMMAT(MAT,OK,DIST) !----------------------------
    !
    ! Test if the supplied MAT indeed corresponds to a symmetry 
    ! operation, in which case each atom in X1 = XCUR has a partner 
    ! in X2=MAT*X1 within a global distance tolerance TOL, 
    ! i.e.: X1(i) <~> X2(perm(i))
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT (IN) :: MAT(3,3)
    LOGICAL, INTENT(OUT) :: OK
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    !
    LOGICAL :: ASSIGNED(NAT)
    INTEGER :: J1,J2,J3
    DOUBLE PRECISION :: DIST2,TOL2
    !
    ! Get XTMP_ij = MAT_ik * XCUR_kj, where i,k=1,2,3; j=1,..,NAT
    CALL PRODMM(3,3,NAT,MAT,XCUR,XTMP)
    !
    TOL2 = TOL(2)*TOL(2)
    ASSIGNED(1:NAT) = .FALSE.
    !
    L1: DO J1=1,NAT
       DIST=10.0D10
       OK=.FALSE.
       L2: DO J2=1,NAT
          IF(ASSIGNED(J2)) CYCLE L2
          IF(LABELS(J1) /= LABELS(J2)) CYCLE L2
          DIST2=0.0D0
          DO J3=1,3
             DIST2 = DIST2 + (XCUR(J3,J1)-XTMP(J3,J2))**2
          ENDDO
          IF(DIST > DIST2) DIST = DIST2
          !WRITE(MYUNIT, *) 'test_symmat> DIST2, TOL2:', &
          !     DIST2, TOL2
          IF(DIST2 < TOL2) THEN
             OK=.TRUE.
             ASSIGNED(J2)=.TRUE.
             EXIT L2
          ENDIF
       ENDDO L2
       IF (.NOT. OK) EXIT L1
    END DO L1
    !
    DIST = DSQRT(DIST)
    !
    RETURN
    !
  END SUBROUTINE TEST_SYMMAT !-------------------------------------
  !
  SUBROUTINE ADD_SYMMAT(MAT) !-------------------------------------
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: MAT(3,3)
    INTEGER :: J1,J2
    !
    IF(NSYMOPS < NSYMOPS_MAX) THEN
       NSYMOPS = NSYMOPS + 1
       SYMOPS(NSYMOPS,1:3,1:3) = MAT(1:3,1:3)
    ELSE
       WRITE(MYUNIT,'(A)') &
            'add_symmat> ERROR: NSYMOPS > NSYMOPS_MAX.'
       STOP
    ENDIF
    !
    IF(DEBUG) THEN
      WRITE(MYUNIT, '(A)') 'add_symmat> Symm. op. added!'
       DO J1=1,3
          WRITE(MYUNIT,'(3(1X,F20.10))') &
               ( MAT(J1,J2), J2=1,3 )        
       ENDDO
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE ADD_SYMMAT !--------------------------------------
  !
  SUBROUTINE ADD_ROTAXIS(ORDER, AXIS) !----------------------------
    !
    ! Add AXIS and its ORDER to global lists ROTAXIS
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ORDER
    DOUBLE PRECISION, INTENT(IN) :: AXIS(3)
    !
    INTEGER :: NOLD
    DOUBLE PRECISION :: TMP(3)
    !
    IF(NROTAXES < NROTAXES_MAX) THEN
       NOLD = NROTAXES
       NROTAXES = NROTAXES + 1
       IF( NOLD == 0 ) THEN
          ROTAXES(1:3,NROTAXES) = AXIS(1:3)
          ROTAXESORDER(NROTAXES) = ORDER          
       ELSEIF( ORDER > ROTAXESORDER(NOLD) ) THEN
          ROTAXES(1:3,NROTAXES) = AXIS(1:3)
          ROTAXESORDER(NROTAXES) = ORDER
       ELSE
          ROTAXES(1:3,NROTAXES) = ROTAXES(1:3,NOLD)
          ROTAXESORDER(NROTAXES) = ROTAXESORDER(NOLD)
          ROTAXES(1:3,NOLD) = AXIS(1:3)
          ROTAXESORDER(NOLD) = ORDER          
       ENDIF
    ELSE
       WRITE(MYUNIT, '(A)') &
            'pgsym> NROTAXES > NROTAXES_MAX, stopping!'
       STOP
    ENDIF
    !
  END SUBROUTINE ADD_ROTAXIS !-------------------------------------
  !
  SUBROUTINE TRANSFORM_SYMOPS() !----------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL :: OK, WARN
    INTEGER :: I,J,K
    DOUBLE PRECISION :: RMAT1(3,3), RMAT2(3,3), RMAT3(3,3), DUMMY
    !
    ! Find the matrix RMAT2 that maps XREF onto XCUR.
    CALL FIND_RMAT(NAT,XCUR,XREF,RMAT2)
    !
    WARN=.FALSE.
    !
    ! Reset XCUR to XREF
    XCUR(:,:) = XREF(:,:)
    !
    DO I=1,NSYMOPS
       !
       RMAT1(1:3,1:3) = SYMOPS(I,1:3,1:3)
       !
       ! Make RMAT1 consistent with the reference state:
       ! 1) Compute product RMAT1 * RMAT2 -> RMAT3
       CALL PRODMM(3,3,3,RMAT1,RMAT2,RMAT3)
       ! 2) Compute product RMAT2^t * RMAT3 -> RMAT1
       CALL PRODMTM(3,3,3,RMAT2,RMAT3,RMAT1)
       !
       ! Test the sym. ops. in the reference frame of XREF.
       CALL TEST_SYMMAT(RMAT1, OK, DUMMY)
       IF(.NOT. OK) THEN
          WARN=.TRUE.
          WRITE(MYUNIT, '(A,I3,A,F8.5)') &
               'transform_symops> Sym. op. ',I,&
               ' failed with mismatch ', DUMMY
       ENDIF
       !
       SYMOPS(I,1:3,1:3) = RMAT1(1:3,1:3)
       !
    END DO
    !
    IF(WARN) THEN
       WRITE(MYUNIT,'(A, A)') &
            'transform_symops> Reorientation failed,', &
            ' discarding all sym. ops.'
       NSYMOPS = 0
       SYMOPS(:,:,:) = 0.0D0
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE TRANSFORM_SYMOPS !--------------------------------
  !
  SUBROUTINE COMPLETE_SYMOPS() !-----------------------------------
    !
    ! permute through all possible combinations of the initially 
    ! supplied symmetry operations to arrive at a complete set of
    ! operations mapping a single atom to all other equivalent 
    ! atoms in the point group. This assumes that the initial 
    ! number already uniquely identifies all operations.
    !
    IMPLICIT NONE
    !
    LOGICAL :: SAME
    INTEGER :: I, J, K
    DOUBLE PRECISION :: MAT1(3,3), MAT2(3,3), MAT3(3,3), DUMMY
    !
    I=0
    DO WHILE (I < NSYMOPS)
       I=I+1
       MAT1(:,:) = SYMOPS(I,:,:)
       J=0
       DO WHILE(J < NSYMOPS)
          J=J+1
          MAT2(:,:) = SYMOPS(J,:,:)
          CALL PRODMM(3,3,3,MAT1,MAT2,MAT3)
          TEST: DO K=1,NSYMOPS
             MAT2(:,:) = SYMOPS(K,:,:)
             CALL COMPARE_SYMMATS(MAT3,MAT2,TOL(3),SAME)
             IF(SAME) EXIT TEST
          ENDDO TEST
          IF(.NOT. SAME) THEN
             CALL TEST_SYMMAT(MAT3,SAME,DUMMY)
             IF(SAME) CALL ADD_SYMMAT(MAT3)
          ENDIF
       ENDDO
    ENDDO
    !
    IF(NSYMOPS == 0) THEN
       MAT3(1:3,1:3) = 0
       DO I=1,3
          MAT3(I,I) = 1.0D0
       ENDDO
       CALL ADD_SYMMAT(MAT3)
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE COMPLETE_SYMOPS !---------------------------------
  !
  !---! Batch of miscellaneous routine !----------------------------
  !
  SUBROUTINE ALLOCATE_ALLX_AND_INI() !-----------------------------
    !
    IMPLICIT NONE
    !
    ALLOCATE(LABELS(NAT),MAJORITY(0:NAT),XREF(3,NAT),XCUR(3,NAT),&
         XCUR_NORM(3,NAT),XTMP(3,NAT)) 
    !
    NSYMOPS = 0
    SYMOPS(1:NSYMOPS_MAX,1:3,1:3) = 0.0D0
    NROTAXES = 0
    ROTAXESORDER(1:NROTAXES_MAX) = 0 
    ROTAXES(1:3,1:NROTAXES_MAX) = 0.0D0
    LINEAR = .FALSE. 
    FPG='   '
    XREF(1:3,1:NAT) = 0.0D0
    XCUR(1:3,1:NAT) = 0.0D0
    XTMP(1:3,1:NAT) = 0.0D0
    NUIT(1:3,1:3) = 0.0D0
    EVECS(1:3,1:3) = 0.0D0
    !
    RETURN
    !
  END SUBROUTINE ALLOCATE_ALLX_AND_INI !---------------------------
  !
  SUBROUTINE DEALLOCATE_ALLX() !-----------------------------------
    !
    IMPLICIT NONE
    !
    DEALLOCATE(XREF,XCUR,XCUR_NORM,XTMP,LABELS,MAJORITY)
    !
    RETURN
    !
  END SUBROUTINE DEALLOCATE_ALLX !---------------------------------
  !
END MODULE PGSYMMOD
