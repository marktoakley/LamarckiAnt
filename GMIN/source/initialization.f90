!   GMIN: A program for finding global minima
!   Copyright (C) 1999- David J. Wales
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

!   js850> This subroutine will be used to do all initializations, e.g. prepare
!   variables/arrays used in potentials, take step routines, etc.  It is assumed
!   that at this point the keywords have been read and the coords have been
!   loaded.

SUBROUTINE INITIALIZATIONS()
USE COMMONS
USE QMODULE
USE CLASS_OVERLAP
USE BGUPMOD ! In BGupta.f90
USE POT_PARAMS 
USE POLIRMOD, ONLY: POLIRINIT
USE HOMOREFMOD ! in homoref_addons.f90
IMPLICIT NONE
INTEGER LUNIT, GETUNIT
integer j1, j2, j6
DOUBLE PRECISION TEMPCOORDS(3*NATOMS)
DOUBLE PRECISION GRAD(NATOMS*3), ENERGY, DUMMY
LOGICAL GTEST

IF (ONEDAPBCT.OR.ONEDPBCT) THEN
   ALLOCATE(XYPHI(3*NATOMS))
   LUNIT=GETUNIT()
   OPEN (LUNIT,FILE='ONED.phi',STATUS='OLD')
   READ(LUNIT,*) XYPHI(1:3*NATOMS)
   CLOSE(LUNIT)
   WRITE(MYUNIT,'(A)') 'initialization> phi values'
   WRITE(MYUNIT,'(3G20.10)') XYPHI(1:3*NATOMS)
ENDIF

IF (TWODAPBCT.OR.TWODPBCT) THEN
   ALLOCATE(XYPHI(2*(NONEDAPBC)*(NONEDAPBC)))
   LUNIT=GETUNIT()
   OPEN (LUNIT,FILE='TWOD.phi',STATUS='OLD')
   READ(LUNIT,*) XYPHI(1:(2*NONEDAPBC*NONEDAPBC))
   CLOSE(LUNIT)
   WRITE(MYUNIT,'(A,I5,I5)') 'initialization> phi values, NONED, NATOMS', NONEDAPBC, NATOMS
   WRITE(MYUNIT,'(3G20.10)') XYPHI(1:(2*NONEDAPBC*NONEDAPBC))
ENDIF

IF (THREEDAPBCT.OR.THREEDPBCT) THEN
   ALLOCATE(XYPHI(3*NONEDAPBC**3))
   LUNIT=GETUNIT()
   OPEN (LUNIT,FILE='THREED.phi',STATUS='OLD')
   READ(LUNIT,*) XYPHI(1:3*NONEDAPBC**3)
   CLOSE(LUNIT)
   WRITE(MYUNIT,'(A)') 'initialization> phi values'
   WRITE(MYUNIT,'(3G20.10)') XYPHI(1:3*NONEDAPBC**3)
ENDIF




! js850> if PERIODIC, then make the default for MAXERISE much larger.  This is
! to avoid failed quenches because of discontinuities due to the existence of
! multiple images
IF ( PERIODIC .AND. .NOT. MAXERISE_SET ) MAXERISE = 1.D-3

! js850> set HARMONICR0 to be the initial coords
IF ( HARMONICF ) THEN
   HARMONICR0(1:3*NATOMS) = COORDS(1:3*NATOMS,1)
ENDIF

! js850> create FROZENLIST from FROZEN
! FROZENLIST holds the sorted list of frozen atoms between 1 and NFREEZE.
! Between NFREEZE+1 and N it holds the sorted list of unfrozen particles.
IF ( FREEZE ) THEN
   ALLOCATE(FROZENLIST(NATOMS))
   J1=0
   J2=NFREEZE
   DO J6=1,NATOMS
      IF ( (BINARY .OR. SOFT_SPHERE) .AND. J6 .EQ. NTYPEA ) THEN
         IF ( FROZEN(J6) ) THEN
            NFREEZETYPEA = J1+1
         ELSE
            NFREEZETYPEA = J1
         ENDIF
      ENDIF
      IF ( FROZEN(J6) ) THEN
         J1 = J1+1
         FROZENLIST(J1) = J6
      ELSE
         J2 = J2+1
         FROZENLIST(J2) = J6
      ENDIF
   END DO
ENDIF

! js850> initialize overlap
IF (OVERLAPK) THEN
   IF (OVERLAP_IMPORT) THEN
      WRITE(MYUNIT,*) " reading overlap comparison file from overlap.input"
      LUNIT = GETUNIT()
      OPEN (LUNIT,FILE='overlap.input',STATUS='OLD')
      DO J1=1,NATOMS
         J2=3*(J1-1)
         READ(LUNIT,*) TEMPCOORDS(J2+1), TEMPCOORDS(J2+2), TEMPCOORDS(J2+3)
      ENDDO
      CLOSE (LUNIT)
   ELSE
      TEMPCOORDS(:) = COORDS(:,1)
   ENDIF
   CALL OVERLAP_SETUP( TEMPCOORDS, OVERLAP_DR)
   !sanity check
   !CALL OVERLAP_get_overlap( COORDS(1:3*NATOMS,1), dummy)
   !write(*,*) "initial overlap", dummy
ENDIF

! csw34> initialize hydrogen-bond grouping
! Maximum of 1000 groups by default in keyword.f
IF (HBONDMATRIX) THEN
   ALLOCATE(HBONDGROUPS(MAXHBONDGROUPS,HBONDNRES,HBONDNRES))
   ALLOCATE(HBONDMAT(HBONDNRES,HBONDNRES))
   ALLOCATE(HBONDSOFTMAT(HBONDNRES,HBONDNRES))
   ALLOCATE(HBONDGROUPPOP(MAXHBONDGROUPS))
   ALLOCATE(HBONDBEST(MAXHBONDGROUPS))
   ALLOCATE(HBONDMARKOV(MAXHBONDGROUPS))
   ALLOCATE(HBONDMAXE(MAXHBONDGROUPS))
   ALLOCATE(HBONDQUE(MCSTEPS(1)))
   ALLOCATE(HBONDGROUPIDS(MCSTEPS(1)))
   ALLOCATE(HBONDBESTCOORDS(MAXHBONDGROUPS,3*NATOMS))
   HBONDGROUPIDS(:)=0
   HBONDGROUPPOP(:)=0
   HBONDGROUPS(:,:,:)=0
   HBONDMAT(:,:)=0
   HBONDSOFTMAT(:,:)=0
   HBONDMARKOV(:)=HUGE(1.0D0)
   HBONDBEST(:)=HUGE(1.0D0)
   HBONDQUE(:)=HUGE(1.0D0)
   HBONDBESTCOORDS(:,:)=0.0D0
   HBONDMAXE(:)=-HUGE(1.0D0)
   HBONDESAVE=HUGE(1.0D0)
   HBONDQUZEROE=HUGE(1.0D0)
ENDIF


!ds656> If NSPECIES is not allocated, presume homogeneneous system
IF(.NOT. ALLOCATED(NSPECIES)) THEN 
   ALLOCATE(NSPECIES(0:1))
   NSPECIES(0)=1          ! one species
   NSPECIES(1)=NATOMS     ! all atoms are of the same species
ELSEIF(NSPECIES(0)==1) THEN
   NSPECIES(1)=NATOMS     ! all atoms are of the same species
ENDIF
NSPECIES_INI(:) = NSPECIES(:)
!
!ds656> Compute the local neighboruhood size in permutation space
QALCS_NBRHD = 0
DO J1=1,NSPECIES(0)-1
   DO J2=J1+1,NSPECIES(0)
      QALCS_NBRHD = QALCS_NBRHD + NSPECIES(J1)*NSPECIES(J2)
   ENDDO
ENDDO
!
!ds656> Allocate the stress tensor if need be
IF(STRESST) THEN
   ALLOCATE(STRESS(0:NATOMS,1:3,1:3))
ENDIF
!
!ds656> Initialise the 5 Binary Gupta parameters for each atom type.
IF (BGUPTAT) THEN
   !AA
   AARRAY2(1,1)=2.0d0*AAA
   AARRAY(1,1)=AAA
   AARRAY_MOD(1,1)=AAB
   PARRAY(1,1)=PAA
   PARRAY_MOD(1,1)=PAB
   QARRAY(1,1)=QAA
   QARRAY_MOD(1,1)=QAB
   ZARRAY(1,1)=ZAA
   ZARRAY_MOD(1,1)=ZAB
   R0ARRAY(1,1)=R0AA
   R0ARRAY_MOD(1,1)=R0AB
   !BB
   AARRAY2(2,2)=2.0d0*ABB
   AARRAY(2,2)=ABB
   AARRAY_MOD(2,2)=AAB
   PARRAY(2,2)=PBB
   PARRAY_MOD(2,2)=PAB
   QARRAY(2,2)=QBB
   QARRAY_MOD(2,2)=QAB
   ZARRAY(2,2)=ZBB
   ZARRAY_MOD(2,2)=ZAB
   R0ARRAY(2,2)=R0BB
   R0ARRAY_MOD(2,2)=R0AB
   !BA
   AARRAY2(2,1)=2.0d0*AAB
   AARRAY(2,1)=AAB
   AARRAY_MOD(2,1)=AAA
   PARRAY(2,1)=PAB
   PARRAY_MOD(2,1)=PAA
   QARRAY(2,1)=QAB
   QARRAY_MOD(2,1)=QAA
   ZARRAY(2,1)=ZAB
   ZARRAY_MOD(2,1)=ZAA
   R0ARRAY(2,1)=R0AB
   R0ARRAY_MOD(2,1)=R0AA
   !AB
   AARRAY2(1,2)=2.0d0*AAB
   AARRAY(1,2)=AAB
   AARRAY_MOD(1,2)=ABB
   PARRAY(1,2)=PAB
   PARRAY_MOD(1,2)=PBB
   QARRAY(1,2)=QAB
   QARRAY_MOD(1,2)=QBB
   ZARRAY(1,2)=ZAB
   ZARRAY_MOD(1,2)=ZBB
   R0ARRAY(1,2)=R0AB
   R0ARRAY_MOD(1,2)=R0BB
   !
ENDIF !<ds656
!
! ds656> allocate some arrays for homotop enumeration:
IF(ENPERMST) THEN
   ALLOCATE(LEHMER_LIST(NATOMS,NPAR), LEHMER_COORDS(3*NATOMS,NPAR), &
        LEHMER_LAST(NPAR),LEHMER_ILASTB(NPAR))
   LEHMER_COORDS(1:3*NATOMS,1:NPAR) = COORDS(1:3*NATOMS,1:NPAR)
   LEHMER_LAST(1:NPAR) = 'B'
   LEHMER_ILASTB(1:NPAR) = NATOMS
   LEHMER_LIST(1:NTYPEA,1:NPAR) = 'A'
   LEHMER_LIST(NTYPEA+1:NATOMS,1:NPAR) = 'B'
ENDIF
!
IF(MSCT .AND. CUTT) THEN 
   !
   ! Allocate and initialise some type-dependent constants
   ! for smooth truncation.
   IF(ALLOCATED(CUTA_REP)) DEALLOCATE(CUTA_REP)
   IF(ALLOCATED(CUTB_REP)) DEALLOCATE(CUTB_REP)
   IF(ALLOCATED(CUTC_REP)) DEALLOCATE(CUTC_REP)
   IF(ALLOCATED(CUTA_ATT)) DEALLOCATE(CUTA_ATT)
   IF(ALLOCATED(CUTB_ATT)) DEALLOCATE(CUTB_ATT)
   IF(ALLOCATED(CUTC_ATT)) DEALLOCATE(CUTC_ATT)
   ALLOCATE(CUTA_REP(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(CUTB_REP(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(CUTC_REP(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(CUTA_ATT(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(CUTB_ATT(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(CUTC_ATT(NSPECIES(0),NSPECIES(0)))
   !
   CUTA_REP(:,:) = 0.0D0
   CUTB_REP(:,:) = 0.0D0
   CUTC_REP(:,:) = 0.0D0
   CUTA_ATT(:,:) = 0.0D0
   CUTB_ATT(:,:) = 0.0D0
   CUTC_ATT(:,:) = 0.0D0
   !
   DO J1=1,NSPECIES(0)
      DO J2=1,J1
         !
         !WRITE(*,*) "INITIALISATION> CUTOFF:", CUTOFF
         !
         ! Truncation coefficients for pairwise additive repulsion.
         ! EPSILON is absorbed here.
         DUMMY = MSC_EPS(J1,J2)*(MSC_A(J1,J2)/CUTOFF)**MSC_N(J1,J2)
         CUTA_REP(J1,J2) = -DUMMY
         CUTA_REP(J2,J1) = CUTA_REP(J1,J2)
         !
         DUMMY = DUMMY/CUTOFF
         CUTB_REP(J1,J2) = DBLE(MSC_N(J1,J2))*DUMMY
         CUTB_REP(J2,J1) = CUTB_REP(J1,J2)
         !
         DUMMY = DUMMY/CUTOFF
         CUTC_REP(J1,J2) = -DBLE(MSC_N(J1,J2)*(MSC_N(J1,J2)+1))*DUMMY
         CUTC_REP(J2,J1) = CUTC_REP(J1,J2)
         !
         ! Truncation coefficinets for many-body attraction.
         DUMMY = (MSC_A(J1,J2)/CUTOFF)**MSC_M(J1,J2)
         CUTA_ATT(J1,J2) = -DUMMY
         CUTA_ATT(J2,J1) = CUTA_ATT(J1,J2)
         !
         DUMMY = DUMMY/CUTOFF
         CUTB_ATT(J1,J2) = DBLE(MSC_M(J1,J2))*DUMMY
         CUTB_ATT(J2,J1) = CUTB_ATT(J1,J2)
         !
         DUMMY = DUMMY/CUTOFF
         CUTC_ATT(J1,J2) = -DBLE(MSC_M(J1,J2)*(MSC_M(J1,J2)+1))*DUMMY
         CUTC_ATT(J2,J1) = CUTC_ATT(J1,J2)
         !
      ENDDO
   ENDDO
   !
ENDIF
!
IF(HOMOREFT) THEN
   ALLOCATE(NNBRS(NATOMS,NSPECIES(0),0:NNMAX))
   ALLOCATE(ANNHIST(NSPECIES(0),NSPECIES(0),-1:NNMAX))
   ALLOCATE(NNHIST(NSPECIES(0),NSPECIES(0),-1:NNMAX))
   ALLOCATE(ANNHIST_MEAN(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(ANNHIST_VAR(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(NNBOND_WEIGHTS(NSPECIES(0),NSPECIES(0)))
   ALLOCATE(IFLIPE(1:NATOMS))
ENDIF
!
! ds656 > Allocate book-keeping arrays for listing
! atoms and their neighbours.  
!
ALLOCATE(ATOMLISTS(NSPECIES(0),3,0:NATOMS))
ALLOCATE(INVATOMLISTS(NATOMS,3))
!
IF(QALCS_SURFT .OR. QALCS_SYMT) THEN
   NNLIST_MAX = 12
   ALLOCATE(NNLISTS(1:NATOMS,0:NNLIST_MAX))
   IF(QALCS_SURFT) THEN
      NVSITES_MAX = 6*NATOMS !6*NINT(DBLE(NATOMS)**(2.0D0/3.0D0))
      ALLOCATE(VSITES(NVSITES_MAX,3))
   ENDIF
ENDIF
!
! set all atoms as mobile and partition by type
CALL RESET_ATOMLISTS(1) 
!
! ds656> If atomic labels are not to be permuted, we set QMINT 
! here and it will remain fixed throughout.
IF(.NOT. QALCST) THEN
   DO J1=1,NSAVE
      QMINT(J1,1:NATOMS) = INVATOMLISTS(1:NATOMS,1)
   END DO
ENDIF
!
! ds656> If species labels have been specified, the number of 
!        labels ought to correspond to the number of species!
IF(SPECLABELST) THEN
   IF(SIZE(SPECLABELS) /= NSPECIES(0)) THEN
      WRITE(MYUNIT,'(A)') &
           'initialization> SPECLABELS count != species count.'
      STOP
   ENDIF
ENDIF
!
IF(KEEPLISTS) THEN
   !
   ALLOCATE(NBRLISTS(2,NATOMS,0:NATOMS))
   ! Initialise neighbour lists to zero
   CALL RESET_NBRLISTS()
   !
ENDIF
!
IF (POLIRT) THEN
   CALL POLIRINIT
ENDIF
!
IF (MIEFT) THEN
   CALL MIEF_INI()
ENDIF
!
RETURN
!
END SUBROUTINE INITIALIZATIONS
