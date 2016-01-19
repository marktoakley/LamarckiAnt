!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales 
!   This file is part of GMIN.                                      
!
!   GMIN is free software; you can redistribute it and/or modIFy
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
!   along with this program; IF not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-130 US
!
!=============================================================
!   All routines in this file were implemented by
!   Dmitri Schebarchov (ds656).
!=============================================================
!
MODULE QALCS_SYM_MOD
  !
  ! ds656> This module contains all the variables and subroutines
  !        specific to QALCS_SYM.
  !
  USE COMMONS, ONLY : ATOMLISTS, NATOMS, NSPECIES, MYUNIT, DEBUG, &
       COORDS, PGSYMTOLS, QALCS_SYM_MINCORESIZE, INVATOMLISTS
  USE PGSYMMOD, ONLY : NSYMOPS_MAX
  !
  IMPLICIT NONE
  !
  INTEGER :: IPROC, INI_CORE_SIZE, NSITES, NORBITS, NSYMOPS_CORE
  INTEGER, ALLOCATABLE :: LIST_CORE(:), LIST_FREE(:)
  DOUBLE PRECISION :: TOLS(3), CNTR_CORE(3)
  DOUBLE PRECISION, ALLOCATABLE :: SYMOPS_CORE(:,:,:)
  !
  DOUBLE PRECISION, ALLOCATABLE :: LAT_COORDS(:,:)
  INTEGER, ALLOCATABLE :: LAT_LISTS(:,:), ORB_POP_COUNTS(:,:)
  !
CONTAINS
  !
  !==========================================================
  !ds656> First come subroutines called from QALCS_SYM
  !==========================================================
  !
  SUBROUTINE ALLOCATE_ALL()
    !
    IMPLICIT NONE
    !
    ALLOCATE(LIST_CORE(0:NATOMS))
    LIST_CORE(:) = 0
    ALLOCATE(LIST_FREE(0:NATOMS))
    LIST_FREE(:) = 0
    ALLOCATE(SYMOPS_CORE(NSYMOPS_MAX,3,3))
    SYMOPS_CORE(:,:,:) = 0.0D0
    ALLOCATE(LAT_COORDS(NATOMS,3))
    LAT_COORDS(:,:) = 0.0D0
    ALLOCATE(LAT_LISTS(NATOMS,0:NSYMOPS_MAX))
    LAT_LISTS(:,:) = 0
    ALLOCATE(ORB_POP_COUNTS(NATOMS,0:NSPECIES(0)))
    !
    RETURN
    !
  END SUBROUTINE ALLOCATE_ALL
  !
  SUBROUTINE DEALLOCATE_ALL()
    !
    IMPLICIT NONE
    !
    DEALLOCATE(LIST_CORE,LIST_FREE,SYMOPS_CORE,&
         LAT_COORDS,LAT_LISTS,ORB_POP_COUNTS)
    !
    RETURN
    !
  END SUBROUTINE DEALLOCATE_ALL
  !
  SUBROUTINE FIND_SYM_CORE(PROGRESS)
    !
    ! Routine for identifying a symmetric core in global COORDS.
    ! It also partitions non-core atom indices into radial shells.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: PROGRESS
    INTEGER :: I,J,K,LIST1(0:NATOMS), LIST2(0:NATOMS), NSYMOPS, &
         LABELS(1:NATOMS), NSPECIES_CORE(NSPECIES(0))
    DOUBLE PRECISION :: CNTR(3), X(3*NATOMS), &
         SYMOPS(NSYMOPS_MAX,3,3)
    CHARACTER(LEN=4) :: FPGRP
    !
    TOLS(1:3) = (/0.01D0, 0.2D0, 0.1D0/) !PGSYMTOLS(1:3) 
    PROGRESS = .FALSE.
    !
    NSYMOPS_CORE = 0
    INI_CORE_SIZE = LIST_CORE(0)
    !
    LIST1(0) = NATOMS; LIST2(0) = 0
    DO I=1,NATOMS
       LIST1(I) = I
       LIST2(I) = 0
    ENDDO
    !
    core_search: DO WHILE(LIST1(0) >= QALCS_SYM_MINCORESIZE)
       !
       ! Extract coordinates for current "core" list
       CALL COORDS_SUBSET(LIST1,COORDS(1:3*NATOMS,IPROC),X)
       !
       ! Create a list of labels for current "core" list
       LABELS(:) = 0
       NSPECIES_CORE(:) = 0
       DO I=1,LIST1(0)
          LABELS(I) = INVATOMLISTS(LIST1(I),1) 
          NSPECIES_CORE(LABELS(I)) = NSPECIES_CORE(LABELS(I)) + 1
       ENDDO
       !
       ! Determine majority label/species foc current "core" list
       J=1 
       DO I=2,NSPECIES(0)
          IF(NSPECIES_CORE(I) > NSPECIES_CORE(I-1)) J = I
       ENDDO
       !
       CALL PGSYM(LIST1(0),X(1:3*LIST1(0)),LABELS(1:LIST1(0)),&
            TOLS,2,J,CNTR,NSYMOPS,SYMOPS,FPGRP)
       WRITE(MYUNIT,'(A,I2,I5)') &
            'find_sym_core> NSYMOPS=',NSYMOPS
       !
       IF(NSYMOPS_CORE < NSYMOPS) THEN ! keep
          NSYMOPS_CORE = NSYMOPS
          SYMOPS_CORE(1:NSYMOPS,:,:) = SYMOPS(1:NSYMOPS,:,:)
          LIST_CORE(:) = LIST1(:)
          LIST_FREE(:) = LIST2(:)
          CNTR_CORE(1:3) = CNTR(1:3)
       ENDIF
       !       
       CALL STRIP_OUTER_SHELL(COORDS(1:3*NATOMS,IPROC),CNTR,&
            TOLS(2)/10.D0,LIST1,LIST2)
       !
    ENDDO core_search
    !
    IF(NSYMOPS_CORE > 0) THEN
       WRITE(MYUNIT,'(A,I5,A,I3,A)') &
            'find_sym_core> Found ',LIST_CORE(0), &
            ' -atom core with ',NSYMOPS_CORE,' sym. ops.'
       IF(NSYMOPS_CORE > 1 .AND. INI_CORE_SIZE < LIST_CORE(0)) &
            PROGRESS = .TRUE.
       IF(LIST_FREE(0) > 0) THEN 
          ! The list is back-to-front, reverse it now.
          ! is this step actually needed?
          LIST2(1:LIST_FREE(0)) = &
               LIST_FREE(LIST_FREE(0):1:-1)
          LIST_FREE(1:LIST_FREE(0)) = &
               LIST2(1:LIST_FREE(0))
       END IF
    ELSE
       WRITE(MYUNIT,'(A)') &
            'find_sym_core> Failed to find symmetry.'
       ! Put all atoms in core with zero symmetry elements.
       LIST_CORE(0) = NATOMS
       DO I=1,NATOMS
          LIST_CORE(I) = I
       ENDDO
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE FIND_SYM_CORE
  !
  SUBROUTINE BUILD_LATTICE(COMMENSURATE)
    !
    ! Generate a lattice by applying the core symmetry operations
    ! to non-core atoms.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: COMMENSURATE
    !
    LOGICAL :: SITE_NEW, SITE_EMPTY
    INTEGER :: I,J,K,L,LIST_SINGLE(0:NATOMS),L_IN_LIST(NATOMS),&
         ORB_SIZE, ORB_POP(NSPECIES(0),0:NATOMS)
    DOUBLE PRECISION :: X_IN_LIST(NATOMS,3), X(3), DX2, &
         ORB_COORDS(NSYMOPS_MAX,3), TOL2
    !
    ! 1) Re-sort LIST_FREE by distance from core centre.
    !    Initialise LIST_CLASH, lattice lists and coordinates.
    ! 2) Use sym. ops. to generate the orbit for each atom, 
    !    also checking if a site overlaps with another atom,
    !    the atom is added to the orbit's population.
    !    (NB: Need a different tolerance for overlap here!)
    ! 4) Add the orbit to the LATTICE data structure.
    !
    TOL2 = 1.0
    !
    CALL SORT_ATOMLIST_BY_DISTANCE_FROM_POINT(LIST_FREE,CNTR_CORE)
    !
    ! Translate shell atoms so that the origin is at core centroid
    DO I=1,LIST_FREE(0)
       L_IN_LIST(I) = INVATOMLISTS(LIST_FREE(I),1)
       K=3*(LIST_FREE(I)-1)
       DO J=1,3
          X_IN_LIST(I,J) = COORDS(K+J,IPROC) - CNTR_CORE(J)
       ENDDO
    ENDDO
    !
    NORBITS = 0
    NSITES = 0
    !
    I=0
    DO WHILE(I<LIST_FREE(0))
       !
       I=I+1
       ORB_SIZE = 1
       ORB_POP(:,:) = 0
       ORB_POP(L_IN_LIST(I),0) = 1 
       ORB_POP(L_IN_LIST(I),1) = LIST_FREE(I)
       ORB_COORDS(ORB_SIZE,1:3) = X_IN_LIST(I,1:3)
       !
       DO K=1,NSYMOPS_CORE ! Loop over sym.ops.
          !
          ! Apply sym. op. on current atom I.
          CALL PRODMM(3,3,1,SYMOPS_CORE(K,1:3,1:3),&
               X_IN_LIST(I,1:3),X(1:3))
          !
          SITE_NEW = .TRUE.; SITE_EMPTY = .TRUE.
          !
          ! Chek if atom I has been mapped onto itself or
          ! another atom index in LIST_FREE.
          J=I
          TEST1: DO WHILE(J<=LIST_FREE(0))
             DX2=0.0D0
             DO L=1,3
                DX2 = DX2 + (X(L) - X_IN_LIST(J,L))**2
             ENDDO
             IF(DX2 < TOL2) THEN
                IF(J == I) THEN
                   ! Atom I mapped onto itself
                   SITE_NEW = .FALSE.; SITE_EMPTY = .FALSE.
                ELSE
                   ! New symmetry-equivalent site overlaps with 
                   ! another atom.
                   SITE_EMPTY = .FALSE.
                   ORB_POP(L_IN_LIST(J),0) = &
                        ORB_POP(L_IN_LIST(J),0) + 1
                   ORB_POP(L_IN_LIST(J),ORB_POP(L_IN_LIST(J),0)) = &
                        LIST_FREE(J)
                   ! Move atom J out of (future) contention.
                   DO L=J,LIST_FREE(0)-1
                      LIST_FREE(L) = LIST_FREE(L+1)
                      X_IN_LIST(L,1:3) = X_IN_LIST(L+1,1:3)
                      L_IN_LIST(L) = L_IN_LIST(L+1)
                   ENDDO
                   LIST_FREE(LIST_FREE(0)) = 0
                   X_IN_LIST(LIST_FREE(0),1:3) = 0.0D0
                   L_IN_LIST(LIST_FREE(0)) = 0
                   LIST_FREE(0) = LIST_FREE(0) - 1
                ENDIF
                EXIT TEST1
             ENDIF
             J=J+1
          ENDDO TEST1
          !
          IF(SITE_NEW .AND. SITE_EMPTY) THEN
             ! Does the site clash with another?
             TEST2: DO J=1,ORB_SIZE
                DX2 = 0.0D0
                DO L=1,3
                   DX2 = DX2 + (X(L) - ORB_COORDS(J,L))**2
                ENDDO
                IF(DX2 < TOL2) THEN
                   SITE_NEW = .FALSE.
                   EXIT TEST2
                ENDIF
             ENDDO TEST2
          ENDIF
          !
          IF(SITE_NEW) THEN
             ORB_SIZE = ORB_SIZE + 1
             ORB_COORDS(ORB_SIZE,1:3) = X(1:3)
          ENDIF
          !
       ENDDO ! Looped over sym. ops. (i.e. created an orbit)
       !
       WRITE(MYUNIT,'(A,I3,A)',ADVANCE='NO') &
            'build_lattice> Orbit of size',&
            ORB_SIZE,' and population:'
       L = 0
       DO J=1,NSPECIES(0)
          L = L + ORB_POP(J,0)
          IF(ORB_POP(J,0) > 0) THEN
             WRITE(MYUNIT,'(I3,A,I1,A)',ADVANCE='NO') &
                  ORB_POP(J,0),' (T',J,'); total= '
          ENDIF
       ENDDO
       WRITE(MYUNIT,'(I3)') L
       !
       ! Now add current orbit to lattice... 
       NORBITS = NORBITS + 1
       ORB_POP_COUNTS(NORBITS,0) = L
       ORB_POP_COUNTS(NORBITS,1:NSPECIES(0)) = ORB_POP(1:NSPECIES(0),0)
       DO J=1,ORB_SIZE
          NSITES = NSITES + 1
          LAT_COORDS(NSITES,1:3) = ORB_COORDS(J,1:3)
          LAT_LISTS(NORBITS,0) = LAT_LISTS(NORBITS,0) + 1
          LAT_LISTS(NORBITS,LAT_LISTS(NORBITS,0)) = NSITES
       ENDDO
       !
    ENDDO ! Looped over free atoms.
    !
    RETURN
    !
  END SUBROUTINE BUILD_LATTICE
  !
  !==========================================================
  !ds656> Now come auxiliary subroutines (the nitty-gritty!)
  !==========================================================
  !
  SUBROUTINE SORT_ATOMLIST_BY_DISTANCE_FROM_POINT(LIST,CNTR)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT) :: LIST(0:NATOMS)
    DOUBLE PRECISION, INTENT(IN) :: CNTR(3)
    !
    INTEGER :: I,J,K
    DOUBLE PRECISION :: DISTS(0:NATOMS), X(3*NATOMS), R(3), DIST
    !
    CALL COORDS_SUBSET(LIST,COORDS(1:3*NATOMS,IPROC),X)
    !
    DISTS(:) = 0.0D0
    !
    ! Re-sort LIST_FREE by distance from CNTR_CORE.
    ! First -> shortest distance; Last -> longest distance.
    DO I=1,LIST(0)
       K=3*(LIST(I)-1)
       DO J=1,3
          R(J) = X(K+J) - CNTR(J)
       ENDDO
       ! Append atom index and sistance to the lists.
       DISTS(I) = DSQRT(DOT_PRODUCT(R,R))
       ! Now shuffle down the list, relying on DISTS(0)=0.0d0.
       J=I
       DO WHILE(DISTS(J) < DISTS(J-1))
          K=LIST(J); LIST(J)=LIST(J-1); LIST(J-1)=K
          DIST=DISTS(J); DISTS(J)=DISTS(J-1); DISTS(J-1)=DIST
          J=J-1
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE SORT_ATOMLIST_BY_DISTANCE_FROM_POINT
  !
  SUBROUTINE STRIP_OUTER_SHELL(X,CNTR,TOL,LIST1,LIST2)
    !
    ! From LIST1 remove a subset of (indices for) atoms that are
    ! furthest from the provided CNTR. Output the reduced LIST1,
    ! which will be sorted by distance! Also output LIST2, which
    ! will absorb entries removed from LIST1. Note that entries
    ! appended to LIST2 will also be sorted by distance, but in
    ! reverse order to LIST1.
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS), CNTR(3), TOL
    INTEGER, INTENT(INOUT) :: LIST1(0:NATOMS), LIST2(0:NATOMS)
    !
    INTEGER :: I, J, K
    DOUBLE PRECISION :: R(3), DIST, DISTS(0:LIST1(0))
    !
    DISTS(:) = 0.0D0
    !
    ! Sort LIST1 by distance
    DO I=1,LIST1(0)
       K=3*(LIST1(I)-1)
       DO J=1,3
          R(J) = X(K+J) - CNTR(J)
       ENDDO
       ! Append atom index and sistance to the lists.
       DISTS(I) = DSQRT(DOT_PRODUCT(R,R))
       ! Now shuffle down the list, relying on DISTS(0)=0.0d0.
       J=I
       DO WHILE(DISTS(J) < DISTS(J-1))
          K=LIST1(J); LIST1(J)=LIST1(J-1); LIST1(J-1)=K
          DIST=DISTS(J); DISTS(J)=DISTS(J-1); DISTS(J-1)=DIST
          J=J-1
       ENDDO
    ENDDO
    !
    ! Sequentially move tail-end entries from LIST1 to LIST2,
    ! until a 'gap' is found.
    K=0
    gap: DO I=LIST1(0)-1,1,-1
       K=K+1
       J=I+1
       LIST2(0) = LIST2(0) + 1
       LIST2(LIST2(0)) = LIST1(J)
       LIST1(J) = 0; LIST1(0) = I
       IF(DISTS(J)-DISTS(I) > TOL) EXIT gap
    ENDDO gap
    !
    WRITE(MYUNIT,'(A,I5)') 'strip_outer_shell> Stripped atoms:', K
    !
    RETURN
    !
  END SUBROUTINE STRIP_OUTER_SHELL
  !
  ! ds656> --- Test by visualising core and orbits ---
  !WRITE(*,'(A)') 'Core and orbits before merging:'
  !DO I = 1,LIST_CORE(0)
  !   J=3*(LIST_CORE(I)-1)
  !   WRITE(*,'(A,3(1X,F20.10))') 'O0', &
  !        (COORDS(J+K,IPROC), K=1,3)
  !ENDDO
  !DO I=1,NORBITS
  !   DO J = 1, ORB_SIZE(I)
  !      WRITE(*,'(A,I1,3(1X,F20.10))') 'O',I, &
  !           (ORB_COORDS(I,J,K) + CNTR_CORE(K), K=1,3)
  !   ENDDO
  !ENDDO
  ! <ds656 -------------------------------------------
  !
  SUBROUTINE COORDS_SUBSET(LIST,X0,X) 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: LIST(0:NATOMS)
    DOUBLE PRECISION, INTENT(IN) :: X0(1:3*NATOMS)
    DOUBLE PRECISION, INTENT(OUT) :: X(1:3*NATOMS)
    !
    INTEGER :: I,J,K,L
    !
    X(1:3*NATOMS) = 0.0D0
    DO I=1,LIST(0)
       J=3*(LIST(I)-1)
       K=3*(I-1)
       DO L=1,3
          X(K+L) = X0(J+L)
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE COORDS_SUBSET
  !
END MODULE QALCS_SYM_MOD
