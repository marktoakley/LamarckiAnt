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
SUBROUTINE QALCS_SURF(NP,ITER,TIME,BRUN,QDONE,SCREENC)
  !
  ! Quench-Assisted Local Combinatorial Search for surface vacancies.
  ! Atoms are systematically swapped with "vacancies" (type 0).
  !
  USE COMMONS, ONLY : NATOMS, VSITES, NNLISTS, COORDS, &
       NQ, ECONV, MYUNIT, TSTART, QALCSV, BOXCENTROIDT
  !
  IMPLICIT NONE
  !
  ! Parse passed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  LOGICAL ::  COMPLETE
  INTEGER :: I1,I2,I3,I33,K,NQTOT,ATOMS_SORTED_BY_NN(0:NATOMS), &
       NN_LOWEST, NVSITES, VSITE_WEIGHT
  DOUBLE PRECISION :: X(3*NATOMS),E,POTEL,NN_AVE, NNDIST(3)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  X(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
  E = POTEL
  COMPLETE = .FALSE.
  !
  WRITE(MYUNIT,'(A,F20.10,A,F11.1)') &
       'QALCS_surf> Initial E= ',POTEL,' t= ',TIME-TSTART
  !
  DO WHILE (.NOT. COMPLETE)
     !
     COMPLETE = .TRUE.
     !
     CALL BUILD_NNLISTS(NP, ATOMS_SORTED_BY_NN, NNDIST)
     CALL GEN_VSITES(NP, ATOMS_SORTED_BY_NN, NNDIST, &
          NVSITES, VSITE_WEIGHT) 
     !
     ! Consider only atoms of coordination 6 and below
     ! NN_AVE is apparently too large as an upper bound... 
     ! Perhaps try NN_LOWEST? (By analogy with MAX_WEIGHT for voids.)
     NN_LOWEST = NNLISTS(ATOMS_SORTED_BY_NN(1),0)
     !write(*,*) NN_LOWEST
     IF(NN_LOWEST <= VSITE_WEIGHT) THEN
        REDUCE: DO I2=1,ATOMS_SORTED_BY_NN(0) 
           I3=ATOMS_SORTED_BY_NN(I2) ! actual atom index
           IF(NNLISTS(I3,0) > NN_LOWEST) THEN
              ATOMS_SORTED_BY_NN(0) = I2-1
              EXIT REDUCE
           ENDIF
        ENDDO REDUCE
     ELSE
        ATOMS_SORTED_BY_NN(0) = 0
     ENDIF
     !
     WRITE(MYUNIT,'(A,I3,A,I3,A)') &
          'QALCS_surf> Swapping ',ATOMS_SORTED_BY_NN(0), &
          ' atom(s) with ',NVSITES,' v-site(s).'
     !
     SITES: DO I1=1,NVSITES ! loop thru best vacant sites
        !
        DO I2=1,ATOMS_SORTED_BY_NN(0) ! loop thru worst atoms
           !
           I3=ATOMS_SORTED_BY_NN(I2) ! actual atom index
           !
           I33=3*(I3-1)
           !
           DO K=1,3
              COORDS(I33+K,NP) = VSITES(I1,K)
           ENDDO
           !
           IF(BOXCENTROIDT) CALL BOXCENTROID(COORDS(:,NP))
           !
           NQTOT = NQTOT + 1
           NQ(NP) = NQ(NP) + 1
           CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
           IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '  ')
           !
           IF(POTEL < E - ECONV) THEN ! improvemenet
              COMPLETE = .FALSE.
              E = POTEL
              X(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
              WRITE(MYUNIT,'(A,F20.10)') &
                   'QALCS_surf> Found lower E= ', E
              EXIT SITES
           ELSE ! revert
              POTEL = E
              COORDS(1:3*NATOMS,NP) = X(1:3*NATOMS)
           ENDIF
           !
        ENDDO
     ENDDO SITES
     !
  END DO
  !
  WRITE(MYUNIT,'(A,F20.10,A,F11.1)') &
       'QALCS_surf> Final E= ',POTEL,' t= ',TIME-TSTART  
  !
  RETURN
  !
END SUBROUTINE QALCS_SURF
!
!=================================================================
!
SUBROUTINE BUILD_NNLISTS(NP, LIST, NNDIST)
  !
  ! This non-standard way of builiding NNLISTS is tailored
  ! for generating vacancies on a surface. One noteworthy side-
  ! effect is that, for a pair of atoms (I and J), it can happen 
  ! that I is listed as a neighbour of J, but J is not listed as a 
  ! neighbour of I.
  !
  ! Although the procedure does not distinguish between different
  ! atomic species/types, it does respect atom groups:
  ! the nearest-neighboour analysis is restricted to atoms in groups
  ! 1 and 2; atoms in group 3 are completely ignored. Further,
  ! only group-1 atoms will be listed in LIST on output; and 
  ! the minimum, average and maximum NN distances stored in NNDIST
  ! will be computed exclusively from the local environment of 
  ! atoms in group 1.
  !
  USE COMMONS, ONLY : NNLISTS,ATOMLISTS,INVATOMLISTS,NATOMS, & 
       NSPECIES, COORDS, MYUNIT, NNLIST_MAX
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(OUT) :: LIST(0:NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: NNDIST(3)
  !
  LOGICAL :: GROWN, GAP, SELF
  INTEGER :: T1,T2,G1,G2,J1,J2,I1,I2,I13,I23,K,L,M,N, IT
  INTEGER :: LJ1,UJ1,LJ2,UJ2,G2START,IPAIR1(2),IPAIR2(2),NN_SURF
  DOUBLE PRECISION :: DUMMY,DIST2,GAPTOL,RNGTOL,NN_AVE, &
       NNDISTS(NATOMS,NNLIST_MAX)
  
  !
  NNLISTS(1:NATOMS,0:NNLIST_MAX) = 0
  NNDISTS(1:NATOMS,1:NNLIST_MAX) = 0.0D0
  !
  ! Set the tolerance for defining a gap in the distribution of
  ! nearest-neighbour distances. GAPTOL = 1.1
  ! maximum deviation from minimal NN distance: RNGTOL=1.457~((1+SQRT(2))/2)^2
  !GAPTOL = 1.457D0 ! 15% jump?
  !RNGTOL = 2.0D0 ! 1+0.75*fccgap
  GAPTOL = 1.333D0
  RNGTOL = 1.457d0
  NNDIST(1) = 1.0D9 ! initialise minimal value squared NN distance
  !
  ! 1st pass: fill and sort lists
  !
  DO T1=1,NSPECIES(0) ! Double loop over ALL types
     DO T2=1,T1
        !
        DO G1=1,2 ! Double loop over groups 1 & 2 (not 3)
           !
           IF(T1 == T2) THEN
              G2START = G1
           ELSE
              G2START = 1
           ENDIF
           !
           DO G2 = G2START,2
              !
              ! Set bounds for outer loop over atoms
              LJ1=1
              UJ1=ATOMLISTS(T1,G1,0)
              IF(T1==T2.AND.G1==G2) THEN
                 SELF = .TRUE.
                 UJ1 = UJ1 - 1 ! Exclude distance to itself
              ELSE
                 SELF = .FALSE.
              ENDIF
              !
              DO J1=LJ1,UJ1 ! Double loop over atoms
                 I1=ATOMLISTS(T1,G1,J1)
                 I13 = 3*(I1-1)
                 !
                 ! Set bounds for inner loop over atoms
                 IF(SELF) THEN
                    LJ2=J1+1
                    UJ2=UJ1+1
                 ELSE
                    LJ2=1
                    UJ2=ATOMLISTS(T2,G2,0)
                 ENDIF
                 !
                 DO J2=LJ2,UJ2
                    I2=ATOMLISTS(T2,G2,J2)
                    I23 = 3*(I2-1)
                    !
                    IPAIR1=(/I1,I2/) ! pair to be swapped
                    IPAIR2=(/I2,I1/) ! useful later
                    !
                    DIST2=0.0D0
                    DO K=1,3
                       DIST2 = DIST2 + &
                            (COORDS(I23+K,NP)-COORDS(I13+K,NP))**2
                    ENDDO
                    IF(DIST2<NNDIST(1)) NNDIST(1) = DIST2
                    !
                    DO K=1,2 ! loop through the pair of lists
                       ! Check if list is still growing
                       IF(NNLISTS(IPAIR1(K),0) < NNLIST_MAX) THEN
                          NNLISTS(IPAIR1(K),0) = &
                               NNLISTS(IPAIR1(K),0) + 1
                          GROWN = .TRUE.
                       ELSEIF( DIST2 < NNDISTS( IPAIR1(K), &
                            NNLISTS(IPAIR1(K),0) ) ) THEN
                          GROWN = .TRUE.
                       ELSE
                          GROWN = .FALSE.
                       ENDIF
                       ! Update the last entry if required
                       IF( GROWN ) THEN
                          M = NNLISTS(IPAIR1(K),0)
                          NNLISTS(IPAIR1(K),M) = IPAIR2(K)
                          NNDISTS(IPAIR1(K),M) = DIST2
                          !
                          !ds656> testing...
                          !write(MYUNIT,*) &
                          !     'build_nnlists> list for atom',&
                          !     IPAIR1(K),' : ',&
                          !     (NNLISTS(IPAIR1(K),IT), &
                          !     IT=1,NNLISTS(IPAIR1(K),0) )
                          !write(MYUNIT,*) &
                          !     'build_nnlists> dists: ',&
                          !     (NNDISTS(IPAIR1(K),IT), &
                          !     IT=1,NNLISTS(IPAIR1(K),0) )
                          !<ds656 ...testing
                          ! Shuffle the last entry if possible
                          SHUFFLE: DO L=NNLISTS(IPAIR1(K),0),2,-1
                             M=L-1
                             IF(NNDISTS(IPAIR1(K),L) < &
                                  NNDISTS(IPAIR1(K),M)) THEN
                                ! Swap atom indices in list
                                N = NNLISTS(IPAIR1(K),M)
                                NNLISTS(IPAIR1(K),M) = &
                                     NNLISTS(IPAIR1(K),L) 
                                NNLISTS(IPAIR1(K),L) = N
                                ! Swap distances
                                DUMMY = NNDISTS(IPAIR1(K),M)
                                NNDISTS(IPAIR1(K),M) = &
                                     NNDISTS(IPAIR1(K),L) 
                                NNDISTS(IPAIR1(K),L) = DUMMY
                             ELSE
                                EXIT SHUFFLE
                             ENDIF
                          ENDDO SHUFFLE
                          !
                       ENDIF ! check for change in list
                       !
                    ENDDO ! Looped over two lists
                    !
                 ENDDO ! Double loop over atom indices
              ENDDO
              !
           ENDDO ! Double loop over groups
        ENDDO
        !
     ENDDO ! Double loop over types
  ENDDO
  !
  !ds656> testing...
  !WRITE(MYUNIT,*) 'gen_nnlists> NN counts after 1st pass:', &
  !     (NNLISTS(J1,0), J1 = 1, NATOMS)
  !DO J1=1,NATOMS
  !   write(MYUNIT,'(A,I3,A)',ADVANCE='NO') &
  !        'build_nnlists> List for ', J1,':'
  !   DO J2=1,NNLISTS(J1,0)
  !      write(MYUNIT,'(I3,A,F6.4,A)',ADVANCE='NO')  &
  !           NNLISTS(J1,J2), '(', NNDISTS(J1,J2),')'
  !   ENDDO
  !   WRITE(MYUNIT,*)
  !ENDDO
  !<ds656 ...testing
  !
  ! 2nd pass: remove comparatively distant neighbours from lists 
  L=0
  M=0
  LIST(0:NATOMS) = 0 ! initialise the sorted list of atom labels
  NNDIST(2:3) = 0.0D0 ! initialise average and maximium values
  DO J1=1,NATOMS
     IF(INVATOMLISTS(J1,2) == 1) THEN ! consider group 1 only
        ! Check that the current atom's nearest neighbour 
        ! distance is not much greater than the overall 
        ! smallest neighbour distance. 
        IF(NNDISTS(J1,1)/NNDIST(1) < RNGTOL) THEN
           GAP=.FALSE.
           DUMMY = NNDISTS(J1,1)
           GAPSEARCH: DO J2=2,NNLISTS(J1,0)
              ! Identify a gap by spanning the ordered list of 
              ! NNDISTs (*squared* distances) and locating the 
              ! first gap greater than GAPTOL
              ! [ Note: RNGTOL = 1.457 ~ ((1+SQRT(2))/2)^2 ]
              !IF(NNDISTS(J1,J2)/NNDISTS(J1,J2-1) < GAPTOL .AND. &
              !     NNDISTS(J1,J2)/NNDISTS(J1,1) < RNGTOL) THEN   
              IF(NNDISTS(J1,J2)/NNDISTS(J1,1) < RNGTOL) THEN     
                 DUMMY = DUMMY + NNDISTS(J1,J2)
                 IF(NNDISTS(J1,J2) > NNDIST(3)) &
                      NNDIST(3) = NNDISTS(J1,J2)
              ELSE
                 GAP=.TRUE.
                 K = J2
                 EXIT GAPSEARCH
              ENDIF
           ENDDO GAPSEARCH
        ELSE
           GAP=.TRUE.
           K=1
        ENDIF
        IF(GAP) THEN
           DO J2=K,NNLISTS(J1,0)
              NNLISTS(J1,J2) = 0
              NNDISTS(J1,J2) = 0.0D0
           ENDDO
           K = K - 1
           ! Update count of neighbours for atom
           NNLISTS(J1,0) = K
        ELSE
           K = NNLISTS(J1,0)
        ENDIF
        ! Accumulate average coordination
        L = L + K
        ! Accumulate average distance
        NNDIST(2) = NNDIST(2) + DUMMY
        ! Add atom index in LIST, sorted by coordination
        ! from lowest to highest.
        M = M + 1
        LIST(M) = J1
        ORDER: DO J2=M,2,-1
           IF(NNLISTS(LIST(J2),0) < NNLISTS(LIST(J2-1),0)) THEN
              K = LIST(J2)
              LIST(J2) = LIST(J2-1)
              LIST(J2-1) = K
           ELSE
              EXIT ORDER
           ENDIF
        ENDDO ORDER
     !ELSE ! Set to zero and ignore in the future
     !   NNLISTS(J1,0:NNLIST_MAX) = 0
     !   NNDISTS(J1,1:NNLIST_MAX) = 0.0D0
     ENDIF ! Group check
  ENDDO
  !
  NNDIST(2) = NNDIST(2)/DBLE(L)
  NN_AVE = DBLE(L)/DBLE(M)
  !
  !ds656> testing...
  ! WRITE(MYUNIT,'(A)') 'build_nnlists> Final NN counts:'
  ! DO K=1,NATOMS
  !      WRITE(MYUNIT,'(I3)',ADVANCE='NO') NNLISTS(K,0)
  ! ENDDO
  ! WRITE(MYUNIT,'(A)') ''
  !WRITE(MYUNIT,'(A)') 'gen_nnlists> Atoms sorted by NN count:'
  !do J1=1,M ! M is the length of LIST set in the previous loop
  !   j2= list(j1)
  !   write(myunit,'(A,3(I5))') 'gen_nnlists> ',j1,j2, nnlists(j2,0)
  !enddo
  !<ds656 ...testing
  !
  ! Choose cutoff for defining 'surface' atoms.
  !NN_SURF = NINT(NN_AVE)
  NN_SURF = 11 !?! Or NNLIST_MAX-1 ?!?
  LIST(0) = M
  SURFACE: DO K=1,M ! M <= NATOMS is the length of LIST
     J1=LIST(K) ! Actual atom index
     IF(NNLISTS(J1,0) > NN_SURF) THEN
        LIST(0) = K-1 ! store count for 'surface' atoms
        EXIT SURFACE 
     ENDIF      
  ENDDO SURFACE
  !
  WRITE(MYUNIT,'(A,I3,1X,F6.3)') &
       'build_nnlists> MIN and AVE NN-counts:', &
       NNLISTS(LIST(1),0), NN_AVE
  WRITE(MYUNIT,'(A,3(1X,F8.5))') &
       'build_nnlists> MIN, AVE, MAX NN-distances: ', &
       (SQRT(NNDIST(K)), K=1,3)
  !
  RETURN
  !
END SUBROUTINE BUILD_NNLISTS
!
!=================================================================
!
SUBROUTINE GEN_VSITES(NP, LIST, NNDIST, NVSITES, VSITE_WEIGHT)
  !
  ! Generate vacant sites from NNLISTS
  ! USE COMMONS, ONLY : NATOMS, NNLISTS, COORDS, NNLIST_MAX, &
  !     VSITES, NVSITES_MAX, MYUNIT, ATOMLISTS, NSPECIES
  USE COMMONS, ONLY : NATOMS, NNLISTS, COORDS, NNLIST_MAX, &
      VSITES, NVSITES_MAX, MYUNIT, MIEFT, MIEF_NSITES, &
      MIEF_SITES, MIEF_SIG, ATOMLISTS, NSPECIES
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(IN) :: LIST(0:NATOMS)
  DOUBLE PRECISION, INTENT(IN) :: NNDIST(3)
  INTEGER, INTENT(OUT) :: NVSITES, VSITE_WEIGHT
  INTEGER :: I0,I1,I2,I3,I13,I23,I33,J1,J2,J3,K, &
       VSITE_WEIGHTS(NVSITES_MAX), NNLIST(0:NNLIST_MAX), &
       VSITE_NNLIST(0:NNLIST_MAX), n0, n1, n2, n3
  LOGICAL :: BAD, ADD
  DOUBLE PRECISION :: VSITE(3), DUMMY, TOL1, TOL2, TOL3, TOL4
  !
  VSITES(1:NVSITES_MAX,1:3) = 0.0D0
  NVSITES=0
  VSITE_WEIGHTS(1:NVSITES_MAX) = 0
  !
  ! Thresholds for cluster-cluster atoms
  !TOL1=0.375d0*NNDIST(1) ! squared rad. of circumsphere of tetrahedron (for overalap)
  TOL1=NNDIST(1)*(0.5d0*(1.0D0+DSQRT(3.0D0/8.0D0)))**2 
  TOL2=NNDIST(3) ! Coordination determined by the largest NN distance
  !
  ! Thresholds for cluster-substrate atoms
  IF(MIEFT) THEN
     TOL3=0.64D0*MINVAL(MIEF_SIG)**2 ! adequate?!?
     TOL4=1.44D0*MAXVAL(MIEF_SIG)**2
  ENDIF
  !
  DO I0=1,LIST(0) ! loop over atoms in list
     !
     I1=LIST(I0) ! actual atom index
     I13=3*(I1-1)
     !
     DO J1=1,NNLISTS(I1,0) ! max = NNLIST_MAX
        I2=NNLISTS(I1,J1)
        I23=3*(I2-1)
        !
        !WRITE(MYUNIT,'(A)') &
        !     'gen_sites> New candidate v-site'
        !
        ! Generate candidate vsite
        DO K=1,3 ! XYZ coordinates for potential vsite
           VSITE(K) = &
                2.0D0*COORDS(I13+K,NP) - COORDS(I23+K,NP)
        ENDDO
        VSITE_NNLIST(0) = 1
        VSITE_NNLIST(1) = I1
        VSITE_NNLIST(2:NNLIST_MAX) = 0
        !
        BAD=.FALSE.
        !
        ! Check distance from vsite to other neighbouring atoms
        BADNESS1: DO J2=1,NNLISTS(I1,0)
           IF(J2 /= J1) THEN
              I2=NNLISTS(I1,J2)
              I23=3*(I2-1)
              DUMMY = 0.0D0
              DO K=1,3
                 DUMMY = DUMMY + &
                      (COORDS(I23+K,NP)-VSITE(K))**2
              ENDDO
              IF(DUMMY < TOL1) THEN ! overlap
                 BAD = .TRUE.
                 EXIT BADNESS1
              ELSEIF(DUMMY < TOL2) THEN ! add as neighbour
                 IF(VSITE_NNLIST(0) < NNLIST_MAX) THEN 
                    VSITE_NNLIST(0) = VSITE_NNLIST(0) + 1
                    VSITE_NNLIST(VSITE_NNLIST(0)) = I2
                    !WRITE(MYUNIT,'(A,I3,A,f10.5)') &
                    !     'gen_sites> Added atom ',I2,' as nbro.',dsqrt(dummy)
                 ELSE
                    !WRITE(MYUNIT,'(A)') &
                    !     'gen_sites> WARNING(1): exceeded max nbro count'
                    !WRITE(MYUNIT,'(A,I3,A,f10.5)') &
                    !     'gen_sites> Not added atom ',I2,' as nbro.',dsqrt(dummy)
                 ENDIF
                 ! --------------------------------------------------
                 ! Also consider nbros of I2 that are not nbros of I1
                 ! --------------------------------------------------
                 CALL REDUCE_NNLIST(NNLIST_MAX, &
                      NNLISTS(I1,0:NNLIST_MAX), &
                      NNLISTS(I2,0:NNLIST_MAX), &
                      NNLIST(0:NNLIST_MAX) )
                 DO J3=1,NNLIST(0)
                    I3=NNLIST(J3)
                    I33=3*(I3-1)
                    DUMMY = 0.0D0
                    DO K=1,3
                       DUMMY = DUMMY + &
                            (COORDS(I33+K,NP)-VSITE(K))**2
                    ENDDO
                    IF(DUMMY < TOL1) THEN ! overlap!
                       ! This block turns out to be important, especially when
                       ! distinction between 1st and 2nd nearest neighbours is
                       ! ambiguous!
                       BAD=.TRUE.
                       EXIT BADNESS1
                    ELSEIF( DUMMY < TOL2 ) THEN
                       ADD=.TRUE.
                       LK: DO K=1,VSITE_NNLIST(0)
                          IF(VSITE_NNLIST(K) == I3) THEN
                             ADD = .FALSE.
                             EXIT LK
                          ENDIF
                       ENDDO LK
                       IF(ADD) THEN 
                          IF (VSITE_NNLIST(0) < NNLIST_MAX) THEN
                             VSITE_NNLIST(0) = VSITE_NNLIST(0) + 1
                             VSITE_NNLIST(VSITE_NNLIST(0)) = I3
                             !WRITE(MYUNIT,'(A,I3,A,f10.5)') &
                             !     'gen_sites> Added atom ',I3,' as nbro.',dsqrt(dummy)
                          !ELSE
                             !WRITE(MYUNIT,'(A)') &
                             !     'gen_sites> WARNING(2): exceeded max nbro count'
                             !WRITE(MYUNIT,'(A,I3,A,f10.5)') &
                             !     'gen_sites> Not added atom ',I3,' as nbro.',dsqrt(dummy)

                          ENDIF
                       ENDIF
                    ENDIF
                 ENDDO
                 ! --------------------------------------------------
              ENDIF
              !
           ENDIF
        ENDDO BADNESS1
        !
        ! Check distance from VSITE to Mie site(s)
        VSITE_WEIGHT = 0
        IF(MIEFT .AND. .NOT. BAD) THEN
           BADNESS2: DO J2=1,MIEF_NSITES
              DUMMY = 0.0D0
              DO K=1,3
                 DUMMY = DUMMY + &
                      (MIEF_SITES(J2,K)-VSITE(K))**2
              ENDDO
              IF(DUMMY < TOL3) THEN ! overlap
                 BAD = .TRUE.
                 EXIT BADNESS2
              ELSEIF(DUMMY < TOL4) THEN ! add extra weight 
                 VSITE_WEIGHT = VSITE_WEIGHT + 1
              ENDIF
           ENDDO BADNESS2
        ENDIF
        !
        IF(.NOT. BAD) THEN
           IF(NVSITES < NVSITES_MAX) THEN
              NVSITES=NVSITES+1
              VSITES(NVSITES,1:3) = VSITE(1:3)
              !test
              !WRITE(MYUNIT,'(A,3(1X,F10.5))') 'gen_sites> NEW V-site:',&
              !     VSITE(1:3) 
              VSITE_WEIGHTS(NVSITES) = VSITE_NNLIST(0) + VSITE_WEIGHT
           ELSE
              WRITE(MYUNIT,'(A,I6)') &
                   'gen_sites> NVSITES exceeded ', NVSITES_MAX
              STOP
           ENDIF
           IF(VSITE_WEIGHTS(NVSITES).GE.12) THEN
              WRITE(MYUNIT,'(A,I3)') 'gen_vsites> ***WARNING*** v-site coordination',&
                   VSITE_WEIGHTS(NVSITES)
              ! ds656> testing ..
              ! WRITE(MYUNIT,'(A)', advance='no') 'gen_vsites> v-site nbros:'
!               do n0=1,vsite_nnlist(0)
!                  WRITE(MYUNIT,'(1x,i3)', advance='no') vsite_nnlist(n0)
!               enddo
!               write(myunit,*)

!               WRITE(MYUNIT,'(A,I3)') 'gen_vsites> v-site generated from nbro of', I1
!               WRITE(MYUNIT,'(A,I3)') 'gen_vsites> the nbro is atom',nnlists(i1,j1)
!               WRITE(MYUNIT,'(A)', advance='no') 'gen_vsites> all nbros:'
!               do n0=1,nnlists(i1,0)
!                  WRITE(MYUNIT,'(1x,i3)', advance='no') nnlists(i1,n0)
!               enddo
!               write(myunit,*)

!               write(*,'(I5)') natoms+1
!               write(*,'(A)') 'Atoms and vacancies before merging'
!               DO n0=1,NSPECIES(0) !ds656> Should not exceed 10           
!                  DO n1=1,ATOMLISTS(n0,1,0)
!                     n2=ATOMLISTS(n0,1,n1)
!                     WRITE(*,'(I2,3(1X,F20.10))') n0,&
!                          (COORDS(n3,NP), n3=3*(n2-1)+1,3*n2)
!                  ENDDO
!               ENDDO
!               WRITE(*,'(A,3(1X,F20.10))') ' 0',(VSITE(J2), J2=1,3)
!              stop
              ! <ds656 .. testing
           ENDIF
        ENDIF
        !
     ENDDO ! loop over neighbour list
     !
  ENDDO ! loop over atoms
  !
  !ds656> testing...
  !WRITE(MYUNIT,'(A)') 'gen_vsites> unsorted Vsites:'  
  !DO I1=1,NVSITES
  !   WRITE(MYUNIT,'(A,I6,I3,3(1X,F10.5))') 'gen_sites>', &
  !        I1, VSITE_WEIGHTS(I1), VSITES(I1,1:3)
  !ENDDO
  !<ds656 ...testing
  !
  !WRITE(MYUNIT,'(A,I6)') 'gen_sites> V-site candidates: ', NVSITES
  !
  ! Sort VSITES by weight from biggest to smallest (CHECK!)
  CALL SORT3_V2(NVSITES,NVSITES_MAX,VSITE_WEIGHTS,VSITES)
  ! weed out vacancies whose avergage coordination is below 
  !
  !ds656> testing...
  !WRITE(MYUNIT,'(A)') 'gen_vsites> Vsites sorted by occupancy:'  
  !DO I1=1,NVSITES
  !   WRITE(MYUNIT,'(A,I6,I3,3(1X,F10.5))') 'gen_sites>', &
  !        I1, VSITE_WEIGHTS(I1), VSITES(I1,1:3)
  !ENDDO
  ! write(*,'(I5)') natoms+nvsites
  ! write(*,'(A)') 'Atoms and all 1st generation vacancies'
  ! DO I0=1,NSPECIES(0) !ds656> Should not exceed 10           
  !    DO I1=1,ATOMLISTS(I0,1,0)
  !       J1=ATOMLISTS(I0,1,I1)
  !       WRITE(*,'(A,I1,3(1X,F20.10))') 'A',I0,&
  !            (COORDS(J2,NP), J2=3*(J1-1)+1,3*J1)
  !    ENDDO
  ! ENDDO
  ! DO I1=1,NVSITES
  !    WRITE(*,'(I2,3(1X,F20.10))') VSITE_WEIGHTS(I1),&
  !         (VSITES(I1, J2), J2=1,3)
  ! ENDDO
  !<ds656 ...testing
  !
  I2=NVSITES
  IF(VSITE_WEIGHTS(1) > 1) THEN ! or 1?
     ! Weed out all but the largest WEIGHTS (highest coordination)
     WEED: DO I1=I2,1,-1
        IF(VSITE_WEIGHTS(I1) < VSITE_WEIGHTS(1)) THEN
           NVSITES = NVSITES-1
           VSITE_WEIGHTS(I1) = 0
           VSITES(I1,1:3) = 0.0D0
        ELSE
           EXIT WEED
        ENDIF
     ENDDO WEED
     VSITE_WEIGHT = VSITE_WEIGHTS(1)
     !WRITE(MYUNIT,'(A,I2)') &
     !     'gen_sites> Highest v-site weight: ',VSITE_WEIGHT
  ELSE
     NVSITES = 0
  ENDIF
  !
  ! ds656> testing ..
  !if(vsite_weights(1) .ge. 12) then
     ! write(*,'(I5)') natoms+nvsites
     ! write(*,'(A)') 'Atoms and vacancies before merging'
     ! DO I0=1,NSPECIES(0) !ds656> Should not exceed 10           
     !    DO I1=1,ATOMLISTS(I0,1,0)
     !       J1=ATOMLISTS(I0,1,I1)
     !       WRITE(*,'(I2,3(1X,F20.10))') I0,&
     !            (COORDS(J2,NP), J2=3*(J1-1)+1,3*J1)
     !    ENDDO
     ! ENDDO
     ! DO I1=1,NVSITES
     !    WRITE(*,'(A,3(1X,F20.10))') ' 0',(VSITES(I1, J2), J2=1,3)
     ! ENDDO
     !stop
  !endif
  ! <ds656 .. testing
  !
  I1 = 1 
  DO WHILE(I1 < NVSITES)
     I2 = I1+1
     DO WHILE(I2 <= NVSITES)
        DUMMY = 0.0D0
        DO K=1,3
           DUMMY = DUMMY + &
                (VSITES(I1,K)-VSITES(I2,K))**2
        ENDDO
        IF(DUMMY < TOL1) THEN ! Absorb I2 into I2 and decrement NVSITES
           VSITES(I1,:) = 0.5D0*( VSITES(I1,:) + VSITES(I2,:) )
           VSITES(I2,:) = VSITES(NVSITES,:)
           VSITE_WEIGHTS(I2) = VSITE_WEIGHTS(NVSITES)
           VSITES(NVSITES,:) = 0.0D0
           VSITE_WEIGHTS(NVSITES) = 0
           NVSITES = NVSITES - 1
        ELSE ! update I2
           I2 = I2 + 1
        ENDIF
     ENDDO
     I1 = I1 + 1
  ENDDO
  !
  IF(NVSITES > 0) THEN
     WRITE(MYUNIT,'(A,I3,A,I3)') &
          'gen_sites> Generated ', NVSITES,' v-site(s) of weight',&
          VSITE_WEIGHT
  ELSE
     WRITE(MYUNIT,'(A,I6)') &
          'gen_sites> Failed to generate adequate v-sites!'     
  ENDIF

  !ds656> testing ..
  !if(vsite_weight == 12) then
  ! write(*,'(I5)') natoms+nvsites
  ! write(*,'(A)') 'Final atoms and vacancies after merging'
  ! DO I0=1,NSPECIES(0) !ds656> Should not exceed 10           
  !    DO I1=1,ATOMLISTS(I0,1,0)
  !       J1=ATOMLISTS(I0,1,I1)
  !       WRITE(*,'(I2,3(1X,F20.10))') I0,&
  !            (COORDS(J2,NP), J2=3*(J1-1)+1,3*J1)
  !    ENDDO
  ! ENDDO
  ! DO I1=1,NVSITES
  !    WRITE(*,'(A,3(1X,F20.10))') '0',(VSITES(I1, J2), J2=1,3)
  ! ENDDO
  ! stop
  !endif
  !<ds656 .. testing
  !
  RETURN
  !
END SUBROUTINE GEN_VSITES
!
SUBROUTINE REDUCE_NNLIST(NNLIST_MAX,LISTA,LISTB,LISTC)
  !
  ! Put in LISTC all entries from LISTB that are not in LISTA.
  ! Lists are not assumed to be sorted, so the procedure is slow.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NNLIST_MAX
  INTEGER, INTENT(IN) :: LISTA(0:NNLIST_MAX), LISTB(0:NNLIST_MAX)
  INTEGER, INTENT(OUT) :: LISTC(0:NNLIST_MAX)
  INTEGER :: IA, IB 
  LOGICAL :: SHARED
  !
  LISTC(0:NNLIST_MAX) = 0
  !
  LB: DO IB=1,LISTB(0)
     SHARED = .FALSE.
     LA: DO IA = 1, LISTA(0)
        IF(LISTB(IB) == LISTA(IA)) THEN
           SHARED = .TRUE.
           EXIT LA
        ENDIF
     ENDDO LA
     IF(.NOT. SHARED) THEN
        LISTC(0) = LISTC(0) + 1
        LISTC(LISTC(0)) = LISTB(IB)
     ENDIF
  ENDDO LB
  !
  RETURN
  !
END SUBROUTINE REDUCE_NNLIST
