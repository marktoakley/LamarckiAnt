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
SUBROUTINE QALCS(NP,ITER,TIME,BRUN,QDONE,SCREENC)
  !
  ! Quench-Assisted Local Combinatorial Search
  !
  USE COMMONS, ONLY : NATOMS,COORDS,ECONV,NSPECIES,TSTART,MYUNIT,&
                      INVATOMLISTS,QALCS_SURFT,QALCST,QALCSMODE,&
                      QALCS_SYMT, SEQLENGTH, SAVEMULTIMINONLY,&
                      LFLIPS_RESET,NPAR,RMS,NQ
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  LOGICAL :: DONE, STEP
  INTEGER :: ALPHA,BETA,AB(2),NQ0,NQTOT,NDUDS,NBLOCKS,L0(NATOMS)
  DOUBLE PRECISION :: POTEL, E0
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.  
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  IF(LFLIPS_RESET) THEN
     WRITE(MYUNIT,'(A)') 'QALCSearch> Resetting stoichiometry...'
     CALL RESET_STOICHIOMETRY()
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF (NPAR.GT.1) THEN
        WRITE(MYUNIT,'(A,I1,A,I10,A,F20.10,A,I5,A,G12.5,A,F11.1)') &
             '[',NP,']Qu ',NQ(NP),' E=',POTEL,' steps=',ITER, &
             ' RMS=',RMS,' t=',TIME
     ELSE
        WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,F11.1)') &
             'Qu ',NQ(NP),' E=',POTEL,' steps=',ITER, &
             ' RMS=',RMS,' t=',TIME
     ENDIF
  ENDIF
  !
  CALL MYCPU_TIME(TIME)
  WRITE(MYUNIT,'(A,F20.10,A,F11.1)') &
       'QALCS> Initial E= ',POTEL,' t= ',TIME-TSTART
  !
  ! Count species with non-zero population
  BETA=0
  DO ALPHA=1,NSPECIES(0)
     IF(NSPECIES(ALPHA) > 0) BETA = BETA+1
  ENDDO
  NBLOCKS = BETA*(BETA-1)/2
  !
  DONE = .FALSE.
  NQ0 = NQTOT
  !
  ! The outer while loop is for QALCS_SYM or QALCS_SURF.
  DO WHILE(.NOT. DONE) 
     !
     DONE = .TRUE.
     !
     L0(1:NATOMS) = INVATOMLISTS(1:NATOMS,1)
     SEQLENGTH = 0
     !
     IF(QALCST) THEN ! Loop over inter-species swap types.
        !
        IF(QALCSMODE == 0) THEN ! Steepest descent (slow)
           !
           STEP = .TRUE.
           DO WHILE(STEP)
              CALL SPAN_SWAPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
           ENDDO
           !
        ELSEIF(QALCSMODE >= 1 .AND. QALCSMODE < 4) THEN 
           !
           ! More efficient descents with each interspecies
           ! swap types are treated separately. There are NBLOCKS
           ! different swap types, and we converge when all of 
           ! them fail to yield an improvement.
           !
           NDUDS = 0
           DO WHILE(NDUDS < NBLOCKS)
              NDUDS = 1
              DO ALPHA=1,NSPECIES(0)-1
                 IF(NSPECIES(ALPHA) == 0) CYCLE
                 AB(1) = ALPHA
                 DO BETA=ALPHA+1,NSPECIES(0)
                    IF(NSPECIES(BETA) == 0) CYCLE
                    AB(2) = BETA
                    E0 = POTEL
                    CALL QALCS_AB(NP,ITER,TIME,BRUN,QDONE,SCREENC,AB)
                    IF(POTEL < E0 - ECONV) THEN
                       NDUDS = 1 ! reset the dud streak
                    ELSE ! Increment the count for dud blocks.
                       NDUDS = NDUDS + 1
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
           !
        ELSEIF(QALCSMODE >= 4 .AND. QALCSMODE < 6) THEN
           !
           ! In modes 4 and 5 all interspecies swap types
           ! are lumped into a single neighbourhood, which 
           ! is then ranked by unquenched swap-gain.
           STEP=.TRUE.
           DO WHILE(STEP)
              CALL SCAN_SWAP_NBRHD(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
           ENDDO
           !
        ELSEIF(QALCSMODE == 6) THEN
           !
           ! Steepest descent FLIP sequence..
           STEP = .TRUE.
           DO WHILE(STEP)
              CALL SPAN_FLIPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
           ENDDO
           !
        ELSEIF(QALCSMODE.EQ.7 .OR. QALCSMODE.EQ.8) THEN
           !
           ! In modes 7 and 8 all interspecies flips
           ! are lumped into a single neighbourhood, which 
           ! is then ranked by unquenched swap-gain (8) or 
           ! scanned in random order (7).
           !
           STEP=.TRUE.
           DO WHILE(STEP)
              CALL SCAN_FLIP_NBRHD(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
           ENDDO
           !
        ELSE
           !
           WRITE(MYUNIT,'(A)') 'QALCSearch> Bad QALCSMODE!'
           STOP
           !
        ENDIF
        !
     ENDIF
     !
     ! Save multminimum
     IF(SAVEMULTIMINONLY) THEN
        CALL GSAVEIT_MC(POTEL,COORDS(1:3*NATOMS,NP), &
             INVATOMLISTS(1:NATOMS,1), NP)
     ENDIF
     !
     !ds656> Testing
     ! WRITE(MYUNIT,'(A)') 'QALCS> Final permutation:'
     ! BETA=0
     ! DO ALPHA=1,NATOMS
     !    WRITE(MYUNIT,'(I3)',ADVANCE='NO') INVATOMLISTS(ALPHA,1)
     !    IF(INVATOMLISTS(ALPHA,1) == 1) BETA=BETA+1
     ! ENDDO
     ! WRITE(MYUNIT,*)
     ! IF(BETA /= NSPECIES(1)) THEN
     !    WRITE(MYUNIT,'(A,I4)') 'QALCS_ab> Bad count for type-1:',BETA
     ! ENDIF
     !<ds656 End of testing.
     !
     CALL CALC_HAMMING_DISTANCE(NATOMS,L0(1:NATOMS),&
          INVATOMLISTS(1:NATOMS,1),ALPHA)
     !
     CALL MYCPU_TIME(TIME)
     WRITE(MYUNIT,'(A,F20.10,A,I9,A,F11.1)') &
          'QALCS> Biminimum E= ',POTEL, &
          ' after ',NQTOT-NQ0,' quenches t= ',TIME-TSTART
     WRITE(MYUNIT,'(A,I5,A,I5)') &
          'QALCS> Swap-sequence length: ',SEQLENGTH,&
          ' Hamming distance: ',ALPHA
          
     !
     E0 = POTEL
     IF(QALCS_SYMT) THEN ! Symmetrisation scheme
        CALL QALCS_SYM(NP,ITER,TIME,BRUN,QDONE,SCREENC)
     ENDIF
     IF(QALCS_SURFT) THEN ! Perform DLS-like surface optimisation.
        CALL QALCS_SURF(NP,ITER,TIME,BRUN,QDONE,SCREENC)
     ENDIF
     IF( POTEL < E0 - ECONV .AND. QALCST) DONE = .FALSE.
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE QALCS
!
!=============================================================
!
SUBROUTINE QALCS_AB(NP,ITER,TIME,BRUN,QDONE,SCREENC,AB)
  !
  ! QALCS for a particular AB-type swap
  !
  USE COMMONS, ONLY : NATOMS, MYUNIT, ATOMLISTS, INVATOMLISTS, &
       NSPECIES, QALCSV, QALCSMODE
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP,AB(2)
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH      
  !
  LOGICAL :: UNCONVERGED, SCAN_ALL
  INTEGER :: LISTS_AB(1:2,0:NATOMS),I_AB(1:2),I,J
  DOUBLE PRECISION :: FLIPS_AB(1:2,1:NATOMS)
  !
  WRITE(MYUNIT,'(A,I2,A,I2)') &
       'QALCS_ab> Treating species ',AB(1),' and ',AB(2)
  !
  LISTS_AB(1,0:NATOMS) = ATOMLISTS(AB(1),1,0:NATOMS)
  LISTS_AB(2,0:NATOMS) = ATOMLISTS(AB(2),1,0:NATOMS)
  !
  UNCONVERGED = .TRUE.
  SCAN_ALL = .TRUE.
  !
  DO WHILE(UNCONVERGED)
     !
     !ds656> Testing...
     IF(QALCSV) THEN
        DO I=1,2
           WRITE(MYUNIT,'(A,I2,A)',ADVANCE='NO') &
                'QALCS_ab> Current auxiliary list for type ', &
                AB(I),':'
           DO J=1,LISTS_AB(I,0)
              WRITE(MYUNIT,'(1X,I4)',ADVANCE='NO') LISTS_AB(I,J)
           ENDDO
           WRITE(MYUNIT,'(A)')'.'
        ENDDO
     ENDIF
     !<ds656 ...testing
     !
     IF(QALCSMODE == 1) THEN
        CALL CALC_FLIPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,AB, &
             LISTS_AB,FLIPS_AB)
        CALL AB_SEARCH(NP,ITER,TIME,BRUN,QDONE,SCREENC,&
             LISTS_AB,FLIPS_AB,I_AB,SCAN_ALL)
     ELSEIF(QALCSMODE == 2 .OR. QALCSMODE == 3) THEN
        CALL AB_SEARCH2(NP,ITER,TIME,BRUN,QDONE,SCREENC,&
             LISTS_AB,I_AB)
     ELSE
        WRITE(MYUNIT,'(A)') 'QALCS_ab> Bad QALCSMODE!'
        STOP
     ENDIF
     !
     IF(I_AB(1)==0 .AND. I_AB(2)==0) THEN ! no good swaps found
        !
        !WRITE(MYUNIT,'(A)') 'QALCS_AB> No good swaps candidates.'
        !
        IF( SCAN_ALL ) THEN
           UNCONVERGED = .FALSE.
        ELSE
           ! Reset auxiliary lists
           LISTS_AB(1,0:NATOMS) = ATOMLISTS(AB(1),1,0:NATOMS)
           LISTS_AB(2,0:NATOMS) = ATOMLISTS(AB(2),1,0:NATOMS)
           SCAN_ALL = .TRUE.
           !
           IF(QALCSV) WRITE(MYUNIT,'(A)') &
                'QALCS_ab> Auxiliary list reset.'
           !
        ENDIF
        !
     ELSEIF(QALCSMODE == 1) THEN ! decrement auxiliary lists
        !
        IF(QALCSV) WRITE(MYUNIT,'(A,2(1X,I3),A)') &
             'QALCS_ab> Removed atoms', &
             I_AB(1),I_AB(2),' from auxiliary list.'
        !
        DO I=1,2 ! span types A and B
           DO J=1,LISTS_AB(I,0) ! scan each list
              IF(LISTS_AB(I,J)==I_AB(1) .OR. &
                   LISTS_AB(I,J)==I_AB(2)) THEN
                 LISTS_AB(I,J) = LISTS_AB(I,LISTS_AB(I,0))
                 LISTS_AB(I,LISTS_AB(I,0)) = 0
                 LISTS_AB(I,0) = LISTS_AB(I,0) - 1
              ENDIF
           ENDDO
        ENDDO
        !
        SCAN_ALL = .FALSE.
        !
     ELSE
        ! Update lists for QALCSMODE 2-5
        LISTS_AB(1,0:NATOMS) = ATOMLISTS(AB(1),1,0:NATOMS)
        LISTS_AB(2,0:NATOMS) = ATOMLISTS(AB(2),1,0:NATOMS)
     ENDIF
     !
     !ds656> Testing...
     IF(QALCSV) THEN
        WRITE(MYUNIT,'(A)',ADVANCE='NO') &
             'QALCS_ab> Current label assignment:'
        J=0
        DO I=1,NATOMS
           WRITE(MYUNIT,'(I3)',ADVANCE='NO') INVATOMLISTS(I,1)
           IF(INVATOMLISTS(I,1) == 1) J=J+1
        ENDDO
        WRITE(MYUNIT,*)
        IF(J /= NSPECIES(1)) THEN
           WRITE(MYUNIT,'(A,I4)') &
                'QALCS_ab> Bad count for type-1:',J
        ENDIF
     ENDIF
     !ds656> ...testing
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE QALCS_AB
!
!=============================================================
!
SUBROUTINE AB_SEARCH2(NP,ITER,TIME,BRUN,QDONE,SCREENC,LISTS_AB,I_AB)
  !
  ! Randomly search for negative swap gain in AB-block.
  !
  USE PORFUNCS
  USE COMMONS, ONLY : NATOMS, COORDS, QALCSV, NQ, ECONV, &
       MYUNIT, QALCSMODE, SEQLENGTH
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP,LISTS_AB(1:2,0:NATOMS)
  INTEGER, INTENT(OUT) :: I_AB(1:2)
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  INTEGER :: NQTOT, N, NTOT, SWAPS(LISTS_AB(1,0)*LISTS_AB(2,0),2),&
       I, J, K
  DOUBLE PRECISION :: POTEL, E0, X0(3*NATOMS), DPRAND, &
       EST(LISTS_AB(1,0)*LISTS_AB(2,0))
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  !=== Build (un)sorted list of non-degenerate swaps ============
  !
  ! Here X0 and I_AB are used as dummy variables!
  !
  N = 0
  DO I=1, LISTS_AB(1,0)
     DO J=1, LISTS_AB(2,0)
        !
        ! Append to list of swaps only if swap is non-degenerate
        !
        N = N + 1
        SWAPS(N,1) = LISTS_AB(1,I)
        SWAPS(N,2) = LISTS_AB(2,J)
        !
        IF(QALCSMODE==3) THEN ! Sort list
           ! Perform swap and compute unquenched energy.
           CALL SWAP_LABELS(LISTS_AB(1,I),LISTS_AB(2,J))
           CALL POTENTIAL(COORDS(:,NP),X0(:),EST(N),.FALSE.,.FALSE.)
           CALL SWAP_LABELS(LISTS_AB(1,I),LISTS_AB(2,J))
           ! Sort by unquenched energy
           SORT_EST: DO K=N,2,-1
              IF(EST(K) < EST(K-1)) THEN
                 E0 = EST(K); EST(K) = EST(K-1); EST(K-1) = E0
                 I_AB(1:2) = SWAPS(K,1:2)
                 SWAPS(K,1:2) = SWAPS(K-1,1:2)
                 SWAPS(K-1,1:2) = I_AB(1:2)
              ELSE
                 EXIT SORT_EST
              ENDIF
           ENDDO SORT_EST
        ENDIF
        !
     ENDDO
  ENDDO
  !==============================================================
  !
  ! Now E0 and I_AB are no longer dummies!
  !
  E0 = POTEL
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  I_AB(1:2) = 0
  !
  NTOT = N
  DO I=1,NTOT
     !
     IF(QALCSMODE==2) THEN    
        J = INT(DPRAND()*DBLE(N)) + 1
     ELSEIF(QALCSMODE==3) THEN
        J=I
     ENDIF
     !
     CALL SWAP_LABELS(SWAPS(J,1),SWAPS(J,2))
     !                                                            
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(QALCSV) THEN
        CALL PRINT_QUENCH(NP, ITER, '  ')
        IF(QALCSMODE==3) THEN
           WRITE(MYUNIT,'(2(A,F15.10))') &
                'ab_search2> dE*= ',EST(J)-E0,' dE= ',POTEL-E0 
        ENDIF
     ENDIF
     !
     IF(POTEL < E0 - ECONV) THEN
        SEQLENGTH = SEQLENGTH + 1
        WRITE(MYUNIT,'(A,I4,A,I4,A,F15.10)') &
             'ab_search2> ',SWAPS(J,1),' <-> ',SWAPS(J,2),&
             ' => dE= ', POTEL-E0
        I_AB(1:2) = SWAPS(J,1:2)
        RETURN
     ELSE 
        ! Undo the swap and restore the configuratiion
        CALL SWAP_LABELS(SWAPS(J,1),SWAPS(J,2))
        COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
        POTEL = E0
        !
        IF(QALCSMODE==2) THEN
           ! Remove the rejected swap from list
           SWAPS(J,1:2) = SWAPS(N,1:2)
           SWAPS(N,1:2) = 0
           N = N - 1
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  WRITE(MYUNIT,'(A)') 'ab_search2> No swaps with -ve gain!'
  !
  RETURN
  !
END SUBROUTINE AB_SEARCH2
!
!=============================================================
!
SUBROUTINE AB_SEARCH(NP,ITER,TIME,BRUN,QDONE,SCREENC, & 
     LISTS_AB,FLIPS_AB,I_AB,SCAN_ALL)
  !
  ! Search for negative swap gain in AB-(sub)block.
  !
  USE PORFUNCS
  USE COMMONS, ONLY : NATOMS, COORDS, QALCSV, NQ, ECONV, &
       MYUNIT, SEQLENGTH
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  LOGICAL, INTENT(IN) :: SCAN_ALL
  INTEGER, INTENT(IN) :: NP,LISTS_AB(1:2,0:NATOMS)
  INTEGER, INTENT(OUT) :: I_AB(1:2)
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  DOUBLE PRECISION, INTENT(IN) :: FLIPS_AB(1:2,1:NATOMS)
  !
  LOGICAL :: COMPLETE, STUCK
  INTEGER :: SUBLISTS_AB(1:2,0:NATOMS), NQTOT, &
       IPICK1, IPICK2, TPICK1, TPICK2, IPICK1BEST, I, J
  DOUBLE PRECISION :: X0(3*NATOMS),XBEST(3*NATOMS),E0,EBEST
  DOUBLE PRECISION :: SUBFLIPS_AB(1:2, 1:NATOMS), DUMMY, POTEL
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  ! Sublfips and sublists must be shuffled syncronously!
  SUBLISTS_AB(1:2, 0:NATOMS) = LISTS_AB(1:2, 0:NATOMS)
  SUBFLIPS_AB(1:2, 1:NATOMS) = FLIPS_AB(1:2, 1:NATOMS)
  !
  I_AB(1:2) = 0
  COMPLETE = .FALSE.
  !
  list: DO WHILE(.NOT. COMPLETE) ! loop over swaps
     !
     !WRITE(MYUNIT,'(A)') 'ab_search> (re)started alternating search'
     !
     X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
     XBEST(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
     E0=POTEL
     EBEST=POTEL
     !
     ! Find atom index with lowest flip gain
     TPICK1 = 1; IPICK1 = 1
     DUMMY = SUBFLIPS_AB(1,1)
     DO J=1,2
        DO I=1,SUBLISTS_AB(J,0)
           IF(SUBFLIPS_AB(J,I) < DUMMY) THEN
              DUMMY = SUBFLIPS_AB(J,I)
              TPICK1 = J
              IPICK1 = I
           ENDIF
        ENDDO
     ENDDO
     !
     WRITE(MYUNIT,'(A,I3)') 'ab_search> Best flipper: ', &
          SUBLISTS_AB(TPICK1, IPICK1)
     !
     STUCK = .FALSE.
     !
     sublist: DO WHILE(.NOT. STUCK)
        !
        STUCK=.TRUE.
        !
        IPICK2 = IPICK1
        TPICK2 = TPICK1
        !
        IF(TPICK1 == 1) THEN
           TPICK1 = 2
        ELSEIF(TPICK1 == 2) THEN
           TPICK1 = 1 
        ELSE
           STOP
        ENDIF
        !
        IPICK1BEST = 0
        !
        !ds656> Testing...
        IF(QALCSV) THEN
           WRITE(MYUNIT,'(A, I3, A)',ADVANCE='NO') &
                'ab_search> Scanning pairs for atom ', &
                SUBLISTS_AB(TPICK2,IPICK2),':'
           DO IPICK1=1,SUBLISTS_AB(TPICK1,0)
              WRITE(MYUNIT,'(1X,I3)',ADVANCE='NO') &
                   SUBLISTS_AB(TPICK1,IPICK1)
           ENDDO
           WRITE(MYUNIT,'(A)') '.'
        ENDIF
        !<ds656 ...testing
        !
        scan: DO IPICK1=1,SUBLISTS_AB(TPICK1,0)
           ! 
           !WRITE(MYUNIT,*) 'AB_SWAP_SEARCH> TPICK1,IPICK1,I,:', &
           !     TPICK1, IPICK1, SUBLISTS_AB(TPICK1, IPICK1)
           !WRITE(MYUNIT,*) 'AB_SWAP_SEARCH> TPICK2,IPICK2,J,:', &
           !     TPICK2, IPICK2, SUBLISTS_AB(TPICK2, IPICK2)
           !
           CALL SWAP_LABELS( SUBLISTS_AB(TPICK1,IPICK1), &
                SUBLISTS_AB(TPICK2,IPICK2) )
           !
           NQTOT = NQTOT + 1
           NQ(NP) = NQ(NP) + 1
           CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
           IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '  ')
           !
           IF(POTEL < EBEST - ECONV) THEN
              COMPLETE = .TRUE.
              STUCK = .FALSE.
              ! Store best-encountered swap candidates
              I_AB(1) = SUBLISTS_AB( TPICK1, IPICK1)
              I_AB(2) = SUBLISTS_AB( TPICK2, IPICK2)
              IPICK1BEST = IPICK1
              EBEST = POTEL
              XBEST(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
              WRITE(MYUNIT,'(A,I4,A,I4,A,F15.10)') &
                   'ab_search> ', I_AB(1),' <-> ', I_AB(2), &
                   ' => dE= ',EBEST-E0
              CALL FLUSH(MYUNIT)
           ENDIF
           !
           ! Repeat the swap (to undo it) and reset the coordinates
           CALL SWAP_LABELS( SUBLISTS_AB(TPICK1,IPICK1), &
                SUBLISTS_AB(TPICK2,IPICK2) )
           !SCREENC(1:3*NATOMS) = X0(1:3*NATOMS) !?
           COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
           POTEL = E0
           !
        ENDDO scan
        !WRITE(MYUNIT,'(A)') 'ab_search> ...finished row/column!'
        !
        IPICK1 = IPICK1BEST
        !
        ! Remove TPICK2, IPICK2 from SUBLISTS_AB and SUBFLIPS_AB
        SUBLISTS_AB(TPICK2,IPICK2) = &
             SUBLISTS_AB(TPICK2,SUBLISTS_AB(TPICK2,0))
        SUBLISTS_AB(TPICK2,SUBLISTS_AB(TPICK2,0)) = 0
        SUBFLIPS_AB(TPICK2,IPICK2) = &
             SUBFLIPS_AB(TPICK2,SUBLISTS_AB(TPICK2,0))
        SUBFLIPS_AB(TPICK2,SUBLISTS_AB(TPICK2,0)) = 0
        SUBLISTS_AB(TPICK2,0) = SUBLISTS_AB(TPICK2,0) - 1
        !
        IF(SUBLISTS_AB(TPICK2,0) == 0) THEN ! list depleted!
           STUCK = .TRUE. ! Exit loop 'sublist'
           COMPLETE = .TRUE. ! Exit loop 'list'
           WRITE(MYUNIT,'(A)') 'ab_search> List depleted!'
        ENDIF
        !
     ENDDO sublist
     !
     IF(.NOT. SCAN_ALL) COMPLETE = .TRUE.
     !
  ENDDO list
  !
  IF(EBEST < E0 - ECONV) THEN
     SEQLENGTH = SEQLENGTH + 1
     CALL SWAP_LABELS( I_AB(1), I_AB(2) )
     !SCREENC(1:3*NATOMS) = XBEST(1:3*NATOMS)
     COORDS(1:3*NATOMS, NP) = XBEST(1:3*NATOMS)
     POTEL = EBEST
     WRITE(MYUNIT,'(A,I4,A,I4,A,F20.10)') &
          'ab_search> Swapped ',I_AB(1),' and ',I_AB(2), &
          ' to get E=',POTEL
     CALL FLUSH(MYUNIT)
  !ELSE
  !   WRITE(MYUNIT,'(A)') 'ab_search> Found no good swaps.'
  ENDIF
  !
  RETURN
  !
END SUBROUTINE AB_SEARCH
!
!=============================================================
!
SUBROUTINE CALC_FLIPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,AB, &
     LISTS_AB, FLIPS_AB)
  !
  ! USE COMMONS, ONLY : NATOMS, COORDS, NQ, QALCSV, MYUNIT
  USE COMMONS, ONLY : NATOMS, COORDS, NQ, QALCSV, SAVEMULTIMINONLY
  !
  IMPLICIT NONE 
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP,LISTS_AB(1:2,0:NATOMS),AB(1:2)
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  DOUBLE PRECISION, INTENT(OUT) :: FLIPS_AB(1:2,1:NATOMS)
  !
  LOGICAL :: TEMPSAVE
  INTEGER :: BA(2), T, I, J, NQTOT
  DOUBLE PRECISION :: X0(1:3*NATOMS), E0, POTEL
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  BA(1) = AB(2)
  BA(2) = AB(1)
  !
  FLIPS_AB(1:2,1:NATOMS) = 0.0D0
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  E0 = POTEL
  !
  DO T=1,2
     DO I=1,LISTS_AB(T,0)
        !
        ! Flip the label of current atom index and quench
        J=LISTS_AB(T,I) ! permanent atom index
        CALL FLIP_LABEL(J,BA(T))
        !write(MYUNIT,'(A,4(1X,I2))') 'calc_flips> T,I,J,BA(T)=',T,I,J,BA(T)
        !
        TEMPSAVE = SAVEMULTIMINONLY
        SAVEMULTIMINONLY = .TRUE. ! never save minima with 1 atom flipped!
        NQTOT = NQTOT + 1
        NQ(NP) = NQ(NP) + 1
        CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
        IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '* ')
        SAVEMULTIMINONLY = TEMPSAVE
        !
        ! Store the new energy after quenched flip
        FLIPS_AB(T,I) = POTEL
        !
        ! Undo the flip and revert the configuration
        CALL FLIP_LABEL(J,AB(T))
        COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
        POTEL = E0
        !
     ENDDO
  ENDDO
  !
  !STOP ! TESTING!
  !
  RETURN
  !
END SUBROUTINE CALC_FLIPS
!
!=============================================================
!
SUBROUTINE SWAP_LABELS(I,J)
  !
  ! Swap the labels/types of two particles identified by their
  ! permanent indices I and J.
  !
  ! USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, MYUNIT, ATMASS, READMASST
  USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, ATMASS, &
       READMASST, SPECMASST
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I,J
  INTEGER :: ITYPE,JTYPE,IPLACE,JPLACE,IGROUP,JGROUP,K
  DOUBLE PRECISION :: DUMMY
  !
  !WRITE(MYUNIT,*) 'SWAPPING I,J=', I,J
  !
  ITYPE = INVATOMLISTS(I,1)
  IGROUP = INVATOMLISTS(I,2)
  IPLACE = INVATOMLISTS(I,3)
  !
  JTYPE = INVATOMLISTS(J,1)
  JGROUP = INVATOMLISTS(J,2)
  JPLACE = INVATOMLISTS(J,3)
  !
  !WRITE(MYUNIT,*) 'ITYPE, IGROUP, IPLACE=', ITYPE, IGROUP, IPLACE
  !WRITE(MYUNIT,*) 'JTYPE, JGROUP, JPLACE=', JTYPE, JGROUP, JPLACE
  !
  !WRITE(MYUNIT,*) 'INITIAL ITYPE ( ', ITYPE, ' ) LIST:', &
  !     (ATOMLISTS(ITYPE,IGROUP,K), K=1,ATOMLISTS(ITYPE,IGROUP,0))
  !WRITE(MYUNIT,*) 'INITIAL JTYPE ( ', JTYPE, ' ) LIST:', &
  !     (ATOMLISTS(JTYPE,JGROUP,K), K=1,ATOMLISTS(JTYPE,JGROUP,0))
  !
  ! Swap labels
  ATOMLISTS(ITYPE, IGROUP, IPLACE) = J
  ATOMLISTS(JTYPE, JGROUP, JPLACE) = I
  !
  IF(READMASST.OR.SPECMASST) THEN   ! Swap masses
     DUMMY = ATMASS(I)
     ATMASS(I) = ATMASS(J)
     ATMASS(J) = DUMMY
  ENDIF
  !
  INVATOMLISTS(I,1) = JTYPE
  INVATOMLISTS(I,2) = JGROUP
  INVATOMLISTS(I,3) = JPLACE
  !
  INVATOMLISTS(J,1) = ITYPE
  INVATOMLISTS(J,2) = IGROUP
  INVATOMLISTS(J,3) = IPLACE
  !
  !WRITE(MYUNIT,*) 'FINAL ITYPE (', ITYPE, ') LIST:', &
  !     (ATOMLISTS(ITYPE,IGROUP,K), K=1,ATOMLISTS(ITYPE,IGROUP,0))
  !WRITE(MYUNIT,*) 'FINLA JTYPE (', JTYPE, ') LIST:', &
  !     (ATOMLISTS(JTYPE,JGROUP,K), K=1,ATOMLISTS(JTYPE,JGROUP,0))
  RETURN
  !
END SUBROUTINE SWAP_LABELS
!
!=============================================================
!
SUBROUTINE FLIP_LABEL(I,NEWTYPE)
  !
  ! Change the current label/type of a particle (identified by its
  ! permanent index I) to NEWTYPE.
  !
  ! USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, MYUNIT
  USE COMMONS, ONLY: ATOMLISTS, INVATOMLISTS, NSPECIES, ATMASS, &
       READMASST, SPECMASST, SPECMASS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I, NEWTYPE
  INTEGER :: ITYPE,IPLACE,IGROUP,J
  !
  ITYPE = INVATOMLISTS(I,1)
  IGROUP = INVATOMLISTS(I,2)
  IPLACE = INVATOMLISTS(I,3)
  !
  IF(IPLACE < ATOMLISTS(ITYPE,IGROUP,0)) THEN
     ! Fill the hole at IPLACE with the last entry for ITYPE
     J = ATOMLISTS(ITYPE,IGROUP,ATOMLISTS(ITYPE,IGROUP,0))
     ATOMLISTS(ITYPE,IGROUP,IPLACE) = J
     INVATOMLISTS(J,3) = IPLACE       
  ENDIF
  ATOMLISTS(ITYPE,IGROUP,ATOMLISTS(ITYPE,IGROUP,0)) = 0
  ATOMLISTS(ITYPE,IGROUP,0) = ATOMLISTS(ITYPE,IGROUP,0) - 1
  ! Is there a problem when IPLACE is the last entry, i.e. J = I?
  !
  ! Now grow append I as the last entry for NEWTYPE
  ATOMLISTS(NEWTYPE,IGROUP,0) = ATOMLISTS(NEWTYPE,IGROUP,0) + 1
  ATOMLISTS(NEWTYPE,IGROUP,ATOMLISTS(NEWTYPE,IGROUP,0)) = I
  INVATOMLISTS(I,1) = NEWTYPE
  INVATOMLISTS(I,3) = ATOMLISTS(NEWTYPE,IGROUP,0)
  !
  NSPECIES(ITYPE) = NSPECIES(ITYPE) - 1
  NSPECIES(NEWTYPE) = NSPECIES(NEWTYPE) + 1
  !
  IF(READMASST.OR.SPECMASST) ATMASS(I) = SPECMASS(NEWTYPE)
  !
  !WRITE(MYUNIT,*) 'TYPE', ITYPE, 'LIST:', &
  !     ( ATOMLISTS(ITYPE,IGROUP,J), J=1,ATOMLISTS(ITYPE,IGROUP,0) )
  !WRITE(MYUNIT,*) 'TYPE', NEWTYPE, 'LIST:', &
  !     ( ATOMLISTS(NEWTYPE,IGROUP,J), J=1,ATOMLISTS(NEWTYPE,IGROUP,0) )
  
  !
END SUBROUTINE FLIP_LABEL
!
!=============================================================
!
SUBROUTINE PRINT_QUENCH(NP, ITER, TAG)
  !
  ! Standard post-quench print statement.
  !
  USE COMMONS, ONLY : RMS, NPAR, NQ, MYUNIT
  !
  INTEGER, INTENT(IN) :: NP, ITER
  CHARACTER(LEN=2), INTENT(IN) :: TAG
  DOUBLE PRECISION :: POTEL
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  !
  IF (NPAR.GT.1) THEN
     WRITE(MYUNIT,'(A,I1,A,A,I10,A,F20.10,A,I5,A,G12.5,A)') &
          '[',NP,']Qu', TAG, NQ(NP),' E=',POTEL,' steps=',ITER, &
          ' RMS=',RMS
  ELSE
     WRITE(MYUNIT,'(A,A,I10,A,F20.10,A,I5,A,G12.5)') &
          'Qu', TAG, NQ(NP),' E=',POTEL,' steps=',ITER, &
          ' RMS=',RMS
  ENDIF
  !
END SUBROUTINE PRINT_QUENCH
!
!=============================================================
!
SUBROUTINE SPAN_SWAPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
  !
  ! Calculated all swap gains and count the number of 
  ! negative ones. If STEP=.true. on input, then find
  ! and execute the lowest flip gain; and in the absence
  ! of negative flip gains set STEP=.false on output.
  ! Determine if the configuration is a minimum or a
  ! saddle (or neither) in permutation space.
  !
  USE COMMONS, ONLY : NSPECIES,ATOMLISTS,NATOMS,COORDS,NQ, &
       QALCSV,ECONV,MYUNIT,SEQLENGTH
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
  ! 
  LOGICAL, INTENT(INOUT) :: STEP
  !
  INTEGER :: TA, TB, GA, GB, IA, IB, I, J, NTOT, NNEG, NQTOT, &
       BESTSWAP(2,2)
  DOUBLE PRECISION :: POTEL, E0, X0(3*NATOMS), E0LO, &
       EMIN(2), XMIN(3*NATOMS)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  E0 = POTEL
  E0LO = E0 - ECONV
  EMIN(1:2) = POTEL
  BESTSWAP(1:2,1:2) = 0
  !
  NNEG=0 ! counter for -ve swap gains
  !NTOT=0
  !
  DO TA=1,NSPECIES(0)-1 ! Double loop over species
     DO TB=TA+1,NSPECIES(0)
        !
        DO GA = 1,3 ! Double loop over all groups
           DO GB = 1,3 
              !
              DO IA=1,ATOMLISTS(TA,GA,0) ! Double loop over atoms
                 I=ATOMLISTS(TA,GA,IA)
                 DO IB=1,ATOMLISTS(TB,GB,0)
                    J=ATOMLISTS(TB,GB,IB)
                    !
                    CALL SWAP_LABELS(I,J)
                    !
                    NQTOT = NQTOT + 1
                    NQ(NP) = NQ(NP) + 1
                    CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN, &
                         QDONE,SCREENC)
                    IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '  ')
                    !
                    !NTOT = NTOT+1
                    IF(POTEL < E0LO) THEN
                       NNEG=NNEG+1                       
                       IF(POTEL < EMIN(2) - ECONV) THEN
                          EMIN(2) = POTEL
                          BESTSWAP(2,1:2) = (/I,J/)
                          IF(POTEL < EMIN(1) - ECONV) THEN
                             EMIN(2) = EMIN(1)
                             BESTSWAP(2,1:2) = BESTSWAP(1,1:2)
                             EMIN(1) = POTEL
                             BESTSWAP(1,1:2) = (/I,J/)
                             IF(STEP) XMIN(:) = COORDS(:,NP)
                          ENDIF
                       ENDIF
                    ENDIF
                    !
                    CALL SWAP_LABELS(I,J)
                    POTEL = E0
                    COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
                    !
                 ENDDO
              ENDDO
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  WRITE(MYUNIT, '(A, I6, A)', ADVANCE='NO') &
       'span_swaps> ',NNEG,' -ve swap-gain(s)'
  !
  IF(NNEG == 0) THEN ! permutational (bi)minimum
     WRITE(MYUNIT, '(A)') ' => local optimum.'
  ELSEIF(NNEG == 2) THEN
     WRITE(MYUNIT,'(A,I3,A,I3,A,I3,A,I3)', ADVANCE='NO') ': ', &
          BESTSWAP(1,1),' <-> ', BESTSWAP(1,2), ' and ', &
          BESTSWAP(2,1),' <-> ', BESTSWAP(2,2)
     IF(  BESTSWAP(1,1) /= BESTSWAP(2,1) .AND. &
          BESTSWAP(1,2) /= BESTSWAP(2,2) ) THEN ! 'permutational saddle'
        WRITE(MYUNIT, '(A)') ' => saddle.'
     ELSE
        WRITE(MYUNIT, '(A)') ' => ?'
     ENDIF
  ELSE
     WRITE(MYUNIT,*)
  ENDIF
  !
  IF(STEP) THEN
     IF(NNEG > 0) THEN
        SEQLENGTH = SEQLENGTH + 1
        POTEL = EMIN(1)
        COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
        CALL SWAP_LABELS(BESTSWAP(1,1),BESTSWAP(1,2))
        WRITE(MYUNIT, '(A, G20.10)') &
             'span_swaps> Executed best swap with dE=', &
             EMIN(1)-E0
     ELSE
        STEP=.FALSE.
     ENDIF
  ELSEIF(NNEG > 0) THEN
     STEP=.TRUE.
  ENDIF
  !
  RETURN
  !
END SUBROUTINE SPAN_SWAPS
!
SUBROUTINE SPAN_FLIPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
  !
  ! Calculated all flip gains and count the number of 
  ! negative ones. If STEP=.true. on input, then find
  ! and execute the lowest flip gain; and in the absence
  ! of negative flip gains set STEP=.false on output.
  ! Determine if the configuration is a minimum or a
  ! saddle (or neither) in permutation space.
  !
  USE COMMONS, ONLY : NSPECIES,ATOMLISTS,NATOMS,COORDS,NQ, &
       QALCSV,ECONV,MYUNIT,SEQLENGTH
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
  ! 
  LOGICAL, INTENT(INOUT) :: STEP
  !
  INTEGER :: TA, TB, GA, IA, I, NTOT, NNEG, NQTOT, &
       BESTFLIP(2,2)
  DOUBLE PRECISION :: POTEL, E0, X0(3*NATOMS), E0LO, &
       EMIN(2), XMIN(3*NATOMS)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  WRITE(MYUNIT,'(A)') &
       '============================================================'
  !
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  E0 = POTEL
  E0LO = E0 - ECONV
  EMIN(1:2) = POTEL
  BESTFLIP(1:2,1:2) = 0
  !
  NNEG=0 ! counter for -ve flip gains
  !NTOT=0
  !
  DO TA=1,NSPECIES(0) ! loop over species
     DO GA = 1,2 ! loop over groups
        DO IA=1,ATOMLISTS(TA,GA,0) ! loop over atoms
           I=ATOMLISTS(TA,GA,IA) ! actual atom index
           DO TB=1,NSPECIES(0)
              IF(TB.NE.TA) THEN ! attempt flip
                 !
                 CALL FLIP_LABEL(I,TB)
                 !
                 NQTOT = NQTOT + 1
                 NQ(NP) = NQ(NP) + 1
                 CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN, &
                      QDONE,SCREENC)
                 IF(QALCSV) CALL PRINT_QUENCH(NP, ITER, '  ')
                 !
                 !NTOT = NTOT+1
                 IF(POTEL < E0LO) THEN
                    NNEG=NNEG+1                       
                    IF(POTEL < EMIN(2) - ECONV) THEN
                       EMIN(2) = POTEL
                       BESTFLIP(2,1:2) = (/I,TB/)
                       IF(POTEL < EMIN(1) - ECONV) THEN
                          EMIN(2) = EMIN(1)
                          BESTFLIP(2,1:2) = BESTFLIP(1,1:2)
                          EMIN(1) = POTEL
                          BESTFLIP(1,1:2) = (/I,TB/)
                          IF(STEP) XMIN(:) = COORDS(:,NP)
                       ENDIF
                    ENDIF
                 ENDIF
                 !
                 CALL FLIP_LABEL(I,TA)
                 POTEL = E0
                 COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
                 !
              ENDIF
           ENDDO
        ENDDO
        !
     ENDDO
  ENDDO
  !
  WRITE(MYUNIT, '(A, I6, A)', ADVANCE='NO') &
       'span_flips> ',NNEG,' -ve flip-gain(s)'
  !
  IF(NNEG == 0) THEN ! permutational (bi)minimum
     WRITE(MYUNIT, '(A)') ' => local optimum.'
  !ELSEIF(NNEG == 2) THEN
  !   WRITE(MYUNIT,'(A,I3,A,I3,A,I3,A,I3)', ADVANCE='NO') ': ', &
  !        BESTSWAP(1,1),' <-> ', BESTSWAP(1,2), ' and ', &
  !        BESTSWAP(2,1),' <-> ', BESTSWAP(2,2)
  !   IF(  BESTSWAP(1,1) /= BESTSWAP(2,1) .AND. &
  !        BESTSWAP(1,2) /= BESTSWAP(2,2) ) THEN ! 'permutational saddle'
  !      WRITE(MYUNIT, '(A)') ' => saddle.'
  !   ELSE
  !      WRITE(MYUNIT, '(A)') ' => ?'
  !   ENDIF
  ELSE
     WRITE(MYUNIT,*)
  ENDIF
  !
  IF(STEP) THEN
     IF(NNEG > 0) THEN
        SEQLENGTH = SEQLENGTH + 1
        POTEL = EMIN(1)
        COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
        CALL FLIP_LABEL(BESTFLIP(1,1),BESTFLIP(1,2))
        WRITE(MYUNIT, '(A, G20.10)') &
             'span_flips> Executed best flip with dE=', &
             EMIN(1)-E0
     ELSE
        STEP=.FALSE.
     ENDIF
  ELSEIF(NNEG > 0) THEN
     STEP=.TRUE.
  ENDIF
  !
  RETURN
  !
END SUBROUTINE SPAN_FLIPS
!
!=============================================================
!
SUBROUTINE RANDMULTIPERM(NP, IGROUP)
  !
  ! ds656> 6/1/2015
  ! Randomly permute atomic labels in a multicomponent system.
  ! (More general than RANPERM.) Permute labels only for atoms
  ! in group IGROUP.
  !
  USE COMMONS, ONLY : NATOMS, NSPECIES, ATOMLISTS, INVATOMLISTS, &
       MYUNIT
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP, IGROUP
  INTEGER :: I,I1,I2,LABELS(NATOMS),LIST(NATOMS), N, L, NN
  DOUBLE PRECISION :: DPRAND, X
  !
  LABELS(:) = 0
  LIST(:) = 0
  N=0
  DO I1=1,NSPECIES(0)
     DO I2=1,ATOMLISTS(I1,IGROUP,0)
        N = N + 1
        LIST(N) = ATOMLISTS(I1,IGROUP,I2)
        LABELS(N) = I1
     ENDDO
     ATOMLISTS(I1,IGROUP,:) = 0 ! Reset global lists to zero.
  ENDDO
  !
  NN = N ! Since N will be decremented in the following loop.
  !
  DO I1=1,NN ! Loop over all atom indices in IGROUP.
     !
     ! Randomly choose a label for current atom index.
     I=LIST(I1) ! Actual atom index.
     I2=INT(DPRAND()*DBLE(N))+1 ! Randomly choose index for LABELS.
     L=LABELS(I2) ! Store chosen label.
     !
     ! Remove chosen label from contention by swapping with
     ! the Nth label in LIST and then decrementing N
     LABELS(I2) = LABELS(N) 
     LABELS(N) = L 
     N = N - 1 
     !
     ! Update global (ionverse) lists and for current atom.
     ATOMLISTS(L,IGROUP,0) = ATOMLISTS(L,IGROUP,0) + 1
     ATOMLISTS(L,IGROUP,ATOMLISTS(L,IGROUP,0)) = I
     INVATOMLISTS(I,1) = L
     INVATOMLISTS(I,2) = IGROUP
     INVATOMLISTS(I,3) = ATOMLISTS(L,IGROUP,0)
     !
  ENDDO
  !
  WRITE(MYUNIT,'(A, I2)') &
       'randmultiperm> Permuted labels of atoms in group ',IGROUP
  !
  RETURN
  !
END SUBROUTINE RANDMULTIPERM
!
!=============================================================
!
SUBROUTINE CALC_HAMMING_DISTANCE(N,LIST1,LIST2,NDIFF)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N, LIST1(N),LIST2(N)
  INTEGER, INTENT(OUT) :: NDIFF
  !
  INTEGER :: I
  !
  NDIFF=0
  DO I=1,N
     IF(LIST1(I).NE.LIST2(I)) NDIFF=NDIFF+1
  ENDDO
  !
  RETURN
  !
END SUBROUTINE CALC_HAMMING_DISTANCE
!
SUBROUTINE SCAN_SWAP_NBRHD(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
  !
  USE PORFUNCS
  USE COMMONS, ONLY : NATOMS, COORDS, QALCSV, NQ, ECONV, &
       MYUNIT, QALCSMODE, SEQLENGTH, ATOMLISTS, &
       NSPECIES, QALCS_NBRHD, QALCS_PARAM
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  LOGICAL, INTENT(OUT) :: STEP
  !
  INTEGER :: NQTOT, N, NTOT, SWAPS(QALCS_NBRHD,2), I, J, K, &
       TA, TB, GA, GB, IA, IB, I_AB(1:2)
  DOUBLE PRECISION :: POTEL, E0, X0(3*NATOMS), DPRAND, &
       EST(QALCS_NBRHD)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  !==============================================================
  ! Firs loop over the entire neighbourhood to compute the
  ! sorted list of pair-swaps.
  !
  N = 0
  DO TA=1,NSPECIES(0)-1 ! Double loop over species
     DO TB=TA+1,NSPECIES(0)
        !
        DO GA = 1,3 ! Double loop over all groups
           DO GB = 1,3 
              !
              DO IA=1,ATOMLISTS(TA,GA,0) ! Double loop over atoms
                 I=ATOMLISTS(TA,GA,IA)
                 DO IB=1,ATOMLISTS(TB,GB,0)
                    J=ATOMLISTS(TB,GB,IB)
                    !
                    N = N + 1
                    SWAPS(N,1) = I
                    SWAPS(N,2) = J
                    !
                    IF(QALCSMODE == 5) THEN
                       !
                       CALL SWAP_LABELS(I,J)
                       CALL POTENTIAL(COORDS(:,NP),X0(:),EST(N),.FALSE.,.FALSE.)
                       CALL SWAP_LABELS(I,J)
                       !
                       SORT_EST2: DO K=N,2,-1
                          IF(EST(K) < EST(K-1)) THEN
                             E0=EST(K); EST(K)=EST(K-1); EST(K-1)=E0
                             I_AB(1:2) = SWAPS(K,1:2)
                             SWAPS(K,1:2) = SWAPS(K-1,1:2)
                             SWAPS(K-1,1:2) = I_AB(1:2)
                          ELSE
                             EXIT SORT_EST2
                          ENDIF
                       ENDDO SORT_EST2
                       !
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
           ENDDO
        ENDDO
        !
     ENDDO
  ENDDO
  !
  IF(QALCS_NBRHD /= N) THEN
     WRITE(MYUNIT,'(A)') 'scan_nbrhd> Inconsistency!'
     STOP
  ENDIF
  !
  ! Now E0 and I_AB are no longer dummies!
  !
  E0 = POTEL
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  !
  IF(QALCS_PARAM > 0 .AND. QALCS_PARAM < N) THEN
     NTOT = QALCS_PARAM
  ELSE
     NTOT = N
  ENDIF
  !
  DO I=1,NTOT
     !
     IF(QALCSMODE==4) THEN    
        J = INT(DPRAND()*DBLE(N)) + 1
     ELSEIF(QALCSMODE==5) THEN
        J=I
     ENDIF
     !
     CALL SWAP_LABELS(SWAPS(J,1),SWAPS(J,2))
     !                                                            
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(QALCSV) THEN
        CALL PRINT_QUENCH(NP, ITER, '  ')
        IF(QALCSMODE==5) THEN
           WRITE(MYUNIT,'(2(A,F15.10))') &
                'scan_nbrhd> dE*= ',EST(J)-E0,' dE= ',POTEL-E0 
        ENDIF
     ENDIF
     !
     IF(POTEL < E0 - ECONV) THEN
        SEQLENGTH = SEQLENGTH + 1
        WRITE(MYUNIT,'(A,I4,A,I4,A,I4,A,F15.10)') &
             'scan_nbrhd> Try ',I,' : ',&
             SWAPS(J,1),' <-> ',SWAPS(J,2),&
             ' => dE= ', POTEL-E0
        STEP = .TRUE.
        RETURN
     ELSE 
        ! Undo the swap and restore the configuration
        CALL SWAP_LABELS(SWAPS(J,1),SWAPS(J,2))
        COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
        POTEL = E0
        !
        IF(QALCSMODE==4) THEN
           ! Remove the rejected swap from list
           SWAPS(J,1:2) = SWAPS(N,1:2)
           SWAPS(N,1:2) = 0
           N = N - 1
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  IF(NTOT == QALCS_PARAM) THEN
     WRITE(MYUNIT,'(A,I6,A)') &
          'scan_nbrhd> Failed to find improvement after ',&
          NTOT,' trial swaps!'
  ELSE
     WRITE(MYUNIT,'(A)') 'scan_nbrhd> No swaps with -ve gain!'
  ENDIF
  STEP = .FALSE.
  !
  RETURN
  !
END SUBROUTINE SCAN_SWAP_NBRHD
!
SUBROUTINE SCAN_FLIP_NBRHD(NP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
  !
  USE PORFUNCS
  USE COMMONS, ONLY : NATOMS, COORDS, QALCSV, NQ, ECONV, &
       MYUNIT,QALCSMODE,SEQLENGTH,ATOMLISTS,INVATOMLISTS,&
       NSPECIES,QALCS_NBRHD,QALCS_PARAM,SEMIGRAND_MUT,&
       SEMIGRAND_MU
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  LOGICAL, INTENT(OUT) :: STEP
  !
  INTEGER :: NQTOT, N, NTOT, FLIPS(NATOMS,2), I, J, K, &
       TA, TB, GA, IA, I_AB(1:2)
  DOUBLE PRECISION :: POTEL, E0, X0(3*NATOMS), DPRAND, &
       EST(QALCS_NBRHD), DUMMY
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  !==============================================================
  ! Firs loop over the entire neighbourhood to compute the
  ! sorted list of pair-swaps.
  !
  N = 0
  DO TA=1,NSPECIES(0) ! loop over species
     DO GA = 1,2 ! loop over groups 1 and 2
        DO IA=1,ATOMLISTS(TA,GA,0) ! loop over atoms
           I=ATOMLISTS(TA,GA,IA)
           DO TB=1,NSPECIES(0) ! loop over species
              IF(TB.NE.TA) THEN
                 N = N + 1
                 FLIPS(N,1) = I; FLIPS(N,2) = TB                    
                 !
                 IF(QALCSMODE == 8) THEN
                    !
                    CALL FLIP_LABEL(I,TB)
                    CALL POTENTIAL(COORDS(:,NP),X0(:),EST(N),.FALSE.,.FALSE.)
                    IF (SEMIGRAND_MUT) THEN
                       DUMMY=0.0D0
                       DO J=2, NSPECIES(0)
                          DUMMY = DUMMY + NSPECIES(J)*SEMIGRAND_MU(J)
                       ENDDO
                       EST(N) = EST(N) - DUMMY
                    ENDIF
                    CALL FLIP_LABEL(I,TA)
                    !
                    SORT_EST3: DO K=N,2,-1
                       IF(EST(K) < EST(K-1)) THEN
                          E0=EST(K); EST(K)=EST(K-1); EST(K-1)=E0
                          I_AB(1:2) = FLIPS(K,1:2)
                          FLIPS(K,1:2) = FLIPS(K-1,1:2)
                          FLIPS(K-1,1:2) = I_AB(1:2)
                       ELSE
                          EXIT SORT_EST3
                       ENDIF
                    ENDDO SORT_EST3
                    !
                 ENDIF
              ENDIF
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  IF(N /= NATOMS*(NSPECIES(0)-1)) THEN
     WRITE(MYUNIT,'(A)') 'scan_flip_nbrhd> Inconsistency!'
     STOP
  ENDIF
  !
  ! Now E0 and I_AB are no longer dummies!
  !
  E0 = POTEL
  X0(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  !
  IF(QALCS_PARAM > 0 .AND. QALCS_PARAM < N) THEN
     NTOT = QALCS_PARAM
  ELSE
     NTOT = N
  ENDIF
  !
  DO I=1,NTOT
     !
     IF(QALCSMODE==7) THEN    
        J = INT(DPRAND()*DBLE(N)) + 1
     ELSEIF(QALCSMODE==8) THEN
        J=I
     ENDIF
     !
     TA=INVATOMLISTS(FLIPS(J,1),1) ! Current label
     TB=FLIPS(J,2) ! new label
     CALL FLIP_LABEL(FLIPS(J,1),TB)
     !                                                            
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(QALCSV) THEN
        CALL PRINT_QUENCH(NP, ITER, '  ')
        IF(QALCSMODE==8) THEN
           WRITE(MYUNIT,'(2(A,F15.10))') &
                'scan_flip_nbrhd> dE*= ',EST(J)-E0,' dE= ',POTEL-E0 
        ENDIF
     ENDIF
     !
     IF(POTEL < E0 - ECONV) THEN
        SEQLENGTH = SEQLENGTH + 1
        WRITE(MYUNIT,'(A,I4,A,I4,A,I4,A,F15.10)') &
             'scan_flip_nbrhd> Trial ',I,' : ',&
             FLIPS(J,1),' ~> ',FLIPS(J,2),&
             ' => dE= ', POTEL-E0
        STEP = .TRUE.
        RETURN
     ELSE 
        ! Undo the flip and restore the configuration
        CALL FLIP_LABEL(FLIPS(J,1),TA)
        COORDS(1:3*NATOMS, NP) = X0(1:3*NATOMS)
        POTEL = E0
        !
        IF(QALCSMODE==7) THEN
           ! Remove the rejected swap from list
           FLIPS(J,1:2) = FLIPS(N,1:2)
           FLIPS(N,1:2) = 0
           N = N - 1
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  IF(NTOT == QALCS_PARAM) THEN
     WRITE(MYUNIT,'(A,I6,A)') &
          'scan_flip_nbrhd> Failed to find improvement after ',&
          NTOT,' trial swaps!'
  ELSE
     WRITE(MYUNIT,'(A)') 'scan_flip_nbrhd> No flips with -ve gain!'
  ENDIF
  STEP = .FALSE.
  !
  RETURN
  !
END SUBROUTINE SCAN_FLIP_NBRHD

! SUBROUTINE BH_SWAPS(NP, ITER, TIME, BRUN, QDONE, SCREENC)
!   !
!   USE COMMONS, ONLY : NSPECIES,ATOMLISTS,NATOMS,COORDS,NQ, &
!        ECONV,MYUNIT
!   !
!   IMPLICIT NONE
!   !
!   ! Parsed variables
!   INTEGER, INTENT(IN) :: NP
!   INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
!   DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
!   !
!   INTEGER :: N
!   !
!   !SWAPS: DO N=1,QALCS_NBHSWAPS
     
!   !ENDDO SWAPS
!   !
! END SUBROUTINE BH_SWAPS
