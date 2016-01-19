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
!   All procedures in this file were implemented by 
!   Dmitri Schebarchov (ds656). 
!
SUBROUTINE HOMOREF(NP, ITER, TIME, BRUN, QDONE, SCREENC)
  !
  ! ds656> 30/10/2013
  ! Iterative search for optimal permutation of current coordinates,
  ! guided entirely by flip gains, and guaranteeing local optimality
  ! when the flip gains are exact.
  !
  ! Global variables
  USE COMMONS, ONLY : NTYPEA, NATOMS, COORDS, NQ, NPAR, &
       HOMOREF_LSMODE,ECONV, TSTART, MYUNIT, HOMOREFTEST, &
       HIT, RMS, HOMOREF_NCYCLES, HOMOREF_AUXT, HOMOREF_FGMODE,&
       HOMOREF_AUX_NSWAPS, MCSTEPS, PRTFRQ, HOMOREFCHECK,&
       HOMOREF_BHT
  USE HOMOREFMOD, ONLY : NBONDS
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  ! NOTE:
  ! SCREENC serves no purpose, but is parsed back and forth throughout
  ! GMIN, starting from the main program. It is parsed by QUENCH (but
  ! not really used!), which is the only reason why it appears in HOMOREF.
  !
  ! Local variables
  LOGICAL :: COMPLETE_KL, GREEDY, STEEP, IMPROVED
  INTEGER :: N, LI, UI, NQTOT, NQINI
  INTEGER :: NSWAPS, NSWAPSTOT, NPASSES, NPASSESTOT, NCYCLES
  DOUBLE PRECISION :: POTEL, ESTART
  DOUBLE PRECISION :: ELOWEST, XLOWEST(3*NATOMS)
  DOUBLE PRECISION :: EMIN, XMIN(3*NATOMS)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  CALL MYCPU_TIME(TIME)
  WRITE(MYUNIT,'(A,F11.1)') 'homoref> Refining... t=',TIME-TSTART
  !
  IF(HOMOREF_LSMODE == 0 .OR. HOMOREF_LSMODE ==1) THEN
     GREEDY = .TRUE.
  ELSE
     GREEDY = .FALSE.
  ENDIF
  !
  IF(HOMOREF_LSMODE == 0) THEN
     STEEP = .TRUE.
  ELSE
     STEEP = .FALSE.
  ENDIF
  !
  ESTART = POTEL
  XLOWEST(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
  ELOWEST = POTEL ! This should be the energy of current COORDS.
  !
  NSWAPSTOT = 0
  NPASSESTOT = 0
  NCYCLES = 0
  !
  ! Reset NN bond count
  NBONDS(0:1) = 0
  !
  cycles: DO WHILE(NCYCLES < HOMOREF_NCYCLES)
     !
     NCYCLES = NCYCLES + 1
     IF(HOMOREF_LSMODE == -1) THEN
        COMPLETE_KL = .TRUE.
     ELSE
        COMPLETE_KL = .FALSE.
     ENDIF
     !                               
     XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
     EMIN = POTEL
     !
     NSWAPS = 0
     NPASSES = 0
     !
     NQINI = NQTOT
     !
     passes: DO WHILE( .NOT. COMPLETE_KL) ! Skipped for LSMODE=-1
        !
        COMPLETE_KL = .TRUE.
        NPASSES = NPASSES + 1
        !
        swaps: DO N=1, MIN(NTYPEA,NATOMS-NTYPEA)
           !
           NSWAPS = NSWAPS + 1
           !
           IF(STEEP) THEN
              ! Do not remove swapped atoms out of contention.
              LI = 1
              UI = NATOMS
           ELSE
              ! Set lower (LI) and upper (UI) indices bracketing
              ! atom candidates for swapping. 1 <= LI < UI <= NATOMS.
              ! This removes atoms from contention for the current pass
              LI = N
              UI = NATOMS - N + 1
           ENDIF
           !
           IF(HOMOREF_FGMODE==2) THEN
              CALL SCAN_SWAPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI,IMPROVED)
           ELSEIF(HOMOREF_FGMODE==1) THEN
              CALL FLIP_SEQ(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI)
           ELSEIF(HOMOREF_FGMODE==0) THEN
              CALL AFLIP_SEQ(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI)
           ELSE
              WRITE(MYUNIT,'(A, I2)') "homoref> Bad HOMOREF_FGMODE:", &
                   HOMOREF_FGMODE
              STOP
           ENDIF
           !
           ! Check if the swap has produced a beter minimum
           IF(POTEL < EMIN - ECONV) THEN
              !
              COMPLETE_KL = .FALSE.
              EMIN = POTEL
              XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
              !
              IF(HOMOREFTEST) THEN
                 WRITE(MYUNIT,'(A,F20.10,A,I4,A,I6,A,I3,A)') &
                      'homoref> New EMIN= ', &
                      EMIN, ' in cycle', NCYCLES, &
                      ' ( swap', NSWAPS,' pass', NPASSES, ' )'
              ENDIF
              !
           ELSEIF(GREEDY) THEN
              !
              ! Revert the last swap with +ve gain
              COORDS(1:3*NATOMS,NP) = XMIN(1:3*NATOMS)
              POTEL = EMIN
              !
              EXIT swaps ! also exit 'passes' if N=1, since COMPLETE_KL=TRUE
              !
           ENDIF
           !
           ! Terminate if the number of quenches exceeds MCSTEPS
           IF(NQ(NP) .GE. MCSTEPS(1)) THEN
              CALL MYCPU_TIME(TIME)
              WRITE(MYUNIT,'(A,I9,A,I9,A)') &
                   'homoref> NQ=', NQ(NP), &
                   ' exceeded MCSTEPS=', MCSTEPS(1),' quitting...'
              RETURN
           ENDIF
           !
        ENDDO swaps
        !
     ENDDO passes
     !
     N=0
     IF(HOMOREFCHECK) THEN
        WRITE(MYUNIT, '(A)') &
             'homoref> Checking and enforcing twofold optimality...'
        POTEL = EMIN
        COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
        IMPROVED=.TRUE.
        DO WHILE(IMPROVED)
           N=N+1
           XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
           EMIN = POTEL
           CALL SCAN_SWAPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,1,NATOMS,IMPROVED)
        ENDDO
        WRITE(MYUNIT,'(A,I4,A,F20.10)') &
             'homoref> Check passed after ',N,' scans, EMIN= ',EMIN
     ENDIF
     !
     NSWAPSTOT = NSWAPSTOT + NSWAPS + N
     NPASSESTOT = NPASSESTOT + NPASSES
     !
     !IF(HOMOREFTEST) THEN
     WRITE(MYUNIT,'(A,F20.10,A,I4,A,I8,A,I6,A,I2,A)') &
          'homoref> Converged to EMIN= ', &
          EMIN, ' in cycle', NCYCLES, &
          ' (',NQTOT-NQINI,' quenches', &
          NSWAPS,' swaps,', NPASSES, ' passes)'
     !ENDIF
     !
     ! Do nearest-neighbour analysis of current XMIN
     IF(NCYCLES == 1) THEN
        CALL ANAL_NNBROS(XMIN(1:3*NATOMS),.TRUE.,.TRUE.,.FALSE.)
     ELSE
        CALL ANAL_NNBROS(XMIN(1:3*NATOMS),.FALSE.,.TRUE.,.FALSE.)
     ENDIF
     !
     IF(EMIN < ELOWEST - ECONV) THEN
        ELOWEST = EMIN
        XLOWEST(1:3*NATOMS) = XMIN(1:3*NATOMS)
        IF(HOMOREFTEST) THEN
           WRITE(MYUNIT,'(A,F20.10,A,I4)') &
                'homoref> New ELOWEST= ', &
                EMIN, ' in cycle', NCYCLES
        ENDIF
     ENDIF
     !
     ! Reinstate current XLOWEST and ELOWEST
     POTEL = ELOWEST
     COORDS(1:3*NATOMS, NP) = XLOWEST(1:3*NATOMS)
     !
     ! Save XMIN and check if target is hit
     CALL GSAVEIT(EMIN, XMIN, NP)
     IF(HIT) EXIT cycles
     !
     IF(NCYCLES < HOMOREF_NCYCLES) THEN
        ! 
        CALL RANDPERM(NP)
        !
        NQTOT = NQTOT + 1
        NQ(NP) = NQ(NP) + 1
        CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
        IF(MOD(NQ(NP),PRTFRQ).EQ.0) THEN
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
        IF(HOMOREF_AUXT) THEN
           IF(HOMOREF_AUX_NSWAPS > 0) THEN
              CALL HOMOREF_AUX(NP, ITER, TIME, BRUN, QDONE, SCREENC)
           ENDIF
        ELSEIF(HOMOREF_BHT) THEN
           CALL HOMOREF_BH(NP, ITER, TIME, BRUN, QDONE, SCREENC)
        ENDIF
        !
     ENDIF
     !
  ENDDO cycles
  !
  CALL MYCPU_TIME(TIME)
  WRITE(MYUNIT,'(A,I5,A,I9,A,F15.10,A,F11.1)') &
       'homoref> Did', NCYCLES, ' cycle(s) with', &
       NSWAPSTOT, ' swaps (total). Final dE= ', &
       ELOWEST-ESTART, ' t= ',TIME-TSTART
  !
  !IF((NCYCLES > 0) .AND. (NCYCLES < HOMOREF_NCYCLES)) THEN 
     ! Print average NBRO table if the "cycles" loop exited early,
     ! which happens when a target is hit.
     CALL ANAL_NNBROS(XMIN(1:3*NATOMS),.FALSE.,.FALSE.,.TRUE.)
     CALL UPDATE_NNBOND_WEIGHTS()
  !ENDIF
  !
  RETURN
  !
END SUBROUTINE HOMOREF
!
!=====================================================================
!
SUBROUTINE FLIP_SEQ(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI)
  !
  ! ds656> 30/10/2013
  ! Perform a sequence of flips, guided by exact flip gains,
  ! stopping when a negative swap gain is encountered or all
  ! flip indices are depleted.
  !
  USE COMMONS, ONLY : NATOMS,NTYPEA,NSPECIES,NTYPEA_FIX,COORDS,&
       MYUNIT, HOMOREFTEST, NQ, RMS, NPAR, ECONV, &
       HOMOREF_NFMAX, PRTFRQ, HOMOREF_NSMAX
  USE HOMOREFMOD, ONLY : NUNIQFLIPS,IFLIPE
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP,LI,UI
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  LOGICAL :: COMPLETE, COMPLETE2, IMPROVED
  INTEGER :: IPICK,I,J,LINDEX,UINDEX,NFLIPS,NQTOT,LI2,UI2,NFLIPS_TOT,&
       NFLIPSA,NFLIPSB,NSEQ,FLIPHIST(0:HOMOREF_NFMAX)
  DOUBLE PRECISION :: DUMMY,POSTFLIPE(NATOMS),GAIN, ESTAR, &
       POTEL,E0,MEMORY(0:1),BESTSWAP_COORDS(3*NATOMS),&
       COORDS0(3*NATOMS),COORDS1(3*NATOMS),GRAD(3*NATOMS)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  ! Initial coordinates that remain unpermuted
  COORDS0(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
  ! Initial coordinates with same-species entries to be be permuted?
  COORDS1(1:3*NATOMS) = COORDS0(1:3*NATOMS)
  E0 = POTEL
  !
  BESTSWAP_COORDS(:) = COORDS1(:)
  !
  LI2=LI
  UI2=UI
  !
  NUNIQFLIPS = 0
  NFLIPS_TOT = 0
  NSEQ = 0
  COMPLETE2= .FALSE.
  IMPROVED= .FALSE.
  !
  FRES: DO WHILE(.NOT.COMPLETE2) 
     !
     FLIPHIST(0:HOMOREF_NFMAX)=0
     !
     NFLIPS = 0
     !
     LINDEX = LI2
     UINDEX = UI2
     !
     NSEQ = NSEQ + 1
     !
     COMPLETE = .FALSE.
     !
     FSEQ: DO WHILE(.NOT.COMPLETE)
        !
        COORDS(1:3*NATOMS,NP) = COORDS1(1:3*NATOMS)
        NFLIPS = NFLIPS + 1
        !
        IF(NFLIPS > 1 .OR. NFLIPS_TOT==0) THEN
           CALL CALC_POSTFLIPE(NP,ITER,TIME,BRUN,QDONE,SCREENC,&
                LINDEX,UINDEX,POSTFLIPE)
           IF(NFLIPS_TOT == 0) IFLIPE(:) = POSTFLIPE(:)
        ELSE
           POSTFLIPE(:) = IFLIPE(:)
        ENDIF        
        !
        ! Find atom index between LINDEX and UINDEX
        ! with the best (approximate) flip gain
        IPICK = LINDEX
        GAIN = POSTFLIPE(LINDEX)
        DO I = LINDEX+1, UINDEX
           DUMMY = POSTFLIPE(I)
           IF( DUMMY < GAIN ) THEN ! Update best gain
              GAIN = DUMMY
              IPICK = I
           ENDIF
        ENDDO
        !
        FLIPHIST(0) = FLIPHIST(0)+1
        FLIPHIST(FLIPHIST(0)) = IPICK
        !
        ! Check if atom IPICK is of type A or B, set I as adjacent
        ! to the partition, then update NTYPEA and the candidate 
        ! range for next flip. This is crafty, and it actually works.
        ! 
        IF(IPICK > NTYPEA_FIX) THEN ! type B, soon to be A
           I = NTYPEA_FIX + 1  ! New index for the ex-B now A-type atom
           NTYPEA = NTYPEA_FIX + 1
           LINDEX = LI2 ! Next time A-type atom will flip
           UINDEX = NTYPEA_FIX
        ELSE ! type A, soon to be B
           I = NTYPEA_FIX  ! New index for the ex-A now B-type atom
           NTYPEA = NTYPEA_FIX - 1 
           LINDEX = NTYPEA_FIX + 1 ! Next time B-type atom will flip
           UINDEX = UI2 
        ENDIF
        NSPECIES(1) = NTYPEA
        NSPECIES(2) = NATOMS-NTYPEA
        !
        ! If after two flips we encounter a best candidate already
        ! adjacent to the partition, then we have converged on the
        ! best pair swap candidates: NTYPEA_FIX and NTYPEA_FIX + 1
        ! If not, then continue the flip sequence. 
        !
        IF(NFLIPS > 2 .AND. I == IPICK) THEN ! Converged!!!
           !
           COMPLETE = .TRUE.
           !
        ELSE ! Continue flipping
           !
           ! Swap coordinates IPICK and I, with IPICK changing type.
           IF(I /= IPICK) THEN
              CALL SWAP_COORDS(NP,IPICK,I,.FALSE.)
              !CALL SWAP_COORDS_V2(COORDS0(1:3*NATOMS),IPICK,I)
              CALL SWAP_COORDS_V2(COORDS1(1:3*NATOMS),IPICK,I)
              CALL SWAP_IFLIPE(IPICK,I)
           ENDIF
           !
           IF(NFLIPS == 2) MEMORY(0) = GAIN  ! initial value              
           IF(NFLIPS >= 2) THEN
              MEMORY(1) = GAIN
              IF(MEMORY(1) <= MEMORY(0)) THEN
                 MEMORY(0) = MEMORY(1)
                 IF(MEMORY(0) < E0 - ECONV) THEN
                    IMPROVED = .TRUE.
                    BESTSWAP_COORDS(:) = COORDS1(:)
                 ENDIF
              ENDIF
              IF(NFLIPS >= MAX(HOMOREF_NFMAX,2)) THEN
                 IF(HOMOREF_NFMAX > 2) THEN
                    WRITE(MYUNIT,'(2A)') 'flip_seq> WARNING: ',&
                         'flip sequence did not converge'
                 ENDIF
                 COMPLETE = .TRUE. ! prepare for termination
              ENDIF
           ENDIF
           !
        ENDIF
        !
     ENDDO FSEQ
     !
     NFLIPS_TOT = NFLIPS_TOT + NFLIPS
     !
     COMPLETE2=.TRUE.
     !
     IF( .NOT. IMPROVED .AND. (NSEQ .LT. HOMOREF_NSMAX) &
          .AND. (LI==1) .AND. (UI==NATOMS) ) THEN 
        ! Process FLIPHIST and prepare for new sequence 
        ! First revert to NTYPE_FIX 
        NTYPEA = NTYPEA_FIX
        NSPECIES(1) = NTYPEA
        NSPECIES(2) = NATOMS-NTYPEA
        COORDS(1:3*NATOMS,NP) = COORDS0(1:3*NATOMS) ! revert COORDS
        CALL PROCESS_FLIPHIST(FLIPHIST,LI2,UI2,NP) ! Shuffle and adjust bounds
        COORDS0(:) = COORDS(:,NP) 
        COORDS1(:) = COORDS0(:)
        ! Adjust ATOMLISTS accordingly 
        CALL RESET_ATOMLISTS(1)
        !
        !IF(LI2 <= NTYPEA .AND. UI2 > NTYPEA) THEN
        IF(LI2 < NTYPEA .AND. UI2 > NTYPEA+1) THEN
           COMPLETE2=.FALSE. ! A new flip sequence will initiate
        ELSE
           WRITE(MYUNIT,'(A, 3(1X,I5))') &
                'flip_seq> Flipping indices depleted. LI2, UI2, NUNIQFLIPS=', &
                LI2, UI2, NUNIQFLIPS
        ENDIF
        !
     ELSE
        NUNIQFLIPS = 0
     ENDIF
     !
     IF(HOMOREFTEST) THEN
        WRITE(MYUNIT,'(A,I3,A,F15.10,A,I3,A)') &
             'flip_seq> ', NSEQ, &
             ' converged to swap gain ', MEMORY(0)-E0, &
             ' after ', NFLIPS, ' flips.'
        WRITE(MYUNIT,'(A)', ADVANCE='NO') 'flip_seq> FLIPHIST:'
        DO I=1,FLIPHIST(0)
           WRITE(MYUNIT,'(I5)', ADVANCE='NO') FLIPHIST(I)
        ENDDO
        WRITE(MYUNIT,'(I5)', ADVANCE='YES')
     ENDIF
     !
  ENDDO FRES
  !
  NTYPEA = NTYPEA_FIX
  NSPECIES(1) = NTYPEA
  NSPECIES(2) = NATOMS-NTYPEA
  COORDS(:,NP) = BESTSWAP_COORDS(:) ! revert to best swap
  !
  IF(IMPROVED) THEN
     CALL SWAP_COORDS(NP,NTYPEA,NTYPEA+1,.FALSE.)
     ! Push the swapped coordinates away from the partition boundary
     CALL SWAP_COORDS(NP,LI,NTYPEA,.FALSE.)
     CALL SWAP_COORDS(NP,NTYPEA+1,UI,.FALSE.)
  ENDIF
  !
  CALL RESET_ATOMLISTS(1)
  !
  ! Quenching after the swap will update COORDS and POTEL, but
  ! NOT check for placement in the saved set.
  NQTOT = NQTOT + 1
  NQ(NP) = NQ(NP) + 1
  CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
  IF(MOD(NQ(NP),PRTFRQ).EQ.0) THEN
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
  IF(HOMOREFTEST) THEN
     WRITE(MYUNIT,'(A,F15.10,A,I3,A)') &
          "flip_seq> Best swap gain of", &
          POTEL-E0, " after ", NFLIPS_TOT, " flips."
  ENDIF
  !
  RETURN
  !
END SUBROUTINE FLIP_SEQ
!
!======================================================================
!
SUBROUTINE CALC_POSTFLIPE(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI,POSTFLIPE)
  !
  ! ds656> 30/10/2013
  ! Calculate flip gains, allowing for local geometry relaxation.
  ! Indices LI and UI specifiy the range of atoms to be considered.
  !
  USE COMMONS, ONLY : NATOMS,NTYPEA,COORDS,NQ,MYUNIT,NPAR,RMS,&
       KEEPLISTS,NBRLISTS,NSPECIES,HOMOREFTEST,NTYPEA_FIX
  !
  IMPLICIT NONE
  !
  ! crap for QUENCH      
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE 
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) 
  !
  ! indices specifying processor ID and atom range to be "flipped"
  INTEGER, INTENT(IN) :: NP,LI,UI
  DOUBLE PRECISION, INTENT(OUT) :: POSTFLIPE(NATOMS)
  !
  INTEGER :: NBRLISTS0(2,NATOMS,0:NATOMS)
  INTEGER :: I,J,K,L1,L2,INBR,NQTOT,NTYPEA0
  DOUBLE PRECISION :: COORDS0(3*NATOMS),POTEL, &
       POTEL0,POT1,POT2,RC(2),GRAD(3*NATOMS)
  CHARACTER(LEN=2) :: ATYPE,TAG
  !   
  ! Energy of COORDS from last quench. Common block in QUENCH.   
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  COORDS0(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
  POTEL0 = POTEL
  NTYPEA0 = NTYPEA
  !
  IF(KEEPLISTS) THEN
     CALL BUILD_NBRLISTS(COORDS0)
     NBRLISTS0(1:2,1:NATOMS,0:NATOMS) = NBRLISTS(1:2,1:NATOMS,0:NATOMS)
  ENDIF
  !
  POSTFLIPE(:) = 0.0D0
  !
  DO I=LI,UI
     !
     ! Check if atom I is of type A or B and pick J near partition
     IF(I > NTYPEA0) THEN ! type B, soon to be A
        J = NTYPEA0 + 1
        NTYPEA = J
     ELSE ! type A, soon to be B
        J = NTYPEA0
        NTYPEA = NTYPEA0 - 1
     ENDIF
     NSPECIES(1) = NTYPEA
     NSPECIES(2) = NATOMS-NTYPEA
     !
     ! Now swap coordinates and NBR lists of I and J
     CALL SWAP_COORDS(NP,I,J,.TRUE.)
     !
     IF(KEEPLISTS) THEN
        !
        ! Put all atoms in group 3 and  sort by type with new NTYPEA.
        CALL RESET_ATOMLISTS(3) 
        !
        ! Make J and its 1st neighbour shell dynamic, 
        ! 2nd neighbour shell as frozen, and ignore the rest
        CALL SET_ATOMLISTS_BY_NBROS(J)
        !
     ELSE
        CALL RESET_ATOMLISTS(1) ! since NTYPEA changes
     ENDIF
     !
     ! Quench state with new NTYPEA (and maybe some frozen atoms).
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     ! 
     IF(HOMOREFTEST) THEN
        IF(NTYPEA == NTYPEA_FIX) THEN
           TAG='  '
        ELSE
           TAG='* '
        ENDIF
        IF (NPAR.GT.1) THEN
           WRITE(MYUNIT,'(A,I1,A,A,I10,A,F20.10,A,I5,A,G12.5,A)') &
                '[',NP,']Qu',TAG,NQ(NP),' E=',POTEL,' steps=',ITER, &
                ' RMS=',RMS
        ELSE
           WRITE(MYUNIT,'(A,A,I10,A,F20.10,A,I5,A,G12.5)') &
                'Qu',TAG,NQ(NP),' E=',POTEL,' steps=',ITER, &
                ' RMS=',RMS
        ENDIF
     ENDIF
     !
     ! Store the energy after each flip
     IF(KEEPLISTS) THEN
        !
        !CALL RESET_ATOMLISTS(1)
        CALL POTENTIAL(COORDS(1:3*NATOMS,NP),GRAD(1:3*NATOMS),POT2,.FALSE.,.FALSE.)
        POSTFLIPE(I) = POT2 ! - POT1
        !
     ELSE
        POSTFLIPE(I) = POTEL ! - POTEL0
        !
     ENDIF
     !
     ! Revert to original state by repeating the swap...
     CALL SWAP_COORDS(NP,I,J,.TRUE.)
     ! .. and restoring the coordinates to avoid wandering.
     COORDS(1:3*NATOMS,NP) = COORDS0(1:3*NATOMS) 
     !
  END DO
  !
  ! Revert to original energies
  NTYPEA = NTYPEA0
  NSPECIES(1) = NTYPEA
  NSPECIES(2) = NATOMS-NTYPEA
  POTEL = POTEL0
  !
  ! Unfreeze all atoms and sort by type with current NTYPEA.
  CALL RESET_ATOMLISTS(1) 
  !
  !WRITE(*,*) "CALC_POSTFLIPE> FINISH"
  !
  RETURN
  !      
END SUBROUTINE CALC_POSTFLIPE
!
!=====================================================================
SUBROUTINE AFLIP_SEQ(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI)
  !
  ! ds656> 30/10/2013
  ! Perform a sequence of flips guided by approximate flip gains.
  !
  USE COMMONS, ONLY : NATOMS,NTYPEA,NSPECIES,NTYPEA_FIX,COORDS,&
       MYUNIT, HOMOREFTEST, NQ, RMS, NPAR, ECONV, &
       HOMOREF_NFMAX, PRTFRQ, HOMOREF_NSMAX
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP,LI,UI
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) ! for QUENCH
  !
  LOGICAL :: COMPLETE, COMPLETE2
  INTEGER :: IPICK,I,J,LINDEX,UINDEX,NFLIPS,NQTOT,LI2,UI2,NFLIPS_TOT,&
       NFLIPSA,NFLIPSB,NSEQ,FLIPHIST(0:HOMOREF_NFMAX)
  DOUBLE PRECISION :: DUMMY,FLIPGAINS(NATOMS),GAIN, ESTAR, &
       POTEL,E0,MEMORY(0:1),BESTSWAP_COORDS(3*NATOMS),&
       COORDS0(3*NATOMS)
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  COORDS0(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)  
  E0 = POTEL
  !
  LI2=LI
  UI2=UI
  !
  NFLIPS_TOT = 0
  NSEQ = 0
  COMPLETE2= .FALSE.
  !
  FRES: DO WHILE(.NOT.COMPLETE2) 
     !
     FLIPHIST(0:HOMOREF_NFMAX)=0
     !
     NFLIPS = 0
     MEMORY(0:1) = 0.0d0 ! No real need to initialiuse here
     !
     LINDEX = LI2
     UINDEX = UI2
     !
     NSEQ = NSEQ + 1
     !
     COMPLETE = .FALSE.
     !
     FSEQ: DO WHILE(.NOT.COMPLETE)
        !
        NFLIPS = NFLIPS + 1
        !
        CALL BPOTENTIAL(COORDS(1:3*NATOMS,NP),DUMMY,FLIPGAINS)
        !
        ! Find atom index between LINDEX and UINDEX
        ! with the best (approximate) flip gain
        IPICK = LINDEX
        GAIN = FLIPGAINS(LINDEX)
        DO I = LINDEX+1, UINDEX
           DUMMY = FLIPGAINS(I)
           IF( DUMMY < GAIN ) THEN ! Update best gain
              GAIN = DUMMY
              IPICK = I
           ENDIF
        ENDDO
        !
        FLIPHIST(0) = FLIPHIST(0)+1
        FLIPHIST(FLIPHIST(0)) = IPICK
        !
        ! Check if atom IPICK is of type A or B, set I as adjacent
        ! to the partition, then update NTYPEA and the candidate 
        ! range for next flip. This is crafty, and it actually works.
        ! 
        IF(IPICK > NTYPEA_FIX) THEN ! type B, soon to be A
           I = NTYPEA_FIX + 1  ! New index for the ex-B now A-type atom
           NTYPEA = NTYPEA_FIX + 1
           LINDEX = LI2 ! Next time A-type atom will flip
           UINDEX = NTYPEA_FIX
        ELSE ! type A, soon to be B
           I = NTYPEA_FIX  ! New index for the ex-A now B-type atom
           NTYPEA = NTYPEA_FIX - 1 
           LINDEX = NTYPEA_FIX + 1 ! Next time B-type atom will flip
           UINDEX = UI2 
        ENDIF
        NSPECIES(1) = NTYPEA
        NSPECIES(2) = NATOMS-NTYPEA
        !
        ! If after two flips we encounter a best candidate already
        ! adjacent to the partition, then we have converged on the
        ! best pair swap candidates: NTYPEA_FIX and NTYPEA_FIX + 1
        ! If not, then continue the flip sequence. 
        !
        IF(NFLIPS > 2 .AND. I == IPICK) THEN ! Converged!!!
           !
           COMPLETE = .TRUE.
           !
        ELSE ! Continue flipping
           !
           ! Swap coordinates IPICK and I, with IPICK changing type.
           IF(I /= IPICK) CALL SWAP_COORDS(NP,IPICK,I,.FALSE.)
           !
           IF(NFLIPS == 2) THEN
              MEMORY(0) = ESTAR + GAIN  ! initial value
              MEMORY(1) = ESTAR + GAIN  ! initial value
              BESTSWAP_COORDS(:) = COORDS(:,NP)
           ENDIF
           IF(NFLIPS >= 2) THEN
              MEMORY(1) = ESTAR + GAIN 
              IF(MEMORY(1) < MEMORY(0)) THEN
                 MEMORY(0) = MEMORY(1)
                 BESTSWAP_COORDS(:) = COORDS(:,NP)
              ENDIF
              IF(NFLIPS >= MAX(HOMOREF_NFMAX,2)) THEN
                 IF(HOMOREF_NFMAX > 2) THEN
                    WRITE(MYUNIT,'(2A)') 'aflip_seq> WARNING: ',&
                         'flip sequence did not converge'
                 ENDIF
                 !COORDS(:,NP) = BESTSWAP_COORDS(:) ! revert to best swap
                 COMPLETE = .TRUE. ! prepare for termination
              ENDIF
           ENDIF
           !
           IF(.NOT. COMPLETE ) THEN
              ! Set ATOMLIST to all mobile and sorted by type
              CALL RESET_ATOMLISTS(1)
              IF(HOMOREF_NFMAX > 1) THEN
                 ! Quenching will update COORDS and POTEL
                 NQTOT = NQTOT + 1
                 NQ(NP) = NQ(NP) + 1
                 CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
                 IF(HOMOREFTEST) THEN
                    IF (NPAR.GT.1) THEN
                       WRITE(MYUNIT, &
                            '(A,I1,A,I10,A,F20.10,A,I5,A,G12.5,A,F11.1)') &
                            '[',NP,']Qu* ',NQ(NP),' E=',POTEL,' steps=',ITER, &
                            ' RMS=',RMS,' t=',TIME
                    ELSE
                       WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,F11.1)') &
                            'Qu* ',NQ(NP),' E=',POTEL,' steps=',ITER, &
                            ' RMS=',RMS,' t=',TIME
                    ENDIF
                 ENDIF
                 ESTAR=POTEL
              ENDIF
           ENDIF
           !
        ENDIF
        !
     ENDDO FSEQ
     !
     NFLIPS_TOT = NFLIPS_TOT + NFLIPS
     !
     COMPLETE2=.TRUE.
     !
     IF( (MEMORY(0)-E0 >= -ECONV) .AND. (NSEQ .LT. HOMOREF_NSMAX) ) THEN 
        ! Process FLIPHIST and prepare for new sequence 
        ! First revert to NTYPE_FIX 
        NTYPEA = NTYPEA_FIX
        NSPECIES(1) = NTYPEA
        NSPECIES(2) = NATOMS-NTYPEA
        COORDS(:,NP) = COORDS0(:) ! revert COORDS
        CALL PROCESS_FLIPHIST(FLIPHIST,LI2,UI2,NP)
        ! Adjust ATOMLISTS accordingly 
        CALL RESET_ATOMLISTS(1)
        !
        IF(LI2 >= NTYPEA .OR. UI2 <= NTYPEA+1) THEN
           WRITE(MYUNIT,'(A)') &
                'aflip_seq> Flipping indices depleted.'
        ELSE
           COMPLETE2=.FALSE. ! A new flip sequence will initiate
        ENDIF
        !
     ENDIF
     !
     IF(HOMOREFTEST) THEN
        WRITE(MYUNIT,'(A,I2,A,F15.10)') &
             'aflip_seq> ', NSEQ, &
             ' converged = ', MEMORY(0)-E0 
        WRITE(MYUNIT,'(A)', ADVANCE='NO') 'aflip_seq> FLIPHIST:'
        DO I=1,FLIPHIST(0)
           WRITE(MYUNIT,'(I5)', ADVANCE='NO') FLIPHIST(I)
        ENDDO
        WRITE(MYUNIT,'(I5)', ADVANCE='YES')
     ENDIF
     !
  ENDDO FRES
  !
  NTYPEA = NTYPEA_FIX
  NSPECIES(1) = NTYPEA
  NSPECIES(2) = NATOMS-NTYPEA
  COORDS(:,NP) = BESTSWAP_COORDS(:) ! revert to best swap
  CALL SWAP_COORDS(NP,NTYPEA,NTYPEA+1,.FALSE.)
  !
  ! Push the swapped coordinates away from the partition boundary
  CALL SWAP_COORDS(NP,LI,NTYPEA,.FALSE.)
  CALL SWAP_COORDS(NP,NTYPEA+1,UI,.FALSE.)
  !
  CALL RESET_ATOMLISTS(1)
  !
  ! Quenching after the swap will update COORDS and POTEL, but
  ! NOT check for placement in the saved set.
  NQTOT = NQTOT + 1
  NQ(NP) = NQ(NP) + 1
  CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
  IF(MOD(NQ(NP),PRTFRQ).EQ.0) THEN
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
  IF(HOMOREFTEST) THEN
     WRITE(MYUNIT,'(A,F13.8,A,I2,A)') &
          "aflip_seq> Best swap gain of", &
          POTEL-E0, " after ", NFLIPS_TOT, " flips."
  ENDIF
  !
  RETURN
  !
END SUBROUTINE AFLIP_SEQ
!
!=====================================================================
!
SUBROUTINE PROCESS_FLIPHIST(FLIPHIST,LI, UI, NP)
  ! ds656> Process FLIPHIST accumulated by FLIP_SEQ and do
  !        the necessary adjustament for another flip sequence for
  !        a smaller subset of atom indices.
  USE COMMONS, ONLY :  HOMOREF_NFMAX, MYUNIT, &
       NTYPEA, NTYPEA_FIX, HOMOREFTEST
  USE HOMOREFMOD, ONLY : NUNIQFLIPS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP,FLIPHIST(0:HOMOREF_NFMAX)
  INTEGER, INTENT(INOUT) :: LI,UI
  INTEGER :: M,N,I,J,FLIPS(0:HOMOREF_NFMAX+2),UNIQ(0:HOMOREF_NFMAX+2),&
       LI0,UI0,NA,NB
  LOGICAL :: DUPE, ADDNA, ADDNAP1
  !
  IF(NTYPEA /= NTYPEA_FIX ) THEN
     WRITE(MYUNIT, '(A)') 'process_fliphist> NTYPEA /= NTYPEA_FIX, stopping'
     STOP
  ENDIF
  !
  IF(HOMOREFTEST) THEN 
     WRITE(MYUNIT,'(A,2(1X,I3))') &
          'process_fliphist> Started with LI0,UI0:',LI,UI
     WRITE(MYUNIT,*) "process_fliphist> FLIPHIST:",(FLIPHIST(J),J=1,FLIPHIST(0))
  ENDIF
  !
  FLIPS(0:HOMOREF_NFMAX) = FLIPHIST(0:HOMOREF_NFMAX)
  FLIPS(HOMOREF_NFMAX+1:HOMOREF_NFMAX+2) = 0
  !
  UNIQ(0:HOMOREF_NFMAX+2) = 0
  ADDNA=.FALSE.
  ADDNAP1=.FALSE.
  !
  ! NTYPEA or NTYPE+1 should be considered
  ! if it appears in the first two entries of FLIPHIST.. 
  DO N=1,2
     IF(FLIPS(N) == NTYPEA) THEN
        ADDNA=.TRUE.
     ELSEIF(FLIPS(N) == NTYPEA+1) THEN
        ADDNAP1=.TRUE.
     ENDIF
  ENDDO
  !
  ! ..or (ii) if one or both of the first two entries have duplicates.
  L0: DO N=3,FLIPS(0)
     IF((FLIPS(N)==FLIPS(1)).OR.(FLIPS(N)==FLIPS(2))) THEN
        IF(FLIPS(N) <= NTYPEA) THEN 
           ADDNA=.TRUE.
        ELSEIF(FLIPS(N) > NTYPEA) THEN
           ADDNAP1=.TRUE.
        ENDIF
     ENDIF
  ENDDO L0
  !
  ! Add (potentially redundant) NTYPEA and NTYPEA+1 to FLIPHIST
  IF(ADDNA) THEN
     FLIPS(0)=FLIPS(0)+1
     FLIPS(FLIPS(0)) = NTYPEA
  ENDIF
  IF(ADDNAP1) THEN
     FLIPS(0)=FLIPS(0)+1
     FLIPS(FLIPS(0)) = NTYPEA+1
  ENDIF
  !
  NA=0
  NB=0
  ! Now sift through all entries 
  ! of FLIPS and ignore duplicates
  !L1: DO N=2,FLIPS(0)
  DO N=1,FLIPS(0)
     !
     I=FLIPS(N)
     DUPE=.FALSE.
     !
     IF( ((I.NE.NTYPEA).OR.ADDNA) .AND. ((I.NE.NTYPEA+1).OR.ADDNAP1) ) THEN
        !
        ! Ignore all instances of I=NTYPEA (unless ADDNA is true) and 
        ! all instances of I=NTYPEA+1 (unless ADDNAP1 is true).
        !
        L1: DO M=1,UNIQ(0)
           IF(I==UNIQ(M)) THEN
              DUPE=.TRUE.
              EXIT L1
           ENDIF
        ENDDO L1
        !
        IF(.NOT. DUPE) THEN ! STORE
           UNIQ(0) = UNIQ(0)+1
           UNIQ(UNIQ(0)) = I
           IF(I>NTYPEA) THEN
              NB=NB+1
           ELSE
              NA=NA+1
           ENDIF
           IF(UNIQ(0)>1) THEN ! SORT
              L2: DO M=UNIQ(0),2,-1
                 IF(UNIQ(M) < UNIQ(M-1)) THEN
                    J=UNIQ(M)
                    UNIQ(M)=UNIQ(M-1)
                    UNIQ(M-1)=J
                 ELSE
                    EXIT L2
                 ENDIF
              ENDDO L2
           ENDIF
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  IF(HOMOREFTEST) THEN 
     WRITE(MYUNIT,*) "process_fliphist> UNIQ:",(UNIQ(J),J=1,UNIQ(0))
     WRITE(MYUNIT,*) "process_fliphist> NA,NB:", NA,NB
  ENDIF
  NUNIQFLIPS = NUNIQFLIPS + UNIQ(0)
  !
  LI0=LI
  UI0=UI
  !
  DO N=1,NA
     I=UNIQ(N)
     IF(I <= NTYPEA .AND. LI <= NTYPEA .AND. LI <= I) THEN
        CALL SWAP_COORDS(NP,I,LI,.FALSE.)
        CALL SWAP_IFLIPE(I,LI)
        LI=LI+1        
     ELSE
        WRITE(MYUNIT,'(A,5(1X,I3))') &
             'process_fliphist> bad LI0 LI, I, UI, UI0:',&
             LI0,LI,I,UI,UI0
        WRITE(MYUNIT,*) "process_fliphist> UNIQ:",(UNIQ(J),J=1,UNIQ(0))
        WRITE(MYUNIT,*) "process_fliphist> FLIPS:",(FLIPHIST(J),J=1,FLIPHIST(0))
        STOP
     ENDIF     
  ENDDO
  !
  DO N=0,NB-1
     I=UNIQ(UNIQ(0)-N) ! Reverse order
     IF(I > NTYPEA .AND. UI > NTYPEA .AND. I <= UI) THEN
        CALL SWAP_COORDS(NP,I,UI,.FALSE.)
        CALL SWAP_IFLIPE(I,UI)
        UI=UI-1
     ELSE
        WRITE(MYUNIT,'(A,5(1X,I3))') &
             'process_fliphist> bad LI0 LI, I, UI, UI0:',&
             LI0,LI,I,UI,UI0
        WRITE(MYUNIT,*) "process_fliphist> UNIQ:",(UNIQ(J),J=1,UNIQ(0))
        WRITE(MYUNIT,*) "process_fliphist> FLIPS:",(FLIPHIST(J),J=1,FLIPHIST(0))
        STOP
     ENDIF
     
  ENDDO
  !
  !DO N=UNIQ(0),1,-1
  !   I=UNIQ(N)
  !   IF(I <= NTYPEA .AND. LI <= NTYPEA .AND. LI <= I) THEN
  !      CALL SWAP_COORDS(NP,I,LI,.FALSE.)
  !      LI=LI+1
  !   ELSEIF(I > NTYPEA .AND. UI > NTYPEA .AND. I <= UI) THEN
  !      CALL SWAP_COORDS(NP,I,UI,.FALSE.)
  !      UI=UI-1
  !   ELSE
  !      WRITE(MYUNIT,'(A,5(1X,I3))') &
  !           'process_fliphist> bad LI0 LI, I, UI, UI0:',&
  !           LI0,LI,I,UI,UI0
  !      WRITE(MYUNIT,*) "process_fliphist> UNIQ:",(UNIQ(J),J=1,UNIQ(0))
  !      WRITE(MYUNIT,*) "process_fliphist> FLIPS:",(FLIPHIST(J),J=1,FLIPHIST(0))
  !      STOP
  !   ENDIF
  !ENDDO
  !
  RETURN
  !
END SUBROUTINE PROCESS_FLIPHIST
!
!=====================================================================
! ds656> Sum over the indices of POTARRAY listed in LISTARRAY
!        and store the result in POTSUM.
!
SUBROUTINE SUM_VT_ENTRIES(POTARRAY,LISTARRAY,POTSUM)
  !
  USE COMMONS, ONLY : NATOMS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: LISTARRAY(0:NATOMS)
  DOUBLE PRECISION, INTENT(IN) :: POTARRAY(NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: POTSUM ! Sum over listed indices
  INTEGER :: I,IINDX
  !
  POTSUM=0.0D0
  !
  DO I=1,LISTARRAY(0)
     IINDX = LISTARRAY(I)
     POTSUM = POTSUM + POTARRAY(IINDX)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE SUM_VT_ENTRIES
!=====================================================================
!
! ds656> Set ATOMSLIST so that atom index I and all its neighbours 
!        are dynamic, and all other atoms are fixed. 
!        13/11/13
!
SUBROUTINE SET_ATOMLISTS_BY_NBROS(I)
  !
  USE COMMONS, ONLY : NTYPEA,NBRLISTS,ATOMLISTS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I ! Atom index
  INTEGER :: K,L,M,ITYPE,IGROUP,INBR,IINDX
  LOGICAL :: CHECK
  !
  !
  DO M = 1,2 ! neighbour shells 
     !
     DO K = 1,NBRLISTS(M,I,0) ! loop through neighbours of I
        !
        INBR = NBRLISTS(M,I,K)
        !
        IF(INBR > NTYPEA) THEN
           ITYPE = 2
        ELSE
           ITYPE = 1
        ENDIF
        !
        ! Sanity check:
        IF(INBR == 0) THEN
           WRITE(*,*) &
                "set_atomlists_by_nbros> INBR should not be zero:", &
                INBR
           STOP
        END IF
        !
        ! Grow group M.
        ATOMLISTS(ITYPE,M,0) = ATOMLISTS(ITYPE,M,0) + 1 
        ATOMLISTS(ITYPE,M,ATOMLISTS(ITYPE,M,0)) = INBR 
        !
     ENDDO
     !
  ENDDO
  !
  ! Sanity check
  IF(ATOMLISTS(1,1,0)+ATOMLISTS(2,1,0)+&
       ATOMLISTS(1,2,0)+ATOMLISTS(2,2,0) /= &
       NBRLISTS(1,I,0) + NBRLISTS(2,I,0)) THEN
     write(*,*) "Too many entries in ATOMLISTS:", &
          ATOMLISTS(1,1,0), ATOMLISTS(2,1,0),&
          ATOMLISTS(1,2,0), ATOMLISTS(2,2,0)
     WRITE(*,*) "index I:",I
     WRITE(*,*) "NBRLIST1:", (NBRLISTS(1,I,K),K=1,NBRLISTS(1,I,0))
     WRITE(*,*) "NBRLIST2:", (NBRLISTS(2,I,K),K=1,NBRLISTS(2,I,0))
     DO ITYPE=1,2
        DO IGROUP=1,2
           DO IINDX=1,ATOMLISTS(ITYPE,IGROUP,0)
              WRITE(*,*) "TYPE, GROUP, INDX, LABELS:", &
                   ITYPE,IGROUP,IINDX,ATOMLISTS(ITYPE,IGROUP,IINDX)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  ! Sanity check:
  ! CHECK=.TRUE.
  ! DO ITYPE=1,2 ! Type
  !    DO IGROUP=1,2 ! Group
  !       DO IINDX = 1,ATOMLISTS(ITYPE,IGROUP,0)
  !          IF (ATOMLISTS(ITYPE,IGROUP,IINDX) == 0) THEN
  !             WRITE(*,*) &
  !                  "set_atomlists_by_nbros> Last value is zero", &
  !                  ITYPE,IGROUP,IINDX,ATOMLISTS(ITYPE,IGROUP,IINDX)
  !             CHECK = .FALSE.
  !          ENDIF
  !       ENDDO
  !    ENDDO
  ! ENDDO
  ! IF(.NOT. CHECK) STOP
  !
  RETURN
  !
END SUBROUTINE SET_ATOMLISTS_BY_NBROS
!
!=====================================================================
!
SUBROUTINE BUILD_NBRLISTS(X)
  !
  USE COMMONS, ONLY : NATOMS,NBRLISTS,NBRCUT1,NBRCUT2
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  !
  INTEGER :: J1,J13,J2,J23,I
  DOUBLE PRECISION :: DIST,RC1,RC2
  !
  CALL RESET_NBRLISTS()
  !
  RC1 = NBRCUT1*NBRCUT1
  RC2 = NBRCUT2*NBRCUT2
  !
  DO J1 = 1,NATOMS-1
     J13 = 3*(J1-1)
     DO J2 = J1+1,NATOMS
        J23 = 3*(J2-1)
        DIST = 0.0D0
        DO I=1,3
           DIST = DIST + (X(J23+I)-X(J13+I))**2
        ENDDO
        !
        IF(DIST < RC2) THEN
           !
           IF(DIST < RC1) THEN
              NBRLISTS(1,J1,0) = NBRLISTS(1,J1,0) + 1
              NBRLISTS(1,J1,NBRLISTS(1,J1,0)) = J2
              NBRLISTS(1,J2,0) = NBRLISTS(1,J2,0) + 1
              NBRLISTS(1,J2,NBRLISTS(1,J2,0)) = J1
           ELSE
              !WRITE(*,*) NBRLISTS(2,J1,0), NBRLISTS(2,J2,0) 
              NBRLISTS(2,J1,0) = NBRLISTS(2,J1,0) + 1
              NBRLISTS(2,J1,NBRLISTS(2,J1,0)) = J2
              NBRLISTS(2,J2,0) = NBRLISTS(2,J2,0) + 1
              NBRLISTS(2,J2,NBRLISTS(2,J2,0)) = J1
           ENDIF
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE BUILD_NBRLISTS
!=====================================================================
!
SUBROUTINE SWAP_COORDS(NP,I,J,SWAPLISTS)
  !
  USE COMMONS, ONLY : COORDS, KEEPLISTS, NBRLISTS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP,I,J
  LOGICAL, INTENT(IN) :: SWAPLISTS
  INTEGER :: K,L,M
  DOUBLE PRECISION :: DUMMY
  !
  K = 3*(I-1)
  L = 3*(J-1)
  DO M = 1,3
     DUMMY = COORDS(K+M,NP)
     COORDS(K+M,NP) = COORDS(L+M,NP)
     COORDS(L+M,NP) = DUMMY
  ENDDO
  !
  ! Update NBRLISTS accordingly. MESSY!!!!
  !
  ! WRITE(*,*) "1st neighbours of",I," before coordinate swap"
  ! WRITE(*,*) (NBRLISTS(1,I,K),K=1,NBRLISTS(1,I,0))
  ! WRITE(*,*) "2nd neighbours of",I," before coordinate swap"
  ! WRITE(*,*) (NBRLISTS(2,I,K),K=1,NBRLISTS(2,I,0))
  
  ! WRITE(*,*) "1st neighbours of",J," before coordinate swap"
  ! WRITE(*,*) (NBRLISTS(1,J,K),K=1,NBRLISTS(1,J,0))
  ! WRITE(*,*) "2nd neighbours of",J," before coordinate swap"
  ! WRITE(*,*) (NBRLISTS(2,J,K),K=1,NBRLISTS(2,J,0))
  !
  IF(KEEPLISTS .AND. SWAPLISTS) THEN
     !
     DO M=1,2 ! loop over the 1st and 2nd neighbour shells
        DO K=0,MAX(NBRLISTS(M,I,0),NBRLISTS(M,J,0))
           !
           ! Exchange I and J lists, but NOT the 1st entry for M=1
           IF(M /= 1 .OR. K /= 1) THEN ! had .AND. here before... d'oh!
              !
              L=NBRLISTS(M,I,K)
              !
              ! Had 1 instead of 0 here... d'oh!
              IF(K > 0 .AND. NBRLISTS(M,J,K) == I) THEN 
                 NBRLISTS(M,I,K) = J 
              ELSE
                 NBRLISTS(M,I,K) = NBRLISTS(M,J,K)
              ENDIF
              !
              IF(K > 0 .AND. L == J) THEN
                 NBRLISTS(M,J,K) = I
              ELSE
                 NBRLISTS(M,J,K) = L
              ENDIF
              !
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! WRITE(*,*) "1st neighbours of",I," after coordinate swap"
  ! WRITE(*,*) (NBRLISTS(1,I,K),K=1,NBRLISTS(1,I,0))
  ! WRITE(*,*) "2nd neighbours of",I," after coordinate swap"
  ! WRITE(*,*) (NBRLISTS(2,I,K),K=1,NBRLISTS(2,I,0))
  ! !
  ! WRITE(*,*) "1st neighbours of",J," after coordinate swap"
  ! WRITE(*,*) (NBRLISTS(1,J,K),K=1,NBRLISTS(1,J,0))
  ! WRITE(*,*) "2nd neighbours of",J," after coordinate swap"
  ! WRITE(*,*) (NBRLISTS(2,J,K),K=1,NBRLISTS(2,J,0))
  !
  !STOP 
  RETURN
  !
END SUBROUTINE SWAP_COORDS
!
!=====================================================================
!
SUBROUTINE SWAP_COORDS_V2(X,I,J)
  !
  USE COMMONS, ONLY : NATOMS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I,J
  DOUBLE PRECISION, INTENT(INOUT) :: X(1:3*NATOMS)
  INTEGER :: K,L,M
  DOUBLE PRECISION :: DUMMY
  !
  K = 3*(I-1)
  L = 3*(J-1)
  DO M = 1,3
     DUMMY = X(K+M)
     X(K+M) = X(L+M)
     X(L+M) = DUMMY
  ENDDO
  !
END SUBROUTINE SWAP_COORDS_V2
!
!=====================================================================
!
SUBROUTINE SWAP_IFLIPE(I,J)
  !
  USE HOMOREFMOD, ONLY : IFLIPE
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I,J
  DOUBLE PRECISION :: DUMMY
  !
  DUMMY = IFLIPE(I)
  IFLIPE(I) = IFLIPE(J)
  IFLIPE(J) = DUMMY
  !
END SUBROUTINE SWAP_IFLIPE
!
!=====================================================================
!
SUBROUTINE BPOTENTIAL(X,E,ALFGAIN)
  !
  USE COMMONS, ONLY : NATOMS, BGUPTAT, BLJCLUSTER_NOCUT, MYUNIT, VT, &
       GLJT, NSPECIES, NPCALL
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: E, ALFGAIN(NATOMS)
  DOUBLE PRECISION :: GRAD(3*NATOMS) ! not used but must be parsed
  !
  !
  NPCALL=NPCALL+1 ! Update function calls
  !
  IF(BLJCLUSTER_NOCUT .OR. (GLJT .AND. NSPECIES(0)==2)) THEN
     CALL BLJ_CLUST(X,GRAD,E,.FALSE.) ! updates global per-atom energies
     CALL BLJ_CLUST2(X,ALFGAIN)
  ELSEIF(BGUPTAT) THEN
     CALL BGUPTA(X,GRAD,E,.FALSE.) ! updates global per-atom energies
     CALL BGUPTA2(X,ALFGAIN)
  ELSE
     WRITE(MYUNIT, '(A)') "homoref> ERROR: BPOTENTIAL either dislikes current potential, or there are more than 2 species!"
     STOP
  ENDIF
  !
  ALFGAIN(:) = ALFGAIN(:) - VT(:)
  !
  RETURN
  !
END SUBROUTINE BPOTENTIAL
!======================================================================
!
SUBROUTINE SCAN_SWAPS(NP,ITER,TIME,BRUN,QDONE,SCREENC,LI,UI,IMP)
  !
  ! ds656> 19/3/2014
  ! Evaluation of all swap gains and search for the lowest one.
  ! Indices LI and UI specifiy the range of atoms to be considered.
  !
  USE COMMONS, ONLY : NATOMS,NTYPEA,COORDS,NQ,MYUNIT,NPAR,RMS,ECONV,&
       HOMOREFTEST
  !
  IMPLICIT NONE
  !
  ! crap for QUENCH      
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE 
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS) 
  !
  ! indices specifying processor ID and atom range to be considered for swap
  INTEGER, INTENT(IN) :: NP,LI,UI
  LOGICAL, INTENT(OUT) :: IMP
  !
  INTEGER :: I,J,II,JJ,NSWAPS,NQTOT
  DOUBLE PRECISION :: COORDS0(3*NATOMS),POTEL,POTEL0,POTEL_LOWEST,&
       COORDS_LOWEST(3*NATOMS)
  !   
  ! Energy of COORDS from last quench. Common block in QUENCH.   
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  COORDS0(:) = COORDS(:,NP)
  POTEL0 = POTEL
  !
  NSWAPS=0
  DO I=LI,NTYPEA
     DO J=NTYPEA+1,UI
        !
        NSWAPS = NSWAPS + 1
        !
        ! Now swap coordinates and NBR lists of I and J
        CALL SWAP_COORDS(NP,I,J,.TRUE.)
        !CALL RESET_ATOMLISTS(1) ! Is this needed?
        !
        ! Quench state with new NTYPEA (and maybe some frozen atoms).
        NQTOT = NQTOT + 1
        NQ(NP) = NQ(NP) + 1
        CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
        !
        IF(HOMOREFTEST) THEN
           IF (NPAR.GT.1) THEN
              WRITE(MYUNIT,'(A,I1,A,I10,A,F20.10,A,I5,A,G12.5)') &
                   '[',NP,']Qu',NQ(NP),' E=',POTEL,' steps=',ITER, &
                   ' RMS=',RMS
           ELSE
              WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5)') &
                   'Qu',NQ(NP),' E=',POTEL,' steps=',ITER, &
                   ' RMS=',RMS
           ENDIF
        ENDIF
        !
        IF((NSWAPS == 1).OR.(POTEL < POTEL_LOWEST-ECONV)) THEN
           POTEL_LOWEST = POTEL
           COORDS_LOWEST(:) = COORDS(:,NP)
           II=I
           JJ=J
        ENDIF
        !
        ! Revert to original state energies
        COORDS(:,NP) = COORDS0(:)
        POTEL = POTEL0
        !CALL RESET_ATOMLISTS(1) ! is this needed?
        !
     END DO
  END DO
  !
  IF(HOMOREFTEST) THEN
     WRITE(MYUNIT,'(A,F15.10,A,2(1X,I4))') &
          'scan_swaps> Best swap gain of', POTEL_LOWEST-POTEL0,&
          'involving ', II, JJ
     
  ENDIF
  IF(POTEL_LOWEST<POTEL0-ECONV) THEN
     IMP=.TRUE.
  ELSE
     IMP=.FALSE.
  ENDIF
  !
  COORDS(:,NP) = COORDS_LOWEST(:)
  POTEL = POTEL_LOWEST
  ! Push the swapped coordinates away from the partition boundary
  CALL SWAP_COORDS(NP,LI,II,.FALSE.)
  CALL SWAP_COORDS(NP,JJ,UI,.FALSE.)
  !
  RETURN
  !      
END SUBROUTINE SCAN_SWAPS
