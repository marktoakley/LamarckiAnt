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
!====================================================================
!   This file contains supplementary routines for homotop refinement.
!   Dmitri Schebarchov (ds656).
!====================================================================
!
MODULE HOMOREFMOD
  !
  USE COMMONS, ONLY : NATOMS, NSPECIES
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: NNMAX=12
  INTEGER :: NBONDS(0:1), NUNIQFLIPS
  DOUBLE PRECISION :: NBONDS_AVE
  !
  ! 1st index: 1..NATMOS 
  ! 2nd index: 1..NSPECIES(0)
  ! 3rd index: 0..NNMAX
  INTEGER, ALLOCATABLE :: NNBRS(:,:,:)
  !
  ! 1st index: 1..NSPECIES(0)
  ! 2nd index: 1..NSPECIES(0)
  ! 3rd index:-1..NNMAX
  INTEGER, ALLOCATABLE :: NNHIST(:,:,:)
  INTEGER, ALLOCATABLE :: ANNHIST(:,:,:)
  !
  ! 1st index: 1..NSPECIES(0)
  ! 2nd index: 1..NSPECIES(0)
  DOUBLE PRECISION, ALLOCATABLE :: ANNHIST_MEAN(:,:), ANNHIST_VAR(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: NNBOND_WEIGHTS(:,:)
  !
  ! 1..NATOMS
  DOUBLE PRECISION, ALLOCATABLE :: IFLIPE(:)
  !
END MODULE HOMOREFMOD
!
!====================================================================
!
SUBROUTINE ANAL_NNBROS(X,INI,HIST,FIN) ! N-species
  !
  ! Nearest neightbour analysis of coordinates in X.
  ! INI -> initialise ANNHIST (Y/N)
  ! HIST -> update ANNHIST (Y/N)
  ! FIN -> print final ANNHIST (Y/N)
  !
  USE COMMONS, ONLY : NATOMS, NSPECIES, MYUNIT,&
       HOMOREF_AUX_NNCUT
  USE HOMOREFMOD, ONLY : NNMAX, NNHIST, ANNHIST, NNBRS, &
       ANNHIST_MEAN, ANNHIST_VAR, NBONDS
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(1:3*NATOMS)
  LOGICAL, INTENT(IN) :: INI, FIN, HIST
  !
  LOGICAL :: STORE
  INTEGER :: J1,J13,J2,J23,I,K,T1,T2,J1REM, J2REM
  DOUBLE PRECISION :: DIST2, CUTOFF2, RIJ2(NATOMS,NATOMS)
  !
  ! Initialise average histogram
  IF(INI) THEN
     ANNHIST(1:NSPECIES(0),1:NSPECIES(0),-1:NNMAX) = 0
  ENDIF
  !
  CUTOFF2 = HOMOREF_AUX_NNCUT**2
  RIJ2(:,:) = 0.0D0
  !
  ! Reset nearest-nehgbour lists
  NNBRS(:,:,:) = 0
  !  
  NBONDS(1)=0
  DO J1=1,NATOMS-1
     J13=3*(J1-1)
     DO J2=J1+1,NATOMS
        J23=3*(J2-1)
        !
        DIST2=0.0D0
        DO I=1,3
           DIST2 = DIST2 + (X(J23+I)-X(J13+I))**2
        ENDDO
        RIJ2(J1,J2) = DIST2
        RIJ2(J2,J1) = DIST2
        !
        IF(DIST2 < CUTOFF2) THEN
           !
           ! Since atoms have no labels, infer
           ! their type from atom indices...
           I=0
           L1: DO K=1,NSPECIES(0)
              IF(J1 <= NSPECIES(K) + I) THEN
                 T1 = K
                 EXIT L1
              ENDIF
              I=I+NSPECIES(K)
           ENDDO L1
           I=0
           L2: DO K=T1,NSPECIES(0) ! Know that T1 <= T2
              IF(J2 <= NSPECIES(K) + I) THEN
                 T2 = K
                 EXIT L2
              ENDIF
              I=I+NSPECIES(K)
           ENDDO L2
           !
           STORE=.TRUE.
           J1REM=0
           J2REM=0
           !
           IF(NNBRS(J1,T2,0) < NNMAX) THEN
              CONTINUE ! do nothing
           ELSEIF(NNBRS(J1,T2,0) == NNMAX) THEN
              IF( DIST2 < RIJ2(J1,NNBRS(J1,T2,NNMAX)) ) THEN
                 J1REM=NNBRS(J1,T2,NNMAX)
              ELSE
                 STORE=.FALSE.
              ENDIF
           ELSEIF(NNBRS(J1,T2,0) > NNMAX) THEN
              WRITE(MYUNIT,'(A)') 'anal_nbros> NN count exceeds NNMAX!'
              STOP
           ENDIF
           !
           IF(NNBRS(J2,T1,0) < NNMAX) THEN
              CONTINUE ! do nothing
           ELSEIF(NNBRS(J2,T1,0) == NNMAX) THEN
              IF(DIST2 < RIJ2(J2,NNBRS(J2,T1,NNMAX)) ) THEN
                 J2REM=NNBRS(J2,T1,NNMAX)
              ELSE
                 STORE=.FALSE.
              ENDIF
           ELSEIF(NNBRS(J1,T2,0) > NNMAX) THEN
              WRITE(MYUNIT,'(A)') 'anal_nbros> NN count exceeds NNMAX!'
              STOP   
           ENDIF
           !
           IF(STORE) THEN
              !
              IF(J1REM == 0) THEN ! simply grow NNBRS(J1)
                 NBONDS(1) = NBONDS(1) + 1
                 NNBRS(J1,T2,0) = NNBRS(J1,T2,0) + 1
                 NNBRS(J1,T2,NNBRS(J1,T2,0)) = J2
              ELSE ! Substitute the last value of NNBRS(J1)..
                 NNBRS(J1,T2,NNBRS(J1,T2,0)) = J2
                 ! ..and remove atom J1 from NNBRS(J1REM)
                 DO I=1,NNBRS(J1REM,T1,0)
                    IF(NNBRS(J1REM,T1,I) == J1) EXIT
                 ENDDO
                 ! Shuffle other entries in NNBRS(J1REM) to fill void
                 DO K=I,NNBRS(J1REM,T1,0)-1
                    NNBRS(J1REM,T1,K) = NNBRS(J1REM,T1,K+1)
                 ENDDO
                 ! Decrement NNBRS(J1REM) count
                 NNBRS(J1REM,T1,NNBRS(J1REM,T1,0))=0
                 NNBRS(J1REM,T1,0) = NNBRS(J1REM,T1,0) - 1
                 NBONDS(1) = NBONDS(1) - 1
              ENDIF
              !
              IF(J2REM == 0) THEN ! simply grow NNBRS(J2)
                 NBONDS(1) = NBONDS(1) + 1
                 NNBRS(J2,T1,0) = NNBRS(J2,T1,0) + 1
                 NNBRS(J2,T1,NNBRS(J2,T1,0)) = J1
              ELSE ! Substitute the last value of NNBRS(J2)
                 NNBRS(J2,T1,NNBRS(J2,T1,0)) = J1
                 ! ..and remove atom J2 from NNBRS(J2REM)
                 DO I=1,NNBRS(J2REM,T2,0)
                    IF(NNBRS(J2REM,T2,I) == J2) EXIT
                 ENDDO
                 ! Shuffle other entries in NNBRS(J2REM) to fill void
                 DO K=I,NNBRS(J2REM,T2,0)-1
                    NNBRS(J2REM,T2,K) = NNBRS(J2REM,T2,K+1)
                 ENDDO
                 ! Decrement NNBRS(J2REM)
                 NNBRS(J2REM,T2,NNBRS(J2REM,T2,0))=0
                 NNBRS(J2REM,T2,0) = NNBRS(J2REM,T2,0) - 1
                 NBONDS(1) = NBONDS(1) - 1
              ENDIF
              !
              ! Shuffle down the newly appended values in 
              ! NNBRS(J1) and NNBRS(J2). This makes the
              ! two lists ORDERED by increasing distance.
              !
              NNJ1: DO I=NNBRS(J1,T2,0),2,-1
                 IF( RIJ2(J1,NNBRS(J1,T2,I)) < &
                      RIJ2(J1,NNBRS(J1,T2,I-1)) ) THEN
                    K=NNBRS(J1,T2,I)
                    NNBRS(J1,T2,I) = NNBRS(J1,T2,I-1)
                    NNBRS(J1,T2,I-1) = K
                 ELSE
                    EXIT NNJ1
                 ENDIF
              ENDDO NNJ1
              !
              NNJ2: DO I=NNBRS(J2,T1,0),2,-1
                 IF( RIJ2(J2,NNBRS(J2,T1,I)) < &
                      RIJ2(J2,NNBRS(J2,T1,I-1)) ) THEN
                    K=NNBRS(J2,T1,I)
                    NNBRS(J2,T1,I) = NNBRS(J2,T1,I-1)
                    NNBRS(J2,T1,I-1) = K
                 ELSE
                    EXIT NNJ2
                 ENDIF
              ENDDO NNJ2
              !
           ENDIF
           !
           ! Sanity check
           IF( (NNBRS(J1,T2,0) > NNMAX) .OR. &
                (NNBRS(J2,T1,0) > NNMAX) ) THEN
              WRITE(MYUNIT,'(A)') 'anal_nbros> NN count exceeds NNMAX!'
              STOP
           ENDIF
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  IF(NBONDS(0) > 0 .AND. NBONDS(0) /= NBONDS(1)) THEN
     WRITE(MYUNIT,'(A, I5, A, I5)') &
          'anal_nbros> WARNING: 2xNBONDS changed from', &
          NBONDS(0),' to', NBONDS(1)
  ENDIF
  NBONDS(0) = NBONDS(1)
  !
  IF(HIST) THEN ! Accumulate bond-count histograms     
     !
     NNHIST(:,:,:) = 0
     !
     I=0
     DO T1=1,NSPECIES(0)
        DO J1=I+1,I+NSPECIES(T1)
           DO T2=1,NSPECIES(0)
              NNHIST(T1,T2,NNBRS(J1,T2,0)) = &
                   NNHIST(T1,T2,NNBRS(J1,T2,0)) + 1 ! T1-T2 NNs for T1
              NNHIST(T1,T2,-1)=NNHIST(T1,T2,-1) + 1 ! total 
           ENDDO
        ENDDO
        I=I+NSPECIES(T1)
     ENDDO
     !
     WRITE(MYUNIT,'(A,I7)') 'anal_nbros> 2xNBONDS= ', NBONDS(1) 
     WRITE(MYUNIT,'(3(A))')'anal_nbros> Count:  ',&
          '|  0  |  1  |  2  |  3  |  4  |  5  |  6  ',&
          '|  7  |  8  |  9  |  10 | 11  | 12  | TOT '
     DO T1=1,NSPECIES(0)
        DO T2=1,NSPECIES(0)
           WRITE(MYUNIT,'(A,I1,A,I1,A,14(1X,I5))') 'anal_nbros> T', &
                T1,'_T',T2,':',(NNHIST(T1,T2,I), I=0,NNMAX),&
                NNHIST(T1,T2,-1)
        ENDDO
     ENDDO
     !
     ! Update cumulative average histograms
     DO T1=1,NSPECIES(0)
        DO T2=1,NSPECIES(0)
           DO I=0,NNMAX ! each coordination
              ANNHIST(T1,T2,I)=ANNHIST(T1,T2,I)+NNHIST(T1,T2,I)
           ENDDO
           ANNHIST(T1,T2,-1)=ANNHIST(T1,T2,-1)+NNHIST(T1,T2,-1) ! Norm
        ENDDO
     ENDDO
     !
     ! Update means and variances
     DO J1=1,NSPECIES(0) 
        DO J2=1,NSPECIES(0)
           !
           ANNHIST_MEAN(J1,J2) = 0.0D0
           DO I=1,NNMAX ! Exclude I=0 since the contribution is 0
              ANNHIST_MEAN(J1,J2) = ANNHIST_MEAN(J1,J2) + &
                   DBLE(I*ANNHIST(J1,J2,I))
           ENDDO
           ANNHIST_MEAN(J1,J2) = ANNHIST_MEAN(J1,J2) / &
                DBLE(ANNHIST(J1,J2,-1))
           !
           ANNHIST_VAR(J1,J2) = 0.0D0
           DO I=0,NNMAX
              ANNHIST_VAR(J1,J2) = ANNHIST_VAR(J1,J2) + &
                   (DBLE(I) - ANNHIST_MEAN(J1,J2))**2 * &
                   DBLE(ANNHIST(J1,J2,I))                
           ENDDO
           ANNHIST_VAR(J1,J2) = ANNHIST_VAR(J1,J2) / &
                DBLE(ANNHIST(J1,J2,-1))
        !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  IF(FIN) THEN ! Print average histograms
     WRITE(MYUNIT,'(3(A))')'anal_nbros> AVDIS:',&
          '   0     1     2     3     4     5     6  ',&
          '   7     8     9     10   11    12  |  MEAN  |  VAR'
     DO T1=1,NSPECIES(0)
        DO T2=1,NSPECIES(0)
           WRITE(MYUNIT,'(A,I1,A,I1,A,13(1X,F5.3),2(1X,F7.4))') &
                'anal_nbros> T',T1,'_T',T2,':', &
                (DBLE(ANNHIST(T1,T2,I))/DBLE(ANNHIST(T1,T2,-1)), &
                I=0,NNMAX), &
                ANNHIST_MEAN(T1,T2), ANNHIST_VAR(T1,T2)
        ENDDO
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE ANAL_NNBROS
!
!=====================================================================
! Auxilliary basin-hopping with random (unweighted) label swaps
!
SUBROUTINE HOMOREF_BH(NP, ITER, TIME, BRUN, QDONE, SCREENC)
  !
  USE COMMONS, ONLY : NATOMS, NTYPEA, MYUNIT, COORDS, NQ, NPAR, RMS,&
       HOMOREF_BH_NSWAPMAX, HOMOREF_BH_NDUDMAX, HOMOREF_BH_TEMP, &
       HOMOREF_BH_FACTOR, ECONV, PRTFRQ, HOMOREFTEST
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT):: ITER,BRUN,QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME,SCREENC(3*NATOMS) ! QUENCH 
  !
  LOGICAL :: COMPLETE
  INTEGER :: NSWAPS,NDUDS,IA,IB,J,K,L,NQTOT
  DOUBLE PRECISION :: DPRAND,R,COORDS0(1:3*NATOMS),E0,POTEL,&
       XMIN(1:3*NATOMS),EMIN, TEMP
  !
  ! Energy of COORDS from last quench. Common block in QUENCH. 
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  EMIN = POTEL
  !
  WRITE(MYUNIT,'(A, F16.8)') &
       'homoref_bh> Initial E= ',EMIN
  !
  TEMP = HOMOREF_BH_TEMP
  NSWAPS=0
  NDUDS=0
  COMPLETE = .FALSE.
  !
  SWAPS: DO WHILE (.NOT. COMPLETE) 
     !
     NSWAPS=NSWAPS+1
     NDUDS=NDUDS+1
     !
     COORDS0(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
     E0 = POTEL
     !
     IA=INT(DPRAND()*DBLE(NTYPEA)) + 1
     IB=INT(DPRAND()*DBLE(NATOMS-NTYPEA)) + 1 + NTYPEA
     !
     ! Swap coords of IPICKA and IPICKB
     J=3*(IA-1)
     K=3*(IB-1)
     DO L=1,3
        R=COORDS(J+L,NP)
        COORDS(J+L,NP) = COORDS(K+L,NP)
        COORDS(K+L,NP) = R
     ENDDO
     !
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(HOMOREFTEST) THEN
        IF (MOD(NSWAPS-1,PRTFRQ).EQ.0) THEN
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
     ENDIF
     !
     R=DEXP((E0-POTEL)/TEMP)
     IF(DPRAND() < R) THEN ! Accept
        IF(POTEL < EMIN - ECONV) THEN
           EMIN=POTEL
           XMIN(1:3*NATOMS)=COORDS(1:3*NATOMS, NP)
           WRITE(MYUNIT,'(A,F10.5,A,I5)') 'homoref_bh> New EMIN:', &
                EMIN,' on swap', NSWAPS
           NDUDS = 0
        ENDIF
     ELSE ! Reject and revert
        COORDS(1:3*NATOMS,NP) = COORDS0(1:3*NATOMS)
        POTEL = E0
     ENDIF
     !
     TEMP = TEMP*HOMOREF_BH_FACTOR
     !
     IF(  (NDUDS .GE. HOMOREF_BH_NDUDMAX) .OR. &
          (NSWAPS .GE. HOMOREF_BH_NSWAPMAX)  ) &
          COMPLETE = .TRUE.
     !
  END DO SWAPS
  !
  POTEL = EMIN
  COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
  !
  WRITE(MYUNIT,'(A, F16.8,A,I5,A)') &
       'homoref_bh> Final E=',EMIN,&
       ' after', NSWAPS,' swaps'
  !
  RETURN
  !
END SUBROUTINE HOMOREF_BH
!
!=====================================================================
!
SUBROUTINE HOMOREF_AUX(NP, ITER, TIME, BRUN, QDONE, SCREENC) ! binary
  !
  ! Auxuliary search by rejection-free BH, with each atom weighted
  ! by a Boltzman factor with a change in weighted bond count as the
  ! the argument.
  !
  USE HOMOREFMOD, ONLY : NNMAX
  USE COMMONS, ONLY : NATOMS, NTYPEA, MYUNIT, COORDS, NQ, NPAR, RMS,&
       HOMOREF_AUX_NSWAPS, HOMOREF_AUX_TEMP, HOMOREF_AUX_FACTOR
  !
  IMPLICIT NONE
  !
  ! Parsed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT):: ITER,BRUN,QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME,SCREENC(3*NATOMS) ! QUENCH
  !
  LOGICAL :: SUCCESS
  INTEGER :: I,J,K,L,N,IPICKA,IPICKB,NQTOT, &
       UNNA(1:2,0:NNMAX),UNNB(1:2,0:NNMAX)
  DOUBLE PRECISION :: PSUMS(0:NATOMS+1),TWEIGHT,R,DPRAND,POTEL,&
       XMIN(3*NATOMS),EMIN,WMAX,GAIN,TEMP
  !
  ! Energy of COORDS from last quench. Common block in QUENCH. 
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  XMIN(1:3*NATOMS)=COORDS(1:3*NATOMS, NP)
  EMIN=POTEL
  !
  ! Make gloabal NNBRS and NNHIST consistent with COORDS
  CALL ANAL_NNBROS(COORDS(1:3*NATOMS,NP),.FALSE.,.FALSE.,.FALSE.)
  ! Update bond weights consistent with current ANNHIST
  CALL UPDATE_NNBOND_WEIGHTS()
  ! Calculate total sum of bond weights
  CALL CALC_TWEIGHT(TWEIGHT)  
  WMAX = TWEIGHT
  !
  WRITE(MYUNIT,'(A, F16.8,A,F10.5)') &
       'homoref_aux> Initial E=',EMIN,' W=',WMAX
  !
  TEMP = HOMOREF_AUX_TEMP
  CALL CALC_ALL_PSUMS(TEMP, PSUMS)
  !
  SWAPS: DO N=1,HOMOREF_AUX_NSWAPS
     !
     ! Pick type A atom
     R=DPRAND()*PSUMS(0)
     IPICKA=1
     DO WHILE (R > PSUMS(IPICKA))
        IPICKA = IPICKA + 1
     ENDDO
     !
     ! Pick type B atom
     R=DPRAND()*PSUMS(NATOMS+1)
     IPICKB=NTYPEA+1
     DO WHILE (R > PSUMS(IPICKB))
        IPICKB = IPICKB + 1
     ENDDO
     !
     ! seek best swap gain... appears ineffective
     !CALL FIND_BEST_SWAP(1,NATOMS,IPICKA,IPICKB,GAIN)
     !IF(GAIN <= 0.0D0) THEN
     !   EXIT
     !ENDIF
     !
     ! Sanity check:
     IF(IPICKA>NTYPEA .OR. IPICKB<=NTYPEA .OR. IPICKB>NATOMS) THEN
        WRITE(MYUNIT,'(A)') "homoref_aux> bad IPICKA and IPICKB:"
        WRITE(MYUNIT,*) IPICKA, IPICKB
        STOP
     ENDIF
     !
     !----------------------------------------------------------
     ! There seems to be a bug in the way I update NNBRO swap...
     !---------------------------------------------------------- 
     ! CALL REM_INTERSECTION(NNBRS(IPICKA,1:2,0:NNMAX), &
     !      NNBRS(IPICKB,1:2,0:NNMAX), UNNA(1:2,0:NNMAX), &
     !      UNNB(1:2,0:NNMAX))
     ! !
     ! DO I=1,2 ! loop over A- and B-type NN lists of IPICKA and IPICKB
     !    !
     !    ! Update NNBRS of IPICKA that are not neighbours of IPICKB
     !    NNA: DO J=1,UNNA(I,0) ! loop over uncommon NNs of IPICKA
     !       K=UNNA(I,J) ! Actual atom index
     !       SUCCESS=.FALSE.
     !       IDA: DO L=1,NNBRS(K,1,0) ! scan for IPICKA
     !          IF(NNBRS(K,1,L) == IPICKA) THEN
     !             ! Remove IPICKA from A-list
     !             NNBRS(K,1,L) = NNBRS(K,1,NNBRS(K,1,0))
     !             NNBRS(K,1,NNBRS(K,1,0)) = 0
     !             NNBRS(K,1,0) = NNBRS(K,1,0) - 1
     !             SUCCESS = .TRUE.
     !             EXIT IDA
     !          ENDIF
     !       ENDDO IDA
     !       IF(.NOT. SUCCESS) THEN
     !          WRITE(MYUNIT,'(A)') "homoref_aux> ERROR1:"
     !          WRITE(MYUNIT,*) IPICKA,IPICKB,I,J,K
     !          STOP
     !       ENDIF
     !       ! Add IPICKB to B-list of K, but NOT if K = IPICKA !!!
     !       IF(K/=IPICKA) THEN
     !          NNBRS(K,2,0) = NNBRS(K,2,0) + 1
     !          NNBRS(K,2,NNBRS(K,2,0)) = IPICKB
     !       ELSE
     !          WRITE(MYUNIT,'(A)') 'homoref_aux> ERROR1a'
     !          STOP
     !       END IF
     !    ENDDO NNA
     !    !
     !    ! Update NNBRS of IPICKB that are not neighbours of IPICKA
     !    NNB: DO J=1,UNNB(I,0) ! loop over uncommon NNs of IPICKB
     !       K=UNNB(I,J) ! actual atom index
     !       SUCCESS=.FALSE.
     !       IDB: DO L=1,NNBRS(K,2,0)
     !          IF(NNBRS(K,2,L) == IPICKB) THEN
     !             ! Remove IPICKB from B-lists
     !             NNBRS(K,2,L) = NNBRS(K,2,NNBRS(K,2,0))
     !             NNBRS(K,2,NNBRS(K,2,0)) = 0
     !             NNBRS(K,2,0) = NNBRS(K,2,0) - 1
     !             SUCCESS = .TRUE.
     !             EXIT IDB
     !          ENDIF
     !       ENDDO IDB
     !       IF(.NOT. SUCCESS) THEN
     !          WRITE(MYUNIT,'(A)') "homoref_aux> ERROR2"
     !          WRITE(MYUNIT,*) IPICKA,IPICKB,I,J,K
     !          STOP
     !       ENDIF
     !       ! Add IPICKA to A-lists
     !       IF(K/=IPICKB) THEN
     !          NNBRS(K,1,0) = NNBRS(K,1,0) + 1
     !          NNBRS(K,1,NNBRS(K,1,0)) = IPICKA
     !       ELSE
     !          WRITE(MYUNIT,'(A)') 'homoref_aux> ERROR2a'
     !          STOP
     !       ENDIF
     !    ENDDO NNB
     !    !
     ! ENDDO
     !
     ! Swap coords of IPICKA and IPICKB
     J=3*(IPICKA-1)
     K=3*(IPICKB-1)
     DO L=1,3
        R=COORDS(J+L,NP)
        COORDS(J+L,NP) = COORDS(K+L,NP)
        COORDS(K+L,NP) = R
     ENDDO
     ! Swap NNBRS of IPICKA and IPICKB
     ! DO J=1,2
     !    DO K=0,NNMAX
     !       L=NNBRS(IPICKA,J,K)
     !       NNBRS(IPICKA,J,K) = NNBRS(IPICKB,J,K)
     !       NNBRS(IPICKB,J,K) = L
     !    ENDDO
     ! ENDDO
     !
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
     !
     ! Could be useful for testing
     !CALL GSAVEIT(POTEL,COORDS(1:3*NATOMS,NP),NP)
     !
     IF(POTEL<EMIN) THEN
        EMIN=POTEL
        XMIN(1:3*NATOMS)=COORDS(1:3*NATOMS, NP)
        WRITE(MYUNIT,'(A,F10.5,A,I5)') 'homoref_aux> New EMIN:', &
             EMIN,' on swap', N
     ENDIF
     !
     ! Make global NNBRS and NNHIST consistent with COORDS
     ! Could devise a more efficient UPDATE scheme, as opposed to
     ! calling ANAL_NBROS every time.
     CALL ANAL_NNBROS(COORDS(1:3*NATOMS,NP),.FALSE.,.FALSE.,.FALSE.)
     !
     TEMP = TEMP*HOMOREF_AUX_FACTOR
     CALL CALC_ALL_PSUMS(TEMP, PSUMS)
     CALL CALC_TWEIGHT(TWEIGHT)
     !
     IF(TWEIGHT>WMAX) THEN
        WMAX=TWEIGHT
        WRITE(MYUNIT,'(A,F10.5,A,I5)') 'homoref_aux> New WMAX:', &
             WMAX,' on swap', N
     ENDIF
     !
  ENDDO SWAPS
  !
  POTEL = EMIN
  COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
  !
  WRITE(MYUNIT,'(A, F16.8,A,F10.5,A,I5,A)') &
       'homoref_aux> Final E=',EMIN,' W=',WMAX,&
       ' after', HOMOREF_AUX_NSWAPS,' swaps'
  !
  RETURN
  !
END SUBROUTINE HOMOREF_AUX
!
!====================================================================
!
SUBROUTINE UPDATE_NNBOND_WEIGHTS() ! N-species
  !
  USE COMMONS, ONLY : NSPECIES, MYUNIT,NATOMS, HOMOREFTEST
  USE HOMOREFMOD, ONLY : ANNHIST_MEAN, NNBOND_WEIGHTS,&
       NBONDS, NBONDS_AVE
  !
  IMPLICIT NONE
  !
  INTEGER :: T1,T2,DELTA,I
  DOUBLE PRECISION :: DUMMY, BONDCOUNT3
  !
  NBONDS_AVE = 0.0
  DO T1=1,NSPECIES(0)
     DUMMY=0.0
     DO T2 = 1,NSPECIES(0)
        DUMMY = DUMMY + ANNHIST_MEAN(T1,T2)             
     ENDDO
     NBONDS_AVE = NBONDS_AVE + DBLE(NSPECIES(T1))*DUMMY
  ENDDO
  NBONDS_AVE = 0.5D0*NBONDS_AVE
  !
  !
  ! !
  ! BONDCOUNT3=0.0D0
  ! DO T1=1,NSPECIES(0)
  !    DUMMY=0.0D0
  !    DO T2=1,T1
  !       IF(T1==T2) THEN
  !          DELTA=1
  !       ELSE
  !          DELTA=0
  !       ENDIF
  !       DUMMY = DUMMY + ANNHIST_MEAN(T1,T2)/DBLE(1+DELTA)
  !    ENDDO
  !    BONDCOUNT3 = BONDCOUNT3 + DBLE(NSPECIES(T1))*DUMMY
  ! ENDDO
  ! WRITE(MYUNIT,'(A,F10.5)') 'upd_nnb_weights> BONDCOUNT3=',BONDCOUNT3
  !
  IF(HOMOREFTEST) THEN
     WRITE(MYUNIT,'(A,F8.1)') 'upd_nnb_weights> NBONDS_AVE=',NBONDS_AVE
     WRITE(MYUNIT,'(A,I6)') 'upd_nnb_weights> 2xNBONDS=',NBONDS(1)
  ENDIF
  !
  DO T1=1,NSPECIES(0)
     !
     WRITE(MYUNIT,'(A)', ADVANCE='NO') 'upd_nnb_weights>'
     !
     DO T2 = 1,NSPECIES(0)
        !
        IF(T1==T2) THEN
           DELTA=1
        ELSE
           DELTA=0
        ENDIF
        !
        IF(NBONDS(1) > 0 .AND. NSPECIES(T2) > 1) THEN
           NNBOND_WEIGHTS(T1,T2) = DBLE(NATOMS*(NATOMS-1))* &
                ANNHIST_MEAN(T1,T2) / &
                (NBONDS_AVE*DBLE(2*(NSPECIES(T2)-DELTA))) - 1.0D0
        ELSE
           NNBOND_WEIGHTS(T1,T2) = 0.0D0
        ENDIF
        !
        WRITE(MYUNIT,'(F8.4)', ADVANCE='NO') NNBOND_WEIGHTS(T1,T2)
        !
     ENDDO
     !
     WRITE(MYUNIT,'(A)', ADVANCE='YES') ''
     !
  ENDDO
  !
END SUBROUTINE UPDATE_NNBOND_WEIGHTS
!
!====================================================================
!
SUBROUTINE CALC_TWEIGHT(TWEIGHT) ! N-species
  !
  USE COMMONS, ONLY : NSPECIES
  USE HOMOREFMOD, ONLY : NNBOND_WEIGHTS,NNBRS
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(OUT) :: TWEIGHT
  INTEGER :: T1, T2, I, J
  !
  TWEIGHT = 0.0D0
  I=0
  DO T1=1,NSPECIES(0) 
     DO J=I+1,I+NSPECIES(T1) ! span T1 atoms
        DO T2=1,NSPECIES(0) 
           TWEIGHT = TWEIGHT + &
                NNBRS(J,T2,0)*NNBOND_WEIGHTS(T1,T2)
        ENDDO
     ENDDO
     I=I+NSPECIES(T1)
  ENDDO
  TWEIGHT = 0.5D0*TWEIGHT
  !
  RETURN
  !
END SUBROUTINE CALC_TWEIGHT
!
!====================================================================
!
! SUBROUTINE FIND_BEST_SWAP(LI,UI,IA,IB,GMAX) ! 2-species
!   !
!   USE COMMONS, ONLY : NATOMS, NTYPEA, MYUNIT
!   USE HOMOREFMOD, ONLY : NNBRS, NNMAX, NNBOND_WEIGHTS
!   !
!   IMPLICIT NONE
!   !
!   INTEGER, INTENT(IN) :: LI, UI
!   INTEGER, INTENT(OUT) :: IA,IB
!   DOUBLE PRECISION, INTENT(OUT) :: GMAX 
!   !
!   INTEGER :: I,J,K,UNNI(2,0:NNMAX),UNNJ(2,0:NNMAX)
!   DOUBLE PRECISION :: G
!   !
!   GMAX=-9.9D9
!   !
!   DO I=LI,NTYPEA
!      DO J=NTYPEA+1, NATOMS
!         CALL REM_INTERSECTION(NNBRS(I,1:2,0:NNMAX), &
!              NNBRS(J,1:2,0:NNMAX),UNNI(1:2,0:NNMAX),&
!              UNNJ(1:2,0:NNMAX))
!         G = DBLE(UNNI(1,0)-UNNJ(1,0))* &
!              (NNBOND_WEIGHTS(1,2)-NNBOND_WEIGHTS(1,1)) + &
!              DBLE(UNNI(2,0)-UNNJ(2,0))* &
!              (NNBOND_WEIGHTS(2,2)-NNBOND_WEIGHTS(2,1))
!         IF(GMAX < G) THEN
!            GMAX=G
!            IA=I
!            IB=J
!         ENDIF
!      ENDDO
!   ENDDO
!   !
!   WRITE(MYUNIT,'(A,F10.5)') 'find_best_swap> GMAX=', GMAX
!   !
!   RETURN
!   !
! END SUBROUTINE FIND_BEST_SWAP
!
!====================================================================
!
SUBROUTINE CALC_ALL_PSUMS(TEMP, PSUMS) ! 2-species
  !
  USE COMMONS, ONLY : NATOMS, NTYPEA
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: TEMP
  DOUBLE PRECISION, INTENT(OUT) :: PSUMS(0:NATOMS+1)
  DOUBLE PRECISION :: EXPDW
  INTEGER :: I
  !
  PSUMS(0:NATOMS+1) = 0.0D0
  !
  CALL CALC_ONE_EXPDW(1,TEMP,EXPDW)
  PSUMS(1) = EXPDW
  PSUMS(0) = EXPDW
  DO I=2,NTYPEA
     CALL CALC_ONE_EXPDW(I,TEMP,EXPDW)
     PSUMS(I) = PSUMS(I-1)+EXPDW
     PSUMS(0) = PSUMS(0) + EXPDW     
  ENDDO
  !
  CALL CALC_ONE_EXPDW(NTYPEA+1,TEMP,EXPDW)
  PSUMS(NTYPEA+1) = EXPDW
  PSUMS(NATOMS+1) = EXPDW
  DO I=NTYPEA+2,NATOMS
     CALL CALC_ONE_EXPDW(I,TEMP,EXPDW)
     PSUMS(I) = PSUMS(I-1)+EXPDW
     PSUMS(NATOMS+1) = PSUMS(NATOMS+1) + EXPDW
  ENDDO
  !
  RETURN
  !
END SUBROUTINE CALC_ALL_PSUMS
!
!====================================================================
!
SUBROUTINE CALC_ONE_EXPDW(I,TEMP,EXPDW) ! binary only...
  !
  USE HOMOREFMOD, ONLY : NNBRS, NNBRS, NNBOND_WEIGHTS
  USE COMMONS, ONLY : NATOMS, NTYPEA, MYUNIT
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: I
  DOUBLE PRECISION, INTENT(IN) :: TEMP
  DOUBLE PRECISION, INTENT(OUT) :: EXPDW
  !
  INTEGER :: T1, T2, T3
  !
  IF(I < 1 .OR. I > NATOMS) THEN
     WRITE(MYUNIT,'(A)') "calc_one_weight> Atom index out of bounds!"
     STOP
  ENDIF
  !
  IF(I>NTYPEA) THEN
     T1=2 ! INITIAL
     T2=1 ! FINAL
  ELSE
     T1=1 ! INITIAL
     T2=2 ! FINAL
  ENDIF
  !
  EXPDW = 0.0
  ! 
  DO T3=1,2
     EXPDW = EXPDW + NNBRS(I,T3,0)* &
          (NNBOND_WEIGHTS(T2,T3)-NNBOND_WEIGHTS(T1,T3))
  ENDDO
  !
  EXPDW = DEXP(EXPDW/TEMP)
  !
  RETURN
  !
END SUBROUTINE CALC_ONE_EXPDW
!
!====================================================================
!
! SUBROUTINE REM_INTERSECTION(V1IN,V2IN,V1OUT,V2OUT) ! 2-species
!   !
!   USE HOMOREFMOD, ONLY : NNMAX
!   !
!   IMPLICIT NONE
!   !
!   INTEGER, INTENT(IN) :: V1IN(1:2,0:NNMAX), V2IN(1:2,0:NNMAX)
!   INTEGER, INTENT(OUT) :: V1OUT(1:2,0:NNMAX), V2OUT(1:2,0:NNMAX)
!   !
!   LOGICAL :: KEEP
!   INTEGER :: I,J,K,ISEC(0:NNMAX)
!   !
!   DO K=1,2
!      ! First find intersection
!      ISEC(0:NNMAX) = 0
!      DO I=1,V1IN(K,0)
!         DO J=1,V2IN(K,0)
!            IF(V1IN(K,I)==V2IN(K,J)) THEN
!               ISEC(0) = ISEC(0)+1
!               ISEC(ISEC(0)) = V1IN(K,I)
!            ENDIF
!         ENDDO
!      ENDDO
!      !
!      ! Remove intersection from V1IN and store in V1OUT
!      V1OUT(K,0:NNMAX) = 0
!      DO I=1,V1IN(K,0)
!         KEEP=.TRUE.
!         DO J=1,ISEC(0)
!            IF(V1IN(K,I)==ISEC(J)) THEN
!               KEEP = .FALSE.
!            ENDIF
!         ENDDO
!         IF(KEEP) THEN
!            V1OUT(K,0) = V1OUT(K,0)+1
!            V1OUT(K,V1OUT(K,0)) = V1IN(K,I)
!         ENDIF
!      ENDDO
!      !
!      ! Remove intersection from V2IN and store in V2OUT
!      V2OUT(K,0:NNMAX) = 0
!      DO I=1,V2IN(K,0)
!         KEEP=.TRUE.
!         DO J=1,ISEC(0)
!            IF(V2IN(K,I)==ISEC(J)) THEN
!               KEEP = .FALSE.
!            ENDIF
!         ENDDO
!         IF(KEEP) THEN
!            V2OUT(K,0) = V2OUT(K,0)+1
!            V2OUT(K,V2OUT(K,0)) = V2IN(K,I)
!         ENDIF
!      ENDDO
!      !
!   ENDDO
!   !
!   RETURN
!   !
! END SUBROUTINE REM_INTERSECTION
!
!====================================================================
! What follows is a set of routines for "random permutation" moves
! The total number of pair swaps is randomly chosen from the 
! interval [0.5*MIN(NA,NB), MIN(NA,NB)].
! The routine SWAP_COORDS is in file homoref.f90
!====================================================================
!
SUBROUTINE RANDPERM(NP)
  !
  ! ds656> 30/10/2013
  ! Swap coordinates of randomly-picked unlike atom pairs in a
  ! binary system. 
  !
  USE COMMONS, ONLY : NATOMS, NTYPEA, MYUNIT
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  INTEGER :: N, NMAX,LI, UI, I1,J1,I2,J2
  DOUBLE PRECISION :: DPRAND, X
  !
  ! randomly pick the total swap count
  NMAX = INT((0.5D0 + 0.5D0*DPRAND())*DBLE(MIN(NTYPEA,NATOMS-NTYPEA))) + 1
  !NMAX = NINT((0.5D0 + 0.5D0*DPRAND())*DBLE(MIN(NTYPEA,NATOMS-NTYPEA)))
  !
  DO N = 1, NMAX
     !
     ! Set lower (LI) and upper (UI) indices bracketing
     ! atom candidates for flipping. 1 <= LI < UI <= NATOMS.
     LI = N
     UI = NATOMS - N + 1
     !
     ! pick atom at random
     I1 = INT(DPRAND()*DBLE(UI-LI+1)) + 1
     !
     ! determine the type of picked atom and update bounds
     IF(I1 > NTYPEA) THEN
        J1 = UI
        J2 = LI
        UI = NTYPEA
     ELSE
        J1 = LI
        J2 = UI
        LI = NTYPEA + 1
     ENDIF
     !
     ! pick 2nd atom of different type
     I2 = INT(DPRAND()*DBLE(UI-LI+1)) + 1
     !
     ! swap and move out of bounds for next iteration
     CALL SWAP_COORDS(NP,I1,I2,.FALSE.)
     CALL SWAP_COORDS(NP,I1,J1,.FALSE.)
     CALL SWAP_COORDS(NP,I2,J2,.FALSE.)
     !
  END DO
  !
  WRITE(MYUNIT, '(A,I3,A)') "randperm> Swapped ",NMAX,' random pairs.'
  !
  RETURN
  !
END SUBROUTINE RANDPERM
