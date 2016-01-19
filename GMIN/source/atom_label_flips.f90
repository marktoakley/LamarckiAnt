SUBROUTINE FLIPSEQ(NP, ITER, TIME, BRUN, QDONE, SCREENC)
  !
  USE COMMONS, ONLY : NATOMS, NSPECIES, ATOMLISTS, INVATOMLISTS, &
       MYUNIT, COORDS, NQ, NPAR, RMS, ECONV, PRTFRQ, HIT, &
       LFLIPS_N2, LFLIPS_TEMP, LFLIPS_TFACTOR, LFLIPS_NUP, &
       LFLIPS_RESET
  !
  IMPLICIT NONE
  !
  ! Parse passed variables
  INTEGER, INTENT(IN) :: NP
  INTEGER, INTENT(INOUT):: ITER,BRUN,QDONE ! for QUENCH
  DOUBLE PRECISION, INTENT(INOUT) :: TIME,SCREENC(3*NATOMS) ! QUENCH
  !
  INTEGER :: NQTOT, ATOMLISTS_MIN(NSPECIES(0),3,0:NATOMS), &
       INVATOMLISTS_MIN(NATOMS,3),I,LIST(0:NATOMS),J,K,L,N, &
       NACC,NREJ,L_OLD,L_NEW,NSPECIES_MIN(0:NSPECIES(0))
  DOUBLE PRECISION :: DPRAND, R, POTEL, X0(1:3*NATOMS), E0, &
       XMIN(1:3*NATOMS), EMIN, TEMP
  !
  ! Energy of COORDS from last quench. Common block in QUENCH. 
  COMMON /MYPOT/ POTEL
  ! Total quench count. Commom block in MC.
  COMMON /TOT/ NQTOT
  !
  WRITE(MYUNIT,'(A)') &
       '============================================================'
  !
  IF(LFLIPS_RESET) THEN
     WRITE(MYUNIT,'(A)') 'flipseq> Resetting stoichiometry...'
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
  NSPECIES_MIN(:) = NSPECIES(:)
  ATOMLISTS_MIN(:,:,:) = ATOMLISTS(:,:,:)
  INVATOMLISTS_MIN(:,:) = INVATOMLISTS(:,:)
  XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS, NP)
  EMIN = POTEL
  TEMP = LFLIPS_TEMP
  !
  WRITE(MYUNIT,'(A,F20.10,A,F10.8)',ADVANCE='NO') &
       'flipseq> Initial E= ', POTEL,' T= ', TEMP
  DO I=1,NSPECIES(0)
     WRITE(MYUNIT,'(1X,A,I1,A,I3)',ADVANCE='NO') &
          'N_',I,'= ', NSPECIES(I)
  ENDDO
  WRITE(MYUNIT,*)
  !
  NACC = 0
  NREJ = 0
  DO N=1,LFLIPS_N2
     !
     X0(1:3*NATOMS) = COORDS(1:3*NATOMS,NP)
     E0 = POTEL
     !
     LIST(0:NATOMS) = 0
     DO J=1,NSPECIES(0)
        DO K=1,ATOMLISTS(J,1,0) ! only group 1
           LIST(0) = LIST(0) + 1
           LIST(LIST(0)) = ATOMLISTS(J,1,K) 
        ENDDO
     ENDDO
     I = INT(DPRAND()*DBLE(LIST(0))) + 1 ! Forgot +1
     I = LIST(I) ! This is the index of atom to be transmuted
     !
     L_OLD = INVATOMLISTS(I,1) ! Store old label of I
     !
     K = INT(DPRAND()*DBLE(NSPECIES(0)-1)) + 1
     L=0
     SPEC_LOOP: DO J=1,NSPECIES(0)
        IF(J==L_OLD) CYCLE
        L=L+1
        IF(L==K) THEN
           L_NEW = J
           EXIT SPEC_LOOP
        ENDIF
     ENDDO SPEC_LOOP
     !
     CALL FLIP_LABEL(I,L_NEW)
     !
     NQTOT = NQTOT + 1
     NQ(NP) = NQ(NP) + 1
     CALL QUENCH(.FALSE.,NP,ITER,TIME,BRUN,QDONE,SCREENC)
     IF (MOD(N-1,PRTFRQ).EQ.0) THEN
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
     R=DEXP( (E0-POTEL)/TEMP )
     IF(DPRAND() < R) THEN ! Accept
        NACC = NACC + 1
        IF(POTEL < EMIN - ECONV) THEN ! Store labels and coordinates
           EMIN=POTEL
           XMIN(1:3*NATOMS)=COORDS(1:3*NATOMS, NP)
           ATOMLISTS_MIN(:,:,:) = ATOMLISTS(:,:,:)
           INVATOMLISTS_MIN(:,:) = INVATOMLISTS(:,:)
           NSPECIES_MIN(:) = NSPECIES(:)
           WRITE(MYUNIT,'(A,F20.10,A,I8,A,F10.8)',ADVANCE='NO') &
                'flipseq> New EMIN= ',EMIN,' on flip', N,&
                ' at T=', TEMP
           DO K=1,NSPECIES(0)
              WRITE(MYUNIT,'(1X,A,I1,A,I3)',ADVANCE='NO') &
                   'N_',K,'= ', NSPECIES_MIN(K)
           ENDDO
           WRITE(MYUNIT,*)
           IF(HIT) RETURN
        ENDIF
     ELSE ! Reject. Revert labels and coordinates.
        NREJ = NREJ + 1
        CALL FLIP_LABEL(I,L_OLD)
        COORDS(1:3*NATOMS,NP) = X0(1:3*NATOMS)
        POTEL = E0
     ENDIF
     !
     IF(MOD(N,LFLIPS_NUP)==0) THEN
        IF(LFLIPS_TFACTOR <= 1.0D0) THEN 
           ! Simulated annealing: always decrease temperature. 
           TEMP = TEMP*LFLIPS_TFACTOR
        ELSE ! Strive for ACC/REJ ratio of 0.5.
           ! This will have to be recoded for desired ACC/REJ 
           ! ratios other than 0.5!
           IF(NACC < NREJ) THEN
              TEMP = TEMP*LFLIPS_TFACTOR
           ELSEIF(NACC > NREJ) THEN
              TEMP = TEMP/LFLIPS_TFACTOR
           ENDIF
           WRITE(MYUNIT,'(A,F7.5,A,F7.5)') &
                'flipseq> Accepted ratio ', &
                DBLE(NACC)/DBLE(NACC+NREJ),' new T= ', TEMP
           NACC = 0; NREJ = 0
        ENDIF
     ENDIF
     !
  ENDDO
  !
  NSPECIES(:) = NSPECIES_MIN(:)
  ATOMLISTS(:,:,:) = ATOMLISTS_MIN(:,:,:)
  INVATOMLISTS(:,:) = INVATOMLISTS_MIN(:,:)
  COORDS(1:3*NATOMS, NP) = XMIN(1:3*NATOMS)
  POTEL = EMIN
  !
  WRITE(MYUNIT,'(A,F20.10,A,F8.6)',ADVANCE='NO') &
       'flipseq> Final E= ', POTEL,' T= ', TEMP
  DO I=1,NSPECIES(0)
     WRITE(MYUNIT,'(1X,A,I1,A,I3)',ADVANCE='NO') &
          'N_',I,'= ', NSPECIES(I)
  ENDDO
  WRITE(MYUNIT,*)
  !
  RETURN
  !
END SUBROUTINE FLIPSEQ
