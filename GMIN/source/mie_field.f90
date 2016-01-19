SUBROUTINE MIEF_INI()
  !
  USE COMMONS, ONLY : MIEF_FILENAME, MIEF_EPS, MIEF_SIG, &
       MIEF_NSITES, MIEF_N, MIEF_M, NSPECIES, MIEF_SITES, &
       MIEF_CUTT, MIEF_RCUT, MIEF_U_RCUT, MIEF_DUDR_RCUT, MYUNIT
  !
  IMPLICIT NONE
  !
  LOGICAL :: YESNO
  INTEGER :: I, J, EOF, LUNIT, GETUNIT
  DOUBLE PRECISION :: PREF, DUMMY, REP_TERM, ATT_TERM
  !
  ALLOCATE(MIEF_SIG(NSPECIES(0)))
  ALLOCATE(MIEF_EPS(NSPECIES(0)))
  !
  ! Process the file with Mie parameters and site coordinates
  !
  INQUIRE(FILE=MIEF_FILENAME, EXIST=YESNO)
  !
  IF(YESNO) THEN
     LUNIT=GETUNIT()
     OPEN (LUNIT, FILE=MIEF_FILENAME, STATUS='OLD')
     READ(LUNIT,*) MIEF_NSITES, MIEF_N, MIEF_M, &
          (MIEF_EPS(I), I=1,NSPECIES(0)), &
          (MIEF_SIG(I), I=1,NSPECIES(0))
     !
     ALLOCATE(MIEF_SITES(MIEF_NSITES,3))
     !
     WRITE(MYUNIT,'(A,I2,A,I1,A,I3,A)') 'mief_ini> ',&
          MIEF_N,'-',MIEF_M,' Mie surface with ', &
          MIEF_NSITES,' sites.'
     WRITE(MYUNIT,'(A)', ADVANCE='NO') 'mief_ini> Epsilon value(s): '
     DO I=1,NSPECIES(0)
        WRITE(MYUNIT,'(F8.4)', ADVANCE='NO') MIEF_EPS(I)
     ENDDO
     WRITE(MYUNIT,*)
     WRITE(MYUNIT,'(A)', ADVANCE='NO') 'mief_ini> Sigma value(s): '
     DO I=1,NSPECIES(0)
        WRITE(MYUNIT,'(F8.4)', ADVANCE='NO') MIEF_SIG(I)
     ENDDO
     WRITE(MYUNIT,*)
     !
     DO I=1,MIEF_NSITES
        READ(LUNIT,*,IOSTAT=EOF) ( MIEF_SITES(I,J), J=1,3 )
        IF(EOF.LT.0) THEN
           WRITE(MYUNIT,*) "mief_ini> Premature end of file ", &
                TRIM(MIEF_FILENAME)
           STOP
        ENDIF
     ENDDO
     !
     CLOSE(LUNIT)
     !
  ELSE
     WRITE(MYUNIT,*) "mief_ini> Missing file ", TRIM(MIEF_FILENAME)
     STOP
  ENDIF
  !
  IF (MIEF_CUTT) THEN
     ALLOCATE(MIEF_U_RCUT(NSPECIES(0)))
     ALLOCATE(MIEF_DUDR_RCUT(NSPECIES(0)))
     PREF = DBLE(MIEF_N)/DBLE(MIEF_N-MIEF_M)*&
          (DBLE(MIEF_N)/DBLE(MIEF_M))**&
          (DBLE(MIEF_M)/DBLE(MIEF_N-MIEF_M))
     DO I=1,NSPECIES(0)
        DUMMY = MIEF_SIG(I)/MIEF_RCUT
        REP_TERM=DUMMY**MIEF_N
        ATT_TERM=DUMMY**MIEF_M
        ! ds656> Had a fucking mistake here for ages!!!!! 
        !MIEF_U_RCUT=PREF*MIEF_EPS(I)*(REP_TERM-ATT_TERM)
        !MIEF_DUDR_RCUT=PREF*MIEF_EPS(I)*( -DBLE(MIEF_N)*REP_TERM + &
        !     DBLE(MIEF_M)*ATT_TERM ) / MIEF_RCUT
        ! ds656> Spot the difference... why was the compiler quiet?!?
        MIEF_U_RCUT(I)=PREF*MIEF_EPS(I)*(REP_TERM-ATT_TERM)
        MIEF_DUDR_RCUT(I)=PREF*MIEF_EPS(I)*( -DBLE(MIEF_N)*REP_TERM + &
             DBLE(MIEF_M)*ATT_TERM ) / MIEF_RCUT
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MIEF_INI
!
SUBROUTINE MIEF(X,GRAD,EREAL,GRADT, STRESST)
  !
  USE COMMONS, ONLY : MIEF_N, MIEF_M, MIEF_NSITES, MIEF_SITES, &
       MIEF_EPS, MIEF_SIG, ATOMLISTS, INVATOMLISTS, VT, NATOMS, &
       NSPECIES, MIEF_PBCT, MIEF_BOX, MIEF_CUTT, MIEF_RCUT, &
       MIEF_U_RCUT, MIEF_DUDR_RCUT, STRESS
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(INOUT) :: GRAD(3*NATOMS)
  DOUBLE PRECISION, INTENT(INOUT) :: EREAL
  LOGICAL, INTENT(IN) :: GRADT, STRESST
  !
  INTEGER :: I, J, K, L, I3, ITYPE
  DOUBLE PRECISION :: PREF, DX(3),DIST2,IDIST,SIG,EPS,&
       DUMMY, ATT_TERM, REP_TERM, DIST
  !
  PREF = DBLE(MIEF_N)/DBLE(MIEF_N-MIEF_M)*&
       (DBLE(MIEF_N)/DBLE(MIEF_M))**&
       (DBLE(MIEF_M)/DBLE(MIEF_N-MIEF_M))
  !
  DO I=1,NATOMS ! Loop over ALL atoms
     !
     I3=3*(I-1)
     ITYPE=INVATOMLISTS(I,1)
     EPS=MIEF_EPS(ITYPE)*PREF
     SIG=MIEF_SIG(ITYPE)
     !
     DO J=1,MIEF_NSITES ! Loope over Mie sites
        !
        DIST2=0.0
        DO K=1,3
           DX(K) = X(I3+K) - MIEF_SITES(J,K)
           IF (MIEF_PBCT) DX(K) = DX(K) - MIEF_BOX(K)*&
                NINT(DX(K)/MIEF_BOX(K))
           DIST2 = DIST2 + DX(K)*DX(K)
        ENDDO
        DIST = DSQRT(DIST2)
        !
        IF(DIST < MIEF_RCUT) THEN
           !
           IDIST = SIG / DIST
           REP_TERM=IDIST**MIEF_N
           ATT_TERM=IDIST**MIEF_M
           !
           DUMMY = EPS*(REP_TERM-ATT_TERM)
           IF(MIEF_CUTT) THEN
              DUMMY = DUMMY - MIEF_U_RCUT(ITYPE) - &
                   (DIST-MIEF_RCUT)*MIEF_DUDR_RCUT(ITYPE)
           ENDIF
           VT(I) = VT(I) + DUMMY
           EREAL = EREAL + DUMMY
           !
           IF(GRADT) THEN
              ! Compute 1/R*dU/dR
              DUMMY = EPS*( -DBLE(MIEF_N)*REP_TERM + &
                   DBLE(MIEF_M)*ATT_TERM ) / DIST2
              IF(MIEF_CUTT) DUMMY = DUMMY - MIEF_DUDR_RCUT(ITYPE)/DIST
              DO K=1,3
                 GRAD(I3+K) = GRAD(I3+K) + DUMMY*DX(K)
                 IF(STRESST) THEN ! Accumulate local stresses
                    DO L=K,3
                       STRESS(I,K,L) = STRESS(I,K,L) + &
                            DUMMY*DX(K)*DX(L)
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF
           !
        ENDIF ! Cutoff check
        !
     ENDDO
  ENDDO
  !
  IF(GRADT) THEN
     DO I=1,NSPECIES(0) ! Atom types
        DO J=1,ATOMLISTS(I,2,0) ! group 2 only
           K=3*(ATOMLISTS(I,2,J)-1)
           GRAD(K+1) = 0.0D0 ! Reset to zero
           GRAD(K+2) = 0.0D0
           GRAD(K+3) = 0.0D0
        ENDDO
     END DO
  ENDIF
  !
  IF(STRESST) THEN
     DO I=1,NATOMS
        DO K=1,3
           DO L=K,3
              STRESS(0,K,L) = STRESS(0,K,L) + STRESS(I,K,L)
              STRESS(I,K,L) = STRESS(I,L,K) ! Impose symmetry
           ENDDO
        ENDDO
     ENDDO
     DO K=1,3
        DO L=K+1,3
           STRESS(0,K,L) = STRESS(0,L,K) ! Impose symmetry
        ENDDO
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MIEF
