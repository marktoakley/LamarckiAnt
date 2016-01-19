
!  GMIN: A program for finding global minima
!  CopyrIGht (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

MODULE VGWSP
     IMPLICIT NONE
     INTEGER, SAVE :: N_ATOM, N3_ATOM, CONPOT, NSTEP
     INTEGER, PARAMETER :: N_DIM=3, NGAUSS=3
     DOUBLE PRECISION :: E_ZERO,BL,ATOL
     DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
     DOUBLE PRECISION, PARAMETER :: ATOMICMASS=0.020614788876D0
     DOUBLE PRECISION, SAVE:: SIGMA, EPSILON, LJG(NGAUSS), LJK(NGAUSS), LJG2(NGAUSS)
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: MASSARRY(:), TAU(:)
     DOUBLE PRECISION, SAVE ::  TAUMAX, CGK(3), CGG(3), CGG2(3), CGCONST
!    EXTERNAL DLSODE
!    EXTERNAL JAC
END MODULE VGWSP

MODULE VGW
     IMPLICIT NONE
     INTEGER, SAVE :: N_ATOM, N3_ATOM, CONPOT
     DOUBLE PRECISION, SAVE :: LJK(6), LJG(6), LJG2(6), E_ZERO, SIGMA, SIG2, U, BL
     INTEGER, PARAMETER :: N_DIM=3, NUM_G=6
     DOUBLE PRECISION, PARAMETER :: LJC1=31279960.65933084D0, LJC2=1668963.963961670D0
     DOUBLE PRECISION, PARAMETER :: LJC3=91092.34069670191D0, LJC4=3354.805129558428D0 
     DOUBLE PRECISION, PARAMETER :: LJC5=-8.46844309983970D0, LJC6=-0.38418467585210D0 
     DOUBLE PRECISION, PARAMETER :: LJA1=35.14249661727566D0, LJA2=21.73050942017830D0 
     DOUBLE PRECISION, PARAMETER :: LJA3=13.25329843520143D0, LJA4=7.609820703336350D0 
     DOUBLE PRECISION, PARAMETER :: LJA5=1.671802581756990D0, LJA6=0.502618140953350D0
     DOUBLE PRECISION, PARAMETER :: ATOMICMASS=0.020614788876D0
     DOUBLE PRECISION, SAVE :: TAUMAX, CGK(3), CGG(3), CGG2(3), CGCONST, EPSILON
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: MASSARRY(:)
!    EXTERNAL DLSODE
!    EXTERNAL JAC
END MODULE VGW
  
SUBROUTINE CLEANUP_VGWSP
  USE VGWSP
  DEALLOCATE(MASSARRY,TAU)
END SUBROUTINE

SUBROUTINE CLEANUP_VGW
  USE VGW
  DEALLOCATE(MASSARRY)
END SUBROUTINE

SUBROUTINE JAC
END SUBROUTINE

SUBROUTINE INITIALIZE_VGW(N, LJS, LJE, TMAX, CPSIGMA, CP)
  USE VGW
  IMPLICIT NONE
  INTEGER :: I, N, NSAVE, CP
  DOUBLE PRECISION :: LJS, LJE, TMAX, RADIUS, LAM, COEF, CPA, CPSIGMA, CPS, CPD
  DOUBLE PRECISION, PARAMETER :: C1=1.83365D0, C2=0.238678D0, C3=1.0029D0

  N_ATOM=N  
  N3_ATOM=3*N_ATOM
  EPSILON=LJE
  BL=-1D0               ! BL=BOXLENGTH*SIGMA

  ALLOCATE(MASSARRY(3*N_ATOM))

  OPEN(UNIT=7,FILE='vgwdata',STATUS='OLD')
    DO I=0,N_ATOM-1
      READ(7,*) MASSARRY(3*I+1)
      MASSARRY(3*I+1)=1/(ATOMICMASS*MASSARRY(3*I+1))
      MASSARRY(3*I+2)=MASSARRY(3*I+1)
      MASSARRY(3*I+3)=MASSARRY(3*I+1)
    ENDDO
  CLOSE(7)

  SIGMA=LJS
  SIG2=SIGMA*SIGMA  

  CONPOT=CP

  CPA = 1000
  CPD = 0.5
  CPS = CPSIGMA*0.00001D0

  LJK(1)=LJC1*LJE
  LJK(2)=LJC2*LJE
  LJK(3)=LJC3*LJE 
  LJK(4)=LJC4*LJE
  LJK(5)=LJC5*LJE
  LJK(6)=LJC6*LJE

  LJG(1)=LJA1/SIG2
  LJG(2)=LJA2/SIG2
  LJG(3)=LJA3/SIG2
  LJG(4)=LJA4/SIG2
  LJG(5)=LJA5/SIG2
  LJG(6)=LJA6/SIG2

  LJG2(1)=LJG(1)*LJG(1)
  LJG2(2)=LJG(2)*LJG(2)
  LJG2(3)=LJG(3)*LJG(3)        
  LJG2(4)=LJG(4)*LJG(4)
  LJG2(5)=LJG(5)*LJG(5)
  LJG2(6)=LJG(6)*LJG(6)
   
  COEF=-CPA/(CPD**2)
  CGK(1)= COEF*(1/(1+CPD))
  CGK(2)=COEF*(1/(1-CPD))
  CGK(3)=-2*COEF

  CGG(1)=CPS*(1+CPD)
  CGG(2)=CPS*(1-CPD)
  CGG(3)=CPS

  CGG2(1)=CGG(1)**2
  CGG2(2)=CGG(2)**2
  CGG2(3)=CGG(3)**2

  CGCONST = 2*CPA/(1-CPD**2)

  TAUMAX=TMAX

  IF(TAUMAX.LE.1.0D-6) THEN
    PRINT *, "TAUMAX TOO SMALL, CHANGING TO 1.0"
    TAUMAX=1.0D0
  ENDIF
  
 ! LAM = (0.42485D0)/SQRT(MI)
 ! CPS = CPS/(C1*LAM**2 + C2*LAM + C3)

END SUBROUTINE INITIALIZE_VGW

SUBROUTINE INITIALIZE_VGWSP(N,SIG,EPS,TMAX,CPSIGMA,CP,AT)
  USE VGWSP
  IMPLICIT NONE
  DOUBLE PRECISION :: SIG2,SIG,EPS,TMAX,TMIN,TSTEP,CPA,CPSIGMA,CPS,CPD,COEF,MI,LAM,AT
  DOUBLE PRECISION, PARAMETER :: C1=1.83365,C2=0.238678,C3=1.0029
  INTEGER :: I,N,NSAVE,CP
 
  N_ATOM=N
  SIGMA = SIG
  EPSILON=EPS
  N3_ATOM = 3*N_ATOM
  SIG2=SIGMA*SIGMA
  BL=-1D0
  ATOL=AT
  ALLOCATE(MASSARRY(3*N_ATOM))

  OPEN(UNIT=7,FILE='vgwdata',STATUS='OLD')
    DO I=0,N_ATOM-1
      READ(7,*) MASSARRY(3*I+1)
      MASSARRY(3*I+1)=1/(ATOMICMASS*MASSARRY(3*I+1))
      MASSARRY(3*I+2)=MASSARRY(3*I+1)
      MASSARRY(3*I+3)=MASSARRY(3*I+1)
    ENDDO
  CLOSE(7)

  CPA=1000.0D0
  CPD=0.5D0
  CPS=CPSIGMA*0.00001D0

  IF(CPS.EQ.0) THEN
    CPS=0.00001
  ENDIF
  
  LJK(1)=460.0D0*EPS*4.0D0        !**** LJ GAUSSIAN PARAMETERS ****
  LJK(2)=-0.37D0*EPS*4.0D0
  LJK(3)=-5.8D0*EPS*4.0D0
  LJG(1)=6.65D0/SIG2;
  LJG(2)=0.79D0/SIG2
  LJG(3)=2.6D0/SIG2
  LJG2(1)=LJG(1)*LJG(1)
  LJG2(2)=LJG(2)*LJG(2) 
  LJG2(3)=LJG(3)*LJG(3)

  COEF=-CPA/(CPD**2)              !**** CONSTRAINING POTENTIAL GUASSIAN PARAMETERS ****
  CGK(1)= COEF*(1/(1+CPD))
  CGK(2)=COEF*(1/(1-CPD))
  CGK(3)=-2*COEF

  CGG(1)=CPS*(1+CPD)
  CGG(2)=CPS*(1-CPD)
  CGG(3)=CPS

  CGG2(1)=CGG(1)**2
  CGG2(2)=CGG(2)**2
  CGG2(3)=CGG(3)**2

  CGCONST = 2*CPA/(1-CPD**2)

  TAUMAX=TMAX

  IF(TAUMAX.LE.1.0D-6) THEN
    PRINT *, "TAUMAX TOO SMALL, CHANGING TO 1.0"
    TAUMAX=1.0D0
  ENDIF

  CONPOT=CP
  
 ! LAM = (0.42485D0)/SQRT(MI)
 ! CPS=CPS/(C1*LAM**2 + C2*LAM + C3)

  TMIN=0.001D0
  TSTEP=0.1D0

  IF(TAUMAX.LE.TMIN) THEN
    TAUMAX = 10*TMIN
    PRINT *, "TAU MAX TOO SMALL, CHANGING TO 0.01"
  ENDIF

  NSTEP=((TAUMAX-TMIN)/TSTEP)+1
 
  ALLOCATE(TAU(NSTEP))
 
  DO I=0,NSTEP-1
    TAU(I+1)=0.5D0*(TMIN+I*TSTEP)                    ! ARRAY OF IMAGINARY TIME POINTS (FOR PROPAGATION)
  ENDDO

END SUBROUTINE INITIALIZE_VGWSP

SUBROUTINE VGWQUENCH(QCNFG, ENRG, CFLAG)
  USE VGW
  IMPLICIT NONE
  INTEGER :: ITASK, ISTATE, IOPT, ITOL, LRW, NEQ, I, J, K
  DOUBLE PRECISION :: QCNFG(N3_ATOM),QN(N3_ATOM), QC(3), ENRG, ULJ
  DOUBLE PRECISION :: T,TOUT,LNZ(2,2), LOGZ
  INTEGER, ALLOCATABLE, DIMENSION (:), SAVE :: IWORK
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: RWORK, Y
  DOUBLE PRECISION, PARAMETER :: TAUI=0.000001D0, RTOL=0.0D0
  DOUBLE PRECISION, PARAMETER :: ATOL=3.0D-8
  INTEGER, PARAMETER :: MF=10, LIW=20
  LOGICAL CFLAG
  EXTERNAL RHS
  EXTERNAL JAC

  CFLAG=.TRUE.

  NEQ=1+N3_ATOM+N3_ATOM*(N3_ATOM+1)/2;
  LRW=20+16*NEQ
  LNZ=0.0D0
  ITOL=1

  ALLOCATE(RWORK(LRW), IWORK(LIW), Y(NEQ))

  ITASK=1
  IOPT=1
  ISTATE=1
  IWORK(5)=4
  IWORK(6)=100000         
  IWORK(7)=0
  IWORK(8)=0
  IWORK(9)=0
  RWORK(5)=0.0D0
  RWORK(6)=0.0D0
  RWORK(7)=0.0D0
  RWORK(8)=0.0D0
  RWORK(9)=0.0D0
  RWORK(10)=0.0D0
 
  E_ZERO=ENRG
 
  DO I=1,N_DIM
    QC(I)=0.0D0
  ENDDO

  DO I=1,N3_ATOM
    QN(I)=QCNFG(I)*SIGMA
  ENDDO
  
  CALL GAUSSENERGYPB(QN, N_ATOM, SIGMA, EPSILON, BL, ULJ)
 
  DO I=0,N_ATOM-1
    DO J=1,3
      QC(J)=QC(J)+QN(3*I+J)                            ! CENTER OF MASS VECTOR
    ENDDO
  ENDDO

  DO I=1,3
    QC(I)=QC(I)/N_ATOM                                 ! CENTER OF MASS
  ENDDO

  DO I=0,N_ATOM-1
    DO J=1,3
      QN(3*I+J)=QN(3*I+J)-QC(J)                        ! SHIFT CENTER OF MASS TO ORIGIN
    ENDDO
  ENDDO

  T=TAUI                                               ! T (CURRENT TAU PARAMTER) = INITITAL TAU
  Y(1)=-TAUI*ULJ                                       ! BEGIN INITIALIZE Y: Y(1) = GAMMA, Y(2...3N+1) = COORDS, 
  J=2 

  DO I=1,N3_ATOM
    Y(J)= QN(I)
    J=J+1
  ENDDO

  DO I=1,N3_ATOM
    Y(J)=TAUI*MASSARRY(I)
    J=J+1       
    DO K=I+1,N3_ATOM
      Y(J)=0.0D0
      J=J+1
    ENDDO
  ENDDO

  TOUT=0.5D0*TAUMAX 

  CALL DLSODE(RHS, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
 
  CALL LNP(NEQ, Y, LOGZ)                                        ! CALCULATE LN OF PARTITION FUNCTION
  
  LNZ(1,1)=2*TOUT                                               ! STORE LAST TWO POINTS OF PARITION FUNCTION
  LNZ(2,1)=LOGZ+1.5*LOG(TOUT)    
  TOUT=0.5D0*(TAUMAX+0.1)

  CALL DLSODE(RHS, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
  CALL LNP(NEQ, Y, LOGZ)

  LNZ(1,2)=2*TOUT
  LNZ(2,2)=LOGZ+1.5*LOG(TOUT)         
      
  ENRG=(E_ZERO-(LNZ(2,2)-LNZ(2,1))/(LNZ(1,2)-LNZ(1,1)))         ! E = -d/dB*lnZ  (CALCULATE ENERGY)

  DO I=1,N3_ATOM
    QCNFG(I) = Y(I+1)/SIGMA
  ENDDO
 
  DEALLOCATE(RWORK, IWORK, Y)

END SUBROUTINE VGWQUENCH

SUBROUTINE VGWQUENCHSP(QCNFG, ENRG, CFLAG)
  USE VGWSP
  IMPLICIT NONE
  DOUBLE PRECISION :: QCNFG(N3_ATOM), ENRG, RAD, QTEMP(N3_ATOM), Q(N3_ATOM)
  DOUBLE PRECISION :: LNZ(2,2), LOGZ, T, GAMMA, QC(3), TOUT, ULJ
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: Y, RWORK
  DOUBLE PRECISION, PARAMETER :: TAUI=0.000001D0, CTOL=0.065D0,RTOL=0.0D0
  INTEGER :: NEQ, LRW, N_STEP, ITASK, ISTATE, IOPT, ITOL, I, J, K, CNT
  INTEGER, ALLOCATABLE, DIMENSION (:) :: IWORK
  INTEGER, PARAMETER :: LIW=20, MF=10
  LOGICAL CFLAG
  EXTERNAL RHSS
  EXTERNAL JAC

  CFLAG=.TRUE.
    
  NEQ=1+N_ATOM*(N_DIM+N_DIM*(N_DIM+1)/2)

  ITOL=1
  LRW=20+16*NEQ

  ALLOCATE(Y(NEQ), RWORK(LRW), IWORK(LIW))
  
  ITASK=1
  IOPT=1
  ISTATE=1
  IWORK(6)=100000
  IWORK(5)=4
  IWORK(7)=0
  IWORK(8)=0
  IWORK(9)=0
  IWORK(10)=0
  RWORK(5)=0.0D0
  RWORK(6)=0.0D0
  RWORK(7)=0.0D0
  RWORK(8)=0.0D0
  RWORK(9)=0.0D0
  RWORK(10)=0.0D0

  LNZ=0

  Y(1)=0

  DO I=1,N3_ATOM
    Y(I+1)=QCNFG(I)*SIGMA                       ! COPY INITITAL COORDS TO Y VECTOR
    Q(I)=Y(I+1)
  ENDDO
 
  CALL GAUSSENERGYPB(Q, N_ATOM, SIGMA, EPSILON, BL, ULJ)

  QC=0
  Y(1)=-TAUI*ULJ
 
  DO I=0,N_ATOM-1                             ! DETERMINE CENTER OF MASS VECTOR
    DO J=1,N_DIM
      QC(J)=QC(J)+Y(1+3*I+J)                    
    ENDDO
  ENDDO

  DO I=1,N_DIM
    QC(I)=QC(I)/N_ATOM
  ENDDO

  DO I=0,N_ATOM-1
    DO J=1,N_DIM
      Y(1+3*I+J)=Y(1+3*I+J)-QC(J)                 ! ZERO CENTER OF MASS
    ENDDO
  ENDDO

  CNT=2+N3_ATOM

  DO I=0,N_ATOM-1                             ! INITIALIZE Y VECTOR WITH INITITIAL CONDITIONS
    DO J=1,N_DIM
      Y(CNT)=TAUI*MASSARRY(3*I+J)
      CNT=CNT+1
      DO K=J+1,N_DIM
        Y(CNT)=0.0D0
        CNT=CNT+1
      ENDDO
    ENDDO
  ENDDO

  E_ZERO = ENRG

  T=TAUI

  DO I=1,NSTEP
    TOUT = TAU(I)
    CALL DLSODE(RHSS, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
    CALL LNPS(NEQ, Y, LOGZ)
    LNZ(1,1)=LNZ(1,2)
    LNZ(2,1)=LNZ(2,2)
    LNZ(1,2)=2*TAU(I)
    LNZ(2,2) = LOGZ+1.5D0*LOG(TAU(I))
  ENDDO
  
  ENRG =(E_ZERO - (LNZ(2,2)-LNZ(2,1))/(LNZ(1,2)-LNZ(1,1)))

  DO I=1,N3_ATOM
    QCNFG(I)=Y(I+1)/SIGMA
  ENDDO
 
  DEALLOCATE(Y, RWORK, IWORK)

END SUBROUTINE

SUBROUTINE LNP(NEQ, Y, LOGZ)                                 ! LOG OF DENSITY MATRIX
  USE VGW
  INTEGER :: I, J, K, NEQ
  DOUBLE PRECISION LDET, GAMMA, LOGZ, Q(N3_ATOM), C(N3_ATOM, N3_ATOM), Y(NEQ)

  GAMMA=Y(1)
  K=2

  DO I=1,N3_ATOM
    Q(I)=Y(K)                                                ! COPY COORDINATES FROM ARRAY Y()
    K=K+1
  ENDDO

  DO I=1,N3_ATOM
    DO J=I,N3_ATOM
      C(I,J) = Y(K)                                          ! GET GAUSSIAN WIDTH MATRIX FROM Y()
      C(J,I) = Y(K)
      K=K+1
    ENDDO
  ENDDO

  CALL LNDET(LDET, C)
  LOGZ=2.0D0*GAMMA-0.5D0*LDET     

END SUBROUTINE LNP

SUBROUTINE LNPS(NEQ, Y, LOGZ)
  USE VGWSP
  IMPLICIT NONE
  INTEGER :: I, J, K, NEQ, CNT
  DOUBLE PRECISION :: Y(NEQ), C(N_DIM,N_DIM), M(N_DIM,N_DIM), LOGZ, GAMMA, DETI, DET

  GAMMA=Y(1)
  CNT=2+N3_ATOM

  DET=0.0D0

  DO I=1,N_ATOM
    DO J=1,N_DIM
      C(J,J) = Y(CNT)
      CNT=CNT+1
      DO K=J+1,N_DIM
        C(J,K)= Y(CNT)
        C(K,J)= C(J,K)
        CNT=CNT+1
      ENDDO
    ENDDO
    CALL INVDET(C,M,DETI)
    DET = DET + LOG(DETI)
  ENDDO

  LOGZ = 2.0D0*GAMMA - 0.5D0*DET

END SUBROUTINE 

SUBROUTINE LNDET(LDET, C)
  USE VGW
  INTEGER :: I, J, K 
  DOUBLE PRECISION :: LDET, SUML, C(N3_ATOM, N3_ATOM), P(N3_ATOM)

  LDET=0.0D0 
     
  DO I=1,N3_ATOM
    DO J=I, N3_ATOM
      SUML=C(I,J)
      K=I-1
      DO WHILE(K.GT.0)
        SUML=SUML-(C(I,K)*C(J,K))
        K=K-1
      ENDDO
      IF(I.EQ.J) THEN
        P(I)=SQRT(SUML)
        LDET=LDET+LOG(SUML)
      ELSE 
        C(J,I)=SUML/P(I)
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE LNDET

SUBROUTINE INVDET(A, M, DET)
  IMPLICIT NONE
  DOUBLE PRECISION :: DET, A(3,3), M(3,3)

  M(1,1) = A(2,2)*A(3,3)-A(2,3)*A(2,3)
  M(2,2) = A(1,1)*A(3,3)-A(1,3)*A(1,3)
  M(3,3) = A(1,1)*A(2,2)-A(1,2)*A(1,2)
  M(1,2) = -A(1,2)*A(3,3)+A(1,3)*A(2,3)
  M(1,3) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
  M(2,3) = -A(1,1)*A(2,3)+A(1,3)*A(1,2)
  DET = M(1,1)*A(1,1)+M(1,2)*A(1,2)+M(1,3)*A(1,3)

END SUBROUTINE INVDET  

SUBROUTINE RHSS(NEQ, T, Y, YPRIME)
  USE VGWSP
  IMPLICIT NONE
  INTEGER :: I, J, K, I1, I2, IG, NEQ, CNT
  DOUBLE PRECISION :: COE, FACTOR, U, UX, UXX, QZQ, EXPAV, GG, GG2, TRMG, DETS, DETI
  DOUBLE PRECISION :: M(N_DIM,N_DIM), A(N_DIM,N_DIM), R(N_DIM), Z(N_DIM,N_DIM)
  DOUBLE PRECISION :: DETA, DETAG, GUG, QP, Q(3*N_ATOM), BLKC(N3_ATOM,N_DIM)      ! ARRAY of N 3x3 BLOCKS (BLOCK DIAGONAL OF C)
  DOUBLE PRECISION :: UPV(N3_ATOM), UPM(N3_ATOM, N_DIM), Q2, Q12(N_DIM), TRUXXGI
  DOUBLE PRECISION :: T, AG(N_DIM,N_DIM), GU(N_DIM,N_DIM), Y(NEQ), YPRIME(NEQ)

  CNT=3*N_ATOM+2

  DO I=1,3*N_ATOM
    Q(I)=Y(I+1)
  ENDDO

  TRMG=0.0D0
  
  DO I=0,N_ATOM-1  
    DO J=1,N_DIM
      DO K=J,N_DIM
        BLKC(3*I+J,k)= Y(CNT)
        BLKC(3*I+K,J)= Y(CNT) 
        CNT=CNT+1
      ENDDO
    ENDDO
  ENDDO

  U=0.0D0
  UPV=0.0D0
  UPM=0.0D0
   
  DO I1=0,N_ATOM-1
    DO I2=I1+1,N_ATOM-1
      Q2=0.0D0          
      DO I=1,N_DIM
        Q12(I)=Q(3*I1+I)-Q(3*I2+I)
        Q2 = Q2 + Q12(I)*Q12(I)
      ENDDO
      
      DO J=1,N_DIM
        DO K=J,N_DIM
          A(J,K)=BLKC(3*I1+J,K)+BLKC(3*I2+J,K)
        ENDDO
      ENDDO

      CALL INVDET(A,M,DETI)       
      DETA=1.0D0/DETI
           
      DO J=1,N_DIM
        DO K=J,N_DIM
          A(J,K)=M(J,K)*DETA     
        ENDDO
      ENDDO
                       
      DO IG=1,NGAUSS                                  ! BEGIN SUMMATION OVER GAUSSIANS
        COE=LJK(IG)
        GG=LJG(IG)
        GG2=LJG2(IG)
        DO J=1,N_DIM
          AG(J,J)=GG+A(J,J)
          DO K=J+1,N_DIM
            AG(J,K)=A(J,K)
          ENDDO
        ENDDO
     
        CALL INVDET(AG,M,DETAG)
        FACTOR=GG2/DETAG  
 
        DO I=1,N_DIM
          Z(I,I)=GG-FACTOR*M(I,I)
          DO J=I+1,N_DIM
            Z(I,J)=-FACTOR*M(I,J)
            Z(J,I)=Z(I,J)
          ENDDO
        ENDDO          
     
        QZQ=0.0D0

        DO I=1,N_DIM
          R(I)=0.0D0
          DO J=1,N_DIM
            R(I)=R(I)+Z(I,J)*Q12(J)
          ENDDO
          QZQ=QZQ+R(I)*Q12(I)
          R(I)=-2.0D0*R(I)
        ENDDO
     
        EXPAV=SQRT(DETA/DETAG)*EXP(-QZQ)
        U=U+EXPAV*COE

        DO I=1,N_DIM
          UX=EXPAV*R(I)*COE
          UPV(3*I1+I)=UPV(3*I1+I)+UX                          ! COMPUTE GRADIENT (UPV = GRADIENT VECTOR)
          UPV(3*I2+I)=UPV(3*I2+I)-UX
          DO J=I,N_DIM
            UXX=(R(I)*R(J)-2*Z(I,J))*EXPAV*COE                ! HESSIAN MATRIX (BLOCK OF I,JTH PARTICLE)
            UPM(3*I1+I,J)=UPM(3*I1+I,J)+UXX                   ! UPM = BLOCK DIAGONAL OF 3Nx3N HESSIAN
            UPM(3*I2+I,J)=UPM(3*I2+I,J)+UXX
            IF(I.NE.J) THEN                                   ! FILL LOWER HALF OF MATRIX
              UPM(3*I1+J,I)=UPM(3*I1+J,I)+UXX
              UPM(3*I2+J,I)=UPM(3*I2+J,I)+UXX
            ENDIF
          ENDDO
        ENDDO
      ENDDO                                                   ! END SUMMATION OVER GUASSIANS
    ENDDO                                                     ! END OF I1 LOOP (ITH PARTICLE)     
  ENDDO                                                       ! END OF I2 LOOP (JTH PARTICLE)
  
  IF(CONPOT.EQ.1) THEN                                        ! CONSTRAINING POTENTIAL
    DO I1=0,N_ATOM-1
      Q2=0.0D0          
      DO I=1,N_DIM
        Q12(I)=Q(3*I1+I)
        Q2 = Q2 + Q12(I)*Q12(I)
      ENDDO
      
      DO J=1,N_DIM
        DO K=J,N_DIM
          A(J,K)=BLKC(3*I1+J,K)
        ENDDO
      ENDDO

      CALL INVDET(A,M,DETI)       
      DETA=1.0D0/DETI
           
      DO J=1,N_DIM
        DO K=J,N_DIM
          A(J,K)=M(J,K)*DETA     
        ENDDO
      ENDDO
                       
      DO IG=1,3                                  ! BEGIN SUMMATION OVER CONSTRAINING GAUSSIANS
        COE=CGK(IG)
        GG=CGG(IG)
        GG2=CGG2(IG)
        DO J=1,N_DIM
          AG(J,J)=GG+A(J,J)
          DO K=J+1,N_DIM
            AG(J,K)=A(J,K)
          ENDDO
        ENDDO
     
        CALL INVDET(AG,M,DETAG)
        FACTOR=GG2/DETAG  
 
        DO I=1,N_DIM
          Z(I,I)=GG-FACTOR*M(I,I)
          DO J=I+1,N_DIM
            Z(I,J)=-FACTOR*M(I,J)
            Z(J,I)=Z(I,J)
          ENDDO
        ENDDO          
     
        QZQ=0.0D0

        DO I=1,3
          R(I)=0.0D0
          DO J=1,3
            R(I)=R(I)+Z(I,J)*Q12(J)
          ENDDO
          QZQ=QZQ+R(I)*Q12(I)
          R(I)=-2.0D0*R(I)
        ENDDO
     
        EXPAV=SQRT(DETA/DETAG)*EXP(-QZQ)
        U=U+EXPAV*COE

        DO I=1,3
          UX=EXPAV*R(I)*COE
          UPV(3*I1+I)=UPV(3*I1+I)+UX                          ! COMPUTE GRADIENT (UPV = GRADIENT VECTOR)
          DO J=I,3
            UXX=(R(I)*R(J)-2*Z(I,J))*EXPAV*COE                ! HESSIAN MATRIX (BLOCK OF I,JTH PARTICLE)
            UPM(3*I1+I,J)=UPM(3*I1+I,J)+UXX                   ! UPM = BLOCK DIAGONAL OF 3Nx3N HESSIAN
            IF(I.NE.J) THEN                                   ! FILL LOWER HALF OF MATRIX
              UPM(3*I1+J,I)=UPM(3*I1+J,I)+UXX
            ENDIF
          ENDDO
        ENDDO
      ENDDO                                                   ! END SUMMATION OVER GAUSSIANS  
      U=U+CGCONST                                           
    ENDDO                                                     ! END OF I1 LOOP (ITH PARTICLE)     
  ENDIF
  TRUXXGI=0.0D0

  DO I1=0,N_ATOM-1
    DO I=1,N_DIM
      TRUXXGI=TRUXXGI+UPM(3*I1+I,I)*BLKC(3*I1+I,I)
      DO J=I+1,N_DIM
        TRUXXGI=TRUXXGI+2.0D0*UPM(3*I1+I,J)*BLKC(3*I1+I,J)
      ENDDO
    ENDDO
  ENDDO

  CNT=1

  YPRIME(CNT)=-0.25D0*TRUXXGI-U+E_ZERO                 ! *(pnt++)=0.25*TrUxxGI-0.5*TrMG-U+E_zero+1.5*N_atom/(*T);
  CNT=CNT+1  

  DO I1=0,N_ATOM-1
    DO I=1,N_DIM
      QP=0.0D0
      DO J=1,N_DIM
        QP = QP-BLKC(3*I1+I,J)*UPV(3*I1+J)
      ENDDO
      YPRIME(CNT)=QP
      CNT=CNT+1
    ENDDO
  ENDDO

  DO I1=0,N_ATOM-1
    DO I=1,N_DIM
      DO J=1,N_DIM
        GU(I,J)=0.0D0
        DO K=1,N_DIM
          GU(I,J)=GU(I,J)+BLKC(3*I1+I,K)*UPM(3*I1+K,J)
        ENDDO
      ENDDO
    ENDDO
          
    DO I=1,N_DIM
      DO J=I, N_DIM
        GUG=0.0D0
        DO K=1,N_DIM
          GUG=GUG-GU(I,K)*BLKC(3*I1+K,J)
        ENDDO
        IF(I.EQ.J) THEN
          GUG=GUG+MASSARRY(3*I1+J)
        ENDIF
        YPRIME(CNT) = GUG
        CNT=CNT+1
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE RHS(NEQ, T, Y, YPRIME)
  USE VGW
  IMPLICIT NONE
  INTEGER :: NEQ, I, J, K, IG, X1, X2, I1, I2, CNT
  DOUBLE PRECISION :: T, Y(NEQ), YPRIME(NEQ), UXP, UXXP, Q(N3_ATOM), C(N3_ATOM, N3_ATOM)
  DOUBLE PRECISION :: A(3,3), AG(3,3), M(3,3), Z(3,3), Q12(3), R(3), CU(N3_ATOM), UX(N3_ATOM)
  DOUBLE PRECISION :: UXX(N3_ATOM, N3_ATOM), DETA, DETAG
  DOUBLE PRECISION :: COE, TRCUXX, CUC, EXPAV, QZQ, GG, GG2, FACTOR, QP, GP

  CNT=N3_ATOM+2

  DO I=1,N3_ATOM
    Q(I)=Y(I+1)
  ENDDO
 
  DO I=1,N3_ATOM
    DO J=I,N3_ATOM
      C(I,J)=Y(CNT)
      C(J,I)=Y(CNT)
      CNT=CNT+1
    ENDDO
  ENDDO

  U=0.0D0                                           ! INITIALIZE POTENTIAL ENERGY
  UX=0.0D0                                          ! INITITALIZE GRADIENT VECTOR
  UXX=0.0D0                                         ! INITIALIZE HESSIAN MATRIX

  DO I1=0,N_ATOM-1                                  ! CALCULATE ALL Rij
    X1=3*I1
    DO I2=I1+1,N_ATOM-1
      X2=3*I2
      DO I=1,3
        Q12(I)=Q(X1+I)-Q(X2+I)                      ! ALL Rij VECTOR
      ENDDO

      DO I=1,3
        DO J=I,3
          A(I,J)=C(X1+I,X1+J)-C(X1+I,X2+J)-C(X1+J,X2+I)+C(X2+I,X2+J)                ! CONSTRUCT MAXTRIX "A"
        ENDDO
      ENDDO

      CALL INVDET(A,M,DETA)                          ! COMPUTE INVERSE DETERMINANT OF "A"
      DETA=1.0D0/DETA                                ! TAKE INVERSE OF INVERSE DET. OF "A" TO GET DET. OF "A"

      DO I=1,N_DIM
        DO J=I,N_DIM
          A(I,J)=M(I,J)*DETA
        ENDDO
      ENDDO
                                                    ! ******* BEGIN SUMMATION OVER GUASSIANS *******
      DO IG=1,NUM_G
        COE=LJK(IG)                                 ! GET COEFFICIENT Ci OF iTH GAUSSIAN TERM IN POTENTIAL
        GG=LJG(IG)                                  ! GET EXPONENT COEFFICIENT (WIDTH) ALPHA OF iTH GAUSSIAN
        GG2=LJG2(IG)                                ! COMPUTE ALPHA SQUARED
        DO I=1,3
          AG(I,I)=GG+A(I,I)                         ! MATRIX "A" PLUS GAUSSIAN COEFFICEINT ALPHA
          DO J=I+1,N_DIM
            AG(I,J)=A(I,J)
          ENDDO
        ENDDO 

        CALL INVDET(AG,M,DETAG)                             ! INVERSE DETERMIANT OF MATRIX "A" + ALPHA
        FACTOR=GG2/DETAG

        DO I=1,3
          Z(I,I)=GG-FACTOR*M(I,I)                   ! CONSTRUCT MATRIX Z (Zij = ALPHA - ALPHA^2(ALPHA + Aij)^-1)
          DO J=I+1,3
            Z(I,J)=-(FACTOR*M(I,J))
            Z(J,I)=Z(I,J)
          ENDDO
        ENDDO

        QZQ=0

        DO I=1,3
          R(I)=0.0D0
          DO J=1,3
            R(I)=R(I)+Z(I,J)*Q12(J)
          ENDDO
          QZQ=QZQ+R(I)*Q12(I)                       ! CALCULATE ENTIRE GAUSSIAN ARGUMENT (-Qij(T)*Zij*Qij)
          R(I)=-(2*R(I))
        ENDDO

        EXPAV=SQRT(DETA/DETAG)*EXP(-QZQ)
        U=U+(EXPAV*COE);
      
        DO I=1,3
          UXP=EXPAV*R(I)*COE;                              ! COMPUTE GRADIENT (FORCES) FOR ijTH PARTICLE        
          UX(X1+I)=UX(X1+I)+UXP
          UX(X2+I)=UX(X2+I)-UXP
          DO J=I,3
            UXXP=(R(I)*R(J)-2*Z(I,J))*EXPAV*COE;
            UXX(X1+I,X1+J)=UXX(X1+I,X1+J)+UXXP             ! COMPUTE HESSIAN MATRIX
            UXX(X2+I,X2+J)=UXX(X2+I,X2+J)+UXXP       
            UXX(X1+J,X2+I)=UXX(X1+J,X2+I)-UXXP     
            IF(I.NE.J) THEN
              UXX(X1+I,X2+J)=UXX(X1+I,X2+J)-UXXP
            ENDIF
          ENDDO
        ENDDO
      ENDDO                                                ! END OF SUMMATION OVER GAUSSIANS
    ENDDO                                                  ! END OF I2 LOOP (jTH PARTICLE)
  ENDDO  
                                                  ! END OF I1 LOOP (iTH PARTICLE)
  IF(CONPOT.EQ.1) THEN
  DO I1=0,N_ATOM-1                                  ! CONSTRAINING POTENTIAL
    X1=3*I1
    DO I=1,3
      Q12(I)=Q(X1+I)
    ENDDO

    DO I=1,3
      DO J=I,3
        A(I,J)=C(X1+I,X1+J)
      ENDDO
    ENDDO

    CALL INVDET(A,M,DETA)                          ! COMPUTE INVERSE DETERMINANT OF "A"
    DETA=1.0D0/DETA                                ! TAKE INVERSE OF INVERSE DET. OF "A" TO GET DET. OF "A"

      DO I=1,N_DIM
        DO J=I,N_DIM
          A(I,J)=M(I,J)*DETA
        ENDDO
      ENDDO
                                                    ! ******* BEGIN SUMMATION OVER CP GAUSSIANS *******
      DO IG=1,3
        COE=CGK(IG)                                 ! GET COEFFICIENT Ci OF iTH GAUSSIAN TERM IN POTENTIAL
        GG=CGG(IG)                                  ! GET EXPONENT COEFFICIENT (WIDTH) ALPHA OF iTH GAUSSIAN
        GG2=CGG2(IG)                                ! COMPUTE ALPHA SQUARED
        DO I=1,3
          AG(I,I)=GG+A(I,I)                         ! MATRIX "A" PLUS GAUSSIAN COEFFICEINT ALPHA
          DO J=I+1,N_DIM
            AG(I,J)=A(I,J)
          ENDDO
        ENDDO 

        CALL INVDET(A,M,DETAG)                             ! INVERSE DETERMIANT OF MATRIX "A" + ALPHA
        FACTOR=GG2/DETAG

        DO I=1,3
          Z(I,I)=GG-FACTOR*M(I,I)                   ! CONSTRUCT MATRIX Z (Zij = ALPHA - ALPHA^2(ALPHA + Aij)^-1)
          DO J=I+1,3
            Z(I,J)=-(FACTOR*M(I,J))
            Z(J,I)=Z(I,J)
          ENDDO
        ENDDO

        QZQ=0

        DO I=1,3
          R(I)=0.0D0
          DO J=1,3
            R(I)=R(I)+Z(I,J)*Q12(J)
          ENDDO
          QZQ=QZQ+R(I)*Q12(I)                       ! CALCULATE ENTIRE GAUSSIAN ARGUMENT (-Qij(T)*Zij*Qij)
          R(I)=-(2*R(I))
        ENDDO

        EXPAV=SQRT(DETA/DETAG)*EXP(-QZQ)
        U=U+(EXPAV*COE);
      
        DO I=1,3
          UXP=EXPAV*R(I)*COE;                              ! COMPUTE GRADIENT (FORCES) FOR ijTH PARTICLE        
          UX(X1+I)=UX(X1+I)+UXP
          DO J=I,3
            UXXP=(R(I)*R(J)-2*Z(I,J))*EXPAV*COE;
            UXX(X1+I,X1+J)=UXX(X1+I,X1+J)+UXXP             ! COMPUTE HESSIAN MATRIX
          ENDDO
        ENDDO
      ENDDO                                                ! END OF SUMMATION OVER GAUSSIANS
      U=U+CGCONST
    ENDDO  
    ENDIF
  
  TRCUXX=0.0D0                                             ! ******* BEGIN CALCULATION OF RIGHT HAND SIDE OF DIFF. EQ. *******
  CNT=2                                                    ! C'(TAU) = -C<Del*Del U>C + HBAR^2*M^-1
  
  DO I=1,N3_ATOM
    QP=0.0D0
    DO J=1,N3_ATOM
      QP=QP-(C(I,J)*UX(J))
    ENDDO
    YPRIME(CNT)=QP                                         ! COMPUTE NEW VELOCTIES
    CNT=CNT+1
  ENDDO
    
  DO I=1,N3_ATOM                                           ! RHS OF GAUSSIAN WIDTH MATRIX
    DO J=1, N3_ATOM
      CU(J)=0.0D0
      DO K=1,J-1
        CU(J)=CU(J)+(C(I,K)*UXX(K,J))
      ENDDO
      DO K=J,N3_ATOM
        CU(J)=CU(J)+C(I,K)*UXX(J,K)
      ENDDO
    ENDDO

    TRCUXX=TRCUXX+CU(I) 
   
    DO J=I,N3_ATOM
      CUC=0.0D0
      DO K=1,N3_ATOM
        CUC=CUC+(CU(K)*C(K,J))
      ENDDO
      GP=-CUC
      IF(I.EQ.J) THEN
        GP=GP+MASSARRY(I)
      ENDIF
      YPRIME(CNT)=GP
      CNT=CNT+1
    ENDDO
  ENDDO
  
  YPRIME(1)=-0.25D0*TRCUXX-U+E_ZERO                        ! CALCULATE NEW GAMMA PRIME

END SUBROUTINE RHS

SUBROUTINE GAUSSENERGYPB(Q, N, LJS, EPS, BL, U)
  IMPLICIT NONE
  INTEGER :: I, J, K, N
  DOUBLE PRECISION :: Q(3,N),U,LJS,EPS,RSQ,SIG2,BL,BL2,QIJ
  DOUBLE PRECISION :: LJA(3), LJC(3)

  BL2=BL/2
  SIG2=LJS**2
  U=0

  LJC(1)=1840D0
  LJC(2)=-1.48D0
  LJC(3)=-23.2D0
  LJA(1)=-6.65D0/SIG2
  LJA(2)=-0.79D0/SIG2
  LJA(3)=-2.60D0/SIG2

  DO I=1,N-1
    DO J=I+1,N
      RSQ=0
      DO K=1,3
        QIJ=Q(K,I)-Q(K,J)
        IF((ABS(QIJ) > BL2).AND.(BL.GT.-1)) THEN             ! PERIODIC BOUNDARY CONDITIONS
          QIJ=QIJ-BL*(QIJ/ABS(QIJ))
        ENDIF
        RSQ=RSQ+QIJ**2
      ENDDO
      DO K=1,3
        U=U+LJC(K)*EXP(LJA(K)*RSQ)
      ENDDO
    ENDDO
  ENDDO
  U=EPS*U

ENDSUBROUTINE


SUBROUTINE CLRADIUS(Q, N, RAD)
  IMPLICIT NONE
  INTEGER :: N, I, J
  DOUBLE PRECISION :: Q(3*N), RCM(3), CV(3), CR, RAD

  RAD=0
  RCM=0
  
  DO I=0,N-1
    DO J=1,3
      RCM(J)= RCM(J)+Q(3*I+J)
    ENDDO
  ENDDO

  DO I=1,3
    RCM(I)=RCM(I)/N
  ENDDO

  DO I=0,N-1
    CR=0
    DO J=1,3
      CV(J)=Q(3*I+J)-RCM(J)
      CR=CR+CV(J)**2
    ENDDO
    CR = SQRT(CR)
    IF(CR.GT.RAD) THEN
      RAD=CR
    ENDIF
  ENDDO

END SUBROUTINE
