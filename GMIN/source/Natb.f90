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
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!       TB potential for Na, Ag, Li clusters
!
! khs26> Some comments were originally in French. I have translated them
! as best I can.

SUBROUTINE NATB(NSIZE, P, GRAD, E0, GTEST, GUIDET)
   IMPLICIT NONE
   INTEGER ION, NSIZE
   DOUBLE PRECISION, PARAMETER :: EPSILON=1.0D-10
   LOGICAL GTEST, STEST, GUIDET
   DOUBLE PRECISION P(3*NSIZE), GRAD(3*NSIZE), E0
   CHARACTER(LEN=2) METAL
   COMMON /ION/ ION
   COMMON /METAL/ METAL

   METAL='NA'
   ION=0

!       IF (metal == 'NA') then
!          WRITE(*, *) 'Parameters optimized for SODIUM'
!       elseIF (metal == 'AG') then
!          WRITE(*, *) 'Parameters optimized for SILVER'
!       elseIF (metal == 'LI') then
!          WRITE(*, *) 'Parameters optimized for LITHIUM'
!       END IF

!       IF (ion == 1) then
!          WRITE(*, *) 'CATION'
!       elseIF (ion == 0) then
!          WRITE(*, *) 'NEUTRAL'
!       else
!          WRITE(*, *) 'ANION'
!       END IF

!       WRITE(*, *) 'GUIDET=', GUIDET
   IF (GUIDET) P(1:3*NSIZE)=P(1:3*NSIZE)*6.02D0
   CALL ENTOTS(P, GRAD, NSIZE, E0)
   IF (GUIDET) THEN
      P(1:3*NSIZE)=P(1:3*NSIZE)/6.02D0
      GRAD(1:3*NSIZE)=GRAD(1:3*NSIZE)*6.02D0
   END IF

END SUBROUTINE NATB

SUBROUTINE ENTOTS(P, DP, NBAT, E0)
   IMPLICIT NONE
   DOUBLE PRECISION ANM, RKB, EREPMAX
   PARAMETER(ANM=41901.052, RKB=3.165D-6, EREPMAX=100.D0)

   INTEGER NBAT
   DOUBLE PRECISION P(3*NBAT), E0
   DOUBLE PRECISION GRAD(3*NBAT), DP(3*NBAT)

   DOUBLE PRECISION HDER(3*NBAT, NBAT*(NBAT+1)/2), R(NBAT*(NBAT+1)/2)
   DOUBLE PRECISION HA(NBAT*(NBAT+1)/2)
   DOUBLE PRECISION VECT(NBAT*(NBAT+1)/2), RINV(NBAT*(NBAT+1)/2)
   DOUBLE PRECISION H(NBAT*(NBAT+1)/2), DIF(3*NBAT*(NBAT+1)/2)
   DOUBLE PRECISION DBSZ(3*NBAT*(NBAT+1)/2), UNIT(3*NBAT*(NBAT+1)/2)
   DOUBLE PRECISION VALP(NBAT), VT(NBAT, NBAT), NUM(NBAT)
   DOUBLE PRECISION HMAT(NBAT, NBAT)
   DOUBLE PRECISION AIJK, APOL, HIJ, HII, HH, ASZ, BIJK, DENOM, DER
   DOUBLE PRECISION DERQ, DI, DII, DISTOR, DJ, DJJ, DK, EN, EREP, ESP
   DOUBLE PRECISION FTOL, GRD, HIJK, HSS
   DOUBLE PRECISION RIJ, RIJ14, RIJ2, RIJK, RIK, RMAX, RMIN, SCAL
   DOUBLE PRECISION UIK, UJK, VAL, VGR, XIJ, YIJ, ZIJ
   INTEGER II, III, IDIAG, I, IDER, IJX, IJY, IJZ, IJQL, IJ, IJQ, IK
   INTEGER IKQ, IKQL, IND, IPARAM, IPRTG, IPRTH, ISI, ISJ, ITMAX
   INTEGER J, JJ, JJJ, JK, JKQ, JKQL, K, KKK, L, M, ML, NALK, NBLEC, NDIMH
   INTEGER NDOC, NGR, NOC
   INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NBAT)
   INTEGER IWORK(33*3*NBAT), INFO, ISTAT
   DOUBLE PRECISION WORK(33*3*NBAT), ABSTOL, DLAMCH
   INTEGER ION
   COMMON /ION/ ION
   CHARACTER(LEN=2) METAL
   COMMON /METAL/ METAL

   DIMENSION DI(3), DJ(3)

   LOGICAL MONOC

!        PARAMETERS

   LWORK=33*3*NBAT
   ILWORK=33*3*NBAT

   IF (METAL == 'NA') THEN
      ESP=0.077307
      RMIN=3.0D0
      RMAX=100.0D0
   ELSE IF (METAL == 'AG') THEN
      ESP=0.227857
      RMIN=4.0D0
      RMAX=20.0D0
   ELSE IF (METAL == 'LI') THEN
      ESP=0.067908D0
      RMIN=0.0D0
      RMAX=1000.0D0
   END IF

   NALK=NBAT
   IDIAG=1
   DISTOR=0.0D0
   HH=1.0D-6
   FTOL=1.0D-4
   ITMAX=100
   IPARAM=1
   ASZ=1.0D0
   NGR=0
   IPRTH=0
   IPRTG=0
   NBLEC=NBAT-ION
   APOL=0.0D0

   NDIMH=NBAT ! CHANGED BY DJW

   H(1:NBAT)=0.0D0
   HDER(1:3*NBAT, 1:NBAT)=0.0D0
   GRAD(1:3*NBAT)=0.0D0
   DENOM=1.0D0/ESP 

   IJ=0
   IJQ=0
   DO I=1, NBAT
      DO J=1, I-1
         IJ=IJ+1
         IJX=IJQ+1
         IJY=IJQ+2
         IJZ=IJQ+3 
         DIF(IJX)=P(3*J-2)-P(3*I-2)
         DIF(IJY)=P(3*J-1)-P(3*I-1)
         DIF(IJZ)=P(3*J)-P(3*I)
         RIJ=DIF(IJX)*DIF(IJX)+DIF(IJY)*DIF(IJY)+DIF(IJZ)*DIF(IJZ)
         RIJ=DSQRT(RIJ)
         R(IJ)=RIJ
         RINV(IJ)=1.D0/RIJ
         UNIT(IJX)=DIF(IJX)*RINV(IJ)
         UNIT(IJY)=DIF(IJY)*RINV(IJ)
         UNIT(IJZ)=DIF(IJZ)*RINV(IJ)
         IJQ=IJQ+3
      END DO
      IJ=IJ+1
      IJQ=IJQ+3
   END DO

! Initialisation
   HDER(1:3*NBAT, 1:(NBAT*(NBAT+1)/2))=0.0D0

! Storage of perturbation terms
   IJQ=0
   IJ=0
   DO I=1, NALK
      DO J=1, I-1
         IJ=IJ+1
         RIJ=R(IJ)
         IF ((RIJ < RMAX) .AND. (RIJ >= RMIN)) THEN
            CALL FTSZ(RIJ, VAL, DER)
         ELSE IF (RIJ < RMIN) THEN
            IF (METAL == 'NA' .OR. METAL == 'AG' .OR. METAL == 'LI') THEN
               VAL=-3.724D-10*EXP(3.32*RIJ)
               DER=3.32*VAL
            END IF
         ELSE
            VAL=0.D0
            DER=0.D0
         END IF
         VECT(IJ)=VAL
         RINV(IJ)=1.D0/RIJ 
         DO L=1, 3
            IJQL=IJQ+L
            DBSZ(IJQL)=DER*UNIT(IJQL)
         END DO
         IJQ=IJQ+3
      END DO
      IJQ=IJQ+3
      IJ=IJ+1
   END DO

! Hamiltonian and gradient
   IJ=0
   III=0
   IJQ=0
   DO I=1, NALK
      JJJ=0
! Non-diagonal elements (beta effective)
      DO J=1, I-1
         IJ=IJ+1
         RIJ=R(IJ)
         IF (RIJ >= RMAX) THEN
            HSS=0.D0
            DER=0.D0
         ELSE
            IF (RIJ >= RMIN) THEN
               CALL FTSS(RIJ, VAL, DER)
            ELSE
               IF (METAL == 'NA') THEN
                  VAL=-3.743D-8*DEXP(2.544*RIJ)
                  DER=2.544*VAL
               ELSEIF (METAL == 'AG') THEN
                  VAL=-3.743D-8*DEXP(2.544*RIJ)
                  DER=2.544*VAL
               ELSEIF (METAL == 'LI') THEN
                  VAL=-3.743D-8*DEXP(2.544*RIJ)
                  DER=2.544*VAL
               END IF
            END IF
            HSS=VAL
         END IF
         IDER=0
         DO L=1, 3 
            IJQL=IJQ+L
            DERQ=DER*UNIT(IJQL)
            HDER(IDER+I, IJ)=-DERQ
            HDER(IDER+J, IJ)=DERQ
            IDER=IDER+NBAT
         END DO
         DI(1:3)=0.0D0
         DJ(1:3)=0.0D0
         KKK=0  
! Perturbed Hamiltonian (3-body terms)
         HIJ=0.D0
         DO K=1, NALK
            IF ((K == J).OR.(K == I)) THEN
               KKK=KKK + K
               CYCLE
            END IF
            IF (I > K) THEN
               ISI=1
               IK=III+K
            ELSE
               IK=KKK+I
               ISI=-1
            END IF
            IF (J > K) THEN
               ISJ=1
               JK=JJJ+K
            ELSE
               JK=KKK+J
               ISJ=-1
            END IF
            IKQ=IK+IK+IK-3
            JKQ=JK+JK+JK-3
            SCAL=UNIT(IKQ+1)*UNIT(JKQ+1)+UNIT(IKQ+2)*UNIT(JKQ+2) &
                +UNIT(IKQ+3)*UNIT(JKQ+3)
            SCAL=SCAL*ISI*ISJ
            BIJK=VECT(IK)*VECT(JK)
            HIJK=BIJK*SCAL
            HIJ=HIJ-HIJK
            RIJK=RINV(IK)*RINV(JK)
! Gradient of SP terms
            IDER=0
            DO L=1, 3
               IKQL=IKQ+L
               JKQL=JKQ+L
               UIK=ISI*UNIT(IKQL)
               UJK=ISJ*UNIT(JKQL)
               AIJK=BIJK*RIJK
               DII=-DBSZ(IKQL)*VECT(JK)*ISI*SCAL &
                  -ISJ*DIF (JKQL)*AIJK &
                  +HIJK*UIK*RINV(IK)
               DJJ=-DBSZ(JKQL)*VECT(IK)*ISJ*SCAL &
                  -ISI*DIF (IKQL)*AIJK &
                  +HIJK*UJK*RINV(JK)
               DK=-DII-DJJ
               DI(L)=DI(L)+DII
               DJ(L)=DJ(L)+DJJ
               HDER(IDER+K, IJ)=-DK*DENOM
               IDER=IDER+NBAT
            END DO
         END DO

         H(IJ)=HSS+HIJ*DENOM
         IDER=0 
         DO L=1, 3
            HDER(IDER+I, IJ)=HDER(IDER+I, IJ)-DI(L)*DENOM
            HDER(IDER+J, IJ)=HDER(IDER+J, IJ)-DJ(L)*DENOM
            IDER=IDER+NBAT
         END DO
         JJJ=JJJ+J
         IJQ=IJQ+3
      END DO
! Diagonal elements (repulsion)
      IJ=IJ+1
      KKK=0
      HII=0.D0 
      DO K=1, NALK
         IF (K == I) THEN
            KKK=KKK + K
            CYCLE
         END IF
         IF (I > K) THEN
            IK=III+K
            ISI=1
         ELSE
            IK=KKK+I
            ISI=-1
         END IF
         IKQ=IK+IK+IK-3
! Hamiltonian
         VAL=0.D0 
         DER=0.D0
         RIK=R(IK)
         IF (RIK <= RMAX) THEN
            IF (RIK >= RMIN) THEN
               CALL FRSS(RIK, VAL, DER)
            ELSE
               VAL=1.D0
               DER=0.D0
            END IF
            HII=HII+VAL
         END IF
! Gradient
         IDER=0
         DO L=1, 3
            IKQL=IKQ+L
            DERQ=ISI*DER*UNIT(IKQL)
            HDER(K+IDER, IJ)=DERQ
            HDER(I+IDER, IJ)=HDER(I+IDER, IJ)-DERQ
            IDER=IDER+NBAT
         END DO
         KKK=KKK+K
      END DO
      H(IJ)=HII
      III=III+I
      IJQ=IJQ+3
   END DO

! Diagonalisation of H

! Number of occupied orbitals
   NDOC=NBLEC/2
   NOC=NDOC
   MONOC=.FALSE.
   IF ((2*NDOC).NE.NBLEC) THEN
      MONOC=.TRUE. 
      NOC=NDOC+1
   END IF
   DO I=1, NALK
      NUM(I)=I
   END DO
! Save the Hamiltonian in HA    
   DO IJ=1, NALK*(NALK+1)/2
      HA(IJ)=H(IJ)
   END DO
   IF (NALK == 1) THEN
      VALP(1)=-H(1)
      VECT(1)=1.D0
   END IF
   IF (NALK > 1) THEN
      IF (IDIAG == 1) THEN
! HA contains the lower-triangular matrix
! H contains the upper-triangular matrix
! H is destroyed after diagonalisation       
         II=0
         IJ=0
         DO I=1, NALK
            JJ=II
            DO J=I, NALK
               IJ=IJ+1
               H(IJ)=HA(JJ+I)
               JJ=JJ+J
            END DO
            II=II+I
         END DO

! start new DJW
         HMAT(1:NBAT, 1:NBAT)=0.0D0
         II=0
         DO I=1, NALK
            DO J=I, NALK
               II=II+1
               HMAT(I, J)=H(II)
            END DO
         END DO

         ABSTOL=DLAMCH('SAFE  MINIMUM')
         CALL DSYEVR('V','I','U', NALK, HMAT, NBAT, 0.0D0, 1.0D0, 1, NOC, ABSTOL, NFOUND, VALP, &
                     VT, NBAT, ISUPPZ, WORK, LWORK, IWORK, ILWORK, INFO )
         IF (INFO.NE.0) WRITE(*, *) 'WARNING - INFO=', INFO,' IN DSYEVR'
! end new DJW

         VALP(1:NOC)=-VALP(1:NOC)
      END IF
   END IF
   
   
   VALP(1:NOC)=-VALP(1:NOC)
! Total energy
   VGR=0.0D0
   EN=SUM(VALP(1:NDOC)) * 2.0D0
   IF (MONOC) EN=EN+VALP(NOC)
! Electron density
   IND=0
   E0=EN
! Calculation of the energy gradient
! At the end GRAD contains the derivatives in order:
! (DE/DXI, I=1, NBAT),(DE/DYI, I=1, NBAT),(DE/DZI, I=1, NBAT)
   IF (NALK > 0) then
      DO K=1, 3*NBAT
         GRD=0.D0
         DO I=1, NDOC
            DO L=1, NALK
               DO M=1, NALK
                  IF (L > M) THEN
                     ML=L*(L-1)/2 + M
                  ELSE
                     ML=M*(M-1)/2 + L
                  END IF
                  GRD=GRD+2.D0*VT(L, I)*VT(M, I)*HDER(K, ML)
               END DO
            END DO
         END DO
         IF (MONOC) THEN 
            DO L=1, NALK
               DO M=1, NALK
                  IF (L > M) THEN
                     ML=L*(L-1)/2+M
                  ELSE
                     ML=M*(M-1)/2+L
                  END IF
                  GRD=GRD + VT(L, NOC)*VT(M, NOC)*HDER(K, ML)
               END DO
            END DO
         END IF
         GRAD(K)=GRD
      END DO
   END IF

   IJ=0
   EREP=0.0
   DO I=1, NBAT
      DO J=1, I-1
         IJ=IJ+1
         EREP=EREP + 1.0D0/(R(IJ)**12)
      END DO
      IJ=IJ+1
   END DO

   EREP=EREPMAX*EREP

   DO I=1, NBAT
      DP(3*I-2)=GRAD(I)
      DP(3*I-1)=GRAD(I+NBAT)
      DP(3*I  )=GRAD(I+2*NBAT)
   END DO

   DO I=1, NBAT
      DO J=1, NBAT
         IF (I.NE.J) THEN
            XIJ=P(3*I-2)-P(3*J-2)
            YIJ=P(3*I-1)-P(3*J-1)
            ZIJ=P(3*I)-P(3*J)
            RIJ2=XIJ**2+YIJ**2+ZIJ**2
            RIJ14=RIJ2**7
            DP(3*I-2)=DP(3*I-2)-EREPMAX*12.*XIJ/RIJ14
            DP(3*I-1)=DP(3*I-1)-EREPMAX*12.*YIJ/RIJ14
            DP(3*I  )=DP(3*I  )-EREPMAX*12.*ZIJ/RIJ14
         END IF
      END DO
   END DO

   E0=E0+EREP

END SUBROUTINE ENTOTS

SUBROUTINE FTSS(R, F, DF)
   IMPLICIT NONE
   DOUBLE PRECISION A(10), B(10), C(10)
   DOUBLE PRECISION AA, BB, R, FEXP, FI, DFI, F, DF, XN
   INTEGER I, NT
   CHARACTER(LEN=2) METAL
   COMMON/METAL/METAL

   IF (METAL == 'NA') THEN
      A(1)=-0.0000600345
      A(2)=0.7366633997
      A(3)=-21.2536277450
      B(1)=7.8227026687
      B(2)=5.0784652072
      B(3)=-5.7014389921
      C(1)=1.4334140314
      C(2)=2.7260385637
      C(3)=0.0D0
      NT=3
      F=0.0D0
      DF=0.0D0
      DO I=1, NT
         FI=A(I)*(R**B(I))*EXP(-R*C(I))
         DFI=(-C(I)+B(I)/R)*FI
         F=F+FI
         DF=DF+DFI
      END DO
   ELSE IF (METAL == 'AG') THEN
      A(1)=-0.0148920093/27.21
      B(1)=9.4301889740
      C(1)=2.2599701079
      A(2)=-11.8760873508/27.21
      B(2)=-1.0347495105
      C(2)=0.5509616893
      NT=2
      F=0.0D0
      DF=0.0D0
      DO I=1, NT
         FI=A(I)*(R**B(I))*EXP(-R*C(I))
         DFI=(-C(I)+B(I)/R)*FI
         F=F+FI
         DF=DF+DFI
      END DO
   ELSE IF (METAL == 'LI') THEN
      AA=-0.000195D0
      BB=1.65D0
      XN=8.0D0
      FEXP=AA*EXP(-BB*R)
      F=(R**XN)*FEXP
      DF=(XN*R**(XN-1.D0)-BB*(R**XN))*FEXP
   END IF

END SUBROUTINE FTSS

SUBROUTINE FTSZ(R, F, DF)
   IMPLICIT NONE
   DOUBLE PRECISION A(10), B(10), C(10)
   DOUBLE PRECISION R, F, DF, AA, BB, DFI, FEXP, FI
   DOUBLE PRECISION XN
   INTEGER I, NT
   CHARACTER(LEN=2) METAL
   COMMON/METAL/METAL

   IF (METAL == 'NA') THEN
      A(1)=0.0002818237
      A(2)=-0.0025316896
      A(3)=67.9894048289
      B(1)=9.2162096989
      B(2)= 4.0958902363
      B(3)= -5.2076863938
      C(1)= 2.2554055362
      C(2)=.9408906266
      C(3)=.6157151927
      NT=3
      F=0.D0
      DF=0.D0
      DO I=1, NT
         FI=A(I)*(R**B(I))*EXP(-R*C(I))
         DFI=(-C(I)+B(I)/R)*FI
         F=F+FI
         DF=DF+DFI
      END DO
   ELSE IF (METAL == 'AG') THEN
      NT=2
      A(1)= -.0058542/27.21
      B(1)= 9.5986911821
      C(1)= 2.0836237596
      A(2)= -17.9409106/27.219
      B(2)=-3.5128659875
      C(2)=-.0000000007
      A(1)=-.0043951845/27.21
      A(2)=-13.469637022/27.219
      F=0.D0
      DF=0.D0
      DO I=1, NT
         FI=A(I)*(R**B(I))*EXP(-R*C(I))
         DFI=(-C(I)+B(I)/R)*FI
         F=F+FI
         DF=DF+DFI
      END DO
   ELSE IF (METAL == 'LI') THEN
      AA=0.0000091D0
      BB=1.80D0
      XN=10.D0
      FEXP=AA*EXP(-BB*R)
      F=(R**XN)*FEXP
      DF=(XN*R**(XN-1.D0)-BB*(R**XN))*FEXP
   END IF

END SUBROUTINE FTSZ

SUBROUTINE FRSS(R, F, DF)
   IMPLICIT NONE
   DOUBLE PRECISION R, F, DF
   DOUBLE PRECISION A, C, RE, C6, X6, FREP, DREP, FVDW, DVDW, ARG
   DOUBLE PRECISION COUP, DCOUP, DARG
   CHARACTER(LEN=2) METAL
   COMMON/METAL/METAL

   IF (METAL == 'NA') THEN
      A=1.6822175236
      C= 1.3007271463
      RE=5.5425217957
      C6= 10.4276970395
      X6= 5.8999900152
      FREP=A*EXP(-C*R)
      DREP=-C*FREP
      FVDW=C6*R**(-X6)
      DVDW=-X6*FVDW/R
      ARG= (RE/R-1.D0)
      IF (R > RE) THEN
         COUP=1.D0
         DCOUP=0.D0
      ELSE
         COUP=EXP(-ARG**2)
         DARG=-RE/R**2
         DCOUP=-2.D0*COUP*ARG*DARG
      END IF
      F=FREP-COUP*FVDW
      DF=DREP-DVDW*COUP-DCOUP*FVDW
   ELSE IF (METAL == 'AG') THEN
      A=8245.404377/27.21
      C=2.292537
      RE=16.726372
      C6=10679.191179/27.21
      X6=6. 
      FREP=A*EXP(-C*R)
      DREP=-C*FREP
      FVDW=C6*R**(-X6)
      DVDW=-X6*FVDW/R
      ARG= (RE/R-1.D0)
      IF (R > RE) THEN
         COUP=1.D0
         DCOUP=0.D0
      ELSE
         COUP=EXP(-ARG**2)
         DARG=-RE/R**2
         DCOUP=-2.D0*COUP*ARG*DARG
      END IF
      F=FREP-COUP*FVDW
      DF=DREP-DVDW*COUP-DCOUP*FVDW
   ELSE IF (METAL == 'LI') THEN
      A=1.6822175236D0
      C= 1.3007271463D0*5.80D0/5.05D0
      RE=5.5425217957D0*5.05D0/5.80D0
      X6= 5.8999900152D0
      C6= 10.4276970395D0*(5.05D0/5.80D0)**X6
      FREP=A*EXP(-C*R)
      DREP=-C*FREP
      FVDW=C6*R**(-X6)
      DVDW=-X6*FVDW/R
      ARG= (RE/R-1.D0)
      IF (R > RE) THEN
         COUP=1.D0
         DCOUP=0.D0
      ELSE
         COUP=EXP(-ARG**2)
         DARG=-RE/R**2
         DCOUP=-2.D0*COUP*ARG*DARG
      END IF
      F=1.43D0*(FREP-COUP*FVDW)
      DF=1.43D0*(DREP-DVDW*COUP-DCOUP*FVDW)
   END IF

END SUBROUTINE FRSS
