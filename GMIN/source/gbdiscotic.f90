      SUBROUTINE GBDISCOTIC (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY : NATOMS, GBSIGNOT, GBEPSNOT, GBMU, GBNU, SIGMAF, INVKAP, &
                          GBCHI, GBCHIPRM

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: SCSIG, SCSIG3, ENERGY
      DOUBLE PRECISION :: EPS1, EPS2, EPS
      DOUBLE PRECISION :: FCT1, FCT2, FCT3, FCT4, FCT5, FCT6, FCT7, FCT8
      DOUBLE PRECISION :: FCT9, FCT10, FCT11, FCT12, FCT13, FCT14, FCT15
      DOUBLE PRECISION :: FCT16, FCT17, FCT18, FCT19, FCT20
      DOUBLE PRECISION :: FCT3P4, FCT3M4, FCT7P8, FCT7M8
      DOUBLE PRECISION :: ALP, BET, GAM, APB, AMB
      DOUBLE PRECISION :: RIJSQ, R2, ABSRIJ, INVR, SCR 
      DOUBLE PRECISION :: SRM1, SRM2, SRM6, SRM7, SRM12, SRM13, SR12M6
      DOUBLE PRECISION :: VR, VA, VB, VG, FIJ(3), FIJN, FIJEI, FIJEJ 
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), P(3), EI(3), EJ(3), OST(3), RMI(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST
      DOUBLE PRECISION :: F, DF

      OST        = (/0.D0, 0.D0, 1.D0/)
      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      ENERGY = 0.D0
      IF (GTEST) G(:)   = 0.D0

      DO J1 = 1, REALNATOMS 

         J3      = 3*J1
         J5      = OFFSET + J3
         RI(:)   = X(J3-2:J3)
         P(:)    = X(J5-2:J5)

!     ROTATION MATRIX AND ITS DERIVATIVES

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         E(J1,:)     = MATMUL(RMI(:,:),OST(:))
         IF (GTEST) THEN
            DE1(J1,:)   = MATMUL(DRMI1(:,:),OST(:))
            DE2(J1,:)   = MATMUL(DRMI2(:,:),OST(:))
            DE3(J1,:)   = MATMUL(DRMI3(:,:),OST(:))
         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS - 1

         J3      = 3*J1
         J5      = OFFSET + J3
         RI(:)   = X(J3-2:J3)

         DO J2 = J1 + 1, REALNATOMS 

            J4    = 3*J2
            J6    = OFFSET + J4
            RJ(:) = X(J4-2:J4) 
            RIJ(:)  = RI(:) - RJ(:)
            RIJSQ   = DOT_PRODUCT(RIJ,RIJ)
            ABSRIJ  = DSQRT(RIJSQ)
            NR(:)   = RIJ(:)/ABSRIJ
            SCR     = ABSRIJ/SIGMAF
            INVR    = 1.D0/ABSRIJ
            R2      = 1.D0/RIJSQ

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

            EI  = E(J1,:)
            EJ  = E(J2,:)
            ALP = DOT_PRODUCT(NR(:),EI(:)) 
            BET = DOT_PRODUCT(NR(:),EJ(:))
            GAM = DOT_PRODUCT(EI(:),EJ(:))

!     CALCULATE USEFUL QUANTITIES

            APB    = ALP+BET
            AMB    = ALP-BET

            FCT1   = 1.D0/(1.D0+GBCHI*GAM)
            FCT2   = 1.D0/(1.D0-GBCHI*GAM)
            FCT3   = APB*FCT1
            FCT4   = AMB*FCT2
            FCT3P4 = FCT3+FCT4
            FCT3M4 = FCT3-FCT4

            FCT5   = 1.D0/(1.D0+GBCHIPRM*GAM)
            FCT6   = 1.D0/(1.D0-GBCHIPRM*GAM)
            FCT7   = (ALP+BET)*FCT5
            FCT8   = (ALP-BET)*FCT6
            FCT7P8 = FCT7+FCT8
            FCT7M8 = FCT7-FCT8

!     CALCULATE $\epsilon$

            EPS1   = DSQRT(FCT1*FCT2)
            EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
            EPS    = GBEPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

            SCSIG  = 1.D0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
            SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE DEL(V)/DEL(R)

            SRM1   = 1.D0/(SCR - SCSIG*INVKAP +1.D0)
            SRM2   = SRM1*SRM1
            SRM6   = SRM2*SRM2*SRM2
            SRM12  = SRM6*SRM6
            SR12M6 = SRM12-SRM6
            FCT9   = 2.D0*SRM13-SRM7
!            VR     = -(24.D0/SIGMAF)*EPS*FCT9

!     CALCULATE ENERGY

            ENERGY = ENERGY + EPS*SR12M6

            IF (GTEST) THEN
               
!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)
               
               SRM7   = SRM6*SRM1
               SRM13  = SRM12*SRM1
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/SIGMAF)*EPS*FCT9   
               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(GBSIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*INVKAP*SCSIG3
               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*INVKAP*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9
               VG     = FCT19 - FCT20

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FIJ  = FIJN*NR + FIJEI*EI + FIJEJ*EJ

               DADPI1 = DOT_PRODUCT(NR(:),DE1(J1,:))
               DADPI2 = DOT_PRODUCT(NR(:),DE2(J1,:))
               DADPI3 = DOT_PRODUCT(NR(:),DE3(J1,:))

               DBDPJ1 = DOT_PRODUCT(NR(:),DE1(J2,:))
               DBDPJ2 = DOT_PRODUCT(NR(:),DE2(J2,:))
               DBDPJ3 = DOT_PRODUCT(NR(:),DE3(J2,:))

               DGDPI1 = DOT_PRODUCT(DE1(J1,:),EJ(:))
               DGDPI2 = DOT_PRODUCT(DE2(J1,:),EJ(:))
               DGDPI3 = DOT_PRODUCT(DE3(J1,:),EJ(:))

               DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J2,:))
               DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J2,:))
               DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J2,:))

               G(J3-2:J3) = G(J3-2:J3) + FIJ 
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

            ENDIF

         ENDDO

      ENDDO

      ENERGY = 4.D0*ENERGY

      END SUBROUTINE GBDISCOTIC 
