      SUBROUTINE GBDDP (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY : NATOMS, GBSIGNOT, GBEPSNOT, GBMU, GBNU, SIGMAF, INVKAP, GBCHI, GBCHIPRM, GBDPFCT 

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET
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
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), P(3), OST(3), EI(3), EJ(3), RMI(3,3), RMJ(3,3), DU(3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DEDPI1(3), DEDPI2(3), DEDPI3(3), DEDPJ1(3), DEDPJ2(3), DEDPJ3(3)
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST

      OST     = (/0.D0, 0.D0, 1.D0/)
      DU      = (/1.D0, 1.D0, 1.D0/)/DSQRT(3.D0) 

!     FROM INPUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      IF (GTEST) THEN

         ENERGY = 0.D0
         G(:)   = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST) 

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT(P, RMJ, DRMJ1, DRMJ2, DRMJ3, GTEST) 
               
!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ
               SCR    = ABSRIJ/SIGMAF
               INVR   = 1.D0/ABSRIJ
               R2     = 1.D0/RIJSQ

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

               EI  = MATMUL(RMI,OST)
               EJ  = MATMUL(RMJ,OST)
               ALP = DOT_PRODUCT(NR,EI) 
               BET = DOT_PRODUCT(NR,EJ)
               GAM = DOT_PRODUCT(EI,EJ)

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

               SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
               SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE DEL(V)/DEL(R)

               SRM1   = 1.D0/(SCR - SCSIG*INVKAP +1.D0)
               SRM2   = SRM1*SRM1
               SRM6   = SRM2*SRM2*SRM2
               SRM7   = SRM6*SRM1
               SRM12  = SRM6*SRM6
               SRM13  = SRM12*SRM1
               SR12M6 = SRM12-SRM6
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/SIGMAF)*EPS*FCT9

!     CALCULATE ENERGY

               ENERGY = ENERGY + 4.D0*EPS*SR12M6

!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)

               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(GBSIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*INVKAP*SCSIG3

               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
!                       /EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*INVKAP*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9

               VG     = FCT19 - FCT20

!     CALCULATE CONTRIBUTION TO FORCES

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FIJ  = FIJN*NR + FIJEI*EI + FIJEJ*EJ

               DEDPI1 = MATMUL(DRMI1,OST)
               DEDPI2 = MATMUL(DRMI2,OST)
               DEDPI3 = MATMUL(DRMI3,OST)

               DEDPJ1 = MATMUL(DRMJ1,OST)
               DEDPJ2 = MATMUL(DRMJ2,OST)
               DEDPJ3 = MATMUL(DRMJ3,OST)
 
               DADPI1 = DOT_PRODUCT(NR,DEDPI1)  
               DADPI2 = DOT_PRODUCT(NR,DEDPI2)  
               DADPI3 = DOT_PRODUCT(NR,DEDPI3)  
                                                           
               DBDPJ1 = DOT_PRODUCT(NR,DEDPJ1)  
               DBDPJ2 = DOT_PRODUCT(NR,DEDPJ2)  
               DBDPJ3 = DOT_PRODUCT(NR,DEDPJ3)

               DGDPI1 = DOT_PRODUCT(DEDPI1,EJ)
               DGDPI2 = DOT_PRODUCT(DEDPI2,EJ)
               DGDPI3 = DOT_PRODUCT(DEDPI3,EJ)

               DGDPJ1 = DOT_PRODUCT(EI,DEDPJ1)
               DGDPJ2 = DOT_PRODUCT(EI,DEDPJ2)
               DGDPJ3 = DOT_PRODUCT(EI,DEDPJ3)  

               G(J3-2:J3) = G(J3-2:J3) + FIJ 
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

!     CONTRIBUTION FROM DIPOLAR INTERACTIONS

               EI = MATMUL(RMI,DU)
               EJ = MATMUL(RMJ,DU)
               ALP = DOT_PRODUCT(NR,EI)
               BET = DOT_PRODUCT(NR,EJ)
               GAM = DOT_PRODUCT(EI,EJ)

               VR  = -GBDPFCT*R2*R2*(GAM - 3.D0*ALP*BET)
               VA  = -GBDPFCT*BET*R2*INVR
               VB  = -GBDPFCT*ALP*R2*INVR
               VG  = GBDPFCT*R2*INVR/3.D0

               ENERGY = ENERGY + GBDPFCT*INVR*R2*(GAM/3.D0 - ALP*BET)

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR
               FIJ    = FIJN*NR + FIJEI*EI + FIJEJ*EJ

               DEDPI1 = MATMUL(DRMI1,DU)
               DEDPI2 = MATMUL(DRMI2,DU)
               DEDPI3 = MATMUL(DRMI3,DU)

               DEDPJ1 = MATMUL(DRMJ1,DU)
               DEDPJ2 = MATMUL(DRMJ2,DU)
               DEDPJ3 = MATMUL(DRMJ3,DU)

               DADPI1 = DOT_PRODUCT(NR,DEDPI1)
               DADPI2 = DOT_PRODUCT(NR,DEDPI2)
               DADPI3 = DOT_PRODUCT(NR,DEDPI3)

               DBDPJ1 = DOT_PRODUCT(NR,DEDPJ1)
               DBDPJ2 = DOT_PRODUCT(NR,DEDPJ2)
               DBDPJ3 = DOT_PRODUCT(NR,DEDPJ3)

               DGDPI1 = DOT_PRODUCT(DEDPI1,EJ)
               DGDPI2 = DOT_PRODUCT(DEDPI2,EJ)
               DGDPI3 = DOT_PRODUCT(DEDPI3,EJ)

               DGDPJ1 = DOT_PRODUCT(EI,DEDPJ1)
               DGDPJ2 = DOT_PRODUCT(EI,DEDPJ2)
               DGDPJ3 = DOT_PRODUCT(EI,DEDPJ3)

               G(J3-2:J3) = G(J3-2:J3) + FIJ
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

            ENDDO

         ENDDO

      ENDIF

      END SUBROUTINE GBDDP 
