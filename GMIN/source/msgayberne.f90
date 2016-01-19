      SUBROUTINE MSGAYBERNE (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY : NATOMS, NGBSITE, GBMU, GBNU, GBCHI, GBCHIPRM, SIGNOT, EPSNOT

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J4, J5, J6, INDX, INDXI, INDXJ, K1, K2, NMOL, OFFSET
      DOUBLE PRECISION :: SCSIG, SCSIG3, EPS1, EPS2, EPS
      DOUBLE PRECISION :: FCT1, FCT2, FCT3, FCT4, FCT5, FCT6, FCT7, FCT8
      DOUBLE PRECISION :: FCT9, FCT10, FCT11, FCT12, FCT13, FCT14, FCT15
      DOUBLE PRECISION :: FCT16, FCT17, FCT18, FCT19, FCT20 
      DOUBLE PRECISION :: FCT3P4, FCT3M4, FCT7P8, FCT7M8
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: PI(3), PJ(3), THETA, THETA2
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: RMJ(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: E(3,3), DE1(3,3), DE2(3,3), DE3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, APB, AMB 
      DOUBLE PRECISION :: R, R2, INVR, SCR 
      DOUBLE PRECISION :: RI(3), RJ(3), NR(3), EI(3), EJ(3), DSS(3) 
      DOUBLE PRECISION :: SRM1, SRM2, SRM6, SRM7, SRM12, SRM13, SR12M6
      DOUBLE PRECISION :: VR, VA, VB, VG, ENERGY
      DOUBLE PRECISION :: FXIJ, FYIJ, FZIJ, FIAJBN, FIAJBI, FIAJBJ
      DOUBLE PRECISION :: D1PST(3), D2PST(3), D3PST(3), D1OST(3), D2OST(3), D3OST(3)
      DOUBLE PRECISION :: D1R, D1A, D1B, D1G, D2R, D2A, D2B, D2G, D3R, D3A, D3B, D3G
      DOUBLE PRECISION :: PST(NGBSITE,3), OST(NGBSITE,3)
      LOGICAL          :: GTEST

      NMOL       = NATOMS/2
      OFFSET     = 3*NMOL

      CALL DEFMOL(PST, OST)

      ENERGY     = 0.D0

      IF (GTEST) THEN

         G(:)  = 0.D0

         DO J1 = 1, NMOL - 1

            J3 = 3*J1
            J5 = OFFSET + J3
            RI = X(J3-2:J3)
            PI = X(J5-2:J5)

!           CALL ROTMAT (PI, RMI, DRMI1, DRMI2, DRMI3) ! wrong number of arguments. Commented by DJW 5/5/13
       
            DO I = 1, NGBSITE
        
               EI(:)  = MATMUL(RMI,OST(I,:))
         
               DO J2 =  J1 + 1, NMOL
 
                  J4      = 3*J2
                  J6      = OFFSET + J4
                  RJ      = X(J4-2:J4) 
                  PJ      = X(J6-2:J6)

!                 CALL ROTMAT (PI, RMI, DRMJ1, DRMJ2, DRMJ3) ! wrong number of arguments. Commented by DJW 5/5/13

                  DO J = 1, NGBSITE

                     EJ(:)  = MATMUL(RMJ,OST(J,:)) 
 
                     DSS(:)  = RI(:) - RJ(:) + MATMUL(RMI,PST(I,:)) - MATMUL(RMJ,PST(J,:))
                     R2      = DOT_PRODUCT(DSS,DSS)
                     R       = DSQRT(R2)
                     SCR     = R/SIGNOT
                     INVR    = 1.D0/R
                     R2      = 1.D0/R2

!     NORMILIZE THE SITE-SITE SEPRATION VECTOR

                     NR(:)   = DSS(:)*INVR

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

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
                     EPS    = EPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

                     SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
                     SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE $\PARTIAL(V)/\PARTIAL(R)$

                     SRM1   = 1.D0/(SCR-SCSIG+1.D0)
                     SRM2   = SRM1*SRM1
                     SRM6   = SRM2*SRM2*SRM2
                     SRM7   = SRM6*SRM1
                     SRM12  = SRM6*SRM6
                     SRM13  = SRM12*SRM1
                     SR12M6 = SRM12-SRM6
                     FCT9   = 2.D0*SRM13-SRM7
                     VR     = -(24.D0/SIGNOT)*EPS*FCT9

!     CALCULATE ENERGY

                     ENERGY = ENERGY + EPS*SR12M6

!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)

                     FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(SIGNOT*EPS2)
                     FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
                     FCT12  = 12.D0*EPS*GBCHI*SCSIG3

                     VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
                     VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

                     FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
                     FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8) / EPS2
                     FCT15  = FCT13 + FCT14
                     FCT16  = FCT3*FCT3 - FCT4*FCT4
                     FCT17  = 4.D0*EPS*FCT15
                     FCT18  = 6.D0*EPS*GBCHI*GBCHI*SCSIG3*FCT16
                     FCT19  = FCT17*SR12M6
                     FCT20  = FCT18*FCT9

                     VG     = FCT19 - FCT20

!     CALCULATE CONTRIBUTION TO FORCES

                     FIAJBN  = VR - (VA*ALP+VB*BET)*INVR
                     FIAJBI  = VA*INVR
                     FIAJBJ  = VB*INVR

                     FXIJ = FIAJBN*NR(1) + FIAJBI*EI(1) + FIAJBJ*EJ(1)
                     FYIJ = FIAJBN*NR(2) + FIAJBI*EI(2) + FIAJBJ*EJ(2)
                     FZIJ = FIAJBN*NR(3) + FIAJBI*EI(3) + FIAJBJ*EJ(3)

                     G(J3-2) = G(J3-2) + FXIJ
                     G(J3-1) = G(J3-1) + FYIJ
                     G(J3)   = G(J3)   + FZIJ

                     G(J4-2) = G(J4-2) - FXIJ
                     G(J4-1) = G(J4-1) - FYIJ
                     G(J4)   = G(J4)   - FZIJ

!     DERIVATIVES WITH RESPECT TO THE ANGULAR COORDINATES OF J1-TH MOLECULE            

                     D1PST = MATMUL(DRMI1,PST(I,:))
                     D2PST = MATMUL(DRMI2,PST(I,:))
                     D3PST = MATMUL(DRMI3,PST(I,:))

                     D1OST = MATMUL(DRMI1,OST(I,:))
                     D2OST = MATMUL(DRMI2,OST(I,:))
                     D3OST = MATMUL(DRMI3,OST(I,:))

                     D1R = DOT_PRODUCT(NR,D1PST)
                     D1A = (DOT_PRODUCT(EI,D1PST) - ALP*D1R)*INVR + DOT_PRODUCT(NR,D1OST)
                     D1B = (DOT_PRODUCT(EJ,D1PST) - BET*D1R)*INVR
                     D1G = DOT_PRODUCT(EJ,D1OST)

                     D2R = DOT_PRODUCT(NR,D2PST)
                     D2A = (DOT_PRODUCT(EI,D2PST) - ALP*D2R)*INVR + DOT_PRODUCT(NR,D2OST)
                     D2B = (DOT_PRODUCT(EJ,D2PST) - BET*D2R)*INVR
                     D2G = DOT_PRODUCT(EJ,D2OST)

                     D3R = DOT_PRODUCT(NR,D3PST)
                     D3A = (DOT_PRODUCT(EI,D3PST) - ALP*D3R)*INVR + DOT_PRODUCT(NR,D3OST)
                     D3B = (DOT_PRODUCT(EJ,D3PST) - BET*D3R)*INVR
                     D3G = DOT_PRODUCT(EJ,D3OST)

                     G(J5-2) = G(J5-2) + VR*D1R + VA*D1A + VB*D1B + VG*D1G
                     G(J5-1) = G(J5-1) + VR*D2R + VA*D2A + VB*D2B + VG*D2G
                     G(J5)   = G(J5)   + VR*D3R + VA*D3A + VB*D3B + VG*D3G

!     DERIVATIVES WITH RESPECT TO THE ANGULAR COORDINATES OF J2-TH MOLECULE

                     D1PST = -MATMUL(DRMJ1,PST(J,:))
                     D2PST = -MATMUL(DRMJ2,PST(J,:))
                     D3PST = -MATMUL(DRMJ3,PST(J,:))

                     D1OST = MATMUL(DRMJ1,OST(J,:))
                     D2OST = MATMUL(DRMJ2,OST(J,:))
                     D3OST = MATMUL(DRMJ3,OST(J,:))

                     D1R = DOT_PRODUCT(NR,D1PST)
                     D1A = (DOT_PRODUCT(EI,D1PST) - ALP*D1R)*INVR 
                     D1B = (DOT_PRODUCT(EJ,D1PST) - BET*D1R)*INVR + DOT_PRODUCT(NR,D1OST)
                     D1G = DOT_PRODUCT(EI,D1OST)

                     D2R = DOT_PRODUCT(NR,D2PST)
                     D2A = (DOT_PRODUCT(EI,D2PST) - ALP*D2R)*INVR 
                     D2B = (DOT_PRODUCT(EJ,D2PST) - BET*D2R)*INVR + DOT_PRODUCT(NR,D2OST)
                     D2G = DOT_PRODUCT(EI,D2OST)

                     D3R = DOT_PRODUCT(NR,D3PST)
                     D3A = (DOT_PRODUCT(EI,D3PST) - ALP*D3R)*INVR 
                     D3B = (DOT_PRODUCT(EJ,D3PST) - BET*D3R)*INVR + DOT_PRODUCT(NR,D3OST)
                     D3G = DOT_PRODUCT(EI,D3OST)

                     G(J6-2) = G(J6-2) + VR*D1R + VA*D1A + VB*D1B + VG*D1G
                     G(J6-1) = G(J6-1) + VR*D2R + VA*D2A + VB*D2B + VG*D2G
                     G(J6)   = G(J6)   + VR*D3R + VA*D3A + VB*D3B + VG*D3G

                  ENDDO

               ENDDO

            ENDDO

         ENDDO

         ENERGY = 4.D0 * ENERGY

      ENDIF

      RETURN
      END SUBROUTINE MSGAYBERNE 

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE DEFMOL(PST, OST)

      USE COMMONS, ONLY : NGBSITE, SIGNOT, GBKAPPA

      IMPLICIT NONE

      DOUBLE PRECISION :: PST(NGBSITE,3), OST(NGBSITE,3)

      PST(1,1) = - 0.25D0*SIGNOT*GBKAPPA
      PST(1,2) = 0.D0
      PST(1,3) = 0.D0

      PST(2,1) = 0.25D0*SIGNOT*GBKAPPA
      PST(2,2) = 0.D0
      PST(2,3) = 0.D0

      OST(1,1) = 1.D0
      OST(1,2) = 0.D0
      OST(1,3) = 0.D0
  
      OST(2,1) = 0.D0
      OST(2,2) = 1.D0
      OST(2,3) = 0.D0

      RETURN 
      END SUBROUTINE DEFMOL
