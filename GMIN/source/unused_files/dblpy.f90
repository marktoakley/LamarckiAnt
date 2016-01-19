      SUBROUTINE DBLPY (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, PYA1, PYA2, PYSIGNOT, PYEPSNOT, RADIFT, CLJ6, CLJ12 

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: AEZR1(3,3), AEZR2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: DEZ(3,3), D1E(3,3), D2E(3,3), D3E(3,3) 
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: SE1(3*NATOMS/2,3), SD1E1(3*NATOMS/2,3), SD2E1(3*NATOMS/2,3), SD3E1(3*NATOMS/2,3)
      DOUBLE PRECISION :: SE2(3*NATOMS/2,3), SD1E2(3*NATOMS/2,3), SD2E2(3*NATOMS/2,3), SD3E2(3*NATOMS/2,3)
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, EI(3), EJ(3), DU(3)
      LOGICAL          :: GTEST

      ENERGY = 0.D0
      IF (GTEST) G(:)   = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      CALL DEFDUMBBELL(CLJ6, CLJ12)

      AEZR1(:,:) = 0.D0
      AEZR2(:,:) = 0.D0

      DO J1 = 1, 3

         AEZR1(J1,J1) = 1.D0/(PYA1(J1)*PYA1(J1))
         AEZR2(J1,J1) = 1.D0/(PYA2(J1)*PYA2(J1))

      ENDDO

      DO J1 = 1, REALNATOMS

         J3      = 3*J1
         J5      = OFFSET + J3
         RI      = X(J3-2:J3)
         P       = X(J5-2:J5)

!     ROTATION MATRIX AND ITS DERIVATIVES

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))
            IF (GTEST) THEN
               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))
            ENDIF

         ENDDO

         SE1(J3-2:J3,:) = MATMUL(RMI(:,:),(MATMUL(AEZR1(:,:),(TRANSPOSE(RMI(:,:))))))

         IF (RADIFT) THEN

            SE2(J3-2:J3,:)   = MATMUL(RMI(:,:),(MATMUL(AEZR2(:,:),(TRANSPOSE(RMI(:,:))))))

         ENDIF

         IF (GTEST) THEN

            DEZ(:,:)         = MATMUL(DRMI1(:,:),AEZR1(:,:))
            SD1E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI2(:,:),AEZR1(:,:))
            SD2E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI3(:,:),AEZR1(:,:))
            SD3E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))

            IF (RADIFT) THEN

            DEZ(:,:)         = MATMUL(DRMI1(:,:),AEZR2(:,:))
            SD1E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI2(:,:),AEZR2(:,:))
            SD2E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI3(:,:),AEZR2(:,:))
            SD3E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
 
            ENDIF

         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS - 1

         J3    = 3*J1
         J5    = OFFSET + J3
         RI    = X(J3-2:J3)

         AE1(:,:) = SE1(J3-2:J3,:)

         IF (RADIFT) THEN

            AE2(:,:) = SE2(J3-2:J3,:)         

         ENDIF

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4
            RJ = X(J4-2:J4) 

!     SUM OVER LJ SITES
            DO I = 1, NRBSITES

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8       = NRBSITES*(J2-1) + J
                  RIJSS(:) = R(J7,:) - R(J8,:)
                  DSS2     = DOT_PRODUCT(RIJSS,RIJSS)
                  R2       = 1.D0/DSS2
                  R6       = R2*R2*R2

                  ENERGY = ENERGY + (CLJ12*R6 - CLJ6)*R6

                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     DVDR   = -6.D0*(2.D0*CLJ12*R6 - CLJ6)*R6*R2

                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RIJSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RIJSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RIJSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RIJSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RIJSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RIJSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

            BE1(:,:) = SE1(J4-2:J4,:)

            IF (RADIFT) THEN
   
               BE2(:,:) = SE2(J4-2:J4,:)

            ENDIF

            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ

!     CALCULATE ECF

            CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

            FCNT1   = -FMIN
            SRTFI1  = 1.D0 / DSQRT(FCNT1)

            RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
            RHO1SQ = RHO1*RHO1
            RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
            RHO112 = RHO16 * RHO16

            IF (RADIFT) THEN

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

               FCNT2    = -FMIN
               SRTFI2   = 1.D0/DSQRT(FCNT2)
               RHO2     = PYSIGNOT/(ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
               RHO2SQ   = RHO2*RHO2
               RHO26    = RHO2SQ*RHO2SQ*RHO2SQ

            ELSE

               RHO2   = RHO1
               RHO26  = RHO16

            ENDIF

            ENERGY = ENERGY + 4.D0*PYEPSNOT*(RHO112 - RHO26)

            IF (GTEST) THEN

               APB(:,:) = LAMDAC1*AE1(:,:) + (1.D0-LAMDAC1)*BE1(:,:)

               CALL MTRXIN (APB(:,:), APBINV(:,:))

               ARIBRJ(:) = LAMDAC1*MATMUL(AE1(:,:),RI(:)) + (1.D0-LAMDAC1)*MATMUL(BE1(:,:),RJ(:))
               XC(:)     = MATMUL(APBINV(:,:), ARIBRJ(:))
               XCMRI(:)  = XC(:) - RI(:)
               XCMRJ(:)  = XC(:) - RJ(:)
               DF1DR(:)  = -2.D0*LAMDAC1*MATMUL(AE1(:,:),XCMRI(:))
               FCTR1     = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR(:)  = (1.D0-SRTFI1)*NR(:)/PYSIGNOT + FCTR1*DF1DR(:)
               DVDF1     = -2.D0*RHO112*RHO1*FCTR1
                
               D1E(:,:)  = SD1E1(J3-2:J3,:)
               D2E(:,:)  = SD2E1(J3-2:J3,:)
               D3E(:,:)  = SD3E1(J3-2:J3,:)

               DF1PI1    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D1E(:,:),XCMRI(:)))
               DF1PI2    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D2E(:,:),XCMRI(:)))
               DF1PI3    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D3E(:,:),XCMRI(:)))

               D1E(:,:)  = SD1E1(J4-2:J4,:)
               D2E(:,:)  = SD2E1(J4-2:J4,:)
               D3E(:,:)  = SD3E1(J4-2:J4,:)
 
               DF1PJ1    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D1E(:,:),XCMRJ(:)))
               DF1PJ2    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D2E(:,:),XCMRJ(:)))
               DF1PJ3    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D3E(:,:),XCMRJ(:)))

               IF (RADIFT) THEN

                  APB(:,:) = LAMDAC2*AE2(:,:) + (1.D0-LAMDAC2)*BE2(:,:)

                  CALL MTRXIN (APB(:,:), APBINV(:,:))

                  ARIBRJ(:) = LAMDAC2*MATMUL(AE2(:,:),RI(:)) + (1.D0-LAMDAC2)*MATMUL(BE2(:,:),RJ(:))
                  XC(:)     = MATMUL(APBINV(:,:), ARIBRJ(:))
                  XCMRI(:)  = XC(:) - RI(:)
                  XCMRJ(:)  = XC(:) - RJ(:)
                  DF2DR(:)  = -2.D0*LAMDAC2*MATMUL(AE2(:,:),XCMRI(:))
                 
                  FCTR2     = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR(:)  = (1.D0-SRTFI2)*NR(:)/PYSIGNOT + FCTR2*DF2DR(:)
                  DVDF2     = RHO26*RHO2*FCTR2

                  D1E(:,:)  = SD1E2(J3-2:J3,:)
                  D2E(:,:)  = SD2E2(J3-2:J3,:)
                  D3E(:,:)  = SD3E2(J3-2:J3,:)

                  DF2PI1    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D1E(:,:),XCMRI(:)))
                  DF2PI2    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D2E(:,:),XCMRI(:)))
                  DF2PI3    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D3E(:,:),XCMRI(:)))

                  D1E(:,:)  = SD1E2(J4-2:J4,:)
                  D2E(:,:)  = SD2E2(J4-2:J4,:)
                  D3E(:,:)  = SD3E2(J4-2:J4,:)

                  DF2PJ1    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D1E(:,:),XCMRJ(:)))
                  DF2PJ2    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D2E(:,:),XCMRJ(:)))
                  DF2PJ3    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D3E(:,:),XCMRJ(:)))

               ELSE

                  DG2DR(:)  = DG1DR(:)
                  DVDF2     = RHO26*RHO2*FCTR1
                  DF2PI1    = DF1PI1
                  DF2PI2    = DF1PI2
                  DF2PI3    = DF1PI3
                  DF2PJ1    = DF1PJ1
                  DF2PJ2    = DF1PJ2
                  DF2PJ3    = DF1PJ3

               ENDIF

!     CALCULATE GRADIENT

               FIJ(:) = 2.D0*RHO112*RHO1*DG1DR(:) - RHO26*RHO2*DG2DR(:)
               TIJ(1) = DVDF1*DF1PI1 + DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + DVDF2*DF2PJ3

               G(J3-2:J3) = G(J3-2:J3) - 24.D0*PYEPSNOT*FIJ(:)
               G(J4-2:J4) = G(J4-2:J4) + 24.D0*PYEPSNOT*FIJ(:)
               G(J5-2:J5) = G(J5-2:J5) + 24.D0*PYEPSNOT*TIJ(:)
               G(J6-2:J6) = G(J6-2:J6) + 24.D0*PYEPSNOT*TJI(:)

            ENDIF

         ENDDO

      ENDDO

      END SUBROUTINE DBLPY 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDUMBBELL(CLJ6, CLJ12)

      USE COMMONS, ONLY: NRBSITES, SITE, PYA1

      IMPLICIT NONE

      DOUBLE PRECISION :: CLJ6, CLJ12

      SITE(:,:) = 0.D0
      SITE(1,3) = PYA1(3)
      SITE(2,3) =-PYA1(3)

      CLJ6  = 4.D0*1.D0**6.D0
      CLJ12 = 4.D0*1.D0**12.D0

      END SUBROUTINE DEFDUMBBELL
