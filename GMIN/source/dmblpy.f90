      SUBROUTINE DUMBBELLPOLARYUKAWA (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, DBPMU, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, DIST, DR, R2, ABSRIJ, RIJSQ, DPFCT, SIGMA, EXPFCT
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RSS(3), NR(3), P(3), DU(3), EI(3), EJ(3), R(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: EPS, KAPPA, SIGAA, SIGBB, SIGAB
      DOUBLE PRECISION :: GLJEXP, CLJNAA, CLJ2NAA, CLJNBB, CLJ2NBB, CLJNAB, CLJ2NAB, LJRN, LJR2N, CLJN, CLJ2N
      LOGICAL          :: GTEST
     
      CALL DEFDMBLYUKAWA(DU, DPFCT, EPS, KAPPA, SIGAA, SIGBB, SIGAB) 

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         E(J1,:)   = MATMUL(RMI(:,:),DU(:))
         DE1(J1,:) = MATMUL(DRMI1(:,:),DU(:))
         DE2(J1,:) = MATMUL(DRMI2(:,:),DU(:))
         DE3(J1,:) = MATMUL(DRMI3(:,:),DU(:))

         DO J2 = 1, NRBSITES - 1

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))

            IF (GTEST) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS 

         J3 = 3*J1
         J5 = OFFSET + J3
 
         RI(:)  = X(J3-2:J3)
         EI(:)  = E(J1,:)

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     LJ CONTRIBUTION

            DO I = 1, NRBSITES - 1

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES - 1

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  R2     = 1.D0/R2 
                  IF (I == 1 .AND. J == 1) THEN
                     SIGMA = SIGAA
                  ELSEIF (I == 2 .AND. J == 2) THEN
                     SIGMA = SIGBB
                  ELSE
                     SIGMA = SIGAB
                  ENDIF
 
                  EXPFCT = EXP(-KAPPA*(ABSRIJ - SIGMA))
                  ENERGY = ENERGY + EPS*SIGMA*EXPFCT/ABSRIJ 

                  IF (GTEST) THEN

                     DVDR =- EPS*SIGMA*EXPFCT*R2*(KAPPA + 1.D0/ABSRIJ)

                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

!     CONTRIBUTION FROM DIPOLAR INTERACTIONS

            RJ(:)  = X(J4-2:J4)
            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ
            EJ(:)  = E(J2,:)
            ALP    = DOT_PRODUCT(NR(:),EI(:))
            BET    = DOT_PRODUCT(NR(:),EJ(:))
            GAM    = DOT_PRODUCT(EI(:),EJ(:))

            ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

            IF (GTEST) THEN

               VR  = -DPFCT*R2*R2*(GAM - 3.D0*ALP*BET)
               VA  = -DPFCT*BET*R2/ABSRIJ
               VB  = -DPFCT*ALP*R2/ABSRIJ
               VG  =  DPFCT*R2/(3.D0*ABSRIJ)

               FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
               FIJEI  = VA/ABSRIJ
               FIJEJ  = VB/ABSRIJ
               FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

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

               G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

            ENDIF

         ENDDO

         ENERGY = ENERGY - DBPMU*EFIELD*EI(3)

         IF (GTEST) THEN

            G(J5-2) = G(J5-2) - DBPMU*EFIELD*DE1(J1,3)
            G(J5-1) = G(J5-1) - DBPMU*EFIELD*DE2(J1,3)
            G(J5)   = G(J5)   - DBPMU*EFIELD*DE3(J1,3)

         ENDIF  

      ENDDO

      END SUBROUTINE DUMBBELLPOLARYUKAWA 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDMBLYUKAWA(DU, DPFCT, EPS, KAPPA, SIGAA, SIGBB, SIGAB)
    
      USE COMMONS, ONLY: SITE, DBSIGBB, YKAPPA, YEPS, DBPMU  
     
      IMPLICIT NONE

      DOUBLE PRECISION :: EPS, KAPPA, SIGAA, SIGBB, SIGAB, DPFCT, DU(3), GLJEXP, CLJNAA, CLJ2NAA, CLJNBB, CLJ2NBB, CLJNAB, CLJ2NAB

      EPS   = YEPS
      KAPPA = YKAPPA         
      SIGAA = 1.D0
      SIGBB = DBSIGBB
      SIGAB = 0.5D0*(SIGAA + SIGBB)        

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.5D0*SIGAA

      SITE(2,1) = 0.D0
      SITE(2,2) = 0.D0
      SITE(2,3) = -0.5D0*DBSIGBB

      SITE(3,1) = 0.D0
      SITE(3,2) = 0.D0
      SITE(3,3) = 0.D0

      DU(:) = (/0.D0, 1.D0, 0.D0/)
      DPFCT = 3.D0*DBPMU*DBPMU*SIGAA**3.D0

      END SUBROUTINE DEFDMBLYUKAWA
