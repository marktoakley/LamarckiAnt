      SUBROUTINE MSTBIN(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, NRBSITES1, NPS, SITE, RBUV, STOCKMU, EFIELDT, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, NRBS1, NRBS2 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DVDR, DPFCT
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3)
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3), E(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: DE1(NATOMS*NRBSITES/2,3), DE2(NATOMS*NRBSITES/2,3), DE3(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      LOGICAL          :: GTEST

      DPFCT   = 3.D0*STOCKMU*STOCKMU

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

         IF (J1 <= NPS) THEN

            DO J2 = 1, NRBSITES1

               J4        = NRBSITES1*(J1-1) + J2
               R(J4,:)   = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
               E(J4,:)   = MATMUL(RMI(:,:),RBUV(J2,:))

               IF (GTEST) THEN

                 DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
                 DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
                 DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

                 DE1(J4,:) = MATMUL(DRMI1(:,:),RBUV(J2,:))
                 DE2(J4,:) = MATMUL(DRMI2(:,:),RBUV(J2,:))
                 DE3(J4,:) = MATMUL(DRMI3(:,:),RBUV(J2,:))

               ENDIF

            ENDDO

         ELSE
  
            DO J2 = NRBSITES1 + 1, NRBSITES

               J4        = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J1-NPS-1) + J2 - NRBSITES1
               R(J4,:)   = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
               E(J4,:)   = MATMUL(RMI(:,:),RBUV(J2,:))

               IF (GTEST) THEN

                  DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
                  DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
                  DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

                  DE1(J4,:) = MATMUL(DRMI1(:,:),RBUV(J2,:))
                  DE2(J4,:) = MATMUL(DRMI2(:,:),RBUV(J2,:))
                  DE3(J4,:) = MATMUL(DRMI3(:,:),RBUV(J2,:))

               ENDIF
            
            ENDDO

         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS  

         J3 = 3*J1
         J5 = OFFSET + J3

         RI(:)  = X(J3-2:J3)

         IF (J1 <= NPS) THEN
            NRBS1 = 1
            NRBS2 = NRBSITES1
         ELSE 
            NRBS1 = NRBSITES1 + 1
            NRBS2 = NRBSITES
         ENDIF
 
         DO I = NRBS1, NRBS2

            IF (J1 <= NPS) THEN
               J7 = NRBSITES1*(J1-1) + I
            ELSE
               J7 = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J1-NPS-1) + I - NRBSITES1
            ENDIF
            EI(:) = E(J7,:)

            DO J2 = J1 + 1, REALNATOMS

               IF (J2 <= NPS) THEN
                  NRBS1 = 1
                  NRBS2 = NRBSITES1
               ELSE
                  NRBS1 = NRBSITES1 + 1
                  NRBS2 = NRBSITES
               ENDIF

               J4 = 3*J2
               J6 = OFFSET + J4

!     LJ CONTRIBUTION

               DO J = NRBS1, NRBS2

                  IF (J2 <= NPS) THEN
                     J8 = NRBSITES1*(J2-1) + J
                  ELSE
                     J8 = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J2-NPS-1) + J - NRBSITES1
                  ENDIF

                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2
                  R6     = R2*R2*R2
                  R12    = R6*R6
            
                  ENERGY = ENERGY + 4.D0*(R12 - R6)
            
                  IF (GTEST) THEN

                     DVDR = 4.D0*(-12.D0*R12 + 6.D0*R6)*R2

                     G(J3-2:J3)  = G(J3-2:J3) + DVDR*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR*RSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))

                  ENDIF

!     DIPOLAR CONTRIBUTION

                  R4     = R2*R2
                  EJ(:)  = E(J8,:)
                  ALP    = DOT_PRODUCT(NR(:),EI(:))
                  BET    = DOT_PRODUCT(NR(:),EJ(:))
                  GAM    = DOT_PRODUCT(EI(:),EJ(:))

                  ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

                  IF (GTEST) THEN

                     VR     = -DPFCT*R4*(GAM - 3.D0*ALP*BET)
                     VA     = -DPFCT*BET*R2/ABSRIJ
                     VB     = -DPFCT*ALP*R2/ABSRIJ
                     VG     =  DPFCT*R2/(3.D0*ABSRIJ)

                     FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
                     FIJEI  = VA/ABSRIJ
                     FIJEJ  = VB/ABSRIJ
                     FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

                     G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
                     G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

                     G(J5-2) = G(J5-2) + FIJN*DOT_PRODUCT(NR(:),DR1(J7,:))      &
                             + VA*(DOT_PRODUCT(NR(:),DE1(J7,:)) + DOT_PRODUCT(EI(:),DR1(J7,:))/ABSRIJ) &
                             + VB*(DOT_PRODUCT(EJ(:),DR1(J7,:))/ABSRIJ) + VG*DOT_PRODUCT(DE1(J7,:),EJ(:))
                     G(J5-1) = G(J5-1) + FIJN*DOT_PRODUCT(NR(:),DR2(J7,:))      &
                             + VA*(DOT_PRODUCT(NR(:),DE2(J7,:)) + DOT_PRODUCT(EI(:),DR2(J7,:))/ABSRIJ) &
                             + VB*(DOT_PRODUCT(EJ(:),DR2(J7,:))/ABSRIJ) + VG*DOT_PRODUCT(DE2(J7,:),EJ(:))
                     G(J5)   = G(J5) + FIJN*DOT_PRODUCT(NR(:),DR3(J7,:))        &
                             + VA*(DOT_PRODUCT(NR(:),DE3(J7,:)) + DOT_PRODUCT(EI(:),DR3(J7,:))/ABSRIJ) &
                             + VB*(DOT_PRODUCT(EJ(:),DR3(J7,:))/ABSRIJ) + VG*DOT_PRODUCT(DE3(J7,:),EJ(:))

                     G(J6-2) = G(J6-2) - FIJN*DOT_PRODUCT(NR(:),DR1(J8,:))      &
                             - VA*DOT_PRODUCT(EI(:),DR1(J8,:))/ABSRIJ + VB*(DOT_PRODUCT(NR(:),DE1(J8,:)) &
                             - DOT_PRODUCT(EJ(:),DR1(J8,:))/ABSRIJ) + VG*DOT_PRODUCT(EI(:),DE1(J8,:))
                     G(J6-1) = G(J6-1) - FIJN*DOT_PRODUCT(NR(:),DR2(J8,:))      &
                             - VA*DOT_PRODUCT(EI(:),DR2(J8,:))/ABSRIJ + VB*(DOT_PRODUCT(NR(:),DE2(J8,:)) &
                             - DOT_PRODUCT(EJ(:),DR2(J8,:))/ABSRIJ) + VG*DOT_PRODUCT(EI(:),DE2(J8,:))
                     G(J6)   = G(J6) - FIJN*DOT_PRODUCT(NR(:),DR3(J8,:))        &
                             - VA*DOT_PRODUCT(EI(:),DR3(J8,:))/ABSRIJ + VB*(DOT_PRODUCT(NR(:),DE3(J8,:)) &
                             - DOT_PRODUCT(EJ(:),DR3(J8,:))/ABSRIJ) + VG*DOT_PRODUCT(EI(:),DE3(J8,:))

                  ENDIF 

               ENDDO

            ENDDO
 
            IF (EFIELDT) THEN

            
               ENERGY = ENERGY - STOCKMU*EFIELD*EI(3)

               IF (GTEST) THEN

                  G(J5-2) = G(J5-2) - STOCKMU*EFIELD*DE1(J7,3)
                  G(J5-1) = G(J5-1) - STOCKMU*EFIELD*DE2(J7,3)
                  G(J5)   = G(J5)   - STOCKMU*EFIELD*DE3(J7,3)

               ENDIF

            ENDIF
       
         ENDDO

      ENDDO

      END SUBROUTINE MSTBIN
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMSTBIN()

      USE COMMONS

      IMPLICIT NONE
      INTEGER          :: J1, J2
      DOUBLE PRECISION :: SIN60, COS60, SIN36, COS36, PI, R

      PI    = 4.D0*DATAN(1.D0)
      SIN36 = DSIN(PI/5.D0)       ! DSQRT(10.D0 - 2.D0*DSQRT(5.D0))/4.D0
      COS36 = DCOS(PI/5.D0)       ! DSQRT((5.D0) + 1.D0)/4.D0
      R     = 0.5D0/SIN36

      DO J1 = 1, NRBSITES1

         SITE(J1,1) = R*DSIN(0.4D0*PI*FLOAT(J1-1))
         SITE(J1,2) = -R*DCOS(0.4D0*PI*FLOAT(J1-1))

      ENDDO

      SITE(:,3) = 0.D0

      RBUV(1,:) = SITE(3,:) - SITE(4,:)
      RBUV(2,:) = SITE(4,:) - SITE(5,:)
      RBUV(3,:) = SITE(5,:) - SITE(1,:)
      RBUV(4,:) = SITE(1,:) - SITE(2,:)
      RBUV(5,:) = SITE(2,:) - SITE(3,:)

      SIN60 = DSIN(PI/3.D0)       
      COS60 = DCOS(PI/3.D0)       
      R     = 1.D0

      DO J1 = NRBSITES1 + 1, NRBSITES

         SITE(J1,1) = R*DSIN(PI*FLOAT(J1-NRBSITES1-1)/3.D0)
         SITE(J1,2) = -R*DCOS(PI*FLOAT(J1-NRBSITES1-1)/3.D0)

      ENDDO

      RBUV(6,:)  = (/-1.D0,   0.D0,  0.D0/)  
      RBUV(7,:)  = (/-COS60, -SIN60, 0.D0/)
      RBUV(8,:)  = (/ COS60, -SIN60, 0.D0/)
      RBUV(9,:)  = (/ 1.D0,   0.D0,  0.D0/)
      RBUV(10,:) = (/ COS60,  SIN60, 0.D0/) 
      RBUV(11,:) = (/-COS60,  SIN60, 0.D0/)

      DO J1 = 1, NRBSITES

         RBUV(J1,:) = RBUV(J1,:)/DSQRT(DOT_PRODUCT(RBUV(J1,:),RBUV(J1,:)))

      ENDDO

      END SUBROUTINE DEFMSTBIN
