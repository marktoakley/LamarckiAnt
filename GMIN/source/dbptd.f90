      SUBROUTINE DMBLTD (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, DBPMU, EFIELDT, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), SITETD(4,3), SIGTD(4), SIGDB(2)
      DOUBLE PRECISION :: ENERGY, DVDR, R2, R6, R12, ABSRIJ, RIJSQ, DPFCT, SIGMA
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RSS(3), NR(3), P(3), DU(3), EI(3), EJ(3), R(NATOMS+4,3) 
      DOUBLE PRECISION :: DR1(NATOMS+4,3), DR2(NATOMS+4,3), DR3(NATOMS+4,3) 
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS, CLJ12, CLJ6
      LOGICAL          :: GTEST

      CALL DEFDUM(DU, DPFCT, CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS)
      CALL DEFTD(SITETD,SIGTD,SIGDB)

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS-1

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

            J4        = (NRBSITES-1)*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))

            IF (GTEST) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS - 1 

         J3 = 3*J1
         J5 = OFFSET + J3
 
         RI(:)  = X(J3-2:J3)
         EI(:)  = E(J1,:)

         DO J2 = J1 + 1, REALNATOMS - 1

            J4 = 3*J2
            J6 = OFFSET + J4

!     LJ CONTRIBUTION

            DO I = 1, NRBSITES - 1

               J7 = (NRBSITES-1)*(J1-1) + I

                  DO J = 1, NRBSITES - 1

                  J8     = (NRBSITES-1)*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)

                  R2     = DOT_PRODUCT(RSS(:),RSS(:))

                  R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
                  R6     = R2*R2*R2
                  R12    = R6*R6
 
                  IF (I == 1 .AND. J == 1) THEN
                     ENERGY = ENERGY + (CLJ12BB*R6 - CLJ6BB)*R6
                  ELSEIF (I == 2 .AND. J == 2) THEN
                     ENERGY = ENERGY + (CLJ12SS*R6 - CLJ6SS)*R6
                  ELSE
                     ENERGY = ENERGY + (CLJ12BS*R6 - CLJ6BS)*R6
                  ENDIF
 
                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     IF (I == 1 .AND. J == 1) THEN
                        DVDR   = -6.D0*(2.D0*CLJ12BB*R6 - CLJ6BB)*R6*R2
                     ELSEIF (I == 2 .AND. J == 2) THEN
                        DVDR   = -6.D0*(2.D0*CLJ12SS*R6 - CLJ6SS)*R6*R2
                     ELSE
                        DVDR   = -6.D0*(2.D0*CLJ12BS*R6 - CLJ6BS)*R6*R2
                     ENDIF

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

         IF (EFIELDT) THEN

            ENERGY = ENERGY - DBPMU*EFIELD*EI(3)

            IF (GTEST) THEN

               G(J5-2) = G(J5-2) - DBPMU*EFIELD*DE1(J1,3)
               G(J5-1) = G(J5-1) - DBPMU*EFIELD*DE2(J1,3)
               G(J5)   = G(J5)   - DBPMU*EFIELD*DE3(J1,3)

            ENDIF  

         ENDIF

      ENDDO

!     LJ CONTRIBUTION WITH TD

      J4 = 3*REALNATOMS
      J6 = 3*NATOMS
      RI = X(J4-2:J4)
      P  = X(J6-2:J6)

      CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)
         
      DO J = 1, 4 

         J8        = (NRBSITES-1)*(REALNATOMS-1) + J
         R(J8,:)   = RI(:) + MATMUL(RMI,SITETD(J,:))

         IF (GTEST) THEN

            DR1(J8,:) = MATMUL(DRMI1,SITETD(J,:))
            DR2(J8,:) = MATMUL(DRMI2,SITETD(J,:))
            DR3(J8,:) = MATMUL(DRMI3,SITETD(J,:))

         ENDIF

      ENDDO 

      DO J1 = 1, REALNATOMS - 1 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO I = 1, NRBSITES - 1

            J7 = (NRBSITES-1)*(J1-1) + I

            DO J = 1, 4

               J8     = (NRBSITES-1)*(REALNATOMS-1) + J 
         
               RSS(:) = R(J7,:) - R(J8,:)

               SIGMA  = 0.5D0*(SIGDB(I)+SIGTD(J))
               CLJ12  = 4.D0*SIGMA**12
               CLJ6   = 4.D0*SIGMA**6
               R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
               R6     = R2*R2*R2
               R12    = R6*R6

               ENERGY = ENERGY + (CLJ12*R6 - CLJ6)*R6

               IF (GTEST) THEN
!     DVDR = DVDR/R
                  DVDR   = -6.D0*(2.D0*CLJ12*R6 - CLJ6)*R6*R2

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
 
      ENDDO

      END SUBROUTINE DMBLTD 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTD(SITETD,SIGTD,SIGDB)

      USE COMMONS

      IMPLICIT NONE

      DOUBLE PRECISION :: SITETD(4,3), SIGTD(4), SIGDB(2), FCTR


      FCTR        = 1.D0/DSQRT(8.D0)
      SITETD(1,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      SITETD(2,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      SITETD(3,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)
      SITETD(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)
      SIGTD(:)    = (/1.D0, 0.9D0, 0.8D0, 0.7D0/)
      SIGDB(:)    = (/1.D0, 0.4D0/)
   
      END SUBROUTINE DEFTD

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWDMBLTD()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3), SITETD(4,3), FCTR
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='dmbltd.xyz', STATUS='UNKNOWN')

      FCTR        = 1.D0/DSQRT(8.D0)
      SITETD(1,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      SITETD(2,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      SITETD(3,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)
      SITETD(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)

      GTEST = .FALSE.

      DO J1 = 1, NSAVE

         WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, QMINNATOMS(J1)/2

            J5   = 3*J3
            J7   = 3*QMINNATOMS(J1)/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            IF (J3 < QMINNATOMS(J1)/2) THEN 

                DO J2 = 1, NRBSITES - 1

                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
                  IF (J2 == 1) THEN
                     WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF

               ENDDO

            ELSE

                DO J2 = 1, 4

                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITETD(J2,:))
                  IF (J2 == 1) THEN
                     WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSEIF (J2 == 2) THEN
                     WRITE(26,'(A4,3F20.10)') 'N', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSEIF (J2 == 3) THEN
                     WRITE(26,'(A4,3F20.10)') 'B', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'S', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF

               ENDDO

            ENDIF  

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWDMBLTD
