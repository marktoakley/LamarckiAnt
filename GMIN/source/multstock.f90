      SUBROUTINE MULTSTOCK (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, RBUV, DPMU, EFIELDT, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DVDR
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3)
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3), E(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: DE1(NATOMS*NRBSITES/2,3), DE2(NATOMS*NRBSITES/2,3), DE3(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3), DPFCT
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DADPJ1, DADPJ2, DADPJ3
      DOUBLE PRECISION :: DBDPI1, DBDPI2, DBDPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST

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

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
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

      ENDDO

      DO J1 = 1, REALNATOMS  

         J3 = 3*J1
         J5 = OFFSET + J3

         RI(:)  = X(J3-2:J3)

         DO I = 1, NRBSITES

            J7    = NRBSITES*(J1-1) + I
            EI(:) = E(J7,:)

            DO J2 = J1 + 1, REALNATOMS

               J4 = 3*J2
               J6 = OFFSET + J4

!     LJ CONTRIBUTION

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2
                  R6     = R2*R2*R2
                  R12    = R6*R6
            
!                  IF (I <= 3 .OR. J <= 3) THEN

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

!                  ENDIF

                  ENDIF

!     DIPOLAR CONTRIBUTION

                  R4     = R2*R2
                  EJ(:)  = E(J8,:)
                  ALP    = DOT_PRODUCT(NR(:),EI(:))
                  BET    = DOT_PRODUCT(NR(:),EJ(:))
                  GAM    = DOT_PRODUCT(EI(:),EJ(:))

                  DPFCT  = 3.D0*DPMU(I)*DPMU(J)
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

                     DOTI1 = DOT_PRODUCT(RSS,DR1(J7,:))
                     DOTI2 = DOT_PRODUCT(RSS,DR2(J7,:))
                     DOTI3 = DOT_PRODUCT(RSS,DR3(J7,:))

                     DOTJ1 =-DOT_PRODUCT(RSS,DR1(J8,:))
                     DOTJ2 =-DOT_PRODUCT(RSS,DR2(J8,:))
                     DOTJ3 =-DOT_PRODUCT(RSS,DR3(J8,:))

                     DADPI1 = DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI1 + DOT_PRODUCT(NR(:),DE1(J7,:))
                     DADPI2 = DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI2 + DOT_PRODUCT(NR(:),DE2(J7,:))
                     DADPI3 = DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI3 + DOT_PRODUCT(NR(:),DE3(J7,:))
                                   
                     DADPJ1 =-DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ1
                     DADPJ2 =-DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ2
                     DADPJ3 =-DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ3

                     DBDPI1 = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI1
                     DBDPI2 = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI2
                     DBDPI3 = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI3

                     DBDPJ1 =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ1 + DOT_PRODUCT(NR(:),DE1(J8,:))
                     DBDPJ2 =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ2 + DOT_PRODUCT(NR(:),DE2(J8,:))
                     DBDPJ3 =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ3 + DOT_PRODUCT(NR(:),DE3(J8,:))
                                   
                     DGDPI1 = DOT_PRODUCT(DE1(J7,:),EJ(:))
                     DGDPI2 = DOT_PRODUCT(DE2(J7,:),EJ(:))
                     DGDPI3 = DOT_PRODUCT(DE3(J7,:),EJ(:))

                     DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J8,:))
                     DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J8,:))
                     DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J8,:))

                     G(J5-2) = G(J5-2) + VR*DOTI1/ABSRIJ + VA*DADPI1 + VB*DBDPI1 + VG*DGDPI1
                     G(J5-1) = G(J5-1) + VR*DOTI2/ABSRIJ + VA*DADPI2 + VB*DBDPI2 + VG*DGDPI2
                     G(J5)   = G(J5)   + VR*DOTI3/ABSRIJ + VA*DADPI3 + VB*DBDPI3 + VG*DGDPI3

                     G(J6-2) = G(J6-2) + VR*DOTJ1/ABSRIJ + VA*DADPJ1 + VB*DBDPJ1 + VG*DGDPJ1
                     G(J6-1) = G(J6-1) + VR*DOTJ2/ABSRIJ + VA*DADPJ2 + VB*DBDPJ2 + VG*DGDPJ2
                     G(J6)   = G(J6)   + VR*DOTJ3/ABSRIJ + VA*DADPJ3 + VB*DBDPJ3 + VG*DGDPJ3

                  ENDIF 

               ENDDO

            ENDDO
 
            IF (EFIELDT) THEN

               ENERGY = ENERGY - DPMU(I)*EFIELD*EI(3)

               IF (GTEST) THEN

                  G(J5-2) = G(J5-2) - DPMU(I)*EFIELD*DE1(J7,3)
                  G(J5-1) = G(J5-1) - DPMU(I)*EFIELD*DE2(J7,3)
                  G(J5)   = G(J5)   - DPMU(I)*EFIELD*DE3(J7,3)

               ENDIF

            ENDIF
       
         ENDDO

      ENDDO

      END SUBROUTINE MULTSTOCK
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULT2STOCK()

      USE COMMONS

      IMPLICIT NONE

      SITE(1,:) = (/0.D0, 0.D0, 0.5D0/)
      SITE(2,:) = (/0.D0, 0.D0,-0.5D0/)
      RBUV(1,:) = (/0.D0, 1.D0, 0.D0/)
      RBUV(2,:) = (/0.D0, 1.D0, 0.D0/)

      END SUBROUTINE DEFMULT2STOCK

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULT3STOCK()

      USE COMMONS

      IMPLICIT NONE

      SITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
      SITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
      SITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)

      RBUV(1,:) = (/-1.D0, 0.D0, 0.D0/)
      RBUV(2,:) = (/0.5D0, 0.5D0*DSQRT(3.D0), 0.D0/)
      RBUV(3,:) = (/0.5D0, -0.5D0*DSQRT(3.D0), 0.D0/)
      
      END SUBROUTINE DEFMULT3STOCK

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULT4STOCK()

      USE COMMONS

      IMPLICIT NONE

!      DOUBLE PRECISION :: INVRT2

!      INVRT2    = 1.D0/DSQRT(2.D0)

!      SITE(1,:) = (/ 0.5D0,  0.5D0, 0.D0/)
!      SITE(2,:) = (/-0.5D0,  0.5D0, 0.D0/)
!      SITE(3,:) = (/-0.5D0, -0.5D0, 0.D0/)
!      SITE(4,:) = (/ 0.5D0, -0.5D0, 0.D0/)
!      RBUV(1,:) = (/ 1.D0,  -1.D0,  0.D0/)*INVRT2
!      RBUV(2,:) = (/ 1.D0,   1.D0,  0.D0/)*INVRT2
!      RBUV(3,:) = (/-1.D0,   1.D0,  0.D0/)*INVRT2
!      RBUV(4,:) = (/-1.D0,  -1.D0,  0.D0/)*INVRT2

      SITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
      SITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
      SITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
      SITE(4,:) = (/ 0.D0, 0.D0, 0.D0/)

      RBUV(1,:) = (/-1.D0, 0.D0, 0.D0/)
      RBUV(2,:) = (/0.5D0, 0.5D0*DSQRT(3.D0), 0.D0/)
      RBUV(3,:) = (/0.5D0, -0.5D0*DSQRT(3.D0), 0.D0/)
      RBUV(4,:) = (/0.D0, 0.D0, 1.D0/)

      END SUBROUTINE DEFMULT4STOCK

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULT5STOCK()

      USE COMMONS

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: PI, SIN36, RD

      PI    = 4.D0*DATAN(1.D0)
      SIN36 = DSIN(PI/5.D0)       ! DSQRT(10.D0 - 2.D0*DSQRT(5.D0))/4.D0
      RD    = 0.5D0/SIN36

      DO J1 = 1, 5

         SITE(J1,1) =  RD*DSIN(0.4D0*PI*FLOAT(J1-1))
         SITE(J1,2) = -RD*DCOS(0.4D0*PI*FLOAT(J1-1))

      ENDDO

      SITE(:,3) = 0.D0

      RBUV(1,:) = SITE(3,:) - SITE(4,:)
      RBUV(2,:) = SITE(4,:) - SITE(5,:)
      RBUV(3,:) = SITE(5,:) - SITE(1,:)
      RBUV(4,:) = SITE(1,:) - SITE(2,:)
      RBUV(5,:) = SITE(2,:) - SITE(3,:)

      DO J1 = 1, 5

         RBUV(J1,:) = RBUV(J1,:)/DSQRT(DOT_PRODUCT(RBUV(J1,:),RBUV(J1,:)))

      ENDDO

      END SUBROUTINE DEFMULT5STOCK
