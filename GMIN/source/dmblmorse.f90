      SUBROUTINE DMBLMORSE(X, G, ENERGY, GTEST)

!     each rigidbody is a dumbbell consisting of 2 Morse sites

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, EPS11, EPS22, EPS12, MRHO11, MRHO22, MRHO12, REQ11, REQ22, REQ12,      &
                         DBPMU, EFIELDT, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, R2, ABSRIJ, EPS, RHO, REQ, FCTR1, FCTR2, RIJSQ
      DOUBLE PRECISION :: RSS(3), NR(3), P(3), EI(3), EJ(3), DU(3), RI(3), RJ(3), RIJ(3)
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3), E(NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3), DPFCT
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DADPJ1, DADPJ2, DADPJ3
      DOUBLE PRECISION :: DBDPI1, DBDPI2, DBDPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST

      DPFCT = 3.D0*DBPMU*DBPMU
      DU(:) = (/0.D0, 1.D0, 0.D0/)
      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = X(J3-2:J3) + MATMUL(RMI(:,:),SITE(J2,:))

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

            ENDIF

         ENDDO

         E(J1,:)   = MATMUL(RMI(:,:),DU(:)) 

         IF (GTEST) THEN

            DE1(J1,:) = MATMUL(DRMI1(:,:),DU(:))
            DE2(J1,:) = MATMUL(DRMI2(:,:),DU(:))
            DE3(J1,:) = MATMUL(DRMI3(:,:),DU(:))

         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS  

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

            DO I = 1, NRBSITES

               J7    = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2

                  IF (I == 1 .AND. J == 1) THEN
                     EPS = EPS11; RHO = MRHO11; REQ = REQ11
                  ELSEIF (I == 2 .AND. J == 2) THEN
                     EPS = EPS22; RHO = MRHO22; REQ = REQ22
                  ELSE
                     EPS = EPS12; RHO = MRHO12; REQ = REQ12
                  ENDIF

                  FCTR1  = EXP(RHO*(1.D0 - ABSRIJ/REQ))
                  FCTR2  = (1.D0 - FCTR1) 
                  ENERGY = ENERGY + EPS*(FCTR2*FCTR2 - 1.D0)
!     DVDR = DVDR/R
                  IF (GTEST) THEN

                     DVDR       = 2.D0*EPS*RHO*FCTR1*FCTR2/(REQ*ABSRIJ)
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

            IF (DBPMU == 0.D0) GO TO 10 
!     DIPOLAR CONTRIBUTION
            RI(:)  = X(J3-2:J3)
            RJ(:)  = X(J4-2:J4)
            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ
            EI(:)  = E(J1,:)
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

10         ENDDO
 
         IF (EFIELDT) THEN

            ENERGY = ENERGY - DBPMU*EFIELD*EI(3)

            IF (GTEST) THEN

               G(J5-2) = G(J5-2) - DBPMU*EFIELD*DE1(J1,3)
               G(J5-1) = G(J5-1) - DBPMU*EFIELD*DE2(J1,3)
               G(J5)   = G(J5)   - DBPMU*EFIELD*DE3(J1,3)

            ENDIF

         ENDIF

      ENDDO

      END SUBROUTINE DMBLMORSE
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDMBL()

      USE COMMONS, ONLY: SITE, EPS11, EPS22, EPS12, MRHO11, MRHO22, MRHO12, REQ11, REQ22, REQ12

      IMPLICIT NONE
      
      SITE(1,:) = (/ 0.D0, 0.D0, 0.5D0*REQ11/)
      SITE(2,:) = (/ 0.D0, 0.D0,-0.5D0*REQ22/)

      EPS12  = DSQRT(EPS11*EPS22)
      MRHO12 = 0.5D0*(MRHO11+MRHO22)
      REQ12  = 0.5D0*(REQ11+REQ22)

      END SUBROUTINE DEFDMBL

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWDMBL()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='dmbl.xyz', STATUS='UNKNOWN')

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

            DO J2 = 1, NRBSITES 

               RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (J2 == 1) THEN
                  WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWDMBL
