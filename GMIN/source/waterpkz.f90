      SUBROUTINE WATERPKZ (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, K1, K2, ITR, REALNATOMS, OFFSET
      INTEGER, PARAMETER :: ITRMAX = 1000 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: CHARGE(NRBSITES), C12, C6, CC, POLAR
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DR(3), RSS(3), PI(3), PJ(3)
      DOUBLE PRECISION :: DSS2, R2, R6, R12, ABSR, R3, R5, DVDR, ENERGY, EP0, EP
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: T(3*NATOMS/2, 3*NATOMS/2)
      DOUBLE PRECISION :: EFIELD(3*NATOMS/2), DFIELD(3*NATOMS/2,3*NATOMS), EFIELD0(3*NATOMS/2), DM(3*NATOMS/2)
      DOUBLE PRECISION :: EX(3), EY(3), EZ(3), DUMMY(3), SDUMMY
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      DOUBLE PRECISION :: DT1(3,3), DT2(3,3), DT3(3,3), DT4(3,3), DT5(3,3), DT6(3,3)
      LOGICAL          :: GTEST, SCF
     
      CALL DEFWATERPKZ(CHARGE, C6, C12, CC, POLAR)
 
      SCF = .TRUE.

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      ENERGY     = 0.D0
      EFIELD0(:) = 0.D0 
      T(:,:)     = 0.D0

      IF (GTEST) THEN
         G(:) = 0.D0
         IF (SCF) THEN
            EX(:) = (/1.D0, 0.D0, 0.D0/)
            EY(:) = (/0.D0, 1.D0, 0.D0/)
            EZ(:) = (/0.D0, 0.D0, 1.D0/)
            DFIELD(:,:) = 0.D0
         ENDIF
      ENDIF

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         PI = X(J5-2:J5)

         CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
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

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     M-M LJ CONTRIBUTION

            J7     = NRBSITES*(J1-1) + 4 
            J8     = NRBSITES*(J2-1) + 4
            RSS(:) = R(J7,:) - R(J8,:)
            R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
            R6     = R2 * R2 * R2
            R12    = R6 ** 2
            ENERGY = ENERGY + C12*R12 - C6*R6

            IF (GTEST) THEN
!     DVDR = DVDR/R
               DVDR       = -6.D0*(2.D0*C12*R12 - C6*R6)*R2

               G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:)

               G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
               G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
               G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))

               G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
               G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
               G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))

            ENDIF

!     SUM OVER CHARGED SITES

            DO I = 1, NRBSITES

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(RSS(:), RSS(:))
                  R2     = 1.D0/DSS2
                  ABSR   = SQRT(DSS2)
                  ENERGY = ENERGY + CC*CHARGE(I)*CHARGE(J)/ABSR

                  IF (SCF) THEN

                     IF (I == 4) THEN
                        EFIELD0(J3-2:J3) = EFIELD0(J3-2:J3) + CHARGE(J)*RSS(:)*R2/ABSR
                     ENDIF
                     IF (J == 4) THEN
                        EFIELD0(J4-2:J4) = EFIELD0(J4-2:J2) - CHARGE(I)*RSS(:)*R2/ABSR
                     ENDIF

                     IF (I == 4 .AND. J == 4) THEN

                        R3 = R2/ABSR
                        R5 = R2*R3

                        T(J3-2,J4-2) = 3.D0*RSS(1)*RSS(1)*R5 - R3
                        T(J3-1,J4-1) = 3.D0*RSS(2)*RSS(2)*R5 - R3
                        T(J3,J4)     = 3.D0*RSS(3)*RSS(3)*R5 - R3  
                        T(J3-2,J4-1) = 3.D0*RSS(1)*RSS(2)*R5 
                        T(J3-1,J4)   = 3.D0*RSS(2)*RSS(3)*R5
                        T(J3,J4-2)   = 3.D0*RSS(3)*RSS(1)*R5
                        T(J3-1,J4-2) = T(J3-2,J4-1)
                        T(J3,J4-1)   = T(J3-1,J4)
                        T(J3-2,J4)   = T(J3,J4-2)
                        
                        T(J4-2,J3-2) = T(J3-2,J4-2)
                        T(J4-1,J3-1) = T(J3-1,J4-1) 
                        T(J4,J3)     = T(J3,J4)
                        T(J4-2,J3-1) = T(J3-2,J4-1)
                        T(J4-1,J3)   = T(J3-1,J4)
                        T(J4,J3-2)   = T(J3,J4-2) 
                        T(J4-1,J3-2) = T(J3-1,J4-2)
                        T(J4,J3-1)   = T(J3,J4-1)
                        T(J4-2,J3)   = T(J3-2,J4)

                     ENDIF

                  ENDIF

                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     R3   = R2/ABSR
                     DVDR = -CC*CHARGE(I)*CHARGE(J)*R3 

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

         ENDDO

      ENDDO

      IF (SCF) THEN

         EFIELD(:) = EFIELD0(:)
         EP0       = DOT_PRODUCT(EFIELD(:),EFIELD(:))
         EP0       =-0.5D0*CC*POLAR*EP0

         ITR = 0

         DO WHILE (ITR < ITRMAX)

            ITR       = ITR + 1
            DM(:)     = POLAR*EFIELD(:)
            EP        = 0.D0
            EFIELD(:) = 0.D0
            DO J1 = 1, REALNATOMS
               J3 = 3*J1
               DO J2 = 1, REALNATOMS
                  IF (J1 == J2) CYCLE
                  J4 = 3*J2
                  EFIELD(J3-2:J3) = EFIELD(J3-2:J3) + MATMUL(T(J3-2:J3,J4-2:J4),DM(J4-2:J4))
               ENDDO
               EFIELD(J3-2:J3) = EFIELD(J3-2:J3) + EFIELD0(J3-2:J3)
            ENDDO
            EP = DOT_PRODUCT(EFIELD0(:),EFIELD(:))
            EP = -0.5D0*CC*POLAR*EP
            IF (ABS(EP - EP0) < 1.D-08) EXIT
            IF (ITR == ITRMAX) THEN
               WRITE(*,*) 'POLARISATION ENERGY IS NOT CONVERGED'
               STOP
            ENDIF
            EP0 = EP

         ENDDO

         DM(:) = POLAR*EFIELD(:)

      ENDIF

      ENERGY = ENERGY + EP

      IF (SCF .AND. GTEST) THEN

         DO J1 = 1, REALNATOMS   

            J3 = 3*J1
            J5 = OFFSET + J3
            I  = 4
            J7 = NRBSITES*(J1-1) + I

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2 
               J6 = OFFSET + J4

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(RSS(:), RSS(:))
                  R2     = 1.D0/DSS2
                  R3     = R2/SQRT(DSS2)
                  R5     = R3*R2

!     DERIVATIVES OF THE ELECTRIC FIELD

                  DUMMY(:)             = CHARGE(J)*R3*(-3.D0*R2*RSS(1)*RSS(:) + EX(:))
                  DFIELD(J3-2:J3,J3-2) = DFIELD(J3-2:J3,J3-2) + DUMMY(:) 
                  DFIELD(J3-2:J3,J4-2) = DFIELD(J3-2:J3,J4-2) - DUMMY(:)

                  DUMMY(:)             = CHARGE(J)*R3*(-3.D0*R2*RSS(2)*RSS(:) + EY(:))
                  DFIELD(J3-2:J3,J3-1) = DFIELD(J3-2:J3,J3-1) + DUMMY(:)
                  DFIELD(J3-2:J3,J4-1) = DFIELD(J3-2:J3,J4-1) - DUMMY(:)

                  DUMMY(:)             = CHARGE(J)*R3*(-3.D0*R2*RSS(3)*RSS(:) + EZ(:))
                  DFIELD(J3-2:J3,J3)   = DFIELD(J3-2:J3,J3) + DUMMY(:)
                  DFIELD(J3-2:J3,J4)   = DFIELD(J3-2:J3,J4) - DUMMY(:)

                  DOTI1   = DOT_PRODUCT(RSS(:),DR1(J7,:))
                  DOTI2   = DOT_PRODUCT(RSS(:),DR2(J7,:))
                  DOTI3   = DOT_PRODUCT(RSS(:),DR3(J7,:))

                  DOTJ1   = DOT_PRODUCT(RSS(:),DR1(J8,:))
                  DOTJ2   = DOT_PRODUCT(RSS(:),DR2(J8,:))
                  DOTJ3   = DOT_PRODUCT(RSS(:),DR3(J8,:))

                  DFIELD(J3-2:J3,J5-2) = DFIELD(J3-2:J3,J5-2) + CHARGE(J)*R3*(-3.D0*R2*DOTI1*RSS(:) + DR1(J7,:))
                  DFIELD(J3-2:J3,J6-2) = DFIELD(J3-2:J3,J6-2) - CHARGE(J)*R3*(-3.D0*R2*DOTJ1*RSS(:) + DR1(J8,:))  

                  DFIELD(J3-2:J3,J5-1) = DFIELD(J3-2:J3,J5-1) + CHARGE(J)*R3*(-3.D0*R2*DOTI2*RSS(:) + DR2(J7,:))
                  DFIELD(J3-2:J3,J6-1) = DFIELD(J3-2:J3,J6-1) - CHARGE(J)*R3*(-3.D0*R2*DOTJ2*RSS(:) + DR2(J8,:))
 
                  DFIELD(J3-2:J3,J5)   = DFIELD(J3-2:J3,J5)   + CHARGE(J)*R3*(-3.D0*R2*DOTI3*RSS(:) + DR3(J7,:))
                  DFIELD(J3-2:J3,J6)   = DFIELD(J3-2:J3,J6)   - CHARGE(J)*R3*(-3.D0*R2*DOTJ3*RSS(:) + DR3(J8,:))

               ENDDO

            ENDDO

         ENDDO

      ENDIF

      IF (GTEST .AND. SCF) THEN

         DO J1 = 1, 3*NATOMS

            G(J1) = G(J1) - CC*DOT_PRODUCT(DFIELD(:,J1),DM(:))

         ENDDO

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3
            I  = 4
            J7 = NRBSITES*(J1-1) + I

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4
               J  = 4
               J8 = NRBSITES*(J2-1) + J

               RSS(:) = R(J7,:) - R(J8,:)
               DSS2   = DOT_PRODUCT(RSS(:), RSS(:))
               R2     = 1.D0/DSS2
               R3     = R2/SQRT(DSS2)
               R5     = R3*R2

!     DERIVATIVES OF THE T MATRIX

               DT1(1,1) = -3.D0*R5*RSS(1)*(5.D0*R2*RSS(1)*RSS(1) - 3.D0)
               DT1(2,2) = -3.D0*R5*RSS(1)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0)
               DT1(3,3) = -3.D0*R5*RSS(1)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0)
               DT1(1,2) = -3.D0*R5*RSS(2)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0)
               DT1(2,3) = -15.D0*R5*R2*RSS(1)*RSS(2)*RSS(3) 
               DT1(3,1) = -3.D0*R5*RSS(3)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0)
               DT1(2,1) = DT1(1,2)
               DT1(1,3) = DT1(3,1)
               DT1(3,2) = DT1(2,3)

               DT2(1,1) = -3.D0*R5*RSS(2)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0)
               DT2(2,2) = -3.D0*R5*RSS(2)*(5.D0*R2*RSS(2)*RSS(2) - 3.D0)
               DT2(3,3) = -3.D0*R5*RSS(2)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0)
               DT2(1,2) = -3.D0*R5*RSS(1)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0)
               DT2(2,3) = -3.D0*R5*RSS(3)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0)
               DT2(3,1) = -15.D0*R5*R2*RSS(1)*RSS(2)*RSS(3) 
               DT2(2,1) = DT2(1,2)
               DT2(1,3) = DT2(3,1)
               DT2(3,2) = DT2(2,3)

               DT3(1,1) = -3.D0*R5*RSS(3)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0)
               DT3(2,2) = -3.D0*R5*RSS(3)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0)
               DT3(3,3) = -3.D0*R5*RSS(3)*(5.D0*R2*RSS(3)*RSS(3) - 3.D0)
               DT3(1,2) = -15.D0*R5*R2*RSS(1)*RSS(2)*RSS(3)
               DT3(2,3) = -3.D0*R5*RSS(2)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0)
               DT3(3,1) = -3.D0*R5*RSS(1)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0)
               DT3(2,1) = DT3(1,2)
               DT3(1,3) = DT3(3,1)
               DT3(3,2) = DT3(2,3)

               DOTI1   = DOT_PRODUCT(RSS(:),DR1(J7,:))
               DOTI2   = DOT_PRODUCT(RSS(:),DR2(J7,:))
               DOTI3   = DOT_PRODUCT(RSS(:),DR3(J7,:))

               DOTJ1   = DOT_PRODUCT(RSS(:),DR1(J8,:))
               DOTJ2   = DOT_PRODUCT(RSS(:),DR2(J8,:))
               DOTJ3   = DOT_PRODUCT(RSS(:),DR3(J8,:))

               DT4(1,1) = -3.D0*R5*DOTI1*(5.D0*R2*RSS(1)*RSS(1) - 1.D0) + 6.D0*R5*RSS(1)*DR1(J7,1)
               DT4(2,2) = -3.D0*R5*DOTI1*(5.D0*R2*RSS(2)*RSS(2) - 1.D0) + 6.D0*R5*RSS(2)*DR1(J7,2)
               DT4(3,3) = -3.D0*R5*DOTI1*(5.D0*R2*RSS(3)*RSS(3) - 1.D0) + 6.D0*R5*RSS(3)*DR1(J7,3)
               DT4(1,2) = -15.D0*R5*R2*RSS(1)*RSS(2)*DOTI1 + 3.D0*R5*(RSS(1)*DR1(J7,2) + RSS(2)*DR1(J7,1))
               DT4(2,3) = -15.D0*R5*R2*RSS(2)*RSS(3)*DOTI1 + 3.D0*R5*(RSS(2)*DR1(J7,3) + RSS(3)*DR1(J7,2))
               DT4(3,1) = -15.D0*R5*R2*RSS(3)*RSS(1)*DOTI1 + 3.D0*R5*(RSS(3)*DR1(J7,1) + RSS(1)*DR1(J7,3))
               DT4(2,1) = DT4(1,2)
               DT4(1,3) = DT4(3,1)
               DT4(3,2) = DT4(2,3)

               DT5(1,1) = -3.D0*R5*DOTI2*(1)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0) + 6.D0*R5*RSS(1)*DR2(J7,1)
               DT5(2,2) = -3.D0*R5*DOTI2*(1)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0) + 6.D0*R5*RSS(2)*DR2(J7,2)
               DT5(3,3) = -3.D0*R5*DOTI2*(1)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0) + 6.D0*R5*RSS(3)*DR2(J7,3)
               DT5(1,2) = -15.D0*R5*R2*RSS(1)*RSS(2)*DOTI2 + 3.D0*R5*(RSS(1)*DR2(J7,2) + RSS(2)*DR2(J7,1))
               DT5(2,3) = -15.D0*R5*R2*RSS(2)*RSS(3)*DOTI2 + 3.D0*R5*(RSS(2)*DR2(J7,3) + RSS(3)*DR2(J7,2))
               DT5(3,1) = -15.D0*R5*R2*RSS(3)*RSS(1)*DOTI2 + 3.D0*R5*(RSS(3)*DR2(J7,1) + RSS(1)*DR2(J7,3))
               DT5(2,1) = DT5(1,2)
               DT5(1,3) = DT5(3,1)
               DT5(3,2) = DT5(2,3)

               DT6(1,1) = -3.D0*R5*DOTI3*(1)*(5.D0*R2*RSS(1)*RSS(1) - 1.D0) + 6.D0*R5*RSS(1)*DR3(J7,1)
               DT6(2,2) = -3.D0*R5*DOTI3*(1)*(5.D0*R2*RSS(2)*RSS(2) - 1.D0) + 6.D0*R5*RSS(2)*DR3(J7,2)
               DT6(3,3) = -3.D0*R5*DOTI3*(1)*(5.D0*R2*RSS(3)*RSS(3) - 1.D0) + 6.D0*R5*RSS(3)*DR3(J7,3)
               DT6(1,2) = -15.D0*R5*R2*RSS(1)*RSS(2)*DOTI3 + 3.D0*R5*(RSS(1)*DR3(J7,2) + RSS(2)*DR3(J7,1))
               DT6(2,3) = -15.D0*R5*R2*RSS(2)*RSS(3)*DOTI3 + 3.D0*R5*(RSS(2)*DR3(J7,3) + RSS(3)*DR3(J7,2))
               DT6(3,1) = -15.D0*R5*R2*RSS(3)*RSS(1)*DOTI3 + 3.D0*R5*(RSS(3)*DR3(J7,1) + RSS(1)*DR3(J7,3))
               DT6(2,1) = DT6(1,2)
               DT6(1,3) = DT6(3,1)
               DT6(3,2) = DT6(2,3)

               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT1(K1,K2)
                  ENDDO
               ENDDO
               G(J3-2) = G(J3-2) - CC*SDUMMY
  
               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT2(K1,K2)
                  ENDDO
               ENDDO
               G(J3-1) = G(J3-1) - CC*SDUMMY

               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT3(K1,K2)
                  ENDDO
               ENDDO
               G(J3) = G(J3) - CC*SDUMMY

               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT4(K1,K2)
                  ENDDO
               ENDDO
               G(J5-2) = G(J5-2) - CC*SDUMMY

               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT5(K1,K2)
                  ENDDO
               ENDDO
               G(J5-1) = G(J5-1) - CC*SDUMMY

               SDUMMY = 0.D0
               DO K1 = 1, 3
                  DO K2 = 1, 3
                     SDUMMY = SDUMMY + DM(J3-3+K1)*DM(J4-3+K2)*DT6(K1,K2)
                  ENDDO
               ENDDO
               G(J5) = G(J5) - CC*SDUMMY

            ENDDO

         ENDDO

      ENDIF

      END SUBROUTINE  WATERPKZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFWATERPKZ(CHARGE, C6, C12, CC, POLAR)
    
      USE COMMONS, ONLY: NRBSITES, SITE  
      
      IMPLICIT NONE

      DOUBLE PRECISION :: CHARGE(NRBSITES), C6, C12, CC, POLAR 

      SITE(1,1) = 0.0D0
      SITE(1,2) = 0.0D0
      SITE(1,3) = 0.0D0
      SITE(2,1) = 0.756690D0
      SITE(2,2) = 0.0D0
      SITE(2,3) =-0.585892D0
      SITE(3,1) =-0.756690D0
      SITE(3,2) = 0.0D0
      SITE(3,3) =-0.585892D0
      SITE(4,1) = 0.0D0
      SITE(4,2) = 0.0D0
      SITE(4,3) =-0.138D0

      C6        = 1704.76651D0 ! LJ coefficients in kcal/mol Angstrom**6 or Angstrom**12
      C12       = 1729898.17D0
      CHARGE(1) = 1.2456D0
      CHARGE(2) = 0.6228D0
      CHARGE(3) = 0.6228D0
      CHARGE(4) =-2.4912D0
      CC        = 332.06378D0
      POLAR     = 1.47D0

      END SUBROUTINE DEFWATERPKZ
