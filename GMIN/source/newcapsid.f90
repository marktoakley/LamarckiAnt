      SUBROUTINE NEWCAPSID (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, RHO
 
      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DSS(3), PI(3)
      DOUBLE PRECISION :: DSS2, RSS, DVDR, ENERGY, SIGMASQ, VSS
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: CAPEPS2, CAPRAD, CAPHEIGHT
      LOGICAL          :: GTEST
      COMMON / CAPS /  CAPEPS2, CAPRAD, CAPHEIGHT

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      
      OFFSET  = 3*NATOMS/2
      CAPRAD  = 2.D0*CAPRAD*SIN(0.8D0*DATAN(1.D0))
      SIGMASQ = (1.D0 + CAPRAD*SQRT(0.5D0*(5.D0 + SQRT(5.D0))))**2
    
      DO J1 = 1, NATOMS/2
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

      DO J1 = 1, NATOMS/2 - 1 
         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, NATOMS/2
            J4 = 3*J2
            J6 = OFFSET + J4
!     Sum over Morse sites
            DO I = 1, NRBSITES - 2
               J7 = NRBSITES*(J1-1) + I
               DO J = 1, NRBSITES - 2
                  J8 = NRBSITES*(J2-1) + J
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(DSS(:),DSS(:))
                  RSS    = SQRT(DSS2)
                  VSS    = EXP(RHO*(1.D0 - RSS))
                  ENERGY = ENERGY + (1.D0 - VSS)*(1.D0 - VSS) - 1.D0 

                  IF (GTEST) THEN
                     DVDR       = 2.D0*RHO*(-VSS*VSS + VSS)/RSS
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))
                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))
                  ENDIF
               ENDDO
            ENDDO
!     Sum over apex sites
            DO I = 7, 8
               DO J = 7, 8
                  IF (I == 8 .AND. J == 8) CYCLE
                  J7     = NRBSITES*(J1-1) + I
                  J8     = NRBSITES*(J2-1) + J
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(DSS,DSS)
                  VSS    = CAPEPS2*(SIGMASQ/DSS2)**6
                  ENERGY = ENERGY + VSS
                  IF (GTEST) THEN
                     DVDR   = -12.D0*VSS/DSS2
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))
                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE NEWCAPSID 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCAPSID(CAPRAD,CAPHEIGHT)
    
      USE COMMONS, ONLY: SITE, NRBSITES  
      
      IMPLICIT NONE
  
      INTEGER :: J1 
      DOUBLE PRECISION :: CAPRAD, CAPHEIGHT
      DOUBLE PRECISION :: PI, SIN36

!     The parametric equation:
!     x = r * cos(theta), y = r * sin(theta), where 0 <= theta < 2 * PI
!     can generate the vertices of an n-side regular polygon by assigning n different values to
!     theta according to 2 * PI * i / n, where 0 <= i < n.

      PI    = 4.D0*DATAN(1.D0)
      SIN36 = DSIN(PI/5.D0)
      SITE(:,3) = 0.D0
      CAPRAD = 2.D0*CAPRAD*SIN36

      DO J1 = 1, NRBSITES - 2
!         SITE(J1,1) = CAPRAD*DSIN(0.4D0*PI*FLOAT(J1-1))
!         SITE(J1,2) =-CAPRAD*DCOS(0.4D0*PI*FLOAT(J1-1))
         SITE(J1,1) = CAPRAD*DSIN(PI*FLOAT(J1-1)/3.D0)
         SITE(J1,2) =-CAPRAD*DCOS(PI*FLOAT(J1-1)/3.D0)
      ENDDO
      SITE(NRBSITES-1,:) = (/0.D0, 0.D0, CAPHEIGHT/)
      SITE(NRBSITES,:)   = (/0.D0, 0.D0,-CAPHEIGHT/)

      END SUBROUTINE DEFCAPSID

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWNEWCAPSID()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='newcapsid.xyz', STATUS='UNKNOWN')

      GTEST = .FALSE.

      DO J1 = 1, NSAVE

         WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, NATOMS/2

            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            DO J2 = 1, NRBSITES

               RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (J2 <= 5) THEN
                  WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWNEWCAPSID
