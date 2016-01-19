      SUBROUTINE CAPBIN(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NPS, NRBSITES, NRBSITES1, NTSITES, SITE, CAPEPS2, CAPRAD, CAPRHO

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, NRBSI1, NRBSI2, NRBSJ1, NRBSJ2 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R12, ABSR, DVDR, VSS, CAPRADP, CAPRADH, SIGMASQ, SIGMASQP, SIGMASQH, SIGMASQM
      DOUBLE PRECISION :: RI(3), RSS(3), P(3)
      DOUBLE PRECISION :: R(NTSITES,3)
      DOUBLE PRECISION :: DR1(NTSITES,3), DR2(NTSITES,3), DR3(NTSITES,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      LOGICAL          :: GTEST

!     Input parameters
!     CAPEPS2    : Scales the repulsive LJ contribution relative to the Morse interactions
!     CAPRHO     : The Morse range parameter
!     CAPRAD     : Radius at the base of a capsomer
!     CAPHEIGHT1 : Height from the basal plane for the primary axial site
!     CAPHEIGHT2 : Height from the basal plane for the secondary axial site

      CAPRADP  = CAPRAD
      CAPRADH  = 2.D0*CAPRAD*DSIN(0.8D0*DATAN(1.D0)) ! CAPRADH = CAPRADP*2.D0*SIN(PI/5)

!     LJ range parameter: calculated from the input parameter in the units of R_{e}
      SIGMASQP = (1.D0 + CAPRADP*SQRT(0.5D0*(5.D0 + SQRT(5.D0))))**2
      SIGMASQH = (1.D0 + CAPRADH*SQRT(0.5D0*(5.D0 + SQRT(5.D0))))**2
      SIGMASQM = (0.5D0*(SQRT(SIGMASQP) + SQRT(SIGMASQH)))**2

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS


      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI(:) = X(J3-2:J3)
         P(:)  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         IF (J1 <= NPS) THEN                     ! NPS -> Number of pentamers
!     Site coordinates for the pentamers
            DO J2 = 1, NRBSITES1                 ! NRBSITES1 -> Number of sites defining pentamers
               J4        = NRBSITES1*(J1-1) + J2
               R(J4,:)   = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (GTEST) THEN
                 DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
                 DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
                 DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))
               ENDIF
            ENDDO
         ELSE
!     Site coordinates for the pentamers
!     NRBSITES -> Number of sites defining pentamers + Number of sites defining hexamers
            DO J2 = NRBSITES1 + 1, NRBSITES     
               J4        = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J1-NPS-1) + J2 - NRBSITES1
               R(J4,:)   = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (GTEST) THEN
                  DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
                  DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
                  DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO J1 = 1, REALNATOMS - 1 
         J3 = 3*J1
         J5 = OFFSET + J3
         IF (J1 <= NPS) THEN
            NRBSI1 = 1
            NRBSI2 = NRBSITES1 
         ELSE 
            NRBSI1 = NRBSITES1 + 1
            NRBSI2 = NRBSITES 
         ENDIF
         DO I = NRBSI1, NRBSI2
            IF (J1 <= NPS) THEN
               J7 = NRBSITES1*(J1-1) + I
            ELSE
               J7 = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J1-NPS-1) + I - NRBSITES1
            ENDIF
            DO J2 = J1 + 1, REALNATOMS
               IF (J2 <= NPS) THEN
                  NRBSJ1 = 1
                  NRBSJ2 = NRBSITES1 
               ELSE
                  NRBSJ1 = NRBSITES1 + 1
                  NRBSJ2 = NRBSITES
               ENDIF
               J4 = 3*J2
               J6 = OFFSET + J4
               DO J = NRBSJ1, NRBSJ2
                  IF (J2 <= NPS) THEN
                     J8 = NRBSITES1*(J2-1) + J
                  ELSE
                     J8 = NRBSITES1*NPS + (NRBSITES - NRBSITES1)*(J2-NPS-1) + J - NRBSITES1
                  ENDIF
                  RSS(:) = R(J7,:) - R(J8,:)
                  IF ((I > NRBSI2-2) .AND. (J > NRBSJ2-2)) THEN  
!     Repulsive contribution between two apex sites: primary-primary and primary-secondary, but not between secondary sites
                     IF (I == NRBSI2 .AND. J == NRBSJ2) CYCLE
                     R2     = 1.D0/DOT_PRODUCT(RSS,RSS)
                     IF (J1 <= NPS .AND. J2 <= NPS) THEN
                        SIGMASQ = SIGMASQP
                     ELSEIF (J1 > NPS .AND. J2 > NPS) THEN
                        SIGMASQ = SIGMASQH
                     ELSE
                        SIGMASQ = SIGMASQM
                     ENDIF
                     VSS    = CAPEPS2*(SIGMASQ*R2)**6
                     ENERGY = ENERGY + VSS

                     IF (GTEST) THEN
                        DVDR =-12.D0*VSS*R2
                        G(J3-2:J3)  = G(J3-2:J3) + DVDR*RSS(:)
                        G(J4-2:J4)  = G(J4-2:J4) - DVDR*RSS(:)

                        G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
                        G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
                        G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))
                        G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
                        G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
                        G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))
                     ENDIF
!     Morse interaction between basal sites          
                  ELSE IF ((I < NRBSI2-1) .AND. (J < NRBSJ2-1)) THEN     ! change here: < to <= in both cases for the previous model
                        R2     = DOT_PRODUCT(RSS,RSS)
                        ABSR   = SQRT(R2)
                        R2     = 1.D0/R2
                        VSS    = EXP(CAPRHO*(1.D0 - ABSR))
                        ENERGY = ENERGY + (1.D0 - VSS)*(1.D0 - VSS) - 1.D0
                        IF (GTEST) THEN
                           DVDR = 2.D0*CAPRHO*(1.D0 - VSS)*VSS/ABSR
                           G(J3-2:J3)  = G(J3-2:J3) + DVDR*RSS(:)
                           G(J4-2:J4)  = G(J4-2:J4) - DVDR*RSS(:)

                           G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
                           G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
                           G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))
                           G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
                           G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
                           G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))
                        ENDIF
                     ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE CAPBIN
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCAPBIN()

      USE COMMONS, ONLY: SITE, NRBSITES, NRBSITES1, CAPRAD, CAPHEIGHT1, CAPHEIGHT2

      IMPLICIT NONE
      INTEGER          :: J1, J2
      DOUBLE PRECISION :: PI, SIN36, CAPRADH 

!     The parametric equation: 
!     x = r * cos(theta), y = r * sin(theta), where 0 <= theta < 2 * PI 
!     can generate the vertices of an n-side regular polygon by assigning n different values to 
!     theta according to 2 * PI * i / n, where 0 <= i < n.  

      PI    = 4.D0*DATAN(1.D0)
      SIN36 = DSIN(PI/5.D0)
      SITE(:,3) = 0.D0

      DO J1 = 1, NRBSITES1 - 2
         SITE(J1,1) = CAPRAD*DSIN(0.4D0*PI*FLOAT(J1-1))
         SITE(J1,2) =-CAPRAD*DCOS(0.4D0*PI*FLOAT(J1-1))
      ENDDO
      SITE(NRBSITES1-1,:) = (/0.D0, 0.D0, CAPHEIGHT1/)
      SITE(NRBSITES1,:)   = (/0.D0, 0.D0,-CAPHEIGHT2/)

      IF (NRBSITES == 7) RETURN

      CAPRADH = 2.D0*CAPRAD*SIN36

      DO J1 = NRBSITES1 + 1, NRBSITES - 2
         SITE(J1,1) = CAPRADH*DSIN(PI*FLOAT(J1-NRBSITES1-1)/3.D0)
         SITE(J1,2) =-CAPRADH*DCOS(PI*FLOAT(J1-NRBSITES1-1)/3.D0)
      ENDDO
      SITE(NRBSITES-1,:) = (/0.D0, 0.D0, CAPHEIGHT1/)
      SITE(NRBSITES,:)   = (/0.D0, 0.D0,-CAPHEIGHT2/)

      END SUBROUTINE DEFCAPBIN

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWCAPBIN()

      USE COMMONS, ONLY: NPS, NRBSITES, NRBSITES1, NTSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='capbin.xyz', STATUS='UNKNOWN')
      GTEST = .FALSE.

      DO J1 = 1, NSAVE
         WRITE(26,'(I6)') NTSITES
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, QMINNATOMS(J1)/2
            J5   = 3*J3
            J7   = 3*QMINNATOMS(J1)/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            IF (J3 <= NPS) THEN 
               DO J2 = 1, NRBSITES1
                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
                  IF (J2 <= NRBSITES1 - 2) THEN
                     WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE IF (J2 == NRBSITES1 - 1) THEN
                     WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF
               ENDDO
            ELSE 
               DO J2 = NRBSITES1 + 1, NRBSITES
                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
                  IF (J2 <= NRBSITES - 2) THEN
                     WRITE(26,'(A4,3F20.10)') 'N', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE IF (J2 == NRBSITES - 1) THEN
                     WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWCAPBIN
