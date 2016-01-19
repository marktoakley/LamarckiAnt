      SUBROUTINE  LWOTP (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY : NATOMS, NRBSITES, BOXLX, BOXLY, BOXLZ, LWCNSTA, LWCNSTB, LWRCUTSQ, LWRCUT2SQ, SITE

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RI(3), RJ(3), DR(3), RIJSS(3), PI(3)
      DOUBLE PRECISION :: D2CMCM, DSS2, R2, R6, R12, VSS, DVSSDR, ENERGY, BOXL(3) 
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      LOGICAL          :: GTEST

!      BOXL = (/7.577639751552796D0, 7.483229813664598D0, 7.284472049689441D0/)
       
!     FROM INPUT PARAMETERS

      BOXL(:)    = (/BOXLX, BOXLY, BOXLZ/)
      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      ENERGY = 0.D0
      IF (GTEST) G(:) = 0.D0

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         PI = X(J5-2:J5)

         CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = MATMUL(RMI,SITE(J2,:))
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
         RI = X(J3-2:J3)

         DO J2 =  J1 + 1, REALNATOMS 

            J4      = 3*J2
            J6      = OFFSET + J4
            RJ      = X(J4-2:J4)
            
            DR(:)  = RI(:) - RJ(:) - ANINT((RI(:) - RJ(:))/BOXL(:))*BOXL(:)
!            DR(:)  = RI(:) - RJ(:) - ANINT((RI(:) - RJ(:))/BOXLX)*BOXLX

            D2CMCM  = DOT_PRODUCT(DR(:),DR(:))

            IF (D2CMCM < LWRCUT2SQ) THEN

!            IF (D2CMCM < 0.25D0*BOXL(3)*BOXL(3)) THEN

               DO I = 1, NRBSITES

                  J7 = NRBSITES*(J1-1) + I
     
                  DO J = 1, NRBSITES

                     J8       = NRBSITES*(J2-1) + J
                     RIJSS(:) = DR(:) + R(J7,:) - R(J8,:)  
                     DSS2     = DOT_PRODUCT(RIJSS(:),RIJSS(:))

                     IF (DSS2 < LWRCUTSQ) THEN

                        R2     = 1.D0 / DSS2
                        R6     = R2 * R2 * R2
                        R12    = R6 ** 2
                        VSS    = R12 - R6 + LWCNSTA * DSS2 + LWCNSTB
                        ENERGY = ENERGY + VSS

                        IF (GTEST) THEN

                           DVSSDR = -(6.D0 * R12 - 3.D0 * R6)/DSS2 + LWCNSTA

                           G(J3-2:J3) = G(J3-2:J3) + DVSSDR*RIJSS(:)
                           G(J4-2:J4) = G(J4-2:J4) - DVSSDR*RIJSS(:)

                           G(J5-2) = G(J5-2) + DVSSDR*DOT_PRODUCT(RIJSS,DR1(J7,:))
                           G(J5-1) = G(J5-1) + DVSSDR*DOT_PRODUCT(RIJSS,DR2(J7,:))
                           G(J5)   = G(J5)   + DVSSDR*DOT_PRODUCT(RIJSS,DR3(J7,:))

                           G(J6-2) = G(J6-2) - DVSSDR*DOT_PRODUCT(RIJSS,DR1(J8,:))
                           G(J6-1) = G(J6-1) - DVSSDR*DOT_PRODUCT(RIJSS,DR2(J8,:))
                           G(J6)   = G(J6)   - DVSSDR*DOT_PRODUCT(RIJSS,DR3(J8,:))

                        ENDIF

                     ENDIF

                  ENDDO

               ENDDO

            ENDIF

         ENDDO

      ENDDO

      ENERGY = 4.D0*ENERGY  
      IF (GTEST) G(:) = 8.D0*G(:)

      END SUBROUTINE LWOTP 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFLWOTP()
    
      USE COMMONS !, ONLY : LWRCUT, LWCNSTA, LWCNSTB, LWRCUTSQ, LWRCUT2SQ, SITE(:,:)   
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: PI, RCUT6, RCUT12, DIST(NRBSITES), RCUT2

      LWRCUTSQ = LWRCUT * LWRCUT
      RCUT6    = (1.D0 / LWRCUT) ** 6
      RCUT12   = RCUT6 * RCUT6
      LWCNSTA  = (6.D0 * RCUT12 - 3.D0 * RCUT6) / LWRCUTSQ
      LWCNSTB  = 4.D0 * RCUT6 - 7.D0 * RCUT12

      PI       = 4.D0 * DATAN(1.D0)

      SITE(1,1) = 0.D0
      SITE(1,2) = - 2.D0 * DSIN(7.D0*PI/24.D0) / 3.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = DCOS(7.D0*PI/24.D0)
      SITE(2,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      SITE(2,3) = 0.D0

      SITE(3,1) = - DCOS(7.D0*PI/24.D0)
      SITE(3,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      SITE(3,3) = 0.D0

      DO I = 1, NRBSITES

         DIST(I) = DSQRT(DOT_PRODUCT(SITE(I,:),SITE(I,:)))

      ENDDO

      DELRC      = 2.D0 * MAXVAL(DIST)
      RCUT2      = LWRCUT + DELRC
      LWRCUT2SQ  = RCUT2 * RCUT2

      END SUBROUTINE DEFLWOTP

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE  LWOTPC (X, G, ENERGY, GTEST)

      USE COMMONS !, ONLY : NATOMS, NRBSITES, SITE(:,:)

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, R2, R6, R12
      DOUBLE PRECISION :: RI(3), RSS(3), P(3), R(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
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

            DO I = 1, NRBSITES 

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES 

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
                  R6     = R2*R2*R2
                  R12    = R6*R6
                  ENERGY = ENERGY + R12 - R6 

                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     DVDR       = -(6.D0*R12 - 3.D0*R6)*R2
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

      ENDDO

      ENERGY = 4.D0*ENERGY
      G(:)   = 8.D0*G(:)

      END SUBROUTINE LWOTPC 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE TAKESTEPLWOTP (NP)

!     THIS ROUTINE TAKES STEP FOR LWOTP MOLECULE ENSURING NO OVERLAP BETWEEN LJSITES

      USE COMMONS 

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2, BOXL(3)
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI, THETA2, FRPISQ
      DOUBLE PRECISION :: RI(3), RJ(3), RISITE(3), RJSITE(3), DR(3), DSS(3)
      DOUBLE PRECISION :: P1(3), P2(3), RIJSQ, SSCSQ, RMI(3,3), RMJ(3,3) 

      BOXL = (/7.577639751552796D0, 7.483229813664598D0, 7.284472049689441D0/)

      PI         = 4.D0*ATAN(1.0D0)
      FRPISQ     = 4.D0*PI*PI
      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS
      SSCSQ      = (1.D0 + DELRC) * (1.D0 + DELRC)

      DO J1 = 1,REALNATOMS

         J2     = 3*J1
         DUMMY2 = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2
         IF (DUMMY2 .GT. RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), &
                                       COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'initial coordinate outside container -- increase container radius'
            STOP
         END IF

      END DO

      DO J1 = 1,3*NATOMS
         COORDSO(J1,NP) = COORDS(J1,NP)
      END DO

      DO J1 = 1,NATOMS/2
         VATO(J1,NP) = VAT(J1,NP)
      END DO

!     FIND THE CENTRE OF MASS

      XMASS = 0.0D0
      YMASS = 0.0D0
      ZMASS = 0.0D0

      DO J1 = 1,NATOMS/2

         XMASS = XMASS + COORDS(3*(J1-1)+1,NP)
         YMASS = YMASS + COORDS(3*(J1-1)+2,NP)
         ZMASS = ZMASS + COORDS(3*(J1-1)+3,NP)

      ENDDO

      XMASS = XMASS/(REALNATOMS)
      YMASS = YMASS/(REALNATOMS)
      ZMASS = ZMASS/(REALNATOMS)

!     Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!     and the pair energy of the most tightly bound atom, VMIN. An angular step is
!     taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!     DMAX (or CMMAX from CM of the cluster).

      DMAX  =  1.0D0
      VMAX  = -1.0D3
      VMAX2 = -1.0D3
      VMIN  =  0.0D0
      CMMAX =  1.0D0

      DO J1 = 1, REALNATOMS

         J2 = 3*J1
         DIST(J1)   = DSQRT( COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2)
         CMDIST(J1) = DSQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1) .GT. CMMAX) CMMAX = CMDIST(J1)
         IF (DIST(J1) .GT. DMAX) DMAX = DIST(J1)
         IF (VAT(J1,NP) .GT. VMAX) THEN
            VMAX = VAT(J1,NP)
            JMAX = J1
         ELSE IF ((VAT(J1,NP).LT. VMAX) .AND. (VAT(J1,NP) .GT. VMAX2)) THEN
              VMAX2 = VAT(J1,NP)
              JMAX2 = J1
         ENDIF
         IF (VAT(J1,NP) .LT. VMIN) VMIN = VAT(J1,NP)

      ENDDO

      IF (VAT(JMAX,NP) > (ASTEP(NP)*VMIN) .AND. (.NOT.NORESET)) THEN

         J2 = 3*JMAX
         THETA           = DPRAND()*PI
         PHI             = DPRAND()*PI*2.0D0
         COORDS(J2-2,NP) = XMASS + (CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
         COORDS(J2-1,NP) = YMASS + (CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
         COORDS(J2,NP)   = ZMASS + (CMMAX+1.0D0)*DCOS(THETA)
         DUMMY           = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY > RADIUS) THEN

            DUMMY           = DSQRT(RADIUS*0.99D0/DUMMY)
            COORDS(J2-2,NP) = COORDS(J2-2,NP)*DUMMY
            COORDS(J2-1,NP) = COORDS(J2-1,NP)*DUMMY
            COORDS(J2,NP)   = COORDS(J2,NP)*DUMMY

         END IF

      ENDIF

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3

!     CHECK FOR OVERLAP

         OVRLPT = .TRUE.

95       DO WHILE (OVRLPT)

            LOCALSTEP = 0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2,NP) = COORDS(J3-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1,NP) = COORDS(J3-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3,NP)   = COORDS(J3,NP) + LOCALSTEP*RANDOM

            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-2,NP) = COORDS(J5-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-1,NP) = COORDS(J5-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5,NP)   = COORDS(J5,NP) + LOCALSTEP*RANDOM
          
            OVRLPT = .FALSE.

            RI(:)  = COORDS(J3-2:J3,NP)
            P1(:)  = COORDS(J5-2:J5,NP)

!     ROTATION MATRIX

            CALL ROTMAT (P1, RMI)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE
 
               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4,NP)
               P2(:) = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

               CALL ROTMAT (P2, RMJ)

               IF (PERIODIC) THEN

                  DR(:)  = RI(:) - RJ(:) - ANINT((RI(:) - RJ(:))/BOXL(:))*BOXL(:)

               ELSE

                  DR(:)  = RI(:) - RJ(:)

               ENDIF

               RIJSQ  = DOT_PRODUCT(DR,DR)
                        
               IF (RIJSQ <= SSCSQ) THEN

                  DO I = 1, NRBSITES

                     RISITE = MATMUL(RMI,SITE(I,:)) 

                     DO J = 1, NRBSITES

                        RJSITE = MATMUL(RMJ,SITE(J,:))

                        DSS = DR + RISITE - RJSITE

                        IF (DOT_PRODUCT(DSS,DSS) < 1.D0) THEN

                           OVRLPT = .TRUE.
                           GO TO 95
                        ENDIF

                     ENDDO

                  ENDDO

               ENDIF

            ENDDO  ! END LOOP OVER J2

         ENDDO  ! END WHILE
 
      ENDDO  ! END LOOP OVER J1   

      DO I = REALNATOMS + 1, NATOMS

         J       = 3*I
         THETA2  = DOT_PRODUCT(COORDS(J-2:J,NP),COORDS(J-2:J,NP))
         IF (THETA2 > FRPISQ) THEN
            THETA2   = DSQRT(THETA2)
            THETA    = THETA2 - INT(THETA2/(2.D0*PI))*2.D0*PI
            COORDS(J-2:J,NP) = COORDS(J-2:J,NP)/THETA2 * THETA
         ENDIF

      ENDDO
      
      END SUBROUTINE TAKESTEPLWOTP

!     ----------------------------------------------------------------------------------------------

! jdf43>        Subroutine ROTMAT now lives in rigidbaa.f90 30/01/12

!     SUBROUTINE ROTMAT (P, RM)

!     IMPLICIT NONE

!     DOUBLE PRECISION :: E(3,3), I3(3,3), RM(3,3), P(3), THETA, THETA2 

!     I3(:,:) = 0.D0
!     I3(1,1) = 1.D0
!     I3(2,2) = 1.D0
!     I3(3,3) = 1.D0
       
!     THETA   = DSQRT(DOT_PRODUCT(P,P))

!     IF (THETA == 0.D0) THEN

!        RM = I3

!     ELSE

!        THETA2  = THETA * THETA
!        E(:,:)  = 0.D0
!        E(1,2)  = -P(3)
!        E(1,3)  =  P(2)
!        E(2,3)  = -P(1)
!        E(2,1)  = -E(1,2)
!        E(3,1)  = -E(1,3)
!        E(3,2)  = -E(2,3)
!        E       = E/THETA

!        RM = I3 + (1.D0-COS(THETA))*MATMUL(E,E) + E*SIN(THETA)

!     ENDIF

!     END SUBROUTINE ROTMAT 
