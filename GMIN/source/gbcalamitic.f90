      SUBROUTINE GBCALAMITIC (X, G, ENERGY, GTEST) 

      USE COMMONS, ONLY : NATOMS, GBSIGNOT, GBEPSNOT, GBMU, GBNU, GBCHI, GBCHIPRM

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: SCSIG, SCSIG3, ENERGY
      DOUBLE PRECISION :: EPS1, EPS2, EPS
      DOUBLE PRECISION :: FCT1, FCT2, FCT3, FCT4, FCT5, FCT6, FCT7, FCT8
      DOUBLE PRECISION :: FCT9, FCT10, FCT11, FCT12, FCT13, FCT14, FCT15
      DOUBLE PRECISION :: FCT16, FCT17, FCT18, FCT19, FCT20
      DOUBLE PRECISION :: FCT3P4, FCT3M4, FCT7P8, FCT7M8
      DOUBLE PRECISION :: ALP, BET, GAM, APB, AMB
      DOUBLE PRECISION :: RIJSQ, R2, ABSRIJ, INVR, SCR 
      DOUBLE PRECISION :: SRM1, SRM2, SRM6, SRM7, SRM12, SRM13, SR12M6
      DOUBLE PRECISION :: VR, VA, VB, VG, FIJ(3), FIJN, FIJEI, FIJEJ 
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), P(3), EI(3), EJ(3), OST(3), RMI(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST

      OST        = (/0.D0, 0.D0, 1.D0/)
      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      ENERGY = 0.D0
      IF (GTEST) G(:)   = 0.D0

      DO J1 = 1, REALNATOMS 

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX AND ITS DERIVATIVES

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST) 

            E(J1,:)     = MATMUL(RMI(:,:),OST(:)) 
            IF (GTEST) THEN
               DE1(J1,:)   = MATMUL(DRMI1(:,:),OST(:))
               DE2(J1,:)   = MATMUL(DRMI2(:,:),OST(:))
               DE3(J1,:)   = MATMUL(DRMI3(:,:),OST(:))
            ENDIF
           
      ENDDO

      DO J1 = 1, REALNATOMS - 1

         J3      = 3*J1
         J5      = OFFSET + J3
         RI      = X(J3-2:J3)

!     BEGIN INNER LOOP OVER PARTICLES

         DO J2 = J1 + 1, REALNATOMS

            J4     = 3*J2
            J6     = OFFSET + J4
            RJ     = X(J4-2:J4) 

!     CALCULATE SEPARATION

            RIJ    = RI - RJ
            RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
            ABSRIJ = DSQRT(RIJSQ)
            NR     = RIJ / ABSRIJ
            SCR    = ABSRIJ/GBSIGNOT
            INVR   = 1.D0/ABSRIJ
            R2     = 1.D0/RIJSQ

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

            EI  = E(J1,:) 
            EJ  = E(J2,:)
            ALP = DOT_PRODUCT(NR,EI) 
            BET = DOT_PRODUCT(NR,EJ)
            GAM = DOT_PRODUCT(EI,EJ)

!     CALCULATE USEFUL QUANTITIES

            APB    = ALP+BET
            AMB    = ALP-BET

            FCT1   = 1.D0/(1.D0+GBCHI*GAM)
            FCT2   = 1.D0/(1.D0-GBCHI*GAM)
            FCT3   = APB*FCT1
            FCT4   = AMB*FCT2
            FCT3P4 = FCT3+FCT4
            FCT3M4 = FCT3-FCT4
            FCT5   = 1.D0/(1.D0+GBCHIPRM*GAM)
            FCT6   = 1.D0/(1.D0-GBCHIPRM*GAM)
            FCT7   = (ALP+BET)*FCT5
            FCT8   = (ALP-BET)*FCT6
            FCT7P8 = FCT7+FCT8
            FCT7M8 = FCT7-FCT8

!     CALCULATE $\epsilon$

            EPS1   = DSQRT(FCT1*FCT2)
            EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
            EPS    = GBEPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

            SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
            SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE DEL(V)/DEL(R)

            SRM1   = 1.D0/(SCR-SCSIG+1.D0)
            SRM2   = SRM1*SRM1
            SRM6   = SRM2*SRM2*SRM2
            SRM12  = SRM6*SRM6
            SR12M6 = SRM12-SRM6

!     CALCULATE ENERGY

            ENERGY = ENERGY + EPS*SR12M6

            IF (GTEST) THEN

               SRM7   = SRM6*SRM1
               SRM13  = SRM12*SRM1
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/GBSIGNOT)*EPS*FCT9

!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)

               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(GBSIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*SCSIG3

               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9

               VG     = FCT19 - FCT20

!     CALCULATE CONTRIBUTION TO FORCES

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FIJ  = FIJN*NR + FIJEI*EI + FIJEJ*EJ

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

               G(J3-2:J3) = G(J3-2:J3) + FIJ 
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

            ENDIF

         ENDDO

      ENDDO

!     CALCULATE THE GAY-BERNE POTENTIAL

      ENERGY = 4.D0*ENERGY

      END SUBROUTINE GBCALAMITIC 

!     --------------------------------------------------------------------------

      SUBROUTINE TKSTDCELPSD (NP)

!     THIS ROUTINE TAKES STEP FOR SINGLE-SITE ELLIPSOIDAL BODIES ENSURING NO OVERLAP

      USE COMMONS

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET, PTINDX
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI, THETA2, FRPISQ
      DOUBLE PRECISION :: AE(3,3), BE(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3)
      DOUBLE PRECISION :: P1(3), P2(3), ABSRIJ, RCUT, XMIN, FMIN, ECFVAL

      PI         = 4.D0*ATAN(1.0D0)
      FRPISQ     = 4.D0*PI*PI

      IF (GBT .OR. GBDT .OR. GBDPT) THEN

         RCUT       = GBEPSNOT*MAX(GBKAPPA,1.D0)

      ENDIF

      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS

      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP)) ! COORDS might have been shifted by symmetry

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

            CALL ELPSDRTMT (P1, ESA, AE)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE

               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4,NP)
               P2(:) = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

               CALL ELPSDRTMT (P2, ESA, BE)

               RIJ    = RI - RJ
               ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))

               IF (ABSRIJ < RCUT) THEN

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                  ECFVAL = - FMIN

                  IF (ECFVAL < 1.D0) THEN
                     OVRLPT = .TRUE.
                     GO TO 95
                  ENDIF

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

      END SUBROUTINE TKSTDCELPSD
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ELPSDRTMT (P, ESA, A)

      IMPLICIT NONE

      INTEGER          :: K
      DOUBLE PRECISION :: E(3,3), I3(3,3), ROT2(3,3), ROTMAT(3,3), P(3), THETA, THETA2 
      DOUBLE PRECISION :: AEZR1(3,3), A(3,3), ESA(3)

      I3(:,:)    = 0.D0
      AEZR1(:,:) = 0.D0
      I3(1,1)    = 1.D0
      I3(2,2)    = 1.D0
      I3(3,3)    = 1.D0
       
      THETA   = DSQRT(DOT_PRODUCT(P,P))

      IF (THETA == 0.D0) THEN

         ROTMAT = I3

      ELSE

         THETA2  = THETA * THETA
         E(:,:)  = 0.D0
         E(1,2)  = -P(3)
         E(1,3)  =  P(2)
         E(2,3)  = -P(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         E       = E/THETA

         ROTMAT = I3 + (1.D0-COS(THETA))*MATMUL(E,E) + E*SIN(THETA)

      ENDIF

      DO K = 1, 3
         AEZR1(K,K) = 1.D0/(ESA(K)*ESA(K))
      ENDDO

      A = MATMUL(ROTMAT,(MATMUL(AEZR1,(TRANSPOSE(ROTMAT)))))

      END SUBROUTINE ELPSDRTMT 
