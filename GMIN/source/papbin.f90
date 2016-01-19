!     ----------------------------------------------------------------------------------------------

      SUBROUTINE PAPBIN(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, PAPEPS, PAPS, PAPANG1, PAPANG2, YKAPPA 

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, NMOL, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, RIJSQ, ABSRIJ, EXPFCT, PI 
      DOUBLE PRECISION :: RIJ(3), NR(3), P(3), WP, DELR, ANGFAC(2), PAPCD(2), DWPDR
      DOUBLE PRECISION :: DOTI, DOTJ, ARGI, ARGJ, PHII, PHIJ
      DOUBLE PRECISION :: E(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DE1(NATOMS*NRBSITES/2,3), DE2(NATOMS*NRBSITES/2,3), DE3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: DDOTIDR(3), DDOTJDR(3), DPHIIDR(3), DPHIJDR(3), DPHIIDPI(3), DPHIJDPJ(3)
      DOUBLE PRECISION :: LAMBDA, INVS
      LOGICAL          :: GTEST
      DOUBLE PRECISION :: EPS, SIGMA, KAPPA     ! potential parameters

      EPS = 1.D0; SIGMA = 1.D0; KAPPA = YKAPPA     ! Yukawa parameters
      PI     = 4.D0*DATAN(1.D0)
      PAPCD(1) = COS(PI/180.D0*PAPANG1)
      PAPCD(2) = COS(PI/180.D0*PAPANG2)
 
      CALL DEFPAPBIN()

      ANGFAC(1) = PI/(1.D0 - PAPCD(1))
      ANGFAC(2) = PI/(1.D0 - PAPCD(2))
      INVS    = 1.D0/PAPS
      LAMBDA  = SIGMA + 1.D0/KAPPA 
      NMOL    = NATOMS/2
      OFFSET  = 3*NMOL

      ENERGY  = 0.D0
      IF (GTEST) THEN
        G(:) = 0.D0
      ENDIF

      DO J1 = 1, NMOL

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO I = 1, NRBSITES

            J7 = NRBSITES*(J1-1) + I
            E(J7,:) = MATMUL(RMI(:,:),RBSTLA(I,:))

            IF (GTEST) THEN
               DE1(J7,:) = MATMUL(DRMI1(:,:),RBSTLA(I,:))
               DE2(J7,:) = MATMUL(DRMI2(:,:),RBSTLA(I,:))
               DE3(J7,:) = MATMUL(DRMI3(:,:),RBSTLA(I,:))
            ENDIF 

         ENDDO
 
      ENDDO

      DO J1 = 1, NMOL - 1  

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, NMOL

            J4 = 3*J2
            J6 = OFFSET + J4

!     Calculate electrostatic repulsion via Yukawa potential
            RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSRIJ = DSQRT(RIJSQ)
            EXPFCT = DEXP(-KAPPA*(ABSRIJ - SIGMA)) 
            ENERGY = ENERGY + EPS*SIGMA*EXPFCT/ABSRIJ
            
            IF (GTEST) THEN
              DVDR =-EPS*SIGMA*EXPFCT/RIJSQ*(KAPPA + 1.D0/ABSRIJ) 
              G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJ(:)
              G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJ(:)
            ENDIF

!      Calculation of the factor, WP, which depends on ABSRIJ and is independent of the sites

!     Find normalised vector between I and J, and distance parameter DELR
            NR(:)  = RIJ(:)/ABSRIJ
            DELR   = ABSRIJ - LAMBDA

!     Set specific cases for the patch-antipatch potential and its gradient
            IF (DELR > INVS) THEN
               RETURN
            ELSE IF (DELR < 0.D0) THEN
!     Distance less than LAMBDA, so full attraction
               WP =-1.D0
               DWPDR = 0.D0
            ELSE
!     Distance is between LAMBDA and LAMBDA + INVS, so potential varies from -1 to 0 
               WP    =-0.5D0*(1.D0 + COS(PI*DELR*PAPS))
               DWPDR = 0.5D0*PI*PAPS*SIN(PI*DELR*PAPS)
            ENDIF

            DO I = 1, NRBSITES         ! Loop over patches/antipatches on body J1

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8 = NRBSITES*(J2-1) + J
 
!     Find the angles between the orientation of (anti)patch I and NR, and between the orientation of (anti)patch J and NR
!     Negative for DOTI, since RIJ is displacement from J2 to J1
                  DOTI =-DOT_PRODUCT(E(J7,:),NR(:))
                  DOTJ = DOT_PRODUCT(E(J8,:),NR(:))

!     Calculate the values of PHI for the angular potential
                  IF (DOTI < PAPCD(I)) THEN
!     Angle greater than patch width, so no attraction
                     PHII =-1.D0
                  ELSE
!     Angle less than patch width, so some attraction
                     ARGI = ANGFAC(I)*(DOTI - PAPCD(I))
                     PHII =-COS(ARGI) 
                  ENDIF
 
                  IF (DOTJ < PAPCD(J)) THEN
!     Angle greater than patch width, so no attraction
                     PHIJ =-1.D0
                  ELSE
!     Angle less than patch width, so some attraction
                     ARGJ = ANGFAC(J)*(DOTJ - PAPCD(J))
                     PHIJ =-COS(ARGJ) 
                  ENDIF 

!     Add the patch-antipatch attraction to the potential energy

                  ENERGY = ENERGY + 0.25D0*PAPEPS*(1.D0 + PHII)*(1.D0 + PHIJ)*WP

                  IF (GTEST) THEN
!     Need to find the derivatives of the patch-antipatch attraction wrt translational and rotational coordinates
                     IF ((DOTI >= PAPCD(I)) .AND. (DOTJ >= PAPCD(J))) THEN

                        DDOTIDR(:) =-DOTI*RIJ(:)/RIJSQ - E(J7,:)/ABSRIJ
                        DDOTJDR(:) =-DOTJ*RIJ(:)/RIJSQ + E(J8,:)/ABSRIJ

                        DPHIIDR(:) = ANGFAC(I)*SIN(ARGI)*DDOTIDR(:) 
                        DPHIJDR(:) = ANGFAC(J)*SIN(ARGJ)*DDOTJDR(:) 

                        DPHIIDPI(1) =-ANGFAC(I)*SIN(ARGI)*DOT_PRODUCT(NR(:),DE1(J7,:))
                        DPHIIDPI(2) =-ANGFAC(I)*SIN(ARGI)*DOT_PRODUCT(NR(:),DE2(J7,:))
                        DPHIIDPI(3) =-ANGFAC(I)*SIN(ARGI)*DOT_PRODUCT(NR(:),DE3(J7,:))

                        DPHIJDPJ(1) = ANGFAC(J)*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE1(J8,:))
                        DPHIJDPJ(2) = ANGFAC(J)*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE2(J8,:))
                        DPHIJDPJ(3) = ANGFAC(J)*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE3(J8,:))

                        G(J3-2:J3) = G(J3-2:J3) + 0.25D0*PAPEPS*((1.D0 + PHIJ)*WP*DPHIIDR(:)                  &
     &                             + (1.D0 + PHII)*WP*DPHIJDR(:) +(1.D0+PHII)*(1.D0+PHIJ)*DWPDR*NR(:))
                        G(J4-2:J4) = G(J4-2:J4) - 0.25D0*PAPEPS*((1.D0 + PHIJ)*WP*DPHIIDR(:)                  &
     &                             + (1.D0 + PHII)*WP*DPHIJDR(:) +(1.D0+PHII)*(1.D0+PHIJ)*DWPDR*NR(:))

                        G(J5-2:J5) = G(J5-2:J5) + 0.25D0*PAPEPS*(1.D0 + PHIJ)*WP*DPHIIDPI(:)
                        G(J6-2:J6) = G(J6-2:J6) + 0.25D0*PAPEPS*(1.D0 + PHII)*WP*DPHIJDPJ(:)

                     ENDIF

                  ENDIF 

               ENDDO ! End loop over patches in J2

            ENDDO ! End loop over patches in J1

         ENDDO 

      ENDDO 

      END SUBROUTINE PAPBIN

!     ----------------------------------------------------------------------------------------------
 
!     DEFPAPBIN           
                                                                           
!     ----------------------------------------------------------------------------------------------                                                                                         
      SUBROUTINE DEFPAPBIN()

      USE COMMONS

      IMPLICIT NONE

      RBSTLA(1,:)= (/  0.D0,  0.D0, 1.D0/)
      RBSTLA(2,:)= (/  0.D0,  0.D0,-1.D0/)
      
      END SUBROUTINE DEFPAPBIN

!     ----------------------------------------------------------------------------------------------


      SUBROUTINE VIEWPAPBIN()

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='papbin.xyz', STATUS='UNKNOWN')

      GTEST = .FALSE.

      DO J1 = 1, NSAVE

         WRITE(26,'(I4)') (NATOMS/2)*(NRBSITES+1)
         WRITE(26,'(A)') ' '

         DO J3 = 1, NATOMS/2

            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            DO J2 = 1, NRBSITES+1
               IF (J2 <= NRBSITES/2) THEN
                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + 0.5D0*MATMUL(RMI(:,:),RBSTLA(J2,:))
                  WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSEIF (J2 <= NRBSITES) THEN
                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + 0.5D0*MATMUL(RMI(:,:),RBSTLA(J2,:))
                  WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  RBCOORDS(1:3) = QMINP(J1,J5-2:J5) 
                  WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF
            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWPAPBIN
