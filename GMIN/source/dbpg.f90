      SUBROUTINE DUMBBELLP (X, G, ENERGY, GTEST)

!      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, DBEPSBB, DBEPSAB, DBSIGBB, DBSIGAB, DBPMU, EFIELDT, EFIELD, DBYUKAWAT, LAMBDAYAA, LAMBDAYBB, LAMBDAYAB, YEPSFAC, RESTRAINLT, RESTRAINLK, RESTRAINLDIST
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, DBEPSBB, DBEPSAB, DBSIGBB, DBSIGAB, DBPMU, EFIELDT, EFIELD, DBYUKAWAT, LAMBDAYAA, LAMBDAYBB, LAMBDAYAB, YEPSFAC



      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, DIST, DR, R2, R6, R12, ABSRIJ, RIJSQ, DPFCT
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RSS(3), NR(3), P(3), DU(3), EI(3), EJ(3), R(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS
      DOUBLE PRECISION :: DBSIGAA, DBEPSAA
      LOGICAL          :: GTEST
      
      IF (DBYUKAWAT .EQV. .TRUE. ) THEN
         CALL DEFDUMYUKAWA(DU, DPFCT) 
         DBSIGAA  = 1.0D0
         DBEPSAA  = 1.0D0     
         DBEPSAB = SQRT(DBEPSAA*DBEPSBB)        
         DBSIGAB = 0.5D0*(DBSIGAA+DBSIGBB)        
      ELSE
         CALL DEFDUM(DU, DPFCT, CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS)
      ENDIF

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

         E(J1,:)   = MATMUL(RMI(:,:),DU(:))
         DE1(J1,:) = MATMUL(DRMI1(:,:),DU(:))
         DE2(J1,:) = MATMUL(DRMI2(:,:),DU(:))
         DE3(J1,:) = MATMUL(DRMI3(:,:),DU(:))

         DO J2 = 1, NRBSITES - 1

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))

            IF (GTEST) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS !- 1 

         J3     = 3*J1
         J5     = OFFSET + J3
         RI(:)  = X(J3-2:J3) 
         EI(:)  = E(J1,:)

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     LJ CONTRIBUTION

            DO I = 1, NRBSITES - 1

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES - 1

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
                  R6     = R2*R2*R2
                  R12    = R6*R6

                  DIST = 1.0D0/SQRT(R2)

                  IF (I == 1 .AND. J == 1) THEN
                     IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                        ENERGY = ENERGY + YEPSFAC * DBEPSAA * 1/DIST * DBSIGAA * EXP(-1.0D0*(DIST-DBSIGAA)/LAMBDAYAA)  
                     ELSE
                        ENERGY = ENERGY + (CLJ12BB*R6 - CLJ6BB)*R6
                     ENDIF
                  ELSEIF (I == 2 .AND. J == 2) THEN
                     IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                        ENERGY = ENERGY + YEPSFAC * DBEPSBB * 1/DIST * DBSIGBB * EXP(-1.0D0*(DIST-DBSIGBB)/LAMBDAYBB)  
                     ELSE
                        ENERGY = ENERGY + (CLJ12SS*R6 - CLJ6SS)*R6
                     ENDIF
                  ELSE
                     IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                        ENERGY = ENERGY + YEPSFAC * DBEPSAB * 1/DIST * DBSIGAB * EXP(-1.0D0*(DIST-DBSIGAB)/LAMBDAYAB)  
                     ELSE
                        ENERGY = ENERGY + (CLJ12BS*R6 - CLJ6BS)*R6
                     ENDIF
                  ENDIF
 
                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     IF (I == 1 .AND. J == 1) THEN
                        IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                           DVDR   = - YEPSFAC * DBEPSAA * DBSIGAA * (1.0D0/DIST**3 + 1.0D0/LAMBDAYAA/DIST**2) * EXP(-1.0D0*(DIST-DBSIGAA)/LAMBDAYAA) 
                        ELSE
                           DVDR   = -6.D0*(2.D0*CLJ12BB*R6 - CLJ6BB)*R6*R2
                        ENDIF
                     ELSEIF (I == 2 .AND. J == 2) THEN
                        IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                           DVDR   = - YEPSFAC * DBEPSBB * DBSIGBB * (1.0D0/DIST**3 + 1.0D0/LAMBDAYBB/DIST**2) * EXP(-1.0D0*(DIST-DBSIGBB)/LAMBDAYBB) 
                        ELSE
                           DVDR   = -6.D0*(2.D0*CLJ12SS*R6 - CLJ6SS)*R6*R2
                        ENDIF
                     ELSE
                        IF (DBYUKAWAT .EQV. .TRUE. ) THEN
                           DVDR   = - YEPSFAC * DBEPSAB * DBSIGAB * (1.0D0/DIST**3 + 1.0D0/LAMBDAYAB/DIST**2) * EXP(-1.0D0*(DIST-DBSIGAB)/LAMBDAYAB) 
                        ELSE
                           DVDR   = -6.D0*(2.D0*CLJ12BS*R6 - CLJ6BS)*R6*R2
                        ENDIF
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

!            RI(:)  = X(J3-2:J3) 
            RJ(:)  = X(J4-2:J4)
            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ
!            EI(:)  = E(J1,:)
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


! hk286
!      IF (RESTRAINLT) THEN
!         J1 = 1
!         J2 = 2
!         DR = (X(3*J1-2) - X(3*J2-2))**2 + (X(3*J1-1) - X(3*J2-1))**2 + (X(3*J1) - X(3*J2))**2 
!         DR = DSQRT(DR)
!         G(3*J1)   = G(3*J1)   + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J1)   - X(3*J2)  ) * RESTRAINLK     
!         G(3*J1-1) = G(3*J1-1) + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J1-1) - X(3*J2-1)) * RESTRAINLK
!         G(3*J1-2) = G(3*J1-2) + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J1-2) - X(3*J2-2)) * RESTRAINLK
         
!         G(3*J2)   = G(3*J2)   + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J2)   - X(3*J1)  ) * RESTRAINLK      
!         G(3*J2-1) = G(3*J2-1) + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J2-1) - X(3*J1-1)) * RESTRAINLK
!         G(3*J2-2) = G(3*J2-2) + 1.0 * (DR - RESTRAINLDIST)/DR * (X(3*J2-2) - X(3*J1-2)) * RESTRAINLK
         
!         ENERGY = ENERGY + 0.5D0* (DR - RESTRAINLDIST)**2 * RESTRAINLK
!      END IF


      END SUBROUTINE DUMBBELLP 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDUM(DU, DPFCT, CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS)
    
      USE COMMONS, ONLY: SITE, DBEPSBB, DBSIGBB, DBPMU  
     
      IMPLICIT NONE

      DOUBLE PRECISION :: DBEPSAA, DBEPSAB, DBSIGAA, DBSIGAB, DPFCT, DU(3)
      DOUBLE PRECISION :: CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS

      DBSIGAA  = 1.D0
      DBEPSAA  = 1.D0
      DBEPSAB = SQRT(DBEPSAA*DBEPSBB)        
      DBSIGAB = 0.5D0*(DBSIGAA+DBSIGBB)        

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.5D0*DBSIGAA

      SITE(2,1) = 0.D0
      SITE(2,2) = 0.D0
      SITE(2,3) = -0.5D0*DBSIGBB

      SITE(3,1) = 0.D0
      SITE(3,2) = 0.D0
      SITE(3,3) = 0.D0

      DU(:) = (/0.D0, 1.D0, 0.D0/)
      DPFCT = 3.D0*DBPMU*DBPMU*DBSIGAA**3.D0

      CLJ6BB  = 4.D0*DBEPSAA*DBSIGAA**6 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      CLJ12BB = 4.D0*DBEPSAA*DBSIGAA**12
      CLJ6SS  = 4.D0*DBEPSBB*DBSIGBB**6
      CLJ12SS = 4.D0*DBEPSBB*DBSIGBB**12
      CLJ6BS  = 4.D0*DBEPSAB*DBSIGAB**6
      CLJ12BS = 4.D0*DBEPSAB*DBSIGAB**12

      END SUBROUTINE DEFDUM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDUMYUKAWA(DU, DPFCT) 
    
      USE COMMONS, ONLY: SITE, DBEPSBB, DBSIGBB, DBPMU  
     
      IMPLICIT NONE

      DOUBLE PRECISION :: DBEPSAA, DBEPSAB, DBSIGAA, DBSIGAB, DPFCT, DU(3)

      DBSIGAA  = 1.D0
      DBEPSAA  = 1.D0
      DBEPSAB = SQRT(DBEPSAA*DBEPSBB)        
      DBSIGAB = 0.5D0*(DBSIGAA+DBSIGBB)        

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.5D0*DBSIGAA

      SITE(2,1) = 0.D0
      SITE(2,2) = 0.D0
      SITE(2,3) = -0.5D0*DBSIGBB

      SITE(3,1) = 0.D0
      SITE(3,2) = 0.D0
      SITE(3,3) = 0.D0

      DU(:) = (/0.D0, 1.D0, 0.D0/)
      DPFCT = 3.D0*DBPMU*DBPMU*DBSIGAA**3.D0

      END SUBROUTINE DEFDUMYUKAWA


!     ----------------------------------------------------------------------------------------------

!      SUBROUTINE TAKESTEPDB (NP)

!     THIS ROUTINE TAKES STEP FOR LWOTP MOLECULE ENSURING NO OVERLAP BETWEEN LJSITES

!      USE COMMONS, ONLY: DBYUKAWAT, LAMBDAYAA, LAMBDAYBB, LAMBDAYAB

!      IF (DBYUKAWAT .EQ. .TRUE.) THEN
!         LAMBDAYAA = LAMBDAYAA + 0.001
!         LAMBDAYAB = LAMBDAYAB + 0.001
!         LAMBDAYBB = LAMBDAYBB + 0.001
!      ENDIF
!      PRINT *,LAMBDAYAA, LAMBDAYAB, LAMBDAYBB
      
!      END SUBROUTINE TAKESTEPDB


      SUBROUTINE TAKESTEPDB (NP)

!     THIS ROUTINE TAKES STEP FOR LWOTP MOLECULE ENSURING NO OVERLAP BETWEEN LJSITES

      USE COMMONS 

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI, THETA2, FRPISQ
      DOUBLE PRECISION :: RI(3), RJ(3), RISITE(3), RJSITE(3), DR(3), DSS(3)
      DOUBLE PRECISION :: P1(3), P2(3), RIJSQ, SSCSQ, RMI(3,3), RMJ(3,3), DRM(3,3) 

      PI         = 4.D0*ATAN(1.0D0)
      FRPISQ     = 4.D0*PI*PI
      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS
      SSCSQ      = (1.D0 + DELRC)*(1.D0 + DELRC)

      IF (CENT) CALL CENTRE2(COORDS(1:3*NATOMS,NP)) 

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
!         write(*,*) j1, ovrlpt

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

            CALL RMDRVT (P1, RMI, DRM, DRM, DRM, .FALSE.)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE
 
               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4,NP)
               P2(:) = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

               CALL RMDRVT (P2, RMJ, DRM, DRM, DRM, .FALSE.)

               DR(:)  = RI(:) - RJ(:)

               RIJSQ  = DOT_PRODUCT(DR,DR)
                        
               IF (RIJSQ <= 4.D0) THEN

                  DO I = 1, NRBSITES - 1

                     RISITE = MATMUL(RMI,SITE(I,:)) 

                     DO J = 1, NRBSITES - 1

                        RJSITE = MATMUL(RMJ,SITE(J,:))

                        DSS = DR + RISITE - RJSITE

                        IF (I == 1 .AND. J == 1) THEN
  
                           IF (DOT_PRODUCT(DSS,DSS) < 1.D0) THEN

                           OVRLPT = .TRUE.
                           GO TO 95

                           ENDIF

                        ELSEIF	(I == 2 .AND. J == 2) THEN

                           IF (DOT_PRODUCT(DSS,DSS) < DBSIGBB*DBSIGBB) THEN

                           OVRLPT = .TRUE.
                           GO TO 95

                           ENDIF

                        ELSE

                           IF (DOT_PRODUCT(DSS,DSS) < DBSIGAB*DBSIGAB) THEN

                           OVRLPT = .TRUE.
                           GO TO 95

                           ENDIF

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
      
      END SUBROUTINE TAKESTEPDB
