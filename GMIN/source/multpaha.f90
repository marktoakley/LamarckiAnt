      SUBROUTINE MULTPAHA (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, TPAHA, NCMP, SITE, RBSTLA, STCHRG, RHOCC0, RHOCC10, RHOCC20, &
     &                   RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20,      &
     &                   ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J4, J5, J6, J7, J8, J9, J10, REALNATOMS, OFFSET, FCT(6) 
      INTEGER          :: J1INT, J1FIN, J2INT, J2FIN, TTLST, SUMNC, SUMST
      INTEGER, ALLOCATABLE :: NCRBN(:), NRBST(:), MINDXI(:), MINDXF(:), STINDX(:)
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DVDR, ENERGY1, ENERGY2, ENERGY3
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3), FRIJ(3), TIJ(3), TJI(3) 
      DOUBLE PRECISION, ALLOCATABLE :: R(:,:), E(:,:), CHARGE(:)
      DOUBLE PRECISION, ALLOCATABLE :: DR1(:,:), DR2(:,:), DR3(:,:) 
      DOUBLE PRECISION, ALLOCATABLE :: DE1(:,:), DE2(:,:), DE3(:,:)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DCADR(3), DCBDR(3)
      DOUBLE PRECISION :: RHOCC, RHOHH, RHOCH, COSTA, COSTB, DMPFCT, DDMPDR, EXPFCT 
      DOUBLE PRECISION :: DRIJDPI(3), DRIJDPJ(3), DCADPI(3), DCBDPI(3), DCADPJ(3), DCBDPJ(3)
      DOUBLE PRECISION, PARAMETER :: B = 1.6485D0
      LOGICAL          :: GTEST

!     TPAHA is the total number of different PAH molecules that can be considered; currently TPAHA = 13

      ALLOCATE(MINDXI(TPAHA+1))
      ALLOCATE(MINDXF(TPAHA))
      ALLOCATE(NCRBN(TPAHA))
      ALLOCATE(NRBST(TPAHA))
      ALLOCATE(STINDX(TPAHA))

!     NCMP(i) is the number of i-th PAH molecules present in the mixture; NCMP(i) = 0 if no i-th PAH is present 

!      ASSIGN number of sites for each PAH molecule

      NRBST(1) = 12; NRBST(2) = 18; NRBST(3) = 24; NRBST(4) = 26 
      NRBST(5) = 24; NRBST(6) = 32; NRBST(7) = 34; NRBST(8) = 36
      NRBST(9) = 42; NRBST(10) = 46; NRBST(11) = 60; NRBST(12) = 64
      NRBST(13) = 72
      NCRBN(1) = 6; NCRBN(2) = 10; NCRBN(3) = 14; NCRBN(4) = 16
      NCRBN(5) = 14; NCRBN(6) = 20; NCRBN(7) = 22; NCRBN(8) = 24
      NCRBN(9) = 28; NCRBN(10) = 32; NCRBN(11) = 42; NCRBN(12) = 46
      NCRBN(13) = 54

      SUMNC = 0; SUMST = 0

      DO I = 1, TPAHA

         IF (NCMP (I) == 0) THEN
            MINDXI(I) = SUMNC
            STINDX(I) = SUMST
          ELSE
            MINDXI(I) = SUMNC + 1
            STINDX(I) = SUMST + 1
            SUMNC     = SUMNC + NCMP(I)
            SUMST     = SUMST + NCMP(I)*NRBST(I)
         ENDIF

         MINDXF(I) = SUMNC

      ENDDO

      IF (NCMP(TPAHA) == 0) THEN 
         MINDXI(TPAHA+1) = NATOMS/2
      ELSE
         MINDXI(TPAHA+1) = NATOMS/2 + 1
      ENDIF

!      PRINT *, STINDX

      TTLST = SUMST

!      PRINT *, 'MINDXI'
!      PRINT *, MINDXI
!      PRINT *, 'MINDXF'
!      PRINT *, MINDXF
!      PRINT *, 'TTLST=', TTLST


      ALLOCATE(R(TTLST,3)); ALLOCATE(E(TTLST,3)); ALLOCATE(CHARGE(TTLST))
      ALLOCATE(DR1(TTLST,3)); ALLOCATE(DR2(TTLST,3)); ALLOCATE(DR3(TTLST,3))
      ALLOCATE(DE1(TTLST,3)); ALLOCATE(DE2(TTLST,3)); ALLOCATE(DE3(TTLST,3))

      CALL DEFPAHA()

      FCT(1) = 1; FCT(2) = 2; FCT(3) = 6; FCT(4) = 24; FCT(5) = 120; FCT(6) = 720
      ENERGY = 0.D0; ENERGY1 = 0.D0; ENERGY2 = 0.D0; ENERGY3 = 0.D0

      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      J4    = 0
       
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         J9 = 1    ! J9 is the component label within this DO loop
 
         DO WHILE (MINDXF(J9) < J1)
            J9 = J9 + 1
         ENDDO

!         PRINT *, J1, J9
        
         IF (J9 == 1) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFBENZENE()
         ELSEIF (J9 == 2) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFNAPHTHALENE()
         ELSEIF (J9 == 3) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFANTHRACENE()
         ELSEIF (J9 == 4) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFPYRENE()
         ELSEIF (J9 == 5) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFPHENANTHRENE()
         ELSEIF (J9 == 6) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFPERYLENE()
         ELSEIF (J9 == 7) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFBENZOPERYLENE()
         ELSEIF (J9 == 8) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFCORONENE()
         ELSEIF (J9 == 9) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFBISANTHENE()
         ELSEIF (J9 == 10) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFOVALENE()
         ELSEIF (J9 == 11) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFHEXABENZOCORONENE()
         ELSEIF (J9 == 12) THEN
            NRBSITES = NRBST(J9)
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFOCTABENZOCORONENE()
         ELSEIF (J9 == 13) THEN
            NRBSITES = NRBST(J9) 
            ALLOCATE(SITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            CALL DEFCIRCUMCORONENE()
         ENDIF
 
!         J4 = J4 + NRBST(J9) ! J4 sums up to the total number of sites

         DO J2 = 1, NRBST(J9)

            J4 = J4 + 1

            R(J4,:)    = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
            CHARGE(J4) = STCHRG(J2)
!            PRINT *, J4, R(J4,:)

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

            ENDIF
            
            E(J4,:)  = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST) THEN

               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF

         ENDDO

         IF (ALLOCATED(SITE)) DEALLOCATE(SITE)
         IF (ALLOCATED(RBSTLA)) DEALLOCATE(RBSTLA)
         IF (ALLOCATED(STCHRG)) DEALLOCATE(STCHRG)

      ENDDO

      DO J9 = 1, TPAHA

         IF (NCMP(J9) == 0) CYCLE

         DO J10 = J9, TPAHA 

            IF (NCMP(J10) == 0) CYCLE

            IF (J9 == J10) THEN

               J1INT = MINDXI(J9)
               J1FIN = MINDXF(J9) - 1
!               J2INT = J1INT + 1
               J2FIN = MINDXF(J9)

            ELSE

               J1INT = MINDXI(J9)
               J1FIN = MINDXF(J9)
!               J2INT = MINDXI(J10) 
               J2FIN = MINDXF(J10)

            ENDIF


            DO J1 = J1INT, J1FIN 

               J3 = 3*J1
               J5 = OFFSET + J3

               RI(:)  = X(J3-2:J3)

               IF (J9 == J10) THEN
                  J2INT = J1 + 1
               ELSE
                 J2INT = MINDXI(J10)
               ENDIF
               
!               PRINT *, 'J910', J9, J10
!               PRINT *, J1INT, J1FIN, J2INT, J2FIN

               DO I = 1, NRBST(J9)

                  J7 = STINDX(J9) - 1 + (J1-MINDXI(J9))*NRBST(J9) + I
                  EI(:) = E(J7,:)

                  DO J2 = J2INT, J2FIN

                     J4 = 3*J2
                     J6 = OFFSET + J4

                     DO J = 1, NRBST(J10)

                        J8 = STINDX(J10) - 1 + (J2-MINDXI(J10))*NRBST(J10) + J
!                        PRINT *, J7, J8
                        EJ(:)  = E(J8,:)
                        RSS(:) = R(J7,:) - R(J8,:)
                        R2     = DOT_PRODUCT(RSS(:),RSS(:))
                        ABSRIJ = DSQRT(R2)
                        NR(:)  = RSS(:)/ABSRIJ
                        R2     = 1.D0/R2
                        R6     = R2*R2*R2

!                        PRINT *, R(J7,:)
!                        PRINT *, R(J8,:)

!     CALCULATE THE DAMPING FACTOR

                        DMPFCT = 1.D0
                        DDMPDR = B

                        DO K = 1, 6

                           DMPFCT = DMPFCT + (B*ABSRIJ)**K/FLOAT(FCT(K))
                           IF (K > 1) DDMPDR = DDMPDR + (B**K)*(ABSRIJ)**(K-1)/FLOAT(FCT(K-1))

                        END DO

                        EXPFCT = DEXP(-B*ABSRIJ)
                        DDMPDR = (B*EXPFCT*DMPFCT - EXPFCT*DDMPDR)/ABSRIJ
                        DMPFCT = 1.D0 - EXPFCT*DMPFCT

!     NOW CALCULATE RHOAB

                        COSTA      =-DOT_PRODUCT(NR(:),EI(:))
                        COSTB      = DOT_PRODUCT(NR(:),EJ(:))

                        IF (GTEST) THEN

                           DCADR(:)   =-EI(:)/ABSRIJ - COSTA*R2*RSS(:)
                           DCBDR(:)   = EJ(:)/ABSRIJ - COSTB*R2*RSS(:)

                           DRIJDPI(1) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                           DRIJDPI(2) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                           DRIJDPI(3) = DOT_PRODUCT(RSS(:),DR3(J7,:))

                           DRIJDPJ(1) =-DOT_PRODUCT(RSS(:),DR1(J8,:))
                           DRIJDPJ(2) =-DOT_PRODUCT(RSS(:),DR2(J8,:))
                           DRIJDPJ(3) =-DOT_PRODUCT(RSS(:),DR3(J8,:))

                           DCADPI(1)  =-DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE1(J7,:)) & 
                                       - COSTA*R2*DRIJDPI(1)
                           DCADPI(2)  =-DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE2(J7,:)) &
                                       - COSTA*R2*DRIJDPI(2)
                           DCADPI(3)  =-DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE3(J7,:)) &
                                       - COSTA*R2*DRIJDPI(3)
                           DCBDPI(1)  = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(1)
                           DCBDPI(2)  = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(2)
                           DCBDPI(3)  = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(3)
                
                           DCADPJ(1)  = DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(1)
                           DCADPJ(2)  = DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(2)
                           DCADPJ(3)  = DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(3)

                           DCBDPJ(1)  =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE1(J8,:)) &
                                       - COSTB*R2*DRIJDPJ(1)
                           DCBDPJ(2)  =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE2(J8,:)) &
                                       - COSTB*R2*DRIJDPJ(2)
                           DCBDPJ(3)  =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE3(J8,:)) &
                                       - COSTB*R2*DRIJDPJ(3)

                        ENDIF
    
                        IF (I <= NCRBN(J9) .AND. J <= NCRBN(J10)) THEN

                           RHOCC   = RHOCC0 + RHOCC10*(COSTA + COSTB) + RHOCC20*(1.5D0*COSTA*COSTA & 
                                    + 1.5D0*COSTB*COSTB - 1.D0)
                           EXPFCT  = KKJ*DEXP(-ALPHACC*(ABSRIJ - RHOCC))
                           ENERGY1 = ENERGY1 + EXPFCT                   
                           ENERGY2 = ENERGY2 - DC6CC*DMPFCT*R6

!                           PRINT *, KKJ, ALPHACC, ABSRIJ, RHOCC
!                           STOP

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6CC*R6*R2*DMPFCT - DC6CC*R6*DDMPDR 
                              FRIJ(:) = ALPHACC*EXPFCT*(-NR(:) + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADR(:) &
                                       + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHACC*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCC10  &
                                       + 3.D0*RHOCC20*COSTA)*DCADPI(:) + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHACC*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCC10  &
                                       + 3.D0*RHOCC20*COSTA)*DCADPJ(:) + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPJ(:))

                           ENDIF

                        ELSEIF (I > NCRBN(J9) .AND. J > NCRBN(J10)) THEN

                           RHOHH  = RHOHH0 + RHOHH10*(COSTA + COSTB) + RHOHH20*(1.5D0*COSTA*COSTA      &
                                   + 1.5D0*COSTB*COSTB - 1.D0) 
                           EXPFCT  = KKJ*DEXP(-ALPHAHH*(ABSRIJ - RHOHH))
                           ENERGY1 = ENERGY1 + EXPFCT
                           ENERGY2 = ENERGY2 - DC6HH*DMPFCT*R6

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6HH*R6*R2*DMPFCT - DC6HH*R6*DDMPDR 
                              FRIJ(:) = ALPHAHH*EXPFCT*(-NR(:) + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADR(:) &
                                       + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHAHH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOHH10  &
                                       + 3.D0*RHOHH20*COSTA)*DCADPI(:) + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHAHH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOHH10  &
                                       + 3.D0*RHOHH20*COSTA)*DCADPJ(:) + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPJ(:))

                           ENDIF

                        ELSE IF (I <= NCRBN(J9) .AND. J > NCRBN(J10)) THEN 

                           RHOCH  = RHOCH0 + RHOC10H*COSTA + RHOCH10*COSTB + RHOC20H*(1.5D0*COSTA*COSTA &
                                   - 0.5D0) + RHOCH20*(1.5D0*COSTB*COSTB - 0.5D0)
                           EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                           ENERGY1 = ENERGY1 + EXPFCT
                           ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6

                           IF (GTEST) THEN
                  
                              DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                              FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADR(:) &
                                       + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOC10H  &
                                       + 3.D0*RHOC20H*COSTA)*DCADPI(:) + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOC10H  &
                                       + 3.D0*RHOC20H*COSTA)*DCADPJ(:) + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPJ(:))

                           ENDIF

                        ELSE ! IF(I > NCRBN(J9) .AND. J <= NCRBN(J10)) THEN

                           RHOCH  = RHOCH0 + RHOCH10*COSTA + RHOC10H*COSTB + RHOCH20*(1.5D0*COSTA*COSTA &
                                   - 0.5D0) + RHOC20H*(1.5D0*COSTB*COSTB - 0.5D0)
                           EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                           ENERGY1 = ENERGY1 + EXPFCT
                           ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                              FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADR(:) &
                                       + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCH10  &
                                       + 3.D0*RHOCH20*COSTA)*DCADPI(:) + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCH10  &
                                       + 3.D0*RHOCH20*COSTA)*DCADPJ(:) + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPJ(:))

                           ENDIF

                        ENDIF

                        ENERGY3   = ENERGY3 + CCKJ*CHARGE(J7)*CHARGE(J8)/ABSRIJ
                         
                        IF (GTEST) THEN

                           DVDR   = DVDR - CCKJ*CHARGE(J7)*CHARGE(J8)*R2/ABSRIJ

                           G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:) + FRIJ(:)
                           G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:) - FRIJ(:)

                           G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPI(:) + TIJ(:)
                           G(J6-2:J6) = G(J6-2:J6) + DVDR*DRIJDPJ(:) + TJI(:)

                        ENDIF

                     ENDDO

                  ENDDO
 
               ENDDO

            ENDDO

         ENDDO

      ENDDO

      ENERGY = (ENERGY1 + ENERGY2 + ENERGY3)*2625.499D0 
      IF (GTEST) G(:) = G(:)*2625.499D0

      IF (ALLOCATED(MINDXI)) DEALLOCATE(MINDXI)
      IF (ALLOCATED(MINDXF)) DEALLOCATE(MINDXF)
      IF (ALLOCATED(NCRBN)) DEALLOCATE(NCRBN)
      IF (ALLOCATED(NRBST)) DEALLOCATE(NRBST)
      IF (ALLOCATED(STINDX)) DEALLOCATE(STINDX)

      IF (ALLOCATED(R)) DEALLOCATE(R); IF (ALLOCATED(E)) DEALLOCATE(E); IF (ALLOCATED(CHARGE)) DEALLOCATE(CHARGE)
      IF (ALLOCATED(DR1)) DEALLOCATE(DR1); IF (ALLOCATED(DR2)) DEALLOCATE(DR2); IF (ALLOCATED(DR3)) DEALLOCATE(DR3)
      IF (ALLOCATED(DE1)) DEALLOCATE(DE1); IF (ALLOCATED(DE2)) DEALLOCATE(DE2); IF (ALLOCATED(DE3)) DEALLOCATE(DE3)

      END SUBROUTINE MULTPAHA


!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWMULTPAHA() 

      USE COMMONS, ONLY: NATOMS, NRBSITES, TPAHA, NCMP, SITE, RBSTLA, STCHRG, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7, J9
      INTEGER          :: SUMNC, TTLST
      INTEGER, ALLOCATABLE :: NCRBN(:), NRBST(:), MINDXI(:), MINDXF(:)
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST
 
      OPEN(UNIT=26, FILE='multpaha.xyz', STATUS='UNKNOWN')

      ALLOCATE(MINDXF(TPAHA))
      ALLOCATE(NCRBN(TPAHA))
      ALLOCATE(NRBST(TPAHA))

      NRBST(1) = 12; NRBST(2) = 18; NRBST(3) = 24; NRBST(4) = 26
      NRBST(5) = 24; NRBST(6) = 32; NRBST(7) = 34; NRBST(8) = 36
      NRBST(9) = 42; NRBST(10) = 46; NRBST(11) = 60; NRBST(12) = 64
      NRBST(13) = 72
      NCRBN(1) = 6; NCRBN(2) = 10; NCRBN(3) = 14; NCRBN(4) = 16
      NCRBN(5) = 14; NCRBN(6) = 20; NCRBN(7) = 22; NCRBN(8) = 24
      NCRBN(9) = 28; NCRBN(10) = 32; NCRBN(11) = 42; NCRBN(12) = 46
      NCRBN(13) = 54

      SUMNC = 0; TTLST = 0

      DO I = 1, TPAHA

         SUMNC     = SUMNC + NCMP(I)
         TTLST     = TTLST + NCMP(I)*NRBST(I)
         MINDXF(I) = SUMNC

      ENDDO

      DO J1 = 1, NSAVE

         WRITE(26,'(I6)') TTLST 
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, NATOMS/2

            J5 = 3*J3
            J7 = 3*NATOMS/2 + J5
            P  = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            J9 = 1    ! J9 is the component label within this DO loop

            DO WHILE (MINDXF(J9) < J3)
               J9 = J9 + 1
            ENDDO

            IF (J9 == 1) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFBENZENE()
            ELSEIF (J9 == 2) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFNAPHTHALENE()
            ELSEIF (J9 == 3) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFANTHRACENE()
            ELSEIF (J9 == 4) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFPYRENE()
            ELSEIF (J9 == 5) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFPHENANTHRENE()
            ELSEIF (J9 == 6) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFPERYLENE()
            ELSEIF (J9 == 7) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFBENZOPERYLENE()
            ELSEIF (J9 == 8) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFCORONENE()
            ELSEIF (J9 == 9) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFBISANTHENE()
            ELSEIF (J9 == 10) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFOVALENE()
            ELSEIF (J9 == 11) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFHEXABENZOCORONENE()
            ELSEIF (J9 == 12) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFOCTABENZOCORONENE()
            ELSEIF (J9 == 13) THEN
               NRBSITES = NRBST(J9)
               ALLOCATE(SITE(NRBSITES,3))
               ALLOCATE(RBSTLA(NRBSITES,3))
               ALLOCATE(STCHRG(NRBSITES))
               CALL DEFCIRCUMCORONENE()
            ENDIF

            DO J2 = 1, NRBST(J9)

               RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (J2 <= NCRBN(J9)) THEN
                  WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

            IF (ALLOCATED(SITE)) DEALLOCATE(SITE)
            IF (ALLOCATED(RBSTLA)) DEALLOCATE(RBSTLA)
            IF (ALLOCATED(STCHRG)) DEALLOCATE(STCHRG)

         ENDDO

      ENDDO

      END SUBROUTINE VIEWMULTPAHA
