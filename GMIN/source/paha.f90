      SUBROUTINE PAHA (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, NCARBON, SITE, RBSTLA, STCHRG, RHOCC0, RHOCC10, RHOCC20, &
     &                   RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20,      &
     &                   ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, FCT(6) 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DVDR, ENERGY1, ENERGY2, ENERGY3
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3), FRIJ(3), TIJ(3), TJI(3) 
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3), E(3*NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: DE1(3*NATOMS*NRBSITES/2,3), DE2(3*NATOMS*NRBSITES/2,3), DE3(3*NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DCADR(3), DCBDR(3)
      DOUBLE PRECISION :: RHOCC, RHOHH, RHOCH, COSTA, COSTB, DMPFCT, DDMPDR, EXPFCT 
      DOUBLE PRECISION :: DRIJDPI(3), DRIJDPJ(3), DCADPI(3), DCBDPI(3), DCADPJ(3), DCBDPJ(3)
      DOUBLE PRECISION, PARAMETER :: B = 1.6485D0, EB = 2.0D0
      LOGICAL          :: GTEST

      FCT(1) = 1; FCT(2) = 2; FCT(3) = 6; FCT(4) = 24; FCT(5) = 120; FCT(6) = 720
      ENERGY = 0.D0; ENERGY1 = 0.D0; ENERGY2 = 0.D0; ENERGY3 = 0.D0

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

            J4      = NRBSITES*(J1-1) + J2
            R(J4,:) = RI(:) + MATMUL(RMI(:,:),SITE(J2,:))
            E(J4,:) = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS  - 1

         J3 = 3*J1
         J5 = OFFSET + J3

         RI(:)  = X(J3-2:J3)

         DO I = 1, NRBSITES

            J7    = NRBSITES*(J1-1) + I
            EI(:) = E(J7,:)
               
            DO J2 = J1 + 1, REALNATOMS

               J4 = 3*J2
               J6 = OFFSET + J4

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  EJ(:)  = E(J8,:)
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2
                  R6     = R2*R2*R2

!     CALCULATE THE DISPERSION DAMPING FACTOR

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
    
                  IF (I <= NCARBON .AND. J <= NCARBON) THEN

                     RHOCC   = RHOCC0 + RHOCC10*(COSTA + COSTB) + RHOCC20*(1.5D0*COSTA*COSTA & 
                             + 1.5D0*COSTB*COSTB - 1.D0)
                     EXPFCT  = KKJ*DEXP(-ALPHACC*(ABSRIJ - RHOCC))
                     ENERGY1 = ENERGY1 + EXPFCT                                
                     ENERGY2 = ENERGY2 - DC6CC*DMPFCT*R6

                     IF (GTEST) THEN

                        DVDR    = 6.D0*DC6CC*R6*R2*DMPFCT - DC6CC*R6*DDMPDR 
                        FRIJ(:) = ALPHACC*EXPFCT*(-NR(:) + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADR(:) &
                                + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDR(:))
                        TIJ(:)  = ALPHACC*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPI(:) &
                                + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPI(:))
                        TJI(:)  = ALPHACC*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPJ(:) &
                                + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPJ(:)) 

                     ENDIF

                  ELSEIF (I > NCARBON .AND. J > NCARBON) THEN

                     RHOHH  = RHOHH0 + RHOHH10*(COSTA + COSTB) + RHOHH20*(1.5D0*COSTA*COSTA      &
                            + 1.5D0*COSTB*COSTB - 1.D0) 
                     EXPFCT  = KKJ*DEXP(-ALPHAHH*(ABSRIJ - RHOHH))
                     ENERGY1 = ENERGY1 + EXPFCT                                
                     ENERGY2 = ENERGY2 - DC6HH*DMPFCT*R6

                     IF (GTEST) THEN

                        DVDR    = 6.D0*DC6HH*R6*R2*DMPFCT - DC6HH*R6*DDMPDR 
                        FRIJ(:) = ALPHAHH*EXPFCT*(-NR(:) + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADR(:) &
                                + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDR(:))
                        TIJ(:)  = ALPHAHH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPI(:) &
                                + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPI(:))
                        TJI(:)  = ALPHAHH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPJ(:) &
                                + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPJ(:))

                     ENDIF

                  ELSE IF (I <= NCARBON .AND. J > NCARBON) THEN 

                     RHOCH  = RHOCH0 + RHOC10H*COSTA + RHOCH10*COSTB + RHOC20H*(1.5D0*COSTA*COSTA &
                            - 0.5D0) + RHOCH20*(1.5D0*COSTB*COSTB - 0.5D0)
                     EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                     ENERGY1 = ENERGY1 + EXPFCT                              
                     ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6

                     IF (GTEST) THEN
                  
                        DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                        FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADR(:) &
                                + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDR(:))
                        TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPI(:) &
                                + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPI(:))
                        TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPJ(:) &
                                + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPJ(:))

                     ENDIF

                  ELSE !IF(I > NCARBON .AND. J <= NCARBON) THEN

                     RHOCH  = RHOCH0 + RHOCH10*COSTA + RHOC10H*COSTB + RHOCH20*(1.5D0*COSTA*COSTA &
                            - 0.5D0) + RHOC20H*(1.5D0*COSTB*COSTB - 0.5D0)
                     EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                     ENERGY1 = ENERGY1 + EXPFCT                         
                     ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6

                     IF (GTEST) THEN

                        DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                        FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADR(:) &
                                + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDR(:))
                        TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPI(:) &
                                + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPI(:))
                        TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPJ(:) &
                                + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPJ(:))

                     ENDIF

                  ENDIF

                  ENERGY3   = ENERGY3 + CCKJ*STCHRG(I)*STCHRG(J)/ABSRIJ
 
                  IF (GTEST) THEN


                     DVDR   = DVDR - CCKJ*STCHRG(I)*STCHRG(J)*R2/ABSRIJ


                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:) + FRIJ(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:) - FRIJ(:)

                     G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPI(:) + TIJ(:)
                     G(J6-2:J6) = G(J6-2:J6) + DVDR*DRIJDPJ(:) + TJI(:)

                  ENDIF

               ENDDO

            ENDDO
 
         ENDDO

      ENDDO

      ENERGY = (ENERGY1 + ENERGY2 + ENERGY3)*2625.499D0 
      IF (GTEST) G(:) = G(:)*2625.499D0

      END SUBROUTINE PAHA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPAHA()

      USE COMMONS, ONLY: RHOCC0, RHOCC10, RHOCC20,  RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20, &
                         ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ

      IMPLICIT NONE
 
      ALPHACC = 1.861500D0
      ALPHAHH = 1.431200D0
      ALPHACH = 1.775600D0

      DC6CC    = 30.469D0
      DC6HH    = 5.359D0
      DC6CH    = 12.840D0

      RHOCC0  = 5.814700D0
      RHOCC10 = 0.021700D0
      RHOCC20 =-0.220800D0

      RHOHH0  = 4.486200D0
      RHOHH10 =-0.271800D0
      RHOHH20 = 0.0D0

      RHOCH0  = 5.150500D0
      RHOC10H = 0.021700D0
      RHOCH10 =-0.271800D0
      RHOC20H =-0.220800D0
      RHOCH20 = 0.0D0

      KKJ     = 1.D-03
      CCKJ    = 1.D0   !1389.354848D0

      END SUBROUTINE DEFPAHA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBENZENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG 

      IMPLICIT NONE
 
      INTEGER :: J1

!     C6H6

!     D6h reference geometry: C-C: 1.397 angstrom; C-H: 1.087 angstrom

!      SITE(1,:)  = (/-2.63923430843701,   0.00000000000000,   0.00000000000000/)
!      SITE(2,:)  = (/ 2.63923430843701,   0.00000000000000,   0.00000000000000/)
!      SITE(3,:)  = (/-1.31961715421850,  -2.28564395764590,   0.00000000000000/)
!      SITE(4,:)  = (/ 1.31961715421850,  -2.28564395764590,   0.00000000000000/)
!      SITE(5,:)  = (/-1.31961715421850,   2.28564395764590,   0.00000000000000/)
!      SITE(6,:)  = (/ 1.31961715421850,   2.28564395764590,   0.00000000000000/)
!      SITE(7,:)  = (/-4.69338981379532,   0.00000000000000,   0.00000000000000/)
!      SITE(8,:)  = (/ 4.69338981379532,   0.00000000000000,   0.00000000000000/)
!      SITE(9,:)  = (/ 2.34669490689766,   4.06459480860986,   0.00000000000000/)
!      SITE(10,:) = (/-2.34669490689766,   4.06459480860986,   0.00000000000000/)
!      SITE(11,:) = (/ 2.34669490689766,  -4.06459480860986,   0.00000000000000/)
!      SITE(12,:) = (/-2.34669490689766,  -4.06459480860986,   0.00000000000000/)

      SITE(1,:)  = (/ 2.63923430843701,   0.00000000000000,   0.00000000000000/)
      SITE(2,:)  = (/ 1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      SITE(3,:)  = (/-1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      SITE(4,:)  = (/-2.63923430843701,   0.00000000000000,   0.00000000000000/)
      SITE(5,:)  = (/-1.31961715421850,   2.28564395764590,   0.00000000000000/)
      SITE(6,:)  = (/ 1.31961715421850,   2.28564395764590,   0.00000000000000/)
      SITE(7,:)  = (/ 4.69338981379532,   0.00000000000000,   0.00000000000000/)
      SITE(8,:)  = (/ 2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      SITE(9,:)  = (/-2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      SITE(10,:) = (/-4.69338981379532,   0.00000000000000,   0.00000000000000/)
      SITE(11,:) = (/-2.34669490689766,   4.06459480860986,   0.00000000000000/)
      SITE(12,:) = (/ 2.34669490689766,   4.06459480860986,   0.00000000000000/)

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO H3
      RBSTLA(4,:)  = SITE(10,:) - SITE(4,:)                 ! Z FROM C4 TO H4
      RBSTLA(5,:)  = SITE(11,:) - SITE(5,:)                 ! Z FROM C5 TO H5
      RBSTLA(6,:)  = SITE(12,:) - SITE(6,:)                 ! Z FROM C6 TO H6
      RBSTLA(7,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1!
      RBSTLA(8,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2!
      RBSTLA(9,:)  = SITE(9,:) -  SITE(3,:)                 ! Z FROM C3 TO H3!
      RBSTLA(10,:) = SITE(10,:) - SITE(4,:)                 ! Z FROM C4 TO H4!
      RBSTLA(11,:) = SITE(11,:) - SITE(5,:)                 ! Z FROM C5 TO H5!
      RBSTLA(12,:) = SITE(12,:) - SITE(6,:)                 ! Z FROM C6 TO H6!

!      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1
!      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2
!      RBSTLA(3,:)  = SITE(12,:) - SITE(3,:)                 ! Z FROM C3 TO H6
!      RBSTLA(4,:)  = SITE(11,:) - SITE(4,:)                 ! Z FROM C4 TO H5
!      RBSTLA(5,:)  = SITE(10,:) - SITE(5,:)                 ! Z FROM C5 TO H4
!      RBSTLA(6,:)  = SITE(9,:)  - SITE(6,:)                 ! Z FROM C6 TO H3
!      RBSTLA(7,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1
!      RBSTLA(8,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2
!      RBSTLA(9,:)  = SITE(9,:)  - SITE(6,:)                 ! Z FROM C6 TO H3
!      RBSTLA(10,:) = SITE(10,:) - SITE(5,:)                 ! Z FROM C5 TO H4
!      RBSTLA(11,:) = SITE(11,:) - SITE(4,:)                 ! Z FROM C4 TO H5
!      RBSTLA(12,:) = SITE(12,:) - SITE(3,:)                 ! Z FROM C3 TO H6
      
      DO J1 = 1, NRBSITES
 
         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      STCHRG(1:6)  = -0.11114D0
      STCHRG(7:12) =  0.11114D0

      END SUBROUTINE DEFBENZENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFNAPHTHALENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE
 
      INTEGER :: J1

!     C10H8

      SITE(1,:)  = (/-1.33862D0, -4.59918D0, 0.D0/)   ! C1
      SITE(2,:)  = (/-2.65019D0, -2.35249D0, 0.D0/)   ! C2
      SITE(3,:)  = (/-1.35523D0, 0.D0, 0.D0/)         ! C3
      SITE(4,:)  = (/ 1.35523D0, 0.D0, 0.D0/)         ! C4
      SITE(5,:)  = (/ 2.65019D0,-2.35249D0, 0.D0/)    ! C5
      SITE(6,:)  = (/ 1.33862D0,-4.59918D0, 0.D0/)    ! C6
      SITE(7,:)  = (/-2.65019D0, 2.35249D0, 0.D0/)    ! C7
      SITE(8,:)  = (/ 2.65019D0, 2.35249D0, 0.D0/)    ! C8
      SITE(9,:)  = (/ 1.33862D0, 4.59918D0, 0.D0/)    ! C9
      SITE(10,:) = (/-1.33862D0, 4.59918D0, 0.D0/)    ! C10
      SITE(11,:) = (/-4.70575D0, 2.34799D0, 0.D0/)    ! H1
      SITE(12,:) = (/-2.35493D0,-6.38388D0, 0.D0/)    ! H2
      SITE(13,:) = (/-4.70575D0,-2.34799D0, 0.D0/)    ! H3
      SITE(14,:) = (/ 4.70575D0,-2.34799D0, 0.D0/)    ! H4
      SITE(15,:) = (/ 2.35493D0,-6.38388D0, 0.D0/)    ! H5
      SITE(16,:) = (/ 4.70575D0, 2.34799D0, 0.D0/)    ! H6
      SITE(17,:) = (/ 2.35493D0, 6.38388D0, 0.D0/)    ! H7
      SITE(18,:) = (/-2.35493D0, 6.38388D0, 0.D0/)    ! H8

      STCHRG(1)  =-0.10048D0 
      STCHRG(2)  =-0.29796D0
      STCHRG(3)  = 0.24018D0
      STCHRG(4)  = 0.24018D0
      STCHRG(5)  =-0.29796D0
      STCHRG(6)  =-0.10048D0
      STCHRG(7)  =-0.29796D0
      STCHRG(8)  =-0.29796D0
      STCHRG(9)  =-0.10048D0 
      STCHRG(10) =-0.10048D0
      STCHRG(11) = 0.15530D0
      STCHRG(12) = 0.12304D0
      STCHRG(13) = 0.15530D0
      STCHRG(14) = 0.15530D0
      STCHRG(15) = 0.12304D0
      STCHRG(16) = 0.15530D0
      STCHRG(17) = 0.12304D0
      STCHRG(18) = 0.12304D0

!      SITE(:,:) =  SITE(:,:)*0.529177D0

      RBSTLA(1,:)  = SITE(12,:) - SITE(1,:)                 ! Z FROM C1 TO H2
      RBSTLA(2,:)  = SITE(13,:) - SITE(2,:)                 ! Z FROM C2 TO H3
      RBSTLA(3,:)  = SITE(4,:)  - SITE(3,:)                 ! Z FROM C3 TO C4
      RBSTLA(4,:)  = SITE(3,:)  - SITE(4,:)                 ! Z FROM C4 TO C3
      RBSTLA(5,:)  = SITE(14,:) - SITE(5,:)                 ! Z FROM C5 TO H4
      RBSTLA(6,:)  = SITE(15,:) - SITE(6,:)                 ! Z FROM C6 TO H5
      RBSTLA(7,:)  = SITE(11,:) - SITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(8,:)  = SITE(16,:) - SITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(9,:)  = SITE(17,:) - SITE(9,:)                 ! Z FROM C9 TO H7
      RBSTLA(10,:) = SITE(18,:) - SITE(10,:)                ! Z FROM C10 TO H8
      RBSTLA(11,:) = SITE(11,:) - SITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(12,:) = SITE(12,:) - SITE(1,:)                 ! Z FROM C1 TO H2
      RBSTLA(13,:) = SITE(13,:) - SITE(2,:)                 ! Z FROM C2 TO H3
      RBSTLA(14,:) = SITE(14,:) - SITE(5,:)                 ! Z FROM C5 TO H4
      RBSTLA(15,:) = SITE(15,:) - SITE(6,:)                 ! Z FROM C6 TO H5
      RBSTLA(16,:) = SITE(16,:) - SITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(17,:) = SITE(17,:) - SITE(9,:)                 ! Z FROM C9 TO H7
      RBSTLA(18,:) = SITE(18,:) - SITE(10,:)                ! Z FROM C10 TO H8
      
      DO J1 = 1, NRBSITES
 
         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFNAPHTHALENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFANTHRACENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C14H10

      SITE(1,:)  = (/ 1.36540D0, 2.31298D0, 0.D0/)    ! C1
      SITE(2,:)  = (/-1.36540D0, 2.31298D0, 0.D0/)    ! C2     
      SITE(3,:)  = (/ 1.36540D0,-2.31298D0, 0.D0/)    ! C3
      SITE(4,:)  = (/-1.36540D0,-2.31298D0, 0.D0/)    ! C4
      SITE(5,:)  = (/-2.65253D0, 0.D0, 0.D0/)         ! C5
      SITE(6,:)  = (/ 2.65253D0, 0.D0, 0.D0/)         ! C6
      SITE(7,:)  = (/ 2.65927D0, 4.68538D0, 0.D0/)    ! C7
      SITE(8,:)  = (/-2.65927D0, 4.68538D0, 0.D0/)    ! C8
      SITE(9,:)  = (/ 2.65927D0,-4.68538D0, 0.D0/)    ! C9
      SITE(10,:) = (/-2.65927D0,-4.68538D0, 0.D0/)    ! C10
      SITE(11,:) = (/ 1.34762D0,-6.91760D0, 0.D0/)    ! C11
      SITE(12,:) = (/-1.34762D0,-6.91760D0, 0.D0/)    ! C12
      SITE(13,:) = (/ 1.34762D0, 6.91760D0, 0.D0/)    ! C13
      SITE(14,:) = (/-1.34762D0, 6.91760D0, 0.D0/)    ! C14
      SITE(15,:) = (/ 4.71450D0,-4.67888D0, 0.D0/)    ! H1
      SITE(16,:) = (/ 2.35428D0,-8.70751D0, 0.D0/)    ! H2
      SITE(17,:) = (/-2.35428D0,-8.70751D0, 0.D0/)    ! H3
      SITE(18,:) = (/-4.71450D0,-4.67888D0, 0.D0/)    ! H4
      SITE(19,:) = (/ 4.71450D0, 4.67888D0, 0.D0/)    ! H5
      SITE(20,:) = (/ 2.35428D0, 8.70751D0, 0.D0/)    ! H6
      SITE(21,:) = (/-2.35428D0, 8.70751D0, 0.D0/)    ! H7
      SITE(22,:) = (/-4.71450D0, 4.67888D0, 0.D0/)    ! H8
      SITE(23,:) = (/-4.70918D0, 0.D0, 0.D0/)         ! H9
      SITE(24,:) = (/ 4.70918D0, 0.D0, 0.D0/)         ! H10

      STCHRG(1)  = 0.23448D0 
      STCHRG(2)  = 0.23448D0
      STCHRG(3)  = 0.23448D0
      STCHRG(4)  = 0.23448D0
      STCHRG(5)  =-0.47174D0
      STCHRG(6)  =-0.47174D0
      STCHRG(7)  =-0.25252D0
      STCHRG(8)  =-0.25252D0
      STCHRG(9)  =-0.25252D0
      STCHRG(10) =-0.25252D0
      STCHRG(11) =-0.11389D0
      STCHRG(12) =-0.11389D0
      STCHRG(13) =-0.11389D0
      STCHRG(14) =-0.11389D0
      STCHRG(15) = 0.14291D0
      STCHRG(16) = 0.12531D0
      STCHRG(17) = 0.12531D0
      STCHRG(18) = 0.14291D0
      STCHRG(19) = 0.14291D0
      STCHRG(20) = 0.12531D0
      STCHRG(21) = 0.12531D0
      STCHRG(22) = 0.14291D0
      STCHRG(23) = 0.19915D0
      STCHRG(24) = 0.19915D0

      RBSTLA(1,:)  = SITE(2,:)  - SITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = SITE(1,:)  - SITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = SITE(4,:)  - SITE(3,:)                 ! Z FROM C3 TO C4
      RBSTLA(4,:)  = SITE(3,:)  - SITE(4,:)                 ! Z FROM C4 TO C3
      RBSTLA(5,:)  = SITE(23,:) - SITE(5,:)                 ! Z FROM C5 TO H9
      RBSTLA(6,:)  = SITE(24,:) - SITE(6,:)                 ! Z FROM C6 TO H10
      RBSTLA(7,:)  = SITE(19,:) - SITE(7,:)                 ! Z FROM C7 TO H5
      RBSTLA(8,:)  = SITE(22,:) - SITE(8,:)                 ! Z FROM C8 TO H8
      RBSTLA(9,:)  = SITE(15,:) - SITE(9,:)                 ! Z FROM C9 TO H1
      RBSTLA(10,:) = SITE(18,:) - SITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(11,:) = SITE(16,:) - SITE(11,:)                ! Z FROM C11 TO H2
      RBSTLA(12,:) = SITE(17,:) - SITE(12,:)                ! Z FROM C12 TO H3
      RBSTLA(13,:) = SITE(20,:) - SITE(13,:)                ! Z FROM C13 TO H6
      RBSTLA(14,:) = SITE(21,:) - SITE(14,:)                ! Z FROM C14 TO H7
      RBSTLA(15,:) = SITE(15,:) - SITE(9,:)                 ! Z FROM C9 TO H1
      RBSTLA(16,:) = SITE(16,:) - SITE(11,:)                ! Z FROM C11 TO H2
      RBSTLA(17,:) = SITE(17,:) - SITE(12,:)                ! Z FROM C12 TO H3
      RBSTLA(18,:) = SITE(18,:) - SITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(19,:) = SITE(19,:) - SITE(7,:)                 ! Z FROM C7 TO H5
      RBSTLA(20,:) = SITE(20,:) - SITE(13,:)                ! Z FROM C13 TO H6
      RBSTLA(21,:) = SITE(21,:) - SITE(14,:)                ! Z FROM C14 TO H7
      RBSTLA(22,:) = SITE(22,:) - SITE(8,:)                 ! Z FROM C8 TO H8
      RBSTLA(23,:) = SITE(23,:) - SITE(5,:)                 ! Z FROM C5 TO H9
      RBSTLA(24,:) = SITE(24,:) - SITE(6,:)                 ! Z FROM C6 TO H10

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFANTHRACENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPYRENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C16H10

      SITE(1,:)  = (/-1.34794D0, 0.D0, 0.D0/)         ! C1
      SITE(2,:)  = (/ 1.34794D0, 0.D0, 0.D0/)         ! C2
      SITE(3,:)  = (/ 2.70059D0, 2.33625D0, 0.D0/)    ! C3
      SITE(4,:)  = (/ 2.70059d0,-2.33625D0, 0.D0/)    ! C4
      SITE(5,:)  = (/-2.70059D0,-2.33625D0, 0.D0/)    ! C5
      SITE(6,:)  = (/-2.70059D0, 2.33625D0, 0.D0/)    ! C6
      SITE(7,:)  = (/ 1.28651D0, 4.65603D0, 0.D0/)    ! C7
      SITE(8,:)  = (/ 5.35355D0, 2.28771D0, 0.D0/)    ! C8
      SITE(9,:)  = (/ 1.28651D0,-4.65603D0, 0.D0/)    ! C9
      SITE(10,:) = (/ 5.35355D0,-2.28771D0, 0.D0/)    ! C10
      SITE(11,:) = (/-1.28651D0,-4.65603D0, 0.D0/)    ! C11
      SITE(12,:) = (/-5.35355D0,-2.28771D0, 0.D0/)    ! C12
      SITE(13,:) = (/-1.28651D0, 4.65603D0, 0.D0/)    ! C13
      SITE(14,:) = (/-5.35355D0, 2.28771D0, 0.D0/)    ! C14
      SITE(15,:) = (/ 6.65929D0, 0.D0, 0.D0/)         ! C15
      SITE(16,:) = (/-6.65929D0, 0.D0, 0.D0/)         ! C16
      SITE(17,:) = (/ 2.32543D0, 6.42907D0, 0.D0/)    ! H1
      SITE(18,:) = (/ 6.38694D0, 4.06382D0, 0.D0/)    ! H2
      SITE(19,:) = (/ 2.32543D0,-6.42907D0, 0.D0/)    ! H3
      SITE(20,:) = (/ 6.38694D0,-4.06382D0, 0.D0/)    ! H4
      SITE(21,:) = (/-2.32543D0,-6.42907D0, 0.D0/)    ! H5
      SITE(22,:) = (/-6.38694D0,-4.06382D0, 0.D0/)    ! H6
      SITE(23,:) = (/-2.32543D0, 6.42907D0, 0.D0/)    ! H7
      SITE(24,:) = (/-6.38694D0, 4.06382D0, 0.D0/)    ! H8
      SITE(25,:) = (/ 8.71284D0, 0.D0, 0.D0/)         ! H9
      SITE(26,:) = (/-8.71284D0, 0.D0, 0.D0/)         ! H10

      STCHRG(1)  =-0.04275D0 
      STCHRG(2)  =-0.04275D0
      STCHRG(3)  = 0.22339D0
      STCHRG(4)  = 0.22339D0
      STCHRG(5)  = 0.22339D0
      STCHRG(6)  = 0.22339D0
      STCHRG(7)  =-0.24782D0
      STCHRG(8)  =-0.29542D0
      STCHRG(9)  =-0.24782D0
      STCHRG(10) =-0.29542D0
      STCHRG(11) =-0.24782D0
      STCHRG(12) =-0.29542D0
      STCHRG(13) =-0.24782D0
      STCHRG(14) =-0.29542D0
      STCHRG(15) =-0.05466D0
      STCHRG(16) =-0.05466D0
      STCHRG(17) = 0.15533D0
      STCHRG(18) = 0.15109D0
      STCHRG(19) = 0.15533D0
      STCHRG(20) = 0.15109D0
      STCHRG(21) = 0.15533D0
      STCHRG(22) = 0.15109D0
      STCHRG(23) = 0.15533D0
      STCHRG(24) = 0.15109D0
      STCHRG(25) = 0.12425D0
      STCHRG(26) = 0.12425D0

      RBSTLA(1,:)  = SITE(2,:)  - SITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = SITE(1,:)  - SITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = SITE(2,:)  - SITE(3,:)                 ! Z FROM C3 TO C2
      RBSTLA(4,:)  = SITE(2,:)  - SITE(4,:)                 ! Z FROM C4 TO C2
      RBSTLA(5,:)  = SITE(1,:)  - SITE(5,:)                 ! Z FROM C5 TO C1
      RBSTLA(6,:)  = SITE(1,:)  - SITE(6,:)                 ! Z FROM C6 TO C1
      RBSTLA(7,:)  = SITE(17,:) - SITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(8,:)  = SITE(18,:) - SITE(8,:)                 ! Z FROM C8 TO H2
      RBSTLA(9,:)  = SITE(19,:) - SITE(9,:)                 ! Z FROM C9 TO H3
      RBSTLA(10,:) = SITE(20,:) - SITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(11,:) = SITE(21,:) - SITE(11,:)                ! Z FROM C11 TO H5
      RBSTLA(12,:) = SITE(22,:) - SITE(12,:)                ! Z FROM C12 TO H6
      RBSTLA(13,:) = SITE(23,:) - SITE(13,:)                ! Z FROM C13 TO H7
      RBSTLA(14,:) = SITE(24,:) - SITE(14,:)                ! Z FROM C14 TO H8
      RBSTLA(15,:) = SITE(25,:) - SITE(15,:)                ! Z FROM C15 TO H9
      RBSTLA(16,:) = SITE(26,:) - SITE(16,:)                ! Z FROM C16 TO H10
      RBSTLA(17,:) = SITE(17,:) - SITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(18,:) = SITE(18,:) - SITE(8,:)                 ! Z FROM C8 TO H2
      RBSTLA(19,:) = SITE(19,:) - SITE(9,:)                 ! Z FROM C9 TO H3
      RBSTLA(20,:) = SITE(20,:) - SITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(21,:) = SITE(21,:) - SITE(11,:)                ! Z FROM C11 TO H5
      RBSTLA(22,:) = SITE(22,:) - SITE(12,:)                ! Z FROM C12 TO H6
      RBSTLA(23,:) = SITE(23,:) - SITE(13,:)                ! Z FROM C13 TO H7
      RBSTLA(24,:) = SITE(24,:) - SITE(14,:)                ! Z FROM C14 TO H8
      RBSTLA(25,:) = SITE(25,:) - SITE(15,:)                ! Z FROM C15 TO H9
      RBSTLA(26,:) = SITE(26,:) - SITE(16,:)                ! Z FROM C16 TO H10

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFPYRENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPHENANTHRENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C14H10

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.728950D0,  1.277200D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/-0.728950D0,  1.277200D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/0.679800D0,  -1.197690D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/-0.679800D0,  -1.197690D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-1.423030D0,  0.030020D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/1.423030D0,  0.030020D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/1.500740D0,  2.462850D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/-1.500740D0,  2.462850D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/2.883480D0,  2.425000D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-2.883480D0,  2.425000D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/2.837310D0,  0.016730D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.837310D0,  0.016730D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/3.562000D0,  1.192060D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/-3.562000D0,  1.192060D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/1.232290D0,  -2.134300D0,  0.000000D0/)   !H1
      SITE(16,:)  = (/-1.232290D0,  -2.134300D0,  0.000000D0/)   !H2
      SITE(17,:)  = (/1.007910D0,  3.429160D0,  0.000000D0/)   !H3
      SITE(18,:)  = (/3.447220D0,  3.354020D0,  0.000000D0/)   !H4
      SITE(19,:)  = (/3.347890D0,  -0.943600D0,  0.000000D0/)   !H5
      SITE(20,:)  = (/-1.007910D0,  3.429160D0,  0.000000D0/)   !H6
      SITE(21,:)  = (/-3.447220D0,  3.354020D0,  0.000000D0/)   !H7
      SITE(22,:)  = (/-3.347890D0,  -0.943600D0,  0.000000D0/)   !H8
      SITE(23,:)  = (/4.648270D0,  1.167030D0,  0.000000D0/)   !H9
      SITE(24,:)  = (/-4.648270D0,  1.167030D0,  0.000000D0/)   !H10

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = 0.0153630D0
      STCHRG(2)  = 0.0153630D0
      STCHRG(3)  = -0.2774610D0
      STCHRG(4)  = -0.2774610D0
      STCHRG(5)  = 0.2396090D0
      STCHRG(6)  = 0.2396090D0
      STCHRG(7)  = -0.1930180D0
      STCHRG(8)  = -0.1930180D0
      STCHRG(9)  = -0.1423640D0
      STCHRG(10)  = -0.1423640D0
      STCHRG(11)  = -0.2681230D0
      STCHRG(12)  = -0.2681230D0
      STCHRG(13)  = -0.086620D0
      STCHRG(14)  = -0.086620D0
      STCHRG(15)  = 0.1696350D0
      STCHRG(16)  = 0.1696350D0
      STCHRG(17)  = 0.1409570D0
      STCHRG(18)  = 0.1312140D0
      STCHRG(19)  = 0.1446370D0
      STCHRG(20)  = 0.1409570D0
      STCHRG(21)  = 0.1312140D0
      STCHRG(22)  = 0.1446370D0
      STCHRG(23)  = 0.1261710D0
      STCHRG(24)  = 0.1261710D0

      RBSTLA(1,:)  = SITE(2,:)  - SITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = SITE(1,:)  - SITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = SITE(15,:)  - SITE(3,:)                 ! Z FROM C3 TO H1
      RBSTLA(4,:)  = SITE(16,:)  - SITE(4,:)                 ! Z FROM C4 TO H2
      RBSTLA(5,:)  = SITE(2,:)  - SITE(5,:)                 ! Z FROM C5 TO C2
      RBSTLA(6,:)  = SITE(1,:)  - SITE(6,:)                 ! Z FROM C6 TO C1
      RBSTLA(7,:)  = SITE(17,:)  - SITE(7,:)                 ! Z FROM C7 TO H3
      RBSTLA(8,:)  = SITE(20,:)  - SITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(9,:)  = SITE(18,:)  - SITE(9,:)                 ! Z FROM C9 TO H4
      RBSTLA(10,:)  = SITE(21,:)  - SITE(10,:)                 ! Z FROM C10 TO H7
      RBSTLA(11,:)  = SITE(19,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(12,:)  = SITE(22,:)  - SITE(12,:)                 ! Z FROM C12 TO H8
      RBSTLA(13,:)  = SITE(23,:)  - SITE(13,:)                 ! Z FROM C13 TO H9
      RBSTLA(14,:)  = SITE(24,:)  - SITE(14,:)                 ! Z FROM C14 TO H10
      RBSTLA(15,:)  = SITE(15,:)  - SITE(3,:)                 ! Z FROM C3 TO H1
      RBSTLA(16,:)  = SITE(16,:)  - SITE(4,:)                 ! Z FROM C4 TO H2
      RBSTLA(17,:)  = SITE(17,:)  - SITE(7,:)                 ! Z FROM C7 TO H3
      RBSTLA(18,:)  = SITE(18,:)  - SITE(9,:)                 ! Z FROM C9 TO H4
      RBSTLA(19,:)  = SITE(19,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(20,:)  = SITE(20,:)  - SITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(21,:)  = SITE(21,:)  - SITE(10,:)                 ! Z FROM C10 TO H7
      RBSTLA(22,:)  = SITE(22,:)  - SITE(12,:)                 ! Z FROM C12 TO H8
      RBSTLA(23,:)  = SITE(23,:)  - SITE(13,:)                 ! Z FROM C13 TO H9
      RBSTLA(24,:)  = SITE(24,:)  - SITE(14,:)                 ! Z FROM C14 TO H10

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFPHENANTHRENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPERYLENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C20H12

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.000000D0,  -1.439410D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/0.000000D0,  1.439410D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/1.249970D0,  -0.738310D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/-1.249970D0,  0.738310D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/1.249970D0,  0.738310D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-1.249970D0,  -0.738310D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/2.427620D0,  -1.479640D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/-2.427620D0,  1.479640D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/2.427620D0,  1.479640D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-2.427620D0,  -1.479640D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/2.422750D0,  -2.886150D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.422750D0,  2.886150D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/2.422750D0,  2.886150D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/-2.422750D0,  -2.886150D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/1.232560D0,  -3.575640D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/-1.232560D0,  3.575640D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/1.232560D0,  3.575640D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/-1.232560D0,  -3.575640D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/0.000000D0,  -2.874710D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/0.000000D0,  2.874710D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/1.217830D0,  -4.662620D0,  0.000000D0/)   !H1
      SITE(22,:)  = (/-1.217830D0,  4.662620D0,  0.000000D0/)   !H2
      SITE(23,:)  = (/1.217830D0,  4.662620D0,  0.000000D0/)   !H3
      SITE(24,:)  = (/-1.217830D0,  -4.662620D0,  0.000000D0/)   !H4
      SITE(25,:)  = (/3.368390D0,  -3.421400D0,  0.000000D0/)   !H5
      SITE(26,:)  = (/-3.368390D0,  3.421400D0,  0.000000D0/)   !H6
      SITE(27,:)  = (/3.368390D0,  3.421400D0,  0.000000D0/)   !H7
      SITE(28,:)  = (/-3.368390D0,  -3.421400D0,  0.000000D0/)   !H8
      SITE(29,:)  = (/3.388320D0,  -0.977500D0,  0.000000D0/)   !H9
      SITE(30,:)  = (/-3.388320D0,  0.977500D0,  0.000000D0/)   !H10
      SITE(31,:)  = (/3.388320D0,  0.977500D0,  0.000000D0/)   !H11
      SITE(32,:)  = (/-3.388320D0,  -0.977500D0,  0.000000D0/)   !H12

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = 0.0235520D0
      STCHRG(2)  = 0.0235520D0
      STCHRG(3)  = 0.0294130D0
      STCHRG(4)  = 0.0294130D0
      STCHRG(5)  = 0.0294130D0
      STCHRG(6)  = 0.0294130D0
      STCHRG(7)  = -0.1999210D0
      STCHRG(8)  = -0.1999210D0
      STCHRG(9)  = -0.1999210D0
      STCHRG(10)  = -0.1999210D0
      STCHRG(11)  = -0.0757740D0
      STCHRG(12)  = -0.0757740D0
      STCHRG(13)  = -0.0757740D0
      STCHRG(14)  = -0.0757740D0
      STCHRG(15)  = -0.3459910D0
      STCHRG(16)  = -0.3459910D0
      STCHRG(17)  = -0.3459910D0
      STCHRG(18)  = -0.3459910D0
      STCHRG(19)  = 0.2862460D0
      STCHRG(20)  = 0.2862460D0
      STCHRG(21)  = 0.1680010D0
      STCHRG(22)  = 0.1680010D0
      STCHRG(23)  = 0.1680010D0
      STCHRG(24)  = 0.1680010D0
      STCHRG(25)  = 0.1286140D0
      STCHRG(26)  = 0.1286140D0
      STCHRG(27)  = 0.1286140D0
      STCHRG(28)  = 0.1286140D0
      STCHRG(29)  = 0.1407590D0
      STCHRG(30)  = 0.1407590D0
      STCHRG(31)  = 0.1407590D0
      STCHRG(32)  = 0.1407590D0

      RBSTLA(1,:)  = SITE(19,:)  - SITE(1,:)                 ! Z FROM C1 TO C19
      RBSTLA(2,:)  = SITE(20,:)  - SITE(2,:)                 ! Z FROM C2 TO C20
      RBSTLA(3,:)  = SITE(1,:)  - SITE(3,:)                 ! Z FROM C3 TO C1
      RBSTLA(4,:)  = SITE(2,:)  - SITE(4,:)                 ! Z FROM C4 TO C2
      RBSTLA(5,:)  = SITE(2,:)  - SITE(5,:)                 ! Z FROM C5 TO C2
      RBSTLA(6,:)  = SITE(1,:)  - SITE(6,:)                 ! Z FROM C6 TO C1
      RBSTLA(7,:)  = SITE(29,:)  - SITE(7,:)                 ! Z FROM C7 TO H9
      RBSTLA(8,:)  = SITE(30,:)  - SITE(8,:)                 ! Z FROM C8 TO H10
      RBSTLA(9,:)  = SITE(31,:)  - SITE(9,:)                 ! Z FROM C9 TO H11
      RBSTLA(10,:)  = SITE(32,:)  - SITE(10,:)                 ! Z FROM C10 TO H12
      RBSTLA(11,:)  = SITE(25,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(12,:)  = SITE(26,:)  - SITE(12,:)                 ! Z FROM C12 TO H6
      RBSTLA(13,:)  = SITE(27,:)  - SITE(13,:)                 ! Z FROM C13 TO H7
      RBSTLA(14,:)  = SITE(28,:)  - SITE(14,:)                 ! Z FROM C14 TO H8
      RBSTLA(15,:)  = SITE(21,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(16,:)  = SITE(22,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(17,:)  = SITE(23,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(18,:)  = SITE(24,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(19,:)  = SITE(21,:)  - SITE(19,:)                 ! Z FROM C19 TO H1
      RBSTLA(20,:)  = SITE(22,:)  - SITE(20,:)                 ! Z FROM C20 TO H2
      RBSTLA(21,:)  = SITE(21,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(22,:)  = SITE(22,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(23,:)  = SITE(23,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(24,:)  = SITE(24,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(25,:)  = SITE(25,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(26,:)  = SITE(26,:)  - SITE(12,:)                 ! Z FROM C12 TO H6
      RBSTLA(27,:)  = SITE(27,:)  - SITE(13,:)                 ! Z FROM C13 TO H7
      RBSTLA(28,:)  = SITE(28,:)  - SITE(14,:)                 ! Z FROM C14 TO H8
      RBSTLA(29,:)  = SITE(29,:)  - SITE(7,:)                 ! Z FROM C7 TO H9
      RBSTLA(30,:)  = SITE(30,:)  - SITE(8,:)                 ! Z FROM C8 TO H10
      RBSTLA(31,:)  = SITE(31,:)  - SITE(9,:)                 ! Z FROM C9 TO H11
      RBSTLA(32,:)  = SITE(32,:)  - SITE(10,:)                 ! Z FROM C10 TO H12

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFPERYLENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBENZOPERYLENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C22H12

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/-1.431340D0,  -0.344590D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/1.431340D0,  -0.344590D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/-0.713810D0,  0.895350D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/0.734840D0,  -1.591900D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/0.713810D0,  0.895350D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-0.734840D0,  -1.591900D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/-1.418430D0,  2.128550D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/1.490630D0,  -2.770590D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/1.418430D0,  2.128550D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-1.490630D0,  -2.770590D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-2.850180D0,  2.113800D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/2.888170D0,  -2.744870D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/2.850180D0,  2.113800D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/-2.888170D0,  -2.744870D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/-3.540760D0,  0.939410D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/3.570020D0,  -1.538170D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/3.540760D0,  0.939410D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/-3.570020D0,  -1.538170D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-0.687470D0,  3.343090D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/0.687470D0,  3.343090D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/-2.860950D0,  -0.323210D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/2.860950D0,  -0.323210D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/-4.628030D0,  0.938640D0,  0.000000D0/)   !H1
      SITE(24,:)  = (/4.657020D0,  -1.517650D0,  0.000000D0/)   !H2
      SITE(25,:)  = (/4.628030D0,  0.938640D0,  0.000000D0/)   !H3
      SITE(26,:)  = (/-4.657020D0,  -1.517650D0,  0.000000D0/)   !H4
      SITE(27,:)  = (/-3.378060D0,  3.064450D0,  0.000000D0/)   !H5
      SITE(28,:)  = (/3.438960D0,  -3.681560D0,  0.000000D0/)   !H6
      SITE(29,:)  = (/3.378060D0,  3.064450D0,  0.000000D0/)   !H7
      SITE(30,:)  = (/-3.438960D0,  -3.681560D0,  0.000000D0/)   !H8
      SITE(31,:)  = (/0.995680D0,  -3.735330D0,  0.000000D0/)   !H9
      SITE(32,:)  = (/-0.995680D0,  -3.735330D0,  0.000000D0/)   !H10
      SITE(33,:)  = (/-1.234990D0,  4.282540D0,  0.000000D0/)   !H11
      SITE(34,:)  = (/1.234990D0,  4.282540D0,  0.000000D0/)   !H12

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0032350D0
      STCHRG(2)  = -0.0032350D0
      STCHRG(3)  = -0.0485040D0
      STCHRG(4)  = 0.0531290D0
      STCHRG(5)  = -0.0485040D0
      STCHRG(6)  = 0.0531290D0
      STCHRG(7)  = 0.1924370D0
      STCHRG(8)  = -0.2123090D0
      STCHRG(9)  = 0.1924370D0
      STCHRG(10)  = -0.2123090D0
      STCHRG(11)  = -0.2288940D0
      STCHRG(12)  = -0.081160D0
      STCHRG(13)  = -0.2288940D0
      STCHRG(14)  = -0.081160D0
      STCHRG(15)  = -0.2808820D0
      STCHRG(16)  = -0.2932090D0
      STCHRG(17)  = -0.2808820D0
      STCHRG(18)  = -0.2932090D0
      STCHRG(19)  = -0.2320860D0
      STCHRG(20)  = -0.2320860D0
      STCHRG(21)  = 0.2334920D0
      STCHRG(22)  = 0.2334920D0
      STCHRG(23)  = 0.1662180D0
      STCHRG(24)  = 0.1553260D0
      STCHRG(25)  = 0.1662180D0
      STCHRG(26)  = 0.1553260D0
      STCHRG(27)  = 0.1545260D0
      STCHRG(28)  = 0.1290550D0
      STCHRG(29)  = 0.1545260D0
      STCHRG(30)  = 0.1290550D0
      STCHRG(31)  = 0.1437480D0
      STCHRG(32)  = 0.1437480D0
      STCHRG(33)  = 0.1523490D0
      STCHRG(34)  = 0.1523490D0

      RBSTLA(1,:)  = SITE(3,:)  - SITE(1,:)                 ! Z FROM C1 TO C3
      RBSTLA(2,:)  = SITE(5,:)  - SITE(2,:)                 ! Z FROM C2 TO C5
      RBSTLA(3,:)  = SITE(7,:)  - SITE(3,:)                 ! Z FROM C3 TO C7
      RBSTLA(4,:)  = SITE(2,:)  - SITE(4,:)                 ! Z FROM C4 TO C2
      RBSTLA(5,:)  = SITE(9,:)  - SITE(5,:)                 ! Z FROM C5 TO C9
      RBSTLA(6,:)  = SITE(1,:)  - SITE(6,:)                 ! Z FROM C6 TO C1
      RBSTLA(7,:)  = SITE(3,:)  - SITE(7,:)                 ! Z FROM C7 TO C3
      RBSTLA(8,:)  = SITE(31,:)  - SITE(8,:)                 ! Z FROM C8 TO H9
      RBSTLA(9,:)  = SITE(5,:)  - SITE(9,:)                 ! Z FROM C9 TO C5
      RBSTLA(10,:)  = SITE(32,:)  - SITE(10,:)                 ! Z FROM C10 TO H10
      RBSTLA(11,:)  = SITE(27,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(12,:)  = SITE(28,:)  - SITE(12,:)                 ! Z FROM C12 TO H6
      RBSTLA(13,:)  = SITE(29,:)  - SITE(13,:)                 ! Z FROM C13 TO H7
      RBSTLA(14,:)  = SITE(30,:)  - SITE(14,:)                 ! Z FROM C14 TO H8
      RBSTLA(15,:)  = SITE(23,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(16,:)  = SITE(24,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(17,:)  = SITE(25,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(18,:)  = SITE(26,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(19,:)  = SITE(33,:)  - SITE(19,:)                 ! Z FROM C19 TO H11
      RBSTLA(20,:)  = SITE(34,:)  - SITE(20,:)                 ! Z FROM C20 TO H12
      RBSTLA(21,:)  = SITE(1,:)  - SITE(21,:)                 ! Z FROM C21 TO C1
      RBSTLA(22,:)  = SITE(2,:)  - SITE(22,:)                 ! Z FROM C22 TO C2
      RBSTLA(23,:)  = SITE(23,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(24,:)  = SITE(24,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(25,:)  = SITE(25,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(26,:)  = SITE(26,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(27,:)  = SITE(27,:)  - SITE(11,:)                 ! Z FROM C11 TO H5
      RBSTLA(28,:)  = SITE(28,:)  - SITE(12,:)                 ! Z FROM C12 TO H6
      RBSTLA(29,:)  = SITE(29,:)  - SITE(13,:)                 ! Z FROM C13 TO H7
      RBSTLA(30,:)  = SITE(30,:)  - SITE(14,:)                 ! Z FROM C14 TO H8
      RBSTLA(31,:)  = SITE(31,:)  - SITE(8,:)                 ! Z FROM C8 TO H9
      RBSTLA(32,:)  = SITE(32,:)  - SITE(10,:)                 ! Z FROM C10 TO H10
      RBSTLA(33,:)  = SITE(33,:)  - SITE(19,:)                 ! Z FROM C19 TO H11
      RBSTLA(34,:)  = SITE(34,:)  - SITE(20,:)                 ! Z FROM C20 TO H12

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFBENZOPERYLENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCORONENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C24H12

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/1.427510D0,  0.000000D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/0.713760D0,  1.236260D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/-0.713750D0,  1.236260D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/-1.427500D0,  0.000000D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-0.713750D0,  -1.236260D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/0.713760D0,  -1.236260D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/2.849060D0,  0.000000D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/1.424540D0,  2.467360D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/-1.424530D0,  2.467360D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-2.849050D0,  0.000000D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-1.424530D0,  -2.467360D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/1.424540D0,  -2.467360D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/3.534380D0,  1.248270D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/2.848230D0,  2.436730D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/0.686160D0,  3.685010D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/-0.686150D0,  3.685010D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/-2.848230D0,  2.436730D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/-3.534380D0,  1.248270D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-3.534380D0,  -1.248270D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-2.848230D0,  -2.436730D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/-0.686150D0,  -3.685010D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/0.686160D0,  -3.685010D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/2.848230D0,  -2.436730D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/3.534380D0,  -1.248270D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/4.621790D0,  1.245370D0,  0.000000D0/)   !H1
      SITE(26,:)  = (/3.389430D0,  3.379890D0,  0.000000D0/)   !H2
      SITE(27,:)  = (/1.232370D0,  4.625280D0,  0.000000D0/)   !H3
      SITE(28,:)  = (/-1.232360D0,  4.625280D0,  0.000000D0/)   !H4
      SITE(29,:)  = (/-3.389420D0,  3.379890D0,  0.000000D0/)   !H5
      SITE(30,:)  = (/-4.621780D0,  1.245370D0,  0.000000D0/)   !H6
      SITE(31,:)  = (/-4.621780D0,  -1.245370D0,  0.000000D0/)   !H7
      SITE(32,:)  = (/-3.389420D0,  -3.379890D0,  0.000000D0/)   !H8
      SITE(33,:)  = (/-1.232360D0,  -4.625280D0,  0.000000D0/)   !H9
      SITE(34,:)  = (/1.232370D0,  -4.625280D0,  0.000000D0/)   !H10
      SITE(35,:)  = (/3.389430D0,  -3.379890D0,  0.000000D0/)   !H11
      SITE(36,:)  = (/4.621790D0,  -1.245370D0,  0.000000D0/)   !H12

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0110220D0
      STCHRG(2)  = -0.0110220D0
      STCHRG(3)  = -0.0110220D0
      STCHRG(4)  = -0.0110220D0
      STCHRG(5)  = -0.0110220D0
      STCHRG(6)  = -0.0110220D0
      STCHRG(7)  = 0.1817100D0
      STCHRG(8)  = 0.1817100D0
      STCHRG(9)  = 0.1817100D0
      STCHRG(10)  = 0.1817100D0
      STCHRG(11)  = 0.1817100D0
      STCHRG(12)  = 0.1817100D0
      STCHRG(13)  = -0.2397720D0
      STCHRG(14)  = -0.2397720D0
      STCHRG(15)  = -0.2397720D0
      STCHRG(16)  = -0.2397720D0
      STCHRG(17)  = -0.2397720D0
      STCHRG(18)  = -0.2397720D0
      STCHRG(19)  = -0.2397720D0
      STCHRG(20)  = -0.2397720D0
      STCHRG(21)  = -0.2397720D0
      STCHRG(22)  = -0.2397720D0
      STCHRG(23)  = -0.2397720D0
      STCHRG(24)  = -0.2397720D0
      STCHRG(25)  = 0.1544280D0
      STCHRG(26)  = 0.1544280D0
      STCHRG(27)  = 0.1544280D0
      STCHRG(28)  = 0.1544280D0
      STCHRG(29)  = 0.1544280D0
      STCHRG(30)  = 0.1544280D0
      STCHRG(31)  = 0.1544280D0
      STCHRG(32)  = 0.1544280D0
      STCHRG(33)  = 0.1544280D0
      STCHRG(34)  = 0.1544280D0
      STCHRG(35)  = 0.1544280D0
      STCHRG(36)  = 0.1544280D0

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO C7
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO C8
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO C9
      RBSTLA(4,:)  = SITE(10,:)  - SITE(4,:)                 ! Z FROM C4 TO C10
      RBSTLA(5,:)  = SITE(11,:)  - SITE(5,:)                 ! Z FROM C5 TO C11
      RBSTLA(6,:)  = SITE(12,:)  - SITE(6,:)                 ! Z FROM C6 TO C12
      RBSTLA(7,:)  = SITE(1,:)  - SITE(7,:)                 ! Z FROM C7 TO C1
      RBSTLA(8,:)  = SITE(2,:)  - SITE(8,:)                 ! Z FROM C8 TO C2
      RBSTLA(9,:)  = SITE(3,:)  - SITE(9,:)                 ! Z FROM C9 TO C3
      RBSTLA(10,:)  = SITE(4,:)  - SITE(10,:)                 ! Z FROM C10 TO C4
      RBSTLA(11,:)  = SITE(5,:)  - SITE(11,:)                 ! Z FROM C11 TO C5
      RBSTLA(12,:)  = SITE(6,:)  - SITE(12,:)                 ! Z FROM C12 TO C6
      RBSTLA(13,:)  = SITE(25,:)  - SITE(13,:)                 ! Z FROM C13 TO H1
      RBSTLA(14,:)  = SITE(26,:)  - SITE(14,:)                 ! Z FROM C14 TO H2
      RBSTLA(15,:)  = SITE(27,:)  - SITE(15,:)                 ! Z FROM C15 TO H3
      RBSTLA(16,:)  = SITE(28,:)  - SITE(16,:)                 ! Z FROM C16 TO H4
      RBSTLA(17,:)  = SITE(29,:)  - SITE(17,:)                 ! Z FROM C17 TO H5
      RBSTLA(18,:)  = SITE(30,:)  - SITE(18,:)                 ! Z FROM C18 TO H6
      RBSTLA(19,:)  = SITE(31,:)  - SITE(19,:)                 ! Z FROM C19 TO H7
      RBSTLA(20,:)  = SITE(32,:)  - SITE(20,:)                 ! Z FROM C20 TO H8
      RBSTLA(21,:)  = SITE(33,:)  - SITE(21,:)                 ! Z FROM C21 TO H9
      RBSTLA(22,:)  = SITE(34,:)  - SITE(22,:)                 ! Z FROM C22 TO H10
      RBSTLA(23,:)  = SITE(35,:)  - SITE(23,:)                 ! Z FROM C23 TO H11
      RBSTLA(24,:)  = SITE(36,:)  - SITE(24,:)                 ! Z FROM C24 TO H12
      RBSTLA(25,:)  = SITE(25,:)  - SITE(13,:)                 ! Z FROM C13 TO H1
      RBSTLA(26,:)  = SITE(26,:)  - SITE(14,:)                 ! Z FROM C14 TO H2
      RBSTLA(27,:)  = SITE(27,:)  - SITE(15,:)                 ! Z FROM C15 TO H3
      RBSTLA(28,:)  = SITE(28,:)  - SITE(16,:)                 ! Z FROM C16 TO H4
      RBSTLA(29,:)  = SITE(29,:)  - SITE(17,:)                 ! Z FROM C17 TO H5
      RBSTLA(30,:)  = SITE(30,:)  - SITE(18,:)                 ! Z FROM C18 TO H6
      RBSTLA(31,:)  = SITE(31,:)  - SITE(19,:)                 ! Z FROM C19 TO H7
      RBSTLA(32,:)  = SITE(32,:)  - SITE(20,:)                 ! Z FROM C20 TO H8
      RBSTLA(33,:)  = SITE(33,:)  - SITE(21,:)                 ! Z FROM C21 TO H9
      RBSTLA(34,:)  = SITE(34,:)  - SITE(22,:)                 ! Z FROM C22 TO H10
      RBSTLA(35,:)  = SITE(35,:)  - SITE(23,:)                 ! Z FROM C23 TO H11
      RBSTLA(36,:)  = SITE(36,:)  - SITE(24,:)                 ! Z FROM C24 TO H12

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFCORONENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBISANTHENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C28H14

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.726050D0,  0.000000D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/-0.726050D0,  0.000000D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/1.431660D0,  -1.230140D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/-1.431660D0,  -1.230140D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/1.431660D0,  1.230140D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-1.431660D0,  1.230140D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/0.736110D0,  -2.487820D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/-0.736110D0,  -2.487820D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/0.736110D0,  2.487820D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-0.736110D0,  2.487820D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/2.872800D0,  -1.222330D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.872800D0,  -1.222330D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/2.872800D0,  1.222330D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/-2.872800D0,  1.222330D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/3.552990D0,  0.000000D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/-3.552990D0,  0.000000D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/3.579590D0,  -2.458470D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/-3.579590D0,  -2.458470D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/3.579590D0,  2.458470D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-3.579590D0,  2.458470D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/2.894490D0,  -3.647210D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/-2.894500D0,  -3.647210D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/2.894500D0,  3.647210D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/-2.894500D0,  3.647210D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/1.483700D0,  -3.659590D0,  0.000000D0/)   !C25
      SITE(26,:)  = (/-1.483700D0,  -3.659590D0,  0.000000D0/)   !C26
      SITE(27,:)  = (/1.483700D0,  3.659590D0,  0.000000D0/)   !C27
      SITE(28,:)  = (/-1.483700D0,  3.659590D0,  0.000000D0/)   !C28
      SITE(29,:)  = (/4.640690D0,  0.000000D0,  0.000000D0/)   !H1
      SITE(30,:)  = (/-4.640690D0,  0.000000D0,  0.000000D0/)   !H2
      SITE(31,:)  = (/4.666360D0,  -2.438830D0,  0.000000D0/)   !H3
      SITE(32,:)  = (/-4.666360D0,  -2.438830D0,  0.000000D0/)   !H4
      SITE(33,:)  = (/4.666360D0,  2.438830D0,  0.000000D0/)   !H5
      SITE(34,:)  = (/-4.666360D0,  2.438830D0,  0.000000D0/)   !H6
      SITE(35,:)  = (/3.432380D0,  -4.591340D0,  0.000000D0/)   !H7
      SITE(36,:)  = (/-3.432380D0,  -4.591340D0,  0.000000D0/)   !H8
      SITE(37,:)  = (/3.432380D0,  4.591340D0,  0.000000D0/)   !H9
      SITE(38,:)  = (/-3.432380D0,  4.591340D0,  0.000000D0/)   !H10
      SITE(39,:)  = (/0.986740D0,  -4.623150D0,  0.000000D0/)   !H11
      SITE(40,:)  = (/-0.986740D0,  -4.623150D0,  0.000000D0/)   !H12
      SITE(41,:)  = (/0.986740D0,  4.623150D0,  0.000000D0/)   !H13
      SITE(42,:)  = (/-0.986740D0,  4.623150D0,  0.000000D0/)   !H14

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0461580D0
      STCHRG(2)  = -0.0461580D0
      STCHRG(3)  = -0.0059470D0
      STCHRG(4)  = -0.0059470D0
      STCHRG(5)  = -0.0059470D0
      STCHRG(6)  = -0.0059470D0
      STCHRG(7)  = 0.0554270D0
      STCHRG(8)  = 0.0554270D0
      STCHRG(9)  = 0.0554270D0
      STCHRG(10)  = 0.0554270D0
      STCHRG(11)  = 0.2774230D0
      STCHRG(12)  = 0.2774230D0
      STCHRG(13)  = 0.2774230D0
      STCHRG(14)  = 0.2774230D0
      STCHRG(15)  = -0.4862500D0
      STCHRG(16)  = -0.4862500D0
      STCHRG(17)  = -0.2731650D0
      STCHRG(18)  = -0.2731650D0
      STCHRG(19)  = -0.2731650D0
      STCHRG(20)  = -0.2731650D0
      STCHRG(21)  = -0.1022400D0
      STCHRG(22)  = -0.1022400D0
      STCHRG(23)  = -0.1022400D0
      STCHRG(24)  = -0.1022400D0
      STCHRG(25)  = -0.2224090D0
      STCHRG(26)  = -0.2224090D0
      STCHRG(27)  = -0.2224090D0
      STCHRG(28)  = -0.2224090D0
      STCHRG(29)  = 0.1992400D0
      STCHRG(30)  = 0.1992400D0
      STCHRG(31)  = 0.1538890D0
      STCHRG(32)  = 0.1538890D0
      STCHRG(33)  = 0.1538890D0
      STCHRG(34)  = 0.1538890D0
      STCHRG(35)  = 0.1352520D0
      STCHRG(36)  = 0.1352520D0
      STCHRG(37)  = 0.1352520D0
      STCHRG(38)  = 0.1352520D0
      STCHRG(39)  = 0.1483550D0
      STCHRG(40)  = 0.1483550D0
      STCHRG(41)  = 0.1483550D0
      STCHRG(42)  = 0.1483550D0

      RBSTLA(1,:)  = SITE(2,:)  - SITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = SITE(1,:)  - SITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = SITE(1,:)  - SITE(3,:)                 ! Z FROM C3 TO C1
      RBSTLA(4,:)  = SITE(2,:)  - SITE(4,:)                 ! Z FROM C4 TO C2
      RBSTLA(5,:)  = SITE(1,:)  - SITE(5,:)                 ! Z FROM C5 TO C1
      RBSTLA(6,:)  = SITE(2,:)  - SITE(6,:)                 ! Z FROM C6 TO C2
      RBSTLA(7,:)  = SITE(8,:)  - SITE(7,:)                 ! Z FROM C7 TO C8
      RBSTLA(8,:)  = SITE(7,:)  - SITE(8,:)                 ! Z FROM C8 TO C7
      RBSTLA(9,:)  = SITE(10,:)  - SITE(9,:)                 ! Z FROM C9 TO C10
      RBSTLA(10,:)  = SITE(9,:)  - SITE(10,:)                 ! Z FROM C10 TO C9
      RBSTLA(11,:)  = SITE(3,:)  - SITE(11,:)                 ! Z FROM C11 TO C3
      RBSTLA(12,:)  = SITE(4,:)  - SITE(12,:)                 ! Z FROM C12 TO C4
      RBSTLA(13,:)  = SITE(5,:)  - SITE(13,:)                 ! Z FROM C13 TO C5
      RBSTLA(14,:)  = SITE(6,:)  - SITE(14,:)                 ! Z FROM C14 TO C6
      RBSTLA(15,:)  = SITE(29,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(16,:)  = SITE(30,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(17,:)  = SITE(31,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(18,:)  = SITE(32,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(19,:)  = SITE(33,:)  - SITE(19,:)                 ! Z FROM C19 TO H5
      RBSTLA(20,:)  = SITE(34,:)  - SITE(20,:)                 ! Z FROM C20 TO H6
      RBSTLA(21,:)  = SITE(35,:)  - SITE(21,:)                 ! Z FROM C21 TO H7
      RBSTLA(22,:)  = SITE(36,:)  - SITE(22,:)                 ! Z FROM C22 TO H8
      RBSTLA(23,:)  = SITE(37,:)  - SITE(23,:)                 ! Z FROM C23 TO H9
      RBSTLA(24,:)  = SITE(38,:)  - SITE(24,:)                 ! Z FROM C24 TO H10
      RBSTLA(25,:)  = SITE(39,:)  - SITE(25,:)                 ! Z FROM C25 TO H11
      RBSTLA(26,:)  = SITE(40,:)  - SITE(26,:)                 ! Z FROM C26 TO H12
      RBSTLA(27,:)  = SITE(41,:)  - SITE(27,:)                 ! Z FROM C27 TO H13
      RBSTLA(28,:)  = SITE(42,:)  - SITE(28,:)                 ! Z FROM C28 TO H14
      RBSTLA(29,:)  = SITE(29,:)  - SITE(15,:)                 ! Z FROM C15 TO H1
      RBSTLA(30,:)  = SITE(30,:)  - SITE(16,:)                 ! Z FROM C16 TO H2
      RBSTLA(31,:)  = SITE(31,:)  - SITE(17,:)                 ! Z FROM C17 TO H3
      RBSTLA(32,:)  = SITE(32,:)  - SITE(18,:)                 ! Z FROM C18 TO H4
      RBSTLA(33,:)  = SITE(33,:)  - SITE(19,:)                 ! Z FROM C19 TO H5
      RBSTLA(34,:)  = SITE(34,:)  - SITE(20,:)                 ! Z FROM C20 TO H6
      RBSTLA(35,:)  = SITE(35,:)  - SITE(21,:)                 ! Z FROM C21 TO H7
      RBSTLA(36,:)  = SITE(36,:)  - SITE(22,:)                 ! Z FROM C22 TO H8
      RBSTLA(37,:)  = SITE(37,:)  - SITE(23,:)                 ! Z FROM C23 TO H9
      RBSTLA(38,:)  = SITE(38,:)  - SITE(24,:)                 ! Z FROM C24 TO H10
      RBSTLA(39,:)  = SITE(39,:)  - SITE(25,:)                 ! Z FROM C25 TO H11
      RBSTLA(40,:)  = SITE(40,:)  - SITE(26,:)                 ! Z FROM C26 TO H12
      RBSTLA(41,:)  = SITE(41,:)  - SITE(27,:)                 ! Z FROM C27 TO H13
      RBSTLA(42,:)  = SITE(42,:)  - SITE(28,:)                 ! Z FROM C28 TO H14

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFBISANTHENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFOVALENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C32H14

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/-1.228120D0,  -1.426290D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/0.000000D0,  -0.717510D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/0.000000D0,  0.717510D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/-1.228120D0,  1.426290D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-2.464490D0,  0.713290D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-2.464490D0,  -0.713290D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/-1.224540D0,  -2.858360D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/1.228120D0,  -1.426290D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/1.228120D0,  1.426290D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/-1.224540D0,  2.858360D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-3.697210D0,  1.423190D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-3.697210D0,  -1.423190D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/1.224530D0,  -2.858360D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/2.464490D0,  -0.713290D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/2.464490D0,  0.713290D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/1.224530D0,  2.858360D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/-3.665710D0,  2.855940D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/-4.907220D0,  0.690170D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-4.907220D0,  -0.690170D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-3.665710D0,  -2.855940D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/0.000000D0,  -3.539310D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/3.697210D0,  -1.423190D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/3.697210D0,  1.423190D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/0.000000D0,  3.539310D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/-2.485870D0,  3.542240D0,  0.000000D0/)   !C25
      SITE(26,:)  = (/-2.485870D0,  -3.542240D0,  0.000000D0/)   !C26
      SITE(27,:)  = (/2.485870D0,  -3.542240D0,  0.000000D0/)   !C27
      SITE(28,:)  = (/3.665710D0,  -2.855940D0,  0.000000D0/)   !C28
      SITE(29,:)  = (/4.907220D0,  -0.690170D0,  0.000000D0/)   !C29
      SITE(30,:)  = (/4.907220D0,  0.690170D0,  0.000000D0/)   !C30
      SITE(31,:)  = (/3.665710D0,  2.855940D0,  0.000000D0/)   !C31
      SITE(32,:)  = (/2.485870D0,  3.542240D0,  0.000000D0/)   !C32
      SITE(33,:)  = (/0.000000D0,  -4.627370D0,  0.000000D0/)   !H1
      SITE(34,:)  = (/0.000000D0,  4.627370D0,  0.000000D0/)   !H2
      SITE(35,:)  = (/-2.483580D0,  4.629530D0,  0.000000D0/)   !H3
      SITE(36,:)  = (/-2.483580D0,  -4.629530D0,  0.000000D0/)   !H4
      SITE(37,:)  = (/-5.849040D0,  1.233670D0,  0.000000D0/)   !H5
      SITE(38,:)  = (/-5.849040D0,  -1.233670D0,  0.000000D0/)   !H6
      SITE(39,:)  = (/-4.611080D0,  3.393330D0,  0.000000D0/)   !H7
      SITE(40,:)  = (/-4.611080D0,  -3.393330D0,  0.000000D0/)   !H8
      SITE(41,:)  = (/2.483580D0,  -4.629530D0,  0.000000D0/)   !H9
      SITE(42,:)  = (/2.483580D0,  4.629530D0,  0.000000D0/)   !H10
      SITE(43,:)  = (/4.611080D0,  -3.393330D0,  0.000000D0/)   !H11
      SITE(44,:)  = (/5.849040D0,  -1.233670D0,  0.000000D0/)   !H12
      SITE(45,:)  = (/5.849040D0,  1.233670D0,  0.000000D0/)   !H13
      SITE(46,:)  = (/4.611080D0,  3.393330D0,  0.000000D0/)   !H14

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0564090D0
      STCHRG(2)  = 0.0160840D0
      STCHRG(3)  = 0.0160840D0
      STCHRG(4)  = -0.0564090D0
      STCHRG(5)  = 0.0101710D0
      STCHRG(6)  = 0.0101710D0
      STCHRG(7)  = 0.2426000D0
      STCHRG(8)  = -0.0564090D0
      STCHRG(9)  = -0.0564090D0
      STCHRG(10)  = 0.2426000D0
      STCHRG(11)  = 0.1609920D0
      STCHRG(12)  = 0.1609920D0
      STCHRG(13)  = 0.2426000D0
      STCHRG(14)  = 0.0101710D0
      STCHRG(15)  = 0.0101710D0
      STCHRG(16)  = 0.2426000D0
      STCHRG(17)  = -0.2434940D0
      STCHRG(18)  = -0.2338890D0
      STCHRG(19)  = -0.2338890D0
      STCHRG(20)  = -0.2434940D0
      STCHRG(21)  = -0.4310520D0
      STCHRG(22)  = 0.1609920D0
      STCHRG(23)  = 0.1609920D0
      STCHRG(24)  = -0.4310520D0
      STCHRG(25)  = -0.2319160D0
      STCHRG(26)  = -0.2319160D0
      STCHRG(27)  = -0.2319160D0
      STCHRG(28)  = -0.2434940D0
      STCHRG(29)  = -0.2338890D0
      STCHRG(30)  = -0.2338890D0
      STCHRG(31)  = -0.2434940D0
      STCHRG(32)  = -0.2319160D0
      STCHRG(33)  = 0.1797660D0
      STCHRG(34)  = 0.1797660D0
      STCHRG(35)  = 0.1566500D0
      STCHRG(36)  = 0.1566500D0
      STCHRG(37)  = 0.1565740D0
      STCHRG(38)  = 0.1565740D0
      STCHRG(39)  = 0.1563220D0
      STCHRG(40)  = 0.1563220D0
      STCHRG(41)  = 0.1566500D0
      STCHRG(42)  = 0.1566500D0
      STCHRG(43)  = 0.1563220D0
      STCHRG(44)  = 0.1565740D0
      STCHRG(45)  = 0.1565740D0
      STCHRG(46)  = 0.1563220D0

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO C7
      RBSTLA(2,:)  = SITE(3,:)  - SITE(2,:)                 ! Z FROM C2 TO C3
      RBSTLA(3,:)  = SITE(2,:)  - SITE(3,:)                 ! Z FROM C3 TO C2
      RBSTLA(4,:)  = SITE(10,:)  - SITE(4,:)                 ! Z FROM C4 TO C10
      RBSTLA(5,:)  = SITE(11,:)  - SITE(5,:)                 ! Z FROM C5 TO C11
      RBSTLA(6,:)  = SITE(12,:)  - SITE(6,:)                 ! Z FROM C6 TO C12
      RBSTLA(7,:)  = SITE(1,:)  - SITE(7,:)                 ! Z FROM C7 TO C1
      RBSTLA(8,:)  = SITE(13,:)  - SITE(8,:)                 ! Z FROM C8 TO C13
      RBSTLA(9,:)  = SITE(16,:)  - SITE(9,:)                 ! Z FROM C9 TO C16
      RBSTLA(10,:)  = SITE(4,:)  - SITE(10,:)                 ! Z FROM C10 TO C4
      RBSTLA(11,:)  = SITE(5,:)  - SITE(11,:)                 ! Z FROM C11 TO C5
      RBSTLA(12,:)  = SITE(6,:)  - SITE(12,:)                 ! Z FROM C12 TO C6
      RBSTLA(13,:)  = SITE(8,:)  - SITE(13,:)                 ! Z FROM C13 TO C8
      RBSTLA(14,:)  = SITE(22,:)  - SITE(14,:)                 ! Z FROM C14 TO C22
      RBSTLA(15,:)  = SITE(23,:)  - SITE(15,:)                 ! Z FROM C15 TO C23
      RBSTLA(16,:)  = SITE(9,:)  - SITE(16,:)                 ! Z FROM C16 TO C9
      RBSTLA(17,:)  = SITE(39,:)  - SITE(17,:)                 ! Z FROM C17 TO H7
      RBSTLA(18,:)  = SITE(37,:)  - SITE(18,:)                 ! Z FROM C18 TO H5
      RBSTLA(19,:)  = SITE(38,:)  - SITE(19,:)                 ! Z FROM C19 TO H6
      RBSTLA(20,:)  = SITE(40,:)  - SITE(20,:)                 ! Z FROM C20 TO H8
      RBSTLA(21,:)  = SITE(33,:)  - SITE(21,:)                 ! Z FROM C21 TO H1
      RBSTLA(22,:)  = SITE(14,:)  - SITE(22,:)                 ! Z FROM C22 TO C14
      RBSTLA(23,:)  = SITE(15,:)  - SITE(23,:)                 ! Z FROM C23 TO C15
      RBSTLA(24,:)  = SITE(34,:)  - SITE(24,:)                 ! Z FROM C24 TO H2
      RBSTLA(25,:)  = SITE(35,:)  - SITE(25,:)                 ! Z FROM C25 TO H3
      RBSTLA(26,:)  = SITE(36,:)  - SITE(26,:)                 ! Z FROM C26 TO H4
      RBSTLA(27,:)  = SITE(41,:)  - SITE(27,:)                 ! Z FROM C27 TO H9
      RBSTLA(28,:)  = SITE(43,:)  - SITE(28,:)                 ! Z FROM C28 TO H11
      RBSTLA(29,:)  = SITE(44,:)  - SITE(29,:)                 ! Z FROM C29 TO H12
      RBSTLA(30,:)  = SITE(45,:)  - SITE(30,:)                 ! Z FROM C30 TO H13
      RBSTLA(31,:)  = SITE(46,:)  - SITE(31,:)                 ! Z FROM C31 TO H14
      RBSTLA(32,:)  = SITE(42,:)  - SITE(32,:)                 ! Z FROM C32 TO H10
      RBSTLA(33,:)  = SITE(33,:)  - SITE(21,:)                 ! Z FROM C21 TO H1
      RBSTLA(34,:)  = SITE(34,:)  - SITE(24,:)                 ! Z FROM C24 TO H2
      RBSTLA(35,:)  = SITE(35,:)  - SITE(25,:)                 ! Z FROM C25 TO H3
      RBSTLA(36,:)  = SITE(36,:)  - SITE(26,:)                 ! Z FROM C26 TO H4
      RBSTLA(37,:)  = SITE(37,:)  - SITE(18,:)                 ! Z FROM C18 TO H5
      RBSTLA(38,:)  = SITE(38,:)  - SITE(19,:)                 ! Z FROM C19 TO H6
      RBSTLA(39,:)  = SITE(39,:)  - SITE(17,:)                 ! Z FROM C17 TO H7
      RBSTLA(40,:)  = SITE(40,:)  - SITE(20,:)                 ! Z FROM C20 TO H8
      RBSTLA(41,:)  = SITE(41,:)  - SITE(27,:)                 ! Z FROM C27 TO H9
      RBSTLA(42,:)  = SITE(42,:)  - SITE(32,:)                 ! Z FROM C32 TO H10
      RBSTLA(43,:)  = SITE(43,:)  - SITE(28,:)                 ! Z FROM C28 TO H11
      RBSTLA(44,:)  = SITE(44,:)  - SITE(29,:)                 ! Z FROM C29 TO H12
      RBSTLA(45,:)  = SITE(45,:)  - SITE(30,:)                 ! Z FROM C30 TO H13
      RBSTLA(46,:)  = SITE(46,:)  - SITE(31,:)                 ! Z FROM C31 TO H14

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFOVALENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFHEXABENZOCORONENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C42H18

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.000000D0,  -1.423760D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/1.233020D0,  -0.711880D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/1.233020D0,  0.711880D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/0.000000D0,  1.423760D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-1.233020D0,  0.711880D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-1.233020D0,  -0.711880D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/0.000000D0,  -2.871600D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/2.486890D0,  -1.435800D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/2.486890D0,  1.435800D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/0.000000D0,  2.871600D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-2.486890D0,  1.435800D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.486890D0,  -1.435800D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/1.200840D0,  -4.993810D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/2.495340D0,  -2.860590D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/3.725020D0,  -0.730730D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/3.725020D0,  0.730730D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/2.495340D0,  2.860590D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/1.200840D0,  4.993810D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-1.200840D0,  4.993810D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-2.495340D0,  2.860590D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/-3.725020D0,  0.730730D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/-3.725020D0,  -0.730730D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/-2.495340D0,  -2.860590D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/-1.200840D0,  -4.993810D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/1.229670D0,  -3.591320D0,  0.000000D0/)   !C25
      SITE(26,:)  = (/4.925190D0,  -1.456950D0,  0.000000D0/)   !C26
      SITE(27,:)  = (/4.925190D0,  1.456950D0,  0.000000D0/)   !C27
      SITE(28,:)  = (/1.229670D0,  3.591320D0,  0.000000D0/)   !C28
      SITE(29,:)  = (/-1.229670D0,  3.591320D0,  0.000000D0/)   !C29
      SITE(30,:)  = (/-4.925190D0,  1.456950D0,  0.000000D0/)   !C30
      SITE(31,:)  = (/-4.925190D0,  -1.456950D0,  0.000000D0/)   !C31
      SITE(32,:)  = (/-1.229670D0,  -3.591320D0,  0.000000D0/)   !C32
      SITE(33,:)  = (/3.724340D0,  -3.536870D0,  0.000000D0/)   !C33
      SITE(34,:)  = (/4.925490D0,  -2.843740D0,  0.000000D0/)   !C34
      SITE(35,:)  = (/4.925490D0,  2.843740D0,  0.000000D0/)   !C35
      SITE(36,:)  = (/3.724340D0,  3.536870D0,  0.000000D0/)   !C36
      SITE(37,:)  = (/0.000000D0,  5.687460D0,  0.000000D0/)   !C37
      SITE(38,:)  = (/-3.724340D0,  3.536870D0,  0.000000D0/)   !C38
      SITE(39,:)  = (/-4.925490D0,  2.843740D0,  0.000000D0/)   !C39
      SITE(40,:)  = (/-4.925490D0,  -2.843740D0,  0.000000D0/)   !C40
      SITE(41,:)  = (/-3.724340D0,  -3.536870D0,  0.000000D0/)   !C41
      SITE(42,:)  = (/0.000000D0,  -5.687460D0,  0.000000D0/)   !C42
      SITE(43,:)  = (/-2.121850D0,  5.562640D0,  0.000000D0/)   !H1
      SITE(44,:)  = (/2.121850D0,  -5.562640D0,  0.000000D0/)   !H2
      SITE(45,:)  = (/3.756480D0,  -4.618910D0,  0.000000D0/)   !H3
      SITE(46,:)  = (/3.756480D0,  4.618910D0,  0.000000D0/)   !H4
      SITE(47,:)  = (/-3.756480D0,  4.618910D0,  0.000000D0/)   !H5
      SITE(48,:)  = (/-3.756480D0,  -4.618910D0,  0.000000D0/)   !H6
      SITE(49,:)  = (/5.878330D0,  -0.943780D0,  0.000000D0/)   !H7
      SITE(50,:)  = (/5.878330D0,  0.943780D0,  0.000000D0/)   !H8
      SITE(51,:)  = (/2.121850D0,  5.562640D0,  0.000000D0/)   !H9
      SITE(52,:)  = (/-5.878330D0,  0.943780D0,  0.000000D0/)   !H10
      SITE(53,:)  = (/-5.878330D0,  -0.943780D0,  0.000000D0/)   !H11
      SITE(54,:)  = (/-2.121850D0,  -5.562640D0,  0.000000D0/)   !H12
      SITE(55,:)  = (/5.866390D0,  -3.386970D0,  0.000000D0/)   !H13
      SITE(56,:)  = (/5.866390D0,  3.386970D0,  0.000000D0/)   !H14
      SITE(57,:)  = (/0.000000D0,  6.773920D0,  0.000000D0/)   !H15
      SITE(58,:)  = (/-5.866390D0,  3.386970D0,  0.000000D0/)   !H16
      SITE(59,:)  = (/-5.866390D0,  -3.386970D0,  0.000000D0/)   !H17
      SITE(60,:)  = (/0.000000D0,  -6.773920D0,  0.000000D0/)   !H18

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0239090D0
      STCHRG(2)  = -0.0239090D0
      STCHRG(3)  = -0.0239090D0
      STCHRG(4)  = -0.0239090D0
      STCHRG(5)  = -0.0239090D0
      STCHRG(6)  = -0.0239090D0
      STCHRG(7)  = 0.0287860D0
      STCHRG(8)  = 0.0287860D0
      STCHRG(9)  = 0.0287860D0
      STCHRG(10)  = 0.0287860D0
      STCHRG(11)  = 0.0287860D0
      STCHRG(12)  = 0.0287860D0
      STCHRG(13)  = -0.2113620D0
      STCHRG(14)  = 0.0526040D0
      STCHRG(15)  = 0.0526040D0
      STCHRG(16)  = 0.0526040D0
      STCHRG(17)  = 0.0526040D0
      STCHRG(18)  = -0.2113620D0
      STCHRG(19)  = -0.2113620D0
      STCHRG(20)  = 0.0526040D0
      STCHRG(21)  = 0.0526040D0
      STCHRG(22)  = 0.0526040D0
      STCHRG(23)  = 0.0526040D0
      STCHRG(24)  = -0.2113620D0
      STCHRG(25)  = 0.0526040D0
      STCHRG(26)  = -0.2113620D0
      STCHRG(27)  = -0.2113620D0
      STCHRG(28)  = 0.0526040D0
      STCHRG(29)  = 0.0526040D0
      STCHRG(30)  = -0.2113620D0
      STCHRG(31)  = -0.2113620D0
      STCHRG(32)  = 0.0526040D0
      STCHRG(33)  = -0.2113620D0
      STCHRG(34)  = -0.1107910D0
      STCHRG(35)  = -0.1107910D0
      STCHRG(36)  = -0.2113620D0
      STCHRG(37)  = -0.1107910D0
      STCHRG(38)  = -0.2113620D0
      STCHRG(39)  = -0.1107910D0
      STCHRG(40)  = -0.1107910D0
      STCHRG(41)  = -0.2113620D0
      STCHRG(42)  = -0.1107910D0
      STCHRG(43)  = 0.1456340D0
      STCHRG(44)  = 0.1456340D0
      STCHRG(45)  = 0.1456340D0
      STCHRG(46)  = 0.1456340D0
      STCHRG(47)  = 0.1456340D0
      STCHRG(48)  = 0.1456340D0
      STCHRG(49)  = 0.1456340D0
      STCHRG(50)  = 0.1456340D0
      STCHRG(51)  = 0.1456340D0
      STCHRG(52)  = 0.1456340D0
      STCHRG(53)  = 0.1456340D0
      STCHRG(54)  = 0.1456340D0
      STCHRG(55)  = 0.1321610D0
      STCHRG(56)  = 0.1321610D0
      STCHRG(57)  = 0.1321610D0
      STCHRG(58)  = 0.1321610D0
      STCHRG(59)  = 0.1321610D0
      STCHRG(60)  = 0.1321610D0

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO C7
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO C8
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO C9
      RBSTLA(4,:)  = SITE(10,:)  - SITE(4,:)                 ! Z FROM C4 TO C10
      RBSTLA(5,:)  = SITE(11,:)  - SITE(5,:)                 ! Z FROM C5 TO C11
      RBSTLA(6,:)  = SITE(12,:)  - SITE(6,:)                 ! Z FROM C6 TO C12
      RBSTLA(7,:)  = SITE(1,:)  - SITE(7,:)                 ! Z FROM C7 TO C1
      RBSTLA(8,:)  = SITE(2,:)  - SITE(8,:)                 ! Z FROM C8 TO C2
      RBSTLA(9,:)  = SITE(3,:)  - SITE(9,:)                 ! Z FROM C9 TO C3
      RBSTLA(10,:)  = SITE(4,:)  - SITE(10,:)                 ! Z FROM C10 TO C4
      RBSTLA(11,:)  = SITE(5,:)  - SITE(11,:)                 ! Z FROM C11 TO C5
      RBSTLA(12,:)  = SITE(6,:)  - SITE(12,:)                 ! Z FROM C12 TO C6
      RBSTLA(13,:)  = SITE(44,:)  - SITE(13,:)                 ! Z FROM C13 TO H2
      RBSTLA(14,:)  = SITE(25,:)  - SITE(14,:)                 ! Z FROM C14 TO C25
      RBSTLA(15,:)  = SITE(16,:)  - SITE(15,:)                 ! Z FROM C15 TO C16
      RBSTLA(16,:)  = SITE(15,:)  - SITE(16,:)                 ! Z FROM C16 TO C15
      RBSTLA(17,:)  = SITE(28,:)  - SITE(17,:)                 ! Z FROM C17 TO C28
      RBSTLA(18,:)  = SITE(51,:)  - SITE(18,:)                 ! Z FROM C18 TO H9
      RBSTLA(19,:)  = SITE(43,:)  - SITE(19,:)                 ! Z FROM C19 TO H1
      RBSTLA(20,:)  = SITE(29,:)  - SITE(20,:)                 ! Z FROM C20 TO C29
      RBSTLA(21,:)  = SITE(22,:)  - SITE(21,:)                 ! Z FROM C21 TO C22
      RBSTLA(22,:)  = SITE(21,:)  - SITE(22,:)                 ! Z FROM C22 TO C21
      RBSTLA(23,:)  = SITE(32,:)  - SITE(23,:)                 ! Z FROM C23 TO C32
      RBSTLA(24,:)  = SITE(54,:)  - SITE(24,:)                 ! Z FROM C24 TO H12
      RBSTLA(25,:)  = SITE(14,:)  - SITE(25,:)                 ! Z FROM C25 TO C14
      RBSTLA(26,:)  = SITE(49,:)  - SITE(26,:)                 ! Z FROM C26 TO H7
      RBSTLA(27,:)  = SITE(50,:)  - SITE(27,:)                 ! Z FROM C27 TO H8
      RBSTLA(28,:)  = SITE(17,:)  - SITE(28,:)                 ! Z FROM C28 TO C17
      RBSTLA(29,:)  = SITE(20,:)  - SITE(29,:)                 ! Z FROM C29 TO C20
      RBSTLA(30,:)  = SITE(52,:)  - SITE(30,:)                 ! Z FROM C30 TO H10
      RBSTLA(31,:)  = SITE(53,:)  - SITE(31,:)                 ! Z FROM C31 TO H11
      RBSTLA(32,:)  = SITE(23,:)  - SITE(32,:)                 ! Z FROM C32 TO C23
      RBSTLA(33,:)  = SITE(45,:)  - SITE(33,:)                 ! Z FROM C33 TO H3
      RBSTLA(34,:)  = SITE(55,:)  - SITE(34,:)                 ! Z FROM C34 TO H13
      RBSTLA(35,:)  = SITE(56,:)  - SITE(35,:)                 ! Z FROM C35 TO H14
      RBSTLA(36,:)  = SITE(46,:)  - SITE(36,:)                 ! Z FROM C36 TO H4
      RBSTLA(37,:)  = SITE(57,:)  - SITE(37,:)                 ! Z FROM C37 TO H15
      RBSTLA(38,:)  = SITE(47,:)  - SITE(38,:)                 ! Z FROM C38 TO H5
      RBSTLA(39,:)  = SITE(58,:)  - SITE(39,:)                 ! Z FROM C39 TO H16
      RBSTLA(40,:)  = SITE(59,:)  - SITE(40,:)                 ! Z FROM C40 TO H17
      RBSTLA(41,:)  = SITE(48,:)  - SITE(41,:)                 ! Z FROM C41 TO H6
      RBSTLA(42,:)  = SITE(60,:)  - SITE(42,:)                 ! Z FROM C42 TO H18
      RBSTLA(43,:)  = SITE(43,:)  - SITE(19,:)                 ! Z FROM C19 TO H1
      RBSTLA(44,:)  = SITE(44,:)  - SITE(13,:)                 ! Z FROM C13 TO H2
      RBSTLA(45,:)  = SITE(45,:)  - SITE(33,:)                 ! Z FROM C33 TO H3
      RBSTLA(46,:)  = SITE(46,:)  - SITE(36,:)                 ! Z FROM C36 TO H4
      RBSTLA(47,:)  = SITE(47,:)  - SITE(38,:)                 ! Z FROM C38 TO H5
      RBSTLA(48,:)  = SITE(48,:)  - SITE(41,:)                 ! Z FROM C41 TO H6
      RBSTLA(49,:)  = SITE(49,:)  - SITE(26,:)                 ! Z FROM C26 TO H7
      RBSTLA(50,:)  = SITE(50,:)  - SITE(27,:)                 ! Z FROM C27 TO H8
      RBSTLA(51,:)  = SITE(51,:)  - SITE(18,:)                 ! Z FROM C18 TO H9
      RBSTLA(52,:)  = SITE(52,:)  - SITE(30,:)                 ! Z FROM C30 TO H10
      RBSTLA(53,:)  = SITE(53,:)  - SITE(31,:)                 ! Z FROM C31 TO H11
      RBSTLA(54,:)  = SITE(54,:)  - SITE(24,:)                 ! Z FROM C24 TO H12
      RBSTLA(55,:)  = SITE(55,:)  - SITE(34,:)                 ! Z FROM C34 TO H13
      RBSTLA(56,:)  = SITE(56,:)  - SITE(35,:)                 ! Z FROM C35 TO H14
      RBSTLA(57,:)  = SITE(57,:)  - SITE(37,:)                 ! Z FROM C37 TO H15
      RBSTLA(58,:)  = SITE(58,:)  - SITE(39,:)                 ! Z FROM C39 TO H16
      RBSTLA(59,:)  = SITE(59,:)  - SITE(40,:)                 ! Z FROM C40 TO H17
      RBSTLA(60,:)  = SITE(60,:)  - SITE(42,:)                 ! Z FROM C42 TO H18

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFHEXABENZOCORONENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFOCTABENZOCORONENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C46H18

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.000000D0,  -1.425060D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/1.231180D0,  -0.710700D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/1.231180D0,  0.710700D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/0.000000D0,  1.425060D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-1.231180D0,  0.710700D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-1.231180D0,  -0.710700D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/0.000000D0,  -2.865070D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/2.478690D0,  -1.434890D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/2.478690D0,  1.434890D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/0.000000D0,  2.865070D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-2.478690D0,  1.434890D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.478690D0,  -1.434890D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/1.203620D0,  -4.987300D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/2.493820D0,  -2.849580D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/3.712050D0,  -0.714060D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/3.712050D0,  0.714060D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/2.493820D0,  2.849580D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/1.203620D0,  4.987300D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-1.203620D0,  4.987300D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-2.493820D0,  2.849580D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/-3.712050D0,  0.714060D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/-3.712050D0,  -0.714060D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/-2.493820D0,  -2.849580D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/-1.203620D0,  -4.987300D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/1.235310D0,  -3.583840D0,  0.000000D0/)   !C25
      SITE(26,:)  = (/4.951530D0,  -1.413420D0,  0.000000D0/)   !C26
      SITE(27,:)  = (/4.951530D0,  1.413420D0,  0.000000D0/)   !C27
      SITE(28,:)  = (/1.235310D0,  3.583840D0,  0.000000D0/)   !C28
      SITE(29,:)  = (/-1.235310D0,  3.583840D0,  0.000000D0/)   !C29
      SITE(30,:)  = (/-4.951530D0,  1.413420D0,  0.000000D0/)   !C30
      SITE(31,:)  = (/-4.951530D0,  -1.413420D0,  0.000000D0/)   !C31
      SITE(32,:)  = (/-1.235310D0,  -3.583840D0,  0.000000D0/)   !C32
      SITE(33,:)  = (/3.745490D0,  -3.515110D0,  0.000000D0/)   !C33
      SITE(34,:)  = (/4.933170D0,  -2.824900D0,  0.000000D0/)   !C34
      SITE(35,:)  = (/6.177280D0,  -0.682150D0,  0.000000D0/)   !C35
      SITE(36,:)  = (/6.177280D0,  0.682150D0,  0.000000D0/)   !C36
      SITE(37,:)  = (/4.933170D0,  2.824900D0,  0.000000D0/)   !C37
      SITE(38,:)  = (/3.745490D0,  3.515110D0,  0.000000D0/)   !C38
      SITE(39,:)  = (/0.000000D0,  5.678280D0,  0.000000D0/)   !C39
      SITE(40,:)  = (/-3.745490D0,  3.515110D0,  0.000000D0/)   !C40
      SITE(41,:)  = (/-4.933170D0,  2.824900D0,  0.000000D0/)   !C41
      SITE(42,:)  = (/-6.177280D0,  0.682150D0,  0.000000D0/)   !C42
      SITE(43,:)  = (/-6.177280D0,  -0.682150D0,  0.000000D0/)   !C43
      SITE(44,:)  = (/-4.933170D0,  -2.824900D0,  0.000000D0/)   !C44
      SITE(45,:)  = (/-3.745490D0,  -3.515110D0,  0.000000D0/)   !C45
      SITE(46,:)  = (/0.000000D0,  -5.678280D0,  0.000000D0/)   !C46
      SITE(47,:)  = (/3.783390D0,  -4.597500D0,  0.000000D0/)   !H1
      SITE(48,:)  = (/3.783390D0,  4.597500D0,  0.000000D0/)   !H2
      SITE(49,:)  = (/-3.783390D0,  4.597500D0,  0.000000D0/)   !H3
      SITE(50,:)  = (/-3.783390D0,  -4.597500D0,  0.000000D0/)   !H4
      SITE(51,:)  = (/2.124570D0,  -5.556960D0,  0.000000D0/)   !H5
      SITE(52,:)  = (/2.124570D0,  5.556960D0,  0.000000D0/)   !H6
      SITE(53,:)  = (/-2.124570D0,  5.556960D0,  0.000000D0/)   !H7
      SITE(54,:)  = (/-2.124570D0,  -5.556960D0,  0.000000D0/)   !H8
      SITE(55,:)  = (/5.876920D0,  -3.364720D0,  0.000000D0/)   !H9
      SITE(56,:)  = (/7.112310D0,  -1.236870D0,  0.000000D0/)   !H10
      SITE(57,:)  = (/7.112310D0,  1.236870D0,  0.000000D0/)   !H11
      SITE(58,:)  = (/5.876920D0,  3.364720D0,  0.000000D0/)   !H12
      SITE(59,:)  = (/0.000000D0,  6.764810D0,  0.000000D0/)   !H13
      SITE(60,:)  = (/-5.876920D0,  3.364720D0,  0.000000D0/)   !H14
      SITE(61,:)  = (/-7.112310D0,  1.236870D0,  0.000000D0/)   !H15
      SITE(62,:)  = (/-7.112310D0,  -1.236870D0,  0.000000D0/)   !H16
      SITE(63,:)  = (/-5.876920D0,  -3.364720D0,  0.000000D0/)   !H17
      SITE(64,:)  = (/0.000000D0,  -6.764810D0,  0.000000D0/)   !H18

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = -0.0628330D0
      STCHRG(2)  = 0.0094960D0
      STCHRG(3)  = 0.0094960D0
      STCHRG(4)  = -0.0628330D0
      STCHRG(5)  = 0.0094960D0
      STCHRG(6)  = 0.0094960D0
      STCHRG(7)  = 0.0231300D0
      STCHRG(8)  = 0.0110570D0
      STCHRG(9)  = 0.0110570D0
      STCHRG(10)  = 0.0231300D0
      STCHRG(11)  = 0.0110570D0
      STCHRG(12)  = 0.0110570D0
      STCHRG(13)  = -0.2250820D0
      STCHRG(14)  = 0.0126350D0
      STCHRG(15)  = -0.0151280D0
      STCHRG(16)  = -0.0151280D0
      STCHRG(17)  = 0.0126350D0
      STCHRG(18)  = -0.2250820D0
      STCHRG(19)  = -0.2250820D0
      STCHRG(20)  = 0.0126350D0
      STCHRG(21)  = -0.0151280D0
      STCHRG(22)  = -0.0151280D0
      STCHRG(23)  = 0.0126350D0
      STCHRG(24)  = -0.2250820D0
      STCHRG(25)  = 0.0861370D0
      STCHRG(26)  = 0.1967670D0
      STCHRG(27)  = 0.1967670D0
      STCHRG(28)  = 0.0861370D0
      STCHRG(29)  = 0.0861370D0
      STCHRG(30)  = 0.1967670D0
      STCHRG(31)  = 0.1967670D0
      STCHRG(32)  = 0.0861370D0
      STCHRG(33)  = -0.1762830D0
      STCHRG(34)  = -0.2579940D0
      STCHRG(35)  = -0.2432570D0
      STCHRG(36)  = -0.2432570D0
      STCHRG(37)  = -0.2579940D0
      STCHRG(38)  = -0.1762830D0
      STCHRG(39)  = -0.1104240D0
      STCHRG(40)  = -0.1762830D0
      STCHRG(41)  = -0.2579940D0
      STCHRG(42)  = -0.2432570D0
      STCHRG(43)  = -0.2432570D0
      STCHRG(44)  = -0.2579940D0
      STCHRG(45)  = -0.1762830D0
      STCHRG(46)  = -0.1104240D0
      STCHRG(47)  = 0.1519720D0
      STCHRG(48)  = 0.1519720D0
      STCHRG(49)  = 0.1519720D0
      STCHRG(50)  = 0.1519720D0
      STCHRG(51)  = 0.1447300D0
      STCHRG(52)  = 0.1447300D0
      STCHRG(53)  = 0.1447300D0
      STCHRG(54)  = 0.1447300D0
      STCHRG(55)  = 0.1529060D0
      STCHRG(56)  = 0.1593090D0
      STCHRG(57)  = 0.1593090D0
      STCHRG(58)  = 0.1529060D0
      STCHRG(59)  = 0.1356020D0
      STCHRG(60)  = 0.1529060D0
      STCHRG(61)  = 0.1593090D0
      STCHRG(62)  = 0.1593090D0
      STCHRG(63)  = 0.1529060D0
      STCHRG(64)  = 0.1356020D0

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO C7
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO C8
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO C9
      RBSTLA(4,:)  = SITE(10,:)  - SITE(4,:)                 ! Z FROM C4 TO C10
      RBSTLA(5,:)  = SITE(11,:)  - SITE(5,:)                 ! Z FROM C5 TO C11
      RBSTLA(6,:)  = SITE(12,:)  - SITE(6,:)                 ! Z FROM C6 TO C12
      RBSTLA(7,:)  = SITE(1,:)  - SITE(7,:)                 ! Z FROM C7 TO C1
      RBSTLA(8,:)  = SITE(2,:)  - SITE(8,:)                 ! Z FROM C8 TO C2
      RBSTLA(9,:)  = SITE(3,:)  - SITE(9,:)                 ! Z FROM C9 TO C3
      RBSTLA(10,:)  = SITE(4,:)  - SITE(10,:)                 ! Z FROM C10 TO C4
      RBSTLA(11,:)  = SITE(5,:)  - SITE(11,:)                 ! Z FROM C11 TO C5
      RBSTLA(12,:)  = SITE(6,:)  - SITE(12,:)                 ! Z FROM C12 TO C6
      RBSTLA(13,:)  = SITE(51,:)  - SITE(13,:)                 ! Z FROM C13 TO H5
      RBSTLA(14,:)  = SITE(25,:)  - SITE(14,:)                 ! Z FROM C14 TO C25
      RBSTLA(15,:)  = SITE(26,:)  - SITE(15,:)                 ! Z FROM C15 TO C26
      RBSTLA(16,:)  = SITE(27,:)  - SITE(16,:)                 ! Z FROM C16 TO C27
      RBSTLA(17,:)  = SITE(28,:)  - SITE(17,:)                 ! Z FROM C17 TO C28
      RBSTLA(18,:)  = SITE(52,:)  - SITE(18,:)                 ! Z FROM C18 TO H6
      RBSTLA(19,:)  = SITE(53,:)  - SITE(19,:)                 ! Z FROM C19 TO H7
      RBSTLA(20,:)  = SITE(29,:)  - SITE(20,:)                 ! Z FROM C20 TO C29
      RBSTLA(21,:)  = SITE(30,:)  - SITE(21,:)                 ! Z FROM C21 TO C30
      RBSTLA(22,:)  = SITE(31,:)  - SITE(22,:)                 ! Z FROM C22 TO C31
      RBSTLA(23,:)  = SITE(32,:)  - SITE(23,:)                 ! Z FROM C23 TO C32
      RBSTLA(24,:)  = SITE(54,:)  - SITE(24,:)                 ! Z FROM C24 TO H8
      RBSTLA(25,:)  = SITE(14,:)  - SITE(25,:)                 ! Z FROM C25 TO C14
      RBSTLA(26,:)  = SITE(15,:)  - SITE(26,:)                 ! Z FROM C26 TO C15
      RBSTLA(27,:)  = SITE(16,:)  - SITE(27,:)                 ! Z FROM C27 TO C16
      RBSTLA(28,:)  = SITE(17,:)  - SITE(28,:)                 ! Z FROM C28 TO C17
      RBSTLA(29,:)  = SITE(20,:)  - SITE(29,:)                 ! Z FROM C29 TO C20
      RBSTLA(30,:)  = SITE(21,:)  - SITE(30,:)                 ! Z FROM C30 TO C21
      RBSTLA(31,:)  = SITE(22,:)  - SITE(31,:)                 ! Z FROM C31 TO C22
      RBSTLA(32,:)  = SITE(23,:)  - SITE(32,:)                 ! Z FROM C32 TO C23
      RBSTLA(33,:)  = SITE(47,:)  - SITE(33,:)                 ! Z FROM C33 TO H1
      RBSTLA(34,:)  = SITE(55,:)  - SITE(34,:)                 ! Z FROM C34 TO H9
      RBSTLA(35,:)  = SITE(56,:)  - SITE(35,:)                 ! Z FROM C35 TO H10
      RBSTLA(36,:)  = SITE(57,:)  - SITE(36,:)                 ! Z FROM C36 TO H11
      RBSTLA(37,:)  = SITE(58,:)  - SITE(37,:)                 ! Z FROM C37 TO H12
      RBSTLA(38,:)  = SITE(48,:)  - SITE(38,:)                 ! Z FROM C38 TO H2
      RBSTLA(39,:)  = SITE(59,:)  - SITE(39,:)                 ! Z FROM C39 TO H13
      RBSTLA(40,:)  = SITE(49,:)  - SITE(40,:)                 ! Z FROM C40 TO H3
      RBSTLA(41,:)  = SITE(60,:)  - SITE(41,:)                 ! Z FROM C41 TO H14
      RBSTLA(42,:)  = SITE(61,:)  - SITE(42,:)                 ! Z FROM C42 TO H15
      RBSTLA(43,:)  = SITE(62,:)  - SITE(43,:)                 ! Z FROM C43 TO H16
      RBSTLA(44,:)  = SITE(63,:)  - SITE(44,:)                 ! Z FROM C44 TO H17
      RBSTLA(45,:)  = SITE(50,:)  - SITE(45,:)                 ! Z FROM C45 TO H4
      RBSTLA(46,:)  = SITE(64,:)  - SITE(46,:)                 ! Z FROM C46 TO H18
      RBSTLA(47,:)  = SITE(47,:)  - SITE(33,:)                 ! Z FROM C33 TO H1
      RBSTLA(48,:)  = SITE(48,:)  - SITE(38,:)                 ! Z FROM C38 TO H2
      RBSTLA(49,:)  = SITE(49,:)  - SITE(40,:)                 ! Z FROM C40 TO H3
      RBSTLA(50,:)  = SITE(50,:)  - SITE(45,:)                 ! Z FROM C45 TO H4
      RBSTLA(51,:)  = SITE(51,:)  - SITE(13,:)                 ! Z FROM C13 TO H5
      RBSTLA(52,:)  = SITE(52,:)  - SITE(18,:)                 ! Z FROM C18 TO H6
      RBSTLA(53,:)  = SITE(53,:)  - SITE(19,:)                 ! Z FROM C19 TO H7
      RBSTLA(54,:)  = SITE(54,:)  - SITE(24,:)                 ! Z FROM C24 TO H8
      RBSTLA(55,:)  = SITE(55,:)  - SITE(34,:)                 ! Z FROM C34 TO H9
      RBSTLA(56,:)  = SITE(56,:)  - SITE(35,:)                 ! Z FROM C35 TO H10
      RBSTLA(57,:)  = SITE(57,:)  - SITE(36,:)                 ! Z FROM C36 TO H11
      RBSTLA(58,:)  = SITE(58,:)  - SITE(37,:)                 ! Z FROM C37 TO H12
      RBSTLA(59,:)  = SITE(59,:)  - SITE(39,:)                 ! Z FROM C39 TO H13
      RBSTLA(60,:)  = SITE(60,:)  - SITE(41,:)                 ! Z FROM C41 TO H14
      RBSTLA(61,:)  = SITE(61,:)  - SITE(42,:)                 ! Z FROM C42 TO H15
      RBSTLA(62,:)  = SITE(62,:)  - SITE(43,:)                 ! Z FROM C43 TO H16
      RBSTLA(63,:)  = SITE(63,:)  - SITE(44,:)                 ! Z FROM C44 TO H17
      RBSTLA(64,:)  = SITE(64,:)  - SITE(46,:)                 ! Z FROM C46 TO H18

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFOCTABENZOCORONENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCIRCUMCORONENE()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C54H18

!     UNITS INITIALLY IN ANGSTROM

      SITE(1,:)  = (/0.000000D0,  -1.419570D0,  0.000000D0/)   !C1
      SITE(2,:)  = (/1.229390D0,  -0.709790D0,  0.000000D0/)   !C2
      SITE(3,:)  = (/1.229390D0,  0.709790D0,  0.000000D0/)   !C3
      SITE(4,:)  = (/0.000000D0,  1.419570D0,  0.000000D0/)   !C4
      SITE(5,:)  = (/-1.229390D0,  0.709790D0,  0.000000D0/)   !C5
      SITE(6,:)  = (/-1.229390D0,  -0.709790D0,  0.000000D0/)   !C6
      SITE(7,:)  = (/0.000000D0,  -2.848790D0,  0.000000D0/)   !C7
      SITE(8,:)  = (/2.467130D0,  -1.424390D0,  0.000000D0/)   !C8
      SITE(9,:)  = (/2.467130D0,  1.424390D0,  0.000000D0/)   !C9
      SITE(10,:)  = (/0.000000D0,  2.848790D0,  0.000000D0/)   !C10
      SITE(11,:)  = (/-2.467130D0,  1.424390D0,  0.000000D0/)   !C11
      SITE(12,:)  = (/-2.467130D0,  -1.424390D0,  0.000000D0/)   !C12
      SITE(13,:)  = (/1.225320D0,  -4.988360D0,  0.000000D0/)   !C13
      SITE(14,:)  = (/2.465830D0,  -2.844750D0,  0.000000D0/)   !C14
      SITE(15,:)  = (/3.696550D0,  -0.713100D0,  0.000000D0/)   !C15
      SITE(16,:)  = (/3.696550D0,  0.713100D0,  0.000000D0/)   !C16
      SITE(17,:)  = (/2.465830D0,  2.844750D0,  0.000000D0/)   !C17
      SITE(18,:)  = (/1.225320D0,  4.988360D0,  0.000000D0/)   !C18
      SITE(19,:)  = (/-1.225320D0,  4.988360D0,  0.000000D0/)   !C19
      SITE(20,:)  = (/-2.465830D0,  2.844750D0,  0.000000D0/)   !C20
      SITE(21,:)  = (/-3.696550D0,  0.713100D0,  0.000000D0/)   !C21
      SITE(22,:)  = (/-3.696550D0,  -0.713100D0,  0.000000D0/)   !C22
      SITE(23,:)  = (/-2.465830D0,  -2.844750D0,  0.000000D0/)   !C23
      SITE(24,:)  = (/-1.225320D0,  -4.988360D0,  0.000000D0/)   !C24
      SITE(25,:)  = (/1.230710D0,  -3.557840D0,  0.000000D0/)   !C25
      SITE(26,:)  = (/3.669170D0,  -4.992040D0,  0.000000D0/)   !C26
      SITE(27,:)  = (/4.932720D0,  -1.433020D0,  0.000000D0/)   !C27
      SITE(28,:)  = (/4.932720D0,  1.433020D0,  0.000000D0/)   !C28
      SITE(29,:)  = (/3.669170D0,  4.992040D0,  0.000000D0/)   !C29
      SITE(30,:)  = (/1.230710D0,  3.557840D0,  0.000000D0/)   !C30
      SITE(31,:)  = (/-1.230710D0,  3.557840D0,  0.000000D0/)   !C31
      SITE(32,:)  = (/-3.669170D0,  4.992040D0,  0.000000D0/)   !C32
      SITE(33,:)  = (/-4.932720D0,  1.433020D0,  0.000000D0/)   !C33
      SITE(34,:)  = (/-4.932720D0,  -1.433020D0,  0.000000D0/)   !C34
      SITE(35,:)  = (/-3.669170D0,  -4.992040D0,  0.000000D0/)   !C35
      SITE(36,:)  = (/-1.230710D0,  -3.557840D0,  0.000000D0/)   !C36
      SITE(37,:)  = (/2.488640D0,  -5.673610D0,  0.000000D0/)   !C37
      SITE(38,:)  = (/3.707390D0,  -3.555340D0,  0.000000D0/)   !C38
      SITE(39,:)  = (/4.908860D0,  -2.834130D0,  0.000000D0/)   !C39
      SITE(40,:)  = (/6.157820D0,  -0.681580D0,  0.000000D0/)   !C40
      SITE(41,:)  = (/6.157820D0,  0.681580D0,  0.000000D0/)   !C41
      SITE(42,:)  = (/4.908860D0,  2.834130D0,  0.000000D0/)   !C42
      SITE(43,:)  = (/3.707390D0,  3.555340D0,  0.000000D0/)   !C43
      SITE(44,:)  = (/2.488640D0,  5.673610D0,  0.000000D0/)   !C44
      SITE(45,:)  = (/0.000000D0,  5.668240D0,  0.000000D0/)   !C45
      SITE(46,:)  = (/-2.488640D0,  5.673610D0,  0.000000D0/)   !C46
      SITE(47,:)  = (/-3.707390D0,  3.555340D0,  0.000000D0/)   !C47
      SITE(48,:)  = (/-4.908860D0,  2.834130D0,  0.000000D0/)   !C48
      SITE(49,:)  = (/-6.157820D0,  0.681580D0,  0.000000D0/)   !C49
      SITE(50,:)  = (/-6.157820D0,  -0.681580D0,  0.000000D0/)   !C50
      SITE(51,:)  = (/-4.908860D0,  -2.834130D0,  0.000000D0/)   !C51
      SITE(52,:)  = (/-3.707390D0,  -3.555340D0,  0.000000D0/)   !C52
      SITE(53,:)  = (/-2.488640D0,  -5.673610D0,  0.000000D0/)   !C53
      SITE(54,:)  = (/0.000000D0,  -5.668240D0,  0.000000D0/)   !C54
      SITE(55,:)  = (/2.484360D0,  -6.760930D0,  0.000000D0/)   !H1
      SITE(56,:)  = (/4.612950D0,  -5.532010D0,  0.000000D0/)   !H2
      SITE(57,:)  = (/5.851170D0,  -3.378180D0,  0.000000D0/)   !H3
      SITE(58,:)  = (/7.097340D0,  -1.228930D0,  0.000000D0/)   !H4
      SITE(59,:)  = (/7.097340D0,  1.228930D0,  0.000000D0/)   !H5
      SITE(60,:)  = (/5.851170D0,  3.378180D0,  0.000000D0/)   !H6
      SITE(61,:)  = (/4.612950D0,  5.532010D0,  0.000000D0/)   !H7
      SITE(62,:)  = (/2.484360D0,  6.760930D0,  0.000000D0/)   !H8
      SITE(63,:)  = (/0.000000D0,  6.756340D0,  0.000000D0/)   !H9
      SITE(64,:)  = (/-2.484360D0,  6.760930D0,  0.000000D0/)   !H10
      SITE(65,:)  = (/-4.612950D0,  5.532010D0,  0.000000D0/)   !H11
      SITE(66,:)  = (/-5.851170D0,  3.378180D0,  0.000000D0/)   !H12
      SITE(67,:)  = (/-7.097340D0,  1.228930D0,  0.000000D0/)   !H13
      SITE(68,:)  = (/-7.097340D0,  -1.228930D0,  0.000000D0/)   !H14
      SITE(69,:)  = (/-5.851170D0,  -3.378180D0,  0.000000D0/)   !H15
      SITE(70,:)  = (/-4.612950D0,  -5.532010D0,  0.000000D0/)   !H16
      SITE(71,:)  = (/-2.484360D0,  -6.760930D0,  0.000000D0/)   !H17
      SITE(72,:)  = (/0.000000D0,  -6.756340D0,  0.000000D0/)   !H18

      SITE(:,:) =  SITE(:,:)/0.5291770D0

      STCHRG(1)  = 0.0029610D0
      STCHRG(2)  = 0.0029610D0
      STCHRG(3)  = 0.0029610D0
      STCHRG(4)  = 0.0029610D0
      STCHRG(5)  = 0.0029610D0
      STCHRG(6)  = 0.0029610D0
      STCHRG(7)  = 0.0008790D0
      STCHRG(8)  = 0.0008790D0
      STCHRG(9)  = 0.0008790D0
      STCHRG(10)  = 0.0008790D0
      STCHRG(11)  = 0.0008790D0
      STCHRG(12)  = 0.0008790D0
      STCHRG(13)  = 0.2298540D0
      STCHRG(14)  = -0.0308620D0
      STCHRG(15)  = -0.0308620D0
      STCHRG(16)  = -0.0308620D0
      STCHRG(17)  = -0.0308620D0
      STCHRG(18)  = 0.2298540D0
      STCHRG(19)  = 0.2298540D0
      STCHRG(20)  = -0.0308620D0
      STCHRG(21)  = -0.0308620D0
      STCHRG(22)  = -0.0308620D0
      STCHRG(23)  = -0.0308620D0
      STCHRG(24)  = 0.2298540D0
      STCHRG(25)  = -0.0308620D0
      STCHRG(26)  = -0.2404810D0
      STCHRG(27)  = 0.2298540D0
      STCHRG(28)  = 0.2298540D0
      STCHRG(29)  = -0.2404810D0
      STCHRG(30)  = -0.0308620D0
      STCHRG(31)  = -0.0308620D0
      STCHRG(32)  = -0.2404810D0
      STCHRG(33)  = 0.2298540D0
      STCHRG(34)  = 0.2298540D0
      STCHRG(35)  = -0.2404810D0
      STCHRG(36)  = -0.0308620D0
      STCHRG(37)  = -0.2404810D0
      STCHRG(38)  = 0.2298540D0
      STCHRG(39)  = -0.4184350D0
      STCHRG(40)  = -0.2404810D0
      STCHRG(41)  = -0.2404810D0
      STCHRG(42)  = -0.4184350D0
      STCHRG(43)  = 0.2298540D0
      STCHRG(44)  = -0.2404810D0
      STCHRG(45)  = -0.4184350D0
      STCHRG(46)  = -0.2404810D0
      STCHRG(47)  = 0.2298540D0
      STCHRG(48)  = -0.4184350D0
      STCHRG(49)  = -0.2404810D0
      STCHRG(50)  = -0.2404810D0
      STCHRG(51)  = -0.4184350D0
      STCHRG(52)  = 0.2298540D0
      STCHRG(53)  = -0.2404810D0
      STCHRG(54)  = -0.4184350D0
      STCHRG(55)  = 0.1585010D0
      STCHRG(56)  = 0.1585010D0
      STCHRG(57)  = 0.1805690D0
      STCHRG(58)  = 0.1585010D0
      STCHRG(59)  = 0.1585010D0
      STCHRG(60)  = 0.1805690D0
      STCHRG(61)  = 0.1585010D0
      STCHRG(62)  = 0.1585010D0
      STCHRG(63)  = 0.1805690D0
      STCHRG(64)  = 0.1585010D0
      STCHRG(65)  = 0.1585010D0
      STCHRG(66)  = 0.1805690D0
      STCHRG(67)  = 0.1585010D0
      STCHRG(68)  = 0.1585010D0
      STCHRG(69)  = 0.1805690D0
      STCHRG(70)  = 0.1585010D0
      STCHRG(71)  = 0.1585010D0
      STCHRG(72)  = 0.1805690D0

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO C7
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO C8
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO C9
      RBSTLA(4,:)  = SITE(10,:)  - SITE(4,:)                 ! Z FROM C4 TO C10
      RBSTLA(5,:)  = SITE(11,:)  - SITE(5,:)                 ! Z FROM C5 TO C11
      RBSTLA(6,:)  = SITE(12,:)  - SITE(6,:)                 ! Z FROM C6 TO C12
      RBSTLA(7,:)  = SITE(1,:)  - SITE(7,:)                 ! Z FROM C7 TO C1
      RBSTLA(8,:)  = SITE(2,:)  - SITE(8,:)                 ! Z FROM C8 TO C2
      RBSTLA(9,:)  = SITE(3,:)  - SITE(9,:)                 ! Z FROM C9 TO C3
      RBSTLA(10,:)  = SITE(4,:)  - SITE(10,:)                 ! Z FROM C10 TO C4
      RBSTLA(11,:)  = SITE(5,:)  - SITE(11,:)                 ! Z FROM C11 TO C5
      RBSTLA(12,:)  = SITE(6,:)  - SITE(12,:)                 ! Z FROM C12 TO C6
      RBSTLA(13,:)  = SITE(25,:)  - SITE(13,:)                 ! Z FROM C13 TO C25
      RBSTLA(14,:)  = SITE(38,:)  - SITE(14,:)                 ! Z FROM C14 TO C38
      RBSTLA(15,:)  = SITE(27,:)  - SITE(15,:)                 ! Z FROM C15 TO C27
      RBSTLA(16,:)  = SITE(28,:)  - SITE(16,:)                 ! Z FROM C16 TO C28
      RBSTLA(17,:)  = SITE(43,:)  - SITE(17,:)                 ! Z FROM C17 TO C43
      RBSTLA(18,:)  = SITE(30,:)  - SITE(18,:)                 ! Z FROM C18 TO C30
      RBSTLA(19,:)  = SITE(31,:)  - SITE(19,:)                 ! Z FROM C19 TO C31
      RBSTLA(20,:)  = SITE(47,:)  - SITE(20,:)                 ! Z FROM C20 TO C47
      RBSTLA(21,:)  = SITE(33,:)  - SITE(21,:)                 ! Z FROM C21 TO C33
      RBSTLA(22,:)  = SITE(34,:)  - SITE(22,:)                 ! Z FROM C22 TO C34
      RBSTLA(23,:)  = SITE(52,:)  - SITE(23,:)                 ! Z FROM C23 TO C52
      RBSTLA(24,:)  = SITE(36,:)  - SITE(24,:)                 ! Z FROM C24 TO C36
      RBSTLA(25,:)  = SITE(13,:)  - SITE(25,:)                 ! Z FROM C25 TO C13
      RBSTLA(26,:)  = SITE(56,:)  - SITE(26,:)                 ! Z FROM C26 TO H2
      RBSTLA(27,:)  = SITE(15,:)  - SITE(27,:)                 ! Z FROM C27 TO C15
      RBSTLA(28,:)  = SITE(16,:)  - SITE(28,:)                 ! Z FROM C28 TO C16
      RBSTLA(29,:)  = SITE(61,:)  - SITE(29,:)                 ! Z FROM C29 TO H7
      RBSTLA(30,:)  = SITE(18,:)  - SITE(30,:)                 ! Z FROM C30 TO C18
      RBSTLA(31,:)  = SITE(19,:)  - SITE(31,:)                 ! Z FROM C31 TO C19
      RBSTLA(32,:)  = SITE(65,:)  - SITE(32,:)                 ! Z FROM C32 TO H11
      RBSTLA(33,:)  = SITE(21,:)  - SITE(33,:)                 ! Z FROM C33 TO C21
      RBSTLA(34,:)  = SITE(22,:)  - SITE(34,:)                 ! Z FROM C34 TO C22
      RBSTLA(35,:)  = SITE(70,:)  - SITE(35,:)                 ! Z FROM C35 TO H16
      RBSTLA(36,:)  = SITE(24,:)  - SITE(36,:)                 ! Z FROM C36 TO C24
      RBSTLA(37,:)  = SITE(55,:)  - SITE(37,:)                 ! Z FROM C37 TO H1
      RBSTLA(38,:)  = SITE(14,:)  - SITE(38,:)                 ! Z FROM C38 TO C14
      RBSTLA(39,:)  = SITE(57,:)  - SITE(39,:)                 ! Z FROM C39 TO H3
      RBSTLA(40,:)  = SITE(58,:)  - SITE(40,:)                 ! Z FROM C40 TO H4
      RBSTLA(41,:)  = SITE(59,:)  - SITE(41,:)                 ! Z FROM C41 TO H5
      RBSTLA(42,:)  = SITE(60,:)  - SITE(42,:)                 ! Z FROM C42 TO H6
      RBSTLA(43,:)  = SITE(17,:)  - SITE(43,:)                 ! Z FROM C43 TO C17
      RBSTLA(44,:)  = SITE(62,:)  - SITE(44,:)                 ! Z FROM C44 TO H8
      RBSTLA(45,:)  = SITE(63,:)  - SITE(45,:)                 ! Z FROM C45 TO H9
      RBSTLA(46,:)  = SITE(64,:)  - SITE(46,:)                 ! Z FROM C46 TO H10
      RBSTLA(47,:)  = SITE(20,:)  - SITE(47,:)                 ! Z FROM C47 TO C20
      RBSTLA(48,:)  = SITE(66,:)  - SITE(48,:)                 ! Z FROM C48 TO H12
      RBSTLA(49,:)  = SITE(67,:)  - SITE(49,:)                 ! Z FROM C49 TO H13
      RBSTLA(50,:)  = SITE(68,:)  - SITE(50,:)                 ! Z FROM C50 TO H14
      RBSTLA(51,:)  = SITE(69,:)  - SITE(51,:)                 ! Z FROM C51 TO H15
      RBSTLA(52,:)  = SITE(23,:)  - SITE(52,:)                 ! Z FROM C52 TO C23
      RBSTLA(53,:)  = SITE(71,:)  - SITE(53,:)                 ! Z FROM C53 TO H17
      RBSTLA(54,:)  = SITE(72,:)  - SITE(54,:)                 ! Z FROM C54 TO H18
      RBSTLA(55,:)  = SITE(55,:)  - SITE(37,:)                 ! Z FROM C37 TO H1
      RBSTLA(56,:)  = SITE(56,:)  - SITE(26,:)                 ! Z FROM C26 TO H2
      RBSTLA(57,:)  = SITE(57,:)  - SITE(39,:)                 ! Z FROM C39 TO H3
      RBSTLA(58,:)  = SITE(58,:)  - SITE(40,:)                 ! Z FROM C40 TO H4
      RBSTLA(59,:)  = SITE(59,:)  - SITE(41,:)                 ! Z FROM C41 TO H5
      RBSTLA(60,:)  = SITE(60,:)  - SITE(42,:)                 ! Z FROM C42 TO H6
      RBSTLA(61,:)  = SITE(61,:)  - SITE(29,:)                 ! Z FROM C29 TO H7
      RBSTLA(62,:)  = SITE(62,:)  - SITE(44,:)                 ! Z FROM C44 TO H8
      RBSTLA(63,:)  = SITE(63,:)  - SITE(45,:)                 ! Z FROM C45 TO H9
      RBSTLA(64,:)  = SITE(64,:)  - SITE(46,:)                 ! Z FROM C46 TO H10
      RBSTLA(65,:)  = SITE(65,:)  - SITE(32,:)                 ! Z FROM C32 TO H11
      RBSTLA(66,:)  = SITE(66,:)  - SITE(48,:)                 ! Z FROM C48 TO H12
      RBSTLA(67,:)  = SITE(67,:)  - SITE(49,:)                 ! Z FROM C49 TO H13
      RBSTLA(68,:)  = SITE(68,:)  - SITE(50,:)                 ! Z FROM C50 TO H14
      RBSTLA(69,:)  = SITE(69,:)  - SITE(51,:)                 ! Z FROM C51 TO H15
      RBSTLA(70,:)  = SITE(70,:)  - SITE(35,:)                 ! Z FROM C35 TO H16
      RBSTLA(71,:)  = SITE(71,:)  - SITE(53,:)                 ! Z FROM C53 TO H17
      RBSTLA(72,:)  = SITE(72,:)  - SITE(54,:)                 ! Z FROM C54 TO H18

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFCIRCUMCORONENE

!     ----------------------------------------------------------------------------------------------
