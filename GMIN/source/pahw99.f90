      SUBROUTINE PAHW99(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, PAHID, SITE, NCARBON

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, NCPHST
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), CHARGE(NRBSITES)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJSS(3), PI(3)
      DOUBLE PRECISION :: DSS2, R2, R6, ABSR, DVDR, ENERGY, EXPFCT
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: CWRCC, CWRHH, CWRCH, CWACC, CWAHH, CWACH, CWECC, CWEHH, CWECH, CCKJ
      LOGICAL          :: GTEST

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      IF (PAHID == 1) THEN
          CALL DEFBNZNW99(CHARGE, CWRCC, CWRHH, CWRCH, CWACC, CWAHH, CWACH, CWECC, CWEHH, CWECH, CCKJ)
          NCARBON = 6
      ENDIF

!     Number carbon plus hydrogen sites

      NCPHST = NCARBON + (NRBSITES - NCARBON)/2

      ENERGY    = 0.D0
      G(:)      = 0.D0

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

      DO J1 = 1, REALNATOMS - 1 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     Sum over sites to calculate repulsion and dispersion interactions

            DO I = 1, NCPHST

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NCPHST

                  J8       = NRBSITES*(J2-1) + J
                  RIJSS(:) = R(J7,:) - R(J8,:)
                  DSS2     = DOT_PRODUCT(RIJSS,RIJSS)
                  R2       = 1.D0/DSS2
                  R6       = R2*R2*R2
                  ABSR     = SQRT(DSS2)
                 
                  IF (I <= NCARBON .AND. J <= NCARBON) THEN
                     EXPFCT = EXP(-CWECC*ABSR)
                     ENERGY = ENERGY + CWRCC*EXPFCT - CWACC*R6 
                  ELSEIF (I > NCARBON .AND. J > NCARBON) THEN
                     EXPFCT = EXP(-CWEHH*ABSR)
                     ENERGY = ENERGY + CWRHH*EXPFCT - CWAHH*R6
                  ELSE
                     EXPFCT = EXP(-CWECH*ABSR)
                     ENERGY = ENERGY + CWRCH*EXPFCT - CWACH*R6
                  ENDIF 

                  ABSR     = SQRT(DSS2)
                  
                  IF (GTEST) THEN
!     DVDR = DVDR/R
                     IF (I <= NCARBON .AND. J <= NCARBON) THEN
                        DVDR   = - CWRCC*CWECC*EXPFCT/ABSR + 6.D0*CWACC*R6*R2
                     ELSEIF (I > NCARBON .AND. J > NCARBON) THEN
                        DVDR   = - CWRHH*CWEHH*EXPFCT/ABSR + 6.D0*CWAHH*R6*R2
                     ELSE
                        DVDR   = - CWRCH*CWECH*EXPFCT/ABSR + 6.D0*CWACH*R6*R2
                     ENDIF

                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RIJSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RIJSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RIJSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RIJSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RIJSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RIJSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

!     Sum over sites to calculate coulomb interactions

            DO I = 1, NRBSITES 

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  IF ((I > NCARBON .AND. I < NCPHST) .OR. (J > NCARBON .AND. J < NCPHST)) CYCLE  

                  J8       = NRBSITES*(J2-1) + J
                  RIJSS(:) = R(J7,:) - R(J8,:)
                  DSS2     = DOT_PRODUCT(RIJSS,RIJSS)
                  R2       = 1.D0/DSS2
                  R6       = R2*R2*R2
                  ABSR     = SQRT(DSS2)
                  ENERGY   = ENERGY + CCKJ*CHARGE(I)*CHARGE(J)/ABSR

                  IF (GTEST) THEN

                     DVDR       =-CCKJ*CHARGE(I)*CHARGE(J)*R2/ABSR
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RIJSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RIJSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RIJSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RIJSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RIJSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RIJSS,DR3(J8,:))

                  ENDIF
              
               ENDDO

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE PAHW99 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBNZNW99(CHARGE, CWRCC, CWRHH, CWRCH, CWACC, CWAHH, CWACH, CWECC, CWEHH, CWECH, CCKJ)

      USE COMMONS, ONLY: NRBSITES, SITE

      IMPLICIT NONE

      DOUBLE PRECISION :: CHARGE(NRBSITES), CWRCC, CWRHH, CWRCH, CWACC, CWAHH, CWACH, CWECC, CWEHH, CWECH, CCKJ

!     C6H6
!     in angstrom; C: 1-6, H: 7-12(foreshortened), H: 13-18 (only for Coulomb)

      SITE(1,:)  = (/ 1.39662264416316D0, 0.00000000000000D0, 0.0D0/)
      SITE(2,:)  = (/ 0.69831132208158D0, 1.20951068934589D0, 0.0D0/)
      SITE(3,:)  = (/-0.69831132208158D0, 1.20951068934589D0, 0.0D0/)
      SITE(4,:)  = (/-1.39662264416316D0, 0.00000000000000D0, 0.0D0/)
      SITE(5,:)  = (/-0.69831132208158D0,-1.20951068934589D0, 0.0D0/)
      SITE(6,:)  = (/ 0.69831132208158D0,-1.20951068934589D0, 0.0D0/)
      SITE(7,:)  = (/ 2.38363492050584D0, 0.00000000000000D0, 0.0D0/)
      SITE(8,:)  = (/ 1.19181746025292D0, 2.06428839450576D0, 0.0D0/)
      SITE(9,:)  = (/-1.19181746025292D0, 2.06428839450576D0, 0.0D0/)
      SITE(10,:) = (/-2.38363492050584D0, 0.00000000000000D0, 0.0D0/)
      SITE(11,:) = (/-1.19181746025292D0,-2.06428839450576D0, 0.0D0/)
      SITE(12,:) = (/ 1.19181746025292D0,-2.06428839450576D0, 0.0D0/)
      SITE(13,:) = (/ 2.48363492050584D0, 0.00000000000000D0, 0.0D0/)
      SITE(14,:) = (/ 1.24181746025292D0, 2.15089093488420D0, 0.0D0/)
      SITE(15,:) = (/-1.24181746025292D0, 2.15089093488420D0, 0.0D0/)
      SITE(16,:) = (/-2.48363492050584D0, 0.00000000000000D0, 0.0D0/)
      SITE(17,:) = (/-1.24181746025292d0,-2.15089093488420d0, 0.0D0/)
      SITE(18,:) = (/ 1.24181746025292D0,-2.15089093488420D0, 0.0D0/)

      CWRCC  = 270363.D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      CWRHH  = 12680.D0
      CWRCH  = SQRT(CWRCC*CWRHH)
      CWACC  = 1701.73D0
      CWAHH  = 278.37D0
      CWACH  = SQRT(CWACC*CWAHH)
      CWECC  = 3.60D0
      CWEHH  = 3.56D0
      CWECH  = 0.5D0*(CWECC + CWEHH)
      CCKJ   = 1389.354848D0

      CHARGE(1:6)   =-0.11114D0
      CHARGE(7:12)  = 0.0D0
      CHARGE(13:18) = 0.11114D0

      END SUBROUTINE DEFBNZNW99
