      SUBROUTINE PTSTST(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, PAPCD, PAPEPS, YKAPPA 

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, JL, JU, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, RIJSQ, R2, RSS(3), RSSSQ, FCTR, ABSRIJ, EXPFCT, KAPPA
      DOUBLE PRECISION :: RIJ(3), P(3)
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      LOGICAL          :: GTEST
      DOUBLE PRECISION :: RLJN, R2LJN, LJSIG2, LJN
      DOUBLE PRECISION :: EPS, SIGMA, SIGSQ, ALPHA     ! potential parameters

      EPS = PAPEPS; SIGMA = PAPCD; KAPPA = YKAPPA 
      LJSIG2 = 1.D0; LJN = 23.D0
      SIGSQ = SIGMA*SIGMA

      OFFSET  = 3*NATOMS/2

      ENERGY  = 0.D0
      IF (GTEST) THEN
        G(:) = 0.D0
      ENDIF

      DO J1 = 1, NATOMS/2

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = X(J3-2:J3) + MATMUL(RMI(:,:),SITE(J2,:))
            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

            ENDIF

         ENDDO

 
      ENDDO

      DO J1 = 1, NATOMS/2 - 1  

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, NATOMS/2

            J4 = 3*J2
            J6 = OFFSET + J4

            RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            RLJN   = (LJSIG2*R2)**(LJN/2.D0)
            R2LJN  = RLJN**2
            ENERGY = ENERGY + 4.D0*(R2LJN-RLJN)   

!            ABSRIJ = DSQRT(RIJSQ)
!            EXPFCT = EXP(-KAPPA*(ABSRIJ - 1.D0))
!            ENERGY = EXPFCT/ABSRIJ

            IF (GTEST) THEN
               DVDR = 4.D0*LJN*R2*(RLJN - 2.D0*R2LJN) 
!               DVDR   =-EXPFCT/RIJSQ*(KAPPA + 1.D0/ABSRIJ)
               G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJ(:)
            ENDIF
           
            DO I = 1, NRBSITES
               J7   = NRBSITES*(J1-1) + I
               IF (I <= NRBSITES/2) THEN
                  JL = NRBSITES/2 + 1
                  JU = NRBSITES
               ELSE
                  JL = 1
                  JU = NRBSITES/2
               ENDIF
               DO J = JL, JU !1, NRBSITES
                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  RSSSQ  = DOT_PRODUCT(RSS(:),RSS(:))
                  FCTR   = EXP(-RSSSQ/(2.D0*SIGSQ))
                  ENERGY = ENERGY - EPS*FCTR
!     DVDR = DVDR/R
                  IF (GTEST) THEN
                     DVDR       = EPS*FCTR/SIGSQ
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

      END SUBROUTINE PTSTST

!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE DEFPTSTST()

      USE COMMONS, ONLY: SITE, NRBSITES

      IMPLICIT NONE

      NRBSITES = 2
      ALLOCATE(SITE(NRBSITES,3))

      SITE(1,:)= 0.5D0*(/ -0.096225D0, 0.301232D0, 0.948682D0/)
      SITE(2,:)= 0.5D0*(/ -0.096225D0,-0.301232D0,-0.948682D0/)

!      SITE(1,:)= 0.5D0/SQRT(3.D0)*(/  1.D0,  1.D0,  1.D0/)
!      SITE(2,:)= 0.5D0/SQRT(3.D0)*(/ -1.D0, -1.D0,  1.D0/)
!      SITE(3,:)= 0.5D0/SQRT(3.D0)*(/ -1.D0,  1.D0, -1.D0/)
!      SITE(4,:)= 0.5D0/SQRT(3.D0)*(/  1.D0, -1.D0, -1.D0/)

      END SUBROUTINE DEFPTSTST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWPTSTST()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(NRBSITES*3), P3(3,3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='ptstst.xyz', STATUS='UNKNOWN')

      GTEST = .FALSE.

      DO J1 = 1, NSAVE

         WRITE(26,'(I6)') (NATOMS/2)*(NRBSITES+1)
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, NATOMS/2

            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)


            WRITE(26,'(A4,3F20.10)') 'O', QMINP(J1,J5-2), QMINP(J1,J5-1), QMINP(J1,J5)
            DO J2 = 1, NRBSITES
               RBCOORDS(3*J2-2:3*J2) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (J2 <= NRBSITES/2) THEN
               WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(3*J2-2), RBCOORDS(3*J2-1), RBCOORDS(3*J2)
               ELSE
               WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(3*J2-2), RBCOORDS(3*J2-1), RBCOORDS(3*J2)
               ENDIF
            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWPTSTST

