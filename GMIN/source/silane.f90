      SUBROUTINE SILANE (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, EPS11, EPS22, EPS12, REQ11, REQ22, REQ12, MRHO11, MRHO22, MRHO12 

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, DSS, FCTR, RHO, REQ, EPS
      DOUBLE PRECISION :: RSS(3), NR(3), P(3)
      DOUBLE PRECISION :: R(NATOMS*NRBSITES/2,3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      LOGICAL          :: GTEST

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = X(J3-2:J3) + MATMUL(RMI(:,:),SITE(J2,:))
!            PRINT *, R(J4,:)

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS - 1 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO I = 1, NRBSITES

            J7    = NRBSITES*(J1-1) + I

            DO J2 = J1 + 1, REALNATOMS

               J4 = 3*J2
               J6 = OFFSET + J4

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  DSS    = DSQRT(DOT_PRODUCT(RSS(:),RSS(:)))

                  IF ((I == 1) .AND. (J == 1))  THEN
                     RHO = MRHO11; REQ = REQ11; EPS = EPS11
                  ELSEIF ((I > 1) .AND. (J > 1)) THEN
                     RHO = MRHO22; REQ = REQ22; EPS = EPS22
                  ELSE
                     RHO = MRHO12; REQ = REQ12; EPS = EPS12
                  ENDIF
                  FCTR   = EXP(RHO*(REQ - DSS)) 
                  ENERGY = ENERGY + EPS*((1.D0 - FCTR)*(1.D0 - FCTR) - 1.D0)

!     DVDR = DVDR/R
                  IF (GTEST) THEN

                     DVDR       = EPS*2.D0*RHO*(-EXP(2.D0*RHO*(REQ - DSS)) + FCTR)/DSS
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

      END SUBROUTINE SILANE
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFSILANE()

      USE COMMONS

      IMPLICIT NONE
      DOUBLE PRECISION :: FCTR
!     SiH4 : Si - 1, H - 2

      FCTR      = 0.85391048D0
      SITE(1,:) = (/ 0.0D0, 0.0D0, 0.0D0/)
      SITE(2,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      SITE(3,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      SITE(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)
      SITE(5,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)

      EPS11  = 1.06D-6 
      EPS12  = 0.042D0         ! kcal/mol
      EPS22  = 0.085D0
      MRHO11 = 1.352d0
      MRHO12 = 8.190d0         ! A-1
      MRHO22 = 1.216d0
      REQ11  = 8.743d0
      REQ12  = 2.836d0         ! A
      REQ22  = 3.348d0

      END SUBROUTINE DEFSILANE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWSILANE()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(NRBSITES*3), P3(3,3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='silane.xyz', STATUS='UNKNOWN')

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

               IF (J2 == 1) THEN
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'Si', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE                                  
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWSILANE
