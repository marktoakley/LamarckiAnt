      SUBROUTINE TETRAHEDRA  (X, G, ENERGY, GTEST)

!     each rigidbody consistes of 4 Morse sites in a regular tetrahedral arrangement plus a repulsive site at the centre

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, RHO, MREQ, EPSR

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), R(NATOMS*NRBSITES/2,3), RSS(3), P(3)
      DOUBLE PRECISION :: DR1(NATOMS*NRBSITES/2,3), DR2(NATOMS*NRBSITES/2,3), DR3(NATOMS*NRBSITES/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ENERGY, DVDR, DSS, FCTR, R2, R6, R12
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

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS  

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

            DO I = 1, NRBSITES

               J7   = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  DSS    = DSQRT(DOT_PRODUCT(RSS(:),RSS(:)))
                  FCTR   = EXP(RHO*(1.D0 - DSS/MREQ)) 
                  ENERGY = ENERGY + (1.D0 - FCTR)*(1.D0 - FCTR) - 1.D0
!     DVDR = DVDR/R
                  IF (GTEST) THEN

                     DVDR       = 2.D0*RHO*(-EXP(2.D0*RHO*(1.D0 - DSS/MREQ)) + FCTR)/(MREQ*DSS)
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

            RSS(:) = X(J3-2:J3) - X(J4-2:J4)
            R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
            R6     = R2*R2*R2
            R12    = R6**2
            ENERGY = ENERGY + EPSR*R12

            IF (GTEST) THEN

!     DVDR = DVDR/R
               DVDR       = -12.D0*EPSR*R12*R2
               G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:)
 
            ENDIF

         ENDDO

      ENDDO

!      G(1:3*NATOMS/2) = 0

      END SUBROUTINE TETRAHEDRA
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTDHD()

      USE COMMONS, ONLY: SITE

      IMPLICIT NONE

      DOUBLE PRECISION :: FCTR

      FCTR      = 1.D0/DSQRT(8.D0)
      SITE(1,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      SITE(2,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      SITE(3,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)
      SITE(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)

      END SUBROUTINE DEFTDHD

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWTDHD()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J7 
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(NRBSITES*3), P3(3,3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='tetrahedra.xyz', STATUS='UNKNOWN')

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

               RBCOORDS(3*J2-2:3*J2) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))

               IF (J2 > 1) P3(J2-1,:) = RBCOORDS(3*J2-2:3*J2) - RBCOORDS(1:3)

            ENDDO

            WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10,2X,A12,2X,3F20.10,2X,A12,2X,3F20.10)')          &
     &      'LA', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3), 'atom_vector', P3(1,1), P3(1,2), P3(1,3),  &
     &      'atom_vector', P3(2,1), P3(2,2), P3(2,3), 'atom_vector', P3(3,1), P3(3,2), P3(3,3)

             DO J2 = 2, NRBSITES

                  J4 = J2 + 1
                  IF (J2 == NRBSITES) J4 = 2
                  P(:) = RBCOORDS(3*J4-2:3*J4) - RBCOORDS(3*J2-2:3*J2)
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')                                        &
     &            'LA', RBCOORDS(3*J2-2), RBCOORDS(3*J2-1), RBCOORDS(3*J2), 'atom_vector', P(1), P(2), P(3)

             ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWTDHD
