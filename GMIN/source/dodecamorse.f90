      SUBROUTINE DODECAMORSE  (X, G, ENERGY, GTEST)
!
!     modified tetrahedra.f90
!     morse sites on the vertices of a dodecahedron
!
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, RHO, MREQ, BEPS, DDMCUT

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

            RSS(:) = X(J4-2:J4) - X(J3-2:J3)
            DSS    = DSQRT(DOT_PRODUCT(RSS(:),RSS(:)))

            IF (DSS.LT.(DDMCUT+1.D0)) THEN

               DO I = 1, NRBSITES
                  J7   = NRBSITES*(J1-1) + I

                  DO J = 1, NRBSITES
                     J8     = NRBSITES*(J2-1) + J
                     RSS(:) = R(J7,:) - R(J8,:)
                     DSS    = DSQRT(DOT_PRODUCT(RSS(:),RSS(:)))

                     IF (DSS.LT.DDMCUT) THEN
                        FCTR   = EXP(RHO*(MREQ-DSS))
                        ENERGY = ENERGY + BEPS*((1.D0 - FCTR)*(1.D0 - FCTR) - 1.D0)

                        IF (GTEST) THEN
                           DVDR       = 2.D0*BEPS*RHO*FCTR*(1.D0-FCTR)
                           G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:)
                           G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:)

                           G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(RSS(:),DR1(J7,:))
                           G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(RSS(:),DR2(J7,:))
                           G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(RSS(:),DR3(J7,:))

                           G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(RSS(:),DR1(J8,:))
                           G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(RSS(:),DR2(J8,:))
                           G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(RSS(:),DR3(J8,:))
                        ENDIF

                     ENDIF

                  ENDDO

               ENDDO

            ENDIF

         ENDDO

      ENDDO


      END SUBROUTINE DODECAMORSE
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDDM()
!
!     sites definition for colloidal rocks
!     the effective hard spheres, diameter HSEFF, are circumscribed by a unit
!     sphere
!
      USE COMMONS, ONLY: SITE, HSEFF
      IMPLICIT NONE
      DOUBLE PRECISION :: PHI, PHI1, FCTR
      INTEGER          :: I

      PHI=(1.D0+DSQRT(5.D0))/2.D0
      PHI1=1.D0/PHI
      FCTR=1.D0/(2.D0*DSQRT(3.D0))

!     vertices of a dodecahedron on a unit sphere

      SITE(1,:) =(/ 0.D0, PHI1, PHI /)*FCTR
      SITE(2,:) =(/ 0.D0,-PHI1, PHI /)*FCTR
      SITE(3,:) =(/ 1.D0, 1.D0, 1.D0/)*FCTR
      SITE(4,:) =(/-1.D0, 1.D0, 1.D0/)*FCTR
      SITE(5,:) =(/-1.D0,-1.D0, 1.D0/)*FCTR
      SITE(6,:) =(/ 1.D0,-1.D0, 1.D0/)*FCTR
      SITE(7,:) =(/ PHI , 0.D0, PHI1/)*FCTR
      SITE(8,:) =(/-PHI , 0.D0, PHI1/)*FCTR
      SITE(9,:) =(/ PHI1, PHI , 0.D0/)*FCTR
      SITE(10,:)=(/-PHI1, PHI , 0.D0/)*FCTR
      SITE(11,:)=(/-PHI1,-PHI , 0.D0/)*FCTR
      SITE(12,:)=(/ PHI1,-PHI , 0.D0/)*FCTR
      SITE(13,:)=(/ PHI , 0.D0,-PHI1/)*FCTR
      SITE(14,:)=(/-PHI , 0.D0,-PHI1/)*FCTR
      SITE(15,:)=(/ 1.D0, 1.D0,-1.D0/)*FCTR
      SITE(16,:)=(/-1.D0, 1.D0,-1.D0/)*FCTR
      SITE(17,:)=(/-1.D0,-1.D0,-1.D0/)*FCTR
      SITE(18,:)=(/ 1.D0,-1.D0,-1.D0/)*FCTR
      SITE(19,:)=(/ 0.D0, PHI1,-PHI /)*FCTR
      SITE(20,:)=(/ 0.D0,-PHI1,-PHI /)*FCTR

!     the sites are inside a unit sphere and meet the sphere at the effective
!     hard sphere radius

      SITE=SITE*(1.D0-HSEFF)

      END SUBROUTINE DEFDDM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWDDM()
!
!     converts to sites representation and writes to dodecamorse.xyz
!
      USE COMMONS, ONLY: NRBSITES, SITE, NSAVE, HSEFF
      USE QMODULE
      IMPLICIT NONE
      INTEGER          :: J1, J2, J3, J4, J5, J7 
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(NRBSITES*3), P3(3,3), RCC(3), DCC
      LOGICAL          :: GTEST, GOODSTRUCT

      OPEN(UNIT=26, FILE='dodecamorse.xyz', STATUS='UNKNOWN')

      GTEST = .FALSE.

      DO J1 = 1, NSAVE
         WRITE(26,'(I6,F8.3)') (QMINNATOMS(J1)/2)*(NRBSITES+1),HSEFF
         WRITE(26,*) J1, QMIN(J1), FF(J1)
  
         DO J3 = 1, QMINNATOMS(J1)/2
            J5   = 3*J3
            J7   = 3*QMINNATOMS(J1)/2 + J5
            P(:) = QMINP(J1,J7-2:J7)
            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)
            WRITE(26,'(A4,3F20.10)')'O',10.D0*QMINP(J1,J5-2:J5)

            DO J2=1,NRBSITES
               RBCOORDS(3*J2-2:3*J2) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               WRITE(26,'(A4,3F20.10)')'LA',10.D0*RBCOORDS(3*J2-2:3*J2)
            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWDDM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DDMCONDENSE(COORDS)
!
!     finds the least distance between rock i and rocks j and normalises Rij
!     condenses loose structures and helps to repair coalesced ones
!     called after step-taking and before quenching
!
      USE COMMONS, ONLY: NATOMS,MREQ,HSEFF
      IMPLICIT NONE
      INTEGER           :: J1,J2,J3,J4
      DOUBLE PRECISION  :: D,MIND
      DOUBLE PRECISION  :: COORDS(3*NATOMS),R(3),MINR(3)

!     IF (.TRUE.) RETURN

      DO J1=1,NATOMS/2
         J3=3*J1
         MIND=1.D10
         MINR(:)=0.D0
         DO J2=J1+1,NATOMS/2
            J4=3*J2
            R(:)=COORDS(J4-2:J4)-COORDS(J3-2:J3)
            D=SQRT(DOT_PRODUCT(R,R))
            IF (D<MIND) THEN
               MIND=D
               MINR(:)=R(:)
            ENDIF
         ENDDO
         COORDS(J3-2:J3)=COORDS(J3-2:J3)+(MIND-1.D0+(HSEFF-MREQ))*MINR(:)/MIND
      ENDDO

      END SUBROUTINE DDMCONDENSE

!     ----------------------------------------------------------------------------------------------
