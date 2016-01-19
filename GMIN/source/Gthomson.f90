!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

      SUBROUTINE GTHOMSON(X,V,ETHOMSON,GTEST)

      USE COMMONS, ONLY: GTHOMMET, GTHOMPOT, GThomsonRho, GThomsonSigma, NATOMS
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT)  :: X(*)
      DOUBLE PRECISION, INTENT(OUT) :: V(*), ETHOMSON
      LOGICAL, INTENT(IN) :: GTEST
      INTEGER J1, J2
      DOUBLE PRECISION DIST
      DOUBLE PRECISION :: TMPCOORDS(3*NATOMS), DR(3), Gradient11(3), Gradient12(3), Gradient21(3), Gradient22(3)

! jwrm2> Reset the range of the polar coordinates
      CALL GTHOMSONWRAPPOLAR (X(1:3*NATOMS), NATOMS)

      CALL GTHOMSONANGTOC(TMPCOORDS(1:3*NATOMS), X(1:3*NATOMS), NATOMS)

      ETHOMSON=0.0D0
      V(1:3*NATOMS)=0.0D0

      DO J1=1,NATOMS-1

         IF (GTEST) THEN
            IF ( GTHOMMET .EQ. 1) THEN
               CALL GRADMETRICCYLINDER (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
            ELSEIF ( GTHOMMET .EQ. 2) THEN
               CALL GRADMETRICCATENOID (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
            ELSEIF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN
               CALL GRADMETRICUNDULOID (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
            ELSEIF ( GTHOMMET .EQ. 5 ) THEN
               CALL GRADMETRICSPHERE (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
            ELSEIF ( GTHOMMET .EQ. 6 ) THEN
               CALL GRADMETRICMOBIUS (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
            ENDIF
         ENDIF

         DO J2=J1+1,NATOMS
            DR(1) = TMPCOORDS(3*J1-2)-TMPCOORDS(3*J2-2)
            DR(2) = TMPCOORDS(3*J1-1)-TMPCOORDS(3*J2-1)
            DR(3) = TMPCOORDS(3*J1  )-TMPCOORDS(3*J2  )
            DIST  = SQRT(DOT_PRODUCT(DR, DR))

            IF (GTHOMPOT .EQ. 1) THEN
               ETHOMSON = ETHOMSON + 1/DIST ! Coulomb
            ELSEIF (GTHOMPOT .EQ. 2) THEN 
               ETHOMSON = ETHOMSON + 1/DIST**3 ! 1/R^3
            ELSEIF (GTHOMPOT .EQ. 3) THEN
               ETHOMSON = ETHOMSON + 1/DIST*EXP(-1.0D0*DIST/GThomsonSigma)                !Yukawa
            ELSEIF (GTHOMPOT .EQ. 4) THEN
               ETHOMSON = ETHOMSON + (GThomsonSigma/DIST)**12 - (GThomsonSigma/DIST)**6  ! LJ
            ELSEIF (GTHOMPOT .EQ. 5) THEN
               ETHOMSON = ETHOMSON + (GThomsonSigma/DIST)**12                            !Repulsive LJ 
            ELSEIF (GTHOMPOT .EQ. 6) THEN
               ETHOMSON = ETHOMSON + (1.0D0 - EXP(GThomsonRho*(GThomsonSigma-DIST)))**2 -1.0D0 ! Morse
            ENDIF


            IF (GTEST) THEN

               IF (GTHOMPOT .EQ. 1) THEN
                  V(3*J1-2) = V(3*J1-2) - 1.0D0/DIST**3 * DOT_PRODUCT(DR,Gradient11)
                  V(3*J1-1) = V(3*J1-1) - 1.0D0/DIST**3 * DOT_PRODUCT(DR,Gradient12)               
               ELSEIF (GTHOMPOT .EQ. 2) THEN 
                  V(3*J1-2) = V(3*J1-2) - 3.0D0/DIST**5 * DOT_PRODUCT(DR,Gradient11)
                  V(3*J1-1) = V(3*J1-1) - 3.0D0/DIST**5 * DOT_PRODUCT(DR,Gradient12)               
               ELSEIF (GTHOMPOT .EQ. 3) THEN
                  V(3*J1-2) = V(3*J1-2) - (1.0D0/DIST**3 + 1.0D0/GThomsonSigma/DIST**2) * EXP(-1.0D0*DIST/GThomsonSigma) &
                       * DOT_PRODUCT(DR,Gradient11)
                  V(3*J1-1) = V(3*J1-1) - (1.0D0/DIST**3 + 1.0D0/GThomsonSigma/DIST**2) * EXP(-1.0D0*DIST/GThomsonSigma) &
                       * DOT_PRODUCT(DR,Gradient12)               
               ELSEIF (GTHOMPOT .EQ. 4) THEN
                  V(3*J1-2) = V(3*J1-2) - ((12.0D0*GThomsonSigma**12/DIST**14) - (6.0D0*GThomsonSigma**6/DIST**8)) &
                       * DOT_PRODUCT(DR,Gradient11)
                  V(3*J1-1) = V(3*J1-1) - ((12.0D0*GThomsonSigma**12/DIST**14) - (6.0D0*GThomsonSigma**6/DIST**8)) & 
                       * DOT_PRODUCT(DR,Gradient12)
               ELSEIF (GTHOMPOT .EQ. 5) THEN
                  V(3*J1-2) = V(3*J1-2) - (12.0D0*GThomsonSigma**12/DIST**14)* DOT_PRODUCT(DR,Gradient11) ! Repulsive LJ
                  V(3*J1-1) = V(3*J1-1) - (12.0D0*GThomsonSigma**12/DIST**14)* DOT_PRODUCT(DR,Gradient12) ! Replusive LJ
               ELSEIF (GTHOMPOT .EQ. 6) THEN
                  V(3*J1-2) = V(3*J1-2) + 2.0D0 * (1.0D0 - EXP(GThomsonRho*(GThomsonSigma-DIST))) * & 
                       EXP(GThomsonRho*(GThomsonSigma-DIST)) * GThomsonRho / DIST * DOT_PRODUCT(DR,Gradient11) ! Morse
                  V(3*J1-1) = V(3*J1-1) + 2.0D0 * (1.0D0 - EXP(GThomsonRho*(GThomsonSigma-DIST))) * &
                       EXP(GThomsonRho*(GThomsonSigma-DIST)) * GThomsonRho / DIST * DOT_PRODUCT(DR,Gradient12) ! Morse
               ENDIF
 
               IF ( GTHOMMET .EQ. 1) THEN
                  CALL GRADMETRICCYLINDER (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
               ELSEIF ( GTHOMMET .EQ. 2) THEN
                  CALL GRADMETRICCATENOID (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
               ELSEIF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN
                  CALL GRADMETRICUNDULOID (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
               ELSEIF ( GTHOMMET .EQ. 5 ) THEN
                  CALL GRADMETRICSPHERE (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
               ELSEIF ( GTHOMMET .EQ. 6 ) THEN
                  CALL GRADMETRICMOBIUS (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
               ENDIF

               IF (GTHOMPOT .EQ. 1) THEN
                  V(3*J2-2) = V(3*J2-2) + 1.0D0/DIST**3 * DOT_PRODUCT(DR,Gradient21)               
                  V(3*J2-1) = V(3*J2-1) + 1.0D0/DIST**3 * DOT_PRODUCT(DR,Gradient22)
               ELSEIF (GTHOMPOT .EQ. 2) THEN 
                  V(3*J2-2) = V(3*J2-2) + 3.0D0/DIST**5 * DOT_PRODUCT(DR,Gradient21)               
                  V(3*J2-1) = V(3*J2-1) + 3.0D0/DIST**5 * DOT_PRODUCT(DR,Gradient22)
               ELSEIF (GTHOMPOT .EQ. 3) THEN
                  V(3*J2-2) = V(3*J2-2) + (1.0D0/DIST**3 + 1.0D0/GThomsonSigma/DIST**2) * EXP(-1.0D0*DIST/GThomsonSigma) &
                       * DOT_PRODUCT(DR,Gradient21)               
                  V(3*J2-1) = V(3*J2-1) + (1.0D0/DIST**3 + 1.0D0/GThomsonSigma/DIST**2) * EXP(-1.0D0*DIST/GThomsonSigma) &
                       * DOT_PRODUCT(DR,Gradient22)
               ELSEIF (GTHOMPOT .EQ. 4) THEN
                  V(3*J2-2) = V(3*J2-2) + ((12.0D0*GThomsonSigma**12/DIST**14) - (6.0D0*GThomsonSigma**6/DIST**8))&
                       * DOT_PRODUCT(DR,Gradient21)
                  V(3*J2-1) = V(3*J2-1) + ((12.0D0*GThomsonSigma**12/DIST**14) - (6.0D0*GThomsonSigma**6/DIST**8))&
                       * DOT_PRODUCT(DR,Gradient22)
               ELSEIF (GTHOMPOT .EQ. 5) THEN
                  V(3*J2-2) = V(3*J2-2) + (12.0D0*GThomsonSigma**12/DIST**14)* DOT_PRODUCT(DR,Gradient21) ! Repulsive LJ
                  V(3*J2-1) = V(3*J2-1) + (12.0D0*GThomsonSigma**12/DIST**14)* DOT_PRODUCT(DR,Gradient22) ! Repulsive LJ
               ELSEIF (GTHOMPOT .EQ. 6) THEN
                  V(3*J2-2) = V(3*J2-2) - 2.0D0 * (1.0D0 - EXP(GThomsonRho*(GThomsonSigma-DIST))) &
                       * EXP(GThomsonRho*(GThomsonSigma-DIST)) * GThomsonRho / DIST * DOT_PRODUCT(DR,Gradient21)
                  V(3*J2-1) = V(3*J2-1) - 2.0D0 * (1.0D0 - EXP(GThomsonRho*(GThomsonSigma-DIST))) &
                       * EXP(GThomsonRho*(GThomsonSigma-DIST)) * GThomsonRho / DIST * DOT_PRODUCT(DR,Gradient22) 
               ENDIF
               !PRINT *, 1.0D0/DIST**3 * DOT_PRODUCT(DR,Gradient2), (1/DIST2 - 1/DIST)/0.00001

            ENDIF

         ENDDO
      ENDDO

      RETURN
    END SUBROUTINE GTHOMSON

!----------------------------------------------------------------------------------------------------------------------------------
!  Subroutine to convert Cartesians to theta, phi.
!  Order of output coordinates is phi1, theta1, r1, phi2, theta2, r2,...
!  where phi is the azimuthal angle and theta the vertical angle.
!  r is set to zero for all atoms and recalculated from the surface parameterisation when converting back to Cartesians
!
!  COORDS: the input coordinates, which should be Cartesians
!  P: the output coordinates in polars
!  NATOMS: the number of bodies

      SUBROUTINE GTHOMSONCTOANG(COORDS, P, NATOMS)
      USE COMMONS, ONLY: DEBUG, GTHOMSONZ, GTHOMMET, GTREFU, GTMU, GTM, GTN, MYUNIT
      IMPLICIT NONE
      INTEGER J1
      DOUBLE PRECISION, INTENT(INOUT)  :: COORDS(*)
      DOUBLE PRECISION, INTENT(OUT) :: P(*)
      INTEGER, INTENT(IN) :: NATOMS
      DOUBLE PRECISION :: PHI, U
      DOUBLE PRECISION RADIUS
      DOUBLE PRECISION, PARAMETER :: TOLERANCE = 1.D-5, PI = 4.0D0*ATAN(1.0D0), HALFPI = 2.0D0*ATAN(1.0D0)
      LOGICAL POLAR

! jwrm2> Attempt to detect if the coordintes are already polar. We'll assume they are if all the x coordinates
!        are between 0 and 2pi, all the y coordinates are between 0 and pi and all the z coordinates are zero, 
!        This is not guaranteed to mean they're polar, but should be ok since it shouldn't be a 2D system.
      POLAR = .TRUE.
      DO J1 = 3, 3*NATOMS, 3
        IF (COORDS(J1-2) .GT. 2*PI   .OR. COORDS(J1-2) .LT. 0.D0) POLAR = .FALSE. ! x coord
        IF (COORDS(J1-1) .GT. PI     .OR. COORDS(J1-1) .LT. 0.D0) POLAR = .FALSE. ! y coord
        IF (COORDS(J1) .GT. TOLERANCE .OR. COORDS(J1) .LT. -TOLERANCE) POLAR = .FALSE. ! z coord
      END DO
      IF (POLAR) THEN
        IF (DEBUG) WRITE (MYUNIT,*) 'GTHOMSONCTOANG> warning - incoming coordinates are already polar, skipping conversion'
        P(1:3*NATOMS) = COORDS(1:3*NATOMS)
        RETURN
      END IF

      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF

      DO J1 = 1, NATOMS

         IF ( (COORDS(3*J1-2) .GE. 0.0D0) .AND. (COORDS(3*J1-1) .GE. 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-2)) < 1.0D-5) THEN
               P(3*J1-2) = HALFPI
            ELSE IF ( ABS(COORDS(3*J1-1)) < 1.0D-5) THEN
               P(3*J1-2) = 0.0D0
            ELSE
               P(3*J1-2) = ATAN(COORDS(3*J1-1)/COORDS(3*J1-2))
            ENDIF
         ELSEIF ( (COORDS(3*J1-2) < 0.0D0) .AND. (COORDS(3*J1-1) .GE. 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-1)) < 1.0D-5) THEN
               P(3*J1-2) = 2*HALFPI
            ELSE
               P(3*J1-2) = 2*HALFPI - ATAN(COORDS(3*J1-1)/(-COORDS(3*J1-2)))
            ENDIF
         ELSEIF ( (COORDS(3*J1-2) < 0.0D0) .AND. (COORDS(3*J1-1) < 0.0D0) ) THEN
               P(3*J1-2) = 2*HALFPI + ATAN(COORDS(3*J1-1)/COORDS(3*J1-2))
         ELSEIF ( (COORDS(3*J1-2) .GE. 0.0D0) .AND. (COORDS(3*J1-1) < 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-2)) < 1.0D-5) THEN
               P(3*J1-2) = 3*HALFPI
            ELSE
               P(3*J1-2) = 4*HALFPI - ATAN(-COORDS(3*J1-1)/COORDS(3*J1-2))
            ENDIF
         ENDIF

         IF ( (GTHOMMET .EQ. 1) .OR. (GTHOMMET .EQ. 2) ) THEN

            IF (COORDS(3*J1) > GThomsonZ) COORDS(3*J1) = GThomsonZ
            IF (COORDS(3*J1) <-GThomsonZ) COORDS(3*J1) = -GThomsonZ
            P(3*J1-1) = ACOS(COORDS(3*J1)/GThomsonZ)

         ELSE IF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) )  THEN

            IF ( ABS(COS(P(3*J1-2))) < 1.0D-5 ) THEN
               u = ASIN( (COORDS(3*J1-1)**2- GTn)/GTm )/GTmu
            ELSE
               u = ASIN( (COORDS(3*J1-2)/COS(P(3*J1-2)))**2/GTm - GTn/GTm )/GTmu
            ENDIF

            IF ( GTHOMMET .EQ. 3 ) THEN
               IF (COORDS(3*J1) .GE. 0.0D0) u = u+2*pi/GTmu
               IF (COORDS(3*J1) <    0.0D0) u = pi/GTmu-u
            ELSEIF ( GTHOMMET .EQ. 4 ) THEN
               IF (COORDS(3*J1) .GE. 0.0D0) u = pi/GTmu-u
            ENDIF

            phi = (u - GTrefU/GTmu)*GTmu/pi/GThomsonZ
            if (phi > 1.0D0) phi = 1.0D0
            if (phi <-1.0D0) phi =-1.0D0
            P(3*J1-1) = ACOS(phi)

         ELSE IF ( GTHOMMET .EQ. 5 )  THEN
            
            IF ( COORDS(3*J1) < 0.0D0 ) THEN               
               P(3*J1-1) = 2*HALFPI - ACOS(-COORDS(3*J1)/RADIUS)
            ELSE
               IF ( COORDS(3*J1)/RADIUS .GE. 1.0D0 ) THEN
                  P(3*J1-1) = 0.0D0
               ELSE IF ( COORDS(3*J1)/RADIUS .LE. -1.0D0 ) THEN
                  P(3*J1-1) = 2*HALFPI
               ELSE
                  P(3*J1-1) = ACOS(COORDS(3*J1)/RADIUS)
               ENDIF
            ENDIF

         ENDIF
         
         IF ( P(3*J1-2) > 4*HALFPI ) P(3*J1-2) = P(3*J1-2) - 4*HALFPI 
         IF ( P(3*J1-2) < 0.0D0    ) P(3*J1-2) = P(3*J1-2) + 4*HALFPI 
         IF ( P(3*J1-1) > 4*HALFPI ) P(3*J1-1) = P(3*J1-1) - 4*HALFPI 
         IF ( P(3*J1-1) < 0.0D0    ) P(3*J1-1) = P(3*J1-1) + 4*HALFPI 

      ENDDO

      DO J1 = 3, 3*NATOMS, 3
         P(J1) = 0.0D0
      END DO

      RETURN
      
    END SUBROUTINE GTHOMSONCTOANG

!----------------------------------------------------------------------------------------------------------------------------------
!  Subroutine to convert theta, phi to Cartesians.
!  Coordinates come in as phi1, theta1, r1, phi2, theta2, r2,...
!  where phi is the azimuthal angle and theta the vertical angle.
!  All rs are zero and the radius is calculated here from the surface parameterisation
!
!  COORDS: output coordinates, in Cartesians
!  P: input coordinates in polars
!  NATOMS: the number of bodies

      SUBROUTINE GTHOMSONANGTOC(COORDS,P,NATOMS)
      USE COMMONS, ONLY: DEBUG, GTHOMSONC, GTHOMSONZ, GTHOMMET, GTREFU, GTREFZ, GTMU, GTK, GTM, GTN, GTA, GTC
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS
      INTEGER J1
      DOUBLE PRECISION, INTENT(OUT)  :: COORDS(*)
      DOUBLE PRECISION, INTENT(IN)   :: P(*)
      DOUBLE PRECISION :: FELINT, SELINT, u
      DOUBLE PRECISION :: ARG1, COSHARG
      DOUBLE PRECISION :: RADIUS
      DOUBLE PRECISION, PARAMETER :: TOLERANCE = 1.D-5, PI = 4.0D0*ATAN(1.0D0)
      LOGICAL POLAR

! jwrm2> Detect if the coordintes are already Cartesian. If the r coordinates are not zero, then
!        the coordinates are already Cartesian. They could still be 
!        Cartesian if they are all zero, but this it unlikely since it shouldn't be a flat 2D system.
      POLAR = .TRUE.
      DO J1 = 3, 3*NATOMS, 3
        IF (P(J1) .GT. TOLERANCE .OR. P(J1) .LT. -TOLERANCE) POLAR = .FALSE. ! r coord
      END DO
      IF (.NOT. POLAR) THEN
        IF (DEBUG) PRINT *, 'GTHOMSONANGTOC> warning - incoming coordinates are already Cartesian, skipping conversion'
        COORDS(1:3*NATOMS) = P(1:3*NATOMS)
        RETURN
      END IF

      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF

      DO J1 = 1, NATOMS
         IF ( GTHOMMET .EQ. 1) THEN
            COORDS(3*(J1-1)+1)= GThomsonC * COS(P(3*J1-2))
            COORDS(3*(J1-1)+2)= GThomsonC * SIN(P(3*J1-2))
            COORDS(3*(J1-1)+3)= GThomsonZ * COS(P(3*J1-1))

         ELSE IF ( GTHOMMET .EQ. 2) THEN
            COSHARG = GThomsonC * COSH(GThomsonZ/GThomsonC*COS(P(3*J1-1)))
            COORDS(3*(J1-1)+1)= COSHARG * COS(P(3*J1-2))
            COORDS(3*(J1-1)+2)= COSHARG * SIN(P(3*J1-2))
            COORDS(3*(J1-1)+3)= GThomsonZ * COS(P(3*J1-1))

         ELSEIF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN

            u = GTrefU/GTmu + GThomsonZ*pi/GTmu * COS(P(3*J1-1))
            CALL EllipIntegral(Felint, Selint, (GTmu*u/2.0D0-pi/4.0D0), GTk)
            ARG1 = SQRT(GTm*SIN(GTmu*u)+GTn)
            COORDS(3*(J1-1)+1) = ARG1 * COS(P(3*J1-2))
            COORDS(3*(J1-1)+2) = ARG1 * SIN(P(3*J1-2))
            COORDS(3*(J1-1)+3) = GTa*Felint + GTc*Selint

            IF (GTHOMMET .EQ. 3) THEN
               COORDS(3*(J1-1)+3) = COORDS(3*(J1-1)+3) - GTrefZ
            ENDIF

         ELSE IF ( GTHOMMET .EQ. 5) THEN

            ARG1 = RADIUS * SIN(P(3*J1-1))
            COORDS(3*J1-2)= ARG1 * COS(P(3*J1-2))
            COORDS(3*J1-1)= ARG1 * SIN(P(3*J1-2))
            COORDS(3*J1  )= RADIUS * COS(P(3*J1-1))

         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE GTHOMSONANGTOC

!----------------------------------------------------------------
! Make sure polar coordinates lie within the correct ranges
! Coordinates come in as phi1, theta1, r1, phi2, theta2, r2,...
!
! P: The polar coordinates
      SUBROUTINE GTHOMSONWRAPPOLAR (P, NATOMS)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: NATOMS
        DOUBLE PRECISION, INTENT(INOUT) :: P(*)
        DOUBLE PRECISION, PARAMETER :: TOLERANCE = 1.D-5, PI = 4.0D0*ATAN(1.0D0)
        INTEGER :: I
        LOGICAL :: POLAR

! jwrm2> Detect if the coordintes are not polar. If the r coordinates are not zero, then
!        the coordinates are already Cartesian. They could still be 
!        Cartesian if they are all zero, but this it unlikely since it shouldn't be a flat 2D system.
        POLAR = .TRUE.
        DO I = 3, 3*NATOMS, 3
          IF (P(I) .GT. TOLERANCE .OR. P(I) .LT. -TOLERANCE) POLAR = .FALSE. ! r coord
        END DO
        IF (.NOT. POLAR) THEN
          PRINT *, 'GTHOMSONWRAPPOLAR> warning - incoming coordinates are Cartesian, skipping'
          RETURN
        END IF

        DO I = 1, 3*NATOMS, 3
! Theta coordinate
          IF (P(I+1) .GT. PI) THEN
            P(I+1) = 2*PI - P(I+1)
            P(I)   = PI + P(I)
          END IF
          IF (P(I+1) .LT. 0.D0) THEN
            P(I+1) = - P(I+1)
            P(I)   = PI + P(I)
          END IF
! Phi coordinate
          DO WHILE (P(I) .GT. 2*PI)
            P(I) = P(I) - 2*PI
          END DO 
          DO WHILE (P(I) .LT. 0.D0)
            P(I) = P(I) + 2*PI
          END DO
! Reset r coordiantes to exactly zero
          P(I+2) = 0.D0
        END DO 

      END SUBROUTINE GTHOMSONWRAPPOLAR

!----------------------------------------------------------------
! take step routine      

      SUBROUTINE TAKESTEPGTHOMSON ()

      USE COMMONS, ONLY : COORDS, NATOMS, TMOVE, OMOVE, STEP, OSTEP, GTHOMMET, GTHOMPOT, GThomsonSigma
      USE ROTATIONS
      IMPLICIT NONE
      DOUBLE PRECISION P(3*NATOMS), LOCALSTEP, RANDOM, DPRAND, CARTCOORDS(3*NATOMS)
      DOUBLE PRECISION RMI(3,3)
      INTEGER NP, J1
      LOGICAL PERCT

      NP = 1
      LOCALSTEP = STEP(NP)

      IF ( GTHOMMET .EQ. 5) THEN

         CALL GTHOMSONANGTOC(CARTCOORDS(1:3*NATOMS), COORDS(1:3*NATOMS, NP), NATOMS)

         DO J1 = 1, NATOMS
            RMI = rot_aa2mx(rot_small_random_aa(LOCALSTEP))
            CARTCOORDS(3*J1-2:3*J1) = MATMUL(RMI(:,:), CARTCOORDS(3*J1-2:3*J1))
         ENDDO

         CALL GTHOMSONCTOANG(CARTCOORDS(1:3*NATOMS), COORDS(1:3*NATOMS, NP), NATOMS)

      ELSE

         CALL GTHOMSONCTOANG(COORDS(:, NP), P, NATOMS)
         DO J1 = 1, NATOMS
            LOCALSTEP = 0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)
            RANDOM            = (DPRAND() - 0.5D0)*2.0D0
            P(3*J1-2) = P(3*J1-2) + LOCALSTEP*RANDOM
            
            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)
            RANDOM            = (DPRAND() - 0.5D0)*2.0D0
            P(3*J1-1)   = P(3*J1-1) + LOCALSTEP*RANDOM   

            P(3*J1) = 0.D0         
         ENDDO
         CALL GTHOMSONANGTOC(COORDS(:, NP), P, NATOMS)
      
      ENDIF

!      DO J1 = 1, NATOMS
!         LOCALSTEP = 0.0D0
!         IF (TMOVE(NP)) LOCALSTEP = STEP(NP)
!         RANDOM            = (DPRAND() - 0.5D0)*2.0D0
!         P(3*J1-2) = P(3*J1-2) + LOCALSTEP*RANDOM         
!      ENDDO      
!      CALL GTHOMSONANGTOC(COORDS, P,NATOMS)
!      COORDS(1:3*NATOMS,NP) = P(1:3*NATOMS)
!      DO J1 = 1, NATOMS
!         LOCALSTEP = 0.0D0
!         IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)
!         RANDOM            = (DPRAND() - 0.5D0)*2.0D0
!         COORDS(3*J1,NP)   = COORDS(3*J1,NP) + LOCALSTEP*RANDOM          
!      ENDDO


! hk286 - 12/5/2013
      IF ( (GTHOMMET .EQ. 5) .AND. ((GTHOMPOT .EQ. 4) .OR. (GTHOMPOT .EQ. 6)) ) THEN
! jwrm2> PERCSPHERE expects Cartesian coordinates, but COORDS is storing polars
         CALL GTHOMSONANGTOC(CARTCOORDS(1:3*NATOMS), COORDS(1:3*NATOMS, NP), NATOMS)
         CALL PERCSPHERE(CARTCOORDS(1:3*NATOMS), NATOMS, 1.5D0*GThomsonSigma, PERCT)
         CALL GTHOMSONCTOANG(CARTCOORDS(1:3*NATOMS), COORDS(1:3*NATOMS, NP), NATOMS)
      ENDIF

    END SUBROUTINE TAKESTEPGTHOMSON

!----------------------------------------------------------------
! newrestart routine      

      SUBROUTINE NEWRESTARTGTHOMSON ()

      USE COMMONS, ONLY : COORDS, NATOMS
      IMPLICIT NONE
      DOUBLE PRECISION P(3*NATOMS), LOCALSTEP, RANDOM, DPRAND
      DOUBLE PRECISION EREAL, GRAD(3*NATOMS)
      INTEGER NP, J1
      
      PRINT *, "NEWRESTART IS CALLED"
      NP = 1
      CALL GTHOMSONCTOANG(COORDS,P,NATOMS)

      DO J1 = 1, NATOMS

         LOCALSTEP = 0.5D0
         RANDOM            = (DPRAND() - 0.5D0)*2.0D0
         P(3*J1-2) = P(3*J1-2) + LOCALSTEP*RANDOM
         
         LOCALSTEP = 0.5D0
         RANDOM            = (DPRAND() - 0.5D0)*2.0D0
         P(3*J1-1)   = P(3*J1-1) + LOCALSTEP*RANDOM

         P(3*J1) = 0.D0
          
      ENDDO
      
      CALL GTHOMSONANGTOC(COORDS(1:3*NATOMS, NP), P, NATOMS)
      CALL GTHOMSON(COORDS(1:3*NATOMS, NP), GRAD, EREAL, .FALSE.)
      
      PRINT *, "NEWRESTART DONE"

    END SUBROUTINE NEWRESTARTGTHOMSON

!----------------------------------------------------------------
! Initialization of random coordinates
      SUBROUTINE INIGTHOMSON ()

      USE COMMONS, ONLY : NATOMS, GTHOMMET, GTHOMSONZ, GTHOMSONC2, GTHOMSONC, GTrefU, GTrefZ, GTmu, GTk, GTm, GTn, GTa, GTc
      USE ROTATIONS
      USE VEC3
      IMPLICIT NONE
      DOUBLE PRECISION P(3*NATOMS), RANDOM, DPRAND, pi, c, a, Felint, Selint
      INTEGER NP, J1
      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),itime1
      DOUBLE PRECISION :: RADIUS, COORDS(3*NATOMS)

      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF
      NP = 1

      IF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN

         pi = 4.0D0*ATAN(1.0D0)
         IF (GTHOMMET .EQ. 3) THEN
            GTrefU = 1.50D0*pi
         ELSEIF (GTHOMMET .EQ. 4) THEN
            GTrefU = 0.50D0*pi
         ENDIF
         IF (GThomsonC > GThomsonC2) THEN
            c = GThomsonC
            a = GThomsonC2
         ELSE
            c = GThomsonC2
            a = GThomsonC
         ENDIF
         GTa = a
         GTc = c
         GTmu = 2.0D0/(a+c)
         GTk = SQRT(1-(a/c)**2)
         GTm = (c**2-a**2)/2.0D0
         GTn = (c**2+a**2)/2.0D0      
      
      ENDIF
      IF (GTHOMMET .EQ. 3) THEN
         CALL EllipIntegral(Felint, Selint, (GTrefU/2.0D0-pi/4.0D0), GTk)
         GTrefZ = GTa*Felint + GTc*Selint
      ENDIF

      CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
      itime1= values(7)*39 + values(8)
      CALL SDPRND(itime1+NP)

!      OPEN(UNIT = 28, FILE = 'coords')            
!      DO J1 = 1, NATOMS
!         READ(28, *) COORDS(3*J1-2:3*J1)
!      ENDDO
!      CLOSE (UNIT=28)
!      PRINT *, COORDS(1:3)
!      CALL GTHOMSONCTOANG(COORDS,P,NATOMS)
!      PRINT *, P(1:2)
!      CALL GTHOMSONANGTOC(P,NATOMS)
!      PRINT *, P(1:3)
     
      IF ( GTHOMMET .EQ. 5 ) THEN

         DO J1 = 1, NATOMS
            COORDS(3*J1-2:3*J1) = RADIUS * vec_random()
         ENDDO

         OPEN(UNIT = 28, FILE = 'coordsini', STATUS = 'REPLACE')            
         DO J1 = 1, NATOMS
            WRITE(28, *) COORDS(3*J1-2:3*J1)
         ENDDO
         CLOSE (UNIT=28)

      ELSE

         DO J1 = 1, NATOMS
            RANDOM  = DPRAND()*8*ATAN(1.0D0)
            P(3*J1-2) = RANDOM
            RANDOM  = DPRAND()*4*ATAN(1.0D0)
            P(3*J1-1) = RANDOM 
            P(3*J1) = 0.D0         
         ENDDO

         CALL GTHOMSONANGTOC(COORDS, P, NATOMS)
         OPEN(UNIT = 28, FILE = 'coordsini', STATUS = 'REPLACE')            
         DO J1 = 1, NATOMS
            WRITE(28, *) COORDS(3*J1-2:3*J1)
         ENDDO
         CLOSE (UNIT=28)

      ENDIF
         
    END SUBROUTINE INIGTHOMSON


!----------------------------------------------------------------
! CYLINDER     
      SUBROUTINE GRADMETRICCYLINDER (X, Gradient1, Gradient2)

      USE COMMONS, ONLY : GThomsonC, GThomsonZ
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: Gradient1(3), Gradient2(3)

      Gradient1(1) = -GThomsonC * SIN(X(1))
      Gradient1(2) =  GThomsonC * COS(X(1))
      Gradient1(3) =  0.0D0

      Gradient2(1) = 0.0D0
      Gradient2(2) = 0.0D0
      Gradient2(3) = -GThomsonZ * SIN(X(2))

      END SUBROUTINE GRADMETRICCYLINDER

!----------------------------------------------------------------
! catenoid      
      SUBROUTINE GRADMETRICCATENOID (X, Gradient1, Gradient2)

      USE COMMONS, ONLY : GThomsonC, GThomsonZ
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: Gradient1(3), Gradient2(3)
      DOUBLE PRECISION ARGHY, COSHA, SINHA, SINX2, SINX1, COSX1
      
      ARGHY = GThomsonZ/GThomsonC*COS(X(2))
      COSHA = COSH(ARGHY)
      SINHA = SINH(ARGHY)
      SINX2 = SIN(X(2))
      SINX1 = SIN(X(1))
      COSX1 = COS(X(1))
      Gradient1(1) = -GThomsonC * COSHA * SINX1
      Gradient1(2) =  GThomsonC * COSHA * COSX1
      Gradient1(3) =  0.0D0

      Gradient2(1) = -GThomsonZ * SINHA * COSX1 * SINX2
      Gradient2(2) = -GThomsonZ * SINHA * SINX1 * SINX2
      Gradient2(3) = -GThomsonZ * SINX2

    END SUBROUTINE GRADMETRICCATENOID

!----------------------------------------------------------------
! unduloid      
      SUBROUTINE GRADMETRICUNDULOID (X, Gradient1, Gradient2)

      USE COMMONS, ONLY : GThomsonZ, GTrefU, GTmu, GTk, GTm, GTn, GTa, GTc
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: Gradient1(3), Gradient2(3)
      DOUBLE PRECISION :: pi, u
      DOUBLE PRECISION :: ARG1, ARG2, SINMU, CSINX2

      pi = 4.0D0*ATAN(1.0D0)

      u = GTrefU/GTmu + GThomsonZ*pi/GTmu * COS(X(2))
      
      ARG1 = SQRT(GTm*SIN(GTmu*u)+GTn)
      ARG2 = 0.5D0 / ARG1 * GTm * GTmu * COS(GTmu*u)
      SINMU = SIN(GTmu*u/2.0D0-pi/4.0D0)
      Gradient2(1) = ARG2* COS(X(1))
      Gradient2(2) = ARG2* SIN(X(1))
      Gradient2(3) = 1.0D0*(GTmu/2.0D0)*(GTa/SQRT(1-GTk**2*(SINMU)**2) + GTc*SQRT(1-GTk**2*(SINMU)**2))

      CSINX2 = GThomsonZ*pi/GTmu * SIN(X(2))
      Gradient2(1) = -Gradient2(1) *  CSINX2
      Gradient2(2) = -Gradient2(2) *  CSINX2
      Gradient2(3) = -Gradient2(3) *  CSINX2

      Gradient1(1) =-ARG1 * SIN(X(1))
      Gradient1(2) = ARG1 * COS(X(1))
      Gradient1(3) = 0.0D0


    END SUBROUTINE GRADMETRICUNDULOID

!----------------------------------------------------------------
! SPHERE
      SUBROUTINE GRADMETRICSPHERE (X, Gradient1, Gradient2)
        
      USE COMMONS, ONLY: GTHOMSONZ, GTHOMMET
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: Gradient1(3), Gradient2(3)
      DOUBLE PRECISION ARG1, ARG2, ARG3, ARG4
      DOUBLE PRECISION RADIUS
      
      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF

      ARG1 = SIN(X(1))
      ARG2 = COS(X(1))
      ARG3 = SIN(X(2))
      ARG4 = COS(X(2))
      
      Gradient1(1) = -ARG3 * ARG1 * RADIUS
      Gradient1(2) =  ARG3 * ARG2 * RADIUS
      Gradient1(3) =  0.0D0

      Gradient2(1) = ARG4 * ARG2 * RADIUS
      Gradient2(2) = ARG4 * ARG1 * RADIUS
      Gradient2(3) = - ARG3 * RADIUS
      
    END SUBROUTINE GRADMETRICSPHERE


!----------------------------------------------------------------
! Mobius Strip     
      SUBROUTINE GRADMETRICMOBIUS (X, Gradient1, Gradient2)

      USE COMMONS, ONLY : GThomsonC, GThomsonZ
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: Gradient1(3), Gradient2(3)
      DOUBLE PRECISION SINQ2, COSQ2, COSQ12, COSQ1, SINQ1, SINQ12
      DOUBLE PRECISION ARG1, ARG2, ARG3

      SINQ2 = SIN(X(2))
      COSQ2 = COS(X(2))
      COSQ1 = COS(X(1))
      SINQ1 = SIN(X(1))
      SINQ12 = SIN(X(1)/2.0D0)
      COSQ12 = COS(X(1)/2.0D0)
      ARG1 = GThomsonC + GThomsonZ * COSQ2 * COSQ12
      ARG2 = GThomsonZ * COSQ2 * SINQ12 / 2.0D0
      ARG3 = GThomsonZ * SINQ2 * COSQ12

      Gradient1(1) = -ARG1*SINQ1 - ARG2 * COSQ1
      Gradient1(2) =  ARG1*COSQ1 - ARG2 * SINQ1
      Gradient1(3) =  GThomsonZ * COSQ2 * COSQ12 /2.0D0

      Gradient2(1) = -ARG3 * COSQ1
      Gradient2(2) = -ARG3 * SINQ1
      Gradient2(3) = -GThomsonZ * SINQ2 * SINQ12

    END SUBROUTINE GRADMETRICMOBIUS
      

!----------------------------------------------------------------
! Check limits of the unduloids      
      SUBROUTINE FINDNGZ()

      USE COMMONS, ONLY : GThomsonC, GThomsonC2, GThomsonZ, GTHOMMET
      IMPLICIT NONE
      DOUBLE PRECISION :: c, a, pi, mu, k, Felint, Selint, u
      DOUBLE PRECISION :: refZ, refU, TEMPZ, GZMIN, GZMAX, ZPOS
      LOGICAL :: TEST

      TEMPZ = GThomsonZ
      GZMAX = 1.0D0
      GZMIN = 0.0D0
      TEST = .FALSE.

      pi = 4.0D0*ATAN(1.0D0)
      IF (GTHOMMET .EQ. 3) THEN
         refU = 1.50D0*pi
      ELSEIF (GTHOMMET .EQ. 4) THEN
         refU = 0.50D0*pi
      ENDIF

      IF (GThomsonC > GThomsonC2) THEN
         c = GThomsonC
         a = GThomsonC2
      ELSE
         c = GThomsonC2
         a = GThomsonC
      ENDIF
      mu = 2.0D0/(a+c)
      k = SQRT(1-(a/c)**2)
!      m = (c**2-a**2)/2.0D0
!      n = (c**2+a**2)/2.0D0      

      u = refU/mu + GZMAX*pi/mu
      CALL EllipIntegral(Felint, Selint, (mu*u/2.0D0-pi/4.0D0), k)
      ZPOS = a*Felint + c*Selint
      IF (GTHOMMET .EQ. 3) THEN
         CALL EllipIntegral(Felint, Selint, (refU/2.0D0-pi/4.0D0), k)
         refZ = a*Felint + c*Selint
         ZPOS = ZPOS -refZ
      ENDIF

      IF ( (ZPOS < TEMPZ) .OR. (ZPOS < 0) ) THEN
         PRINT *, "ERROR in the definition of unduloids"
         STOP
      ENDIF

      DO WHILE ( (.NOT. TEST) )
         GThomsonZ = (GZMAX + GZMIN)/2.0D0
         u    = refU/mu + GThomsonZ*pi/mu
         CALL EllipIntegral(Felint, Selint, (mu*u/2.0D0-pi/4.0D0), k)
         ZPOS = a*Felint + c*Selint
         IF (GTHOMMET .EQ. 3) THEN
            CALL EllipIntegral(Felint, Selint, (refU/2.0D0-pi/4.0D0), k)
            refZ = a*Felint + c*Selint
            ZPOS = ZPOS -refZ
         ENDIF

         IF (TEMPZ > ZPOS) GZMIN = GThomsonZ
         IF (TEMPZ < ZPOS) GZMAX = GThomsonZ

         IF ( (TEMPZ - ZPOS < 0.00001) .AND. (TEMPZ - ZPOS > -0.00001) ) THEN
            TEST = .TRUE.
            !PRINT *, GThomsonZ
         ENDIF
      ENDDO 
     
    END SUBROUTINE FINDNGZ


!----------------------------------------------------------------
! Check limits of the unduloids      
      SUBROUTINE CONVERTUNDULOIDPARAMETERS(V)

      USE COMMONS, ONLY : GThomsonC, GThomsonC2, GThomsonZ, GTHOMMET
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: V
      DOUBLE PRECISION :: c, a, pi, mu, k, m, n, Felint, Selint, u
      DOUBLE PRECISION :: refZ, refU, TEMPZ, GZMIN, GZMAX, ZPOS
      DOUBLE PRECISION :: amax, amin, cmax, cmin
      DOUBLE PRECISION :: Volume
      LOGICAL :: TEST, TEST2


      TEMPZ = GThomsonZ
      TEST2 = .FALSE.

      pi = 4.0D0*ATAN(1.0D0)
      IF (GTHOMMET .EQ. 3) THEN
         refU = 1.50D0*pi
      ELSEIF (GTHOMMET .EQ. 4) THEN
         refU = 0.50D0*pi
      ENDIF

      IF (GThomsonC > GThomsonC2) THEN
         c = GThomsonC
         a = GThomsonC2
      ELSE
         c = GThomsonC2
         a = GThomsonC
      ENDIF
      cmax = c
      cmin = a
      amax = c
      amin = a
 
      DO WHILE ( (.NOT. TEST2) )
         
         IF (GTHOMMET .EQ. 3) THEN
            c = (cmax + cmin)/2.0D0
         ELSE IF (GTHOMMET .EQ. 4) THEN
            a = (amax + amin)/2.0D0
         ENDIF

         GZMAX = 1.0D0
         GZMIN = 0.0D0
         TEST = .FALSE.
         mu = 2.0D0/(a+c)
         k = SQRT(1-(a/c)**2)
         m = (c**2-a**2)/2.0D0
         n = (c**2+a**2)/2.0D0      
                  
         DO WHILE ( (.NOT. TEST) )
            GThomsonZ = (GZMAX + GZMIN)/2.0D0
            u    = refU/mu + GThomsonZ*pi/mu
            CALL EllipIntegral(Felint, Selint, (mu*u/2.0D0-pi/4.0D0), k)
            ZPOS = a*Felint + c*Selint
            IF (GTHOMMET .EQ. 3) THEN
               CALL EllipIntegral(Felint, Selint, (refU/2.0D0-pi/4.0D0), k)
               refZ = a*Felint + c*Selint
               ZPOS = ZPOS -refZ
            ENDIF
            
            IF (TEMPZ > ZPOS) GZMIN = GThomsonZ
            IF (TEMPZ < ZPOS) GZMAX = GThomsonZ
            
            IF ( (TEMPZ - ZPOS < 0.00001) .AND. (TEMPZ - ZPOS > -0.00001) ) THEN
               TEST = .TRUE.
            ENDIF

         ENDDO

         IF (GTHOMMET .EQ. 3) THEN

            CALL EllipIntegral(Felint, Selint, (mu*u/2.0D0-pi/4.0D0), k)            
            Volume = pi * (2.0D0*(a**2+c**2)*c+3.0D0*a*c**2)/3.0D0 * Selint - pi*a**2*c/3.0D0*Felint &
                 - pi*(c**2-a**2)/6.0D0*SQRT(m*SIN(mu*u)+n)*COS(mu*u)
            CALL EllipIntegral(Felint, Selint, (refU/2.0D0-pi/4.0D0), k)
            Volume = Volume - ( pi * (2.0D0*(a**2+c**2)*c+3.0D0*a*c**2)/3.0D0 * Selint - pi*a**2*c/3.0D0*Felint &
                 - pi*(c**2-a**2)/6.0D0*SQRT(m*SIN(refU)+n)*COS(refU) )            
            IF (Volume > V) THEN 
               cmax = c
            ELSEIF (Volume < V) THEN
               cmin = c
            ENDIF

         ELSEIF (GTHOMMET .EQ. 4) THEN

            CALL EllipIntegral(Felint, Selint, (mu*u/2.0D0-pi/4.0D0), k)            
            Volume = pi * (2.0D0*(a**2+c**2)*c+3.0D0*a*c**2)/3.0D0 * Selint - pi*a**2*c/3.0D0*Felint &
                 - pi*(c**2-a**2)/6.0D0*SQRT(m*SIN(mu*u)+n)*COS(mu*u)
            CALL EllipIntegral(Felint, Selint, (refU/2.0D0-pi/4.0D0), k)
            Volume = Volume - ( pi * (2.0D0*(a**2+c**2)*c+3.0D0*a*c**2)/3.0D0 * Selint - pi*a**2*c/3.0D0*Felint &
                 - pi*(c**2-a**2)/6.0D0*SQRT(m*SIN(refU)+n)*COS(refU) )            

            IF (Volume > V) THEN 
               amax = a
            ELSEIF (Volume < V) THEN
               amin = a
            ENDIF

         ENDIF

         IF ( (Volume - V < 0.00001) .AND. (Volume - V > -0.00001) ) THEN
            TEST2 = .TRUE.
         ENDIF
         
      ENDDO

      PRINT *, "GTHOMSONC ", c
      PRINT *, "GTHOMSONA ", a
      PRINT *, "umax ", GTHOMSONZ

      GThomsonC = c
      GThomsonC2 = a
      GThomsonZ = TEMPZ
         
      PRINT *, "zmax ", GTHOMSONZ

    END SUBROUTINE CONVERTUNDULOIDPARAMETERS


!----------------------------------------------------------------

SUBROUTINE GTHOMSONHESSIAN(X,HESS)

  USE commons, ONLY: GTHOMMET, GTHOMPOT, GThomsonSigma, NATOMS
  IMPLICIT NONE

  INTEGER J1, J2, J3, J4
  DOUBLE PRECISION, INTENT(IN)  :: X(*)
  DOUBLE PRECISION, INTENT(OUT) :: HESS(3*NATOMS, 3*NATOMS)
  DOUBLE PRECISION DIST, D2VDR2, DVDR
  DOUBLE PRECISION :: TMPCOORDS(3*NATOMS), DR(3), Gradient11(3), Gradient12(3), Gradient21(3), Gradient22(3) 
  DOUBLE PRECISION :: HESS111(3), HESS112(3), HESS122(3)
  DOUBLE PRECISION :: HESS211(3), HESS212(3), HESS222(3)

  CALL GTHOMSONANGTOC(TMPCOORDS(1:3*NATOMS), X(1:3*NATOMS), NATOMS)      
  HESS(1:3*NATOMS,1:3*NATOMS)=0.0D0
  
  DO J1=1,NATOMS
     J3=3*J1

     IF ( GTHOMMET .EQ. 1) THEN
        PRINT *, "METRIC NOT IMPLEMENTED IN HESSIAN"
        STOP
     ELSEIF ( GTHOMMET .EQ. 2) THEN
        CALL GRADMETRICCATENOID (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
        CALL HESSMETRICCATENOID (X(3*J1-2:3*J1-1), HESS111, HESS112, HESS122)
     ELSEIF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN
        PRINT *, "METRIC NOT IMPLEMENTED IN HESSIAN"
        STOP
     ELSEIF ( GTHOMMET .EQ. 5 ) THEN
        CALL GRADMETRICSPHERE (X(3*J1-2:3*J1-1), Gradient11, Gradient12)
        CALL HESSMETRICSPHERE (X(3*J1-2:3*J1-1), HESS111, HESS112, HESS122)
     ENDIF
     
     DO J2=J1+1,NATOMS
        J4=3*J2
        DR(1) = TMPCOORDS(J3-2)-TMPCOORDS(J4-2)
        DR(2) = TMPCOORDS(J3-1)-TMPCOORDS(J4-1)
        DR(3) = TMPCOORDS(J3  )-TMPCOORDS(J4  )
        DIST  = SQRT(DOT_PRODUCT(DR, DR))

        IF (GTHOMPOT .EQ. 1) THEN
           D2VDR2 =   3/DIST**5 
           DVDR   = - 1/DIST**3
        ELSEIF (GTHOMPOT .EQ. 2) THEN 
           PRINT *, "POTENTIAL NOT IMPLEMENTED IN HESSIAN"
           STOP
        ELSEIF (GTHOMPOT .EQ. 3) THEN
           D2VDR2 =  ( 3/DIST**5 + 3/GThomsonSigma/DIST**4 + 1/GThomsonSigma**2/DIST**3 ) * EXP(-1.0D0*DIST/GThomsonSigma)
           DVDR   = -( 1/DIST**3 + 1/GThomsonSigma/DIST**2 ) * EXP(-1.0D0*DIST/GThomsonSigma)
        ELSEIF (GTHOMPOT .EQ. 4) THEN
           PRINT *, "POTENTIAL NOT IMPLEMENTED IN HESSIAN"
           STOP
        ELSEIF (GTHOMPOT .EQ. 5) THEN
           PRINT *, "POTENTIAL NOT IMPLEMENTED IN HESSIAN"
           STOP
        ELSEIF (GTHOMPOT .EQ. 6) THEN
           PRINT *, "POTENTIAL NOT IMPLEMENTED IN HESSIAN"
           STOP
        ENDIF
        
        HESS(3*J1-2,3*J1-2) = HESS(3*J1-2,3*J1-2) +  D2VDR2 * (DOT_PRODUCT(DR,Gradient11))**2 &
             + DVDR * (DOT_PRODUCT(Gradient11,Gradient11)+ DOT_PRODUCT(DR,HESS111)) 
        HESS(3*J1-1,3*J1-1) = HESS(3*J1-1,3*J1-1) +  D2VDR2 * (DOT_PRODUCT(DR,Gradient12))**2 &
             + DVDR * (DOT_PRODUCT(Gradient12,Gradient12)+ DOT_PRODUCT(DR,HESS122)) 
        HESS(3*J1-2,3*J1-1) = HESS(3*J1-2,3*J1-1) +  D2VDR2 * DOT_PRODUCT(DR,Gradient11) * DOT_PRODUCT(DR,Gradient12) &
             + DVDR * (DOT_PRODUCT(Gradient11,Gradient12) + DOT_PRODUCT(DR,HESS112)) 

        IF ( GTHOMMET .EQ. 1) THEN
        ELSEIF ( GTHOMMET .EQ. 2) THEN
           CALL GRADMETRICCATENOID (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
           CALL HESSMETRICCATENOID (X(3*J2-2:3*J2-1), HESS211, HESS212, HESS222)
        ELSEIF ( (GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4) ) THEN
        ELSEIF ( GTHOMMET .EQ. 5 ) THEN
           CALL GRADMETRICSPHERE (X(3*J2-2:3*J2-1), Gradient21, Gradient22)
           CALL HESSMETRICSPHERE (X(3*J2-2:3*J2-1), HESS211, HESS212, HESS222)
        ENDIF

        HESS(3*J2-2,3*J2-2) = HESS(3*J2-2,3*J2-2) +  D2VDR2 * (DOT_PRODUCT(DR,Gradient21))**2 &
             + DVDR * (DOT_PRODUCT(Gradient21,Gradient21) - DOT_PRODUCT(DR,HESS211)) 
        HESS(3*J2-1,3*J2-1) = HESS(3*J2-1,3*J2-1) +  D2VDR2 * (DOT_PRODUCT(DR,Gradient22))**2 &
             + DVDR * (DOT_PRODUCT(Gradient22,Gradient22) - DOT_PRODUCT(DR,HESS222)) 
        HESS(3*J2-2,3*J2-1) = HESS(3*J2-2,3*J2-1) +  D2VDR2 * DOT_PRODUCT(DR,Gradient21) * DOT_PRODUCT(DR,Gradient22) &
             + DVDR * (DOT_PRODUCT(Gradient11,Gradient12) - DOT_PRODUCT(DR,HESS212)) 
        
        HESS(3*J1-2,3*J2-2) =  HESS(3*J1-2,3*J2-2) -D2VDR2 * DOT_PRODUCT(DR,Gradient11) * DOT_PRODUCT(DR,Gradient21) &
                               - DVDR * DOT_PRODUCT(Gradient11,Gradient21)          
        HESS(3*J1-2,3*J2-1) =  HESS(3*J1-2,3*J2-1) -D2VDR2 * DOT_PRODUCT(DR,Gradient11) * DOT_PRODUCT(DR,Gradient22) &
                               - DVDR * DOT_PRODUCT(Gradient11,Gradient22)
        HESS(3*J1-1,3*J2-2) =  HESS(3*J1-1,3*J2-2) -D2VDR2 * DOT_PRODUCT(DR,Gradient12) * DOT_PRODUCT(DR,Gradient21) &
                               - DVDR * DOT_PRODUCT(Gradient12,Gradient21)
        HESS(3*J1-1,3*J2-1) =  HESS(3*J1-1,3*J2-1) -D2VDR2 * DOT_PRODUCT(DR,Gradient12) * DOT_PRODUCT(DR,Gradient22) &
                               - DVDR * DOT_PRODUCT(Gradient12,Gradient22)

     ENDDO
  ENDDO

  DO J1=1,3*NATOMS
     DO J2 = J1+1,3*NATOMS
        HESS(J2,J1) = HESS(J1,J2)
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE GTHOMSONHESSIAN

!----------------------------------------------------------------

SUBROUTINE GTHOMSONNUMHESSIAN(X,HESS)

  USE commons, ONLY : NATOMS
  IMPLICIT NONE

  INTEGER J1, J2
  DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: HESS(3*NATOMS,3*NATOMS)
  DOUBLE PRECISION :: X1(3*NATOMS), X2(3*NATOMS), DX, ETHOMSON
  DOUBLE PRECISION :: V1(3*NATOMS), V2(3*NATOMS)
  
  DX = 0.00001

  DO J1=1,3*NATOMS
     X1(:) = X(:)
     X1(J1) = X1(J1) - DX
     CALL GTHOMSON(X1,V1,ETHOMSON,.TRUE.)
     X2(:) = X(:)
     X2(J1) = X2(J1) + DX
     CALL GTHOMSON(X2,V2,ETHOMSON,.TRUE.)
     
     DO J2=J1,3*NATOMS
        HESS(J1,J2) = (V2(J2)-V1(J2))/(2.0D0*DX)
        HESS(J2,J1) = HESS(J1,J2)
     END DO

  END DO
  
  RETURN
  
END SUBROUTINE GTHOMSONNUMHESSIAN

!----------------------------------------------------------------
! catenoid      
      SUBROUTINE HESSMETRICCATENOID (X, HESS11, HESS12, HESS22)

      USE COMMONS, ONLY : GThomsonC, GThomsonZ
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: X(2)
      DOUBLE PRECISION, INTENT(OUT) :: HESS11(3), HESS12(3), HESS22(3)
      DOUBLE PRECISION ARGHY, COSHA, SINHA, SINX2, SINX1, COSX2, COSX1

      COSX2 = COS(X(2))
      ARGHY = GThomsonZ/GThomsonC*COSX2
      COSHA = COSH(ARGHY)
      SINHA = SINH(ARGHY)
      SINX2 = SIN(X(2))
      SINX1 = SIN(X(1))
      COSX1 = COS(X(1))
      
      HESS11(1) = -GThomsonC * COSHA * COSX1
      HESS11(2) = -GThomsonC * COSHA * SINX1
      HESS11(3) =  0.0D0

      HESS12(1) =  GThomsonZ * SINHA * SINX1 * SINX2
      HESS12(2) = -GThomsonZ * SINHA * COSX1 * SINX2
      HESS12(3) =  0.0D0

      HESS22(1) = -GThomsonZ * COSX1 * (COSX2 * SINHA - GThomsonZ/GThomsonC *SINX2**2 *COSHA)
      HESS22(2) = -GThomsonZ * SINX1 * (COSX2 * SINHA - GThomsonZ/GThomsonC *SINX2**2 *COSHA)
      HESS22(3) = -GThomsonC * ARGHY

    END SUBROUTINE HESSMETRICCATENOID


!----------------------------------------------------------------
! sphere
      SUBROUTINE HESSMETRICSPHERE (X, HESS11, HESS12, HESS22)

        USE COMMONS, ONLY : GThomsonZ, GTHOMMET
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN)  :: X(2)
        DOUBLE PRECISION, INTENT(OUT) :: HESS11(3), HESS12(3), HESS22(3)
        DOUBLE PRECISION ARG1, ARG2, ARG3, ARG4
        DOUBLE PRECISION RADIUS
        
        IF ( GTHOMMET .EQ. 5 ) THEN
           RADIUS = GTHOMSONZ
        ENDIF

        ARG1 = SIN(X(1))
        ARG2 = COS(X(1))
        ARG3 = SIN(X(2))
        ARG4 = COS(X(2))

        HESS11(1) = -ARG3 * ARG2 * RADIUS
        HESS11(2) = -ARG3 * ARG1 * RADIUS
        HESS11(3) = 0.0

        HESS12(1) = -ARG4 * ARG1 * RADIUS
        HESS12(2) =  ARG4 * ARG2 * RADIUS
        HESS12(3) = 0.0
        
        HESS22(1) = -ARG3 * ARG1 * RADIUS
        HESS22(2) = -ARG3 * ARG2 * RADIUS
        HESS22(3) = -ARG4 * RADIUS
        
      END SUBROUTINE HESSMETRICSPHERE

!----------------------------------------------------------------------------------------------------------------

SUBROUTINE HKMINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RMATBEST)
USE COMMONS, ONLY: BESTPERM, EFIELDT, GEOMDIFFTOL, GTHOMMET, NFREEZE, NPERMGROUP, NPERMSIZE, &
                 & NSETS, PERMDIST, PERMGROUP, PULLT, SETS, STOCKT, MYUNIT
IMPLICIT NONE

INTEGER :: MAXIMUMTRIES=20
INTEGER NATOMS, NPERM, PATOMS, NTRIES
INTEGER J3, INVERT, NORBIT1, NORBIT2, NCHOOSE2, NDUMMY, LPERM(NATOMS), J1, J2, NCHOOSE1, NROTDONE, NORBITB1, NORBITB2, &
  &     NCHOOSEB1, NCHOOSEB2
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), &
  &              DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX, BOXLY, BOXLZ, WORSTRAD, RMAT(3,3), DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3)
DOUBLE PRECISION RMATCUMUL(3,3)
DOUBLE PRECISION REFXZ(3,3)
LOGICAL DEBUG, TWOD, BULKT
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)
CHARACTER(LEN=5) ZSYMSAVE
COMMON /SYS/ ZSYMSAVE

NROTDONE=-1
!
REFXZ(1:3,1:3)=0.0D0
REFXZ(1,1)=1.0D0; REFXZ(2,2)=-1.0D0; REFXZ(3,3)=1.0D0
!
! It is possible for the standard orientation to result in a distance that is worse than
! the starting distance. Hence we need to set XBEST here.
!
DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)
DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)
DBEST = 1.0D100

XDUMMY = 0.0D0
XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
DO J1 = 1,NATOMS
   XDUMMY = XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
&                  (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
&                  (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
ENDDO
DBEST = XDUMMY

DO J1=1,NATOMS
   BESTPERM(1:NATOMS)=J1
ENDDO
RMATBEST(1:3,1:3)=RMAT(1:3,1:3)
ROTINVBBEST(1:3,1:3)=0.0D0
ROTINVBBEST(1,1)=1.0D0;ROTINVBBEST(2,2)=1.0D0;ROTINVBBEST(3,3)=1.0D0;
ROTABEST(1:3,1:3)=0.0D0
ROTABEST(1,1)=1.0D0;ROTABEST(2,2)=1.0D0;ROTABEST(3,3)=1.0D0;

!
! End of XBEST associated initialisation.
!
IF (.NOT. PERMDIST) THEN
   CALL GTHOMSONNEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,.FALSE.,DEBUG,RMAT)
   DISTANCE = DISTANCE**2
   IF (DISTANCE .LT. DBEST) THEN
     IF (DEBUG) PRINT *, 'Reducing best distance squared from ', DBEST, ' to ', DISTANCE 
     DBEST=DISTANCE
     XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
     RMATBEST(1:3,1:3)=RMAT(1:3,1:3)
     RMATBEST=MATMUL(RMATBEST,ROTABEST)
   END IF
   GOTO 50
END IF


NROTDONE=-1
11 CONTINUE
NROTDONE=NROTDONE+1

INVERT=1
60 CONTINUE ! jump back here if INVERT changes sign.
   NCHOOSEB1=0
66 NCHOOSEB1=NCHOOSEB1+1
   NCHOOSEB2=0
31 NCHOOSEB2=NCHOOSEB2+1
   NCHOOSE1=0
65 NCHOOSE1=NCHOOSE1+1
   NCHOOSE2=0
30 NCHOOSE2=NCHOOSE2+1
DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)

DO J1=1,NATOMS
   ALLPERM(J1)=J1
ENDDO

! The optimal alignment returned by minpermdist is a local minimum, but may not
! be the global minimum. Calling MYORIENT first should put permutational isomers
! into a standard alignment and spot the global minimum zero distance in one
! go. However, we also need to cycle over equivalent atoms in orbits using NCHOOSE2.
!
! Problems can occur if we don't use all the atoms specified by NORBIT1 and NORBIT2
! because of the numerical cutoffs employed in MYORIENT. We could miss the
! right orientation! 
!
! If we use MYORIENT to produce particular orientations then we end up aligning 
! COORDSA not with COORDSB but with the standard orientation of COORDSB in DUMMYB.
! We now deal with this by tracking the complete transformation, including the
! contribution of MYORIENT using ROTB and ROTINVB.
!

DISTANCE=0.0D0
IF (NFREEZE.LE.0) THEN
      IF ((PULLT.OR.EFIELDT.OR.TWOD.OR.(GTHOMMET < 5)).AND.(INVERT.EQ.-1)) THEN ! reflect in xz plane
         DO J1=1,NATOMS
            DUMMYC(3*(J1-1)+1)=DUMMYA(3*(J1-1)+1)
            DUMMYC(3*(J1-1)+2)=-DUMMYA(3*(J1-1)+2)
            DUMMYC(3*(J1-1)+3)=DUMMYA(3*(J1-1)+3)
         ENDDO
      ELSE
         DUMMYC(1:3*NATOMS)=INVERT*DUMMYA(1:3*NATOMS)
      ENDIF 

      CALL MYORIENT(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)
      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      CALL MYORIENT(DUMMYB,DUMMY,NORBITB1,NCHOOSEB1,NORBITB2,NCHOOSEB2,NATOMS,DEBUG,ROTB,ROTINVB,STOCKT)
      DUMMYB(1:3*NATOMS)=DUMMY(1:3*NATOMS)

      DISTANCE=0.0D0
      DO J1=1,3*NATOMS
         DISTANCE=DISTANCE+(DUMMYA(J1)-DUMMYB(J1))**2
      ENDDO
!      IF (DEBUG) PRINT *, ' minpermdist> after initial call to MYORIENT distance squared = ', DISTANCE, ' for ',NATOMS,' atoms'
!  IF (DEBUG) PRINT '(A,6I8)',' minpermdist> size of orbits, selected atoms, random rotations, invert: ', &
! &       NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2,NROTDONE,INVERT
ELSE
   NORBIT1=1; NORBIT2=1; NORBITB1=1; NORBITB2=1
   CALL GTHOMSONNEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,.FALSE.,DEBUG,RMAT)
   IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> after initial call to NEWMINDIST distance=',DISTANCE
   DISTANCE=DISTANCE**2
ENDIF

!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case.
!  We return to label 10 after every round of permutational/orientational alignment
!  unless we have converged to the identity permutation.
!
!  Atoms are not allowed to appear in more than one group.
!  The maximum number of pair exchanges associated with a group is two.
!
NTRIES=0
!
!  RMATCUMUL contains the accumulated rotation matrix that relates the original 
!  DUMMYA obtained from COORDSA to the final one.
!
RMATCUMUL(1:3,1:3)=0.0D0
RMATCUMUL(1,1)=1.0D0; RMATCUMUL(2,2)=1.0D0; RMATCUMUL(3,3)=1.0D0
10 CONTINUE

NTRIES=NTRIES+1

NDUMMY=1
DO J1=1,NATOMS
   NEWPERM(J1)=J1
ENDDO

!
! ALLPERM saves the permutation from the previous cycle.
! NEWPERM contains the permutation for this cycle, relative to the identity.
! SAVEPERM is temporary storage for NEWPERM.
! NEWPERM must be applied to ALLPERM after the loop over NPERMGROUP and
! corresponding swaps.
!
! New version allows for overlapping atoms in NPERMGROUP, so that atoms
! can appear in more than one group. This was needed for flexible water potentials.
!
DO J1=1,NPERMGROUP
   PATOMS=NPERMSIZE(J1)
   DO J2=1,PATOMS
      PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
      PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
   ENDDO
!
! All permutations within this group of size NPERMSIZE(J1) are now tried.
!
   CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)
   SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
   DO J2=1,PATOMS
      SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1))
   ENDDO
!
! Update permutation of associated atoms, if any. 
! We must do this as we go along, because these atoms could move in more than
! one permutational group now.
!
   IF (NSETS(J1).GT.0) THEN
      DO J2=1,PATOMS
         DO J3=1,NSETS(J1)
            SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J3)
         ENDDO
      ENDDO
   ENDIF
   NDUMMY=NDUMMY+NPERMSIZE(J1)
   NEWPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
! PRINT '(A)',' minpermdist> NEWPERM:'
! PRINT '(20I6)',NEWPERM(1:NATOMS)

ENDDO

!
! Update the overall permutation here.
! The latest NEWPERM(J1) tells us which position moves to J1 in the latest
! permutation, relative to the identity.
! ALLPERM(J2) tells us which atom has moved to position J2.
! So, the new overall permutation, i.e. the atoms that moves to position J1
! After ALLPERM followed by NEWPERM is ALLPERM(NEWPERM(J1))
!
DO J1=1,NATOMS
!  SAVEPERM(ALLPERM(J1))=ALLPERM(NEWPERM(J1)) !!! BUG 12/9/11 DJW
   SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
ENDDO
ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
! PRINT '(A)',' minpermdist> ALLPERM:'
! PRINT '(20I6)',ALLPERM(1:NATOMS)

DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
DISTANCE=0.0D0
!
! Update coordinates in DUMMYA to overall permutation using NEWPERM.
!
DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)

   IF (J3.NE.NEWPERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' minpermdist> move position ',NEWPERM(J3),' to ',J3
      NPERM=NPERM+1
   ENDIF
   DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))**2 &
  &                    +(DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))**2 &
  &                    +(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))**2
ENDDO
!
!  Optimal alignment. Coordinates in DUMMYA are reset by NEWMINDIST (second argument).
!  Must allow at least one call to NEWMINDIST in case the MYORIENT result is terrible
!  but gives zero permutations!
!  
 
! PRINT '(A,I6,2G20.10)','NPERM,DBEST,DISTANCE=',NPERM,DBEST,DISTANCE

IF ((NPERM.NE.0).OR.(NTRIES.EQ.1)) THEN 
   CALL GTHOMSONNEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,.FALSE.,DEBUG,RMAT)
   RMATCUMUL=MATMUL(RMAT,RMATCUMUL)
   DISTANCE=DISTANCE**2 ! we are using DISTANCE^2 further down
!  IF (DEBUG) WRITE(*,'(A,G20.10)') ' minpermdist> distance after NEWMINDIST=                     ', &
! &                                    SQRT(DISTANCE) 
   IF (NTRIES.LT.MAXIMUMTRIES) THEN
      GOTO 10
   ELSE ! prevent infinite loop
      IF(DEBUG) PRINT '(A)',' minpermdistGTHOMSON> WARNING - number of tries exceeded, giving up'
   ENDIF
ENDIF

IF (DISTANCE.LT.DBEST) THEN
   IF (DEBUG) PRINT *, 'Reducing best distance squared from ', DBEST, ' to ', DISTANCE 
   DBEST=DISTANCE
   XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
   BESTPERM(1:NATOMS)=ALLPERM(1:NATOMS)
   RMATBEST(1:3,1:3)=RMATCUMUL(1:3,1:3)
   ROTINVBBEST(1:3,1:3)=ROTINVB(1:3,1:3) 
   ROTABEST(1:3,1:3)=ROTA(1:3,1:3)      
   RMATBEST=MATMUL(RMATBEST,ROTABEST)
   IF (INVERT.EQ.-1) THEN
      IF (PULLT.OR.EFIELDT.OR.TWOD.OR.(GTHOMMET < 5)) THEN ! reflect in xz plane rather than invert!
         RMATBEST(1:3,1:3)=MATMUL(RMATBEST,REFXZ)
      ELSE
         RMATBEST(1:3,1:3)=-RMATBEST(1:3,1:3)
      ENDIF
   ENDIF

ENDIF

IF (DISTANCE .LT. GEOMDIFFTOL ) THEN
   GOTO 50 ! hk286
ENDIF

IF (NCHOOSE2.LT.NORBIT2) GOTO 30
IF (NCHOOSE1.LT.NORBIT1) GOTO 65
IF (NCHOOSEB2.LT.NORBITB2) GOTO 31
IF (NCHOOSEB1.LT.NORBITB1) GOTO 66
!
!  Now try the enantiomer (or xz reflected structure for PULLT.OR.EFIELDT.OR.TWOD).
!
!  GO TO 50
IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
! don't try inversion for bulk or charmm or amber or frozen atoms
!
   INVERT=-1
   GOTO 60
ENDIF

50 DISTANCE=DBEST
!
!  XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!  Rotate XBEST by ROTINVBBEST to put in best correspondence with COORDSB, 
!  undoing the reorientation to DUMMYB from MYORIENT. 
!  We should get the same result for ROTINVBBEST * RMATBEST * (COORDSA-CMA) 
!  where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
!  (aside from a possible permutation of the atom ordering)
!

   XDUMMY=0.0D0
   DO J1=1,NATOMS
      XBEST(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROTINVBBEST,XBEST(3*(J1-1)+1:3*(J1-1)+3))
      XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
   ENDDO
   IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
      PRINT '(2(A,F20.10))',' minpermdistGTHOMSON> ERROR *** distance between transformed XBEST and COORDSB=',SQRT(XDUMMY), &
  &                         ' should be ',SQRT(DISTANCE)
      STOP
   ENDIF

   RMATBEST=MATMUL(ROTINVBBEST,RMATBEST)
   COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!
   DISTANCE=SQRT(DISTANCE) ! now changed to return distance, not distance^2 22/11/10 DJW
   
RETURN
END SUBROUTINE HKMINPERMDIST

! -----------------------------------------------------------

SUBROUTINE HKMINDIST(RA,RB,NATOMS,DIST,BULKT,PRESERVET)
! hk286 - add GTHOMSONT,  NGTHORI
  USE COMMONS, ONLY : GTHOMSONT, GTHOMMET
  IMPLICIT NONE
  INTEGER J1, ITER, J2, NCOUNT
  DOUBLE PRECISION DPRAND
  DOUBLE PRECISION P(3), DIST, DIST0, DISTFUNC, & 
       &                 MYROTMAT(3,3), &
       &                 OMEGATOT(3,3), RA(*), RB(*)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R, R0, R1, R1SAVE
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RBSAVE
  INTEGER NSIZE, NATOMS
  LOGICAL BULKT, MFLAG, PRESERVET
  COMMON /GEOM/ NSIZE
  COMMON /MINDOM/ MYROTMAT, OMEGATOT

  ALLOCATE(RBSAVE(3*NATOMS))
  RBSAVE(1:3*NATOMS)=RB(1:3*NATOMS)
  !
  !  Initialise accumulated rotation matrix to the identity.
  !
  DO J1=1,3
     DO J2=1,3
        OMEGATOT(J2,J1)=0.0D0
     ENDDO
     OMEGATOT(J1,J1)=1.0D0
  ENDDO
  
  ALLOCATE(R(3,NATOMS),R0(3,NATOMS),R1(3,NATOMS),R1SAVE(3,NATOMS))
  NSIZE=NATOMS
  DO J1=1,NSIZE
     R(1,J1)=RA(3*(J1-1)+1)
     R(2,J1)=RA(3*(J1-1)+2)
     R(3,J1)=RA(3*(J1-1)+3)
     R0(1,J1)=RB(3*(J1-1)+1)
     R0(2,J1)=RB(3*(J1-1)+2)
     R0(3,J1)=RB(3*(J1-1)+3)
  ENDDO
  !
  !     initial angles
  !
  P(1)=0.0D0
  P(2)=0.0D0
  P(3)=0.0D0
  !     IF (TWOD) P(1)=0.0D0
  !     IF (TWOD) P(2)=0.0D0
  !
  !     calculate initial distance
  !
  NCOUNT=0
10 DIST0=DISTFUNC(P,R,R0,R1)
  DIST0=SQRT(DIST0)
  IF (BULKT) THEN
     DIST=DIST0
     IF (PRESERVET) RB(1:3*NATOMS)=RBSAVE(1:3*NATOMS)
     DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)
     RETURN
  ENDIF
  
  CALL MMYLBFGS(P,1.0D-7,MFLAG,DIST,500,ITER,R,R0,R1)
  DIST=DISTFUNC(P,R,R0,R1)
  DIST=SQRT(DIST)
  IF (MFLAG) THEN
     !        WRITE(*,'(A,2F15.5,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER
     !        PRINT*
  ELSE
     NCOUNT=NCOUNT+1
     IF (NCOUNT.GT.0) THEN 
        PRINT*,'convergence failure in mind'
        !           STOP
     ELSE
        !           WRITE(*,'(A,2F15.5,A,I6,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER,' NCOUNT=',NCOUNT
        DO J1=1,NSIZE
           R0(1,J1)=R1(1,J1)
           R0(2,J1)=R1(2,J1)
           R0(3,J1)=R1(3,J1)
        ENDDO
        ! hk286
        IF (.NOT. (GTHOMSONT .AND. (GTHOMMET < 5))) P(1)=2*(DPRAND()-0.5D0)/10.0D0
        IF (.NOT. (GTHOMSONT .AND. (GTHOMMET < 5))) P(2)=2*(DPRAND()-0.5D0)/10.0D0
        P(3)=2*(DPRAND()-0.5D0)/10.0D0
        GOTO 10
     ENDIF
  ENDIF
  DO J1=1,NSIZE
     RB(3*(J1-1)+1)=R1(1,J1)
     RB(3*(J1-1)+2)=R1(2,J1)
     RB(3*(J1-1)+3)=R1(3,J1)
     RA(3*(J1-1)+1)=R(1,J1)
     RA(3*(J1-1)+2)=R(2,J1)
     RA(3*(J1-1)+3)=R(3,J1)
  ENDDO
  
  DO J1=1,3
     DO J2=1,3
        MYROTMAT(J2,J1)=OMEGATOT(J2,J1)
     ENDDO
  ENDDO
  
  DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)
  
  RETURN
END SUBROUTINE HKMINDIST

! -------------------------------------------------------------

SUBROUTINE GTHOMSONNEWMINDIST(RA,RB,NATOMS,DIST,BULKT,PRESERVET,DEBUG,RMAT)
  USE COMMONS, ONLY : GTHOMSONT, GTHOMMET
  
  IMPLICIT NONE
  INTEGER J1, NATOMS, NSIZE, INFO, JMIN
  DOUBLE PRECISION RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), XM, YM, ZM, XP, YP, ZP, &
       &              DIAG(4), TEMPA(9*NATOMS), RMAT(3,3), MINV, Q1, Q2, Q3, Q4, CMXA, CMYA, CMZA, &
       &              MYROTMAT(3,3), OMEGATOT(3,3)
  DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
  LOGICAL BULKT, PRESERVET, DEBUG
  COMMON /MINDOM/ MYROTMAT, OMEGATOT

! jwrm2> Convert coordinates from polar to Cartesian
!  CALL GTHOMSONANGTOC(RA,COORDSA,NATOMS)
!  CALL GTHOMSONANGTOC(RB,COORDSB,NATOMS)

  IF (GTHOMSONT .AND. (GTHOMMET < 5)) THEN
  !  ALLOCATE(XA(3*(NATOMS/2)*number of sites,XB(3*(NATOMS/2)*number of sites))
  !  NSIZE=(NATOMS/2)*number of sites
  !  PRINT '(A)',' newmindist> New quaternion procedure not yet coded for flatland'
  ! There is one unknown angle, so this should be trivial!'
     CALL HKMINDIST(RA,RB,NATOMS,DIST,BULKT,PRESERVET)
     RMAT(1:3,1:3)=OMEGATOT(1:3,1:3)
     RETURN
     !  STOP
  ELSE
     ALLOCATE(XA(3*NATOMS),XB(3*NATOMS))
     NSIZE=NATOMS
     XA(1:3*NATOMS)=RA(1:3*NATOMS)
     XB(1:3*NATOMS)=RB(1:3*NATOMS)
  ENDIF  
  !
  ! Move centre of coordinates of XA and XB to the origin.
  !
  CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
  
  !
  !  The formula below is not invariant to overall translation because XP, YP, ZP
  !  involve a sum of coordinates! We need to have XA and XB coordinate centres both 
  !  at the origin!!
  !
  QMAT(1:4,1:4)=0.0D0
  DO J1=1,NSIZE
     XM=XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
     YM=XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
     ZM=XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
     XP=XA(3*(J1-1)+1)+XB(3*(J1-1)+1)
     YP=XA(3*(J1-1)+2)+XB(3*(J1-1)+2)
     ZP=XA(3*(J1-1)+3)+XB(3*(J1-1)+3)
     QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
     QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
     QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
     QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
     QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
     QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
     QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
     QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
     QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
     QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
  ENDDO
  QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3); QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)
  
  CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)
  IF (INFO.NE.0) PRINT '(A,I6,A)',' gthomsonnewmindist> WARNING - INFO=',INFO,' in DSYEV'
  
  MINV=1.0D100
  DO J1=1,4
     !     PRINT '(A,I8,G20.10)','newmindist> J1,DIAG=',J1,DIAG(J1)
     IF (DIAG(J1).LT.MINV) THEN
        JMIN=J1
        MINV=DIAG(J1)
     ENDIF
  ENDDO
  IF (MINV.LT.0.0D0) THEN
     IF (ABS(MINV).LT.1.0D-6) THEN
        MINV=0.0D0
     ELSE
        PRINT '(A,G20.10,A)',' gthomsonnewmindist> WARNING MINV is ',MINV,' change to absolute value'
        MINV=-MINV
     ENDIF
  ENDIF
  DIST=SQRT(MINV)
  
  IF (DEBUG) PRINT '(A,G20.10,A,I6)',' newmindist> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
  Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
  !
  ! RMAT will contain the matrix that maps XB onto the best correspondence with XA
  !
  RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
  RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
  RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
  RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
  RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
  RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
  RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
  RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
  RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2

IF (.NOT.PRESERVET) THEN
   CALL NEWROTGEOM(NSIZE,RB,RMAT,CMXA,CMYA,CMZA)
ENDIF
  
  DEALLOCATE(XA,XB)
  
  RETURN

END SUBROUTINE GTHOMSONNEWMINDIST

!----------------------------------------------------------------------------------------------------------------------

SUBROUTINE GTHOMSONMINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,DISTANCE,DIST2,RMATBEST)
  
  USE COMMONS, ONLY : GTHOMMET
  IMPLICIT NONE
  INTEGER NATOMS, J1
  DOUBLE PRECISION XCOORDSA (3*NATOMS), XCOORDSB(3*NATOMS), DD, DD2
  DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE
  DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,RMATBEST(3,3),RMATBEST2(3,3),REFXY(3,3)
  LOGICAL DEBUG, BULKT

  REFXY(:,:) = 0.0D0
  REFXY(1,1) = 1.0D0; REFXY(2,2) = 1.0D0; REFXY(3,3) = -1.0D0

  IF (GTHOMMET .EQ. 5) THEN
    CALL GTHOMSONANGTOC(XCOORDSA,COORDSA,NATOMS)
    CALL GTHOMSONANGTOC(XCOORDSB,COORDSB,NATOMS)
  ELSE
    XCOORDSA = COORDSA
    XCOORDSB = COORDSB
  END IF

  DISTANCE = 1.0D10
     
  CALL HKMINPERMDIST(XCOORDSB,XCOORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,.FALSE.,DD,DD2,RMATBEST)

!  IF (DD < DISTANCE .AND. DD < GEOMDIFFTOL) THEN
  IF (DD < DISTANCE) THEN
     DISTANCE = DD
     DIST2 = DD2
     CALL GTHOMSONCTOANG(XCOORDSA,COORDSA,NATOMS)
     CALL GTHOMSONCTOANG(XCOORDSB,COORDSB,NATOMS)
  ENDIF

  IF (GTHOMMET < 5) THEN
     DO J1 = 1, NATOMS
        XCOORDSA(3*J1) = -XCOORDSA(3*J1)
     ENDDO
     CALL HKMINPERMDIST(XCOORDSB,XCOORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,.FALSE.,DD,DD2,RMATBEST2)
     !  IF (DD < DISTANCE .AND. DD < GEOMDIFFTOL) THEN
     IF (DD < DISTANCE) THEN
        DISTANCE = DD
        DIST2 = DD2
        CALL GTHOMSONCTOANG(XCOORDSA,COORDSA,NATOMS)
        CALL GTHOMSONCTOANG(XCOORDSB,COORDSB,NATOMS)
        RMATBEST(:,:)=MATMUL(REFXY,RMATBEST)
        RMATBEST(:,:)=MATMUL(RMATBEST2,RMATBEST)     
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE GTHOMSONMINPERMDIST
