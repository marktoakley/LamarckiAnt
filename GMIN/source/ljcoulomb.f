C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C  Energy and gradient for LJ plus coulomb repulsion on some particles.
C  Mark Miller and Marie-Pierre Gaigeot
C
      SUBROUTINE LJCOUL(X,V,ENERGY,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), R1,
     1                 R6, ENERGY, DUMMYX, DUMMYY, DUMMYZ, DUMMY,
     2                 COULQ2, COULQ2_ENERGY, R8

C     Divide Coulomb contribution to energy by 4.  This factor is cancelled at the
C     end of the energy loop by the 4 prefactor in the LJ part of the potential.
      COULQ2 = COULQ * COULQ
      COULQ2_ENERGY = COULQ2 / 4.0D0
  
      ENERGY=0.0D0
      DO J1=1,3*NATOMS 
         V(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DUMMYX = X(J3-2)-X(J4-2)
               DUMMYY = X(J3-1)-X(J4-1)
               DUMMYZ = X(J3)-X(J4)
               DIST = DUMMYX**2 + DUMMYY**2 + DUMMYZ**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               ENERGY=ENERGY+DUMMY
               R8=DIST*R6
               DUMMY = -24.0D0*(2.0D0*R6-1.0D0)*R8
               IF ( (J1.LE.COULN).AND.(J2.LE.COULN) ) THEN
                  R1 = SQRT(DIST)
                  ENERGY = ENERGY + COULQ2_ENERGY*R1
                  DUMMY = DUMMY - COULQ2*R1*DIST
               END IF
               V(J3-2) = V(J3-2) + DUMMY * DUMMYX
               V(J3-1) = V(J3-1) + DUMMY * DUMMYY
               V(J3) = V(J3) + DUMMY * DUMMYZ
               V(J4-2) = V(J4-2) - DUMMY * DUMMYX
               V(J4-1) = V(J4-1) - DUMMY * DUMMYY
               V(J4) = V(J4) - DUMMY * DUMMYZ
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               ENERGY=ENERGY+DUMMY
               IF ( (J1.LE.COULN).AND.(J2.LE.COULN) ) THEN
                  ENERGY = ENERGY + COULQ2_ENERGY*SQRT(DIST)
               END IF
            ENDDO
         ENDDO
      ENDIF

      ENERGY=ENERGY*4.0D0

      RETURN
      END

C--------------------------------------------------------------------------------

C     Attempt swap between charged and neutral particles in the LJ+Coulomb
C     potential.

      SUBROUTINE TAKESTEPLJC(NP)

      USE COMMONS, ONLY : NATOMS, COORDS, COULN
      IMPLICIT NONE
      INTEGER NP

      DOUBLE PRECISION DPRAND

      DOUBLE PRECISION X, Y, Z
      INTEGER J1, J2

      IF (COULN==0 .OR. COULN==NATOMS) RETURN

C     Choose one charged particle (1 to COULN) and one neutral (COULN+1 to NATOMS)
      J1 = 3 * (INT(DPRAND()*COULN))
      J2 = 3 * (INT(DPRAND()*(NATOMS-COULN)) + COULN)

      X = COORDS(J1+1, NP)
      Y = COORDS(J1+2, NP)
      Z = COORDS(J1+3, NP)

      COORDS(J1+1, NP) = COORDS(J2+1, NP)
      COORDS(J1+2, NP) = COORDS(J2+2, NP)
      COORDS(J1+3, NP) = COORDS(J2+3, NP)

      COORDS(J2+1, NP) = X
      COORDS(J2+2, NP) = Y
      COORDS(J2+3, NP) = Z

      RETURN
      END
