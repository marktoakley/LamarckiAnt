C
C   GMIN: A program for optimizing geometries and calculating reaction pathways
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
      SUBROUTINE MYINERTIA(COORDS,ITDET)
      USE COMMONS, ONLY: NATOMS, ATMASS
      IMPLICIT NONE
      DOUBLE PRECISION ITDET,COORDS(3*NATOMS),DUMQ(3*NATOMS)
      INTEGER J1, J2, J3
      DOUBLE PRECISION IT(3,3), CMX, CMY, CMZ, MASST, IV(3,3)

      DUMQ(1:3*NATOMS)=COORDS(1:3*NATOMS)
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      MASST=0.0D0
      DO J1=1,NATOMS
            CMX=CMX+DUMQ(3*(J1-1)+1)*ATMASS(J1)
            CMY=CMY+DUMQ(3*(J1-1)+2)*ATMASS(J1)
            CMZ=CMZ+DUMQ(3*(J1-1)+3)*ATMASS(J1)
            MASST=MASST+ATMASS(J1)
      ENDDO
      CMX=CMX/MASST
      CMY=CMY/MASST
      CMZ=CMZ/MASST
      DO J1=1,NATOMS
         DUMQ(3*(J1-1)+1)=DUMQ(3*(J1-1)+1)-CMX
         DUMQ(3*(J1-1)+2)=DUMQ(3*(J1-1)+2)-CMY
         DUMQ(3*(J1-1)+3)=DUMQ(3*(J1-1)+3)-CMZ
      ENDDO

      DO J1=1,3
         DO J2=1,3
            IT(J1,J2)=0.0D0
            j3loop1: DO J3=1,NATOMS
               IT(J1,J2)=IT(J1,J2)-DUMQ(3*(J3-1)+J1)*DUMQ(3*(J3-1)+J2)*ATMASS(J3)
            ENDDO j3loop1
            IF (J1.EQ.J2) THEN
               j3loop2: DO J3=1,NATOMS
                  IT(J1,J2)=IT(J1,J2)+(DUMQ(3*(J3-1)+1)**2+DUMQ(3*(J3-1)+2)**2+DUMQ(3*(J3-1)+3)**2)*ATMASS(J3)
               ENDDO j3loop2
            ENDIF
         ENDDO
      ENDDO

C
C Diagonalize inertia tensor. The 0 flag reorders the e/values and
C e/vectors from smallest to largest.
C
      CALL EIG(IT,IV,3,3,0)

      ITDET=IT(1,1)*IT(2,2)*IT(3,3)
      IF (NATOMS.EQ.2) ITDET=IT(2,2)*IT(3,3)

      RETURN
      END

