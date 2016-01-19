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
C  Add energy and gradient correction terms for compression
C
      SUBROUTINE COMPRESS(X,V,ENERGY,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3, NMOL
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), ENERGY, XMASS, YMASS, ZMASS

C jwrm2> for rigid body angle axis system, only the first half of X is the body position
C        Note, as written this is a harmonic potential acting on the centre of the rigid
C        body, NOT on each site of the rigid body 
      NMOL = NATOMS
      IF (RIGID) NMOL = NATOMS/2

C Find centre of mass
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NMOL
         XMASS=XMASS+X(3*(J1-1)+1)
         YMASS=YMASS+X(3*(J1-1)+2)
         ZMASS=ZMASS+X(3*(J1-1)+3)
      ENDDO
      XMASS=XMASS/NMOL
      YMASS=YMASS/NMOL
      ZMASS=ZMASS/NMOL

      IF (PERIODIC) RETURN

      IF (DEBUG) WRITE (MYUNIT,*) 'compress> Compressing'

C Loop over all bodies
      DO J1=1, NMOL
         J3=3*J1
         DIST=(X(J3-2)-XMASS)**2+(X(J3-1)-YMASS)**2+(X(J3)-ZMASS)**2

C Add energy
         ENERGY=ENERGY+K_COMP*DIST/2.0D0

         IF (GTEST) THEN
C Add gradients
            V(J3-2)=V(J3-2)+K_COMP*(X(J3-2)-XMASS)
            V(J3-1)=V(J3-1)+K_COMP*(X(J3-1)-YMASS)
            V(J3)=  V(J3)  +K_COMP*(X(J3)  -ZMASS)
         ENDIF
      ENDDO

      RETURN
      END
