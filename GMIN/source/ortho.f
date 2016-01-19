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
      SUBROUTINE ORTHO(VEC1,COORDS,NOPT,NDEC,NATOMS)
      IMPLICIT NONE
      INTEGER J2,NOPT,NDEC,NATOMS,J3
      DOUBLE PRECISION DUMMY1, DUMMY2, VEC1(NDEC), COORDS(NDEC)
C
C  X trans
C
      DUMMY1=0.0D0
      DO J2=1,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/(1.0D0*NATOMS)
      DO J2=1,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
C
C  Y trans
C
      DUMMY1=0.0D0
      DO J2=2,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/(1.0D0*NATOMS)
      DO J2=2,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
C
C  Z trans
C
      DUMMY1=0.0D0
      DO J2=3,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/(1.0D0*NATOMS)
      DO J2=3,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
C
C  For a periodic system we should now return
C

C
C  X rot
C
      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-1)*COORDS(J3)-VEC1(J3)*COORDS(J3-1)
         DUMMY2=DUMMY2+COORDS(J3)**2+COORDS(J3-1)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-1)=VEC1(J3-1)-DUMMY2*COORDS(J3)
            VEC1(J3)=VEC1(J3)+DUMMY2*COORDS(J3-1)
         ENDDO
      ENDIF
C
C  Y rot
C
      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1-VEC1(J3-2)*COORDS(J3)+VEC1(J3)*COORDS(J3-2)
         DUMMY2=DUMMY2+COORDS(J3)**2+COORDS(J3-2)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)+DUMMY2*COORDS(J3)
            VEC1(J3)=VEC1(J3)-DUMMY2*COORDS(J3-2)
         ENDDO
      ENDIF
C
C  Z rot
C
      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-2)*COORDS(J3-1)-VEC1(J3-1)*COORDS(J3-2)
         DUMMY2=DUMMY2+COORDS(J3)**2+COORDS(J3-2)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)-DUMMY2*COORDS(J3-1)
            VEC1(J3-1)=VEC1(J3-1)+DUMMY2*COORDS(J3-2)
         ENDDO
      ENDIF

      RETURN
      END
