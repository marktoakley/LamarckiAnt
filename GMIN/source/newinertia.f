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
      SUBROUTINE NEWINERTIA(X,MXATMS,NATOMS,SUM)
      IMPLICIT NONE
      INTEGER NATOMS, J1, J3, MXATMS
      DOUBLE PRECISION X(3*MXATMS), XM, YM, ZM, SUM

      XM=0.0D0
      YM=0.0D0
      ZM=0.0D0
      DO J1=1,NATOMS
         XM=XM+X(3*(J1-1)+1)
         YM=YM+X(3*(J1-1)+2)
         ZM=ZM+X(3*(J1-1)+3)
      ENDDO
      XM=XM/NATOMS
      YM=YM/NATOMS
      ZM=ZM/NATOMS
     
C     DUMMY=0.0D0
C     DO J3=1,NATOMS
C        DUMMY=DUMMY+(X(3*(J3-1)+2)-YM)**2+(X(3*(J3-1)+3)-ZM)**2
C     ENDDO
C     XIN=DUMMY
C     DUMMY=0.0D0
C     DO J3=1,NATOMS
C        DUMMY=DUMMY+(X(3*(J3-1)+1)-XM)**2+(X(3*(J3-1)+3)-ZM)**2
C     ENDDO
C     YIN=DUMMY
C     DUMMY=0.0D0
C     DO J3=1,NATOMS
C        DUMMY=DUMMY+(X(3*(J3-1)+1)-XM)**2+(X(3*(J3-1)+2)-YM)**2
C     ENDDO
C     ZIN=DUMMY
C
C  SUM of eigenvalues is the trace of the matrix.

C     SUM=XIN+YIN+ZIN
      SUM=0.0D0
      DO J3=1,NATOMS
         SUM=SUM+(X(3*(J3-1)+1)-XM)**2+(X(3*(J3-1)+2)-YM)**2+(X(3*(J3-1)+3)-ZM)**2
      ENDDO
      SUM=2.0D0*SUM
C

      RETURN
      END
