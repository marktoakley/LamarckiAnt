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
C the function F1DIM has to accompany LINMIN
C
      FUNCTION F1DIM(XX)
      USE commons
      use f1com
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION GRAD(3*NATOMS), XT(3*NATOMS)
      INTEGER J

      DO 10 J=1,NCOM
         XT(J)=PCOM(J)+XX*XICOM(J)
10    CONTINUE
      CALL POTENTIAL(XT,GRAD,POTEL,.FALSE.,.FALSE.)
      F1DIM=POTEL
C     PRINT*,'XX, energy in F1DIM=',XX,POTEL
C     CALL FLUSH(6)
      RETURN
      END
