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
C     This subprogram performs a sort on the input data and
C     arranges it from smallest to biggest. The exchange-sort
C     algorithm is used.
C
      SUBROUTINE SORT4(N,NATOMS,A,F)
      IMPLICIT NONE
      INTEGER J1, L, N, J2, NATOMS, F(NATOMS), NTEMP
      DOUBLE PRECISION A(NATOMS), TEMP
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=F(L)
         F(L)=F(J1)
         F(J1)=NTEMP
20    CONTINUE
      RETURN
      END

      SUBROUTINE SORT4_INT(N,A)
      !same as sort4, but with only one array of integers
      !i.e. simply sort an array of integers of length N
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(INOUT) :: A(N)
      INTEGER J1, L, J2, NTEMP
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
10       CONTINUE
         NTEMP=A(L)
         A(L)=A(J1)
         A(J1)=NTEMP
20    CONTINUE
      RETURN
      END
