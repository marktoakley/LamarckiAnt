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
      SUBROUTINE SORT3(N,J3,A,B)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2
      DOUBLE PRECISION TEMP, A(J3), B(3*J3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         DO J2=0,2
            TEMP=B(3*L-J2)
            B(3*L-J2)=B(3*J1-J2)
            B(3*J1-J2)=TEMP
         ENDDO
20    CONTINUE
      RETURN
      END

C ds656> Sort by entries in A(1:N) from biggest to smallest!
      SUBROUTINE SORT3_V2(N,J3,A,X)
      IMPLICIT NONE 
      INTEGER J1, L, N, J2, J3, A(J3), ITEMP
      DOUBLE PRECISION DTEMP, X(J3,3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         ITEMP=A(L)
         A(L)=A(J1)
         A(J1)=ITEMP
         DO J2=1,3
            DTEMP=X(L,J2)
            X(L,J2)=X(J1,J2)
            X(J1,J2)=DTEMP
         ENDDO
20    CONTINUE
      RETURN
      END
