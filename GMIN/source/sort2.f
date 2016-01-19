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
      SUBROUTINE GSORT2(N,NATOMS)
      USE COMMONS, ONLY : QMINAV, QMINPCSMAV, CSMT, FEBHT
      USE QMODULE
      IMPLICIT NONE
      INTEGER NATOMS
      INTEGER J1, L, N, J3, J2, NTEMP
      DOUBLE PRECISION TEMP, C
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         TEMP=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=TEMP
         NTEMP=QMINNATOMS(L)
         QMINNATOMS(L)=QMINNATOMS(J1)
         QMINNATOMS(J1)=NTEMP
         IF (CSMT) THEN
            TEMP=QMINAV(L)
            QMINAV(L)=QMINAV(J1)
            QMINAV(J1)=TEMP
         ENDIF
C Swap the first found (FF) and number of potential calls when first found (NPCALL_QMIN) array elements to match
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         NTEMP=NPCALL_QMIN(L)
         NPCALL_QMIN(L)=NPCALL_QMIN(J1)
         NPCALL_QMIN(J1)=NTEMP
C Swap the coordinates to match
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
         IF (CSMT) THEN
            DO J2=1,3*NATOMS
               C=QMINPCSMAV(L,J2)
               QMINPCSMAV(L,J2)=QMINPCSMAV(J1,J2)
               QMINPCSMAV(J1,J2)=C
            ENDDO
         ENDIF
20    CONTINUE

      IF (FEBHT) THEN
         DO J1 = N, 1, -1
            QMIN(J1) = QMIN(J1) - QMIN(1)
         ENDDO
      ENDIF

      RETURN
      END
