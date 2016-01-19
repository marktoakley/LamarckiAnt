C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      PROGRAM MAKEXYZ
      LOGICAL YESNO

C     READ(*,*)
      INQUIRE(FILE='vector.dump',EXIST=YESNO)
      IF (YESNO) OPEN(UNIT=8,FILE='vector.dump',STATUS='OLD')
C     BTOA=0.529177D0
      BTOA=1.0D0
      DO J1=1,16000
         REWIND(8)
         READ(*,*,END=20) X, Y, Z
         IF (YESNO) READ(8,*) XV
         IF (YESNO) READ(8,*) XV, YV, ZV
         PRINT*,'38'
         PRINT*
         DO J2=1,38
            WRITE(*,10) 'LA ',X*BTOA,Y*BTOA,Z*BTOA,XV, YV, ZV
10          FORMAT(A2,6F20.10)
            IF (J2.NE.38) THEN
               READ(*,*,END=20) X, Y, Z
               IF (YESNO) READ(8,*) XV, YV, ZV
            ENDIF
         ENDDO
      ENDDO

20    STOP
      END
