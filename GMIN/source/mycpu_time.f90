!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MYCPU_TIME(POO)
USE COMMONS
IMPLICIT NONE
DOUBLE PRECISION MYTIME, POO

CALL CPU_TIME(MYTIME)
POO=MYTIME ! without this extra assignment NAG f95 returns a random number!
           ! saving TSTART is necessary for PG compiler, where initial time is
           ! not zero!
! IF (MYTIME-TSTART.GT.TIMELIMIT) THEN
!    PRINT '(A,F12.1,A,F12.1)','CPU time limit exceeded: time=',MYTIME-TSTART,' limit=',TIMELIMIT
!    IF (DUMPSP) CALL DODUMPSP
!    STOP 1
! ENDIF

RETURN
END SUBROUTINE MYCPU_TIME
