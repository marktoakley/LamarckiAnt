!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

SUBROUTINE PERC(P,NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NATOMS, MYUNIT
  LOGICAL, INTENT(IN)  :: DEBUG, RIGID
  LOGICAL, INTENT(OUT) :: PERCT
  DOUBLE PRECISION, INTENT(IN)  :: P(3*NATOMS), PERCCUT
  INTEGER NSITES, J1, J2, NCYCLE, DMIN1, DMAX1, NUNCON1
  INTEGER, ALLOCATABLE :: NDIST1(:)
  DOUBLE PRECISION DUMMY
  LOGICAL CHANGED
  LOGICAL, ALLOCATABLE :: CON(:,:)

  NSITES = NATOMS
  IF (RIGID) THEN
    NSITES = NSITES/2
  ENDIF

  ALLOCATE(NDIST1(NSITES), CON(NSITES,NSITES))

  CON(1:NSITES,1:NSITES)=.FALSE.
  DO J1=1,NSITES
    DO J2=J1+1,NSITES
      DUMMY=(P(3*(J2-1)+1)-P(3*(J1-1)+1))**2+(P(3*(J2-1)+2)-P(3*(J1-1)+2))**2+(P(3*(J2-1)+3)-P(3*(J1-1)+3))**2
      IF (DUMMY.LT.PERCCUT) THEN
        CON(J2,J1)=.TRUE.
        CON(J1,J2)=.TRUE.
!       IF (DEBUG) WRITE(MYUNIT,'(A,2I8)') 'perc> connecting atoms ',J1,J2
      ENDIF
    ENDDO
  ENDDO

! 
! Check that we have a percolating constraint network.
!
  NDIST1(1:NSITES)=1000000
  NDIST1(1)=0
  NCYCLE=0
5 CHANGED=.FALSE.
  NCYCLE=NCYCLE+1
  DMIN1=100000
  DMAX1=0
  NUNCON1=0
  DO J1=1, NSITES
    IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
    DO J2=1, NSITES
      IF (CON(J2,J1)) THEN
        IF (NDIST1(J2)+1.LT.NDIST1(J1)) THEN
          CHANGED=.TRUE.
          NDIST1(J1)=NDIST1(J2)+1
        ENDIF
        IF (NDIST1(J1)+1.LT.NDIST1(J2)) THEN
          CHANGED=.TRUE.
          NDIST1(J2)=NDIST1(J1)+1
        ENDIF
      ENDIF
    ENDDO
    IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
    IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
    IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
  ENDDO
!  PRINT *,'DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED=',DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED
  IF (CHANGED) GOTO 5
  IF (DEBUG) WRITE(MYUNIT,'(3(A,I8))') 'perc> steps to atom 1 converged in ',NCYCLE-1, &
     &                    ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
  PERCT=.TRUE.
  IF (NUNCON1.GT.0) THEN
     PERCT=.FALSE.
     IF (DEBUG) WRITE(MYUNIT,'(3G20.10)') P(1:3*NATOMS)
  ENDIF

END SUBROUTINE PERC
