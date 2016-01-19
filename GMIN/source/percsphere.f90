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

SUBROUTINE PERCSPHERE(P,NATOMS,PERCCUT,PERCT)

  IMPLICIT NONE

! NATOMS: the number of bodies
! P: the coordinates of those bodies
! PERCCUT: the square root of the maximum distance for bodies to be considered connected
! PERCT: whether the whole structure is connected
  INTEGER, INTENT(IN) :: NATOMS
  DOUBLE PRECISION, INTENT(INOUT) :: P(3*NATOMS)
  DOUBLE PRECISION, INTENT(IN) :: PERCCUT
  LOGICAL, INTENT(OUT) :: PERCT

  INTEGER NSITES, J1, J2, J3, JCLUST, JR1, JR2
  INTEGER, ALLOCATABLE :: NDIST1(:), JCENTRE(:)
  DOUBLE PRECISION, ALLOCATABLE :: COM(:,:)
  DOUBLE PRECISION DUMMY, DMIN, THETA
  DOUBLE PRECISION ROTVEC(3), RMI(3,3), DRMI(3,3), TEMP(3)
  LOGICAL CHANGED, NEWT, SINGLET
! CON: whether any pair of bodies are connected
  LOGICAL, ALLOCATABLE :: CON(:,:)

  NSITES = NATOMS

  ALLOCATE(NDIST1(NSITES), JCENTRE(NSITES), CON(NSITES,NSITES))
  ALLOCATE(COM(NSITES,3))
  NDIST1(:) = 0
  JCENTRE(:) = 0
  COM(:,:) = 0.0D0

! First, simply scan through all pairs of bodies and call them connected if they lie within PERCCUT squared
  CON(1:NSITES,1:NSITES)=.FALSE.
  DO J1=1,NSITES
    DO J2=J1+1,NSITES
      DUMMY=(P(3*(J2-1)+1)-P(3*(J1-1)+1))**2+(P(3*(J2-1)+2)-P(3*(J1-1)+2))**2+(P(3*(J2-1)+3)-P(3*(J1-1)+3))**2
      IF (DUMMY.LT.PERCCUT**2) THEN
        CON(J2,J1)=.TRUE.
        CON(J1,J2)=.TRUE.
      ENDIF
    ENDDO
  ENDDO

  J1 = 1
  JCLUST = 0

10 JCLUST = JCLUST + 1

! Now we say, if x is connected to y and y to z, set x connected to z
5 CHANGED=.FALSE.
  DO J2=1, NSITES
     IF (CON(J1,J2)) THEN
        DO J3 = 1, NSITES
           IF (CON(J2,J3)) THEN
              IF (CON(J1,J3).EQV..FALSE.) THEN
                 CON(J1,J3) = .TRUE.
                 CON(J3,J1) = .TRUE.
                 CHANGED = .TRUE.
              ENDIF
           ENDIF
        ENDDO
     ENDIF
  ENDDO
! Need to rerun the above if a new connection was made
  IF (CHANGED) GOTO 5

! Now we count the number of isolated clusters
! NDIST1 stores the number of bodies in each cluster
! JCENTRE stores the body used for the start of the search for each cluster
  NDIST1(JCLUST) = 0
  SINGLET = .TRUE.
  DO J2 = 1, NSITES
     IF(CON(J1,J2)) THEN
        NDIST1(JCLUST) = NDIST1(JCLUST) + 1
        SINGLET = .FALSE.
     ENDIF
  ENDDO
  IF (SINGLET) THEN
! This means that there's an isolated body
     NDIST1(JCLUST) = NDIST1(JCLUST) + 1
     CON(J1,J1) = .TRUE.
  ENDIF
  JCENTRE(JCLUST) = J1

! Now check whether there are any bodies not connected to the current cluster
  NEWT = .FALSE.
  DO J3 = 1, NSITES
     IF (CON(J1,J3) .EQV. .FALSE.) THEN
        NEWT = .TRUE.
        DO J2 = 1, JCLUST
! or connected to any previous cluster
           IF (CON(J3,JCENTRE(J2))) THEN
              NEWT = .FALSE.
           ENDIF
        ENDDO
        IF (NEWT) THEN
           J1 = J3
           IF (JCLUST < NSITES) THEN
! If we've reached here, we've found a body not connected to any cluster,
! so use it as the start of a new connection search
              GOTO 10
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  IF (JCLUST .EQ. 1) THEN
     PERCT = .TRUE.
  ELSE
! The goal is now to rotate clusters so they all connect
     PERCT = .FALSE.
     DO J3 = 2, JCLUST

! First, find the pair of bodies from each cluster with the minimum distance
! DMIN stores the minimum distance; JR1 is the index of the body from cluster 1;
! JR2 is the index of the body from cluster J3
        DMIN = 1.0D10
        DO J1 = 1, NSITES
           DO J2 = 1, NSITES
              IF (CON(JCENTRE(1),J1) .AND. CON(JCENTRE(J3),J2) ) THEN
                 DUMMY=(P(3*(J2-1)+1)-P(3*(J1-1)+1))**2+(P(3*(J2-1)+2)-P(3*(J1-1)+2))**2+(P(3*(J2-1)+3)-P(3*(J1-1)+3))**2
                 IF (DUMMY < DMIN) THEN
                    JR1 = J1
                    JR2 = J2
                 ENDIF
              ENDIF
           ENDDO
        ENDDO

! Theta becomes the angle subtended at the origin between JR1 and JR2
        THETA = DOT_PRODUCT(P(3*JR1-2:3*JR1),P(3*JR2-2:3*JR2))
        DUMMY = DOT_PRODUCT(P(3*JR1-2:3*JR1),P(3*JR1-2:3*JR1))
        THETA = THETA / DSQRT(DUMMY)
        DUMMY = DOT_PRODUCT(P(3*JR2-2:3*JR2),P(3*JR2-2:3*JR2))
        THETA = THETA / DSQRT(DUMMY)
        THETA = ACOS(THETA)

! Set up ROTVEC as the angle-axis rotation vector to take JR2 to JR1
        ROTVEC(1) = P(3*JR1-1) * P(3*JR2)   - P(3*JR2-1) * P(3*JR1)
        ROTVEC(2) = P(3*JR1)   * P(3*JR2-2) - P(3*JR2)   * P(3*JR1-2)
        ROTVEC(3) = P(3*JR1-2) * P(3*JR2-1) - P(3*JR2-2) * P(3*JR1-1)
        DUMMY = DOT_PRODUCT(ROTVEC,ROTVEC)
        ROTVEC = THETA * ROTVEC / DSQRT(DUMMY)

! Test rotation of JR2 onto JR1
        CALL RMDRVT(ROTVEC, RMI, DRMI, DRMI, DRMI, .FALSE.)
        TEMP = MATMUL(RMI,P(3*JR2-2:3*JR2))
        DUMMY=(TEMP(3*JR2-2)-P(3*JR1-2))**2+(TEMP(3*JR2-1)-P(3*JR1-1))**2+(TEMP(3*JR2)-P(3*JR1))**2
        IF (DUMMY > PERCCUT) THEN
! Reverse the angle if we've rotated the wrong way
           ROTVEC = -1.0D0 * ROTVEC
        ENDIF

! Reduce the rotation a bit, so we don't put bodies on top of each other
        ROTVEC = ROTVEC * 0.975D0
        CALL RMDRVT(ROTVEC, RMI, DRMI, DRMI, DRMI, .FALSE.)

! Rotate all the bodies connected to the root of cluster J3
        DO J2 = 1, NSITES
           IF ( CON(JCENTRE(J3),J2) ) THEN
              P(3*J2-2:3*J2) = MATMUL(RMI,P(3*J2-2:3*J2))
           ENDIF
        ENDDO
     ENDDO
  ENDIF

  DEALLOCATE(NDIST1, JCENTRE, CON)
  DEALLOCATE(COM)

END SUBROUTINE PERCSPHERE
