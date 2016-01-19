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
C*************************************************************************
C
C  Subroutine MORSE calculates the cartesian gradient analytically for 
C  the Morse potential.
C
C*************************************************************************
C
      SUBROUTINE MORSE(X,V,EMORSE,GTEST)
      USE commons
      IMPLICIT NONE 
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), DIST, R, DUMMY,
     1                 RR(NATOMS,NATOMS), EMORSE, DUMMYX, 
     2                 DUMMYY, DUMMYZ, XMUL2
!     LOGICAL EVAP, evapreject
!     COMMON /EV/ EVAP, evapreject

!     EVAP=.FALSE.
      EMORSE=0.0D0
      DO J1=1,NATOMS
         VT(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            RR(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=MAX(DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
     1                  + (X(J3)-X(J4))**2),1.0D-5)
!              DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
!    1                  + (X(J3)-X(J4))**2)
               R=DEXP(RHO-RHO*DIST)
               RR(J2,J1)=2.0D0*R*(R-1.0D0)/DIST
               RR(J1,J2)=RR(J2,J1)
               DUMMY=R*(R-2.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
     1                  + (X(J3)-X(J4))**2)
               R=DEXP(RHO-RHO*DIST)
               DUMMY=R*(R-2.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) RETURN
 
      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
            V(J3-2)=0.0D0
            V(J3-1)=0.0D0
            V(J3)=0.0D0
         ELSE
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=RR(J4,J1)
               DUMMYX=DUMMYX-(X(J3-2)-X(J2-2))*XMUL2
               DUMMYY=DUMMYY-(X(J3-1)-X(J2-1))*XMUL2
               DUMMYZ=DUMMYZ-(X(J3)-X(J2))*XMUL2
            ENDDO
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
            IF (DIST.GT.RADIUS) THEN
               DUMMYX=DUMMYX+4.0D0*(DIST-RADIUS)*X(J3-2)
               DUMMYY=DUMMYY+4.0D0*(DIST-RADIUS)*X(J3-1)
               DUMMYZ=DUMMYZ+4.0D0*(DIST-RADIUS)*X(J3)
            ENDIF
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

      RETURN
      END
