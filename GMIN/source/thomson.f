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
C  Energy and gradient for Thomson problem in theta, phi coordinates.
C
      SUBROUTINE THOMSON(X,V,ETHOMSON,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(*), DIST, V(*), ETHOMSON, DUMMY, CT1, ST1, CT2, ST2, CPDIFF, SPDIFF
      DOUBLE PRECISION COST(NATOMS), SINT(NATOMS), COSP(NATOMS), SINP(NATOMS)
      DOUBLE PRECISION, PARAMETER :: SR2=1.4142135623730950488D0

      DO J1=1,NATOMS
         J3=2*J1
         COST(J1)=COS(X(J3-1))
         SINT(J1)=SIN(X(J3-1))
         COSP(J1)=COS(X(J3))
         SINP(J1)=SIN(X(J3))
      ENDDO

      ETHOMSON=0.0D0
      VT(1:NATOMS)=0.0D0
      V(1:2*NATOMS)=0.0D0
      DO J1=1,NATOMS
         J3=2*J1
         CT1=COST(J1)
         ST1=SINT(J1)
         DO J2=J1+1,NATOMS
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
C           DIST=1.0D0/(SR2*SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2))
            DIST=1.0D0/SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            ETHOMSON=ETHOMSON+DIST
            DIST=DIST**3
            DUMMY=SPDIFF*ST1*ST2*DIST
            V(J3-1)=V(J3-1)+(CPDIFF*CT1*ST2-CT2*ST1)*DIST
            V(J3)=V(J3)    -DUMMY
            V(J4-1)=V(J4-1)+(CPDIFF*CT2*ST1-CT1*ST2)*DIST
            V(J4)=V(J4)    +DUMMY
         ENDDO
      ENDDO
      ETHOMSON=ETHOMSON/SR2
      DO J1=1,2*NATOMS
         V(J1)=V(J1)/SR2
      ENDDO

      RETURN
      END
C
C  Energy and gradient for Thomson problem in theta, phi coordinates.
C  In this case there is one odd charge with magnitude different from unity.
C
      SUBROUTINE ODDTHOMSON(X,V,ETHOMSON,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(*), DIST, V(*), ETHOMSON, DUMMY, CT1, ST1, CT2, ST2, CPDIFF, SPDIFF
      DOUBLE PRECISION COST(NATOMS), SINT(NATOMS), COSP(NATOMS), SINP(NATOMS)
      DOUBLE PRECISION, PARAMETER :: SR2=1.4142135623730950488D0

      DO J1=1,NATOMS
         J3=2*J1
         COST(J1)=COS(X(J3-1))
         SINT(J1)=SIN(X(J3-1))
         COSP(J1)=COS(X(J3))
         SINP(J1)=SIN(X(J3))
      ENDDO

      ETHOMSON=0.0D0
      VT(1:NATOMS)=0.0D0
      V(1:2*NATOMS)=0.0D0
C
C  Terms involving the odd charge.
C
      CT1=COST(1)
      ST1=SINT(1)
      DO J2=2,NATOMS
         J4=2*J2
         CT2=COST(J2)
         ST2=SINT(J2)
         CPDIFF=COSP(1)*COSP(J2)+SINP(1)*SINP(J2)
         SPDIFF=COSP(J2)*SINP(1)-SINP(J2)*COSP(1)
         DIST=1.0D0/SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
         ETHOMSON=ETHOMSON+DIST
         DIST=DIST**3
         DIST=DIST*ODDCHARGE ! put the odd charge into the gradient; energy corrected below
         DUMMY=SPDIFF*ST1*ST2*DIST
         V(1)=V(1)+(CPDIFF*CT1*ST2-CT2*ST1)*DIST
         V(2)=V(2)    -DUMMY
         V(J4-1)=V(J4-1)+(CPDIFF*CT2*ST1-CT1*ST2)*DIST
         V(J4)=V(J4)    +DUMMY
      ENDDO
      ETHOMSON=ETHOMSON*ODDCHARGE

      DO J1=2,NATOMS
         J3=2*J1
         CT1=COST(J1)
         ST1=SINT(J1)
         DO J2=J1+1,NATOMS
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
            DIST=1.0D0/SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            ETHOMSON=ETHOMSON+DIST
            DIST=DIST**3
            DUMMY=SPDIFF*ST1*ST2*DIST
            V(J3-1)=V(J3-1)+(CPDIFF*CT1*ST2-CT2*ST1)*DIST
            V(J3)=V(J3)    -DUMMY
            V(J4-1)=V(J4-1)+(CPDIFF*CT2*ST1-CT1*ST2)*DIST
            V(J4)=V(J4)    +DUMMY
         ENDDO
      ENDDO
      ETHOMSON=ETHOMSON/SR2
      DO J1=1,2*NATOMS
         V(J1)=V(J1)/SR2
      ENDDO

      RETURN
      END
!
!  Subroutine to convert Cartesians to theta, phi.
!
      SUBROUTINE THOMSONCTOANG(COORDS,P,NATOMS,MYUNIT)
      IMPLICIT NONE
      INTEGER NATOMS, J1, MYUNIT
      DOUBLE PRECISION COORDS(3*NATOMS), P(3*NATOMS), DIST
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0

      DO J1=1,NATOMS
         DIST=SQRT(COORDS(3*(J1-1)+1)**2+COORDS(3*(J1-1)+2)**2+COORDS(3*(J1-1)+3)**2)
         COORDS(3*(J1-1)+1)=COORDS(3*(J1-1)+1)/DIST
         COORDS(3*(J1-1)+2)=COORDS(3*(J1-1)+2)/DIST
         COORDS(3*(J1-1)+3)=COORDS(3*(J1-1)+3)/DIST
         P(2*(J1-1)+1)=ACOS(COORDS(3*(J1-1)+3))
         IF (ABS(COORDS(3*(J1-1)+3)-COS(P(2*(J1-1)+1))).GT.1.0D-10) THEN
             WRITE(MYUNIT, '(A)') 'inconsistent conversion for z'
            STOP
         ENDIF
         IF (COORDS(3*(J1-1)+1).EQ.0.0D0) THEN
            IF (COORDS(3*(J1-1)+2).GT.0.0D0) THEN
               P(2*(J1-1)+2)=HALFPI
            ELSE 
               P(2*(J1-1)+2)=-HALFPI
            ENDIF
         ELSE IF (COORDS(3*(J1-1)+2).EQ.0.0D0) THEN
            IF (COORDS(3*(J1-1)+1).GT.0.0D0) THEN
               P(2*(J1-1)+2)=0.0D0
            ELSE 
               P(2*(J1-1)+2)=2*HALFPI
            ENDIF
         ELSE
            P(2*(J1-1)+2)=ATAN(COORDS(3*(J1-1)+2)/COORDS(3*(J1-1)+1))
         ENDIF
         IF (ABS(COORDS(3*(J1-1)+1)-SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2))).GT.1.0D-5) THEN
            P(2*(J1-1)+2)=P(2*(J1-1)+2)+2*HALFPI
            IF (ABS(COORDS(3*(J1-1)+1)-SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                WRITE(MYUNIT, '(A)') 'inconsistent conversion for x'
               STOP
            ENDIF
         ENDIF
         IF (ABS(COORDS(3*(J1-1)+2)-SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2))).GT.1.0D-5) THEN
            P(2*(J1-1)+2)=-P(2*(J1-1)+2)
            IF (ABS(COORDS(3*(J1-1)+2)-SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                WRITE(MYUNIT, '(A)') 'inconsistent conversion for y'
                WRITE(MYUNIT, '(A,3G20.10)') 'x,y,z:      ',COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
                WRITE(MYUNIT, '(A,3G20.10)') 'theta,phi: ',P(2*(J1-1)+1),P(2*(J1-1)+2)
                WRITE(MYUNIT, '(A,3G20.10)') 'x,y,z calc: ',SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2)),
     &                                            SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2)),
     &                                            COS(P(2*(J1-1)+1))
               STOP
            ENDIF
         ENDIF
      ENDDO
      RETURN

      END
!
!  Subroutine to convert theta, phi to Cartesians.
!
      SUBROUTINE THOMSONANGTOC(P,NATOMS)
      IMPLICIT NONE
      INTEGER NATOMS, J1
      DOUBLE PRECISION TMPCOORDS(3*NATOMS), P(3*NATOMS)

      TMPCOORDS(1:2*NATOMS)=P(1:2*NATOMS) 
      DO J1=1,NATOMS
         P(3*(J1-1)+1)=SIN(TMPCOORDS(2*(J1-1)+1))*COS(TMPCOORDS(2*(J1-1)+2))
         P(3*(J1-1)+2)=SIN(TMPCOORDS(2*(J1-1)+1))*SIN(TMPCOORDS(2*(J1-1)+2))
         P(3*(J1-1)+3)=COS(TMPCOORDS(2*(J1-1)+1))
      ENDDO

      RETURN
      END

