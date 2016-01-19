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
C  Energy and gradient for LJ.
C
      SUBROUTINE OLAMI(X,V,EOLAMI,GTEST,SECT)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*MXATMS), DIST, V(3*MXATMS), G(MXATMS,MXATMS), 
     1                 R6, EOLAMI, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 A, B, B2, C, D, V0, V1, V2
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      IF (SECT) THEN
         PRINT*,'ERROR no second derivatives'
         STOP
      ENDIF
      EVAP=.FALSE.
      EOLAMI=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS 
         V(J1)=0.0D0
      ENDDO
      A=1.0D0
      B=1.1D0*A
      C=1.5D0*A
      D=1.618D0*A
      V0=-1.0D0
      V1=0.75D0
      V2=-0.75D0
      B2=B+(B-A)*V1/ABS(V0)
C     WRITE(*,'(A,5F15.5)') 'A,B,B2,C,D=',A,B,B2,C,D
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         IF (DIST.GT.RADIUS) THEN
C           IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
C    1                 'Atom ',J1,' at ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
            EVAP=.TRUE.
            EOLAMI=EOLAMI+10.0D2*(DIST-RADIUS)**2
         ENDIF
         DO J2=J1+1,NATOMS
            J4=3*J2
            DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
C           WRITE(*,'(A,2I5,G20.10)') 'J1,J2,dist=',J1,J2,SQRT(DIST)
            IF (DIST.LT.1.0D0-1.0D-4) THEN
               VT(J1)=VT(J1)+1.0D2
               VT(J2)=VT(J2)+1.0D2
               EOLAMI=EOLAMI+1.0D2
               V(3*(J1-1)+1)=V(3*(J1-1)+1)-(X(3*(J1-1)+1)-X(3*(J2-1)+1))*10.0D0
               V(3*(J1-1)+2)=V(3*(J1-1)+2)-(X(3*(J1-1)+2)-X(3*(J2-1)+2))*10.0D0
               V(3*(J1-1)+3)=V(3*(J1-1)+3)-(X(3*(J1-1)+3)-X(3*(J2-1)+3))*10.0D0
               V(3*(J2-1)+1)=V(3*(J2-1)+1)-(X(3*(J2-1)+1)-X(3*(J1-1)+1))*10.0D0
               V(3*(J2-1)+2)=V(3*(J2-1)+2)-(X(3*(J2-1)+2)-X(3*(J1-1)+2))*10.0D0
               V(3*(J2-1)+3)=V(3*(J2-1)+3)-(X(3*(J2-1)+3)-X(3*(J1-1)+3))*10.0D0
            ELSE IF (DIST.LT.1.0D0) THEN
               VT(J1)=VT(J1)+V0
               VT(J2)=VT(J2)+V0
               EOLAMI=EOLAMI+V0
            ELSE IF (DIST.LT.B2**2) THEN
               VT(J1)=VT(J1)+V0+(V1-V0)*(SQRT(DIST)-A)/(B2-A)
               VT(J2)=VT(J2)+V0+(V1-V0)*(SQRT(DIST)-A)/(B2-A)
               EOLAMI=EOLAMI+V0+(V1-V0)*(SQRT(DIST)-A)/(B2-A)
               V(3*(J1-1)+1)=V(3*(J1-1)+1)+(V1-V0)*(X(3*(J1-1)+1)-X(3*(J2-1)+1))/((B2-A)*SQRT(DIST))
               V(3*(J1-1)+2)=V(3*(J1-1)+2)+(V1-V0)*(X(3*(J1-1)+2)-X(3*(J2-1)+2))/((B2-A)*SQRT(DIST))
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+(V1-V0)*(X(3*(J1-1)+3)-X(3*(J2-1)+3))/((B2-A)*SQRT(DIST))
               V(3*(J2-1)+1)=V(3*(J2-1)+1)+(V1-V0)*(X(3*(J2-1)+1)-X(3*(J1-1)+1))/((B2-A)*SQRT(DIST))
               V(3*(J2-1)+2)=V(3*(J2-1)+2)+(V1-V0)*(X(3*(J2-1)+2)-X(3*(J1-1)+2))/((B2-A)*SQRT(DIST))
               V(3*(J2-1)+3)=V(3*(J2-1)+3)+(V1-V0)*(X(3*(J2-1)+3)-X(3*(J1-1)+3))/((B2-A)*SQRT(DIST))
            ELSE IF (DIST.LT.C**2) THEN
               VT(J1)=VT(J1)+V1
               VT(J2)=VT(J2)+V1
               EOLAMI=EOLAMI+V1
            ELSE IF (DIST.LT.D**2) THEN
               VT(J1)=VT(J1)+V2
               VT(J2)=VT(J2)+V2
               EOLAMI=EOLAMI+V2
            ENDIF
         ENDDO
      ENDDO
      WRITE(*,'(A,G20.10)') 'energy=',EOLAMI
      WRITE(*,'(I5,G20.10)') (J1,V(J1),J1=1,3*NATOMS)
      IF (DEBUG.AND.EVAP) THEN
         WRITE(*,'(A)') 'An atom has evaporated - dumping coordinates'
C        WRITE(40,'(I4)') NATOMS
C        WRITE(40,'(A,I4,A,F15.5)') 'energy after evap=',EOLAMI
C        WRITE(40,'(A2,3F20.10)') ('LJ ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),J1=1,NATOMS)
      ENDIF

      RETURN
      END
