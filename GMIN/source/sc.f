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
C      Sutton-Chen potentials
C
      SUBROUTINE SC(X,V,PSC,GRADT)
      USE commons
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, I, J
      LOGICAL GRADT
      DOUBLE PRECISION X(3*NATOMS), DMM(NATOMS,NATOMS), POTA, POTB, SIGMM, DNN(NATOMS,NATOMS),
     1                 V(3*NATOMS), SCRHO(NATOMS), DUMMY, PSC, TPOTA, TPOTB, SIGNN, D2(NATOMS,NATOMS),
     2                 DUMMY2, CPMDFAC
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      EVAP=.FALSE.
      SIGMM=SIG**MM
      SIGNN=SIG**NN
      PSC=0.0D0
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO
      CPMDFAC=1.0D0
C      IF (CPMD) CPMDFAC=0.5291772D0
      DO J1=1,NATOMS
         VT(J1)=0.0D0
         D2(J1,J1)=0.0D0
         DMM(J1,J1)=0.0D0
         DNN(J1,J1)=0.0D0
         J4=3*(J1-1)
         DUMMY=X(J4+1)**2+X(J4+2)**2+X(J4+3)**2
         DUMMY=DUMMY*CPMDFAC**2
C        IF (DUMMY.GT.RADIUS) THEN
C           EVAP=.TRUE.
C           PSC=PSC+(DUMMY-RADIUS)**2
C           IF (GRADT) THEN
C              V(J4+3)=4.0D0*(DUMMY-RADIUS)*X(J4+3)*CPMDFAC
C              V(J4+2)=4.0D0*(DUMMY-RADIUS)*X(J4+2)*CPMDFAC
C              V(J4+1)=4.0D0*(DUMMY-RADIUS)*X(J4+1)*CPMDFAC
C           ENDIF
C        ENDIF
         DO J2=1,J1-1
            J3=3*(J2-1)
            DUMMY=(X(J4+1)-X(J3+1))**2+(X(J4+2)-X(J3+2))**2+(X(J4+3)-X(J3+3))**2
            DUMMY=DUMMY*CPMDFAC**2
            DUMMY=1.0D0/DUMMY
            D2(J2,J1)=DUMMY
            D2(J1,J2)=DUMMY
            DUMMY=DSQRT(DUMMY)
            DMM(J2,J1)=SIGMM*DUMMY**MM
            DMM(J1,J2)=DMM(J2,J1)
            DNN(J2,J1)=SIGNN*DUMMY**NN
            DNN(J1,J2)=DNN(J2,J1)
         ENDDO
      ENDDO
C
C Store density matrix.
C
      DO I=1,NATOMS
         DUMMY=0.0D00
         DO J=1,NATOMS
            DUMMY=DUMMY + DMM(J,I)
         ENDDO
         SCRHO(I)=DSQRT(DUMMY)
      ENDDO
C
C First calculate the potential energy. Choosing SIG=sqrt(2) makes the
C unit of length the bulk nearest-neighbour distance. Choosing SIG to 
C be the bulk lattice constant (=sqrt(2)x nn-distance) gives the original
C SC form. Can also choose sigma=eps=1 and then scale the energies by eps
C and the distances by sig.
C metal  NN  MM   SCEPS/eV    SCC   SIG
C  Cu     9   6   0.012382  39.432
C  Ni     9   6   0.015707  39.432
C  Ag    12   6   0.002542  144.41
C  Au    10   8   0.012793  34.408
C  Pt    10   8   0.019833  34.408
C  Pd    12   7   0.004179  108.27
C
      POTA=0.0D0
      POTB=0.0D0
      DO I=1,NATOMS
         TPOTA=0.0D0
         DO J=1,NATOMS
            TPOTA=TPOTA + DNN(J,I)
         ENDDO
         TPOTA=TPOTA
         POTA=POTA+TPOTA
         TPOTB=SCRHO(I)
         POTB=POTB + TPOTB
         VT(I)=SCEPS*(TPOTA/2.0D0 - TPOTB*SCC)
      ENDDO
      PSC=PSC + SCEPS*(POTA/2.0D0 - POTB*SCC)
C
C     PRINT*,'Sutton-Chen n=', NN,' m=', MM
C     PRINT*,'Two-body contribution=',POTA,' eV'
C     PRINT*,'Many-body contribution=',-POTB,' eV'
C     PRINT*,'Total Energy for last step=', PSC,' eV'
C
      IF (.NOT.GRADT) RETURN
C
C Now calculate the gradient analytically.
C
      DO J1=1,NATOMS
         DUMMY2=1.0D0/SCRHO(J1)
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,NATOMS
               DUMMY=DUMMY+( MM*SCC*DMM(J4,J1)*(DUMMY2 + 1.0D0/SCRHO(J4))/2.0D0
     2                      -NN*DNN(J4,J1) ) * (X(J3)-X(3*(J4-1)+J2))*D2(J4,J1)*CPMDFAC
            ENDDO
            V(J3)=V(J3)+SCEPS*DUMMY
         ENDDO
      ENDDO
C     PRINT*,'energy=',PSC

      RETURN
      END
