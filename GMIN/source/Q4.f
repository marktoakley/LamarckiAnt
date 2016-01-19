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
C******************************************************************************
C
C  Subroutine order calculates Q4 only
C
C******************************************************************************

      SUBROUTINE ORDERQ4(NATOMS,POINTS,CURQ4)

      IMPLICIT DOUBLE PRECISION(A-H,P-Z)
      INTEGER J, K, NB, NATOMS
      DOUBLE PRECISION DX,DY,DZ,DIST,POINTS(3*NATOMS),
     1                 Q40,Q4R(4),Q4I(4), CURQ4

      Q40=0.0d0
      do 10 j=1,4
        Q4R(j)=0.0d0
        Q4I(j)=0.0d0
10    continue

      NB=0
      DO 30 J=1,NATOMS
         DO 40 K=J+1,NATOMS
            DX=POINTS(3*(J-1)+1)-POINTS(3*(K-1)+1)
            DY=POINTS(3*(J-1)+2)-POINTS(3*(K-1)+2)
            DZ=POINTS(3*(J-1)+3)-POINTS(3*(K-1)+3)
            DIST=DSQRT(DX*DX+DY*DY+DZ*DZ)
            IF (DIST.LT.1.3909d0) THEN
               NB=NB+1
               CALL EVASH4(DX/DIST,DY/DIST,DZ/DIST,Q40,Q4R,Q4I)
            ENDIF
40       CONTINUE
30    CONTINUE

      IF (NB.EQ.0) PRINT*,'*** WARNING - NB=0'

      CURQ4 = Q4(NB,Q40,Q4R,Q4I)

      RETURN
      END

C
C**************************************************************
C     
      SUBROUTINE SHINIT
 
C *** Calculate coefficients for Q4, Q6, W4 and W6
C     Notebook JvD A92
 
      COMMON/Q4COEF/ Q4COEF
 
      INTEGER    M, I
      DOUBLE PRECISION       FAC, FACS
 
      DOUBLE PRECISION      Q4COEF(3,0:4)
 
C *** It would be nice indeed to know the exact formula!
 
C     DATA      Q4COEF /   4.375,  -3.75,  0.375,
C    :                   -17.5,     7.5,   0.0,
C    :                   -52.5,    60.0,  -7.5,
C    :                  -105.0,     0.0,   0.0,
C    :                   105.0,     0.0,   0.0 /
 
C *** Now multiply the coefficients with sqrt(2.0* (l-m)!/(l+m)! )

      Q4COEF(1,0)=4.375D0
      Q4COEF(2,0)=-3.75D0
      Q4COEF(3,0)=0.375D0
      Q4COEF(1,1)=-17.5D0
      Q4COEF(2,1)=7.5D0
      Q4COEF(3,1)=0.0D0
      Q4COEF(1,2)=-52.5D0
      Q4COEF(2,2)=60.0D0
      Q4COEF(3,2)=-7.5D0
      Q4COEF(1,3)=-105.0D0
      Q4COEF(2,3)=0.0D0
      Q4COEF(3,3)=0.0D0
      Q4COEF(1,4)=105.0D0
      Q4COEF(2,4)=0.0D0
      Q4COEF(3,4)=0.0D0
 
      DO 99 M=1, 4
         DO 98 I=1, 3
            FACS=DSQRT(2.0*FAC(4-M)/FAC(4+M))
            Q4COEF(I,M)=Q4COEF(I,M)*FACS
98       CONTINUE
99    CONTINUE
 
 
      DO 89 M=1, 6
         DO 88 I=1, 4
            FACS=DSQRT(2.0*FAC(6-M)/FAC(6+M))
C           Q6COEF(I,M)=Q6COEF(I,M)*FACS
88       CONTINUE
89    CONTINUE
 
      RETURN
      END
 
      SUBROUTINE EVASH4(X, Y, Z, Q0, QR, QI)
 
C *** Evaluate the spherical harmonics of degree 4 ********************
C
C     REAL*8   X,Y,Z    A vector on the unit-sphere to be processed
C     DOUBLE PRECISION
C            Q0       Accumulator of P40(cos(THETA))
C            QR(4)    Accumulator of P4m(cos(THETA))*cos(PHI)
C            QI(4)    Accumulator of P4m(cos(THETA))*sin(PHI)
C
C     Note that the actual spherical harmonics Y4m differ from Q0
C     and Q (= QR + i*QI) by a factor
C                     sqrt((9/8*pi)
C     This will be taken care of in FUNCTION Q4
C
C     Angles THETA and PHI are defined in the regular way:
C
C         X = cos(PHI)sin(THETA)
C         Y = sin(PHI)sin(THETA)
C         Z = cos(THETA)
C
C     Notebook JvD A38 and A85
C
C *********************************************************************
 
      COMMON/Q4COEF/ Q4COEF
      DOUBLE PRECISION      Q4COEF(3,0:4)
 
      DOUBLE PRECISION      X, Y, Z
      DOUBLE PRECISION  Q0, QR(4), QI(4)
 
      DOUBLE PRECISION COSTH, COSTH2, COSTH4, SINTH, SINTH2, SCTH
      DOUBLE PRECISION TWOCPH, COSPHI(4), SINPHI(4), TEMP
 
      COSTH  = Z
      COSTH2 = COSTH  * COSTH
      COSTH4 = COSTH2 * COSTH2
      SINTH2 = 1. - COSTH2
      SINTH  = DSQRT(SINTH2)
      SCTH   = SINTH  * COSTH
 
 
C *** Is THETA = 0 ? Then PHI is irrelevant *****************
 
      IF (SINTH .EQ. 0.) THEN
         COSPHI(1) = 1.
         SINPHI(1) = 1.
      ELSE
         COSPHI(1) = X/SINTH
         SINPHI(1) = Y/SINTH
      ENDIF
      TWOCPH    = 2.*COSPHI(1)
      COSPHI(2) = TWOCPH*COSPHI(1)-1.
      SINPHI(2) = TWOCPH*SINPHI(1)
      COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
      SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
      COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
      SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)
 
C *** Now the spherical harmonics are calculated and summed with Q
C
C     This part of the subroutine would have been easier to understand
C     when complex numbers would have been used. However, COMPLEX*16
C     variables are not included in the FORTRAN 77 standard.
 
      Q0    = Q0 + Q4COEF(1,0)*COSTH4+Q4COEF(2,0)*COSTH2
     :             +Q4COEF(3,0)
 
      TEMP  = (Q4COEF(1,1) * COSTH2 + Q4COEF(2,1)) * SCTH
      QR(1) = QR(1) + COSPHI(1)*TEMP
      QI(1) = QI(1) + SINPHI(1)*TEMP
 
      TEMP  = Q4COEF(1,2)*COSTH4+Q4COEF(2,2)*COSTH2+Q4COEF(3,2)
      QR(2) = QR(2) + COSPHI(2)*TEMP
      QI(2) = QI(2) + SINPHI(2)*TEMP
 
      TEMP  = Q4COEF(1,3) * SCTH * SINTH2
      QR(3) = QR(3) + COSPHI(3)*TEMP
      QI(3) = QI(3) + SINPHI(3)*TEMP
 
      TEMP  = Q4COEF(1,4) * SINTH2 * SINTH2
      QR(4) = QR(4) + COSPHI(4)*TEMP
      QI(4) = QI(4) + SINPHI(4)*TEMP
 
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION Q4(NB, Q0, QR, QI)
 
C *** The order parameter Q4 is calculated from Q0 and Q *************
C
C     INTEGER  NB      Number of bonds accumulated
C     DOUBLE PRECISION
C              Q0      Accumulator of P40(cos(THETA))
C              QR(4)   Accumulator of P4m(cos(THETA))*cos(PHI)
C              QI(4)   Accumulator of P4m(cos(THETA))*sin(PHI)
C
C     In the expression for Q4, the factor pi*9/4 disappears
C     The summation has to be performed over the squares of Q4m, with
C     m running from -4 to +4. However, Q4m and Q4-m are equal when
C     squared. So the summation is done for positive m and a factor
C     2 is introduced.
C
C ********************************************************************
 
      INTEGER   NB
      DOUBLE PRECISION  Q0, QR(4), QI(4)
 
      Q4 = DSQRT(    Q0 * Q0
     :          + ( QR(1)*QR(1) + QI(1)*QI(1) )
     :          + ( QR(2)*QR(2) + QI(2)*QI(2) )
     :          + ( QR(3)*QR(3) + QI(3)*QI(3) )
     :          + ( QR(4)*QR(4) + QI(4)*QI(4) ) ) / NB
 
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION FAC(N)
 
      INTEGER I, N
 
      FAC = 1.
      IF(N.LT.2)RETURN
      DO 100 I = 2, N
         FAC = FAC*I
100   CONTINUE
      RETURN
      END
