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
C MNBRAK is called by LINMIN. Given distinct initial points AX and BX,
C it searches in downhill dirction and returns new points AX, BX, CX
C which bracket a minimum of the function.
C function values of the new points: FA, FB, FC
C
      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,F1DIM)
      USE commons
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      EXTERNAL F1DIM
      PARAMETER (GOLD=1.618034D00, GLIMIT=100.D0, TINY=1.D-8)
      LOGICAL FTEST
      COMMON /FAIL/ FTEST

      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=F1DIM(CX)
      IF (FTEST) RETURN
10    IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.0D0*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=F1DIM(U)
          IF (FTEST) RETURN
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 10
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 10
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=F1DIM(U)
          IF (FTEST) RETURN
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=F1DIM(U)
          IF (FTEST) RETURN
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=F1DIM(U)
            IF (FTEST) RETURN
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=F1DIM(U)
          IF (FTEST) RETURN
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=F1DIM(U)
          IF (FTEST) RETURN
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 10
      ENDIF
      RETURN
      END
