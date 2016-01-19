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
C Given the function F1DIM, its derivative function DF1DIM
C and a bracketing triplet XA, XB, XC, this routine
C isolates the minimum and returns the abcissa of the minimum as XMIN,
C and the minimum function value is returned as DBRENT.
C
      FUNCTION DBRENT(AX,BX,CX,TOLB,XMIN)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION OLDE
      PARAMETER (ITMAX=100,ZEPS=1.0D-10)
      LOGICAL OK1,OK2
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      XX=V
      E=0.0D0
      CALL ZWISCHEN(XX,DF1DIM,POTEL1)
      FX=POTEL1
      FV=FX
      FW=FX
      DX=DF1DIM
      DV=DX
      DW=DX
      DO 30 ITER=1,ITMAX
        XM=0.50D0*(A+B)
        QMAX1=TOLB*DABS(XX)+ZEPS
        QMAX2=2.0D0*QMAX1
        IF(DABS(XX-XM).LE.(QMAX2-0.50D0*(B-A))) GOTO 40
        IF(DABS(E).GT.QMAX1) THEN
          D1=2.0D0*(B-A)
          D2=D1
          IF(DW.NE.DX) D1=(W-XX)*DX/(DX-DW)
          IF(DV.NE.DX) D2=(V-XX)*DX/(DX-DV)
          U1=XX+D1
          U2=XX+D2
          OK1=((A-U1)*(U1-B).GT.0.).AND.(DX*D1.LE.0.)
          OK2=((A-U2)*(U2-B).GT.0.).AND.(DX*D2.LE.0.)
          OLDE=E
          E=D
          IF(.NOT.(OK1.OR.OK2))THEN
            GO TO 10
          ELSE IF (OK1.AND.OK2)THEN
            IF(DABS(D1).LT.ABS(D2))THEN
              D=D1
            ELSE
              D=D2
            ENDIF
          ELSE IF (OK1)THEN
            D=D1
          ELSE
            D=D2
          ENDIF
          IF(DABS(D).GT.DABS(0.5D0*OLDE))GO TO 10
          U=XX+D
          IF(U-A.LT.QMAX2 .OR. B-U.LT.QMAX2) D=DSIGN(QMAX1,XM-XX)
          GOTO 20
        ENDIF
10      IF(DX.GE.0.) THEN
          E=A-XX
        ELSE
          E=B-XX
        ENDIF
        D=0.5D0*E
20      IF(DABS(D).GE.QMAX1) THEN
          U=XX+D
          CALL ZWISCHEN(U,DF1DIM,POTEL1)
          FU=POTEL1
        ELSE
          U=XX+DSIGN(QMAX1,D)
          CALL ZWISCHEN(U,DF1DIM,POTEL1)
          FU=POTEL1
          IF(FU.GT.FX)GO TO 40
        ENDIF
        DU=DF1DIM
        IF(FU.LE.FX) THEN
          IF(U.GE.XX) THEN
            A=XX
          ELSE
            B=XX
          ENDIF
          V=W
          FV=FW
          DV=DW
          W=XX
          FW=FX
          DW=DX
          XX=U
          FX=FU
          DX=DU
        ELSE
          IF(U.LT.XX) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.XX) THEN
            V=W
            FV=FW
            DV=DW
            W=U
            FW=FU
            DW=DU
          ELSE IF(FU.LE.FV .OR. V.EQ.XX .OR. V.EQ.W) THEN
            V=U
            FV=FU
            DV=DU
          ENDIF
        ENDIF
30    CONTINUE
      PRINT*,'DBRENT exceeded maximum iterations.'
40    XMIN=XX
C     PRINT*,'DBRENT did ',ITER,' iterations.'
      DBRENT=FX
      RETURN
      END
