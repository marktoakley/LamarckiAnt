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
      FUNCTION BRENT(AX,BX,CX,TOLB,XMIN)
      IMPLICIT NONE
      INTEGER ITER,ITMAX
      DOUBLE PRECISION AX,BX,CX,TOLB,XMIN,CGOLD,ZEPS
      PARAMETER (ITMAX=100,CGOLD=0.3819660D0,ZEPS=1.0D-10)
      DOUBLE PRECISION A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XM,F1DIM,BRENT
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.0D0
      FX=F1DIM(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5D0*(A+B)
        TOL1=TOLB*ABS(X)+ZEPS
        TOL2=2.0D0*TOL1
        IF(ABS(X-XM).LE.(TOL2-0.5D0*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) then
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.0d0*(Q-R)
          IF(Q.GT.0.0D0) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=d
          IF(ABS(P).GE.ABS(.5D0*Q*ETEMP).OR.P.LE.Q*(A-X).OR.P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F1DIM(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIf
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      PRINT*, 'brent exceed maximum iterations'
3     XMIN=X
      BRENT=FX
      RETURN
      END
