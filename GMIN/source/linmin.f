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
C The subroutine LINMIN takes the 3*N dimensional point P and the 3*N dimensio-
C nal direction GRAD, moves and resets P to where the function FUNC(P) takes on
C a minimum in the direction GRAD from P, and replaces GRAD by the actual vector
C displacement that P was moved.
C The value of FUNC at the returned location P is returned as FRET.
C
C
      SUBROUTINE LINMIN(ITER,P,GRAD,N,FRET)
      USE commons
      use f1com
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      INTEGER ITER, N, J, NCOUNT
      DOUBLE PRECISION F1DIM
      LOGICAL STUCK, FTEST
      COMMON /ST/ STUCK
      COMMON /FAIL/ FTEST
      EXTERNAL F1DIM
      DIMENSION GRAD(3*NATOMS),P(3*NATOMS)
      SAVE BX

      NCOM=3*N
      DO J=1,3*N
        PCOM(J)=P(J)
        XICOM(J)=GRAD(J)
        DUMMY=DUMMY+GRAD(J)**2
      ENDDO

      IF (STUCK) BX=1.0D0

      NCOUNT=0
      AX=0.0D0
C
C  We already know FA!
C
C     FA=F1DIM(AX)
C
      FA=FRET

      IF (FTEST) RETURN

      IF (BX.EQ.0.0D0) BX=1.0D0
C     BX=MIN(DABS(BX),STEP(1)/(DSQRT(3.0D0*NATOMS)*RMS))
      BX=MIN(DABS(BX),0.1D0/(DSQRT(3.0D0*NATOMS)*RMS))

C     IF ((BX.EQ.0.0D0).OR.(ITER.EQ.1)) BX=1.0D0/(SQRT(3.0D0*NATOMS)*RMS)

15    FB=F1DIM(BX)
C     PRINT*,'AX,BX,FA,FB=',AX,BX,FA,FB
      IF (FTEST) RETURN
      IF ((FB.GT.FA+ECONV).AND.(NCOUNT.LT.10)) THEN
         CX=BX
         FC=FB
         BX=BX/10.0D0
         NCOUNT=NCOUNT+1
         GOTO 15
      ENDIF

      EMIN1=FB-FA
C     WRITE(*,'(A,E20.10)') 'Initial energy lowering in LINMIN=',EMIN1
      DUMMY=FB

C
C  If NCOUNT > 0 we have already bracketed the minimum!
C
      IF ((NCOUNT.EQ.0).OR.(EMIN1.GE.0.0D0)) THEN
         CALL MNBRAK(AX,BX,CX,FA,FB,FC,F1DIM)
C        WRITE(*,'(A,6F15.7)') 'after MNBRACK AX,BX,CX,FA,FB,FC=',AX,BX,CX,FA,FB,FC
         EMIN2=FB-DUMMY
C        WRITE(*,'(A,E20.10)')    'Energy lowering by MNBRAK=        ',EMIN2
      ENDIF
      DUMMY=FB

      IF (DBRENTT) THEN
         FRET=DBRENT(AX,BX,CX,TOLB,XMIN)
         EMIN3=FRET-DUMMY
C        WRITE(*,'(A,E20.10)') 'Energy lowering by DBRENT=        ',EMIN3
      ELSE
         FRET=BRENT(AX,BX,CX,TOLB,XMIN)
         EMIN3=FRET-DUMMY
C        WRITE(*,'(A,E20.10)') 'Energy lowering by BRENT=         ',EMIN3
      ENDIF

      DO 20 J=1,3*N
        GRAD(J)=XMIN*GRAD(J)
        P(J)=P(J)+GRAD(J)
20    CONTINUE

      BX=BX*10.0D0

      RETURN
      END
