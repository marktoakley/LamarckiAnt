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
C  Conjugate gradient minimization.
C
      SUBROUTINE CGMIN(MAXI,P,CFLAG,ITER,EREAL,NP)
      USE commons
      IMPLICIT NONE


      INTEGER I, J, MAXI, NP, J1
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.D-6)
      DOUBLE PRECISION G(3*NATOMS),H(3*NATOMS),GRAD(3*NATOMS),P(3*NATOMS),FP,
     1                 PGSUM,PPGSUM,EDIFF,POTEL,QMAX,
     2                 GSUM,EREAL,VAR,FRET,GG,DGG,GAM,QSTART,QFINISH
      INTEGER ITER, MYPSAVE

      LOGICAL STUCK, FTEST, NOTCALLED, CTEST, CFLAG
      COMMON /ST/ STUCK
      COMMON /MYPOT/ POTEL
      COMMON /FAIL/ FTEST
      COMMON /Q4C/ QSTART, QFINISH

      NOTCALLED=.TRUE.

C     CALL ORDERQ4(NATOMS,P,QSTART)

      IF (DEBUG.AND.DUMPT) THEN
         IF (ARNO) THEN
            WRITE(40,'(I4)') NATOMS+2
            WRITE(40,10) NP,NQ(NP)
            WRITE(40,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(40,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(40,65) (P(I),I=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE
            WRITE(40,'(I4)') NATOMS
            WRITE(40,10) NQ(NP)
10          FORMAT(1X,'QUENCH NUMBER ',I6,' initial points in cgmin')
            WRITE(40,'(A2,3F20.10)') ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(40,'(A2,3F20.10)') ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
         ENDIF
      ENDIF

      CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.FALSE.)
      POTEL=EREAL

      IF (FTEST) THEN
         CFLAG=.FALSE.
         RETURN
      ENDIF

      FP=EREAL
      GSUM=RMS
      PGSUM=0.0D0
      PPGSUM=0.0D0
      DO J=1,3*NATOMS
        G(J)=-GRAD(J)
        H(J)=G(J)
        GRAD(J)=H(J)
      ENDDO

      DO ITER=1,MAXI

C        WRITE(40,'(I4)') NATOMS
C        WRITE(40,'(A,I4,A,F15.5)') 'At step number ',ITER,' energy=',EREAL
C        WRITE(40,'(A2,3F20.10)') ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NSEED)

C        PRINT*,'variables:'
C        WRITE(*,'(5F20.10)') (P(J1),J1=1,3*NATOMS)
C        PRINT*,'gradient:'
C        WRITE(*,'(5F20.10)') (GRAD(J1),J1=1,3*NATOMS)

         STUCK=.FALSE.

         IF ((DABS(GSUM-PGSUM)/GSUM.LT.1.0D-5).AND.
     1       (DABS(GSUM-PPGSUM)/GSUM.LT.1.0D-5).AND.(MOD(ITER,3).EQ.0)) THEN
            PRINT*,'STUCK'
            STUCK=.TRUE.
            IF (AMBER.AND.(RMS.GT.1.0D0)) THEN
               EREAL=0.0D0
               POTEL=0.0D0
               RMS=1.0D0
               WRITE(*,'(A)') ' Stuck - step discarded'
               RETURN
C
C  Short range Morse potentials can get stuck through atoms getting too far out
C  of the core. Try contracting.
C
            ELSE IF (MORSET.AND.RHO.GT.6.0D0) THEN
               DO I=1,NATOMS-NSEED
                  P(3*(I-1)+1)=P(3*(I-1)+1)*0.9D0   
                  P(3*(I-1)+2)=P(3*(I-1)+2)*0.9D0   
                  P(3*(I-1)+3)=P(3*(I-1)+3)*0.9D0   
               ENDDO
            ELSE IF (CENT.AND.(.NOT.SEEDT)) THEN
               CALL CENTRE2(P)
            ELSE IF (FIXCOM.AND.(.NOT.SEEDT)) THEN
               CALL CENTRECOM(P)
            ELSE
               CFLAG=.FALSE.
               RETURN
            ENDIF
         ENDIF

         IF (DEBUG) WRITE(6,'(A,G20.10,A,G15.5,A,I4,A)') 
     1                    ' Potential energy=',EREAL,' RMS force=',GSUM,' after ',ITER-1,' CG steps'

         !FIXIMAGE=.TRUE.
         FRET=EREAL
C        CALL MYLINMIN(P,GRAD,NATOMS,FRET)
         CALL LINMIN(ITER,P,GRAD,NATOMS,FRET)
         FIXIMAGE=.FALSE.
         IF (ITER.EQ.1) MYPSAVE=MYPOWER
         IF (FTEST) THEN
            CFLAG=.FALSE.
            RETURN
         ENDIF
         QMAX=1.0D-10
         VAR=QMAX*(ABS(FRET)+ABS(FP)+EPS)
         EDIFF=DABS(FRET-FP)
         IF (NATOMS-NSEED.EQ.1) THEN
            IF (GSUM.LT.GMAX/10.0D0) THEN
               CFLAG=.TRUE.
               RETURN
            ENDIF
         ELSE
            IF (GSUM.LT.GMAX) THEN
               CFLAG=.TRUE.
               RETURN
            ENDIF
         ENDIF
         CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.FALSE.)
         POTEL=EREAL
         IF (EREAL.GT.FP) THEN
            EREAL=0.0D0
            POTEL=0.0D0
            RMS=1.0D0
            WRITE(*,'(A)') ' Energy increased in quench - step discarded'
            RETURN
         ENDIF

C
C  Catch cold fusion for ionic potentials and discard.
C
C  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
C  for systems with positive energies. - khs26 26/11/09
C
         IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO).AND.(EREAL.LT.-1.0D3)) THEN
            EREAL=0.0D6
            POTEL=0.0D6
            RMS=1.0D0
            WRITE(*,'(A)') ' Cold fusion diagnosed - step discarded'
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
            COLDFUSION=.TRUE.
            RETURN
         ENDIF
!        IF ((AMBER.AND.NOTCALLED).AND.(RMS.LT.1.0D0)) THEN
!           CALL CHIRALTEST(CTEST,P)
!           IF (CTEST) THEN
!              WRITE(*,'(A)') ' Change in chirality detected - step rejected'
!              POTEL=1.0D6
!              EREAL=1.0D6
!              RMS=1.0D0
!              RETURN
!           ENDIF
!           NOTCALLED=.FALSE.
!        ENDIF
         FP=EREAL
         GG=0.0D0
         DGG=0.0D0
         DO J=1,3*NATOMS
            GG=GG+G(J)**2
C
C  This is Fletcher-Reeves
C
C           DGG=DGG+GRAD(J)**2
C
C  This is Polak-Ribiere
C
            DGG=DGG+(GRAD(J)+G(J))*GRAD(J)
         ENDDO
         GAM=DGG/GG
         PPGSUM=PGSUM
         PGSUM=GSUM
         GSUM=RMS
         DO J=1,3*NATOMS
            G(J)=-GRAD(J)
            H(J)=G(J)+GAM*H(J)
            GRAD(J)=H(J)
         ENDDO
      ENDDO
 
      CFLAG=.FALSE.

      RETURN
      END
