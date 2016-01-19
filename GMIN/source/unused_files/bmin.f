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
      SUBROUTINE BMIN(MAXI,P,CFLAG,ITER,EREAL,NP)
      USE COMMONS
      USE MODAMBER
      IMPLICIT NONE

      INTEGER NMAX, MAXI, NP, IFLAG, J2, J3, INFO, J1, NOPT
      DOUBLE PRECISION EPS
      PARAMETER (NMAX=3*MXATMS,EPS=1.D-6)
      DOUBLE PRECISION GRAD(NMAX),P(NMAX),POTEL,EREAL,QSTART,QFINISH,XTOL,BVEC(3*MXATMS)
      DOUBLE PRECISION PP(3*MXATMS), GP(3*MXATMS), DUMMY, R12, R13, TH, EVEC(3), FVEC(3)
      INTEGER ITER
      LOGICAL ITINC, BFGSTST
      DOUBLE PRECISION DIAG(3*MXATMS),BASE,DIGITS,WORK(3*MXATMS*(2*MUPDATE+1)+2*MUPDATE)
      EXTERNAL LB2
      LOGICAL DZT
      COMMON /DZ/ DZT

      LOGICAL FTEST, NOTCALLED, CTEST, CFLAG
      COMMON /MYPOT/ POTEL
      COMMON /FAIL/ FTEST
      COMMON /Q4C/ QSTART, QFINISH

      NOTCALLED=.TRUE.
      BFGSTST=.FALSE.

C     BASE=DLAMCH('B')
C     DIGITS=DLAMCH('N')
C     XTOL=BASE**(1 - DIGITS)
      BASE=2
      DIGITS=53
      XTOL=2.2204460492503D-16
      IF (DEBUG) IPRINT(1)=1
     
      IF (TNM6) THEN
         DUMMY=SQRT((P(4)-P(1))**2+(P(5)-P(2))**2+(P(6)-P(3))**2)
         EVEC(1)=(P(4)-P(1))/DUMMY
         EVEC(2)=(P(5)-P(2))/DUMMY
         EVEC(3)=(P(6)-P(3))/DUMMY
         DUMMY=SQRT((-(P(3)**2*P(5)) + P(3)*P(5)*P(6) - P(4)*P(5)*P(7) + P(3)**2*P(8) + P(4)**2*P(8) - 
     -      2*P(3)*P(6)*P(8) + P(6)**2*P(8) + P(1)**2*(-P(5) + P(8)) + 
     -      P(1)*(P(2)*(P(4) - P(7)) + P(5)*(P(4) + P(7)) - 2*P(4)*P(8)) - 
     -      P(2)*(P(4)**2 - P(4)*P(7) - (P(3) - P(6))*(P(6) - P(9))) + P(5)*(P(3) - P(6))*P(9))**2 + 
     -   (-(P(1)*P(5)**2) - P(1)*P(6)**2 + P(5)**2*P(7) + P(6)**2*P(7) + P(2)**2*(-P(4) + P(7)) + 
     -      P(3)**2*(-P(4) + P(7)) + P(1)*P(5)*P(8) - P(4)*P(5)*P(8) + 
     -      P(2)*(P(5)*(P(1) + P(4) - 2*P(7)) - (P(1) - P(4))*P(8)) + (P(1) - P(4))*P(6)*P(9) + 
     -      P(3)*(P(6)*(P(1) + P(4) - 2*P(7)) - (P(1) - P(4))*P(9)))**2 + 
     -   (-(P(6)*(P(4)*P(7) + (P(2) - P(5))*(P(2) - P(8)))) - 
     -      P(3)*(P(4)**2 - P(4)*P(7) - (P(2) - P(5))*(P(5) - P(8))) + (P(4)**2 + (P(2) - P(5))**2)*P(9) + 
     -      P(1)**2*(-P(6) + P(9)) + P(1)*(P(3)*(P(4) - P(7)) + P(6)*(P(4) + P(7)) - 2*P(4)*P(9)))**2)
         FVEC(1)=(-(P(1)*P(5)**2) - P(1)*P(6)**2 + P(5)**2*P(7) + P(6)**2*P(7) + P(2)**2*(-P(4) + P(7)) + 
     -              P(3)**2*(-P(4) + P(7)) + P(1)*P(5)*P(8) - P(4)*P(5)*P(8) + 
     -              P(2)*(P(5)*(P(1) + P(4) - 2*P(7)) - (P(1) - P(4))*P(8)) + (P(1) - P(4))*P(6)*P(9) + 
     -              P(3)*(P(6)*(P(1) + P(4) - 2*P(7)) - (P(1) - P(4))*P(9)))/DUMMY
         FVEC(2)=(-(P(3)**2*P(5)) + P(3)*P(5)*P(6) - P(4)*P(5)*P(7) + P(3)**2*P(8) + P(4)**2*P(8) - 
     -            2*P(3)*P(6)*P(8) + P(6)**2*P(8) + P(1)**2*(-P(5) + P(8)) + 
     -              P(1)*(P(2)*(P(4) - P(7)) + P(5)*(P(4) + P(7)) - 2*P(4)*P(8)) - 
     -              P(2)*(P(4)**2 - P(4)*P(7) - (P(3) - P(6))*(P(6) - P(9))) + P(5)*(P(3) - P(6))*P(9))/DUMMY
         FVEC(3)=(-(P(6)*(P(4)*P(7) + (P(2) - P(5))*(P(2) - P(8)))) - 
     -              P(3)*(P(4)**2 - P(4)*P(7) - (P(2) - P(5))*(P(5) - P(8))) + (P(4)**2 + (P(2) - P(5))**2)*P(9) + 
     -              P(1)**2*(-P(6) + P(9)) + P(1)*(P(3)*(P(4) - P(7)) + P(6)*(P(4) + P(7)) - 2*P(4)*P(9)))/DUMMY
      ENDIF
C
C  Subsequent minimisations seem to work better without reinitialisations,
C  despite the IFLAG<0 messages. Not true for AMBER, though!
C
C     IFLAG=0

      IF (AMBER) THEN
         IFLAG=0
         INFO=0
      ENDIF

C     CALL ORDERQ4(NATOMS,P,QSTART)

      IF (DEBUG.AND.DUMPT) THEN
         IF (ARNO) THEN
            WRITE(40,'(I4)') NATOMS+2
            WRITE(40,10) NP,NQ(NP)
            WRITE(40,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(40,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(40,65) (P(J1),J1=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (AMBER) THEN
            WRITE(40,'(I4)') NATOMS
            WRITE(40,10) NP,NQ(NP)
            DO J2=1,NATOMS
               WRITE(40,'(A,3F20.10)') typech(J2)(1:1),(P(3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            WRITE(40,'(I4)') NATOMS
            WRITE(40,10) NQ(NP)
10          FORMAT(1X,'QUENCH NUMBER ',I6,' initial points in bmin')
            WRITE(40,'(A2,3F20.10)') ('LA ',P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(40,'(A2,3F20.10)') ('LB',P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
         ENDIF
      ENDIF

      DO ITER=1,MAXI
11       CONTINUE !FIXIMAGE=.TRUE.
         ITINC=.FALSE.
C        IF (DEBUG.AND.DUMPT) THEN
C           WRITE(40,'(I4)') NATOMS
C           WRITE(40,'(A,I4,A,F15.5,A,F15.5)') 'At step number ',ITER,' energy=',EREAL,' RMS=',RMS
C           WRITE(40,'(A2,3F20.10)') ('LA ',P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3),J1=1,NATOMS-NS)
C           IF (NS.GT.0) WRITE(40,'(A2,3F20.10)') ('LB',P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
C        ENDIF
         CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.FALSE.)
         IF (TNM6) THEN
            NOPT=3*NATOMS-6
            R12=SQRT((P(4)-P(1))**2+(P(5)-P(2))**2+(P(6)-P(3))**2)
            R13=SQRT((P(7)-P(1))**2+(P(8)-P(2))**2+(P(9)-P(3))**2)
            TH=ACOS(( (P(4)-P(1))*(P(7)-P(1))+(P(5)-P(2))*(P(8)-P(2))+(P(6)-P(3))*(P(9)-P(3)) )/(R12*R13))
            PP(1)=R12
            PP(2)=R13
            PP(3)=TH
            GP(1)=GRAD(4)*EVEC(1)+GRAD(5)*EVEC(2)+GRAD(6)*EVEC(3)
            GP(2)=GRAD(7)*(EVEC(1)*COS(TH)+FVEC(1)*SIN(TH))
     1           +GRAD(8)*(EVEC(2)*COS(TH)+FVEC(2)*SIN(TH))
     2           +GRAD(9)*(EVEC(3)*COS(TH)+FVEC(3)*SIN(TH))
            GP(3)=GRAD(7)*R13*(-EVEC(1)*SIN(TH)+FVEC(1)*COS(TH))
     1           +GRAD(8)*R13*(-EVEC(2)*SIN(TH)+FVEC(2)*COS(TH))
     2           +GRAD(9)*R13*(-EVEC(3)*SIN(TH)+FVEC(3)*COS(TH))
            RMS=GP(1)**2+GP(2)**2+GP(3)**2
            DO J1=10,3*NATOMS
               PP(J1-6)=P(J1)
               GP(J1-6)=GRAD(J1)
               RMS=RMS+GP(J1)**2
            ENDDO
            RMS=SQRT(RMS/NOPT)
         ELSE
            NOPT=3*NATOMS
            DO J1=1,3*NATOMS
               PP(J1)=P(J1)
               GP(J1)=GRAD(J1)
            ENDDO
         ENDIF
         
         POTEL=EREAL
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
         IF ((AMBER.AND.NOTCALLED).AND.(RMS.LT.1.0D0)) THEN
            CALL CHIRALTEST(CTEST,P)
            IF (CTEST) THEN
               WRITE(*,'(A)') ' Change in chirality detected - step rejected'
               POTEL=1.0D6
               EREAL=1.0D6
               RMS=1.0D0
               RETURN
            ENDIF
            NOTCALLED=.FALSE.
         ENDIF

         CALL LBFGS(NOPT,MUPDATE,PP,EREAL,GP,.FALSE.,DIAG,IPRINT,GMAX,XTOL,WORK,IFLAG,
     1              DGUESS,ITINC,BFGSTST,BVEC,FTOL,INFO,AMBER)
C        PRINT*,'After LBFGS INFO, IFLAG=',INFO, IFLAG
         IF (TNM6) THEN
            P(4)=P(1)+PP(1)*EVEC(1)
            P(5)=P(2)+PP(1)*EVEC(2)
            P(6)=P(3)+PP(1)*EVEC(3)
            P(7)=P(1)+PP(2)*(SIN(PP(3))*FVEC(1)+COS(PP(3))*EVEC(1))
            P(8)=P(2)+PP(2)*(SIN(PP(3))*FVEC(2)+COS(PP(3))*EVEC(2))
            P(9)=P(3)+PP(2)*(SIN(PP(3))*FVEC(3)+COS(PP(3))*EVEC(3))
            DO J1=10,3*NATOMS
               P(J1)=PP(J1-6)
            ENDDO
         ELSE
            DO J1=1,3*NATOMS
               P(J1)=PP(J1)
            ENDDO
         ENDIF

         IF (.NOT.ITINC.AND.(.NOT.(IFLAG.LE.0))) THEN
            GOTO 11
         ENDIF
         FIXIMAGE=.FALSE.

         IF (IFLAG.EQ.0.AND.(.NOT.DZT)) THEN
            CFLAG=.TRUE.
            IF (AMBER.AND.DEBUG) THEN
               PRINT*,'DIAG:'
               WRITE(*,'(12F11.3)') (DIAG(J2),J2=1,3*NATOMS)
            ENDIF
            RETURN
         ENDIF

         IF ((IFLAG.LT.0).OR.(INFO.EQ.6)) THEN
            PRINT*,'WARNING - after LBFGS IFLAG=',IFLAG,' energy=',EREAL
            IFLAG=0
            INFO=0
            IF (AMBER.AND.(RMS.GT.1.0D0)) THEN
C              EREAL=1.0D6
C              POTEL=1.0D6
C              RMS=1.0D0
C              WRITE(*,'(A)') ' Stuck - step discarded'
C              RETURN
C
C  Short range Morse potentials can get stuck through atoms getting too far out
C  of the core. Try contracting.
C
            ELSE IF (MORSET.AND.RHO.GT.6.0D0) THEN
               DO J1=1,NATOMS-NSEED
                  P(3*(J1-1)+1)=P(3*(J1-1)+1)*0.9D0   
                  P(3*(J1-1)+2)=P(3*(J1-1)+2)*0.9D0   
                  P(3*(J1-1)+3)=P(3*(J1-1)+3)*0.9D0   
               ENDDO
            ELSE IF (CENT.AND.(.NOT.SEEDT)) THEN
               CALL CENTRE2(P)
            ELSE IF (FIXCOM.AND.(.NOT.SEEDT)) THEN
               CALL CENTRECOM(P)
            ELSE
               CFLAG=.FALSE.
C              RETURN
            ENDIF
         ENDIF

         IF (AMBER.AND.DEBUG.AND.DUMPT) THEN
            WRITE(40,'(I4)') NATOMS
            WRITE(40,'(A,I4,A,F15.5)') 'At step number ',ITER,' energy=',EREAL
            DO J2=1,NATOMS
               WRITE(40,'(A,3F20.10)') typech(J2)(1:1),(P(3*(J2-1)+J3),J3=1,3)
            ENDDO
         ENDIF
         CALL FLUSH(6)
      ENDDO

      CFLAG=.FALSE.

      RETURN
      END
