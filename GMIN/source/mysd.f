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
      SUBROUTINE MYSD(ITMAX,VARS,MFLAG,NSTP,ENERGY)
      USE COMMONS,ONLY : NATOMS, MYUNIT, GMAX, FIXIMAGE, RMS, SDTOL, DEBUG
      IMPLICIT NONE
      INTEGER NSTP, ITMAX, J1, NOPT
      LOGICAL MFLAG
      DOUBLE PRECISION VARS(3*NATOMS), GRAD(3*NATOMS), ENERGY, STEP(3*NATOMS), NEWVARS(3*NATOMS),
     &                 EPRED, ENEW, NEWGRAD(3*NATOMS), SFAC, GNORM, SLENGTH, DDOT, PERROR

      SFAC=1.1D0
      NOPT=3*NATOMS
      NSTP=1
      CALL POTENTIAL(VARS,GRAD,ENERGY,.TRUE.,.FALSE.)

10    GNORM=DSQRT(DDOT(NOPT,GRAD,1,GRAD,1))
      IF (GNORM.LE.0.0D0) THEN
         WRITE(MYUNIT,'(A)') 'mysd> ERROR - GNORM is zero'
         STOP
      ENDIF
      IF (NSTP.EQ.1) THEN
         SLENGTH=MIN(GNORM,1.0D0/GNORM)
      ENDIF

      !FIXIMAGE=.TRUE.
!
!  Check for convergence on RMS force.
!
      IF (RMS.LT.GMAX) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            RETURN
         ENDIF
      ENDIF
!
!  Check if we have exceeded maximum allowed iterations.
!
      IF (NSTP.EQ.ITMAX) THEN
         MFLAG=.FALSE.
         FIXIMAGE=.FALSE.
         RETURN
      ENDIF

11    DO J1=1,NOPT
         STEP(J1)=-SLENGTH*GRAD(J1)/GNORM
         NEWVARS(J1)=VARS(J1)+STEP(J1)
      ENDDO
      EPRED=DDOT(NOPT,GRAD,1,STEP,1)
 
      CALL POTENTIAL(NEWVARS,NEWGRAD,ENEW,.TRUE.,.FALSE.)
      IF (ENEW-ENERGY.EQ.0.0D0) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'mysd> WARNING - ENEW=ENERGY=',ENERGY
         PERROR=0.0D0
      ELSE 
         PERROR=(EPRED-(ENEW-ENERGY))*100.0D0/(ENEW-ENERGY)
      ENDIF
      IF (PERROR.GT.SDTOL) THEN
         SLENGTH=SLENGTH/SFAC
         IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,G20.10)') 'mysd> % error=',PERROR,
     &     ' exceeds tolerance, decreasing step length to',SLENGTH
         GOTO 11
      ELSEIF (PERROR.GT.SDTOL*0.66D0) THEN
         SLENGTH=SLENGTH/SFAC
         IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,G20.10)') 'mysd> % error=',PERROR,
     &     ' decreasing step length to',SLENGTH
      ELSE
         SLENGTH=SLENGTH*SFAC
         IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,G20.10)') 'mysd> % error=',PERROR,
     &     ' increasing step length to',SLENGTH
      ENDIF
      ENERGY=ENEW
      GRAD(1:NOPT)=NEWGRAD(1:NOPT)
      VARS(1:NOPT)=NEWVARS(1:NOPT)

      IF (DEBUG) WRITE(MYUNIT,'(A,2F20.10,A,I6,A,G15.5)') 'mysd> E and RMS=',ENERGY,RMS,' after ',
     &        NSTP,' SD cycles, step=',SLENGTH

      IF (SLENGTH.LT.1.0D-200) THEN
         WRITE(MYUNIT,'(A)') 'mysd> Step size underflow - quit'
         STOP
      ENDIF
      MFLAG=.FALSE.
      FIXIMAGE=.FALSE.
      NSTP=NSTP+1
      GOTO 10

      RETURN
      END
