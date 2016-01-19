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
      SUBROUTINE DFPMIN(ITMAX,P,N,GTOL,ITER,FRET,BFSUCCESS)
      USE commons
C      use modamber
      IMPLICIT NONE
      INTEGER ITER,N,ITMAX,ICRAP
      DOUBLE PRECISION FRET,GTOL,P(N),EPS,STPMX,TOLX,prms
      PARAMETER (STPMX=1.d1,EPS=3.0D-8,TOLX=4.0D0*EPS)
      LOGICAL BFSUCCESS,BSTUCK,CRAP,CTEST,NOTCALLED,DFLAG
      INTEGER I,ITS,J
      LOGICAL CHECK
      DOUBLE PRECISION FAC,FAD,FAE,FP,STPMAX,SUM,SUMDG,SUMXI,DG(3*NATOMS),
     1                 G(3*NATOMS),HDG(3*NATOMS),HESSIN(3*NATOMS,3*NATOMS),PNEW(3*NATOMS),XI(3*NATOMS)

C     FP=F1DIM(P)

      ICRAP=1
      CRAP=.FALSE.
      NOTCALLED=.TRUE.
      CTEST=.FALSE.
      DFLAG=DEBUG
C
C  Get initial gradient and energy.
C
      CALL POTENTIAL(P,G,FP,.TRUE.,.FALSE.)
      SUM=0.0D0
C
C  Initialize Hessian to unit matrix. Initial search direction, XI,
C  is antiparallel to the gradient.
C
      DO I=1,N
         DO J=1,N
            HESSIN(J,I)=0.0D0
         ENDDO
         HESSIN(I,I)=1.d0
         XI(I)=-G(I)
         SUM=SUM+P(I)**2
      ENDDO
C     STPMAX=STPMX*MAX(SQRT(SUM), DBLE(N))
      STPMAX=MAXBFGS

      DO ITS=1,ITMAX
!        IF (MOD(ITS,20).EQ.0) CALL MAKELIST ! Paul Mortenson's amber specific!
         ITER=ITS
         CALL LNSRCH(N,P,FP,G,XI,PNEW,FRET,STPMAX,CHECK,BSTUCK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        DO I=1,N
C           PNEW(I)=P(I)
C        ENDDO
C        CALL LINMIN(ITS,PNEW,G,NATOMS,FRET)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (BSTUCK) THEN
            WRITE (*,*) ' *** Stuck in lnsrch - taking 5 CG steps'
            DEBUG=.TRUE.
            CALL CGMIN(5,P,CRAP,ICRAP,FRET,1)
            IF (.NOT.(DFLAG)) DEBUG=.FALSE.
         ENDIF
         FP=FRET
         DO I=1,N
           XI(I)=PNEW(I)-P(I)
           P(I)=PNEW(I)
           DG(I)=G(I)
         ENDDO
C
C  Get new gradient and Hessian.
C
         CALL POTENTIAL(P,G,FP,.TRUE.,.FALSE.)
         PRMS=0.0D0
         DO I=1,N
            PRMS=PRMS+G(I)**2
         ENDDO
         PRMS=SQRT(PRMS/N)
         IF (DEBUG) WRITE (*,*) ' Iteration ',ITS,' energy= ',FP,' RMS force= ',PRMS
C        IF (DEBUG) CALL FLUSH(6)
         IF (PRMS.LT.GTOL) THEN
            BFSUCCESS=.TRUE.
            RETURN
         ENDIF
C
C  Catch cold fusion for ionic potentials and discard.
C
C  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
C  for systems with positive energies. - khs26 26/11/09
C
         IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.PACHECO).AND.(FRET.LT.-1.0D4)) THEN
            FRET=1.0D4
            RMS=1.0D0
            WRITE(*,'(A)') ' Cold fusion diagnosed - step discarded'
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
            COLDFUSION=.TRUE.
            RETURN
         ENDIF
!        IF ((AMBER.AND.NOTCALLED).AND.(PRMS.LT.1.0D0).AND.(.NOT.(FIX))) THEN
!           CALL CHIRALTEST(CTEST,P)
!           IF (CTEST) THEN
!              WRITE(*,'(A)') ' Change in chirality detected - step rejected'
!              FRET=1.0D6
!              RMS=1.0D0
!              RETURN
!           ENDIF
!           NOTCALLED=.FALSE.
!        ENDIF

         DO I=1,N
            DG(I)=G(I)-DG(I)
         ENDDO
         DO I=1,N
            HDG(I)=0.0D0
            DO J=1,N
               HDG(I)=HDG(I)+HESSIN(J,I)*DG(J)
            ENDDO
         ENDDO
C
C  Hessian update.
C
         FAC=0.0D0
         FAE=0.0D0
         SUMDG=0.0D0
         SUMXI=0.0D0
         DO I=1,N
            FAC=FAC+DG(I)*XI(I)
            FAE=FAE+DG(I)*HDG(I)
            SUMDG=SUMDG+DG(I)**2
            SUMXI=SUMXI+XI(I)**2
         ENDDO
         IF (FAC**2.GT.EPS*SUMDG*SUMXI) THEN
            FAC=1.0D0/FAC
            FAD=1.0D0/FAE
            DO I=1,N
               DG(I)=FAC*XI(I)-FAD*HDG(I)
            ENDDO
            DO I=1,N
               DO J=I,N
                  HESSIN(J,I)=HESSIN(J,I)+FAC*XI(I)*XI(J)-FAD*HDG(I)*HDG(J)+FAE*DG(I)*DG(J)
                  HESSIN(I,J)=HESSIN(J,I)
               ENDDO
            ENDDO
         ENDIF
C
C  Search direction update.
C
         DO I=1,N
            XI(I)=0.0D0
            DO J=1,N
              XI(I)=XI(I)-HESSIN(J,I)*G(J)
            ENDDO
         ENDDO
      ENDDO

      BFSUCCESS=.FALSE.
      RETURN
      END
