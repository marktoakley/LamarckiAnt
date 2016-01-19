!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
!        Line search removed plus small modifications, DJW 2001
!
      SUBROUTINE MYMYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)
      USE COMMONS
      USE GENRIGID
      USE MODAMBER
      USE MODAMBER9, ONLY : STEEREDMINT, LOCALSTEEREDMINT, SMINK, SMINKINC, SMINKCURRENT
      USE MODCHARMM
      USE PORFUNCS
      IMPLICIT NONE

!hk286 - introduce XRIGIDCOORDS
      DOUBLE PRECISION XRIGIDCOORDS(3*NATOMS)

      INTEGER N,M,J1,ITMAX,ITDONE,NP,J2,J3,NFAIL,NDECREASE,NGUESS,NDUMMY
      DOUBLE PRECISION XCOORDS(3*NATOMS),GRAD(3*NATOMS),SLENGTH,DDOT,EPLUS,EMINUS,DIFF,DUMMY,WTEMP(3*NATOMS)
      DOUBLE PRECISION TMPANG(3*NATOMS), TMPCOORDS(3*NATOMS), RIGIDDUMPCOORDS(3*NATOMS)
      DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(3*NATOMS),OVERLAP,OLDX(3*NATOMS),OLDOLDX(3*NATOMS),VGUESS(3),
     1                 X1, Y1, Z1, X2, Y2, Z2, TRY(3*NATOMS), D1, D2, RBCOORDS(18), DUMMY2, DIST, DIST1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
!     DOUBLE PRECISION DIAG(N), W(N*(2*M+1)+2*M)
      LOGICAL DIAGCO, YESNO, RESET, NOTCALLED, CTEST, MFLAG, LRESET
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,POTEL,QSTART,QFINISH
      DOUBLE PRECISION OLDCART(3*NATOMS), DELTAQ(N),DELTACART(3*NATOMS),LEPSILON,DOT1,DOT2
      DOUBLE PRECISION LCART(3*NATOMS),OLDQ(N),NEWQ(N),OLDGINT(N),GINT(N),XINT(N),XSAVE(N),SMINKCURRENTP
      DOUBLE PRECISION, ALLOCATABLE :: FRAMES(:,:), PE(:), MODGRAD(:)
      LOGICAL NOCOOR, FAILED, COREDONE
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,INMC,IYCN,ISCN
      INTEGER KD, NNZ, GETUNIT
      COMMON /MYPOT/ POTEL
      COMMON /Q4C/ QSTART, QFINISH
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, SMINKCHANGET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAP, evapreject
      LOGICAL FTEST
      COMMON /FAIL/ FTEST
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT
      LOGICAL RESETLBFGS
      COMMON /CRESETLBFGS/ RESETLBFGS
!     SAVE 

      LOGICAL CELL_HAS_CHANGED ! vr274> test for cell reduction

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for pgi and ifort
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for pgi and ifort
!     IF (QUENCHDOS) ALLOCATE(FRAMES(N,ITMAX), PE(ITMAX), MODGRAD(ITMAX))
      LRESET=.FALSE.
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         IF (GCBHT) THEN
            DEALLOCATE(W)
            ALLOCATE(W(N*(2*M+1)+2*M))
            DEALLOCATE(DIAG)
            ALLOCATE(DIAG(N))  
!
! We cannot change RESET itself - it gives a segmentation fault because MYLBFGS is called
! with fixed values in this argument, not variables.
!
            LRESET=.TRUE.
         ELSE
            WRITE(MYUNIT, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
            call EXIT(10)
         ENDIF
      ENDIF
      COREDONE=.FALSE.
      SMINKCHANGET=.FALSE.
      SMINKCURRENT=0.0D0
      SMINKCURRENTP=0.0D0
      LOCALSTEEREDMINT=.FALSE.
      IF (STEEREDMINT) LOCALSTEEREDMINT=.TRUE.

! for CHARMM: update nonbonded list at the start of each minimization
! hk286 - revised since XCOORDS is actually rigid body coordinates here
      IF(CHRMMT) THEN
         IF (ATOMRIGIDCOORDT) THEN
            CALL UPDATENBONDS(XCOORDS)
         ELSE
            XRIGIDCOORDS(1:DEGFREEDOMS) = XCOORDS(1:DEGFREEDOMS)
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
            CALL UPDATENBONDS(XCOORDS)  
            CALL TRANSFORMCTORIGID(XCOORDS, XRIGIDCOORDS)
            XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
            XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
         ENDIF
      ENDIF
!hk286

      IF(DMACRYST) CELL_HAS_CHANGED=.FALSE.

      NFAIL=0
      IF (GUIDECHANGET) ITER=0
      IF (RESET.OR.LRESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.
      IF (DEBUG) THEN
         IF (LRESET.OR.RESET.OR.GUIDECHANGET) WRITE(MYUNIT,'(A)') 'mylbfgs> Resetting LBFGS minimiser'
         IF (.NOT.(LRESET.OR.RESET.OR.GUIDECHANGET)) WRITE(MYUNIT,'(A)') 'mylbfgs> Not resetting LBFGS minimiser'
!fh301>{{{
         IF (CHEMSHIFT) THEN
             OPEN(UNIT=137,FILE='rmsshifts',STATUS='old',ACTION='write',
     &       position='append')
             WRITE(137,'(I6,6(F15.10))') ITDONE,RMSSHIFT(:)
             WRITE(137,'(A,I6)') 'Quench',NQ(NP)
             WRITE(137,'(A)') '  Qu        HA             CA             HN             N              C              CB'
             WRITE(137,'(A)') '------------------------------------------------------------------------------------------------'
             CLOSE(137)
         ENDIF
!fh301>}}}
      ENDIF

      IF (CALCQT) CALL ORDERQ4(NATOMS,XCOORDS,QSTART)

      IF (DUMPT) THEN
         IF (ARNO) THEN
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS+2
            WRITE(DUMPXYZUNIT(NP),11) NP,NQ(NP)
            WRITE(DUMPXYZUNIT(NP),'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(DUMPXYZUNIT(NP),'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(DUMPXYZUNIT(NP),65) (XCOORDS(J1),J1=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TIP) THEN
            WRITE(DUMPXYZUNIT(NP),'(I6)') (NATOMS/2)*3
            WRITE(DUMPXYZUNIT(NP),'(A,I5)') 'LBFGS iteration ',ITER
            DO J2=1,NATOMS/2
               CALL TIPIO(XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),
     1              XCOORDS(3*(NATOMS/2+J2-1)+1),XCOORDS(3*(NATOMS/2+J2-1)+2),XCOORDS(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ELSEIF (NCORE(NP).GT.0) THEN
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
            WRITE(DUMPXYZUNIT(NP),11) NQ(NP)
            WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') ('LB',XCOORDS(3*(I-1)+1),XCOORDS(3*(I-1)+2),XCOORDS(3*(I-1)+3),
     &                                               I=1,NATOMS-NCORE(NP))
            IF (NCORE(NP).GT.0) WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     &            ('LA ',XCOORDS(3*(I-1)+1),XCOORDS(3*(I-1)+2),XCOORDS(3*(I-1)+3),I=NATOMS-NCORE(NP)+1,NATOMS)
         ELSE
! csw34> Transform rigid coordinates to atomistic before dumping so we can visualise
            IF (.NOT.ATOMRIGIDCOORDT) THEN
               RIGIDDUMPCOORDS(1:3*NATOMS) = 0.0D0
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, RIGIDDUMPCOORDS, XCOORDS)
               WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
               WRITE(DUMPXYZUNIT(NP),11) NQ(NP)
               WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') ('LA ',RIGIDDUMPCOORDS(3*(J1-1)+1),
     &             RIGIDDUMPCOORDS(3*(J1-1)+2),RIGIDDUMPCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
               IF (NS.GT.0) WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     1             ('LB',RIGIDDUMPCOORDS(3*(J1-1)+1),RIGIDDUMPCOORDS(3*(J1-1)+2),RIGIDDUMPCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
            ELSE
               WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
               WRITE(DUMPXYZUNIT(NP),11) NQ(NP)
               WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),
     &             XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
               IF (NS.GT.0) WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     1             ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
11             FORMAT(1X,'QUENCH NUMBER ',I6,' initial points in mymylbfgs')
            ENDIF
         ENDIF
      ENDIF
      CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
!
!  Catch cold fusion for ionic potentials and discard.
!
!  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
!  for systems with positive energies. - khs26 26/11/09
!
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT.OR.SILANET 
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.PYT.OR.JMT.OR.DMACRYST.OR.TTM3T.OR.AMBER12T.OR.MBPOLT)
     &   .AND.((ENERGY.LT.COLDFUSIONLIMIT).OR.PYCOLDFUSION)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENERGY,COLDFUSIONLIMIT
         ENERGY=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
         COLDFUSION=.TRUE.
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
         RETURN
      ENDIF
      IF (FTEST) THEN
         ENERGY=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Diagonalisation failure - step discarded'
         RETURN
      ENDIF

!     IF (QUENCHDOS) THEN
!        MODGRAD(1)=DSQRT(DDOT(N,GRAD,1,GRAD,1))
!        PE(1)=ENERGY
!        FRAMES(1:N,1)=XCOORDS(1:N)
!     ENDIF
!
! Stop the core from changing morphology easily, but allow it to relax
!
!     IF ((NCORE(NP).GT.0).AND.(.NOT.COREDONE)) THEN
!        DUMMY2=0.0D0
!        DO J1=1,3*(NATOMS-NCORE(NP))
!           DUMMY2=DUMMY2+GRAD(J1)**2
!        ENDDO
!        DUMMY2=DSQRT(DUMMY2/(3*NATOMS))
!        IF (DUMMY2.GT.EPS*2.0D0) THEN
!           GRAD(3*(NATOMS-NCORE(NP))+1:3*NATOMS)=0.0D0
!        ELSE
!           COREDONE=.TRUE.
!        ENDIF
!     ENDIF

      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - geometry discarded.'
!         IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
          RETURN
      ENDIF

!
!  If INTMINT and CHRMMT need to transform to internal coordinates
!  See COPTIM.2.3 for switching to internals from Cartesians using LIMINCUT.
!
      IF (INTMINT) THEN
         OLDCART(1:3*NATOMS)=XCOORDS(1:3*NATOMS) ! store cartesians in OLDCART for both CHARMM and UNRES
!         IF (UNRST) THEN
!C
!C store internals (in OLDQ) and update X to contain internals
!C
!            CALL geom_to_var(N,OLDQ)
!            XCOORDS(1:N)=OLDQ(1:N)
!         ELSE IF (CHRMMT) THEN
            CALL GETKD(KD) ! get width of sparse band in G matrix KD
            CALL GETNNZ(NNZ) ! get number of non-zero elements in B-matrix
            NOCOOR=.FALSE. ! calculate internals therefore NOCOOR is false
            GINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            XINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            CALL TRANSFORM(XCOORDS,GRAD,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,KD)
            OLDQ(1:N)=XINT(1:N)    ! store internals
            OLDGINT(1:N)=GINT(1:N) ! store gradient in internals
!         ENDIF
      ENDIF
!
!  for CHRMMT:
!  XCOORDS contains current Cartesians
!  GRAD    contains current gradient
!  XINT    contains current internals
!  GINT    contains current gradient in internals
!  OLDQ    contains internals for initial geometry
!  OLDGINT contains gradient in internals for initial geometry
!  OLDCART contains Cartesian coordinates for initial geometry
!
      IF (EVAPREJECT) RETURN
      POTEL=ENERGY

      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

!
!  Catch cold fusion for ionic potentials and discard.
!
      IF ((DBPT.OR.DBPTDT.OR.DMBLPYT.OR.MSTBINT.OR.MSSTOCKT.OR.MULTPAHAT.OR.NPAHT.OR.PAHW99T.OR.TDHDT.OR.SILANET)
     &    .AND.(ENERGY.LT.-5.0D4)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENERGY,COLDFUSIONLIMIT
         ENERGY= 1.D06
         POTEL = 1.D06
         RMS   = 1.0D0
         COLDFUSION = .TRUE.
         RETURN
      ENDIF

!
!  Termination test. 
!
10    CALL FLUSH(MYUNIT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN ! vr274> added cell change test for DMACRYS
         IF (CHRMMT.AND.ACESOLV) THEN
            NCHENCALLS=ACEUPSTEP-1
            CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
            PRINT *, 'HERE'
! hk286
!            PRINT *, "ENERGY ", ENERGY

            IF (DEBUG) WRITE(*,'(A,2G20.10,A)') ' mylbfgs> Energy and RMS force=',ENERGY,RMS,' after ACE update'
            IF (RMS.LE.EPS) MFLAG=.TRUE.
         ELSE
            MFLAG=.TRUE.
         ENDIF
         ! vr274>
        IF(DMACRYST) THEN
          CALL DMACRYS_REDUCE_CELL_RIGID(NATOMS, XCOORDS,CELL_HAS_CHANGED)
          IF(CELL_HAS_CHANGED) THEN
              MFLAG=.false.
              write(MYUNIT,*) "cell was reduced"
          ENDIF
        ENDIF


         IF (EVAP) MFLAG=.FALSE. ! do not allow convergence if we happen to have a small RMS and EVAP is true'
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

!             IF (QUENCHDOS) THEN
!                DO J1=1,ITDONE+1
!                   DIST=0.0D0
!                   DO J2=1,N
!                      DIST=DIST+(FRAMES(J2,J1)-FRAMES(J2,ITDONE+1))**2
!                   ENDDO
!                   DIST=SQRT(DIST)
!                   IF (J1.EQ.1) DIST1=DIST
! !                   DOSSTATS(J1,1)=PE(J1)
!                   IF ((MODGRAD(J1).GT.0.0D0).AND.(DIST.GT.0.0D0)) THEN
! !                    DOSSTATS(J1,2)=(MODGRAD(1)/MODGRAD(J1))*(DIST/DIST1)**(N-1)
! !                    DOSSTATS(J1,2)=(N-1)*LOG(DIST)-LOG(MODGRAD(J1))
!                      DOSSTATS(J1,2)=DIST**(N-1)/MODGRAD(J1)
!                   ELSE
!                      DOSSTATS(J1,2)=0.0D0
!                   ENDIF
! !                 WRITE(MYUNIT,'(A,I6,4G18.8)') 'lbfgs> J1,MODGRAD,DIST,DOSSTATS(J1,2),DOSSTATS(J1,1)=',
! !    &                      J1,MODGRAD(J1),DIST,DOSSTATS(J1,2),DOSSTATS(J1,1)
!                ENDDO
!                DEALLOCATE(FRAMES, PE, MODGRAD)
!             ENDIF
            
            RETURN
         ENDIF
      ENDIF
      

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) FIXIMAGE=.FALSE.
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
!fh301>{{{
         IF (CHEMSHIFT2) THEN
           CHEMSHIFT=.TRUE.
!           CHEMSHIFTITER=0
         ENDIF
!fh301>}}}
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
         RETURN
      ENDIF


      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(MYUNIT,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(MYUNIT,235) J1
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
!           INQUIRE(FILE='diag',EXIST=YESNO)
!           IF (YESNO) THEN
!              OPEN(UNIT=34,FILE='diag',STATUS='OLD')
!              READ(34,*) (DIAG(J1),J1=1,N)
!              PRINT*,'diag read in LBFGS'
!              WRITE(*,'(6F15.5)') (DIAG(J1),J1=1,N)
!           ELSE
            DO J1=1,N
               DIAG(J1)=DGUESS
            ENDDO
         ENDIF

!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY INFORMATION.
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!
         ISPT= N+2*M    ! index for storage of search steps
         IYPT= ISPT+N*M ! index for storage of gradient differences
!
!  NR step for diagonal inverse Hessian
!
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(ISPT+I)= -GINT(I)*DIAG(I)
               W(I)= -GINT(I)*DIAG(I)
            ENDDO
            GNORM= DSQRT(DDOT(N,GINT,1,GINT,1))
         ELSE
            DO J1=1,N
               DUMMY=-GRAD(J1)*DIAG(J1)
               W(ISPT+J1)=DUMMY
               W(J1)=DUMMY
            ENDDO
            GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))
         ENDIF
!
!  Make the first guess for the step length cautious.
!
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!
!  Update estimate of diagonal inverse Hessian elements
!
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (ABS(YY).LT.1.0D-20) THEN
               YY=SIGN(1.0D-20,YY)
               IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A)') 'WARNING, resetting YY to ',YY,' in mymylbfgs'
            ENDIF
            IF (ABS(YS).LT.1.0D-20) THEN
               YS=SIGN(1.0D-20,YS)
               IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A)') 'WARNING, resetting YS to ',YS,' in mymylbfgs'
            ENDIF
!           WRITE(MYUNIT,'(A,2G20.10)') 'YS/YY=',YS/YY
            DO J1=1,N
!              DIAG(J1)= ABS(YS/YY) ! messes up after step reversals!
               DIAG(J1)= YS/YY
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(MYUNIT,235) J1
                  STOP
               ENDIF
            ENDDO
         ENDIF
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!

         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(I)= -GINT(I)
            ENDDO
         ELSE
            DO J1=1,N
               W(J1)= -GRAD(J1)
            ENDDO
         ENDIF
         CP= POINT
         DO J1= 1,BOUND
            CP=CP-1
            IF (CP.EQ.-1) CP=M-1
            SQ=DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)=W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO J1=1,N
            W(J1)=DIAG(J1)*W(J1)
            IF (ABS(W(J1)).GT.1.0D2) W(J1)=SIGN(1.0D2,W(J1)) ! DJW
         ENDDO

         DO J1=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0  
      ENDIF
!
!  Store the new search direction
!
      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

      IF (CHRMMT.AND.INTMINT) THEN
         DOT1=SQRT(DDOT(N,GINT,1,GINT,1))
      ELSE
         DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
      ENDIF
!
!  Overflow has occasionally occurred here.
!  We only need the sign of the overlap, so use a temporary array with
!  reduced elements.
!
      DUMMY=1.0D0
      DO J1=1,N
         IF (ABS(W(J1)).GT.DUMMY) DUMMY=ABS(W(J1))
      ENDDO
      DO J1=1,N
         WTEMP(J1)=W(J1)/DUMMY
      ENDDO
      DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            OVERLAP=DDOT(N,GINT,1,WTEMP,1)/(DOT1*DOT2)
         ELSE
            OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
        ENDIF
      ENDIF
!     WRITE(MYUNIT,'(A,2G20.10)') 'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
!     WRITE(MYUNIT,'(A,2G20.10)')'GRAD . GRAD=',DDOT(N,GRAD,1,GRAD,1)
!     WRITE(MYUNIT,'(A,2G20.10)') 'WTEMP . WTEMP=',DOT2
      IF (OVERLAP.GT.0.0D0) THEN
!        IF (DEBUG) PRINT*,'Search direction has positive projection onto gradient - resetting'
!        ITER=0
!        GOTO 10
         IF (DEBUG) WRITE(MYUNIT,'(A)') 'Search direction has positive projection onto gradient - reversing step'
         DO J1=1,N
            W(ISPT+POINT*N+J1)= -W(J1)  !!! DJW, reverses step
            WTEMP(J1)= -W(J1)/DUMMY  !!! DJW, reverses step
         ENDDO
      ENDIF

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(I)=GINT(I)
         ENDDO
      ELSE
         DO J1=1,N
            W(J1)=GRAD(J1)
         ENDDO
      ENDIF
      SLENGTH=0.0D0
      DO J1=1,N
!        SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
         SLENGTH=SLENGTH+WTEMP(J1)**2
      ENDDO
!
! Taking out the magnitude of the largest element can prevent overflow.
!
      SLENGTH=SQRT(SLENGTH)*DUMMY
!     WRITE(MYUNIT,'(A,2G20.10)') 'SLENGTH=',SLENGTH
      IF (STP*SLENGTH.GT.MAXBFGS) THEN
!         WRITE(MYUNIT,'(A,3G30.15)') 'STP,SLENGTH,MAXBFGS=',STP,SLENGTH,MAXBFGS
         STP=MAXBFGS/SLENGTH
      ENDIF


!
!  We now have the proposed step.
!
      IF (CHRMMT.AND.INTMINT) THEN
         DO J1=1,N
            XINT(J1)=XINT(J1)+STP*W(ISPT+POINT*N+J1)
            DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)
         ENDDO
      ELSE
!
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
         XSAVE(1:N)=XCOORDS(1:N) 
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
!
! For Thomson try projection for the geometry after the step.
!
         IF (PROJIT) THEN
            IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJI(TMPCOORDS,NATOMS)
               CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJI(TMPCOORDS,NATOMS)
               CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL PROJI(TMPCOORDS,NATOMS)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
         IF (PROJIHT) THEN
            IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJIH(TMPCOORDS,NATOMS)
               CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJIH(TMPCOORDS,NATOMS)
               CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL PROJIH(TMPCOORDS,NATOMS)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
      ENDIF
!
!  For charmm internals must transform and back-transform!
!
      NDECREASE=0
      LEPSILON=1.0D-6

20    IF (INTMINT) THEN
         IF (CHRMMT) THEN
            NEWQ(1:N)=OLDQ(1:N)
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)
!
! Need to keep OLDQ constant for repeated back-transformations if first step size fails.
! Therefore pass dummy array newq that can change.
! Similarly with LCART and OLDCART.
!
!           CALL TRANSBACK(XINT,NEWQ,LCART,N,3*NATOMS,NNZ,KD)
            CALL TRANSBACKDELTA(DELTAQ,DELTACART,LCART,N,3*NATOMS,NNZ,KD,FAILED,.FALSE.,LEPSILON) ! transform step to Cartesians
            IF (FAILED) THEN
!              NCOUNT=NCOUNT+1
!              IF (NCOUNT.GT.1) STOP
!              LEPSILON=1.0D-5*DPRAND()
!              GOTO 21
! or
               MFLAG=.FALSE.
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
            ENDIF
!
! now add DELTACART to LCART to get new cartesians. Put these in X.
!
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
            XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
!
!  for CHRMMT:
!  LCART    contains new Cartesians (after step)
!  XCOORDS contains new Cartesians (after step)
!  XINT    contains new internals (after step)
!  GRAD    contains old gradient
!  GINT    contains old gradient in internals
!  OLDQ    contains old internals
!  OLDGINT contains old gradient in internals for the last successful geometry
!  NEWQ    contains old internals for the last successful geometry
!  OLDCART contains old Cartesians for the last successful geometry
!
!         ELSEIF (UNRST) THEN
!            NEWQ(1:N)=X(1:N) ! store new internals in NEWQ
!C
!C need a temporary array NEWQ here as argument to var_to_geom to keep X unchanged in case we need to
!C modify the step below.
!C
!            CALL var_to_geom(N,NEWQ) ! update internals
!            CALL chainbuild ! get cartesians
         ENDIF
      ENDIF

! csw34> INCREMENT THE FORCE CONSTANT FOR STEERED MINIMISATION
      SMINKCHANGET=.FALSE.
      IF (LOCALSTEEREDMINT) THEN
         SMINKCURRENT=MIN(SMINKCURRENT+SMINKINC,SMINK)
         IF (SMINKCURRENT.NE.SMINKCURRENTP) SMINKCHANGET=.TRUE.
! a bit of useful debug printing         
        IF (DEBUG) WRITE(MYUNIT,'(A,2F20.10,L5)') 'SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET=',SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET
      ENDIF

      CALL POTENTIAL(XCOORDS,GNEW,ENEW,.TRUE.,.FALSE.)
      IF(RESETLBFGS) THEN
         RESETLBFGS=.FALSE.
         IF(DEBUG) WRITE(MYUNIT,*) "Reset l-bfgs when Charmm updating nonbond list"
         ITER=0
         ENERGY=ENEW
      ENDIF
      IF (FTEST) THEN
         ENERGY=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Diagonalisation failure - step discarded'
         RETURN
      ENDIF
! hk286
      !PRINT*, "ENEW ", ENEW
      !PRINT*, XCOORDS
      !PRINT*, GNEW
      !READ*, DUMMY2

!     IF (QUENCHDOS) THEN
!        MODGRAD(ITDONE+2)=DSQRT(DDOT(N,GNEW,1,GNEW,1))
!        PE(ITDONE+2)=ENEW
!        FRAMES(1:N,ITDONE+2)=XCOORDS(1:N)
!     ENDIF
!
! Stop the core from changing morphology easily, but allow it to relax
!
!     IF ((NCORE(NP).GT.0).AND.(.NOT.COREDONE)) THEN
!        DUMMY2=0.0D0
!        DO J1=1,3*(NATOMS-NCORE(NP))
!           DUMMY2=DUMMY2+GNEW(J1)**2
!        ENDDO
!        DUMMY2=DSQRT(DUMMY2/(3*NATOMS))
!        IF (DUMMY2.GT.EPS*2.0D0) THEN
!           GNEW(3*(NATOMS-NCORE(NP))+1:3*NATOMS)=0.0D0
!        ELSE
!           COREDONE=.TRUE.
!        ENDIF
!     ENDIF

      IF (EVAPREJECT) return
      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - step discarded.'
!         IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
          RETURN
      ENDIF
!
!  Catch cold fusion for ionic potentials and discard.
!
!  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
!  for systems with positive energies. - khs26 26/11/09
!
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT.OR.SILANET
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.PYT.OR.JMT.OR.DMACRYST.OR.TTM3T.OR.AMBER12T.OR.MBPOLT)
     &   .AND.((ENEW.LT.COLDFUSIONLIMIT).OR.PYCOLDFUSION)) THEN
         WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENEW,COLDFUSIONLIMIT
         ENERGY=1.0D6
         ENEW=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
         COLDFUSION=.TRUE.
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
         RETURN
      ENDIF
      IF ((DBPT.OR.DBPTDT.OR.DMBLPYT.OR.MSTBINT.OR.MSSTOCKT.OR.MULTPAHAT.OR.NPAHT.OR.PAHW99T.OR.TDHDT.OR.SILANET)
     &    .AND.(ENEW.LT.-5.0D4)) THEN
         WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENEW,COLDFUSIONLIMIT
         ENERGY=0.0D0
         ENEW=0.0D0
         POTEL=0.0D0
         RMS=1.0D0
         RETURN
      ENDIF


!
!  We need to transform the newly obtained Cartesian gradient for CHARMM and internals.
!  NOCOOR is true because we dont need to transform the coordinates.
!
      IF (CHRMMT.AND.INTMINT) THEN
         NOCOOR=.TRUE.
         CALL TRANSFORM(XCOORDS,GNEW,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,KD)
      ENDIF

!     IF (TIP) THEN
!           WRITE(DUMPXYZUNIT(NP),'(I6)') (NATOMS/2)*3
!           WRITE(DUMPXYZUNIT(NP),'(A,I5,A,F20.10)') 'LBFGS iteration ',ITER,' energy =',ENEW
!           DO J2=1,NATOMS/2
!              CALL TIPIO(XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),
!    1              XCOORDS(3*(NATOMS/2+J2-1)+1),XCOORDS(3*(NATOMS/2+J2-1)+2),XCOORDS(3*(NATOMS/2+J2-1)+3),RBCOORDS)
!              WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
!              WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
!              WRITE(DUMPXYZUNIT(NP),'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
!           ENDDO
!     ENDIF

!     WRITE(*,'(A,F20.10)') 'ENEW=',ENEW
!     WRITE(*,'(I6,F20.10)') (J1,GNEW(J1),J1=1,N)

!
! csw34 Force acceptance of step if FIXDIHEFLAG is TRUE
!
      IF (FIXDIHEFLAG) ENERGY=ENEW

      IF(DMACRYST.AND.CELL_HAS_CHANGED) THEN ! force acceptance after lattice reduction
          ENERGY=ENEW
          CELL_HAS_CHANGED=.FALSE.
      ENDIF

      IF (((ENEW-ENERGY.LE.MAXERISE).OR.EVAP.OR.GUIDECHANGET.OR.SMINKCHANGET).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,3*NATOMS
            GRAD(J1)=GNEW(J1)
         ENDDO
!fh301>{{{
!         IF (CHEMSHIFT2.EQV. .TRUE.) THEN
!           IF (CHEMSHIFT.EQV. .FALSE.) CHEMSHIFTITER=CHEMSHIFTITER-1
!           IF (CHEMSHIFTITER .EQ. 0) CHEMSHIFT=.TRUE.
!         ENDIF
!fh301>}}}
         IF (DEBUG) THEN
!fh301>{{{
           IF (CHEMSHIFT) THEN
             WRITE(MYUNIT,'(A,F20.10,F20.10,G20.10,A,I6,A,F13.10)') ' Energy, Camshiftenergy and RMS force=',ENERGY,
     1       ENERGYCAMSHIFT,RMS,' after ',ITDONE,' LBFGS steps, step:',STP*SLENGTH
           ELSE
!fh301>}}}
             WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1             ' LBFGS steps, step:',STP*SLENGTH
           ENDIF
!fh301>{{{
         ENDIF
!fh301>}}}

!
!  Step finished so can reset OLDQ to new XINT, OLDCART to new LCART,
!  as well as the Cartesian and internal gradients.
!
         IF (CHRMMT.AND.INTMINT) THEN
            OLDGINT(1:N)=GINT(1:N)
            OLDCART(1:3*NATOMS)=LCART(1:3*NATOMS)
!
!  Need to remake XINT because step was only projected in Cartesians?
!  Actually, just setting OLDQ=XINT without this correction seems to
!  be OK. Due to numerical imprecision, it might still be possible
!  for X and XINT to get out of register. Perhaps this doesn't matter
!  because the energy and gradient are always calculated in Cartesians.
!
!           IF (BFGSTST) CALL TRANSDELTA(DELTACART,DELTAQ,LCART,N,3*NATOMS,NNZ,KD)
!           OLDQ(1:N)=OLDQ(1:N)+DELTAQ(1:N)
            OLDQ(1:N)=XINT(1:N)
!         ELSEIF (UNRST) THEN
!!           TEST1(1:N)=X(1:N)
!            CALL geom_to_var(N,X(1:N)) ! testing!!! - to put X back into register with the common block internals (and g)
!!           CALL geom_to_var(N,TEST1(1:N))
!!           do j1=1,N
!!           if (abs((TEST1(j1)-x(j1))/x(j1))*100.0d0.gt.1.0D-6) print *,'hello coords ',J1
!!           enddo
         ENDIF
!
!  Try to take an extra step using the two previous geometries.
! 
!          GOTO 112
!          IF (MOD(ITDONE,3).EQ.0) THEN
!             NGUESS=0
! 111         CONTINUE
!             DO J1=1,NATOMS
!                X1=OLDX(3*(J1-1)+1)-OLDOLDX(3*(J1-1)+1)
!                Y1=OLDX(3*(J1-1)+2)-OLDOLDX(3*(J1-1)+2)
!                Z1=OLDX(3*(J1-1)+3)-OLDOLDX(3*(J1-1)+3)
!                X2=XCOORDS(3*(J1-1)+1)-OLDX(3*(J1-1)+1)
!                Y2=XCOORDS(3*(J1-1)+2)-OLDX(3*(J1-1)+2)
!                Z2=XCOORDS(3*(J1-1)+3)-OLDX(3*(J1-1)+3)
!                VGUESS(1)=(x2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
!      -  ((x2*(y1*y2 + z1*z2) - x1*(y2**2 + z2**2))*
!      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
!      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
!                VGUESS(2)=(y2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
!      -  ((-(x2**2*y1) + x1*x2*y2 + z2*(y2*z1 - y1*z2))*
!      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
!      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
!                VGUESS(3)=(z2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
!      -  ((-(x2**2*z1) + x1*x2*z2 + y2*(-(y2*z1) + y1*z2))*
!      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
!      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
!                D1=SQRT(VGUESS(1)**2+VGUESS(2)**2+VGUESS(3)**2)
!                IF (D1.LT.0.1) THEN
!                   TRY(3*(J1-1)+1)=XCOORDS(3*(J1-1)+1)+VGUESS(1)*1.0D0
!                   TRY(3*(J1-1)+2)=XCOORDS(3*(J1-1)+2)+VGUESS(2)*1.0D0
!                   TRY(3*(J1-1)+3)=XCOORDS(3*(J1-1)+3)+VGUESS(3)*1.0D0
!                ENDIF
!             ENDDO
!             CALL POTENTIAL(TRY,GNEW,EGUESS,.FALSE.,.FALSE.)
!             WRITE(*,'(A,3G20.10)') 'ENEW,EGUESS,change=',ENEW,EGUESS,EGUESS-ENEW
!             IF (EGUESS-ENEW.LT.0.0D0) THEN
!                NGUESS=NGUESS+1
!                ENEW=EGUESS
!                DO J1=1,N
!                   OLDOLDX(J1)=OLDX(J1)
!                   OLDX(J1)=XCOORDS(J1)
!                   XCOORDS(J1)=TRY(J1)
!                ENDDO
!                IF (NGUESS.LT.6) GOTO 111
!             ENDIF
!          ENDIF
!          DO J1=1,N
!             OLDOLDX(J1)=OLDX(J1)
!             OLDX(J1)=XCOORDS(J1)
!          ENDDO
! 
! 112      CONTINUE
!
!  May want to prevent the PE from falling too much if we are trying to visit all the
!  PE bins. Halve the step size until the energy decrease is in range.
!
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
!
!  Energy decreased too much - try again with a smaller step size
!
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
!              XINT=OLDQ ! should be the same as subtracting the step
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
!
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point

!              DO J1=1,N
!                 XCOORDS(J1)=XCOORDS(J1)-STP*W(ISPT+POINT*N+J1)
!              ENDDO

            ENDIF
            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
            ENDIF
            GOTO 30
         ENDIF
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.5*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.5D0
            ENDDO
         ELSE
!
! Resetting to XSAVE and adding half the step should be the same as subtracting 
! half the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.5*STP*W(ISPT+POINT*N+J1)
            ENDDO 
!           DO J1=1,N
!              XCOORDS(J1)=XCOORDS(J1)-0.5*STP*W(ISPT+POINT*N+J1)
!           ENDDO 
!
! For Thomson try projection for the geometry after the step.
!        
            IF (PROJIT) THEN
               IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF

         ENDIF
         STP=STP/2.0D0
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(MYUNIT,'(A,F19.10,A,F16.10,A,F15.8)') 
     1                      ' energy decreased too much from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         
         !FIXIMAGE=.TRUE.
         GOTO 20
      ELSE
!
!  Energy increased - try again with a smaller step size
!
         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
!              XINT=OLDQ ! should be the same as subtracting the step
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
!
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
!              DO J1=1,N
!                 GRAD(J1)=GNEW(J1) ! GRAD contains the gradient at the lowest energy point
!                 XCOORDS(J1)=XCOORDS(J1)-STP*W(ISPT+POINT*N+J1)
!              ENDDO
            ENDIF
            ITER=0   !  try resetting
!            IF (NFAIL.GT.20) THEN
! bs360: smaller NFAIL 
             IF (NFAIL.GT.5) THEN         
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
             ENDIF
            GOTO 30
         ENDIF
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.9*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.1D0
            ENDDO
         ELSE
!
! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
! 0.9 of the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
            ENDDO 

!           DO J1=1,N
!              XCOORDS(J1)=XCOORDS(J1)-0.9*STP*W(ISPT+POINT*N+J1)
!           ENDDO 
!
! For Thomson try projection for the geometry after the step.
!        
            IF (PROJIT) THEN
               IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT .AND. (.NOT. GTHOMSONT) ) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSEIF (THOMSONT.AND.GTHOMSONT.AND.(GTHOMMET .EQ. 5)) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL GTHOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL GTHOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
         ENDIF
         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
!fh301>{{{
!         IF (CHEMSHIFT2.eqv..FALSE.) THEN
!fh301>}}}
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,A,F20.10,A,F20.10)') 
     1                      ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
!fh301>{{{
!         ENDIF
!         IF (CHEMSHIFT2 .EQV. .TRUE.) THEN
!
!           IF (CHEMSHIFT .EQV. .TRUE.) THEN
!             CHEMSHIFTITER=CHEMSHIFTITER+1
!             IF (CHEMSHIFTITER .EQ. 1000) THEN
!               CHEMSHIFTITER=0
!               RETURN
!               CHEMSHIFT=.FALSE.
!               CHEMSHIFTITER=10
!             ENDIF
!           ELSE
!             CHEMSHIFTITER=CHEMSHIFTITER+1
!             IF (CHEMSHIFTITER .GE. 100) THEN
!               CHEMSHIFTITER=0
!               CHEMSHIFT=.TRUE.
!               RETURN
!             ENDIF
!           ENDIF
!         ENDIF
!fh301>}}}
         !FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
!
!     Compute the new step and gradient change
!
30    NPT=POINT*N

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
            W(IYPT+NPT+I)= GINT(I)-W(I)
         ENDDO
      ELSE
         DO J1=1,N
            W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
            W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
!
!  Prevent overflow of dot product for large gradients.
!
            IF (ABS(W(IYPT+NPT+J1)).GT.1.0D10) W(IYPT+NPT+J1)=SIGN(1.0D10,W(IYPT+NPT+J1))
         ENDDO
      ENDIF
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF (DUMPT.AND.DEBUG) THEN
         IF (.NOT.ATOMRIGIDCOORDT) THEN
            RIGIDDUMPCOORDS(1:3*NATOMS) = 0.0D0
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, RIGIDDUMPCOORDS, XCOORDS)
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
            WRITE(DUMPXYZUNIT(NP),'(A,I8,A,G20.10)') 'at step ',ITER,' energy=',ENERGY
            WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     1         ('LA ',RIGIDDUMPCOORDS(3*(J1-1)+1),RIGIDDUMPCOORDS(3*(J1-1)+2),RIGIDDUMPCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     1         ('LB',RIGIDDUMPCOORDS(3*(J1-1)+1),RIGIDDUMPCOORDS(3*(J1-1)+2),RIGIDDUMPCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
         ELSE
            WRITE(DUMPXYZUNIT(NP),'(I4)') NATOMS
            WRITE(DUMPXYZUNIT(NP),'(A,I8,A,G20.10)') 'at step ',ITER,' energy=',ENERGY
            WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT(NP),'(A2,3F20.10)') 
     1         ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
         ENDIF
      ENDIF

      IF (CENT) CALL CENTRE2(XCOORDS)

      SMINKCURRENTP=SMINKCURRENT
      GOTO 10
!     IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)

      RETURN
      END
