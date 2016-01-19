!   Copyright (C) 2003-2010 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE INTLBFGS(QSTART,QFINISH)
USE PORFUNCS
USE COMMONS, ONLY : FREEZENODEST, FREEZETOL, MAXBFGS, CHRMMT, MYUNIT, CQMAX, &
     & INTRMSTOL, INTIMAGE, NREPMAX, NREPULSIVE, MUPDATE, INTDGUESS, &
     & NCONSTRAINT, CONI, CONJ, CONDISTREF, INTCONMAX, &
     & INTCONSTRAINREPCUT, REPCON, INTCONSTRAINTREP, INTREPSEP, NREPI, NREPJ, &
     & CONDISTREFLOCAL, INTCONFRAC, CONACTIVE, REPI, &
     & REPJ, NREPMAX, ATOMACTIVE, NCONSTRAINTON, CONION, CONJON, CONDISTREFLOCALON, CONDISTREFON, &
     & NREPCUT, REPCUT, CHECKCONINT, INTCONSTEPS, INTRELSTEPS, MAXCONE, COLDFUSIONLIMIT, &
     & INTSTEPS1, DUMPINTXYZ, DUMPINTXYZFREQ, DUMPINTEOS, DUMPINTEOSFREQ, &
     & IMSEPMIN, IMSEPMAX, MAXINTIMAGE, INTFREEZET, INTFREEZETOL, FREEZE, &
     & INTFROZEN, CHECKREPINTERVAL, NNREPULSIVE, INTFREEZEMIN, INTIMAGECHECK, &
     & CONCUT, CONCUTLOCAL, NATOMS, DEBUG, STEP, MCSTEPS

IMPLICIT NONE 

DOUBLE PRECISION, INTENT(IN) :: QSTART(3*NATOMS), QFINISH(3*NATOMS)  ! The two end points
INTEGER D, U
DOUBLE PRECISION DMAX, DF, DMIN, LOCALSTEP
INTEGER NDECREASE, NFAIL, NMAXINT, NMININT, JMAX, JMIN, INTIMAGESAVE, NOFF, J1, J2, NQDONE
LOGICAL KNOWE, KNOWG, KNOWH, ADDATOM, ADDREP(NATOMS)
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

DOUBLE PRECISION DUMMY, DPRAND
INTEGER POINT,NPT,J3,J4,NIMAGEFREEZE,NACTIVE,NBEST,NEWATOM
INTEGER TURNONORDER(NATOMS),NBACKTRACK,NQCIFREEZE
INTEGER NDUMMY, NLASTGOODE, NSTEPSMAX
INTEGER NTRIES(NATOMS), NITERDONE, EXITSTATUS, DLIST(NATOMS)
DOUBLE PRECISION :: DDOT,STPMIN, ETOTALTMP, RMSTMP, USEFRAC, STIME, FTIME, &
  &                 ETOTAL, LASTGOODE, RMS, STEPTOT, LINTCONSTRAINTTOL, LXYZ(2*3*NATOMS), &
  &                 BESTWORST, WORST
DOUBLE PRECISION, DIMENSION(MUPDATE)     :: RHO1,ALPHA
DOUBLE PRECISION :: EOLD, DMOVED(NATOMS)
LOGICAL SWITCHED
DOUBLE PRECISION, POINTER :: X(:), G(:)
!
! efk: for freezenodes
!
DOUBLE PRECISION :: TESTG, TOTGNORM
INTEGER :: IM
!
! Dimensions involving INTIMAGE
!
DOUBLE PRECISION, ALLOCATABLE :: TRUEEE(:), &
  &              EEETMP(:), MYGTMP(:), EEE(:), STEPIMAGE(:), &
  &              GTMP(:), DIAG(:), STP(:), SEARCHSTEP(:,:), GDIF(:,:), GLAST(:), XSAVE(:)
DOUBLE PRECISION, ALLOCATABLE, TARGET :: XYZ(:), GGG(:), DPTMP(:), D2TMP(:,:)
! saved interpolation
DOUBLE PRECISION, ALLOCATABLE :: BESTXYZ(:), BESTEEE(:)
INTEGER BESTINTIMAGE, NSTEPS, NITERUSE
LOGICAL, ALLOCATABLE :: CHECKG(:), IMGFREEZE(:)

ALLOCATE(TRUEEE(INTIMAGE+2), &
  &      EEETMP(INTIMAGE+2), MYGTMP(3*NATOMS*INTIMAGE), &
  &      GTMP(3*NATOMS*INTIMAGE), &
  &      DIAG(3*NATOMS*INTIMAGE), STP(3*NATOMS*INTIMAGE), SEARCHSTEP(0:MUPDATE,(3*NATOMS)*INTIMAGE), &
  &      GDIF(0:MUPDATE,(3*NATOMS)*INTIMAGE),GLAST((3*NATOMS)*INTIMAGE), XSAVE((3*NATOMS)*INTIMAGE), &
  &      XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), CHECKG((3*NATOMS)*INTIMAGE), IMGFREEZE(INTIMAGE), &
  &      EEE(INTIMAGE+2), STEPIMAGE(INTIMAGE))
ALLOCATE(BESTXYZ((3*NATOMS)*(INTIMAGE+2)),BESTEEE(INTIMAGE+2))

SWITCHED=.FALSE.
INTIMAGESAVE=INTIMAGE
NBACKTRACK=1
CALL MYCPU_TIME(STIME,.FALSE.)
WRITE(MYUNIT,'(A,I6)') ' intlbfgs> Maximum number of steps for constraint potential phase is ',INTSTEPS1
WRITE(MYUNIT,'(A,I6)') ' intlbfgs> Updates: ',MUPDATE
ADDATOM=.FALSE.
NFAIL=0
IMGFREEZE(1:INTIMAGE)=.FALSE.
D=(3*NATOMS)*INTIMAGE
U=MUPDATE
NITERDONE=1
NITERUSE=1
NQDONE=0

IF ( D<=0 ) THEN
   WRITE(MYUNIT,*) 'd is not positive, d=',d
   STOP
ENDIF
IF ( U<=0 ) THEN
   WRITE(MYUNIT,*) 'u is not positive, u=',u
   STOP
ENDIF
IF (INTSTEPS1 < 0) THEN
   WRITE(MYUNIT,'(1x,a)') 'Maximal number of iterations is less than zero! Stop.'
   STOP
ENDIF
!
! XYZ, GGG, EEE include the end point images
! X, G do not.
!
IF (.NOT.ALLOCATED(CONI)) THEN 
   ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX))
   ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
ENDIF
X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
!
! Initialise XYZ
!
XYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
XYZ((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=QFINISH(1:(3*NATOMS))
DO J1=1,INTIMAGE+2
   XYZ((J1-1)*(3*NATOMS)+1:J1*(3*NATOMS))=((INTIMAGE+2-J1)*QSTART(1:(3*NATOMS))+(J1-1)*QFINISH(1:(3*NATOMS)))/(INTIMAGE+1)
ENDDO
      WRITE(MYUNIT,'(A)') 'intlbfgs> here Z'
      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(398-1)+1:3*(398-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(398-1)+1:(INTIMAGE+1)*3*NATOMS+3*(398-1)+3)
      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(400-1)+1:3*(400-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(400-1)+1:(INTIMAGE+1)*3*NATOMS+3*(400-1)+3)
      WRITE(MYUNIT,'(6G20.10)') QSTART(1:6)
      WRITE(MYUNIT,'(6G20.10)') QFINISH(1:6)

NQCIFREEZE=0
IF (FREEZE) THEN
   WRITE(MYUNIT,'(A)') ' intlbfgs> ERROR *** QCI has not been coded for frozen atoms yet'
   STOP     
ENDIF
IF (ALLOCATED(INTFROZEN)) DEALLOCATE(INTFROZEN)
ALLOCATE(INTFROZEN(NATOMS))
INTFROZEN(1:NATOMS)=.FALSE.
DLIST(1:NATOMS)=-1
DMOVED(1:NATOMS)=1.0D100
IF (INTFREEZET) THEN
   DUMMY=INTFREEZETOL**2
   WRITE(MYUNIT,'(A,6G20.10)') ' intlbfgs> INTFREEZETOL,DUMMY=',INTFREEZETOL,DUMMY
   DO J1=1,NATOMS
      DF=(XYZ(3*(J1-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &     +(XYZ(3*(J1-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &     +(XYZ(3*(J1-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2
      IF (DF.LT.DUMMY) THEN
         NQCIFREEZE=NQCIFREEZE+1
         INTFROZEN(J1)=.TRUE.
         IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,F12.6,A,I6)') &
  &            ' intlbfgs> atom ',J1,' moves less than threshold: dist^2=',DF,' total=',NQCIFREEZE
      ENDIF
      sortd: DO J2=1,J1
         IF (DF.LT.DMOVED(J2)) THEN
            DO J3=J1,J2+1,-1
               DMOVED(J3)=DMOVED(J3-1)
               DLIST(J3)=DLIST(J3-1)
            ENDDO
            DMOVED(J2)=DF
            DLIST(J2)=J1
            EXIT sortd
         ENDIF
      ENDDO sortd
   ENDDO
   WRITE(MYUNIT,'(A,I6,A,F12.6,A,I6)') ' intlbfgs> Total number of atoms moving less than threshold=',NQCIFREEZE
ENDIF

      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(398-1)+1:3*(398-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(398-1)+1:(INTIMAGE+1)*3*NATOMS+3*(398-1)+3)
      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(400-1)+1:3*(400-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(400-1)+1:(INTIMAGE+1)*3*NATOMS+3*(400-1)+3)

IF (NATOMS-NQCIFREEZE.LT.INTFREEZEMIN) THEN
   DO J1=NATOMS,NATOMS-INTFREEZEMIN+1,-1
      INTFROZEN(DLIST(J1))=.FALSE.
   ENDDO
   NQCIFREEZE=NATOMS-INTFREEZEMIN
   WRITE(MYUNIT,'(A,I6,A)') ' intlbfgs> Freezing ',NQCIFREEZE,' atoms'
ENDIF

NLASTGOODE=0
LASTGOODE=1.0D100

!
! Constraints are collected in a list and activated via the CONACTIVE(J1)
! logical array. There will generally be of order NATOMS. However, the
! repulsions will scale as NATOMS**2 and are treated differently. The
! active repulsions are stored sequentially as atoms are added to the
! growing list. This is done even if we have congeom or congeom.dat files
! available. In this case we use the fixed list of possible constraints
! via CHECKPERC, but the list of repulsions and cutoffs is recreated on
! the fly. The fixed lists are used in make_conpot, since this is called
! for pairs of minima with all atoms active to obtain an interpolation
! metric.
!
! Perhaps we should use the fixed list to activate the repulsions below?
! A neighbour list for repulsions is maintained to make the constraint
! potential evaluation scale as order N.
!
IF (NQCIFREEZE.LT.NATOMS) THEN
   LXYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   LXYZ((3*NATOMS)+1:2*(3*NATOMS))=QFINISH(1:(3*NATOMS))
   CALL CHECKPERC(LXYZ,LINTCONSTRAINTTOL,NQCIFREEZE,2)
ELSE
   IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
   NCONSTRAINT=0
   WRITE(MYUNIT,'(A)') ' intlbfgs> All atoms move less than threshold - skip to linear interpolation for end points'
   INTIMAGE=0
   XYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   XYZ((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=QFINISH(1:(3*NATOMS))
   DO J1=1,INTIMAGE+2
      XYZ((J1-1)*(3*NATOMS)+1:J1*(3*NATOMS))=((INTIMAGE+2-J1)*QSTART(1:(3*NATOMS))+(J1-1)*QFINISH(1:(3*NATOMS)))/(INTIMAGE+1)
   ENDDO
   GOTO 678
ENDIF

NACTIVE=0
ATOMACTIVE(1:NATOMS)=.FALSE.
IF (INTFREEZET) THEN
   DO J1=1,NATOMS
      IF (INTFROZEN(J1)) THEN
! 
! linear interpolation 
! 
         DO J2=2,INTIMAGE+1
            XYZ((J2-1)*3*NATOMS+3*(J1-1)+1:(J2-1)*3*NATOMS+3*(J1-1)+3)= &
  &            (INTIMAGE-J2+2)*XYZ(3*(J1-1)+1:3*(J1-1)+3)/(INTIMAGE+1) &
  &           +(J2-1)*XYZ(3*NATOMS*(INTIMAGE+1)+3*(J1-1)+1:3*NATOMS*(INTIMAGE+1)+3*(J1-1)+3)/(INTIMAGE+1)
         ENDDO
         ATOMACTIVE(J1)=.TRUE.
         NACTIVE=NACTIVE+1
         TURNONORDER(NACTIVE)=J1
         NTRIES(J1)=1
      ENDIF
   ENDDO
ENDIF

REPCON=-INTCONSTRAINTREP/INTCONSTRAINREPCUT**6 ! also needed for congrad.f90 potential
IF (ALLOCATED(CONDISTREFLOCAL)) DEALLOCATE(CONDISTREFLOCAL)
IF (ALLOCATED(CONCUTLOCAL)) DEALLOCATE(CONCUTLOCAL)
ALLOCATE(CONDISTREFLOCAL(NCONSTRAINT))
ALLOCATE(CONCUTLOCAL(NCONSTRAINT))
IF (ALLOCATED(CONDISTREFLOCALON)) DEALLOCATE(CONDISTREFLOCALON)
IF (ALLOCATED(CONDISTREFON)) DEALLOCATE(CONDISTREFON)
IF (ALLOCATED(CONION)) DEALLOCATE(CONION)
IF (ALLOCATED(CONJON)) DEALLOCATE(CONJON)
ALLOCATE(CONDISTREFLOCALON(NCONSTRAINT),CONDISTREFON(NCONSTRAINT),CONION(NCONSTRAINT),CONJON(NCONSTRAINT))
CONDISTREFLOCAL(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
CONCUTLOCAL(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)
DUMMY=1.0D100
IF (NCONSTRAINT.EQ.0) THEN
   NACTIVE=NATOMS
   EOLD=ETOTAL
   SWITCHED=.TRUE.
   USEFRAC=1.0D0
   NREPULSIVE=0
   NNREPULSIVE=0
   GLAST(1:D)=G(1:D)
   XSAVE(1:D)=X(1:D)
   GOTO 567
ENDIF
DO J1=1,NCONSTRAINT
   DF=SQRT((XYZ(3*(CONI(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+3))**2)&
  &  +SQRT((XYZ(3*(CONJ(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+3))**2)
   IF (J1.EQ.3505) THEN
      WRITE(MYUNIT,'(A,3I10)') 'intlbfgs> J1,CONI(J1),CONJ(J1)=',J1,CONI(J1),CONJ(J1)
      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(CONI(J1)-1)+1:3*(CONI(J1)-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+1:(INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+3)
      WRITE(MYUNIT,'(6G20.10)') XYZ(3*(CONJ(J1)-1)+1:3*(CONJ(J1)-1)+3), &
  &                             XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+1:(INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+3)
   ENDIF
   IF (DF.LT.DUMMY) THEN
      NBEST=J1
      DUMMY=DF
   ENDIF
ENDDO
IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6,A,F15.5)') ' intlbfgs> Smallest overall motion for constraint ',NBEST,' atoms ', &
  &                           CONI(NBEST),CONJ(NBEST),' distance=',DUMMY

TURNONORDER(1:NATOMS)=0
NTRIES(1:NATOMS)=1
IF (ALLOCATED(CONACTIVE)) DEALLOCATE(CONACTIVE)
ALLOCATE(CONACTIVE(NCONSTRAINT))
CONACTIVE(1:NCONSTRAINT)=.FALSE.
CONACTIVE(NBEST)=.TRUE.
ATOMACTIVE(CONI(NBEST))=.TRUE.
ATOMACTIVE(CONJ(NBEST))=.TRUE.
IF (.NOT.INTFROZEN(CONI(NBEST))) THEN
   TURNONORDER(NACTIVE+1)=CONI(NBEST)
   NACTIVE=NACTIVE+1
ENDIF
IF (.NOT.INTFROZEN(CONJ(NBEST))) THEN
   TURNONORDER(NACTIVE+2)=CONJ(NBEST)
   NACTIVE=NACTIVE+1
ENDIF
NTRIES(CONI(NBEST))=1
NTRIES(CONJ(NBEST))=1
NREPULSIVE=0
NCONSTRAINTON=1
CONDISTREFLOCALON(1)=CONDISTREFLOCAL(NBEST)
CONDISTREFON(1)=CONDISTREF(NBEST)
CONION(1)=CONI(NBEST)
CONJON(1)=CONJ(NBEST)
IF (DEBUG) WRITE(MYUNIT,'(A,I6)') ' intlbfgs> Number of active atoms is now ',NACTIVE
!
! If INTFREEZET is true we need to add constraints and replusions to the frozen atoms.
!
IF (INTFREEZET) THEN
DO J1=1,NCONSTRAINT
   IF (CONACTIVE(J1)) CYCLE
   IF ((CONI(J1).EQ.CONI(NBEST)).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.CONI(NBEST)).AND.(ATOMACTIVE(CONI(J1)))) THEN
      CONACTIVE(J1)=.TRUE.
      IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
   ENDIF
   IF ((CONI(J1).EQ.CONJ(NBEST)).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.CONJ(NBEST)).AND.(ATOMACTIVE(CONI(J1)))) THEN
      CONACTIVE(J1)=.TRUE.
      IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
   ENDIF
ENDDO

DO J1=1,NATOMS
   IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
   IF (ABS(J1-CONI(NBEST)).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
   IF (INTFROZEN(J1).AND.INTFROZEN(CONI(NBEST))) CYCLE
   DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
      IF (.NOT.CONACTIVE(J2)) CYCLE ! identify active constraints
      IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.CONI(NBEST))).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.CONI(NBEST)))) GOTO 545
   ENDDO
   DMIN=1.0D100
   DMAX=-1.0D0
   DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
      DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
      IF (DF.GT.DMAX) DMAX=DF
      IF (DF.LT.DMIN) DMIN=DF
   ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
   DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
   NREPULSIVE=NREPULSIVE+1
   IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
   REPI(NREPULSIVE)=J1
   REPJ(NREPULSIVE)=CONI(NBEST)
   REPCUT(NREPULSIVE)=DMIN
!  IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',CONI(NBEST),' with atom ',J1, &
! &                                          ' cutoff=',DMIN
545 CONTINUE
ENDDO

DO J1=1,NATOMS
   IF (ABS(J1-CONJ(NBEST)).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
   IF (INTFROZEN(J1).AND.INTFROZEN(CONJ(NBEST))) CYCLE
   DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
      IF (.NOT.CONACTIVE(J2)) CYCLE ! identify active constraints
      IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.CONJ(NBEST))).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.CONJ(NBEST)))) GOTO 541
   ENDDO
   DMIN=1.0D100
   DMAX=-1.0D0
   DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
      DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
      IF (DF.GT.DMAX) DMAX=DF
      IF (DF.LT.DMIN) DMIN=DF
   ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
   DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
   NREPULSIVE=NREPULSIVE+1
   IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
   REPI(NREPULSIVE)=J1
   REPJ(NREPULSIVE)=CONJ(NBEST)
   REPCUT(NREPULSIVE)=DMIN
!  IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',CONJ(NBEST),' with atom ',J1, &
! &                                          ' cutoff=',DMIN
541 CONTINUE
ENDDO
ENDIF
CALL MYCPU_TIME(FTIME,.FALSE.)
WRITE(MYUNIT,'(A,F10.1)') ' intlbfgs> constrained potential finished, time=',FTIME-STIME
STIME=FTIME
NSTEPSMAX=INTSTEPS1
!
! Don;t want to redistribute images before even taking a step, so don;t call CHECKSEP.
! Must call CHECKREP to initialise NNREULSIVE, NREPI, NREPJ, etc. SEGV otherwise on second cycle!
!
! To take BH-type steps in the QCI space, jump back here. Leave SWITCHED true.
!
BESTWORST=1.0D100
9876 CONTINUE
CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
IF (CHECKCONINT) THEN
   CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ELSE
   CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ENDIF
EOLD=ETOTAL
GLAST(1:D)=G(1:D)
XSAVE(1:D)=X(1:D)

IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
   WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=', &
  &                       ETOTAL/INTIMAGE,COLDFUSIONLIMIT
   DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT)
   DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
   INTIMAGE=INTIMAGESAVE
   RETURN
ENDIF

! IF (DEBUG) WRITE(*,'(A6,A20,A20,A9,A9)') 'Iter','Energy per image','RMS Force','Step'

567 CONTINUE

DO ! Main do loop with counter NITERDONE, initially set to one
!
!  Add next atom to active set if ADDATOM is true. 
!  Constraints to atoms already in the active set are turned on
!  and short-range repulsions to active atoms that are not distance constrained are turned on.
!  *** OLD Find nearest atom to active set attached by a constraint
!  *** NEW Find atom with most constraints to active set
!  Turn on constraint terms for this atom with all previous members of the active set
!  Add repulsions to non-constrained atoms in this set
!  NTOADD is the number of atoms to add to the active set in each pass. 1 seems best!
!
   IF (ADDATOM.AND.(NACTIVE.LT.NATOMS)) THEN
      CALL DOADDATOM(NCONSTRAINT,NTRIES,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE)
      NLASTGOODE=NITERDONE
      LASTGOODE=ETOTAL
   ENDIF
   GTMP(1:D)=0.0D0
   CALL MAKESTEP(NITERUSE,POINT,DIAG,INTIMAGE,SEARCHSTEP,G,GTMP,STP,GDIF,NPT,D,RHO1,ALPHA)
!
! If the number of images has changed since G was declared then G is not the same
! size as Gtmp and Dot_Product cannot be used.
!
!  IF (Dot_Product(G,Gtmp)/SQRT( Dot_Product(G,G)*Dot_Product(Gtmp,Gtmp) ) > 0.0D0) THEN
!
!  Separate sqrt;s to avoid overflow.
!
   IF (DDOT(D,G,1,GTMP,1)/MAX(1.0D-100,SQRT( DDOT(D,G,1,G,1))*SQRT(DDOT(D,GTMP,1,GTMP,1)) ) > 0.0D0) THEN
        IF (DEBUG) WRITE(MYUNIT,*) 'Search direction has positive projection onto gradient - reversing step'
        GTMP(1:D)=-GTMP(1:D)
        SEARCHSTEP(POINT,1:D)=GTMP(1:D)
   ENDIF
   GTMP(1:D)=G(1:D)

!  We should apply the maximum LBFGS step to each image separately.
!  However, using different scale factors for different images leads to huge
!  discontinuities! Now take the minimum scale factor for all images. DJW 26/11/07

   STPMIN=1.0D0
   DO J2=1,INTIMAGE
      STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,(3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2), &
  &                                    SEARCHSTEP(POINT,(3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2)))
      DUMMY=STEPIMAGE(J2)
      IF (STEPIMAGE(J2) > MAXBFGS) THEN
           STP((3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2) = MAXBFGS/STEPIMAGE(J2)
           STPMIN=MIN(STPMIN,STP((3*NATOMS)*(J2-1)+1))
      ENDIF
!     WRITE(MYUNIT,'(A,I8,3G20.10)') ' image,initial step size,STP,prod=',J2,DUMMY,STP(3*NATOMS*(J2-1)+1), &
! &                                   STEPIMAGE(J2)*STP(3*NATOMS*(J2-1)+1)   
   ENDDO
   STP(1:D)=STPMIN
! EFK: decide whether to freeze some nodes
   IF (FREEZENODEST) THEN
      TOTGNORM=SQRT(DOT_PRODUCT(G(1:(3*NATOMS)*INTIMAGE),G(1:(3*NATOMS)*INTIMAGE))/INTIMAGE)
      NIMAGEFREEZE=0
      DO IM=1,INTIMAGE
         TESTG=SQRT(DOT_PRODUCT(G((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM),G((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM)))
         IMGFREEZE(IM)=.FALSE.
         IF (TOTGNORM.NE.0.0D0) THEN
!           IF (TESTG/TOTGNORM.LT.FREEZETOL) THEN
            IF (TESTG/SQRT(3.0D0*NATOMS).LT.FREEZETOL) THEN
!              IF (DEBUG) PRINT '(A,I6,3G20.10)', ' intlbfgs> Freezing image: ',IM,TESTG,FREEZETOL,TOTGNORM
               IMGFREEZE(IM)=.TRUE.
               STEPIMAGE(IM)=0.0D0
               NIMAGEFREEZE=NIMAGEFREEZE+1
               STP((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM)=0.0D0
            ENDIF
         ENDIF
      ENDDO
      IF (DEBUG) PRINT '(2(A,I6))', ' intlbfgs> Number of frozen images=',NIMAGEFREEZE,' / ',INTIMAGE
   ENDIF
   !  We now have the proposed step - update geometry and calculate new gradient
   NDECREASE=0
20 X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)

!  IF (.NOT.SWITCHED) THEN
   IF (.TRUE.) THEN
!     IF ((RMS.LT.INTRMSTOL*1.0D10).AND.(MOD(NITERDONE,10).EQ.0).AND.(NSTEPSMAX-NITERDONE.GT.100)) &
! &               CALL CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,(3*NATOMS),NATOMS)
      IF (MOD(NITERDONE,INTIMAGECHECK).EQ.0) THEN
864      CONTINUE ! for adding more than one image at a time
         DMAX=0.0D0
         DMIN=HUGE(1.0D0)
         DO J1=1,INTIMAGE+1
            DUMMY=0.0D0
            DO J2=1,3*NATOMS
               IF (ATOMACTIVE((J2-1)/3+1)) THEN
                  DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
               ENDIF
            ENDDO
            DUMMY=SQRT(DUMMY)
            IF (DUMMY.GT.DMAX) THEN
               DMAX=DUMMY
               JMAX=J1
            ENDIF
            IF (DUMMY.LT.DMIN) THEN
               DMIN=DUMMY
               JMIN=J1
            ENDIF
            IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)')' intlbfgs> distance between images ', &
  &                                                  J1,' and ',J1+1,' is ',DUMMY
         ENDDO
         IF ((DMAX.GT.IMSEPMAX).AND.(INTIMAGE.LT.MAXINTIMAGE)) THEN
            WRITE(MYUNIT,'(A,I6,A,I6,A,I6)') ' intlbfgs> Add an image between ',JMAX,' and ',JMAX+1,' INTIMAGE=',INTIMAGE
            NITERUSE=0
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+3)))
            XYZ(1:3*NATOMS*JMAX)=DPTMP(1:3*NATOMS*JMAX)
            XYZ(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1))=(DPTMP(3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX) &
  &                                               + DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1)))/2.0D0
            XYZ(3*NATOMS*(JMAX+1)+1:3*NATOMS*(INTIMAGE+3))=DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+2))
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to MUPDATE over memories and
! 1:(3*NATOMS)*INTIMAGE over only the variable images.
!
            DEALLOCATE(DPTMP)
            ALLOCATE(D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE))
            D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)=SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE+1)))
            DO J1=0,MUPDATE
               IF (JMAX.GT.1) SEARCHSTEP(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) SEARCHSTEP(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               SEARCHSTEP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                             D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)=GDIF(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE+1)))
            DO J1=0,MUPDATE
               IF (JMAX.GT.1) GDIF(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) GDIF(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               GDIF(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                       D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            GDIF(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DEALLOCATE(D2TMP)

            DEALLOCATE(TRUEEE,EEETMP,MYGTMP,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(TRUEEE(INTIMAGE+3), &
  &                  EEETMP(INTIMAGE+3), MYGTMP(3*NATOMS*(INTIMAGE+1)), &
  &                  GTMP(3*NATOMS*(INTIMAGE+1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE+1)), STP(3*NATOMS*(INTIMAGE+1)), &
  &                  GLAST((3*NATOMS)*(INTIMAGE+1)), &
  &                  XSAVE((3*NATOMS)*(INTIMAGE+1)), CHECKG((3*NATOMS)*(INTIMAGE+1)), IMGFREEZE(INTIMAGE+1), &
  &                  EEE(INTIMAGE+3), STEPIMAGE(INTIMAGE+1), GGG(3*NATOMS*(INTIMAGE+3)))
            GGG(1:3*NATOMS*(INTIMAGE+3))=0.0D0
            TRUEEE(1:INTIMAGE+3)=0.0D0
            EEETMP(1:INTIMAGE+3)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            GLAST(1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
            XSAVE(1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
            CHECKG(1:(3*NATOMS)*(INTIMAGE+1))=.FALSE.
            IMGFREEZE(1:INTIMAGE+1)=.FALSE.
            EEE(1:INTIMAGE+3)=0.0D0
            STEPIMAGE(1:INTIMAGE+1)=0.0D0

            X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+2))
            G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+2))
            INTIMAGE=INTIMAGE+1
            D=(3*NATOMS)*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
!           GOTO 864
         ELSEIF ((DMIN.LT.IMSEPMIN).AND.(INTIMAGE.GT.1)) THEN
            IF (JMIN.EQ.1) JMIN=2
            WRITE(MYUNIT,'(A,I6,A,I6)') ' intlbfgs> Remove image ',JMIN
            NITERUSE=0
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+1)))
            XYZ(1:3*NATOMS*(JMIN-1))=DPTMP(1:3*NATOMS*(JMIN-1))
            XYZ(3*NATOMS*(JMIN-1)+1:3*NATOMS*(INTIMAGE+1))=DPTMP(3*NATOMS*JMIN+1:3*NATOMS*(INTIMAGE+2))

            DEALLOCATE(DPTMP)
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to MUPDATE over memories and
! 1:(3*NATOMS)*INTIMAGE over only the variable images.
!
            ALLOCATE(D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE))
            D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)=SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE-1)))
            DO J1=0,MUPDATE
               SEARCHSTEP(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               SEARCHSTEP(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SEARCHSTEP(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            D2TMP(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)=GDIF(0:MUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE-1)))
            DO J1=0,MUPDATE
               GDIF(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               GDIF(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            GDIF(0:MUPDATE,1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DEALLOCATE(D2TMP)

            DEALLOCATE(TRUEEE,EEETMP,MYGTMP,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(TRUEEE(INTIMAGE+1),&
  &                  EEETMP(INTIMAGE+1), MYGTMP(3*NATOMS*(INTIMAGE-1)), &
  &                  GTMP(3*NATOMS*(INTIMAGE-1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE-1)), STP(3*NATOMS*(INTIMAGE-1)), &
  &                  GLAST((3*NATOMS)*(INTIMAGE-1)), &
  &                  XSAVE((3*NATOMS)*(INTIMAGE-1)), CHECKG((3*NATOMS)*(INTIMAGE-1)), IMGFREEZE(INTIMAGE-1), &
  &                  EEE(INTIMAGE+1), STEPIMAGE(INTIMAGE-1), GGG(3*NATOMS*(INTIMAGE+1)))
            GGG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            TRUEEE(1:INTIMAGE+1)=0.0D0
            EEETMP(1:INTIMAGE+1)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            GLAST(1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
            XSAVE(1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
            CHECKG(1:(3*NATOMS)*(INTIMAGE-1))=.FALSE.
            IMGFREEZE(1:INTIMAGE-1)=.FALSE.
            EEE(1:INTIMAGE+1)=0.0D0
            STEPIMAGE(1:INTIMAGE-1)=0.0D0

            X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE))
            G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE))
            INTIMAGE=INTIMAGE-1
            D=(3*NATOMS)*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
!           GOTO 864
         ENDIF
      ENDIF
   ENDIF
!
! End of add/subtract images block.
!
   IF (.NOT.SWITCHED) THEN
      IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      IF ((ETOTAL-EOLD.LT.1.0D100).OR.ADDATOM) THEN ! MAXERISE effectively set to 1.0D100 here
         EOLD=ETOTAL
         GLAST(1:D)=G(1:D)
         XSAVE(1:D)=X(1:D)
      ELSE
         NDECREASE=NDECREASE+1
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(*,'(A,I6)') ' intlbfgs> WARNING *** in lbfgs cannot find a lower energy, NFAIL=',NFAIL
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
         ELSE
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
            STP(1:D)=STP(1:D)/10.0D0
            WRITE(*,'(A,G25.15,A,G25.15,A)') ' intlbfgs> energy increased from ',EOLD,' to ',ETOTAL, &
     &          ' decreasing step size'
            GOTO 20
         ENDIF
      ENDIF
      ADDATOM=.FALSE.
   ELSE ! combine constraint and true potentials
      IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
      ETOTALTMP=0.0D0
      IF (INTCONFRAC.NE.0.0D0) THEN
         DO J4=2,INTIMAGE+1
            IF (CHRMMT) CALL UPDATENBONDS(XYZ((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4))
            CALL POTENTIAL(XYZ((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4),GGG((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4),EEE(J4), &
  &                                    .TRUE.,.FALSE.)
            ETOTALTMP=ETOTALTMP+EEE(J4)
         ENDDO
      ENDIF
      EEETMP(1:INTIMAGE+2)=EEE(1:INTIMAGE+2)
      MYGTMP(1:D)=G(1:D)
      IF (USEFRAC.LT.1.0D0) THEN
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ELSE
         ETOTAL=0.0D0
         G(1:D)=0.0D0
      ENDIF
      ETOTAL=USEFRAC*ETOTALTMP+(1.0D0-USEFRAC)*ETOTAL
      G(1:D)=USEFRAC*MYGTMP(1:D)+(1.0D0-USEFRAC)*G(1:D)
      RMS=SUM(G(1:D)**2)
      RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
      EEE(1:INTIMAGE+2)=USEFRAC*EEETMP(1:INTIMAGE+2)+(1.0D0-USEFRAC)*EEE(1:INTIMAGE+2)
      WORST=-1.0D100
      DO J4=2,INTIMAGE+1
         IF (EEE(J4).GT.WORST) WORST=EEE(J4)
      ENDDO
      WRITE(MYUNIT,'(A,G20.10,A,I8)') 'intlbfgs> Highest QCI image energy=',WORST,' images=',INTIMAGE
   ENDIF
   IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
      WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ETOTAL/INTIMAGE,COLDFUSIONLIMIT
      DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT)
      DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &              DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
      INTIMAGE=INTIMAGESAVE
      RETURN
   ENDIF

   STEPTOT = SUM(STEPIMAGE)/INTIMAGE

   IF (DEBUG) THEN
      WRITE(MYUNIT,'(A,I6,2G20.10,G20.10,I8)') ' intlbfgs> steps: ',NITERDONE,ETOTAL/INTIMAGE,RMS,STEPTOT,NACTIVE
      CALL FLUSH(6)
   ENDIF

   IF (.NOT.SWITCHED) THEN
!     IF ((NITERDONE-NLASTGOODE.GT.INTRELSTEPS).AND.((ETOTAL.GT.LASTGOODE).OR.(ETOTAL/INTIMAGE.GT.MAXCONE*1.0D8))) THEN
      IF (.FALSE.) THEN ! no backtracking
         WRITE(MYUNIT,'(2(A,I6))') ' intlbfgs> Backtracking ',NBACKTRACK,' steps, current active atoms=',NACTIVE
         NTRIES(NEWATOM)=NTRIES(NEWATOM)+1
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
!
! Backtrack by removing the last NBACKTRACK atoms along with their active constraints and
! repulsions.
!
         NOFF=0
         DO J1=1,NBACKTRACK
            NDUMMY=TURNONORDER(NACTIVE-J1+1)
            IF (INTFROZEN(NDUMMY)) THEN
               IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Not turning off frozen active atom ',NDUMMY
               CYCLE
            ENDIF
            IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning off active atom ',NDUMMY
            DO J2=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J2)) CYCLE 
               IF ((CONI(J2).EQ.NDUMMY).OR.(CONJ(J2).EQ.NDUMMY)) THEN
                  CONACTIVE(J2)=.FALSE.
                  IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning off constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
            ATOMACTIVE(NDUMMY)=.FALSE.
            NOFF=NOFF+1
         ENDDO
         NACTIVE=NACTIVE-NOFF
         NDUMMY=1
         NREPULSIVE=0
         DO J1=1,NATOMS
! 
! Make a list of repelling atoms here and then use it
! CONI(J2) is always less than CONJ(J2) so we only need to
! cycle over a given range of constraints and continue from
! where we left off for the next atom j1
!  
            ADDREP(1:J1+INTREPSEP)=.FALSE.
            ADDREP(J1+INTREPSEP+1:NATOMS)=.TRUE. ! no repulsion for atoms too close in sequence
            IF (INTFROZEN(J1)) THEN
               DO J2=J1+INTREPSEP+1,NATOMS
                  IF (INTFROZEN(J2)) ADDREP(J2)=.FALSE.
               ENDDO
            ENDIF
            addloop: DO J2=NDUMMY,NCONSTRAINT
               IF (CONI(J2).EQ.J1) THEN
                  ADDREP(CONJ(J2))=.FALSE.
               ELSE
                  NDUMMY=J2 ! for next atom
                  EXIT addloop
               ENDIF
            ENDDO addloop
            rep2: DO J2=J1+INTREPSEP+1,NATOMS

               IF (.NOT.ADDREP(J2)) CYCLE

               DMIN=1.0D100
               DO J3=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
                  DF=SQRT((XYZ((J3-1)*3*NATOMS+3*(J2-1)+1)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+1))**2+ &
    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+2)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+2))**2+ &
    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+3)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+3))**2)
                  IF (DF.LT.DMIN) DMIN=DF
               ENDDO

               NREPULSIVE=NREPULSIVE+1
               IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
               REPI(NREPULSIVE)=J1
               REPJ(NREPULSIVE)=J2
! 
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
               REPCUT(NREPULSIVE)=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
            ENDDO rep2
         ENDDO


         NBACKTRACK=MAX(MIN(MIN(1.0D0*(NBACKTRACK+1),1.0D0*50),1.0D0*(NACTIVE-2-NQCIFREEZE)),1.0D0)
!        IF (DEBUG) WRITE(MYUNIT,'(A,I6)') ' intlbfgs> Number of atoms to backtrack is now ',NBACKTRACK
         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            WRITE(MYUNIT,'(A,I6)') ' intlbfgs> ERROR *** inconsistency in number of active atoms. ',NDUMMY,' should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) WRITE(MYUNIT,'(A,I6)') ' active atom ',J1
            ENDDO
            STOP
         ENDIF
         ADDATOM=.TRUE.

         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ENDIF
      LASTGOODE=ETOTAL
   ENDIF
   EXITSTATUS=0
   INTDGUESS=DIAG(1) ! should be ok for subsequent runs of the same system DJW
   IF ((.NOT.SWITCHED).AND.(RMS<=INTRMSTOL).AND.NITERDONE>1) EXITSTATUS=1 
   IF (SWITCHED.AND.(RMS<=CQMAX).AND.NITERDONE>1) EXITSTATUS=1 
   IF (NITERDONE==NSTEPSMAX) EXITSTATUS=2

   IF (EXITSTATUS > 0) THEN  
      IF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.1)) THEN ! add active atom or restart with true potential on
         IF (ETOTAL/INTIMAGE.GT.MAXCONE) GOTO 777
         IF (NACTIVE.LT.NATOMS) THEN 
            ADDATOM=.TRUE.
            GOTO 777
         ENDIF
         CALL MYCPU_TIME(FTIME,.FALSE.)
         WRITE(MYUNIT,'(A,I6,A,F12.6,A,I6,A,F10.1)') ' intlbfgs> switch on true potential at step ',NITERDONE, &
  &                                     ' fraction=',INTCONFRAC,' images=',INTIMAGE,' time=',FTIME-STIME
         IF (DEBUG) CALL RWG(NITERDONE,INTIMAGE,XYZ)
         IF (DEBUG) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE,MYUNIT)
         WRITE(MYUNIT,'(A,I6,A,F15.6)') ' intlbfgs> Allowing ',INTCONSTEPS,' further optimization steps'
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) THEN
               WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> ERROR *** number of active atoms=',NACTIVE,' but atom ',J1,' is not active'
            ENDIF
         ENDDO
         NSTEPSMAX=NITERDONE+INTCONSTEPS
         SWITCHED=.TRUE.
         RMS=INTRMSTOL*10.0D0 ! to prevent premature convergence
         G(1:(3*NATOMS)*INTIMAGE)=INTRMSTOL*10.0D0
         USEFRAC=INTCONFRAC
         GOTO 777
      ELSEIF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.2)) THEN 
         WRITE(MYUNIT,'(A,I6)') ' intlbfgs> ERROR *** number of active atoms at final step=',NACTIVE
         CALL FLUSH(6)
         RETURN
      ELSEIF (DEBUG) THEN
         WRITE(MYUNIT,'(A,I6,A,I6)') 'intlbfgs> energies for images:'
         WRITE(MYUNIT,'(I6,F20.10)') (J2,EEE(J2),J2=1,INTIMAGE+2)
      ENDIF
      EXIT
   ENDIF
   777 CONTINUE
!
! Compute the new step and gradient change
!
   NPT=POINT*D
   SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
   GDIF(POINT,:)=G-GTMP
   
   POINT=POINT+1; IF (POINT==MUPDATE) POINT=0

   IF (DUMPINTXYZ.AND.MOD(NITERDONE,DUMPINTXYZFREQ)==0) CALL RWG(NITERDONE,INTIMAGE,XYZ)
   IF (DUMPINTEOS.AND.MOD(NITERDONE,DUMPINTEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE,MYUNIT)

   NITERDONE=NITERDONE+1
   NITERUSE=NITERUSE+1

   IF (NITERDONE.GT.NSTEPSMAX) EXIT
   IF (NACTIVE.EQ.NATOMS) THEN
      IF (.NOT.SWITCHED) THEN
         CALL MYCPU_TIME(FTIME,.FALSE.)
         WRITE(MYUNIT,'(A,I6,A,F12.6,A,I6,A,F10.1)') ' intlbfgs> switch on true potential at step ',NITERDONE, &
  &                                     ' fraction=',INTCONFRAC,' images=',INTIMAGE,' time=',FTIME-STIME
         WRITE(MYUNIT,'(A,I6,A,F15.6)') ' intlbfgs> Allowing ',INTCONSTEPS,' further optimization steps'
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) THEN
               WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> ERROR *** number of active atoms=',NACTIVE,' but atom ',J1,' is not active'
            ENDIF
         ENDDO
         NSTEPSMAX=NITERDONE+INTCONSTEPS
         SWITCHED=.TRUE.
         IF (FREEZENODEST) THEN
            IMGFREEZE(1:INTIMAGE)=.FALSE.
         ENDIF
         RMS=INTRMSTOL*10.0D0 ! to prevent premature convergence
         USEFRAC=INTCONFRAC
      ENDIF
   ENDIF

ENDDO ! end of main do loop over counter NITERDONE

      CALL FLUSH(6)

IF (.NOT.SWITCHED) THEN 
   WRITE(MYUNIT,'(A,I6,A)') ' intlbfgs> ERROR *** number of active atoms at final step=',NACTIVE,' no potential switch'
   STOP
ENDIF
IF (EXITSTATUS.EQ.1) THEN
   WRITE(MYUNIT,'(A,I6,A,G20.10,A,G15.8,A,I4)') ' intlbfgs> Converged after ',NITERDONE,' steps, energy/image=',ETOTAL/INTIMAGE, &
  &                               ' RMS=',RMS,' images=',INTIMAGE
ELSEIF (EXITSTATUS.EQ.2) THEN
   WRITE(MYUNIT,'(A,I6,A,G20.10,A,G15.8,A,I4)') ' intlbfgs> After ',NITERDONE,' steps, energy/image=',ETOTAL/INTIMAGE, &
  &                               ' RMS=',RMS,' images=',INTIMAGE
ENDIF
!
! Linear interpolation for constraint potential and real potential separately.
! Constraint potential need not be flat if we have done some steps with both
! potentials turned on.
!
678 CONTINUE

CALL RWG(NITERDONE,INTIMAGE,XYZ)
CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE,MYUNIT)
NQDONE=NQDONE+1

WRITE(MYUNIT,'(A,G20.10)') 'intlbfgs> WORST=',WORST
WRITE(MYUNIT,'(A,2I8)') 'intlbfgs> NQDONE,MCSTEPS=',NQDONE,MCSTEPS(1)
IF (WORST.EQ.0.0D0) GOTO 8765
IF (NQDONE.EQ.MCSTEPS(1)) GOTO 8765

!
! Accept/reject this QCI set until BH steps exceeded, or worst energy is zero.
! Reset everything necessary and go back to 9876 CONTINUE - start of QCI optimisation.
!
IF (WORST.GT.BESTWORST) THEN
   WRITE(MYUNIT,'(A)') 'intlbfgs> rejecting step - resetting'
   DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &         DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
   NULLIFY(X,G)
   INTIMAGE=BESTINTIMAGE
   D=(3*NATOMS)*INTIMAGE

   ALLOCATE(TRUEEE(INTIMAGE+2), &
  &      EEETMP(INTIMAGE+2), MYGTMP(3*NATOMS*INTIMAGE), &
  &      GTMP(3*NATOMS*INTIMAGE), &
  &      DIAG(3*NATOMS*INTIMAGE), STP(3*NATOMS*INTIMAGE), SEARCHSTEP(0:MUPDATE,(3*NATOMS)*INTIMAGE), &
  &      GDIF(0:MUPDATE,(3*NATOMS)*INTIMAGE),GLAST((3*NATOMS)*INTIMAGE), XSAVE((3*NATOMS)*INTIMAGE), &
  &      XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), CHECKG((3*NATOMS)*INTIMAGE), IMGFREEZE(INTIMAGE), &
  &      EEE(INTIMAGE+2), STEPIMAGE(INTIMAGE))
   XYZ(1:(3*NATOMS)*(INTIMAGE+2))=BESTXYZ(1:(3*NATOMS)*(INTIMAGE+2))
   EEE(1:INTIMAGE+2)=BESTEEE(1:INTIMAGE+2)
   X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
   G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
   WRITE(MYUNIT,'(A,I8,A,G20.10)') 'intlbfgs> resetting QCI to ',INTIMAGE,' images, highest energy=',BESTWORST
ELSE
   BESTWORST=WORST
   BESTINTIMAGE=INTIMAGE
   IF (ALLOCATED(BESTXYZ)) DEALLOCATE(BESTXYZ,BESTEEE)
   ALLOCATE(BESTXYZ((3*NATOMS)*(INTIMAGE+2)),BESTEEE(INTIMAGE+2))
   BESTXYZ(1:(3*NATOMS)*(INTIMAGE+2))=XYZ(1:(3*NATOMS)*(INTIMAGE+2))
   BESTEEE(1:INTIMAGE+2)=EEE(1:INTIMAGE+2)
   WRITE(MYUNIT,'(A,I8,A,G20.10)') 'intlbfgs> retaining ',INTIMAGE,' QCI images, highest energy=',BESTWORST
ENDIF
!
! Take a step by perturbing all atoms in all images. Perhaps not images with zero energy?
!
WRITE(MYUNIT,'(A,3I8)') 'intlbfgs> NSTEPSMAX,INTCONSTEPS,NITERDONE=',NSTEPSMAX,INTCONSTEPS,NITERDONE
NSTEPSMAX=INTCONSTEPS
NITERDONE=1
LOCALSTEP=STEP(1)
DO J4=2,INTIMAGE+1
   IF (EEE(J4).GT.0.0D0) THEN
      IF (DEBUG) WRITE(MYUNIT,'(A,I8,A,G20.10)') 'intlbfgs> Perturbing image ',J4,' energy=',EEE(J4)
      DO J1=1,NATOMS
         J2=3*NATOMS*(J4-1)+3*J1
         XYZ(J2-2)=XYZ(J2-2)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
         XYZ(J2-1)=XYZ(J2-1)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
         XYZ(J2)=  XYZ(J2)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
      ENDDO
   ELSE
      IF (DEBUG) WRITE(MYUNIT,'(A,I8,A,G20.10)') 'intlbfgs> Not perturbing image ',J4,' energy=',EEE(J4)
   ENDIF
ENDDO
GOTO 9876

8765 CONTINUE ! jump here if all images have zero energy

CALL RWG(0,INTIMAGE,XYZ)
CALL WRITEPROFILE(0,EEE,INTIMAGE,MYUNIT)

DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT)
DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
INTIMAGE=INTIMAGESAVE

STOP

END SUBROUTINE INTLBFGS
!
! Neighbour list for repulsions to reduce cost of constraint potential.
!
SUBROUTINE CHECKREP(INTIMAGE,XYZ,NOPT,NNSTART,NSTART)
USE COMMONS,ONLY : NREPI, NREPJ, NREPCUT, NNREPULSIVE, NREPULSIVE, REPI, REPJ, REPCUT, CHECKREPCUTOFF, DEBUG, MYUNIT 
USE PORFUNCS
IMPLICIT NONE
INTEGER JJ, KK, NI1, NJ1, NI2, NJ2, INTIMAGE, NOPT, NI, NJ, NNSTART, NSTART
DOUBLE PRECISION LDIST, XYZ(NOPT*(INTIMAGE+2)),COMPARE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN
LOGICAL NOINT

NNREPULSIVE=NNSTART
DO JJ=NSTART,NREPULSIVE
   COMPARE=(CHECKREPCUTOFF*REPCUT(JJ))**2
   NI=REPI(JJ)
   NJ=REPJ(JJ)
   DO KK=1,INTIMAGE+2 ! first check for standard distances within threshold
      LDIST=(XYZ((KK-1)*NOPT+3*(NI-1)+1)-XYZ((KK-1)*NOPT+3*(NJ-1)+1))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+2)-XYZ((KK-1)*NOPT+3*(NJ-1)+2))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+3)-XYZ((KK-1)*NOPT+3*(NJ-1)+3))**2
      IF (LDIST.LT.COMPARE) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=NI
         NREPJ(NNREPULSIVE)=NJ
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
         GOTO 246
      ENDIF
   ENDDO 
   COMPARE=CHECKREPCUTOFF*REPCUT(JJ)
   DO KK=2,INTIMAGE+2 ! now check internal minima within threshold
      DMIN=1.0D10
      NI2=NOPT*(KK-2)+3*(NI-1)
      NI1=NOPT*(KK-1)+3*(NI-1)
      NJ2=NOPT*(KK-2)+3*(NJ-1)
      NJ1=NOPT*(KK-1)+3*(NJ-1)
      R1AX=XYZ(NI2+1); R1AY=XYZ(NI2+2); R1AZ=XYZ(NI2+3)
      R1BX=XYZ(NJ2+1); R1BY=XYZ(NJ2+2); R1BZ=XYZ(NJ2+3)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      CALL INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN,NOINT)

      IF (NOINT) CYCLE
      IF (DMIN.LT.COMPARE) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=NI
         NREPJ(NNREPULSIVE)=NJ
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
         GOTO 246
      ENDIF
   ENDDO 
246 CONTINUE
ENDDO
IF (DEBUG) WRITE(MYUNIT,'(A,2I8)') ' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE

END SUBROUTINE CHECKREP

SUBROUTINE RWG(NITER,INTIMAGE,XYZ)
USE PORFUNCS
USE COMMONS,ONLY: STOCKT,STOCKAAT, RBAAT, ZSYM, NATOMS, MYUNIT
IMPLICIT NONE
CHARACTER(LEN=10) :: XYZFILE   = 'int.xyz   '
INTEGER,INTENT(IN) :: NITER
INTEGER :: J1,J2,INTIMAGE
CHARACTER(LEN=80) :: FILENAME,DUMMYS
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2))

FILENAME=XYZFILE

IF (NITER.GT.0) THEN
   WRITE(DUMMYS,'(I8)') NITER
   FILENAME='int.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz' ! so that vmd recognises the file type!
ENDIF
OPEN(UNIT=993,FILE=FILENAME,STATUS='replace')
DO J2=1,INTIMAGE+2
   WRITE(993,'(I4/)') NATOMS
!  WRITE(993,'(A5,1X,3F20.10)') (ZSYM((J1+2)/3),xyz( (j2-1)*(3*NATOMS)+j1),&
   WRITE(993,'(A5,1X,3F20.10)') ('LA ',XYZ( (j2-1)*(3*NATOMS)+J1), &
 &     XYZ((J2-1)*(3*NATOMS)+J1+1), XYZ((J2-1)*(3*NATOMS)+J1+2),J1=1,(3*NATOMS),3)
ENDDO

WRITE(MYUNIT,*) 'rwg> Interpolated image coordinates were saved to xyz file "'//TRIM(FILENAME)//'"'

CLOSE(UNIT=993)
END SUBROUTINE RWG

SUBROUTINE WRITEPROFILE(NITER,EEE,INTIMAGE,MYUNIT)
IMPLICIT NONE 
INTEGER,INTENT(IN) :: NITER, INTIMAGE
INTEGER :: I,UNIT,MYUNIT
DOUBLE PRECISION :: EEE(INTIMAGE+2)
CHARACTER(LEN=20) :: FILENAME

UNIT=992
IF (NITER.GT.0) THEN
   WRITE(FILENAME,'(I8)') NITER
   FILENAME='int.EofS.' // TRIM(ADJUSTL(FILENAME))
ELSE   
   FILENAME='int.EofS'
ENDIF
OPEN(UNIT=UNIT,FILE=FILENAME,STATUS='replace')

WRITE(UNIT=UNIT,FMT='(2g24.13)') EEE(1)
DO I=2,INTIMAGE+1
   WRITE(UNIT=UNIT,FMT='(2G24.13)') EEE(I)
ENDDO
WRITE(UNIT=UNIT,FMT='(2G24.13)') EEE(INTIMAGE+2)

CLOSE(UNIT)
WRITE(MYUNIT,'(A)') ' writeprofile> Interpolated energy profile was saved to file "'//trim(filename)//'"'

END SUBROUTINE WRITEPROFILE

SUBROUTINE DOADDATOM(NCONSTRAINT,NTRIES,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE)
USE COMMONS, ONLY : CONACTIVE, CONI, CONJ, ATOMACTIVE, CONDISTREF, REPI, REPJ, REPCUT, INTREPSEP,  &
  &             INTCONSTRAINREPCUT, NREPULSIVE, NREPMAX, MAXCONUSE, CHECKCONINT, &
  &             FREEZENODEST, NNREPULSIVE, NATOMS, DEBUG, MYUNIT
IMPLICIT NONE
INTEGER INTIMAGE
INTEGER NBEST, NCONTOACTIVE(NATOMS),  NCONSTRAINT, J2, NTRIES(NATOMS), NEWATOM,  CONLIST(NATOMS), N1, N2, N3, &
  &     NTOADD, NADDED, NMININT, NMAXINT, TURNONORDER(NATOMS), NDUMMY, J1, J3, NITERDONE, NCONFORNEWATOM, NACTIVE
DOUBLE PRECISION DUMMY, DUMMY2, DPRAND, RANDOM, CONDIST(NATOMS), DMIN
INTEGER NDFORNEWATOM, BESTPRESERVEDN(NATOMS)
DOUBLE PRECISION BESTPRESERVEDD(NATOMS), BESTCLOSESTD(NATOMS), INVDTOACTIVE(NATOMS)
LOGICAL IMGFREEZE(INTIMAGE)
DOUBLE PRECISION C1, C2, C3, VEC1(3), VEC2(3), VEC3(3), ESAVED, ESAVEC, ESAVE0
INTEGER NCFORNEWATOM, BESTCLOSESTN(NATOMS), NNREPSAVE, NREPSAVE
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), XSAVED(3,INTIMAGE+2), XSAVEC(3,INTIMAGE+2), XSAVE0(3,INTIMAGE+2),FRAC,RAN1, &
  &              RMS,EEE(INTIMAGE+2),GGG((3*NATOMS)*(INTIMAGE+2)),ETOTAL,DS,DF

NTOADD=1
!  NTOADD=NATOMS-2  !!!! DJW
NADDED=0

!
! Save current number of repulsions and number that are active to speed up the
! calls to CHECKREP
!
NNREPSAVE=NNREPULSIVE
NREPSAVE=NREPULSIVE
542   CONTINUE
!     DUMMY=1.0D100
      NBEST=0
      NCONTOACTIVE(1:NATOMS)=0
      INVDTOACTIVE(1:NATOMS)=0.0D0
      DO J2=1,NCONSTRAINT
         IF (CONACTIVE(J2)) CYCLE   ! count new, inactive constraints
         IF (ATOMACTIVE(CONI(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONJ(J2))) THEN
               NCONTOACTIVE(CONJ(J2))=NCONTOACTIVE(CONJ(J2))+1
               INVDTOACTIVE(CONJ(J2))=INVDTOACTIVE(CONJ(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (ATOMACTIVE(CONJ(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONI(J2))) THEN
               NCONTOACTIVE(CONI(J2))=NCONTOACTIVE(CONI(J2))+1
               INVDTOACTIVE(CONI(J2))=INVDTOACTIVE(CONI(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (NCONTOACTIVE(CONI(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONI(J2))
         ENDIF
         IF (NCONTOACTIVE(CONJ(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONJ(J2))
         ENDIF
!        WRITE(MYUNIT,'(A,7I6)') 'J2,NCONTOACTIVEI,NCONTOACTOVEJ,CONI,CONJ,NEWATOM,NBEST=', &
! &                             J2,NCONTOACTIVE(CONI(J2)),NCONTOACTIVE(CONJ(J2)),CONI(J2),CONJ(J2),NEWATOM,NBEST

      ENDDO
!
!  Choose NEWATOM stochastically. Bias towards atoms with the maximum constraints.
!  Use a normalised probability and generate a random number between 0 and 1.
!
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         IF (NCONTOACTIVE(J2).EQ.0) CYCLE
         IF (ATOMACTIVE(J2)) CYCLE
!        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
         DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4 
!        WRITE(MYUNIT,'(A,I6,A,G20.10)') ' intlbfgs> Unnormalised probability for choosing atom ',J2,' is ', &
! &                ((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4
      ENDDO

      RANDOM=DUMMY2*DPRAND()
      DUMMY2=0.0D0
      choosenew: DO J2=1,NATOMS
         IF (NCONTOACTIVE(J2).EQ.0) CYCLE
         IF (ATOMACTIVE(J2)) CYCLE
!        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
         DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4 
         IF (DUMMY2.GE.RANDOM) THEN
            NEWATOM=J2
            IF (DEBUG) WRITE(MYUNIT,'(3(A,I6))') ' intlbfgs> Choosing new active atom ',NEWATOM,' new constraints=', &
  &                                       NCONTOACTIVE(J2),' maximum=',NBEST
            EXIT choosenew
         ENDIF
      ENDDO choosenew
          
      IF (NEWATOM*NBEST.EQ.0) THEN ! sanity check
         WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> ERROR *** new active atom not set'
         STOP
      ELSE
!
!  We need a sorted list of up to 3 active atoms, sorted according to how well the
!  end point distance is preserved, even if they don't satisfy the constraint 
!  condition. We want three atoms to use for a local axis system in the interpolation.
!
!  Try sorting on the shortest average distances in the endpoint structures instead, to avoid
!  problems with distant atoms acidentally having a well-preserved distance.
!
         NDFORNEWATOM=0
         BESTPRESERVEDD(1:NATOMS)=1.0D100
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
            DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
            DUMMY=ABS(DS-DF)
            NDFORNEWATOM=NDFORNEWATOM+1
            DO J2=1,NDFORNEWATOM 
               IF (DUMMY.LT.BESTPRESERVEDD(J2)) THEN
!                 WRITE(MYUNIT,'(A,I6,G12.4,I6,G12.4)') 'J1,DUMMY < J2,BESTPRESERVEDD: ',J1,DUMMY,J2,BESTPRESERVEDD(J2)
                  DO J3=NDFORNEWATOM,J2+1,-1 
!                    WRITE(MYUNIT,'(A,I6,A,I6,A,G12.4)') ' moving diff and list from ',J3-1,' to ',J3, &
!&                                               ' DIFF=',BESTPRESERVEDD(J3-1)
                     BESTPRESERVEDD(J3)=BESTPRESERVEDD(J3-1)
                     BESTPRESERVEDN(J3)=BESTPRESERVEDN(J3-1)
                  ENDDO
                  BESTPRESERVEDD(J2)=DUMMY
!                 WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting BESTPRESERVEDD element ',J2,' to ',DUMMY
                  BESTPRESERVEDN(J2)=J1
!                 WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting BESTPRESERVEDN element ',J2,' to ',J1
                  GOTO 653
               ENDIF
            ENDDO
653         CONTINUE
         ENDDO
         IF (DEBUG) THEN
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' best preserved distances:'
            WRITE(MYUNIT,'(20I6)') BESTPRESERVEDN(1:MIN(10,NDFORNEWATOM))
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> sorted differences:'
            WRITE(MYUNIT,'(10G12.4)') BESTPRESERVEDD(1:MIN(10,NDFORNEWATOM))
         ENDIF
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.

         NCFORNEWATOM=0
         BESTCLOSESTD(1:NATOMS)=1.0D100
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
            DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
            DUMMY=(DS+DF)/2.0D0
            NCFORNEWATOM=NCFORNEWATOM+1
            DO J2=1,NCFORNEWATOM
               IF (DUMMY.LT.BESTCLOSESTD(J2)) THEN
!                 WRITE(MYUNIT,'(A,I6,G12.4,I6,G12.4)') 'J1,DUMMY < J2,BESTCLOSESTD: ',J1,DUMMY,J2,BESTCLOSESTD(J2)
                  DO J3=NCFORNEWATOM,J2+1,-1
!                    WRITE(MYUNIT,'(A,I6,A,I6,A,G12.4)') ' moving diff and list from ',J3-1,' to ',J3, &
!&                                               ' DIFF=',BESTCLOSESTD(J3-1)
                     BESTCLOSESTD(J3)=BESTCLOSESTD(J3-1)
                     BESTCLOSESTN(J3)=BESTCLOSESTN(J3-1)
                  ENDDO
                  BESTCLOSESTD(J2)=DUMMY
!                 WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting BESTCLOSESTD element ',J2,' to ',DUMMY
                  BESTCLOSESTN(J2)=J1
!                 PRINT '(A,I6,A,G12.4)',' setting BESTCLOSESTN element ',J2,' to ',J1
                  GOTO 659
               ENDIF
            ENDDO
659         CONTINUE
         ENDDO
         IF (DEBUG) THEN
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' shortest average distances in endpoints:'
            WRITE(MYUNIT,'(20I6)') BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> sorted differences:'
            WRITE(MYUNIT,'(10G12.4)') BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
         ENDIF
!
!  Maintain a sorted list of active atoms that are constrained to the new atom, sorted
!  according to their distance.
!
         NCONFORNEWATOM=0
         CONDIST(1:NATOMS)=1.0D100
         IF (DEBUG) WRITE(MYUNIT,'(3(A,I6))') ' intlbfgs> New active atom is number ',NEWATOM,' total=',NACTIVE+1, &
 &                        ' steps=',NITERDONE
         DO J1=1,NCONSTRAINT
            IF (CONACTIVE(J1)) CYCLE
            IF ((CONI(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONI(J1)))) THEN  
                 NCONFORNEWATOM=NCONFORNEWATOM+1
!                CONACTIVE(J1)=.TRUE.
!                NITSTART(J1)=NITERDONE
!                NCONSTRAINTON=NCONSTRAINTON+1
! !
! ! The ...ON variables are not actually used in congrad.f90.
! !
!                CONDISTREFLOCALON(NCONSTRAINTON)=CONDISTREFLOCAL(J1)
!                CONDISTREFON(NCONSTRAINTON)=CONDISTREF(J1)
!                CONION(NCONSTRAINTON)=CONI(J1)
!                CONJON(NCONSTRAINTON)=CONJ(J1)
! 
!                IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
               IF (NCONFORNEWATOM.EQ.1) THEN
                  CONDIST(1)=CONDISTREF(J1)
                  IF (CONI(J1).EQ.NEWATOM) CONLIST(1)=CONJ(J1)
                  IF (CONJ(J1).EQ.NEWATOM) CONLIST(1)=CONI(J1)
               ENDIF
               DO J2=1,NCONFORNEWATOM-1
                  IF (CONDISTREF(J1).LT.CONDIST(J2)) THEN
!                    WRITE(MYUNIT,'(A,I6,G12.4,I6,G12.4)') 'J1,CONDISTREF < J2,CONDIST: ',J1,CONDISTREF(J1),J2,CONDIST(J2)
                     DO J3=NCONFORNEWATOM,J2+1,-1
!                       WRITE(MYUNIT,'(A,I6,A,I6,A,G12.4)') ' moving dist and list from ',J3-1,' to ',J3,' CONDIST=',CONDIST(J3-1)
                        CONDIST(J3)=CONDIST(J3-1)
                        CONLIST(J3)=CONLIST(J3-1)
                     ENDDO
                     CONDIST(J2)=CONDISTREF(J1)
!                    WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting condist element ',J2,' to ',CONDISTREF(J1)
                     IF (CONI(J1).EQ.NEWATOM) CONLIST(J2)=CONJ(J1)
                     IF (CONJ(J1).EQ.NEWATOM) CONLIST(J2)=CONI(J1)
!                    WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting conlist element ',J2,' to ',CONLIST(J2)
                     GOTO 654
                  ENDIF
               ENDDO 
               CONDIST(NCONFORNEWATOM)=CONDISTREF(J1)
!              WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting condist element ',NCONFORNEWATOM,' to ',CONDISTREF(J1)
               IF (CONI(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONJ(J1)
               IF (CONJ(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONI(J1)
!              WRITE(MYUNIT,'(A,I6,A,G12.4)') ' setting conlist element ',NCONFORNEWATOM,' to ',CONLIST(NCONFORNEWATOM)
654          CONTINUE
            ENDIF
         ENDDO 
         IF (DEBUG) THEN
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' is constrained to ',NCONFORNEWATOM, &
  &                                       ' other active atoms:'
            WRITE(MYUNIT,'(20I6)') CONLIST(1:NCONFORNEWATOM)
            WRITE(MYUNIT,'(A,I6,A,I6,A)') ' intlbfgs> sorted distances:'
            WRITE(MYUNIT,'(10G12.4)') CONDIST(1:NCONFORNEWATOM)
         ENDIF
         DO J1=1,MIN(MAXCONUSE,NCONFORNEWATOM)
            DO J2=1,NCONSTRAINT
               IF ((CONI(J2).EQ.NEWATOM).AND.(CONJ(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ELSE IF ((CONJ(J2).EQ.NEWATOM).AND.(CONI(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
            IF (ABS(J1-NEWATOM).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
            DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
               IF (.NOT.CONACTIVE(J2)) CYCLE ! identify active constraints 
               IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.NEWATOM)).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.NEWATOM))) GOTO 543
            ENDDO
            DMIN=1.0D100
            DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
               DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
               IF (DF.LT.DMIN) DMIN=DF
            ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
            DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
            NREPULSIVE=NREPULSIVE+1
            IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
            REPI(NREPULSIVE)=J1
            REPJ(NREPULSIVE)=NEWATOM
            REPCUT(NREPULSIVE)=DMIN
!           IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',NEWATOM,' with atom ',J1, &
! &                                                   ' cutoff=',DMIN
543         CONTINUE
         ENDDO
         ATOMACTIVE(NEWATOM)=.TRUE.
         NACTIVE=NACTIVE+1

         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            WRITE(MYUNIT,'(A,I6)') ' intlbfgs> ERROR *** inconsistency in number of active atoms. ',NDUMMY,' should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) WRITE(MYUNIT,'(A,I6)') ' active atom ',J1
            ENDDO
            STOP
         ENDIF

         TURNONORDER(NACTIVE)=NEWATOM
!
! Initial guess for new active atom position. This is crucial for success in INTCONSTRAINT schemes!
!
         ESAVED=1.0D100
         ESAVE0=1.0D100
         ESAVEC=1.0D100
         IF (NCONFORNEWATOM.GE.3) THEN
!
! Move the new atom consistently in the local environment of its three nearest actively constrained atoms.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
            IF (DEBUG) WRITE(MYUNIT,'(A)') ' intlbfgs> initial guess from closest three constrained active atoms'
            VEC1(1:3)=XYZ(3*(CONLIST(2)-1)+1:3*(CONLIST(2)-1)+3)-XYZ(3*(CONLIST(1)-1)+1:3*(CONLIST(1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(CONLIST(3)-1)+1:3*(CONLIST(3)-1)+3)-XYZ(3*(CONLIST(1)-1)+1:3*(CONLIST(1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)+C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVE0=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVE0(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF
         IF (NDFORNEWATOM.GE.3) THEN
!
! Choose three atoms from the BESTPRESERVEDN list at random with bias towards the 
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
            DUMMY=0.0D0
            DO J1=1,NDFORNEWATOM
!              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
!              DUMMY=DUMMY+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
               DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
            ENDDO
            N1=0; N2=0; N3=0
            DO WHILE (N3.EQ.0)
               DUMMY2=0.0D0
               RAN1=DPRAND()*DUMMY
               DO J1=1,NDFORNEWATOM
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
                  DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
                  IF (DUMMY2.GE.RAN1) THEN
                     IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
                     IF (N1.EQ.0) THEN
                        N1=J1
                        EXIT
                     ENDIF
                     IF (N2.EQ.0) THEN
                        N2=J1
                        EXIT
                     ENDIF
                     N3=J1
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF (DEBUG) WRITE(MYUNIT,'(A,3I6,A)') ' intlbfgs> choosing positions ',N1,N2,N3,' in best preserved list'
            IF (DEBUG) WRITE(MYUNIT,'(A,3I6)') ' intlbfgs> atoms are ',BESTPRESERVEDN(N1),BESTPRESERVEDN(N2),BESTPRESERVEDN(N3)
!           IF (DEBUG) WRITE(MYUNIT,'(A,3I6,A)') ' intlbfgs> full list has length ',NDFORNEWATOM
!           IF (DEBUG) WRITE(MYUNIT,'(20I6)') BESTPRESERVEDN(1:NDFORNEWATOM)

!
! Move the new atom consistently in the local environment of the three active atoms with the
! best preserved absolute distances or the shortest average distances in the end points.
! Check the energies and compare linear interpolation as well, then choose the interpolation
! with the lowest energy.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
            VEC1(1:3)=XYZ(3*(BESTPRESERVEDN(N2)-1)+1:3*(BESTPRESERVEDN(N2)-1)+3) &
  &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(BESTPRESERVEDN(N3)-1)+1:3*(BESTPRESERVEDN(N3)-1)+3) &
  &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)+ &
  &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO

            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVED=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVED(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF

         IF (NCFORNEWATOM.GE.3) THEN
!
! Choose three atoms from the BESTCLOSEST list at random with bias towards the
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
            DUMMY=0.0D0
            DO J1=1,NCFORNEWATOM
!              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
!              DUMMY=DUMMY+1.0D0/(1.0D0*BESTCLOSESTD(J1))
               DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
            ENDDO
            N1=0; N2=0; N3=0
            DO WHILE (N3.EQ.0)
               DUMMY2=0.0D0
               RAN1=DPRAND()*DUMMY
               DO J1=1,NCFORNEWATOM
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTCLOSESTD(J1))
                  DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
                  IF (DUMMY2.GE.RAN1) THEN
                     IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
                     IF (N1.EQ.0) THEN
                        N1=J1
                        EXIT
                     ENDIF
                     IF (N2.EQ.0) THEN
                        N2=J1
                        EXIT
                     ENDIF
                     N3=J1
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF (DEBUG) WRITE(MYUNIT,'(A,3I6,A)') ' intlbfgs> choosing positions ',N1,N2,N3,' in closest list'

            VEC1(1:3)=XYZ(3*(BESTCLOSESTN(N2)-1)+1:3*(BESTCLOSESTN(N2)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(BESTCLOSESTN(N3)-1)+1:3*(BESTCLOSESTN(N3)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)+ &
  &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO

            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVEC=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVEC(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF
!
! Standard linear interpolation, with constraint distance scaled by FRAC.
! Works for FRAC as small as 0.1 with repulsion turned off.
! We use an appropriately weighted displacement from atom CONLIST(1) using the displacements
! in the two end points.
!
         FRAC=1.0D0
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+1))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+2)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+2))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+3))/(INTIMAGE+1)
         ENDDO
         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         IF (DEBUG) WRITE(MYUNIT,'(A,4G15.5)') ' intlbfgs> energies for constrained, preserved, closest, and linear schemes=', &
  &                 ESAVE0,ESAVED,ESAVEC,ETOTAL
         IF ((ETOTAL.LT.ESAVEC).AND.(ETOTAL.LT.ESAVED).AND.(ETOTAL.LT.ESAVE0)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> lowest energy from linear interpolation'
         ELSE IF ((ESAVEC.LT.ESAVED).AND.(ESAVEC.LT.ESAVE0)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using closest atoms'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVEC(1:3,J1)
            ENDDO
            ETOTAL=ESAVEC
         ELSE IF (ESAVED.LT.ESAVE0) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using preserved distances'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVED(1:3,J1)
            ENDDO
            ETOTAL=ESAVED
         ELSE 
            IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using closest constraints'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVE0(1:3,J1)
            ENDDO
            ETOTAL=ESAVE0
         ENDIF
      ENDIF
      NADDED=NADDED+1
      IF (NADDED.LT.NTOADD) GOTO 542
!
! Turn frozen images off for new added atom.
!
!     IF (DEBUG) WRITE(MYUNIT,'(A)') ' intlbfgs> turning off frozen images'
!     IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
      CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
!
! need a new gradient since the active atom has changed !
!
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF

END SUBROUTINE DOADDATOM

SUBROUTINE CHECKPERC(LXYZ,LINTCONSTRAINTTOL,NQCIFREEZE,NCPFIT)
USE COMMONS, ONLY : ATOMACTIVE, NCONSTRAINT, INTFROZEN, CONI, CONJ, CONDISTREF, INTCONMAX, INTCONSTRAINTTOL, &
  &             INTCONSEP, NCONGEOM, CONGEOM, CONIFIX, CONJFIX, CONDISTREFFIX, MYUNIT, &
  &             NCONSTRAINTFIX, PERIODIC, TWOD, RIGID, CONDATT, CONCUT, CONCUTFIX, NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ
IMPLICIT NONE
INTEGER NDIST1(NATOMS), NCYCLE, DMIN1, DMAX1, NUNCON1, J1, J2, J3, NQCIFREEZE, J4, NCPFIT
DOUBLE PRECISION LINTCONSTRAINTTOL, MAXCONDIST, MINCONDIST, DS, DF, LXYZ((3*NATOMS)*2)
DOUBLE PRECISION DSMIN, DSMAX, DSMEAN, D, DIST2, RMAT(3,3)
LOGICAL CHANGED
LOGICAL :: CALLED=.FALSE.
SAVE CALLED

LINTCONSTRAINTTOL=INTCONSTRAINTTOL

IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
!
! Fixed constraints based on congeom file entries
! Just need to adjust the list based on any frozen atoms. We
! want to exclude any constraints between two frozen atoms 
! from the list, because subsequent code depends on this.
!

IF (NCONGEOM.GE.2) THEN
   IF (CALLED.OR.CONDATT) THEN
      J2=0
      DO J1=1,NCONSTRAINTFIX
!
! If called with two minima check that CONCUTFIX is large enough to
! accommodate the separation of the two atoms in both minima.
!
         IF (NCPFIT.EQ.2) THEN
            DF=MAX(ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ(3*(CONIFIX(J1)-1)+1)-LXYZ(3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+2)-LXYZ(3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+3)-LXYZ(3*(CONJFIX(J1)-1)+3))**2)),&
                   ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+1)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+2)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+3)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+3))**2)))
            IF (DF.GT.CONCUTFIX(J1)) THEN
               IF (ABS(DF-CONCUTFIX(J1)).GT.1.0D-6) &
  &                WRITE(MYUNIT,'(A,2I5,3(A,G15.5))') ' checkperc> Increasing con cutoff atoms ', &
  &                CONIFIX(J1),CONJFIX(J1),' from ',CONCUTFIX(J1),' to ',DF,' ref=',CONDISTREFFIX(J1)
               CONCUTFIX(J1)=DF
            ENDIF
         ENDIF
         IF (INTFROZEN(CONIFIX(J1)).AND.INTFROZEN(CONJFIX(J1))) CYCLE
         J2=J2+1
         CONI(J2)=CONIFIX(J1)
         CONJ(J2)=CONJFIX(J1)
         CONDISTREF(J2)=CONDISTREFFIX(J1)
         CONCUT(J2)=CONCUTFIX(J1)
      ENDDO
      NCONSTRAINT=J2
      WRITE(MYUNIT,'(A,I6,A)') ' checkperc> After allowing for frozen atoms there are ',NCONSTRAINT,' constraints'
      RETURN 
   ELSE
!
! Put reference minima in optimal permutational alignment with reference minimum one.
!
      DO J2=2,NCONGEOM
         CALL MINPERMDIST(CONGEOM(1,1:3*NATOMS),CONGEOM(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                       BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,D,DIST2,RIGID,RMAT)
      ENDDO
   ENDIF
   ALLOCATE(CONIFIX(INTCONMAX),CONJFIX(INTCONMAX),CONCUTFIX(INTCONMAX),CONDISTREFFIX(INTCONMAX))
ENDIF

51   NCONSTRAINT=0 
MAXCONDIST=-1.0D0
MINCONDIST=1.0D100
IF (NCONGEOM.LT.2) THEN 
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS

         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
         DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2) 
         IF (DS.GT.5.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
!        IF (DS.GT.15.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
         DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2) 
         IF (DF.GT.5.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
!        IF (DF.GT.15.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
!        IF (2.0D0*ABS(DS-DF)/(DS+DF).LT.LINTCONSTRAINTTOL) THEN
         IF (ABS(DS-DF).LT.LINTCONSTRAINTTOL) THEN
!
!  Add constraint for this distance to the list.
!
            NCONSTRAINT=NCONSTRAINT+1
!           PRINT '(A,2I6,A,I6)','checkperc> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
            IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
            CONI(NCONSTRAINT)=J2
            CONJ(NCONSTRAINT)=J3
            CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
            CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
            IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
            IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
!           IF (DEBUG) PRINT '(A,2I6,A,2F12.2,A,F12.4,A,I8)',' checkperc> constrain distance for atoms ',CONI(NCONSTRAINT), &
! &                 CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
! &                ' # constraints=',NCONSTRAINT
         ENDIF
      ENDDO
   ENDDO
   IF (DEBUG) WRITE(MYUNIT,'(A,I6,2(A,F15.5))') ' checkperc> Total distance constraints=',NCONSTRAINT, &
  &                                     ' shortest=',MINCONDIST,' longest=',MAXCONDIST
ELSE
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS
         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         DSMIN=1.0D100
         DSMAX=-1.0D100
         DSMEAN=0.0D0
         DO J4=1,NCONGEOM
            DS=SQRT((CONGEOM(J4,3*(J2-1)+1)-CONGEOM(J4,3*(J3-1)+1))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+2)-CONGEOM(J4,3*(J3-1)+2))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+3)-CONGEOM(J4,3*(J3-1)+3))**2) 
            IF (DS.GT.DSMAX) DSMAX=DS
            IF (DS.LT.DSMIN) DSMIN=DS
            IF ((J4.GT.1).AND.(ABS(DSMIN-DSMAX).GT.LINTCONSTRAINTTOL)) GOTO 753 ! unconstrained
            DSMEAN=DSMEAN+DS
         ENDDO
!
!  Add constraint for this distance to the list if we make it to here.
!
         NCONSTRAINT=NCONSTRAINT+1
!        PRINT '(A,2I6,A,I6)','checkperc> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DSMAX+DSMIN)/2.0D0 
         CONCUT(NCONSTRAINT)=(DSMAX-DSMIN)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) WRITE(MYUNIT,'(A,2I5,A,2F10.4,A,F12.4,A,I8)') &
  &                       ' checkperc> constrain atoms ',CONI(NCONSTRAINT), &
  &                       CONJ(NCONSTRAINT),' max, min ',DSMAX,DSMIN, &
  &                       ' cutoff=',CONCUT(NCONSTRAINT),' constraints=',NCONSTRAINT
753      CONTINUE
      ENDDO
   ENDDO
   CONIFIX(1:NCONSTRAINT)=CONI(1:NCONSTRAINT)
   CONJFIX(1:NCONSTRAINT)=CONJ(1:NCONSTRAINT)
   CONDISTREFFIX(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
   CONCUTFIX(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)
   NCONSTRAINTFIX=NCONSTRAINT
ENDIF
!
! Check that we have a percolating constraint network. If not, increase the tolerance and try again!
! Calculate minimum number of steps of each atom from number 1 or any frozen atom.
!
NDIST1(1:NATOMS)=1000000
IF (NQCIFREEZE.EQ.0) THEN
   NDIST1(1)=0
ELSE
   DO J1=1,NATOMS
      IF (INTFROZEN(J1)) NDIST1(J1)=0
   ENDDO
ENDIF
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN1=100000
DMAX1=0
NUNCON1=0
DO J1=1,NATOMS
   IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
   DO J2=1,NCONSTRAINT
      IF (CONI(J2).EQ.J1) THEN
         IF (NDIST1(CONJ(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONJ(J2))+1
         ENDIF
      ELSE IF (CONJ(J2).EQ.J1) THEN
         IF (NDIST1(CONI(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONI(J2))+1
         ENDIF
      ENDIF
   ENDDO
   IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
   IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
   IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
ENDDO
IF (CHANGED) GOTO 5
  IF (DEBUG) WRITE(MYUNIT,'(3(A,I8))') ' checkperc> steps to atom 1 converged in ',NCYCLE-1, &
    &               ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
IF (NUNCON1.GT.0) THEN
   LINTCONSTRAINTTOL=LINTCONSTRAINTTOL*1.1D0
   IF (DEBUG) WRITE(MYUNIT,'(A,F15.5)') ' checkperc> increasing the local constraint tolerance parameter to ',LINTCONSTRAINTTOL
   IF (LINTCONSTRAINTTOL.GT.100.0D0) THEN
      WRITE(MYUNIT,'(A,G20.10)') 'checkperc> likely ERROR *** LINTCONSTRAINTTOL=',LINTCONSTRAINTTOL
      STOP
   ENDIF
   GOTO 51
ENDIF
! IF (DEBUG) WRITE(MYUNIT,'(A,F15.5)') ' checkperc> Final constraint tolerance parameter ',LINTCONSTRAINTTOL

! WRITE(MYUNIT,'(A,I6,3(A,F15.5))') ' checkperc> Total distance constraints=',NCONSTRAINT, &
!   &                    ' shortest=',MINCONDIST,' longest=',MAXCONDIST,' tolerance=',LINTCONSTRAINTTOL

CALLED=.TRUE.

END SUBROUTINE CHECKPERC

SUBROUTINE MAKESTEP(NITERDONE,POINT,DIAG,INTIMAGE,SEARCHSTEP,G,GTMP,STP,GDIF,NPT,D,RHO1,ALPHA)
USE COMMONS, ONLY : MUPDATE, INTDGUESS, NATOMS, MYUNIT
IMPLICIT NONE
INTEGER NITERDONE, POINT, BOUND, NPT, D, CP, INTIMAGE, I
DOUBLE PRECISION DIAG(3*NATOMS*INTIMAGE),SEARCHSTEP(0:MUPDATE,(3*NATOMS)*INTIMAGE),G((3*NATOMS)*INTIMAGE), &
  &  GTMP(3*NATOMS*INTIMAGE), GNORM, STP(3*NATOMS*INTIMAGE), YS, GDIF(0:MUPDATE,(3*NATOMS)*INTIMAGE), YY, &
  &  SQ, YR, BETA
DOUBLE PRECISION, DIMENSION(MUPDATE)     :: RHO1,ALPHA
LOGICAL CHANGEIMAGE
SAVE

! WRITE(MYUNIT,'(A,I8)') 'makestep> NITERDONE=',NITERDONE
! WRITE(MYUNIT,'(A,3G20.10)') 'makestep> DIAG: ',DIAG(118:120)
! WRITE(MYUNIT,'(A,3G20.10)') 'makestep> INTDGUESS: ',INTDGUESS
! WRITE(MYUNIT,'(A,3G20.10)') 'makestep> G: ',G(118:120)
! WRITE(MYUNIT,'(A,4I8)') 'NPT,D,NPT/D,POINT=',NPT,D,NPT/D,POINT
! WRITE(MYUNIT,'(A,3G20.10)') 'makestep> GDIF: ',GDIF(NPT/D,118:120)
! WRITE(MYUNIT,'(A,3G20.10)') 'makestep> SEARCHSTEP: ',SEARCHSTEP(NPT/D,118:120)

MAIN: IF (NITERDONE==1) THEN
     POINT = 0
     DIAG(1:D)=INTDGUESS
     SEARCHSTEP(0,1:D)= -G(1:D)*INTDGUESS            ! NR STEP FOR DIAGONAL INVERSE HESSIAN
     GTMP(1:D)        = SEARCHSTEP(0,1:D)
     GNORM            = MAX(SQRT(DOT_PRODUCT(G(1:D),G(1:D))),1.0D-100)
     STP(1:D)         = MIN(1.0D0/GNORM, GNORM) ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
ELSE MAIN
     BOUND=NITERDONE-1
     IF (NITERDONE.GT.MUPDATE) BOUND=MUPDATE
     YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
     IF (YS==0.0D0) YS=1.0D0
    
! Update estimate of diagonal inverse Hessian elements.
! We divide by both YS and YY at different points, so they had better not be zero!

     YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
     IF (YY==0.0D0) YY=1.0D0
!    DIAG = ABS(YS/YY)
     DIAG(1) = YS/YY
      
! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, 
! "Updating quasi-Newton matrices with limited storage",
! Mathematics of Computation, Vol.35, No.151, pp. 773-782

     CP= POINT; IF (POINT==0) CP = MUPDATE
     RHO1(CP)=1.0D0/YS
     GTMP(1:D) = -G(1:D)
     CP= POINT 
                   
     DO I= 1,BOUND 
          CP = CP - 1; IF (CP == -1) CP = MUPDATE - 1
          SQ= DOT_PRODUCT( SEARCHSTEP(CP,1:D),GTMP(1:D) )
          ALPHA(CP+1) = RHO1(CP+1) * SQ
          GTMP(1:D)        = -ALPHA(CP+1)*GDIF(CP,1:D) + GTMP(1:D)
     ENDDO
              
     GTMP(1:D)=DIAG(1)*GTMP(1:D)

     DO I=1,BOUND
          YR= DOT_PRODUCT( GDIF(CP,1:D) , GTMP )
          BETA= RHO1(CP+1)*YR
          BETA= ALPHA(CP+1)-BETA
!         WRITE(MYUNIT,'(A,I8,4G20.10)') 'makestep> I,YR,BETA,RHO1,ALPHA=',I,YR,BETA,RHO1(CP+1),ALPHA(CP+1)
          GTMP(1:D) = BETA*SEARCHSTEP(CP,1:D) + GTMP(1:D)
          CP=CP+1
!         IF (CP==M) CP=0
          IF (CP==MUPDATE) CP=0
     ENDDO
              
     STP(1:D) = 1.0D0
ENDIF MAIN

!  Store the new search direction
IF (NITERDONE.GT.1) SEARCHSTEP(POINT,1:D)=GTMP(1:D)

END SUBROUTINE MAKESTEP
