!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! Objective - to symmetrise a set of input coordinates using a
! continuous symmetry measure and constructing the closest structure with
! the point group in question. In fact, this averaged structure may not have 
! the exact point group because we do not insist that the permutations
! corresponding to minimum distances are themselves organised symmetrically. 
!
SUBROUTINE SYMMETRYCSM(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL, &
  &         LNQUENCH,EBEST,BESTCOORDS,JBEST)

USE COMMONS
USE porfuncs
IMPLICIT NONE
DOUBLE PRECISION :: LCOORDS(3*NATOMS), T0, SCREENC(3*NATOMS), LEBEST, QBEST(3*NATOMS), VATBEST(NATOMS), VATORIG(NATOMS), POO
DOUBLE PRECISION :: TIME, EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), POTEL, CSMBEST, OLCOORDS(3*NATOMS), CSMSTEP, EREAL
DOUBLE PRECISION :: DUMMY(3*NATOMS), CSMIMAGESSAVE(3*NATOMS*CSMGPINDEX), CSMPMATSAVE(3,3), CSMEBEST
DOUBLE PRECISION :: XTEMP(3*NATOMS), CSMVALUE, RANDOM, ANGLE, COST, SINT, TY, TZ, TX, DPRAND, CSMRMS, RMAT(3,3), DIST2
DOUBLE PRECISION :: QORIG(3*NATOMS), ASTEPSAVE, CSMAVNEW(3*NATOMS), GSAVE(3*NATOMS), AA(3)
INTEGER :: ISTAT, LNQUENCH, JP, BRUN,QDONE,JBEST(NPAR), ITERATIONS, J1, J2, CSMIT, J3, NDONE
LOGICAL :: CHANGEDE, LDEBUG
INTEGER :: NQTOT, NQTOTSAVE, NSYMCALL
COMMON /MYPOT/ POTEL
COMMON /TOT/ NQTOT

LDEBUG=DEBUG
LDEBUG=.TRUE.
NSYMCALL=NSYMCALL+1   ! NSYMCALL should be the number of consecutive calls to symmetry with this minimum
! PRINT '(A,5I)','NSYMCALL,NORDER,NORBIT,LASTORBIT,NORBIT-LASTORBIT+1=',NSYMCALL,NORDER,NORBIT,LASTORBIT,NORBIT-LASTORBIT+1
! IF (NSYMCALL.GT.1) THEN
!    IF (NSYMCALL.GT.1) THEN
!       IF (LDEBUG) WRITE(MYUNIT, '(A)') 'maximum consecutive calls to symmetry reached for this minimum'
!       RETURN 
!    ENDIF
! ENDIF
NSYMREM=0
CHANGEDE=.FALSE.
CSMEBEST=EPREV(JP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Accept the CSM minimised geometry whatever!
! CSMEBEST=1.0D100
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
QBEST(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
VATBEST(1:NATOMS)=VAT(1:NATOMS,JP)
NQTOTSAVE=NQTOT
LCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
VATORIG(1:NATOMS)=VAT(1:NATOMS,JP)
QORIG(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
!
!  Since the structure in COORDS is changing we need to recalculate the 
!  normalisation.
!

CSMNORM=0.0D0
DO J1=1,NATOMS
   CSMNORM=CSMNORM+LCOORDS(3*(J1-1)+1)**2+LCOORDS(3*(J1-1)+2)**2+LCOORDS(3*(J1-1)+3)**2
ENDDO
CSMNORM=2*CSMGPINDEX*CSMNORM

!
!  (1) First do CSMSTEPS of basin-hopping for the CSM value of the given point group.
!      Accept downhill moves only.
!  (2) Then do CSMRELAX basin-hopping steps to relax the averaged structure, which
!      ideally has the given point group symmetry (but probably doesn;t).
!  (3) Accept/reject lowest energy structure found according to usual Metropolis condition.
!

CSMBEST=1.0D100
OLCOORDS(1:3*NATOMS)=LCOORDS(3*NATOMS)
CSMSTEP=0.5D0
DO J2=1,CSMSTEPS
   IF (PERMDIST) THEN
      DO J1=1,CSMGPINDEX
         XTEMP(1:3*NATOMS)=LCOORDS(1:3*NATOMS)
         CALL CSMROT(XTEMP,DUMMY,1,J1)
         CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,EREAL,DIST2,RIGID,RMAT)
         CALL CSMROT(DUMMY,XTEMP,-1,J1) ! need to rotate the permuted rotated images back to the reference orientation
         CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
      ENDDO
   ELSE
      DO J1=1,CSMGPINDEX
         CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=LCOORDS(1:3*NATOMS)
      ENDDO
   ENDIF
!
! At this point we have the best permutations for each point group operation,
! with CSMIMAGES containing the permuted coordinates of X rotated back into
! its reference orientation, i.e. the point group mapping has been reversed.
!
   CALL CSMMIN(LCOORDS,CSMVALUE,CSMRMS,CSMIT)
   IF (CSMVALUE.LT.CSMBEST) THEN
      CSMBEST=CSMVALUE
      IF (LDEBUG) WRITE(MYUNIT,'(A,I6,A,G20.10,A,I6)') 'symmetrizecsm> accepted CSM geometry ',J2,' for CSM=', &
 &                                                     CSMVALUE,' iterations=',CSMIT
      OLCOORDS(1:3*NATOMS)=LCOORDS(1:3*NATOMS)
      CSMIMAGESSAVE(1:3*NATOMS*CSMGPINDEX)=CSMIMAGES(1:3*NATOMS*CSMGPINDEX)
      CSMPMATSAVE(1:3,1:3)=CSMPMAT(1:3,1:3)
   ELSE
      IF (LDEBUG) WRITE(MYUNIT,'(A,I6,A,G20.10,A,I6)') 'symmetrizecsm> rejected CSM geometry ',J2,' for CSM=', &
 &                                                     CSMVALUE,' iterations=',CSMIT
      LCOORDS(1:3*NATOMS)=OLCOORDS(1:3*NATOMS)
   ENDIF
   IF (CSMVALUE.LT.1.0D-6) EXIT ! we can;t do better than exact symmetry!
   RANDOM=(DPRAND()-0.5D0)*2.0D0
   ANGLE=RANDOM*CSMSTEP ! in radians
   COST=COS(ANGLE)
   SINT=SIN(ANGLE)
!
! Random rotation about the x axis.
!
   DO J1=1,NATOMS
      TY= COST*LCOORDS(3*(J1-1)+2)+SINT*LCOORDS(3*(J1-1)+3)
      TZ=-SINT*LCOORDS(3*(J1-1)+2)+COST*LCOORDS(3*(J1-1)+3)
      LCOORDS(3*(J1-1)+2)=TY
      LCOORDS(3*(J1-1)+3)=TZ
   ENDDO
   RANDOM=(DPRAND()-0.5D0)*2.0D0
   ANGLE=RANDOM*CSMSTEP ! in radians
   COST=COS(ANGLE)
   SINT=SIN(ANGLE)
!
! Random rotation about the y axis.
!
   DO J1=1,NATOMS
      TX= COST*LCOORDS(3*(J1-1)+1)+SINT*LCOORDS(3*(J1-1)+3)
      TZ=-SINT*LCOORDS(3*(J1-1)+1)+COST*LCOORDS(3*(J1-1)+3)
      LCOORDS(3*(J1-1)+1)=TX
      LCOORDS(3*(J1-1)+3)=TZ
   ENDDO
!
! Random rotation about the z axis.
! 
   RANDOM=(DPRAND()-0.5D0)*2.0D0
   ANGLE=RANDOM*CSMSTEP ! in radians
   COST=COS(ANGLE)
   SINT=SIN(ANGLE)
   DO J1=1,NATOMS
      TX= COST*LCOORDS(3*(J1-1)+1)+SINT*LCOORDS(3*(J1-1)+2)
      TY=-SINT*LCOORDS(3*(J1-1)+1)+COST*LCOORDS(3*(J1-1)+2)
      LCOORDS(3*(J1-1)+1)=TX
      LCOORDS(3*(J1-1)+2)=TY
   ENDDO
ENDDO
!
! Construct averaged geometry from the orientation with the lowest CSM.
!
CSMIMAGES(1:3*NATOMS*CSMGPINDEX)=CSMIMAGESSAVE(1:3*NATOMS*CSMGPINDEX)
CSMPMAT(1:3,1:3)=CSMPMATSAVE(1:3,1:3)
CSMAV(1:3*NATOMS)=0.0D0
DO J2=1,CSMGPINDEX
! 
! rotate permuted image to best orientation with CSMPMAT
! apply point group operation J2
!
   DO J3=1,NATOMS
      XTEMP(3*(J3-1)+1)=CSMPMAT(1,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
  &                    +CSMPMAT(1,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
  &                    +CSMPMAT(1,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
      XTEMP(3*(J3-1)+2)=CSMPMAT(2,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
  &                    +CSMPMAT(2,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
  &                    +CSMPMAT(2,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
      XTEMP(3*(J3-1)+3)=CSMPMAT(3,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
  &                    +CSMPMAT(3,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
  &                    +CSMPMAT(3,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
   ENDDO
   CALL CSMROT(XTEMP,DUMMY,1,J2)
   CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)+DUMMY(1:3*NATOMS)
ENDDO
CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)/CSMGPINDEX
! IF (LDEBUG) THEN
!     OPEN(UNIT=1,FILE='CSMav_structures.xyz',STATUS='UNKNOWN')
!    WRITE(1,'(I6)') NATOMS
!    WRITE(1,'(A,I6)') 'coordinates after quench ',NQ(JP)-1
!    WRITE(1,'(A3,3G20.10)') ('LA ',LCOORDS(3*(J1-1)+1),LCOORDS(3*(J1-1)+2),LCOORDS(3*(J1-1)+3),J1=1,NATOMS)
!    WRITE(1,'(I6)') NATOMS
!    WRITE(1,'(A)') 'averaged coordinates '
!    WRITE(1,'(A3,3G20.10)') ('LA ',CSMAV(3*(J1-1)+1),CSMAV(3*(J1-1)+2),CSMAV(3*(J1-1)+3),J1=1,NATOMS)
!     CLOSE(1)
! ENDIF
OLCOORDS(1:3*NATOMS)=CSMAV(1:3*NATOMS)
COORDS(1:3*NATOMS,JP)=CSMAV(1:3*NATOMS)

! !
! !  Repeat averaging for the averaged structure!
! !
! NDONE=0
! 5 CONTINUE
! CALL CENTRE2(CSMAV)
! NDONE=NDONE+1
! IF (PERMDIST) THEN
!    CSMVALUE=0.0D0
!    DO J1=1,CSMGPINDEX
!       XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
!       CALL CSMROT(XTEMP,DUMMY,1,J1)
!       CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,EREAL,DIST2,RIGID,RMAT)
!       CSMVALUE=CSMVALUE+EREAL**2
!       CALL CSMROT(DUMMY,XTEMP,-1,J1) ! need to rotate the permuted rotated images back to the reference orientation
!       CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
!    ENDDO
! ELSE
!    DO J1=1,CSMGPINDEX
!       CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=CSMAV(1:3*NATOMS)
!    ENDDO
! ENDIF
! CSMNORM=0.0D0
! DO J1=1,NATOMS
!    CSMNORM=CSMNORM+CSMAV(3*(J1-1)+1)**2+CSMAV(3*(J1-1)+2)**2+CSMAV(3*(J1-1)+3)**2
! ENDDO
! CSMNORM=2*CSMGPINDEX*CSMNORM
! CSMVALUE=CSMVALUE/CSMNORM
! !
! ! Checked that CSMPOTGRAD value for CSM agrees with above.
! !
! ! AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0
! ! CALL CSMPOTGRAD(CSMAV,AA,POO,.FALSE.,GSAVE)
! 
! CSMAVNEW(1:3*NATOMS)=0.0D0
! DO J2=1,CSMGPINDEX
!    XTEMP(1:3*NATOMS)=CSMIMAGES(1+3*NATOMS*(J2-1):3*NATOMS*J2)
!    CALL CSMROT(XTEMP,DUMMY,1,J2)
!    CSMAVNEW(1:3*NATOMS)=CSMAVNEW(1:3*NATOMS)+DUMMY(1:3*NATOMS)
! ENDDO
! CSMAVNEW(1:3*NATOMS)=CSMAVNEW(1:3*NATOMS)/CSMGPINDEX
! ! WRITE(3,'(I6)') NATOMS
! ! WRITE(3,'(A,2G20.10)') 'REaveraged coordinates, CSM=',CSMVALUE
! ! WRITE(MYUNIT,'(A,3G20.10)') 'REaveraged coordinates, CSM=',CSMVALUE
! ! WRITE(3,'(A3,3G20.10)') ('LA ',CSMAVNEW(3*(J1-1)+1),CSMAVNEW(3*(J1-1)+2),CSMAVNEW(3*(J1-1)+3),J1=1,NATOMS)
! CSMAV(1:3*NATOMS)=CSMAVNEW(1:3*NATOMS)
! IF (NDONE.LT.10) GOTO 5
! ! 
! ! STOP
! WRITE(3,'(I6)') NATOMS
! WRITE(3,'(A,2G20.10)') 'REaveraged coordinates, CSM=',CSMVALUE
! WRITE(3,'(A3,3G20.10)') ('LA ',CSMAV(3*(J1-1)+1),CSMAV(3*(J1-1)+2),CSMAV(3*(J1-1)+3),J1=1,NATOMS)
! OLCOORDS(1:3*NATOMS)=CSMAV(1:3*NATOMS)
! COORDS(1:3*NATOMS,JP)=CSMAV(1:3*NATOMS)

ASTEPSAVE=ASTEP(JP)
ASTEP(JP)=0.0D0 ! turn off angular steps for CSM relaxation
DO J1=1,CSMQUENCHES
   NQTOT=NQTOT+1
   NQ(JP)=NQ(JP)+1
   CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
   CALL MYCPU_TIME(T0)

   IF (NPAR.GT.1) THEN
      WRITE(MYUNIT,'(A,I1,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,I5)') '[',JP,']Qu ',NQ(JP),' E=', &
  &        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',T0-TSTART,' in CSM relaxation ',J1
   ELSE
      WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,I5)') 'Qu ',NQ(JP),' E=', &
  &        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',T0-TSTART,' in CSM relaxation ',J1
   ENDIF
!
! Accept downhill moves.
!
   IF (CSMEBEST-POTEL.GT.0.0D0) THEN
      CSMEBEST=POTEL
      OLCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
      QBEST(1:3*NATOMS)=OLCOORDS(1:3*NATOMS)
      VATBEST(1:NATOMS)=VAT(1:NATOMS,JP)
   ELSE
      COORDS(1:3*NATOMS,JP)=OLCOORDS(1:3*NATOMS)
   ENDIF
   IF (HIT) EXIT 
   CALL TAKESTEP(JP)
ENDDO
ASTEP(JP)=ASTEPSAVE ! restore angular step

IF (EPREV(JP)-CSMEBEST.GT.ECONV) THEN
   WRITE(MYUNIT,'(A,G20.10)') 'symmetrisecsm>     changing current minimum to energy ',CSMEBEST
   CHANGEDE=.TRUE. 
   EPREV(JP)=CSMEBEST
   POTEL=CSMEBEST
   VAT(1:NATOMS,JP)=VATBEST(1:NATOMS)
   COORDS(1:3*NATOMS,JP)=QBEST(1:3*NATOMS)
   IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
   COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
   VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
ELSE
   WRITE(MYUNIT,'(A,G20.10)') 'symmetrisecsm> not changing current minimum energy     ',EPREV(JP)
   COORDS(1:3*NATOMS,JP)=QORIG(1:3*NATOMS)
ENDIF

RETURN

END SUBROUTINE SYMMETRYCSM

