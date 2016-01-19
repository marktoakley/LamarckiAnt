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
! Objective - to symmetrise a set of input coordinates.
!
SUBROUTINE SYMMETRY2(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL)
USE COMMONS
USE porfuncs
IMPLICIT NONE
DOUBLE PRECISION :: LCOORDS(3*NATOMS), CMDIST(NATOMS), ORBDIST(NATOMS), T0
DOUBLE PRECISION :: EBEST, QBEST(3*NATOMS), VATBEST(NATOMS)
DOUBLE PRECISION :: CMX(2), CMY(2), CMZ(2), CM(3), CMSAVE(3), COREVT(NATOMS), OTHERVT(NATOMS)
DOUBLE PRECISION :: DUMMY, XMASS, YMASS, ZMASS
DOUBLE PRECISION :: ORIGIN(3)
DOUBLE PRECISION :: TIME, LTOLD
INTEGER :: J1, J2, NDUMMY, I, NORBIT, J, ISTART, LARGESIZE, MYUNIT2, MYUNIT3, HPTGRPSAVE
INTEGER :: JP, NATOMSCORE, NORBITSCORE
INTEGER :: NFLOAT
INTEGER :: NINDEX(NATOMS), ORBSIZE(NATOMS)
DOUBLE PRECISION :: SCOORDS(3*NATOMS)
! DOUBLE PRECISION :: ORBSORT(3*NATOMS), LCOORDSSORTED(3*NATOMS), RCOORDSSORTED(3*NATOMS)
DOUBLE PRECISION :: CORECOORDS(3*NATOMS), OTHERCOORDS(3*NATOMS)
DOUBLE PRECISION :: DENOM
LOGICAL :: CHANGEDE, LDEBUG
DOUBLE PRECISION :: POTEL, SCREENC(3*NATOMS), GENMAT(100,3,3), GENMATSAVE(100,3,3), X(3*NATOMS)
DOUBLE PRECISION :: DUMMYE, DUMMYGRAD(3*NATOMS)
INTEGER :: NCORESAVE
INTEGER :: NQTOT, QDONE, BRUN, ITERATIONS, NQTOTSAVE, IGEN, IGENSAVE, ORBSYM
INTEGER :: IPRNT
INTEGER :: NSYMCALL
CHARACTER(LEN=4) FPGRP, POINTGROUP
CHARACTER(LEN=6) JPSTRING
DOUBLE PRECISION, PARAMETER :: EPS=1.0D-10
COMMON /MYPOT/ POTEL
COMMON /TOT/ NQTOT
SAVE NORBIT

! DO J1=1,NPAR
  ! WRITE(MYUNIT,*) 'symmetry2> J1,SHELLMOVES,PTGROUP',J1,SHELLMOVES(J1),PTGROUP(J1)
! ENDDO
NCORESAVE=NCORE(JP)
WRITE(JPSTRING,'(I6)') JP
MYUNIT2=NPAR+MYUNIT
MYUNIT3=2*NPAR+MYUNIT+1
IPRNT=0
LDEBUG=DEBUG
NSYMCALL=NSYMCALL+1   ! NSYMCALL should be the number of consecutive calls to symmetry with this minimum
IF (NSYMCALL.GT.2) THEN
   IF (LDEBUG) WRITE(MYUNIT, '(A)') 'symmetry2> maximum consecutive calls to symmetry reached for this minimum'
   WRITE(MYUNIT, '(A)') 'symmetry2> maximum consecutive calls to symmetry reached for this minimum'
   RETURN 
ENDIF
NSYMREM=0
CHANGEDE=.FALSE.
IF (LDEBUG) IPRNT=11
CALL MYCPU_TIME(T0)
EBEST=EPREV(JP)
QBEST(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
VATBEST(1:NATOMS)=VAT(1:NATOMS,JP)
NQTOTSAVE=NQTOT
LCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)

!   DO J2=1,3*NATOMS
!      X(J2)=COORDS(J2,JP)
!   ENDDO
!   CALL POTENTIAL(X,DUMMYGRAD,DUMMYE,.FALSE.,.FALSE.)
!   DO J2=1,NATOMS
!      IF (VT(J2).NE.0.0D0) THEN
!         IF (ABS((VT(J2)-VAT(J2,JP))/VT(J2)).GT.0.01D0) THEN
!            WRITE(MYUNIT,'(A,I8,2F15.5)') 'symmetry2> A J2,VAT,VT(J2)=',J2,VAT(J2,JP),VT(J2)
!            STOP
!         ENDIF
!      ENDIF
!   ENDDO


LTOLD=SYMTOL2
ORIGIN(1:3)=1.0D0
cmloop: DO J1=1,200
   XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
   DENOM=0.0D0
   DO I=1,NATOMS
      IF (J1.GT.1) THEN
         DUMMY=EXP(-DISTFAC*CMDIST(I))
      ELSE
         DUMMY=1.0D0
      ENDIF
      XMASS=XMASS+LCOORDS(3*(I-1)+1)*DUMMY
      YMASS=YMASS+LCOORDS(3*(I-1)+2)*DUMMY
      ZMASS=ZMASS+LCOORDS(3*(I-1)+3)*DUMMY
      DENOM=DENOM+DUMMY
   ENDDO
   DUMMY=SQRT( (ORIGIN(1)-XMASS/DENOM)**2+(ORIGIN(2)-YMASS/DENOM)**2+(ORIGIN(3)-ZMASS/DENOM)**2)
   ORIGIN(1)=XMASS/DENOM;  ORIGIN(2)=YMASS/DENOM; ORIGIN(3)=ZMASS/DENOM
   DO J=1,NATOMS
      CMDIST(J)=SQRT((LCOORDS(3*(J-1)+1)-ORIGIN(1))**2+ &
                     (LCOORDS(3*(J-1)+2)-ORIGIN(2))**2+ &
                     (LCOORDS(3*(J-1)+3)-ORIGIN(3))**2)
   ENDDO
   IF (LDEBUG) WRITE(MYUNIT, '(A,I5,3G20.10)') 'symmetry2> cycle, origin:',J1,ORIGIN(1:3)
   IF ((J1.GT.1).AND.(DUMMY.LT.1.0D-4)) EXIT cmloop 
ENDDO cmloop

DO J=1,NATOMS
   NINDEX(J)=J
ENDDO

IF (LDEBUG) WRITE(MYUNIT, '(A,3F15.5)') 'symmetry2> initial centre of mass: ',ORIGIN(1:3)
SCOORDS(1:3*NATOMS)=LCOORDS(1:3*NATOMS)
CALL PIKSR2(NATOMS,CMDIST,NINDEX) ! sorts CMDIST and NINDEX
! CALL GETORBITS(NORBIT,ORBDIST,CMDIST,NATOMS,ORBSIZE,SYMTOL1,LARGESIZE,NINDEX,LCOORDS,LDEBUG)
! Find the largest gap between CM distances and use the centre-of-mass of all the atoms up to that point
CMX(1)=0.0D0; CMY(1)=0.0D0; CMZ(1)=0.0D0
CMX(2)=LCOORDS(3*(NINDEX(1)-1)+1); CMY(2)=LCOORDS(3*(NINDEX(1)-1)+2); CMZ(2)=LCOORDS(3*(NINDEX(1)-1)+3)
! Store sum of centres of mass from previous largest gap to current in CMx(2).
GMAX=-1.0D0
DO J=2,NATOMS
   IF (ABS(CMDIST(J)-CMDIST(J-1)).GT.GMAX) THEN
      GMAX=ABS(CMDIST(J)-CMDIST(J-1))
      NDUMMY=J-1
      CMX(1)=CMX(1)+CMX(2); CMY(1)=CMY(1)+CMY(2); CMZ(1)=CMZ(1)+CMZ(2)
      CMX(2)=LCOORDS(3*(NINDEX(J)-1)+1); CMY(2)=LCOORDS(3*(NINDEX(J)-1)+2); CMZ(2)=LCOORDS(3*(NINDEX(J)-1)+3)
   ELSE
      CMX(2)=CMX(2)+LCOORDS(3*(NINDEX(J)-1)+1); CMY(2)=CMY(2)+LCOORDS(3*(NINDEX(J)-1)+2); CMZ(2)=CMZ(2)+LCOORDS(3*(NINDEX(J)-1)+3)
   ENDIF
ENDDO
CMX(1)=CMX(1)/NDUMMY; CMY(1)=CMY(1)/NDUMMY; CMZ(1)=CMZ(1)/NDUMMY

IF (LDEBUG) WRITE(MYUNIT, '(A,3F15.5)') 'symmetry2> origin will be moved to centre of mass for atoms up to biggest gap: ', &
  &                           CMX(1),CMY(1),CMZ(1)
IF (LDEBUG) WRITE(MYUNIT, '(A,I6,A,F15.5,A)') 'symmetry2> number of atoms=',NDUMMY
ORBSIZE(1)=NDUMMY
IF (LDEBUG) WRITE(MYUNIT, '(12I5)')  (NINDEX(J1),J1=1,ORBSIZE(1))
ORIGIN(1)=CMX(1)
ORIGIN(2)=CMY(1)
ORIGIN(3)=CMZ(1)
DO J=1,NATOMS
   NINDEX(J)=J
   LCOORDS(3*(J-1)+1)=LCOORDS(3*(J-1)+1)-ORIGIN(1)
   LCOORDS(3*(J-1)+2)=LCOORDS(3*(J-1)+2)-ORIGIN(2)
   LCOORDS(3*(J-1)+3)=LCOORDS(3*(J-1)+3)-ORIGIN(3)
   CMDIST(J)=SQRT(LCOORDS(3*(J-1)+1)**2+LCOORDS(3*(J-1)+2)**2+LCOORDS(3*(J-1)+3)**2)
ENDDO
CALL PIKSR2(NATOMS,CMDIST,NINDEX) ! sorts CMDIST and NINDEX again
CALL GETORBITS(NORBIT,ORBDIST,CMDIST,NATOMS,ORBSIZE,SYMTOL1,LARGESIZE,NINDEX,LCOORDS,LDEBUG)
IF (LDEBUG) THEN
   ISTART=1
   OPEN(UNIT=MYUNIT2,FILE='orbits.'//TRIM(ADJUSTL(JPSTRING)) // '.xyz',STATUS='UNKNOWN')
   DO J1=1,NORBIT
      WRITE(MYUNIT2,'(I5)') SUM(ORBSIZE(1:J1))
      WRITE(MYUNIT2,'(A)') ' '
      DO J2=ISTART,ISTART+ORBSIZE(J1)-1
         WRITE(MYUNIT2,'(A2,2X,3G20.10)') 'LA',LCOORDS(3*(NINDEX(J2)-1)+1:3*(NINDEX(J2)-1)+3)
      ENDDO
      DO J2=1,ISTART-1
         WRITE(MYUNIT2,'(A2,2X,3G20.10)') 'LB',LCOORDS(3*(NINDEX(J2)-1)+1:3*(NINDEX(J2)-1)+3)
      ENDDO
      ISTART=ISTART+ORBSIZE(J1)
   ENDDO
   CLOSE(MYUNIT2)
ENDIF

IF (LARGESIZE.EQ.1) THEN
   IF (LDEBUG) THEN
      CALL MYCPU_TIME(TIME)
      WRITE(MYUNIT, '(A,F15.2)') 'symmetry2> No nontrivial orbits - return from symmetry2, time taken=',TIME-T0
   ENDIF
   RETURN
ENDIF

! Examine how the point group changes as we include more orbits.
 
ISTART=1
NATOMSCORE=0
NORBITSCORE=0
IGENSAVE=0
HPTGRPSAVE=0
symcore: DO J1=1,NORBIT
   DO J2=ISTART,ISTART+ORBSIZE(J1)-1
      CORECOORDS(3*(J2-1)+1:3*(J2-1)+3)=LCOORDS(3*(NINDEX(J2)-1)+1:3*(NINDEX(J2)-1)+3)
   ENDDO
   ISTART=ISTART+ORBSIZE(J1)
   NATOMSCORE=NATOMSCORE+ORBSIZE(J1)
   SCOORDS(1:3*NATOMSCORE)=CORECOORDS(1:3*NATOMSCORE) 
   IF (NATOMSCORE.GT.3) THEN
      DUMMY=MATDIFF
      CALL PTGRP(SCOORDS,NATOMSCORE,LDEBUG,SYMTOL1,SYMTOL2,SYMTOL3,GENMAT,IGEN,FPGRP,CM,DUMMY) ! SCOORDS is changed!
      IF (LDEBUG) WRITE(MYUNIT, '(A,I5,A,I5,A,I4,A,A4,A,I6)') 'symmetry2> number of orbits=',J1,' number of atoms=',NATOMSCORE, &
                          ' number of generators=',IGEN,' point group= ',FPGRP,' order=',HPTGRP
      IF (HPTGRP.GE.HPTGRPSAVE) THEN
         HPTGRPSAVE=HPTGRP
         ORBSYM=J1
         CMSAVE(1:3)=CM(1:3)
         NCORE(JP)=NATOMSCORE
         NORBITSCORE=J1
         IGENSAVE=IGEN
         POINTGROUP=FPGRP
         GENMATSAVE(1:IGEN,1:3,1:3)=GENMAT(1:IGEN,1:3,1:3)
         IF ((HPTGRP.EQ.1).AND.(NCORE(JP).GE.20)) THEN
            NCORE(JP)=0
            EXIT symcore ! looks like there is no symmetry
         ENDIF
      ELSEIF (HPTGRP.GE.1) THEN ! symmetry has decreased
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I1,3A,I6,A,I6,A)') '[',JP,']symmetry2> highest symmetry ',POINTGROUP,' for ',NORBITSCORE, &
  &                        ' orbits and ',NCORE(JP),' atoms'
         ELSE
            WRITE(MYUNIT,'(3A,I6,A,I6,A)') 'symmetry2> highest symmetry ',POINTGROUP,' for ',NORBITSCORE, &
  &                        ' orbits and ',NCORE(JP),' atoms'
         ENDIF
         EXIT symcore
      ENDIF
   ELSE
      IF (LDEBUG) WRITE(MYUNIT, '(A,I5,A,I5)') 'symmetry2> number of orbits=',NORBITSCORE,' number of atoms=',NATOMSCORE
   ENDIF
ENDDO symcore

IF (HPTGRPSAVE.EQ.1) THEN
   IF (LDEBUG) WRITE(MYUNIT, '(A)') 'symmetry2> no symmetry detected'
   NCORE(JP)=NCORESAVE
   RETURN
ELSE
   IF (NCORE(JP).EQ.NATOMS) THEN
      RETURN
   ENDIF
!  WRITE(MYUNIT, '(A,A,A,I4)') 'symmetry2> symmetry analysis for point group ',TRIM(ADJUSTL(POINTGROUP)), &
!   &                   ' number of generators=',IGENSAVE
ENDIF

! Move the CORECOORDS centre of mass to the origin.

DO J1=1,NATOMS
   LCOORDS(3*(J1-1)+1)=LCOORDS(3*(J1-1)+1)-CMSAVE(1)
   LCOORDS(3*(J1-1)+2)=LCOORDS(3*(J1-1)+2)-CMSAVE(2)
   LCOORDS(3*(J1-1)+3)=LCOORDS(3*(J1-1)+3)-CMSAVE(3)
ENDDO
DO J1=1,NCORE(JP)
   CORECOORDS(3*(J1-1)+1)=CORECOORDS(3*(J1-1)+1)-CMSAVE(1)
   CORECOORDS(3*(J1-1)+2)=CORECOORDS(3*(J1-1)+2)-CMSAVE(2)
   CORECOORDS(3*(J1-1)+3)=CORECOORDS(3*(J1-1)+3)-CMSAVE(3)
ENDDO

! Move coordinates into CORE and OTHER vectors. Define NCORE(JP), the number of core atoms.

ISTART=1
DO J1=1,NORBITSCORE
   DO J2=ISTART,ISTART+ORBSIZE(J1)-1
      CORECOORDS(3*(J2-1)+1:3*(J2-1)+3)=LCOORDS(3*(NINDEX(J2)-1)+1:3*(NINDEX(J2)-1)+3)
      COREVT(J2)=VAT(NINDEX(J2),JP)
   ENDDO
   ISTART=ISTART+ORBSIZE(J1)
ENDDO
NCORE(JP)=ISTART-1
DO J1=NORBITSCORE+1,NORBIT
   DO J2=ISTART,ISTART+ORBSIZE(J1)-1
      OTHERCOORDS(3*(J2-NCORE(JP)-1)+1:3*(J2-NCORE(JP)-1)+3)=LCOORDS(3*(NINDEX(J2)-1)+1:3*(NINDEX(J2)-1)+3)
      OTHERVT(J2-NCORE(JP))=VAT(NINDEX(J2),JP)
   ENDDO
   ISTART=ISTART+ORBSIZE(J1)
ENDDO

NFLOAT=NATOMS-NCORE(JP)

IF (LDEBUG) THEN
   OPEN(MYUNIT2,FILE='coreplusother.' // TRIM(ADJUSTL(JPSTRING)) // '.xyz',STATUS='UNKNOWN')
   WRITE(MYUNIT2,*) NATOMS
   WRITE(MYUNIT2,'(A)') ' '
   IF (NCORE(JP).GT.0) WRITE(MYUNIT2,'(A2,3X,3F20.10)') ('LA',CORECOORDS(3*(J2-1)+1:3*(J2-1)+3),J2=1,NCORE(JP))
   DO J1=1,NFLOAT
      WRITE(MYUNIT2,'(A2,3X,3F20.10)') 'LB',OTHERCOORDS(3*(J1-1)+1:3*(J1-1)+3)
      NDUMMY=NDUMMY+1
   ENDDO
   CLOSE(MYUNIT2)
ENDIF
IF (NCORE(JP).EQ.NATOMS) RETURN
! IF (LDEBUG) WRITE(MYUNIT,'(A,3I8)') 'symmetry2> NCORE(JP),NCORESAVE,NFLOAT=',NCORE(JP),NCORESAVE,NFLOAT
IF (NCORE(JP).LE.NCORESAVE) THEN
   NCORE(JP)=NCORESAVE
   RETURN
ENDIF

DO J1=1,NPAR
   IF (J1.EQ.JP) CYCLE
!  WRITE(MYUNIT,*) 'symmetry2> J1,SHELLMOVES,PTGROUP,POINTGROUP=',J1,SHELLMOVES(J1),PTGROUP(J1),POINTGROUP
   IF (.NOT.SHELLMOVES(J1)) CYCLE
   IF (POINTGROUP.EQ.PTGROUP(J1)) THEN
      IF (NPAR.GT.1) THEN
         WRITE(MYUNIT,'(A,I1,3A,I6)') '[',JP,']symmetry2> point group ',POINTGROUP,' coincides with run ',J1
      ELSE
         WRITE(MYUNIT,'(3A,I6)') 'symmetry2> point group ',POINTGROUP,' coincides with run ',J1
      ENDIF
      NCORE(JP)=0
      NSURFMOVES(JP)=0
      SHELLMOVES(JP)=.FALSE. ! in case this was true and the symmetry has just increased to that of another run
      RETURN
   ENDIF
ENDDO
NSURFMOVES(JP)=0
SHELLMOVES(JP)=.TRUE.
PTGROUP(JP)=POINTGROUP

COORDS(3*NFLOAT+1:3*NATOMS,JP)=CORECOORDS(1:3*NCORE(JP))
VAT(NFLOAT+1:NATOMS,JP)=COREVT(1:NCORE(JP))
DO J1=1,NFLOAT
   COORDS(3*(J1-1)+1:3*(J1-1)+3,JP)=OTHERCOORDS(3*(J1-1)+1:3*(J1-1)+3)
ENDDO
VAT(1:NFLOAT,JP)=OTHERVT(1:NFLOAT)

IF (NPAR.GT.1) THEN
   WRITE(MYUNIT,'(A,I1,A,I8,A,I8,3A,I8,A)') '[',JP,']symmetry2> starting a block of ',SHELLMOVEMAX,' moves with ', &
  &      NCORE(JP),' atoms unperturbed, core symmetry ',POINTGROUP,' after ',NQ(JP),' quenches'
ELSE
   WRITE(MYUNIT,'(A,I8,A,I8,3A,I8,A)') 'symmetry2> starting a block of ',SHELLMOVEMAX,' moves with ', &
  &      NCORE(JP),' atoms unperturbed, core symmetry ',POINTGROUP,' after ',NQ(JP),' quenches'
ENDIF
IF (DEBUG) THEN ! check for overlaps
   DO J1=1,NATOMS ! check overlaps with current orbit
      DO J2=J1+1,NATOMS 
         DUMMY=SQRT((COORDS(3*(J1-1)+1,JP)-COORDS(3*(J2-1)+1,JP))**2 &
                   +(COORDS(3*(J1-1)+2,JP)-COORDS(3*(J2-1)+2,JP))**2 &
                   +(COORDS(3*(J1-1)+3,JP)-COORDS(3*(J2-1)+3,JP))**2)
         IF (DUMMY.LT.SYMTOL5) THEN
            WRITE(MYUNIT,'(A,I8,A,I8,A,G20.10)') 'symmetry2> ERROR - atoms ',J1,' and ',J2,' are separated by only ',DUMMY
      WRITE(41,'(I4)') NATOMS
      WRITE(41,*) ' ' 
      WRITE(41,'(A2,3F20.10)') ('LA ',COORDS(3*(I-1)+1,JP),COORDS(3*(I-1)+2,JP),COORDS(3*(I-1)+3,JP),I=1,NCORE(JP))
      WRITE(41,'(A2,3F20.10)') ('LB',COORDS(3*(I-1)+1,JP),COORDS(3*(I-1)+2,JP),COORDS(3*(I-1)+3,JP),I=NCORE(JP)+1,NATOMS)
            STOP
         ENDIF
      ENDDO
   ENDDO
ENDIF

!   DO J2=1,3*NATOMS
!      X(J2)=COORDS(J2,JP)
!   ENDDO
!   CALL POTENTIAL(X,DUMMYGRAD,DUMMYE,.FALSE.,.FALSE.)
!   DO J2=1,NATOMS
!      IF (VT(J2).NE.0.0D0) THEN
!         IF (ABS((VT(J2)-VAT(J2,JP))/VT(J2)).GT.0.01D0) THEN
!            WRITE(MYUNIT,'(A,I8,2F15.5)') 'symmetry2> A J2,VAT,VT(J2)=',J2,VAT(J2,JP),VT(J2)
!            STOP
!         ENDIF
!      ENDIF
!   ENDDO

!           WRITE(77,'(I4)') NATOMS
!           WRITE(77,*) 'coords at the end of symmetry2'
!           WRITE(77,'(A2,3F20.10)') ('LB',X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS-NCORE(JP))
!           WRITE(77,'(A2,3F20.10)') ('LA ',X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=NATOMS-NCORE(JP)+1,NATOMS)

RETURN

END SUBROUTINE SYMMETRY2
