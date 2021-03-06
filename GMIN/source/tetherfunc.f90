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


! Functions module to serve tethered WL subroutine (Tetyana Bogdan)
!---======================================---
      module tetherfunc 
         implicit none
      contains

!---======================================---`
      function Energy2Index(CurrentPointEnergy, BinLabelBottom)
      use Commons
      implicit none

      integer Energy2Index 

      real(8) CurrentPointEnergy, BinLabelBottom, HISTINT

      HISTINT=(HISTMAX-HISTMIN)/HBINS
      if (nint((CurrentPointEnergy-(BinLabelBottom-histint/2.0d0))/HistInt)<0) then
         Energy2Index=-1
      else
         Energy2Index=nint((CurrentPointEnergy-(BinLabelBottom-histint/2.0d0))/HistInt)+1
      endif
!     WRITE(MYUNIT,'(A,3G20.10,I6)') 'Energy2Index> CurrentPointEnergy,BinLabelBottom,histint,Energy2Index=', &
!    &                                              CurrentPointEnergy,BinLabelBottom,histint,Energy2Index

      end function Energy2Index


!---======================================---`
      real(8) function GetDisplacement(step)
           implicit none
           real(8),intent(in) :: step
           real(8) :: harvest
           call random_number(harvest)
           GetDisplacement=2.0d0*harvest*step-step
      end function GetDisplacement

!---======================================---`
      function PerturbGeometry(CurrentPointCoordinates)
           use Commons, only: Radius,Natoms,step
           implicit none

           real(8),dimension(3*Natoms,1) :: PerturbGeometry

           real(8),intent(in) :: CurrentPointCoordinates(3*Natoms,1)
           !real(8) GetDisplacement
           integer i, j

           PerturbGeometry=CurrentPointCoordinates

           do j=1, Natoms
                do
                 print *, 'pert'
                     PerturbGeometry(3*(j-1)+1,1) = PerturbGeometry(3*(j-1)+1,1) + GetDisplacement(step(1))
                     PerturbGeometry(3*(j-1)+2,1) = PerturbGeometry(3*(j-1)+2,1) + GetDisplacement(step(1))
                     PerturbGeometry(3*(j-1)+3,1) = PerturbGeometry(3*(j-1)+3,1) + GetDisplacement(step(1))
                     if ( PerturbGeometry(3*(j-1)+1, 1)**2 &
                      & + PerturbGeometry(3*(j-1)+2, 1)**2 &
                      & + PerturbGeometry(3*(j-1)+3, 1)**2 < Radius ) then
                          exit
                     else
                          PerturbGeometry(3*(j-1)+1,1) = CurrentPointCoordinates(3*(j-1)+1,1)
                          PerturbGeometry(3*(j-1)+2,1) = CurrentPointCoordinates(3*(j-1)+2,1)
                          PerturbGeometry(3*(j-1)+3,1) = CurrentPointCoordinates(3*(j-1)+3,1)
                     endif
		enddo
           enddo
      end function PerturbGeometry

!---======================================---`

      function CalculatedDistance(CurrentPointCoordinates, PerturbedCoordinates)
      use Commons
      implicit none

      integer i
      real(8) CalculatedDistance, CurrentPointCoordinates(3*Natoms,1), PerturbedCoordinates(3*Natoms,1), &
              & x1(Natoms, 1), y1(Natoms, 1), z1(Natoms, 1), x2(Natoms, 1), y2(Natoms, 1), z2(Natoms, 1)


      do i=1,Natoms
         x1(i,1)=CurrentPointCoordinates(3*(i-1)+1,1)
         y1(i,1)=CurrentPointCoordinates(3*(i-1)+2,1)
         z1(i,1)=CurrentPointCoordinates(3*(i-1)+3,1)
         x2(i,1)=PerturbedCoordinates(3*(i-1)+1,1)
         y2(i,1)=PerturbedCoordinates(3*(i-1)+2,1)
         z2(i,1)=PerturbedCoordinates(3*(i-1)+3,1)
      enddo

      CalculatedDistance=dsqrt(sum ( (x1(:,1)-x2(:,1))**2 + (y1(:,1)-y2(:,1))**2 + (z1(:,1)-z2(:,1))**2 ) )

      !print *, 'Distance' , CalculatedDistance

      !if (CalculatedDistance>0.8) then
      !    call PrintXyz(Natoms,PerturbedCoordinates(:,1))
      !    call PrintXyz(Natoms,CurrentPointCoordinates(:,1))
      !endif
      end function CalculatedDistance

!---======================================---`
      subroutine PrintXyz(Natoms,X)
          implicit none

          integer,intent(in) :: Natoms
          real(8),intent(in) :: X(3*Natoms)

          integer :: i

          print *, Natoms
          print *
          do i=1,Natoms
             write(*,'(a4,3f20.10)') 'C   ',X(3*(i-1)+1),X(3*(i-1)+2),X(3*(i-1)+3)
          enddo
      end subroutine PrintXyz
!---======================================---`
      function CheckFlatness(lVisits_S,lnModFac,nWL)
      use Commons, only: lhbins, sampledbins
      implicit none

      integer lVisits(lhbins), n, i, VisitsNonzero(lhbins),Hmin, Hmax, deltaHk, nWL, lVisits_S(sampledbins)
      logical CheckFlatness, FlatBin(lhbins), FlatBin_S(sampledbins)
      real(8) HDev, lnModFac
      character (len =256) filename
      character (len= 10)  istr


      CheckFlatness=.false.


!      Hmax=maxval(lVisits)
!      do i=1, lHbins
!         if (lVisits(i).ne.0) then
!            VisitsNonzero(i)=lVisits(i)
!         else
!            VisitsNonzero(i)=huge(lhbins)
!         endif
!      enddo
!      Hmin=minval(VisitsNonzero)
!      if (debug) print *, 'Min and max', Hmin, Hmax
!
!      HDev=(1.0d0*Hmax-Hmin)/(Hmax+Hmin)
    
!  
!     (abs(Hmin-Hmax)>HPercent)) needed to prevent the run from exiting right after the equilibration

!     FLATNESS CONDITION
!     if ((HDev<HPercent).and.(abs(1.0d0*Hmin-Hmax)>HPercent)) CheckFlatness=.true.


!     VisitProp:minimal number of visits should be proportional to 1/sqrt(ln(f))

         CheckFlatness=.false.
         FlatBin_S=.false.
 
         do i=1, sampledbins
            if (lVisits_S(i) > 1.0d0/sqrt(lnModFac)) FlatBin_S(i)=.True.
         enddo
 
         if (All(FlatBin_S)) CheckFlatness=.true.

      end function CheckFlatness

!---======================================---
    SUBROUTINE SAVEBINSTRUCTURES(CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES, BININDEX, MINIMANUMBER, WRITESTRUCT)
        USE COMMONS, ONLY:NATOMS, HBINS, CHRMMT, ZSYM, ECONV
        USE MODCHARMM, ONLY: OSASAT, RPRO, ODIHET
        IMPLICIT NONE

        INTEGER BININDEX, JB, JC, MINIMANUMBER(HBINS)
        INTEGER, SAVE :: BINENSAVED=0
        REAL(8) CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES(3*NATOMS, 1), Q4ORDERPARAM
        REAL(8) DIHEORDERPARAM,SASAORDERPARAM
        LOGICAL NEWBINENERGY, YESNO, WRITESTRUCT
        CHARACTER (LEN =256)  BINNAME
        CHARACTER (LEN= 10)  ISTR
        REAL(8),ALLOCATABLE, SAVE :: BINENERGIES(:)
        INTEGER,ALLOCATABLE, SAVE :: BINENIMPORTANCEINDEX(:)
        
         IF (BINENSAVED.EQ.0) THEN
           ALLOCATE(BINENERGIES(100000000))
           BINENERGIES=0.0
           ALLOCATE(BINENIMPORTANCEINDEX(100000000))
           BINENIMPORTANCEINDEX=0
         ENDIF
      
         IF (BINENSAVED.EQ.SIZE(BINENERGIES)) THEN
            RETURN
         ENDIF

         INQUIRE(FILE='binenergies', exist=yesno)

         NEWBINENERGY=.TRUE.
         DO JB=1, BINENSAVED
!            IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.BQMAX) THEN
             IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.ECONV) THEN
                NEWBINENERGY=.FALSE.
             ENDIF
         ENDDO

         IF (NEWBINENERGY) THEN
            BINENSAVED=BINENSAVED+1
            !call ORDERQ4(Natoms, CurrentPointCoordinates, Q4orderparam)
            IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,CURRENTPOINTCOORDINATES(1:3*NATOMS,1))
            IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,CURRENTPOINTCOORDINATES(1:3*NATOMS:3,1), &
                      & CURRENTPOINTCOORDINATES(2:3*NATOMS:3,1),CURRENTPOINTCOORDINATES(3:3*NATOMS:3,1))
            BINENERGIES(BINENSAVED)=CURRENTPOINTENERGY
            IF (WRITESTRUCT)  THEN
                IF ((BINENSAVED.EQ.1).AND.(YESNO)) THEN
                   OPEN(UNIT=1979,FILE='binenergies',status='old')
                ELSE
                   OPEN(UNIT=1979,FILE='binenergies',status='unknown',position='append')
                ENDIF
! jmc                write(1979,'(2G20.10)') CurrentPointEnergy, Q4orderparam
                IF (ODIHET.AND.OSASAT) THEN
                   WRITE(1979,'(3G20.10)') CurrentPointEnergy, diheorderparam, SASAorderparam
                ELSEIF (OSASAT) THEN
                   WRITE(1979,'(2G20.10)') CurrentPointEnergy, SASAorderparam
                ELSEIF (ODIHET) THEN
                   WRITE(1979,'(2G20.10)') CurrentPointEnergy, diheorderparam
                ELSE
                   WRITE(1979,'(G20.10)') CurrentPointEnergy
                ENDIF
                CLOSE (1979)
                WRITE (ISTR, '(i10)') binindex
                BINNAME="binstructures."//trim(adjustl(istr))
                IF ((BINENSAVED.EQ.1).AND.(YESNO)) THEN
                   OPEN(UNIT=1979,FILE=BINNAME, STATUS="old", form="formatted")
                ELSE
                   OPEN(UNIT=1979,FILE=BINNAME, STATUS="unknown", form="formatted", position="append")
                ENDIF
                WRITE(1979,'(I10)') natoms
                WRITE(1979,'(A,1G20.10)') '     Structure with energy ', CurrentPointEnergy
                IF (CHRMMT) THEN
                   DO JB=1,NATOMS
                      WRITE(1979,'(A,1X,3F20.10)') ZSYM(JB)(1:1),(CurrentPointCoordinates(3*(JB-1)+JC,1),JC=1,3)
                   ENDDO
                ELSE
                   WRITE(1979,30) (CURRENTPOINTCOORDINATES(JB,1),JB=1,3*NATOMS)
30                 FORMAT('LA ',3F20.10)
                ENDIF
                CLOSE (1979)
             ENDIF
             MINIMANUMBER(BININDEX)=MINIMANUMBER(BININDEX)+1
         ENDIF

! Importance Index 

 	 IF ((.NOT.NEWBINENERGY).AND.(WRITESTRUCT)) THEN 
		OPEN(UNIT=1979,FILE="BinEnImportanceIndex", status="unknown", form="formatted")
		DO JB=1, BINENSAVED
!                  IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.BQMAX) THEN
	           IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.ECONV) THEN
	              BINENIMPORTANCEINDEX(JB)=BINENIMPORTANCEINDEX(JB)+1
             	   ENDIF
         	ENDDO
		DO JB=1, BINENSAVED
		WRITE(1979,'(1G20.10, I10)') binenergies(jb), binenImportanceIndex(jb)
		ENDDO
		CLOSE (1979)
	 ENDIF
 
         IF (BINENSAVED.EQ.SIZE(BINENERGIES)) THEN
            IF (ALLOCATED(BINENERGIES)) DEALLOCATE(BINENERGIES)
         ENDIF

         END SUBROUTINE SAVEBINSTRUCTURES

!---======================================---
    SUBROUTINE SAVEBINSTRUCTURESMPI(CURRENTPOINTENERGY,CURRENTPOINTCOORDINATES,BININDEX,WRITESTRUCT,NODE,NEWBINENERGY,MYUNIT)

        USE COMMONS, ONLY : NATOMS, HBINS, CHRMMT, ZSYM, NPAR, ECONV
        USE MODCHARMM, ONLY : OSASAT, RPRO, ODIHET
        IMPLICIT NONE

        INTEGER BININDEX, JB, JC, NODE, I, MINIMANUMBER(HBINS,0:NPAR-1), MYUNIT, GETUNIT
        INTEGER, SAVE :: BINENSAVED=0
        INTEGER, SAVE :: BINENMAX=1000
        INTEGER, SAVE :: LUNIT
        REAL(8) CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES(3*NATOMS, 1), Q4ORDERPARAM
        REAL(8) DIHEORDERPARAM,SASAORDERPARAM
        LOGICAL YESNO, WRITESTRUCT, NEWBINENERGY
        CHARACTER (LEN=256)  BINNAME, FILENAME, BINNAME2
        CHARACTER (LEN=10)  ISTR, JSTR

        REAL(8),ALLOCATABLE, SAVE :: BINENERGIES(:), ETEMP(:)
!       INTEGER,ALLOCATABLE, SAVE :: BINENIMPORTANCEINDEX(:,:)
        
         IF (BINENSAVED.EQ.0) THEN
           ALLOCATE(BINENERGIES(BINENMAX))
           BINENERGIES=0.0D0
           LUNIT=GETUNIT()
           WRITE (JSTR,'(I10)') NODE+1
           FILENAME=TRIM(ADJUSTL(JSTR)) // '/quench_structures.xyz'
           OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='UNKNOWN')
!          ALLOCATE(BINENIMPORTANCEINDEX(10000,0:NPAR-1))
!          BINENIMPORTANCEINDEX(:,NODE)=0
         ENDIF
      
!        WRITE (ISTR, '(i10)') NODE+1
!        FILENAME=trim(adjustl(istr)) // "/binenergies"
!        INQUIRE(FILE='binenergies', exist=yesno)

         NEWBINENERGY=.TRUE.
         DO JB=1, BINENSAVED
             IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.ECONV) THEN
                NEWBINENERGY=.FALSE.
                RETURN
             ENDIF
         ENDDO

         BINENSAVED=BINENSAVED+1
         IF (BINENSAVED.GT.BINENMAX) THEN
            ALLOCATE(ETEMP(BINENMAX))
            ETEMP(1:BINENMAX)=BINENERGIES(1:BINENMAX)
            DEALLOCATE(BINENERGIES)
            ALLOCATE(BINENERGIES(2*BINENMAX))
            BINENERGIES(1:BINENMAX)=ETEMP(1:BINENMAX)
            BINENMAX=2*BINENMAX
            DEALLOCATE(ETEMP)
         ENDIF

         !call ORDERQ4(Natoms, CurrentPointCoordinates, Q4orderparam)
         IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,CURRENTPOINTCOORDINATES(1:3*NATOMS,1))
         IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,CURRENTPOINTCOORDINATES(1:3*NATOMS:3,1), &
                   & CURRENTPOINTCOORDINATES(2:3*NATOMS:3,1),CURRENTPOINTCOORDINATES(3:3*NATOMS:3,1))
         BINENERGIES(BINENSAVED)=CURRENTPOINTENERGY
         IF (WRITESTRUCT)  THEN
             WRITE(LUNIT,'(I10)') NATOMS
             WRITE(LUNIT,'(A,G20.10)') 'energy ', CURRENTPOINTENERGY
             IF (CHRMMT) THEN
                DO JB=1,NATOMS
                   WRITE(LUNIT,'(A,1X,3F20.10)') ZSYM(JB)(1:1),(CURRENTPOINTCOORDINATES(3*(JB-1)+JC,1),JC=1,3)
                ENDDO
             ELSE
                WRITE(LUNIT,30) (CURRENTPOINTCOORDINATES(JB,1),JB=1,3*NATOMS)
30              FORMAT('LA ',3F20.10)
             ENDIF
         ENDIF
!
! Not currently used.
!
!        MINIMANUMBER(BININDEX, NODE)=MINIMANUMBER(BININDEX, NODE)+1
!
! Importance Index 
!
!	 IF ((.NOT.NEWBINENERGY).AND.(WRITESTRUCT)) THEN 
!           BINNAME=trim(adjustl(istr)) // "/BinEnImportanceIndex"
!           OPEN(UNIT=1980+node,FILE=binname, status="unknown", form="formatted")
!           DO JB=1, BINENSAVED
!              IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.ECONV) THEN
!	          BINENIMPORTANCEINDEX(JB,node)=BINENIMPORTANCEINDEX(JB,node)+1
!              ENDIF
!           ENDDO
!	    DO JB=1, BINENSAVED
!	       WRITE(1980+node,'(1G20.10, I10)') binenergies(jb), binenImportanceIndex(jb,node)
!	    ENDDO
!	    CLOSE (1980+node)
!        ENDIF
 
         END SUBROUTINE SAVEBINSTRUCTURESMPI
   
    end module tetherfunc
