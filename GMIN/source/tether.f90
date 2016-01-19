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
!
! Tether performs a conventional WL run for selected local minima that have the highest BinImportanceIndex
! with the additional constraint that the distance between the given minimum and any structure that is accepted during the 
! sample is less than the average distance d_jmin for the BS bin to which the minimum belongs. 
! It is possible to split the energy range into a number of windows 'hwindows'. Setting 'hwindows' to 1 
! using a keyword produces a WL
! run of the whole range. (Tetyana Bogdan)
! T.V. Bogdan, D.J. Wales and F. Calvo, J. Chem. Phys., 124, 044102 (2006).
!
!---======================================---
      subroutine TetheredWL
      use Commons, only: Natoms,  coords, histmin, histMAX, hbins, histfac, TargetWL, &
                         & histfacmul, debug, Equil, hdistconstraint, &
                         & DumpEveryNthQuench, hwindows, lhbins, sampledbins, lnHarmFreq, CHRMMT, FixedEndMoveT
      use TetherFunc

      implicit none
    

      integer nsteps, i, j, k, Visits(Hbins), VisitsTotal(Hbins), lbfgs_iterations, BinIndexold, BinIndexnew, NQuenches, &
              & NQuenchesSinceFlat, NQuenchesSuccess, nWL, NQuenchesSinceLastUpdate, EscapedFromTether, NQuenchesMax, &
              & lVisits(lhbins,hwindows), lVisitsTotal(lhbins,hwindows), ndummy, Converged, &
              & lVisits_S(sampledbins,hwindows), lVisitsTotal_S(sampledbins,hwindows) 
      real(8) potel, Weight(Hbins) , Distance(Hbins), dummy, CurrentPointEnergy, CurrentPointCoordinates(3*Natoms, 1), &
              & Distanceold, lnModFac, PerturbedCoordinates(3*Natoms, 1), PerturbedCoordinatesSave(3*Natoms, 1), lnRatio, &
              & harvest, BinLabel(HBins), Norm, MinDistance(Hbins), MinDistanceold, grad(3*Natoms), &
              & SourceMinimum(3*Natoms, 1), HISTINT, &
              & CheckTether, lnharvest, OldEnergy, delta, &
              &  window_coords(3*Natoms, hwindows), & 
              & lBinLabel(HBins/hwindows,hwindows ), l_lnWeight(lhbins,hwindows), tot_lnWeight(hbins-hwindows+1), ScaleFac, &
              & lBinLabel_S(sampledbins,hwindows ), l_lnWeight_S(sampledbins,hwindows), ExtrapolFac, ShiftFac, BottomWeight

!!!! THINK ABOUT THIS:, tot_lnWeight_S(hbins-hwindows+1)


      real(8) lnWeight(Hbins)

      character (len =256)  histfilename
      character (len= 10)  istr

      logical yesno, Flat, evap, evapreject, AcceptMove, allfound

      common /mypot/ potel

      common /ev/ evap, evapreject

      HISTINT=(HISTMAX-HISTMIN)/HBINS
      window_coords=0.0d0
      print *, 'Tethering in ', hwindows, ' windows started. Sampling for coordinates:'
      print *, 'Energy range split as follows: window, Emin, Emax:'
      do j=1, hwindows
          print *, j ! this is a rough not exact estimation of energy ranges
          print *, histmin+(j-1)*histint*(hbins/hwindows), histmin+(j)*histint*(hbins/hwindows)
      enddo


      CurrentPointCoordinates(:,1)=coords(:,1)
      PerturbedCoordinatesSave(:,1)=coords(:,1)
      print *, 'Calculating initial energy'
      call quench(.false.,1,lbfgs_iterations,dummy,ndummy,Converged,CurrentPointCoordinates)
      CurrentPointEnergy=potel
      SourceMinimum(:,1)=coords(:,1) 
      print *, 'Initial energy', CurrentPointEnergy

      do
        allfound=.false.
        do j=1, hwindows
           do i=1, 3*natoms
             if (window_coords(i,j).ne.0.0d0) then
               allfound=.true.
             else
               allfound=.false.
             endif
            enddo
        enddo
        
        if (allfound) exit       
    
        CurrentPointCoordinates(:,1)=coords(:,1)
        
        do i=1, 3*Natoms
             PerturbedCoordinatesSave(i,1)=coords(i,1)
        enddo

        IF (CHRMMT) THEN
            IF (FixedEndMoveT) THEN ! actually basin-hopping-type step in internal coordinates
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
!!!             ELSE
               coords(:3*NATOMS,1) = PerturbedCoordinatesSave(:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
!              CALL FILLICT(1)
               CALL TAKESTEPCH(1)
               ! new geometry is now in coords
               PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
            ELSE ! try taking step in Cartesians
               PerturbedCoordinates = PerturbGeometry(PerturbedCoordinatesSave)
               coords(1:3*NATOMS,1) = PerturbedCoordinates(1:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
            ENDIF
        ELSE
           PerturbedCoordinates=PerturbGeometry(PerturbedCoordinatesSave)
           coords(:,1)=PerturbedCoordinates(:,1)
        ENDIF

   
        !call PrintXyz(Natoms,coords(:,1)) 
        call potential(coords, grad, potel, .true., .false.)
        CurrentPointEnergy=potel
        !call PrintXyz(Natoms,coords(:,1)) 
        CurrentPointCoordinates(:,1)=coords(:,1)

        do j=1, hwindows

           if ((CurrentPointEnergy.gt.histmin+(lhbins-sampledbins)*HistInt+(j-1)*histint*(hbins/hwindows)).and. &
              & (CurrentPointEnergy.lt.histmin+(lhbins-sampledbins)*HistInt+(j)*histint*(hbins/hwindows)).and. & 
              & window_coords(1,j).eq.0.d0) then

              do i=1, 3*Natoms
                 window_coords(i,j)=CurrentPointCoordinates(i,1)
              enddo
              print *, window_coords(1,j), 'e', CurrentPointEnergy
            endif
        enddo

      enddo

      print *, 'Seeding windows with starting geometries successful'

      print *, 'Sampling requested for ', sampledbins, 'bins'


    
!     Global density of states variables: 
      lnWeight=0.0d0
      Visits=0
      VisitsTotal=0
      do i=1, Hbins
         BinLabel(i)=HistMin + HistInt*(i-0.5)
      enddo

      

!****  Variables for each window:

      print *, 'Each window is divided into ', lhbins, 'bins.'

      
   do j=1, hwindows

      print *, 'WL simulation for window' , j


      do i=1, lhbins
         lVisits(i,j)=0
         lVisitsTotal(i,j)=0       
      enddo

      do i=1, sampledbins 
         lVisits_S(i,j)=0
         lVisitsTotal_S(i,j)=0       
      enddo

      do i=1, lhbins
         lBinLabel(i,j)=(HistMin+(j-1)*histint*lhbins) + HistInt*(i-0.5)
      enddo

      do i=1, sampledbins 
         lBinLabel_S(i,j)=(HistMin+(lhbins-sampledbins)*HistInt+(j-1)*histint*lhbins) + HistInt*(i-0.5)
      enddo
      
      if (j.ne.1) then
        do i=1, lhbins
           lBinLabel(i,j)=lBinLabel(i,j)-histint*(j-1)
        enddo
      endif


      !do i=1, sampledbins
      !   print *, lBinLabel_S(i,j), lVisits_S(i,j), l_lnWeight_S(i,j)
      !enddo


 !****    seeding w(x) with harmonic density of states for each window ****

      print *, 'Harmonic Frequency for the Minimum: ', lnHarmFreq
      do i=1, lhbins
         l_lnWeight(i,j)=log((lBinLabel(i,j)-HistMin)**((3*Natoms/2)-1))-lnHarmFreq
      enddo

      l_lnWeight_S=l_lnWeight(lhbins-sampledbins+1:lhbins,:)

!*** OR NOT:
      !l_lnWeight_S=0.d0

       open(unit=422, file='HarmVibDOS.unweighted', status='unknown')
       do i=1, lhbins
          write(422, '(2G20.10)') lBinLabel(i,j), l_lnWeight(i,j)
       enddo
       close(422)

       open(unit=422, file='HarmVibDOS_S.unweighted', status='unknown')
       do i=1, sampledbins 
          write(422, '(2G20.10)') lBinLabel_S(i,j), l_lnWeight_S(i,j)
       enddo
       close(422)

      print *, 'Seeding with starting geometry:'

      do      
        
        if (j==1) then
           IF (CHRMMT) THEN
               IF (FixedEndMoveT) THEN ! actually basin-hopping-type steps in internal coordinates
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
!!!             ELSE
                  coords(:3*NATOMS,1) = window_coords(:,j)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
                  CALL TAKESTEPCH(1)
                  ! new geometry is now in coords
                  PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
               ELSE ! try taking step in Cartesians
                  PerturbedCoordinates=PerturbGeometry(window_coords(:,j))
                  coords(1:3*NATOMS,1) = PerturbedCoordinates(1:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
               ENDIF
           ELSE
              PerturbedCoordinates=PerturbGeometry(window_coords(:,j))
              coords(:,1)=PerturbedCoordinates(:,1)
           ENDIF

        else
           IF (CHRMMT) THEN
               IF (FixedEndMoveT) THEN ! actually basin-hopping-type steps in internal coordinates
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
!!!             ELSE
                  coords(:3*NATOMS,1) = window_coords(:,(j-1))
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
                  CALL TAKESTEPCH(1)
                  ! new geometry is now in coords
                  PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
               ELSE ! try taking step in Cartesians
                  PerturbedCoordinates=PerturbGeometry(window_coords(:,(j-1)))
                  coords(1:3*NATOMS,1) = PerturbedCoordinates(1:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
               ENDIF
           ELSE
              PerturbedCoordinates=PerturbGeometry(window_coords(:,(j-1)))
              coords(:,1)=PerturbedCoordinates(:,1)
           ENDIF
        endif

        call potential(coords, grad, potel, .true., .false.)
        CurrentPointEnergy=potel
        CurrentPointCoordinates(:,1)=coords(:,1)
        BinIndexold=Energy2Index(CurrentPointEnergy, lBinLabel_S(1,j))
        print *, CurrentPointEnergy, BinIndexold, lBinLabel_S(1,j), lBinLabel_S(sampledbins,j)
        if ((BinIndexold.ge.1 ).and.(BinIndexold.le.sampledbins)) then
          do i=1, 3*Natoms
             window_coords(i,j)=CurrentPointCoordinates(i,1)
          enddo
          print *, window_coords(1,j), 'e', CurrentPointEnergy
          exit
        else
          do i=1, 3*Natoms
             window_coords(i,(j-1))=CurrentPointCoordinates(i,1)
          enddo
        endif

      enddo


      lnModFac=log(HistFac)
      nWL=0

      do i=1, 3*Natoms
         CurrentPointCoordinates(i,1)=window_coords(i,j)
      enddo

      do i=1, 3*Natoms
         PerturbedCoordinatesSave(i,1)=window_coords(i,j)
      enddo
      coords=CurrentPointCoordinates
      call potential(coords, grad, potel, .true., .false.)
      !call quench(.false.,1,lbfgs_iterations,dummy,dummy,Converged,coords)
      CurrentPointEnergy=potel
      OldEnergy=potel
      print *, 'Initial energy', CurrentPointEnergy, 'for window', j



      BinIndexold=Energy2Index(CurrentPointEnergy, lBinLabel_S(1,j))

      !if (debug) print *, CurrentPointEnergy, BinIndexold, 'here'

      if ((BinIndexold < 1 ).or.(BinIndexold > sampledbins)) then
         print *, 'Starting geometry is outside requested range. Exiting.'
         return
      endif 

      NQuenches=0
      NQuenchesSinceFlat=0
      NQuenchesSuccess=0
      EscapedFromTether=0
      NQuenchesSinceLastUpdate=0


!Repeat for the requested number of Wang-Landau iterations. 


      do
         if ( nWL==TargetWL )  then
          print *, 'WL run for window', j, 'complete.'
          exit
         endif

         IF (CHRMMT) THEN
             IF (FixedEndMoveT) THEN ! actually basin-hopping-type step in internal coordinates
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
!!!             ELSE
                coords(:3*NATOMS,1) = PerturbedCoordinatesSave(:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
                CALL TAKESTEPCH(1)
                ! new geometry is now in coords
                PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
             ELSE ! try taking step in Cartesians
                PerturbedCoordinates = PerturbGeometry(PerturbedCoordinatesSave)
                coords(1:3*NATOMS,1) = PerturbedCoordinates(1:3*NATOMS,1)
                WRITE(*,'(A)') 'tether> please fix the call to FILLICT to compile this branch of GMIN'
                STOP 
!               CALL FILLICT(1)
             ENDIF
         ELSE
            PerturbedCoordinates=PerturbGeometry(PerturbedCoordinatesSave)
            coords(:,1)=PerturbedCoordinates(:,1)
         ENDIF

         CheckTether=CalculatedDistance(SourceMinimum, PerturbedCoordinates)
         !print *, 'Distance2',CheckTether,hdistconstraint
         call potential(coords, grad, potel, .true., .false.)
         !print *, 'Distance3', CheckTether, hdistconstraint
         if ((evapreject).or.(CheckTether.gt.hdistconstraint)) then
            EscapedFromTether=EscapedFromTether+1
            cycle
         endif
         oldenergy=CurrentPointEnergy
         CurrentPointEnergy=potel
         CurrentPointCoordinates=coords


         if (.not.evapreject) then
            NQuenchesSuccess=NQuenchesSuccess+1
            BinIndexnew=Energy2Index(CurrentPointEnergy, lBinLabel_S(1,j))


            if (BinIndexnew < 1 .or.BinIndexnew>sampledbins) then
               !print *, 'Structure outside energy range. Rejecting move'
                AcceptMove=.false. 
            else 
               lnRatio=l_lnWeight_S(BinIndexold,j)-l_lnWeight_S(BinIndexnew,j)
               call random_number(harvest)
               lnharvest=log(harvest)

                  if ((lnRatio>0.d0).or.(lnRatio>lnharvest)) then
                     AcceptMove=.true.
                  else
                     AcceptMove=.false.
                  endif
            endif

            !print *, 'Eold, Enew, Iold, Inew, Wold, Wnew, AcceptMove , Ratio'
            !print *, oldenergy, CurrentPointEnergy, BinIndexold, BinIndexnew, l_lnWeight(BinIndexold,j) , &
            !         & l_lnWeight(BinIndexnew,j), AcceptMove, lnRatio


            if (AcceptMove) then 
               ! move accepted
                 lVisits_S(BinIndexnew,j)=lVisits_S(BinIndexnew,j)+1
                 lVisitsTotal_S(BinIndexnew,j)=lVisitsTotal_S(BinIndexnew,j)+1
                 l_lnWeight_S(BinIndexnew,j)=l_lnWeight_S(BinIndexnew,j)+lnModFac
               PerturbedCoordinatesSave=PerturbedCoordinates
               BinIndexold=BinIndexnew
               OldEnergy=CurrentPointEnergy
            else
               ! move rejected
               !print *, 'Staying in bin', BinIndexold
                  lVisits_S(BinIndexold,j)=lVisits_S(BinIndexold,j)+1
                  lVisitsTotal_S(BinIndexold,j)=lVisitsTotal_S(BinIndexold,j)+1
                  l_lnWeight_S(BinIndexold,j)=l_lnWeight_S(BinIndexold,j)+lnModFac
            endif
         else
         if (debug) print *, 'Optimisation unsuccessful'
         endif

         if ((mod(NQuenchesSuccess,DumpEveryNthQuench).eq.0).and.(NQuenchesSinceLastUpdate > Equil)) then 
            ! record statistics
            ! print *, NQuenchesSuccess
             write (istr, '(i10)') j 

             histfilename="lnWeight.his."//trim(adjustl(istr))
             open(unit=421, file=histfilename, status='unknown')
             do i=1, sampledbins
                write(421, '(2G20.10)')  lBinLabel_S(i,j), l_lnWeight_S(i,j)
             enddo
             close(421)
           
             histfilename="Visits.his."//trim(adjustl(istr))
             open(unit=422, file=histfilename, status='unknown')
             do i=1, sampledbins
                write(422, '(2G20.10)')  lBinLabel_S(i,j), lVisits_S(i,j)
             enddo
             close(422)
      
             histfilename="VisitsTotal.his."//trim(adjustl(istr))
             open(unit=422, file=histfilename, status='unknown')
             do i=1, sampledbins
                write(422, '(2G20.10)')  lBinLabel_S(i,j), lVisitsTotal_S(i,j)
             enddo
             close(422)
           
             open(unit=422, file='nWL.restart', status='unknown')
             write(422, '(4G20.10)') nWL, nint(1.0d0/lnModFac), j, NQuenches
             close(422)
      
             open(unit=422, file='EscapedFromTether', status='unknown')
             write(422, '(G20.10)') EscapedFromTether
             close(422)
 

             Flat=CheckFlatness(lVisits_S(:,j), lnModFac, nWL) 
             if (Flat) then
                NQuenchesSinceLastUpdate=0
                print *, '--===Visits histogram satisfied flatness criterion===--'
                lnModFac=HistFacMul*lnModFac
                print *, 'Updating modfac to', lnModFac, 'next WL iteration:', nWL+1, 'Visits needed:', nint(1.0d0/lnModFac)
                do i=1, sampledbins
                   lVisits_S(i,j)=0
                enddo
                nWL=nWL+1
                endif
      
         endif

      NQuenches=NQuenches+1
      NQuenchesSinceLastUpdate=NQuenchesSinceLastUpdate+1
      if (NQuenches.gt.200000) exit
      enddo

    enddo


    print *, 'Vibrational density of states generation in windows complete. Starting scaling:'

    tot_lnWeight=0.d0


    if (hwindows.gt.1) then
       ScaleFac=1.0d0
       do j=1, hwindows
           if (j.eq.1) then
              do i=1, lhbins
                 lnWeight(i+(j-1)*lhbins)=l_lnWeight(i, j)
              enddo
           else
              ScaleFac=l_lnWeight(lhbins, j-1)/l_lnWeight(1, j)
              print *, l_lnWeight(lhbins, j-1), l_lnWeight(1, j), ScaleFac, l_lnWeight(1, j)*ScaleFac 
              do i=1, lhbins
                 l_lnWeight(i, j)=l_lnWeight(i, j)*ScaleFac
              enddo
              write (istr, '(i10)') j 
              histfilename="lnWeight2.his."//trim(adjustl(istr))
              open(unit=421, file=histfilename, status='unknown')
              do i=1, lHbins
                 write(421, '(2G20.10)')  lBinLabel(i,j), l_lnWeight(i,j)
              enddo
              close(421)
              do i=1, lhbins
                 lnWeight(i+(j-1)*lhbins)=l_lnWeight(i, j)
              enddo
           endif
        enddo
      else
       do j=1, hwindows

         print *, lBinLabel(lhbins-sampledbins+1,j), lBinLabel_S(1,j) 
         print *, l_lnWeight(lhbins-sampledbins+1,j), l_lnWeight_S(1,j) 

         ExtrapolFac=abs(l_lnWeight_S(1,j)-l_lnWeight((lhbins-sampledbins+1),j))
         print *, 'ExtrapolFac', ExtrapolFac 

        do i=1, sampledbins
           l_lnWeight_S(i,j)=l_lnWeight_S(i,j)-ExtrapolFac
        enddo
        l_lnWeight(lhbins-sampledbins+1:lhbins,:)=l_lnWeight_S

!        ShiftFac=(log((lBinLabel(i,j)-HistMin)**((3*Natoms/2)-1))-lnHarmFreq)-l_lnWeight(1,j)
!        print *, ShiftFac
!
!        do i=1, lhbins
!              l_lnWeight(i,j)=l_lnWeight(i,j)+ShiftFac
!        enddo

    
        lnWeight=l_lnWeight(:,j) 
       enddo
      endif

           

    histfilename="lnWeight.his"
    open(unit=422, file=histfilename, status='unknown')
    !do i=1, hbins-hwindows+1
    do i=1, hbins
        write(422, '(2G20.10)')  BinLabel(i), lnWeight(i)
    enddo
    close(422)
    end subroutine TetheredWL

