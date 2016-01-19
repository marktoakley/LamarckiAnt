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
! The `basin-sampling' algorithm combines `basin-hopping' and Wang-Landau
! sampling technique to study the thermodynamics of the transformed PES.
! It provides a direct temperature-independent estimate of the density of states
! along with thermodynamic properties such as the free energy and entropy via ensemble averages using
! samples of local minima rather than instantaneous configurations. (Tetyana Bogdan)
! T.V. Bogdan, D.J. Wales and F. Calvo, J. Chem. Phys., 124, 044102 (2006).
!                                            
!---======================================---
      SUBROUTINE BASINSAMPLING

      USE MODCHARMM
      USE COMMONS, ONLY: NATOMS,  COORDS, BSPTRESTART, HISTMIN, HISTMAX, HBINS, HISTFAC, TARGETWL, &
                         & HISTFACMUL, DEBUG, EQUIL , PERIODIC, TWOD, BINSTRUCTURES, SAVENTH, DUMPEVERYNTHQUENCH, &
                         & FIXEDENDMOVET, CHRMMT, RIGID, FIXCOM
      IMPLICIT NONE
    

      INTEGER NSTEPS, I, J, VISITS(HBINS), VISITSTOTAL(HBINS), LBFGS_ITERATIONS, BININDEXOLD, BININDEXNEW, NQUENCHES, &
              & NQUENCHESSINCEFLAT, NQUENCHESSUCCESS, NWL, NQUENCHESSINCELASTUPDATE, MINIMANUMBER(HBINS), NDUMMY, CONVERGED, &
              & FROZENVISITS(HBINS)
      REAL(8) POTEL, DUMMY, CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES(3*NATOMS, 1), &
              & DISTANCEOLD, LNMODFAC, PERTURBEDCOORDINATES(3*NATOMS, 1), PERTURBEDCOORDINATESSAVE(3*NATOMS, 1), LNRATIO, &
              & HARVEST, BINLABEL(HBINS), OLDENERGY, NORM, MINDISTANCE(HBINS), MINDISTANCEOLD, LNHARVEST, QWEIGHT(HBINS,2), &
              & CURRENTQ, Q4ORDERPARAM, KNOWNENERGIES(61836), RMAT(3,3), HISTINT

      !real(16) lnWeight(Hbins), Distance(Hbins) 
      REAL(8) LNWEIGHT(HBINS), DISTANCE(HBINS) 
         
      LOGICAL YESNO, FLAT, EVAP, EVAPREJECT, ACCEPTMOVE

      COMMON /MYPOT/ POTEL

      COMMON /EV/ EVAP, EVAPREJECT

! If the BS run is requested to be restarted the following files will be read in:
      
      IF (BSPTRESTART) THEN
         OPEN(UNIT=421, FILE='lnWeight.his', status='old')
         DO I=1, HBINS
            READ(421, '(2G20.10)') dummy, lnWeight(i)
         ENDDO
         PRINT *, 'Following ln(Weight) histogram was read:'
         PRINT *, LNWEIGHT
         CLOSE(421)
         OPEN(UNIT=421, FILE='Distance.his', status='old')
         DO I=1, HBINS
            READ(421, '(2G20.10)') dummy, Distance(i)
         ENDDO
         PRINT *, 'Following Distance histogram was read:'
         PRINT *, DISTANCE
         CLOSE(421)
         OPEN(UNIT=421, FILE='MinDistance.his', status='old')
         DO I=1, HBINS
            READ(421, '(2G20.10)') dummy, MinDistance(i)
         ENDDO
         PRINT *, 'Following Minimized Distance histogram was read:'
         PRINT *, MINDISTANCE
         CLOSE(421)
      ELSE
! note that weight accumulation is is done in logarithms
         LNWEIGHT=0.0D0 
         QWEIGHT=0.0D0
         DISTANCE=0.0D0
         MINDISTANCE=0.0D0
      ENDIF

      IF (BSPTRESTART) THEN
         OPEN(UNIT=422, FILE='Visits.his', status='old')
         DO I=1, HBINS
            READ(422, '(2G20.10)') dummy, Visits(i)
         ENDDO
         PRINT *, 'Following Visits histogram was read:'
         PRINT *, VISITS
         CLOSE(422)
         OPEN(UNIT=422, FILE='VisitsTotal.his', status='old')
         DO I=1, HBINS
            READ(422, '(2G20.10)') dummy, VisitsTotal(i)
         ENDDO
         PRINT *, 'Following VisitsTotal histogram was read:'
         PRINT *, VISITSTOTAL
         CLOSE(422)
         OPEN(UNIT=422, FILE='lnModfac.restart', status='old')
         READ(422, '(G20.10)') lnModFac
         PRINT *, 'Restarting from modification factor', lnModFac
         CLOSE(422)
         OPEN(UNIT=422, FILE='nWL.restart', status='old')
         READ(422, '(G20.10)') nWL
         PRINT *, 'Number of Wang Landau iterations already completed: ', nWL
         CLOSE(422)
      ELSE
         VISITS=0
         VISITSTOTAL=0
         MINIMANUMBER=0
         LNMODFAC=LOG(HISTFAC) ! REFER TO THE TABLE
         PRINT *, 'lnModfac, exp', lnModFac, exp(lnModFac)
         NWL=0
      ENDIF

      HISTINT=(HISTMAX-HISTMIN)/HBINS
      DO I=1, HBINS
         BINLABEL(I)=HISTMIN + HISTINT*(I-0.5)
      ENDDO


! for VISITPROP convergence scenario 
! if one wants to use previously accumulated number of visits, set Equil to -1 in data and copy 
! old VisitsTotal.his to Frozen.his
!
      IF (EQUIL.EQ.-1) THEN
         OPEN(UNIT=422, FILE='Frozen.his', status='old')
         DO I=1, HBINS
            READ(422, '(2G20.10)') dummy, FrozenVisits(i)
         ENDDO
       ENDIF

      
      IF (FIXCOM) CALL CENTRECOM(COORDS(1:3*NATOMS,1))
      CURRENTPOINTCOORDINATES=COORDS
      PERTURBEDCOORDINATESSAVE=COORDS
      PRINT *, 'Calculating initial energy'
      CALL QUENCH(.FALSE.,1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,CURRENTPOINTCOORDINATES)
      CURRENTPOINTENERGY=POTEL
      PRINT *, 'Initial energy', CurrentPointEnergy


      BININDEXOLD=ENERGY2INDEX(CURRENTPOINTENERGY)

      IF ((BININDEXOLD < 1 ).OR.(BININDEXOLD > HBINS)) THEN
         PRINT *, 'Starting geometry is outside requested range. Exiting.'
         RETURN
      ENDIF 
      DISTANCEOLD=0.0D0       

      NQUENCHES=0
      NQUENCHESSINCEFLAT=0
      NQUENCHESSUCCESS=0
      CONVERGED=0
      NQUENCHESSINCELASTUPDATE=0

! Repeat for the requested number of Wang-Landau iterations. 
! The iteration is complete when the visits histogram satisfies the flatness criterion.

      DO
         IF ( NWL==TARGETWL ) EXIT
! jmc do something different here if charmm
         IF (CHRMMT) THEN

! Fixed end move scheme
             IF (FIXEDENDMOVET) THEN ! STEALING FEM KEYWORD FOR BASIN-HOPPING GMIN IC MOVES...
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(1:3*NATOMS,1) = coords(1:3*NATOMS,1)
!!!             ELSE
! step in internal coordinates
                COORDS(1:3*NATOMS,1) = PERTURBEDCOORDINATESSAVE(1:3*NATOMS,1)
                IF(CHRIGIDTRANST) CALL MKRIGIDTRANS(1)
                IF(CHRIGIDROTT) CALL MKRIGIDROT(1)
                CALL TAKESTEPCH(1)
                ! new geometry is now in coords
                IF (FIXCOM) CALL CENTRECOM(COORDS(1:3*NATOMS,1))
                PERTURBEDCOORDINATES(1:3*NATOMS,1) = COORDS(1:3*NATOMS,1)
             ELSE
! try taking step in Cartesians; need to preserve detailed balance and explore all of configuration space.
                PERTURBEDCOORDINATES = PERTURBGEOMETRY(PERTURBEDCOORDINATESSAVE)
                IF (FIXCOM) CALL CENTRECOM(PERTURBEDCOORDINATES(1:3*NATOMS,1))
                COORDS(1:3*NATOMS,1) = PERTURBEDCOORDINATES(1:3*NATOMS,1)
                WRITE(*,'(A)') 'basinsampling> please fix the call to FILLICT to compile this branch of GMIN'
                STOP
!               CALL FILLICT(1)
             ENDIF
         ELSE
            PERTURBEDCOORDINATES=PERTURBGEOMETRY(PERTURBEDCOORDINATESSAVE) 
            IF (FIXCOM) CALL CENTRECOM(PERTURBEDCOORDINATES(1:3*NATOMS,1))
            COORDS=PERTURBEDCOORDINATES
         ENDIF

         CALL QUENCH(.FALSE.,1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS)
         IF (EVAPREJECT) THEN
            IF (DEBUG) PRINT *, 'Evaporation during minimisation'
            IF (DEBUG) PRINT *, 'oldenergy, CurrentPointEnergy, potel=',oldenergy, CurrentPointEnergy, potel
            CYCLE
         ENDIF
         OLDENERGY=CURRENTPOINTENERGY
         CURRENTPOINTENERGY=POTEL
         CURRENTPOINTCOORDINATES=COORDS


         IF ((CONVERGED.EQ.1).AND.(.NOT.EVAPREJECT)) THEN
            NQUENCHESSUCCESS=NQUENCHESSUCCESS+1
            BININDEXNEW=ENERGY2INDEX(CURRENTPOINTENERGY)
            !call ORDERQ4(Natoms, CurrentPointCoordinates, Q4orderparam)
            !CurrentQ=Q4orderparam
            !print *, BinIndexnew, CurrentPointEnergy, CurrentQ, Q4orderparam


            IF (BININDEXNEW < 1 .OR.BININDEXNEW>HBINS) THEN
               IF (DEBUG) PRINT *, 'Structure outside energy range. Rejecting move'
               ACCEPTMOVE=.FALSE.
            ELSE 
         	IF ((BINSTRUCTURES).AND.(MOD(NQUENCHES, SAVENTH).EQ.0)) THEN
             	   CALL SAVEBINSTRUCTURES(CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES, BININDEXNEW, MINIMANUMBER, .TRUE.)
           	 ENDIF
               CALL SAVEBINSTRUCTURES(CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES, BININDEXNEW, MINIMANUMBER, .FALSE.)
               LNRATIO=LNWEIGHT(BININDEXOLD)-LNWEIGHT(BININDEXNEW)
               CALL RANDOM_NUMBER(HARVEST)
               LNHARVEST=LOG(HARVEST)
               IF ((LNRATIO>0.D0).OR.(LNRATIO>LNHARVEST)) THEN
                  ACCEPTMOVE=.TRUE.
               ELSE
                  ACCEPTMOVE=.FALSE.
               ENDIF
            ENDIF

            IF (DEBUG) PRINT *, 'Eold, Enew, Iold, Inew, Wold, Wnew, Converged , Ratio'
            IF (DEBUG) PRINT *, OLDENERGY, CURRENTPOINTENERGY, BININDEXOLD, BININDEXNEW, LNWEIGHT(BININDEXOLD) , &
                     & LNWEIGHT(BININDEXNEW), CONVERGED, LNRATIO

            IF (ACCEPTMOVE) THEN 
               ! move accepted
               !print *, 'Moving from bin', BinIndexold, 'to bin', BinIndexnew

               VISITS(BININDEXNEW)=VISITS(BININDEXNEW)+1
               VISITSTOTAL(BININDEXNEW)=VISITSTOTAL(BININDEXNEW)+1
               LNWEIGHT(BININDEXNEW)=LNWEIGHT(BININDEXNEW)+LNMODFAC
               DISTANCEOLD=CALCULATEDDISTANCE(CURRENTPOINTCOORDINATES, PERTURBEDCOORDINATES)
!              CALL MINDGMIN(CURRENTPOINTCOORDINATES,PERTURBEDCOORDINATES,NATOMS,MINDISTANCEOLD,PERIODIC,TWOD)
               CALL NEWMINDIST(CURRENTPOINTCOORDINATES,PERTURBEDCOORDINATES,NATOMS,MINDISTANCEOLD,PERIODIC,TWOD, &
    &                            'AX    ',.FALSE.,RIGID,DEBUG,RMAT)

               DISTANCE(BININDEXNEW)=DISTANCE(BININDEXNEW)+DISTANCEOLD
               MINDISTANCE(BININDEXNEW)=MINDISTANCE(BININDEXNEW)+MINDISTANCEOLD
               BININDEXOLD=BININDEXNEW
               PERTURBEDCOORDINATESSAVE=PERTURBEDCOORDINATES
            ELSE
               ! move rejected
               !print *, 'Staying in bin', BinIndexold
               VISITS(BININDEXOLD)=VISITS(BININDEXOLD)+1
               VISITSTOTAL(BININDEXOLD)=VISITSTOTAL(BININDEXOLD)+1
               LNWEIGHT(BININDEXOLD)=LNWEIGHT(BININDEXOLD)+LNMODFAC
               DISTANCE(BININDEXOLD)=DISTANCE(BININDEXOLD)+DISTANCEOLD
               MINDISTANCE(BININDEXOLD)=MINDISTANCE(BININDEXOLD)+MINDISTANCEOLD
            ENDIF
         ELSE
         IF (DEBUG) PRINT *, 'Optimisation unsuccessful'
         ENDIF

         IF ((MOD(NQUENCHESSUCCESS,DUMPEVERYNTHQUENCH).EQ.0).AND.(NQUENCHESSINCELASTUPDATE > EQUIL)) THEN  
            ! record statistics
             OPEN(UNIT=421, FILE='lnWeight.his', status='unknown')
             DO I=1, HBINS
                WRITE(421, '(2G20.10)')  BinLabel(i), lnWeight(i)
             ENDDO
             CLOSE(421)
           
             OPEN(UNIT=421, FILE='Distance.his', status='unknown')
             DO I=1, HBINS
                WRITE(421, '(2G20.10)')  BinLabel(i), Distance(i)
             ENDDO
             CLOSE(421)
           
             OPEN(UNIT=421, FILE='MinDistance.his', status='unknown')
             DO I=1, HBINS
                WRITE(421, '(2G20.10)')  BinLabel(i), MinDistance(i)
             ENDDO
             CLOSE(421)
      
             OPEN(UNIT=422, FILE='Visits.his', status='unknown')
             DO I=1, HBINS
                WRITE(422, '(2G20.10)')  BinLabel(i), Visits(i)
             ENDDO
             CLOSE(422)
      
             OPEN(UNIT=422, FILE='VisitsTotal.his', status='unknown')
             DO I=1, HBINS
                WRITE(422, '(2G20.10)')  BinLabel(i), VisitsTotal(i)
             ENDDO
             CLOSE(422)
           
             OPEN(UNIT=422, FILE='lnModfac.restart', status='unknown')
             WRITE(422, '(G20.10)') lnModFac
             CLOSE(422)
           
             OPEN(UNIT=422, FILE='nWL.restart', status='unknown')
             WRITE(422, '(G20.10)') nWL
             CLOSE(422)
      
             OPEN(UNIT=421, FILE='MinimaNumber.his', status='unknown')
             DO I=1, HBINS
                WRITE(421, '(2G20.10)')  BinLabel(i), MinimaNumber(i)
             ENDDO
             CLOSE(421)


             CALL RECORD_STAT(LNWEIGHT,DISTANCE, MINDISTANCE, VISITSTOTAL, BINLABEL)

             ! checking histogram for convergence

             FLAT=CHECKFLATNESS(VISITS,LNMODFAC,NWL, FROZENVISITS)
             IF (FLAT) THEN
                NQUENCHESSINCELASTUPDATE=0
                PRINT *, '--===Visits histogram satisfied flatness criterion===--'
                LNMODFAC=HISTFACMUL*LNMODFAC
                PRINT *, 'lnModfac, exp', lnModFac, exp(lnModFac)
                PRINT *, 'Updating modification factor to', lnModFac
                FROZENVISITS=VISITS
                VISITS=0
                NWL=NWL+1
                
                OPEN(UNIT=422, FILE='lnModfac.restart', status='unknown')
                WRITE(422, '(G20.10)') lnModFac
                CLOSE(422)
              
                OPEN(UNIT=422, FILE='nWL.restart', status='unknown')
                WRITE(422, '(G20.10)') nWL
                CLOSE(422)
             ENDIF
      

              PRINT *, NQUENCHES, ' quenches completed'
         ENDIF

      NQUENCHES=NQUENCHES+1
      NQUENCHESSINCELASTUPDATE=NQUENCHESSINCELASTUPDATE+1
      ENDDO

      CONTAINS

!---======================================---`
      FUNCTION ENERGY2INDEX(CURRENTPOINTENERGY)
      USE COMMONS
      IMPLICIT NONE

      INTEGER ENERGY2INDEX

      REAL(8) CURRENTPOINTENERGY


      IF (NINT((CURRENTPOINTENERGY-HISTMIN)/HISTINT)<0) THEN
         ENERGY2INDEX=-1
      ELSE
         ENERGY2INDEX=NINT((CURRENTPOINTENERGY-HISTMIN)/HISTINT)+1
      ENDIF

      END FUNCTION ENERGY2INDEX


!---======================================---`
      REAL(8) FUNCTION GETDISPLACEMENT(STEP)
           IMPLICIT NONE
           REAL(8),INTENT(IN) :: STEP
           REAL(8) :: HARVEST
           CALL RANDOM_NUMBER(HARVEST)
           GETDISPLACEMENT=2.0D0*HARVEST*STEP-STEP
      END FUNCTION GETDISPLACEMENT

!---======================================---`
      FUNCTION PERTURBGEOMETRY(CURRENTPOINTCOORDINATES)
           USE COMMONS, ONLY: RADIUS,NATOMS,STEP
           IMPLICIT NONE

           REAL(8),DIMENSION(3*NATOMS,1) :: PERTURBGEOMETRY

           REAL(8),INTENT(IN) :: CURRENTPOINTCOORDINATES(3*NATOMS,1)
           INTEGER I, J

           PERTURBGEOMETRY=CURRENTPOINTCOORDINATES

           DO J=1, NATOMS
                DO
                     PERTURBGEOMETRY(3*(J-1)+1,1) = PERTURBGEOMETRY(3*(J-1)+1,1) + GETDISPLACEMENT(STEP(1))
                     PERTURBGEOMETRY(3*(J-1)+2,1) = PERTURBGEOMETRY(3*(J-1)+2,1) + GETDISPLACEMENT(STEP(1))
                     PERTURBGEOMETRY(3*(J-1)+3,1) = PERTURBGEOMETRY(3*(J-1)+3,1) + GETDISPLACEMENT(STEP(1))
                     IF ( PERTURBGEOMETRY(3*(J-1)+1, 1)**2 &
                      & + PERTURBGEOMETRY(3*(J-1)+2, 1)**2 &
                      & + PERTURBGEOMETRY(3*(J-1)+3, 1)**2 < RADIUS ) THEN
                          EXIT
                     ELSE
                          PERTURBGEOMETRY(3*(J-1)+1,1) = CURRENTPOINTCOORDINATES(3*(J-1)+1,1)
                          PERTURBGEOMETRY(3*(J-1)+2,1) = CURRENTPOINTCOORDINATES(3*(J-1)+2,1)
                          PERTURBGEOMETRY(3*(J-1)+3,1) = CURRENTPOINTCOORDINATES(3*(J-1)+3,1)
                     ENDIF
		ENDDO
           ENDDO
      END FUNCTION PERTURBGEOMETRY

!---======================================---`

      FUNCTION CALCULATEDDISTANCE(CURRENTPOINTCOORDINATES, PERTURBEDCOORDINATES)
      USE COMMONS
      IMPLICIT NONE

      INTEGER I
      REAL(8) CALCULATEDDISTANCE, CURRENTPOINTCOORDINATES(3*NATOMS,1), PERTURBEDCOORDINATES(3*NATOMS,1), &
              & X1(NATOMS, 1), Y1(NATOMS, 1), Z1(NATOMS, 1), X2(NATOMS, 1), Y2(NATOMS, 1), Z2(NATOMS, 1)

      
      DO I=1,NATOMS
         X1(I,1)=CURRENTPOINTCOORDINATES(3*(I-1)+1,1)
         Y1(I,1)=CURRENTPOINTCOORDINATES(3*(I-1)+2,1)
         Z1(I,1)=CURRENTPOINTCOORDINATES(3*(I-1)+3,1)
         X2(I,1)=PERTURBEDCOORDINATES(3*(I-1)+1,1)
         Y2(I,1)=PERTURBEDCOORDINATES(3*(I-1)+2,1)
         Z2(I,1)=PERTURBEDCOORDINATES(3*(I-1)+3,1)
      ENDDO

      CALCULATEDDISTANCE=DSQRT(SUM ( (X1(:,1)-X2(:,1))**2 + (Y1(:,1)-Y2(:,1))**2 + (Z1(:,1)-Z2(:,1))**2 ) )

      IF (DEBUG) PRINT *, 'Distance' , CalculatedDistance

      END FUNCTION CALCULATEDDISTANCE

!---======================================---`
      FUNCTION CHECKFLATNESS(VISITS,LNMODFAC,NWL, FROZENVISITS)
      USE COMMONS
      IMPLICIT NONE

      INTEGER VISITS(HBINS), N, I, VISITSNONZERO(HBINS),HMIN, HMAX, DELTAHK, NWL, FROZENVISITS(HBINS)
      LOGICAL CHECKFLATNESS, FLATBIN(HBINS)
      REAL(8) HDEV, LNMODFAC
      CHARACTER (LEN =256) FILENAME
      CHARACTER (LEN= 10)  ISTR


      CHECKFLATNESS=.FALSE.
      DELTAHK=0


      HMAX=MAXVAL(VISITS)
      DO I=1, HBINS
         IF (VISITS(I).NE.0) THEN
            VISITSNONZERO(I)=VISITS(I)
         ELSE
            VISITSNONZERO(I)=HUGE(HBINS)
         ENDIF
      ENDDO
      HMIN=MINVAL(VISITSNONZERO)
      IF (DEBUG) PRINT *, 'Min and max', Hmin, Hmax

      HDEV=(1.0D0*HMAX-HMIN)/(HMAX+HMIN)
    
      OPEN(UNIT=421, FILE='HistFluctuations', status='unknown', position='append')
      WRITE(421, '(G20.10)') HDev
      CLOSE(421)

      IF (VISITPROP) THEN
          !minimal number of visits should be proportional to 1/sqrt(ln(f))
 
          CHECKFLATNESS=.FALSE.
          FLATBIN=.FALSE.
 
          IF ((EQUIL.NE.-1).AND.(NWL.EQ.0)) THEN
             IF ((HDEV<HPERCENT).AND.(ABS(1.0D0*HMIN-HMAX)>HPERCENT)) CHECKFLATNESS=.TRUE.
          ELSE
             DO I=1, HBINS
                IF ((VISITS(I) > 1.0D0/SQRT(LNMODFAC)).OR.((VISITS(I).EQ.0).AND.(FROZENVISITS(I).EQ.0))) FLATBIN(I)=.TRUE.
             ENDDO
             IF (DEBUG)  PRINT *, 'Flatness by bin', FlatBin 
             IF (ALL(FLATBIN)) CHECKFLATNESS=.TRUE.
          ENDIF
      ELSE
         ! histogram flatness criterion: should be used for disconnected PES in the initial run
         ! where 1/lnf criterion will never be satisfied.
         IF ((HDEV<HPERCENT).AND.(ABS(1.0D0*HMIN-HMAX)>HPERCENT)) CHECKFLATNESS=.TRUE.
      ENDIF

      DO I=1, HBINS
         IF (VISITS(I).NE.0) THEN
            DELTAHK=DELTAHK+(VISITS(I)-HMIN)
         ELSE
            DELTAHK=DELTAHK+0
         ENDIF
      ENDDO
      OPEN(UNIT=421, FILE='HistSaturation', status='unknown', position='append')
      WRITE(421, '(G20.10)') deltaHk
      CLOSE(421)

      WRITE (ISTR, '(i10)') nWL
      FILENAME="HistSaturation."//trim(adjustl(istr))
      OPEN(UNIT=421,FILE=FILENAME, STATUS="unknown", form="formatted", position="append")
      WRITE(421, '(G20.10)') deltaHk
      CLOSE(421)

      IF (CHECKFLATNESS) THEN
         OPEN(UNIT=421, FILE='deltaHk.vs.lnModfac', status='unknown', position='append')
         WRITE(421, '(2G20.10)') lnModFac,deltaHk
         CLOSE(421)
      ENDIF

      END FUNCTION CHECKFLATNESS

!---======================================---`
      SUBROUTINE RECORD_STAT(LNG_J,DISTANCE, MINDISTANCE, VISITSTOTAL, BINLABEL)
      USE COMMONS
      IMPLICIT NONE

      INTEGER VISITSTOTAL(HBINS), I, T, FULLBINS

      REAL(8)  MINDISTANCE(HBINS), KB, BINLABEL(HBINS),  &
              & A, B, AP, BP, AVLNG_J, TOTLNG_J, MINIMAL_D_J, &
              & DELTAT, SUMD_JNONMIN, SUMD_JMIN, &
              & AVD_JNONMIN, AVD_JMIN, SUMD2_JNONMIN, SUMD2_JMIN, AVD2_JNONMIN, AVD2_JMIN, &
	      & STD_DEVD_JNONMIN, STD_DEVD_JMIN, STD_DEVSMOOTHD_J, AVSMOOTHD2_J, AVSMOOTHD_J, &
              & TOTLNP_J, AVLNP_J

      REAL(8) LNP_JMIN(HBINS), P_JNONMIN(HBINS), LNP_JNONMIN(HBINS), LNG_J(HBINS), KAPPA, &
              & D_JNONMIN(HBINS), DISTANCE(HBINS), D_JMIN(HBINS), SCALEDLNG_JNONMIN(HBINS),&
              & P_JNORM(HBINS), SMOOTHD_JNONMIN(HBINS) , G_J(HBINS), G_JNORM(HBINS), P_J(HBINS), SCALEDLNP_JNONMIN(HBINS)


      IF (PERIODIC) THEN
         KAPPA=3*NATOMS-3
      ELSE
         KAPPA=3*NATOMS-6
      ENDIF

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            D_JNONMIN(I)=DISTANCE(I)/VISITSTOTAL(I)
         ELSE
            D_JNONMIN(I)=0.0D0
         ENDIF
      ENDDO 

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            D_JMIN(I)=MINDISTANCE(I)/VISITSTOTAL(I)
         ELSE
            D_JMIN(I)=0.0D0
         ENDIF
      ENDDO 

      MINIMAL_D_J=1.0D100
      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
           IF (D_JNONMIN(I) < MINIMAL_D_J) MINIMAL_D_J=D_JNONMIN(I)
         ENDIF
      ENDDO

     ! calculate density of minima in the BS approximation: P_j=[G_j*(dm/d_j)**kappa]/Norm*deltaV 
     ! now we will actually be storing ln(p_j)

      FULLBINS=0
      TOTLNG_J=0.D0
      AVLNG_J=0.D0

      DO I=1, HBINS
         TOTLNG_J=TOTLNG_J+LNG_J(I)
         FULLBINS=FULLBINS+1
      ENDDO

      AVLNG_J=TOTLNG_J/FULLBINS

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            SCALEDLNG_JNONMIN(I)=LNG_J(I)-AVLNG_J
         ELSE
            SCALEDLNG_JNONMIN(I)=0.D0
         ENDIF
      ENDDO
 
      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            G_J(I)=EXP(SCALEDLNG_JNONMIN(I))
         ELSE
            G_J(I)=0.D0
         ENDIF
      ENDDO

	DO I=1,HBINS
	    IF (VISITSTOTAL(I).NE.0) THEN
	       G_JNORM(I)=G_J(I)/(SUM(G_J)*HISTINT)
	    ELSE
	       G_JNORM(I)=0.D0
            ENDIF
        ENDDO

         
      FULLBINS=0
      TOTLNP_J=0.D0
      AVLNP_J=0.D0

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            LNP_JNONMIN(I)=LNG_J(I)-(KAPPA+3)*LOG(D_JNONMIN(I))
            TOTLNP_J=TOTLNP_J+LNP_JNONMIN(I)
            FULLBINS=FULLBINS+1 
         ELSE 
            LNP_JNONMIN(I)=-HUGE(LNG_J(I))
         ENDIF
      ENDDO

      AVLNP_J=TOTLNP_J/FULLBINS

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            SCALEDLNP_JNONMIN(I)=LNP_JNONMIN(I)-AVLNP_J
         ELSE
            SCALEDLNP_JNONMIN(I)=0.D0
         ENDIF
      ENDDO

      DO I=1, HBINS
         IF (VISITSTOTAL(I).NE.0) THEN
            P_JNONMIN(I)=EXP(SCALEDLNP_JNONMIN(I))
         ELSE
            P_JNONMIN(I)=0.D0
         ENDIF
      ENDDO

        DO I=1,HBINS
            IF (VISITSTOTAL(I).NE.0) THEN
               P_JNORM(I)=P_JNONMIN(I)/(SUM(P_JNONMIN)*HISTINT)
            ELSE
               P_JNORM(I)=0.D0
            ENDIF
        ENDDO


      OPEN(UNIT=421, FILE='BL.Pjnorm.lnGj.Djnm.Djm.VT.his', status='unknown')
      DO I=1, HBINS
         WRITE(421, '(6G20.10)')  BinLabel(i), p_jnorm(i), lng_j(i), d_jnonmin(i), d_jmin(i), VisitsTotal(i)
      ENDDO
         CLOSE(421)

      END SUBROUTINE RECORD_STAT

!---======================================---
    SUBROUTINE SAVEBINSTRUCTURES(CURRENTPOINTENERGY, CURRENTPOINTCOORDINATES, BININDEX, MINIMANUMBER, WRITESTRUCT)
        USE COMMONS, ONLY:NATOMS, HBINS,BQMAX, CHRMMT, ZSYM
        USE MODCHARMM, ONLY: OSASAT, RPRO, ODIHET
        IMPLICIT NONE

        INTEGER BININDEX, JB,  MINIMANUMBER(HBINS), JC
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
             IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.BQMAX) THEN
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
	             IF (ABS(BINENERGIES(JB)-CURRENTPOINTENERGY).LE.BQMAX) THEN
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



      END SUBROUTINE BASINSAMPLING

