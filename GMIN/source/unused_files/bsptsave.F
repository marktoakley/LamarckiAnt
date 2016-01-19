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
!---======================================---
      SUBROUTINE PTBASINSAMPLING

      USE MODCHARMM
      USE COMMONS
      USE TETHERFUNC

      IMPLICIT NONE
#ifdef MPI

      INCLUDE 'mpif.h'
       

      
      INTEGER, PARAMETER :: NHIST=100, NHISTE=1000
      INTEGER :: IACCEPT(0:NTRAJ-1), MPIERR, I,J,K, IS(MPI_STATUS_SIZE)
      INTEGER NHISTQ4(NHIST,0:NTRAJ-1), NDUMMY, NTYPEA, NTOT, NH, IMESG, IQE, IQ4, IQ6, J1, J2,
     1        NHISTQ6(NHIST,0:NTRAJ-1), NHISTQE(NHISTE, 0:NTRAJ-1), MSTEP, NSTAT, IENR, 
     2        NACCEPTPT(0:NTRAJ-1), NOUT(0:NTRAJ-1), ITRAJ, ITRAJO,NEACCEPT, RNDSEED, NUPDATE,
     3        CONVERGED,LBFGS_ITERATIONS, JD, BININDEX, MINIMANUMBER(HBINS,0:NTRAJ-1),
     4        NHISTALLQ(NHIST, NHIST, 0:NTRAJ-1),
     5        NHISTGRAND(NHISTE, NHIST, NHIST, 0:NTRAJ-1),
     6        Visits(hbins, 0:NTRAJ-1)

      REAL(8) V(NATOMS), VO(NATOMS), TEMPTRAJ(0:NTRAJ-1), H(0:NTRAJ-1), BETA(0:NTRAJ-1), 
     1        EAV(0:NTRAJ-1), EAV2(0:NTRAJ-1), Q(3,NATOMS), Q4AV(0:NTRAJ-1), Q4AV2(0:NTRAJ-1), 
     2        Q6AV(0:NTRAJ-1), Q6AV2(0:NTRAJ-1),VMIN(0:NTRAJ-1), VMAX(0:NTRAJ-1), VENR(NENRPER), 
     3        HINIT(0:NTRAJ-1), EPSAB, EPSBB, SIGAB, SIGBB, X(NATOMS), Y(NATOMS), Z(NATOMS), 
     4        CTE, T, VINIT,VFINAL, POTEL, GRAD(3*NATOMS), Q4, Q6, RANDOM, DPRAND, Q4MAX, Q6MAX , 
     5	      DQ4, DQ6, DHISTE, ENUL, XO(NATOMS), YO(NATOMS), ZO(NATOMS), DDX, DDY, DDZ, DE, 
     6        W, WCOMP, WAC, E, ER, DBETA, DELTA, CV, FQ4, FQ6, RMAX, DDXN, DDYN, DDZN, R2, 
     7        DUMMY,BINLABEL(HBINS), VFINALSAVE, DIHEORDERPARAM,SASAORDERPARAM, 
     8        diheorderparam_av(0:NTRAJ-1), SASAorderparam_av(0:NTRAJ-1), MINCOORDS(3*NATOMS,NTRAJ), interimCOORDS(3*NATOMS,1),
     9        interimCOORDS2(3*NATOMS,1), DOSSTATS(MAXIT,2)

              

       REAL(8)  diheorderparam_min, diheorderparam_max, SASAorderparam_min, SASAorderparam_max, 
     1          ddihe, dsasa, rgyr_min, rgyr_max, eint_min, eint_max, drgyr, deint,
     2          rgyr_av(0:NTRAJ-1), eint_av(0:NTRAJ-1), rgyr_av2(0:NTRAJ-1), eint_av2(0:NTRAJ-1),
     3          orderparam1, orderparam2, Distance(HBINS, 0:NTRAJ-1),
     4          Distance_av(HBINS,0:NTRAJ-1), DISTANCEOLD, norm_pj(0:NTRAJ-1), Pj(HBINS,0:NTRAJ-1),
     5          VMINPREV, VMINNEW, DELTA_EAV(0:NTRAJ-1),
     6        diheorderparam_av2(0:NTRAJ-1), SASAorderparam_av2(0:NTRAJ-1), RGYR, EINT



      CHARACTER (LEN =256)  FILENAME, FILENAME2,FILENAME3,FILENAME4,FILENAME5,FILENAME6,
     1                      FILENAME7,FILENAME8,FILENAME9, FILENAME10, FILENAME11, FILENAME12, 
     2                      FILENAME100,FILENAME101,FILENAME102, FILENAME103, FILENAME104, FILENAME105 
      CHARACTER (LEN= 10)  ISTR, CDUMMY
      LOGICAL EXCHANGE, EXCHANGEACCEPT, FITS, NEWENERGY, EVAP, EVAPREJECT, ACCEPTED


      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /MYPOT/ POTEL
      COMMON /EV/ EVAP, EVAPREJECT
     

!      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NDUMMY,MPIERR)
      IF ((NDUMMY.NE.NPAR).OR.(NDUMMY.NE.NTRAJ)) THEN
         WRITE(MYUNIT,'(A)')  'Number of temperature trajectories does not correspond to the number  
     1           of processors. Stop.'
         RETURN
      ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)

      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Processor", mynode, " of", NPAR, " speaking:"
      WRITE(MYUNIT, '(A,I10)') 'Number of atoms', natoms
      IF (PERIODIC) THEN
         WRITE(MYUNIT, '(A,6G20.10)') 'Binary data', ntypea, epsab, epsbb, sigab, sigbb, cutoff
         WRITE(MYUNIT, '(A,3G20.10)') 'Box data', boxlx, boxly, boxlz
      ELSEIF(CHRMMT) THEN
         WRITE(MYUNIT, '(A)') 'CHARMM job'
         WRITE(MYUNIT, '(A,G20.10)') 'Starting stepsize ', STEP(MYNODE+1)
      ELSE
         WRITE(MYUNIT, '(A,G20.10)') 'Radius**2', radius
      ENDIF

      IF (CHRMMT) THEN
         IF (HISTRESTART) THEN
           WRITE (ISTR, '(i10)') MYNODE+1
           FILENAME12="xyz_out."//trim(adjustl(istr))
           OPEN(UNIT=1986+MYNODE,FILE=FILENAME12, STATUS="old", form="formatted")
               READ(1986+MYNODE, '(A,G20.10)') cdummy, dummy
            DO J1=1,NATOMS
               J2=3*(J1-1)
               READ(1986+MYNODE,*) COORDS(J2+1,MYNODE+1), COORDS(J2+2,MYNODE+1), COORDS(J2+3,MYNODE+1)
            ENDDO
            CLOSE(1986+MYNODE)
         ELSE
             DO J1=1,NATOMS
                COORDS(3*(J1-1)+1,mynode+1)=COORDS(3*(J1-1)+1,1)
                COORDS(3*(J1-1)+2,mynode+1)=COORDS(3*(J1-1)+2,1)
                COORDS(3*(J1-1)+3,mynode+1)=COORDS(3*(J1-1)+3,1)
             ENDDO
         ENDIF
      ENDIF

      DO I=1, HBINS
         BINLABEL(I)=HISTMIN + HISTINT*(I-0.5)
      ENDDO


      ITRAJ=MYNODE
      NEACCEPT=0
      MINIMANUMBER=0

      ! Initialisation

      DO I=1,NATOMS
         X(I)=COORDS(3*(I-1)+1,mynode+1)
         Y(I)=COORDS(3*(I-1)+2,mynode+1)
         Z(I)=COORDS(3*(I-1)+3,mynode+1)
      ENDDO
      DO I=1,NATOMS
         Q(1,I)=X(I)
         Q(2,I)=Y(I)
         Q(3,I)=Z(I)
      ENDDO

      WRITE(MYUNIT, '(A,3G20.10)') 'xyz', x(1), y(1), z(1) 
 
      CTE=(LOG(PTTMAX/PTTMIN))/(NTRAJ-1)
      CTE=EXP(CTE)

      DO I=0, NTRAJ-1
         TEMPTRAJ(I)=PTTMIN*CTE**I
         T=TEMPTRAJ(I)
         BETA(I)=1.0D0/T

         IF (T.LT.1.0d0) THEN
            H(I)=(16./(1.-SQRT(T)))**(1./6.)-(16./(1.+SQRT(T)))**(1./6.)
            H(I)=0.46*H(I)
         ELSE
            H(I)=1.
         ENDIF
         H(I)=H(I)/2.
      !   IF (BINSTRUCTURES) H(I)=H(I)*3.0
      ENDDO

      WRITE(MYUNIT, '(A,2G20.10)') "Temp, H", TEMPTRAJ(MYNODE), H(MYNODE) 

      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
      VINIT=POTEL

      IF (PERIODIC) THEN
         CALL QORDER_BLJ(Q,Q4,Q6)
      ELSE IF (CHRMMT) THEN
           IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
           IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
           IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.false.) 
           IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
      ELSE
         CALL QORDER_LJ(Q,Q4,Q6)
      ENDIF

      WRITE(MYUNIT, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NTRAJ-1)
      WRITE(MYUNIT, '(A, G20.10)') 'This temperature trajectory=', TEMPTRAJ(MYNODE)
      WRITE(MYUNIT, '(A, G20.10)') 'Starting E=', VINIT
      VMINPREV=VINIT
      IF (PERIODIC) WRITE(MYUNIT, '(A, 2G20.10)') 'Starting Q4, Q6=', Q4, Q6
      IF (CHRMMT) THEN
        IF (OSASAT)  WRITE(MYUNIT, '(A, 2G20.10)') 'Starting Met-enk order params=', 
     1                                         diheorderparam, SASAorderparam 
        IF (ORGYT)  WRITE(MYUNIT, '(A, G20.10)') 'Starting Rgyr order params=', 
     1                                          RGYR 
        IF (OEINTT) WRITE(MYUNIT, '(A, G20.10)') 'Starting Eint order params=',
     1                                          EINT
      ENDIF

      ! Initialisation complete

      RNDSEED=2002+MYNODE
      CALL SDPRND(RNDSEED)
      RANDOM=DPRAND()
      WRITE(MYUNIT, '(A, G20.10)') 'Starting random number=', RANDOM

      NTOT=0

      DO I=0, NTRAJ-1
         HINIT(I)=H(I)
      ENDDO 

      NACCEPTPT(MYNODE)=0
      IACCEPT(MYNODE)=0
      EAV(MYNODE)=0.
      EAV2(MYNODE)=0.
      Q4AV(MYNODE)=0.
      Q6AV(MYNODE)=0.
      Q4AV2(MYNODE)=0.
      Q6AV2(MYNODE)=0.
      diheorderparam_av(mynode)=0.
      SASAorderparam_av(mynode)=0.
      diheorderparam_av2(mynode)=0.
      SASAorderparam_av2(mynode)=0.
      rgyr_av(mynode)=0.
      eint_av(mynode)=0.
      rgyr_av2(mynode)=0.
      eint_av2(mynode)=0.
      DO I=1,NHIST
         NHISTQ4(I,MYNODE)=0
         NHISTQ6(I,MYNODE)=0
      ENDDO
      DO I=1,NHISTE
         NHISTQE(I,MYNODE)=0
      ENDDO
      NHISTALLQ(:,:,MYNODE)=0
      NHISTGRAND(:,:,:,MYNODE)=0
      Visits(:,MYNODE)=0
      Distance(:,MYNODE)=0
      Distance_av(:,MYNODE)=0
      DELTA_EAV(MYNODE)=0.0d0

      Q4MAX=0.1
      Q6MAX=0.5 ! will have to change as is rather system-specific
      diheorderparam_min=0.4d0
      diheorderparam_max=1.0d0
      SASAorderparam_min=300.0d0
      SASAorderparam_max=800.0d0
      rgyr_min=4.0d0
      rgyr_max=12.0d0
      eint_min=-35.0
      eint_max=10.0
      ddihe=(diheorderparam_max-diheorderparam_min)/(NHIST-1)
      dsasa=(SASAorderparam_max-SASAorderparam_min)/(NHIST-1)
      drgyr=(rgyr_max-rgyr_min)/(NHIST-1)
      deint=(eint_max-eint_min)/(NHIST-1)
      DQ4=Q4MAX/(NHIST-1)
      DQ6=Q6MAX/(NHIST-1)
      DHISTE=(PTEMAX-PTEMIN)/(NHISTE-1)

      VMIN(MYNODE)=-VINIT
      VMAX(MYNODE)=VINIT
      ENUL=VINIT

      IF (MYNODE.EQ.0) THEN
         OPEN(UNIT=10, FILE='distributions.header', FORM='formatted')
         WRITE(10,*) NATOMS
         WRITE(10,*) NTRAJ, PTSTEPS/NENRPER, NENRPER
         WRITE(10,*) ENUL
         CLOSE(10)
      ENDIF 
      
      WRITE (ISTR, '(i10)') MYNODE+1
      IF (.NOT.BINSTRUCTURES) THEN
         FILENAME2="config."//trim(adjustl(istr))
         OPEN(UNIT=11+MYNODE,FILE=FILENAME2, form="formatted")
      ENDIF

      NUPDATE=100
      MSTEP=NATOMS ! is it a general rule?
      NSTAT=MSTEP*(PTSTEPS+NEQUIL)

      DO I=1, PTSTEPS+NEQUIL
         DO J=1,MSTEP
         ! Propagate

            IF (CHRMMT) THEN
                DO K=1,NATOMS
                  COORDS(3*(K-1)+1,mynode+1)=X(K)
                  COORDS(3*(K-1)+2,mynode+1)=Y(K)
                  COORDS(3*(K-1)+3,mynode+1)=Z(K)
                ENDDO
                DO K=1, NATOMS
                   XO(K)=X(K)
                   YO(K)=Y(K)
                   ZO(K)=Z(K)
                ENDDO
                  IF (FIXEDENDMOVET) THEN 
               ! step in internal coordinates
                   IF(CHRIGIDROTT) CALL MKRIGIDROT(mynode+1)
                   IF(CHRIGIDTRANST) CALL MKRIGIDTRANS(mynode+1)
                   CALL TAKESTEPCH(mynode+1)
                   ! new geometry is now in coords
                   IF (FIXCOM) CALL CENTRECOM(COORDS(1:3*NATOMS,mynode+1))
                   DO K=1,NATOMS
                      X(K)=COORDS(3*(K-1)+1,mynode+1)
                      Y(K)=COORDS(3*(K-1)+2,mynode+1)
                      Z(K)=COORDS(3*(K-1)+3,mynode+1)
                   ENDDO
                   FITS=.true.
                  ENDIF
            ELSE IF (BINSTRUCTURES) THEN
                DO K=1, NATOMS
                   XO(K)=X(K)
                   YO(K)=Y(K)
                   ZO(K)=Z(K)
                ENDDO
                   DO K=1,NATOMS
                      interimCOORDS(3*(K-1)+1,1)=X(K)
                      interimCOORDS(3*(K-1)+2,1)=Y(K)
                      interimCOORDS(3*(K-1)+3,1)=Z(K)
                   ENDDO
               interimCoords2=PerturbGeometry(interimCoords)
!               WRITE(MYUNIT, '(I10, A)') j, ' step pert'
!                  CALL FLUSH(MYUNIT)
                   DO K=1,NATOMS
                      X(K)=interimCOORDS2(3*(K-1)+1,1)
                      Y(K)=interimCOORDS2(3*(K-1)+2,1)
                      Z(K)=interimCOORDS2(3*(K-1)+3,1)
                   ENDDO

                 DO K=1,NATOMS
                   COORDS(3*(K-1)+1,mynode+1)=X(K)
                   COORDS(3*(K-1)+2,mynode+1)=Y(K)
                   COORDS(3*(K-1)+3,mynode+1)=Z(K)
                 ENDDO
            ELSE
                DO K=1, NATOMS
                   XO(K)=X(K)
                   YO(K)=Y(K)
                   ZO(K)=Z(K)
                ENDDO
7777            RANDOM=DPRAND()
                H=INT(1+NATOMS*RANDOM) ! checked that it covers the whole range
                IF (NH.LE.0.OR.NH.GT.NATOMS) GOTO 7777
             
                ! move randomely selected particle - will have to be generalized
                RANDOM=DPRAND()
                DDX=H(MYNODE)*(RANDOM-0.5)
                RANDOM=DPRAND()
                DDY=H(MYNODE)*(RANDOM-0.5)
                RANDOM=DPRAND()
                DDZ=H(MYNODE)*(RANDOM-0.5)
             
                IF (PERIODIC) THEN
                   X(NH)=X(NH)+DDX
                   Y(NH)=Y(NH)+DDY
                   Z(NH)=Z(NH)+DDZ
                   X(NH)=X(NH)-BOXLX*ANINT(X(NH)/BOXLX) 
                   Y(NH)=Y(NH)-BOXLY*ANINT(Y(NH)/BOXLY) 
                   Z(NH)=Z(NH)-BOXLZ*ANINT(Z(NH)/BOXLZ) 
                   DO K=1,NATOMS
                      COORDS(3*(K-1)+1,mynode+1)=X(K)
                      COORDS(3*(K-1)+2,mynode+1)=Y(K)
                      COORDS(3*(K-1)+3,mynode+1)=Z(K)
                   ENDDO
                   FITS=.TRUE.
                 ELSE ! e.g. LJ 
                   DDXN=DDX/NATOMS
                   DDYN=DDY/NATOMS
                   DDZN=DDZ/NATOMS
             
                   RMAX=0.
                   DO K=1,NATOMS
                      IF (K.NE.NH) THEN
                         X(K)=X(K)-DDXN
                         Y(K)=Y(K)-DDYN
                         Z(K)=Z(K)-DDZN
                      ELSE
                         X(K)=X(K)-DDXN+DDX
                         Y(K)=Y(K)-DDYN+DDY
                         Z(K)=Z(K)-DDZN+DDZ
                      ENDIF
             
                      R2=X(K)**2+Y(K)**2+Z(K)**2
                      RMAX=MAX(R2,RMAX)
                    ENDDO
                    DO K=1,NATOMS
                      COORDS(3*(K-1)+1,mynode+1)=X(K)
                      COORDS(3*(K-1)+2,mynode+1)=Y(K)
                      COORDS(3*(K-1)+3,mynode+1)=Z(K)
                    ENDDO
                    IF (RMAX.GT.RADIUS) THEN
                       FITS=.FALSE.
                    ELSE
                       FITS=.TRUE.
                    ENDIF
                 ENDIF
            ENDIF


            VFINAL=0.D0
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
            VFINAL=POTEL

       IF (.NOT.BINSTRUCTURES) THEN
            DELTA=VFINAL-VINIT
            WCOMP=DELTA*BETA(MYNODE)
            W=MIN(1.0d0,EXP(-WCOMP))
            RANDOM=DPRAND()
            IF ((RANDOM.GT.W).OR.(.NOT.FITS)) THEN !reject move
               DO K=1, NATOMS
                  X(K)=XO(K)
                  Y(K)=YO(K)
                  Z(K)=ZO(K)
               ENDDO
               DO K=1,NATOMS
                 COORDS(3*(K-1)+1,mynode+1)=X(K)
                 COORDS(3*(K-1)+2,mynode+1)=Y(K)
                 COORDS(3*(K-1)+3,mynode+1)=Z(K)
               ENDDO
               VFINAL=VINIT
               IF (.NOT.FITS) NOUT(MYNODE)=NOUT(MYNODE)+1
            ELSE
               NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1
               IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
            ENDIF
         ENDIF
            
            IF (I.GT.NEQUIL) THEN
               VMIN(MYNODE)=MIN(VMIN(MYNODE),VINIT)
               VMAX(MYNODE)=MAX(VMAX(MYNODE), VINIT)
            ENDIF
!Quenching part
            IF ((BINSTRUCTURES).AND.(J.eq.1)) THEN 
               COORDSO(:,MYNODE+1)= COORDS(:,MYNODE+1)
               CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS, 
     1           DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1),DOSSTATS)
             IF (CONVERGED.NE.1) THEN
                WRITE(MYUNIT, '(A)') 'QUENCH DID NOT CONVERGE' 
                DO K=1, NATOMS
                    X(K)=XO(K)
                    Y(K)=YO(K)
                    Z(K)=ZO(K)
                ENDDO
             ENDIF
             IF ((CONVERGED.EQ.1).AND.(.NOT.EVAPREJECT)) THEN
             MINCOORDS(:,MYNODE+1)= COORDS(:,MYNODE+1)
                VMINNEW=POTEL
                DELTA=VMINNEW-VMINPREV
                WCOMP=DELTA*BETA(MYNODE)
                W=MIN(1.0d0,EXP(-WCOMP))
                RANDOM=DPRAND()
               WRITE(MYUNIT, '(I10, A, A, G20.10, A, G20.10)') i, ' step ',
     1          ' En= ', VMINNEW, ' Eo= ',VMINPREV
                  CALL FLUSH(MYUNIT)
                IF (RANDOM.GT.W) THEN !reject move
                   ACCEPTED=.false.
                   DO K=1, NATOMS
                      X(K)=XO(K)
                      Y(K)=YO(K)
                      Z(K)=ZO(K)
                   ENDDO
                   DO K=1,NATOMS
                     COORDS(3*(K-1)+1,mynode+1)=X(K)
                     COORDS(3*(K-1)+2,mynode+1)=Y(K)
                     COORDS(3*(K-1)+3,mynode+1)=Z(K)
                   ENDDO
                   VMINNEW=VMINPREV
                   IF (.NOT.FITS) NOUT(MYNODE)=NOUT(MYNODE)+1
                ELSE
                   ACCEPTED=.true.
                   NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1
                   IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
                ENDIF
               VMINPREV=VMINNEW 
                BININDEX=Energy2Index(VMINNEW, BINLABEL(1))
                IF ((MOD(I, SAVENTH).EQ.0).AND.ACCEPTED) THEN 
                       CALL SAVEBINSTRUCTURESMPI(POTEL,MINCOORDS(:,MYNODE+1),BININDEX,.true., MYNODE+1, 
     1                    NEWENERGY)
                       IF (NEWENERGY) THEN
                           MINIMANUMBER(BININDEX, MYNODE)=MINIMANUMBER(BININDEX, MYNODE)+1
                       ENDIF
                ENDIF
                ! IF (I.GT.NEQUIL) THEN
                 !  IF (BININDEX.GT.0.AND.BININDEX.LE.HBINS) then
                       VISITS(BININDEX,MYNODE)=VISITS(BININDEX,MYNODE)+1
                       DISTANCEOLD=CALCULATEDDISTANCE(COORDS(:,MYNODE+1), COORDSO(:,MYNODE+1))
                       DISTANCE(BININDEX,MYNODE)=DISTANCE(BININDEX,MYNODE)+DISTANCEOLD
                 !  endif
                ! ENDIF
                 DO K=1,NATOMS ! coords have to be resaved because we have not to use noreset to be able to calc distance 
                    COORDS(3*(K-1)+1,mynode+1)=X(K)
                    COORDS(3*(K-1)+2,mynode+1)=Y(K)
                    COORDS(3*(K-1)+3,mynode+1)=Z(K)
                 ENDDO
            ENDIF ! this part executed only if quench converged
            ENDIF ! binstructures endif

            VINIT=VFINAL
         ENDDO

       
         IF (MOD(I,NUPDATE).EQ.0) THEN ! update mc steps
            WAC=1.0*IACCEPT(MYNODE)/NUPDATE/MSTEP
            IF (WAC.LT.0.4) THEN
               H(MYNODE)=H(MYNODE)*0.9
            ENDIF
            IF (WAC.GT.0.6) THEN
               H(MYNODE)=H(MYNODE)*1.1
            ENDIF
!           IF (CHRMMT) THEN
!               IF ( (1.0*IACCEPT(MYNODE)/NUPDATE/MSTEP).GT.0.50d0) THEN 
!                  STEP(MYNODE+1)=STEP(MYNODE+1)*1.05
!                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=TRANSMAX*1.05
!                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=ROTMAX*1.05
!               ELSE
!                  STEP(MYNODE+1)=STEP(MYNODE+1)*(1.0D0/1.05)
!                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=TRANSMAX*(1.0D0/1.05)
!                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=ROTMAX*(1.0D0/1.05)
!               ENDIF
!           ENDIF
           IACCEPT(MYNODE)=0
         ENDIF
         
         IF (.NOT.BINSTRUCTURES) E=VFINAL
         IF (BINSTRUCTURES) E=VMINNEW
         IF(MYNODE.EQ.0) THEN
            RANDOM=DPRAND()
            J=(NTRAJ-1)*RANDOM
            RANDOM=DPRAND()
            IF (RANDOM.GT.EXCHPROB) THEN ! 0.1 probability of exchange 
               J=-2
               EXCHANGE=.FALSE.
            ELSE
               EXCHANGE=.TRUE.
            ENDIF
          ENDIF
          CALL MPI_BCAST(J,1,MPI_INTEGER,0,MPI_COMM_WORLD, MPIERR)
          IF (MYNODE.EQ.(J+1)) THEN
!             CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J+1,0,MPI_COMM_WORLD,IS,MPIERR)
             CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J,0,MPI_COMM_WORLD,MPIERR)
          ENDIF
          IF (MYNODE.EQ.J) THEN
             CALL MPI_RECV(ER,1,MPI_DOUBLE_PRECISION,J+1,0,MPI_COMM_WORLD,IS,MPIERR)

             DBETA=BETA(J)-BETA(J+1)
             DELTA=E-ER
             W=MIN(1.0D0,DEXP(DELTA*DBETA))
             NTOT=NTOT+1 
             RANDOM=DPRAND()
             IF (W.GT.RANDOM) THEN
                EXCHANGEACCEPT=.TRUE.
                IMESG=1
                CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_RECV(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(X,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(Y,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(Z,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
                E=ER
                NEACCEPT=NEACCEPT+1
              ELSE
                EXCHANGEACCEPT=.FALSE.
                IMESG=0
                CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
              ENDIF
          ENDIF
          IF (MYNODE.EQ.(J+1)) THEN
             CALL MPI_RECV(IMESG,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
             NTOT=NTOT+1
             IF (IMESG.EQ.1) THEN
                CALL MPI_RECV(ITRAJO,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(XO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(YO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(ZO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(E,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
                DO K=1, NATOMS
                   X(K)=XO(K)
                   Y(K)=YO(K)
                   Z(K)=ZO(K)
                ENDDO
                ITRAJ=ITRAJO
                NEACCEPT=NEACCEPT+1
              ENDIF
           ENDIF

           IF (.NOT.BINSTRUCTURES) VINIT=E
           IF (BINSTRUCTURES) VMINPREV=E

           IF (I.GT.NEQUIL) THEN 
              EAV(MYNODE)=EAV(MYNODE)+E
              EAV2(MYNODE)=EAV2(MYNODE)+E**2
              IQE=INT((E-PTEMIN)/DHISTE+1)
              IF (IQE.GT.0.AND.IQE.LT.NHISTE) THEN
                 NHISTQE(IQE,MYNODE)=NHISTQE(IQE,MYNODE)+1
              ENDIF
              DO K=1,NATOMS
                 Q(1,K)=X(K)
                 Q(2,K)=Y(K)
                 Q(3,K)=Z(K)
              ENDDO
              IF (PERIODIC) THEN
                 CALL QORDER_BLJ(Q,Q4,Q6)
              ELSE IF (CHRMMT) THEN
                   IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
                   IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                     COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
                   IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                     COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.false.) 
                   IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                     COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
              ELSE
                 CALL QORDER_LJ(Q,Q4,Q6)
              ENDIF
              Q4AV(MYNODE)=Q4AV(MYNODE)+Q4
              Q4AV2(MYNODE)=Q4AV2(MYNODE)+Q4**2
              Q6AV(MYNODE)=Q6AV(MYNODE)+Q6
              Q6AV2(MYNODE)=Q6AV2(MYNODE)+Q6**2


              DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)+DIHEORDERPARAM
              DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)+DIHEORDERPARAM**2
              SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)+SASAORDERPARAM
              SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)+SASAORDERPARAM**2
              RGYR_AV(MYNODE)=RGYR_AV(MYNODE)+RGYR
              RGYR_AV2(MYNODE)=RGYR_AV2(MYNODE)+RGYR**2
              EINT_AV(MYNODE)=EINT_AV(MYNODE)+EINT
              EINT_AV2(MYNODE)=EINT_AV2(MYNODE)+EINT**2
              IF (CHRMMT) THEN
                  IF (ODIHET) THEN
                     IQ4=INT((DIHEORDERPARAM-DIHEORDERPARAM_min)/ddihe+1)
                     IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                        NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                     ENDIF
                     IQ6=INT((SASAORDERPARAM-SASAORDERPARAM_min)/dsasa+1)             
                     IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                        NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                     ENDIF
                  ELSE IF (ORGYT) THEN
                     IQ4=INT((rgyr-rgyr_min)/drgyr+1)
                     IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                        NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                     ENDIF
                     IQ6=INT((eint-eint_min)/deint+1)             
                     IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                        NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                     ENDIF
                  ENDIF 
              ELSE
                  IQ4=INT(Q4/DQ4+1)
                  IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                     NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                  ENDIF
                  IQ6=INT(Q6/DQ6+1)             
                  IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                     NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                  ENDIF
              ENDIF

! Free energy statistic
              IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                 NHISTALLQ(IQ4,IQ6,MYNODE)=NHISTALLQ(IQ4,IQ6,MYNODE)+1
               ENDIF
               IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST.AND.
     1            IQE.GT.0.AND.IQE.LT.NHISTE) THEN
                  NHISTGRAND(IQE,IQ4,IQ6,MYNODE)=NHISTGRAND(IQE,IQ4,IQ6,MYNODE)+1
               ENDIF

              
              IENR=IENR+1
              IF (.NOT.BINSTRUCTURES) THEN
                 IF (IENR.EQ.NENRPER) THEN
                    IENR=0
                    CALL FLUSH(11+MYNODE)
                    WRITE(11+MYNODE,*) '#E=',E
                    DO K=1,NATOMS
                       WRITE(11+MYNODE,*) X(K),Y(K),Z(K)
                    ENDDO
                  ENDIF 
               ENDIF
           ENDIF                   
      ENDDO

!      CALL MPI_FINALIZE(MPIERR)
 
      ! computing the averages
      EAV(MYNODE)=EAV(MYNODE)/PTSTEPS
      EAV2(MYNODE)=EAV2(MYNODE)/PTSTEPS
      CV=(EAV2(MYNODE)-EAV(MYNODE)**2)*BETA(MYNODE)**2

      WRITE (ISTR, '(i10)') MYNODE+1
      FILENAME3="T.Ev.Cv.Ev2.Steps."//trim(adjustl(istr))
      OPEN(UNIT=41+MYNODE,FILE=FILENAME3, STATUS="unknown", form="formatted")
      IF (CHRMMT) THEN
         WRITE(41+MYNODE,'(6G20.10)') TEMPTRAJ(MYNODE),TEMPTRAJ(MYNODE)/0.001987, 
     1                                EAV(MYNODE), CV, EAV2(MYNODE),PTSTEPS
      ELSE
         WRITE(41+MYNODE,'(5G20.10)') TEMPTRAJ(MYNODE),EAV(MYNODE),CV, EAV2(MYNODE),PTSTEPS
      ENDIF
      CALL FLUSH(41+MYNODE)
      CLOSE(41+MYNODE)

      Q4AV(MYNODE)=Q4AV(MYNODE)/PTSTEPS
      Q6AV(MYNODE)=Q6AV(MYNODE)/PTSTEPS
      Q4AV2(MYNODE)=Q4AV2(MYNODE)/PTSTEPS
      Q6AV2(MYNODE)=Q6AV2(MYNODE)/PTSTEPS
      DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)/PTSTEPS
      SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)/PTSTEPS
      DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)/PTSTEPS
      SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)/PTSTEPS
      rgyr_AV(MYNODE)=rgyr_AV(MYNODE)/PTSTEPS
      eint_AV(MYNODE)=eint_AV(MYNODE)/PTSTEPS
      rgyr_AV2(MYNODE)=rgyr_AV2(MYNODE)/PTSTEPS
      eint_AV2(MYNODE)=eint_AV2(MYNODE)/PTSTEPS
      FQ4=SQRT(Q4AV2(MYNODE)-Q4AV(MYNODE)**2)
      FQ6=SQRT(Q6AV2(MYNODE)-Q6AV(MYNODE)**2)
      IF (BINARY) THEN  
         FILENAME4="T.Q4Av.Q6Av.Q4Av2.Q6Av2.Steps."//trim(adjustl(istr))
         OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
         WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), Q4AV(MYNODE), Q6AV(MYNODE), Q4AV2(MYNODE),  
     1                        Q6AV2(MYNODE), PTSTEPS
         CLOSE(1980)
      ENDIF
      IF (CHRMMT) THEN 
         IF (ODIHET) THEN 
            FILENAME4="T.DiheAv.SasaAv.DiheAv2.SasaAv2.Steps."//trim(adjustl(istr))
            OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
            WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), DIHEORDERPARAM_AV(MYNODE), 
     1              SASAORDERPARAM_AV(MYNODE), DIHEORDERPARAM_AV2(MYNODE),  
     2              SASAORDERPARAM_AV2(MYNODE), PTSTEPS
            CLOSE(1980)
            ELSE IF (ORGYT) THEN
            FILENAME4="T.RgyrAv.EintAv.RgyrAv2.Eint2.Steps."//trim(adjustl(istr))
            OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
            WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), RGYR_AV(MYNODE), 
     1              EINT_AV(MYNODE), RGYR_AV2(MYNODE),  
     2              EINT_AV2(MYNODE), PTSTEPS
            CLOSE(1980)
         ENDIF
      ENDIF

      FILENAME5="profile_E."//trim(adjustl(istr)) 
      OPEN(UNIT=1981,FILE=FILENAME5, STATUS="unknown", form="formatted")
      DO K=1,NHISTE
         WRITE(1981,'(2G20.10)') PTEMIN+(K-1)*DHISTE,NHISTQE(K,MYNODE)
      ENDDO
      CLOSE(1981)

      IF (BINARY) THEN
         FILENAME6="profile_Q4."//trim(adjustl(istr))
         FILENAME7="profile_Q6."//trim(adjustl(istr))
         OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
         OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
         DO K=1,NHIST
            WRITE(1982,'(2G20.10)') (K-1)*DQ4,NHISTQ4(K,MYNODE)
            WRITE(1983,'(2G20.10)') (K-1)*DQ6,NHISTQ6(K,MYNODE)
         ENDDO
         CLOSE(1982)
         CLOSE(1983)
       ELSE IF (CHRMMT) THEN
         IF (ODIHET) THEN
            FILENAME6="profile_dihe."//trim(adjustl(istr))
            FILENAME7="profile_sasa."//trim(adjustl(istr))
            OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
            OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
            DO K=1,NHIST
               WRITE(1982,'(2G20.10)') diheorderparam_min+(K-1)*Ddihe,NHISTQ4(K,MYNODE)
               WRITE(1983,'(2G20.10)') SASAorderparam_min+(K-1)*Dsasa,NHISTQ6(K,MYNODE)
            ENDDO
            CLOSE(1982)
            CLOSE(1983)
         ELSE IF (ORGYT) THEN
            FILENAME6="profile_rgyr."//trim(adjustl(istr))
            FILENAME7="profile_eint."//trim(adjustl(istr))
            OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
            OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
            DO K=1,NHIST
               WRITE(1982,'(2G20.10)') rgyr_min+(K-1)*Drgyr,NHISTQ4(K,MYNODE)
               WRITE(1983,'(2G20.10)') eint_min+(K-1)*Deint,NHISTQ6(K,MYNODE)
            ENDDO
            CLOSE(1982)
            CLOSE(1983)
            ENDIF
       ENDIF

      WRITE (ISTR, '(i10)') MYNODE+1
      FILENAME8="I.Vmin.Vmax."//trim(adjustl(istr)) 
      OPEN(UNIT=1984,FILE=FILENAME8, STATUS="unknown", form="formatted")
      WRITE(1984,'(3G20.10)') MYNODE, VMIN(MYNODE),VMAX(MYNODE)
      CLOSE(1984)

      WRITE(MYUNIT, '(G20.10,A,G20.10)') NACCEPTPT(MYNODE), ' steps accepted out of ', NSTAT
      WRITE(MYUNIT, '(A,G20.10)') 'Final stepsize ', STEP(MYNODE+1)
      WRITE(MYUNIT, '(G20.10,A,G20.10)') NEACCEPT, ' exchanges accepted out of ', NTOT

      DO K=1,NATOMS
         Q(1,K)=X(K)
         Q(2,K)=Y(K)
         Q(3,K)=Z(K)
      ENDDO
      IF (PERIODIC) THEN
         CALL QORDER_BLJ(Q,Q4,Q6)
      ELSE
         CALL QORDER_LJ(Q,Q4,Q6)
      ENDIF

      IF (.NOT.BINSTRUCTURES) THEN
         WRITE (ISTR, '(i10)') MYNODE+1
         FILENAME5="FreeEnStatHistO1.HistO2."//trim(adjustl(istr)) 
         OPEN(UNIT=747+MYNODE,FILE=FILENAME5, STATUS="unknown", form="formatted")
            DO J=1,NHIST
               DO I=1,NHIST
                   WRITE(747+MYNODE,'(G20.10)') NHISTALLQ(i,j,MYNODE)
               ENDDO
            ENDDO
         CLOSE(747+MYNODE)
         WRITE (ISTR, '(i10)') MYNODE+1
         FILENAME12="FreeEnStatHistE.HistO1.HistO2."//trim(adjustl(istr))
         OPEN(UNIT=447+MYNODE,FILE=FILENAME12, STATUS="unknown", form="formatted")
         DO K=1,NHISTE
            DO J=1,NHIST
               DO I=1,NHIST
                   IF (NHISTGRAND(k,i,j,MYNODE).ne.0) then
                   WRITE(447+MYNODE,'(4I10)') k,i,j,NHISTGRAND(k,i,j,MYNODE)
                   endif
               ENDDO
            ENDDO
         ENDDO
         CLOSE(447+MYNODE)
      ENDIF


      WRITE (ISTR, '(i10)') MYNODE 
      FILENAME101="Visits.his."//trim(adjustl(istr))
      OPEN(UNIT=1986+MYNODE,FILE=FILENAME101, STATUS="unknown", form="formatted")
      DO k=1, HBINS
           WRITE(1986+MYNODE, '(2G20.10)') BinLabel(k), Visits(k,MYNODE)
      ENDDO
      CLOSE(1986+MYNODE)


      WRITE (ISTR, '(i10)') MYNODE 
      FILENAME102="Distance_av.his."//trim(adjustl(istr))
      OPEN(UNIT=1987+MYNODE,FILE=FILENAME102, STATUS="unknown", form="formatted")
      DO k=1, HBINS
           IF (Visits(k,MYNODE).EQ.0) THEN
              WRITE(1987+MYNODE, '(2G20.10)') BinLabel(k), 0.0D0
           ELSE
              WRITE(1987+MYNODE, '(2G20.10)') BinLabel(k), Distance(k,MYNODE)/Visits(k,MYNODE)
           ENDIF
      ENDDO
      CLOSE(1987+MYNODE)

      WRITE (ISTR, '(i10)') MYNODE
      FILENAME103="MinimaNumber.his."//trim(adjustl(istr))
      OPEN(UNIT=1988+MYNODE,FILE=FILENAME103, STATUS="unknown", form="formatted")
      DO k=1, HBINS
           WRITE(1988+MYNODE, '(2G20.10)') BinLabel(k), MinimaNumber(k,MYNODE)
      ENDDO
      CLOSE(1985)

      IF (HISTRESTART) THEN
         WRITE (ISTR, '(i10)') MYNODE+1
         FILENAME9="xyz_out."//trim(adjustl(istr))
      ELSE
         WRITE (ISTR, '(i10)') MYNODE+1
         FILENAME9="xyz_out."//trim(adjustl(istr))
      ENDIF
      OPEN(UNIT=1985+MYNODE,FILE=FILENAME9, STATUS="unknown", form="formatted")
      WRITE(1985+MYNODE,*) '#E=', E
      DO I=1,NATOMS
         WRITE(1985+MYNODE,*) X(I),Y(I),Z(I)
      ENDDO
      CLOSE(1985+MYNODE)
      CLOSE(11+MYNODE)



#else
      RETURN
#endif
      END SUBROUTINE PTBASINSAMPLING

                !CALL FLUSH(MYUNIT)
                !RETURN
