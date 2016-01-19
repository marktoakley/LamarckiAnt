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


      INCLUDE 'mpif.h'
       

      
      INTEGER :: NHIST=100, NHISTE=2000, IACCEPT(0:NTRAJ-1), MPIERR, I,J,K, IS(MPI_STATUS_SIZE)
      INTEGER NHISTQ4(NHIST,0:NTRAJ-1), NDUMMY, NTOT, NH, IMESG, IQE, IQ4, IQ6, 
     1        NHISTQ6(NHIST,0:NTRAJ-1), NHISTQE(NHISTE, 0:NTRAJ-1), MSTEP, NSTAT, IENR, 
     2        NACCEPTPT(0:NTRAJ-1), NOUT(0:NTRAJ-1), ITRAJ, ITRAJO,NEACCEPT, RNDSEED, NUPDATE,
     3        CONVERGED,LBFGS_ITERATIONS, JD
      REAL(8) V(NATOMS), VO(NATOMS), TEMPTRAJ(0:NTRAJ-1), H(0:NTRAJ-1), BETA(0:NTRAJ-1), 
     1        EAV(0:NTRAJ-1), EAV2(0:NTRAJ-1), Q(3,NATOMS), Q4AV(0:NTRAJ-1), Q4AV2(0:NTRAJ-1), 
     2        Q6AV(0:NTRAJ-1), Q6AV2(0:NTRAJ-1),VMIN(0:NTRAJ-1), VMAX(0:NTRAJ-1), VENR(NENRPER), 
     3        HINIT(0:NTRAJ-1), EPSAB, EPSBB, SIGAB, SIGBB, X(NATOMS), Y(NATOMS), Z(NATOMS), 
     4        CTE, T, VINIT,VFINAL, POTEL, GRAD(3*NATOMS), Q4, Q6, RANDOM, DPRAND, Q4MAX, Q6MAX , 
     5	      DQ4, DQ6, DHISTE, ENUL, XO(NATOMS), YO(NATOMS), ZO(NATOMS), DDX, DDY, DDZ, DE, 
     6        W, WCOMP, WAC, E, ER, DBETA, DELTA, CV, FQ4, FQ6, RMAX, DDXN, DDYN, DDZN, R2, 
     7        DUMMY, BINLABEL(HBINS), BININDEX

      CHARACTER (LEN =256)  FILENAME, FILENAME2,FILENAME3,FILENAME4,FILENAME5,FILENAME6,
     1                      FILENAME7,FILENAME8,FILENAME9, FILENAME10 
      CHARACTER (LEN= 10)  ISTR
      LOGICAL EXCHANGE, EXCHANGEACCEPT, FITS

      COMMON /MYPOT/ POTEL
     

      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NDUMMY,MPIERR)
      IF ((NDUMMY.NE.NPAR).OR.(NDUMMY.NE.NTRAJ)) THEN
         PRINT *, 'Number of temperature trajectories does not correspond to the number  
     1           of processors. Stop.'
         RETURN
      ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)

      WRITE (ISTR, '(i10)') MYNODE
      FILENAME="output."//trim(adjustl(istr))
      OPEN(UNIT=1979,FILE=FILENAME, STATUS="unknown", form="formatted")
      WRITE(1979, '(A,I10,A,I10,A)') "Processor", mynode, " of", NPAR, " speaking:"
      WRITE(1979, '(A,I10)') 'Number of atoms', natoms
      IF (PERIODIC) THEN
         WRITE(1979, '(A,6G20.10)') 'Binary data', ntypea, epsab, epsbb, sigab, sigbb, cutoff
         WRITE(1979, '(A,3G20.10)') 'Box data', boxlx, boxly, boxlz
      ELSEIF(CHRMMT) THEN
         WRITE(1979, '(A)') 'CHARRRMMM'
      ELSE
         WRITE(1979, '(A,G20.10)') 'Radius**2', radius
      ENDIF
      CLOSE(1979)

      DO I=1, HBINS
         BINLABEL(I)=HISTMIN + HISTINT*(I-0.5)
      ENDDO


      ITRAJ=MYNODE
      NEACCEPT=0

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

      OPEN(UNIT=1979,FILE=FILENAME, STATUS="old", form="formatted", position="append")
      WRITE(1979, '(A,3G20.10)') 'xyz', x(1), y(1), z(1) 
      CLOSE(1979)
 
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
      ENDDO

      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
      VINIT=POTEL

      IF (PERIODIC) THEN
         CALL QORDER_BLJ(Q,Q4,Q6)
      ELSE
         CALL QORDER_LJ(Q,Q4,Q6)
      ENDIF

      OPEN(UNIT=1979,FILE=FILENAME, STATUS="old", form="formatted", position="append")
      WRITE(1979, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NTRAJ-1)
      WRITE(1979, '(A, G20.10)') 'This temperature trajectory=', TEMPTRAJ(MYNODE)
      WRITE(1979, '(A, G20.10)') 'Starting E=', VINIT
      WRITE(1979, '(A, 2G20.10)') 'Starting Q4, Q6=', Q4, Q6
      CLOSE(1979)

      ! Initialisation complete

      RNDSEED=2002+MYNODE
      CALL SDPRND(RNDSEED)
      RANDOM=DPRAND()
      OPEN(UNIT=1979,FILE=FILENAME, STATUS="old", form="formatted", position="append")
      WRITE(1979, '(A, G20.10)') 'Starting random number=', RANDOM
      CLOSE(1979)

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
      DO I=1,NHIST
         NHISTQ4(I,MYNODE)=0
         NHISTQ6(I,MYNODE)=0
      ENDDO
      DO I=1,NHISTE
         NHISTQE(I,MYNODE)=0
      ENDDO

      Q4MAX=0.1
      Q6MAX=0.5 ! will have to change as is rather system-specific
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
      

      FILENAME2="config."//trim(adjustl(istr))
      OPEN(UNIT=11+MYNODE,FILE=FILENAME2, form="formatted")

      NUPDATE=100
      MSTEP=NATOMS ! is it a general rule?
      NSTAT=MSTEP*(PTSTEPS+NEQUIL)

      DO I=1, PTSTEPS+NEQUIL
         DO J=1,MSTEP
         ! Propagate

            DO K=1, NATOMS
               XO(K)=X(K)
               YO(K)=Y(K)
               ZO(K)=Z(K)
            ENDDO
7777        RANDOM=DPRAND()
            NH=INT(1+NATOMS*RANDOM) ! checked that it covers the whole range
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


            VFINAL=0.D0
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
            VFINAL=POTEL

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
               VFINAL=VINIT
               IF (.NOT.FITS) NOUT(MYNODE)=NOUT(MYNODE)+1
            ELSE
               NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1
               IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
            ENDIF
            
            IF (I.GT.NEQUIL) THEN
               VMIN(MYNODE)=MIN(VMIN(MYNODE),VINIT)
               VMAX(MYNODE)=MAX(VMAX(MYNODE), VINIT)
            ENDIF
!Quenching part
!            CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS, &
!            & DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1))
!            IF (POTEL.LT.ENUL) then ! we have found a lower configuration
!               FILENAME10="lowest."//trim(adjustl(istr))
!               OPEN(UNIT=555+MYNODE,FILE=FILENAME10, STATUS="unknown", form="formatted")
!               WRITE(555+MYNODE,'(A,G20.10,A,I10,A,G20.10)') 'Energy', POTEL, 'steps ',i,'Temp', &
!               & TEMPTRAJ(MYNODE)
!               WRITE(555+MYNODE,30) (COORDS(JD,MYNODE+1),JD=1,3*NATOMS)
!30             FORMAT('LA ',3F20.10)
!               CALL FLUSH(555+MYNODE)
!               CLOSE(555+MYNODE)
!               ENUL=POTEL
!            ENDIF

            
            VINIT=VFINAL
            BININDEX=Energy2Index(VFINAL, BINLABEL(1))
            OPEN(UNIT=1979,FILE=FILENAME, STATUS="old", form="formatted", position="append")
            WRITE(1979, '(4G20.10)') I, J, BININDEX, VFINAL
            CLOSE(1979)
         ENDDO
       

         IF (MOD(I,NUPDATE).EQ.0) THEN ! update mc steps
            WAC=1.0*IACCEPT(MYNODE)/NUPDATE/MSTEP
            IF (WAC.LT.0.4) THEN
               H(MYNODE)=H(MYNODE)*0.9
            ENDIF
            IF (WAC.GT.0.6) THEN
               H(MYNODE)=H(MYNODE)*1.1
            ENDIF
            IACCEPT(MYNODE)=0
         ENDIF
         
         E=VFINAL
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

9999       continue

           VINIT=E

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
              ELSE
                 CALL QORDER_LJ(Q,Q4,Q6)
              ENDIF
              Q4AV(MYNODE)=Q4AV(MYNODE)+Q4
              Q4AV2(MYNODE)=Q4AV2(MYNODE)+Q4**2
              Q6AV(MYNODE)=Q6AV(MYNODE)+Q6
              Q6AV2(MYNODE)=Q6AV2(MYNODE)+Q6**2
              IQ4=INT(Q4/DQ4+1)
              IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                 NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
              ENDIF
              IQ6=INT(Q6/DQ6+1)             
              IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                 NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
              ENDIF
              IENR=IENR+1
              IF (IENR.EQ.NENRPER) THEN
                 IENR=0
                 CALL FLUSH(11+MYNODE)
                 WRITE(11+MYNODE,*) '#E=',E
                 DO K=1,NATOMS
                    WRITE(11+MYNODE,*) X(K),Y(K),Z(K)
                 ENDDO
               ENDIF 
           ENDIF                   
      ENDDO

      CALL MPI_FINALIZE(MPIERR)
 
      ! computing the averages
      EAV(MYNODE)=EAV(MYNODE)/PTSTEPS
      EAV2(MYNODE)=EAV2(MYNODE)/PTSTEPS
      CV=(EAV2(MYNODE)-EAV(MYNODE)**2)*BETA(MYNODE)**2

      FILENAME3="T.Ev.Cv.Ev2.Steps."//trim(adjustl(istr))
      OPEN(UNIT=41+MYNODE,FILE=FILENAME3, STATUS="unknown", form="formatted")
      WRITE(41+MYNODE,'(5G20.10)') TEMPTRAJ(MYNODE),EAV(MYNODE),CV, EAV2(MYNODE),PTSTEPS
      CALL FLUSH(41+MYNODE)
      CLOSE(41+MYNODE)

      Q4AV(MYNODE)=Q4AV(MYNODE)/PTSTEPS
      Q6AV(MYNODE)=Q6AV(MYNODE)/PTSTEPS
      Q4AV2(MYNODE)=Q4AV2(MYNODE)/PTSTEPS
      Q6AV2(MYNODE)=Q6AV2(MYNODE)/PTSTEPS
      FQ4=SQRT(Q4AV2(MYNODE)-Q4AV(MYNODE)**2)
      FQ6=SQRT(Q6AV2(MYNODE)-Q6AV(MYNODE)**2) 
      FILENAME4="T.Q4Av.Q6Av.Q4Av2.Q6Av2.Steps."//trim(adjustl(istr))
      OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
      WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), Q4AV(MYNODE), Q6AV(MYNODE), Q4AV2(MYNODE),  
     1                        Q6AV2(MYNODE), PTSTEPS
      CLOSE(1980)

      FILENAME5="profile_E."//trim(adjustl(istr)) 
      OPEN(UNIT=1981,FILE=FILENAME5, STATUS="unknown", form="formatted")
      DO K=1,NHISTE
         WRITE(1981,'(2G20.10)') PTEMIN+(K-1)*DHISTE,NHISTQE(K,MYNODE)
      ENDDO
      CLOSE(1981)

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

      FILENAME8="I.Vmin.Vmax."//trim(adjustl(istr)) 
      OPEN(UNIT=1984,FILE=FILENAME8, STATUS="unknown", form="formatted")
      WRITE(1984,'(3G20.10)') MYNODE, VMIN(MYNODE),VMAX(MYNODE)
      CLOSE(1984)

      OPEN(UNIT=1979,FILE=FILENAME, STATUS="old", form="formatted", position="append")
      WRITE(1979, '(G20.10,A,G20.10)') NACCEPTPT(MYNODE), ' steps accepted out of ', NSTAT
      WRITE(1979, '(G20.10,A,G20.10)') NEACCEPT, ' exchanges accepted out of ', NTOT
      CLOSE(1979)

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



      FILENAME9="xyz_out."//trim(adjustl(istr))
      OPEN(UNIT=1985,FILE=FILENAME9, STATUS="unknown", form="formatted")
      WRITE(1985, '(G20.10)') NTYPEA 
      WRITE(1985, '(4G20.10)') epsab, epsbb, sigab, sigbb
      WRITE(1985, '(3G20.10)') boxlx, boxly, boxlz
      WRITE(1985, '(G20.10)') cutoff
      DO I=1,NATOMS
         WRITE(1985,*) X(I),Y(I),Z(I)
      ENDDO
      WRITE(1985,*) '#E=', E
      WRITE(1985,*) '#Q4=',Q4 
      WRITE(1985,*) '#Q6=', Q6
      CLOSE(1985)
      CLOSE(11+MYNODE)




      END SUBROUTINE PTBASINSAMPLING

