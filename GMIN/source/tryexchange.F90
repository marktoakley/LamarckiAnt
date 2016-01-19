!  GMIN: A program for finding global minima
!  Copyright (C) 1999- David J. Wales
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
SUBROUTINE TRYEXCHANGE(E,X,Y,Z,XO,YO,ZO,VOLD,EXCHANGEACCEPT,JLOW, &
  &                    VNEW,GRAD,VNEWSAVE,VMINOLD,VMINNEW,BETA,ITRAJ,ITRAJO,NTOT, &
  &                    LBFGS_ITERATIONS,NEACCEPT,LBFGS_ITERATIONSO,QV,XDUMMY,PEINT,NCHOSEN,LASTEXDOWN, MCSTEP)
USE COMMONS
IMPLICIT NONE
DOUBLE PRECISION RANDOM, DPRAND, E, X(NATOMS), Y(NATOMS), Z(NATOMS), XO(NATOMS), YO(NATOMS), ZO(NATOMS), &
  &              POTEL, VOLD, VNEW, GRAD(3*NATOMS), &
  &              DELTA, W, ER, DBETA, VNEWSAVE, VMINOLD, VMINNEW, BETA(0:NPAR-1), LESAVE, PERTCOORDS(3*NATOMS), &
  &              PERTCOORDSO(3*NATOMS), DUMMY, NEWE, QV(NENRPER), XDUMMY, PEINT, OLDPE
INTEGER KAPPA
DOUBLE PRECISION, intent(in)    :: MCSTEP ! the monte carlo step. Used for printing 
INTEGER, INTENT(INOUT)    :: NEACCEPT, NTOT
!DOUBLE PRECISION    :: MCSTEP ! the monte carlo step. Used for printing 
LOGICAL, intent(OUT) :: EXCHANGEACCEPT
COMMON /MYPOT/ POTEL
INTEGER JLOW, MPIERR, ITRAJ, ITRAJO, K, I, LBFGS_ITERATIONS, IMESG, &
  &     LBFGS_ITERATIONSO, J2, NDUMMY, J1, NCHOSEN, MINCHOSEN
DOUBLE PRECISION LASTEXDOWN(0:NPAR-1)
double precision mycoords(3*natoms)

!ab2111>
DOUBLE PRECISION POTEL_HSA,CONFIG(3*NATOMS)
DOUBLE PRECISION WORKF,WORKR,WORK,Pdummy
INTEGER WELLLOW,WELLphys,ITRAJLIST(NPAR),ITRAJLISTO(NPAR),HEADNODE
INTEGER GETUNIT
LOGICAL WRITELIST, RESERVOIR_SERIAL
LOGICAL, SAVE :: FIRST = .TRUE.

#ifdef MPI
INCLUDE 'mpif.h'
INTEGER IS(MPI_STATUS_SIZE)


EXCHANGEACCEPT=.FALSE.

IF (PTINTERVAL) THEN
   !js850> this is a temporary implementation, eventually it should be 
   !done in a more cohesive way with all the other exchange options
   IF (INT(MOD(MCSTEP+1,EXCHINT)) .NE. 0) RETURN
ENDIF

IF (PTEX_INDEP) THEN
   MYCOORDS(1:3*NATOMS-2:3) = X(:)
   MYCOORDS(2:3*NATOMS-1:3) = Y(:)
   MYCOORDS(3:3*NATOMS-0:3) = Z(:)
   CALL TRYEXCHANGE_INDEPENDENCE(MYNODE, E, MYCOORDS, NPAR, &
         beta, EXCHANGEACCEPT, MCSTEP, myunit)
   NTOT=NTOT+1 
   if (exchangeaccept) then
      X(:) = MYCOORDS(1:3*NATOMS-2:3)
      Y(:) = MYCOORDS(2:3*NATOMS-1:3)
      Z(:) = MYCOORDS(3:3*NATOMS-0:3)
      VNEW = E
      VOLD = E
      NEACCEPT = NEACCEPT + 1
   endif
   return
ENDIF

! attempt reservoir exchange move
IF (RESERVOIRT.AND.(MYNODE.EQ.USERES)) THEN
   RANDOM=DPRAND()
   IF(RANDOM.LT.RES_PSWAP) THEN
      CALL RESERVOIR_EXCHANGE(X,Y,Z,XO,YO,ZO,VOLD, &                                           
                              VNEW,GRAD,BETA, &
                              LBFGS_ITERATIONS,MCSTEP,ITRAJ)
   ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Decision time
! jdf43> if RESERVOIR, JLOW must be set on node rank USERES+1, so to keep things
!        simple, determine JLOW and broadcast on node rank USERES+1, regardless
!        of whether a reservoir is used.
!IF (MYNODE.EQ.(USERES+1)) THEN
!   RANDOM=DPRAND()
!   JLOW=(NPAR-1)*RANDOM
!   IF (.NOT. PTINTERVAL) THEN !js850> this should be temporary
!      RANDOM=DPRAND()
!      IF (RANDOM.GT.EXCHPROB) JLOW=-2
!   ENDIF
!ENDIF
!! Decides destination for RANDOM/INTERVAL and SINGLE/SETS exchanges 
CALL BHPT_GET_DESTINATION(MCSTEP,JLOW)
IF (JLOW==-1) JLOW=-2
IF (MYNODE==(JLOW-1)) JLOW=MYNODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! No exchange for reservoirs unless this is the top replica using the reservoir.
! Also no exchange attempted within the equilibration period specified by EXEQ.
!
!IF (RESERVOIRT) THEN
!   IF (JLOW.LT.USERES) JLOW=-2
!ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Node 0 broadcasts, the others receive the value of JLOW.
!
!CALL MPI_BCAST(JLOW,1,MPI_INTEGER,USERES+1,MPI_COMM_WORLD,MPIERR)
!
! JLOW is set to -2 for no exchange. 
! We try to exchange with trajectory JLOW+1. JLOW starts from 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!DO J1=1,NPAR
!ITRAJLIST(J1) = 0
!ENDDO
HEADNODE=0

IF(JLOW.EQ.-2) THEN
!   CALL MPI_GATHER(ITRAJ,1,MPI_INTEGER,ITRAJLIST,1,MPI_INTEGER,HEADNODE,MPI_COMM_WORLD,MPIERR)
!   RETURN
ENDIF

! disable deprecated reservoir steps below (where reservoir has its own node)
RESERVOIR_SERIAL = .FALSE.

IF (MYNODE.EQ.(JLOW+1)) THEN
   IF (RESERVOIR_SERIAL.AND.((MYNODE-1).EQ.USERES)) THEN
      ! store local copy of configuration
      !CONFIG(:) = COORDS(:,MYNODE+1)

      DO J1=1,NATOMS
         CONFIG(3*(J1-1)+1) = X(J1)
         CONFIG(3*(J1-1)+2) = Y(J1)
         CONFIG(3*(J1-1)+3) = Z(J1)
         !WRITE(*,'(A,I6,G20.10)') "try> config ", 3*(J1-1)+1,CONFIG(3*(J1-1)+1) 
      ENDDO

      CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.TRUE.,.FALSE.)
      !WRITE(*,'(A,G20.10)') "try> potential ", POTEL

      CALL RESERVOIR_CALC_EHSA(X,Y,Z,WELLLOW,POTEL,POTEL_HSA,BETA(:),LBFGS_ITERATIONS)
      WORK = BETA(MYNODE-1)*POTEL_HSA - BETA(MYNODE)*POTEL
      CALL MPI_SEND(WELLLOW,1,MPI_INTEGER,JLOW,0,MPI_COMM_WORLD,MPIERR)

   ELSEIF (RESERVOIR_SERIAL.AND.(MYNODE.EQ.USERES)) THEN
      CALL RESERVOIR_SAMPLE_HSA(X,Y,Z,WELLLOW,POTEL,POTEL_HSA,BETA,LBFGS_ITERATIONS)

      DO J1=1,NATOMS
         CONFIG(3*(J1-1)+1) = X(J1)
         CONFIG(3*(J1-1)+2) = Y(J1)
         CONFIG(3*(J1-1)+3) = Z(J1)
         !WRITE(*,'(A,I,G)') "try> config ", 3*(J1-1)+1,CONFIG(3*(J1-1)+1) 
      ENDDO

      CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.TRUE.,.FALSE.)

      WORK = BETA(MYNODE-1)*POTEL - BETA(MYNODE)*POTEL_HSA
      CALL MPI_SEND(WELLLOW,1,MPI_INTEGER,JLOW,0,MPI_COMM_WORLD,MPIERR)

   ! Standard work calculation for parallel tempering and non-reservoir replica exchanges
   ELSEIF ((RESERVOIR_SERIAL.EQV..FALSE.).OR.(JLOW.GT.USERES)) THEN
      DO J1=1,NATOMS
         CONFIG(3*(J1-1)+1) = X(J1)
         CONFIG(3*(J1-1)+2) = Y(J1)
         CONFIG(3*(J1-1)+3) = Z(J1)
         !WRITE(*,'(A,I,G)') "try> config ", 3*(J1-1)+1,CONFIG(3*(J1-1)+1) 
      ENDDO
      CALL POTENTIAL(CONFIG(:),GRAD,E,.TRUE.,.FALSE.)
      WORK = BETA(MYNODE-1)*VNEW - BETA(MYNODE)*VNEW
      !WORK = BETA(MYNODE-1)*E - BETA(MYNODE)*E
   ENDIF
      

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Send work value to JLOW for acception / rejection criterion
   !
   ! Instead of energies, more general "work" values are now passed between replicas, which
   ! facilitates hamiltonian exchange as well as parallel tempering.
   !  Work values are defined as 
   !     WORK = BETA(MYNODE+1) H_{MYNODE+1}(x) - BETA(MYNODE) H_MYNODE(x)
   !   where H_i is energy function of replica i
   !
   CALL MPI_SEND(WORK,1,MPI_DOUBLE_PRECISION,JLOW,0,MPI_COMM_WORLD,MPIERR)
  !CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,JLOW,0,MPI_COMM_WORLD,MPIERR)
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDIF

IF (MYNODE.EQ.JLOW) THEN
   IF (RESERVOIR_SERIAL.AND.(MYNODE.EQ.USERES)) THEN
      CALL RESERVOIR_SAMPLE_HSA(X,Y,Z,WELLLOW,POTEL,POTEL_HSA,BETA,LBFGS_ITERATIONS)
      !DO J1=1,NATOMS
      !   CONFIG(3*(J1-1)+1) = X(J1)
      !   CONFIG(3*(J1-1)+2) = Y(J1)
      !   CONFIG(3*(J1-1)+3) = Z(J1)
      !ENDDO
      !CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.FALSE.,.FALSE.)
      CALL MPI_RECV(WELLphys,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,IS,MPIERR)
      WORKF = BETA(MYNODE+1)*POTEL - BETA(MYNODE)*POTEL_HSA
   ! Standard work calculation for non-reservoir replicas 
   ELSEIF (RESERVOIR_SERIAL.AND.((MYNODE+1).EQ.USERES)) THEN
      ! store local copy of configuration
      !CONFIG(:) = COORDS(:,MYNODE+1)

      DO J1=1,NATOMS
         CONFIG(3*(J1-1)+1) = X(J1)
         CONFIG(3*(J1-1)+2) = Y(J1)
         CONFIG(3*(J1-1)+3) = Z(J1)
         !WRITE(*,'(A,I,G)') "try> config ", 3*(J1-1)+1,CONFIG(3*(J1-1)+1) 
      ENDDO

      CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.TRUE.,.FALSE.)
      !WRITE(*,'(A,G)') "try> potential ", POTEL

      CALL RESERVOIR_CALC_EHSA(X,Y,Z,WELLLOW,POTEL,POTEL_HSA,BETA(:),LBFGS_ITERATIONS)
      CALL MPI_RECV(WELLphys,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,IS,MPIERR)
      WORKF = BETA(MYNODE+1)*POTEL_HSA - BETA(MYNODE)*POTEL

   ELSE
      !ab2111> debug:
      IF (DEBUG) THEN
         DO J1=1,NATOMS
            CONFIG(3*(J1-1)+1) = X(J1)
            CONFIG(3*(J1-1)+2) = Y(J1)
            CONFIG(3*(J1-1)+3) = Z(J1)
         ENDDO
         CALL POTENTIAL(CONFIG(:),GRAD,Pdummy,.TRUE.,.FALSE.)
      !WRITE(MYUNIT,'(A,5F20.10)') "sanity check: ",E,Pdummy,E-Pdummy,VNEW,VNEW-Pdummy
      ENDIF
       !WORKF = BETA(MYNODE+1)*E - BETA(MYNODE)*E
       WORKF = BETA(MYNODE+1)*VNEW - BETA(MYNODE)*VNEW
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Receive work value from JLOW+1 for acception / rejection criterion
   !
   CALL MPI_RECV(WORKR,1,MPI_DOUBLE_PRECISION,JLOW+1,0,MPI_COMM_WORLD,IS,MPIERR)
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Acceptance / Rejection decision
   !
   ! Acceptance criterion, written in terms of Work values:
   ! P_acc = min(1,e^{-WORKF-WORKR})
   !
   ! where WORKF+WORKR = BETA(i+1) H_{i+1}(x) - BETA(i) H_i(x)
   !                    + BETA(i) H_{i}(y) - BETA(i+1) H_{i+1}(y)
 
   WORK = WORKF + WORKR
   W=MIN(1.0D0,DEXP(-WORK))
   ! check if NaN
   IF(WORK.EQ.(WORK+1.0)) W=0.0

   !CALL MPI_RECV(ER,1,MPI_DOUBLE_PRECISION,JLOW+1,0,MPI_COMM_WORLD,IS,MPIERR)
   !DBETA=BETA(JLOW)-BETA(JLOW+1)
   !DELTA=E-ER
   !W=MIN(1.0D0,DEXP(DELTA*DBETA))

   NTOT=NTOT+1 
   RANDOM=DPRAND()
   LESAVE=E

   !ab2111> debug: (nores implementation with reservoir keyword turned on)
   !IF(JLOW.EQ.0) THEN
   !   W=0.0
   !ENDIF
   !W=1.0

   IF (MOD(MCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) THEN
      WRITE(MYUNIT,'(A,4G14.7,2I6,F20.1,I6)') 'tryexchange> work:  ',WORKF,WORKR,WORK,W,WELLLOW,WELLphys,MCSTEP,ITRAJ
   ENDIF

   !IF (MOD(MCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) THEN
   !   WRITE(MYUNIT,'(A,4G14.7,2I,F20.1)') 'tryexchange> work:',WORKF,WORKR,WORK,W,WELLLOW,WELLphys,MCSTEP
   !ENDIF
   !IF (RESERVOIR_SERIAL) THEN 
   !   WRITE(MYUNIT,'(A,4G20.10,2I6)') 'te> work:',WORKF,WORKR,WORK,W,WELLLOW,WELLphys
   !ENDIF

   IF (W.GT.RANDOM) THEN
      EXCHANGEACCEPT=.TRUE.
      IMESG=1
      
      CALL MPI_SEND(IMESG,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(VMINNEW,1,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(VNEW,1,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(LBFGS_ITERATIONS,1,MPI_INTEGER,JLOW+1,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_RECV(ITRAJ,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(X,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(Y,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(Z,NATOMS,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(VMINNEW,1,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(VNEW,1,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(LBFGS_ITERATIONS,1,MPI_DOUBLE_PRECISION,JLOW+1,1,MPI_COMM_WORLD,IS,MPIERR)
      ER=WORKR * 1./(BETA(JLOW)-BETA(JLOW+1))
      E=ER
      NEACCEPT=NEACCEPT+1
      LBFGS_ITERATIONSO=LBFGS_ITERATIONS
      VMINOLD=VMINNEW ! these variables are the same at this point in the program
      VOLD=VNEW
   ELSE
      EXCHANGEACCEPT=.FALSE.
      IMESG=0
      CALL MPI_SEND(IMESG,1,MPI_INTEGER,JLOW+1,0,MPI_COMM_WORLD,MPIERR)
   ENDIF
   IF (EXCHANGEACCEPT) THEN
      ER=WORKR * 1./(BETA(JLOW)-BETA(JLOW+1))
!     IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(MYUNIT,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') &
      IF(DEBUG) THEN
         WRITE(MYUNIT,'(A,I5,A,G16.6,A,I5,A,G16.6,A,F20.1)') &
           &               'tryexchange> Exchange from this replica ',JLOW,' energy ',LESAVE,&
           &                                           ' to replica ',JLOW+1,' energy ',ER, ' ACC STEP', MCSTEP
      ENDIF
      VNEWSAVE=VOLD
      IF (DEBUG) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'tryexchange> Instantaneous PE is now ',VOLD

         !DO K=1,NATOMS
         !   COORDS(3*(K-1)+1,MYNODE+1)=X(K)
         !   COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
         !   COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         !ENDDO
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
         WRITE(MYUNIT,'(2(A,G20.10))') 'tryexchange> Recalculated PE is ',POTEL,' VNEW=',VNEW
      ENDIF

   ELSE
      ! if Hamiltonian exchange (i.e. default RESERVOIR implementation), then
      ! print \Delta E = H_RESERVOIR - H_Physical for ER
      IF (BETA(JLOW).EQ.BETA(JLOW+1)) THEN
         ER=WORKR * 1./BETA(JLOW)
      ELSE
         ER=WORKR * 1./(BETA(JLOW)-BETA(JLOW+1))
      ENDIF
!     IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0.0D0) WRITE(MYUNIT,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') &
      IF (DEBUG) WRITE(MYUNIT,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') &
     &           'tryexchange> Exchange from replica ',JLOW,' energy ',E,' to replica ',&
     &                  JLOW+1,' energy ',ER,' REJ'
   ENDIF
ENDIF
CALL FLUSH(MYUNIT)
IF (MYNODE.EQ.(JLOW+1)) THEN
   CALL MPI_RECV(IMESG,1,MPI_INTEGER,JLOW,0,MPI_COMM_WORLD,IS,MPIERR)
   NTOT=NTOT+1
   ER=E
   IF (IMESG.EQ.1) THEN
!
!  Here we receive first, so we need to save variables that would otherwise be overwritten
!  before they are sent to the other replica. Or use XO, YO, ZO in the case of X, Y, Z.
!
      LBFGS_ITERATIONSO=LBFGS_ITERATIONS
      CALL MPI_RECV(ITRAJO,1,MPI_INTEGER,JLOW,0,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(XO,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(YO,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(ZO,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(VMINOLD,1,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(VOLD,1,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_RECV(LBFGS_ITERATIONS,1,MPI_INTEGER,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,JLOW,0,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,MPIERR)
      CALL MPI_SEND(VMINNEW,1,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_SEND(VNEW,1,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      CALL MPI_SEND(LBFGS_ITERATIONSO,1,MPI_DOUBLE_PRECISION,JLOW,1,MPI_COMM_WORLD,IS,MPIERR)
      DO K=1, NATOMS
         X(K)=XO(K)
         Y(K)=YO(K)
         Z(K)=ZO(K)
      ENDDO
      ITRAJ=ITRAJO
      NEACCEPT=NEACCEPT+1
      EXCHANGEACCEPT=.TRUE.
      LBFGS_ITERATIONSO=LBFGS_ITERATIONS
     !IF (MOD(MCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) THEN
     !   WRITE(MYUNIT,'(A,I5,A,G16.6,A,I5,A,G16.6,A,F20.1)') &
     !&            'tryexchange> Exchange from replica ',MYNODE-1,' energy ',VOLD,&
     !&                                                ' to this replica ',MYNODE,' energy ',VNEW, ' ACC STEP', MCSTEP
     !ENDIF
      VMINNEW=VMINOLD ! these variables are the same at this point in the program
      VNEW=VOLD
      ! RESERVOIR exchange: make sure VNEW = Hphysical(X,Y,Z), not Hreservoir(X,Y,Z)
      IF(RESERVOIR_SERIAL.AND.(JLOW.EQ.USERES)) THEN
         DO K=1,NATOMS
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, VNEW, .TRUE., .FALSE.)
         VOLD=VNEW
      ENDIF
      IF (DEBUG) THEN
         WRITE(MYUNIT,'(2(A,G20.10))') 'tryexchange> Instantaneous PE is now ',VOLD,' and VNEW=',VNEW
         DO K=1,NATOMS
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
         WRITE(MYUNIT,'(A,G20.10)') 'tryexchange> Recalculated PE is ',POTEL
      ENDIF

   ELSE
      EXCHANGEACCEPT=.FALSE.
!     IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0.0D0) WRITE(MYUNIT,'(A,I5,A,I5,A,G16.6,A)') &
      IF (DEBUG) WRITE(MYUNIT,'(A,I5,A,I5,A,G16.6,A)') &
     &          'tryexchange> Exchange from      replica ',MYNODE-1,&
     &          ' energy not available to this replica ',MYNODE,' energy ',VOLD,' REJ'
   ENDIF
ENDIF

!open permutations file if hasn't been done so
IF (FIRST.AND.(MYNODE.EQ.HEADNODE)) THEN
   FIRST = .FALSE.
   ! open the file if it hasn't been done yet
   PTEX_INDEP_UNIT = GETUNIT()
   OPEN(UNIT=PTEX_INDEP_UNIT,FILE="exchanges.permutations", STATUS="unknown", form="formatted")
   !warning, this unit is never closed.  I don't think it matters because it
   !should stay open until the program exits anyway
ENDIF

!ab2111> send exchange pairs to HEADNODE and write replica trajectory changes to permutations file
ITRAJLISTO(:)=ITRAJLIST(:)
!CALL MPI_GATHER(ITRAJ,1,MPI_INTEGER,ITRAJLIST,1,MPI_INTEGER,HEADNODE,MPI_COMM_WORLD,MPIERR)
IF (MOD(MCSTEP-1.0D0,1.0D0*PRTFRQ).EQ.0) THEN
   CALL MPI_GATHER(ITRAJ,1,MPI_INTEGER,ITRAJLIST,1,MPI_INTEGER,HEADNODE,MPI_COMM_WORLD,MPIERR)
   IF(MYNODE.EQ.HEADNODE) THEN
   !WRITELIST=.FALSE.
   !DO J1=1,NPAR
   !   IF (ITRAJLISTO(J1).NE.ITRAJLIST(J1)) THEN
   !      WRITELIST=.TRUE.
   !      EXIT 
   !   ENDIF
   !ENDDO
   !IF(WRITELIST) THEN
      WRITE(PTEX_INDEP_UNIT,'(F20.10,A)', advance="no") MCSTEP, " "
      DO J1=1,NPAR
         WRITE(PTEX_INDEP_UNIT,'(I3)',advance="no") ITRAJLIST(J1)
      ENDDO
      WRITE(PTEX_INDEP_UNIT,'(A)') ""
      CALL FLUSH(PTEX_INDEP_UNIT)
   ENDIF
ENDIF

!ab2111> make sure all copies of coordinates are consistent:
DO K=1,NATOMS
    COORDS(3*(K-1)+1,MYNODE+1)=X(K)
    COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
    COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
ENDDO

#else
      RETURN
#endif

END SUBROUTINE TRYEXCHANGE
