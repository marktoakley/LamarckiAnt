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

! 
! This version of tryexchange tries to redistribute all replicas according
! to the appropriate probability distribution.  In order to implement this
! the head node gets the energies from all the replicas, and does a short 
! Monte Carlo simulation.  This consists of doing N exchanges between randomly
! chosen replicas where N is large enough to reach equilibration.  See
! J. D. Chodera, and M. R. Shirts, JCP 2011, http://dx.doi.org/10.1063/1.3660669
! In that paper, they suggest that from their own experience N should be from
! nreps**3 to nreps**5.  
! 
! Then, once the new distribution of replicas is chosen, the actual exchanges
! are done via communication with the other nodes
! 

module tryexchange_monte_carlo_mod
!this module implements the mini Monte Carlo routine redistribute
!the replicas independently
implicit none
private

integer nreps
!integer ncount
!double precision time_mc, time_indep, time_wait


public :: tryexchange_monte_carlo

contains

subroutine get_pair(j1, j2)
   implicit none
   integer, intent(out) :: j1, j2
   double precision dprand
   j1 = 1
   j2 = j1
   do while (j1 .eq. j2) 
      !here we choose j1 and j2 randomly, but it might be more efficient to
      !exchange neighbors because the probability to accept a long distance
      !exchange will probably be very small.
      j1 = int(dprand() * nreps) + 1 ! 1 to nreps inclusive
      j2 = int(dprand() * nreps) + 1 ! 1 to nreps inclusive
   enddo
end subroutine get_pair

subroutine tryexchange_monte_carlo(nsteps, nreps_i, energies, beta, newpositions)
   implicit none
   integer, intent(in) :: nsteps, nreps_i
   double precision, intent(in) :: beta(nreps_i)
   double precision, intent(in) :: energies(nreps_i)
   integer, intent(out) :: newpositions(nreps_i)
   double precision e1, e2, beta1, beta2, w, random
   double precision DPRAND
   integer j1, j2, n, i1, i2, jtemp
   !
   !newpositions holds the permutation of the replicas.
   !newpositions(i) is the label of the replica which is currently at i.  E.g.
   !after the monte carlo run.  Replica newposition(i) should be moved to
   !position i
   !
   nreps = nreps_i

   !initialize newpositions
   do j1=1,nreps
      newpositions(j1) = j1
   enddo

   !do monte carlo loop
   do n=1,nsteps
      call get_pair(j1, j2)

      !try to exchange j1, j2
      !i1 and i2 are the labels of the replicas currently at j1 and
      !j2
      i1 = newpositions(j1) 
      i2 = newpositions(j2)
      e1 = energies(i1) !use the original configurations
      e2 = energies(i2)
      beta1 = beta(j1) !use the beta from the current position
      beta2 = beta(j2)
      W=MIN(1.0D0,DEXP((e1-e2)*(beta1-beta2)))
      RANDOM=DPRAND()
      !write(*,*) "trying exchange", j1, j2, e1, e2
      IF (W.GT.RANDOM) THEN
         !accept exchange
         !write(*,*) "   accepting"
         jtemp = newpositions(j1)
         newpositions(j1) = newpositions(j2)
         newpositions(j2) = jtemp
      ENDIF
   enddo

END subroutine TRYEXCHANGE_MONTE_CARLO
END MODULE TRYEXCHANGE_MONTE_CARLO_MOD

SUBROUTINE TRYEXCHANGE_INDEPENDENCE(MYNODE, MYENERGY, MYCOORDS, nreps, beta, EXCHANGEACCEPT, &
      MCSTEP, myunit)
#ifdef MPI
  USE COMMONS, only : natoms, ptex_indep_nsteps, PTEX_INDEP_UNIT, debug
  USE TRYEXCHANGE_MONTE_CARLO_MOD, ONLY : TRYEXCHANGE_MONTE_CARLO
#else
  USE COMMONS, only: natoms
#endif
IMPLICIT NONE
INTEGER, INTENT(IN) :: MYNODE, nreps, myunit
DOUBLE PRECISION, INTENT(INOUT) :: MYCOORDS(3*NATOMS), MYENERGY
DOUBLE PRECISION, INTENT(IN) :: beta(nreps), mcstep
LOGICAL, INTENT(OUT) :: EXCHANGEACCEPT
integer headnode, nmcsteps, node, tempinteger
integer nperm, permutation(nreps)
integer newpositions(nreps), oldpositions(nreps), newnode, oldnode
DOUBLE PRECISION ENERGIES(NREPS), newenergy
double precision allcoords(3*natoms, nreps) !this is large
double precision tempcoords(3*natoms)
integer j1
integer MPITAG
INTEGER MPIERR
INTEGER getunit
LOGICAL, SAVE :: FIRST = .TRUE.
double precision timestart, timeend, time2, time3
integer, save :: ncount
double precision, save :: time_indep, time_mc, time_wait
#ifdef MPI
INCLUDE 'mpif.h'
INTEGER MPISTATUS(MPI_STATUS_SIZE)


CALL MYCPU_TIME(TIMESTART)

MPITAG = 1
HEADNODE = 0
IF (PTEX_INDEP_NSTEPS  .GT. 0) THEN
   NMCSTEPS = PTEX_INDEP_NSTEPS
ELSE
   NMCSTEPS = NREPS**3 
ENDIF

EXCHANGEACCEPT = .FALSE.

!IF (MYNODE .EQ. HEADNODE) THEN
   !write(*,*) ""
   !write(*,*) "in tryexchange_independence"
!endif
!write(*,'(A,I5,A,4F20.10)') "I'm node", mynode, " old", myenergy, mycoords(3*natoms-2:3*natoms)

IF (MYNODE .EQ. HEADNODE) THEN
   !write(*,*) ""
   !write(*,*) "in tryexchange_independence"
   !
   !get energies from all nodes
   !
   ENERGIES(1) = MYENERGY
   !write(*,*) "getting energy from nodes"
   DO NODE=2,NREPS
      !write(*,*) "   getting energy from node", node-1
      CALL MPI_RECV(NEWENERGY, 1, MPI_DOUBLE_PRECISION, NODE-1, &
         MPITAG, MPI_COMM_WORLD, MPISTATUS, MPIERR)
      ENERGIES(NODE) = NEWENERGY
   ENDDO
   CALL MYCPU_TIME(TIME2)


   !
   !do monte carlo simulation
   !
   !write(*,*) "doing tryexchange monte carlo"
   CALL TRYEXCHANGE_MONTE_CARLO(NMCSTEPS, NREPS, ENERGIES, BETA, PERMUTATION)
   CALL MYCPU_TIME(TIME3)

   !
   !print the permutation
   !
   IF (FIRST) THEN
      FIRST = .FALSE.
      ncount = 0
      time_mc = 0. !this overwrites the first one... oh well
      time_indep = 0.
      ! open the file if it hasn't been done yet
      PTEX_INDEP_UNIT = GETUNIT()
      OPEN(UNIT=PTEX_INDEP_UNIT,FILE="exchanges.permutations", STATUS="unknown", form="formatted")
      !warning, this unit is never closed.  I don't think it matters because it
      !should stay open until the program exits anyway
   ENDIF
   !write to file exchanges.permutations
   write(PTEX_INDEP_UNIT,'(F20.1,A)', advance="no") mcstep, " "
   do j1=1,nreps
      write(PTEX_INDEP_UNIT,'(I5)', advance="no") permutation(j1)
   enddo
   write(PTEX_INDEP_UNIT,'(A)') ""
   CALL FLUSH(PTEX_INDEP_UNIT)

   !write(*,'(A,F20.1)') "tryexchange_independence> exchange permutation", mcstep
   !do j1=1,nreps
   !   write(*,'(I5)', advance="no") permutation(j1)
   !enddo
   !write(*,*) ""

   !write(myunit,'(A,F20.1)') "replica permutation", mcstep
   !do j1=1,nreps
   !   write(myunit,'(I5)', advance="no") permutation(j1)
   !enddo
   !write(myunit,*) ""

   !do j1=1,nreps
      !write(*,'(F20.10)', advance="no") energies(j1)
   !enddo
   !write(*,*) ""
   

   !
   !now we do exchanges:  Actually, it's a permutation, which could be thought of
   !as a series of exchanges, but it would be more complicated.  So we copy all the 
   !data to the head node and then redistribute it to the appropriate place
   !

   !
   !Get list of nodes that need to change and a list of where they are going
   !Also, tell all the nodes whether they need to do an exchange
   !
   !write(*,*) "telling the nodes whether they need to exchange"
   NPERM = 0
   DO J1=1,NREPS
      IF (PERMUTATION(J1) .NE. J1) THEN
         NPERM = NPERM + 1
         OLDPOSITIONS(NPERM) = PERMUTATION(J1)
         NEWPOSITIONS(NPERM) = J1
         TEMPINTEGER = 1
      ELSE
         TEMPINTEGER = 0
      ENDIF
      IF (J1-1 .NE. HEADNODE) THEN
         CALL MPI_SEND(TEMPINTEGER, 1, MPI_INTEGER, J1-1, &
            MPITAG, MPI_COMM_WORLD, MPIERR)
      ENDIF
   ENDDO

   !
   !copy all relevant data (energy, coords, ect) from the nodes to this, the head
   !node.  For PTMC we only need energy and coords.  We already have energy, so we
   !really just need coords
   !
   !write(*,*) "copying data from nodes"
   DO J1=1,NPERM
      OLDNODE = OLDPOSITIONS(J1) 
      IF (OLDNODE-1 .EQ. HEADNODE) THEN
         ALLCOORDS(1:3*NATOMS,OLDNODE) = MYCOORDS(:)
      ELSE
         CALL MPI_RECV(TEMPCOORDS, 3*NATOMS, MPI_DOUBLE_PRECISION, OLDNODE-1, &
            MPITAG, MPI_COMM_WORLD, MPISTATUS, MPIERR)
         ALLCOORDS(1:3*NATOMS,OLDNODE) = TEMPCOORDS(:)  !THIS DOUBLE COPY MAY NOT BE NECESSARY
      endif
   ENDDO

   !
   !redistribute the data (energy, coords, ect) to their new homes
   !
   !write(*,*) "sending data to nodes"
   DO J1=1,NPERM
      OLDNODE = OLDPOSITIONS(J1) 
      NEWNODE = NEWPOSITIONS(J1) 
      IF (NEWNODE-1 .EQ. HEADNODE) THEN
         MYENERGY = ENERGIES(OLDNODE)
         MYCOORDS(:) = ALLCOORDS(1:3*NATOMS, OLDNODE)
      ELSE
         !write(*,'(A,4I5)') "   sending to", newnode-1, oldnode-1, newnode, oldnode
         CALL MPI_SEND(ENERGIES(OLDNODE), 1, MPI_DOUBLE_PRECISION, NEWNODE-1, &
            MPITAG, MPI_COMM_WORLD, MPIERR)
         TEMPCOORDS(:) = ALLCOORDS(1:3*NATOMS, OLDNODE) !THIS DOUBLE COPY MAY NOT BE NECESSARY
         CALL MPI_SEND(TEMPCOORDS, 3*NATOMS, MPI_DOUBLE_PRECISION, NEWNODE-1, &
            MPITAG, MPI_COMM_WORLD, MPIERR)
      ENDIF
   ENDDO

   !
   !here do anything else that the head node needs to do 
   !
   IF ( PERMUTATION(1) .EQ. 1 ) THEN
      EXCHANGEACCEPT = .FALSE.
   else
      EXCHANGEACCEPT = .TRUE.
   ENDIF
   !write(*,*) "headnode done"
   !write(*,*) ""

   ! print some timeing info
   CALL MYCPU_TIME(TIMEEND)
   time_indep = time_indep + timeend - timestart
   time_mc = time_mc + time3 - time2
   time_wait = time_wait + time2 - timestart
   ncount = ncount + 1
   if (mod(ncount,10) .eq. 0) then
      write(myunit,'(3(A,G15.5),I15)') "tryexchange_independence> times: ", time_indep, " waiting ", time_wait, " in mc ", time_mc, ncount
   endif

ELSE
   !
   ! This is run by all nodes except the head node
   ! Here goes all the calls for communication with the head node
   !

   !
   ! return energy
   !
   !write(*,*) "i'm node", mynode, "sending energy"
   CALL MPI_SEND(MYENERGY, 1, MPI_DOUBLE_PRECISION, HEADNODE, &
      MPITAG, MPI_COMM_WORLD, MPIERR)

   !
   ! recieve from the head node if this node will be doing an exchange.
   !
   CALL MPI_RECV(TEMPINTEGER, 1, MPI_INTEGER, HEADNODE, &
      MPITAG, MPI_COMM_WORLD, MPISTATUS,MPIERR)
   EXCHANGEACCEPT = (.NOT. TEMPINTEGER .EQ. 0)
   !write(*,*) "I'm node", mynode, "EXCHANGEACCEPT", EXCHANGEACCEPT


   IF (EXCHANGEACCEPT) THEN
      !
      ! return energy, coords, etc
      !
      CALL MPI_SEND(MYCOORDS, 3*NATOMS, MPI_DOUBLE_PRECISION, HEADNODE, &
         MPITAG, MPI_COMM_WORLD, MPIERR)

      !
      ! recieve new coords, energy etc
      !
      !write(*,*) "I'm node", mynode, "recieving energy"
      CALL MPI_RECV(MYENERGY, 1, MPI_DOUBLE_PRECISION, HEADNODE, &
         MPITAG, MPI_COMM_WORLD, MPISTATUS, MPIERR)
      !write(*,*) "I'm node", mynode, "got energy", myenergy
      !write(*,*) "I'm node", mynode, "recieving coords"
      CALL MPI_RECV(MYCOORDS, 3*NATOMS, MPI_DOUBLE_PRECISION, HEADNODE, &
         MPITAG, MPI_COMM_WORLD, MPISTATUS, MPIERR)
      !write(*,'(A,I5,A,3F20.10)') "I'm node", mynode, " got coords", mycoords(1:3)
   ENDIF

   !write(*,'(A,I5,A,4F20.10)') "I'm node", mynode, " new", myenergy, mycoords(3*natoms-2:3*natoms)

ENDIF


if (exchangeaccept) then
   write(myunit,*) "exchanging"
else
   if (debug) write(myunit,*) "not exchanging"
endif

#else
      RETURN
#endif

END SUBROUTINE TRYEXCHANGE_INDEPENDENCE

