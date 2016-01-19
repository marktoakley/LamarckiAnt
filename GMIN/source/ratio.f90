! RATIO routines
!----------------------------------------------------------------------------------------------------------------------------------!
! NULLMOVE counts the steps that don't move to new minima. Identity is
! identified when the distance between minima is less than GEOMDIFFTOL. In this
! instance, the NULLMOVES counter is incremented.
!
! NULLMOVE also acts as the Metropolis test.
SUBROUTINE NULLMOVE(NEWCOORDS,OLDCOORDS,ATEST,NULLMOVES,JP)
USE COMMONS, ONLY: NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, PERIODIC, GEOMDIFFTOL, RIGID, TWOD, TEMP, MYUNIT, DEBUG, NPAR
IMPLICIT NONE
DOUBLE PRECISION  :: NEWCOORDS(3*NATOMS), OLDCOORDS(3*NATOMS), OLDCOORDSCP(3*NATOMS)
DOUBLE PRECISION  :: RMAT(3,3), DISTANCE, DIST2, ENEW, EOLD, GRAD(3*NATOMS), DPRAND
INTEGER           :: NULLMOVES(NPAR), JP
LOGICAL           :: ATEST 

OLDCOORDSCP(:)=OLDCOORDS(:) ! make local copy to avoid returning modified markov coordinates
CALL MINPERMDIST(NEWCOORDS,OLDCOORDSCP,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
CALL POTENTIAL(OLDCOORDSCP,GRAD,EOLD,.FALSE.,.FALSE.) ! possibly
CALL POTENTIAL(NEWCOORDS,GRAD,ENEW,.FALSE.,.FALSE.) ! redundant

IF (DISTANCE<GEOMDIFFTOL) THEN ! you didn't move to a new minimum
   ATEST=.FALSE.
   NULLMOVES(JP)=NULLMOVES(JP)+1
   IF (DEBUG) WRITE(MYUNIT,*)ATEST,"DISTFAIL, DISTANCE=",DISTANCE
ELSE IF ((EXP(-(ENEW-EOLD)/MAX(TEMP(JP),1.D-100))<DPRAND()).OR.(ENEW<-1.0D6)) THEN ! you failed the Metropolis test
   ATEST=.FALSE.
   IF (DEBUG) WRITE(MYUNIT,*)ATEST,"EFAIL, EDIFF=",ENEW-EOLD
ELSE ! you moved somewhere and passed the metropolis test 
   ! ATEST=.TRUE.
   IF (DEBUG) WRITE(MYUNIT,*)ATEST,"SUCCESS, DISTANCE=",DISTANCE,", DELTAE=",ENEW-EOLD
ENDIF

END SUBROUTINE NULLMOVE

!----------------------------------------------------------------------------------------------------------------------------------!
! FIXRATIO determines whether to increase or decrease the step size(s) or
! temperature. It also prints helpful statistics.
!
! P0 is the probability that moves end in go to new minima. If P0 is less than
! the step ratio, SRATIO, the step size(s) are increased. Otherwise, conversely.
!
! P1 is the probability that a move passes the Metropolis test GIVEN that move
! has ended in a different minimum. If P1 is less than the temperature ratio,
! TRATIO, the temperature is increased. Otherwise, conversely.
SUBROUTINE FIXRATIO(NSTEPS,NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
USE COMMONS, ONLY: STEP,OSTEP,TEMP,NACCEPT,SRATIO,TRATIO,NPAR,MYUNIT,SUMTEMP,SUMSTEP,SUMOSTEP,SUPPRESST,TMOVE,OMOVE, &
                 & GROUPROTT,GR_SCALEPROB,GR_SCALEROT 
USE GROUPROTMOD
IMPLICIT NONE
DOUBLE PRECISION  :: P0, P1, P0T, P1T
INTEGER           :: JP, NULLMOVES(NPAR), NSUCCESS(NPAR), NULLMOVEST(NPAR), NSUCCESST(NPAR), NSTEPS, LOCALSTEPS

LOCALSTEPS=MOD(NSTEPS,NACCEPT) ! these lines are necessary for
IF (LOCALSTEPS==0) LOCALSTEPS=NACCEPT ! the final summary output.

P0=1.D0-1.D0*NULLMOVES(JP)/(1.D0*LOCALSTEPS) ! probability moves go to new minima
P1=1.D0*NSUCCESS(JP)/(1.D0*(LOCALSTEPS-NULLMOVES(JP))) ! probability, given new minimum, a move is accepted

SUMSTEP=SUMSTEP+STEP(JP)*LOCALSTEPS/NACCEPT ! sum of translation step sizes
SUMOSTEP=SUMOSTEP+OSTEP(JP)*LOCALSTEPS/NACCEPT ! sum of rotation step sizes
SUMTEMP=SUMTEMP+TEMP(JP)*LOCALSTEPS/NACCEPT ! sum of temperatures

IF (SRATIO>0) THEN
   IF (P0>SRATIO) THEN ! decrease step size(s)
      IF (TMOVE(JP)) STEP(JP)=STEP(JP)/1.01
      IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/1.01
! GROUPROTATION step scaling
      IF(GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,1.0D0/1.01)
   ELSE ! increase step size(s)
      IF (TMOVE(JP)) STEP(JP)=STEP(JP)*1.01
      IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)*1.01
! GROUPROTATION step scaling
      IF(GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,1.0D0*1.01)
   ENDIF
ENDIF
IF (TRATIO>0) THEN ! decrease temperature
   IF (P1>TRATIO) THEN
      TEMP(JP)=TEMP(JP)/1.1
   ELSE ! increase temperature
      TEMP(JP)=TEMP(JP)*1.1
   ENDIF
ENDIF

NULLMOVEST(JP)=NULLMOVEST(JP)+NULLMOVES(JP) ! total number of moves that stay in the same minimum
NULLMOVES(JP)=0 ! reset local counter
NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP) ! total number of accepted moves
NSUCCESS(JP)=0 ! reset local counter

P0T=1.D0-1.D0*NULLMOVEST(JP)/(1.D0*NSTEPS) ! simulation-wide probability of moving to a new basin
P1T=1.D0*NSUCCESST(JP)/(1.D0*(NSTEPS-NULLMOVEST(JP))) ! simulation wide probablity of accepting a step

IF (.NOT.SUPPRESST) THEN ! print some useful statistics
   WRITE(MYUNIT,'(4(A,F8.5),A)')'FIXRATIO>Step ratio (mean) = ',P0,' (',P0T,'); Temperature ratio (mean) = ',P1,' (',P1T,')'
   WRITE(MYUNIT,'(A,2F10.5,A,F10.5)')'FIXRATIO>Step is now ',STEP(JP),OSTEP(JP),'; Temperature is now ',TEMP(JP)
   IF (JP==NPAR) WRITE(MYUNIT,'(A,2F10.5,A,F10.5)')'FIXRATIO>Average step =',SUMSTEP/NPAR/NSTEPS*NACCEPT,&
                       SUMOSTEP/NPAR/NSTEPS*NACCEPT,'; Average temperature =',SUMTEMP/NPAR/NSTEPS*NACCEPT
ENDIF

END SUBROUTINE FIXRATIO

!----------------------------------------------------------------------------------------------------------------------------------!
! FINALRATIO prints summary statistics at the end of the simulation.
SUBROUTINE FINALRATIO(NULLMOVEST,NSUCCESST,JP,NSTEPS)
USE COMMONS, ONLY: MYUNIT,NPAR,SUMTEMP,SUMSTEP,SUMOSTEP,NACCEPT
IMPLICIT NONE
DOUBLE PRECISION  :: P0, P1
INTEGER           :: JP, NULLMOVEST(NPAR), NSUCCESST(NPAR), NSTEPS

IF (JP<NPAR) RETURN ! what?
 
P0=1.D0-1.D0*SUM(NULLMOVEST(:))/(1.D0*NSTEPS)
P1=1.D0*SUM(NSUCCESST(:))/(1.D0*(NSTEPS-SUM(NULLMOVEST(:))))
SUMSTEP=SUMSTEP/NSTEPS*NACCEPT
SUMOSTEP=SUMOSTEP/NSTEPS*NACCEPT
SUMTEMP=SUMTEMP/NSTEPS*NACCEPT
 
WRITE(MYUNIT,'(3(A,I0))')'FINALRATIO>Total:',NSTEPS,' Null:',SUM(NULLMOVEST(:)),' Success:',SUM(NSUCCESST(:))
WRITE(MYUNIT,'(2(A,F8.5))')'FINALRATIO>Step ratio = ',P0,'; Temperature ratio = ',P1
WRITE(MYUNIT,'(A,2F10.5,A,F10.5)')'FINALRATIO>Average step =',SUMSTEP,SUMOSTEP,&
                                  '; Average temperature =',SUMTEMP

END SUBROUTINE FINALRATIO

!----------------------------------------------------------------------------------------------------------------------------------!
