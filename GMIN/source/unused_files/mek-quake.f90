!
!  Calculate Cv curve from histograms of quench probabilities 
!  for replicas at different temperatures REPT(J1).
!  Minimise chi^2 statistic to extract best probabilities.
!
PROGRAM CALCCV
IMPLICIT NONE
DOUBLE PRECISION TMIN, TMAX, TINT, Z0, Z1, Z2, TEMPERATURE, DIFF, CHI2, RMS, VTOTAL2, BESTRMS
DOUBLE PRECISION DUMMY2, DUMMY3
DOUBLE PRECISION, ALLOCATABLE :: QENERGY(:), IENERGY(:), IVISITS(:,:),WEIGHTIVISITS(:,:), QVISITS(:,:)
DOUBLE PRECISION, ALLOCATABLE :: VAR(:), WIJ(:,:), BESTVAR(:), REPT(:), VREPTOTAL(:), RELWEIGHTS(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: VAR2(:), BESTVAR2(:)
DOUBLE PRECISION, ALLOCATABLE :: GRAD(:), DGRAD(:)
DOUBLE PRECISION CHI2PLUS, CHI2MINUS, TOL, VARETOTAL, MAXBHSTEP, DPRAND, TSTAR, QFAC, Z0MAX, DUMMY
INTEGER J1, NTEMP, J2, KAPPA, NREP, NVAR, NEVAR, NCOUNT, MUPDATES, ITMAX, ITDONE, NHOPS, RNDSEED, J3, &
  &     NDUMMY, BSPTDUMPFRQ, NQBINS, NIBINS, NAPPEND, BSPTSTART, NREPVARS, NEVAR2, NCOUNT2, NLAST
CHARACTER(LEN=80) FSTRING
CHARACTER(LEN=80) IDSTRING
CHARACTER(LEN=8) DUMMYSTRING
LOGICAL, ALLOCATABLE :: ALLZERO(:), ALLZERO2(:)
LOGICAL CONVERGED, YESNO, FIRSTGO

OPEN(UNIT=10,FILE='Cv.data.doubleT',STATUS='OLD')
READ(10,*) NREP, NQBINS, NIBINS, TMIN, TMAX, NTEMP, KAPPA, TSTAR, BSPTDUMPFRQ, BSPTSTART
PRINT '(A,3I8,2G20.10,2I10,G20.10,2(I8,1X))','nrep,# q bins, # i bins,Tmin,Tmax,#T,kappa,T*,inc,start=', &
  &    NREP, NQBINS, NIBINS, TMIN, TMAX, NTEMP, KAPPA, TSTAR, BSPTDUMPFRQ, BSPTSTART

ALLOCATE(REPT(NREP),VREPTOTAL(NREP)) ! temperature of replicas

ALLOCATE(QVISITS(NQBINS,NREP),QENERGY(NQBINS),IENERGY(NIBINS),IVISITS(NIBINS,NREP), &
&  WEIGHTIVISITS(NIBINS,NREP),ALLZERO(NIBINS),RELWEIGHTS(NIBINS,NIBINS,NREP))

IDSTRING=''
NAPPEND=BSPTSTART-BSPTDUMPFRQ

! Loop over Visits.His.<n> intermediate dumps if required
! Start with the final dump, which may be the only one

FIRSTGO=.TRUE.
222 CONTINUE ! loop over Visits.His.<n> intermediate dumps if required

IVISITS(1:NIBINS,1:NREP)=0.0D0
WEIGHTIVISITS(1:NIBINS,1:NREP)=0.0D0
RELWEIGHTS(1:NIBINS,1:NIBINS,1:NREP)=0.0D0
DO J1=1,NREP
   IF (J1-1.GE.10) THEN
      WRITE(FSTRING,'(A11,I2)') 'Visits.his.',J1-1
   ELSE
      WRITE(FSTRING,'(A11,I1)') 'Visits.his.',J1-1
   ENDIF
   FSTRING=TRIM(ADJUSTL(FSTRING)) // TRIM(ADJUSTL(IDSTRING))
   INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',FSTRING,' not found - quit'
      STOP
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD')
   READ(1,*) REPT(J1)
   VREPTOTAL(J1)=0.0D0
   DO J2=1,NQBINS
      READ(1,*) QENERGY(J2),QVISITS(J2,J1)
      VREPTOTAL(J1)=VREPTOTAL(J1)+QVISITS(J2,J1)
   ENDDO
   DO J2=1,NQBINS
      READ(1,*) FSTRING ! This is the dummy "Instantaneous PE histogram for quench bin energy" line
      QFAC=(QENERGY(J2)-QENERGY(1))*(1.0D0/REPT(J1)-1.0D0/TSTAR)
      DO J3=1,NIBINS
         READ(1,*) IENERGY(J3), NDUMMY
         IVISITS(J3,J1)=IVISITS(J3,J1)+NDUMMY
         IF (NDUMMY.GT.0.0D0) THEN ! Can prevent Inf and NaN if exponent overflows but NDUMMY is zero
            WEIGHTIVISITS(J3,J1)=WEIGHTIVISITS(J3,J1)+NDUMMY*EXP((IENERGY(J3)-IENERGY(1))/TSTAR+QFAC)
         ENDIF
      ENDDO
   ENDDO
   CLOSE(1)
   WRITE(*,'(A,I5,2(A,G20.10))') 'temperature for replica ',J1,' is ',REPT(J1),' total visits=',VREPTOTAL(J1)
   PRINT '(A)',' '
   IF (J1-1.GE.10) THEN
      WRITE(FSTRING,'(A12,I2)') 'Visits2.his.',J1-1
   ELSE
      WRITE(FSTRING,'(A12,I1)') 'Visits2.his.',J1-1
   ENDIF
   FSTRING=TRIM(ADJUSTL(FSTRING)) // TRIM(ADJUSTL(IDSTRING))
   INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',FSTRING,' not found - quit'
      STOP
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD')
   READ(1,*) DUMMY ! should be REPT(J1) again
   DUMMY=0.0D0
   DO J2=1,NIBINS
      READ(1,*) DUMMY2, DUMMY3
      IF (DUMMY3.NE.IVISITS(J2,J1)) THEN
         PRINT '(2(A,I8),2(A,G20.10))','ERROR - visits to PE bin ',J2,' in replica ',J1,' should be ',IVISITS(J2,J1), &
  &                                    ' not ',DUMMY3
         STOP
      ENDIF
   ENDDO
   DO J2=1,NIBINS
      READ(1,*) DUMMYSTRING
      DO J3=1,NIBINS
         READ(1,*) DUMMY, DUMMY2, DUMMY3
         IF (DUMMY3.GT.0.0D0) RELWEIGHTS(J3,J2,J1)=DUMMY2/DUMMY3 ! J3 labels the inferred weight bin
                                                                 ! J2 labels the instantaneous PE bin
                                                                 ! J1 labels the replica
      ENDDO
   ENDDO
   CLOSE(1)
ENDDO

NEVAR=0
VTOTAL2=0.0D0
eloop1: DO J1=1,NIBINS
   ALLZERO(J1)=.TRUE.
   DO J2=1,NREP
      IF (IVISITS(J1,J2).GT.0.0D0) THEN
         ALLZERO(J1)=.FALSE.
         NEVAR=NEVAR+1
         VTOTAL2=VTOTAL2+IVISITS(J1,J2)
         CYCLE eloop1
      ENDIF
   ENDDO
ENDDO eloop1
NVAR=NREP+NEVAR
PRINT '(3(A,I5))','number of V^I bin weight variables=',NEVAR,' number of Z variables=',NREP,' total=',NVAR
!
! Calculate the fixed quantities 
!
IF (ALLOCATED(WIJ)) DEALLOCATE(WIJ)
ALLOCATE(WIJ(NEVAR,NREP))
NCOUNT=0
eloop2: DO J1=1,NIBINS
   IF (ALLZERO(J1)) CYCLE eloop2
   NCOUNT=NCOUNT+1
   DO J2=1,NREP
      IF (WEIGHTIVISITS(J1,J2).GT.0.0D0) THEN
         WIJ(NCOUNT,J2)=LOG(WEIGHTIVISITS(J1,J2))
      ENDIF
   ENDDO
ENDDO eloop2

IF (ALLOCATED(VAR)) DEALLOCATE(VAR,BESTVAR)
ALLOCATE(VAR(NEVAR+NREP),BESTVAR(NEVAR+NREP))
!
! Initialise variables. First NEVAR of VAR are the quench energy probabilities and 
! the last NREP are the partition functions.
! VAR(1) is frozen as 0.0D0.
!
DO J1=1,NEVAR+NREP
   VAR(J1)=J1-1
ENDDO
!
! Test analytic derivatives
!
IF (ALLOCATED(GRAD)) DEALLOCATE(GRAD,DGRAD)
ALLOCATE(GRAD(NVAR),DGRAD(NVAR))
IF (.FALSE.) THEN
   CALL GETCHI2(NIBINS,NEVAR,NREP,WIJ,VAR,CHI2,GRAD,ALLZERO,IVISITS,RMS,VTOTAL2)
   DIFF=0.001
   DO J1=1,NVAR
      VAR(J1)=VAR(J1)+DIFF
      CALL GETCHI2(NIBINS,NEVAR,NREP,WIJ,VAR,CHI2PLUS,DGRAD,ALLZERO,IVISITS,RMS,VTOTAL2)
      VAR(J1)=VAR(J1)-2.0D0*DIFF
      CALL GETCHI2(NIBINS,NEVAR,NREP,WIJ,VAR,CHI2MINUS,DGRAD,ALLZERO,IVISITS,RMS,VTOTAL2)
      VAR(J1)=VAR(J1)+DIFF
      IF (GRAD(J1).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
      ENDIF
   ENDDO
ENDIF

MUPDATES=100
TOL=1.0D-7
ITMAX=1000000
NHOPS=1
BESTRMS=1.0D100
RNDSEED=1
MAXBHSTEP=1.0D0
CALL SDPRND(RNDSEED)

DO J1=1,NHOPS

   CALL MYLBFGS(NVAR,MUPDATES,VAR,TOL,ITDONE,ITMAX,CHI2,CONVERGED,NIBINS,NEVAR,NREP,WIJ,ALLZERO,IVISITS,VTOTAL2,RMS)

   IF (RMS.LT.BESTRMS) THEN
      BESTRMS=RMS
      BESTVAR(1:NVAR)=VAR(1:NVAR)
   ENDIF

   DO J2=1,NVAR
      VAR(J2)=VAR(J2)+VAR(J2)*(2.0D0*DPRAND()-1.0D0)*MAXBHSTEP
   ENDDO

ENDDO

VAR(1:NVAR)=EXP(BESTVAR(1:NVAR))
PRINT '(A)','Final bin weights and relative values:'
NCOUNT=0
VARETOTAL=0.0D0
DO J1=1,NIBINS
   IF (ALLZERO(J1)) THEN
      PRINT '(3E20.10)',IENERGY(J1),0.0D0,0.0D0
   ELSE 
      NCOUNT=NCOUNT+1
      PRINT '(3E20.10)',IENERGY(J1),VAR(NCOUNT),VAR(NCOUNT)/VAR(1)
      VARETOTAL=VARETOTAL+VAR(NCOUNT)
   ENDIF
ENDDO
!DO J1=1,NREP
!   ZCALC(J1)=0.0D0
!   NCOUNT=0
!   DO J2=1,NIBINS
!      IF (.NOT.ALLZERO(J2)) THEN
!         NCOUNT=NCOUNT+1
!         ZCALC(J1)=ZCALC(J1)+VAR(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/REPT(J1))
!      ENDIF
!   ENDDO
!ENDDO
!DO J1=1,NREP
!   PRINT '(A,I5,A,G20.10)','Actual and calculated non-zero weights for replica ',J1,' T=',REPT(J1)
!   NCOUNT=0
!   VNORM=0.0D0
!   DO J2=1,NIBINS
!      IF (.NOT.ALLZERO(J2)) THEN
!         NCOUNT=NCOUNT+1
!         VNORM=VNORM+IVISITS(J2,J1)*EXP((IENERGY(J2)-IENERGY(1))/REPT(J1))
!      ENDIF
!   ENDDO
!   NCOUNT=0
!   DO J2=1,NIBINS
!      IF (.NOT.ALLZERO(J2)) THEN
!         NCOUNT=NCOUNT+1
!         VTOTAL(J1)=VTOTAL(J1)+IVISITS(J2,J1)
!         WRITE(*,'(4G20.10)') IVISITS(J2,J1)*EXP((IENERGY(J2)-IENERGY(1))/REPT(J1))/VNORM,VAR(NCOUNT)/VARETOTAL, &
!  &                           IVISITS(J2,J1)*EXP((IENERGY(J2)-IENERGY(1))/REPT(J1))*VARETOTAL/(VNORM*VAR(NCOUNT))
!      ENDIF
!   ENDDO
!ENDDO

!PRINT '(A)','Final partition functions, relative values, recalculated values, ratio:'
!DO J1=1,NREP
!   PRINT '(5E20.10)',REPT(J1),VAR(NEVAR+J1),VAR(NEVAR+J1)/VAR(NEVAR+1),ZCALC(J1),VAR(NEVAR+J1)/ZCALC(J1)
!ENDDO

!PRINT '(A)','Final partition functions and relative values:'
!DO J1=1,NREP
!   PRINT '(5E20.10)',REPT(J1),VAR(NEVAR+J1),VAR(NEVAR+J1)/VAR(NEVAR+1)
!ENDDO

TINT=(TMAX-TMIN)/(NTEMP-1)
!
!  Calculate Z0, Z1 and Z2 over the required T range. Omit factors of (kT/h)^kappa,
!  which cancel.
!
FSTRING='Cv.out' // TRIM(ADJUSTL(IDSTRING))
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
! PRINT '(A)','Most probable energy as a function of T:'
DO J1=1,NTEMP
   Z0=0.0D0
   Z1=0.0D0
   Z2=0.0D0
   TEMPERATURE=TMIN+(J1-1)*TINT
   NCOUNT=0
   Z0MAX=-1.0D0
   DO J2=1,NIBINS
      IF (ALLZERO(J2)) CYCLE 
      NCOUNT=NCOUNT+1
      DUMMY=VAR(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)
      Z0=Z0+DUMMY
      Z1=Z1+VAR(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)*(IENERGY(J2)-IENERGY(1))
      Z2=Z2+VAR(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)*(IENERGY(J2)-IENERGY(1))**2
      IF (DUMMY.GT.Z0MAX) THEN
         Z0MAX=DUMMY
         NDUMMY=J2
      ENDIF
   ENDDO
!  PRINT '(A,F15.4,A,I7,2(A,G20.10))','For T=',TEMPERATURE,' bin ',NDUMMY,' energy=',IENERGY(NDUMMY),' Z=',Z0MAX
   WRITE(1,'(6G20.10)') TEMPERATURE, Z0, Z1, Z2, KAPPA*TEMPERATURE/2.0D0 + Z1/Z0, &
  &                     KAPPA/2.0D0-(Z1/(Z0*TEMPERATURE))**2 + Z2/(Z0*TEMPERATURE**2)
ENDDO
CLOSE(1)

! GOTO 444 to skip the quench DoS phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Second pass using quench DoS inferrence !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The number of omega bin weight variables can now increase to include any instantaneous pe bin for
! which there is a non-zero relative weight.
! The number of normalisation variables is equal to the number of non-zero visits to all the 
! instantaneous pe bins summed over replicas.
!
IF (ALLOCATED(ALLZERO2)) DEALLOCATE(ALLZERO2)
ALLOCATE(ALLZERO2(NIBINS))
ALLZERO2(1:NIBINS)=.TRUE.
VTOTAL2=0.0D0
DO J1=1,NREP
   jbin2 : DO J2=1,NIBINS
      IF (IVISITS(J2,J1).EQ.0.0D0) CYCLE jbin2
      VTOTAL2=VTOTAL2+IVISITS(J2,J1)
      DO J3=1,J2
         IF (RELWEIGHTS(J3,J2,J1).GT.0.0D0) THEN
            ALLZERO2(J3)=.FALSE.
         ENDIF
      ENDDO
   ENDDO jbin2
ENDDO

NEVAR2=0
DO J1=1,NIBINS
   IF (.NOT.ALLZERO2(J1)) NEVAR2=NEVAR2+1
ENDDO

NREPVARS=0
DO J1=1,NREP
   DO J2=1,NIBINS
      IF (IVISITS(J2,J1).GT.0.0D0) NREPVARS=NREPVARS+1
   ENDDO
ENDDO
NVAR=NEVAR2+NREPVARS
PRINT '(3(A,I5))','number of V^I bin weight variables=',NEVAR2,' number of normalisation variables=',NREPVARS,' total=',NVAR

IF (ALLOCATED(VAR2)) DEALLOCATE(VAR2,BESTVAR2)
ALLOCATE(VAR2(NVAR),BESTVAR2(NVAR))
!
! Initialise variables on first round only. First NEVAR2 of VAR are the omega values and
! the last NREPVARS are the multiplicative factors, one for each non-zero visits instantaneous
! PE bin in each replica.
! VAR(1) is frozen as 0.0D0. Other values are initialised according to the values obtained
! without the quench visits data.
!
! IF (.NOT.FIRSTGO) GOTO 111
NCOUNT=0
NCOUNT2=0
NLAST=0
VAR2(1:NVAR)=0.0D0
DO J2=1,NIBINS
   IF (.NOT.ALLZERO2(J2)) NCOUNT2=NCOUNT2+1
   IF (.NOT.ALLZERO(J2)) THEN
      NCOUNT=NCOUNT+1 ! this is a variable omega in the first scheme
      IF (VAR(NCOUNT).GE.0.0D0) VAR2(NLAST+1:NCOUNT2)=LOG(VAR(NCOUNT))
!     DO J3=NLAST+1,NCOUNT2
!        PRINT '(A,I6,A,G20.10,A,I6)','Initialising var2 element ',J3,' to ',VAR(NCOUNT),' NCOUNT=',NCOUNT
!     ENDDO
      NLAST=NCOUNT2
   ENDIF
ENDDO 

NCOUNT=0
DO J1=1,NREP
   DO J2=1,NIBINS
      IF (IVISITS(J2,J1).GT.0.0D0) THEN
         NCOUNT=NCOUNT+1
         IF (VAR(NEVAR+J1).GE.0.0D0) VAR2(NEVAR2+NCOUNT)=LOG(VAR(NEVAR+J1))
!        PRINT '(A,I6,A,G20.10,A,I6)','Initialising var2 element ',NEVAR2+NCOUNT,' to ',VAR(NEVAR+J1),' NCOUNT=',NCOUNT
      ENDIF
   ENDDO
ENDDO

111 CONTINUE
!
! Test analytic derivatives
!
IF (.FALSE.) THEN
   IF (ALLOCATED(GRAD)) DEALLOCATE(GRAD,DGRAD)
   ALLOCATE(GRAD(NVAR),DGRAD(NVAR))
   DO J1=1,NVAR
      VAR2(J1)=DPRAND()
   ENDDO
   CALL GETCHIX2(NIBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,VAR2,CHI2,GRAD,ALLZERO2,IVISITS,RMS,VTOTAL2,RELWEIGHTS)
   DIFF=0.001
   DO J1=1,NREPVARS
      VAR2(J1)=VAR2(J1)+DIFF
      CALL GETCHIX2(NIBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,VAR2,CHI2PLUS,DGRAD,ALLZERO2,IVISITS,RMS,VTOTAL2,RELWEIGHTS)
      VAR2(J1)=VAR2(J1)-2.0D0*DIFF
      CALL GETCHIX2(NIBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,VAR2,CHI2MINUS,DGRAD,ALLZERO2,IVISITS,RMS,VTOTAL2,RELWEIGHTS)
      VAR2(J1)=VAR2(J1)+DIFF
      IF (GRAD(J1).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
      ELSE
         PRINT '(A,I5,2G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1)
      ENDIF
   ENDDO
ENDIF

MUPDATES=100
TOL=1.0D-5
ITMAX=1000000
CALL MYLBFGS2(NVAR,MUPDATES,VAR2,TOL,ITDONE,ITMAX,CHI2,CONVERGED,NIBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,ALLZERO2,IVISITS, &
  &           VTOTAL2,RMS,RELWEIGHTS)

PRINT '(A)','Final bin weights and relative values:'
NCOUNT=0
VAR2(1:NVAR)=EXP(VAR2(1:NVAR))
DO J1=1,NIBINS
   IF (ALLZERO2(J1)) THEN
      PRINT '(3E20.10)',IENERGY(J1),0.0D0,0.0D0
   ELSE 
      NCOUNT=NCOUNT+1
      PRINT '(3E20.10)',IENERGY(J1),VAR2(NCOUNT),VAR2(NCOUNT)/VAR2(1)
   ENDIF
ENDDO

TINT=(TMAX-TMIN)/(NTEMP-1)
!
!  Calculate Z0, Z1 and Z2 over the required T range. Omit factors of (kT/h)^kappa,
!  which cancel.
!
FSTRING='Cv.out.qd' // TRIM(ADJUSTL(IDSTRING))
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
PRINT '(A)','Most probable energy as a function of T:'
DO J1=1,NTEMP
   Z0=0.0D0
   Z1=0.0D0
   Z2=0.0D0
   TEMPERATURE=TMIN+(J1-1)*TINT
   NCOUNT=0
   Z0MAX=-1.0D0
   DO J2=1,NIBINS
      IF (ALLZERO2(J2)) CYCLE 
      NCOUNT=NCOUNT+1
      DUMMY=VAR2(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)
      Z0=Z0+DUMMY
      Z1=Z1+VAR2(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)*(IENERGY(J2)-IENERGY(1))
      Z2=Z2+VAR2(NCOUNT)*EXP(-(IENERGY(J2)-IENERGY(1))/TEMPERATURE)*(IENERGY(J2)-IENERGY(1))**2
      IF (DUMMY.GT.Z0MAX) THEN
         Z0MAX=DUMMY
         NDUMMY=J2
      ENDIF
   ENDDO
   PRINT '(A,F15.4,A,I7,2(A,G20.10))','For T=',TEMPERATURE,' bin ',NDUMMY,' energy=',IENERGY(NDUMMY),' Z=',Z0MAX
   WRITE(1,'(6G20.10)') TEMPERATURE, Z0, Z1, Z2, KAPPA*TEMPERATURE/2.0D0 + Z1/Z0, &
  &                     KAPPA/2.0D0-(Z1/(Z0*TEMPERATURE))**2 + Z2/(Z0*TEMPERATURE**2)
ENDDO
CLOSE(1)

444 CONTINUE

!
! Are there intermediate Visits.His.<n> to do?
!
NAPPEND=NAPPEND+BSPTDUMPFRQ
IF ((NAPPEND.LT.100000).OR.(NAPPEND.GE.100000000)) THEN
   PRINT '(A)','NAPPEND is out of allowed range - quit'
   GOTO 333
ENDIF
IF (NAPPEND.GE.10000000) THEN
   WRITE(IDSTRING,'(A1,I8)') '.', NAPPEND
ELSEIF (NAPPEND.GE.1000000) THEN
   WRITE(IDSTRING,'(A1,I7)') '.', NAPPEND
ELSEIF (NAPPEND.GE.100000) THEN
   WRITE(IDSTRING,'(A1,I6)') '.', NAPPEND
ELSEIF (NAPPEND.GE.10000) THEN
   WRITE(IDSTRING,'(A1,I5)') '.', NAPPEND
ENDIF
FIRSTGO=.FALSE.
GOTO 222 ! do the Visits histograms corresponding to NAPPEND steps

333 CONTINUE

END PROGRAM CALCCV

SUBROUTINE GETCHI2(NBINS,NEVAR,NREP,WIJ,VAR,CHI2,GRAD,ALLZERO,VISITS,RMS,VTOTAL2)
IMPLICIT NONE
INTEGER NEVAR, NREP, NBINS, NCOUNT, J1, J2
DOUBLE PRECISION GRAD(NEVAR+NREP), WIJ(NEVAR,NREP), VAR(NEVAR+NREP), CHI2, RIJ, VISITS(NBINS,NREP), RMS, VTOTAL2
LOGICAL ALLZERO(NBINS)

CHI2=0.0D0
GRAD(1:NEVAR+NREP)=0.0D0
NCOUNT=0
eloop: DO J1=1,NBINS
   IF (ALLZERO(J1)) CYCLE eloop
   NCOUNT=NCOUNT+1
   DO J2=1,NREP
      IF (VISITS(J1,J2).GT.0.0D0) THEN
         RIJ=WIJ(NCOUNT,J2)-(VAR(NCOUNT)-VAR(NEVAR+J2))
         CHI2=CHI2+VISITS(J1,J2)*RIJ**2
         GRAD(NCOUNT)=GRAD(NCOUNT)    -VISITS(J1,J2)*RIJ
         GRAD(NEVAR+J2)=GRAD(NEVAR+J2)+VISITS(J1,J2)*RIJ
      ENDIF
   ENDDO
ENDDO eloop
CHI2=CHI2/VTOTAL2
GRAD(1:NEVAR+NREP)=2*GRAD(1:NEVAR+NREP)/VTOTAL2
GRAD(1)=0.0D0 ! freeze first weight
RMS=0.0D0
DO J1=1,NEVAR+NREP
   RMS=RMS+GRAD(J1)**2
ENDDO
RMS=SQRT(RMS/(NEVAR+NREP))

END SUBROUTINE GETCHI2

SUBROUTINE GETCHIX2(NBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,VAR2,CHI2,GRAD,ALLZERO2,VISITS,RMS,VTOTAL2,RELWEIGHTS)
IMPLICIT NONE
INTEGER NEVAR2, NREP, NBINS, NCOUNT, J1, J2, NREPVARS, J3, NCOUNT2, NCOUNT3, NEVAR
DOUBLE PRECISION GRAD(NEVAR2+NREPVARS), WIJ(NEVAR,NREP), VAR2(NEVAR2+NREPVARS), CHI2, RIJ, VISITS(NBINS,NREP), RMS, VTOTAL2
DOUBLE PRECISION RELWEIGHTS(NBINS,NBINS,NREP)
LOGICAL ALLZERO2(NBINS)

CHI2=0.0D0
GRAD(1:NEVAR2+NREPVARS)=0.0D0

NCOUNT2=0
DO J1=1,NREP
   NCOUNT3=0
   jbin : DO J2=1,NBINS
      IF (VISITS(J2,J1).GT.0.0D0) THEN
         NCOUNT2=NCOUNT2+1 ! counts bins with non-zero visits over all replicas
         NCOUNT3=NCOUNT3+1 ! counts bins with non-zero visits for replica J1
      ELSE
         CYCLE jbin
      ENDIF
      NCOUNT=0 
      ibin : DO J3=1,J2 ! weights must vanish above J2
         IF (ALLZERO2(J3)) CYCLE ibin 
         NCOUNT=NCOUNT+1 ! counts omega values that aren't identically zero, i.e. variables
         IF (RELWEIGHTS(J3,J2,J1).GT.0.0D0) THEN
            RIJ=WIJ(NCOUNT3,J1)+LOG(RELWEIGHTS(J3,J2,J1))-(VAR2(NCOUNT)-VAR2(NEVAR2+NCOUNT2))
            CHI2=CHI2+VISITS(J2,J1)*RIJ**2
            GRAD(NCOUNT)=GRAD(NCOUNT)    -VISITS(J2,J1)*RIJ
            GRAD(NEVAR2+NCOUNT2)=GRAD(NEVAR2+NCOUNT2)+VISITS(J2,J1)*RIJ
         ENDIF
      ENDDO ibin
   ENDDO jbin
ENDDO 

CHI2=CHI2/VTOTAL2
GRAD(1:NEVAR+NREPVARS)=2*GRAD(1:NEVAR+NREPVARS)/VTOTAL2
GRAD(1)=0.0D0 ! freeze first weight
RMS=0.0D0
DO J1=1,NEVAR+NREPVARS
   RMS=RMS+GRAD(J1)**2
ENDDO
RMS=SQRT(RMS/(NEVAR+NREPVARS))

END SUBROUTINE GETCHIX2

SUBROUTINE MYLBFGS(N,M,X,EPS,ITDONE,ITMAX,ENERGY,CONVERGED,NBINS,NEVAR,NREP,WIJ,ALLZERO,VISITS,VTOTAL2,RMS)
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL
DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT
DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),RMS
LOGICAL CONVERGED
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,OVERLAP,DOT1,DOT2,DGUESS,MAXBFGS
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
LOGICAL MFLAG, FAILED
INTEGER NBINS, NEVAR, NREP
DOUBLE PRECISION WIJ(NEVAR,NREP), VISITS(NBINS,NREP), VTOTAL2
LOGICAL ALLZERO(NBINS)

DGUESS=1.0D0
MAXBFGS=100.0D0
ITER=0
ITDONE=0
NFAIL=0
FAILED=.FALSE.
CALL GETCHI2(NBINS,NEVAR,NREP,WIJ,X,ENERGY,G,ALLZERO,VISITS,RMS,VTOTAL2)

! WRITE(*,'(A,2G20.10,A,I6,A)') 'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
   IF (MFLAG) THEN
      WRITE(*,'(A,I8,A,F15.7,A,F15.7)') 'mylbfgs> chi^2 converged in ',ITDONE, &
  &                                                ' steps. Value=',ENERGY,' RMS force=',RMS
      IF (ITER.GT.0) WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
      CONVERGED=.TRUE.
      RETURN
   ENDIF
ENDIF

IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
   WRITE(*,'(A,G15.7,A,G15.7)') 'mylbfgs> **WARNING - chi^2 did not converge, value=',ENERGY,' RMS force=',RMS
   WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
   CONVERGED=.FALSE.
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240   FORMAT('xmylbfgs> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   DO I=1,N
      DIAG(I)=DGUESS
   ENDDO
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   IF (GNORM.EQ.0.0D0) THEN
      GNORM=1.0D0 ! exact zero is presumably wrong!
      PRINT '(A)','WARNING - GNORM was zero in xmylbfgs, resetting to one'
   ENDIF
   STP=MIN(GNORM,1.0D0/GNORM)
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
   IF (YY.EQ.0.0D0) YY=1.0D0
   DO I=1,N
      DIAG(I)= YS/YY
!     DIAG(I)= ABS(YS/YY)
   ENDDO
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
   DO I=1,N
      W(I)= -G(I)
   ENDDO
   CP= POINT
   DO I= 1,BOUND
      CP=CP-1
      IF (CP.EQ. -1)CP=M-1
      SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
      INMC=N+M+CP+1
      IYCN=IYPT+CP*N
      W(INMC)= W(N+CP+1)*SQ
      CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
   ENDDO
  
   DO I=1,N
      W(I)=DIAG(I)*W(I)
   ENDDO

   DO I=1,BOUND
      YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
      BETA= W(N+CP+1)*YR
      INMC=N+M+CP+1
      BETA= W(INMC)-BETA
      ISCN=ISPT+CP*N
      CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
      CP=CP+1
      IF (CP.EQ.M) CP=0
   ENDDO
   STP=1.0D0
ENDIF
!
!  Store the new search direction
!
W(1)=0.0D0 ! freeze first weight
DO I=1,N
   W(ISPT+POINT*N+I)= W(I)
ENDDO

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

IF (OVERLAP.GT.0.0D0) THEN
   PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
ENDIF
      
DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
ENDDO 
NDECREASE=0
20    CONTINUE
CALL GETCHI2(NBINS,NEVAR,NREP,WIJ,X,ENEW,GNEW,ALLZERO,VISITS,RMS,VTOTAL2)

IF (ENEW.EQ.0.0D0) ENEW=1.0D-100 ! to prevent divide by zero
IF (((ENEW-ENERGY)/ABS(ENEW).LE.1.0D-10)) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   DO J1=1,N
      G(J1)=GNEW(J1)
   ENDDO
   WRITE(*,'(A,2G20.10,A,I6,A,F13.10)') &
  &                'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
ELSE
!
!  chi^2 increased - try again with a smaller step size?
!
   IF (NDECREASE.GT.5) THEN  
      NFAIL=NFAIL+1
      PRINT*,' LBFGS step cannot find a lower non-zero eigenvalue, NFAIL=',NFAIL
      ITER=0  !  try resetting
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
      ENDDO 
!     IF (NFAIL.GT.1) FAILED=.TRUE.
      GOTO 10
   ENDIF
   DO J1=1,N
      X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   NDECREASE=NDECREASE+1
   STP=STP/10.0D0
!  WRITE(*,'(A,G20.10,A,G20.10,A,F15.8)') &
! &                         'mylbfgs> chi^2 increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
   GOTO 20
ENDIF
!
!  Compute the new step and gradient change
!
NPT=POINT*N
DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO

POINT=POINT+1
IF (POINT.EQ.M) POINT=0
GOTO 10

RETURN

END SUBROUTINE MYLBFGS

SUBROUTINE MYLBFGS2(N,M,X,EPS,ITDONE,ITMAX,ENERGY,CONVERGED,NBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,ALLZERO2,VISITS, &
  &                 VTOTAL2,RMS,RELWEIGHTS)
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL
DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT
DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),RMS
LOGICAL CONVERGED
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,OVERLAP,DOT1,DOT2,DGUESS,MAXBFGS
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
LOGICAL MFLAG, FAILED
INTEGER NBINS, NEVAR, NREP, NEVAR2, NREPVARS
DOUBLE PRECISION WIJ(NEVAR,NREP), VISITS(NBINS,NREP), VTOTAL2, RELWEIGHTS(NBINS,NBINS,NREP)
LOGICAL ALLZERO2(NBINS)

DGUESS=1.0D0
MAXBFGS=100.0D0
ITER=0
ITDONE=0
NFAIL=0
FAILED=.FALSE.
CALL GETCHIX2(NBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,X,ENERGY,G,ALLZERO2,VISITS,RMS,VTOTAL2,RELWEIGHTS)

WRITE(*,'(A,2G20.10,A,I6,A)') 'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
   IF (MFLAG) THEN
      WRITE(*,'(A,I8,A,F15.7,A,F15.7)') 'mylbfgs> chi^2 converged in ',ITDONE, &
  &                                                ' steps. Value=',ENERGY,' RMS force=',RMS
      IF (ITER.GT.0) WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
      CONVERGED=.TRUE.
      RETURN
   ENDIF
ENDIF

IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
   WRITE(*,'(A,G15.7,A,G15.7)') 'mylbfgs> **WARNING - chi^2 did not converge, value=',ENERGY,' RMS force=',RMS
   WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
   CONVERGED=.FALSE.
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240   FORMAT('xmylbfgs> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   DO I=1,N
      DIAG(I)=DGUESS
   ENDDO
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   IF (GNORM.EQ.0.0D0) THEN
      GNORM=1.0D0 ! exact zero is presumably wrong!
      PRINT '(A)','WARNING - GNORM was zero in xmylbfgs, resetting to one'
   ENDIF
   STP=MIN(GNORM,1.0D0/GNORM)
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
   IF (YY.EQ.0.0D0) YY=1.0D0
   DO I=1,N
      DIAG(I)= YS/YY
!     DIAG(I)= ABS(YS/YY)
   ENDDO
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
   DO I=1,N
      W(I)= -G(I)
   ENDDO
   CP= POINT
   DO I= 1,BOUND
      CP=CP-1
      IF (CP.EQ. -1)CP=M-1
      SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
      INMC=N+M+CP+1
      IYCN=IYPT+CP*N
      W(INMC)= W(N+CP+1)*SQ
      CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
   ENDDO
  
   DO I=1,N
      W(I)=DIAG(I)*W(I)
   ENDDO

   DO I=1,BOUND
      YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
      BETA= W(N+CP+1)*YR
      INMC=N+M+CP+1
      BETA= W(INMC)-BETA
      ISCN=ISPT+CP*N
      CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
      CP=CP+1
      IF (CP.EQ.M) CP=0
   ENDDO
   STP=1.0D0
ENDIF
!
!  Store the new search direction
!
W(1)=0.0D0 ! freeze first weight
DO I=1,N
   W(ISPT+POINT*N+I)= W(I)
ENDDO

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

IF (OVERLAP.GT.0.0D0) THEN
   PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
ENDIF
      
DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
ENDDO 
NDECREASE=0
20    CONTINUE
CALL GETCHIX2(NBINS,NEVAR,NEVAR2,NREP,NREPVARS,WIJ,X,ENEW,GNEW,ALLZERO2,VISITS,RMS,VTOTAL2,RELWEIGHTS)

IF (ENEW.EQ.0.0D0) ENEW=1.0D-100 ! to prevent divide by zero
IF (((ENEW-ENERGY)/ABS(ENEW).LE.1.0D-10)) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   DO J1=1,N
      G(J1)=GNEW(J1)
   ENDDO
   WRITE(*,'(A,2G20.10,A,I6,A,F13.10)') &
  &                'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
ELSE
!
!  chi^2 increased - try again with a smaller step size?
!
   IF (NDECREASE.GT.5) THEN  
      NFAIL=NFAIL+1
      PRINT*,' LBFGS step cannot find a lower non-zero eigenvalue, NFAIL=',NFAIL
      ITER=0  !  try resetting
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
      ENDDO 
!     IF (NFAIL.GT.1) FAILED=.TRUE.
      GOTO 10
   ENDIF
   DO J1=1,N
      X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   NDECREASE=NDECREASE+1
   STP=STP/10.0D0
   WRITE(*,'(A,G20.10,A,G20.10,A,F15.8)') &
  &                         'mylbfgs> chi^2 increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
   GOTO 20
ENDIF
!
!  Compute the new step and gradient change
!
NPT=POINT*N
DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO

POINT=POINT+1
IF (POINT.EQ.M) POINT=0
GOTO 10

RETURN

END SUBROUTINE MYLBFGS2

      double precision function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

!   Copyright (C) 1992  N.M. Maclaren
!   Copyright (C) 1992  The University of Cambridge

!   This software may be reproduced and used freely, provided that all
!   users of it agree that the copyright holders are not liable for any
!   damage or injury caused by use of this software and that this
!   condition is passed onto all subsequent recipients of the software,
!   whether modified or not.



        SUBROUTINE SDPRND (ISEED)
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
!
!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.
!
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
!
!   INDEX must be initialised to an integer between 1 and 101
!   inclusive, POLY(1...N) to integers between 0 and 1000009710
!   inclusive (not all 0), and OTHER to a non-negative proper fraction
!   with denominator 33554432.  It uses the Wichmann-Hill generator to
!   do this.
!
        IX = MOD(ABS(ISEED),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER = AINT(YMOD*X)/YMOD
        OFFSET = 1.0D0/YMOD
        INDEX = 1
        END

        DOUBLE PRECISION FUNCTION DPRAND()
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101), OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0, XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0, TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
!
!   This returns a uniform (0,1) random number, with extremely good
!   uniformity properties.  It assumes that double precision provides
!   at least 33 bits of accuracy, and uses a power of two base.
!
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
!
!   See [Knuth] for why this implements the algorithm described in
!   the paper.  Note that this code is tuned for machines with fast
!   double precision, but slow multiply and divide; many, many other
!   options are possible.
!
        N = INDEX-64
        IF (N .LE. 0) N = N+101
        X = POLY(INDEX)+POLY(INDEX)
        X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY(INDEX) = X
        INDEX = INDEX+1
        IF (INDEX .GT. 101) INDEX = INDEX-101
!
!   Add in the second generator modulo 1, and force to be non-zero.
!   The restricted ranges largely cancel themselves out.
!
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        END
