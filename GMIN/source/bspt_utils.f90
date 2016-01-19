subroutine bspt_takestep_amber(X, Y, Z, IMCSTEP, &
      jumpt, peint, peqv, vnew, exab_count,BETA,deltalnJ )
! choose which takestep routine to use and take the step
! note: MYNODE is in commons.  is it safe to use here?
   use commons
   use grouprotmod
   use modamber9, only : mdstept
   !use grouprotation_ab
   implicit none

   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP
   logical, intent(IN) :: jumpt
   DOUBLE PRECISION, INTENT(IN) :: peint
   DOUBLE PRECISION, INTENT(IN) :: peqv(NENRPER)

   DOUBLE PRECISION, INTENT(INOUT) :: X(NATOMS), Y(NATOMS), Z(NATOMS)
   DOUBLE PRECISION, INTENT(INOUT) :: VNEW
   DOUBLE PRECISION, INTENT(OUT) :: deltalnJ

   integer, intent(INOUT) :: exab_count

   integer j1, j2, j3, j4, k, nchosen
   double precision dummy, dprand,random,Xo(NATOMS),Yo(NATOMS),Zo(NATOMS),BETA(0:NPAR-1)

   !ab2111> debug
   double precision GRAD(3*NATOMS),potel,potelold,point(3*NATOMS),Cpert(3*NATOMS)

   LOGICAL RANDOM_PERT

   CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
   CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   !WRITE(MYUNIT,'(A,F20.10)'), '-1 bspt_utils_amber> before random pert ',POTEL

   IF (AMBERT.EQV..FALSE.) THEN
      WRITE(*,'(A)') "Amber not turned on for this takestep routine!"
      STOP
   ENDIF

   RANDOM_PERT=.TRUE.
   ! random perturbations
   IF (RANDOM_PERT) THEN
      CALL RANDOM_STEP(X,Y,Z)
   ENDIF

   !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )

   CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
   CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   POTELold = POTEL

   !WRITE(MYUNIT,'(A,F20.10)'), '0 bspt_utils_amber> after random pert ',POTEL

   ! grouprotation moves
   IF (GROUPROTT.AND.MOD(INT(IMCSTEP),GROUPROTFREQ).EQ.0) THEN
      CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
      !CALL GROUPROTSTEPab(COORDS(:,MYNODE+1))
      !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )
      CALL GROUPROTSTEP(MYNODE+1)
      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   !   WRITE(MYUNIT,'(A,F20.10)'), '0 bspt_utils_amber> after grouprotation ',POTEL
      !CALL CtoX(COORDS(:,MYNODE+1),X,Y,Z)
      !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )
   ENDIF

   IF (DIHEDRALROTT.AND.MOD(INT(IMCSTEP),DIHEDRALROTFREQ).EQ.0) THEN
      !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )
      !WRITE(MYUNIT,'(A,F20.10,A,F20.10)'), 'a bspt_utils_amber> CA ',COORDS(3*(317-1)+1,MYNODE+1),"C ",COORDS(3*(325-1)+1,MYNODE+1)
      CALL DIHEDRALROTSTEP(deltalnJ,MYNODE+1,BETA,IMCSTEP)
      !WRITE(MYUNIT,'(A,F20.10,A,F20.10)'), 'b bspt_utils_amber> CA ',COORDS(3*(317-1)+1,MYNODE+1),"C ",COORDS(3*(325-1)+1,MYNODE+1)
      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   !   WRITE(MYUNIT,'(A,2F20.10)'), '1 bspt_utils_amber> after dihedral rotation ',POTEL,deltalnJ
      !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )
   ENDIF

   !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )

   !CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
   !CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   !WRITE(MYUNIT,'(A,F20.10)'), '1 bspt_utils_amber> after rotation ',POTEL
   ! ambermd steps
   IF (MDSTEPT) THEN
      !CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
      !WRITE(MYUNIT,'(A,F)'), 'bspt_utils_amber> Coords', X(1)
      CALL TAKESTEPAMBER(MYNODE+1,COORDS(:,MYNODE+1),MOVABLEATOMLIST,NMOVABLEATOMS,LIGMOVET,MDSTEPT,RANDOMSEEDT, &
     &      BLOCKMOVET,NBLOCKS,ATOMSINBLOCK)
      CALL CtoX(COORDS(:,MYNODE+1),X,Y,Z)
   ENDIF

   CALL CtoX(COORDS(:,MYNODE+1),X,Y,Z)
   !IF (DPRAND().LT.0.4) THEN
   !   DO J1=1,10
   !      Xo(:) = X(:)
   !      Yo(:) = Y(:)
   !      Zo(:) = Z(:)
   !      CALL RANDOM_STEP(X,Y,Z)
   !      call accrej_me(X,Y,Z,xo,yo,zo,BETA)
   !   ENDDO
   !ENDIF

   !CALL XtoC(COORDS(:,MYNODE+1),X,Y,Z)
   CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD,POTEL,.TRUE.,.FALSE.)
   !CALL CtoX(COORDS(:,MYNODE+1),X,Y,Z)
   !WRITE(MYUNIT,'(A,3F)'), '2 bspt_utils_amber> ',POTELold,POTEL,POTEL-POTELold

   !CALL ptmc_dumpstruct_dump(COORDS(:,MYNODE+1), POTEL, imcstep )

end subroutine bspt_takestep_amber

subroutine random_step(X,Y,Z)

   USE COMMONS
   IMPLICIT NONE
   INTEGER K
   DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),xo(NATOMS),yo(NATOMS),zo(NATOMS),RANDOM,DPRAND
   
   DO K=1,NATOMS
      RANDOM=DPRAND()
      X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
      RANDOM=DPRAND()
      Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
      RANDOM=DPRAND()
      Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
   ENDDO

end subroutine random_step

subroutine accrej_me(x,y,z,xo,yo,zo,BETA)

   USE COMMONS
   IMPLICIT NONE
   DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),xo(NATOMS),yo(NATOMS),zo(NATOMS),RANDOM,C(3*NATOMS),Co(3*NATOMS),GRAD(3*NATOMS)
   DOUBLE PRECISION V,Vo,Weight,DPRAND,BETA(0:NPAR-1)

   CALL XtoC(C,x,y,z)
   CALL XtoC(Co,xo,yo,zo)
   CALL POTENTIAL(Co,GRAD,Vo,.TRUE.,.FALSE.)
   CALL POTENTIAL(C,GRAD,V,.TRUE.,.FALSE.)

   RANDOM = DPRAND()
   Weight=MIN(1.0D0,EXP(-BETA(MYNODE)*(V-Vo)))

   !WRITE(MYUNIT, '(A,4G20.10)') 'accrej_me> VOLD,VNEW,W,RANDOM,RECOUNT=',Vo,V,Weight,RANDOM
   IF (Weight.LT.RANDOM) THEN
      ! rej
      X(:) = Xo(:)
      Y(:) = Yo(:)
      Z(:) = Zo(:)
   ELSE
      ! accept
      Xo(:) = X(:)
      Yo(:) = Y(:)
      Zo(:) = Z(:)
      !Ehsa = Ehsapert
      !Ephys = Ephyspert
   ENDIF

end subroutine accrej_me

subroutine CtoX(C,X,Y,Z)

   use commons
   implicit none

   INTEGER I
   DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),C(3*NATOMS)

   DO I=1,NATOMS
      X(I) = C(3*I-2)
      Y(I) = C(3*I-1)
      Z(I) = C(3*I)
   ENDDO

end subroutine CtoX

subroutine XtoC(C,X,Y,Z)

   use commons
   implicit none

   INTEGER I
   DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),C(3*NATOMS)

   DO I=1,NATOMS
      C(3*I-2) = X(I)
      C(3*I-1) = Y(I)
      C(3*I) = Z(I)
   ENDDO

end subroutine XtoC

subroutine bspt_takestep(X, Y, Z, IMCSTEP, &
      jumpt, peint, peqv, vnew, exab_count )
! choose which takestep routine to use and take the step
! note: MYNODE is in commons.  is it safe to use here?
   use commons
   implicit none

   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP
   logical, intent(IN) :: jumpt
   DOUBLE PRECISION, INTENT(IN) :: peint
   DOUBLE PRECISION, INTENT(IN) :: peqv(NENRPER)

   DOUBLE PRECISION, INTENT(INOUT) :: X(NATOMS), Y(NATOMS), Z(NATOMS)
   DOUBLE PRECISION, INTENT(INOUT) :: VNEW
   integer, intent(INOUT) :: exab_count

   integer j1, j2, j3, j4, k, nchosen
   double precision dummy, dprand, random

   IF (CHRMMT) THEN
      ! Random cartesian move for CHARMM
      DO K=1,NATOMS
         RANDOM=DPRAND()
         X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
         RANDOM=DPRAND()
         Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
         RANDOM=DPRAND()
         Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
      ENDDO
!
! Try jumps with no perturbations.
!
   ELSEIF (.NOT.JUMPT) THEN
      IF ( BINARY_EXAB .AND. (MOD(IMCSTEP-1.0D0,BINARY_EXAB_FRQ*1.0D0).EQ.0.0D0) ) THEN
         !js850> try to exchange a type A and type B particle as the MC step
         EXAB_COUNT = EXAB_COUNT + 1
         !choose an A and a B particle
         IF ( FREEZE ) THEN
           !get a mobile type B particle
           !NMOBILETYPEB = NTYPEB - NFREEZETYPEB = (NATOMS - NTYPEA) - (NFREEZE - NFREEZETYPEA)
           J3 = (NATOMS - NTYPEA) - (NFREEZE - NFREEZETYPEA) !NMOBILETYPEB
     !      write(*,*) NATOMS , NTYPEA,  NFREEZE , NFREEZETYPEA
           !J4=       *( NMOBILETYPEB ) +1 + NFREEZE + (NMOBILETYPEA)
           RANDOM=DPRAND()
           J4 = RANDOM*( J3           ) +1 + NFREEZE + (NTYPEA-NFREEZETYPEA)
           J2 = FROZENLIST(J4)
           !get a mobile type A particle
           !NMOBILETYPEA = NTYPEA - NFREEZETYPEA
           !J1=       *( NMOBILETYPEA ) +1 + NFREEZE
           RANDOM=DPRAND()
           J4 = RANDOM*( NTYPEA-NFREEZETYPEA ) + 1 + NFREEZE
           J1 = FROZENLIST(J4)
         ELSE
           RANDOM=DPRAND()
           J1 = RANDOM*NTYPEA + 1
           RANDOM=DPRAND()
           J2 = RANDOM*(NATOMS-NTYPEA) + 1 + NTYPEA
         ENDIF
         !WRITE(MYUNIT,*) "trying exchange ", J1, J2, IMCSTEP
         !IF ( FROZEN(J1) .OR. FROZEN(J2) .OR. J1 .GT. NTYPEA .OR. J2 .LE.  NTYPEA ) THEN
           !WRITE(*,*) "problem with binary_exab", j1, j2, frozen(j1), frozen(j2), nfreeze
           !CALL EXIT()
         !endif
         !exchange the xyz coordinates of J1 and J2
         DUMMY = X(J1)
         X(J1) = X(J2)
         X(J2) = DUMMY
         DUMMY = Y(J1)
         Y(J1) = Y(J2)
         Y(J2) = DUMMY
         DUMMY = Z(J1)
         Z(J1) = Z(J2)
         Z(J2) = DUMMY
      !ELSEIF (.NOT.((MYNODE.EQ.USERES).AND.RESERVOIRT)) THEN
      ELSE
         ! Random cartesian move if 
         !       CHRMMT=RESERVOIRT=JUMPT=FALSE & MYNODE>=USERES
         IF (DEBUGss2029) THEN  
            WRITE(MYUNIT, '(A,G20.10)') "bspt> Random cartesian move; MCSTEP = ", IMCSTEP 
         ENDIF 

         DO K=1,NATOMS
            RANDOM=DPRAND()
            X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            RANDOM=DPRAND()
            Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            RANDOM=DPRAND()
            Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
         ENDDO
      ENDIF
   ENDIF   ! closes IF (CHRMMT) THEN 

end subroutine bspt_takestep


subroutine bspt_check_configuration(X, Y, Z, BAD_CONFIG, &
      imcstep, xout, outside )
   ! do the configuration checks on the new coordinates.
   ! the result will be saved in BAD_CONFIG
   ! if any of the checks fail then BAD_CONFIG will be TRUE
   ! if they all pass, then BAD_CONFIG will be FALSE
   use commons
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP
   logical, INTENT(OUT) :: BAD_CONFIG
   DOUBLE PRECISION, INTENT(INOUT) :: X(NATOMS), Y(NATOMS), Z(NATOMS)
   DOUBLE PRECISION, INTENT(INOUT) :: XOUT
   logical, intent(INOUT) :: outside

   integer k, j1
   double precision rrx, rry, rrz, rrr
   double precision dist

   BAD_CONFIG = .FALSE.

   OUTSIDE=.FALSE.

   IF (.NOT.((MYNODE.LE.USERES).AND.RESERVOIRT)) THEN
      ! reservoir stuff 
      IF (.NOT.(CHRMMT.OR.MODEL1T.OR.PERCOLATET.OR.PERIODIC)) THEN
         cloop: DO K=1,NATOMS
            DIST=X(K)**2+Y(K)**2+Z(K)**2
            IF (DIST.GT.RADIUS) THEN
               IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) WRITE(MYUNIT,'(A,I6,A)') 'bspt> Perturbed atom ',K, &
                    ' outside container, recount previous configuration results'
               XOUT=XOUT+1.0D0
               OUTSIDE=.TRUE.
               BAD_CONFIG=.TRUE. ! The right way to deal with rejected steps!
               EXIT cloop
            ENDIF
         ENDDO cloop
      ENDIF

      IF (PERCOLATET) THEN
         CALL PERC(COORDS(1:3*NATOMS,MYNODE+1),NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)
         IF (.NOT.PERCT) THEN
            IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) WRITE(MYUNIT,'(A,I6,A)') &
                 'bspt> After step system is not a percolating network, recount previous configuration results'
            XOUT=XOUT+1.0D0
            OUTSIDE=.TRUE.
            BAD_CONFIG=.TRUE. ! The right way to deal with rejected steps!
         ENDIF
      ENDIF
!
! js850> Check to see if RESTRICTREGION constraints have been violated - if they have, reject the step
!
      IF ( RESTRICTREGION ) THEN
        !RESTRICTREGIONTEST=.FALSE.
        DO J1=1,NATOMS
          IF ( .NOT. FROZEN(J1) .AND. .NOT. HARMONICFLIST(J1) &
       .AND. .NOT. DONTMOVE(J1) ) THEN
            RRX = ( X(J1)-RESTRICTREGIONX0 )
            RRY = ( Y(J1)-RESTRICTREGIONY0 )
            RRZ = ( Z(J1)-RESTRICTREGIONZ0 )
            RRX = RRX - ANINT(RRX/BOXLX)*BOXLX
            RRY = RRY - ANINT(RRY/BOXLY)*BOXLY
            IF ( RESTRICTCYL ) THEN
              RRZ = 0
              RRR = DSQRT(RRX**2+RRY**2 )
            ELSE
              RRZ = RRZ - ANINT(RRZ/BOXLZ)*BOXLZ
              RRR = DSQRT(RRX**2+RRY**2+RRZ**2 )
            ENDIF
            IF ( RRR > RESTRICTREGIONRADIUS ) THEN
              BAD_CONFIG=.TRUE.
              WRITE(MYUNIT, *) 'bspt> restrictregion> rejecting step ', IMCSTEP
              !GOODSTRUCTURE=.FALSE.
              !RESTRICTREGIONTEST=.TRUE.
              EXIT !exit loop
            ENDIF
          ENDIF
        ENDDO
      ENDIF
   ENDIF ! closes IF (.NOT.((MYNODE.LE.USERES).AND.RESERVOIRT))
end subroutine bspt_check_configuration


SUBROUTINE BSPT_READ_BSPTRESTART( XYZ, PREVSTEPS, VOLD, VMINOLD, &
      NACCEPTPT, NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN )
   ! For restart we need to get the current configuration, its pe, the pe of the minimum it quenched to,
   ! if applicable, the number of steps already done, the maximum step size, and the Visits and Visits2
   ! histograms. If we dump using BSPTDUMPFRQ then we can restore from the last such file. We can work
   ! out what the last dump was once we know how many steps have been done!
   use commons, only : mynode, natoms, nq, step, npar, myunit
   implicit none
   double precision, intent(OUT) :: xyz(3*natoms), PREVSTEPS, vold, vminold, NACCEPTPT(0:NPAR-1)
   integer, intent(OUT) :: NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN 
   CHARACTER (LEN=256)  ISTR, SDUMMY, filename12
   integer GETUNIT
   integer lunit, j1, j2

   !
   !read file [mynode+1]/bsptrestart
   !
   WRITE (ISTR, '(I10)') MYNODE+1
   FILENAME12=TRIM(ADJUSTL(ISTR)) // "/bsptrestart"
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME12, STATUS="old", form="formatted")
   READ(LUNIT,*) PREVSTEPS,VOLD,VMINOLD,STEP(MYNODE+1), &
      NACCEPTPT(MYNODE),NEACCEPT,NTOT,NOUTQBIN,NOUTPEBIN,NQ(MYNODE+1)
   WRITE(MYUNIT,'(A,A)') 'bspt> Reading restart information from ',TRIM(ADJUSTL(FILENAME12))
   WRITE (MYUNIT,'(A,F20.1)')      'bspt> Number of previous steps=     ',PREVSTEPS
   WRITE (MYUNIT,'(A,2F20.10)')    'bspt> Quench and instantaneous PE=',VMINOLD,VOLD
   WRITE (MYUNIT,'(A,F20.10)')     'bspt> Step size=',STEP(MYNODE+1)
   WRITE (MYUNIT,'(A,F15.1,2I10)') 'bspt> Accepted MC steps, PT steps and total PT steps=',NACCEPTPT(MYNODE),NEACCEPT,NTOT
   WRITE (MYUNIT,'(A,3I10)')       'bspt> Quenches and instantaneous energies outside range=',NOUTQBIN,NOUTPEBIN
   WRITE (MYUNIT,'(A,I10)')        'bspt> Total quenches=',NQ(MYNODE+1)

   DO J1=1,NATOMS
      J2=3*(J1-1)
      READ(LUNIT,*) xyz(J2+1), xyz(J2+2), xyz(J2+3)
   ENDDO
   CLOSE(LUNIT)

end subroutine bspt_read_bsptrestart

subroutine bspt_read_vists_his( PREVSTEPS, pevisits, qvisits, T )
   ! For restart we need to get the current configuration, its pe, the pe of the minimum it quenched to,
   ! if applicable, the number of steps already done, the maximum step size, and the Visits and Visits2
   ! histograms. If we dump using BSPTDUMPFRQ then we can restore from the last such file. We can work
   ! out what the last dump was once we know how many steps have been done!
   USE PORFUNCS
   USE COMMONS, ONLY : MYNODE, BSPTDUMPFRQ, NENRPER, HBINS, NPAR, MYUNIT, BSPT
   implicit none
   double precision, intent(IN) :: PREVSTEPS
   DOUBLE PRECISION, intent(out) :: QVISITS(HBINS, 0:NPAR-1), PEVISITS(NENRPER,0:NPAR-1)
   double precision, intent(out) :: T
   CHARACTER (LEN=256)  ISTR, SDUMMY
   CHARACTER (LEN=256)  FILENAME101, FILENAME1, FILENAME2, FILENAME3, FILENAMEPREFIX
   integer lunit, j1, j2, k
   double precision dummy
   logical FILETEST
   integer GETUNIT

   !
   ! read Visits.his file:  it could have various names.  First look for
   ! FILENAME1 = [MYNODE+1]/Visits.his.[IMCSTEP]
   ! if that doesn't exists, then look for
   ! FILENAME2 = [mynode+1]/Visits.his.restart
   ! if that doesn't exists, then look for
   ! FILENAME3 = [mynode+1]/Visits.his
   !
   ! the final file name will be in FILENAME101
   !
   WRITE (ISTR, '(I10)') MYNODE+1
   filenameprefix=TRIM(ADJUSTL(ISTR)) // "/Visits.his"
   FILENAME2=TRIM(ADJUSTL(filenameprefix)) // ".restart"
   FILENAME3=TRIM(ADJUSTL(filenameprefix))
   IF (BSPTDUMPFRQ.GT.0) THEN
      DUMMY=INT(PREVSTEPS/(1.0D0*BSPTDUMPFRQ))*1.0D0
      DUMMY=DUMMY*BSPTDUMPFRQ
      WRITE (SDUMMY, '(F15.1)') DUMMY
      FILENAME101=TRIM(ADJUSTL(filenameprefix)) // '.' // TRIM(ADJUSTL(SDUMMY))
      !
      !js850> if file #/Visits.his.#####.0 doesn't exist, then check if
      !          file #/Visits.his.restart  exists and use that instead.  This is
      !          a way to change BSPTDUMPFRQ after a restart
      !

      WRITE(MYUNIT,'(A,A)') 'bspt> Looking for file ',TRIM(ADJUSTL(FILENAME101))
      INQUIRE(FILE=TRIM(ADJUSTL(FILENAME101)), EXIST=FILETEST)
   ELSE
      FILETEST = .FALSE.
   ENDIF
   IF ( .NOT. FILETEST ) THEN
      WRITE(MYUNIT,'(A,A)') 'bspt> Looking for file ',TRIM(ADJUSTL(FILENAME2))
      INQUIRE(FILE=FILENAME2, EXIST=FILETEST)
      IF ( FILETEST ) THEN 
         FILENAME101 = TRIM(ADJUSTL(FILENAME2)) 
      ELSE
         WRITE(MYUNIT,'(A,A)') 'bspt> Looking for file ',TRIM(ADJUSTL(FILENAME3))
         FILENAME101 = TRIM(ADJUSTL(FILENAME3)) 
      ENDIF
   ENDIF

   !
   ! Now we have the file name in FILENAME101.  Load the data
   !
   WRITE(MYUNIT,'(A,A)') 'bspt> Reading restart information from ',TRIM(ADJUSTL(FILENAME101))
   CALL FLUSH(MYUNIT)
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME101, STATUS="unknown", form="formatted")
   READ(LUNIT, '(G20.10)') T !this is the temperature
   READ(LUNIT, '(A)') SDUMMY
   DO K=1, NENRPER
      READ(LUNIT,*) DUMMY, PEVISITS(K,MYNODE)
   ENDDO
   IF (BSPT) THEN
      READ(LUNIT, '(A)') SDUMMY
      DO K=1, HBINS
         WRITE(LUNIT,*) DUMMY, QVISITS(K,MYNODE)
      ENDDO
   ELSE
      QVISITS(:,MYNODE) = 0.
   ENDIF
   CLOSE(LUNIT)
end subroutine bspt_read_vists_his

subroutine bspt_read_vists_his2( pevisits2 )
   ! read pevisits2 from file 
   ! [mynode+1]/Visits2.his
   use commons, only : mynode, nenrper, hbins, npar, myunit
   implicit none
   double precision, intent(out) :: PEVISITS2(NENRPER, HBINS, 0:NPAR-1)
   CHARACTER (LEN=256)  ISTR
   CHARACTER (LEN=256)  FILENAME101
   integer lunit
   integer GETUNIT
   FILENAME101=TRIM(ADJUSTL(ISTR)) // "/Visits2.his"
   WRITE(MYUNIT,'(A,A)') 'bspt> Reading restart information from ',TRIM(ADJUSTL(FILENAME101))
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME101,STATUS="unknown",FORM="UNFORMATTED")
   READ(LUNIT) PEVISITS2(1:NENRPER,1:HBINS,MYNODE)
   CLOSE(LUNIT) 
end subroutine bspt_read_vists_his2

SUBROUTINE PTMC_DUMP_HISTOGRAM(PTEMIN, PEINT, PEVISITS, binlabel, qvisits, &
      TEMPTRAJ, IMCSTEP, FINALT)
   ! write the visits histogram of energies and quench energies to file n/Visits.his.#
   ! where n is the node number and # is the current step number.
   ! Note: the step number is written as a float, so it includes
   ! the first decimal place.  e.g. 10000.0 
   !
   ! if FINALT is true, then print to file n/Visits.his
   !
   ! ss2029 > note: suffix of Visits file is IMCSTEP, but it
   ! includes only post-equilibration data 
   USE COMMONS, ONLY : MYNODE, NPAR, NENRPER, HBINS, BSPT
   implicit none
   double precision, intent(IN) :: imcstep
   LOGICAL, INTENT(IN) :: FINALT
   DOUBLE PRECISION, INTENT(IN) :: PTEMIN, PEINT, PEVISITS(NENRPER,0:NPAR-1)
   DOUBLE PRECISION, INTENT(IN) :: TEMPTRAJ(0:NPAR-1)
   DOUBLE PRECISION, INTENT(IN) :: binlabel(hbins), QVISITS(HBINS,0:NPAR-1)
   CHARACTER (LEN=256)  ISTR, SDUMMY, filename101
   integer lunit, k
   integer getunit

   ! set up the file name to print to
   WRITE (ISTR, '(I2,A1)') MYNODE+1
   if (finalt) then
      ! [mynode+1]/Visits.his
      FILENAME101=trim(adjustl(istr)) // "/Visits.his"
   else
      ! [mynode+1]/Visits.his.[imcstep]
      WRITE (SDUMMY, '(F15.1)') IMCSTEP
      ISTR=TRIM(ADJUSTL(ISTR)) // '/Visits.his.' // TRIM(ADJUSTL(SDUMMY))
      FILENAME101=TRIM(ADJUSTL(ISTR))
   endif

   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME101, STATUS="unknown", form="formatted")
   WRITE(LUNIT, '(G20.10)') TEMPTRAJ(MYNODE)
   WRITE(LUNIT, '(A)') 'Visits to instantaneous PE bins without quench contributions'
   DO K=1, NENRPER
      WRITE(LUNIT, '(G20.10,F20.1)') PTEMIN+(K-1)*PEINT,PEVISITS(K,MYNODE)
   ENDDO
   !js850> this should be hidden behind bspt. but we must be careful
   !       when reading visits.his files on restart
   IF (BSPT) THEN
      WRITE(LUNIT, '(A)') 'Visits to quench bins'
      DO K=1, HBINS
         WRITE(LUNIT, '(G20.10,F20.1)') BINLABEL(K), QVISITS(K,MYNODE)
      ENDDO
   ENDIF
   CLOSE(LUNIT)

end subroutine ptmc_dump_histogram

subroutine bspt_dump_histogram2(pevisits2, imcstep, finalt)
   ! dump 2d-array PEVISITS2 to [mynode+1]/Visits2.his.[imcstep].0
   use commons, only : nenrper, hbins, npar, mynode
   implicit none
   integer, intent(IN) :: PEVISITS2(NENRPER, HBINS, 0:NPAR-1)
   double precision, intent(IN) :: imcstep
   LOGICAL, INTENT(IN) :: FINALT
   CHARACTER (LEN=256)  ISTR, SDUMMY, filename101
   integer lunit
   integer getunit

   ! set up the file name to print to
   WRITE (ISTR, '(I2,A1)') MYNODE+1
   if (finalt) then
      ! [mynode+1]/Visits2.his
      FILENAME101=trim(adjustl(istr)) // "/Visits2.his"
   else
      ! [mynode+1]/Visits2.his.[imcstep]
      WRITE (SDUMMY, '(F15.1)') IMCSTEP
      ISTR=TRIM(ADJUSTL(ISTR)) // '/Visits2.his.' // TRIM(ADJUSTL(SDUMMY))
      FILENAME101=TRIM(ADJUSTL(ISTR))
   endif

   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME101,STATUS="unknown",FORM="UNFORMATTED")
   WRITE(LUNIT) PEVISITS2(1:NENRPER,1:HBINS,MYNODE)
   CLOSE(LUNIT)
end subroutine bspt_dump_histogram2

SUBROUTINE PTMC_DUMP_RESTART_INFO(X, Y, Z, IMCSTEP, VOLD, VMINOLD, &
      NACCEPTPT, NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN )
   ! dump restart information to file n/bsptrestart
   ! where n is the node number.
   ! but first copy n/bsptrestart to n/bsptrestart.save
   USE PORFUNCS
   USE COMMONS, ONLY : NATOMS, MYNODE, NPAR, NQ, STEP, PTMC, BSPT
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NEACCEPT, NTOT, NOUTQBIN, NOUTPEBIN
   double precision, INTENT(IN) :: NACCEPTPT(0:NPAR-1)
   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP
   DOUBLE PRECISION, INTENT(INOUT) :: X(NATOMS), Y(NATOMS), Z(NATOMS)
   DOUBLE PRECISION, INTENT(INOUT) :: VOLD, VMINOLD
   CHARACTER (LEN=256)  ISTR, SDUMMY, filename9, FILENAME10
   integer lunit, j1
   integer getunit

   WRITE (ISTR,'(I2)') MYNODE+1
   FILENAME9=TRIM(ADJUSTL(ISTR)) // "/bsptrestart"
   !js850> copy #/bsptrestart to #/bsptrestart.save
   FILENAME10=TRIM(ADJUSTL(ISTR)) // "/bsptrestart.save"

   ! js850> note: if the " &> /dev/null " is included in the "cp"
   ! statement below then sometimes the file is copied before
   ! the shell finished writing from the previous time this
   ! subroutine was called.  I don't know why removing the " &> ! /dev/null" 
   ! helps, but it does.
   SDUMMY="cp -p "//TRIM(ADJUSTL(FILENAME9))//" "//TRIM(ADJUSTL(FILENAME10))!//" &> /dev/null"
   CALL SYSTEM(SDUMMY)
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=FILENAME9, STATUS="unknown", form="formatted")
   IF (PTMC .AND. .NOT. BSPT) VMINOLD = 0.D0 !TO PREVENT IT FROM PRINTING AS ***********
   WRITE(LUNIT,'(F20.1,3G20.10,F15.1,5I15)') IMCSTEP,VOLD,VMINOLD,STEP(MYNODE+1),&
                                       NACCEPTPT(MYNODE),NEACCEPT,NTOT,NOUTQBIN,NOUTPEBIN,NQ(MYNODE+1)
   DO J1=1,NATOMS
      WRITE(LUNIT,'(3G25.15)') X(J1),Y(J1),Z(J1)
   ENDDO
   CALL FLUSH(LUNIT)
   CLOSE(LUNIT)
end subroutine ptmc_dump_restart_info

subroutine ptmc_dumpstruct_dump(xyz, vnew, imcstep )
   use commons, only : bsptrestart, natoms, mynode
   implicit none
   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP, VNEW, XYZ(3*NATOMS)
   LOGICAL, SAVE :: first_call = .true.
   integer tempunit, GETUNIT, j1
   CHARACTER (LEN=256)  SDUMMY, filename106

   ! create dumpstruct.repname 
   WRITE(FILENAME106,*) 'dumpstruct.'
   WRITE(SDUMMY,'(I3)') MYNODE+1 
   FILENAME106=TRIM(ADJUSTL(FILENAME106)) // TRIM(ADJUSTL(SDUMMY))

   !ss2029> create replicaNumber/dumpstruct instead of dumpstruct.replicaNumber 
   !
   !WRITE(SDUMMY,'(I3)') MYNODE+1 
   !FILENAME106=TRIM(ADJUSTL(SDUMMY)) // '/dumpstruct' 

   !
   !  ss2029> dumpstruct.# files are being opened and closed within this main loop. Might be better to 
   !       open and close only once outside the loop as done for
   !       FILENAME_ENER above. todo 
   !
   TEMPUNIT=GETUNIT()
   IF ( first_call .AND. (.not.BSPTRESTART) ) THEN 
   OPEN(UNIT=TEMPUNIT,FILE=FILENAME106, STATUS="unknown", form="formatted")
   first_call = .FALSE.
   ELSE
   OPEN(UNIT=TEMPUNIT,FILE=FILENAME106, STATUS="unknown", form="formatted", POSITION="APPEND")
   END IF
   WRITE(TEMPUNIT,*) NATOMS
   WRITE(TEMPUNIT,*) VNEW, IMCSTEP
   DO J1=1,NATOMS
   WRITE(TEMPUNIT,'(A,3F20.10)') 'LA  ', XYZ( 3*(J1-1)+1), &
     XYZ( 3*(J1-1)+2), XYZ( 3*(J1-1)+3)
   END DO
   CLOSE(TEMPUNIT)
end subroutine ptmc_dumpstruct_dump

subroutine bspt_do_quenching(x, y, z, xo, yo, zo, mincoords, imcstep, vminnew,&
   vminold, LBFGS_ITERATIONS, recount, histint, NOUTQBIN, beta, outside )
   ! do the quenching part of basin sampling
   ! 
   ! 1. quench the current monte carlo chain configuration. Note that if this
   ! step was rejected we need to quench the structure in xO, yO, zO.
   ! 
   ! 2. if the quench energy is out of bounds (of this histogram), then reject 
   ! the step by setting RECOUNT to .TRUE.
   !
   ! the quenched structure will be in MINCOORDS
   use commons
   !use commons, only: bspt, mindensityT, RESERVOIRT, quenchfrq, mynode, natoms,&
      !nq, debug, npar, histmin, hbins, CHRMMT, myunit
   USE MODCHARMM
   implicit none
   DOUBLE PRECISION XYZ(3*NATOMS)
   DOUBLE PRECISION, INTENT(IN) :: X(NATOMS), Y(NATOMS), Z(NATOMS)
   DOUBLE PRECISION, INTENT(IN) :: XO(NATOMS), YO(NATOMS), ZO(NATOMS)
   DOUBLE PRECISION, INTENT(IN) :: BETA(0:NPAR-1)
   DOUBLE PRECISION, INTENT(IN) :: IMCSTEP
   logical, intent(IN) :: outside
   double precision, intent(IN) :: histint

   DOUBLE PRECISION, INTENT(INOUT) :: MINCOORDS(3*NATOMS,NPAR)
   DOUBLE PRECISION, INTENT(INOUT) :: vminnew, vminold

   integer, intent(INOUT) :: NOUTQBIN, LBFGS_ITERATIONS
   LOGICAL, INTENT(INOUT) :: RECOUNT

   LOGICAL EVAPREJECT, QUENCHEDT
   integer converged
   integer binindex, ndummy
   double precision dummy, random, dprand, wcomp, w

   DOUBLE PRECISION potel
   COMMON /MYPOT/ POTEL

   ! js850> these are not always initialized below
   CONVERGED = 1
   EVAPREJECT = .FALSE.
   QUENCHEDT = .FALSE.

   ! Quenching part if required.

   ! COORDSO saves the perturbed coordinates before the quench in order to calculate
   ! a quench distance. Should no longer be needed. COORDS are used as scratch for quenches.
   !COORDSO(:,MYNODE+1)=COORDS(:,MYNODE+1) 

   !js850> what is this for?
   IF ((RECOUNT.AND.BSPT).OR.(IMCSTEP.LE.NEQUIL+PTSTEPS)) THEN
      VMINNEW=0.0D0
      LBFGS_ITERATIONS=0
   ENDIF

   ! ss2029> advanced BSPT to outermost IF
   IF (BSPT.AND.(MINDENSITYT.OR.(IMCSTEP.GT.NEQUIL+PTSTEPS))) THEN
      
      IF (.NOT.RECOUNT) THEN 
         IF (MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.1.0D0*0) THEN 
            !
            ! do the quenching
            !
            IF ((MYNODE.LE.USERES).AND.RESERVOIRT) THEN
               WRITE(MYUNIT,'(A)') 'bspt> Reservoir not coded yet for BSPT - quit'
               VMINNEW=0.0D0
               LBFGS_ITERATIONS=0
               EVAPREJECT=.FALSE.
               CONVERGED=1
               POTEL=HUGE(1.0D0)
               STOP
            ELSE
               ! copy X, Y, Z into XYZ
               XYZ(1:3*NATOMS-2:3) = X(:)
               XYZ(2:3*NATOMS-1:3) = Y(:)
               XYZ(3:3*NATOMS-0:3) = Z(:)
               !WRITE(MYUNIT, '(A)') 'bspt> calling quench'
               IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
               CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,XYZ)
               IF (CONVERGED.NE.1) WRITE(MYUNIT, '(A)') 'bspt> WARNING - quench did not converge' 
               QUENCHEDT = .TRUE.
               VMINNEW=POTEL
               NQ(MYNODE+1)=NQ(MYNODE+1)+1
               MINCOORDS(:,MYNODE+1)=XYZ(:) ! MINCOORDS contains quench coords
            ENDIF
         ELSE
            !
            ! not doing the quenching this step
            !
            VMINNEW=0.0D0
            !WRITE(MYUNIT, '(A)') 'bspt> not quenching'
            LBFGS_ITERATIONS=0
            EVAPREJECT=.FALSE.
            CONVERGED=1
            POTEL=HUGE(1.0D0)
         ENDIF
      ELSEIF (RECOUNT.AND.(VMINOLD.EQ.0.0D0).AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.0.0D0)) THEN
         !
         ! We might not have quenched at the step to be recounted, so we have to do so here.
         !
         IF (DEBUG) WRITE(MYUNIT, '(A,G20.10,A)') 'bspt> recounting step for previous configuration with VMINOLD=',&
                                                   VMINOLD,' need to call quench'
         ! copy X, Y, Z into XYZ
         XYZ(1:3*NATOMS-2:3) = XO(:)
         XYZ(2:3*NATOMS-1:3) = YO(:)
         XYZ(3:3*NATOMS-0:3) = ZO(:)
         !WRITE(MYUNIT, '(A)') 'bspt> calling quench for recount'
         IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,XYZ)
         QUENCHEDT = .TRUE.
         NQ(MYNODE+1)=NQ(MYNODE+1)+1
         MINCOORDS(:,MYNODE+1)=XYZ(:) ! MINCOORDS contains quench coords
         VMINOLD=POTEL
         VMINNEW = POTEL ! js850> I added this.  is it supposed to be here?
         IF (CONVERGED.NE.1) WRITE(MYUNIT, '(A)') 'bspt> WARNING - quench did not converge' 
         IF (DEBUG) WRITE(MYUNIT, '(A,G20.10)') 'bspt> VMINOLD set to ',VMINOLD
      ENDIF
   ENDIF ! closes IF (MINDENSITYT.OR.(IMCSTEP.GT.NEQUIL+PTSTEPS))
!
! If either the PE (VNEW) or quench energy (VMINNEW) is out of range of the histogram
! then we may reject the step and recount the previous one.
!
! BININDEX for quench bin
! IBININDEX for instantaneous bin
!
   IF (QUENCHEDT) THEN
      IF (BSPT.AND.((IMCSTEP.GT.NEQUIL+PTSTEPS).OR.(MINDENSITYT)).AND.(MOD(IMCSTEP,1.0D0*QUENCHFRQ).EQ.0.0D0)) THEN
         BININDEX=INT((VMINNEW-HISTMIN)/HISTINT)+1
         IF ((BININDEX.GT.HBINS).OR.(BININDEX.LT.1)) THEN
            RECOUNT=.TRUE.
            ! I think this flag OUTSIDE is wrong.  if it was outside then RECOUNT = .true. and the quench would have been done on the previous configuration, one that was not outside
            IF (.NOT.OUTSIDE) NOUTQBIN=NOUTQBIN+1
            !WRITE(MYUNIT, '(A,G20.10)') 'bspt> WARNING: quench energy out of bounds ', VMINNEW, histmin, histint, binindex, hbins
            WRITE(MYUNIT, *) 'bspt> WARNING: quench energy out of bounds ', VMINNEW, histmin, histint, binindex, hbins
         ENDIF
      ENDIF
      IF (MINDENSITYT) THEN
         WCOMP=(VMINNEW-VMINOLD)*BETA(MYNODE) ! use difference in quench energies for MINDENSITYT
         W=MIN(1.0D0,EXP(-WCOMP))
         RANDOM=DPRAND()
         IF (RANDOM.GT.W) RECOUNT=.TRUE. 
      ENDIF

      !IF (DEBUGss2029) THEN
          !IF (MYNODE.EQ.0) print *,'bspt> ss2029: three MC step = ', IMCSTEP,' RECOUNT = ', RECOUNT 
      !ENDIF

      IF ((CONVERGED.NE.1).OR.EVAPREJECT) RECOUNT=.TRUE. ! reject and recount
   endif
end subroutine bspt_do_quenching


