SUBROUTINE DUMPSTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
!USE COMMONS,ONLY : NATOMS, NPAR, STEP, ASTEP, OSTEP, TEMP, NMSBSAVE, MSBE, MSBCOORDS, &
!  &                NSAVE, MAXSAVE, MYUNIT, A9INTET, NSAVEINTE, COORDSO, &
!  &                DUMPSTRUCTURES, DUMPMINT, NPCALL, NATOMSALLOC !jdf43>
USE COMMONS,ONLY : NATOMS, NPAR, STEP, ASTEP, OSTEP, TEMP, NMSBSAVE, MSBE, MSBCOORDS, &
  &                NSAVE, MAXSAVE, MYUNIT, A9INTET, NSAVEINTE, &
  &                DUMPSTRUCTURES, DUMPMINT, NPCALL, NATOMSALLOC, & !jdf43>
  &                GCBHT, GCNATOMS, QENERGIES, QPE, QCOORDINATES, AVOIDNATOMS

USE QMODULE
USE PORFUNCS
USE OUTPUT, ONLY : WRITE_MARKOV_COORDS
IMPLICIT NONE
INTEGER NDONE, JBEST(NPAR), JP, J1, MYUNIT2, GETUNIT, LUNIT
DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMSALLOC,NPAR)
DOUBLE PRECISION  DUMGRAD(3*NATOMSALLOC), OPOTEL, SAVEP
CHARACTER(LEN=20) :: ISTR
! ss2029> added variable ISTR1 
CHARACTER(LEN=20) :: ISTR1
LOGICAL LOPEN
DOUBLE PRECISION POTEL
COMMON /MYPOT/ POTEL

IF (NPAR.GT.1) THEN
   WRITE (ISTR, '(i10)') JP
   MYUNIT2=GETUNIT()
   CALL FLUSH(MYUNIT)
!   CALL SYSTEM('cp "GMIN.dump."//trim(adjustl(istr)) "GMIN.dump."//trim(adjustl(istr))//".save"')
   OPEN(MYUNIT2,FILE="GMIN.dump."//trim(adjustl(istr)),STATUS='UNKNOWN')
ELSE
   MYUNIT2=GETUNIT()
   WRITE (ISTR, '(A)') ''
   CALL SYSTEM('cp GMIN.dump GMIN.dump.save')
   OPEN(MYUNIT2,FILE='GMIN.dump',STATUS='UNKNOWN')
ENDIF

  WRITE(MYUNIT2, '(A)') 'steps completed NQ(JP) in mc'
  WRITE(MYUNIT2, '(I8)') NDONE
  WRITE(MYUNIT2, '(I8)') NSAVE
! Number of potential energy calls done so far
  WRITE(MYUNIT2, '(I20)') NPCALL 
  WRITE(MYUNIT2, '(A)') 'COORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  WRITE(MYUNIT2, '(I8)') NATOMS
  CALL WRITE_MARKOV_COORDS(MYUNIT2, '(3F25.15)', JP)
!  WRITE(MYUNIT2, '(3F25.15)') COORDSO(1:3*NATOMS,JP)
  WRITE(MYUNIT2, '(A)') 'STEP, ASTEP, OSTEP, TEMP:'
  WRITE(MYUNIT2, '(4F25.15)') STEP(JP),ASTEP(JP),OSTEP(JP),TEMP(JP)
  WRITE(MYUNIT2, '(A)') 'QMIN and QMINP'
! Dump info for each saved minimum
  DO J1=1,NSAVE
     WRITE(MYUNIT2, '(A,I8)') 'saved minimum ',J1
! Energy of saved minimum J1
     WRITE(MYUNIT2, '(G25.15)') QMIN(J1)
! Number of potential energy calls taken when this minimum was first encountered
     WRITE(MYUNIT2, '(I20)') NPCALL_QMIN(J1)
! Coordinates of saved minimum J1
     WRITE(MYUNIT2, '(I20)') QMINNATOMS(J1)
     WRITE(MYUNIT2, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1))
! ss2029> Write out minima from each replica (coded with BHPT output in mind) 
!         in xyz format. 27/02/2012 
     IF (QMIN(J1).LT.1.0D10) THEN   
        WRITE (ISTR1, '(i10)') J1 
        LUNIT=GETUNIT()
        IF ((NPAR.GT.1).AND.(DUMPMINT)) THEN
          OPEN(LUNIT,FILE="dumpmin."//TRIM(ADJUSTL(ISTR))//"."//TRIM(ADJUSTL(ISTR1)),STATUS='UNKNOWN')
          WRITE(LUNIT, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1)) 
          CLOSE(LUNIT)
!        ELSE !jdf43> 
        ELSEIF (DUMPSTRUCTURES.AND.DUMPMINT) THEN !jdf43>
          OPEN(LUNIT,FILE="dumpmin."//TRIM(ADJUSTL(ISTR1)),STATUS='UNKNOWN')
          WRITE(LUNIT, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1)) 
          CLOSE(LUNIT) 
        ENDIF
     ENDIF  
  ENDDO
  IF (GCBHT) THEN
     WRITE(MYUNIT2, '(A)') 'for GCBH QENERGIES, QPE, QCOORDINATES'
     DO J1=1,GCNATOMS
        WRITE(MYUNIT2, '(A,I8)') 'number of atoms',J1
        WRITE(MYUNIT2, '(2G25.15)') QENERGIES(J1),QPE(J1)
        WRITE(MYUNIT2, '(3F25.15)') QCOORDINATES(J1,1:3*J1)
     ENDDO
  ENDIF
  WRITE(MYUNIT2,'(A)') 'new restart procedure - JBEST, EBEST, BESTCOORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  WRITE(MYUNIT2, '(I8,F25.15)') JBEST(JP), EBEST(JP)
  WRITE(MYUNIT2, '(I8)') QMINNATOMS(1)
  WRITE(MYUNIT2, '(3F25.15)') BESTCOORDS(1:3*QMINNATOMS(1),JP)
! bs360: the following is not tested for MPI (July 2010)
  WRITE(MYUNIT2, '(A)') 'new restart procedure - saved energies and coordinates in MSB list'
  WRITE(MYUNIT2, '(I8)') NMSBSAVE
  DO J1=1,MIN(NMSBSAVE,MAXSAVE)
     WRITE(MYUNIT2, '(A,I8)') 'structure number ',J1
     WRITE(MYUNIT2, '(G25.15)') MSBE(J1)
     WRITE(MYUNIT2, '(I8)') AVOIDNATOMS(J1)
     WRITE(MYUNIT2, '(3F25.15)') MSBCOORDS(1:3*AVOIDNATOMS(J1),JP)
  ENDDO
! WRITE(MYUNIT2, '(A)','old taboo procedure: energy and inertia determinant'
! DO J1=1,NPAR
   ! WRITE(MYUNIT2, '(A,I8)') 'taboo structure ',J1
   ! WRITE(MYUNIT2, '(2G25.15)') ESAVE(J1), XINSAVE(J1)
! ENDDO
!
  CLOSE(MYUNIT2)


!   csw34> Add in dumping of the lowest interaction energies if A9INTE is specified
IF (NPAR.EQ.1) THEN
  IF (A9INTET) THEN
     CALL SYSTEM('cp GMINinte.dump GMINinte.dump.save')

     LUNIT=GETUNIT()
     OPEN(UNIT=LUNIT,FILE='GMINinte.dump',STATUS='UNKNOWN')
     WRITE(LUNIT, '(A)') 'SAVEINTE'
     WRITE(LUNIT, '(I8)') NSAVEINTE
     WRITE(LUNIT, '(A)') 'steps completed J1 in mc'
     WRITE(LUNIT, '(I8)') NDONE
     WRITE(LUNIT, '(A)') 'INTEQMIN and INTEQMINP'
     DO J1=1,NSAVEINTE
        WRITE(LUNIT, '(A,I8)') 'saved interaction enthalpy minimum ',J1
        WRITE(LUNIT, '(G25.15)') INTEQMIN(J1)
        WRITE(LUNIT, '(3F25.15)') INTEQMINP(J1,1:3*NATOMS)
     ENDDO
     CLOSE(UNIT=LUNIT)
  ENDIF
ENDIF

RETURN
END SUBROUTINE DUMPSTATE

SUBROUTINE RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
USE COMMONS,ONLY : NATOMS, COORDS, NPAR, STEP, ASTEP, OSTEP, TEMP, NMSBSAVE, MSBE, MSBCOORDS, &
  &                NSAVE, MAXSAVE, MYUNIT, AMHT, SEEDT, A9INTET, NSAVEINTE, INTEDUMPFILE, INTERESTORE, NPCALL, GCBHT, &
  &                GCNATOMS, QENERGIES, QPE, QCOORDINATES, AVOIDNATOMS
USE QMODULE
IMPLICIT NONE
INTEGER NDONE, JBEST(NPAR), JP, J1, OLDSAVE, MYUNIT2, LUNIT, GETUNIT, NDUMMY
DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), DUMMYARRAY(3*NATOMS)
CHARACTER(LEN=1) DUMMYS
CHARACTER(LEN=20) :: ISTR

IF (NPAR.GT.1) THEN
   WRITE (ISTR, '(i10)') JP
   MYUNIT2=GETUNIT()
   OPEN(MYUNIT2,FILE="GMIN.dump."//trim(adjustl(istr)),STATUS='OLD')
ELSE
   MYUNIT2=1
   OPEN(MYUNIT2,FILE='GMIN.dump',STATUS='OLD')
ENDIF

   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,'(I8)') NDONE
   READ(MYUNIT2,'(I8)') OLDSAVE
   READ(MYUNIT2,'(I20)') NPCALL 
   READ(MYUNIT2,*) DUMMYS
! Read in last minimum, which is also the last minimum in the Markov chain.
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) NATOMS
   READ(MYUNIT2,*) COORDS(1:3*NATOMS,JP)
   DO J1=1,NATOMS*(NPAR-1)
      READ(MYUNIT2,*)
   ENDDO
   IF ((.NOT.SEEDT).AND.(.NOT.AMHT)) THEN
      WRITE(MYUNIT,'(A,I4)') 'Initial coordinates: process',JP
      WRITE(MYUNIT,'(3F20.10)') (COORDS(J1,JP),J1=1,3*NATOMS)
   ENDIF
! step and temperature information
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) STEP(JP),ASTEP(JP),OSTEP(JP),TEMP(JP)
! best OLDSAVE minima
  READ(MYUNIT2,*) DUMMYS
!       csw34> this is the loop that has been edited to allow SAVE to vary between
!       runs. There are three cases which must be considered. 1) NSAVE=OLDSAVE.
!       In this case, nothing new needs to be done. 2) NSAVE<OLDSAVE. Now, the
!       extra minima that are saved in the dump file need to be cycled over
!       (ignored). 3) NSAVE>OLDSAVE. The array must be allocated with
!       dimensions suitable for the new number of minima, and the dump file read
!       in as far as possible. The rest of the array must then be padded with
!       zeros. The QMIN and QMINP arrays are already allocated in main.F with
!       the correct dimension so that need not worry us here :)
        IF (NSAVE.EQ.OLDSAVE) THEN
           DO J1=1,NSAVE
              READ(MYUNIT2,*) DUMMYS
              READ(MYUNIT2,*) QMIN(J1)
              READ(MYUNIT2,*) NPCALL_QMIN(J1)
              READ(MYUNIT2,*) QMINNATOMS(J1)
              READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
           ENDDO
        ELSEIF (NSAVE.LT.OLDSAVE) THEN
           PRINT *,'NSAVE<OLDSAVE - truncating read in'
           DO J1=1,OLDSAVE
              IF (J1.LE.NSAVE) THEN
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) QMIN(J1)
                 READ(MYUNIT2,*) NPCALL_QMIN(J1)
                 READ(MYUNIT2,*) QMINNATOMS(J1)
                 READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
              ELSE 
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) QMINNATOMS(J1)
                 READ(MYUNIT2,*) DUMMYARRAY(1:3*QMINNATOMS(J1))
              ENDIF
           ENDDO
        ELSEIF (NSAVE.GT.OLDSAVE) THEN
           PRINT *,'NSAVE>OLDSAVE - padding QMIN and QMINP with zeros'
           DO J1=1,NSAVE
              IF (J1.LE.OLDSAVE) THEN
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) QMIN(J1)
                 READ(MYUNIT2,*) NPCALL_QMIN(J1)
                 READ(MYUNIT2,*) QMINNATOMS(J1)
                 READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
              ELSE 
                 QMIN(J1)=0.0D0
                 NPCALL_QMIN(J1)=0
                 QMINNATOMS(J1)=NATOMS
                 QMINP(J1,1:3*QMINNATOMS(J1))=0.0D0
              ENDIF
           ENDDO
        ENDIF
        IF (GCBHT) THEN
           READ(MYUNIT2, *) DUMMYS
           DO J1=1,GCNATOMS
              READ(MYUNIT2, *) DUMMYS
              READ(MYUNIT2, *) QENERGIES(J1),QPE(J1)
              READ(MYUNIT2, *) QCOORDINATES(J1,1:3*J1)
           ENDDO
        ENDIF

! read in EBEST and BESTCOORDS
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) JBEST(JP), EBEST(JP)
   READ(MYUNIT2,*) NDUMMY
   READ(MYUNIT2,*) BESTCOORDS(1:3*NDUMMY,JP)
! bs360: the following is not tested for MPI (July 2010)
  READ(MYUNIT2,*) DUMMYS
  READ(MYUNIT2,*) NMSBSAVE
  DO J1=1,MIN(NMSBSAVE,MAXSAVE)
     READ(MYUNIT2,*) DUMMYS
     READ(MYUNIT2,*) MSBE(J1)
     READ(MYUNIT2,*) AVOIDNATOMS(J1)
     READ(MYUNIT2,*) MSBCOORDS(1:3*AVOIDNATOMS(J1),J1)
  ENDDO

  CLOSE(MYUNIT2)


!   csw34> Add in restoring of the lowest interaction energies if A9INTE is specified
!          INTERESTORE is .TRUE. if the RESTORE keyword has two arguements. The second is
!          the name of the interaction energy dump file!
IF (NPAR.EQ.1) THEN
IF (A9INTET.AND.INTERESTORE) THEN
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(INTEDUMPFILE)),STATUS='OLD')
   READ(LUNIT,*) DUMMYS
   READ(LUNIT,'(I8)') OLDSAVE
   READ(LUNIT,*) DUMMYS
   READ(LUNIT,'(I8)') NDONE
   READ(LUNIT,*) DUMMYS
!       csw34> loop as above to allow SAVE to vary between runs
   IF (NSAVEINTE.EQ.OLDSAVE) THEN
           DO J1=1,NSAVEINTE
                   READ(LUNIT,*) DUMMYS
                   READ(LUNIT,*) INTEQMIN(J1)
                   READ(LUNIT,*) INTEQMINP(J1,1:3*NATOMS)
           ENDDO
!       csw34> NSAVEINTE<OLDSAVE
   ELSE IF (NSAVEINTE.LT.OLDSAVE) THEN
           PRINT *,'NSAVEINTE<OLDSAVE - truncating read in'
           DO J1=1,OLDSAVE
                   IF (J1.LE.NSAVEINTE) THEN
                           READ(LUNIT,*) DUMMYS
                           READ(LUNIT,*) INTEQMIN(J1)
                           READ(LUNIT,*) INTEQMINP(J1,1:3*NATOMS)
                   ELSE
                           READ(LUNIT,*) DUMMYS
                           READ(LUNIT,*) DUMMYS
!       csw34> This is the key bit, you've got to skip the correct number of
!       lines and that includes those for the coordinates. The new array
!       DUMMYARRAY allows this in a very simple way.
                           READ(LUNIT,*) DUMMYARRAY(1:3*NATOMS)
                  ENDIF
           ENDDO
!       csw34> NSAVEINTE>OLDSAVE
   ELSE IF (NSAVEINTE.GT.OLDSAVE) THEN
           PRINT *,'NSAVEINTE>OLDSAVE - padding INTEQMIN and INTEQMINP with zeros'
           DO J1=1,NSAVEINTE
                   IF (J1.LE.OLDSAVE) THEN
                           READ(LUNIT,*) DUMMYS
                           READ(LUNIT,*) INTEQMIN(J1)
                           READ(LUNIT,*) INTEQMINP(J1,1:3*NATOMS)
                   ELSE
                           INTEQMIN(J1)=0.0D0
                           INTEQMINP(J1,1:3*NATOMS)=0.0D0
                   ENDIF
           ENDDO
   ENDIF
CLOSE(UNIT=LUNIT)
ENDIF
ENDIF

RETURN
END SUBROUTINE RESTORESTATE
