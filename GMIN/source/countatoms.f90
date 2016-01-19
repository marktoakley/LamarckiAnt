! GPL License Info {{{
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
!}}}

MODULE NOA
  USE MODAMBER9 , ONLY : INPCRD
  USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_SETUP
  IMPLICIT NONE
  SAVE
  INTEGER :: NUMBER_OF_ATOMS
  
CONTAINS
  
  SUBROUTINE COUNTATOMS(MYUNIT, NPAR, GCBHT, GCMU, GCNATOMS)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NPAR
    INTEGER :: EOF,NRES,SEQ(500),I_RES,NOGLY,GLY, MYUNIT,J
    LOGICAL :: YESNO, YESNOA, YESNOC, YESNOAMH, YESNOA9
    CHARACTER(LEN=5) TARFL
    CHARACTER(LEN=10)  CHECK
    CHARACTER(LEN=80) MYLINE,INPCRD1
    LOGICAL GCBHT, END
    DOUBLE PRECISION GCMU
    INTEGER GCNATOMS, LUNIT, GETUNIT
    CHARACTER(LEN=16) WORD
    
    YESNO=.FALSE.
    YESNOA=.FALSE.
    YESNOAMH=.FALSE.
    YESNOA9=.FALSE.
    !
    ! If the current working directory contains more than one of 
    ! these files then the precedence is coords, then input.crd, 
    ! then coords.amber 
    ! OPTIM does this a bit better by calling getparams first to 
    ! see if we are actually doing AMBER or CHARMM. 
    !
    INQUIRE(FILE='pro.list',EXIST=YESNOAMH)
    INQUIRE(FILE='coords',EXIST=YESNO)
    INQUIRE(FILE='coords.amber',EXIST=YESNOA)
    INQUIRE(FILE='input.crd',EXIST=YESNOC)
    INQUIRE(FILE='coords.inpcrd',EXIST=YESNOA9)
    
    NUMBER_OF_ATOMS=0
    
    IF (YESNO) THEN
       OPEN(UNIT=7,FILE='coords',STATUS='OLD')
       DO
          READ(7,*,IOSTAT=EOF)
          IF (EOF==0) THEN
             NUMBER_OF_ATOMS = NUMBER_OF_ATOMS + 1
          ELSE
             EXIT
          ENDIF
       ENDDO
       !js850> if a parallel run, then file coords is NPAR*NATOMS long.
       !       NPAR must be known at the time of this call
       IF ( NPAR .GT. 1 ) THEN 
          IF (MOD(NUMBER_OF_ATOMS,NPAR).NE.0) THEN
             WRITE(MYUNIT,'(A,I8,A,I8)') &
                  'Number of atoms in coords file=', &
                  NUMBER_OF_ATOMS, &
                  ' is not divisible by number of runs=',NPAR
             STOP
          ENDIF
          NUMBER_OF_ATOMS=NUMBER_OF_ATOMS/NPAR
       ENDIF
    ELSEIF (YESNOAMH) THEN
       open(unit=30,file='pro.list',status='old',form='formatted')
       read (30,1000)tarfl
1000   format(a5)
       close(30)
       
       open(30,file='proteins/'//tarfl,status='old')
       read(30,*)
       read(30,*)nres
       if (nres.gt.500) then
          write(6,*) 'failure nres gr than 500 countatoms'
          stop
       endif
       read (30,25)(seq(i_res),i_res=1,nres)
       !            write(6,25)(seq(i_res),i_res=1,nres)
25     format(25(i2,1x))
       close(30)
       
       NOGLY = 0
       GLY = 0
       
       do i_res=1,nres
          if (seq(i_res).ne.8) NOGLY = NOGLY +1
          if (seq(i_res).eq.8) GLY = GLY +1
       enddo
       
       Number_of_Atoms = NOGLY*3 + GLY*2
    ELSE IF (YESNOA9) THEN
       ! OPEN(UNIT=7,FILE='coords.gayberne',STATUS='OLD')
       ! PRINT '(A)','reading coordinates from file coords.gayberne'
       ! inpcrd1='coords.inpcrd'
       ! inpcrd1=trim(adjustl(inpcrd1))
       ! call amberinterface(Number_of_Atoms,1,inpcrd1,MYUNIT)
       CALL AMBER12_SETUP(NUMBER_OF_ATOMS)
       IF (NUMBER_OF_ATOMS .EQ. 0) THEN
          NUMBER_OF_ATOMS = COUNT_ATOMS_FROM_INPCRD('coords.inpcrd')
       END IF
       
    ELSEIF (YESNOC) THEN
       OPEN(UNIT=7,FILE='input.crd',STATUS='OLD')
       do
          read(7,*) myline
          if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
             cycle
          else
             read(myline,*) Number_of_Atoms
             exit
          endif
       enddo
       
       ! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
       ! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.
       
       CALL GETMAXAIM
       WRITE(MYUNIT,'(A,I8)') 'countatoms> Number_of_Atoms=',Number_of_Atoms
    ELSEIF (YESNOA) THEN
       OPEN(UNIT=7,FILE='coords.amber',STATUS='OLD')
       do
          read(7,'(A3)',iostat=eof) check
          if (eof.LT.0) then
             PRINT *,'End of file before all information specified'
             STOP
          ENDIF
          IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') THEN
             CLOSE(7)
             EXIT
          ENDIF
          Number_of_Atoms = Number_of_Atoms + 1
       enddo
    ELSE
! ------------------------
! vr274 (DMACRYS)
! now try to parse the data file in advance to check what's going on.
! DMACRYS doen's have a coordiate file, that's why we initialize the dmacrys stuff here in
! a firts pass reading the data file
! -------------------------
       if(.not.countatoms_from_datafile()) then
          PRINT '(A)','ERROR - no coords, input.crd, coords.inpcrd or coords.amber &
               &file to determine number of atoms. Also could not be guessed from &
               &data file (only works for DMACRYS)'
          STOP
       endif
    ENDIF
    
    CLOSE(7)
    
    ! ds656> Determine the number of species from 'data'
    CALL COUNT_SPECIES_FROM_DATAFILE()

    LUNIT=GETUNIT()
    OPEN(UNIT=LUNIT,FILE='data',STATUS='OLD')
    
190 CALL INPUT(END,LUNIT)
    IF (.NOT. END) THEN
       CALL READU(WORD)
    ENDIF
    IF (END .OR. WORD .EQ. 'STOP') GOTO 191
!
! Grand canonical basin-hopping.
!
    IF (WORD.EQ.'GCBH') THEN
       GCBHT=.TRUE.
       CALL READF(GCMU) ! chemical potential
       CALL READI(GCNATOMS) ! maximum number of atoms
       WRITE(MYUNIT,'(A,G20.10,A,I6)') &
            'countatoms> Grand canonical basin-hopping, chemical potential=', &
            GCMU,' maximum atoms=',GCNATOMS
       GOTO 191
    ENDIF
    
    GOTO 190
191 CONTINUE
    CLOSE(LUNIT)
    
    !print *, Number_of_Atoms, ' atoms in the system'    

    RETURN
    !
  END SUBROUTINE COUNTATOMS
  !
  SUBROUTINE COUNT_SPECIES_FROM_DATAFILE()
    !
    USE COMMONS, ONLY : NSPECIES, NSPECIES_INI
    !
    IMPLICIT NONE
    !
    LOGICAL :: END
    INTEGER :: LUNIT, J, GETUNIT
    CHARACTER(LEN=16) :: WORD
    !
    ! The following variables are for the common block
    ! containing NITEMS, which is required for GLJ keyword.
    LOGICAL :: SKIPBL,CLEAR,ECHO,CAT
    INTEGER :: ITEM,NITEMS,LOC,LINE,NCR,NERROR,LAST
    COMMON /BUFINF/ ITEM,NITEMS,LOC(80),LINE,SKIPBL, &
         CLEAR,NCR,NERROR,ECHO,LAST,CAT
    !
    LUNIT=GETUNIT()
    J=1
    !
    !ds656> Avoiding GOTO...
    OPEN(UNIT=LUNIT,FILE='data',STATUS='OLD')
    DATALOOP: DO ! Loop through lines in data file
       CALL INPUT(END,LUNIT)
       IF(END) EXIT DATALOOP
       CALL READU(WORD)
       IF(WORD.EQ.'STOP') EXIT DATALOOP
       IF(WORD.EQ.'BLJCLUSTER_NOCUT'.OR.WORD.EQ.'BGUPTAT') THEN
          J=2
          EXIT DATALOOP
       ELSEIF(WORD.EQ.'MLJ'.OR.WORD.EQ.'MGUPTA'.OR.&
            WORD.EQ.'MSC') THEN
          CALL READI(J)
          EXIT DATALOOP
       ELSEIF(WORD.EQ.'GLJ') THEN
          J=NITEMS-1
          EXIT DATALOOP
       ENDIF
    ENDDO DATALOOP
    CLOSE(LUNIT)
    !
    ! Now allocate and initialise NSPECIES.
    ALLOCATE(NSPECIES(0:J), NSPECIES_INI(0:J))
    NSPECIES(0) = J; NSPECIES(1:J) = 0
    NSPECIES_INI(:) = NSPECIES(:)
    !
    RETURN
    !
  END SUBROUTINE COUNT_SPECIES_FROM_DATAFILE
!
! ------------------------
! vr274 (DMACRYS)
! parses the datafile to check whether to run DMACRYS
! returns true if number of atoms could be determined
! -------------------------
      function countatoms_from_datafile() result(found)
         implicit none

         logical done, found
         character word*16
         ! we need these variables for common block of parser - change to module!!
         logical SKIPBL, CLEAR, ECHO, CAT
         INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, LAST, DATA_UNIT
         COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,&
                     NERROR, ECHO, LAST, CAT

         ! so far we didn't find any information on atoms
         found=.false.
         print *,"COUNTATOMS> Trying to determine atoms from data file"
!         open (5,FILE='data',STATUS='OLD')
         CALL FILE_OPEN('data', DATA_UNIT, .FALSE.)
         do
            call INPUT(done, DATA_UNIT)
            if(done) exit
            call READU(word)
            ! user wants to run dmacrys, so read in file here
            if(word.eq.'DMACRYS') then
               print *,"COUNTATOMS> Found DMACRYS, initializing"
               call DMACRYS_SETUP
               call dmacrys_get_natoms(NUMBER_OF_ATOMS)
               NUMBER_OF_ATOMS = NUMBER_OF_ATOMS + 2 ! for the strain matrix
               found=.true.
! vr274
! we cannot rely on a file for a in general defined user potential, therefore call
! routines from there to get number of atom informations.
! TODO: We should change number of atoms to degrees of freedom.
! -------------------------

            else if(word.eq.'USERPOT') then
               print *,"COUNTATOMS> Using User defined potential, initializing"
               call USERPOT_INIT
               call USERPOT_GET_NATOMS(NUMBER_OF_ATOMS)
               found=.true.
            endif

         end do
         close(DATA_UNIT)
      end function

      FUNCTION COUNT_ATOMS_FROM_INPCRD(INPCRD_FILE) RESULT(NUM_ATOMS)
         IMPLICIT NONE
! Arguments and return value
         CHARACTER(LEN=*), INTENT(IN)  :: INPCRD_FILE
         INTEGER                       :: NUM_ATOMS
! Local variables         
         INTEGER                       :: NUM_COORDS, NUM_ATOMS_FROM_FILE, &
                                        & UNIT_NUMBER, I, READ_STATUS, &
                                        & COORD_START, COORD_END
         CHARACTER(LEN=80)             :: CURRENT_LINE
         LOGICAL                       :: UNIT_IN_USE

! Check which unit numbers are available for use at the moment.
         DO I=5001, 5999
            INQUIRE(I, OPENED=UNIT_IN_USE)
            IF (.NOT. UNIT_IN_USE) THEN
               UNIT_NUMBER = I
               EXIT
            END IF
         END DO
! Open the input coordinates file and find the number of atoms reported on the 
! second line of the file.
         OPEN(UNIT = UNIT_NUMBER, FILE = INPCRD_FILE, STATUS = 'OLD')
         READ(UNIT_NUMBER, '(A80)') CURRENT_LINE
         READ(UNIT_NUMBER,*) NUM_ATOMS_FROM_FILE
! Now read the coordinates and count how many there are (should be 
! 3 * num_atoms)
         READ_STATUS = 0
         NUM_COORDS  = 0
         DO WHILE (READ_STATUS .EQ. 0)
            CURRENT_LINE = ''
            READ(UNIT_NUMBER, '(A80)', IOSTAT = READ_STATUS) CURRENT_LINE
            DO I = 0,5
               COORD_START = I * 12 + 1
               COORD_END   = I * 12 + 12
               IF ((CURRENT_LINE .NE. '') .AND. &
                  &(CURRENT_LINE(COORD_START:COORD_END) .NE. '            '))&
                  & THEN
                  NUM_COORDS = NUM_COORDS + 1
               END IF
            END DO
         END DO
! Close the input coordinates file.
         CLOSE(UNIT = UNIT_NUMBER)
! If there's a discrepancy between the number of atoms the file claims that 
! there are and the number of coordinates 
! (3 * number of atoms (-1 for periodic)), stop the program and throw an error.
         IF (NUM_ATOMS_FROM_FILE .NE. NUM_COORDS / 3 .AND. &
             NUM_ATOMS_FROM_FILE .NE. NUM_COORDS / 3 - 1) THEN
            STOP 'There is a discrepancy between the number of atoms listed &
                 &in the input file and the number of coordinates. Please &
                 &check your input.'
         END IF
! Otherwise, return the number of atoms.
         NUM_ATOMS = NUM_ATOMS_FROM_FILE
      END FUNCTION COUNT_ATOMS_FROM_INPCRD

END MODULE NOA
