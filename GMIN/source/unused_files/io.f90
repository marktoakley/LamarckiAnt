! GMIN I/O ROUTINES
!
! CONTENTS
! --------
!
! FUNCTIONS:
! 
! - FILEEXISTS(FILENAME) - returns .TRUE. if FILENAME exists
! - FILELENGTH(FILENAME) - returns the number of lines in FILENAME
! - WRITEPDB(NATOMS,MOLINFO,COORDS,FILENAME) <UNDER CONSTRUCTION>
! - READRST(NATOMS,FILENAME) - returns coordinates in 1D array from AMBER FILENAME.rst file
! - WRITERST(NATOMS,COORDS,FILENAME) - writes coordinates into AMBER FILENAME.rst file
! - WRITECRD(NATOMS,MOLINFO,COORDS,FILENAME) <UNDER CONSTRUCTION>
!
! SUBROUTINES:
!

FUNCTION FILEEXISTS(FILENAME)
! Takes a file name and returns .TRUE. if the file exists, and .FALSE. if it does not
IMPLICIT NONE

CHARACTER(LEN=40), INTENT(IN) :: FILENAME
LOGICAL, INTENT(OUT) :: EXISTS

EXISTS=.FALSE.
INQUIRE(FILE=TRIM(ADJUSTL(FILENAME),EXIST=EXISTS)

RETURN EXISTS 
END FUNCTION FILEEXISTS



FUNCTION FILELENGTH(FILENAME)
! Returns the length of a file as an integer
IMPLICIT NONE

CHARACTER(LEN=40), INTENT(IN) :: FILENAME
INTEGER, INTENT(OUT)  :: LENGTH
INTEGER :: IOSTATUS

LENGTH=0
! Open the file
OPEN(UNIT=20,FILE=TRIM(ADJUSTL(FILENAME),STATUS='OLD')
DO
   READ(UNIT=20,IOSTAT=IOSTATUS,FMT='(A6)') LINE
! If it reads past the end of the file, IOSTATUS=-1
   IF (IOSTATUS<0) THEN
      CLOSE(20)
      EXIT
   ELSE 
      LENGTH=LENGTH+1
   ENDIF
END DO

RETURN LENGTH
END FUNCTION FILELENGTH



FUNCTION WRITEPDB(NATOMS,MOLINFO,COORDS,FILENAME)
! Write COORDS to FILENAME.pdb in PDB format
! MOLINFO(NATOMS) is a derived type array 
! Each element (i.e. MOLINFO(1)) contains info about that atom:
! MOLINFO(1)%name = atom name
! MOLINFO(1)%resid = residue number
! MOLINFO(1)%resname = residue name
! MOLINFO(1)%chainid = chain identifier
! MOLINFO(1)%segname = segment name
! MOLINFO(1)%element = element name
! MOLINFO(1)%occupancy = site occupancy
! MOLINFO(1)%bfactor = crystalographic temperature/b factor
! MOLINFO(1)%charge = atomic charge
! MOLINFO(1)%weight = CHARMM weighing
IMPLICIT NONE
CHARACTER(LEN=40), INTENT(IN) :: FILENAME
INTEGER, INTENT(IN) :: NATOMS 
INTEGER :: I
TYPE(atom), INTENT(IN) :: MOLINFO(NATOMS)

OPEN(UNIT=20,FILE=TRIM(ADJUSTR(FILENAME))//TRIM(ADJUSTL('.pdb')),STATUS='UNKNOWN')
DO I=1,NATOMS
   WRITE(20,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)')
     &        'ATOM',I,MOLINFO(I)%name,' ',MOLINFO(I)%resname,MOLINFO(I)%chainid
     &        ,MOLINFO(I)%resid,' ',COORDS(3*I-2:3*I),MOLINFO(I)%occupancy
     &        ,MOLINFO(I)%bfactor,MOLINFO(I)%element,MOLINFO(I)%charge
ENDDO
CLOSE(20)
RETURN
END FUNCTION WRITEPDB



FUNCTION READRST(FILENAME)
! Reads COORDS from FILENAME.rst (AMBER restart file)
IMPLICIT NONE

CHARACTER(LEN=40), INTENT(IN) :: FILENAME
CHARACTER(LEN=1) :: DUMMY
INTEGER :: NATOMS, I
DOUBLE PRECISION, INTENT(OUT) :: COORDS(3*NATOMS)

NATOMS=0
! Open the file
OPEN(UNIT=20,FILE=TRIM(ADJUSTR(FILENAME))//TRIM(ADJUSTL('.rst')),STATUS='UNKNOWN')
! .rst files have a title line - read it into DUMMY
READ(20,'(20A4)') DUMMY
! Read the number of atoms from the .rst file
READ(20,'(I5)') NATOMS
! Read coordinates
READ(20,'(6F12.7)') (COORDS(I),I=1,3*NATOMS)
! Close the file
CLOSE(20)

! Return the coordinates in place of the function call
RETURN COORDS
END FUNCTION READRST



FUNCTION WRITERST(NATOMS,COORDS,FILENAME)
! Writes COORDS to FILENAME.rst in AMBER restart format
IMPLICIT NONE

CHARACTER(LEN=40), INTENT(IN) :: FILENAME
INTEGER, INTENT(IN) :: NATOMS, I
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)

! Open the file
OPEN(UNIT=20,FILE=TRIM(ADJUSTR(FILENAME))//TRIM(ADJUSTL('.rst')),STATUS='UNKNOWN')
! Write the .rst title line. We make it the filename!
WRITE(20,'(20A4)') FILENAME
! Write the number of atoms
WRITE(20,'(I5)') NATOMS
! Write the coordinates
WRITE(20,'(6F12.7)') (COORDS(I),I=1,3*NATOMS)
! Close the file
CLOSE(20)

RETURN
END FUNCTION WRITERST



FUNCTION WRITECRD(NATOMS,MOLINFO,COORDS,FILENAME)
! Write COORDS to FILENAME.crd in CHARMM card format
! MOLINFO(NATOMS) is a derived type array 
! Each element (i.e. MOLINFO(1)) contains info about that atom:
! MOLINFO(1)%name = atom name
! MOLINFO(1)%resid = residue number
! MOLINFO(1)%resname = residue name
! MOLINFO(1)%chainid = chain identifier
! MOLINFO(1)%segname = segment name
! MOLINFO(1)%element = element name
! MOLINFO(1)%occupancy = site occupancy
! MOLINFO(1)%bfactor = crystalographic temperature/b factor
! MOLINFO(1)%charge = atomic charge
! MOLINFO(1)%weight = CHARMM weighing
IMPLICIT NONE

CHARACTER(LEN=40), INTENT(IN) :: FILENAME
INTEGER, INTENT(IN) :: NATOMS, I
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)

OPEN(UNIT=20,FILE=TRIM(ADJUSTR(FILENAME))//TRIM(ADJUSTL('.rst')),STATUS='UNKNOWN')
! CHARMM .crd files start with the number of atoms
WRITE(20,'(I5)') NATOMS
WRITE(20,'(2I5,1X,A4,1X,A4,3F17.12,2X,A4,X,A2)')
! FROM CHARMM: (info needs to be in MOLINFO)
     &        I,IRES,RES(IRES),TYPE(I),X(I),Y(I),Z(I),SEGID(ISEG),RID
CLOSE(20)
RETURN
END FUNCTION WRITECRD
