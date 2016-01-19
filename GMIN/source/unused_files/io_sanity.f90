! GMIN I/O FILE OPENING MODULE
!
! CONTENTS
! --------
!
! FUNCTIONS:
! 
! - OPENINPUT(FILENAME) <UNDER CONSTRUCTION>
! - OPENDATA(FILENAME) <UNDER CONSTRUCTION>
! - READTOBUFFER(FILEUNIT) - returns a BUFFER(100,NLINES), containing each line
!   of FILEUNIT as a 100 character wide array
! - FILELENGTH(FILEUNIT) - returns the number of lines in FILEUNIT
! - FILEEXISTS(FILEUNIT) - returns .TRUE. if FILEUNIT exists
!
! SUBROUTINES:
!

FUNCTION OPENINPUT(FILENAME)
! Open the input file and perform sanity checks, returning molecular information
! for suitable files
    IMPLICIT NONE

    CHARACTER, DIMENSION(:), INTENT(IN) :: FILENAME
    INTEGER :: FILEUNIT
    CHARACTER, DIMENSION(:, :) :: BUFFER
    CHARACTER, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: MOLINFO    
   
    
 
    RETURN MOLINFO
END FUNCTION OPENINPUT

FUNCTION OPENDATA(FILENAME)
! Open the data file and perform sanity checks, setting the appropriate flags
! and settings
    IMPLICIT NONE

    CHARACTER, DIMENSION(:), INTENT(IN) :: FILENAME
    INTEGER :: FILEUNIT

END FUNCTION OPENDATA

FUNCTION READTOBUFFER(FILEUNIT)
! Read the file into a buffer and return the buffered file
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: FILEUNIT
    CHARACTER, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: BUFFER
    INTEGER :: FILELENGTH, IOSTATUS, LINENUMBER

    ! Find the length of the file and allocate the buffer accordingly
    FILELENGTH = FILELENGTH(FILEUNIT)
    ALLOCATE(BUFFER(100, FILELENGTH))

    ! Initialise IOSTATUS and LINENUMBER to 0 and then loop through the file
    IOSTATUS = 0
    LINENUMBER = 1
    DO WHILE (IOSTATUS .NE. -1)
        READ(UNIT = FILEUNIT, FMT = '(A100)', IOSTAT = IOSTATUS) &
        & BUFFER(:, LINENUMBER)
        LINENUMBER = LINENUMBER + 1
    END DO

    RETURN BUFFER
END FUNCTION READTOBUFFER

FUNCTION FILELENGTH(FILEUNIT)
! Returns the length of a file as an integer
    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: FILEUNIT
    INTEGER, INTENT(OUT)  :: LENGTH
    INTEGER :: IOSTATUS

    LENGTH = 0
    IOSTATUS = 0

    ! If it reads past the end of the file, IOSTATUS=-1
    DO WHILE (IOSTATUS .NE. -1)
        READ(UNIT = FILEUNIT, IOSTAT = IOSTATUS, FMT = *)
        LENGTH = LENGTH + 1
    END DO

    RETURN LENGTH
END FUNCTION FILELENGTH

FUNCTION FILEEXISTS(FILEUNIT)
! Check whether the file exists
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: FILEUNIT
    LOGICAL, INTENT(OUT) :: EXISTS

    EXISTS = .FALSE.
    INQUIRE(UNIT = FILEUNIT, EXIST = EXISTS)

    RETURN EXISTS 
END FUNCTION FILEEXISTS

! GMIN I/O ROUTINES
!
! CONTENTS
! --------
!
! FUNCTIONS:
! 
! - READPDB(BUFFER) - returns MOLINFO and COORDS from BUFFER <UNDER CONSTRUCTION>
! - WRITEPDB(NATOMS,MOLINFO,COORDS,FILENAME) <UNDER CONSTRUCTION>
! - READRST(BUFFER) - returns COORDS 1D array from BUFFER <UNDER CONSTRUCTION>
! - WRITERST(NATOMS,COORDS,FILENAME) - writes coordinates into AMBER FILENAME.rst file
! - READCRD(BUFFER) <UNDER CONSTRUCTION>
! - WRITECRD(NATOMS,MOLINFO,COORDS,FILENAME) <UNDER CONSTRUCTION>
! - READXYZ(BUFFER) <UNDER CONSTRUCTION>
!
! SUBROUTINES:
!

FUNCTION WRITEPDB(NATOMS,MOLINFO,COORDS,FILENAME)
! Write COORDS to FILENAME.pdb in PDB format
! MOLINFO(NATOMS) is a derived type array 
! Each element (i.e. MOLINFO(1)) contains info about that atom:
! MOLINFO(1)%index = atom index in coords array
! MOLINFO(1)%name = atom name
! MOLINFO(1)%resid = residue number
! MOLINFO(1)%resname = residue name
! MOLINFO(1)%chainid = chain identifier
! MOLINFO(1)%segname = segment name
! MOLINFO(1)%element = element name
! MOLINFO(1)%occupancy = site occupancy
! MOLINFO(1)%bfactor = crystallographic temperature/b factor
! MOLINFO(1)%charge = atomic charge
! MOLINFO(1)%weight = CHARMM weighting
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
