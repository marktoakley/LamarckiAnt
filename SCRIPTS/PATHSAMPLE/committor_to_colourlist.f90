! This little utility program takes a file commit.data from PATHSAMPLE and 
! picks out the minima with committor probabilities within the range specified 
! via the command line arguments.  Output is a file which contains a list of 
! indices of minima, one per line, with committors in that range. This file is 
! appropriate for use in disconnectionDPS to colour a tree etc.
! The commit.data file must be in the working directory where this program is run.

! The three command line arguments to this program are:
! 1. lower bound on the committor probability range (double precision variable)
! 2. upper bound on the committor probability range (double precision variable)
! 3. output file name (80 characters or less)

! JMC March 2011

PROGRAM commit_to_colourlist
IMPLICIT NONE

INTEGER :: ierr, i
CHARACTER(LEN=80) :: FILENAME, buffer
DOUBLE PRECISION :: commit, low, high

CALL GETARG(1,buffer)
READ(buffer,*) low
CALL GETARG(2,buffer)
READ(buffer,*) high
CALL GETARG(3,FILENAME)
OPEN(UNIT=10,FILE=FILENAME,STATUS='UNKNOWN')

PRINT *,'List of minima with committor probabilities between ',low,' and ',high
PRINT *,'will be written to file ',FILENAME

OPEN(UNIT=11,FILE='commit.data',STATUS='OLD')

DO i=1,huge(1)
   READ(11,*,IOSTAT=ierr) commit
   IF (ierr /= 0) EXIT
   IF (commit < high .AND. commit >= low) WRITE(10,*) i
END DO

CLOSE(10)
CLOSE(11)
END PROGRAM commit_to_colourlist
