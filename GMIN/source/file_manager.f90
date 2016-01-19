MODULE FILE_MANAGER

CONTAINS

   FUNCTION CHECK_FILE(FILE_UNIT, FOR_READ, FOR_WRITE)
   IMPLICIT NONE
   ! Function
   LOGICAL             :: CHECK_FILE
   ! Arguments
   INTEGER, INTENT(IN) :: FILE_UNIT
   LOGICAL, INTENT(IN) :: FOR_READ, FOR_WRITE
   ! Local variables
   CHARACTER (LEN=7)   :: FILE_READABLE, FILE_WRITABLE
   LOGICAL             :: FILE_EXISTS, FILE_OPEN
   
   ! Initialise to .FALSE. (i.e. file isn't good for use)
   CHECK_FILE = .FALSE.
   
   INQUIRE(FILE_UNIT, EXIST=FILE_EXISTS, OPENED=FILE_OPEN, READ=FILE_READABLE, WRITE=FILE_WRITABLE)
   ! Did the file open correctly and connect to the unit?
   IF (.NOT. (FILE_EXISTS .AND. FILE_OPEN)) THEN
      RETURN
   END IF
   ! Can we read, if read access has been requested?
   IF (FOR_READ) THEN
      IF (FILE_READABLE .NE. 'YES') THEN
         RETURN
      END IF
   END IF
   ! Can we write, if write access has been requested?
   IF (FOR_WRITE) THEN
      IF (FILE_WRITABLE .NE. 'YES') THEN
         RETURN
      END IF
   END IF
   ! If all the other checks succeed, set CHECK_FILE to .TRUE., the file is ok
   CHECK_FILE = .TRUE.
   
   RETURN
   END FUNCTION CHECK_FILE

END MODULE FILE_MANAGER
