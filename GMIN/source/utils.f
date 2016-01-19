!   GMIN: A program for finding global minima
!   Copyright (C) 1999- David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!      SUBROUTINE TOKENISE(STRING, TOKENS)
!         IMPLICIT NONE 
!         CHARACTER(LEN=*), INTENT(IN)  :: STRING
!         
!      END SUBROUTINE TOKENISE
      
      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
         INTEGER                       :: GETUNIT
         
         FILE_UNIT = GETUNIT()
         IF (APPEND) THEN
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN', POSITION='APPEND')
         ELSE
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN')
         END IF
      
      END SUBROUTINE FILE_OPEN

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
      !
      ! start checking for available units > 103, to avoid system default units
      ! 100, 101 and 102 are stdin, stdout and stderr respectively.
      ! 
         INTEGER :: UNITNUM

         INUSE=.TRUE.
         UNITNUM=103

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT

      INTEGER FUNCTION FILE_LENGTH(FILE_NAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         INTEGER                       :: FILE_UNIT
         CHARACTER(LEN=80)             :: DUMMY_LINE
         INTEGER                       :: IO_STATUS

         CALL FILE_OPEN(FILE_NAME, FILE_UNIT, .FALSE.)
         FILE_LENGTH = 0
         DO
            READ(FILE_UNIT, *, IOSTAT=IO_STATUS) DUMMY_LINE
            IF (IO_STATUS /= 0) EXIT
            FILE_LENGTH = FILE_LENGTH + 1
         ENDDO

         CLOSE(FILE_UNIT)
      END FUNCTION FILE_LENGTH 

      SUBROUTINE VECNORM(VECTOR,NOPT)
      !
      ! normalises VECTOR (adapted from OPTIM)
      !
         IMPLICIT NONE
         INTEGER :: J1 
         DOUBLE PRECISION :: DUMMY 
         INTEGER, INTENT(IN) :: NOPT
         DOUBLE PRECISION, INTENT(INOUT) :: VECTOR(NOPT)

         DUMMY=0.0D0
         DO J1=1,NOPT
            DUMMY=DUMMY+VECTOR(J1)**2
         ENDDO

         IF (DUMMY.GT.0.0D0) THEN
            DUMMY=1.0D0/DSQRT(DUMMY)
            DO J1=1,NOPT
               VECTOR(J1)=VECTOR(J1)*DUMMY
            ENDDO
         ELSE
            PRINT *,'WARNING: zero size vector passed to VECNORM - STOP'
            STOP
         ENDIF

      END SUBROUTINE VECNORM

