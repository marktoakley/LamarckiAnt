      PROGRAM ADD_GLY
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXATOMS=10000
      INTEGER, PARAMETER :: MAXRES=1700

      REAL TGCORD(MAXRES,3,3)

      CHARACTER(LEN=2) :: SDUMMY
      DOUBLE PRECISION :: DUMMY, X, Y, Z
      INTEGER NDUMMY, NFRAMES, NATOMS, J1, I1, NRES
      INTEGER I_RES, SEQ(MAXRES), NOGLY, GLY

!     This program adds missing GLY coordinates from files such as int*xyz
!     so that they can be read back into OPTIM.
!     A few pieces of information are hardwired and will need changing. 
!     The output is printed to STDOUT.
!     1.  The sequence is read from the directory ./proteins/<name> 
!         The name should be changed to the 5 character name of the 
!         protein which can be found in the first line of the pro.list 
!         file for this system.
!     2.  The xyz file is read from standard input

      OPEN(30,FILE='./proteins/1uamd',STATUS='OLD')
      READ(30,*)
      READ(30,*)NRES
      WRITE(6,*)NRES
      IF (NRES.GT.500) THEN
         WRITE(6,*)'FAILURE NRES GR THAN 500 COUNTATOMS'
         STOP
      ENDIF
      READ (30,25)(SEQ(I_RES),I_RES=1,NRES)
25    FORMAT(25(I2,1X))
      CLOSE(30)

      NOGLY = 0
      GLY = 0

      DO I_RES=1,NRES
         IF (SEQ(I_RES).NE.8) NOGLY = NOGLY +1
         IF (SEQ(I_RES).EQ.8) GLY = GLY +1
      ENDDO

      NATOMS = NOGLY*3 + GLY*2

      WRITE(6,334)NRES,3,1,1
334   FORMAT(4(I8,1X),' NMRES NMCRD NUMPRO NMSNAP')
      WRITE(6,683)1,1,1,1.0,1
683   FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

      DO
         READ(*,*,END=666) NDUMMY
         READ(*,*) 
         DO J1=1,NRES
            IF (SEQ(J1).NE.8) THEN
               READ(*,*)SDUMMY,X,Y,Z 
               TGCORD(J1,1,1)=REAL(X)
               TGCORD(J1,2,1)=REAL(Y)
               TGCORD(J1,3,1)=REAL(Z)

               READ(*,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,2)=REAL(X)
               TGCORD(J1,2,2)=REAL(Y)
               TGCORD(J1,3,2)=REAL(Z)

               READ(*,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,3)=REAL(X) 
               TGCORD(J1,2,3)=REAL(Y)
               TGCORD(J1,3,3)=REAL(Z)
        
!              WRITE(6,332) (TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,2),I1=1,3),(TGCORD(J1,I1,3),I1=1,3)
332            FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1X))
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,1),I1=1,3)
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,2),I1=1,3)
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,3),I1=1,3)
            ELSE

               READ(*,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,1)=REAL(X)
               TGCORD(J1,2,1)=REAL(Y)
               TGCORD(J1,3,1)=REAL(Z)

               READ(*,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,3)=REAL(X)
               TGCORD(J1,2,3)=REAL(Y)
               TGCORD(J1,3,3)=REAL(Z)
!              WRITE(6,332) (TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,3),I1=1,3)
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,1),I1=1,3)
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,1),I1=1,3)
               WRITE(*,'(3F20.10)') (TGCORD(J1,I1,3),I1=1,3)
            ENDIF
         ENDDO
      ENDDO

666   CONTINUE

      END PROGRAM ADD_GLY
