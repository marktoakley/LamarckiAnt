        PROGRAM XYZ2MOVIESEG
        IMPLICIT NONE

        INTEGER, PARAMETER :: MAXATOMS=10000
        INTEGER, PARAMETER :: MAXRES=1700

        REAL TGCORD(MAXRES,3,3)

        CHARACTER(LEN=2) :: SDUMMY
        DOUBLE PRECISION :: DUMMY, X, Y, Z
        INTEGER NDUMMY, NFRAMES, NATOMS, J1, I1, NRES
        INTEGER I_RES, SEQ(MAXRES), NOGLY, GLY

!     This program converts xyz files into a movieseg file for the AMH. 
!     A few pieces of information are hardwired and will need changing. 
!     The output is printed to STDOUT.
!     1.  The sequence is read from Mike's amh/proteins directory.
!         The name may need to be changed to the 5 character name of the 
!         protein which can be found in the first line of the pro.list 
!         file for this system.
!     2.  The xyz file is called target.xyz.
!     1.  path to amh/proteins

!        OPEN(30,FILE='/home/mp466/amh/proteins/1tjba',STATUS='OLD')
        OPEN(30,FILE='/home/mprentis/amh/proteins/1tjba',STATUS='OLD')
           READ(30,*)
           READ(30,*)NRES
           WRITE(6,*)NRES
            IF (NRES.GT.500) THEN
               WRITE(6,*)'FAILURE NRES GR THAN 500 COUNTATOMS'
               STOP
            ENDIF
           READ (30,25)(SEQ(I_RES),I_RES=1,NRES)
!           WRITE(6,25)(SEQ(I_RES),I_RES=1,NRES)
25         FORMAT(25(I2,1X))
        CLOSE(30)

        NOGLY = 0
        GLY = 0

        DO I_RES=1,NRES
           IF (SEQ(I_RES).NE.8) NOGLY = NOGLY +1
           IF (SEQ(I_RES).EQ.8) GLY = GLY +1
        ENDDO

        NATOMS = NOGLY*3 + GLY*2

        OPEN(UNIT=3,FILE='target.xyx',STATUS='OLD')

        WRITE(6,334)NRES,3,1,1
334     FORMAT(4(I8,1X),' NMRES NMCRD NUMPRO NMSNAP')
        WRITE(6,683)1,1,1,1.0,1
683     FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

        DO
           READ(3,*)NDUMMY
           READ(3,*) 
         DO J1=1,NRES
             IF (SEQ(J1).NE.8) THEN

               READ(3,*)SDUMMY,X,Y,Z 
               TGCORD(J1,1,1)=REAL(X)
               TGCORD(J1,2,1)=REAL(Y)
               TGCORD(J1,3,1)=REAL(Z)

               READ(3,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,2)=REAL(X)
               TGCORD(J1,2,2)=REAL(Y)
               TGCORD(J1,3,2)=REAL(Z)

               READ(3,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,3)=REAL(X) 
               TGCORD(J1,2,3)=REAL(Y)
               TGCORD(J1,3,3)=REAL(Z)
        
               WRITE(6,332)(TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,2),I1=1,3),(TGCORD(J1,I1,3),I1=1,3)
332            FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1X))
          ELSE

               READ(3,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,1)=REAL(X)
               TGCORD(J1,2,1)=REAL(Y)
               TGCORD(J1,3,1)=REAL(Z)

               READ(3,*)SDUMMY,X,Y,Z
               TGCORD(J1,1,3)=REAL(X)
               TGCORD(J1,2,3)=REAL(Y)
               TGCORD(J1,3,3)=REAL(Z)
               WRITE(6,332)(TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,1),I1=1,3),(TGCORD(J1,I1,3),I1=1,3)

         ENDIF
        ENDDO
       ENDDO

      END PROGRAM XYZ2MOVIESEG
