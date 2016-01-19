       PROGRAM XYZ2PDB
       IMPLICIT NONE

       INTEGER, PARAMETER :: MAXATOMS=10000
       CHARACTER(LEN=30) :: S1(MAXATOMS)
       CHARACTER(LEN=10) :: S2(MAXATOMS)
       CHARACTER(LEN=2) :: SDUMMY
       DOUBLE PRECISION :: R1(MAXATOMS), R2(MAXATOMS), DUMMY, X, Y, Z
       INTEGER NDUMMY, NFRAMES, J1, NATOMS


!     This program converts a single xyz file into a pdb file
!     A few pieces of information 
!     are hardwired and will need changing. The output is printed to 
!     out.pdb.
!     1.  The sequence is read from Mike's amh/proteins directory.
!         The name may need to be changed to the 5 character name of the 
!         protein which can be found in the first line of the pro.list 
!         file for this system.
!     2.  The path to the target xyz is set to the local directory.


       OPEN (UNIT=1,FILE='1r69z.pdb',STATUS='OLD')
       NATOMS=0
       DO 
           NATOMS=NATOMS+1
           READ(1,'(A30,3F8.3,2F6.2,A10)',END=10) S1(NATOMS),DUMMY,DUMMY,DUMMY,R1(NATOMS),R2(NATOMS),S2(NATOMS)
       ENDDO
       10 NATOMS=NATOMS-1
       CLOSE(1)
       PRINT '(I6,A)',NATOMS,' coordinates read from dummy pdb file'

       OPEN (UNIT=2,FILE='out.pdb',STATUS='UNKNOWN')
       OPEN (UNIT=3,FILE='target.xyz',STATUS='OLD')

        NFRAMES=0
       DO
         READ(3,*,END=20) NDUMMY
         READ(3,*) 
         WRITE(2,'(A6)') 'HEADER'
       DO J1=1,NATOMS
         READ(3,*) SDUMMY,X,Y,Z 
         WRITE(2,'(A30,3F8.3,2F6.2,A10)') S1(J1),X,Y,Z,R1(J1),R2(J1),S2(J1)
       ENDDO
       WRITE(2,'(A3)') 'END'
       NFRAMES=NFRAMES+1
       ENDDO
       20 PRINT '(I6,A)',NFRAMES,' frames converted to pdb format and written to out.pdb'

       END PROGRAM XYZ2PDB
