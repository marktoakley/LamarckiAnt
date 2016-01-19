      PROGRAM XYZ2PDB_AMH
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXATOMS=10000
      CHARACTER(LEN=2) :: DUMMY
      CHARACTER(LEN=3) :: AANAME(20),RES_TYPE
      CHARACTER(LEN=5) :: MEMNAME 
      INTEGER :: SEQ(500),GLY,NRES_AMH, JTGRES

      DOUBLE PRECISION :: STRUCT(MAXATOMS)
      DOUBLE PRECISION :: CPRCORD_X(MAXATOMS),CPRCORD_Y(MAXATOMS),CPRCORD_Z(MAXATOMS)
      DOUBLE PRECISION :: NITCORD_X(MAXATOMS),NITCORD_Y(MAXATOMS),NITCORD_Z(MAXATOMS)

      INTEGER :: NDUMMY, NFRAMES, J1, NATOMS, J2, J3, NRES, NCOUNT, I_RES
      INTEGER :: GLYC, COUNTER, OPEN_STATUS

      DATA AANAME /'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS', &
         'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'/

!     This program converts path.xyz files into a set of pdb files
!     for generating movies for the AMH.  A few pieces of information 
!     are hardwired and will need changing. The output is printed to 
!     STDOUT 
!     1.  The sequence is read from Mike's amh/proteins directory.
!     2.  The path to the path.xyz is set to the local directory.

!     1.  path to amh/proteins

      OPEN(8,FILE='pro.list',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_STATUS)
         IF (OPEN_STATUS.NE.0) THEN
           WRITE(6,*) 'FAILURE TO OPEN PRO.LIST FILE'
           WRITE(6,*) 'NEEDS TO BE IN LOCAL DIRECTORY'
           WRITE(6,*) 'ERROR NUMBER ',OPEN_STATUS
           STOP
         ENDIF
         READ(8,2002)MEMNAME
2002     FORMAT (A5)
      CLOSE(8)

      OPEN (UNIT=1,FILE='./proteins/'//memname,STATUS='OLD')
        READ(1,*)
        READ(1,*)NRES
           IF (NRES.GT.500) THEN
              WRITE(6,*) 'FAILURE NRES_AMH GR THAN 500 CONNECTODATA'
              STOP
           ENDIF
        READ (1,25)(SEQ(I_RES),I_RES=1,NRES)
!        WRITE (6,25)(SEQ(I_RES),I_RES=1,NRES)
25      FORMAT(25(I2,1X))
      CLOSE(1)

       CALL SYSTEM('cat path.xyz | grep  "Ener" | wc -l > temptemp') 

       OPEN(UNIT=1,FILE='temptemp',STATUS='OLD')
         READ(1,*)NFRAMES
       CLOSE(1)

       CALL SYSTEM('rm temptemp') 

          DO J3 = 1,NRES
             CPRCORD_X(J3)=0.D0
             CPRCORD_Y(J3)=0.D0
             CPRCORD_Z(J3)=0.D0
             NITCORD_X(J3)=0.D0
             NITCORD_Y(J3)=0.D0
             NITCORD_Z(J3)=0.D0
          ENDDO
! 2.  path to the path.xys is set to the local directory

      OPEN (UNIT=99,FILE='path.xyz',STATUS='OLD')
      DO J3=1,NFRAMES
       READ(99,*)
       READ(99,*) 
        GLYC=0
       DO J2=1,NRES
!        WRITE(6,*)'J2= ',J2

         IF (SEQ(J2).EQ.8) THEN
          READ(99,*)DUMMY,STRUCT(9*(J2-1)+1),STRUCT(9*(J2-1)+2),STRUCT(9*(J2-1)+3)
          STRUCT(9*(J2-1)+4) = STRUCT(9*(J2-1)+1)
          STRUCT(9*(J2-1)+5) = STRUCT(9*(J2-1)+2)
          STRUCT(9*(J2-1)+6) = STRUCT(9*(J2-1)+3)

          READ(99,*)DUMMY,STRUCT(9*(J2-1)+7),STRUCT(9*(J2-1)+8),STRUCT(9*(J2-1)+9)
         ELSE
          READ(99,*)DUMMY,STRUCT(9*(J2-1)+1),STRUCT(9*(J2-1)+2),STRUCT(9*(J2-1)+3)
          READ(99,*)DUMMY,STRUCT(9*(J2-1)+4),STRUCT(9*(J2-1)+5),STRUCT(9*(J2-1)+6)
          READ(99,*)DUMMY,STRUCT(9*(J2-1)+7),STRUCT(9*(J2-1)+8),STRUCT(9*(J2-1)+9)
         ENDIF
        ENDDO

!    This part completes the backbone by calculate in the
!    C-PRIME and NITROGEN suing ideal stero-chemistry for
!    a trans-peptide bond.

         DO JTGRES = 1,NRES      ! CARBON-PRIME
          CPRCORD_X(JTGRES)=0.4436538D0*STRUCT(9*(JTGRES-1)+1) &
                           +0.2352006D0*STRUCT(9*(JTGRES)+1) &
                           +0.3211456D0*STRUCT(9*(JTGRES-1)+7)

          CPRCORD_Y(JTGRES)=0.4436538D0*STRUCT(9*(JTGRES-1)+2) &
                           +0.2352006D0*STRUCT(9*(JTGRES)+2) &
                           +0.3211456D0*STRUCT(9*(JTGRES-1)+8)

          CPRCORD_Z(JTGRES)=0.4436538D0*STRUCT(9*(JTGRES-1)+3) &
                           +0.2352006D0*STRUCT(9*(JTGRES)+3) &
                           +0.3211456D0*STRUCT(9*(JTGRES-1)+9)
         ENDDO

         DO JTGRES = 1,NRES      ! NITROGEN
          NITCORD_X(JTGRES+1)=0.4831806D0*STRUCT(9*(JTGRES-1)+1) &
                             +0.7032820D0*STRUCT(9*(JTGRES)+1) &
                             -0.1864626D0*STRUCT(9*(JTGRES-1)+7)

          NITCORD_Y(JTGRES+1)=0.4831806D0*STRUCT(9*(JTGRES-1)+2) &
                             +0.7032820D0*STRUCT(9*(JTGRES)+2) &
                             -0.1864626D0*STRUCT(9*(JTGRES-1)+8)

          NITCORD_Z(JTGRES+1)=0.4831806D0*STRUCT(9*(JTGRES-1)+3) &
                             +0.7032820D0*STRUCT(9*(JTGRES)+3) &
                             -0.1864626D0*STRUCT(9*(JTGRES-1)+9)
         ENDDO


!     The PDB files are printed to STDOUT here.
          COUNTER=1
          DO JTGRES = 1,NRES
             RES_TYPE=AANAME(SEQ(JTGRES))

       IF (JTGRES .NE. 1 ) THEN
        WRITE(6,56)COUNTER,RES_TYPE,JTGRES, &
            NITCORD_X(JTGRES),NITCORD_Y(JTGRES),NITCORD_Z(JTGRES),JTGRES
56       FORMAT('ATOM',4X,I3,2X,'N ',2X,A3,1X,'A',1X,I3,4X,F8.3,F8.3,F8.3, &
                2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)
         COUNTER=COUNTER+1
       ENDIF

         WRITE(6,52)COUNTER,RES_TYPE,JTGRES,STRUCT(9*(JTGRES-1)+1), &
           STRUCT(9*(JTGRES-1)+2),STRUCT(9*(JTGRES-1)+3),JTGRES
52       FORMAT('ATOM',4X,I3,2X,'CA',2X,A3,1X,'A',1X,I3,4X,F8.3,F8.3,F8.3, &
              2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)
         COUNTER=COUNTER+1

         IF (SEQ(JTGRES).NE.8) THEN
         WRITE(6,53)COUNTER,RES_TYPE,JTGRES,STRUCT(9*(JTGRES-1)+4), &
       STRUCT(9*(JTGRES-1)+5),STRUCT(9*(JTGRES-1)+6),JTGRES
53       FORMAT('ATOM',4X,I3,2X,'CB',2X,A3,1X,'A',1X,I3,4X,F8.3,F8.3,F8.3, &
                2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)
         COUNTER=COUNTER+1
          ENDIF

       IF (JTGRES .NE. NRES)THEN
        WRITE(6,55)COUNTER,RES_TYPE,JTGRES, &
                CPRCORD_X(JTGRES),CPRCORD_Y(JTGRES),CPRCORD_Z(JTGRES),JTGRES
55       FORMAT('ATOM',4X,I3,2X,'C ',2X,A3,1X,'A',1X,I3,4X,F8.3,F8.3,F8.3, &
                2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)
         COUNTER=COUNTER+1
       ENDIF

        WRITE(6,54)COUNTER,RES_TYPE,JTGRES,STRUCT(9*(JTGRES-1)+7), &
       STRUCT(9*(JTGRES-1)+8),STRUCT(9*(JTGRES-1)+9),JTGRES
54      FORMAT('ATOM',4X,I3,2X,'O ',2X,A3,1X,'A',1X,I3,4X,F8.3,F8.3,F8.3, &
                2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)
        COUNTER=COUNTER+1

        ENDDO
        WRITE(6,'(a3)')'END'
      
       ENDDO

      END PROGRAM 
