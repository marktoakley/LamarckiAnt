      PROGRAM EXTRACTTS_CHAINCROSSING
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXATOMS=10000
      CHARACTER(LEN=2) :: DUMMY
      CHARACTER(LEN=3) :: AANAME(20),RES_TYPE
      CHARACTER(LEN=5) :: MEMNAME 
      INTEGER :: SEQ(500),GLY,NRES_AMH, JTGRES

      DOUBLE PRECISION :: STRUCT(MAXATOMS),DDUMMY,DIST
      DOUBLE PRECISION :: CA(MAXATOMS,MAXATOMS)

      INTEGER :: NDUMMY, NFRAMES, J1, NATOMS, J2, J3, NRES, NCOUNT, I_RES
      INTEGER :: GLYC, COUNTER, OPEN_STATUS

      DATA AANAME /'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS', &
         'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'/

!     This program calculated potential chaincrossing events.
!     A few pieces of information are hardwired and may need changing. 
!     The output is printed to STDOUT 
!     1.  The sequence is read from Mike's amh/proteins directory.
!     2.  The path to the extractedts is set to the local directory.

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

      OPEN (UNIT=1,FILE='proteins/'//memname,STATUS='OLD')
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


      CALL SYSTEM("cat pathsample_out | grep 'number of ts extracted=' | awk '{print$6}'  > temptemp") 

       OPEN(1,FILE='temptemp',STATUS='OLD',ACTION='READ',IOSTAT=OPEN_STATUS)
         IF (OPEN_STATUS.NE.0) THEN
           WRITE(6,*) 'FAILURE TO OPEN TEMPTEMP FILE'
           WRITE(6,*) 'NEEDS TO BE IN LOCAL DIRECTORY'
           WRITE(6,*) 'ERROR NUMBER ',OPEN_STATUS
           STOP
         ENDIF
         READ(1,*)NFRAMES
       CLOSE(1)

         IF (NFRAMES.LT.1) THEN
           WRITE(6,*) 'FAILURE IN TEMPTEMP FILE'
           WRITE(6,*) 'NFRAMES LT 1 '
           STOP
         ENDIF

!           WRITE(6,*) 'NFRAMES LT 1 ' , NFRAMES
!           WRITE(6,*) 'NRES ', NRES
       CALL SYSTEM('rm temptemp') 

! 2.  path to the extractedts is set to the local directory

      OPEN (UNIT=99,FILE='extractedts',STATUS='OLD')
      DO J3=1,NFRAMES
       DO J2=1,NRES
!        WRITE(6,*)'J2= ',J2

         IF (SEQ(J2).EQ.8) THEN
          READ(99,*)STRUCT(9*(J2-1)+1),STRUCT(9*(J2-1)+2),STRUCT(9*(J2-1)+3)

          CA(9*(J2-1)+1,J3) = STRUCT(9*(J2-1)+1)
          CA(9*(J2-1)+2,J3) = STRUCT(9*(J2-1)+2)
          CA(9*(J2-1)+3,J3) = STRUCT(9*(J2-1)+3)

!         The C-prime and N postions don't depend on C-beta positions
          READ(99,*)DDUMMY,DDUMMY,DDUMMY

          READ(99,*)STRUCT(9*(J2-1)+7),STRUCT(9*(J2-1)+8),STRUCT(9*(J2-1)+9)
         ELSE
          READ(99,*)STRUCT(9*(J2-1)+1),STRUCT(9*(J2-1)+2),STRUCT(9*(J2-1)+3)
          CA(9*(J2-1)+1,J3) = STRUCT(9*(J2-1)+1)
          CA(9*(J2-1)+2,J3) = STRUCT(9*(J2-1)+2)
          CA(9*(J2-1)+3,J3) = STRUCT(9*(J2-1)+3)

          READ(99,*)STRUCT(9*(J2-1)+4),STRUCT(9*(J2-1)+5),STRUCT(9*(J2-1)+6)
          READ(99,*)STRUCT(9*(J2-1)+7),STRUCT(9*(J2-1)+8),STRUCT(9*(J2-1)+9)
         ENDIF
        ENDDO
       ENDDO

      DO J3=1,NFRAMES
       DO J2=1,NRES - 2
        DO J1=J2 +2 ,NRES
          dist=dsqrt(     (CA(9*(J2-1)+1,J3)       &
                        -  CA(9*(J1-1)+1,J3))**2  &
                        + (CA(9*(J2-1)+2,J3)       &
                        -  CA(9*(J1-1)+2,J3))**2  &
                        + (CA(9*(J2-1)+3,J3)       &
                        -  CA(9*(J1-1)+3,J3))**2  )
          if (dist .lt. 2.0D0)WRITE(6,*)'POTENTIAL CHAIN CROSSSING'
          if (dist .lt. 2.0D0)WRITE(6,*)'FRAME ', J3
          if (dist .lt. 2.0D0)WRITE(6,*)'RESIDUES ', J2, J1
          if (dist .lt. 2.0D0)WRITE(6,*)'DIST ', dist
         ENDDO
        ENDDO
       ENDDO

      END PROGRAM 
