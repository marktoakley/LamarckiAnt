!
! You need this to exclude pairwise interactions in AMBER
!
PROGRAM MAIN
  IMPLICIT NONE
  INTEGER TNATOMS 
  INTEGER, ALLOCATABLE :: EXCLUDEDNO(:), ATOMINBODY(:), ATOMLIST(:,:), GROUPMEMBERS(:,:)
  INTEGER, ALLOCATABLE :: EXCLUDEDLIST(:), EXCLUDEDPERATOM(:,:)
  INTEGER J1, J2, J3, J4, IOSTATUS
  CHARACTER(LEN=5) RECORDNAME
  LOGICAL ADDITIONALEXCLUSION

  GROUPMEMBERS(1:TNATOMS,1:TNATOMS) = 0

  PRINT *, "ENTER THE NUMBER OF ATOMS IN THE SYSTEM"
  READ (UNIT = 5, *) TNATOMS
  ALLOCATE (EXCLUDEDNO(TNATOMS))
  ALLOCATE (ATOMINBODY(TNATOMS))
  ALLOCATE (ATOMLIST(TNATOMS,TNATOMS))
  ALLOCATE (GROUPMEMBERS(TNATOMS,TNATOMS))
  ALLOCATE (EXCLUDEDLIST(TNATOMS*TNATOMS))
  ALLOCATE (EXCLUDEDPERATOM(TNATOMS,TNATOMS))
  
! hk286 - input
  OPEN(UNIT=1,FILE='excluded_list_no',STATUS='UNKNOWN')
  DO J1 = 1, TNATOMS/10
     READ(1,*,IOSTAT=iostatus) EXCLUDEDNO((J1-1)*10+1:J1*10)
  ENDDO
  READ(1,*,IOSTAT=iostatus) EXCLUDEDNO( (TNATOMS/10)*10+1 : TNATOMS )

  PRINT *, "FINISH READING"
  
  J2 = 0
  DO J1 = 1, TNATOMS
     J2 = J2 + EXCLUDEDNO(J1)
  ENDDO
  CLOSE(1)
 
  OPEN(UNIT=1,FILE='excluded_list',STATUS='UNKNOWN')
  DO J1 = 1, J2/10
     READ(1,*,IOSTAT=iostatus) EXCLUDEDLIST((J1-1)*10+1:J1*10)
  ENDDO
  READ(1,*,IOSTAT=iostatus) EXCLUDEDLIST((J2/10)*10+1:J2)

  PRINT *, "FINISH READING"

  J4 = 0
  DO J1 = 1, TNATOMS
     DO J3 = 1, EXCLUDEDNO(J1)
        J4 = J4 + 1
        EXCLUDEDPERATOM(J1,J3) = EXCLUDEDLIST(J4)
     ENDDO
  ENDDO
  CLOSE(1)

  OPEN(UNIT=1,FILE='rbodyconfig',STATUS='UNKNOWN')
  READ (1,*,IOSTAT=iostatus) RECORDNAME, ATOMINBODY(1)
  DO J1 = 1, ATOMINBODY(1)
     READ(1,*,IOSTAT=iostatus) ATOMLIST(1,J1)
  ENDDO
  
  DO J1 = 1, ATOMINBODY(1)
     DO J4 = J1+1, ATOMINBODY(1)
        ADDITIONALEXCLUSION = .TRUE.
        DO J3 = 1, EXCLUDEDNO(ATOMLIST(1,J1))
           IF (EXCLUDEDPERATOM(ATOMLIST(1, J1),J3).EQ. ATOMLIST(1,J4)) THEN
              ADDITIONALEXCLUSION = .FALSE.
           ENDIF
        ENDDO
        IF (ADDITIONALEXCLUSION .EQ. .TRUE.) THEN
           EXCLUDEDNO(ATOMLIST(1,J1)) = EXCLUDEDNO(ATOMLIST(1,J1)) + 1
           DO J3 = 1, EXCLUDEDNO(ATOMLIST(1,J1))-1
              IF (ATOMLIST(1,J4) < EXCLUDEDPERATOM(ATOMLIST(1, J1),EXCLUDEDNO(ATOMLIST(1,J1))-J3)) THEN
                 EXCLUDEDPERATOM(ATOMLIST(1,J1),EXCLUDEDNO(ATOMLIST(1,J1))-J3+1) = EXCLUDEDPERATOM(ATOMLIST(1,J1),EXCLUDEDNO(ATOMLIST(1,J1))-J3)
                 EXCLUDEDPERATOM(ATOMLIST(1,J1),EXCLUDEDNO(ATOMLIST(1,J1))-J3) = ATOMLIST(1,J4)
              ELSEIF (J3 .EQ. 1) THEN
                 EXCLUDEDPERATOM(ATOMLIST(1,J1),EXCLUDEDNO(ATOMLIST(1,J1))-J3+1) = ATOMLIST(1,J4)
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  CLOSE(1)

  PRINT *, "NOW WRITING"

  
  OPEN(UNIT=1,FILE='excluded_list_no_new',STATUS='UNKNOWN')
  DO J1 = 1, TNATOMS/10
     WRITE(1,'(10(I8))') EXCLUDEDNO((J1-1)*10+1:J1*10)
  ENDDO
  WRITE(1,'(10(I8))') EXCLUDEDNO((TNATOMS/10)*10+1 : TNATOMS )
  CLOSE(1)  

  OPEN(UNIT=1,FILE='excluded_list_new',STATUS='UNKNOWN')
  J4 = 0
  DO J1 = 1, TNATOMS
     DO J3 = 1, EXCLUDEDNO(J1)
        J4 = J4 + 1
        EXCLUDEDLIST(J4) = EXCLUDEDPERATOM(J1,J3)
     ENDDO
  ENDDO

  DO J1 = 1, J4/10
     WRITE(1,'(10(I8))') EXCLUDEDLIST((J1-1)*10+1:J1*10)
  ENDDO
  WRITE(1,'(10(I8))') EXCLUDEDLIST((J4/10)*10+1:J4)
  CLOSE(1)

END PROGRAM MAIN