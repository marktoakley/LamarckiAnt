     PROGRAM GETMETRICORDERED 
     IMPLICIT NONE
     DOUBLE PRECISION METRIC1(1:510000), METRIC2(1:510000),METRIC3(1:510000) 
     DOUBLE PRECISION M1(1:510000), M2(1:510000), M3(1:510000) 
     INTEGER  MNUM(1:510000) 
     DOUBLE PRECISION, ALLOCATABLE :: M1b(:), M2b(:), M3b(:) 
     INTEGER, ALLOCATABLE ::  MNUMb(:) 
     INTEGER I, J, NUMMETRIC, NUM, NMIN, NBASIN, MINE, NRECS
     LOGICAL BASINONLY
     CHARACTER*12 U, W 

! Read parameters required from command line
     CALL GETARG(1,U) ! Number of metrics, 1, 2 or 3
     CALL GETARG(2,W) ! Basin number, 0 to ignore basins
 
     READ(U, *, END=53) NUMMETRIC
     READ(W, *, END=53) NBASIN
     
     IF (NBASIN.LE.0) THEN
       BASINONLY=.FALSE.
     ELSE
       BASINONLY=.TRUE.
     ENDIF
!    PRINT*, 'Number of metrics= ', NUMMETRIC
!    PRINT*, 'Basin number= ', NBASIN  

     METRIC1(:)=HUGE(1.0D0) 
     METRIC2(:)=HUGE(1.0D0)
     METRIC3(:)=HUGE(1.0D0)
     M1(:)=HUGE(1.0D0) 
     M2(:)=HUGE(1.0D0)
     M3(:)=HUGE(1.0D0)
     MNUM(:)=0   
 
     IF (BASINONLY) OPEN(UNIT=2, FILE='basins',status='old')
     OPEN(UNIT=4, FILE='metric',status='old')
     
     I=1
     DO
       IF (I.GT.510000) THEN
         PRINT*, 'Oops - database bigger than default array size'
         PRINT*, 'Increase array size before proceeding'
         STOP
       ENDIF 
       IF (NUMMETRIC.EQ.3) THEN
        READ(4,*, END=22) METRIC1(I), NUM, METRIC2(I), METRIC3(I)
       ELSE IF (NUMMETRIC.EQ.2) THEN 
        READ(4,*, END=22) METRIC1(I), NUM, METRIC2(I)
       ELSE  
        READ(4,*, END=22) METRIC1(I), NUM
       ENDIF 
       I=I+1
     ENDDO

22   NMIN=I-1 
 
     NRECS=0
     DO J=1, NMIN 
      IF (BASINONLY) THEN
        READ(2,*, END=23) MINE, NUM
        IF (MINE.GE.NMIN) EXIT
        IF (NUM.EQ.NBASIN) THEN 
          NRECS=NRECS+1 
          M1(NRECS)=METRIC1(MINE)
          MNUM(NRECS)=MINE
          IF (NUMMETRIC.GE.2) M2(NRECS)=METRIC2(MINE) 
          IF (NUMMETRIC.EQ.3) M3(NRECS)=METRIC3(MINE) 
        ENDIF
      ELSE
          NRECS=NRECS+1 
          M1(NRECS)=METRIC1(J)
          MNUM(NRECS)=J
          IF (NUMMETRIC.GE.2) M2(NRECS)=METRIC2(J) 
          IF (NUMMETRIC.EQ.3) M3(NRECS)=METRIC3(J) 
      ENDIF
     ENDDO
23   ALLOCATE(M1b(1:NRECS))
     ALLOCATE(M2b(1:NRECS))
     ALLOCATE(M3b(1:NRECS))
     ALLOCATE(MNUMb(1:NRECS))
     M1b(1:NRECS)=M1(1:NRECS)
     M2b(1:NRECS)=M2(1:NRECS)
     M3b(1:NRECS)=M3(1:NRECS)
     MNUMb(1:NRECS)=MNUM(1:NRECS)
     
!    DO I=1, NRECS
!     PRINT*, M1b(I), MNUMb(I), M2b(I) 
!    ENDDO

     IF (NUMMETRIC.EQ.3) CALL SORT(M3b, NRECS, MNUMb, M1b, M2b)
!    DO I=1, NRECS
!     IF (NUMMETRIC.EQ.3) WRITE(6,'(F25.10, I12, 2F25.10)') M1b(I), MNUMb(I), M2b(I), M3b(I) 
!    ENDDO
     IF (NUMMETRIC.GE.2) CALL SORT(M2b, NRECS, MNUMb, M1b, M3b)
!    DO I=1, NRECS
!     IF (NUMMETRIC.EQ.2) WRITE(6,'(F25.10, I12, F25.10)') M1b(I), MNUMb(I), M2b(I) 
!     IF (NUMMETRIC.EQ.3) WRITE(6,'(F25.10, I12, 2F25.10)') M1b(I), MNUMb(I), M2b(I), M3b(I) 
!    ENDDO
     CALL SORT(M1b, NRECS, MNUMb, M2b, M3b)
     DO I=1, NRECS
      IF (NUMMETRIC.EQ.3) WRITE(6,'(F25.10, I12, 2F25.10)') M1b(I), MNUMb(I), M2b(I), M3b(I) 
      IF (NUMMETRIC.EQ.2) WRITE(6,'(F25.10, I12, F25.10)') M1b(I), MNUMb(I), M2b(I) 
      IF (NUMMETRIC.EQ.1) WRITE(6,'(F25.10, I12)')  M1b(I), MNUMb(I) 
     ENDDO
     STOP 

53   PRINT*, 'Two arguments needed: [number of metrics], [basin number]'
     END

  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Sorting subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sort(XNUMS, N, DUM1, DUM2, DUM3)
      implicit none
      doubleprecision XNUMS(N)
      integer, allocatable :: ITEMP(:)   
      integer DUM1(N), i, N, minind
      doubleprecision DUM2(N)
      doubleprecision DUM3(N)
      doubleprecision, allocatable :: DTEMP(:)

      Allocate(DTEMP(1:N))
      Allocate(ITEMP(1:N))
!
      do 300 i=1, N
!        IF (i.gt.1) PRINT*, XNUMS(i-1)
!
!        find minind, the index of minimum elements in 
!        XNUMS(i), XNUMS(i+1), ..., XNUMS(n)
!
         call fmin(XNUMS, i, N, minind)
!
!        exchange XNUMS(minind) with XNUMS(i) if necessary
!
         if(XNUMS(minind).LT. XNUMS(i)) then
!            call swap(XNUMS(i), XNUMS(minind))
!            call intswap(DUM1(i), DUM1(minind))
             ITEMP(i)=DUM1(minind)
             IF (I.GT.1) ITEMP(1:i-1)=DUM1(1:i-1)
             ITEMP(i+1:minind)=DUM1(i:minind-1)
             ITEMP(minind+1:N)=DUM1(minind+1:N)
             DUM1(1:N)=ITEMP(1:N)
             DTEMP(i)=DUM2(minind)
             IF (I.GT.1) DTEMP(1:i-1)=DUM2(1:i-1)
             DTEMP(i+1:minind)=DUM2(i:minind-1)
             DTEMP(minind+1:N)=DUM2(minind+1:N)
             DUM2(1:N)=DTEMP(1:N)
!            call swap(DUM2(i), DUM2(minind))
             DTEMP(i)=DUM3(minind)
             IF (I.GT.1) DTEMP(1:i-1)=DUM3(1:i-1)
             DTEMP(i+1:minind)=DUM3(i:minind-1)
             DTEMP(minind+1:N)=DUM3(minind+1:N)
             DUM3(1:N)=DTEMP(1:N)
!            call swap(DUM3(i), DUM3(minind))
            DTEMP(1:N)=XNUMS(1:N)
            IF (i.GT.1) DTEMP(1:i-1)=XNUMS(1:i-1)
            DTEMP(i)=XNUMS(minind)
            DTEMP(i+1:minind)=XNUMS(i:minind-1)
            DTEMP(minind+1:N)=XNUMS(minind+1:N)
            XNUMS(1:N)=DTEMP(1:N)
         endif
 300  continue
      return
      end

!
!     find the index of the minmum elements in 
!     sub-array XNUMS(i), XNUMS(i+1), ..., XNUMS(n)
!
!
      subroutine fmin(XNUMS, i, N, minind)
      integer i, N, minind
      doubleprecision XNUMS(N)
!
      doubleprecision xmin

      xmin   = XNUMS(i)
      minind = i
      do 200 j=i+1, N
         if(XNUMS(j) .LT. xmin) then
            xmin = XNUMS(j)
            minind = j
         endif
 200  continue
!
      return
      end

!
!     x <--> y
!
      subroutine swap(x, y)
      doubleprecision x, y, temp
!
      temp = x
      x    = y
      y    = temp
!
      return
      end
!
!     x <--> y
!
      subroutine intswap(x, y)
      integer x, y, temp
!
      temp = x
      x    = y
      y    = temp
!
      return
      end
 
