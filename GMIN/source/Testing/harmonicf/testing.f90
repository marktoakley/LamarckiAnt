SUBROUTINE RUN_TESTS_AFTER_INIT()
      USE COMMONS, ONLY : BOXLX, BOXLY, BOXLZ, HARMONICSTR, NATOMS, COORDS
      IMPLICIT NONE
    
      INTEGER J1, J2
      DOUBLE PRECISION R(3*NATOMS), V(3*NATOMS), E(NATOMS)
      DOUBLE PRECISION accuracy, err, maxerr, temp, X(3), X0(3), &
     &                 Etrue(Natoms), Vtrue(3*NATOMS)
      accuracy=1.0D-10

      !set parameters
      HARMONICSTR=1.5D0
      BOXLX = 100.D0
      BOXLY = 100.D0
      BOXLZ = 100.D0
      x0(1) = 0.D0
      x0(2) = 0.D0
      x0(3) = 0.D0
      X(1) = 0.1D0
      X(2) = 0.2D0
      X(3) = 0.3D0

      ! calculate energy only
      E(:)=0.D0
      V(:)=0.D0
      DO J1=1,NATOMS
        CALL HARMONICFIELD(COORDS((J1-1)*3+1,1), X0, V((J1-1)*3+1), E(J1))
        R(J1) = SQRT( &
     &               (COORDS(3*(J1-1)+1,1)-X0(1))**2+&
     &               (COORDS(3*(J1-1)+2,1)-X0(2))**2+&
     &               (COORDS(3*(J1-1)+3,1)-X0(3))**2)
      ENDDO

      OPEN(17,FILE='E',status='unknown')
      DO J1=1,NATOMS
        WRITE(17,*) R(J1), E(J1)
      ENDDO
      CLOSE(17)
      OPEN(17,FILE='Etrue',status='unknown')
      DO J1=1,NATOMS
        READ(17,*) TEMP, ETRUE(J1)
      ENDDO
      CLOSE(17)
      call maxdiff(natoms, e, etrue, maxerr)
      IF ( maxerr.eq.0.d0 ) then
        write(*,*) "energy exact"
      ENDIF
      IF ( maxerr.lt.accuracy ) then
        write(*,*) "energy approximate"
      else
        write(*,*) "fail: maxerr ", maxerr
      ENDIF


      ! test gradient, but only every 10th atom to save space in the repository
      !write V
      OPEN(17,FILE='V',status='unknown')
      DO J1=1,NATOMS
        !WRITE(17,*) R(j1), 
!    &  coords(3*(J1-1)+1,1), coords(3*(J1-1)+2,1), coords(3*(J1-1)+3,1),&
        WRITE(17,*) V(3*(J1-1)+1), V(3*(J1-1)+2), V(3*(J1-1)+3)
      ENDDO
      CLOSE(17)
      !read correct V
      OPEN(17,FILE='Vtrue',status='unknown')
      DO J1=1,NATOMS
        READ(17,*) Vtrue(3*(J1-1)+1), Vtrue(3*(J1-1)+2), Vtrue(3*(J1-1)+3)
      ENDDO
      CLOSE(17)
      call maxdiff(3*NATOMs, V, Vtrue, maxerr)
      IF ( maxerr.eq.0.d0 ) then
        write(*,*) "potential exact"
      ENDIF
      IF ( maxerr.lt.accuracy ) then
        write(*,*) "potential approximate"
      ENDIF


      CALL EXIT(1)
END SUBROUTINE

SUBROUTINE MAXDIFF( N, V, VTRUE, MAXERR )
      implicit none
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: V(N), VTRUE(N)
      double precision, intent(out) :: maxerr
      double precision err
      integer j1
      maxerr=0.D0
      DO J1=1,N
        ERR = ABS(V(J1) - VTRUE(J1))
        IF ( ERR .GT. MAXERR) THEN
          MAXERR = ERR
        ENDIF
      END DO
END SUBROUTINE
