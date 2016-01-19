SUBROUTINE RUN_TESTS_AFTER_INIT()
      USE COMMONS
      USE CLASS_OVERLAP
      IMPLICIT NONE
    
      DOUBLE PRECISION POTEL, GRAD(3*NATOMS), GRADtrue(3*NATOMS), Etrue
      DOUBLE PRECISION COORDSNEW(3*NATOMS)
      INTEGER J1, J2
      LOGICAL PASS
      DOUBLE PRECISION accuracy, err, maxerr
      double precision qa, qb, qab
      accuracy=1.0D-10


      CALL OVERLAP_GET_OVERLAP( COORDS(1:3*NATOMS, 1), qa, qb, qab)
      !OPEN(17,FILE='q',status='unknown')
      write(*,*) qa, qb, qab
      !CLOSE(17)

      CALL OVERLAP_GET_OVERLAP2( COORDS(1:3*NATOMS, 1), qa, qb, qab, 0.2D0)
      write(*,*) qa, qb, qab

      !
      !load data from file coords2
      OPEN(17,FILE='coords2',STATUS='UNKNOWN')
      DO J1=1,NATOMS
        READ(17,*) COORDSNEW(3*J1-2), COORDSNEW(3*J1-1), COORDSNEW(3*J1-0)
      ENDDO
      CLOSE(17)

      maxerr = 0.D0

      !check first type of overlap
      CALL OVERLAP_GET_OVERLAP( COORDSnew(1:3*NATOMS), qa, qb, qab)
      write(*,*) qa, qb, qab
      err = abs(qa - 0.9D0)
      if ( err .gt. maxerr ) maxerr = err
      err = abs(qb - 0.78571428571428570D0)
      if ( err .gt. maxerr ) maxerr = err

      !check second type of overlap
      CALL OVERLAP_GET_OVERLAP2( COORDSnew(1:3*NATOMS), qa, qb, qab, 0.2D0)
      write(*,*) qa, qb, qab
      err = abs(qa - 0.96965948569901439D0)
      if ( err .gt. maxerr ) maxerr = err
      err = abs(qb - 0.85540858769456618D0)
      if ( err .gt. maxerr ) maxerr = err

      !check third type of overlap
      CALL OVERLAP_GET_OVERLAP2_R( COORDSnew(1:3*NATOMS), qa, qb, qab, 0.2D0)
      write(*,*) qa, qb, qab
      err = abs(qa - 0.97303401115958743D0)
      if ( err .gt. maxerr ) maxerr = err
      err = abs(qb - 0.66159544348207544D0)
      if ( err .gt. maxerr ) maxerr = err

      write(*,*) "maxerr = ", maxerr
      IF ( maxerr.lt.accuracy ) then
        write(*,*) "overlap approximate"
      ENDIF



      CALL EXIT(1)
END SUBROUTINE
