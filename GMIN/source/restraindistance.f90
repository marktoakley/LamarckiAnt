SUBROUTINE READRESTRAINL()
  USE commons
  IMPLICIT NONE

  CHARACTER(LEN=100) CHECK1
  INTEGER :: iostatus, J1


! hk286 > determine no of pairs
  RESTRAINLNOPAIR = 0
  OPEN(UNIT=222,FILE='restraindistances',status='old')
  DO
     READ(222,*,IOSTAT=iostatus) CHECK1
     IF (iostatus<0) THEN
        CLOSE(222)
        EXIT
     ELSE 
        RESTRAINLNOPAIR = RESTRAINLNOPAIR + 1
     ENDIF
  END DO
  CLOSE(222)

  ALLOCATE (RESTRAINLDIST(RESTRAINLNOPAIR))
  ALLOCATE (RESTRAINLPAIRS(RESTRAINLNOPAIR,2))

  OPEN(UNIT=222,FILE='restraindistances',status='old')
  DO J1 = 1, RESTRAINLNOPAIR 
     READ(222,*) RESTRAINLPAIRS(J1,1), RESTRAINLPAIRS(J1,2), RESTRAINLDIST(J1)
  END DO
  CLOSE(222)


END SUBROUTINE READRESTRAINL


SUBROUTINE RESTRAINLPOTENTIAL(X,GRAD,EREAL,GRADT,SECT)

  USE commons
  IMPLICIT NONE
  
  DOUBLE PRECISION :: X(3*NATOMS), GRAD(3*NATOMS), EREAL, DX(3), DIST
  LOGICAL :: GRADT, SECT
  INTEGER :: J1, JP1, JP2

  DO J1 = 1, RESTRAINLNOPAIR 
 
     JP1 = RESTRAINLPAIRS(J1,1)
     JP2 = RESTRAINLPAIRS(J1,2)
     DX = X(3*JP1-2:3*JP1) - X(3*JP2-2:3*JP2)
     DIST = DOT_PRODUCT(DX,DX)
     EREAL = EREAl + 0.5D0 * RESTRAINLK * (DIST - RESTRAINLDIST(J1)**2)**2
     IF (GRADT .EQV. .TRUE.) THEN
        GRAD(3*JP1-2:3*JP1) = GRAD(3*JP1-2:3*JP1) +  2.0D0 * RESTRAINLK * (DIST - RESTRAINLDIST(J1)**2) * DX
        GRAD(3*JP2-2:3*JP2) = GRAD(3*JP2-2:3*JP2) -  2.0D0 * RESTRAINLK * (DIST - RESTRAINLDIST(J1)**2) * DX
     END IF
     IF (SECT .EQV. .TRUE.) THEN
        PRINT *, "restraindistance > Hessians not yet implemented"
        STOP
     ENDIF
  ENDDO
  
END SUBROUTINE RESTRAINLPOTENTIAL
