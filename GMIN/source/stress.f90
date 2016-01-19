SUBROUTINE CALC_STRESS(X,EPOT)
  !
  USE COMMONS, ONLY : STRESST, STRESS, NATOMS, MYUNIT, &
       MGUPTAT, MIEFT
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: EPOT
  !
  INTEGER :: I,J
  DOUBLE PRECISION :: GRAD(3*NATOMS)
  !
  IF(.NOT.STRESST) THEN
     WRITE(MYUNIT,'(A)') &
          'calc_stress> Should not be here!'
     RETURN
  ELSE
     STRESS(:,:,:) = 0.0D0 ! Initialise
  ENDIF
  !
  IF(MGUPTAT) THEN
     CALL MGUPTA(X, GRAD, EPOT, .TRUE., .FALSE., .TRUE.)
  ELSE
     WRITE(MYUNIT,'(A)') &
          'calc_stress> Stress calculation not implemented for current potential.'
  ENDIF
  !
  IF(MIEFT) CALL MIEF(X,GRAD,EPOT,.TRUE.,.TRUE.)
  !
  WRITE(MYUNIT,'(A)') 'stress> Overall stress tensor:'
  DO I=1,3
     WRITE(MYUNIT,'(3(1X,E15.8))') (STRESS(0,I,J),J=1,3)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE CALC_STRESS
