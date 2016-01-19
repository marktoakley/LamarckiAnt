SUBROUTINE RUN_TESTS_AFTER_INIT()
      USE COMMONS
      IMPLICIT NONE
    
      DOUBLE PRECISION POTEL, GRAD(3*NATOMS), GRADtrue(3*NATOMS), Etrue
      INTEGER J1, J2, J3, NWRITE
      LOGICAL PASS
      DOUBLE PRECISION accuracy, err, maxerr, v1, v2, v3
      DOUBLE PRECISION TIMESTART, TIMEEND
      LOGICAL ETEST, VTEST
      LOGICAL Epass, Vpass
      epass = .true.
      vpass = .true.


      !
      ! do test with FREEZE and RESTRICTREGION
      !
      write(*,*) ""
      write(*,*) "do test with FREEZE and RESTRICTREGION"
      CALL TEST_EV(etest, vtest)
      if (.not. etest) epass = .false.
      if (.not. vtest) vpass = .false.




      !
      !turn RESTRICTREGION off and repeat
      !
      write(*,*) ""
      write(*,*) "turn RESTRICTREGION off and repeat"
      RESTRICTREGION=.FALSE.
      CALL TEST_EV(etest, vtest)
      if (.not. etest) epass = .false.
      if (.not. vtest) vpass = .false.



      !
      !turn FREEZE off and repeat
      !
      write(*,*) ""
      write(*,*) "turn FREEZE off and repeat"
      FREEZE=.FALSE.
      CALL TEST_EV(etest, vtest)
      if (.not. etest) epass = .false.
      if (.not. vtest) vpass = .false.


      WRITE(*,*) ""
      IF (EPASS) THEN
         WRITE(*,*) "ENERGY CALCULATION PASS"
      ELSE
         WRITE(*,*) "ENERGY CALCULATION FAIL"
      ENDIF
      IF (VPASS) THEN
         WRITE(*,*) "GRADIENT CALCULATION PASS"
      ELSE
         WRITE(*,*) "GRADIENT CALCULATION FAIL"
      ENDIF
      WRITE(*,*) ""






!      ! test how long it takes to call ljpshift_old
!      CALL MYCPU_TIME(TIMESTART)
!      do j1=1,1000
!         CALL LJPSHIFT_old(COORDS, GRAD, POTEL, .true., .false.)  
!      enddo
!      CALL MYCPU_TIME(TIMEEND)
!      WRITE(*,*) "TIME ljpshift ", TIMEEND-TIMESTART
!
!      ! test how long it takes to call ljpshift
!      CALL MYCPU_TIME(TIMESTART)
!      do j1=1,1000
!         CALL LJPSHIFT(COORDS, GRAD, POTEL, .true., .false.)  
!      enddo
!      CALL MYCPU_TIME(TIMEEND)
!      WRITE(*,*) "TIME ljpshift ", TIMEEND-TIMESTART
!
!      ! test how long it takes to call soft_sphere
!      CALL MYCPU_TIME(TIMESTART)
!      do j1=1,1000
!         CALL soft_sphere_pot(COORDS, GRAD, POTEL, .true., .false.)  
!      enddo
!      CALL MYCPU_TIME(TIMEEND)
!      WRITE(*,*) "TIME soft_sphere ", TIMEEND-TIMESTART

      CALL EXIT(1)
END SUBROUTINE

SUBROUTINE TEST_EV(etest, vtest) 
      USE COMMONS
      IMPLICIT NONE
    
      DOUBLE PRECISION POTEL, GRAD(3*NATOMS), GRADtrue(3*NATOMS), Etrue
      INTEGER J1, J2, J3, NWRITE
      LOGICAL PASS
      DOUBLE PRECISION err, Eerr, Vmaxerr, v1, v2, v3, accuracy
      logical etest, vtest

      accuracy=1.0D-10

      ! calculate energy only
      CALL LJPSHIFT(COORDS, GRAD, POTEL, .false., .false.)  

      OPEN(17,FILE='E',status='unknown')
        WRITE(17,*) POTEL
      CLOSE(17)
      OPEN(17,FILE='Etrue',status='unknown')
        read(17,*) etrue
      CLOSE(17)
      Eerr = POTEL - ETRUE

      ! calculate energy and gradient
      CALL LJPSHIFT(COORDS, GRAD, POTEL, .true., .false.)  

      ! test gradient, but only every 10th atom to save space in the repository
      !write V
      OPEN(17,FILE='V',status='unknown')
      nwrite = 0
      DO J1=1,NATOMS
         IF (.NOT. FROZEN(J1) ) THEN
            WRITE(17,*) J1, GRAD(3*J1-2), GRAD(3*J1-1), GRAD(3*J1-0)
            NWRITE = NWRITE + 1
         ENDIF
      END DO
      CLOSE(17)
      !read correct V
      Vmaxerr=0.d0
      OPEN(17,FILE='Vtrue',status='unknown')
      DO J1=1,NWRITE
        READ(17,*) J2, V1, V2, V3
        GRADTRUE(3*J2-2) = V1
        GRADTRUE(3*J2-1) = V2
        GRADTRUE(3*J2-0) = V3
        DO J3=1,3
          ERR = ABS(GRAD(3*(J2-1)+J3) - GRADTRUE(3*(J2-1)+J3))
          IF ( ERR .GT. VMAXERR) THEN
            VMAXERR = ERR
          ENDIF
        ENDDO
      END DO
      CLOSE(17)

      IF ( Eerr.eq. 0.D0 ) then
         write(*,*) "   energy calculation exact"
      else 
         IF ( ABS(Eerr).lt.accuracy ) then
            write(*,*) "   energy calculation approximate"
         ENDIF
         write (*,*) "      energy error ", potel-ETRUE
      ENDIF

      IF ( Vmaxerr.eq.0.d0 ) then
         write(*,*) "   gradient exact"
      else
         IF ( Vmaxerr.lt.accuracy ) then
            write(*,*) "   gradient approximate"
         ENDIF
         write (*,*) "      gradient maxerr ", Vmaxerr
      ENDIF

      etest = err .lt. accuracy
      vtest = ( Vmaxerr.lt.accuracy )
END SUBROUTINE TEST_EV
