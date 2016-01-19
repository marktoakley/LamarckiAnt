SUBROUTINE RUN_TESTS_AFTER_INIT()
    ! this space should be used to test a given subroutine or function.
    !
    ! e.g.
    !
    ! CALL MY_POTENTIAL(MYCOORDS, ENERGY)
    !
    ! IF ( ENERGY .EQ. ENERGY_TRUE ) THEN

    print *,"Test successful"

    ! ENDIF
    CALL EXIT(1)
END SUBROUTINE
