message("Adding harmonicf test for GMIN")

GMIN_test_binary(harmonicf_prog testing.f90)


add_test(
    NAME harmonicf_E_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND harmonicf_prog 
)

set_tests_properties(harmonicf_E_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "energy approximate"
)


add_test(
    NAME harmonicf_V_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND harmonicf_prog 
)

set_tests_properties(harmonicf_V_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "potential approximate"
)
