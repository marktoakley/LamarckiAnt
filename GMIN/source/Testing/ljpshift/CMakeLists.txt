message("Adding ljpshift test for GMIN")

GMIN_test_binary(ljpshift_prog testing.f90)


add_test(
    NAME ljpshift_E_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ljpshift_prog 
)

set_tests_properties(ljpshift_E_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
  PASS_REGULAR_EXPRESSION "ENERGY CALCULATION PASS"
)

add_test(
    NAME ljpshift_V_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ljpshift_prog 
)

set_tests_properties(ljpshift_V_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "GRADIENT CALCULATION PASS"
)
