message("Adding soft_sphere test for GMIN")

GMIN_test_binary(soft_sphere_prog testing.f90)


add_test(
    NAME soft_sphere_E_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND soft_sphere_prog 
)

set_tests_properties(soft_sphere_E_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
  PASS_REGULAR_EXPRESSION "ENERGY CALCULATION PASS"
)


add_test(
    NAME soft_sphere_V_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND soft_sphere_prog 
)

set_tests_properties(soft_sphere_V_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "GRADIENT CALCULATION PASS"
)
