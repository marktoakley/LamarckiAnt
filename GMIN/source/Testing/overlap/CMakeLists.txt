message("Adding overlap tests for GMIN")

GMIN_test_binary(overlap_prog testing.f90)

add_test(
    NAME overlap_approx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND overlap_prog 
)

set_tests_properties(overlap_approx PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "overlap approximate"
)
