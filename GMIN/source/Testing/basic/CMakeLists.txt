message("Adding example test for GMIN")

#include_directories(${CMAKE_BINARY_DIR})
add_executable(test_rotations testrot.f90 
  ${SVN_ROOT}/GMIN/source/rotations.f90 ${SVN_ROOT}/GMIN/source/dprand.f 
  ${SVN_ROOT}/GMIN/source/vec3.f90 )

add_test(
    NAME  rotation_tests 
    COMMAND test_rotations
)

set_tests_properties(rotation_tests PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "Test failed"
    PASS_REGULAR_EXPRESSION "All tests successful"
)

