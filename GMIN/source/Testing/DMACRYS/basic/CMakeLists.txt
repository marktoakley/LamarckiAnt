message("Adding example test for GMIN")

#include_directories(${CMAKE_BINARY_DIR})
add_executable(mimimum_reduction_tests reduction.f90 ${SVN_ROOT}/GMIN/source/DMACRYSinterface/cell_reduction.f90 ${SVN_ROOT}/GMIN/source/dprand.f ${SVN_ROOT}/GMIN/source/vec3.f90)

add_test(
    NAME  mimimum_reduction_tests 
    COMMAND mimimum_reduction_tests
)

set_tests_properties(mimimum_reduction_tests PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "Test failed"
    PASS_REGULAR_EXPRESSION "All tests successful"
)

