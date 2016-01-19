message("Adding example test for GMIN")

#this tells cmake to compile a new binary called GMIN_test_example_prog using
#the subroutine in testing.f90 the binary name must be unique
GMIN_test_binary(GMIN_test_example_prog testing.f90)

#this tells cmake that you're making a new test
#the test name must be unique
add_test(
    NAME GMIN_test_example 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND GMIN_test_example_prog 
)

#This tells cmake how to tell if your test is successful.  cmake will search
#standard out for the string "Test successful".  If found, cmake will
#assume the test was successful.  See testing.f90
set_tests_properties(GMIN_test_example PROPERTIES 
#    FAIL_REGULAR_EXPRESSION "failed"
    PASS_REGULAR_EXPRESSION "Test successful"
)

