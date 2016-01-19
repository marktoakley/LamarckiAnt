message("Generating tests")

# defines a macro to run & check a basic gmin quench test
macro (gmin_quench_test testname gmincmd workdir energy)
   add_test(
      NAME ${testname} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${workdir} 
      COMMAND ${gmincmd}
   )

   add_test (
      NAME ${testname}_rounded 
      WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/${workdir} 
      COMMAND ${CMAKE_SOURCE_DIR}/Testing/scripts/check_energy.sh ${energy} 1e-4 
   )

   set_tests_properties(${testname}_rounded PROPERTIES 
      PASS_REGULAR_EXPRESSION "1"
   )

#  add_test (
#     NAME ${testname}_exact 
#     WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/${workdir} 
#     COMMAND ${CMAKE_SOURCE_DIR}/Testing/scripts/check_energy.sh ${energy} 0 
#  )
#  
#  set_tests_properties(${testname}_exact PROPERTIES 
#     PASS_REGULAR_EXPRESSION "1"
#  )
endmacro (gmin_quench_test)

include_directories(${CMAKE_BINARY_DIR})

macro (GMIN_test_binary binname testfile)
  message("Adding test binary: ${binname}")
  add_executable(${binname} ${GMIN_MAIN} ${DUMMY_CHARMM} ${DUMMY_AMH} ${DUMMY_AMBER9} ${DUMMY_DMACRYS} ${DUMMY_USERPOT} ${testfile})
  target_link_libraries(${binname} gminlib ${MYLAPACK_LIBS})
endmacro(GMIN_test_binary)

# add all subdirectories which contain CMakeLists.txt
file(GLOB all_valid_subdirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*/CMakeLists.txt")
foreach(dir ${all_valid_subdirs})
    message(STATUS "path = ${dir}")
    if(${dir} MATCHES "^([^/]*)//CMakeLists.txt")
        string(REGEX REPLACE "^([^/]*)//CMakeLists.txt" "\\1" dir_trimmed ${dir})
        add_subdirectory(${dir_trimmed})
    endif()
endforeach(dir)
