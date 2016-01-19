include_directories(${CMAKE_BINARY_DIR})
file(GLOB OPEP_SOURCES *.f90 *.f *.F *.F90)
add_library(OPEP ${OPEP_SOURCES})
add_dependencies(OPEP gminlib)
