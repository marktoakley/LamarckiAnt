include_directories(${CMAKE_BINARY_DIR})
file(GLOB AMH_SOURCES *.f *.f90)
add_library(amh ${AMH_SOURCES})
ADD_DEPENDENCIES(amh gminlib)
