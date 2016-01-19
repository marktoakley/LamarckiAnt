include_directories(${CMAKE_BINARY_DIR})
file(GLOB DIPOLE_SOURCES *.f90 *.f)
add_library(SPIN ${DIPOLE_SOURCES})
add_dependencies(SPIN gminlib)
