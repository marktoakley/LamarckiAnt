include_directories(${CMAKE_BINARY_DIR})

file(GLOB DMACRYSINTERFACE_SOURCES *.f *.f90 ${SVN_ROOT}/DMACRYS/*.f90 ${SVN_ROOT}/DMACRYS/*.c)
add_library(dmacrysinterface ${DMACRYSINTERFACE_SOURCES})
ADD_DEPENDENCIES(dmacrysinterface gminlib)
