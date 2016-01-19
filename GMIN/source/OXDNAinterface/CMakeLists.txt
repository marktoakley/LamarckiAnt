#cmake_minimum_required(VERSION 2.8)

#project(GMIN)

#set(PROJECT_VERSION "0.0")
#set(PROJECT_CONTACT "todo@email")

# if not specified by user, the standard build type is release
#if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
#  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
#endif(NOT CMAKE_BUILD_TYPE)

#enable_language(CXX)

file(GLOB OXDNA_SOURCES *.cpp *.c *.f90)
file(GLOB NOT_OXDNA_SOURCES test.cpp)

list(REMOVE_ITEM OXDNA_SOURCES ${NOT_OXDNA_SOURCES})

include_directories(${CMAKE_BINARY_DIR})
add_library(OXDNA ${OXDNA_SOURCES})
set_property(TARGET OXDNA PROPERTY COMPILE_DEFINITIONS NOCUDA)
add_dependencies(OXDNA gminlib)

#add_executable(oxdna_test test.cpp)
#target_link_libraries(oxdna_test OXDNA)

