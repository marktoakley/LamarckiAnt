cmake_minimum_required(VERSION 2.6)

project(MYLAPACK)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "dw34@cam.ac.uk")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../CMakeModules)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

#enable_language(CXX)
find_package(FORTRANCOMPILER)

file(GLOB MYLAPACK_SOURCES *.f)

add_library(mylapack ${MYLAPACK_SOURCES})
