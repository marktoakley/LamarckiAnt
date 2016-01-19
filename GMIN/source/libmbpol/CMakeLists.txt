cmake_minimum_required(VERSION 2.6)

project(libmbpol)
# Francesco Paesani's MBPOL flexible water potential
set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "jdf43@cam.ac.uk")

# optimisation of polynomials takes forever
# would be better as a CMAKE variable
set(CMAKE_CXX_FLAGS_RELEASE -O0)
set(CMAKE_CXX_FLAGS_RELEASE -fPIC)

file(GLOB LIBMBPOL_SOURCES *.cpp *.h)
add_library(libmbpol ${LIBMBPOL_SOURCES})
