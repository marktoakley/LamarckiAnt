###########################################################
#
# MAIN CMAKE FILE FOR GMIN 
# ========================
#
# Usage: 
#   
#   $ mkdir ~/svn/GMIN/builds/pgf90
#   $ cd ~/svn/GMIN/builds/pgf90
#   $ FC=pgf90 cmake ../../source 
#   $ make -j8 <target>    
#              -- if <target> not specified, all active ones are built  
#
# Good to know: 
#   1. VERBOSE=1 make   
#         -- gives the full compilation command for each file, useful for
#         debugging  
#
#   2. To compile with MPI 
#         FC=mpif90 cmake ../source -DCOMPILER_SWITCH=pgi 
#         make GMIN 
#       
#   3. GUI for CMakeCache.txt  
#         ccmake .    (after running cmake command) 
#         e.g.: 
#              set WITH_MPI to ON -> Press [c] -> Press [e] -> Press [g]  
#           equivalent command line version
#            $ cmake . -DWITH_MPI=yes 
# 
#   4. Out-of-source build only 
#         do not run cmake in ~svn/GMIN/source    
#
#   5. Tested with 
#        $ cmake --version 
#        cmake version 2.8.7.20120416-g31e7
#
#   6. To do a clean build, first delete CMakeCache.txt, and if that doesn't
#   help, delete the whole directory!
# 
############################################################

cmake_minimum_required(VERSION 2.8)
# 2.8 required particularly for the charmm-interfaced targets.

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(GMIN)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "dw34@cam.ac.uk")
get_filename_component(SVN_ROOT ${GMIN_SOURCE_DIR}/../.. ABSOLUTE)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${SVN_ROOT}/CMakeModules)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Release Debug Debug_slow RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
set(CMAKE_Fortran_MODULE_DIRECTORY ${GMIN_BINARY_DIR}/modules)
#include_directories(${CMAKE_CURRENT_BINARY_DIR}/modules)
#include_directories(${CMAKE_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/PythonGMIN.cmake)

file(GLOB DUMMY_AMH     amhdummy.f90)
file(GLOB DUMMY_CHARMM  chdummy.f)
file(GLOB DUMMY_AMBER9  amber9dummy.f90)
file(GLOB DUMMY_AMBER12 ${SVN_ROOT}/AMBER12/amber12_interface.F90)
file(GLOB DUMMY_CUDA    modcudalbfgs.F90)
file(GLOB DUMMY_USERPOT dummy_userpot.f90)
file(GLOB DUMMY_DMACRYS dummy_dmacrys.f90)
file(GLOB DUMMY_TESTING dummy_testing.f90)
file(GLOB DUMMY_OPTIM   dummyoptim.f90)

set(ALL_DUMMIES ${DUMMY_AMH}
                ${DUMMY_CHARMM}
                ${DUMMY_AMBER9}
                ${DUMMY_AMBER12}
                ${DUMMY_CUDA}
                ${DUMMY_DMACRYS}
                ${DUMMY_USERPOT}
                ${DUMMY_TESTING} 
                ${DUMMY_OPTIM} )

# Glob all the sources
file(GLOB GMIN_LIB_SOURCES *.f
                           *.f90
                           *.F
                           *.F90 
                           *.c )

file(GLOB NOT_GMIN_SOURCES myblas.f
                           mylapack.f
                           # These files are generated later
                           display_version.f90
                           porfuncs.f90 )

# Due to a compiler bug in ifort 13.1.3, we can't use -O3 for genrigid.f90
# Investigations continue...
if( ${COMPILER_SWITCH} STREQUAL "ifort" )
  SET_SOURCE_FILES_PROPERTIES( genrigid.f90 PROPERTIES COMPILE_FLAGS -O2 )
endif ( ${COMPILER_SWITCH} STREQUAL "ifort" )


list(REMOVE_ITEM GMIN_LIB_SOURCES ${NOT_GMIN_SOURCES}
                                  ${ALL_DUMMIES} )

# Add generated source files 
set(GMIN_LIB_SOURCES ${GMIN_LIB_SOURCES}
                     display_version.f90
                     porfuncs.f90 )

# ---- porfuncs.f90 TODO: John is writing a replacement
add_custom_command (
   OUTPUT porfuncs.f90 
   COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/porfuncs.sh ${COMPILER_SWITCH} > porfuncs.f90)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES porfuncs.f90)

# ---- display_version.f90 
add_custom_command (
  OUTPUT display_version.f90 
  COMMAND bash ${SVN_ROOT}/SCRIPTS/all/display_version.sh > display_version.f90
)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES display_version.f90)

# Collect all the extra source files into one place. It's just dummies at the moment.
set(ALL_EXTRA_SOURCES ${ALL_DUMMIES})

########################
# LIBRARIES START HERE #
########################

# Make an extras library for dummies, this is necessary because some of the
# files in gminlib use modules made by extralib
add_library(extralib ${ALL_EXTRA_SOURCES})
set_module_dir(extralib)
set_target_properties(extralib PROPERTIES LINKER_LANGUAGE "Fortran") 
set_target_properties(extralib PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")

# Make a gmin library
add_library(gminlib ${GMIN_LIB_SOURCES})
set_module_dir(gminlib)
set_module_depends(gminlib extralib)
set_target_properties(gminlib PROPERTIES LINKER_LANGUAGE "Fortran") 

# If you try to optimise this code your computer will rot
add_subdirectory(libmbpol)
target_link_libraries(gminlib libmbpol)

find_package(MYBLAS)
find_package(MYLAPACK)

# MPI flag 
option(WITH_MPI "Enable MPI compilation" OFF)
if(WITH_MPI)
  # append -DMPI to compile definitions
  add_definitions(-DMPI) 
  message("building with MPI: \n  ${CMAKE_Fortran_FLAGS}")   
endif(WITH_MPI)

add_definitions(-D_SVN_ROOT_='${SVN_ROOT}')

##########################
# EXECUTABLES START HERE #
##########################

# standard GMIN
add_executable(GMIN main.F)
set_target_properties(GMIN PROPERTIES LINKER_LANGUAGE "Fortran")
set_module_depends(GMIN gminlib
                        extralib )
target_link_libraries(GMIN gminlib
                           extralib
                           ${MYLAPACK_LIBS})

# AMHGMIN
option(WITH_AMH "Enable AMHGMIN compilation" OFF)
if(WITH_AMH)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
# AMH has issues with the e_write subroutine. The amh_interfaces_mod has the interface for e_write, but
# does not define it and so causes problems in linking. 
# At the moment, I just link in the dummy module instead.
#  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_AMH})
  add_subdirectory(${GMIN_SOURCE_DIR}/AMH)
  add_executable(AMHGMIN ${EXTRA_SOURCES})
  set_module_dir(AMHGMIN)
  set_module_depends(AMHGMIN gminlib)
  set_module_dir(amh)
  set_module_depends(amh gminlib
                         extralib )
  set_target_properties(AMHGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(AMHGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(AMHGMIN gminlib 
                                amh
                                ${MYLAPACK_LIBS})
endif(WITH_AMH)

# A9GMIN  
option(WITH_AMBER9 "Enable Amber GMIN with AMBER9 compilation" OFF)
if(WITH_AMBER9)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_AMBER9})
  find_package(NAB)
  add_subdirectory(AMBER) # temporary folder to copy files from sander;  
  add_executable(A9GMIN main.F
                        ${EXTRA_SOURCES} )
  set_module_dir(A9GMIN)
  set_module_depends(A9GMIN gminlib)
  set_module_dir(AMBER_LIB)
  set_module_depends(AMBER_LIB gminlib)
  set_target_properties(A9GMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(A9GMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(A9GMIN gminlib
                               AMBER_LIB
                               NAB_LIB
                               ${MYLAPACK_LIBS} )
endif(WITH_AMBER9)

# For A12GMIN or CUDAGMIN, we need to build and add the AMBER 12 interface.
if(WITH_AMBER12 OR WITH_CUDA)
  add_subdirectory(${SVN_ROOT}/AMBER12 AMBER12_interface)
endif()

# A12GMIN  
option(WITH_AMBER12 "Enable Amber GMIN with AMBER12 compilation" OFF)
# A12GMIN and CUDAGMIN both use the same Amber interface. Thus, we can only build the CPU version or the
# GPU version at any given time.
if(WITH_AMBER12 AND NOT WITH_CUDA)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  add_executable(A12GMIN main.F
                         ${EXTRA_SOURCES})
  set_module_dir(A12GMIN)
  set_module_depends(A12GMIN gminlib
                             amber12_base )
  set_target_properties(A12GMIN PROPERTIES LINKER_LANGUAGE "Fortran") 
  set_target_properties(A12GMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_CUDA") 
  target_link_libraries(A12GMIN gminlib
                                amber12_base
                                ${MYLAPACK_LIBS} )
endif(WITH_AMBER12 AND NOT WITH_CUDA)

# CUDAGMIN
option(WITH_CUDA "Enable CUDAGMIN compilation" OFF)
if(WITH_CUDA)
  if(NOT ${COMPILER_SWITCH} STREQUAL "ifort" AND NOT ${COMPILER_SWITCH} STREQUAL "gfortran")
    message( FATAL_ERROR "The cuda version of AMBER can only be compiled with gnu or Intel compilers. " )
  endif (NOT ${COMPILER_SWITCH} STREQUAL "ifort" AND NOT ${COMPILER_SWITCH} STREQUAL "gfortran") 
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  add_subdirectory(AMBER12CUDAinterface)
  add_subdirectory(${SVN_ROOT}/cuda_lbfgs cuda_lbfgs)
  add_executable(CUDAGMIN main.F
                          ${EXTRA_SOURCES})
  set_module_dir(CUDAGMIN)
  set_module_depends(CUDAGMIN gminlib 
                              amber12_base )
  set_target_properties(CUDAGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(CUDAGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};CUDA")
  target_link_libraries(CUDAGMIN gminlib 
                                 amber12cuda
                                 cuda_lbfgs
                                 amber12_base
                                 ${MYLAPACK_LIBS} )
endif(WITH_CUDA)

#DMACRYS
option(WITH_DMACRYS "Enable DMAGMIN compilation (DMACRYS needs to be present!)" OFF)
if(WITH_DMACRYS)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_DMACRYS})
  add_subdirectory(DMACRYSinterface)
  add_executable(DMAGMIN main.F
                         ${EXTRA_SOURCES} )
  set_module_dir(DMAGMIN)
  set_module_depends(DMAGMIN gminlib)
  set_module_dir(dmacrysinterface)
  set_module_depends(dmacrysinterface gminlib)
  set_target_properties(DMAGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(DMAGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(DMAGMIN gminlib
                                dmacrysinterface
                                ${MYLAPACK_LIBS})
endif(WITH_DMACRYS)

#OXDNA
option(WITH_OXDNA "Enable OXDNAGMIN compilation" OFF)
if(WITH_OXDNA)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_USERPOT})
  add_subdirectory(OXDNAinterface)
  add_executable(OXDNAGMIN main.F
                           ${EXTRA_SOURCES})
  set_module_dir(OXDNAGMIN)
  set_module_depends(OXDNAGMIN gminlib)
  set_module_dir(OXDNA)
  set_module_depends(OXDNA gminlib)
  set_target_properties(OXDNAGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(OXDNAGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(OXDNAGMIN gminlib
                                  OXDNA 
                                  ${MYLAPACK_LIBS} )
endif(WITH_OXDNA)

#SPIN interface
option(WITH_SPIN "Enable SPINGMIN compilation" OFF)
if(WITH_SPIN)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_USERPOT})
  add_subdirectory(SPINinterface)
  add_executable(SPINGMIN main.F
                          ${EXTRA_SOURCES})
  set_module_dir(SPINGMIN)
  set_module_depends(SPINGMIN gminlib)
  set_module_dir(SPIN)
  set_module_depends(SPIN gminlib)
  set_target_properties(SPINGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(SPINGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(SPINGMIN gminlib
                                 SPIN
                                 ${MYLAPACK_LIBS} )
endif(WITH_SPIN)

#OPEP interface
option(WITH_OPEP "Enable OPEPGMIN compilation" OFF)
if(WITH_OPEP)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_USERPOT})
  add_subdirectory(OPEPinterface)
  add_executable(OPEPGMIN main.F
                          ${EXTRA_SOURCES} )
  set_module_dir(OPEPGMIN)
  set_module_depends(OPEPGMIN gminlib)
  set_module_dir(OPEP)
  set_module_depends(OPEP gminlib)
  set_target_properties(OPEPGMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(OPEPGMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(OPEPGMIN gminlib
                                 OPEP
                                 ${MYLAPACK_LIBS})
endif(WITH_OPEP)

# C35GMIN 
option(WITH_CHARMM35 "Enable C35GMIN compilation" OFF)
if(WITH_CHARMM35)
  if(WITH_CHARMM36)
    message( FATAL_ERROR "cannot build CHARMM35 and CHARMM36 at the same time" )
  endif(WITH_CHARMM36)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_CHARMM})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CHARMMboth)
  add_executable(C35GMIN main.F
                         myblas.f
                         mylapack.f
                         ${EXTRA_SOURCES})
  set_module_dir(C35GMIN)
  set_module_depends(C35GMIN gminlib)
  set_module_dir(charmm35_interface_lib)
  set_module_depends(charmm35_interface_lib gminlib)
  set_target_properties(C35GMIN PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(C35GMIN PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  # add the CHARMM library as a dependency
  add_dependencies(C35GMIN CHARMM_WAS_BUILT)
  # Link in the CHARMM libraries and the CHARMM35 GMIN interface
  target_link_libraries(C35GMIN gminlib
                                mylapack
                                myblascharmm
                                charmm35_interface_lib
                                ${CHARMM_LIBS} )
endif(WITH_CHARMM35)

# C36GMIN, doesn't do anything just yet, just here as a template 
option(WITH_CHARMM36 "Enable C36GMIN compilation" OFF)
if(WITH_CHARMM36)
  message( FATAL_ERROR "CHARMM36 not yet implemented with CMake" )
#  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CHARMMboth)
#  add_executable(C36GMIN ${GMIN_MAIN} ${GMIN_PREPROCESS} myblas.f mylapack.f ${DUMMY_AMBER9} ${DUMMY_AMH}
#${DUMMY_DMACRYS} ${DUMMY_USERPOT} ${DUMMY_TESTING} ${DUMMY_AMBER12} ${DUMMY_CUDA})
#  target_link_libraries(C36GMIN gminlib mylapack myblascharmm charmm ${CHARMM_LIBS})
endif(WITH_CHARMM36)

# Example section for USERPOT entry, please do not commit CMAkelists.txt if you make changes here!
# I'm working on a more general solutioin where you can give the path to the userpotenital source
# directory
# add_subdirectory(path/to/userpot/library)
# add_executable(USERGMIN main.F.f ${DUMMY_CHARMM} ${DUMMY_AMH} ${DUMMY_AMBER9} ${DUMMY_DMACRYS})
# target_link_libraries(GMIN gminlib ${MYLAPACK_LIBS} USERPOTLIBRARY <whatever else is needed>)

# Testing framework
option(WITH_TESTING "Enable automatic generation of tests (Testing folder has to be present!)" OFF)
if(WITH_TESTING)
  enable_testing()
  add_subdirectory(Testing)
endif(WITH_TESTING)

gmin_python_interface()

# Test for duplicates.
# In the standard build we first merge the object files together in a library
# (gminlib) to avoid duplicate compilation for multiple targets.  A
# drawback of this is that duplicate symbols (duplicate names of subroutines,
# functions, etc.) are simply ignored.  This is a separate target which builds
# GMIN directly from the sources, specifically to check for duplicate
# symbols.
option(WITH_DUPLICATE_TESTING "Compile GMIN in such a way that duplicate symbols are flagged" OFF)
if(WITH_DUPLICATE_TESTING)
  add_executable(GMIN_TEST_DUPLICATES main.F
                                      ${GMIN_SOURCES}
                                      ${ALL_EXTRA_SOURCES} )
  set_module_dir(GMIN_TEST_DUPLICATES)
  set_module_depends(GMIN_TEST_DUPLICATES gminlib)
  set_target_properties(GMIN_TEST_DUPLICATES PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(GMIN_TEST_DUPLICATES PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(GMIN_TEST_DUPLICATES ${MYLAPACK_LIBS}
                                             libmbpol )
  add_dependencies(GMIN_TEST_DUPLICATES gminlib) # so I don't have to manually specify all the other dependencies
endif(WITH_DUPLICATE_TESTING)
