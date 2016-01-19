option(WITH_PYTHON "Generate python interface for GMIN" OFF)

if(WITH_PYTHON)
  find_package(PythonInterp)
  find_package(PythonLibs)
#   set(Boost_USE_STATIC_LIBS        ON)
#   set(Boost_USE_MULTITHREADED      ON)
#   set(Boost_USE_STATIC_RUNTIME    OFF)

  find_package( Boost 1.36.0 COMPONENTS python )
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Error, boost not found. Please install boost or disable the Python interface with -DWITH_PYTHON=no")
    include_directories(${Boost_INCLUDE_DIRS})
  endif(NOT Boost_FOUND)

  include_directories(${PYTHON_INCLUDE_PATH})

  add_custom_command(
    OUTPUT pygmin.F
    DEPENDS ${CMAKE_SOURCE_DIR}/main.F
    COMMAND sed -ne '1,/RUN_TESTS_AFTER_INIT/p' ${CMAKE_SOURCE_DIR}/main.F > pygmin.F
    COMMAND sed -e '/RUN_TESTS_AFTER_INIT/d' -e 's/PROGRAM GMIN/SUBROUTINE PYGMIN_INIT/' -i pygmin.F
    COMMAND echo "      END SUBROUTINE" >> pygmin.F
  )
endif(WITH_PYTHON)

macro(gmin_python_interface)
  if(WITH_PYTHON)
    add_library(gmin_ SHARED pygmin.F python/pygmin.cpp python/wrapper.f90 ${DUMMY_CHARMM} ${DUMMY_AMH} ${DUMMY_AMBER9} ${DUMMY_DMACRYS} ${DUMMY_USERPOT} ${DUMMY_TESTING})
    target_link_libraries(gmin_ gminlib ${MYLAPACK_LIBS} ${Boost_PYTHON_LIBRARY})
    set_target_properties(gmin_ PROPERTIES PREFIX "")
    #set_target_properties($(PY_MODULE} PROPERTIES LINKER_LANGUAGE CXX)
    install_targets(/share/gmin/python gmin_)
    set_target_properties(gmin_ PROPERTIES COMPILE_FLAGS "-DGMIN_PYTHON_MODULE=1")

    if(WITH_DMACRYS)
      include_directories(${CMAKE_BINARY_DIR}/DMACRYSinterface)
      add_library(dmagmin_ SHARED pygmin.F python/pydmagmin.cpp python/dmagmin_wrapper.f90 ${DUMMY_CHARMM} ${DUMMY_AMBER9} ${DUMMY_AMH} ${DUMMY_USERPOT} ${DUMMY_TESTING})
      target_link_libraries(dmagmin_ gminlib dmacrysinterface gminlib ${MYLAPACK_LIBS} ${Boost_PYTHON_LIBRARY})
      set_target_properties(dmagmin_ PROPERTIES PREFIX "")
      #set_target_properties($(PY_MODULE} PROPERTIES LINKER_LANGUAGE CXX)
      install_targets(/share/gmin/python dmagmin_)
      set_target_properties(dmagmin_ PROPERTIES COMPILE_FLAGS "-DGMIN_PYTHON_MODULE=2")
    endif(WITH_DMACRYS)

    if(WITH_AMBER9)
      add_library(ambgmin_ SHARED pygmin.F python/pygmin.cpp python/wrapper.f90 ${DUMMY_CHARMM} ${DUMMY_AMH} ${DUMMY_DMACRYS} ${DUMMY_USERPOT} ${DUMMY_TESTING})
      target_link_libraries(ambgmin_ gminlib AMBER_LIB NAB_LIB dummylib ${MYLAPACK_LIBS} ${Boost_PYTHON_LIBRARY})
      set_target_properties(ambgmin_ PROPERTIES PREFIX "")
      #set_target_properties($(PY_MODULE} PROPERTIES LINKER_LANGUAGE CXX)
      install_targets(/share/gmin/python ambgmin_)
      set_target_properties(ambgmin_ PROPERTIES COMPILE_FLAGS "-DGMIN_PYTHON_MODULE=3")
    endif(WITH_AMBER9)

    if(WITH_OXDNA)
      add_library(oxdnagmin_ SHARED pygmin.F python/pygmin.cpp python/wrapper.f90 ${DUMMY_CHARMM} ${DUMMY_AMBER9} ${DUMMY_AMH} ${DUMMY_DMACRYS} ${DUMMY_TESTING})
      target_link_libraries(oxdnagmin_ gminlib OXDNA ${MYLAPACK_LIBS} ${Boost_PYTHON_LIBRARY})
      set_target_properties(oxdnagmin_ PROPERTIES PREFIX "")
      #set_target_properties($(PY_MODULE} PROPERTIES LINKER_LANGUAGE CXX)
      install_targets(/share/gmin/python oxdnagmin_)
      set_target_properties(oxdnagmin_ PROPERTIES COMPILE_FLAGS "-DGMIN_PYTHON_MODULE=4")
    endif(WITH_OXDNA)
  endif(WITH_PYTHON)

endmacro(gmin_python_interface)
