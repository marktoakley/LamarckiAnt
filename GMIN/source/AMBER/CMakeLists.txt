###########################################################
#
# CMAKE FILE TO GENERATE libamber.a 
# 
############################################################

include_directories(${CMAKE_BINARY_DIR} ${SVN_ROOT}/AMBER/src/sander)

file( GLOB AMBER_SOURCES ${SVN_ROOT}/AMBER/src/sander/*.f* )

# excluded files 
file(GLOB NOT_AMBER_SOURCES 
	${SVN_ROOT}/AMBER/src/sander/egb_backup.f 
	${SVN_ROOT}/AMBER/src/sander/bond.f
	${SVN_ROOT}/AMBER/src/sander/ephi.f
	${SVN_ROOT}/AMBER/src/sander/setbon.f
       	${SVN_ROOT}/AMBER/src/sander/setdih.f
	${SVN_ROOT}/AMBER/src/sander/setvar.f
	${SVN_ROOT}/AMBER/src/sander/dihdup.f
	${SVN_ROOT}/AMBER/src/sander/movecm.f
	${SVN_ROOT}/AMBER/src/sander/angl.f
	)

# exclude files 
list(REMOVE_ITEM AMBER_SOURCES ${NOT_AMBER_SOURCES} )

# we need fortran free-from flag
SET_SOURCE_FILES_PROPERTIES(${AMBER_SOURCES}
  PROPERTIES
  COMPILE_FLAGS "${FORTRAN_FREEFORM_FLAG}"
)

# preprocess 
#message("preprocessing AMBER/nxtsec.f") 
add_custom_command (
	OUTPUT _nxtsec.f 
	DEPENDS ${SVN_ROOT}/AMBER/src/lib/nxtsec.f 
	COMMAND cpp -traditional -P ${SVN_ROOT}/AMBER/src/lib/nxtsec.f > _nxtsec.f
)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES _nxtsec.f )
SET(AMBER_SOURCES ${AMBER_SOURCES} _nxtsec.f)

# build library of local sources (libamber.a)  
add_library(AMBER_LIB ${AMBER_SOURCES} ${HEADERS_COPIED})
ADD_DEPENDENCIES(AMBER_LIB gminlib)

