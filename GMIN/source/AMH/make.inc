####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################
#
# SAT: specify list of object files to be excluded from the BLASLIB archive;
# 
#  Modify the FC and FFLAGS definitions to refer to the
#  compiler and desired compiler options for your machine.
#
# FC = ifort
# FFLAGS= -132 -C -g -traceback
#  FFLAGS= -132 -Vaxlib -O3 # for ifc
#  FFLAGS= -132 -O0 # for ifc
#FC = f95
#FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O0
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH = ar
ARCHFLAGS = cr
ARCHDELFLAGS = d
RANLIB = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
AMHLIB = ~/svn/GMIN/source/libamh.a
