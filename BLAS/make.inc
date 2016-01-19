####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################
#
# SAT: specify list of object files to be excluded from the BLASLIB archive;
# 
BLAS_EXCLUDE_LIST = 
#
#  Modify the FC and FFLAGS definitions to refer to the
#  compiler and desired compiler options for your machine.
#
FC = 
FFLAGS =
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
BLASLIBsingle = libmyblasS.a
BLASLIBdouble = libmyblas.a
BLASLIBcomplex = libmyblasC.a
BLASLIBcomplex16 = libmyblasC16.a
