####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################
#
#  Modify the FC and FFLAGS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.
#
FC = ifc
FFLAGS = 
NOOPT = -O0
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
LAPACKLIBsingle = libmylapackS.a
LAPACKLIBdouble = libmylapackD.a
LAPACKLIBcomplex = libmylapackC.a
LAPACKLIBcomplex16 = libmylapackC16.a
LAPACKLIBselection = libmylapack.a
