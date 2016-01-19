#  
# Default target gives GMIN
# binaries are now produced in the bin directory!
# 'make chgmin' produces CGMIN
# 'make AMBGMIN' produces AMBGMIN
# 'make clean' removed all object and library files, and any binaries in the bin directory
#
#  DONT FORGET - mpirun
#              - add MPI keyword to data to run with mpirun
#              - change CHARMM library to the appropriate version
#
#  64 bit CGMIN openmpi executable works with BSPT LJ31. 
#  
#  Tried linking to CHARMM compiled with 64 bit openmpi. Reduced optimisation options.
#  BSPT metenk run works fine on one node and more than one node when linked to openmpi charmm build.
#  Also works when linked to 64 bit charmm library compiled with pgf90.
#  -C flag cannot be used with any CHARMM build due to CHARMM intentional out of bounds code.
#
#  pgf90 7.1.6 works for 64 bit lam and openmpi
#  pgf90 7.1.6 works for 64 bit openmpi with CHARMM library ~/svn/CHARMM35
#        for -Mextend -O3 -Munroll -Mnoframe flags in 
#        ~/svn/CHARMM35/build/UNX/Makefile_gnu and modules
#  pgi/64/7.1/6 and mpi/openmpi/64/pgi71/1.2.6
#
#  pgf90 7.1.6 64 bit lam cannot be tested due to a dependency
#
#  openmpi works with ifort on mek-quake:
#  ifort_em64t/10.0.025 mpi/openmpi/64/intel10/1.2.6 icc_em64t/10.0.023
#
#  openmpi fails on clust with module:
#
# ifort/64/10.1/015 icc/64/10.1/015 mpi/openmpi/64/intel10/1.2.6
#

#include build_config

SVN_ROOT = $(realpath $(PWD)/../..)
#SHELL = /bin/bash -O expand_aliases
#SVN_REV_NO = "$(shell module load subversion/1.7.7; svn --version --quiet)"

GPROG =	../bin/GMIN
CPROG = ../bin/CGMIN
AMHPROG = ../bin/AMHGMIN
AMBPROG = ../bin/AMBGMIN
CUDAPROG = ../bin/CUDAGMIN
#QUIPPROG = ../bin/QUIPGMIN
AMB12PROG = ../bin/AMB12GMIN

OBJS1 = porfuncs.o commons.o countatoms.o modamb.o modcharmm.o modmxatms.o modhess.o modamber9.o grouprotation.o dampedmove.o moves.o sanity.o pot_params.o \
	operations.o overlap.o polirmod.o ps_intra_polir_b.o swmod.o chirality.o grafpack.o mbpol.o precision.o libmbpol.a gauss.o rad.o dihedral.o twist.o \
	rotamer_move.o
OBJS2 =	centre.o genrigid.o rotations.o finalio.o modconsts_trans_97.o modconsts.o dist.o tryexchange.o \
	io1.o keywords.o  main.o mc.o mcruns.o morse.o \
	potential.o quench.o dprand.o saveit.o seed.o \
	sort.o sort2.o sort3.o sort4.o takestep.o mycpu_time.o trans.o \
	finalq.o symmetry.o symmetrycsm.o ptgrp.o eigsrt.o SiSW.o taboo.o reseed.o newinertia.o supermc.o pgsym.o pgsym_mod.o\
	tosifumi.o ortho.o compress.o mylbfgs.o mymylbfgs.o input.o ddfpmin.o dlnsrch.o cgmin.o linmin.o \
	brent.o mnbrak.o dbrent.o f1dim.o zwischen.o hsmove.o PachecoC60.o AT.o EAMLJ_sub.o \
	Pbglue.o wenzel.o odesd.o capsid.o rigidfuncs.o tip.o pah.o strand.o \
	SW.o qmod.o ljpbin.o fdm.o dftb.o ljpshift.o dzugutov.o ljcoulomb.o \
	binary_id_swaps.o homoref.o homoref_addons.o enperms.o multiperm.o mie_field.o boxcentroid.o \
	QALCSearch.o QALCS_surface.o QALCS_symmetry.o QALCS_symmetry_mod.o atom_label_swaps.o atom_label_flips.o stress.o \
	fd.o fedor.o welch.o glj_yukawa.o BGupta.o BLJcluster.o BLJcluster_nocut.o stock.o Farkas.o getorbits.o \
	sc.o MSC.o MGupta.o MLJ.o Zetterling.o MSorig.o MSorigc.o MStrans.97.o convert.o \
	frausi.o p46merdiff.o g46merdiff.o lj.o modperm.o modf1com.o mododesd.o EAMal.o Alglue.o Mgglue.o \
	Gupta.o orient.o Natb.o sticky.o enumerate.o minperm.o minpermdist.o LB2.o \
	dgetrf.o dgetri.o reorient.o thomson.o Q4.o basinsampling.o tether.o tetherfunc.o BLN.o \
	newmindist.o centrecom.o qorderparam_blj.o qorderparam_lj.o bspt.o GMINdump.o quad.o \
	rigidbaa.o checkd.o capbin.o dbpg.o dbptd.o dmblmorse.o dmblpy.o gbcalamitic.o gbdiscotic.o gem.o gbddp.o linrod.o \
        lwotp.o	newcapsid.o newpah.o msgayberne.o mstbin.o multstock.o paha.o pahw99.o pap.o ptstst.o papbin.o papjanus.o silane.o \
        multpaha.o newtip.o patchy.o asaoos.o stockaa.o tetrahedra.o waterpdc.o waterpkz.o takestepmsgb.o dipolarmorse.o \
	gay-berne.o JM.o projI.o projIh.o model1.o FS.o vgw.o \
	mysd.o qdtest.o qdtest2.o MB.o dummyoptim.o Ackland_metals.o Ackland_wrapper.o DF1.o \
	read_cmd_args.o display_version.o file_manager.o output.o chaperonin.o bulkmindist.o perc.o \
	ga_modules.o ga_main.o ga_select.o ga_bln.o ga_cluster.o dummy_dmacrys.o dummy_userpot.o vec3.o \
	dummy_testing.o harmonic.o utils.o minpermdistrbcom.o rbperm.o msaveit.o soft_sphere.o bhpt.o \
	initialization.o chiro.o sandbox.o hbondmatrix.o py.o Energy_1d_APBC.o atomlists.o \
	Energy_1d_PBC.o Energy_2d_APBC.o Energy_2d_PBC.o Energy_3d_APBC.o Energy_3d_PBC.o \
	cell_lists.o cell_lists_binary.o \
	neighbor_list.o neighbor_list_bin.o bspt_one_atom.o ratio.o neighbor_list_moveone.o dodecamorse.o \
	tryexchange_independence.o dprand_universal.o dlsode.o ttm3f.o restraindistance.o mW_wrapper.o mindist.o \
	bspt_utils.o selectmove.o CSMinit.o Gthomson.o percsphere.o ellipticintegral.o reservoir_utils.o biased_torsional_sampling.o random_normal.o \
	dftbC.o detsym.o symutils.o ljgh.o morsegh.o GLJ.o phi4model.o dihedralrotation.o parse_pot_params.o inertia.o mlowest.o qcipot.o \
   make_conpot.o congrad.o intlbfgs.o


GMIN_DIR = ${SVN_ROOT}/GMIN/source
OBJS = ${OBJS1} ${OBJS2}
GENF90FILES = porfuncs.f90
#QUIPINTERFACE = gmin_quip_wrapper.o
CUDAINTERFACE = modcudalbfgs.o
CHDUM = chdummy.o
AMHDUM = amhdummy.o
AMB12DUM = ${SVN_ROOT}/AMBER12/amber12_interface.o
AMB9DUM = amber9dummy.o
QUIPDUM = gmin_quipdummy.o
CUDADUM = modcudalbfgs.o
VPATH = .:AMH:.
#LDFLAGS = -L.
DEFS =
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P

INCLUDE_COMMAND = -I${GMIN_DIR} -I${GMIN_DIR}/AMH -I${SVN_ROOT}/SHARED_SOURCES -I${SVN_ROOT}/AMBER12

#QUIP library files
#QUIPLIBS = libfullquip.a

#CUDA stuff
NVCC = nvcc
NVCCFLAGS = -O3 -arch compute_35
CUDASRC = /home/${USER}/svn/GMIN/source/cuda_lbfgs
CUDAFILES = $(CUDASRC)/lbfgs/timer_cuda.o $(CUDASRC)/projects/setup_cuda.o $(CUDASRC)/projects/lj_cuda.o $(CUDASRC)/lbfgs/lbfgs_cuda.o
CUDAOBJS = timer_cuda.o setup_cuda.o lj_cuda.o lbfgs_cuda.o
CUDAINCDIRS = -I/usr/local/cuda/include -I$(CUDASRC)/include -I$(CUDASRC)/include/CudaLBFGS -I$(CUDASRC)/lbfgs -I/usr/local/cuda/include
CUDALIBS = -L/usr/local/cuda-5.5/targets/x86_64-linux/lib -lcudart  -L/usr/lib -lcuda -L/usr/local/cuda-5.5/targets/x86_64-linux/lib -lcublas

#
# location of the source for the AMBER9 interface
#
##### AMBER9 ##################### 
AMB9SRC = ../../AMBER
SANDERSRC = ../../AMBER/src/sander
NABSRC = ../../NAB/src

# Define which CHARMM version you want to use with GMIN. Just uncomment the option you want. 
# for CHARMM 31
#  CTYPE=C31
# for CHARMM 35 (recommended)
#  CTYPE=C35
# do you want to use SCC-DFTB with CHARMM35 (will not work for CHARMM31!)?
#  QMTYPE=DFTB
# No CHARMM
CTYPE=C0
#
# set CHARMM build environment variable and compilation command
# Only used for C31 and C35
#
# For DFTB compilation, need to add 'T' to this line to replace QUANTUM with SCC-DFTB
ifeq (${QMTYPE},DFTB)
    chcommand = ./install.com gnu large ${FC} T keepo keepf OPT
else
    chcommand = ./install.com gnu medium ${FC} keepo keepf OPT
endif

ifeq (${CTYPE},C31)
#### start CHARMM 31 
    C31SRC = /home/${USER}/svn/CHARMM31
##  C31SRC = /home/${USER}/svn/charmm31.mpi
    BLAS_EXCLUDE_LIST = dnrm2.o daxpy.o dcopy.o ddot.o
    FCMDIR = ${C31SRC}/source/fcm
    SRC31 = charmm_main.src # energy.src
    EXTRAS = myblas.o mylapack.o
    PREFLX = ${C31SRC}/tool/prefx_gnu
#
# The following line should probably never be uncommented!
##  PREFLX = ${C31SRC}/tool/preflx
    PREFDIR = ${C31SRC}/build/gnu
##### CHARMM 31 block end 
    LIBDIR31=${C31SRC}/lib/gnu
    CHOBJS31 = $(LIBDIR31)/help.o $(LIBDIR31)/iniall.o $(LIBDIR31)/miscom.o $(LIBDIR31)/usersb.o
    CHLIBS31 = $(LIBDIR31)/adumb.a \
	        $(LIBDIR31)/flucq.a $(LIBDIR31)/cadint.a $(LIBDIR31)/cheq.a $(LIBDIR31)/cff.a $(LIBDIR31)/correl.a $(LIBDIR31)/dimb.a \
	        $(LIBDIR31)/emap.a $(LIBDIR31)/dynamc.a $(LIBDIR31)/energy.a $(LIBDIR31)/gamint.a $(LIBDIR31)/gukint.a \
		   $(LIBDIR31)/gener.a $(LIBDIR31)/image.a $(LIBDIR31)/io.a $(LIBDIR31)/machdep.a $(LIBDIR31)/manip.a $(LIBDIR31)/mbond.a \
		   $(LIBDIR31)/mc.a $(LIBDIR31)/minmiz.a $(LIBDIR31)/misc.a $(LIBDIR31)/mmff.a $(LIBDIR31)/molvib.a $(LIBDIR31)/nbonds.a \
		   $(LIBDIR31)/pert.a $(LIBDIR31)/quantum.a $(LIBDIR31)/rxncor.a $(LIBDIR31)/shapes.a $(LIBDIR31)/solvation.a \
		   $(LIBDIR31)/util.a $(LIBDIR31)/vibran.a libcharmm.a
#### end CHARMM 31 
endif

ifeq (${CTYPE},C35)
#### start CHARMM 35 
    C35SRC = /home/${USER}/svn/CHARMM35
    BLAS_EXCLUDE_LIST = dnrm2.o daxpy.o dcopy.o ddot.o
    FCMDIR = ${C35SRC}/source/fcm
    SRC35 = charmm_main.src # energy.src
    EXTRAS = myblas.o mylapack.o
    PREFLX = ${C35SRC}/tool/prefx_gnu
    PREFDIR = ${C35SRC}/build/gnu
    LIBDIR35=${C35SRC}/lib/gnu
    CHOBJS35 = $(LIBDIR35)/help.o $(LIBDIR35)/iniall.o $(LIBDIR35)/miscom.o $(LIBDIR35)/usersb.o
ifeq (${QMTYPE},DFTB)
    CHLIBS35 = $(LIBDIR35)/adumb.a \
	        $(LIBDIR35)/flucq.a $(LIBDIR35)/cadint.a $(LIBDIR35)/cheq.a $(LIBDIR35)/cff.a $(LIBDIR35)/correl.a $(LIBDIR35)/dimb.a \
	        $(LIBDIR35)/emap.a $(LIBDIR35)/dynamc.a $(LIBDIR35)/energy.a $(LIBDIR35)/gamint.a $(LIBDIR35)/gukint.a \
		   $(LIBDIR35)/gener.a $(LIBDIR35)/image.a $(LIBDIR35)/io.a $(LIBDIR35)/machdep.a $(LIBDIR35)/manip.a $(LIBDIR35)/mbond.a \
		   $(LIBDIR35)/mc.a $(LIBDIR35)/minmiz.a $(LIBDIR35)/misc.a $(LIBDIR35)/mmff.a $(LIBDIR35)/molvib.a $(LIBDIR35)/nbonds.a \
		   $(LIBDIR35)/pert.a $(LIBDIR35)/quantum.a $(LIBDIR35)/rxncor.a $(LIBDIR35)/shapes.a $(LIBDIR35)/solvation.a \
		   $(LIBDIR35)/util.a $(LIBDIR35)/vibran.a $(LIBDIR35)/sccdftb.a $(LIBDIR35)/sccdftbint.a libcharmm.a
else
    CHLIBS35 = $(LIBDIR35)/adumb.a \
	        $(LIBDIR35)/flucq.a $(LIBDIR35)/cadint.a $(LIBDIR35)/cheq.a $(LIBDIR35)/cff.a $(LIBDIR35)/correl.a $(LIBDIR35)/dimb.a \
	        $(LIBDIR35)/emap.a $(LIBDIR35)/dynamc.a $(LIBDIR35)/energy.a $(LIBDIR35)/gamint.a $(LIBDIR35)/gukint.a \
		   $(LIBDIR35)/gener.a $(LIBDIR35)/image.a $(LIBDIR35)/io.a $(LIBDIR35)/machdep.a $(LIBDIR35)/manip.a $(LIBDIR35)/mbond.a \
		   $(LIBDIR35)/mc.a $(LIBDIR35)/minmiz.a $(LIBDIR35)/misc.a $(LIBDIR35)/mmff.a $(LIBDIR35)/molvib.a $(LIBDIR35)/nbonds.a \
		   $(LIBDIR35)/pert.a $(LIBDIR35)/quantum.a $(LIBDIR35)/rxncor.a $(LIBDIR35)/shapes.a $(LIBDIR35)/solvation.a \
		   $(LIBDIR35)/util.a $(LIBDIR35)/vibran.a libcharmm.a
endif
#### end CHARMM 35 

endif

###################################### COMPILERS AND COMPILER FLAGS ###################################### 
#
######### The Portland Group Compiler Technology Fortran 90 compiler {{
# FC = pgf90
## FC = mpif77  # for lam - don't forget to uncomment MPI!
## FC = mpif90  # for mpich and openmpi - don't forget to uncomment MPI!
## DEFS = -DMPI
## The usual flags for AMBGMIN:
## FFLAGS= -Mextend -O3 -Mvect=assoc,cachesize:1024000,recog,transform
## These are the CHARMM31 flags for mpif90 64 bit library.
## It is ESSENTIAL to use the same flags as for the CHARMM build!!!!
##
## FFLAGS= -Mextend -O3 -Munroll -Mnoframe # used this for normal stuff
# FFLAGS= -Mextend -O0 -Mnoframe
## FFLAGS= -Mextend -O0 -C -g -traceback
## FFLAGS= -Mextend -O0 -g -traceback -pg ## used this for debugging with pgdbg swo24
## FFLAGS= -Mextend -g -traceback ## used this for debugging with pgdbg swo24
## Debugging flags
## FFLAGS= -Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Melf -Mpgicoff -traceback
## FFLAGS= -Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mpgicoff -traceback
## flags for AMBER9
##
## FREEFORMAT_FLAG= -Mfree
## EXTRA_FLAGS=-module
## MOD_LOCATION= -module $(@D)
##
# NOOPT = -O0 -Mextend
# SEARCH_PATH =  -I..
## the double reference to -lblas seems to be needed!
# LIBS = libmyblas.a libmylapack.a libmyblas.a #libmbpol.a
# LDFLAGS= -L. -lstdc++
# SWITCH=pgi
## }}}
####### end of The Portland Group Compiler Technology Fortran 90 compiler }}}
#
#### NAGWare Fortran 
## use this compiler for nag/64/5.1
## FC = f95
## DEFS+=-DNAG
## use this compiler for nag/64/5.2
#FC = nagfor
#FFLAGS = -132 -Bdynamic  -maxcontin=3000 -kind=byte -mismatch_all -O4 
## this line is for garden variety debugging 
## FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all -ieee=stop
## this line is for thorough but slow debugging 
## FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte -ieee=stop
#NOOPT= -O0 -132  -kind=byte
#SEARCH_PATH = -I..
#LDFLAGS= -L. -lstdc++
##MOD_LOCATION= -mdir $(@D)
#SWITCH=nag
##
## the double reference to -lblas seems to be needed!
##
#LIBS = libmylapack.a libmyblas.a 
##
## flags for AMBER9
##
##FREEFORMAT_FLAG= -free
##EXTRA_FLAGS=-I
##
####### end of NAGWare Fortran 95 compiler flags 
#
###########################################################
# Intel compilers 
#
# GMIN + CHARMM35 - tested with the options which contain "##" in front 
# 
 FC = ifort
# FC = mpif77 
# FC = mpif90  
# DEFS = -DMPI 
# DEFS = -DCUDA 	# CUDAGMIN compiled with the Makefile currently cannot be used with AMBER - in the meantime use CMake
##### ifort debugging flags 
# FFLAGS= -132 -C -g -traceback -debug full
# FFLAGS= -132 -O0 -g -traceback -fpe:0 -check all
# FFLAGS= -132 -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
# FFLAGS= -132 -g -debug all -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds \
#      -O0 -traceback -debug full -check uninit -stand f03 -assume realloc_lhs -check all,noarg_temp_created \
#      -traceback -warn all -fstack-protector -assume protect_parens -implicitnone
##### ifort optimization flags
 FFLAGS= -132 -O3
# FFLAGS= -132 -Vaxlib -O3 # for ifc
# FFLAGS= -132 -O3 -ip -static # -ipo or -fast doesnt work 
# FFLAGS= -132 -O4
## FFLAGS= -O2 -extend_source -g -traceback
 NOOPT= -132 -O0
 SWITCH=ifort
 SEARCH_PATH= -I.
 SEARCH_PATH= -I..
 LDFLAGS= -L. -lstdc++
 LIBS = libmyblas.a libmylapack.a libmyblas.a  
 FREEFORMAT_FLAG= -free
 EXTRA_FLAGS=-I
 MOD_LOCATION= -module $(@D)
#
# End of intel compilers 
## Gfortran  
#
# FC = gfortran
# FC = mpif90
# DEFS = -DMPI
#DEFS = -DCUDA		# CUDAGMIN compiled with the Makefile currently cannot be used with AMBER - in the meantime use CMake
# FFLAGS= -ffixed-line-length-265 -O3 -ffree-line-length-0
# FFLAGS= -ffixed-line-length-132 -O3 -ftree-vectorize -ffree-line-length-0
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -ffree-line-length-0
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic -ffree-line-length-0
# NOOPT= -O0 -ffixed-line-length-132
# SWITCH=gfortran
# SEARCH_PATH =  -I..
# LIBS = libmylapack.a libmyblas.a
# LDFLAGS = -LBLAS -LLAPACK
# LDFLAGS = -L. -lstdc++
# FREEFORMAT_FLAG= -ffree-form
# EXTRA_FLAGS=-I
#
#  End Gfortran 
# 
###############################################################################################
# Pass the subversion revision number into the code to print it in the output
DEFS+=-D_SVN_ROOT_="'$(SVN_ROOT)'"
DEFS+=-DDUMMY_CUDA
DEFS+=-DDUMMY_AMBER12
#DEFS+=-D_SVN_REV_NO_="'$(SVN_REV_NO)'"
   #DEFS+=-DSVNVERSION="`./svn_revision.sh`"
###################################### RULES AND TARGETS ###################################### 
.SUFFIXES:
.SUFFIXES: .o .f .F .f90 .F90

.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

.F90.f90:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

.f90.o:
#	$(FC) $(FFLAGS) ${SEARCH_PATH} -module . -o $@ -c $<
	$(FC) $(FFLAGS) ${SEARCH_PATH} ${INCLUDE_COMMAND} ${MOD_LOCATION} -o $@ -c $<

.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} ${INCLUDE_COMMAND} ${MOD_LOCATION} -o $@ -c $<

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) ${CUDAINCDIRS} -c $< 

default: $(GPROG)

# vr274 not needed anymore since there is a rule for .F90->.f90 now
#tryexchange.o:
#	$(CPP) $(CPFLAGS) $(DEFS) tryexchange.F > tryexchange.f90
#	$(FC) $(FFLAGS) ${SEARCH_PATH} -c tryexchange.f90

GMIN: $(GPROG)
gmin: $(GPROG)

CHGMIN: $(CPROG)
chgmin: $(CPROG)

$(GPROG): $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(OBJS) $(EXTRAS) 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS) $(OBJS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(LDFLAGS) $(LIBS)

ifeq (${CTYPE},C31)
$(CPROG): $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(OBJS) $(EXTRAS) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ ${CHOBJS31} $(EXTRAS) $(OBJS) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(LDFLAGS) $(LIBS) \
        ${CHLIBS31} ${CHLIBS31} ${CHLIBS31}
endif
ifeq (${CTYPE},C35)
$(CPROG): $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(OBJS) $(EXTRAS) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ ${CHOBJS35} $(EXTRAS) $(OBJS) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(LDFLAGS) $(LIBS) \
        ${CHLIBS35} ${CHLIBS35} ${CHLIBS35}
endif

AMBGMIN: $(AMBPROG) 
ambgmin: $(AMBPROG) 

$(AMBPROG): $(OBJS) $(EXTRAS) $(CHDUM) $(AMHDUM) $(AMB12DUM) $(CUDADUM) $(QUIPDUM) libamber.a libnab.a
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(OBJS) $(EXTRAS) $(CHDUM) $(AMHDUM) $(AMB12DUM) $(CUDADUM) $(QUIPDUM) $(LDFLAGS) libamber.a libnab.a $(LIBS) 

CUDAGMIN: $(CUDAPROG)
cudagmin: $(CUDAPROG)

$(CUDAPROG): $(CUDAFILES) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(CUDAINTERFACE) $(OBJS) $(EXTRAS)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS)  $(CUDAOBJS) $(CUDAINTERFACE) $(OBJS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPDUM) $(LDFLAGS) $(LIBS) $(CUDALIBS)

#QUIPGMIN: $(QUIPPROG)
#quipgmin: $(QUIPPROG)

#$(QUIPPROG): $(QUIPLIBS) $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPINTERFACE)
#	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(AMHDUM) $(QUIPINTERFACE) $(LDFLAGS) $(QUIPLIBS) -static -llapack -lopenblas 


AMHGMIN: $(AMHPROG)
amhgmin: $(AMHPROG)

$(AMHPROG): $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(QUIPDUM) libamh.a
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) $(AMB12DUM) $(CUDADUM) $(QUIPDUM) $(LDFLAGS) libamh.a $(LIBS) 

#  no agressive optimizations for selected files to cut down on compile time
keywords.o: keywords.f
	${FC} ${NOOPT} ${SEARCH_PATH} ${INCLUDE_COMMAND} -c keywords.f

clean:
	rm -f ${GPROG} ${CPROG} $(AMHPROG) $(AMBPROG) $(OBJS) $(CUDAPROG) *.mod $(EXTRAS) *.lst *.o pref.dat prefx.msg nag.f90 FOR021.DAT porfuncs.f90 display_version.f90 mc.f bspt.f bspt_one_atom.f main.f tryexchange.f90 bhpt.f90 *.a tryexchange_independence.f90 $(AMB12DUM) $(SVN_ROOT)/AMBER12/*.mod
#	rm -f ${GPROG} ${CPROG} $(AMHPROG) $(AMBPROG) $(OBJS) $(CUDAPROG) $(QUIPPROG) *.mod $(EXTRAS) *.lst *.o pref.dat prefx.msg nag.f90 FOR021.DAT porfuncs.f90 display_version.f90 mc.f ptmc.f bspt.f bspt_one_atom.f main.f tryexchange.f90 bhpt.f90 *.a tryexchange_independence.f90 io1.f mylbfgs.f chirality.f90 potential.f keywords.f
	if test -d ../../BLAS ;  then cd ../../BLAS ; make clean ; fi
	if test -d ../../LAPACK ;  then cd ../../LAPACK ; make clean ; fi
	if test -d CHARMMinterface ;  then cd CHARMMinterface ; make clean ; fi
	if test -d CHARMMinterface35 ;  then cd CHARMMinterface35 ; make clean ; fi
	if test -d CHARMMboth ;  then cd CHARMMboth ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi
	if test -d $(SANDERSRC) ; then cd $(SANDERSRC) ; make clean ; fi
#	if test -d ../../QUIP ; then cd ../../QUIP ; export QUIP_ARCH=linux_x86_64_ifort_icc; make clean ; fi
	cd libmbpol; make clean; cd ../

cleanexe:
	rm -f $(CPROG) $(GPROG) $(AMHPROG) $(AMBPROG) $(CUDAPROG)
#	rm -f $(CPROG) $(GPROG) $(AMHPROG) $(AMBPROG) $(CUDAPROG) $(QUIPPROG)

timing:
	rm -f GMIN.2.0
	$(FC) $(FFLAGS) -p $(OBJS) $(EXTRAS) -o $(GPROG) $(LIBS)

feedback:
	$(FC) $(FFLAGS) -xprofile=use:gmin.profile *.f -o $(GPROG)

rebuild:
	make clean
	make

rebuildamh:
	make clean
	make AMHGMIN

libamber.a:
	export SRCDIR=$(CURDIR);cd ${SANDERSRC}; make lib1 FC="${FC}" FFLAGS="${FFLAGS}" \
	FREEFORMAT_FLAG="${FREEFORMAT_FLAG}" EXTRA_FLAGS="${EXTRA_FLAGS}"
libnab.a:
	cd $(NABSRC)/../ucpp-1.3/; make install; cd $(NABSRC)/../byacc/; make install; cd $(NABSRC); make install NABHOME="$(CURDIR)/../../NAB" LIBDIR="$(CURDIR)/../../NAB/lib";cp $(CURDIR)/../../NAB/lib/libnab.a $(CURDIR)/libnab.a
libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 
libmyblas.a: SAT-Ghost
	cd ../../BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
	cp libmyblas.a ../GMIN/source/
libmylapack.a: SAT-Ghost
	cd ../../LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
	cp libmylapack.a ../GMIN/source/

#libfullquip.a:
#	cd ../../QUIP; export QUIP_ARCH=linux_x86_64_ifort_icc; make all
#	mkdir ../../QUIP/QUIPOBJS
#	cd ../../QUIP; find . -name "*.o" -exec cp {} QUIPOBJS/ ";"
#	ar rv libfullquip.a ../../QUIP/QUIPOBJS/*.o
#	rm -r ../../QUIP/QUIPOBJS

ifeq (${CTYPE},C31)
libcharmm.a: SAT-Ghost
	echo ${chcommand}
	cd ${C31SRC}; ${chcommand}
	echo "Makefile> tail of" ${PREFDIR}/gnu.log
	tail -2 ${PREFDIR}/gnu.log
	cd CHARMMboth; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH} -I${PREFDIR}" PREFLX="${PREFLX}" \
	PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR="${FCMDIR}" CHSRC="${C31SRC}" SRCCH="${SRC31}"
endif
ifeq (${CTYPE},C35)
libcharmm.a: SAT-Ghost
	echo ${chcommand}
	cd ${C35SRC} ; ${chcommand}
	echo "Makefile> tail of" ${PREFDIR}/gnu.log
	tail -2 ${PREFDIR}/gnu.log
	cd CHARMMboth; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH} -I${PREFDIR}" PREFLX="${PREFLX}" \
	PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR="${FCMDIR}" CHSRC="${C35SRC}" SRCCH="${SRC35}"
endif

SAT-Ghost:

porfuncs.f90: porfuncs.sh
	./porfuncs.sh ${SWITCH} > porfuncs.f90

libmbpol.a:
	cd libmbpol; make; cp libmbpol.a ../; cd ../ 

# include ../../include/Makefile.gmin.display_version 

DISPLAY_VERSION=../../SCRIPTS/all/display_version.sh

dvopts=fflags "${FFLAGS}" prog GMIN fc_full_name "${FULL_COMPILER_NAME}" fc_exec "${FC}" make_opts "${MAKE_OPTS}" 

display_version.f90: ${DISPLAY_VERSION} 
	${DISPLAY_VERSION} ${dvopts} > display_version.f90

#read_cmd_args.f90: read_cmd_args.F90

# 
###################################### DEPENDENCIES ###################################### 
${GPROG}: libmyblas.a libmylapack.a 
${CPROG}: libmyblas.a libmylapack.a libcharmm.a
${AMHPROG}: commons.o libmyblas.a libmylapack.a libamh.a
${AMBPROG}: libmyblas.a libmylapack.a libamber.a porfuncs.o
${CUDAPROG}: libmyblas.a libmylapack.a
#${QUIPPROG}: libfullquip.a 
libcharmm.a: commons.o modcharmm.o modmxatms.o
libamber.a: commons.o modamber9.o porfuncs.o grouprotation.o
#${OBJS2}: ${OBJS1} 


Alglue.o:      commons.o
atomlists.o:   commons.o
atom_label_swaps.o: commons.o
atom_label_flips.o: commons.o
boxcentroid.o: commons.o	
binary_id_swaps.o: commons.o
BLJcluster.o:  commons.o
GLJ.o:  commons.o
BLJcluster_nocut.o: commons.o
BGupta.o:      commons.o
detsym.o:      commons.o
homoref.o:     commons.o
homoref_addons.o: commons.o
enperms.o:     commons.o
multiperm.o:   commons.o porfuncs.o
EAMLJ_sub.o:   commons.o
EAMal.o:       commons.o
Gupta.o:       commons.o
glj_yukawa.o:  commons.o
MSorig.o:      commons.o modconsts.o dist.o
MSorigc.o:     commons.o modconsts.o dist.o
MStrans.97.o:  commons.o modconsts_trans_97.o dist.o
Mgglue.o:      commons.o
PachecoC60.o:  commons.o
Pbglue.o:      commons.o
SW.o:          commons.o
SiSW.o:        commons.o
Zetterling.o:  commons.o
amber.o:       commons.o modamb.o
capbin.o:      commons.o
capsid.o:      commons.o
centre.o:      commons.o genrigid.o
genrigid.o:    commons.o rotations.o vec3.o ${AMB12DUM}
Gthomson.o:    commons.o rotations.o vec3.o
centrecom.o:   commons.o precision.o
cgmin.o:       commons.o
checkd.o:      commons.o
chiro.o:       commons.o rotations.o qmod.o vec3.o
compress.o:    commons.o
countatoms.o:  modamber9.o
overlap.o:     commons.o porfuncs.o
dbpg.o:        commons.o
dbptd.o:       commons.o
dmblmorse.o:   commons.o
dmblpy.o:      commons.o
ddfpmin.o:     commons.o
dftb.o:        commons.o
dftbC.o:        commons.o
dzugutov.o:    commons.o
evstep.o:      commons.o
f1dim.o:       commons.o modf1com.o
fd.o:          commons.o
fdm.o:         commons.o
finalio.o:     commons.o modamb.o qmod.o modcharmm.o modamber9.o BGupta.o modperm.o 
finalq.o:      commons.o qmod.o
frausi.o:      dist.o
gbcalamitic.o: commons.o
gbdiscotic.o:  commons.o
gbddp.o:       commons.o
gem.o:         commons.o
grnd.o:        commons.o
harmonic.o:    commons.o
hmat1n.o:      commons.o
hmatd_.o:      commons.o
hsmove.o:      commons.o
io1.o:         commons.o modamb.o modperm.o qmod.o modcharmm.o porfuncs.o BGupta.o sandbox.o py.o
keywords.o:    commons.o modamb.o modcharmm.o porfuncs.o ga_modules.o BGupta.o glj_yukawa.o chiro.o sandbox.o restraindistance.o polirmod.o swmod.o chirality.o mbpol.o parse_pot_params.o ${AMB12DUM} twist.o rotamer_move.o
linmin.o:      commons.o modf1com.o
linrod.o:      commons.o
ljcoulomb.o:   commons.o
lj.o:          commons.o 
ljgh.o:        commons.o modhess.o
ljpbin.o:      commons.o
ljpshift.o:    commons.o neighbor_list.o neighbor_list_bin.o cell_lists.o cell_lists_binary.o neighbor_list_moveone.o
#gmin_quip_wrapper.o: commons.o
lwotp.o:       commons.o rbperm.o
newpah.o:      commons.o
main.f: main.F
main.o:        commons.o modf1com.o countatoms.o modperm.o qmod.o modamb.o modmxatms.o porfuncs.o initialization.o homoref_addons.o ${AMB12DUM} twist.o
mcruns.o:      commons.o
mf.o:          commons.o
mnbrak.o:      commons.o
modcudadummy.o: commons.o
modcudalbfgs.o: commons.o
moves.o:       vec3.o sanity.o
dummyoptim.o:  commons.o
mindist.o:      commons.o porfuncs.o
morse.o:       commons.o
morsegh.o:       commons.o modhess.o
mylbfgs.o:     commons.o modamb.o porfuncs.o ${CUDADUM} rad.o
mymylbfgs.o:     commons.o modamb.o porfuncs.o
mstbin.o:      commons.o
multstock.o:   commons.o
newtip.o:      commons.o
newcapsid.o:   commons.o
patchy.o:      commons.o
asaoos.o:      commons.o
odesd.o:       commons.o mododesd.o
olami.o:       commons.o
otp.o:         commons.o
BLN.o:         commons.o
pah.o:         commons.o
multpaha.o:    commons.o qmod.o
pahw99.o:      commons.o
paha.o:        commons.o
pap.o:         commons.o
papbin.o:      commons.o qmod.o
papjanus.o:    commons.o qmod.o
percsphere.o:  commons.o
parse_pot_params.o: commons.o pot_params.o
ptstst.o:      commons.o
potential.o:   commons.o modperm.o qmod.o modcharmm.o restraindistance.o polirmod.o swmod.o mbpol.o ${AMB12DUM} ${CUDADUM} rad.o twist.o precision.o porfuncs.o
polirmod.o:    commons.o ps_intra_polir_b.o
ptgrp.o:       commons.o
pgsym_mod.o:   commons.o
pgsym.o:       commons.o pgsym_mod.o
py.o:          commons.o vec3.o qmod.o rotations.o
phi4model.o:	commons.o
quench.o:      commons.o qmod.o porfuncs.o 
QALCSearch.o:  commons.o
QALCS_surface.o:  commons.o
QALCS_symmetry.o:  commons.o QALCS_symmetry_mod.o
QALCS_symmetry_mod.o: commons.o pgsym_mod.o
rad.o:         commons.o precision.o
ratio.o:       commons.o
rotations.o:   commons.o
rdpot.o:       commons.o
sandbox.o:     commons.o rotations.o qmod.o vec3.o
saveit.o:      commons.o qmod.o
sc.o:          commons.o
mie_field.o:   commons.o
MGupta.o:      commons.o pot_params.o 
MLJ.o:         commons.o pot_params.o 
MSC.o:         commons.o pot_params.o 
seed.o:        commons.o
silane.o:      commons.o
stockaa.o:     commons.o
strand.o:      commons.o
stress.o:      commons.o
supermc.o:     commons.o
swmod.o:       commons.o
symmetry.o:    commons.o porfuncs.o
symmetrycsm.o:    commons.o porfuncs.o
taboo.o:       commons.o
takestep.o:    commons.o rotations.o vec3.o genrigid.o
tetrahedra.o:  commons.o
tip.o:         commons.o
tosifumi.o:    commons.o modhess.o
welch.o:       commons.o
waterpdc.o:    commons.o
waterpkz.o:    commons.o
zwischen.o:    commons.o modf1com.o precision.o
stock.o:       commons.o
sticky.o:      commons.o
tether.o:      tetherfunc.o
mycpu_time.o:  commons.o
gauss.o:       commons.o precision.o
mc.o:          qmod.o modcharmm.o porfuncs.o mc.f mc.F commons.o grouprotation.o operations.o utils.o rotations.o vec3.o bhpt.o hbondmatrix.o Gthomson.o  chirality.o mbpol.o ${AMB12DUM} rotamer_move.o
mc.f: mc.F
reservoir_utils.o: commons.o random_normal.o
biased_torsional_sampling.o: commons.o random_normal.o
bspt_utils.o:      commons.o grouprotation.o mymylbfgs.o dprand.o perc.o porfuncs.o
bspt.f: bspt.F
bspt.o:          qmod.o modcharmm.o porfuncs.o tetherfunc.o bspt.f bspt.F utils.o overlap.o tryexchange.o bspt_utils.o
bspt_one_atom.f: bspt_one_atom.F
bspt_one_atom.o:          qmod.o modcharmm.o porfuncs.o tetherfunc.o bspt_one_atom.f bspt_one_atom.F utils.o overlap.o ljpshift.o tryexchange.o bspt_utils.o
GMINdump.o: commons.o qmod.o porfuncs.o output.o
output.o: commons.o file_manager.o
quad.o:        commons.o
enumerate.o:  commons.o
convert.o:  porfuncs.o
getorbits.o: commons.o
mysd.o: commons.o
mW_wrapper.o: commons.o
qdtest.o: commons.o
MB.o: commons.o 
sort2.o: qmod.o commons.o
model1.o: commons.o
FS.o: commons.o
Ackland_wrapper.o: commons.o
Ackland_metals.o:      ackland_sma.h ackland_mishin_cu.h ackland_mendelev_fe.h commons.o
DF1.o: commons.o
CSMinit.o: commons.o
minpermdist.o: commons.o porfuncs.o rbperm.o
grouprotation.o: commons.o moves.o
hbondmatrix.o : porfuncs.o commons.o ${AMB12DUM}
dampedmove.o: commons.o
finalio.o: gay-berne.o sandbox.o hbondmatrix.o py.o ${AMB12DUM}
chaperonin.o: commons.o
bulkmindist.o: commons.o
read_cmd_args.o: porfuncs.o display_version.o read_cmd_args.f90 commons.o
ga_main.o: ga_modules.o commons.o
ga_select: ga_modules.o commons.o
ga_bln.o: ga_modules.o commons.o
ga_cluster.o: ga_modules.o commons.o
tryexchange.f90: tryexchange.F90
tryexchange.o: commons.o tryexchange.f90 tryexchange.F90 tryexchange_independence.o reservoir_utils.o
tryexchange_independence.f90: tryexchange_independence.F90
tryexchange_independence.o: commons.o tryexchange_independence.f90
minpermdistrbcom.o: commons.o
rbperm.o: commons.o minpermdistrbcom.o
rotations.o: vec3.o
msaveit.o: commons.o modamber9.o porfuncs.o rotations.o
bhpt.o: bhpt.f90 bhpt.F90 commons.o dprand.o
bhpt.f90: bhpt.F90
initialization.o: commons.o overlap.o polirmod.o BGupta.o homoref_addons.o qmod.o
py.o: commons.o rotations.o vec3.o qmod.o gay-berne.o
Energy_1d_APBC.o: commons.o
Energy_2d_PBC.o: commons.o
Energy_3d_APBC.o: commons.o
Energy_3d_PBC.o: commons.o
orient.o: commons.o
cell_lists.o: cell_lists.f90
cell_lists_binary.o: cell_lists.f90
neighbor_list.o: neighbor_list.f90 cell_lists.o
neighbor_list_bin.o: cell_lists.o
neighbor_list_moveone.o: 
selectmove.o: commons.o
dihedralrotation.o: commons.o biased_torsional_sampling.o
dipolarmorse.o: commons.o
chirality.o: grafpack.o ${AMB12DUM}
mbpol.o: commons.o libmbpol.a
inertia.o: commons.o 
dlnsrch.o: precision.o
reseed.o: precision.o
trans.o: precision.o
mlowest.o: commons.o 
twist.o: precision.o dihedral.o
dihedral.o: precision.o
