include ./make.inc
# SAT: This version of BLAS was adapted from LAPACK distribution (version 3.0)
#
#######################################################################
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       blas.a
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------

FMODULES_SOURCES = amhglobals.f altpot_interfaces.f90  amh_interfaces.f90  globals_alt.f90 
FMODULES_OBJECTS = altpot_interfaces.o  amh_interfaces.o  globals_alt.o amhglobals.o
FMODULES = altpot_interfaces.mod altpot_interfaces.mod  amhglobals.mod globals_alt.mod 
F_RAW_SOURCES = $(wildcard *.f)
F_SOURCES = $(filter-out $(FMODULES_SOURCES),$(F_RAW_SOURCES))
F90_RAW_SOURCES = $(wildcard *.f90)
F90_SOURCES = $(filter-out $(FMODULES_SOURCES),$(F90_RAW_SOURCES))
F_OBJECTS = $(F_SOURCES:%.f=%.o)
F90_OBJECTS = $(F90_SOURCES:%.f90=%.o)

all: single
single: $(AMHLIB)

${AMHLIB}: $(FMODULES_OBJECTS) $(F_OBJECTS) $(F90_OBJECTS)
	$(ARCH) $(ARCHFLAGS) $(AMHLIB) $(F_OBJECTS) $(F90_OBJECTS) $^

$(F_OBJECTS),$(F90_OBJECTS): $(FMODULES)

$(FMODULES_OBJECTS): $(FMODULES) 

amhglobals.mod: amhglobals.f
	$(FC) $(FFLAGS)  -c $<  

amh_interfaces.mod: amh_interfaces.f90 amhglobals.mod
	$(FC) $(FFLAGS)  -c $< 

globals_alt.mod: globals_alt.f90 amhglobals.mod
	$(FC) $(FFLAGS)  -c $<  

altpot_interfaces.mod: altpot_interfaces.f90 amhglobals.mod globals_alt.mod
	$(FC) $(FFLAGS)  -c $< 

%.o: %.f
	$(FC) $(FFLAGS) -c $< 

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -f *.o ${AMHLIB} ifc?????? *.lst *mod

################### DEPENDENCIES ######################################

#walesamh_finalio.o: amhglobals.o 
additive_ev.o: amhglobals.o
dssp_hdrgn_eastwood.o:   amhglobals.o
dssp_hforce.o:  amhglobals.o
force.o:       amhglobals.o  amh_interfaces.o  globals_alt.o
gaspot.o:      amhglobals.o
gengrd.o:      amhglobals.o
gentab.o:      amhglobals.o  amh_interfaces.o
getmem.o:      amhglobals.o
gomb.o:        amhglobals.o
force.o:       amhglobals.o
initil.o:      amhglobals.o
initv.o:       amhglobals.o
intstr.o:      amhglobals.o
lookup.o:      amhglobals.o  amh_interfaces.o
walesamh_initial.o:     amhglobals.o amh_interfaces.o ../commons.o
oxy.o:         amhglobals.o  amh_interfaces.o
potind.o:      amhglobals.o
Q_bias_seg_a.o:      amhglobals.o
Q_bias_seg_b.o:      amhglobals.o
qchrgmk.o:      amhglobals.o
rama.o:        amhglobals.o
restrt.o:      amhglobals.o
rndcol.o:      amhglobals.o
savstr.o:      amhglobals.o
scltab.o:      amhglobals.o
walesamh_interface.o: amhglobals.o 
zero_amh.o:        amhglobals.o
altpot_io.o:      globals_alt.o   altpot_interfaces.o  amhglobals.o
altpot_onebody.o:      altpot_interfaces.o  amhglobals.o  globals_alt.o
altpot_pairpotential.o:      amhglobals.o globals_alt.o
altpot_subroutines.o:      amhglobals.o  globals_alt.o
altpot_utilities.o:      amhglobals.o  globals_alt.o
contact_P_AP.o:      amhglobals.o
ev_gamma.o:      amhglobals.o
ev_set_up.o:     amhglobals.o
E_write.o:       amhglobals.o
hdrgn.o:         amhglobals.o
harm_spring:	 amhglobals.o
move.o:          amhglobals.o  amh_interfaces.o  globals_alt.o
non_add_contact.o:       amhglobals.o
num_to_char.o:   amhglobals.o
rep_bias.o:      amhglobals.o
rep_contact.o:   amhglobals.o
Rg_bias.o:       amhglobals.o
