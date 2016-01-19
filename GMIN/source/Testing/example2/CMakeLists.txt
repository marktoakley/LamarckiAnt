message("Adding test for GMIN")

#this calls a macro which sets up a test called example2 which runs
#the program GMIN in directory run1.  It then runs a script which tests
#if the energy in GMIN_out agrees with the correct value
gmin_quench_test(example2 GMIN run1 -405.8736755850)
