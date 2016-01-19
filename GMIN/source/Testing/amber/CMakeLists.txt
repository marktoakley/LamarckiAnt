if(WITH_AMBER9)
   message("Adding AMBER test for AMBGMIN")

   #this calls a macro which sets up a test called example2 which runs
   #the program GMIN in directory run1.  It then runs a script which tests
   #if the energy in GMIN_out agrees with the correct value
   gmin_quench_test(amber_singlequench AMBGMIN singlequench -199.1192693581)
endif(WITH_AMBER9)
