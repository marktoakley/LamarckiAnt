if(WITH_AMBER9)
    message("Adding binary test for AMBGMIN")

    #this calls a macro which sets up a test called example2 which runs
    #the program GMIN in directory run1.  It then runs a script which tests
    #if the energy in GMIN_out agrees with the correct value

#    gmin_quench_test(AMBGMINenery AMBGMIN run1 -1.9089938038)
    gmin_quench_test(AMBGMINenery AMBGMIN run1 -2.9089938038)
endif(WITH_AMBER9)
