      module globals_alt

        implicit none

        save

        ! PARAMETERS TO BE SET BEFORE COMPILATION:
        ! n_altgamtypes = # of different gammas for a particular well
        ! currently this can only be set to 2 
        ! max_well_alt = max wells for alternative gamma potential
        ! Alimits set the densitylimits for defining the wells of the 1-body potential

        integer max_well_alt,n_altgamtypes,max_letters,debugflag
        integer outfile1_alt(10),n_OB_wells,timings_file_alt
        parameter(n_altgamtypes=2)
        parameter(max_well_alt=2)
        parameter(max_letters=20)
        parameter(debugflag=0)
        parameter(n_OB_wells=3)


        character, dimension(20) ::  aminoacids*3 = (/"ALA",  "ARG",   "ASN",   "ASP",   "CYS", &
             "GLN",   "GLU",   "GLY",   "HIS",   "ILE",   "LEU", &
             "LYS",   "MET",   "PHE",   "PRO",   "SER",   "THR", &
             "TRP",   "TYR",   "VAL" /) 
         double precision, dimension(20) :: hp_scale = (/ 0.45D0,  0.08D0,  0.15D0,  0.12D0, &      ! residue hydrophobicity scale
             0.40D0,  0.11D0,  0.11D0,  0.35D0,  0.16D0,  1.00D0,  0.61D0,  0.00D0,  0.18D0, &
             0.41D0,  0.23D0,  0.25D0,  0.23D0, 0.39D0,  0.25D0,  0.70D0 /)


        ! PARAMETERS TO BE READ IN FROM INPUTFILE:
        ! altpotflag decides wether or not to use alternative anisotropic potential
        ! (which is manybody and therefore probably a bit slowing down)

        integer altpotflag(max_well_alt)                             ! ALTPOT
        double precision kappa_alt,treshold_alt,kappa_well
        integer onebodyflag,onebody_type                          ! ONEBODY
        double precision kappa_OB
        logical debug_numer_forces

        ! VARIABLES TO BE SET DURING EXECUTION:

        integer output_stepsize_alt , output_step_alt ,count_alt                 ! ALTPOT
        double precision altgamma(max_letters,max_letters,n_altgamtypes,max_well_alt),T_alt
        double precision onebody_gamma(max_letters,n_OB_wells)     ! ONEBODY
        double precision Alimits_OB(2,n_OB_wells) 
        logical do_send_output


        !  DIAGNOSTIC VARIABLES (NOT NECESSARY FOR EXECUTION)
        !  OBSERVE!  n_OB_wells must be 3 currently !!!
        double precision OB_density(max_letters,n_OB_wells),OB_dns_count
         double precision, dimension(20) :: accumulated_time=0.0D0


      end module globals_alt
