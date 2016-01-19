      Module amhglobals
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old arraysize.inc include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
        save

c     maxsiz is the maximum number of amino acid residues

         integer maxsiz
         parameter (maxsiz=80)

c     maxcnt is the maximum number of (i,j)
c     interactions

        integer maxcnt
        parameter(maxcnt=(maxsiz**2-maxsiz)/2)
c       parameter(maxcnt=31250)
c       parameter(maxcnt=7022)
        integer max_well
        parameter(max_well=10)

      integer max_letters
      parameter (max_letters=21)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old head include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer maxmem
      parameter(maxmem=40)

      integer maxpro
      parameter(maxpro=1)

c     maxres is the maximum number of amino acid residues
c     for a database protein

      integer maxres
      parameter(maxres=650)

c     maxcrd is the maximum number of coordinate types;

      integer maxcrd
      parameter(maxcrd=3)

c     maxtab is the maximum number of tables, i.e.,
c     the number of pairwise ixns;

      integer maxtab
      parameter(maxtab=4)

c     multiple seq  + partial constraining
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical cyclic,ss_a,ss_b
      integer num_foldon_a,num_foldon_b,const_mode, totalssnorm
      integer foldstrt_min_a(maxres),foldstrt_max_a(maxres)
      integer foldstrt_min_b(maxres),foldstrt_max_b(maxres)
      integer ss_pattern_a(maxres), ss_pattern_b(maxres)
!         integer ss_alpha(maxres), ss_beta(maxres)
!         integer sa_alpha(maxres), sa_beta(maxres)
 
      double precision ss_dist(500,500,maxtab)

c     switch for averaging force table over sequences

      logical ave_seq,targ_cons, mem_cons,ave_seq_hb
      logical ssweight,ave_seq_amw, ave_seq_amc
      integer maxseq,numseq,numseq_hb,numseq_amw
      integer numseq_amc
      parameter (maxseq=2)
      
      integer tgsequences(maxres,maxseq)
      integer tgsequences_hb(maxres,maxseq)
      integer tgsequences_amw(maxres,maxseq)
      integer tgsequences_amc(maxres,maxseq)

      double precision aps(maxres,maxres)

       character homlfl(maxseq)*5
cccc  MCP

      DOUBLE PRECISION x_int(maxres*3*3)
      DOUBLE PRECISION x_mcp(maxres*3*3) 
      double precision E_temp, E_harm_springs

cccc   go contact

      logical go_con
      double precision    go_con_dist

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old additions include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical alpha_c_of_n,ab_c_of_n_old,ab_c_of_n_new,quench
      logical beta_c_of_n,dist_cut
      integer nquench

c     oxscl is scaling constant for oxygen potential

      double precision oxscl

c     hbscl is the scale factor for hbonds

        double precision hbscl(17)
        double precision sigma_h,sigma_NO,ho_zero,NO_zero
        double precision :: anti_NHB(20,20,2),anti_HB(20,20,2),
     *                  para_HB(20,20,2),para_one(20),anti_one(20)

c     lamdaR is a scale factor for incorrect 
c     memories

        double precision lambdaR

c     chrlscl is the scaling parameter for chiral forces

      double precision chrlscl

c     scaling parameters for ramachondrian plots

      double precision ramascl

c     scaling parameters for oxygen excluded volume

      double precision oexcldscl

c     C-C excluded volume distance array

      double precision ccev_dist(1:maxcnt,1:maxtab),ooev_dist(1:maxcnt)

c       flag for sec-setruc dependendent gamma is ss class

        logical four_plus
c     other flags controlling encoding used
c     i_ignore_chain_dirn says encoding doesn't depend on chain direction 
c     mismatch says encoding doesn't depend on what memory protein pair
c     is if the identities are not identical to the target pair (ie they are mismatched)

        logical i_ignore_chain_dirn,mismatch

c     flag and scle factor for biasing the  funnel
c temp old 
!      logical i_Qbias,i_bias_native
!      integer n_Qbias,i_Q_format
!      double precision Qbiaspoly(1:100),
!      double precision          width_Qexp,Q0,Q_weight,Q_clip
!      integer n_Q_anneal


      integer n_divs_max
      parameter (n_divs_max=1000)
      logical ibiasgauss,ibiaspoly,i_Qbias_a,i_bias_native_a
      logical con_local_a, con_local_b
      integer con_local_a_cut, con_local_b_cut
      integer nbiaspoly,n_Qbias_a,i_Q_format_a
      double precision bias_weight,bias_av,bias_var,bias_prefactor,
     *           biaspoly(1:100),Qbiaspoly_a(1:100),
     *           width_Qexp_a,Q0_a,Q_weight_a,Q_clip_a,
     *           Q0_safe_a, Q0_safe_b,
     *           Q0_inc_a,Q0_inc_b,Qvalue_a,Qvalue_b
      integer n_Q_anneal_a,n_Q_anneal_b

       double precision del_r_a(1:maxsiz),
     *           Q_ij_a(0:n_divs_max,maxsiz),
     *          dQ_dr_ij_a(1:n_divs_max,1:maxsiz)

      logical i_Qbias_b,i_bias_native_b
      integer n_Qbias_b,i_Q_format_b
      double precision   Qbiaspoly_b(1:100)
      double precision        width_Qexp_b, Q0_b,Q_weight_b,Q_clip_b

      integer seglist_a(maxres), seglist_b(maxres),
     *        numconst_a,i_ixn_Qbias_b(maxres,maxres),
     *        numconst_b,i_ixn_Qbias_a(maxres,maxres)

       double precision del_r_b(1:maxsiz),
     *           Q_ij_b(0:n_divs_max,maxsiz),
     *          dQ_dr_ij_b(1:n_divs_max,1:maxsiz)

c     flag and parameters for |i-j| dependent 'contact order' interaction
      logical i_contact_order
      integer n_contact_order_terms
      integer, parameter:: max_contact_order_terms=5
      integer i_contact_order_min(max_contact_order_terms),
     *        i_contact_order_max(max_contact_order_terms)
      double precision r_min_contact_order(max_contact_order_terms),
     *     r_max_contact_order(max_contact_order_terms)
      double precision gamma_contact_order(max_contact_order_terms,3)
c     GAP, 03/23/04, The minimum separation + 1 for the interactions in the contact class
      integer, parameter:: CUTOFF_CONT_LOW=13

c     flag and parameters for radius of gyration force

      logical i_rg_corey
      double precision rg_bounds(2)
      logical i_Rg_bias
      integer n_Rg_bias
      double precision Rg_biaspoly(1:100)
      double precision rg_shift
      double precision rg_scl

      logical i_rg_garyk, i_rg_first
      double precision D_rg, T_rg, delR_rg, M_rg, kappa_rg

      double precision targ_dist(maxcnt,maxtab)
!      double precision dQ_dr_ij(1:n_divs_max,1:maxsiz)
!    e double precision del_r(1:maxsiz), Q_ij(0:n_divs_max,maxsiz)

      double precision ngomb
      logical igomb

c iscltab is switch for scltab subroutine
      logical iscltab
c iclassweight allows (rather sneaky) scaling of gammas by class
      logical iclassweight

c some variables for non-additive contact potential

      logical :: i_non_add_contact
      integer, parameter:: ngamma_non_add=32
      integer class_of_res_2(maxsiz),sort_non_add(2,2,2,2,2),
     *   class_2(20)
      double precision :: gamma_non_add(ngamma_non_add)

c variables for Parallel/Anti-Parallel contact potential

      logical:: i_con_P_AP
      integer:: i_diff_P_AP,i_atom_P_AP
       double precision :: weight_P_AP(3),r_cut_P_AP

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old potent include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     parameter file for potential

c     nummem is the number of protein memories 
c     used in the potential

      integer nummem

c     numpro is the number of ensemble proteins 
 
      integer numpro
c     nmres is the actual protein size used

      integer nmres

c     i_rama is the flag for rama potential
c     (0=off 1=on (new) 2=on (old)

       integer i_rama

c     hbond is the hbond flag for hbonds

        logical hbond,i_hbond_eastwood,i_hbond_czong

c     oxexcldv is the flag for oxygen-oxygen excluded volume

        logical oxexcldv
        double precision o_exvmin,o_exvminS(1:4)   !sep where O-O ev comes on

c     numcrd is the number of coordinate types

      integer numcrd,numtab

c     charge switches

      logical i_alt_prox,i_3rd_is_contact,gly_con
      integer n_letters,min_seq_sep,num_well,n_letters_con
      double precision alt_prox_cut

c     maxr is the grid resolution of the force lookup table

      integer maxr,maxs
      parameter(maxr=400)

c     srplr is  rel weights of short-med-long range energy if iclassweight=true

      double precision srplr(3)

c     eres is the target energy/residue

      double precision eres

c     minmr and maxmr partition the protein
c     into short- med- and long-range ixns
c     the medium range class goes from minmr--maxmr inclusive

      integer minmr,maxmr

c     delta is used in determining the range
c     of the upper and lower bounds for the memory
c     distances (the square-well half width)

      double precision delta,delte

c     rcut is the maximum ixn distance
c     sr cut is division between sr/lr ints (in space) 

      double precision rcut,rcutAMH,srcut

c     if iexcld, then turn on excluded-volume penalty

      logical iexcld

      logical iexcld_gamma,iexcld_beta
      integer i_type_ev_gamma
      double precision pexcld_gamma
      double precision exvmin_gamma(maxsiz,maxsiz)
      double precision c_of_m_dist(maxsiz)

c     rexcld is the minimum allowed distance 
c     between residues
c     pexcld is the penalty associated w/ violating the self-
c     avoiding cutoff, rexcld

      double precision pexcld

c     exvmin is where C-C excluded vol starts (in Angstroms)
c     and exvminS is where C-C excluded vol starts (in Angstroms) for |i-j|<5

      double precision exvmin,exvminS,exvminS_beta(1:4)

c     if welscl, then scale well-widths by |i-j|;
c     otherwise, scale by rij

      logical welscl

c     maxsec is the maximum number of secondary
c     structures of any given type for a given protein

      integer maxsec
      parameter(maxsec=30)

c     maxmov is the maximum number of quench structures

      integer maxmov
      parameter(maxmov=1)

c     timstp is the timestep for the verlet algorithm

      double precision timstp

c     tolshk is the maximum deviation allowed for each constraint

      double precision tolshk

c Protein name  protnm
      character protnm(0:maxmem)*6

c     flags for random interactions

      logical ran_force,ran_file,allow_neg
      double precision  r_ran,s_ran
      integer ran_min_seq_dist

      integer ilong(maxcnt,2,maxtab),
     *        numlng(0:maxsiz,maxtab),ires(maxsiz),
     *        sa(maxres),
     *        tarpre(maxres,maxtab),mempre(maxres,maxtab),
     *        jres(maxres),
     *        crdixn(maxtab,2)

      integer ixn_from_site(maxsiz,maxsiz,maxtab)

 
      double precision vpotnt(0:maxr+1,maxcnt,maxtab),
     *     amh_gse,
     *     rinc(maxcnt,maxtab),
     *     rincinv(maxcnt,maxtab),
     *     rincsq(maxcnt,maxtab),
     *     forse(0:maxr+1,maxcnt,maxtab),
     *     deltz(maxcnt,maxtab),shydro(maxres,2,maxtab),
     *     ywork(maxres,3,maxcrd),
     *     hydseq(maxres,2,maxtab),
     *     qchrg(21,21,20,20,max_well+2),
     *     qchrg2(21,21,20,20,2,2)

c++ Johan
!          double precision altgamma(20,20,max_well)
c-- Johan
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old utility include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        character*10 name_a(100000)

c     include file for seeds,
c     input_amh-output files, do loop

c     --- work space ---

      integer iwork(maxcnt)

      double precision work1(maxcnt),work2(maxcnt),work3(maxcnt),
     *     work4(maxcnt),work6(maxcnt),work8(maxcnt)

c     --- i/o tapes ---

      integer imemri,iprolst,iprolstscl,iran,oarchv,oconv,oecorr,imat,
     *        iss_bias,ihbond,imoviein,igamma,itargseqs,ihomol,
     *        omovi,omoviseg,SO,input_amh,ievgamma,ievbeta,icon,imem_cons,
     *        ohdrgn,orama,ooxy,ochiral,oamh,orep,
     *        oamhsr,oamhlr,oran,occev,ooev,oamhmr,
     *        obias_Rg,oKE,ibiasfile,ohdrgn_s,
     *        ohdrgn_m,ohdrgn_l,ohdrgn_seq,ononadd,ocon_P_AP,pdb,
     *        oPE_plus_KE,oPE_no_bias,oPE_with_bias,oPE_backbone,oPE_backbone_norama,
     *        oobiassega,oobiassegb,oharmspring,
     *        ovgaspot,ovscltab,iphi,icontacts,
     *        iss_struct,itarg_seq

c     --- seeds ---

c     iseed_amh(4) are the seeds used for the random number generators

      integer iseed_amh(4)

c     --- do loop indicies ---

      integer i502,i503,i505,i507,
     *        i508,i509,i511,i512,i513,i514,i515,
     *        i516,i518,i521,i523,i525,i526,i527,i528,
     *        i540
c     --- implied do loop indicies ---

      integer i1,i2

c     --- diagnostic flags ---

      logical idigns

c     --- amino acid abbreviations --

      character aminoa(21)*3

      double precision hydscl(0:21,2),hrdrad(0:21,2),eqdist(5)

c     set up mapping from sequence id to amino acid residue

      data (aminoa(i1),i1=1,21)/'ALA','ARG','ASN','ASP',
     *                          'CYS','GLN','GLU','GLY',
     *                          'HIS','ILE','LEU','LYS',
     *                          'MET','PHE','PRO','SER',
     *                          'THR','TRP','TYR','VAL',
     *                          'INI'/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             old param.s include file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     incmov is the increment for the number of 
c     structures to be saved fixed T

      integer incmov

c     nmdif is greater than mxtemp/incmov

      integer nmdif
      parameter(nmdif=2001)

c     t1-t5 are the cutoff points in the temperature-annealing
c     schedule

      double precision t1,t2,t3,t4,t5
 
c     it1-it5 are then number of T grid points allocated to each
c     of the intervals [t1,t2], [t2,t3],...

      integer it1,it2,it3,it4,it5

c     mxtemp is the number of temperature grid points

      integer mxtemp,nmtemp
      parameter(mxtemp=12001)

      integer nmstep

c     set hierarchy parameters:

c     mxlevl is the maximum number of hierarchy levels
c     mxsct is the maximum number of divisions per level       

      integer mxlevl,mxsct
      parameter(mxlevl=5,mxsct=5)

c     lenfun is the maximum protein segment length 
c     considered per level

      integer lenfun
      parameter(lenfun=155)

c     if ictemp, then the minimization is performed at a 
c     constant temperature, ctemp

      logical ictemp
      double precision ctemp

!     if i_V_test print out test potentials between 2 sites
!     in gaspot, and again in scltab

      logical i_V_test
      integer test_site(2)

c     if iresrt, then read in previously generated
c     proteins used for restarting program or analysis

      logical iresrt



c     if movanal, then analyse movie structures only

      logical movanal

      logical known

c     i_Etim=T if want to output energies every i_Etim_step timesteps
      logical i_Etim
      integer i_Etim_step


      integer itgrd(4)
 
      double precision temtur(mxtemp),temtur_quench(maxmov,maxmov),temgrd(5),
     *     bondln(maxsiz,maxcrd),target(maxsiz,3,maxcrd),
     *     prcord(maxsiz,3,maxpro,maxcrd),
     *     zrcord(maxsiz,3,maxpro,maxcrd),
     *     avep(maxpro,2,50),
     *     quench_crd(maxsiz,3,maxpro,maxcrd,maxmov)

      double precision avepp(maxpro,2,50)

! distance cut-offs for contact potential
! r(1) is closest range well
      double precision r_min(1:max_well),r_max(1:max_well)

! variables for replica force

      logical i_rep,i_rep_lambda_uniform
      double precision rep_phi_exp(maxsiz),rep_lambda(maxsiz),
     *   rep_cut_off,rep_tol
      integer n_rep_con(maxsiz),rep_con_2_res(maxsiz,maxsiz)

      end module amhglobals
