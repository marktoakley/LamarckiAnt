

c     --------------------- initil ----------------------

      subroutine initil

c     INITIL open required files, read in input_amh 
c            parameters and print header

c     ---------------------------------------------------

      use amhglobals,  only:SO,iprolst,iprolstscl,input_amh,oamh,oamhsr,
     * oamhlr,nummem,maxmem,numpro,maxpro,nmres,maxsiz,numcrd,
     * maxcrd,numtab,maxtab,crdixn,i505,hrdrad,igomb,
     * ngomb,n_letters,max_letters,i_alt_prox,
     * alt_prox_cut,iexcld,allow_neg,lambdaR,eres,minmr,maxmr,maxs,
     * maxr,i_rama,hbond,hbscl,chrlscl,delta,delte,
     * ibiasgauss,bias_weight,i_hbond_eastwood,i_hbond_czong,
     * bias_av,bias_var,bias_prefactor,ibiaspoly,nbiaspoly,
     * biaspoly,i_Qbias_a,n_Qbias_a,Qbiaspoly_a,nmstep,iseed_amh,
     * it1,it2,it3,
     * it4,it5,ictemp,ctemp,nmtemp,mxtemp,incmov,nmdif,iresrt,
     * iscltab,movanal,itgrd,i1,i2,i503,hydscl,
     * aminoa,idigns,oarchv,omovi,
     *  orep,imemri,ires,
     *  oPE_no_bias,oPE_with_bias,oPE_plus_KE,oPE_backbone,oPE_backbone_norama,
     * ohdrgn,orama,ooxy,ochiral,oKE,occev,ooev,
     * obias_Rg,oamhmr,ohdrgn_s,ohdrgn_m,ohdrgn_l,ohdrgn_seq,
     * ononadd,ocon_P_AP,oobiassega,oobiassegb,oharmspring,
     * omoviseg,pexcld,oxexcldv,o_exvmin,o_exvminS,
     * oexcldscl,ran_force,ran_file,
     * r_ran,s_ran,ran_min_seq_dist,srplr,oxscl,ramascl,timstp,
     * tolshk,rcut,rcutAMH,srcut,welscl,width_Qexp_a,t1,t2,
     *  t3,t4,t5,temgrd,qchrg,qchrg2,i_Rg_bias,n_Rg_bias,
     *  i_rg_garyk,i_rg_first, D_rg, T_rg, delR_rg, M_rg, kappa_rg,
     *  rg_shift,rg_scl,Rg_biaspoly,i_3rd_is_contact,
     * min_seq_sep,iclassweight,exvmin,exvminS,exvminS_beta,
     *  r_min,r_max,gly_con,width_Qexp_b,
     * i_V_test,test_site,rg_bounds,i_rg_corey,
     * alpha_c_of_n,ab_c_of_n_old,ab_c_of_n_new,num_well,sigma_h,
     * sigma_NO,ho_zero,maxseq,homlfl,
     * ave_seq,ave_seq_hb,ave_seq_amw,ave_seq_amc,
     * numseq,numseq_hb,numseq_amw,numseq_amc,
     * NO_zero,four_plus,i_non_add_contact,i_Etim_step,i_Etim,quench,
     * known,pexcld_gamma,iexcld_gamma,i_type_ev_gamma,iexcld_beta,
     * i_contact_order,i_contact_order_min,i_contact_order_max,
     * gamma_contact_order,r_min_contact_order,r_max_contact_order,
     * n_contact_order_terms,max_contact_order_terms,
     * i_Qbias_b,n_Qbias_b,Qbiaspoly_b,
     * num_foldon_a,foldstrt_min_a, foldstrt_max_a,cyclic,
     * num_foldon_b,foldstrt_min_b, foldstrt_max_b,const_mode,
     * targ_cons,mem_cons,Q0_safe_a,Q0_safe_b,iss_bias,
     * i_Q_format_a, Q0_a, Q_weight_a,n_Q_anneal_a,Q0_inc_a, Q_clip_a,
     * seglist_a,numconst_a,seglist_b,numconst_b,i_bias_native_a,
     * i_Q_format_b, Q0_b, Q_weight_b,n_Q_anneal_b,Q0_inc_b, Q_clip_b,
     * ss_a,ss_b,
     * con_local_a, con_local_a_cut , con_local_b, con_local_b_cut,
     * i_bias_native_b,totalssnorm,
     * i_con_P_AP,i_atom_P_AP,weight_P_AP,r_cut_P_AP,i_diff_P_AP,
     * i_ignore_chain_dirn,mismatch,
     * go_con, go_con_dist,ssweight,aps,ss_pattern_a,ss_pattern_b,
     * i_rep,rep_cut_off,rep_tol,
     * i_rep_lambda_uniform,rep_lambda,n_letters_con,dist_cut,
     * itarg_seq,tgsequences,maxres,tgsequences_hb,tgsequences_amw,
     * tgsequences_amc 

      implicit none

c     internal variables:

         integer iunit,i,i_homlist,iseq,j,k,gap(maxsiz),
     *              tempgap, switch, ss_pat_iter_a,
     *              ss_pat_iter_b,temp_a, id3,nmrss,open_status
c     getval

         double precision rnorm 

         character*5 proflinit

c     list required subroutines

         external opnfil,qchrgmk,read_altgamma,read_input_alt
      
         data i_homlist/21/

	iunit = 998
	rnorm = 0.0D0

c     --------------------- begin -----------------------


c      write(6,*) 'here we go initil'

c      write(*,*) 'initil.f: here we go initil'
      call read_input_alt()  ! For alternative potential

c      write(SO,*) 'call opnfil'

        call  opnfil(iprolst,iprolstscl,input_amh,oarchv,
     *             omovi,ohdrgn,ohdrgn_s,ohdrgn_m,ohdrgn_l,ohdrgn_seq,
     *             orama,ooxy,ochiral,oamh,oamhsr,orep,
     *             oPE_no_bias,oPE_with_bias,oPE_backbone,
     *             oPE_backbone_norama,oamhlr,oamhmr,ononadd,
     *             ocon_P_AP,occev,ooev,obias_Rg,omoviseg,
     *             oobiassega,oobiassegb,oharmspring)


c             write(*,*) 'exit opnfil'

c       write(SO,*)'iprolst ',iprolst
c       write(SO,*)'input_amh ',input_amh
c       write(SO,*)'oarchv ',oarchv
c       write(SO,*)'omovi ',omovi
c       write(SO,*)'ohdrgn ',ohdrgn
c       write(SO,*)'ohdrgn_s ',ohdrgn_s
c       write(SO,*)'ohdrgn_m ',ohdrgn_m
c       write(SO,*)'ohdrgn_l ',ohdrgn_l
c       write(SO,*)'ohdrgn_seq ',ohdrgn_seq
c       write(SO,*)'orama ',orama
c       write(SO,*)'ooxy ',ooxy
c       write(SO,*)'ochiral ',ochiral
c       write(SO,*)'oamh ',oamh
c       write(SO,*)'ototal ',ototal
c       write(SO,*)'oamhsr ',oamhsr
c       write(SO,*)'oPE ',oPE
c       write(SO,*)'oKE ',oKE
c       write(SO,*)'oamhlr ',oamhlr
c       write(SO,*)'oamhmr ',oamhmr
c       write(SO,*)'ononadd ',ononadd
c       write(SO,*)'oamhlr ',oamhlr
c       write(SO,*)'oamhmr ',oamhmr
c       write(SO,*)'ononadd ',ononadd
c       write(SO,*)'ocon_P_AP',ocon_P_AP
c       write(SO,*)'obias_Rg',obias_Rg
c      write(SO,*)'omoviseg',omoviseg
c      write(SO,*)'ocrap',ocrap
c      write(SO,*)'oobiassega',oobiassega
c      write(SO,*)'oobiassegb',oobiassegb

c   read in various parameters and check that they are consistent with 
c   the maximum space allowed
c   nummem is the number of protein memories used in the potential

      read(input_amh,*)nummem
      if( nummem.gt.maxmem )then
         write(oarchv,720)nummem,maxmem
  720    format('nummem too large',2(1x,i3))
         stop
      endif

c     numpro is the number of ensemble proteins 
 
      read(input_amh,*)numpro

      write(SO,*)'numpro=',numpro

      if( numpro.gt.maxpro )then
         write(oarchv,721)numpro,maxpro
  721    format('numpro too large',2(1x,i3))
         stop
      endif

c     maxsiz is the maximum number of amino acid residues;
c     nmres is the actual protein size used

       read(input_amh,*)nmres

       write(SO,*)'nmres=',nmres

      if( nmres.gt.maxsiz )then
         write(oarchv,722)nmres,maxsiz
         write(SO,722)nmres,maxsiz
  722    format('nmres too large',2(1x,i3))
         stop
      endif

      read(input_amh,*) known

c     assign number of coordinate types

      numcrd=3

c     check that the number of specified
c     coordinate types does not exceed the maximum

      if( numcrd.gt.maxcrd )then
         write(oarchv,236)numcrd,maxcrd
  236    format(/'number of specified coordinate types ',i3,
     *          ' is greater than maximum ',i3)
         stop
      endif

c     assign number of different types of ixns

      numtab=4
     
c     check that the number of specified
c     interactions does not exceed the maximum

      if( numtab.gt.maxtab )then
         write(oarchv,237)numtab,maxtab
  237    format(/'number of specified tables ',i3,
     *          ' is greater than maximum ',i3)
         stop
      endif

c     assign coordinate types for each interaction

      crdixn(1,1)=1
      crdixn(1,2)=1

      crdixn(2,1)=1
      crdixn(2,2)=2

      crdixn(3,1)=2
      crdixn(3,2)=1

      crdixn(4,1)=2
      crdixn(4,2)=2
      
c     read in coordinate type to be used
c     (beta or center of mass)

c     sequence charge switches and their weights

      read(input_amh,*)igomb
      read(input_amh,*)ngomb

      write(SO,*) 'ngomb=',ngomb

      if (igomb) then
        write(oarchv,1724) ngomb
1724    format('Many-body factor ngomb=',f6.3)
       endif

c      degeneracy of our code

       read(input_amh,*) 
       read(input_amh,*)dist_cut
       read(input_amh,*)
       read(input_amh,*)ave_seq_amc
       read(input_amh,*)ave_seq_hb
       read(input_amh,*)ave_seq_amw

c       write(6,*)'ave_seq_amc ', ave_seq_amc
c       write(6,*)'ave_seq_hb ', ave_seq_hb
c       write(6,*)'ave_seq_amw ', ave_seq_amw

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read in sequences to average force over, if avg_seq is on
         do iseq = 1,maxseq
           do id3 = 1,maxres
             tgsequences(id3,iseq)=0
             tgsequences_hb(id3,iseq)=0
             tgsequences_amw(id3,iseq)=0
              tgsequences_amc(id3,iseq)=0
           enddo
         enddo

         rewind iprolst
         read (iprolst,1000)proflinit
1000     format(a5)
         write(SO,*) 'target in initil',proflinit
         rewind iprolst
            open(imemri,file='proteins/'//proflinit,
     *                    status='old',iostat=open_status)
         if (open_status.ne.0) then
           write(SO,*) 'failure to open file in initil',proflinit
           stop
         endif

         read(imemri,*)
         read(imemri,101)nmrss

c         write(6,*)'nmrss'
c         write(6,101)nmrss
101      format(i5)
         read(imemri,106)(ires(i1),i1=1,nmrss)
106      format(25(i2,1x))
         close(imemri)
         write(SO,*)'in initil ires'
         write(SO,106)(ires(i1),i1=1,nmrss)

          if (ave_seq_amc .or. ave_seq_hb .or. ave_seq_amw) then
           write(6,*)'target_sequences itarg_seq '
           open(itarg_seq,file='target_sequences',status='old')
           read(itarg_seq,*)numseq
454          format(i5)
c           write(6,*)'numseq ',numseq
             do iseq = 1,numseq
                  read(itarg_seq,999)(tgsequences(id3,iseq),id3=1,nmres)
c                  write(6,999)(tgsequences(id3,iseq),id3=1,nmres)
999               format(25(i2,1x))

           do id3 = 1,nmres
                if( (tgsequences(id3,iseq) .lt. 1) .or.
     *                    (tgsequences(id3,iseq) .gt. 20) ) then
                    write(6,*)'target_sequence  out of seq range '
                    write(6,*)'resnum   ', id3
                    write(6,*)'seqnum   ', iseq
                    write(6,*)'target_sequence  ', tgsequences(id3,iseq)
                    stop
                endif
           enddo

             enddo
            close(itarg_seq)
c            write(6,*)'itarg_seq closed ',itarg_seq
          endif

          if  (numseq.gt.maxseq) then
            write(6,*) 'numseq greater than maxseq'
            stop
          endif

         if ( ave_seq_amc ) then
            numseq_amc=numseq
          else
             numseq_amc=1
            tgsequences_amc(1:nmres,1)=ires(1:nmres)
          endif
c         write(6,*)'ave_seq_amc numseq_amc ', ave_seq_amc, numseq_amc

         if ( ave_seq_hb ) then
            numseq_hb=numseq
          else
            numseq_hb=1
            tgsequences_hb(1:nmres,1)=ires(1:nmres)
          endif
c         write(6,*)'ave_seq_hb numseq_hb ', ave_seq_hb, numseq_hb

       if ( ave_seq_amw) then
        numseq_amw=numseq
        write(6,*)'tgsequences_amw sequences'
         do iseq = 1,numseq_amw
            do id3 = 1,nmres
              tgsequences_amw(id3,iseq) = tgsequences(id3,iseq)
c              write(6,*)'tgseq_amw ave_seq_amw ',ave_seq_amw
c              write(6,*)'numseq_amw ', numseq_amw
c              write(6,106)(tgsequences_amw(i1,iseq),i1=1,nmrss)
            enddo
         enddo
        else
           numseq_amw=1
           tgsequences_amw(1:nmres,1)=ires(1:nmres)
c           write(6,*)'tgseq_amw ave_seq_amw ',ave_seq_amw
c           write(6,*)'numseq_amw ', numseq_amw
c           write(6,*)'numseq ', numseq
c           write(6,106)(tgsequences_amw(i1,1),i1=1,nmrss)
        endif
c        write(6,*)'ave_seq_amw numseq_amw ', ave_seq_amw, numseq_amw

c        if (ave_seq_amw) then
c             write(6,*)'tgsequences_amw sequences'
c         do iseq = 1,numseq_amw
c          do id3 = 1,nmres
c             tgsequences_amw(id3,iseq) = tgsequences(id3,iseq)
c             write(6,*)'tgsequences_amw sequences'
c             write(6,106)(tgsequences_amw(i1,1),i1=1,nmrss)
c          enddo
c         enddo
c        endif

       if (ave_seq_hb) then
             write(6,*)'tgsequences_hb sequences'
         do iseq = 1,numseq_hb
          do id3 = 1,nmres
             tgsequences_hb(id3,iseq) = tgsequences(id3,iseq)
             write(6,*)'tgsequences_hb sequences'
             write(6,106)(tgsequences_hb(i1,1),i1=1,nmrss)
          enddo
         enddo
        endif

        if ( ave_seq_amc) then
             write(6,*)'tgsequences_amc sequences'
         do iseq = 1,numseq_amc
          do id3 = 1,nmres
             tgsequences_amc(id3,iseq) = tgsequences(id3,iseq)
             write(6,*)'tgsequences_amc sequences'
             write(6,106)(tgsequences_amc(i1,1),i1=1,nmrss)
          enddo
         enddo
        endif
                       
c       write(6,106)(tgsequences_amc(i1,1),i1=1,nmrss)
                                      
       read(input_amh,*) targ_cons
       read(input_amh,*) mem_cons
       read(input_amh,*)ssweight
       read(input_amh,*)

c           write(6,*) 'targ_cons = ',targ_cons
c           write(6,*) 'mem_cons = ',mem_cons
c           write(6,*) 'ssweight = ',ssweight
          if((ave_seq) .and. (numseq.gt.maxseq) ) then
            write(6,*) 'numseq greater than maxseq'
            stop
          endif

        if (ssweight) then

        open(iss_bias,file='rama_ss_bias',status='old')

        do i = 1,nmres
        read(iss_bias,*)aps(i,5),aps(i,6)
        aps(i,1)=1.0D0
        aps(i,2)=1.0D0
        aps(i,3)=1.0D0
        aps(i,4)=1.0D0
        enddo

        close(iss_bias)
        else

        aps(:,1)=1.0D0
        aps(:,2)=1.0D0
        aps(:,3)=1.0D0
        aps(:,4)=1.0D0
        aps(:,5)=0.0D0
        aps(:,6)=0.0D0
        endif

c   gamma averaging

c      if (ave_seq) then
c            open(i_homlist,file='target_sequences',status='old')
c         do iseq = 1,numseq
c             read (i_homlist,1000)homlfl(iseq)
c1000         format(a5)
c             write(SO,*) 'homologue ',homlfl(iseq)
c         enddo
c         close (i_homlist)
c       endif  !ave_seq  avg_gammas


       read(input_amh,*) n_letters
c       write(SO,*)'n_letters ',n_letters
       if ((n_letters.ne.2).and.(n_letters.ne.4)) then
         write(SO,*) 'must be 2 or 4 letter code, not',n_letters
         stop
       endif
       read(input_amh,*) four_plus
       write(SO,*)'four_plus ',four_plus
       read(input_amh,*)go_con
       read(input_amh,*)go_con_dist
       write(SO,*)'go_con ',go_con
       write(SO,*)'go_con_dist ',go_con_dist
       read(input_amh,*)i_ignore_chain_dirn
       read(input_amh,*)mismatch
       write(SO,*)'mismatch ',mismatch

       read(input_amh,*) i_non_add_contact
       if (four_plus.and.(n_letters.ne.4)) then
         write(SO,*) 'can only extend 4 letter code if using it'
         stop
       endif

       if (n_letters.gt.max_letters) then
          write(SO,*) 'degeneracy too high'
          stop
       endif

       read(input_amh,*) min_seq_sep
       write(SO,*) ' min_seq_sep',  min_seq_sep
       read(input_amh,*) i_3rd_is_contact
       write(SO,*) 'i_3rd_is_contact' , i_3rd_is_contact
       if (i_3rd_is_contact.and.(n_letters.ne.4)) then
           write(SO,*) 'contact potential for 4 letter code only'
           stop
       endif

       read(input_amh,*) n_letters_con
       if (n_letters_con.ne.4.and.n_letters_con.ne.20) then
         write(SO,*) 'only 2,4 or 20 letter contact int supported',
     *                                                 n_letters_con
         stop
       endif

       read(input_amh,*) num_well
       write(SO,*)'num_well=',num_well
       if (num_well.ne.2.and.num_well.ne.3.and.num_well.ne.5
     *                                      .and.num_well.ne.10) then
         write(SO,*) 'only 3 or 10 well contact pot supported',num_well
         stop
       endif

       read(input_amh,*) gly_con
       do i=1,num_well
         read(input_amh,*) r_min(i),r_max(i)
         write(SO,*)'r_min(i),r_max(i)=',r_min(i),r_max(i)
       enddo
       read(input_amh,*) i_alt_prox
       write(SO,*)'i_alt_prox=',i_alt_prox
       read(input_amh,*) alt_prox_cut
       write(SO,*)'alt_prox_cut=',alt_prox_cut

       if (i_3rd_is_contact.and.(.not.i_alt_prox)) then
           write(SO,*) 'conatct potential only supported'
           write(SO,*) 'in Coreys alternative prox class'
           stop
       endif



c     excluded volume parameters

      read(input_amh,*)
      read(input_amh,*)iexcld
      read(input_amh,*)iexcld_beta
      read(input_amh,*)iexcld_gamma,i_type_ev_gamma
      if (i_type_ev_gamma.ne.1.and.i_type_ev_gamma.ne.2.and.i_type_ev_gamma.ne.3) then
        write(SO,*) 'i_type_ev_gamma cock-up',i_type_ev_gamma
        stop
      endif
c  pexcld  penalty for violating excluded volume
      read(input_amh,*)pexcld
       write(SO,*)'pexcld=',pexcld
      read(input_amh,*)exvmin
       write(SO,*)'exvmin=',exvmin
      read(input_amh,*)exvminS
       write(SO,*)'exvminS=',exvminS
      read(input_amh,*)exvminS_beta(1:4)
       write(SO,*) 'exvminS_beta(1:4)=',exvminS_beta(1:4)
      read(input_amh,*)oxexcldv
       write(SO,*)'oxexcldv=',oxexcldv
      read(input_amh,*) o_exvmin
       write(SO,*)'o_exvmin=',o_exvmin
      read(input_amh,*) o_exvminS(1:4)
       write(SO,*) 'o_exvminS(1:4)',o_exvminS(1:4)
      read(input_amh,*)oexcldscl
      write(SO,*) 'oexcldscl',oexcldscl

      pexcld_gamma=pexcld

         hrdrad(0,1)=0.0D0
         hrdrad(0,2)=0.0D0
         do 505 i505=1,21
            hrdrad(i505,1)=exvmin*0.5D0
            hrdrad(i505,2)=exvmin*0.5D0
  505    continue
ccccccc

      read(input_amh,*)

c     ran_force is flag for random force
      read(input_amh,*) ran_force
      write(SO,*) 'ran_force ',ran_force
      if (ran_force  .and. nummem .gt. 1) then
        write(SO,*) 'random force and >1 memory, so stopping'
        stop
      endif

c     use previous random data if ran_file true

      read(input_amh,*) ran_file

c     r_ran  is eqm dist for random force

      read(input_amh,*) r_ran

c     s_ran  is sd of randon int

      read(input_amh,*) s_ran

c     ran_min_seq_dist is minimum seq-sep for random ints

      read(input_amh,*) ran_min_seq_dist

c     allow_neg allows negative prefactors to random potential if true

      read(input_amh,*) allow_neg

c     lambdaR is the scale factor given to
c     incorrect memories or random ints

      read(input_amh,*) lambdaR

      write(SO,*) 'lambdaR',lambdaR

c     srplr is ratio of total short-med-long range (in seq) energy
c     in native state if iclassweight=.true.

      read(input_amh,*)
      read(input_amh,*)srplr

      write(SO,*) 'srplr',srplr

c     eres is the desired potential energy/residue

      read(input_amh,*)eres

c     short-range ixn cutoff
      write(SO,*) 'bottom 2a'

      read(input_amh,*)minmr,maxmr
      write(SO,*) 'minmr, maxmr',minmr,maxmr

c     maxs is the r-grid resolution

      read(input_amh,*)maxs

c     check that maxs is less than or equal to the
c     maximum number of r-grid points

      if( maxs.gt.maxr )then
         write(oarchv,726)maxs,maxr
  726    format('maxs too large',2(1x,i3))
         stop
      endif

c     read switch for rama potential


      read(input_amh,*)
      read(input_amh,*) i_rama

      if (i_rama.ne.0 .and. i_rama.ne.1 .and. i_rama.ne.2 
     *  .and. i_rama.ne.3) then
         write(SO,*) 'i_rama set incorrectly',i_rama
         stop
      endif

c     read hbond flag
      write(SO,*) 'reading hbond'

      read(input_amh,*)hbond
      read(input_amh,*)i_hbond_eastwood ! old dssp_hbond
      read(input_amh,*)i_hbond_czong

       if((hbond.and.(i_hbond_eastwood.or.i_hbond_czong)).or.
     * (hbond .and. (i_hbond_eastwood .or. i_hbond_czong))) then
         write(SO,*) 'cannot have both sorts of H-bond'
         write(SO,*) 'for a start they use same avep array'
         stop
       endif

c     input_amh in switch for oxygen potential and scaling factor
c     and hbscl

      write(SO,*) 'reading input_amh in switch'
      read(input_amh,*)NO_zero
      read(input_amh,*)sigma_NO
      read(input_amh,*)ho_zero
      read(input_amh,*)sigma_h
      read(input_amh,*)hbscl(1:17)


c     input_amh in switch for oxygen potential and scaling factor
c     and hbscl

      write(SO,*) 'reading input_amh in switch'

      read(input_amh,*)oxscl
      write(SO,*) 'oxscl',oxscl
      read(input_amh,*)ramascl
      write(SO,*) 'ramascl', ramascl
c     chiral forces
      read(input_amh,*)chrlscl
      write(SO,*)'chrscl',chrlscl

! read in parameters for replica force
      read(input_amh,*) 
      read(input_amh,*) i_rep
      write(SO,*) i_rep
      read(input_amh,*) rep_cut_off
      write(SO,*) rep_cut_off
      read(input_amh,*) rep_tol
      write(SO,*) rep_tol
      read(input_amh,*) i_rep_lambda_uniform,rep_lambda(1)
      write(SO,*) i_rep_lambda_uniform,rep_lambda(1)

c     central potential switch and scaling factor
c     for central potential

      write(SO,*) 'bottom 2b'

c     MD timestep

      read(input_amh,*)
      read(input_amh,*)timstp
      write(SO,*) 'timstp',timstp

c     set tolerance for SHAKE routines
c
c
c Changing this temporarily
c tolshk=0.2

c        tolshk=200000000*timstp
       tolshk=0.2*timstp
c     minimal well-width

      read(input_amh,*)delta

c     exponent for well-width as function |i-j-1|

      read(input_amh,*)delte

c     maximum interaction distance

      read(input_amh,*)rcut
      read(input_amh,*)rcutAMH
      read(input_amh,*)srcut
      read(input_amh,*)welscl
      write(SO,*) 'rcut,rcutAMH,srcut,welscl',rcut
     *                              ,rcutAMH,welscl


c     read parameters/flags for biasing functions

      read(input_amh,*)
      read(input_amh,*) ibiasgauss           ! first biasong funcs that depend on E_amh
      read(input_amh,*) bias_weight
      read(input_amh,*) bias_av
      read(input_amh,*) bias_var
      bias_prefactor=bias_weight/sqrt(2.0*3.1415926*
     *                                        bias_var)
      if (.not.ibiasgauss) bias_prefactor=0.0D0
      read(input_amh,*) ibiaspoly
      read(input_amh,*) nbiaspoly
      read(input_amh,*) (biaspoly(i),i=1,nbiaspoly)
      if (nbiaspoly.gt.100) then
         write(SO,*) 'nbias too high'
         stop
      endif
      if (ibiasgauss.and.ibiaspoly) then
        write(SO,*) 'polynomial + gaussian bias cannot
     *                 be on at same time'
         stop
      endif
      do i=nbiaspoly+1,100
        biaspoly(i)=0.0D0
      enddo
ccc  part for partial contstraint
ccc  mode 1  for generate list for f of qs
ccc  mode 2  for two seperate bits with flexible part in between

       read(input_amh,*)
       read(input_amh,*)cyclic
       write(SO,*)'cyclic = ',cyclic
       read(input_amh,*)const_mode
       write(SO,*)'const_mode = ',const_mode
       read(input_amh,*) num_foldon_a
       write(SO,*)'num_foldon_a',num_foldon_a
c        if (num_foldon_a .eq. 0) goto 75
         if (num_foldon_a .gt. 30) STOP 'Too many Foldons: Max = 30'
       do i = 1,num_foldon_a
          read(input_amh,*) foldstrt_min_a(i),foldstrt_max_a(i)
          write(SO,*) foldstrt_min_a(i),foldstrt_max_a(i),i
          if( foldstrt_max_a(i).gt.nmres )then
          write(oarchv,724)foldstrt_max_a,nmres
724       format('foldstrt_max_a too large',2(1x,i3))
          stop
          endif
       enddo

       if( const_mode .eq. 2 )then
              read(input_amh,*) num_foldon_b
              write(SO,*) num_foldon_b
            if (num_foldon_b .gt. 20) STOP 'Too many Foldons: Max = 20'
                do i = 1,num_foldon_b
                   read(input_amh,*) foldstrt_min_b(i),foldstrt_max_b(i)
                   write(SO,*) foldstrt_min_b(i),foldstrt_max_b(i)
                if( foldstrt_max_b(i).gt.nmres )then
                    write(oarchv,725)foldstrt_max_b,nmres
725                 format('foldstrt_max_b too large',2(1x,i3))
                     stop
                 endif
            enddo
        endif ! const_mode

      read(input_amh,*)
      read(input_amh,*) i_Qbias_a          ! now biasing funcs that dep on Q
      write(SO,*)'iQbias_a ',i_Qbias_a
      read(input_amh,*)n_Q_anneal_a, Q0_inc_a
      write(SO,*)'n_Q_anneal_a Q0_inc_a', n_Q_anneal_a, Q0_inc_a
      read(input_amh,*)con_local_a, con_local_a_cut
      write(SO,*)'con_local_a  con_local_a_cut',
     *                             con_local_a,con_local_a_cut
      read(input_amh,*) i_bias_native_a
      write(SO,*)'i_bias_native_a',  i_bias_native_a
      read(input_amh,*) ss_a
      write(SO,*)'ss_a ', ss_a
      write(SO,*)'simulate s.s. units ss_a ', ss_a
      if (.not. ss_a) ss_pat_iter_a = 1
      if (ss_a) ss_pat_iter_a = num_foldon_a 
      write(SO,*)'ss_pat_iter_a ',ss_pat_iter_a
      read(input_amh,*)(ss_pattern_a(i),i=1,ss_pat_iter_a)
       do i = 1, ss_pat_iter_a
            write(SO,*)'ss_pattern_a ',  ss_pattern_a(i)
       enddo

      read(input_amh,*) width_Qexp_a
      write(SO,*)'width_Qexp_a',  width_Qexp_a
      read(input_amh,*) n_Qbias_a
      write(SO,*)'n_Qbias_a', n_Qbias_a
      read(input_amh,*) i_Q_format_a
      write(SO,*)'i_Q_format_a', i_Q_format_a 

      read(input_amh,*) (Qbiaspoly_a(i),i=1,n_Qbias_a)
      if (n_Qbias_a.gt.100) then
         write(SO,*) 'n_Qbias too high'
         stop
      endif
      read(input_amh,*) Q0_safe_a,Q_weight_a
      write(SO,*)'Q0_safe_a  ', Q0_safe_a
      read(input_amh,*) Q_clip_a
      write(SO,*)'Q_clip_a  ', Q_clip_a

       Q0_a =  Q0_safe_a

      do i=n_Qbias_a+1,100
        Qbiaspoly_a(i)=0.0D0
      enddo

c  count number of contacts for ss_a
          if(ss_a)then
            totalssnorm = 0
          do i = 1,num_foldon_a
              temp_a = 0
           do j = foldstrt_min_a(i), foldstrt_max_a(i)
              temp_a  = temp_a + 1
           enddo
              totalssnorm = totalssnorm + (temp_a-1)*(temp_a-2)
          enddo
              write(SO,*)'totalssnorm',totalssnorm
          endif ! ss_a

ccccccccccccccc   NEWNEWNEWNEWNEWNEW

      read(input_amh,*)
      read(input_amh,*) i_Qbias_b         ! now biasing funcs that dep on Q
      write(SO,*)'i_Qbias_b ',i_Qbias_b
      read(input_amh,*)n_Q_anneal_b, Q0_inc_b
      write(SO,*)'n_Q_anneal_b, Q0_inc_b ', n_Q_anneal_b, Q0_inc_b
      read(input_amh,*)con_local_b, con_local_b_cut
      write(SO,*)'con_local_b  con_local_b_cut',
     *                                  con_local_b,con_local_b_cut
      read(input_amh,*) i_bias_native_b
      write(SO,*)'i_bias_native_b ', i_bias_native_b
      read(input_amh,*) ss_b
      write(SO,*)'simulate s.s. units  ss_b ', ss_b 
      if (.not. ss_b) ss_pat_iter_b = 1
      if (ss_b) ss_pat_iter_b = num_foldon_b
      write(SO,*)'ss_pat_iter_b ',ss_pat_iter_b
      read(input_amh,*)(ss_pattern_b(i),i=1,ss_pat_iter_b)
      read(input_amh,*) width_Qexp_b
      read(input_amh,*) n_Qbias_b
      read(input_amh,*) i_Q_format_b
      read(input_amh,*) (Qbiaspoly_b(i),i=1,n_Qbias_b)

      if (n_Qbias_b .gt. 100) then
         write(SO,*) 'n_Qbias_b too high'
         stop
      endif
      read(input_amh,*) Q0_safe_b,Q_weight_b
      read(input_amh,*) Q_clip_b

      Q0_b = Q0_safe_b


      do i=n_Qbias_b+1,100
        Qbiaspoly_b(i)=0.0D0
      enddo

c make constraint lists

          do i = 1,  maxsiz
             gap(i) = 0
          enddo

          do i = 1,  nmres
	                seglist_a(i) = 0
          enddo

          if ( num_foldon_a .ne. 1 .and. num_foldon_a .ne. 0) then
             do i = 1,num_foldon_a - 1
               gap(i) = foldstrt_min_a(i+1) - foldstrt_max_a(i) - 1
             enddo
          endif

          do i = 1,num_foldon_a
          if ( i .eq. 1 .and. num_foldon_a .ne. 0 ) then
            do j = foldstrt_min_a(i), foldstrt_max_a(i)
                seglist_a(j - (foldstrt_min_a(i)-1)) = j
                numconst_a = j - (foldstrt_min_a(i)-1)
            enddo
          endif
          tempgap = 0
          if ( i .ne. 1 .and. num_foldon_a .ne. 0 ) then
              do k = 1 , i - 1
                  tempgap = tempgap + gap(k)
              enddo
              do j = foldstrt_min_a(i), foldstrt_max_a(i)
                  seglist_a(j - tempgap -  (foldstrt_min_a(1)-1)) = j
                  numconst_a = j - tempgap  - (foldstrt_min_a(1)-1)
              enddo
          endif   ! if ( i .ne. 1) then
          enddo

             write(SO,*)'numconst_a',numconst_a
             do i = 1, numconst_a
               write(SO,*)'seglist_a ',  seglist_a(i)
             enddo


cc construct  second list

       if (const_mode .eq. 1) then
             do i = 1,  nmres
                       seglist_b(i)= 0
             enddo
             numconst_b  = 0
             do  266 k = 1 ,  nmres
                 switch = 0
              if (num_foldon_a .eq. 0) numconst_a = nmres
                  do 255 j = 1, numconst_b
                     if  (seglist_b(j) .eq. k) switch = 1
255            enddo
                 if (switch .eq. 0) then
                    numconst_b  = numconst_b  + 1
                    seglist_b(numconst_b) =  k
                 endif
266           enddo
      endif  !    listmode 1

c   lists for manual assignment of two lists for q constraints

       if (const_mode .eq. 2) then
             do i = 1,  nmres
                       seglist_b(i)= 0
             enddo
                numconst_b  = 0

          if ( num_foldon_b .ne. 1 .and. num_foldon_b .ne. 0) then
             do i = 1,num_foldon_b - 1
                  gap(i) = foldstrt_min_b(i+1) - foldstrt_max_b(i) - 1
             enddo
          endif

          do i = 1,num_foldon_b
          if ( i .eq. 1 .and. num_foldon_b .ne. 0 ) then
            do j = foldstrt_min_b(i), foldstrt_max_b(i)
                seglist_b(j - (foldstrt_min_b(i)-1)) = j
                numconst_b = j - (foldstrt_min_b(i)-1)
            enddo
          endif
          tempgap = 0
          if ( i .ne. 1 .and. num_foldon_b .ne. 0 ) then
              do k = 1 , i - 1
                  tempgap = tempgap + gap(k)
              enddo
              do j = foldstrt_min_b(i), foldstrt_max_b(i)
                  seglist_b(j - tempgap -  (foldstrt_min_b(1)-1)) = j
                  numconst_b = j - tempgap  - (foldstrt_min_b(1)-1)
              enddo
          endif   ! if ( i .ne. 1) then
          enddo

           write(SO,*)'numconst_b',numconst_b
                  do i = 1,  numconst_b
             write(SO,*)'seglist_b ',  seglist_b(i)
                  enddo

        endif  !    listmode 2

c       if (con_local_a) then
c          j = 0 
c          do i= 1, numconst_a
c           if((seglist_a(i)+1-seglist_a(i)).lt.con_local_a_cut)then
c                j = j  + 1
c           endif
c          enddo
c         numconst_a = j
c       endif   ! if (con_local_a)  then
c
c
c       if ( con_local_b) then
c        j = 0 
c         do i= 1, numconst_b
c           if((seglist_b(i)+1 -seglist_b(i)) .lt. con_local_b_cut)then 
c             j = j  + 1 
c           endif 
c         enddo
c        numconst_b = j 
c       endif   ! if (con_local_b)  then
c
cccccccccccccccccc NEWNEWNEWNEWNEWNEW

      read(input_amh,*)
      read(input_amh,*) i_contact_order
      write(SO,*)'i_contact_order ',i_contact_order
      read(input_amh,*) n_contact_order_terms
      if (n_contact_order_terms.gt.max_contact_order_terms) then
         write(SO,*) 'n_contact_order_terms too big'
         stop
      endif
      do i=1,n_contact_order_terms
        read(input_amh,*) i_contact_order_min(i),i_contact_order_max(i)
        read(input_amh,*) r_min_contact_order(i),r_max_contact_order(i)
        read(input_amh,*) gamma_contact_order(i,1:3)
      enddo

ccccccccccccccccccccHydrogoen bond parameters

      read(input_amh,*)
      read(input_amh,*) i_con_P_AP,i_atom_P_AP
      write(SO,*)'i_con_P_AP   i_atom_P_AP ',i_con_P_AP, i_atom_P_AP
      read(input_amh,*) r_cut_P_AP,i_diff_P_AP
      write(SO,*)'r_cut_P_AP  i_diff_P_AP ',r_cut_P_AP, i_diff_P_AP
      read(input_amh,*) weight_P_AP(1:3)

      read(input_amh,*)
      read(input_amh,*) alpha_c_of_n
      read(input_amh,*) ab_c_of_n_new
      read(input_amh,*) ab_c_of_n_old
      if ( (alpha_c_of_n.and.ab_c_of_n_new) .or.
     *     (alpha_c_of_n.and.ab_c_of_n_old) .or.
     *     (ab_c_of_n_new.and.ab_c_of_n_old) ) then
         write(SO,*) 'c of n problem',alpha_c_of_n
         write(SO,*) 'c of n problem',ab_c_of_n_new
         write(SO,*) 'c of n problem',ab_c_of_n_old
         stop
      endif

      read(input_amh,*)
      read(input_amh,*) i_Rg_bias          ! now biasing funcs that dep on Rg
      read(input_amh,*) i_rg_corey
      read(input_amh,*) rg_bounds(1)
      read(input_amh,*) rg_bounds(2)
      read(input_amh,*) rg_shift, rg_scl
      read(input_amh,*) i_rg_garyk
      read(input_amh,*) D_rg, T_rg, delR_rg, M_rg, kappa_rg
      read(input_amh,*) n_Rg_bias
      read(input_amh,*) (Rg_biaspoly(i),i=1,n_Rg_bias)
      if (n_Rg_bias.gt.100) then
         write(SO,*) 'n_Rg_bias too high'
         stop
      endif
                                                                                                     
      if (i_rg_corey.and.i_rg_garyk) then
         write(SO,*)
     *     'Rg_bias: you must choose between i_rg_corey and i_rg_garyk'
         stop
      endif
                                                                                                     
      i_rg_first=.true.

      do i=n_Rg_bias+1,100
        Rg_biaspoly(i)=0.0D0
      enddo

c     nmstep is the number of time steps per temperature

      read(input_amh,*)nmstep

c     read seed values for random number generators

      read(input_amh,*)iseed_amh
      write(SO,*) 'iseed_amh',iseed_amh

c     it1-it5 are then number of T grid points allocated to each
c     of the intervals [t1,t2], [t2,t3],...

c     t1-t5 are the cutoff points in the temperature-annealing
c     schedule

      read(input_amh,*)
      read(input_amh,*)
      read(input_amh,*)it1,t1
      read(input_amh,*)it2,t2
      read(input_amh,*)it3,t3
      read(input_amh,*)it4,t4
      read(input_amh,*)it5,t5

      write(SO,*) it1,it2,it3,it4,it5
      write(SO,*) t1,t2,t3,t4,t5

c     if ictemp, then the minimization is performed at a
c     constant temperature, ctemp

      read(input_amh,*)ictemp
      read(input_amh,*)ctemp

      write(SO,*)'ictemp ',ictemp
      write(SO,*)'ctemp  ',ctemp

c     nmtemp is the number of temperature grid points

      nmtemp=it1 + it2 + it3 + it4
      if( mxtemp.lt.nmtemp )then
         write(oarchv,981)mxtemp,nmtemp
  981    format(/'mxtemp ',i6,' too small ',i6)
         stop
      endif

       write(SO,*) 'bottom 3b'

c     incmov is the increment for the number of
c     structures to be saved fixed T

      read(input_amh,*)
      read(input_amh,*)incmov
      rnorm=float(nmtemp)/float(incmov)
      if( int(rnorm).gt.nmdif )then
         write(oarchv,728)incmov
  728    format(/'incmov too small ',i3)
         stop
      endif

      read(input_amh,*)

      read(input_amh,*) i_V_test
      read(input_amh,*) test_site(1),test_site(2)

      read(input_amh,*)


c     if iresrt, then read in previously generated
c     proteins used for restarting program or analysis

      read(input_amh,*)iresrt

c     if iscltab go into scltab

      read(input_amh,*) iscltab

c     if iclassweight scale gammas according to srplr
      read(input_amh,*) iclassweight
      if (iclassweight.and.igomb) then
         write(SO,*) 'cannot divide energies into s/lr if igomb=true'
         stop
      endif

c      if movanal, then perform analysis of movie strs only

      read(input_amh,*) movanal
      read(input_amh,*) quench
      write(SO,*) 'quench=',quench 
      read(input_amh,*) i_Etim,i_Etim_step

c     check nmdif large enough -- rolling D average


      rnorm=float(nmstep)/float(incmov)
      if( int(rnorm).gt.nmdif )then
         write(oarchv,251)nmstep,incmov,nmdif
  251    format(/'nmdif too small:nmstep ',i5,
     *           ' incmov ',i5,' nmdif ',i5)
         stop
      endif


c     set annealing grid parameters

      itgrd(1)=it1
      itgrd(2)=it2
      itgrd(3)=it3
      itgrd(4)=it4
  
      temgrd(1)=t1
      temgrd(2)=t2
      temgrd(3)=t3
      temgrd(4)=t4
      temgrd(5)=t5

c     set charge weights



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     print header for archive file

      write(SO,*) 'about to print header'
 
            iunit=oarchv

c        echo relevant parameters

         if (movanal) write(iunit,*) 'movie str analysis only'

c         write(iunit,145)
c  145    format('PARAMETER LIST'/)

c         write(iunit,102)nmres
c  102    format('protein length ',i3)

c         write(iunit,103)numpro
c  103    format('ensemble size ',i3)

c        potential parameters

        write(iunit,151)
  151    format(/'POTENTIAL PARAMETERS'/)

         write(iunit,141)nummem
  141    format('actual # protein memories ',i4)

         write(iunit,270)numcrd
  270    format(/'number of coordinate types ',i2,
     *          ' beta-coordinates ',l1)

         write(iunit,271)numtab,
     *                  ((crdixn(i1,i2),i2=1,2),i1=1,numtab)
  271    format('Ixns: ',i2,' coordinates ',4('(',i2,',',i2,')'))

         write(iunit,107)timstp
  107    format(/'timstp ',1pe10.3)
       
         write(iunit,378)oxscl,ramascl
  378    format('ox pot on: ',' oxy scl ',f8.3,' rma scl ',f8.3)

           write(iunit,376)chrlscl
  376      format('chiral pot on: ',' scaling factor ',f8.3)

c         write(iunit,11379) lambdaR
c11379    format('lambdaR: ',f8.3)

c            write (iunit,123)
c  123       format('charges=hydrophobicity ')

c         write(iunit,148)delta,delte
c  148    format(/'well-width ',1pe10.3,' exponent ',
c     *          1pe10.3)

c         write(iunit,817)iexcld,pexcld
c  817    format(/'Excluded volume ',l1,
c     *          ' harmonic coefficient ',1pe10.3)

c	  write(iunit,818)oxexcldv,oexcldscl
c  818	  format(/'Ox. Exlcuded Volume',l1,
c     *          'Scale factor ', 1pe10.3)


c       print seeds for random number generators

c         write(iunit,844)
c  844    format(/'Seeds for Random Number Generators')
c         write(iunit,843)iseed_amh

c  843    format('iseed_amhs ',4(i10,2x))


c         write(iunit,795)welscl
c  795    format('Well scale ',l1)

c        annealing parameters

c         write(iunit,140)
c  140    format(/'ANNEALING PARAMETERS'/)

c         write(iunit,104)nmtemp,nmstep
c  104    format('# Temperature sets ',i8,
c     *          ' # time steps/T ',i6)

c         write(iunit,113)ictemp,ctemp
c  113    format('fixed T ',l1,' @T=',1pe12.5)

c         write(iunit,114)
c  114    format('temperature-annealing grid')

c         write(iunit,115)(itgrd(i1),i1=1,4),nmtemp
c  115    format(5(i12,1x))

c         write(iunit,116)(temgrd(i1),i1=1,5)
c  116    format(5(1pe12.5,1x))


c        print h scale

         do 503 i503=1,20
            hydscl(i503,2)=hydscl(i503,1)
  503    continue
         hydscl(8,2)=0.0D0

         write(iunit,118)
  118    format(/'Hard-Sphere Radii'/)

c         write(oarchv,868)
c  868    format('Reduced hard-sphere radii on Beta-carbons'/)

         write(iunit,389)
  389    format('alpha-carbons')
         write(iunit,119)(aminoa(i1),hrdrad(i1,1),i1=1,20)
  119    format(5(a4,1x,f5.2,2x))
         write(iunit,*)
         write(iunit,379)
  379    format('beta-carbons')
         write(iunit,119)(aminoa(i1),hrdrad(i1,2),i1=1,20)

c     end header
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     create qchrg
      call qchrgmk(qchrg,qchrg2,oarchv,n_letters)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     create alternative gamma (or qchrg if you will)
c      call read_altgamma()


c     set diagnostic flag

      idigns=.true.
c     ---------------------- done -----------------------
      return
      end
