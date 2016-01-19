c   --------------------- gaspot ----------------------

      subroutine gaspot(maxsd,target,bondln,curtab,nmcrd,cdtyp1,cdtyp2,passi,passf)

c     ---------------------------------------------------

c     GASPOT find ungeneralized potential

ccccccccccccccccccccccccccccccccccccccccccccccc
c

c     arguments:

c        maxsd - maximum protein length (i)
c        target- target protein (i)
c        bondln- target protein's bond lengths (o)
c        curtab- index for current table to be
c                constructed (i)
c        nmcrd- number of coordinate types (i)
c        cdtyp1- tpye id for coordinate set 1 (i)
c        cdtyp2- tpye id for coordinate set 2 (i)
c        passi - true on first pass through 
c                subroutine; otherwise false (i)[
c        passf - true on last pass through 
c                subroutine; otherwise false(i)
c
c      Variable names lifted from gaussv.f
c        gaussp- gaussian as a function of the r-grid
c                in rgrid (o)
c        rgrid - grid of r points for which the gaussian
c                is to be computed (i)        
c      -----------------------------------------------
c
c
c
c        dist_ij        distance i to j in the memory
c        id1
c     ---------------------------------------------------
c
c  NB:
c      Part of the reason for the length, and confusion in this subroutine is
c     historical in nature.  Before homology matching of residues, each 
c     residue was compared to other residues in a sliding window,
c     from -n11 to n11, and from -n33 to n33, were compared.
c
c
c Program Logic
c -------------
c       rtlist = 0
c
c
c509        Loop over memories (25%)
c                Call getmem (29%)
c               Call hprofl (hydscl
c                Call hprofl (hydscl
c
c                If first memory (target)  -->
c                       Call pttarg
c                       Continue to next memory (35%)
c
c                If first memory protien;
c                        Compare names protnm(1) and protnm(0)
c
c               Read in match file
c
c519                Loop over each interaction (43%)
c
c
c                        Determine category of interaction (55%)
c                        Call charge
c                       If charge = 0, to to 519 (? next ixn)
c                        Call gaussv to create gaussian curve around dist_ij
c
c                        Continue to next interaction (80%)
c
c      Continue to next memory (80%)
c
      use amhglobals,  only:SO, maxs,maxcnt,maxmem,numlng,nmres,
     *      iprolst,vpotnt,icon,forse,dist_cut,
     *      idigns,nummem,protnm,imemri,maxres,jres,numcrd,
     *      mempre,oarchv,sa,iwork,aminoa,hydseq,
     *      maxsiz,ires,shydro,tarpre,numtab,ran_force,
     *      ran_file,iran,iseed_amh,s_ran,lambdaR,allow_neg,oran,ilong,
     *      i_alt_prox,rinc,srcut,alt_prox_cut,ran_min_seq_dist,
     *      qchrg,qchrg2,deltz,welscl,r_ran,
     *      rincinv,rincsq,work2,four_plus,
     *      ywork,work8,hydscl,min_seq_sep,i_3rd_is_contact,
     *      n_letters,r_min,r_max,gly_con,i_V_test,test_site,
     *      ixn_from_site,minmr,maxmr,max_well,alpha_c_of_n,
     *      ab_c_of_n_old,ab_c_of_n_new,num_well,i_contact_order,
     *      i_contact_order_min,i_contact_order_max,r_min_contact_order,
     *      r_max_contact_order,gamma_contact_order,imat,
     *      n_contact_order_terms,numseq,ave_seq,itarg_seq,
     *      ave_seq_amc,numseq_amc,tgsequences_amc,
     *      ovgaspot,maxseq,tgsequences,targ_cons,go_con,go_con_dist

      implicit none


c     argument declarations:

         logical passi,passf,iprox_ran

         integer maxsd,curtab,cdtyp1,cdtyp2,nmcrd,iprox,helix,i_contact_term
     
         double precision target(maxsd,3,nmcrd),bondln(maxsd,nmcrd)

c     internal variables:

         logical passt,helix_A,helix_B
         double precision rgrid(maxs), gaussp(maxs),gaussq(maxs),ran_num(1:maxcnt),delt_safe,theta,rnmres,
     *        theta_dot,force_term(6),vpotint_test(1:maxs+1),t_min,t_max,r_temp

         integer indxt,nmrss(0:maxmem),maxsz,maxs2, match(maxres),j,i_r,ii,
     *           i,id1, id2, id3, id4,idummy,i_test_tab,
     *           i_test_ixn,tab_for_contact,open_status,count,i_mem, i_ixn, i_res,i_well

        integer i514,i515,i544,i1,i2,i516,i507,i522,i588,i504,i523,i589,i503,i532,i518

         integer iseq,connmres,temp_mcp_pass

         double precision rnorm,delt2,totchg(max_well),long_nfactor(max_well),dist_ij,m,c,dist_prox

         character*5 profl,tarfl
         character*42 confile
         character*33 profile
         character*36 matfile
         character*10 ccount1,ccount2,ccount3
         character*30 blah
c     required subroutines

         external getmem,gaussv,hprofl,pttarg,gaussw,num_to_char,SLARNV

c     --------------------- begin -----------------------
c     set various variables

      indxt=numlng(nmres,curtab)

c     rewind data base file
      rewind iprolst

c     set number points to be analyzed or set the table size

       maxsz=maxs
       maxs2=maxs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   set long_nfactor which is some fudge factor for the
c   contact potentials (it partially allows for the effect of 
c   size -- which the fraction of contacts at a certain distance
c   depends on). It should be the same as included in the
c   optimisation. In corey's code it is in qchrgmk.f instead.
c   Johan: added long_nfactor for 2-wells 

        long_nfactor=1.0D0

c   start if alpha_c_of_n 
        if ( alpha_c_of_n ) then

       do ii=1,max_well,1
         long_nfactor(ii)=1.0D0
       enddo
       rnmres=real(nmres)

        if (num_well .eq. 2) then
         long_nfactor(2)=1.0D0/(rnmres*0.012D0 + 0.87D0)
        elseif (num_well .eq. 3) then
         long_nfactor(2)=1.0D0/(rnmres*0.0065D0 + 0.87D0)
         long_nfactor(3)=1.0D0/(rnmres*0.04187261D0 + 0.1256658D0)
        elseif (num_well .eq. 10) then
         long_nfactor(1)=1.0D0/(rnmres*0.0008D0 + 0.09D0)
         long_nfactor(2)=1.0D0/(rnmres*0.0009D0 + 0.16D0)
         long_nfactor(3)=1.0D0/(rnmres*0.001D0 + 0.20D0)
         long_nfactor(4)=1.0D0/(rnmres*0.003D0 + 0.29D0)
         long_nfactor(5)=1.0D0/(rnmres*0.004D0 + 0.53D0)
         long_nfactor(6)=1.0D0/(rnmres*0.004D0 + 0.76D0)
         long_nfactor(7)=1.0D0/(rnmres*0.005D0 + 0.77D0)
         long_nfactor(8)=1.0D0/(rnmres*0.005D0 + 0.94D0)
         long_nfactor(9)=1.0D0/(rnmres*0.006D0 + 1.18D0)
         long_nfactor(10)=1.0D0/(rnmres*0.013D0 + 1.8D0)
        else
           write(SO,*) 'unsupported no of wells (alpha_c_of_n)',num_well
           stop
        endif

       endif
c   end if alpha_c_of_n 

       if (ab_c_of_n_new) then
       rnmres=real(nmres)
        if (num_well .eq. 2) then
          long_nfactor(1)= ( 0.035D0*rnmres)/(rnmres*0.043D0 + 1.0D0)
          long_nfactor(2)=( 0.07D0*rnmres)/(rnmres*0.023D0 + 1.0D0)
           long_nfactor(1)=1.0D0/(long_nfactor(1))
           long_nfactor(2)=1.0D0/(long_nfactor(2))
        elseif (num_well .eq. 3) then
          long_nfactor(1)= ( 0.0843467D0*rnmres)/(rnmres*0.0453928D0 + 1.0D0)
          long_nfactor(2)=( 0.0669808D0*rnmres)/(rnmres*0.025112D0 + 1.0D0)
          long_nfactor(3)=( 0.18665D0*rnmres) /(rnmres*0.0107983D0 + 1.0D0)
           long_nfactor(1)=1.0D0/(long_nfactor(1))
           long_nfactor(2)=1.0D0/(long_nfactor(2))
           long_nfactor(3)=1.0D0/(long_nfactor(3))
        elseif (num_well .eq. 5) then

         long_nfactor(1)=( 0.0297375D0*rnmres) /(rnmres*0.02977935D0 + 1.0D0)
         long_nfactor(2)=( 0.0389704D0*rnmres) /(rnmres*0.021101D0 + 1.0D0)
         long_nfactor(3)=( 0.0596751D0*rnmres) /(rnmres*0.0133269D0  + 1.0D0)
         long_nfactor(4)=( 0.0681322D0*rnmres) /(rnmres*0.0100256D0 + 1.0D0)
         long_nfactor(5)=( 0.0729201D0*rnmres) /(rnmres*0.00347563D0 + 1.0D0)

         long_nfactor(1)=1.0D0/(long_nfactor(1))
         long_nfactor(2)=1.0D0/(long_nfactor(2))
         long_nfactor(3)=1.0D0/(long_nfactor(3))
         long_nfactor(4)=1.0D0/(long_nfactor(4))
         long_nfactor(5)=1.0D0/(long_nfactor(5))

       elseif (num_well .eq. 10) then

        long_nfactor(1)=( 0.0785047D0*rnmres) /(rnmres*0.245032D0 + 1.0D0)
        long_nfactor(2)=( 0.0152761D0*rnmres) /(rnmres*0.0283803D0 + 1.0D0)
        long_nfactor(3)=( 0.205481D0*rnmres) /(rnmres*0.385813D0  + 1.0D0)
        long_nfactor(4)=( 0.0174765D0*rnmres) /(rnmres*0.0174638D0 + 1.0D0)
        long_nfactor(5)=( 0.0352685D0*rnmres) /(rnmres*0.0269838D0 + 1.0D0)
        long_nfactor(6)=( 0.0474026D0*rnmres) /(rnmres*0.0249249D0 + 1.0D0)
        long_nfactor(7)=( 0.18665D0 *rnmres) /(rnmres*0.0107983D0 + 1.0D0)
        long_nfactor(8)=( 0.0390303D0*rnmres) /(rnmres*0.0140943D0 + 1.0D0)
        long_nfactor(9)=( 0.0327411D0*rnmres) /(rnmres*0.00812347D0 + 1.0D0)
        long_nfactor(10)=( 0.0561461D0*rnmres) /(rnmres*0.00743991D0 + 1.0D0)

           long_nfactor(1)=1.0D0/(long_nfactor(1))
           long_nfactor(2)=1.0D0/(long_nfactor(2))
           long_nfactor(3)=1.0D0/(long_nfactor(3))
           long_nfactor(4)=1.0D0/(long_nfactor(4))
           long_nfactor(5)=1.0D0/(long_nfactor(5))
           long_nfactor(6)=1.0D0/(long_nfactor(6))
           long_nfactor(7)=1.0D0/(long_nfactor(7))
           long_nfactor(8)=1.0D0/(long_nfactor(8))
           long_nfactor(9)=1.0D0/(long_nfactor(9))
           long_nfactor(10)=1.0D0/(long_nfactor(10))

        else
            write(SO,*) 'num_well problem for ab_c_of_n_new',num_well
            stop 
        endif
       endif

       if (ab_c_of_n_old) then
       rnmres=real(nmres)
       long_nfactor(1)=1.0D0/(rnmres*0.0015D0 + 1.94D0)
       long_nfactor(2)=1.0D0/(rnmres*0.0032D0 + 1.83D0)
       long_nfactor(3)=1.0D0/(rnmres*0.022D0 + 7.77D0)

       if (num_well.ne.3) then 
         write(SO,*) 'numwell must be 3 for ab_c_of_n_*old*',num_well
         stop
       endif

       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize potential and force table

        do 514 i514=1,indxt
          do 515 i515=0,maxs2+1
            vpotnt(i515,i514,curtab)=0.0D0
            forse(i515,i514,curtab)=0.0D0
  515     continue
  514   continue

      idigns=.false.
      if( passi )then
         passt=.true.
      else
         passt=.false.
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   read target (which is considered 0th memory)
         read (iprolst,1000)tarfl
1000     format(a5)
c         profile='/home/mprentis/amh/proteins/'//tarfl
CCC           profile='/home/mprentis/amh/proteins/'//tarfl

         protnm(0)=tarfl
         matfile='match/'//trim(tarfl)//'/'//tarfl

          open(imemri,file='proteins/'//tarfl,
     *                          status='old',iostat=open_status)

         if (open_status.ne.0) then
           write(6,*) 'failure to open profile 0th mem file '
           write(6,*) 'error number ',open_status
           stop
         endif

c        read in memory protein coordinates, and primary sequence and crystal
c        secondary structure

!         write(6,*)'call getmem '

         call getmem(protnm(0),nmrss(0),maxres,jres,imemri,numcrd,ywork,
     *               mempre(1,curtab),oarchv,sa,iwork,passt,0)

         close(imemri)

         work8(1)=float(nmrss(0))

c        set hydrophobicity profile

         call hprofl(maxres,nmrss(0),jres,hydseq(1,1,curtab),hydscl(0,1))
         call hprofl(maxres,nmrss(0),jres,hydseq(1,2,curtab),hydscl(0,2))

c   place coordinates in target
c       and sequence profile in tres

               call pttarg(maxsiz,nmrss(0),numcrd,target,maxres,ywork,ires,jres,
     *                     shydro(1,1,curtab),hydseq(1,1,curtab),tarpre(1,curtab),
     *                     mempre(1,curtab),bondln,oarchv,passi)


c  end of reading target info
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     contact order term

       if (i_contact_order.and.curtab.eq.4) then
         do i_ixn=1,indxt
           id3=ilong(i_ixn,1,curtab)
           id4=ilong(i_ixn,2,curtab)
           do i_contact_term=1,n_contact_order_terms
           if (  (id4-id3.ge.i_contact_order_min(i_contact_term)) .and.
     *        (id4-id3.le.i_contact_order_max(i_contact_term)) .and.
     *                    ires(id3).ne.8 .and.
     *                        ires(id4).ne.8 )  then
              do i_r=1,maxsz
                 r_temp=rinc(i_ixn,curtab)*real(i_r)
          t_min=tanh(7.0D0*(r_temp-r_min_contact_order(i_contact_term)))
          t_max=tanh(7.0D0*(r_max_contact_order(i_contact_term)-r_temp))
                 theta = (gamma_contact_order(i_contact_term,1)
     *              +gamma_contact_order(i_contact_term,2)*
     *           ((id4-id3)**gamma_contact_order(i_contact_term,3)) )*
     *              0.25D0*( 1.0D0+t_min )*( 1.0+t_max )
                 theta_dot=7.0D0*theta*(t_max-t_min)

                 vpotnt(i_r,i_ixn,curtab)=vpotnt(i_r,i_ixn,curtab)+theta
                 forse(i_r,i_ixn,curtab)=forse(i_r,i_ixn,curtab)-theta_dot/r_temp
              enddo

           endif
           enddo
         enddo
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read in sequences to average force over, if avg_seq is on
c      if (ave_seq) then
c          write(6,*)'target_sequences'
c        open(itarg_seq,file='target_sequences',status='old')
c           read(itarg_seq,*)numseq
c           write(6,*)'numseq'
c        do iseq = 1,numseq
c          read(itarg_seq,999)(tgsequences(id3,iseq),id3=1,nmres)
c          write(6,999)(tgsequences(id3,iseq),id3=1,nmres)
c999       format(25(i2,1x))
c        enddo
c        close(itarg_seq)
c      else
c        numseq = 1
c!        tgsequences(1:nmres,1)=ires(1:nmres)
c      endif
c          if((ave_seq) .and. (numseq.gt.maxseq) ) then
c            write(6,*) 'numseq greater than maxseq'
c            stop
c          endif



         if (targ_cons) then
       confile='/home/mprentis/amh/md_input/targcons/'//tarfl
c               123456789012345678901234567890123456789012345678
          write(6,*)'targ_cons directory = ' , confile
          open(icon,file=confile,status='old',iostat=open_status)
               if (open_status.ne.0) then
                 write(6,*) 'failure to open file ',profile
                 write(6,*) 'error number ',open_status
                 stop
               endif
               read(icon,*)blah
               read(icon,102)connmres
               read(icon,999)(tgsequences(i1,numseq),i1=1,nmres)
999            format(25(i2,1x))
              write(6,*)'targ_cons sequence'
               write(6,999)(tgsequences(i1,numseq),i1=1,nmres)
                                                                                
102            format (i5)
                write(6,*)'connmres =',connmres
                if( connmres .gt.maxres )then
C                   write(oarchv,611)protnm,nmrss,maxres
611                format('targ_cons ',a5,' too large ',i4,
     *                   ' reserved space ',i4)
                     write(6,611)protnm,nmrss,maxres
                     stop
                endif
                if (connmres .ne. nmres )then
c                   write(oarchv,612)connmres,nmres
612                format('connmres',i5,' .ne. nmres ',i5,' 612 gaspot')
                   write(6,612)connmres,nmres
                   stop
                endif
          close(icon)
        endif ! if (targ_cons)

         if (curtab.eq.1) then
c                   write(6,*)'target  sequence   ',targ_cons
c               do 909 a2=1,numseq
c                    write(6,999)(tgsequences(i1,a2),i1=1,nmres)
c909            continue
         endif  
                                                                          
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     loop over each memory

      do 509 i_mem=1,nummem
         read (iprolst,1000)profl
                                                                                
         protnm(i_mem)=profl
CCC         matfile='/home/mprentis/amh/match/'//trim(tarfl)//'/'//profl
c insert targ_cons seq part here for memories
c  /home/mprentis/amh/consens_seq

c           write(6,*) 'matfile ',matfile
c           write(6,*) 'profile ',profile 

          open(imemri,file='proteins/'//profl,status='old',iostat=open_status)

         if (open_status.ne.0) then
           write(SO,*) ' proteins ',profl 
           write(SO,*) 'failure to open protein 543 file '
           write(SO,*) 'error number ',open_status
           stop
         endif
 
c        read in memory protein coordinates, and primary sequence and crystal 
c        secondary structure

         call getmem(protnm(i_mem),nmrss(i_mem),maxres,jres,imemri,numcrd,ywork,
     *               mempre(1,curtab),oarchv,sa,iwork,passt,i_mem)

         close(imemri)
 
         work8(i_mem+1)=float(nmrss(i_mem))

c        set hydrophobicity profile

         call hprofl(maxres,nmrss(i_mem),jres,hydseq(1,1,curtab),hydscl(0,1))
         call hprofl(maxres,nmrss(i_mem),jres,hydseq(1,2,curtab),hydscl(0,2))

c        read in and label secondary structure based on data bank lables

c000000000000000000000000000000000000000000000000000000000000000000


         if( passf.and.(i_mem.eq.1) )then
c            write(oarchv,905)protnm(0),protnm(i_mem)
  905       format(/'2nd structures for ',a5,2x,a5)
            do 544 i544=1,nmrss(i_mem)
               write(oarchv,904)i544,aminoa(jres(i544)),
     *                          (tarpre(i544,i1),
     *                          (int(shydro(i544,i2,i1)),
     *                           i2=1,2),i1=1,numtab)
  904          format(i4,1x,a3,2x,4(3(i2,1x),1x))
  544       continue
         endif

cMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
c  Read in match file
c
c          read in alignments
           do 670 i1=1,maxres
             match(i1)=0
670           continue
         open(imat,file='match/'//trim(tarfl)//'/'//profl,
     *                     status='old',iostat=open_status)

         if (open_status.ne.0) then
           write(SO,*) 'match ', profl 
           write(SO,*) 'failure to open file match 594 '
           write(SO,*) 'error number ',open_status
           stop
         endif

           read (imat,*)
           read (imat,1002)(match(i_res), i_res=1,nmres)
1002       format(10(i4))
          close (imat)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (ran_force .and. curtab.eq.1) then
           if (ran_file) then
             open(iran,file='random_in',status='old')
             read(iran,*) idummy
              if (idummy .ne. indxt) then
                write(SO,*) 'wrong no of interactions in random_in'
                write(SO,*) idummy,indxt
                stop
              endif
              read(iran,*) (ran_num(i), i=1,indxt)
           else
             temp_mcp_pass=3
             call SLARNV(temp_mcp_pass,iseed_amh(1),indxt,ran_num) !should be gaussian dist mean 0, sd 1  
             do j=1,indxt
               ran_num(j) = ran_num(j)*s_ran + lambdaR
               if (ran_num(j).lt.0.0D0.and.(.not.allow_neg)) ran_num(j)=0.0D0
             enddo
           endif
           open(oran,file='random_out',status='new')
             write(oran,*) indxt
             write(oran,*) (ran_num(i), i=1,indxt)
           close(oran)
         endif


c        loop over each constraint

         do 519 i_ixn=1,indxt

c            set r-grid
               
             do 516 i516=1,maxsz
               rgrid(i516)=rinc(i_ixn,curtab)*float(i516)
516          continue

c      id3, id4 are the i and j indices (respectively) for the target protein

            id3=ilong(i_ixn,1,curtab)
            id4=ilong(i_ixn,2,curtab)

c      id1, id2 are the i and j indices (respectively) in the memory protein
c      that id3 and id4 are aligned to

            id1=match(id3)
            id2=match(id4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c use memory secondary structural info: currently in some encodings the interactions
c are given a different weight (ie gamma value) if *either* one of the residues
c involved is aligned to a helical residue in the memory
c I made this bit of the code a bit more long-winded than before to avoid looking up
c mempre(id1,curtab) when id1=0 (non-aligned) because this is outside the array.

              if (id1.eq.0) then  !first residue not aligned => not aligned to helix
                helix_A=.false.
              else
                if (mempre(id1,curtab) .eq. 1) then 
                   helix_A=.true.    !first residue aligned to helix
                else
                   helix_A=.false.    !first residue *NOT* aligned to helix
                endif
              endif
              if (id2.eq.0) then  !second residue not aligned => not aligned to helix
                helix_B=.false.
              else
                if (mempre(id2,curtab) .eq. 1) then 
                   helix_B=.true.    !second residue aligned to helix
                else
                   helix_B=.false.    !second residue *NOT* aligned to helix
                endif
              endif
              
              if (helix_A.or.helix_B) then
                 helix = 1
              else
                 helix = 2
              endif

              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc             check if match file corrupted   ccccccccc

              if( (id1.gt.nmrss(i_mem)).or.(id1.lt.0) ) then
                  write(SO,*) 'match file corrupted'
                  write(SO,*) 'memory number',i_mem,'name',protnm(i_mem),'residue',id3
                  write(SO,*) 'id1=',id1
                  write(SO,*) 'id3=',id3
                  write(SO,*) (match(i),i=1,nmres)
                  stop
              endif
              if( (id2.gt.nmrss(i_mem)).or.(id2.lt.0).or.
     *              ( (id2-id1.le.0) .and. (id2*id1.ne.0) ) )then
                  write(SO,*) 'match file corrupted'
                  write(SO,*) 'memory number',i_mem,'name',protnm(i_mem),'residue',id4
                  write(SO,*) 'id2=',id2
                  write(SO,*) 'id1=',id1
                  stop
              endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


            if ( ires(id3).eq.8 .and. ires(id4).eq.8 ) then
               tab_for_contact=1
            elseif ( ires(id3).eq.8 .and. ires(id4).ne.8 ) then
               tab_for_contact=2
            elseif ( ires(id3).ne.8 .and. ires(id4).eq.8 ) then
               tab_for_contact=3
            else
               tab_for_contact=4
            endif                 !decide table for contact potential   
                                  !--usually beta-beta (ie 4) but must
                                  ! allow for glycines
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc reasons to exit loop over ints (consider first only if not contact-int)

            if ( .not. ( id4-id3.gt.12 .and. i_mem.eq.1 .and.
     *                  curtab.eq.tab_for_contact .and. n_letters.eq.4
     *                  .and. i_3rd_is_contact ) ) then


!             only interactions between residues of certain separation
              if ((id4-id3).lt.min_seq_sep) goto 519

c             no potential for beta-carbon on glycine 
               if( (cdtyp1.eq.2).and.(ires(id3).eq.8) )go to 519
               if( (cdtyp2.eq.2).and.(ires(id4).eq.8) )go to 519


!              if not aligned, go to next interaction
              if ( (id1.eq.0) .or. (id2.eq.0) ) goto 519

!             if align to glycine and on beta carbon then skip
              if( (cdtyp1.eq.2).and.(jres(id1).eq.8) )go to 519
              if( (cdtyp2.eq.2).and.(jres(id2).eq.8) )go to 519

! exit loop if contact interaction and have a glycine and gly_con is false
            elseif (.not.gly_con) then 
              if ( (ires(id3).eq.8).or.(ires(id4).eq.8) ) goto 519
            endif

ccccccccccc      end of reasons to exit loop over ints
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


               dist_ij = dsqrt((ywork(id1,1,cdtyp1) -
     *                            ywork(id2,1,cdtyp2))**2 +
     *                           (ywork(id1,2,cdtyp1) -
     *                            ywork(id2,2,cdtyp2))**2 +
     *                           (ywork(id1,3,cdtyp1) -
     *                            ywork(id2,3,cdtyp2))**2 )

              dist_prox = dsqrt((ywork(id1,1,1) -
     *                            ywork(id2,1,1))**2 +
     *                           (ywork(id1,2,1) -
     *                            ywork(id2,2,1))**2 +
     *                           (ywork(id1,3,1) -
     *                            ywork(id2,3,1))**2 ) 

c              if ( (curtab.eq.4).and.(i_ixn.eq.101) ) then
c                write(SO,*) 'distance is',dist_ij
c                write(SO,*) 'cdtyp1 and 2 are',cdtyp1,cdtyp2
c              endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set category of interaction based on distance in seq and in space.

              if (.not.i_alt_prox) then

                   if (dist_ij.lt.srcut) then
                     if ((id4-id3).lt.minmr) then
                       iprox=1
                     else
                       iprox=2
                     endif
                   else
                     iprox=3
                   endif

c      Mike changed this so short range in seq pairs actually
c      go into the long range in seq and space category if they
c      are separated by more than srcut Angstroms
c      note also, that I appear to set this for each table, while
c      corey choses class based only on alpha-alpha dist (+prox)
c      This is probably not a big deal, but when running using
c      his gammas, switch i_alt_prox flag to true and will get 
c      the following proximity classes


              elseif (.not.i_3rd_is_contact) then
 
                   
                   iprox = 4
                   if ((id4-id3).lt.minmr) iprox=1
                   if ( ((id4-id3).ge.minmr)
     *             .and. ((id4-id3).le.maxmr)) iprox=2
                   if ( ((id4-id3) .gt. maxmr) .and.
     *             (dist_prox .lt. alt_prox_cut) ) iprox = 3

              else

                   iprox = 4
                   if ((id4-id3).lt.minmr) iprox=1
                   if ( ((id4-id3).ge.minmr)
     *             .and. ((id4-id3).le.maxmr)) iprox=2
                   if ( ((id4-id3) .gt. maxmr).and.(i_mem.eq.1)
     *               .and. curtab.eq.tab_for_contact ) iprox = 3

              endif

c      random interactions are for sequence separations of
c      ran_min_seq_dist and above

              if ((id4-id3).lt.ran_min_seq_dist) then
                iprox_ran=.false.
              else
                iprox_ran=.true.
              endif   

              if ( i_alt_prox.and.(iprox.eq.4) ) goto 519



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               find individual charges and total charge contribution
c     used to be done in subroutine charge, but now done here
c     because so short 

               totchg=0.0D0
          if (i_3rd_is_contact.and.iprox.eq.3) then
           do i_well=1,num_well
            do iseq=1,numseq_amc
              totchg(i_well)= totchg(i_well)-
     *        qchrg(tgsequences_amc(id3,iseq),tgsequences_amc(id4,iseq),
     *        jres(id1),jres(id2),i_well+2)
                   enddo
              totchg(i_well)=totchg(i_well)*
     *                          long_nfactor(i_well)
                 enddo
               else
              if (four_plus .and. iprox .le. 2) then
          do iseq=1,numseq_amc
          totchg(1)=totchg(1)-
     *     qchrg2(tgsequences_amc(id3,iseq),tgsequences_amc(id4,iseq),jres(id1),jres(id2),helix,iprox)
          enddo
                 else
          do iseq=1,numseq_amc
          totchg(1)=totchg(1)-
     *     qchrg(tgsequences_amc(id3,iseq),tgsequences_amc(id4,iseq),jres(id1),jres(id2),iprox)
          enddo
                 endif
                 do i_well=2,num_well
                   totchg(i_well)=0.0D0
                 enddo
               endif

	       if (dist_cut) then
                if ( (iprox .eq. 2) .and.(dist_ij .gt. 6.5)) then
                  totchg=0.0D0
                endif
                endif

               totchg=totchg/real(numseq_amc)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              find r-dependence of gaussian potential
               delt_safe=deltz(i_ixn,curtab)
               if( welscl )then
                  call gaussv(maxsz,dist_ij,deltz(i_ixn,curtab),gaussp,rgrid)
         if ( (go_con) .and. (dist_ij .gt. go_con_dist)) then
                   gaussp=0.0D0
          endif

               else
                  call gaussw(maxsz,dist_ij,0.25D0*dist_ij,gaussp,rgrid)
               endif

               if (ran_force.and. (curtab.eq.1).and.(iprox_ran)) then
                  call gaussv(maxsz,r_ran,delt_safe,gaussq,rgrid)
               endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              combine sequence-structure terms

               do 507 i507=1,maxsz

                 if ( i_3rd_is_contact .and.(iprox .eq. 3) ) then
                   gaussp(i507) = 0.0D0
                   do i_well=1,num_well
                     gaussp(i507)= gaussp(i507)+totchg(i_well)*
     *               0.25D0*((1.0D0 + tanh(7.0D0*(rgrid(i507)-r_min(i_well) )))
     *               *(1.0D0+tanh(7.0D0*(r_max(i_well)-rgrid(i507)))))
                   enddo
                 else
                   gaussp(i507)=totchg(1)*gaussp(i507)
                 endif
 
  507          continue

c                 set up potential
     
                    do 522 i522=1,maxsz
                        vpotnt(i522,i_ixn,curtab)=vpotnt(i522,i_ixn,curtab) + gaussp(i522)
  522               continue

                  if ((curtab.eq.1) .and. ran_force.and. 
     *                                   (iprox_ran)) then
                    do 588 i588=1,maxs
                        gaussq(i588)=gaussq(i588)*ran_num(i_ixn)
                        vpotnt(i588,i_ixn,1)=vpotnt(i588,i_ixn,1) +gaussq(i588)
 588               enddo
                  endif                       

c                 include r-dependent force term=
c                 -(1/r)*dV(r)/dr

                  delt2=2.0D0*deltz(i_ixn,curtab)
                  do 504 i504=1,maxsz

                    if (i_3rd_is_contact.and.iprox .eq. 3) then

                      gaussp(i504)=0.0D0
                      do i_well=1,num_well
                        force_term(1)=
     *                    tanh(7.0D0*(r_max(i_well)-rgrid(i504)))**2  ! G
                        force_term(2)=
     *                    -tanh(7.0D0*(rgrid(i504)-r_min(i_well)))**2   ! U
                        force_term(3)= 
     *                    tanh(7.0D0*(r_max(i_well)-rgrid(i504)))    !G
                        force_term(4)= 
     *                    -tanh(7.0D0*(r_max(i_well)-rgrid(i504)))*    !GU^2
     *                      tanh(7.0D0*(rgrid(i504)-r_min(i_well)))**2
                        force_term(5)= 
     *                    -tanh(7.0D0*(rgrid(i504)-r_min(i_well)))    !U
                        force_term(6)= 
     *                    tanh(7.0D0*(rgrid(i504)-r_min(i_well)))*  !UG^2
     *                    tanh(7.0D0*(r_max(i_well)-rgrid(i504)))**2


                       gaussp(i504)=gaussp(i504)
     *                      -7.0D0*(totchg(i_well)/4.0)* (
     *                      force_term(1) + force_term(2)
     *                      + force_term(3) + force_term(4) +
     *                      force_term(5) + force_term(6) )/rgrid(i504)
                      enddo

                    else
                      gaussp(i504)=delt2*gaussp(i504)*
     *                         (1.0D0 - (dist_ij/
     *                          rgrid(i504)) )

                    endif

  504             continue

c                 set up force term
    
                    do 523 i523=1,maxsz
                      forse(i523,i_ixn,curtab)=forse(i523,i_ixn,curtab) + gaussp(i523)
  523               continue

                  if (curtab.eq.1 .and. ran_force
     *                          .and.(iprox_ran)) then
                    do 589 i589=1,maxs
                      gaussq(i589)=delt2*gaussq(i589)*
     *                           (1.0D0 - (r_ran/rgrid(i589)) )

                        forse(i589,i_ixn,1)=forse(i589,i_ixn,1) + gaussq(i589)
                      
 589                enddo
                  endif                       

c            end loop over insertions and deletions

  519    continue                        ! End of loop over interactions
  509 continue                                ! End of loop over memories

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         prints test potential + force (+integrated force)
 
       if (i_V_test) then

         i_test_tab=curtab
         i_test_ixn=ixn_from_site(test_site(1),test_site(2),i_test_tab)
         if ( (test_site(1).ne.ilong(i_test_ixn,1,i_test_tab))
     *    .or.(test_site(2).ne.ilong(i_test_ixn,2,i_test_tab)) )
     *      then
            write(SO,*) 'error in gaspot'
            write(SO,*) 'test sites for table',i_test_tab,'wrong'
            write(SO,*) 'should be',test_site(1),test_site(2)
            write(SO,*) 'find from looking up interaction they are'
            write(SO,*) ilong(i_test_ixn,1,i_test_tab),
     *                    ilong(i_test_ixn,1,i_test_tab)
            goto 1111
          endif

         vpotint_test=0.0D0

         count=test_site(1)
         call num_to_char(count,ccount1)
         count=test_site(2)
         call num_to_char(count,ccount2)
         count=i_test_tab
         call num_to_char(count,ccount3)
         open(ovgaspot,file='V_gaspot_'//trim(ccount1)//'_'//
     *    trim(ccount2)//'.'//trim(ccount3),status='unknown',
     *    action='write')

       write(ovgaspot,*) 'r,vpotnt,forse, int force for ixn num,table'
     *                                  ,i_test_ixn,i_test_tab
         write(ovgaspot,*) 'sites',ilong(i_test_ixn,1,i_test_tab)
     *                               ,ilong(i_test_ixn,2,i_test_tab)
         do i=maxs,1,-1
           vpotint_test(i)=vpotint_test(i+1)+
     *                         (1.0D0*rinc(i_test_ixn,i_test_tab))*0.5D0*
     *                  (forse(i+1,i_test_ixn,i_test_tab)*float(i+1)+
     *                    forse(i,i_test_ixn,i_test_tab)*
     *                      float(i))*rinc(i_test_ixn,i_test_tab)
         enddo
         do i=1,maxs+1
           write(ovgaspot,*) float(i)*rinc(i_test_ixn,i_test_tab),
     *                   vpotnt(i,i_test_ixn,i_test_tab),
     *     (forse(i,i_test_ixn,i_test_tab)*rinc(i_test_ixn,i_test_tab))
     *                          *float(i),vpotint_test(i)
         enddo    
         close(ovgaspot)
1111   endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do i_ixn=1,indxt
         vpotnt(maxs+1,i_ixn,curtab)=0.0D0


       do i=maxs,0,-1

         m=(1.0D0/rinc(i_ixn,curtab))*(forse(i+1,i_ixn,curtab)-forse(i,i_ixn,curtab))

         c=forse(i,i_ixn,curtab)*float(i+1)-forse(i+1,i_ixn,curtab)*float(i)

         vpotnt(i,i_ixn,curtab)=vpotnt(i+1,i_ixn,curtab)+
     *     m*rinc(i_ixn,curtab)*rinc(i_ixn,curtab)*rinc(i_ixn,curtab)*
     *      (float( (i+1)*(i+1)-i)-2.0D0/3.0D0)  +
     *     c*rinc(i_ixn,curtab)*rinc(i_ixn,curtab)*(float(i)+0.5D0)

       enddo
       enddo

c     write out proteins used in generating potential

      if( passi )then
c         write(oarchv,177)protnm(0),int(work8(1))
  177    format(/'target protein ',a5,1x,i5)
c         write(oarchv,111)nummem
  111    format('# proteins used in generating potential ',i3)
c         write(oarchv,114)(protnm(i1),int(work8(i1+1)),i1=1,nummem)
  114    format(5(a5,1x,i5,2x))
      endif

c     invert r-grid increment

      do 518 i518=1,indxt
         rincinv(i518,curtab)=1.0D0/rinc(i518,curtab)
         rincsq(i518,curtab)=rinc(i518,curtab)*
     *                              rinc(i518,curtab)
  518 continue

c     --- diagnostics ---


         do 534 i_ixn=1,indxt

            id1=ilong(i_ixn,1,curtab)
            id2=ilong(i_ixn,2,curtab)

c           check if endpoints of force small enough

c            if( abs(forse(maxs,i_ixn,curtab)).gt.5.0e-5 )
c     *         write(oarchv,121)curtab,i_ixn,id1,id2,
c     *               forse(maxs,i_ixn,curtab)
c  121          format('F > 5.0e-05:table ',i2,' ixn ',i5,
c     *                ' between ',2(i3,1x),' F(maxs) ',
c     *                1pe10.3)

            if( (id1.eq.0).and.(id2.le.12).and.(cdtyp1.eq.2)  )then

c              check if potential has been correctly
c              generated

c              generate r-grid

               do 503 i503=1,maxs
                 rgrid(i503)=float(i503)*rinc(i_ixn,curtab)
  503          continue

c               write(oarchv,110)curtab,cdtyp1,cdtyp2,id1,id2,
c     *                          rincinv(i_ixn,curtab)
  110          format(/'test force: table number ',i3,
     *                ' coord1 ',i1,' coord2 ',i1,
     *                ' sites',2(1x,i3),' rincinv ',1pe10.3)
               do 532 i532=2,maxs-1

                  work2(i532)=0.5D0*(vpotnt(i532-1,i_ixn,curtab)-
     *                             vpotnt(i532+1,i_ixn,curtab))/
     *                            (rgrid(i532-1)-rgrid(i532))

                  if( forse(i532,i_ixn,curtab)*rgrid(i532)
     *                    .ne.0.0D0 )
     *               rnorm=-work2(i532)/
     *                  (forse(i532,i_ixn,curtab)*rgrid(i532))

  112             format(i3,5(1x,1pe10.3))
  532          continue

            endif

  534    continue

c     --- end diagnostics ---

          write(SO,*) 'end of gaspot '
c     ---------------------- done -----------------------

      return
      end
