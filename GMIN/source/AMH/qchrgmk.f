c     --------------------- qchrgmk ----------------------

      subroutine qchrgmk(qchrg,qchrg2,oarchv,n_letters)

c     ---------------------------------------------------
c     QCHRGMK makes table qchrg
c     ---------------------------------------------------

      use amhglobals,  only: SO,max_letters,i_3rd_is_contact,four_plus,
     *  i_non_add_contact,sort_non_add,
     *  gamma_non_add,ngamma_non_add,class_2,max_well,num_well,
     *  para_HB,anti_HB,anti_NHB,i_ignore_chain_dirn,mismatch,
     *  para_one,anti_one,n_letters_con, srplr,igamma,ihbond

      use globals_alt, only:altpotflag ! For alternative potential

      implicit none

c
c	i_ranges	3; near in seq.   -->     near in space
c                          far  in seq.   and     near in space
c                          far  in seq.   and     far  in space

         integer n_letters

c     internal variables:
                                                                                                     
         double precision qchrg(21,21,20,20,max_well+2),optgamma(512),qchrg2(21,21,20,20,2,2),gamma_temp

         integer i500,i501,i502,i503,i504,i_ranges,ngamma,ngamma_check,open_status,
     *           read_status,isit1,isit2,n_gammas_per_con_well,ic1,oarchv                               

c     tshen
      integer myj,myip,myjp,myk,mykp,mywel
      double precision mychch, mychpo

        integer class(20),class_con(20),icl,jcl,ipcl,jpcl,
     *       gammapt,sort_2(2,2,2,2),class_4(20),class_20(20),
     *       sort(max_letters,max_letters,max_letters,max_letters),
     *       gammapt_contact,sort_contact(max_letters,max_letters),
     *       n_ranges,n_gammas_in_range(max_well+2),i_well,i_well2,
     *       i_type,icl_con,jcl_con,
     *  sort_sym(max_letters,max_letters,max_letters,max_letters),
     *  sort_plus(max_letters,max_letters,max_letters,max_letters,2)

         integer i, j, ip, jp,count

        data class_4  /1,3,2,2,4,2,2,1,3,4,4,3,4,4,1,1,1,4,4,4/
        data class_20 /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
        data sort_2 /1,5,6,16,7,3,14,12,8,15,2,11,13,10,9,4/
        class_2=(/2,1,1,1,2,1,1,2,1,2,2,1,2,2,1,1,1,2,2,2/)

c  set classes for different types of code
c  due to history, gammas are grouped differently
c  for the two-letter code compared to multi-letter
c  codes. The two-letter code is thus treated on a
c  slightly different footing.

c          write(6,*)'n_letters ',n_letters 

          if (n_letters.eq.2) then
            do i=1,20                  !two-letter code
              class(i)=class_2(i)
            enddo
          elseif (n_letters.eq.4) then
            do i=1,20                  !four-letter code
              class(i)=class_4(i)
            enddo
          endif
                                                                                                     
         if (n_letters_con.eq.2) then  !classes may be different for contact interactions
             class_con(1:20)=class_2(1:20)
          elseif (n_letters_con.eq.4) then
             class_con(1:20)=class_4(1:20)
          elseif (n_letters_con.eq.20) then
             class_con(1:20)=class_20(1:20)
          endif
        if (n_letters.eq.2) then

          do i=1,2                  ! two-letter code
          do j=1,2
          do ip=1,2
          do jp=1,2
          sort(i,j,ip,jp)=sort_2(i,j,ip,jp)
          enddo
          enddo
          enddo
          enddo
 
        else

c        sort_sym is for the case where four_plus is false and chain direction is ignored
                                                                                                     
          count = 1
             do i = 1,n_letters
             do j = i,n_letters   ! note j starts at i (and jp at ip) so looping over all
             do ip = 1,n_letters  ! distinguishable pairs in both target and memory
             do jp = ip,n_letters  !(a total of n_letters+(n_letters -choose-2) in each
               sort_sym(i,j,ip,jp) = count
               sort_sym(j,i,jp,ip) = count  ! swapping both pairs (ie (ij-i'j' -> ji-j'i')
                                            ! gives same gamma if chain direction considered irrelevant
               if ( (i.ne.j) .and. (ip.ne.jp) ) count=count+1 !if both pairs different then swapping identities in just
               sort_sym(j,i,ip,jp) = count            !one pair will give a new gamma, so increment
               sort_sym(i,j,jp,ip) = count
               count = count +1
             enddo
             enddo
             enddo
             enddo


          count = 1                 !multi-letter code
          do 2 i = 1,n_letters
          do 4 j = 1,n_letters
          do 6 ip = 1,n_letters
          do 8 jp = 1,n_letters

          sort(i,j,ip,jp) = count
          count = count +1
8         continue
6         continue
4         continue
2         continue

           if (four_plus) then           !include secondary structural information

           if (i_ignore_chain_dirn.and.(.not.mismatch)) then  ! here consider the order of residues in chain
                                                              !  not to effect their interaction energy
             count = 1      
             do i = 1,n_letters
             do j = i,n_letters   ! note j starts at i (and jp at ip) so looping over all 
             do ip = 1,n_letters  ! distinguishable pairs in both target and memory
             do jp = ip,n_letters  !(a total of n_letters+(n_letters -choose-2) in each 
             do ic1  = 1,2
               sort_plus(i,j,ip,jp,ic1) = count
               sort_plus(j,i,jp,ip,ic1) = count  ! swapping both pairs (ie (ij-i'j' -> ji-j'i') 
                                                    ! gives same gamma if chain direction considered irrelevant
               if ( (i.ne.j) .and. (ip.ne.jp) ) count=count+1 !if both pairs different then swapping identities in just
               sort_plus(j,i,ip,jp,ic1) = count            !one pair will give a new gamma, so increment
               sort_plus(i,j,jp,ip,ic1) = count
               count = count +1
             enddo
             enddo
             enddo
             enddo
             enddo


           elseif (mismatch.and.i_ignore_chain_dirn) then !here ignore not only chain direction but also
                                                          !assume pairs that for each kind of target pair have only two
                                                          !kinds of alignment: exactly right and not exactly right (ie mismatched)
                 sort_plus=0
                 count = 1
                 do i = 1,n_letters
                 do j = i,n_letters
                 do ic1 = 1,2
                 sort_plus(i,j,i,j,ic1)=count
                 sort_plus(j,i,j,i,ic1)=count
                 count = count + 1
                 enddo
                 enddo
                 enddo

                 count = 21

                 do i = 1,n_letters
                 do j = i,n_letters
                 do ic1 = 1,2
                        do ip = 1,4
                        do jp = 1,4
                        if (ip .ne. i .or. jp .ne. j) then
                         sort_plus(i,j,ip,jp,ic1)=count
                         sort_plus(j,i,jp,ip,ic1)=count
                        endif
                        enddo
                        enddo
                  count = count + 1
                  enddo
                  enddo
                  enddo

           else

          count = 1                 !multi-letter code
          do i = 1,n_letters
          do j = 1,n_letters
          do ip = 1,n_letters
          do jp = 1,n_letters
          do ic1  = 1,2
          sort_plus(i,j,ip,jp,ic1) = count
          count = count +1
          enddo
          enddo
          enddo
          enddo
          enddo

           endif
           endif  ! four_plus

          if (i_3rd_is_contact) then
            count = 1
            do i=1,n_letters_con
            do j=i,n_letters_con
            sort_contact(i,j)=count 
            sort_contact(j,i)=count
            count=count+1
            enddo
            enddo
          endif

        endif
        
        if (i_non_add_contact) then
          count=1
          do i=1,2
          do j=1,2
          do i_well=1,2
          do i_type=1,2
          do i_well2=1,2
            sort_non_add(i,j,i_well,i_type,i_well2)= count
            count=count+1
          enddo
          enddo
          enddo
          enddo
          enddo
        endif
 
       n_gammas_per_con_well=n_letters_con*(n_letters_con+1)/2

        if (n_letters.eq.2) then
          ngamma=48
        elseif ((n_letters.eq.4).and.i_3rd_is_contact) then
          if (four_plus) then
            if (i_ignore_chain_dirn.and.(.not.mismatch)) then
              ngamma=544+num_well*n_gammas_per_con_well
            elseif (i_ignore_chain_dirn.and.mismatch) then
              ngamma=80+num_well*n_gammas_per_con_well
            else
              ngamma=1024+num_well*n_gammas_per_con_well
            endif
          else
            if (i_ignore_chain_dirn.and.(.not.mismatch)) then
            ngamma=272+num_well*n_gammas_per_con_well
            else
            ngamma=512+num_well*n_gammas_per_con_well
            endif
          endif
        elseif ((n_letters.eq.4).and.(.not.i_3rd_is_contact)) then
          ngamma=768
        endif

        if (i_non_add_contact) ngamma=ngamma+ngamma_non_add
    
        open(igamma,file='gamma.dat',action='read',status='old',iostat=open_status)
        if (open_status.ne.0) then
          write(SO,*) 'failure to open gamma file'
          write(SO,*) 'error number ',open_status
          stop
        endif

        read(igamma,*) ngamma_check
        if (ngamma_check.ne.ngamma) then
          write(SO,*) 'gamma cock-up',ngamma_check,ngamma
          stop
        endif

           write (oarchv,*)
           write (oarchv,*) 'Gammas Used'
           write (oarchv,*) '-----------'


        if (i_3rd_is_contact.and.(n_letters.eq.4)) then
           n_ranges=2+num_well
           if (four_plus) then
             if (i_ignore_chain_dirn.and.(.not.mismatch)) then
               n_gammas_in_range(1)=272
               n_gammas_in_range(2)=272
             elseif (i_ignore_chain_dirn.and.mismatch) then
               n_gammas_in_range(1)=40
               n_gammas_in_range(2)=40
             else
               n_gammas_in_range(1)=512
               n_gammas_in_range(2)=512
             endif
           else

             if (i_ignore_chain_dirn) then
               n_gammas_in_range(1)= 136
               n_gammas_in_range(2)= 136
             else
               n_gammas_in_range(1)= 256
               n_gammas_in_range(2)= 256
             endif

           endif

          do i = 1,num_well
             n_gammas_in_range(2+i)=n_gammas_per_con_well
           enddo

        else 
           n_ranges=3
           do i=1,n_ranges
             n_gammas_in_range(i)=4**n_letters ! eg 1 to 16 for 2-letter code
           enddo
        endif


        do 510 i_ranges=1,n_ranges

        do 504 i504=1,n_gammas_in_range(i_ranges)  
           read (igamma,*,iostat=read_status,err=99) optgamma(i504)
99         if (read_status.ne.0) then
             write(SO,*)  'failure reading from gamma file (1st read)'
             write(SO,*) 'error number ',read_status
             write(SO,*) optgamma(i504)
             stop
           endif
           write (oarchv,444) optgamma(i504)
444        format(f8.4)

504    continue

        do 500 i500=1,20
           icl=class(i500)
           icl_con=class_con(i500)

           do 501 i501=1,20
             jcl=class(i501)
             jcl_con=class_con(i501)
              gammapt_contact=sort_contact(icl_con,jcl_con)

              do 502 i502=1,20
                 ipcl=class(i502)
                 do 503 i503=1,20
                   jpcl=class(i503)

              if ( i_ranges .le. 2) then
                 if (four_plus) then
                   do ic1 = 1,2
                   gammapt=sort_plus(icl,jcl,ipcl,jpcl,ic1)
                   qchrg2(i500,i501,i502,i503,ic1,i_ranges)=optgamma(gammapt)
                   enddo
                 else ! four_plus false
                      if (i_ignore_chain_dirn) then
                       gammapt=sort_sym(icl,jcl,ipcl,jpcl)
                       gamma_temp=optgamma(gammapt)
                      else
                       gammapt=sort(icl,jcl,ipcl,jpcl)
                       gamma_temp=optgamma(gammapt)
                      endif
                 endif ! four_plus
               endif ! i_ranges 2
                                                                                                     
                   if (i_3rd_is_contact.and.(n_letters.eq.4))then
                    if(i_ranges.gt.2) then
                       gamma_temp=optgamma(gammapt_contact)
                    endif
                   endif
c                  third distance class Go contacts
              if ((.not.i_3rd_is_contact).and.(n_letters.eq.4))then
                 if(i_ranges.gt.2) then
                     gamma_temp=optgamma(gammapt)
                 endif

              endif
                                                                                  
c   If alternative potential qchrg for the corresponding well
c   should be set to 0
                  if((altpotflag(1).gt.0).and.(i_ranges.eq.3))then
                      qchrg(i500,i501,i502,i503,i_ranges)=0.0D0
                  elseif((altpotflag(2).gt.0).and.(i_ranges.eq.4))then
                      qchrg(i500,i501,i502,i503,i_ranges)=0.0D0
c   these are added by Garegin Papoian, 4/8/04
c   Modifying medium-range AMH interactions
                  elseif((altpotflag(2).gt.0).and.(i_ranges.eq.2))then
                      qchrg(i500,i501,i502,i503,i_ranges)=gamma_temp*srplr(2)
c   Modifying short-range AMH interactions
                  elseif((altpotflag(2).gt.0).and.(i_ranges.eq.1))then
                      qchrg(i500,i501,i502,i503,i_ranges)=gamma_temp*srplr(1)
c   end of GAP
                  else
                   qchrg(i500,i501,i502,i503,i_ranges)=gamma_temp


                  endif
                                                                                       
503              continue
502            continue

501      continue
500     continue

510     continue

c----------tshen add the 21th one, just rely on the relations
c need to gamma(21,"1 - 20","1 - 20","1 - 20","1 - max_well+2")  and its oppo
c and gamma2 but what in gamma2??
c first try mychch=1.2 and mychpo=1.1
      mychch=1.2D0*1.2D0
      mychpo=1.2D0
                
      do myj=1,20
         do myip=1,20
            do myjp=1,20
               do mywel=1, max_well+2
                  if((class(myj).eq.4).or.(class(myip).eq.4).or.(class(myjp).eq.4)) then
                     qchrg(21,myj,myip,myjp,mywel)=qchrg(7,myj,myip,myjp,mywel)
                     qchrg(myj,21,myip,myjp,mywel)=qchrg(myj,7,myip,myjp,mywel)
                  elseif ((class(myj).eq.1).or.(class(myip).eq.1).or.(class(myjp).eq.1)) then
                     qchrg(21,myj,myip,myjp,mywel)=mychpo*qchrg(7,myj,myip,myjp,mywel)
                     qchrg(myj,21,myip,myjp,mywel)=mychpo*qchrg(myj,7,myip,myjp,mywel)
                  else
                     qchrg(21,myj,myip,myjp,mywel)=mychch*qchrg(7,myj,myip,myjp,mywel)
                     qchrg(myj,21,myip,myjp,mywel)=mychch*qchrg(myj,7,myip,myjp,mywel)
                  endif
               enddo
            enddo
         enddo
      enddo


      do myj=1,20
         do myip=1,20
            do myjp=1,20
               do myk=1, 2
                  do mykp=1,2
                     if((class(myj).eq.4).or.(class(myip).eq.4).or.(class(myjp).eq.4)) then
                        qchrg2(21,myj,myip,myjp,myk,mykp)=qchrg2(7,myj,myip,myjp,myk,mykp)
                        qchrg2(myj,21,myip,myjp,myk,mykp)=qchrg2(myj,7,myip,myjp,myk,mykp)
                     elseif ((class(myj).eq.1).or.(class(myip).eq.1).or.(class(myjp).eq.1)) then
                        qchrg2(21,myj,myip,myjp,myk,mykp)=mychpo*qchrg2(7,myj,myip,myjp,myk,mykp)
                        qchrg2(myj,21,myip,myjp,myk,mykp)=mychpo*qchrg2(myj,7,myip,myjp,myk,mykp)
                     else
                        qchrg2(21,myj,myip,myjp,myk,mykp)=mychch*qchrg2(7,myj,myip,myjp,myk,mykp)
                        qchrg2(myj,21,myip,myjp,myk,mykp)=mychch*qchrg2(myj,7,myip,myjp,myk,mykp)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      do myip=1,20
            do myjp=1,20
               do mywel=1, max_well+2
                  if((class(myip).eq.4).or.(class(myjp).eq.4)) then
                     qchrg(21,21,myip,myjp,mywel)=qchrg(7,21,myip,myjp,mywel)
                     qchrg(21,21,myip,myjp,mywel)=qchrg(21,7,myip,myjp,mywel)
                  elseif ((class(myip).eq.1).or.(class(myjp).eq.1)) then
                     qchrg(21,21,myip,myjp,mywel)=mychpo*qchrg(7,21,myip,myjp,mywel)
                     qchrg(21,21,myip,myjp,mywel)=mychpo*qchrg(21,7,myip,myjp,mywel)
                  else
                     qchrg(21,21,myip,myjp,mywel)=mychch*qchrg(7,21,myip,myjp,mywel)
                     qchrg(21,21,myip,myjp,mywel)=mychch*qchrg(21,7,myip,myjp,mywel)
                  endif
               enddo

               do myk=1, 2
                  do mykp=1,2
                     if((class(myip).eq.4).or.(class(myjp).eq.4)) then
                        qchrg2(21,21,myip,myjp,myk,mykp)=qchrg2(7,21,myip,myjp,myk,mykp)
                        qchrg2(21,21,myip,myjp,myk,mykp)=qchrg2(21,7,myip,myjp,myk,mykp)
                     elseif ((class(myip).eq.1).or.(class(myjp).eq.1)) then
                        qchrg2(21,21,myip,myjp,myk,mykp)=mychpo*qchrg2(7,21,myip,myjp,myk,mykp)
                        qchrg2(21,21,myip,myjp,myk,mykp)=mychpo*qchrg2(21,7,myip,myjp,myk,mykp)
                     else
                        qchrg2(21,21,myip,myjp,myk,mykp)=mychch*qchrg2(7,21,myip,myjp,myk,mykp)
                        qchrg2(21,21,myip,myjp,myk,mykp)=mychch*qchrg2(21,7,myip,myjp,myk,mykp)
                     endif
                  enddo
               enddo
            enddo
         enddo

c-------end of tshen modi-----------

        if (i_non_add_contact) then
          do i=1,ngamma_non_add
            read (igamma,*,iostat=read_status,err=199) gamma_non_add(i)
199            if (read_status.ne.0) then
              write(SO,*) 'iailure reading from gamma file (2nd read)'
              write(SO,*)'error number ',read_status
              stop
            endif
          enddo
        endif

        close (igamma)

        open(ihbond,file='params/anti_HB',status='old',iostat=open_status)
        if (open_status.ne.0) then
          write(SO,*)'failure to open anti_HB file'
          write(SO,*) 'error number ',open_status
          stop
        endif

        do isit1 = 1,20
        read(ihbond,77)(anti_HB(isit1,isit2,1),isit2=1,20)
        enddo
        read(ihbond,*)
        do isit1 = 1,20
        read(ihbond,77)(anti_HB(isit1,isit2,2),isit2=1,20)
        enddo
        close(ihbond)

        open(ihbond,file='params/anti_NHB',status='old',iostat=open_status)
        if (open_status.ne.0) then
          write(SO,*)'failure to open anti_NHB file'
          write(SO,*) 'error number ',open_status
          stop
        endif

        do isit1 = 1,20
        read(ihbond,77)(anti_NHB(isit1,isit2,1),isit2=1,20)
        enddo
        read(ihbond,*)
        do isit1 = 1,20
        read(ihbond,77)(anti_NHB(isit1,isit2,2),isit2=1,20)
        enddo
        close(ihbond)

        open(ihbond,file='params/para_one',status='old')
        do isit1 = 1,20
        read(ihbond,*)para_one(isit1)
        enddo
        close(ihbond)

        open(ihbond,file='params/anti_one',status='old')
        do isit1 = 1,20
        read(ihbond,*)anti_one(isit1)
        enddo
        close(ihbond)

        open(ihbond,file='params/para_HB',
     *  status='old',iostat=open_status)
        if (open_status.ne.0) then
          write(SO,*) 'failure to open para_HB file'
          write(SO,*) 'error number ',open_status
          stop
        endif

        do isit1 = 1,20
        read(ihbond,77)(para_HB(isit1,isit2,1),isit2=1,20)
        enddo
        read(ihbond,*)
        do isit1 = 1,20
        read(ihbond,77)(para_HB(isit1,isit2,2),isit2=1,20)
        enddo
        close(ihbond)

77        format(20(f8.5,1x))

        return
      end
