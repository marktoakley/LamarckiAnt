
c     --------------------- gentab ----------------------

      subroutine gentab

c     ---------------------------------------------------

c     GENTAB is driver routine for generating required tables

c     ---------------------------------------------------
c
c  Subroutine Logic
c
c     * Call time
c     * Loop over number of tables
c       * Determine the interacting subset of residues
c         [ this may be historical ]
c       * Call potind
c
c       * Call gengrd
c 
c       * If Table = 1, passi = true
c       * If Table = numtab, passf = true
c
c       * Call Gaspot
c
c     * Call time
c     * Construct hbond tables
c
c  Subroutine Variables
c
c      rsep

      use amhglobals,  only:SO,numtab,crdixn,nmres,
     *  idigns,maxsiz,maxcnt,numlng,ilong,oarchv,maxs,work6,work3,
     *  delta,delte,deltz,rinc,rcutAMH,target,bondln,numcrd,
     *  n_divs_max,targ_dist,maxcrd,ibiasfile,i_rep,maxtab,iss_struct,
     *  i_Qbias_a,width_Qexp_a,del_r_a,Q_ij_a,dQ_dr_ij_a,i_bias_native_a,
     *  i_Qbias_b,i_bias_native_b,totalssnorm,numconst_a, numconst_b,
     *  del_r_b,Q_ij_b,dQ_dr_ij_b,width_Qexp_b,ss_a,ss_b,ss_dist

      use amh_interfaces, only:rep_contact

      implicit none

c     internal variables:

         logical passi,passf

         integer jbegn,jend, i_tab, i_dist,i_diff,tab,i_ixn,isit1,isit2,nmres_check,idummy,
     *           i_res,i_res1, i_res2

         double precision rsep, max_cont,max_rij,
     *         small,xdiff,ydiff,zdiff,
     *         bias_cord(maxsiz,3,maxcrd),
     *         xdifftemp,ydifftemp,zdifftemp

        integer tempNres,seq_a(500),seq_b(500),icoord,sa_alpha(500),ss_alpha(500),
     *          sa_beta(500),ss_beta(500)

        double precision Q_var_a,del_Q_a,Q_var_b,del_Q_b,cords(500,3,3)
        
        character blah*38

c     required subroutines

        external potind,gengrd,gaspot,ev_set_up

c     --------------------- begin -----------------------
c        write(SO,*) 'gentab'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     generate potential from memory proteins
c     construct potential and force tables
c     find tables for alpha-alpha coordinates

      do 500 i_tab=1,numtab

c        determine which residues will interact given
c        the protein length and 'mutation' region

         if( (crdixn(i_tab,1).eq.crdixn(i_tab,2)).and.
     *       (crdixn(i_tab,1).eq.1) )then
            jbegn=2
            jend=nmres
            rsep=5.0D0
         elseif( crdixn(i_tab,1).ne.crdixn(i_tab,2) )then
            jbegn=1
            jend=nmres
             rsep=5.5D0
         else
            jbegn=1
            jend=nmres
            rsep=10.5D0
         endif

         idigns=.false.
         call potind(maxsiz,maxcnt,numlng(0,i_tab),ilong(1,1,i_tab),jbegn,jend,nmres,i_tab)

c        set analysis parameters to false;
c        full tables are to be generated

c        construct r-grid for potential

c         if( i_tab.eq.1 )idigns=.false.
         call gengrd(maxcnt,ilong(1,1,i_tab),
     *               numlng(0,i_tab),nmres,maxs,
     *               work6,work3,delta,delte,
     *               maxsiz,deltz(1,i_tab),
     *               rinc(1,i_tab),idigns,
     *               oarchv,rsep,rcutAMH)
         idigns=.false.

c        set passi true if this is the first call
c        to gaspot; otherwise set passi to false

         if( i_tab.eq.1 )then
            passi=.true.
         else
            passi=.false.
         endif

c        set passf true if this is the las tcall
c        to gaspot; otherwise set passf to false

         if( i_tab.eq.numtab )then
            passf=.true.
         else
            passf=.false.
         endif

c         if( i_tab.eq.2 )idigns=.true.
c         write(SO,*) 'calling gaspot'         
        call gaspot(maxsiz,target, bondln,i_tab,numcrd,crdixn(i_tab,1),crdixn(i_tab,2),passi,passf)

         idigns=.false.

500     continue


!*****************************************************
! set up excluded volumes

      call ev_set_up

! set up replica interaction if required

      if (i_rep) call rep_contact(target)
        
        open(iss_struct,file='params/alphz',status='old')
                                                                                                         
        read(iss_struct,991)blah
991     format(a38)
                                                                                                         
        if (blah.eq."This is a non-continuous protein chain")then
           write(6,*) 'I found a non-continutous protein'
           stop
        endif
        read(iss_struct,199)tempNres
199     format(i5)
        read (iss_struct,25)(seq_a(i_res),i_res=1,tempNres)
25      format(25(i2,1x))
        
        do 100 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,1),i_res=1,tempNres)
30          format(8(f8.3,1x))
100     continue
        do 101 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,2),i_res=1,tempNres)
101     continue
        do 102 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,3),i_res=1,tempNres)
102     continue
        read(iss_struct,25)(sa_alpha(i_res),i_res=1,tempNres)
        read(iss_struct,25)(ss_alpha(i_res),i_res=1,tempNres)
        close (iss_struct)

c    targ_dist used in Q_bias_seg alpha

         do 2501 i_res1=1, tempNres-1
          do  2502 i_res2=i_res1 +1, tempNres
            xdifftemp=cords(i_res1,1,1) - cords(i_res2,1,1)
            ydifftemp=cords(i_res1,2,1) - cords(i_res2,2,1)
            zdifftemp=cords(i_res1,3,1) - cords(i_res2,3,1)
            ss_dist(i_res1,i_res2,1)=dsqrt(xdifftemp**2 +
     *                        ydifftemp**2 + zdifftemp**2)
2502       enddo
2501      enddo

        open(iss_struct,file='params/betaz',status='old')
        read(iss_struct,991)blah
        if (blah.eq."This is a non-continuous protein chain")then
           write(6,*) 'I found a non-continutous protein'
           stop
        endif
        read(iss_struct,199)tempNres
        read (iss_struct,25)(seq_b(i_res),i_res=1,tempNres)
        do 200 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,1),i_res=1,tempNres)
200     continue
        do 201 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,2),i_res=1,tempNres)
201     continue
        do 202 icoord=1,3
            read (iss_struct,30)(cords(i_res,icoord,3),i_res=1,tempNres)
202     continue
        read(iss_struct,25)(sa_beta(i_res),i_res=1,tempNres)
        read(iss_struct,25)(ss_beta(i_res),i_res=1,tempNres)
        close (iss_struct)

c    targ_dist used in Q_bias_seg beta
                                                                                                         
         do 3501 i_res1=1, tempNres-1
          do  3502 i_res2=i_res1 +1, tempNres
                                                                                                         
            xdifftemp=cords(i_res1,1,1) - cords(i_res2,1,1)
            ydifftemp=cords(i_res1,2,1) - cords(i_res2,2,1)
            zdifftemp=cords(i_res1,3,1) - cords(i_res2,3,1)
            ss_dist(i_res1,i_res2,2)=dsqrt(xdifftemp**2 + ydifftemp**2 + zdifftemp**2)
3502       enddo
3501      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Make table for calculation of Q, and dQ/dr<ij>
c      both scaled by 1/ ( 0.5 (N-1)(N-2) )
c      I have notes on this 'A specific Q-dependent potential'

        max_cont=1E-8    !max cont of Q<ij> to Q before ignore it

       if (i_Qbias_a)then
       do  i_diff = 1,  nmres
         Q_var_a=(float(i_diff)**(2.0D0*width_Qexp_a))
         max_rij=dsqrt( -2.0D0*log(max_cont)*Q_var_a)
         del_r_a(i_diff)=max_rij/n_divs_max
           do i_dist=0,n_divs_max-1
       if (ss_a) then
       Q_ij_a(i_dist,i_diff)=
     *     dexp(-0.5D0*(float(i_dist)*del_r_a(i_diff))**2/
     *     Q_var_a)*(2.0D0/(float(totalssnorm)))
       endif

      if (.not. ss_a) then
         Q_ij_a(i_dist,i_diff)=
     *     dexp(-0.5D0*(float(i_dist)*del_r_a(i_diff))**2/
     *    Q_var_a)* (2.0D0/(float((numconst_a-1)*(numconst_a-2))))
      endif
       
       enddo ! do   i_dist=0,n_divs_max-1
!       write(6,*)'numconst_a  in gentab ',numconst_a
         
         Q_ij_a(n_divs_max,i_diff)=0.0D0
!     use expansion to calc gradient, for accuracy
         do i_dist=1,n_divs_max
           small=-(0.5D0/Q_var_a)*(del_r_a(i_diff)**2)
     *                     *(2.0D0*float(i_dist-1)+1.0D0)
         del_Q_a=Q_ij_a(i_dist-1,i_diff)*(
     *          +small+(1.0D0/2.0D0)*small**2+(1.0D0/6.0D0)*small**3+
     *          (1.0D0/24.0D0)*small**4+(1.0D0/120.0D0)*small**5)
         dQ_dr_ij_a(i_dist,i_diff)= (del_Q_a/del_r_a(i_diff))

c         if (i_diff.eq.20) then   ! uncomment to check Q_ij_a and dQ_ij_a/dr_ij
c             write(6,*) float(i_dist)*del_r_a(i_diff),
c    *         Q_ij_a(i_dist-1,i_diff),dQ_dr_ij(i_dist,i_diff),
c    *        (Q_ij_a(i_dist,i_diff)-Q_ij_a(i_dist-1,i_diff))/del_r_a(i_diff)
c          endif
           enddo  ! do i_dist=1,n_divs_max
       enddo   ! i
       endif  ! if (i_Qbias_a)then

         if (i_Qbias_b)then

c  normailization of q constraint of secondary structure units
c  over right length of residues for proper Q normalization

       do  i_diff = 1,  nmres
         Q_var_b=(float(i_diff)**(2.0D0*width_Qexp_b))
         max_rij=dsqrt( -2.0D0*log(max_cont)*Q_var_b)
         del_r_b(i_diff)=max_rij/n_divs_max
         do i_dist=0,n_divs_max-1
        
         if (.not. ss_b) then
       Q_ij_b(i_dist,i_diff)=
     * dexp(-0.5D0*(float(i_dist)*del_r_b(i_diff))**2/
     * Q_var_b ) *(2.0D0/(float((numconst_b-1)*(numconst_b-2))))
        endif

      if (ss_b) then
         Q_ij_b(i_dist,i_diff)=dexp(-0.5D0*(float(i_dist)*del_r_b(i_diff))**2/
     *Q_var_b)*(2.0D0/(float((totalssnorm-1)*(totalssnorm-2))))
        endif
                                                                                                         
           enddo ! do i_dist=0,n_divs_max-1
                                                                                                         
         Q_ij_b(n_divs_max,i_diff)=0.0D0
                                                                                                         
!     use expansion to calc gradient, for accuracy
                                                                                                         
         do i_dist=1,n_divs_max
                                                                                                         
           small=-(0.5D0/Q_var_b)*(del_r_b(i_diff)**2)*(2.0D0*float(i_dist-1)+1.0D0)
                                                                                                         
           del_Q_b=Q_ij_b(i_dist-1,i_diff)*(
     *            +small+(1.0D0/2.0D0)*small**2+(1.0D0/6.0D0)*small**3+
     *            (1.0D0/24.0D0)*small**4+(1.0D0/120.0D0)*small**5)


           dQ_dr_ij_b(i_dist,i_diff)= (del_Q_b/del_r_b(i_diff))
                                                                                                         
c         if (i_diff.eq.20) then    ! uncomment to check Q_ij and dQ_ij/dr_ij
c         write(6,*) float(i_dist)*del_r_b(i_diff),
c     *     Q_ij_b(i_dist-1,i_diff),dQ_dr_ij_b(i_dist,i_diff),
c     *   (Q_ij_b(i_dist,i_diff)-Q_ij_b(i_dist-1,i_diff))/del_r_b(i_diff)
c          endif
                                                                                                         
               enddo  ! do i_dist=1,n_divs_max
           enddo   ! i
           endif ! if (i_Qbias_b)then

c       also calculate separations in the target, which will
c       need later in subroutine Q_bias
c       (transferred there via globals)
c       --or alternatively read in structure to bias to

           if ((.not.i_bias_native_a .and.  i_Qbias_a ).or.
     *         (.not.i_bias_native_b .and.  i_Qbias_b )) then

c         if (.not.i_bias_native_a .and.  i_Qbias_a )then

           open(ibiasfile,file='movieseg_bias',status='unknown',
     *                                            action='read')
           read(ibiasfile,1040) nmres_check,idummy,idummy,idummy
1040       format(4(i8,1x))
           if (nmres_check.ne.nmres) then
             write(SO,*) 'movieseg_bias file wrong',nmres_check,nmres
             stop
           endif
           read(ibiasfile,*)
1050       format(3x,3(1x,f8.4),2(4x,3(1x,f8.4)))
           do i_res=1,nmres
             read(ibiasfile,1050)  bias_cord(i_res,1,1),
     *       bias_cord(i_res,2,1),bias_cord(i_res,3,1),   
     *       bias_cord(i_res,1,2),bias_cord(i_res,2,2),   
     *       bias_cord(i_res,3,2),bias_cord(i_res,1,3),   
     *       bias_cord(i_res,2,3),bias_cord(i_res,3,3)
           enddo
           close(ibiasfile)

          endif 

        if ((i_bias_native_a .and.  i_Qbias_a ).or.
     *      (i_bias_native_b .and.  i_Qbias_b )) then

             bias_cord=target
             write(6,*) 'target  used for constraint'

         do 1001 tab=1,maxtab
         do 501 i_ixn=1,numlng(nmres,tab)

            isit1=ilong(i_ixn,1,tab)
            isit2=ilong(i_ixn,2,tab)

            xdiff=bias_cord(isit1,1,crdixn(tab,1)) -
     *                  bias_cord(isit2,1,crdixn(tab,2))
            ydiff=bias_cord(isit1,2,crdixn(tab,1)) -
     *                  bias_cord(isit2,2,crdixn(tab,2))
            zdiff=bias_cord(isit1,3,crdixn(tab,1)) -
     *                  bias_cord(isit2,3,crdixn(tab,2))

            targ_dist(i_ixn,tab)=dsqrt(xdiff**2 + ydiff**2 + zdiff**2)

501      enddo
1001     enddo

       endif
c     ---------------------- done -----------------------
c        write(SO,*) 'leaving gentab'      return
      end
