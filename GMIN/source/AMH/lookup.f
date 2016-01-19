
c     --------------------- lookup ----------------------

      subroutine lookup(maxr,crdixn,vpotnt,forse,pro_cord,f_cord,
     *                  trgeng,numlng,nmres,rincinv,rincsq,
     *                  igomb,ngomb,tempav,maxtab,ilong,
     *                  E,ibiasgauss,bias_av,bias_var,
     *                  bias_prefactor,ibiaspoly,nbiaspoly,
     *                  biaspoly,i_Qbias_a,i_Qbias_b,ccev_dist,pexcld)

c     ---------------------------------------------------

c     LOOKP  performs lookup of potntlial or force
c            for interacting sites given the distances
c            between them

c     arguments:

c        tempav- if true, then find potntlial energy
c
c        i_ixn   Counter which ranges over all of the interactions.  The 
c                actual number depends on if nearest neighbor interactions are
c                included or not, but is aprox.  N(N-1)/2
c
c     ---------------------------------------------------

      use amhglobals,  only:maxcnt,maxsiz,i_Rg_bias,iexcld,minmr,maxmr,
     *      iexcld_gamma,maxcrd,SO
      use amh_interfaces, only:gomb,additive_ev,Rg_bias,ev_gamma,
     *                      Q_bias_seg_a,Q_bias_seg_b
      
      implicit none
      
      integer, intent(in):: maxr,maxtab
      integer, intent(in):: crdixn(maxtab,2),
     *  numlng(0:maxsiz,maxtab),nmres,ilong(maxcnt,2,maxtab),
     *  nbiaspoly
       double precision, intent(in):: vpotnt(0:maxr+1,maxcnt,maxtab),
     *  forse(0:maxr+1,maxcnt,maxtab),pro_cord(maxsiz,3,maxcrd),
     *  rincinv(maxcnt,maxtab),
     *  rincsq(maxcnt,maxtab),ngomb,bias_av,bias_var,bias_prefactor,
     *  biaspoly(1:100),ccev_dist(maxcnt,maxtab),pexcld
       double precision, intent(out):: f_cord(maxsiz,3,maxcrd),
     *  E(:,:),trgeng(maxtab,3)
      logical, intent(in):: igomb,tempav,ibiasgauss,ibiaspoly,
     *           i_Qbias_a,i_Qbias_b


         integer index(maxcnt,maxtab),isit1,isit2,tab


         double precision distne(maxcnt,maxtab),tforce(maxcnt,maxtab),
     *        potntl(maxcnt),f_cord_temp(maxsiz,3,maxcrd),
     *        eta(maxcnt,maxtab),A_to_nmo(1:maxsiz,1:2),
     *        xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab),
     *        zdiff(maxcnt,maxtab),E_temp(size(E,1),size(E,2)),bias_F


c     internal variables:

         integer i_ixn,ia,ib,i515,indx 

          double precision m,c,eta_prime

c -comment back in if want to use these variables
c        double precision forsego1,forsego2,arsetot,deltaV,
c    *         oldamhE,olddist(maxcnt,maxtab),arsetot2,
c    *         deltadist(maxcnt,maxtab),deltaVbias,
c    *         oldbiasE,
c    *         oldtforce(maxcnt,maxtab),oldbias_F,
c        double precision V240,V240old,deltaV240
c        double precision f_cord_old(maxsiz,3,maxcrd)
c        integer i1,i2

c     --------------------- begin -----------------------

!     zero force and energy

      f_cord=0.0D0
      E=0.0D0


c        find (x_i-x_j), (y_i-y_j), (z_i-z_j) for
c        each pair of interacting sites i and j

         do 1001 tab=1,maxtab
         do 501 i_ixn=1,numlng(nmres,tab)

            isit1=ilong(i_ixn,1,tab)
            isit2=ilong(i_ixn,2,tab)

            xdiff(i_ixn,tab)=pro_cord(isit1,1,crdixn(tab,1)) -
     *                  pro_cord(isit2,1,crdixn(tab,2))


            ydiff(i_ixn,tab)=pro_cord(isit1,2,crdixn(tab,1)) -
     *                  pro_cord(isit2,2,crdixn(tab,2))

            zdiff(i_ixn,tab)=pro_cord(isit1,3,crdixn(tab,1)) -
     *                  pro_cord(isit2,3,crdixn(tab,2))

            distne(i_ixn,tab)=dsqrt(xdiff(i_ixn,tab)**2 +
     *                        ydiff(i_ixn,tab)**2 +
     *                        zdiff(i_ixn,tab)**2)

c     find potntlial index; distances should be in array distne

         index(i_ixn,tab)=int( distne(i_ixn,tab)*rincinv(i_ixn,tab) )
         eta(i_ixn,tab)=distne(i_ixn,tab)*rincinv(i_ixn,tab)
     *          -dble(index(i_ixn,tab))
         tforce(i_ixn,tab)=0.0D0



501    continue
1001   continue

!         write(6,*)'i_Qbias'
c        if (i_Qbias) then
c           call Q_bias(distne,f_cord_temp,nmres,E_temp,
c    *                 numlng,ilong,
c    *                 xdiff,ydiff,zdiff)
c          f_cord=f_cord+f_cord_temp
c          E=E+E_temp
c       endif

        if (i_Qbias_a) then
           call Q_bias_seg_a(distne,f_cord_temp,nmres,E_temp,xdiff,ydiff,zdiff)
           f_cord=f_cord+f_cord_temp
           E=E+E_temp
        endif

        if (i_Qbias_b) then
           call Q_bias_seg_b(distne,f_cord_temp,nmres,E_temp,xdiff,ydiff,zdiff)
           f_cord=f_cord+f_cord_temp
           E=E+E_temp
        endif

        if (i_Rg_bias) then
           call Rg_bias(pro_cord,f_cord_temp,E_temp,tempav)
           f_cord=f_cord+f_cord_temp
           E=E+E_temp
        endif

        if (iexcld) then
            call additive_ev(distne,f_cord_temp,nmres,E_temp,
     *                 numlng,ilong,tempav,crdixn,
     *                 xdiff,ydiff,zdiff,ccev_dist,pexcld)
           f_cord=f_cord+f_cord_temp
           E=E+E_temp
        endif

        if (iexcld_gamma) then
           call ev_gamma(pro_cord,f_cord_temp,E_temp,tempav)
           f_cord=f_cord+f_cord_temp
           E=E+E_temp
        endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find potntlial for each of the interacting sites
c     note that interpolation used is to ensure that 
c     energy = integrated force (that's why don't 
c     linearly interpolate for energy like we do for
c     F(r)/r
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if( tempav.and.(.not.igomb) )then
             do 606 tab=1,maxtab

             trgeng(tab,1)=0.0D0
             trgeng(tab,2)=0.0D0
             trgeng(tab,3)=0.0D0

             do 505 i_ixn=1,numlng(nmres,tab)
               isit1=ilong(i_ixn,1,tab)
               isit2=ilong(i_ixn,2,tab)

               indx=index(i_ixn,tab)
               if (indx.gt.maxr) cycle !outside range of force table
         
               m=(forse(indx+1,i_ixn,tab)-
     *                       forse(indx,i_ixn,tab))
c note that m is therefore what i had in my notes times rinc
c `c', on the other hand is the same:  c=forse(indx)-m*float(indx)
c which simplifies to ...

         c=forse(indx,i_ixn,tab)*dble(indx+1)
     *           -forse(indx+1,i_ixn,tab)*dble(indx)

           eta_prime=1.0D0-eta(i_ixn,tab) 

c the following formula is in my notes (remember m here
c not quite identical to there (see above)

               potntl(i_ixn)=vpotnt(indx+1,i_ixn,tab)+
     *      rincsq(i_ixn,tab)*
     *      eta_prime* (m* 
     *      (dble(indx+1)*(dble(indx+1)-eta_prime)+
     *       eta_prime*eta_prime/3.0D0)  +
     *      c* (dble(indx+1)-0.5D0*eta_prime)  )
               

c    calc contribution to energy due to 1st mem
c    and that due to mems without excluded vol
c    this is rough in that proper interpolation
c    as above is not used. Thus will apparently  get small
c    residual excluded volume contribution.


c              if (i_ixn.eq.240.and.tab.eq.1)  then
c                      V240=potntl(i_ixn)
c              endif


               E(1,5) = E(1,5) + potntl(i_ixn)

c these lines were used to check out contributions from a
c particular gamma value (we zeroed out all other gamma
c values in gamma.dat
c              if (abs(potntl(i_ixn)).gt.0.001) then
c                 write(SO,*) '****************************'
c                 write(SO,*) 'lookup reports:'
c                 write(SO,*) 'tab=',tab
c                 write(SO,*) 'isit1,isit2,r,contrib==',isit1,isit2,
c    *                distne(i_ixn,tab),potntl(i_ixn)
c                 write(SO,*) 'r increment=',1.0/rincinv(i_ixn,tab)
c                 write(SO,*) 'grid points =',
c    *                  real(indx)*(1.0/rincinv(i_ixn,tab)),
c    *                  real(indx+1)*(1.0/rincinv(i_ixn,tab))
c              endif 

               if ((isit2-isit1) .lt. minmr) then
                    E(1,7) = E(1,7) + potntl(i_ixn)
                    trgeng(tab,1)=trgeng(tab,1)+potntl(i_ixn)
               elseif ((isit2-isit1) .gt. maxmr) then
                    E(1,8) = E(1,8) + potntl(i_ixn)
                    trgeng(tab,3)=trgeng(tab,3)+potntl(i_ixn)
               else
                    E(1,13) = E(1,13) + potntl(i_ixn)
                    trgeng(tab,2)=trgeng(tab,2)+potntl(i_ixn)
               endif
  505        continue
  606        continue
              
        endif                        ! end of     if ( tempav and (not igomb)  )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     find force for each of the interacting sites
          
         if (.not.igomb) then
           do 1595 tab=1,4
           do 595 i_ixn=1,numlng(nmres,tab)
            isit1=ilong(i_ixn,1,tab)
            isit2=ilong(i_ixn,2,tab)
            indx=index(i_ixn,tab)
            if (indx.gt.maxr) cycle !outside range of force table
            tforce(i_ixn,tab)=
     *            (1.0D0-eta(i_ixn,tab))*forse(indx,i_ixn,tab)
     *             +eta(i_ixn,tab)*forse(indx+1,i_ixn,tab)

  595     continue
 1595     continue

         else
           call gomb(   eta,index,tempav,ngomb,A_to_nmo,E_temp,ibiasgauss,bias_av,
     *                  bias_var,bias_prefactor,bias_F,ibiaspoly,nbiaspoly,biaspoly)
           E=E+E_temp

          do 777 tab=1,maxtab
          if (tab.eq.1.or.tab.eq.2) then
            ia=1
          else
            ia=2
          endif
          if (tab.eq.1.or.tab.eq.3) then
            ib=1
          else
            ib=2
          endif
          do 506 i_ixn=1,numlng(nmres,tab)
            isit1=ilong(i_ixn,1,tab)
            isit2=ilong(i_ixn,2,tab)
            indx=index(i_ixn,tab)
            if (indx.gt.maxr) cycle !outside range of force table


            tforce(i_ixn,tab)=
     *       ((1.0D0-eta(i_ixn,tab))*forse(indx,i_ixn,tab)
     *        +eta(i_ixn,tab)*forse(indx+1,i_ixn,tab)  )
     *         *0.5D0*ngomb*(A_to_nmo(isit1,ia)+A_to_nmo(isit2,ib))
 
            tforce(i_ixn,tab)=tforce(i_ixn,tab)*(1.0D0+bias_F)

  506     continue
  777     continue

         endif



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        the below lines may be (selectively) uncommented to 
c        check that energy corresponds
c        to the integrated force. Also uncomment the
c        'V240' line in E-cal above if required.
c      
c        
c        open(77,file='temp.dat',status='old')
c        read(77,*) oldamhE
c        if (oldamhE.eq.99.0) goto 123
c        read(77,*) oldbiasE
c        read(77,*) V240old
c        read(77,*) oldbias_F
c        deltaV=E(1,5)-oldamhE
c        deltaVbias=E(1,10)-oldbiasE
c        deltaV240=V240-V240old
c        arsetot=0d0
c        arsetot2=0d0
c        do tab=1,maxtab
c        do  i_ixn=1,numlng(nmres,tab)
c          read(77,*) olddist(i_ixn,tab)
c          read(77,*) oldtforce(i_ixn,tab)
c          deltadist(i_ixn,tab)=distne(i_ixn,tab)-olddist(i_ixn,tab)
c          arsetot=arsetot+0.5*
c    +         (   tforce(i_ixn,tab)*distne(i_ixn,tab)/(1.0+bias_F)
c    +      +oldtforce(i_ixn,tab)*olddist(i_ixn,tab)/(1.0+oldbias_F) ) 
c    +                   *deltadist(i_ixn,tab)
c          arsetot2=arsetot2+0.5*
c    +         (   tforce(i_ixn,tab)*distne(i_ixn,tab)
c    +                                     *bias_F/(1.0+bias_F)
c    +      +oldtforce(i_ixn,tab)*olddist(i_ixn,tab)
c    +                             *oldbias_F/(1.0+oldbias_F)  )
c    +                   *deltadist(i_ixn,tab)
c        enddo
c        enddo
c        write(SO,*) 'amhE (inc excld V)',E(1,5)
c        write(SO,*) 'sum of -Fxdelta_r,and delta_V',-arsetot,deltaV
c        write(SO,*)  'ratio',-arsetot/deltaV
c        write(SO,*) 'bias sum of -Fxdelta_r,and delta_V',
c    +                                  -arsetot2,deltaVbias
c        write(SO,*) 'bias_F,oldbias_F',bias_F,oldbias_F
c123      close(77)

c        open(77,file='temp.dat',status='unknown')
c        write(77,*) E(1,5)
c        write(77,*) E(1,10)
c        write(77,*) V240
c        write(77,*) bias_F
c        do tab=1,maxtab
c        do  i_ixn=1,numlng(nmres,tab)
c          write(77,*) distne(i_ixn,tab)
c          write(77,*) tforce(i_ixn,tab)
c        enddo
c        enddo
c        close(77)

c         write(SO,*) '240 F, Fold, Fav',
c    +                tforce(240,1)*distne(240,1)/(1.0+bias_F),
c    +             oldtforce(240,1)*olddist(240,1)/(1.0+oldbias_F),
c    +             0.5*(tforce(240,1)*distne(240,1)/(1.0+bias_F)+
c    +             oldtforce(240,1)*olddist(240,1)/(1.0+oldbias_F))          
c         write(SO,*) ' 240 r, rold, delr',distne(240,1),olddist(240,1),
c    +                       deltadist(240,1)
c         write(SO,*) '240 V, Vold, delV',V240,V240old,V240-V240old
c         write(SO,*) '240  -F*delta_dist ',-deltadist(240,1)*
c    +             0.5*(tforce(240,1)*distne(240,1)/(1.0+bias_F)+
c    +             oldtforce(240,1)*olddist(240,1)/(1.0+oldbias_F))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*         if (tempav) write(SO,*) tab,E(1,5)

c        find force for each interaction

         do 616 tab=1,maxtab
         do 515 i515=1,numlng(nmres,tab)
            xdiff(i515,tab)=tforce(i515,tab)*xdiff(i515,tab)
            ydiff(i515,tab)=tforce(i515,tab)*ydiff(i515,tab)
            zdiff(i515,tab)=tforce(i515,tab)*zdiff(i515,tab)
  515    continue
  616    continue

         do 613 tab=1,maxtab

         do i_ixn=1,numlng(nmres,tab)
            isit1=ilong(i_ixn,1,tab)
            isit2=ilong(i_ixn,2,tab)
            
               f_cord(isit1,1,crdixn(tab,1))=
     *        f_cord(isit1,1,crdixn(tab,1)) + xdiff(i_ixn,tab)

               f_cord(isit1,2,crdixn(tab,1))=
     *        f_cord(isit1,2,crdixn(tab,1)) + ydiff(i_ixn,tab)

               f_cord(isit1,3,crdixn(tab,1))=
     *        f_cord(isit1,3,crdixn(tab,1)) + zdiff(i_ixn,tab)

               f_cord(isit2,1,crdixn(tab,2))=
     *        f_cord(isit2,1,crdixn(tab,2)) - xdiff(i_ixn,tab)

               f_cord(isit2,2,crdixn(tab,2))=
     *        f_cord(isit2,2,crdixn(tab,2)) - ydiff(i_ixn,tab)


               f_cord(isit2,3,crdixn(tab,2))=
     *        f_cord(isit2,3,crdixn(tab,2)) - zdiff(i_ixn,tab)

           enddo

  613      continue 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------- done -----------------------

      return
      end
