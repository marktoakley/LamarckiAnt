
c     --------------------- force ----------------------

      subroutine force(numpro,prcord,zrcord,avep,tempav,maxr,vpotnt,forse,
     *          trgeng,pexcld,numlng,nmres,rincinv,rincsq,ilong,crdixn,ires,
     *                 eqdist,hbond,oxexcldv,scl_call)
 
c     ---------------------------------------------------

c     FORCE  finds the forces and potential for the configurations in prcord

c     arguments:

c        maxsiz- maximum protein length (i)
c        maxpro- maximum number of proteins (i)
c        numpro- actual ensemble number (i)
c        prcord- current set of configurations (i)
c        maxcnt- maximum number of pair ixns (i)
c        zrcord- total force on each residue (o)
c        avep  - array used to record 1st and 2nd 
c                moments of the potential energy
c        tempav- if tempav, then compute total forces
c                and potential for averaging
c        maxr  - maximum r-grid length (i)
c        maxres- maximum number of residues (i)
c        ires  - coded amino acid sequence (i)
         
         use amhglobals,  only: maxtab,maxsiz,maxcnt,maxpro,maxcrd,
     *         oxscl,chrlscl,ramascl,igomb,ngomb,ibiasgauss,bias_av,bias_var,
     *         bias_prefactor,ibiaspoly,nbiaspoly,biaspoly,
     *         i_Qbias_a,i_Qbias_b,ccev_dist,i_non_add_contact,i_con_P_AP,
     *         i_rep,i_hbond_eastwood,SO

         use amh_interfaces, only:oxy,lookup,
     *    contact_P_AP,non_add_contact,rep_bias,dssp_hdrgn_eastwood,hdrgn

         use globals_alt, only:altpotflag

         implicit none

c     argument declarations:

         logical tempav,scl_call,hbond,oxexcldv

         integer numpro,crdixn(maxtab,2),maxr,numlng(0:maxsiz,maxtab),nmres,
     *           ilong(maxcnt,2,maxtab),ires(maxsiz)
     
         double precision prcord(maxsiz,3,maxpro,maxcrd),
     *        zrcord(maxsiz,3,maxpro,maxcrd),f_cord(maxsiz,3,maxcrd),pro_cord(maxsiz,3,maxcrd),
     *        avep(maxpro,2,50),rincinv(maxcnt,maxtab),
     *        rincsq(maxcnt,maxtab),vpotnt(0:maxr+1,maxcnt,maxtab),
     *        forse(0:maxr+1,maxcnt,maxtab),trgeng(maxtab,3),pexcld,E(2,50),eqdist(5),
     *        frcord(maxsiz,3,maxpro,maxcrd)

c     internal variables:
c        --- do loop indices ---

         integer  i_pro

          external  altpot_master

c     --------------------- begin -----------------------

c             write(SO,*) 'begin force'
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize storage for potential energy
c     and initialize zrcord

      avep=0.0D0
      zrcord=0.0D0
      pro_cord=0.0D0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (i_rep) then
        call rep_bias(prcord,frcord,E(1,43))
        zrcord=zrcord+frcord
        if (tempav) avep(1,:,:)=avep(1,:,:)+E
      endif
               
      do 500 i_pro=1,numpro
        pro_cord=prcord(:,:,i_pro,:)

c       find hbond force and potential
c       if ( .not.scl_call ) then
         if ( hbond ) then
            call hdrgn(pro_cord,f_cord,tempav,E)
            zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
            if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
         endif

         if ( i_hbond_eastwood ) then
            call dssp_hdrgn_eastwood(pro_cord,f_cord,tempav,E)
            zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
            if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
c           write(6,*)'E hydrogen bond ', E
       endif
 
c       if( .not.scl_call ) then
           call oxy(maxtab,ires,1,nmres,pro_cord,f_cord,
     *              eqdist,oxscl,chrlscl,ramascl,
     *              oxexcldv,numlng,ilong,nmres,E)

            zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
            if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
c        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c        find potential and force due to Ca/Cb via lookup 

c          write(SO,*)'in lookup'

         call lookup(maxr,crdixn,vpotnt,forse,pro_cord,f_cord,
     *      trgeng,numlng,nmres,rincinv,rincsq,igomb,ngomb,tempav,maxtab,ilong,
     *      E,ibiasgauss,bias_av,bias_var,bias_prefactor,ibiaspoly,nbiaspoly,
     *               biaspoly,i_Qbias_a,i_Qbias_b,ccev_dist,pexcld)
         zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
         if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       find non-additive contact contribution
c                  write(SO,*) 'in non_add'

          if (i_non_add_contact) then
             E=0.0D0
             call non_add_contact(pro_cord,nmres,tempav,f_cord,E(1,17))
             zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
             if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
          endif

         if((altpotflag(1).gt.0).or.(altpotflag(2).gt.0))then ! Anisotropically varying potential:
            E=0.0D0
            call altpot_master(pro_cord,f_cord,E,tempav,nmres)
             zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
             if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
         endif

         if (i_con_P_AP) then 
             E=0.0D0
             call contact_P_AP(pro_cord,nmres,f_cord,E(1,40:42))
             zrcord(:,:,i_pro,:)=zrcord(:,:,i_pro,:)+f_cord
             if (tempav) avep(i_pro,:,:)=avep(i_pro,:,:)+E
         endif
  500 continue                        ! End of loop over proteins

c   energy 
c   avep(num_pro , ???? , energy term , nmdif)

c force 
c zrcord(numres,i_axis,numpro,cord type ),

c     ---------------------- done -----------------------
      return
      end
