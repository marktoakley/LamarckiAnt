
c     --------------------- gomb ----------------------

      subroutine gomb(  eta,index,
     *                  tempav,
     *                  ngomb,
     *                  A_to_nmo,E,ibiasgauss,bias_av,
     *                  bias_var,bias_prefactor,bias_F,
     *                  ibiaspoly,nbiaspoly,biaspoly)

c     ---------------------------------------------------

c      gomb (`Go many-body') evaluates the A-factors 
c      required for the gomb implementation
c
c     ---------------------------------------------------

      use amhglobals,  only: maxcnt,maxsiz,nmres,numlng,ilong,
     *  forse,vpotnt,rincsq,amh_gse,eres,maxr

      implicit none

       double precision, intent(in):: eta(1:maxcnt,1:4),ngomb,bias_av,bias_var,
     *  bias_prefactor,biaspoly(1:100)
       double precision, intent(out):: A_to_nmo(1:maxsiz,1:2),E(:,:),bias_F
      integer, intent(in):: index(1:maxcnt,1:4),nbiaspoly
      logical, intent(in):: tempav,ibiasgauss,ibiaspoly

c     internal variables:
         double precision A_to_n(1:maxsiz,1:2),m,c,E_addition,eta_prime
         double precision Agomb1(1:maxsiz,1:4),Agomb2(1:maxsiz,1:4),q
         integer i_ixn,isit1,isit2,i,j,i_tab,i_res,indx

c     --------------------- begin -----------------------

!     zero energy and force factor

      E=0.0D0
      bias_F=0.0D0

c     initialise --- i labels residue, j table

      do 504 i=1,nmres
        do 501 j=1,4
          Agomb1(i,j)=0.0D0
          Agomb2(i,j)=0.0D0
501     enddo
504   enddo

         do 506 i_tab=1,4
           do 505 i_ixn=1,numlng(nmres,i_tab)
c  want to work out Agomb(i)=sum_j(theta_ij) where theta_ij
c  are the pair energies. This is slightly complicated by
c  the fact that i labels an atom (Ca or Cb) and things
c  are annoyingly stored in 4 tables. Agomb is the quantity
c  that we take mod of and ^(ngomb-1) in the i509 loop. The
c  second index in A_to_nmo labels alpha/beta (1/2)

             isit1=ilong(i_ixn,1,i_tab)
             isit2=ilong(i_ixn,2,i_tab)

             indx=index(i_ixn,i_tab)
             if (indx.gt.maxr) cycle !outside range of force table

c  for accuracy choose interpolate potential so it exactly 
c  (hopefully) corresponds to integrated force. This
c  naturally depends on method of interpolating the force
         
             m=     (forse(indx+1,i_ixn,i_tab)-
     *                       forse(indx,i_ixn,i_tab))

             c=forse(indx,i_ixn,i_tab)
     *                                *float(indx+1)
     *           -forse(indx+1,i_ixn,i_tab)
     *                                   *float(indx)

             eta_prime=1.0D0-eta(i_ixn,i_tab) 

             E_addition=vpotnt(indx+1,i_ixn,i_tab)+
     *      rincsq(i_ixn,i_tab)*
     *      eta_prime* (m*
     *      (float(indx+1)*(float(indx+1)-eta_prime)+
     *       eta_prime*eta_prime/3.0D0)  +
     *      c* (float(indx+1)-0.5D0*eta_prime)  )

             Agomb1(isit1,i_tab)=Agomb1(isit1,i_tab) + 
     *         E_addition

             Agomb2(isit2,i_tab)=Agomb2(isit2,i_tab) + 
     *         E_addition          

  505     continue
506     continue

       do 509 i_res=1,nmres
        A_to_nmo(i_res,1)  = abs( ( Agomb1(i_res,1) + Agomb1(i_res,2)
     +        + Agomb2(i_res,1) + Agomb2(i_res,3) ) )**(ngomb-1.0D0)
        A_to_nmo(i_res,2)  = abs( ( Agomb1(i_res,3) + Agomb1(i_res,4)
     +        + Agomb2(i_res,2) + Agomb2(i_res,4) ) )**(ngomb-1.0D0)
509    continue


      if (tempav.or.ibiasgauss.or.ibiaspoly) then
        do 510 i_res=1,nmres
          A_to_n(i_res,1)  =  A_to_nmo(i_res,1) * ( Agomb1(i_res,1)
     +      + Agomb1(i_res,2) + Agomb2(i_res,1) + Agomb2(i_res,3) )
          A_to_n(i_res,2)  =  A_to_nmo(i_res,2) * ( Agomb1(i_res,3)
     +      + Agomb1(i_res,4) + Agomb2(i_res,2) + Agomb2(i_res,4) )
        E(1,5)=E(1,5) + A_to_n(i_res,1) + A_to_n(i_res,2)
!       write(SO,*) i_res,A_to_n(i_res,1) + A_to_n(i_res,2)
510     continue
        E(1,5)=E(1,5)*0.5D0
        amh_gse=-4.0D0*eres*float(nmres)

c    These below are backwards ways of adding an
c    unmbrelling potential

        if (ibiasgauss) then
          E(1,10)=bias_prefactor*exp(-(0.5D0/bias_var)*
     *                    (E(1,5)/amh_gse-bias_av)**2)
          bias_F=-(E(1,10)/bias_var)*
     *                ((E(1,5)/amh_gse)-bias_av)/amh_gse
        elseif (ibiaspoly) then
          q=E(1,5)/amh_gse
          E(1,10)=0.0D0
          do i=1,nbiaspoly
            E(1,10)=E(1,10)+biaspoly(i)*q**i
            bias_F=bias_F+(1.0D0/amh_gse)*(float(i))
     *                             *biaspoly(i)*q**(i-1) 
          enddo
        endif
      endif

c     ---------------------- done -----------------------

      return
      end
