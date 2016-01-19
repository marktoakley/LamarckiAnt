
c     --------------------- gengrd ----------------------

      subroutine gengrd(maxcnt,ilong,numlng,nmres,maxs,maxrij,xwork,delta,delte,
     *                  maxsiz,deltz,rinc,idigns,oarchv,rsep,rcutAMH)

c     ---------------------------------------------------

c     GENGRD generates the r-grid

c     arguments:

c        maxcnt- maximum number of interactions (i)
c        ilong - list of interacting sites (i)
c        numlng- breakdown of structure of ilong (i)
c        nmres - number of target sites (i)
c        maxs  - number of r-grid points (i)
c        maxrij- maximum grid value for each  (w)
c        xwork - work array (w)
c        delta - Gaussian well-width (i)
c        delte - exponent for scaling of Gaussian as 
c                function of i-j (i)
c        maxsiz- maximum protein length (i)
c        deltz - denominator in Gaussian for each interaction
c                pair (o)
c        rinc  - r-grid increment for each interaction pair (o)


c       indxt  - total number of interactions.
c       idiff  - distance in number of residues
c
c     ---------------------------------------------------

      use amhglobals,  only: n_letters

c     set required parameters

      implicit none

c     argument declarations:

         logical idigns

         integer maxsiz,maxcnt,ilong(maxcnt,2),numlng(0:maxsiz),nmres,maxs,oarchv

         double precision deltz(maxcnt),rinc(maxcnt),xwork(maxcnt),
     *        maxrij(maxcnt),delta,delte,rsep,rcutAMH
     
c     internal variables:

         integer indxt,idiff

c        --- do loop indices ---

         integer i_indx
 
c        --- implied do loop indices ---


         double precision deltc,scalr


c     --------------------- begin -----------------------

c     --- diagnostics ---

c      write(oarchv,145)maxcnt,numlng(nmres),nmres,maxs,
c     *                 maxsiz,delta,delte
c  145 format(/'Gengrd: maxcnt ',i5,' numlng ',i5,
c     *        ' nmres ',i3,' maxs ',i4,' maxsiz ',i4,
c     *        ' delta ',1pe10.3,' delte ',1pe10.3)

c     --- end diagnostics ---

c     initialize number of interactions

      indxt=numlng(nmres)
      numlng(0)=0

c     initialize maximum r_{ij} 

      do 500 i_indx=1,indxt
         idiff=ilong(i_indx,2) - ilong(i_indx,1)    ! distance in sequence space
         maxrij(i_indx)=min( rsep*float(idiff), rcutAMH )
  500 continue
c     set endpoints so that gaussian potential is continuous,
c     i.e., V and F at endpoints are 'small'

      deltc=0.5D0/delta**2
      scalr=(-delte*2)

      do 501 i_indx=1,indxt
         deltz(i_indx)=deltc*(float(ilong(i_indx,2) - ilong(i_indx,1) ))**scalr
        if ( (ilong(i_indx,2) - ilong(i_indx,1) .gt. 4 ) .and.
     *  (ilong(i_indx,2) - ilong(i_indx,1) .lt. 13).and.n_letters.eq.4) then

         deltz(i_indx)=2.0D0*(float(ilong(i_indx,2) - ilong(i_indx,1) ))**(-0.60D0) 
        endif

  501 continue

c     set upper bound for r-grid

      deltc=-log(1.0e-08)

c     set addend for determining maximum r-grid value

      do 502 i_indx=1,indxt
         xwork(i_indx)=sqrt(deltc/deltz(i_indx))
  502 continue
c     set endpoint for each interaction pair

      do 503 i_indx=1,indxt
         maxrij(i_indx)=maxrij(i_indx) + xwork(i_indx)
  503 continue
c     initialize grid increments based on min and max
c     distance values for each constraint

      do 504 i_indx=1,indxt
         rinc(i_indx)=maxrij(i_indx)/float(maxs-1)
  504 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- diagnostics ---
      if( idigns )then
c        list spanning interval for each | i - j |
c         write(oarchv,113)indxt
c  113    format(/'b1 and deltz:indx ',i8)
c         do 517 i517=1,min( 53,nmres )
cc            if( (ilong(i517,1).eq.4).and.
cc     *          (ilong(i517,2).eq.88) )then 
cc            do 518 i518=numlng(i517-1)+1,
cc     *              min(numlng(i517-1)+10,numlng(i517))
c             do 519 i519=1,2
c               if( i519.eq.1 )then
c                  i518=numlng(i517-1)+1
c               else
c                  i518=numlng(i517)
c               endif
c               write(oarchv,114)ilong(i518,1),ilong(i518,2),
c     *                          maxrij(i518)-xwork(i518),
c     *                          xwork(i518),deltz(i518),
c     *                          rinc(i518)
c  114          format(2(i3,1x),4(1x,1pe10.3))
c  519       continue
cc  518       continue
cc            endif
c  517    continue
      endif
c     --- end diagnostics ---
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ---------------------- done -----------------------

      return
      end
