      subroutine additive_ev(distne,f_cord,nmres,E,
     *                  numlng,ilong,tempav,crdixn,
     *                  xdiff,ydiff,zdiff,ccev_dist,pexcld)


c     this routine works out the excluded volume interaction between
c     carbon atoms in the case when it is not already included in the AMH
c     potential, ie when we have non-additive amh part, but want additive
c     excluded volume contribution

      use amhglobals,  only: maxsiz,maxtab,maxcnt,maxcrd,ires


c     argument declarations

      implicit none

       double precision, intent(in):: distne(maxcnt,maxtab),pexcld,
     * xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab),
     * zdiff(maxcnt,maxtab),ccev_dist(maxcnt,maxtab)

       double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
      integer, intent(in):: nmres,numlng(0:maxsiz,maxtab),
     *        ilong(maxcnt,2,maxtab),crdixn(maxtab,2)
      logical, intent(in):: tempav

c     internal variables

      integer i_ixn,i_tab,isit1,isit2
      double precision r_dev,ev_force,xi

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     zero force and energy

      f_cord=0.0D0
      E=0.0D0
      xi=0.5D0

      do i_tab=1,4
      do i_ixn=1,numlng(nmres,i_tab)
        isit1=ilong(i_ixn,1,i_tab)
        isit2=ilong(i_ixn,2,i_tab)
        if (ires(isit1).eq.8 .and. (i_tab.eq.3.or. i_tab.eq.4) ) cycle
        if (ires(isit2).eq.8 .and. (i_tab.eq.2.or. i_tab.eq.4) ) cycle
        if( (isit2-isit1).eq.1 ) cycle
!no beta-beta term on glycines

        r_dev=ccev_dist(i_ixn,i_tab)-distne(i_ixn,i_tab)
        E(1,9)=E(1,9)+0.5D0*pexcld*(1.0D0+tanh((r_dev)/xi))
        ev_force=(pexcld/(2.0D0*xi))*(cosh((r_dev)/xi)**(-2))/distne(i_ixn,i_tab)

          f_cord(isit1,1,crdixn(i_tab,1))=
     *          f_cord(isit1,1,crdixn(i_tab,1)) + ev_force*xdiff(i_ixn,i_tab)
          f_cord(isit1,2,crdixn(i_tab,1))=
     *          f_cord(isit1,2,crdixn(i_tab,1)) + ev_force*ydiff(i_ixn,i_tab)
          f_cord(isit1,3,crdixn(i_tab,1))=
     *          f_cord(isit1,3,crdixn(i_tab,1)) + ev_force*zdiff(i_ixn,i_tab)

          f_cord(isit2,1,crdixn(i_tab,2))=
     *          f_cord(isit2,1,crdixn(i_tab,2)) - ev_force*xdiff(i_ixn,i_tab)
          f_cord(isit2,2,crdixn(i_tab,2))=
     *          f_cord(isit2,2,crdixn(i_tab,2)) - ev_force*ydiff(i_ixn,i_tab)
          f_cord(isit2,3,crdixn(i_tab,2))=
     *          f_cord(isit2,3,crdixn(i_tab,2)) - ev_force*zdiff(i_ixn,i_tab)

       enddo
      enddo

      end
