c     ------------------- oxy ----------------------

      subroutine oxy(maxtab,ires, jstrt,jfins,pro_cord,f_cord_oxy,
     *               eqdist,oxscl,chrlscl,ramascl,
     *               oxexcldv,numlng,ilong,nmres,E)
 
c     --------------------------------------------------

c     OXY includes oxygen coordinates 
c     arguments:

c        jstrt - first modified site (i)
c        jfins - last modifies site (i)
c        f_cord_oxy- work array of forces (o)
c        AMHmaxsiz- maximum number of residues (i)

c     ---------------------------------------------------

      use amhglobals, only:maxsiz,maxcnt,i_rama,ooev_dist,oexcldscl,maxcrd
      use amh_interfaces, only:rama

      implicit none
      double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd),eqdist(:),oxscl,chrlscl,ramascl
      double precision, intent(out) :: f_cord_oxy(maxsiz,3,maxcrd),E(:,:)

      integer, intent(in):: maxtab,ires(maxsiz),jstrt,jfins,numlng(0:maxsiz,maxtab),ilong(maxcnt,2,4),nmres
      logical, intent(in):: oxexcldv
 
c     internal variables:
      double precision a(3), b(3), c(3),anorm, bnorm, cnorm,amnr(3),bmnr(3),cmnr(3),
     * aprl(maxsiz,3),bprl(maxsiz,3),cprl(maxsiz,3),norm,eqprod,f_cord_oxy_rama(maxsiz,3,maxcrd),
     * nitcord(maxsiz,3),cprcord(maxsiz,3),delta(0:maxsiz+1,6,3),dist(maxsiz,6),
     * fac1(0:maxsiz+1,6), prod(maxsiz),ramapot(maxsiz),oxf,r_dev,oxdiff(maxcnt),
     * oydiff(maxcnt),ozdiff(maxcnt),oxdist(maxcnt),cval1,cval2,cval3,nval1,nval2,nval3,xi,exvmin_mcp,factor(3),diff(3),
     * distmcp

c        --- do loop indices ---
         integer i501,i503,i504,i505,i506,i507,i_ixn,i_axis,i_res,isit1,isit2
         integer i,j

c     --------------------- begin -----------------------

c     zero force and energy
      f_cord_oxy_rama=0.D0
      f_cord_oxy=0.0D0
      E(1,3)=0.0D0
      E=0.0D0

       do 610 i505=jstrt,jfins
        do 606 i506=1,3
          fac1(i505,i506)=0.0D0
606     continue
610     continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Calculate Chain Connectivity
c        between N-Cb, Cprime-Cb, N-Cprime 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     loop over all amino acids in protein 

       do 500 i_res=jstrt+1,jfins-1

c     calculate nitrogen, cprime position

         cval1=0.444D0
         cval2=0.235D0
         cval3=0.321D0

       cprcord(i_res,1)=cval1*pro_cord(i_res,1,1)+cval2*pro_cord(i_res+1,1,1)+cval3*pro_cord(i_res,1,3)
       cprcord(i_res,2)=cval1*pro_cord(i_res,2,1)+cval2*pro_cord(i_res+1,2,1)+cval3*pro_cord(i_res,2,3)
       cprcord(i_res,3)=cval1*pro_cord(i_res,3,1)+cval2*pro_cord(i_res+1,3,1)+cval3*pro_cord(i_res,3,3)

          nval1=0.483D0
          nval2=0.703D0
          nval3=-0.186D0

       nitcord(i_res,1)=nval1*pro_cord(i_res-1,1,1)+nval2*pro_cord(i_res,1,1)+nval3*pro_cord(i_res-1,1,3)
       nitcord(i_res,2)=nval1*pro_cord(i_res-1,2,1)+nval2*pro_cord(i_res,2,1)+nval3*pro_cord(i_res-1,2,3)
       nitcord(i_res,3)=nval1*pro_cord(i_res-1,3,1)+nval2*pro_cord(i_res,3,1)+nval3*pro_cord(i_res-1,3,3)
          
c     calculate distances        

          delta(i_res,1,1)=pro_cord(i_res,1,2)-nitcord(i_res,1)
          delta(i_res,1,2)=pro_cord(i_res,2,2)-nitcord(i_res,2)
          delta(i_res,1,3)=pro_cord(i_res,3,2)-nitcord(i_res,3)
          delta(i_res,2,1)=pro_cord(i_res,1,2)-cprcord(i_res,1)
          delta(i_res,2,2)=pro_cord(i_res,2,2)-cprcord(i_res,2)
          delta(i_res,2,3)=pro_cord(i_res,3,2)-cprcord(i_res,3)
          delta(i_res,3,1)=cprcord(i_res,1)-nitcord(i_res,1)
          delta(i_res,3,2)=cprcord(i_res,2)-nitcord(i_res,2)
          delta(i_res,3,3)=cprcord(i_res,3)-nitcord(i_res,3)

          dist(i_res,1)=dsqrt( delta(i_res,1,1)**2+delta(i_res,1,2)**2+delta(i_res,1,3)**2 )
          dist(i_res,2)=dsqrt( delta(i_res,2,1)**2+delta(i_res,2,2)**2+delta(i_res,2,3)**2 )
          dist(i_res,3)=dsqrt( delta(i_res,3,1)**2+delta(i_res,3,2)**2+delta(i_res,3,3)**2 )

          fac1(i_res,1)=2.0D0*oxscl*(eqdist(1)-dist(i_res,1))/dist(i_res,1)
          fac1(i_res,2)=2.0D0*oxscl*(eqdist(2)-dist(i_res,2))/dist(i_res,2)
          fac1(i_res,3)=2.0D0*oxscl*(eqdist(3)-dist(i_res,3))/dist(i_res,3)

          if (ires(i_res).eq.8) then
            fac1(i_res,1)=0.0D0
            fac1(i_res,2)=0.0D0
          endif
500    continue

       do 504 i504=1,3
        cprcord(jstrt,i504)=cval1*pro_cord(jstrt,i504,1)+cval2*pro_cord(jstrt+1,i504,1)+cval3*pro_cord(jstrt,i504,3)
        nitcord(jfins,i504)=nval1*pro_cord(jfins-1,i504,1)+nval2*pro_cord(jfins,i504,1)+nval3*pro_cord(jfins-1,i504,3)
 
          delta(jstrt-1,1,i504)=0.0D0
          delta(jstrt-1,2,i504)=0.0D0
          delta(jstrt-1,3,i504)=0.0D0

          delta(jstrt,1,i504)=0.0D0
          delta(jstrt,2,i504)=pro_cord(jstrt,i504,2)-cprcord(jstrt,i504)
          delta(jstrt,3,i504)=0.0D0
          delta(jfins,1,i504)=pro_cord(jfins,i504,2)-nitcord(jfins,i504)
          delta(jfins,2,i504)=0.0D0
          delta(jfins,3,i504)=0.0D0

          delta(jfins+1,1,i504)=0.0D0
          delta(jfins+1,2,i504)=0.0D0
          delta(jfins+1,3,i504)=0.0D0
504     continue

        do 506 i506=1,3
          fac1(jstrt-1,i506)=0.0D0
          fac1(jfins+1,i506)=0.0D0
506     continue

         dist(jstrt,2)=dsqrt( delta(jstrt,2,1)**2+delta(jstrt,2,2)**2+delta(jstrt,2,3)**2 )
         dist(jfins,1)=dsqrt( delta(jfins,1,1)**2+delta(jfins,1,2)**2+delta(jfins,1,3)**2 )
         fac1(jstrt,1)=0.0D0
         fac1(jstrt,2)=2.0D0*oxscl*(eqdist(2)-dist(jstrt,2))/dist(jstrt,2)
         fac1(jstrt,3)=0.0D0

         fac1(jfins,1)=2.0D0*oxscl*(eqdist(1)-dist(jfins,1))/dist(jfins,1)
         fac1(jfins,2)=0.0D0
         fac1(jfins,3)=0.0D0

         do 666 i_res = 1, jfins
          if (ires(i_res).eq.8) then
            fac1(i_res,1)=0.0D0
            fac1(i_res,2)=0.0D0
          endif
666      continue

c  Add forces

        do 503 i503=1,3
        do 501 i501=jstrt,jfins

          f_cord_oxy(i501,i503,1)=f_cord_oxy(i501,i503,1)
     *          +fac1(i501,1)*delta(i501,1,i503)*(-nval2)               
     *          +fac1(i501+1,1)*delta(i501+1,1,i503)*(-nval1)
     *          +fac1(i501,2)*delta(i501,2,i503)*(-cval1)               
     *          +fac1(i501-1,2)*delta(i501-1,2,i503)*(-cval2)
     *          +fac1(i501+1,3)*delta(i501+1,3,i503)*(-nval1)     
     *          +fac1(i501,3)*delta(i501,3,i503)*(-nval2+cval1)
     *          +fac1(i501-1,3)*delta(i501-1,3,i503)*(cval2)

          f_cord_oxy(i501,i503,2)=f_cord_oxy(i501,i503,2)
     *          +fac1(i501,1)*delta(i501,1,i503)               
     *          +fac1(i501,2)*delta(i501,2,i503)               

          f_cord_oxy(i501,i503,3)= f_cord_oxy(i501,i503,3)+
     *           fac1(i501+1,1)*delta(i501+1,1,i503)*(-nval3)
     *          +fac1(i501,2)*delta(i501,2,i503)*(-cval3)               
     *          +fac1(i501+1,3)*delta(i501+1,3,i503)*(-nval3)           
     *          +fac1(i501,3)*delta(i501,3,i503)*(cval3)

501           continue
503        continue

c----------------------------------------------------
c  E stores potentials for output at incmov
c---------------------------------------------------

      if (ires(jstrt).ne.8)  then
        E(1,3)=E(1,3)+oxscl*(eqdist(2)-dist(jstrt,2))**2
      endif

      if (ires(jfins).ne.8)  then
       E(1,3)=E(1,3)+oxscl*(eqdist(1)-dist(jfins,1))**2
      endif

       do 505 i505=jstrt+1,jfins-1
         if (ires(i505).ne.8) then
           E(1,3)=E(1,3)+oxscl*((eqdist(1)-dist(i505,1))**2+(eqdist(2)-dist(i505,2))**2+(eqdist(3)-dist(i505,3))**2 )
         else
           E(1,3)=E(1,3)+oxscl*((eqdist(3)-dist(i505,3))**2)
         endif
505    continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Calculate Oxygen Excluded Volume if Needed
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if ( oxexcldv ) then

        xi=0.5D0

         do i_ixn= 1,numlng(nmres,2)
            i = ilong(i_ixn,1,2)
            j = ilong(i_ixn,2,2)
            if((j-i).eq.1) cycle
          diff(:)=pro_cord(i,:,3)-pro_cord(j,:,3)  ! X,Y,Z
          distmcp=dsqrt( diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3) )

          factor(:)=(oexcldscl/(2.0D0*xi))*(cosh((ooev_dist(i_ixn)-distmcp)/xi)**(-2)) *diff(:)/distmcp

          f_cord_oxy(i,:,3)=f_cord_oxy(i,:,3)+factor(:)
          f_cord_oxy(j,:,3)=f_cord_oxy(j,:,3)-factor(:)

          E(1,11)=E(1,11)+0.5D0*oexcldscl*(1.0D0+tanh((ooev_dist(i_ixn)-distmcp)/xi))
c         write(6,*)'ij E',j-i,0.5D0*oexcldscl*(1.0D0+tanh((ooev_dist(i_ixn)-distmcp)/xi)) 
        enddo

       endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     call rama to get ramachandrian potentials

          ramapot=0.0D0

          if (i_rama.eq.1) then
        call rama(ires,jstrt,jfins,pro_cord,f_cord_oxy_rama,ramascl,ramapot,nitcord,cprcord)
          endif
          E(1,2)=E(1,2)+sum(ramapot(jstrt:jfins))
          f_cord_oxy=f_cord_oxy+f_cord_oxy_rama

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute chiral forces
c        eqprod= -2.5      ! was erroneously eqprod=-0.7698

         eqprod= -0.83D0

       do 507 i507=jstrt+1,jfins-1
         if (ires(i507).ne.8) then
c     calculate vectors
            a(1)=pro_cord(i507,1,1)-nitcord(i507,1)
            b(1)=cprcord(i507,1)-pro_cord(i507,1,1)
            c(1)=pro_cord(i507,1,2)-pro_cord(i507,1,1)
            a(2)=pro_cord(i507,2,1)-nitcord(i507,2)
            b(2)=cprcord(i507,2)-pro_cord(i507,2,1)
            c(2)=pro_cord(i507,2,2)-pro_cord(i507,2,1)
            a(3)=pro_cord(i507,3,1)-nitcord(i507,3)
            b(3)=cprcord(i507,3)-pro_cord(i507,3,1)
            c(3)=pro_cord(i507,3,2)-pro_cord(i507,3,1)
c     calculate minors
            amnr(1)=b(2)*c(3)-b(3)*c(2)
            amnr(2)=b(3)*c(1)-b(1)*c(3)
            amnr(3)=b(1)*c(2)-b(2)*c(1)
            bmnr(1)=c(2)*a(3)-c(3)*a(2)
            bmnr(2)=c(3)*a(1)-c(1)*a(3)
            bmnr(3)=c(1)*a(2)-c(2)*a(1)
            cmnr(1)=a(2)*b(3)-a(3)*b(2)
            cmnr(2)=a(3)*b(1)-a(1)*b(3)
            cmnr(3)=a(1)*b(2)-a(2)*b(1)
c     calculate products
            anorm=1.0D0/(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
            bnorm=1.0D0/(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
            cnorm=1.0D0/(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
            norm=dsqrt(anorm*bnorm*cnorm)
            prod(i507)=(a(1)*amnr(1)+a(2)*amnr(2)+a(3)*amnr(3))*norm
c     calculate partial derivitives
            aprl(i507,1)=norm*amnr(1)-a(1)*prod(i507)*anorm
            aprl(i507,2)=norm*amnr(2)-a(2)*prod(i507)*anorm
            aprl(i507,3)=norm*amnr(3)-a(3)*prod(i507)*anorm
            bprl(i507,1)=norm*bmnr(1)-b(1)*prod(i507)*bnorm
            bprl(i507,2)=norm*bmnr(2)-b(2)*prod(i507)*bnorm
            bprl(i507,3)=norm*bmnr(3)-b(3)*prod(i507)*bnorm
            cprl(i507,1)=norm*cmnr(1)-c(1)*prod(i507)*cnorm
            cprl(i507,2)=norm*cmnr(2)-c(2)*prod(i507)*cnorm
            cprl(i507,3)=norm*cmnr(3)-c(3)*prod(i507)*cnorm
          end if
507     continue

c     add potential to accumulator if averages are to be computed

        do 518 i_res=jstrt+1,jfins-1
        if (ires(i_res).ne.8) then
            E(1,4)=E(1,4)+chrlscl*(prod(i_res)-eqprod)**2
        endif
518     continue

c     add terms for each coordinate
 
      do 509 i_axis=1,3
        do 510 i_res=jstrt+1,jfins-1
        if (ires(i_res).ne.8) then
        
          f_cord_oxy(i_res,i_axis,1)=f_cord_oxy(i_res,i_axis,1)
     *        - 2.0D0*chrlscl*(prod(i_res)-eqprod) * ( (nval1+nval3)*aprl(i_res,i_axis)
     *         -(cval2+cval3)*bprl(i_res,i_axis)-cprl(i_res,i_axis) )

          f_cord_oxy(i_res+1,i_axis,1)=f_cord_oxy(i_res+1,i_axis,1)
     *        - 2.0D0*chrlscl*(prod(i_res)-eqprod) * cval2*bprl(i_res,i_axis)

          f_cord_oxy(i_res-1,i_axis,1)=f_cord_oxy(i_res-1,i_axis,1)
     *        + 2.0D0*chrlscl*(prod(i_res)-eqprod) * nval1*aprl(i_res,i_axis)

          f_cord_oxy(i_res,i_axis,2)=f_cord_oxy(i_res,i_axis,2)
     *        - 2.0D0*chrlscl*(prod(i_res)-eqprod) * cprl(i_res,i_axis)

          f_cord_oxy(i_res,i_axis,3)=f_cord_oxy(i_res,i_axis,3)
     *        - 2.0D0*chrlscl*(prod(i_res)-eqprod) * cval3*bprl(i_res,i_axis)

          f_cord_oxy(i_res-1,i_axis,3)=f_cord_oxy(i_res-1,i_axis,3)
     *        - 2.0D0*chrlscl*(prod(i_res)-eqprod)*(-nval3)*aprl(i_res,i_axis) 

        endif

510             continue
509          continue

c     ---------------------- done -----------------------

      return

      end
