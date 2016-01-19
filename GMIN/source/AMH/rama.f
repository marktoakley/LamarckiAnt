c     ------------------- rama ----------------------
      subroutine rama(ires,jstrt,jfins,pro_cord,f_cord,
     *                ramascl,ramapot,nitcord,cprcord)

c     --------------------------------------------------

c     RAMA figures out chiral forces 

      use amhglobals,  only:maxsiz,maxcrd,aps

      implicit none

      double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd),ramascl,
     *          nitcord(maxsiz,3),cprcord(maxsiz,3)
      double precision, intent(out) :: f_cord(maxsiz,3,maxcrd),ramapot(maxsiz)
      integer, intent(in):: ires(maxsiz),jstrt,jfins
 
c
c
c        a        vector from CA (pro_cord) to C'
c        b        vector from N to CA
c        c        vector from C' to N
c
c        adb        a dot b
c        bdc        b dot c
c        adc        a dot c
c
c        bxa        b cross a
c
c
c
c     sign : Returns the magnitude of its first argument with the sign of its
c            second argument.
c
c     argument declarations:

             double precision 
     *        a(3),b(3),c(3),pa(3),pb(3),pc(3),
     *              pa2,pb2,pc2,padb,padc,pbdc,pbxa(3),
     *              pcdbxa,pchi,pquad,pfac,
     *        adir(0:maxsiz+1,3,2),
     *        bdir(0:maxsiz+1,3,2),w(6),
     *        cdir(0:maxsiz+1,3,2),
     *        a2,b2,c2,adb,adc,bdc,bxa(3),cdbxa,
     *        chi,quad,phi(maxsiz),psi(maxsiz),
     *        slope(2),prob(maxsiz),
     *        pkappa1,pkappa2,
     *        fac,kappa1,kappa2,kappa,pkappa,
     *        cval1,cval2,cval3,nval1,nval2,nval3

c        variables for new potential

            double precision phid(6),psid(6),f(6),sigma(6),s(6)
        integer ips

c     internal variables:
c        --- do loop indices ---
         integer i_axis,i_res,i1,i2,i3

             double precision pi

        data pi /3.141592653589793238462643383279502884197D0/
        data w /1.3149D0,1.17016D0,1.29264D0,1.78596D0,1.0D0,1.0D0/
        data sigma /15.398D0,100.521D0,49.0954D0,419.123D0,0.0D0,0.0D0/
        data phid /-2.051D0,-1.353D0,-1.265D0,-1.265D0,0.0D0,0.0D0/
        data psid /2.138D0,2.4D0,-0.218D0,-0.929D0,0.0D0,0.0D0/

         cval1=0.444D0
         cval2=0.235D0
         cval3=0.321D0

         nval1=0.483D0
         nval2=0.703D0
         nval3=-0.186D0

c     --------------------- begin -----------------------
c     testing for CHRIS

!  zero force and energy
      ramapot(:)=0.0D0
      f_cord=0.0D0

        do 600 i1=1,maxsiz
          phi(i1)=0.0D0
          psi(i1)=0.0D0
600        continue

        do 610 i1=0,maxsiz
          do 620 i2=1,3
            do 630 i3=1,2
              adir(i1,i2,i3)=0.0D0
              bdir(i1,i2,i3)=0.0D0
              cdir(i1,i2,i3)=0.0D0
630            continue
620          continue
610        continue
            
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHI
ccccccccccccccccccccccccccccccccccccccccccccccccc

        do 501 i_res=jstrt+1,jfins-1
c       calculate vectors        

            pa(1)= cprcord(i_res,1)-pro_cord(i_res,1,1)
            pa(2)= cprcord(i_res,2)-pro_cord(i_res,2,1)
            pa(3)= cprcord(i_res,3)-pro_cord(i_res,3,1)

            pb(1)= pro_cord(i_res,1,1)-nitcord(i_res,1)
            pb(2)= pro_cord(i_res,2,1)-nitcord(i_res,2)
            pb(3)= pro_cord(i_res,3,1)-nitcord(i_res,3)

            pc(1)= nitcord(i_res,1)-cprcord(i_res-1,1)
            pc(2)= nitcord(i_res,2)-cprcord(i_res-1,2)
            pc(3)= nitcord(i_res,3)-cprcord(i_res-1,3)

c     calculate dot products 

        padb=pa(1)*pb(1)+pa(2)*pb(2)+pa(3)*pb(3)
        pbdc=pb(1)*pc(1)+pb(2)*pc(2)+pb(3)*pc(3)
        padc=pa(1)*pc(1)+pa(2)*pc(2)+pa(3)*pc(3)

        pbxa(1)=pb(2)*pa(3)-pb(3)*pa(2)
        pbxa(2)=pb(3)*pa(1)-pb(1)*pa(3)
        pbxa(3)=pb(1)*pa(2)-pb(2)*pa(1)

        pcdbxa=pc(1)*pbxa(1)+pc(2)*pbxa(2)+pc(3)*pbxa(3)

        pa2=pa(1)*pa(1)+pa(2)*pa(2)+pa(3)*pa(3)
        pb2=pb(1)*pb(1)+pb(2)*pb(2)+pb(3)*pb(3)
        pc2=pc(1)*pc(1)+pc(2)*pc(2)+pc(3)*pc(3)

        pkappa1=dsqrt( pa2*pb2-padb*padb )
        pkappa2=dsqrt( pb2*pc2-pbdc*pbdc )
        pkappa=pkappa1*pkappa2

        pchi=(padb*pbdc-padc*pb2)/pkappa
c        pchi=DMAX1(DMIN1(pchi,0.99999999D9),-0.99999999D9)

c        if (pchi.gt.1d0) pchi=1d0
c        if (pchi.lt.-1d0) pchi=-1d0

        if (pchi.ge.1.0d0)pchi=0.9999
        if (pchi.le.-1.0d0)pchi=-0.9999

        pquad=DSIGN(1.0D0,pcdbxa)
        phi(i_res)=pquad*DACOS(pchi)

ccccccccccccccccccccccccccccccccccccccccccccccccc
c  PSI
ccccccccccccccccccccccccccccccccccccccccccccccccc

            a(1)=  nitcord(i_res+1,1)-cprcord(i_res,1)
            a(2)=  nitcord(i_res+1,2)-cprcord(i_res,2)
            a(3)=  nitcord(i_res+1,3)-cprcord(i_res,3)

            b(1)=  cprcord(i_res,1)-pro_cord(i_res,1,1)
            b(2)=  cprcord(i_res,2)-pro_cord(i_res,2,1)
            b(3)=  cprcord(i_res,3)-pro_cord(i_res,3,1)

            c(1)=  pro_cord(i_res,1,1)-nitcord(i_res,1)
            c(2)=  pro_cord(i_res,2,1)-nitcord(i_res,2)
            c(3)=  pro_cord(i_res,3,1)-nitcord(i_res,3)

c     calculate dot products

        adb=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
        bdc=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
        adc=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)

        bxa(1)=b(2)*a(3)-b(3)*a(2)
        bxa(2)=b(3)*a(1)-b(1)*a(3)
        bxa(3)=b(1)*a(2)-b(2)*a(1)

        cdbxa=c(1)*bxa(1)+c(2)*bxa(2)+c(3)*bxa(3)

        a2=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
        b2=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
        c2=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)

        kappa1=dsqrt( a2*b2-adb*adb )
        kappa2=dsqrt( b2*c2-bdc*bdc )
        kappa=kappa1*kappa2

        chi=(adb*bdc-adc*b2)/kappa
c        chi=DMAX1(DMIN1(chi,0.999999999999999D9),-0.999999999999999999D9)

c        if (chi.gt.1d0) chi=1d0
c        if (chi.lt.-1d0) chi=-1d0

        if (chi.ge.1.0d0)chi=0.9999
        if (chi.le.-1.0d0)chi=-0.9999

        quad=DSIGN(1.0D0,cdbxa)
        psi(i_res)=quad*DACOS(chi)

        if (psi(i_res).ge. pi) psi(i_res)=-psi(i_res)
        if (phi(i_res).ge. pi) phi(i_res)=-phi(i_res)
        if (ires(i_res).ne.8) then
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEW POTENTIAL%%%%%%%%%%%%%%%%%%

        prob(i_res) = 0.0D0

        do  ips = 1,6
          f(ips) = DEXP(-sigma(ips)*(DCOS(phi(i_res) -
     *     phid(ips)) -1.0D0)**2)
           f(ips) = f(ips)
           s(ips) = DEXP(-sigma(ips)*(DCOS(psi(i_res) -
     *     psid(ips)) -1.0D0)**2)
           s(ips) = s(ips)
        prob(i_res) = prob(i_res) + w(ips)*f(ips)*s(ips)*aps(i_res,ips)

        enddo   ! ips

        ramapot(i_res) = - ramascl*prob(i_res)
        slope(1) = 0.0D0
        slope(2) = 0.0D0
        do  ips = 1,6
           slope(1) = slope(1) +aps(i_res,ips)*
     *     w(ips)*ramascl*2.0D0*sigma(ips)*f(ips)*s(ips)*
     *      (DCOS(phi(i_res)-phid(ips)) -1.0D0)*DSIN(phi(i_res) - phid(ips))
           slope(2) = slope(2) + aps(i_res,ips)*
     *     w(ips)*ramascl*2.0D0*sigma(ips)*f(ips)*s(ips)*
     *      (DCOS(psi(i_res)-psid(ips)) -1.0D0)*DSIN(psi(i_res) - psid(ips))

        enddo  ! ips

        pfac=-pquad*slope(1)/dsqrt(1.0D0-pchi*pchi)
       fac=-quad*slope(2)/dsqrt(1.0D0-chi*chi)

        adir(i_res,1,1)=pfac*((pbdc*pb(1)-pb2*pc(1))/pkappa
     *         -pchi*(pb2*pa(1)-padb*pb(1))/pkappa1**2)
        adir(i_res,2,1)=pfac*((pbdc*pb(2)-pb2*pc(2))/pkappa
     *         -pchi*(pb2*pa(2)-padb*pb(2))/pkappa1**2)
        adir(i_res,3,1)=pfac*((pbdc*pb(3)-pb2*pc(3))/pkappa
     *         -pchi*(pb2*pa(3)-padb*pb(3))/pkappa1**2)

        bdir(i_res,1,1)=pfac
     *         *((padb*pc(1)+pbdc*pa(1)-2.0D0*padc*pb(1))/pkappa
     *         -pchi*(pa2*pb(1)-padb*pa(1))/pkappa1**2
     *         -pchi*(pc2*pb(1)-pbdc*pc(1))/pkappa2**2)
        bdir(i_res,2,1)=pfac
     *         *((padb*pc(2)+pbdc*pa(2)-2.0D0*padc*pb(2))/pkappa
     *         -pchi*(pa2*pb(2)-padb*pa(2))/pkappa1**2
     *         -pchi*(pc2*pb(2)-pbdc*pc(2))/pkappa2**2)
        bdir(i_res,3,1)=pfac
     *         *((padb*pc(3)+pbdc*pa(3)-2.0D0*padc*pb(3))/pkappa
     *         -pchi*(pa2*pb(3)-padb*pa(3))/pkappa1**2
     *         -pchi*(pc2*pb(3)-pbdc*pc(3))/pkappa2**2)

        cdir(i_res,1,1)=pfac*((padb*pb(1)-pb2*pa(1))/pkappa
     *         -pchi*(pb2*pc(1)-pbdc*pb(1))/pkappa2**2)
        cdir(i_res,2,1)=pfac*((padb*pb(2)-pb2*pa(2))/pkappa
     *         -pchi*(pb2*pc(2)-pbdc*pb(2))/pkappa2**2)
        cdir(i_res,3,1)=pfac*((padb*pb(3)-pb2*pa(3))/pkappa
     *         -pchi*(pb2*pc(3)-pbdc*pb(3))/pkappa2**2)

        adir(i_res,1,2)=fac*((bdc*b(1)-b2*c(1))/kappa
     *         -chi*(b2*a(1)-adb*b(1))/kappa1**2)
        adir(i_res,2,2)=fac*((bdc*b(2)-b2*c(2))/kappa
     *         -chi*(b2*a(2)-adb*b(2))/kappa1**2)
        adir(i_res,3,2)=fac*((bdc*b(3)-b2*c(3))/kappa
     *         -chi*(b2*a(3)-adb*b(3))/kappa1**2)

        bdir(i_res,1,2)=fac
     *         *((adb*c(1)+bdc*a(1)-2.0D0*adc*b(1))/kappa
     *         -chi*(a2*b(1)-adb*a(1))/kappa1**2
     *         -chi*(c2*b(1)-bdc*c(1))/kappa2**2)
        bdir(i_res,2,2)=fac
     *         *((adb*c(2)+bdc*a(2)-2.0D0*adc*b(2))/kappa
     *         -chi*(a2*b(2)-adb*a(2))/kappa1**2
     *         -chi*(c2*b(2)-bdc*c(2))/kappa2**2)
        bdir(i_res,3,2)=fac
     *         *((adb*c(3)+bdc*a(3)-2.0D0*adc*b(3))/kappa
     *         -chi*(a2*b(3)-adb*a(3))/kappa1**2
     *         -chi*(c2*b(3)-bdc*c(3))/kappa2**2)

        cdir(i_res,1,2)=fac*((adb*b(1)-b2*a(1))/kappa
     *         -chi*(b2*c(1)-bdc*b(1))/kappa2**2)
        cdir(i_res,2,2)=fac*((adb*b(2)-b2*a(2))/kappa
     *         -chi*(b2*c(2)-bdc*b(2))/kappa2**2)
        cdir(i_res,3,2)=fac*((adb*b(3)-b2*a(3))/kappa
     *         -chi*(b2*c(3)-bdc*b(3))/kappa2**2)

             endif
501        continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
c       compute forces
ccccccccccccccccccccccccccccccccccccccccccccccccc

        do 503 i_axis=1,3
          do 504 i_res=jstrt,jfins

      f_cord(i_res,i_axis,1)=f_cord(i_res,i_axis,1)
     * +cval2*adir(i_res-1,i_axis,1)
     * -(cval2+cval3)*adir(i_res,i_axis,1)
     * +(nval1+nval3)*bdir(i_res,i_axis,1)
     * -(nval1)*bdir(i_res+1,i_axis,1)
     * +(nval2-cval2)*cdir(i_res,i_axis,1)
     * +(nval1-cval1)*cdir(i_res+1,i_axis,1)
     * +(nval2-cval2)*adir(i_res-1,i_axis,2)
     * +(nval1-cval1)*adir(i_res,i_axis,2)
     * +cval2*bdir(i_res-1,i_axis,2)
     * -(cval2+cval3)*bdir(i_res,i_axis,2)
     * +(nval1+nval3)*cdir(i_res,i_axis,2)
     * -nval1*cdir(i_res+1,i_axis,2)

       f_cord(i_res,i_axis,3)=f_cord(i_res,i_axis,3)
     * +cval3*adir(i_res,i_axis,1)
     * +(-nval3)*bdir(i_res+1,i_axis,1)
     * +(-cval3+nval3)*cdir(i_res+1,i_axis,1)
     * +(-cval3+nval3)*adir(i_res,i_axis,2)
     * +cval3*bdir(i_res,i_axis,2)
     * +(-nval3)*cdir(i_res+1,i_axis,2)

504      continue
503    continue

c     ---------------------- done -----------------------

      return
      end
