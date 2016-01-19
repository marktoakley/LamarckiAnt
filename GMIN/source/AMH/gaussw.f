
c     --------------------- gaussw ----------------------

      subroutine gaussw(maxs,dist,delta,gaussp,rgrid)

c     ---------------------------------------------------

c     GAUSSV generates a guassian curve centered about 
c            dist with a 'standard deviation' of deltr

c     arguments:

c        maxs  - number of points the grid is to be
c                evaluated at (i)
c        dist  - center or average of gaussian; 
c                it is the distance between residues 
c                i and j (i)
c        delta - well width (i)
c        gaussp- gaussian as a function of the r-grid
c                in rgrid (o)
c        rgrid - grid of r points for which the gaussian
c                is to be computed (i)        

c     ---------------------------------------------------

      implicit none

c     argument declarations:

         integer maxs

         double precision dist,gaussp(maxs),rgrid(maxs),delta

c     internal variables:
 
         integer i500,i501,i502,i503

         double precision deltr

c     --------------------- begin -----------------------

c     --- diagnostics ---

c     echo input
 
c      write(oarchv,100)maxs,dist,deltr
c  100 format('gaussw:maxs ',i3,' dist and deltr',
c    *        2(1x,1pe10.3))

c     --- end diagnostics ---

c     find (r - r(mem))**2

      do 500 i500=1,maxs
         gaussp(i500)=-( rgrid(i500) - dist )**2
  500 continue

c     divide by 'variance' of well-width
            
      deltr=0.5/(delta**2)
      do 501 i501=1,maxs
         gaussp(i501)=gaussp(i501)*deltr
  501 continue

c     make sure underflow doesn't occur

      do 502 i502=1,maxs
         gaussp(i502)=max( gaussp(i502),-60.0D0 )
  502 continue

c     compute gaussian

      do 503 i503=1,maxs
         gaussp(i503)=exp( gaussp(i503) )
  503 continue

c     --- diagnostics ---

c     echo gaussian grid
 
c      write(oarchv,101)
c  101 format(/'gaussian grid ')
c      do 505 i505=1,maxs
c         write(oarchv,102)i505,gaussp(i505),rgrid(i505)
c  102    format(i3,2(1x,1pe12.5))
c  505 continue
c      stop

c     --- end diagnostics ---

c     ---------------------- done -----------------------
 
      return
      end
