
c     --------------------- intstr ----------------------

      subroutine intstr

c     ---------------------------------------------------

c     INTSTR generates initial structures and performs
c            a quick analysis on the Rg and lr/sr

c     ---------------------------------------------------

      use amhglobals,  only: iresrt,movanal,maxsiz,nmres,maxpro,
     *  numpro,maxcrd,numcrd,prcord,oarchv,ires,iseed_amh,quench

      implicit none

c     internal variables:

         integer procnt

c     required subroutines

         external rndcol,restrt    !   restrt

c     --------------------- begin -----------------------

c     procnt tracks the number of initial structures
c     generated

      procnt=0

      if( ( iresrt.and.(.not.movanal) ) .or. quench )then

         call restrt(procnt,maxsiz,nmres,maxpro,
     *               numpro,maxcrd,numcrd,prcord,
     *               oarchv)

      endif

      if( procnt.lt.numpro )then

c        if number of initial configurations found
c        up to this point is less than the required
c        number, then find random coils for
c        the remainder
         call rndcol(nmres,maxpro,numpro,
     *                ires,procnt,numcrd,prcord,
     *                iseed_amh)

      endif


c     ---------------------- done -----------------------

      return
      end
