
c     ---------------------  ynorm ----------------------

      subroutine ynorm(maxsiz,jstrt,jfins,observ,rnorm)

c     ---------------------------------------------------

c     YNORM  normalizes the array observ by rnorm

c     arguments:

c        maxsiz- maximum array length
c        jstrt - first site to be included
c        jfins - last site to be included
c        observ- set of observables to be normalized
c        rnorm - normalizing factor

c     ---------------------------------------------------

      implicit none

c      include 'utility'  ! utility file

c     argument declarations:

         integer maxsiz,jstrt,jfins
     
         double precision observ(maxsiz),rnorm

c     internal variables:
          integer i_res


c     --------------------- begin -----------------------

c     normalize observ by rnorm

      do 515 i_res=jstrt,jfins
         observ(i_res)=rnorm*observ(i_res)
  515 continue

c     --------------------- done  -----------------------
      
      return
      end
