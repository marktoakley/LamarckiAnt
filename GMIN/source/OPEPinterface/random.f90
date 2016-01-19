MODULE random
  ! Random number generator (from "Numerical Recipes").
  ! Returns a uniform random deviate between 0.0 and 1.0.
  ! Set idum to any negative value to initialize or  
  ! reinitialize the sequence.                      

  ! Shared variables
  save
  integer :: idum, inext, inextp
  integer :: iff = 0
  integer, dimension(55) :: ma

  integer :: gaussian_flag
  double precision :: gaussian_number2
  
end module random

double precision function ran3()
      use RANDOM
      implicit none

      integer, parameter :: mbig=1000000000
      integer, parameter :: mseed=161803398
      integer, parameter :: mz=0
      double precision, parameter :: fac=1./mbig

      integer :: i,mj, mk, ii, k

      ! Any large mbig, and any smaller (but still large) mseed can be
      !  substituted for the above values.

      if(idum.lt.0.or.iff.eq.0)then
           iff=1
           mj=mseed-iabs(idum)
           mj=mod(mj,mbig)
           ma(55)=mj
           mk=1
           do i=1,54
             ii=mod(21*i,55)
             ma(ii)=mk
             mk=mj-mk
             if(mk.lt.mz)mk=mk+mbig
             mj=ma(ii)
           enddo
           do k=1,4
             do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz)ma(i)=ma(i)+mbig
             enddo
           enddo
           inext=0
           inextp=31
           idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
   end function ran3


!! Generates a Gaussian distribution - taken from Numerical Recipes 
!double precision function gaussian_number()
!  use random
!  implicit none
!  double precision :: u, v, x, y, q
!  double precision :: ran3
!
!  do
!    u = ran3()
!    v = 1.7156*(ran3()-0.5)
!    x = u - 0.449871
!    y = abs(v) + 0.386595
!    q = x**2 + y*(0.19600*y - 0.25472*x)
!    if (q < 0.27597) then
!      exit
!    else if (q > 0.27846) then
!      continue
!    else if (v**2 < -4*log(u)*u**2) then
!      exit
!    endif
!  enddo
!
!  gaussian_number = v/u
!  return
!
!end function gaussian_number

! Generates a Gaussian distribution - taken from Numerical Recipes 
double precision function gaussian_number()
  use random
  implicit none
  double precision :: fac, rsq, v1, v2
  double precision :: ran3

  if (gaussian_flag == 0 ) then
     do
        v1 = 2.0d0 * (ran3() - 0.50d0)
        v2 = 2.0d0 * (ran3() - 0.50d0)
        rsq = v1*v1 + v2*v2
        if ( (rsq .lt. 1.0d0) .and. (rsq .ne. 0.0d0)) exit
     end do
     fac = dsqrt( -2.0 * dlog(rsq) / rsq)
     gaussian_number2 = v1 * fac
     gaussian_flag = 1
     gaussian_number = v2*fac
     return
  else
     gaussian_flag = 0
     gaussian_number = gaussian_number2
     return
  endif

  gaussian_number = -9999
  return

end function gaussian_number


