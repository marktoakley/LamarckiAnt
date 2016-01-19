!Module for calculating random numbers

module marsaglias
 	
 	implicit none
	double precision, dimension(97) :: u
	double precision :: c, cd, cm
	integer :: i97, j97
	logical :: set
	data set /.false./

	contains

		subroutine amrset(iseed)
 			
 			implicit none
!   initialization routine for marsaglias random number generator
!   testing:  call amrset with iseed = 54185253.  this should result
!   in is1 = 1802 and is2 = 9373.  call amrand 20000 times, then six
!   more times, printing the six random numbers * 2**24 (ie, 4096*4096)
!   they should be: (6f12.1)
!   6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
!
!   input:

			integer, intent(in) :: iseed
! ... integer seed greater than zero
!
!     internal:
!
			integer :: is1, is2
! ... the two internal seeds used in marsaglia algorithm
			integer :: is1max
! ... max value of first seed (is1), 31328
			integer :: is2max
! ... max value of second seed (is2), 30081
			integer :: i, j, k, l, m
! ... used in generation of u()
			double precision :: s, t
! ... used in generation of u()
			integer :: ii, jj
! ... loop indices
			data is1max, is2max /31328, 30081/
! 
!     --- construct two internal seeds from single unbound seed ---
! 
!         is1 and is2 are quotient and remainder of iseed/is2max.  we add
!         one to keep zero and one results from both mapping to one.
!         max and min functions keep is1 and is2 in required bounds.    
			is1 = max((iseed / is2max)+1, 1)
			is1 = min(is1, is1max)
			is2 = max(1, mod(iseed, is2max)+1)
			is2 = min(is2, is2max)
!
			i = mod(is1/177, 177) + 2
			j = mod(is1    , 177) + 2
			k = mod(is2/169, 178) + 1
			l = mod(is2    , 169)
			do ii = 1, 97
				s = 0.0d0
				t = 0.5d0
				do  jj = 1, 24
					m = mod(mod(i*j, 179)*k, 179)
					i = j
					j = k
					k = m
					l = mod(53*l+1, 169)
					if (mod(l*m, 64) .ge. 32) s = s + t
					t = 0.5d0 * t
				end do
				u(ii) = s
			end do
			c  = 362436.0d0   / 16777216.0d0
			cd = 7654321.0d0  / 16777216.0d0
			cm = 16777213.0d0 / 16777216.0d0
			i97 = 97
			j97 = 33

			set = .true.
			return
		
		end subroutine amrset

  
		double precision function amrand ( ) 
			
			implicit none
!     portable random number generator by george marsaglia
!
!     output:
!
!     amrand:  a random number between 0.0 and 1.0
!     internal:
!
			double precision :: uni
! ... working var. for random number
!
!     this random number generator originally appeared in *toward a universal
!     random number generator* by george marsaglia and arif zaman.  florida
!     state university report: fsu-scri-87-50 (1987)
!
!     it was later modified by f. james and published in *a review of pseudo-
!     random number generators*
!
!     this is claimed to be the best known random number generator available.
!     it passes all of the tests for random number generators and has a
!     period of 2^144, is completely portable (gives bit identical results on
!     all machines with at least 24-bit mantissas in the floating point
!     representation).
!
!     the algorithm is a combination of a fibonacci sequence (with lags of 97
!     and 33, and operation "subtraction plus one, modulo one") and an
!     "arithmetic sequence" (using subtraction).

			if ( .not. set ) then

			endif
			uni = u(i97) - u(j97)
			if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
			u(i97) = uni
			i97 = i97 - 1
			if (i97 .eq. 0) i97 = 97
			j97 = j97 - 1
			if (j97 .eq. 0) j97 = 97
			c = c - cd
			if ( c .lt. 0.0d0 ) c = c + cm
			uni = uni - c
			if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
			amrand = uni
    
			return
		
		end function amrand



		double precision function  gauss(am,sd)

			implicit none

!    this is a version of amrand() that adds the constraint of
!    a gaussian distribution, with mean "am" and standard deviation "sd".
!    output is to variable "gauss". it also requires amrset to have
!    been called first, and "uses up" the same sequence that
!    amrand() does.
			double precision, intent(in)  :: am, sd
			double precision :: a, uni, zero, six
			integer :: i
			data zero,six /0.0d0,6.0d0/

			if ( .not. set ) then
			
			endif
			a = zero
			do  i = 1, 12
				uni = u(i97) - u(j97)
				if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
				u(i97) = uni
				i97 = i97 - 1
				if (i97 .eq. 0) i97 = 97
				j97 = j97 - 1
				if (j97 .eq. 0) j97 = 97
				c = c - cd
				if ( c .lt. 0.0d0 ) c = c + cm
				uni = uni - c
				if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
				a = a + uni
			end do
			gauss=(a-six)*sd+am

			return

		end function gauss


		subroutine gauss_vec(am,sd,v,n)
    
			implicit none
!   this is a version of amrand() that adds the constraint of
!   a gaussian distribution, with mean "am" and standard deviation "sd".
!   output is to variable "v". it also requires amrset to have
!   been called first, and "uses up" the same sequence that
!   amrand() does.
			double precision, dimension(n), intent(out) :: v
			double precision, intent(in) :: sd, am
			double precision :: zero, six, a, uni
			integer, intent(in) :: n
			integer :: iter, i
			data zero,six /0.0d0,6.0d0/
    
			if ( .not. set ) then
			endif

			do iter=1,n
				a = zero
				do i = 1, 12
					uni = u(i97) - u(j97)
					if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
					u(i97) = uni
					i97 = i97 - 1
					if (i97 .eq. 0) i97 = 97
					j97 = j97 - 1
					if (j97 .eq. 0) j97 = 97
					c = c - cd
					if ( c .lt. 0.0d0 ) c = c + cm
					uni = uni - c
					if ( uni .lt. 0.0d0 ) uni = uni + 1.0d0
					a = a + uni
				end do
				v(iter) = (a-six)*sd+am
			end do
			
			return

		end subroutine gauss_vec

end module marsaglias