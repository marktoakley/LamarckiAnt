!Module for mathematical operations. Belonging to the TWIRL-Program for analyzing
!the energy landscape of magnetic nanosystems.

module math_mod

	use main_mod
	use list_mod
	implicit none
		
	contains

!--------------------------------------------------------------------------
	
		double precision function math_abs(x)
			
			implicit none
			!Formal Variables
			double precision, dimension(:), intent(in) :: x					
			!Local Variables
			integer :: dims,i
				
		!Beginning of the function:
			math_abs = 0.d0
			dims = size(x)				
			do i=1,dims
				math_abs = math_abs + x(i)*x(i)
			enddo				
			math_abs = sqrt(math_abs)
				
			return
				
		end function math_abs

!--------------------------------------------------------------------------

		double precision function math_dis(x,y)
		
			implicit none
			!Formal Variables
			double precision, intent(in), dimension(:) :: x,y
			!Local Variables
			double precision :: s1,s2
			integer :: i,j
					
		!Start of Function
				
			i = size(x)
			j = size(y)
			math_dis = 0.d0
					
			if (i .ne. j) then
				print *,'VECTORS MUST HAVE THE SAME SIZE IN POLDIS FUNCTION!'
				return
			else				
				do i=1,n
					s1 = abs(x(i)-y(i))
					s2 = 2.d0*pi - abs(x(i) - y(i))					
					math_dis = math_dis + min(s1,s2)*min(s1,s2)
				enddo
			endif
				
			math_dis = sqrt(math_dis)
				
			return
				
		end function math_dis
		
!-------------------------------------------------------------------------------		
		
		subroutine math_rescale(x)
		
			implicit none
			!Formal Variables
			double precision, dimension(:), intent(inout) :: x
			!Local Variables
			integer :: dims,i
				
		!Beginning of the routine:
									
			do i=1,dims
				do while(x(i) .gt. pi)
					x(i) = x(i) - 2.d0*pi
				enddo
				do while(x(i) .lt. -pi)
					x(i) = x(i) + 2.d0*pi
				enddo
			enddo
				
			return
				
		end subroutine math_rescale
		
!--------------------------------------------------------------------------
		
		double precision function math_rms(grad)
			
			implicit none
			!Formal Variables
			double precision, dimension(:), intent(in) :: grad
			!Local Variables
			integer :: dims
				
		!Beginning of the function:
			
			dims = size(grad)
			math_rms = math_abs(grad)
			math_rms = sqrt(math_rms/dble(dims))
				
			return
				
		end function math_rms
		
!-------------------------------------------------------------------------------	
		
		subroutine math_magnet(x,mag)
		
			implicit none			
			!Formal Variables
			double precision, dimension(:), intent(in) :: x
			double precision, dimension(2), intent(out):: mag
			!Local Variables
			integer :: dims,i
					
		!Beginning of the routine:
					
			dims = size(x)
			mag = 0.d0	
				
			do i=1,dims
				mag(1) = mag(1) + spin(i)*cos(x(i))
				mag(2) = mag(2) + spin(i)*sin(x(i))
			enddo
				
			return
				
		end subroutine math_magnet
		
!-------------------------------------------------------------------------------

		double precision function math_round(number,m_prec)
	
			implicit none
			!Formal Variables
			double precision, intent(in) :: number
			integer, intent(in) :: m_prec
			!Local Variables
			character(len=30)	:: m_number
			character(len=10)	:: m_format
				
		!Beginning of the function:
				
			math_round = 0.d0
			
			select case(m_prec)						
				case(1)
					m_format = '(F12.1)'
				case(2)
					m_format = '(F12.2)'
				case(3)
					m_format = '(F12.3)'
				case(4)
					m_format = '(F12.4)'
				case(5)
					m_format = '(F12.5)'
				case(6)
					m_format = '(F12.6)'
				case(7)
					m_format = '(F12.7)'
				case(8)
					m_format = '(F12.8)'
				case(9)
					m_format = '(F15.9)'
				case(10)
					m_format = '(F15.10)'
				case(11)
					m_format = '(F15.11)'
				case(12)
					m_format = '(F15.12)'			
			end select

			write(m_number,m_format) number
			read(m_number,*) math_round

			return
		
		end function math_round
		
!--------------------------------------------------------------------------
					
end module math_mod