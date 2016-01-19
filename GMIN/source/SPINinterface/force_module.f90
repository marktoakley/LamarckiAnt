!Module for the potential and forces!

module force_mod
	
	use main_mod
	use math_mod
	implicit none
	
	contains
	
!--------------------------------------------------------------------------
		subroutine get_gradient(x,grad)

			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, dimension(n), intent(out) :: grad
		
		!Beginning of the routine:
		
			call dd_gradient(x,grad)
		
			return
		
		end subroutine get_gradient

!--------------------------------------------------------------------------	
		subroutine get_energy(x,enrg)

			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, intent(out) :: enrg
		
		!Beginning of the routine:
			
			call dd_energy(x,enrg)
		
			return
		
		end subroutine get_energy

!--------------------------------------------------------------------------	
		
		subroutine get_hessian(x,hessian)

			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, dimension(n,n), intent(out) :: hessian
		
		!Beginning of the routine:
			
			call dd_hessian(x,hessian)
		
			return
		
		end subroutine get_hessian
		
!--------------------------------------------------------------------------

		subroutine dd_gradient(x,grad)
	
			use list_mod
			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, dimension(n), intent(out) :: grad
			!Local Variables
			double precision, dimension(2) :: g1,g2,s2
			double precision :: yy2xx,xx2yy,xy3,term1,term2
			integer :: i,j
			
		!Beginning of the routine:
		
			sinx = sin(x)
			cosx = cos(x)
		
			do i=1,n
		
				g1(1) = (-1.d0)*spin(i)*sinx(i)
				g1(2) = spin(i)*cosx(i)
		
				term1 = 0.d0
				term2 = 0.d0
		
				do j=1,n
				
					s2(1) = spin(j)*cosx(j)
					s2(2) = spin(j)*sinx(j)
				
					yy2xx = yy(i,j) - 2.d0*xx(i,j)
					xx2yy = xx(i,j) - 2.d0*yy(i,j)
					xy3 = 3.d0 * xy(i,j) * (-1.d0)
      		
					term1 = term1 + yy2xx * s2(1) + xy3 * s2(2)
					term2 = term2 + xx2yy * s2(2) + xy3 * s2(1)
		
				enddo
			
				grad(i) = dd_cons * ( term1 * g1(1) + term2 * g1(2) ) &
						- h_ext * dot_product(g1,h_vec)
      
			enddo
    
			return
    
		end subroutine dd_gradient
	
!--------------------------------------------------------------------------

		subroutine dd_energy(x,enrg)

			use list_mod	
			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, intent(out) :: enrg
			!Local Variables
			double precision, dimension(2) :: s1,s2
			double precision :: yy2xx,xx2yy,xy3,term1,term2
			integer :: i,j
			
		!Beginning of the routine:
		
			sinx = sin(x)
			cosx = cos(x)
			enrg = 0.d0
		
			do i=1,n
		
				s1(1) = spin(i)*cosx(i)
				s1(2) = spin(i)*sinx(i)
		
				term1 = 0.d0
				term2 = 0.d0
		
				do j=1,n
				
					s2(1) = spin(j)*cosx(j) 
					s2(2) = spin(j)*sinx(j)
				
					yy2xx = yy(i,j) - 2.d0*xx(i,j)
					xx2yy = xx(i,j) - 2.d0*yy(i,j)
					xy3 = 3.d0 * xy(i,j) * (-1.d0)
      		
					term1 = term1 + yy2xx * s2(1) + xy3 * s2(2)
					term2 = term2 + xx2yy * s2(2) + xy3 * s2(1)

				enddo
			
				enrg = enrg + (dd_cons * ( term1 * s1(1) + term2 * s1(2) ) / 2.d0) &
					 - h_ext * dot_product(s1,h_vec)
				      
			enddo
    
			return
    
		end subroutine dd_energy

!---------------------------------------------------------------------------------------------------
		
		subroutine dd_hessian(x,hessian)

			use list_mod
			implicit none
			!Formal Variables
			double precision, dimension(n), intent(in) :: x
			double precision, dimension(n,n), intent(out) :: hessian
			!Local Variables
			double precision, dimension(2) :: g1,g2,s2,h1
			double precision :: yy2xx,xx2yy,xy3,term1,term2
			integer :: i,j
			
		!Beginning of the routine:

			sinx = sin(x)
			cosx = cos(x)
		
			do i=1,n
			
				h1(1) = (-1.d0)*spin(i)*cosx(i)
				h1(2) = (-1.d0)*spin(i)*sinx(i)
			
				term1 = 0.d0	
				term2 = 0.d0
		
				do j=1,n
				
					s2(1) = spin(j)*cosx(j)
					s2(2) = spin(j)*sinx(j)
				
					yy2xx = yy(i,j) - 2.d0*xx(i,j)
					xx2yy = xx(i,j) - 2.d0*yy(i,j)
					xy3	= 3.d0 * xy(i,j) * (-1.d0)
      		
					term1 = term1 + yy2xx * s2(1) + xy3 * s2(2)
					term2 = term2 + xx2yy * s2(2) + xy3 * s2(1)
	
				enddo
			
				hessian(i,i) = dd_cons * ( term1 * h1(1) + term2 * h1(2) ) &
							 - h_ext * dot_product(h1,h_vec)
      
			enddo
   
			do i=1,n
   	
				g1(1) = (-1.d0)*spin(i)*sinx(i)
				g1(2) = spin(i)*cosx(i)
			
				do j=1,n
				
					if (j .eq. i) then
						cycle
					else
						g2(1) = (-1.d0)*spin(j)*sinx(j)
						g2(2) = spin(j)*cosx(j)
					
						yy2xx = yy(i,j) - 2.d0*xx(i,j)
						xx2yy = xx(i,j) - 2.d0*yy(i,j)
						xy3	= 3.d0 * xy(i,j) * (-1.d0)
      		
						term1 = yy2xx * g2(1) + xy3 * g2(2)
						term2 = xx2yy * g2(2) + xy3 * g2(1)
      		
						hessian(i,j) = dd_cons * (term1*g1(1) + term2*g1(2)) 
					endif
      	
				enddo
     
			enddo
				
			return
	
		end subroutine dd_hessian		
		
!--------------------------------------------------------------------------

end module force_mod
			