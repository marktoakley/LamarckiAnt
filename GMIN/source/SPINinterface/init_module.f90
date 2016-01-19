!Module for initializing the lattice. Belonging to the TWIRL-Program for analyzing
!the energy landscape of magnetic nanosystems.

module init_mod

	use main_mod
	use list_mod
	use marsaglias
	implicit none
	
	contains
	
!--------------------------------------------------------------------------

		subroutine init_placeatoms(crystal)
		
			implicit none
			!Formval Variables
			character(len=4) :: crystal
			!Local Variables
			integer :: i
			
		!Beginning of the routine:
		
			if (crystal == 'hexa')  then
	
				do i=1,n	
					if (mod((i-1)/lcons,2) == 0) then
						lattice(i,1) = dble(mod(i-1,lcons))*xdist 		
					elseif (mod((i-1)/lcons,2) == 1) then
						lattice(i,1) = dble(mod(i-1,lcons))*xdist + xdist/2.d0 		
					endif
					lattice(i,2) = dble((i-1)/lcons)*ydist 							
				enddo		

			elseif (crystal == 'quad') then
		
				do i=1,n	
					lattice(i,1) = dble(mod(i-1,lcons))*xdist 							
					lattice(i,2) = dble((i-1)/lcons)*ydist 									
				enddo
			
			endif
		
			if (dorder > 0.d0) then
		
				do i=1,n
					lattice(i,1) = lattice(i,1) + (0.5d0 - gauss(0.5d0,dorder))*xdist
					lattice(i,2) = lattice(i,2) + (0.5d0 - gauss(0.5d0,dorder))*ydist
				enddo
		
			endif
						
			return
	
		end subroutine init_placeatoms

!--------------------------------------------------------------------------

		subroutine init_placespins

			implicit none
			!Loca Variables
			double precision :: number
			integer :: i
			
		!Beginning of the routine:
		
			do i=1,n
		
				call random_number(number)
				if (amrand() > dmiss) then
					spin(i) = 1.d0
				else
					spin(i) = 0.d0
				endif
	
			enddo

			do i=1,n		
				spin(i) = spin(i) + (0.5d0 - gauss(0.5d0,dspin))*spin(i)	
			enddo 

			return

		end subroutine init_placespins

!--------------------------------------------------------------------------

		subroutine init_shiftspins(x)

			implicit none
			!Formal Variables
			double precision, dimension(n), intent(inout) :: x
			!Local Variables
			integer :: i
		
		!Beginning of the routine:
	
			do i=1,n
				x(i) = amrand()*2.d0*pi - pi			
			enddo
	
			return

		end subroutine init_shiftspins	

!--------------------------------------------------------------------------

		subroutine init_createlists

			implicit none
			!Local Variables
			double precision, dimension(6) :: asum
			double precision :: rx,ry,lx,ly
			integer :: i,j,nx,ny
			
		!Beginning of the routine:
		
			xx = 0.d0
			yy = 0.d0
			xy = 0.d0

			if (ewaldsum .eqv. .false.) then
	
				do nx=-xsize,xsize
		
					do ny =-ysize,ysize
			
						lx = nx * xdist * lcons
						ly = ny * ydist * lcons
			
						do i=1,n
				
							do j=1,n
					
								rx = lattice(j,1) + lx - lattice(i,1)
								ry = lattice(j,2) + ly - lattice(i,2)
						
								if ((nx==0).and.(ny==0)) then
						
									if (j==i) then
							
										xx(i,j) = xx(i,j)
										yy(i,j) = yy(i,j)
										xy(i,j) = xy(i,j)
								
									else
									
										xx(i,j) = xx(i,j) + rx * rx / ((sqrt(rx*rx + ry*ry))**5)
										yy(i,j) = yy(i,j) + ry * ry / ((sqrt(rx*rx + ry*ry))**5)
										xy(i,j) = xy(i,j) + rx * ry / ((sqrt(rx*rx + ry*ry))**5)
								
									endif

								else
							
									xx(i,j) = xx(i,j) + rx * rx / ((sqrt(rx*rx + ry*ry))**5)
									yy(i,j) = yy(i,j) + ry * ry / ((sqrt(rx*rx + ry*ry))**5)
									xy(i,j) = xy(i,j) + rx * ry / ((sqrt(rx*rx + ry*ry))**5)
								
								endif
										
							enddo
					
						enddo
				
					enddo
			
				enddo
			
			elseif (ewaldsum .eqv. .true.) then

				do i=1,n
		
					do j=1,n
			
						rx = lattice(j,1) - lattice(i,1)
						ry = lattice(j,2) - lattice(i,2)
			
						call dipsum((lcons)*xdist,(lcons)*ydist,rx,ry,0.d0,asum)
				
						xx(i,j) = asum(2)
						yy(i,j) = asum(3)
						xy(i,j) = asum(4)
			
					enddo
			
				enddo
	
			endif		
		
			return

		end subroutine init_createlists

!--------------------------------------------------------------------------

		subroutine init_rdmnumbers(rlyrnd)
	
			implicit none
			!Formal Variables
			logical, intent(in) :: rlyrnd
			!Local Variables
			double precision :: dptime
			character(len=10) :: time
			integer :: i,seed
			
		!Beginning of the routine:
		
			if (rlyrnd) then				
				call date_and_time(TIME=time)
				read(time,'(F8.0)') dptime
				seed = int(dptime)
				call amrset(seed)
			else
				call amrset(1000)
			endif
		
			return

		end subroutine init_rdmnumbers
		
!--------------------------------------------------------------------------

		subroutine init_lattice(rrdm,lat)
		
			implicit none
			!Formal Variables
			logical, intent(in) :: rrdm
			character(len=4), intent(in) :: lat

		!Beginning of the routine:
					
			call init_rdmnumbers(rrdm)
			call init_placeatoms(lat)
			call init_placespins
			call init_createlists
				
			return
				
		end subroutine init_lattice
		
!--------------------------------------------------------------------------

end module init_mod