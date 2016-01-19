!Module for the a lot of important parameters. Belonging to the TWIRL-Program for analyzing
!the energy landscape of magnetic nanosystems.

module main_mod

	implicit none
	save
	
!Main-Parameters for the TWIRL-Code
	
	!Mathematical constants
	
		double precision, parameter :: pi = 3.141592653589793d0 		
		double precision, parameter :: sqtwo = 1.414213562373095d0		
	
	!Field constants
	
		double precision, parameter :: dd_cons = 214.72588d0 			
		double precision ::	h_ext = 0.d0		 						 
		double precision, dimension(2) :: h_vec = (/1.d0 ,0.d0/)					
				
	!Lattice constants

		integer, parameter :: lcons = 6
		integer, parameter :: n = 36
		integer, parameter :: nats = 12
		integer, parameter :: xsize = 0
		integer, parameter :: ysize = 0
		logical, parameter :: ewaldsum = .true.								
		
		double precision, parameter :: xdist= 5.d0										
		double precision, parameter :: ydist= 5.d0										
		
		double precision, parameter :: spinsize = 2000.d0							
		double precision, parameter :: dorder = 0.05d0
																							
		double precision, parameter :: dmiss = 0.d0										
		double precision, parameter :: dspin = 0.d0								
		
		double precision, dimension(n) :: radius(n)
				
	!Choose the potential here
	
		character(len=6) :: pot = 'dipole'

end module main_mod

!--------------------------------------------------------------------------
!module for all needed vectors

module list_mod

	use main_mod
	implicit none
	save

		double precision, dimension(n) :: sinx 				!sinus-array for each angle
		double precision, dimension(n) :: cosx 				!cosinus-array for each angle
		double precision, dimension(n) :: spin				!spin-value of each particle
		double precision, dimension(n,2) :: lattice			!coordinates of each particle
		double precision, dimension(n,n) :: xx				!values for dd-interaction
		double precision, dimension(n,n) :: yy				!values for dd-interaction
		double precision, dimension(n,n) :: xy				!values for dd-interaction

end module list_mod

