module md_initialise
  use md_defs
  use md_utils
  use restart_module
contains
  ! Initialisation of the simulation
  subroutine initialise(conformation)
    implicit none
    
    logical :: flag
    character(len=20) :: dummy
    integer :: ierror, i, i3, j3, it
    double precision :: temperature
    type(t_conformations), intent(inout) :: conformation
    double precision :: redum
    double precision :: m1, m2

    ! We first initialise the various parameters
    dt=TIMESTEP*timeunit;
    dt1=0.5d0*dt;
    dt2=0.5d0*dt*dt;
    
    ! We then initialise the potential
    call initialise_potential()

    ! We also include the inverse mass
    if (.not.allocated(invmass))  allocate(invmass(vecsize))
    if (.not.allocated(dt1_mass)) allocate(dt1_mass(vecsize))
    if (.not.allocated(dt2_mass)) allocate(dt2_mass(vecsize))


    invmass = 1.0d0/mass
    dt2_mass = dt2 * invmass
    dt1_mass = dt1 * invmass
    
    call center_of_mass(natoms,pos,mass)


    ! We get the information on the bonds and their ideal length
    call readtop1(nbondh, nbonda,nbonds)
    nbondt = nbondh + nbonda

    ! Defines the various vectors necessary
    if (.not.allocated(ia)) allocate(ia(nbondt), ib(nbondt))
    if (.not.allocated(redu1_mass))  &
    allocate(redu1_mass(nbondt), redu2_mass(nbondt), reduA(nbondt), reduB(nbondt), reduAA(nbondt), reduBB(nbondt))
    if (.not.allocated(beq)) allocate(beq(nbondt), beq2(nbondt), ibeq2(nbondt), rij1(3,nbondt), rij2(nbondt))

    call readtop2(nbondt,ia,ib,beq,beq2) !redu1_mass,redu2_mass)

    ! Set up the variables whether rattle is used or not
    if ((.not. rattle) .or. (nrattle<=0)) then
      degree_freedom = dble(vecsize) - 6.0d0
    else  
      degree_freedom = dble(vecsize) - 6.0d0 - dble(nbondt)
    end if  


    ! Compute the various reduced masses needed for rattle 
    do it = 1, nbondt
      i3 = ia(it)
      j3 = ib(it)
      m1 = mass(i3)
      m2 = mass(j3)
      
      ! If we use a Deuterium mass
#ifdef DMASS  
      if ((m1/m2)>=10.0d0) then
        m2 = m2*2.0d0
      else if ((m2/m1)>=10.0d0) then
        m1 = m1*2.0d0
      end if
#endif

      redum = (m1*m2) / (m1+m2) ! Reduced masses
      reduA(it) = redum/(dt*m1) 
      reduB(it) = redum/(dt*m2) 
  
      redu1_mass(it) = redum / (2.0d0 * dt)
      redu2_mass(it) = redum / beq2(it) 
  
      reduAA(it) = redu2_mass(it)/m1
      reduBB(it) = redu2_mass(it)/m2
        
      ibeq2(it) = 1.0d0 / beq2(it)
    end do
      
    ! Generate the list of temperatures for the thermalization
    ! Each new temperature is 1.5 as high as the previous one
    if (.not.allocated(temperatures)) allocate(temperatures(n_steps_thermalization))
    
    temperature = conformation%temperature


    do i=1, n_steps_thermalization
       temperatures(n_steps_thermalization-i+1) = temperature
       temperature = temperature / 1.5
    end do
    
    ! we then take care of the velocities
    if (restart .eq. 'new') then            ! New set of velocities
       call initialise_velocities(temperature)
       restart_simulation_time = 0
    else if (restart .eq. 'restart') then  ! Reuse the previous velocities and position
       write(501,*) 'On relance'
       call read_restart(restart_simulation_time,conformation)
       n_steps_thermalization = 0         ! Make sure that there is no thermalization
    else if (restart .eq. 'new_velocities') then  ! Generates a new set of velocities
       call read_restart(restart_simulation_time,conformation)
       call initialise_velocities(temperature)
    else
       stop 'Wrong choice for restart'
    endif

    if (thermostat .eq. 'langevin') then
       ! Calculate the random force variance for the langevin thermostat
       langevin_scalev = exp(-friction_coef*dt)
    endif


    ! Check whether or not the filecounter file exists or not. If yes, reads it;
    ! if not, give a basic value to the counter. 
    if (restart .ne. 'restart') then
      inquire(file=COUNTER, exist=flag) 
      if (flag) then
        open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
        read(FCOUNTER,'(A13,I6)') dummy, mincounter
        close(FCOUNTER)
      else
        mincounter = 0
      endif
    endif

    return
  end subroutine initialise

  subroutine initialise_potential()
    implicit none
    
    integer :: i
    double precision, dimension(vecsize) :: xpos    ! Positions ( (x1,y1,z1),(x2,y2,z2)...
    double precision, dimension(natoms)  :: amass   ! 1 over the masses
  
    ! We now call the protein part to get the positions and force
  
    inquire(file="conf_initiale_RNA.pdb", exist=RNA_simulation)
    if (RNA_simulation) then
      call initialise_RNA(PBC,BL,C_M,natoms,xpos,amass,atomic_type,force_scaling_factor, use_qbug)
    else
      call initialise_protein(PBC,BL,C_M,natoms,xpos,amass,atomic_type,force_scaling_factor, use_qbug)
    endif

    ! We reorder the positions so that they are consistent with the rest of the program
    do i = 1, natoms
       x(i) = xpos(3*i-2)
       y(i) = xpos(3*i-1)
       z(i) = xpos(3*i)
    end do
    pos = xpos
   
!    mass(1:natoms) = 1.0d0/amass
    mass(1:3*natoms:3) = 1.0d0/amass
    mass(2:3*natoms:3) = 1.0d0/amass
    mass(3:3*natoms:3) = 1.0d0/amass
  
    mass = mass * 2390.0 ! From amu to kcal/mol fs^2 / A^2 
  end subroutine initialise_potential

end module md_initialise
