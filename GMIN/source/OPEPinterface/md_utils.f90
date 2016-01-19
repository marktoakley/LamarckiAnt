! This set of routines serves for the initialisation 
!
! Copyright Normand Mousseau January 2006

!
! This method computes the force on the Calpha atoms based on set of harmonic constraints
! Defined for a number of fragements
!
module constraints
  use defs
  
  save
  integer :: n_constrained_atoms
  logical, dimension(:), allocatable          :: constrained_atoms

  double precision                                     :: energy_c
  double precision, dimension(:), allocatable, target  :: fc, pc
  double precision, dimension(:), pointer              :: cx, cy, cz
  double precision, dimension(:), pointer              :: fcx, fcy, fcz

  
contains
  subroutine set_harmonic_constraints(nat,posa)
    implicit none
    integer :: i, j, iatom

    integer, intent(in) :: nat
    double precision, dimension(3*nat), target, intent(in) :: posa
    double precision, dimension(:), pointer      :: xx, yy, zz

    xx   => posa(1:3*natoms:3)
    yy   => posa(2:3*natoms:3)
    zz   => posa(3:3*natoms:3)


    if (.not.allocated(constrained_atoms)) allocate(constrained_atoms(natoms))
    if (.not.allocated(fc))  allocate(fc(vecsize))
    if (.not.allocated(pc))  allocate(pc(vecsize))

    cx   => pc(1:3*natoms:3)
    cy   => pc(2:3*natoms:3)
    cz   => pc(3:3*natoms:3)
    
    fcx   => fc(1:3*natoms:3)
    fcy   => fc(2:3*natoms:3)
    fcz   => fc(3:3*natoms:3)

    iatom = 0
    do i = 1, nfrag
      do j =1, list_fragments(i,2)
        iatom = iatom + 1
        
        if ((list_fragments(i,3) .eq. 1)    &
          &    .and.(atomic_type(iatom) .eq. '  CA ')) then
          constrained_atoms(iatom) = .true.
          cx(iatom) = xx(iatom)
          cy(iatom) = yy(iatom)
          cz(iatom) = zz(iatom)
        else
          constrained_atoms(iatom) = .false.
          cx(iatom)=0.0d0
          cy(iatom)=0.0d0
          cz(iatom)=0.0d0
        endif
      end do
    end do
    
  end subroutine set_harmonic_constraints
  
  subroutine harmonic_constraints(nat,posa,forca,t_energy)
    implicit none
  
    integer :: i
    double precision :: r2, dx, dy, dz

    integer, intent(in) :: nat
    double precision, dimension(3*nat), target, intent(in)    :: posa
    double precision, dimension(3*nat), intent(inout) :: forca
    double precision, intent(inout)                   :: t_energy

    double precision, dimension(:), pointer      :: xx, yy, zz

    xx   => posa(1:3*natoms:3)
    yy   => posa(2:3*natoms:3)
    zz   => posa(3:3*natoms:3)

    fc(:) = 0.0d0

    do i=1, natoms
       if (constrained_atoms(i)) then
          ! We do not apply periodic boundary conditions as we want to enforce the constraint
          dx = xx(i) - cx(i)
          dy = yy(i) - cy(i)
          dz = zz(i) - cz(i)
          
          r2 = dx*dx + dy*dy + dz*dz
          energy_c = energy_c + r2
          fcx(i) = dx
          fcy(i) = dy
          fcz(i) = dz
       end if
    end do
    
    fc(:) = -2.0 * k_spring * fc(:)
    energy_c = k_spring * energy_c

    t_energy = energy_c   !  t_energy + energy_c

    forca(:) = forca(:) + fc(:)

  end subroutine harmonic_constraints
end module constraints


module md_utils
  use md_defs
  use geometric_corrections
  use constraints

  implicit none
  
  save
  contains

  subroutine integrate(id,energyscale,constrained_energy)
    integer :: id
    double precision, intent(in) :: energyscale
    double precision :: constrained_energy

    call vv_1()
    ! Compute the force (with constraints, if needed) at this new position
    call calcforce(energyscale,pos,force,total_energy)
    call vv_integrate(constrained_energy)
  end subroutine integrate

  subroutine integrate_log(id,energyscale,constrained_energy,logener,&
  simuPercent)
    integer :: id
    double precision, intent(in) :: energyscale
    double precision :: constrained_energy
    logical :: logener
    double precision :: simuPercent

    call vv_1()
    ! Compute the force (with constraints, if needed) at this new position
    call calcforce_log(energyscale,pos,force,total_energy,logener, simuPercent)
    call vv_integrate(constrained_energy)
  end subroutine integrate_log

  subroutine vv_1()
    vel = vel + force*dt1_mass
    
    ! Call the first part of RATTLE 
    if (rattle) then 
     call rattles_1()
    else
     pos = pos + vel*dt
    endif
  end subroutine vv_1


  ! This routine applies velocity Verlet with the possibility of rattle
  subroutine vv_integrate(constrained_energy)
    double precision :: constrained_energy

    ! begin PLUMED 
!!!RL     if(meta_on.and.meta_start) then
!!!RL       call meta_force_calculation(pos,meta_force)
!!!RL       force(:) = force(:) + meta_force(:) 
!!!RL     endif
    ! end PLUMED 
 
    ! If we have set-up a sphere, apply the confining condition
    if  (confining_sphere)  call confine_in_sphere()

    ! If we constrain the fragments, then add harmonic constraints
    if (constrained_fragments.or.restrained_fragments) then
      call harmonic_constraints(natoms,pos,force,constrained_energy)
    end if
    
    ! Second half of velocity Verlet
     vel = vel + force*dt1_mass

    ! and of RATTLE 
    if (rattle) call rattles_2() 

  end subroutine vv_integrate

  ! Initialise the velocities for the MD
  subroutine initialise_velocities(temperature) 
    implicit none
    
    integer :: i
    double precision, intent(in) :: temperature
    double precision :: rescale, current_temperature
    double precision :: gaussian_number
    
    ! First, generate the velocity according to a Gaussian distribution
    do i=1, vecsize
       vel(i) = gaussian_number()
    enddo

    ! Remove the total momentum and rotation     
    call remove_total_momentum(natoms,vel,mass)
    call remove_rotation(natoms,pos,vel,mass)
    
    ! Then renormalize the velocities to start at the correct temperature
    current_temperature = sum(mass*vel*vel) / degree_freedom ! (VECSIZE -6.0d0)  ! 3 rotations, 3 translations 
    rescale = dsqrt( temperature/current_temperature)
    vel = vel*rescale
    
    return
  end subroutine initialise_velocities
  
  
  !****************************************************************************
  !* This routine rescales the velocities using either a Berendsen bath or
  !* straight renormalisation, depending on the value of "thermal_bath"
  !* 
  !* It also checks on the hydrogen velocities. If they are too large, they are
  !* rescaled to the value (amplitude and direction) of their neighbouring N
  !****************************************************************************
  subroutine scale_velocities(thermal_bath, temperature, nt)
    implicit none
    
    integer :: i, nt
    double precision, intent(in) :: temperature
    double precision :: twokinetic, current_temperature, rescale
    character(len=20), intent(in) :: thermal_bath

    integer :: it, i1, i2, xii
    double precision :: gmaxbf, sdots, g0
    double precision :: vh2, vh2max, vn2
    double precision, dimension(3) :: qij, rij, sij
    double precision pbc_mic
    double precision :: c2
    double precision, external :: gaussian_number
 

    ! If control on hydrogen is one, then we go through the following steps
    if (control_hydrogen) then 
       twokinetic = sum(mass*vel*vel)

       ! (VECSIZE -6.0d0)  ! 3 rotations, 3 translations
       current_temperature = twokinetic / degree_freedom 
       
       ! Defines the maximum velocity acceptable for the hydrogen
       vh2max = 0.03d0*beq2(1)/(dt*dt)
       gmaxbf = 0.0d0
       do it = 1, nbondt
          i1 = ia(it)
          i2 = ib(it)
          
          rij = (/  x(i1),  y(i1),  z(i1) /) - (/  x(i2),  y(i2),  z(i2) /)

          ! Apply periodic boundary conditions, if necessary
          if(PBC)then    
             do xii = 1, 3
                rij(xii) = pbc_mic(rij(xii))
             end do
          endif
          
          rij1(:,it) = rij
          rij2(it) = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
          
          if (it<=nbondh) then
             vh2 = vx(i2)*vx(i2) + vy(i2)*vy(i2) + vz(i2)*vz(i2)
             if (vh2>=vh2max) then
                vn2 = vx(i1)*vx(i1) + vy(i1)*vy(i1) + vz(i1)*vz(i1)
                write(501, "(a4, 3i10, 3f12.6)") " nt ", it, i2, i1, vh2, vn2, &
                     temperature/current_temperature
                
                vx(i2) = vx(i1)
                vy(i2) = vy(i1)
                vz(i2) = vz(i1)
             end if
          end if
          
          qij = (/ vx(i1), vy(i1), vz(i1) /) - (/ vx(i2), vy(i2), vz(i2) /)
          sij = rij + dt*qij
          
          sdots = sij(1)*sij(1) + sij(2)*sij(2) + sij(3)*sij(3)
          g0 = (sdots - beq2(it))*ibeq2(it)
          if ( abs(g0) > abs(gmaxbf) ) then
             gmaxbf = g0
          end if
       end do
    endif
    
    ! Compute the current kinetic energy and temperature
    twokinetic = sum(mass*vel*vel)
    current_temperature = twokinetic / degree_freedom ! (VECSIZE -6.0d0)  ! 3 rotations, 3 translations

    if (thermal_bath .eq. 'berendsen') then
       rescale = dsqrt ( 1.0d0 - (TIMESTEP/dble(BERENDSEN_TIME)) *&
            ( 1.0d0 - (temperature/current_temperature)))
       vel = vel*rescale
    else if (thermal_bath .eq. 'rescale_vel') then
       rescale = dsqrt( temperature/current_temperature)
       vel = vel*rescale
    else if (thermal_bath .eq. 'langevin') then
       ! langevin_scalev may be quite small, so we don't use a square here
       c2  =  (1.d0-langevin_scalev)*(1.d0+langevin_scalev)*temperature
       ! generate the random component according to a Gaussian distribution
       do i=1, vecsize
          vel(i) = vel(i)*langevin_scalev + sqrt(c2*invmass(i))*gaussian_number()
       enddo
!!!RL        c2  = sqrt((1.d0-c1)*temperature)
!!!RL        do i=1,size(mass)
!!!RL         vel(i) = c1*vel(i) + c2*sqrt((c1+1.d0)/mass(i))*gaussian_number()
!!!RL        enddo
    else
       write(*,*) 'Error in scale_velocities - wrong definition of thermal_bath: ', thermal_bath
       write(FLOG,*) 'Error in scale_velocities - wrong definition of thermal_bath: ', thermal_bath
       stop
    endif
    
    
    return
  end subroutine scale_velocities
  
  !****************************************************************************
  !*  Forces the atoms to stay within a sphere by bouncing all atoms that cross
  !****************************************************************************
  subroutine confine_in_sphere()
    implicit none
    integer :: i
    double precision :: r2, xx, yy, zz, inverse_r, parallel_component
    
    
    ! We go through all atoms and identify those that are outside the radius
    do i = 1, natoms
      r2 = x(i)*x(i) + y(i)*y(i) + z(i)*z(i)
      if (r2 .ge. radius_confining_sphere2) then
        ! We compute the inverse radius to normalise the vector from the center
        inverse_r = 1.0d0/dsqrt(r2)
        xx = x(i)*inverse_r
        yy = y(i)*inverse_r
        zz = z(i)*inverse_r
        ! The bouncing will be parallel to the radius
        parallel_component = xx*vx(i) + yy*vy(i) + zz * vz(i)
        
        ! We substract twice the component to reverse the velocity but keep the total momentum
        if(parallel_component > 0) then
          vx(i) = vx(i) - 2.0d0*parallel_component*xx
          vy(i) = vy(i) - 2.0d0*parallel_component*yy
          vz(i) = vz(i) - 2.0d0*parallel_component*zz
        endif
        
        ! Indicates when an atom is bouncing
!!!        write(*,*) 'Bouncing atom: ', i
      endif
    enddo
    return
    
  end subroutine confine_in_sphere


  !****************************************************************************
  !*  First step of rattle
  !****************************************************************************
  subroutine rattles_1()
    implicit none

    integer :: i,it, ic, i1, i2, xii
    double precision :: g0, rij2, sdotr
    double precision :: pbc_mic, RMA, RMB, GAB
    double precision, dimension(3) :: rij, DX
    double precision, dimension(:,:),  allocatable :: sij 
    logical :: done
    logical, dimension (:),   allocatable :: moving, moved

    allocate(moving(natoms), moved(natoms), sij(nbondt,3))

!! initializing stuff
    moving = .false.
    moved  = .true.
    ic     = 0
    done   = .false.

!! saving current positions
    do i=1,nbondt
     i1 = ia(i)
     i2 = ib(i)
     sij(i,1)=x(i1)-x(i2)
     sij(i,2)=y(i1)-y(i2)
     sij(i,3)=z(i1)-z(i2)
     if(PBC)then
        sij(i,1) = pbc_mic(sij(i,1))
        sij(i,2) = pbc_mic(sij(i,2))
        sij(i,3) = pbc_mic(sij(i,3)) 
     endif
    enddo
 
! Evolve positions
    pos = pos + vel*dt

    do while ((.not.done).and.(ic.le.nrattle))

     done = .true. 

      do it = 1, nbondt

        i1 = ia(it)
        i2 = ib(it)

        if(moved(i1).or.moved(i2)) then
       
         rij = (/  x(i1),  y(i1),  z(i1) /) - (/ x(i2),  y(i2),  z(i2) /) 
 
         ! We apply periodic boundary conditions if necessary
         if(PBC)then
            do xii = 1, 3
               rij(xii) = pbc_mic(rij(xii))
            end do
         endif
         
         rij2 = sum(rij(:)**2)
         g0 = beq2(it)-rij2

         if (abs(g0) > epspos) then

          sdotr = sum(sij(it,:)*rij(:))

          RMA=1.d0 / mass(i1)
          RMB=1.d0 / mass(i2)
      
          GAB = g0 / (2.d0 * (RMA + RMB) * sdotr)

          DX(:)  = sij(it,:) * GAB

          x(i1)   = x(i1)  + RMA * DX(1)
          y(i1)   = y(i1)  + RMA * DX(2)
          z(i1)   = z(i1)  + RMA * DX(3)
 
          x(i2)   = x(i2)  - RMB * DX(1)
          y(i2)   = y(i2)  - RMB * DX(2)
          z(i2)   = z(i2)  - RMB * DX(3)

          DX(:) = DX(:)/dt
           
          vx(i1) = vx(i1) + RMA * DX(1)
          vy(i1) = vy(i1) + RMA * DX(2)
          vz(i1) = vz(i1) + RMA * DX(3)

          vx(i2) = vx(i2) - RMB * DX(1)
          vy(i2) = vy(i2) - RMB * DX(2)
          vz(i2) = vz(i2) - RMB * DX(3)

          moving(i1) = .true.
          moving(i2) = .true.
          done       = .false. 

        end if
       endif
      enddo ! end loop on constraint

      do i=1, natoms
       moved(i)  = moving(i)
       moving(i) = .false.
      enddo
    
      ic = ic + 1
     enddo  ! end of iteration

     if(.not.done) then
      write(*,*) "Try to increase MAX_NUM_RATTLE_ITER"
      stop
     endif

     deallocate(moving, moved, sij)

  end subroutine rattles_1

  !****************************************************************************
  !* Second step for rattle
  !****************************************************************************
  subroutine rattles_2()
    implicit none

    integer :: i,it, ic, i1, i2
    double precision :: RMA, RMB, GAB, rdotv
    double precision pbc_mic
    double precision, dimension(3) :: vij, DX
    double precision, dimension(:,:),  allocatable :: sij
    logical :: done
    logical, dimension (:), allocatable :: moving, moved



    allocate(moving(natoms), moved(natoms), sij(nbondt,3))

!! initialization
    moving = .false.
    moved  = .true.
    ic     = 0
    done   = .false.

!! saving current positions
    do i=1,nbondt
     i1 = ia(i)
     i2 = ib(i)
     sij(i,1)=x(i1)-x(i2)
     sij(i,2)=y(i1)-y(i2)
     sij(i,3)=z(i1)-z(i2)
     if(PBC)then
        sij(i,1) = pbc_mic(sij(i,1))
        sij(i,2) = pbc_mic(sij(i,2))
        sij(i,3) = pbc_mic(sij(i,3))
     endif
    enddo

   do while ((.not.done).and.(ic.le.nrattle))  

    done = .true.

    do it = 1, nbondt

      i1 = ia(it)
      i2 = ib(it)

      if(moved(i1).or.moved(i2)) then

       vij = (/ vx(i1), vy(i1), vz(i1) /) - (/ vx(i2), vy(i2), vz(i2) /)

       rdotv = sum(sij(it,:)*vij(:))
       RMA   = 1.d0/mass(i1)
       RMB   = 1.d0/mass(i2)
       GAB   = -rdotv / (( RMA + RMB ) * beq2(it))

       if(abs(rdotv) .gt. epsvel) then

        DX(:) = sij(it,:) * GAB

        vx(i1) = vx(i1) + RMA * DX(1)
        vy(i1) = vy(i1) + RMA * DX(2)
        vz(i1) = vz(i1) + RMA * DX(3)

        vx(i2) = vx(i2) - RMB * DX(1) 
        vy(i2) = vy(i2) - RMB * DX(2)
        vz(i2) = vz(i2) - RMB * DX(3)
        
        moving(i1) = .true.
        moving(i2) = .true.
        done       = .false. 

       endif
      endif

     end do ! end cycle on constraint

     do i=1,natoms
       moved(i)  = moving(i)
       moving(i) = .false.
     enddo  

     ic = ic+1

    enddo ! end of iteration

   if(.not.done) then
      write(*,*) "Try to increase MAX_NUM_RATTLE_ITER"
      stop
   endif

   deallocate(moving, moved, sij)
  
  end subroutine rattles_2

end module md_utils

