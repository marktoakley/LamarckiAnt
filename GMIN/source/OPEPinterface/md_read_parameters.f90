subroutine read_parameters_md()
    use md_defs
    implicit none
    
    integer :: ierror
    character(len=20) :: dummy, dummy2
      
    ! Read the time step (in fs)
    call Get_Environment_Variable('timestep',dummy)
    if (dummy .eq. '') then
       write(*,*) 'Error : timestep is not defined !'
       stop
    else
       read(dummy,*) timestep
    endif
    
    ! Reads the time unit (in fs)
    call Get_Environment_Variable('timeunit',dummy)
    if (dummy .eq. '') then
       write(*,*) 'Error : timeunit is not defined !'
       stop
    else
       read(dummy,*) timeunit
    endif
    
    ! Reads the number of equilibration temperatures
    call Get_Environment_Variable('n_equilibration',dummy)
    if (dummy .eq. '') then
       write(*,*) 'Error : n_equilibration is not defined !'
       stop
    else
       read(dummy,*) n_equilibration
    endif
    
    ! Reads the number of steps between rescaling velocity during equilibration
    call Get_Environment_Variable('n_rescale_v_equil',dummy)
    if (dummy .eq. '') then
       n_rescale_v_equil = 10
    else
       read(dummy,*) n_rescale_v_equil
    endif
    
    ! Simulation time in fs or in number of steps
    call Get_Environment_Variable('simulation_time',dummy)
    if (dummy .eq. '') then
       
       call Get_Environment_Variable('n_production',dummy2)
       if (dummy2 .eq. '') then
          write(*,*) 'Error : Neither Simulation_time nor n_production is not defined !'
          stop
       else
          read(dummy2,*) n_production
          simulation_time = timeunit*timestep*n_production/1.0e6
       endif
    else
       read(dummy,*) simulation_time
       n_production = int ( simulation_time / (timeunit*timestep) * 1.0e6 )
    endif
 
    !Corrects for rotation of the molecule every x steps
    call Get_Environment_Variable('n_correct_rotation',dummy)
    if (dummy .eq. '') then
       n_correct_rotation = 0
    else
       read(dummy,*) n_correct_rotation
    endif

    if (PBC .and. n_correct_rotation .ne. 0 )then
      write(*,*) '#############       WARNING     #######################'
      write(*,*) ' WARNING ! you have periodic boundary condition and    '
      write(*,*) ' correct_rotation at the same time!                    '
      write(*,*) '#######################################################'
    endif

    
    ! Same thing of the center of mass
    call Get_Environment_Variable('n_correct_center_of_mass',dummy)
    if (dummy .eq. '') then
       n_correct_com = 0
    else
       read(dummy,*) n_correct_com
    endif

    if (PBC .and. n_correct_com .ne. 0 )then
      write(*,*) '#############       WARNING     #######################'
      write(*,*) '  WARNING ! you have periodic boundary condition and   '
      write(*,*) '  center of mass correction at the same time!          '
      write(*,*) '#######################################################'
    endif
    
    ! Accumulates statistics every n steps
    call Get_Environment_Variable('n_statistics',dummy)
    if (dummy .eq. '') then
       write(*,*) 'Error : n_statistics is not defined !'
       stop
    else
       read(dummy,*) n_stats
    endif
    
    ! And computes the average every n steps
    call Get_Environment_Variable('n_stats_average',dummy)
    if (dummy .eq. '') then
       write(*,*) 'Error : n_stats_average is not defined !'
       stop
    else
       read(dummy,*) n_stats_average
    endif
    
    ! Defines the thermostat
    thermo = .false.
    call Get_Environment_Variable('thermostat',dummy)
    if (dummy .eq. '') then
       thermostat = 'none'
    else
       read(dummy,*) thermostat
       thermo = (thermostat=='berendsen' .or. thermostat=='Nose-Hoover' .or. thermostat=="langevin")
       if (.not. thermo) thermostat = 'none'
    endif
    
    ! If Berendsen is chosen,  reads the Berendsen time. 
    if (thermostat.eq.'berendsen') then
       call Get_Environment_Variable('berendsen_time',dummy)
       if (dummy .eq. '') then
          write(*,*) 'Error : berendsen_time (in steps) is not defined !'
          stop
       else
          read(dummy,*) berendsen_time
       endif
    endif
    
    ! If langevin is chosen,  reads the friction coefficient.
    if (thermostat.eq.'langevin') then
       call Get_Environment_Variable('friction_coef',dummy)
       if (dummy .eq. '') then
          write(*,*) 'Error : Langevin thermostat selected but friction_coef is not defined !'
          stop
         !friction_coef = 0.01
       else
          read(dummy,*) friction_coef
       endif
    endif
    
    ! Writes  the configurations every n steps
    call Get_Environment_Variable('n_save_configs',dummy)
    if (dummy .eq. '') then
       n_save_configs = 1000
    else
       read(dummy,*) n_save_configs
    endif
    
    ! And save the restart information
    call Get_Environment_Variable('n_save_restarts',dummy)
    if (dummy .eq. '') then
       n_save_restarts = n_save_configs
    else
       read(dummy,*) n_save_restarts
    endif

    ! By default, save the structures after removing rotation. It is possible to save them
    ! unrotated, however
    call Get_Environment_Variable('save_unrotated_structures',dummy)
    if (dummy .eq. '') then
      save_unrotated_structures = .false.
    else
      read(dummy,*) save_unrotated_structures
    endif
    
    ! How many steps during each temperature of thermalization
    call Get_Environment_Variable('n_steps_thermalization',dummy)
    if (dummy .eq. '') then
       n_steps_thermalization = 1
    else
       read(dummy,*) n_steps_thermalization
    endif
    
    ! The following three options define various file names
    call Get_Environment_Variable('FINAL', dummy)
    if (dummy .eq. '') then
       FINAL  = 'min'
    else
       read(dummy,*) FINAL
    endif
    
    call Get_Environment_Variable('FILECOUNTER', dummy)
    if (dummy .eq. '') then
       COUNTER  = 'filecounter'
    else
       read(dummy,*) COUNTER
    endif
    
    call Get_Environment_Variable('RESTART_FILE', dummy)
    if (dummy .eq. '') then
       RESTARTFILE = 'restart.dat'
    else
       read(dummy,*) RESTARTFILE
    endif
    
    ! Confines the proteins to a sphere with elastic boundaries
    call Get_Environment_Variable('confining_sphere',dummy)
    if (dummy .eq. '') then
        CONFINING_SPHERE = .false.
    else  
        read(dummy,*) CONFINING_SPHERE
    endif
    
    if (confining_sphere) then
      call Get_Environment_Variable('radius_confining_sphere',dummy)
      if (dummy .eq. '') then
          stop "Must provide a value for the radius of the confining sphere"
      else  
          read(dummy,*) RADIUS_CONFINING_SPHERE
          radius_confining_sphere2 = radius_confining_sphere*radius_confining_sphere
      endif

      call Get_Environment_Variable('confining_force_constant',dummy)
      if (dummy .eq. '') then
          k_confine = 1.0d0
      else
          read(dummy,*) k_confine
      endif
    endif        

    ! Control on hydrogen atomns (if velocity too large, restrict motion)
    call Get_Environment_Variable('Control_hydrogen',dummy)
    if (dummy .eq. '') then 
      control_hydrogen = .true.
    else
      read(dummy,*) control_hydrogen
    endif 

    ! envs for rattle
    call Get_Environment_Variable('RATTLE',dummy)
    if (dummy .eq. '') then
      rattle = .false.
    else 
      read(dummy, *) rattle
      if (rattle) then
        call Get_Environment_Variable('RATTLE_Bond_Length_Tolerance', dummy)
        if (dummy .eq. '') then
          epspos = 1.0d-6
        else
          read(dummy, *) epspos
        end if
        call Get_Environment_Variable('RATTLE_Velocity_Tolerance', dummy)
        if (dummy .eq. '') then
          epsvel = 1.0d-12
        else
          read(dummy, *) epsvel
        end if
        call Get_Environment_Variable('MAX_NUM_RATTLE_ITER', dummy)
        if (dummy .eq. '') then
          nrattle = 200
        else
          read(dummy, *) nrattle
        end if  
      end if
    end if  

    
    ! We write down the various parameters for the simulation
    if(taskid .eq. 0 ) then 
      open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
      
      write(FLOG,*) '**************************************************************************'
      write(FLOG,*) ' '
      write(FLOG,'(1X,A39,A20  )') ' Restart the simulation              : ', restart 
      write(FLOG,'(1X,A39,F12.4)') ' Time unit (in femtoseconds)         : ', timeunit
      write(FLOG,'(1X,A39,F12.4)') ' Basic time step (in time units)     : ', timestep
      write(FLOG,'(1X,A39,F12.4)') ' Simulation time (in ns)             : ', simulation_time
      !write(FLOG,'(1X,A39,I12  )') ' Number of production steps          : ', n_production
      write(FLOG,'(1X,A39,I12  )') ' Number of thermalization steps      : ', n_steps_thermalization
      write(FLOG,'(1X,A39,I12  )') ' Number of equilibration steps       : ', n_equilibration
      if (n_equilibration .gt. 0 ) then
         write(FLOG,'(1X,A39,I12  )') ' Number of equilibration steps       : ', n_rescale_v_equil
      endif
      write(FLOG,'(1X,A39,A12  )') ' Thermostat:                         : ', thermostat
      if (thermostat .eq. 'berendsen') then
         write(FLOG,'(1X,A39,I12  )') ' Berendsen time (in steps)           : ', berendsen_time
      endif
      
      write(flog,*) ' '
      write(FLOG,'(1X,A39,L12    )') ' Use a confining sphere              : ', CONFINING_SPHERE
      if (confining_sphere) then
        write(FLOG,'(1X,A39,F12.4  )') ' Radius of the confining sphere      : ', RADIUS_CONFINING_SPHERE
      endif
      write(flog,*) ' '      
      write(FLOG,'(A39,L12)')   ' Bond constraint with RATTLE         : ', rattle
      write(FLOG,'(A34,e12.6)') ' RATTLE Bond Length Tolerance        : ', epspos
      write(FLOG,'(A34,e12.6)') ' RATTLE Velocity Tolerance           : ', epsvel
      write(FLOG,'(A34,I12)')   ' MAX NUM of Iterations per RATTLE    : ', nrattle
      write(FLOG,'(A39,L12)')   ' Use restricted motion for H         : ', control_hydrogen
      write(flog,*) ' '
      write(FLOG,'(1X,A39,I12  )') ' Correct center-of-mass every n steps: ', n_correct_com
      write(FLOG,'(1X,A39,I12  )') ' Correct rotation every n steps      : ', n_correct_rotation
      write(FLOG,'(1X,A39,I12  )') ' Accumulate statistics every n steps : ', n_stats
      write(FLOG,'(1X,A39,I12  )') ' Average statistics every n steps    : ', n_stats_average
      write(FLOG,'(1X,A39,I12  )') ' Write configurations every n steps  : ', n_save_configs
      write(FLOG,'(1X,A39,I12  )') ' Save restart configs every n steps  : ', n_save_restarts
      write(FLOG,'(1X,A39,l12  )') ' Save unrotated structures           : ', save_unrotated_structures
      write(flog,*) ' '
      write(flog,*) 'Input / Output '
      write(flog,*) '********************* '
      write(flog,'(A39,A16  )')  ' Name of log file                     : ', master_logfile
      write(flog,'(A39,A16  )')  ' Prefix for minima (file)             : ', FINAL
      write(flog,'(A39,A16  )')  ' File with filecounter                : ', counter
      write(flog,'(A39,A16  )')  ' Restart file                         : ', RESTARTFILE
      close(FLOG)
    endif
 
    if (.not.allocated(vel)) allocate(vel(vecsize))
    vx => vel(1:3*natoms:3)
    vy => vel(2:3*natoms:3)
    vz => vel(3:3*natoms:3)
    
end subroutine read_parameters_md
