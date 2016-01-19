! This set of routines serves for the initialisation 
!
! Copyright Normand Mousseau January 2006


module restart_module
  use defs
  use fileio

  implicit none
  save
  contains
 
  subroutine save_restart(current_simulation_time,conformation)
    implicit none  
    integer :: i,ierror
    type(t_conformations), intent(in):: conformation
    integer, intent(in) :: current_simulation_time
    double precision, dimension(NATOMS) :: x,y,z,vx,vy,vz,xp,yp,zp
    character(len=100) fname
    
    x  = conformation%pos(1:3*natoms:3)
    y  = conformation%pos(2:3*natoms:3)
    z  = conformation%pos(3:3*natoms:3)

    xp = conformation%posref(1:3*natoms:3)
    yp = conformation%posref(2:3*natoms:3)
    zp = conformation%posref(3:3*natoms:3)

    vx = conformation%vel(1:3*natoms:3)
    vy = conformation%vel(2:3*natoms:3)
    vz = conformation%vel(3:3*natoms:3)
    
    fname = trim(conformation%path) // trim(RESTARTFILE)
    open(unit=FRESTART,file=fname,status='unknown',action='write',position='rewind',iostat=ierror)
    write(FRESTART,*) current_simulation_time
    write(FRESTART,*)  conformation%id, conformation%counter
    write(FRESTART,*)  conformation%energy, conformation%temperature
    write(FRESTART,'(1x,3(2x,f14.6))') (x(i),y(i),z(i),i=1,natoms)
    if (SIMULATION_METHOD .eq. 'MD') then
      write(FRESTART,'(1x,3(2x,f14.6))') (vx(i),vy(i),vz(i),i=1,natoms)
    end if
    write(FRESTART,'(1x,3(2x,f14.6))') (xp(i),yp(i),zp(i),i=1,natoms)
    close(FRESTART)
    return
  end subroutine save_restart
  
  subroutine read_restart(current_simulation_time, conformation)
    implicit none
    type(t_conformations), intent(inout) :: conformation
    integer :: i, ierror  
    double precision :: current_temperature
    integer, intent(out) :: current_simulation_time
    double precision, dimension(NATOMS) :: x,y,z,vx,vy,vz,xp,yp,zp
    character(len=100) fname

    fname = trim(conformation%path) // trim(RESTARTFILE)
    call testfile(fname, ex=.true.)
    open(unit=FRESTART,file=fname,status='old',action='read',position='rewind',iostat=ierror)
    read(FRESTART,*) current_simulation_time
    read(FRESTART,*) conformation%id, conformation%counter
    read(FRESTART,*) conformation%energy, current_temperature
    read(FRESTART,'(1x,3(2x,f14.6))') (x(i),y(i),z(i),i=1,natoms)
    if (SIMULATION_METHOD .eq. 'MD') then
      read(FRESTART,'(1x,3(2x,f14.6))') (vx(i),vy(i),vz(i),i=1,natoms)
    end if
    read(FRESTART,'(1x,3(2x,f14.6))') (xp(i),yp(i),zp(i),i=1,natoms)
    close(FRESTART)

    conformation%pos(1:3*natoms:3) = x
    conformation%pos(2:3*natoms:3) = y
    conformation%pos(3:3*natoms:3) = z

    conformation%posref(1:3*natoms:3) = xp
    conformation%posref(2:3*natoms:3) = yp
    conformation%posref(3:3*natoms:3) = zp

    conformation%vel(1:3*natoms:3) = vx
    conformation%vel(2:3*natoms:3) = vy
    conformation%vel(3:3*natoms:3) = vz

    if (abs(current_temperature - conformation%temperature) > 1e-9) then
      write(6,*) 'ERROR : Temperature of replica ', conformation%path, 'does not match run parameters'
      write(6,*) 'Restart T:', current_temperature, '  T_Param:',conformation%temperature
    endif
    
    open(unit=FLOG,file=conformation%logfile,status='unknown',action='write',position='append',iostat=ierror)
    write(flog,*) 'Restart at time :', current_simulation_time
    close(flog)
    
    return
  end subroutine read_restart
  
  subroutine save_master_restart(current_simulation_time)
    
    implicit none  
    integer ::ierror
    integer, intent(in) :: current_simulation_time
    
    open(unit=FRESTART,file=RESTARTFILE,status='unknown',action='write',position='rewind',iostat=ierror)
    write(FRESTART,*) current_simulation_time
    close(FRESTART)
    
    return
  end subroutine save_master_restart

  subroutine read_master_restart(current_simulation_time)
    
    implicit none  
    integer ::ierror
    integer, intent(out) :: current_simulation_time

    call testfile(RESTARTFILE, ex=.true.)
    open(unit=FRESTART,file=RESTARTFILE,status='old',action='read',position='rewind',iostat=ierror)
    read(FRESTART,*) current_simulation_time
    close(FRESTART)
    
    return
  end subroutine read_master_restart
  
end module restart_module
