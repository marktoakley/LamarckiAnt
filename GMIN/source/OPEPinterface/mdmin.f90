! Subroutine mdmin  !! Identical to Genevieve version Dec04
!
! Minimizes the energy. It first uses a steepest descent algorithm and after 
! several hundred steps then uses a conjugate gradient method.

! This minimization is done with only a minimal knowledge of the physics
! of the problem so that it is portable
!
! The minimization uses a simple steepest descent with variable step size.
!
! This file contains 3 routines
!

MODULE minimization
! This module defines a number of parameters used during the minimization
  implicit none
  save
  
  integer :: MAXSTEP_MIN
  double precision :: MINF 
  double precision :: STEP_DMD
  double precision, parameter :: FRIC = 0.5  !! 
  double precision, parameter :: reps = 1.0e-12 
  double precision  :: MINF2 
 
 contains

    subroutine calcforce(scale,xa,fa,etot)
    use defs
    use calcforces
    double precision :: scale
    double precision, intent(out) :: etot
    double precision, dimension(vecsize) :: xpos, xfor
    double precision, dimension(vecsize), intent(inout) :: xa
    double precision, dimension(vecsize), intent(out) :: fa
    

    ! Call the force field
    if (RNA_simulation) then
      call calcforce_RNA(scale,xa,fa,etot,.true., 0.0d0)
    else
      call calcforce_protein(scale,xpos,xfor,etot)
    endif
    
  end subroutine calcforce
  

END MODULE minimization

subroutine initialise_minimization()
  use defs
  use minimization
  use geometric_corrections

  implicit none
  integer :: ierror
  character(len=20) dummy
  
  call Get_Environment_Variable('Max_Num_Iter_Minimization', dummy)
  if (dummy .eq. '') then
     MAXSTEP_MIN = 45000
  else  
     read(dummy,*) MAXSTEP_MIN
  endif
  
  call Get_Environment_Variable('Force_Threshold_Minimization', dummy)
  if (dummy .eq. '') then
     MINF = 0.01
  else  
     read(dummy,*) MINF
  endif
  MINF2 = MINF*MINF  

  call Get_Environment_Variable('Time_Step_DampedMD', dummy)
  if (dummy .eq. '') then
     STEP_DMD = 0.1
  else  
     read(dummy,*) STEP_DMD
  endif
 
  if(taskid .eq. 0 ) then
    open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror) 
    write(flog,*) ' '
    write(flog,*) '****** Damped MD Minimization  parameters **************************************'
    write(flog,'(1X,A51,I8)')    ' Maximum number of iteration (minimization)      : ', MAXSTEP_MIN
    write(flog,'(1X,A51,F8.5)')    ' Force threshold for minimization                : ', MINF
    write(flog,'(1X,A51,F8.5)')    ' Time step for the damped MD minimization        : ', STEP_DMD
    write(flog,*) ' '
    close(flog)
  endif
end subroutine initialise_minimization


subroutine min_converge(success, logfile, relax_name)
  use defs
  use minimization
  use geometric_corrections
  implicit none
  
  character(len=20), intent(in) :: LOGFILE
  logical, intent(out) :: success
  integer :: iter, i, npart, ierror
  double precision :: ftot, ftot2, VERSTEP, delr, rmsd
  double precision :: dotvf, v2, frac, ref_rmsd, dif
  double precision :: scale
  double precision, dimension(VECSIZE) :: vel

  character(len=100),intent(in),optional :: relax_name
  character(len=100) :: fname
  if(present(relax_name)) then
    fname = trim(relax_name)
  else
    fname = 'relaxation.pdb'
  endif

!  type (t_conformations) :: conformation

!  scale = conformation%energyscale

!  write(*,*) "scale---mdmmin", scale
 
  scale = 1.0
 
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)

  evalf_number = 0   ! Set the number of force evaluations to zero
  

  call calcforce(scale,pos,force,total_energy)
  evalf_number = evalf_number + 1

  ftot2 = dot_product(force,force)
  
  vel = 0.0  ! Vectorial operation
  VERSTEP = STEP_DMD
  do iter=1, MAXSTEP_MIN   
     
     call minimize_rmsd(natoms,posref, pos, delr, rmsd, npart)
     
     ref_rmsd = rmsd
     ! v(t+ 1/2 * dt) = v(t)+(dt/2)*F(t)/m   
     do i=1, vecsize
        vel(i)= vel(i) + (0.5*VERSTEP/MASS(i)) * force(i)
     enddo
     
     ! x(t+dt) = x(t) + dt * v(t+ 1/2* dt)     
     v2 = 0.0
     do i=1, vecsize
        v2 = v2 + vel(i) * vel(i) 
     enddo
     frac = exp(-VERSTEP * FRIC * sqrt(v2))
     
     if (frac .le. 0.0000001 .OR. frac .le. 0.008) then   !! PHIL 14 SEPT 04
        frac = 0.9 
     endif   !!! END NEW PHIL 14 SEPT 04
     
     do i=1, vecsize
        vel(i) = vel(i) * frac
     enddo
     do i=1, vecsize
        pos(i) = pos(i) + VERSTEP * vel(i)/(1.0 + sqrt(ftot2))
     enddo
     
     ! Calculating F(t+dt)
     call calcforce(scale,pos,force,total_energy)
!     write(flog,*) 'total energy for F(t+dt)', total_energy
     evalf_number = evalf_number + 1
     ftot2 = 0.0
     do i=1, vecsize
        ftot2 = ftot2 + force(i) * force(i)
        if (force(i) > 10) then
!      write(98,*) i, force(i)
        endif
     enddo
     
     call minimize_rmsd(natoms, posref, pos, delr, rmsd, npart)
     
     dif = abs(rmsd - ref_rmsd)
!     if (dif < reps ) then
!     endif
     
     if (mod(iter, 50) == 0) then
       if (debug_status .eq. 'debug') write(*, "(' ','iter: ',i5,' energy: ',f12.4,' frac: ', f8.4,       &
             & ' ftot: ', f12.4, ' dotvf: ', f12.6, ' vel: ', f12.4,  '  delr: ',  &       
             & f7.4, ' rmsd: ',f7.3,' npart: ',i4,' verstep: ', f6.3)")  &
             & iter, total_energy, frac, sqrt(ftot2), dotvf, sqrt(v2),   &
             & delr, rmsd, npart, VERSTEP 
        write(flog, "(' ','iter: ',i5,' energy: ',f12.4,' frac: ', f8.4,       &
             & ' ftot: ', f12.4, ' dotvf: ', f12.6, ' vel: ', f12.4,  '  delr: ',  &       
             & f7.4, ' rmsd: ',f7.3,' npart: ',i4,' verstep: ', f6.3)")  &
             & iter, total_energy, frac, sqrt(ftot2), dotvf, sqrt(v2),   &
             & delr, rmsd, npart, VERSTEP 
     endif
     
     
     if ( (ftot2 < MINF2 .or. dif < reps) .and. (iter .gt. 5) ) then
        write(flog,*) sqrt(ftot2), sqrt(MINF2), dif, rmsd, ref_rmsd,reps, iter
        exit
     endif
     
     ! let velocities to be zero if (v dot F) is negative
     dotvf = 0.0
     do i=1, vecsize
        dotvf = dotvf + vel(i) * force(i)
     enddo
     if (dotvf < 0.0) then
        do i=1, vecsize
           vel(i) = 0.0
        enddo
        VERSTEP = 0.95 * VERSTEP
     endif
     
     if (VERSTEP < 0.01) VERSTEP = 1.8 * VERSTEP
     
     ! v(t+dt) = v(t+1/2*dt) + dt/2 * F(t+dt)/m 
     do i=1, vecsize  
        vel(i) = vel(i) + (0.5 * VERSTEP / MASS(i)) * force(i)
     enddo
     if (mod(iter,50) .eq. 0) then
       !call write_to_file('relaxe',0,natoms,pos,posref,ndigits_filename,fname,iter,0.0d0, & 
         !total_energy,.true.,singlefile)
     endif
  enddo
  
  ftot = sqrt(ftot2)
  if (ftot < MINF  .and. iter > 1) then
     success = .true.
     if (debug_status .eq. 'debug') then 
        write(*,*) 'Minimization successful   ftot : ', ftot
        write(*,*) 'Minimization Energy : ', total_energy 
     endif
     write(FLOG,*) 'Minimization successful   ftot : ', ftot
     write(FLOG,*) 'Minimization Energy : ', total_energy 
  else
     success = .false.
     if (debug_status .eq. 'debug') then 
        write(*,*) 'Minimization failed   ftot : ', ftot
        write(*,*) 'Minimization Energy : ', total_energy
     endif
     write(FLOG,*) 'Minimization failed   ftot : ', ftot
     write(FLOG,*) 'Minimization Energy : ', total_energy
  endif
  
  close(flog)
end subroutine min_converge
