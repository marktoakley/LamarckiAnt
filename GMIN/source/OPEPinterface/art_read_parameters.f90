subroutine read_parameters_art()
  use defs
  use saddlepoint
  implicit none

  integer :: ierror
  logical :: flag
  character(len=20) :: dummy

  ! Read the maximum number of events
  call Get_Environment_Variable('Max_Number_Events', dummy)
  if (dummy .eq. '') then   
     NUMBER_EVENTS = 1
  else
     read(dummy,*) NUMBER_EVENTS
  endif

  ! Info regarding initial displacement
  call Get_Environment_Variable('Initial_Step_Size', dummy)
  if (dummy .eq. '') then
     INITSTEPSIZE = 0.001             ! Size of initial displacement in Ang. 
  else
     read(dummy,*) INITSTEPSIZE
  endif

  ! Maximum number of iteration in the activation 
  call Get_Environment_Variable('Max_Iter_Activation',dummy)
  if (dummy .eq. '') then
     MAXITER = 10000                
  else
     read(dummy,*) MAXITER
  endif

  ! Maximum number of iteration in the activation 
  call Get_Environment_Variable('Max_Iter_Basin',dummy)
  if (dummy .eq. '') then
     MAXKTER = 100                
  else
     read(dummy,*) MAXKTER
  endif


  ! Number of relaxation perpendicular moves - basin and activation
  call Get_Environment_Variable('Max_Perp_Moves_Basin',dummy)
  if (dummy .eq. '') then
     MAXKPERP = 2       
  else
     read(dummy,*) MAXKPERP
  endif

  call Get_Environment_Variable('Max_Perp_Moves_Activ',dummy)
  if (dummy .eq. '') then
     MAXIPERP = 12       
  else
     read(dummy,*) MAXIPERP
  endif

  ! Increment size - overall scaling (in angstroems) 
  call Get_Environment_Variable('Increment_Size',dummy)
  if (dummy .eq. '') then
     INCREMENT = 0.01
  else
     read(dummy,*) INCREMENT
  endif

  ! Eigenvalue threshold
  call Get_Environment_Variable('Eigenvalue_Threshold',dummy)
  if (dummy .eq. '') then
     write(*,*) 'Error : No eigenvalue threshold provided  (Eigenvalue_Threshold)'
     stop
  else
     read(dummy,*) EIGEN_THRESH
  endif

  ! Force threshold for the perpendicular relaxation
  call Get_Environment_Variable('Force_Threshold_Perp_Rel',dummy)
  if (dummy .eq. '') then
     FTHRESHOLD = 1.0
  else
     read(dummy,*) FTHRESHOLD
  endif
  FTHRESH2 = FTHRESHOLD * FTHRESHOLD

  ! Force threshold for the convergence at the saddle point
  call Get_Environment_Variable('Exit_Force_Threshold',dummy)
  if (dummy .eq. '') then
     EXITTHRESH = 1.0
  else
     read(dummy,*) EXITTHRESH
  endif

  ! Initial move for the event
  call Get_Environment_Variable('Type_Initial_Move',dummy)
  if (dummy .eq. '') then
     INIT_MOVE = 'global'
  else
     if ( (dummy .ne. 'global') .and. (dummy .ne. 'local')) then
        write(*,*) 'Error : Type_Initial_Move can only be global or local - current value : ', dummy
        stop
     endif     
     read(dummy,*) INIT_MOVE 
  endif

  call Get_Environment_Variable('Fraction_Fragment_Initally_Moved', dummy)
  if (dummy .eq. '') then
     FRACTIONFRAGMENT = 1.0           ! % du fragment est déplacé aléatoirement
  else
     read(dummy,*) FRACTIONFRAGMENT
     if (INIT_MOVE .eq. 'global') then
        write(*,*) 'Warning: The Fraction Fragment is not used in global moves'
     endif
  endif


  ! Force threshold for the convergence at the saddle point
  call Get_Environment_Variable('Number_Lanczos_Vectors',dummy)
  if (dummy .eq. '') then
     NVECTOR_LANCZOS = 16
  else
     read(dummy,*) NVECTOR_LANCZOS
  endif

  ! Check whether or not the filecounter file exists or not. If yes, reads it;
  ! if not, give a basic value to the counter. 

  call Get_Environment_Variable('FILECOUNTER', dummy)
  if (dummy .eq. '') then
     COUNTER  = 'filecounter'
  else
     read(dummy,*) COUNTER
  endif

  inquire(file=COUNTER, exist=flag)
  if (flag) then
    open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
    read(FCOUNTER,'(A13,I6)') dummy, mincounter
    close(FCOUNTER)
  else
    mincounter = 0
  endif

  ! File names
     
  call Get_Environment_Variable('EVENTSLIST', dummy)
  if (dummy .eq. '') then
     EVENTSLIST   = 'events.list'
  else 
     read(dummy,*) EVENTSLIST
  endif
     
  call Get_Environment_Variable('REFCONFIG', dummy)
  if (dummy .eq. '') then
     REFCONFIG   = 'refconfig.dat'
  else
     read(dummy,*) REFCONFIG
  endif
     
  call Get_Environment_Variable('FINAL', dummy)
  if (dummy .eq. '') then
     FINAL  = 'min'
  else
     read(dummy,*) FINAL
  endif

  call Get_Environment_Variable('SADDLE', dummy)
  if (dummy .eq. '') then
     SADDLE  = 'sad'
  else
     read(dummy,*) SADDLE
  endif

  call Get_Environment_Variable('RESTART_FILE', dummy)
  if (dummy .eq. '') then
     RESTARTFILE = 'restart.dat'
  else
     read(dummy,*) RESTARTFILE
  endif


  ! We write down the various parameters for the simulation
  open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) '**************************************************************************'
  write(FLOG,'(1X,A39,I12)')   ' - Number of events                  : ', NUMBER_EVENTS
  write(FLOG,'(1X,A39,F12.4)') ' - Eigenvalue threshold              : ', EIGEN_THRESH
  write(FLOG,'(1X,A39,F12.4)') ' - Total force threshold (saddle)    : ', EXITTHRESH

  write(FLOG,'(1X,A39,A12)')   ' - Type of initial move              : ', INIT_MOVE
  write(FLOG,'(1X,A39,F12.4)') ' - Initial step size                 : ', INITSTEPSIZE
  if (INIT_MOVE .eq. 'local') then 
     write(FLOG,'(1X,A39,F12.4)') ' - Initial fraction of fragment moved: ', FRACTIONFRAGMENT
  endif
  write(FLOG,'(1X,A39,F12.4)') ' - Increment size                    : ', INCREMENT
  write(FLOG,'(1X,A39,I12)')   ' - Number of vectors computed by Lanzcos         : ', NVECTOR_LANCZOS
  write(FLOG,'(1X,A51,I8)')    ' - Maximum number of iteractions (basin -kter)   : ', MAXKTER
  write(FLOG,'(1X,A51,I8)')    ' - Maximum number of iteractions (activation)    : ', MAXITER
  write(FLOG,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (basin) : ', MAXKPERP
  write(FLOG,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (activ) : ', MAXIPERP
  write(FLOG,'(1X,A51,F8.4)')  ' - Force threshold for perpendicular relaxation  : ', FTHRESHOLD

  write(FLOG,*)
  write(FLOG,*) '**************************************************************************'

  write(flog,*) ' '
  write(flog,*) 'Input / Output '
  write(flog,*) '********************* '
  write(flog,'(A39,A16  )')  ' Name of log file                     : ', master_logfile
  write(flog,'(A39,A16  )')  ' Liste of events                      : ', eventslist
  write(flog,'(A39,A16  )')  ' Reference configuration              : ', refconfig
  write(flog,'(A39,A16  )')  ' Restart file                         : ', restartfile
  write(flog,'(A39,A16  )')  ' Prefix for minima (file)             : ', FINAL
  write(flog,'(A39,A16  )')  ' Prefix for saddle points             : ', SADDLE
  write(flog,'(A39,A16  )')  ' File with filecounter                : ', counter

  close(flog)



  allocate(atom_displaced(NATOMS))
  allocate(initial_direction(VECSIZE))

  call initialise_lanczos()


end subroutine read_parameters_art
