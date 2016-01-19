module lanczos_defs
  use defs
  implicit none
  save

  logical :: first_time = .true.
  real(8) :: lanczos_step
  real(8) :: eigenvalue
  real(8), dimension(30) :: eigenvals

  ! Projection direction based on lanczos computations of lowest eigenvalues
  real(8), dimension(:), allocatable :: old_projection, projection
end module lanczos_defs
  
subroutine initialise_lanczos()
  use lanczos_defs

  integer :: ierror
  character(len=20) :: dummy
  
  call Get_Environment_Variable('LANCZOS_STEP', dummy)
  if (dummy .eq. '' ) then
    LANCZOS_STEP = 0.01
  else
    read(dummy,*) LANCZOS_STEP
  endif


  if(taskid .eq. 0 ) then
    open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
    write(flog,*) ' '
    write(flog,*) '****** Lanczos parameters ************************************************'
    write(flog,'(1X,A51,F8.5)')    ' Size of Lanczos step (numerical derivative)     : ', LANCZOS_STEP
    close(flog)
  endif 
  allocate(old_projection(VECSIZE))
  allocate(projection(VECSIZE))
end subroutine


subroutine lanczos(maxvec,new_projection)
  use defs
  use lanczos_defs
  use saddlepoint
  use geometric_corrections
  use random
  implicit none
 
  integer, intent(in) :: maxvec, new_projection
  real(8), dimension( 2 * maxvec -1 ) :: scratcha
  real(8), dimension(maxvec) :: diag
  real(8), dimension(maxvec-1) :: offdiag
  real(8), dimension(maxvec, maxvec) :: vector
  real(8), dimension(VECSIZE, maxvec), target :: lanc
  real(8), dimension(vecsize) :: fohig

  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
 
  integer :: i,k, i_err,ivec
  real(8) :: a1,a0,b2,b1
  real(8) :: excited_energy
  real(8) :: sum2, invsum
  real(8), dimension(VECSIZE) :: newpos,newforce,ref_force
  real(8) :: ran3, scale

  ! We now take the current position as the reference point and will make 
  ! a displacement in a random direction or using the previous direction as
  ! the starting point.

  scale = 1.0

!  type (t_conformations) :: conformation
  
!  scale = conformation%energyscale

!  write(*,*) "scale---lanczos = " , scale

  call calcforce(scale,pos,ref_force,total_energy,fohig)  
  evalf_number = evalf_number + 1
  z0 => lanc(:,1)

  if(new_projection > 0 ) then
    z0 = projection             ! Vectorial operation
    old_projection = projection ! Vectorial operation
  else
    do i=1, VECSIZE
      z0(i) = 0.5d0 - ran3()
    end do

    call geometric_center(natoms,z0)

    if(first_time) then
      old_projection = z0   ! Vectorial operation
      first_time = .false.
    else
      old_projection = projection
    endif 
  endif

  ! We normalize the displacement to 1 total
  sum2 = dot_product(z0,z0)
  invsum = 1.0/sqrt(sum2)
  z0 = z0 * invsum

  newpos = pos + z0 * LANCZOS_STEP   ! Vectorial operation
  call calcforce(scale,newpos,newforce,excited_energy,fohig)
  evalf_number = evalf_number + 1
  
  ! We extract lanczos(1)
  newforce = newforce - ref_force  

  ! We get a0
  a0 = dot_product(z0,newforce)
  diag(1) = a0

  z1 => lanc(:,2)
  z1 = newforce - a0 * z0    ! Vectorial operation

  b1 = dot_product(z1,z1)
  offdiag(1) = sqrt(b1)

  invsum = 1.0d0 / sqrt ( b1 )
  z1 = z1 * invsum           ! Vectorial operation
  
  ! We can now repeat this game for the next vectors
  do ivec = 2, maxvec-1
    z1 => lanc(:,ivec)
    newpos = pos + z1 * LANCZOS_STEP
    call calcforce(scale,newpos,newforce,excited_energy,fohig) 
    evalf_number = evalf_number + 1
    newforce = newforce - ref_force  

    a1 = dot_product(z1,newforce)
    diag(ivec) = a1

    b1 = offdiag(ivec-1)
    z0 => lanc(:,ivec-1)
    z2 => lanc(:,ivec+1)
    z2 = newforce - a1*z1 -b1*z0

    b2 = dot_product(z2,z2)
    offdiag(ivec) = sqrt(b2)
    
    invsum = 1.0/sqrt(b2)
    z2 = z2 * invsum
  end do

  ! We now consider the last line of our matrix
  ivec = maxvec
  z1 => lanc(:,maxvec)
  newpos = pos + z1 * LANCZOS_STEP    ! Vectorial operation
  call calcforce(scale,newpos,newforce,excited_energy,fohig)
  evalf_number = evalf_number + 1
  newforce = newforce - ref_force

  a1 = dot_product(z1,newforce)
  diag(maxvec) = a1

  ! We now have everything we need in order to diagonalise and find the
  ! eigenvectors.

  diag = -1.0d0 * diag
  offdiag = -1.0d0 * offdiag

  ! We now need the routines from Lapack. We define a few values
  i_err = 0


  ! We now reconstruct the eigenvectors in the real space
  ! Of course, we need only the first maxvec elements of vec

  projection = 0.0d0    ! Vectorial operation
  do k=1, maxvec
    z1 => lanc(:,k)
    a1 = vector(k,1)
    projection = projection + a1 * z1   ! Vectorial operation
  end do 
  
  ! The following lines are probably not needed.
  newpos = pos + projection * LANCZOS_STEP   ! Vectorial operation
  call calcforce(scale,newpos,newforce,excited_energy,fohig)
  evalf_number = evalf_number + 1
  newforce = newforce - ref_force

  eigenvalue=diag(1)/LANCZOS_STEP
  do i=1, 4
    eigenvals(i) = diag(i) / LANCZOS_STEP
  end do

  a1 = dot_product(old_projection,projection)
  if(a1<0.0d0) projection = -1.0d0 * projection 
  
  ! write(*,*) 'eigenvalue : ', eigenvalue, 'Lanczos step: ', lanczos_step
end subroutine lanczos
