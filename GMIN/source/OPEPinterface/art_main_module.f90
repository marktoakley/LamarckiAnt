module saddlepoint
  implicit none
  save

  real(8) :: INITSTEPSIZE       ! Size of initial displacement in Ang. 
  real(8) :: FRACTIONFRAGMENT   ! Portion du fragment est déplacé aléatoirement
  integer :: MAXITER, MAXKTER 
  integer :: MAXIPERP, MAXKPERP 
  integer :: NVECTOR_LANCZOS
  real(8) :: INCREMENT 
  real(8) :: FTHRESHOLD, FTHRESH2 
  real(8) :: EIGEN_THRESH
  real(8) :: EXITTHRESH 
  character(len=20) :: INIT_MOVE
  character(len=20) :: LOGFILE

  integer, dimension(:), allocatable :: atom_displaced    ! Id of local atoms displaced
  integer                            :: natom_displaced   ! Number of local atoms displaced
  real(8), dimension(:), allocatable :: initial_direction ! Initial move for leaving harmonic
                                                          ! well.
contains

  subroutine calcforce(scale,xa,fa,etot,fohig)
    use defs
    use calcforces

    integer :: i
    real(8) :: scale
    real(8), intent(out) :: etot
    real(8), dimension(vecsize), intent(inout) :: xa
    real(8), dimension(vecsize), intent(out) :: fa, fohig
    real(8), dimension(vecsize) :: xfor
  
!    scale  = 1.0
 

    call calcforce_protein(scale,xa,xfor,etot)
    
    do i=1, natoms
       fa    = xfor
    end do
    fohig = 0

  end subroutine calcforce
end module saddlepoint



