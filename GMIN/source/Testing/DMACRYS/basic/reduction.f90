subroutine test_cell_parameters(m)
  use cell_reduction
  implicit none
  double precision m(6),p(3)
  double precision lattice(3,3)
  logical execute_tests
  
  print '(A,6F8.3)', 'Testing cell parameters ', m
  call cell_set_angles(m)
  if(execute_tests()) then
    print *,"Test successful"
  else
    print *,"Test failed"
    stop
  endif
  return
end subroutine

function execute_tests()
  use cell_reduction
  use vec3
  implicit none
  logical execute_tests
  double precision lattice1(3,3), dprand
  double precision lattice1inv(3,3)
  double precision lattice2(3,3)
  double precision r(3,3)
  double precision v1(3,3), b1(3,3)
  double precision v2(3,3), b2(3,3)
  integer i, j, n

  execute_tests = .false.

  call cell_get_lattice(lattice1)
  call cell_minimum_reduction()
  call cell_get_lattice(lattice2)
!  print *,lattice1
!  print *,lattice2
  !  print *, "The transformation matrix is"
  !  print *,r_inv_
  r = cell_get_transformation()
  !  print *, "The transformation matrix r is"
  !  print *,r
  do n=1,100
    do i=1,3
      b1(:,i) = (/dprand(), dprand(), dprand()/)
      v1(:,i) = matmul(lattice1,b1(:,i))
      b2(:,i) = matmul(r,b1(:,i))
      !      print *,"b1: ", b1(:,i)
      !      print *,"b2: ", b2(:,i)
      v2(:,i) = matmul(lattice2,b2(:,i))
    end do

    ! now check vector length
    do i=1,3
      do j=i,3
        if(.not.(abs(dot_product(v1(:,i),v1(:,j)) - dot_product(v2(:,i),v2(:,j))) < 1e-6)) then

          return
        endif
      end do
    end do
  end do

  execute_tests = .true.
  return
end function


program test_minimum_reduction
  use cell_reduction
  implicit none
  double precision pi
  double precision dprand
  double precision lattice(3,3)
  parameter (pi=3.141592654d0)
  logical execute_tests
  integer n,i,j

  call SDPRND(100)
  call test_cell_parameters((/1d0,1d0,1d0,90d0,90d0,90d0/))
  call test_cell_parameters((/2d0,3d0,5d0,80d0,90d0,100d0/))
  call test_cell_parameters((/5d0,3d0,2d0,50d0,120d0,130d0/))
  call test_cell_parameters((/0.98d0,0.99*dsqrt(2d0),1d0,45d0,90d0,90d0/))
  call test_cell_parameters((/6.8847d0,9.6634d0,7.3173d0,90d0,90d0,90d0/))

  call test_cell_parameters((/7.3554d0,9.6526d0,6.8981d0,90d0,90d0,93.5056d0/))
  call test_cell_parameters((/5.8624d0,6.3141d0,6.8046d0,89.9745d0,89.9979d0,75.5530d0/))

  lattice(:,1) = (/5.154516032580123d0, 1.8552065440687974D-009,2.108593398352645d0/)
  lattice(:,2) = (/0d0,5.634946737466266d0, -1.2430171254899343D-009/)
  lattice(:,3) = (/0d0,0d0,8.227665973975959d0/)
  call cell_set_lattice(lattice)
  if(.not.execute_tests()) then
      print *,"Test failed"
      stop
  endif

  print '(A)',"Testing 100 random cells"
  do n=1,100
    do i=1,3
      do j=1,3
        lattice(i,j) = dprand()-0.5d0
      end do
    end do
    call cell_set_lattice(lattice)
    if(.not.execute_tests()) then
      print *,"Test failed"
      stop
    endif
  end do
  print *,"Test successful"
  print *,"All tests successful"
end program

