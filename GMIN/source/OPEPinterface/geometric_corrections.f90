module geometric_corrections
  
contains
  ! Euclidean norm of a vector
  function euc_norm(v)
    double precision :: euc_norm
    double precision, intent(in) :: v(3)

    euc_norm = dsqrt(dot_product(v,v))
  end function euc_norm
  
  function crossproduct(u, v)
    double precision :: crossproduct(3)
    double precision, intent (in) :: u(3), v(3)
  
    crossproduct(1) = u(2)*v(3) - u(3)*v(2)
    crossproduct(2) = u(3)*v(1) - u(1)*v(3)
    crossproduct(3) = u(1)*v(2) - u(2)*v(1)
  
  end function crossproduct


  subroutine cross(v, v1, v2) 
    implicit none
    double precision, dimension(3), intent(in)  :: v1, v2
    double precision, dimension(3), intent(out) :: v
    
    v(1) = v1(2)* v2(3) - v1(3)*v2(2)
    v(2) = v1(3)* v2(1) - v1(1)*v2(3)
    v(3) = v1(1)* v2(2) - v1(2)*v2(1)
  end subroutine cross
  
  subroutine tensor_product(t, v1, v2, f )
    implicit none
    double precision, intent(in) :: f
    double precision, dimension(3), intent(in) :: v1, v2
    double precision, dimension(3,3), intent(out) :: t
    integer :: i, j
    
    do i=1, 3
       do j= 1, 3
          t(i,j) = f * v1(i) * v2(j)
       end do
    end do
  end subroutine tensor_product
  
  !******************************************************************************************
  ! This subroutine get the rotation matrix using iterative method 
  subroutine rotation_matrix(natoms,posa,posb,M)
    implicit none
    
    integer, intent(in) :: natoms
    double precision, dimension(3*natoms), target :: posa,posb
    double precision, dimension(3,3)  :: cpp,cpp1,cpp2,cpp3    ! Matrix C"  
    double precision, dimension(3, 3) :: M,M1,M2,M3            ! Matrix M 
    double precision, dimension(3, 3) :: matrc                 ! Matrix C 
    double precision, parameter :: eps = 1.0E-12
    double precision :: alpha, beta, gama, s, t, st2, st
    integer :: k, p, q, r, l, iter, maxiter
    
    ! We choose the identity matrix as the initial value of M   
    M = 0.0d0
    
    M(1,1)=1.d0
    M(2,2)=1.d0
    M(3,3)=1.d0
    
    ! We can get the initial value of cpp by cpp=M*matrc, matrc is an constant, the initial value 
    ! of M is an identity matrix 
    
    call matrix_c(natoms,posa,posb,matrc)   ! Calculating the value of matrix C (matrc)
    cpp = matmul(M,matrc)       ! Calculating the value of matrix C" (cpp), cpp=M*matrc
    
    l=1                             ! The initial value of iteration number l 
    ! We make the iteration using "while loop" statement
    ! We change M iteratively by applying the best rotation as the order x-y-z, x-y-z, and so on
    ! until the absolute value of alpha, beta, gama all have become less than a preset tolerance 

    ! Adding an iteration counter to this, since there are cases where this
    ! numeric procedure will not converge
    iter = 0
    maxiter = 1000
    xyz_rotation: do  
       
       ! Read in next value of l
       ! We make the best rotation about x axis once
       p=1
       q=2
       r=3
       
       M1 = M
       cpp1 = cpp
       
       s=cpp1(r,q)-cpp1(q,r)
       t=cpp1(q,q)+cpp1(r,r)
       st2=s*s+t*t
       st=sqrt(st2)
       alpha=atan(s/t)
       do k=1,3
          M(p,k)=M1(p,k)
          M(q,k)=(t*M1(q,k)+s*M1(r,k))/st
          M(r,k)=(-s*M1(q,k)+t*M1(r,k))/st
          cpp(p,k)=cpp1(p,k)
          cpp(q,k)=(t*cpp1(q,k)+s*cpp1(r,k))/st
          cpp(r,k)=(-s*cpp1(q,k)+t*cpp1(r,k))/st  
       end do
       ! We finished the best rotation about x axis once 
       
       ! We make the best rotation about y axis once
       p=2
       q=3
       r=1
       M2 = M
       cpp2 = cpp
       
       s=cpp2(r,q)-cpp2(q,r)
       t=cpp2(q,q)+cpp2(r,r)
       st2=s*s+t*t
       st=sqrt(st2)
       beta=atan(s/t)
       do k=1,3
          M(p,k)=M2(p,k)
          M(q,k)=(t*M2(q,k)+s*M2(r,k))/st
          M(r,k)=(-s*M2(q,k)+t*M2(r,k))/st        
          cpp(p,k)=cpp2(p,k)
          cpp(q,k)=(t*cpp2(q,k)+s*cpp2(r,k))/st
          cpp(r,k)=(-s*cpp2(q,k)+t*cpp2(r,k))/st  
       end do
       ! We finished the best rotation about y-axis once
       
       ! We make the best rotation about z-axis once 
       p=3
       q=1
       r=2
       M3 = M
       cpp3 = cpp
       
       s=cpp3(r,q)-cpp3(q,r)
       t=cpp3(q,q)+cpp3(r,r)
       st2=s*s+t*t
       st=sqrt(st2)
       gama=atan(s/t)
       do k=1,3
          M(p,k)=M3(p,k)
          M(q,k)=(t*M3(q,k)+s*M3(r,k))/st
          M(r,k)=(-s*M3(q,k)+t*M3(r,k))/st        
          cpp(p,k)=cpp3(p,k)
          cpp(q,k)=(t*cpp3(q,k)+s*cpp3(r,k))/st
          cpp(r,k)=(-s*cpp3(q,k)+t*cpp3(r,k))/st  
       end do
       ! We finished the best rotation about y axis once  
       ! At the same time, we finished the best rotation as the order x, y, and z axis once
       
       if ((abs(alpha).le.eps) .and. (abs(beta).le.eps) .and. (abs(gama).le.eps)) exit
       l=l+1

       iter = iter+1
       if (iter .ge. maxiter) then
!         write(*,*) "Exiting after ", maxiter, "iterations."
         ! If this did not converge, there may be junk (i.e, NaN) in the matrix,
         ! so we reset it
         M = 0.0d0
         M(1,1)=1.d0
         M(2,2)=1.d0
         M(3,3)=1.d0

         exit
       endif
       
    end do xyz_rotation
    
  end subroutine rotation_matrix
  
  subroutine minimize_rmsd(nat,posref,pos,delr,rmsd,npart,calpha_in,atomic_type_in)
    ! This is the subroutine of Rotation which gets the rotation matrix using iterative method
    implicit none
    
    double precision, parameter :: THRESHOLD = 0.1          ! In Angstroems
    integer, intent(in) :: nat
    integer ::  natoms

    ! The following two variables are optional if not present, the minimization takes place
    ! over the full set of atoms
    logical, intent(in), optional :: calpha_in
    character(len=5), dimension(nat), intent(in), optional :: atomic_type_in

    logical :: calpha
    character(len=5), dimension(nat) :: atomic_type

    double precision, dimension(3*nat), intent(inout), target :: posref, pos
    double precision, dimension(:), allocatable, target :: posa, posb, pospb

    double precision, dimension(:), pointer :: xa, ya, za   ! Pointers for working position
    double precision, dimension(:), pointer :: xb, yb, zb   ! Pointers for reference position
    double precision, dimension(:), pointer :: xpb, ypb, zpb
    double precision, dimension(3,3) :: M
    double precision :: dr2,dr,delr2, delx,dely,delz,delr, rmsd
    integer :: i,j,npart
    
    if (present(calpha_in)) then 
       calpha = calpha_in
    else
       calpha = .false.
    endif
    if (present(atomic_type_in)) atomic_type(:) = atomic_type_in(:)

    ! We first set up pointers for the x, y, z components 
    if (calpha) then
       natoms = 0
       do i=1, nat
          if (atomic_type(i) .eq. '  CA ') natoms = natoms +1
       enddo
    else
       natoms = nat
    endif

    allocate(posa(3*natoms))
    allocate(posb(3*natoms))
    allocate(pospb(3*natoms))

    if (calpha) then
       j = 0
       do i=1, nat
          if (atomic_type(i) .eq. '  CA ') then
             j = j +1
             posa(3*j-2) =   posref(3*i-2) 
             posa(3*j-1) =   posref(3*i-1) 
             posa(3*j  ) =   posref(3*i  ) 
             posb(3*j-2) =   pos(3*i-2) 
             posb(3*j-1) =   pos(3*i-1) 
             posb(3*j  ) =   pos(3*i  ) 
          endif
       enddo
    else
       posa(:) = posref(:)
       posb(:) = pos(:)
    endif

    xa => posa(1:3*NATOMS:3)
    ya => posa(2:3*NATOMS:3) 
    za => posa(3:3*NATOMS:3)
    
    xb => posb(1:3*NATOMS:3)
    yb => posb(2:3*NATOMS:3)
    zb => posb(3:3*NATOMS:3)
    
    xpb => pospb(1:3*NATOMS:3)
    ypb => pospb(2:3*NATOMS:3)
    zpb => pospb(3:3*NATOMS:3)

    call geometric_center(natoms,posa)
    call geometric_center(natoms,posb)
    
    call rotation_matrix(natoms,posa,posb,M)
    
    do i=1, NATOMS
       xpb(i)=M(1,1)*xb(i) + M(1,2)*yb(i) + M(1,3)*zb(i)
       ypb(i)=M(2,1)*xb(i) + M(2,2)*yb(i) + M(2,3)*zb(i)
       zpb(i)=M(3,1)*xb(i) + M(3,2)*yb(i) + M(3,3)*zb(i)
    end do
    
    !
    ! If compute on the whole protein, then send back the corrected positions
    !
    if (.not.calpha) then
      pos(:) = pospb(:)
    endif

    delr2 = 0.0d0
    npart = 0
    
    do i=1, NATOMS
       delx = (xa(i) - xpb(i))
       dely = (ya(i) - ypb(i))
       delz = (za(i) - zpb(i))
       
       dr2   = delx*delx + dely*dely + delz*delz
       delr2 = delr2 + dr2
       dr = sqrt(dr2) 
       if(dr > THRESHOLD) npart = npart + 1
    end do
    
    delr = sqrt(delr2)
    rmsd = sqrt(delr2/natoms)

    deallocate(posa)
    deallocate(posb)
    deallocate(pospb)
    
  end subroutine minimize_rmsd


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! The subroutine "geometric_center" places the geometric center of a 3D vector
  ! at (0,0,0)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine geometric_center(natoms,vector)
    implicit none
    integer :: natoms
    double precision, dimension(3*natoms), intent(inout), target :: vector
    double precision, dimension(3) :: centrd  ! centroid or center of mass     
    
    double precision, dimension(:), pointer :: x, y, z     ! Pointers for coordinates
    
    ! We first set-up pointers for the x, y, z components 
    x => vector(1:3*natoms:3)
    y => vector(2:3*natoms:3)
    z => vector(3:3*natoms:3)
 
    centrd(1) = sum(x)/natoms
    centrd(2) = sum(y)/natoms
    centrd(3) = sum(z)/natoms
    
    
    x = x - centrd(1)
    y = y - centrd(2)
    z = z - centrd(3)
    
  end subroutine geometric_center


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  Removes the total momentum of a velocity vector
  ! 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine remove_total_momentum(natoms,vec,mass)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3*natoms), intent(inout), target :: vec
    double precision, dimension(3*natoms), intent(in),    target :: mass 
    double precision, dimension(:), pointer :: vx, vy, vz

    double precision :: pxtot, pytot, pztot, mcm
 
    vx => vec(1:3*natoms:3)
    vy => vec(2:3*natoms:3)
    vz => vec(3:3*natoms:3)
   
    mcm = sum(mass(1:natoms))

    pxtot = dot_product(mass(1:natoms),vx)
    pytot = dot_product(mass(1:natoms),vy)
    pztot = dot_product(mass(1:natoms),vz)

!    pxtot = dot_product(mass(1:natoms),vx)/natoms
!    pytot = dot_product(mass(1:natoms),vy)/natoms
!    pztot = dot_product(mass(1:natoms),vz)/natoms

    vx = vx - pxtot/mcm
    vy = vy - pytot/mcm
    vz = vz - pztot/mcm
  
!    vx = vx - pxtot/mass(1:natoms)
!    vy = vy - pytot/mass(1:natoms)
!    vz = vz - pztot/mass(1:natoms)   

 
    return
  end subroutine remove_total_momentum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  Removes the center of mass of a position vector
  ! 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine center_of_mass(natoms,vec,mass)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3*natoms), intent(inout), target :: vec
    double precision, dimension(3*natoms), intent(in),    target :: mass 
    double precision, dimension(:), pointer :: x, y, z

    double precision :: cm1,cm2,cm3,total_mass
    
    x => vec(1:3*natoms:3)
    y => vec(2:3*natoms:3)
    z => vec(3:3*natoms:3)

    total_mass = sum(mass(1:natoms))
    cm1        = dot_product(mass(1:natoms),x)
    cm2        = dot_product(mass(1:natoms),y)
    cm3        = dot_product(mass(1:natoms),z)

    cm1 = cm1 / total_mass
    cm2 = cm2 / total_mass
    cm3 = cm3 / total_mass

    x = x - cm1
    y = y - cm2
    z = z - cm3

    return
  end subroutine center_of_mass
    
  subroutine remove_rotation(natoms,pos,vel,mass) 
    implicit none
    
    integer, intent(in) :: natoms
    double precision, dimension(3*natoms), intent(in), target    :: pos, mass
    double precision, dimension(3*natoms), intent(inout), target :: vel
    
    integer :: i
    double precision, dimension(:), pointer :: x, y, z, vx, vy, vz
    double precision :: total_mass, trace
    double precision, dimension(3)   :: cm, l, temp, v, l1, v1, o
    double precision, dimension(3,3) :: inertia, i1
 
    x => pos(1:3*natoms:3)
    y => pos(2:3*natoms:3)
    z => pos(3:3*natoms:3)

    vx => vel(1:3*natoms:3)
    vy => vel(2:3*natoms:3)
    vz => vel(3:3*natoms:3)
   
!    cm(1) = dot_product(mass(1:natoms),x)
!    cm(2) = dot_product(mass(1:natoms),y)
!    cm(3) = dot_product(mass(1:natoms),z)
!    total_mass = sum(mass(1:natoms))
!    cm = cm/total_mass
 
    inertia = 0.0d0
    l = 0.0d0
    do i=1, natoms
       temp(1) = x(i)  ! - cm(1)
       temp(2) = y(i)  ! - cm(2)
       temp(3) = z(i)  ! - cm(3)
       v(1) = vx(i)
       v(2) = vy(i)
       v(3) = vz(i)
       call cross(l1, temp, v)
       l1 = l1 * mass(i)    ! Vectorial operation
       l = l + l1           ! Vectorial operation

       call tensor_product(i1, temp, temp, mass(i));
       inertia = inertia - i1    ! Vectorial operation
    end do
    
    trace = inertia(1,1) + inertia(2,2) + inertia(3,3)
    inertia(1,1) = inertia(1,1) - trace
    inertia(2,2) = inertia(2,2) - trace
    inertia(3,3) = inertia(3,3) - trace
    
    call solve_3x3(inertia, l, o);
    
    do i=1, natoms
       temp(1) = x(i)  ! - cm(1)
       temp(2) = y(i)  ! - cm(2)
       temp(3) = z(i)  ! - cm(3)
       call cross(v1, o, temp)
       vx(i) = vx(i) - v1(1)
       vy(i) = vy(i) - v1(2)
       vz(i) = vz(i) - v1(3)
    end do
    
  end subroutine remove_rotation
  
   subroutine solve_3x3(U,V,W)
    implicit none
    double precision, dimension(3,3), intent(in) :: U
    double precision, dimension(3), intent(in)   :: V
    double precision, dimension(3), intent(out)  :: W

    double precision :: a, b, c, d, e, f, p, q, o
    double precision :: af_de, aq_eo, ab_dd, ac_ee
    double precision :: x, y, z

    a = U(1,1)
    b = U(2,2)
    c = U(3,3)
    d = U(1,2)
    e = U(1,3)
    f = U(2,3)
    o = V(1)
    p = V(2)
    q = V(3)

    af_de = a*f-d*e
    aq_eo = a*q-e*o
    ab_dd = a*b-d*d
    ac_ee = a*c-e*e

    z = (af_de*(a*p-d*o)-ab_dd*aq_eo) / (af_de*af_de-ab_dd*ac_ee)
    y = (aq_eo - z*ac_ee)/af_de
    x = (o - d*y - e*z)/a

    W(1) = x
    W(2) = y
    W(3) = z

    return
  end subroutine solve_3x3

  ! This subroutine computes the elements of matrix mtr_c, here mtr_c is a 3*3 matrix. 
  ! mtr_c(1,1) is the sum of xa(i)*xb(i), mtr_c(1,2) is the sum of ya(i)*xb(i), mtr_c(1,3) 
  ! is the sum of za(i)*xb(i) and so on 
  subroutine matrix_c(natoms,posa,posb,mtr_c)
    implicit none
    
    integer, intent(in) :: natoms
    double precision, dimension(3*natoms), intent(in), target :: posa, posb
    
    double precision, dimension(3,3), intent(out) :: mtr_c
    double precision, dimension(:), pointer :: xa, ya, za, xb, yb, zb
    integer :: i
    
    ! We first set-up pointers for the x, y, z components for posref and pos
    xa => posa(1:3*NATOMS:3)
    ya => posa(2:3*NATOMS:3)
    za => posa(3:3*NATOMS:3)
    
    xb => posb(1:3*NATOMS:3)
    yb => posb(2:3*NATOMS:3)
    zb => posb(3:3*NATOMS:3)
    
    mtr_c = 0.0d0
    
    ! we get the matrix C(mtr_c), the matrix element mtr_c(k,j)=sum(i)((x_ji*y_ki))   
    
    do i=1, NATOMS
       mtr_c(1,1) = mtr_c(1,1)+xa(i)*xb(i)
       mtr_c(1,2) = mtr_c(1,2)+ya(i)*xb(i)
       mtr_c(1,3) = mtr_c(1,3)+za(i)*xb(i)
       
       mtr_c(2,1) = mtr_c(2,1)+xa(i)*yb(i)
       mtr_c(2,2) = mtr_c(2,2)+ya(i)*yb(i)
       mtr_c(2,3) = mtr_c(2,3)+za(i)*yb(i)  
       
       mtr_c(3,1) = mtr_c(3,1)+xa(i)*zb(i)
       mtr_c(3,2) = mtr_c(3,2)+ya(i)*zb(i)
       mtr_c(3,3) = mtr_c(3,3)+za(i)*zb(i)
    end do
    
  end subroutine matrix_c

  subroutine get_momentum_rotation(nstep,natoms,pos,vel,mass)
    implicit none
    integer, intent(in) :: nstep,natoms
    double precision, dimension(3*natoms), intent(in), target :: pos, vel, mass
    double precision, dimension(:), pointer :: x, y, z, vx, vy, vz, ms

    double precision, dimension(3) :: com, momentum
    double precision :: total_mass
    double precision :: trace
    double precision, dimension(3) :: temp, L, V, L1, O
    double precision, dimension(3,3) :: inertia, i1

    integer :: ia

    x => pos(1:3*natoms:3); vx => vel(1:3*natoms:3)
    y => pos(2:3*natoms:3); vy => vel(2:3*natoms:3)
    z => pos(3:3*natoms:3); vz => vel(3:3*natoms:3)

    ms => mass(1:natoms); total_mass = sum(ms)

    com = (/ dot_product(ms,x), dot_product(ms,y), dot_product(ms,z) /) /total_mass
    momentum = (/ dot_product(ms,vx), dot_product(ms,vy), dot_product(ms,vz) /) /natoms
    inertia = 0.0d0
    L = 0.0d0
    do ia = 1, natoms
      temp = (/ x(ia), y(ia), z(ia) /) -com
      V = (/ vx(ia), vy(ia), vz(ia) /)
      call cross(L1, temp, V)
      L = L + L1
      call tensor_product(i1, temp, temp, ms(ia))
      inertia = inertia - i1
    end do  
    trace = inertia(1,1) + inertia(2,2) + inertia(3,3)
    inertia(1,1) = inertia(1,1) - trace
    inertia(2,2) = inertia(2,2) - trace
    inertia(3,3) = inertia(3,3) - trace
    
    call solve_3x3(inertia, L, o)

    write(99, "(i10,9f12.6)") nstep,com,momentum,O
    
  end subroutine get_momentum_rotation
  
end module geometric_corrections
