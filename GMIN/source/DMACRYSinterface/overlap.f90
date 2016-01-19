module dma_overlap

double precision cutoff
parameter(cutoff = 1.5)

contains
! calcualte the minimal connection between vectors using triclinic pbc
subroutine mindist(r_i,r_j,r_ij,mlattice)
    use vec3
    use cell_reduction
    implicit none
    double precision, intent(in) :: r_i(3), r_j(3), mlattice(3,3)
    double precision, intent(out) :: r_ij(3)
    double precision x1(3),x2(3)
    double precision a(3),b(3),c(3)
    double precision r_tp(3), r_dp(3), r_sp(3)
    a = mlattice(:,1)
    b = mlattice(:,2)
    c = mlattice(:,3)
    r_tp = r_j - r_i;
    r_dp = r_tp - c*dnint(r_tp(3)/c(3))
    r_sp = r_dp - b*dnint(r_dp(2)/b(2))
    r_ij = r_sp - a*dnint(r_sp(1)/a(1))
end subroutine

function overlap_potential(natoms, x, g, mlattice) result(e)
    use vec3
    implicit none
    integer, intent(in) :: natoms
    double precision e
    double precision, intent(in) :: x(3*natoms), mlattice(3,3)
    double precision, intent(out) :: g(3*natoms)
    double precision k

    double precision rij(3), minr,r
    integer i,j
    double precision a(3),b(3),c(3), mi(3,3)
    a = mlattice(:,1)
    b = mlattice(:,2)
    c = mlattice(:,3)
    k=1000

    e=0
    minr=1000.
    call invert3x3(mlattice, mi)
    do i=1,natoms
        do j=i+1,natoms
            call mindist(x(3*i-2:3*i), x(3*j-2:3*j), rij, mlattice)
            r = vec_len(rij)
            minr = min(r,minr)
            if(r<cutoff) then
                e=e+k*(cutoff-r)**2
                g(3*i-2:3*i)=g(3*i-2:3*i)+k*(cutoff-r)*rij/r
                g(3*j-2:3*j)=g(3*j-2:3*j)-k*(cutoff-r)*rij/r
            endif
        end do
    end do
end function


function overlap_potential_wrapper(x,grad, xrigid,mlattice, gradt) result(ereal)
    use commons, only: natoms
    use genrigid, only: transformgrad, degfreedoms
    implicit none
    double precision, intent(in)  :: x(1:3*NATOMS)
    double precision, intent(in)  :: xrigid(degfreedoms)
    double precision, intent(in)  :: mlattice(3,3)
    double precision ereal
    double precision, intent(out) :: grad(1:3*NATOMS)
    logical gradt

    double precision :: grad_atom(1:3*NATOMS)

!    print *, "mindist", minimum_distance(x)

    grad(:)=0.0d0
    grad_atom(:)=0.0d0
    ereal = overlap_potential(natoms-2, x(1:3*natoms-6), grad_atom(1:3*natoms-6), mlattice)
    grad_atom(3*natoms-5:)=0.0d0
    if(gradt) call transformgrad(grad_atom, xrigid, grad)

end function

function minimum_distance(x) result(minr)
    use commons, only: natoms
    use vec3
    use genrigid, only: get_lattice_matrix

    implicit none
    double precision x(3*natoms)
    double precision mlattice(3,3)
    double precision rij(3), minr,r
    integer i,j
    double precision a(3),b(3),c(3), mi(3,3)

    call get_lattice_matrix(x(3*natoms-5:3*natoms), mlattice)


    a = mlattice(:,1)
    b = mlattice(:,2)
    c = mlattice(:,3)

    minr=1000.
    call invert3x3(mlattice, mi)
    do i=1,natoms
        do j=i+1,natoms
            call mindist(x(3*i), x(3*j), rij, mlattice)
            r = vec_len(rij)
            minr = min(r,minr)
        end do
    end do

end function

end module
