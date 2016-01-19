SUBROUTINE RMDRVT(P, RM)

    IMPLICIT NONE
    DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)

    !     DEk = derivate of E with respect to the kth component of P
    DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
    LOGICAL          :: GTEST

    !     Set the values of the idenity matrix I3
    I3(:,:) = 0.D0
    I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

    !     Calculate the value of THETA2 as the square modulus of P
    THETA2  = DOT_PRODUCT(P,P)

    IF (THETA2 < 1.0D-12) THEN
        !        Execute if the angle of rotation is zero
        !        In this case the rotation matrix is the identity matrix
        RM(:,:) = I3(:,:)

        ! vr274> first order corrections to rotation matrix
        RM(1,2) = -P(3)
        RM(2,1) = P(3)
        RM(1,3) = P(2)
        RM(3,1) = -P(2)
        RM(2,3) = -P(1)
        RM(3,2) = P(1)

    ELSE
        !       Execute for the general case, where THETA dos not equal zero
        !       Find values of THETA, CT, ST and THETA3
        THETA   = SQRT(THETA2)
        CT      = COS(THETA)
        ST      = SIN(THETA)
        THETA3  = 1.D0/(THETA2*THETA)

        !        Set THETA to 1/THETA purely for convenience
        THETA   = 1.D0/THETA

        !        Normalise P and construct the skew-symmetric matrix E
        !        ESQ is calculated as the square of E
        PN(:)   = THETA*P(:)
        E(:,:)  = 0.D0
        E(1,2)  = -PN(3)
        E(1,3)  =  PN(2)
        E(2,3)  = -PN(1)
        E(2,1)  = -E(1,2)
        E(3,1)  = -E(1,3)
        E(3,2)  = -E(2,3)
        ESQ     = MATMUL(E,E)

        !        RM is calculated from Rodrigues' rotation formula (equation (1)
        !        in the paper)
        RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)
    ENDIF
END SUBROUTINE RMDRVT


function matrix_equals(m1,m2) result(matches)
    implicit none
    logical matches
    integer i,j
    double precision m1(3,3), m2(3,3)
    matches=.false.
    do i=1,3
        do j=1,3
            if(abs(m1(i,j)-m2(i,j)) > 1e-8) return
        enddo
    enddo
    matches=.true.
end function

subroutine test_conversions(pin)
    use vec3
    use rotations
    implicit none
    double precision pi
    parameter (pi=3.141592654d0)
    logical matrix_equals

    double precision p(3), m(3,3), p2(3), q(4),pin(3)
    double precision x(3,3), xref(3,3)
    double precision mref(3,3)
    integer i

    p = pin
    call rmdrvt(p, mref)

    ! angle axis to quaternion to matrix
    m = rot_q2mx(rot_aa2q(p))
    if(.not.matrix_equals(m,mref)) then
        print *,"Conversion angle axis to quaternion to matrix failed"
        print *,pin
        stop
    endif

    ! quaternion to angle axis
    q = rot_aa2q(p)
    p2 = rot_q2aa(q)
    call rmdrvt(p2, m)
    if(.not.matrix_equals(m,mref)) then
        print *,"Conversion quaternion to angle axis failed"
        print *,pin,vec_len(pin)
        stop
    endif

    ! matrix to quaternion
    m = rot_q2mx(rot_mx2q(mref))
    if(.not.matrix_equals(m,mref)) then
        print *,"Conversion matrix to quaternion failed"
        print *,pin
        stop
    endif
  
    ! matrix to angle axis
    p =rot_mx2aa(mref)
    call RMDRVT(p,m)
    if(.not.matrix_equals(m,mref)) then
        print *,"Conversion matrix to angle axis failed"
        print *,pin
        print *,p
        stop
    endif
  
    ! now testing calculationg orientation
    do i=1,3
        xref(:,i) = vec_random()
        x(:,i) = matmul(mref,xref(:,i))
    end do

    p = rot_get_orientation_aa(x, xref)
    call rmdrvt(p, m)
    if(.not.matrix_equals(m,mref)) then
        print *,"Calculating orientation from 3 points failed"
        print *,pin
        print *,p
        print *,"mref"
        print *,mref
        stop
    endif

    return
end subroutine

subroutine test_vector_proximity(p)
    use rotations
    implicit none
    double precision, intent(in) :: p(3)
    double precision pn(3)
    integer i
    call test_conversions(p)
    do i=1,10000
        pn = p
        call rot_takestep_aa(pn, 1d-5)
        call test_conversions(pn)
    end do
end subroutine

subroutine test_vector_angle(p, a)
    implicit none
    double precision, intent(in) :: p(3), a
    call test_vector_proximity(a*p)
end subroutine

subroutine test_multiplication
    use rotations
    implicit none
    double precision p1(3), p2(3), m1(3,3), m2(3,3), mref(3,3), p(3)
    integer i
    logical matrix_equals

    do i=1,100000
        p1 = rot_random_aa()
        p2 = rot_random_aa()
        call RMDRVT(p1, m1)
        call RMDRVT(p2, m2)
        mref = matmul(m2,m1)

        p = rot_rotate_aa(p1, p2)
        call RMDRVT(p, m1)

        if(.not.matrix_equals(m1,mref)) then
            print *,"Combining 2 angle axis vectors failed"
            print *,p1
            print *,p2
            stop
        endif
    end do

end subroutine

program test_rotations
    use vec3
    use rotations
    implicit none
    double precision pi, dprand
    integer i
    parameter (pi=3.141592654d0)

    call test_vector_angle((/1d0,0d0,0d0/), pi/4.0d0)

    do i=0,10
        call test_vector_angle((/1d0,0d0,0d0/), dble(i)*pi/4.0d0)
        call test_vector_angle((/0d0,1d0,0d0/), dble(i)*pi/4.0d0)
        call test_vector_angle((/1d0,0d0,1d0/), dble(i)*pi/4.0d0)
        call test_vector_angle(vec_random(), dble(i)*pi/4.0d0)
    end do
    do i=1,100000
        call test_conversions(rot_random_aa())
        call test_conversions((/10d0*dprand() - 5d0,10d0*dprand() - 5d0,10d0*dprand() - 5d0/))
    end do

    call test_multiplication

    print *,"All tests successful"
end program

