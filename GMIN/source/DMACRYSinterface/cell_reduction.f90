! Minimum cell algorithm
!
! This module is an implementation of the numerical stablilized mimimum cell algorithm
! as given by Grosse-Kunstleve, Sauter & Adams
! "Numerically stable algorithms for the computation of reduced unit cells"
! Acta Crystallographica Section A Foundations of Crystallography 60, 1, 2003
!
! This implementation was adapted from the c++ reference implentation given in the
! cctbx framework (cctbx.sourceforge.net), rewritten and adjusted for Fortran
!
! vr274>
! TODO: the current implementation is not thread safe! One has to
!        use an TYPE for the cell setup which is passed in every function
!
module cell_reduction
  ! these are the public functions to access
  ! funcionality of cell reduction algorithm
  public cell_set_grubernotation, cell_set_lattice, cell_set_angles
  public cell_minimum_reduction, cell_get_lattice, cell_get_transformation
  public cell_has_changed, entier, cell_get_lattice_for_cif

private
  ! the lattice matrix in gruber notation
  double precision a_,b_,c_,d_,e_,f_
  ! inverse transformation to reduced cell
  ! will be constructed on minimum reduction
  double precision r_inv_(3,3)

  ! some control variables from the original algorithm
  double precision :: multiplier_significant_change_test_
  integer n_no_significant_change_, min_n_no_significant_change_
  integer n_iterations_, iteration_limit_

  double precision last_abc_significant_change_test_(3)

  parameter(multiplier_significant_change_test_ = 16d0)
  parameter(iteration_limit_=100)
  parameter(min_n_no_significant_change=2)

  logical verbose
  parameter(verbose = .false.)
  logical do_acute_cell
  parameter(do_acute_cell = .false.)

contains
  ! set cell dimensions in gruber notation
  ! (a*a, b*b, c*c, 2*b*c, 2*a*c, 2*a*b)
  subroutine cell_set_grubernotation(cell)
    implicit none
    double precision cell(6)
    a_ = cell(1)
    b_ = cell(2)
    c_ = cell(3)
    d_ = cell(4)
    e_ = cell(5)
    f_ = cell(6)
  end subroutine

  ! set cell dimensions in lattice vectors
  ! lattice(:,1), lattice(:,2), lattice(:,3)
  subroutine cell_set_lattice(lattice)
    implicit none
    double precision lattice(3,3)
    double precision la(3), lb(3), lc(3)

    la = lattice(:,1)
    lb = lattice(:,2)
    lc = lattice(:,3)

    a_ = dot_product(la,la)
    b_ = dot_product(lb,lb)
    c_ = dot_product(lc,lc)
    d_ = 2d0*dot_product(lb,lc)
    e_ = 2d0*dot_product(la,lc)
    f_ = 2d0*dot_product(la,lb)
  end subroutine

  ! set cell dimenstion with angle notation
  ! (a, b, c, alpha, beta, gamma)
  subroutine cell_set_angles(m)
    implicit none
    double precision m(6)
    double precision pi
    parameter (pi=3.141592654d0)

    a_ = m(1)*m(1)
    b_ = m(2)*m(2)
    c_ = m(3)*m(3)
    d_ = 2d0*cos(m(4)/180d0*pi)*m(2)*m(3)
    e_ = 2d0*cos(m(5)/180d0*pi)*m(1)*m(3)
    f_ = 2d0*cos(m(6)/180d0*pi)*m(1)*m(2)
  end subroutine

  ! execute the minimum reduction
  subroutine cell_minimum_reduction()
    implicit none
    n_iterations_ = 0
    n_no_significant_change_ = 0
    r_inv_=0d0
    r_inv_(1,1)=1d0
    r_inv_(2,2)=1d0
    r_inv_(3,3)=1d0

    last_abc_significant_change_test_ = (/-a_,-b_,-c_/)

    do while (.not.step())
    enddo
  end subroutine

  function cell_has_changed()
    use vec3
    implicit none
    logical cell_has_changed
    double precision id(3,3), n
    call identity3x3(id)
    cell_has_changed = .true.
    id = r_inv_ - id
    n =  dot_product(id(:,1),id(:,1)) + dot_product(id(:,2),id(:,2)) + dot_product(id(:,3),id(:,3))
    if(n < 1e-8) cell_has_changed = .false.
  end function

function cell_get_niter()
	integer cell_get_niter
    cell_get_niter = n_iterations_
end function

  ! get the lattice dimension as lattice vectors
  ! lattice(:,1), lattice(:,2), lattice(:,3)
  subroutine cell_get_lattice(cell)
    implicit none
    double precision cell(3,3)
    double precision a, b, c, alpha, beta, gam
    a = sqrt(a_)
    b = sqrt(b_)
    c = sqrt(c_)
    alpha = acos(0.5d0 * d_ / b / c)
    beta = acos(0.5d0 * e_ / a / c)
    gam = acos(0.5d0 * f_ / a / b)
    cell=0
    cell(3,3)=c
    cell(2,2)=b*sin(alpha)
    cell(3,2)=b*cos(alpha)
    cell(1,1)=a/sin(alpha)*sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gam)**2  + 2*cos(alpha)*cos(beta)*cos(gam))
    cell(2,1)=a*(cos(gam) - cos(beta)*cos(alpha))/sin(alpha)
    cell(3,1)=a*cos(beta)
  end subroutine

  subroutine cell_get_lattice_for_cif(cell)
    double precision cell(3,3)
    double precision a, b, c, alpha, beta, gamma
    double precision alphaa, betaa, gammaa, v, aa, bb, cc

    a = sqrt(a_)
    b = sqrt(b_)
    c = sqrt(c_)
    alpha = acos(0.5d0 * d_ / b / c)
    beta = acos(0.5d0 * e_ / a / c)
    gamma = acos(0.5d0 * f_ / a / b)

    v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)+2*cos(alpha)*cos(beta)*cos(gamma))

    aa=sin(alpha)/a/v
    bb=sin(beta )/b/v
    cc=sin(gamma)/c/v

    alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) )
    betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) )
    gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) )

    cell(1,1)=a
    cell(1,2)=b*cos(gamma)
    cell(1,3)=c*cos(beta)

    cell(2,1)=0
    cell(2,2)=b*sin(gamma)
    cell(2,3)=-c*sin(beta)*cos(alphaa)

    cell(3,1)=0
    cell(3,2)=0
    cell(3,3)=1/cc
  end subroutine

  ! get the transformation matrix from cell to reduced cell
  function cell_get_transformation() result(r)
    use vec3
    double precision r(3,3)
    double precision xx(3,3)
    call invert3x3(r_inv_,r)
    call invert3x3(r_inv_,xx)
!    print *,"The transformation is"
!    print *,r_inv_
  end function

  ! greatest integer which is not greater than x
  function entier(x)
    implicit none
    double precision, intent(in) :: x
    double precision entier
    entier = dint(x)
    if (x-entier < 0d0) entier = entier - 1d0
    if (x-entier >= 1d0) entier = entier + 1d0 ! work around rounding errors
  end function


  subroutine def_test(n_zero, n_positive)
    implicit none
    integer, intent(out) ::  n_zero, n_positive
    n_zero = 0
    n_positive = 0

    if (0d0 < d_) then
      n_positive = n_positive + 1
    else if (.not.(d_ < 0d0)) then
      n_zero = n_zero + 1
    endif

    if (0d0 < e_) then
      n_positive = n_positive + 1
    else if (.not.(e_ < 0d0)) then
      n_zero = n_zero + 1
    endif

    if (0d0 < f_) then
      n_positive = n_positive + 1
    else if (.not.(f_ < 0d0)) then
      n_zero = n_zero + 1
    endif
  end subroutine

  function def_gt_0()
    implicit none
    integer n_zero, n_positive
    logical def_gt_0
    call def_test(n_zero, n_positive);
    def_gt_0 =  (n_positive==3).or.((n_zero == 0).and.(n_positive==1))
  end function


  function significant_change_testd(new_value, i) result(out)
    implicit none
    logical out
    integer i
    double precision :: new_value
    double precision ::  m_new, diff, m_new_plus_diff
    double precision ::  m_new_plus_diff_minus_m_new

    m_new = multiplier_significant_change_test_ * new_value;
    diff = new_value - last_abc_significant_change_test_(i);

    ! TODO: what is this
    ! m_new_plus_diff = spoil_optimization(m_new + diff);
    m_new_plus_diff = m_new + diff;

    m_new_plus_diff_minus_m_new = m_new_plus_diff - m_new;
    out = abs(m_new_plus_diff_minus_m_new)  > 1d-10
  end function

  function significant_change_test() result(out)
    implicit none
    logical out

    if ( significant_change_testd(a_, 1).or.significant_change_testd(b_, 2) &
    .or.significant_change_testd(c_, 3)) then
      n_no_significant_change_ = 0
    else
      n_no_significant_change_ = n_no_significant_change_ +1
      if (n_no_significant_change_ == min_n_no_significant_change_) then
        out = .false.
        return
      endif
    endif
    last_abc_significant_change_test_ = (/a_,b_,c_/)
    out = .true.
    return
  end function

  function step() result(done)
    implicit none
    logical done

    !    print *,"Do step"
    !    print *,a_,b_,c_
    !    print *,d_,e_,f_

    ! N1
    if (b_ < a_) then
      call n1_action
    endif
    ! N2
    if (c_ < b_) then
      call n2_action
      done=.false.
      return;
    endif
    
    ! N3
    if (def_gt_0().or.do_acute_cell) then
      call n3_true_action
    else
      call n3_false_action
      if (.not.significant_change_test()) then
        done=.true.
        return
      endif
    endif
    
    if (b2_action()) then
      done=.false.;
      return
    endif
    if (b3_action()) then
      done=.false.
      return
    endif
    if (b4_action()) then
      done=.false.
      return
    endif
    if (b5_action()) then
      done=.false.
      return
    endif
    done=.true.
  end function


  subroutine cb_update(mt)
    use vec3
    implicit none
    double precision mt(9)
    double precision m(3,3)
    if (n_iterations_ == iteration_limit_) then
      print *, "ERROR: minimum reduction failed, iteration limit exceeded";
      stop
    endif

    m(1,1:3) = mt(1:3)
    m(2,1:3) = mt(4:6)
    m(3,1:3) = mt(7:9)
    if(abs(det3x3(m)) < 0.1d0) then
        print *,"ERROR in cell_reduction::cb_update, invalid transformation"
        print *,m
        print *,mt
        stop
    endif
    r_inv_ = matmul(r_inv_, m)
    n_iterations_ = n_iterations_ + 1
  end subroutine

  subroutine swap(a, b)
    implicit none
    double precision a,b,tmp
    tmp=a
    a=b
    b=tmp
  end subroutine

  ! swap a and b
  subroutine n1_action
    implicit none
    if(verbose) print *, "N1 action"
    call cb_update((/0d0,-1d0,0d0, -1d0,0d0,0d0, 0d0,0d0,-1d0/))
    call swap(a_, b_)
    call swap(d_, e_)
  end subroutine

  ! swap b and c
  subroutine n2_action
    implicit none
    if(verbose) print *, "N2 action"
    call cb_update((/-1d0,0d0,0d0, 0d0,0d0,-1d0, 0d0,-1d0,0d0/))
    call swap(b_, c_)
    call swap(e_, f_)
  end subroutine

  subroutine n3_true_action
    use vec3
    implicit none
    double precision m(9)

    m=0
    m(1) = 1d0
    m(5) = 1d0
    m(9) = 1d0

    if (d_ < 0d0) m(1) = -1d0;
    if (e_ < 0d0) m(5) = -1d0;
    if (f_ < 0d0) m(9) = -1d0;
    if(verbose) print *, "N3 true action"
    call cb_update(m)
    d_ = abs(d_)
    e_ = abs(e_)
    f_ = abs(f_)
  end subroutine

  subroutine n3_false_action
    implicit none
    double precision m(9)
    integer z
    m=0
    m(1) = 1d0
    m(5) = 1d0
    m(9) = 1d0
    z = -1
    if (0d0 < d_) then
      m(1) = -1d0
    else if (.not.(d_ < 0d0)) then
      z = 1
    endif
    if (0d0 < e_) then
      m(5) = -1d0
    else if (.not.(e_ < 0d0)) then
      z = 5
    endif
    if (0d0 < f_) then
      m(9) = -1d0
    else if (.not.(f_ < 0d0)) then
      z = 9
    endif
    if (m(1)*m(5)*m(9) < 0d0) then
      if(z == -1) then
        print *,"error in n3_false_action"
        stop
      endif
      m(z) = -1d0
    endif
    if(verbose) print *, "N3 false"

    call cb_update(m)
    d_ = -abs(d_)
    e_ = -abs(e_)
    f_ = -abs(f_)

  end subroutine

  function b2_action() result(out)
    implicit none
    logical out
    double precision j


    out=.false.
    if (.not.(b_ < abs(d_))) return
    j = entier((d_+b_)/(2d0*b_))
    if (int(j) == 0) return
    if(verbose) print *, "B2 action"
    call cb_update((/1d0,0d0,0d0,0d0,1d0,-j,0d0,0d0,1d0/))
    c_ = c_ + j*j*b_ - j*d_;
    d_ = d_ - 2d0*j*b_;
    e_ = e_ - j*f_;
    if (.not.(0d0 < c_)) then
      print *,"ERROR: degenerate unit cell"
      stop
    endif
    out=.true.
  end function

  function b3_action() result(out)
    implicit none
    logical out
    double precision j
    out = .false.
    if (.not.(a_ < abs(e_))) return
    j = entier((e_+a_)/(2d0*a_));
    if (int(j) == 0) return
    if(verbose) print *, "B3 action"
    call cb_update((/1d0,0d0,-j,0d0,1d0,0d0,0d0,0d0,1d0/))
    c_ = c_ + j*j*a_ - j*e_;
    d_ = d_ - j*f_;
    e_ = e_ - 2d0*j*a_;
    if (.not.(0d0 < c_)) then
      print *,"ERROR: degenerate unit cell",c_
      stop
    endif
    out=.true.
    if(verbose) print *, "N3 false"

  end function

  function b4_action() result(out)
    implicit none

    logical out
    double precision j

    out = .false.
    if (.not.(a_ < abs(f_))) return
    j = entier((f_+a_)/(2*a_));
    if (int(j) == 0) return
    if(verbose) print *, "B4 action"
    call cb_update((/1d0,-j,0d0,0d0,1d0,0d0,0d0,0d0,1d0/))
    b_ = b_ + j*j*a_ - j*f_;
    d_ = d_ - j*e_;
    f_ = f_ - 2d0*j*a_;
    if (.not.(0d0 < b_)) then
      print *,"ERROR: degenerate unit cell"
      stop
    endif
    out = .true.
  end function

  function b5_action() result(out)
    implicit none
    double precision j, de, fab

    logical out
    out=.false.

    de = d_ + e_;
    fab = f_ + a_ + b_;
    if (.not.(de+fab < 0d0)) return
    j = entier((de+fab)/(2*fab));
    if (int(j) == 0) return
    if(verbose) print *, "B5 action"
    call cb_update((/1d0,0d0,-j,0d0,1d0,-j,0d0,0d0,1d0/))
    c_ = c_ + j*j*fab-j*de;
    d_ = d_ - j*(2d0*b_+f_);
    e_ = e_ - j*(2d0*a_+f_);
    if (.not.(0d0 < c_)) then
      print *,"ERROR: degenerate unit cell"
      stop
    endif
    out=.true.
  end function

end module cell_reduction
