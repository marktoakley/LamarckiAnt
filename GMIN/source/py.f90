module py_module
    ! swo24> equation references are to:
    !  CW = Chakrabarti & Wales, PCCP 11 1970-1976 (2009)
    !  PY = Paramonov & Yaliraki, J Chem Phys 123 194111 (2005)
    ! It would be wise to read and understand these papers before trying to work through this code!

    ! General naming conventions:
    !  x_p2 = x to power 2 = x ** 2
    !  x_pm2 = x to power minus 2 = x ** (-2)
    !  x_mod = modulus of vector x = |x|
    !  dx_dy = derivative of x with respect to y
    !  mat_inv = inverse of matrix 'mat'

    implicit none

    type py_ellipsoid
        type(py_site), pointer :: site  ! pointer to parent site
        double precision :: r(3)    ! site position (can vary when using PBCs)
        double precision :: semiaxes(3), long_semiaxis, short_semiaxis  ! semiaxes of ellipsoid; longest semiaxis
        ! shape matrix in body frame (i.e., diagonal), in molecule frame, in lab frame, and its inverse (CW eq 23)
        double precision :: shapemat_bodyframe(3, 3), shapemat_molframe(3, 3), shapemat(3, 3), shapemat_inv(3, 3)
        double precision :: shapemat_det, shapemat_inv_det  ! determinant and inverse det of shape matrix
        double precision :: shapemat_inv_cofactors(3, 3)    ! cofactor matrix of inverse shape matrix
        double precision :: mat_dot_r(3)    ! shape matrix dot position (PY eq 5, second term)
        double precision :: dF_dr(3), dF_dp(3), dr_dp(3)    ! (CW eq 28); (CW eq 32); (CW eq 12 if at molecule center)
        double precision :: dmat_dp(3, 3, 3)  ! shape matrix with respect to angle axis (DW eq 33)
        double precision :: strength    ! coefficient to determiante interaction strength
    end type py_ellipsoid

    type py_site
        type(py_molecule), pointer :: molecule  ! pointer to parent molecule
        type(py_ellipsoid) :: rep, att  ! repulsive and attractive ellipsoids
        double precision :: r(3), r_molframe(3), p(3), rotmat(3, 3) ! lab and molecule frame positions; rotations
        logical :: ells_same    ! are the two ellipsoids the same shape and orentiation?
    end type py_site

    type py_molecule
        type(py_site), allocatable :: sites(:)  ! container for sites inside the molecule
        double precision :: r(3), p(3)  ! position and orientation
        double precision :: rotmat(3, 3), rotmat_deriv(3, 3, 3) ! rotation matrix and derivative
        integer :: r_index, p_index ! indices of r and p inside the big x array
    end type py_molecule

    type(py_molecule), target, allocatable :: molecules(:)  ! container for all the molecules
    integer :: total_num_sites  ! total number of sites in the whole cluster

    ! gcut = value of g (CW eq 25) above which the interaction is cut off
    double precision :: gcut, gcut_pm2, gcut_pm6, gcut_pm7, gcut_pm8, gcut_pm12, gcut_pm14

    logical :: above_cutoff ! signals between subroutines if ellipsoids are too far apart

contains
    subroutine initialize_py_potential
        ! compute values that pertain to the whole PY potential

        use commons, only: paramonovcutoff, pcutoff, pysignot, vt

        if(.not. allocated(vt)) allocate(vt(size(molecules)))   ! allocate array for pairwise energies

        if(paramonovcutoff) then
            gcut = (pcutoff + pysignot) / pysignot
            gcut_pm2 = 1.d0 / (gcut ** 2)
            gcut_pm6 = gcut_pm2 ** 3
            gcut_pm7 = gcut_pm6 / gcut
            gcut_pm8 = gcut_pm6 * gcut_pm2
            gcut_pm12 = gcut_pm6 ** 2
            gcut_pm14 = gcut_pm12 * gcut_pm2
        end if

    end subroutine initialize_py_potential

    subroutine initialize_py_ellipsoid(site, ell)
        ! set up the py ellipsoid values for the first time

        implicit none
        type(py_site), target, intent(in) :: site
        type(py_ellipsoid), intent(inout) :: ell

        double precision :: dummy(3, 3)
        integer :: k

        ell%site => site    ! set up the parent pointer
        ell%long_semiaxis = maxval(ell%semiaxes)
        ell%short_semiaxis = minval(ell%semiaxes)

        ell%shapemat_bodyframe(:, :) = 0.d0
        ell%shapemat_det = 1.d0
        do k = 1, 3
            ell%shapemat_bodyframe(k, k) = 1.d0 / (ell%semiaxes(k) ** 2)
            ell%shapemat_det = ell%shapemat_det * ell%shapemat_bodyframe(k, k) 
        end do

        ell%shapemat_inv_det = 1.d0 / ell%shapemat_det
        call rmdrvt(site%p, site%rotmat, dummy, dummy, dummy, .false.)
        ell%shapemat_molframe(:, :) = matmul(site%rotmat, matmul(ell%shapemat_bodyframe, transpose(site%rotmat)))

    end subroutine initialize_py_ellipsoid

    subroutine pairwise_py(sitei, sitej, grad_contrib, energy_contrib, gtest)
        ! compute energy and gradient due to a pair of sites

        use commons, only: pysignot, paramonovcutoff, pyepsnot

        implicit none
        type(py_site), target, intent(inout) :: sitei, sitej
        logical, intent(in) :: gtest
        double precision, intent(out) :: energy_contrib, grad_contrib(12)

        type(py_ellipsoid), pointer :: elli, ellj   ! pointers to ellipsoids inside the sites
        double precision :: g, g_pm1, g_pm2, g_pm6  ! powers of g
        double precision :: F_pm1_2, dg_dF, dg_dr(3), dg_drmod  ! F_pm1_2 = F ** (-1/2)
        double precision :: e_rep, g_rep(12), e_att, g_att(12)
        double precision :: g_pm7, g_pm12, g_pm13, g_pm14
        double precision :: dU_dF, dU_drmod

        ! repulsive part
        elli => sitei%rep
        ellj => sitej%rep

        call compute_py_values(elli, ellj, gtest, g, g_pm1, g_pm2, g_pm6, F_pm1_2, dg_dF, dg_dr, dg_drmod)
        if(paramonovcutoff .and. above_cutoff) then
            e_rep = 0.d0
            g_rep(:) = 0.d0
        else ! not over the cutoff
            g_pm12 = g_pm6 * g_pm6
            g_pm13 = g_pm12 * g_pm1

            e_rep = g_pm12
            if(paramonovcutoff) e_rep = e_rep + gcut_pm12 * (6.d0 * gcut_pm2 / g_pm2 - 7.d0)

            if(gtest) then
                dU_dF = -2.d0 * g_pm13 * dg_dF ! (CW eq 30) (need 2.d0 for 0.5d0 in dg_dF)
                dU_drmod = -2.d0 * g_pm13 * (1.d0 - F_pm1_2) / pysignot ! derivate with respect to scalar r_ij

                g_rep(1: 3) = -2.d0 * g_pm13 * dg_dr(:)   ! (CW eq 26, but note the error in the paper's sign!)
                g_rep(7: 9) = dU_dF * elli%dF_dp(:) + dU_drmod * elli%dr_dp(:) ! (CW eq 29, with multisite term)
                g_rep(10: 12) = dU_dF * ellj%dF_dp(:) + dU_drmod * ellj%dr_dp(:)

                if(paramonovcutoff) then
                    ! add splines to achieve cutoff
                    g_rep(1: 3) = g_rep(1: 3) + 2.d0 * gcut_pm14 * g * dg_dr(:)
                    g_rep(7: 9) = g_rep(7: 9) + 2.d0 * gcut_pm14 * g * (dg_dF * elli%dF_dp(:) + dg_drmod * elli%dr_dp(:))
                    g_rep(10: 12) = g_rep(10: 12) + 2.d0 * gcut_pm14 * g * (dg_dF * ellj%dF_dp(:) + dg_drmod * ellj%dr_dp(:))
                end if

                g_rep(4: 6) = -g_rep(1: 3)  ! (CW eq 34)
            end if
        end if

        ! attractive part: if ellipsoids are the same, keep the old values and ellipsoid pointers
        if(.not.(sitei%ells_same .and. sitej%ells_same)) then
            elli => sitei%att
            ellj => sitej%att
            call compute_py_values(elli, ellj, gtest, g, g_pm1, g_pm2, g_pm6, F_pm1_2, dg_dF, dg_dr, dg_drmod)
        end if

        if(paramonovcutoff .and. above_cutoff) then
            e_att = 0.d0
            g_att(:) = 0.d0
            
        else ! not over cutoff
            g_pm12 = g_pm6 * g_pm6
            g_pm13 = g_pm12 * g_pm1

            e_att = -g_pm6
            if(paramonovcutoff) e_att = e_att + gcut_pm6 * (-3.d0 * gcut_pm2 / g_pm2 + 4.d0)

            if(gtest) then
                g_pm7 = g_pm6 * g_pm1
                dU_dF = g_pm7 * dg_dF   ! (CW eq 31, but with factors of 24 and 0.5d0)
                dU_drmod = g_pm7 * (1.d0 - F_pm1_2) / pysignot

                g_att(1: 3) = g_pm7 * dg_dr(:)
                g_att(7: 9) = dU_dF * elli%dF_dp(:) + dU_drmod * elli%dr_dp(:)
                g_att(10: 12) = dU_dF * ellj%dF_dp(:) + dU_drmod * ellj%dr_dp(:)

                if(paramonovcutoff) then
                    g_att(1: 3) = g_att(1: 3) - gcut_pm8 * g * dg_dr(:)
                    g_att(7: 9) = g_att(7: 9) - gcut_pm8 * g * (dg_dF * elli%dF_dp(:) + dg_drmod * elli%dr_dp(:))
                    g_att(10: 12) = g_att(10: 12) - gcut_pm8 * g * (dg_dF * ellj%dF_dp(:) + dg_drmod * ellj%dr_dp(:))
                end if

                g_att(4: 6) = -g_att(1: 3)
            end if
        end if

        ! add on final factors
        energy_contrib = 4.d0 * pyepsnot * (sitei%rep%strength * sitej%rep%strength * e_rep &
            & + sitei%att%strength * sitej%att%strength * e_att)
        if(gtest) grad_contrib(:) = 24.d0 * pyepsnot * (sitei%rep%strength * sitej%rep%strength * g_rep(:) &
            & + sitei%att%strength * sitej%att%strength * g_att(:))

    end subroutine pairwise_py

    subroutine compute_py_values(elli, ellj, gtest, g, g_pm1, g_pm2, g_pm6, F_pm1_2, dg_dF, dg_dr, dg_drmod)
        use commons, only: pysignot, paramonovcutoff
        implicit none

        type(py_ellipsoid), target, intent(inout) :: elli, ellj
        logical, intent(in) :: gtest
        double precision, intent(out) :: g, g_pm1, g_pm2, g_pm6, F_pm1_2, dg_dF, dg_dr(3), dg_drmod

        type(py_site), pointer :: sitei, sitej
        type(py_molecule), pointer :: moli, molj
        double precision :: lambda_c, F, rij(3), rij_p2, rijmod, rij_hat(3)
        double precision :: x_c(3), xc_minus_ri(3), xc_minus_rj(3)  ! PY eq 9; last term in CW eq 28
        double precision :: xlambdac_term1(3, 3), xlambdac_term2(3) ! PY eq 5 with lambda_c from after PY eq 7
        integer :: i

        sitei => elli%site
        sitej => ellj%site
        moli => sitei%molecule
        molj => sitej%molecule

        call compute_contact(elli, ellj, lambda_c, F)

        rij(:) = elli%r(:) - ellj%r(:)
        rij_p2 = dot_product(rij, rij)
        rijmod = sqrt(rij_p2)
        rij_hat(:) = rij(:) / rijmod
        F_pm1_2 = 1.d0 / sqrt(F)
        g = (rijmod * (1.d0 - F_pm1_2) + pysignot) / pysignot    ! (CW eq 25)

        if(paramonovcutoff .and. (above_cutoff .or. g > gcut)) then
            ! if over cutoff, do nothing
            above_cutoff = .true.
        else
            g_pm1 = 1.d0 / g
            g_pm2 = g_pm1 * g_pm1
            g_pm6 = g_pm2 * g_pm2 * g_pm2

            if(gtest) then
                xlambdac_term1(:, :) = inverse(lambda_c * elli%shapemat + (1.d0 - lambda_c) * ellj%shapemat)
                xlambdac_term2(:) = lambda_c * elli%mat_dot_r + (1.d0 - lambda_c) * ellj%mat_dot_r
                x_c(:) = matmul(xlambdac_term1, xlambdac_term2)
                xc_minus_ri(:) = x_c(:) - elli%r(:)
                xc_minus_rj(:) = x_c(:) - ellj%r(:)

                do i = 1, 3
                    elli%dF_dp(i) = lambda_c * (dot_product(xc_minus_ri, matmul(elli%dmat_dp(i, :, :), xc_minus_ri)) &
                        &- 2.d0 * dot_product(xc_minus_ri, matmul(elli%shapemat, &
                        & matmul(moli%rotmat_deriv(i, :, :), sitei%r_molframe(:)))))    ! CW eq 32
                   ellj%dF_dp(i) = (1.d0 - lambda_c) * (dot_product(xc_minus_rj, matmul(ellj%dmat_dp(i, :, :), xc_minus_rj)) &
                        &- 2.d0 * dot_product(xc_minus_rj, matmul(ellj%shapemat, &
                        & matmul(molj%rotmat_deriv(i, :, :), sitej%r_molframe(:)))))    ! CW eq 35
                end do
                elli%dF_dr(:) = -2.d0 * lambda_c * matmul(elli%shapemat, xc_minus_ri)    ! PY eq C2
                ellj%dF_dr(:) = -elli%dF_dr(:)  ! PY eq C3

                dg_dF = 0.5d0 * rijmod * F_pm1_2 / (F * pysignot)
                dg_dr(:) = (1.d0 - F_pm1_2) * rij_hat(:) / pysignot + dg_dF * elli%dF_dr(:) ! CW eq 17 first term

                if(paramonovcutoff) dg_drmod = (1.d0 - F_pm1_2) / pysignot

                do i = 1, 3
                    elli%dr_dp(i) = dot_product(rij_hat(:), matmul(moli%rotmat_deriv(i, :, :), sitei%r_molframe(:)))
                    ellj%dr_dp(i) = dot_product(-rij_hat(:), matmul(molj%rotmat_deriv(i, :, :), sitej%r_molframe(:)))
                end do
            end if
        end if
    end subroutine compute_py_values

    subroutine compute_contact(elli, ellj, lambda_c, F)
        use commons, only: paramonovpbcx, paramonovpbcy, paramonovpbcz, boxlx, boxly, boxlz, &
            & paramonovcutoff, pcutoff, pycfthresh, pycoldfusion

        implicit none
        type(py_ellipsoid), intent(inout) :: elli, ellj
        double precision, intent(out) :: lambda_c, F

        integer :: image_list_offset, xyz, o1, o2, image_num
        double precision :: pbc_values(3), image_rj(8, 3), rij(3), rij_p2, rijmod
        double precision :: this_lambda_c, d_R, this_d_R, this_F

        ! make sure images are initially set properly
        above_cutoff = .false.
        elli%r(:) = elli%site%r(:)
        ellj%r(:) = ellj%site%r(:)
        rij(:) = elli%r(:) - ellj%r(:)
        ellj%mat_dot_r(:) = matmul(ellj%shapemat, ellj%r)

        if(paramonovpbcx .or. paramonovpbcy .or. paramonovpbcz) then
            ! set up the pbc matrix
            pbc_values(:) = 0.d0
            if(paramonovpbcx) pbc_values(1) = boxlx
            if(paramonovpbcy) pbc_values(2) = boxly
            if(paramonovpbcz) pbc_values(3) = boxlz

            do xyz = 1, 3
                if(pbc_values(xyz) /= 0.d0 .and. abs(rij(xyz)) > pbc_values(xyz) / 2.d0) &
                    ellj%r(xyz) = ellj%r(xyz) + sign(pbc_values(xyz), rij(xyz))
            end do

            ! redo dependent values, but don't change any orientational stuff
            rij(:) = elli%r(:) - ellj%r(:)
            ellj%mat_dot_r(:) = matmul(ellj%shapemat, ellj%r)

        end if

        rijmod = sqrt(dot_product(rij, rij))

        if(paramonovcutoff .and. rijmod - elli%long_semiaxis - ellj%long_semiaxis > pcutoff) then
            ! above cutoff
            lambda_c = 0.d0
            F = -1.d0
            above_cutoff = .true.
        else
            call polyecf(elli%shapemat_inv, ellj%shapemat_inv, elli%shapemat_inv_cofactors, &
                & ellj%shapemat_inv_cofactors, elli%shapemat_inv_det, ellj%shapemat_inv_det, & 
                & rij, lambda_c, F)
                if(F < pycfthresh) pycoldfusion = .true.
        end if
        ! if ECF falls below given value, then flag for cold fusion of the ellipsoids

    end subroutine compute_contact

    subroutine update_py_molecule(mol, r, p, gtest)
        ! changes the r and p of a molecule, then updates all the appropriate mol, site, and ell values

        use commons, only: paramonovpbcx, paramonovpbcy, paramonovpbcz, boxlx, boxly, boxlz ! PBC parameters
        implicit none
        double precision, parameter :: pi = 3.14159265358979323846
        type(py_molecule), target, intent(inout) :: mol
        type(py_site), pointer :: site
        double precision, intent(inout) :: r(3), p(3)   ! position and orientation vectors
        logical, intent(in) :: gtest    ! if gradients are required

        integer :: i
        double precision :: pmod    ! modulus of p

        ! put molecule back in PBC box
        if(paramonovpbcx) r(1) = r(1) - boxlx * nint(r(1) / boxlx)
        if(paramonovpbcy) r(2) = r(2) - boxly * nint(r(2) / boxly)
        if(paramonovpbcz) r(3) = r(3) - boxlz * nint(r(3) / boxlz)
        mol%r(:) = r(:)

        ! make sure that 0 < |p| < 2*pi
        pmod = sqrt(dot_product(p, p))
        if(pmod > 2 * pi) p(:) = p(:) / pmod * mod(pmod, 2 * pi)
        mol%p(:) = p(:)

        ! recompute the rotation matrices and derivatives
        call rmdrvt(mol%p, mol%rotmat, &
            & mol%rotmat_deriv(1, :, :), mol%rotmat_deriv(2, :, :), mol%rotmat_deriv(3, :, :), gtest)

        ! loop over sites
        do i = 1, size(mol%sites)
            site => mol%sites(i)

            site%r(:) = mol%r(:) + matmul(mol%rotmat, site%r_molframe(:))   ! reset site position
            call update_py_ellipsoid(site%rep)  ! update rep ellipsoid
            if(.not. site%ells_same) call update_py_ellipsoid(site%att) ! update att ell if it's different
        end do

    end subroutine update_py_molecule

    subroutine update_py_ellipsoid(ell)
        ! update values specific to a ellipsoid

        implicit none
        type(py_ellipsoid), intent(inout) :: ell

        type(py_site), pointer :: site
        type(py_molecule), pointer :: mol
        integer :: k
        double precision :: dummy(3, 3)

        ! identify parent site and molecule
        site => ell%site
        mol => site%molecule

        ell%r(:) = site%r(:)    ! put ellipsoid at site position
        ell%shapemat(:, :) = matmul(mol%rotmat, matmul(ell%shapemat_molframe, transpose(mol%rotmat)))
        ell%shapemat_inv = adjugate(ell%shapemat) / ell%shapemat_det    ! compute inverse
        ell%shapemat_inv_cofactors = transpose(adjugate(ell%shapemat_inv)) ! cofactor matrix is transpose of adjugate matrix
        ell%mat_dot_r(:) = matmul(ell%shapemat, ell%r)

        do k = 1, 3
            dummy(:, :) = matmul(mol%rotmat_deriv(k, :, :), ell%shapemat_molframe)
            ell%dmat_dp(k, :, :) = matmul(dummy, transpose(mol%rotmat)) + matmul(mol%rotmat, transpose(dummy))
        end do

    end subroutine update_py_ellipsoid

    function py_overlap(elli, ellj)
        ! return true if two ellipsoids overlap

        use commons, only: pyoverlapthresh, pycoldfusion
        use vec3
        type(py_ellipsoid), intent(inout) :: elli, ellj
        logical :: py_overlap

        double precision :: rij(3), rijmod, lambda_c, F

        py_overlap = .false.
        pycoldfusion = .false. 
        rij(:) = elli%r(:) - ellj%r(:)
        rijmod = vec_len(rij(:))

        ! if ellipsoids could overlap for some orientation
        if(rijmod < elli%short_semiaxis + ellj%short_semiaxis) then
            py_overlap = .true.
        elseif(rijmod < elli%long_semiaxis + ellj%long_semiaxis) then
            call compute_contact(elli, ellj, lambda_c, F)
            if(F < pyoverlapthresh) py_overlap = .true.
        end if

    end function py_overlap

    subroutine polyecf(a, b, at, bt, deta, detb, x, x1, smax)
    ! compute ECF using polynomial method and Newton's method to locate the root

        implicit none

        ! a = inverse shape matrix of one ellipsoid; at = cofactor matrix of a; deta = det of A
        ! b = other matrix
        ! x[] = separation between ellipsoid centers
        ! x1 = x_n in newton's method, and output value for lambda
        ! smax = s(lambda_c)
        double precision, intent(in) :: a(3,3), at(3,3), b(3,3), bt(3,3)
        double precision, intent(in) :: deta, detb, x(3)
        double precision, intent(out) :: x1, smax

        ! as = a*, etc
        double precision :: as, bs, cs, ds, es

        ! x_p2 = x(i) * x(j) terms
        double precision :: x_p2(3, 3)

        ! f1 = f_1, coefficient of numerator of s(lambda)
        ! g0 = g_0, coefficient of denominator of s(l)
        ! h0 = h_0, coefficient of numerator of s'(l)
        double precision :: f1, f2, f3, f4, g0, g1, g2, g3
        double precision :: h0, h1, h2, h3, h4, h5, h6

        ! x0 = x_(n-1) in newton's method, test values for lambda
        ! xx = placeholder for powers of x0 (aka lambda)
        ! p  = value of polynomial h at lambda=x0
        ! dp = value of derivative of h
        ! n = numerator of s
        ! d = denominator of s
        double precision :: x0, xx, p, dp, n, d

        ! i, j = indices
        ! converged = if x1-x0 got sufficiently small
        integer :: i, j
        logical :: converged

        ! compute x ^ 2
        forall(i = 1:3)
            forall(j = i:3)
                x_p2(i, j) = x(i) * x(j)
            end forall
        end forall

        ! compute a* = sum_ij x_i x_j at_ij (noting that at is symmetric)
        as =        x_p2(1, 1)*at(1,1) + x_p2(2, 2)*at(2,2) + x_p2(3, 3)*at(3,3) + &
            & 2.d0*(x_p2(1, 2)*at(1,2) + x_p2(1, 3)*at(1,3) + x_p2(2, 3)*at(2,3)) 

        ! compute b* = sum_ij x_i x_j bt_ij (noting that bt is symmetric)
        bs =        x_p2(1, 1)*bt(1,1) + x_p2(2, 2)*bt(2,2) + x_p2(3, 3)*bt(3,3) + &
            & 2.d0*(x_p2(1, 2)*bt(1,2) + x_p2(1, 3)*bt(1,3) + x_p2(2, 3)*bt(2,3))

        ! compute c* = sum_ij x_i x_j ct_ij (noting that ct is symmetric)
        cs =       x_p2(1, 1)*(a(2,2)*b(3,3) - a(2,3)*b(3,2) + b(2,2)*a(3,3) - b(2,3)*a(3,2)) + &
           &       x_p2(2, 2)*(a(3,3)*b(1,1) - a(3,1)*b(1,3) + b(3,3)*a(1,1) - b(3,1)*a(1,3)) + &
           &       x_p2(3, 3)*(a(1,1)*b(2,2) - a(1,2)*b(2,1) + b(1,1)*a(2,2) - b(1,2)*a(2,1)) + &
           & 2.d0*(x_p2(1, 2)*(a(2,3)*b(3,1) - a(2,1)*b(3,3) + b(2,3)*a(3,1) - b(2,1)*a(3,3)) + &
           &       x_p2(1, 3)*(a(2,1)*b(3,2) - a(2,2)*b(3,1) + b(2,1)*a(3,2) - b(2,2)*a(3,1)) + &
           &       x_p2(2, 3)*(a(3,1)*b(1,2) - a(3,2)*b(1,1) + b(3,1)*a(1,2) - b(3,2)*a(1,1)))

        ! compute d* = sum_ij a_ij bt_ij (noting symmetry)
        ds =       a(1,1)*bt(1,1) + a(2,2)*bt(2,2) + a(3,3)*bt(3,3) + &
           & 2.d0*(a(1,2)*bt(1,2) + a(1,3)*bt(1,3) + a(2,3)*bt(2,3))

        ! compute e* = sum_ij b_ij at_ij (noting symmetry)
        es =        b(1,1)*at(1,1) + b(2,2)*at(2,2) + b(3,3)*at(3,3) + &
            & 2.d0*(b(1,2)*at(1,2) + b(1,3)*at(1,3) + b(2,3)*at(2,3))

        ! compute coefficients of f = numerator of s
        f1 = as
        f2 = -3.d0*as + cs
        f3 = 3.d0*as + bs - 2.d0*cs
        f4 = - as - bs + cs

        ! compute coefficients of g = denominator of s
        g0 = deta
        g1 = -3.d0*deta + es
        g2 = 3.d0*deta + ds - 2.d0*es
        g3 = -deta + detb - ds + es

        ! compute coefficients of h = numerator of s'
        h0 = f1*g0
        h1 = 2.d0*f2*g0
        h2 = 3.d0*f3*g0 + f2*g1 - f1*g2
        h3 = 4.d0*f4*g0 + 2.d0*f3*g1 - 2.d0*f1*g3
        h4 = 3.d0*f4*g1 + f3*g2 - f2*g3
        h5 = 2.d0*f4*g2
        h6 = f4*g3

        ! begin iteration loop. guess initial position.
        x0 = 0.51d0
        converged = .false.
        do i = 1, 25
            ! compute polynomial and its derivative: first, add the constant parts. each
            ! step adds the next power onto both, and xx increases to a new power of x0
            ! each time.
            p  = h0
            dp = h1
            xx = x0
            p  = p  + h1*xx
            dp = dp + 2.d0*h2*xx
            xx = xx*x0
            p  = p  + h2*xx
            dp = dp + 3.d0*h3*xx
            xx = xx*x0
            p  = p  + h3*xx
            dp = dp + 4.d0*h4*xx
            xx = xx*x0
            p  = p  + h4*xx
            dp = dp + 5.d0*h5*xx
            xx = xx*x0
            p  = p  + h5*xx
            dp = dp + 6.d0*h6*xx
            xx = xx*x0
            p  = p  + h6*xx

            ! compute new position. if dp=0, we are at a critical point and need to 
            ! do some extra work.
            if(dp .ne. 0.d0) then
                ! take the normal newton-rapheson step
                x1 = x0 - p / dp

                ! if x1 ended up outside [0,1], then pretend that x0 was the first point
                ! computed in a bisection method. so if p>0, go right, and if p<0 go left.
                if (x1 .lt. 0.d0 .or. x1 .gt. 1.d0) then
                    if(p .gt. 0.d0) x1 = (x0 + 1.d0)/2.d0
                    if(p .lt. 0.d0) x1 =  x0/2.d0

                ! if the new position is sufficiently close to the old one, we say we found
                ! the root.
                elseif (abs(x1 - x0) .lt. 1d-8) then
                    converged = .true.
                    exit ! leave do loop
                endif

            ! if dp=0, take a bisection step.
            else
                if(p .gt. 0.d0) x1 = (x0 + 1.d0)/2.d0
                if(p .lt. 0.d0) x1 =  x0/2.d0
            endif

            ! update old position
            x0 = x1
        enddo

        ! compute numerator and denominator of s(lambda_c)
        d  = g0
        xx = x1
        n  = f1*xx
        d  = d + g1*xx
        xx = xx*x1
        n  = n + f2*xx
        d  = d + g2*xx
        xx = xx*x1
        n  = n + f3*xx
        d  = d + g3*xx
        xx = xx*x1
        n  = n + f4*xx

        ! compute s(lambda_c), which gets returned as a parameter
        smax = n / d

        ! if we did not find a root, notify the user. if this happens,
        ! something has gone very wrong!
!        if (.not. converged) print *, 'polyecf> newton did not converge'

    end subroutine polyecf

    function max_distance(mols) result(d)
        ! compute the maximum distance of a molecule from the origin

        use vec3, only: vec_len

        implicit none
        type(py_molecule), intent(in) :: mols(:)
        double precision :: d

        double precision :: this_d
        integer :: i

        d = 0.d0

        do i = 1, size(mols)
            this_d = vec_len(mols(i)%r)
            if(this_d > d) d = this_d
        end do

    end function max_distance

    function inverse(m) result(inv)
        ! compute inverse of 3x3 matrix

        implicit none
        double precision, intent(in) :: m(:, :)
        double precision :: inv(3, 3)

        double precision :: adj(3, 3)

        adj(:, :) = adjugate(m)
        inv = adj / (m(1, 1) * adj(1, 1) + m(1, 2) * adj(2, 1) + m(1, 3) * adj(3, 1))

    end function inverse

    function adjugate(m) result(a)
        ! return the adjugate (aka adjoint) of a 3x3 matrix

        ! the adjugate matrix has the properties for 3x3 matrices A:
        !  inverse(A) = adj(A) / det(A)
        !  det(A) = sum_i A_ki adj(A)_ik for k = 1, 2, or 3

        double precision, intent(in) :: m(3, 3)
        double precision :: a(3, 3)

        a(1, 1) = m(2, 2) * m(3, 3) - m(3, 2) * m(2, 3)
        a(1, 2) = m(1, 3) * m(3, 2) - m(3, 3) * m(1, 2)
        a(1, 3) = m(1, 2) * m(2, 3) - m(2, 2) * m(1, 3)
        a(2, 1) = m(2, 3) * m(3, 1) - m(3, 3) * m(2, 1)
        a(2, 2) = m(1, 1) * m(3, 3) - m(3, 1) * m(1, 3)
        a(2, 3) = m(1, 3) * m(2, 1) - m(2, 3) * m(1, 1)
        a(3, 1) = m(2, 1) * m(3, 2) - m(3, 1) * m(2, 2)
        a(3, 2) = m(1, 2) * m(3, 1) - m(3, 2) * m(1, 1)
        a(3, 3) = m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)

    end function adjugate

end module py_module

subroutine py_input
    use commons, only: natoms, site, nrbsites, ntsites
    use py_module
    implicit none
    integer, parameter :: py_unit = 299

    type(py_molecule), pointer :: mol
    type(py_site), pointer :: this_site
    integer :: num_mols, npysites, this_npysite, arity, num_these_mols
    integer :: mol_index, mol_type_index, site_index
    double precision :: number_fraction
    character(len=5) :: label1
    character(len=7) :: label2
    character(len=11) :: label3

    call initialize_py_potential

    num_mols = natoms / 2
    total_num_sites = 0
    allocate(molecules(num_mols))

    !--- parse input
    open(unit=py_unit, file='pysites.xyz', status='old')
    read(py_unit, *) npysites
    if(npysites > 0) then
        ! read in identical particles
        arity = 1
        this_npysite = npysites
        num_these_mols = num_mols

        ! allocate the SITE array to store the site information for use with rigid body routines
        nrbsites = npysites
    elseif(npysites == 0) then
        ! read in an n-ary mixture
        read(py_unit, *) arity

        ! pretend that each molecule has only one site, since the rigid body routines can't handle n-ary mixtures
        nrbsites = 1
    end if

    ! allocate the SITE array for other rigid body routines
    allocate(site(nrbsites,3))
    ntsites = num_mols * nrbsites

    total_num_sites = 0
    mol_index = 1
    do mol_type_index = 1, arity
        read(py_unit, *) ! blank line
        if(npysites == 0) then
            read(py_unit, *) number_fraction   ! read in number fraction for n-ary mixture
            read(py_unit, *) this_npysite      ! read number of site here

            num_these_mols = nint(num_mols * number_fraction)    ! compute number of this type
        end if

        if(num_these_mols == 0) then
            read(py_unit, *) ! there will be no molecules of this type, so just toss the information
        else
            ! read in all the sites to the first molecule in this set
            allocate(molecules(mol_index)%sites(this_npysite))
            do site_index = 1, this_npysite
                this_site => molecules(mol_index)%sites(site_index)
                read(py_unit, *) label1, this_site%r_molframe(1), this_site%r_molframe(2), this_site%r_molframe(3), &
                    & label2, this_site%rep%strength, this_site%att%strength, &
                    & this_site%rep%semiaxes(1), this_site%rep%semiaxes(2), this_site%rep%semiaxes(3), &
                    & this_site%att%semiaxes(1), this_site%att%semiaxes(2), this_site%att%semiaxes(3), &
                    & label3, this_site%p(1), this_site%p(2), this_site%p(3)

                ! if there is only one type of molecule, then store the site position information in SITE too
                if(arity == 1) site(site_index, :) = this_site%r_molframe(:)
            end do
            mol_index = mol_index + 1

            ! copy the first molecule into the rest of the molecules of this type
            do mol_index = mol_index, mol_index + (num_these_mols - 2)
                allocate(molecules(mol_index)%sites(this_npysite))
                molecules(mol_index) = molecules(mol_index - 1)
            end do

            ! update the total site count
            total_num_sites = total_num_sites + num_these_mols * this_npysite
        end if
    end do
    !--- end parsing input

    ! initialize
    do mol_index = 1, size(molecules)
        mol => molecules(mol_index)

        ! initialize molecule values
        mol%r_index = 3 * mol_index - 2
        mol%p_index = mol%r_index + 3 * num_mols

        do site_index = 1, size(mol%sites)
            this_site => mol%sites(site_index)
            this_site%molecule => mol    ! give all sites a pointer to their parent

            call initialize_py_ellipsoid(this_site, this_site%rep)

            if(all(this_site%att%semiaxes == this_site%rep%semiaxes) &
              & .and. (this_site%att%strength == this_site%rep%strength)) then
                this_site%ells_same = .true.
                this_site%att = this_site%rep
            else
                this_site%ells_same = .false.
                call initialize_py_ellipsoid(this_site, this_site%att)
            end if
                
        end do
    end do

end subroutine py_input

subroutine py(x, grad, energy, gtest)
    use commons, only: natoms, vt, frozen
    use py_module
    implicit none
    logical, intent(in) :: gtest
    double precision, intent(inout) :: x(3 * natoms)
    double precision, intent(out) :: energy, grad(3 * natoms)

    type(py_molecule), pointer :: moli, molj
    type(py_site), pointer :: sitei, sitej
    integer :: num_mols
    integer :: moli_index, molj_index, sitei_index, sitej_index
    double precision :: energy_contrib, grad_contrib(12)

    num_mols = size(molecules)

    energy = 0.d0
    vt(:) = 0.d0
    if(gtest) grad(:) = 0.d0

    ! update the values for all the molecules
    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        call update_py_molecule(moli, x(moli%r_index: moli%r_index + 2), &
            & x(moli%p_index: moli%p_index + 2), gtest)
    end do

    ! outer loop over molecules
    do moli_index = 1, num_mols - 1
        moli => molecules(moli_index)

        ! inner loop over molecules
        do molj_index = moli_index + 1, num_mols
            molj => molecules(molj_index)

            ! loop over sites in outer molecule
            do sitei_index = 1, size(moli%sites)
                sitei => moli%sites(sitei_index)

                ! loop over sites in inner molecule
                do sitej_index = 1, size(molj%sites)
                    sitej => molj%sites(sitej_index)

                    energy_contrib = 0.d0
                    grad_contrib(:) = 0.d0

                    ! compute the energy and gradient for this pair of sites
                    call pairwise_py(sitei, sitej, grad_contrib, energy_contrib, gtest)

                    ! add the energy to total and pairwise
                    energy = energy + energy_contrib
                    vt(moli_index) = vt(moli_index) + energy_contrib
                    vt(molj_index) = vt(molj_index) + energy_contrib

                    if(gtest) then  ! add gradient contributions
                        grad(moli%r_index: moli%r_index + 2) = grad(moli%r_index: moli%r_index + 2) + grad_contrib(1: 3)
                        grad(molj%r_index: molj%r_index + 2) = grad(molj%r_index: molj%r_index + 2) + grad_contrib(4: 6)
                        grad(moli%p_index: moli%p_index + 2) = grad(moli%p_index: moli%p_index + 2) + grad_contrib(7: 9)
                        grad(molj%p_index: molj%p_index + 2) = grad(molj%p_index: molj%p_index + 2) + grad_contrib(10: 12)
                    end if

                end do
            end do
        end do
    end do

    ! freeze keyword: set gradients to zero
    if(gtest .and. any(frozen)) then
        do moli_index = 1, num_mols
            moli => molecules(moli_index)
            if(frozen(moli_index)) grad(moli%r_index: moli%r_index + 2) = 0.d0
            if(frozen(moli_index + num_mols)) grad(moli%p_index: moli%p_index + 2) = 0.d0
        end do
    end if

end subroutine py

subroutine py_output(output_file)
    ! writes ellipsoid.xyz for visualization

    ! nsave = number of saved minima
    ! qmin[min] = energy of saved minimum
    ! qminp[min][xyz] = r and p coordinates of saved minima (as per x array)
    ! ff[min] = first step when minimum was found
    use commons, only: nsave, mpit, mynode
    use qmodule, only: qmin, qminp, ff
    use py_module
    implicit none
    integer, intent(in) :: output_file  ! unit number for output

    integer :: min_num, mol_num, site_num   ! indices for loops
    type(py_molecule), pointer :: mol       ! pointers for loops
    type(py_site), pointer :: site
    type(py_ellipsoid), pointer :: ell
    double precision :: mat_out(3, 3)   ! rotation matrix used in the output
    character(len=25) :: filename, node

    ! open file for writing
    if(mpit) then
        write(node, *) mynode + 1
        filename = 'ellipsoid.' // trim(adjustl(node)) // '.xyz'
    else
        filename = 'ellipsoid.xyz'
    end if
    open(unit=output_file, file=filename, status='unknown')

    ! loop over saved minima
    do min_num = 1, nsave
        ! put number of ellipsoids and comment line
        write(output_file, *) total_num_sites
        write(output_file, '(A,I10,A,F20.10,A,I10)') 'Energy of minimum', min_num, ' =', qmin(min_num), ' first found at step', ff(min_num)

        ! loop over molecules
        do mol_num = 1, size(molecules)
            mol => molecules(mol_num)

            ! use the saved r and p values to recreate the saved configuration
            call update_py_molecule(mol, qminp(min_num, mol%r_index: mol%r_index + 2), &
                & qminp(min_num, mol%p_index: mol%p_index + 2), .false.)

            ! loop over sites
            do site_num = 1, size(mol%sites)
                site => mol%sites(site_num)
                ell => site%rep     ! use the repulsive ellipsoid for visualization

                mat_out = matmul(mol%rotmat, site%rotmat)

                write(output_file, '(a5, 2x, 3f20.10, 2x, a8, 12f15.8, 2x, a11, 3f15.8)') &
                    & 'O', site%r(1), site%r(2), site%r(3), &
                    & 'ellipse', 2.d0 * ell%semiaxes(1), 2.d0 * ell%semiaxes(2), 2.d0 * ell%semiaxes(3), &
                    & mat_out(1, 1), mat_out(1, 2), mat_out(1, 3), &
                    & mat_out(2, 1), mat_out(2, 2), mat_out(2, 3), &
                    & mat_out(3, 1), mat_out(3, 2), mat_out(3, 3), &
                    & 'atom_vector', mol%p(1), mol%p(2), mol%p(3) 

            end do
        end do
    end do

    close(output_file)

end subroutine py_output

subroutine py_takestep(np)
    use commons, only: myunit, coords, radius, percolatet, perccut, tmove, omove, step, ostep, astep, vt, frozen
    use py_module
    use vec3, only: vec_len, vec_random
    use rotations, only: rot_takestep_aa

    implicit none
    integer, intent(in) :: np

    double precision, target :: x(size(coords(:, np)))
    double precision, pointer :: r(:), p(:)
    double precision :: step_size, min_vt, dice_roll
    type(py_molecule), pointer :: moli, molj
    type(py_site), pointer :: sitei, sitej
    integer :: num_mols, moli_index, molj_index, sitei_index, sitej_index, move_fails
    logical :: stop_pulling ! for angular moves

    num_mols = size(molecules)
    x(:) = coords(:, np)
    min_vt = minval(vt)

    ! update ellipsoids before taking any steps
    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        r => x(moli%r_index: moli%r_index + 2)
        p => x(moli%p_index: moli%p_index + 2)

        if(.not. percolatet) then
            if(vec_len(moli%r) > radius) then
                write(myunit, *) 'py_takestep> initial coord outside container, brining in'
                r(:) = r(:) - sqrt(radius) * nint(r(:) / sqrt(radius))
            end if
        end if

        call update_py_molecule(moli, r, p, .false.)
    end do

    do moli_index = 1, num_mols
        moli => molecules(moli_index)
        r => x(moli%r_index: moli%r_index + 2)
        p => x(moli%p_index: moli%p_index + 2)

        ! initialize fail counter. provide landing spot for overlap aborts. if
        ! the number of fails is too high, just leave this rigid body where it
        ! is.
        move_fails = 0
95      continue
        if(move_fails > 500) cycle

        ! angular move: accept step with probability 1 - exp(-beta * delta_energy), where beta = 1/kT = astep
        ! if xyz or angular parts are frozen, do not move
        if(.not.(frozen(moli_index) .or. frozen(moli_index + num_mols))) then
            if(astep(np) > 0.d0) then
                call random_number(dice_roll)
                if(dice_roll > exp(-astep(np) * (vt(moli_index) - min_vt))) then
                    r = max_distance(molecules) * vec_random()

                    ! if using percolate, then pull the molecule in until it is within percolation distance of another
                    if(percolatet) then
                        stop_pulling = .false.
                        do
                            do molj_index = 1, num_mols
                                if(moli_index == molj_index) cycle
                                molj => molecules(molj_index)

                                if(vec_len(moli%r - molj%r) < perccut) then
                                    stop_pulling = .true.
                                    exit
                                end if
                            end do

                            if(stop_pulling) exit
                            r(:) = 0.95 * r(:)
                        end do
                    else ! if using radius, pull in until the molecule is inside the container
                        do
                            if(vec_len(r) < radius) exit
                            r(:) = 0.95 * r(:)
                        end do
                    end if

                end if
            end if
        end if

        ! translational move: uniformly distributed in a sphere of radius step(np) (only if xyz not frozen)
        if(.not. frozen(moli_index)) then
            if(tmove(np)) then
                call random_number(step_size)
                step_size = sqrt(step(np) * step_size)  ! need to take square root to get uniform volume distribution
                r(:) = r(:) + step_size * vec_random()
            end if
        end if

        ! orientational move (only if angle part not frozen)
        if(.not. frozen(moli_index + num_mols)) then
            if(omove(np)) then
                call random_number(step_size)
                call rot_takestep_aa(p(:), step_size * ostep(np))
            end if
        end if

        ! update molecule, then check for overlap with other molecules
        call update_py_molecule(moli, r, p, .false.)

        do molj_index = 1, num_mols
            if(moli_index == molj_index) cycle
            molj => molecules(molj_index)

            do sitei_index = 1, size(moli%sites)
                sitei => moli%sites(sitei_index)

                do sitej_index = 1, size(molj%sites)
                    sitej => molj%sites(sitej_index)

                    ! if there is overlap, reset the coordinates, increment the
                    ! fail counter and try again
                    if(py_overlap(sitei%rep, sitej%rep)) then
                        r = coords(moli%r_index: moli%r_index + 2, np)
                        p = coords(moli%p_index: moli%p_index + 2, np)
                        move_fails = move_fails + 1
                        go to 95
                    end if

                end do
            end do
        end do

!        print *, 'succeeded with fails: ', move_fails
    end do

    ! copy local results back
    coords(:, np) = x(:)

end subroutine py_takestep


subroutine rigidcontour
 use commons, only : natoms, contourbounds, gridsize 
implicit none
   
   double precision :: x(3*natoms), grad(3*natoms), stepsize(3)
   double precision :: energymatrix(gridsize,gridsize,gridsize),ereal
   integer          :: ix, iy, iz, myunit3, myunit4, getunit
   logical          :: gradt
    
    x(:)=0.0D0
    x(9)=6.283185D0
    x(10)= 0.80045706105333891411
    x(11)=-1.41840695812026718059 
    x(12)=-0.80045807397906165725
!    x(12)=x(9)
    stepsize(:)=(contourbounds(:,2)-contourbounds(:,1))/gridsize
    myunit3=getunit()
    myunit4=getunit()+1
    open(unit=myunit3,file='energymatrix',status='unknown')
    open(unit=myunit4,file='energycurve',status='unknown')
    ix=1
    iy=1
    iz=1
    do ix=1,gridsize
      do iz=1,gridsize
        x(4)=contourbounds(1,1)+ix*stepsize(1)
        x(6)=contourbounds(3,1)+iz*stepsize(3)
        call py(x, grad, ereal, gradt)
        energymatrix(ix,iy,iz)=ereal
       write(myunit3,*) x(4), x(6), energymatrix(ix,iy,iz)
        if(x(4)==x(6).and.x(4)>0.0D0) write(myunit4,*) sqrt(x(4)**2+x(6)**2), energymatrix(ix,iy,iz) 
!       write(myunit3,*) x(:), energymatrix(ix,iy,iz)
      end do
       write(myunit3,*)
    end do
end subroutine rigidcontour

