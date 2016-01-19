module chiro_module
    
    implicit none

    double precision, parameter :: z_hat(3) = (/ 0.0, 0.0, 1.0 /)

    type chiro_molecule
        double precision :: r(3), p(3), rm(3, 3), rmd(3, 3, 3)  ! pos, orient, rot mat, rot mat derivs
        double precision :: mu_hat0(3), mu_hat(3), mu_hat_deriv(3, 3)
        integer :: rind, pind   ! indices of values in the big X array
        ! could put an LJ sites array in here
    end type chiro_molecule

    integer :: nmols, x_size    ! number of molecules, size of x array

    type(chiro_molecule), allocatable, target :: molecules(:)

    double precision :: chiro_sigma, chiro_mu, chiro_gamma, chiro_l
    double precision :: mu_p2, sin_gamma, cos_gamma

contains
    subroutine update_chiro_molecule(mol, r, p, gtest)
        use commons, only: periodic, boxlx
        implicit none
        type(chiro_molecule), target, intent(inout) :: mol
        double precision, intent(inout) :: r(3), p(3)
        logical, intent(in) :: gtest

        double precision, parameter :: pi = 3.14159265359
        double precision :: pmod
        integer :: k
        double precision :: dummy(3, 3)

        mol%r(:) = r(:)

        ! make sure that 0 < |p| < 2*pi
        pmod = sqrt(dot_product(p, p))
        if(pmod > 2 * pi) p(:) = p(:) / pmod * mod(pmod, 2 * pi)
        mol%p(:) = p(:)

        ! recompute the rotation matrices and derivatives
        call rmdrvt(mol%p, mol%rm, mol%rmd(1, :, :), mol%rmd(2, :, :), mol%rmd(3, :, :), gtest)

        mol%mu_hat = matmul(mol%rm, mol%mu_hat0)

        ! dmu/dp^I_k = dR^I/dp^I_k . mu_hat0
        do k = 1, 3
            mol%mu_hat_deriv(k, :) = matmul(mol%rmd(k, :, :), mol%mu_hat0)
        end do

        ! check for periodic box condition
        if(periodic .and. abs(mol%r(1)) > boxlx / 2) then
            mol%r(1) = mol%r(1) - sign(boxlx, mol%r(1))
        end if

    end subroutine update_chiro_molecule

    subroutine rod_mindist2(ri, rj, mu_hati, mu_hatj, mui_deriv, muj_deriv, li, lj, dij_p2, xi, xj, grad, gtest)
        implicit none
        ! ri = position of center of first rod
        ! mu_hati = unit vector pointing along first rod
        ! li = half the length of the first rod
        double precision, intent(in) :: ri(3), rj(3), mu_hati(3), mu_hatj(3), mui_deriv(3, 3), muj_deriv(3, 3), li, lj
        logical, intent(in) :: gtest

        ! dij_p2 = square of minimum contact distance
        ! xi = point of closest approach on first rod
        double precision, intent(out) :: dij_p2, xi(3), xj(3), grad(12)

        double precision :: rij(3), rij_modp2   ! |rij|^2
        double precision :: r_dot_mui, r_dot_muj, mu_dot_mu
        double precision :: cc
        double precision :: lambdai, lambdaj ! position of point of contact as measured along rod from the center
        double precision :: offi, offj ! position of naive point of contact as measured from end of rod
        double precision :: xij(3)

        double precision :: dd2_dri(3), dd2_dpi(3), dd2_dpj(3), dd2_dl
        integer :: k

        logical :: inti, intj

        ! assume our points are in the interior
        inti = .true.
        intj = .true.

        rij(:) = ri - rj
        rij_modp2 = dot_product(rij, rij)
        r_dot_mui = dot_product(rij, mu_hati)
        r_dot_muj = dot_product(rij, mu_hatj)
        mu_dot_mu = dot_product(mu_hati, mu_hatj)
        cc = 1.d0 - mu_dot_mu * mu_dot_mu

        ! compute naive values for the lambdas
        lambdai = (r_dot_mui - mu_dot_mu * r_dot_muj) / cc
        lambdaj = (-r_dot_muj + mu_dot_mu * r_dot_mui) / cc

        if(cc == 0.d0) then
            ! mui is parallel to muj
            inti = .false.
            if(r_dot_mui == 0.d0) then
                ! rods are side by side: there is no displacement along the
                ! rods. choose the contact point as the middle
                lambdai = 0.d0
                lambdaj = 0.d0
                intj = .false.
            else
                ! two rods are not side by side. see below for comments.
                lambdai = sign(li, r_dot_mui)
                lambdaj = lambdai * mu_dot_mu - r_dot_muj
                if(abs(lambdaj) > lj) then
                    lambdaj = sign(lj, lambdaj)
                    intj = .false.
                end if
            end if
        else
            if(.not. (abs(lambdai) < li .and. abs(lambdaj) < lj)) then
                ! lambdas are not in the interior of the rods. compute how far off the end of the rods
                ! the naive points lie.
                offi = abs(lambdai) - li
                offj = abs(lambdaj) - lj

                if(offi > offj) then
                    ! first point of contact is farther off than second, so put first point of contact at
                    ! end of the first rod.
                    lambdai = sign(li, lambdai)
                    inti = .false.

                    ! recompute second point of contact
                    lambdaj = lambdai * mu_dot_mu - r_dot_muj

                    if(abs(lambdaj) > lj) then
                        ! second point of contact is off the end of the rod, so put it at end of rod
                        lambdaj = sign(lj, lambdaj)
                        intj = .false.
                    end if
                else
                    ! do above, but reverse i and j
                    lambdaj = sign(lj, lambdaj)
                    intj = .false.

                    lambdai = lambdaj * mu_dot_mu + r_dot_mui
                    if(abs(lambdai) > li) then
                        lambdai = sign(li, lambdai)
                        inti = .false.
                    end if
                end if
            end if
        end if

        ! Lagos & Sevilla compute dij_p2 this way:
        !dij_p2 = rij_modp2 + lambdai * lambdai + lambdaj * lambdaj - 2.d0 * lambdai * lambdaj * mu_dot_mu &
        !    & + 2.d0 * lambdaj * r_dot_muj - 2.d0 * lambdai * r_dot_mui
        
        xi = ri - lambdai * mu_hati
        xj = rj - lambdaj * mu_hatj
        xij = xi - xj
        dij_p2 = dot_product(xij, xij)

        if(gtest) then
            ! D(d^2)/D(r^i) = d(d^2)/d(r^i) + d(d^2)/dlambda^i . d(lambda^i)/dr^i + d(d^2)/dlambda^j . d(lambda^j)/dr^i
            ! D = derivative, d = partial derivative

            ! add on first term
            dd2_dri(:) = 2.d0 * xij(:)
            do k = 1, 3
                dd2_dpi(k) = dot_product(mui_deriv(k, :), -2.d0 * lambdai * xij(:))
                dd2_dpj(k) = dot_product(muj_deriv(k, :), 2.d0 * lambdaj * xij(:))
            end do

            if(inti) then
                ! compute d(d^2)/d(lambda^i)
                dd2_dl = 2.d0 * (lambdai - lambdaj * mu_dot_mu - r_dot_mui)

                ! if at end point, then d(lambda^i)/d(r^i) = 0. if in interior, add on second term.
                dd2_dri(:) = dd2_dri(:) + dd2_dl * (mu_hati - mu_dot_mu * mu_hatj) / cc
                do k = 1, 3
                    dd2_dpi(k) = dd2_dpi(k) + dd2_dl * (2.d0 * mu_dot_mu * dot_product(mui_deriv(k, :), mu_hatj) * &
                        & (r_dot_mui - mu_dot_mu * r_dot_muj) / (cc * cc) + (dot_product(rij, mui_deriv(k, :)) - &
                        & dot_product(mui_deriv(k, :), mu_hatj) * r_dot_muj) / cc )
                    dd2_dpj(k) = dd2_dpj(k) + dd2_dl * (2.d0 * mu_dot_mu * dot_product(mu_hati, muj_deriv(k, :)) * &
                        & (r_dot_mui - mu_dot_mu * r_dot_muj) / (cc * cc) - (dot_product(mu_hati, muj_deriv(k,:)) * &
                        & r_dot_muj + mu_dot_mu * dot_product(rij, muj_deriv(k, :))) / cc )
                end do
            end if

            if(intj) then
                ! procedure as above, but for lambda^j
                dd2_dl = 2.d0 * (lambdaj - lambdai * mu_dot_mu + r_dot_muj)
                dd2_dri(:) = dd2_dri(:) + dd2_dl * (-mu_hatj + mu_dot_mu * mu_hati) / cc
                do k = 1, 3
                    dd2_dpi(k) = dd2_dpi(k) + dd2_dl * (2.d0 * mu_dot_mu * dot_product(mui_deriv(k, :), mu_hatj) * &
                        & (-r_dot_muj + mu_dot_mu * r_dot_mui) / (cc * cc) + (dot_product(mui_deriv(k, :) , mu_hatj) * &
                        & r_dot_mui + mu_dot_mu * dot_product(rij, mui_deriv(k, :))) / cc )
                    dd2_dpj(k) = dd2_dpj(k) + dd2_dl * (2.d0 * mu_dot_mu * dot_product(mu_hati, muj_deriv(k, :)) * &
                        & (-r_dot_muj + mu_dot_mu * r_dot_mui) / (cc * cc) + (-dot_product(rij, muj_deriv(k, :)) + &
                        & dot_product(mu_hati, muj_deriv(k, :)) * r_dot_mui) / cc)
                end do
            end if

            grad(1:3) = dd2_dri
            grad(4:6) = -dd2_dri
            grad(7:9) = dd2_dpi
            grad(10:12) = dd2_dpj
        end if

        return
    end subroutine rod_mindist2

    function max_distance(mols) result(d)
        ! compute the maximum distance of a molecule from the origin

        use vec3, only: vec_len

        implicit none
        type(chiro_molecule), intent(in) :: mols(:)
        double precision :: d

        double precision :: this_d
        integer :: i

        d = 0.0

        do i = 1, size(mols)
            this_d = vec_len(mols(i)%r)
            if(this_d > d) d = this_d
        end do

    end function max_distance

end module chiro_module

    subroutine initialize_chiro(out_unit)
        ! for n-ary mixture or multisite models, this subroutine would have to read data from a file
        ! and do the appropriate allocations and assignments to the molecules array
        use chiro_module
        use commons, only: natoms
        implicit none
        
        integer, intent(in) :: out_unit ! for GMIN_out

        type(chiro_molecule), pointer :: mol
        integer :: mol_ind

        nmols = natoms / 2
        x_size = 3 * natoms
        write(out_unit, *) nmols, 'chiropole molecules'

        ! allocate array of molecules
        allocate(molecules(nmols))

        ! put all the indices in the right places
        do mol_ind = 1, nmols
            mol => molecules(mol_ind)
            mol%rind = 3 * mol_ind - 2
            mol%pind = mol%rind + 3 * nmols

            ! set all molecules mu0 to z_hat
            mol%mu_hat0(:) = z_hat(:)

            ! set up for LJ sites array would go here
        end do

        ! process the input parameters
        mu_p2 = chiro_mu * chiro_mu
        sin_gamma = sin(chiro_gamma)
        cos_gamma = cos(chiro_gamma)

    end subroutine initialize_chiro

subroutine chiro(x, grad, energy, gtest)
    use commons, only: vt, twod, periodic, boxlx
    use vec3, only: vec_len, vec_cross
    use chiro_module
    implicit none
    double precision, intent(inout) :: x(x_size)
    logical, intent(in) :: gtest
    double precision, intent(out) :: energy, grad(x_size)

    type(chiro_molecule), pointer :: moli, molj
    integer moli_ind, molj_ind
    double precision :: grad_pair(12)
    double precision :: energy_contrib, grad_contrib(12)

    double precision :: rij(3), rijmod, rij_hat(3)
    double precision :: r_pm1, r_pm3, r_pm4, r_pm6, r_pm7, r_pm12, r_pm13 ! r is in reduced units
    double precision :: mu_dot, mu_cross(3)  ! mu_hati . mu_hatj and mu_hati x mu_hatj
    integer :: k

    double precision :: d_pm1, d_pm2, d_pm6, d_pm12
    double precision :: d2, d2_pm1, d2_pm3, d2_pm6
    double precision :: xi(3), xj(3)

    double precision :: image_moljr(3)
    integer :: lshift, ushift, jshift
    logical :: newpair

    energy = 0.0
    vt(:) = 0.0
    if(gtest) grad(:) = 0.0

    ! update the values for all the molecules
    do moli_ind = 1, nmols
        moli => molecules(moli_ind)
        call update_chiro_molecule(moli, x(moli%rind: moli%rind + 2), x(moli%pind: moli%pind + 2), gtest)
    end do

    if(periodic) then
        lshift = -1
        ushift = 1
    else
        lshift = 0
        ushift = 0
    end if

    ! outer loop over molecules
    do moli_ind = 1, nmols - 1
        moli => molecules(moli_ind)

        ! inner loop over molecules
        do molj_ind = moli_ind + 1, nmols
            molj => molecules(molj_ind)

            ! do pbc shift for inner molecule
            newpair = .true.
            do jshift = lshift, ushift

            ! --- chiropole contributions ---
            ! compute angular parts only once per set of images
            if(newpair) then
                mu_dot = dot_product(moli%mu_hat, molj%mu_hat)
                mu_cross(:) = vec_cross(moli%mu_hat, molj%mu_hat)
                newpair = .false.
            end if

            energy_contrib = 0.0
            grad_contrib(:) = 0.0

            ! move image according to current pbc shift
            image_moljr = molj%r
            if(periodic) image_moljr(1) = image_moljr(1) + jshift * boxlx
            rij(:) = moli%r - image_moljr

            rijmod = vec_len(rij(:))
            if(rijmod == 0.0) then
                rij_hat(:) = 0.0
            else
                rij_hat(:) = rij(:) / rijmod
            end if

            r_pm1 = chiro_sigma / rijmod
            r_pm3 = r_pm1 * r_pm1 * r_pm1

            energy_contrib = energy_contrib &
                & - mu_p2 * r_pm3 * (cos_gamma * mu_dot + sin_gamma * dot_product(mu_cross, rij_hat))

            if(gtest) then
                grad_pair(:) = 0.0
                r_pm4 = r_pm3 * r_pm1

                grad_pair(1: 3) = -mu_p2 * r_pm4 / chiro_sigma * (cos_gamma * (-3.0 * mu_dot * rij_hat(:)) + &
                    & sin_gamma * (mu_cross(:) - 4.0 * dot_product(mu_cross, rij_hat) * rij_hat(:)))
                grad_pair(4: 6) = -grad_pair(1: 3)

                do k = 1, 3
                    grad_pair(6 + k) = -mu_p2 * r_pm3 * (cos_gamma * dot_product(moli%mu_hat_deriv(k, :), molj%mu_hat) & 
                        & + sin_gamma * dot_product(vec_cross(moli%mu_hat_deriv(k, :), molj%mu_hat), rij_hat))

                    grad_pair(9 + k) = -mu_p2 * r_pm3 * (cos_gamma * dot_product(moli%mu_hat, molj%mu_hat_deriv(k, :)) & 
                        & + sin_gamma * dot_product(vec_cross(moli%mu_hat, molj%mu_hat_deriv(k, :)), rij_hat))
                end do

                ! if sites were off-center, would need to do torque terms here

                grad_contrib(:) = grad_contrib(:) + grad_pair(:)
            end if

            if(chiro_l == 0.d0) then
                ! --- lj contributions ---
                d_pm2 = 1.d0 / dot_product(rij, rij)
                d_pm1 = sqrt(d_pm2)
                d_pm6 = d_pm2 ** 3
                d_pm12 = d_pm6 ** 2
                energy_contrib = energy_contrib + 4.d0 * (d_pm12 - d_pm6)

                if(gtest) then
                    grad_pair(1: 3) = 4.d0 * (-12.d0 * d_pm12 + 6.d0 * d_pm6) * d_pm1 * rij_hat(:)
                    grad_pair(4: 6) = -grad_pair(1: 3)
                    grad_pair(7: 12) = 0.d0
                    grad_contrib = grad_contrib + grad_pair
                end if

            elseif(chiro_l > 0.d0) then
                ! --- rod contribution ---
                call rod_mindist2(moli%r, image_moljr, moli%mu_hat, molj%mu_hat, moli%mu_hat_deriv, molj%mu_hat_deriv, &
                    & chiro_l, chiro_l, d2, xi, xj, grad_pair, gtest)
                d2_pm1 = 1.d0 / d2
                d2_pm3 = d2_pm1 * d2_pm1 * d2_pm1
                d2_pm6 = d2_pm3 * d2_pm3

                ! U = 4*(d^-12 - d^-6) = 4*((d^2)^-6 - (d^2)^-3)
                energy_contrib = energy_contrib + 4.d0 * (d2_pm6 - d2_pm3)
                if(gtest) then
                    ! dU/dr^i = 4(-6(d^2)^-7 + 3(d^2)^-4) . d(d^2)/dr^i
                    grad_pair(:) = grad_pair(:) * 4.d0 * (-6.d0 * d2_pm6 + 3.d0 * d2_pm3) * d2_pm1
                    grad_contrib = grad_contrib + grad_pair
                end if

            end if

            ! add energy contribuions
            energy = energy + energy_contrib
            vt(moli_ind) = vt(moli_ind) + energy_contrib
            vt(molj_ind) = vt(molj_ind) + energy_contrib

            ! add gradient contributions
            if(gtest) then
                grad(moli%rind: moli%rind + 2) = grad(moli%rind: moli%rind + 2) + grad_contrib(1: 3)
                grad(molj%rind: molj%rind + 2) = grad(molj%rind: molj%rind + 2) + grad_contrib(4: 6)
                grad(moli%pind: moli%pind + 2) = grad(moli%pind: moli%pind + 2) + grad_contrib(7: 9)
                grad(molj%pind: molj%pind + 2) = grad(molj%pind: molj%pind + 2) + grad_contrib(10: 12)
            end if

            ! increment pbc shift for inner loop
            end do
        end do

    end do

    ! 2D: set all z gradients to zero
    if(twod) then
        do moli_ind = 1, nmols
            grad(molecules(moli_ind)%rind + 2) = 0.0
        end do
    end if

end subroutine chiro

subroutine chiro_output
    use commons, only: natoms, nsave, mpit, mynode
    use qmodule, only: qmin, qminp, ff
    use chiro_module
    implicit none

    integer :: out_unit, coords_unit
    integer :: min_num, mol_num, site_num, atom_num, atom_index
    type(chiro_molecule), pointer :: mol
    character(len=25) :: out_name, coords_name, min_name, node
    double precision :: rotmat(3, 3), dummy(3, 3)

    integer :: getunit

    out_unit = getunit()

    ! open sandout file for writing
    if(mpit) then
        write(node, *) mynode + 1
        out_name = 'chiro.' // trim(adjustl(node)) // '.xyz'
    else
        out_name = 'chiro.xyz'
    end if
    open(unit=out_unit, file=out_name, status='unknown')

    ! loop over saved minima
    do min_num = 1, nsave
        ! put number of atoms and comment line. now we one O for each molecule
        write(out_unit, '(i8)') nmols
        write(out_unit, '(a, i4, a, f20.8, a, i10, a, f6.2)') 'energy of minimum ', min_num, ' = ', qmin(min_num), ' first found at step ', ff(min_num), &
            & ' rod:', chiro_l

        ! loop over molecules
        do mol_num = 1, size(molecules)
            mol => molecules(mol_num)

            call rmdrvt(qminp(min_num, mol%pind: mol%pind + 2), rotmat, dummy, dummy, dummy, .false.)

            ! put an O (red) at the north pole and N (blue) at south pole
            write(out_unit, '(a, f15.6, f15.6, f15.6, a, f15.6, f15.6, f15.6)') 'O ', &
                & qminp(min_num, mol%rind: mol%rind + 2), '      atom_vector ', matmul(rotmat, z_hat)
        end do

        ! output coords. files
        coords_unit = getunit()
        write(min_name, '(i3)') min_num

        if(mpit) then
            coords_name = 'coords.' // trim(adjustl(min_name)) // '.' // trim(adjustl(node))
        else
            coords_name = 'coords.' // trim(adjustl(min_name))
        end if

        open(coords_unit, file=coords_name, status='unknown')
        do atom_num = 1, natoms
            atom_index = 3 * (atom_num - 1) + 1
            write(coords_unit, *) qminp(min_num, atom_index), qminp(min_num, atom_index + 1), qminp(min_num, atom_index + 2)
        end do
        close(coords_unit)

    end do

    close(out_unit)
            
end subroutine chiro_output

subroutine chiro_takestep(np)
    use commons, only: myunit, coords, radius, percolatet, perccut, tmove, omove, step, ostep, astep, vt, twod
    use vec3, only: vec_len, vec_random
    use rotations, only: rot_takestep_aa
    use chiro_module

    implicit none
    integer, intent(in) :: np

    double precision, target :: x(size(coords(:, np)))
    double precision, pointer :: r(:), p(:)
    double precision :: step_size, min_vt, dice_roll, twod_step(3)
    type(chiro_molecule), pointer :: moli, molj
    integer :: moli_ind, molj_ind
    logical :: stop_pulling ! for angular moves

    x(:) = coords(:, np)
    min_vt = minval(vt)

    ! update molecules before taking any steps
    do moli_ind = 1, nmols
        moli => molecules(moli_ind)
        r => x(moli%rind: moli%rind + 2)
        p => x(moli%pind: moli%pind + 2)

        if(.not. percolatet) then
            if(vec_len(moli%r) > radius) then
                write(myunit, *) 'chiro_takestep> initial coord outside container, brining in'
                r(:) = r(:) - sqrt(radius) * nint(r(:) / sqrt(radius))
            end if
        end if

        call update_chiro_molecule(moli, r, p, .false.)
    end do

    do moli_ind = 1, nmols
        moli => molecules(moli_ind)
        r => x(moli%rind: moli%rind + 2)
        p => x(moli%pind: moli%pind + 2)

        if(twod) then
            ! translational step
            do
                call random_number(twod_step(1))
                call random_number(twod_step(2))
                twod_step(3) = 0.0
                if(vec_len(twod_step) <= 1.0) exit
            end do
            r(:) = r(:) + (step(np) * twod_step(:))
        else
            ! angular move: accept step with probability 1 - exp(-beta * delta_energy), where beta = 1/kT = astep
            if(astep(np) > 0.0) then
                call random_number(dice_roll)
                if(dice_roll > exp(-astep(np) * (vt(moli_ind) - min_vt))) then
                    r = max_distance(molecules) * vec_random()

                    ! if using percolate, then pull the molecule in until it is within percolation distance of another
                    if(percolatet) then
                        stop_pulling = .false.
                        do
                            do molj_ind = 1, nmols
                                if(moli_ind == molj_ind) cycle
                                molj => molecules(molj_ind)

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

            ! translational move: uniformly distributed in a sphere of radius step(np)
            if(tmove(np)) then
                call random_number(step_size)
                step_size = sqrt(step(np) * step_size)  ! need to take square root to get uniform volume distribution
                r(:) = r(:) + step_size * vec_random()
            end if
        end if ! 2d

        ! orientational move
        if(omove(np)) then
            call random_number(step_size)
            call rot_takestep_aa(p(:), step_size * ostep(np))
        end if

    end do

    ! copy local results back
    coords(:, np) = x(:)

end subroutine chiro_takestep
