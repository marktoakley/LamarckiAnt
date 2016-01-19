C   GPL License info {{{
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C }}}
C
C Doxygen: g46merdiff {{{
C
C> \mainpage 
C> \name g46merdiff

C
C> \brief Calculate the energy, gradient, and second derivatives matrix for the Go-like BLN model \n
C> \author John Rose
C>
C> A particle configuration and number of particles is passed to the subroutine and
C> the energy, gradient, and matrix of second derivatives is returned.
C>
C> \param N          number of particles
C> \param QO         array of cartesian particle coordinates
C> \param GRAD       array of gradients
C> \param ENERGY     energy
C }}}
        subroutine g46merdiff(qo, n, grad, energy, gtest)
C {{{ 
C declarations {{{
        USE MODHESS
        IMPLICIT NONE
        logical gtest, stest
        INTEGER ntype(46), N
        DOUBLE PRECISION QO(3*N), GRAD(3*N), ENERGY
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N),D_PARAM(N),
     1                   c_param(n), rk_theta, rk_r, epsilon, sigma, theta_0, delta, rmass
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        DOUBLE PRECISION X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n)
C }}}
C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    1  d_param(n),c_param(n)

        STEST=.FALSE.

        call gparam_array(a_param,b_param,c_param,d_param,n)
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        call calc_energy(qo,energy,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        call calc_gradient(qo,grad,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

C commented section {{{
C        DIF=1.0D-4
C        DO J1=1,3*N
C           TEMP1=QO(J1)
C           QO(J1)=QO(J1)+DIF
C           call calc_int_coords(qo,n)
C           call calc_energy(qo,V1,n)
C           QO(J1)=QO(J1)-2.0D0*DIF
C           call calc_int_coords(qo,n)
C           call calc_energy(qo,V2,n)
C           tgrad(J1)=(V1-V2)/(2.0D0*DIF)
C           QO(J1)=TEMP1
C        ENDDO
C        call calc_int_coords(qo,n)

C       PRINT*,'Analytical/Numerical first derivatives:'
C       WRITE(*,'(3G20.10)') (GRAD(J1)/TGRAD(J1),J1=1,3*N)
C }}}

        IF (.NOT.STEST) RETURN
        call calc_dyn(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

        return
        end
C }}}
C
C Doxygen: gparam_array {{{
C>
C> \brief Fill the parameter arrays which specify interaction potentials
C> \param N INTEGER  - number of particles
C> \param a_param \param b_param - LJ interaction between non-bonded particles 
C> \param c_param \param d_param - dihedral angle potential
C>
C}}}
        subroutine gparam_array(a_param,b_param,c_param,d_param,n)
C {{{
C Declarations {{{
        IMPLICIT NONE
        logical connect(46,46)
        INTEGER J, ICOUNT, I, J2, J1, N
        DOUBLE PRECISION NTYPE(46), A_PARAM(N,N), B_PARAM(N,N)
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N), EPSILON
        parameter (epsilon = 0.0100570)
C }}}
C Specify amino acid types by filling in the array ntype(:) {{{

        ntype(1) = 1
        ntype(2) = 1
        ntype(3) = 1
        ntype(4) = 1
        ntype(5) = 1
        ntype(6) = 1
        ntype(7) = 1
        ntype(8) = 1
        ntype(9) = 1
        ntype(10) = 3
        ntype(11) = 3
        ntype(12) = 3
        ntype(13) = 2
        ntype(14) = 1
        ntype(15) = 2
        ntype(16) = 1
        ntype(17) = 2
        ntype(18) = 1
        ntype(19) = 2
        ntype(20) = 1
        ntype(21) = 3
        ntype(22) = 3
        ntype(23) = 3
        ntype(24) = 1
        ntype(25) = 1
        ntype(26) = 1
        ntype(27) = 1
        ntype(28) = 1
        ntype(29) = 1
        ntype(30) = 1
        ntype(31) = 1
        ntype(32) = 1
        ntype(33) = 3
        ntype(34) = 3
        ntype(35) = 3
        ntype(36) = 2
        ntype(37) = 1
        ntype(38) = 2
        ntype(39) = 1
        ntype(40) = 2
        ntype(41) = 1
        ntype(42) = 2
        ntype(43) = 1
        ntype(44) = 2
        ntype(45) = 1
        ntype(46) = 2
     
C }}}
C Go-like model connectivities: fill in array CONNECT(:,:) {{{
C
        DO J1=1,46
           DO J2=J1,46
              CONNECT(J2,J1)=.FALSE.
           ENDDO
        ENDDO
        CONNECT(20, 1)=.TRUE.
        CONNECT(24, 1)=.TRUE.
        CONNECT(45, 1)=.TRUE.
        CONNECT(24, 2)=.TRUE.
        CONNECT(43, 2)=.TRUE.
        CONNECT(45, 2)=.TRUE.
        CONNECT(18, 3)=.TRUE.
        CONNECT(20, 3)=.TRUE.
        CONNECT(25, 3)=.TRUE.
        CONNECT(26, 3)=.TRUE.
        CONNECT(43, 3)=.TRUE.
        CONNECT(26, 4)=.TRUE.
        CONNECT(41, 4)=.TRUE.
        CONNECT(16, 5)=.TRUE.
        CONNECT(18, 5)=.TRUE.
        CONNECT(26, 5)=.TRUE.
        CONNECT(27, 5)=.TRUE.
        CONNECT(28, 5)=.TRUE.
        CONNECT(41, 5)=.TRUE.
        CONNECT(28, 6)=.TRUE.
        CONNECT(39, 6)=.TRUE.
        CONNECT(16, 7)=.TRUE.
        CONNECT(28, 7)=.TRUE.
        CONNECT(29, 7)=.TRUE.
        CONNECT(30, 7)=.TRUE.
        CONNECT(39, 7)=.TRUE.
        CONNECT(30, 8)=.TRUE.
        CONNECT(37, 8)=.TRUE.
        CONNECT(14, 9)=.TRUE.
        CONNECT(30, 9)=.TRUE.
        CONNECT(31, 9)=.TRUE.
        CONNECT(32, 9)=.TRUE.
        CONNECT(37, 9)=.TRUE.
        CONNECT(30, 14)=.TRUE.
        CONNECT(31, 14)=.TRUE.
        CONNECT(28, 16)=.TRUE.
        CONNECT(29, 16)=.TRUE.
        CONNECT(26, 18)=.TRUE.
        CONNECT(24, 20)=.TRUE.
        CONNECT(25, 20)=.TRUE.
        CONNECT(45, 24)=.TRUE.
        CONNECT(41, 26)=.TRUE.
        CONNECT(43, 26)=.TRUE.
        CONNECT(39, 28)=.TRUE.
        CONNECT(41, 28)=.TRUE.
        CONNECT(39, 30)=.TRUE.
        CONNECT(37, 32)=.TRUE.
        CONNECT(1, 20)=.TRUE.
        CONNECT(1, 24)=.TRUE.
        CONNECT(1, 45)=.TRUE.
        CONNECT(2, 24)=.TRUE.
        CONNECT(2, 43)=.TRUE.
        CONNECT(2, 45)=.TRUE.
        CONNECT(3, 18)=.TRUE.
        CONNECT(3, 20)=.TRUE.
        CONNECT(3, 25)=.TRUE.
        CONNECT(3, 26)=.TRUE.
        CONNECT(3, 43)=.TRUE.
        CONNECT(4, 26)=.TRUE.
        CONNECT(4, 41)=.TRUE.
        CONNECT(5, 16)=.TRUE.
        CONNECT(5, 18)=.TRUE.
        CONNECT(5, 26)=.TRUE.
        CONNECT(5, 27)=.TRUE.
        CONNECT(5, 28)=.TRUE.
        CONNECT(5, 41)=.TRUE.
        CONNECT(6, 28)=.TRUE.
        CONNECT(6, 39)=.TRUE.
        CONNECT(7, 16)=.TRUE.
        CONNECT(7, 28)=.TRUE.
        CONNECT(7, 29)=.TRUE.
        CONNECT(7, 30)=.TRUE.
        CONNECT(7, 39)=.TRUE.
        CONNECT(8, 30)=.TRUE.
        CONNECT(8, 37)=.TRUE.
        CONNECT(9, 14)=.TRUE.
        CONNECT(9, 30)=.TRUE.
        CONNECT(9, 31)=.TRUE.
        CONNECT(9, 32)=.TRUE.
        CONNECT(9, 37)=.TRUE.
        CONNECT(14, 30)=.TRUE.
        CONNECT(14, 31)=.TRUE.
        CONNECT(16, 28)=.TRUE.
        CONNECT(16, 29)=.TRUE.
        CONNECT(18, 26)=.TRUE.
        CONNECT(20, 24)=.TRUE.
        CONNECT(20, 25)=.TRUE.
        CONNECT(24, 45)=.TRUE.
        CONNECT(26, 41)=.TRUE.
        CONNECT(26, 43)=.TRUE.
        CONNECT(28, 39)=.TRUE.
        CONNECT(28, 41)=.TRUE.
        CONNECT(30, 39)=.TRUE.
        CONNECT(32, 37)=.TRUE.

C }}}
C Parameters for the dihedral angle potential: fill in arrays c_param(:,:), d_param(:,:) {{{

        do i = 1, n-3
        icount = 0

        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo

        if(icount .ge. 2)then
        c_param(i+1) = 0.0
        d_param(i+1) = 0.2*epsilon
        else
        c_param(i+1) = 1.2*epsilon
        d_param(i+1) = 1.2*epsilon
        endif

        icount = 0

        enddo
C }}}
C Parameters for the L-J interaction between non-bonded particles:
C arrays a_param(:,:), b_param(:,:)
C {{{

        do i = 1, n-1
           do j = i+1, n

           if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3) then
             a_param(i,j) = 1.0*epsilon 
             b_param(i,j) = 0.0 
             a_param(j,i) = 1.0*epsilon 
             b_param(j,i) = 0.0
           elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
             a_param(i,j) =  epsilon
             a_param(j,i) =  epsilon
             IF (CONNECT(I,J)) THEN
                b_param(i,j) = -epsilon 
                b_param(j,i) = -epsilon
             ELSE
                b_param(i,j) = 0.0D0
                b_param(j,i) = 0.0D0
             ENDIF
           else
             a_param(i,j) = epsilon*2.0/3.0 
             b_param(i,j) = epsilon*2.0/3.0 
             a_param(j,i) = epsilon*2.0/3.0 
             b_param(j,i) = epsilon*2.0/3.0 
           endif
   
           enddo
        enddo
C }}}
        return
        end
C }}}
