C   GPL License Info {{{
C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C }}}

C> Calculate the energy and gradient
C> for a given configuration of a BLN polymer chain.

      SUBROUTINE CHAPERONIN(QO,GRAD,ENERGY,GRADT)
C{{{
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION QO(3*NATOMS), GRAD(3*NATOMS)
      LOGICAL GRADT
      INTEGER N
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), YR(NATOMS,NATOMS), 
     &                 ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3),
     &                 X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS), 
     &                 ENERGY, DFAC(NATOMS), SINBOND(NATOMS)

      N=NATOMS
!
! Without these initialisations the NAG compiler fills in random numbers for
! unassigned elements with optimisation turned on.
!
      BOND_ANGLE(1:NATOMS)=0.0D0
      TOR_ANGLE(1:NATOMS)=0.0D0
      DFAC(1:NATOMS)=0.0D0
      SINBOND(1:NATOMS)=0.0D0
      DOT_PROD(1:NATOMS,1:3)=0.0D0
      X_PROD(1:NATOMS)=0.0D0
      RADII(1:NATOMS,1:NATOMS)=0.0D0
      XR(1:NATOMS,1:NATOMS)=0.0D0
      YR(1:NATOMS,1:NATOMS)=0.0D0
      ZR(1:NATOMS,1:NATOMS)=0.0D0

      CALL CALC_INT_COORDSCHAPERONIN(QO,N,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD
     $     ,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,SINBOND,A_BLN,B_BLN
     $     ,C_BLN,D_BLN,DFAC)
      CALL CALC_ENERGYCHAPERONIN(QO,ENERGY,N,LJREP_BLN,LJATT_BLN,A_BLN
     $     ,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE
     $     ,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,radius_container
     $     ,HYDRO_BLN)
      IF (.NOT.GRADT) RETURN
 
      CALL CALC_GRADIENTCHAPERONIN(QO,GRAD,N,LJREP_BLN,LJATT_BLN,A_BLN
     $     ,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE
     $     ,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,DFAC,SINBOND
     $     ,RADIUS_CONTAINER,HYDRO_BLN)

      RETURN
      END
C }}}
C
C> Calculate the internal coordinates
C
      SUBROUTINE CALC_INT_COORDSCHAPERONIN(QO,N,X,Y,Z,XR,YR,ZR,DOT_PROD
     $     ,X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,SINBOND ,A_BLN
     $     ,B_BLN,C_BLN,D_BLN,DFAC)
C {{{
C Declarations {{{
      implicit NONE
      INTEGER I, N, J, NATOMS
      DOUBLE PRECISION COS_PHI, COS_THETA, DUMMY, DUMMY2, duma, dumb
      DOUBLE PRECISION QO(3*NATOMS)
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.283185307179586477D0
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS
     $     ), DFAC(NATOMS),YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS),
     $     DOT_PROD(NATOMS,3), D_BLN(NATOMS)
     $     ,A_BLN(NATOMS), B_BLN(NATOMS), X_PROD(NATOMS),
     $     BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS),RADII(NATOMS,NATOMS),
     $     SINBOND(NATOMS), C_BLN(NATOMS)
C }}}
      do i = 1, n
         j = (i-1)*3
         x(i) = qo((i-1)*3+1)
         y(i) = qo((i-1)*3+2)
         z(i) = qo((i-1)*3+3)
      enddo
C
C Inter-particle distances {{{
C
      do i = 1, n-1
         do j = i+1, n
C        do j = 1, n
            xr(i,j) = x(j) - x(i)
            yr(i,j) = y(j) - y(i)
            zr(i,j) = z(j) - z(i)
            radii(i,j) = sqrt(xr(i,j)*xr(i,j) + yr(i,j)*yr(i,j) + zr(i,j)*zr(i,j))
            radii(j,i) = radii(i,j)
         enddo
      enddo
C }}}
C
C Dot products between bond vectors {{{
C

      do i = 1, n-3
         dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + zr(i,i+1)*zr(i,i+1)
         dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ zr(i,i+1)*zr(i+1,i+2)
         dot_prod(i,3) = xr(i,i+1)*xr(i+2,i+3)+yr(i,i+1)*yr(i+2,i+3)+ zr(i,i+1)*zr(i+2,i+3)
      enddo

!     i = n-2
      dot_prod(n-2,1) = xr(n-2,n-2+1)*xr(n-2,n-2+1) + yr(n-2,n-2+1)*yr(n-2,n-2+1) + zr(n-2,n-2+1)*zr(n-2,n-2+1)
      dot_prod(n-2,2) = xr(n-2,n-2+1)*xr(n-2+1,n-2+2)+yr(n-2,n-2+1)*yr(n-2+1,n-2+2)+ zr(n-2,n-2+1)*zr(n-2+1,n-2+2)
!     i = n-1
      dot_prod(n-1,1) = xr(n-1,n-1+1)*xr(n-1,n-1+1) + yr(n-1,n-1+1)*yr(n-1,n-1+1) + zr(n-1,n-1+1)*zr(n-1,n-1+1)
C }}}
C Cross-products between adjacent bond vectors {{{

      do i = 1, n-2
         x_prod(i) = dot_prod(i,1)*dot_prod(i+1,1) - dot_prod(i,2)*dot_prod(i,2)   
      enddo
C }}}
C Bond angles {{{

      do i = 1, n-2
         cos_theta=-dot_prod(i,2)/(sqrt(dot_prod(i,1)*dot_prod(i+1,1)))
         bond_angle(i+1) = dacos(cos_theta)
         SINBOND(i+1)=SIN(bond_angle(i+1))*SQRT(dot_prod(i,1)*dot_prod(i+1,1))
      enddo
C }}}
C Torsional angles {{{

      do i = 1, n-3
         cos_phi = (dot_prod(i,2)*dot_prod(i+1,2) - dot_prod(i,3)
     $        *dot_prod(i+1,1))/sqrt(x_prod(i)*x_prod(i+1))
         IF (ABS(cos_phi).GT.1.0D0) cos_phi=SIGN(1.0D0,cos_phi)
         tor_angle(i+1) = dacos(cos_phi)
C }}}
C
C tor_angle is returned in the range 0 to Pi.
C dummy should take the opposite sign from the dihedral angle. Negative
C values of the dihedral should be subtracted from 2*pi.
C This is only necessary when the potential contains cos(phi+pi/4) terms because
C the gradient is discontinuous when phi goes through pi for such terms if phi
C is restricted to 0 < phi < pi.
C {{{
C        dummy=xr(i+2,i+3)*(yr(i+1,i)*zr(i+1,i+2)-yr(i+1,i+2)*zr(i+1,i))+
C    &         yr(i+2,i+3)*(xr(i+1,i+2)*zr(i+1,i)-xr(i+1,i)*zr(i+1,i+2))+
C    &         zr(i+2,i+3)*(xr(i+1,i)*yr(i+1,i+2)-xr(i+1,i+2)*yr(i+1,i))
         dummy=xr(i+2,i+3)*(-yr(i,i+1)*zr(i+1,i+2)+yr(i+1,i+2)*zr(i,i+1))+
     &         yr(i+2,i+3)*(-xr(i+1,i+2)*zr(i,i+1)+xr(i,i+1)*zr(i+1,i+2))+
     &         zr(i+2,i+3)*(-xr(i,i+1)*yr(i+1,i+2)+xr(i+1,i+2)*yr(i,i+1))
         IF (DUMMY.GT.0.0D0) tor_angle(i+1)=TWOPI-tor_angle(i+1)
C }}}
C  This is an ugly hack to prevent division by zero. There will be a loss of precision
C  if a dihedral is 0, PI, 2*PI, if D_BLN is non-zero.
C {{{
         IF (TAN(tor_angle(i+1)).EQ.0.0D0) THEN
            PRINT '(A,I8,A,G20.10)','WARNING in BLN, dihedral angle ',i+1,' is ',tor_angle(i+1)
            tor_angle(i+1)=tor_angle(i+1)+1.0D-10
            PRINT '(A,G20.10)','WARNING in BLN, TAN perturbed angle=',tor_angle(i+1)
         ENDIF
c         DUMMY2=TAN(tor_angle(i+1))
c
c         DFAC(i+1)=(A_BLN(i+1)+D_BLN(i+1)*( 1.0D0+1.0D0/DUMMY2 )
c     $        *0.7071067811865475244D0-B_BLN(i+1)+C_BLN(i+1)*(12.0
c     $        *costor(i+1)**2-3.0))/sqrt(x_prod(i+1)*x_prod(i))

         duma=A_BLN(i+1)*cos((tor_angle(i+1)-C_BLN(i+1))/2.D0)**5
     $        * sin((tor_angle(i+1)-C_BLN(i+1))/2.D0)
         dumb=B_BLN(i+1)*cos((tor_angle(i+1)-D_BLN(i+1))/2.D0)**5
     $        * sin((tor_angle(i+1)-D_BLN(i+1))/2.D0)

         DFAC(i+1)=-3.D0*(duma+dumb)/sin(tor_angle(i+1))/sqrt(x_prod(i+1
     $        )*x_prod(i))

C }}}
      enddo

      return
      end
C }}}
C> Calculate the energy

      subroutine calc_energyCHAPERONIN(qo,energy,n,LJREP_BLN,LJATT_BLN
     $     ,A_BLN,B_BLN,C_BLN,D_BLN,x,y,z,xr,yr,zr,dot_prod,x_prod
     $     ,bond_angle,tor_angle,radii,natoms,rk_r,rk_theta
     $     ,radius_container,HYDRO_BLN)
C {{{
      implicit NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, parameter :: theta_0 = 1.8326D0, PI4=0.7853981633974483096D0 ! 1.8326 radians is 105 degrees
      DOUBLE PRECISION qo(3*NATOMS), ENERGY
      DOUBLE PRECISION x(NATOMS), y(NATOMS), z(NATOMS), xr(NATOMS,NATOMS
     $     ),yr(NATOMS,NATOMS), zr(NATOMS,NATOMS), dot_prod(NATOMS,3),
     $     x_prod(NATOMS), bond_angle(NATOMS), tor_angle(NATOMS),
     $     radii(NATOMS,NATOMS)
      DOUBLE PRECISION RK_R, RK_THETA, e_nbond, e_bond, e_bangle,
     $     e_tangle, rad6
      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS)
     $     ,C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)
      DOUBLE PRECISION HYDRO_BLN(NATOMS)
      DOUBLE PRECISION e_container, radius_container,pi
      DOUBLE PRECISION radi,dum10p,dum10m,dum4p,dum4m
      PARAMETER(PI=3.141592653589793D0)

      e_nbond=0.0D0
      e_bond=0.0D0
      e_bangle=0.0D0
      e_tangle=0.0D0

      e_container=0.0D0

      do i = 1, n-2
         do j = i+2, n
            rad6 = radii(i,j)**6
            e_nbond = e_nbond + (LJREP_BLN(i,j)/rad6 - LJATT_BLN(i,j))/rad6
         enddo
      enddo
      e_nbond=e_nbond*4.0D0

      do i = 1, n-1
         e_bond = e_bond + (radii(i,i+1)-1.0D0)**2
      enddo
      e_bond=e_bond*rk_r/2.0D0

      do i = 2, n-1
         e_bangle = e_bangle + (bond_angle(i)-theta_0)**2
      enddo
      e_bangle=e_bangle*rk_theta/2.0D0

      do i = 2, n-2
         e_tangle = e_tangle - A_BLN(i)*(cos((tor_angle(i)-C_BLN(i))/2.D0))**6
     &                       - B_BLN(i)*(cos((tor_angle(i)-D_BLN(i))/2.D0))**6
      enddo

      if (RADIUS_CONTAINER.GE.0.D0) THEN
c     container contribution

         do i = 1, n

            radi=sqrt(x(i)**2+y(i)**2+z(i)**2)
            dum10p=1.D0/(radi+RADIUS_CONTAINER)**10
            dum10m=1.D0/(radi-RADIUS_CONTAINER)**10
            dum4p=1.D0/(radi+RADIUS_CONTAINER)**4
            dum4m=1.D0/(radi-RADIUS_CONTAINER)**4
            e_container = e_container +4.D0*pi*RADIUS_CONTAINER*(((dum10m
     $           -dum10p)/5.D0 - hydro_bln(i)*(dum4m-dum4p)/2.D0))/radi
            
         enddo

      else

         e_container=0.D0

      endif

      energy = e_nbond + e_bond + e_bangle + e_tangle + e_container

C     write(*,'(A,4F20.10)') 'nbond,bond,bangle,tangle=',e_nbond,e_bond,e_bangle,e_tangle

      return
      end
C }}}
C
C> Calculate the gradients
C
      SUBROUTINE CALC_GRADIENTCHAPERONIN(QO,FQ,N,LJREP_BLN,LJATT_BLN
     $     ,A_BLN,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD
     $     ,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA
     $     ,DFAC,SINBOND,RADIUS_CONTAINER,HYDRO_BLN)
C {{{
C Declarations {{{
!      USE COMMONS,ONLY : MYUNIT
      IMPLICIT NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, parameter :: theta_0 = 1.8326D0
      DOUBLE PRECISION qo(3*NATOMS),fq(3*NATOMS),fx(NATOMS),fy(NATOMS),
     $     fz(NATOMS), DFAC(NATOMS), SINBOND(NATOMS)
      DOUBLE PRECISION fnb_x(NATOMS),fnb_y(NATOMS),fnb_z(NATOMS),
     $     fb_x(NATOMS),fb_y(NATOMS)
      DOUBLE PRECISION fb_z(NATOMS),fba_x(NATOMS),fba_y(NATOMS)
     $     ,fba_z(NATOMS)
      DOUBLE PRECISION fta_x(NATOMS),fta_y(NATOMS),fta_z(NATOMS), a3,
     $     coef, coef1, coef2, coef3, a4
      DOUBLE PRECISION fc_x(NATOMS),fc_y(NATOMS),fc_z(NATOMS)

      DOUBLE PRECISION RK_R, RK_THETA, rad7, rad14, df, fxx, fzz, fyy,
     $     rvar, den, rnum, den1, a1, a2, den2

      DOUBLE PRECISION x(NATOMS), y(NATOMS), z(NATOMS), xr(NATOMS,NATOMS
     $     ),yr(NATOMS,NATOMS), zr(NATOMS,NATOMS), dot_prod(NATOMS,3)
     $     ,x_prod(NATOMS), bond_angle(NATOMS), tor_angle(NATOMS),
     $     radii(NATOMS,NATOMS)

      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS)
     $     ,C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)

      DOUBLE PRECISION HYDRO_BLN(NATOMS),radius_container
      double precision radi,dum10p,dum10m,dum4p,dum4m
      double precision pi,dummy1,dummy2,dummy3,dummy4
      pi=acos(-1.D0)
C }}}
C
C Gradients of potential
C {{{
      do i = 1,n

         fnb_x(i) = 0.0  
         fnb_y(i) = 0.0 
         fnb_z(i) = 0.0 
   
         fb_x(i)  = 0.0 
         fb_y(i)  = 0.0 
         fb_z(i)  = 0.0 
   
         fba_x(i) = 0.0 
         fba_y(i) = 0.0 
         fba_z(i) = 0.0 
   
         fta_x(i) = 0.0 
         fta_y(i) = 0.0 
         fta_z(i) = 0.0 

         fc_x(i) = 0.0 
         fc_y(i) = 0.0 
         fc_z(i) = 0.0 
   
         fx(i)= 0.0 
         fy(i)= 0.0 
         fz(i)= 0.0 

      enddo
C }}}
C ..... Non-bonded interaction forces ..... 
C {{{
      do i = 1, n-2
         do j = i+2, n

            rad7 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)   
            rad14 = rad7*rad7 

            df = -24.0*((2.0*LJREP_BLN(i,j)/rad14) - (LJATT_BLN(i,j)/(rad7*radii(i,j))))

            fxx = df*xr(i,j) 
            fyy = df*yr(i,j) 
            fzz = df*zr(i,j) 

            fnb_x(i) = fxx + fnb_x(i)
            fnb_y(i) = fyy + fnb_y(i)
            fnb_z(i) = fzz + fnb_z(i)

            fnb_x(j) = -fxx + fnb_x(j)
            fnb_y(j) = -fyy + fnb_y(j)
            fnb_z(j) = -fzz + fnb_z(j)

         enddo
      enddo
C }}}
C ... Bond interaction forces ... 
C {{{
      do i = 1, n-1

         rvar = 1.0D0/radii(i,i+1) 

         df = rk_r*(1.0 - rvar) 
         fxx = df*xr(i,i+1) 
         fyy = df*yr(i,i+1) 
         fzz = df*zr(i,i+1) 

         fb_x(i) = fxx + fb_x(i)
         fb_y(i) = fyy + fb_y(i)
         fb_z(i) = fzz + fb_z(i)

         fb_x(i+1) = -fxx + fb_x(i+1)
         fb_y(i+1) = -fyy + fb_y(i+1)
         fb_z(i+1) = -fzz + fb_z(i+1)

      enddo
C }}}
C bond angle forces  particle 1
C particles 1,2,n-1, and n done outside of the loop
C {{{
!     i = 1
      den = sinbond(1+1)
      rnum = rk_theta*(bond_angle(1+1) - theta_0)

      fba_x(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*xr(1,1+1) - xr(1+1,1+2))/den

      fba_y(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*yr(1,1+1) - yr(1+1,1+2))/den

      fba_z(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*zr(1,1+1) - zr(1+1,1+2))/den

C }}}
C particle 2
C {{{
!     i = 2
      den = sinbond(2)
      den1 = sinbond(3)

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*xr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *xr(2-1,2) + xr(2,2+1) - xr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*xr(2,2+1) - xr(2+1,2+2))/den1

      fba_x(2) = a1 + a2 

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*yr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *yr(2-1,2) + yr(2,2+1) - yr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*yr(2,2+1) - yr(2+1,2+2))/den1

      fba_y(2) = a1 + a2 

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*zr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *zr(2-1,2) + zr(2,2+1) - zr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*zr(2,2+1) - zr(2+1,2+2))/den1

      fba_z(2) = a1 + a2 

C }}}
C particles 3 thru n-2 
C {{{
      do i = 3, n-2

         den = sinbond(i)
         den1 = sinbond(i+1)
         den2 = sinbond(i-1)

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

         fba_x(i) = a1 + a2 + a3 

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1   dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1      dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1      dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

         fba_y(i) = a1 + a2 + a3 

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1   dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1      dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1      dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

         fba_z(i) = a1 + a2 + a3 

      enddo
C }}}
C particle n-1 
C {{{
!     i = n-1
      den = sinbond(n-1)
      den1 = sinbond(n-1-1)

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1   dot_prod(n-1,1))*xr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *xr(n-1-1,n-1) + xr(n-1,n-1+1) - xr(n-1-1,n-1))/den

      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*xr(n-1-1,n-1) - xr(n-1-2,n-1-1))/den1

      fba_x(n-1) = a1 + a2

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1   dot_prod(n-1,1))*yr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *yr(n-1-1,n-1) + yr(n-1,n-1+1) - yr(n-1-1,n-1))/den
   
      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*yr(n-1-1,n-1) - yr(n-1-2,n-1-1))/den1

      fba_y(n-1) = a1 + a2

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1  dot_prod(n-1,1))*zr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *zr(n-1-1,n-1) + zr(n-1,n-1+1) - zr(n-1-1,n-1))/den

      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*zr(n-1-1,n-1) - zr(n-1-2,n-1-1))/den1

      fba_z(n-1) = a1 + a2
C }}}
C particle n
C {{{
!     i = n
      den = sinbond(n-1)

      fba_x(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*xr(n-1,n) 
     1      - xr(n-2,n-1))/den

      fba_y(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*yr(n-1,n) 
     1      - yr(n-2,n-1))/den

      fba_z(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*zr(n-1,n) 
     1      - zr(n-2,n-1))/den

C }}}
C Torsional angle forces
C particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
C particle 1
C {{{
!     i = 1
      coef =DFAC(1+1)

      fta_x(1) = -coef*(-dot_prod(1+1,2)*xr(1+1,1+2) +
     1       dot_prod(1+1,1)*xr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*xr(1,1+1) +
     1      dot_prod(1,2)*xr(1+1,1+2))) 

      fta_y(1) = -coef*(-dot_prod(1+1,2)*yr(1+1,1+2) +
     1      dot_prod(1+1,1)*yr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*yr(1,1+1) +
     1      dot_prod(1,2)*yr(1+1,1+2))) 

      fta_z(1) = -coef*(-dot_prod(1+1,2)*zr(1+1,1+2) +
     1      dot_prod(1+1,1)*zr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*zr(1,1+1) +
     1      dot_prod(1,2)*zr(1+1,1+2))) 

C }}}
C particle 2
C {{{
!     i = 2
      coef =DFAC(2+1)

      coef1 = DFAC(2)

      a1 =  -coef*(-dot_prod(2+1,2)*xr(2+1,2+2) +
     1      dot_prod(2+1,1)*xr(2+2,2+3) -
     1      (1.0/x_prod(2))*(dot_prod(2+1,2)*dot_prod(2,2) -
     1      dot_prod(2,3)*dot_prod(2+1,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     1      dot_prod(2,2)*xr(2+1,2+2))) 

      a2 = -coef1*(-dot_prod(2-1,2)*xr(2+1,2+2) +
     1      dot_prod(2,2)*xr(2,2+1) - dot_prod(2,2)*xr(2-1,2) -
     1      dot_prod(2,1)*xr(2+1,2+2) + 2.0*dot_prod(2-1,3)*xr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*xr(2-1,2) -
     1      dot_prod(2-1,1)*xr(2,2+1) - dot_prod(2-1,2)*xr(2,2+1) +
     1      dot_prod(2-1,2)*xr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     1      dot_prod(2,2)*xr(2+1,2+2))) 

      fta_x(2) = a1 + a2 

      a1 = -dot_prod(3,2)*yr(3,4) + dot_prod(3,1)*yr(4,5) 
      a1=a1-(dot_prod(3,2)*dot_prod(2,2) - dot_prod(2,3)*dot_prod(3,1))*
     &     (-dot_prod(3,1)*yr(2,3) + dot_prod(2,2)*yr(3,4))/x_prod(2)
      a1=-coef*a1

      a2 = -coef1*(-dot_prod(2-1,2)*yr(2+1,2+2) +
     1      dot_prod(2,2)*yr(2,2+1) - dot_prod(2,2)*yr(2-1,2) -
     1      dot_prod(2,1)*yr(2+1,2+2) + 2.0*dot_prod(2-1,3)*yr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*yr(2-1,2) -
     1      dot_prod(2-1,1)*yr(2,2+1) - dot_prod(2-1,2)*yr(2,2+1) +
     1      dot_prod(2-1,2)*yr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*yr(2,2+1) +
     1      dot_prod(2,2)*yr(2+1,2+2)))

      fta_y(2) = a1 + a2 
!     WRITE(MYUNIT,'(A,I6,5G15.5)') 'y i,a1,a2,coef,coef1,fta_y=',2,a1,a2,coef,coef1,fta_y(2)
      
      a1 = -coef*(-dot_prod(3,2)*zr(3,4) +
     1      dot_prod(3,1)*zr(4,5) -
     1      (1.0/x_prod(2))*(dot_prod(3,2)*dot_prod(2,2) -
     1      dot_prod(2,3)*dot_prod(3,1))*(-dot_prod(3,1)*zr(2,3) +
     1      dot_prod(2,2)*zr(3,4))) 

      a2 = -coef1*(-dot_prod(2-1,2)*zr(2+1,2+2) +
     1      dot_prod(2,2)*zr(2,2+1) - dot_prod(2,2)*zr(2-1,2) -
     1      dot_prod(2,1)*zr(2+1,2+2) + 2.0*dot_prod(2-1,3)*zr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*zr(2-1,2) -
     1      dot_prod(2-1,1)*zr(2,2+1) - dot_prod(2-1,2)*zr(2,2+1) +
     1      dot_prod(2-1,2)*zr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*zr(2,2+1) +
     1      dot_prod(2,2)*zr(2+1,2+2))) 

      fta_z(2) = a1 + a2 
!     WRITE(MYUNIT,'(A,I6,4G15.5)') 'z i,a1,a2,coef,coef1=',2,a1,a2,coef,coef1
C }}}
C particle 3
C {{{
!     i = 3
      coef=DFAC(3+1)

      coef1=DFAC(3)

      coef2=DFAC(3-1)

      a1 = -coef*(-dot_prod(3+1,2)*xr(3+1,3+2) +
     1      dot_prod(3+1,1)*xr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     1      dot_prod(3,2)*xr(3+1,3+2))) 

      a2 = -coef1*(-dot_prod(3-1,2)*xr(3+1,3+2) +
     1      dot_prod(3,2)*xr(3,3+1) - dot_prod(3,2)*xr(3-1,3) -
     1      dot_prod(3,1)*xr(3+1,3+2) + 2.0*dot_prod(3-1,3)*xr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*xr(3-1,3) -
     1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     1      dot_prod(3-1,2)*xr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     1      dot_prod(3,2)*xr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*xr(3,3+1) -
     1      dot_prod(3-2,2)*xr(3-1,3) + dot_prod(3-1,2)*xr(3-2,3-1) +
     1      dot_prod(3-1,1)*xr(3-2,3-1) - 2.0*dot_prod(3-2,3)*xr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*xr(3-1,3) -
     1      dot_prod(3-2,2)*xr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*xr(3-1,3) -
     1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     1      dot_prod(3-1,2)*xr(3-1,3))) 

      fta_x(3) = a1 + a2 + a3 
 
      a1 = -coef*(-dot_prod(3+1,2)*yr(3+1,3+2) +
     1      dot_prod(3+1,1)*yr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     1      dot_prod(3,2)*yr(3+1,3+2))) 
      
      a2 = -coef1*(-dot_prod(3-1,2)*yr(3+1,3+2) +
     1      dot_prod(3,2)*yr(3,3+1) - dot_prod(3,2)*yr(3-1,3) -
     1      dot_prod(3,1)*yr(3+1,3+2) + 2.0*dot_prod(3-1,3)*yr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*yr(3-1,3) -
     1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     1      dot_prod(3-1,2)*yr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     1      dot_prod(3,2)*yr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*yr(3,3+1) -
     1      dot_prod(3-2,2)*yr(3-1,3) + dot_prod(3-1,2)*yr(3-2,3-1) +
     1      dot_prod(3-1,1)*yr(3-2,3-1) - 2.0*dot_prod(3-2,3)*yr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*yr(3-1,3) -
     1      dot_prod(3-2,2)*yr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*yr(3-1,3) -
     1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     1      dot_prod(3-1,2)*yr(3-1,3)))

      fta_y(3) = a1 + a2 + a3 
 
      a1 = -coef*(-dot_prod(3+1,2)*zr(3+1,3+2) +
     1      dot_prod(3+1,1)*zr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     1      dot_prod(3,2)*zr(3+1,3+2))) 

      a2 =  -coef1*(-dot_prod(3-1,2)*zr(3+1,3+2) +
     1      dot_prod(3,2)*zr(3,3+1) - dot_prod(3,2)*zr(3-1,3) -
     1      dot_prod(3,1)*zr(3+1,3+2) + 2.0*dot_prod(3-1,3)*zr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*zr(3-1,3) -
     1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     1      dot_prod(3-1,2)*zr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     1      dot_prod(3,2)*zr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*zr(3,3+1) -
     1      dot_prod(3-2,2)*zr(3-1,3) + dot_prod(3-1,2)*zr(3-2,3-1) +
     1      dot_prod(3-1,1)*zr(3-2,3-1) - 2.0*dot_prod(3-2,3)*zr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*zr(3-1,3) -
     1      dot_prod(3-2,2)*zr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*zr(3-1,3) -
     1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     1      dot_prod(3-1,2)*zr(3-1,3))) 

      fta_z(3) = a1 + a2 + a3 
C }}}
C particles 4 to n-3
C {{{
      do i = 4, n-3

         coef=DFAC(i+1)

         coef1=DFAC(i)

         coef2=DFAC(i-1)

         coef3=DFAC(i-2)

         a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1      dot_prod(i+1,1)*xr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1      dot_prod(i,2)*xr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1      dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1      dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1      dot_prod(i-1,2)*xr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1      dot_prod(i,2)*xr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1      dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1      dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1      dot_prod(i-2,2)*xr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1      dot_prod(i-1,2)*xr(i-1,i))) 

         a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1      dot_prod(i-2,1)*xr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1      dot_prod(i-2,2)*xr(i-2,i-1))) 

         fta_x(i) = a1 + a2 + a3 + a4 

         a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1      dot_prod(i+1,1)*yr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1      dot_prod(i,2)*yr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1      dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1      dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1      dot_prod(i-1,2)*yr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1      dot_prod(i,2)*yr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1      dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1      dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1      dot_prod(i-2,2)*yr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1      dot_prod(i-1,2)*yr(i-1,i))) 

         a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1      dot_prod(i-2,1)*yr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1      dot_prod(i-2,2)*yr(i-2,i-1))) 

         fta_y(i) = a1 + a2 + a3 + a4 

         a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1      dot_prod(i+1,1)*zr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1      dot_prod(i,2)*zr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1      dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1      dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1      dot_prod(i-1,2)*zr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1      dot_prod(i,2)*zr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1      dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1      dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1      dot_prod(i-2,2)*zr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1      dot_prod(i-1,2)*zr(i-1,i))) 
      
         a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1      dot_prod(i-2,1)*zr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1      dot_prod(i-2,2)*zr(i-2,i-1))) 

         fta_z(i) = a1 + a2 + a3 + a4 

      enddo
C }}}
C particle n-2
C {{{
!     i = n-2
      coef1=DFAC(n-2)

      coef2=DFAC(n-2-1)

      coef3=DFAC(n-2-2)

      a1 = -coef1*(-dot_prod(n-2-1,2)*xr(n-2+1,n-2+2) + 
     1      dot_prod(n-2,2)*xr(n-2,n-2+1) - dot_prod(n-2,2)*xr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*xr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*xr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*xr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*xr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*xr(n-2+1,n-2+2)))


      a2 = -coef2*(dot_prod(n-2-2,2)*xr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*xr(n-2-1,n-2) + dot_prod(n-2-1,2)*xr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*xr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*xr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*xr(n-2-1,n-2))) 
      
      a3 = -coef3*(dot_prod(n-2-3,2)*xr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*xr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1))) 

      fta_x(n-2) = a1 + a2 + a3 

      a1 =  -coef1*(-dot_prod(n-2-1,2)*yr(n-2+1,n-2+2) +  
     1      dot_prod(n-2,2)*yr(n-2,n-2+1) - dot_prod(n-2,2)*yr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*yr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*yr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*yr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*yr(n-2+1,n-2+2))) 

      a2 =  -coef2*(dot_prod(n-2-2,2)*yr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*yr(n-2-1,n-2) + dot_prod(n-2-1,2)*yr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*yr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*yr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)))

      a3 = -coef3*(dot_prod(n-2-3,2)*yr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*yr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1))) 

      fta_y(n-2) = a1 + a2 + a3 
 
      a1 = -coef1*(-dot_prod(n-2-1,2)*zr(n-2+1,n-2+2) +  
     1      dot_prod(n-2,2)*zr(n-2,n-2+1) - dot_prod(n-2,2)*zr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*zr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*zr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*zr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*zr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*zr(n-2+1,n-2+2))) 

      a2 = -coef2*(dot_prod(n-2-2,2)*zr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*zr(n-2-1,n-2) + dot_prod(n-2-1,2)*zr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*zr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*zr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*zr(n-2-1,n-2))) 

      a3 = -coef3*(dot_prod(n-2-3,2)*zr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*zr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1))) 

      fta_z(n-2) = a1 + a2 + a3 
C }}}
C particle n-1
C {{{
!     i = n-1
      coef2=DFAC(n-1-1)

      coef3=DFAC(n-1-2)

      a1 = -coef2*(dot_prod(n-1-2,2)*xr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*xr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*xr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*xr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*xr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-1,1)*xr(n-1,n-1+1) - dot_prod(n-1-1,2)*xr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*xr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*xr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*xr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1))) 

      fta_x(n-1) = a1 + a2  

      a1 = -coef2*(dot_prod(n-1-2,2)*yr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*yr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*yr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*yr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*yr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*yr(n-1-1,n-1) -
     1  dot_prod(n-1-1,1)*yr(n-1,n-1+1) - dot_prod(n-1-1,2)*yr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*yr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*yr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*yr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1))) 

      fta_y(n-1) = a1 + a2  

      a1 = -coef2*(dot_prod(n-1-2,2)*zr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*zr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*zr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*zr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*zr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-1,1)*zr(n-1,n-1+1) - dot_prod(n-1-1,2)*zr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*zr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*zr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*zr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1))) 

      fta_z(n-1) = a1 + a2 
C }}} 
C particle n
C {{{
!     i = n
      coef3=DFAC(n-2)

      fta_x(n) = -coef3*(dot_prod(n-3,2)*xr(n-2,n-1) 
     1      - dot_prod(n-2,1)*xr(n-3,n-2) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-1,n) -
     1      dot_prod(n-2,2)*xr(n-2,n-1))) 

      fta_y(n) = -coef3*(dot_prod(n-3,2)*yr(n-2,n-1) - 
     1      dot_prod(n-2,1)*yr(n-3,n-2) 
     1      - (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) 
     1      - dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-1,n) -
     1      dot_prod(n-2,2)*yr(n-2,n-1))) 

      fta_z(n) = -coef3*(dot_prod(n-3,2)*zr(n-2,n-1) - 
     1      dot_prod(n-2,1)*zr(n-3,n-2) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-1,n) -
     1      dot_prod(n-2,2)*zr(n-2,n-1))) 

      if (RADIUS_CONTAINER.GE.0.D0) THEN

c     container contribution

         do i = 1, n

            radi=sqrt(x(i)**2+y(i)**2+z(i)**2)
            dum10p=1.D0/(radi+RADIUS_CONTAINER)**10
            dum10m=1.D0/(radi-RADIUS_CONTAINER)**10
            dum4p=1.D0/(radi+RADIUS_CONTAINER)**4
            dum4m=1.D0/(radi-RADIUS_CONTAINER)**4
            
            dummy1=(dum10m-dum10p)/5.D0-hydro_bln(i)*(dum4m-dum4p)/2.D0
            dummy2=-dum10m/(radi-RADIUS_CONTAINER)
     &           +dum10p/(radi+RADIUS_CONTAINER)
            dummy3=-dum4m/(radi-RADIUS_CONTAINER)
     &           +dum4p/(radi+RADIUS_CONTAINER)
            
            dummy4=-dummy1/radi+2.D0*dummy2-2.D0*hydro_BLN(i)*dummy3
            
            dummy4=-4.D0*pi*RADIUS_CONTAINER*dummy4/radi**2
            
            fc_x(i)=dummy4*x(i)
            fc_y(i)=dummy4*y(i)
            fc_z(i)=dummy4*z(i)
            
         enddo

      endif
C }}}
C Total up the gradients
C {{{
      do i = 1, n
!        IF (I.EQ.2) THEN
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbx,fbx,fbax,ftax=',i,fnb_x(i),fb_x(i),fba_x(i),fta_x(i)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnby,fby,fbay,ftay=',i,fnb_y(i),fb_y(i),fba_y(i),fta_y(i)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbz,fbz,fbaz,ftaz=',i,fnb_z(i),fb_z(i),fba_z(i),fta_z(i)
!        ENDIF
         fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i) + fc_x(i) 
         fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i) + fc_y(i) 
         fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i) + fc_z(i) 
      enddo

      do i = 1, n
         j = (i-1)*3
         fq(j+1) = -fx(i)
         fq(j+2) = -fy(i)
         fq(j+3) = -fz(i)
      enddo
C }}}
      return
      end
C }}}
C
C>  Fill the parameter arrays
C
      subroutine param_arrayCHAPERONIN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN
     $     ,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,LJREPBB, LJATTBB, LJREPLL,
     $     LJATTLL, LJREPNN, LJATTNN,LJREPBL, LJATTBL, LJREPBN, LJATTBN,
     $     LJREPLN, LJATTLN,HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN,
     $     ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, HYDROPHOBIC,
     $     HYDRO_BLN, n)
C {{{
C Declarations {{{
      implicit NONE
      INTEGER N
      INTEGER ntype(n), J1, i, j, icount
      DOUBLE PRECISION LJREP_BLN(n,n), LJATT_BLN(n,n), A_BLN(n), B_BLN(n
     $     ), C_BLN(N), D_BLN(N),LJREPBB, LJATTBB, LJREPLL, LJATTLL,
     $     LJREPNN, LJATTNN,LJREPBL, LJATTBL, LJREPBN, LJATTBN, LJREPLN,
     $     LJATTLN,HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN,
     $     EDBLN, TABLN, TBBLN, TCBLN, TDBLN, HYDROPHOBIC, HYDRO_BLN(n)
      CHARACTER(LEN=1) BEADLETTER(N), BLNSSTRUCT(N)
      DOUBLE PRECISION PI
      PARAMETER(PI=3.141592653589793D0)

C }}}
C Amino Acid types {{{
C B=1 L=2 N=3
C
      DO J1=1,N
         IF (BEADLETTER(J1).EQ.'B') THEN
            ntype(J1)=1
            hydro_bln(j1)=hydrophobic
         ELSEIF (BEADLETTER(J1).EQ.'L') THEN
            ntype(J1)=2
            hydro_bln(j1)=0.D0
         ELSEIF (BEADLETTER(J1).EQ.'N') THEN
            ntype(J1)=3
            hydro_bln(j1)=0.D0
         ELSE
            PRINT '(A,A1)','ERROR in param_arrayBLN, unrecognised bead type: ',BEADLETTER(J1)
            STOP
         ENDIF
      ENDDO

C }}}

C Parameters for the dihedral angle potential 
C The end bonds have no dihedral term, so the total number of terms
C is N-3 (e.g. 4 atoms, 1 dihedral). Non-zero terms for
C 2 to 3, 3 to 4, 4 to 5, ... , N-3 to N-2, N-2 to N-1. The
C H, E and T parameters are defined for the first bead of each edge,
C i.e. for 2, 3, 4, ..., N-2.
C {{{

      A_BLN(1:N)=0.0D0
      B_BLN(1:N)=0.0D0
      C_BLN(1:N)=0.0D0
      D_BLN(1:N)=0.0D0
      DO I=1,N-3
         IF (BLNSSTRUCT(I).EQ.'H') THEN
            A_BLN(I+1)=HABLN
            B_BLN(I+1)=HBBLN
            C_BLN(I+1)=HCBLN
            D_BLN(I+1)=HDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'E') THEN
            A_BLN(I+1)=EABLN
            B_BLN(I+1)=EBBLN
            C_BLN(I+1)=ECBLN
            D_BLN(I+1)=EDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'T') THEN
            A_BLN(I+1)=TABLN
            B_BLN(I+1)=TBBLN
            C_BLN(I+1)=TCBLN
            D_BLN(I+1)=TDBLN
         ELSE
            PRINT '(A,A1)','ERROR in param_arrayCHAPERONIN, unrecognised SS type: ',BLNSSTRUCT(J1)
            STOP
         ENDIF
C        PRINT '(A,I6,A,A1,A,4F12.4)','I+1=',I+1,' symbol=',BLNSSTRUCT(I),' A,B,C,D=',A_BLN(I+1),B_BLN(I+1),C_BLN(I+1),D_BLN(I+1)
      ENDDO
C }}}
C  Parameters for the L-J interaction between non-bonded particles {{{

      do i = 1, n-1
         do j = i+1, n
            if (ntype(i) .eq. 1) then
               if (ntype(j) .eq. 1) then  ! B-B
                  LJREP_BLN(i,j) = LJREPBB
                  LJATT_BLN(i,j) = LJATTBB
                  LJREP_BLN(j,i) = LJREPBB
                  LJATT_BLN(j,i) = LJATTBB
               elseif (ntype(j) .eq. 2) then ! B-L
                  LJREP_BLN(i,j) = LJREPBL
                  LJATT_BLN(i,j) = LJATTBL
                  LJREP_BLN(j,i) = LJREPBL
                  LJATT_BLN(j,i) = LJATTBL
               elseif (ntype(j) .eq. 3) then ! B-N
                  LJREP_BLN(i,j) = LJREPBN
                  LJATT_BLN(i,j) = LJATTBN
                  LJREP_BLN(j,i) = LJREPBN
                  LJATT_BLN(j,i) = LJATTBN
               endif
            elseif (ntype(i) .eq. 2) then
               if (ntype(j) .eq. 1) then  ! L-B
                  LJREP_BLN(i,j) = LJREPBL
                  LJATT_BLN(i,j) = LJATTBL
                  LJREP_BLN(j,i) = LJREPBL
                  LJATT_BLN(j,i) = LJATTBL
               elseif (ntype(j) .eq. 2) then ! L-L
                  LJREP_BLN(i,j) = LJREPLL
                  LJATT_BLN(i,j) = LJATTLL
                  LJREP_BLN(j,i) = LJREPLL
                  LJATT_BLN(j,i) = LJATTLL
               elseif (ntype(j) .eq. 3) then ! L-N
                  LJREP_BLN(i,j) = LJREPLN
                  LJATT_BLN(i,j) = LJATTLN
                  LJREP_BLN(j,i) = LJREPLN
                  LJATT_BLN(j,i) = LJATTLN
               endif
            elseif (ntype(i) .eq. 3) then
               if (ntype(j) .eq. 1) then  ! N-B
                  LJREP_BLN(i,j) = LJREPBN
                  LJATT_BLN(i,j) = LJATTBN
                  LJREP_BLN(j,i) = LJREPBN
                  LJATT_BLN(j,i) = LJATTBN
               elseif (ntype(j) .eq. 2) then ! N-L
                  LJREP_BLN(i,j) = LJREPLN
                  LJATT_BLN(i,j) = LJATTLN
                  LJREP_BLN(j,i) = LJREPLN
                  LJATT_BLN(j,i) = LJATTLN
               elseif (ntype(j) .eq. 3) then ! N-N
                  LJREP_BLN(i,j) = LJREPNN
                  LJATT_BLN(i,j) = LJATTNN
                  LJREP_BLN(j,i) = LJREPNN
                  LJATT_BLN(j,i) = LJATTNN
               endif
            endif
         enddo
      enddo
C }}}

      return
      end
C }}}
