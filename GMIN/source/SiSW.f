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
C**********************************************************************
C
C   Subroutines for calculating the force (minus the gradient) and the
C   energy of a Stillinger-Weber system.
C
C   as mentioned in Barkema+Mousseau (PRL 77, 4358 (1996)) the 3-body
C   strength of the SW potential is too weak to lead to good amorphous
C   structures. It is strengthed here by a factor 1.5 (see LAMBDA).
C
C   The subroutine presupposes periodic boundary coundition in a box
C   centered around -0.5*box and +0.5*box (the same for the y and z
C   direction)
C
C   In input:  pos: positions in angstroems (as a 3*NATOMS vector)
C              nei: list of neighbours 
C              numnei: number of neighbours for atom i
C              box,boy,boz: size of the box along x,y and z
C                           directions. 
C
C  As output:  the total energy (in both calcforce and calcenergy
C              force:  a 3*NATOMS vector of the force on the atoms.
C
C   Written by Normand Mousseau 
C              13 Oct. 1997
C
C********************************************************************

      subroutine SISW(pos,GRAD,totenergy,GTEST)
      USE commons
      implicit none
      LOGICAL GTEST
      integer MAXNEI,N
      double precision SIGMA,A,ALPHA,BETA,EPSILON,SWGAMMA,LAMBDA,ONE_THIRD,A_EPS,P,RCUT,
     1                 box, boy, boz
      parameter (MAXNEI=36)
      parameter (SIGMA=2.0951,A=7.049556277,ALPHA=1.8)
      parameter (BETA=0.6022245584,EPSILON=2.16823,P=4)
      parameter (SWGAMMA=1.2,LAMBDA=(21.0*1.5),RCUT= (ALPHA*SIGMA))
      parameter (ONE_THIRD=(1.0/3.0),A_EPS=(A*EPSILON))
      integer numnei(NATOMS),nei(NATOMS,MAXNEI)
      double precision pos(3*NATOMS),GRAD(3*NATOMS),totenergy
 
      integer  i,j,k,indj,indk
      double precision  xi,yi,zi,xij,yij,zij,rij,rij2
      double precision  rhoij,cos_x_ij,cos_y_ij,cos_z_ij,
     +     one_o_a_ij,expo,gam_o_a_ij,exp_gam_ij,gam_o_a_ik,
     +     cos_jik,cos_x_ik,cos_y_ik,cos_z_ik,cos_p_1o3,
     +     xik,yik,zik,rik,rik2,rhoik,ffx,ffy,ffz,
     +     r_to_minusp,one_o_a_ik,term_ik,fact,
     +     term1,term2,fact3_ij,term_ij,dhdcos_ik,
     +     dhdcos_ij, one_o_a2,twobody,threebody
 
      double precision dcosdxj,dcosdyj,dcosdzj,dcosdxk,dcosdyk,dcosdzk
      double precision invsig,invrij,invrik
      double precision rcut2
      double precision x(NATOMS),y(NATOMS),z(NATOMS),
     +                 fx(NATOMS),fy(NATOMS),fz(NATOMS)

      N=3*NATOMS
      rcut2=RCUT*RCUT

      box=2.0D0*boxlx
      boz=2.0D0*boxlx
      boy=2.0D0*boxlx
      do i=1,NATOMS
        x(i)= pos(3*(I-1)+1)
        y(i)= pos(3*(I-1)+2)
        z(i)= pos(3*(I-1)+3)
      end do

      call neighbour(NATOMS,N,pos,numnei,nei,box,boy,boz,x,y,z)

      IF (GTEST) THEN
         do i=1,NATOMS
           fx(i)= 0.0d0
           fy(i)= 0.0d0
           fz(i)= 0.0d0
         end do
      ENDIF

      twobody=0.0d0
      threebody=0.0d0
      totenergy=0.0d0
      invsig=1.0/SIGMA

      do i=1,NATOMS
        xi=x(i)
        yi=y(i)
        zi=z(i)
        do indj=1,numnei(i)
          j=nei(i,indj)

C  Pair interactions of i and j  with periodic boundary counditions
C  Distance 
          xij=x(j)-xi
          yij=y(j)-yi
          zij=z(j)-zi
          xij=xij-boxlx*anint(xij/boxlx)
          yij=yij-boxly*anint(yij/boxly)
          zij=zij-boxlz*anint(zij/boxlz)
          rij2=xij*xij+yij*yij+zij*zij

C   Check whether the distance is too large 
          if (rij2.lt.rcut2) then
            rij=dsqrt(rij2)
            invrij = 1.0/rij
            rhoij=rij*invsig
            cos_x_ij=xij*invrij
            cos_y_ij=yij*invrij
            cos_z_ij=zij*invrij
    
C  Some useful quantities 
            one_o_a_ij =1.0d0/(rhoij-ALPHA)
            expo=dexp(one_o_a_ij)
            gam_o_a_ij=SWGAMMA*one_o_a_ij
            exp_gam_ij=dexp(gam_o_a_ij)
            r_to_minusp=rhoij**(-1.0*P)
    
C  Two body energy and force 
            term1=A_EPS*(BETA*r_to_minusp-1.0)*expo
            one_o_a2 = one_o_a_ij*one_o_a_ij
            term2=(one_o_a2*term1
     +            +A_EPS*P*BETA*r_to_minusp*expo/rhoij)*invsig
    
            twobody=twobody+0.5*term1
            totenergy=totenergy+0.5*term1
    
            IF (GTEST) THEN
               fx(i)=fx(i)-term2*cos_x_ij
               fy(i)=fy(i)-term2*cos_y_ij
               fz(i)=fz(i)-term2*cos_z_ij
            ENDIF

C  Prepare for the three body term  
            fact3_ij=gam_o_a_ij*one_o_a_ij*invsig
    
            do indk=indj+1,numnei(i)
              k=nei(i,indk)
C  Distance  and PRB
              xik=x(k)-xi
              yik=y(k)-yi
              zik=z(k)-zi
              xik=xik-boxlx*anint(xik/boxlx)
              yik=yik-boxly*anint(yik/boxly)
              zik=zik-boxlz*anint(zik/boxlz)
              rik2=xik*xik+yik*yik+zik*zik
    
              if (rik2.le.rcut2) then
                rik=dsqrt(rik2)
                invrik=1.0/rik
                rhoik=rik*invsig
                cos_x_ik=xik*invrik
                cos_y_ik=yik*invrik
                cos_z_ik=zik*invrik
    
C          Some useful quantities 
                one_o_a_ik=1.0/(rhoik-ALPHA)
                gam_o_a_ik=SWGAMMA*one_o_a_ik
                fact     =EPSILON*LAMBDA*exp(gam_o_a_ik)*exp_gam_ij
                cos_jik=cos_x_ij*cos_x_ik+cos_y_ij*cos_y_ik+ 
     +                  cos_z_ij*cos_z_ik
                cos_p_1o3=cos_jik+ONE_THIRD
    
C          Energy (added only to central atom)  
                totenergy=totenergy+fact*cos_p_1o3*cos_p_1o3
                threebody=threebody+fact*cos_p_1o3*cos_p_1o3
    
    
C          Force 
                IF (GTEST) THEN
                   term_ij=fact*fact3_ij*cos_p_1o3*cos_p_1o3
                   dhdcos_ij=2*fact*cos_p_1o3
                   term_ik=fact*gam_o_a_ik*one_o_a_ik*cos_p_1o3*cos_p_1o3/SIGMA
                   dhdcos_ik=2*fact*cos_p_1o3
    
                   dcosdxj=(cos_x_ik-cos_jik*cos_x_ij)*invrij
                   dcosdyj=(cos_y_ik-cos_jik*cos_y_ij)*invrij
                   dcosdzj=(cos_z_ik-cos_jik*cos_z_ij)*invrij
    
                   dcosdxk=(cos_x_ij-cos_jik*cos_x_ik)*invrik
                   dcosdyk=(cos_y_ij-cos_jik*cos_y_ik)*invrik
                   dcosdzk=(cos_z_ij-cos_jik*cos_z_ik)*invrik
    
                   ffx=term_ij*cos_x_ij-dhdcos_ij*dcosdxj
                   ffy=term_ij*cos_y_ij-dhdcos_ij*dcosdyj
                   ffz=term_ij*cos_z_ij-dhdcos_ij*dcosdzj
                   fx(j)=fx(j)+ffx
                   fy(j)=fy(j)+ffy
                   fz(j)=fz(j)+ffz
                   fx(i)=fx(i)-ffx
                   fy(i)=fy(i)-ffy
                   fz(i)=fz(i)-ffz
                   ffx=term_ik*cos_x_ik-dhdcos_ik*dcosdxk
                   ffy=term_ik*cos_y_ik-dhdcos_ik*dcosdyk
                   ffz=term_ik*cos_z_ik-dhdcos_ik*dcosdzk
                   fx(k)=fx(k)+ffx
                   fy(k)=fy(k)+ffy
                   fz(k)=fz(k)+ffz
                   fx(i)=fx(i)-ffx
                   fy(i)=fy(i)-ffy
                   fz(i)=fz(i)-ffz
                 ENDIF
              endif
            enddo
          endif
        enddo
      enddo

      IF (GTEST) THEN
         do i=1,NATOMS
C          force(i)=fx(i)
C          force(i+NATOMS)=fy(i)
C          force(i+2*NATOMS)=fz(i)

           GRAD(3*(I-1)+1)=-fx(i)
           GRAD(3*(I-1)+2)=-fy(i)
           GRAD(3*(I-1)+3)=-fz(i)

         end do
      ENDIF

C     print *,'twobody: ',twobody,'  threebody',threebody
      return
      end
    
      subroutine neighbour(NATOMS,N,pos,numnei,nei,box,boy,boz,x,y,z)
      implicit none
      integer NATOMS,MAXNEI,N
      double precision SIGMA,ALPHA,RCUT
      parameter (MAXNEI=36)
      parameter (SIGMA=2.0951,ALPHA=1.8,RCUT= (ALPHA*SIGMA))
      integer numnei(NATOMS),nei(NATOMS,MAXNEI)
      double precision pos(N),box,boy,boz
      integer i,j
      double precision xi,yi,zi,xij,yij,zij,rcut2,rij2
      double precision x(NATOMS),y(NATOMS),z(NATOMS)

C  Note that the neighbour list uses a slightly enlarged cut-off
C  so that one does not need to recompute the neighbour list at
C  each iteration

      rcut2=1.4*RCUT*RCUT

      do i=1,NATOMS
        numnei(i)=0
        do j=1,MAXNEI
          nei(i,j)=0
        end do
      end do 

      do i=1,NATOMS
        xi=x(i)
        yi=y(i)
        zi=z(i)
        do j=i+1,NATOMS

C  Pair interactions of i and j  with periodic boundary counditions
C  Distance
          xij=x(j)-xi
          yij=y(j)-yi
          zij=z(j)-zi
          xij=xij-box*anint(xij/box)
          yij=yij-boy*anint(yij/boy)
          zij=zij-boz*anint(zij/boz)
          rij2=xij*xij+yij*yij+zij*zij

C   Check whether the distance is too large
          if (rij2.lt.rcut2) then
            numnei(i)=numnei(i)+1
            numnei(j)=numnei(j)+1
            if ((numnei(i).gt.MAXNEI).OR.(numnei(j).gt.MAXNEI))  then
               print *, 'WARNING - Too many neighbours'
            ELSE
               nei(i,numnei(i))=j
               nei(j,numnei(j))=i
            endif
          endif
        end do
      end do

      return
      end
