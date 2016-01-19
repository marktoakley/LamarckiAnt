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
C*************************************************************************
C
C  Here we calculate the EAMLJ potential and gradient
C  (Baskes, PRL 27, 2592 (1999)). It has three parameters A0, beta and Z0. 
C                                        
C*************************************************************************

      SUBROUTINE EAMLJ(XALL,V,energy,GRADT)
      USE commons
      IMPLICIT NONE 
      INTEGER NSIZE, i, J
      DOUBLE PRECISION XALL(3*NATOMS), energy, x(NATOMS), y(NATOMS), z(NATOMS),
     +                 A0, beta, Z0, rhol(NATOMS), r(NATOMS,NATOMS), ww, v1, v2, 
     +                 V(3*NATOMS), fx(NATOMS), fy(NATOMS), fz(NATOMS), br, func, 
     +                 rij, r6, wexp, wLJ, wtot 
      LOGICAL GRADT
      common /EAMLJCOMM/ A0,beta,Z0

      NSIZE=NATOMS

      do i=1,nsize
        x(i)=xall(3*(i-1)+1)
        y(i)=xall(3*(i-1)+2)
        z(i)=xall(3*(i-1)+3)
      enddo

c       CALCULATE ENERGY

      do i=1,Nsize
         do j=i+1,Nsize
            ww=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
            ww=dsqrt(ww)
            r(i,j)=ww
            r(j,i)=r(i,j)
         enddo
      enddo
      
      v1=0.D0
      v2=0.D0

      do i=1,nsize
         fx(i)=0.D0
         fy(i)=0.D0
         fz(i)=0.D0
         rhol(i)=0.D0
         do j=1,Nsize
            if (i.ne.j) then
             rhol(i)=rhol(i)+dexp(-beta*(r(i,j)-1.D0))
            endif
         enddo
         rhol(i)=rhol(i)/Z0
         v1=v1+func(rhol(i))

         do j=i+1,Nsize
            r6=1.D0/r(i,j)**6
            v2=v2+r6*(r6-2.D0)
            ww=dexp(-beta*(r(i,j)-1.D0))
            v2=v2-2.D0*func(ww)/Z0
         enddo

      enddo

      energy=v1+v2

C        print*,energy
 
        IF (.NOT.GRADT) RETURN

      do i=1,Nsize
         do j=1,Nsize
            if (j.ne.i) then
             rij=r(j,i)

             br=beta*(rij-1.D0)
             wexp=dexp(-br)

             wtot=-0.5D0*(dlog(rhol(j)*rhol(i)))
             wtot=wtot-br
             wtot=A0*beta*wexp*wtot
             wlj=12.D0/(rij**7)
             wlj=wlj*(1.D0-1.D0/rij**6)
             wtot=wtot+wlj
C             wtot=-wtot
             fx(i)=fx(i)+(x(i)-x(j))*wtot/rij
             fy(i)=fy(i)+(y(i)-y(j))*wtot/rij
             fz(i)=fz(i)+(z(i)-z(j))*wtot/rij
            endif
         enddo
      enddo

        do i=1,nsize
c          write(*,'(3F20.10)') fx(i), fy(i), fz(i)
          v(3*(i-1)+1)=fx(i)
          v(3*(i-1)+2)=fy(i)
          v(3*(i-1)+3)=fz(i)
        enddo 

      return
      end

c__________________________________________________________________________

      function func(x)
        implicit none
        double precision a0, beta, z0, x, func
      common/EAMLJCOMM/A0,beta,Z0
      func=A0*Z0*x
      func=func*(dlog(x)-1.D0)
      func=0.5D0*func
      return
      end
