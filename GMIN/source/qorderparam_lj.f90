!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! 
!---======================================---
      SUBROUTINE QORDER_LJ(Q,Q4,Q6)
      USE COMMONS, ONLY: NATOMS

      IMPLICIT NONE

      REAL(8) Q(3,NATOMS), Q4, Q6
        INTEGER J, K, NB, I,m
        DOUBLE PRECISION DISTBOND
        PARAMETER (DISTBOND=1.3909)
        DOUBLE PRECISION DX,DY,DZ,DIST,phi,costheta,coef2,arg,pi
        COMPLEX Y2,Q4bar(0:4),Q6bar(0:6)


	NB=0
	pi=dacos(-1d0)
	do m=0,4
	   Q4bar(m)=(0,0d0)
	enddo
	do m=0,6
	   Q6bar(m)=(0,0d0)
	enddo
	do J=1,NATOMS-1
	   do K=J+1,NATOMS
	      DX=Q(1,J)-Q(1,K)
	      DY=Q(2,J)-Q(2,K)
	      DZ=Q(3,J)-Q(3,K)

	      DIST=DSQRT(DX*DX+DY*DY+DZ*DZ)
	      IF (DIST.lt.DISTBOND) THEN
		 NB=NB+1
		 costheta =DZ/DIST
		 phi = datan(DY/DX)
		 if(DX.lt.0) phi=phi+pi
		 do m=0,4
		    Q4bar(m)=Q4bar(m)+Y2(4,m,costheta,phi)
		 enddo
		 do m=0,6
		    Q6bar(m)=Q6bar(m)+Y2(6,m,costheta,phi)
		 enddo
	      ENDIF
	   enddo
	enddo
	Q4=0
	do m=0,4
	   Q4=Q4+coef2(4,m)*Q4bar(m)*conjg(Q4bar(m))
	enddo
	Q4=dsqrt(Q4)/(3*Nb)
	Q6=0
	do m=0,6
	   Q6=Q6+coef2(6,m)*Q6bar(m)*conjg(Q6bar(m))
	enddo
	Q6=dsqrt(Q6/13.D0)/Nb
	return
	end
	
	DOUBLE PRECISION function coef2(l,m)
	integer l,m,k
	if(m.eq.0) then
	   coef2=2*l+1d0
	else
	   coef2=4*l+2.
	   do k=l-m+1,l+m
	      coef2=coef2/k
	   enddo
	endif
	return
	end
	
	COMPLEX function Y2(l,m,costheta,phi)
!Computes the the spherical harmonic Y2(l.m)*sqrt(4*pi)  
	implicit none
	INTEGER l,m
	DOUBLE PRECISION plgndr2,costheta,phi
	if(m.lt.0) stop 'm<0 in Y2'
	if(m.eq.0) then
	   Y2=plgndr2(l,m,costheta)
	else 
	   Y2=plgndr2(l,m,costheta)*exp(cmplx(0d0,m*phi))
	endif
	return 
	END
	
	DOUBLE PRECISION FUNCTION plgndr2(l,m,x) 
	implicit none
	INTEGER l,m 
	DOUBLE PRECISION x
!Computes the associated Legendre polynomial Pml (x).
	INTEGER i,ll 
	DOUBLE PRECISION fact,pll,pmm,pmmp1,somx2 
	if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) &
        & stop 'bad arguments in plgndr2'
	pmm=1.  
	if(m.gt.0) then 
	   somx2=sqrt((1.-x)*(1.+x)) 
	   fact=1. 
	   do i=1,m 
	      pmm=-pmm*fact*somx2 
	      fact=fact+2. 
	   enddo 
	endif 
	if(l.eq.m) then 
	   plgndr2=pmm 
	else 
	   pmmp1=x*(2*m+1)*pmm 
	   if(l.eq.m+1) then 
	      plgndr2=pmmp1 
	   else 
	      do ll=m+2,l
		 pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m) 
		 pmm=pmmp1 
		 pmmp1=pll 
	      enddo  
	      plgndr2=pll 
	   endif 
	endif 
	return 
	END



