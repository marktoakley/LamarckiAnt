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
MODULE distance 
      IMPLICIT NONE
      SAVE

!   Some of the distances etc...that I will need in common blocks for
!   deivatives program (derivs.1 etc...)

!INTEGER N


DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG   
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: deriv1st  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R, S2, LAMBDA, REPULSE, BONDAGE, KAPPANEW  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Horig  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SVtype  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DIRCOS 

!      INTEGER COUNTER
!      DOUBLE PRECISION Sorig(4*MXATMS,4*MXATMS),UBOND,UREP 
!   Should have separate common blocks for integers and double 
!   precision thingamies. 
!      COMMON /MISC/ R,DIRCOS,S2,LAMBDA,SVtype,Horig,DIAG,REPULSE,BONDAGE,&
!     &              KAPPANEW


END MODULE distance
