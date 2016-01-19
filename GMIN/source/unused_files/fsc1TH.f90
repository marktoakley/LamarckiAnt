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
subroutine fsc1TH ( kEX, coso, xv, fsc, dfscCJ )

USE nrtype
use stats
IMPLICIT NONE
include 'NP.H'

integer, parameter                 :: lmax = NP-1
real(sp)                           :: kEX, coso, cosi, intTH
real(sp), dimension(NP)            :: xv 
integer(i4b)                       :: il
complex(dpc)                       :: fsc
complex(dpc),dimension(0:lmax)     :: dfscCJ
real(sp), dimension(0:lmax)        :: LegendreP, deltaTH

stats_fsc = stats_fsc + 1.

do il = 0, lmax
   deltaTH (il) = xv (il+1)
enddo

call PLegendre ( coso, LegendreP, lmax )

fsc = (0.0_dp, 0.0_dp)
do il = 0, lmax
   fsc = fsc + LegendreP(il) * (2*il+1) * ( Exp ( 2.*CI*deltaTH(il) ) - 1. )
enddo
fsc = fsc / ( 2.*CI*kEX)

dfscCJ = (0.0_dp,0.0_dp)
do il = 0, lmax
   dfscCJ(il) = &
            LegendreP(il) * (2*il+1) *  Exp ( - 2.*CI*xv(il+1) ) / kEX
enddo   ! il

return
END subroutine fsc1TH
