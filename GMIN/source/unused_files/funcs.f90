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
FUNCTION funcs ( xv )
! Una funcion escalar (argumento vectorial), 
! f ( \vec x )

! Su gradiente (componente a componente), 
! grad_{j} = \sum_{i} \frac{\partial f_{i}}{\partial x_{j}}.

USE nrtype
IMPLICIT NONE
include 'NP.H'
REAL(SP), DIMENSION(:), INTENT(IN) :: xv
REAL(SP)                           :: funcs
complex(dpc)                       :: fsc
complex(dpc),dimension(NP)         :: dfscCJ
real(sp), parameter                :: favg = 1.e6
REAL(SP)                           :: coso, cosi 
integer(i4b)                       :: io, ii, ix

real(sp)                           :: kEX
common /ENE/ kEX
real(sp), dimension (ncoso)        :: cosoEX
real(sp), dimension (ncosi)        :: cosiEX
real(sp), dimension (ncoso,ncosi)  :: intEX
common /EXP/ cosoEX, cosiEX, intEX

funcs = 0.0_sp
do ii = 1, ncosi
   cosi = cosiEX (ii) 
   do io = 1, ncoso
      coso = cosoEX (io) 
      call fsc1TH ( kEX, coso, xv, fsc, dfscCJ )
      funcs = funcs + ( ( Abs (fsc) )**2  - intEX(io,ii)  )**2
   enddo
enddo

do ix = 1, NP

   if ( xv(ix) .gt.  PI ) then
      funcs  = funcs * ( 1. + favg * ( xv(ix) - PI )**2 )
   endif
   if ( xv(ix) .lt.  0. ) then
      funcs  = funcs * ( 1. + favg * ( xv(ix) + 0. )**2 )
   endif

enddo


END FUNCTION funcs
