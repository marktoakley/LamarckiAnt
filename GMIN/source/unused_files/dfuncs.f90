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
FUNCTION dfuncs ( xv )
USE nrtype
IMPLICIT NONE
include 'NP.H'

INTERFACE
!  s = \sum_{i} f_{i} (\vec)
   FUNCTION funcs(xv)
      USE nrtype
      REAL(SP), DIMENSION(:), INTENT(IN)       :: xv
      REAL(SP)                                 :: funcs
   END FUNCTION funcs
END INTERFACE

REAL(SP), DIMENSION(:), INTENT(IN)    :: xv
REAL(SP), DIMENSION(size(xv))         :: dfuncs
complex(dpc)                          :: fsc
complex(dpc),dimension(NP)            :: dfscCJ
REAL(SP)                              :: coso, cosi, dIdx
real(sp), parameter                   :: favg = 1.e6
integer(i4b)                          :: ix, io, ii

real(sp)                           :: kEX
common /ENE/ kEX
real(sp), dimension (ncoso)        :: cosoEX
real(sp), dimension (ncosi)        :: cosiEX
real(sp), dimension (ncoso,ncosi)  :: intEX
common /EXP/ cosoEX, cosiEX, intEX

do ix = 1, size(xv)

   dfuncs (ix) = 0.0
   do ii = 1, ncosi
      cosi = cosiEX (ii) 
      do io = 1, ncoso
         coso = cosoEX (io) 
         call fsc1TH ( kEX, coso, xv, fsc, dfscCJ )

         dIdx = 2. * Real ( fsc * dfscCJ(ix) )
         dfuncs (ix) =  &
                        2. * ( ( Abs (fsc) )**2  - intEX(io,ii)  ) * dIdx

      enddo   ! io
   enddo   ! ii

   if ( xv(ix) .gt.  PI ) then
      dfuncs (ix) = dfuncs (ix) * ( 1. + favg * ( xv(ix) - PI )**2 ) + &
                    2.*favg* ( xv(ix) - PI ) * funcs (xv)
   endif
   if ( xv(ix) .lt.  0. ) then
      dfuncs (ix) = dfuncs (ix) * ( 1. + favg * ( xv(ix) + 0. )**2 ) + &
                    2.*favg* ( xv(ix) + 0. ) * funcs (xv)
   endif


enddo  ! ix
END FUNCTION dfuncs
