!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!  Energy and gradient for 1-D step landscape with three minima.
!
SUBROUTINE MODEL1(X,V,ENERGY,QE,QX)
USE commons
IMPLICIT NONE
DOUBLE PRECISION ENERGY, X(3*NATOMS), DIST, V(3*NATOMS), QE, QX
DOUBLE PRECISION, PARAMETER :: EDGE=1.0D3
  
VT(1)=0.0D0
IF (X(1).LE.MSTART) THEN
   ENERGY=EDGE
   QE=ME1
   QX=(MBSTART1+MSTART)/2.0D0
ELSEIF (X(1).LE.MBSTART1) THEN
   ENERGY=ME1
   QE=ME1
   QX=(MBSTART1+MSTART)/2.0D0
ELSEIF (X(1).LE.MBFINISH1) THEN
   ENERGY=ME2+MBHEIGHT1
   IF (X(1).LE.(MBFINISH1+MBSTART1)/2.0D0) THEN
      QE=ME1
      QX=(MBSTART1+MSTART)/2.0D0
   ELSE
      QE=ME2
      QX=(MBSTART2+MBFINISH1)/2.0D0
   ENDIF
ELSEIF (X(1).LE.MBSTART2) THEN
   ENERGY=ME2
   QE=ME2
   QX=(MBSTART2+MBFINISH1)/2.0D0
ELSEIF (X(1).LE.MBFINISH2) THEN
   ENERGY=ME3+MBHEIGHT2
   IF (X(1).LE.(MBFINISH1+MBSTART1)/2.0D0) THEN
      QE=ME2
      QX=(MBSTART2+MBFINISH1)/2.0D0
   ELSE
      QE=ME3
      QX=(MFINISH+MBFINISH2)/2.0D0
   ENDIF
ELSEIF (X(1).LE.MFINISH) THEN
   ENERGY=ME3
   QE=ME3
   QX=(MFINISH+MBFINISH2)/2.0D0
ELSE
   ENERGY=EDGE
   QE=ME3
   QX=(MFINISH+MBFINISH2)/2.0D0
ENDIF

! PRINT '(A,4F10.1)','X,E,QE,QX=',X(1),ENERGY,QE,QX

RETURN
END
