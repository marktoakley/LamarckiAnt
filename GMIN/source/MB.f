C   OPTIM: A program for optimizing geometries and calculating reaction pathways
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
      SUBROUTINE MB(Z,VNEW,POTEL,GTEST,STEST)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      DOUBLE PRECISION Z(3*NATOMS), POTEL
      DOUBLE PRECISION VNEW(3*NATOMS), X, Y

      X=Z(1)
      Y=Z(2)

C     POTEL=-170.0D0*DEXP(-49.0D0-46.0D0*X-13.0D0*X**2+50.0D0*Y+22.0D0*X*Y-13.0D0*Y**2)/2.0D0)
C    1      -200.0D0*DEXP(-1.0D0+2.0D0*X-X**2-10.0D0*Y**2)
C    2      -100.0D0*DEXP(-5.0D0/2.0D0-X+10.0D0*Y-10.0D0*Y**2)+
C    3       +15.0D0*DEXP((8.0D0+8.0D0*X+7.0D0*X**2-8.0D0*Y+6.0D0*X*Y+7.0D0*Y**2)/10.0D0)
      POTEL=
     1      -200.0D0*DEXP(-1.0D0*(X-1.0D0)**2+ 0.0D0*(X-1.0D0)*(Y-0.0D0)-10.0D0*(Y-0.0D0)**2)
     2      -100.0D0*DEXP(-1.0D0*(X+0.0D0)**2+ 0.0D0*(X+0.0D0)*(Y-0.5D0)-10.0D0*(Y-0.5D0)**2)
     3      -170.0D0*DEXP(-6.5D0*(X+0.5D0)**2+11.0D0*(X+0.5D0)*(Y-1.5D0) -6.5D0*(Y-1.5D0)**2)
     4       +15.0D0*DEXP( 0.7D0*(X+1.0D0)**2+ 0.6D0*(X+1.0D0)*(Y-1.0D0) +0.7D0*(Y-1.0D0)**2)

      IF (GTEST) THEN
      
      VNEW(1)=400*DEXP(-1.0D0+2.0D0*X-X**2-10.0D0*Y**2)*(-1.0D0+X)+   
     1  200.0D0*DEXP(-5.0D0/2.0D0-X**2+10.0D0*Y-10.0D0*Y**2)*X+
     2  3.0D0*DEXP((8.0D0+8.0D0*X+7.0D0*X**2-8.0D0*Y+6.0D0*X*Y
     3  +7.0D0*Y**2)/10.0D0)*(4.0D0+7.0D0*X+3.0D0*Y)- 
     4  170.0D0*DEXP((-49.0D0-46.0D0*X-13.0D0*X**2+50.0D0*Y+22.0D0*X*Y
     5  -13.0D0*Y**2)/2.0D0)*(-23.0D0-13.0D0*X+11.0D0*Y)
      VNEW(2)=-170.0D0*DEXP((-49.0D0-46.0D0*X-13.0D0*X**2+50.0D0*Y
     1  +22.0D0*X*Y-13.0D0*Y**2)/2.0D0)*
     2   (25.0D0+11.0D0*X-13.0D0*Y)+2000.0D0*DEXP(-5.0D0/2.0D0-X**2
     3  +10.0D0*Y-10.0D0*Y**2)*
     4   (-1.0D0/2.0D0+Y)+4000.0D0*DEXP(-1.0D0+2.0D0*X-X**2-10.0D0*Y**2)
     5  *Y+3.0D0*DEXP((8.0D0+8.0D0*X+7.0D0*X**2-8.0D0*Y+6.0D0*X*Y
     6  +7.0D0*Y**2)/10.0D0)*(-4.0D0+3.0D0*X+7.0D0*Y)
      VNEW(3)=0.0D0
      ENDIF

      IF (STEST) THEN
         STOP
!      HESS(1,1)=2210*DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0)+
!     1  400*DEXP(-1+2*X-X**2-10*Y**2)+
!     2  200*DEXP(-5.0D0/2.0D0-X**2+10*Y-10*Y**2)+
!     3  21*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)-
!     4  800*DEXP(-1+2*X-X**2-10*Y**2)*(-1+X)**2-
!     5  400*DEXP(-5.0D0/2.0D0-X**2+10*Y-10*Y**2)*X**2+
!     6  3*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)*
!     7    (4+7*X+3*Y)**2/5-
!     8  170*DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0)*
!     9   (-23-13*X+11*Y)**2
!      HESS(1,2)=-1870*DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0) 
!     1  +9*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)+
!     2  2000*DEXP(-5.0D0/2.0D0-X**2+10*Y-10*Y**2)*X*(1-2*Y)+
!     3  8000*DEXP(-1+2*X-X**2-10*Y**2)*(1-X)*Y+
!     4  3*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)*(4+7*X+3*Y)*
!     5    (-4+3*X+7*Y)/5.0D0-170*
!     6   DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0)*
!     7   (25+11*X-13*Y)*(-23-13*X+11*Y)
!      HESS(2,1)=HESS(1,2)
!      HESS(2,2)=2210*DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0)+ 
!     1  4000*DEXP(-1+2*X-X**2-10*Y**2)+
!     2  2000*DEXP(-5.0D0/2.0D0-X**2+10*Y-10*Y**2)+
!     3  21*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)-
!     4  170*DEXP((-49-46*X-13*X**2+50*Y+22*X*Y-13*Y**2)/2.0D0)*
!     5   (25+11*X-13*Y)**2-80000*DEXP(-1+2*X-X**2-10*Y**2)*Y**2-
!     6  10000*DEXP(-5.0D0/2.0D0-X**2+10*Y-10*Y**2)*(-1+2*Y)**2+
!     7  3*DEXP((8+8*X+7*X**2-8*Y+6*X*Y+7*Y**2)/10.0D0)*
!     8    (-4+3*X+7*Y)**2/5
      ENDIF

!     PRINT '(A,4F20.10)','X,Y,GX,GY=',X,Y,VNEW(1),VNEW(2)

      RETURN
      END
