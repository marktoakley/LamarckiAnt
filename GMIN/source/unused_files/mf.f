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
C
C  Energy and Gradient for the Mean Field Potential.
C
      SUBROUTINE MF(X,V,EMF,GTEST)
      USE commons


      LOGICAL GTEST
      INTEGER J1, J2, J3, J4, IJ, I, I2
      DOUBLE PRECISION X(3*MXATMS), DIST, V(3*MXATMS), G(MXATMS,MXATMS) 
      DOUBLE PRECISION R6,EMF,DUMMYX,DUMMYY,DUMMYZ,XMUL2,DUMMY,RMIN,RMAX
      DOUBLE PRECISION R(MXATMS),RR,RR2,RC,RC2,SUM,DD,FF
      DOUBLE PRECISION VMF,D_R_VMF,D_RC_VMF,Y,Y2,DDY,GDR,VLJ
      DOUBLE PRECISION INT1,INT2,INT3,INT4,INT5,FINT
      DOUBLE PRECISION INFTY
      PARAMETER (INFTY=1.E30)
      EXTERNAL FUNC1,FUNC2,FUNC3,FUNC4
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      COMMON /RRR/ RR,RC2
 
        RMAX=0.0D0
        EMF =0.0D0

         DO J1=1,NATOMS

           J3=3*J1
           DIST =X(J3-2)**2+X(J3-1)**2+X(J3)**2
           R(J1)=DSQRT(DIST)

           IF(DIST.GT.RMAX) THEN
             RMAX=DIST
             IJ  =J1
           ENDIF

         ENDDO

          RC    = DSQRT(RMAX)+LAMBDA
          RC2   = RC*RC
          RADIUS= RC2
          SUM   = 0.0D0
         
         DO J1=1,NATOMS
           
           RR  =R(J1)
           RR2 =RR*RR

           RMIN=RC-RR
           RMAX=RC+RR

C CALCULATION OF POTENTIAL AND DERIVATIVES

       call qromb(func1,RMIN,RMAX,INT1)
       call qromb(func2,RMIN,RMAX,INT2)
       call qromb(func3,RMIN,RMAX,INT4)
       call qromb(func4,RMIN,RMAX,INT5)

C CALCULATION OF THE INTEGRAL I3 (TO INFINITY)

      call qromo(func1,RMAX,INFTY,INT3)

        VMF      = FACT * (INT1-INT2+2.0D0*INT3)
        D_R_VMF  = FACT * INT4
        D_RC_VMF = -FACT * RC * INT5

           EMF   = EMF+VMF

           VT(J1)= VT(J1)+VMF/4.0D0

             DD  = D_R_VMF/RR

           J3     =3*J1 
           V(J3)  =DD*X(J3)  
           V(J3-1)=DD*X(J3-1)
           V(J3-2)=DD*X(J3-2)

           SUM = SUM + D_RC_VMF
         
          ENDDO

C************************
C PARTICULAR CASE RR=RMAX
C************************

          J3=3*IJ
          RR=R(IJ)
          DD=SUM/RR

           V(J3)  =V(J3)  +DD*X(J3)  
           V(J3-1)=V(J3-1)+DD*X(J3-1)
           V(J3-2)=V(J3-2)+DD*X(J3-2)

      RETURN
      END

