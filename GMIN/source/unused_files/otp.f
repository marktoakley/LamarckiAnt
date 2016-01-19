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
C  Energy and gradient for OTP.
C
      SUBROUTINE OTP(X,V,EOTP,GTEST,SECT)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2
      INTEGER IA,IM
      DOUBLE PRECISION X(3*MXATMS), V(3*MXATMS), X1, X2, Y1, Y2, Z1, Z2, FOTP, DFOTP, PHI,
     1                 EOTP, DUMMY, RALPHA12, RALPHA22, RDIST, SITE(6,3),
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, RR2, CX, CY, CZ, SUM,
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, C2A2, C2A1,
     5                 D1, D2, D3, D4, D5, D6, D10, D11, D12, DIST, DUMMY2,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, XDUMM,
     7                 DCM1, DCM2, DCM3, DCM4, DCM5, DCM6
      DOUBLE PRECISION CMDIST(MXATMS*3)
      DOUBLE PRECISION RMIN,INFTY,RC,RC2,RR
      DOUBLE PRECISION INT1,INT2,INT3,INT4,INT5
      DOUBLE PRECISION RMAX,EMF,VMF,D_R_VMF,D_RC_VMF
      DOUBLE PRECISION x3,y3,z3
      PARAMETER(INFTY=1.0D30)
      LOGICAL EVAP, evapreject
      EXTERNAL FUNC1,FUNC2,FUNC3,FUNC4
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      COMMON /EV/ EVAP, evapreject

      COMMON /RRR/RR,RC2
C
C  Attractive and repulsive function definitions and their derivatives.
C
       FOTP(XDUMM)=4.0D0*(1.0D0/XDUMM**12-1.0D0/XDUMM**6)
       DFOTP(XDUMM)=24.0D0/XDUMM**7*(-2.0D0/XDUMM**6+1.0D0)
C
C Derivatives of R(i,j) with respect to rigid body coordinates.
C P(i) = (p1x,p1y,p1z), P(j)=(p2x,p2y,p2z) are the site cordinates in
C the reference geometry for the molecule at the origin.
C
      D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @    rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + 
     @    m2*p2z*s2 + x1 - x2)

      D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @    rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + n2*p2x*s2 - 
     @    l2*p2z*s2 + y1 - y2)

      D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @    rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + 
     @    l2*p2y*s2 + z1 - z2)

      D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @      rdist*((c3a1*(2*l1*p1x*(-1 + l12*ralpha12) + 
     @            (m1*p1y + n1*p1z)*(-1 + 2*l12*ralpha12)) + 
     @         (l1*l12*p1x + l12*m1*p1y - l1*n1*p1y + l1*m1*p1z + l12*n1*p1z)*
     @          ralpha12*s1 + l1*(c2a1*(n1*p1y - m1*p1z)*ralpha12 - p1x*s1))*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     @         c3a1*m1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-(l1*p1y) + p1z + l1* (l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + m1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (-(c2a1*l1*(-(m1*p1x) + l1*p1y)*ralpha12) + 
     @         c3a1*n1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-p1y - l1*p1z + l1*(-(m1*p1x) + l1*n1*p1x + l1*p1y + m1*n1*p1y + n1**2*p1z)*ralpha12)*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @ rdist*((-(c2a1*m1*(-(n1*p1y) + m1*p1z)*ralpha12) + 
     @         c3a1*l1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-(m1*p1x) - p1z + m1* (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (-(c3a1*(l1*p1x + 2*m1*p1y + n1*p1z)) + 
     @         c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @         (-(m1*p1y) + (l1*m12*p1x + m1*n1*p1x + m1*m12*p1y - l1*m1*p1z + 
     @               m12*n1*p1z)*ralpha12)*s1)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     @         c3a1*n1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (p1x - m1*(m1 - l1*n1)*p1x*ralpha12 + 
     @            m1*((l1 + m1*n1)*p1y*ralpha12 + p1z*(-1 + n1**2*ralpha12)))*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @ rdist*((c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     @         c3a1*l1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + (-(n1*p1x) + p1y + n1*
     @             (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (-(c2a1*n1*(n1*p1x - l1*p1z)*ralpha12) + 
     @         c3a1*m1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-p1x - n1*p1y + n1*(l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + 
     @               m1*n1*p1z)*ralpha12)*s1)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (-(c3a1*(l1*p1x + m1*p1y + 2*n1*p1z)) + 
     @         c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     @         2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @         (-(n1*p1z) + (-(m1*n1*p1x) + l1*n12*p1x + l1*n1*p1y + m1*n12*p1y + 
     @               n1*n12*p1z)*ralpha12)*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      D10(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @         rdist*((c3a2*(2*l2*p2x + m2*p2y + n2*p2z - 
     @            2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
     @         (l2*l22*p2x + l22*m2*p2y - l2*n2*p2y + l2*m2*p2z + l22*n2*p2z)*
     @          ralpha22*s2 + l2*(-(c2a2*n2*p2y*ralpha22) + c2a2*m2*p2z*ralpha22 + 
     @            p2x*s2))*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - 
     @         c2a2*p2x + c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + 
     @         (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (-(c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22) + 
     @         c3a2*m2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - (-(l2*p2y) + p2z + l2*
     @             (l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*
     @          s2)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (c2a2*l2*(-(m2*p2x) + l2*p2y)*ralpha22 + 
     @         c3a2*n2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @    (p2y + l2*p2z - l2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + n2**2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      D11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @ rdist*((c2a2*m2*(-(n2*p2y) + m2*p2z)*ralpha22 + 
     @         c3a2*l2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (m2*p2x + p2z - m2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c3a2*(l2*p2x + 2*m2*p2y + n2*p2z - 
     @            2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
     @         (l2*m22*p2x + m2*n2*p2x + m2*m22*p2y - l2*m2*p2z + m22*n2*p2z)*
     @          ralpha22*s2 + m2*(c2a2*(n2*p2x - l2*p2z)*ralpha22 + p2y*s2))*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (-(c2a2*m2*(m2*p2x - l2*p2y)*ralpha22) + 
     @         c3a2*n2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (-p2x + m2*p2z - m2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + 
     @               n2**2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      D12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
     @ rdist*((-(c2a2*n2*(n2*p2y - m2*p2z)*ralpha22) + 
     @         c3a2*l2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (n2*p2x - p2y - n2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c2a2*n2*(n2*p2x - l2*p2z)*ralpha22 + 
     @         c3a2*m2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (p2x + n2*p2y - n2*(l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c3a2*(l2*p2x + m2*p2y + 2*n2*p2z - 
     @            2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (n2*(m2*p2x - l2*p2y) - n22*(l2*p2x + m2*p2y + n2*p2z))*ralpha22*
     @          s2 + n2*(-(c2a2*m2*p2x*ralpha22) + c2a2*l2*p2y*ralpha22 + p2z*s2))*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))
c***************************
c CALCULATION MF INTERACTION
****************************
c********************************************************************** 
c
      DCM1(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) = 
     @ rdist*(-CX + c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) + 
     @ n1*p1y*s1 - m1*p1z*s1 + x1)
      DCM2(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) = 
     @ rdist*(-CY + c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - 
     @ n1*p1x*s1 + l1*p1z*s1 + y1)
      DCM3(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) = 
     @ rdist*(-CZ + c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) + 
     @ m1*p1x*s1 - l1*p1y*s1 + z1)
      DCM4(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) = 
     @ rdist*((c3a1*(2*l1*p1x*(-1 + l12*ralpha12) + 
     @      (m1*p1y + n1*p1z)*(-1 + 2*l12*ralpha12)) + 
     @   (l1*l12*p1x + l12*m1*p1y - l1*n1*p1y + l1*m1*p1z + l12*n1*p1z)*
     @    ralpha12*s1 + l1*(c2a1*(n1*p1y - m1*p1z)*ralpha12 - p1x*s1))*
     @ (-CX + c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) + n1*p1y*s1 - 
     @   m1*p1z*s1 + x1) + (c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     @   c3a1*m1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (-(l1*p1y) + p1z + l1*
     @       (l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + m1*n1*p1z)*ralpha12)*
     @    s1)*(-CY + c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - 
     @   n1*p1x*s1 + l1*p1z*s1 + y1) + 
     @ (-(c2a1*l1*(-(m1*p1x) + l1*p1y)*ralpha12) + 
     @   c3a1*n1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (-p1y - l1*p1z + l1*(-(m1*p1x) + l1*n1*p1x + l1*p1y + m1*n1*p1y + 
     @         n1**2*p1z)*ralpha12)*s1)*
     @ (-CZ + c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) + m1*p1x*s1 - 
     @   l1*p1y*s1 + z1))
        DCM5(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) = 
     @  rdist*((-(c2a1*m1*(-(n1*p1y) + m1*p1z)*ralpha12) + 
     @   c3a1*l1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (-(m1*p1x) - p1z + m1*
     @       (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @    s1)*(-CX + c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) + 
     @   n1*p1y*s1 - m1*p1z*s1 + x1) + 
     @ (-(c3a1*(l1*p1x + 2*m1*p1y + n1*p1z)) + 
     @   c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     @   2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @   (-(m1*p1y) + (l1*m12*p1x + m1*n1*p1x + m1*m12*p1y - l1*m1*p1z + 
     @         m12*n1*p1z)*ralpha12)*s1)*
     @ (-CY + c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - n1*p1x*s1 + 
     @   l1*p1z*s1 + y1) + (c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     @   c3a1*n1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (p1x - m1*(m1 - l1*n1)*p1x*ralpha12 + 
     @      m1*((l1 + m1*n1)*p1y*ralpha12 + p1z*(-1 + n1**2*ralpha12)))*s1)*
     @ (-CZ + c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) + m1*p1x*s1 - 
     @   l1*p1y*s1 + z1))
        DCM6(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist) =
     @   rdist*((c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     @   c3a1*l1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (-(n1*p1x) + p1y + n1*
     @       (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @    s1)*(-CX + c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) + 
     @   n1*p1y*s1 - m1*p1z*s1 + x1) + 
     @ (-(c2a1*n1*(n1*p1x - l1*p1z)*ralpha12) + 
     @   c3a1*m1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @   (-p1x - n1*p1y + n1*(l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + 
     @         m1*n1*p1z)*ralpha12)*s1)*
     @ (-CY + c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - n1*p1x*s1 + 
     @   l1*p1z*s1 + y1) + (-(c3a1*(l1*p1x + m1*p1y + 2*n1*p1z)) + 
     @   c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     @   2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @   (-(n1*p1z) + (-(m1*n1*p1x) + l1*n12*p1x + l1*n1*p1y + m1*n12*p1y + 
     @         n1*n12*p1z)*ralpha12)*s1)*
     @ (-CZ + c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) + m1*p1x*s1 - 
     @   l1*p1y*s1 + z1))
C************************************************************
C
C  The six reference site positions per capped pentagon. These need to
C  be multiplied by RAD.
C
c******************** 
c      OTP Molecule Geometry
c
       Phi=37.5D0
       Phi=Phi/180.0D0*3.1415927D0
 
       x1=0.0D0
       y1=0.0D0
       z1=0.0D0
 
       x2=dsin(Phi)
       y2=dcos(Phi)
       z2=0.0D0
 
       x3=2.0D0*dsin(Phi)
       y3=0.0D0
       z3=0.0D0
 
       cx=(x1+x2+x3)/3.0D0
       cy=(y1+y2+y3)/3.0D0
       cz=(z1+z2+z3)/3.0D0
 
       x1=x1-cx
       y1=y1-cy
       z1=z1-cz
       x2=x2-cx
       y2=y2-cy
       z2=z2-cz
       x3=x3-cx
       y3=y3-cy
       z3=z3-cz
  
       SITE(1,1)=X1
       SITE(1,2)=Y1
       SITE(1,3)=Z1
       SITE(2,1)=X2
       SITE(2,2)=Y2
       SITE(2,3)=Z2
       SITE(3,1)=X3
       SITE(3,2)=Y3
       SITE(3,3)=Z3

       DUMMY=0.0D0
       GX1=0.0D0
       GY1=0.0D0
       GZ1=0.0D0
       GL1=0.0D0
       GM1=0.0D0
       GN1=0.0D0
       GL2=0.0D0
       GM2=0.0D0
       GN2=0.0D0

c****************

      IF (SECT) THEN
         PRINT*,'ERROR - no OTP second derivatives!'
         STOP
      ENDIF
      EVAP=.FALSE.
      EOTP=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO
      
      CX=0.0D0
      CY=0.0D0
      CZ=0.0D0

      DO J1=1,NATOMS/2
         CX=CX+X(6*(J1-1)+1)
         CY=CY+X(6*(J1-1)+2)
         CZ=CZ+X(6*(J1-1)+3)
      ENDDO

      CX=CX*2/NATOMS
      CY=CY*2/NATOMS
      CZ=CZ*2/NATOMS
      CX=0.0D0
      CY=0.0D0
      CZ=0.0D0
      PRINT*,'CX,CY,CZ=',CX,CY,CZ
C
C  Potential energy first.
C
      DO J1=1,NATOMS/2

         X1=X(6*(J1-1)+1)
         Y1=X(6*(J1-1)+2)
         Z1=X(6*(J1-1)+3)
         L1=X(6*(J1-1)+4)
         M1=X(6*(J1-1)+5)
         N1=X(6*(J1-1)+6)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-10)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0001D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF
C        WRITE(*,'(A,6F15.5)') 'ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1=',ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1

         DO J2=J1+1,NATOMS/2
            X2=X(6*(J2-1)+1)
            Y2=X(6*(J2-1)+2)
            Z2=X(6*(J2-1)+3)
            L2=X(6*(J2-1)+4)
            M2=X(6*(J2-1)+5)
            N2=X(6*(J2-1)+6)
            L22=L2**2
            M22=M2**2
            N22=N2**2
            ALPHA2=SQRT(L22+M22+N22)
            RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-10)
            CA2=COS(ALPHA2)
            C2A2=CA2
            IF (ALPHA2.LT.0.0001D0) THEN
C              C3A2=-ALPHA2/2+ALPHA2**3/24
               C3A2=-0.5D0+ALPHA2**2/24.0D0
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  Sum over LJ sites.
C
            DO K1=1,3

               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)

               CMDIST(3*(J1-1)+K1)=
     @ Sqrt((-CX + c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) + 
     @ (n1*p1y - m1*p1z)*s1 + x1)**2 + 
     @ (-CY + c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) + (-(n1*p1x) 
     @ + l1*p1z)*s1 + y1)**2 + 
     @ (-CZ + c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) + (m1*p1x - 
     @ l1*p1y)*s1 + z1)**2)
C              PRINT*,'J1,K1,3*(J1-1)+K1,CMDIST=',J1,K1,3*(J1-1)+K1,CMDIST(3*(J1-1)+K1)

           DO K2=1,3

                  P2X=SITE(K2,1)
                  P2Y=SITE(K2,2)
                  P2Z=SITE(K2,3)

               CMDIST(3*(J2-1)+K2)=
     @ Sqrt((-CX + c2a2*p2x - c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + 
     @ (n2*p2y - m2*p2z)*s2 + x2)**2 + 
     @ (-CY + c2a2*p2y - c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n2*p2x) 
     @ + l2*p2z)*s2 + y2)**2 + 
     @ (-CZ + c2a2*p2z - c3a2*n1*(l2*p2x + m2*p2y + n2*p2z) + (m2*p2x - 
     @ l2*p2y)*s2 + z2)**2)
C              PRINT*,'J2,K2,3*(J2-1)+K2,CMDIST=',J2,K2,3*(J2-1)+K2,CMDIST(3*(J2-1)+K2)

                  DIST=
     @  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     @      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     @   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     @      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     @   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     @      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)

                  DUMMY=DUMMY+FOTP(DIST)
C                 WRITE(*,'(A,2I5,3G20.10)') 'K1,K2,DIST,FOTP,DUMMY=',K1,K2,DIST,FOTP(DIST),DUMMY

                  IF (GTEST) THEN
                     DUMMY2=DFOTP(DIST)
                     RDIST =1.0D0/DIST

                     GX1=GX1+DUMMY2*D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GY1=GY1+DUMMY2*D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GZ1=GZ1+DUMMY2*D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GL1=GL1+DUMMY2*D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GM1=GM1+DUMMY2*D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GN1=GN1+DUMMY2*D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GL2=GL2+DUMMY2*D10(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GM2=GM2+DUMMY2*D11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                     GN2=GN2+DUMMY2*D12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)
                  ENDIF
               ENDDO
            ENDDO
   
            EOTP=EOTP+DUMMY
            VT(J1)=VT(J1)+DUMMY
            VT(J2)=VT(J2)+DUMMY

            V(6*(J1-1)+1)=V(6*(J1-1)+1)+GX1
            V(6*(J1-1)+2)=V(6*(J1-1)+2)+GY1
            V(6*(J1-1)+3)=V(6*(J1-1)+3)+GZ1
            V(6*(J1-1)+4)=V(6*(J1-1)+4)+GL1
            V(6*(J1-1)+5)=V(6*(J1-1)+5)+GM1
            V(6*(J1-1)+6)=V(6*(J1-1)+6)+GN1
            V(6*(J2-1)+1)=V(6*(J2-1)+1)-GX1
            V(6*(J2-1)+2)=V(6*(J2-1)+2)-GY1
            V(6*(J2-1)+3)=V(6*(J2-1)+3)-GZ1
            V(6*(J2-1)+4)=V(6*(J2-1)+4)+GL2
            V(6*(J2-1)+5)=V(6*(J2-1)+5)+GM2
            V(6*(J2-1)+6)=V(6*(J2-1)+6)+GN2
         ENDDO
      ENDDO

C     WRITE(*,'(A,G20.10)') 'energy=',EOTP
C     PRINT*,'coords:'
C     WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'gradient:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)
C     WRITE(*,'(A,2G20.10)') 'ALPHA1,ALPHA2=',ALPHA1,ALPHA2

        RMAX=0.0D0
        EMF =0.0D0

         DO J1=1,NATOMS/2

            DO K1=1,3

           J2=3*(J1-1)+K1
           DIST =CMDIST(J2)
           
           IF(DIST.GT.RMAX) THEN
             RMAX=DIST
             IA =K1
             IM =J1 
           ENDIF

           ENDDO

         END DO
C        PRINT*,'RMAX=',RMAX

          RC    = DSQRT(RMAX)+LAMBDA
          RC2   = RC*RC
          RADIUS= RC2
          SUM   = 0.0D0
C         PRINT*,'RADIUS=',RADIUS
           
C***********************************

         DO J1=1,NATOMS/2

         X1=X(6*(J1-1)+1)
         Y1=X(6*(J1-1)+2)
         Z1=X(6*(J1-1)+3)
         L1=X(6*(J1-1)+4)
         M1=X(6*(J1-1)+5)
         N1=X(6*(J1-1)+6)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-10)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0001D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF
c***************************

           GX1=0.0D0
           GY1=0.0D0
           GZ1=0.0D0
           GL1=0.0D0
           GM1=0.0D0
           GN1=0.0D0

           DO K1=1,3

           J2=3*(J1-1)+K1
           
           RR=CMDIST(J2)
C          PRINT*,'J2,CMDIST=',J2,CMDIST(J2)
           RR2 =RR*RR

           RMIN=RC-RR
           RMAX=RC+RR
C          PRINT*,'RMIN,RMAX=',RMIN,RMAX

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

           VT(J1)= VT(J1)+VMF

               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)

               RDIST=1.0D0/RR

           GX1=GX1+D_R_VMF*DCM1(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GY1=GY1+D_R_VMF*DCM2(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GZ1=GZ1+D_R_VMF*DCM3(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GL1=GL1+D_R_VMF*DCM4(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GM1=GM1+D_R_VMF*DCM5(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GN1=GN1+D_R_VMF*DCM6(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)

           SUM = SUM + D_RC_VMF
         
          ENDDO ! MOLECULE SITES

C           V(6*(J1-1)+1)=V(6*(J1-1)+1)+GX1
C           V(6*(J1-1)+2)=V(6*(J1-1)+2)+GY1
C           V(6*(J1-1)+3)=V(6*(J1-1)+3)+GZ1
C           V(6*(J1-1)+4)=V(6*(J1-1)+4)+GL1
C           V(6*(J1-1)+5)=V(6*(J1-1)+5)+GM1
C           V(6*(J1-1)+6)=V(6*(J1-1)+6)+GN1
            V(6*(J1-1)+1)=GX1
            V(6*(J1-1)+2)=GY1
            V(6*(J1-1)+3)=GZ1
            V(6*(J1-1)+4)=GL1
            V(6*(J1-1)+5)=GM1
            V(6*(J1-1)+6)=GN1

          ENDDO ! MOLECULES


C*************
C Add energies
C*************

C         EOTP=EOTP+EMF
          EOTP=EMF

C************************
C PARTICULAR CASE RR=RMAX
C************************

         J1=IM
         
         X1=X(6*(J1-1)+1)
         Y1=X(6*(J1-1)+2)
         Z1=X(6*(J1-1)+3)
         L1=X(6*(J1-1)+4)
         M1=X(6*(J1-1)+5)
         N1=X(6*(J1-1)+6)

         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-10)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0001D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF

               P1X=SITE(IA,1)
               P1Y=SITE(IA,2)
               P1Z=SITE(IA,3)

           J2   = 3*(IM-1)+IA
           DIST = CMDIST(J2)
           RDIST= 1.0D0/DIST

           GX1=SUM*DCM1(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GY1=SUM*DCM2(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GZ1=SUM*DCM3(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GL1=SUM*DCM4(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GM1=SUM*DCM5(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)
           GN1=SUM*DCM6(p1x,p1y,p1z,x1,y1,z1,n1,l1,m1,CX,CY,CZ,rdist)

            V(6*(IM-1)+1)=V(6*(IM-1)+1)+GX1
            V(6*(IM-1)+2)=V(6*(IM-1)+2)+GY1
            V(6*(IM-1)+3)=V(6*(IM-1)+3)+GZ1
            V(6*(IM-1)+4)=V(6*(IM-1)+4)+GL1
            V(6*(IM-1)+5)=V(6*(IM-1)+5)+GM1
            V(6*(IM-1)+6)=V(6*(IM-1)+6)+GN1

            
      RETURN
      END


c***************************************
      DOUBLE PRECISION FUNCTION func1(y)
c***************************************
      USE commons
      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RR,RC2
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU
      
        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     #  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        func1   =Y2*GDR*VLJ
                     
      END

c***************************************
      DOUBLE PRECISION FUNCTION func2(y)
c***************************************

      USE commons
      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC2,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      COMMON /RRR/ RR,RC2
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     #  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ
        func2 = Y2*FF*(RC2-RR*RR-Y2)/(2.0D0*RR*Y)

      END

c***************************************
      DOUBLE PRECISION FUNCTION func3(y)
c***************************************

      USE commons
      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC2,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      COMMON /RRR/ RR,RC2
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     #  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ
        func3 = Y2*FF*(RC2+RR*RR-Y2)/(2.0D0*RR*RR*Y)

      END

c***************************************
      DOUBLE PRECISION FUNCTION func4(y)
c***************************************

      USE commons
      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC2,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      COMMON /RRR/ RR,RC2
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     #  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ
        func4 = Y*FF/RR

      END

      SUBROUTINE OTPPARAMMF
      IMPLICIT NONE
      USE commons
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592654)
      DOUBLE PRECISION H,M,GD,L,A,T,DY,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      open(unit=15,file='parameters.dat',status='old')

      READ(15,*)H
      READ(15,*)M
      READ(15,*)GD
      READ(15,*)L
      READ(15,*)A
      READ(15,*)B
      READ(15,*)T
      READ(15,*)EPSILON
      READ(15,*)RHOHAT
      READ(15,*)LAMBDA
      READ(15,*)MU

      FACT=MU*PI*RHOHAT

      CLOSE(15)

      RETURN
      END
