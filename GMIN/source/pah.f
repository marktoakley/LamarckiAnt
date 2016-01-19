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
C  Energy and gradient for rigid body PAH molecule using new rigid body
C  derivative functions etc.
C
      SUBROUTINE PAH(X,V,EPAH,GTEST,SECT)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2, NAT2
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, FC, DFC,
     1                 EPAH, DUMMY, RALPHA12, RALPHA22, RDIST, RM6, RM1,
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, C6cc, C12cc, C6hh, C12hh,
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, C2A2, C2A1,
     5                 D1, D2, D3, D4, D5, D6, DA, DB, DC, DIST, DUMMY2,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, Q1, Q2, RDIST6,
     7           D1S, D2S, D3S, D4S, D5S, D6S, DAS, DBS, DCS, CHARGE(36), HTOKJ
      PARAMETER (HTOKJ=1389.354848D0) ! conversion factor for coulomb energy and gradients
      DOUBLE PRECISION FLJcc, FLJhh, FLJch, DFLJcc, DFLJhh, DFLJch, C6ch, C12ch
C
C  Statement functions.
C  Site-site energy terms - Coulombic and LJ. 
C
      FLJcc(RM6)=(C12cc*RM6-C6cc)*RM6
      DFLJcc(RM1,RM6)=-6.0D0*(-C6cc+2.0D0*C12cc*RM6)*RM6*RM1
      FLJhh(RM6)=(C12hh*RM6-C6hh)*RM6
      DFLJhh(RM1,RM6)=-6.0D0*(-C6hh+2.0D0*C12hh*RM6)*RM6*RM1
      FLJch(RM6)=(C12ch*RM6-C6ch)*RM6
      DFLJch(RM1,RM6)=-6.0D0*(-C6ch+2.0D0*C12ch*RM6)*RM1*RM6

      FC(Q1,Q2,RM1)=Q1*Q2*RM1
      DFC(Q1,Q2,RM1)=-Q1*Q2*RM1**2
C
      DO J1=1,NRBSITES
         SITE(J1,3)=0.0D0
      ENDDO

      SITE( 1,1)=1.42d0
      SITE( 1,2)=0.0d0
      CHARGE( 1)=-0.0066d0
      SITE( 2,1)=0.71d0
      SITE( 2,2)=1.229755d0
      CHARGE( 2)=-0.0066d0
      SITE( 3,1)=-0.71d0
      SITE( 3,2)=1.229755d0
      CHARGE( 3)=-0.0066d0
      SITE( 4,1)=-1.42d0
      SITE( 4,2)=0.0d0
      CHARGE( 4)=-0.0066D0
      SITE( 5,1)=-0.71d0
      SITE( 5,2)=-1.229755d0
      CHARGE( 5)=-0.0066D0
      SITE( 6,1)=0.71d0
      SITE( 6,2)=-1.229755d0
      CHARGE( 6)=-0.0066D0
      SITE( 7,1)=2.838d0
      SITE( 7,2)=0.0d0
      CHARGE( 7)=0.0158d0
      SITE( 8,1)=1.419d0
      SITE( 8,2)=2.457779d0
      CHARGE( 8)=0.0158d0
      SITE( 9,1)=-1.419d0
      SITE( 9,2)=2.457779d0
      CHARGE( 9)=0.0158d0
      SITE(10,1)=-2.838d0
      SITE(10,2)=0.0d0
      CHARGE(10)=0.0158d0
      SITE(11,1)=-1.419d0
      SITE(11,2)=-2.45779d0
      CHARGE(11)=0.0158d0
      SITE(12,1)=1.419d0
      SITE(12,2)=-2.45779d0
      CHARGE(12)=0.0158d0
      SITE(13,1)=3.523366d0
      SITE(13,2)=1.248219d0
      CHARGE(13)=-0.0986d0
      SITE(14,1)=2.842673d0
      SITE(14,2)=2.427211d0
      CHARGE(14)=-0.0986d0
      SITE(15,1)=0.680706d0
      SITE(15,2)=3.675430d0
      CHARGE(15)=-0.0986d0
      SITE(16,1)=-0.680677d0
      SITE(16,2)=3.675437d0
      CHARGE(16)=-0.0986d0
      SITE(17,1)=-2.842673d0
      SITE(17,2)=2.427211d0
      CHARGE(17)=-0.0986D0
      SITE(18,1)=-3.523366D0
      SITE(18,2)=1.248219D0
      CHARGE(18)=-0.0986D0
      SITE(19,1)=-3.523366D0
      SITE(19,2)=-1.248219D0
      CHARGE(19)=-0.0986D0
      SITE(20,1)=-2.842673D0
      SITE(20,2)=-2.427211D0
      CHARGE(20)=-0.0986D0
      SITE(21,1)=-0.680677D0
      SITE(21,2)=-3.675437D0
      CHARGE(21)=-0.0986D0
      SITE(22,1)=0.680706D0
      SITE(22,2)=-3.675437D0
      CHARGE(22)=-0.0986D0
      SITE(23,1)=2.842673D0
      SITE(23,2)=-2.427211D0
      CHARGE(23)=-0.0986D0
      SITE(24,1)=3.523366D0
      SITE(24,2)=-1.248219D0
      CHARGE(24)=-0.0986D0
      SITE(25,1)=4.602048D0
      SITE(25,2)=1.274393D0
      CHARGE(25)=0.094D0
      SITE(26,1)=3.404682D0
      SITE(26,2)=3.348290D0
      CHARGE(26)=0.094D0
      SITE(27,1)=1.197383D0
      SITE(27,2)=4.622681D0
      CHARGE(27)=0.094D0
      SITE(28,1)=-1.197347D0
      SITE(28,2)=4.622693D0
      CHARGE(28)=0.094D0
      SITE(29,1)=-3.404682D0
      SITE(29,2)=3.348290D0
      CHARGE(29)=0.094D0
      SITE(30,1)=-4.602048D0
      SITE(30,2)=1.274393D0
      CHARGE(30)=0.094D0
      SITE(31,1)=-4.602048D0
      SITE(31,2)=-1.274393D0
      CHARGE(31)=0.094D0
      SITE(32,1)=-3.404682D0
      SITE(32,2)=-3.348290D0
      CHARGE(32)=0.094D0
      SITE(33,1)=-1.197347D0
      SITE(33,2)=-4.622693D0
      CHARGE(33)=0.094D0
      SITE(34,1)=1.197383D0
      SITE(34,2)=-4.622681D0
      CHARGE(34)=0.094D0
      SITE(35,1)=3.404682D0
      SITE(35,2)=-3.348290D0
      CHARGE(35)=0.094D0
      SITE(36,1)=4.602048D0
      SITE(36,2)=-1.274393D0
      CHARGE(36)=0.094D0

      C6cc=0.3926*12.075625**3 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12cc=0.3926*12.075625**6
      C6hh=0.0543*8.625969**3
      C12hh=0.0543*8.625969**6
      C6ch=0.1435*10.291264**3
      C12ch=0.1435*10.291264**6

      NAT2=NATOMS/2
      EPAH=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO

C
C  Potential energy first.
C
      DO J1=1,NAT2-1
         X1=X(3*(J1-1)+1)
         Y1=X(3*(J1-1)+2)
         Z1=X(3*(J1-1)+3)
         L1=X(3*(NAT2+J1-1)+1)
         M1=X(3*(NAT2+J1-1)+2)
         N1=X(3*(NAT2+J1-1)+3)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-20)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0010D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF

c        WRITE(*,'(A,6F15.5)') 'ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1=',ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1
         DO J2=J1+1,NAT2
            X2=X(3*(J2-1)+1)
            Y2=X(3*(J2-1)+2)
            Z2=X(3*(J2-1)+3)
            L2=X(3*(NAT2+J2-1)+1)
            M2=X(3*(NAT2+J2-1)+2)
            N2=X(3*(NAT2+J2-1)+3)
            L22=L2**2
            M22=M2**2
            N22=N2**2
            ALPHA2=SQRT(L22+M22+N22)
            RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-20)
            CA2=COS(ALPHA2)
            C2A2=CA2
            IF (ALPHA2.LT.0.00001D0) THEN
C              C3A2=-ALPHA2/2+ALPHA2**3/24
               C3A2=-0.5D0+ALPHA2**2/24.0D0
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  LJ and Coulomb contributions.
C
            DO K1=1,36
               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)
               DO K2=1,36
                  P2X=SITE(K2,1)
                  P2Y=SITE(K2,2)
                  P2Z=SITE(K2,3)
                  DIST= 
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
C           PRINT*,'coordinates of molecule pair:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
                  RDIST=1.0D0/DIST
                  RDIST6=RDIST**6
C
C  LJ term
C
                  IF (K1.LE.24.AND.K2.LE.24) THEN     ! C-C
                     DUMMY=FLJcc(RDIST6)
                  ELSEIF (K1.GT.24.AND.K2.GT.24) THEN ! H-H
                     DUMMY=FLJhh(RDIST6)
                  ELSE                                ! C-H
                     DUMMY=FLJch(RDIST6)
                  ENDIF
C                 PRINT*,'distance,DUMMY=',DIST,DUMMY
C
C  Coulomb term
C
                  DUMMY=DUMMY+HTOKJ*FC(CHARGE(K1),CHARGE(K2),RDIST)
                  IF (GTEST.OR.SECT) THEN
                     IF (K1.LE.24.AND.K2.LE.24) THEN     ! C-C
                        DUMMY2=DFLJcc(RDIST,RDIST6)
                     ELSEIF (K1.GT.24.AND.K2.GT.24) THEN ! H-H
                        DUMMY2=DFLJhh(RDIST,RDIST6)
                     ELSE                                ! C-H
                        DUMMY2=DFLJch(RDIST,RDIST6)
                     ENDIF
                     D1S=D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D2S=D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D3S=D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D4S=D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D5S=D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D6S=D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DAS=DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DBS=DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DCS=DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
C
C  Add Coul;omb contribution
C
                     DUMMY2=DUMMY2+HTOKJ*DFC(CHARGE(K1),CHARGE(K2),RDIST)
                     GX1=DUMMY2*D1S
                     GY1=DUMMY2*D2S
                     GZ1=DUMMY2*D3S
                     GL1=DUMMY2*D4S
                     GM1=DUMMY2*D5S
                     GN1=DUMMY2*D6S
                     GL2=DUMMY2*DAS
                     GM2=DUMMY2*DBS
                     GN2=DUMMY2*DCS

                     EPAH=EPAH+DUMMY
                     VT(J1)=VT(J1)+DUMMY
                     VT(J2)=VT(J2)+DUMMY

                     IF (GTEST) THEN
                        V(3*(J1-1)+1)=V(3*(J1-1)+1)+GX1
                        V(3*(J1-1)+2)=V(3*(J1-1)+2)+GY1
                        V(3*(J1-1)+3)=V(3*(J1-1)+3)+GZ1
                        V(3*(NAT2+J1-1)+1)=V(3*(NAT2+J1-1)+1)+GL1
                        V(3*(NAT2+J1-1)+2)=V(3*(NAT2+J1-1)+2)+GM1
                        V(3*(NAT2+J1-1)+3)=V(3*(NAT2+J1-1)+3)+GN1
                        V(3*(J2-1)+1)=V(3*(J2-1)+1)-GX1
                        V(3*(J2-1)+2)=V(3*(J2-1)+2)-GY1
                        V(3*(J2-1)+3)=V(3*(J2-1)+3)-GZ1
                        V(3*(NAT2+J2-1)+1)=V(3*(NAT2+J2-1)+1)+GL2
                        V(3*(NAT2+J2-1)+2)=V(3*(NAT2+J2-1)+2)+GM2
                        V(3*(NAT2+J2-1)+3)=V(3*(NAT2+J2-1)+3)+GN2
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     WRITE(*,'(A,G20.10)') 'energy=',EPAH
C     PRINT*,'coords:'
C     WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'gradient:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

      RETURN
      END
