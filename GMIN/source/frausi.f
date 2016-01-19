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
C   This program is supposed to calculate the energy, 
C   first derivatives and second derivatives for Silicon 
C   clusters using a Tight Binding (TB) approximation.  
C   This will then be used in conjunction with OPTIM,
C   as it is now known, to solve the problems of the universe
C   and in particular those of the silicon cluster world.
C   The parameters etc..are taken from PRB 47 12754 1993
C   and PRB 50 5645 1994 in case you were wondering.

      SUBROUTINE FRAUSI(N,XS,deriv1st,ENERGY,GTEST,ANGSTROM,natoms)
      use distance
      IMPLICIT NONE

      DOUBLE PRECISION EPSs,EPSp,C1, C2, C3, C4, C5,
     3                 C6, C7, C8, C9, C10,A1,B1,Ssssig,Sspsig,
     4                 Sppsig,Spppi,Hsssig,Hspsig,Hppsig,Hpppi,
     5                 C11,C12, Ang2Bohr, H2eV, HBohr2eVAng
 
      PARAMETER (EPSs=-0.3878175D0, EPSp=-0.166954D0,
     /     A1=1.0D0,B1=5.05D0,
     /     C1=20.9904791D0,C2=-18.6227459D0,
     /    C3=13.0417919D0, C4=-7.1832245D0, C5=3.0080981D0,
     /    C6=-0.8713337D0, C7=0.1321772D0, C8=-3.6D-06,
     /    C9=-4.5D-06, C10=1.71D-05,C11=-9.0D-06,C12=-2.17D-05,
     /    Ang2Bohr=1.889725989D0, H2eV=27.2113961D0, HBohr2eVAng=51.4220824D0)

      INTEGER N,natoms
      DOUBLE PRECISION Q(NATOMS,NATOMS),deriv1st(3*NATOMS) 


      LOGICAL GTEST,ANGSTROM
      INTEGER I, J, K, AI, AJ, NDIM, NMAX, I2, J2, K2,  INFO
      DOUBLE PRECISION DIST(NATOMS,NATOMS),
     1                 WORK(24*NATOMS),
     2                 S(4*NATOMS,4*NATOMS),
     3                 ENERGY,XS(3*NATOMS),
     4                 H(4*NATOMS,4*NATOMS)
      DOUBLE PRECISION UREP,REP
      LOGICAL FTEST
      COMMON /FAIL/ FTEST
      allocate(R(NATOMS,NATOMS),DIRCOS(NATOMS,NATOMS,3),Horig(4*NATOMS,4*NATOMS),DIAG(4*NATOMS))


      

C   Calculate distance matrix

      DO I=1,N
       I2=3*(I-1)
       R(I,I)=0.0D0
       DO J=I+1,N
            J2=3*(J-1)
            DIST(I,J)=(XS(I2+1)-XS(J2+1))**2+(XS(I2+2)-XS(J2+2))**2
     1                +(XS(I2+3)-XS(J2+3))**2
          R(I,J)=SQRT(DIST(I,J))

C   Now the direction cosines
C   Which incidentally are the projections of a unit vector onto 
C   the axes and NOT the cosines.
C   Note we have not calculated DIRCOS(I,I)

            DO K=1,3
               DIRCOS(I,J,K)=(XS(J2+K)-XS(I2+K))/R(I,J)
               DIRCOS(J,I,K)=-DIRCOS(I,J,K)
               DIRCOS(I,I,K)=0.0D0
            END DO

C   Convert distance given in Angstroms to bohr
            IF (ANGSTROM) THEN
             R(I,J)=Ang2Bohr*R(I,J)
            END IF
             R(J,I)=R(I,J)
         END DO
      END DO 


C   Calculate Repulsive Free Energy

      UREP=0.0D0
      DO I=1,N
         DO J=I+1,N

C CALCUALTIONS ARE IN BOHR
    
      Q(J,I)=(R(J,I)-((B1+A1)/2.0))/((B1-A1)/2.0)

      REP=C1 + C2*Q(J,I) + C3*(-1.0D0 + 2.0D0*Q(J,I)**2) + C4*(-3.0D0*Q(J,I) + 4.0D0*Q(J,I)**3) +
     /    C5*(1 -8.0D0*Q(J,I)**2 + 8.0D0*Q(J,I)**4) + C6*(5.0D0*Q(J,I) -20.0D0*Q(J,I)**3 + 16.0D0*Q(J,I)**5)
     /     + C7*(-1.0D0+18.0D0*Q(J,I)**2 -48.0D0*Q(J,I)**4 + 32.0D0*Q(J,I)**6)
     /   + C8*(-7.0D0*Q(J,I) + 56.0D0*Q(J,I)**3 - 112.0D0*Q(J,I)**5 + 64.0D0*Q(J,I)**7)
     /    + C9*(1.0D0 - 32.0D0*Q(J,I)**2 + 160.0D0*Q(J,I)**4 - 256.0D0*Q(J,I)**6
     /   +128.0D0*Q(J,I)**8) + C10*(9.0D0*Q(J,I) - 120.0D0*Q(J,I)**3 + 432.0D0*Q(J,I)**5
     /   - 576.0D0*Q(J,I)**7 + 256.0D0*Q(J,I)**9) 
     /    + C11*(-1.0D0 + 50.0D0*Q(J,I)**2 - 400.0D0*Q(J,I)**4 + 1120.0D0*Q(J,I)**6
     /    - 1280.0D0*Q(J,I)**8 + 512.0D0*Q(J,I)**10)
     /    + C12*(-11.0D0*Q(J,I) + 220.0D0*Q(J,I)**3 - 1232.0D0*Q(J,I)**5 + 2816.0D0*Q(J,I)**7
     /    - 2816.0D0*Q(J,I)**9 + 1024.0D0*Q(J,I)**11)
     /     - C1/2.0D0 

         IF (ABS(R(I,J)).GT.5.053328858168296D0) THEN
            REP=0.0D0
         ELSE
            REP=REP+0.00002119467638195263D0
         ENDIF

        UREP=REP+UREP

        ENDDO
      ENDDO

C10    FORMAT(A8, F20.17)

C   Only include overlap or interaction with different atoms here.
               
      DO I=1,N
         AI=4*(I-1)
C         AI2=2*AI
         DO J=I+1,N
            AJ=4*(J-1)

C   Dependence of V on distance between atoms R(I,J)
C   Similarly, S must scale with R(I,J)

C   Calculate the overlap matrix S(I,J), related to these parameters
C   by the Slater-Koster scheme

           CALL OVERL(N,I,J,Ssssig, Sspsig, Sppsig, Spppi,natoms)


            S((AI+1),(AJ+1))=Ssssig
            S((AI+1),(AJ+2))=DIRCOS(I,J,1)*Sspsig
            S((AI+1),(AJ+3))=DIRCOS(I,J,2)*Sspsig
            S((AI+1),(AJ+4))=DIRCOS(I,J,3)*Sspsig

            S((AI+2),(AJ+1))=-S((AI+1),(AJ+2))
            S((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,1)**2)*Spppi
            S((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Sppsig
     1                       -Spppi)
            S((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Sppsig
     1                       -Spppi)

            S((AI+3),(AJ+1))=-S((AI+1),(AJ+3))
            S((AI+3),(AJ+2))=S((AI+2),(AJ+3))
            S((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,2)**2)*Spppi
            S((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Sppsig
     1                       -Spppi)

            S((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Sspsig
            S((AI+4),(AJ+2))=S((AI+2),(AJ+4))
            S((AI+4),(AJ+3))=S((AI+3),(AJ+4))
            S((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,3)**2)*Spppi

            
C   Hamiltonian is found using the Slater-Koster scheme
C   as shown in their paper and Harrison's book p.481
C   which relates H to the parameter V using direction cosines.
      
           CALL HAMIL(N,I,J,Hsssig, Hspsig, Hppsig, Hpppi,natoms)


            H((AI+1),(AJ+1))=Hsssig
            H((AI+1),(AJ+2))=DIRCOS(I,J,1)*Hspsig
            H((AI+1),(AJ+3))=DIRCOS(I,J,2)*Hspsig
            H((AI+1),(AJ+4))=DIRCOS(I,J,3)*Hspsig

            H((AI+2),(AJ+1))=-H((AI+1),(AJ+2))
            H((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Hppsig+(1.0D0-
     1           DIRCOS(I,J,1)**2)*Hpppi
            H((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Hppsig
     1           -Hpppi)
            H((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Hppsig
     1           -Hpppi)

            H((AI+3),(AJ+1))=-H((AI+1),(AJ+3))
            H((AI+3),(AJ+2))=H((AI+2),(AJ+3))
            H((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Hppsig+(1.0D0-
     1           DIRCOS(I,J,2)**2)*Hpppi
            H((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Hppsig
     1           -Hpppi)

            H((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Hspsig
            H((AI+4),(AJ+2))=H((AI+2),(AJ+4))
            H((AI+4),(AJ+3))=H((AI+3),(AJ+4))
            H((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Hppsig+(1.0D0-
     1           DIRCOS(I,J,3)**2)*Hpppi

       END DO

C   Interaction/Overlap of different AOs on the same atom
         
          DO K=1,4
           DO K2=K+1,4
             S((AI+K),(AI+K2))=0.0D0
             H((AI+K),(AI+K2))=0.0D0
           END DO
          END DO

C   Interaction of same AO on the same atom
       H((AI+1),(AI+1))=EPSs
        H((AI+2),(AI+2))=EPSp
        H((AI+3),(AI+3))=EPSp
        H((AI+4),(AI+4))=EPSp

      END DO

C   Now we have to make sure the matrix S and H are symmetric.
C   Ergo, have to get the upper diagonal 'cos must have complete
C   matrix for Cholesky decomposition.
      
      DO I=1,4*N
        DO J=I+1,4*N
           S(J,I)=S(I,J)
           H(J,I)=H(I,J)
          END DO
    

C   Interaction/overlap of the same AO on the same atom
         S(I,I)=1.0D0

      END DO

C50    FORMAT (10(5X,20F8.5,5X))
C50     FORMAT (10F8.4)

      DO I=1,4*N
         DO J=1,4*N
            Horig(I,J)=H(I,J)
         END DO
      END DO


C      DO I=1,4*N
C        DO J=1,4*N
C           Sorig(I,J)=S(I,J)
C        END DO
C      END DO

      NDIM=4*N
      NMAX=4*NATOMS

      FTEST=.FALSE.

      CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, 24*NATOMS, INFO )
C
C  This sort isn't needed for the SGI lapack.f routine which
C  incorporates it.
C
      CALL EIGSRT(DIAG,H,NDIM,NMAX)

 
      IF (INFO.NE.0) THEN
         FTEST=.TRUE.
         ENERGY=1.0D6
         PRINT*,'DSYGV failed - INFO=',INFO
         RETURN
      ENDIF

C      WRITE(6,*)
C      WRITE(6,*)  'The eigenvalues of the Hessian matrix:'
C      DO I=1,4*N
C         WRITE(6,*)  DIAG(I)
C      END DO
C      WRITE(6,*)

      ENERGY=0.0D0

      DO J=4*N,(2*N+1),-1
         ENERGY=ENERGY+DIAG(J)
      END DO
C     WRITE(*,'(A,2F20.10)') 'HOMO and LUMO at ',DIAG(2*N+1),DIAG(2*N)

C   For single atom there should be no contribution from the
C   bond counting term
C     WRITE(*,'(A,2F20.10)') 'ENERGY and UREP=',2*ENERGY,UREP
      IF (N .EQ. 1) THEN
         ENERGY=2*ENERGY
      ELSE
          ENERGY=2*ENERGY+UREP
      END IF
C      WRITE(6,*)
C      WRITE(6,*) 'ENERGY=',ENERGY

C   Now call the all new originally crafted subroutine to calculate
C   the first derivatives.
C   Convert R(I,J) from Bohr to Angstrom and ENERGY from hartree to eV
C   derivs1st should NOW be in eV/Angstrom as returned from subroutine.

      IF (GTEST) THEN
          CALL FRAUSIDER(N,deriv1st,H,ANGSTROM,natoms)
      END IF
      IF (ANGSTROM) THEN
      DO I=1,N
         DO J=1,N
             R(I,J)=R(I,J)/Ang2Bohr
         END DO
      END DO
      ENERGY=ENERGY*H2eV
      END IF

      RETURN
      END

C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE OVERL(N,K,L,SSSs,SSPs,SPPs,SPPp,natoms)
      use distance
      IMPLICIT NONE
      INTEGER  K,L,natoms

      INTEGER N
      DOUBLE PRECISION Q(NATOMS,NATOMS)

      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     /     D1,D2,D3,
     7     D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,
     8     D21,D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36,
     9     D37,D38,D39,D40,A3,B3,D41,D42,D43,D44,D45,D46,D47,D48

      PARAMETER ( D1=0.3942975, D2=-0.3092186, D3=0.1395490,
     /     D4=-0.0176746, D5=-0.0169834,
     /    D6=8.3055D-03,D7=1.4080D-03,D8=-2.9477D-03,D9=1.3781D-03,
     /    D10=-2.380D-04, D11=-7.16D-05,D12=5.26D-05,
     /    D13=-0.3568323,D14=0.2467967,D15=-0.0456893,
     /    D16=-0.0605347,D17=0.0497093,D18=-0.0102684,
     /    D19=-6.4619D-03,D20=5.4908D-03,D21=-1.6044D-03,
     /    D22=-1.030D-04,D23=2.934D-04,D24=-1.581D-04,
     /    D25=-0.1400088,D26=0.0192343,D27=0.1646371,
     /    D28=-0.181142,D29=0.0754221,D30=3.0624D-03,
     /    D31=-0.0183958,D32=8.7078D-03,D33=-1.1229D-03,
     /    D34=-8.468D-04,D35=5.924D-04,D36=-1.579D-04,
     /    D37=0.3722275,D38=-0.3063175,D39=0.1654376,
     /    D40=-0.0484825,D41=-2.4093D-03,D42=9.0576D-03,
     /    D43=-3.7347D-03,D44=3.162D-04,D45=3.926D-04,
     /    D46=-2.436D-04,D47=7.48D-05,D48=1.01D-05)

       A3=1.5D0
       B3=9.5D0
       Q(K,L)=(R(K,L)-((B3+A3)/2.0D0))/((B3-A3)/2.0D0)

C       PRINT*,'r of ij in OVERL is',R(K,L)

C  Ssssig
      SSSs=D1 + D2*Q(K,L) + D3*(-1.0D0 + 2.0D0*Q(K,L)**2) 
     /    + D4*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    D5*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4) 
     /    + D6*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + D7*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + D8*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + D9*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + D10*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +D11*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +D12*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - D1/2.0D0

      IF (R(K,L).GT.9.489340860704047D0) THEN
         SSSs=0.0D0
      ELSE
         SSSs=SSSs-0.0007078111506558793D0
      ENDIF

C   Sspsig
      SSPs=D13 + D14*Q(K,L) + D15*(-1.0D0 + 2.0D0*Q(K,L)**2)
     /    + D16*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    D17*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + D18*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + D19*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + D20*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + D21*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) 
     /    + D22*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9) 
     /   +D23*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +D24*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - D13/2.0D0

      IF (R(K,L).GT.9.360442319381175D0) THEN
         SSPs=0.0D0
      ELSE
         SSPs=SSPs+0.0008518797271876711D0
      ENDIF

C Sppsig
      SPPs=D25 + D26*Q(K,L) + D27*(-1.0D0 + 2.0D0*Q(K,L)**2)        
     /    + D28*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    D29*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + D30*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + D31*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + D32*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + D33*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + D34*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +D35*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +D36*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /      - D25/2.0D0

      IF (R(K,L).GT.9.096495577461203D0) THEN
         SPPs=0.0D0
      ELSE
         SPPs=SPPs-0.0002134373553017432D0
      ENDIF

C   Spppi
      SPPp=D37 + D38*Q(K,L) + D39*(-1.0D0 + 2.0D0*Q(K,L)**2)        
     /    + D40*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    D41*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + D42*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + D43*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + D44*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + D45*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + D46*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +D47*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +D48*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - D37/2.0D0

      IF (R(K,L).GT.9.358643583694771D0) THEN
         SPPp=0.0D0
      ELSE
         SPPp=SPPp-0.0001838066680916461D0
      ENDIF

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE HAMIL(N,K,L,SSSs,SSPs,SPPs,SPPp,natoms)
      use distance
      IMPLICIT NONE
      INTEGER  K,L,natoms

      INTEGER N
      DOUBLE PRECISION Q(NATOMS,NATOMS)

      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     /     G1,G2,G3,
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2,G41,G42,G43,G44,G45,G46,G47,G48

C    /    G13=0.1777735,G14=-0.1094317,G15=-0.0071041,
C  Needed to adjust G14 since the function doesn't actually have a turning
C  point around the cutoff region otherwise!
C
C
C    /    G25=0.0510829,G26=0.0092781,G27=-0.0894841,
C  Needed to adjust G26 since the function doesn't actually have a turning
C  point around the cutoff region otherwise!
C
       PARAMETER (G1=-0.2601372,G2=0.1958030,G3=-0.0716540,
     /    G4=-0.0084811,G5=0.0229926,G6=-9.8866D-03,
     /    G7=-8.281D-04,G8=3.5950D-03,G9=-2.6770D-03,
     /    G10=1.3587D-03,G11=-5.617D-04,G12=2.165D-04,
     /    G13=0.1777735,G14=-0.1082,G15=-0.0071041,
     /    G16=0.0557698,G17=-0.0376445,G18=0.0088910,
     /    G19=0.0041357,G20=-0.0052914,G21=0.0031425,
     /    G22=-0.0014231,G23=5.843D-04,G24=-2.103D-04,
     /    G25=0.0510829,G26=0.0145,G27=-0.0894841,
     /    G28=0.0856069,G29=-0.0355870,G30=1.3150D-03,
     /    G31=7.8269D-03,G32=-6.1665D-03,G33=3.2455D-03,
     /    G34=-1.4904D-03,G35=6.809D-04,G36=-2.655D-04,
     /    G37=-0.1737379,G38=0.1403235,G39=-0.0716452,
     /    G40=0.0185100,G41=0.0027857, G42=-5.0867D-03,
     /    G43=2.5525D-03, G44=-6.749D-04, G45=-2.12D-05,
     /    G46=1.537D-04, G47=-1.269D-04, G48=7.84D-05 )

       A2=1.5D0
       B2=9.5D0
       Q(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C  Ssssig
      SSSs=G1 + G2*Q(K,L) + G3*(-1.0D0 + 2.0D0*Q(K,L)**2) 
     /    + G4*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    G5*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4) 
     /    + G6*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + G7*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + G8*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + G9*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + G10*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +G11*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +G12*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - G1/2.0D0

         IF (R(K,L).GT.9.074040365621258D0) THEN
            SSSs=0.0D0
         ELSE
            SSSs=SSSs+0.0002808028388218697D0
         ENDIF

C   Sspsig
      SSPs=G13 + G14*Q(K,L) + G15*(-1.0D0 + 2.0D0*Q(K,L)**2)
     /    + G16*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    G17*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + G18*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + G19*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + G20*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + G21*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) 
     /    + G22*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9) 
     /   +G23*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +G24*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - G13/2.0D0
       IF (R(K,L).GT.9.246570341781329D0) THEN
          SSPs=0.0D0
       ELSE
          SSPs=SSPs-0.00164172705840912D0
       ENDIF

C Sppsig
      SPPs=G25 + G26*Q(K,L) + G27*(-1.0D0 + 2.0D0*Q(K,L)**2)
     /    + G28*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    G29*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + G30*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + G31*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + G32*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + G33*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + G34*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +G35*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +G36*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /      - G25/2.0D0
       IF (R(K,L).GT.9.273572893686612D0) THEN
          SPPs=0.0D0
       ELSE
          SPPs=SPPs-0.005828116913694732D0
       ENDIF


C   Spppi
      SPPp=G37 + G38*Q(K,L) + G39*(-1.0D0 + 2.0D0*Q(K,L)**2)        
     /    + G40*(-3.0D0*Q(K,L) + 4.0D0*Q(K,L)**3) +
     /    G41*(1 -8.0D0*Q(K,L)**2 + 8.0D0*Q(K,L)**4)
     /    + G42*(5.0D0*Q(K,L) -20.0D0*Q(K,L)**3 + 16.0D0*Q(K,L)**5)
     /     + G43*(-1.0D0+18.0D0*Q(K,L)**2 -48.0D0*Q(K,L)**4 + 32.0D0*Q(K,L)**6)
     /   + G44*(-7.0D0*Q(K,L) + 56.0D0*Q(K,L)**3 - 112.0D0*Q(K,L)**5 + 64.0D0*Q(K,L)**7)
     /    + G45*(1.0D0 - 32.0D0*Q(K,L)**2 + 160.0D0*Q(K,L)**4 - 256.0D0*Q(K,L)**6
     /   +128.0D0*Q(K,L)**8) + G46*(9.0D0*Q(K,L) - 120.0D0*Q(K,L)**3 + 432.0D0*Q(K,L)**5
     /   - 576.0D0*Q(K,L)**7 + 256.0D0*Q(K,L)**9)
     /   +G47*(-1.0D0 + 50.0D0*Q(K,L)**2 - 400.0D0*Q(K,L)**4 + 1120.0D0*Q(K,L)**6
     /    - 1280.0D0*Q(K,L)**8 + 512.0D0*Q(K,L)**10)
     /   +G48*(-11.0D0*Q(K,L) + 220.0D0*Q(K,L)**3 - 1232.0D0*Q(K,L)**5 + 2816.0D0*Q(K,L)**7
     /    - 2816.0D0*Q(K,L)**9 + 1024.0D0*Q(K,L)**11)
     /    - G37/2.0D0

       IF (R(K,L).GT.9.10748558948335D0) THEN
          SPPp=0.0D0
       ELSE
          SPPp=SPPp+0.00009727014254751198D0
       ENDIF

      RETURN
      END
C   This is phase two where we calculate the first derivatives
C   of the silicon potential coded in TB.Si.f. It is in tight
C   binding parameterised form.
C   Off I go...
C   This is derivs1.f.jloop.betterish in the resrves directory
C   for the record

      SUBROUTINE FRAUSIDER(N,deriv1st,EIGVEC,ANGSTROM,natoms)
      use distance
      IMPLICIT NONE 

      DOUBLE PRECISION EPSs,EPSp,C1, C2, C3, C4, C5,
     3                 C6, C7, C8, C9, C10,A1,B1,Ssssig,Sspsig,
     4                 Sppsig,Spppi,Hsssig,Hspsig,Hppsig,Hpppi,
     5                 C11,C12, Ang2Bohr, H2eV, HBohr2eVAng
 
      PARAMETER (EPSs=-0.3878175D0, EPSp=-0.166954D0,
     /     A1=1.0D0,B1=5.05D0,
     /     C1=20.9904791D0,C2=-18.6227459D0,
     /    C3=13.0417919D0, C4=-7.1832245D0, C5=3.0080981D0,
     /    C6=-0.8713337D0, C7=0.1321772D0, C8=-3.6D-06,
     /    C9=-4.5D-06, C10=1.71D-05,C11=-9.0D-06,C12=-2.17D-05,
     /    Ang2Bohr=1.889725989D0, H2eV=27.2113961D0, HBohr2eVAng=51.4220824D0)

      INTEGER N,natoms
      DOUBLE PRECISION Q(NATOMS,NATOMS),deriv1st(3*NATOMS)

      INTEGER I, J, K , AI, AJ, LEV, NUM1, NUM2, 
     1        I2,J2
      DOUBLE PRECISION REP1st(NATOMS,3),RR, D1, D2, D3,
     1                 ELEC1st(NATOMS,3), E1, E2, DUMMY, E3,
     2                 diffSV(2,4,4), SDC1, SQC1, SQC2, SQC3,
     3                 R1st(NATOMS,NATOMS),EIGVEC(4*NATOMS,4*NATOMS),
     4                 UREP1st
      DOUBLE PRECISION SVSUB,
     1                 SQDIRCOS(3*NATOMS,NATOMS),SUBDIRCOS(3*NATOMS,NATOMS),
     2                 DELD1, DELD12, DELD2, DELD22, DELD3, DELD32, DELRDX,DRSss,DRSsp,
     3                 DRSpps,DRSppp,DQDR,DSVSUB
      LOGICAL ANGSTROM


C
C  Take the transpose of EIGVEC - makes loops more efficient
C
      DO J=1,4*N
         DO I=J+1,4*N
            DUMMY=EIGVEC(J,I)
            EIGVEC(J,I)=EIGVEC(I,J)
            EIGVEC(I,J)=DUMMY
         ENDDO
      ENDDO
C
C   Can put arrays into a block data format instead of parameter
C   so can just put in once in a file to be included

C   Can't start an array with a no. ie. 1stderiv...syntax error

C   First of all calculate derivatives of repulsive energy term
C   This is called REP1st.

      DQDR=1.0D0/((B1-A1)/2.0D0)

      DO K=1,3
         DO I=1,N
            DO J=I+1,N

      Q(J,I)=(R(J,I)-((B1+A1)/2.0))/((B1-A1)/2.0)

       
         UREP1st=C2*DQDR + C3*4.0D0*Q(J,I)*DQDR
     /   +C4*(-3.0D0*DQDR + 12.0D0*(Q(J,I)**2)*DQDR)
     /   + C5*(-16.0D0*Q(J,I)*DQDR + 32.0D0*(Q(J,I)**3)*DQDR) +C6*(5.0D0*DQDR -60.0D0*(Q(J,I)**2)*DQDR
     /   + 80.0D0*(Q(J,I)**4)*DQDR) +C7*(36.0D0*Q(J,I)*DQDR - 192.0D0*(Q(J,I)**3)*DQDR +192.0D0*(Q(J,I)**5)*DQDR)
     /   +C8*(-7.0D0*DQDR +168.0D0*(Q(J,I)**2)*DQDR - 560.0D0*(Q(J,I)**4)*DQDR +448.0D0*(Q(J,I)**6)*DQDR)
     /   +C9*(-64.0D0*Q(J,I)*DQDR +640.0D0*(Q(J,I)**3)*DQDR -1536.0D0*(Q(J,I)**5)*DQDR +1024.0D0*(Q(J,I)**7)*DQDR)
     /   + C10*(9.0D0*DQDR -360.0D0*(Q(J,I)**2)*DQDR +2160*(Q(J,I)**4)*DQDR
     /   -4032.0D0*(Q(J,I)**6)*DQDR +2304.0D0*(Q(J,I)**8)*DQDR)
     /   +C11*(100.0D0*Q(J,I)*DQDR - 1600*(Q(J,I)**3)*DQDR + 6720*(Q(J,I)**5)*DQDR - 10240*(Q(J,I)**7)*DQDR
     /       +5120*(Q(J,I)**9)*DQDR)
     /   +C12*(-11.0D0*DQDR + 660.0D0*(Q(J,I)**2)*DQDR - 6160*(Q(J,I)**4)*DQDR +19712*(Q(J,I)**6)*DQDR
     /   -25344*(Q(J,I)**8)*DQDR + 11264*(Q(J,I)**10)*DQDR)


C           IF (ABS(R(J,I)).GT.B1) UREP1st=0.0D0
            IF (ABS(R(I,J)).GT.5.053328858168296D0) THEN
               UREP1st=0.0D0
            ENDIF

               R1st(J,I)=-DIRCOS(I,J,K)*UREP1st
               R1st(I,J)=-R1st(J,I)

            END DO

            R1st(I,I)=0.0D0
            D1=0.0D0
            DO J=1,N
              D1=D1+R1st(J,I)
            END DO

            REP1st(I,K)=D1
         END DO
      END DO

C   CO=1,2 or 3 represents the x,y or z coordinate and thus differentiation
C   that coordinate. Cripes, I hope my work is more compact than my comment
C   statements.
      DO I=1,N
        DO J=I+1,N
           J2=3*(J-1)

               DO NUM1=1,3
                  SQDIRCOS(J2+NUM1,I)=DIRCOS(J,I,NUM1)**2
                  SUBDIRCOS(J2+NUM1,I)=SQDIRCOS(J2+NUM1,I)-1.0D0
               END DO

        END DO
        END DO



C        PRINT*,'Entering the differentiation of electronic energy zone'
        DO K=1,N
          ELEC1st(K,1)=0.0D0
        ENDDO

        DO I=1,N
         AI=4*(I-1)

C   Diff of interactions of AOs on same atom
C   Don't need to include for diffSV since not used - loop
C   is only for J=I+1,N

         DO J=I+1,N
           RR=1.0D0/R(J,I)
           AJ=4*(J-1)
           J2=3*(J-1)

C   A good idea appears to be simplify the expressions in the
C   differentiating section and reduce the number of indices


           D1=DIRCOS(J,I,1)
           D2=DIRCOS(J,I,2)
           D3=DIRCOS(J,I,3)

           SDC1=SUBDIRCOS(J2+1,I)
C           SDC2=SUBDIRCOS(J2+2,I)
C           SDC3=SUBDIRCOS(J2+3,I)
           SQC1=SQDIRCOS(J2+1,I)
           SQC2=SQDIRCOS(J2+2,I)
           SQC3=SQDIRCOS(J2+3,I)


C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

C               DO SV=1,2
              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

               
C           IF (SV.EQ.1) THEN

           CALL OVERL(N,I,J,Ssssig, Sspsig, Sppsig, Spppi,natoms)

           CALL DEROV(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)

             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp

                  diffSV(1,1,1)=DRSss*DELRDX
                  diffSV(1,1,2)=Sspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(1,1,3)=Sspsig*DELD2 - D2*DRSsp*DELRDX
                  diffSV(1,1,4)=Sspsig*DELD3 - D3*DRSsp*DELRDX


                  diffSV(1,2,1)=-diffSV(1,1,2)
                  diffSV(1,2,2)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
                  diffSV(1,2,3)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
                  diffSV(1,2,4)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              
                  diffSV(1,3,1)=-diffSV(1,1,3)
                  diffSV(1,3,2)= diffSV(1,2,3)
                  diffSV(1,3,3)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
                  diffSV(1,3,4)= D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

                  diffSV(1,4,1)=-diffSV(1,1,4)
                  diffSV(1,4,2)= diffSV(1,2,4)
                  diffSV(1,4,3)= diffSV(1,3,4)
                  diffSV(1,4,4)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ELSE

               CALL HAMIL(N,I,J,Hsssig, Hspsig, Hppsig, Hpppi,natoms)


               CALL DERHAM(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)

             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(2,1,1)=DRSss*DELRDX
              diffSV(2,1,2)=Hspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(2,1,3)=Hspsig*DELD2 - D2*DRSsp*DELRDX
              diffSV(2,1,4)=Hspsig*DELD3 - D3*DRSsp*DELRDX

              diffSV(2,2,1)=-diffSV(2,1,2)
              diffSV(2,2,2)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(2,2,3)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSV(2,2,4)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(2,3,1)=-diffSV(2,1,3)
              diffSV(2,3,2)= diffSV(2,2,3)
              diffSV(2,3,3)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(2,3,4)= D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSV(2,4,1)=-diffSV(2,1,4)
              diffSV(2,4,2)= diffSV(2,2,4)
              diffSV(2,4,3)= diffSV(2,3,4)
              diffSV(2,4,4)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ENDIF

C               END DO ! SV

           E1=0.0D0


           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)
              E3=diffSV(2,NUM1,NUM2)


              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(E3-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,1)=ELEC1st(I,1)+2.0D0*E1
           ELEC1st(J,1)=ELEC1st(J,1)-2.0D0*E1

         END DO  !J (1st time)
        END DO  !I(1st time)

        DO K=1,N
          ELEC1st(K,2)=0.0D0
        ENDDO


        DO I=1,N
         AI=4*(I-1)

C   Diff of interactions of AOs on same atom
C   Don't need to include for diffSV since not used - loop
C   is only for J=I+1,N

         DO J=I+1,N

           RR=1.0D0/R(J,I)
           AJ=4*(J-1)
           J2=3*(J-1)

           D1=DIRCOS(J,I,2)
           D2=DIRCOS(J,I,3)
           D3=DIRCOS(J,I,1)

           SDC1=SUBDIRCOS(J2+2,I)
C           SDC2=SUBDIRCOS(J2+3,I)
C           SDC3=SUBDIRCOS(J2+1,I)
           SQC1=SQDIRCOS(J2+2,I)
           SQC2=SQDIRCOS(J2+3,I)
           SQC3=SQDIRCOS(J2+1,I)

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

C          PRINT*,'YTERMS'

              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

C               DO SV=1,2
C           IF (SV.EQ.1) THEN

           CALL OVERL(N,I,J,Ssssig, Sspsig, Sppsig, Spppi,natoms)

           CALL DEROV(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)

             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp


              diffSV(1,1,1)=DRSss*DELRDX
              diffSV(1,1,3)=Sspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(1,1,4)=Sspsig*DELD2 - D2*DRSsp*DELRDX
              diffSV(1,1,2)=Sspsig*DELD3 - D3*DRSsp*DELRDX

              diffSV(1,2,1)=-diffSV(1,1,2)
              diffSV(1,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(1,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSV(1,2,4)=D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2

              diffSV(1,3,1)=-diffSV(1,1,3)
              diffSV(1,3,2)= diffSV(1,2,3)
              diffSV(1,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(1,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSV(1,4,1)=-diffSV(1,1,4)
              diffSV(1,4,2)= diffSV(1,2,4)
              diffSV(1,4,3)= diffSV(1,3,4)
              diffSV(1,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ELSE

               CALL HAMIL(N,I,J,Hsssig, Hspsig, Hppsig, Hpppi,natoms)

               CALL DERHAM(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)


             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(2,1,1)=DRSss*DELRDX
              diffSV(2,1,3)=Hspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(2,1,4)=Hspsig*DELD2 - D2*DRSsp*DELRDX
              diffSV(2,1,2)=Hspsig*DELD3 - D3*DRSsp*DELRDX

              diffSV(2,2,1)=-diffSV(2,1,2)
              diffSV(2,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(2,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSV(2,2,4)=D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2

              diffSV(2,3,1)=-diffSV(2,1,3)
              diffSV(2,3,2)= diffSV(2,2,3)
              diffSV(2,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(2,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSV(2,4,1)=-diffSV(2,1,4)
              diffSV(2,4,2)= diffSV(2,2,4)
              diffSV(2,4,3)= diffSV(2,3,4)
              diffSV(2,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ENDIF

C               END DO ! SV

           E1=0.0D0

           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)
              E3=diffSV(2,NUM1,NUM2)

              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(E3-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,2)=ELEC1st(I,2)+2.0D0*E1
           ELEC1st(J,2)=ELEC1st(J,2)-2.0D0*E1

         END DO  !J (1st time)
        END DO  !I(1st time)

        DO K=1,N
          ELEC1st(K,3)=0.0D0
        ENDDO
        DO I=1,N
         AI=4*(I-1)

         DO J=I+1,N

           RR=1.0D0/R(J,I)
           AJ=4*(J-1)
           J2=3*(J-1)

           D1=DIRCOS(J,I,3)
           D2=DIRCOS(J,I,1)
           D3=DIRCOS(J,I,2)


           SDC1=SUBDIRCOS(J2+3,I)
C           SDC2=SUBDIRCOS(J2+1,I)
C           SDC3=SUBDIRCOS(J2+2,I)
           SQC1=SQDIRCOS(J2+3,I)
           SQC2=SQDIRCOS(J2+1,I)
           SQC3=SQDIRCOS(J2+2,I)


C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

C               DO SV=1,2
              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

           CALL OVERL(N,I,J,Ssssig, Sspsig, Sppsig, Spppi,natoms)


C              IF (SV.EQ.1) THEN

               CALL DEROV(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)

             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp


              diffSV(1,1,1)=DRSss*DELRDX
              diffSV(1,1,4)=Sspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(1,1,2)=Sspsig*DELD2 - D2*DRSsp*DELRDX
              diffSV(1,1,3)=Sspsig*DELD3 - D3*DRSsp*DELRDX

              diffSV(1,2,1)=-diffSV(1,1,2)
              diffSV(1,4,4)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(1,2,3)=D3*D2*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2
              diffSV(1,2,4)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSV(1,3,1)=-diffSV(1,1,3)
              diffSV(1,3,2)= diffSV(1,2,3)
              diffSV(1,2,2)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(1,3,4)= D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(1,4,1)=-diffSV(1,1,4)
              diffSV(1,4,2)= diffSV(1,2,4)
              diffSV(1,4,3)= diffSV(1,3,4)
              diffSV(1,3,3)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ELSE

               CALL HAMIL(N,I,J,Hsssig, Hspsig, Hppsig, Hpppi,natoms)


               CALL DERHAM(N,I,J,DRSss,DRSsp,DRSpps,DRSppp,natoms)

             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(2,1,1)=DRSss*DELRDX
              diffSV(2,1,4)=Hspsig*DELD1 - D1*DRSsp*DELRDX
              diffSV(2,1,2)=Hspsig*DELD2 - D2*DRSsp*DELRDX
              diffSV(2,1,3)=Hspsig*DELD3 - D3*DRSsp*DELRDX

              diffSV(2,2,1)=-diffSV(2,1,2)
              diffSV(2,4,4)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(2,2,3)=D3*D2*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2
              diffSV(2,2,4)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSV(2,3,1)=-diffSV(2,1,3)
              diffSV(2,3,2)= diffSV(2,2,3)
              diffSV(2,2,2)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(2,3,4)= D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(2,4,1)=-diffSV(2,1,4)
              diffSV(2,4,2)= diffSV(2,2,4)
              diffSV(2,4,3)= diffSV(2,3,4)
              diffSV(2,3,3)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

C              ENDIF
C               END DO ! SV

          E1=0.0D0

           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)
              E3=diffSV(2,NUM1,NUM2)

              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(E3-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,3)=ELEC1st(I,3)+2.0D0*E1
           ELEC1st(J,3)=ELEC1st(J,3)-2.0D0*E1


         END DO  !J (1st time)
        END DO  !I(1st time)

C   In order to differentiate by the coordinate of atom K, you only get 
C   a value for the derivative if one of the two atoms involved in the
C   pair interaction is K otherwise it is zero. 
C   Have to keep it in a loop where I and J are the atom number rather
C   than the AO.

C   Must be aware that although none of the diagonal elements are defined
C   they are equal to zero. Not deemed necessary to calc. them since
C   the loop never reaches I=J.  Good to use J=I+1 since uses symmetry to 
C   cut down on calculations.

C   Include factor of 2 for ELEC1st since two electrons fill each MO
C   Include Hbohr2eVAng conversion factor to return in units of eV/Angstrom 
C   that I know and love

      DO K=1,3
         DO I=1,N
            I2=3*(I-1)
            IF (ANGSTROM) THEN
             deriv1st(I2+K)=(2.0D0*ELEC1st(I,K)+REP1st(I,K))*HBohr2eVAng
            ELSE
             deriv1st(I2+K)=(2.0D0*ELEC1st(I,K)+REP1st(I,K))
            END IF
         END DO
      END DO
       
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE DEROV(N,K,L,SSSs,SSPs,SPPs,SPPp,natoms)
      use distance
      IMPLICIT NONE
      INTEGER  K,L,natoms

      INTEGER N
      DOUBLE PRECISION Q(NATOMS,NATOMS)

      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,DQDR,
     /     D1,D2,D3,
     7     D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,
     8     D21,D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36,
     9     D37,D38,D39,D40,A3,B3,D41,D42,D43,D44,D45,D46,D47,D48

      PARAMETER ( D1=0.3942975, D2=-0.3092186, D3=0.1395490,
     /     D4=-0.0176746, D5=-0.0169834,
     /    D6=8.3055D-03,D7=1.4080D-03,D8=-2.9477D-03,D9=1.3781D-03,
     /    D10=-2.380D-04, D11=-7.16D-05,D12=5.26D-05,
     /    D13=-0.3568323,D14=0.2467967,D15=-0.0456893,
     /    D16=-0.0605347,D17=0.0497093,D18=-0.0102684,
     /    D19=-6.4619D-03,D20=5.4908D-03,D21=-1.6044D-03,
     /    D22=-1.030D-04,D23=2.934D-04,D24=-1.581D-04,
     /    D25=-0.1400088,D26=0.0192343,D27=0.1646371,
     /    D28=-0.181142,D29=0.0754221,D30=3.0624D-03,
     /    D31=-0.0183958,D32=8.7078D-03,D33=-1.1229D-03,
     /    D34=-8.468D-04,D35=5.924D-04,D36=-1.579D-04,
     /    D37=0.3722275,D38=-0.3063175,D39=0.1654376,
     /    D40=-0.0484825,D41=-2.4093D-03,D42=9.0576D-03,
     /    D43=-3.7347D-03,D44=3.162D-04,D45=3.926D-04,
     /    D46=-2.436D-04,D47=7.48D-05,D48=1.01D-05)

       A3=1.5D0
       B3=9.5D0

       Q(K,L)=(R(K,L)-((B3+A3)/2.0D0))/((B3-A3)/2.0D0)

C         PRINT*,'r IN DEROV is',R(K,L)

      DQDR=1.0D0 /((B3-A3)/2.0D0)

C  Ssssig

       SSSs=D2*DQDR + D3*4.0D0*Q(K,L)*DQDR
     /   +D4*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + D5*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +D6*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +D7*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +D8*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +D9*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +D10*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +D11*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +D12*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

      IF (R(K,L).GT.9.489340860704047D0) SSSs=0.0D0

C   Sspsig
      SSPs=D14*DQDR + D15*4.0D0*Q(K,L)*DQDR
     /   +D16*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + D17*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +D18*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +D19*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +D20*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +D21*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +D22*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +D23*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +D24*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

      IF (R(K,L).GT.9.360442319381175D0) SSPs=0.0D0

C Sppsig
      SPPs=D26*DQDR + D27*4.0D0*Q(K,L)*DQDR
     /   +D28*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + D29*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +D30*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +D31*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +D32*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +D33*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +D34*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +D35*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +D36*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)


      IF (R(K,L).GT.9.096495577461203D0) SPPs=0.0D0

C   Spppi
      SPPp=D38*DQDR + D39*4.0D0*Q(K,L)*DQDR
     /   +D40*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + D41*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +D42*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +D43*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +D44*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +D45*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +D46*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +D47*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +D48*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

      IF (R(K,L).GT.9.358643583694771) SPPp=0.0D0

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   


      SUBROUTINE DERHAM(N,K,L,SSSs,SSPs,SPPs,SPPp,natoms)
      use distance
      IMPLICIT NONE
      INTEGER  K,L,natoms

      INTEGER N
      DOUBLE PRECISION Q(NATOMS,NATOMS)

      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,DQDR,
     /     G1,G2,G3,
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2,G41,G42,G43,G44,G45,G46,G47,G48

C
C N.B. G14 and G26 have been altered to give functions with turning points
C near the cutoff.
C
       PARAMETER (G1=-0.2601372,G2=0.1958030,G3=-0.0716540,
     /    G4=-0.0084811,G5=0.0229926,G6=-9.8866D-03,
     /    G7=-8.281D-04,G8=3.5950D-03,G9=-2.6770D-03,
     /    G10=1.3587D-03,G11=-5.617D-04,G12=2.165D-04,
     /    G13=0.1777735,G14=-0.1082D0 ,G15=-0.0071041,
     /    G16=0.0557698,G17=-0.0376445,G18=0.0088910,
     /    G19=0.0041357,G20=-0.0052914,G21=0.0031425,
     /    G22=-0.0014231,G23=5.843D-04,G24=-2.103D-04,
     /    G25=0.0510829,G26=0.0145D0,G27=-0.0894841,
     /    G28=0.0856069,G29=-0.0355870,G30=1.3150D-03,
     /    G31=7.8269D-03,G32=-6.1665D-03,G33=3.2455D-03,
     /    G34=-1.4904D-03,G35=6.809D-04,G36=-2.655D-04,
     /    G37=-0.1737379,G38=0.1403235,G39=-0.0716452,
     /    G40=0.0185100,G41=0.0027857, G42=-5.0867D-03,
     /    G43=2.5525D-03, G44=-6.749D-04, G45=-2.12D-05,
     /    G46=1.537D-04, G47=-1.269D-04, G48=7.84D-05 )

       A2=1.5D0
       B2=9.5D0

       Q(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

      DQDR=1.0D0/((B2-A2)/2.0D0)

C  Ssssig

       SSSs=G2*DQDR + G3*4.0D0*Q(K,L)*DQDR
     /   +G4*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + G5*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +G6*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +G7*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +G8*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +G9*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +G10*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +G11*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +G12*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

       IF (R(K,L).GT.9.074040365621258D0) SSSs=0.0D0

C   Sspsig
      SSPs=G14*DQDR + G15*4.0D0*Q(K,L)*DQDR
     /   +G16*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + G17*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +G18*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +G19*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +G20*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +G21*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +G22*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +G23*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +G24*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

       IF (R(K,L).GT.9.246570341781329D0) SSPs=0.0D0

C Sppsig
      SPPs=G26*DQDR + G27*4.0D0*Q(K,L)*DQDR
     /   +G28*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + G29*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +G30*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +G31*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +G32*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +G33*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +G34*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +G35*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +G36*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

       IF (R(K,L).GT.9.273572893686612D0) SPPs=0.0D0

C   Spppi
      SPPp=G38*DQDR + G39*4.0D0*Q(K,L)*DQDR
     /   +G40*(-3.0D0*DQDR + 12.0D0*(Q(K,L)**2)*DQDR)
     /   + G41*(-16.0D0*Q(K,L)*DQDR + 32.0D0*(Q(K,L)**3)*DQDR) +G42*(5.0D0*DQDR -60.0D0*(Q(K,L)**2)*DQDR
     /   + 80.0D0*(Q(K,L)**4)*DQDR) +G43*(36.0D0*Q(K,L)*DQDR - 192.0D0*(Q(K,L)**3)*DQDR +192.0D0*(Q(K,L)**5)*DQDR)
     /   +G44*(-7.0D0*DQDR +168.0D0*(Q(K,L)**2)*DQDR - 560.0D0*(Q(K,L)**4)*DQDR +448.0D0*(Q(K,L)**6)*DQDR)
     /   +G45*(-64.0D0*Q(K,L)*DQDR +640.0D0*(Q(K,L)**3)*DQDR -1536.0D0*(Q(K,L)**5)*DQDR +1024.0D0*(Q(K,L)**7)*DQDR)
     /   +G46*(9.0D0*DQDR -360.0D0*(Q(K,L)**2)*DQDR +2160*(Q(K,L)**4)*DQDR
     /   -4032.0D0*(Q(K,L)**6)*DQDR +2304.0D0*(Q(K,L)**8)*DQDR)
     /   +G47*(100.0D0*Q(K,L)*DQDR - 1600*(Q(K,L)**3)*DQDR + 6720*(Q(K,L)**5)*DQDR - 10240*(Q(K,L)**7)*DQDR
     /       +5120*(Q(K,L)**9)*DQDR)
     /   +G48*(-11.0D0*DQDR + 660.0D0*(Q(K,L)**2)*DQDR - 6160*(Q(K,L)**4)*DQDR +19712*(Q(K,L)**6)*DQDR 
     /   -25344*(Q(K,L)**8)*DQDR + 11264*(Q(K,L)**10)*DQDR)

       IF (R(K,L).GT.9.10748558948335D0) SPPp=0.0D0

      RETURN
      END
