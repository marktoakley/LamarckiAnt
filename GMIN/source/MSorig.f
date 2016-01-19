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

      SUBROUTINE MSORIG(N,XS,deriv1st,ENERGY,GTEST)
      USE commons
      USE consts
      USE distance
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER I, J, K, N, AI, AJ, NDIM, NMAX, I2, J2, K2, AI2, INFO, 
     1        MX, MY, MZ, J1
      DOUBLE PRECISION DIST,ZETA,NB,VECT(4),EPS(4),WORK(24*NATOMS),
     1                 S(4*NATOMS,4*NATOMS),ENERGY,XS(3*NATOMS),
     2                 H(4*NATOMS,4*NATOMS),VECTOR(3), deriv1st(3*NATOMS)
      DOUBLE PRECISION UREP,UBOND,ELEMENT
      LOGICAL FTEST
      COMMON /FAIL/ FTEST
      COMMON /NEW/ ELEMENT
      allocate(DIAG(4*natoms), R(natoms,natoms), S2(natoms,natoms), LAMBDA(natoms,natoms),
     &         REPULSE(natoms,natoms),BONDAGE(natoms,natoms), KAPPANEW(natoms,natoms),
     &         Horig(4*natoms,4*natoms), SVtype(8*natoms,natoms), DIRCOS(natoms,natoms,3))
 

C
C************ START PERIODIC BOUNDARY CONDITIONS **************************
C  Deal with atoms leaving the box:
C
      IF (PERIODIC) THEN
C        PRINT*,'Points before:'
C        WRITE(*,'(3F15.7)') (XS(J1),J1=1,3*N)
         DO J1=1,N
            J2=3*(J1-1)
            XS(J2+1)=XS(J2+1) - BOXLX*DNINT(XS(J2+1)/BOXLX) 
            XS(J2+2)=XS(J2+2) - BOXLY*DNINT(XS(J2+2)/BOXLY)
            XS(J2+3)=XS(J2+3) - BOXLZ*DNINT(XS(J2+3)/BOXLZ)
         ENDDO
C        PRINT*,'Points after:'
C        WRITE(*,'(3F15.7)') (XS(J1),J1=1,3*N)
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is used:
C
         DO I=1,N
            I2=3*(I-1)
            R(I,I)=0.0D0
            DO J=I+1,N
               J2=3*(J-1)
               VECTOR(1)=XS(J2+1)-XS(I2+1)
               VECTOR(2)=XS(J2+2)-XS(I2+2)
               VECTOR(3)=XS(J2+3)-XS(I2+3)
               MX=NINT(VECTOR(1)/BOXLX)
               MY=NINT(VECTOR(2)/BOXLY)
               MZ=NINT(VECTOR(3)/BOXLZ)
               VECTOR(1)= VECTOR(1) - BOXLX * MX
               VECTOR(2)= VECTOR(2) - BOXLY * MY
               VECTOR(3)= VECTOR(3) - BOXLZ * MZ
               DIST=VECTOR(1)**2+VECTOR(2)**2+VECTOR(3)**2
               R(I,J)=SQRT(DIST)
               R(J,I)=R(I,J)
               DO K=1,3
                  DIRCOS(I,J,K)=VECTOR(K)/R(I,J)
                  DIRCOS(J,I,K)=-DIRCOS(I,J,K)
               ENDDO
C              WRITE(*,'(2I4,3F15.7)') (I,J,(DIRCOS(I,J,K),K=1,3))
            ENDDO
            DO K=1,3
               DIRCOS(I,I,K)=0.0D0
            ENDDO
         ENDDO
C************ END PERIODIC BOUNDARY CONDITIONS ****************************
      ELSE

C   Calculate distance matrix

         DO I=1,N
          I2=3*(I-1)
          R(I,I)=0.0D0
             DO J=I+1,N
               J2=3*(J-1)
               DIST=(XS(I2+1)-XS(J2+1))**2+(XS(I2+2)-XS(J2+2))**2
     1                   +(XS(I2+3)-XS(J2+3))**2
                R(I,J)=SQRT(DIST)
                R(J,I)=R(I,J)

C   Now the direction cosines
C   Which incidentally are the projections of a unit vector onto 
C   the axes and NOT the cosines.
C   Note we have not calculated DIRCOS(I,I)

               DO K=1,3
                  DIRCOS(I,J,K)=(XS(J2+K)-XS(I2+K))/R(I,J)
                  DIRCOS(J,I,K)=-DIRCOS(I,J,K)
               ENDDO
            ENDDO
            DO K=1,3
               DIRCOS(I,I,K)=0.0D0
            ENDDO
         ENDDO 
      ENDIF

C   Calculate Repulsive Free Energy
C   CHI(R) is the repulsive pair potential, the sum of which is 
C   equal to the repulsive energy (UREP) in eV.
C   BETA is 4/r_0, where r_0 is one half the dimer bond length.
C   D is the bond length for the crystal in Angstroms.

      UREP=0.0D0
      DO I=1,N
         DO J=I+1,N
             REPULSE(J,I)=CHI*DEXP(-BETA*(R(J,I)-D))
             UREP=UREP+REPULSE(J,I)
         END DO
      END DO 
C     WRITE(6,10) 'UREP=',UREP
C10    FORMAT(A8, F20.17)

C   Calculate ad hoc bond counting term UBOND in eV.
C   This is using the a and b parameters determined using the PRB
C   50 5645 1994. It is just a correction to the cohesive energy
C   for clusters of different sizes in order to obtain a good agreement
C   with ab initio values.
C   Sum is over all bonds.

      NB=0.0D0

      DO I=1,N
         DO J=I+1,N
             BONDAGE(J,I)=DEXP((R(J,I)-RC)/DELTA)
             NB=NB+1.0D0/(BONDAGE(J,I)+1.0D0)
         END DO
      END DO
      UBOND=-DBLE(N)*(A*(NB/DBLE(N))+B)
C     WRITE(6,10) 'UBOND=',UBOND

C   Calculate electronic energy UEL
C   Not sure how the hell to do this but I'll see what I come up with.
C   First find the H(I,J) matrix.
C   The interaction of the 3s and 3p atomic orbitals on one atom
C   with those on another atom leads to 8 different interactions
C   Obvious right? 'cos there are 8AOs involved.
C   Still got to include the fact that only the nearest neighbour
C   Interactions should be included...but I am working on it.



C   Suppose it would help if I found the universal parameters 
C   I am assuming that the ones quoted in PRB 50 12754 1993
C   are at the crystal bond length (D0=2.35)
C   VECT(1) is Vsssig, 2 Vspsig, 3 Vppsig and 4 Vpppi evaluated at D0
C   Thus, for every atom there are 4 atomic orbitals and the matrices
C   will be (4N, 4N) in size.

      VECT(1)=Vsssig
      VECT(2)=Vspsig
      VECT(3)=Vppsig
      VECT(4)=Vpppi

      EPS(1)=EPSs+EPSs
      EPS(2)=EPSs+EPSp
      EPS(3)=EPSp+EPSp
      EPS(4)=EPSp+EPSp
    
      ZETA=2.0D0*SQRT(3.0D0)
  
C   Only include overlap or interaction with different atoms here.
               
      DO I=1,N
         AI=4*(I-1)
         AI2=2*AI
         DO J=I+1,N
            AJ=4*(J-1)

C   Dependence of V on distance between atoms R(I,J)
C   Similarly, S must scale with R(I,J)

            DO K=1,4
               SVtype(AI2+4+K,J)=VECT(K)*DEXP(-ALPHA*(R(I,J)-D))
               SVtype(AI2+K,J)=2.0D0*SVtype(AI2+4+K,J)/(KAPPA*EPS(K))
            END DO

C   Set DIRCOS2 to specific value of DIRCOS(I,J,K) so that
C   it is easier to read the equations

C           DO K=1,3
C              DIRCOS2(K)=DIRCOS(I,J,K)
C           END DO
C   Calculate the overlap matrix S(I,J), related to these parameters
C   by the Slater-Koster scheme

            S((AI+1),(AJ+1))=SVtype(AI2+1,J)
            S((AI+1),(AJ+2))=DIRCOS(I,J,1)*SVtype(AI2+2,J)
            S((AI+1),(AJ+3))=DIRCOS(I,J,2)*SVtype(AI2+2,J)
            S((AI+1),(AJ+4))=DIRCOS(I,J,3)*SVtype(AI2+2,J)

            S((AI+2),(AJ+1))=-S((AI+1),(AJ+2))
            S((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*SVtype(AI2+3,J)+(1.0D0-
     1                       DIRCOS(I,J,1)**2)*SVtype(AI2+4,J)
            S((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))
            S((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))

            S((AI+3),(AJ+1))=-S((AI+1),(AJ+3))
            S((AI+3),(AJ+2))=S((AI+2),(AJ+3))
            S((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*SVtype(AI2+3,J)+(1.0D0-
     1                       DIRCOS(I,J,2)**2)*SVtype(AI2+4,J)
            S((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))
            
            S((AI+4),(AJ+1))=-S((AI+1),(AJ+4))
            S((AI+4),(AJ+2))=S((AI+2),(AJ+4))
            S((AI+4),(AJ+3))=S((AI+3),(AJ+4))
            S((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*SVtype(AI2+3,J)+(1.0D0-
     1                       DIRCOS(I,J,3)**2)*SVtype(AI2+4,J)
            
            S2(I,J)=(SVtype(AI2+1,J)-ZETA*SVtype(AI2+2,J)-3.0D0*
     1                        SVtype(AI2+3,J))/4.0D0

C   Hamiltonian is found using the Slater-Koster scheme
C   as shown in their paper and Harrison's book p.481
C   which relates H to the parameter V using direction cosines.
C   LAMBDA changes for each I and J
      
      LAMBDA(I,J)=1.0D0+1.0D0/KAPPA-S2(I,J)**2

            H((AI+1),(AJ+1))=SVtype(AI2+5,J)*LAMBDA(I,J)
            H((AI+1),(AJ+2))=DIRCOS(I,J,1)*SVtype(AI2+6,J)*LAMBDA(I,J)
            H((AI+1),(AJ+3))=DIRCOS(I,J,2)*SVtype(AI2+6,J)*LAMBDA(I,J)
            H((AI+1),(AJ+4))=DIRCOS(I,J,3)*SVtype(AI2+6,J)*LAMBDA(I,J)

            H((AI+2),(AJ+1))=-H((AI+1),(AJ+2))
            H((AI+2),(AJ+2))=(DIRCOS(I,J,1)**2*SVtype(AI2+7,J)+(1.0D0-
     1           DIRCOS(I,J,1)**2)*SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            
            H((AI+3),(AJ+1))=-H((AI+1),(AJ+3))
            H((AI+3),(AJ+2))=H((AI+2),(AJ+3))
            H((AI+3),(AJ+3))=(DIRCOS(I,J,2)**2*SVtype(AI2+7,J)+(1.0D0-
     1           DIRCOS(I,J,2)**2)*SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            
            H((AI+4),(AJ+1))=-H((AI+1),(AJ+4))
            H((AI+4),(AJ+2))=H((AI+2),(AJ+4))
            H((AI+4),(AJ+3))=H((AI+3),(AJ+4))
            H((AI+4),(AJ+4))=(DIRCOS(I,J,3)**2*SVtype(AI2+7,J)+(1.0D0-
     1           DIRCOS(I,J,3)**2)*SVtype(AI2+8,J))*LAMBDA(I,J)


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

C     ELEMENT=S(1,6)
      
      DO I=1,4*N
          DO J=I+1,4*N
             S(J,I)=S(I,J)
             H(J,I)=H(I,J)
          END DO
    

C   Interaction/overlap of the same AO on the same atom
         S(I,I)=1.0D0

      END DO
      DO I=1,4*N
         DO J=1,4*N
            Horig(I,J)=H(I,J)
         END DO
      END DO
      
      NDIM=4*N
      NMAX=4*NATOMS

      FTEST=.FALSE.
 
C     CALL LANCZOS(NDIM,H,DIAG,1.0D0,-1.0D0,NFOUND)
C     CALL CGMIN(H,S,DIAG)
C     INFO=0

      CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, 24*NATOMS, INFO )

C
C  This sort isn't needed for the SGI lapack.f routine which
C  incorporates it. 
C
      CALL EIGSRT(DIAG,H,NDIM,NMAX)
C     CALL EIGSRT(DIAG,H,NFOUND,NMAX)
C     WRITE(*,'(6F15.7)') (DIAG(J1),J1=1,NDIM)

      IF (INFO.NE.0) THEN
         FTEST=.TRUE.
         ENERGY=1.0D6
         PRINT*,'DSYGV failed - INFO=',INFO
         RETURN
      ENDIF

      ENERGY=0.0D0
      DO J=4*N,(2*N+1),-1
         ENERGY=ENERGY+DIAG(J)
      END DO

C   For single atom there should be no contribution from the 
C   bond counting term      

      IF (N .EQ. 1) THEN
         ENERGY=2*ENERGY
      ELSE
         ENERGY=UREP+UBOND+2*ENERGY
      END IF

C   Now call the all new originally crafted subroutine to calculate
C   the first derivatives.

      IF (GTEST) THEN
         CALL DERIVS1(N,deriv1st,H,natoms)
C        CALL DERIVS1(N)
      ELSE
         RETURN
      ENDIF

      RETURN
      END

C   This is phase two where we calculate the first derivatives
C   of the silicon potential coded in TB.Si.f. It is in tight
C   binding parameterised form.
C   Off I go...
C   This is derivs1.f.jloop.betterish in the resrves directory
C   for the record

      SUBROUTINE DERIVS1(N,deriv1st,EIGVEC,natoms)
      USE consts
      USE distance
      IMPLICIT NONE 
      INTEGER I, J, K, N, AI, AJ, LEV, SV, ASV, NUM1, NUM2, I2, J2,natoms
      DOUBLE PRECISION REP1st(NATOMS,3),BOND1st(NATOMS,3), RR, D1, D2, D3,
     1                 ELEC1st(NATOMS,3), diffH, E1, E2, AS1, DUMMY,
     2                 diffSV(2,4,4), SDC1, SDC2, SDC3, SQC1, SQC2, SQC3,
     3                 diffS2,R1st(NATOMS,NATOMS),EIGVEC(4*NATOMS,4*NATOMS),
     4                 B1st(NATOMS,NATOMS), deriv1st(3*NATOMS)
      DOUBLE PRECISION LMC(NATOMS,NATOMS),ALPHASUM1(NATOMS,NATOMS),ALPHASUM2,SVSUB,
     1                 SQDIRCOS(3*NATOMS,NATOMS),SUBDIRCOS(3*NATOMS,NATOMS),LONG,OPTION(3)

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

      DO K=1,3
         DO I=1,N
            REP1st(I,K)=0.0D0
                   DO J=I+1,N
               R1st(J,I)=DIRCOS(I,J,K)*REPULSE(J,I)
               R1st(I,J)=-R1st(J,I)
            END DO
            R1st(I,I)=0.0D0
            D1=REP1st(I,K)
            DO J=1,N
               D1=D1+BETA*R1st(J,I)
            END DO
            REP1st(I,K)=D1
         END DO
      END DO

C   Easy isn't it.
C   Now do the bond counting term derivative.

      DO K=1,3
         DO I=1,N
            BOND1st(I,K)=0.0D0
            DO J=I+1,N
               B1st(J,I)=DIRCOS(I,J,K)*BONDAGE(J,I)/(BONDAGE(J,I)+1.0D0)**2
               B1st(I,J)=-B1st(J,I)
            END DO
            B1st(I,I)=0.0D0
            D1=BOND1st(I,K)
            DO J=1,N
               D1=D1-A*B1st(J,I)/DELTA
            END DO
            BOND1st(I,K)=D1
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

                 LMC(J,I)=DIRCOS(J,I,1)*DIRCOS(J,I,2)*DIRCOS(J,I,3)

                 ALPHASUM1(J,I)=(ALPHA+2.0D0/R(J,I))
        END DO
      END DO

C     Entering the differentiation of electronic energy zone
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


C   EXT is for when CO+X is greater than 3 for CO=3 and 2 ie. Z and Y
C   EXT2 is for when this only the case for CO=3
C   It is just a method to permute the X, Y and Z values etc.

           D1=DIRCOS(J,I,1)
           D2=DIRCOS(J,I,2)
           D3=DIRCOS(J,I,3)

           AS1=ALPHASUM1(J,I)
           SDC1=SUBDIRCOS(J2+1,I)
           SDC2=SUBDIRCOS(J2+2,I)
           SDC3=SUBDIRCOS(J2+3,I)
           SQC1=SQDIRCOS(J2+1,I)
           SQC2=SQDIRCOS(J2+2,I)
           SQC3=SQDIRCOS(J2+3,I)

                 LONG=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

                 OPTION(1)=-D2*LONG
                 OPTION(2)=-D3*LONG
                 OPTION(3)=-LMC(J,I)*AS1

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

                 DO SV=1,2
                    ASV=2*AI+4*(SV-1)
                    ALPHASUM2=D1*SVtype(ASV+2,J)*(ALPHA+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-ALPHA*D1*SVtype(ASV+1,J)
C                 diffSV(SV,1,2)=-SVtype(ASV+2,J)*(RR*SDC1+ALPHA*SQC1)
C             diffSV(SV,1,3)=-D2*ALPHASUM2
C                 diffSV(SV,1,4)=-D3*ALPHASUM2
                  diffSV(SV,1,2)= SVtype(ASV+2,J)*(RR*SDC1+ALPHA*SQC1)
              diffSV(SV,1,3)= D2*ALPHASUM2
                  diffSV(SV,1,4)= D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*(ALPHA*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1)
                  diffSV(SV,2,3)=SVSUB*OPTION(1)
                  diffSV(SV,2,4)=SVSUB*OPTION(2)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*SQC2*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC2-2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=SVSUB*OPTION(3)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*SQC3*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC3-2.0D0*SQC3*RR))
               ENDDO

           diffS2=ALPHA*DIRCOS(I,J,1)*S2(I,J)
           E1=0.0D0
           DO NUM2=1,4   
             DO NUM1=1,4

C             diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-2.0D0*S2(I,J)*Horig((AJ+NUM1),(AI+NUM2))/LAMBDA(I,J)*diffS2
              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-2.0D0*S2(I,J)*Horig((AJ+NUM2),(AI+NUM1))/LAMBDA(I,J)*diffS2
              E2=diffSV(1,NUM1,NUM2)
C             WRITE(*,'(A,4I3,2F15.7)') 'I,J,NUM1,NUM2,diffH,diffS=',I,J,NUM1,NUM2,diffH,E2


              DO LEV=4*N,2*N+1,-1
C              E1=E1+EIGVEC(LEV,AI+NUM2)*EIGVEC(LEV,AJ+NUM1)*(diffH-DIAG(LEV)*E2)
               E1=E1+EIGVEC(LEV,AJ+NUM2)*EIGVEC(LEV,AI+NUM1)*(diffH-DIAG(LEV)*E2)
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

           AS1=ALPHASUM1(J,I)
           SDC1=SUBDIRCOS(J2+2,I)
           SDC2=SUBDIRCOS(J2+3,I)
           SDC3=SUBDIRCOS(J2+1,I)
           SQC1=SQDIRCOS(J2+2,I)
           SQC2=SQDIRCOS(J2+3,I)
           SQC3=SQDIRCOS(J2+1,I)

               LONG=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

               OPTION(1)=-D2*LONG
               OPTION(2)=-D3*LONG
               OPTION(3)=-LMC(J,I)*AS1

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

               DO SV=1,2
                  ASV=2*AI+4*(SV-1)
                  ALPHASUM2=D1*SVtype(ASV+2,J)*(ALPHA+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-ALPHA*D1*SVtype(ASV+1,J)
                  diffSV(SV,1,3)= SVtype(ASV+2,J)*(RR*SDC1+ALPHA*SQC1)
              diffSV(SV,1,4)= D2*ALPHASUM2
                  diffSV(SV,1,2)= D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*(ALPHA*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1)
                  diffSV(SV,2,3)=SVSUB*OPTION(2)
                  diffSV(SV,2,4)=SVSUB*OPTION(3)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*SQC2*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC2-2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=SVSUB*OPTION(1)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*SQC3*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC3-2.0D0*SQC3*RR))

               END DO ! SV

           diffS2=ALPHA*DIRCOS(I,J,2)*S2(I,J)
           E1=0.0D0
           DO NUM2=1,4   
             DO NUM1=1,4

              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-2.0D0*S2(I,J)*Horig((AJ+NUM2),(AI+NUM1))/LAMBDA(I,J)*diffS2

              E2=diffSV(1,NUM1,NUM2)
              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AJ+NUM2)*EIGVEC(LEV,AI+NUM1)*(diffH-DIAG(LEV)*E2)
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

           AS1=ALPHASUM1(J,I)
           SDC1=SUBDIRCOS(J2+3,I)
           SDC2=SUBDIRCOS(J2+1,I)
           SDC3=SUBDIRCOS(J2+2,I)
           SQC1=SQDIRCOS(J2+3,I)
           SQC2=SQDIRCOS(J2+1,I)
           SQC3=SQDIRCOS(J2+2,I)

               LONG=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

               OPTION(1)=-D2*LONG
               OPTION(2)=-D3*LONG
               OPTION(3)=-LMC(J,I)*AS1

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

               DO SV=1,2
                  ASV=2*AI+4*(SV-1)
                  ALPHASUM2=D1*SVtype(ASV+2,J)*(ALPHA+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-ALPHA*D1*SVtype(ASV+1,J)
                  diffSV(SV,1,4)= SVtype(ASV+2,J)*(RR*SDC1+ALPHA*SQC1)
              diffSV(SV,1,2)= D2*ALPHASUM2
                  diffSV(SV,1,3)= D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*(ALPHA*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1)
                  diffSV(SV,2,3)=SVSUB*OPTION(3)
                  diffSV(SV,2,4)=SVSUB*OPTION(1)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*SQC2*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC2-2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=SVSUB*OPTION(2)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*SQC3*AS1+SVtype(ASV+4,J)*(-ALPHA*SDC3-2.0D0*SQC3*RR))

               END DO ! SV

           diffS2=ALPHA*DIRCOS(I,J,3)*S2(I,J)
           E1=0.0D0
           DO NUM2=1,4   
             DO NUM1=1,4

              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-2.0D0*S2(I,J)*Horig((AJ+NUM2),(AI+NUM1))/LAMBDA(I,J)*diffS2

              E2=diffSV(1,NUM1,NUM2)
              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AJ+NUM2)*EIGVEC(LEV,AI+NUM1)*(diffH-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,3)=ELEC1st(I,3)+2.0D0*E1
           ELEC1st(J,3)=ELEC1st(J,3)-2.0D0*E1

         END DO  !J (1st time)
        END DO  !I(1st time)

C   In order to differentiate by the coordinate of atom K, you only get 
C   a value for the derivative if one of the two atoms involved in the
C   pair interaction is K otherwise it is zero. This is achieved by use
C   of the FACTOR term. I'm sure this is highly inefficient BUT if it
C   works then that is worthwhile. I can improve it later....when I'm 
C   drawing my pension.            
C   Have to keep it in a loop where I and J are the atom number rather
C   than the AO.

C   Recent changes involve replacing derivV with its specific value 
C   FACTOR*diffSV(2,X,Y) and perhaps similarly for derivS.
C   diffSV and diffH will be zero.
C   Must be aware that although none of the diagonal elements are defined
C   they are equal to zero. Not deemed necessary to calc. them since
C   the loop never reaches I=J.  Good to use J=I+1 since uses symmetry to 
C   cut down on calculations.

C   Keep diffS2 since it is only dependent on which atoms are interacting 
C   not the AOs

C   Include factor of 2 for ELEC1st since two electrons fill each MO

      DO K=1,3
         DO I=1,N
            I2=3*(I-1)
            deriv1st(I2+K)=2.0D0*ELEC1st(I,K)+BOND1st(I,K)+REP1st(I,K)
         END DO
      END DO

      END
