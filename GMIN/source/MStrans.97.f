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
C   This program is supposed to calculate the energy 
C   and first derivatives of Silicon 
C   clusters using a Tight Binding (TB) approximation.  
C   This will then be used in conjunction with GMIN
C   to solve the problems of the universe
C   and in particular those of the silicon cluster world.

C   This is a subsequent adaptation to the original potential
C   and uses Menon and Subbaswamy's newest transferable potential 
C   (PRB 55 9231 1997) 

C   Things to note about THIS version of the transferable potential:
C   Gets the same energy as M+S, compatible with ab initio
C   No longer uses coordination dependent term as in previous versions of their potential
C   - apparently this is key to transferability.

C   This is the UPDATED TRANSFERABLE version from M+S's April 1997 PRB paper
C   it makes S2 as in the original potential AND K(r) an exponential
C   dependence.
C   Still no bond counting term.

C   NOT really MSORIG - purely for the purpose of easy testing
C   it is IN fact MSTRANS97

      SUBROUTINE MSTRANS97(N,XS,deriv1st,ENERGY,GTEST)
      USE commons
      USE consts_trans_97
      USE distance
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER I, J, K, N, AI, AJ, NDIM, NMAX, I2, J2, K2, AI2,  INFO,
     1        MX, MY, MZ, J1, NFOUND
      DOUBLE PRECISION DIST,ZETA,WORK(24*NATOMS), deriv1st(3*NATOMS),
     1                 S(4*NATOMS,4*NATOMS),ENERGY,XS(3*NATOMS),COHESIVE,
     2                 H(4*NATOMS,4*NATOMS),D1,D2,D3,VECTOR(3)
C    5                 VECT(4),EPS(4),NB
      DOUBLE PRECISION UREP,UBOND
      LOGICAL FTEST
      COMMON /FAIL/ FTEST
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
     1                +(XS(I2+3)-XS(J2+3))**2
           R(I,J)=SQRT(DIST)
             R(J,I)=R(I,J)

C   Now the direction cosines
C   Which incidentally are the projections of a unit vector onto 
C   the axes and NOT the cosines.
C   Note we have not calculated DIRCOS(I,I)

            DO K=1,3
               DIRCOS(I,J,K)=(XS(J2+K)-XS(I2+K))/R(I,J)
               DIRCOS(J,I,K)=-DIRCOS(I,J,K)
            END DO

         END DO
      END DO 
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
C      WRITE(6,10) 'UREP=',UREP
C10    FORMAT(A8, F20.17)

C   Calculate ad hoc bond counting term UBOND in eV.
C   This is using the a and b parameters determined using the PRB
C   50 5645 1994. It is just a correction to the cohesive energy
C   for clusters of different sizes in order to obtain a good agreement
C   with ab initio values.
C   Sum is over all bonds.

C   To attain something like the values in PRB 50 11577 1994 then ENERGY
C   is just made up of UREP and electronic ENERGY. ie. does NOT contain
C   UBOND.
C   Ramakrishna and Bahel, JCP 104 9833 1996 use M+S's transferable potential
C   but also include bond counting term in their investigation of Si_12. Not sure 
C   how this affects things...may try both ways. 

C   NO BOND COUNTING TERM in this version

C      NB=0.0D0
C
C      DO I=1,N
C         DO J=I+1,N
C             BONDAGE(J,I)=DEXP((R(J,I)-RC)/DELTA)
C           NB=NB+1.0D0/(BONDAGE(J,I)+1.0D0)
C         END DO
C      END DO
C      UBOND=-DBLE(N)*(A*(NB/DBLE(N))+B)
C      WRITE(6,10) 'UBOND=',UBOND

C   Calculate electronic energy (ENERGY)
C   Not sure how the hell to do this but I'll see what I come up with.
C   First find the H(I,J) matrix.
C   The interaction of the 3s and 3p atomic orbitals on one atom
C   with those on another atom leads to 8 different interactions
C   Obvious right? 'cos there are 8AOs involved.

C   UNRESOLVED FACT
C   Still got to include the fact that only the nearest neighbour
C   Interactions should be included...but I am working on it.
C   Could be the reason why M+S potential is not so good.
C   Some potentials use a decay factor (Khan and Broughton,
C   PRB 43 11754 1991)

C   UNRESOLVED FACT
C   Suppose it would help if I found the universal parameters 
C   I am assuming that the ones quoted in PRB 50 12754 1993
C   are at the crystal bond length (D0=2.35)
C   VECT(1) is Vsssig, 2 Vspsig, 3 Vppsig and 4 Vpppi evaluated at D0.
C   STILL do NOT know how M+S found these values
C   Thus, for every atom there are 4 atomic orbitals and the matrices
C   will be (4N, 4N) in size.
C   These are now defined by a DATA line in consts.f90
    
      ZETA=2.0D0*SQRT(3.0D0)
  
C   Only include overlap or interaction with different atoms here.
               
C      PRINT*,'coords:'
C      WRITE(*,'(3F20.10)') (XS(I),I=1,3*N)
      DO I=1,N
         AI=4*(I-1)
         AI2=2*AI
         DO J=I+1,N
            AJ=4*(J-1)

C   Transferable TB potential introduces a distance potential to the
C   non-orthogonality coefficient K (1/KAPPANEW)
C   Problems arise when R is large and therfore DEXP(R) is HUGE
C   The program fails due to an overflow
C   Therefore, since K is only ever in calculations as the reciprocal
C   I have calculated that instead.
            WRITE(6,'(A,2I2,2F20.10)') 'I,J,R(J,I),arg=',I,J,R(J,I),SIGMA*(R(I,J)-D)**2
C            KAPPANEW(J,I)=KAPPA*DEXP(SIGMA*(R(I,J)-D)**2)
            KAPPANEW(J,I)=DEXP(-SIGMA*(R(I,J)-D)**2)/KAPPA

C   Dependence of V on distance between atoms R(I,J)
C   Similarly, S must scale with R(I,J)

            DO K=1,4
               SVtype(AI2+4+K,J)=VECT(K)*DEXP(-ALPHA*(R(I,J)-D))
C               SVtype(AI2+K,J)=2.0D0*SVtype(AI2+4+K,J)/(KAPPANEW(J,I)*EPS(K))
               SVtype(AI2+K,J)=2.0D0*SVtype(AI2+4+K,J)*KAPPANEW(J,I)/(EPS(K))
            END DO

C   Set variable with no arguments equal to DIRCOS(I,J,K) 

            D1=DIRCOS(I,J,1)
            D2=DIRCOS(I,J,2)
            D3=DIRCOS(I,J,3)

C   Calculate the overlap matrix S(I,J), related to these parameters
C   by the Slater-Koster scheme

            S((AI+1),(AJ+1))=SVtype(AI2+1,J)
            S((AI+1),(AJ+2))=D1*SVtype(AI2+2,J)
            S((AI+1),(AJ+3))=D2*SVtype(AI2+2,J)
            S((AI+1),(AJ+4))=D3*SVtype(AI2+2,J)

            S((AI+2),(AJ+1))=-S((AI+1),(AJ+2))
            S((AI+2),(AJ+2))=D1**2*SVtype(AI2+3,J)+(1.0D0-
     1                       D1**2)*SVtype(AI2+4,J)
            S((AI+2),(AJ+3))=D1*D2*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))
            S((AI+2),(AJ+4))=D1*D3*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))

            S((AI+3),(AJ+1))=-S((AI+1),(AJ+3))
            S((AI+3),(AJ+2))=S((AI+2),(AJ+3))
            S((AI+3),(AJ+3))=D2**2*SVtype(AI2+3,J)+(1.0D0-
     1                       D2**2)*SVtype(AI2+4,J)
            S((AI+3),(AJ+4))=D2*D3*(SVtype(AI2+3,J)
     1                       -SVtype(AI2+4,J))
            
            S((AI+4),(AJ+1))=-S((AI+1),(AJ+4))
            S((AI+4),(AJ+2))=S((AI+2),(AJ+4))
            S((AI+4),(AJ+3))=S((AI+3),(AJ+4))
            S((AI+4),(AJ+4))=D3**2*SVtype(AI2+3,J)+(1.0D0-
     1                       D3**2)*SVtype(AI2+4,J)
           
C   Transferable potential introduces the dependence of S2 on SVtype(4) reducing the
C   value of CHI.
 
            S2(I,J)=(SVtype(AI2+1,J)-ZETA*SVtype(AI2+2,J)-3.0D0*
     1                        SVtype(AI2+3,J))/4.0D0

C   Hamiltonian is found using the Slater-Koster scheme
C   as shown in their paper and Harrison's book p.481
C   which relates H to the parameter V using direction cosines.
C   LAMBDA changes for each I and J. Notice additional
C   dependence on position of atoms through KAPPANEW
      
C      LAMBDA(I,J)=1.0D0+1.0D0/KAPPANEW(J,I)-S2(I,J)**2
      LAMBDA(I,J)=1.0D0+KAPPANEW(J,I)-S2(I,J)**2

            H((AI+1),(AJ+1))=SVtype(AI2+5,J)*LAMBDA(I,J)
            H((AI+1),(AJ+2))=D1*SVtype(AI2+6,J)*LAMBDA(I,J)
            H((AI+1),(AJ+3))=D2*SVtype(AI2+6,J)*LAMBDA(I,J)
            H((AI+1),(AJ+4))=D3*SVtype(AI2+6,J)*LAMBDA(I,J)

            H((AI+2),(AJ+1))=-H((AI+1),(AJ+2))
            H((AI+2),(AJ+2))=(D1**2*SVtype(AI2+7,J)+(1.0D0-
     1           D1**2)*SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+2),(AJ+3))=D1*D2*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+2),(AJ+4))=D1*D3*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            
            H((AI+3),(AJ+1))=-H((AI+1),(AJ+3))
            H((AI+3),(AJ+2))=H((AI+2),(AJ+3))
            H((AI+3),(AJ+3))=(D2**2*SVtype(AI2+7,J)+(1.0D0-
     1           D2**2)*SVtype(AI2+8,J))*LAMBDA(I,J)
            H((AI+3),(AJ+4))=D2*D3*(SVtype(AI2+7,J)
     1           -SVtype(AI2+8,J))*LAMBDA(I,J)
            
            H((AI+4),(AJ+1))=-H((AI+1),(AJ+4))
            H((AI+4),(AJ+2))=H((AI+2),(AJ+4))
            H((AI+4),(AJ+3))=H((AI+3),(AJ+4))
            H((AI+4),(AJ+4))=(D3**2*SVtype(AI2+7,J)+(1.0D0-
     1           D3**2)*SVtype(AI2+8,J))*LAMBDA(I,J)


        END DO !Of J loop

C   Interaction/Overlap of different AOs on the same atom
C   Since require whole of S and H matrices for diagonalisation of matrix routines
         
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

      END DO !End of I

C   Now we have to make sure the matrix S and H are symmetric.
C   PERHAPS NOT: the lapack routine takes 'U' ie. the upper matrix
C   Similarly Horig is only used in J>I routines.      
C      WRITE(6,*) 'I,J,S(I,J)'
      DO I=1,4*N
        DO J=I+1,4*N
            Horig(J,I)=H(I,J)
C            Horig(I,J)=H(I,J)
          END DO

C   Interaction/overlap of the same AO on the same atom
         S(I,I)=1.0D0

      END DO
      
      NDIM=4*N
      NMAX=4*NATOMS

      FTEST=.FALSE.

C   This calls the lapack routine to diagonalize the matrix:
C   1 - specifies type of problem to be solved
C   V - calculate eigenvalues and eigenvectors
C   U - upper triangles of next two matrices are stored
C   NDIM - order of next two matrices
C   H - matrix of eigenvectors. Leading NxN upper triangle part of H contains 
C       upper triangle part of H
C   NMAX - leading dimension of array H/S
C   S - cholesky decompostion
C   DIAG - Eigenvalues in ascending order
C   WORK - returns optimal LWORK=24*NATOMS
C   LWORK =24*NATOMS - length of array work
C   INFO - =0 successful exit, otherwise error of some kind.

      CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, 24*NATOMS, INFO )

C   This routine is not required for SGI since lapack routine includes it.

      CALL EIGSRT(DIAG,H,NDIM,NMAX)


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
         ENERGY=2.0D0*ENERGY
      ELSE IF(PERIODIC) THEN
         ENERGY=UREP+2.0D0*ENERGY
      ELSE
C         ENERGY=UREP+UBOND+2.0D0*ENERGY
C   To get the same results as in PRB 50 11577 1994 then do not include bond counting
C   term but do include N*1.
         ENERGY=UREP+2.0D0*ENERGY+FUDGE*DBLE(N)
         COHESIVE=(ENERGY-N*2.0D0*(EPSs+EPSp))/N
      END IF
C      WRITE(6,*) 'ENERGY in potential=',ENERGY
C      WRITE(6,*) 'COHESIVE ENERGY=',COHESIVE

C   Now call the all new originally crafted subroutine to calculate
C   the first derivatives.

      IF (GTEST) THEN
         CALL DERIVS1T97(N,deriv1st,H,natoms)
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

C   This is the IMPROVED transferable potential April 1997 PRB.
C   Changes to S2 and K(r) and the appropriate derivatives.

      SUBROUTINE DERIVS1T97(N,deriv1st,EIGVEC,natoms)
      USE consts_trans_97
      USE distance
      IMPLICIT NONE 
      INTEGER I, J, K , N, AI, AJ, LEV, SV, ASV, NUM1, NUM2, 
     1        I2,J2,natoms
      DOUBLE PRECISION REP1st(NATOMS,3), RR, D1, D2, D3, deriv1st(3*NATOMS),
     1                 ELEC1st(NATOMS,3), diffH, E1, E2, AS1(2), DUMMY,
     2                 diffSV(2,4,4), SDC1, SDC2, SDC3, SQC1, SQC2, SQC3,
     3                 diffS2,R1st(NATOMS,NATOMS),EIGVEC(4*NATOMS,4*NATOMS),
     4                 ALPHAADD(NATOMS,NATOMS),diffK,AADD(2)
C    5                 B1st(NATOMS,NATOMS),BOND1st(NATOMS,3)
      DOUBLE PRECISION LMC(NATOMS,NATOMS),ALPHASUM1S(NATOMS,NATOMS),ALPHASUM2,SVSUB,
     1                 ALPHASUM1V(NATOMS,NATOMS),
     2                 SQDIRCOS(3*NATOMS,NATOMS),SUBDIRCOS(3*NATOMS,NATOMS),LONG(2),OPTION(3)

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
                 DO J=I+1,N
               R1st(J,I)=DIRCOS(I,J,K)*REPULSE(J,I)
               R1st(I,J)=-R1st(J,I)
            END DO
            R1st(I,I)=0.0D0
            D1=0.0D0
            DO J=1,N
               D1=D1+BETA*R1st(J,I)
            END DO
            REP1st(I,K)=D1
         END DO
      END DO

C      WRITE(6,*)'First going to do derivative of repulsive energy:'
C      DO K=1,3
C         DO I=1,N
C            WRITE(6,*) (REP1st(I,K))
C         END DO
C         WRITE(6,*)
C      END DO


C   Easy isn't it.
C   Now do the bond counting term derivative.

C NO BOND COUNTING TERM IN THIS VERSION

C      DO K=1,3
C         DO I=1,N
C            DO J=I+1,N
C               B1st(J,I)=DIRCOS(I,J,K)*BONDAGE(J,I)/(BONDAGE(J,I)+1.0D0)**2
C               B1st(I,J)=-B1st(J,I)
C            END DO
C            B1st(I,I)=0.0D0
C            D1=0.0D0
C            DO J=1,N
C               D1=D1-A*B1st(J,I)/DELTA
C            END DO
C            BOND1st(I,K)=D1
C         END DO
C      END DO
C      WRITE(6,*)'Then going to do derivative of bond energy:'
C      DO K=1,3
C         DO I=1,N
C            WRITE(6,*) BOND1st(I,K)
C         END DO
C         WRITE(6,*)
C      END DO

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

C   Transferable potential changes value of S2 and Kappa...therefore changes
C   derivatives by a factor which is always multiplied by alpha.
C   IMPROVED version has NEW ALPHAADD and diffK, the latter of which has to be changed
C   in x, y and z.

           ALPHAADD(J,I)=ALPHA+2.0D0*SIGMA*(R(J,I)-D)

               ALPHASUM1S(J,I)=(ALPHAADD(J,I)+2.0D0/R(J,I))
               ALPHASUM1V(J,I)=(ALPHA+2.0D0/R(J,I))
        END DO
      END DO

C     WRITE(6,*)'Entering the differentiation of electronic energy zone'
        DO K=1,N
          ELEC1st(K,1)=0.0D0
          ELEC1st(K,2)=0.0D0
          ELEC1st(K,3)=0.0D0
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

           AADD(1)=ALPHAADD(J,I)
           AADD(2)=ALPHA

           AS1(1)=ALPHASUM1S(J,I)
           AS1(2)=ALPHASUM1V(J,I)
           SDC1=SUBDIRCOS(J2+1,I)
           SDC2=SUBDIRCOS(J2+2,I)
           SDC3=SUBDIRCOS(J2+3,I)
           SQC1=SQDIRCOS(J2+1,I)
           SQC2=SQDIRCOS(J2+2,I)
           SQC3=SQDIRCOS(J2+3,I)

               LONG(1)=RR*(2.0D0*SQC1-1.0D0)+SQC1*AADD(1)
               LONG(2)=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

C               OPTION(1)=-D2*LONG
C               OPTION(2)=-D3*LONG
C               OPTION(3)=-LMC(J,I)*AS1

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

               DO SV=1,2
                  ASV=2*AI+4*(SV-1)
                  ALPHASUM2=D1*SVtype(ASV+2,J)*(AADD(SV)+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-AADD(SV)*D1*SVtype(ASV+1,J)
                  diffSV(SV,1,2)=SVtype(ASV+2,J)*(RR*SDC1+AADD(SV)*SQC1)
              diffSV(SV,1,3)=D2*ALPHASUM2
                  diffSV(SV,1,4)=D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*(AADD(SV)*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1(SV))
                  diffSV(SV,2,3)=-SVSUB*D2*LONG(SV)
                  diffSV(SV,2,4)=-SVSUB*D3*LONG(SV)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*SQC2*AS1(SV)-SVtype(ASV+4,J)*(AADD(SV)*SDC2+2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=-SVSUB*LMC(J,I)*AS1(SV)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*SQC3*AS1(SV)-SVtype(ASV+4,J)*(AADD(SV)*SDC3+2.0D0*SQC3*RR))

               END DO ! SV
           diffS2=-AADD(1)*D1*S2(I,J)
C   New value of diffK in IMPROVED version
C           KAPSQ=(1.0D0/KAPPANEW(J,I))**2
C           diffK=2.0D0*SIGMA*D1*(R(J,I)-D)*KAPPANEW(J,I)
C   Newest version where KAPPANEW is reciprocal of previous value, the new 
C   diffK=diffK(old)*KAPSQ
           diffK=2.0D0*SIGMA*D1*(R(J,I)-D)*KAPPANEW(J,I)
           E1=0.0D0
           DO NUM1=1,4
             DO NUM2=1,4   

              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-Horig((AJ+NUM2),(AI+NUM1))
C     1             /LAMBDA(I,J)*(diffK*KAPSQ+2.0D0*S2(I,J)*diffS2)
     1             /LAMBDA(I,J)*(diffK+2.0D0*S2(I,J)*diffS2)

              E2=diffSV(1,NUM1,NUM2)
              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(diffH-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,1)=ELEC1st(I,1)+2.0D0*E1
           ELEC1st(J,1)=ELEC1st(J,1)-2.0D0*E1

         END DO  !J (1st time)
        END DO  !I(1st time)

C   Differentiating by y

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

           AADD(1)=ALPHAADD(J,I)
           AADD(2)=ALPHA

           AS1(1)=ALPHASUM1S(J,I)
           AS1(2)=ALPHASUM1V(J,I)
           SDC1=SUBDIRCOS(J2+2,I)
           SDC2=SUBDIRCOS(J2+3,I)
           SDC3=SUBDIRCOS(J2+1,I)
           SQC1=SQDIRCOS(J2+2,I)
           SQC2=SQDIRCOS(J2+3,I)
           SQC3=SQDIRCOS(J2+1,I)

               LONG(1)=RR*(2.0D0*SQC1-1.0D0)+SQC1*AADD(1)
               LONG(2)=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

C               OPTION(1)=-D2*LONG(SV)
C               OPTION(2)=-D3*LONG(SV)
C               OPTION(3)=-LMC(J,I)*AS1(SV)

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

               DO SV=1,2
                  ASV=2*AI+4*(SV-1)
                  ALPHASUM2=D1*SVtype(ASV+2,J)*(AADD(SV)+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-AADD(SV)*D1*SVtype(ASV+1,J)
                  diffSV(SV,1,3)=SVtype(ASV+2,J)*(RR*SDC1+AADD(SV)*SQC1)
              diffSV(SV,1,4)=D2*ALPHASUM2
                  diffSV(SV,1,2)=D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*(AADD(SV)*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1(SV))
                  diffSV(SV,2,3)=-SVSUB*D3*LONG(SV)
                  diffSV(SV,2,4)=-SVSUB*LMC(J,I)*AS1(SV)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*SQC2*AS1(SV)+SVtype(ASV+4,J)*(-AADD(SV)*SDC2-2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=-SVSUB*D2*LONG(SV)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*SQC3*AS1(SV)+SVtype(ASV+4,J)*(-AADD(SV)*SDC3-2.0D0*SQC3*RR))

               END DO ! SV

           diffS2=-AADD(1)*D1*S2(I,J)
C           KAPSQ=(1.0D0/KAPPANEW(J,I))**2
C   New value of diffK in IMPROVED version
           diffK=2.0D0*SIGMA*D1*(R(J,I)-D)*KAPPANEW(J,I)

           E1=0.0D0
           DO NUM2=1,4   
             DO NUM1=1,4

              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-Horig((AJ+NUM2),(AI+NUM1))
     1             /LAMBDA(I,J)*(diffK+2.0D0*S2(I,J)*diffS2)

              E2=diffSV(1,NUM1,NUM2)
              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(diffH-DIAG(LEV)*E2)
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
           ELEC1st(I,2)=ELEC1st(I,2)+2.0D0*E1
           ELEC1st(J,2)=ELEC1st(J,2)-2.0D0*E1

         END DO  !J (1st time)
        END DO  !I(1st time)

C   Now differentiating by z

        DO I=1,N
         AI=4*(I-1)

         DO J=I+1,N
           RR=1.0D0/R(J,I)
           AJ=4*(J-1)
           J2=3*(J-1)

           D1=DIRCOS(J,I,3)
           D2=DIRCOS(J,I,1)
           D3=DIRCOS(J,I,2)

           AADD(1)=ALPHAADD(J,I)
           AADD(2)=ALPHA

           AS1(1)=ALPHASUM1S(J,I)
           AS1(2)=ALPHASUM1V(J,I)
           SDC1=SUBDIRCOS(J2+3,I)
           SDC2=SUBDIRCOS(J2+1,I)
           SDC3=SUBDIRCOS(J2+2,I)
           SQC1=SQDIRCOS(J2+3,I)
           SQC2=SQDIRCOS(J2+1,I)
           SQC3=SQDIRCOS(J2+2,I)

               LONG(1)=RR*(2.0D0*SQC1-1.0D0)+SQC1*AADD(1)
               LONG(2)=RR*(2.0D0*SQC1-1.0D0)+SQC1*ALPHA

C               OPTION(1)=-D2*LONG(SV)
C               OPTION(2)=-D3*LONG(SV)
C               OPTION(3)=-LMC(J,I)*AS1(SV)

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

               DO SV=1,2
                  ASV=2*AI+4*(SV-1)
                  ALPHASUM2=D1*SVtype(ASV+2,J)*(AADD(SV)+RR)
                  SVSUB=SVtype(ASV+3,J)-SVtype(ASV+4,J)

                  diffSV(SV,1,1)=-AADD(SV)*D1*SVtype(ASV+1,J)
                  diffSV(SV,1,4)=SVtype(ASV+2,J)*(RR*SDC1+AADD(SV)*SQC1)
              diffSV(SV,1,2)=D2*ALPHASUM2
                  diffSV(SV,1,3)=D3*ALPHASUM2

                  diffSV(SV,2,1)=-diffSV(SV,1,2)
                  diffSV(SV,4,4)=-D1*(SVtype(ASV+3,J)*(AADD(SV)*SQC1+2.0D0*RR*SDC1)-SVtype(ASV+4,J)*SDC1*AS1(SV))
                  diffSV(SV,2,3)=-SVSUB*LMC(J,I)*AS1(SV)
                  diffSV(SV,2,4)=-SVSUB*D2*LONG(SV)
              
                  diffSV(SV,3,1)=-diffSV(SV,1,3)
                  diffSV(SV,3,2)= diffSV(SV,2,3)
                  diffSV(SV,2,2)=-D1*(SVtype(ASV+3,J)*SQC2*AS1(SV)+SVtype(ASV+4,J)*(-AADD(SV)*SDC2-2.0D0*SQC2*RR))
                  diffSV(SV,3,4)=-SVSUB*D3*LONG(SV)

                  diffSV(SV,4,1)=-diffSV(SV,1,4)
                  diffSV(SV,4,2)= diffSV(SV,2,4)
                  diffSV(SV,4,3)= diffSV(SV,3,4)
                  diffSV(SV,3,3)=-D1*(SVtype(ASV+3,J)*SQC3*AS1(SV)+SVtype(ASV+4,J)*(-AADD(SV)*SDC3-2.0D0*SQC3*RR))

               END DO ! SV

           diffS2=-AADD(1)*D1*S2(I,J)
C           KAPSQ=(1.0D0/KAPPANEW(J,I))**2
C   New value of diffK in IMPROVED version
           diffK=2.0D0*SIGMA*D1*(R(J,I)-D)*KAPPANEW(J,I)

           E1=0.0D0
           DO NUM2=1,4   
             DO NUM1=1,4

              diffH=diffSV(2,NUM1,NUM2)*LAMBDA(I,J)-Horig((AJ+NUM2),(AI+NUM1))
     1             /LAMBDA(I,J)*(diffK+2.0D0*S2(I,J)*diffS2)

              E2=diffSV(1,NUM1,NUM2)
              DO LEV=4*N,2*N+1,-1
               E1=E1+EIGVEC(LEV,AI+NUM1)*EIGVEC(LEV,AJ+NUM2)*(diffH-DIAG(LEV)*E2)
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

C      DO I=1,3
C        WRITE(6,*) 'Electronic derivative'
C        DO K=1,N
C          WRITE(6,*) ELEC1st(K,I)
C        END DO
C       WRITE(6,*)
C      END DO


C   Include factor of 2 for ELEC1st since two electrons fill each MO
      DO K=1,3
         DO I=1,N
            I2=3*(I-1)
C            deriv1st(I2+K)=2.0D0*ELEC1st(I,K)+BOND1st(I,K)+REP1st(I,K)
            deriv1st(I2+K)=2.0D0*ELEC1st(I,K)+REP1st(I,K)
         END DO
      END DO

C      DO I=1,3*N
C            WRITE(6,*) 'Thus the first derivative of the total energy'
C     1                 ,deriv1st(I)
C      END DO

      END
