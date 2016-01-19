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

      SUBROUTINE DFTB(N,XS,deriv1st,ENERGY,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      INTEGER I, J, K, AI, AJ, NDIM,NMAX, I2, J2, K2,  INFO,
     /        NELEC,SUB(NATOMS),LEVA,LEVB,J1,N
      DOUBLE PRECISION DIST,XMULREAL,WORK(24*natoms),
     2                 S(4*NATOMS,4*NATOMS),Sorig(4*NATOMS,4*NATOMS),
     3                 ENERGY,XS(3*NATOMS),HB(4*NATOMS,4*NATOMS),
     4                 H(4*NATOMS,4*NATOMS),HA(4*NATOMS,4*NATOMS),
     /                 Horiga(4*NATOMS,4*NATOMS),Horigb(4*NATOMS,4*NATOMS),
     /                 VA(4*NATOMS,4*NATOMS),VB(4*NATOMS,4*NATOMS)

      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),DIAG(4*NATOMS),
     3                 deriv1st(3*NATOMS),DIAGA(4*NATOMS),DIAGB(4*NATOMS)

      DOUBLE PRECISION Ssssig,Sspsig,Sppsig,Spppi,
     1                 EPSsac,EPSsbc,EPSpac,EPSpbc,Hsssiga,Hspsiga,
     2        Hppsiga,Hpppia,Hsssigb,Hspsigb,Hppsigb,Hpppib

       PARAMETER(EPSsac=-0.519289D0, EPSpac=-0.216472D0,
     1    EPSsbc=-0.434382D0, EPSpbc=-0.216472D0)


      DOUBLE PRECISION UREP,REP
      LOGICAL DIAGTEST
      COMMON /FAIL/ DIAGTEST

      XMUL=INT(XMULREAL)

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
               DIRCOS(I,I,K)=0.0D0
            END DO
C            PRINT*,'D1 IS',DIRCOS(I,J,1)
C            PRINT*,'D2 IS',DIRCOS(I,J,2)
C            PRINT*,'D3 IS',DIRCOS(I,J,3)

         END DO
      END DO 


C   Calculate Repulsive Free Energy

      UREP=0.0D0

      DO I=1,N
         DO J=I+1,N

CC          PRINT*,'js and atnos are',I,J,IATNUM(I),IATNUM(J)
C REMEMBER  : THESE CALCUALTIONS ARE IN BOHR, NOT ANGSTROM

              CALL AREPCC(N,XS,I,J,REP,R, DIRCOS, DIAGA, DIAGB,natoms)

C        PRINT*,'rep in h is',REP

        UREP=REP+UREP

         REP=0.0D0

        ENDDO
      ENDDO

C      PRINT*,' NEW urep in h is',UREP

C10    FORMAT(A8, F20.17)

C   Only include overlap or interaction with different atoms here.
               
      DO I=1,N
         AI=4*(I-1)
         DO J=I+1,N
            AJ=4*(J-1)

C   Calculate the overlap matrix S(I,J), related to these parameters
C   by the Slater-Koster scheme

          CALL OVECC(N,XS,I,J,Ssssig, Sspsig, Sppsig, Spppi,R, DIRCOS, DIAGA, DIAGB,natoms)


           S((AI+1),(AJ+2))=DIRCOS(I,J,1)*Sspsig
           S((AI+1),(AJ+3))=DIRCOS(I,J,2)*Sspsig
           S((AI+1),(AJ+4))=DIRCOS(I,J,3)*Sspsig

           S((AI+2),(AJ+1))=-DIRCOS(I,J,1)*Sspsig
           S((AI+3),(AJ+1))=-DIRCOS(I,J,2)*Sspsig
           S((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Sspsig


            S((AI+1),(AJ+1))=Ssssig
 
            S((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,1)**2)*Spppi
            S((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Sppsig
     1                       -Spppi)
            S((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Sppsig
     1                       -Spppi)

            S((AI+3),(AJ+2))=S((AI+2),(AJ+3))
            S((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,2)**2)*Spppi
            S((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Sppsig
     1                       -Spppi)

            S((AI+4),(AJ+2))=S((AI+2),(AJ+4))
            S((AI+4),(AJ+3))=S((AI+3),(AJ+4))
            S((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Sppsig+(1.0D0-
     1                       DIRCOS(I,J,3)**2)*Spppi

            
C   Hamiltonian is found using the Slater-Koster scheme
C   as shown in their paper and Harrison's book p.481
C   which relates H to the parameter V using direction cosines.
C   LAMBDA changes for each I and J
      
C my change means constructing  a and b type interaction matrices
C as an uhf analogue

          CALL HAMCCA(N,XS,I,J,Hsssiga, Hspsiga, Hppsiga, Hpppia,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL HAMCCB(N,XS,I,J,Hsssigb, Hspsigb, Hppsigb, Hpppib,R, DIRCOS, DIAGA, DIAGB,natoms)

          HA((AI+1),(AJ+2))=DIRCOS(I,J,1)*Hspsiga
          HA((AI+1),(AJ+3))=DIRCOS(I,J,2)*Hspsiga
          HA((AI+1),(AJ+4))=DIRCOS(I,J,3)*Hspsiga

          HA((AI+2),(AJ+1))=-DIRCOS(I,J,1)*Hspsiga
          HA((AI+3),(AJ+1))=-DIRCOS(I,J,2)*Hspsiga
          HA((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Hspsiga


          HB((AI+1),(AJ+2))=DIRCOS(I,J,1)*Hspsigb
          HB((AI+1),(AJ+3))=DIRCOS(I,J,2)*Hspsigb
          HB((AI+1),(AJ+4))=DIRCOS(I,J,3)*Hspsigb

          HB((AI+2),(AJ+1))=-DIRCOS(I,J,1)*Hspsigb
          HB((AI+3),(AJ+1))=-DIRCOS(I,J,2)*Hspsigb
          HB((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Hspsigb


C  transform the Ha into SK elements....

            HA((AI+1),(AJ+1))=Hsssiga

            HA((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Hppsiga+(1.0D0-
     1           DIRCOS(I,J,1)**2)*Hpppia
            HA((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Hppsiga
     1           -Hpppia)
            HA((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Hppsiga
     1           -Hpppia)

            HA((AI+3),(AJ+2))=HA((AI+2),(AJ+3))
            HA((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Hppsiga+(1.0D0-
     1           DIRCOS(I,J,2)**2)*Hpppia
            HA((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Hppsiga
     1           -Hpppia)

            HA((AI+4),(AJ+2))=HA((AI+2),(AJ+4))
            HA((AI+4),(AJ+3))=HA((AI+3),(AJ+4))
            HA((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Hppsiga+(1.0D0-
     1           DIRCOS(I,J,3)**2)*Hpppia

C now for the beta type...

            HB((AI+1),(AJ+1))=Hsssigb

            HB((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Hppsigb+(1.0D0-
     1           DIRCOS(I,J,1)**2)*Hpppib
            HB((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Hppsigb
     1           -Hpppib)
            HB((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Hppsigb
     1           -Hpppib)

            HB((AI+3),(AJ+2))=HB((AI+2),(AJ+3))
            HB((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Hppsigb+(1.0D0-
     1           DIRCOS(I,J,2)**2)*Hpppib
            HB((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Hppsigb
     1           -Hpppib)

            HB((AI+4),(AJ+2))=HB((AI+2),(AJ+4))
            HB((AI+4),(AJ+3))=HB((AI+3),(AJ+4))
            HB((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Hppsigb+(1.0D0-
     1           DIRCOS(I,J,3)**2)*Hpppib


         END DO

C   Interaction/Overlap of different AOs on the same atom
         
          DO K=1,4
           DO K2=K+1,4
             S((AI+K),(AI+K2))=0.0D0
C             H((AI+K),(AI+K2))=0.0D0

             HA((AI+K),(AI+K2))=0.0D0
             HB((AI+K),(AI+K2))=0.0D0

           END DO
          END DO

C   Interaction of same AO on the same atom

          HA((AI+1),(AI+1))=EPSsac
          HA((AI+2),(AI+2))=EPSpac
          HA((AI+3),(AI+3))=EPSpac
          HA((AI+4),(AI+4))=EPSpac

          HB((AI+1),(AI+1))=EPSsbc
          HB((AI+2),(AI+2))=EPSpbc
          HB((AI+3),(AI+3))=EPSpbc
          HB((AI+4),(AI+4))=EPSpbc

      END DO

C   Now we have to make sure the matrix S and H are symmetric.
C   Ergo, have to get the upper diagonal 'cos must have complete
C   matrix for Cholesky decomposition.
      
      DO I=1,4*N
            DO J=I+1,4*N
            S(J,I)=S(I,J)
C            H(J,I)=H(I,J)
            HA(J,I)=HA(I,J)
            HB(J,I)=HB(I,J)
          END DO
    

C   Interaction/overlap of the same AO on the same atom
         S(I,I)=1.0D0

      END DO

C       PRINT*,'ha is'
C       DO J=1,4*N
C       WRITE(*,50) (HA(I,J),I=1,4*N)
C       ENDDO

C       PRINT*,'hb is'
C       DO J=1,4*N
C       WRITE(*,50) (HB(I,J),I=1,4*N)
C       ENDDO

C           PRINT*,'s is'
C          DO J=1,4*N
C          WRITE(*,50) (S(I,J),I=1,4*N)
C          ENDDO

C50    FORMAT (10(5X,20F8.5,5X))
C50     FORMAT (10F8.4)

      DO I=1,4*N
         DO J=1,4*N
C            Horig(I,J)=H(I,J)
           Horiga(I,J)=HA(I,J)
           Horigb(I,J)=HB(I,J)
         END DO
      END DO

      DO I=1,4*N
        DO J=1,4*N
           Sorig(I,J)=S(I,J)
        END DO
      END DO

C      PRINT*,'before dsy, ha is'
      DO I=1,4*N
         DO J=1,4*N
            H(I,J)=Horiga(I,J)
         END DO
      END DO

C       DO J=1,4*N
C       WRITE(*,50) (H(I,J),I=1,4*N)
C       ENDDO

C      PRINT*,'before dsy, S is'

      DO I=1,4*N
        DO J=1,4*N
           S(I,J)=Sorig(I,J)
        END DO
      END DO

C       DO J=1,4*N
C       WRITE(*,50) (S(I,J),I=1,4*N)
C       ENDDO

      NDIM=4*N
      NMAX=4*NATOMS

      DIAGTEST=.FALSE.

      CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, 24*natoms, INFO )
C
C   This routine is not required for SGI since lapack routine includes it.
C
      IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGSRT(DIAG,H,NDIM,NMAX)
 
      IF (INFO.NE.0) THEN
         DIAGTEST=.TRUE.
         ENERGY=1.0D6
         PRINT*,'alpha fucking DFTB DSYGV failed - INFO=',INFO
         PRINT*,'coords are'
         WRITE(*,80) (XS(J1),J1=1,3*N)
80       FORMAT(3F15.5)
C         CALL DFTB3(XMUL,IATNUM,N,XS,deriv1st,ENERGY,GTEST)
C         RETURN
         STOP
      ENDIF

C      PRINT*,'now the returned alpha evectors are'

C       DO J=1,4*N
C       WRITE(*,50) (H(I,J),I=1,4*N)
C       ENDDO

C writing the evectors to HA
      DO J=1,4*N
         DO I=1,4*N
            VA(J,I)=H(J,I)
         ENDDO
      ENDDO

      DO I=1,4*N
      DIAGA(I)=DIAG(I)
      END DO

C      PRINT*,'done the alpha bit'

C      PRINT*,'before dsy, S is'

      DO I=1,4*N
        DO J=1,4*N
           S(I,J)=Sorig(I,J)
        END DO
      END DO

C       DO J=1,4*N
C       WRITE(*,50) (S(I,J),I=1,4*N)
C       ENDDO

C      PRINT*,'before dsy, hb is'
      DO I=1,4*N
         DO J=1,4*N
            H(I,J)=Horigb(I,J)
         END DO
      END DO
C       DO J=1,4*N
C       WRITE(*,50) (H(I,J),I=1,4*N)
C       ENDDO



      DIAGTEST=.FALSE.

       NDIM=4*N
       NMAX=4*NATOMS


      CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, 24*natoms, INFO )
C
C   This routine is not required for SGI since lapack routine includes it.
C
      IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGSRT(DIAG,H,NDIM,NMAX)

      IF (INFO.NE.0) THEN
         DIAGTEST=.TRUE.
         ENERGY=1.0D6
         PRINT*,'DSYGV failed - INFO=',INFO
         RETURN
      ENDIF

C      PRINT*,'now the returned beta evectors are'

C       DO J=1,4*N
C       WRITE(*,50) (H(I,J),I=1,4*N)
C       ENDDO

C writing the evectors to HB
      DO J=1,4*N
         DO I=1,4*N
            VB(J,I)=H(J,I)
         ENDDO
      ENDDO

      DO I=1,4*N
      DIAGB(I)=DIAG(I)
      ENDDO


C      WRITE(6,*)
C      WRITE(6,*)  'The alpha eigenvalues of the Hessian matrix:'
C      DO I=1,4*N
CC         WRITE(6,*) I,DIAGA(I)
C      END DO
C      WRITE(6,*)

C      WRITE(6,*)
C      WRITE(6,*)  'The beta eigenvalues of the Hessian matrix:'
C      DO I=1,4*N
C         WRITE(6,*)  DIAGB(I)
C      END DO
C      WRITE(6,*)

      ENERGY=0.0D0

C figuring out the occupation numbers
C       NELEC=0
C       PRINT*,'xmul is',XMUL
C       DO J=1,N

C       IF (IATNUM(J).NE.1) THEN
C       SUB(J)=IATNUM(J)-2
C       ELSE
C       SUB(J)=IATNUM(J)
C       ENDIF

C       NELEC=SUB(J)+NELEC 
C       ENDDO

C       LEVA=(NELEC-XMUL+1)/2
C       LEVB=LEVA+(XMUL-1)
C       PRINT*,'NELEC and LEVA/B are',NELEC,LEVA,LEVB
       
      DO J=4*N,((4*N)-120+1),-1
C        DO J=1,LEVA
C       PRINT*,'J is',J
C         PRINT*,'diag of j is',DIAGA(J)
         ENERGY=ENERGY+DIAGA(J)
      END DO
      DO J=4*N,((4*N)-120+1),-1
C        DO J=1,LEVB
C       PRINT*,'J is',J
C         PRINT*,'diag of j is',DIAGB(J)
         ENERGY=ENERGY+DIAGB(J)
      END DO

C        PRINT*,'elec energy is',ENERGY

C   For single atom there should be no contribution from the
C   bond counting term
C      IF (N .EQ. 1) THEN
C         ENERGY=2*ENERGY
C      ELSE
C altered with no urep to figure out what my urep should be...
          ENERGY=ENERGY+UREP
C           ENERGY=UREP+0.0D0
C testing the rep part of the potential (derivs)

C      END IF
C      WRITE(6,*)
C      WRITE(6,*) 'ENERGY=',ENERGY

C   Now call the all new originally crafted subroutine to calculate
C   the first derivatives.

      IF (GTEST) CALL DFDERIV1(XMUL,IATNUM,N,XS,deriv1st,VA,VB,R, DIRCOS, DIAGA, DIAGB,natoms)

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   
      SUBROUTINE OVECC(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),XS(3*NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     3                  D1,D2,D3,D4,D5,D6,A3,B3,
     4                 D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,
     5                 D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,
     6             D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40
C old param :D2= -0.3593880000 ,
C old param: D22=  0.0272058000 ,
      PARAMETER (   D1=  0.4635980000 ,D2=-0.354943D0,
     /  D3=  0.1597240000 ,D4= -0.0236036000
     / ,D5= -0.0160748000 ,D6=  0.0101939000 ,D7= -0.0010466800 ,D8= -0.0013805000
     / ,D9=  0.0007822820,D10= -0.0001806570,
     /  D11= -0.3627180000 ,D12=  0.2504088000 ,D13= -0.0459648000 ,D14= -0.0572412000
     / ,D15=  0.0499554000 ,D16= -0.0161908000 ,D17= -0.0007504350 ,D18=  0.0030270400
     / ,D19= -0.0013407800,D20=  0.0002213130,
     /  D21= -0.1381954000 , D22=0.021232D0,
     /  D23=  0.1362280000 ,D24= -0.1495750000
     / ,D25=  0.0722040000 ,D26= -0.0144966000 ,D27= -0.0040281700 ,D28=  0.0047257800
     / ,D29= -0.0021491900,D30=  0.0004729620,
     /  D31=  0.3653820000 ,D32= -0.3031560000 ,D33=  0.1686080000 ,D34= -0.0580734000
     / ,D35=  0.0075940400 ,D36=  0.0040539700 ,D37= -0.0031628400 ,D38=  0.0011715700
     / ,D39= -0.0002414690,D40= -0.0000131366)

       A3=1.0D0
       B3=7.0D0

       QTB(K,L)=(R(K,L)-((B3+A3)/2.0D0))/((B3-A3)/2.0D0)

C       PRINT*,'r of ij in OVERL is',R(K,L)

C  Ssssig

      IF (R(K,L).LT.6.858869592109767D0) THEN

      SSSs=D1 + D2*QTB(K,L) + D3*(-1.0D0 + 2.0D0*QTB(K,L)**2) 
     /    + D4*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    D5*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4) 
     /    + D6*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + D7*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + D8*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + D9*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + D10*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - D1/2.0D0 - 0.005280155435851808D0
     
      ELSE
       
      SSSs=0.0D0

       ENDIF

C   Sspsig

            IF (R(K,L).LT.6.968036930355212D0) THEN

      SSPs= D11 + D12*QTB(K,L) + D13*(-1.0D0 + 2.0D0*QTB(K,L)**2)
     /    + D14*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    D15*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + D16*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + D17*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + D18*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + D19*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + D20*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - D11/2.0D0 - 0.0007691252996453457D0

      ELSE

      SSPs=0.0D0

       ENDIF

C *****
C when using Fs params, swap pps and pppi
C *****
C Sppsig

            IF (R(K,L).LT.7.031215111280781D0) THEN

      SPPs= D21 + D22*QTB(K,L) + D23*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + D24*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    D25*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + D26*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + D27*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + D28*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + D29*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + D30*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - D21/2.0D0 + 0.004483419251163735D0

          ELSE

        SPPs=0.0D0

         ENDIF

C   Spppi

            IF (R(K,L).LT.6.4433341605186D0) THEN

      SPPp= D31 + D32*QTB(K,L) + D33*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + D34*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    D35*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + D36*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + D37*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + D38*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + D39*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + D40*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - D31/2.0D0 + 0.0006053764295787245D0

          ELSE

      SPPp=0.0D0

          ENDIF

C       PRINT*,'in overlapcc, elements are',SSSs,SSPs,SPPs,SPPp
C       PRINT*,'B3 is',B3

C         IF (R(K,L).GT.B3) THEN
C          SSSs=0.0D0
C          SSPs=0.0D0
C          SPPs=0.0D0
C          SPPp=0.0D0
C         ENDIF

       RETURN
       END
   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   
      SUBROUTINE HAMCCA(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),XS(3*NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     /     G1,G2,G3,
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

C old param G22= -0.1770288000 ,
       PARAMETER (  G1= -0.4546540000 ,G2=  0.3501960000 ,G3= -0.1513000000 ,G4=  0.0159980000
     / ,G5=  0.0247918000 ,G6= -0.0192764000 ,G7=  0.0071180300 ,G8= -0.0002604590
     / ,G9= -0.0014905100,G10=  0.0013685600,
     /  G11=  0.3813740000 ,G12= -0.2675544000 ,G13=  0.0587541000 ,G14=  0.0551394000
     / ,G15= -0.0624710000 ,G16=  0.0366550000 ,G17= -0.0161735000 ,G18=  0.0047799300
     / ,G19= -0.0000637944,G20= -0.0010983800,
     /  G21=  0.2788480000 , G22=-0.1669651D0,
     /  G23=  0.0113421000 ,G24=  0.0677288000
     / ,G25= -0.0683438000 ,G26=  0.0436374000 ,G27= -0.0211787000 ,G28=  0.0074277200
     / ,G29= -0.0013746800,G30= -0.0004145400,
     /  G31= -0.3851400000 ,G32=  0.3337080000 ,G33= -0.2141850000 ,G34=  0.1061410000
     / ,G35= -0.0414461000 ,G36=  0.0121444000 ,G37= -0.0017954400 ,G38= -0.0007822330
     / ,G39=  0.0008151550,G40= -0.0004982710)


       A2=1.0D0
       B2=7.0D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C  Ssssig

       IF (R(K,L).LT.6.501645139158383D0) THEN

      SSSs=G1 + G2*QTB(K,L) + G3*(-1.0D0 + 2.0D0*QTB(K,L)**2) 
     /    + G4*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G5*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4) 
     /    + G6*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G7*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G8*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G9*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G10*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G1/2.0D0 + 0.001465818855739154D0 

       ELSE
 
       SSSs=0.0D0

       ENDIF

C   Sspsig

       IF (R(K,L).LT.6.514708000907687D0) THEN

      SSPs=G11 + G12*QTB(K,L) + G13*(-1.0D0 + 2.0D0*QTB(K,L)**2)
     /    + G14*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G15*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G16*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G17*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G18*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G19*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G20*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G11/2.0D0 + 0.0001542094201995026D0

       ELSE
 
      SSPs=0.0D0

       ENDIF

C Sppsig

       IF (R(K,L).LT.6.699591659971658D0) THEN

      SPPs=G21 + G22*QTB(K,L) + G23*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + G24*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G25*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G26*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G27*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G28*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G29*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G30*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G21/2.0D0 - 0.01228748543456073D0

       ELSE

      SPPs=0.0D0

        ENDIF


C   Spppi

       IF (R(K,L).LT.6.841936337171249D0) THEN

      SPPp=G31 + G32*QTB(K,L) + G33*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + G34*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G35*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G36*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G37*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G38*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G39*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G40*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G31/2.0D0 - 0.001936846016333094D0 

         ELSE

            SPPp=0.0D0

         ENDIF


C         IF (R(K,L).GT.B2) THEN
C          SSSs=0.0D0
C          SSPs=0.0D0
C          SPPs=0.0D0
C          SPPp=0.0D0
C         ENDIF

      RETURN
      END

C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE HAMCCB(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     /     G1,G2,G3,XS(3*NATOMS),
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

       PARAMETER (   G1= -0.3789880000 , G2=0.287963,
C old param
C G2=  0.2886558000 ,
     /  G3= -0.1201660000 ,G4=  0.0100538000
     / ,G5=  0.0203983000 ,G6= -0.0149822000 ,G7=  0.0056378700 ,G8= -0.0006071250
     / ,G9= -0.0008101620,G10=  0.0006987440,
     /  G11=  0.2956740000 ,G12= -0.2010048000 ,G13=  0.0327500000 ,G14=  0.0514607000
     / ,G15= -0.0508434000 ,G16=  0.0288404000 ,G17= -0.0133668000 ,G18=  0.0046855100
     / ,G19=  0.0006011540,G20= -0.0000060000,
     /  G21=  0.2314380000 ,G22= -0.1469886000 ,G23=  0.0151875000 ,G24=  0.0440803000
     / ,G25= -0.0458782000 ,G26=  0.0312287000 ,G27= -0.0171949000 ,G28=  0.0079412400
     / ,G29= -0.0028890500,G30=  0.0012923000,
     /  G31= -0.3325740000 ,G32=  0.2847522000 ,G33= -0.1816530000 ,G34=  0.0901204000
     / ,G35= -0.0360385000 ,G36=  0.0113274000 ,G37= -0.0021334300 ,G38= -0.0006478810
     / ,G39=  0.0009400840,G40= -0.0009030410)

       A2=1.0D0
       B2=7.0D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C  Ssssig

       IF (R(K,L).LT.6.695259436126274D0) THEN

      SSSs=G1 + G2*QTB(K,L) + G3*(-1.0D0 + 2.0D0*QTB(K,L)**2) 
     /    + G4*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G5*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4) 
     /    + G6*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G7*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G8*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G9*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G10*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G1/2.0D0 + 0.002109343015707699D0

      ELSE

      SSSs=0.0D0

      ENDIF


C   Sspsig

       IF (R(K,L).LT.6.796735011037814D0) THEN

      SSPs=G11 + G12*QTB(K,L) + G13*(-1.0D0 + 2.0D0*QTB(K,L)**2)
     /    + G14*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G15*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G16*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G17*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G18*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G19*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G20*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G11/2.0D0 + 0.001551075520518252D0

         ELSE

        SSPs=0.0D0

         ENDIF


C Sppsig

       IF (R(K,L).LT.6.880535598840217D0) THEN

      SPPs=G21 + G22*QTB(K,L) + G23*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + G24*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G25*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G26*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G27*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G28*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G29*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G30*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G21/2.0D0 - 0.001666042658392236D0

         ELSE

         SPPs=0.0D0

         ENDIF


C   Spppi

       IF (R(K,L).LT.5.944816716454667D0) THEN

      SPPp=G31 + G32*QTB(K,L) + G33*(-1.0D0 + 2.0D0*QTB(K,L)**2)        
     /    + G34*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    G35*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4)
     /    + G36*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + G37*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + G38*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + G39*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + G40*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - G31/2.0D0 + 0.0007476927053454984D0

        ELSE

          SPPp=0.0D0

        ENDIF


C         IF (R(K,L).GT.B2) THEN
C          SSSs=0.0D0
C          SSPs=0.0D0
C          SPPs=0.0D0
C          SPPp=0.0D0
C         ENDIF

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE AREPCC(N,XS,K,L,rep,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION rep,
     /     C1,C2,C3,XS(3*NATOMS),
     7     C4,C5,C6,C7,C8,C9,C10,
     9     A2,B2

       PARAMETER (   C1=  2.5911600000 ,C2= -2.0924020000 ,C3=  1.1228900000 ,C4= -0.4302620000
     / ,C5=  0.1257390000 ,C6= -0.0215445000 ,C7=  0.0000000000 ,C8=  0.0000000000
     / ,C9=  0.0000000000,C10=  0.0000000000) 

       A2=1.0D0
       B2=3.3D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C  Ssssig

      rep=C1 + C2*QTB(K,L) + C3*(-1.0D0 + 2.0D0*QTB(K,L)**2) 
     /    + C4*(-3.0D0*QTB(K,L) + 4.0D0*QTB(K,L)**3) +
     /    C5*(1 -8.0D0*QTB(K,L)**2 + 8.0D0*QTB(K,L)**4) 
     /    + C6*(5.0D0*QTB(K,L) -20.0D0*QTB(K,L)**3 + 16.0D0*QTB(K,L)**5)
     /     + C7*(-1.0D0+18.0D0*QTB(K,L)**2 -48.0D0*QTB(K,L)**4 + 32.0D0*QTB(K,L)**6)
     /   + C8*(-7.0D0*QTB(K,L) + 56.0D0*QTB(K,L)**3 - 112.0D0*QTB(K,L)**5 + 64.0D0*QTB(K,L)**7)
     /    + C9*(1.0D0 - 32.0D0*QTB(K,L)**2 + 160.0D0*QTB(K,L)**4 - 256.0D0*QTB(K,L)**6
     /   +128.0D0*QTB(K,L)**8) + C10*(9.0D0*QTB(K,L) - 120.0D0*QTB(K,L)**3 + 432.0D0*QTB(K,L)**5
     /   - 576.0D0*QTB(K,L)**7 + 256.0D0*QTB(K,L)**9) - C1/2.0D0




         IF (R(K,L).GT.3.251838574097556D0) THEN
            rep=0.0D0
         ELSE
           rep=rep+0.00002383611355998103D0
         ENDIF

      RETURN
      END
C   This is phase two where we calculate the first derivatives
C   of the silicon potential coded in TB.Si.f. It is in tight
C   binding parameterised form.
C   Off I go...
C   This is derivs1.f.jloop.betterish in the resrves directory
C   for the record

      SUBROUTINE DFDERIV1(XMUL,IATNUM,N,XS,deriv1st,EIGVECA,EIGVECB,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE 
      INTEGER N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),DERIV1ST(3*NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      INTEGER I, J, K , AI, AJ, LEV, SV, NUM1, NUM2, 
     1        I2,J2,IATNUM(NATOMS),XMUL,NELEC,SUB(NATOMS),NOCCB,NOCCA
      DOUBLE PRECISION REP1st(NATOMS,3),RR, D1, D2, D3,
     1                 E2, AS1, DUMMY,
     2                 diffSV(3,4,4), SDC1, SDC2, SDC3, SQC1, SQC2, SQC3,
     3                 R1st(NATOMS,NATOMS),
C made urep first a scalar
     4                 UREP1st,XS(3*NATOMS),E1A,E1B,
     5                 ELEC1sta(NATOMS,3),ELEC1stb(NATOMS,3),EIGVECA(4*NATOMS,4*NATOMS),
     6                 EIGVECB(4*NATOMS,4*NATOMS)
      DOUBLE PRECISION LMC(NATOMS,NATOMS),ALPHASUM1(NATOMS,NATOMS),SVSUB,
     1                 SQDIRCOS(3*NATOMS,NATOMS),SUBDIRCOS(3*NATOMS,NATOMS),LONG,OPTION(3),
     2                 DELD1, DELD12, DELD2, DELD22, DELD3, DELD32, DELRDX,DRSss,DRSsp,
     3                 DRSpps,DRSppp,DSVSUB, Ssssig, Sspsig, Sppsig, Spppi,
     4                 Hsssig, Hspsig, Hppsig, Hpppi

C       PRINT*,'makes it to dfderiv'
C
C  Take the transpose of EIGVEC - makes loops more efficient
C
      DO J=1,4*N
         DO I=J+1,4*N
            DUMMY=EIGVECA(J,I)
            EIGVECA(J,I)=EIGVECA(I,J)
            EIGVECA(I,J)=DUMMY

            DUMMY=EIGVECB(J,I)
            EIGVECB(J,I)=EIGVECB(I,J)
            EIGVECB(I,J)=DUMMY
         ENDDO
      ENDDO

C50     FORMAT (10F12.8)

C       PRINT*,'the a vectors in dfderiv are'
C       DO J=1,4*N
C       WRITE(*,50) (EIGVECA(I,J),I=1,4*N)
C       ENDDO

C       PRINT*,'the b vectors in dfderiv are'
C       DO J=1,4*N
C       WRITE(*,50) (EIGVECB(I,J),I=1,4*N)
C       ENDDO

C occupation numbers

C       NELEC=0
C       PRINT*,'xmul is',XMUL
C       DO J=1,N

C       IF (IATNUM(J).NE.1) THEN
C       SUB(J)=IATNUM(J)-2
C       ELSE
C       SUB(J)=IATNUM(J)
C       ENDIF

C       NELEC=SUB(J)+NELEC
C       ENDDO

C       NOCCA=(NELEC-XMUL+1)/2
C       NOCCB=NOCCA+(XMUL-1)
C       PRINT*,'in derivs NELEC and LEVA/B are',NELEC,NOCCA,NOCCB

C
C   Can put arrays into a block data format instead of parameter
C   so can just put in once in a file to be included

C   Can't start an array with a no. ie. 1stderiv...syntax error

C   First of all calculate derivatives of repulsive energy term
C   This is called REP1st.

C      WRITE(6,*)'First going to do derivative of repulsive energy:'

      DO K=1,3
         DO I=1,N
            DO J=I+1,N

              CALL DREPCC(N,XS,I,J,UREP1st,natoms,R, DIRCOS, DIAGA, DIAGB)

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
          ELEC1sta(K,1)=0.0D0
          ELEC1stb(K,1)=0.0D0
        ENDDO

        DO I=1,N
         AI=4*(I-1)

C   Diff of interactions of AOs on same atom
C   Don't need to include for diffSV since not used - loop
C   is only for J=I+1,N

         DO J=I+1,N
           RR=1.0D0/R(J,I)
              IF ((I.EQ.1).AND.(J.EQ.26)) THEN
              ENDIF

              IF ((I.EQ.1).AND.(J.EQ.27)) THEN
              ENDIF

           AJ=4*(J-1)
           J2=3*(J-1)

C   A good idea appears to be simplify the expressions in the
C   differentiating section and reduce the number of indices


           D1=DIRCOS(J,I,1)
           D2=DIRCOS(J,I,2)
           D3=DIRCOS(J,I,3)

           SDC1=SUBDIRCOS(J2+1,I)
           SDC2=SUBDIRCOS(J2+2,I)
           SDC3=SUBDIRCOS(J2+3,I)
           SQC1=SQDIRCOS(J2+1,I)
           SQC2=SQDIRCOS(J2+2,I)
           SQC3=SQDIRCOS(J2+3,I)


C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

           DO SV=1,3
              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

               
C              PRINT*,'XTERMS'

           IF (SV.EQ.1) THEN

          CALL OVECC(N,XS,I,J,Ssssig, Sspsig, Sppsig, Spppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DOVECC(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)

 
           diffSV(SV,1,2)=(Sspsig*DELD1 - D1*DRSsp*DELRDX)
            diffSV(SV,1,3)=(Sspsig*DELD2 - D2*DRSsp*DELRDX)
           diffSV(SV,1,4)=(Sspsig*DELD3 - D3*DRSsp*DELRDX)

           diffSV(SV,2,1)=-(Sspsig*DELD1 - D1*DRSsp*DELRDX)
           diffSV(SV,3,1)=-(Sspsig*DELD2 - D2*DRSsp*DELRDX)
           diffSV(SV,4,1)=-(Sspsig*DELD3 - D3*DRSsp*DELRDX)

             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp

              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,2,2)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSV(SV,2,4)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
          
              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,3,3)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,4,4)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX


              ELSE IF (SV.EQ.2) THEN

          CALL HAMCCA(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCA(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)

          diffSV(SV,1,2)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,1,3)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,1,4)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,2,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,3,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,4,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)


             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,2,2)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSV(SV,2,4)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,3,3)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,4,4)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

            ELSE


          CALL HAMCCB(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCB(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


          diffSV(SV,1,2)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,1,3)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,1,4)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,2,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,3,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,4,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)



             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,2,2)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSV(SV,2,4)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,3,3)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,4,4)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX


              ENDIF

           END DO ! SV

C      WRITE(6,*) 'X doves'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(1,NUM1,NUM2)
C        END DO
C        WRITE(6,*)
C      END DO
 
C      WRITE(6,*) 'XA dhams'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(2,NUM1,NUM2)
C         END DO
C         WRITE(6,*)
C      END DO
 
C      WRITE(6,*) 'XB dhams'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(3,NUM1,NUM2)
C         END DO
C        WRITE(6,*)
C      END DO

           E1A=0.0D0
           E1B=0.0D0

           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)

C you'll have to change levels for diff systems until you get a system going...
C this is for CH

              DO LEV=4*N,(4*N-120+1),-1
C               LEV=4*N
C               DO LEV=1,NOCCA
               E1A=E1A+EIGVECA(LEV,AI+NUM1)*EIGVECA(LEV,AJ+NUM2)*(diffSV(2,NUM1,NUM2)-DIAGA(LEV)*E2)
              ENDDO !LEV

              DO LEV=4*N,(4*N-120+1),-1
C               DO LEV=1,NOCCB
               E1B=E1B+EIGVECB(LEV,AI+NUM1)*EIGVECB(LEV,AJ+NUM2)*(diffSV(3,NUM1,NUM2)-DIAGB(LEV)*E2)
              END DO !LEV

             END DO !NUM2
           END DO !NUM1
C later on you'll have to alter the 2 factor to allow for spin up/down occupation for Halpha/Hbeta

           ELEC1sta(I,1)=ELEC1sta(I,1)+2.0D0*E1A
           ELEC1sta(J,1)=ELEC1sta(J,1)-2.0D0*E1A

           ELEC1stb(I,1)=ELEC1stb(I,1)+2.0D0*E1B
           ELEC1stb(J,1)=ELEC1stb(J,1)-2.0D0*E1B

         END DO  !J (1st time)
        END DO  !I(1st time)

        DO K=1,N
          ELEC1sta(K,2)=0.0D0
          ELEC1stb(K,2)=0.0D0
        ENDDO


        DO I=1,N
         AI=4*(I-1)
         I2=3*(I-1)


C   Diff of interactions of AOs on same atom
C   Don't need to include for diffSV since not used - loop
C   is only for J=I+1,N

         DO J=I+1,N
           J2=3*(J-1)

           RR=1.0D0/R(J,I)

           AJ=4*(J-1)

           D1=DIRCOS(J,I,2)
           D2=DIRCOS(J,I,3)
           D3=DIRCOS(J,I,1)

           SDC1=SUBDIRCOS(J2+2,I)
           SDC2=SUBDIRCOS(J2+3,I)
           SDC3=SUBDIRCOS(J2+1,I)
           SQC1=SQDIRCOS(J2+2,I)
           SQC2=SQDIRCOS(J2+3,I)
           SQC3=SQDIRCOS(J2+1,I)


C   Differentiation of overlap or interaction between orbitals
C   on different atoms

C          PRINT*,'YTERMS'
           DO SV=1,3

              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1


           IF (SV.EQ.1) THEN


          CALL OVECC(N,XS,I,J,Ssssig, Sspsig, Sppsig, Spppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DOVECC(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


          diffSV(SV,1,3)=(Sspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,1,4)=(Sspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,1,2)=(Sspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,3,1)=-(Sspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,4,1)=-(Sspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,2,1)=-(Sspsig*DELD3 - D3*DRSsp*DELRDX)


             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSV(SV,2,4)=D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX


              ELSE  IF (SV.EQ.2) THEN

          CALL HAMCCA(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCA(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


            diffSV(SV,1,3)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
            diffSV(SV,1,4)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
            diffSV(SV,1,2)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

            diffSV(SV,3,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
            diffSV(SV,4,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
            diffSV(SV,2,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)


             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSV(SV,2,4)=D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

              ELSE

          CALL HAMCCB(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCB(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


          diffSV(SV,1,3)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
           diffSV(SV,1,4)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,1,2)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,3,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,4,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,2,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)

             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSV(SV,2,4)=D2*D3*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

              ENDIF

           END DO ! SV

C      WRITE(6,*) 'Y doves'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(1,NUM1,NUM2)
C         END DO
C         WRITE(6,*)
C      END DO

C      WRITE(6,*) 'Y dhams'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(2,NUM1,NUM2)
C         END DO
C         WRITE(6,*)
C      END DO

           E1A=0.0D0
           E1B=0.0D0

           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)

C you'll have to change levels for diff systems until you get a system going...
C altered for CH

              DO LEV=4*N,(4*N-120+1),-1
C               LEV=4*N
C                DO LEV=1,NOCCA
               E1A=E1A+EIGVECA(LEV,AI+NUM1)*EIGVECA(LEV,AJ+NUM2)*(diffSV(2,NUM1,NUM2)-DIAGA(LEV)*E2)
              ENDDO ! LEV
               DO LEV=4*N,(4*N-120+1),-1
C               DO LEV=1,NOCCB
               E1B=E1B+EIGVECB(LEV,AI+NUM1)*EIGVECB(LEV,AJ+NUM2)*(diffSV(3,NUM1,NUM2)-DIAGB(LEV)*E2) 
              END DO !LEV
             END DO !NUM2
           END DO !NUM1
C later on you'll have to alter the 2 factor to allow for spin up/down occupation for Halpha/Hbeta


           ELEC1sta(I,2)=ELEC1sta(I,2)+2.0D0*E1A
           ELEC1sta(J,2)=ELEC1sta(J,2)-2.0D0*E1A

           ELEC1stb(I,2)=ELEC1stb(I,2)+2.0D0*E1B
           ELEC1stb(J,2)=ELEC1stb(J,2)-2.0D0*E1B

C           PRINT*,'elec1sta for y is',ELEC1sta(I,2)
C           PRINT*,'elec1stb for y is',ELEC1stb(I,2)

         END DO  !J (1st time)
        END DO  !I(1st time)

        DO K=1,N
          ELEC1sta(K,3)=0.0D0
          ELEC1stb(K,3)=0.0D0
        ENDDO

        DO I=1,N
         AI=4*(I-1)
         I2=3*(I-1)


C   Diff of interactions of AOs on same atom
C   Don't need to include for diffSV since not used - loop
C   is only for J=I+1,N

         DO J=I+1,N
           J2=3*(J-1)

           RR=1.0D0/R(J,I)
           AJ=4*(J-1)

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

           OPTION(1)=-D2*LONG
           OPTION(2)=-D3*LONG
           OPTION(3)=-LMC(J,I)*AS1

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

           DO SV=1,3
              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1


              IF (SV.EQ.1) THEN

          CALL OVECC(N,XS,I,J,Ssssig, Sspsig, Sppsig, Spppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DOVECC(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)



           diffSV(SV,1,4)=(Sspsig*DELD1 - D1*DRSsp*DELRDX)
           diffSV(SV,1,2)=(Sspsig*DELD2 - D2*DRSsp*DELRDX)
           diffSV(SV,1,3)=(Sspsig*DELD3 - D3*DRSsp*DELRDX)

           diffSV(SV,4,1)=-(Sspsig*DELD1 - D1*DRSsp*DELRDX)
           diffSV(SV,2,1)=-(Sspsig*DELD2 - D2*DRSsp*DELRDX)
           diffSV(SV,3,1)=-(Sspsig*DELD3 - D3*DRSsp*DELRDX)


             SVSUB=Sppsig - Spppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,4,4)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D3*D2*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2
              diffSV(SV,2,4)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,2,2)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,3,3)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

              ELSE IF (SV.EQ.2) THEN

          CALL HAMCCA(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCA(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


          diffSV(SV,1,4)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,1,2)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,1,3)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,4,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,2,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,3,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)

             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp


              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,4,4)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D3*D2*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2
              diffSV(SV,2,4)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,2,2)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,3,3)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

              ELSE


          CALL HAMCCB(N,XS,I,J,Hsssig, Hspsig, Hppsig, Hpppi,R, DIRCOS, DIAGA, DIAGB,natoms)
          CALL DAMCCB(N,XS,I,J,DRSss,DRSsp,DRSpps,DRSppp,R, DIRCOS, DIAGA, DIAGB,natoms)


          diffSV(SV,1,4)=(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,1,2)=(Hspsig*DELD2 - D2*DRSsp*DELRDX)
           diffSV(SV,1,3)=(Hspsig*DELD3 - D3*DRSsp*DELRDX)

          diffSV(SV,4,1)=-(Hspsig*DELD1 - D1*DRSsp*DELRDX)
          diffSV(SV,2,1)=-(Hspsig*DELD2 - D2*DRSsp*DELRDX)
          diffSV(SV,3,1)=-(Hspsig*DELD3 - D3*DRSsp*DELRDX)

             SVSUB=Hppsig - Hpppi
             DSVSUB=DRSpps - DRSppp

              diffSV(SV,1,1)=DRSss*DELRDX

              diffSV(SV,4,4)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSppp*DELRDX
              diffSV(SV,2,3)=D3*D2*DELRDX*DSVSUB - SVSUB*D2*DELD3 - SVSUB*D3*DELD2
              diffSV(SV,2,4)=D1*D2*DELRDX*DSVSUB - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSV(SV,3,2)= diffSV(SV,2,3)
              diffSV(SV,2,2)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSppp*DELRDX
              diffSV(SV,3,4)= D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSV(SV,4,2)= diffSV(SV,2,4)
              diffSV(SV,4,3)= diffSV(SV,3,4)
              diffSV(SV,3,3)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSppp*DELRDX

              ENDIF

           END DO ! SV

C printing out the matrix

C      WRITE(6,*) 'Z doves'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(1,NUM1,NUM2)
C         END DO
C        WRITE(6,*)
C      END DO

C      WRITE(6,*) 'Z alpha dhams'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(2,NUM1,NUM2)
C         END DO
C         WRITE(6,*)
C      END DO

C      WRITE(6,*) 'Z beta dhams'
C      DO NUM1 =1,4
C         DO NUM2=1,4
C            WRITE(6,*)diffSV(3,NUM1,NUM2)
C         END DO
C         WRITE(6,*)
C      END DO

          E1A=0.0D0
          E1B=0.0D0

           DO NUM1=1,4
             DO NUM2=1,4

              E2=diffSV(1,NUM1,NUM2)

C you'll have to change levels for diff systems until you get a system going...

              DO LEV=4*N,(4*N-120+1),-1
C               LEV=4*N
C               DO LEV=1,NOCCA
               E1A=E1A+EIGVECA(LEV,AI+NUM1)*EIGVECA(LEV,AJ+NUM2)*(diffSV(2,NUM1,NUM2)-DIAGA(LEV)*E2)
              ENDDO ! LEV

               DO LEV=4*N,(4*N-120+1),-1 
C                DO LEV=1,NOCCB
               E1B=E1B+EIGVECB(LEV,AI+NUM1)*EIGVECB(LEV,AJ+NUM2)*(diffSV(3,NUM1,NUM2)-DIAGB(LEV)*E2)
              END DO !LEV

             END DO !NUM2
           END DO !NUM1
C later on you'll have to alter the 2 factor to allow for spin up/down occupation for Halpha/Hbeta

           ELEC1sta(I,3)=ELEC1sta(I,3)+2.0D0*E1A
           ELEC1sta(J,3)=ELEC1sta(J,3)-2.0D0*E1A

           ELEC1stb(I,3)=ELEC1stb(I,3)+2.0D0*E1B
           ELEC1stb(J,3)=ELEC1stb(J,3)-2.0D0*E1B

C           PRINT*,'elec1st for at 1 z is',ELEC1st(I,3)

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
C            deriv1st(I2+K)=2.0D0*ELEC1st(I,K)+REP1st(I,K)
C             deriv1st(I2+K)=2.0D0*ELEC1st(I,K)
              deriv1st(I2+K)=ELEC1sta(I,K)+ELEC1stb(I,K)+REP1st(I,K)
C              deriv1st(I2+K)=REP1st(I,K)
C              deriv1st(I2+K)=ELEC1sta(I,K)+ELEC1stb(I,K)
C            PRINT*,'i,k and der1st',I,K,deriv1st(I2+K)
         END DO
      END DO
       
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE DOVECC(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,
     3                  D1,D2,D3,D4,D5,D6,A3,B3,XS(3*NATOMS),
     4                 D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,
     5                 D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,
     6             D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40,DQDR

      PARAMETER ( D1=  0.4635980000 , D2=-0.354943D0,
C old params D2= -0.3593880000 ,
C old params D22=  0.0272058000 ,
     /  D3=  0.1597240000 ,D4= -0.0236036000
     / ,D5= -0.0160748000 ,D6=  0.0101939000 ,D7= -0.0010466800 ,D8= -0.0013805000
     / ,D9=  0.0007822820,D10= -0.0001806570,
     /  D11= -0.3627180000 ,D12=  0.2504088000 ,D13= -0.0459648000 ,D14= -0.0572412000
     / ,D15=  0.0499554000 ,D16= -0.0161908000 ,D17= -0.0007504350 ,D18=  0.0030270400
     / ,D19= -0.0013407800,D20=  0.0002213130,
     /  D21= -0.1381954000 , D22=0.021232D0,
     /  D23=  0.1362280000 ,D24= -0.1495750000
     / ,D25=  0.0722040000 ,D26= -0.0144966000 ,D27= -0.0040281700 ,D28=  0.0047257800
     / ,D29= -0.0021491900,D30=  0.0004729620,
     /  D31=  0.3653820000 ,D32= -0.3031560000 ,D33=  0.1686080000 ,D34= -0.0580734000
     / ,D35=  0.0075940400 ,D36=  0.0040539700 ,D37= -0.0031628400 ,D38=  0.0011715700
     / ,D39= -0.0002414690,D40= -0.0000131366)

       A3=1.0D0
       B3=7.0D0

       QTB(K,L)=(R(K,L)-((B3+A3)/2.0D0))/((B3-A3)/2.0D0)

C         PRINT*,'r IN DEROV is',R(K,L)

      DQDR=1.0D0 /((B3-A3)/2.0D0)

C  Ssssig

            IF (R(K,L).LT.6.858869592109767D0) THEN

        SSSs=(D2 + D3*4.0D0*QTB(K,L)
     /   +D4*(-3.0D0 + 12.0D0*(QTB(K,L)**2))
     /   + D5*(-16.0D0*QTB(K,L) + 32.0D0*(QTB(K,L)**3)) +D6*(5.0D0 -60.0D0*(QTB(K,L)**2)
     /   + 80.0D0*(QTB(K,L)**4)) +D7*(36.0D0*QTB(K,L) - 192.0D0*(QTB(K,L)**3) +192.0D0*(QTB(K,L)**5))
     /   +D8*(-7.0D0 +168.0D0*(QTB(K,L)**2) - 560.0D0*(QTB(K,L)**4) +448.0D0*(QTB(K,L)**6))
     /   +D9*(-64.0D0*QTB(K,L) +640.0D0*(QTB(K,L)**3) -1536.0D0*(QTB(K,L)**5) +1024.0D0*(QTB(K,L)**7))
     / +D10*(9.0D0 -360.0D0*(QTB(K,L)**2) +2160*(QTB(K,L)**4)
     /   -4032.0D0*(QTB(K,L)**6) +2304.0D0*(QTB(K,L)**8)))*DQDR

         ELSE
       
      SSSs=0.0D0

       ENDIF

C   Sspsig

            IF (R(K,L).LT.6.968036930355212D0) THEN

      SSPs=D12*DQDR + D13*4.0D0*QTB(K,L)*DQDR
     /   +D14*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + D15*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +D16*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +D17*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +D18*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +D19*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +D20*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR - 4032*(QTB(K,L)**6)*DQDR
     /   +2304*(QTB(K,L)**8)*DQDR)

        ELSE

      SSPs=0.0D0

       ENDIF

C *****
C when using Fs params, swap pps and pppi
C *****

C Sppsig

            IF (R(K,L).LT.7.031215111280781D0) THEN

      SPPs=(D22 + D23*4.0D0*QTB(K,L)
     /   +D24*(-3.0D0 + 12.0D0*(QTB(K,L)**2))
     /   + D25*(-16.0D0*QTB(K,L) + 32.0D0*(QTB(K,L)**3)) +D26*(5.0D0 -60.0D0*(QTB(K,L)**2)
     /   + 80.0D0*(QTB(K,L)**4)) +D27*(36.0D0*QTB(K,L) - 192.0D0*(QTB(K,L)**3) +192.0D0*(QTB(K,L)**5))
     /   +D28*(-7.0D0 +168.0D0*(QTB(K,L)**2) - 560.0D0*(QTB(K,L)**4) +448.0D0*(QTB(K,L)**6))
     /   +D29*(-64.0D0*QTB(K,L) +640.0D0*(QTB(K,L)**3) -1536.0D0*(QTB(K,L)**5) +1024.0D0*(QTB(K,L)**7))
     /   +D30*(9.0D0 -360.0D0*(QTB(K,L)**2) +2160*(QTB(K,L)**4)
     /   -4032.0D0*(QTB(K,L)**6) +2304.0D0*(QTB(K,L)**8)))*DQDR

          ELSE

        SPPs=0.0D0

         ENDIF

C   Spppi

            IF (R(K,L).LT.6.4433341605186D0) THEN

      SPPp=(D32 + D33*4.0D0*QTB(K,L)
     /   +D34*(-3.0D0 + 12.0D0*(QTB(K,L)**2))
     /   + D35*(-16.0D0*QTB(K,L) + 32.0D0*(QTB(K,L)**3)) +D36*(5.0D0 -60.0D0*(QTB(K,L)**2)
     /   + 80.0D0*(QTB(K,L)**4)) +D37*(36.0D0*QTB(K,L) - 192.0D0*(QTB(K,L)**3) +192.0D0*(QTB(K,L)**5))
     /   +D38*(-7.0D0 +168.0D0*(QTB(K,L)**2) - 560.0D0*(QTB(K,L)**4) +448.0D0*(QTB(K,L)**6))
     /   +D39*(-64.0D0*QTB(K,L) +640.0D0*(QTB(K,L)**3) -1536.0D0*(QTB(K,L)**5) +1024.0D0*(QTB(K,L)**7))
     /   +D40*(9.0D0 -360.0D0*(QTB(K,L)**2) +2160*(QTB(K,L)**4)
     /   -4032.0D0*(QTB(K,L)**6) +2304.0D0*(QTB(K,L)**8)))*DQDR


          ELSE

        SPPp=0.0D0

          ENDIF

C   IF (R(K,L).GT.B3) THEN
C          SSSs=0.0D0
C          SSPs=0.0D0
C          SPPs=0.0D0
c          SPPp=0.0D0
c         ENDIF

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE DAMCCA(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,DQDR,
     /     G1,G2,G3,XS(3*NATOMS),
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

C old param: G22= -0.1770288000 ,
       PARAMETER ( G1= -0.4546540000 ,G2=  0.3501960000 ,G3= -0.1513000000 ,G4=  0.0159980000
     / ,G5=  0.0247918000 ,G6= -0.0192764000 ,G7=  0.0071180300 ,G8= -0.0002604590
     / ,G9= -0.0014905100,G10=  0.0013685600,
     /  G11=  0.3813740000 ,G12= -0.2675544000 ,G13=  0.0587541000 ,G14=  0.0551394000
     / ,G15= -0.0624710000 ,G16=  0.0366550000 ,G17= -0.0161735000 ,G18=  0.0047799300
     / ,G19= -0.0000637944,G20= -0.0010983800,
     /  G21=  0.2788480000 ,G22=-0.1669651D0,
     /  G23=  0.0113421000 ,G24=  0.0677288000
     / ,G25= -0.0683438000 ,G26=  0.0436374000 ,G27= -0.0211787000 ,G28=  0.0074277200
     / ,G29= -0.0013746800,G30= -0.0004145400,
     /  G31= -0.3851400000 ,G32=  0.3337080000 ,G33= -0.2141850000 ,G34=  0.1061410000
     / ,G35= -0.0414461000 ,G36=  0.0121444000 ,G37= -0.0017954400 ,G38= -0.0007822330
     / ,G39=  0.0008151550,G40= -0.0004982710)

       A2=1.0D0
       B2=7.0D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C      PRINT*,'r IN DERHAM is',R(K,L)

      DQDR=1.0D0/((B2-A2)/2.0D0)

C  Ssssig

       IF (R(K,L).LT.6.501645139158383D0) THEN

       SSSs=G2*DQDR + G3*4.0D0*QTB(K,L)*DQDR
     /   +G4*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G5*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G6*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G7*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G8*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G9*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G10*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

       ELSE
 
       SSSs=0.0D0

       ENDIF


C   Sspsig

       IF (R(K,L).LT.6.514708000907687D0) THEN

      SSPs=G12*DQDR + G13*4.0D0*QTB(K,L)*DQDR
     /   +G14*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G15*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G16*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G17*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G18*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G19*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G20*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

       ELSE
 
      SSPs=0.0D0

       ENDIF


C swapped these two, due to fuckup in PRB

C Sppsig

       IF (R(K,L).LT.6.699591659971658D0) THEN

      SPPs=G22*DQDR + G23*4.0D0*QTB(K,L)*DQDR
     /   +G24*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G25*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G26*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G27*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G28*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G29*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G30*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

       ELSE

      SPPs=0.0D0

        ENDIF


C   Spppi

       IF (R(K,L).LT.6.841936337171249D0) THEN


      SPPp=G32*DQDR + G33*4.0D0*QTB(K,L)*DQDR
     /   +G34*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G35*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G36*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G37*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G38*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G39*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G40*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

         ELSE

            SPPp=0.0D0

         ENDIF

c         IF (R(K,L).GT.B2) THEN
c          SSSs=0.0D0
c          SSPs=0.0D0
c          SPPs=0.0D0
c          SPPp=0.0D0
c         ENDIF


      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   


      SUBROUTINE DAMCCB(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB,natoms)
      IMPLICIT NONE
      INTEGER  K,L,N, natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION SSSs,SSPs,SPPs,SPPp,DQDR,
     /     G1,G2,G3,XS(3*NATOMS),
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

       PARAMETER (  G1= -0.3789880000 ,G2=0.287963,
C old param; G2=  0.2886558000 ,
     /  G3= -0.1201660000 ,G4=  0.0100538000
     / ,G5=  0.0203983000 ,G6= -0.0149822000 ,G7=  0.0056378700 ,G8= -0.0006071250
     / ,G9= -0.0008101620,G10=  0.0006987440,
     /  G11=  0.2956740000 ,G12= -0.2010048000 ,G13=  0.0327500000 ,G14=  0.0514607000
     / ,G15= -0.0508434000 ,G16=  0.0288404000 ,G17= -0.0133668000 ,G18=  0.0046855100
     / ,G19=  0.0006011540,G20= -0.0000060000,
     /  G21=  0.2314380000 ,G22= -0.1469886000 ,G23=  0.0151875000 ,G24=  0.0440803000
     / ,G25= -0.0458782000 ,G26=  0.0312287000 ,G27= -0.0171949000 ,G28=  0.0079412400
     / ,G29= -0.0028890500,G30=  0.0012923000,
     /  G31= -0.3325740000 ,G32=  0.2847522000 ,G33= -0.1816530000 ,G34=  0.0901204000
     / ,G35= -0.0360385000 ,G36=  0.0113274000 ,G37= -0.0021334300 ,G38= -0.0006478810
     / ,G39=  0.0009400840,G40= -0.0009030410)

       A2=1.0D0
       B2=7.0D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

C      PRINT*,'r IN DERHAM is',R(K,L)

      DQDR=1.0D0/((B2-A2)/2.0D0)

C  Ssssig

       IF (R(K,L).LT.6.695259436126274D0) THEN

       SSSs=G2*DQDR + G3*4.0D0*QTB(K,L)*DQDR
     /   +G4*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G5*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G6*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G7*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G8*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G9*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G10*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)


      ELSE

      SSSs=0.0D0

      ENDIF

C   Sspsig

       IF (R(K,L).LT.6.796735011037814D0) THEN


      SSPs=G12*DQDR + G13*4.0D0*QTB(K,L)*DQDR
     /   +G14*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G15*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G16*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G17*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G18*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G19*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G20*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

         ELSE

        SSPs=0.0D0

         ENDIF



C swapped these two, due to fuckup in PRB

C Sppsig

       IF (R(K,L).LT.6.880535598840217D0) THEN

      SPPs=G22*DQDR + G23*4.0D0*QTB(K,L)*DQDR
     /   +G24*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G25*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G26*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G27*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G28*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G29*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G30*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

         ELSE

         SPPs=0.0D0

         ENDIF



C   Spppi

       IF (R(K,L).LT.5.944816716454667D0) THEN


      SPPp=G32*DQDR + G33*4.0D0*QTB(K,L)*DQDR
     /   +G34*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + G35*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +G36*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +G37*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +G38*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +G39*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     /   +G40*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)

        ELSE

          SPPp=0.0D0

        ENDIF


c         IF (R(K,L).GT.B2) THEN
c          SSSs=0.0D0
c          SSPs=0.0D0
c          SPPs=0.0D0
c          SPPp=0.0D0
c         ENDIF

      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE DREPCC(N,XS,K,L,rep,NATOMS,R, DIRCOS, DIAGA, DIAGB)
      IMPLICIT NONE
      INTEGER  K,L,N,natoms
      DOUBLE PRECISION R(NATOMS,NATOMS),
     1                 DIRCOS(NATOMS,NATOMS,3),
     3                 QTB(NATOMS,NATOMS),
     6                 DIAGA(4*NATOMS),DIAGB(4*NATOMS)


      DOUBLE PRECISION rep,
     /     C1,C2,C3,XS(3*NATOMS),
     7     C4,C5,C6,C7,C8,C9,C10,
     9     A2,B2,DQDR

       PARAMETER (   C1=  2.5911600000 ,C2= -2.0924020000 ,C3=  1.1228900000 ,C4= -0.4302620000
     / ,C5=  0.1257390000 ,C6= -0.0215445000 ,C7=  0.0000000000 ,C8=  0.0000000000
     / ,C9=  0.0000000000,C10=  0.0000000000)

       A2=1.0D0
       B2=3.3D0
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)

      DQDR=1.0D0/((B2-A2)/2.0D0)


C  Ssssig
      rep=C2*DQDR + C3*4.0D0*QTB(K,L)*DQDR
     /   +C4*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + C5*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +C6*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +C7*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +C8*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +C9*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     / +C10*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)


c         IF (R(K,L).GT.B2) rep=0.0D0

         IF (R(K,L).GT.3.251838574097556D0) THEN
            rep=0.0D0
         ELSE
c           rep=rep+0.00002383611355998103D0
         ENDIF


C       PRINT*,'in drepcc, rep is',rep
      RETURN
      END
