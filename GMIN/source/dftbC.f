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
C
C  Tiffany s TB routine for C.
C
      SUBROUTINE DFTBC(N,XS,deriv1st,ENERGY,GTEST)
      USE COMMONS, ONLY : MYUNIT
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER I, J, K, AI, AJ, NDIM,NMAX, I2, J2, K2, LWORK, INFO, J1, N, NDIAG
      DOUBLE PRECISION DIST,WORK(320*N),DPRAND,S(4*N,4*N),Sorig(4*N,4*N),
     &                 ENERGY,XS(3*N),HB(4*N,4*N),HA(4*N,4*N),VA(4*N,4*N),VB(4*N,4*N)

      DOUBLE PRECISION R(N,N), ABSTOL,DIRCOS(N,N,3),DIAG(4*N),deriv1st(3*N),DIAGA(4*N),DIAGB(4*N)

      DOUBLE PRECISION Ssssig,Sspsig,Sppsig,Spppi,EPSsac,EPSsbc,EPSpac,EPSpbc,Hsssiga,Hspsiga,
     &        Hppsiga,Hpppia,Hsssigb,Hspsigb,Hppsigb,Hpppib

      PARAMETER(EPSsac=-0.519289D0, EPSpac=-0.216472D0, EPSsbc=-0.434382D0, EPSpbc=-0.216472D0)
   
      INTEGER IWORK(20*N), IFAIL(4*N) 

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB


      DOUBLE PRECISION UREP,REP
      LOGICAL FTEST
      COMMON /FAIL/ FTEST

      LWORK=320*N

C   Calculate distance matrix

      DO I=1,N
       I2=3*(I-1)
       R(I,I)=0.0D0
          DO J=I+1,N
            J2=3*(J-1)
            DIST=(XS(I2+1)-XS(J2+1))**2+(XS(I2+2)-XS(J2+2))**2+(XS(I2+3)-XS(J2+3))**2
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

              CALL AREPCCC(N,I,J,REP,R)

C        PRINT*,'rep in h is',REP

        UREP=REP+UREP

         REP=0.0D0

        ENDDO
      ENDDO

C      PRINT*,' NEW urep in h is',UREP

C   Only include overlap or interaction with different atoms here.
               
      DO I=1,N
         AI=4*(I-1)
         DO J=I+1,N
            AJ=4*(J-1)

C   Calculate the overlap matrix S(I,J), related to these parameters
C   by the Slater-Koster scheme

           CALL OVECCC(N,I,J,Ssssig, Sspsig, Sppsig, Spppi, R)

           S((AI+1),(AJ+2))=DIRCOS(I,J,1)*Sspsig
           S((AI+1),(AJ+3))=DIRCOS(I,J,2)*Sspsig
           S((AI+1),(AJ+4))=DIRCOS(I,J,3)*Sspsig

           S((AI+2),(AJ+1))=-DIRCOS(I,J,1)*Sspsig
           S((AI+3),(AJ+1))=-DIRCOS(I,J,2)*Sspsig
           S((AI+4),(AJ+1))=-DIRCOS(I,J,3)*Sspsig

           S((AI+1),(AJ+1))=Ssssig
 
           S((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Sppsig+(1.0D0-DIRCOS(I,J,1)**2)*Spppi
           S((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Sppsig-Spppi)
           S((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Sppsig-Spppi)

           S((AI+3),(AJ+2))=S((AI+2),(AJ+3))
           S((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Sppsig+(1.0D0-DIRCOS(I,J,2)**2)*Spppi
           S((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Sppsig-Spppi)

           S((AI+4),(AJ+2))=S((AI+2),(AJ+4))
           S((AI+4),(AJ+3))=S((AI+3),(AJ+4))
           S((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Sppsig+(1.0D0-DIRCOS(I,J,3)**2)*Spppi
            
C   Hamiltonian is found using the Slater-Koster scheme
C   as shown in their paper and Harrison s book p.481
C   which relates H to the parameter V using direction cosines.
C   LAMBDA changes for each I and J
      
C my change means constructing  a and b type interaction matrices
C as an uhf analogue

           CALL HAMCCAC(N,I,J,Hsssiga, Hspsiga, Hppsiga, Hpppia, R)
           CALL HAMCCBC(N,I,J,Hsssigb, Hspsigb, Hppsigb, Hpppib, R)

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

           HA((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Hppsiga+(1.0D0-DIRCOS(I,J,1)**2)*Hpppia
           HA((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Hppsiga-Hpppia)
           HA((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Hppsiga-Hpppia)

           HA((AI+3),(AJ+2))=HA((AI+2),(AJ+3))
           HA((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Hppsiga+(1.0D0-DIRCOS(I,J,2)**2)*Hpppia
           HA((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Hppsiga-Hpppia)

           HA((AI+4),(AJ+2))=HA((AI+2),(AJ+4))
           HA((AI+4),(AJ+3))=HA((AI+3),(AJ+4))
           HA((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Hppsiga+(1.0D0-DIRCOS(I,J,3)**2)*Hpppia

C now for the beta type...

           HB((AI+1),(AJ+1))=Hsssigb

           HB((AI+2),(AJ+2))=DIRCOS(I,J,1)**2*Hppsigb+(1.0D0-DIRCOS(I,J,1)**2)*Hpppib
           HB((AI+2),(AJ+3))=DIRCOS(I,J,1)*DIRCOS(I,J,2)*(Hppsigb-Hpppib)
           HB((AI+2),(AJ+4))=DIRCOS(I,J,1)*DIRCOS(I,J,3)*(Hppsigb-Hpppib)

           HB((AI+3),(AJ+2))=HB((AI+2),(AJ+3))
           HB((AI+3),(AJ+3))=DIRCOS(I,J,2)**2*Hppsigb+(1.0D0-DIRCOS(I,J,2)**2)*Hpppib
           HB((AI+3),(AJ+4))=DIRCOS(I,J,2)*DIRCOS(I,J,3)*(Hppsigb-Hpppib)

           HB((AI+4),(AJ+2))=HB((AI+2),(AJ+4))
           HB((AI+4),(AJ+3))=HB((AI+3),(AJ+4))
           HB((AI+4),(AJ+4))=DIRCOS(I,J,3)**2*Hppsigb+(1.0D0-DIRCOS(I,J,3)**2)*Hpppib

        ENDDO

C   Interaction/Overlap of different AOs on the same atom
         
        DO K=1,4
           DO K2=K+1,4
              S((AI+K),(AI+K2))=0.0D0
C             H((AI+K),(AI+K2))=0.0D0

              HA((AI+K),(AI+K2))=0.0D0
              HB((AI+K),(AI+K2))=0.0D0
           ENDDO
        ENDDO

C   Interaction of same AO on the same atom

         HA((AI+1),(AI+1))=EPSsac
         HA((AI+2),(AI+2))=EPSpac
         HA((AI+3),(AI+3))=EPSpac
         HA((AI+4),(AI+4))=EPSpac

         HB((AI+1),(AI+1))=EPSsbc
         HB((AI+2),(AI+2))=EPSpbc
         HB((AI+3),(AI+3))=EPSpbc
         HB((AI+4),(AI+4))=EPSpbc

      ENDDO

      DO I=1,4*N
         S(I,I)=1.0D0
      ENDDO

      DO I=1,4*N
        DO J=1,4*N
           Sorig(I,J)=S(I,J)
        ENDDO
      ENDDO

      NDIM=4*N
      NMAX=4*N

      FTEST=.FALSE.

      ABSTOL=-1.0D0
      DO J1=1,4*N
         DIAG(J1)=1.0D100
      ENDDO
      CALL DSYGVX(1, 'V', 'I', 'U', NDIM, HA, NMAX, S, NMAX, 0.0D0,
     1     0.0D0, 1, 120, ABSTOL, NDIAG, DIAG, VA, NMAX, WORK, LWORK, IWORK, IFAIL, INFO)
C
C     CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, LWORK, INFO )
C
C     IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGENSORT_VAL_ASC(DIAG,VA,NDIM,NMAX)
      IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGSRT(DIAG,VA,NDIM,NMAX)
 
      IF (INFO.NE.0) THEN
         FTEST=.TRUE.
         ENERGY=1.0D6
         WRITE(MYUNIT,*) 'alpha DFTB diagonalisation failed - INFO=',INFO
         RETURN
      ENDIF

      DO I=1,4*N
         DIAGA(I)=DIAG(I)
      ENDDO

      DO I=1,4*N
        DO J=1,4*N
           S(I,J)=Sorig(I,J)
        END DO
      END DO

      FTEST=.FALSE.

      DO J1=1,4*N
         DIAG(J1)=1.0D100
      ENDDO

      CALL DSYGVX(1, 'V', 'A', 'U', NDIM, HB, NMAX, S, NMAX, 0.0D0,
     1     0.0D0, 1, 120, ABSTOL, NDIAG, DIAG, VB, NMAX, WORK, LWORK, IWORK, IFAIL, INFO)

C     CALL DSYGV( 1, 'V', 'U', NDIM, H, NMAX, S, NMAX, DIAG, WORK, LWORK, INFO )
C     IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGENSORT_VAL_ASC(DIAG,VB,NDIM,NMAX)
      IF (DIAG(1).LT.DIAG(NDIM)) CALL EIGSRT(DIAG,VB,NDIM,NMAX)

      IF (INFO.NE.0) THEN
         FTEST=.TRUE.
         ENERGY=1.0D6
         WRITE(MYUNIT,*) 'DSYGV failed - INFO=',INFO
         RETURN
      ENDIF

      DO I=1,4*N
         DIAGB(I)=DIAG(I)
      ENDDO

      ENERGY=0.0D0

      DO J=4*N-119,4*N
         ENERGY=ENERGY+DIAGA(J)+DIAGB(J)
      ENDDO

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

C   Now call the all new originally crafted subroutine to calculate
C   the first derivatives.

      IF (GTEST) CALL DFDERIV1C(N,XS,deriv1st,VA,VB,R, DIRCOS, DIAGA, DIAGB)


      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   
      SUBROUTINE OVECCC(N,K,L,SSSs,SSPs,SPPs,SPPp,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,
     3                  D1,D2,D3,D4,D5,D6,A3,B3,
     4                 D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,
     5                 D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,
     6             D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40

      PARAMETER (   D1=  0.4635980000 ,D2=-0.354943D0,
C old param :D2= -0.3593880000 ,
     /  D3=  0.1597240000 ,D4= -0.0236036000
     / ,D5= -0.0160748000 ,D6=  0.0101939000 ,D7= -0.0010466800 ,D8= -0.0013805000
     / ,D9=  0.0007822820,D10= -0.0001806570,

     /  D11= -0.3627180000 ,D12=  0.2504088000 ,D13= -0.0459648000 ,D14= -0.0572412000
     / ,D15=  0.0499554000 ,D16= -0.0161908000 ,D17= -0.0007504350 ,D18=  0.0030270400
     / ,D19= -0.0013407800,D20=  0.0002213130,

     /  D21= -0.1381954000 , D22=0.021232D0,
C old param: D22=  0.0272058000 ,
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
      SUBROUTINE HAMCCAC(N,K,L,SSSs,SSPs,SPPs,SPPp,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,
     /     G1,G2,G3,
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

       PARAMETER (  G1= -0.4546540000 ,G2=  0.3501960000 ,G3= -0.1513000000 ,G4=  0.0159980000
     / ,G5=  0.0247918000 ,G6= -0.0192764000 ,G7=  0.0071180300 ,G8= -0.0002604590
     / ,G9= -0.0014905100,G10=  0.0013685600,

     /  G11=  0.3813740000 ,G12= -0.2675544000 ,G13=  0.0587541000 ,G14=  0.0551394000
     / ,G15= -0.0624710000 ,G16=  0.0366550000 ,G17= -0.0161735000 ,G18=  0.0047799300
     / ,G19= -0.0000637944,G20= -0.0010983800,

     /  G21=  0.2788480000 , G22=-0.1669651D0,
C old param G22= -0.1770288000 ,
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

      SUBROUTINE HAMCCBC(N,K,L,SSSs,SSPs,SPPs,SPPp,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,
     /     G1,G2,G3,
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

      SUBROUTINE AREPCCC(N,K,L,rep,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

      DOUBLE PRECISION REP,
     /     C1,C2,C3,
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

      SUBROUTINE DFDERIV1C(N,XS,deriv1st,EIGVECA,EIGVECB,R, DIRCOS, DIAGA, DIAGB)
      IMPLICIT NONE 
      INTEGER N
      DOUBLE PRECISION R(N,N),DERIV1ST(3*N),
     1                 DIRCOS(N,N,3),
     6                 DIAGA(4*N),DIAGB(4*N)

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB

      INTEGER I, J, K , AI, AJ, LEV, NUM1, NUM2, I2,J2
      DOUBLE PRECISION REP1ST(N,3),RR, D1, D2, D3,
     1                 DUMMY, DUMMY1,
     2                 diffSVc(3,4,4), SDC1, SQC1, SQC2, SQC3,
     3                 R1st(N,N), diffSVa(3,4,4), diffSVb(3,4,4),
C made urep first a scalar
     4                 UREP1st,XS(3*N),E1Aa,E1Ba,E1Ab,E1Bb,E1Ac,E1Bc,
     5                 ELEC1sta(N,3),ELEC1stb(N,3),EIGVECA(4*N,4*N),DUMVECA(4*N),DUMVECB(4*N),
     6                 EIGVECB(4*N,4*N)
      DOUBLE PRECISION SVSUB,TEMP1A,TEMP2A,TEMP1B,TEMP2B,
     1                 SQDIRCOS(3*N,N),SUBDIRCOS(3*N,N),
     2                 DELD1, DELD12, DELD2, DELD22, DELD3, DELD32, DELRDX, DSVSUB,
     3                 DRSssa,DRSspa,DRSppsa,DRSpppa,
     3                 DRSssb,DRSspb,DRSppsb,DRSpppb,
     4                 Ssssig, Sspsig, Sppsig, Spppi,
     5                 Hsssiga, Hspsiga, Hppsiga, Hpppia,
     5                 Hsssigb, Hspsigb, Hppsigb, Hpppib,
     5                 DORSss,DORSsp,DORSpps,DORSppp
      LOGICAL TEST1,TEST2,TEST3

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

C occupation numbers

C
C   Can put arrays into a block data format instead of parameter
C   so can just put in once in a file to be included

C   Can t start an array with a no. ie. 1stderiv...syntax error

C   First of all calculate derivatives of repulsive energy term
C   This is called REP1st.

C      WRITE(6,*)'First going to do derivative of repulsive energy:'

      DO K=1,3
         DO I=1,N
            DO J=I+1,N

               CALL DREPCCC(N,I,J,UREP1st,R)

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
            ENDDO
         ENDDO
      ENDDO

C     PRINT*,'Entering the differentiation of electronic energy zone'
        DO K=1,N
           ELEC1sta(K,1)=0.0D0
           ELEC1stb(K,1)=0.0D0
           ELEC1sta(K,2)=0.0D0
           ELEC1stb(K,2)=0.0D0
           ELEC1sta(K,3)=0.0D0
           ELEC1stb(K,3)=0.0D0
        ENDDO

        DO I=1,N
           AI=4*(I-1)
           I2=3*(I-1)


           DO J=I+1,N
              IF (R(J,I).GT.7.031215111280781D0) GOTO 666
              RR=1.0D0/R(J,I)
              AJ=4*(J-1)
              J2=3*(J-1)

              D1=DIRCOS(J,I,1)
              D2=DIRCOS(J,I,2)
              D3=DIRCOS(J,I,3)

              SDC1=SUBDIRCOS(J2+1,I)
              SQC1=SQDIRCOS(J2+1,I)
              SQC2=SQDIRCOS(J2+2,I)
              SQC3=SQDIRCOS(J2+3,I)

C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

              CALL OVECCC(N,I,J,Ssssig, Sspsig, Sppsig, Spppi,R)
              CALL DOVECCC(N,I,J,DORSss,DORSsp,DORSpps,DORSppp,R)

              DUMMY1=DORSsp*DELRDX

              diffSVa(1,1,2)=Sspsig*DELD1 - D1*DUMMY1
              diffSVa(1,1,3)=Sspsig*DELD2 - D2*DUMMY1
              diffSVa(1,1,4)=Sspsig*DELD3 - D3*DUMMY1

              diffSVa(1,2,1)=-diffSVa(1,1,2)
              diffSVa(1,3,1)=-diffSVa(1,1,3)
              diffSVa(1,4,1)=-diffSVa(1,1,4)

              SVSUB=Sppsig - Spppi
              DSVSUB=DORSpps - DORSppp

              diffSVa(1,1,1)=DORSss*DELRDX
              DUMMY1=DELRDX*DSVSUB

              diffSVa(1,2,2)=SQC1*DUMMY1  + SVSUB*DELD12   + DORSppp*DELRDX
              diffSVa(1,2,3)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSVa(1,2,4)=D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
                
              diffSVa(1,3,2)= diffSVa(1,2,3)
              diffSVa(1,3,3)= SQC2*DUMMY1 + SVSUB*DELD22 + DORSppp*DELRDX
              diffSVa(1,3,4)= D2*D3*DUMMY1 - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSVa(1,4,2)= diffSVa(1,2,4)
              diffSVa(1,4,3)= diffSVa(1,3,4)
              diffSVa(1,4,4)= SQC3*DUMMY1 + SVSUB*DELD32 + DORSppp*DELRDX

              CALL HAMCCAC(N,I,J,Hsssiga, Hspsiga, Hppsiga, Hpppia,R)
              CALL DAMCCAC(N,I,J,DRSssa,DRSspa,DRSppsa,DRSpppa,R)

              DUMMY1=DRSspa*DELRDX

              diffSVa(2,1,2)=Hspsiga*DELD1 - D1*DUMMY1
              diffSVa(2,1,3)=Hspsiga*DELD2 - D2*DUMMY1
              diffSVa(2,1,4)=Hspsiga*DELD3 - D3*DUMMY1

              diffSVa(2,2,1)=-diffSVa(2,1,2)
              diffSVa(2,3,1)=-diffSVa(2,1,3)
              diffSVa(2,4,1)=-diffSVa(2,1,4)

              SVSUB=Hppsiga - Hpppia
              DSVSUB=DRSppsa - DRSpppa

              diffSVa(2,1,1)=DRSssa*DELRDX
              DUMMY1=DELRDX*DSVSUB

              diffSVa(2,2,2)=SQC1*DUMMY1  + SVSUB*DELD12   + DRSpppa*DELRDX
              diffSVa(2,2,3)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSVa(2,2,4)=D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSVa(2,3,2)= diffSVa(2,2,3)
              diffSVa(2,3,3)= SQC2*DUMMY1  + SVSUB*DELD22    + DRSpppa*DELRDX
              diffSVa(2,3,4)= D2*D3*DUMMY1 - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSVa(2,4,2)= diffSVa(2,2,4)
              diffSVa(2,4,3)= diffSVa(2,3,4)
              diffSVa(2,4,4)= SQC3*DUMMY1 + SVSUB*DELD32 + DRSpppa*DELRDX

              CALL HAMCCBC(N,I,J,Hsssigb, Hspsigb, Hppsigb, Hpppib,R)
              CALL DAMCCBC(N,XS,I,J,DRSssb,DRSspb,DRSppsb,DRSpppb,R, DIRCOS, DIAGA, DIAGB)

              DUMMY1=DRSspb*DELRDX
      
              diffSVa(3,1,2)=Hspsigb*DELD1 - D1*DUMMY1
              diffSVa(3,1,3)=Hspsigb*DELD2 - D2*DUMMY1
              diffSVa(3,1,4)=Hspsigb*DELD3 - D3*DUMMY1
     
              diffSVa(3,2,1)=-diffSVa(3,1,2)
              diffSVa(3,3,1)=-diffSVa(3,1,3)
              diffSVa(3,4,1)=-diffSVa(3,1,4)

              SVSUB=Hppsigb - Hpppib
              DSVSUB=DRSppsb - DRSpppb

              diffSVa(3,1,1)=DRSssb*DELRDX
              DUMMY1=DELRDX*DSVSUB

              diffSVa(3,2,2)=SQC1*DUMMY1 + SVSUB*DELD12 + DRSpppb*DELRDX
              diffSVa(3,2,3)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1
              diffSVa(3,2,4)=D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSVa(3,3,2)= diffSVa(3,2,3)
              diffSVa(3,3,3)= SQC2*DUMMY1 + SVSUB*DELD22 + DRSpppb*DELRDX
              diffSVa(3,3,4)= D2*D3*DUMMY1 - SVSUB*D2*DELD3  - SVSUB*D3*DELD2

              diffSVa(3,4,2)= diffSVa(3,2,4)
              diffSVa(3,4,3)= diffSVa(3,3,4)
              diffSVa(3,4,4)= SQC3*DUMMY1 + SVSUB*DELD32 + DRSpppb*DELRDX

Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
C   Differentiation of overlap or interaction between orbitals
C   on different atoms

              D1=DIRCOS(J,I,2)
              D2=DIRCOS(J,I,3)
              D3=DIRCOS(J,I,1)

              SDC1=SUBDIRCOS(J2+2,I)
              SQC1=SQDIRCOS(J2+2,I)
              SQC2=SQDIRCOS(J2+3,I)
              SQC3=SQDIRCOS(J2+1,I)

              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

              diffSVb(1,1,3)=(Sspsig*DELD1 - D1*DORSsp*DELRDX)
              diffSVb(1,1,4)=(Sspsig*DELD2 - D2*DORSsp*DELRDX)
              diffSVb(1,1,2)=diffSVa(1,1,3)

              diffSVb(1,3,1)=-diffSVb(1,1,3)
              diffSVb(1,4,1)=-diffSVb(1,1,4)
              diffSVb(1,2,1)=-diffSVb(1,1,2)

              SVSUB=Sppsig - Spppi
              DSVSUB=DORSpps - DORSppp

              diffSVb(1,1,1)=DORSss*DELRDX

              diffSVb(1,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DORSppp*DELRDX
              diffSVb(1,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSVb(1,2,4)=diffSVa(1,3,4)

              diffSVb(1,3,2)= diffSVb(1,2,3)
              diffSVb(1,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DORSppp*DELRDX
              diffSVb(1,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSVb(1,4,2)= diffSVb(1,2,4)
              diffSVb(1,4,3)= diffSVb(1,3,4)
              diffSVb(1,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DORSppp*DELRDX

              diffSVb(2,1,3)=(Hspsiga*DELD1 - D1*DRSspa*DELRDX)
              diffSVb(2,1,4)=(Hspsiga*DELD2 - D2*DRSspa*DELRDX)
              diffSVb(2,1,2)=diffSVa(2,1,3)

              diffSVb(2,3,1)=-diffSVb(2,1,3)
              diffSVb(2,4,1)=-diffSVb(2,1,4)
              diffSVb(2,2,1)=-diffSVb(2,1,2)

              SVSUB=Hppsiga - Hpppia
              DSVSUB=DRSppsa - DRSpppa

              diffSVb(2,1,1)=DRSssa*DELRDX

              diffSVb(2,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSpppa*DELRDX
              diffSVb(2,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSVb(2,2,4)=diffSVa(2,3,4)

              diffSVb(2,3,2)= diffSVb(2,2,3)
              diffSVb(2,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSpppa*DELRDX
              diffSVb(2,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSVb(2,4,2)= diffSVb(2,2,4)
              diffSVb(2,4,3)= diffSVb(2,3,4)
              diffSVb(2,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSpppa*DELRDX

              diffSVb(3,1,3)=(Hspsigb*DELD1 - D1*DRSspb*DELRDX)
              diffSVb(3,1,4)=(Hspsigb*DELD2 - D2*DRSspb*DELRDX)
              diffSVb(3,1,2)=diffSVa(3,1,3)

              diffSVb(3,3,1)=-diffSVb(3,1,3)
              diffSVb(3,4,1)=-diffSVb(3,1,4)
              diffSVb(3,2,1)=-diffSVb(3,1,2)

              SVSUB=Hppsigb - Hpppib
              DSVSUB=DRSppsb - DRSpppb

              diffSVb(3,1,1)=DRSssb*DELRDX

              diffSVb(3,3,3)=SQC1*DELRDX*DSVSUB + SVSUB*DELD12 + DRSpppb*DELRDX
              diffSVb(3,2,3)=D1*D3*DELRDX*DSVSUB - SVSUB*D1*DELD3 - SVSUB*D3*DELD1
              diffSVb(3,2,4)=diffSVa(3,3,4)

              diffSVb(3,3,2)= diffSVb(3,2,3)
              diffSVb(3,4,4)= SQC2*DELRDX*DSVSUB + SVSUB*DELD22 + DRSpppb*DELRDX
              diffSVb(3,3,4)= D2*D1*DELRDX*DSVSUB - SVSUB*D2*DELD1 - SVSUB*D1*DELD2

              diffSVb(3,4,2)= diffSVb(3,2,4)
              diffSVb(3,4,3)= diffSVb(3,3,4)
              diffSVb(3,2,2)= SQC3*DELRDX*DSVSUB + SVSUB*DELD32 + DRSpppb*DELRDX

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C   diffSV(diff. of S or V,AO no.,AO no.)

C   Differentiation of overlap or interaction between orbitals
C   on different atoms

              D1=DIRCOS(J,I,3)
              D2=DIRCOS(J,I,1)
              D3=DIRCOS(J,I,2)

              SDC1=SUBDIRCOS(J2+3,I)
              SQC1=SQDIRCOS(J2+3,I)
              SQC2=SQDIRCOS(J2+1,I)
              SQC3=SQDIRCOS(J2+2,I)

              DELD1=RR*SDC1
              DELD12=-2*RR*D1*SDC1
              DELD2=D2*D1*RR
              DELD22=-2*SQC2*D1*RR
              DELD3=D3*D1*RR
              DELD32=-2*SQC3*D1*RR
              DELRDX=D1

              DUMMY1=DORSsp*DELRDX

              diffSVc(1,1,4)=Sspsig*DELD1 - D1*DUMMY1
              diffSVc(1,1,2)=diffSVa(1,1,4)
              diffSVc(1,1,3)=diffSVb(1,1,4)

              diffSVc(1,4,1)=-diffSVc(1,1,4)
              diffSVc(1,2,1)=-diffSVc(1,1,2)
              diffSVc(1,3,1)=-diffSVc(1,1,3)

              SVSUB=Sppsig - Spppi
              DSVSUB=DORSpps - DORSppp

              diffSVc(1,1,1)=DORSss*DELRDX

              DUMMY1=DELRDX*DSVSUB

              diffSVc(1,4,4)=SQC1*DUMMY1 + SVSUB*DELD12 + DORSppp*DELRDX
              diffSVc(1,2,3)=diffSVa(1,3,4)
              diffSVc(1,2,4)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSVc(1,3,2)= diffSVc(1,2,3)
              diffSVc(1,2,2)= SQC2*DUMMY1 + SVSUB*DELD22 + DORSppp*DELRDX
              diffSVc(1,3,4)= D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSVc(1,4,2)= diffSVc(1,2,4)
              diffSVc(1,4,3)= diffSVc(1,3,4)
              diffSVc(1,3,3)= SQC3*DUMMY1 + SVSUB*DELD32 + DORSppp*DELRDX

              DUMMY1=DRSspa*DELRDX

              diffSVc(2,1,4)=Hspsiga*DELD1 - D1*DUMMY1
              diffSVc(2,1,2)=diffSVa(2,1,4)
              diffSVc(2,1,3)=diffSVb(2,1,4)

              diffSVc(2,4,1)=-diffSVc(2,1,4)
              diffSVc(2,2,1)=-diffSVc(2,1,2)
              diffSVc(2,3,1)=-diffSVc(2,1,3)

              SVSUB=Hppsiga - Hpppia
              DSVSUB=DRSppsa - DRSpppa

              diffSVc(2,1,1)=DRSssa*DELRDX

              DUMMY1=DELRDX*DSVSUB

              diffSVc(2,4,4)=SQC1*DUMMY1 + SVSUB*DELD12 + DRSpppa*DELRDX
              diffSVc(2,2,3)=diffSVa(2,3,4)
              diffSVc(2,2,4)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSVc(2,3,2)= diffSVc(2,2,3)
              diffSVc(2,2,2)= SQC2*DUMMY1 + SVSUB*DELD22 + DRSpppa*DELRDX
              diffSVc(2,3,4)= D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSVc(2,4,2)= diffSVc(2,2,4)
              diffSVc(2,4,3)= diffSVc(2,3,4)
              diffSVc(2,3,3)= SQC3*DUMMY1 + SVSUB*DELD32 + DRSpppa*DELRDX

              DUMMY1=DRSspb*DELRDX

              diffSVc(3,1,4)=Hspsigb*DELD1 - D1*DUMMY1
              diffSVc(3,1,2)=diffSVa(3,1,4)
              diffSVc(3,1,3)=diffSVb(3,1,4)

              diffSVc(3,4,1)=-diffSVc(3,1,4)
              diffSVc(3,2,1)=-diffSVc(3,1,2)
              diffSVc(3,3,1)=-diffSVc(3,1,3)

              SVSUB=Hppsigb - Hpppib
              DSVSUB=DRSppsb - DRSpppb

              diffSVc(3,1,1)=DRSssb*DELRDX

              DUMMY1=DELRDX*DSVSUB

              diffSVc(3,4,4)=SQC1*DUMMY1 + SVSUB*DELD12 + DRSpppb*DELRDX
              diffSVc(3,2,3)=diffSVa(3,3,4)
              diffSVc(3,2,4)=D1*D2*DUMMY1 - SVSUB*D1*DELD2 - SVSUB*D2*DELD1

              diffSVc(3,3,2)= diffSVc(3,2,3)
              diffSVc(3,2,2)= SQC2*DUMMY1 + SVSUB*DELD22 + DRSpppb*DELRDX
              diffSVc(3,3,4)= D1*D3*DUMMY1 - SVSUB*D1*DELD3 - SVSUB*D3*DELD1

              diffSVc(3,4,2)= diffSVc(3,2,4)
              diffSVc(3,4,3)= diffSVc(3,3,4)
              diffSVc(3,3,3)= SQC3*DUMMY1 + SVSUB*DELD32 + DRSpppb*DELRDX

              E1Aa=0.0D0
              E1Ba=0.0D0
              E1Ab=0.0D0
              E1Bb=0.0D0
              E1Ac=0.0D0
              E1Bc=0.0D0

              DO NUM1=1,4
                 DO LEV=4*N-119,4*N
                    DUMVECA(LEV)=EIGVECA(LEV,AI+NUM1)
                    DUMVECB(LEV)=EIGVECB(LEV,AI+NUM1)
                 ENDDO
                 DO NUM2=1,4
 
                    TEST1=((diffSVa(1,NUM1,NUM2).EQ.0.0D0).AND.
     1                     (diffSVb(1,NUM1,NUM2).EQ.0.0D0).AND.
     2                     (diffSVc(1,NUM1,NUM2).EQ.0.0D0))
                    TEST2=((diffSVa(2,NUM1,NUM2).EQ.0.0D0).AND.
     1                     (diffSVb(2,NUM1,NUM2).EQ.0.0D0).AND.
     2                     (diffSVc(2,NUM1,NUM2).EQ.0.0D0))
                    TEST3=((diffSVa(3,NUM1,NUM2).EQ.0.0D0).AND.
     1                     (diffSVb(3,NUM1,NUM2).EQ.0.0D0).AND.
     2                     (diffSVc(3,NUM1,NUM2).EQ.0.0D0))
C                   WRITE(*,'(A,4I6,3L)') 'AI,AJ,NUM1,NUM2,T1,T2,T3=',AI,AJ,NUM1,NUM2,TEST1,TEST2,TEST3
                    IF (TEST1.AND.TEST2.AND.TEST3) GOTO 555
                    TEMP1A=0.0D0
                    TEMP2A=0.0D0
                    TEMP1B=0.0D0
                    TEMP2B=0.0D0
                    IF (TEST1) THEN
                       DO LEV=4*N-119,4*N
                          TEMP1A=TEMP1A+DUMVECA(LEV)*EIGVECA(LEV,AJ+NUM2)
                          TEMP1B=TEMP1B+DUMVECB(LEV)*EIGVECB(LEV,AJ+NUM2)
                       ENDDO !LEV
                    ELSE IF (TEST2.AND.TEST3) THEN
                       DO LEV=4*N-119,4*N
                          TEMP2A=TEMP2A-DUMVECA(LEV)*EIGVECA(LEV,AJ+NUM2)*DIAGA(LEV)
                          TEMP2B=TEMP2B-DUMVECB(LEV)*EIGVECB(LEV,AJ+NUM2)*DIAGB(LEV)
                       ENDDO !LEV
                    ELSE IF (TEST2) THEN
                       DO LEV=4*N-119,4*N
                          TEMP2A=TEMP2A-DUMVECA(LEV)*EIGVECA(LEV,AJ+NUM2)*DIAGA(LEV)
                          DUMMY=DUMVECB(LEV)*EIGVECB(LEV,AJ+NUM2)
                          TEMP1B=TEMP1B+DUMMY
                          TEMP2B=TEMP2B-DUMMY*DIAGB(LEV)
                       ENDDO !LEV
                    ELSE IF (TEST3) THEN
                       DO LEV=4*N-119,4*N
                          DUMMY=DUMVECA(LEV)*EIGVECA(LEV,AJ+NUM2)
                          TEMP1A=TEMP1A+DUMMY
                          TEMP2A=TEMP2A-DUMMY*DIAGA(LEV)
                          TEMP2B=TEMP2B-DUMVECB(LEV)*EIGVECB(LEV,AJ+NUM2)*DIAGB(LEV)
                       ENDDO !LEV
                    ELSE
                       DO LEV=4*N-119,4*N
                          DUMMY=DUMVECA(LEV)*EIGVECA(LEV,AJ+NUM2)
                          TEMP1A=TEMP1A+DUMMY
                          TEMP2A=TEMP2A-DUMMY*DIAGA(LEV)
                          DUMMY=DUMVECB(LEV)*EIGVECB(LEV,AJ+NUM2)
                          TEMP1B=TEMP1B+DUMMY
                          TEMP2B=TEMP2B-DUMMY*DIAGB(LEV)
                       ENDDO !LEV
                    ENDIF
                    E1Aa=E1Aa+TEMP1A*diffSVa(2,NUM1,NUM2)+TEMP2A*diffSVa(1,NUM1,NUM2)
                    E1Ab=E1Ab+TEMP1A*diffSVb(2,NUM1,NUM2)+TEMP2A*diffSVb(1,NUM1,NUM2)
                    E1Ac=E1Ac+TEMP1A*diffSVc(2,NUM1,NUM2)+TEMP2A*diffSVc(1,NUM1,NUM2)
                    E1Ba=E1Ba+TEMP1B*diffSVa(3,NUM1,NUM2)+TEMP2B*diffSVa(1,NUM1,NUM2)
                    E1Bb=E1Bb+TEMP1B*diffSVb(3,NUM1,NUM2)+TEMP2B*diffSVb(1,NUM1,NUM2)
                    E1Bc=E1Bc+TEMP1B*diffSVc(3,NUM1,NUM2)+TEMP2B*diffSVc(1,NUM1,NUM2)
555                 CONTINUE 
                 ENDDO !NUM2
              ENDDO !NUM1
 
              ELEC1sta(I,1)=ELEC1sta(I,1)+2.0D0*E1Aa
              ELEC1sta(J,1)=ELEC1sta(J,1)-2.0D0*E1Aa
              ELEC1stb(I,1)=ELEC1stb(I,1)+2.0D0*E1Ba
              ELEC1stb(J,1)=ELEC1stb(J,1)-2.0D0*E1Ba
              ELEC1sta(I,2)=ELEC1sta(I,2)+2.0D0*E1Ab
              ELEC1sta(J,2)=ELEC1sta(J,2)-2.0D0*E1Ab
              ELEC1stb(I,2)=ELEC1stb(I,2)+2.0D0*E1Bb
              ELEC1stb(J,2)=ELEC1stb(J,2)-2.0D0*E1Bb
              ELEC1sta(I,3)=ELEC1sta(I,3)+2.0D0*E1Ac
              ELEC1sta(J,3)=ELEC1sta(J,3)-2.0D0*E1Ac
              ELEC1stb(I,3)=ELEC1stb(I,3)+2.0D0*E1Bc
              ELEC1stb(J,3)=ELEC1stb(J,3)-2.0D0*E1Bc

666           CONTINUE

           ENDDO  !J 
        ENDDO  !

C   In order to differentiate by the coordinate of atom K, you only get 
C   a value for the derivative if one of the two atoms involved in the
C   pair interaction is K otherwise it is zero. This is achieved by use
C   of the FACTOR term. I m sure this is highly inefficient BUT if it
C   works then that is worthwhile. I can improve it later....when I m 
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
            deriv1st(I2+K)=ELEC1sta(I,K)+ELEC1stb(I,K)+REP1st(I,K)
         ENDDO
      ENDDO
       
      RETURN
      END
C   
C  subroutine to expand the polynomial sum for the overlap matrix elements
C  this is between same element....same element/diff element will be more
C   

      SUBROUTINE DOVECCC(N,K,L,SSSs,SSPs,SPPs,SPPp,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,
     3                  D1,D2,D3,D4,D5,D6,A3,B3,
     4                 D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,
     5                 D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,
     6             D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40,DQDR

      PARAMETER ( D1=  0.4635980000 , D2=-0.354943D0,
C old params D2= -0.3593880000 ,
     /  D3=  0.1597240000 ,D4= -0.0236036000
     / ,D5= -0.0160748000 ,D6=  0.0101939000 ,D7= -0.0010466800 ,D8= -0.0013805000
     / ,D9=  0.0007822820,D10= -0.0001806570,

     /  D11= -0.3627180000 ,D12=  0.2504088000 ,D13= -0.0459648000 ,D14= -0.0572412000
     / ,D15=  0.0499554000 ,D16= -0.0161908000 ,D17= -0.0007504350 ,D18=  0.0030270400
     / ,D19= -0.0013407800,D20=  0.0002213130,

     /  D21= -0.1381954000 , D22=0.021232D0,
C old params D22=  0.0272058000 ,
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

      SUBROUTINE DAMCCAC(N,K,L,SSSs,SSPs,SPPs,SPPp,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,DQDR,
     /     G1,G2,G3,
     7     G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,
     8     G21,G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G33,G34,G35,G36,
     9     G37,G38,G39,G40,A2,B2

       PARAMETER ( G1= -0.4546540000 ,G2=  0.3501960000 ,G3= -0.1513000000 ,G4=  0.0159980000
     / ,G5=  0.0247918000 ,G6= -0.0192764000 ,G7=  0.0071180300 ,G8= -0.0002604590
     / ,G9= -0.0014905100,G10=  0.0013685600,

     /  G11=  0.3813740000 ,G12= -0.2675544000 ,G13=  0.0587541000 ,G14=  0.0551394000
     / ,G15= -0.0624710000 ,G16=  0.0366550000 ,G17= -0.0161735000 ,G18=  0.0047799300
     / ,G19= -0.0000637944,G20= -0.0010983800,

     /  G21=  0.2788480000 ,G22=-0.1669651D0,
C old param: G22= -0.1770288000 ,
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


      SUBROUTINE DAMCCBC(N,XS,K,L,SSSs,SSPs,SPPs,SPPp,R, DIRCOS, DIAGA, DIAGB)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N),
     1                 DIRCOS(N,N,3),
     3                 QTB(N,N),
     6                 DIAGA(4*N),DIAGB(4*N)

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB

      DOUBLE PRECISION SSSS,SSPS,SPPS,SPPP,DQDR,
     /     G1,G2,G3,XS(3*N),
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

      SUBROUTINE DREPCCC(N,K,L,rep,R)
      IMPLICIT NONE
      INTEGER  K,L,N
      DOUBLE PRECISION R(N,N), QTB(N,N)

C     COMMON /RDIST/ R, DIRCOS, DIAGA, DIAGB

      DOUBLE PRECISION REP,
     /     C1,C2,C3,
     7     C4,C5,C6,C7,C8,C9,C10,
     9     A2,B2,DQDR

       PARAMETER (   C1=  2.5911600000 ,C2= -2.0924020000 ,C3=  1.1228900000 ,C4= -0.4302620000
     / ,C5=  0.1257390000 ,C6= -0.0215445000 ,C7=  0.0000000000 ,C8=  0.0000000000
     / ,C9=  0.0000000000,C10=  0.0000000000, A2=1.0D0, B2=3.3D0)

C      A2=1.0D0
C      B2=3.3D0

C         IF (R(K,L).GT.B2) rep=0.0D0

         IF (R(K,L).GT.3.251838574097556D0) THEN
            rep=0.0D0
         ELSE
C  Ssssig
       QTB(K,L)=(R(K,L)-((B2+A2)/2.0D0))/((B2-A2)/2.0D0)
       DQDR=1.0D0/((B2-A2)/2.0D0)
      rep=C2*DQDR + C3*4.0D0*QTB(K,L)*DQDR
     /   +C4*(-3.0D0*DQDR + 12.0D0*(QTB(K,L)**2)*DQDR)
     /   + C5*(-16.0D0*QTB(K,L)*DQDR + 32.0D0*(QTB(K,L)**3)*DQDR) +C6*(5.0D0*DQDR -60.0D0*(QTB(K,L)**2)*DQDR
     /   + 80.0D0*(QTB(K,L)**4)*DQDR) +C7*(36.0D0*QTB(K,L)*DQDR - 192.0D0*(QTB(K,L)**3)*DQDR +192.0D0*(QTB(K,L)**5)*DQDR)
     /   +C8*(-7.0D0*DQDR +168.0D0*(QTB(K,L)**2)*DQDR - 560.0D0*(QTB(K,L)**4)*DQDR +448.0D0*(QTB(K,L)**6)*DQDR)
     /   +C9*(-64.0D0*QTB(K,L)*DQDR +640.0D0*(QTB(K,L)**3)*DQDR -1536.0D0*(QTB(K,L)**5)*DQDR +1024.0D0*(QTB(K,L)**7)*DQDR)
     / +C10*(9.0D0*DQDR -360.0D0*(QTB(K,L)**2)*DQDR +2160*(QTB(K,L)**4)*DQDR
     /   -4032.0D0*(QTB(K,L)**6)*DQDR +2304.0D0*(QTB(K,L)**8)*DQDR)
C           rep=rep+0.00002383611355998103D0
         ENDIF


C       PRINT*,'in drepcc, rep is',rep
      RETURN
      END
