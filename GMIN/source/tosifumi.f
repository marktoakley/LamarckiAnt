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

C  Energy and derivatives of a general Born-Mayer potential with
C  Tosi-Fumi parameters specified after the keyword TOSI in the
C  odata file. Energy is in hartree, length in Bohr.
C
      SUBROUTINE TOSIFUMI(X, V, POTEL, GTEST, STEST)
      USE MODHESS
      USE commons
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2, K, L, J3, J4
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), POTEL, 
     1                 DUM, AC(NATOMS,NATOMS), Q(NATOMS), RD(NATOMS,NATOMS), RRD(NATOMS,NATOMS)
      LOGICAL GTEST, STEST

      N=NATOMS
      DO I=1,N
         IF (ZSYM(I).EQ.'MI') THEN
            Q(I)=-1.0D0
         ELSE IF (ZSYM(I).EQ.'PL') THEN
            Q(I)=1.0D0
         ELSE
            WRITE(*,'(A)') ' All atoms must be type PL or MI for Tosi-Fumi'
            STOP
         ENDIF
      ENDDO

       DO I=1,N
          DO J=I,N
             IF (ZSYM(I).EQ.'PL') THEN 
                IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APP
                   AC(J,I)=APP
                ELSE IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ELSE IF (ZSYM(I).EQ.'MI') THEN
                IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=AMM
                   AC(J,I)=AMM
                ELSE IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ENDIF 
         ENDDO
      ENDDO
C
C  Calculate the interparticle distances.
C 
      DO I=1,N
         K=1+3*(I-1)         
         RD(I,I)=0.0D0
         RRD(I,I)=0.0D0
         DO J=I+1,N
            L=K+3*(J-I)
            RD(I,J)= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                     (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                     (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RRD(I,J)=1.0D0/RD(I,J)
            RD(J,I)=RD(I,J)
            RRD(J,I)=RRD(I,J)
         ENDDO
      ENDDO

      POTEL=0.0D0
      DO I=1,N
         DO J=I+1,N
            POTEL=POTEL+Q(I)*Q(J)*RRD(J,I)+AC(J,I)*DEXP(-RD(J,I)/RHO)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  Gradient.
C
      DO J1=1,N
         DO J2=1,3
            DUM=0.0D0
            J4=3*(J1-1)+J2
            DO J3=1,N
               DUM=DUM-(Q(J3)*Q(J1)*RRD(J3,J1)**2+AC(J3,J1)*DEXP(-RD(J3,J1)/RHO)/RHO)*RRD(J3,J1)*(X(J4)-X(3*(J3-1)+J2))
            ENDDO
            V(J4)=DUM
         ENDDO
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  Hessian.
C
      DO I=1,3*N
         DO J=1,3*N
            HESS(I,J)=0.0D0
         ENDDO
      ENDDO
      DO I=1,3*N
         DO J=I,3*N
C 
C  Determine the cartesian coordinate of X(I) and X(J).
C
            K=MOD((I-1),3)
            L=MOD((J-1),3)
C 
C  This IF statement calls SHESS if I and J are the same
C  Cartesian coordinate and calls DHESS if the are different
C  Cartesian coordinates
C
            IF(K .EQ. L)THEN
               CALL NSHESS(I,J,K,L,X,N,AC,Q,RD)
            ELSE
               CALL NDHESS(I,J,K,L,X,N,AC,Q,RD)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C....................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE THE SAME
 
      SUBROUTINE NSHESS(I,J,K,L,CFG,NUM,AC,Q,RD)
      USE MODHESS
      USE commons
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM
      DOUBLE PRECISION CFG(3*NATOMS), RW, R, QQ, AB,
     1                 EP, DR, DF,AC(NATOMS,NATOMS),Q(NATOMS),RD(NATOMS,NATOMS)

      RW=RHO
 
c.....I1 and J1 are the particle numbers for X(I) and X(J)

      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
c.....if the particles are the same do the next two loops
c.....there are two loops so that we can skip particle I1
c.....without resorting to an IF statement

      IF(I1 .EQ. J1)THEN
      
      DO M = 1, I1-1
        
         N = 3*(M-1) + K + 1
 
         R = RD(I1,M)
         QQ = Q(I1)*Q(M)
         AB = AC(I1,M)
         EP = DEXP(-R/RW)
         DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
         DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1        + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2        - (AB*EP/(RW*R))
 
         HESS(I,I) = HESS(I,I) + DF
      
      ENDDO
 
      DO 20 M = I1+1, NUM
         
      N = 3*(M-1) + K + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
      DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1     + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2     - (AB*EP/(RW*R))
 
      HESS(I,I) = HESS(I,I) + DF
      
   20 CONTINUE
 
      ELSE
 
c.....do this if the particles are not the same
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(J))*(CFG(I)-CFG(J))
      
      HESS(I,J) = (QQ/(R*R*R)) - (3.0D0*QQ*DR/(R*R*R*R*R))
     1     - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
     2     + (AB*EP/(RW*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
C..........................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE NOT THE SAME
 
      SUBROUTINE NDHESS(I,J,K,L,CFG,NUM,AC,Q,RD)
      USE MODHESS
      USE commons
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM,N1
      DOUBLE PRECISION CFG(3*NATOMS),RW,R,QQ,AB,
     1                 EP,DR,DF,AC(NATOMS,NATOMS),Q(NATOMS),RD(NATOMS,NATOMS)
      RW=RHO
      
      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
c.....do the next to loops if the particles are the same

      IF(I1 .EQ. J1)THEN
 
      DO 10 M = 1, I1-1
 
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   10 CONTINUE
 
      DO 20 M = I1+1 ,NUM
         
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   20 CONTINUE
      HESS(J,I) = HESS(I,J)
 
      ELSE
 
c.....do this if the particles are not the same
      N = K + 3*(J1-1) + 1
      N1 = L + 3*(I1-1) + 1
 
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(N1)-CFG(J))
      
      HESS(I,J) = -(3.0D0*QQ*DR/(R*R*R*R*R))
     1         - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
