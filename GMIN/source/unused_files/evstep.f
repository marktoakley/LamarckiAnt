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
      SUBROUTINE EVSTEP(NP,VECTORS)
      USE commons
      IMPLICIT NONE
      
      DOUBLE PRECISION ENERGY,
     1                 EVALMIN, EVALMAX, X(3*MXATMS), GRAD(3*MXATMS),
     2                 VECTORS(3*MXATMS,20), DUMMY
      INTEGER J1, NP, J2
C
C  Save coordinates and atomic energies in case the step is reversed.
C
      DO J1=1,3*(NATOMS-NSEED)
         COORDSO(J1,NP)=COORDS(J1,NP)
         X(J1)=COORDS(J1,NP)
      ENDDO
      DO J1=1,NATOMS
         VATO(J1,NP)=VAT(J1,NP)
      ENDDO

      CALL POTENTIAL(X,GRAD,ENERGY,.FALSE.,.TRUE.)
      CALL MYITEIG(NP,VECTORS,EVALMIN,EVALMAX)
C
C  Renormalize to give the same likely normalization as for vectors with components
C  randomly sampled from the uniform distribution. Expected sum of squares is
C  3*NATOMS*2/3=2*NATOMS. Hence we need to multiply by SQRT(2*NATOMS).
C
      DUMMY=DSQRT(2.0D0*NATOMS)
      DO J1=1,NVECTORS
         DO J2=1,3*NATOMS
            VECTORS(J2,J1)=VECTORS(J2,J1)*DUMMY
         ENDDO
      ENDDO
   
      RETURN 
      END
C
C***************************************************************************
C
C  Estimate the smallest non-zero eigenvalue and the associated
C  eigenvector by repeated iteration.
C
      SUBROUTINE MYITEIG(NP,VECO,EVALMIN,EVALMAX)
      USE commons
      IMPLICIT NONE

      INTEGER NMAX, J1, J2, NP, NOPT, JV
      PARAMETER (NMAX=3*MXATMS)
      DOUBLE PRECISION DUMMY1,VEC1(3*MXATMS),PERCENT,
     1                 VEC2(3*MXATMS),DUMMY2,EVALMAX,EVALMIN,
     2                 VECO(3*MXATMS,20)
      SAVE

      NOPT=3*NATOMS
C
C  First estimate the largest eigenvalue of HESS. 
C
      IF (VECO(1,NVECTORS).EQ.0.0D0) THEN
         DO J1=1,NOPT
            VEC1(J1)=0.0D0
         ENDDO
         VEC1(1)=1.0D0
      ELSE
         DO J1=1,NOPT
            VEC1(J1)=VECO(J1,NVECTORS)
         ENDDO
      ENDIF

      DUMMY1=-1.0D10
      DO J1=1,NEVL
         CALL ORTHOG(VEC1,NP,VECO,0,NOPT,NATOMS,NMAX,.TRUE.)
         CALL DSYMV('U',NOPT,1.0D0,HESS,NMAX,VEC1,1,0.0D0,VEC2,1)
         DUMMY2=0.0D0
         DO J2=1,NOPT
            DUMMY2=DUMMY2+VEC2(J2)**2
         ENDDO
         EVALMAX=DSQRT(DUMMY2)
         DO J2=1,NOPT
            VEC2(J2)=VEC2(J2)/EVALMAX
         ENDDO
         PERCENT=100.0D0*DABS((DUMMY1-EVALMAX)/DUMMY1)
         IF (100.0D0*DABS((DUMMY1-EVALMAX)/DUMMY1).LT.CEIG) THEN
            WRITE(*,'(A,I4,A,F15.7)') ' Largest  eigenvalue converged in ',J1,' steps. Eigenvalue=',EVALMAX
            GOTO 10
         ENDIF
         DUMMY1=EVALMAX

         DO J2=1,NOPT
            VEC1(J2)=VEC2(J2)
         ENDDO
      ENDDO

      WRITE(*,'(A,F15.7)') ' ****WARNING - Largest  eigenvalue did not converge, value=',EVALMAX
10    DO J2=1,NOPT
         VECO(J2,NVECTORS)=VEC2(J2)
      ENDDO
C
C  Shift all the eigenvalues according to the size of EVALMAX.
C
      DO J1=1,NOPT
         HESS(J1,J1)=HESS(J1,J1)-EVALMAX*0.55D0
      ENDDO
C
C  Now estimate the new largest magnitude eigenvalue of HESS.
C
      DO JV=1,NVECTORS-1
         IF (VECO(1,JV).EQ.0.0D0) THEN
            DO J1=1,NOPT
               VEC1(J1)=0.0D0
            ENDDO
            VEC1(1)=1.0D0
         ELSE
            DO J1=1,NOPT
               VEC1(J1)=VECO(J1,JV)
            ENDDO
         ENDIF
         DUMMY1=-1.0D10
         DO J1=1,NEVS
            CALL ORTHOG(VEC1,NP,VECO,JV-1,NOPT,NATOMS,NMAX,.TRUE.)
            CALL DSYMV('U',NOPT,1.0D0,HESS,NMAX,VEC1,1,0.0D0,VEC2,1)
   
            DUMMY2=0.0D0
            DO J2=1,NOPT
               DUMMY2=DUMMY2+VEC2(J2)**2
            ENDDO
            EVALMIN=DSQRT(DUMMY2)
C
C  Normalize. Changing the phase may be necessary for extrapolation to work.
C
            DO J2=1,NOPT
               VEC2(J2)=-VEC2(J2)/EVALMIN
            ENDDO
            EVALMIN=-EVALMIN+EVALMAX*0.55D0
            PERCENT=100.0D0*DABS((DUMMY1-EVALMIN)/DUMMY1)
            IF (PERCENT.LT.CEIG) THEN
               WRITE(*,'(A,I4,A,F15.7)') ' Smallest eigenvalue converged in ',J1,' steps. Eigenvalue=',EVALMIN
               GOTO 20
            ENDIF
            DUMMY1=EVALMIN
            DO J2=1,NOPT
               VEC1(J2)=VEC2(J2)
            ENDDO
         ENDDO
         WRITE(*,'(A,F15.7)') ' ****WARNING - Smallest eigenvalue did not converge, value=',EVALMIN


20       DO J2=1,NOPT
            VECO(J2,JV)=VEC2(J2)
         ENDDO
      ENDDO


      RETURN
      END

C**********************************************************************************
C
C  Orthogonalise VEC1 to overall translations and rotations.
C
      SUBROUTINE ORTHOG(VEC1,NP,VECO,NVECO,NOPT,NREAL,NMAX,OTEST) 
      USE commons
      IMPLICIT NONE
      
      INTEGER J2, J3, NOPT, NREAL, NMAX, NVECO, NP
      DOUBLE PRECISION VEC1(NMAX), DUMMY1, DUMMY2, VECO(3*MXATMS,20)
      LOGICAL OTEST

      DO J2=1,NVECO
         DUMMY1=0.0D0
         DO J3=1,NOPT
         DUMMY1=DUMMY1+VECO(J3,J2)*VEC1(J3)
         ENDDO
         DUMMY2=0.0D0
         DO J3=1,NOPT
            VEC1(J3)=VEC1(J3)-DUMMY1*VECO(J3,J2)
            DUMMY2=DUMMY2+VEC1(J3)**2
         ENDDO
         DUMMY2=1.0D0/DSQRT(DUMMY2)
         DO J3=1,NOPT
            VEC1(J3)=VEC1(J3)*DUMMY2
         ENDDO
      ENDDO
      CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DO J2=1,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/DFLOAT(NREAL)
      DO J2=1,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DO J2=2,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/DFLOAT(NREAL)
      DO J2=2,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DO J2=3,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/DFLOAT(NREAL)
      DO J2=3,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      IF (PERIODIC) RETURN

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NREAL
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-1)*COORDS(J3,NP)-VEC1(J3)*COORDS(J3-1,NP)
         DUMMY2=DUMMY2+COORDS(J3,NP)**2+COORDS(J3-1,NP)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NREAL
            J3=3*J2
            VEC1(J3-1)=VEC1(J3-1)-DUMMY2*COORDS(J3,NP)
            VEC1(J3)=VEC1(J3)+DUMMY2*COORDS(J3-1,NP)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NREAL
         J3=3*J2
         DUMMY1=DUMMY1-VEC1(J3-2)*COORDS(J3,NP)+VEC1(J3)*COORDS(J3-2,NP)
         DUMMY2=DUMMY2+COORDS(J3,NP)**2+COORDS(J3-2,NP)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NREAL
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)+DUMMY2*COORDS(J3,NP)
            VEC1(J3)=VEC1(J3)-DUMMY2*COORDS(J3-2,NP)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NREAL
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-2)*COORDS(J3-1,NP)-VEC1(J3-1)*COORDS(J3-2,NP)
         DUMMY2=DUMMY2+COORDS(J3,NP)**2+COORDS(J3-2,NP)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NREAL
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)-DUMMY2*COORDS(J3-1,NP)
            VEC1(J3-1)=VEC1(J3-1)+DUMMY2*COORDS(J3-2,NP)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      RETURN
      END

C
C  Normalize vector VEC1
C
      SUBROUTINE  VECNORM(VEC1,NOPT)
      IMPLICIT NONE
      INTEGER J2, NOPT
      DOUBLE PRECISION DUMMY2, VEC1(NOPT)

      DUMMY2=0.0D0
      DO J2=1,NOPT
         DUMMY2=DUMMY2+VEC1(J2)**2
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=1.0D0/DSQRT(DUMMY2)
         DO J2=1,NOPT
            VEC1(J2)=VEC1(J2)*DUMMY2
         ENDDO
      ENDIF

      RETURN
      END
