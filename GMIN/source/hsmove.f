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
      SUBROUTINE HSMOVE(LCOORDSN,NP,NMOVE)
      USE commons
      IMPLICIT NONE
      INTEGER J1, J2, PARTNR(NATOMS), J3, NUMBERS(NATOMS), NPAIR, NP, NMOVE
      DOUBLE PRECISION FIXDIR(3*NATOMS), BOXL, COLTIM(NATOMS), TIMBIG, RX12, RY12, RZ12, MINIM, 
     1                 VX12, VY12, VZ12, B12, R12SQ, V12SQ, DISCR, SIGSQ, T12, COLTIM2(NATOMS), T122, DNEW, 
     2                 LCOORDSN(3*NATOMS,NPAR), DPRAND

      BOXL=BOXLX
      DO J1=1,3*NATOMS
         FIXDIR(J1)=2*(DPRAND()-0.5D0)
      ENDDO

      TIMBIG=1.0D100
      IF (.NOT.BINARY) SIGSQ=1.0D0

      DO J1=1,NATOMS
         COLTIM(J1)=TIMBIG
      ENDDO

C     DO J1=1,3*NATOMS
C        LCOORDSN(J1,NP)=LCOORDS(J1,NP)
C     ENDDO

      DO J1=1,NATOMS
         DO J2=J1+1,NATOMS

            IF (PERIODIC) THEN
               RX12=MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)
               RY12=MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)
               RZ12=MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)
            ELSE
               RX12=LCOORDSN(3*(J1-1)+1,NP)-LCOORDSN(3*(J2-1)+1,NP)
               RY12=LCOORDSN(3*(J1-1)+2,NP)-LCOORDSN(3*(J2-1)+2,NP)
               RZ12=LCOORDSN(3*(J1-1)+3,NP)-LCOORDSN(3*(J2-1)+3,NP)
            ENDIF

            VX12=FIXDIR(3*(J1-1)+1)-FIXDIR(3*(J2-1)+1)
            VY12=FIXDIR(3*(J1-1)+2)-FIXDIR(3*(J2-1)+2)
            VZ12=FIXDIR(3*(J1-1)+3)-FIXDIR(3*(J2-1)+3)

            B12=RX12*VX12+RY12*VY12+RZ12*VZ12

            IF (B12.LT.0.0D0) THEN
               R12SQ=RX12**2+RY12**2+RZ12**2
               V12SQ=VX12**2+VY12**2+VZ12**2
               IF (BINARY .OR. SOFT_SPHERE) THEN
                  IF (J1.LE.NTYPEA) THEN
                     IF (J2.LE.NTYPEA) THEN
                        DISCR=B12**2-V12SQ*(R12SQ-1.0D0)
                     ELSE
                        DISCR=B12**2-V12SQ*(R12SQ-SIGAB**2)
                     ENDIF
                  ELSE
                     DISCR=B12**2-V12SQ*(R12SQ-SIGBB**2)
                  ENDIF
               ELSE
                  DISCR=B12**2-V12SQ*(R12SQ-SIGSQ)
               ENDIF
               IF (DISCR.GT.0.0D0) THEN
                  T12=(-B12-SQRT(DISCR))/V12SQ
   
                  IF (T12.LT.COLTIM(J1)) THEN
                     COLTIM(J1)=T12
                     COLTIM2(J1)=(-B12+SQRT(DISCR))/V12SQ
                     PARTNR(J1)=J2
                  ENDIF
   
                  IF (T12.LT.COLTIM(J2)) THEN
                     COLTIM(J2)=T12
                     COLTIM2(J2)=(-B12+SQRT(DISCR))/V12SQ
                     PARTNR(J2)=J1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO J1=1,NATOMS
         NUMBERS(J1)=J1
      ENDDO
      CALL SORTN(NATOMS,COLTIM,NUMBERS)
C     T12=TIMBIG
C     DO J1=1,NATOMS
C        IF (COLTIM(J1).LT.T12) THEN
C           T12=COLTIM(J1)
C           T122=COLTIM2(J1)
C           J2=J1
C        ENDIF
C     ENDDO

      DO NPAIR=1,NMOVE

      J2=NUMBERS(NPAIR)
      J1=PARTNR(J2)
      T12=COLTIM(NPAIR)
      T122=COLTIM2(J2)

      IF (DEBUG) WRITE(*,'(A,I4,A,I4,A,I4,A,G15.6)') ' collision ',NPAIR,' is between atoms ',J1,' and ',J2,' at time ',T12
C
C  Advance all positions to time T12*T12FAC. 
C
      IF (T12FAC.LE.1.0D0) THEN
         DO J3=1,NATOMS
            COLTIM(J3)=COLTIM(J3)-T12*T12FAC
            LCOORDSN(3*(J3-1)+1,NP)=LCOORDSN(3*(J3-1)+1,NP)+FIXDIR(3*(J3-1)+1)*T12*T12FAC
            LCOORDSN(3*(J3-1)+2,NP)=LCOORDSN(3*(J3-1)+2,NP)+FIXDIR(3*(J3-1)+2)*T12*T12FAC
            LCOORDSN(3*(J3-1)+3,NP)=LCOORDSN(3*(J3-1)+3,NP)+FIXDIR(3*(J3-1)+3)*T12*T12FAC
C
C  Could put atoms leaving the primary supercell back in the box here
C  if desired.
C
         ENDDO
C
C  Advance positions of the colliding pair only to time T12. 
C
C        LCOORDSN(3*(J1-1)+1,NP)=LCOORDSN(3*(J1-1)+1,NP)+FIXDIR(3*(J1-1)+1)*T12
C        LCOORDSN(3*(J1-1)+2,NP)=LCOORDSN(3*(J1-1)+2,NP)+FIXDIR(3*(J1-1)+2)*T12
C        LCOORDSN(3*(J1-1)+3,NP)=LCOORDSN(3*(J1-1)+3,NP)+FIXDIR(3*(J1-1)+3)*T12
C        LCOORDSN(3*(J2-1)+1,NP)=LCOORDSN(3*(J2-1)+1,NP)+FIXDIR(3*(J2-1)+1)*T12
C        LCOORDSN(3*(J2-1)+2,NP)=LCOORDSN(3*(J2-1)+2,NP)+FIXDIR(3*(J2-1)+2)*T12
C        LCOORDSN(3*(J2-1)+3,NP)=LCOORDSN(3*(J2-1)+3,NP)+FIXDIR(3*(J2-1)+3)*T12

C        IF (PERIODIC) THEN
C           PRINT*,'separation between these atoms is now ',SQRT(MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)**2
C    1                                      +MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)**2
C    2                                      +MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)**2)
C        ELSE
C           PRINT*,'separation between these atoms is now ',SQRT((LCOORDSN(3*(J1-1)+1,NP)-LCOORDSN(3*(J2-1)+1,NP))**2
C    1                                      +(LCOORDSN(3*(J1-1)+2,NP)-LCOORDSN(3*(J2-1)+2,NP))**2
C    2                                      +(LCOORDSN(3*(J1-1)+3,NP)-LCOORDSN(3*(J2-1)+3,NP))**2)
C        ENDIF
      ELSE 
C
C  Put the first colliding pair half way between their entrance/exit positions
C

         LCOORDSN(3*(J1-1)+1,NP)=LCOORDSN(3*(J1-1)+1,NP)+FIXDIR(3*(J1-1)+1)*(T122-T12)/2
         LCOORDSN(3*(J1-1)+2,NP)=LCOORDSN(3*(J1-1)+2,NP)+FIXDIR(3*(J1-1)+2)*(T122-T12)/2
         LCOORDSN(3*(J1-1)+3,NP)=LCOORDSN(3*(J1-1)+3,NP)+FIXDIR(3*(J1-1)+3)*(T122-T12)/2
         LCOORDSN(3*(J2-1)+1,NP)=LCOORDSN(3*(J2-1)+1,NP)+FIXDIR(3*(J2-1)+1)*(T122-T12)/2
         LCOORDSN(3*(J2-1)+2,NP)=LCOORDSN(3*(J2-1)+2,NP)+FIXDIR(3*(J2-1)+2)*(T122-T12)/2
         LCOORDSN(3*(J2-1)+3,NP)=LCOORDSN(3*(J2-1)+3,NP)+FIXDIR(3*(J2-1)+3)*(T122-T12)/2

         IF (PERIODIC) THEN
            DNEW=SQRT(MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)**2
     1            +MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)**2
     2            +MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)**2)
         ELSE
            DNEW=SQRT((LCOORDSN(3*(J1-1)+1,NP)-LCOORDSN(3*(J2-1)+1,NP))**2
     1            +(LCOORDSN(3*(J1-1)+2,NP)-LCOORDSN(3*(J2-1)+2,NP))**2
     2            +(LCOORDSN(3*(J1-1)+3,NP)-LCOORDSN(3*(J2-1)+3,NP))**2)
         ENDIF
C        PRINT*,'separation half way between entrance and exit=',DNEW

         IF (BINARY) THEN
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  SIGSQ=1.0D0
               ELSE
                  SIGSQ=SIGAB**2
               ENDIF
            ELSE
               SIGSQ=SIGBB**2
            ENDIF
         ENDIF
C
C  Rescale the distance between the first colliding pair
C
         IF (PERIODIC) THEN
            RX12=MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)*SQRT(SIGSQ)/DNEW/2
            RY12=MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)*SQRT(SIGSQ)/DNEW/2
            RZ12=MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)*SQRT(SIGSQ)/DNEW/2

            VX12=MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)
            VY12=MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)
            VZ12=MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)

            LCOORDSN(3*(J1-1)+1,NP)=LCOORDSN(3*(J1-1)+1,NP)-VX12/2+RX12
            LCOORDSN(3*(J1-1)+2,NP)=LCOORDSN(3*(J1-1)+2,NP)-VY12/2+RY12
            LCOORDSN(3*(J1-1)+3,NP)=LCOORDSN(3*(J1-1)+3,NP)-VZ12/2+RZ12
            LCOORDSN(3*(J2-1)+1,NP)=LCOORDSN(3*(J2-1)+1,NP)+VX12/2-RX12
            LCOORDSN(3*(J2-1)+2,NP)=LCOORDSN(3*(J2-1)+2,NP)+VY12/2-RY12
            LCOORDSN(3*(J2-1)+3,NP)=LCOORDSN(3*(J2-1)+3,NP)+VZ12/2-RZ12

            DNEW=SQRT(MINIM(LCOORDSN(3*(J1-1)+1,NP),LCOORDSN(3*(J2-1)+1,NP),BOXL)**2
     1            +MINIM(LCOORDSN(3*(J1-1)+2,NP),LCOORDSN(3*(J2-1)+2,NP),BOXL)**2
     2            +MINIM(LCOORDSN(3*(J1-1)+3,NP),LCOORDSN(3*(J2-1)+3,NP),BOXL)**2)
         ELSE
            RX12=(LCOORDSN(3*(J1-1)+1,NP)-LCOORDSN(3*(J2-1)+1,NP))*SQRT(SIGSQ)/DNEW/2
            RY12=(LCOORDSN(3*(J1-1)+2,NP)-LCOORDSN(3*(J2-1)+2,NP))*SQRT(SIGSQ)/DNEW/2
            RZ12=(LCOORDSN(3*(J1-1)+3,NP)-LCOORDSN(3*(J2-1)+3,NP))*SQRT(SIGSQ)/DNEW/2
            VX12=(LCOORDSN(3*(J1-1)+1,NP)+LCOORDSN(3*(J2-1)+1,NP))/2
            VY12=(LCOORDSN(3*(J1-1)+2,NP)+LCOORDSN(3*(J2-1)+2,NP))/2
            VZ12=(LCOORDSN(3*(J1-1)+3,NP)+LCOORDSN(3*(J2-1)+3,NP))/2
            LCOORDSN(3*(J1-1)+1,NP)=VX12+RX12
            LCOORDSN(3*(J1-1)+2,NP)=VY12+RY12
            LCOORDSN(3*(J1-1)+3,NP)=VZ12+RZ12

            LCOORDSN(3*(J2-1)+1,NP)=VX12-RX12
            LCOORDSN(3*(J2-1)+2,NP)=VY12-RY12
            LCOORDSN(3*(J2-1)+3,NP)=VZ12-RZ12

            DNEW=SQRT((LCOORDSN(3*(J1-1)+1,NP)-LCOORDSN(3*(J2-1)+1,NP))**2
     1            +(LCOORDSN(3*(J1-1)+2,NP)-LCOORDSN(3*(J2-1)+2,NP))**2
     2            +(LCOORDSN(3*(J1-1)+3,NP)-LCOORDSN(3*(J2-1)+3,NP))**2)
         ENDIF
C        PRINT*,'separation reset to ',DNEW
      ENDIF

C     LSTEP=0.0D0
C     DO J3=1,3*NATOMS
C        FIXDIR(J3)=LCOORDSN(J3,NP)-LCOORDS(J3,NP)
C        LSTEP=LSTEP+(LCOORDSN(J3,NP)-LCOORDS(J3,NP))**2
C     ENDDO
C     LSTEP=SQRT(LSTEP)
C     CALL ORTHOG(FIXDIR,LCOORDS,NATOMS,3*NATOMS,.TRUE.)

      ENDDO
     
      IF (DEBUG) WRITE(*,'(I4,A)') NMOVE,' hard sphere type collisions completed'
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION MINIM(A,B,BOXL)
      DOUBLE PRECISION A, B, BOXL

      MINIM=A-B-BOXL*ANINT((A-B)/BOXL)

      RETURN
      END

      SUBROUTINE SORTN(N,A,NA)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      INTEGER J1, L, N, J2, NA(N), NTEMP
      DOUBLE PRECISION TEMP, A(N)
C
      DO J1=1,N-1
         L=J1
         DO J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
         ENDDO
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
      ENDDO
      RETURN
      END
