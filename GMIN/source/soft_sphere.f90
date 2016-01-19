!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*************************************************************************
!
!  Subroutine SOFT_SPHERE calculates the energy, cartesian gradient and second
!  derivative matrix analytically for the soft sphere potential for a binary
!  liquid in reduced units using a shifted, truncated potential.  
!
!*************************************************************************
!
      MODULE SOFT_SPHERE_CLASS
         !this module defines the potential SOFT_SPHERE
         IMPLICIT NONE
         PRIVATE

         TYPE INTERACTION_DEF
           DOUBLE PRECISION SIG, RC, A, CIJ, BIJ, SIG2, ISIG2, A2, RC2
         END TYPE INTERACTION_DEF
         TYPE(INTERACTION_DEF) :: AA, BB, AB
         DOUBLE PRECISION SIZERATIO
         INTEGER N

         PUBLIC :: SOFT_SPHERE_CLASS_SETUP, &
     &   SOFT_SPHERE_UPDATE_E_AA, SOFT_SPHERE_UPDATE_E_BB, SOFT_SPHERE_UPDATE_E_AB, &
     &   SOFT_SPHERE_UPDATE_EG_AA, SOFT_SPHERE_UPDATE_EG_BB, SOFT_SPHERE_UPDATE_EG_AB, &
     &   SOFT_SPHERE_GET_CUT

         CONTAINS

         !define the potential
         FUNCTION VIJ( R, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: R
            DOUBLE PRECISION VAL
            IF (R < T%RC) THEN
              VAL = (T%SIG/R)**12 - T%CIJ
            ELSE IF (R < T%A) THEN
              VAL = T%BIJ*(T%A-R)**3 
            ELSE 
              VAL = 0.D0
            ENDIF
         END FUNCTION VIJ

         !define the first derivative of the potential
         !To save some computation, return the gradient / R
         FUNCTION DVIJ( R, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: R
            DOUBLE PRECISION VAL
            IF (R < T%RC) THEN
              VAL = -12.D0 * (T%SIG/R)**14/T%SIG**2
            ELSE IF (R < T%A) THEN
              VAL = -3.D0 * T%BIJ*(T%A-R)**2 /R
            ELSE 
              VAL = 0.D0
            ENDIF
         END FUNCTION DVIJ

         !define the second derivative of the potential
         FUNCTION DDVIJ( R, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: R
            DOUBLE PRECISION VAL
            IF (R < T%RC) THEN
              VAL = 156.D0 * (T%SIG/R)**14/(T%SIG**2)
            ELSE IF (R < T%A) THEN
              VAL = 6.D0 * T%BIJ*(T%A-R)
            ELSE 
              VAL = 0.D0
            ENDIF
         END FUNCTION DDVIJ



         SUBROUTINE GET_R2(X, J1, J2, R2, XVEC)
            USE COMMONS, ONLY : BOXLX, BOXLY, BOXLZ
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(OUT) :: R2, XVEC(3)
            DOUBLE PRECISION, INTENT(IN) :: X(N)
            INTEGER, INTENT(IN) :: J1, J2
            INTEGER J3, J4
            J3=3*(J1-1)
            J4=3*(J2-1)
            !calculate atom separation
            XVEC(1)=X(J3+1)-X(J4+1)
            XVEC(2)=X(J3+2)-X(J4+2)
            XVEC(3)=X(J3+3)-X(J4+3)
            XVEC(1)=XVEC(1)-BOXLX*NINT(XVEC(1)/BOXLX)
            XVEC(2)=XVEC(2)-BOXLY*NINT(XVEC(2)/BOXLX)
            XVEC(3)=XVEC(3)-BOXLZ*NINT(XVEC(3)/BOXLX)
            R2=(XVEC(1)**2 + XVEC(2)**2 + XVEC(3)**2)
            !if (abs(xvec(1)) .gt. 0.5*boxlx .or.(j1.eq.1 .and. j2.eq.111)) then
            !if ((j1.eq.1 .and. j2.eq.111)) then
               !write(*,*) j1, j2, xvec(1), fiximage, x(j3+1), x(j4+1), NINT((X(J3+1)-X(J4+1))/BOXLX)
            !endif
         END SUBROUTINE GET_R2

         SUBROUTINE UPDATE_POTENTIAL( X, J1, J2, POTEL, T )
            IMPLICIT NONE
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            DOUBLE PRECISION R2, XVEC(3)
            CALL GET_R2(X, J1, J2, R2, XVEC)
            IF (R2 .GE. T%A2) THEN
               RETURN
            ELSE IF (R2 .GE. T%RC2) THEN
               POTEL = POTEL + T%BIJ * (T%A - SQRT(R2))**3 
            ELSE 
               POTEL = POTEL + (T%SIG2/R2)**6 - T%CIJ
            ENDIF
            !E = VIJ(sqrt(R2), T)  !avoid the square root
            !if (j2 .lt. 204 .and. j1 .lt. 204) THEN
               !write(*,*) sqrt(R2), E, T%SIG, j1, j2, AA%sig, BB%sig, AB%sig
            !endif
            !POTEL = POTEL + E
         END SUBROUTINE UPDATE_POTENTIAL

         SUBROUTINE UPDATE_POTENTIAL_GRADIENT( X, J1, J2, POTEL, V, T )
            IMPLICIT NONE
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            DOUBLE PRECISION R, R2, XVEC(3), G
            INTEGER J5
            CALL GET_R2(X, J1, J2, R2, XVEC)
            IF (R2 .GE. T%A2) THEN
               RETURN
            ELSE IF (R2 .GE. T%RC2) THEN
               POTEL = POTEL + T%BIJ * (T%A - SQRT(R2))**3 
               R = SQRT(R2)
               G = -3.D0 * T%BIJ*(T%A-R)**2 /R
            ELSE 
               POTEL = POTEL + (T%SIG2/R2)**6 - T%CIJ
               G = -12.D0 * (T%SIG2/R2)**7 * T%ISIG2
            ENDIF
            !update gradient
            !G = DVIJ(sqrt(R2), T)  !avoid the square root
            DO J5=1,3
              V(3*(J1-1)+J5)=V(3*(J1-1)+J5)+G*XVEC(J5)
              V(3*(J2-1)+J5)=V(3*(J2-1)+J5)-G*XVEC(J5)
            END DO
            !if (sqrt(R2) < T%A) write(*,*) sqrt(R2), E, G
         END SUBROUTINE UPDATE_POTENTIAL_GRADIENT

         SUBROUTINE SET_INTERACTION_DEF( T, sig1, sig2, rc )
            IMPLICIT NONE
            !careful, with INTENT(OUT) here, any element of T not assigned in this subroutine could become undefined
            TYPE(INTERACTION_DEF), INTENT(OUT) ::  T
            DOUBLE PRECISION, INTENT(IN) ::  SIG1, SIG2, RC
            T%SIG = SIG1 + SIG2
            T%RC = RC !SHOULD I MULTIPLY BY SIG HERE?
            T%A = T%RC * (1.D0+2.D0/13)
            T%BIJ = 13.D0**2 * T%SIG**12 / T%RC**15
            T%CIJ = (T%SIG/T%RC)**12 * (1.D0-8.D0/13)

            T%A2 = T%A**2
            T%RC2 = T%RC**2
            T%SIG2 = T%SIG**2
            T%ISIG2 = 1.D0 / T%SIG2
         END SUBROUTINE SET_INTERACTION_DEF

         !Above this will be only general subroutines and functions.
         !Below this will be subroutines which depend on the specific type of
         !interaction, AA, AB, BB.

         SUBROUTINE SOFT_SPHERE_GET_CUT( CUTAA, CUTBB, CUTAB )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(OUT) :: CUTAA, CUTBB, CUTAB 
            CUTAA = AA%A
            CUTBB = BB%A
            CUTAB = AB%A
         END SUBROUTINE SOFT_SPHERE_GET_CUT

         SUBROUTINE SOFT_SPHERE_CLASS_SETUP(CUTOFF)
            USE COMMONS, ONLY : NATOMS
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: CUTOFF
            DOUBLE PRECISION SIGA, SIGB
            LOGICAL :: FIRST_THROUGH = .TRUE. 
            !FIRST_THROUGH will be defined as .TRUE. on first initialization, but will retain value on subsequent calls

            IF ( .NOT. FIRST_THROUGH ) RETURN
            FIRST_THROUGH = .FALSE.

            N = NATOMS

            !In principle SIG and rc and be user defined for each interaction
            !type.  At the moment, this is all predefined

            SIZERATIO = 1.2D0 !SIGB = SIGA*SIZERATIO
            SIGB = ( 4.D0 / ( 8.D0 +2.D0*(1.D0+SIZERATIO)**3 + 8.D0*SIZERATIO**3 ))**(1./3)
            SIGA = SIZERATIO * SIGB

            CALL SET_INTERACTION_DEF(AA, SIGA, SIGA, CUTOFF)
            CALL SET_INTERACTION_DEF(BB, SIGB, SIGB, CUTOFF)
            CALL SET_INTERACTION_DEF(AB, SIGA, SIGB, CUTOFF)

            write(*,*) "#SIG ", AA%SIG, BB%SIG, AB%SIG
            write(*,*) "#A ", AA%A, BB%A, AB%A
            write(*,*) "#Bij ", AA%Bij, BB%Bij, AB%Bij
            write(*,*) "#Cij ", AA%Cij, BB%Cij, AB%Cij
         END SUBROUTINE SOFT_SPHERE_CLASS_SETUP

         !3 subroutines for calculating the potentials
         SUBROUTINE SOFT_SPHERE_UPDATE_E_AA( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, AA)
         END SUBROUTINE SOFT_SPHERE_UPDATE_E_AA
         SUBROUTINE SOFT_SPHERE_UPDATE_E_BB( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, BB)
         END SUBROUTINE SOFT_SPHERE_UPDATE_E_BB
         SUBROUTINE SOFT_SPHERE_UPDATE_E_AB( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, AB)
         END SUBROUTINE SOFT_SPHERE_UPDATE_E_AB

         !3 subroutines for calculating the potential and gradients
         SUBROUTINE SOFT_SPHERE_UPDATE_EG_AA( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, AA)
         END SUBROUTINE SOFT_SPHERE_UPDATE_EG_AA
         SUBROUTINE SOFT_SPHERE_UPDATE_EG_BB( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, BB)
         END SUBROUTINE SOFT_SPHERE_UPDATE_EG_BB
         SUBROUTINE SOFT_SPHERE_UPDATE_EG_AB( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, AB)
         END SUBROUTINE SOFT_SPHERE_UPDATE_EG_AB

      END MODULE SOFT_SPHERE_CLASS

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

      SUBROUTINE SOFT_SPHERE_POT(X, V, POTEL, GTEST, STEST)
      USE COMMONS, ONLY : NATOMS, BOXLX, BOXLY, BOXLZ, CUTOFF, FIXIMAGE, NORESET, &
     &    FREEZE, FREEZESAVE, FREEZESAVEE, NFREEZE, NFREEZETYPEA, FROZEN, FROZENLIST, &
     &    FREEZEIL, FREEZEIL_USE, FREEZEIL_NAA, FREEZEIL_NAB, FREEZEIL_NBB, FREEZEIL_E, &
     &    RESTRICTREGIONRADIUS, RESTRICTREGIONX0, RESTRICTREGIONY0, RESTRICTREGIONZ0, &
     &    RESTRICTREGION, SOFT_SPHERE_NTYPEA
      USE SOFT_SPHERE_CLASS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, NNOFREEZETYPEA
      DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
      DOUBLE PRECISION CUTAA, CUTAB, CUTBB, POTELDUM
      LOGICAL, INTENT(IN) :: GTEST, STEST
      INTEGER AASTART, ABSTART, BBSTART
      DOUBLE PRECISION RRX, RRY, RRZ, RRR1, RRR2, CUT
      LOGICAL USEIL
      double precision x111

      NTYPEA = SOFT_SPHERE_NTYPEA
      CALL SOFT_SPHERE_CLASS_SETUP( CUTOFF )
      CALL SOFT_SPHERE_GET_CUT( CUTAA, CUTBB, CUTAB )

      !WRITE(*,*) "FIXIMAGE ", FIXIMAGE, "SETTING FALSE"
      !FIXIMAGE = .FALSE.


      N = NATOMS
      NNOFREEZETYPEA = NTYPEA - NFREEZETYPEA
!
!  Deal with any atoms that have left the box.
!
      !x111 = x(111*3-2)
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
      !if ( abs(x111 - x(111*3-2)) .gt. 1.) then
         !write (*,*) "changing x111", x111, "->", x(111*3-2)
      !endif
      !write(*,*) X(1), X(2), X(3)
!
!  Calculate interatomic vectors using the minimum image convention.
!  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
!  atoms i and j.
!

      POTEL=0.0D0
      IF (STEST) THEN
        write(*,*) "warning: calculation of the hessian is not implimented in GMIN"
        ! js850> it could easily be implimented, but there is no way to return
        ! the calculated matrix.  OPTIM usis a module which is not implimented
        ! (and not needed) in GMIN
      ENDIF

      IF (.True.) THEN !bug testing only
      IF (FREEZE .AND. RESTRICTREGION) THEN
        !choose AASTART, ABSTART, BBSTART so that all interactions can fit in
        !the array FREEZEIL.
        AASTART = 1;
        BBSTART = AASTART+NTYPEA*(NTYPEA-1);
        ABSTART = BBSTART+(NATOMS-NTYPEA)*(NATOMS-NTYPEA-1);

        IF ( .NOT. FREEZEIL_USE) THEN
          !allocate memory here.
          ALLOCATE(FREEZEIL(NATOMS*(NATOMS+1),2))

          !
          !set up FREEZEIL
          !
          FREEZEIL_USE = .TRUE.
          FREEZEIL_NAA = 0
          FREEZEIL_NBB = 0
          FREEZEIL_NAB = 0
          FREEZEIL(:,:) = 0
          DO J1=1,NATOMS
            DO J2=J1+1,NATOMS
              IF (J1.LE.NTYPEA .AND. J2.LE.NTYPEA) THEN
                CUT = CUTAA
              ELSE IF (J1.GT.NTYPEA .AND. J2.GT.NTYPEA) THEN
                CUT = CUTBB
              ELSE
                CUT = CUTAB
              ENDIF
              RRX = ( X(3*(J1-1)+1)-RESTRICTREGIONX0 )
              RRY = ( X(3*(J1-1)+2)-RESTRICTREGIONY0 )
              RRZ = ( X(3*(J1-1)+3)-RESTRICTREGIONZ0 )
              RRR1 = DSQRT(RRX**2+RRY**2+RRZ**2 )
              RRX = ( X(3*(J2-1)+1)-RESTRICTREGIONX0 )
              RRY = ( X(3*(J2-1)+2)-RESTRICTREGIONY0 )
              RRZ = ( X(3*(J2-1)+3)-RESTRICTREGIONZ0 )
              RRR2 = DSQRT(RRX**2+RRY**2+RRZ**2 )
              USEIL = .TRUE.
              IF (RRR1.GT.RESTRICtREGIONRADIUS+CUT) USEIL=.FALSE.
              IF (RRR2.GT.RESTRICtREGIONRADIUS+CUT) USEIL=.FALSE.
              IF (FROZEN(J1) .AND. FROZEN(J2)) USEIL=.FALSE.
              IF (USEIL) THEN
                !add to FREEZEIL
                IF (J1.LE.NTYPEA .AND. J2.LE.NTYPEA) THEN
                  FREEZEIL( AASTART + FREEZEIL_NAA ,1) = J1
                  FREEZEIL( AASTART + FREEZEIL_NAA ,2) = J2
                  FREEZEIL_NAA = FREEZEIL_NAA + 1
                ELSE IF (J1.GT.NTYPEA .AND. J2.GT.NTYPEA) THEN
                  FREEZEIL( BBSTART + FREEZEIL_NBB ,1) = J1
                  FREEZEIL( BBSTART + FREEZEIL_NBB ,2) = J2
                  FREEZEIL_NBB = FREEZEIL_NBB + 1
                else
                  FREEZEIL( ABSTART + FREEZEIL_NAB ,1) = J1
                  FREEZEIL( ABSTART + FREEZEIL_NAB ,2) = J2
                  FREEZEIL_NAB = FREEZEIL_NAB + 1
                endif

              ELSE
                !calc energy of j1,j2 interaction and add to FREEZEIL_E
                IF (J1.LE.NTYPEA .AND. J2.LE.NTYPEA) THEN
                  CALL SOFT_SPHERE_UPDATE_E_AA(X, J1, J2, FREEZEIL_E)
                ELSE IF (J1.GT.NTYPEA .AND. J2.GT.NTYPEA) THEN
                  CALL SOFT_SPHERE_UPDATE_E_BB(X, J1, J2, FREEZEIL_E)
                else
                  CALL SOFT_SPHERE_UPDATE_E_AB(X, J1, J2, FREEZEIL_E)
                endif
              ENDIF
            ENDDO
          ENDDO
        endif

        !write(*,*) FREEZEIL_NAA, FREEZEIL_NBB, FREEZEIL_NAB, FREEZEIL_NAA+ FREEZEIL_NBB+ FREEZEIL_NAB
        IF ((GTEST .OR. STEST)) THEN
          V(1:3*NATOMS)=0.D0
          do J3=AASTART,AASTART+FREEZEIL_NAA-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            !write (*,*) "AA", j1, j2, frozen(j1), frozen(j2)
            CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTEL, V)
          enddo
          do J3=BBSTART,BBSTART+FREEZEIL_NBB-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTEL, V)
          enddo
          do J3=ABSTART,ABSTART+FREEZEIL_NAB-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
          enddo
        ELSE
          do J3=AASTART,AASTART+FREEZEIL_NAA-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            !write (*,*) "AA", j1, j2, frozen(j1), frozen(j2)
            CALL SOFT_SPHERE_UPDATE_E_AA(X, J1, J2, POTEL)
          enddo
          do J3=BBSTART,BBSTART+FREEZEIL_NBB-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            CALL SOFT_SPHERE_UPDATE_E_BB(X, J1, J2, POTEL)
          enddo
          do J3=ABSTART,ABSTART+FREEZEIL_NAB-1
            j1=FREEZEIL(j3,1)
            j2=FREEZEIL(j3,2)
            CALL SOFT_SPHERE_UPDATE_E_AB(X, J1, J2, POTEL)
          enddo
        ENDIF
          POTEL = POTEL + FREEZEIL_E


      ELSEIF ((GTEST .OR. STEST) .AND. FREEZE) THEN
         !if GTEST, update V in adition to POTEL

         DO J1=1,3*NATOMS
           V(J1) = 0.D0
         END DO

         IF ( .NOT. FREEZESAVE ) THEN

           !calculate interactions between mobile particles
           !AA
           DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
             J1=FROZENLIST(J3)
             DO J4=J3+1,NFREEZE+NNOFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !BB
           DO J3=NFREEZE+NNOFREEZETYPEA+1,N
             J1=FROZENLIST(J3)
             DO J4=J3+1,N
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !AB
           DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
             J1=FROZENLIST(J3)
             DO J4=NFREEZE+NNOFREEZETYPEA+1,N
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !calculate interactions between mobile and imobile atoms
           !AA
           DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
             J1=FROZENLIST(J3)
             DO J4=1,NFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !BB
           DO J3=NFREEZE+NNOFREEZETYPEA+1,N
             J1=FROZENLIST(J3)
             DO J4=NFREEZETYPEA+1,NFREEZE
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !AB mobile-imobile
           DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
             J1=FROZENLIST(J3)
             DO J4=NFREEZETYPEA+1,NFREEZE
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO
           !BA mobile-imobile
           DO J3=NFREEZE+NNOFREEZETYPEA+1,N
             J1=FROZENLIST(J3)
             DO J4=1,NFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
             ENDDO
           ENDDO

           POTEL=POTEL+FREEZESAVEE

         ELSE

           FREEZESAVEE = 0.D0
           DO J1=1,NTYPEA
             DO J2=J1+1,NTYPEA
               POTELDUM=0.D0
               CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
             ENDDO
           ENDDO
           DO J1=1,NTYPEA
             DO J2=NTYPEA+1,N
               POTELDUM=0.D0
               CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
             ENDDO
           ENDDO
           DO J1=NTYPEA+1,N
             DO J2=J1+1,N
               POTELDUM=0.D0
               CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
             ENDDO
           ENDDO

           FREEZESAVE = .FALSE.

         ENDIF

      ELSEIF ( GTEST .OR. STEST ) THEN

        DO J1=1,3*NATOMS
          V(J1) = 0.D0
        END DO

        DO J1=1,NTYPEA
          DO J2=J1+1,NTYPEA
            CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=1,NTYPEA
          DO J2=NTYPEA+1,N
            CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=NTYPEA+1,N
          DO J2=J1+1,N
            CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO

      ELSE
         !only update POTEL
         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL SOFT_SPHERE_UPDATE_E_AA(X, J1, J2, POTEL)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1, N
               CALL SOFT_SPHERE_UPDATE_E_AB(X, J1, J2, POTEL)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, N
            DO J2=J1+1, N
               CALL SOFT_SPHERE_UPDATE_E_BB(X, J1, J2, POTEL)
            ENDDO
         ENDDO
      ENDIF
      !POTEL = POTEL * 4.D0
      ELSE  !for bugtesting only
        !write(*,*) "HERE"
        DO J1=1,3*NATOMS
          V(J1) = 0.D0
        END DO
        POTEL = 0.D0

        DO J1=1,NTYPEA
          DO J2=J1+1,NTYPEA
            CALL SOFT_SPHERE_UPDATE_EG_AA(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=1,NTYPEA
          DO J2=NTYPEA+1,N
            CALL SOFT_SPHERE_UPDATE_EG_AB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=NTYPEA+1,N
          DO J2=J1+1,N
            CALL SOFT_SPHERE_UPDATE_EG_BB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO

      ENDIF  !for bugtesting only

      RETURN
      END SUBROUTINE SOFT_SPHERE_POT

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
      
!*****************************************************************************

!*****************************************************************************

!*****************************************************************************

!*****************************************************************************


