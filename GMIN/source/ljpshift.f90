!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful, !   but WITHOUT ANY WARRANTY; without even the implied warranty of
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
!  Subroutine LJPSHIFT calculates the energy, cartesian gradient and second
!  derivative matrix analytically for Lennard-Jones in reduced units
!  (epsilon=sigma=1) using a shifted, truncated potential.
!
!  Adapted for the binary LJ glass described by Sastry, Debenetti and
!  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
!  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
!  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
!  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
!  of 5.975206 for 256 atoms. 
!
!
!*************************************************************************
!

      MODULE LJPSHIFT_CLASS
         !this module defines the potential LJPSHIFT
         IMPLICIT NONE
         PRIVATE

         TYPE INTERACTION_DEF
           DOUBLE PRECISION EPS, SIG, CONST, RCONST, RCUT, IRCUT2, SIG6, SIG12
         END TYPE INTERACTION_DEF
         TYPE(INTERACTION_DEF) :: AA, BB, AB
         DOUBLE PRECISION SIZERATIO
         INTEGER N
         DOUBLE PRECISION BOXLX, BOXLY, BOXLZ, BOXLMIN_HALF_2
         DOUBLE PRECISION IBOXLX, IBOXLY, IBOXLZ
         DOUBLE PRECISION BOXLVEC(3), IBOXLVEC(3)

         PUBLIC :: LJPSHIFT_CLASS_SETUP, &
     &   LJPSHIFT_UPDATE_E_AA, LJPSHIFT_UPDATE_E_BB, LJPSHIFT_UPDATE_E_AB, &
     &   LJPSHIFT_UPDATE_EG_AA, LJPSHIFT_UPDATE_EG_BB, LJPSHIFT_UPDATE_EG_AB, &
     &   LJPSHIFT_GET_CUT, &
     &   LJPSHIFT_UPDATE_E_AA_offset, LJPSHIFT_UPDATE_E_BB_offset, &
          LJPSHIFT_UPDATE_E_AB_offset

         CONTAINS

         !define the potential
         !assume R is not more than RCUT
         FUNCTION VIJ( R2, IR6, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: R2, IR6
            DOUBLE PRECISION VAL
            VAL = 4.D0*T%EPS*(T%SIG6*IR6*(T%SIG6*IR6-1.0D0) + T%RCONST*R2 + T%CONST)
         END FUNCTION VIJ

         !define the first derivative of the potential
         !assume R is not more than RCUT
         !To save some computation, return the gradient / R
         FUNCTION DVIJ( IR8, IR14, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: IR14, IR8
            DOUBLE PRECISION VAL
            VAL = -8.0D0*T%EPS*(3.0D0*(2.0D0*IR14*(T%SIG12)-IR8*T%SIG6)-T%RCONST)
         END FUNCTION DVIJ


         !define the second derivative of the potential
         FUNCTION DDVIJ( R, T ) RESULT(VAL)
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: R
            DOUBLE PRECISION VAL
            VAL = 0.D0 !not implimented
         END FUNCTION DDVIJ



         SUBROUTINE GET_R2(X, J1, J2, R2, XVEC)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(OUT) :: R2, XVEC(3)
            DOUBLE PRECISION, INTENT(IN) :: X(N)
            INTEGER, INTENT(IN) :: J1, J2
            INTEGER J3, J4
            J3=3*(J1-1)
            J4=3*(J2-1)
            !calculate atom separation
            XVEC(:) = X(J3+1:J3+3) - X(J4+1:J4+3)
            XVEC(:) = XVEC(:) - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
            R2 = sum(XVEC**2)

            !XVEC(:) = X(3*(J1-1)+1:3*(J1-1)+3) - X(3*(J2-1)+1:3*(J2-1)+3)
            !R2 = sum(XVEC**2)
            !IF (R2 .GE. BOXLMIN_HALF_2) THEN
               !XVEC(:) = XVEC(:) - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
               !R2 = SUM(XVEC**2)
            !ENDIF
         END SUBROUTINE GET_R2

         SUBROUTINE UPDATE_POTENTIAL( X, J1, J2, POTEL, T )
            IMPLICIT NONE
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            DOUBLE PRECISION R2, E, XVEC(3), IR2, IR6

            CALL GET_R2(X, J1, J2, R2, XVEC)
            !XVEC(:) = X(3*(J1-1)+1:3*(J1-1)+3) - X(3*(J2-1)+1:3*(J2-1)+3)
            !R2 = sum(XVEC**2)
            !IF (R2 .GE. BOXLMIN_HALF_2) THEN
            !   XVEC(:) = XVEC(:) - BOXLVEC(:) * NINT(XVEC(:)*IBOXLVEC(:))
            !   R2 = SUM(XVEC**2)
            !ENDIF
            
            IR2 = 1.D0/R2
            IF (IR2.GT.T%IRCUT2) THEN
               IR6=IR2**3
               E = VIJ(R2, IR6, T)
               !E = T%EPS*(T%SIG6*IR6*(T%SIG6*IR6-1.0D0) + T%RCONST/IR2 + T%CONST)
               POTEL = POTEL + E
            ENDIF
         END SUBROUTINE UPDATE_POTENTIAL

         SUBROUTINE UPDATE_POTENTIAL_offset( X, J1, J2, POTEL, T, xoffset )
            !POTENTIAL WHERE THE PERIODIC OFFSET XOFFSET IS SPECIFIED DIRECTLY
            !use NL_BIN_MOVEONE, only : nl_myunit
            IMPLICIT NONE
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: X(3*N), xoffset(3)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            DOUBLE PRECISION R2, E, IR2, IR6, r2test, drtest(3)

            !CALL GET_R2(X, J1, J2, R2, XVEC)
            R2 = sum((X(3*(J1-1)+1:3*(J1-1)+3) - X(3*(J2-1)+1:3*(J2-1)+3) &
                  +xoffset)**2)
            
            IR2 = 1.D0/R2
            IF (IR2.GT.T%IRCUT2) THEN
               !if (.true.) then
                  !CALL GET_R2(X, J1, J2, R2test, drtest)
                  !if (abs(r2-r2test) .gt.1d-6) then
                     !write(nl_myunit,*) "ERROR R2 incorrect"
                     !write(nl_myunit,*) "xoffset", sqrt(r2), sqrt(r2test)
                     !write(nl_myunit,*) "       ", j1, j2
                     !write(nl_myunit,*) X(3*(J1-1)+1:3*(J1-1)+3)
                     !write(nl_myunit,*) X(3*(J2-1)+1:3*(J2-1)+3)
                     !write(nl_myunit,*) X(3*(J1-1)+1:3*(J1-1)+3) - X(3*(J2-1)+1:3*(J2-1)+3)
                     !write(nl_myunit,*) xoffset
                  !endif
               !endif
               IR6=IR2**3
               E = VIJ(R2, IR6, T)
               !E = T%EPS*(T%SIG6*IR6*(T%SIG6*IR6-1.0D0) + T%RCONST/IR2 + T%CONST)
               POTEL = POTEL + E
            ENDIF
         END SUBROUTINE UPDATE_POTENTIAL_offset

         SUBROUTINE UPDATE_POTENTIAL_GRADIENT( X, J1, J2, POTEL, V, T )
            IMPLICIT NONE
            TYPE(INTERACTION_DEF), INTENT(IN) ::  T
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            DOUBLE PRECISION R2, E, XVEC(3), G, IR2, IR6, IR8, IR14
            INTEGER J5
            CALL GET_R2(X, J1, J2, R2, XVEC)
            IR2 = 1.D0/R2
            IF (IR2.GT.T%IRCUT2) THEN
               !update potential
               IR6=IR2**3
               E = VIJ(R2, IR6, T)
               !E = T%EPS*(T%SIG6*IR6*(T%SIG6*IR6-1.0D0) + T%RCONST/IR2 + T%CONST)
               POTEL = POTEL + E
               !update gradient
               IR8=IR6*IR2
               IR14=IR8*IR6
               G = DVIJ(IR8, IR14, T)
               !G = -8.0D0*T%EPS*(3.0D0*(2.0D0*IR14*(T%SIG12)-IR8*T%SIG6)-T%RCONST)
               DO J5=1,3
                  V(3*(J1-1)+J5)=V(3*(J1-1)+J5)+G*XVEC(J5)
                  V(3*(J2-1)+J5)=V(3*(J2-1)+J5)-G*XVEC(J5)
               END DO
            ENDIF
         END SUBROUTINE UPDATE_POTENTIAL_GRADIENT

         SUBROUTINE SET_INTERACTION_DEF( T, EPS, SIG, RC )
            IMPLICIT NONE
            !careful, with INTENT(OUT) here, any element of T not assigned in this subroutine could become undefined
            TYPE(INTERACTION_DEF), INTENT(OUT) ::  T
            DOUBLE PRECISION, INTENT(IN) ::  eps, SIG, RC
            DOUBLE PRECISION SIGRC6, SIGRC12
            T%EPS = EPS
            T%SIG = SIG
            T%RCUT = RC * T%SIG
            T%IRCUT2 = 1.D0/T%RCUT**2
            T%SIG6 = T%SIG**6
            T%SIG12 = T%SIG6**2

            SIGRC6 = T%SIG6/T%RCUT**6
            SIGRC12 = SIGRC6**2

            T%CONST = 4.0D0*(SIGRC6)-7.0D0*SIGRC12
            T%RCONST = (6.0D0*SIGRC12-3.0D0*SIGRC6)/T%RCUT**2
         END SUBROUTINE SET_INTERACTION_DEF

         !Above this will be only general subroutines and functions.
         !Below this will be subroutines which depend on the specific type of
         !interaction, AA, AB, BB.

         SUBROUTINE LJPSHIFT_GET_CUT( CUTAA, CUTBB, CUTAB )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(OUT) :: CUTAA, CUTBB, CUTAB 
            CUTAA = AA%RCUT
            CUTBB = BB%RCUT
            CUTAB = AB%RCUT
         END SUBROUTINE LJPSHIFT_GET_CUT

         SUBROUTINE LJPSHIFT_CLASS_SETUP(CUTOFF, EPSAA, EPSBB, EPSAB, SIGAA, SIGBB, SIGAB, NATOMS, BOXLX_I, BOXLY_I, BOXLZ_I)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: CUTOFF, EPSAA, EPSBB, EPSAB, SIGAA, SIGBB, SIGAB
            DOUBLE PRECISION, INTENT(IN) :: BOXLX_I, BOXLY_I, BOXLZ_I
            LOGICAL :: FIRST_THROUGH = .TRUE. !will be defined as .TRUE. on first initialization, but will retain value on subsequent calls
            INTEGER, INTENT(IN) :: NATOMS

            IF ( .NOT. FIRST_THROUGH ) RETURN
            FIRST_THROUGH = .FALSE.

            N = NATOMS
            BOXLX = BOXLX_I
            BOXLY = BOXLY_I
            BOXLZ = BOXLZ_I
            IBOXLX = 1.d0/BOXLX_I
            IBOXLY = 1.d0/BOXLY_I
            IBOXLZ = 1.d0/BOXLZ_I
            BOXLVEC(1) = BOXLX
            BOXLVEC(2) = BOXLY
            BOXLVEC(3) = BOXLZ
            IBOXLVEC(:) = 1.d0 / BOXLVEC(:)
            BOXLMIN_HALF_2 = (MINVAL(BOXLVEC)/2.D0)**2

            CALL SET_INTERACTION_DEF(AA, EPSAA, SIGAA, CUTOFF)
            CALL SET_INTERACTION_DEF(BB, EPSBB, SIGBB, CUTOFF)
            CALL SET_INTERACTION_DEF(AB, EPSAB, SIGAB, CUTOFF)

         END SUBROUTINE LJPSHIFT_CLASS_SETUP

         !3 subroutines for calculating the potentials
         SUBROUTINE LJPSHIFT_UPDATE_E_AA( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, AA)
         END SUBROUTINE LJPSHIFT_UPDATE_E_AA
         SUBROUTINE LJPSHIFT_UPDATE_E_BB( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, BB)
         END SUBROUTINE LJPSHIFT_UPDATE_E_BB
         SUBROUTINE LJPSHIFT_UPDATE_E_AB( X, J1, J2, POTEL )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL(X, J1, J2, POTEL, AB)
         END SUBROUTINE LJPSHIFT_UPDATE_E_AB

         !3 subroutines for calculating the potential and gradients
         SUBROUTINE LJPSHIFT_UPDATE_EG_AA( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, AA)
         END SUBROUTINE LJPSHIFT_UPDATE_EG_AA
         SUBROUTINE LJPSHIFT_UPDATE_EG_BB( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, BB)
         END SUBROUTINE LJPSHIFT_UPDATE_EG_BB
         SUBROUTINE LJPSHIFT_UPDATE_EG_AB( X, J1, J2, POTEL, V )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL, V(3*N)
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_GRADIENT(X, J1, J2, POTEL, V, AB)
         END SUBROUTINE LJPSHIFT_UPDATE_EG_AB

         !3 subroutines for calculating the potentials with offset
         SUBROUTINE LJPSHIFT_UPDATE_E_AA_offset( X, J1, J2, POTEL, XOFFSET )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N), XOFFSET(3)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_OFFSET(X, J1, J2, POTEL, AA, XOFFSET)
         END SUBROUTINE LJPSHIFT_UPDATE_E_AA_offset
         SUBROUTINE LJPSHIFT_UPDATE_E_BB_offset( X, J1, J2, POTEL, XOFFSET )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N), XOFFSET(3)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_OFFSET(X, J1, J2, POTEL, BB, XOFFSET)
         END SUBROUTINE LJPSHIFT_UPDATE_E_BB_offset
         SUBROUTINE LJPSHIFT_UPDATE_E_AB_offset( X, J1, J2, POTEL, XOFFSET )
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: X(3*N), XOFFSET(3)
            DOUBLE PRECISION, INTENT(OUT) :: POTEL
            INTEGER, INTENT(IN) :: J1, J2
            CALL UPDATE_POTENTIAL_OFFSET(X, J1, J2, POTEL, AB, XOFFSET)
         END SUBROUTINE LJPSHIFT_UPDATE_E_AB_offset

      END MODULE LJPSHIFT_CLASS

!*************************************************************************
!*************************************************************************

MODULE FREEZE_NL_MOD
   ! REMOVE THE FROZEN-FROZEN INTERACTIONS FROM THE neighbro LISTS.
   ! Calculate the energy for all frozen-frozen interactions and record it in EOFFSET
   IMPLICIT NONE
   PRIVATE
   INTEGER, ALLOCATABLE :: FREEZE_NL_AALIST(:,:)
   INTEGER, ALLOCATABLE :: FREEZE_NL_BBLIST(:,:)
   INTEGER, ALLOCATABLE :: FREEZE_NL_ABLIST(:,:)
   INTEGER FREEZE_NL_NAA, FREEZE_NL_NBB, FREEZE_NL_NAB
   INTEGER NATOMS, NTYPEA, NTYPEB
   DOUBLE PRECISION FREEZE_NL_EOFFSET

   LOGICAL FREEZE
   LOGICAL, ALLOCATABLE :: FROZEN(:)

   PUBLIC :: FREEZE_NL_SETUP, FREEZE_NL_UPDATE, FREEZE_NL_AALIST, FREEZE_NL_BBLIST, &
      FREEZE_NL_ABLIST, FREEZE_NL_NAA, FREEZE_NL_NBB, FREEZE_NL_NAB, FREEZE_NL_EOFFSET

   CONTAINS

   SUBROUTINE FREEZE_NL_SETUP( NATOMS_I, NTYPEA_I, FREEZE_I, FROZEN_I)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS_I, NTYPEA_I
      LOGICAL, INTENT(IN) :: FROZEN_I(NATOMS_I), FREEZE_I
      LOGICAL :: FIRST = .TRUE.

      IF (.NOT. FIRST) RETURN
      FIRST = .FALSE.

      NATOMS = NATOMS_I
      NTYPEA = NTYPEA_I
      NTYPEB = NATOMS - NTYPEA
      ALLOCATE( FREEZE_NL_AALIST(2,NTYPEA * (NTYPEA-1)/2) )
      ALLOCATE( FREEZE_NL_BBLIST(2,NTYPEB * (NTYPEB-1)/2) )
      ALLOCATE( FREEZE_NL_ABLIST(2,NTYPEA * NTYPEB) )

      FREEZE = FREEZE_I
      IF (FREEZE) THEN
         ALLOCATE( FROZEN(NATOMS) )
         FROZEN(:) = FROZEN_I(:)
      ENDIF
   END SUBROUTINE FREEZE_NL_SETUP

   SUBROUTINE FREEZE_NL_UPDATE ( COORDS, AALIST_I, NAA_I, BBLIST_I, NBB_I, ABLIST_I, NAB_I )
      !MAKE NEW LISTS FROM THE OLD ONES EXCLUDING FROZEN-FROZEN INTERACTIONS.  
      !CALCULATE THE ENERGY FROM THE FROZEN-FROZEN INTERACTIONS AND SAVE THE SUM
      !IN EOFFSET
      USE LJPSHIFT_CLASS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NAA_I, NBB_I, NAB_I
      INTEGER, INTENT(IN) :: AALIST_I(2,NAA_I), BBLIST_I(2,NBB_I), ABLIST_I(2,NAB_I)
      DOUBLE PRECISION, INTENT(IN) :: COORDS(NATOMS*3)
      INTEGER J1, J2, J3

      FREEZE_NL_EOFFSET = 0.D0
      FREEZE_NL_NAA = 0
      FREEZE_NL_NBB = 0
      FREEZE_NL_NAB = 0
      DO J3=1,NAA_I
         J1=AALIST_I(1,J3)
         J2=AALIST_I(2,J3)
         IF (FROZEN(J1) .AND. FROZEN(J2)) THEN
            CALL LJPSHIFT_UPDATE_E_AA(COORDS, J1, J2, FREEZE_NL_EOFFSET)
         ELSE
            FREEZE_NL_NAA = FREEZE_NL_NAA + 1
            FREEZE_NL_AALIST(1,FREEZE_NL_NAA) = J1
            FREEZE_NL_AALIST(2,FREEZE_NL_NAA) = J2
         ENDIF
      ENDDO
      DO J3=1,NBB_I
         J1=BBLIST_I(1,J3)
         J2=BBLIST_I(2,J3)
         IF (FROZEN(J1) .AND. FROZEN(J2)) THEN
            CALL LJPSHIFT_UPDATE_E_BB(COORDS, J1, J2, FREEZE_NL_EOFFSET)
            !PRINT *, J1, J2, FROZEN(J1), FROZEN(J2)
         ELSE
            !PRINT *, "UNFROZEN", J1, J2, FROZEN(J1), FROZEN(J2)
            FREEZE_NL_NBB = FREEZE_NL_NBB + 1
            FREEZE_NL_BBLIST(1,FREEZE_NL_NBB) = J1
            FREEZE_NL_BBLIST(2,FREEZE_NL_NBB) = J2
         ENDIF
      ENDDO
      DO J3=1,NAB_I
         J1=ABLIST_I(1,J3)
         J2=ABLIST_I(2,J3)
         IF (FROZEN(J1) .AND. FROZEN(J2)) THEN
            CALL LJPSHIFT_UPDATE_E_AB(COORDS, J1, J2, FREEZE_NL_EOFFSET)
         ELSE
            FREEZE_NL_NAB = FREEZE_NL_NAB + 1
            FREEZE_NL_ABLIST(1,FREEZE_NL_NAB) = J1
            FREEZE_NL_ABLIST(2,FREEZE_NL_NAB) = J2
         ENDIF
      ENDDO
      !WRITE(*,*) "FIL UPDATE", NAA_I, NBB_I, NAB_I
      !WRITE(*,*) "          ", FREEZE_NL_NAA, FREEZE_NL_NBB, FREEZE_NL_NAB, FREEZE_NL_EOFFSET

   END SUBROUTINE FREEZE_NL_UPDATE 

END MODULE FREEZE_NL_MOD

!*************************************************************************
!*************************************************************************

      SUBROUTINE LJPSHIFT(X, V, POTEL, GTEST, STEST)
      !This subroutine calculates the binary lennard jones potential with cutoff
      !
      !This is essentially a wrapper function.  The actual work is done in one
      !of the other LJPSHIFT subroutines
      !
      USE COMMONS, ONLY : NATOMS, CUTOFF, FIXIMAGE, NORESET, BOXLX, BOXLY, BOXLZ, &
     &    FREEZE, RESTRICTREGION, EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA, &
     &    ONE_ATOM_TAKESTEP
      USE LJPSHIFT_CLASS
      USE FREEZE_NL_MOD
      USE NEIGHBOR_LIST_MOD
      USE BIN_NL_MOD
      USE CELL_LISTS_BINARY_MOD
      IMPLICIT NONE
      INTEGER J1, J2
      DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
      LOGICAL, INTENT(IN) :: GTEST, STEST
      LOGICAL IL_CHANGED
      LOGICAL, SAVE :: USE_NEIGHBOR_LISTS = .TRUE. !this should be passable
      LOGICAL :: FIRST = .TRUE.
      integer, save :: num_calls = 0
      double precision :: oldpotel
      !double precision :: time0, time1, time2, time01=0.d0, time12=0.d0, time02=0.d0
      num_calls = num_calls + 1

      !CALL MYCPU_TIME(TIME0)
      IF (FIRST) THEN
         CALL LJPSHIFT_CLASS_SETUP( CUTOFF, 1.D0, EPSBB, EPSAB, 1.D0, SIGBB, SIGAB, NATOMS, BOXLX, BOXLY, BOXLZ)
      ENDIF

      FIRST = .FALSE.

      !CALL MYCPU_TIME(TIME1)
!
!  Calculate interatomic vectors using the minimum image convention.
!
      POTEL=0.0D0
      IF (GTEST .OR. STEST) THEN
         V(1:3*NATOMS)=0.D0
      ENDIF

      IF (STEST) THEN
        write(*,*) "warning: calculation of the hessian is not implimented in GMIN"
        ! js850> it could easily be implimented, but there is no way to return
        ! the calculated matrix.  OPTIM uses a module which is not implimented
        ! (and not needed) in GMIN
      ENDIF


      !IF (ONE_ATOM_TAKESTEP .AND. ONE_ATOM_USE .AND. .NOT. GTEST )  THEN
      IF (ONE_ATOM_TAKESTEP .AND. .NOT. GTEST )  THEN
         CALL LJPSHIFT_ONE_ATOM2(X, V, POTEL, .false.)
         IF (.false. .and. mod(num_calls,100000).eq.1) THEN
            OLDPOTEL = 0.D0
            CALL LJPSHIFT_NEIGHBOR_LIST( X, V, oldpotel, .FALSE., .FALSE.)
            WRITE(*,'(A,3G27.12)') "ljpshift> potel: ", OLDPOTEL-POTEL, POTEL, OLDPOTEL
         ENDIF
         RETURN
      ENDIF

!
!  Deal with any atoms that have left the box.
!
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET).and.(.not.ONE_ATOM_TAKESTEP)) THEN
         DO J1=1,NATOMS
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF

      IF (USE_NEIGHBOR_LISTS) THEN
         CALL LJPSHIFT_NEIGHBOR_LIST( X, V, POTEL, GTEST, STEST)
         return 
      ENDIF



      IF (FREEZE .AND. RESTRICTREGION) THEN
         CALL LJPSHIFT_INTERACTION_LIST (X, V, POTEL, GTEST, STEST)
      ELSEIF ((GTEST .OR. STEST) .AND. FREEZE) THEN
         CALL LJPSHIFT_FROZENLIST (X, V, POTEL, GTEST, STEST)
      ELSEIF ( GTEST .OR. STEST ) THEN
        !update POTEL and V
        DO J1=1,NTYPEA
          DO J2=J1+1,NTYPEA
            CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=1,NTYPEA
          DO J2=NTYPEA+1,NATOMS
            CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
        DO J1=NTYPEA+1,NATOMS
          DO J2=J1+1,NATOMS
            CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
          ENDDO
        ENDDO
      ELSE
         !only update POTEL
         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_E_AA(X, J1, J2, POTEL)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1, NATOMS
               CALL LJPSHIFT_UPDATE_E_AB(X, J1, J2, POTEL)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, NATOMS
            DO J2=J1+1, NATOMS
               CALL LJPSHIFT_UPDATE_E_BB(X, J1, J2, POTEL)
            ENDDO
         ENDDO
      ENDIF


      !CALL MYCPU_TIME(TIME2)
      !Time01 = time01 + TIME1 - TIME0
      !Time12 = time12 + TIME2 - TIME1
      !Time02 = time02 + TIME2 - TIME0
      !write(*,*) "ljpshift> times", time01, time12, time02
      RETURN
      END SUBROUTINE LJPSHIFT

!****************************************************************************
!****************************************************************************

      SUBROUTINE LJPSHIFT_INTERACTION_LIST (X, V, POTEL, GTEST, STEST)
      !This subroutine calculates the potential using interaction lists.  The
      !interaction lists are used so that atom pairs with unchanging
      !interactions are calculated only once.  It is assumed that there are frozen
      !particles and that restrictregion is set
      !
      !This subroutine is outdated.  It's probably better to use
      !LJPSHIFT_NEIGHBOR_LIST
      USE COMMONS, ONLY : NATOMS, &
     &    FROZEN, &
     &    FREEZEIL, FREEZEIL_USE, FREEZEIL_NAA, FREEZEIL_NAB, FREEZEIL_NBB, FREEZEIL_E, &
     &    RESTRICTREGIONRADIUS, RESTRICTREGIONX0, RESTRICTREGIONY0, RESTRICTREGIONZ0, &
     &    NTYPEA
      USE LJPSHIFT_CLASS
      IMPLICIT NONE
      INTEGER J1, J2, J3
      DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
      DOUBLE PRECISION  CUTAA, CUTBB, CUTAB
      LOGICAL, INTENT(IN) :: GTEST, STEST
      INTEGER AASTART, ABSTART, BBSTART
      DOUBLE PRECISION RRX, RRY, RRZ, RRR1, RRR2, CUT
      LOGICAL USEIL
      !choose AASTART, ABSTART, BBSTART so that all interactions can fit in
      !the array FREEZEIL.
      AASTART = 1;
      BBSTART = AASTART+NTYPEA*(NTYPEA-1);
      ABSTART = BBSTART+(NATOMS-NTYPEA)*(NATOMS-NTYPEA-1);

      !write(*,*) "using interactionlist"
      CALL LJPSHIFT_GET_CUT( CUTAA, CUTBB, CUTAB )

      !IF FREEZEIL_USE is false then we must set up the interactionlist
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
               IF (RRR1.GT.RESTRICTREGIONRADIUS+CUT) USEIL=.FALSE.
               IF (RRR2.GT.RESTRICTREGIONRADIUS+CUT) USEIL=.FALSE.
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
                  ELSE
                     FREEZEIL( ABSTART + FREEZEIL_NAB ,1) = J1
                     FREEZEIL( ABSTART + FREEZEIL_NAB ,2) = J2
                     FREEZEIL_NAB = FREEZEIL_NAB + 1
                  ENDIF

               ELSE
                  !calc energy of j1,j2 interaction and add to FREEZEIL_E
                  IF (J1.LE.NTYPEA .AND. J2.LE.NTYPEA) THEN
                     CALL LJPSHIFT_UPDATE_E_AA(X, J1, J2, FREEZEIL_E)
                  ELSE IF (J1.GT.NTYPEA .AND. J2.GT.NTYPEA) THEN
                     CALL LJPSHIFT_UPDATE_E_BB(X, J1, J2, FREEZEIL_E)
                  ELSE
                     CALL LJPSHIFT_UPDATE_E_AB(X, J1, J2, FREEZEIL_E)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      !write(*,*) FREEZEIL_NAA, FREEZEIL_NBB, FREEZEIL_NAB, FREEZEIL_NAA+ FREEZEIL_NBB+ FREEZEIL_NAB
      !
      !calculate the potential using the interaction list
      !
      IF ((GTEST .OR. STEST)) THEN
         !calculate both potential and gradient.
         DO J3=AASTART,AASTART+FREEZEIL_NAA-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            !WRITE (*,*) "AA", J1, J2, FROZEN(J1), FROZEN(J2)
            CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
         ENDDO
         DO J3=BBSTART,BBSTART+FREEZEIL_NBB-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
         ENDDO
         DO J3=ABSTART,ABSTART+FREEZEIL_NAB-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
         ENDDO
      ELSE
         !calculate only the potential
         DO J3=AASTART,AASTART+FREEZEIL_NAA-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            CALL LJPSHIFT_UPDATE_E_AA(X, J1, J2, POTEL)
         ENDDO
         DO J3=BBSTART,BBSTART+FREEZEIL_NBB-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            CALL LJPSHIFT_UPDATE_E_BB(X, J1, J2, POTEL)
         ENDDO
         DO J3=ABSTART,ABSTART+FREEZEIL_NAB-1
            J1=FREEZEIL(J3,1)
            J2=FREEZEIL(J3,2)
            CALL LJPSHIFT_UPDATE_E_AB(X, J1, J2, POTEL)
         ENDDO
      ENDIF
      POTEL = POTEL + FREEZEIL_E
      END SUBROUTINE LJPSHIFT_INTERACTION_LIST

!****************************************************************************
!****************************************************************************

      SUBROUTINE LJPSHIFT_FROZENLIST (X, V, POTEL, GTEST, STEST)
      !This subroutine calculates the ljpshift potential using FROZENLIST to
      !ensure that interactions between frozen atoms are calculated only once.
      !
      !This subroutine is outdated.  It's probably better to use
      !LJPSHIFT_NEIGHBOR_LIST
      !
      USE COMMONS, ONLY : NATOMS, &
     &    FREEZESAVE, FREEZESAVEE, NFREEZE, NFREEZETYPEA, FROZEN, FROZENLIST, NTYPEA
      USE LJPSHIFT_CLASS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NNOFREEZETYPEA
      DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
      DOUBLE PRECISION POTELDUM
      LOGICAL, INTENT(IN) :: GTEST, STEST

      NNOFREEZETYPEA = NTYPEA - NFREEZETYPEA

      IF ( .NOT. FREEZESAVE ) THEN

         !calculate interactions between mobile particles
         !AA
         DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
            J1=FROZENLIST(J3)
            DO J4=J3+1,NFREEZE+NNOFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !BB
         DO J3=NFREEZE+NNOFREEZETYPEA+1,NATOMS
            J1=FROZENLIST(J3)
            DO J4=J3+1,NATOMS
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !AB
         DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
            J1=FROZENLIST(J3)
            DO J4=NFREEZE+NNOFREEZETYPEA+1,NATOMS
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !calculate interactions between mobile and imobile atoms
         !AA
         DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
            J1=FROZENLIST(J3)
            DO J4=1,NFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !BB
         DO J3=NFREEZE+NNOFREEZETYPEA+1,NATOMS
            J1=FROZENLIST(J3)
            DO J4=NFREEZETYPEA+1,NFREEZE
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !AB mobile-imobile
         DO J3=NFREEZE+1,NFREEZE+NNOFREEZETYPEA
            J1=FROZENLIST(J3)
            DO J4=NFREEZETYPEA+1,NFREEZE
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO
         !BA mobile-imobile
         DO J3=NFREEZE+NNOFREEZETYPEA+1,NATOMS
            J1=FROZENLIST(J3)
            DO J4=1,NFREEZETYPEA
               J2=FROZENLIST(J4)
               CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
            ENDDO
         ENDDO

         POTEL=POTEL+FREEZESAVEE

      ELSE

         !The first time through we must calculate the contribution to the
         !energy from the frozen frozen interactions.  This energy will be
         !stored in FREEZESAVEE for later use
         FREEZESAVEE = 0.D0
         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               POTELDUM=0.D0
               CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,NATOMS
               POTELDUM=0.D0
               CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,NATOMS
            DO J2=J1+1,NATOMS
               POTELDUM=0.D0
               CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTELDUM, V)
               POTEL=POTEL+POTELDUM
               IF ( FROZEN(J1) .AND. FROZEN(J2) ) FREEZESAVEE=FREEZESAVEE+POTELDUM
            ENDDO
         ENDDO

         FREEZESAVE = .FALSE.

      ENDIF


      END SUBROUTINE LJPSHIFT_FROZENLIST 

SUBROUTINE LJPSHIFT_NEIGHBOR_LIST(X, V, POTEL, GTEST, STEST)
   !This subroutine calculates the potential using neighbor lists
   !The neighbor lists are created and updated using modules FREEZE_NL_MOD,
   !NEIGHBOR_LIST_MOD, and BIN_NL_MOD
   !
   ! NEIGHBOR_LISTS_MOD is the top level module, it creates and maintains the
   !     neighborlists.  The neighbor list is only updated if at least one atom
   !     has moved more than a certain amount (related to the skin depth).
   !
   ! BIN_NL_MOD converts the full neighbor list into 3 neighbor lists, one for
   !     each interaction type
   !
   ! FREEZE_NL_MOD removes the frozen-frozen interactions from the 3 neighbor
   !     lists of BIN_NL_MOD.  It also records the energy of the frozen-frozen
   !     interactions in FREEZE_NL_EOFFSET
   !
   USE COMMONS, ONLY : NATOMS, CUTOFF, BOXLX, BOXLY, BOXLZ, &
   &    NTYPEA, FREEZE, FROZEN
   USE FREEZE_NL_MOD
   USE NEIGHBOR_LIST_MOD
   USE BIN_NL_MOD
   IMPLICIT NONE
   INTEGER J1, J2
   DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   LOGICAL, INTENT(IN) :: GTEST, STEST
   LOGICAL IL_CHANGED
   LOGICAL, SAVE :: FIRST = .TRUE.

   !CALL MYCPU_TIME(TIME0)
   IF (FIRST) THEN
      CALL NL_SETUP( NATOMS, X, CUTOFF, BOXLX, BOXLY, BOXLZ )
      CALL BIN_NL_SETUP( NATOMS, NTYPEA )
      IF (FREEZE) CALL FREEZE_NL_SETUP( NATOMS, NTYPEA, FREEZE, FROZEN )
   ENDIF

   CALL NL_UPDATE( X, IL_CHANGED )

   IF ( (IL_CHANGED .OR. FIRST)) THEN
      CALL BIN_NL_UPDATE ( NL_LIST, NL_NLIST )
      if (FREEZE) CALL FREEZE_NL_UPDATE ( X, BIN_NL_AALIST, BIN_NL_NAA, BIN_NL_BBLIST, BIN_NL_NBB, BIN_NL_ABLIST, BIN_NL_NAB )
   ENDIF

   FIRST = .FALSE.

   IF (FREEZE) THEN
      CALL LJPSHIFT_INTERACTION_LIST2 (X, V, POTEL, GTEST, STEST, NATOMS, FREEZE_NL_AALIST, &
         FREEZE_NL_NAA, FREEZE_NL_BBLIST, FREEZE_NL_NBB, FREEZE_NL_ABLIST, FREEZE_NL_NAB, FREEZE_NL_EOFFSET)
   ELSE
      CALL LJPSHIFT_INTERACTION_LIST2 (X, V, POTEL, GTEST, STEST, NATOMS, BIN_NL_AALIST, &
         BIN_NL_NAA, BIN_NL_BBLIST, BIN_NL_NBB, BIN_NL_ABLIST, BIN_NL_NAB, 0.D0)
   ENDIF
END SUBROUTINE LJPSHIFT_NEIGHBOR_LIST

SUBROUTINE LJPSHIFT_INTERACTION_LIST2 (X, V, POTEL, GTEST, STEST, NATOMS, AALIST, &
NAA, BBLIST, NBB, ABLIST, NAB, EOFFSET)
   !THIS SUBROUTINE CALCULATES THE POTENTIAL USING INTERACTION LISTS.
   !Calculate the potential energy of all atom pairs in the neighbor lists
   !AALIST, BBLIST, ABLIST
   USE LJPSHIFT_CLASS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NAA, NBB, NAB, NATOMS
   DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   DOUBLE PRECISION, INTENT(IN) :: EOFFSET
   INTEGER, INTENT(IN) :: AALIST(2,NAA), BBLIST(2,NBB), ABLIST(2,NAB)
   LOGICAL, INTENT(IN) :: GTEST, STEST
   INTEGER J1, J2, J3

   !
   !CALCULATE THE POTENTIAL USING THE INTERACTION LIST
   !

   IF ((GTEST .OR. STEST)) THEN
      !CALCULATE BOTH POTENTIAL AND GRADIENT.
      DO J3=1,NAA
         J1=AALIST(1,J3)
         J2=AALIST(2,J3)
         !WRITE (*,*) "AA", J1, J2, FROZEN(J1), FROZEN(J2)
         CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
      ENDDO
      DO J3=1,NBB
         J1=BBLIST(1,J3)
         J2=BBLIST(2,J3)
         CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
      ENDDO
      DO J3=1,NAB
         J1=ABLIST(1,J3)
         J2=ABLIST(2,J3)
         CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
      ENDDO
   ELSE
      !CALCULATE ONLY THE POTENTIAL
      DO J3=1,NAA
         J1=AALIST(1,J3)
         J2=AALIST(2,J3)
         CALL LJPSHIFT_UPDATE_E_AA(X, J1, J2, POTEL)
      ENDDO
      DO J3=1,NBB
         J1=BBLIST(1,J3)
         J2=BBLIST(2,J3)
         CALL LJPSHIFT_UPDATE_E_BB(X, J1, J2, POTEL)
      ENDDO
      DO J3=1,NAB
         J1=ABLIST(1,J3)
         J2=ABLIST(2,J3)
         CALL LJPSHIFT_UPDATE_E_AB(X, J1, J2, POTEL)
      ENDDO
   ENDIF
   POTEL = POTEL + EOFFSET
END SUBROUTINE LJPSHIFT_INTERACTION_LIST2

SUBROUTINE LJPSHIFT_INTERACTION_LIST3 (X, V, POTEL, GTEST, STEST, NATOMS, AALIST, &
   NAA, BBLIST, NBB, ABLIST, NAB, EOFFSET, ATOMI)
   !THIS SUBROUTINE CALCULATES THE POTENTIAL USING INTERACTION LISTS.
   !Calculate the potential energy of atomi with all the atoms in the lists
   !AALIST, BBLIST, ABLIST
   USE LJPSHIFT_CLASS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NAA, NBB, NAB, NATOMS, ATOMI
   DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   DOUBLE PRECISION, INTENT(IN) :: EOFFSET
   INTEGER, INTENT(IN) :: AALIST(NAA), BBLIST(NBB), ABLIST(NAB)
   LOGICAL, INTENT(IN) :: GTEST, STEST
   INTEGER J1, J2, J3

   !write(*,*) "LJPSHIFT_INTERACTION_LIST3>", atomi, naa, nbb, nab

   !
   !CALCULATE THE POTENTIAL USING THE INTERACTION LIST
   !

   J1 = ATOMI
   IF ((GTEST .OR. STEST)) THEN
      !CALCULATE BOTH POTENTIAL AND GRADIENT.
      DO J3=1,NAA
         J2=AALIST(J3)
         CALL LJPSHIFT_UPDATE_EG_AA(X, J1, J2, POTEL, V)
      ENDDO
      DO J3=1,NBB
         J2=BBLIST(J3)
         CALL LJPSHIFT_UPDATE_EG_BB(X, J1, J2, POTEL, V)
      ENDDO
      DO J3=1,NAB
         J2=ABLIST(J3)
         CALL LJPSHIFT_UPDATE_EG_AB(X, J1, J2, POTEL, V)
      ENDDO
   ELSE
      !CALCULATE ONLY THE POTENTIAL
      DO J3=1,NAA
         J2=AALIST(J3)
         CALL LJPSHIFT_UPDATE_E_AA(X, J1, J2, POTEL)
      ENDDO
      DO J3=1,NBB
         J2=BBLIST(J3)
         CALL LJPSHIFT_UPDATE_E_BB(X, J1, J2, POTEL)
      ENDDO
      DO J3=1,NAB
         J2=ABLIST(J3)
         CALL LJPSHIFT_UPDATE_E_AB(X, J1, J2, POTEL)
      ENDDO
   ENDIF
   POTEL = POTEL + EOFFSET
END SUBROUTINE LJPSHIFT_INTERACTION_LIST3

SUBROUTINE LJPSHIFT_INTERACTION_LIST4 (X, V, POTEL, GTEST, STEST, NATOMS, list, &
   nlist)
   !THIS SUBROUTINE CALCULATES THE POTENTIAL between all the atoms in atomlist
   !this is not a very fast method for many atoms.  it is intended to be used
   !for very short lists
   USE LJPSHIFT_CLASS
   use commons, only : ntypea
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NATOMS, nlist
   DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   INTEGER, INTENT(IN) :: list(nlist)
   LOGICAL, INTENT(IN) :: GTEST, STEST
   INTEGER J1, J2, J3, j4

   !
   !CALCULATE THE POTENTIAL USING THE INTERACTION LIST
   !

   IF ((GTEST .OR. STEST)) THEN
      !CALCULATE BOTH POTENTIAL AND GRADIENT.
      DO J1=1,NLIST
         J3 = LIST(J1)
         DO J2=1,J1-1
            J4 = LIST(J2)
            !AA
            IF ((J3.LE.NTYPEA) .AND. (J4.LE.NTYPEA)) THEN
               CALL LJPSHIFT_UPDATE_EG_AA(X, J3, J4, POTEL, V)
            ELSEIF ((J3.GT.NTYPEA) .AND. (J4.GT.NTYPEA)) THEN
               CALL LJPSHIFT_UPDATE_EG_BB(X, J3, J4, POTEL, V)
            ELSE
               CALL LJPSHIFT_UPDATE_EG_AB(X, J3, J4, POTEL, V)
            ENDIF
         ENDDO
      ENDDO
   ELSE
      !CALCULATE ONLY THE POTENTIAL
      DO J1=1,NLIST
         J3 = LIST(J1)
         DO J2=1,J1-1
            J4 = LIST(J2)
            !write(*,*) "double count", j3, j4
            !AA
            IF ((J3.LE.NTYPEA) .AND. (J4.LE.NTYPEA)) THEN
               CALL LJPSHIFT_UPDATE_E_AA(X, J3, J4, POTEL)
            ELSEIF ((J3.GT.NTYPEA) .AND. (J4.GT.NTYPEA)) THEN
               CALL LJPSHIFT_UPDATE_E_BB(X, J3, J4, POTEL)
            ELSE
               CALL LJPSHIFT_UPDATE_E_AB(X, J3, J4, POTEL)
            ENDIF
         ENDDO
      ENDDO
   ENDIF
END SUBROUTINE LJPSHIFT_INTERACTION_LIST4

SUBROUTINE LJPSHIFT_INTERACTION_LIST5 (X, V, POTEL, GTEST, &
   STEST, NATOMS, &
   AALIST, NAA, xoffsetAA, &
   BBLIST, NBB, xoffsetBB, &
   ABLIST, NAB, xoffsetAB, &
   EOFFSET, ATOMI)
   !THIS SUBROUTINE CALCULATES THE POTENTIAL USING INTERACTION LISTS in the case
   !where the periodic offset has been predetermined
   !Calculate the potential energy of atomi with all the atoms in the lists
   !AALIST, BBLIST, ABLIST
   USE LJPSHIFT_CLASS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NAA, NBB, NAB, NATOMS, ATOMI
   DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   DOUBLE PRECISION, INTENT(IN) :: EOFFSET
   DOUBLE PRECISION, INTENT(IN) :: XOFFSETAA(3,NAA)
   DOUBLE PRECISION, INTENT(IN) :: XOFFSETBB(3,NAA)
   DOUBLE PRECISION, INTENT(IN) :: XOFFSETAB(3,NAA)
   INTEGER, INTENT(IN) :: AALIST(NAA), BBLIST(NBB), ABLIST(NAB)
   LOGICAL, INTENT(IN) :: GTEST, STEST
   INTEGER J1, J2, J3

   !write(*,*) "LJPSHIFT_INTERACTION_LIST5>", atomi, naa, nbb, nab

   !
   !CALCULATE THE POTENTIAL USING THE INTERACTION LIST
   !

   J1 = ATOMI
   IF ((GTEST .OR. STEST)) THEN
      !CALCULATE BOTH POTENTIAL AND GRADIENT.
      write(*,*) "error LJPSHIFT_INTERACTION_LIST5 not implemented with gradient"
   ELSE
      !CALCULATE ONLY THE POTENTIAL
      DO J3=1,NAA
         J2=AALIST(J3)
         CALL LJPSHIFT_UPDATE_E_AA_OFFSET(X, J1, J2, POTEL, XOFFSETAA(:,J3))
      ENDDO
      DO J3=1,NBB
         J2=BBLIST(J3)
         CALL LJPSHIFT_UPDATE_E_BB_OFFSET(X, J1, J2, POTEL, XOFFSETBB(:,J3))
      ENDDO
      DO J3=1,NAB
         J2=ABLIST(J3)
         CALL LJPSHIFT_UPDATE_E_AB_OFFSET(X, J1, J2, POTEL, XOFFSETAB(:,J3))
      ENDDO
   ENDIF
   POTEL = POTEL + EOFFSET
END SUBROUTINE LJPSHIFT_INTERACTION_LIST5

SUBROUTINE LJPSHIFT_ONE_ATOM_REBUILD(coords, v, energy)
   !assume all atoms have moved.  rebuild everything
   USE CELL_LISTS_BINARY_MOD
   USE COMMONS, ONLY : NATOMS
   IMPLICIT NONE
   DOUBLE PRECISION :: coords(3*NATOMS) , v(3*natoms)
   DOUBLE PRECISION, INTENT(out) :: energy
   integer i, j1, j2
   DOUBLE PRECISION IBOXLX, IBOXLY, IBOXLZ
   call cl_rebuild_cell_lists(coords)
   energy = 0.d0
   CALL LJPSHIFT_NEIGHBOR_LIST( COORDS, V, energy, .FALSE., .FALSE.)
END SUBROUTINE LJPSHIFT_ONE_ATOM_REBUILD

SUBROUTINE LJPSHIFT_ONE_ATOM(COORDS, V, POTEL, GTEST)
   !This subroutine calculates the binary lennard jones potential with cutoff
   !This subroutine is designed for potential updates for one atoms moves.  
   !
   !Use LJPSHIFT_SET_MOVED to communicate with this potential which atoms have
   !moved.
   !
   !If atom i has moved, the new potential energy is
   !
   ! Enew = Eold + sum_j ( Enew_{ij} - Eold_{ij} )
   !
   !Therefore we need to be able to calculate Eold_{ij} for each neighbor of i.
   !To do this, we save the a copy of coords in CL_oldcoords.  
   !
   !Cell lists are used to avoid calculation of atoms pairs that are very far
   !apart.  The module CELL_LISTS_BINARY_MOD is used for this purpose.  We also
   !use CELL_LISTS_BINARY_MOD to store things like the list of moved atoms, and
   !cl_oldcoords.  These perhaps should be moved out of that module
   USE COMMONS, ONLY : NATOMS, CUTOFF, BOXLX, BOXLY, BOXLZ, &
   &    EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
   USE CELL_LISTS_BINARY_MOD
   USE LJPSHIFT_CLASS, ONLY : LJPSHIFT_CLASS_SETUP
   IMPLICIT NONE
   INTEGER J1, J2
   DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   LOGICAL, INTENT(IN) :: GTEST
   LOGICAL IL_CHANGED
   LOGICAL, SAVE :: FIRST = .TRUE.
   integer :: i, imove, MYNODE, IERR, n
   DOUBLE PRECISION newenergy_i
   DOUBLE PRECISION, SAVE :: OLD_POTEL
   DOUBLE PRECISION :: emoved_old, emoved_new
   !double precision :: time0, time1, time2, time01=0.d0, time12=0.d0, time02=0.d0

   if (first) then
      !initialize the modules we use
      first = .false.
      CALL LJPSHIFT_CLASS_SETUP( CUTOFF, 1.D0, EPSBB, EPSAB, 1.D0, SIGBB, SIGAB, NATOMS, BOXLX, BOXLY, BOXLZ)
      CALL CELL_LISTS_BINARY_SETUP(COORDS, NATOMS, BOXLX, CUTOFF, NTYPEA)

      cl_nmoved = natoms 
   endif

   IF (2*CL_NMOVED .gt. natoms) THEN
      !if many particles have moved, then it's faster to just rebuild the whole
      !thing
      call LJPSHIFT_ONE_ATOM_REBUILD(COORDS, v, potel)
      !write(*,*) "initial potel ", sum(ONE_ATOM_ENERGIES)/2.d0
      cl_oldcoords(:) = coords(:)
      cl_nmoved = 0
      old_potel = potel
      !write(*,*) "energy from rebuild", potel
      return
   endif


!      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
!         do i=1,cl_nmoved
!            imove = cl_moved_atoms(i)
!            J2 = 3*(imove-1)
!            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
!            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
!            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
!         enddo
!      endif
   
   !
   !get the contribution to the old energy of all the moved particles
   !
   Emoved_old = 0.d0
   !write(*,*) "moved", cl_moved_atoms(1:cl_nmoved)
   call cl_get_neighbors_several(cl_oldcoords, cl_nmoved, cl_moved_atoms)
   CALL LJPSHIFT_INTERACTION_LIST2 (cl_oldCOORDS, V, Emoved_old, .false., .false., NATOMS, CL_LISTAA, &
         CL_NLISTAA, CL_LISTBB, CL_NLISTBB, CL_LISTAB, CL_NLISTAB, 0.d0)
   !write(*,"(A,F25.12,3I8)") "Emoved_old", Emoved_old, cl_nlistaa, cl_nlistbb, cl_nlistab


   !now update the cell lists for all the moved particles
   call CL_update_cell_lists(COORDS, cl_nmoved, cl_moved_atoms)

   !
   !now calculate the new energy contribution of the moved particles
   !
   Emoved_new = 0.d0
   call cl_get_neighbors_several(coords, cl_nmoved, cl_moved_atoms)
   CALL LJPSHIFT_INTERACTION_LIST2 (COORDS, V, Emoved_new, .false., .false., NATOMS, CL_LISTAA, &
         CL_NLISTAA, CL_LISTBB, CL_NLISTBB, CL_LISTAB, CL_NLISTAB, 0.d0)
   !write(*,"(A,F25.12,3I8)") "Emoved_new", Emoved_new, cl_nlistaa, cl_nlistbb, cl_nlistab
   !write(*,"(A,2F25.12,I9)") "Emoved_old ,new", Emoved_old, emoved_new, cl_nmoved

   !
   !update potential
   !
   potel = old_potel + (emoved_new - emoved_old)
   !write(*,*) "energy one atom", old_potel, potel


   !now update oldcoords and old_potel
   do n=1,cl_nmoved
      i = cl_moved_atoms(n)
      cl_oldcoords(3*(i-1)+1 : 3*(i-1)+3) = coords(3*(i-1)+1 : 3*(i-1)+3)
   enddo
   old_potel = potel




   !reset cl_nmoved
   cl_nmoved = 0

   

END SUBROUTINE LJPSHIFT_ONE_ATOM

SUBROUTINE LJPSHIFT_SET_MOVED(NMOVED, MOVED_ATOMS, MOVEALL)
   !This subroutine is used to communicate with the potential which atoms have
   !been moved.
   !It should be called from the monte carlo loop anytime an atom is moved.
   use commons, only : natoms
   use cell_lists_binary_mod
   implicit none
   integer, intent(in) :: nmoved, moved_atoms(natoms)
   logical, intent(in) :: moveall
   logical, save :: first = .true.
   integer :: i1, i2, j1, j2
   logical :: isnew
   !write(*,*) "in LJPSHIFT_SET_MOVED", nmoved
   if (first) then
      !wait until  the class has been fully set up
      first = .false.
      cl_nmoved = natoms
      return
   endif

   if (moveall .or. cl_nmoved.ge.natoms) then
      cl_nmoved = natoms
      return
   endif

   !WRITE(*,*) "nmoved", nmoved
   !don't assume cl_nmoved is 0
   !CL_NMOVED = NMOVED
   !CL_MOVED_ATOMS(1:NMOVED) = MOVED_ATOMS(1:NMOVED)

   !MAKE SURE THE NEW LIST WON'T BE TOO LONG
   IF (CL_NMOVED + NMOVED .GT. 10) THEN
      WRITE(*,*) "LJPSHIFT_SET_MOVED> WARNING: CL_NMOVED is very large", CL_NMOVED + NMOVED
      IF ((CL_NMOVED + NMOVED) .GT. NATOMS) THEN
         WRITE(*,*) "LJPSHIFT_SET_MOVED> ERROR: CL_NMOVED is greater than NATOMS"
         cl_nmoved = natoms
         return
      ENDIF
   ENDIF

   !add the new moved atoms to the list
   !removing duplicates
   do j1=1,nmoved
      i1 = moved_atoms(j1)
      !check if the new one is already on the list
      isnew = .true.
      do j2=1,cl_nmoved
         i2 = cl_moved_atoms(j2)
         if (i1 .eq. i2) then
            isnew = .false.
            exit !exit loop
         endif
      enddo
      if (isnew) then
         cl_nmoved = cl_nmoved + 1
         cl_moved_atoms(cl_nmoved) = i1
      endif
   enddo


   !eliminate duplicates
END SUBROUTINE LJPSHIFT_SET_MOVED

SUBROUTINE LJPSHIFT_SET_MOVED2(NMOVED, MOVED_ATOMS, MOVEALL, myunit)
   !This subroutine is used to communicate with the potential which atoms have
   !been moved.
   !It should be called from the monte carlo loop anytime an atom is moved.
   use commons, only : natoms
   use NL_BIN_MOVEONE, only : nl_bin_moveone_nmoved, nl_bin_moveone_moved_atoms, &
   nl_myunit
   implicit none
   integer, intent(in) :: nmoved, moved_atoms(natoms), myunit
   logical, intent(in) :: moveall
   logical, save :: first = .true.
   integer :: i1, i2, j1, j2
   logical :: isnew
   !write(myunit,*) "in LJPSHIFT_SET_MOVED", nl_bin_moveone_nmoved, nmoved, moveall
   if (first) then
      !wait until  the class has been fully set up
      first = .false.
      nl_bin_moveone_nmoved = natoms
      nl_myunit = myunit
      return
   endif

   !WRITE(myunit,*) "LJPSHIFT_SET_MOVED>", nmoved, moveall
   if (moveall .or. nmoved.ge.natoms) then
      !WRITE(myunit,*) "LJPSHIFT_SET_MOVED> moveall set", nl_bin_moveone_nmoved
      nl_bin_moveone_nmoved = natoms
      return
   endif

   !WRITE(*,*) "nmoved", nmoved
   !don't assume cl_nmoved is 0
   !CL_NMOVED = NMOVED
   !CL_MOVED_ATOMS(1:NMOVED) = MOVED_ATOMS(1:NMOVED)

   !MAKE SURE THE NEW LIST WON'T BE TOO LONG
   IF (nl_bin_moveone_NMOVED + NMOVED .GT. natoms/2) then
      nl_bin_moveone_nmoved = natoms
      return
   ENDIF
   !IF (nl_bin_moveone_NMOVED + NMOVED .GT. 20 .and. &
         !nl_bin_moveone_NMOVED.ne.natoms .and. nmoved .ne. natoms ) THEN
      !WRITE(myunit,*) "LJPSHIFT_SET_MOVED> WARNING: NMOVED is very large", nl_bin_moveone_NMOVED , NMOVED
      !IF ((nl_bin_moveone_NMOVED + NMOVED) .GT. NATOMS) THEN
         !WRITE(myunit,*) "LJPSHIFT_SET_MOVED> ERROR: NMOVED is greater than NATOMS"
         !nl_bin_moveone_nmoved = natoms
         !return
      !ENDIF
   !ENDIF

   !add the new moved atoms to the list
   !removing duplicates
   do j1=1,nmoved
      i1 = moved_atoms(j1)
      !check if the new one is already on the list
      isnew = .true.
      do j2=1,nl_bin_moveone_nmoved
         i2 = nl_bin_moveone_moved_atoms(j2)
         if (i1 .eq. i2) then
            isnew = .false.
            exit !exit loop
         endif
      enddo
      if (isnew) then
         nl_bin_moveone_nmoved = nl_bin_moveone_nmoved + 1
         nl_bin_moveone_moved_atoms(nl_bin_moveone_nmoved) = i1
      endif
   enddo


   !eliminate duplicates
END SUBROUTINE LJPSHIFT_SET_MOVED2

SUBROUTINE LJPSHIFT_ONE_ATOM_REBUILD2(coords, v, energy)
   !assume all atoms have moved.  rebuild everything
   USE NL_BIN_MOVEONE, ONLY : NL_BIN_MOVEONE_UPDATE_LISTS
   USE COMMONS, ONLY : NATOMS
   IMPLICIT NONE
   DOUBLE PRECISION :: coords(3*NATOMS) , v(3*natoms)
   DOUBLE PRECISION, INTENT(out) :: energy
   integer i, j1, j2, moved_atoms(natoms)
   DOUBLE PRECISION IBOXLX, IBOXLY, IBOXLZ
   !do j1=1,natoms
      !moved_atoms(j1) = j1
   !enddo
   call NL_BIN_MOVEONE_update_lists(coords, natoms, moved_atoms)
   energy = 0.d0
   CALL LJPSHIFT_NEIGHBOR_LIST( COORDS, V, energy, .FALSE., .FALSE.)
END SUBROUTINE LJPSHIFT_ONE_ATOM_REBUILD2

SUBROUTINE LJPSHIFT_ONE_ATOM2(COORDS, V, POTEL, GTEST)
   !This subroutine calculates the binary lennard jones potential with cutoff
   !This subroutine is designed for potential updates for one atoms moves.  
   !
   !Use LJPSHIFT_SET_MOVED2 to communicate with this potential which atoms have
   !moved.
   !
   !If atom i has moved, the new potential energy is
   !
   ! Enew = Eold + sum_j ( Enew_{ij} - Eold_{ij} )
   !
   !Therefore we need to be able to calculate Eold_{ij} for each neighbor of i.
   !To do this, we save the a copy of coords in oldcoords.  
   !
   !neighbor lists are used to avoid having to recompute which atoms are
   !neighbors.  periodic corrections are also saved for every neighbor pair
   !to avoid having to recalculate NINT(x / boxL).  Because of this the atoms
   !cannot be shifted to another periodic cell without rebuilding the lists.
   !i.e. you can't put the atoms back in the box
   USE COMMONS, ONLY : NATOMS, CUTOFF, BOXLX, BOXLY, BOXLZ, &
   &    EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
   USE NL_BIN_MOVEONE
   USE LJPSHIFT_CLASS, ONLY : LJPSHIFT_CLASS_SETUP
   IMPLICIT NONE
   INTEGER J1, J2
   DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATOMS) 
   DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POTEL 
   LOGICAL, INTENT(IN) :: GTEST
   LOGICAL IL_CHANGED
   LOGICAL, SAVE :: FIRST = .TRUE.
   integer :: i, imove, MYNODE, IERR, n
   integer nmoved, moved_atoms(natoms)
   DOUBLE PRECISION newenergy_i
   DOUBLE PRECISION, SAVE :: OLD_POTEL
   DOUBLE PRECISION :: emoved_old, emoved_new
   DOUBLE PRECISION :: edouble_count_old, edouble_count_new
   !double precision :: time0, time1, time2, time01=0.d0, time12=0.d0, time02=0.d0

   nmoved = nl_bin_moveone_NMOVED
   !if (nmoved .gt. 3) then
      !write(nl_myunit,*) "ljpshift> nmoved", nmoved
   !endif

   if (first) then
      !initialize the modules we use
      first = .false.
      CALL LJPSHIFT_CLASS_SETUP( CUTOFF, 1.D0, EPSBB, EPSAB, 1.D0, SIGBB, SIGAB, NATOMS, BOXLX, BOXLY, BOXLZ)
      CALL NL_BIN_MOVEONE_SETUP(COORDS, NATOMS, BOXLX, CUTOFF, NTYPEA)

      nmoved = natoms 
   endif

   IF (2*NMOVED .gt. natoms) THEN
      !if many particles have moved, then it's faster to just rebuild the whole
      !thing.  
      !Take this opportunity to put the atoms back in the box
      !
      !do j1=1,natoms
         !j2 = 3*(j1-1)
         !coords(j2+1) = coords(j2+1) - boxlx * nint(coords(j2+1) / boxlx)
         !coords(j2+2) = coords(j2+2) - boxly * nint(coords(j2+2) / boxly)
         !coords(j2+3) = coords(j2+3) - boxlz * nint(coords(j2+3) / boxlz)
      !enddo
      call LJPSHIFT_ONE_ATOM_REBUILD2(COORDS, v, potel)
      !write(*,*) "initial potel ", sum(ONE_ATOM_ENERGIES)/2.d0
      !nl_bin_moveone_oldcoords(:) = coords(:)
      nl_bin_moveone_nmoved = 0
      old_potel = potel
      !write(*,*) "energy from rebuild", potel
      return
   endif

   moved_atoms(1:nmoved) = nl_bin_moveone_moved_atoms(1:nmoved) !only for simplicity


!      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
!         do i=1,nl_bin_moveone_nmoved
!            imove = cl_moved_atoms(i)
!            J2 = 3*(imove-1)
!            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
!            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
!            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
!         enddo
!      endif
   
   !
   !get the contribution to the old energy of all the moved particles
   !
   Emoved_old = 0.d0
   !write(*,*) "moved", cl_moved_atoms(1:nl_bin_moveone_nmoved)
   do j1=1,nmoved
      j2 = moved_atoms(j1)
      !CALL LJPSHIFT_INTERACTION_LIST3 (nl_bin_moveone_oldCOORDS, V, &
         !Emoved_old, .false., .false., NATOMS,&
         !NL_BIN_MOVEONE_LISTSAA(:,J2), NLISTSAA(J2), &
         !NL_BIN_MOVEONE_LISTSBB(:,J2), NLISTSBB(J2), &
         !NL_BIN_MOVEONE_LISTSAB(:,J2), NLISTSAB(J2), &
         !0.D0, J2)
      CALL LJPSHIFT_INTERACTION_LIST5(NL_BIN_MOVEONE_OLDCOORDS, V, &
         EMOVED_OLD, .FALSE., .FALSE., NATOMS, &
         NL_BIN_MOVEONE_LISTSAA(:,J2), NLISTSAA(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSAA(:,:,j2), &
         NL_BIN_MOVEONE_LISTSBB(:,J2), NLISTSBB(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSBB(:,:,j2), &
         NL_BIN_MOVEONE_LISTSAB(:,J2), NLISTSAB(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSAB(:,:,j2), &
         0.D0, J2)
   enddo
   !
   !at this point we've counted the interactions between the atoms
   !in moved list twice.  We must correct for this
   !
   Edouble_count_old = 0.d0
   CALL LJPSHIFT_INTERACTION_LIST4 (NL_BIN_MOVEONE_OLDCOORDS, V, &
      Edouble_count_old, .false., .false., NATOMS,&
      moved_atoms, nmoved)
   !WRITE(*,*) "EMOVED_OLD", EMOVED_OLD, EDOUBLE_COUNT_OLD
   EMOVED_OLD = EMOVED_OLD - EDOUBLE_COUNT_OLD


   !
   !now update the lists for all the moved particles
   !also update oldcoords
   call nl_bin_moveone_update_lists(COORDS, nl_bin_moveone_nmoved, &
      nl_bin_moveone_moved_atoms)
   !


   !
   !now calculate the new energy contribution of the moved particles
   !
   Emoved_new = 0.d0
   do j1=1,nmoved
      j2 = moved_atoms(j1)
      !CALL LJPSHIFT_INTERACTION_LIST3 (nl_bin_moveone_oldCOORDS, V, &
         !Emoved_new, .false., .false., NATOMS,&
         !NL_BIN_MOVEONE_LISTSAA(:,J2), NLISTSAA(J2), &
         !NL_BIN_MOVEONE_LISTSBB(:,J2), NLISTSBB(J2), &
         !NL_BIN_MOVEONE_LISTSAB(:,J2), NLISTSAB(J2), &
         !0.D0, J2)
      CALL LJPSHIFT_INTERACTION_LIST5 (nl_bin_moveone_oldCOORDS, V, &
         Emoved_new, .false., .false., NATOMS,&
         NL_BIN_MOVEONE_LISTSAA(:,J2), NLISTSAA(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSAA(:,:,j2), &
         NL_BIN_MOVEONE_LISTSBB(:,J2), NLISTSBB(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSBB(:,:,j2), &
         NL_BIN_MOVEONE_LISTSAB(:,J2), NLISTSAB(J2), &
            NL_BIN_MOVEONE_XOFFSET_LISTSAB(:,:,j2), &
         0.D0, J2)
   enddo
   !
   !at this point we've counted the interactions between the atoms
   !in moved list twice.  We must correct for this
   !
   Edouble_count_new = 0.d0
   CALL LJPSHIFT_INTERACTION_LIST4 (NL_BIN_MOVEONE_OLDCOORDS, V, &
      Edouble_count_new, .false., .false., NATOMS,&
      moved_atoms, nmoved)
   !WRITE(*,*) "EMOVED_new", EMOVED_new, EDOUBLE_COUNT_new
   EMOVED_new = EMOVED_new - EDOUBLE_COUNT_new

   !
   !update potential
   !
   potel = old_potel + (emoved_new - emoved_old)
   !write(*,*) "energy one atom", old_potel, potel


   !now update old_potel
   old_potel = potel


   !reset cl_nmoved
   nl_bin_moveone_nmoved = 0

   

END SUBROUTINE LJPSHIFT_ONE_ATOM2
