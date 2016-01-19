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
!=================================================================
! This module contains subroutines that parse arguments for 
! keywords specifying a potential with many parameters. The
! idea is to remove all the potential-specific clutter from 
! keywords.F and contain it here.
! -ds656 (Dec 2014)
!=================================================================
!
MODULE PARSE_POT_PARAMS
  !
  USE COMMONS, ONLY : NATOMS, NSPECIES, MYUNIT
  !
  IMPLICIT NONE
  !
  !-------------------------------------------------------------
  ! To use NITEMS from the common block BUFINF, we must declare 
  ! lots of other crap that is included in BUFINF.
  LOGICAL :: SKIPBL, CLEAR, ECHO, CAT
  INTEGER :: ITEM, NITEMS, LOC, LINE, NCR, NERROR, LAST
  COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, &
       NCR, NERROR, ECHO, LAST, CAT
  !-------------------------------------------------------------
  !
  LOGICAL :: END 
  CHARACTER(LEN=16) :: WORD  
  !
CONTAINS
  !
  SUBROUTINE PARSE_MLJ_PARAMS(DATA_UNIT)
    !
    ! Parse parameters for Multicomponent Lennard-Jones system.
    !
    USE POT_PARAMS, ONLY : MLJ_EPS, MLJ_SIG
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: DATA_UNIT
    !
    INTEGER :: N0,N1,N2,N3,NONTYPEA,NTYPEA
    !
    ! First check the number of species from the first argument
    CALL READI(N0) 
    IF (NSPECIES(0) /= N0) THEN
       WRITE(MYUNIT,'(A)') &
            'parse_pot_params> Inconsistent species count for MLJ!'
       STOP
    ENDIF
    !
    IF(ALLOCATED(MLJ_EPS)) DEALLOCATE(MLJ_EPS)
    IF(ALLOCATED(MLJ_SIG)) DEALLOCATE(MLJ_SIG)
    !
    ALLOCATE(MLJ_EPS(N0,N0),MLJ_SIG(N0,N0))
    !
    CALL READF(MLJ_EPS(1,1))
    CALL READF(MLJ_SIG(1,1))
    !
    NONTYPEA=0
    DO N1 = 1,NSPECIES(0) 
       !
       IF (N1==1) THEN ! skip first line
          N0=2
       ELSE
          N0=N1
       ENDIF
       !
       DO N2 = N0,NSPECIES(0)
          !
          CALL INPUT(END, DATA_UNIT)
          IF (.NOT. END) THEN
             CALL READU(WORD)
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MLJ_params> Bad line count in ''data''!'
             STOP
          ENDIF
          !
          IF(WORD .EQ. 'MLJ') THEN
             !
             IF(N1==N2) THEN
                N3 = 4
                CALL READI(NSPECIES(N1))
                NONTYPEA = NONTYPEA + NSPECIES(N1)
             ELSE
                N3 = 3
             ENDIF
             !
             IF(NITEMS < N3) THEN
                WRITE(MYUNIT,'(A)') &
                   'parse_MLJ_param> Insufficient param count!'
                STOP
             ENDIF
             !
             CALL READF(MLJ_EPS(N1,N2))
             MLJ_EPS(N2,N1) = MLJ_EPS(N1,N2) ! Impose symmetry
             CALL READF(MLJ_SIG(N1,N2))
             MLJ_SIG(N2,N1) = MLJ_SIG(N1,N2)
             !
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MLJ_params> Missing ''MLJ'' header!'
             STOP
          ENDIF
          !
       ENDDO
       !
    ENDDO
    !
    ! Atom number for species 1 has not yet been specified.        
    NTYPEA = NATOMS - NONTYPEA
    NSPECIES(1) = NTYPEA
    !                                                              
    WRITE(MYUNIT,'(A,I4,A)') &
         'parse_MLJ_params> Lennard-Jones system with', &
         NSPECIES(0),' species.' 
    WRITE(MYUNIT, '(A)') &
         'parse_MLJ_params> Atom count for each species:'
    DO N1=1,NSPECIES(0)
       WRITE(MYUNIT,'(I6)',ADVANCE='NO') NSPECIES(N1)
    ENDDO
    WRITE(MYUNIT,'(A)') ' '
    WRITE(MYUNIT,'(A)') &
         'parse_MLJ_params> SPEC_I, SPEC_J, EPS, SIG:'
    DO N1 = 1,NSPECIES(0)
       DO N2 = N1,NSPECIES(0)
          !
          WRITE(MYUNIT,'(2(1X,I3),2(1X,1pE12.6E1))') &
               N1,N2,MLJ_EPS(N1,N2),MLJ_SIG(N1,N2)
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE PARSE_MLJ_PARAMS
  !
  SUBROUTINE PARSE_MGUPTA_PARAMS(DATA_UNIT)
    !
    ! Parse parameters for Multicomponent Gupta system.
    !
    USE POT_PARAMS, ONLY : MGUPTA_A, MGUPTA_XI, MGUPTA_P, &
         MGUPTA_Q, MGUPTA_R0, MGUPTA_M2Q, MGUPTA_XI_SQ, &
         MGUPTA_MP_DIVBY_R0, MGUPTA_M2Q_DIVBY_R0
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: DATA_UNIT
    !
    INTEGER :: N0,N1,N2,N3,NONTYPEA,NTYPEA
    !
    ! First check the number of species from the first argument
    CALL READI(N0)
    IF (NSPECIES(0) /= N0) THEN
       WRITE(MYUNIT,'(A)') &
            'parse_pot_params> Inconsistent species count for MGUPTA!'
       STOP
    ENDIF
    !
    IF(ALLOCATED(MGUPTA_A)) DEALLOCATE(MGUPTA_A)
    IF(ALLOCATED(MGUPTA_P)) DEALLOCATE(MGUPTA_P)
    IF(ALLOCATED(MGUPTA_Q)) DEALLOCATE(MGUPTA_Q)
    IF(ALLOCATED(MGUPTA_XI)) DEALLOCATE(MGUPTA_XI)
    IF(ALLOCATED(MGUPTA_R0)) DEALLOCATE(MGUPTA_R0)
    !
    IF(ALLOCATED(MGUPTA_M2Q)) DEALLOCATE(MGUPTA_M2Q)
    IF(ALLOCATED(MGUPTA_XI_SQ)) DEALLOCATE(MGUPTA_XI_SQ)
    IF(ALLOCATED(MGUPTA_MP_DIVBY_R0)) DEALLOCATE(MGUPTA_MP_DIVBY_R0)
    IF(ALLOCATED(MGUPTA_M2Q_DIVBY_R0)) &
         DEALLOCATE(MGUPTA_M2Q_DIVBY_R0)
    !
    ALLOCATE(MGUPTA_A(N0,N0),MGUPTA_XI(N0,N0))
    ALLOCATE(MGUPTA_P(N0,N0),MGUPTA_Q(N0,N0),MGUPTA_R0(N0,N0))
    ALLOCATE(MGUPTA_M2Q(N0,N0),MGUPTA_XI_SQ(N0,N0))
    ALLOCATE(MGUPTA_MP_DIVBY_R0(N0,N0))
    ALLOCATE(MGUPTA_M2Q_DIVBY_R0(N0,N0))
    !
    CALL READF(MGUPTA_A(1,1))
    CALL READF(MGUPTA_P(1,1))
    CALL READF(MGUPTA_Q(1,1))
    CALL READF(MGUPTA_XI(1,1))
    CALL READF(MGUPTA_R0(1,1))
    !
    NONTYPEA=0
    DO N1 = 1,NSPECIES(0) 
       !
       IF (N1==1) THEN ! skip first line
          N0=2
       ELSE
          N0=N1
       ENDIF
       !
       DO N2 = N0,NSPECIES(0)
          !
          CALL INPUT(END, DATA_UNIT)
          IF (.NOT. END) THEN
             CALL READU(WORD)
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MGupta_params> Bad line count in ''data''!'
             STOP
          ENDIF
          !
          IF(WORD .EQ. 'MGUPTA') THEN
             !
             IF(N1==N2) THEN
                N3 = 7
                CALL READI(NSPECIES(N1))
                NONTYPEA = NONTYPEA + NSPECIES(N1)
             ELSE
                N3 = 6
             ENDIF
             !
             IF(NITEMS < N3) THEN
                WRITE(MYUNIT,'(A)') &
                   'parse_MGupta_param> Insufficient param count!'
                STOP
             ENDIF
             !
             CALL READF(MGUPTA_A(N1,N2))
             MGUPTA_A(N2,N1) = MGUPTA_A(N1,N2) ! Impose symmetry
             CALL READF(MGUPTA_P(N1,N2))
             MGUPTA_P(N2,N1) = MGUPTA_P(N1,N2)
             CALL READF(MGUPTA_Q(N1,N2))
             MGUPTA_Q(N2,N1) = MGUPTA_Q(N1,N2)
             CALL READF(MGUPTA_XI(N1,N2))
             MGUPTA_XI(N2,N1) = MGUPTA_XI(N1,N2)
             CALL READF(MGUPTA_R0(N1,N2))
             MGUPTA_R0(N2,N1) = MGUPTA_R0(N1,N2)
             !
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MGupta_params> Missing ''MGUPTA'' header!'
             STOP
          ENDIF
          !
       ENDDO
       !
    ENDDO
    !
    ! Also pre-compute some auxilliary constants for speed
    DO N1 = 1,NSPECIES(0) 
       DO N2 = N1,NSPECIES(0)
          !
          MGUPTA_XI_SQ(N1,N2) = MGUPTA_XI(N1,N2)**2
          MGUPTA_XI_SQ(N2,N1) = MGUPTA_XI_SQ(N1,N2)
          !
          MGUPTA_MP_DIVBY_R0(N1,N2) = &
               -MGUPTA_P(N1,N2)/MGUPTA_R0(N1,N2)
          MGUPTA_MP_DIVBY_R0(N2,N1) = MGUPTA_MP_DIVBY_R0(N1,N2)
          !
          MGUPTA_M2Q(N1,N2) = -2.0D0*MGUPTA_Q(N1,N2)
          MGUPTA_M2Q(N2,N1) = MGUPTA_M2Q(N1,N2)
          !
          MGUPTA_M2Q_DIVBY_R0(N1,N2) = &
               MGUPTA_M2Q(N1,N2)/MGUPTA_R0(N1,N2)
          MGUPTA_M2Q_DIVBY_R0(N2,N1) = MGUPTA_M2Q_DIVBY_R0(N1,N2)
          !
       ENDDO
    ENDDO
    !
    ! Atom number for species 1 has not yet been specified.        
    NTYPEA = NATOMS - NONTYPEA
    NSPECIES(1) = NTYPEA
    !                                                              
    WRITE(MYUNIT,'(A,I4,A)') &
         'parse_MGupta_params> Gupta system with', &
         NSPECIES(0),' species.' 
    WRITE(MYUNIT, '(A)') &
         'parse_MGupta_params> Atom count for each species:'
    DO N1=1,NSPECIES(0)
       WRITE(MYUNIT,'(I6)',ADVANCE='NO') NSPECIES(N1)
    ENDDO
    WRITE(MYUNIT,'(A)') ' '
    WRITE(MYUNIT,'(A)') &
         'parse_MGupta_params> SPEC_I, SPEC_J, A, P, Q, XI, R0:'
    DO N1 = 1,NSPECIES(0)
       DO N2 = N1,NSPECIES(0)
          !
          WRITE(MYUNIT,'(2(1X,I3),5(1X,1pE12.6E1))') &
               N1,N2,MGUPTA_A(N1,N2),MGUPTA_P(N1,N2), &
               MGUPTA_Q(N1,N2),MGUPTA_XI(N1,N2),MGUPTA_R0(N1,N2)
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE PARSE_MGUPTA_PARAMS
  !
  SUBROUTINE PARSE_MSC_PARAMS(DATA_UNIT)
    !
    USE POT_PARAMS, ONLY : MSC_N, MSC_M, MSC_EPS, MSC_A, MSC_C

    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: DATA_UNIT
    INTEGER :: J1,J2,J3,J4,NONTYPEA,NTYPEA
    !
    CALL READI(J1) ! First argument gives the number of species
    IF (NSPECIES(0) /= J1) THEN
       WRITE(MYUNIT,'(A)') &
            'parse_pot_params> Inconsistent species count for MSC!'
       STOP
    ENDIF
    !
    IF(ALLOCATED(MSC_N)) DEALLOCATE(MSC_N)
    IF(ALLOCATED(MSC_M)) DEALLOCATE(MSC_M)
    IF(ALLOCATED(MSC_EPS)) DEALLOCATE(MSC_EPS)
    IF(ALLOCATED(MSC_C)) DEALLOCATE(MSC_C)
    IF(ALLOCATED(MSC_A)) DEALLOCATE(MSC_A)
    ALLOCATE(MSC_N(J1,J1),MSC_M(J1,J1))
    ALLOCATE(MSC_EPS(J1,J1),MSC_A(J1,J1),MSC_C(J1))
    !
    CALL READI(MSC_N(1,1))
    CALL READI(MSC_M(1,1))
    CALL READF(MSC_A(1,1))
    CALL READF(MSC_EPS(1,1))
    CALL READF(MSC_C(1))
    !
    NONTYPEA=0
    DO J1 = 1,NSPECIES(0) ! read all the parameters
       IF(J1==1) THEN
          J3 = 2
       ELSE
          J3 = J1
       ENDIF
       DO J2 = J3,NSPECIES(0)
          CALL INPUT(END, DATA_UNIT)
          IF (.NOT. END) THEN
             CALL READU(WORD)
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MSC_params> Bad line count in ''data''!'
             STOP
          ENDIF
          IF(WORD .EQ. 'MSC') THEN
             IF(J2==J1) THEN
                J4 = 7
                CALL READI(NSPECIES(J1))
                NONTYPEA = NONTYPEA + NSPECIES(J1)
             ELSE
                J4 = 5
             ENDIF
             IF(NITEMS < J4) THEN
                WRITE(MYUNIT,'(A)') &
                     'parse_MSC_params> Insufficient param count!'
                STOP
             ENDIF
             CALL READI(MSC_N(J1,J2))
             MSC_N(J2,J1) = MSC_N(J1,J2) ! Impose symmetry
             CALL READI(MSC_M(J1,J2))
             MSC_M(J2,J1) = MSC_M(J1,J2)
             CALL READF(MSC_A(J1,J2))
             MSC_A(J2,J1) = MSC_A(J1,J2)
             CALL READF(MSC_EPS(J1,J2))
             MSC_EPS(J2,J1) = MSC_EPS(J1,J2)
             IF(J1.EQ.J2) CALL READF(MSC_C(J1))
          ELSE
             WRITE(MYUNIT,'(A)') &
                  'parse_MSC_params> Missing ''MSC'' header!'
             STOP
          ENDIF
       ENDDO
    ENDDO
    !
    ! Atom number for species 1 has not yet been specified.
    NTYPEA = NATOMS - NONTYPEA
    NSPECIES(1) = NTYPEA
    !
    WRITE(MYUNIT,'(A,I4,A)') &
         'parse_MSC_params> Sutton-Chen system with ',NSPECIES(0),&
         ' species. Atom count for each species:'
    DO J1=1,NSPECIES(0)
       WRITE(MYUNIT,'(I6)',ADVANCE='NO') NSPECIES(J1)
    ENDDO
    WRITE(MYUNIT,'(A)') ' '
    WRITE(MYUNIT,'(A)') &
         'parse_MSC_params> SPEC_I, SPEC_J, N, M, A, EPS, (C):'
    DO J1 = 1,NSPECIES(0)
       DO J2 = J1,NSPECIES(0)
          IF(J1.EQ.J2) THEN
             WRITE(MYUNIT,'(4(1X,I3),3(1X,1pE12.6E1))') &
                  J1,J2,MSC_N(J1,J2),MSC_M(J1,J2),MSC_A(J1,J2),&
                  MSC_EPS(J1,J2),MSC_C(J1)
          ELSE
             WRITE(MYUNIT,'(4(1X,I3),2(1X,1pE12.6E1))') &
                  J1,J2,MSC_N(J1,J2),MSC_M(J1,J2),MSC_A(J1,J2), &
                  MSC_EPS(J1,J2)
          ENDIF
       ENDDO
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE PARSE_MSC_PARAMS
  !
END MODULE PARSE_POT_PARAMS
