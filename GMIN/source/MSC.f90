!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modIFy
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
!   along with this program; IF not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!******************************************
!   Multicomponent Sutton-Chen potentials
!   ds656> 28/11/2013
!******************************************
!
SUBROUTINE MSC (X,V,PG,GRADT)
  !
  USE COMMONS, ONLY : NATOMS,VT,ATOMLISTS,NSPECIES,&
       CUTT,CUTOFF,PERIODIC,BOX3D
  USE POT_PARAMS, ONLY : MSC_N, MSC_M, MSC_EPS, MSC_A, MSC_C, &
       CUTA_REP, CUTB_REP, CUTC_REP, CUTA_ATT, CUTB_ATT, CUTC_ATT
  !
  IMPLICIT NONE
  !
  ! Parsed variables -----------------------------------------
  LOGICAL, INTENT(IN) :: GRADT
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), PG
  !
  ! Local variables -----------------------------------------
  LOGICAL :: SELF, ADD
  INTEGER :: I,J1,J2,J13,J23,T1,T2,G1,G2,K1,K2,LI1,UI1,LI2,UI2,G2START
  DOUBLE PRECISION :: DIST,DX,RHO(NATOMS),DUMMY,ATT,REP,IBOX3D(3),&
       IDIST,RIJ(NATOMS,NATOMS),VTEMP(NATOMS,3),C1,C2
  !
  IF(PERIODIC) IBOX3D(1:3) = 1.0D0/BOX3D(1:3)
  !
  ! Do some initialisations
  RIJ(1:NATOMS,1:NATOMS) = 0.0D0
  RHO(1:NATOMS) = 0.0D0
  VT(1:NATOMS) = 0.0D0
  !
  ! Double loop over atom types
  LT1: DO T1=1,NSPECIES(0)
     LT2: DO T2=1,T1 
        !
        ! Double loop over mobile / frozen groups.
        LG1: DO G1 = 1,2 ! 1 -> mobile; 2 -> frozen
           !
           IF(T1 == T2) THEN
              G2START = G1
           ELSE
              G2START = 1
           ENDIF
           !
           LG2: DO G2 = G2START,2 ! include frozen-frozen interactions
              !
              LI1=1
              UI1=ATOMLISTS(T1,G1,0)
              !
              IF(T1==T2.AND.G1==G2) THEN
                 SELF = .TRUE.
                 UI1 = UI1 - 1 ! Exclude self-interaction
              ELSE
                 SELF = .FALSE.
              ENDIF
              !
              LA1: DO K1=LI1,UI1
                 !
                 J1 = ATOMLISTS(T1,G1,K1) ! actual atom index
                 J13 = 3*(J1-1) 
                 !
                 ! Atom loop bounds depend on type and group.
                 ! Can this IF block be moved outside of loop LA1?
                 IF(SELF) THEN
                    LI2=K1+1
                    UI2=UI1+1
                 ELSE
                    LI2=1
                    UI2=ATOMLISTS(T2,G2,0)
                 ENDIF
                 !
                 LA2: DO K2=LI2,UI2
                    !
                    J2 = ATOMLISTS(T2,G2,K2) ! actual atom 2 index
                    J23 = 3*(J2-1)
                    !
                    DIST=0.0D0
                    DO I=1,3
                       DX = X(J13+I)-X(J23+I)
                       ! apply PBCs with nearest image convention
                       IF(PERIODIC) THEN 
                          DX = DX - BOX3D(I)*NINT(DX*IBOX3D(I))
                       ENDIF
                       DIST = DIST + DX*DX
                    ENDDO
                    DIST = DSQRT(DIST)
                    !
                    ADD = .TRUE.
                    IF(CUTT) THEN ! check for cutoff
                       IF(DIST.GE.CUTOFF) THEN
                          ADD = .FALSE. ! do not add
                          DIST = -1.0D0 ! flag by non-physical value
                       ENDIF
                    ENDIF
                    !
                    IF(GRADT) THEN
                       RIJ(J2,J1)=DIST ! store for later
                       RIJ(J1,J2)=DIST ! impose symmetry
                    ENDIF
                    !
                    IF(ADD) THEN
                       !
                       ! compute scaled inverse distance
                       DUMMY = MSC_A(T1,T2)/DIST 
                       !
                       ! accumulate the pairwise additive repulsion 
                       ! and contribution to the many-body attraction
                       REP = MSC_EPS(T1,T2)*DUMMY**MSC_N(T1,T2)
                       ATT = DUMMY**MSC_M(T1,T2)
                       !
                       IF(CUTT) THEN ! Add smooth truncation
                          !
                          DUMMY = DIST - CUTOFF 
                          REP = REP + CUTA_REP(T1,T2) + &
                               CUTB_REP(T1,T2)*DUMMY
                          ATT = ATT + CUTA_ATT(T1,T2) + &
                               CUTB_ATT(T1,T2)*DUMMY
                          !
                          DUMMY = 0.5D0*DUMMY*DUMMY
                          REP = REP + CUTC_REP(T1,T2)*DUMMY
                          ATT = ATT + CUTC_ATT(T1,T2)*DUMMY
                          !
                       ENDIF
                       VT(J1) = VT(J1) + REP
                       VT(J2) = VT(J2) + REP
                       RHO(J1) = RHO(J1) + ATT
                       RHO(J2) = RHO(J2) + ATT
                    ENDIF
                    !
                 ENDDO LA2 ! Double loop over atoms
              ENDDO LA1
              !
           ENDDO LG2 ! Double loop over groups (1 and 2 only)
        ENDDO LG1
        !
     ENDDO LT2 ! Double loop over atom types
  ENDDO LT1
  !
  ! Store per-atom energies in array VT for atom groups 1 and 2
  PG=0.0D0
  DO T1=1,NSPECIES(0) ! Span all tom types / species
     C1 = -MSC_C(T1)*MSC_EPS(T1,T1) ! type-dependent coefficient
     C2 = 0.5D0*C1 ! type-dependent coefficient
     DO G1=1,2 ! Span groups 1 and 2
        DO J2=1,ATOMLISTS(T1,G1,0)! loop through atoms
           J1=ATOMLISTS(T1,G1,J2) ! actual atom index
           RHO(J1)=DSQRT(RHO(J1)) ! square root density
           VT(J1) = 0.5D0*VT(J1) + C1*RHO(J1) ! potential per atom
           PG = PG + VT(J1) ! accumulate total potential
           IF(GRADT) THEN ! keep recpirocal of RHO for the gradient
              IF(RHO(J1) > 0.0D0) RHO(J1) = C2 / RHO(J1) 
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  ! Calculate the gradient if required
  !
  IF (GRADT) THEN
     ! 
     ! Do some initialisations
     V(1:3*NATOMS) = 0.0D0
     VTEMP(1:NATOMS,1:3) = 0.0D0
     !
     ! Double loop over atom types
     LT3: DO T1=1,NSPECIES(0)
        LT4: DO T2=1,T1
           !
           ! Precompute some auxiliary type-dependent coefficients
           C1 = -1.0D0/MSC_A(T1,T2)
           C2 = DBLE(MSC_M(T1,T2))*C1 ! Attractive
           C1 = DBLE(MSC_N(T1,T2))*C1*MSC_EPS(T1,T2) ! Repulsive
           !
           ! Double loop over free / frozen groups.
           LG3: DO G1 = 1,2 ! loop over goups 1 and 2
              !
              IF(T1 == T2) THEN
                 G2START = G1
              ELSE
                 G2START = 1
              ENDIF
              !
              LG4: DO G2 = G2START,2 ! loop over groups 1 and 2
                 !
                 LI1=1
                 UI1=ATOMLISTS(T1,G1,0)
                 !
                 IF(T1==T2.AND.G1==G2) THEN
                    SELF = .TRUE.
                    UI1 = UI1 - 1 ! Exclude self-interaction
                 ELSE
                    SELF = .FALSE.
                 ENDIF
                 !
                 LA3: DO K1=LI1,UI1
                    !
                    J1 = ATOMLISTS(T1,G1,K1) ! actual atom 1 index
                    J13 = 3*(J1-1)
                    !
                    ! 2nd atom loop bounds depend on group.
                    ! Can this IF block be moved outside of loop A1?
                    IF(SELF) THEN
                       LI2=K1+1 
                       UI2=UI1+1
                    ELSE
                       LI2=1
                       UI2=ATOMLISTS(T2,G2,0)
                    ENDIF
                    !
                    LA4: DO K2=LI2,UI2
                       !
                       J2 = ATOMLISTS(T2,G2,K2) ! actual atom 2 index
                       J23 = 3*(J2-1)
                       !
                       ! compute DUMMY=dU/dRij within cutoff
                       DIST=RIJ(J1,J2)
                       !
                       IF(DIST.GT.0.0D0) THEN 
                          ! 
                          IDIST = 1.0D0/DIST
                          ! accumulate attractive and repulsive parts
                          DUMMY = MSC_A(T1,T2)*IDIST
                          REP = C1*DUMMY**(MSC_N(T1,T2)+1)
                          ATT = C2*DUMMY**(MSC_M(T1,T2)+1)
                          
                          IF(CUTT) THEN ! apply smooth truncation
                             DUMMY = DIST - CUTOFF
                             REP = REP + CUTB_REP(T1,T2) + &
                                  DUMMY*CUTC_REP(T1,T2)
                             ATT = ATT + CUTB_ATT(T1,T2) + &
                                  DUMMY*CUTC_ATT(T1,T2)
                          ENDIF
                          
                          ! compute (dU/drij)/RIJ. 
                          DUMMY = (REP + ATT*(RHO(J2)+RHO(J1)))*IDIST
                          ! accumulate cartesian components
                          DO I=1,3
                             DX = X(J13+I) - X(J23+I)
                             IF(PERIODIC) THEN 
                                DX = DX-BOX3D(I)*NINT(DX*IBOX3D(I))
                             ENDIF
                             VTEMP(J1,I) = VTEMP(J1,I) + DUMMY*DX
                             VTEMP(J2,I) = VTEMP(J2,I) - DUMMY*DX
                          ENDDO
                          !
                       ENDIF
                       !
                    ENDDO LA4
                 ENDDO LA3
                 !
              ENDDO LG4
           ENDDO LG3
           !
        ENDDO LT4
     ENDDO LT3
     !
     ! Finally, sum the gradient for group 1 atoms only
     DO T1=1,NSPECIES(0) ! Atom types
        DO J2=1,ATOMLISTS(T1,1,0) 
           J1=ATOMLISTS(T1,1,J2)
           J13=3*(J1-1)
           DO I=1,3
              V(J13+I)=V(J13+I)+VTEMP(J1,I)
           ENDDO
        ENDDO
     ENDDO
     !
     !write(*,*) "----------------------------------------------"
     ! TEST PERATOM ENERGIES AND GRADIENTS
     !DO J1 = 1,NATOMS
     !   J2=3*(J1-1)
     !   WRITE(*,'(A,4(1X,F12.6))') "BGUPTA", VT(J1), &
     !        V(J2+1), V(J2+2), V(J2+3)
     !ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MSC
