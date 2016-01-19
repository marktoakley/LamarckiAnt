!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!   Loop structure recoded by J.A. Elliott 2009
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
!*************************************************************************
!
!   BINARY GUPTA POTENTIAL
!   CJH 12/2011
!
!   Loop structure changes 
!   ds656 07/2013
!
!*************************************************************************
!
MODULE BGUPMOD
  !
  DOUBLE PRECISION :: AAA,AAB,ABB,PAA,PAB,PBB,QAA,QAB,QBB,ZAA,ZAB,ZBB,&
       R0AA,R0AB,R0BB
  DOUBLE PRECISION :: AARRAY(2,2),AARRAY2(2,2),ZARRAY(2,2),&
       PARRAY(2,2), QARRAY(2,2), R0ARRAY(2,2)
  ! Generate parameter arrays with atom types reversed
  DOUBLE PRECISION :: AARRAY_MOD(2,2), ZARRAY_MOD(2,2), &
       PARRAY_MOD(2,2), QARRAY_MOD(2,2), R0ARRAY_MOD(2,2)  
  CHARACTER(LEN=2) :: BGUPTANAME1,BGUPTANAME2
  !
END MODULE BGUPMOD
!
SUBROUTINE BGUPTA (X,V,PG,GRADT)
  !
  USE COMMONS, ONLY : NATOMS,VT,ATOMLISTS,NSPECIES
  USE BGUPMOD
  !
  IMPLICIT NONE
  !
  ! Parsed variables -----------------------------------------
  LOGICAL, INTENT(IN) :: GRADT
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), PG
  !
  ! Local variables -----------------------------------------
  LOGICAL :: SELF
  INTEGER :: I,J1,J2,J13,J23,ATOMTYPE,MBI,T1,T2,G1,G2,K1,K2,&
       LI1,UI1,LI2,UI2,G2START
  DOUBLE PRECISION :: DIST, DX, DY, DZ, &
       GRHO(NATOMS), VTEMP, DUMMY, DISTANCE_MATRIX(NATOMS,NATOMS), &
       VTEMP1(NATOMS), VTEMP2(NATOMS), VTEMP3(NATOMS), PWR1, PWR2
  !
  ! Do some initialisations
  DISTANCE_MATRIX(1:NATOMS,1:NATOMS)=0.0D0
  GRHO(1:NATOMS)=0.0D0
  VT(1:NATOMS) = 0.0D0
  !
  MBI = 2 ! Many-Body Index (2 for Gupta; 1 for Morse-like pair pot.)
  PWR1 = 1.0D0/DBLE(MBI)
  PWR2 = PWR1-1.0D0
  !
  ! Double loop over atom types
  LT1: DO T1=1,NSPECIES(0)
     LT2: DO T2=1,T1 ! T1,NSPECIES(0)
        !
        ! Double loop over mobile and frozen atom groups
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
                 J1 = ATOMLISTS(T1,G1,K1) ! actual atom 1 index
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
                    !TESTCOUNT=TESTCOUNT+1
                    !WRITE(*,*) "BGupta> TEST:", J1,J2,T1,T2
                    !
                    DIST=0.0D0
                    DO I=1,3
                       DIST = DIST + (X(J23+I)-X(J13+I))**2
                    ENDDO
                    DIST = DSQRT(DIST)
                    !
                    DISTANCE_MATRIX(J2,J1)=DIST ! store distance
                    DISTANCE_MATRIX(J1,J2)=DIST ! impose symmetry
                    !
                    DIST=DIST/R0ARRAY(T1,T2)
                    !
                    DUMMY = AARRAY(T1,T2)*DEXP(PARRAY(T1,T2)*(1-DIST))
                    VT(J1) = VT(J1) + DUMMY
                    VT(J2) = VT(J2) + DUMMY
                    !
                    ! calculate many-body potential term:
                    DUMMY = DEXP(2.0D0*QARRAY(T1,T2)*(1-DIST)) 
                    DUMMY = DUMMY * ZARRAY(T1,T2)**MBI
                    GRHO(J1) = GRHO(J1) + DUMMY
                    GRHO(J2) = GRHO(J2) + DUMMY
                    !
                 ENDDO LA2 ! Double loop over atoms
              ENDDO LA1
              !
           ENDDO LG2 ! Span dynamic (1) and frozen (2) groups only
        ENDDO LG1
        !
     ENDDO LT2 ! Double loop over atom types
  ENDDO LT1
  !
  ! Now store per-atom energies in array VT and sum over all atoms
  !
  PG=0.0D0
  DO G1=1,2 ! Groups 1 and 2
     DO T1=1,NSPECIES(0) ! Atom types
        DO J2=1,ATOMLISTS(T1,G1,0) 
           J1=ATOMLISTS(T1,G1,J2)
           !GRHO(J1)=DSQRT(GRHO(J1)) ! square root density
           VT(J1) = VT(J1) - GRHO(J1)**PWR1 ! accummulate per-atom energy
           PG = PG + VT(J1)       ! accumulate potential energy
        ENDDO
     ENDDO
  ENDDO
  !
  ! Calculate gradient terms, if required
  !
  IF (GRADT) THEN
     ! 
     V(:) = 0.0d0
     VTEMP1(:)=0.0D0
     VTEMP2(:)=0.0D0
     VTEMP3(:)=0.0D0
     !
     ! Double loop over atom types
     LT3: DO T1=1,NSPECIES(0)
        LT4: DO T2=1,T1 
           !
           ! Double loop over free / frozen groups
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
                    ! store reciprocal of RHO for atom J1
                    DUMMY=GRHO(J1)**PWR2 
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
                       DIST=DISTANCE_MATRIX(J1,J2) ! recall distance 
                       DIST=DIST/R0ARRAY(T1,T2)
                       !
                       VTEMP = 2.0D0*( PWR1*QARRAY(T1,T2)* &
                            (ZARRAY(T1,T2)**MBI)* &
                            (DUMMY+GRHO(J2)**PWR2)* &
                            DEXP(2.0D0*QARRAY(T1,T2)* &
                            (1-DIST)) - AARRAY(T1,T2)*PARRAY(T1,T2)* &
                            DEXP(PARRAY(T1,T2)*(1-DIST))) &
                            /(R0ARRAY(T1,T2)**2*DIST) ! ARRAY2 vs ARRAY!!
                       ! was ARRAY2, but with MBI need ARRAY1 !
                       !
                       ! Cartesian components of distance
                       DX=(X(J13+1)-X(J23+1)) 
                       DY=(X(J13+2)-X(J23+2))
                       DZ=(X(J13+3)-X(J23+3))
                       !
                       VTEMP1(J1)=VTEMP1(J1)+VTEMP*DX 
                       VTEMP2(J1)=VTEMP2(J1)+VTEMP*DY
                       VTEMP3(J1)=VTEMP3(J1)+VTEMP*DZ
                       !
                       ! impose symmetry on gradient
                       VTEMP1(J2)=VTEMP1(J2)-VTEMP*DX 
                       VTEMP2(J2)=VTEMP2(J2)-VTEMP*DY
                       VTEMP3(J2)=VTEMP3(J2)-VTEMP*DZ
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
           V(J13+1)=V(J13+1)+VTEMP1(J1)
           V(J13+2)=V(J13+2)+VTEMP2(J1)
           V(J13+3)=V(J13+3)+VTEMP3(J1)
        ENDDO
     ENDDO
     !
     ! TEST PERATOM ENERGIES AND GRADIENTS
     ! DO T1=1,NSPECIES(0)
     !    DO G1=1,2
     !       WRITE(*,'(A,3(1X,I5))') "BGUPTA> T,G,NAT(T,G)=", &
     !            T1, G1, ATOMLISTS(T1,G1,0)
     !       DO J2=1,ATOMLISTS(T1,G1,0)
     !          J1=ATOMLISTS(T1,G1,J2)
     !          J13=3*(J1-1)
     !          WRITE(*,'(A,1X, I4, 4(1X,F12.6))') "BGUPTA>", &
     !               J1, VT(J1), V(J13+1), V(J13+2), V(J13+3)
     !       ENDDO
     !    ENDDO
     ! ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE BGUPTA
!
SUBROUTINE BGUPTA2 (X,VTMOD)
  !
  !     Calculate what the energy of each atom would be if that particular
  !     atom was of different type. The geometry is absolutely fixed.
  !
  USE COMMONS, ONLY : NATOMS, NTYPEA, NSPECIES
  USE BGUPMOD
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) ! coordinates
  DOUBLE PRECISION, INTENT(OUT) :: VTMOD(NATOMS) ! potential
  INTEGER :: I,J1,J2,J13,J23,MBI,LI1,UI1,LI2,UI2,T1,T2,TINDEX(2,2)
  DOUBLE PRECISION :: DIST,DIST1,DIST2,GRHO(NATOMS),DUMMY1,DUMMY2
  DOUBLE PRECISION :: PWR1
  !
  ! Do some initialisations
  GRHO(:)=0.0D0
  VTMOD(1:NATOMS) = 0.0D0
  MBI = 2 ! (2 for Gupta; 1 gives Morse-like pair potential)
  PWR1 = 1.0D0/DBLE(MBI)
  !
  TINDEX(1,1:2) = (/1,NTYPEA/)
  TINDEX(2,1:2) = (/NTYPEA+1,NATOMS/)
  !
  LT1: DO T1=1,NSPECIES(0)
     LT2: DO T2=T1,NSPECIES(0)
        !
        LI1=TINDEX(T1,1)
        UI1=TINDEX(T1,2)
        IF(T1==T2) UI1 = UI1 - 1
        !
        LA1: DO J1=LI1,UI1
           !
           ! Atom loop bounds depend on group.
           ! Can this IF block be moved outise of loop A1?
           IF(T1==T2) THEN
              LI2=J1+1 
              UI2=UI1+1
           ELSE
              LI2=TINDEX(T2,1)
              UI2=TINDEX(T2,2)
           ENDIF
           !
           LA2: DO J2=LI2,UI2
              !
              J13 = 3*(J1-1)
              J23 = 3*(J2-1)
              !
              ! Calc and store distance between atoms J1,J2
              DIST=0.0D0
              DO I=1,3
                 DIST = DIST + (X(J23+I)-X(J13+I))**2
              ENDDO
              DIST = DSQRT(DIST)
              !
              DIST1=DIST/R0ARRAY_MOD(T1,T2)
              DIST2=DIST/R0ARRAY_MOD(T2,T1)
              !
              DUMMY1 = AARRAY_MOD(T1,T2)*DEXP(PARRAY_MOD(T1,T2)*(1-DIST1))
              DUMMY2 = AARRAY_MOD(T2,T1)*DEXP(PARRAY_MOD(T2,T1)*(1-DIST2))
              !VTMOD(J1) = VTMOD(J1) + 2.0D0*DUMMY1 
              !VTMOD(J2) = VTMOD(J2) + 2.0D0*DUMMY2
              VTMOD(J1) = VTMOD(J1) + DUMMY1 
              VTMOD(J2) = VTMOD(J2) + DUMMY2
              !
              ! calculate many-body potential term
              !
              DUMMY1 = DEXP(2.0D0*QARRAY_MOD(T1,T2)*(1-DIST1))*&
                   ZARRAY_MOD(T1,T2)**MBI
              DUMMY2 = DEXP(2.0D0*QARRAY_MOD(T2,T1)*(1-DIST2))*&
                   ZARRAY_MOD(T2,T1)**MBI
              !GRHO(J1) = GRHO(J1) + 2.0D0*DUMMY1
              !GRHO(J2) = GRHO(J2) + 2.0D0*DUMMY2
              GRHO(J1) = GRHO(J1) + DUMMY1
              GRHO(J2) = GRHO(J2) + DUMMY2
              !
           ENDDO LA2
        ENDDO LA1
        !
     ENDDO LT2
  ENDDO LT1
  !
  !     Now store the potential energy of each atom in array VTMOD
  !
  DO J1=1,NATOMS
     !GRHO(J1)=DSQRT(GRHO(J1)) ! square root density
     VTMOD(J1) = VTMOD(J1) - GRHO(J1)**PWR1 ! accummulate per-atom energy
  ENDDO
  !
  RETURN
  !
END SUBROUTINE BGUPTA2

