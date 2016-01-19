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
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-130 US
! =====================================================================
!
! Routines for performing atom identity swaps in bynary systems.
! ds656> 20/06/2013
!
! =====================================================================
!
SUBROUTINE ABSWAP(NP, DE1, DE2)
  !
  ! ds656> Swap absolutely randomly picked AB pair 
  !
  USE COMMONS,  ONLY : NATOMS, NTYPEA, COORDS, MYUNIT, VT, &
       BGUPTAT, BLJCLUSTER_NOCUT, BASWAPTEST
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  DOUBLE PRECISION, INTENT(OUT) :: DE1, DE2
  DOUBLE PRECISION :: DPRAND, X, VTMOD(NATOMS)
  INTEGER ::  I, I1, I2, J1, J2
  !
  IF ( (NTYPEA == 0) .OR. (NTYPEA == NATOMS) ) RETURN
  !
  J1 = INT(DPRAND()*NTYPEA)
  J2 = INT(DPRAND()*(NATOMS-NTYPEA)) + NTYPEA
  !
  IF(BASWAPTEST) THEN
     !
     IF(BGUPTAT) THEN
        CALL BGUPTA2(COORDS(:,NP),VTMOD)
        DE1 = VTMOD(J1+1) - VT(J1+1)
        DE2 = VTMOD(J2+1) - VT(J2+1)
     ELSEIF(BLJCLUSTER_NOCUT) THEN
        CALL BLJ_CLUST2(COORDS(:,NP),VTMOD)
        DE1 = VTMOD(J1+1) - VT(J1+1)
        DE2 = VTMOD(J2+1) - VT(J2+1)
     ELSE
        WRITE(MYUNIT, '(A)') 'ABswap> DE1 and DE2 set to zero.'
        DE1 = 0.0d0
        DE2 = 0.0d0
     ENDIF
     !
  ENDIF
  !
  I1 = 3*J1
  I2 = 3*J2
  !
  ! NOTE: J1 and J2 are between 0 and NATOMS - 1
  DO I=1,3
     X = COORDS(I1+I, NP)
     COORDS(I1+I, NP) = COORDS(I2+I, NP)
     COORDS(I2+I, NP) = X
  ENDDO
  !
  IF(BASWAPTEST) WRITE(MYUNIT, '(A, I6, A, I6)') 'ABswap> Swapped atoms', J1 +  1, ' and', J2 + 1
  !
  RETURN
  !
END SUBROUTINE ABSWAP
!
! =====================================================================
!
SUBROUTINE ABSWAP2(NP, DE1, DE2)
  !
  ! ds656> Swap AB pair picked using a cheap(ish) energy bias.  
  !
  USE COMMONS, ONLY : NATOMS, VT, NTYPEA, COORDS, MYUNIT, &
       TEMP, BGUPTAT, BLJCLUSTER_NOCUT, BASWAPTEST
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NP
  DOUBLE PRECISION, INTENT(OUT) :: DE1, DE2
  !
  DOUBLE PRECISION VTMOD(NATOMS), VTDIFF(NATOMS), X(3*NATOMS), &
       DE, DUMMY, DPRAND
  INTEGER I, IA, IB, JA, JB
  !
  IF ( (NTYPEA == 0) .OR. (NTYPEA == NATOMS) ) RETURN
  !
  X(:) = COORDS(:,NP)
  !
  IF(BGUPTAT) THEN
     CALL BGUPTA2(X,VTMOD)
  ELSEIF(BLJCLUSTER_NOCUT) THEN
     CALL BLJ_CLUST2(X,VTMOD)
  ELSE
     WRITE(MYUNIT, '(A)') "ABswap2> ERROR: Energy-biased swaps work only for BGUPTA and BLJCLUSTER_NOCUT!"
     STOP
  ENDIF
  !
  ! Accummulte Boltzmann-weighted partial sums for A
  DUMMY=0.0D0
  DO I = 1, NTYPEA
     VTDIFF(I) = DEXP((VT(I) - VTMOD(I))/TEMP(NP))
     DUMMY = DUMMY + VTDIFF(I)
  ENDDO
  DE = DPRAND()*DUMMY
  IA=1
  DUMMY = VTDIFF(IA)
  DO WHILE(DUMMY < DE)
     IA = IA + 1
     DUMMY = DUMMY + VTDIFF(IA)
  END DO
  !
  ! Accummulte Boltzmann-weighted partial sums for B
  DUMMY=0.0D0
  DO I = NTYPEA+1, NATOMS
     VTDIFF(I) = DEXP((VT(I) - VTMOD(I))/TEMP(NP))
     DUMMY = DUMMY + VTDIFF(I)
  ENDDO
  DE = DPRAND()*DUMMY
  IB = NTYPEA+1
  DUMMY = VTDIFF(IB)
  DO WHILE(DUMMY < DE)
     IB = IB + 1
     DUMMY = DUMMY + VTDIFF(IB)
  END DO
  !
  ! NOTE: DPRAND should never yield 1.0d0. 
  ! Now swap coordinates.
  !
  JA = 3*(IA-1)
  JB = 3*(IB-1)
  DO I=1,3
     DUMMY = COORDS(JA+I, NP)
     COORDS(JA+I,NP) = COORDS(JB+I,NP)
     COORDS(JB+I,NP) = DUMMY
  ENDDO
  !
  DE1 = VTMOD(IA) - VT(IA)
  DE2 = VTMOD(IB) - VT(IB)
  !
  IF(BASWAPTEST) WRITE(MYUNIT, '(A,I5,A,I5)') "ABswap2> Swapped atoms: ", IA, " and ", IB
  !
  RETURN
  !
END SUBROUTINE ABSWAP2
