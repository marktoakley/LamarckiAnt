!   GMIN: A program for finding global minima
!   Copyright (C) 1999- David J. Wales
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
! ds656> Routine for resetting ATOMLISTS so that atoms are
! listed as either all mobile or all frozen, and grouped by 
! type in order of appearance (arbitrary choice).
SUBROUTINE RESET_ATOMLISTS(IGROUP)
  !
  USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, NATOMS, NSPECIES, &
       SPECMASST, SPECMASS, ATMASS, MYUNIT
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: IGROUP
  !
  INTEGER :: I, J, K
  !
  ! Indices in ATOMLISTS:
  ! 1) first index groups atom types / species: 1 -> A, 2 -> B, ...
  ! 2) second index groups mobile (1) and frozen (2) atoms,
  !    determined by IGROUP
  ! 3) third index spans all the atom indices from COORDS that
  !    are in the group, with the zeroth entry giving 
  !    the total size of each list.
  !
  ATOMLISTS(:,:,:) = 0 ! reset to zero
  INVATOMLISTS(:,:) = 0
  !
  ! Separate atoms by type/species and make them all of group IGROUP
  K = 0
  DO I = 1, NSPECIES(0) ! span all species
     ATOMLISTS(I,IGROUP,0) = NSPECIES(I)
     DO J = 1, NSPECIES(I) ! span all atoms of each species
        K = K+1 ! increment atom index
        ATOMLISTS(I,IGROUP,J) = K
        INVATOMLISTS(K,1) = I
        INVATOMLISTS(K,2) = IGROUP
        INVATOMLISTS(K,3) = J
     ENDDO
  ENDDO
  !
  ! Sanity check
  IF(NATOMS /= K) THEN
     WRITE(MYUNIT, '(A)') 'reset_atomlists> Inconsistent atom count!'
     !WRITE(MYUNIT, *) NATOMS, K, NSPECIES(0)
     STOP
  ENDIF
  !
  IF(SPECMASST) THEN
     WRITE(MYUNIT,'(A)') &
          'initialization> (Re)setting masses in accord with SPECMASS:'
     DO I=1,NATOMS
        ATMASS(I) = SPECMASS(INVATOMLISTS(I,1))
        WRITE(MYUNIT,*) ATMASS(I)
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE RESET_ATOMLISTS
!
! ds656> Routine for resetting ATOMLISTS in accord with provided
! list of atomic labels. 
SUBROUTINE SET_ATOMLISTS(LABELS,IGROUP)
  !
  USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, NATOMS, NSPECIES, &
       SPECMASST, SPECMASS, ATMASS, MYUNIT, LFLIPST, QALCSMODE
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: LABELS(NATOMS), IGROUP
  !
  INTEGER :: I,J,TYPECOUNTS(1:NSPECIES(0))
  !
  ATOMLISTS(:,:,:) = 0
  INVATOMLISTS(:,:)=0
  TYPECOUNTS(:) = 0
  !
  !write(MYUNIT,*) 'set_atomlists>', (LABELS(I), I=1,NATOMS)
  DO I=1,NATOMS
     !
     J=LABELS(I)
     TYPECOUNTS(J) = TYPECOUNTS(J) + 1
     !
     ATOMLISTS(J,IGROUP,0) = ATOMLISTS(J,IGROUP,0) + 1
     ATOMLISTS(J,IGROUP,ATOMLISTS(J,IGROUP,0)) = I
     !
     INVATOMLISTS(I,1) = J
     INVATOMLISTS(I,2) = IGROUP
     INVATOMLISTS(I,3) = ATOMLISTS(J,IGROUP,0)
     !
     IF(SPECMASST) ATMASS(I) = SPECMASS(J)
     !
  END DO
  !
  DO I=1,NSPECIES(0)
     IF(LFLIPST.OR.(QALCSMODE.GE.6)) THEN
        NSPECIES(I) = TYPECOUNTS(I)
     ELSE IF(NSPECIES(I) /= TYPECOUNTS(I)) THEN
        WRITE(MYUNIT,'(A, I2, A, I5, I5)') &
             'set_atomlists> Inconsistent counts for atom-type ', &
             I,' ,these numbers should equal:', NSPECIES(I), &
             TYPECOUNTS(I)
        STOP
     ENDIF
  END DO
  !
END SUBROUTINE SET_ATOMLISTS
!
SUBROUTINE CHANGE_ATOM_GROUP(LIST,IGROUP)
  !
  ! Move all atoms in LIST from their current group to IGROUP 
  !
  USE COMMONS, ONLY : ATOMLISTS, INVATOMLISTS, NATOMS, MYUNIT
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: LIST(0:NATOMS),IGROUP
  INTEGER :: I,IT,IG,IP,J,K
  !
  DO K=1,LIST(0) ! span all atoms in list
     !
     I=LIST(K)
     IF(I<1 .OR. I>NATOMS) THEN
        WRITE(MYUNIT,'(A)') 'put_atom_in_group> Failed sanity check!'
        STOP
     ENDIF
     !
     IT=INVATOMLISTS(I,1)
     IG=INVATOMLISTS(I,2)
     IP=INVATOMLISTS(I,3)
     !
     ! Store the last atom index in ATOMLISTS(IT,IG,:) as J
     J=ATOMLISTS(IT,IG,ATOMLISTS(IT,IG,0))
     ! Overwrite the entry of I with J
     ATOMLISTS(IT,IG,IP) = J
     ! Update the inverse list entry of atom J
     INVATOMLISTS(J,3) = IP
     ! Decrement the reduced list length
     ATOMLISTS(IT,IG,0) = ATOMLISTS(IT,IG,0) - 1
     !
     ! Increment the new enlarged list length and store as J
     J = ATOMLISTS(IT,IGROUP,0) + 1
     ! Update the enlarged list length
     ATOMLISTS(IT,IGROUP,0) = J
     ! Update the last entry in the enlarged list
     ATOMLISTS(IT,IGROUP,J) = I
     ! Update inverse list of atom I
     INVATOMLISTS(I,2) = IGROUP
     INVATOMLISTS(I,3) = J
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE CHANGE_ATOM_GROUP
!
SUBROUTINE RESET_STOICHIOMETRY()
  !
  ! ds656> 31/07/2015
  ! Reset stoichiometry to the initial value, with the labels 
  ! distributed randomly.
  !
  USE COMMONS, ONLY : NATOMS, MYUNIT, NSPECIES_INI
  !
  IMPLICIT NONE
  !
  INTEGER :: I, J, K, LABELS(NATOMS)
  !
  LABELS(:) = 0
  K = 0
  DO I=1, NSPECIES_INI(0)
     DO J= 1,NSPECIES_INI(I)
        K=K+1
        IF(K > NATOMS) THEN
           WRITE(MYUNIT,'(A)') 'set_stoichiometry> WTF?!? Terminatig..'
           STOP
        ELSE
           LABELS(K) = I
        ENDIF
     ENDDO
  ENDDO
  !
  CALL KNUTH_SHUFFLE(NATOMS, LABELS)
  CALL SET_ATOMLISTS(LABELS, 1) ! put all in group 1
  !
  RETURN
  !
END SUBROUTINE RESET_STOICHIOMETRY
!
SUBROUTINE KNUTH_SHUFFLE(N,A)
  !
  ! ds656> Randomly permute integer vector A
  !
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(INOUT) :: A(N)
  INTEGER :: I, RANDPOS, TEMP
  DOUBLE PRECISION :: DPRAND
  !
  DO I = N, 2, -1
     RANDPOS = INT(DPRAND()*DBLE(I)) + 1
     TEMP = A(RANDPOS)
     A(RANDPOS) = A(I)
     A(I) = TEMP
  ENDDO
  !
  RETURN
  !
END SUBROUTINE KNUTH_SHUFFLE
!
! ds656> Routine for resetting NBRLISTS so that each
!        atom is counted as its own neighbour, and all 
!        other entries are set to zero.
SUBROUTINE RESET_NBRLISTS()
  !
  USE COMMONS, ONLY : NBRLISTS, NATOMS
  !
  IMPLICIT NONE
  !
  INTEGER :: I
  !
  ! Set all neighbour indices to zero
  NBRLISTS(1,1:NATOMS,2:NATOMS) = 0
  NBRLISTS(2,1:NATOMS,0:NATOMS) = 0
  !
  ! Count each atom as its own (1st) neighbour
  DO I = 1,NATOMS
     NBRLISTS(1,I,0) = 1
     NBRLISTS(1,I,1) = I
  ENDDO
  !
END SUBROUTINE RESET_NBRLISTS
