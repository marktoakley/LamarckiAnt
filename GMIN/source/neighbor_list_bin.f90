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



MODULE BIN_NL_MOD
   !convert a single neighbor list into 3 different neighbor lists depending on
   !interaction type.  AA, AB, BB
   IMPLICIT NONE
   PRIVATE
   INTEGER, ALLOCATABLE :: BIN_NL_AALIST(:,:)
   INTEGER, ALLOCATABLE :: BIN_NL_BBLIST(:,:)
   INTEGER, ALLOCATABLE :: BIN_NL_ABLIST(:,:)
   INTEGER BIN_NL_NAA, BIN_NL_NBB, BIN_NL_NAB
   INTEGER NATOMS, NTYPEA, NTYPEB


   PUBLIC :: BIN_NL_SETUP, BIN_NL_UPDATE, BIN_NL_AALIST, BIN_NL_BBLIST, &
      BIN_NL_ABLIST, BIN_NL_NAA, BIN_NL_NBB, BIN_NL_NAB

   CONTAINS

   SUBROUTINE BIN_NL_SETUP( NATOMS_I, NTYPEA_I)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS_I, NTYPEA_I
      LOGICAL :: FIRST = .TRUE.

      IF (.NOT. FIRST) RETURN
      FIRST = .FALSE.

      NATOMS = NATOMS_I
      NTYPEA = NTYPEA_I
      NTYPEB = NATOMS - NTYPEA
      ALLOCATE( BIN_NL_AALIST(2,NTYPEA * (NTYPEA-1)/2) )
      ALLOCATE( BIN_NL_BBLIST(2,NTYPEB * (NTYPEB-1)/2) )
      ALLOCATE( BIN_NL_ABLIST(2,NTYPEA * NTYPEB) )

   END SUBROUTINE BIN_NL_SETUP

   SUBROUTINE BIN_NL_UPDATE ( LIST_I, NLIST_I )
      !convert the neighbor list to typed AA, BB, and AB neighbor lists
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NLIST_I
      INTEGER, INTENT(IN) :: LIST_I(2,NLIST_I)
      INTEGER J1, J2, J3

      BIN_NL_NAA = 0
      BIN_NL_NBB = 0
      BIN_NL_NAB = 0
      DO J3=1,NLIST_I
         J1=LIST_I(1,J3)
         J2=LIST_I(2,J3)
         IF (J1 .LE. NTYPEA .AND. J2 .LE. NTYPEA) THEN
            BIN_NL_NAA = BIN_NL_NAA + 1
            BIN_NL_AALIST(1,BIN_NL_NAA) = J1
            BIN_NL_AALIST(2,BIN_NL_NAA) = J2
         ELSEIF (J1 .GT. NTYPEA .AND. J2 .GT. NTYPEA) THEN
            BIN_NL_NBB = BIN_NL_NBB + 1
            BIN_NL_BBLIST(1,BIN_NL_NBB) = J1
            BIN_NL_BBLIST(2,BIN_NL_NBB) = J2
         ELSE
            BIN_NL_NAB = BIN_NL_NAB + 1
            BIN_NL_ABLIST(1,BIN_NL_NAB) = J1
            BIN_NL_ABLIST(2,BIN_NL_NAB) = J2
         ENDIF
      ENDDO
      !WRITE(*,*) "bin UPDATE", NLIST_I
      !WRITE(*,*) "          ", BIN_NL_NAA, BIN_NL_NBB, BIN_NL_NAB, BIN_NL_NAA + BIN_NL_NBB + BIN_NL_NAB

   END SUBROUTINE BIN_NL_UPDATE 

END MODULE BIN_NL_MOD
