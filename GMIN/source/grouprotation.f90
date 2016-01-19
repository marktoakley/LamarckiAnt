!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2010 David J. Wales
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

! GROUPROTMOD module - contains subroutines involved in doing GROUPROTATION moves
MODULE GROUPROTMOD

CONTAINS

! ROUTINE TO SCALE GROUPROTATION MOVES
! GROUPROTATION moves can be scaled up and down by changing either the selection
! probability or the maximum rotation amplitude of all/some groups. In this initial
! implementation, a change will be applied to ALL groups. Selective changes will need
! some additional logging routines before they can be implemented - it's on the list!
! Selecting which groups to scale in the future could be done with a mask array.

! TO-DO: Should minimum values be specified in data rather than with an optional argument?
! This would mean we could safely have an optional GROUPMASK argument to control selective
! changes...


! SCALE_PROB and SCALE_ROT: logicals controling what to scale
! FACTOR: the desired scaling factor
! MIN_PROB and MIN_ROT: OPTIONAL arguments specifying minimum values for probability and scaling
! ATOMGROUPSCALING(NGROUPS) contains the rotation amplitude scaling factors for each group
! ATOMGROUPPSELECT(NGROUPS) contains the group selection probabilities for each group
      SUBROUTINE GROUPROTSCALE(SCALE_PROB,SCALE_ROT,FACTOR,MIN_PROB,MIN_ROT)
      USE commons, only: NGROUPS, ATOMGROUPSCALING, ATOMGROUPPSELECT, MYUNIT, DEBUG
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: SCALE_PROB, SCALE_ROT
      DOUBLE PRECISION, INTENT(IN) :: FACTOR
      DOUBLE PRECISION, OPTIONAL :: MIN_PROB, MIN_ROT
      INTEGER :: I1
! Scale selection probability
      IF(SCALE_PROB) THEN
         ATOMGROUPPSELECT(:)=ATOMGROUPPSELECT(:)*FACTOR
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> all group selection probabilities scaled by ',FACTOR
      ENDIF
! Scale rotation amplitude
      IF(SCALE_ROT) THEN 
         ATOMGROUPSCALING(:)=ATOMGROUPSCALING(:)*FACTOR
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> all group rotation amplitudes scaled by ',FACTOR
      ENDIF
! Sanity checks
      DO I1=1,NGROUPS
! 1. Rotation amplitude and rotation scaling should not exceed 1.0
         ATOMGROUPPSELECT(I1)=MIN(ATOMGROUPPSELECT(I1),1.0D0)
         ATOMGROUPSCALING(I1)=MIN(ATOMGROUPSCALING(I1),1.0D0)
! 2. If provided, check that the minimum values (MINPROB, MINROT) are respected
         IF(PRESENT(MIN_PROB)) ATOMGROUPPSELECT(I1)=MAX(ATOMGROUPPSELECT(I1),MIN_PROB)
         IF(PRESENT(MIN_ROT)) ATOMGROUPSCALING(I1)=MAX(ATOMGROUPSCALING(I1),MIN_ROT) 
      ENDDO
! 3. Warn user if selection probabilities all set to 1 if only they are being scaled
      IF((SUM(ATOMGROUPPSELECT)/NGROUPS.EQ.1.0D0).AND.(SCALE_PROB.AND.(.NOT.SCALE_ROT))) THEN
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> WARNING: selection probability set to 1.0 for all groups'
      ENDIF
! DEBUG PRINTING
      IF(DEBUG) THEN
         PRINT *,"GROUPROTSCALE> SCALE_PROB,SCALE_ROT,FACTOR=",SCALE_PROB,SCALE_ROT,FACTOR
         PRINT *,"GROUPROTSCALE> ATOMGROUPSCALING(:)=",ATOMGROUPSCALING(:)
         PRINT *,"GROUPROTSCALE> ATOMGROUPPSELECT(:)=",ATOMGROUPPSELECT(:)
      ENDIF
 
      END SUBROUTINE GROUPROTSCALE

! ROUTINE TO DRIVE GROUPROTATION MOVES
! JP is the parallel run ID so that only the appropriate coordinates are
! altered during parallel runs.
      SUBROUTINE GROUPROTSTEP(JP)
      USE commons
      USE modcharmm, only: CHNEIGHBOURT
      IMPLICIT NONE
      DOUBLE PRECISION :: DPRAND, PI, TWOPI, GROUPROTANGLE, GROUPROTANGLEDEG, MYRANDOM
      INTEGER :: I1,JP,RESC,RES1,RES2
! Some helpful parameters
      PI=ATAN(1.0D0)*4
      TWOPI=2.0D0*PI
! For each group....      
      IF ( USERPOTT .OR. AMBERT .OR. AMBER12T .OR. (CHRMMT.AND.(.NOT.CHNEIGHBOURT)) ) THEN
         DO I1=1,NGROUPS
            IF (ATOMGROUPPSELECT(I1).GE.DPRAND()) THEN
! Group selected to be rotated - calculate rotation angle
               GROUPROTANGLE=(DPRAND()-0.5)*twopi*ATOMGROUPSCALING(I1)
               GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
! Print some into to GMIN_out for the user
               IF (.NOT. GROUPROT_SUPPRESS) THEN
                  WRITE(MYUNIT,*) 'GROUPROTATION> Rotating group ',TRIM(ADJUSTL(ATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
               END IF
! Call the rotation subroutine
               CALL GROUPROTATION(ATOMGROUPAXIS(I1,1),ATOMGROUPAXIS(I1,2),GROUPROTANGLE,ATOMGROUPS(I1,:),COORDS(:,JP))
            ENDIF
         ENDDO 
      ENDIF
      IF (CHRMMT.AND.CHNEIGHBOURT) THEN
! pick group and number of groups to be rotated (between 2 and 4)
          MYRANDOM=DPRAND()
          RESC=1 + NINT(MYRANDOM*(NGROUPS-1))
          MYRANDOM=DPRAND()
          RES1=RESC - NINT(2*MYRANDOM)
          IF (RES1.LT.1) RES1=1
          MYRANDOM=DPRAND()
          RES2=RESC + NINT(2*MYRANDOM)
          IF (RES2.GT.NGROUPS) RES2=NGROUPS
          DO I1=RES1,RES2
             GROUPROTANGLE=(DPRAND()-0.5)*twopi*ATOMGROUPSCALING(I1)
             GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
! Print some into to GMIN_out for the user
             IF (.NOT. GROUPROT_SUPPRESS) THEN
                WRITE(MYUNIT,*) 'GROUPROTATION> Rotating group ',TRIM(ADJUSTL(ATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
             END IF
! Call the rotation subroutine
             CALL GROUPROTATION(ATOMGROUPAXIS(I1,1),ATOMGROUPAXIS(I1,2),GROUPROTANGLE,ATOMGROUPS(I1,:),COORDS(:,JP))
          ENDDO
      ENDIF

      END SUBROUTINE GROUPROTSTEP

! The GROUPROTATION subroutine allows for almost any rotation of a defined set of atoms.
! The rotation axis is defined by two atoms (BATOMS1 and BATOM2), the group of atoms to 
! rotate is defined by the logical array ATOMINGROUP (if element is .TRUE., that atom is
! in the group), ANGLE is the rotation angle in radians and STEPCOORDS contains the current
! atomic coordinates.
!
      SUBROUTINE GROUPROTATION(BATOM1,BATOM2,ANGLE,ATOMINGROUP,STEPCOORDS)
      USE commons
      USE MOVES
      IMPLICIT NONE
      INTEGER :: BATOM1, BATOM2, I1
      DOUBLE PRECISION :: BVECTOR(3), LENGTH, ANGLE, DUMMYMAT(3,3)=0.0D0, ROTMAT(3,3)
      DOUBLE PRECISION :: GROUPATOM(3), GROUPATOMROT(3), STEPCOORDS(3*NATOMS)
      LOGICAL :: ATOMINGROUP(NATOMS)
! ===============================TESTING ROTATE_ABOUT_AXIS============================
!      DOUBLE PRECISION, DIMENSION(3*NATOMS)  :: COORDS_COPY
!      INTEGER, DIMENSION(:), ALLOCATABLE     :: ATOM_LIST
!      INTEGER                                :: NUM_ROTATING_ATOMS
!      INTEGER                                :: I2
!      DOUBLE PRECISION, DIMENSION(3)         :: AXIS_START, AXIS_END
!      DOUBLE PRECISION                       :: PI, ANGLE_DEGREES
!
!! Take a copy of the coordinates
!      COORDS_COPY(:) = STEPCOORDS(:)
!
!! Work out how many atoms are rotating to allocate ATOM_LIST
!      NUM_ROTATING_ATOMS = 0
!      DO I1 = 1, NATOMS
!         IF (ATOMINGROUP(I1) .EQV. .TRUE.) THEN
!            NUM_ROTATING_ATOMS = NUM_ROTATING_ATOMS + 1
!         END IF
!      END DO
!      ALLOCATE(ATOM_LIST(NUM_ROTATING_ATOMS))
!
!! Assign ATOM_LIST
!      I2 = 1
!      DO I1 = 1, NATOMS
!         IF (ATOMINGROUP(I1) .EQV. .TRUE.) THEN
!            ATOM_LIST(I2) = I1
!            I2 = I2 + 1
!         END IF
!      END DO
!
!! Assign the start and ends of the axis
!      AXIS_START(1) = COORDS_COPY(3 * BATOM2 - 2)
!      AXIS_START(2) = COORDS_COPY(3 * BATOM2 - 1)
!      AXIS_START(3) = COORDS_COPY(3 * BATOM2    )
!      AXIS_END(1)   = COORDS_COPY(3 * BATOM1 - 2)
!      AXIS_END(2)   = COORDS_COPY(3 * BATOM1 - 1)
!      AXIS_END(3)   = COORDS_COPY(3 * BATOM1    )
!
!! Convert the angle from radians to degrees
!      PI = 4.0D0 * ATAN(1.0D0)
!      ANGLE_DEGREES = 180.0D0 * ANGLE / PI
!
!! Call the subroutine
!      CALL ROTATION_ABOUT_AXIS(COORDS_COPY, AXIS_START, AXIS_END, &
!                               ANGLE_DEGREES, ATOM_LIST)
! ===============================TESTING ROTATE_ABOUT_AXIS============================

! STEP 1
! Produce notmalised bond vector corresponding to the rotation axis
! BATOM1 and BATOM2 are the atoms defining this vector
      BVECTOR(1)=STEPCOORDS(3*BATOM1-2)-STEPCOORDS(3*BATOM2-2)
      BVECTOR(2)=STEPCOORDS(3*BATOM1-1)-STEPCOORDS(3*BATOM2-1)
      BVECTOR(3)=STEPCOORDS(3*BATOM1  )-STEPCOORDS(3*BATOM2  )
! Find length   
      LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
      BVECTOR(1)=BVECTOR(1)/LENGTH
      BVECTOR(2)=BVECTOR(2)/LENGTH
      BVECTOR(3)=BVECTOR(3)/LENGTH
! STEP 2
! Scale this vector so its length is the rotation to be done (in radians)
      BVECTOR(1)=BVECTOR(1)*ANGLE
      BVECTOR(2)=BVECTOR(2)*ANGLE
      BVECTOR(3)=BVECTOR(3)*ANGLE
! STEP 3
! Get the rotation matrix for this vector axis from RMDRVT
! Interface:
! SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)
! P is an un-normalised vector you wish to rotate around. Its length equals the desired rotation in radians
! RM will return the 3x3 rotation matrix
! DRM1-3 are derivative matricies, not needed here
! GTEST is also not needed so set to .FALSE.
      CALL RMDRVT(BVECTOR,ROTMAT,DUMMYMAT,DUMMYMAT,DUMMYMAT,.FALSE.)
! STEP 4
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
      DO I1=1,NATOMS
         IF (ATOMINGROUP(I1)) THEN
            GROUPATOM(1)=STEPCOORDS(3*I1-2)-STEPCOORDS(3*BATOM2-2)
            GROUPATOM(2)=STEPCOORDS(3*I1-1)-STEPCOORDS(3*BATOM2-1)
            GROUPATOM(3)=STEPCOORDS(3*I1  )-STEPCOORDS(3*BATOM2  )
! Apply the rotation matrix
            GROUPATOMROT=MATMUL(ROTMAT,GROUPATOM)
! Translate back to the origin and copy to COORDS
            STEPCOORDS(3*I1-2)=GROUPATOMROT(1)+STEPCOORDS(3*BATOM2-2)
            STEPCOORDS(3*I1-1)=GROUPATOMROT(2)+STEPCOORDS(3*BATOM2-1)
            STEPCOORDS(3*I1  )=GROUPATOMROT(3)+STEPCOORDS(3*BATOM2  )
         ENDIF
      ENDDO

! ===============================TESTING ROTATE_ABOUT_AXIS============================
!      OPEN(UNIT=8293, FILE='testing_rotation', POSITION='APPEND')
!      WRITE(8293, '(10I5)') ATOM_LIST
!      DO I1 = 1, NATOMS
!         WRITE(8293, '(I5)') I1
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1-2), STEPCOORDS(3*I1-2), COORDS_COPY(3*I1-2)-STEPCOORDS(3*I1-2)
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1-1), STEPCOORDS(3*I1-1), COORDS_COPY(3*I1-1)-STEPCOORDS(3*I1-1)
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1  ), STEPCOORDS(3*I1  ), COORDS_COPY(3*I1  )-STEPCOORDS(3*I1  )
!      END DO
!      WRITE(8293, '(A40)') '=================================================='
!      CLOSE(8293) 
! ===============================TESTING ROTATE_ABOUT_AXIS============================
      END SUBROUTINE GROUPROTATION

END MODULE GROUPROTMOD
