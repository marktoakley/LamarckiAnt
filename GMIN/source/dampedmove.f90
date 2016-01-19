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

! TO DRIVE DAMPED GROUP MOVES
      SUBROUTINE DAMPEDMOVESTEP(JP)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION :: DPRAND, PI, TWOPI, GROUPROTANGLE, GROUPROTANGLEDEG
      DOUBLE PRECISION :: BVECTOR(3), BVECTOR2(3), BVECTOR3(3), LENGTH, DUMMY
      INTEGER :: I1,JP
! Some helpful parameters
      PI=ATAN(1.0D0)*4
      TWOPI=2.0D0*PI
! For each group....      
      DO I1=1,NDGROUPS
! Rotation along axis of symmetry
         IF(DATOMGROUPPROBA(I1).GE.DPRAND()) THEN
! Group selected to be rotated - calculate rotation angle
            GROUPROTANGLE=(DPRAND()-0.5)*twopi*DATOMGROUPSCALINGA(I1)
!            GROUPROTANGLE= twopi
            GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
! Print some into to GMIN_out for the user
            WRITE(MYUNIT,*) 'DAMPED GROUP MOVE> Rotating along axis group ',TRIM(ADJUSTL(DATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
! Compute the axis of rotation
            BVECTOR(1) = COORDS(3*DATOMGROUPAXIS(I1,1)-2,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP)
            BVECTOR(2) = COORDS(3*DATOMGROUPAXIS(I1,1)-1,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP)
            BVECTOR(3) = COORDS(3*DATOMGROUPAXIS(I1,1)  ,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP)
! Find length   
            LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
            BVECTOR(1)=BVECTOR(1)/LENGTH
            BVECTOR(2)=BVECTOR(2)/LENGTH
            BVECTOR(3)=BVECTOR(3)/LENGTH
! Origin of axis of rotation
            BVECTOR2(1) = COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP)
            BVECTOR2(2) = COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP)
            BVECTOR2(3) = COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP)
! Call the rotation subroutine
            CALL DGROUPROTATION(BVECTOR,BVECTOR2,GROUPROTANGLE,DATOMGROUPS(I1,:),DATOMGROUPATT(I1,:),COORDS(:,JP))
         ENDIF

! Translation along axis of symmetry
         IF(DATOMGROUPPROBC(I1).GE.DPRAND()) THEN
! calculate amount of translation
            GROUPROTANGLE=(DPRAND()-0.5)*DATOMGROUPSCALINGC(I1)
! Print some into to GMIN_out for the user
            WRITE(MYUNIT,*) 'DAMPED GROUP MOVE> Translation along axis group ',TRIM(ADJUSTL(DATOMGROUPNAMES(I1))),' by ',GROUPROTANGLE
! Compute the axis
            BVECTOR(1) = COORDS(3*DATOMGROUPAXIS(I1,1)-2,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP)
            BVECTOR(2) = COORDS(3*DATOMGROUPAXIS(I1,1)-1,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP)
            BVECTOR(3) = COORDS(3*DATOMGROUPAXIS(I1,1)  ,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP)
! Find length   
            LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
            BVECTOR(1)=BVECTOR(1)/LENGTH
            BVECTOR(2)=BVECTOR(2)/LENGTH
            BVECTOR(3)=BVECTOR(3)/LENGTH
! Call the translation subroutine
            CALL DGROUPTRANSLATION(BVECTOR,GROUPROTANGLE,DATOMGROUPS(I1,:),DATOMGROUPATT(I1,:),COORDS(:,JP))
         ENDIF

! Rotation perpendicular to axis of symmetry
         IF(DATOMGROUPPROBB(I1).GE.DPRAND()) THEN
! Group selected to be rotated - calculate rotation angle
            GROUPROTANGLE=(DPRAND()-0.5)*twopi*DATOMGROUPSCALINGB(I1)
            GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
! Print some into to GMIN_out for the user
            WRITE(MYUNIT,*) 'DAMPED GROUP MOVE> Rotating perpendicular to axis group ',TRIM(ADJUSTL(DATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
! Compute the axis of rotation
            BVECTOR(1) = COORDS(3*DATOMGROUPAXIS(I1,1)-2,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP)
            BVECTOR(2) = COORDS(3*DATOMGROUPAXIS(I1,1)-1,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP)
            BVECTOR(3) = COORDS(3*DATOMGROUPAXIS(I1,1)  ,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP)
            LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
            BVECTOR(1)=BVECTOR(1)/LENGTH
            BVECTOR(2)=BVECTOR(2)/LENGTH
            BVECTOR(3)=BVECTOR(3)/LENGTH
            BVECTOR2(1) = DPRAND() - 0.5D0
            BVECTOR2(2) = DPRAND() - 0.5D0
            BVECTOR2(3) = DPRAND() - 0.5D0
            LENGTH=DSQRT(BVECTOR2(1)**2 + BVECTOR2(2)**2 + BVECTOR2(3)**2)
            IF (LENGTH .EQ. 0) THEN
               BVECTOR2(1) = DPRAND() - 0.5D0
               BVECTOR2(2) = DPRAND() - 0.5D0
               BVECTOR2(3) = DPRAND() - 0.5D0
               LENGTH=DSQRT(BVECTOR2(1)**2 + BVECTOR2(2)**2 + BVECTOR2(3)**2)
            ENDIF
            BVECTOR2(1)=BVECTOR2(1)/LENGTH
            BVECTOR2(2)=BVECTOR2(2)/LENGTH
            BVECTOR2(3)=BVECTOR2(3)/LENGTH
            BVECTOR3(1) = BVECTOR(2)*BVECTOR2(3) - BVECTOR2(2)*BVECTOR(3)
            BVECTOR3(2) = BVECTOR(3)*BVECTOR2(1) - BVECTOR2(3)*BVECTOR(1)
            BVECTOR3(3) = BVECTOR(1)*BVECTOR2(2) - BVECTOR2(1)*BVECTOR(2)
            LENGTH=DSQRT(BVECTOR3(1)**2 + BVECTOR3(2)**2 + BVECTOR3(3)**2)
            BVECTOR3(1)=BVECTOR3(1)/LENGTH
            BVECTOR3(2)=BVECTOR3(2)/LENGTH
            BVECTOR3(3)=BVECTOR3(3)/LENGTH
! Origin of axis of rotation
            DUMMY = DPRAND()
            BVECTOR2(1) = (1.0D0 - DUMMY) * COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP) + DUMMY * COORDS(3*DATOMGROUPAXIS(I1,1)-2,JP)
            BVECTOR2(2) = (1.0D0 - DUMMY) * COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP) + DUMMY * COORDS(3*DATOMGROUPAXIS(I1,1)-1,JP)
            BVECTOR2(3) = (1.0D0 - DUMMY) * COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP) + DUMMY * COORDS(3*DATOMGROUPAXIS(I1,1)  ,JP)
! Call the rotation subroutine
            CALL DGROUPROTATION(BVECTOR3,BVECTOR2,GROUPROTANGLE,DATOMGROUPS(I1,:),DATOMGROUPATT(I1,:),COORDS(:,JP))
         ENDIF

! Translation perpendicular to axis of symmetry
         IF(DATOMGROUPPROBD(I1).GE.DPRAND()) THEN
! calculate amount of translation
            GROUPROTANGLE=(DPRAND()-0.5)*DATOMGROUPSCALINGD(I1)
! Print some into to GMIN_out for the user
            WRITE(MYUNIT,*) 'DAMPED GROUP MOVE> Translation perpendicular to axis group ',TRIM(ADJUSTL(DATOMGROUPNAMES(I1))),' by ',GROUPROTANGLE
! Compute the axis
            BVECTOR(1) = COORDS(3*DATOMGROUPAXIS(I1,1)-2,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-2,JP)
            BVECTOR(2) = COORDS(3*DATOMGROUPAXIS(I1,1)-1,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)-1,JP)
            BVECTOR(3) = COORDS(3*DATOMGROUPAXIS(I1,1)  ,JP) - COORDS(3*DATOMGROUPAXIS(I1,2)  ,JP)
            LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
            BVECTOR(1)=BVECTOR(1)/LENGTH
            BVECTOR(2)=BVECTOR(2)/LENGTH
            BVECTOR(3)=BVECTOR(3)/LENGTH
            BVECTOR2(1) = DPRAND() - 0.5D0
            BVECTOR2(2) = DPRAND() - 0.5D0
            BVECTOR2(3) = DPRAND() - 0.5D0
            LENGTH=DSQRT(BVECTOR2(1)**2 + BVECTOR2(2)**2 + BVECTOR2(3)**2)
            IF (LENGTH .EQ. 0) THEN
               BVECTOR2(1) = DPRAND() - 0.5D0
               BVECTOR2(2) = DPRAND() - 0.5D0
               BVECTOR2(3) = DPRAND() - 0.5D0
               LENGTH=DSQRT(BVECTOR2(1)**2 + BVECTOR2(2)**2 + BVECTOR2(3)**2)
            ENDIF
            BVECTOR2(1)=BVECTOR2(1)/LENGTH
            BVECTOR2(2)=BVECTOR2(2)/LENGTH
            BVECTOR2(3)=BVECTOR2(3)/LENGTH
            BVECTOR3(1) = BVECTOR(2)*BVECTOR2(3) - BVECTOR2(2)*BVECTOR(3)
            BVECTOR3(2) = BVECTOR(3)*BVECTOR2(1) - BVECTOR2(3)*BVECTOR(1)
            BVECTOR3(3) = BVECTOR(1)*BVECTOR2(2) - BVECTOR2(1)*BVECTOR(2)
            LENGTH=DSQRT(BVECTOR3(1)**2 + BVECTOR3(2)**2 + BVECTOR3(3)**2)
            BVECTOR3(1)=BVECTOR3(1)/LENGTH
            BVECTOR3(2)=BVECTOR3(2)/LENGTH
            BVECTOR3(3)=BVECTOR3(3)/LENGTH
! Call the translation subroutine
            CALL DGROUPTRANSLATION(BVECTOR3,GROUPROTANGLE,DATOMGROUPS(I1,:),DATOMGROUPATT(I1,:),COORDS(:,JP))
         ENDIF

      ENDDO 
      END SUBROUTINE DAMPEDMOVESTEP

      SUBROUTINE DGROUPROTATION(BVECTOR,BVECTOR2,ANGLE,ATOMINGROUP,ATTENUATION,STEPCOORDS)
      USE commons, ONLY : NATOMS
      IMPLICIT NONE
      INTEGER :: I1
      DOUBLE PRECISION :: BVECTOR(3), BVECTOR2(3), ANGLE, DUMMYMAT(3,3)=0.0D0, ROTMAT(3,3), LOCALATT
      DOUBLE PRECISION :: GROUPATOM(3), GROUPATOMROT(3), STEPCOORDS(3*NATOMS), ATTENUATION(NATOMS)
      LOGICAL :: ATOMINGROUP(NATOMS)

      LOCALATT = 0.0D0
! Scale this vector so its length is the rotation to be done (in radians)
      BVECTOR(1)=BVECTOR(1)*ANGLE
      BVECTOR(2)=BVECTOR(2)*ANGLE
      BVECTOR(3)=BVECTOR(3)*ANGLE

      DO I1=1,NATOMS
         IF (ATOMINGROUP(I1)) THEN
            IF (.NOT.(LOCALATT .EQ. ATTENUATION(I1))) THEN
               LOCALATT = ATTENUATION(I1)
! Get the rotation matrix for this vector axis from RMDRVT
! hk286 - do this only whenever we have different angles, remember that rotation is damped
               CALL RMDRVT(BVECTOR*LOCALATT,ROTMAT,DUMMYMAT,DUMMYMAT,DUMMYMAT,.FALSE.)
            ENDIF
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
            GROUPATOM(1)=STEPCOORDS(3*I1-2)-BVECTOR2(1)
            GROUPATOM(2)=STEPCOORDS(3*I1-1)-BVECTOR2(2)
            GROUPATOM(3)=STEPCOORDS(3*I1  )-BVECTOR2(3)
! Apply the rotation matrix
            GROUPATOMROT=MATMUL(ROTMAT,GROUPATOM)
! Translate back to the origin and copy to COORDS
            STEPCOORDS(3*I1-2)=GROUPATOMROT(1)+BVECTOR2(1)
            STEPCOORDS(3*I1-1)=GROUPATOMROT(2)+BVECTOR2(2)
            STEPCOORDS(3*I1  )=GROUPATOMROT(3)+BVECTOR2(3)
         ENDIF
      ENDDO
      
      END SUBROUTINE DGROUPROTATION

      SUBROUTINE DGROUPTRANSLATION(BVECTOR,STEPSIZE,ATOMINGROUP,ATTENUATION,STEPCOORDS)
      USE commons, ONLY : NATOMS
      IMPLICIT NONE
      INTEGER :: I1
      DOUBLE PRECISION :: BVECTOR(3), STEPSIZE
      DOUBLE PRECISION :: STEPCOORDS(3*NATOMS), ATTENUATION(NATOMS)
      LOGICAL :: ATOMINGROUP(NATOMS)

      DO I1=1,NATOMS
         IF (ATOMINGROUP(I1)) THEN
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
            STEPCOORDS(3*I1-2) = STEPCOORDS(3*I1-2) + BVECTOR(1)*STEPSIZE*ATTENUATION(I1)
            STEPCOORDS(3*I1-1) = STEPCOORDS(3*I1-1) + BVECTOR(2)*STEPSIZE*ATTENUATION(I1)
            STEPCOORDS(3*I1  ) = STEPCOORDS(3*I1  ) + BVECTOR(3)*STEPSIZE*ATTENUATION(I1)
         ENDIF
      ENDDO
      
    END SUBROUTINE DGROUPTRANSLATION
