C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C**************************************************************************
C
C  Subroutine CENTRE moves the centre of mass to the origin.
C
C*********************************************************************
C
      SUBROUTINE CENTRE2(X)
      USE commons
      use genrigid
      IMPLICIT NONE
      DOUBLE PRECISION XMASS, YMASS, ZMASS, DIST, DISTMAX, X(3*NATOMS)
!hk286 - additional variables for generalised rigid body
      DOUBLE PRECISION XRIGIDCOORDS(DEGFREEDOMS)
      INTEGER I, J1
!hk286
         
! hk286 > if generalised rigid body is used, be careful when averaging
! hk286 > no need to shift the rotational degrees of freedom!
      IF ( .NOT. ATOMRIGIDCOORDT ) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NRIGIDBODY
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               XMASS=XMASS+X(6*NRIGIDBODY + 3*J1-2)
               YMASS=YMASS+X(6*NRIGIDBODY + 3*J1-1)
               ZMASS=ZMASS+X(6*NRIGIDBODY + 3*J1  )
            ENDDO
         ENDIF
         XMASS=XMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         YMASS=YMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         ZMASS=ZMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         DO I=1,NRIGIDBODY
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            IF (.NOT.CENTXY) X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               X(6*NRIGIDBODY + 3*J1-2) = X(6*NRIGIDBODY + 3*J1-2) - XMASS
               X(6*NRIGIDBODY + 3*J1-1) = X(6*NRIGIDBODY + 3*J1-1) - YMASS
               IF (.NOT.CENTXY) X(6*NRIGIDBODY+3*J1) = X(6*NRIGIDBODY+3*J1) - ZMASS
            ENDDO
         ENDIF

      ELSE
!hk286 - proceeds as usual if everything is in atom coords

         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         IF (RIGID) THEN
            DO I=1,NATOMS/2
               XMASS=XMASS+X(3*(I-1)+1)
               YMASS=YMASS+X(3*(I-1)+2)
               ZMASS=ZMASS+X(3*(I-1)+3)
            ENDDO
            XMASS=2*XMASS/NATOMS
            YMASS=2*YMASS/NATOMS
            ZMASS=2*ZMASS/NATOMS
C     PRINT*,'initial coordinates in centre:'
C     WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)
            DO I=1,NATOMS/2
               X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
               X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
               IF (.NOT.CENTXY) X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
            ENDDO
         ELSE
            DO I=1,NATOMS
               XMASS=XMASS+X(3*(I-1)+1)
               YMASS=YMASS+X(3*(I-1)+2)
               ZMASS=ZMASS+X(3*(I-1)+3)
            ENDDO
            XMASS=XMASS/NATOMS
            YMASS=YMASS/NATOMS
            ZMASS=ZMASS/NATOMS
            DO I=1,NATOMS
               X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
               X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            ENDDO
            IF (.NOT.CENTXY) THEN
               DO I=1,NATOMS
                  X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
               ENDDO
            ENDIF
         ENDIF
!        IF (DEBUG.AND.(.NOT.CENTXY)) WRITE(MYUNIT,'(A,3G20.10)') 'centre2> centre of mass reset to the origin from ',
!    &        XMASS,YMASS,ZMASS
!        IF (DEBUG.AND.CENTXY) WRITE(MYUNIT,'(A,3G20.10)') 'centre2> centre of mass reset to centre of xy plane from ',
!    &        XMASS,YMASS,ZMASS
C     PRINT*,'final coordinates in centre:'
C     WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)
C
C  Check that all the atoms are in the container. If not then rescale.
C  Must not do this - it could change the system from a minimum!
C
C     IF (RIGID) RETURN
C     DISTMAX=0.0D0
C     DO J1=1,NATOMS
C        DIST=X(3*J1-2)**2+X(3*J1-1)**2+X(3*J1)**2
C        IF (DIST.GT.DISTMAX) DISTMAX=DIST
C     ENDDO
C     IF (DISTMAX.GT.RADIUS) THEN
C        DISTMAX=DSQRT(DISTMAX/RADIUS)*0.99D0
C        DO J1=1,NATOMS
C           X(3*J1-2)=X(3*J1-2)*DISTMAX
C           X(3*J1-1)=X(3*J1-1)*DISTMAX
C           X(3*J1)  =X(3*J1)  *DISTMAX
C        ENDDO
C     ENDIF

      ENDIF
!     hk286 - endif of generalised rigid body if statement

      RETURN
      END

      SUBROUTINE SETCENTRE(X)
      USE commons
      use genrigid

      IMPLICIT NONE
      DOUBLE PRECISION XMASS, YMASS, ZMASS, X(3*NATOMS)
!hk286 - additional variables for generalised rigid body
      DOUBLE PRECISION XRIGIDCOORDS(DEGFREEDOMS)
      INTEGER I, J1

! hk286 > if generalised rigid body is used, be careful when averaging
! hk286 > no need to shift the rotational degrees of freedom!
      IF ( .NOT. ATOMRIGIDCOORDT ) THEN

         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NRIGIDBODY
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               XMASS=XMASS+X(6*NRIGIDBODY + 3*J1-2)
               YMASS=YMASS+X(6*NRIGIDBODY + 3*J1-1)
               ZMASS=ZMASS+X(6*NRIGIDBODY + 3*J1  )
            ENDDO
         ENDIF
         XMASS=XMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         YMASS=YMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         ZMASS=ZMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         DO I=1,NRIGIDBODY
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            IF (.NOT.CENTXY) X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               X(6*NRIGIDBODY + 3*J1-2) = X(6*NRIGIDBODY + 3*J1-2) - XMASS
               X(6*NRIGIDBODY + 3*J1-1) = X(6*NRIGIDBODY + 3*J1-1) - YMASS
               IF (.NOT.CENTXY) X(6*NRIGIDBODY+3*J1) = X(6*NRIGIDBODY+3*J1) - ZMASS
            ENDDO
         ENDIF

      ELSE
!hk286 - proceed as usual if in atom coords

C csw34> XMASS, YMASS and ZMASS are the components of the COM position
C vector
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      IF (RIGID) THEN
C csw34> First need to calculate the position of the COM
         DO I=1,NATOMS/2
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=2*XMASS/NATOMS
         YMASS=2*YMASS/NATOMS
         ZMASS=2*ZMASS/NATOMS
C csw34> Then need to move the COM to the new centre via the origin
         DO I=1,NATOMS/2
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS+CENTX
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS+CENTY
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS+CENTZ
         ENDDO
      ELSE
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
         DO I=1,NATOMS
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS+CENTX
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS+CENTY
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS+CENTZ
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'centre of mass moved to ',CENTX,CENTY,CENTZ

      ENDIF
!     hk286 - endif of generalised rigid body if statement

      RETURN
      END

