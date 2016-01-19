C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C  Finds the minimum distance between two geometries using LBFGS
C  Geometry R0 is reset to R1 after each optimisation step. Geometry
C  in R should not change, so neither should RA. RB is returned as the
C  closest geometry to RA.
C
C jmc As long as zsym isn't 'W' (in which case mind does something special) mind
C doesn't care what atomic symbol we give it
C
      SUBROUTINE MINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
      USE COMMONS,ONLY : FREEZE, PULLT, EFIELDT ! hk286
      IMPLICIT NONE
      INTEGER J1, IG, I, ITER, J2, NCOUNT
      DOUBLE PRECISION DPRAND
      DOUBLE PRECISION P(3), DIST, DIST0, RG, RG0, DISTFUNC, 
     1                 MYROTMAT(3,3), 
     2                 OVEC(3), H1VEC(3), H2VEC(3), DSAVE, OMEGATOT(3,3), RA(*), RB(*)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R, R0, R1, R1SAVE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RBSAVE
      INTEGER NSIZE, NATOMS
      LOGICAL BULKT, MFLAG, TWOD, AGAIN, PRESERVET
      COMMON /GEOM/ NSIZE
      COMMON /MINDOM/ MYROTMAT, OMEGATOT
      CHARACTER(LEN=5) ZUSE

! hk286
!     IF (GTHOMSONT) THEN
!        CALL GTHOMSONMINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
!        RETURN
!     ENDIF

      ALLOCATE(RBSAVE(3*NATOMS))
      RBSAVE(1:3*NATOMS)=RB(1:3*NATOMS)
C
C  Initialise accumulated rotation matrix to the identity.
C
      DO J1=1,3
         DO J2=1,3
            OMEGATOT(J2,J1)=0.0D0
         ENDDO
         OMEGATOT(J1,J1)=1.0D0
      ENDDO

      AGAIN=.FALSE.
      DSAVE=0.0D0
      IF (ZUSE(1:1).EQ.'W') THEN
         ALLOCATE(R(3,3*(NATOMS/2)),R0(3,3*(NATOMS/2)),R1(3,3*(NATOMS/2)),R1SAVE(3,3*(NATOMS/2)))
C        NSIZE=3*NATOMS
         NSIZE=3*(NATOMS/2)
C        DO J1=1,NATOMS
         DO J1=1,NATOMS/2
            CALL CONVERT(RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
     1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
            R(1,(J1-1)*3+1)=OVEC(1)
            R(2,(J1-1)*3+1)=OVEC(2)
            R(3,(J1-1)*3+1)=OVEC(3)
            R(1,(J1-1)*3+2)=H1VEC(1)
            R(2,(J1-1)*3+2)=H1VEC(2)
            R(3,(J1-1)*3+2)=H1VEC(3)
            R(1,(J1-1)*3+3)=H2VEC(1)
            R(2,(J1-1)*3+3)=H2VEC(2)
            R(3,(J1-1)*3+3)=H2VEC(3)
            CALL CONVERT(RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
C    1            RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
     1            RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
            R0(1,(J1-1)*3+1)=OVEC(1)
            R0(2,(J1-1)*3+1)=OVEC(2)
            R0(3,(J1-1)*3+1)=OVEC(3)
            R0(1,(J1-1)*3+2)=H1VEC(1)
            R0(2,(J1-1)*3+2)=H1VEC(2)
            R0(3,(J1-1)*3+2)=H1VEC(3)
            R0(1,(J1-1)*3+3)=H2VEC(1)
            R0(2,(J1-1)*3+3)=H2VEC(2)
            R0(3,(J1-1)*3+3)=H2VEC(3)
         ENDDO
      ELSE
         ALLOCATE(R(3,NATOMS),R0(3,NATOMS),R1(3,NATOMS),R1SAVE(3,NATOMS))
         NSIZE=NATOMS
         DO J1=1,NSIZE
            R(1,J1)=RA(3*(J1-1)+1)
            R(2,J1)=RA(3*(J1-1)+2)
            R(3,J1)=RA(3*(J1-1)+3)
            R0(1,J1)=RB(3*(J1-1)+1)
            R0(2,J1)=RB(3*(J1-1)+2)
            R0(3,J1)=RB(3*(J1-1)+3)
         ENDDO
      ENDIF
C
C     put c.o.m. to origin
C
      IF (.NOT.FREEZE) THEN
         do ig=1,3
            rg=0.d0
            rg0=0.d0
            do i=1,nsize
               rg=rg+R(ig,i)
               rg0=rg0+R0(ig,i)
            enddo
            rg=rg/nsize
            rg0=rg0/nsize
            do i=1,nsize
               R(ig,i)=R(ig,i)-rg
               R0(ig,i)=R0(ig,i)-rg0
            enddo
         enddo
      ENDIF
C
C     initial angles
C
      P(1)=0.0D0
      P(2)=0.0D0
      P(3)=0.0D0
C     IF (TWOD) P(1)=0.0D0
C     IF (TWOD) P(2)=0.0D0
C
C     calculate initial distance
C
      NCOUNT=0
10    DIST0=DISTFUNC(P,R,R0,R1)
      DIST0=SQRT(DIST0)
      IF (BULKT) THEN
         DIST=DIST0
         IF (PRESERVET) RB(1:3*NATOMS)=RBSAVE(1:3*NATOMS)
         DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)
         RETURN
      ENDIF

      CALL MMYLBFGS(P,1.0D-7,MFLAG,DIST,500,ITER,R,R0,R1)
      DIST=DISTFUNC(P,R,R0,R1)
      DIST=SQRT(DIST)
      IF (MFLAG) THEN
C        WRITE(*,'(A,2F15.5,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER
C        PRINT*
      ELSE
         NCOUNT=NCOUNT+1
         IF (NCOUNT.GT.0) THEN 
            PRINT*,'convergence failure in mind'
c           STOP
         ELSE
C           WRITE(*,'(A,2F15.5,A,I6,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER,' NCOUNT=',NCOUNT
            DO J1=1,NSIZE
               R0(1,J1)=R1(1,J1)
               R0(2,J1)=R1(2,J1)
               R0(3,J1)=R1(3,J1)
            ENDDO
            IF (.NOT.(TWOD.OR.PULLT.OR.EFIELDT)) P(1)=2*(DPRAND()-0.5D0)/10.0D0
            IF (.NOT.(TWOD.OR.PULLT.OR.EFIELDT)) P(2)=2*(DPRAND()-0.5D0)/10.0D0
            P(3)=2*(DPRAND()-0.5D0)/10.0D0
            GOTO 10
         ENDIF
      ENDIF
C
C  This block allows the second geometry to rotate out of the xy plane; not allowed!
C
C     IF (TWOD) THEN
C        IF (AGAIN) THEN
C           WRITE(*,'(A,2F20.10)') 'DIST,DSAVE=',DIST,DSAVE
C           IF (DIST.GT.DSAVE) THEN
C              DIST=DSAVE
C              DO J1=1,NSIZE
C                 R1(1,J1)=R1SAVE(1,J1)
C                 R1(2,J1)=R1SAVE(2,J1)
C                 R1(3,J1)=R1SAVE(3,J1)
C                 R0(1,J1)=-R0(1,J1)
C              ENDDO
C              DO J1=1,3
C                 DO J2=1,3
C                    OSAVE(J2,J1)=OSAVESAVE(J2,J1)
C                 ENDDO
C              ENDDO
C           ENDIF
C        ELSE
C           AGAIN=.TRUE.
C           DSAVE=DIST
C           DO J1=1,NSIZE
C              R0(1,J1)=-R0(1,J1)
C              R1SAVE(1,J1)=R1(1,J1)
C              R1SAVE(2,J1)=R1(2,J1)
C              R1SAVE(3,J1)=R1(3,J1)
C           ENDDO
C           DO J1=1,3
C              DO J2=1,3
C                 OSAVE(J2,J1)=0.0D0
C                 OSAVESAVE(J2,J1)=OSAVE(J2,J1)
C              ENDDO
C              OSAVE(J1,J1)=1.0D0
C           ENDDO
C           GOTO 10
C        ENDIF
C     ENDIF

      IF (ZUSE(1:1).EQ.'W') THEN
C        DO J1=1,NATOMS  !  WCOMMENT
         DO J1=1,NATOMS/2
            OVEC(1)=R(1,(J1-1)*3+1)
            OVEC(2)=R(2,(J1-1)*3+1)
            OVEC(3)=R(3,(J1-1)*3+1)
            H1VEC(1)=R(1,(J1-1)*3+2)
            H1VEC(2)=R(2,(J1-1)*3+2)
            H1VEC(3)=R(3,(J1-1)*3+2)
            H2VEC(1)=R(1,(J1-1)*3+3)
            H2VEC(2)=R(2,(J1-1)*3+3)
            H2VEC(3)=R(3,(J1-1)*3+3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,OVEC =',J1,OVEC(1),OVEC(2),OVEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H1VEC=',J1,H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H2VEC=',J1,H2VEC(1),H2VEC(2),H2VEC(3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
     1                 RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3))
C    1                 RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3)) ! WCOMMENT
C           WRITE(*,'(A,I6,6F15.5)') 'J1,RA=',J1,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3)
            OVEC(1)=R1(1,(J1-1)*3+1)
            OVEC(2)=R1(2,(J1-1)*3+1)
            OVEC(3)=R1(3,(J1-1)*3+1)
            H1VEC(1)=R1(1,(J1-1)*3+2)
            H1VEC(2)=R1(2,(J1-1)*3+2)
            H1VEC(3)=R1(3,(J1-1)*3+2)
            H2VEC(1)=R1(1,(J1-1)*3+3)
            H2VEC(2)=R1(2,(J1-1)*3+3)
            H2VEC(3)=R1(3,(J1-1)*3+3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,OVEC =',J1,OVEC(1),OVEC(2),OVEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H1VEC=',J1,H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H2VEC=',J1,H2VEC(1),H2VEC(2),H2VEC(3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
C    1                 RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3))
     1                 RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3))
C           WRITE(*,'(A,I6,6F15.5)') 'J1,RA=',J1,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3)
         ENDDO
      ELSE
         DO J1=1,NSIZE
            RB(3*(J1-1)+1)=R1(1,J1)
            RB(3*(J1-1)+2)=R1(2,J1)
            RB(3*(J1-1)+3)=R1(3,J1)
            RA(3*(J1-1)+1)=R(1,J1)
            RA(3*(J1-1)+2)=R(2,J1)
            RA(3*(J1-1)+3)=R(3,J1)
         ENDDO
      ENDIF

      DO J1=1,3
         DO J2=1,3
            MYROTMAT(J2,J1)=OMEGATOT(J2,J1)
         ENDDO
      ENDDO

      IF (PRESERVET) RB(1:3*NATOMS)=RBSAVE(1:3*NATOMS)

      DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)

      RETURN
      END

c__________________________________________________________________________

        subroutine myinverse(A,B)
        IMPLICIT NONE
        DOUBLE PRECISION A(3,3),B(3,3),C(2,2),D(3,3),DET
        INTEGER I, J, JJ, II

        do i=1,3
           do j=1,3
              do jj=1,j-1
                 do ii=1,i-1
                    C(ii,jj)=A(ii,jj)
                 enddo
                 do ii=i+1,3
                    C(ii-1,jj)=A(ii,jj)
                 enddo
              enddo
              
              do jj=j+1,3
                 do ii=1,i-1
                    C(ii,jj-1)=A(ii,jj)
                 enddo
                 do ii=i+1,3
                    C(ii-1,jj-1)=A(ii,jj)
                    
                 enddo
              enddo
              
              B(i,j)=C(1,1)*C(2,2)-C(1,2)*C(2,1)
              D(i,j)=B(i,j)*(-1)**(i+j)
           enddo
        enddo
        
        det=A(1,1)*D(1,1)+A(1,2)*D(1,2)+A(1,3)*D(1,3)
        IF (DET.EQ.0.0D0) WRITE(*,'(A,G20.10)') 'ERROR, determinant in mind routine inverse=',DET
        
        do i=1,3
           do j=1,3
              B(i,j)=D(j,i)/det
           enddo
        enddo
        
        return
        end
c_____________________________________________________________________

      SUBROUTINE PROD(A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I, J

      do i=1,3
         do j=1,3
            C(i,j)=A(i,1)*B(1,j)+A(i,2)*B(2,j)+A(i,3)*B(3,j)
         enddo
      enddo

      return
      end
c_____________________________________________________________________

      function distfunc(p,r,r0,r1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),A(3,3),R1(3,*),DIST,DISTFUNC
      DOUBLE PRECISION MYROTMAT(3,3),ROTX(3,3),ROTY(3,3),ROTZ(3,3),OMEGATOT(3,3)
      DOUBLE PRECISION P(3)
      integer nsize,I,J
      COMMON /MINDOM/ MYROTMAT, OMEGATOT

      common/geom/nsize

      rotx(1,1)=1.0D0
      rotx(1,2)=0.0D0
      rotx(2,1)=0.0D0
      rotx(1,3)=0.0D0
      rotx(3,1)=0.0D0
      rotx(2,2)=dcos(p(1))
      rotx(3,3)=rotx(2,2)
      rotx(2,3)=dsin(p(1))
      rotx(3,2)=-rotx(2,3)

      roty(1,2)=0.0D0
      roty(2,1)=0.0D0
      roty(2,3)=0.0D0
      roty(3,2)=0.0D0
      roty(2,2)=1.0D0
      roty(1,1)=dcos(p(2))
      roty(3,3)=roty(1,1)
      roty(1,3)=-dsin(p(2))
      roty(3,1)=-roty(1,3)

      rotz(1,3)=0.0D0
      rotz(3,1)=0.0D0
      rotz(2,3)=0.0D0
      rotz(3,2)=0.0D0
      rotz(3,3)=1.0D0
      rotz(1,1)=dcos(p(3))
      rotz(2,2)=rotz(1,1)
      rotz(1,2)=dsin(p(3))
      rotz(2,1)=-rotz(1,2)

      A(1,1)=ROTY(1,1)*ROTZ(1,1)+ROTY(1,2)*ROTZ(2,1)+ROTY(1,3)*ROTZ(3,1)
      A(1,2)=ROTY(1,1)*ROTZ(1,2)+ROTY(1,2)*ROTZ(2,2)+ROTY(1,3)*ROTZ(3,2)
      A(1,3)=ROTY(1,1)*ROTZ(1,3)+ROTY(1,2)*ROTZ(2,3)+ROTY(1,3)*ROTZ(3,3)
      A(2,1)=ROTY(2,1)*ROTZ(1,1)+ROTY(2,2)*ROTZ(2,1)+ROTY(2,3)*ROTZ(3,1)
      A(2,2)=ROTY(2,1)*ROTZ(1,2)+ROTY(2,2)*ROTZ(2,2)+ROTY(2,3)*ROTZ(3,2)
      A(2,3)=ROTY(2,1)*ROTZ(1,3)+ROTY(2,2)*ROTZ(2,3)+ROTY(2,3)*ROTZ(3,3)
      A(3,1)=ROTY(3,1)*ROTZ(1,1)+ROTY(3,2)*ROTZ(2,1)+ROTY(3,3)*ROTZ(3,1)
      A(3,2)=ROTY(3,1)*ROTZ(1,2)+ROTY(3,2)*ROTZ(2,2)+ROTY(3,3)*ROTZ(3,2)
      A(3,3)=ROTY(3,1)*ROTZ(1,3)+ROTY(3,2)*ROTZ(2,3)+ROTY(3,3)*ROTZ(3,3)

      MYROTMAT(1,1)=ROTX(1,1)*A(1,1)+ROTX(1,2)*A(2,1)+ROTX(1,3)*A(3,1)
      MYROTMAT(1,2)=ROTX(1,1)*A(1,2)+ROTX(1,2)*A(2,2)+ROTX(1,3)*A(3,2)
      MYROTMAT(1,3)=ROTX(1,1)*A(1,3)+ROTX(1,2)*A(2,3)+ROTX(1,3)*A(3,3)
      MYROTMAT(2,1)=ROTX(2,1)*A(1,1)+ROTX(2,2)*A(2,1)+ROTX(2,3)*A(3,1)
      MYROTMAT(2,2)=ROTX(2,1)*A(1,2)+ROTX(2,2)*A(2,2)+ROTX(2,3)*A(3,2)
      MYROTMAT(2,3)=ROTX(2,1)*A(1,3)+ROTX(2,2)*A(2,3)+ROTX(2,3)*A(3,3)
      MYROTMAT(3,1)=ROTX(3,1)*A(1,1)+ROTX(3,2)*A(2,1)+ROTX(3,3)*A(3,1)
      MYROTMAT(3,2)=ROTX(3,1)*A(1,2)+ROTX(3,2)*A(2,2)+ROTX(3,3)*A(3,2)
      MYROTMAT(3,3)=ROTX(3,1)*A(1,3)+ROTX(3,2)*A(2,3)+ROTX(3,3)*A(3,3)

C     call prod(roty,rotz,A)
C     call prod(rotx,A,MYROTMAT)

      do i=1,nsize
         do j=1,3
            R1(j,i)=MYROTMAT(j,1)*R0(1,i)+MYROTMAT(j,2)*R0(2,i)+MYROTMAT(j,3)*R0(3,i)
         enddo
      enddo

      dist=0.d0
      do i=1,nsize
         dist=dist+(r(1,i)-r1(1,i))**2+(r(2,i)-r1(2,i))**2+(r(3,i)-r1(3,i))**2
C        WRITE(*,'(A,I5,2G20.10)') 'I,DIST,TOTAL=',I,(r(1,i)-r1(1,i))**2+(r(2,i)-r1(2,i))**2+(r(3,i)-r1(3,i))**2,DIST
      enddo

      DISTFUNC=dist
      
      return
      end
c_____________________________________________________________________

      SUBROUTINE DDISTFUNC(P,F,R,R0,R1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),R1(3,*)
      DOUBLE PRECISION P(3),F(3)
      integer nsize, I

      common/geom/nsize

      F(1)=0.d0
      F(2)=0.d0
      F(3)=0.d0
      do i=1,nsize
         F(1)=F(1)+R(2,i)*R1(3,i)-R(3,i)*R1(2,i)
         F(2)=F(2)+R(3,i)*R1(1,i)-R(1,i)*R1(3,i)
         F(3)=F(3)+R(1,i)*R1(2,i)-R(2,i)*R1(1,i)
      enddo
      F(1)=-2.d0*F(1)
      F(2)=-2.d0*F(2)
      F(3)=-2.d0*F(3)

      return

      end
c_____________________________________________________________________

      subroutine DDDISTFUNC(P,H,R,R0,R1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),R1(3,*)
      DOUBLE PRECISION P(3),H(3,3)
      integer nsize, I, J

      common/geom/nsize

      do i=1,3
         do j=1,3
            H(i,j)=0.d0
            H(i,j)=0.d0
            H(i,j)=0.d0
         enddo
      enddo
      do i=1,nsize
         H(1,1)=H(1,1)+R(2,i)*R1(2,i)+R(3,i)*R1(3,i)
         H(2,2)=H(2,2)+R(1,i)*R1(1,i)+R(3,i)*R1(3,i)
         H(3,3)=H(3,3)+R(1,i)*R1(1,i)+R(2,i)*R1(2,i)
         H(2,1)=H(2,1)-R(2,i)*R1(1,i)
         H(3,1)=H(3,1)-R(3,i)*R1(1,i)
         H(3,2)=H(3,2)-R(3,i)*R1(2,i)
      enddo
      H(1,2)=H(2,1)
      H(1,3)=H(3,1)
      H(2,3)=H(3,2)

      do i=1,3
         do j=1,3
            H(i,j)=2.d0*H(i,j)
         enddo
      enddo

      return

      end

      SUBROUTINE ROTGEOM(NATOMS,COORDS)
      IMPLICIT NONE
      INTEGER I, J, K, NATOMS
      DOUBLE PRECISION COORDS(*), R1, R0(3), MYROTMAT(3,3),OMEGATOT(3,3)
      COMMON /MINDOM/ MYROTMAT, OMEGATOT

      DO I=1,NATOMS
         R0(1)=COORDS(3*(I-1)+1)
         R0(2)=COORDS(3*(I-1)+2)
         R0(3)=COORDS(3*(I-1)+3)
         DO J=1,3
            R1=0.0D0
            DO K=1,3
               R1=R1+MYROTMAT(J,K)*R0(K)
            ENDDO
            COORDS(3*(I-1)+J)=R1
         ENDDO
      ENDDO

      RETURN
      END

C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C
      SUBROUTINE MMYLBFGS(X,EPS,MFLAG,ENERGY,ITMAX,ITDONE,R,R0,R1)
      USE COMMONS
      use porfuncs
      IMPLICIT NONE
      INTEGER N,M,J1,J2,ITMAX,ITDONE,NFAIL,J
      PARAMETER (N=3,M=5)  !  MMUPDATE is actually ignored
      DOUBLE PRECISION X(N),G(3),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,OVERLAP,DISTFUNC
      DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,ALPHA,GSAVE(3)
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,MYROTMAT(3,3),OMEGATOT(3,3),OTEMP(3,3)
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
      LOGICAL MFLAG, PTEST
      INTEGER NSIZE, ISTAT
      DOUBLE PRECISION R(3,NSIZE), R0(3,NSIZE), R1(3,NSIZE), DOT1, DOT2, MAXMBFGS
      COMMON /GEOM/ NSIZE
      COMMON /MINDOM/ MYROTMAT, OMEGATOT
C
C  SGI appears to need this SAVE statement!
C
      SAVE

      PTEST=.TRUE.
      PTEST=.FALSE.
      ALPHA=1.0D0
      NFAIL=0
      ITER=0
      ITDONE=0
      ENERGY=DISTFUNC(X,R,R0,R1)
      CALL DDISTFUNC(X,GSAVE,R,R0,R1)
! hk286 
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) GSAVE(1)=0.0D0
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) GSAVE(2)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO

      IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A)') ' Distance and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

10    CALL FLUSH(6)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
C           WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
C        WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         DO I=1,N
            DIAG(I)=1.0D0
         ENDDO
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR step for diagonal inverse Hessian
C
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))
C
C  Make the first guess for the step length cautious.
C
         STP=MIN(1.0D0/GNORM,GNORM)
C        STP=1.0D0
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  Update estimate of diagonal inverse Hessian elements
C  We divide by both YS and YY at different points, so
C  they had better not be zero!
C
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         IF (YY.EQ.0.0D0) YY=1.0D0
         DUMMY1=YS/YY
C        DUMMY1=ABS(YS/YY)
         DO I=1,N
            DIAG(I)= DUMMY1
         ENDDO
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  Store the new search direction
C
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF
! hk286
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) W(ISPT+POINT*N+1)=0.0D0
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) W(ISPT+POINT*N+2)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) W(ISPT+POINT*N+1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) W(ISPT+POINT*N+2)=0.0D0

C     OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'G . G=',DDOT(N,G,1,G,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reversing step'
C        IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reset'
         DO I=1,N
C           W(I)=-W(I)
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF

      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      MAXMBFGS=0.1D0
      IF (STP*SLENGTH.GT.MAXMBFGS) STP=MAXMBFGS/SLENGTH
C
C  We now have the proposed step.
C
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    ENEW=DISTFUNC(X,R,R0,R1)

      IF (ENEW-ENERGY.LE.1.0D-3) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
      ELSE 
C
C  Energy increased - try again with a smaller step size
C
C        IF (STP*SLENGTH.LT.1.0D-10) THEN
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            PRINT*,' in mmylbfgs LBFGS step cannot find a lower distance, NFAIL=',NFAIL
            ITER=0  !  try resetting
            IF (NFAIL.GT.20) THEN
               PRINT*,' Too many failures - give up'
               RETURN
            ENDIF
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
            ENDDO 
            GOTO 30
         ENDIF
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         NDECREASE=NDECREASE+1
         STP=STP/10.0D0
         IF (PTEST) 
     1    WRITE(*,'(A,F19.10,A,F16.10,A,F15.8)') ' Distance increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         GOTO 20
      ENDIF
C
C  Reset comparison geometry to R1 and rotation angles to zero. Accumulate the overall transformation matrix in 
C  OMEGATOT.
C
      OTEMP(1,1)=MYROTMAT(1,1)*OMEGATOT(1,1)+MYROTMAT(1,2)*OMEGATOT(2,1)+MYROTMAT(1,3)*OMEGATOT(3,1)
      OTEMP(1,2)=MYROTMAT(1,1)*OMEGATOT(1,2)+MYROTMAT(1,2)*OMEGATOT(2,2)+MYROTMAT(1,3)*OMEGATOT(3,2)
      OTEMP(1,3)=MYROTMAT(1,1)*OMEGATOT(1,3)+MYROTMAT(1,2)*OMEGATOT(2,3)+MYROTMAT(1,3)*OMEGATOT(3,3)
      OTEMP(2,1)=MYROTMAT(2,1)*OMEGATOT(1,1)+MYROTMAT(2,2)*OMEGATOT(2,1)+MYROTMAT(2,3)*OMEGATOT(3,1)
      OTEMP(2,2)=MYROTMAT(2,1)*OMEGATOT(1,2)+MYROTMAT(2,2)*OMEGATOT(2,2)+MYROTMAT(2,3)*OMEGATOT(3,2)
      OTEMP(2,3)=MYROTMAT(2,1)*OMEGATOT(1,3)+MYROTMAT(2,2)*OMEGATOT(2,3)+MYROTMAT(2,3)*OMEGATOT(3,3)
      OTEMP(3,1)=MYROTMAT(3,1)*OMEGATOT(1,1)+MYROTMAT(3,2)*OMEGATOT(2,1)+MYROTMAT(3,3)*OMEGATOT(3,1)
      OTEMP(3,2)=MYROTMAT(3,1)*OMEGATOT(1,2)+MYROTMAT(3,2)*OMEGATOT(2,2)+MYROTMAT(3,3)*OMEGATOT(3,2)
      OTEMP(3,3)=MYROTMAT(3,1)*OMEGATOT(1,3)+MYROTMAT(3,2)*OMEGATOT(2,3)+MYROTMAT(3,3)*OMEGATOT(3,3)
      DO J1=1,3
         DO J2=1,3
            OMEGATOT(J2,J1)=OTEMP(J2,J1)
         ENDDO
      ENDDO
      X(1)=0.0D0
      X(2)=0.0D0
      X(3)=0.0D0
      DO I=1,NSIZE
         DO J=1,3
            R0(J,I)=R1(J,I)
         ENDDO
      ENDDO

      CALL DDISTFUNC(X,GSAVE,R,R0,R1)
! hk286
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) GSAVE(1)=0.0D0
!     IF (TWOD.OR.PULLT.OR.EFIELDT.OR.GTHOMSONT) GSAVE(2)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO
      IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A,G15.5)') ' Distance and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1        ' LBFGS steps, step:',STP*SLENGTH
C
C     Compute the new step and gradient change
C
30    NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      GOTO 10

      RETURN
      END
