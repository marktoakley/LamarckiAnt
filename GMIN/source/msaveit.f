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
      SUBROUTINE MSAVEIT(EREAL,P,NP)
      USE COMMONS
      USE MODAMBER9, ONLY : IH,M04
      USE PORFUNCS
      IMPLICIT NONE

      INTEGER J1, J2, J3, J4, NP, MYUNIT2, GETUNIT, I, J, MYSTATUS
      DOUBLE PRECISION EREAL, P(3*NATOMS), LP(3*NATOMS), CMX, CMY, CMZ, MOLWT, CM(3), RANG, IT(3,3), ATMP, NEWLP(3*NATOMS)
      DOUBLE PRECISION P2(3)
      DOUBLE PRECISION VTINT, VTMIN, VTMAX
      DOUBLE PRECISION IV(3,3), LTEMP(3,3), RM(3,3)
      DOUBLE PRECISION :: ATOMCOORDS(3,3)
      CHARACTER(LEN=80) ISTR
      CHARACTER(LEN=2) ZSTRING(NATOMS)
      DOUBLE PRECISION ETheta,EPhi,EPsi,LPS(3*NATOMS)  ! sf344> Euler angles and coordinates for ellipsoids
      INTEGER :: NCALLED=0
      SAVE NCALLED

      NCALLED=NCALLED+1
      IF (.NOT.AMBERT) THEN
         ZSTRING(1:NATOMS)='X1'
         VTMIN=1.0D100
         VTMAX=-1.0D100
         IF (PYGPERIODICT.OR.PYBINARYT.OR.RIGID) NATOMS=NATOMS/2 ! temporary change!
         DO J1=1,NATOMS
            IF (VT(J1).GT.VTMAX) VTMAX=VT(J1)
            IF (VT(J1).LT.VTMIN) VTMIN=VT(J1)
         ENDDO
         VTINT=VTMAX-VTMIN
         IF (VTINT.GT.0.0D0) THEN
            DO J1=1,NATOMS
               IF (VT(J1).GT.VTMIN+2.0D0*VTINT/3.0D0) THEN
                  ZSTRING(J1)='X3'
               ELSE IF (VT(J1).GT.VTMIN+VTINT/3.0D0) THEN
                  ZSTRING(J1)='X2'
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IF (PYGPERIODICT.OR.PYBINARYT.OR.RIGID) NATOMS=NATOMS*2 
      LP(1:3*NATOMS)=P(1:3*NATOMS)
      LPS(1:3*NATOMS)=P(1:3*NATOMS)
C
C  Save the configuration to monitor.n.xyz.
C

      IF (RIGID) NATOMS=NATOMS/2
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      iloop: DO I = 1,NATOMS
         CMX = LP(3*I-2)+CMX
         CMY = LP(3*I-1)+CMY
         CMZ = LP(3*I)+CMZ
         MOLWT = MOLWT+1.0D0
      ENDDO iloop
      CM(1) = CMX/MOLWT
      CM(2) = CMY/MOLWT
      CM(3) = CMZ/MOLWT
      DO I = 1,NATOMS
         DO J=0,2
            LP(3*I-J)=LP(3*I-J)-CM(3-J)
         ENDDO
      ENDDO
C 
C     Build inertia tensor
C
      NEWLP(1:3*NATOMS)=LP(1:3*NATOMS)
      CALL INERTIA(LP,IT,NATOMS)
C
C Diagonalize inertia tensor. The 0 flag reorders the e/values and
C e/vectors from smallest to largest.
C
      CALL EIG(IT,IV,3,3,0)
C
C  Check *now* for degeneracy of eigenvalues -- If present, then see if
C  unique moment of inertia is along x.  If so, change this axis to z
C  by rotating the eigenvector matrix about y.  This guarantees that
C  the unique axis will lie along z.
C
C  The moments could be large, so this should be a different tolerance
C  to the distance checking, i.e. a relative tolerance.
C
      IF (DABS((IT(2,2)-IT(3,3))/(IT(1,1)+IT(2,2)+IT(3,3))).LT.GEOMDIFFTOL) THEN
         RANG=90.D0
         CALL ROTM(2,RANG,1,RM)
         CALL MUL3(LTEMP,IV,RM)
         DO J1=1,3
            DO J2=1,3
               IV(J2,J1)=LTEMP(J2,J1)
            ENDDO
         ENDDO
         ATMP=IT(1,1)
         IT(1,1)=IT(3,3)
         IT(3,3)=ATMP
      ENDIF
!jdf43> cartesian
      CALL MATMULV(NEWLP(1:3*NATOMS),LP(1:3*NATOMS),IV,NATOMS,3,3)
!jdf43> now rotate angle-axis coordinates
      IF (RIGID) THEN
        NATOMS=NATOMS*2
        CALL RBIMPROPERROTATION(NEWLP(:),LP(:),TRANSPOSE(IV(:,:)))
      ENDIF
!jdf43>
      IF (MOD(NCALLED,MONITORINT).EQ.0) THEN
         IF (NPAR.GT.1) THEN
            WRITE (ISTR, '(I10)') NP
            MYUNIT2=GETUNIT()
            OPEN(MYUNIT2,FILE="monitor."//TRIM(ADJUSTL(ISTR))//'.xyz',STATUS='UNKNOWN')
         ELSE
            MYUNIT2=GETUNIT()
            OPEN(MYUNIT2,FILE='monitor.xyz',STATUS='UNKNOWN')
         ENDIF
         IF (PYGPERIODICT.OR.PYBINARYT) THEN
               WRITE(MYUNIT2,*) NATOMS/2 
         ELSEIF (NTIPT) THEN
            WRITE(MYUNIT2,*) 3*NATOMS/2
         ELSE
               WRITE(MYUNIT2,*) NATOMS
         ENDIF
         WRITE(MYUNIT2,'(G20.10)') EREAL
         IF (PYGPERIODICT.OR.PYBINARYT) THEN
            DO J2=1,NATOMS/2
               CALL AAtoEuler(LPS(3*NATOMS/2+3*(J2-1)+1),LPS(3*NATOMS/2+3*(J2-1)+2),LPS(3*NATOMS/2+3*(J2-1)+3),EPhi,EPsi,ETheta)
               WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') ZSTRING(J2),LPS(3*(J2-1)+1),
     &          LPS(3*(J2-1)+2),LPS(3*(J2-1)+3),
     &         'ellipse ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EPhi,EPsi,ETheta,
     &         'atom_vector',LPS(3*NATOMS/2+3*(J2-1)+1),LPS(3*NATOMS/2+3*(J2-1)+2),LPS(3*NATOMS/2+3*(J2-1)+3)
            ENDDO
         ELSEIF (AMBERT) THEN
            WRITE(MYUNIT2,'(A2,3F20.10)') (ih(m0 4+J1-1),NEWLP(3*J1-2),NEWLP(3*J1-1),NEWLP(3*J1),J1=1,NATOMS)
         ELSEIF (NTIPT) THEN
            DO J2=1,NATOMS/2 
               J3=J2*3
               J4=(J2+NATOMS/2)*3
               CALL RIGIDTOTIP(J2,NEWLP(J3-2:J3),NEWLP(J4-2:J4),ATOMCOORDS,MYUNIT2)
            ENDDO
         ELSE
            WRITE(MYUNIT2,'(A,1X,3G20.10)') (ZSTRING(J1),NEWLP(3*(J1-1)+1),NEWLP(3*(J1-1)+2),
     &                                    NEWLP(3*(J1-1)+3),J1=1,NATOMS)
         ENDIF
         CLOSE(MYUNIT2)
      ENDIF
!
! This cp is to try and eliminate vmd load errors, which appear to occur
! if monitor.xyz is being written to by GMIN.
!
!     CALL SYSTEM_SUBR('cp monitor.xyz monitor2.xyz',MYSTATUS)
!
! Save the lowest energy encountered in best.n.xyz
!
      IF (EREAL.LT.LOWESTE) THEN
         LOWESTE=EREAL
         IF (NPAR.GT.1) THEN
            WRITE (ISTR, '(I10)') NP
            MYUNIT2=GETUNIT()
            OPEN(MYUNIT2,FILE="best."//TRIM(ADJUSTL(ISTR))//'.xyz',STATUS='UNKNOWN')
         ELSE
            MYUNIT2=GETUNIT()
            OPEN(MYUNIT2,FILE='best.xyz',STATUS='UNKNOWN')
         ENDIF
         IF (PYGPERIODICT.OR.PYBINARYT) THEN
            WRITE(MYUNIT2,*) NATOMS/2
         ELSEIF (NTIPT) THEN
            WRITE(MYUNIT2,*) 3*NATOMS/2
         ELSE
            WRITE(MYUNIT2,*) NATOMS
         ENDIF
         WRITE(MYUNIT2,'(G20.10)') EREAL
         IF (PYGPERIODICT.OR.PYBINARYT) THEN
            DO J2=1,NATOMS/2
               CALL AAtoEuler(LPS(3*NATOMS/2+3*(J2-1)+1),LPS(3*NATOMS/2+3*(J2-1)+2),LPS(3*NATOMS/2+3*(J2-1)+3),EPhi,EPsi,ETheta)
             WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') ZSTRING(J2),LPS(3*(J2-1)+1),
     &           LPS(3*(J2-1)+2),LPS(3*(J2-1)+3),
     &          'ellipse ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EPhi,EPsi,ETheta,
     &          'atom_vector',LPS(3*NATOMS/2+3*(J2-1)+1),LPS(3*NATOMS/2+3*(J2-1)+2),LPS(3*NATOMS/2+3*(J2-1)+3)
            ENDDO
         ELSEIF (AMBERT) THEN
            WRITE(MYUNIT2,'(A2,3F20.10)') (ih(m 04+J1-1),NEWLP(3*J1-2),NEWLP(3*J1-1),NEWLP(3*J1),J1=1,NATOMS)
         ELSEIF (NTIPT) THEN
            DO J2=1,NATOMS/2
               J3=J2*3
               J4=(J2+NATOMS/2)*3
               CALL RIGIDTOTIP(J2,NEWLP(J3-2:J3),NEWLP(J4-2:J4),ATOMCOORDS,MYUNIT2)
            ENDDO
         ELSE
            WRITE(MYUNIT2,'(A,1X,3G20.10)') (ZSTRING(J1),NEWLP(3*(J1-1)+1),NEWLP(3*(J1-1)+2),
     &                                       NEWLP(3*(J1-1)+3),J1=1,NATOMS)
         ENDIF
         CLOSE(MYUNIT2)
      ENDIF
C
C Dump quench number, energy, and Markov energy in file lowest.n.e.
C
C     IF (NPAR.GT.1) THEN
C        WRITE (ISTR, '(I10)') NP
C        MYUNIT2=GETUNIT()
C        OPEN(MYUNIT2,FILE="monitor."//TRIM(ADJUSTL(ISTR))//'.e',STATUS='UNKNOWN',POSITION='APPEND')
C     ELSE
C        MYUNIT2=GETUNIT()
C        OPEN(MYUNIT2,FILE='monitor.e',STATUS='UNKNOWN',POSITION='APPEND')
C     ENDIF

C     WRITE(MYUNIT2,'(I10,2F25.15)') NQ(NP),EREAL,EPREV(NP)
C     CLOSE(MYUNIT2)

      RETURN
      END

      SUBROUTINE MUL3(TEMP,RM,ROT)
      IMPLICIT NONE

      DOUBLE PRECISION TEMP(3,3), RM(3,3), ROT(3,3)

      TEMP(1,1)=RM(1,1)*ROT(1,1)+RM(1,2)*ROT(2,1)+RM(1,3)*ROT(3,1)
      TEMP(2,1)=RM(2,1)*ROT(1,1)+RM(2,2)*ROT(2,1)+RM(2,3)*ROT(3,1)
      TEMP(3,1)=RM(3,1)*ROT(1,1)+RM(3,2)*ROT(2,1)+RM(3,3)*ROT(3,1)
      TEMP(1,2)=RM(1,1)*ROT(1,2)+RM(1,2)*ROT(2,2)+RM(1,3)*ROT(3,2)
      TEMP(2,2)=RM(2,1)*ROT(1,2)+RM(2,2)*ROT(2,2)+RM(2,3)*ROT(3,2)
      TEMP(3,2)=RM(3,1)*ROT(1,2)+RM(3,2)*ROT(2,2)+RM(3,3)*ROT(3,2)
      TEMP(1,3)=RM(1,1)*ROT(1,3)+RM(1,2)*ROT(2,3)+RM(1,3)*ROT(3,3)
      TEMP(2,3)=RM(2,1)*ROT(1,3)+RM(2,2)*ROT(2,3)+RM(2,3)*ROT(3,3)
      TEMP(3,3)=RM(3,1)*ROT(1,3)+RM(3,2)*ROT(2,3)+RM(3,3)*ROT(3,3)

      RETURN
      END

