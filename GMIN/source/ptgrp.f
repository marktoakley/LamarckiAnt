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
C  SYMMETRY ANALYSIS - save the generators of symmetry
C  operations detected for a subset of atoms (or all of them).
C  To apply these later we need all atoms in the corresponding
C  orientation. NATOMS is the number of atoms at the beginning of
C  Q whose symmetry elements we seek.
C
      SUBROUTINE PTGRP(Q,NATOMS,DEBUG,SYMTOL1,SYMTOL2,SYMTOL3,GENMAT,IGEN,FPGRP,CM,MATDIFF)
      USE COMMONS,ONLY : MYUNIT,HPTGRP,ATMASS
      IMPLICIT NONE
      INTEGER NATOMS
      DOUBLE PRECISION IT(3,3),CM(3),IV(3,3),TOLD,QREF(3*NATOMS),DIST2,MOLWT,
     1                 RM(3,3),QSORT(3*NATOMS),TATB(3), XTEMP(3,3), RMAT(3,3), DUM(3,3), DUM2(3,3),
     2                 SCRATCH(9*NATOMS), Q(3*NATOMS), GENMAT(100,3,3), 
     3                 MATDIFF, DISTANCE(120), QSAVE(3*NATOMS),
     4                 ANG,ANGLE2,ARGU,CMX,CMY,CMZ,RMINV(3,3),DUMMY,DUMMYSINV,WORSTRAD
      DOUBLE PRECISION NEWQ(3*NATOMS), MYDOT, ANGMAG, ANGMG, ANGIAX, DIP, TEST, BILEN, DIST, ZANG2, ZANGS, SYMTOL2, SYMTOL3, TOLO
      DOUBLE PRECISION ARG, ANGL, DPROJ, ZANG, Z, X, ATMP, XXX, RANG, TOLE, DTOR, ORDIS, SYMTOL1
      DOUBLE PRECISION CX, CY, CZ
      INTEGER IGEN, JPLANE, J3, IAXORD, ICLIP, INDEX, LOOK, IBOT2, ISIGMV, IMNAX, ISTART
      INTEGER MINORB, IK, IORBIT, NORBIT, IBOT, IDIR, IXX, IDEGEN, J2, J1, J, I, ISKIP, ICOMPX, JCOMP
      INTEGER JCOMPX, IHIGHX, ISAXIS, ICOMP2, NUMB, ISET, IERR, IIAX, ICOMPQ, ISINV
      INTEGER IDONE, ICOUNT, ICOMP, ILINEAR, IDGRP, IINV, IREF, IROT, IHIGH, IPRNT, IBTOR, PERM(NATOMS), IDUMMY
      LOGICAL NEW, DEBUG
C
C     Symmetry Information
C     FPGRP   Full point group
C     BPGRP   Largest Abelian subgroup
C     PGRP    "Computational" point group
C
      CHARACTER(LEN=4) FPGRP, BPGRP, PGRP
C     INTEGER NORD(2*NATOMS+1)
      INTEGER NORD2(2*NATOMS+1)
      LOGICAL AGAIN
      CHARACTER(LEN=4) STRING, JNKSTR
      CHARACTER GPSTRING*80
 
      IBTOR(I,J) =IOR(I,J)
 
      IPRNT=0
      IF (DEBUG) IPRNT=3
      HPTGRP=1
      FPGRP='   '
      DTOR=ACOS(-1.D0)/180.D0
      AGAIN=.FALSE.
      TOLO=SYMTOL1
      TOLD=SYMTOL2
      TOLE=SYMTOL3
      QSAVE(1:3*NATOMS)=Q(1:3*NATOMS)
C
C Initialize.
C
651   Q(1:3*NATOMS)=QSAVE(1:3*NATOMS) ! otherwise the CofM can be reset differently 
      IHIGH=0
      IGEN=0
      IROT=0
      IREF=0
      IINV=0
      IDGRP=0
      ILINEAR=0
C     ICOMP=0
      ICOUNT=0
      IDONE=0
C     ISINV=0
C     ICOMPQ=0
C     IIAX=0
      ORDIS=0
      IERR=0
      ISET=0
      NUMB=0
C     ICOMP2=0
      IHIGHX=0
      ISAXIS=0
      JCOMPX=0
      JCOMP=0
C     ICOMPX=0
      ISKIP=0
C
C  Translate to center of mass.
C
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      DO I=1,NATOMS
         CMX=Q(3*I-2)*ATMASS(I)+CMX
         CMY=Q(3*I-1)*ATMASS(I)+CMY
         CMZ=Q(3*I)*ATMASS(I)+CMZ
         MOLWT=MOLWT+ATMASS(I)
      ENDDO
      CM(1)=CMX/MOLWT
      CM(2)=CMY/MOLWT
      CM(3)=CMZ/MOLWT
Cds656> testing...
      IF(DEBUG) THEN
         CX=0.D0
         CY=0.D0
         CZ=0.D0
         DO I=1,NATOMS
            CX=Q(3*I-2)+CX
            CY=Q(3*I-1)+CY
            CZ=Q(3*I)+CZ
         ENDDO
         WRITE(MYUNIT,'(A,3(F10.5))')'ptgrp> CoG:',
     2        CX/DBLE(NATOMS), CY/DBLE(NATOMS), CZ/DBLE(NATOMS)
         WRITE(MYUNIT,'(A,3(F10.5))')'ptgrp> CoM:',(CM(I),I=1,3)
      ENDIF
C<ds656 ...testing
      DO I=1,NATOMS
         DO J=0,2
            Q(3*I-J)=Q(3*I-J)-CM(3-J)
         ENDDO
      ENDDO
      QREF(1:3*NATOMS)=Q(1:3*NATOMS)
50    FORMAT(3(2X,F12.6))
C
C     Build inertia tensor
C
      CALL INERTIA(Q,IT,NATOMS)
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
      IF (ABS((IT(2,2)-IT(3,3))/(IT(1,1)+IT(2,2)+IT(3,3))).LT.TOLE) THEN
         RANG=90.0D0
         CALL ROTM(2,RANG,1,RM)
         CALL MATMUL(XTEMP,IV,RM,3,3,3,3,3,3)
         DO J1=1,3
            DO J2=1,3
               IV(J2,J1)=XTEMP(J2,J1)
            ENDDO
         ENDDO
         ATMP=IT(1,1)
         IT(1,1)=IT(3,3)
         IT(3,3)=ATMP
      ENDIF

      IF (IPRNT.GE.3) WRITE(MYUNIT,*)' Diagonalized inertia tensor '
      IF (IPRNT.GE.3) WRITE(MYUNIT,50) ((IT(I,J),J=1,3),I = 1,3)
      CALL MATMULV(NEWQ,Q,IV,NATOMS,3,3)
      IF (IPRNT.GE.4) WRITE(MYUNIT,*) ' Principal axis orientation for molecular system '
      IF (IPRNT.GE.4) WRITE(MYUNIT,50) (NEWQ(I),I=1,3*NATOMS)
      AGAIN=.TRUE.
C
C  Check handedness of inertial axes.
C  If dot product is negative; switch sign of
C  y axis (the choice of axis to change is arbitrary).
C
      CALL MYCROSS(IV(1,1),IV(1,2),TATB,1)
      XXX=MYDOT(TATB,IV(1,3),3)
      IF (XXX.LT.0.D0) THEN
         DO I=2,3*NATOMS-1,3
            NEWQ(I)=-NEWQ(I)
         ENDDO
      ENDIF
C
C  The geometry used to find the
C  symmetry operations changes orientation throughout ptgrp. Use new subroutine reorient to
C  find the matrix that maps the reference configuration in QREF onto the current
C  configuration, and get the generator for the reference configuration from a
C  suitable matrix transformation.
C
C  Generate well-defined sorted coordinate vector now. Used in
C  symmetry evaluation routines which follow.
C  If QSORT changes then the reference orientation has changed, and
C  subsequent generators will be with respect to the new orientation.
C
C     CALL SORTXYZ(NEWQ,QSORT,NORD,TOLD,NATOMS)
      QSORT(1:3*NATOMS)=NEWQ(1:3*NATOMS)
      IF (IPRNT.GE.4) WRITE(MYUNIT,*) ' SORTED COORDINATE VECTOR'
      IF (IPRNT.GE.4) WRITE(MYUNIT,50) (QSORT(J),J=1,NATOMS*3)
C
C  Check to see if point group is Abelian - if so, skip to the Abelian
C  operations block.
C
      X=-1.D0
      IDEGEN=0
      DO I=1,3
         Z=IT(I,I)
         IF (ABS((Z-X)/(IT(1,1)+IT(2,2)+IT(3,3))).LT.TOLE) IDEGEN=IDEGEN+1 
!        WRITE(MYUNIT,'(A,I6,4G20.10,I5)') 'I,X,Z,test,tole,IDEGEN=',I,X,Z,ABS((Z-X)/(IT(1,1)+IT(2,2)+IT(3,3))),TOLE,IDEGEN
         X=Z
      ENDDO
      IF (IDEGEN.EQ.0) THEN
         IF (IPRNT.GE.3) WRITE(MYUNIT,'(A)') ' The molecule belongs to an Abelian group.'
         GOTO 630
      ENDIF
C
C  *** Start of the big IDEGEN=1 block.  ***
C
      IF (IDEGEN.EQ.1) THEN
C
C  Check to make sure that the principal axes
C  are parallel to x, y and z. If not, rotate orientation to this point.
C
C  Use a different tolerance for zero moment of inertia checking.
C
         FPGRP='   '
         IF (IPRNT.GE.3) WRITE(MYUNIT,110)
110      FORMAT(' The molecule belongs to a point group with doubly',
     1        ' degenerate representations.')
         DO I=1,3
            IF (ABS(IT(I,I)).LT.5.0D-2) ILINEAR=1
         ENDDO
         IF (ILINEAR.EQ.1) THEN
            IF (IPRNT.GE.3) WRITE(MYUNIT,'(A)') ' The molecule is linear.'
            GOTO 630
         ENDIF
C
C  Check for rotational axes from C2 to 6. 
C
         DO I=2,6
            ZANG=360.D0/(1.0D0*I)
            DO J=1,3
               CALL ROTM(J,ZANG,1,RM)
               CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
C              CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C              CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
               CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C              CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
C              IF (IPRNT.GE.4) WRITE(MYUNIT,*) ICOMP
C              IF (ICOMP.EQ.0) THEN
               IF (IPRNT.GE.4) WRITE(MYUNIT,'(A,I5,3F12.3)') 'IAXORD,DUMMY,DIST2,TOLD=',I,DUMMY,DIST2,TOLD
               IF (DIST2.LT.TOLD) THEN
                  IHIGH=I
                  IHIGHX=J
                  IGEN=1 ! just save the generator for the principal axis at this point
C                  IF (IGEN.LT.100) THEN
C                     DO J2=1,IGEN
C                        CALL CMPMAT(GENMAT,100,J2,RM,NEW,MATDIFF)
C                        IF (.NOT.NEW) GOTO 10
C                     ENDDO
C                     IGEN=IGEN+1
CC                    PRINT*,'here A IGEN=',IGEN
                     IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'axis ',J,' is a rotation axis of order ',I
C
C   Make the transformation matrix corresponding to the reference state in QREF
C
!                    IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient A'
                     CALL REORIENT(NATOMS,SCRATCH,QREF,RMAT)
                     CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
                     CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
                     GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!                    IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here A IGEN=',IGEN,' new GENMAT:'
!                    IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
10                   CONTINUE
C                  ENDIF
C
C  Stuff ordering into NORD2
C
C                 DO IXX=1,2*NATOMS
C                    NORD2(IXX)=NORD(IXX)
C                 ENDDO
                  DO IXX=1,NATOMS
                     NORD2(IXX)=IXX
                     NORD2(IXX+NATOMS)=PERM(IXX)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
C
C  Catch problem if IHIGH=0 because two of the inertia tensor
C  e/values just happen to be very close - start all over again
C  with smaller TOLE.
C
         IF (IHIGH.EQ.0) THEN
            IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'Accidental degeneracy detected'
            TOLE=TOLE/10.0D0
            IF (TOLE.GT.1.0D-7) THEN
               GOTO 651
            ELSE
               IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) ' The full molecular point group is undetermined'
               WRITE(GPSTRING,653) ' The full molecular point group is undetermined'
653            FORMAT(A80)
               GOTO 652
            ENDIF
         ENDIF

         IF (IPRNT.GE.3) WRITE(MYUNIT,160) IHIGH,IHIGHX
160      FORMAT(' The highest order rotational axis is C',I2,' about ',I2)
C
C  If highest order axis is 2, and there is twofold degeneracy, the
C  group must be D2d. Proceed to find the S4 axis which determines
C  the unique rotational axis. This is subsequently used to rotate
C  molecule to a useful orientation. (NORD2 now contains an effective
C  permutation list).
C
         IF (IHIGH.EQ.2) THEN
            IHIGHX=0
            FPGRP='D2d'
            HPTGRP=8
            ZANG=90.0D0
            DO IDIR=1,3
               CALL ROTM(IDIR,ZANG,1,RM)
               CALL MATMULV(SCRATCH(NATOMS*3+1),NEWQ,RM,NATOMS,3,3)
               CALL REFLECT(SCRATCH(NATOMS*3+1),SCRATCH,NATOMS,IDIR)
C              CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C              IF (ICOMP.EQ.0) THEN
C              CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
               CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C              CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
               IF (DIST2.LT.TOLD) THEN
                  IHIGHX=IDIR
C                 DO I=1,2*NATOMS
C                    NORD2(I)=NORD(I)
C                 ENDDO
                  DO I=1,NATOMS
                     NORD2(I)=I
                     NORD2(I+NATOMS)=PERM(I)
                  ENDDO
               ENDIF
            ENDDO
            IF (IHIGHX.EQ.0) THEN
               IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'Accidental degeneracy detected'
               TOLE=TOLE/10.0D0
               IF (TOLE.GT.1.0D-7) THEN
                  GOTO 651
               ELSE
                  IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) ' The full molecular point group is undetermined'
                  WRITE(GPSTRING,653) ' The full molecular point group is undetermined'
                  GOTO 652
               ENDIF
            ENDIF
            IF (IPRNT.GE.3) WRITE(MYUNIT,190) IHIGHX
190         FORMAT(' The unique rotational axis is ',I1,'.')
         ENDIF
C
C  Now, the unique axis is known, no matter what the group is.
C  Proceed to rotate around unique axis such that the projection
C  of one of the atoms which are equivalent under the rotation lies
C  along one of the other two Cartesian axes.
C
C  Find an atom which is permuted under the rotation.
C  New problem - in a large molecule it is possible to have, say,
C  a 10-orbit, with a principal axis of order only 5. In this case
C  we need to find an atom from an orbit of order equal to the
C  order of the principal axis to look for vertical mirror planes.
C  If we put an atom from the 10-orbit in the xz or yz plane we
C  have not put a vertical mirror plane coincident with the axis.
C  On the other hand, for D_{nd} groups with odd n we actually
C  may not have an orbit of size equal to the order of the principal
C  axis - there is sure to be one of twice the size though. Here
C  the original algorithm would have worked.
C  The principal axis order is IHIGH - this is only needed if
C  IHIGH is odd.
C
         IF (MOD(IHIGH,2).EQ.1) THEN
            IF (IPRNT.GE.4) WRITE(MYUNIT,*)' Before rotation'
            IF (IPRNT.GE.4) WRITE(MYUNIT,50) (NEWQ(J),J=1,NATOMS*3)
            CALL ZERO(SCRATCH,3*NATOMS*3)
            DO J=1,NATOMS
               NORD2(J)=J
               IBOT=3*J-2
               SCRATCH(J)=ATMASS(J)*SQRT(MYDOT(NEWQ(IBOT),NEWQ(IBOT),3))
            ENDDO
            CALL PIKSR2(NATOMS,SCRATCH,NORD2)
            IF (IPRNT.GE.3) WRITE(MYUNIT,440) (SCRATCH(J),J=1,NATOMS)
440         FORMAT(' ORBIT MATRIX:',(F10.6,/))
            NORBIT=1
            IORBIT=1
            DO J=2,NATOMS+1
C              IF (ABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLD) IORBIT=IORBIT+1
C              IF (ABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLD) THEN
               IF (ABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLO) IORBIT=IORBIT+1
               IF (ABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLO) THEN
                  NORD2(NATOMS+NORBIT)=IORBIT
                  SCRATCH(6*NATOMS+NORBIT)=ABS(SCRATCH(J-1))/ATMASS(NORD2(J-1))
                  IF (J.LE.NATOMS) THEN
                     NORBIT=NORBIT+1
                     IORBIT=1
                  ENDIF
               ENDIF
            ENDDO

            IF (IPRNT.GE.130) THEN
                WRITE(MYUNIT,*)' SCRATCH VECTOR '
                WRITE(MYUNIT,*) (SCRATCH(IK),IK=1,NATOMS)
                WRITE(MYUNIT,*)' ORBIT DISTANCES: '
                WRITE(MYUNIT,*) (SCRATCH(6*NATOMS+IK),IK=1,NORBIT)
            ENDIF
            NUMB=0
            IREF=0
C
C  Here we are trying to find an orbit whose size NORD2 is a multiple
C  of IHIGH. If such an orbit is not found then IREF will remain zero 
C  and we must find more self-consistent tolerances.
C
            DO I=1,NORBIT
               NUMB=NORD2(NATOMS+I)+NUMB
C              WRITE(MYUNIT,*) 'I,IHIGH,NORD2=',I,IHIGH,NORD2(NATOMS+I)
               IF (NORD2(NATOMS+I).EQ.IHIGH) THEN
                  MINORB=NORD2(NATOMS+I)
                  IREF=1+NUMB-NORD2(NATOMS+I)
                  GOTO 250
               ENDIF
            ENDDO
250         IF (IREF.NE.0) THEN
               IF (IPRNT.GE.13) WRITE(MYUNIT,260) MINORB, NORD2(IREF) 
260            FORMAT(' The principal axis has order ',i3,/,
     1         ' Atom number ',i3,' belongs to an orbit of this order and', 
     2         ' will be used as a reference.')
               IF (IPRNT.GE.13) WRITE(MYUNIT,520) (NORD2(IREF+J),J=1,MINORB-1)
               ICOUNT=NORD2(IREF)
            ELSE
               NUMB=0
               IREF=0
               DO I=1,NORBIT
                  NUMB=NORD2(NATOMS+I)+NUMB
                  IF (NORD2(NATOMS+I).EQ.2*IHIGH) THEN
                     MINORB=NORD2(NATOMS+I)
                     IREF=1+NUMB-NORD2(NATOMS+I)
                     GOTO 255
                  ENDIF
               ENDDO
255            IF (IREF.NE.0) THEN
                  IF (IPRNT.GE.13) WRITE(MYUNIT,265) MINORB, NORD2(IREF)
265               FORMAT(' The principal axis has order ',i3,/,
     1                   ' Atom number ',i3,' belongs to an orbit of twice this order', 
     2                   ' and will be used as a reference.')
                  IF (IPRNT.GE.13) WRITE(MYUNIT,520) (NORD2(IREF+J),J=1,MINORB-1)
                  ICOUNT=NORD2(IREF)
               ELSE
                  NUMB=0
                  IREF=0
                  DO I=1,NORBIT
                     NUMB=NORD2(NATOMS+I)+NUMB
                     IF (NORD2(NATOMS+I).EQ.4*IHIGH) THEN
                        MINORB=NORD2(NATOMS+I)
                        IREF=1+NUMB-NORD2(NATOMS+I)
                        GOTO 256
                     ENDIF
                  ENDDO
256               IF (IREF.EQ.0) THEN
                     IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(A,G20.10)')
     &   'Orbit size inconsistent with axis order - try decreasing TOLD to ',TOLD/10.0D0
                     TOLD=TOLD/10.0D0
                     IF (TOLD.GT.1.0D-7) THEN
                        GOTO 651
                     ELSE
                        IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) ' The full molecular point group is undetermined'
                        WRITE(GPSTRING,653) ' The full molecular point group is undetermined'
                        GOTO 652
                     ENDIF
                  ENDIF
                  IF (IPRNT.GE.13) WRITE(MYUNIT,266) MINORB, NORD2(IREF)
266               FORMAT(' The principal axis has order ',i3,/,
     1' Atom number ',i3,' belongs to an orbit of four times this order'
     2            ,' and will be used as a reference.')
                  IF (IPRNT.GE.13) WRITE(MYUNIT,520) (NORD2(IREF+J),J=1,MINORB-1)
                  ICOUNT=NORD2(IREF)
               ENDIF
            ENDIF
         ELSE
            DO I=NATOMS,1,-1
C              PRINT*,'I,ICOUNT,NORD2,NORD2(I+NATOMS)=',I,ICOUNT,NORD2(I),NORD2(I+NATOMS)
C              IF (NORD2(I).NE.NORD2(I+NATOMS)) ICOUNT=NORD2(I)
               IF (I.NE.PERM(I)) ICOUNT=I
C              PRINT '(A,3I5)','I,PERM(I),ICOUNT=',I,PERM(I),ICOUNT
            ENDDO
         ENDIF
C
C   Calculate angle between projection of this vector and one of the 
C   Cartesian axes
C
         ISTART=3*(ICOUNT-1)
         IF (ICOUNT.LT.1) THEN
            WRITE(MYUNIT,*) 'ERROR - ICOUNT=',ICOUNT
         ENDIF
         DPROJ=SQRT( NEWQ(ISTART+1)**2 + NEWQ(ISTART+2)**2 +
     1        NEWQ(ISTART+3)**2 - NEWQ(ISTART+IHIGHX)**2 )
         IF (DPROJ.EQ.0.0D0) THEN
            WRITE(MYUNIT,'(A)') 'WARNING in ptgrp, zero projection, this should never happen'
            DPROJ=1.0D-3
         ENDIF
C
C   Find first Cartesian axis which is not the unique axis.
C   (Has to be x or y)
C
         IMNAX=2-IHIGHX/2
         ARG=NEWQ(ISTART+IMNAX)/DPROJ
         IF (ABS(ARG).GT.1.0D0) ARG=1.0D0*ARG/ABS(ARG)
         ANGL=ACOS(ARG)/DTOR
C
C   Skip out of loop if alignment already satisfied
C
         DO I=0,2
            X=ANGL-1.0D0*(I*90)
C           IF (ABS(X).LT.TOLD) ISKIP=1
            IF (ABS(X).LT.TOLO) ISKIP=1
         ENDDO
         IF (ISKIP.EQ.1) GOTO 340
C
C   Get sign of angle needed to rotate molecule into position
C
C   I found a case where this didn't work. There are only two
C   possibilities (+ and -) so if one doesn't work, just use 
C   the other!
C
C        Z0=1
C        DO I=1,3
C           IF (I.NE.IHIGHX) Z0=Z0*NEWQ(ISTART+I)
C        ENDDO
C        IF (Z0.GT.0.D0) ANGL=-ANGL
C
C   Rotate molecule into position - all NATOMS 
C
         CALL ZERO(SCRATCH,NATOMS*2*3)
         CALL VADD(SCRATCH(1),SCRATCH(NATOMS*3+1),NEWQ,NATOMS*3,1) ! puts NEWQ in SCRATCH
         CALL ROTM(IHIGHX,ANGL,1,RM)
         CALL ROTM(IHIGHX,-ANGL,1,RMINV)
C        PRINT*,'NEWQ before transformation:'
C        PRINT '(3G20.10)',NEWQ(1:3*NATOMS)
         CALL MATMULV(NEWQ,SCRATCH,RM,NATOMS,3,3) ! puts rotated structure in NEWQ
         IF (IPRNT.GE.3) WRITE(MYUNIT,310) ANGL,IHIGHX
310      FORMAT(' Molecule rotated through ',F20.10,' degrees ',
     1        'about the ',i1,' axis.')
         IF (IPRNT.GE.4) WRITE(MYUNIT,*)' After rotation:'
         IF (IPRNT.GE.4) WRITE(MYUNIT,50) (NEWQ(J),J=1,NATOMS*3)

         DO 320 J1=1,3
            IF (J1.EQ.IHIGHX) GOTO 320
            IF (J1.EQ.IMNAX) GOTO 320
            IF (ABS(NEWQ(ISTART+J1)).GT.1.0D-6) THEN
               ANGL=-2.0D0*ANGL
               CALL ZERO(SCRATCH,NATOMS*2*3)
               CALL VADD(SCRATCH(1),SCRATCH(NATOMS*3+1),NEWQ,NATOMS*3,1) ! puts NEWQ in SCRATCH
               CALL ROTM(IHIGHX,ANGL,1,RM)
               CALL ROTM(IHIGHX,-ANGL,1,RMINV)
               CALL MATMULV(NEWQ,SCRATCH,RM,NATOMS,3,3) ! puts transformed structure in NEWQ
               IF (IPRNT.GE.3) WRITE(MYUNIT,310) ANGL,IHIGHX
               IF (IPRNT.GE.4) WRITE(MYUNIT,*)' After rotation:'
               IF (IPRNT.GE.4) WRITE(MYUNIT,50) (NEWQ(J),J=1,NATOMS*3)  
            ENDIF
320      CONTINUE
C
C   New reference structure, so we need to generate new sorted vector
C   And update the transformation matrix for the reference structure
C
C        CALL SORTXYZ(NEWQ,QSORT,NORD,TOLD,NATOMS)
         QSORT(1:3*NATOMS)=NEWQ(1:3*NATOMS)
C     PRINT*,'B qsort:'
C     PRINT '(3G20.10)',QSORT(1:3*NATOMS)
C     PRINT*,'newq:'
C     PRINT '(3G20.10)',NEWQ(1:3*NATOMS)
         IF (IPRNT.GE.4) WRITE(MYUNIT,*) ' New sorted coordinates'
         IF (IPRNT.GE.4) WRITE(MYUNIT,50) (QSORT(J),J=1,NATOMS*3)
C
C   Check for sigma(v) planes now.
C
         CALL ZERO(SCRATCH,NATOMS*3*2)
         ISIGMV=0
         DO 330 J=1,3
            IF (J.EQ.IHIGHX) GOTO 330
            CALL REFLECT(NEWQ,SCRATCH,NATOMS,J)
C           CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMPX,TOLD,NATOMS,IPRNT)
C           IF (ICOMPX.EQ.0) ISIGMV=1
C           CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
            CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C           CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
            IF (DIST2.LT.TOLD) ISIGMV=1
330      CONTINUE
C
C   Check for S(2n) axis
C
340      ZANGS=180.D0/(1.0D0*IHIGH)
         CALL ROTM(IHIGHX,ZANGS,1,RM)
         CALL MATMULV(SCRATCH(NATOMS*3+1),NEWQ,RM,NATOMS,3,3)
         CALL REFLECT(SCRATCH(NATOMS*3+1),SCRATCH,NATOMS,IHIGHX)
C        CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C        IF (ICOMP.EQ.0) THEN
C        CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
         CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C        PRINT '(A,G20.10,A,I6,A,G20.10,I6)','DIST2 after check for S(2n) axis=',DIST2,' IHIGH=',IHIGH,' TOLD,IPRNT=',TOLD,IPRNT
C        CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
         IF (DIST2.LT.TOLD) THEN
            ISAXIS=1
            IF (IPRNT.GE.3) THEN
               WRITE(MYUNIT,350) 2*IHIGH
350            FORMAT(' S',I3,' axis exists.')
            ENDIF
            IF (IGEN.LT.100) THEN
               DO J1=1,3
                  RM(J1,IHIGHX)=-RM(J1,IHIGHX)
               ENDDO
C
C   Make the transformation matrix corresponding to the reference state in QREF 
C
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient F'
               CALL REORIENT(NATOMS,SCRATCH,QREF,RMAT) 
               CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
               CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
               DO J2=1,IGEN
                  CALL CMPMAT(GENMAT,100,J2,DUM2,NEW,MATDIFF)
                  IF (.NOT.NEW) GOTO 61
               ENDDO
               IGEN=IGEN+1
               GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here F IGEN=',IGEN,' new GENMAT:'
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
61          CONTINUE
            ENDIF
         ENDIF
C        PRINT*,'dummy after check for S(2n) axis = ',DUMMY,' ISAXIS=',ISAXIS
C
C   Look for perpendicular C2 axes - check along Cartesian axes and
C   along angle bisectors of S2n operations
C
C   first rotate molecule by 1/2 S(2n) angle around IHIGH x and put in
C   SCRATCH(6*NATOMS+1)
C
         ZANG=180.D0
         ZANG2=180.D0/(2.D0*IHIGH)
         CALL ROTM(IHIGHX,ZANG2,1,RM)
         CALL MATMULV(SCRATCH(6*NATOMS+1),NEWQ,RM,NATOMS,3,3)
         CALL ZERO(SCRATCH(3*NATOMS+1),3*NATOMS)
C        CALL SORTXYZ(SCRATCH(6*NATOMS+1),SCRATCH(3*NATOMS+1),NORD,TOLD,NATOMS)
         SCRATCH(3*NATOMS+1:6*NATOMS)=SCRATCH(6*NATOMS+1:9*NATOMS)
         DO 360 I=1,3
            IF (I.EQ.IHIGHX) GOTO 360
C
C   Checking Cartesian axes for both rotated and unrotated molecule
C
            CALL ROTM(I,ZANG,1,RM)
            CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
C           CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C           CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
            CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C           CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
            IF (DIST2.LT.TOLD) JCOMP=JCOMP+1
            CALL MATMULV(SCRATCH,SCRATCH(6*NATOMS+1),RM,NATOMS,3,3)
C           CALL COMPARE2(SCRATCH,SCRATCH(3*NATOMS+1),NORD,ICOMPX,TOLD,NATOMS,IPRNT)
C           CALL MINPERM(NATOMS,SCRATCH(3*NATOMS+1),SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
            CALL TESTSYMOP(NATOMS,SCRATCH(3*NATOMS+1),SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C           CALL BIPARTITE(NATOMS,SCRATCH(3*NATOMS+1),SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
            IF (DIST2.LT.TOLD) JCOMPX=JCOMPX+1
C           IF (ICOMP.EQ.0) JCOMP=JCOMP+1
C           IF (ICOMPX.EQ.0) JCOMPX=JCOMPX+1
360      CONTINUE
         IF (JCOMP.GE.1) THEN
            IF (IPRNT.GE.3) WRITE(MYUNIT,'(A)') ' On-axis perpendicular C2 elements found.'
            IDGRP=1
         ELSEIF (JCOMPX.GE.1) THEN
            IF (IPRNT.GE.3) WRITE(MYUNIT,'(A)') ' Off-axis perpendicular C2 elements found.'
            IDGRP=1
         ENDIF
C
C   Check for sigma(h)
C
         IF (IDGRP.EQ.1) THEN
            CALL REFLECT(NEWQ,SCRATCH,NATOMS,IHIGHX)
C           CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP2,TOLD,NATOMS,IPRNT)
C           IF (ICOMP2.EQ.0) THEN
C           CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
            CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C           CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
            IF (DIST2.LT.TOLD) THEN
               FPGRP='DNh'
               HPTGRP=4*IHIGH
               IF (IPRNT.GE.3) WRITE(MYUNIT,390)
            ENDIF
         ELSE
            CALL REFLECT(NEWQ,SCRATCH,NATOMS,IHIGHX)
390         FORMAT(' Horizontal plane of symmetry found.')
C           CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP2,TOLD,NATOMS,IPRNT)
C           IF (ICOMP2.EQ.0) THEN
C           CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
            CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C           CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
            IF (DIST2.LT.TOLD) THEN
               FPGRP='CNh'
               HPTGRP=2*IHIGH
               IF (IPRNT.GE.3) WRITE(MYUNIT,390)
            ELSE
               FPGRP='CN '
               HPTGRP=IHIGH
               IF (ISAXIS.EQ.1) THEN
                  FPGRP='SN '
                  IHIGH=2*IHIGH
                  HPTGRP=IHIGH
               ENDIF
               DO 400 J=1,3
                  IF (J.EQ.IHIGHX) GOTO 400
                  CALL REFLECT(NEWQ,SCRATCH,NATOMS,J)
C                 CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP2,TOLD,NATOMS,IPRNT)
C                 IF (ICOMP2.EQ.0) FPGRP='CNv'
C                 CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
                  CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C                 CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
                  IF (DIST2.LT.TOLD) THEN
                     FPGRP='CNv'
                     HPTGRP=2*IHIGH
                  ENDIF
400            CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
C  *** End of IDEGEN=1 block. *** 
C
      IF (IDEGEN.EQ.2) THEN
         IF (IPRNT.GE.3) WRITE(MYUNIT,*) ' In cubic loop '
         IF (IPRNT.GE.3) WRITE(MYUNIT,410)
410      FORMAT(' The molecule belongs to a point group with doubly', 
     1          ' and triply degenerate representations.')
C
C   Now we know that the point group is either T, Td, O, Oh, T, I, Th
C   or Ih.  For these groups, however, the inertia tensor
C   is rotationally invariant; therefore different (and significantly
C   more complicated) methods need to be used.  This one seems to
C   be particularly efficient.
C
C   Check for inversion symmetry. This will tell us a great deal.
C   (If inversion center present, then group is either Oh, Th or Ih.
C
         DO I=1,NATOMS*3
            SCRATCH(I)=-NEWQ(I)
         ENDDO
C        CALL COMPARE2(SCRATCH,QSORT,NORD,ISINV,TOLD,NATOMS,IPRNT)
C        CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DUMMYSINV,WORSTRAD)
         CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DUMMYSINV,WORSTRAD)
C        CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DUMMYSINV,WORSTRAD)
C
C   Locate all orbits (groups of equivalent atoms at distance r) -
C   do this by generating and then sorting the c.o.m-outer atom
C   distance matrix (mass-weighted). Use NORD2 to keep track of
C   which atoms are in which orbit.
C
         CALL ZERO(SCRATCH,3*NATOMS*3)
         DO J=1,NATOMS
            NORD2(J)=J
            IBOT=3*J-2
            SCRATCH(J)=ATMASS(J)*SQRT(MYDOT(NEWQ(IBOT),NEWQ(IBOT),3))
         ENDDO
         CALL PIKSR2(NATOMS,SCRATCH,NORD2)
         IF (IPRNT.GE.3) WRITE(MYUNIT,440) (SCRATCH(J),J=1,NATOMS)
C
C   Now count number of orbits; place pertinent info (distances) in top
C   end of SCRATCH array and use top end of NORD2 to hold number
C   of centers/orbit.  This bookkeeping facilitates things down
C   the road.  
C
         NORBIT=1
         IORBIT=1
         DO J=2,NATOMS+1
C           IF (ABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLD) IORBIT=IORBIT+1
C           IF (ABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLD) THEN
            IF (ABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLO) IORBIT=IORBIT+1
            IF (ABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLO) THEN
               NORD2(NATOMS+NORBIT)=IORBIT
               SCRATCH(6*NATOMS+NORBIT)=ABS(SCRATCH(J-1))/ATMASS(NORD2(J-1))
               IF (J.LE.NATOMS) THEN
                  NORBIT=NORBIT+1
                  IORBIT=1
               ENDIF
            ENDIF
         ENDDO

         IF (IPRNT.GE.130) THEN
            WRITE(MYUNIT,*)' SCRATCH VECTOR '
            WRITE(MYUNIT,*)(SCRATCH(IK),IK=1,NATOMS)
            WRITE(MYUNIT,*)' ORBIT DISTANCES: '
            WRITE(MYUNIT,*)(SCRATCH(6*NATOMS+IK),IK=1,NORBIT)
         ENDIF
C
C   Debug print to make sure we've got all the orbits correct.
C
         IF (IPRNT.GE.3) THEN
            WRITE(MYUNIT,470) NORBIT
470         FORMAT(' Number of distinct orbits: ',i3,/,' Summary of orbit sizes:')
            DO I=1,NORBIT
               WRITE(MYUNIT,480) I, NORD2(NATOMS+I)
480            FORMAT(2I3)
            ENDDO
         ENDIF
C
C  The programmer thinks that only the smallest orbit of non-unit
C  magnitude will be needed to identify the point group uniquely
C  by the following method.  Find this orbit and proceed.  
C
C  SCRATCH(6*NATOMS+orbit number) will be the radius of the orbit.
C
         MINORB=NATOMS
         NUMB=0
         IREF=0
         DO 500 I=1,NORBIT
            NUMB=NORD2(NATOMS+I)+NUMB
            IF (SCRATCH(6*NATOMS+I).LT.1.0D-8) GOTO 500
            IF ((NORD2(NATOMS+I).GT.1).AND.(NORD2(NATOMS+I).LE.MINORB)) THEN
               MINORB=NORD2(NATOMS+I)
               ORDIS=SCRATCH(6*NATOMS+I)
               IREF=1+NUMB-NORD2(NATOMS+I)
            ENDIF
500      CONTINUE
         IF (IREF.EQ.0) THEN
            IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'Accidental degeneracy detected'
            TOLE=TOLE/10.0D0
            IF (TOLE.GT.1.0D-7) THEN
               GOTO 651
            ELSE
               IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) ' The full molecular point group is undetermined'
               WRITE(GPSTRING,653) ' The full molecular point group is undetermined'
               GOTO 652
            ENDIF
         ENDIF
         IF (IPRNT.GE.13) WRITE(MYUNIT,510) MINORB,NORD2(IREF)
510      FORMAT(' The smallest non-unit orbit contains ',i3,' members.',/,
     1    ' Atom number ',i3,' belongs to this orbit and will be',
     2    ' used as a reference.')
         IF (IPRNT.GE.13) WRITE(MYUNIT,520)(NORD2(IREF+J),J=1,MINORB-1)  
520      FORMAT(' Other members of this orbit are: ',100I4)
         CALL ZERO(SCRATCH,9*NATOMS)
         IBOT=3*NORD2(IREF)-2
C
C   Reorder NORD2(IREF+J), J=1,MINORB-1 in terms of decreasing distance from
C   atom NORD2(IREF). Coordinates are in NEWQ. Hopefully, this reordering
C   will mean that generators are not missed in the loop below.
C
C   Appears to make no difference at all.
C
CC        PRINT*,'reordering'
CC        PRINT*,'NORD2: ',NORD2(IREF+1:IREF+MINORB-1)
C         DO J=1,MINORB-1
C            DISTANCE(J)=(NEWQ(3*(NORD2(IREF)-1)+1)-NEWQ(3*(NORD2(IREF+J)-1)+1))**2
C     &                 +(NEWQ(3*(NORD2(IREF)-1)+2)-NEWQ(3*(NORD2(IREF+J)-1)+2))**2
C     &                 +(NEWQ(3*(NORD2(IREF)-1)+3)-NEWQ(3*(NORD2(IREF+J)-1)+3))**2
CC           PRINT '(A,2I5,G20.10)','J,NORD2(IREF+J),DISTANCE=',J,NORD2(IREF+J),DISTANCE(J)
C            distloop: DO J1=1,J-1 
C               IF (DISTANCE(J).GT.DISTANCE(J1)) THEN
C                  DUMMY=DISTANCE(J)
C                  IDUMMY=NORD2(IREF+J)
C                  DO J2=J,J1+1,-1
CC                    PRINT*,'changing NORD2 ',IREF+J2,' to NORD2 ',IREF+J2-1
C                     NORD2(IREF+J2)=NORD2(IREF+J2-1)
C                     DISTANCE(J2)=DISTANCE(J2-1)
C                  ENDDO
CC                 PRINT*,'changing NORD2 ',IREF+J1,' to ',IDUMMY
C                  NORD2(IREF+J1)=IDUMMY
C                  DISTANCE(J1)=DUMMY
C                  EXIT distloop
C               ENDIF
C            ENDDO distloop
C         ENDDO
CC        PRINT*,'neworder'
CC        PRINT*,'NORD2: ',NORD2(IREF+1:IREF+MINORB-1)
CC        CALL FLUSH(6)

C
C   Now use reference atom from orbit and loop over all other atoms,
C   this is the second loop (610) below.
C   Outer loop below (620) is fully executed two times
C   if and only if the FPGRP is of the T or I variety and no
C   C2 axis appropriate for orientation fudging was found after
C   identification of the group.
C   First rotate molecule so that the interatomic distance vector
C   is parallel to the z-axis, with bisector along x-axis.
C   Structure now in SCRATCH(1). Here one has to allow for the
C   possibility that atoms I and J are 180 degrees apart; if this
C   is the case, loop through other atoms in orbit and find one
C   that is perpendicular (there has to be one for cubic groups)
C   and use this vector as the "bisector". If
C   it isn't done this way, the orientation fudging gets fudged up.
C
         DO 620 ICOUNT=1,2
            IF (IPRNT.GE.3) WRITE(MYUNIT,*)' PASS ',ICOUNT,' THROUGH ROT. FINDER '
            DO 610 J=1,MINORB-1
               IF (IDONE.EQ.1) GOTO 610
               IF (IPRNT.GE.3) WRITE(MYUNIT,*)' ATOMS ',NORD2(IREF),NORD2(IREF+J)
               IBOT2=3*NORD2(IREF+J)-2
               CALL XVEC(NEWQ(IBOT),NEWQ(IBOT2),SCRATCH(1),0)
               DIST=SQRT(MYDOT(SCRATCH(1),SCRATCH(1),3))
C              CALL VADD(SCRATCH(1),NEWQ(IBOT),NEWQ(IBOT2),3*NATOMS,1) ! bug DJW
               CALL VADD(SCRATCH(1),NEWQ(IBOT),NEWQ(IBOT2),3,1)
C
C   Deal with the linear problem right now.
C
               BILEN=MYDOT(SCRATCH(1),SCRATCH(1),3)
C              IF (BILEN.LT.TOLD) THEN
               IF (BILEN.LT.TOLO) THEN
                  ISET=0
                  DO 530 LOOK=1,MINORB-1
                     IF (ISET.EQ.1) GOTO 530
                     INDEX=3*NORD2(IREF+LOOK)-2
                     TEST=MYDOT(NEWQ(IBOT),NEWQ(INDEX),3)
                     IF (IPRNT.GT.13) WRITE(MYUNIT,*) TEST
C                    IF (ABS(TEST).LT.TOLD) THEN
                     IF (ABS(TEST).LT.TOLO) THEN
                        CALL ZERO(SCRATCH,3)
                        CALL VADD(SCRATCH(1),SCRATCH(1),NEWQ(INDEX),3,1) 
                        ISET=1
                     ENDIF
530               CONTINUE
               ENDIF
               CALL SIAZ(SCRATCH(1),RM,1)
               CALL ZERO(SCRATCH,3*NATOMS*3)
               CALL MATMULV(SCRATCH(3*NATOMS+1),NEWQ,RM,NATOMS,3,3)
C
C   Now you have bisector along x. Rotate about x to bring the two
C   atoms into position parallel to z.
C
               DIP=SQRT(MYDOT(SCRATCH(3*NATOMS+IBOT+1),SCRATCH(3*NATOMS+IBOT+1),2))
               ARGU=SCRATCH(3*NATOMS+IBOT+2)/DIP
               IF (ABS(ARGU).GT.1.D0) ARGU=ARGU/ABS(ARGU)
               ANGLE2=-ACOS(ARGU)/DTOR
               IF (SCRATCH(3*NATOMS+IBOT+1).GT.0.D0)ANGLE2=-ANGLE2
               CALL ROTM(1,ANGLE2,1,RM)
               CALL MATMULV(SCRATCH,SCRATCH(3*NATOMS+1),RM,NATOMS,3,3)
C
C   Now find new vector which bisects the two atoms in question -
C   it had best be x! (Unless linear case where TATB is forced.)
C
               CALL VADD(TATB,SCRATCH(IBOT),SCRATCH(IBOT2),3,1)
               IF (ISET.EQ.1) THEN
                  CALL ZERO(TATB,3)
                  TATB(1)=1.D0
               ENDIF
C
C 1. Test if the two atoms are connected by symmetry axes, looping
C    over possibilities (2,3,4 and 5). T1 contains the matrix transformation
C    required up to this point from NEWQ, and T1T is its transpose (inverse).
C
C    The algorithm uses the fact that if there is an n-fold axis then the
C    edge chosen for the reference atoms is part of an n-sided polygon. ANGMAG
C    is the angle between the vector through the centre of this hypothetical
C    polygon and the vector which bisects the reference edge. SIAZ is then
C    called to determine the rotation needed to map the hypothetical n-fold
C    axis onto the z axis. This transformation is saved in T2, and its
C    transpose in T2T.
C
               DO 600 IAXORD=5,2,-1
                  ICLIP=0
                  ANGIAX=360.D0/(1.0D0*IAXORD)
                  ANGMG=ANGMAG(ORDIS,DIST,IAXORD,IERR)
                  IF (IERR.EQ.1) GOTO 600
570               CALL ROTM(3,ANGMG,1,RM)
                  CALL MATMULV(SCRATCH(3*NATOMS+4),TATB,RM,1,3,3)
                  CALL SIAZ(SCRATCH(3*NATOMS+4),RM,3)
                  CALL ZERO(SCRATCH(3*NATOMS+1),NATOMS*3*2)
                  CALL MATMULV(SCRATCH(3*NATOMS+1),SCRATCH,RM,NATOMS,3,3)  
C                 CALL SORTXYZ(SCRATCH(3*NATOMS+1),QSORT,NORD,TOLD,NATOMS)
                  QSORT(1:3*NATOMS)=SCRATCH(3*NATOMS+1:6*NATOMS)
                  CALL ROTM(3,ANGIAX,1,RM)
                  CALL MATMULV(SCRATCH(6*NATOMS+1),SCRATCH(3*NATOMS+1),RM,NATOMS,3,3)
C                 CALL COMPARE2(SCRATCH(6*NATOMS+1),QSORT,NORD,IIAX,TOLD,NATOMS,IPRNT) 
C                 CALL MINPERM(NATOMS,QSORT,SCRATCH(6*NATOMS+1),0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
                  CALL TESTSYMOP(NATOMS,QSORT,SCRATCH(6*NATOMS+1),PERM,TOLD,DIST2,WORSTRAD)
C                 CALL BIPARTITE(NATOMS,QSORT,SCRATCH(6*NATOMS+1),PERM,DUMMY,DIST2,WORSTRAD)
C
C   Point group determination happens here, as well as orientation
C   fudging.
C
C                 IF (IIAX.EQ.0) THEN
                  IF (IPRNT.GT.4) WRITE(MYUNIT,'(A,I5,3F12.3)') 'IAXORD,DUMMY,DIST2,TOLD=',IAXORD,DUMMY,DIST2,TOLD
C                 PRINT*,'QSORT:'
C                 PRINT '(3G20.10)',QSORT(1:3*NATOMS)
C                 PRINT*,'SCRATCH(6*NATOMS+1):'
C                 PRINT '(3G20.10)',SCRATCH(6*NATOMS+1:9*NATOMS)
                  IF (DIST2.LT.TOLD) THEN
                     IF (IGEN.LT.100) THEN
C
C   Make the transformation matrix corresponding to the reference state in QREF
C
!                       IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient B'
                        CALL REORIENT(NATOMS,SCRATCH(6*NATOMS+1),QREF,RMAT)
                        CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
                        CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
                        DO J2=1,IGEN
                           CALL CMPMAT(GENMAT,100,J2,DUM2,NEW,MATDIFF)
!                          WRITE(MYUNIT,'(A,2I6,G20.10)') 'comparing matrices: J2,IGEN,MATDIFF=',J2,IGEN,MATDIFF
!                          WRITE(MYUNIT,'(A,L5)') 'after cmpmat NEW=',NEW
                           IF (.NOT.NEW) GOTO 20
                        ENDDO
                        IGEN=IGEN+1
                        GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!                       IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here B IGEN=',IGEN,' new GENMAT:'
!                       IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
20                      CONTINUE
                     ENDIF
                     IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,580) IAXORD
580                  FORMAT(' Axis of order ',I1,' identified.')
                     IF (IAXORD.EQ.5) THEN
                        FPGRP='I  '
                        HPTGRP=60
C                       IF (ISINV.EQ.0) FPGRP='I h'
                        IF (DUMMYSINV.LT.TOLD) THEN
                           FPGRP='I h'
                           HPTGRP=120
                        ENDIF
                     ENDIF
                     IF (IAXORD.EQ.4) THEN
                        FPGRP='O  '
                        HPTGRP=24
C                       IF (ISINV.EQ.0) FPGRP='O h'
                        IF (DUMMYSINV.LT.TOLD) THEN
                           FPGRP='O h'
                           HPTGRP=48
                        ENDIF
C
C   Orientation is good. Put structure into NEWQ and leave loop.
C
                        CALL ZERO(NEWQ,3*NATOMS)
                        CALL VADD(NEWQ,NEWQ,SCRATCH(3*NATOMS+1),3*NATOMS,1)
C                       CALL SORTXYZ(NEWQ,QSORT,NORD,TOLD,NATOMS)
                        QSORT(1:3*NATOMS)=NEWQ(1:3*NATOMS)
                        IDONE=1
                        GOTO 610
                     ENDIF
                     IF (IAXORD.EQ.3) THEN
                        IF (FPGRP.EQ.'   ') THEN
                           FPGRP='T  '
                           HPTGRP=12
C                          IF (ISINV.EQ.0) FPGRP='T h'
                           IF (DUMMYSINV.LT.TOLD) THEN
                              FPGRP='T h'
                              HPTGRP=24
                           ENDIF
                        ENDIF
                     ENDIF
                     IF (IAXORD.EQ.2) THEN
                        IF (FPGRP.NE.'   '.AND.FPGRP.NE.'O h'.AND.FPGRP.NE.'O  ') THEN
C
C   Have a C2 axis for I or T group. This is very good. Td and T can be differentiated
C   here.  Do this now if appropriate. Then put structure into NEWQ and
C   leave loop if it is an I group. Otherwise, need to continue through
C   (since there might be a higher order axis down the line.) If
C   pass two, however, then exit for T since it must be the correct group.
C
                           IF (ICOUNT.EQ.2.AND.(FPGRP.EQ.'T  '.OR.FPGRP.EQ.'T h')) THEN
                              DO 590 JPLANE=1,3
                                 CALL REFLECT(SCRATCH(3*NATOMS+1),SCRATCH,NATOMS,JPLANE)
C                                CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMPQ,TOLD,NATOMS,IPRNT)
C                                CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
                                 CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C                                CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
C                                IF (ICOMPQ.EQ.0.AND.ISINV.NE.0) FPGRP='T d'
                                 IF ((DIST2.LT.TOLD).AND.(DUMMYSINV.GT.TOLD)) THEN
                                    FPGRP='T d'
                                    HPTGRP=24
                                 ENDIF
                                 CALL ZERO(NEWQ,3*NATOMS)
                                 CALL VADD(NEWQ,NEWQ,SCRATCH(3*NATOMS+1),3*NATOMS,1)
C                                CALL SORTXYZ(NEWQ,QSORT,NORD,TOLD,NATOMS)
                                 QSORT(1:3*NATOMS)=NEWQ(1:3*NATOMS)
                                 IDONE=1
                                 GOTO 610
590                           CONTINUE
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  ANGMG=-ANGMG
                  IF (ICLIP.EQ.0) THEN
                     ICLIP=1
                     GOTO 570
                  ENDIF
600            CONTINUE
610         CONTINUE
620      CONTINUE
         IF (FPGRP.EQ.'   ') THEN
            IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'Accidental degeneracy detected'
            TOLE=TOLE/10.0D0
            IF (TOLE.GT.1.0D-7) THEN
               GOTO 651
            ELSE
               IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) ' The full molecular point group is undetermined'
               WRITE(GPSTRING,653) ' The full molecular point group is undetermined'
               GOTO 652
            ENDIF
         ENDIF
      ENDIF
C
C   *** This is the end of the IDEGEN=2 block. ***
C
C   Check symmetry operations belonging to abelian groups
C
630   ANG=180.D0
      IREF=0
      DO I=3,1,-1
         CALL REFLECT(NEWQ,SCRATCH,NATOMS,I)
C        CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C        IF (ICOMP.EQ.0) THEN
C        CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
         CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C        CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
         IF (DIST2.LT.TOLD) THEN
            IF (IGEN.LT.100) THEN
               DO J1=1,3
                  DO J2=1,3
                     RM(J2,J1)=0.0D0
                  ENDDO
                  RM(J1,J1)=1.0D0
                  IF (J1.EQ.I) RM(J1,J1)=-1.0D0
               ENDDO
C
C   Make the transformation matrix corresponding to the reference state in QREF 
C
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient C'
               CALL REORIENT(NATOMS,SCRATCH,QREF,RMAT) 
               CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
               CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
               DO J2=1,IGEN
                  CALL CMPMAT(GENMAT,100,J2,DUM2,NEW,MATDIFF)
                  IF (.NOT.NEW) GOTO 30
               ENDDO
               IGEN=IGEN+1
               GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here C IGEN=',IGEN,' new GENMAT:'
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
30             CONTINUE
            ENDIF
            IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,640) I
640         FORMAT(' Reflection in plane ',i2,' is a valid symmetry operation.')
            IREF=IBTOR(IREF,2**I/2)
         ENDIF
         CALL ROTM(I,ANG,1,RM)
         CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
C        CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C        IF (ICOMP.EQ.0) THEN
C        CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
         CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C        CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
         IF (DIST2.LT.TOLD) THEN
            IF (IGEN.LT.100) THEN
C
C   Make the transformation matrix corresponding to the reference state in QREF 
C
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient D'
               CALL REORIENT(NATOMS,SCRATCH,QREF,RMAT) 
               CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
               CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
               DO J2=1,IGEN
                  CALL CMPMAT(GENMAT,100,J2,DUM2,NEW,MATDIFF)
                  IF (.NOT.NEW) GOTO 40
               ENDDO
               IGEN=IGEN+1
               GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here D IGEN=',IGEN,' new GENMAT:'
!              IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
40             CONTINUE
            ENDIF
            IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,650) I
650         FORMAT(' C2 rotation about ',i2,' is a valid symmetry operation ')
            IROT=IBTOR(IROT,2**I/2)
         ENDIF
      ENDDO
      DO I=1,NATOMS*3
C        SCRATCH(I)=-NEWQ(I)
         SCRATCH(I)=-QSORT(I)
      ENDDO
C     CALL COMPARE2(SCRATCH,QSORT,NORD,ICOMP,TOLD,NATOMS,IPRNT)
C     IF (ICOMP.EQ.0) THEN
C     CALL MINPERM(NATOMS,QSORT,SCRATCH,0.0D0,0.0D0,0.0D0,.FALSE.,PERM,DUMMY,DIST2,WORSTRAD)
      CALL TESTSYMOP(NATOMS,QSORT,SCRATCH,PERM,TOLD,DIST2,WORSTRAD)
C     CALL BIPARTITE(NATOMS,QSORT,SCRATCH,PERM,DUMMY,DIST2,WORSTRAD)
C
C  Even DNd do not have the inversion operation - bug fixed DJW 13/5/92
C
C     PRINT '(A,I6,G20.10)','ISAXIS,DIST2=',ISAXIS,DIST2
C     IF ((ISAXIS.EQ.1).OR.(DIST2.LT.TOLD)) THEN
      IF (DIST2.LT.TOLD) THEN
         IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,680)
680      FORMAT(' The molecule possesses an inversion center. ')
         IF (IGEN.LT.100) THEN
            DO J1=1,3
               DO J2=1,3
                  RM(J2,J1)=0.0D0
               ENDDO
               RM(J1,J1)=-1.0D0
            ENDDO
C
C   Make the transformation matrix corresponding to the reference state in QREF 
C
!           IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'calling reorient E'
            CALL REORIENT(NATOMS,SCRATCH,QREF,RMAT) 
            CALL MYMATMUL(DUM,RM,RMAT,3,3,3)
            CALL MATMULV(DUM2,DUM,RMAT,3,3,3)
            DO J2=1,IGEN
               CALL CMPMAT(GENMAT,100,J2,DUM2,NEW,MATDIFF)
               IF (.NOT.NEW) GOTO 60
            ENDDO
            IGEN=IGEN+1
            GENMAT(IGEN,1:3,1:3)=DUM2(1:3,1:3)
!           IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,*) 'here E IGEN=',IGEN,' new GENMAT:'
!           IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,'(9F12.4)') GENMAT(IGEN,1:3,1:3)
60          CONTINUE
         ENDIF
         IF (ILINEAR.EQ.1) THEN
            FPGRP='DXh'
            HPTGRP=2
         ENDIF
         IF (FPGRP.EQ.'   '.AND.IDGRP.EQ.1) THEN
            FPGRP='DNd'
            HPTGRP=4*IHIGH
         ENDIF
         IINV=1
      ELSE
         IF (ILINEAR.EQ.1) THEN
            FPGRP='CXv'
            HPTGRP=1
         ENDIF
         IF (FPGRP.EQ.'   '.AND.IDGRP.EQ.1) THEN
            FPGRP='DN '
            HPTGRP=2*IHIGH
         ENDIF
      ENDIF
C
C   Determine largest abelian subgroup of molecule
C
      IF (IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP='C1 '
      IF (IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.1) PGRP='C i'
      IF (IROT.EQ.0.AND.IREF.NE.0.AND.IINV.EQ.0) PGRP='C s'
      IF (IROT.NE.0.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP='C2 '
      IF (IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.0) PGRP='C2v'
      IF (IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.1) PGRP='C2h'
      IF (IROT.EQ.7.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP='D2 '
      IF (IROT.EQ.7.AND.IREF.EQ.7.AND.IINV.EQ.1) PGRP='D2h'
      IF (IPRNT.GE.4) WRITE(MYUNIT,690) IREF,IROT,IINV
690   FORMAT ('Symmetry bits: ',3(1X,I3))
C
C   Put NEWQ into Q.
C
      BPGRP=PGRP
      CALL ZERO(SCRATCH,NATOMS*3)
      CALL VADD(Q,NEWQ,SCRATCH(1),NATOMS*3,1)

      IF (PGRP.EQ.'C2h') PGRP='C s'

      IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,770)
770   FORMAT(80('*'))
C
C  This should cope with anything we are likely to come up against,
C  until Mark tries a ring with 100 atoms!
C
      IF (IHIGH.EQ.0) THEN 
         WRITE(JNKSTR,'(A3)') FPGRP(1:3)
      ELSE IF (IHIGH.LT.10) THEN
         WRITE(JNKSTR,'(A1,I1,A1)') FPGRP(1:1),IHIGH,FPGRP(3:3)
      ELSE
         WRITE(JNKSTR,'(A1,I2,A1)') FPGRP(1:1),IHIGH,FPGRP(3:3)
      ENDIF
      FPGRP=JNKSTR
      IF (IDEGEN.GT.0) THEN
         IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,750) FPGRP
750      FORMAT(' The full molecular point group is ',A4,'.')
         WRITE(GPSTRING,750) FPGRP
      ENDIF
      IF (IDEGEN.EQ.0) THEN
         IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,750) STRING(BPGRP,IHIGH)
         FPGRP= STRING(BPGRP,IHIGH)
         WRITE(GPSTRING,750) STRING(BPGRP,IHIGH)
      ENDIF
      IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,790) BPGRP
790   FORMAT(' The largest Abelian subgroup of the full molecular point group is ',A4,'.')
      IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,760) TOLD, TOLE, TOLO
760   FORMAT(' Distance tolerance=',F12.5,' Inertia tolerance=',F12.5,' Orbit tolerance=',F12.5)  
      IF (DEBUG.OR.(IPRNT.GE.3)) WRITE(MYUNIT,770)

652   CONTINUE

      RETURN
      END

      SUBROUTINE SIAZ(VEC,ROT,IAXIS)
      IMPLICIT NONE
      INTEGER IAXIS,J1,J2
      DOUBLE PRECISION VEC(3),ROT(3,3),RM(3,3),XTEMP(3,3),ANG1,ARG2,DISP,AZIM,ARG1,DIS,MYDOT
C
C     Write (6,*) 'SIAZ: vec, iaxis'
C     Write (6,'(a,3f12.6,i3)') 'SIAZ: ',vec, iaxis
C 
C DETERMINE LENGTH OF VECTOR.
C
      ROT(1,1)=1.0D0
      ROT(2,2)=1.0D0
      ROT(3,3)=1.0D0
      ROT(1,2)=0.0D0
      ROT(1,3)=0.0D0
      ROT(2,1)=0.0D0
      ROT(2,3)=0.0D0
      ROT(3,1)=0.0D0
      ROT(3,2)=0.0D0
      DIS=SQRT(MYDOT(VEC,VEC,3))
      IF (ABS(DIS).LT.1.0D-10) RETURN
C     Write (6,'(a,f12.6)') 'SIAZ: dis', dis
C
C   RETURNS ROTATION MATRIX NEEDED TO PLACE VECTOR
C   VEC ALONG THE X (IAXIS=1), Y (IAXIS=2) OR Z (IAXIS=3)
C   AXIS.  TAIL OF VECTOR IMPLICITLY ASSUMED TO BE AT 
C   THE ORIGIN.
C
      IF(IAXIS.EQ.3)THEN
        ARG1=VEC(3)/DIS
        IF(ABS(ARG1).GT.1.D0)ARG1=ARG1/ABS(ARG1)
        AZIM=-1.D0*ACOS(ARG1)
        DISP=SQRT(MYDOT(VEC,VEC,2))
        IF (ABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(1)/DISP
        IF(ABS(ARG2).GT.1.D0)ARG2=ARG2/ABS(ARG2)
        ANG1=ACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(2))
        CALL ROTM(3,ANG1,0,RM)
        CALL ROTM(2,AZIM,0,ROT)
        CALL MATMUL(XTEMP,RM,ROT,3,3,3,3,3,3)
        DO 17 J1=1,3
           DO 18 J2=1,3
              ROT(J2,J1)=XTEMP(J2,J1)
18         CONTINUE
17      CONTINUE
      ELSEIF(IAXIS.EQ.1)THEN
        ARG1=VEC(1)/DIS
        IF(ABS(ARG1).GT.1.D0)ARG1=ARG1/ABS(ARG1)
        AZIM=-1.D0*ACOS(ARG1)
        DISP=SQRT(MYDOT(VEC(2),VEC(2),2))
        IF (ABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(2)/DISP
        IF(ABS(ARG2).GT.1.D0)ARG2=ARG2/ABS(ARG2)
        ANG1=ACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(3))
        CALL ROTM(1,ANG1,0,RM)
        CALL ROTM(3,AZIM,0,ROT)
        CALL MATMUL(XTEMP,RM,ROT,3,3,3,3,3,3)
        DO 27 J1=1,3
           DO 28 J2=1,3
              ROT(J2,J1)=XTEMP(J2,J1)
28         CONTINUE
27      CONTINUE
      ELSEIF(IAXIS.EQ.2)THEN
        ARG1=VEC(2)/DIS
        IF(ABS(ARG1).GT.1.D0)ARG1=ARG1/ABS(ARG1)
        AZIM=-1.D0*ACOS(ARG1)
        DISP=SQRT(DIS**2-VEC(2)**2)
        IF (ABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(3)/DISP
        IF(ABS(ARG2).GT.1.D0)ARG2=ARG2/ABS(ARG2)
        ANG1=ACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(1))
        CALL ROTM(2,ANG1,0,RM)
        CALL ROTM(1,AZIM,0,ROT)
        CALL MATMUL(XTEMP,RM,ROT,3,3,3,3,3,3)
        DO 37 J1=1,3
           DO 38 J2=1,3
              ROT(J2,J1)=XTEMP(J2,J1)
38         CONTINUE
37      CONTINUE
       ENDIF
       RETURN
       END

      SUBROUTINE REFLECT(A,SCRATCH,NATOM,IPLANE)
C
C REFLECTS POINTS IN PLANE
C
      IMPLICIT NONE
      INTEGER NATOM,IPLANE,IX,I,J
      DOUBLE PRECISION A(3,NATOM),SCRATCH(3,NATOM)
 
      IX(I,J)=MIN(I,J)/MAX(I,J)
 
      DO I=1,NATOM
         DO J=1,3
            SCRATCH(J,I)=A(J,I)*(1.0D0-2.D0*IX(IPLANE,J))
         ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE ZERO(A,NA)
      IMPLICIT NONE
      INTEGER NA,I
      DOUBLE PRECISION A(NA)

      DO I=1,NA
        A(I)=0.D0
      END DO
      CONTINUE
      RETURN
      END

      SUBROUTINE VADD(A,B,C,N,IP)
      IMPLICIT NONE
      INTEGER N,IP,I
      DOUBLE PRECISION A(N),B(N),C(N)

      DO I=1,N
        A(I)=B(I)+C(I)*IP
      END DO
      RETURN
      END         
      SUBROUTINE COMPARE2(VECS,VECR,NORD,ICOMP,TOL,NATOMS,IPRNT)
C
C     ROBUST EQUIVALENCE CHECK - DO WELL DEFINED SORT ON COORDINATE
C     MATRIX AND COMPARE ELEMENT BY ELEMENT.  SHOULD BE FOOLPROOF.
C
C     VECS     coordinate vector to be checked (modified)
C     VECR     sorted reference coordinate vector (input only)
C     NORD     ???
C     ICOMP    number of coordinates outside of TOL (output only)
C     TOL      tolerance for comparison of coords (input only)
C
      IMPLICIT NONE
      INTEGER I,JAP,NATOMS,ICOMP,IPRNT
      DOUBLE PRECISION VECS(3*NATOMS),VECR(3*NATOMS),Z,TOL
      INTEGER NORD(2*NATOMS+1)
 
      ICOMP=0
      CALL SORTXYZ(VECS,VECS,NORD(NATOMS+1),TOL,NATOMS)
      IF (IPRNT.GT.10) THEN
         WRITE(6,*)'Sorted reference vector (input)'
         WRITE(6,'(3G20.10)')(VECR(JAP),JAP=1,3*NATOMS)
         WRITE(6,*)'Sorted transformed vector'
         WRITE(6,'(3G20.10)')(VECS(JAP),JAP=1,3*NATOMS)
      ENDIF
      DO I=1,NATOMS*3
         Z=ABS(VECR(I)-VECS(I))
         IF (Z.GT.TOL) THEN
C           PRINT '(A,I5,4G20.10)','I,VECR,VECS,diff,TOL=',I,VECR(I),VECS(I),Z,TOL
            ICOMP=ICOMP+1
            RETURN
         ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE SORTXYZ(XX,Y,NORD,TOL,NATOMS)
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR 
C
      IMPLICIT NONE
      INTEGER NATOMS,J,JK,I,NORD(NATOMS*2+1),J1
      DOUBLE PRECISION X(3*NATOMS),XX(3*NATOMS),Y(3*NATOMS),TOL
C
C     Sort on the X - if two X's are equivalent, sort on Y and so on.
C     If the coordinate < the tolerance we should ignore it! However,
C     if the tolerance is sloppy that can lead to the sorting ignoring
C     genuine small differences between coordinates. Sigh. 
C
C     PRINT '(A)','in SORTXYZ'
      DO 10 I=1,3*NATOMS
         X(I)=XX(I)
10    CONTINUE
      JK=1
40    J=1
C     PRINT*,'TOL=',TOL
      outer: DO I=1,3*NATOMS-2,3
C        DO J1=1,NATOMS
C           PRINT '(A,I5,3G20.10)','J1,coords: ',J1,X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3)
C        ENDDO
C        PRINT '(A,2I5,3G20.10)','I,J,X(I),X(J),diff    =',I,J,X(I),X(J),X(I)-X(J)
         IF (X(I)-X(J).GT.TOL) THEN
            J=I          ! J is set to the index of the x coordinate of the
            CYCLE outer  ! the corresponding coordinate is set to -99999.D0, so it won;t be used again.
         ENDIF
         IF (ABS(X(I)-X(J)).LT.TOL) THEN
C           PRINT '(A,2I5,3G20.10)','I,J,X(I+1),X(J+1),diff=',I,J,X(I+1),X(J+1),X(I+1)-X(J+1)
            IF (X(I+1)-X(J+1).GT.TOL) J=I
            IF (ABS(X(I+1)-X(J+1)).LT.TOL) THEN
               IF (X(I+2)-X(J+2).GT.TOL) J=I
            ENDIF
         ENDIF
      ENDDO outer
      DO I=0,2
         Y(3*JK-2+I)=X(J+I)
         X(J)=-99999.D0
      END DO
      NORD(JK)=(J+2)/3
      JK=JK+1
      IF (JK.EQ.NATOMS+1) GOTO 70
      GOTO 40
70    CONTINUE
      RETURN
      END

      SUBROUTINE ROTM(IDIR,ANG,IX,RT)
C
C     FORMS THREE DIMENSIONAL ROTATION MATRIX (IX=0 IF ANG IN RADS,
c     IX=1 IF DEG,  ADD 10 TO IX TO GET TRANSPOSE)
C
      IMPLICIT NONE
      INTEGER IDIR,IX,I,J,IBTAND,IEO
      DOUBLE PRECISION RT(3,3),ANG,ATOR,TANG,INK

      IBTAND(I,J)=IAND(I,J)
      IEO(I)=2*IBTAND(I,1)-1
C
      ATOR=ACOS(-1.D0)/180.D0
      INK=IX/10
      CALL ZERO(RT,9)
      TANG=ANG
      IF(IX.EQ.1)TANG=ANG*ATOR
      DO J=1,3
         RT(J,J)=DCOS(TANG)
         RT(IDIR,J)=IBTAND(IDIR,J)/MAX(IDIR,J)
      END DO
      DO I=1,3
         DO J=1,I
            IF(I.NE.J.AND.I.NE.IDIR.AND.J.NE.IDIR)THEN
               RT(I,J)=IEO(I*J)*DSIN(TANG)
               RT(J,I)=-RT(I,J)
            ENDIF
         END DO
      END DO
      IF(INK.NE.0)CALL MTRANS(RT,RT,3,3,3,3)
      RETURN
      END

      SUBROUTINE INERTIA(Q,IT,NATOMS)
      USE COMMONS,ONLY : ATMASS
C
C BUILD INERTIA TENSOR
C
      IMPLICIT NONE
      INTEGER I,J,IND,K,NATOMS
      DOUBLE PRECISION IT(3,3),SCRATCH(9),Q(3*NATOMS),FDIST
 
      IND(I,J)=3*(I-1)+J
 
      CALL ZERO(SCRATCH,9)
      CALL ZERO(IT,9)
      DO I=1,3
         DO K=1,NATOMS
            IT(I,I)=ATMASS(K)*(FDIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)+IT(I,I)
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,3
            IF (I.NE.J) THEN
               DO K=1,NATOMS
                  IT(I,J)=-ATMASS(K)*(Q(IND(K,I))*Q(IND(K,J)))+IT(I,J)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
    
      SUBROUTINE XVEC(A,B,V,IX)
C
C CALCULATES THE VECTOR V BETWEEN CARTESIAN POINTS A AND B
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IX, I
      DIMENSION A(3),B(3),V(3)
      DO I=1,3
        V(I)=B(I)-A(I)
      END DO
      IF(IX.EQ.1)CALL NORMAL(V,3)
      RETURN
      END 

      SUBROUTINE MYCROSS(A,B,C,IX)
C
C CALCULATES THE (OPTIONALLY) NORMALIZED VECTOR CROSS PRODUCT C=A x B
C
      IMPLICIT NONE
      INTEGER IX
      DOUBLE PRECISION A(3),B(3),C(3)

      C(3)=A(1)*B(2)-B(1)*A(2)
      C(2)=-A(1)*B(3)+A(3)*B(1)
      C(1)=A(2)*B(3)-A(3)*B(2)
      IF(IX.EQ.1)CALL NORMAL(C,3)
      RETURN
      END
C
C********************************************************************
C
      SUBROUTINE EIG(A,B,L,N,N1)
C
C DIAGONALIZATION BY THE JACOBI METHOD.
C A - MATRIX TO BE DIAGONALIZED (eigenvalues returned in diagonal
C        elements of A).  If you want to save A, you must do this before
C        calling EIG.  Set N to the same value as L.
C B - EIGENVECTORS
C L - DIMENSION OF A AND B
C N - SIZE OF SUBMATRIX USED
C N1 - A FLAG INDICATING WHETHER THE EIGENVECTORS AND
C      EIGENVALUES ARE TO BE REORDERED.
C
CSW      1
      IMPLICIT NONE
      INTEGER L,N,N1,JJ,IOFF,I,IM1,MU,MM,JI,II,J
      DOUBLE PRECISION A(L,L),B(L,L),W2,W1,C,T,ALP,ALN,D,SUM,S,DIFF,R,Q,P,TOL2,TOL,ZER,ONE
      DATA ZER/0.D00/,ONE/1.D00/
CSW     2
      TOL=1.D-14
      TOL2=1.D-10
      JJ=0
      IOFF=0
      B(1,1)=ONE
      IF(N.EQ.1) RETURN
      DO I=2,N
        IM1=I-1
        DO J=1,IM1
          B(I,J)=ZER
          B(J,I)=ZER
        END DO
        B(I,I)=ONE
      END DO
C
C FIRST SEE IF MATRIX IS ALREADY DIAGONAL- IF SO THEN
C  TAKE APPROPRIATE ACTION
C
      DO II=1,L
        DO JI=II+1,L
          IF(ABS(A(II,JI)).GT.TOL2)IOFF=IOFF+1
          IF(ABS(A(JI,II)).GT.TOL2)IOFF=IOFF+1
        END DO
      END DO
      IF(IOFF.EQ.0)THEN
          CALL ZERO(B,L*L)
          DO 40 I=1,L
          B(I,I)=ONE
40        CONTINUE
      ELSE
50    P=ZER
      DO 70 I=2,N
      IM1=I-1
      DO 60 J=1,IM1
      Q=A(I,J)
      IF(P.GE. ABS(Q)) GO TO 60
      P= ABS(Q)
      II=I
      JJ=J
60    CONTINUE
70    CONTINUE
      IF(P.EQ.0.) GO TO 140
      P=A(II,II)
      Q=A(II,JJ)
      R=A(JJ,JJ)
      DIFF=0.5D0*(P-R)
      IF( ABS(DIFF).LT. ABS(Q)) GO TO 80
      IF( ABS(Q/DIFF).GT.TOL) GO TO 80
      A(II,JJ)=ZER
      A(JJ,II)=ZER
      GO TO 50
80    S=SQRT(0.250D0*(P-R)**2+Q**2)
      SUM=0.5D0*(P+R)
      D=R*P-Q**2
      IF(SUM.GT.ZER) GO TO 90
      ALN=SUM-S
      ALP=D/ALN
      GO TO 100
90    ALP=SUM+S
      ALN=D/ALP
100   IF(DIFF.GT.ZER) GO TO 110
      T=Q/(DIFF-S)
      A(II,II)=ALN
      A(JJ,JJ)=ALP
      GO TO 120
110   T=Q/(DIFF+S)
      A(II,II)=ALP
      A(JJ,JJ)=ALN
120   C=1.0D0/SQRT(1.0D0+T**2)
      S=T*C
      A(II,JJ)=ZER
      A(JJ,II)=ZER
      DO 130 I=1,N
      P=B(I,II)
      Q=B(I,JJ)
      B(I,II)=C*P+S*Q
      B(I,JJ)=C*Q-S*P
      IF(I.EQ.II.OR.I.EQ.JJ) GO TO 130
      P=A(I,II)
      Q=A(I,JJ)
      R=C*P+S*Q
       A(I,II)=R
      A(II,I)=R
      R=Q*C-P*S
      A(I,JJ)=R
      A(JJ,I)=R
130   CONTINUE
      GO TO 50
      ENDIF
140   IF(N1.EQ.1) RETURN
      MM=N-1
      DO I=1,MM
        II=I+1
        DO J=II,N
          IF(A(I,I)-A(J,J) .LT. 0) THEN
            CYCLE
          ELSE 
            GOTO 150
          END IF
150       W1=A(I,I)
          A(I,I)=A(J,J)
          A(J,J)=W1
          DO MU=1,N
            W2=B(MU,I)
            B(MU,I)=B(MU,J)
            B(MU,J)=W2
          END DO
        END DO
      END DO
      RETURN
      END

C
C  Returns A = B C
C
      SUBROUTINE MATMUL(A,B,C,NA,NB,NC,NTA,NTB,NTC)
      IMPLICIT NONE
      INTEGER NA,NB,NC,NTA,NTB,NTC,I,J1,J,K
      DOUBLE PRECISION B(NA,NB),C(NB,NC),A(NA,NC),T1(NTB),Z

      DO 10 I=1,NTA
         DO 15 J1=1,NTB
            T1(J1)=B(I,J1)
15       CONTINUE
         DO 30 J=1,NTC
            Z=0.0D0
            DO 20 K=1,NTB
               Z=Z+T1(K)*C(K,J)
20          CONTINUE
            A(I,J)=Z
30       CONTINUE
10    CONTINUE
      RETURN
      END     
      
C
C               t
C  Returns A = C B
C
      SUBROUTINE MATMULV(A,B,C,NA,NB,NC)
      IMPLICIT NONE
      INTEGER NA, NB, NC, I, K, J
      DOUBLE PRECISION A(NC,NA),B(NB,NA),C(NB,NC),Z

      DO I=1,NA
         DO K=1,NC
            Z=0.D0
            DO J=1,NB
               Z=Z+B(J,I)*C(J,K)
            ENDDO
            A(K,I)=Z
         ENDDO
      ENDDO
      RETURN 
      END
C                
C  Returns A = B C
C
      SUBROUTINE MYMATMUL(A,B,C,NA,NB,NC)
      IMPLICIT NONE
      INTEGER NA, NB, NC, I, K, J
      DOUBLE PRECISION A(NC,NA),B(NB,NA),C(NB,NC),Z

      DO I=1,NA
         DO K=1,NC
            Z=0.D0
            DO J=1,NB
               Z=Z+B(K,J)*C(J,I)
            ENDDO
            A(K,I)=Z
         ENDDO
      ENDDO
      RETURN 
      END

      DOUBLE PRECISION FUNCTION MYDOT(A,B,N)
      IMPLICIT NONE
      INTEGER N, I
      DOUBLE PRECISION A(N),B(N)

      MYDOT=0.D0
      DO I=1,N
        MYDOT=MYDOT+A(I)*B(I)
      END DO
      RETURN
      END

        SUBROUTINE PIKSR2(N,ARR,NLIST)
C
C SIMPLE SORTER MODIFIED TO MAKE JFS HAPPY.
C
        IMPLICIT NONE
        INTEGER N, J, I
        DOUBLE PRECISION ARR(N), A
        INTEGER NLIST(N), NLST

        DO 12 J=2,N
        A=ARR(J)
        NLST=NLIST(J)
        DO 11 I=J-1,1,-1
         IF(ARR(I).LE.A)GOTO 10
         ARR(I+1)=ARR(I)
         NLIST(I+1)=NLIST(I)
   11   CONTINUE
        I=0
   10   ARR(I+1)=A
        NLIST(I+1)=NLST
   12   CONTINUE
        RETURN
        END

        DOUBLE PRECISION FUNCTION ANGMAG(ORBIT,D,IORDER,IRETURN)
        IMPLICIT NONE
        INTEGER IRETURN, IORDER
        DOUBLE PRECISION DTOR, ORDER, A1, X, TOP, BOT, ORBIT, D

        IRETURN=0
        ANGMAG=0.D0
        DTOR=ACOS(-1.D0)/180.D0
        ORDER=1.0D0*(IORDER)
        A1=180.D0/ORDER
        IF(ABS(A1-90.D0).LT.1.D-3)THEN
          ANGMAG=0.D0
           RETURN
        ENDIF
        X=ORBIT**2-0.25D0*D**2/DSIN(DTOR*A1)**2
        TOP=D/DTAN(DTOR*A1)
        IF (X.LT.-1.D-12) THEN
C
C THE ROTATION IS IMPOSSIBLE
C
           IRETURN=1
           RETURN
        ENDIF
        IF(ABS(X).LT.1.D-6)THEN
          ANGMAG=90.D0
        ELSE
           BOT=2.D0*SQRT(X)
          ANGMAG=DATAN(TOP/BOT)/DTOR
        ENDIF
        RETURN
        END         

      CHARACTER(LEN=4) FUNCTION STRING(NAME,Iord)
      CHARACTER(LEN=4) NAME, ItoA*3
      INTEGER IORD
      If (Name(2:2) .eq. 'N') then
         String(1:1)=Name(1:1)
         String(2: )=ItoA(Iord, 0)
         String(LNBlnk(String)+1:)=Name(3:3)
      Else
         String=Name
      EndIf
      RETURN
      END

      SUBROUTINE MTRANS(A,AT,NR,NC,NTR,NTC)
      IMPLICIT NONE
      INTEGER NR,NC,NTR,NTC, I, J
      DOUBLE PRECISION A(NR,NC),AT(NC,NR)
      DO I=1,NTR
        DO J=1,NTC
          AT(J,I)=A(I,J)
        END DO
      END DO
      RETURN
      END
      Double Precision FUNCTION FDIST(A,B)
C
C CALCULATES THE DISTANCE BETWEEN TWO POINTS IN CARTESIAN SPACE
C
      IMPLICIT NONE
      DOUBLE PRECISION A(3),B(3),Z
      INTEGER I

      Z=0.D0
      DO I=1,3
        Z=Z+(A(I)-B(I))**2
      END DO
      FDIST=SQRT(Z)
      RETURN
      END

      SUBROUTINE NORMAL(X,N)
C
C NORMALIZES VECTOR X
C
      DOUBLE PRECISION X(N), Q, P
      INTEGER N, I
      Q=0.D0
      DO I=1,N
        Q=Q+X(I)**2
      END DO
      P=SQRT(Q)
      IF(P.LT.1D-14)THEN
      WRITE(6,*)' null vector returned from NORMAL'
      RETURN
      ENDIF
      DO I=1,N
        X(I)=X(I)/P
      END DO 
      RETURN
      END 
       CHARACTER*(*) FUNCTION ITOA(NR, FRCPLS)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Convert NR to a left justified string
C
C Arguments:
C     NR       number to be converted (input only)
C     FRCPLS   Force leading '+' if NR positive (input only)
C
C Limitations:
C     May return with incomplete conversion if length of ITOA is too
C     short.  Puts '*' in last position of ITOA to indicate overlow.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Header: itoa.f,v 4.0 89/03/14 01:15:45 bernhold Exp $
C
C $Log:	itoa.f,v $
C Revision 4.0  89/03/14  01:15:45  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:10:22  bernhold
C First working release for VAX
C 
C Revision 2.1  89/01/02  20:36:12  bernhold
C To keep consistent with .u file just checked in.
C 
C     Revision 1.1  88/12/07  13:38:51  bernhold
C     Initial revision
C     
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      INTEGER NR, FRCPLS, NRABS, NDIG, N
      INTEGER I, J
C
C     Clear out the string
C
      DO 10 I=1, LEN(ITOA)
         ITOA=' '
 10   CONTINUE
C
C     Start counting position in string
C
      J=1
      NRABS=ABS (NR)
C
C     Put in sign as appropriate
C
      IF (NR .LT. 0) THEN
         ITOA(J:J)='-'
         J=J + 1
      ENDIF
      IF (FRCPLS .NE. 0 .AND. NR .GT. 0) THEN
         ITOA(J:J)='+'
         J=J + 1
      ENDIF
C
C     Check if we are about to overflow the string
C
      IF (J .GT. LEN(ITOA)) THEN
         ITOA(J-1:J-1)='*'
         RETURN
      ENDIF
C
C     Loop over nr of digits in number
C
      NDIG=INT( LOG10( FLOAT(NRABS)) ) + 1
      DO 100 I=NDIG, 1, -1
         N=MOD ( ( NRABS / (10**(I-1) ) ), 10)
         ITOA(J:J)=CHAR(N + 48)
         J=J + 1
C
C        Check for overflow of the string, but if this is last digit
C        then its okay.
C
         IF (J .GT. LEN(ITOA) .AND. I .GT. 1) THEN
            ITOA(J-1:J-1)='*'
            RETURN
         ENDIF
 100  CONTINUE
C
      RETURN
      END
       INTEGER FUNCTION LNBLNK (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:   Returns the position of the last non-blank character
C
C Arguments: STRING   character string (input only)
C
C Remarks:   All FORTRAN 77 character variables are blank padded on the
C            right.  The intrinsic function LEN returns the dimension
C            of the character object, not the length of the contents.
C            
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Header: lnblnk.f,v 1.1 88/01/11 22:08:15 bernhold Exp $
C
C    Revision 0.0  87/07/24  bernholdt (VAX)
C $Log:	lnblnk.f,v $
C    Revision 1.1  88/01/11  22:08:15  bernhold
C    Initial revision
C    
C
C System:     Standard FORTRAN 77
C
C Copyright 1987 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       INTEGER I
       CHARACTER*(*) STRING
       CHARACTER(LEN=1) BLANK
       PARAMETER (BLANK = ' ')
C
C      Start at the end and work backwards
C
       DO 100 I = LEN(STRING), 1, -1
C         Look for first non-whitespace character
          IF (STRING(I:I) .NE. BLANK) THEN
             LNBLNK = I
             RETURN
          ENDIF
  100  CONTINUE
C
C      If we get this far, the string is empty
       LNBLNK = 0
       RETURN
       END

      SUBROUTINE CMPMAT(XMAT1,D1,N1,XMAT2,TEST,MATDIFF)
!      USE COMMONS,ONLY : MYUNIT
      IMPLICIT NONE
      LOGICAL TEST
      INTEGER D1, N1, J1, J2
      DOUBLE PRECISION XMAT1(D1,3,3), XMAT2(3,3), MATDIFF

      TEST=.FALSE.
      DO J2=1,3
         DO J1=1,3
!           WRITE(MYUNIT,'(A,2I6,3G20.10)') 'J2,J1,xmat1,xmat2,diff=',J2,J1,XMAT1(N1,J2,J1),XMAT2(J2,J1),
!    &                                       ABS(XMAT1(N1,J2,J1)-XMAT2(J2,J1))
            IF (ABS(XMAT1(N1,J2,J1)-XMAT2(J2,J1)).GT.MATDIFF) THEN
               TEST=.TRUE.
               RETURN
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
C
C  Test for whether we have a symmetry operation or not by checking
C  that each atom in vector Q1 has a partner within distance TOLD in
C  vector Q2.
C       Q1(i) <--> Q2(perm(i))
C
      SUBROUTINE TESTSYMOP(NATOMS,Q1,Q2,PERM,TOLD,DIST2,WORSTRAD)
!      USE COMMONS, ONLY : MYUNIT,DEBUG,ATMASS
      USE COMMONS, ONLY : ATMASS
      IMPLICIT NONE
      INTEGER NATOMS
      DOUBLE PRECISION Q1(3*NATOMS), Q2(3*NATOMS), TOLD, DIST2, DIST, TOLDSQ, WORSTRAD, DMIN
      LOGICAL PARTNER, ASSIGNED(NATOMS)
      INTEGER J1, J2, PERM(NATOMS), JWORST
   
      TOLDSQ=TOLD*TOLD
      PERM(1:NATOMS)=0
      ASSIGNED(1:NATOMS)=.FALSE.
      DIST2=-1.0D0
      DO J1=1,NATOMS
         PARTNER=.FALSE.
         DMIN=1.0D100
         innerloop: DO J2=1,NATOMS
            IF (ASSIGNED(J2)) CYCLE innerloop
!
! Don't allow atoms with different masses to be permutable.
!
            IF (ATMASS(J2).NE.ATMASS(J1)) CYCLE innerloop
            DIST=(Q1(3*(J1-1)+1)-Q2(3*(J2-1)+1))**2+ 
     &           (Q1(3*(J1-1)+2)-Q2(3*(J2-1)+2))**2+ 
     &           (Q1(3*(J1-1)+3)-Q2(3*(J2-1)+3))**2
C           IF (DEBUG) WRITE(MYUNIT,'(A,2I6,G20.10)') 'testsymop> J1,J2,DIST=',J1,J2,DIST
            IF (DIST.LT.TOLDSQ) THEN
               PARTNER=.TRUE.
               PERM(J1)=J2
               ASSIGNED(J2)=.TRUE.
               IF (DIST.GT.DIST2) THEN
                  DIST2=DIST
                  JWORST=J1
               ENDIF
C              WRITE(MYUNIT,'(A,2I6,G20.10)') 'testsymop> match J1,J2,DIST=',J1,J2,DIST
               EXIT innerloop
            ENDIF
         ENDDO innerloop
         IF (.NOT.PARTNER) THEN
!           IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'testsymop> no partner for atom ',J1
            DIST2=1.0D3
            RETURN
         ENDIF
      ENDDO 
!
!  If we reach here then the two structures are permutational isomers to within the tolerance TOLD.
!
      DIST2=SQRT(DIST2)
      WORSTRAD=SQRT(Q1(3*(JWORST-1)+1)**2+Q1(3*(JWORST-1)+2)**2+Q1(3*(JWORST-1)+3)**2)

      RETURN
      END
