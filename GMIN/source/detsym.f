C
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
      SUBROUTINE DETSYMMETRY(Q,HORDER,IT,PTEST)
      USE COMMONS, ONLY: DEBUG, NATOMS, ATMASS 

      IMPLICIT NONE
      INTEGER HORDER
      LOGICAL PTEST
      DOUBLE PRECISION Q(*)
      DOUBLE PRECISION IT(3,3)

      INTEGER IPRNT, NHCHECK
      DOUBLE PRECISION TOLD, TOLE
      LOGICAL TWOD

      INTEGER IBOT2,INDEX,LOOK,JJ,JPLANE,ICLIP,IAXORD,JX,IX,ISIGMV,IMNAX,ISTART,IK,IORBIT,
     1        NORBIT,IBOT,IDIR,IXX,J2,I,J,ISKIP,ICOMPX,JCOMP,JCOMPX,ISAXIS,IHIGHX,ICOMP2,
     2        NUMB,ISET,IERR,IIAX,ICOMPQ,ISINV,IDONE,ICOUNT,ICOMP,ILINEAR,IDGRP,IINV,MINORB,IBTOR,IREF,
     3        IROT,IDEGEN,IHIGH,J1,NORD(2*NATOMS+1),NORD2(2*NATOMS+1)
      DOUBLE PRECISION CM(3),IV(3,3),NEWQ(3*NATOMS),DOTOPT,ANGMAG
      DOUBLE PRECISION RM(3,3),QSORT(3*NATOMS),TATB(3), TEMP(3,3), ITX, ITY, ITZ
      DOUBLE PRECISION SCRATCH(9*NATOMS),MOLWT,QSAVE(3*NATOMS),LTOLE
      DOUBLE PRECISION DELX,DELY,DIST,BILEN,TEST,
     1                 DIP,ARGU,ANGLE2,ANGIAX,ANGMG,ANG,ZANG2,ZANGS,CMX,CMY,CMZ,RANG,
     2                 ATMP,XXX,X,Z,ZANG,DPROJ,ARG,ANGL,DTOR,ORDIS
      LOGICAL AGAIN
      CHARACTER(LEN=4) STRING, JNKSTR, SAVEGRP
      CHARACTER(LEN=87) ESTRING 
      CHARACTER(LEN=80) GPSTRING, NSTRING, FSTRING
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      CHARACTER(LEN=4) :: FPGRP, BPGRP, PGRP
      IBTOR(I,J)  = IOR(I,J)
C
      NHCHECK = 6
C     TOLD initial distance tolerance in symmetry subroutine - default 0.0001
      TOLD = 0.0001D0
C     TOLE initial tolerance for the difference in principal moments of inertia divided by the sum 
C          of the principal moments in symmetry subroutine - default 0.0001
      TOLE = 0.0001D0
      TWOD = .FALSE.


      IPRNT=131
      IPRNT=0
      MINORB=1
      LTOLE=TOLE
      FPGRP='   '
      PGRP= '   '
      BPGRP='   '
      DTOR=DACOS(-1.D0)/180.D0
      AGAIN=.FALSE.
C 
C INITIALIZE MILLIONS OF VARIABLES THAT CFT77 SAYS AREN'T INITIALIZED. 
C

      DO 10 J1=1,3*NATOMS
         QSAVE(J1)=Q(J1)
10    CONTINUE
651   IHIGH=0
      IDEGEN=0
      IROT=0
      IREF=0
      IINV=0
      IDGRP=0
      ILINEAR=0
      ICOMP=0
      ICOUNT=0
      IDONE=0
      ISINV=0
      ICOMPQ=0
      IIAX=0
      ORDIS=0
      IERR=0
      ISET=0
      NUMB=0
      ICOMP2=0
      IHIGHX=0
      ISAXIS=0
      JCOMPX=0
      JCOMP=0
      ICOMPX=0
      ISKIP=0
C 
C TRANSLATE TO CENTER OF MASS
C
      IF (IPRNT.GE.3) WRITE(*,20) (ATMASS(J),J = 1,NATOMS)
20    FORMAT(1X,F15.10)
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      iloop: DO I = 1,NATOMS
         CMX = ATMASS(I)*Q(3*I-2)+CMX
         CMY = ATMASS(I)*Q(3*I-1)+CMY
         CMZ = ATMASS(I)*Q(3*I)+CMZ
         MOLWT = MOLWT+ATMASS(I)
      ENDDO iloop
      CM(1) = CMX/MOLWT
      CM(2) = CMY/MOLWT
      CM(3) = CMZ/MOLWT
      DO I = 1,NATOMS
         DO J = 0,2
            Q(3*I-J) = Q(3*I-J)-CM(3-J)
         ENDDO
      ENDDO
      IF (IPRNT .GE. 4)WRITE(*,*)
     1           'After translation to center of mass coordinates '
      IF (IPRNT .GE. 4)WRITE(*,50)(Q(I),I = 1,3*NATOMS)
50    FORMAT(3(2X,F12.6))
C 
C Build, print out, and then diagonalize inertia tensor 
C     Build inertia tensor
C
      CALL INERTIA2(IT,Q)
C
C e/vectors from smallest to largest.
C
      IF (IPRNT.GE.4) WRITE(*,*) 'Inertia tensor'
      IF (IPRNT.GE.4) WRITE(*,50) ((IT(I,J),J = 1,3),I = 1,3)
C
C Diagonalize inertia tensor. The 0 flag reorders the e/values and
C e/vectors from smallest to largest.
C
      CALL EIG(IT,IV,3,3,0)
C 
C Check *now* for degeneracy of eigenvalues 
C
C  Check *now* for degeneracy of eigenvalues -- If present, then see if 
C  unique moment of inertia is along x.  If so, change this axis to z
C  by rotating the eigenvector matrix about y.  This guarantees that
C  the unique axis will lie along z.
C
C  The moments could be large, so this should be a different tolerance
C  to the distance checking, i.e. a relative tolerance.
C
      IF (DABS((IT(2,2)-IT(3,3))/(IT(1,1)+IT(2,2)+IT(3,3))).LT.LTOLE)THEN
         RANG=90.D0
         CALL ROTM(2,RANG,1,RM)
C        CALL MATMUL(TEMP,IV,RM,3,3,3,3,3,3)
         CALL MUL3(TEMP,IV,RM)
         DO 70 J1=1,3
            DO 60 J2=1,3
               IV(J2,J1)=TEMP(J2,J1)
60          CONTINUE
70       CONTINUE
         ATMP=IT(1,1)
         IT(1,1)=IT(3,3)
         IT(3,3)=ATMP
      ENDIF

      IF (IPRNT.GE.4) WRITE (*,*) ' points before reorientation'
      IF (IPRNT.GE.4) WRITE (*,50) (Q(I),I = 1,3*NATOMS)

      IF (IPRNT.GE.4) WRITE(*,*)' Diagonalized inertia tensor '
      IF (IPRNT.GE.4) WRITE(*,50) ((IT(I,J),J = 1,3),I = 1,3)
C     WRITE(*,'(3F20.10)') IT(1,1), IT(2,2), IT(3,3)
      IF (IPRNT.GE.4) WRITE(*,*)' Eigenvectors of inertia tensor ' 
      IF (IPRNT.GE.4) WRITE(*,50)((IV(I,J),J = 1,3),I = 1,3)
      CALL MATMULV(NEWQ,Q,IV,NATOMS,3,3)
      IF (IPRNT.GE.4) WRITE (*,*) ' Principal axis orientation for molecular system '
      IF (IPRNT.GE.4) WRITE (*,50) (NEWQ(I),I = 1,3*NATOMS)
C 
C  Print rotational constants; check handedness of in. axes; generate 
C  sorted coordinate vector 
C  
C
C  Print out the rotational constants from the principal 
C  moments of inertia.
C
      IF (.NOT.AGAIN) CALL ROTCON (IT, 1, IERR, PTEST)
      AGAIN=.TRUE.
C
C  Check handedness of inertial axes - this was wreaking havoc with
C  dihedral angles!!!  If dot product is negative; switch sign of
C  y axis (the choice of axis to change is arbitrary)
C
      CALL CROSSOPT(IV(1,1),IV(1,2),TATB,1)
      XXX=DOTOPT(TATB,IV(1,3),3)
C     WRITE(*,*)' INERTIAL AXIS DOT PRODUCT IS ',XXX
      IF (XXX.LT.0.D0) THEN
         DO 80 I=2,3*NATOMS-1,3
            NEWQ(I)=-NEWQ(I)
80       CONTINUE
      ENDIF
C
C  Want transformation matrix to be the identity in this case.
C  Generate well-defined sorted coordinate vector now.  Used in
C  symmetry evaluation routines which follow.
C
      IF (IPRNT.GE.4) WRITE(*,50) (NEWQ(J),J = 1,NATOMS*3)
      CALL SORTXYZ2(NEWQ,QSORT,NORD,TOLD)
      IF (IPRNT.GE.4) WRITE(*,*) ' SORTED COORDINATE VECTOR'
      IF (IPRNT.GE.4) WRITE(*,50) (QSORT(J),J = 1,NATOMS*3)

C 
C
C Check if point group is Abelian 
C 
C  Check to see if point group is abelian - if not, go through the
C  very hairy routine to determine point group and put in an
C  orientation that aces can deal with.
C
      X = -1.D0
      DO I = 1,3
         Z = IT(I,I)
         IF (DABS((Z-X)/(IT(1,1)+IT(2,2)+IT(3,3))).LT.LTOLE) 
     1                            IDEGEN=IDEGEN+1 
         X = Z
      ENDDO
      IF (IDEGEN.EQ.0) THEN
         IF (IPRNT .GE. 3)WRITE(*,100)
100      FORMAT(' The molecule belongs to an Abelian group.')
         GOTO 630
      ENDIF
C 
C 
C  Start of the big block for point groups containing degenrate IR's   
C
C
      IF (IDEGEN.EQ.1) THEN
C
C  If doubly degenerate, check to make sure that the principal axes
C  are parallel to x,y and z.  If not, rotate orientation to this point.
C
C  Use a different tolerance for zero moment of inertia checking. 
C
         FPGRP = '   '
         IF (IPRNT .GE. 3)WRITE(*,110)
110      FORMAT(' The molecule belongs to a point group with doubly degenerate representations.')
         DO 120 I=1,3
            IF (DABS(IT(I,I)).LT.5.0D-2) ILINEAR=1
120      CONTINUE
         IF (ILINEAR.EQ.1) THEN
            IF (IPRNT .GE. 3)WRITE(*,130)
130         FORMAT(' The molecule is linear.')
C           ILINEAR = 1
            GOTO 630
         ENDIF
C  
C  Determination of unique axis - check for rot. axes from C2 to NHCHECK 
C
         DO I = 2,NHCHECK
            ZANG = 360.D0/FLOAT(I)
            DO J = 1,3
               CALL ROTM(J,ZANG,1,RM)
               CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
               CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
               IF (IPRNT.GE.4) WRITE(*,*) ICOMP
               IF (TWOD.AND.(J.NE.3)) ICOMP=1
               IF (ICOMP.EQ.0) THEN
                  IHIGH = I
                  IHIGHX = J
C 
C  Stuff ordering into NORD2 
C
                  DO 140 IXX = 1,2*NATOMS
                     NORD2(IXX) = NORD(IXX)
140               CONTINUE
               ENDIF
            ENDDO
         ENDDO

C 
C IHIGH=0: accidental degeneracy  
C  Catch problem if IHIGH=0 because two of the inertia tensor
C  e/values just happen to be very close
C
            IF (IHIGH.EQ.0) THEN
               IF (PTEST) PRINT*,'Accidental degeneracy detected'
               LTOLE=LTOLE/10.0D0
               IF (LTOLE.GT.1.0D-7) THEN
                  GOTO 651
               ELSE
                IF (PTEST) WRITE(*,'(A)') ' symmetry> The full molecular point group is undetermined'
                 WRITE(GPSTRING,653) 
     1              ' symmetry> The full molecular point group is undetermined'
653              FORMAT(A80)
                 GOTO 652
               ENDIF
            ENDIF
C 

            IF (IPRNT .GE. 3)WRITE(*,160)IHIGH,IHIGHX
160         FORMAT(' symmetry> The highest order rotational axis is C',I2,' about ',I2)
C
C IF ((IHIGH.EQ.2).AND.(.NOT.TWOD)) THEN: highest order axis=2 + twofold C deg. 
C
C  If highest order axis is 2, and there is twofold degeneracy, the
C  group must be D2d.  Proceed to find the S4 axis which determines
C  the unique rotational axis.  This is subsequently used to rotate
C  molecule to a useful orientation. (NORD2 now contains an effective
C  permutation list).
C
            IF ((IHIGH.EQ.2).AND.(.NOT.TWOD)) THEN
               FPGRP = 'D2d'
               ZANG = 90.D0
               DO 180 IDIR = 1,3
                  CALL ROTM(IDIR,ZANG,1,RM)
                  CALL MATMULV(SCRATCH(NATOMS*3+1),NEWQ,RM,NATOMS,3,3)
                  CALL REFLECT(SCRATCH(NATOMS*3+1),SCRATCH,NATOMS,IDIR)
                  CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
                  IF (ICOMP.EQ.0) THEN
                     IHIGHX = IDIR
                     DO 170 I = 1,2*NATOMS
                        NORD2(I) = NORD(I)
170                  CONTINUE
                  ENDIF
180            CONTINUE
               IF (IPRNT.GE.3)WRITE(*,190)IHIGHX
190            FORMAT(' The unique rotational axis is ',I1,'.')
            ENDIF
C
C
! IF (MOD(IHIGH,2).EQ.1) THEN 
!
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
C  have not put a vertical mirror plane coincident with the axis!
C  I first noticed this with an LJ55 cluster!
C  On the other hand, for D_{nd} groups with odd n we actually
C  may not have an orbit of size equal to the order of the principal
C  axis - there is sure to be one of twice the size though. Here
C  the original algorithm would have worked.
C  The principal axis order is IHIGH - this is only needed if
C  IHIGH is odd.
C
            IF (MOD(IHIGH,2).EQ.1) THEN
               IF (IPRNT.GE.4) WRITE(*,*)' Before rotation'
               IF (IPRNT.GE.4) WRITE(*,50) (NEWQ(J),J = 1,NATOMS*3)
               SCRATCH(1:3*NATOMS*3)=0.0D0
               DO J=1,NATOMS
                  NORD2(J)=J
                  IBOT=3*J-2
                    SCRATCH(J)=ATMASS(J)*DSQRT(DOTOPT(NEWQ(IBOT),NEWQ(IBOT),3))
               ENDDO
               CALL PIKSR2(NATOMS,SCRATCH,NORD2)
               IF (IPRNT.GE.3) WRITE(*,440)(SCRATCH(J),J=1,NATOMS)
               NORBIT=1
               IORBIT=1
               DO 210 J=2,NATOMS+1
                  IF (DABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLD) IORBIT=IORBIT+1
                  IF (DABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLD) THEN
                     NORD2(NATOMS+NORBIT)=IORBIT
                     IF (ATMASS(NORD2(J-1)).LT.0.1D0) THEN
                        SCRATCH(6*NATOMS+NORBIT)=0.0D0
                     ELSE
                          SCRATCH(6*NATOMS+NORBIT)=DABS(SCRATCH(J-1))/ATMASS(NORD2(J-1))
                     ENDIF
                     IF (J.LE.NATOMS) THEN
                        NORBIT=NORBIT+1
                        IORBIT=1
                     ENDIF
                  ENDIF
210            CONTINUE
C
C  NORD2(NATOMS+I) contains the size of orbit I. The atom numbers belonging to 
C  each orbit are stored in order in NORD2 up to element NORD2(NATOMS).
C
               IF (IPRNT.GE.130) THEN
                   WRITE(*,*)' SCRATCH VECTOR '
                   WRITE(*,*) (SCRATCH(IK),IK=1,NATOMS)
                   WRITE(*,*)' ORBIT DISTANCES: '
                   WRITE(*,*) (SCRATCH(6*NATOMS+IK),IK=1,NORBIT)
               ENDIF
               IF (IPRNT.GE.3) THEN
                  WRITE(*,470) NORBIT
                  NUMB=0
                  DO I=1,NORBIT
                     NUMB=NORD2(NATOMS+I)+NUMB
                     IREF=1+NUMB-NORD2(NATOMS+I)
C                    PRINT*,'I,NORBIT,IREF=',I,NORBIT,IREF
                     WRITE(*,220) I,NORD2(NATOMS+I),SCRATCH(6*NATOMS+I),(NORD2(IREF+J),J=0,NORD2(NATOMS+I)-1)
220                  FORMAT(2I3,F15.7,100I3)
                  ENDDO
               ENDIF
               NUMB=0
               IREF=0
               DO I=1,NORBIT
                  NUMB=NORD2(NATOMS+I)+NUMB
                  IF (NORD2(NATOMS+I).EQ.IHIGH) THEN
                     MINORB=NORD2(NATOMS+I)
                     IREF=1+NUMB-NORD2(NATOMS+I)
                     GOTO 250
                  ENDIF
               ENDDO
250            IF (IREF.NE.0) THEN
                  IF (IPRNT.GE.13) WRITE(*,260) MINORB, NORD2(IREF) 
260               FORMAT(' The principal axis has order ',i3,/,
     1            ' Atom number ',i3,' belongs to an orbit of this order and will be used as a reference.')
                  IF (IPRNT.GE.13) WRITE(*,520) (NORD2(IREF+J),J=1,MINORB-1)
                  ICOUNT=NORD2(IREF)
               ELSE
                  NUMB=0
                  IREF=0
                  DO 245 I=1,NORBIT
                     NUMB=NORD2(NATOMS+I)+NUMB
                     IF (NORD2(NATOMS+I).EQ.2*IHIGH) THEN
                        MINORB=NORD2(NATOMS+I)
                        IREF=1+NUMB-NORD2(NATOMS+I)
                        GOTO 255
                     ENDIF
245               CONTINUE
255               IF (IREF.NE.0) THEN
                     IF (IPRNT.GE.13) WRITE(*,265) IHIGH, NORD2(IREF)
265                  FORMAT(' The principal axis has order ',i3,/,
     1              ' Atom number ',i3,' belongs to an orbit of twice this order', 
     2              ' and will be used as a reference.')
                     IF (IPRNT.GE.13) WRITE(*,520) (NORD2(IREF+J),J=1,MINORB-1)
                     ICOUNT=NORD2(IREF)
                  ELSE
                     NUMB=0
                     IREF=0
                     DO 246 I=1,NORBIT
                        NUMB=NORD2(NATOMS+I)+NUMB
                        IF (NORD2(NATOMS+I).EQ.4*IHIGH) THEN
                           MINORB=NORD2(NATOMS+I)
                           IREF=1+NUMB-NORD2(NATOMS+I)
                           GOTO 256
                        ENDIF
246                  CONTINUE
256                  IF (IREF.EQ.0) THEN
                        IF (DEBUG) PRINT '(A,G20.10)','Orbit size inconsistent with axis order - try decreasing TOLD to ',
     &                             TOLD/10.0D0
                        TOLD=TOLD/10.0D0
                        IF (TOLD.GT.1.0D-7) THEN
                           GOTO 651
                        ELSE
                           IF (DEBUG) WRITE(*,'(A)') ' symmetry> The full molecular point group is undetermined'
                           WRITE(GPSTRING,653) ' symmetry> The full molecular point group is undetermined'
                           GOTO 652
                        ENDIF
                     ENDIF
                     IF (IPRNT.GE.13) WRITE(*,266) IHIGH, NORD2(IREF)
266                  FORMAT(' The principal axis has order ',I3,/,
     1                ' Atom number ',i3,' belongs to an orbit of four times this order'
     2               ,' and will be used as a reference.')
                     IF (IPRNT.GE.13) WRITE(*,520)
     1                   (NORD2(IREF+J),J=1,MINORB-1)
                     ICOUNT=NORD2(IREF)
                  ENDIF
               ENDIF
            ELSE
               DO 270 I= NATOMS,1,-1
                  IF (NORD2(I).NE.NORD2(I+NATOMS)) ICOUNT=NORD2(I)
270            CONTINUE
            ENDIF
!

C
C   Calculate angle between projection of this vector and one of the Cartesian axes 
C
            ISTART = 3*(ICOUNT-1)
            DPROJ = DSQRT(NEWQ(ISTART+1)**2+NEWQ(ISTART+2)**2+NEWQ(ISTART+3)**2 - NEWQ(ISTART+IHIGHX)**2)
C
C       
C   Find first Cartesian axis which is not the unique axis.                
C   (Has to be x or y)
C
            IMNAX = 2-IHIGHX/2
            ARG=NEWQ(ISTART+IMNAX)/DPROJ
            IF (DABS(ARG).GT.1.0D0) ARG=1.0D0*ARG/DABS(ARG)
            ANGL=DACOS(ARG)/DTOR
C                                                                       
C
C   Skip out of loop if alignment already satisfied    
C
            DO 280 I = 0,2
               X = ANGL-FLOAT(I*90)
               IF (DABS(X).LT.TOLD)ISKIP = 1
280         CONTINUE
            IF (ISKIP.EQ.1) GOTO 340
C       
C
C   Get sign of angle needed to rotate molecule into position   
C
C   I found a case where this didn't work. There are only two
C   possibilities (+ and -) so if one doesn't work, just use 
C   the other!
C
C           Z0 = 1
C           DO 210 I = 1,3
C              IF (I.EQ.IHIGHX) GOTO 210
C              Z0 = Z0*NEWQ(ISTART+I)
C210         CONTINUE
C           IF (Z0.GT.0.D0) ANGL = -ANGL
C
C   Rotate molecule into position 
C
            SCRATCH(1:NATOMS*2*3)=0.0D0
            CALL VADD(SCRATCH(1),SCRATCH(NATOMS*3+1),NEWQ,NATOMS*3,1)
            CALL ROTM(IHIGHX,ANGL,1,RM)
            CALL MATMULV(NEWQ,SCRATCH,RM,NATOMS,3,3)
            IF (IPRNT .GE. 3)WRITE(*,310)ANGL,IHIGHX
310         FORMAT(' Molecule rotated through ',F20.10,' degrees ',
     1           'about the ',i1,' axis.')
            IF (IPRNT.GE.4) WRITE(*,*)' After rotation:'
            IF (IPRNT.GE.4) WRITE(*,50) (NEWQ(J),J = 1,NATOMS*3)

            DO 320 J1=1,3
               IF (J1.EQ.IHIGHX) GOTO 320
               IF (J1.EQ.IMNAX) GOTO 320
               IF (DABS(NEWQ(ISTART+J1)).GT.1.0D-6) THEN
                  ANGL=-2.0D0*ANGL
                  SCRATCH(1:NATOMS*2*3)=0.0D0
                  CALL VADD(SCRATCH(1),SCRATCH(NATOMS*3+1),NEWQ,NATOMS*3,1)
                  CALL ROTM(IHIGHX,ANGL,1,RM)
                  CALL MATMULV(NEWQ,SCRATCH,RM,NATOMS,3,3)
                  IF (IPRNT.GE.3) WRITE(*,310) ANGL,IHIGHX
                  IF (IPRNT.GE.4) WRITE(*,*)' After rotation:'
                  IF (IPRNT.GE.4) WRITE(*,50) (NEWQ(J),J=1,NATOMS*3)  
               ENDIF
320         CONTINUE
C                                       
C   New reference structure, so we need to generate new sorted vector 
C
C
            CALL SORTXYZ2(NEWQ,QSORT,NORD,TOLD)
            IF (IPRNT.GE.4) WRITE(*,*)' New sorted coordinates'
            IF (IPRNT.GE.4) WRITE(*,50) (QSORT(J),J = 1,NATOMS*3)
C
C                                               
C   Check for sigma(v) planes now.              
C
            SCRATCH(1:NATOMS*3*2)=0.0D0
            ISIGMV=0
            DO 330 J=1,3
               IF (J.EQ.IHIGHX) GOTO 330
               CALL REFLECT(NEWQ,SCRATCH,NATOMS,J)
               CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMPX,TOLD)
               IF ((ICOMPX.EQ.0).AND.(.NOT.TWOD)) ISIGMV=1
330         CONTINUE
            IF (IPRNT.GT.4) PRINT*,'ISIGMV=',ISIGMV
C
C                                               
C   Check for S(2n) axis                        
C
340         ZANGS = 180.D0/FLOAT(IHIGH)
            CALL ROTM(IHIGHX,ZANGS,1,RM)
            CALL MATMULV(SCRATCH(NATOMS*3+1),NEWQ,RM,NATOMS,3,3)
            CALL REFLECT(SCRATCH(NATOMS*3+1),SCRATCH,NATOMS,IHIGHX)
            CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
            IF ((ICOMP.EQ.0).AND.(.NOT.TWOD)) THEN
               ISAXIS = 1
               IF (IPRNT .GE. 3) THEN
                  WRITE(*,350) 2*IHIGH
350               FORMAT(' S',I3,' axis exists.')
               ENDIF
            ENDIF
C 
C
C Look for perpendicular C2 axes            
C
C   Look for perpendicular C2 axes - check along Cartesian axes and
C   along angle bisectors of S2n operations,
C
C   This isn't good enough - in Dnd the C2 axes need not be in either
C   of the above orientations. Try a rotation that puts the reference
C   atom NORD2(IREF) and the next one in the orbit NORD2(IREF+1)
C
C   first rotate molecule by 1/2 S(2n) angle around IHIGH x and put in
C   SCRATCH(6*NATOMS+1)
C  
C
C           PRINT*,'looking for perpendicular C2 axes'
            ZANG = 180.D0
            ZANG2 = 180.D0/(2.D0*IHIGH)
C           PRINT*,'IHIGHX,ZANG2=',IHIGHX,ZANG2
            CALL ROTM(IHIGHX,ZANG2,1,RM)
!           SCRATCH=NEWQ*RM 
            CALL MATMULV(SCRATCH(6*NATOMS+1),NEWQ,RM,NATOMS,3,3)
C           PRINT*,'reference geometry without rotation and without sorting '
C           WRITE(*,'(3F15.7)') (NEWQ,J1=1,3*NATOMS)
C           PRINT*,'reference geometry after rotation through ',ZANG2,' before sorting'
C           WRITE(*,'(3F15.7)') (SCRATCH(6*NATOMS+J1),J1=1,3*NATOMS)
C 
C
C   Sort the reference geometry 
C
            SCRATCH(3*NATOMS+1:3*NATOMS+3*NATOMS)=0.0D0
            CALL SORTXYZ2(SCRATCH(6*NATOMS+1),SCRATCH(3*NATOMS+1),NORD,TOLD)
            DO 360 I = 1,3
               IF (I.EQ.IHIGHX) GOTO 360
C 
C   Checking Cartesian axes for both rotated and unrotated molecule             
C
               CALL ROTM(I,ZANG,1,RM)
               IF (IPRNT.GE.5) WRITE(*,50) ((RM(IX,JX),JX = 1,3),IX = 1,3)
               CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)                 ! unrotated molecule
               CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)        
               CALL MATMULV(SCRATCH,SCRATCH(6*NATOMS+1),RM,NATOMS,3,3)  ! rotated molecule
               CALL COMPARE3(SCRATCH,SCRATCH(3*NATOMS+1),NORD,ICOMPX,TOLD)
               IF ((ICOMP.EQ.0).AND.(.NOT.TWOD)) JCOMP = JCOMP+1
               IF ((ICOMPX.EQ.0).AND.(.NOT.TWOD)) JCOMPX = JCOMPX+1
360         CONTINUE
            IF (JCOMP.GE.1) THEN
               IF (IPRNT .GE. 3) WRITE(*,370)
370            FORMAT(' On-axis perpendicular C2 elements found.')
               IDGRP = 1
            ELSE IF (JCOMPX.GE.1) THEN
               IF (IPRNT .GE. 3) WRITE(*,380)
380            FORMAT(' Off-axis perpendicular C2 elements found.')
               IDGRP = 1
            ENDIF
C                                                                       
C  Another attempt which actually works out where the perpendicular C2 axis should be.
C  Assumes IHIGHX=3 so that the principal axis is the z axis.
C
            IF (IDGRP.NE.1) THEN
               DO J=1,MINORB-1
                  DELX=NEWQ(3*(NORD2(IREF+J)-1)+1)
                  DELY=NEWQ(3*(NORD2(IREF+J)-1)+2)
                  IF (DELX.NE.0.0D0) THEN
                     ZANG2=-ATAN(DELY/DELX)/2
                     ZANG2=ZANG2/DTOR
                     CALL ROTM(IHIGHX,ZANG2,1,RM)
                     CALL MATMULV(SCRATCH(6*NATOMS+1),NEWQ,RM,NATOMS,3,3)
C                    PRINT*,'reference geometry after rotation through ',ZANG2,' before sorting'
C                    WRITE(*,'(3F15.7)') (SCRATCH(6*NATOMS+J1),J1=1,3*NATOMS)
                     SCRATCH(3*NATOMS+1:3*NATOMS+3*NATOMS)=0.0D0
                     CALL SORTXYZ2(SCRATCH(6*NATOMS+1),SCRATCH(3*NATOMS+1),NORD,TOLD)
                     CALL ROTM(1,ZANG,1,RM)
                     CALL MATMULV(SCRATCH,SCRATCH(6*NATOMS+1),RM,NATOMS,3,3)  ! rotated molecule
                     CALL COMPARE3(SCRATCH,SCRATCH(3*NATOMS+1),NORD,ICOMPX,TOLD)
                     IF ((ICOMPX.EQ.0).AND.(.NOT.TWOD)) JCOMPX = JCOMPX+1
                     IF (JCOMPX.GE.1) THEN
                        IF (IPRNT .GE. 3) WRITE(*,380)
                        IDGRP=1
                     ELSE
                        IF (DELY.NE.0.0D0) ZANG2=ATAN(DELX/DELY)/DTOR
                     ENDIF
                  ENDIF 
               ENDDO
            ENDIF
C                                       
C   Check for sigma(h)                  
C
            IF (IDGRP.EQ.1) THEN
               CALL REFLECT(NEWQ,SCRATCH,NATOMS,IHIGHX)
               CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP2,TOLD)
               IF ((ICOMP2.EQ.0).AND.(.NOT.TWOD)) THEN
                  FPGRP = 'DNh'
                  IF (IPRNT .GE. 3)WRITE(*,390)
               ENDIF
            ELSE
               CALL REFLECT(NEWQ,SCRATCH,NATOMS,IHIGHX)
390            FORMAT(' Horizontal plane of symmetry found.')
               CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP2,TOLD)
               IF ((ICOMP2.EQ.0).AND.(.NOT.TWOD)) THEN
                  FPGRP = 'CNh'
                  IF (IPRNT .GE. 3)WRITE(*,390)
               ELSE
                  FPGRP = 'CN '
                  IF (ISAXIS.EQ.1) THEN
                     FPGRP = 'SN '
                     IHIGH = 2*IHIGH
                  ENDIF
                  DO J = 1,3
                     IF (J.EQ.IHIGHX) GOTO 400
                     CALL REFLECT(NEWQ,SCRATCH,NATOMS,J)
                     CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP2,TOLD)
400                  IF ((ICOMP2.EQ.0).AND.(.NOT.TWOD)) FPGRP = 'CNv'
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
         IF (IDEGEN.EQ.2) THEN
            IF (IPRNT.GE.3) WRITE(*,*) ' In td oh loop '
            IF (IPRNT.GE.3) WRITE(*,410)
410         FORMAT(' The molecule belongs to a point group with doubly', 
     1             ' and triply degenerate representations.')
C
C       
C
C   Now we know that the point group is either T, Td, O, Oh, T, I, Th
C   or Ih.  For these groups, however, the inertia tensor
C   is rotationally invariant; therefore different (and significantly
C   more complicated) methods need to be used.  This one seems to
C   be particularly efficient.
C
C   Check for inversion symmetry. This will tell us a great deal. 
C   (If inversion center present, then group is either Oh,Th OR Ih.
C
         DO 420 I = 1,NATOMS*3
            SCRATCH(I) = -NEWQ(I)
420      CONTINUE
         CALL COMPARE3(SCRATCH,QSORT,NORD,ISINV,TOLD)
C
C
C   Locate all orbits (groups of equivalent atoms at distance r) - 
C   do this by generating and then sorting the c.o.m-outer atom
C   distance matrix (mass-weighted).  Use NORD2 to keep track of
C   which atoms are in which orbit.
C
         SCRATCH(1:3*NATOMS*3)=0.0D0
         DO 430 J=1,NATOMS
            NORD2(J)=J
            IBOT=3*J-2
              SCRATCH(J)=ATMASS(J)*DSQRT(DOTOPT(NEWQ(IBOT),NEWQ(IBOT),3))
430      CONTINUE
         CALL PIKSR2(NATOMS,SCRATCH,NORD2)
         IF (IPRNT.GE.3) WRITE(*,440)(SCRATCH(J),J=1,NATOMS)
440      FORMAT(' ORBIT MATRIX:',(F10.6,/))
C
C 
C   Now count number of orbits; place pertinent info (distances) in top 
C   end of SCRATCH array and use top end of NORD2 to hold number
C   of centers/orbit.  This bookkeeping facilitates things down
C   the road.  Don't count dummy atoms!
C
         NORBIT=1
         IORBIT=1
C        CALL VADD(SCRATCH(NATOMS+1),SCRATCH(NATOMS+1),SCRATCH(1),NATOM
         DO 460 J=2,NATOMS+1
            IF (DABS(SCRATCH(J)-SCRATCH(J-1)).LT.TOLD) IORBIT=IORBIT+1
            IF (DABS(SCRATCH(J)-SCRATCH(J-1)).GT.TOLD) THEN
               NORD2(NATOMS+NORBIT)=IORBIT
               IF (ATMASS(NORD2(J-1)).LT.0.1D0) THEN
                  SCRATCH(6*NATOMS+NORBIT)=0.D0
               ELSE
                    SCRATCH(6*NATOMS+NORBIT)=DABS(SCRATCH(J-1))/ATMASS(NORD2(J-1))
               ENDIF
               IF (J.LE.NATOMS) THEN
                  NORBIT=NORBIT+1
                  IORBIT=1
               ENDIF
            ENDIF
460      CONTINUE

         IF (IPRNT.GE.130) THEN
            WRITE(*,*)' SCRATCH VECTOR '
            WRITE(*,*)(SCRATCH(IK),IK=1,NATOMS)
            WRITE(*,*)' ORBIT DISTANCES: '
            WRITE(*,*)(SCRATCH(6*NATOMS+IK),IK=1,NORBIT)
         ENDIF
C        NORBIT=NORBIT-1
C
C 
C   Debug print to make sure we've got all the orbits correct. 
C
         IF (IPRNT.GE.3) THEN
            WRITE(*,470) NORBIT
470         FORMAT(' Number of distinct planetary orbits: ',i3,/,  
     1             ' Summary of orbit sizes:')
            DO 490 I=1,NORBIT
               WRITE(*,480) I, NORD2(NATOMS+I)
480            FORMAT(2I3)
490         CONTINUE
         ENDIF
C 
C
C Find the smallest orbit of non-unit magnitude and proceed 
C
C  The programmer thinks that only the smallest orbit of non-unit
C  magnitude will be needed to identify the point group uniquely
C  by the following method.  Find this orbit and proceed.  Again,
C  be wary of dummy atoms!
C
C  SCRATCH(6*NATOMS+orbit number) will be the radius of the orbit.
C
         MINORB=NATOMS
         NUMB=0
         DO 500 I=1,NORBIT
            NUMB=NORD2(NATOMS+I)+NUMB
            IF (SCRATCH(6*NATOMS+I).LT.1.0D-8) GOTO 500
            IF (NORD2(NATOMS+I).LE.MINORB) THEN
               MINORB=NORD2(NATOMS+I)
               ORDIS=SCRATCH(6*NATOMS+I)
               IREF=1+NUMB-NORD2(NATOMS+I)
C              PRINT*,'I,NUMB,NORD2=',I,NUMB,NORD2(NATOMS+I)
            ENDIF
500      CONTINUE
         IF (IREF.EQ.0) THEN
            IF (DEBUG) PRINT*,'Accidental degeneracy detected'
            TOLE=TOLE/10.0D0
            IF (TOLE.GT.1.0D-7) THEN
               GOTO 651
            ELSE
               IF (DEBUG) WRITE(*,'(A)')' symmetry> The full molecular point group is undetermined'
               WRITE(GPSTRING,653) ' symmetry> The full molecular point group is undetermined'
               GOTO 652
            ENDIF
         ENDIF
         IF (IPRNT.GE.13) WRITE(*,510)MINORB,NORD2(IREF)
510      FORMAT(' The smallest orbit contains ',i3,' members.',/,
     1    ' Atom number ',i3,' belongs to this orbit and will be',
     2    ' used as a reference.')
         IF (IPRNT.GE.13) WRITE(*,520)(NORD2(IREF+J),J=1,MINORB-1)  
520      FORMAT(' Other members of this orbit are: ',100I4)
C
C 
C
C   Now use reference atom from orbit and loop over all other atoms.
C   find highest order rotational axis which permutes the two; this
C   uniquely identifies the highest order rotational axis present in
C   the molecule.  No kidding.
C
         SCRATCH(1:9*NATOMS)=0.0D0
         IBOT=3*NORD2(IREF)-2
C
C   Loop below is the "big one". It is fully executed two times
C   if and only if the FPGRP is of the T or I variety and no
C   C2 axis appropriate for orientation fudging was found after
C   identification of the group.
C
         IF (IPRNT.GE.13) WRITE(*,*)' JUST BEFORE BIG LOOP'
         IF (IPRNT.GE.13) WRITE(*,*)' MINORB SET TO ',MINORB
         DO 620 ICOUNT=1,2
            IF (IPRNT.GE.3)
     1         WRITE(*,*)' PASS ',ICOUNT,' THROUGH ROT. FINDER '
            DO 610 J=1,MINORB-1
               IF (IDONE.EQ.1) GOTO 610
               IF (IPRNT.GE.3)
     1         WRITE(*,*)' ATOMS ',NORD2(IREF),NORD2(IREF+J)
               IBOT2=3*NORD2(IREF+J)-2
C
C   First rotate molecule so that the interatomic distance vector
C   is parallel to the z-axis, with bisector along x-axis.
C   Structure now in SCRATCH(1). Here one has to allow for the
C   possibility that atoms I and J are 180 degrees apart; if this
C   is the case, loop through other atoms in orbit and find one
C   that is perpendicular (there has to be one for cubic groups)
C   and use this vector as the "bisector". If
C   it isn't done this way, the orientation fudging gets fudged up.
C
C              CALL CALCVEC(NEWQ(IBOT),NEWQ(IBOT2),SCRATCH(1),0)
               SCRATCH(1)=NEWQ(IBOT2)-NEWQ(IBOT)
               SCRATCH(2)=NEWQ(IBOT2+1)-NEWQ(IBOT+1)
               SCRATCH(3)=NEWQ(IBOT2+2)-NEWQ(IBOT+2)
               DIST=DSQRT(DOTOPT(SCRATCH(1),SCRATCH(1),3))
C
C  Bug detected via GMIN 1/1/07 DJW
C
C              CALL VADD(SCRATCH(1),NEWQ(IBOT),NEWQ(IBOT2),3*NATOMS,1)
               CALL VADD(SCRATCH(1),NEWQ(IBOT),NEWQ(IBOT2),3,1)
C
C   Deal with the linear problem right now. 
C
               BILEN=DOTOPT(SCRATCH(1),SCRATCH(1),3)
               IF (BILEN.LT.TOLD) THEN
                  ISET=0
                  DO 530 LOOK=1,MINORB-1
                     IF (ISET.EQ.1) GOTO 530
                     INDEX=3*NORD2(IREF+LOOK)-2
                     TEST=DOTOPT(NEWQ(IBOT),NEWQ(INDEX),3)
                     IF (IPRNT.GT.13) WRITE(*,*) TEST
                     IF (DABS(TEST).LT.TOLD) THEN
                        SCRATCH(1:3)=0.0D0
                        CALL VADD(SCRATCH(1),SCRATCH(1),NEWQ(INDEX),3,1) 
                        ISET=1
                     ENDIF
530               CONTINUE
               ENDIF
               CALL SIAZ(SCRATCH(1),RM,1)
               SCRATCH(1:3*NATOMS*3)=0.0D0
               CALL MATMULV(SCRATCH(3*NATOMS+1),NEWQ,RM,NATOMS,3,3)
C 
C   Now you have bisector along x. Rotate about x to bring the two
C   atoms into position parallel to z.
C
C
               DIP=DSQRT(DOTOPT(SCRATCH(3*NATOMS+IBOT+1),SCRATCH
     1            (3*NATOMS+IBOT+1),2))
               ARGU=SCRATCH(3*NATOMS+IBOT+2)/DIP
               IF (DABS(ARGU).GT.1.D0)ARGU=ARGU/DABS(ARGU)
               ANGLE2=-DACOS(ARGU)/DTOR
               IF (SCRATCH(3*NATOMS+IBOT+1).GT.0.D0)ANGLE2=-ANGLE2
               CALL ROTM(1,ANGLE2,1,RM)
               CALL MATMULV(SCRATCH,SCRATCH(3*NATOMS+1),RM,NATOMS,3,3)
C
C 
C   Now find new vector which bisects the two atoms in question - 
C   it had best be x! (Unless linear case where TATB is forced.)
C
               CALL VADD(TATB,SCRATCH(IBOT),SCRATCH(IBOT2),3,1)
               IF (ISET.EQ.1) THEN
                  TATB(1:3)=0.0D0
                  TATB(1)=1.D0
               ENDIF
C 
C 1. Test if the two atoms are connected by symmetry axes, looping over possible options 
C  (2,3,4 and 5).
C
               DO 600 IAXORD=5,2,-1
                  IF (IPRNT.GE.13) WRITE(*,*)' IN BIG CUBIC LOOP'
                  ICLIP=0
                  ANGIAX=360.D0/FLOAT(IAXORD)
                  ANGMG=ANGMAG(ORDIS,DIST,IAXORD,IERR)
                  IF (IPRNT.GE.3) WRITE(*,540)ANGMG
540               FORMAT(' ANGMAG angle is ',f8.4)
                  IF (IPRNT.GE.13) WRITE(*,550) IERR,IAXORD
550               FORMAT(' IERR is ',I2,' for ',I1)
                  IF (IERR.EQ.1) GOTO 600
                  IF (IPRNT.GE.13) THEN
                     WRITE(*,*)
     1                  ' IN ROTATIONAL LOOP, LOOKING FOR ORDER ',IAXORD  
                     WRITE(*,560)NORD2(IREF),NORD2(IREF+J)
560                  FORMAT(' PLAYING WITH ATOMS ',I3,' AND ',I3)
                        WRITE(*,*)' IBOT AND IBOT2 ARE ',IBOT,IBOT2
                     WRITE(*,*)' ORBIT DISTANCE IS ',ORDIS
                     WRITE(*,*)' INTERATOMIC DISTANCE IS ',DIST
                     WRITE(*,*)' MAGIC ANGLE IS ',ANGMG
                     WRITE(*,*)' MAGIC VECTOR IS ',(TATB(JJ),JJ=1,3)
                  ENDIF
570               CALL ROTM(3,ANGMG,1,RM)
                  IF (IPRNT.GE.13) 
     1               WRITE(*,50)(SCRATCH(JJ),JJ=1,3*NATOMS)
                  CALL MATMULV(SCRATCH(3*NATOMS+4),TATB,RM,1,3,3)

C                 IF (IPRNT.GE.13)WRITE(*,80)(SCRATCH(JJ),JJ=1,3*NATOMS)
C                 WRITE(*,*)' ROT. MAGIC VECTOR IS ',
C    1                        (SCRATCH(3*NATOMS+3+JJ),JJ=1,3)

                  CALL SIAZ(SCRATCH(3*NATOMS+4),RM,3)
                  SCRATCH(3*NATOMS+1:3*NATOMS+NATOMS*3*2)=0.0D0
                  CALL MATMULV(SCRATCH(3*NATOMS+1),SCRATCH,RM,NATOMS,3
     1                         ,3)  

C              IF (IPRNT.GE.13) WRITE(*,80) 
C    1             (SCRATCH(3*NATOMS+JJ),JJ=1,3*NATOMS)

                  CALL SORTXYZ2(SCRATCH(3*NATOMS+1),QSORT,NORD,TOLD)
                  CALL ROTM(3,ANGIAX,1,RM)
                  CALL MATMULV(SCRATCH(6*NATOMS+1),SCRATCH(3*NATOMS+1),               
     1                         RM,NATOMS,3,3)
C
C 
C   Point group determination happens here, as well as orientatioN 
C
                  CALL COMPARE3(SCRATCH(6*NATOMS+1),QSORT,NORD,IIAX,
     1                          TOLD) 
C 
C   Point group determination happens here, as well as orientation fudging. 
C
                  IF (IIAX.EQ.0) THEN
                     IF (IPRNT.GE.3) WRITE(*,580) IAXORD
580                  FORMAT(' Axis of order ',I1,' identified.')
                     IF (IAXORD.EQ.5) THEN
                        FPGRP='I  '
                        IF (ISINV.EQ.0)FPGRP='I h'
                     ENDIF
                     IF (IAXORD.EQ.4) THEN
                        FPGRP='O  '
                        IF (ISINV.EQ.0) FPGRP='O h'
C
C 
C   Orientation is good.  Put structure into NEWQ and leave loop. 
C
                        NEWQ(1:3*NATOMS)=0.0D0
                        CALL VADD(NEWQ,NEWQ,SCRATCH(3*NATOMS+1),3*NATOMS
     1                                                               ,1)
                        CALL SORTXYZ2(NEWQ,QSORT,NORD,TOLD)
                        IDONE=1
                        GOTO 610
                     ENDIF
                     IF (IAXORD.EQ.3) THEN
                        IF (FPGRP.EQ.'   ') THEN
                           FPGRP='T  '
                           IF (ISINV.EQ.0)FPGRP='T h'
                        ENDIF
                     ENDIF
                     IF (IAXORD.EQ.2) THEN
                        IF (FPGRP.NE.'   '.AND.FPGRP.NE.'O h'.AND.          
     1                              FPGRP.NE.'O  ') THEN
C 
C
C   Have a C2 axis for I or T group. This is very good. You want this
C   to be the orientation for PITZER. Td and T can be differentiated
C   here.  Do this now if appropriate. Then put structure into NEWQ and
C   leave loop if it is an I group. Otherwise, need to continue through
C   (since there might be a higher order axis down the line.) If
C   pass two, however, then exit for T since it must be the correct group.
C
C 
                           IF (ICOUNT.EQ.2.AND.(FPGRP.EQ.'T  '.OR.FPGRP.EQ.'T h')) THEN
                              DO 590 JPLANE=1,3
                                 CALL REFLECT(SCRATCH(3*NATOMS+1),SCRATCH,NATOMS,JPLANE)
                                 CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMPQ,TOLD)
                                 IF (ICOMPQ.EQ.0.AND.ISINV.NE.0) FPGRP='T d'
                                 NEWQ(1:3*NATOMS)=0.0D0
                                 CALL VADD(NEWQ,NEWQ,SCRATCH(3*NATOMS+1),3*NATOMS,1)
                                 CALL SORTXYZ2(NEWQ,QSORT,NORD,TOLD)
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
            IF (PTEST) PRINT*,'Accidental degeneracy detected'
            LTOLE=LTOLE/10.0D0
            IF (LTOLE.GT.1.0D-7) THEN
               GOTO 651
            ELSE
               IF (PTEST) WRITE(*,'(A)') ' symmetry> The full molecular point group is undetermined'
               WRITE(GPSTRING,653) ' symmetry> The full molecular point group is undetermined'
               GOTO 652
            ENDIF
         ENDIF
      ENDIF
C
C   Cubic point group now determined.
C 
C   *** This is the end of the non-abelian point group block. ***
C
C   Check symmetry operations belonging to abelian groups 
C
630   ANG = 180.D0
      IREF=0
      DO 660 I = 3,1,-1
         CALL REFLECT(NEWQ,SCRATCH,NATOMS,I)
         CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
         IF ((ICOMP.EQ.0).AND.(.NOT.TWOD)) THEN
            IF (IPRNT.GE.3) WRITE(*,640) I
640         FORMAT(' Reflection in plane ',i2,' is a valid ',
     1           ' symmetry operation.')
            IREF=IBTOR(IREF,2**I/2)
         ENDIF
         CALL ROTM(I,ANG,1,RM)
         IF (IPRNT.GE.5) WRITE(*,50) ((RM(IX,JX),JX = 1,3),IX = 1,3)
         CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
         CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
         IF (TWOD.AND.(I.NE.3)) ICOMP=1
         IF (ICOMP.EQ.0) THEN
            IF (IPRNT.GE.3) WRITE(*,650) I
650         FORMAT(' Rotation about ',i2,' is a valid symmetry operation ')
            IROT = IBTOR(IROT,2**I/2)
         ENDIF
660   CONTINUE
      DO 670 I = 1,NATOMS*3
         SCRATCH(I)=-NEWQ(I)
670   CONTINUE
      CALL COMPARE3(SCRATCH,QSORT,NORD,ICOMP,TOLD)
C
C  Even DNd do not have the inversion operation - bug fixed DJW 13/5/92
C
      IF (((ISAXIS.EQ.1).OR.(ICOMP.EQ.0)).AND.(.NOT.TWOD)) THEN
         IF (IPRNT.GE.3) WRITE(*,680)
680      FORMAT(' The molecule possesses an inversion center. ')
         IF (ILINEAR.EQ.1)FPGRP = 'DXh'
         IF (FPGRP.EQ.'   '.AND.IDGRP.EQ.1)FPGRP = 'DNd'
         IINV=1
      ELSE
         IF (ILINEAR.EQ.1)FPGRP = 'CXv'
         IF (FPGRP.EQ.'   '.AND.IDGRP.EQ.1)FPGRP='DN '
      ENDIF
C
C   Determine largest abelian subgroup of molecule 
C
      IF (IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP = 'C1 '
      IF (IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.1) PGRP = 'C i'
      IF (IROT.EQ.0.AND.IREF.NE.0.AND.IINV.EQ.0) PGRP = 'C s'
      IF (IROT.NE.0.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP = 'C2 '
      IF (IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.0) PGRP = 'C2v'
      IF (IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.1) PGRP = 'C2h'
      IF (IROT.EQ.7.AND.IREF.EQ.0.AND.IINV.EQ.0) PGRP = 'D2 '
      IF (IROT.EQ.7.AND.IREF.EQ.7.AND.IINV.EQ.1) PGRP = 'D2h'
      IF (IPRNT.GE.4) WRITE(*,690) IREF,IROT,IINV
690   FORMAT ('Symmetry bits: ',3(1X,I3))
C 
C 
C
C   Rotate to default symmetry frame if needed  
C- completely done by
C   brute force. No cute algorithm used here.
C
C   ACES can not deal with all abelian point groups as sad as thiS
C   may seem. Allow for this now, and hope that someday this isn't
C   necessary.
C
      BPGRP = PGRP
      SCRATCH(1:NATOMS*3)=0.0D0
      IF (PGRP.EQ.'C1 '.OR.PGRP.EQ.'C i'.OR.PGRP.EQ.'C2 '.OR.
     1     PGRP.EQ.'D2 '.OR.PGRP.EQ.'D2h') THEN
         IF (PGRP.NE.'D2h')PGRP = 'C1 '
         CALL VADD(Q,NEWQ,SCRATCH(1),NATOMS*3,1)
         GOTO 740
      ENDIF
      IF (PGRP.EQ.'C2h') THEN
         PGRP = 'C s'
         CONTINUE
      ENDIF
C 
C
C   C2v
C
      IF (PGRP.EQ.'C2v'.OR.PGRP.EQ.'C2h') THEN
C 
C
C   If x is rotation axis - switch to z - and put largest moment of
C   inertia around x - eigenvectors returned from EIG are sorted
C   highest to lowest for D2h, we take standard orientation such
C   that Ix>Iy>Iz
C
         IF (MOD(IROT,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J)
            ENDDO
         ELSE IF (MOD(IROT/2,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J)
            ENDDO
         ELSEIF (MOD(IROT/4,2).EQ.1) THEN
            CALL VADD(Q,NEWQ,SCRATCH(1),NATOMS*3,1)
         ENDIF
C 
      ELSE
C 
         IF (MOD(IREF,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J)
            ENDDO
         ELSEIF (MOD(IREF/2,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J)
            ENDDO
         ELSEIF (MOD(IREF/4,2).EQ.1) THEN
            CALL VADD(Q,NEWQ,SCRATCH(1),3*NATOMS,1)
         ENDIF
C 
      ENDIF
740   CONTINUE

C     IF (PTEST) WRITE(*,770)
770   FORMAT(80('*'))
C
C  This should cope with anything we are likely to come up against,
C  until Mark tries a ring with 100 atoms!
C
C 
      IF (IHIGH.EQ.0) THEN
         WRITE(JNKSTR,777) FPGRP(1:3)
777      FORMAT(A3)
      ELSE IF (IHIGH.LT.10) THEN
         WRITE(JNKSTR,775) FPGRP(1:1),IHIGH,FPGRP(3:3)
775      FORMAT(A1,I1,A1)
      ELSE
         WRITE(JNKSTR,776) FPGRP(1:1),IHIGH,FPGRP(3:3)
776      FORMAT(A1,I2,A1)
      ENDIF
      FPGRP = JNKSTR
      IF (IDEGEN.GT.0) THEN
         IF (TWOD) THEN
            IF (PTEST) WRITE(*,'(A,A4,A)') ' symmetry> The 2D rotational subgroup is ',FPGRP,'.'
            WRITE(GPSTRING,'(A,A4,A)') ' symmetry> The 2D rotational subgroup is ',FPGRP,'.'
         ELSE
            IF (PTEST) WRITE(*,750) FPGRP
            WRITE(GPSTRING,750) FPGRP
         ENDIF
750      FORMAT(' symmetry> The full molecular point group is ',A4,'.')
         SAVEGRP=FPGRP
      ENDIF
      IF (IDEGEN.EQ.0) THEN
         IF (TWOD) THEN
            IF (PTEST) WRITE(*,'(A,A4,A)') ' symmetry> The 2D rotational subgroup is ',STRING(BPGRP,IHIGH),'.'
            WRITE(GPSTRING,'(A,A4,A)') ' symmetry> The 2D rotational subgroup is ',STRING(BPGRP,IHIGH),'.'
         ELSE
            IF (PTEST) WRITE(*,750) STRING(BPGRP,IHIGH)
            WRITE(GPSTRING,750) STRING(BPGRP,IHIGH)
         ENDIF
         SAVEGRP=STRING(BPGRP,IHIGH)
      ENDIF
      IF (PTEST.AND.(.NOT.TWOD)) WRITE(*,790) BPGRP
      IF (PTEST) WRITE(*,760) TOLD, LTOLE
760   FORMAT(' symmetry> Distance tolerance=',F12.5,' Inertia tolerance=',F12.5)  
790   FORMAT(' symmetry> The largest Abelian subgroup of the full molecular point group is ',A4,'.')

C     CALL GEOMOUT
C     CALL ADM(Q)

652   CONTINUE
C
C 
C  Turn off reorientation of system after the first step for INR=5 
C
C      IF (INR.NE.5) THEN
C         DO 800 J1=1,3*NATOMS
C            Q(J1)=QSAVE(J1)
C800      CONTINUE
C      ELSE
C         INR=INR-5
C      ENDIF
! 

!IF (SAVEGRP(1:3).EQ. ...) HORDER=... 
      HORDER=1
      IF (SAVEGRP(1:3).EQ.'C i') HORDER=2
      IF (SAVEGRP(1:3).EQ.'C s') HORDER=2
      IF (SAVEGRP(1:3).EQ.'C2 ') HORDER=2
      IF (SAVEGRP(1:3).EQ.'C3 ') HORDER=3
      IF (SAVEGRP(1:3).EQ.'C4 ') HORDER=4
      IF (SAVEGRP(1:3).EQ.'C5 ') HORDER=5
      IF (SAVEGRP(1:3).EQ.'C6 ') HORDER=6
      IF (SAVEGRP(1:3).EQ.'C7 ') HORDER=7
      IF (SAVEGRP(1:3).EQ.'C8 ') HORDER=8
      IF (SAVEGRP(1:3).EQ.'D2 ') HORDER=4
      IF (SAVEGRP(1:3).EQ.'D3 ') HORDER=6
      IF (SAVEGRP(1:3).EQ.'D4 ') HORDER=8
      IF (SAVEGRP(1:3).EQ.'D5 ') HORDER=10
      IF (SAVEGRP(1:3).EQ.'D6 ') HORDER=12
      IF (SAVEGRP(1:3).EQ.'D7 ') HORDER=14
      IF (SAVEGRP(1:3).EQ.'D8 ') HORDER=16
      IF (SAVEGRP(1:3).EQ.'C2v') HORDER=4
      IF (SAVEGRP(1:3).EQ.'C3v') HORDER=6
      IF (SAVEGRP(1:3).EQ.'C4v') HORDER=8
      IF (SAVEGRP(1:3).EQ.'C5v') HORDER=10
      IF (SAVEGRP(1:3).EQ.'C6v') HORDER=12
      IF (SAVEGRP(1:3).EQ.'C7v') HORDER=14
      IF (SAVEGRP(1:3).EQ.'C8v') HORDER=16
      IF (SAVEGRP(1:3).EQ.'C2h') HORDER=4
      IF (SAVEGRP(1:3).EQ.'C3h') HORDER=6
      IF (SAVEGRP(1:3).EQ.'C4h') HORDER=8
      IF (SAVEGRP(1:3).EQ.'C5h') HORDER=10
      IF (SAVEGRP(1:3).EQ.'C6h') HORDER=12
      IF (SAVEGRP(1:3).EQ.'C7h') HORDER=14
      IF (SAVEGRP(1:3).EQ.'C8h') HORDER=16
      IF (SAVEGRP(1:3).EQ.'D2h') HORDER=8
      IF (SAVEGRP(1:3).EQ.'D3h') HORDER=12
      IF (SAVEGRP(1:3).EQ.'D4h') HORDER=16
      IF (SAVEGRP(1:3).EQ.'D5h') HORDER=20
      IF (SAVEGRP(1:3).EQ.'D6h') HORDER=24
      IF (SAVEGRP(1:3).EQ.'D7h') HORDER=28
      IF (SAVEGRP(1:3).EQ.'D8h') HORDER=32
      IF (SAVEGRP(1:3).EQ.'D2d') HORDER=8
      IF (SAVEGRP(1:3).EQ.'D3d') HORDER=12
      IF (SAVEGRP(1:3).EQ.'D4d') HORDER=16
      IF (SAVEGRP(1:3).EQ.'D5d') HORDER=20
      IF (SAVEGRP(1:3).EQ.'D6d') HORDER=24
      IF (SAVEGRP(1:3).EQ.'D7d') HORDER=28
      IF (SAVEGRP(1:3).EQ.'D8d') HORDER=32
      IF (SAVEGRP(1:3).EQ.'S4 ') HORDER=4
      IF (SAVEGRP(1:3).EQ.'S6 ') HORDER=6
      IF (SAVEGRP(1:3).EQ.'S8 ') HORDER=8
      IF (SAVEGRP(1:3).EQ.'T  ') HORDER=12
      IF (SAVEGRP(1:3).EQ.'T d') HORDER=24
      IF (SAVEGRP(1:3).EQ.'T h') HORDER=24
      IF (SAVEGRP(1:3).EQ.'O  ') HORDER=24
      IF (SAVEGRP(1:3).EQ.'O h') HORDER=48
      IF (SAVEGRP(1:3).EQ.'I  ') HORDER=60
      IF (SAVEGRP(1:3).EQ.'I h') HORDER=120
      IF (SAVEGRP(1:3).EQ.'DXh') HORDER=2
      IF (SAVEGRP(1:3).EQ.'CXv') HORDER=2
! 
      IF (PTEST) THEN
         IF (TWOD) THEN
            WRITE(*,'(A,I6)') ' symmetry> Order of 2D rotational subgroup=',HORDER
         ELSE
            WRITE(*,'(A,I6)') ' symmetry> Order of full point group=',HORDER
         ENDIF
      ENDIF
C     IF (PTEST) WRITE(*,770)

      RETURN
      END
