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
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!=============================================================
!   All routines in this file were implemented by
!   Dmitri Schebarchov (ds656).
!=============================================================
!
SUBROUTINE PGSYM(N_IN,X_IN,L_IN,TOLS_IN,IMODE_IN,ML_IN, &
     CNTR_OUT,NSYMOPS_OUT,SYMOPS_OUT,FPG_OUT)
  !
  ! Detect approximate point-group symmetry for N_IN particles with 
  ! coordinates X_IN and labels L_IN, using 3 tolerances specified 
  ! in TOLS_IN(1:3).
  !
  ! TOLS_IN(1) is used for detecting degenerate principal moments. 
  ! TOLS_IN(2) is a distance threshold for testing sym. ops. 
  ! TOLS_IN(3) is for testing if a sym. op. matrix is new.
  !
  ! IMODE_IN specifies if sym. op. matrices are to be transformed 
  ! to the appropriate reference frame (IMODE_IN=1), and if the 
  ! set of sym. ops. is to be completed (IMODE_IN > 1). 
  ! IMODE_IN<1 is adequate only for classification purposes, 
  ! when the actual sym. ops. are not needed.
  !
  ! Output the centroid (i.e. centre of coordinates) the symmetry 
  ! operations in matrix representation, and a string describing 
  ! the detected point group (FPG).
  !
  ! NOTE: Sym. ops. are diagnosed using just the majority species, 
  ! but candidate sym. ops. are tested on all species. The tests 
  ! are sensitive to the supplied atomic labels (L_IN): atoms with 
  ! different labels cannot be symmetry equivalent!
  !
  ! All core routines are encapsulated in module PGSYMMOD.
  !
  ! The general outline of the procedure is loosely based on the
  ! pymatgen.symmetry.pointgroup module:
  !
  ! 1. Centre the molecule around the majority's center of mass.
  ! 2. Compute the inertia tensor and the e-vals and e-vecs.
  ! 3. Handle the symmetry detection based on eigenvalues.
  !
  !    a. Linear molecules have one zero e-val, amd the two non-
  !       zero ones ought to be degenerate. Possible symmetry
  !       operations are C*v or D*v.
  !    b. Asymetric top molecules have all different e-vals. The
  !       maximum rotational symmetry in such molecules is 2-fold.
  !    c. Symmetric top molecules have 1 unique e-val, giving a
  !       unique rotation axis. All axial point groups are possible
  !       except the cubic groups (T & O) and I.
  !    d. Spherical top molecules have all three eigenvalues equal. 
  !       They have the rare T, O or I point groups.
  !
  ! See http://pymatgen.org/index.html for more details.
  !
  USE COMMONS, ONLY : MYUNIT, DEBUG
  USE PGSYMMOD
  !
  IMPLICIT NONE
  !
  ! Passed variabels
  INTEGER, INTENT(IN) :: N_IN ! Atom count
  DOUBLE PRECISION, INTENT(IN) :: X_IN(3*N_IN) ! flattened coords.
  INTEGER, INTENT(IN) :: L_IN(N_IN) ! Atom labels
  DOUBLE PRECISION, INTENT(IN) :: TOLS_IN(3) ! tolerances
  INTEGER, INTENT(IN) :: IMODE_IN ! mode specifier
  INTEGER, INTENT(IN) :: ML_IN ! The majority label/species
  DOUBLE PRECISION, INTENT(OUT) :: CNTR_OUT(3) ! centroid
  INTEGER, INTENT(OUT) :: NSYMOPS_OUT ! numuber of matrices
  DOUBLE PRECISION, INTENT(OUT) :: SYMOPS_OUT(120,3,3) ! matrices
  CHARACTER(LEN=4), INTENT(OUT) :: FPG_OUT ! Point group label
  !
  ! Variables for internal use
  INTEGER :: I,J,K
  !
  NAT=N_IN
  TOL(1:3) = TOLS_IN(1:3)
  !
  CALL ALLOCATE_ALLX_AND_INI()
  !
  MAJORITY(0:NAT) = 0
  DO I=1,NAT
     LABELS(I) = L_IN(I)
     IF(L_IN(I) == ML_IN) THEN
        MAJORITY(0) = MAJORITY(0) + 1
        MAJORITY(MAJORITY(0)) = I
     ENDIF
  ENDDO
  IF(DEBUG) THEN
     WRITE(MYUNIT,'(A,I6)') &
          'pgsym> Majority species population: ', MAJORITY(0)
  ENDIF
  !
  ! The TRY loop is for potential failures due to loose
  ! tolerances, which are tightened by an order of magnitude
  ! with each iteration of the loop.
  KEEPTRYING = .TRUE.
  TRY: DO WHILE(KEEPTRYING) 
     KEEPTRYING = .FALSE.
     !
     ! Reset all the necessary symmetry elements
     NSYMOPS = 0
     SYMOPS(:,:,:) = 0.0D0
     NROTAXES = 0
     ROTAXES(:,:) = 0.0
     ROTAXESORDER(:) = 0
     !
     ! Find the centroid for the majority species
     CALL CALC_CENTROID(NAT,X_IN,MAJORITY,CNTR_OUT)
     ! Translate X_IN to CNTR_OUT at origin and store XREF as
     ! a set of COLUMN vectors (i.e. 1st index spans XYZ)
     IF(DEBUG) THEN
        WRITE(MYUNIT,'(A,3(1X,F16.8))') &
             'pgsym> Centroid: ', (CNTR_OUT(J), J=1,3)
        WRITE(MYUNIT,'(A,3(1X,F16.8))') &
             'pgsym> Tolerances: ', (TOL(J), J=1,3)
     ENDIF
     DO I=1,NAT
        DO J=1,3
           XREF(J,I) = X_IN(3*(I-1)+J) - CNTR_OUT(J)
        ENDDO
     ENDDO
     !
     CALL PGSYM_DIAGNOSE()
     !
     IF(IDEGEN==0) THEN 
        WRITE(MYUNIT,'(A)') 'pgsym> Asymmetric top diagnosed.'
        CALL PROC_ASYM_TOP()
     ELSEIF(IDEGEN==1) THEN 
        IF(LINEAR) THEN
           WRITE(MYUNIT,'(A)') 'pgsym> Linear structure diagnosed.'
           CALL PROC_LINEAR() 
        ELSE
           WRITE(MYUNIT,'(A)') 'pgsym> Symmetric top diagnosed.'
           CALL PROC_SYM_TOP() 
        ENDIF
     ELSEIF(IDEGEN==2) THEN
        WRITE(MYUNIT,'(A)') 'pgsym> Spherical top diagnosed.'
        CALL PROC_SPH_TOP()
     ELSEIF(IDEGEN==3 .AND. TOL(1) > 1.0D0-6) THEN
        WRITE(MYUNIT,'(A)') &
             'pgsym> Triple degeneracy diagnosed. Dodgy tolerance?'
        WRITE(MYUNIT, '(A)') &
             'pgsym> Tightening TOL(1).'
        TOL(1) = TOL(1)/10.0D0
        KEEPTRYING=.TRUE.
     ELSE
        WRITE(MYUNIT,'(A)') 'pgsym> ERROR: Bad diagnosis.'
        WRITE(MYUNIT,'(A,I2)') 'pgsym> IDEGEN=', IDEGEN
        STOP
     ENDIF
     !
  ENDDO TRY
  !
  WRITE(MYUNIT,'(A, A, A, I3, A)') &
       'pgsym> Identified point group ', FPG, &
       ' from ', NSYMOPS,' sym. ops.'
  !
  IF(IMODE_IN > 0) THEN
     CALL TRANSFORM_SYMOPS()
     IF(IMODE_IN > 1) THEN
        CALL COMPLETE_SYMOPS()
        WRITE(MYUNIT, '(A, I3, A)') &
             'pgsym> Completing the set yields ', &
             NSYMOPS,' sym. ops.'
     ENDIF
  ENDIF
  !
  NSYMOPS_OUT = NSYMOPS
  SYMOPS_OUT(:,:,:) = SYMOPS(:,:,:)
  FPG_OUT = FPG
  !
  CALL DEALLOCATE_ALLX()
  !
  RETURN
  !
END SUBROUTINE PGSYM
!
!==================================================================
!!!!!!!! Now come more self-contained subroutines !!!!!!!!!!!!!!!!!
!==================================================================
!
SUBROUTINE CALC_CENTROID(N,X,LIST,CNTR)
  !
  ! Compute the centroid of 3*N coordinated listed in X.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N, LIST(0:N)
  DOUBLE PRECISION, INTENT(IN) :: X(3*N)
  DOUBLE PRECISION, INTENT(OUT) :: CNTR(3)
  !
  INTEGER :: I,J,K,L
  !
  ! Find the centre of coordinates (centroid) for X
  CNTR(1:3) = 0.0D0
  DO L=1,LIST(0)
     I=LIST(L)
     K=3*(I-1)
     DO J=1,3
        CNTR(J) = CNTR(J) + X(K+J)
     ENDDO
  ENDDO
  CNTR(1:3) = CNTR(1:3)/DBLE(LIST(0))
  !
  RETURN
  !
END SUBROUTINE CALC_CENTROID
!
SUBROUTINE CALC_NUIT(N,X,LIST,NUIT)
  !
  ! Compute the normalised, uniform inertia tensor (UIT) 
  ! for the 3*N coordinates listed in in X. The centroid 
  ! is assumed to be at the origin, and the normalisation 
  ! is achieved via devision by the trace.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N,LIST(0:N)
  DOUBLE PRECISION, INTENT(IN) :: X(3,N)
  DOUBLE PRECISION, INTENT(OUT) :: NUIT(3,3)
  !
  INTEGER :: I,J,K,L
  DOUBLE PRECISION :: X2(3),R2,NORM
  !
  NUIT(:,:) = 0.0D0
  !
  ! Compute elements only in the upper triangle
  DO L=1,LIST(0)
     I=LIST(L)
     X2(1:3) = 0.0D0
     R2=0.0D0
     DO J=1,3
        X2(J) = X(J,I)**2
        R2 = R2 + X2(J)
     ENDDO
     DO J=1,3
        NUIT(J,J) = NUIT(J,J) + R2 - X2(J) ! diagonal
        DO K=J+1,3 ! off-diagonal
           NUIT(J,K) = NUIT(J,K) - X(J,I)*X(K,I)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Normalise and impose symmetry for the remaining elements
  NORM = NUIT(1,1) + NUIT(2,2) + NUIT(3,3) ! Trace (invariant!)
  DO J=1,3
     NUIT(J,J) = NUIT(J,J)/NORM
     DO K=J+1,3
        NUIT(J,K) = NUIT(J,K)/NORM
        NUIT(K,J) = NUIT(J,K)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE CALC_NUIT
!
SUBROUTINE PRODMM(NA,NB,NC,A,B,C) !--------------------------------
  !
  ! C=AB, i.e. C_ij=A_ik*B_kj, where i=[1,NA]; k=[1,NB]; j=[1,NC]
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NA,NB,NC
  DOUBLE PRECISION, INTENT(IN) :: A(NA,NB), B(NB,NC)
  DOUBLE PRECISION, INTENT(OUT) :: C(NA,NC)
  !
  INTEGER :: I,J,K
  DOUBLE PRECISION :: X
  !
  DO I=1,NA
     DO J=1,NC
        X = 0.0D0
        DO K=1,NB
           X = X + A(I,K)*B(K,J)
        ENDDO
        C(I,J) = X
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE PRODMM
!
SUBROUTINE PRODMTM(NA,NB,NC,A,B,C) !-------------------------------
  !    t 
  ! C=A B , i.e. C_ij=A_ki*B_kj
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NA,NB,NC
  DOUBLE PRECISION, INTENT(IN) :: A(NA,NB), B(NB,NC)
  DOUBLE PRECISION, INTENT(OUT) :: C(NA,NC)
  !
  INTEGER :: I,J,K
  DOUBLE PRECISION :: X
  !
  DO I=1,NA
     DO J=1,NC
        X = 0.0D0
        DO K=1,NB
           X = X + A(K,I)*B(K,J)
        ENDDO
        C(I,J) = X
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE PRODMTM
!
SUBROUTINE FIND_RMAT(N,X,X0,RMAT)
  !
  ! Find matrix RMAT that maps a reference structure X0 to 
  ! a structure X, i.e. X = RMAT * X0 for all N
  ! 3-dimensional atomic position vectors of X and X0.
  ! A more concise implementation of that in reorient.f90.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N ! particles
  DOUBLE PRECISION, INTENT(IN) :: X(3,N),X0(3,N) ! 2 sets of XYZ coords
  DOUBLE PRECISION, INTENT(OUT) :: RMAT(3,3) ! The rotation matrix
  !
  INTEGER :: I,J,K,IPIVOT(3),INFO
  DOUBLE PRECISION, PARAMETER :: EPS=1.0D-10 ! for planar systems
  DOUBLE PRECISION :: AMAT(3,3), CMAT(3,3), WORK(3)
  !
  AMAT(1:3,1:3)=0.0D0
  CMAT(1:3,1:3)=0.0D0
  DO I=1,N
     DO J=1,3
        DO K=1,J
           AMAT(J,K) = AMAT(J,K) + X0(J,I)*X0(K,I) + EPS
           CMAT(J,K) = CMAT(J,K) + X(J,I)*X0(K,I)
        ENDDO
        DO K=J+1,3
           CMAT(J,K) = CMAT(J,K) + X(J,I)*X0(K,I)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Symmetrise AMAT
  AMAT(1,2)=AMAT(2,1); AMAT(1,3)=AMAT(3,1); AMAT(2,3)=AMAT(3,2)
  !
  ! Invert AMAT using routines from LAPACK 
  CALL DGETRF(3,3,AMAT,3,IPIVOT,INFO) ! AMAT is modified!!!   
  CALL DGETRI(3,AMAT,3,IPIVOT,WORK,3,INFO)
  !
  IF (INFO.NE.0) THEN
     PRINT '(A,I6)','ERROR - INFO after DGETRI in FIND_RMAT=',INFO
     PRINT '(A)','X:'
     PRINT '(3G20.10)',X
     PRINT '(A)','X0:'
     PRINT '(3G20.10)',X0
     PRINT '(A)','AMAT:'
     PRINT '(3G20.10)',AMAT
     PRINT '(A)','CMAT:'
     PRINT '(3G20.10)',CMAT
     STOP
  ENDIF
  !                           -1
  ! Compute RMAT = CMAT * AMAT   [ least squares fitting ]
  CALL PRODMM(3,3,3,CMAT,AMAT,RMAT)
  !
  RETURN
  !
END SUBROUTINE FIND_RMAT
!
SUBROUTINE VCROSSPROD(A,B,C)
  !
  ! Returns C = A x B (the cross-product of two vectors)
  !
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: A(3), B(3)
  DOUBLE PRECISION, INTENT(OUT) :: C(3)
  !
  C(3)=A(1)*B(2)-B(1)*A(2)
  C(2)=-A(1)*B(3)+A(3)*B(1)
  C(1)=A(2)*B(3)-A(3)*B(2)
  !
  RETURN
  !
END SUBROUTINE VCROSSPROD
!
SUBROUTINE COMPARE_SYMMATS(A,B,TOL,SAME)
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: A(3,3), B(3,3), TOL
  LOGICAL, INTENT(OUT) :: SAME
  !
  INTEGER :: I,J
  !
  SAME = .TRUE.
  DO I=1,3
     DO J=1,3
        IF(ABS(A(I,J)-B(I,J)) > TOL) THEN
           SAME=.FALSE.
           RETURN
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE COMPARE_SYMMATS
