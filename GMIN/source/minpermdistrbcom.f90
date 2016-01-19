!
!     Follows the same procedure as in minpermdist.f90 to align a system of rigid bodies based
!     on the centres of mass only.
!     ----------------------------------------------------------------------------------------------
! jdf43>        Added to GMIN 30/01/12
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE MINPERMDISTRBCOM(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT)
!     DISTANCE returns the squared distance

      USE COMMONS, ONLY : NATOMS, NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, EFIELDT, MYUNIT, GEOMDIFFTOL 

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXIMUMTRIES = 100
      INTEGER            :: NPERM, PATOMS, NTRIES, NSIZE, JMAX, LOCMAX(1), J1, J2, J3, INFO
      INTEGER            :: INVERT, NORBIT1, NORBIT2, PERM(NATOMS), NCHOOSE2, NDUMMY, LPERM(NATOMS), NCHOOSE1
      INTEGER            :: NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)
      DOUBLE PRECISION   :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DISTWP, DIST2, TEMPA(9*NATOMS) 
      DOUBLE PRECISION   :: DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS), DUMMYWP(3*NATOMS)
!      DOUBLE PRECISION   :: XA(3*NATOMS*NRBSITES/2),  XB(3*NATOMS*NRBSITES/2), XBS(3*NATOMS*NRBSITES/2)
      DOUBLE PRECISION   :: XA(3*NATOMS),  XB(3*NATOMS), XBS(3*NATOMS)
      DOUBLE PRECISION   :: XTMP(3*NATOMS)
      DOUBLE PRECISION   :: RMAT(3,3), RMATI(3,3), ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
      DOUBLE PRECISION   :: ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTBINV(3,3), RMATCUMUL(3,3), LMAT(3,3)
      DOUBLE PRECISION   :: ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), RMATWP(3,3)
      DOUBLE PRECISION   :: CMX, CMY, CMZ, CMBX, CMBY, CMBZ
      DOUBLE PRECISION   :: PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, XDUMMY, BOXLX, BOXLY, BOXLZ, WORSTRAD
      DOUBLE PRECISION   :: Q(4), Q1(4), Q2(4), AMAT(4,4), BMAT(4,4), DIAG(4), P(3)
      DOUBLE PRECISION   :: ST, THETA, THETAH, FCT, DUMMYC(3*NATOMS), DUMMYD(3*NATOMS)
      LOGICAL            :: DEBUG, BULKT
      DOUBLE PRECISION   :: BESTA(3*NATOMS), RBDISTANCE, PVEC(3), RTEMP1(3,3), RTEMP2(3,3), SITEDIST
      DOUBLE PRECISION   :: QCUMUL(4), QBEST(4), QA(4), QB(4), QBINV(4), QTMP(4), QI(4)
      DOUBLE PRECISION   :: COORDSAS(3*NATOMS), COORDSBS(3*NATOMS), T(3*NATOMS)

      COORDSAS(1:3*NATOMS) = COORDSA(1:3*NATOMS) ! to trace error, see at the end
      COORDSBS(1:3*NATOMS) = COORDSB(1:3*NATOMS) ! to trace error, see at the end

      CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
      DO J1 = 1, NATOMS
         J2 = 3*J1
         CMBX = CMBX + COORDSB(J2-2)
         CMBY = CMBY + COORDSB(J2-1)
         CMBZ = CMBZ + COORDSB(J2)
       ENDDO
      CMBX = CMBX/NATOMS; CMBY = CMBY/NATOMS; CMBZ = CMBZ/NATOMS
!
!     Bring COORDSB into standard orientation with respect to the site position
!     The standard orientation needs to be done for the sites if we are going to identify
!     permutation-inversion isomers with respect to the sites metric!
!
!     ----------------------------------------------------------------------------------------------
!     !
!     CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))','before ENERGYB =',ENERGY,' RMS=',RMS
!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))',' after ENERGYB =',ENERGY,' RMS=',RMS
!     !
!     ---------------------------------------------------------------------------------------------- 
!
      INVERT = 1

      DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)
      DUMMYC(1:3*NATOMS) = DUMMYB(1:3*NATOMS)
      CALL ORIENTA(DUMMYC,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,QB,DEBUG)
      CALL QROTMAT(QB,ROTB)
      DUMMYB(1:3*NATOMS) = DUMMY(1:3*NATOMS)

      DBEST    = 1.0D100
60    NCHOOSE1 = 0
65    NCHOOSE1 = NCHOOSE1+1
40    NCHOOSE2 = 0
30    NCHOOSE2 = NCHOOSE2+1

      DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)

      DO J1 = 1, NATOMS
         ALLPERM(J1) = J1
      ENDDO
 
!     The optimal alignment returned by minperdist is a local minimum, but may not
!     be the global minimum. Calling RBORIENT first should put permutational isomers
!     into a standard alignment and spot the global minimum zedro distance in one
!     go. However, we also need to cycle over equivalent atoms in orbits using NCHOOSE2.
!
!     Problems can occur if we don't use all the atoms specified by NORBIT1 and NORBIT2
!     because of the numerical cutoffs employed in MYORIENT. We could miss the
!     right orientation! 
!
!     If we use RBORIENT to produce particular orientations then we end up aligning 
!     COORDSA not with COORDSB but with the standard orientation of COORDSB in DUMMYB.
!     We now deal with this by tracking the complete transformation, including the
!     contribution of MYORIENT using ROTB and ROTINVB.
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO J1 = 1, NATOMS
         J2 = 3*J1
         CMX = CMX + DUMMYA(J2-2)
         CMY = CMY + DUMMYA(J2-1)
         CMZ = CMZ + DUMMYA(J2)
      ENDDO

      CMX = CMX/NATOMS; CMY = CMY/NATOMS; CMZ = CMZ/NATOMS
      DO J1 = 1, NATOMS
         J2 = 3*J1
         DUMMYA(J2-2) = DUMMYA(J2-2) - CMX
         DUMMYA(J2-1) = DUMMYA(J2-1) - CMY
         DUMMYA(J2)   = DUMMYA(J2)   - CMZ
      ENDDO

      IF (EFIELDT .AND. INVERT ==-1) THEN
         RMATI(:,:) = 0.D0
         RMATI(1,1) = 1.D0; RMATI(2,2) =-1.D0; RMATI(3,3) = 1.D0
         DO J1 = 1, NATOMS
            J2 = 3*J1
            DUMMYC(J2-2:J2) = MATMUL(RMATI,DUMMYA(J2-2:J2))
         ENDDO
      ELSE
         DUMMYC(1:3*NATOMS) = INVERT*DUMMYA(1:3*NATOMS)
      ENDIF

      CALL ORIENTA(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,QA,DEBUG)
      CALL QROTMAT(QA,ROTA)
      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)

!      PRINT *, 'DUMMYA'
!      PRINT *, DUMMYA
!      PRINT *, 'DUMMYB'
!      PRINT *, DUMMYB
!      STOP
!     ----------------------------------------------------------------------------------------------
!     !
!      CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYBD =',ENERGY,' RMS=',RMS
!      CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYAD =',ENERGY,' RMS=',RMS
!
!      DISTANCE = 0.D0
!      DO J1 = 1, 3*NRB
!         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
!      ENDDO
!     !
!     ----------------------------------------------------------------------------------------------
      
      DISTANCE = 0.D0
      DO J1 = 1, 3*NATOMS
         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
      ENDDO

!     ----------------------------------------------------------------------------------------------
!     !
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'A sites before rotation'
!      WRITE(975,'(A,3F20.10)') ('LA ',XA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'B sites before rotation'
!      WRITE(975,'(A,3F20.10)') ('LA ',XB(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!     !
!     ----------------------------------------------------------------------------------------------

      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)')' minpermdistrbcom> after initial call to RBSITESORIENT distance=',SQRT(DISTANCE)
!
!     Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!     but the coordinates in DUMMYA do. DISTANCE is the distance in this case.
!     We return to label 10 after every round of permutational/orientational alignment
!     unless we have converged to the identity permutation.
!
!     Atoms are not allowed to appear in more than one group.
!     The maximum number of pair exchanges associated with a group is two.
!
      NTRIES = 0
!
!     QCUMUL is a quaternion containing the information of accumulated rotation that relates the  
!     original DUMMYA obtained from COORDSA to the final one. RMATCUMUL is the corresponding
!     rotation matrix. Initialize QCUMUL, so as RMATCUMUL.
!
      QCUMUL(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)
      RMATCUMUL(:,:) = 0.D0; RMATCUMUL(1,1) = 1.D0; RMATCUMUL(2,2) = 1.D0; RMATCUMUL(3,3) = 1.D0
10    CONTINUE

      NTRIES = NTRIES + 1
      NDUMMY = 1

      DO J1 = 1, NATOMS
         NEWPERM(J1) = J1
      ENDDO

!     ALLPERM saves the permutation from the previous cycle.
!     NEWPERM contains the permutation for this cycle, relative to the identity.
!     SAVEPERM is temporary storage for NEWPERM.
!     NEWPERM must be applied to ALLPERM after the loop over NPERMGROUP and corresponding swaps.

      DO J1 = 1, NPERMGROUP

         PATOMS = NPERMSIZE(J1)
          IF (PATOMS.GT.NATOMS) THEN
             WRITE(MYUNIT,'(A,I6,A,I6,A,I6)')' minpermdistrbcom> ERROR *** number of permutable sites in group ',J1,' is ',PATOMS, &
  &                                   ' which exceeds the number of atoms ',NATOMS
             STOP
         ENDIF

         DO J2 = 1, PATOMS
            PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
            PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
            PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
            PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
            PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
            PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
         ENDDO
!
!     All permutations within this group of size NPERMSIZE(J1) are now tried.
!     Note that we are just using a metric based on the rigid body centre of mass coordinates here!
!

         CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)

         SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
         DO J2=1,PATOMS
            SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1))
         ENDDO
!
! Update permutation of associated atoms, if any.
! We must do this as we go along, because these atoms could move in more than
! one permutational group now.
!
         IF (NSETS(J1) > 0) THEN
            DO J2 = 1, PATOMS
               DO J3 = 1, NSETS(J1)
                  SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J3)
               ENDDO
            ENDDO
         ENDIF
         NDUMMY = NDUMMY + NPERMSIZE(J1)
         NEWPERM(1:NATOMS) = SAVEPERM(1:NATOMS)

      ENDDO
!
!     Due to possible swaps above, we need to update the overall permutation here.
!
!     Update the overall permutation here.
!
      DO J1=1,NATOMS
!        SAVEPERM(ALLPERM(J1))=ALLPERM(NEWPERM(J1)) BUG FIX 12.9.11 DJW
         SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
      ENDDO

      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)
      DUMMY(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
      NPERM    = 0
      DISTANCE = 0.0D0
!
!     Now permute the rotational coordinates to correspond to the centre of mass coordinate permutation.
!
!      DO J1 = (NATOMS/2)+1, NATOMS
!         IF (PERM(J1).NE.J1) THEN
!            PRINT '(A,2I8)',' minpermdistrbcom> ERROR - J1,PERM(J1) should be equal', J1,PERM(J1)
!         ENDIF
!         PERM(J1) = PERM(J1-(NATOMS/2)) + (NATOMS/2)
!      ENDDO

      DO J3 = 1, NATOMS

         DUMMYA(3*(J3-1)+1) = DUMMY(3*(NEWPERM(J3)-1)+1)
         DUMMYA(3*(J3-1)+2) = DUMMY(3*(NEWPERM(J3)-1)+2)
         DUMMYA(3*(J3-1)+3) = DUMMY(3*(NEWPERM(J3)-1)+3)

         IF (J3.NE.NEWPERM(J3)) THEN

            IF (DEBUG) WRITE(MYUNIT,'(A,I5,A,I5)') ' minpermdistrbcom> move position ',NEWPERM(J3),' to ',J3
            NPERM = NPERM + 1

         ENDIF

      ENDDO

      DO J1 = 1, 3*NATOMS
         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
      ENDDO
!
!     Further alignment. Coordinates in DUMMYA are reset by RBMINDIST (second argument).
!     Must allow at least one call to RBMINDIST in case the standard orientation result is terrible
!     but gives zero permutations!
!     We try internal symmetry operations first for each rigid body, then minimise the
!     distance further (if possible) .
!  
      IF ((NPERM.NE.0) .OR. (NTRIES.EQ.1)) THEN 
!
!     Now if RBSYMT is .TRUE. we should minimise the distance metric for all the rigid bodies
!     by considering all the allowed symmetry operations for each rigid body in turn.
!     Do this by altering the orientational coordinates in DUMMYA for each rigid body in turn,
!     after saving them, calculate the new distance using the all sites metric, and accept the
!     new orientation if the distance is lower.
!     Could generalise to more than one sort of rigid body as well.
! 
!     This call aligns the overall orientation with respect to the sites metric.
!     Internal symmetry operations of the rigid bodies are not considered.
!
            CALL MINDISTA(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)

            CALL QROTMAT(Q2,RMAT)

            DISTANCE  = DISTANCE*DISTANCE

!     accumulate rotation

         CALL QROTQ(Q2,QCUMUL)

         RMATCUMUL(:,:) = MATMUL(RMAT,RMATCUMUL)  ! also in the form of the rotationa matrix

         IF (NTRIES .LT. MAXIMUMTRIES) THEN
            GOTO 10
         ELSE ! prevent infinite loop
            IF (INVERT == -1) THEN
               IF(DEBUG) WRITE(MYUNIT,'(A)')' minpermdistrbcom> WARNING - number of tries exceeded, giving up'
            ELSE
               WRITE(MYUNIT,'(A)')' minpermdistrbcom> WARNING - number of tries exceeded, giving up'
            ENDIF
         ENDIF

      ENDIF

      IF (DISTANCE .LT. DBEST) THEN

         DBEST                = DISTANCE
         XBEST(1:3*NATOMS)    = DUMMYA(1:3*NATOMS)
         QTMP(1:4)  = QA(1:4)
         CALL QROTQ(QCUMUL,QTMP)
         QBEST(1:4) = QTMP(1:4)
         QBINV(1:4) = (/QB(1), -QB(2:4)/)
         RMATBEST(:,:) = MATMUL(RMATCUMUL,ROTA)
         IF (INVERT == -1) THEN
            IF (EFIELDT) THEN
               RMATBEST(:,:) = MATMUL(RMATBEST,RMATI)
            ELSE
               RMATBEST(:,:) =-RMATBEST(:,:)
            ENDIF
         ENDIF

      ENDIF
!
! If GEOMDIFFTOL is set too sloppy we could miss the best solution by exiting via the line
! below without trying other orbits. Turn off this escape?!
! The turn off seems to be a bug for non-RBPERM runs. LJ38 tests fail with all compilers, producing
! a "Distance is zero: this should not happen" error!!! DJW
! Put the escape back in for now.
!
      IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
      IF (NCHOOSE2.LT.NORBIT2) GOTO 30
      IF (NCHOOSE1.LT.NORBIT1) GOTO 65
      IF (EFIELDT) GOTO 50
!      GOTO 50  !!!! this is a workaround to prevent problems with making the inverted structure. DJW

      IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
!     don't try inversion for bulk or charmm or amber or frozen atoms
!
         IF (DEBUG) WRITE(MYUNIT,'(A)')' minpermdistrbcom> inverting geometry for comparison with target'
         INVERT=-1
         GOTO 60

      ENDIF

!
50    DISTANCE = DBEST       ! squared distance
!
!     XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!     Rotate XBEST by ROTINVBBEST to put in best correspondence with COORDSB, undoing the reorientation to  
!     DUMMYB from RBORIENT. 
!     We should get the same result for ROTINVBBEST * RMATBEST * (COORDSA-CMA), 
!     where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
!     (aside from a possible permutation of the atom ordering)
!
!      CALL QROTMAT(QBINV, ROTBINV)

      ROTBINV = TRANSPOSE(ROTB)
 
      DO J1 = 1, NATOMS
         J2 = 3*J1
         XBEST(J2-2:J2) = MATMUL(ROTBINV,XBEST(J2-2:J2))
      ENDDO

      DO J1 = 1, NATOMS
         XBEST(3*(J1-1)+1) = XBEST(3*(J1-1)+1) + CMBX
         XBEST(3*(J1-1)+2) = XBEST(3*(J1-1)+2) + CMBY
         XBEST(3*(J1-1)+3) = XBEST(3*(J1-1)+3) + CMBZ
      ENDDO

      XDUMMY = 0.D0
      DO J1 = 1, 3*NATOMS
         XDUMMY = XDUMMY + (XBEST(J1) - COORDSB(J1))**2
      ENDDO

      IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
          WRITE(MYUNIT,'(2(A,F20.10))')'minpermdistrbcom> ERROR *** distance between transformed XBEST and COORDSB=',  &
     &    SQRT(XDUMMY),  ' should be ', SQRT(DISTANCE)
          IF(DEBUG) THEN
             WRITE(MYUNIT)'COORDSBCOM'
             DO J1 = 1, NATOMS
                J2 = 3*J1
                WRITE(MYUNIT)COORDSBS(J2-2), COORDSBS(J2-1), COORDSBS(J2)
             ENDDO 
             WRITE(MYUNIT)'COORDSACOM'
!             PRINT *, 'XBEST'
             DO J1 = 1, NATOMS
                J2 = 3*J1
                WRITE(MYUNIT)COORDSAS(J2-2), COORDSAS(J2-1), COORDSAS(J2)
!                PRINT *, XBEST(J2-2), XBEST(J2-1), XBEST(J2)
             ENDDO
          ENDIF 

          STOP

      ENDIF

      CALL QROTQ(QBINV,QBEST)

      RMATBEST(:,:) = MATMUL(ROTBINV,RMATBEST)

      COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!      CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))',' finally A energy=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

!      PRINT *, RMATBEST
!      STOP

      END SUBROUTINE MINPERMDISTRBCOM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ORIENTA(X, T1S, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)

!     This subroutine puts the configuration, X, of an atomic  system into a standard alignment, T1.
!
      USE COMMONS, ONLY: EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
      DOUBLE PRECISION :: X(3*NATOMS), T1(3*NATOMS)
      DOUBLE PRECISION :: XS(3*NATOMS), T1S(3*NATOMS), T2S(3*NATOMS), DIST(NATOMS)
      DOUBLE PRECISION :: AX(3), P(3), Q2(4), ROTM(3,3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, COST, SINT, COSTH, SINTH, ST, FCT
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, DUMMY, PROJ, DMAX2, CUTOFF1, DTEMP
      LOGICAL          :: DEBUG

!      EFIELDT = .FALSE.
      CUTOFF1 = 1.D-03
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS
         J = 3*I
         CMX = CMX + X(J-2)
         CMY = CMY + X(J-1)
         CMZ = CMZ + X(J)
      ENDDO

      CMX = CMX/NATOMS; CMY = CMY/NATOMS; CMZ = CMZ/NATOMS
      DO I = 1, NATOMS
         J = 3*I
         XS(J-2) = X(J-2) - CMX
         XS(J-1) = X(J-1) - CMY
         XS(J)   = X(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN
         T1S(:) = XS(:)
         GOTO 100
      ENDIF

!
!     Find the atom which is at the largest distance from the centre of mass
!
      DO J1 = 1, NATOMS

         J = 3*J1
         DIST(J1) = SQRT(XS(J-2)**2 + XS(J-1)**2 + XS(J)**2)

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            NORBIT1 = NORBIT1 + 1
            IF (NORBIT1 == NCHOOSE1) JMAX1 = J1

         ELSE IF (DIST(J1) > DMAX) THEN

            DMAX = DIST(J1)
            NORBIT1 = 1
            JMAX1 = J1

         ENDIF

      ENDDO

!     For tagged atoms, the choice of the first atom matters if it belongs to an orbit of size > 1.
!
      IF ((ABS(XS(3*JMAX1-2)) < 1.D-08) .AND. (ABS(XS(3*JMAX1-1)) < 1.D-08)) THEN

!
!     atom JMAX1 is already on the z axis!
!
         IF (XS(3*(JMAX1-1)+3) > 0.D0) THEN

            T1S(1:3*NATOMS) = XS(1:3*NATOMS)
            Q2(1:4)   = (/1.D0, 0.D0, 0.D0, 0.D0/)   ! Identity operation       
                  
         ELSE  ! rotate about the x-axis by \pi, DO NOT INVERT!

            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)      ! the corresponding quaternion

            DO J1 = 1, NATOMS
               J       = 3*J1
               T1S(J-2) = XS(J-2)
               T1S(J-1) =-XS(J-1)
               T1S(J)   =-XS(J)
            ENDDO

         ENDIF

      ELSE

!
!     Rotate all atoms through an angle THETA about the axis P(:), so as to rotate the rbsite JMAX1 onto the z axis 
!
!       THETA = DACOS(XS(3*JMAX1)/DMAX)
!     For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!

        THETA = DACOS(XS(3*JMAX1)/DIST(JMAX1))

!     The axis is on the xy-plane and perpendicular to the projection of the translation vector of the rb JMAX1
!     on the xy-plane. The axis is so chosen that a right-handed rotation in the right-handed coordinate system
!     will result in the required transformation.

         AX(3) = 0.D0
         FCT   = DSQRT(XS(3*JMAX1-2)**2 + XS(3*JMAX1-1)**2)
         AX(1) = THETA*XS(3*JMAX1-1)/FCT
         AX(2) =-THETA*XS(3*JMAX1-2)/FCT

         CALL ROTMAT(AX(:), ROTM(:,:))
         THETAH  = 0.5D0*THETA
         AX(1:3) = AX(1:3)/SQRT(DOT_PRODUCT(AX(1:3),AX(1:3)))
         Q2(1:4) = (/COS(THETAH), SIN(THETAH)*AX(1:3)/)          ! the corresponding quaternion 
          
         DO J1 = 1, NATOMS
            J2 = 3*J1
            T1S(J2-2:J2) = MATMUL(ROTM(:,:),XS(J2-2:J2))
         ENDDO   

      ENDIF

!
!     Now find the atom with the largest distance from the z-axis 
!
100   DMAX = -1.0D0

      DO J1 = 1, NATOMS
         J2 = 3*J1 
         DIST(J1) = SQRT(T1S(J2-2)**2 + T1S(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ is the sum of the x components. Use T2S as a dummy in order not to 
!     change T1 until we have decided which atom to put in the xz plane.
!
      DO J1 = 1, NATOMS

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2S(1:3*NATOMS) = T1S(1:3*NATOMS)

            CALL ROTATMXZ(NATOMS, J1, T2S, PROJ, DIST, Q2, .FALSE.)

            IF (ABS(PROJ - DMAX2) < CUTOFF1) THEN
               NORBIT2 = NORBIT2+1
               IF (NORBIT2 == NCHOOSE2) THEN
                  JMAX2 = J1
                  DTEMP = PROJ
               ENDIF
            ELSE IF (PROJ > DMAX2) THEN
               NORBIT2 = 1
               DMAX2   = PROJ
               DTEMP   = PROJ
               JMAX2   = J1
            ENDIF

         ENDIF

      ENDDO

!
!     and now rotate it into the xz plane.
!

      CALL ROTATMXZ(NATOMS, JMAX2, T1S, DTEMP, DIST, Q2, .TRUE.)
!      EFIELDT = .TRUE.

      END SUBROUTINE ORIENTA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTATMXZ(NATOMS, JDO, T1S, PROJ, DIST, Q1, ROTT)

      USE COMMONS, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, JDO, I, J, J2
      DOUBLE PRECISION :: T1S(3*NATOMS), PROJ, DIST(NATOMS), THETA, THETAH, COST, SINT, ST
      DOUBLE PRECISION :: COSTH, SINTH, FCT, P(3), Q2(4), Q1(4), RM(3,3), ROTM(3,3) 
      LOGICAL          :: ROTT

      J2     = 3*JDO

      IF (ABS(T1S(J2-1)) < 1.0D-8) THEN            ! already on the xz plane

         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/) 

         IF (T1S(J2-2) < 0.D0) THEN                ! rotate about the z axis by \pi, DO NOT INVERT!!

            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/) 

            DO I = 1, NATOMS
               J = 3*I
               T1S(J-2) =-T1S(J-2)
               T1S(J-1) =-T1S(J-1)
            ENDDO

         ENDIF

      ELSE

         COST   = T1S(J2-2)/DIST(JDO)       
         SINT   = T1S(J2-1)/DIST(JDO)
         THETA  =-ATAN2(SINT,COST)              ! the negative sign appears as the rotation is reversed  
         THETAH = 0.5D0*THETA
         Q2(1:4) = (/COS(THETAH), 0.D0, 0.D0, SIN(THETAH)/)     ! the quaternion corresponding to R_{z}(-\theta)
         RM(1:3,1:3) = 0.D0                           ! the rotation matrix for R_{z}(-\theta)
         RM(1,1) = COST
         RM(2,2) = COST
         RM(3,3) = 1.D0
         RM(1,2) = SINT
         RM(2,1) =-SINT

         CALL QROTMAT(Q2,RM)

         DO I = 1, NATOMS
            IF (DIST(I) /= 0.0D0) THEN
               J = 3*I
               T1S(J-2:J) = MATMUL(RM(:,:),T1S(J-2:J))
            ENDIF  
         ENDDO

      ENDIF

      IF (ROTT) THEN

         IF (EFIELDT) THEN
            Q1(1:4) = Q2(1:4)
         ELSE
            CALL QROTQ(Q2,Q1)
         ENDIF

      ELSE

         PROJ = 0.D0

         DO I = 1, NATOMS 

            J = 3*I
            IF (T1S(J) > 1.0D-2) PROJ = PROJ + T1S(J-2)

         ENDDO

      ENDIF

      END SUBROUTINE ROTATMXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTATM(T, X, Q2, NATOMS)
!     takes the set of coordinates T for NATOMS number of atoms and returns X after rotation via the 
!     quaternion Q2 about the origin 

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: T(3*NATOMS), X(3*NATOMS), T1(1:3), Q2(4), RM(3,3) 
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = CMX/NATOMS; CMY = CMY/NATOMS; CMZ = CMZ/NATOMS
      DO I = 1, NATOMS
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
         X(J)   = T(J)   - CMZ
      ENDDO

!     extract the rotation matrix corresponding to rotation via Q

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS

         J        = 3*I
         T1(1:3)  = MATMUL(RM(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)  

      ENDDO

      END SUBROUTINE ROTATM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE MINDISTA(RA,RB,NATOMS,DIST,Q2,DEBUG)

!     Returns DIST as the actual distance, rather than the squared distance

      USE COMMONS, ONLY: EFIELDT, MYUNIT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3), THETA, FCT
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG

      NSIZE = NATOMS
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))
      XA(1:3*NSIZE) = RA(1:3*NSIZE); XB(1:3*NSIZE) = RB(1:3*NSIZE)

!     MOVE CENTRE OF COORDINATES OF XA AND XB TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXA = CMXA + XA(J2+1)
         CMYA = CMYA + XA(J2+2)
         CMZA = CMZA + XA(J2+3)
      ENDDO
      CMXA = CMXA/NSIZE; CMYA = CMYA/NSIZE; CMZA = CMZA/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XA(J2+1) = XA(J2+1) - CMXA
         XA(J2+2) = XA(J2+2) - CMYA
         XA(J2+3) = XA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXB = CMXB + XB(J2+1)
         CMYB = CMYB + XB(J2+2)
         CMZB = CMZB + XB(J2+3)
      ENDDO
      CMXB = CMXB/NSIZE; CMYB = CMYB/NSIZE; CMZB = CMZB/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XB(J2+1) = XB(J2+1) - CMXB
         XB(J2+2) = XB(J2+2) - CMYB
         XB(J2+3) = XB(J2+3) - CMZB
      ENDDO

      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)

         IF (EFIELDT) THEN
            QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
            QMAT(1,2) = QMAT(1,2) - XP*YM + XM*YP
            QMAT(2,2) = QMAT(2,2) + XP**2 + YP**2 + ZM**2
         ELSE
            QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
            QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP
            QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
            QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
            QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
            QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP
            QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
            QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
            QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP
            QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
         ENDIF
      ENDDO

      IF (EFIELDT) THEN

!     QMAT IS SYMMETRIC; QMAT(2,1) = QMAT(1,2)

         MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))
         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
         Q2(2) = 0.D0
         Q2(3) = 0.D0
         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))

         IF (MINV < 0.0D0) THEN
            IF (ABS(MINV)< 1.0D-6) THEN
               MINV = 0.0D0
            ELSE
               WRITE(MYUNIT,'(A,G20.10,A)')'newmindist> WARNING MINV is ',MINV,' change to absolute value'
               MINV = -MINV
            ENDIF
         ENDIF
       
      ELSE

        QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4)
        QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)
        CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

        IF (INFO /= 0) WRITE(MYUNIT,'(A,I6,A)')'mindista> WARNING - INFO=',INFO,' in DSYEV'

        MINV = 1.0D100
        DO J1 = 1,4
           IF (DIAG(J1).LT.MINV) THEN
              JMIN = J1
              MINV = DIAG(J1)
           ENDIF
        ENDDO
        IF (MINV < 0.0D0) THEN
           IF (ABS(MINV)< 1.0D-6) THEN
              MINV = 0.0D0
           ELSE
              WRITE(MYUNIT,'(A,G20.10,A)')'mindista> WARNING MINV is ',MINV,' change to absolute value'
              MINV = -MINV
           ENDIF
        ENDIF

        Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      ENDIF

      DIST = SQRT(MINV)

      DO J1 = 1, NATOMS
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL NEWROTGEOMA(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE MINDISTA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE NEWROTGEOMA(NATOMS,COORDS,Q2,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS

         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ

      ENDDO

      END SUBROUTINE NEWROTGEOMA
