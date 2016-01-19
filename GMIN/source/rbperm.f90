!
!     COORDSA becomes the optimal alignment of the optimal permutation(-inversion) isomer, but 
!     without the permutations. DISTANCE is the residual square distance for the best alignment with
!     respect to permutation(-inversion)s as well as orientation and centre of mass.
!
!     RBSITESORIENT is called first for both COORDSA and COORDSB to put them into a standard
!     orientation in DUMMYA and DUMMYB (which both have the centre of coordinates at the origin)
!     with respect to the sites via rotations corresponding to the quaternions QA and QB, 
!     respectively. The objective is to identify permutation-inversion isomers without fail.
!     However, we have to cycle over all equivalent rigid bodies in two particular orbits for DUMMYA
!     to achieve this.
!     We iterate permutations and rbmindist minimisations up to a maximum number or until no more
!     permutations are required for each instance of DUMMYA aligned according to NCHOOSE1 and 
!     NCHOOSE2 by RBSITESORIENT. The cumulative rotation that takes the initial DUMMYA to the one
!     that aligns best with DUMMYB is saved in the quaternion QCUMUL.

!     Then, if we've not going with BULK, we try again for the inverted
!     version of COORDSA. The transformation corresponding to the minimum distance
!     is saved whenever it is improved - the best alignment including permutations
!     is saved in XBEST, and the last step is to rotate this back to coincide best
!     with COORDSB (rather than DUMMYB) using QBINV which is obtained from QB. This gives suitable
!     fixed end points for DNEB.

!     Finally, we transform COORDSA to be in optimal alignment, but without the
!     permutations in XBEST. The overall transformation is
!     COORDSA -> +/- QBINV QCUMUL QA (COORDSA - CMA) 
!
!     The correspondence between COORDSA and DUMMYA after DUMMYA has been aligned by rbmindist is
!     +/- RMATCUMUL ROTA (COORDSA - CMA) = permutation(DUMMYA)
!     where +/- is given by the value of INVERT which can assume +1 or -1.
!     The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
!     centre of coordinates of COORDSA will be the same as for COORDSB.
!
!     ----------------------------------------------------------------------------------------------
! jdf43>        Added to GMIN 30/01/12
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE RBMINPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,SITESB,SITESA)

!     returns DISTANCE as the actual distance, rather than the squared distance

      USE COMMONS, ONLY : NATOMS, NRBSITES, NTSITES, NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, &
     &                    GEOMDIFFTOL, BESTPERM, EFIELDT, RBOPS, NRBGROUP, RBSYMT, DBPTDT, MYUNIT

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXIMUMTRIES = 100
      INTEGER            :: NPERM, PATOMS, NTRIES, NSIZE, JMAX, LOCMAX(1), J1, J2, J3, INFO, NRB
      INTEGER            :: INVERT, NORBIT1, NORBIT2, PERM(NATOMS), NCHOOSE2, NDUMMY, LPERM(NATOMS), NCHOOSE1
      INTEGER            :: NEWPERM(NATOMS), ALLPERM(NATOMS)
      DOUBLE PRECISION   :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DISTWP, DIST2, TEMPA(9*NATOMS) 
      DOUBLE PRECISION   :: DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS), DUMMYWP(3*NATOMS)
      DOUBLE PRECISION   :: XA(3*NTSITES),  XB(3*NTSITES), XBS(3*NTSITES), XTMP(3*NATOMS)
      DOUBLE PRECISION   :: SITESA(3*NTSITES),  SITESB(3*NTSITES)
      DOUBLE PRECISION   :: RMAT(3,3), RMATI(3,3), ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
      DOUBLE PRECISION   :: ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), RMATCUMUL(3,3), LMAT(3,3)
      DOUBLE PRECISION   :: ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), RMATWP(3,3)
      DOUBLE PRECISION   :: CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ
      DOUBLE PRECISION   :: PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, XDUMMY, BOXLX, BOXLY, BOXLZ, WORSTRAD
      DOUBLE PRECISION   :: Q(4), Q1(4), Q2(4), AMAT(4,4), BMAT(4,4), DIAG(4), P(3)
      DOUBLE PRECISION   :: ST, THETA, THETAH, FCT, DUMMYC(3*NATOMS), DUMMYD(3*NATOMS)
      LOGICAL            :: DEBUG, BULKT
      DOUBLE PRECISION   :: BESTA(3*NATOMS), RBDISTANCE, PVEC(3), RTEMP1(3,3), RTEMP2(3,3), SITEDIST
      DOUBLE PRECISION   :: QCUMUL(4), QBEST(4), QA(4), QB(4), QBINV(4), QTMP(4), QI(4)
      DOUBLE PRECISION   :: COORDSAS(3*NATOMS), COORDSBS(3*NATOMS), COORDSCOMA(3*NATOMS/2), COORDSCOMB(3*NATOMS/2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      COORDSCOMA(:) = COORDSA(1:3*NATOMS/2)
!      COORDSCOMB(:) = COORDSB(1:3*NATOMS/2)
!
!      NATOMS = NATOMS/2
!
!      CALL  MINPERMDISTRBCOM(COORDSCOMB,COORDSCOMA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT)
!
!      NATOMS = 2*NATOMS
!
!      IF (SQRT(DISTANCE) <= GEOMDIFFTOL/1.D2) THEN
!         DISTANCE = SQRT(DISTANCE)
!         IF (DEBUG) WRITE(MYUNIT,'(A)')' rbpermdist> minpermdistrbcom suggests identical'
!!         WRITE(MYUNIT,'(A)')' rbpermdist> minpermdistrbcom suggests identical'
!         CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
!         CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
!         DO J1 = 1, NATOMS/2
!            CMAX = CMAX + COORDSA(3*(J1-1)+1)
!            CMAY = CMAY + COORDSA(3*(J1-1)+2)
!            CMAZ = CMAZ + COORDSA(3*(J1-1)+3)
!            CMBX = CMBX + COORDSB(3*(J1-1)+1)
!            CMBY = CMBY + COORDSB(3*(J1-1)+2)
!            CMBZ = CMBZ + COORDSB(3*(J1-1)+3)
!         ENDDO
!         CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
!         CMBX = 2*CMBX/NATOMS; CMBY = 2*CMBY/NATOMS; CMBZ = 2*CMBZ/NATOMS
!         DO J1 = 1, NATOMS/2
!            COORDSA(3*(J1-1)+1) = COORDSA(3*(J1-1)+1) - CMAX
!            COORDSA(3*(J1-1)+2) = COORDSA(3*(J1-1)+2) - CMAY
!            COORDSA(3*(J1-1)+3) = COORDSA(3*(J1-1)+3) - CMAZ
!            COORDSB(3*(J1-1)+1) = COORDSB(3*(J1-1)+1) - CMBX
!            COORDSB(3*(J1-1)+2) = COORDSB(3*(J1-1)+2) - CMBY
!            COORDSB(3*(J1-1)+3) = COORDSB(3*(J1-1)+3) - CMBZ
!         ENDDO
!         CALL SITEPOS(COORDSB,SITESB)
!         CALL SITEPOS(COORDSA,SITESA)
!         DO J1 = 1, NTSITES
!            J2 = 3*J1
!            SITESA(J2-2:J2) = MATMUL(RMATBEST,SITESA(J2-2:J2))
!         ENDDO
!         RETURN
!      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

100   CONTINUE

      COORDSAS(1:3*NATOMS) = COORDSA(1:3*NATOMS) ! to trace error, see at the end
      COORDSBS(1:3*NATOMS) = COORDSB(1:3*NATOMS) ! to trace error, see at the end

      NRB   = (NATOMS/2)
      IF (DBPTDT) THEN
         NSIZE = NTSITES
      ELSE
         NSIZE = NRB*NRBSITES
      ENDIF

      CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
      DO J1 = 1, NATOMS/2
         CMBX = CMBX + COORDSB(3*(J1-1)+1)
         CMBY = CMBY + COORDSB(3*(J1-1)+2)
         CMBZ = CMBZ + COORDSB(3*(J1-1)+3)
       ENDDO
      CMBX = 2*CMBX/NATOMS; CMBY = 2*CMBY/NATOMS; CMBZ = 2*CMBZ/NATOMS
!
!     Bring COORDSB into standard orientation with respect to the site position
!     The standard orientation needs to be done for the sites if we are going to identify
!     permutation-inversion isomers with respect to the sites metric!
!
!     DUMMY(1:3*NATOMS)=DUMMYB(1:3*NATOMS)
!     CALL RBSITESORIENT(DUMMY,DUMMYB,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTB,ROTINVB)
!
!     ----------------------------------------------------------------------------------------------
!     !
!     CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     WRITE(MYUNIT,'(2(A,F25.15))')'before ENERGYB =',ENERGY,' RMS=',RMS
!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     WRITE(MYUNIT,'(2(A,F25.15))')' after ENERGYB =',ENERGY,' RMS=',RMS
!     !
!     ---------------------------------------------------------------------------------------------- 
!
!     CALL SITEPOS(DUMMYB,XB)
!     WRITE(975,'(I6)') NRB*NRBSITES
!     WRITE(975,'(A,F20.10)') 'B sites 1'
!     WRITE(975,'(A,3F20.10)') ('LA ',XB(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)

      INVERT = 1

      DBEST    = 1.0D100
60    NCHOOSE1 = 0
65    NCHOOSE1 = NCHOOSE1+1
40    NCHOOSE2 = 0
30    NCHOOSE2 = NCHOOSE2+1

      DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)
      DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)

      DO J1 = 1, NATOMS
         ALLPERM(J1) = J1
         NEWPERM(J1) = J1
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

!      IF (INVERT.NE.1) THEN
!         WRITE(MYUNIT,'(A)')' rbminpermdist> ERROR *** inversion not yet programmed'
!         STOP
!      ENDIF

      DUMMYC(1:3*NATOMS) = DUMMYA(1:3*NATOMS)

      CALL RBSITESORIENT(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,QA,DEBUG)
!      CALL RBSITESORIENTPCA(DUMMYC,DUMMY,NATOMS,QA)

      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      DUMMYC(1:3*NATOMS)=DUMMYB(1:3*NATOMS)

!     CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     WRITE(MYUNIT,'(2(A,F25.15))')'before DUMMYA =',ENERGY,' RMS=',RMS

      CALL RBSITESORIENT(DUMMYC,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,QB,DEBUG)
!      CALL RBSITESORIENTPCA(DUMMYC,DUMMY,NATOMS,QB)

      DUMMYB(1:3*NATOMS)=DUMMY(1:3*NATOMS)

!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     WRITE(MYUNIT,'(2(A,F25.15))')'before DUMMYB =',ENERGY,' RMS=',RMS

!     ----------------------------------------------------------------------------------------------
!     !
!      CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      WRITE(MYUNIT,'(2(A,F25.15))') ' ENERGYBD =',ENERGY,' RMS=',RMS
!      CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      WRITE(MYUNIT,'(2(A,F25.15))') ' ENERGYAD =',ENERGY,' RMS=',RMS
!
!      DISTANCE = 0.D0
!      DO J1 = 1, 3*NRB
!         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
!      ENDDO
!     !
!     ----------------------------------------------------------------------------------------------
      
!      IF (INVERT == -1) ROTA(:,:) = MATMUL(ROTA(:,:),RMATI(:,:))  
!     RMATI was returned on inversion followed by rotation to bring about best superposition on DUMMYB in one go

      IF (INVERT == -1) THEN
         QTMP(1:4) = QI(1:4)
         CALL QROTQ(QA,QTMP)
         QA(1:4) = QTMP(1:4)
      ENDIF  

      CALL SITEPOS(DUMMYA,XA)
      CALL SITEPOS(DUMMYB,XB)

      DISTANCE = 0.D0
      DO J1 = 1, 3*NSIZE
         DISTANCE = DISTANCE + (XA(J1) - XB(J1))**2
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

      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)')' rbminpermdist> after initial call to RBSITESORIENT distance=',SQRT(DISTANCE)
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
!     original DUMMYA obtained from COORDSA to the final one. Initialize QCUMUL.
!
      QCUMUL(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)

10    CONTINUE

      NTRIES = NTRIES + 1
      NDUMMY = 1

      DO J1 = 1, NATOMS
         PERM(J1) = J1
      ENDDO

      DO J1 = 1, NPERMGROUP

         PATOMS = NPERMSIZE(J1)
          IF (PATOMS.GT.NRB) THEN
             WRITE(MYUNIT,'(A,I6,A,I6,A,I6)')' rbminpermdist> ERROR *** number of permutable sites in group ',J1,' is ',PATOMS, &
  &                                   ' which exceeds the number of rigid bodies ',NRB
             STOP
         ENDIF

         DO J2 = 1, PATOMS
            PDUMMYA(3*(J2-1)+1) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
            PDUMMYA(3*(J2-1)+2) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
            PDUMMYA(3*(J2-1)+3) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
            PDUMMYB(3*(J2-1)+1) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
            PDUMMYB(3*(J2-1)+2) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
            PDUMMYB(3*(J2-1)+3) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
         ENDDO
!
!     All permutations within this group of size NPERMSIZE(J1) are now tried.
!     Note that we are just using a metric based on the rigid body centre of mass coordinates here!
!

         CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)

         DO J2 = 1, PATOMS
            PERM(PERMGROUP(NDUMMY+J2-1)) = PERMGROUP(NDUMMY+LPERM(J2)-1)
         ENDDO
!
!  Will this branch ever be executed for rigid bodies?!
!
         IF (NSETS(J1).GT.0) THEN
            WRITE(MYUNIT,'(A)')' rbperm> ERROR *** associated atom swaps not allowed for rigid bodies'
            STOP
         ENDIF

         NDUMMY = NDUMMY + NPERMSIZE(J1)

      ENDDO
!
!     Due to possible swaps above, we need to update the overall permutation here.
!
      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)
      DUMMY(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
      NPERM    = 0
      DISTANCE = 0.0D0
!
!     Now permute the rotational coordinates to correspond to the centre of mass coordinate permutation.
!
      DO J1 = (NATOMS/2)+1, NATOMS
         IF (PERM(J1).NE.J1) THEN
            WRITE(MYUNIT,'(A,2I8)')' rbminpermdist> ERROR - J1,PERM(J1) should be equal', J1,PERM(J1)
         ENDIF
         PERM(J1) = PERM(J1-(NATOMS/2)) + (NATOMS/2)
      ENDDO

      DO J3 = 1, NATOMS

         DUMMYA(3*(J3-1)+1) = DUMMY(3*(PERM(J3)-1)+1)
         DUMMYA(3*(J3-1)+2) = DUMMY(3*(PERM(J3)-1)+2)
         DUMMYA(3*(J3-1)+3) = DUMMY(3*(PERM(J3)-1)+3)

         IF (J3.NE.PERM(J3)) THEN

            IF (DEBUG) WRITE(MYUNIT,'(A,I5,A,I5)') ' rbpermdist> permute rigid bodies ',J3,' and ',PERM(J3)
            NPERM = NPERM + 1
            NEWPERM(J3) = ALLPERM(PERM(J3))

         ENDIF

      ENDDO

!     Site-site distance

      CALL SITEPOS(DUMMYA,XA)      ! DUMMYB remained unaltered, so as XB.

!     WRITE(975,'(I6)') NRB*NRBSITES
!     WRITE(975,'(A,F20.10)') 'A sites after MINPERM'
!     WRITE(975,'(A,3F20.10)') ('LA ',XA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)

      DO J1 = 1, 3*NSIZE
         DISTANCE = DISTANCE + (XA(J1) - XB(J1))**2
      ENDDO
!
!     Update ALLPERM again. NEWPERM is the same as ALLPERM at the start of the next pass.
!
      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)

      IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,G20.10)') ' rbminpermdist> sites distance after permuting ',NPERM,'    &
     &                                        pairs of rigid bodies=', SQRT(DISTANCE)
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
         IF (RBSYMT) THEN
!
!     This call block does not do any further alignment so as not to break the standard reference geometry.
!
            CALL SITEPOS(DUMMYB,SITESB) 
            RBDISTANCE=1.0D10
            DO J1=1,NRB
               DO J2=1,NRBGROUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!                  CALL SITEPOS(DUMMYA,SITESA)
!                  WRITE(975,'(I6)') NRB*NRBSITES
!                  WRITE(975,'(A,2I6,A,F20.10)') 'A sites before rotation'
!                  WRITE(975,'(A,3F20.10)') ('LA ',SITESA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

                  P(:)      = DUMMYA(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3)
                  IF (J2 /= 1) THEN
                     CALL ROTMAT(P,RTEMP1)
                     PVEC(1:3) = MATMUL(RTEMP1,RBOPS(1:3,J2))
                     PVEC(1:3) = PVEC(1:3)/SQRT(DOT_PRODUCT(PVEC(1:3),PVEC(1:3)))
                     THETAH    = 0.5D0*RBOPS(4,J2) 
                     QTMP(1)   = COS(THETAH)
                     QTMP(2:4) = SIN(THETAH)*PVEC(1:3)
                     CALL QROTAA(QTMP,P)
                  ENDIF
                  DUMMYC(1:6*NRB) = DUMMYA(1:6*NRB)
                  DUMMYC(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3) = P(1:3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!                  CALL POTENTIAL(DUMMYC,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                  WRITE(MYUNIT,'(A,I6,2(A,F25.15))')' rbminpermdist> after RBROT following internal symmetry  & 
!     &                                          operation ',J2,' energy for A=',  ENERGY,' RMS=',RMS
!
!     Calculate new distance based on all sites metric rather than just centre-of-mass rigid body coordinates.
! 
                  CALL SITEPOS(DUMMYC,SITESA)
                  SITEDIST=0.D0
                  DO J3=1,3*NRB*NRBSITES
                     SITEDIST=SITEDIST+(SITESA(J3)-SITESB(J3))**2
                  ENDDO
!                  WRITE(975,'(I6)') NRB*NRBSITES
!                  WRITE(975,'(A,2I6,A,F20.10)') 'C sites for J1,J2=',J1,J2,' distance=',SQRT(SITEDIST)
!                  WRITE(975,'(A,3F20.10)') ('LA ',SITESA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

                  IF (SITEDIST.LT.RBDISTANCE) THEN
                     RBDISTANCE=SITEDIST
                     BESTA(1:6*NRB)=DUMMYC(1:6*NRB)
                  ENDIF
!                  WRITE(MYUNIT,'(A,2F20.10)')' rbpermdist> sites distance, RBDISTANCE, =',SQRT(SITEDIST), SQRT(RBDISTANCE)
               ENDDO
               DUMMYA(1:6*NRB)=BESTA(1:6*NRB)
            ENDDO
!
!     This call aligns the overall orientation with respect to the sites metric. We should already
!     have identified permutation-inversion isomers using the alignment based on centres of mass
!     above, if applicable. We now repeat but with an orientational alignment with respect to sites.
!
            CALL RBMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)
 
            CALL QROTQ(Q2,QCUMUL)

            DO J1=1,NRB
               DO J2=1,NRBGROUP
!
!     MUST NOT USE RMAT HERE - because RMAT is accumulated below.
! 
                  P(:)      = DUMMYA(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3)
                  IF (J2 /= 1) THEN
                     CALL ROTMAT(P,RTEMP1)
                     PVEC(1:3) = MATMUL(RTEMP1,RBOPS(1:3,J2))
                     PVEC(1:3) = PVEC(1:3)/SQRT(DOT_PRODUCT(PVEC(1:3),PVEC(1:3)))
                     THETAH    = 0.5D0*RBOPS(4,J2)
                     QTMP(1)   = COS(THETAH)
                     QTMP(2:4) = SIN(THETAH)*PVEC(1:3)
                     CALL QROTAA(QTMP,P)
                  ENDIF
                  DUMMYC(1:6*NRB) = DUMMYA(1:6*NRB)
                  DUMMYC(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3) = P(1:3)
!
!                 CALL POTENTIAL(DUMMYC,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                 WRITE(MYUNIT,'(A,I6,2(A,F25.15))')' rbminpermdist> after RBROT following internal symmetry 
!                                              operation ',J2,' energy for A=', ENERGY,' RMS=',RMS
!
!     Calculate new distance based on all sites metric rather than just centre-of-mass rigid body coordinates.
! 
                  CALL SITEPOS(DUMMYC,SITESA)
                  SITEDIST=0.D0
                  DO J3=1,3*NRB*NRBSITES
                     SITEDIST=SITEDIST+(SITESA(J3)-SITESB(J3))**2
                  ENDDO
!                 WRITE(MYUNIT,'(A,2F20.10)')' rbminpermdist> sites distance^2,RBDISTANCE=',SITEDIST,RBDISTANCE
                  IF (SITEDIST.LT.RBDISTANCE) THEN
                     RBDISTANCE=SITEDIST
                     BESTA(1:6*NRB)=DUMMYC(1:6*NRB)
                  ENDIF
               ENDDO
               DUMMYA(1:6*NRB)=BESTA(1:6*NRB)
            ENDDO
         ELSE
!
!     This call aligns the overall orientation with respect to the sites metric.
!     Internal symmetry operations of the rigid bodies are not considered.
!
            CALL RBMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)

         ENDIF

         DISTANCE  = DISTANCE*DISTANCE

!     accumulate rotation
         CALL QROTQ(Q2,QCUMUL)

         IF (NTRIES .LT. MAXIMUMTRIES) THEN
            GOTO 10
         ELSE ! prevent infinite loop
            IF (INVERT == -1) THEN
               IF(DEBUG) WRITE(MYUNIT,'(A)')' rbminpermdist> WARNING - number of tries exceeded, giving up'
            ELSE
               WRITE(MYUNIT,'(A)')' rbminpermdist> WARNING - number of tries exceeded, giving up'
            ENDIF
         ENDIF

      ENDIF

      IF (DISTANCE .LT. DBEST) THEN

         DBEST                = DISTANCE
         XBEST(1:3*NATOMS)    = DUMMYA(1:3*NATOMS)
         BESTPERM(1:NATOMS)   = ALLPERM(1:NATOMS)
         QTMP(1:4)  = QA(1:4)
         CALL QROTQ(QCUMUL,QTMP)
         QBEST(1:4) = QTMP(1:4)
         QBINV(1:4) = (/QB(1), -QB(2:4)/)

      ENDIF
!
! If GEOMDIFFTOLL is set too sloppy we could miss the best solution by exiting via the line
! below without trying other orbits. Turn off this escape?!
! The turn off seems to be a bug for non-RBPERM runs. LJ38 tests fail with all compilers, producing
! a "Distance is zero: this should not happen" error!!! DJW
! Put the escape back in for now.
!
      IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
      IF (NCHOOSE2.LT.NORBIT2) GOTO 30
      IF (NCHOOSE1.LT.NORBIT1) GOTO 65
      IF (EFIELDT) GOTO 50
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
      CALL RBROT(XBEST, XTMP, QBINV, NATOMS)
      XBEST(1:3*NATOMS) = XTMP(1:3*NATOMS)

      DO J1 = 1, (NATOMS/2)
         XBEST(3*(J1-1)+1) = XBEST(3*(J1-1)+1) + CMBX
         XBEST(3*(J1-1)+2) = XBEST(3*(J1-1)+2) + CMBY
         XBEST(3*(J1-1)+3) = XBEST(3*(J1-1)+3) + CMBZ
      ENDDO

!     Site-site distance

      CALL SITEPOS(XBEST,XA)

      CALL SITEPOS(COORDSB,XB)

      XDUMMY = 0.D0
      DO J1 = 1, 3*NSIZE
         XDUMMY = XDUMMY + (XA(J1) - XB(J1))**2
      ENDDO

      IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
          WRITE(MYUNIT,'(2(A,F20.10))')'rbminpermdist> ERROR *** distance between transformed XBEST and COORDSB=',  &
     &    SQRT(XDUMMY),  ' should be ', SQRT(DISTANCE)
          WRITE(MYUNIT,*)'COORDSB'
          DO J1 = 1, NATOMS
             J2 = 3*J1
             WRITE(MYUNIT,*) COORDSBS(J2-2), COORDSBS(J2-1), COORDSBS(J2)
          ENDDO 
!          WRITE(MYUNIT,*) 'COORDSA'
!          DO J1 = 1, NATOMS
!             J2 = 3*J1
!             WRITE(MYUNIT,*) COORDSAS(J2-2), COORDSAS(J2-1), COORDSAS(J2)
!          ENDDO 
          WRITE(MYUNIT,*) 'XBEST'
          DO J1 = 1, NATOMS
             J2 = 3*J1
             WRITE(MYUNIT,*) XBEST(J2-2), XBEST(J2-1), XBEST(J2)
          ENDDO 

!          WRITE(MYUNIT,*) 'RMSD=', SQRT(DISTANCE/DFLOAT(NATOMS*NRBSITES/2))
          STOP

      ENDIF

      CALL QROTQ(QBINV,QBEST)

      CALL QROTMAT(QBEST,RMATBEST)

      COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!      CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      WRITE(MYUNIT,'(2(A,F25.15))')' finally A energy=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

      DISTANCE=SQRT(DISTANCE)

      END SUBROUTINE RBMINPERMDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SITEPOS(Q, X)

!     GENERATE SITE POSITIONS FROM THE RIGID-BODY COORDINATES

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NTSITES, DBPTDT

      IMPLICIT NONE

      INTEGER :: J1, J2, J3, J4
      DOUBLE PRECISION :: Q(3*NATOMS), X(3*NTSITES), R(3), P(3), RM(3,3)
      DOUBLE PRECISION :: RBSITETD(4,3), SIGTD(4), SIGDB(2)

      IF (DBPTDT) CALL DEFTD(RBSITETD,SIGTD,SIGDB)

      DO J1 = 1, NATOMS/2

         J2   = 3*J1
         R(:) = Q(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = Q(J2-2:J2)

         CALL ROTMAT(P, RM)

         IF (DBPTDT) THEN
            IF (J1 < NATOMS/2) THEN
               DO J3 = 1, NRBSITES
                  J4          = 3*((J1-1)*NRBSITES + J3)
                  X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), SITE(J3,:))
               ENDDO
            ELSE
               DO J3 = 1, 4
                  J4          = 3*((J1-1)*NRBSITES + J3)
                  X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), RBSITETD(J3,:))
               ENDDO
            ENDIF
         ELSE
            DO J3 = 1, NRBSITES
               J4          = 3*((J1-1)*NRBSITES + J3)
               X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), SITE(J3,:))
            ENDDO
         ENDIF
      ENDDO

      END SUBROUTINE SITEPOS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBSITESORIENT(X, T1, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)

!      USE COMMONS, ONLY: NRBSITES, NTSITES, EFIELDT 
      USE COMMONS, ONLY: NTSITES, EFIELDT 
!                        
!     This subroutine puts permutational isomers of a rigid-body system into a standard orientation
!     with respect to the sites
!
      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, OFFSET, NSIZE
      DOUBLE PRECISION :: X(3*NATOMS), T1(3*NATOMS)
!      DOUBLE PRECISION :: XS(3*NATOMS*NRBSITES/2), T1S(3*NATOMS*NRBSITES/2), T2S(3*NATOMS*NRBSITES/2), DIST(NATOMS*NRBSITES/2)
      DOUBLE PRECISION :: XS(3*NTSITES), T1S(3*NTSITES), T2S(3*NTSITES), DIST(NTSITES)
      DOUBLE PRECISION :: AX(3), P(3), Q2(4), ROTM(3,3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, COST, SINT, COSTH, SINTH, ST, FCT
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, DUMMY, PROJ, DMAX2, CUTOFF1, DTEMP
      LOGICAL          :: DEBUG

      NSIZE = NTSITES
      OFFSET  = 3*NATOMS/2
      CUTOFF1 = 1.D-03

      CALL SITEPOS(X, XS)
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NSIZE
         J = 3*I
         CMX = CMX + XS(J-2)
         CMY = CMY + XS(J-1)
         CMZ = CMZ + XS(J)
      ENDDO

      CMX = CMX/NSIZE; CMY = CMY/NSIZE; CMZ = CMZ/NSIZE
      DO I = 1, NSIZE
         J = 3*I
         XS(J-2) = XS(J-2) - CMX
         XS(J-1) = XS(J-1) - CMY
         XS(J)   = XS(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN
         T1S(:) = XS(:)
         GOTO 100
      ENDIF
!
!     Find the rbsite which is at the largest distance from the centre of mass
!
      DO J1 = 1, NSIZE

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

!     For tagged rbsites, the choice of the first rbsite matters if it belongs to an orbit of size > 1.
!
      IF ((ABS(XS(3*JMAX1-2)) < 1.D-08) .AND. (ABS(XS(3*JMAX1-1)) < 1.D-08)) THEN 
!
!     rb JMAX1 is already on the z axis!
!
         IF (XS(3*(JMAX1-1)+3) > 0.D0) THEN

            T1S(1:3*NSIZE) = XS(1:3*NSIZE)
            Q2(1:4)   = (/1.D0, 0.D0, 0.D0, 0.D0/)   ! Identity operation       
                  
         ELSE  ! rotate about the x-axis by \pi, DO NOT INVERT!

!            P(:) = (/4.D0*DATAN(1.D0), 0.D0, 0.D0/)
!            CALL ROTMAT(P(:), ROTM(:,:))
            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)      ! the corresponding quaternion

            DO J1 = 1, NATOMS/2
               J       = 3*J1
               T1S(J-2) = XS(J-2)
               T1S(J-1) =-XS(J-1)
               T1S(J)   =-XS(J)
            ENDDO

         ENDIF

      ELSE
!
!     Rotate all rb's through an angle THETA about the axis P(:), so as to rotate the rbsite JMAX1 onto the z axis 
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
          
         DO J1 = 1, NSIZE
            J2 = 3*J1
            T1S(J2-2:J2) = MATMUL(ROTM(:,:),XS(J2-2:J2))
         ENDDO   

      ENDIF

!
!     Now find the rbsite with the largest distance from the z-axis 
!
100   DMAX = -1.0D0

      DO J1 = 1, NSIZE
         J2 = 3*J1 
         DIST(J1) = SQRT(T1S(J2-2)**2 + T1S(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ is the sum of the x components. Use T2S as a dummy in order not to 
!     change T1 until we have decided which atom to put in the xz plane.
!
      DO J1 = 1, NSIZE

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2S(1:3*NSIZE) = T1S(1:3*NSIZE)

            CALL RBSITESROTXZ(NSIZE, J1, T2S, PROJ, DIST, Q2, .FALSE.)

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

!      WRITE(MYUNIT,'(A,I6,F20.10)') 'JMAX2,DMAX2=', JMAX2,DMAX2
!
!     and now rotate it into the xz plane.
!
      CALL RBSITESROTXZ(NSIZE, JMAX2, T1S, DTEMP, DIST, Q2, .TRUE.)

!     Now change the coordinates by the rotation matrix corresponding to the net transformation

      CALL RBROT(X, T1, Q2, NATOMS)

      END SUBROUTINE RBSITESORIENT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBSITESROTXZ(NSIZE, JDO, T1S, PROJ, DIST, Q1, ROTT)

      USE COMMONS, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NSIZE, JDO, I, J, J2
      DOUBLE PRECISION :: T1S(3*NSIZE), PROJ, DIST(NSIZE), THETA, THETAH, COST, SINT, ST
      DOUBLE PRECISION :: COSTH, SINTH, FCT, P(3), Q2(4), Q1(4), RM(3,3), ROTM(3,3) 
      LOGICAL          :: ROTT

      J2     = 3*JDO
 
      IF (ABS(T1S(J2-1)) < 1.0D-8) THEN            ! already on the xz plane

         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/) 

         IF (T1S(J2-2) < 0.D0) THEN                ! rotate about the z axis by \pi, DO NOT INVERT!!

            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/) 

            DO I = 1, NSIZE
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

         DO I = 1, NSIZE
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
!            ROTM(:,:) = MATMUL(RM(:,:),ROTM(:,:))
            CALL QROTQ(Q2,Q1)
         ENDIF

      ELSE

         PROJ = 0.D0

         DO I = 1, NSIZE

            J = 3*I
            IF (T1S(J) > 1.0D-2) PROJ = PROJ + T1S(J-2)

         ENDDO

      ENDIF

      END SUBROUTINE RBSITESROTXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE INVMTRX(S, SINV)
!     inverts a 3x3 matrix S

      IMPLICIT NONE

      DOUBLE PRECISION   :: S(3,3), SINV(3,3), A(3,3), DET

!     ADJOINT

      A(1,1) = S(2,2)*S(3,3) - S(2,3)*S(3,2)
      A(1,2) = S(1,3)*S(3,2) - S(1,2)*S(3,3)
      A(1,3) = S(1,2)*S(2,3) - S(1,3)*S(2,2)
      A(2,1) = S(3,1)*S(2,3) - S(2,1)*S(3,3)
      A(2,2) = S(1,1)*S(3,3) - S(1,3)*S(3,1)
      A(2,3) = S(1,3)*S(2,1) - S(1,1)*S(2,3)
      A(3,1) = S(2,1)*S(3,2) - S(3,1)*S(2,2)
      A(3,2) = S(1,2)*S(3,1) - S(1,1)*S(3,2)
      A(3,3) = S(1,1)*S(2,2) - S(1,2)*S(2,1)

      DET    =  S(1,1)* A(1,1) + S(1,2) * A(2,1) + S(1,3) * A(3,1)

      SINV(:,:) = A(:,:)/DET

      END SUBROUTINE INVMTRX

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBROT(T, X, Q2, NATOMS)
!     takes the set of rigid-body coordinates T and returns X after rotation via the quaternion Q2 
!     about the origin

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: T(3*NATOMS), X(3*NATOMS), T1(1:3), Q2(4), P(3), RM(3,3) 
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS/2
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = 2*CMX/NATOMS; CMY = 2*CMY/NATOMS; CMZ = 2*CMZ/NATOMS
      DO I = 1, NATOMS/2
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
         X(J)   = T(J)   - CMZ
      ENDDO

!     extract the rotation matrix corresponding to rotation via Q

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS/2

         J        = 3*I
         T1(1:3)  = MATMUL(RM(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)  
!     CONVERT THE ANGLE-AXIS COORDINATES

         J        = 3*NATOMS/2 + J
         P(:)     = T(J-2:J)
         CALL QROTAA(Q2,P)
         X(J-2:J) = P(1:3)

      ENDDO

      END SUBROUTINE RBROT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTAA(Q2,P)

!     transforms the angle-axis variables P corresponding to rotation via quaternion Q2 and returns 
!     the transformed variables in P

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), Q1(4), Q2(4), P(3), THETA, THETAH, ST, FCT

      THETA   = DSQRT(DOT_PRODUCT(P,P))
      THETAH  = 0.5D0*THETA
      ST      = SIN(THETAH)
      Q1(1)   = COS(THETAH)
      Q1(2:4) = P(:)*ST/THETA

      Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
      Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
      Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
      Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

      THETA  = 2.D0*ACOS(Q(1))

      IF (THETA == 0.D0) THEN
         P (1:3) = 0.D0
      ELSE
         FCT = DSQRT(DOT_PRODUCT(Q(2:4),Q(2:4)))
         P(1:3) = THETA*Q(2:4)/FCT
      ENDIF

      END SUBROUTINE QROTAA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTQ(Q2,Q1)

!     transforms the quaternion Q1 corresponding to rotation via quaternion Q2 and returns the
!     transformed variables in Q1

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), Q1(4), Q2(4)

      Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
      Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
      Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
      Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

      Q1(1:4) = Q(1:4)

      END SUBROUTINE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTMAT(Q,RM)

!     provides the rotation matrix RM corresponding to quaternion Q
!     right-handed rotation in the right-handed coordinate system

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), RM(3,3)

      RM(1,1) = Q(1)**2 + Q(2)**2 - Q(3)**2 - Q(4)**2
      RM(1,2) = 2.D0*(Q(2)*Q(3) - Q(1)*Q(4))
      RM(1,3) = 2.D0*(Q(2)*Q(4) + Q(1)*Q(3))
      RM(2,1) = 2.D0*(Q(2)*Q(3) + Q(1)*Q(4))
      RM(2,2) = Q(1)**2 + Q(3)**2 - Q(2)**2 - Q(4)**2
      RM(2,3) = 2.D0*(Q(3)*Q(4) - Q(1)*Q(2))
      RM(3,1) = 2.D0*(Q(2)*Q(4) - Q(1)*Q(3))
      RM(3,2) = 2.D0*(Q(3)*Q(4) + Q(1)*Q(2))
      RM(3,3) = Q(1)**2 + Q(4)**2 - Q(2)**2 - Q(3)**2

      END SUBROUTINE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBORIENT(X, T1, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)
      
      USE COMMONS, ONLY : EFIELDT 
!                        
!     This subroutine puts permutational isomers of a rigid-body system into a standard orientation.
!
      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), DIST(NATOMS), T1(3*NATOMS), T2(3*NATOMS)
      DOUBLE PRECISION :: RM(3,3), AX(3), P(3), Q1(4), Q2(4), Q(4), ROTM(3,3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, COST, SINT, COSTH, SINTH, ST, FCT
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, DUMMY, PROJ, DMAX2, CUTOFF1, DTEMP
      LOGICAL          :: DEBUG

      OFFSET  = 3*NATOMS/2
      CUTOFF1 = 1.D-03
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS/2
         J = 3*I
         CMX = CMX + X(J-2)
         CMY = CMY + X(J-1)
         CMZ = CMZ + X(J)
      ENDDO

      CMX = 2*CMX/NATOMS; CMY = 2*CMY/NATOMS; CMZ = 2*CMZ/NATOMS
      DO I = 1, NATOMS/2
         J = 3*I
         X(J-2) = X(J-2) - CMX
         X(J-1) = X(J-1) - CMY
         X(J)   = X(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN

         T1(:) = X(:)
         GOTO 100

      ENDIF
!
!     Find the rb which is at the largest distance from the centre of mass
!
      DO J1 = 1, NATOMS/2

         J = 3*J1
         DIST(J1) = SQRT(X(J-2)**2 + X(J-1)**2 + X(J)**2)

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            NORBIT1 = NORBIT1 + 1
            IF (NORBIT1 == NCHOOSE1) JMAX1 = J1

         ELSE IF (DIST(J1) > DMAX) THEN

            DMAX = DIST(J1)
            NORBIT1 = 1
            JMAX1 = J1

         ENDIF

      ENDDO

!     For tagged rbs, the choice of the first rb matters if it belongs to an orbit of size > 1.
!
      IF ((ABS(X(3*JMAX1-2)) < 1.D-08) .AND. (ABS(X(3*JMAX1-1)) < 1.D-08)) THEN 
!
!     rb JMAX1 is already on the z axis!
!
         IF (X(3*(JMAX1-1)+3) > 0.D0) THEN

            T1(1:3*NATOMS) = X(1:3*NATOMS)
            Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)
 
         ELSE  ! rotate about the x-axis by \pi, DO NOT INVERT!

            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)

            DO J1 = 1, NATOMS/2
               J       = 3*J1
               T1(J-2) = X(J-2)
               T1(J-1) =-X(J-1)
               T1(J)   =-X(J)
            ENDDO

         ENDIF

      ELSE
!
!     Rotate all rb's through an angle THETA about the axis P(:), so as to rotate the rb JMAX1 onto the z axis 
!
!     For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!

        THETA = DACOS(X(3*JMAX1)/DIST(JMAX1))

!     The axis is on the xy-plane and perpendicular to the projection of the translation vector of the rb JMAX1
!     on the xy-plane. The axis is so chosen that a right-handed rotation in the right-handed coordinate system
!     will result in the required transformation.

         AX(3) = 0.D0
         FCT   = DSQRT(X(3*JMAX1-2)**2 + X(3*JMAX1-1)**2)
         AX(1) = THETA*X(3*JMAX1-1)/FCT
         AX(2) =-THETA*X(3*JMAX1-2)/FCT

         CALL ROTMAT(AX(:), ROTM(:,:))

         THETAH  = 0.5D0*THETA
         AX(1:3) = AX(1:3)/SQRT(DOT_PRODUCT(AX(1:3),AX(1:3)))
         Q2(1:4) = (/COS(THETAH), SIN(THETAH)*AX(1:3)/)          ! the corresponding quaternion

         CALL RBROT(X, T1, Q2, NATOMS)

      ENDIF

!
!     Now find the rb with the largest distance from the z-axis 
!
100   DMAX = -1.0D0

      DO J1 = 1, NATOMS/2

         J2 = 3*J1 
         DIST(J1) = SQRT(T1(J2-2)**2 + T1(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
         
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ is the sum of the x components. Use T2 as a dummy in order not to 
!     change T1 until we have decided which atom to put in the xz plane.
!
      DO J1 = 1, NATOMS/2

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2(1:3*NATOMS) = T1(1:3*NATOMS)
            CALL RBROTXZ(NATOMS, J1, T2, PROJ, DIST, Q2, .FALSE.)

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
!      CALL RBROTXZ(NATOMS, JMAX2, T1, DMAX2, DIST, ROTM, .TRUE.)

      CALL RBROTXZ(NATOMS, JMAX2, T1, DTEMP, DIST, Q2, .TRUE.)

!      CALL INVMTRX(ROTM(:,:), ROTMINV(:,:))

      END SUBROUTINE RBORIENT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBROTXZ(NATOMS, JDO, T1, PROJ, DIST, Q1, ROTT)

      USE COMMONS, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, JDO, I, J, J2, OFFSET
      DOUBLE PRECISION :: T1(3*NATOMS), X(3*NATOMS), PROJ, DIST(NATOMS/2), THETA, THETAH, COST, SINT, ST
      DOUBLE PRECISION :: COSTH, SINTH, FCT, P(3), Q2(4), Q1(4), RM(3,3), ROTM(3,3), CMX, CMY, CMZ
      LOGICAL          :: ROTT

      J2     = 3*JDO
      OFFSET = 3*NATOMS/2

      X(1:3*NATOMS) = T1(1:3*NATOMS)

      IF (ABS(T1(J2-1)) < 1.0D-8) THEN            ! already on the xz plane

         RM(:,:) = 0.D0
         RM(1,1) = 1.D0; RM(2,2) = 1.D0; RM(3,3) = 1.D0
         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)

         IF (T1(J2-2) < 0.D0) THEN                ! rotate about the z axis by \pi, DO NOT INVERT!!

            P(:) = (/0.D0, 0.D0, 4.D0*DATAN(1.D0)/)
            CALL ROTMAT(P(:), RM(:,:))
            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/)

            DO I = 1, NATOMS/2
               J = 3*I
               X(J-2) = -T1(J-2)
               X(J-1) = -T1(J-1)
            ENDDO

         ENDIF

      ELSE

         COST   = T1(J2-2)/DIST(JDO)
         SINT   = T1(J2-1)/DIST(JDO)
         THETA  =-ATAN2(SINT,COST)              ! the negative sign appears as the rotation is reversed
         THETAH = 0.5D0*THETA
         Q2(1:4) = (/COS(THETAH), 0.D0, 0.D0, SIN(THETAH)/)
         RM(:,:) = 0.D0                           ! the rotation matrix for R_{z}(-\theta)
         RM(1,1) = COST
         RM(2,2) = COST
         RM(3,3) = 1.D0
         RM(1,2) = SINT
         RM(2,1) =-SINT

         DO I = 1, NATOMS/2
            IF (DIST(I) /= 0.0D0) THEN
               J = 3*I
               X(J-2:J) = MATMUL(RM(:,:),T1(J-2:J))
            ENDIF
         ENDDO

      ENDIF

      IF (ROTT) THEN

         IF (EFIELDT) THEN
            ROTM(:,:) = RM(:,:)
         ENDIF

!     Now change the angle-axis coordinates via quaternion multiplication equivalent to rotation by ROTM

         CALL RBROT(T1, X, Q2, NATOMS)
         T1(1:3*NATOMS) = X(1:3*NATOMS)

      ELSE

         PROJ = 0.D0
         DO I = 1, NATOMS/2
            J = 3*I
            IF (X(J) > 1.0D-2) PROJ = PROJ + X(J-2)
         ENDDO

      ENDIF

      T1(1:3*NATOMS) = X(1:3*NATOMS)

      END SUBROUTINE RBROTXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBROTMAT(T, X, ROTMAT, NTSITES)
!     takes the set of rigid-body coordinates T and returns X after rotation via the quaternion Q2
!     about the origin

      IMPLICIT NONE

      INTEGER          :: I, J, NTSITES
      DOUBLE PRECISION :: T(3*NTSITES), X(3*NTSITES), T1(1:3), Q2(4), ROTMAT(3,3)
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NTSITES
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = CMX/NTSITES; CMY = CMY/NTSITES; CMZ = CMZ/NTSITES
      DO I = 1, NTSITES
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
         X(J)   = T(J)   - CMZ
      ENDDO

!     extract the rotation matrix corresponding to rotation via Q

!      CALL QROTMAT(Q2,RM)

      DO I = 1, NTSITES

         J        = 3*I
         T1(1:3)  = MATMUL(ROTMAT(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)

      ENDDO

      END SUBROUTINE RBROTMAT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBREFLECTYZ(X)
! reflects rigidbody coordinates X in the yz plane
      USE COMMONS, ONLY : NATOMS
      USE VEC3
      USE ROTATIONS
      IMPLICIT NONE
      INTEGER          :: J1,J2,J3
      DOUBLE PRECISION :: X(3*NATOMS),M(3,3),I(3,3),A(3),P(3)

      CALL IDENTITY3X3(I(:,:))

      DO J1=1,NATOMS/2
! reflect COM coordinates
         J2=3*J1
         X(J2-2)=(-1.D0)*X(J2-2)
! reflect AA coordinates
         J2=3*(J1+NATOMS/2)
         M(:,:)=ROT_AA2MX(X(J2-2:J2))
         A(:)=VEC_CROSS( MATMUL(M(:,:),I(2,:)), MATMUL(M(:,:),I(3,:)) )
         P(:)=2.D0*ACOS(A(1)/VEC_LEN(A(:)))*(/0.D0,A(3),-1.D0*A(2)/)/SQRT(A(2)**2+A(3)**2)
         X(J2-2:J2)=ROT_ROTATE_AA(X(J2-2:J2),P(:))
      ENDDO

      END SUBROUTINE RBREFLECTYZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBIMPROPERROTATION(Y,X,RM)
! applies the improper rotation RM to rigidbody coordinates X, yielding Y
!    splits RM into a reflection in yz, RFYZ, followed by a pure rotation, RROT
!    reflects the system in yz with RBREFLECTYZ
!    rotates the system with RROT
      USE COMMONS, ONLY : NATOMS
      USE VEC3
      USE ROTATIONS
      IMPLICIT NONE
      INTEGER                         :: J1,J2
      DOUBLE PRECISION                :: RFYZ(3,3), RROT(3,3), P(3)
      DOUBLE PRECISION, INTENT (IN) :: RM(3,3), X(NATOMS*3)
      DOUBLE PRECISION, INTENT(OUT) :: Y(NATOMS*3)

      Y(:)=X(:)
      RROT(:,:)=RM(:,:)
      IF (DET3x3(RM).LT.0.D0) THEN ! RM is an improper rotation
        CALL IDENTITY3x3(RFYZ(:,:))
        RFYZ(1,1)=-1.0D0
        RROT(:,:)=MATMUL(RM(:,:),RFYZ(:,:))
!        RROT(:,:)=MATMUL(RFYZ(:,:),RM(:,:))
        CALL RBREFLECTYZ(Y(:))
      ENDIF
! rotate COM coordinates
      Y(:3*NATOMS/2)=RESHAPE(MATMUL(RROT(:,:),RESHAPE(Y(:3*NATOMS/2),(/ 3, NATOMS/2 /))),(/ 3*NATOMS/2 /))
! rotate AA coordinates
      P(:)=ROT_MX2AA(RROT(:,:))
      DO J1=1,NATOMS/2
         J2=(J1+NATOMS/2)*3
         Y(J2-2:J2)=ROT_ROTATE_AA(Y(J2-2:J2),P(:))
      ENDDO

      END SUBROUTINE RBIMPROPERROTATION

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBSITESORIENTPCA(X,XR,NATOMS,Q)
      USE COMMONS, ONLY: NTSITES
      USE ROTATIONS
      USE VEC3
      IMPLICIT NONE
      INTEGER            :: NATOMS, NSIZE, J1, J2, OFFSET, INFO
      DOUBLE PRECISION   :: X(3*NATOMS), XR(3*NATOMS), XS(3*NTSITES), XST(3,NTSITES), XSTMP(3*NTSITES)
      DOUBLE PRECISION   :: COV(3,3), EVAL(3), WORK(102), RMAT(3,3), P(3), Q(4), EYE(3,3)
      DOUBLE PRECISION   :: CMX, CMY, CMZ

      NSIZE = NTSITES
      OFFSET  = 3*NATOMS/2

      CALL SITEPOS(X, XS)
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO J1 = 1, NSIZE
         J2 = 3*J1
         CMX = CMX + XS(J2-2)
         CMY = CMY + XS(J2-1)
         CMZ = CMZ + XS(J2)
      ENDDO

      DO J1 = 1, NSIZE
         J2 = 3*J1
         XS(J2-2) = XS(J2-2) - CMX/NSIZE
         XS(J2-1) = XS(J2-1) - CMY/NSIZE
         XS(J2)   = XS(J2)   - CMZ/NSIZE
      ENDDO
!
!     get eigenvectors of covariance matrix
!
      XST = RESHAPE(XS,(/3,NTSITES/))
      COV=MATMUL(XST,TRANSPOSE(XST))
      RMAT=COV
      CALL DSYEV('V','L',3,RMAT,3,EVAL,WORK,102,INFO)
      RMAT=TRANSPOSE(RMAT)
!
!     force the transformation to place the cubic centre of mass in the
!     top-front-right (+,+,+) or bottom-back-left (-,-,-) octants (whichever
!     doesn't require a reflection)
!
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0

      DO J1=1,NSIZE
         J2=3*J1
         XSTMP(J2-2:J2)=MATMUL(RMAT(:,:),XS(J2-2:J2))
         CMX=CMX+XSTMP(J2-2)**3
         CMY=CMY+XSTMP(J2-1)**3
         CMZ=CMZ+XSTMP(J2-0)**3
      ENDDO

      CALL IDENTITY3X3(EYE)
      EYE(1,1)=CMX/ABS(CMX)
      EYE(2,2)=CMY/ABS(CMY)
      EYE(3,3)=CMZ/ABS(CMZ)
      RMAT=MATMUL(EYE,RMAT)
      RMAT=DET3X3(RMAT)*RMAT
!
!     generate the corresponding quaternion and angle-axis rotations
!
      Q(:)=ROT_MX2Q(RMAT(:,:))
      P(:)=ROT_Q2AA(Q(:))
!
!     put the transformed coordinates into XR
!
      DO J1=1,NATOMS/2
         J2=3*J1
         XR(J2-2:J2)=MATMUL(RMAT(:,:),X(J2-2:J2))
         J2=3*(J1+NATOMS/2)
         XR(J2-2:J2)=ROT_ROTATE_AA(X(J2-2:J2),P(:))
      ENDDO

      END SUBROUTINE RBSITESORIENTPCA

!     ----------------------------------------------------------------------------------------------
