!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAP                                                                                          !
!                                                                                              !
! Calculates Energy and gradients for the patch-antipatch potential, using the rigid body      !
! angle axis framework described in:                                                           !
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and          !
! David Wales, Phys. Chem. Chem. Phys., 2009, 11, 1970-1976                                    !
! X: the positions and orientations of the bodies, and the three lattice vectors               !
! G: the gradients of the potential energy surface for changes in positions and orientations,  !
!    and the lattice derivtives                                                                !
! ENERGY: the potential energy for the configuration stored in X                               !
! GTEST: logical, true if gradients need to be calculated                                      !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAP_LATTICE(X, G, ENERGY, GTEST)

! NATOMS: twice the number of bodies, one for position and one for orientation
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the directions of the patch and antipatch vectors
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, BOXLX, BOXLY, BOXLZ

      IMPLICIT NONE

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY
      LOGICAL :: GTEST

! I, J: iterators over bodies
! IT, IR: indices of translational and rotational coordinates of I
! JT, JR: indices of translational and rotational coordinates of J
! S, T: indices of patches in body I and J
! S2: iterator over patches in I
! NMOL: the number of bodies
! OFFSET: the indicial offset to the start of orientations in X
      INTEGER :: I, IT, IR, J, JT, JR, S, S2, T, NMOL, OFFSET

! E: the orientations of patches and antipatches in the lab frame
! P: the rotation vector for a body
      DOUBLE PRECISION :: E((NATOMS-1)*NRBSITES/2,3), P(3)

! DEk: matrix product of the derivative of the rotation matrix with respect to kth component of the rotation vector and the
!      position of each site in the reference geometry, as in equation (12) of the paper
      DOUBLE PRECISION :: DE1((NATOMS-1)*NRBSITES/2,3), DE2((NATOMS-1)*NRBSITES/2,3), DE3((NATOMS-1)*NRBSITES/2,3) 

! RMI: the rotation matrix for body I
! DRMIk: the derivative RMI with respect to the kth component of the rotation vector
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

! IJENERGY: energy of interaction between bodies I and J
! GIT: gradients of translational coordinates of I
! GJT: gradients of translational coordinates of J
! GIR: gradients of rotational coordinates of I
! GJR: gradients of rotational coordinates of J
      DOUBLE PRECISION :: IJENERGY, GIT(3), GJT(3), GIR(3), GJR(3)

! PER: loops over coordinate configurations for periodic boundary conditions
! XI: modified coordinates of body I for periodic boundary conditions
! RIJ: vector displacement between bodies I and J
! DSS2: squared distance between I and J
! RCUT: cutoff distance for periodic LJ potential
! RCUTSQ: RCUT squared
      INTEGER          :: PER
      DOUBLE PRECISION :: XI(3), RIJ(3), DSS2, RCUT, RCUTSQ

! BOXLk: the dimensions of the periodic box
     ! DOUBLE PRECISION :: BOXLX, BOXLY, BOXLZ

! Initialise energy and gradients
      ENERGY  = 0.D0
      IF (GTEST) THEN
        G(:) = 0.D0
      ENDIF

! Find NMOL as the actual number of bodies, and OFFSET as the start of rotational coordinates
      NMOL    = (NATOMS-1)/2
      OFFSET  = 3*NMOL
! Set the box lengths
      BOXLX = X(6*NMOL+1)
      BOXLY = X(6*NMOL+2)
      BOXLZ = X(6*NMOL+3)
      
! I Loop over all bodies  
      DO I = 1, NMOL

! Set IT as the index of the third translational coordinate of the Ith body, and IR as the index of the third rotational
! coordinate of the Ith body
! Set P to be the rotation vector of body I
        IT = 3*I
        IR = OFFSET + IT
        P  = X(IR-2:IR)

! Calculate the rotation matrix and derivatives thereof
        CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! S2 Loop over all patches and antipatches for body I
        DO S2 = 1, NRBSITES

! Set S as the index of the S2th patch of body I
          S = NRBSITES*(I-1) + S2
! Calculate E from the rotation matrix acting on the patch orientation in the body frame
          E(S,:) = MATMUL(RMI(:,:),RBSTLA(S2,:))

          IF (GTEST) THEN
! If derivatives are required, calculate the kth derivative of the rotation matrix for body I acting on the orientation of
! the Sth patch in the body frame
            DE1(S,:) = MATMUL(DRMI1(:,:),RBSTLA(S2,:))
            DE2(S,:) = MATMUL(DRMI2(:,:),RBSTLA(S2,:))
            DE3(S,:) = MATMUL(DRMI3(:,:),RBSTLA(S2,:))
          ENDIF ! End IF GTEST
        ENDDO ! End loop over patches
      ENDDO ! End loop over bodies

! Periodic Boundary Conditions
! Calculate the cutoff distance of the potential
      RCUT = BOXLX
      IF (BOXLY .LT. RCUT) RCUT = BOXLY
      IF (BOXLZ .LT. RCUT) RCUT = BOXLZ
      RCUT=5.0 ! 2.702*2.0
      RCUTSQ = (RCUT/2.D0)**2

! Make sure all coordinates lie within the box (which extends from 0 to BOXLk in all 3 directions
      DO I = 1, NMOL
        IT=I*3 ! X(J3) is z comp of position of mol J1
        X(IT-2)=X(IT-2)-BOXLX*IDNINT(X(IT-2)/BOXLX)
        X(IT-1)=X(IT-1)-BOXLY*IDNINT(X(IT-1)/BOXLY)
        X(IT)  =X(IT)  -BOXLZ*IDNINT(X(IT)  /BOXLZ)
      ENDDO

! I Loop over all bodies, except the last, to perform the double sum over bodies without double counting
      DO I = 1, NMOL - 1  

! Set IT as the index of the third translational coordinate of I, and IR as the index of the third rotational
! coordinate of I
        IT = 3*I
        IR = OFFSET + IT

! J Loop over bodies index greater than I, to perform the double sum over bodies without double counting
        DO J = I + 1, NMOL

! Set JT as the index of the third translational coordinate of J, and JR as the index of the third rotational
! coordinate of J
          JT = 3*J
          JR = OFFSET + JT

! Set S and T as appropriate indices of patches in bodies I and J
          S = (I-1)*NRBSITES
          T = (J-1)*NRBSITES

! Loop over different coordinate possibilities
          DO PER = 0, 7
! Set the coordinates
            XI(:) = X(IT-2:IT)
            IF (BTEST(PER,0)) XI(1) = XI(1) + (-1.D0)**IDNINT(XI(1)/BOXLX + 0.5D0)*BOXLX
            IF (BTEST(PER,1)) XI(2) = XI(2) + (-1.D0)**IDNINT(XI(2)/BOXLY + 0.5D0)*BOXLY
            IF (BTEST(PER,2)) XI(3) = XI(3) + (-1.D0)**IDNINT(XI(3)/BOXLZ + 0.5D0)*BOXLZ

! Find the distance between the bodies. Due to the choice of cutoff, only
! one distance is less than the cutoff, so only one of the eight interactions is non-zero.
            RIJ(:) = XI(:) - X(JT-2:JT)
            DSS2   = DOT_PRODUCT(RIJ(:),RIJ(:))

            IF (DSS2 .LE. RCUTSQ) THEN
! We have found the interacting coordinates. Calculate energy and gradients.
              CALL PAPLJPERIODIC(IJENERGY, GIT, GJT, GTEST, RIJ(:), RCUTSQ)
              ENERGY = ENERGY + IJENERGY
              IF (GTEST) THEN
                G(IT-2:IT) = G(IT-2:IT) + GIT(:)
                G(JT-2:JT) = G(JT-2:JT) + GJT(:)
! Lattice derivatives
                G(NMOL*6+1) = G(NMOL*6+1) + GIT(1)*RIJ(1)/BOXLX
                G(NMOL*6+2) = G(NMOL*6+2) + GIT(2)*RIJ(2)/BOXLY
                G(NMOL*6+3) = G(NMOL*6+3) + GIT(3)*RIJ(3)/BOXLZ
              ENDIF
              CALL PAPENERGYPAP(IJENERGY, GIT, GJT, GIR, GJR, GTEST, XI(:), X(JT-2:JT),              &
     &                          E(S+1:S+NRBSITES,:), E(T+1:T+NRBSITES,:), DE1(S+1:S+NRBSITES,:),     &
     &                          DE2(S+1:S+NRBSITES,:), DE3(S+1:S+NRBSITES,:), DE1(T+1:T+NRBSITES,:), &
     &                          DE2(T+1:T+NRBSITES,:), DE3(T+1:T+NRBSITES,:))
              ENERGY = ENERGY + IJENERGY
              IF (GTEST) THEN
                G(IT-2:IT) = G(IT-2:IT) + GIT(:)
                G(JT-2:JT) = G(JT-2:JT) + GJT(:)
                G(IR-2:IR) = G(IR-2:IR) + GIR(:)
                G(JR-2:JR) = G(JR-2:JR) + GJR(:)
! Lattice derivatives
                G(NMOL*6+1) = G(NMOL*6+1) + GIT(1)*RIJ(1)/BOXLX
                G(NMOL*6+2) = G(NMOL*6+2) + GIT(2)*RIJ(2)/BOXLY
                G(NMOL*6+3) = G(NMOL*6+3) + GIT(3)*RIJ(3)/BOXLZ
              ENDIF 

! No point continuing to loop over coordinate configurations now that we have found the interacting one
              EXIT

            ENDIF ! End if found best configuration
          ENDDO ! End loop over coordinate configurations
        ENDDO ! End loop over particles J 
      ENDDO ! End loop over particles I

      END SUBROUTINE PAP_LATTICE
       
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! VIEWPAP                                                                                      !
!                                                                                              !
! Writes positions to file in a format readable by XMakeMol                                    !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE VIEWPAP_LATTICE()

! NATOMS: twice the number of bodies
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the patch and antipatch directions in the body frame
! NSAVE: the number of different configurations to be written
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, NSAVE

! Global variable from QMODULE:
! QMIN: the energies of the stored configurations
! FF: the steps at which the configurations were first found
! QMINP: two dimensional array storing the positions, then orientations, for each configuration
      USE QMODULE

      IMPLICIT NONE

! I, J1-7: indices and counters to iterate over
      INTEGER          :: J1, J2, J3, J5, J7

! RMI: the rotation matrix for a body
! DRMI: dummy variable
! P: the rotation vector for a body
! RBCOORDS: the position of a patch or antipatch for a body
! LFCTR: length factor, distance of patches from centre of body
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3), LFCTR

! BOXLk: box length in kth direction
      DOUBLE PRECISION :: BOXLX, BOXLY, BOXLZ

! GTEST: indicates whether to calculate gradients or not
      LOGICAL          :: GTEST

! Choose vale of LFCTR to look pretty
      LFCTR = 0.4D0

! Open writer to file 'pap.xyz'
      OPEN(UNIT=26, FILE='pap.xyz', STATUS='UNKNOWN')

! Gradients are not required
      GTEST = .FALSE.

! J1 Loop over each configuration
      DO J1 = 1, NSAVE

! Write the number of cores, patches and antipatches 
         WRITE(26,'(I6)') ((NATOMS-1)/2)*(NRBSITES+1)
! Write the energy and step at which the configuration was first found
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

! Set the box lengths
         BOXLX = QMINP(J1, NATOMS*3-2)
         BOXLY = QMINP(J1, NATOMS*3-1)
         BOXLZ = QMINP(J1, NATOMS*3-0)

! J3 Loop over each body
         DO J3 = 1, (NATOMS-1)/2

! Set J5 as the index of the third translational coordinate of body J3, and J7 as the third rotational coordinate
            J5   = 3*J3
            J7   = 3*(NATOMS-1)/2 + J5
! Set P as the rotational coordinates
            P(:) = QMINP(J1,J7-2:J7)

! Find the rotation matrix
! Since GTEST is false, DRMI is a dummy variable
            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

! Write the position of the cores of the bodies
! Display as an Oxygen atom
! If this is the first body, and periodic boundary conditions are used, include the box
            IF (J3 .EQ. 1) THEN
              WRITE(26,'(A4,3F20.10,A10,6F20.10)') 'O', QMINP(J1,J5-2), QMINP(J1,J5-1), QMINP(J1,J5), &
     &        'bbox_xyz', (-BOXLX/2.D0), (BOXLX/2.D0), (-BOXLY/2.D0), (BOXLY/2.D0), (-BOXLZ/2.D0), (BOXLZ/2.D0)
            ELSE
              WRITE(26,'(A4,3F20.10)') 'O', QMINP(J1,J5-2), QMINP(J1,J5-1), QMINP(J1,J5)
            ENDIF

! J2 Loop over patches and antipatches within body J3
            DO J2 = 1, NRBSITES

! Calculate the position of the patch J2 in the lab frame
              RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + LFCTR*MATMUL(RMI(:,:),RBSTLA(J2,:))

! jwrm2> For patches and antipatches different, use this block
               IF (J2 <= NRBSITES/2) THEN
! If J2 is a patch, display as a fluorine atom
                  WRITE(26,'(A4,3F20.10)') 'F', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
! If J2 is an antipatch, display as a carbon atom
                  WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF
! For patches and antipatches the same, use this block
! Display patches and antipatches both as fluorine atoms
!              WRITE(26,'(A4,3F20.10)') 'F', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWPAP_LATTICE

