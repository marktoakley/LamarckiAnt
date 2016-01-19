!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAP                                                                                          !
!                                                                                              !
! Calculates Energy and gradients for the patch-antipatch potential, using the rigid body      !
! angle axis framework described in:                                                           !
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and          !
! David Wales, Phys. Chem. Chem. Phys., 2009, 11, 1970-1976                                    !
! X: the positions and orientations of the bodies                                              !
! G: the gradients of the potential energy surface for changes in positions and orientations   !
! ENERGY: the potential energy for the configuration stored in X                               !
! GTEST: logical, true if gradients need to be calculated                                      !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAP(X, G, ENERGY, GTEST)

! NATOMS: twice the number of bodies, one for position and one for orientation
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the directions of the patch and antipatch vectors
! PERIODIC: true if periodic boundary conditions are to be used
! BOXLX, BOXLY, BOXLZ: dimensions of box for periodic boundary conditions
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, PERIODIC, BOXLX, BOXLY, BOXLZ

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
      DOUBLE PRECISION :: E(NATOMS*NRBSITES/2,3), P(3)

! DEk: matrix product of the derivative of the rotation matrix with respect to kth component of the rotation vector and the
!      position of each site in the reference geometry, as in equation (12) of the paper
      DOUBLE PRECISION :: DE1(NATOMS*NRBSITES/2,3), DE2(NATOMS*NRBSITES/2,3), DE3(NATOMS*NRBSITES/2,3) 

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

! XI, XJ: coordinates of bodies, used for periodic testing

! Initialise energy and gradients
      ENERGY  = 0.D0
      IF (GTEST) THEN
        G(:) = 0.D0
      ENDIF

! Find NMOL as the actual number of bodies, and OFFSET as the start of rotational coordinates
      NMOL    = NATOMS/2
      OFFSET  = 3*NMOL

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

      IF (.NOT. PERIODIC) THEN
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

! Calculate LJ energy
            CALL PAPENERGYLJ(IJENERGY, GIT, GJT, GTEST, X(IT-2:IT), X(JT-2:JT))
            ENERGY = ENERGY + IJENERGY
            IF (GTEST) THEN
              G(IT-2:IT) = G(IT-2:IT) + GIT(:)
              G(JT-2:JT) = G(JT-2:JT) + GJT(:)
            ENDIF

! Set S and T as appropriate indices of patches in bodies I and J
            S = (I-1)*NRBSITES
            T = (J-1)*NRBSITES

! Calculate PAP energy
            CALL PAPENERGYPAP(IJENERGY, GIT, GJT, GIR, GJR, GTEST, X(IT-2:IT), X(JT-2:JT),         &
     &                        E(S+1:S+NRBSITES,:), E(T+1:T+NRBSITES,:), DE1(S+1:S+NRBSITES,:),     &
     &                        DE2(S+1:S+NRBSITES,:), DE3(S+1:S+NRBSITES,:), DE1(T+1:T+NRBSITES,:), &
     &                        DE2(T+1:T+NRBSITES,:), DE3(T+1:T+NRBSITES,:))
! Add on energy contribution
            ENERGY = ENERGY + IJENERGY
            IF (GTEST) THEN
! Add on gradient contributions from LJ and PAP
              G(IT-2:IT) = G(IT-2:IT) + GIT(:)
              G(JT-2:JT) = G(JT-2:JT) + GJT(:)
              G(IR-2:IR) = G(IR-2:IR) + GIR(:)
              G(JR-2:JR) = G(JR-2:JR) + GJR(:)
            ENDIF ! End if GTEST
          ENDDO ! End loop over particles J
        ENDDO ! End loop over particles I

      ELSE ! Periodic Boundary Conditions

! Calculate the cutoff distance of the potential
        RCUT = BOXLX
        IF (BOXLY .LT. RCUT) RCUT = BOXLY
        IF (BOXLZ .LT. RCUT) RCUT = BOXLZ
        RCUTSQ = (RCUT/2.D0)**2

! Make sure all coordinates lie within the box (which extends from -BOXLk/2 to +BOXLk/2 in all 3 directions
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
                ENDIF
                CALL PAPENERGYPAP(IJENERGY, GIT, GJT, GIR, GJR, GTEST, XI(:), X(JT-2:JT),              &
     &                            E(S+1:S+NRBSITES,:), E(T+1:T+NRBSITES,:), DE1(S+1:S+NRBSITES,:),     &
     &                            DE2(S+1:S+NRBSITES,:), DE3(S+1:S+NRBSITES,:), DE1(T+1:T+NRBSITES,:), &
     &                            DE2(T+1:T+NRBSITES,:), DE3(T+1:T+NRBSITES,:))
                ENERGY = ENERGY + IJENERGY
                IF (GTEST) THEN
                  G(IT-2:IT) = G(IT-2:IT) + GIT(:)
                  G(JT-2:JT) = G(JT-2:JT) + GJT(:)
                  G(IR-2:IR) = G(IR-2:IR) + GIR(:)
                  G(JR-2:JR) = G(JR-2:JR) + GJR(:)
                ENDIF 

! No point continuing to loop over coordinate configurations now that we have found the interacting one
                EXIT

              ENDIF ! End if found best configuration
            ENDDO ! End loop ove coordinate configurations
          ENDDO ! End loop over particles J 
        ENDDO ! End loop over particles I
      ENDIF ! End if PERIODIC

      END SUBROUTINE PAP
 
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! DEFPAP                                                                                       !
!                                                                                              !
! Defines pap potential by specifying the positions of patches and antipatches                 !
! The first half of the vectors are patches, the second half anitpatches                       !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE DEFPAP()

! MYUNIT: file handle for GMIN_out
! NRBSITES: the total number of patches and antipatches per body
! RBSTLA: vectors specifying the directions of patches and antipatches
!         these vecotrs MUST be NORMALISED
      USE COMMONS, ONLY: MYUNIT, NRBSITES, RBSTLA

      IMPLICIT NONE

! Loop counter
      INTEGER :: I
      DOUBLE PRECISION :: LENGTH

      IF (NRBSITES == 0) THEN
! Read from input file papsites.xyz
! File format of papsites.xyz:
! [Number of sites]
! [blank line]
! [X patch1] [Y patch1] [Z patch1]
! ...
! [X antipatch1] [Y anitpatch2] [Z antipatch2]
! ...
!
! The first half of the sites will be interpreted as patches and the second half as antipatches

! Open file and read in number of sites
        OPEN(UNIT=394,FILE="papsites.xyz",STATUS="old")
        READ(394,*) NRBSITES 
        READ(394,*)

! Check NRBSITES has been read correctly
        IF (NRBSITES .EQ. 0) THEN
          WRITE(MYUNIT, *) 'DEFPAP> ERROR: NRBSITES not read correctly from papsites.xyz'
          STOP
        ELSE IF (MOD(NRBSITES, 2) .NE. 0) THEN
          WRITE(MYUNIT, *) 'DEFPAP> ERROR: NRBSITES must be even'
          STOP
        ENDIF

! Allocate memory for the sites
        ALLOCATE(RBSTLA(NRBSITES,3))

! Loop over the number of sites
        DO I = 1, NRBSITES
          READ(394, *) RBSTLA(I, 1), RBSTLA(I, 2), RBSTLA(I, 3)
          
! Now normalise the vector
          LENGTH = SQRT(RBSTLA(I, 1)**2 + RBSTLA(I, 2)**2 + RBSTLA(I, 3)**2)
          RBSTLA(I,:) = RBSTLA(I,:)/LENGTH
          
        ENDDO
! Close papsites.xyz
        CLOSE(394)

      ELSE IF (NRBSITES == 2) THEN
! Create a patch along the positive z axis and an antipatch along the negative z axis
        ALLOCATE(RBSTLA(NRBSITES,3))
        RBSTLA(1,:) = (/ 0.D0, 0.D0, 1.D0/)
        RBSTLA(2,:) = (/ 0.D0, 0.D0,-1.D0/)

      ELSE IF (NRBSITES == 4) THEN
! Create a tetrahedral arrangement of patches and antipatches
        ALLOCATE(RBSTLA(NRBSITES,3))
!        RBSTLA(1,:)= 1.D0/SQRT(3.D0)*(/  1.D0,  1.D0,  1.D0/)
!        RBSTLA(2,:)= 1.D0/SQRT(3.D0)*(/ -1.D0, -1.D0,  1.D0/)
!        RBSTLA(3,:)= 1.D0/SQRT(3.D0)*(/ -1.D0,  1.D0, -1.D0/)
!        RBSTLA(4,:)= 1.D0/SQRT(3.D0)*(/  1.D0, -1.D0, -1.D0/)
! Bernal patches 1 and 2 only
        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSTLA(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSTLA(3,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSTLA(4,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)

      ELSE IF (NRBSITES == 6) THEN
! Trying to make Bernel spiral
        ALLOCATE(RBSTLA(NRBSITES,3))
        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSTLA(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSTLA(3,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
        RBSTLA(4,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSTLA(5,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
        RBSTLA(6,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE
        WRITE(MYUNIT,*) 'Number of patches must be 2 (linear), 4 (tetrahedral), ', &
     &                  '6 (Bernal spiral) or 0 (read from papsites.xyz)'

      ENDIF

      END SUBROUTINE DEFPAP

!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAPENERGYLJ                                                                                  !
!                                                                                              !
! Calculates LJ energy interaction and gradients between I and J                               !
! IJENERGY: the energy contribution due to interaction between bodies I and J                  !
! GIT: the gradients for translational coordinates of body I                                   !
! GJT: the gradients for translational coordinates of body J                                   !
! GTEST: true if gradients are required                                                        !
! XI: position of body I                                                                       !
! XJ: position of body J                                                                       !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAPENERGYLJ(IJENERGY, GIT, GJT, GTEST, XI, XJ)

! PAPALP: PAP alpha parameter
      USE COMMONS, ONLY: PAPALP

      IMPLICIT NONE

      DOUBLE PRECISION :: IJENERGY, GIT(3), GJT(3)
      DOUBLE PRECISION :: XI(3), XJ(3)
      LOGICAL          :: GTEST


! RIJ: vector displacement between I and J
! LJSIGMASQ: the sigma value in the LJ potential squared
! LJN: the exponent in the LJ potential
! DSS2: distance between I and J squared
! R2: 1/DSS2, relevant to LJ potential
! RLJN, RLJ2N: terms in LJ potntial
! DVDR: derivative of potential wrt distance, divided by the distance
      DOUBLE PRECISION :: RIJ(3), LJSIGMASQ, LJN, DSS2, R2, RLJN, R2LJN, DVDR

! EXP6: R2D**6, denominator in LJ potential
! EXP3: R2D**3, denomiator in LJ potential
! R2D: 1/(DSS2-1), relevant to LJ potential
!      DOUBLE PRECISION :: EXP3, EXP6, R2D, PREFAC

! Initialise variables
      IJENERGY = 0.D0
      GIT(:) = 0.D0
      GJT(:) = 0.D0 
!      PREFAC  = 4.D0/(PAPALP*PAPALP)
      RIJ(:) = XI(:) - XJ(:)
      DSS2   = DOT_PRODUCT(RIJ(:),RIJ(:))

! Calculate denomiators relevant to the LJ potential
! Add the LJ contribution to the potential energy from the interaction of bodies J1 and J2
! LJ potential from proposed PAP potential description, generates structures with overlapping cores
!      R2D    = 1.D0/(DSS2 - 1.D0)
!      EXP6   = R2D**6
!      EXP3   = R2D**3
!      IJENERGY = PREFAC*(EXP6 - PAPALP*EXP3)

! LJ 2n-n potential, with epsilon equal to one
! Sigma must be adjusted to match the proposed potential above
! LJN chosen to give a reasonable match to the above potential
      R2 = 1.D0/DSS2
      LJN = 23
      LJSIGMASQ = 1.D0 + 1.D0/PAPALP**(1.D0/3.D0)
      RLJN = (LJSIGMASQ*R2)**(LJN/2.D0)
      R2LJN = RLJN**2
      IJENERGY = IJENERGY + 4.D0*(R2LJN-RLJN)

      IF (GTEST) THEN
! LJ potential from proposed PAP potential description, generates structures with overlapping cores
!        DVDR = 2.D0*PREFAC*(-6*EXP6 + 3.D0*PAPALP*EXP3)*R2D

! LJ 2n-n potential, with epsilon equal to one
        DVDR = 4.D0*LJN*R2*(RLJN-2*R2LJN)
        GIT(:) = GIT(:) + DVDR*RIJ(:)
        GJT(:) = GJT(:) - DVDR*RIJ(:)
!        WRITE(*,*) 'RIJ: ', RIJ(1), RIJ(2), RIJ(3)
!        WRITE(*,*) 'DSS2, R2: ', DSS2, R2
!        WRITE(*,*) 'LJN, LJSIGMASQ: ', LJN, LJSIGMASQ
!        WRITE(*,*) 'RLJN, R2LJN: ', RLJN, R2LJN
!        WRITE(*,*) 'GIT, GJT: ', GIT(1), GIT(2), GIT(3), GJT(1), GJT(2), GJT(3)
      ENDIF

      END SUBROUTINE PAPENERGYLJ

!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAPLJPERIODIC                                                                                !
!                                                                                              !
! Calculates LJ energy and gradients between I and J                                           !
! Uses a modified form of the LJ potential which goes precisely and continuously to zero at    !
! a cutoff distance, which should be at least 3*LJSIGMA for good agreement with normal LJ.     !
! See: 'Numerical Experiments on the Stochastic Behaviour of a Lennard-Jones Gas System',      !
! Spotswood D. Stoddard and Joseph Ford, Physical Review A, 1973, 8, 1504-1512                 !
! Particularly equation (2) and figure 1.                                                      !
! Calculates LJ energy interaction and gradients between I and J                               !
! IJENERGY: the energy contribution due to interaction between bodies I and J                  !
! GIT: the gradients for translational coordinates of body I                                   !
! GJT: the gradients for translational coordinates of body J                                   !
! GTEST: true if gradients are required                                                        !
! RIJ: vector displacement between I and J                                                !
! RCUTSQ: cutoff distance for the potential                                                    !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAPLJPERIODIC(IJENERGY, GIT, GJT, GTEST, RIJ, RCUTSQ)

! PAPALP: PAP alpha parameter
      USE COMMONS, ONLY: PAPALP

      IMPLICIT NONE

      DOUBLE PRECISION :: IJENERGY, GIT(3), GJT(3), RIJ(3), RCUTSQ
      LOGICAL          :: GTEST

! LJSIGMASQ: the sigma value in the LJ potential squared
! LJN: the exponent in the LJ potential
! DSS2: squared distance between I and J
! R2: 1/DSS2, relevant to LJ potential
! RLJN, RLJ2N: terms in LJ potntial
! DVDR: derivative of potential wrt distance, divided by the distance
! RCUTSQINV: 1/RCUTSQ
      DOUBLE PRECISION :: LJSIGMASQ, LJN, DSS2, R2, RLJN, R2LJN, DVDR, RCUTSQINV

! RCLJN, RC2LJN: cutoff terms in LJ potential
      DOUBLE PRECISION :: RCLJN, RC2LJN

! Initialise variables
      IJENERGY = 0.D0
      GIT(:) = 0.D0
      GJT(:) = 0.D0 
      DSS2 = DOT_PRODUCT(RIJ(:),RIJ(:))

! LJ 2n-n potential, with epsilon equal to one
! Sigma must be adjusted to match the normal LJ potential above
! LJN chosen to give a reasonable match to the above potential
      RCUTSQINV = 1.D0/RCUTSQ
      R2 = 1.D0/DSS2
      LJN = 23
      LJSIGMASQ = 1.D0 + 1.D0/PAPALP**(1.D0/3.D0)
      RLJN = (LJSIGMASQ*R2)**(LJN/2.D0)
      R2LJN = RLJN**2.D0
      RCLJN = (LJSIGMASQ*RCUTSQINV)**(LJN/2.D0) 
      RC2LJN = (LJSIGMASQ*RCUTSQINV)**(LJN)
      IJENERGY = IJENERGY + 4.D0*(R2LJN-RLJN + (LJN*RC2LJN - (LJN/2.D0)*RCLJN)*(DSS2*RCUTSQINV) &
     &                            - (LJN+1)*RC2LJN + (LJN/2.D0 + 1)*RCLJN)

      IF (GTEST) THEN

! LJ 2n-n potential, with epsilon equal to one
        DVDR = 4.D0*LJN*(R2*(RLJN-2*R2LJN) + RCUTSQINV*(2*RC2LJN-RCLJN))
        GIT(:) = GIT(:) + DVDR*RIJ(:)
        GJT(:) = GJT(:) - DVDR*RIJ(:)
      ENDIF

      END SUBROUTINE PAPLJPERIODIC

!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAPENERGYPAP                                                                                 !
!                                                                                              !
! Calculates the energy interaction between I and J for the pap interaction                    !
! IJENERGY: energy interaction between bodies I and J                                          !
! GIT: the gradients for translational coordinates of body I                                   !
! GJT: the gradients for translational coordinates of body J                                   !
! GIR: the gradients for rotational coordinates of body I                                      !
! GJR: the gradients for rotational coordinates of body J                                      !
! GTEST: true if gradients are required                                                        !
! XI: position of body I                                                                       !
! XJ: position of body J                                                                       !
! EI(NRBSITES:3): E vectors for body I (see PAP())                                             !
! EJ(NRBSITES:3): E vectors for body J                                                         !
! DEkI: derivatives of EI                                                                      !
! DEkJ: derivatives of EJ                                                                      !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAPENERGYPAP(IJENERGY, GIT, GJT, GIR, GJR, GTEST, XI, XJ, EI, EJ, &
     &                        DE1I, DE2I, DE3I, DE1J, DE2J, DE3J) 

! NRBSITES: number of patches plus anti-patches per body
! PAPALP: PAP alpha parameter
! PAPCD: PAP cos delta parameter
! PAPEPS: PAP epsilon parameter
! PAPS: PAP s parameter
      USE COMMONS, ONLY: NRBSITES, PAPALP, PAPCD, PAPEPS, PAPS

      IMPLICIT NONE

      DOUBLE PRECISION :: IJENERGY, GIT(3), GJT(3), GIR(3), GJR(3)
      DOUBLE PRECISION :: XI(3), XJ(3), EI(NRBSITES,3), EJ(NRBSITES,3)
      DOUBLE PRECISION :: DE1I(NRBSITES,3), DE2I(NRBSITES,3), DE3I(NRBSITES,3)
      DOUBLE PRECISION :: DE1J(NRBSITES,3), DE2J(NRBSITES,3), DE3J(NRBSITES,3)
      LOGICAL :: GTEST

! ANGFAC: factor in the angular potential phi, PI/(1-PAPCD)
! INVS: the reciprocal of the pap s parameter
! LAMBDA: longest distance with full PAP attraction
! PI: the mathematical constant
      DOUBLE PRECISION :: ANGFAC, INVS, LAMBDA, PI

! RIJ: vector displacement between I and J
! NR: normalised RIJ
! DSS2: distance between I and J squared
! ABSR: distance between I and J
! DELR: distance parameter for PAP interaction
      DOUBLE PRECISION :: RIJ(3), NR(3), DSS2, ABSR, DELR

! S: patch/anti-patch in body I
! T: patch/anti-patch in body J
! TLOW: lower bound on patch index for looping on body J
! TUP: upper bound on patch index for looping on body J
      INTEGER :: S, T, TLOW, TUP 

! WP: distance part of PAP potential
! DWPDR: derivative of WP wrt distance
! DOTS: angle between orientation of patch S and the vector displacement between two bodies
! DOTT: angle between orientation of patch T and the vector displacement between two bodies
! ARGS: argument for PHIS in the angular potential
! ARGT: argument for PHIT in the angular potential
! PHIS: function in the angular potential
! PHIT: function in the angular potential
      DOUBLE PRECISION :: WP, DWPDR, DOTS, DOTT, ARGS, ARGT, PHIS, PHIT

! DDOTSDR:  derivative of DOTS wrt a change in distance bewteen I and J
! DDOTTDR:  derivative of DOTT wrt a change in distance bewteen I and J
! DPHISDR:  derivative of PHIS wrt a change in distance bewteen I and J
! DPHITDR:  derivative of PHIT wrt a change in distance bewteen I and J
! DPHISDPI: derivative of PHIS wrt a change in orientation of body I (note DPHIIDPJ would be zero)
! PPHITDPJ: derivative of PHIT wrt a change in orientation of body J (note DPHIJDPI would be zero)
      DOUBLE PRECISION :: DDOTSDR(3), DDOTTDR(3), DPHISDR(3), DPHITDR(3), DPHISDPI(3), DPHITDPJ(3)

! Initialise variables
      IJENERGY = 0.D0
      GIT(:)   = 0.D0
      GJT(:)   = 0.D0
      GIR(:)   = 0.D0
      GJR(:)   = 0.D0
      PI       = 4.D0*DATAN(1.D0)
      ANGFAC   = PI/(1.D0 - PAPCD)
      INVS     = 1.D0/PAPS
      LAMBDA   = SQRT(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))

! Find normalised vector between I and J, and distance parameter DELR
      RIJ(:) = XI(:) - XJ(:)
      DSS2   = DOT_PRODUCT(RIJ(:),RIJ(:))
      ABSR   = SQRT(DSS2)
      NR(:)  = RIJ(:)/ABSR
      DELR   = ABSR - LAMBDA

! Set specific cases for the patch-antipatch potential and its gradient
      IF (DELR > INVS) THEN
! Distance is greater than LAMBDA+INVS, so no attraction
! No point proceeding with the calculation
        RETURN
      ELSE IF (DELR < 0.D0) THEN
! Distance less than LAMBDA, so full attraction
        WP =-1.D0 
        DWPDR = 0.D0           
      ELSE
! Distance is between LAMBDA and LAMBDA+INVS, so potential varies from -1 to 0 
        WP =-0.5D0*(1.D0 + COS(PI*DELR*PAPS))
        DWPDR = 0.5D0*PI*PAPS*SIN(PI*DELR*PAPS) 
      ENDIF

! S: Loop over patches in body I
      DO S = 1, NRBSITES

! jwrm2> set bounds on T for patches only, use this for patches only
!        TLOW = 1
!        TUP  = NRBSITES
! jwrm2> set bounds on T for patch-antipatch, use this for patch-antipatch
        IF (S .LE. NRBSITES/2) THEN
          TLOW = NRBSITES/2 + 1
          TUP  = NRBSITES
        ELSE
          TLOW = 1
          TUP = NRBSITES/2
        ENDIF
! jwrm2> select a specific patch, use this for patch-antipatch pairs
!         IF (S .LE.NRBSITES/2) THEN
!           TLOW = S + NRBSITES/2
!         ELSE
!           TLOW = S - NRBSITES/2
!         ENDIF
!         TUP = TLOW
! jwrm2> end comments for patch only/patch-antipatch

! T Loop over patches in body J
        DO T = TLOW, TUP

! Find the angles between the orientation of (anti)patch S and NR, and between the orientation of (anti)patch T and NR
! Negative for DOTI, since RIJ is displacement from J2 to J1
          DOTS =-DOT_PRODUCT(EI(S,:),NR(:))
          DOTT = DOT_PRODUCT(EJ(T,:),NR(:))

! Calculate the values of PHI for the angular potential
          IF (DOTS < PAPCD) THEN
! Angle greater than patch width, so no attraction
            PHIS =-1.D0
          ELSE
! Angle less than patch width, so some attraction
            ARGS = ANGFAC*(DOTS-PAPCD)
            PHIS = -COS(ARGS) ! jwrm2> changed to minus
          ENDIF 
          IF (DOTT < PAPCD) THEN
! Angle greater than patch width, so no attraction
            PHIT =-1.D0
          ELSE
! Angle less than patch width, so some attraction
            ARGT = ANGFAC*(DOTT-PAPCD)
            PHIT = -COS(ARGT) ! jwrm2> changed to minus
          ENDIF 

! Add the patch-antipatch attraction to the potential energy

          IJENERGY = IJENERGY + 0.25D0*PAPEPS*(1.D0 + PHIS)*(1.D0 + PHIT)*WP
          IF (GTEST) THEN
! Need to find the derivatives of the patch-antipatch attraction wrt translational and rotational coordinates
            IF ((DOTS >= PAPCD) .AND. (DOTT >= PAPCD)) THEN
! Calculate the derivates of DOTI and DOTJ wrt a change in distance between I and J
              DDOTSDR(:) = -DOTS*RIJ(:)/DSS2 - EI(S,:)/ABSR
              DDOTTDR(:) = -DOTT*RIJ(:)/DSS2 + EJ(T,:)/ABSR
! Find the derivatives of the PHIs wrt a change in distance between I and J
              DPHISDR(:) = ANGFAC*SIN(ARGS)*DDOTSDR(:) ! changed to plus
              DPHITDR(:) = ANGFAC*SIN(ARGT)*DDOTTDR(:) ! changed to plus

! Add the contribution to the gradient of a change in position of I due to J, and vice versa
              GIT(:) = GIT(:) + 0.25D0*PAPEPS*((1.D0 + PHIT)*WP*DPHISDR(:)                  &
     &                 + (1.D0 + PHIS)*WP*DPHITDR(:) +(1.D0+PHIS)*(1.D0+PHIT)*DWPDR*NR(:))
              GJT(:) = GJT(:) - 0.25D0*PAPEPS*((1.D0 + PHIT)*WP*DPHISDR(:)                  &
     &                 + (1.D0 + PHIS)*WP*DPHITDR(:) +(1.D0+PHIS)*(1.D0+PHIT)*DWPDR*NR(:))
            ENDIF

            IF (DOTS >= PAPCD) THEN
! Find the derivatives of PHII wrt a change in each of the rotational coordinates of body I
              DPHISDPI(1) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE1I(S,:)) ! changed to minus
              DPHISDPI(2) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE2I(S,:)) ! changed to minus
              DPHISDPI(3) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE3I(S,:)) ! changed to minus
! Add the contribution to the gradient of a change in orientation of J1 due to J2
! Simpler because WP does not depend on orientation, and PHIJ does not depend on the orientation of J1
              GIR(:) = GIR(:) + 0.25D0*PAPEPS*(1.D0 + PHIT)*WP*DPHISDPI(:)
            ENDIF

            IF (DOTT >= PAPCD) THEN
! Find the derivatives of PHIJ wrt a change in each of the rotational coordinates of body J
              DPHITDPJ(1) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE1J(T,:)) ! changed to plus
              DPHITDPJ(2) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE2J(T,:)) ! changed to plus
              DPHITDPJ(3) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE3J(T,:)) ! changed to plus
! Add the contribution to the gradient of a change in orientation of J2 due to J1
! Simpler because WP does not depend on orientation, and PHII does not depend on the orientation of J2
              GJR(:) = GJR(:) + 0.25D0*PAPEPS*(1.D0 + PHIS)*WP*DPHITDPJ(:)
            ENDIF
          ENDIF !End GTEST
        ENDDO !End loop over patches in J
      ENDDO ! End loop over patches in I

!      WRITE(*,*) 'XJ: ', XJ

!      WRITE(*,*) 'PAPENERGYPAP: end'

      END SUBROUTINE PAPENERGYPAP
       
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! VIEWPAP                                                                                      !
!                                                                                              !
! Writes positions to file in a format readable by XMakeMol                                    !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE VIEWPAP()

! NATOMS: twice the number of bodies
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the patch and antipatch directions in the body frame
! NSAVE: the number of different configurations to be written
! PERIODIC: true if periodic boundary conditions are used
! BOXLk: Size of periodic box
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, NSAVE, PERIODIC, BOXLX, BOXLY, BOXLZ

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
         WRITE(26,'(I6)') (NATOMS/2)*(NRBSITES+1)
! Write the energy and step at which the configuration was first found
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

! J3 Loop over each body
         DO J3 = 1, NATOMS/2

! Set J5 as the index of the third translational coordinate of body J3, and J7 as the third rotational coordinate
            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
! Set P as the rotational coordinates
            P(:) = QMINP(J1,J7-2:J7)

! Find the rotation matrix
! Since GTEST is false, DRMI is a dummy variable
            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

! Write the position of the cores of the bodies
! Display as an Oxygen atom
! If this is the first body, and periodic boundary conditions are used, include the box
            IF (PERIODIC .AND. J3 .EQ. 1) THEN
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

      END SUBROUTINE VIEWPAP

