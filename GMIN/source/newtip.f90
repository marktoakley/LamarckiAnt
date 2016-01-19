! Calculates the potential energy and gradients for the TIPN water
! potential
! See 'Simulations of rigid bodies in an angle axis framework', Dwaipayan
! Chakrabarti and David Wales, Phys. Chem. Chem. Phys., 2009, 11, 
! 1970-1976 for an explanation of the angle axis scheme

!     X = the positions and orientations of the water molecules
!     G = the gradients of the potential energy surface for changes in
!         the positions and orientations of the water molecules
!     ENERGY = the potential energy of the current configuration
!     GTEST = true if gradients need to be calculated
      SUBROUTINE NEWTIP(X, G, ENERGY, GTEST) 

!     Global variables:
!     NATOMS = the number of lines in coords input file
!     NRBSITES = the number of sites per water molecule
!     SITE = the positions of sites in the refernce geometry
!     TIPID = identifier for the specific TIP potential
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, TIPID, VT, VTVT, WATERDEGREE

      IMPLICIT NONE

!     I, J, J1-J8 = counters to iterate over and site indices
!     OFFSET = the offset to the start of orientations in X array
      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 

!     CHARGE = the charges on each site in the water molecule
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), CHARGE(NRBSITES)
      DOUBLE PRECISION :: CVMAT(2)

!     RISITE = unused
!     RJSITE = unused
!     RI = the translational coordinates of molecule I
!     RJ = unused
!     DR = unused
!     DSS = the vector displacement between two sites
!     PI = the rotational coordinates of molecule I
!     PJ = unused
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DR(3), DSS(3), PI(3), PJ(3)

!     DSS2 = square modulus of the vector displacement between two sites
!            equivalent to r**2, r is the distance between two sites 
!     R2 = DSS2**-1, equivalent to r**-2
!     R6 = DSS2**-3, equivalent to r**-6
!     R12 = DSS2**-6, equivalent to r**-12
!     ABSR = DSS**0.5, equivalent to r
!     DVDR = the derivative of the LJ potential with respect to site-site
!            distance, divided by the distance 
!     C12 = coefficient of r**-12 term in LJ potential in
!           kJ / mol Angstrom**12
!     C6 = coefficient of r**-6 term in LJ potential in
!          kJ / mol Angstrom**6
!     CH2O = energy conversion factor for coulomb energy
      DOUBLE PRECISION :: DSS2, R2, R6, R12, ABSR, DVDR, ENERGY, C12, C6, CH2O

!     R = the position of each site
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)

!     DRk = matrix product of the derivative of the rotation matrix with
!           respect to kth component of the rotation vector and the 
!           position of each site in the reference geometry, as in
!           equation (12) of the paper
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3)

!     RMI = rotation matrix of molecule I
!     DRMIk = derivate of RMI with respect to the kth component of the
!             rotation vector
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

!     DOTIk = unused
!     DOTJk = unused
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      LOGICAL          :: GTEST

!     Set up the specific potentials for different TIPID
!     Subroutines are at the bottom of this file
      IF (TIPID == 1) THEN
         CALL DEFTIP1(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 2) THEN
         CALL DEFTIP2(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 3) THEN
         CALL DEFTIP3(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 4) THEN
         CALL DEFTIP4(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 5) THEN
         CALL DEFTIP5(CHARGE, C12, C6, CH2O)
      ENDIF

!     Initialise energy and gradients
      ENERGY          = 0.D0
      IF (GTEST) G(:) = 0.D0

      VT(:)=0.D0
      VTVT(:,:)=0.D0
      WATERDEGREE(:,:)=0.D0

!     Set offset as 3 times the number of molecules
      OFFSET          = 3*NATOMS/2

!     J1 Loop over all molecules
      DO J1 = 1, NATOMS/2

!        Set J3 as the index of the third translational coordinate of
!        molecule J1 and J5 as the index of the third rotational
!        coordinate of molecule J1
         J3 = 3*J1
         J5 = OFFSET + J3

!        Set RI as the translational coordinates of molecule J1 and PI 
!        as the rotational coordinates of molecule J1
         RI = X(J3-2:J3)
         PI = X(J5-2:J5)

!        Calculate the rotation matrix and derivatives thereof
         CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)

!        J2 loop over sites within each molecule
         DO J2 = 1, NRBSITES

!           Set J4 as the index of the J2th site of molecule J1
            J4        = NRBSITES*(J1-1) + J2

!           Calculate the position of the J2th site as the position of
!           the centre of mass of the molecule plus the rotation matrix
!           acting on the displacement from the centre of mass of the 
!           site in the reference geometry
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))

            IF (GTEST) THEN

!              If derivatives are required, calculate DRMIk.r0i, as in
!              equation (12) of the paper
               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

!     J1 loop over molecules, excluding the final molecule
      DO J1 = 1, NATOMS/2 - 1 

!        Set J3 as the index of the third translational coordinate of
!        molecule J1 and J5 as the index of the third rotational
!        coordinate of molecule J1
         J3 = 3*J1
         J5 = OFFSET + J3

!        J2 loop over molecules with index greater than J1, to perform
!        the double sum J2 > J1
         DO J2 = J1 + 1, NATOMS/2   

!        Set J4 as the index of the third translational coordinate of
!        molecule J2 and J6 as the index of the third rotational
!        coordinate of molecule J2
            J4 = 3*J2
            J6 = OFFSET + J4

!           Oxygen-Oxygen LJ contribution to the energy
!           Set J7 as the index of the oxygen site of molecule J1 and J8
!           as the index of the oxygen site of molecule J2
            J7 = NRBSITES*(J1-1) + 1
            J8 = NRBSITES*(J2-1) + 1

!           Calculate DSS as the vector displacement between the two 
!           oxygens and take the square modulus
            DSS(:) = R(J7,:) - R(J8,:)
            DSS2   = DOT_PRODUCT(DSS,DSS)

!           Use DSS2 to find r**-6 and r**-12 required for LJ potential
            R2     = 1.D0 / DSS2
            R6     = R2 * R2 * R2
            R12    = R6 ** 2

!           Calculate the LJ contribution to the potential energy
            ENERGY = ENERGY + C12*R12 - C6*R6
            VT(J1)=VT(J1) + ( C12*R12 - C6*R6 )
            VT(J2)=VT(J2) + ( C12*R12 - C6*R6 )
            VTVT(J1,J2)=VTVT(J1,J2) + ( C12*R12 - C6*R6 )
            VTVT(J2,J1)=VTVT(J1,J2)

            IF (GTEST) THEN

!              Need to calculate derivatives
!              Find derivative of LJ potential with respect to change in
!              site-site distance
!              Note DVDR is the derivative divided by the distance
               DVDR       = -6.D0*(2.D0*C12*R12 - C6*R6)*R2

!              Find the gradients for a change in position of molecule
!              J1 due to molecule J2, and vice-versa, using equations
!              (10) and (11) in the paper
               G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

!              Find the gradients for a change in orientation of molecule
!              J1 due to molecule J2, and vice-versa, using equations
!              (10) and (12) in the paper  
               G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
               G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
               G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

               G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
               G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
               G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

            ENDIF

            CVMAT=0.D0

!           I loop over each site in molecule J1
            DO I = 1, NRBSITES

!              Set J7 as index of site I in molecule J1
               J7 = NRBSITES*(J1-1) + I

!              J loop over each site in molecule J2
               DO J = 1, NRBSITES

!                 Set J8 as index of site J in molecule J2
                  J8 = NRBSITES*(J2-1) + J 

!                 Calculate DSS as the vector displacement between the
!                 two sites, and take the square modulus
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DSS(1)*DSS(1) + DSS(2)*DSS(2) + DSS(3)*DSS(3)

!                 Use DSS2 to find r**-2 and r**-1 as required for
!                 Coulomb potential and its derivatives
                  R2     = 1.D0 / DSS2
                  ABSR   = SQRT(DSS2)

!                 Calculate the Coulomb contribution to the potential
!                 energy
                  ENERGY = ENERGY + CH2O*CHARGE(I)*CHARGE(J)/ABSR
                  VT(J1)=VT(J1) + ( CH2O*CHARGE(I)*CHARGE(J)/ABSR )
                  VT(J2)=VT(J2) + ( CH2O*CHARGE(I)*CHARGE(J)/ABSR )
                  VTVT(J1,J2)=VTVT(J1,J2) + ( CH2O*CHARGE(I)*CHARGE(J)/ABSR )
                  VTVT(J2,J1)=VTVT(J2,J1) + ( CH2O*CHARGE(I)*CHARGE(J)/ABSR )

                  IF (I==4) CVMAT(1)=CVMAT(1)+CH2O*CHARGE(I)*CHARGE(J)/ABSR
                  IF (J==4) CVMAT(2)=CVMAT(2)+CH2O*CHARGE(I)*CHARGE(J)/ABSR

                  IF (GTEST) THEN

!                    Need to calculate derivatives
!                    Find derivative of LJ potential with respect to
!                    change in site-site distance
!                    Note DVDR is the derivative divided by the distance
                     DVDR       = -CH2O*CHARGE(I)*CHARGE(J)*R2/ABSR 

!                    Find the gradients for a change in position of
!                    molecule J1 due to molecule J2, and vice-versa,
!                    using equations (10) and (11) in the paper
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:) 


!                    Find the gradients for a change in orientation of
!                    molecule J1 due to molecule J2, and vice-versa,
!                    using equations (10) and (12) in the paper  
                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

            IF (CVMAT(1)<CVMAT(2)) THEN
               VTVT(J2,J1)=0.D0
               IF (VTVT(J1,J2)<-12.D0) THEN
                  WATERDEGREE(:,J1)=WATERDEGREE(:,J1)+(/-1, 1 /)
                  WATERDEGREE(:,J2)=WATERDEGREE(:,J2)+(/ 1, 1 /)
               ENDIF
            ENDIF
            IF (CVMAT(2)<CVMAT(1)) THEN
               VTVT(J1,J2)=0.D0
               IF (VTVT(J2,J1)<-12.D0) THEN
                  WATERDEGREE(:,J2)=WATERDEGREE(:,J2)+(/-1, 1 /)
                  WATERDEGREE(:,J1)=WATERDEGREE(:,J1)+(/ 1, 1 /)
               ENDIF
            ENDIF

         ENDDO

      ENDDO

      

      END SUBROUTINE NEWTIP 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP1(CHARGE, C12, C6, CH2O)
!     TIPS water    

      USE COMMONS, ONLY: NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      CHARGE(:) = (/-0.8D0, 0.4D0, 0.4D0/)
      C6        = 2510.4D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2426720.D0
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP1

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP2(CHARGE, C12, C6, CH2O)
!     TIPS2 water    
    
      USE COMMONS, ONLY: NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      SITE(4,1) = 0.D0
      SITE(4,2) = 0.D0
      SITE(4,3) = ROM

      CHARGE(:) = (/0.D0, 0.535D0, 0.535D0, -1.07D0/)
      C6        = 2510.4D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2907880.D0
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP2

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP3(CHARGE, C12, C6, CH2O)
!     TIP3P water    
    
      USE COMMONS, ONLY: NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      CHARGE(:) = (/-0.834D0, 0.417D0, 0.417D0/)
      C6        = 2489.48D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2435088.D0
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP3

!     ----------------------------------------------------------------------------------------------

!     Sets up the TIP4P potential
!     CHARGE = the charge on each site
!     C12 = coefficient of r**-12 term in LJ potential in
!           kJ / mol Angstrom**12
!     C6 = coefficient of r**-6 term in LJ potential in
!          kJ / mol Angstrom**6
!     CH2O = energy conversion factor for coulomb energy
      SUBROUTINE DEFTIP4(CHARGE, C12, C6, CH2O)
    
!     NATOMS = the number of postion and orientation coordinates
!     NRBSITES = the number of sites per molecule
!     SITE = the positions of the sites in the refernce geometry
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, VTVT, WATERDEGREE
      
      IMPLICIT NONE

!     I = counter to iterate over
      INTEGER          :: I

!     M = mass of each site
!     MASS = the mass of a water molecule
!     CM(3) = position of the centre of mass relative to the oxygen
!             position
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O

!     THETA = the H-O-H bond angle
!     ROH = the O-H bond length
!     ROM = displacement of the 4th site from oxygen along the z axis
!     PI = the mathematical constant pi
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      IF (.NOT.ALLOCATED(VTVT)) ALLOCATE(VTVT(NATOMS/2,NATOMS/2))
      IF (.NOT.ALLOCATED(WATERDEGREE)) ALLOCATE(WATERDEGREE(2,NATOMS/2))

!     Define the geometry of the water molecule
!     Convert THETA to radians
      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     Set the positions of the sites relative to oxygen
!     The reference geometry is on the y-z plane

!     Oxygen at the origin
      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

!     Hydrogen 1
      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

!     Hydrogen 2
      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH

!     4th site
      SITE(4,1) = 0.D0
      SITE(4,2) = 0.D0
      SITE(4,3) = ROM

!     Set the charges for O, H1, H2, and the 4th site
      CHARGE(:) = (/0.D0, 0.52D0, 0.52D0, -1.04D0/)

!     Set the LJ coeffeicients
      C6        = 2552.24D0
      C12       = 2510.4D3
      CH2O      = 1389.354848D0

!     Set the mass of each site
      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)

!     Calculate the position of the centre of mass relative to oxygen
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS

!     Convert site coordinates from oxygen origin to centre of mass 
!     origin
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP4

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP5(CHARGE, C12, C6, CH2O)
!     TIP5P water    
    
      USE COMMONS, ONLY: NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETAH, THETAM, ROH, ROM, PI

      PI     = 4.D0*DATAN(1.D0)
      ROH    = 0.9572D0
      ROM    = 0.70D0
      THETAH = 104.52D0
      THETAM = 109.47D0
      THETAH = PI*THETAH/180.D0
      THETAM = PI*THETAM/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETAH)*ROH
      SITE(2,3) = COS(0.5D0*THETAH)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) =-SIN(0.5D0*THETAH)*ROH
      SITE(3,3) = COS(0.5D0*THETAH)*ROH
  
      SITE(4,1) = SIN(0.5D0*THETAM)*ROM 
      SITE(4,2) = 0.D0
      SITE(4,3) =-COS(0.5D0*THETAM)*ROM 

      SITE(5,1) =-SIN(0.5D0*THETAM)*ROM 
      SITE(5,2) = 0.D0
      SITE(5,3) =-COS(0.5D0*THETAM)*ROM
 
      CHARGE(:) = (/0.D0, 0.241D0, 0.241D0, -0.241D0, -0.241D0/)
      C6        = 2470.012857D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2278383.244D0
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP5

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWNEWTIP

      USE COMMONS, ONLY: NATOMS, SITE, NSAVE, WATERDEGREE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7, GETUNIT, NEWTIPUNIT
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      NEWTIPUNIT=GETUNIT()
      OPEN(UNIT=NEWTIPUNIT, FILE='newtip.xyz', STATUS='UNKNOWN')

      GTEST = .FALSE. 

      DO J1 = 1, NSAVE

         WRITE(NEWTIPUNIT,'(I6)') (NATOMS/2)*3
         WRITE(NEWTIPUNIT,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

         DO J3 = 1, NATOMS/2

            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            DO J2 = 1, 3

               RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))

               IF (J2 == 1) THEN
                  WRITE(NEWTIPUNIT,'(A1,I1,2X,3F20.10)') 'O', WATERDEGREE(1,J3)+1, RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  WRITE(NEWTIPUNIT,'(A1,I1,2X,3F20.10)') 'H', WATERDEGREE(1,J3)+1, RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

         ENDDO

      ENDDO

      CLOSE(NEWTIPUNIT)

      END SUBROUTINE VIEWNEWTIP

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RIGIDTOTIP(N,COMCOORDS,AACOORDS,ATOMCOORDS,FUNIT)
      USE COMMONS, ONLY: SITE, WATERDEGREE

      DOUBLE PRECISION, INTENT(IN)  :: COMCOORDS(3),AACOORDS(3)
      DOUBLE PRECISION, INTENT(OUT) :: ATOMCOORDS(3,3)
      INTEGER :: J1, N, NEWTIPUNIT
      INTEGER :: FUNIT
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)

!     IF (PRESENT(FUNIT)) THEN
         NEWTIPUNIT=FUNIT
!     ELSE
!        NEWTIPUNIT=-1
!     ENDIF

      CALL RMDRVT(AACOORDS, RMI, DRMI, DRMI, DRMI, .FALSE.)

      DO J1 = 1, 3

         ATOMCOORDS(:,J1) = COMCOORDS(1:3) + MATMUL(RMI(:,:),SITE(J1,:))

         IF (NEWTIPUNIT>0) THEN
            IF (J1==1) THEN
               WRITE(NEWTIPUNIT,'(A1,I1,2X,3F20.10)') 'O', WATERDEGREE(1,N)+2, ATOMCOORDS(:,J1)
            ELSE
               WRITE(NEWTIPUNIT,'(A1,I1,2X,3F20.10)') 'H', WATERDEGREE(1,N)+2, ATOMCOORDS(:,J1)
            ENDIF
         ENDIF

      ENDDO

      END SUBROUTINE RIGIDTOTIP 
