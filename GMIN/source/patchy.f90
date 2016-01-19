!|gd351>

SUBROUTINE PATCHYPOT (X, G, ENERGY, GTEST)
!SIGMASQ: squared parameter sigma of Lennard-Jones potential
!RANGESQ: squared range of potential
!FACTOR: controls patchwidth
!if SIGMASQ or RANGESQ are changed, parameters of smoothing functions VF1 and VF2 have to be adjusted; if FACTOR is changed, ALPHALEFT, ALPHARIGHT have to be adjusted
  USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, SIGMASQ, RANGESQ, FACTOR

  IMPLICIT NONE
  
  INTEGER          :: I, I1, I2, I3, J, J1, J2, J3, J4, J5, J6, J7, J8, J3T, J4T, REALNATOMS, OFFSET 
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: ENERGY, R, RSQ, RCUB, R2, R6, RHALF, ARG1, ARG2, ALPHA1, ALPHA2, VLJ, VEXP1, VEXP2, VEXP, V, DVLJ, A0, A1, A2
  DOUBLE PRECISION :: DARG1(6), DALPHA1(6), DARG2(6), DALPHA2(6), DVDX(3)
  DOUBLE PRECISION :: RI(NATOMS/2,3), RIJ(3), P(NATOMS/2,3), A(3) 
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: PATCHPOS(NATOMS/2,NRBSITES,3), DPATCHPOS1(NATOMS/2,NRBSITES,3)
  DOUBLE PRECISION :: DPATCHPOS2(NATOMS/2,NRBSITES,3), DPATCHPOS3(NATOMS/2,NRBSITES,3)
  DOUBLE PRECISION :: TEMPX, TEMPLEFT, TEMPRIGHT, TEMPPATCHPOS(NRBSITES,3)
  LOGICAL          :: GTEST
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
  DOUBLE PRECISION :: ALPHAMIN, ALPHA1T, ALPHA2T, VF1, DVF1, VF2, DVF2, DVEXP(1:3), DVF1DX(1:3)
  DOUBLE PRECISION :: S2A(0:3), S1A(0:3), AI(0:3), BI(0:3)
  

  !parameters for smoothing function VF2 (polynomial of grade 3) - has to fulfill
  !VF2(0.9*range) = VLJ(0.9*range)
  !VF2(range) = 0
  !VF2'(0.9*range) = VLJ'(0.9*range)
  !VF2'(range) = 0
  S2A(0:3) =  (/ 4.3196421660790450D1, -7.2976415249550730D1, 4.0919975890836945D1, -7.6195284520433670D0 /)

  ENERGY = 0.D0
  IF (GTEST) G(:) = 0.D0
  REALNATOMS = NATOMS/2
  OFFSET = 3*REALNATOMS
  
  DO J1 = 1, REALNATOMS
    
    J3 = 3*J1
    J5 = OFFSET + J3
    RI(J1,1:3) = X(J3-2:J3)
    P(J1,1:3)  = X(J5-2:J5)

    !calculate actual position of patches
    CALL RMDRVT(P(J1,:), RMI, DRMI1, DRMI2, DRMI3, GTEST)

    DO J2 = 1, NRBSITES
      PATCHPOS(J1,J2,:) = MATMUL(RMI(:,:),SITE(J2,:))
      IF (GTEST) THEN
        DPATCHPOS1(J1,J2,:) = MATMUL(DRMI1(:,:),SITE(J2,:))
        DPATCHPOS2(J1,J2,:) = MATMUL(DRMI2(:,:),SITE(J2,:))
        DPATCHPOS3(J1,J2,:) = MATMUL(DRMI3(:,:),SITE(J2,:))
      END IF
    END DO
       
  END DO


  DO J1 = 1, REALNATOMS

    DO J2 = J1+1, REALNATOMS

      RIJ(:) = RI(J1,:) - RI(J2,:)
      RSQ = DOT_PRODUCT(RIJ(:),RIJ(:))

      IF (RSQ.LE.0.81D0*SIGMASQ) THEN
      !core - LJ repulsion
        R2 = 1.D0/RSQ
        R6 = R2*R2*R2

        ENERGY = ENERGY + (R6 - 1.D0) * R6

        IF (GTEST) THEN
          !calculate derivatives
          DVLJ = -6.D0 * (2.D0 * R6 - 1.D0) * R6 * R2 !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVLJ * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVLJ * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.0.81D0*SIGMASQ).AND.(RSQ.LE.SIGMASQ)) THEN
      !near edge of core - smoothend repulsion
        ALPHAMIN = 2.D0*PI

        R = SQRT(RSQ)
        RCUB = R*RSQ
        R2 = 1.D0 / RSQ
        RHALF = R / 2.D0

        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
        END IF

        DO J3T = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3T,:),-RIJ(:)) / RHALF
          ALPHA1T = ACOS(ARG1)
          DO J4T = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4T,:),RIJ(:)) / RHALF
            ALPHA2T = ACOS(ARG2)
            IF ((ALPHA1T+ALPHA2T).LT.ALPHAMIN) THEN
              ALPHAMIN=ALPHA1T+ALPHA2T
              ALPHA1=ALPHA1T
              ALPHA2=ALPHA2T
              J3=J3T
              J4=J4T
            END IF
          END DO
        END DO

        VEXP1 = EXP(-ALPHA1**2/FACTOR)              
        VEXP2 = EXP(-ALPHA2**2/FACTOR)
        VEXP = VEXP1 * VEXP2
        IF (GTEST) THEN
          !derivatives of ARG1 wrt interparticle vector RIJ
          DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
          !derivatives of ARG1 wrt orientational vector P of particles J1
          DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
          DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
          DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
          DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)
          !derivatives of ARG2 wrt interparticle vector RIJ
          DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
          !derivatives of ARG2 wrt orientational vector P of particle J2
          DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
          DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
          DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
          DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)
        END IF

        !parameters for smoothing function VF1 (polynomial of grade 3) - has to fulfill
        !VF1(0.9*sigma) = VLJ(0.9*sigma)
        !VF1(sigma) = VLJ(sigma)*VANG(THETA) = 0
        !VF1'(0.9*sigma) = VLJ'(0.9*sigma)
        !VF1'(sigma) = VLJ'(sigma)*VANG(THETA)
        AI(0:3) = (/ 299.49098473874074D0,-747.4130927079459D0,596.3532311996697D0,-148.4311232304645D0 /)
        BI(0:3) = (/ 485.999999999995D0,-1565.9999999999836D0,1679.9999999999825D0,-599.9999999999939D0 /)
        S1A(0:3) = AI(0:3) + BI(0:3)*VEXP
        VF1 = S1A(0) + S1A(1) * R + S1A(2) * RSQ + S1A(3) * RCUB
        ENERGY = ENERGY + VF1


        IF (GTEST) THEN
          !calculate derivatives

          DVF1 = S1A(1)/R + 2.D0*S1A(2) + 3.D0*S1A(3)*R
          A0 = -((2.D0*VEXP)/FACTOR)*(BI(0)+BI(1)*R+BI(2)*RSQ+BI(3)*RCUB)
          DVEXP(1:3) = ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)
          DVF1DX(:) = DVF1 * RIJ(:) + A0 * DVEXP(:)
          !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVF1DX(1:3)
          !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVF1DX(1:3)


          A1 = A0*ALPHA1
          A2 = A0*ALPHA2

          !particle 1, derivatives wrt orientational coordinates
          G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
          !particle 2, derivatives wrt orientational coordinates
          G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)
          
        END IF

      ELSEIF ((RSQ.GT.SIGMASQ).AND.(RSQ.LE.0.81D0*RANGESQ)) THEN
      !patchy region - LJ attraction if patches are aligned
        R = SQRT(RSQ)
        R2 = 1.D0 / RSQ
        R6 = R2*R2*R2
        RHALF = R / 2.D0
        VLJ = (R6 - 1.D0) * R6


        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
          DVLJ = -6.D0 * (2.D0 * R6 - 1.D0) * R6 * R2 !factor 1/R from dR/dXi = Xi/R already included
        END IF

        DO J3 = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) / RHALF
          ALPHA1 = ACOS(ARG1)
          !some kind of angle-cutoff could be added..
          VEXP1 = EXP(-ALPHA1**2/FACTOR)

          IF (GTEST) THEN

            !derivatives of ARG1 wrt interparticle vector RIJ
            DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
            !derivatives of ARG1 wrt orientational vector P of particles J1
            DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
            DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)

          END IF 

          DO J4 = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) / RHALF
            ALPHA2 = ACOS(ARG2)
            !some kind of angle-cutoff could be added..
            VEXP2 = EXP(-ALPHA2**2/FACTOR)
            VEXP = VEXP1 * VEXP2
            V = VLJ * VEXP

            ENERGY = ENERGY + V

            IF (GTEST) THEN

              !derivatives of ARG2 wrt interparticle vector RIJ
              DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
              !derivatives of ARG2 wrt orientational vector P of particle J2
              DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
              DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
              DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
              DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)

              A0 = -2.D0 * VLJ / FACTOR
              DVDX(1:3) = VEXP * (DVLJ * RIJ(:) + A0 * (ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)))

              !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
              G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVDX(1:3)
              !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
              G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVDX(1:3)

              A1 = -2.D0 * V * ALPHA1 / FACTOR
              A2 = -2.D0 * V * ALPHA2 / FACTOR

              !particle 1, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
              !particle 2, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)

            END IF 

          END DO
        END DO

      ELSEIF ((RSQ.GT.0.81*RANGESQ).AND.(RSQ.LE.RANGESQ)) THEN
      !patchy region near cutoff - smoothend attraction if patches are aligned
        R = SQRT(RSQ)
        R2 = 1.D0/RSQ
        RHALF = R / 2.D0
        VF2 = S2A(0) + S2A(1) * R + S2A(2) * RSQ + S2A(3) * R * RSQ


        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
          DVF2 = S2A(1)/R  + 2.D0*S2A(2)  + 3.D0*S2A(3) * R !factor 1/R from dR/dXi = Xi/R already included
        END IF

        DO J3 = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) / RHALF
          ALPHA1 = ACOS(ARG1)
          !some kind of angle-cutoff could be added..
          VEXP1 = EXP(-ALPHA1**2/FACTOR)

          IF (GTEST) THEN

            !derivatives of ARG1 wrt interparticle vector RIJ
            DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
            !derivatives of ARG1 wrt orientational vector P of particles J1
            DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
            DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)

          END IF 

          DO J4 = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) / RHALF
            ALPHA2 = ACOS(ARG2)
            !some kind of angle-cutoff could be added..
            VEXP2 = EXP(-ALPHA2**2/FACTOR)
            VEXP = VEXP1 * VEXP2
            V = VF2 * VEXP

            ENERGY = ENERGY + V

            IF (GTEST) THEN

              !derivatives of ARG2 wrt interparticle vector RIJ
              DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
              !derivatives of ARG2 wrt orientational vector P of particle J2
              DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
              DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
              DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
              DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)

              A0 = -2.D0 * VF2 / FACTOR
              DVDX(1:3) = VEXP * (DVF2 * RIJ(:) + A0 * (ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)))

              !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
              G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVDX(1:3)
              !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
              G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVDX(1:3)

              A1 = -2.D0 * V * ALPHA1 / FACTOR
              A2 = -2.D0 * V * ALPHA2 / FACTOR

              !particle 1, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
              !particle 2, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)

            END IF 

          END DO
        END DO

      END IF

    END DO
  END DO

END SUBROUTINE PATCHYPOT


SUBROUTINE PATCHYPAIR (RJ1, RJ2, PATCHPOSJ1, PATCHPOSJ2, ENERGY)
!compute energy of a single pair interaction
!this can be used for calculating numerical derivatives - obsolete when smoothing functions are used

  USE COMMONS, ONLY: NRBSITES, SIGMASQ, RANGESQ, FACTOR

  IMPLICIT NONE
  
  INTEGER ::          J3, J4
  DOUBLE PRECISION :: ENERGY, RSQ, R2, R6, RHALF, ARG1, ARG2, ALPHA1, ALPHA2, VLJ, VEXP1, VEXP2
  DOUBLE PRECISION :: RJ1(3), RJ2(3), RIJ(3) 
  DOUBLE PRECISION :: PATCHPOSJ1(NRBSITES,3), PATCHPOSJ2(NRBSITES,3)
  

  ENERGY = 0.D0

  RIJ(:) = RJ1(:) - RJ2(:)
  RSQ = DOT_PRODUCT(RIJ(:),RIJ(:))

  IF (RSQ.LE.SIGMASQ) THEN
    R2 = 1.D0/RSQ
    R6 = R2*R2*R2
    ENERGY = (R6 - 1.D0) * R6

  ELSEIF ((RSQ.GT.SIGMASQ).AND.(RSQ.LT.RANGESQ)) THEN
    R2 = 1.D0 / RSQ
    R6 = R2*R2*R2
    RHALF = SQRT(RSQ) / 2.D0
    VLJ = (R6 - 1.D0) * R6

    DO J3 = 1, NRBSITES
      ARG1 = DOT_PRODUCT(PATCHPOSJ1(J3,:),-RIJ(:)) / RHALF
      ALPHA1 = ACOS(ARG1)
      !some kind of angle-cutoff could be added..
      VEXP1 = EXP(-ALPHA1**2/FACTOR)
      DO J4 = 1, NRBSITES
        ARG2 = DOT_PRODUCT(PATCHPOSJ2(J4,:),RIJ(:)) / RHALF
        ALPHA2 = ACOS(ARG2)
        !some kind of angle-cutoff could be added..
        VEXP2 = EXP(-ALPHA2**2/FACTOR)
        ENERGY = VLJ * VEXP1 * VEXP2
      END DO
    END DO
  END IF

END SUBROUTINE PATCHYPAIR



SUBROUTINE DEFINE_PATCHES(A)
!distribution of patches on particle surface - this is called from keyword.f
!A: geometry parameter

  USE COMMONS, ONLY: NRBSITES, SITE

  IMPLICIT NONE
  DOUBLE PRECISION :: A
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0



  IF (NRBSITES.EQ.3) THEN 

    SITE(1,1)=5.D-1*COS(0.D0)*SIN(A*PI/12.D0)
    SITE(1,2)=5.D-1*SIN(0.D0)*SIN(A*PI/12.D0)
    SITE(1,3)=5.D-1*COS(A*PI/12.D0)
    SITE(2,1)=5.D-1*COS(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(2,2)=5.D-1*SIN(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(2,3)=5.D-1*COS(A*PI/12.D0)
    SITE(3,1)=5.D-1*COS(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(3,2)=5.D-1*SIN(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(3,3)=5.D-1*COS(A*PI/12.D0)

  ELSEIF (NRBSITES.EQ.4) THEN

!!$    SITE(1,1)=5.D-1*COS(0.D0)*SIN(A*PI/12.D0)
!!$    SITE(1,2)=5.D-1*SIN(0.D0)*SIN(A*PI/12.D0)
!!$    SITE(1,3)=5.D-1*COS(A*PI/12.D0)
!!$    SITE(2,1)=5.D-1*COS(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
!!$    SITE(2,2)=5.D-1*SIN(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
!!$    SITE(2,3)=5.D-1*COS(A*PI/12.D0)
!!$    SITE(3,1)=5.D-1*COS(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
!!$    SITE(3,2)=5.D-1*SIN(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
!!$    SITE(3,3)=5.D-1*COS(A*PI/12.D0)
!!$    SITE(4,1)=5.D-1*COS(0.D0)*SIN(0.D0)
!!$    SITE(4,2)=5.D-1*SIN(0.D0)*SIN(0.D0)
!!$    SITE(4,3)=5.D-1*COS(0.D0)

!direct definition of tetrahedral distribution (equivalent: A~7.29808138D0)

     SITE(1,:)= (5.D-1/SQRT(3.D0))*(/  1.D0,  1.D0,  1.D0/)
     SITE(2,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0, -1.D0,  1.D0/)
     SITE(3,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0,  1.D0, -1.D0/)
     SITE(4,:)= (5.D-1/SQRT(3.D0))*(/  1.D0, -1.D0, -1.D0/)


  ELSEIF (NRBSITES.EQ.5) THEN

    SITE(1,1)=5.D-1*COS(0.D0)*SIN(A*PI/12.D0)
    SITE(1,2)=5.D-1*SIN(0.D0)*SIN(A*PI/12.D0)
    SITE(1,3)=5.D-1*COS(A*PI/12.D0)
    SITE(2,1)=5.D-1*COS(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(2,2)=5.D-1*SIN(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(2,3)=5.D-1*COS(A*PI/12.D0)
    SITE(3,1)=5.D-1*COS(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(3,2)=5.D-1*SIN(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
    SITE(3,3)=5.D-1*COS(A*PI/12.D0)
    SITE(4,1)=5.D-1*COS(0.D0)*SIN(0.D0)
    SITE(4,2)=5.D-1*SIN(0.D0)*SIN(0.D0)
    SITE(4,3)=5.D-1*COS(0.D0)
    SITE(5,1)=5.D-1*COS(0.D0)*SIN(PI)
    SITE(5,2)=5.D-1*SIN(0.D0)*SIN(PI)
    SITE(5,3)=5.D-1*COS(PI)

  ELSEIF (NRBSITES.EQ.6) THEN

    SITE(1,1)=0.D0
    SITE(1,2)=0.D0
    SITE(1,3)=5.D-1
    SITE(2,1)=-5.D-1
    SITE(2,2)=0.D0
    SITE(2,3)=0.D0
    SITE(3,1)=5.D-1
    SITE(3,2)=0.D0
    SITE(3,3)=0.D0
    SITE(4,1)=0.D0
    SITE(4,2)=0.D0
    SITE(4,3)=-5.D-1
    SITE(5,1)=0.D0
    SITE(5,2)=5.D-1
    SITE(5,3)=0.D0
    SITE(6,1)=0.D0
    SITE(6,2)=-5.D-1
    SITE(6,3)=0.D0

  ELSE

    PRINT*,'Patchnumber ',NRBSITES,' not implemented.'
    STOP

  END IF


END SUBROUTINE DEFINE_PATCHES



SUBROUTINE VIEWPATCHY()
!generate output files for visualization (coordinates of 'best' configuration)
!  USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
  USE COMMONS, ONLY: NATOMS, NRBSITES, SITE
  USE QMODULE
  IMPLICIT NONE

  INTEGER          :: I, J1, J2, J3, J5, J7
  DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
  LOGICAL          :: GTEST

  OPEN(UNIT=26, FILE='particles.xyz', STATUS='UNKNOWN')
  OPEN(UNIT=27, FILE='patches.xyz', STATUS='UNKNOWN')

  GTEST = .FALSE. 

  DO J1 = 1, 1!NSAVE

!    WRITE(26,'(I6)') (NATOMS/2)*3
!    WRITE(26,10) J1, QMIN(J1), FF(J1)
!10  FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

     WRITE(26,'(I4)') (NATOMS/2)
     WRITE(26,'(A)') ' '
     WRITE(27,'(I4)') (NATOMS/2)*NRBSITES
     WRITE(27,'(A)') ' '

    DO J3 = 1, NATOMS/2

      J5   = 3*J3
      J7   = 3*NATOMS/2 + J5
      P(:) = QMINP(J1,J7-2:J7)

      CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

      DO J2 = 1, NRBSITES+1
        IF (J2 == 1) THEN
          RBCOORDS(1:3) = QMINP(J1,J5-2:J5)
          WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
        ELSE
          RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2-1,:))
          WRITE(27,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
        ENDIF
      ENDDO

    ENDDO
  ENDDO
  CLOSE (UNIT=26)
  CLOSE (UNIT=27)

END SUBROUTINE VIEWPATCHY



SUBROUTINE VIEWPATCHYCURR(X)
!generate output file for visualization (coordinates explicitly passed to subroutine)
  USE COMMONS, ONLY: NATOMS, NRBSITES, SITE
  USE QMODULE
  IMPLICIT NONE

  INTEGER          :: I, I1, I2, I3, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), RBCOORDS(3)
  DOUBLE PRECISION :: ENERGY, R, RSQ, R2, R6, RHALF, ARG1, ARG2, ALPHA1, ALPHA2, VLJ, VEXP1, VEXP2, VEXP, V, DVLJ, A0, A1, A2
  DOUBLE PRECISION :: DARG1(6), DALPHA1(6), DARG2(6), DALPHA2(6), DVDX(3)
  DOUBLE PRECISION :: RI(NATOMS/2,3), RIJ(3), P(NATOMS/2,3), A(3) 
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: PATCHPOS(NATOMS/2,NRBSITES,3), DPATCHPOS1(NATOMS/2,NRBSITES,3)
  DOUBLE PRECISION :: DPATCHPOS2(NATOMS/2,NRBSITES,3), DPATCHPOS3(NATOMS/2,NRBSITES,3)
  DOUBLE PRECISION :: TEMPX, TEMPLEFT, TEMPRIGHT, TEMPPATCHPOS(NRBSITES,3)
  LOGICAL          :: GTEST
  

  REALNATOMS = NATOMS/2
  OFFSET = 3*REALNATOMS
  
  DO J1 = 1, REALNATOMS
    
    J3 = 3*J1
    J5 = OFFSET + J3
    RI(J1,1:3) = X(J3-2:J3)
    P(J1,1:3)  = X(J5-2:J5)

    !calculate actual position of patches
    CALL RMDRVT(P(J1,:), RMI, DRMI1, DRMI2, DRMI3, .FALSE.)

    DO J2 = 1, NRBSITES
      PATCHPOS(J1,J2,:) = MATMUL(RMI(:,:),SITE(J2,:))
    END DO
       
  END DO


  OPEN(UNIT=26, FILE='particles.xyz', STATUS='UNKNOWN')
  OPEN(UNIT=27, FILE='patches.xyz', STATUS='UNKNOWN')

  GTEST = .FALSE. 

  DO J1 = 1, 1

!    WRITE(26,'(I6)') (NATOMS/2)*3
!    WRITE(26,10) J1, QMIN(J1), FF(J1)
!10  FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

     WRITE(26,'(I4)') (NATOMS/2)
     WRITE(26,'(A)') ' '
     WRITE(27,'(I4)') (NATOMS/2)*NRBSITES
     WRITE(27,'(A)') ' '

    DO J3 = 1, NATOMS/2

      DO J2 = 1, NRBSITES+1
        IF (J2 == 1) THEN
          RBCOORDS(1:3) = RI(J3,:)
          WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
        ELSE
          RBCOORDS(1:3) = RI(J3,:) + PATCHPOS(J3,J2-1,:)
          WRITE(27,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
        ENDIF
      ENDDO

    ENDDO
  ENDDO
  CLOSE (UNIT=26)
  CLOSE (UNIT=27)

END SUBROUTINE VIEWPATCHYCURR



SUBROUTINE DISPLAYGRADIENT(X)
!for bughunting
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: I, J1, J2, J3
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), RI(NATOMS/2,3), RIJ(3), R
  DOUBLE PRECISION :: ENERGY

  CALL PATCHYPOT (X, G, ENERGY, .TRUE.)

  OPEN(9,FILE='gradtest.dat')

  DO J1 = 1,NATOMS/2
    
    J3 = 3*J1
    RI(J1,1:3) = X(J3-2:J3)
       
  END DO

  WRITE(9,*) SQRT(DOT_PRODUCT(G(:),G(:)))

  DO I=1,3*NATOMS
    IF (ABS(G(I)).GT.1.D-2) THEN
      IF (I.LE.3*(NATOMS/2)) THEN
        WRITE(9,*) (I+2)/3,MOD(I,3),G(I)
      ELSE
        WRITE(9,*) (I-3*(NATOMS/2)+2)/3,MOD(I,3)+3,G(I)
      END IF
    END IF
  END DO

  WRITE(9,*) '___'
  
  DO J1=1,NATOMS/2
    DO J2=1,NATOMS/2
      RIJ(:) = RI(J1,:) - RI(J2,:)
      R = SQRT(DOT_PRODUCT(RIJ(:),RIJ(:)))
      IF (R.LE.2.5D0) THEN
        WRITE(9,*) J1,J2,R
      END IF
    END DO
  END DO

  CLOSE(9)

END SUBROUTINE

!<gd351|
