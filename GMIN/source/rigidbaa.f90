! Finds the rotation matrix and its derivatives with respect to 
! components of the rotation vector.
!
! See 'Simulations of rigid bodies in an angle axis framework', Dwaipayan
! Chakrabarti and David Wales, Phys. Chem. Chem. Phys., 2009, 11, 
! 1970-1976
!
!     ----------------------------------------------------------------------------------------------
! jdf43>        Modified for general angle-axis 30/01/12
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTMAT(P, RM)

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, CT, ST, I3(3,3), E(3,3), RM(3,3)

      I3(:,:) = 0.D0
      DO K = 1, 3
         I3(K,K) = 1.D0
      ENDDO

      THETA2 = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN
         RM(:,:) = I3(:,:)
      ELSE
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         RM      = I3(:,:) + (1.D0-CT)*MATMUL(E(:,:),E(:,:)) + ST*E(:,:)
      ENDIF

      END SUBROUTINE ROTMAT

!     --------------------------------------------------------------------------

!     RMDVDT = rotation matrix derivative
!     P(3) = rotation vector 
!     RM(3,3) = rotation matrix
!     DRMk(3,3) = derivative of the rotation matrix with respect to the 
!                 kth component of the rotation vector
!     GTEST = true if derivatives are to be found  

      SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)

      IMPLICIT NONE

!     PN(3) = the unit vector parallel to P
!     THETA = the modulus of the rotation vector P, equivalent to the 
!             angle of rotation
!     THETA2 = THETA squared
!     THETA3 = THETA**-3
!     CT = cos(THETA)
!     ST = sin(THETA)
!     I3(3,3) = 3x3 identity matrix
!     E(3,3) = the skew-symmetric matrix obtained from a unit vector 
!              parallel to P (equation (2) in the paper)
!     ESQ(3,3) = the square of E
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)

!     DEk = derivate of E with respect to the kth component of P
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      LOGICAL          :: GTEST

!     Set the values of the idenity matrix I3
      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!     Calculate the value of THETA2 as the square modulus of P
      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 < 1.0D-12) THEN
!        Execute if the angle of rotation is zero
!        In this case the rotation matrix is the identity matrix
         RM(:,:) = I3(:,:)

         ! vr274> first order corrections to rotation matrix
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

!        If derivatives do not need to found, we're finished
         IF (.NOT. GTEST) RETURN

!        This is the special case described in the paper, where DRMk is
!        equal to E(k), which is the skew-symmetric matrix obtained from
!        P with P(k) equal to 1 and other components equal to zero
!         PN        = (/1.D0, 0.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(2,3)    = -PN(1)
!         E(3,2)    = -E(2,3)
!         DRM1(:,:) = E(:,:)

!         PN        = (/0.D0, 1.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(1,3)    =  PN(2)
!         E(3,1)    = -E(1,3)
!         DRM2(:,:) = E(:,:)

!         PN        = (/0.D0, 0.D0, 1.D0/)
!         E(:,:)    = 0.D0
!         E(1,2)    = -PN(3)
!         E(2,1)    = -E(1,2)
!         DRM3(:,:) = E(:,:)

! hk286 - now up to the linear order in theta
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    = -2.0D0*P(1)
         E(2,3)    = -2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    = -2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    = -2.0D0
         E(3,2)    = P(3)
         E(3,3)    = -2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(3)
         E(1,2)    = -2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    = -2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

      ELSE
!       Execute for the general case, where THETA dos not equal zero
!       Find values of THETA, CT, ST and THETA3
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)

!        Set THETA to 1/THETA purely for convenience
         THETA   = 1.D0/THETA

!        Normalise P and construct the skew-symmetric matrix E
!        ESQ is calculated as the square of E
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         ESQ     = MATMUL(E,E)

!        RM is calculated from Rodrigues' rotation formula (equation (1)
!        in the paper)
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

!        If derivatives do not need to found, we're finished
         IF (.NOT. GTEST) RETURN

!        Set up DEk using the form given in equation (4) in the paper
         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

!        Use equation (3) in the paper to find DRMk
         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

      ENDIF
      END SUBROUTINE RMDRVT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RMDFAS(P, RM, DRM1, DRM2, DRM3, D2RM1, D2RM2, D2RM3, D2RI12, D2RI23, D2RI31, GTEST, STEST)

      IMPLICIT NONE

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, E(3,3), ESQ(3,3), I3(3,3) 
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      DOUBLE PRECISION :: D2E1(3,3), D2E2(3,3), D2E3(3,3), D2E12(3,3), D2E23(3,3), D2E31(3,3)
      DOUBLE PRECISION :: D2RM1(3,3), D2RM2(3,3), D2RM3(3,3)
      DOUBLE PRECISION :: D2RI12(3,3), D2RI23(3,3), D2RI31(3,3)
      DOUBLE PRECISION :: FCTR, FCTRSQ1, FCTRSQ2, FCTRSQ3
      LOGICAL          :: GTEST, STEST

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!      RM(3,3) = 0.D0; DRM1(3,3) = 0.D0; DRM2(3,3) = 0.D0; DRM3(3,3) = 0.D0 
!      D2RM1(3,3) = 0.D0; D2RM2(3,3) = 0.D0; D2RM3(3,3) = 0.D0
!      D2RI12(3,3) = 0.D0; D2RI23(3,3) = 0.D0; D2RI31(3,3) = 0.D0

      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN

         RM(:,:)     = I3(:,:)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

         PN          = (/1.D0, 0.D0, 0.D0/)
         E(:,:)      = 0.D0
         E(2,3)      = -PN(1)
         E(3,2)      = -E(2,3)
         DRM1(:,:)   = E(:,:)

         PN          = (/0.D0, 1.D0, 0.D0/)
         E(:,:)      = 0.D0
         E(1,3)      = PN(2)
         E(3,1)      = -E(1,3)
         DRM2(:,:)   = E(:,:)

         PN          = (/0.D0, 0.D0, 1.D0/)
         E(:,:)      = 0.D0
         E(1,2)      = -PN(3)
         E(2,1)      = -E(1,2)
         DRM3(:,:)   = E(:,:)

         IF (.NOT. STEST) RETURN

         D2RM1(:,:)  = MATMUL(DRM1(:,:),DRM1(:,:))
         D2RM2(:,:)  = MATMUL(DRM2(:,:),DRM2(:,:))
         D2RM3(:,:)  = MATMUL(DRM3(:,:),DRM3(:,:))

         D2RI12(:,:) = MATMUL(DRM1(:,:),DRM2(:,:))
         D2RI23(:,:) = MATMUL(DRM2(:,:),DRM3(:,:))
         D2RI31(:,:) = MATMUL(DRM3(:,:),DRM1(:,:))

      ELSE

         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)

         ESQ     = MATMUL(E(:,:),E(:,:))
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

         IF (.NOT. STEST) RETURN

         FCTRSQ1 = PN(1)*PN(1) 
         FCTRSQ2 = PN(2)*PN(2)
         FCTRSQ3 = PN(3)*PN(3)

         D2E1(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(1)*PN(1))*THETA3
         D2E1(1,2) =  P(3)*FCTR
         D2E1(1,3) = -P(2)*FCTR
         D2E1(2,3) =  P(1)*FCTR + 2.D0*P(1)*THETA3
         D2E1(2,1) = -D2E1(1,2)
         D2E1(3,1) = -D2E1(1,3)
         D2E1(3,2) = -D2E1(2,3)

         D2E2(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(2)*PN(2))*THETA3
         D2E2(1,2) =  P(3)*FCTR
         D2E2(1,3) = -P(2)*FCTR - 2.D0*P(2)*THETA3
         D2E2(2,3) =  P(1)*FCTR 
         D2E2(2,1) = -D2E2(1,2)
         D2E2(3,1) = -D2E2(1,3)
         D2E2(3,2) = -D2E2(2,3)

         D2E3(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(3)*PN(3))*THETA3
         D2E3(1,2) =  P(3)*FCTR + 2.D0*P(3)*THETA3
         D2E3(1,3) = -P(2)*FCTR
         D2E3(2,3) =  P(1)*FCTR 
         D2E3(2,1) = -D2E3(1,2)
         D2E3(3,1) = -D2E3(1,3)
         D2E3(3,2) = -D2E3(2,3)

         D2E12(:,:) =  0.D0
         D2E12(1,2) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2 
         D2E12(1,3) = -PN(1)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E12(2,3) =  PN(2)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E12(2,1) = -D2E12(1,2)
         D2E12(3,1) = -D2E12(1,3)
         D2E12(3,2) = -D2E12(2,3)

         D2E23(:,:) =  0.D0
         D2E23(1,2) =  PN(2)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E23(1,3) = -PN(3)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E23(2,3) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E23(2,1) = -D2E23(1,2)
         D2E23(3,1) = -D2E23(1,3)
         D2E23(3,2) = -D2E23(2,3)

         D2E31(:,:) =  0.D0
         D2E31(1,2) =  PN(1)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E31(1,3) =  3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E31(2,3) =  PN(3)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E31(2,1) = -D2E31(1,2)
         D2E31(3,1) = -D2E31(1,3)
         D2E31(3,2) = -D2E31(2,3)

         D2RM1(:,:)  = ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (CT*FCTRSQ1 - ST*FCTRSQ1*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE1(:,:),DE1(:,:)) + MATMUL(D2E1(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E1(:,:))) + (- ST*FCTRSQ1 - CT*FCTRSQ1*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(1)*DE1(:,:) + ST*D2E1(:,:)
         
         D2RM2(:,:)  = ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (CT*FCTRSQ2 - ST*FCTRSQ2*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE2(:,:),DE2(:,:)) + MATMUL(D2E2(:,:),E(:,:)) & 
                     + MATMUL(E(:,:),D2E2(:,:))) + (- ST*FCTRSQ2 - CT*FCTRSQ2*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(2)*DE2(:,:) + ST*D2E2(:,:)

         D2RM3(:,:)  = ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (CT*FCTRSQ3 - ST*FCTRSQ3*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (1.D0-CT)*(2.D0*MATMUL(DE3(:,:),DE3(:,:)) + MATMUL(D2E3(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E3(:,:))) + (- ST*FCTRSQ3 - CT*FCTRSQ3*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(3)*DE3(:,:) + ST*D2E3(:,:)

         D2RI12(:,:) = ST*PN(1)*(MATMUL(E(:,:),DE2(:,:)) + MATMUL(DE2(:,:),E(:,:))) &
                     + (CT*PN(1)*PN(2) - ST*PN(1)*PN(2)*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E12(:,:),E(:,:)) + MATMUL(DE1(:,:),DE2(:,:)) & 
                     + MATMUL(DE2(:,:),DE1(:,:)) + MATMUL(E(:,:),D2E12(:,:))) & 
                     - (ST*PN(1)*PN(2) + CT*PN(1)*PN(2)*THETA)*E(:,:) + CT*PN(1)*DE2(:,:) + CT*PN(2)*DE1(:,:) &
                     + ST*D2E12(:,:)

         D2RI23(:,:) = ST*PN(2)*(MATMUL(E(:,:),DE3(:,:)) + MATMUL(DE3(:,:),E(:,:))) &
                     + (CT*PN(2)*PN(3) - ST*PN(2)*PN(3)*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E23(:,:),E(:,:)) + MATMUL(DE2(:,:),DE3(:,:)) &
                     + MATMUL(DE3(:,:),DE2(:,:)) + MATMUL(E(:,:),D2E23(:,:))) &
                     - (ST*PN(2)*PN(3) + CT*PN(2)*PN(3)*THETA)*E(:,:) + CT*PN(2)*DE3(:,:) + CT*PN(3)*DE2(:,:) &
                     + ST*D2E23(:,:)

         D2RI31(:,:) = ST*PN(3)*(MATMUL(E(:,:),DE1(:,:)) + MATMUL(DE1(:,:),E(:,:))) &
                     + (CT*PN(3)*PN(1) - ST*PN(3)*PN(1)*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E31(:,:),E(:,:)) + MATMUL(DE3(:,:),DE1(:,:)) &
                     + MATMUL(DE1(:,:),DE3(:,:)) + MATMUL(E(:,:),D2E31(:,:))) &
                     - (ST*PN(3)*PN(1) + CT*PN(3)*PN(1)*THETA)*E(:,:) + CT*PN(3)*DE1(:,:) + CT*PN(1)*DE3(:,:) &
                     + ST*D2E31(:,:)

      ENDIF
 
      END SUBROUTINE RMDFAS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SHIFTRIGID (Q, NATOMS)

!     THIS SUBROUTINE SHIFTS THE 'ZERO' EIGENVALUES CORRESPONDING TO OVERALL TRANSLATION AND
!     ROTATION OF A SYSTEM OF (IDENTICAL) RIGID BODIES WITH C-O-M & ANGLE-AXIS COORDINATES.

      USE COMMONS, ONLY: EFIELDT, PYT, GBT, GBDT, STOCKAAT, PYGPERIODICT, PYBINARYT, SHIFTED, SHIFTV, NZERO
      USE MODHESS

      IMPLICIT NONE

      INTEGER            :: NATOMS, I, J, J1, J2
      DOUBLE PRECISION   :: Q(3*NATOMS), EV(3*NATOMS,6+NATOMS/2), NRMFCT(6+NATOMS/2)
      DOUBLE PRECISION   :: CMX, CMY, CMZ, THETA, THETA2, THETAH, DUMMY
      LOGICAL            :: UNIAXT            

!     INITIALIZE

      EV(:,:)   = 0.D0
      NRMFCT(:) = 0.D0
      CMX       = 0.D0
      CMY       = 0.D0
      CMZ       = 0.D0
      UNIAXT    = .FALSE. !jdf43>

      SHIFTED = .TRUE.
      IF (EFIELDT) THEN
          NZERO = 4
      ELSE
          NZERO = 6
      ENDIF

      DO I = 1, NATOMS/2

         J = 3*I

         CMX = CMX + Q(J-2)
         CMY = CMY + Q(J-1)
         CMZ = CMZ + Q(J)

      ENDDO

      CMX = CMX / FLOAT(NATOMS/2)
      CMY = CMY / FLOAT(NATOMS/2)
      CMZ = CMZ / FLOAT(NATOMS/2)

      DO I = 1, NATOMS/2

         J  = 3*I
         J1 = 3*NATOMS/2 + J

         THETA2 = DOT_PRODUCT(Q(J1-2:J1), Q(J1-2:J1))
         THETA  = DSQRT(THETA2)
         THETAH = 0.5D0*THETA

!     TRANSLATION ALONG X
         EV(J-2,1) = 1.D0
         NRMFCT(1) = NRMFCT(1) + 1.D0

!     TRANSLATION ALONG Y
         EV(J-1,2) = 1.D0
         NRMFCT(2) = NRMFCT(2) + 1.D0

!     TRANSLATION ALONG Z
         EV(J,3) = 1.D0
         NRMFCT(3) = NRMFCT(3) + 1.D0

!     IF EFIELD IS PRESENT, IT HAS TO BE ALONG THE Z-AXIS IN THE LAB FRAME BY CONVENTION.

         IF (THETA == 0.D0) THEN

!     ROTATION ABOUT Z
            EV(J-2,4)  = - Q(J-1) + CMY
            EV(J-1,4)  = Q(J-2) - CMX
            EV(J1,4)   = 1.D0
            NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2

            IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
               EV(J-1,5)  = - Q(J) + CMZ
               EV(J,5)    = Q(J-1) - CMY
               EV(J1-2,5) = 1.D0
               NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
               EV(J-2,6)  = Q(J) - CMZ
               EV(J,6)    = - Q(J-2) + CMX
               EV(J1-1,6) = 1.D0
               NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

            ENDIF

         ELSE

!     ROTATION ABOUT Z
           EV(J-2,4)  = - Q(J-1) + CMY
           EV(J-1,4)  = Q(J-2) - CMX
           EV(J1-2,4) = - 0.5D0*Q(J1-1) + Q(J1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1)*Q(J1-2)/(THETA*TAN(THETAH))
           EV(J1-1,4) = 0.5D0*Q(J1-2) + Q(J1)*Q(J1-1)/THETA2 - 0.5D0*Q(J1)*Q(J1-1)/(THETA*TAN(THETAH))
           EV(J1,4)   = THETAH/TAN(THETAH) + Q(J1)**2/THETA2 - 0.5D0*Q(J1)**2/(THETA*TAN(THETAH))
           NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
             
           IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
              EV(J-1,5)  = - Q(J) + CMZ
              EV(J,5)    = Q(J-1) - CMY
              EV(J1-2,5) = THETAH/TAN(THETAH) + Q(J1-2)**2/THETA2 - 0.5D0*Q(J1-2)**2/(THETA*TAN(THETAH))
              EV(J1-1,5) = - 0.5D0*Q(J1) + Q(J1-2)*Q(J1-1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1-1)/(THETA*TAN(THETAH))
              EV(J1,5)   = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
              NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
              EV(J-2,6)  = Q(J) - CMZ
              EV(J,6)    = - Q(J-2) + CMX
              EV(J1-2,6) = 0.5D0*Q(J1) + Q(J1-1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1-1)*Q(J1-2)/(THETA*TAN(THETAH))
              EV(J1-1,6) = THETAH/TAN(THETAH) + Q(J1-1)**2/THETA2 - 0.5D0*Q(J1-1)**2/(THETA*TAN(THETAH))
              EV(J1,6)   = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
              NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

            ENDIF

         ENDIF

         IF (GBT.OR.GBDT.OR.STOCKAAT.OR.((PYGPERIODICT.OR.PYBINARYT.OR.PYT).AND.UNIAXT)) THEN

!     ROTATION ABOUT THE SYMMETRY AXIS
            IF (THETA == 0.D0) THEN

               EV(J1,NZERO+I)  = 1.D0 
               NRMFCT(NZERO+I) = NRMFCT(NZERO+I) + EV(J1,NZERO+I)**2
            
            ELSE 
               EV(J1-2,NZERO+I) = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1-1,NZERO+I) = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1,NZERO+I)   = THETAH*SIN(THETA) + Q(J1)*Q(J1)/THETA2 &
                                + 0.5D0*(THETA*COS(THETA)-Q(J1)*Q(J1)/THETA)/TAN(THETAH)
               NRMFCT(NZERO+I)  = NRMFCT(NZERO+I) + EV(J1-2,NZERO+I)**2 + EV(J1-1,NZERO+I)**2 + EV(J1,NZERO+I)**2  
            ENDIF
       
         ENDIF

      ENDDO

      IF (GBT.OR.GBDT.OR.STOCKAAT.OR.((PYGPERIODICT.OR.PYBINARYT.OR.PYT).AND.UNIAXT)) THEN
         NZERO = NZERO + NATOMS/2
      ELSE
         NZERO = NZERO
      ENDIF

      DO J = 1, NZERO

         NRMFCT(J) = DSQRT(NRMFCT(J))
         EV(:,J)   = EV(:,J)/NRMFCT(J)

      ENDDO

!     GRAM-SCHMIDT ORTHOGONALISATION TO OBTAIN ORTHONORMAL ROTATIONAL EIGENVECTORS

      DO J = 4, NZERO

         DO J1 = 4, J-1

            EV(:,J) = EV(:,J) - DOT_PRODUCT(EV(:,J),EV(:,J1))*EV(:,J1)

         ENDDO

         EV(:,J) = EV(:,J) / DSQRT(DOT_PRODUCT(EV(:,J),EV(:,J)))

      ENDDO

      DO J1 = 1, 3*NATOMS

         DO J2 = 1, 3*NATOMS

            DO J = 1, NZERO 

               HESS(J2,J1) = HESS(J2,J1) + SHIFTV*EV(J2,J)*EV(J1,J)

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE SHIFTRIGID
 
      SUBROUTINE RBMINDIST(RA,RB,NATOMS,DIST,Q2,DEBUG)

!     Follows the prescription of Kearsley, Acta Cryst. A, 45, 208-210, 1989, making necessary changes 
!     to conform to right-handed rotation in the right-handed coordinate system.
!     Brings RB to the best alignment with RA at the centre of mass of RA
!     Returns DIST as the actual distance, rather than the squared distance 

      USE COMMONS, ONLY: NTSITES, DBPT, DBPTDT, MSSTOCKT, STOCKAAT, EFIELDT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3)
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG

      IF ((DBPT .AND. EFIELDT) .OR. (DBPTDT .AND. EFIELDT) .OR. (MSSTOCKT .AND. EFIELDT) .OR. (STOCKAAT .AND. EFIELDT)) THEN

         CALL FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)
         RETURN

      ENDIF 

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      CALL SITEPOS(RA,XA)

      CALL SITEPOS(RB,XB)

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
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4) 
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)

      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','rbmindist> WARNING - INFO=',INFO,' in DSYEV'

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
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)

      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE RBMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBNEWROTGEOM(NATOMS,COORDS,Q2,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS/2
  
         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ
      
!     CONVERT THE ANGLE-AXIS COORDINATES

         J      = 3*NATOMS/2 + J
         P(:)   = COORDS(J+1:J+3)

         CALL QROTAA(Q2,P)

         COORDS(J+1:J+3) = P(1:3)

      ENDDO

      END SUBROUTINE RBNEWROTGEOM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)

!     returns DIST as the actual distance, rather than the squared distance

      USE COMMONS, ONLY: NTSITES, STOCKAAT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(2,2), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: DEBUG

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      IF (STOCKAAT) THEN

         XA(:) = RA(:)
         XB(:) = RB(:)
   
      ELSE

         CALL SITEPOS(RA,XA)

         CALL SITEPOS(RB,XB)

      ENDIF

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

      QMAT(1:2,1:2) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1) 
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + XP**2 + YP**2 + ZM**2
      ENDDO

!     QMAT IS SYMMETRIC; QMAT(2,1) = QMAT(1,2)

      MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))

!      Q2(1) = SQRT(QMAT(1,2)*QMAT(1,2)/((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)))
      Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
      Q2(2) = 0.D0
      Q2(3) = 0.D0
!      Q2(4) = (MINV - QMAT(1,1))*Q2(1)/QMAT(1,2)
      Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))

      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)
  
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE FLDMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBNEWROTGEOMMYORIENT(NATOMS,COORDS,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

      DOUBLE PRECISION :: r1, diag1, diag2, diag3
      INTEGER          :: u,v,w
!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA
!     sf344> when called from minpermdist, RM is already given from putting the centres of rigid bodies 
!            into standard orientation, so this subroutine simply rotates all coordinates with that rotation
!            matrix

!      RM(1,1) = Q2(1)**2 + Q2(2)**2 - Q2(3)**2 - Q2(4)**2
!      RM(1,2) = 2.D0*(Q2(2)*Q2(3) - Q2(1)*Q2(4))
!      RM(1,3) = 2.D0*(Q2(2)*Q2(4) + Q2(1)*Q2(3))
!      RM(2,1) = 2.D0*(Q2(2)*Q2(3) + Q2(1)*Q2(4))
!      RM(2,2) = Q2(1)**2 + Q2(3)**2 -Q2(2)**2-Q2(4)**2
!      RM(2,3) = 2.D0*(Q2(3)*Q2(4) - Q2(1)*Q2(2))
!      RM(3,1) = 2.D0*(Q2(2)*Q2(4) - Q2(1)*Q2(3))
!      RM(3,2) = 2.D0*(Q2(3)*Q2(4) + Q2(1)*Q2(2))
!      RM(3,3) = Q2(1)**2 +Q2(4)**2 -Q2(2)**2 - Q2(3)**2

      diag1=RM(1,1)
      diag2=RM(2,2)
      diag3=RM(3,3)

       
!     if the rotation matrix is the identity matrix, then return
      IF(abs(diag1-1.0D0)<1.0D-8.AND.abs(diag2-1.0D0)<1.0D-8.AND.abs(diag3-1.0D0)<1.0D-8) THEN
        RETURN
      END IF
!     otherwise figure out the quaternion from the rotation matrix
!      WRITE(*,'(A)') 'coords before rotation'
!      WRITE(*,'(3F13.8)') COORDS(:)

      IF(ABS(RM(1,1))>=ABS(RM(2,2)).AND.ABS(RM(1,1))>=ABS(RM(3,3))) THEN
          diag1=RM(1,1)
          u=1
          v=2
          w=3
      ELSE IF (ABS(RM(2,2))>=ABS(RM(1,1)).AND.ABS(RM(2,2))>=ABS(RM(3,3))) THEN
          diag1=RM(2,2)
          diag2=RM(3,3)
          diag3=RM(1,1)
          u=2
          v=3
          w=1
      ELSE IF (ABS(RM(3,3))>=ABS(RM(1,1)).AND.ABS(RM(3,3))>=ABS(RM(2,2))) THEN
          diag1=RM(3,3)
          diag2=RM(1,1)
          diag3=RM(2,2)
          u=3
          v=1
          w=2
      END IF

      r1=SQRT(1+RM(u,u)-RM(v,v)-RM(w,w))
      
      Q2(1)=(RM(w,v)-RM(v,w))/(2.0D0*r1)
      Q2(u+1)=r1/2.0D0
      Q2(v+1)=(RM(u,v)+RM(v,u))/(2.0D0*r1)
      Q2(w+1)=(RM(w,u)+RM(u,w))/(2.0D0*r1)

      DO I = 1, NATOMS/2
  
         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ
      
!     CONVERT THE ANGLE-AXIS COORDINATES

         J      = 3*NATOMS/2 + J
         P(:)   = COORDS(J+1:J+3)
         THETA  = DSQRT(DOT_PRODUCT(P,P))
         THETAH = 0.5D0*THETA
         ST     = SIN(THETAH)
!        WRITE(*,*) 'st=', THETAH, ST, P(:)
         Q1(1)  = COS(THETAH)
         Q1(2)  = P(1)*ST/THETA
         Q1(3)  = P(2)*ST/THETA
         Q1(4)  = P(3)*ST/THETA

         Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
         Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
         Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
         Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

         THETA  = 2.D0*ACOS(Q(1))

         IF (THETA == 0.D0) THEN
            COORDS (J+1:J+3) = 0.D0
         ELSE
            FCT    = DSQRT(DOT_PRODUCT(Q(2:4),Q(2:4)))
            COORDS(J+1) = THETA*Q(2)/FCT 
            COORDS(J+2) = THETA*Q(3)/FCT
            COORDS(J+3) = THETA*Q(4)/FCT
         ENDIF

      ENDDO
!      WRITE(*,'(A)') 'coords after rotation'
!      WRITE(*,'(3F13.8)') COORDS(:)

      END SUBROUTINE RBNEWROTGEOMMYORIENT

!     ----------------------------------------------------------------------------------------------


SUBROUTINE UNIAXGETPATHLENGTH(RA,RB,TEMP)

USE COMMONS, ONLY : NATOMS, PYA1
implicit none

integer          :: NSIZE, J1, J2, J3, J4, NRBSITES
double precision :: XA(9*NATOMS/2), XB(9*NATOMS/2), RA(3*NATOMS), RB(3*NATOMS), P(3), RM(3,3), RBSITE(3,3), TEMP
double precision :: R(3)

      NRBSITES = 3
      NSIZE = NRBSITES*NATOMS/2

      RBSITE(1,1) = 0.0D0
      RBSITE(1,2) = 0.0D0
      RBSITE(1,3) = 0.0D0
      RBSITE(2,1) = 0.0D0
      RBSITE(2,2) = 0.0D0
      RBSITE(2,3) = PYA1(3)
      RBSITE(3,1) = 0.0D0
      RBSITE(3,2) = 0.0D0
      RBSITE(3,3) =-PYA1(3)

!      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!      write(*,*) 'ra'
!      write(*,*) ra(:)
!     ----------------------------------------------------------------------------------------------
!      CALL LWOTPGH(RB,VNEW,ENERGY,.TRUE.,.FALSE.)
!      WRITE(*,*) ENERGY, SQRT(DOT_PRODUCT(VNEW(:),VNEW(:)))
!      CALL DUMBBELLP(RB,VNEW,ENERGY,.TRUE.,.FALSE.)
!      WRITE(*,*) ENERGY, SQRT(DOT_PRODUCT(VNEW(:),VNEW(:)))
!     ----------------------------------------------------------------------------------------------
!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      DO J1 = 1, NATOMS/2

         J2   = 3*J1
         R(:) = RA(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = RA(J2-2:J2)

         CALL ROTMAT(P, RM)
 
         DO J3 = 1, NRBSITES

            J4          = 3*((J1-1)*NRBSITES + J3)
            XA(J4-2:J4) = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))

         ENDDO

         J2   = 3*J1
         R(:) = RB(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = RB(J2-2:J2)
 
         CALL ROTMAT(P, RM)

         DO J3 = 1, NRBSITES

            J4          = 3*((J1-1)*NRBSITES + J3)
            XB(J4-2:J4) = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))

         ENDDO

      ENDDO

      TEMP=0.0D0
      DO J2=1,NSIZE
            TEMP=TEMP+(XB(J2)-XA(J2))**2
      ENDDO

      END SUBROUTINE UNIAXGETPATHLENGTH

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBCOMMINDIST(RA,RB,NATOMS,DIST,RM,DEBUG)

!     returns squared distance DIST

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG 

      NSIZE = NATOMS/2
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
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4)
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)
      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','newmindist> WARNING - INFO=',INFO,' in DSYEV'

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
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = MINV

      IF (DEBUG) PRINT '(A,G20.10,A,I6)',' rbmindist2> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN

      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE RBCOMMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CHECKDRVTS(X)

      USE COMMONS, ONLY: NATOMS

      IMPLICIT NONE

      INTEGER          :: IVRNO1
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY, RMS, FM, FP, DF, DFN
      LOGICAL          :: GTEST, STEST

      DO IVRNO1 = 1, 3*NATOMS

         WRITE(*, *) IVRNO1

         GTEST    = .FALSE.
         X(IVRNO1) = X(IVRNO1) - 1.D-06
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         FM   = ENERGY
!         WRITE(*, *) FM

         X(IVRNO1) = X(IVRNO1) + 2.D-06
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         FP   = ENERGY
!         WRITE(*, *) FP

         GTEST = .TRUE.
         X(IVRNO1) = X(IVRNO1) - 1.D-06
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         DFN = (FP - FM) / 2.D-06
         DF   = G(IVRNO1)

         WRITE(*, *) DFN
         WRITE(*, *) DF

         IF (ABS(DFN - DF) > 1.D-06) WRITE(*, *) IVRNO1, DFN, DF

      ENDDO

      STOP

      END SUBROUTINE CHECKDRVTS
