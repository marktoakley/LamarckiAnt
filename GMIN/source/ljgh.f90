      SUBROUTINE LJGH(X, G, ENERGY, GTEST, HTEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R6, R12, RIJSQ, DUMMY
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3)
      DOUBLE PRECISION :: DVDR, D2VDR2
      LOGICAL          :: GTEST, HTEST

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (HTEST) HESS(:,:) = 0.D0

      DO J1 = 1, NATOMS  
         J3 = 3*J1
         RI(:)  = X(J3-2:J3)

         DO J2 = J1 + 1, NATOMS
            J4 = 3*J2
            RJ(:)  = X(J4-2:J4)
            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            R6     = R2*R2*R2
            R12    = R6*R6
            
            ENERGY = ENERGY + (R12 - R6)
            IF (GTEST .OR. HTEST) THEN
               DVDR = (-12.D0*R12 + 6.D0*R6)*R2
               G(J3-2:J3)  = G(J3-2:J3) + DVDR*RIJ(:)
               G(J4-2:J4)  = G(J4-2:J4) - DVDR*RIJ(:)
            ENDIF

            IF (HTEST) THEN  
               D2VDR2 = 168.D0*R12*R2*R2 - 48.D0*R6*R2*R2
!     [1] Completely diagonal terms: Same atom, same coordinates
!     xi,xi
               HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2*RIJ(1)*RIJ(1) + DVDR
               HESS(J4-2,J4-2) = HESS(J4-2,J4-2) + D2VDR2*RIJ(1)*RIJ(1) + DVDR
!     yi,yi
               HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2*RIJ(2)*RIJ(2) + DVDR
               HESS(J4-1,J4-1) = HESS(J4-1,J4-1) + D2VDR2*RIJ(2)*RIJ(2) + DVDR
!     zi,zi
               HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2*RIJ(3)*RIJ(3) + DVDR
               HESS(J4,J4)     = HESS(J4,J4)     + D2VDR2*RIJ(3)*RIJ(3) + DVDR

!     [2] Off-diagonal terms on the diagonal block: Same atom, different coordinates
!     xi,yi
               DUMMY = D2VDR2*RIJ(1)*RIJ(2) 
               HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
               HESS(J4-2,J4-1) = HESS(J4-2,J4-1) + DUMMY
!     yi,zi
               DUMMY = D2VDR2*RIJ(2)*RIJ(3) 
               HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
               HESS(J4-1,J4)   = HESS(J4-1,J4) + DUMMY
!     zi,xi
               DUMMY = D2VDR2*RIJ(3)*RIJ(1)
               HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
               HESS(J4-2,J4)   = HESS(J4-2,J4) + DUMMY

!     [3] Diagonal elements on off-diagonal blocks: Different atoms, same coordinate
!     xi,xj
               HESS(J3-2,J4-2) =-D2VDR2*RIJ(1)*RIJ(1) - DVDR
!     yi,yj
               HESS(J3-1,J4-1) =-D2VDR2*RIJ(2)*RIJ(2) - DVDR
!     zi,zj
               HESS(J3,J4)     =-D2VDR2*RIJ(3)*RIJ(3) - DVDR

!     [4] Completely off-diagonal terms: Different atoms and different coordinates
!     xi,yj
               HESS(J3-2,J4-1) =-D2VDR2*RIJ(1)*RIJ(2)
               HESS(J3-1,J4-2) = HESS(J3-2,J4-1) 
!     yi,zj
               HESS(J3-1,J4)   =-D2VDR2*RIJ(2)*RIJ(3)
               HESS(J3,J4-1)   = HESS(J3-1,J4)
!     xi,zj
               HESS(J3-2,J4)   =-D2VDR2*RIJ(3)*RIJ(1)
               HESS(J3,J4-2)   = HESS(J3-2,J4)

            ENDIF
         ENDDO
      ENDDO

      ENERGY    = 4.D0*ENERGY
      G(:)      = 4.D0*G(:)
      IF (.NOT. HTEST) RETURN
      HESS(:,:) = 4.D0*HESS(:,:)

!     Symmetrise Hessian
      DO J1 = 1, 3*NATOMS
         DO J2 = J1+1, 3*NATOMS
            HESS(J2,J1) = HESS(J1,J2)
         ENDDO
      ENDDO

      END SUBROUTINE LJGH
