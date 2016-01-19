      SUBROUTINE GEM (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, GEMRC, BOXLX

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, DVDR, ABSR, R2, R4, RC4, EXPMR4, CA, CB, RCUTSQ
      DOUBLE PRECISION :: RIJ(3) 
      LOGICAL          :: GTEST

      RC4 = GEMRC**4.D0
      CA  = - DEXP(-RC4)*(4.D0*RC4 + 1.D0)
      CB  = 4.D0*GEMRC**3.D0*DEXP(-RC4)
      RCUTSQ = GEMRC*GEMRC

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

      DO I = 1, NATOMS - 1 

         J1 = 3*I

         DO J = I + 1, NATOMS

            J2 = 3*J

            RIJ(:) = X(J1-2:J1) - X(J2-2:J2)
            RIJ(:) = RIJ(:) - ANINT(RIJ(:)/BOXLX)*BOXLX
            R2     = DOT_PRODUCT(RIJ(:),RIJ(:))

            IF (R2 < RCUTSQ) THEN

               R4     = R2*R2
               EXPMR4 = DEXP(-R4)
               ABSR   = DSQRT(R2)
          
               ENERGY = ENERGY + EXPMR4 + CA + CB*ABSR

               IF (R2 == 0.D0) THEN
                  DVDR = 0.D0
               ELSE 
                  DVDR   = -4.D0*R2*DEXP(-R4) + CB/ABSR
               ENDIF
 
               IF (GTEST) THEN

                  G(J1-2:J1) = G(J1-2:J1) + DVDR*RIJ(:)
                  G(J2-2:J2) = G(J2-2:J2) - DVDR*RIJ(:)

               ENDIF

            ENDIF

         ENDDO

      ENDDO

      END SUBROUTINE GEM 
