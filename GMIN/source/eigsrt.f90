SUBROUTINE EIGSRT(DIAG, V, NDIM, NMAX)
! Performs a selection sort on DIAG and V, such that the elements of DIAG
! (usually the eigenvalues) are ordered from highest to lowest.
!
! Arguments
! ---------
!
! DIAG(in/out):   1D array containing the diagonal elements corresponding to V
! V(in/out):      2D array containing arrays corresponding to DIAG (usually the
!                 hessian)
! NDIM(in):       dimension of DIAG
! NMAX(in):       first dimension of V 
!
   IMPLICIT NONE
   ! Arguments
   INTEGER, INTENT(IN)              :: NDIM, NMAX
   DOUBLE PRECISION, INTENT(INOUT)  :: DIAG(NDIM), V(NMAX, NDIM)
   ! Variables
   DOUBLE PRECISION                 :: PIVOT
   INTEGER                          :: I, J, K

   DO I = 1, NDIM - 1
      K     = I
      PIVOT = DIAG(I)
! Find the maximum value in DIAG(J) and set this as the pivot and swap it.
      DO J = I + 1, NDIM
         IF (DIAG(J) >= PIVOT) THEN
            K     = J
            PIVOT = DIAG(J)
         ENDIF
      END DO
! Swap elements of DIAG
      DIAG(K) = DIAG(I)
      DIAG(I) = PIVOT
! Swap columns of V
      DO J = 1, NMAX
         PIVOT   = V(J, I)
         V(J, I) = V(J, K)
         V(J, K) = PIVOT
      END DO
   END DO
END SUBROUTINE EIGSRT
