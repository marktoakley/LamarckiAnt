!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MULTIPERM() 
  !
  ! ds656> Span all the label multipermutations
  !
  USE PORFUNCS
  USE COMMONS, ONLY : NATOMS, NSPECIES, TSTART, NQ, RMS, &
       COORDS, MYUNIT, QALCST, SPANSWAPST, SAVEMULTIMINONLY
  !
  IMPLICIT NONE
  !
  LOGICAL :: MORE, STEP
  INTEGER :: LIST(NATOMS), I, J, JP, K, NQTOT, ITER, QDONE, BRUN
  DOUBLE PRECISION :: POTEL, TIME, SCREENC(3*NATOMS), &
       X0(3*NATOMS)
  !
  COMMON /MYPOT/ POTEL
  COMMON /TOT/ NQTOT
  !
  JP=1 ! single processor...
  !
  ! ============================================================
  ! --- Initialise the list of labels with ascending order -----
  LIST(1:NATOMS) = 0
  K=0
  DO I=1,NSPECIES(0)
     DO J=1,NSPECIES(I)
        K=K+1
        LIST(K) = I
     ENDDO
  ENDDO
  CALL SET_ATOMLISTS(LIST,1)
  ! ============================================================
  !
  ! ============================================================
  ! --- Perform quech zero and store for reference -------------
  IF(.FALSE.) THEN
     WRITE(MYUNIT, '(A)')  'multiperm> Calculating initial energy'
     CALL QUENCH(.FALSE.,JP,ITER,TIME,BRUN,QDONE,SCREENC)
     NQTOT=NQTOT+1
     WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,F11.1)') &
          'Qu ',NQ(JP),' E=', POTEL,' steps=',ITER, &
          ' RMS=',RMS,' t=',TIME-TSTART
     CALL FLUSH(MYUNIT)
     !CALL GSAVEIT_MC(POTEL,COORDS(:,JP),LIST,JP)
  ENDIF
  NQTOT=0
  NQ(JP)=0
  X0(:) = COORDS(:,JP)
  ! ============================================================
  ! 
  ! ============================================================
  ! --- Now cycle through all the multipermutations ------------
  !
  K=0
  MORE = .TRUE.
  !
  DO WHILE(MORE)
     !
     K=K+1
     !
     !IF(K>1) THEN
        CALL QUENCH(.FALSE.,JP,ITER,TIME,BRUN,QDONE,SCREENC)
        NQTOT=NQTOT+1
        NQ(JP)=NQ(JP)+1
        WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,F11.1)') &
             'Qu ',NQ(JP),' E=', POTEL,' steps=',ITER, &
             ' RMS=',RMS,' t=',TIME-TSTART
        CALL FLUSH(MYUNIT)
     !ENDIF
     !
     ! *** Do something with current multipermutation **********
     IF(QALCST) THEN ! Perform biminimisation
        CALL QALCS(JP,ITER,TIME,BRUN,QDONE,SCREENC)
     ELSEIF(SPANSWAPST) THEN ! Count -ve swap gains
        STEP = .FALSE.
        CALL SPAN_SWAPS(JP,ITER,TIME,BRUN,QDONE,SCREENC,STEP)
        IF(SAVEMULTIMINONLY.AND..NOT.STEP) &
             CALL GSAVEIT_MC(POTEL,COORDS(:,JP),LIST,JP)
     ELSE
        CALL GSAVEIT_MC(POTEL,COORDS(:,JP),LIST,JP)    
     ENDIF
     ! *********************************************************
     !
     COORDS(:,JP) = X0(:) ! re-initialise the coordinates
     CALL MULTIPERM_NEXT(NATOMS,LIST,MORE) ! permutation step
     CALL SET_ATOMLISTS(LIST,1) ! update atom lists
     !
  ENDDO
  ! ============================================================
  !
  WRITE(MYUNIT, '(A,I10)')  &
       'multiperm> Finished with permutation count=', K
  CALL FINALQ  
  CALL FINALIO 
  !
  RETURN
  !
END SUBROUTINE MULTIPERM
!
! ================================================================
! ================================================================
! ds656> All the routines that follow have been taken from
! http://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90
! ================================================================
! ================================================================
!
SUBROUTINE MULTIPERM_NEXT( N, A, MORE )
  !
  !***************************************************************
  !
  !! MULTIPERM_NEXT returns the next multipermutation.
  !
  !  Discussion:
  !
  !    A multipermutation is a permutation of objects, 
  !    some of which are identical.
  !
  !    While there are 6 permutations of the distinct 
  !    objects A,B,C, there are only 3 multipermutations of 
  !    the objects A,B,B.
  !
  !    In general, there are N! permutations of N distinct 
  !    objects, but there are N! / ( (M1!) (M2!) ... (MK!) ) 
  !    multipermutations of N objects, in the case where the 
  !    N objects consist of K types, with M1 examples of type 1,
  !    M2 examples of type 2 and so on, and for which objects 
  !    of the same type are indistinguishable.
  !
  !    To begin the computation, the user must set up the first 
  !    multipermutation. To compute ALL possible multipermutations,
  !    this first permutation should list the values in ascending 
  !    order.
  !
  !    The routine will compute, one by one, the next multi-
  !    permutation, in lexicographical order. On the call after 
  !    computing the last multipermutation, the routine will 
  !    return MORE = FALSE (and will reset the multipermutation 
  !    to the FIRST one again.)
  !
  !  Example:
  !
  !    1  1 2 2 3 3
  !    2  1 2 3 2 3
  !    3  1 2 3 3 2
  !    4  1 3 2 2 3
  !    5  1 3 2 3 2
  !    6  1 3 3 2 2
  !    7  2 1 2 3 3
  !    8  2 1 3 2 3
  !    ...
  !   30  3 3 2 2 1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 March 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of items in the 
  !    multipermutation.
  !
  !    Input/output, integer ( kind = 4 ) A(N); on input, the 
  !    current multipermutation. On output, the next 
  !    multipermutation.
  !
  !    Output, logical MORE, is TRUE if the next multipermutation
  !    was computed, or FALSE if no further multipermutations could
  !    be computed.
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(OUT) :: MORE
  INTEGER, INTENT(IN)  :: N 
  INTEGER, INTENT(INOUT) :: A(N)
  INTEGER :: I, M, TEMP
  !
  !  Step 1:
  !  Find M, the last location in A for which A(M) < A(M+1).
  M = 0
  DO I = 1, N-1
     IF(A(I) < A(I+1)) THEN
        M = I
     ENDIF
  ENDDO
  !
  !  Step 2:
  !  If no M was found, we've run out of multipermutations.
  !
  IF(M == 0) THEN
     MORE = .FALSE.
     CALL IVEC_SORT_HEAP_A ( N, A )
     RETURN
  ELSE
     MORE = .TRUE.
  ENDIF
  !
  !  Step 3:
  !  Ascending sort A(M+1:N).
  !
  IF( M + 1 < N ) THEN
     CALL IVEC_SORT_HEAP_A ( N-M, A(M+1:N) )
  ENDIF
  !
  !  Step 4:
  !  Locate the first larger value after A(M).
  !
  I = 1
  DO WHILE ( A(M+I) <= A(M) )
     I = I + 1
  ENDDO
  !
  !  Step 5:
  !  Interchange A(M) and the next larger value.
  !
  TEMP = A(M)
  A(M) = A(M+I)
  A(M+I) = TEMP
  !
  RETURN
  !
END SUBROUTINE MULTIPERM_NEXT
!
SUBROUTINE IVEC_SORT_HEAP_A(N,A)
  !
  !*************************************************************
  !
  !! IVEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
  !
  !  Discussion:
  !
  !    An IVEC is a vector of integer values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Albert Nijenhuis, Herbert Wilf,
  !    Combinatorial Algorithms,
  !    Academic Press, 1978, second edition,
  !    ISBN 0-12-519260-6.
  !
  !  Parameters:
  !
  !    Input, integer N, the number of entries in the array.
  !
  !    Input/output, integer A(N).
  !    On input, the array to be sorted;
  !    On output, the array has been sorted.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(INOUT) :: A(N)
  !
  INTEGER :: N1
  !
  IF ( N <= 1 ) THEN
     RETURN
  END IF
  !
  !  1: Put A into descending heap form.
  !
  CALL IVEC_HEAP_D ( N, A )
  !
  !  2: Sort A.
  !
  !  The largest object in the heap is in A(1).
  !  Move it to position A(N).
  !
  CALL I_SWAP ( A(1), A(N) )
  !
  !  Consider the diminished heap of size N1.
  !
  DO N1 = N - 1, 2, -1
     !
     !  Restore the heap structure of A(1) through A(N1).
     !
    CALL IVEC_HEAP_D ( N1, A )
    !
    !  Take the largest object from A(1) and move it to A(N1).
    !
    CALL I_SWAP ( A(1), A(N1) )
    !
 ENDDO
 !
 RETURN 
 !
END SUBROUTINE IVEC_SORT_HEAP_A
!
SUBROUTINE IVEC_HEAP_D ( N, A )
  !**************************************************************
  !
  !! IVEC_HEAP_D reorders an IVEC into an descending heap.
  !
  !  Discussion:
  !
  !    An IVEC is a vector of integer values.
  !
  !    A descending heap is an array A with the property that, 
  !    for every index J, A(J) >= A(2*J) and A(J) >= A(2*J+1), 
  !    (as long as the indices 2*J and 2*J+1 are legal).
  !
  !                  A(1)
  !                /      \
  !            A(2)         A(3)
  !          /     \        /  \
  !      A(4)       A(5)  A(6) A(7)
  !      /  \       /   \
  !    A(8) A(9) A(10) A(11)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Albert Nijenhuis, Herbert Wilf,
  !    Combinatorial Algorithms,
  !    Academic Press, 1978, second edition,
  !    ISBN 0-12-519260-6.
  !
  !  Parameters:
  !
  !    Input, integer N, the size of the input array.
  !
  !    Input/output, integer A(N).
  !    On input, an unsorted array.
  !    On output, the array has been reordered into a heap.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(INOUT) :: A(N)
  INTEGER :: I, IFREE, KEY, M
  !
  !  Only nodes N/2 down to 1 can be "parent" nodes.
  !
  DO I = N / 2, 1, -1
     !
     !  Copy the value out of the parent node.
     !  Position IFREE is now "open".
     !
     KEY = A(I)
     IFREE = I
     !
     DO
        !
        !  Positions 2*IFREE and 2*IFREE + 1 are the descendants 
        !  of position IFREE. (One or both may not exist because 
        !  they exceed N.)
        !
        M = 2 * IFREE
        !
        !  Does the first position exist?
        !
        IF ( N < M ) THEN
           EXIT
        END IF
        !
        !  Does the second position exist?
        !
        IF ( M + 1 <= N ) THEN
           !
           !  If both positions exist, take the larger of the 
           !  two values, and update M if necessary.
           !
           IF ( A(M) < A(M+1) ) THEN
              M = M + 1
           END IF
           !
        END IF
        !
        !  If the large descendant is larger than KEY, move it up,
        !  and update IFREE, the location of the free position, and
        !  consider the descendants of THIS position.
        !
        IF ( A(M) <= KEY ) THEN
           EXIT
        END IF
        !
        A(IFREE) = A(M)
        IFREE = M
        !
     END DO
     !
     !  Once there is no more shifting to do, KEY moves into the 
     !  free spot IFREE.
     !
     A(IFREE) = KEY
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE IVEC_HEAP_D
!
SUBROUTINE I_SWAP( I, J )
  !
  !********************************************************
  !
  !! I_SWAP switches two I's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 November 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, integer I, J.  On output, the values of 
  !    I and J have been interchanged.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(INOUT) :: I,J
  INTEGER :: K
  !
  K=I
  I=J
  J=K
  !
  RETURN
  !
END SUBROUTINE I_SWAP
