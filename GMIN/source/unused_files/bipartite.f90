!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  Bipartite matching routine to find the closest permutational isomer
!  of one structure from another.
!  The matching is p(i) <--> q(perm(i))
!  Check for "obvious" solution where the minimum distance for each
!  atom defines a matching as well - this will find the answer
!  immediately for aligned permutational isomers.  NOT DONE YET!
!

SUBROUTINE BIPARTITE(N,P,Q,PERM,DIST,WORSTDIST,WORSTRADIUS)
IMPLICIT NONE
INTEGER N, PERM(N), NCONN(2*N+2), NEIGHBOUR(2*N+2,2*N+2), NVERTEX, J1, J2, NPERM, PARENT(2*N+2), J4, OTHERMIN, J3, JMINL, &
        NDIJ, J5, CURRENTMIN, PARENTMIN, LOSE, NDUMMY, NCOUNT, I, MINDIST, MINJ
DOUBLE PRECISION P(3*N), Q(3*N), DIST, WORSTDIST, WORSTRADIUS, WEIGHT(2*N+2,2*N+2), PATHLENGTH(2*N+2), TMPLENGTH, MINLENGTH, &
                 DUMMY
LOGICAL PERMANENT(2*N+2), SIMPLESOL, CHANGE
!
!  Set up initial adjacency lists for the underlying network,
!  which has dimension 2*N+2.
!  Vertex 1 is the source, vertex 2*N+2 is the sink, and the
!  N atoms of structures P and Q occupy vertices 2 to N+1 and
!  N+2 to 2*N+1, respectively.
!
NVERTEX=2*N+2

NCONN(1)=N ! source
DO J1=1,N
   NEIGHBOUR(1,J1)=J1+1
   WEIGHT(1,J1)=0.0D0
ENDDO

SIMPLESOL=.TRUE.
DO J1=1,N  ! structure P
   NCONN(J1+1)=N
   MINDIST=HUGE(MINDIST)
   DO J2=1,N
      NEIGHBOUR(J1+1,J2)=N+J2+1
      WEIGHT(J1+1,J2)=(P(3*(J1-1)+1)-Q(3*(J2-1)+1))**2+(P(3*(J1-1)+2)-Q(3*(J2-1)+2))**2+(P(3*(J1-1)+3)-Q(3*(J2-1)+3))**2
      IF (WEIGHT(J1+1,J2).LT.MINDIST) THEN
         MINDIST=WEIGHT(J1+1,J2)
         MINJ=J2
      ENDIF
   ENDDO
   PERM(J1)=MINJ
   IF (SIMPLESOL) THEN
      DO J2=1,J1-1
         IF (PERM(J1).EQ.PERM(J2)) THEN
            SIMPLESOL=.FALSE.
            EXIT
         ENDIF
      ENDDO
   ENDIF
ENDDO

IF (SIMPLESOL) THEN
   ! PRINT '(A)','simple solution found in bipartite'
   GOTO 100
ENDIF

DO J1=1,N  ! structure Q
   NCONN(J1+N+1)=1
   NEIGHBOUR(J1+N+1,1)=NVERTEX
   WEIGHT(J1+N+1,1)=0.0D0
ENDDO

NCONN(2*N+2)=0 ! sink vertex
!
!  Now run Dijkstra's algorithm to find the shortest path lengths from the source to the
!  other vertices.
!
NDIJ=0
outer: DO
   NDIJ=NDIJ+1
   PATHLENGTH(2:NVERTEX)=HUGE(1.0D0)
   PERMANENT(1:NVERTEX)=.FALSE.
   PARENT(1)=0 
   PATHLENGTH(1)=0.0D0
   NCOUNT=0
   dijkstraloop: DO
      MINLENGTH=HUGE(MINLENGTH)
      DO J1=1,NVERTEX
         IF (PERMANENT(J1)) CYCLE
         IF (PATHLENGTH(J1).LT.MINLENGTH) THEN
            MINLENGTH=PATHLENGTH(J1)
            NDUMMY=J1
         ENDIF
      ENDDO
      PERMANENT(NDUMMY)=.TRUE.
!     CHANGE=.FALSE.
      DO J2=1,NCONN(NDUMMY) ! neighbours of vertex with minimum pathlength
         IF (PATHLENGTH(NEIGHBOUR(NDUMMY,J2)).GT.PATHLENGTH(NDUMMY)+WEIGHT(NDUMMY,J2)) THEN
            PATHLENGTH(NEIGHBOUR(NDUMMY,J2))=PATHLENGTH(NDUMMY)+WEIGHT(NDUMMY,J2)
            PARENT(NEIGHBOUR(NDUMMY,J2))=NDUMMY
!           CHANGE=.TRUE.
         ENDIF
      ENDDO
      NCOUNT=NCOUNT+1
!     PRINT '(A,I6,A,L5)','Dijkstra iteration ',NCOUNT,' change=',CHANGE
      IF (NCOUNT.GE.NVERTEX) EXIT
   ENDDO dijkstraloop
!
! Adjust edge weights
!
   DO J1=1,NVERTEX
      DO J2=1,NCONN(J1)
         WEIGHT(J1,J2)=WEIGHT(J1,J2)+PATHLENGTH(J1)-PATHLENGTH(NEIGHBOUR(J1,J2))
      ENDDO
   ENDDO
!
! Reverse shortest path. NCONN should stay the same for vertices 2 to 2*N+1,
! and change to NATOMS-NDIJ and NDIJ for the source and sink on cycle NDIJ.
!
   CURRENTMIN=NVERTEX
   PARENTMIN=PARENT(CURRENTMIN)
   NCONN(CURRENTMIN)=NCONN(CURRENTMIN)+1
   NEIGHBOUR(CURRENTMIN,NCONN(CURRENTMIN))=PARENTMIN
   DO J1=1,NCONN(PARENTMIN)
      IF (NEIGHBOUR(PARENTMIN,J1).EQ.CURRENTMIN) THEN
         WEIGHT(CURRENTMIN,NCONN(CURRENTMIN))=WEIGHT(PARENTMIN,J1)
         LOSE=J1
         EXIT
      ENDIF
   ENDDO
   DO
     CURRENTMIN=PARENTMIN
     PARENTMIN=PARENT(CURRENTMIN)
     DO J1=1,NCONN(PARENTMIN)
        IF (NEIGHBOUR(PARENTMIN,J1).EQ.CURRENTMIN) THEN
           WEIGHT(CURRENTMIN,LOSE)=WEIGHT(PARENTMIN,J1)
           NEIGHBOUR(CURRENTMIN,LOSE)=PARENTMIN
           LOSE=J1
           EXIT
        ENDIF
     ENDDO
     IF (PARENTMIN.EQ.1) EXIT
   ENDDO
   DO J2=LOSE,NCONN(1)-1
      WEIGHT(1,J2)=WEIGHT(1,J2+1)
      NEIGHBOUR(1,J2)=NEIGHBOUR(1,J2+1)
   ENDDO
   NCONN(1)=NCONN(1)-1
   IF (NCONN(1).EQ.0) EXIT outer
ENDDO outer

DO J1=1,N
   PERM(NEIGHBOUR(N+1+J1,1)-1)=J1
ENDDO

100 CONTINUE
WORSTDIST=-1.0D0
DIST=0.0D0
DO I=1,N
   DUMMY=(p(3*(i-1)+1)-q(3*(perm(i)-1)+1))**2+(p(3*(i-1)+2)-q(3*(perm(i)-1)+2))**2+(p(3*(i-1)+3)-q(3*(perm(i)-1)+3))**2
   DIST=DIST+DUMMY
   IF (DUMMY.GT.WORSTDIST) THEN
      WORSTDIST=DUMMY
      WORSTRADIUS=p(3*(i-1)+1)**2+p(3*(i-1)+2)**2+p(3*(i-1)+3)**2
   ENDIF
ENDDO
WORSTDIST=SQRT(WORSTDIST)
WORSTRADIUS=MAX(SQRT(WORSTRADIUS),1.0D0)

END SUBROUTINE BIPARTITE
