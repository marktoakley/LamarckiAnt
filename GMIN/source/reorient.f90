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
!  This routine finds the matrix r that maps a reference structure
!  q0 into a structure q, i.e.
!  q = r q0 for all N 3-dimensional atomic position vectors of q and q0
!
SUBROUTINE REORIENT(NATOMS,Q,Q0,RMAT)
IMPLICIT NONE
INTEGER NATOMS, J1
DOUBLE PRECISION Q(3*NATOMS), Q0(3*NATOMS), RMAT(3,3), AMAT(3,3), CMAT(3,3), WORK(3), SERROR, RQ0(3)
DOUBLE PRECISION, PARAMETER :: EPS=1.0D-10
INTEGER IPIVOT(3), INFO

AMAT(1:3,1:3)=0.0D0
CMAT(1:3,1:3)=0.0D0
DO J1=1,NATOMS
   AMAT(1,1)=AMAT(1,1)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+1)+EPS ! EPS is designed to avoid a singular matrix
   AMAT(1,2)=AMAT(1,2)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+2)+EPS ! for planar systems
   AMAT(1,3)=AMAT(1,3)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+3)+EPS
   AMAT(2,2)=AMAT(2,2)+Q0(3*(J1-1)+2)*Q0(3*(J1-1)+2)+EPS
   AMAT(2,3)=AMAT(2,3)+Q0(3*(J1-1)+2)*Q0(3*(J1-1)+3)+EPS
   AMAT(3,3)=AMAT(3,3)+Q0(3*(J1-1)+3)*Q0(3*(J1-1)+3)+EPS

   CMAT(1,1)=CMAT(1,1)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+1)
   CMAT(1,2)=CMAT(1,2)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+2)
   CMAT(1,3)=CMAT(1,3)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+3)
   CMAT(2,1)=CMAT(2,1)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+1)
   CMAT(2,2)=CMAT(2,2)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+2)
   CMAT(2,3)=CMAT(2,3)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+3)
   CMAT(3,1)=CMAT(3,1)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+1)
   CMAT(3,2)=CMAT(3,2)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+2)
   CMAT(3,3)=CMAT(3,3)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+3)
ENDDO
! symmetrise AMAT
AMAT(2,1)=AMAT(1,2); AMAT(3,1)=AMAT(1,3); AMAT(3,2)=AMAT(2,3)
! invert AMAT
CALL DGETRF(3,3,AMAT,3,IPIVOT,INFO) ! AMAT is modified!
CALL DGETRI(3,AMAT,3,IPIVOT,WORK,3,INFO)
IF (INFO.NE.0) THEN
   PRINT '(A,I6)','ERROR - INFO after DGETRI in reorient=',INFO
   PRINT '(A)','Q:'
   PRINT '(3G20.10)',Q
   PRINT '(A)','Q0:'
   PRINT '(3G20.10)',Q0
   PRINT '(A)','AMAT:'
   PRINT '(3G20.10)',AMAT
   PRINT '(A)','CMAT:'
   PRINT '(3G20.10)',CMAT
   STOP
ENDIF

RMAT(1,1)=CMAT(1,1)*AMAT(1,1)+CMAT(1,2)*AMAT(2,1)+CMAT(1,3)*AMAT(3,1)
RMAT(1,2)=CMAT(1,1)*AMAT(1,2)+CMAT(1,2)*AMAT(2,2)+CMAT(1,3)*AMAT(3,2)
RMAT(1,3)=CMAT(1,1)*AMAT(1,3)+CMAT(1,2)*AMAT(2,3)+CMAT(1,3)*AMAT(3,3)
RMAT(2,1)=CMAT(2,1)*AMAT(1,1)+CMAT(2,2)*AMAT(2,1)+CMAT(2,3)*AMAT(3,1)
RMAT(2,2)=CMAT(2,1)*AMAT(1,2)+CMAT(2,2)*AMAT(2,2)+CMAT(2,3)*AMAT(3,2)
RMAT(2,3)=CMAT(2,1)*AMAT(1,3)+CMAT(2,2)*AMAT(2,3)+CMAT(2,3)*AMAT(3,3)
RMAT(3,1)=CMAT(3,1)*AMAT(1,1)+CMAT(3,2)*AMAT(2,1)+CMAT(3,3)*AMAT(3,1)
RMAT(3,2)=CMAT(3,1)*AMAT(1,2)+CMAT(3,2)*AMAT(2,2)+CMAT(3,3)*AMAT(3,2)
RMAT(3,3)=CMAT(3,1)*AMAT(1,3)+CMAT(3,2)*AMAT(2,3)+CMAT(3,3)*AMAT(3,3)

RETURN 
!
!  Check squared error - debugging only.
!
SERROR=0.0D0
DO J1=1,NATOMS
   RQ0(1)=RMAT(1,1)*Q0(3*(J1-1)+1)+RMAT(1,2)*Q0(3*(J1-1)+2)+RMAT(1,3)*Q0(3*(J1-1)+3)
   RQ0(2)=RMAT(2,1)*Q0(3*(J1-1)+1)+RMAT(2,2)*Q0(3*(J1-1)+2)+RMAT(2,3)*Q0(3*(J1-1)+3)
   RQ0(3)=RMAT(3,1)*Q0(3*(J1-1)+1)+RMAT(3,2)*Q0(3*(J1-1)+2)+RMAT(3,3)*Q0(3*(J1-1)+3)
   SERROR=SERROR+(Q(3*(J1-1)+1)-RQ0(1))**2+(Q(3*(J1-1)+2)-RQ0(2))**2+(Q(3*(J1-1)+3)-RQ0(3))**2
ENDDO
PRINT '(A,G20.10)','error in reorient=',SERROR
IF (SERROR.GT.0.5D0) THEN
   RMAT(1:3,1:3)=-RMAT(1:3,1:3)
   SERROR=0.0D0
   DO J1=1,NATOMS
      RQ0(1)=RMAT(1,1)*Q0(3*(J1-1)+1)+RMAT(1,2)*Q0(3*(J1-1)+2)+RMAT(1,3)*Q0(3*(J1-1)+3)
      RQ0(2)=RMAT(2,1)*Q0(3*(J1-1)+1)+RMAT(2,2)*Q0(3*(J1-1)+2)+RMAT(2,3)*Q0(3*(J1-1)+3)
      RQ0(3)=RMAT(3,1)*Q0(3*(J1-1)+1)+RMAT(3,2)*Q0(3*(J1-1)+2)+RMAT(3,3)*Q0(3*(J1-1)+3)
      SERROR=SERROR+(Q(3*(J1-1)+1)-RQ0(1))**2+(Q(3*(J1-1)+2)-RQ0(2))**2+(Q(3*(J1-1)+3)-RQ0(3))**2
   ENDDO
   IF (SERROR.GT.0.5D0) THEN
      PRINT '(A,G20.10)','error in reorient for inversion=',SERROR
      STOP
   ENDIF
ENDIF

RETURN
END SUBROUTINE REORIENT
