!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2010 David J. Wales
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
! 
!  These subroutines are taken mostly from GROUPROTATION subroutines. 

SUBROUTINE SAMPLE_TORSIONS(angles,lnJ,torsdelta,a4list,groupsindex,JP)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION :: Gbias(8,8),Gprime(8,8),angles(8),weight,lnJ,torsdelta
      INTEGER :: JP,a4list(3),groupsindex(8),I1,I2
! Some helpful parameters

      ! define list of atoms at boundary of rotation
      a4list(:) = (/DIHEDRALGROUPAXIS(groupsindex(8),2),DIHEDRALGROUPAXIS(groupsindex(8),3),&
     &    DIHEDRALGROUPAXIS(groupsindex(8),3)+1/)

      ! calculate matrix G 
      call calculate_G(Gbias,COORDS(:,JP),groupsindex,a4list,JP)

      ! calculate torsion rotations via cholesky decomp. of biasing matrix
      call sample_rotations(angles,lnJ,Gbias)
       
      ! calculate displacement of vector of angles:
      torsdelta = sqrt(dot_product(angles(:),angles(:)))

END SUBROUTINE SAMPLE_TORSIONS

subroutine calculate_G(Gbias,C,groupsindex,a4list,JP)

   use commons
   implicit none

   double precision :: Gbias(8,8),Gbias2(8,8)
   double precision :: C(3*NATOMS),phi,Mdrdphi(3,8,3),drdphi(3),A(3),B(3),ex(3),ey(3),ez(3),Coordsystem(3,8,9),Rmat(3,3),rhoD
   integer :: a1,a2,a3,a4list(3),groupsindex(8),f1,g1,gi,mygroup,I,atomB3,JP

   DO I=1,3 
      ! calculate rhoD,phi, for atom I
      DO gi=1,8
         mygroup=groupsindex(gi)

         call calc_cylindrical_coords(drdphi(:),DIHEDRALGROUPAXIS(mygroup,1),DIHEDRALGROUPAXIS(mygroup,2),&
     &       DIHEDRALGROUPAXIS(mygroup,3),a4list(I),C)

         ! take care of special cases:
         IF((gi.eq.7).AND.(I.EQ.1)) THEN
            drdphi(1) = 0.0
            drdphi(2) = 0.0
            drdphi(3) = 0.0
         ELSE IF(gi==8.AND.((I.EQ.1).OR.(I.EQ.2))) THEN
            drdphi(1) = 0.0
            drdphi(2) = 0.0
            drdphi(3) = 0.0
         ENDIF

         Mdrdphi(I,gi,1:3) = drdphi(1:3)

      ENDDO
   ENDDO
      
   ! calculate dot product
   do f1=1,8
      do g1=1,8
         Gbias(f1,g1) = 0.0 
         DO I=1,3
            Gbias(f1,g1) = Gbias(f1,g1) + dot_product(Mdrdphi(I,g1,1:3),Mdrdphi(I,f1,1:3))
         ENDDO
      ENDDO
   ENDDO


end subroutine calculate_G

subroutine calc_cylindrical_coords(drdphi,a1,a2,a3,a4,C)

   use commons
   implicit none

   double precision :: C(3*NATOMS),rhoD,phi,n23(3),n34(3),R(3),lenR,theta,ex(3),ey(3),ez(3),drdphi(3)
   integer :: a1,a2,a3,a4

   ! dihedral calculation
   call my_calc_dihedral(phi,a1,a2,a3,a4,C)
   

   ! rhoD calculation
   ! unit vector along bond axis of atoms 2-3 
   n23(:) = C(3*(a3-1)+1:3*(a3-1)+3) - C(3*(a2-1)+1:3*(a2-1)+3)
   n23(:) = n23(:) / sqrt(DOT_PRODUCT(n23,n23))
   R(:) = C(3*(a4-1)+1:3*(a4-1)+3) - C(3*(a3-1)+1:3*(a3-1)+3)
   lenR = sqrt(DOT_PRODUCT(R(:),R(:)))
   n34(:) = R(:) / lenR
   theta = acos(dot_product(n23,n34))

   rhoD = lenR * sin(theta)
   
   ey(:) = n23(:)
   call cross_product(ex(:),ey(:),n34(:))
   !write(*,'(A,4F20.10,2I2.1)') "e> ",ex(1),ey(1),n34(1),lenR,a3,a4
   ex(:) = ex(:) * 1./(sqrt(dot_product(ex(:),ex(:))))
   call cross_product(ez(:),ex(:),ey(:))

   call cross_product(drdphi(:),n34(:),n23)

end subroutine calc_cylindrical_coords

subroutine get_Abias(Abias,Gbias)

   use commons, only: bgsb1,bgsb2
   implicit none

   double precision :: Gbias(8,8),Abias(8,8)
   integer :: a1,a2,a3,a4,a,b,info,Lwork

   !bgsb1 = 30.0
   !b2 = 0.0
   !bgsb2 = 1.0

   Abias(:,:) = 0.5*bgsb1*bgsb2*Gbias(:,:)
   do a1=1,8
      Abias(a1,a1) = Abias(a1,a1) + 0.5*bgsb1
   enddo

end subroutine get_Abias

subroutine sample_rotations(angles,lnJ,Gbias)

   use commons
   USE random_normal_module

   implicit none

   double precision :: angles(8),Gbias(8,8),psi(8),Workev(8),EvalsG(8),Abias(8,8),b1,b2,weight,detAfac,lnJ,lnJtest,dummy1(8),anglesnormal(8)
   double precision :: AP(36),Uinv(8,8)
   integer :: a1,a2,a3,a4,a,b,info,Lwork,I,J,JC

   ! calculate standard normal values psi
   do a=1,8
      angles(a) = sqrt(0.5)*RANDOM_NORMAL()
   enddo
   anglesNormal(:) = angles(:)

   call get_Abias(Abias,Gbias)

   ! do cholesky decomposition of Abias, using "U"pper triangle 
   call dpotrf('U',8,Abias,8,info)
   if (info.ne.0) then
      write(myunit,'(A)') "error: info != 0 "
      STOP
   endif

! calculate weight factor of move
   !calculate determinant factor of Abias: detAfac = (det A)^(1/2)
   ! here Abias below is now the cholesky decomposed matrix L^T , where A = L L^T
   detAfac = 1.0
   DO I=1,8
      detAfac = detAfac * Abias(I,I)
   ENDDO

   lnJ = -dot_product(angles(:),angles(:)) + log(detAfac)

   ! express Upper triangular matrix U of Abias in packed format
     JC = 1                     
     DO 2 J = 1, 8
        DO 1 I = 1, J          
           AP(JC+I-1) = Abias(I,J) 
   1    CONTINUE               
        JC = JC + J            
   2 CONTINUE                   
  
   ! Find inverse of U
   call DTPTRI('U','N',8,AP,info)
   !write(*,'(A,I2.1)') "info ", info
   if (info.ne.0) then
      write(myunit,'(A)') "error: info != 0 "
      STOP
   endif

   ! write inverse in standard (i,j) form
   DO I=1,8
      Do J=1,8
         Uinv(i,j) = 0.0
      ENDDO
   ENDDO

     JC = 1                     
     DO 4 J = 1, 8
        DO 3 I = 1, J          
           Uinv(I,J) = Uinv(I,J) + AP(JC+I-1)
   3    CONTINUE               
        JC = JC + J            
   4 CONTINUE                   
   
   ! solve linear system of equations using U^-1
   angles(:) = MATMUL(Uinv,anglesnormal(:))

   ! sanity check
   !call calculate_weighting_factor(lnJtest,Gbias,angles(:))

end subroutine sample_rotations


subroutine calculate_weighting_factor(lnJ,Gbias,dphi)

   implicit none

   double precision :: weight,Gbias(8,8),Abias(8,8),dphi(8),psi(8),detAfac,lnJ,dummy(8,8),AbiasT(8,8)
   integer :: info,I,J

   call get_Abias(Abias,Gbias)

   !do cholesky decomposition of Abias
   call dpotrf('U',8,Abias,8,info)
   !call F07FDF('L',8,Abias,8,info)
   if (info.ne.0) then
      write(*,'(A)') "error: info != 0 "
      STOP
   endif
   
   DO I=1,8
      DO J=I+1,8
         Abias(J,I) = 0.0
      ENDDO
   ENDDO
   
   !calculate determinant factor of Abias: detAfac = (det A)^(1/2)
   detAfac = 1.0
   DO I=1,8
      detAfac = detAfac * Abias(I,I)
   ENDDO

   !evaluate standard normal variates psi from dphi
   psi(:) = MATMUL(Abias(:,:),dphi(:))

   !calculate exponential factor
   ! exp(ln J) = weight
   weight = exp(-dot_product(psi(:),psi(:)))

   ! bug fix here:
   lnJ = -dot_product(psi(:),psi(:)) + log(detAfac)
   !lnJ = dot_product(psi(:),psi(:)) - log(detAfac)

   weight = weight * detAfac


end subroutine calculate_weighting_factor
