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

! NEW ROUTINE TO DRIVE DIHEDRALROTATION MOVES
! JP is the parallel run ID so that only the appropriate coordinates are
! altered during parallel runs.
SUBROUTINE DIHEDRALROTSTEP(deltalnJ,JP,BETA,IMCSTEP)
      USE commons
      USE modcharmm, only: CHNEIGHBOURT
      USE random_normal_module
      IMPLICIT NONE
      DOUBLE PRECISION :: DPRAND, PI, TWOPI, angle, degreeangle, MYRANDOM, Coordinates(NATOMS),angle0,POTEL,GRAD(3*NATOMS),thetavar,torsions(8),BETA(0:NPAR-1)
      DOUBLE PRECISION :: a4dist,COORDSOLD(3*NATOMS),weight,weight_perturbed,Gbias(8,8),lnJ,lnJpert,torsdelta,deltalnJ,potelold,NDIHEDRALBBGROUPS,imcstep
      INTEGER :: I1,RESC,RES1,RES2,Gstart,Gstop,Numrot,JP,Nresiduemin,Nresiduemax,a,nangle,groupsindex(8),a4list(3),J1
! Some helpful parameters
      PI=ATAN(1.0D0)*4
      TWOPI=2.0D0*PI

! For each group....      

      COORDSOLD(:) = COORDS(:,JP)
      CALL POTENTIAL(COORDS(:,jp),GRAD,POTELold,.TRUE.,.FALSE.)
      !The BGS move is currently selected under the keyword BGSMOVE
      ! currently implemented with probability 0.5, with a dihedralrotation
      ! (pivot) attempted if not

      IF (BGSMOVE.and.(DPRAND().lt.pselectBGS)) THEN
         Nresiduemin = 8
         Nresiduemax = 8
         Numrot = 8
         !NDIHEDRAL_BB_GROUPS=NDIHEDRALGROUPS

         ! choose a random starting residue
         Gstart = 1 + DPRAND()*(NDIHEDRAL_BB_GROUPS-Numrot)
         Gstop = Numrot + Gstart - 1
         !WRITE(MYUNIT,*) 'BGSMOVE> Rotating dihedrals ',Gstart, "to", Gstop,Numrot

         DO I1=1,8
            groupsindex(I1) = Gstart + I1-1
         ENDDO

         ! sample torsion rotations from groups 
         call sample_torsions(torsions,lnJ,torsDelta,a4list,groupsindex,JP)

         !do I1=1,8
         !   write(myunit,'(A,I2.1,F20.10)') "torsions> ", I1, torsions(I1)*180./3.14
         !enddo


         DO I1=Gstart,Gstop
            ! calculate current value of dihedral
            call my_calc_dihedral(angle0,DIHEDRALGROUPAXIS(I1,1),DIHEDRALGROUPAXIS(I1,2),DIHEDRALGROUPAXIS(I1,3),DIHEDRALGROUPAXIS(I1,4),COORDS(:,JP))
            !write(myunit,"(2A,I5.1,F20.10)") "dihedral values> ", TRIM(ADJUSTL(DIHEDRALGROUPNAMES(I1))),I1, angle0

            ! choose an angle to move to
            angle = torsions(I1-Gstart+1)
            degreeangle=angle*(180./PI)

            ! Print some info to GMIN_out for the user
            !WRITE(MYUNIT,'(3A,F20.10,A)') 'Helixmove> Rotating group ',TRIM(ADJUSTL(DIHEDRALGROUPNAMES(I1))),' by ',degreeangle, " degrees"

            ! Call the rotation subroutine
            CALL DIHEDRALROTATION(angle,DIHEDRALGROUPAXIS(I1,1),DIHEDRALGROUPAXIS(I1,2),DIHEDRALGROUPAXIS(I1,3),DIHEDRALGROUPAXIS(I1,4),DIHEDRALGROUPS(I1,:),COORDS(:,JP))

            ! compute final dihedral value
            !call my_calc_dihedral(angle0,DIHEDRALGROUPAXIS(I1,1),DIHEDRALGROUPAXIS(I1,2),DIHEDRALGROUPAXIS(I1,3),DIHEDRALGROUPAXIS(I1,4),COORDS(:,JP))
            !WRITE(MYUNIT,*) 'DIHEDRALROTATION> Final angle is ',angle0*(180./pi),'target',degreeangle
         ENDDO 

         ! calculate jacobian weight of final state
         call calculate_G(Gbias,Coords(:,JP),groupsindex,a4list,JP)
         call calculate_weighting_factor(lnJpert,Gbias,-torsions)

         ! calculate weights of F / R moves
         !write(*,'(A,2F20.10,A,F20.10)') "Jacobian weight of trial move: ",lnJpert,lnJ," wprime/w ", exp(-(lnJpert-lnJ))
         deltalnJ = lnJpert-lnJ

         CALL POTENTIAL(COORDS(:,jp),GRAD,POTEL,.TRUE.,.FALSE.)

         ! compute displacement of a4 coords by rotation
         a4dist = 0.0
         DO J1=1,3
            a4dist = a4dist + (COORDS(3*(a4list(J1)-1)+1,JP)-COORDSOLD(3*(a4list(J1)-1)+1))**2
            a4dist = a4dist + (COORDS(3*(a4list(J1)-1)+2,JP)-COORDSOLD(3*(a4list(J1)-1)+2))**2
            a4dist = a4dist + (coords(3*(a4list(j1)-1)+3,jp)-coordsold(3*(a4list(j1)-1)+3))**2
         enddo
         a4dist = sqrt(a4dist)
              
!            write(myunit,'(A,2I6.1,2F20.10,A,F20.10,A,F20.10,A,3F20.10)') "dihedralroation>",Gstart,Gstop,POTEL,POTELOLD, " atomdist: ",a4dist," stepsize (in torsions) ",&
!     &       torsdelta, " dE/T, -deltalnJ, dE/T-deltalnJ", (POTEL-POTELold)*BETA(MYNODE),-deltalnJ,  (POTEL-POTELold)*BETA(MYNODE)-deltalnJ

         IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) THEN
            write(myunit,'(A,2I6.1,5F20.10)') "dihedralroation>",Gstart,Gstop,POTEL,POTELOLD,a4dist,&
     &       torsdelta, (POTEL-POTELold)*BETA(MYNODE)-deltalnJ
         ENDIF

      ELSE 
         ! do standard independent dihedral rotations 
         DO I1=1,NDIHEDRALGROUPS
            !IF (DIHEDRALGROUPPSELECT(I1).GE.0.0) THEN
            IF (DIHEDRALGROUPPSELECT(I1).GE.DPRAND()) THEN
               angle=twopi * (dihedralgroupscaling(i1) * (dprand()-0.5))
               degreeangle = angle*(180./pi)
               call my_calc_dihedral(angle0,dihedralgroupaxis(i1,1),dihedralgroupaxis(i1,2),dihedralgroupaxis(i1,3),dihedralgroupaxis(i1,4),coords(:,jp))
            ! call the rotation subroutine
               call dihedralrotation(angle,dihedralgroupaxis(i1,1),dihedralgroupaxis(i1,2),dihedralgroupaxis(i1,3),dihedralgroupaxis(i1,4),dihedralgroups(i1,:),coords(:,jp))
            ! print some into to gmin_out for the user
               CALL POTENTIAL(COORDS(:,jp),GRAD,POTEL,.TRUE.,.FALSE.)
               IF (MOD(IMCSTEP-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) THEN
                  write(myunit,'(A,A,A,F20.10,A,F20.10,A,2F20.10,F10.1)') 'dihedralrotation> rotating group ',trim(adjustl(dihedralgroupnames(i1))),' by ',degreeangle,' from '&
     &                  ,angle0*180/pi," new energy ", POTEL,POTELold,IMCSTEP
               ENDIF
            endif
         enddo
      endif


end subroutine dihedralrotstep 

subroutine dihedralrotation(angle,a0,a1,a2,a3,atomingroup,stepcoords)
   USE commons
   IMPLICIT NONE
   INTEGER :: A1, A2, A0,A3,I1
   DOUBLE PRECISION :: BVECTOR(3), LENGTH, ANGLE, DUMMYMAT(3,3)=0.0D0, ROTMAT(3,3)
   DOUBLE PRECISION :: GROUPATOM(3), GROUPATOMROT(3),STEPCOORDS(3*NATOMS),v1(3),v2(3),v3(3),alpha
   LOGICAL :: ATOMINGROUP(NATOMS)

    !WRITE(*,*) "DIHEDRALROTATION> dihedral angle value to rot", angle*180/3.1415
! STEP 1
! Produce notmalised bond vector corresponding to the rotation axis
! A1 and A2 are the atoms defining this vector
   BVECTOR(1)=STEPCOORDS(3*A2-2)-STEPCOORDS(3*A1-2)
   BVECTOR(2)=STEPCOORDS(3*A2-1)-STEPCOORDS(3*A1-1)
   BVECTOR(3)=STEPCOORDS(3*A2  )-STEPCOORDS(3*A1  )
! Find length   
   LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
   BVECTOR(1)=BVECTOR(1)/LENGTH
   BVECTOR(2)=BVECTOR(2)/LENGTH
   BVECTOR(3)=BVECTOR(3)/LENGTH
! STEP 2
! Scale this vector so its length is the rotation to be done (in radians)
   BVECTOR(1)=BVECTOR(1)*ANGLE
   BVECTOR(2)=BVECTOR(2)*ANGLE
   BVECTOR(3)=BVECTOR(3)*ANGLE
! STEP 3
! Get the rotation matrix for this vector axis from RMDRVT
! Interface:
! SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)
! P is an un-normalised vector you wish to rotate around. Its length equals the desired rotation in radians
! RM will return the 3x3 rotation matrix
! DRM1-3 are derivative matricies, not needed here
! GTEST is also not needed so set to .FALSE.
   CALL RMDRVT(BVECTOR,ROTMAT,DUMMYMAT,DUMMYMAT,DUMMYMAT,.FALSE.)
! STEP 4
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
   DO I1=1,NATOMS
      IF (ATOMINGROUP(I1)) THEN
         GROUPATOM(1)=STEPCOORDS(3*I1-2)-STEPCOORDS(3*A1-2)
         GROUPATOM(2)=STEPCOORDS(3*I1-1)-STEPCOORDS(3*A1-1)
         GROUPATOM(3)=STEPCOORDS(3*I1  )-STEPCOORDS(3*A1  )
! Apply the rotation matrix
         GROUPATOMROT=MATMUL(ROTMAT,GROUPATOM)
! Translate back to the origin and copy to COORDS
         STEPCOORDS(3*I1-2)=GROUPATOMROT(1)+STEPCOORDS(3*A1-2)
         STEPCOORDS(3*I1-1)=GROUPATOMROT(2)+STEPCOORDS(3*A1-1)
         STEPCOORDS(3*I1  )=GROUPATOMROT(3)+STEPCOORDS(3*A1  )
      ENDIF
   ENDDO

   !call my_calc_dihedral(alpha,A0,A1,A2,A3,STEPCOORDS(:))
   !WRITE(*,*) "DIHEDRALROTATION> final dihedral angle value  ", alpha*180/3.1415

END SUBROUTINE DIHEDRALROTATION 

SUBROUTINE MY_CALC_DIHEDRAL(alpha,a0,a1,a2,a3,c)

   ! calculates dihedral angle alpha defined by four atoms a0,a1,a2,a3 for
   ! configuration c
   USE commons
   IMPLICIT NONE
   DOUBLE PRECISION v1(3),v2(3),v3(3),r(3),n1(3),n2(3),c(3*NATOMS)
   DOUBLE PRECISION v2norm(3),rdotv2norm,alpha,norm,n1dotn2,normr
   INTEGER a0,a1,a2,a3

   !WRITE(*,*), "my_dihed initial vectors> ", v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),v3(1),v3(2),v3(3)

   ! calculate current dihedral angle corresponding to particular rotation
   v1(1)=c(3*a1-2)-c(3*a0-2)
   v1(2)=c(3*a1-1)-c(3*a0-1)
   v1(3)=c(3*a1)-c(3*a0  )
   v2(1)=c(3*a2-2)-c(3*a1-2)
   v2(2)=c(3*a2-1)-c(3*a1-1)
   v2(3)=c(3*a2)-c(3*a1)
   v3(1)=c(3*a3-2)-c(3*a2-2)
   v3(2)=c(3*a3-1)-c(3*a2-1)
   v3(3)=c(3*a3  )-c(3*a2)

   CALL cross_product(n1,v1,v2)
   CALL cross_product(n2,v2,v3)

   CALL cross_product(r,n1,n2)

   n1dotn2 = DOT_PRODUCT(n1,n2) 
   norm = DOT_PRODUCT(v2,v2) 
   norm = SQRT(norm)

   normr = DOT_PRODUCT(r,r)
   normr = SQRT(normr)
   !v2norm(:) = v2norm(:) / norm
   v2norm(1) = v2(1) / norm
   v2norm(2) = v2(2) / norm
   v2norm(3) = v2(3) / norm
   rdotv2norm = DOT_PRODUCT(r,v2norm)
   alpha = ATAN2(rdotv2norm,n1dotn2)
   !WRITE(*,*), "my_dihed> ", r(1),r(2),r(3),rdotv2norm,normr

   
END SUBROUTINE MY_CALC_DIHEDRAL

SUBROUTINE cross_product(c,u,v)

   ! calculates cross product c of vectors u and v (3 dimensions)
   IMPLICIT NONE

   DOUBLE PRECISION u(3),v(3), c(3)

   c(1) = u(2)*v(3)-u(3)*v(2)
   c(2) = u(3)*v(1)-u(1)*v(3)
   c(3) = u(1)*v(2)-u(2)*v(1)


END SUBROUTINE cross_product

SUBROUTINE sample_angle_value(ngroup,alpha)

   USE commons
   use random_normal_module
   IMPLICIT NONE
   DOUBLE PRECISION mean,sigma,pi,alpha
   INTEGER ngroup

   PI=ATAN(1.0D0)*4

   IF (ANGLETYPE(ngroup).EQ."phi") THEN
      mean = PHI0
      sigma = PHIk
   ELSE IF (ANGLETYPE(ngroup).EQ."psi") THEN
      mean = PSI0
      sigma = PSIk
   ELSE 
      WRITE(*,*) "angletype ", ANGLETYPE(ngroup), "for group",ngroup,"not understood"
      STOP
   ENDIF

   alpha = PI*(mean + sigma*random_normal())

END SUBROUTINE sample_angle_value
