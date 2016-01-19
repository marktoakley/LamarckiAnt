C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C   Finnis-Sinclair potential added by J.A. Elliott 2009
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C*************************************************************************
C
C  Here we calculate the Finnis-Sinclair potential and gradient
C
C*************************************************************************
C
      SUBROUTINE FS (X,V,PG,GRADT)
      USE commons, ONLY : NATOMS, GATOM
      IMPLICIT NONE
      INTEGER J1, J2, j, atomtype
      DOUBLE PRECISION X(3*NATOMS), PG, DIST, PTEMP, V(3*NATOMS),
     1       GA, fsrho(NATOMS), VTEMP, RTEMP, DUMMY,
     2       VTEMP1(NATOMS), VTEMP2(NATOMS), VTEMP3(NATOMS),
     3       d, A, beta, c, c0, c1, c2, c3, dx, dy, dz, min_dist,
     4       cutp, cutr1, cutr2, temp_dist(natoms,natoms)
      LOGICAL GRADT

C Finnis-Sinclar potential parameters
C Phil Mag A 50, 45 (1984)
C parameters for Fe subsequently modified in Erratum
C Phil Mag A 53, 161 (1986)
C (both sets for Fe are included)
C
C bcc metals
C
C d = density cut-off [Angstroms]
C A = binding energy [eV]
C beta = introduces maximum value of phi within the first-nearest-neighbor distance.
C c = two-body interaction cut-off [Angstroms]
C c0,c1,c2 = fitting parameters

      atomtype=gatom

      if (atomtype.eq.1) then

C Vanadium
        d=3.692767d0
        A=2.010637d0
        beta=0.0d0
        c=3.8d0
        c0=-0.8816318d0
        c1=1.4907756d0
        c2=-0.3976370d0

      else if (atomtype.eq.2) then

C Niobium
        d=3.915354d0
        A=3.013789d0
        beta=0.0d0
        c=4.2d0
        c0=-1.5640104d0
        c1=2.0055779d0
        c2=-0.4663764d0

      else if (atomtype.eq.3) then

C Tantalum
        d=4.076980d0
        A=2.591061d0
        beta=0.0d0
        c=4.2d0
        c0=1.2157373d0
        c1=0.0271471d0
        c2=-0.1217350d0

      else if (atomtype.eq.4) then

C Chromium
        d=3.915720d0
        A=1.453418d0
        beta=1.8d0
        c=2.9d0
        c0=29.1429813d0
        c1=-23.3975027d0
        c2=4.7578297d0

      else if (atomtype.eq.5) then

C Molybdenum
        d=4.114825d0
        A=1.887117d0
        beta=0.0d0
        c=3.25d0
        c0=43.4475218d0
        c1=-31.9332978d0
        c2=6.0804249d0

      else if (atomtype.eq.6) then

C Tungsten
        d=4.400224d0
        A=1.896373d0
        beta=0.0d0
        c=3.25d0
        c0=47.1346499d0
        c1=-33.7665655d0
        c2=6.2541999d0

      else if (atomtype.eq.7) then

C Iron - original parameters in Phil Mag A vol. 50 p.45 (1984)
        d=3.699579d0
        A=1.889846d0
        beta=1.8d0
        c=3.4d0
        c0=1.2110601d0
        c1=-0.7510840d0
        c2=0.1380773d0

      else if (atomtype.eq.8) then

C Iron - modified parameters in Erratum, Phil Mag A vol. 53 p.161 (1986)
        d=3.569745d0
        A=1.828905d0
        beta=1.8d0
        c=3.4d0
        c0=1.2371147d0
        c1=-0.3592185d0
        c2=-0.0385607d0

      endif

c If beta is non-zero, set a minimum distance cut-off to avoid rtemp going negative
      if (beta.ne.0d0) then
        min_dist=d*(beta-1.0d0)/beta
      else
        min_dist=0.0d0
      endif

c Calculate the potential and gradient terms

c First, populate density and distance matrices

      do j2=1,natoms                               ! initialise distance and density arrays
        do j1=1,natoms
          temp_dist(j1,j2)=0.0d0
        enddo
        fsrho(j2)=0.0d0
      enddo

      rtemp = 0.0d0                               ! initialise energy accumulators
      ptemp = 0.0d0

      do j1=1,natoms-1                              ! begin outer loop over all atoms except last

        do j2=j1+1,natoms                           ! begin first inner loop over all j2>j1

          temp_dist(j2,j1)=                           ! calc and store distance between atoms j1,j2
     1         dsqrt(( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2         ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3         ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2)

          dist=temp_dist(j2,j1)                       ! store the actual distance used in inner loop
          temp_dist(j1,j2)=dist                       ! impose symmetry on distance matrix

          cutp=dsign(0.5d0,c-dist)+0.5d0                               ! cut-off for 2-body interaction
          ptemp=ptemp
     1          +cutp*(((dist-c)**2)*(c0+dist*(c1+c2*dist)))           ! accumulate two-body potential term

          cutr1=dsign(0.5d0,d-dist)+0.5d0                              ! cut-off for many-body interaction
          cutr2=dsign(0.5d0,dist-min_dist)+0.5d0                       ! minimum distance cut-off
          rtemp=cutr1*cutr2*((dist-d)**2+(beta*(dist-d)**3/d))
          fsrho(j1)=fsrho(j1)+rtemp                   ! accumulate contribution to density matrix
          fsrho(j2)=fsrho(j2)+rtemp                   ! accumulate contribution to density matrix

        enddo                                       ! end inner loop over all j2>j1

      enddo                                         ! end outer loop over all atoms except last

c     Now, sum the potential energy over all atoms

      pg=ptemp                                    ! initialise potential energy

      do j1=1,natoms
        fsrho(j1)=dsqrt(fsrho(j1))                ! square root density
        pg=pg-a*fsrho(j1)                         ! accumulate potential energy
      enddo

c     Calculate gradient terms, if required

      if (gradt) then

        do j1=1,natoms                            ! initialise total gradient terms
          V(3*(J1-1)+1)=0.d0
          V(3*(J1-1)+2)=0.d0
          V(3*(J1-1)+3)=0.d0
          vtemp1(J1)=0.0d0
          vtemp2(J1)=0.0d0
          vtemp3(J1)=0.0d0
        enddo

        do j1=1,natoms-1                                ! begin outer loop over all atoms except last

          dummy=1.0d0/fsrho(j1)                         ! store reciprocal of density element for atom j1

          do j2=j1+1,natoms                             ! begin inner loop over all j2>j1

            dist=temp_dist(j1,j2)                                        ! recall distance from earlier loop

            cutp=dsign(0.5d0,c-dist)+0.5d0                               ! cut-off for 2-body gradient
            vtemp= cutp*(2.0d0*(dist-c)*(c0+dist*(c1+c2*dist))           ! accumulate 2-body gradient term
     1             +((dist-c)**2)*(c1+2.0d0*c2*dist))/dist

            cutr1=dsign(0.5d0,d-dist)+0.5d0                              ! cut-off for many-body gradient
            cutr2=dsign(0.5d0,dist-min_dist)+0.5d0                       ! minimum distance cut-off
            vtemp=vtemp                                                  ! accumulate many-body gradient term
     1            -cutr1*cutr2*((0.5d0*A)*(dummy+1.0d0/fsrho(j2))*
     2            (2.0d0*(dist-d)+3.0d0*(beta/d)*(dist-d)**2))/dist

            dx=(X(3*(J1-1)+1)-X(3*(J2-1)+1))         ! calculate Cartesian components of distance
            dy=(X(3*(J1-1)+2)-X(3*(J2-1)+2))
            dz=(X(3*(J1-1)+3)-X(3*(J2-1)+3))

            vtemp1(j1)=vtemp1(j1)+vtemp*dx           ! accumulate primary gradient components
            vtemp2(j1)=vtemp2(j1)+vtemp*dy
            vtemp3(j1)=vtemp3(j1)+vtemp*dz

            vtemp1(j2)=vtemp1(j2)-vtemp*dx           ! accumulate symmetric gradient components
            vtemp2(j2)=vtemp2(j2)-vtemp*dy
            vtemp3(j2)=vtemp3(j2)-vtemp*dz

          enddo

        enddo

c       Finally, sum the gradient terms over all atoms

        do j1=1,natoms
          V(3*(J1-1)+1)=V(3*(J1-1)+1)+vtemp1(j1)
          V(3*(J1-1)+2)=V(3*(J1-1)+2)+vtemp2(j1)
          V(3*(J1-1)+3)=V(3*(J1-1)+3)+vtemp3(j1)
        enddo

      endif

      RETURN
      END
