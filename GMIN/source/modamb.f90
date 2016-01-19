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
      MODULE MODAMBER 
      IMPLICIT NONE
      SAVE

      DOUBLE PRECISION  top,bottom,ax,bx,cx,ay,by,cy,az,bz,cz,mu,lambda,numer,denom,dx,dy,dz
      DOUBLE PRECISION  totenergy,benergy,tenergy,penergy,answer,impenergy,qenergy,vdwenergy
      INTEGER  a,b,c,ios,atoms,d,e,f,g,h,i,j,k,l,colin,fish,ans,rings,ang,imp,count
      INTEGER  match,bean,bondnumber,t, NDIHEDRALS
!      INTEGER  arraysize
!      PARAMETER  (arraysize=NATOMS)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DATOM1, DATOM2 !mxatms
      INTEGER loopatom(20)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  bondarray   ! (arraysize*4,2)
      DOUBLE PRECISION TINY
      PARAMETER  (TINY=1.0D-12)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x, y, z, q    !(arraysize)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r, vdwa, vdwb   !(arraysize,arraysize)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mass   !(arraysize)
      LOGICAL WATERSTEP
      DOUBLE PRECISION xbar, ybar, zbar

      INTEGER qlistcount, rlistcount, listupdate
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  qlist, rlist   !(arraysize**2,2)
      LOGICAL AMCUT
      DOUBLE PRECISION  QCUTOFF, RCUTOFF, REALQCUTOFF, REALRCUTOFF

      DOUBLE PRECISION gentorsparams(100,6),spectorsparams(100,14),specimpparams(15,7)
      DOUBLE PRECISION midimpparams(4,6),genimpparams(15,5)
  
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  aa1, aa2, aa3   !(arraysize*5)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  theta   !(arraysize*5)

      INTEGER, ALLOCATABLE, DIMENSION(:) ::  da1, da2, da3, da4   !(arraysize*5)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  dphi, did, dvn, dvn2   !(arraysize*5)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  dvn3, ddelta, ddelta2, ddelta3  !(arraysize*5)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  dn, dn2, dn3   !(arraysize*5)

      INTEGER, ALLOCATABLE, DIMENSION(:) ::  ia1, ia2, ia3, ia4  !(arraysize*5)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  iphi, ivn, idelta, in1 !(arraysize*5)

      DOUBLE PRECISION  qx(3),qy(3),qz(3)
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  atnum, bondedto, type ! (arraysize)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bonds   !(arraysize,arraysize)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: one_four, one_three   !(arraysize,arraysize)
      CHARACTER(LEN=2)  typechb,typechc,typechd,typeche
      CHARACTER(LEN=20)  filename
      CHARACTER, ALLOCATABLE, DIMENSION(:) ::  label
      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: typech

! Variables for energy calculation

      CHARACTER(LEN=10)  check
      DOUBLE PRECISION dielec
      PARAMETER  (dielec=1)
      LOGICAL FAKEWATER, MGBWATER
      LOGICAL CART
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: alpha  !(arraysize)
      DOUBLE PRECISION rstar,epsilon,dparam
      DOUBLE PRECISION vdwr(42),vdwe(42)
      DOUBLE PRECISION kr(42,42),ro(42,42)

      DOUBLE PRECISION kt(42,42,42),to(42,42,42),degto(42,42,42)
      DOUBLE PRECISION PK,PHASE,PK2,PK3,PHASE2,PHASE3,IPK,IPHASE

      INTEGER PN,IDIVF,PN2,PN3,IPN

      INTEGER ambercount, ambergcount

! Variables for derivative calculation

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dbondEbydx, dbondEbydy, dbondEbydz, dangEbydx,dangEbydy, dangEbydz,& 
                                                   & dtorsEbydx, dtorsEbydy, dtorsEbydz, dvdwEbydx, dvdwEbydy, dvdwEbydz,&
                                                   & dEbydx, dEbydy, dEbydz, dimpEbydx, dimpEbydy, dimpEbydz,&
                                                   & dqEbydx, dqEbydy, dqEbydz   ! arraysize
      DOUBLE PRECISION arms

      INTEGER nbin1, nbin2, nbin3, nbin4, nbin5, nbin6, nbin7
      LOGICAL BIN
      LOGICAL AMBERSEED, CAP
      DOUBLE PRECISION  seedphi(50), seedpsi(50), actualphi(50), actualpsi(50), actualomega(50)
      INTEGER phiatom1(50), phiatom2(50), psiatom1(50), psiatom2(50), res
     
      INTEGER,ALLOCATABLE,DIMENSION(:) :: chiral
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: chiralarray
 
! Variables for second derivative calculation

!      DOUBLE PRECISION bondhell(3*arraysize,3*arraysize), anglehell(3*arraysize,3*arraysize), torshell(3*arraysize,3*arraysize), 
!     1                 imphell(3*arraysize,3*arraysize), qhell(3*arraysize,3*arraysize), vdwhell(3*arraysize,3*arraysize) 
!     2                 ,hell(3*arraysize,3*arraysize)
!      COMMON /HELL/ bondhell, anglehell, torshell, imphell, qhell, vdwhell,hell
!      INTEGER hellcount
!      COMMON /HELLCO/ hellcount

END MODULE MODAMBER 

