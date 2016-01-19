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
MODULE consts_trans_97 
      IMPLICIT NONE
      SAVE

!   These are the constants taken from PRB 15th April (1997)          
!   the IMPROVED transferable potential of Menon and Subbaswamy.
!   I have also used data lines to fix the values of VECT and EPS
!   Don't know how this affects speed.

      DOUBLE PRECISION CHI,BETA,D,A,B,RC,DELTA,KAPPA,SIGMA,&
     &                 ALPHA,VECT(4),EPS(4),EPSs,EPSp,FUDGE
!     1                 Vsssig,Vspsig,Vppsig,Vpppi

      PARAMETER (CHI=0.41D0, D=2.36D0, A=0.08D0,&
     &    B=-1.4D0, RC=3.5D0, DELTA=0.1D0, SIGMA=2.5D0,&
     &    ALPHA = 1.62D0, KAPPA=1.7D0, BETA=4.0D0*ALPHA,&
     &    EPSs=-13.55D0, EPSp=-6.52D0, FUDGE=1.1D0)
      DATA VECT(1), VECT(2), VECT(3), VECT(4)/ -2.37D0,2.52D0,3.32D0,-1.07D0 /
      DATA EPS(1), EPS(2), EPS(3), EPS(4)/ -27.1D0,-20.07D0,-13.04D0,-13.04D0 /

!   First those for the repulsive energy term
!   CHI is CHI(0), BETA and D are pretty self explanatory
!   although BETA was =4/R(0), where R(0) is half the dimer bond length
!   and D is the bond length of a crystal.
!   However, in the transferable version of Menon and subbaswamy's potential
!   BETA=4*ALPHA

!   SIGMA is the new parameter introduced for the improved 
!   transferable parameter.

!   Then in the bond counting terms
!   A and B are fitted parameters so tha the cohesive energies of
!   several clusters in close agreement with the corresponding ab 
!   initio values.
!   RC is NOT a cutoff distance for interaction apparently but
!   merely a distance for calibrating the energy term.
!   There is no cutoff between atoms in the Menon-Subbaswamy scheme.

!   VECT(1) is the Vsssig parameter allegedly from Harrison's 
!   universal parameter scheme (These figures come from 
!   PRB 50 12754 (1993) and I have assumed that they are at the
!   crystal bond length (D=2.35) although this does not appear to be
!   consistent with the scheme in Harrrison't book "Electronic 
!   Structure and The Properties if Solids".
!   Similarly, VECT(2) is Vspsig, VECT(3) Vppsig and VECT(4) Vpppi.
!   These are assigned in the main program silicon.f
!   KAPPA, ZETA and ALPHA are all pretty self-explanatory.

!   The PARAMETER function will not accept ARRAYS...most useful.

END MODULE consts_trans_97
