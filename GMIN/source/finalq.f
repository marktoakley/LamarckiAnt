C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
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
      SUBROUTINE FINALQ
      USE COMMONS
      USE GENRIGID
      use QMODULE
      USE MODCHARMM, ONLY: ACESOLV,NCHENCALLS,ACEUPSTEP
      IMPLICIT NONE


      INTEGER J1, J2, ITERATIONS, BRUN,QDONE, J3, NP
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMSALLOC), X(3*NATOMS), ENERGY, GRAD(3*NATOMSALLOC), TIME
      DOUBLE PRECISION SAVECSMNORM, DUMMY2, DIST2, RMAT(3,3), AVVAL, CSMGRAD(3), XTEMP(1:3*NATOMSALLOC)
      DOUBLE PRECISION DUMMY(3*NATOMSALLOC), AA(3)

      COMMON /MYPOT/ POTEL
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
 
C
C  Make sure the lowest minima are tightly converged and then sort
C  them just to be on the safe side.
C
      CSMGUIDET=.FALSE.
C
C Turn off guiding potential.
C
      IF (DFTBCT) LJATT=.FALSE.
      GUIDET=.FALSE.
      SHELLMOVES(1:NPAR)=.FALSE.
      IF (CUTT) CUTOFF=FINALCUTOFF
      SAVEQ=.FALSE.
      NQ(:)=0
      IF (FIELDT) FIELDT=.FALSE.
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      IF (SQUEEZET) SQUEEZET=.FALSE.
      MAXIT=MAXIT2
      IF (MPIT) THEN
         NP=MYNODE+1
      ELSE
         NP=1
      ENDIF
! csw34> Adding in clarification of what is being tightly quenched
      WRITE(MYUNIT,'(A)') 'Tightly converging the SAVE lowest energy minima found'
      WRITE(MYUNIT,'(A)') 'NOTE: these may NOT match the other output files - see below for a sorted list of Lowest minima'
      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            NATOMS=QMINNATOMS(J1)
            DO J2=1,3*NATOMS
               COORDS(J2,NP)=QMINP(J1,J2)
            ENDDO
            
            !ds656>
            !write(MYUNIT,*) (QMINT(J1,J2),J2=1,NATOMS)
            IF(NSPECIES(0)>1) CALL SET_ATOMLISTS(QMINT(J1,1:NATOMS),1)
            !<ds656
            
            NQ(NP)=NQ(NP)+1
            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
!fh301>{{{
            STRUC = NQ(1)
            COUT=.TRUE.
!fh301>}}}
! hk286

            IF (RELAXFQ.OR.HYBRIDMINT) THEN
               RIGIDINIT = .FALSE.
            ENDIF
! hk286
            IF(UNFREEZEFINALQ) FROZEN(:)=.FALSE.  ! unfreeze everything before the final quenches
            CALL QUENCH(.TRUE.,NP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)

!cl457
            IF (CHEMSHIFT2) THEN
               WRITE(MYUNIT,'(A,I6,A,G20.10,A,G20.10,A,I5,A,G15.7,A,F12.2)') 'Final Quench ',NQ(NP),' E=',
     1            POTEL,' ECAM=', ENERGYCAMSHIFT,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART
            ELSE
!cl457
               WRITE(MYUNIT,'(A,I6,A,G20.10,A,I5,A,G15.7,A,F12.2)') 'Final Quench ',NQ(NP),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART
            ENDIF

!cl457
!            IF (CHEMSHIFT2) THEN
!               QMIN(J1)=FULLENERGY
!            ELSE
!cl457
               QMIN(J1)=POTEL
!            ENDIF

            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,NP)
            ENDDO

            IF (CSMT) THEN
               CSMAV(1:3*NATOMS)=0.0D0
               DO J2=1,CSMGPINDEX
!
! rotate permuted image to best orientation with CSMPMAT
! apply point group operation J2
! 
                  DO J3=1,NATOMS
                     XTEMP(3*(J3-1)+1)=CSMPMAT(1,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(1,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(1,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                     XTEMP(3*(J3-1)+2)=CSMPMAT(2,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(2,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(2,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                     XTEMP(3*(J3-1)+3)=CSMPMAT(3,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(3,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(3,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  ENDDO
                  CALL CSMROT(XTEMP,DUMMY,1,J2)
                  CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)+DUMMY(1:3*NATOMS)
               ENDDO
               CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)/CSMGPINDEX
!
!  Check the CSM for the averaged structure. It should be zero if this structure has the
!  right point group. Need to reset CSMIMAGES and CSMNORM temporarily. 
!  
               IF (PERMDIST) THEN
                  DO J2=1,CSMGPINDEX
                     XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
                     CALL CSMROT(XTEMP,DUMMY,1,J2)
                     CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY2,DIST2,RIGID,RMAT)
                     CALL CSMROT(DUMMY,XTEMP,-1,J2) ! need to rotate the permuted rotated images back to the reference orientation
                     CSMIMAGES(1+3*NATOMS*(J2-1):3*NATOMS*J2)=XTEMP(1:3*NATOMS)
                  ENDDO
               ELSE
                  DO J2=1,CSMGPINDEX
                     CSMIMAGES(1+3*NATOMS*(J2-1):3*NATOMS*J2)=CSMAV(1:3*NATOMS)
                  ENDDO
               ENDIF
               AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
               SAVECSMNORM=CSMNORM
               CSMNORM=0.0D0
               DO J2=1,NATOMS
                  CSMNORM=CSMNORM+CSMAV(3*(J2-1)+1)**2+CSMAV(3*(J2-1)+2)**2+CSMAV(3*(J2-1)+3)**2
               ENDDO
               CSMNORM=2*CSMGPINDEX*CSMNORM
               CALL CSMPOTGRAD(CSMAV,AA,AVVAL,.TRUE.,CSMGRAD)
               CSMNORM=SAVECSMNORM
               IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'finalq> CSM for averaged structure=',AVVAL
               QMINAV(J1)=AVVAL
               IF (DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'finalq> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
               QMINPCSMAV(J1,1:3*NATOMS)=CSMAV(1:3*NATOMS)
            ENDIF
         ENDIF
      ENDDO
C
C       sf344> sometimes we can have a lower number of minima found than NSAVE. Resetting
C              NSAVE to the number of minima found should get rid of entries with null 
C              coordinates in the file 'lowest' (and other final output files)
C
C  DJW - this may not work because we may not have found enough minima considered 
C        different according to the EDIFF criterion.
C
      IF (DEBUG) WRITE(MYUNIT,'(A,8I6)') 'finalq> NQ=',NQ(1:NPAR)
      IF (DEBUG) WRITE(MYUNIT,'(A,8I6)') 'finalq> NP=',NP
!      WRITE(MYUNIT,'(A,8I6)') 'finalq> NQ=',NQ(1:NPAR)
!      WRITE(MYUNIT,'(A,8I6)') 'finalq> NP=',NP
      NSAVE=MIN(NSAVE,NQ(NP))
! csw34> Re-sort the saved minima now that they have been tightly converged
      CALL GSORT2(NSAVE,NATOMSALLOC)
! csw34> Print a list of sorted minima energies which will 100% match other output files
      WRITE(MYUNIT,'(A)') 'After re-sorting, the lowest found minima are (lowest free energy subtracted if applicable):'
      DO J1=1,NSAVE
         WRITE(MYUNIT,'(A,I6,A,G20.10)') 'Lowest Minimum ',J1,' Energy= ',QMIN(J1)
      ENDDO
C
C  Optionally sort the atoms from most bound to least bound according to VAT. {{{
C
       
      IF (SORTT) THEN
         DO J1=1,NSAVE
            IF (QMIN(J1).LT.0.0D0) THEN
               DO J2=1,3*NATOMS
                  COORDS(J2,1)=QMINP(J1,J2)
                  X(J2)=QMINP(J1,J2)
               ENDDO
               CALL POTENTIAL(X,GRAD,ENERGY,.FALSE.,.FALSE.)
               DO J2=1,NATOMS
                  VAT(J2,1)=VT(J2)
               ENDDO
               CALL SORT3(NATOMS,NATOMS,VAT(1:NATOMS,1),COORDS(1:3*NATOMS,1))
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=COORDS(J2,1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

C     IF (DEBUG) THEN
         IF (TABOOT) THEN
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A)') 'Taboo lists:'
               DO J1=1,NPAR
                  WRITE(MYUNIT,'(A,G20.10)') 'Parallel run ',J1
                  WRITE(MYUNIT,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
               ENDDO
            ELSE
               WRITE(MYUNIT,'(A)') 'Taboo list:'
               WRITE(MYUNIT,'(6F15.7)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
C     ENDIF
      RETURN
      END
