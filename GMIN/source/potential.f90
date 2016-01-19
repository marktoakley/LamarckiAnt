
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
SUBROUTINE POTENTIAL(X, GRAD, EREAL, GRADT, SECT)
! This subroutine acts as a wrapper for all of the potentials implemented in GMIN.
! The potential is chosen in the 'data' file (default is Lennard-Jones particles).
! 
! In addition, this routine also applies fields, if appropriate. 
   USE COMMONS
   USE GENRIGID, ONLY: RIGIDINIT, NRIGIDBODY, DEGFREEDOMS, ATOMRIGIDCOORDT, &
                       AACONVERGENCET, AACONVERGENCE, TRANSFORMRIGIDTOC, &
                       TRANSFORMGRAD, TRANSFORMHESSIAN, TRANSFORMCTORIGID, &
                       GENRIGID_DISTANCECHECK, GENRIGID_COMPRESS
   USE QMODULE, ONLY: QMIN, QMINP 
   USE PERMU, ONLY: FIN
   USE PORFUNCS
   USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_ENERGY_AND_GRADIENT, POT_ENE_REC_C
   USE MODHESS, ONLY: HESS, RBAANORMALMODET 
   USE MODAMBER9, ONLY: ATMASS1
   USE POLIRMOD, ONLY: POLIR
   USE MBPOLMOD, ONLY: MBPOL
   USE SWMOD, ONLY: SWTYPE
   USE MODCUDALBFGS, ONLY: CUDA_ENEGRAD_WRAPPER, CUDA_NUMERICAL_HESS
   USE GAUSS_MOD, ONLY: GFIELD
   USE RAD_MOD, ONLY: RAD, RADR
   USE PREC, ONLY: INT32, REAL64
   USE TWIST_MOD, ONLY: TWISTT, TWIST
   IMPLICIT NONE
! Arguments
! TODO: Work out intents
! TODO: Fix array dimensions?
   REAL(REAL64)               :: X(*)
   REAL(REAL64)               :: GRAD(*)
   REAL(REAL64)               :: EREAL
   LOGICAL, INTENT(IN)        :: GRADT
   LOGICAL, INTENT(IN)        :: SECT
! Variables      
! hk286 > generalised rigid body additions
   REAL(REAL64)               :: XCOORDS(3*NATOMS), GRADATOMS(3*NATOMS), &
                                 XRIGIDCOORDS(DEGFREEDOMS), XRIGIDGRAD(DEGFREEDOMS)
   REAL(REAL64), ALLOCATABLE  :: TEMPHESS(:)
   REAL(REAL64)               :: GRAD1(3*NATOMS)
   INTEGER(INT32)             :: I, J, K
   REAL(REAL64)               :: XRIGIDHESS(DEGFREEDOMS, DEGFREEDOMS)
   LOGICAL                    :: FTEST, EVAP, COMPON, YESNO, GUIDET, EVAPREJECT
   INTEGER(INT32)             :: J1, J2, J3, CSMIT
   CHARACTER                  :: FNAME*80, DUMM*4
   REAL(REAL64)               :: DUMMY2, GEMAX, RMAT(3, 3), XTEMP(3*NATOMS), AVVAL, &
                                 DUMMY(3*NATOMS), DIST2, QE, QX, AA(3), SAVECSMNORM, &
                                 CSMRMS, CSMGRAD(3), SAVECSMPMAT(3, 3), &
                                 SAVECSMIMAGES(3*NATOMS*CSMGPINDEX)
   CHARACTER(LEN=3)           :: CSMGPSAVE
   INTEGER(INT32)             :: CSMGPINDEXSAVE  
   REAL(REAL64)               :: PTGPSAVE(3, 3, 2*CSMGPINDEX), CSMNORMSAVE, ENERGY, VNEW(3*NATOMS)
   INTEGER(INT32)             :: BRUN
   INTEGER(INT32)             :: NQTOT
   LOGICAL                    :: SOCOUPLE, GUIDECHANGET, CSMDOGUIDET, COMTEST
! khs26 > AMBER12 energy decomposition, type defined in AMBER12 interface
   TYPE(POT_ENE_REC_C)        :: AMBER12_ENERGY_DECOMP
! Used only in commented sections
!   REAL(REAL64)               :: XG, YG, ZG, GRADLJ(3*NATOMS), EREALLJ, GRADMF(3*NATOMS), EREALMF, TERMLJ, TERMMF &
!                                 GRADDUM(3*NATOMS), GRADDUM1(3*NATOMS), GRADDUM2(3*NATOMS), EPLUS, EMINUS, DIFF
!                                 EDUM1, EDUM2
!DS656> GRADIENT/HESSIAN TESTING
   DOUBLE PRECISION :: GPLUS(3*NATOMS), GMINUS(3*NATOMS), DIFF, &
        EPLUS, EMINUS
   
! TODO: Delete common blocks
   COMMON /CO/ COMPON
   COMMON /FAIL/ FTEST
   COMMON /EV/ EVAP, EVAPREJECT
   COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
   COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT
   COMMON /TOT/ NQTOT

   GUIDECHANGET=.FALSE.
!
! Test BRUN to see if we should stop if a screen saver is interrupted.
! Need to save a restart file containing:
! Current minimum in the Markov chain. COORDS
! Number of steps done. NQTOT/NPAR should be close enough!
! The current lowest minima. QMIN has the energies, QMINP has the points.
! The current values of the temperature, acceptance ratio and step length,
! TEMP(JP), ACCRAT(JP), STEP(JP), ASTEP(JP)
! which can get changed dynamically.
!
! khs26 > This never seems to be executed, because of BRUN=0
!   BRUN=0
!   IF (BRUN == 1) THEN
!      WRITE(MYUNIT,'(A)' ) 'dumping restart file ssdump'
!      OPEN(UNIT=88, FILE='ssdump', STATUS='UNKNOWN')
!      WRITE(88,'(3G20.10)') ((COORDS(J1, J2), J1=1, 3*NATOMS), J2=1, NPAR)
!      WRITE(88,'(I6)') NQTOT/NPAR, NPCALL
!      WRITE(88,'(G20.10)') (QMIN(J1), J1=1, NSAVE)
!      WRITE(88,'(3G20.10)') ((QMINP(J2, J1), J1=1, 3*NATOMS), J2=1, NSAVE)
!      WRITE(88,'(G20.10)') (TEMP(J1), J1=1, NPAR)
!      WRITE(88,'(G20.10)') (ACCRAT(J1), J1=1, NPAR)
!      WRITE(88,'(G20.10)') (STEP(J1), J1=1, NPAR)
!      WRITE(88,'(G20.10)') (ASTEP(J1), J1=1, NPAR)
!      WRITE(88,'(G20.10)') (OSTEP(J1), J1=1, NPAR)
!      CALL SYSTEM('rm ssave')
!      STOP
!   END IF

! Count the number of potential calls.
   NPCALL=NPCALL+1

10 CONTINUE

   IF (MSORIGT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      IF (CUTT) THEN
         CALL MSORIGC(NATOMS, X, GRAD, EREAL, GRADT)
      ELSE
         CALL MSORIG(NATOMS, X, GRAD, EREAL, GRADT)
      END IF
      IF (FTEST) THEN
         RETURN
      END IF
   ELSE IF (MSTRANST) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL MSTRANS97(NATOMS, X, GRAD, EREAL, GRADT)
      IF (FTEST) THEN
         RETURN
      END IF
   ELSE IF (FRAUSIT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL FRAUSI(NATOMS, X, GRAD, EREAL, GRADT, ANGST, NATOMS)
      IF (FTEST) THEN
         RETURN
      END IF
!
!  DIM Ne^+, Ne*, Ar^+, Ar*
!
   ELSE IF ((NEON .OR. ARGON) .AND. (PLUS .OR. STAR)) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
!      CALL RGNI(NATOMS, X, GRAD, EREAL, GRADT, h0, h1, ee, ev, w, NATOMS)
      CALL RGNI(NATOMS, X, GRAD, EREAL, GRADT)
!
!  DIM Ar^{2+}
!
   ELSE IF (TWOPLUS) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
!      CALL RGNII(NATOMS, X, GRAD, EREAL, GRADT, h0, h1, ee, ev, w, NATOMS)
      CALL RGNII(NATOMS, X, GRAD, EREAL, GRADT)
!
!  DIM Ar and Ne neutrals.
!
   ELSE IF (GROUND) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL GRND(NATOMS, X, GRAD, EREAL, GRADT)
      IF (AXTELL) CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)
   ELSE IF (RGCL2) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL RGNX2(NATOMS, X, GRAD, EREAL, GRADT)
      IF (AXTELL) CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)
   ELSE IF (ARNO) THEN
      SOCOUPLE=.TRUE.
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL RGNXY(NATOMS, X, GRAD, EREAL, GRADT, SOCOUPLE)
!      WRITE(*, *) 'SOCOUPLE=', SOCOUPLE
!      WRITE(*, *) 'POINTS:'
!      WRITE(*,'(3F20.10)') (X(J1), J1=1, 3*NATOMS)
!      WRITE(*, *) 'Energy=', EREAL
!      WRITE(*, *) 'Gradient:'
!      WRITE(*,'(3F20.10)') (GRAD(J1), J1=1, 3*NATOMS)
!      IF (GRADT) THEN
!         DISP=1.0D-4
!         DO J1=1, 3*NATOMS
!            DUMMY2=X(J1)
!            X(J1)=X(J1)+DISP
!            CALL RGNXY(NATOMS, X, GRADNUM, V1,.FALSE., SOCOUPLE)
!            X(J1)=X(J1)-2.0D0*DISP
!            CALL RGNXY(NATOMS, X, GRADNUM, V2,.FALSE., SOCOUPLE)
!            GRADNUM(J1)=(V1-V2)/(2.0D0*DISP)
!            X(J1)=DUMMY2
!         END DO
!      END IF
!      WRITE(*, *) 'Numerical derivatives for displacement ', DISP
!      WRITE(*,'(3F20.10)') (GRADNUM(J1), J1=1, 3*NATOMS)
!      WRITE(*, *) 'Analytic/numerical derivatives'
!      WRITE(*,'(3F20.10)') (GRAD(J1)/GRADNUM(J1), J1=1, 3*NATOMS)

      IF (AXTELL) CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)
   ELSE IF (OTPT) THEN
      WRITE(MYUNIT,'(A)') 'potential> OTP subroutine is commented'
      STOP
!      CALL OTP(X, GRAD, EREAL, GRADT, SECT)
!      DIFF=1.0D-4
!      WRITE(*, *) 'analytic and numerical gradients:'
!      DO J1=1, 3*NATOMS
!         X(J1)=X(J1)+DIFF
!         CALL OTP(X, GRADDUM, EPLUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)-2.0D0*DIFF
!         CALL OTP(X, GRADDUM, EMINUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)+DIFF
!         IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.1.0D0)) THEN
!            WRITE(*,'(I5, 2F20.10)') J1, GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!         END IF
!      END DO

   ELSE IF (LJMFT) THEN
      WRITE(MYUNIT,'(A)') 'potential> LJ/MF subroutines are commented'
      STOP
!      CALL LJ(X, GRADLJ, EREALLJ, GRADT, SECT)
!      CALL MF(X, GRADMF, EREALMF, GRADT)
!      CALL LJ(X, GRADLJ, EREALLJ,.FALSE.,.FALSE.)
!      CALL MF(X, GRADMF, EREALMF,.FALSE.)
!      WRITE(*,'(A, G20.10)') 'radius=', RADIUS
!      WRITE(*, *) 'EREALLJ, EREALMF=', EREALLJ, EREALMF
!      EREAL=EREALLJ+EREALMF
!      TERMLJ=EREALLJ
!      TERMMF=EREALMF
!      DO J1=1, 3*NATOMS
!         WRITE(*, *) 'J1, GRADLJ, GRADMF=', J1, GRADLJ(J1), GRADMF(J1)
!         GRAD(J1)=GRADLJ(J1)+GRADMF(J1)
!      END DO
!      DIFF=1.0D-6
!      WRITE(*, *) 'analytic and numerical gradients:'
! 
!      DO J1=1, 3*NATOMS
!         X(J1)=X(J1)+DIFF
!         CALL LJ(X, GRADDUM1, EDUM1,.FALSE.,.FALSE.)
!         CALL MF(X, GRADDUM2, EDUM2,.FALSE.)
!         EPLUS=EDUM1+EDUM2
! 
!         X(J1)=X(J1)-2.0D0*DIFF
! 
!         CALL LJ(X, GRADDUM1, EDUM1,.FALSE.,.FALSE.)
!         CALL MF(X, GRADDUM2, EDUM2,.FALSE.)
! 
!         EMINUS=EDUM1+EDUM2
! 
!         X(J1)=X(J1)+DIFF
!         GRAD(J1)=(EPLUS-EMINUS)/(2.0D0*DIFF)
!         WRITE(*,'(I5, 4F20.10)') J1, GRAD(J1), EPLUS, EMINUS, ABS(EPLUS-EMINUS)
!     END DO
!

   ELSE IF (MORSET) THEN
!      CALL RAD(X, GRAD, EREAL, GRADT)
!      CALL MORSE(X, GRAD, EREAL, GRADT)
      CALL MORSEGH(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (QCIPOTT) THEN
!     IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREPQCIPOT(X,3*NATOMS,0,1)
      IF (QCIPOT2T) THEN
         CALL QCIPOT2(EREAL,X,GRAD)
      ELSE
         CALL QCIPOT(EREAL,X,GRAD)
      ENDIF

   ELSE IF (TOSI) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      IF (EVAPREJECT) RETURN
      CALL TOSIFUMI(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (WELCH) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL WEL(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (SCT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL SC(X, GRAD, EREAL, GRADT)
      IF (CPMD) EREAL=EREAL+1.0D6

   ELSE IF (MSCT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL MSC(X, GRAD, EREAL, GRADT)

   ELSE IF (ACKLANDT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL ACK(X, GRAD, EREAL, GRADT)

   ELSE IF (FAL.OR.FNI) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL FARKAS(X, GRAD, EREAL, GRADT, NATOMS)

   ELSE IF (LJATT) THEN
      CALL LJ(X, GRAD, EREAL, GRADT, SECT)
      CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)

   ELSE IF (DFTBCT) THEN
      IF (.NOT.PERCOLATET) CALL RAD(X, GRAD, EREAL, GRADT)
      CALL DFTBC(NATOMS, X, GRAD, EREAL, GRADT)
      IF (FTEST) THEN
         RETURN
      END IF

   ELSE IF (DFTBT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
!      CALL SECDIFF(NATOMS, X, W)
      CALL DFTB(NATOMS, X, GRAD, EREAL, GRADT)
      IF (FTEST) THEN
         RETURN
      END IF

   ELSE IF (P46) THEN
!      CALL P46MER(X, GRAD, EREAL, GRADT)
      CALL P46MERDIFF(X, NATOMS, GRAD, EREAL, GRADT)

   ELSE IF (G46) THEN
      CALL G46MERDIFF(X, NATOMS, GRAD, EREAL, GRADT)

   ELSE IF (BLNT) THEN
      CALL BLN(X, GRAD, EREAL, GRADT)

   ELSE IF (CHAPERONINT) THEN
      CALL CHAPERONIN(X, GRAD, EREAL, GRADT)

   ELSE IF (SW) THEN
      IF (PERIODIC) THEN
         CALL SISW(X, GRAD, EREAL, GRADT)
      ELSE
         CALL RAD(X, GRAD, EREAL, GRADT)
         CALL SWTWO(X, GRAD, EREAL, GRADT)
         CALL SWTHREE(GRAD, EREAL, GRADT)
      END IF

   ELSE IF (QUIPT) THEN
      IF (.NOT. PERCOLATET) CALL RAD(X, GRAD, EREAL, GRADT)
      CALL GMIN_QUIP_WRAPPER(X, GRAD, EREAL, QUIPATOMTYPE, QUIPARGSTR)

   ELSE IF (AMHT) THEN
      CALL WALESAMH_INTERFACE(X, GRAD, EREAL)
!      DIFF=1.0D-4
!      SUMMDIFF=0.D0
!      WRITE(*, *) 'analytic and numerical gradients:'
!      DO J1=1, 3*NATOMS
!         WRITE(*,'(F20.10, 2x, F20.10, 2xI5)')X(J1), GRAD(J1), J1
!         X(J1)=X(J1)+DIFF
!         CALL SCLTAB_WALES(X, GRADDUM, EPLUS)
!         X(J1)=X(J1)-2.0D0*DIFF
!         CALL SCLTAB_WALES(X, GRADDUM, EMINUS)
!         X(J1)=X(J1)+DIFF
!         IF (100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1)) > 0.0D0) THEN
!            WRITE(*,'(I5, 3F15.8)') J1, GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), &
!                                   100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!         END IF
!         SUMMDIFF=SUMMDIFF + (GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1)
!         GRAD(J1)=(EPLUS-EMINUS)/(2.0D0*DIFF)
!      END DO
!      WRITE(6,*)'SUMM DIFF ', SUMMDIFF

! khs26> AMBER12 energy and gradient call
   ELSE IF (AMBER12T) THEN
! If coordinates are atomistic (i.e. no rigid bodies) just call the energy and gradient
      IF (ATOMRIGIDCOORDT) THEN
         IF (CUDAT) THEN
! This call copies CPU coordinates to GPU, calculates energy/gradient and copies energy/gradient back to CPU
! Calls to CUDA_ENEGRAD_WRAPPER are for debugging/testing/printing purposes
            CALL CUDA_ENEGRAD_WRAPPER(NATOMS, X, EREAL, GRADATOMS)
            GRAD(1:3*NATOMS)=GRADATOMS(:)
         ELSE
            CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, X, EREAL, GRADATOMS, AMBER12_ENERGY_DECOMP)
            GRAD(1:3*NATOMS)=GRADATOMS(:)
         END IF
! If the coordinates include rigid bodies, transform them back to being atomistic first
      ELSE
         XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
! Now call the AMBER 12 energy and gradient
         IF (CUDAT) THEN
            CALL CUDA_ENEGRAD_WRAPPER(NATOMS, XCOORDS, EREAL, GRADATOMS)
            GRAD(1:3*NATOMS)=GRADATOMS(:)
         ELSE
            CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, XCOORDS, EREAL, GRADATOMS, AMBER12_ENERGY_DECOMP)
            GRAD(1:3*NATOMS)=GRADATOMS(:)
         END IF
! Transform the gradient and coordinates back to the rigid body representation 
         CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
         X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
         X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
         GRAD(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
         GRAD(DEGFREEDOMS+1:3*NATOMS)=0.0D0
      END IF
! AMBER 9 Energy and gradient calls
   ELSE IF (AMBERT) THEN
! hk286 > Generalised rigid body
! hk286 > If rigid coords is used, then need to convert to atom coords first, 
! hk286 > compute energies, then convert gradients & coords back to rigid
      IF (ATOMRIGIDCOORDT) THEN
         CALL AMBERENERGIES(X, GRAD, EREAL, .FALSE., .FALSE.)
         IF (RESTRAINLT) THEN
            CALL RESTRAINLPOTENTIAL(X, GRAD, EREAL, GRADT, SECT) ! hk286
         END IF
      ELSE
         XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
         CALL AMBERENERGIES(XCOORDS, GRADATOMS, EREAL, .FALSE., .FALSE.)
         IF (RESTRAINLT) THEN
            CALL RESTRAINLPOTENTIAL(XCOORDS, GRADATOMS, EREAL, GRADT, SECT) ! hk286
         END IF
         CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
         X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
         X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
         GRAD(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
         GRAD(DEGFREEDOMS+1:3*NATOMS)=0.0D0
! csw34> apply compression if using COMPRESSRIGID and at least one COM-COM distance
!        is greater than RIGIDCOMDIST
         IF (COMPRESSRIGIDT) THEN
            COMTEST=.TRUE.
            CALL GENRIGID_DISTANCECHECK(XRIGIDCOORDS, RIGIDCOMDIST, COMTEST)
! If the threshold distance is exceeded, apply the compression
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL before compression=', EREAL
            IF (.NOT.COMTEST) CALL GENRIGID_COMPRESS(XRIGIDCOORDS, GRAD, EREAL, KCOMP_RIGID) 
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL after compression=', EREAL
! Equate GRAD and XRIGIDGRAD to allow new RMS convergence condition to work
            XRIGIDGRAD(1:DEGFREEDOMS)=GRAD(1:DEGFREEDOMS) 
         END IF  
      END IF

      IF (SECT) THEN
         VNEW=0.0D0
! khs26> Copied analytical second derivatives from OPTIM
         IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
            XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, X, XRIGIDCOORDS)               
         END IF
         IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
         ALLOCATE(HESS(3*NATOMS, 3*NATOMS))
         IF (CUDAT) THEN
            CALL CUDA_NUMERICAL_HESS(NATOMS, X, HESS, DELTA=1.0D-6)
         ELSE 
            IF (.NOT. ALLOCATED(TEMPHESS)) ALLOCATE(TEMPHESS(9*NATOMS*NATOMS))
            TEMPHESS(:)=0.0D0
            CALL MME2WRAPPER(X, ENERGY, VNEW, TEMPHESS, ATMASS1, GRAD1)
            IF (SQRT(ABS(ENERGY)) > 1.0D40) THEN
               EREAL = ENERGY
               WRITE(MYUNIT, *) "Linear dihedral detected in NAB routines (sff2.c)."
            END IF
            K=1
            DO I=1, 3*NATOMS
               DO J=1, 3*NATOMS
                  HESS(I, J)=TEMPHESS(K)
                  K=K + 1
               END DO
            END DO
            DEALLOCATE(TEMPHESS)
         END IF
         IF ( RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
            CALL TRANSFORMHESSIAN(HESS, VNEW, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
            CALL TRANSFORMGRAD(VNEW, XRIGIDCOORDS, XRIGIDGRAD)
            X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
            X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
            VNEW(DEGFREEDOMS+1:3*NATOMS)=0.0D0
            VNEW(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
            HESS(DEGFREEDOMS+1:3*NATOMS,:)=0.0D0
            HESS(:, DEGFREEDOMS+1:3*NATOMS)=0.0D0
            HESS(1:DEGFREEDOMS, 1:DEGFREEDOMS)=XRIGIDHESS(1:DEGFREEDOMS, 1:DEGFREEDOMS)
         END IF
      END IF
! End of second derivative calculation 

! End of AMBERT section

! hk286
   ELSE IF (CHRMMT) THEN
! hk286 > Generalised rigid body
! hk286 > If rigid coords is used, then need to convert to atom coords first, 
! hk286 > compute energies, then convert gradients & coords back to rigid
      IF (ATOMRIGIDCOORDT) THEN
         CALL OCHARMM(X, GRAD, EREAL, GRADT)
      ELSE
         XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
         CALL OCHARMM(XCOORDS, GRADATOMS, EREAL, GRADT)
         CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
         X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
         X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
         GRAD(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
         GRAD(DEGFREEDOMS+1:3*NATOMS)=0.0D0
! csw34> apply compression if using COMPRESSRIGID and at least one COM-COM distance
!        is greater than RIGIDCOMDIST
         IF (COMPRESSRIGIDT) THEN
            COMTEST=.TRUE.
            CALL GENRIGID_DISTANCECHECK(XRIGIDCOORDS, RIGIDCOMDIST, COMTEST)
! If the threshold distance is exceeded, apply the compression
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL before compression=', EREAL
            IF (.NOT.COMTEST) CALL GENRIGID_COMPRESS(XRIGIDCOORDS, GRAD, EREAL, KCOMP_RIGID) 
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL after compression=', EREAL
! Equate GRAD and XRIGIDGRAD to allow new RMS convergence condition to work
            XRIGIDGRAD(1:DEGFREEDOMS)=GRAD(1:DEGFREEDOMS) 
         END IF  
      END IF

! hk286
   ELSE IF (DMACRYST) THEN
      IF (.NOT. ATOMRIGIDCOORDT) THEN
         CALL DMACRYS_POTENTIAL(X, GRAD, EREAL, GRADT)
      ELSE
         CALL TRANSFORMCTORIGID(X, XRIGIDCOORDS)
         CALL DMACRYS_POTENTIAL(XRIGIDCOORDS, GRAD, EREAL, GRADT)
      END IF

   ELSE IF (DZTEST) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL DZPOT(X, GRAD, EREAL, GRADT, SECT)
! vr274: userpot, linkage dependant potential

   ELSE IF (USERPOTT) THEN
! hk286 > Generalised rigid body
! hk286 > If rigid coords is used, then need to convert to atom coords first, 
! hk286 > compute energies, then convert gradients & coords back to rigid
      IF (ATOMRIGIDCOORDT) THEN
         CALL USERPOT_POTENTIAL(3*NATOMS, X, GRAD, EREAL, GRADT)
         IF (RESTRAINLT) THEN
            CALL RESTRAINLPOTENTIAL(X, GRAD, EREAL, GRADT, SECT) ! hk286
         END IF
      ELSE
         XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
         CALL USERPOT_POTENTIAL(3*NATOMS, XCOORDS, GRADATOMS, EREAL, GRADT)
         IF (RESTRAINLT) THEN
            CALL RESTRAINLPOTENTIAL(XCOORDS, GRADATOMS, EREAL, GRADT, SECT) ! hk286
         END IF
         CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
         X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
         X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
         GRAD(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
         GRAD(DEGFREEDOMS+1:3*NATOMS)=0.0D0
! csw34> apply compression if using COMPRESSRIGID and at least one COM-COM distance
!        is greater than RIGIDCOMDIST
         IF (COMPRESSRIGIDT) THEN
            COMTEST=.TRUE.
            CALL GENRIGID_DISTANCECHECK(XRIGIDCOORDS, RIGIDCOMDIST, COMTEST)
! If the threshold distance is exceeded, apply the compression
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL before compression=', EREAL
            IF (.NOT.COMTEST) CALL GENRIGID_COMPRESS(XRIGIDCOORDS, GRAD, EREAL, KCOMP_RIGID) 
            IF (DEBUG) WRITE(MYUNIT,'(A, F20.10)' ) 'EREAL after compression=', EREAL
! Equate GRAD and XRIGIDGRAD to allow new RMS convergence condition to work
            XRIGIDGRAD(1:DEGFREEDOMS)=GRAD(1:DEGFREEDOMS) 
         END IF  
      END IF

      IF (SECT) THEN
         VNEW=0.0D0
! khs26> Copied analytical second derivatives from OPTIM
         IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
            XRIGIDCOORDS(1:DEGFREEDOMS)=X(1:DEGFREEDOMS)
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, X, XRIGIDCOORDS)               
         END IF
         IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
         IF (.NOT. ALLOCATED(TEMPHESS)) ALLOCATE(TEMPHESS(9*NATOMS*NATOMS))
         TEMPHESS(:)=0.0D0
         CALL MME2WRAPPER(X, ENERGY, VNEW, TEMPHESS, ATMASS1, GRAD1)
         IF (SQRT(ABS(ENERGY)) > 1.0D40) THEN
            EREAL = ENERGY
            WRITE(MYUNIT, *) "Linear dihedral detected in NAB routines (sff2.c)."
         END IF
         ALLOCATE(HESS(3*NATOMS, 3*NATOMS))
         K=1
         DO I=1, 3*NATOMS
            DO J=1, 3*NATOMS
               HESS(I, J)=TEMPHESS(K)
               K=K + 1
            END DO
         END DO
         DEALLOCATE(TEMPHESS)
         IF ( RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
            CALL TRANSFORMHESSIAN(HESS, VNEW, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
            CALL TRANSFORMGRAD(VNEW, XRIGIDCOORDS, XRIGIDGRAD)
            X(DEGFREEDOMS+1:3*NATOMS)=0.0D0
            X(1:DEGFREEDOMS)=XRIGIDCOORDS(1:DEGFREEDOMS)
            VNEW(DEGFREEDOMS+1:3*NATOMS)=0.0D0
            VNEW(1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
            HESS(DEGFREEDOMS+1:3*NATOMS,:)=0.0D0
            HESS(:, DEGFREEDOMS+1:3*NATOMS)=0.0D0
            HESS(1:DEGFREEDOMS, 1:DEGFREEDOMS)=XRIGIDHESS(1:DEGFREEDOMS, 1:DEGFREEDOMS)
         END IF
      END IF

   
   ELSE IF (ZETT1) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL Z1(X, GRAD, EREAL, GRADT)

   ELSE IF (ZETT2) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL Z2(X, GRAD, EREAL, GRADT)

   ELSE IF (PACHECO) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL PRC60(X, GRAD, EREAL, GRADT)
      IF (AXTELL) CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)

   ELSE IF (MODEL1T) THEN
      CALL MODEL1(X, GRAD, EREAL, QE, QX)
   
   ELSE IF (EAMLJT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL EAMLJ(X, GRAD, EREAL, GRADT)

   ELSE IF (PBGLUET) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL PBGLUE(X, GRAD, EREAL, GRADT)
   
   ELSE IF (FST) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL FS(X, GRAD, EREAL, GRADT)

   ELSE IF (BGUPTAT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL BGUPTA(X, GRAD, EREAL, GRADT)

   ELSE IF (MGUPTAT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL MGUPTA(X, GRAD, EREAL, GRADT, SECT, .FALSE.)
      ! Testing...
      !ds656> Test gradient and Hessian
      IF (SECT.AND..FALSE.) THEN
         DIFF=1.0D-5
         PRINT*,SECT
         PRINT*,'analytic and numerical gradients:'
         DO J1=1,3*NATOMS
            X(J1)=X(J1)+DIFF
            CALL MGUPTA(X,GPLUS,EPLUS,.FALSE.,.FALSE.,.false.)
            X(J1)=X(J1)-2.0D0*DIFF
            CALL MGUPTA(X,GMINUS,EMINUS,.FALSE.,.FALSE.,.false.)
            X(J1)=X(J1)+DIFF
            IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.1.0D0)) THEN
               WRITE(*,'(I5,2F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ENDIF
         ENDDO
         PRINT*,'analytic and numerical second derivatives:'
         DO J1=1,3*NATOMS
            X(J1)=X(J1)+DIFF
            CALL MGUPTA(X,GPLUS,EPLUS,.TRUE.,.FALSE.,.false.)
            X(J1)=X(J1)-2.0D0*DIFF
            CALL MGUPTA(X,GMINUS,EMINUS,.TRUE.,.FALSE.,.false.)
            X(J1)=X(J1)+DIFF
            DO J2=1,3*NATOMS
               IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND. &
                    (ABS(100.0D0*(HESS(J1,J2)-(GPLUS(J2)-GMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
                  WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(GPLUS(J2)-GMINUS(J2))/(2.0D0*DIFF),'   X'
               ELSE
                  WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(GPLUS(J2)-GMINUS(J2))/(2.0D0*DIFF)
               ENDIF
            ENDDO
         ENDDO
         STOP
      ENDIF
      ! <ds656 ...end of testing.
   ELSE IF (GLJY) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL GLJYPOT(X, GRAD, EREAL, GRADT)

   ELSE IF (GUPTAT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL GUPTA(X, GRAD, EREAL, GRADT)

   ELSE IF (NATBT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL NATB(NATOMS, X, GRAD, EREAL, GRADT, GUIDET)
!      DIFF=1.0D-3
!      WRITE(*, *) 'analytic and numerical gradients:'
!      DO J1=1, 3*NATOMS
!         X(J1)=X(J1)+DIFF
!         CALL NATB(NATOMS, X, GRADDUM, EPLUS,.FALSE.)
!         X(J1)=X(J1)-2.0D0*DIFF
!         CALL NATB(NATOMS, X, GRADDUM, EMINUS,.FALSE.)
!         X(J1)=X(J1)+DIFF
!         IF (GRAD(J1).NE.0.0D0) WRITE(*,'(I5, 3F20.10)') J1, GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), &
!                                100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!      END DO

! dm368 for the phi^4 model
   ELSE IF (PHI4MODELT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL PHI4MODEL(X, GRAD, EREAL, GRADT)

   ELSE IF (ALGLUET) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL ALGLUE(X, GRAD, EREAL, GRADT)

   ELSE IF (MGGLUET) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL MGGLUE(X, GRAD, EREAL, GRADT)
   
   ELSE IF (EAMALT) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL EAMAL(X, GRAD, EREAL, GRADT)

   ELSE IF (WENZEL) THEN
      CALL WEN(X, GRAD, EREAL, GRADT, NATOMS)
   
   ELSE IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
!      IF (DEBUG) OPEN(UNIT=765, FILE='CSMrot.xyz', STATUS='UNKNOWN')
      IF (CSMDOGUIDET) THEN ! we want to use a guiding group
         PTGPSAVE(1:3, 1:3, 1:2*CSMGPINDEX)=PTGP(1:3, 1:3, 1:2*CSMGPINDEX) ! before we change the index!
         CSMNORMSAVE=CSMNORM
         CSMNORM=CSMGUIDENORM
         CSMGPSAVE=CSMGP
         CSMGP=CSMGUIDEGP
         CSMGPINDEXSAVE=CSMGPINDEX
         CSMGPINDEX=CSMGUIDEGPINDEX
         PTGP(1:3, 1:3, 1:2*CSMGPINDEX)=PTGPGUIDE(1:3, 1:3, 1:2*CSMGPINDEX)
      END IF
      IF (PERMDIST) THEN
         DO J1=1, CSMGPINDEX
            XTEMP(1:3*NATOMS)=X(1:3*NATOMS)
            CALL CSMROT(XTEMP, DUMMY, 1, J1)
            CALL MINPERMDIST(XTEMP, DUMMY, NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, PERIODIC, TWOD, EREAL, DIST2, RIGID, RMAT)
            CALL CSMROT(DUMMY, XTEMP,-1, J1) ! need to rotate the permuted rotated images back to the reference orientation
            CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
         END DO
      ELSE
         DO J1=1, CSMGPINDEX
            CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=X(1:3*NATOMS)
         END DO
      END IF
      CALL CSMMIN(X, EREAL, CSMRMS, CSMIT)

      IF (DEBUG) THEN
!
!  Saving CSMIMAGES, CSMPMAT is necessary here because otherwise these quantities will
!  be different in finalq when we need to write out the results to file CSMav.xyz.
!
         CSMAV(1:3*NATOMS)=0.0D0
!         WRITE(MYUNIT,'(A)') 'potential> CSMPMAT:'
!         WRITE(MYUNIT,'(3G20.10)') CSMPMAT(1:3, 1:3)
         SAVECSMIMAGES(1:3*NATOMS*CSMGPINDEX)=CSMIMAGES(1:3*NATOMS*CSMGPINDEX)

         DO J2=1, CSMGPINDEX
!
! Rotate permuted image to best orientation with CSMPMAT
! Then apply point group operation J2
!
            DO J3=1, NATOMS
               XTEMP(3*(J3-1)+1)=CSMPMAT(1, 1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
                                +CSMPMAT(1, 2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
                                +CSMPMAT(1, 3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
               XTEMP(3*(J3-1)+2)=CSMPMAT(2, 1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
                                +CSMPMAT(2, 2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
                                +CSMPMAT(2, 3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
               XTEMP(3*(J3-1)+3)=CSMPMAT(3, 1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
                                +CSMPMAT(3, 2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
                                +CSMPMAT(3, 3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
            END DO
!            WRITE(MYUNIT,'(A)') 'potential> XTEMP:'
!            WRITE(MYUNIT,'(3G20.10)') XTEMP(1:3*NATOMS)
            CALL CSMROT(XTEMP, DUMMY, 1, J2)
            CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)+DUMMY(1:3*NATOMS)
         END DO
         CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)/CSMGPINDEX
!
!  Check the CSM for the averaged structure. It should be zero if this structure has the
!  right point group. Need to reset CSMIMAGES and CSMNORM temporarily.
!
         IF (PERMDIST) THEN
            DO J1=1, CSMGPINDEX
               XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
               CALL CSMROT(XTEMP, DUMMY, 1, J1)
               CALL MINPERMDIST(XTEMP, DUMMY, NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, PERIODIC, TWOD, DUMMY2, DIST2, RIGID, RMAT)
               CALL CSMROT(DUMMY, XTEMP,-1, J1) ! need to rotate the permuted rotated images back to the reference orientation
               CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
            END DO
         ELSE
            DO J1=1, CSMGPINDEX
               CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=CSMAV(1:3*NATOMS)
            END DO
         END IF
         AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
         SAVECSMPMAT(1:3, 1:3)=CSMPMAT(1:3, 1:3)
         SAVECSMNORM=CSMNORM
         CSMNORM=0.0D0
         DO J1=1, NATOMS
            CSMNORM=CSMNORM+CSMAV(3*(J1-1)+1)**2+CSMAV(3*(J1-1)+2)**2+CSMAV(3*(J1-1)+3)**2
         END DO
         CSMNORM=2*CSMGPINDEX*CSMNORM
         CALL CSMPOTGRAD(CSMAV, AA, AVVAL,.TRUE., CSMGRAD)
         CSMNORM=SAVECSMNORM
         CSMPMAT(1:3, 1:3)=SAVECSMPMAT(1:3, 1:3)
         CSMIMAGES(1:3*NATOMS*CSMGPINDEX)=SAVECSMIMAGES(1:3*NATOMS*CSMGPINDEX)
         WRITE(MYUNIT,'(A, 2G20.10)') 'potential> CSM values for reference structure and average=', EREAL, AVVAL
      END IF ! DEBUG
      IF (CSMDOGUIDET) THEN ! undo guiding changes
         CSMGP=CSMGPSAVE
         CSMGPINDEX=CSMGPINDEXSAVE
         PTGP(1:3, 1:3, 1:2*CSMGPINDEX)=PTGPSAVE(1:3, 1:3, 1:2*CSMGPINDEX)
         CSMNORM=CSMNORMSAVE
      END IF

   ELSE IF (PERMOPT .OR. PERMINVOPT .OR. DISTOPT) THEN
!  EREAL is the distance in this case
      CALL MINPERMDIST(FIN, X, NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, PERIODIC, TWOD, EREAL, DIST2, RIGID, RMAT)

   ELSE IF (BLJCLUSTER) THEN
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL LJPSHIFTBINC(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (BLJCLUSTER_NOCUT) THEN
! ds656> Binary LJ without cutoff
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL BLJ_CLUST(X, GRAD, EREAL, GRADT)

   ELSE IF (MLJT) THEN
! ds656> Multicomponent LJ with no cutoff (akin to BLJ_CLUST).
!        Should give same results as GLJ, but different code.
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL MLJ(X, GRAD, EREAL, GRADT)

   ELSE IF (GLJT) THEN
! generalised LJ with no cutoff
      CALL RAD(X, GRAD, EREAL, GRADT)
      CALL GLJ(X, GRAD, EREAL, GRADT)

   ELSE IF (BINARY) THEN
      IF (SHIFTCUT) THEN
         CALL LJPSHIFT(X, GRAD, EREAL, GRADT, SECT)
      ELSE
         CALL LJPBIN(X, GRAD, EREAL, GRADT, SECT)
      END IF

   ELSE IF (SOFT_SPHERE) THEN
      CALL SOFT_SPHERE_POT(X, GRAD, EREAL, GRADT, SECT)
   
   ELSE IF (LB2T) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL LB2(X, GRAD, EREAL, GRADT)

   ELSE IF (LJCOULT) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL LJCOUL(X, GRAD, EREAL, GRADT)

   ELSE IF (STOCKT) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL STOCK(X, GRAD, EREAL, GRADT)

! RBAA potentials

   ELSE IF (CAPBINT) THEN
      CALL CAPBIN (X, GRAD, EREAL, GRADT)

   ELSE IF (CHIROT) THEN
      CALL CHIRO (X, GRAD, EREAL, GRADT)

   ELSE IF (DBPT) THEN
      CALL DUMBBELLP (X, GRAD, EREAL, GRADT)

   ELSE IF (DBPTDT) THEN
      CALL DMBLTD (X, GRAD, EREAL, GRADT)

   ELSE IF (DMBLMT) THEN
      CALL DMBLMORSE (X, GRAD, EREAL, GRADT)

   ELSE IF (DMBLPYT) THEN
      CALL DUMBBELLPOLARYUKAWA (X, GRAD, EREAL, GRADT)

   ELSE IF (LINRODT) THEN
      CALL LINROD (X, GRAD, EREAL, GRADT)

   ELSE IF (LWOTPT) THEN
      IF (PERIODIC) THEN
         CALL LWOTP(X, GRAD, EREAL, GRADT)
      ELSE
         CALL LWOTPC(X, GRAD, EREAL, GRADT)
      END IF

   ELSE IF (NCAPT) THEN
      CALL NEWCAPSID (X, GRAD, EREAL, GRADT)

   ELSE IF (NPAHT) THEN
      CALL NEWPAH (X, GRAD, EREAL, GRADT)

   ELSE IF (NTIPT) THEN
      CALL NEWTIP (X, GRAD, EREAL, GRADT)

   ELSE IF (TTM3T) THEN ! Xantheas' TTM3-F water potential
      EREAL=0.0D0
      GRAD(1:3*NATOMS)=0.0D0
!      CALL TTM3FCALL(NATOMS/3, X, EREAL, GRAD)

   ELSE IF (PATCHY) THEN
      CALL PATCHYPOT (X, GRAD, EREAL, GRADT)

   ELSE IF (ASAOOS) THEN ! this is not anisotropic
      CALL ASAOOSPOT (X, GRAD, EREAL, GRADT)

   ELSE IF (GBT) THEN
      CALL GBCALAMITIC (X, GRAD, EREAL, GRADT)

   ELSE IF (GBDT) THEN
      CALL GBDISCOTIC (X, GRAD, EREAL, GRADT)

   ELSE IF (GBDPT) THEN 
      CALL GBDDP (X, GRAD, EREAL, GRADT)

   ELSE IF (GEMT) THEN
      CALL GEM (X, GRAD, EREAL, GRADT)

   ELSE IF (PAHAT) THEN
      CALL PAHA (X, GRAD, EREAL, GRADT)

   ELSE IF (PAHW99T) THEN
      CALL PAHW99 (X, GRAD, EREAL, GRADT)

   ELSE IF (PAPT) THEN
      CALL PAP (X, GRAD, EREAL, GRADT)

   ELSE IF (PAPBINT) THEN
      CALL PAPBIN (X, GRAD, EREAL, GRADT)

   ELSE IF (PAPJANT) THEN
      CALL PAPJANUS (X, GRAD, EREAL, GRADT)

   ELSE IF (PTSTSTT) THEN
      CALL PTSTST (X, GRAD, EREAL, GRADT)

   ELSE IF (MULTPAHAT) THEN
      CALL MULTPAHA (X, GRAD, EREAL, GRADT)

   ELSE IF (PYT .AND. (.NOT.RIGIDCONTOURT)) THEN
      CALL PY(X, GRAD, EREAL, GRADT)

   ELSE IF (MSTBINT) THEN
      CALL MSTBIN (X, GRAD, EREAL, GRADT)

   ELSE IF (MSSTOCKT) THEN
      CALL MULTSTOCK (X, GRAD, EREAL, GRADT)

   ELSE IF (SANDBOXT) THEN
      CALL SANDBOX (X, GRAD, EREAL, GRADT)

   ELSE IF (SILANET) THEN
      CALL SILANE (X, GRAD, EREAL, GRADT)

   ELSE IF (STOCKAAT) THEN
      CALL STOCKAA (X, GRAD, EREAL, GRADT)

   ELSE IF (MORSEDPT) THEN
      CALL DIPOLARMORSE (X, GRAD, EREAL, GRADT)

   ELSE IF (TDHDT) THEN
      CALL TETRAHEDRA (X, GRAD, EREAL, GRADT)

   ELSE IF (DDMT) THEN
      CALL DODECAMORSE (X, GRAD, EREAL, GRADT)

   ELSE IF (MWPOTT .OR. SWPOTT) THEN
      CALL SWTYPE (X, GRAD, EREAL, GRADT)

   ELSE IF (MWFILMT) THEN
      CALL MWFILM (X, NATOMS, GRAD, EREAL, GRADT, BOXLX, BOXLY, BOXLZ, LAT)

   ELSE IF (POLIRT) THEN
      CALL POLIR (X, GRAD, EREAL, GRADT)

   ELSE IF (MBPOLT) THEN
      CALL MBPOL (X, EREAL, GRAD, GRADT)

   ELSE IF (WATERDCT) THEN
      CALL WATERPDC (X, GRAD, EREAL, GRADT)

   ELSE IF (WATERKZT) THEN
      CALL WATERPKZ (X, GRAD, EREAL, GRADT)

   ELSE IF (GAYBERNET) THEN
      CALL GAYBERNE(X, GRAD, EREAL, GRADT,.FALSE.)

   ELSE IF (PYGPERIODICT .OR. PYBINARYT) THEN
      CALL PYGPERIODIC(X, GRAD, EREAL, GRADT)

   ELSE IF (LJCAPSIDT) THEN
      CALL LJCAPSIDMODEL(X, GRAD, EREAL, GRADT)

   ELSE IF (STICKYT) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL STICKY(X, GRAD, EREAL, GRADT,.FALSE.)
!      DIFF=1.0D-3
!      WRITE(*, *) 'analytic and numerical gradients:'
!      DO J1=1, 3*NATOMS
!         X(J1)=X(J1)+DIFF
!         CALL STICKY(X, GRADDUM, EPLUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)-2.0D0*DIFF
!         CALL STICKY(X, GRADDUM, EMINUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)+DIFF
!         IF (GRAD(J1).NE.0.0D0) WRITE(*,'(I5, 3F20.10)') J1, GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), &
!                                100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!      END DO

   ELSE IF (TIP) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL TIP4P(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (QUADT) THEN
      CALL RADR(X, GRAD, EREAL, GRADT)
      CALL QUAD(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (CAPSID) THEN
!      CALL RADC(X, GRAD, EREAL, GRADT)
      CALL FCAPSID(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (STRANDT) THEN
      CALL STRAND(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (PAHT) THEN
      CALL PAH(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (DIFFRACTT) THEN
      WRITE(MYUNIT,'(A)') 'potential> diffract subroutine is commented'
      STOP
!      CALL DIFFRACT(X, GRAD, EREAL, GRADT, SECT, NATOMS)

!      DIFF=1.0D-4
!      WRITE(*, *) 'analytic and numerical gradients:'
!      DO J1=1, 3*NATOMS
!         X(J1)=X(J1)+DIFF
!         CALL DIFFRACT(X, GRADDUM, EPLUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)-2.0D0*DIFF
!         CALL DIFFRACT(X, GRADDUM, EMINUS,.FALSE.,.FALSE.)
!         X(J1)=X(J1)+DIFF
!         IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.-1.0D0)) THEN
!            WRITE(*,'(I5, 2F20.10)') J1, GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!         END IF
!         WRITE(*, '(A, 3G20.10)') 'function: ', EREAL
!         WRITE(*, '(A, 3G20.10)') 'variables: ', X(1:3)
!         WRITE(*, '(A, 3G20.10)') 'gradient:  ', GRAD(1:3)
!      END DO

   ELSE IF (THOMSONT .AND. (.NOT. GTHOMSONT)) THEN
      IF (ODDCHARGE == 1.0D0) THEN
         CALL THOMSON(X, GRAD, EREAL, GRADT)
      ELSE
         CALL ODDTHOMSON(X, GRAD, EREAL, GRADT)
      END IF

   ELSE IF (THOMSONT .AND. GTHOMSONT) THEN
      CALL GTHOMSON(X, GRAD, EREAL, GRADT)

   ELSE IF (GAUSST) THEN
      CALL GFIELD(X, GRAD, EREAL, GRADT)

   ELSE IF (QDT) THEN
      CALL QDTEST(X, GRAD, EREAL, GRADT)

   ELSE IF (QD2T) THEN
      CALL QDTEST2(X, GRAD, EREAL, GRADT)

   ELSE IF (MULLERBROWNT) THEN
      CALL MB(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (JMT) THEN
      CALL JMEC(NATOMS, X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (DF1T) THEN
      CALL DF1GRAD(X, NATOMS, GRAD, EREAL, GRADT, BOXLX, BOXLY)

   ELSE IF (ONEDAPBCT) THEN
      CALL ENERGY_1D_APBC(X, GRAD, EREAL, GRADT, SECT)

   ELSE IF (ONEDPBCT) THEN
      WRITE(*, *) 'This potential option is commented out.'
      STOP
!      CALL ENERGY_1D_PBC(X, GRAD, EREAL, GRADT, SECT)
   ELSE IF (TWODAPBCT) THEN
      WRITE(*, *) 'This potential option is commented out.'
      STOP
!      CALL ENERGY_2D_APBC(X, GRAD, EREAL, GRADT, SECT)            
   ELSE IF (TWODPBCT) THEN
      WRITE(*, *) 'This potential option is commented out.'
      STOP
!      CALL ENERGY_2D_PBC(X, GRAD, EREAL, GRADT, SECT)
   ELSE IF (THREEDAPBCT) THEN
      WRITE(*, *) 'This potential option is commented out.'
      STOP
!      CALL ENERGY_3D_APBC(X, GRAD, EREAL, GRADT, SECT)
   ELSE IF (THREEDPBCT) THEN
      WRITE(*, *) 'This potential option is commented out.'
      STOP
!      CALL ENERGY_3D_PBC(X, GRAD, EREAL, GRADT, SECT)

   ELSE
! Otherwise, assume Lennard-Jones type particles.

! RAD must be called before the routine that calculates the potential or LBFGS
! will get confused even if EVAP is set .TRUE. correctly.
      CALL RAD(X, GRAD, EREAL, GRADT)
      IF (EVAPREJECT) RETURN
      IF (CUTT) THEN
         CALL LJCUT(X, GRAD, EREAL, GRADT, SECT)
      ELSE
         CALL LJ(X, GRAD, EREAL, GRADT, SECT)
         IF (AXTELL) CALL AXT(NATOMS, X, GRAD, EREAL, GRADT, ZSTAR)
      END IF
   END IF
!
!  --------------- End of possible potentials - now add fields if required ------------------------------
!

   IF (RIGIDCONTOURT) THEN
!       sf344> sanity check, since this keyword works only for a pair of rigid particles        
      IF (NATOMS.NE.4) THEN
         WRITE(MYUNIT,'(A)') 'potential> ERROR - RIGIDCONTOUR works only for a PAIR of rigid bodies, exiting'
         STOP
      ELSE
         CALL RIGIDCONTOUR()
         WRITE(MYUNIT,'(A)') 'potential> finished generating the potential energy matrix, exiting GMIN'
         STOP
      END IF
   END IF 

   IF (PULLT) THEN
      EREAL=EREAL - PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
      GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3) - PFORCE
      GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3) + PFORCE
   END IF

   IF (GCBHT) THEN
! Need a negative sign because the exponential is -(V-N mu)/kT
      EREAL=EREAL - NATOMS*GCMU
   END IF

   IF (FIELDT) THEN
      IF (CENT) THEN
         CALL FDM(X, GRAD, EREAL, GRADT)
      ELSE
         CALL FD(X, GRAD, EREAL, GRADT)
      END IF
   ELSE IF (CPMD .AND. (NPCALL == 1)) THEN
      CALL SYSTEM(' sed -e "s/DUMMY/RESTART WAVEFUNCTION GEOFILE LATEST/" ' //  SYS(1:LSYS) // ' > temp ')
      CALL SYSTEM(' mv temp ' // SYS(1:LSYS) // '.restart')
   ELSE IF (CPMD .AND. (.NOT. SCT)) THEN
      INQUIRE(FILE='RESTART.1', EXIST=YESNO)
      OPEN(UNIT=8, FILE='newgeom', STATUS='UNKNOWN')
      DO J1=1, NATOMS
         WRITE(8,'(6F20.10)') X(3*(J1-1)+1), X(3*(J1-1)+2), X(3*(J1-1)+3), 0.0D0, 0.0D0, 0.0D0
      END DO
      CLOSE(8)
      CALL SYSTEM(' mv newgeom GEOMETRY ')
      CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out >& /dev/null ')
      IF (.NOT. YESNO) THEN
         CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      ELSE
         CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // '.restart > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      END IF
      CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) // &
                  '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
      OPEN (UNIT=7, FILE='temp', STATUS='OLD')
      READ(7,'(A)') FNAME
      WRITE(MYUNIT,'(A)') FNAME
      CLOSE(7)
      OPEN (UNIT=7, FILE='ENERGY', STATUS='OLD')
      READ(7,*) EREAL, GEMAX
      CLOSE(7)
      IF (GEMAX > 1.0D-5) WRITE(MYUNIT,'(A, G15.5, A)') 'WARNING, GEMAX=', GEMAX,' CPMD wavefunction convergence suspect'
      OPEN(UNIT=7, FILE='GEOMETRY', STATUS='OLD')
      DO J1=1, NATOMS
         READ(7,*) GEMAX, GEMAX, GEMAX, GRAD(3*(J1-1)+1), GRAD(3*(J1-1)+2), GRAD(3*(J1-1)+3)
!         WRITE(*,'(6F20.10)') GEMAX, GEMAX, GEMAX, GRAD(3*(J1-1)+1), GRAD(3*(J1-1)+2), GRAD(3*(J1-1)+3)
         GRAD(3*J1-2:3*J1)=-GRAD(3*J1-2:3*J1)
      END DO
      CLOSE(7)
   ELSE IF (DL_POLY) THEN
      CALL SYSTEM('DLPOLY.X > forces')
      OPEN (UNIT=91, FILE='STATIS', STATUS='OLD')
      READ(91,*) DUMM
      READ(91,*) DUMM
      READ(91,*) DUMM
      READ(91,*) EREAL
      WRITE(MYUNIT,'(A, G20.10)' ) 'EREAL=', EREAL
      CLOSE(91)
      OPEN (UNIT=91, FILE='STATIS', STATUS='OLD')
      READ(91,*) (GRAD(J3), J3=1, 3*NATOMS)
      CLOSE(91)
      WRITE(91,'(3F20.10)') (GRAD(J3), J3=1, 3*NATOMS)
   END IF

!   IF (SQUEEZET) CALL SQUEEZE(X, GRAD, EREAL, GRADT)

   IF (K_COMP == 0.0D0) COMPON=.FALSE.
   IF (COMPON) CALL COMPRESS(X, GRAD, EREAL, GRADT)

! Apply twist forces, if appropriate
   IF (TWISTT) THEN
      CALL TWIST(X(1:3*NATOMS), NATOMS, GRAD(1:3*NATOMS), EREAL, GRADT)
   END IF

! ds656> Apply a substrate field for keyword MIE_FIELD
   IF(MIEFT) CALL MIEF(X,GRAD,EREAL,GRADT,.FALSE.)

!
! js850> apply harmonic field for keyword HARMONICF
!
   IF ( HARMONICF ) THEN
      DO J1=1, NATOMS
         IF (HARMONICFLIST(J1)) THEN
!            WRITE(*,'(1x, 6F25.16)') X(1), X(2), X(3), X0(1), X0(2), X0(3)
!            WRITE(*,*) "#V(1) ", GRAD((J1-1)*3+1)
            CALL HARMONICFIELD( X(3*J1-2), HARMONICR0(3*J1-2), GRAD(3*J1-2), EREAL)
!            WRITE(*,*) "#V(1) after ", GRAD((J1-1)*3+1)
         END IF
      END DO
   END IF

   IF (GRADT .OR. CSMT) THEN
      IF (FREEZE) THEN
         DO J1=1, NATOMS
            IF (FROZEN(J1)) THEN
               GRAD(3*J1-2 : 3*J1)=0.0D0
            END IF
         END DO
      END IF

      IF (SEEDT .AND. FREEZECORE) THEN
         GRAD(3*(NATOMS-NSEED)+1 : 3*NATOMS)=0.0D0
      END IF

!  Preserve centre of mass if required.
!      IF (CENT.AND.(.NOT.SEEDT)) THEN
!         XG=0.0D0
!         YG=0.0D0
!         ZG=0.0D0
!         DO J1=1, NATOMS-NSEED
!            J2=3*J1
!            XG=XG+GRAD(J2-2)
!            YG=YG+GRAD(J2-1)
!            ZG=ZG+GRAD(J2)
!         END DO
!         XG=XG/(NATOMS-NSEED)
!         YG=YG/(NATOMS-NSEED)
!         ZG=ZG/(NATOMS-NSEED)
!         DO J1=1, NATOMS-NSEED
!            J2=3*J1
!            GRAD(J2-2)=GRAD(J2-2)-XG
!            GRAD(J2-1)=GRAD(J2-1)-YG
!            GRAD(J2)=  GRAD(J2)-ZG
!         END DO
!      END IF

      IF (CSMT .AND. (.NOT.SYMMETRIZECSM)) THEN
         DUMMY2=0.0D0
         RMS=0.0D0
      ELSE IF (AACONVERGENCET .AND. (ATOMRIGIDCOORDT)) THEN
         DUMMY2=SUM(GRAD(1:DEGFREEDOMS)**2)
         RMS=MAX(SQRT(DUMMY2/DEGFREEDOMS), 1.0D-100)
         IF (RMS < 5.0D0 * BQMAX) THEN
            CALL AACONVERGENCE (GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD, RMS)           
         END IF
      ELSE IF (.NOT.THOMSONT) THEN
         DUMMY2=SUM(GRAD(1:3*NATOMS)**2)
         RMS=MAX(DSQRT(DUMMY2/(3*NATOMS)), 1.0D-100)
      ELSE
         DUMMY2=SUM(GRAD(1:2*NATOMS)**2)
         RMS=MAX(DSQRT(DUMMY2/(2*NATOMS)), 1.0D-100)
      END IF
      IF(DEBUG.AND.(RMS.NE.RMS)) THEN
         WRITE(MYUNIT,'(A)' ) 'potential> WARNING - RMS force is NaN - if using AMBER igb=1, can be due to negative Born radii'
      END IF

!   IF ((CSMRMS < GUIDECUT).AND.(CSMDOGUIDET)) THEN
      IF (CSMDOGUIDET) THEN
         IF (DEBUG) WRITE(MYUNIT,'(A)') 'potential> switching off guiding point group'
         WRITE(MYUNIT,'(A, F20.10, A, I5, A, G12.5, A, G20.10, A, F11.1)') 'Qu guide      E=', &
                                                                  EREAL,' steps=', CSMIT,' RMS=', CSMRMS
         GUIDECHANGET=.TRUE.
         CSMDOGUIDET=.FALSE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
         GOTO 10
      END IF

      IF (RMS < GUIDECUT) THEN
! Guiding potentials are set back to true in quench.f
         IF (DFTBCT .AND. LJATT) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'potential> switching off LJAT - rescaling distances for DFTB carbon'
            LJATT=.FALSE.
            GUIDECHANGET=.TRUE.
            X(1:3*NATOMS)=X(1:3*NATOMS)*LJATTOC
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF (NATBT .AND. GUPTAT) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Gupta'
            GUPTAT=.FALSE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF (CPMD .AND. SCT) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Sutton-Chen'
            SCT=.FALSE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF (WELCH .AND. TOSI) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Tosi'
            TOSI=.FALSE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF ((ZETT1 .OR. ZETT2) .AND. MORSET) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off MORSE'
            MORSET=.FALSE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF (PACHECO .AND. (.NOT.AXTELL)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching on AXTELL'
            AXTELL=.TRUE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         ELSE IF (PERCOLATET .AND. COMPON) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A, G20.10, A, G20.10)' ) 'Switching off compression, GUIDECUT=', &
                                                               GUIDECUT, '; RMS=', RMS
            COMPON=.FALSE.
            GUIDECHANGET=.TRUE.
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
            GOTO 10
         END IF
      END IF
   END IF ! GRADT .OR. CSMT

END SUBROUTINE POTENTIAL
