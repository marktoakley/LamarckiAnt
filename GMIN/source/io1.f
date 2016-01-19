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
      SUBROUTINE IO1
      USE commons
      USE genrigid
      USE modamber
      use qmodule
      use permu
      USE modcharmm
      USE porfuncs
      use bgupmod
      IMPLICIT NONE

      DOUBLE PRECISION VECMN, DUMMY
      LOGICAL END, YESNO, EXISTS
      INTEGER J1, J2, JP, NTYPE(105), ISTAT, NSTART, NFINISH, NDUMMY, J3
      INTEGER MP, LP
      INTEGER Iostatus
      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      CHARACTER FNAME*80, TSTRING*80
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

      INTEGER  NQTOT
      COMMON /TOT/ NQTOT
      DOUBLE PRECISION EPS2, RAD, HEIGHT
      COMMON /CAPS/ EPS2, RAD, HEIGHT
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, LAST
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT

      INTEGER COORDS_UNIT 

      IF (DL_POLY) THEN
         OPEN(UNIT=91,FILE='CONFIG',STATUS='OLD')
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
13       READ(91,'(A1)',END=14) DUMMY
         NATOMS=NATOMS+1
         READ(91,*) COORDS(3*(NATOMS-1)+1,1),COORDS(3*(NATOMS-1)+2,1),COORDS(3*(NATOMS-1)+3,1)
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
!        WRITE(MYUNIT,'(3G20.10)') COORDS(3*(NATOMS-1)+1,1),COORDS(3*(NATOMS-1)+2,1),COORDS(3*(NATOMS-1)+3,1)
         GOTO 13
14       CONTINUE
         CLOSE(91)
      ELSEIF (AMHT) THEN
          write(MYUNIT,'(A)')'DUMMY    '
      ELSEIF (.NOT.(AMBER12T.OR.AMBERT.OR.AMBER.OR.CPMD.OR.CHRMMT.OR.DMACRYST.OR.USERPOTT)) THEN
! jwrm2> Added coordsini for GTHOMSON
         IF (GTHOMSONT) THEN
           CALL FILE_OPEN('coordsini', COORDS_UNIT, .FALSE.)
         ELSE
           CALL FILE_OPEN('coords', COORDS_UNIT, .FALSE.)
         END IF
!         CLOSE(7)
!         OPEN(UNIT=7,FILE='coords',STATUS='OLD')
!js850>  this step (NATOMS=NATOMS/NPAR) is moved to earlier in the program (to
!COUNTATOMS), so that the allocation is done with the correct NATOMS
!        IF (MOD(NATOMS,NPAR).NE.0) THEN
!           WRITE(MYUNIT,'(A,I8,A,I8)') 'Number of atoms in coords file=',NATOMS,
!    &               ' is not divisible by number of runs=',NPAR
!           STOP
!        ENDIF
!        NATOMS=NATOMS/NPAR
         IF (PERMOPT.OR.PERMINVOPT) THEN
            OPEN(UNIT=17,FILE='finish',STATUS='OLD')
            READ(17,*) (FIN(J1),J1=1,3*NATOMS)
            WRITE(MYUNIT,'(A)') 'Target coordinates read from file finish'
         ENDIF
         IF (TSALLIST.AND.(QTSALLIS.EQ.0)) QTSALLIS=1.0D0+1.0D0/(3*NATOMS)
         IF (DFTBT.OR.TOSI.OR.WELCH) THEN
!           REWIND(7) ! this seems to be needed now?
            DO JP=1,NPAR
               DO J1=1,NATOMS
                  CALL INPUT(END, COORDS_UNIT)
                  J2=3*(J1-1)
                  CALL READA(ZSYM(J1))
                  CALL MYUPCASE(ZSYM(J1))
                  IF (ZSYM(J1).EQ.'H ') IATNUM(J1+1)=1
                  IF (ZSYM(J1).EQ.'HE') IATNUM(J1+1)=2
                  IF (ZSYM(J1).EQ.'LI') IATNUM(J1+1)=3
                  IF (ZSYM(J1).EQ.'BE') IATNUM(J1+1)=4
                  IF (ZSYM(J1).EQ.'B ') IATNUM(J1+1)=5
                  IF (ZSYM(J1).EQ.'C ') IATNUM(J1+1)=6
                  IF (ZSYM(J1).EQ.'N ') IATNUM(J1+1)=7
                  IF (ZSYM(J1).EQ.'O ') IATNUM(J1+1)=8
                  IF (ZSYM(J1).EQ.'F ') IATNUM(J1+1)=9
                  IF (ZSYM(J1).EQ.'S ') IATNUM(J1+1)=18
                  CALL READF(COORDS(J2+1,JP))
                  CALL READF(COORDS(J2+2,JP))
                  CALL READF(COORDS(J2+3,JP))
               ENDDO
            ENDDO
         ELSE
            rewind(COORDS_UNIT)
            DO JP=1,NPAR
               DO J1=1,NATOMS
                  J2=3*(J1-1)
                   READ(COORDS_UNIT,*) COORDS(J2+1,JP), COORDS(J2+2,JP), COORDS(J2+3,JP)
               ENDDO
            ENDDO
         ENDIF
         CLOSE(COORDS_UNIT)
!      ELSE IF (AMBERT) THEN
!         DO JP=1,NPAR
!             COORDS(:,JP)=COORDS(:,1)   ! we have already read the coordinates for AMBER runs into this array
!         END DO
      ENDIF

      IF (CPMD) THEN
         FNAME=SYS(1:LSYS)
         WRITE(MYUNIT,'(A,A)') ' Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
111      READ(7,'(A)') FNAME
         IF (FNAME(1:6).EQ.'&ATOMS') THEN
            J1=0
121         READ(7,'(A)') TSTRING
            IF (TSTRING(1:1).EQ.'*') THEN
               J1=J1+1
               READ(7,'(A)') FNAME
               READ(7,*) NTYPE(J1)
               DO J2=1,NTYPE(J1)
                  IATNUM(NATOMS+J2)=1
                  ZSYM(NATOMS+J2)='CP'
                  READ(7,*) COORDS(3*(NATOMS+J2-1)+1,1),COORDS(3*(NATOMS+J2-1)+2,1),COORDS(3*(NATOMS+J2-1)+3,1)
               ENDDO
               NATOMS=NATOMS+NTYPE(J1)
               GOTO 121
            ELSE IF (TSTRING(1:1).EQ.' ') THEN
               GOTO 121
            ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               GOTO 131
            ENDIF
         ELSE
            GOTO 111
         ENDIF

131      CONTINUE
         CLOSE(7)

         CALL SYSTEM(' grep -c ANGSTROM ' // SYS(1:LSYS) // ' > temp')
         OPEN(UNIT=7,FILE='temp',STATUS='OLD')
         READ(7,*) J1
         CLOSE(7)
         IF (J1.EQ.1) THEN
            WRITE(MYUNIT,'(A)') ' Converting initial coordinates from Angstrom to Bohr'
            DO J1=1,3*NATOMS
               COORDS(J1,1)=COORDS(J1,1)*1.889726164D0
            ENDDO
         ENDIF
      ENDIF

      IF (PERMOPT) WRITE(MYUNIT,'(A)') 'Minimising distance with respect to permutations'
      IF (PERMINVOPT) WRITE(MYUNIT,'(A)') 'Minimising distance with respect to permutation-inversions'
      IF (RESIZET) THEN
         WRITE(MYUNIT,'(A,F15.7)') 'Multiplying coordinates by ',RESIZE
         DO JP=1,NPAR
            DO J1=1,3*NATOMS
               COORDS(J1,JP)=COORDS(J1,JP)*RESIZE
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.SEEDT.AND..NOT.AMHT.AND..NOT.SUPPRESST) THEN
         WRITE(MYUNIT,20) 
20       FORMAT('Initial coordinates:')
         IF (MPIT) THEN
            WRITE(MYUNIT,30) (COORDS(J1,MYNODE+1),J1=1,3*NATOMS)
         ELSE 
           DO JP=1,NPAR
               WRITE(MYUNIT,30) (COORDS(J1,JP),J1=1,3*NATOMS)
30             FORMAT(3F20.10)
            ENDDO
         ENDIF
      ENDIF

      IF (MSORIGT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' M and S silicon atoms'
      ELSE IF (MSTRANST) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' M and S transferable silicon atoms'
      ELSE IF (FRAUSIT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Frauenheim silicon atoms'
      ELSE IF (MORSET) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Morse atoms with rho=',RHO
      ELSE IF (LB2T) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' LB2 atoms'
      ELSE IF (NEON) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' neon atoms'
      ELSE IF (ARGON) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' argon atoms'
      ELSE IF (SCT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' SC atoms'
         WRITE(MYUNIT,'(A,I2,A,I2,A,F12.6,A,F12.6,A,F12.6)') 'n=',NN,' m=',MM,' epsilon=',SCEPS,' sigma=',SIG,' C=',SCC
      ELSE IF (DFTBT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Tiffany atoms'
      ELSE IF (SW) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Stillinger-Weber Si atoms'
      ELSE IF (QUIPT) THEN
         WRITE(MYUNIT,'(I4,A,A,A)') NATOMS, ' ', QUIPATOMTYPE, 'atoms'
      ELSE IF (FAL) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Farkas Al atoms'
         CALL ALINIT
      ELSE IF (FNI) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Farkas Ni atoms'
         CALL NIINIT
      ELSE IF (RGCL2) THEN
         WRITE(MYUNIT,'(A,I3,A)') 'DIM potential for Cl_2 with ',NATOMS,' rare gas atoms'
         IF (AXTELL) WRITE(MYUNIT,'(A,F12.5)') 'Including Axilrod-Teller term with Z/hartree Angstrom^9=',ZSTAR
      ELSE IF (ARNO) THEN
         WRITE(MYUNIT,'(A,I3,A)') 'DIM potential for NO with ',NATOMS,' Ar atoms'
         IF (AXTELL) WRITE(MYUNIT,'(A,F12.5)') 'Including Axilrod-Teller term with Z/hartree Angstrom^9=',ZSTAR
      ELSE IF (TOSI) THEN
         WRITE(MYUNIT,'(A,I3,A)') ' Born-Mayer potential for ',NATOMS,' ions with Tosi-Fumi parameters:'
         WRITE(MYUNIT,'(A,F13.8,A,F13.8,A,F13.8,A,F13.8)') ' A++=',APP,' A--=',AMM,' A+-=',APM,' rho=',RHO
      ELSE IF (WELCH) THEN
         WRITE(MYUNIT,'(A,I3,A)') ' Welch potential for ',NATOMS,' ions with parameters:'
         WRITE(MYUNIT,'(A,F13.8,A,F13.8,A,F13.8,A,F13.8)') 
     1           ' A++=',APP,' A--=',AMM,' A+-=',APM,' rho=',RHO
         WRITE(MYUNIT,'(A,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     1           ' Q+=',XQP,' Q-=',XQM,' alpha+=',ALPHAP,' alpha-=',ALPHAM
      ELSE IF (P46) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' 3-colour, 46 bead model polypeptide'
      ELSE IF (BLNT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' bead BLN model'
      ELSE IF (AMBER) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' AMBER atoms'
         IF (AMCUT) WRITE(MYUNIT,'(A,2F12.2)') 'Cutoffs/Angstrom=',RCUTOFF,QCUTOFF
      ELSE IF (AMBERT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' AMBER9 atoms'
!hk286
         IF (RIGIDINIT) THEN
            WRITE(MYUNIT, '(I6,A)') NRIGIDBODY,' RIGID BODIES'
            WRITE(MYUNIT, '(I6,A)') (DEGFREEDOMS-6*NRIGIDBODY)/3, ' SINGLE ATOMS'
         ENDIF
!hk286
      ELSE IF (AMBER12T .AND. (.NOT. CUDAT)) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' AMBER12 atoms'
      ELSE IF (SOFT_SPHERE) THEN !js850>
         IF (NTYPEA.GT.NATOMS) THEN
            WRITE(MYUNIT,'(A, 2G20.10)')'Error: NTYPEA, NATOMS=',NTYPEA, NATOMS
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,I4,A,I4,A)') ' Binary soft sphere solid: ',NTYPEA,' A atoms and ',NATOMS-NTYPEA,' B atoms'
         IF (SHIFTCUT) WRITE(MYUNIT,'(A,F15.5)') ' Shifted, truncated potential, cutoff=',CUTOFF
      ELSE IF (BINARY) THEN
         IF (NTYPEA.GT.NATOMS) THEN
            WRITE(MYUNIT,'(A, 2G20.10)')'Error: NTYPEA, NATOMS=',NTYPEA, NATOMS
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,I4,A,I4,A)') ' Binary Lennard-Jones solid: ',NTYPEA,' A atoms and ',NATOMS-NTYPEA,' B atoms'
         WRITE(MYUNIT,'(A,F10.4,A,F10.4,A,F10.4,A,F10.4)') 
     1     ' eps(AA)=1, sigma(AA)=1, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,' sigma(AB)=',SIGAB,' sigma(BB)=',SIGBB
         IF (SHIFTCUT) WRITE(MYUNIT,'(A,F15.5,A)') ' Shifted, truncated potential, cutoff=',CUTOFF,' in sigma units'
      ELSE IF (BLJCLUSTER) THEN
         IF (NTYPEA.GT.NATOMS) THEN
            WRITE(MYUNIT,'(A, 2G20.10)') 'Error: NTYPEA, NATOMS=',NTYPEA, NATOMS
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,I4,A,I4,A)') ' Binary Lennard-Jones cluster: ',NTYPEA,' A atoms and ',NATOMS-NTYPEA,' B atoms'
         WRITE(MYUNIT,'(A,F10.4,A,F10.4,A,F10.4,A,F10.4)')
     1     ' eps(AA)=1, sigma(AA)=1, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,' sigma(AB)=',SIGAB,' sigma(BB)=',SIGBB
         WRITE(MYUNIT,'(A,F15.5,A)') ' Shifted, truncated potential, cutoff=',CUTOFF,' in sigma units'

      ELSE IF (BLJCLUSTER_NOCUT) THEN
         IF (NTYPEA.GT.NATOMS) THEN
            WRITE(MYUNIT,'(A, 2G20.10)') 'Error: NTYPEA, NATOMS=',NTYPEA, NATOMS
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,I4,A,I4,A)') ' Binary Lennard-Jones cluster: ',NTYPEA,' A atoms and ',NATOMS-NTYPEA,' B atoms'
         WRITE(MYUNIT,'(A,F10.4,A,F10.4,A,F10.4,A,F10.4)')
     1     ' eps(AA)=1, sigma(AA)=1, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,' sigma(AB)=',SIGAB,' sigma(BB)=',SIGBB
         WRITE(MYUNIT,'(A)') ' No cutoff.'
!
! Sanity check for atoms of different elements and masses. If symmetrize is set
! then these must be consistent or we may permute atoms of different elements
! and change the energy in the Markov chain accidentally.
!
         IF(.NOT. SPECMASST .AND. NTYPEA < NATOMS) THEN !ds656> Check when SPECMASST=.FALSE.
            DO J1=2,NTYPEA
               IF (ABS(ATMASS(J1)-ATMASS(1)).GT.1.0D-10) THEN
                  WRITE(MYUNIT,'(A,I6,A,2G20.10)') 'io1> WARNING *** atom masses for 1 and ',J2,' are ',ATMASS(1),ATMASS(J2)
               ENDIF
            ENDDO
            IF (ABS(ATMASS(1)-ATMASS(NTYPEA+1)).LT.1.0D-10) THEN
               WRITE(MYUNIT,'(A,2G20.10)') 'io1> WARNING *** atom masses for A and B are ',ATMASS(1),ATMASS(NTYPEA+1)
               IF (SYMMETRIZE.AND.READMASST) THEN
                  WRITE(MYUNIT,'(A)') 'io1> ERROR *** atom masses are incompatible with SYMMETRIZE keyword'
                  STOP
               ENDIF
            ENDIF
            DO J1=NTYPEA+1,NATOMS
               IF (ABS(ATMASS(J1)-ATMASS(NTYPEA+1)).GT.1.0D-10) THEN
                  WRITE(MYUNIT,'(A,I6,A,2G20.10)') 'io1> WARNING *** atom masses for NTYPEA+1 and ',J2,' are ',
     &                 ATMASS(NTYPEA+1),ATMASS(J2)
               ENDIF
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'Atom masses will be set to species indices.'
         ENDIF !... <ds656
      ELSE IF (GLJT) THEN
! 
! Sanity check for atoms of different elements and masses. If symmetrize is set
! then these must be consistent or we may permute atoms of different elements
! and change the energy in the Markov chain accidentally.
! Not tested yet for GLJ and MLJ! DJW
!                 
         IF(.NOT. SPECMASST) THEN !ds656> check when SPECMASST=.FALSE.
            NSTART=1
            NFINISH=NSPECIES(1)
            DO J2=1,NSPECIES(0)
               DO J1=NSTART+1,NFINISH ! consistency check within set J2
                  IF (ABS(ATMASS(J1)-ATMASS(NSTART)).GT.1.0D-10) THEN
                     WRITE(MYUNIT,'(A,I6,A,2G20.10)') 'io1> WARNING *** atom masses for ',NSTART,' and ',J1,' are ',
     &                    ATMASS(NSTART),ATMASS(J1)
                  ENDIF
               ENDDO
               NDUMMY=NFINISH
               DO J3=J2+1,NSPECIES(0) ! check for mass differences with other species
                  IF (ABS(ATMASS(NSTART)-ATMASS(NDUMMY+1)).LT.1.0D-10) THEN
                     WRITE(MYUNIT,'(A,2G20.10)') 'io1> WARNING *** atom masses for A and B are ',ATMASS(NSTART),ATMASS(NDUMMY+1)
                     IF (SYMMETRIZE.AND.READMASST) THEN
                        WRITE(MYUNIT,'(A)') 'io1> ERROR *** atom masses are incompatible with SYMMETRIZE keyword'
                        STOP
                     ENDIF
                  ENDIF
                  NDUMMY=NDUMMY+NSPECIES(J3)
               ENDDO
               NSTART=NSTART+NSPECIES(J2)
               NFINISH=NFINISH+NSPECIES(J2)
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'Atom masses will be set to species indices.'
         ENDIF !...<ds656
      ELSE IF(BGUPTAT) THEN
         IF (NTYPEA.GT.NATOMS) THEN
            WRITE(MYUNIT,'(A, 2G20.10)') 'Error: NTYPEA, NATOMS=',NTYPEA, NATOMS
            STOP
         ENDIF
         WRITE(MYUNIT,'(A,I4,A,I4,A)') 'Binary Gupta:  ',NTYPEA,'A ATOMS
     1    and ',NATOMS-NTYPEA,' B ATOMS '
         WRITE(MYUNIT,'(A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,
     1         A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,
     2         A,F10.4,A,F10.4)')
     3         'AAA=' ,AAA, 'AAB=' ,AAB, 'ABB=' ,ABB, 'PAA=' ,PAA,
     4         'PAB=' ,PAB, 'PBB=' ,PBB, 'QAA=' ,QAA, 'QAB=' ,QAB,
     5         'QBB=' ,QBB, 'ZAA=' ,ZAA, 'ZAB=' ,ZAB, 'ZBB=' ,ZBB,
     6         'R0AA=' ,R0AA, 'R0AB=' ,R0AB, 'R0BB=' ,R0BB
!
! Sanity check for atoms of different elements and masses. If symmetrize is set
! then these must be consistent or we may permute atoms of different elements
! and change the energy in the Markov chain accidentally.
!
         IF(.NOT. SPECMASST .AND. NTYPEA < NATOMS) THEN !ds656> check when SPECMASST=.FALSE.
            DO J1=2,NTYPEA
               IF (ABS(ATMASS(J1)-ATMASS(1)).GT.1.0D-10) THEN
                  WRITE(MYUNIT,'(A,I6,A,2G20.10)') 'io1> WARNING *** atom masses for 1 and ',J1,' are ',ATMASS(1),ATMASS(J1)
               ENDIF
            ENDDO
            IF (ABS(ATMASS(1)-ATMASS(NTYPEA+1)).LT.1.0D-10) THEN
               WRITE(MYUNIT,'(A,2G20.10)') 'io1> WARNING *** atom masses for A and B are ',ATMASS(1),ATMASS(NTYPEA+1)
               IF (SYMMETRIZE.AND.READMASST) THEN
                  WRITE(MYUNIT,'(A)') 'io1> ERROR *** atom masses are incompatible with SYMMETRIZE keyword'
                  STOP
               ENDIF
            ENDIF
            DO J1=NTYPEA+1,NATOMS
               IF (ABS(ATMASS(J1)-ATMASS(NTYPEA+1)).GT.1.0D-10) THEN
                  WRITE(MYUNIT,'(A,I6,A,2G20.10)') 'io1> WARNING *** atom masses for NTYPEA+1 and ',J1,' are ',
     &                 ATMASS(NTYPEA+1),ATMASS(J1)
               ENDIF
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'Atom masses will be set to species indices.'
         ENDIF !...<ds656
      ELSE IF (BASWAP) THEN
         WRITE(MYUNIT,'(A,I6)') 'Atom identity swaps will be performed. NWAIT= ', BASWAP_NWAIT
         WRITE(MYUNIT,'(A,F10.6)') 'Maximum fraction of swap moves between A and
     1     B type particles:  ', BASWAP_FRAC
      ELSE IF (DZTEST) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Dzugutov atoms'
      ELSE IF (ZETT1) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Zetterling type 1 atoms'
      ELSE IF (ZETT2) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Zetterling type 2 atoms'
      ELSE IF (PACHECO) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Pacheco-Ramalho C60 molecules'
      ELSE IF (MODEL1T) THEN
         WRITE(MYUNIT,'(A)') ' 1D landscape model 1'
      ELSE IF (EAMLJT) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Baskes EAMLJ atoms'
      ELSE IF (PBGLUET) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Lead glue atoms'
      ELSE IF (ACKLANDT) THEN
         IF (ACKLANDID.GT.0) THEN
            WRITE(MYUNIT,'(I4,A,I4)') NATOMS,' Ackland atoms, periodic, Ackland id is ',ACKLANDID
         ELSE
            WRITE(MYUNIT,'(I4,A,I4)') NATOMS,' Ackland atoms, cluster, Ackland id is ',ACKLANDID
         ENDIF
      ELSE IF (ALGLUET) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Aluminium glue atoms'
      ELSE IF (MGGLUET) THEN
         WRITE(MYUNIT,'(I4,A,F10.4)') NATOMS,' Magnesium glue atoms'
      ELSE IF (FST) THEN
         WRITE(MYUNIT,'(I4,A,I2,A)') NATOMS,' Finnis-Sinclair (TYPE ', GATOM,') atoms'
      ELSE IF (GUPTAT) THEN
         WRITE(MYUNIT,'(I4,A,I2,A)') NATOMS,' Gupta (TYPE ', GATOM,') atoms'
         IF (NATBT) WRITE(MYUNIT,'(A)') ' guiding potential for Na tight-binding'
      ELSE IF (QUADT) THEN
         WRITE(*,'(I4,A,I4,A,I1,A)') 3*NATOMS,' rigid body coordinates will be optimised for ',NATOMS/2,' quadrupoles'
      ELSE IF (TIP) THEN
         WRITE(MYUNIT,'(I4,A,I4,A,I1,A)') 3*NATOMS,' rigid body coordinates will be optimised for ',NATOMS/2,
     1        ' TIP',TIPID,' waters'
         IF (SORTT) THEN
            WRITE(MYUNIT,'(A)') 'Turning off SORT option for TIP'
            SORTT=.FALSE.
         ENDIF
      ELSE IF (CAPSID) THEN
         WRITE(MYUNIT,'(I4,A,I4,A)') 3*NATOMS,' rigid body coordinates will be optimised for ',NATOMS/2,' capsid pentamers'
         WRITE(MYUNIT,'(4(A,F12.8))') ' rho=',RHO,' repulsive site epsilon=',EPS2,' pentamer radius=',RAD,
     1                                        ' height=',HEIGHT
         IF (SORTT) THEN
            WRITE(MYUNIT,'(A)') 'Turning off SORT option for capsid'
            SORTT=.FALSE.
         ENDIF
      ELSE IF (STRANDT) THEN
         WRITE(MYUNIT,'(I4,A,I4,A)') 3*NATOMS,' rigid body coordinates will be optimised for ',NATOMS/2,' beta strands'
         IF (SORTT) THEN
            WRITE(MYUNIT,'(A)') 'Turning off SORT option for strands'
            SORTT=.FALSE.
         ENDIF
      ELSE IF (PAHT) THEN
         WRITE(MYUNIT,'(I4,A,I4,A)') 3*NATOMS,' rigid body coordinates will be optimised for ',NATOMS/2,' PAH molecules'
         IF (SORTT) THEN
            WRITE(MYUNIT,'(A)') 'Turning off SORT option for strands'
            SORTT=.FALSE.
         ENDIF
      ELSE IF (STOCKT) THEN
         WRITE(MYUNIT,'(I4,A,F15.5,A,F15.5)') NATOMS/2,' Stockmayer atoms with mu=',STOCKMU,' and lambda=',STOCKLAMBDA
      ELSE IF (LJCOULT) THEN
         WRITE(MYUNIT,'(I4,A,I4,A,F10.6)') NATOMS, ' Lennard-Jones with ', COULN, ' carrying positive charge ', COULQ
         WRITE(MYUNIT,'(A,F10.6)') 'Fraction of swap moves between charged and neutral particles: ', COULSWAP
         WRITE(MYUNIT,'(A,F10.6)') 'MC temperature for swap moves: ', COULTEMP
!
!  anisotropic particles:
!
!     DC430 >

      ELSE IF (LINRODT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' linear rod-like molecules '

      ELSE IF (CAPBINT) THEN

         WRITE(MYUNIT,'(I4,A,I4,A,I4,A)') NATOMS/2,' CAPSID molecules: ', NPS, ' pentamers and ', (NATOMS/2-NPS), ' hexamers'
         WRITE(MYUNIT,'(A,5F15.5)') ' epsilon_rep, RHO, radius, height1, height2 = ', 
     &                                      CAPEPS2, CAPRHO, CAPRAD, CAPHEIGHT1, CAPHEIGHT2
      ELSE IF (DBPT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polar dumbbell molecules '
         WRITE(MYUNIT,'(A,3F15.5)') 'DBEPSAA, DBEPSBB, DBEPSAB =', 1.0, DBEPSBB, SQRT(DBEPSBB)
         WRITE(MYUNIT,'(A,3F15.5)') 'DBSIGAA, DBSIGBB, DBSIGAB =', 1.D0, DBSIGBB, 0.5D0*(1.D0+DBSIGBB)   
         WRITE(MYUNIT,'(A,2F15.5)') 'DBPMU, EFIELD =', DBPMU, EFIELD   

      ELSE IF (DBPTDT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polar dumbbell molecules '
         WRITE(MYUNIT,'(A,4F15.5)') 'DBSIGAA, DBEPSAA, DBEPSBB, DBEPSAB =', 1.0, 1.0, DBEPSBB, DBEPSAB
         WRITE(MYUNIT,'(A,4F15.5)') 'DBSIGBB, DBSIGAB, DBPMU, EFIELD =', DBSIGBB, DBSIGAB, DBPMU, EFIELD   

      ELSE IF (DMBLMT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polar dumbbell molecules with Morse potential'
         WRITE(MYUNIT,'(A,3F15.5)') 'EPS11, EPS22, EPS12 =', EPS11, EPS22, EPS12
         WRITE(MYUNIT,'(A,3F15.5)') 'RHO11, RHO22, RHO12 =', MRHO11, MRHO22, MRHO12
         WRITE(MYUNIT,'(A,3F15.5)') 'REQ11, REQ22, REQ12 =', REQ11, REQ22, REQ12
         WRITE(MYUNIT,'(A,2F15.5)') 'DBPMU, EFIELD =', DBPMU, EFIELD

      ELSE IF (DMBLPYT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polar dumbbell molecules with  Yukawa potential'
         WRITE(MYUNIT,'(A,3F15.5)') 'YEPS, YKAPPA, DBSIGBB =', YEPS, YKAPPA, DBSIGBB
         WRITE(MYUNIT,'(A,2F15.5)') 'DBPMU, EFIELD =', DBPMU, EFIELD

      ELSE IF (LWOTPT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' lwotp molecules '

      ELSE IF (NCAPT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' CAPSID molecules '

      ELSE IF (NPAHT .OR. PAHAT .OR. PAHW99T) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAH molecules '

      ELSE IF (PAPT) THEN

         CALL DEFPAP()
         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAP building blocks with ', NRBSITES, ' patches and anitpatches'
         WRITE(MYUNIT,'(A,4F15.5)') 'ALPHA, S, COSDEL, EPSPA =', PAPALP, PAPS, PAPCD, PAPEPS

      ELSE IF (PAPBINT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAPBIN building blocks with ', NRBSITES, ' patches of two different kinds'
         WRITE(MYUNIT,'(A,4F15.5)') 'PAPEPS, PAPS, PAPANG1, PAPANG2, YKAPPA =', PAPEPS, PAPS, PAPANG1, PAPANG2, YKAPPA

      ELSE IF (PAPJANT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAPJANUS building blocks with ', NRBSITES, ' patches of two different kinds'
         WRITE(MYUNIT,'(A,4F15.5)') 'PAPEPS, PAPS, PAPANG1, YKAPPA =', PAPEPS, PAPS, PAPANG1, YKAPPA

      ELSE IF (NTIPT) THEN

         WRITE(MYUNIT,'(I4,A, I1)') NATOMS/2,' water molecules with TIPID = ', TIPID

!|gd351>
      ELSE IF (PATCHY) THEN

         WRITE(MYUNIT,'(I4,A, I1)') NATOMS/2,' patchy particles with sites = ', NRBSITES

      ELSE IF (ASAOOS) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS,' Asakura-Oosawa particles'

!<gd351|

      ELSE IF (GBT .OR. GBDT .OR. GBDPT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' rigid bodies '

      ELSE IF (GEMT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS,' gem atoms '

      ELSE IF (NPAHT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAH molecules with anistropic potential'
         WRITE(MYUNIT,'(A,I4)') 'PAHID =', PAHID

      ELSE IF (MSGBT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' rigid bodies '

      ELSE IF (MSTBINT) THEN

         WRITE(MYUNIT,'(I4,A)') NPS,' pentagonal polar molecules '
         WRITE(MYUNIT,'(I4,A)') NATOMS/2 - NPS,' hexagonal polar molecules '  
         WRITE(MYUNIT,'(A,4F15.5)') 'STOCKMU =', STOCKMU
         WRITE(MYUNIT,'(A,4F15.5)') 'EFIELD =', EFIELD

      ELSE IF (MSSTOCKT) THEN

         IF (NRBSITES == 3) THEN
            WRITE(MYUNIT,'(I4,A)') NATOMS/2,' triangular polar stockmayer molecules '
         ELSE IF (NRBSITES == 4) THEN
            WRITE(MYUNIT,'(I4,A)') NATOMS/2,' square polar stockmayer molecules '
         ELSE IF (NRBSITES == 5) THEN
            WRITE(MYUNIT,'(I4,A)') NATOMS/2,' pentagonal polar stockmayer molecules '
         ENDIF 
         WRITE(MYUNIT,'(A,6F15.5)') 'DPMU =', (DPMU(J1), J1=1, NRBSITES)
         WRITE(MYUNIT,'(A,4F15.5)') 'EFIELD =', EFIELD

      ELSE IF (MULTPAHAT) THEN

           WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PAHA molecules in a mixture '
           WRITE(MYUNIT,'(A,4I4)')  ' Nos. benzene, napthtalene, anthracene, pyrene, ', NCMP(1), NCMP(2), NCMP(3), NCMP(4)

      ELSE IF (SANDBOXT) THEN

         CALL SANDBOX_INPUT(MYUNIT)

      ELSE IF (SILANET) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' silane  molecules '

      ELSE IF (STOCKAAT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' stockmayer  particle '
         WRITE(MYUNIT,'(A,4F15.5)') 'STOCKMU =', STOCKMU
         WRITE(MYUNIT,'(A,4F15.5)') 'EFIELD =', EFIELD

      ELSE IF (MORSEDPT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' dipolar morse particles '
         WRITE(MYUNIT,'(A,4F15.5)') 'dipole strength =', STOCKMU
         WRITE(MYUNIT,'(A,4F15.5)') 'electric field strength =', EFIELD

      ELSE IF (TDHDT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' tetrahedra '
         WRITE(MYUNIT,'(A,2F15.5)') 'Morse rho, req =', RHO, MREQ

      ELSE IF (WATERDCT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polarisableDC water molecules '

      ELSE IF (WATERKZT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' polarisableKZ water molecules '

      ELSE IF (GAYBERNET) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' Gay-Berne atoms '

      ELSE IF (PYGPERIODICT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' ellipsoids '

      ELSE IF (PYT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' PY building blocks'

      ELSE IF (GAYBERNEDCT) THEN

         WRITE(MYUNIT,'(I4,A)') NATOMS/2,' Gay-Berne atoms '
         DO J2=1,NPAR
            WRITE(MYUNIT,*) COORDS(:,NPAR)
            CALL GBINIT(COORDS(:,NPAR))
         ENDDO
      ELSE IF (STICKYT) THEN
         WRITE(MYUNIT,'(I4,A,F15.5,A,F15.5)') NATOMS/2,' Sticky patch molecules'
      ELSE IF (NATBT) THEN
         WRITE(MYUNIT,'(I4,A,F15.5,A,F15.5)') NATOMS,' Na tight-binding atoms'
         IF (GUPTAT) WRITE(MYUNIT,'(A,I4,A,F12.5)') ' Minimisations will be guided by Gupta potential type ',GATOM,
     1            ' until the RMS force < ',GUIDECUT
      ELSE IF (DIFFRACTT) THEN
         WRITE(MYUNIT,'(A,I6,A)') 'Fit a diffraction pattern with ',3*NATOMS,' variables'
      ELSE IF (CHRMMT) THEN
         IF (CHARMMTYPE.EQ.1) WRITE(MYUNIT, '(I6,A)') NATOMS,' CHARMM atoms for CHARMM22 with ACE capgroup'
         IF (CHARMMTYPE.EQ.2) WRITE(MYUNIT, '(I6,A)') NATOMS,' CHARMM atoms for CHARMM19 with ACE capgroup'
         IF (CHARMMTYPE.EQ.3) WRITE(MYUNIT, '(I6,A)') NATOMS,' CHARMM atoms for CHARMM19 with no capgroup'
!hk286
         IF (RIGIDINIT) THEN
            WRITE(MYUNIT, '(I6,A)') NRIGIDBODY,' RIGID BODIES'
            WRITE(MYUNIT, '(I6,A)') (DEGFREEDOMS-6*NRIGIDBODY)/3, ' SINGLE ATOMS'
         ENDIF
!hk286
         WRITE(MYUNIT, '(A,2F12.4)') 'Maximum and minimum probabilities for twisting are ',CHPMAX,CHPMIN
!bs360         WRITE(MYUNIT, '(A,2F12.4)') 'Maximum and minimum number of twists are ',CHNMAX,CHNMIN
         IF (OMEGAT) WRITE(MYUNIT, '(A)') 'Peptide bonds will be twisted along with all other dihedrals'
      ELSEIF (THOMSONT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' Thomson ions'
         IF (ODDCHARGE.NE.1.0D0) WRITE(MYUNIT,'(A,F20.10)') 'one odd charge of magnitude ',ODDCHARGE
      ELSEIF (AMHT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' AMH Atoms'
      ELSEIF (QDT) THEN
         WRITE(*,'(I4,A)') 3*NATOMS,' quadratic degrees of freedom'
      ELSEIF (QD2T) THEN
         WRITE(*,'(I4,A)') 3*NATOMS,' test degrees of freedom'
      ELSEIF (MULLERBROWNT) THEN
         WRITE(MYUNIT,'(A)') 'Muller-Brown surface'
      ELSEIF (JMT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' atoms for JM potential (with cutoff)'
      ELSEIF (DF1T) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' atoms for DF1 potential'
      ELSEIF (TTM3T) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS/3,' TTM3-F water molecules'
      ELSEIF (DFTBCT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' DFTB carbon atoms'
      ELSEIF (GLJT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' generalised LJ atoms'
      ELSEIF (CUDAT) THEN
          SELECT CASE(CUDAPOT)
              CASE('L')
                  WRITE(MYUNIT,'(A,I4,A)') 'LBFGS minimisation on the GPU with ', NATOMS,' LJ atoms'
              CASE('A')
                  WRITE(MYUNIT,'(A,I4,A)') 'LBFGS minimisation on the GPU with ', NATOMS,' AMBER12 atoms'
              CASE DEFAULT
                  WRITE(MYUNIT,'(A,I4,A)') 'LBFGS minimisation on the GPU with ', NATOMS,' atoms. Potential not recognised. '
          END SELECT
      ELSEIF (QCIPOTT.OR.INTCONSTRAINTT) THEN
         WRITE(MYUNIT,'(I4,A)') NATOMS,' QCI atoms'
      ELSE
         WRITE(MYUNIT,'(I4,A)') NATOMS,' LJ atoms'
      ENDIF
      IF (PYGPERIODICT.OR.PYBINARYT) CALL INITIALISEPYGPERIODIC
      IF (LJCAPSIDT) CALL INITIALISELJCAPSIDMODEL
      IF (PYT) call py_input
      IF (INTMINT)  WRITE(MYUNIT,'(A)') 'Internal coordinate transformation will be used'
      IF (STAR) THEN
         WRITE(MYUNIT,'(A)') 'Excited state'
      ELSE IF (PLUS) THEN
          WRITE(MYUNIT,'(A)') 'Single positive charge'
      ELSE IF (TWOPLUS) THEN
         WRITE(MYUNIT,'(A)') 'Double positive charge'
      ENDIF
      IF (SPARSET) THEN
         WRITE(MYUNIT, '(A,G20.10)') 'io1> Using sparse matrix approach for Hessian, cutoff=',ZERO_THRESH
      ENDIF
      IF (SYMMETRIZE) THEN
         WRITE(MYUNIT, '(A,I6,A)') 'Searching for approximate symmetry elements every ',NSYMINTERVAL,' steps'
         WRITE(MYUNIT, '(A,5F15.5)') 'Distance tolerances: ',SYMTOL1,SYMTOL2,SYMTOL3,SYMTOL4,SYMTOL5
         WRITE(MYUNIT, '(A,F15.5)') 'Threshold for distinguishing transformation matrices: ',MATDIFF
         WRITE(MYUNIT,'(A,F15.5)') 'Exponential factor in core-weighted centre of mass calculation: ',DISTFAC
         WRITE(MYUNIT, '(A,I5)') 'Maximum number of quenches for floater/hole permutations=',NSYMQMAX
         IF (MOVESHELLT) WRITE(MYUNIT,'(A,I8,A,F12.5)') 'Shell moves allowed in blocks of ',SHELLMOVEMAX,' with probability ',
     &                        SHELLPROB
      ENDIF
      IF (DEBUG.OR.CHECKMARKOVT) WRITE(MYUNIT,'(A,I6,A)') 'io1> checking the energy of the saved coordinates in the chain'
      IF (FREEZE) THEN
         WRITE(MYUNIT,'(A,I6,A)') 'io1> ', NFREEZE,' atoms will be frozen:'
         DO J1=1,NATOMS
            IF (FROZEN(J1)) WRITE(MYUNIT,'(I6)') J1
         ENDDO
      ENDIF
      IF (HARMONICF) THEN
         WRITE(MYUNIT,'(A,F12.4)') 'io1> harmonically constrained atoms: strength = ', HARMONICSTR
         DO J1=1,NATOMS
            IF (HARMONICFLIST(J1)) WRITE(MYUNIT,'(I6)') J1
         ENDDO
      ENDIF
      IF (COOP) THEN
         WRITE(MYUNIT,'(A,I5,A,F12.4,A)') 'Cooperative moves will be taken for the ',NCOOP,
     &                             ' nearest neighbours within a cutoff of ',COOPCUT,' for a randomly chosen atom'
      ENDIF
      IF (COMPRESST) THEN
         WRITE(MYUNIT,'(A,G15.8)') 'Initial quenching with a radial compression, force constant=',K_COMP
      ENDIF
      IF (FIXCOM .AND. (.NOT. CHRMMT)) THEN
          INQUIRE(FILE='masses',EXIST=EXISTS)
          IF (EXISTS) THEN
              OPEN(UNIT=1978,FILE="masses")
              DO J2=1,NATOMS
                 READ(1978,'(F12.5)') ATMASS(J2)
              ENDDO
              CLOSE(1978)
          ELSE
             WRITE(MYUNIT,'(A)') 'WARNING, FIXCOM is specified, but "masses" file is not present: setting all masses to unity.'
             ATMASS(1:NATOMS) = 1.0D0
          ENDIF
      ENDIF
      IF (BSPT) THEN
          WRITE(MYUNIT,'(A)') '------------------------------------------------------------------------------------'
          WRITE(MYUNIT,'(A)') 'S T A R T I N G   P T   B A S I N   S A M P L I N G   R U N:'
          WRITE(MYUNIT, '(A)') 'Parallel tempering basin-sampling run'
          WRITE(MYUNIT,'(A,F15.8)') 'Lowest  potential energy', PTEMIN
          WRITE(MYUNIT,'(A,F15.8)') 'Highest potential energy', PTEMAX
          WRITE(MYUNIT,'(A,G15.8)') 'Lowest  temperature', PTTMIN
          WRITE(MYUNIT,'(A,G15.8)') 'Highest temperature', PTTMAX
          WRITE(MYUNIT,'(A,F15.8)') 'Exchange probability', EXCHPROB
          WRITE(MYUNIT,'(A,F15.1)') 'Number of equilibration steps=',NEQUIL
          WRITE(MYUNIT,'(A,I15)')   'Interval between quenches=',QUENCHFRQ
          WRITE(MYUNIT,'(A)') '-------------------------------------------------------------------------------------'
      ENDIF
      IF (BSWL) THEN
          WRITE(MYUNIT,'(A)') '------------------------------------------------------------------------------------'
          WRITE(MYUNIT,'(A)') 'S T A R T I N G   W L  B A S I N   S A M P L I N G   R U N:'
          WRITE(MYUNIT, '(A)') 'A zero-temperature histogram-based MC sampling of the energy density of local minima'
          WRITE(MYUNIT,'(A,F15.8)') 'Energy of the lowest bin ', HISTMIN
          WRITE(MYUNIT,'(A,F15.8)') 'Bin width ', (HISTMAX-HISTMIN)/HBINS
          WRITE(MYUNIT,'(A,I10)') 'Number of bins', HBINS
          WRITE(MYUNIT,'(A,F15.8,A)') 'Flatness criterion',  HPERCENT, ' percent'
          WRITE(MYUNIT,'(A,F15.8)') 'Starting WL modification factor', HISTFAC
          WRITE(MYUNIT,'(A,G20.10)') 'Requested number of Wang-Landau iterations', TargetWL
          IF (EQUIL.GT.0)  WRITE(MYUNIT,'(A,I15,A)') 
     1                     'The accumulation of the density of states will be preceded with ',  EQUIL , ' steps'
          WRITE(MYUNIT,'(A,2G20.10)') 'Temperature range for calculation ensemble averages:', MinimalTemperature, MaximalTemperature
          IF (BINSTRUCTURES) WRITE(MYUNIT,'(A,I15,A)') 'Coordinate of every',SaveNth, 'will be recorded' 

          WRITE(MYUNIT,'(A)') '-------------------------------------------------------------------------------------'

!
!  There is printing in keyword.f if hist.minlist or hist.old are used.
!
      ENDIF
      IF (RADIUS.EQ.0.0D0) THEN
!        RADIUS=1.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0) ! this container is too small for angular moves
         RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
         IF (MORSET) THEN
         ELSE IF (NATBT) THEN
            RADIUS=RADIUS*6.2D0
         ELSE IF (MSORIGT.OR.FRAUSIT.OR.MSTRANST) THEN
            RADIUS=RADIUS*2.5D0
            IF (.NOT.ANGST) RADIUS=RADIUS*1.89D0
         ELSE IF (NEON) THEN
            RADIUS=RADIUS*2.6D0
         ELSE IF (ARGON) THEN
            RADIUS=RADIUS*7.0D0
         ELSE IF (SCT) THEN
            RADIUS=SIG*RADIUS/DSQRT(2.0D0)
         ELSE IF (DFTBT) THEN
             RADIUS=RADIUS*5.0D0
         ELSE IF (P46) THEN
            RADIUS=RADIUS*3.0D0
         ELSE IF (SW) THEN
            RADIUS=RADIUS*2.4D0
         ELSE IF (QUIPT) THEN
            RADIUS=RADIUS*QUIPEQDIST
         ELSE IF (RGCL2) THEN
            RADIUS=RADIUS*3.0D0
         ELSE IF (ARNO) THEN
            RADIUS=RADIUS*3.8D0
         ELSE IF (TOSI.OR.WELCH) THEN
!
!  The extra factor is needed because the natural lattice isnt close-packed.
!
            RADIUS=RADIUS*10.0D0
         ELSE IF (AMBER) THEN
         ELSE IF (CHRMMT) THEN
         ELSE IF (CPMD) THEN
         ELSE IF (DZTEST) THEN
!            RADIUS=1.3D0*RADIUS*2.0D0**(1.0D0/6.0D0)
            RADIUS=2.0D0*RADIUS*2.0D0**(1.0D0/6.0D0)
         ELSE IF (ZETT1.OR.ZETT2) THEN
!           RADIUS=4.0D0*RADIUS*2.0D0**(1.0D0/6.0D0)
            RADIUS=MAX(NATOMS*4.0D0/(2*19.0D0),2.0D0)
         ELSE IF (PACHECO) THEN
            RADIUS=RADIUS*10.0185D0
         ELSE IF (PBGLUET) THEN
            RADIUS=RADIUS*3.4d0
         ELSE IF (ACKLANDT) THEN
            RADIUS=RADIUS*3.0D0
         ELSE IF (EAMALT) THEN
            RADIUS=RADIUS*2.8d0
         ELSE IF (ALGLUET) THEN
            RADIUS=RADIUS*2.85d0
         ELSE IF (MGGLUET) THEN
            RADIUS=RADIUS*3.2d0
         ELSE IF (TIP) THEN
            RADIUS=RADIUS*3.0D0
         ELSE IF (CAPSID) THEN
            RADIUS=RADIUS*RAD*1.5D0
         ELSE IF (LB2T) THEN
            RADIUS=RADIUS*1.3D0
         ELSE IF (CHRMMT) THEN
            RADIUS=1.0D6
         ELSE IF (THOMSONT) THEN
            RADIUS=2.0D0
!           DO J1=1,NPAR
!              STEP(J1)=STEP(J1)*(0.677441D0-0.0037582*NATOMS+9.40318D-6*NATOMS**2-6.21931D-9*NATOMS**3)
!           ENDDO
!            WRITE(MYUNIT, '(A,G15.5,A,G15.5)') 'Maximum step size scaled by estimated nearest neighbour distance of ',
!    &                    0.677441D0-0.0037582*NATOMS+9.40318D-6*NATOMS**2-6.21931D-9*NATOMS**3,' to give ',STEP(1)
         ELSEIF (MULLERBROWNT) THEN 
            RADIUS=100.0D0
         ELSE 
            RADIUS=RADIUS*2.0D0**(1.0D0/6.0D0)
         ENDIF
      ENDIF
      IF ((.NOT.PERIODIC).AND.(.NOT.AMBER).AND.(.NOT.BLNT).AND.(.NOT.MULLERBROWNT).AND.(.NOT.MODEL1T).AND.(.NOT.PERCOLATET) 
     &                    .AND.(.NOT.QCIPOTT).AND.(.NOT.INTCONSTRAINTT)) 
     1                    WRITE(MYUNIT,'(A,F20.10)') 'Container radius=',RADIUS
      RADIUS=RADIUS**2
      IF (PERCOLATET) WRITE(MYUNIT,'(A,F20.10)') 'Checking for percolated structure, cutoff=',PERCCUT
      PERCCUT=PERCCUT**2
      IF (NPAR.GT.1) THEN
         WRITE(MYUNIT,'(I2,A)') NPAR,' parallel runs'
         IF (TABOOT) WRITE(MYUNIT,'(A,I4,A)') 'Taboo lists contain the lowest ',NTAB,' minima'
      ELSE IF (TABOOT) THEN
         WRITE(MYUNIT,'(A,I4,A)') 'Taboo list contains the lowest ',NTAB,' minima'
      ENDIF
      IF (SUPERSTEP) THEN
         WRITE(MYUNIT,'(A,I3,A)') 'Taking supersteps over ',NSUPER,' minima'
         WRITE(MYUNIT,'(A,F15.7)') 'Fixed temperature for supersteps=',TEMPS
         WRITE(MYUNIT,'(A,F12.5,A,I4,A,F12.4)') 
     1   'Initial superstep factor ',SUPSTEP,' will be adjusted every ',NSACCEPT,' steps to give acceptance ratio ',SACCRAT
      ENDIF

      DO J1=1,NPAR
         IF (JUMPMOVE(J1)) WRITE(MYUNIT,'(A,I1,A,I1,A,I4,A)') 
     1                  'Jump moves will be attempted from run ',J1,' to run ',JUMPTO(J1),' every ',JUMPINT(J1),' steps'
      ENDDO
      IF (NEWJUMP) WRITE(MYUNIT,'(A,F12.3)') 
     1   'Jumping based only on current energies (parallel tempering) attempt probability=',PNEWJUMP
      WRITE(MYUNIT,'(A,G15.10)') 'Sloppy quench tolerance for RMS gradient ',BQMAX
      IF ((.NOT.BSPT).AND.(.NOT.PTMC)) THEN
         DO JP=1,NPAR
            IF (RIGID.AND.(BLOCK(NPAR).GT.0)) WRITE(MYUNIT,'(A,I6,A)') 
     1     'Rigid body translations and orientational displacements will be made separately in blocks of ',
     2      BLOCK(NPAR),' steps'
            IF (FIXBOTH(JP)) THEN
               WRITE(MYUNIT,'(A,I3,A,F12.4,A,2F12.4,A)') 
     1                 'In run ',JP,' temperature=',TEMP(JP),' step size and angular threshold=',
     1                  STEP(JP),ASTEP(JP),' all fixed'
               IF (RIGID) WRITE(MYUNIT,'(A,F12.4)') 'Orientational step size for rigid bodies also fixed at ',
     1                                          OSTEP(JP)
               IF (FIXD) WRITE(MYUNIT,'(A,I3,A,I4)') 'In run ',JP,' number of hard sphere collision moves fixed at ',
     1                                          NHSMOVE
            ELSE IF (FIXSTEP(JP)) THEN
               WRITE(MYUNIT,'(A,I3,A,2F12.4)') 'In run ',JP,' step size and angular threshold fixed at ',
     1                                    STEP(JP),ASTEP(JP)
               IF (RIGID) WRITE(MYUNIT,'(A,F12.4)') 'Orientational step size for rigid bodies also fixed at ',
     1                                          OSTEP(JP)
               IF (FIXD) WRITE(MYUNIT,'(A,I3,A,I4)') 'In run ',JP,' number of hard sphere collision moves fixed at ',
     1                                          NHSMOVE
               IF (.NOT.FIXTEMP(JP)) THEN
                  WRITE(MYUNIT,'(A,F12.4,A,F12.4)') 
     1                    'Temperature will be adjusted for acceptance ratio ',ACCRAT(JP),' initial value=',TEMP(JP)
               ELSE
                  WRITE(MYUNIT,'(A,I1,A,G12.4)') 'In run ',JP,' temperature will be fixed - see below for value'
               ENDIF
            ELSE IF (STEPOUT) THEN
               WRITE(MYUNIT,'(A,I3,A,2F12.4,A,2F12.4)') 
     1   'In run ',JP,' step size and angular threshold will be adjusted to escape from basins. Initial values=',
     1                  STEP(JP),ASTEP(JP)
               IF (RIGID) WRITE(MYUNIT,'(A,F12.4)') 'Orientational step size for rigid bodies initial value ',
     1                                          OSTEP(JP)
               IF (FIXD) WRITE(MYUNIT,'(A,I3,A,I4)') 
     1            'In run ',JP,' number of hard sphere collision moves will be adjusted. Initial value=',NHSMOVE
               IF (.NOT.FIXTEMP(JP)) THEN
                  WRITE(MYUNIT,'(A,F12.4,A,F12.4)') 
     1                    'Temperature will be adjusted for acceptance ratio ',ACCRAT(JP),' initial value=',TEMP(JP)
               ELSE
                  WRITE(MYUNIT,'(A,I1,A,G12.4)') 'In run ',JP,' temperature fixed at ',TEMP(JP)
               ENDIF
            ELSE 
               WRITE(MYUNIT,'(A,I3,A,G12.4)') 'In run ',JP,' temperature will be fixed - see below for value'
               WRITE(MYUNIT,'(A,F12.4,A,2F12.4)') 'Step size and angular threshold will be adjusted for acceptance ratio ',
     1                ACCRAT(JP),' initial values=',STEP(JP),ASTEP(JP)
               IF (RIGID) WRITE(MYUNIT,'(A,F12.4)') 'Orientational step size for rigid bodies initial value ',
     1                                          OSTEP(JP)
               IF (FIXD) WRITE(MYUNIT,'(A,I3,A,I4)') 
     1            'In run ',JP,' number of hard sphere collision moves will be adjusted. Initial value=',NHSMOVE
            ENDIF
         ENDDO 
      ENDIF
      IF (EFAC.NE.0.0D0) WRITE(MYUNIT,'(A,F12.4)') 'Exponential factor for proposed steps=',EFAC
      IF (NORESET.OR.BSPT) THEN
         WRITE(MYUNIT,'(A)') 'Configuration will not be reset to quench geometry'
         IF (CENT) THEN
            WRITE(MYUNIT,'(A)') 'WARNING CENTRE can lead to atoms leaving '
            WRITE(MYUNIT,'(A)') 'the container after takestep when the centre of mass is moved.'
!           STOP
         ENDIF
      ELSE
         WRITE(MYUNIT,'(A)') 'Configuration will be reset to quench geometry'
      ENDIF
      IF (CENT .AND. FIXCOM) THEN
          WRITE(MYUNIT,'(A)') 'WARNING: keywords CENTRE (fixing centre of coordinates) and FIXCOM (fixing centre of mass) 
     1                    are incompatible'
          STOP
      ENDIF
      IF (WELCH.OR.CPMD.OR.PACHECO) WRITE(MYUNIT,'(A,F15.5)') 'Guiding function used for RMS>',GUIDECUT
      IF (DFTBCT) WRITE(MYUNIT,'(2(A,F15.5))') 'Guiding function LJAT with Z*=3.0 will be used for RMS>',GUIDECUT,
     &                                         ' rescale factor: ',LJATTOC
      IF (TSALLIST) THEN    
         WRITE(MYUNIT,'(A,F15.5)') 'Sampling with Tsallis statistics, q=',QTSALLIS
         IF (.NOT.SAVEQ) SAVEQ=.TRUE.
         IF (NSAVE.LT.1) NSAVE=1
      ELSE
         WRITE(MYUNIT,'(A)') 'Sampling using Boltzmann weights'
      ENDIF
      IF (PERIODIC) WRITE(MYUNIT,'(A,3F15.7)') 'Periodic boundary conditions, box lengths: ',BOXLX,BOXLY,BOXLZ
      IF (CUTT) WRITE(MYUNIT,'(A,F15.7)') 'Cutoff=',CUTOFF
      IF (BFGS) THEN
         WRITE (MYUNIT,'(A)') 'BFGS minimization'
      ELSE IF (RKMIN) THEN
         WRITE (MYUNIT,'(A)') 'Runge-Kutta steepest-descent minimisation'
      ELSE IF (BSMIN) THEN
         WRITE (MYUNIT,'(A)') 'Bulirsch-Stoer steepest-descent minimisation'
      ELSE IF (LBFGST) THEN
         WRITE (MYUNIT,'(A)') 'Nocedal LBFGS minimisation'
         WRITE(MYUNIT,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
         WRITE(MYUNIT,'(A,F20.10)') 'Maximum step size=',MAXBFGS
         WRITE(MYUNIT,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
      ELSE
         WRITE (MYUNIT,'(A)') 'Conjugate gradient minimization'
      ENDIF
      WRITE(MYUNIT,'(A,G15.10)') 'Final quench tolerance for RMS gradient ',CQMAX
      WRITE(MYUNIT,'(A,F15.10)') 'Energy difference criterion for minima=',ECONV
      WRITE(MYUNIT,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2
      IF (.NOT.BSPT) THEN
         DO J1=1,NRUNS
            WRITE(MYUNIT,120) J1, MCSTEPS(J1), TFAC(J1)
120         FORMAT('Run ',I3,': ',I9,' steps with temperature scaled by ',E15.8)
         ENDDO
      ENDIF
      IF (TUNNELT) THEN
         WRITE(MYUNIT,'(A,F15.5,A)') ' Transforming potential to 1-exp(-(E-E0)*',GAMMA,')'
      ENDIF
      IF (DEBUG) THEN
         WRITE(MYUNIT,160) 
160      FORMAT('Debug printing is on')
      ENDIF

! beh35 > adding stuff for mlowest

      IF (MLOWEST) THEN
         WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(MYUNIT,'(F20.10)',ADVANCE='NO') (OBJ(J1)%SOUGHT,J1=1,MTARGETS)
         WRITE(MYUNIT,'(A)') ' '
      ENDIF


! beh35 > end added stuff


       WRITE(MYUNIT, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE
      IF (OHT)   WRITE(MYUNIT,'(A26,F12.4)') 'Octahedral field strength=',FOH
      IF (IHT)   WRITE(MYUNIT,'(A27,F12.4)') 'Icosahedral field strength=',FIH
      IF (TDT)   WRITE(MYUNIT,'(A27,F12.4)') 'Tetrahedral field strength=',FTD
      IF (D5HT)  WRITE(MYUNIT,'(A19,F12.4)') 'D5h field strength=',FD5H
      IF (TARGET) THEN
         WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(MYUNIT,'(F20.10)',ADVANCE='NO') (TARGETS(J1),J1=1,NTARGETS)
         WRITE(MYUNIT,'(A)') ' '
      ENDIF
      IF (SQUEEZET) THEN
          WRITE(MYUNIT, '(A)') 'squeeze option currently commented'
         STOP
         OPEN (UNIT=7,FILE='vectors',STATUS='OLD')
         VECMN=1.0D6
         DO J1=1,NVEC
            READ(7,*,END=180) VEC(3*J1-2),VEC(3*J1-1),VEC(3*J1)
            DUMMY=VEC(3*J1-2)**2+VEC(3*J1-1)**2+VEC(3*J1)**2
            IF (DUMMY.LT.VECMN) VECMN=DUMMY
         ENDDO
180      VECMN=SQUEEZER/DSQRT(VECMN)
         DO J1=1,NVEC
            VEC(3*J1-2)=VEC(3*J1-2)*VECMN
            VEC(3*J1-1)=VEC(3*J1-1)*VECMN
            VEC(3*J1)=VEC(3*J1)*VECMN
         ENDDO
         WRITE(MYUNIT,190) NVEC
190      FORMAT('System will be squeezed into the shape defined by the '
     1          ,I4,' vectors:')
         WRITE(MYUNIT,200) (VEC(J1),J1=1,3*NVEC)
200      FORMAT(3F20.10)
         WRITE(MYUNIT,210) SQUEEZER, SQUEEZED
210      FORMAT('Initial smallest radius=',F15.5,' with shrink factor='
     1          ,F15.5)
      ENDIF
      IF (RENORM) THEN
         XMOVERENORM=MIN(XMOVERENORM,NATOMS*0.9D0)
         IF (XMOVERENORM.GT.3.0D0) THEN
            WRITE(MYUNIT,'(A,F12.1,A)')  'Large steps of ',XMOVERENORM,' hard sphere type moves will be used:'
         ELSE
            WRITE(MYUNIT,'(A,F15.5,A)')  'Large steps using maximum displacement ',XMOVERENORM,' will be used:'
         ENDIF
         WRITE(MYUNIT,'(A,I6)')    '  Initial interval for large steps is ',NRENORM
         WRITE(MYUNIT,'(A,F12.5)') '  Temperature used in Metropolis test for acceptance of large steps is ',TRENORM
      ENDIF
      IF (RESTORET) THEN
         WRITE(MYUNIT,'(A,A)') 'Restoring GMIN run from file ',TRIM(ADJUSTL(DUMPFILE))
      ENDIF
      IF (NEWRESTART) THEN
         IF (.NOT.AVOIDRESEEDT) THEN
            WRITE(MYUNIT,'(A,F12.5,A,I6,A)') 'Steps will be rejected and taboo list populated if the energy decrease < ',ECONV,
     &                                ' within ',NRELAX,' steps'
         ELSE
            WRITE(MYUNIT,'(A,F12.5,A,I6,A)') 'Runs will be reseeded if the energy does not decrease by at least ',ECONV,
     &                                ' within ',NRELAX,' steps'
         ENDIF
         IF (NHSRESTART.GT.0)  WRITE(MYUNIT, '(I6,A)') NHSRESTART,' hard sphere-type moves will be used to reseed'
      ENDIF
      IF (AVOID) THEN
         IF (AVOIDRESEEDT) THEN
            WRITE(MYUNIT,'(A,F10.2,A,I6,A)') 'Runs will be reseeded if the current minimum comes within ',
     1                 AVOIDDIST,' of up to ',MAXSAVE,' previous minima'
         ELSE
            WRITE(MYUNIT,'(A,F10.2,A,I6,A)') 'Steps will be rejected if the current minimum comes within ',
     1                 AVOIDDIST,' of up to ',MAXSAVE,' previous minima'
         ENDIF
      ENDIF
!
!  Look for the file that contains interrupted screen saver restart information.
!  Current minimum in the Markov chain. COORDS
!  Number of steps done. NQTOT/NPAR should be close enough!
!  The current lowest minima. QMIN has the energies, QMINP has the points.
!  The current values of the temperature, acceptance ratio and step length,
!  TEMP(JP), ACCRAT(JP), STEP(JP), ASTEP(JP) and OSTEP(JP)
!  which can get changed dynamically.
!
      INQUIRE(FILE='ssdump',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=88,FILE='ssdump',STATUS='UNKNOWN')
          WRITE(MYUNIT,'(A)') 'reading dump information from file ssdump'
         READ(88,'(3G20.10)') ((COORDS(J1,J2),J1=1,3*NATOMS),J2=1,NPAR)
         READ(88,'(2I6)') NQTOT, NPCALL
         MCSTEPS(1)=MAX(MCSTEPS(1)-NQTOT*NPAR,1)
         NQTOT=NQTOT*NPAR
         READ(88,'(G20.10)') (QMIN(J1),J1=1,NSAVE)
         READ(88,'(3G20.10)') ((QMINP(J2,J1),J1=1,3*NATOMS),J2=1,NSAVE)
         READ(88,'(G20.10)') (TEMP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (ACCRAT(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (STEP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (ASTEP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (OSTEP(J1),J1=1,NPAR)
         CLOSE(88)
         YESNO=.FALSE.
      ENDIF


      RETURN
      END
