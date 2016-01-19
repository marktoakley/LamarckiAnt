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
      SUBROUTINE GSAVEIT(EREAL,P,NP)
      USE COMMONS
      use QMODULE
      IMPLICIT NONE


      INTEGER J1, J2, J3, NP, NQTOT, CSMIT 
      DOUBLE PRECISION EREAL,P(3*NATOMS), AVVAL, CSMRMS

      COMMON /TOT/ NQTOT
     
!
! Save the lowest minimum for each size with GCBH.
!
      IF (GCBHT) THEN
!        IF (EREAL.LT.QENERGIES(NATOMS)) THEN
         IF (FEBH_POT_ENE+NATOMS*GCMU.LT.QPE(NATOMS)) THEN
            QENERGIES(NATOMS)=EREAL
            QPE(NATOMS)=FEBH_POT_ENE+NATOMS*GCMU
            DO J2=1,3*NATOMS
               QCOORDINATES(NATOMS,J2)=P(J2)
            ENDDO
         ENDIF
      ENDIF
C
C  Save the lowest NSAVE distinguishable configurations.
C
C     WRITE(*,'(A,12E15.7)') 'EREAL,ECONV,QMIN',EREAL,ECONV,(QMIN(J1),J1=1,NSAVE)
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
C
C  These are probably the same - but just to make sure we save the 
C  lowest.
C
            IF (EREAL.LT.QMIN(J1)) THEN
               QMINNATOMS(J1)=NATOMS
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               QMINNATOMS(J2)=QMINNATOMS(J2-1)
               FF(J2)=FF(J2-1)
               NPCALL_QMIN(J2)=NPCALL_QMIN(J2-1)
               DO J3=1,3*QMINNATOMS(J2-1)
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            QMINNATOMS(J1)=NATOMS
            FF(J1)=NQ(NP)
            NPCALL_QMIN(J1)=NPCALL
            DO J2=1,3*QMINNATOMS(J1)
               QMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE


! beh35 > adding stuff to accommodate for mlowest

! here is how it works normally when you specifiy targets without mlowest modification
      
      IF(TARGET) THEN ! beh35
      DO J1=1,NTARGETS
         IF (EREAL-TARGETS(J1).LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(MYUNIT,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(MYUNIT,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',TARGETS(J1)
!            ELSE
!               WRITE(MYUNIT,'(A,I0.2,A,I2,2(A,I15),A)') '[',NP,']Target hit in parallel run ',NP,
!     1                    ' after ',NQTOT,' quenches ',NPCALL,' function calls.'
!               WRITE(MYUNIT,'(A,I0.2,2(A,F20.10))') '[',NP,']Energy=',EREAL,' target=',TARGETS(J1)
            ENDIF
            HIT=.TRUE.
         ENDIF
      ENDDO

! and here is what happens with mlowest
      
      ELSE IF(MLOWEST) THEN
      DO J1=1,MTARGETS
         IF (EREAL-OBJ(J1)%SOUGHT.LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(MYUNIT,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(MYUNIT,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',OBJ(J1)%SOUGHT

               !* maybe go back and look at this later to deal with parallel processes
            ENDIF
            OBJ(J1)%FOUND=.TRUE.
         ENDIF
      ENDDO

! now check and see if we've hit all of them

      CALL decide_to_quit()

      ENDIF

      RETURN
      END
C
C  csw34> new subroutine to save the lowest NSAVEINTE configurations indexed by interaction enthalpy between protein and ligand
C
      SUBROUTINE A9INTESAVEIT(INTEREAL,P,NP)
      USE commons
      use qmodule
      IMPLICIT NONE


      INTEGER J1, J2, J3, NP, NQTOT
      DOUBLE PRECISION INTEREAL,P(3*NATOMS)

      COMMON /TOT/ NQTOT
C
C  Save the lowest NSAVEINTE distinguishable configurations.
C
C     WRITE(*,'(A,12E15.7)') 'INTEREAL,ECONV,QMIN',INTEREAL,ECONV,(INTEQMIN(J1),J1=1,NSAVEINTE)
      DO J1=1,NSAVEINTE
         IF (DABS(INTEREAL-INTEQMIN(J1)).LT.ECONV) THEN
C
C  These are probably the same - but just to make sure we save the 
C  lowest.
C
            IF (INTEREAL.LT.INTEQMIN(J1)) THEN
               INTEQMIN(J1)=INTEREAL
               DO J2=1,3*NATOMS
                  INTEQMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (INTEREAL.LT.INTEQMIN(J1)) THEN

            J2=NSAVEINTE
20          CONTINUE
      
            IF (NSAVEINTE.GT.1) THEN
               INTEQMIN(J2)=INTEQMIN(J2-1)
               INTEFF(J2)=INTEFF(J2-1)
               DO J3=1,3*NATOMS
                  INTEQMINP(J2,J3)=INTEQMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            INTEQMIN(J1)=INTEREAL
            INTEFF(J1)=NQ(NP)
            DO J2=1,3*NATOMS
               INTEQMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      RETURN

      END SUBROUTINE A9INTESAVEIT
C
C ds656> Routine for MultiComponent systems with particle labels
C     
      SUBROUTINE GSAVEIT_MC(EREAL,P,LOT,NP)
      USE commons
      use qmodule
      IMPLICIT NONE
      
      INTEGER LOT(NATOMS)
      INTEGER J1, J2, J3, NP, NQTOT, CSMIT
      DOUBLE PRECISION EREAL,P(3*NATOMS), AVVAL, CSMRMS

      COMMON /TOT/ NQTOT
C
C  Save the lowest NSAVE distinguishable configurations.
C
C     WRITE(*,'(A,12E15.7)') 'EREAL,ECONV,QMIN',EREAL,ECONV,(QMIN(J1),J1=1,NSAVE)

C ds656> Testing
C      WRITE(MYUNIT,'(A)') 'saveit_mc> List of labels:'
c      J2=0
c      DO J1=1,NATOMS
c         WRITE(MYUNIT,'(I3)',ADVANCE='NO') LOT(J1)
c         IF(LOT(J1) == 1) J2=J2+1
c      ENDDO
c      WRITE(MYUNIT,*)
c      IF(J2 /= NSPECIES(1)) THEN
c         WRITE(MYUNIT,'(A,I4)') 'saveit_mc> Bad count for type-1:',J2
c      ENDIF
C <ds656 End of testing.
      
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
C
C  These are probably the same - but just to make sure we save the 
C  lowest.
C
            IF (EREAL.LT.QMIN(J1)) THEN
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
               DO J2=1,NATOMS
                  QMINT(J1,J2)=LOT(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               FF(J2)=FF(J2-1)
               NPCALL_QMIN(J2)=NPCALL_QMIN(J2-1)
               DO J3=1,3*NATOMS
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO
               DO J3=1,NATOMS
                  QMINT(J2,J3)=QMINT(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            FF(J1)=NQ(NP)
            NPCALL_QMIN(J1)=NPCALL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=P(J2)
            ENDDO
            DO J2=1,NATOMS
               QMINT(J1,J2)=LOT(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      IF(TARGET) THEN ! beh35

! normal target procedure 

      DO J1=1,NTARGETS
         IF (EREAL-TARGETS(J1).LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(MYUNIT,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(MYUNIT,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',TARGETS(J1)
!            ELSE
!               WRITE(MYUNIT,'(A,I0.2,A,I2,2(A,I15),A)') '[',NP,']Target hit in parallel run ',NP,
!     1                    ' after ',NQTOT,' quenches ',NPCALL,' function calls.'
!               WRITE(MYUNIT,'(A,I0.2,2(A,F20.10))') '[',NP,']Energy=',EREAL,' target=',TARGETS(J1)
            ENDIF
            HIT=.TRUE.
         ENDIF
      ENDDO

! and here is what happens with mlowest

      ELSE IF(MLOWEST) THEN
      DO J1=1,MTARGETS
         IF (EREAL-OBJ(J1)%SOUGHT.LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(MYUNIT,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(MYUNIT,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',OBJ(J1)%SOUGHT

               !* maybe go back and look at this later to deal with parallel processes
            ENDIF
            OBJ(J1)%FOUND=.TRUE.
         ENDIF
      ENDDO

! now check and see if we've hit all of them

      CALL decide_to_quit()

      ENDIF


      RETURN
      END

