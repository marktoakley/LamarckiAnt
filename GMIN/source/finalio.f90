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
SUBROUTINE FINALIO
    USE COMMONS
    USE GENRIGID, ONLY : RIGIDINIT
    USE MODAMBER
    USE MODAMBER9, ONLY : COORDS1,IH,M04,AMBFINALIO_NODE
    USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_FINISH, AMBER12_WRITE_RESTART, AMBER12_WRITE_PDB, &
                                      AMBER12_WRITE_XYZ
    USE QMODULE
    USE MODCHARMM
    USE AMHGLOBALS, ONLY:NMRES,IRES
    USE BGUPMOD
    USE PERMU

    IMPLICIT NONE

    !   MCP
    INTEGER III, I3,  GLY_COUNT, ID, NUMCRD, NUMPRO, NCPHST, GETUNIT, AMHUNIT1, AMHUNIT2, LUNIT, NSYMOPS
    INTEGER J1, J2, J3, J4, J5, MYUNIT2, I1, NDUMMY, MYUNIT3, NC, NRBS1, NRBS2, MYUNIT4
    DOUBLE PRECISION RBCOORDS(NRBSITES*3), DCOORDS(3*NATOMS), EDUMMY
    DOUBLE PRECISION P3(3,3), P(3), DU(3), RMI(3,3), DRMI(3,3), PI, PHI, THT, CHI
    !DOUBLE PRECISION, ALLOCATABLE :: XCOORDS(:), YCOORDS(:)
    CHARACTER(LEN=13) J1CHAR,J1CHAR2                  !for gay-berne output files
    CHARACTER(LEN=6) ZSTR
    CHARACTER(LEN=2)  DUMMYSTR
    DOUBLE PRECISION EulerPhi,EulerPsi,EulerTheta,EulerThetadeg,EulerPhiDeg,EulerPsiDeg  ! Euler angles for ellipsoids of revolution

    DOUBLE PRECISION EPS2, RAD, HEIGHT,sumx,sumy,sumz,CM(3),SYMOPS(120,3,3)
    LOGICAL :: GTEST
    COMMON /CAPS/ EPS2, RAD, HEIGHT

    CHARACTER(LEN=4) :: FPGRP
    CHARACTER(LEN=6) :: CRMS
    CHARACTER(LEN=20) :: MYFILENAME2, ISTR, DBNUM, MYFILENAME3
    CHARACTER(LEN=15), ALLOCATABLE :: DBNAME(:)

    CHARACTER(LEN=6) :: J1_STRING, MYNODE_STRING

    !  AMH
    CHARACTER(LEN=3) :: RES_TYPE
    CHARACTER(LEN=2) :: ATOM_TYPE
    CHARACTER*1 COUNTTT
    INTEGER COUNTT
    DOUBLE PRECISION  PPPCORD(NMRES*3*3,3,3,5)
    EXTERNAL NUM_TO_CHAR
    DOUBLE PRECISION TEND
    INTEGER ITERATIONS, BRUN,QDONE, NP
    DOUBLE PRECISION SCREENC(3*NATOMSALLOC), TIME

    PI = 4.D0*DATAN(1.D0)
    
    !ds656> test
    !write(*,*) "finalio> START, PI=", PI
    
    IF (PERMOPT) THEN
       NSAVE=2
       QMIN(2)=0.0D0
       QMINP(2,:)=FIN(:)
       FF(:)=0
    ENDIF

    NUMPRO = 1
    NUMCRD = 3

    ALLOCATE(DBNAME(NSAVE))
     
    IF (DEBUG) WRITE(MYUNIT,'(A,3I6)') ' in finalio MYNODE,MYUNIT,NSAVE=',MYNODE,MYUNIT,NSAVE !jdf43>
    DO J1=1,NSAVE
        WRITE(DBNUM,*) J1
        DBNAME(J1)='dbase.'//TRIM(ADJUSTL(DBNUM))
    ENDDO
    DCOORDS(1:3*NATOMS)=0.0D0
    !      IF (AMH) THEN
    !         CALL WALESAMH_FINALIO
    !         STOP
    !      ENDIF

    IF (AMHT) THEN
        AMHUNIT1=GETUNIT()
        OPEN(UNIT=AMHUNIT1,FILE='movie_gmin',STATUS='UNKNOWN')
        WRITE(AMHUNIT1,334)NMRES,NUMCRD,NUMPRO,NSAVE
334     FORMAT(4(I8,1X),' NMRES NMCRD NUMPRO NMSNAP')
    ENDIF

    IF (DMACRYST) THEN
        CALL DMACRYS_DUMP_LOWEST
    ENDIF

    IF(USERPOTT) THEN
        CALL USERPOT_DUMP_LOWEST
    ENDIF


    IF (MPIT) THEN
        WRITE (ISTR, '(I10)') MYNODE+1
        MYUNIT2=GETUNIT()
        MYFILENAME2="lowest."//TRIM(ADJUSTL(ISTR))
!        WRITE(MYUNIT,'(A,I6,A)') ' in finalio MYUNIT2,MYFILENAME2=',MYUNIT2,TRIM(ADJUSTL(MYFILENAME2))
        OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
        IF (CHRMMT) THEN
            DO J1=1,NSAVE
                WRITE(DBNUM,*) J1
                DBNAME(J1)='dbase.'//TRIM(ADJUSTL(ISTR))//'.'//TRIM(ADJUSTL(DBNUM))
            ENDDO
        ENDIF
        IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=GETUNIT()
            MYFILENAME3="CSMav."//trim(adjustl(istr))
        ENDIF
    ELSE
        MYUNIT2=GETUNIT()
        ! hk286
        IF (RIGIDINIT) THEN
            OPEN(MYUNIT2,FILE='GRlowest',STATUS='UNKNOWN')
        ELSE
            OPEN(MYUNIT2,FILE='lowest',STATUS='UNKNOWN')
        ENDIF
        ! hk286
        IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=26 
            OPEN(MYUNIT3,FILE='CSMav.xyz',STATUS='UNKNOWN')
        ENDIF
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! debug DJW
!
!        MYUNIT3=GETUNIT()
!        OPEN(MYUNIT3,FILE='GCBHlowest',STATUS='UNKNOWN')
!        MYUNIT4=GETUNIT()
!        OPEN(MYUNIT4,FILE='GCBHlowest.new',STATUS='UNKNOWN')
! 777    READ(MYUNIT3,*,END=888) J1
!        NATOMS=J1
!        READ(MYUNIT3,*) DUMMYSTR
!        NP=1
!        DO J2=1,J1
!           READ(MYUNIT3,'(A2,1X,3G20.10)') DUMMYSTR,COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)
!        ENDDO
! 
!        CALL QUENCH(.TRUE.,NP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!        WRITE(MYUNIT4,'(I10)') J1
!        WRITE(MYUNIT4,'(A,2G20.10)') 'energy=',QPE(J1),QENERGIES(J1)
!        DO J2=1,J1
!           WRITE(MYUNIT4,'(A2,1X,3G20.10)') 'AX',QCOORDINATES(J1,3*(J2-1)+1:3*(J2-1)+3)
!        ENDDO
! 
!        GOTO 777
! 888    CLOSE(MYUNIT3)
!        CLOSE(MYUNIT4)
!        STOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (GCBHT) THEN
       MYUNIT3=GETUNIT()
       OPEN(MYUNIT3,FILE='GCBHlowest',STATUS='UNKNOWN')
       DO J1=1,GCNATOMS
          IF (QENERGIES(J1).LT.HUGE(1.0D0)/10.0) THEN
             WRITE(MYUNIT3,'(I10)') J1
             WRITE(MYUNIT3,'(A,2G20.10)') 'energy=',QPE(J1),QENERGIES(J1)
             DO J2=1,J1
                WRITE(MYUNIT3,'(A2,1X,3G20.10)') 'AX',QCOORDINATES(J1,3*(J2-1)+1:3*(J2-1)+3)
             ENDDO
          ENDIF
       ENDDO
       CLOSE(MYUNIT3)
    ENDIF

    DO J1=1,NSAVE
        NATOMS=QMINNATOMS(J1)
        IF (AMHT) THEN
            COUNTT=J1
            CALL NUM_TO_CHAR(COUNTT,COUNTTT)
            OPEN(UNIT=27,FILE='movie_gmin.'//COUNTTT//'.pdb',STATUS='UNKNOWN')
        ENDIF

        IF (RGCL2.OR.ARNO) THEN
            WRITE(MYUNIT2,*) NATOMS+2
        ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
            WRITE(MYUNIT2,*) NATOMS/2
        ELSE IF (AMHT) THEN
            WRITE(MYUNIT2,*) NMRES*3
        ELSE
            WRITE(MYUNIT2,*) NATOMS
        ENDIF
        !        IF (CSMT.AND.DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'finalio> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
        !
        ! ds656> print point-group symmetry if required
        IF(PRINT_PTGRP) THEN
           WRITE(MYUNIT,'(A,I4)') 'finalio> Analysing point group of minimum', J1
           CM(1:3) = 0.0D0; NSYMOPS=0; SYMOPS(:,:,:) = 0.0D0
           !ds656> TEST
           J3=1 ! Determine majority label/species
           DO J2=2,NSPECIES(0)
              IF(NSPECIES(J2) > NSPECIES(J3)) J3 = J2
           ENDDO
           ! Determine the overall point-group (for ALL species!).
           CALL PGSYM( NATOMS, QMINP(J1,1:3*NATOMS), &
                QMINT(J1,1:NATOMS), PGSYMTOLS, 2, J3,&
                CM, NSYMOPS, SYMOPS, FPGRP )
           WRITE(MYUNIT2,99) J1, QMIN(J1), FF(J1), NPCALL_QMIN(J1), &
                FPGRP, NSYMOPS
99         FORMAT('Energy of minimum ',I6,'=',G20.10, &
                ' first found at step ',I8,' after ', &
                I20,' function calls. Point group ',A4, &
                ' of order ', I3,'.')
        ELSE ! <ds656
           WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1), NPCALL_QMIN(J1)
10         FORMAT('Energy of minimum ',I6,'=',G20.10, &
                ' first found at step ',I8,' after ',I20,' function calls')
        ENDIF
        !
        IF (MSORIGT.OR.FRAUSIT) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
20          FORMAT('Si',3F20.10)
        ELSE IF (MSTRANST) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
        ELSE IF (RGCL2) THEN
            WRITE(MYUNIT2,'(A,F20.10)') 'Cl 0.0 0.0 ', 0.995D0
            WRITE(MYUNIT2,'(A,F20.10)') 'Cl 0.0 0.0 ',-0.995D0
            WRITE(MYUNIT2,60) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
60          FORMAT('AR ',3F20.10)
        ELSE IF (AMHT) THEN
            !
            !   OUTPUT CORDS FOR LOWEST IN X,Y,Z FORMAT
            !   OUTPUT CORDS FOR MOVIE_GMIN IN MOVIESEG FORMAT
            !
            WRITE(AMHUNIT1,683)NUMCRD,j1,NUMCRD,REAL(NUMCRD),NUMCRD
683         FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

            GLY_COUNT = 0

            DO 1964 III = 1,NMRES
                IF (IRES(III).EQ.8) THEN
                    !!                pppcord(residue, xyz, numpro, atom types
                    PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CA X
                    PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                    PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
                    !    SWAP  CA for CB
                    PPPCORD(III, 1, 1, 2) = REAL(QMINP(j1,9*(III-1)+1- GLY_COUNT*3)) !  CB X
                    PPPCORD(III, 2, 1, 2) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CB Y
                    PPPCORD(III, 3, 1, 2) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CB Z
                    PPPCORD(III, 1, 1, 3) = REAL(QMINP(J1,9*(III-1)+4- GLY_COUNT*3)) !  O X
                    PPPCORD(III, 2, 1, 3) = REAL(QMINP(J1,9*(III-1)+5- GLY_COUNT*3)) !  O Y
                    PPPCORD(III, 3, 1, 3) = REAL(QMINP(J1,9*(III-1)+6- GLY_COUNT*3)) !  O Z

                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3), &
                    QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3), &
                    QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+4-GLY_COUNT*3),QMINP(J1,9*(III-1)+5-GLY_COUNT*3), &
                    QMINP(J1,9*(III-1)+6-GLY_COUNT*3)
31                  FORMAT('AM',3G25.15)

                    GLY_COUNT = GLY_COUNT +1
                ELSE
                    PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(iii-1)+1- GLY_COUNT*3)) !  CA X
                    PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                    PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
                    PPPCORD(III, 1, 1, 2) = REAL(QMINP(J1,9*(III-1)+4- GLY_COUNT*3)) !  CB X
                    PPPCORD(III, 2, 1, 2) = REAL(QMINP(J1,9*(III-1)+5- GLY_COUNT*3)) !  CB Y
                    PPPCORD(III, 3, 1, 2) = REAL(QMINP(J1,9*(III-1)+6- GLY_COUNT*3)) !  CB Z
                    PPPCORD(III, 1, 1, 3) = REAL(QMINP(J1,9*(III-1)+7- GLY_COUNT*3)) !  O X
                    PPPCORD(III, 2, 1, 3) = REAL(QMINP(J1,9*(III-1)+8- GLY_COUNT*3)) !  O Y
                    PPPCORD(III, 3, 1, 3) = REAL(QMINP(J1,9*(III-1)+9- GLY_COUNT*3)) !  O Z
 
                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3),&
                    QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+4-GLY_COUNT*3),QMINP(J1,9*(III-1)+5-GLY_COUNT*3),&
                    QMINP(J1,9*(III-1)+6-GLY_COUNT*3)
                    WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+7-GLY_COUNT*3),QMINP(J1,9*(III-1)+8-GLY_COUNT*3),&
                    QMINP(J1,9*(III-1)+9-GLY_COUNT*3)
                ENDIF
1964        CONTINUE

            DO 526 III=1,NMRES
                WRITE(AMHUNIT1,632)(PPPCORD(III,I3,1,1),I3=1,3),(PPPCORD(III,I3,1,2),I3=1,3),(PPPCORD(III,I3,1,3),I3=1,3)
632             FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1X))
            !632           FORMAT('CA: ',3(f25.15,1x),'CB: ',3(f25.15,1x),'Ox: ', 3(f25.15,1x))
526         CONTINUE

            DO III = 1+1, NMRES
                PPPCORD(III,1,1,4) = &
                0.4831806D0*PPPCORD(III-1,1,1,1) + 0.7032820D0*PPPCORD(III,1,1,1) - 0.1864626D0*PPPCORD(III-1,1,1,3)
                PPPCORD(III,2,1,4) = &
                0.4831806D0*PPPCORD(III-1,2,1,1) + 0.7032820D0*PPPCORD(III,2,1,1) - 0.1864626D0*PPPCORD(III-1,2,1,3)
                PPPCORD(III,3,1,4) = &
                0.4831806D0*PPPCORD(III-1,3,1,1) + 0.7032820d0*PPPCORD(III,3,1,1) - 0.1864626D0*PPPCORD(III-1,3,1,3)
            ENDDO

            DO III = 1, NMRES-1
                PPPCORD(III,1,1,5) = &
                0.4436538d0*PPPCORD(III,1,1,1)+0.2352006D0*PPPCORD(III+1,1,1,1)+0.3211455D0*PPPCORD(III,1,1,3)
                PPPCORD(III,2,1,5) = &
                0.4436538D0*PPPCORD(III,2,1,1)+0.2352006D0*PPPCORD(III+1,2,1,1)+0.3211455D0*PPPCORD(III,2,1,3)
                PPPCORD(III,3,1,5) = &
                0.4436538d0*PPPCORD(III,3,1,1)+0.2352006D0*PPPCORD(III+1,3,1,1)+0.3211455D0*PPPCORD(III,3,1,3)
            ENDDO

            DO III = 1, NMRES
                RES_TYPE = AMINOA(IRES(III))

                IF (III .NE. 1) THEN
                    ATOM_TYPE='N '
                    ID = 4
                    WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
61                  FORMAT('ATOM',4X,i3,2X,A2,2X,A3,3X,i3,4X,F8.3,F8.3,F8.3,2X,'1.00',2X,'0.00',6X,'TPDB',1x,I3)

                ENDIF

                ATOM_TYPE='CA'
                ID = 1
                WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III

                IF (RES_TYPE .NE. 'gly') THEN
                    ATOM_TYPE='CB'
                    ID = 2
                    WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
                ENDIF

                IF (III .NE. NMRES) THEN
                    ATOM_TYPE='C '
                    ID = 5
                    WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
                ENDIF

                ATOM_TYPE='O '
                ID = 3
                WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III

            ENDDO
            CLOSE(27)

        ELSE IF (ARNO) THEN
            WRITE(MYUNIT2,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(MYUNIT2,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(MYUNIT2,65) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
        ELSE IF (TOSI.OR.WELCH) THEN
            DO J2=1,NATOMS
                IF (ZSYM(J2).EQ.'PL') WRITE(MYUNIT2,'(A,3F20.10)') 'Na  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
                IF (ZSYM(J2).EQ.'MI') WRITE(MYUNIT2,'(A,3F20.10)') 'Cl  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
! hk286 - generalised Thomson problem, convert to Cartesians before printing
        ELSE IF (GTHOMSONT) THEN
            CALL GTHOMSONANGTOC(DCOORDS(1:3*NATOMS), QMINP(J1, 1:3*NATOMS), NATOMS) 
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(A,3F20.10)') 'C  ',(DCOORDS(3*(J2-1)+J3),J3=1,3)
            ENDDO
        ELSE IF ((BLJCLUSTER.OR.BLJCLUSTER_NOCUT).AND..NOT.QALCST) THEN
            DO J2=1,NATOMS
                IF (J2.LE.NTYPEA) THEN
                    WRITE(MYUNIT2,'(A,3F20.10)') 'LA  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
                ELSE
                    WRITE(MYUNIT2,'(A,3F20.10)') 'LB  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
                ENDIF
            ENDDO
        ELSE IF (GLJT) THEN
            J3=1
            J4=1
            DO J2=1,NATOMS
               !
               DO J5=1,3
                  P(J5) = QMINP(J1,3*(J2-1)+J5)                  
                  IF(PERIODIC) THEN ! wrap back into the box
                     P(J5)=P(J5)-BOX3D(J5)*ANINT(P(J5)/BOX3D(J5))
                  ENDIF
               ENDDO
               !
               WRITE(ZSTR,'(I6)') J4
               ZSTR='L' // TRIM(ADJUSTL(ZSTR))
               WRITE(MYUNIT2,'(A2,1X,3F20.10)') ZSTR,(P(J5),J5=1,3)
               J3=J3+1
               IF (J3.GT.NSPECIES(J4)) THEN
                  J3=1
                  J4=J4+1
               ENDIF
            ENDDO
            !
         ELSE IF (BGUPTAT .AND. .NOT. QALCST) THEN
            DO J2=1,NATOMS
                IF (J2.LE.NTYPEA) THEN
                    WRITE(MYUNIT2,'(A,3F20.10)') BGUPTANAME1, (QMINP(J1,3*(J2-1)+J3),J3=1,3)
                ELSE
                    WRITE(MYUNIT2,'(A,3F20.10)') BGUPTANAME2, (QMINP(J1,3*(J2-1)+J3),J3=1,3)
                ENDIF
            END DO


        ELSE IF (AMBER.OR.MOLECULART) THEN
            DO J2=1,NATOMS
                WRITE(MYUNIT2,'(A,3F20.10)') typech(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
        ELSE IF (AMBER12T) THEN
            ! Create a string for J1
            WRITE(J1_STRING,'(I6)') J1
            IF (DUMPSTRUCTURES) THEN
               IF (MPIT) THEN
                  ! If we're writing output for multiple parallel MPI runs we need to use MYNODE to
                  ! distinguish the outputs.
                  ! Create a string for the node too.
                  WRITE(MYNODE_STRING,'(I6)') MYNODE
                  CALL AMBER12_WRITE_RESTART(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                               &'.'//TRIM(ADJUSTL(MYNODE_STRING))//&
                                                               &'.rst', &
                                                       &  LEN('coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                                   &'.'//TRIM(ADJUSTL(MYNODE_STRING))//&
                                                                   &'.rst'))
                  CALL AMBER12_WRITE_PDB(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                           &'.'//TRIM(ADJUSTL(MYNODE_STRING))//&
                                                           &'.pdb', &
                                                    & LEN('coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                               &'.'//TRIM(ADJUSTL(MYNODE_STRING))//&
                                                               &'.pdb'))
               ELSE
                  ! Otherwise, just use one output without the MYNODE suffix.
                  CALL AMBER12_WRITE_RESTART(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                               &'.rst', &
                                                        & LEN('coords.'//TRIM(ADJUSTL(J1_STRING))//'.rst'))
                  CALL AMBER12_WRITE_PDB(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(J1_STRING))//&
                                                           &'.pdb', &
                                                    & LEN('coords.'//TRIM(ADJUSTL(J1_STRING))//'.pdb'))
               END IF
            END IF
            IF (MPIT) THEN
               WRITE(MYNODE_STRING,'(I6)') MYNODE
               ! Pass the name, length of the string and whether or not to write a header.
               CALL AMBER12_WRITE_XYZ(QMINP(J1,:), 'lowest.'//TRIM(ADJUSTL(MYNODE_STRING)), &
                                                   LEN('lowest.'//TRIM(ADJUSTL(MYNODE_STRING))), &
                                                   LOGICAL(.FALSE.,KIND=1))
            ELSE 
               ! Pass the name, length of the string and whether or not to write a header.
               CALL AMBER12_WRITE_XYZ(QMINP(J1,:), 'lowest', LEN('lowest'), LOGICAL(.FALSE.,KIND=1))
            END IF
        ELSE IF (AMBERT) THEN
            ! sf344> write out coordinates
            COORDS1(1:3*NATOMS) = QMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
                CALL AMBERFINALIO(j1,MYUNIT2,MYNODE+1,'0',0,COORDS1(1:3*NATOMS))
                WRITE(J1CHAR2,'(I3)') J1
             IF (MPIT) THEN
                WRITE (ISTR, '(I10)') MYNODE+1
                MYUNIT3=GETUNIT()
                MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))//'.'//trim(adjustl(istr))
                OPEN(MYUNIT3,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
             ELSE
                MYUNIT3=GETUNIT()
                MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))
                OPEN(MYUNIT3,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
             END IF
                DO J2=1,NATOMS
                    WRITE(MYUNIT3,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
                ENDDO
                CLOSE(MYUNIT3)
            ELSE
                DO I1=1,NATOMS
                    WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
                ENDDO
            ENDIF

        ELSE IF (CHIROT) THEN
            CALL CHIRO_OUTPUT

        ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
            !
            ! determine the centre of coordinates, then centre them
            !
            SUMX=0.0D0
            SUMY=0.0d0
            SUMZ=0.0d0

            DO J2=1,NATOMS/2
                SUMX=SUMX+QMINP(J1,3*(J2-1)+1)
                SUMY=SUMY+QMINP(J1,3*(J2-1)+2)
                SUMZ=SUMZ+QMINP(J1,3*(J2-1)+3)
            ENDDO
            SUMX=2*SUMX/NATOMS
            SUMY=2*SUMY/NATOMS
            SUMZ=2*SUMZ/NATOMS
            DO J2=1,NATOMS/2
                QMINP(J1,3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)-SUMX
                QMINP(J1,3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)-SUMY
                QMINP(J1,3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)-SUMZ
            ENDDO

            DO j2=1,NATOMS/2
                IF (PARAMONOVPBCX) THEN
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
                    QMINP(J1,3*J2-2)=QMINP(J1,3*J2-2)-BOXLX*NINT(QMINP(J1,3*J2-2)/BOXLX)
                ENDIF
                IF (PARAMONOVPBCY) THEN
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
                    QMINP(J1,3*J2-1)=QMINP(J1,3*J2-1)-BOXLY*NINT(QMINP(J1,3*J2-1)/BOXLY)
                ENDIF
                IF (PARAMONOVPBCZ) THEN
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
                    QMINP(J1,3*J2  )=QMINP(J1,3*J2  )-BOXLZ*NINT(QMINP(J1,3*J2  )/BOXLZ)
                ENDIF
            ENDDO

            DO J2=1,NATOMS/2
                WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a11,3f20.10)') 'H',QMINP(J1,3*(J2-1)+1), &
                QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3), &
                'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1), &
                QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
            ENDDO

        ELSE IF (CHRMMT) THEN
            DO J2=1,NATOMS
                WRITE(MYUNIT2,'(A,1X,3F20.10)') ZSYM(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            !       csw34> This DO loop appeared to be be missing on 30/9/08 which would easily
            !       explain the output problems!
            DO J2=1,NATOMS
                DCOORDS(3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)
                DCOORDS(3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)
                DCOORDS(3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)
            ENDDO
            CALL CHARMMDUMP(DCOORDS,DBNAME(J1))

        !    DC430 >
        !    |gd351> added patchy

        ELSE IF (CAPBINT.OR.DBPT.OR.DBPTDT.OR.DMBLMT.OR.DMBLPYT.OR.LINRODT.OR.LWOTPT.OR.MSTBINT.OR.MSSTOCKT.OR.NCAPT.OR.NPAHT &
        .OR. NTIPT .OR. STOCKAAT .OR. PAHAT .OR. PAHW99T .OR. TDHDT .OR. WATERDCT .OR. WATERKZT .OR. PATCHY .OR. PAPT   &
        .OR. PAPBINT .OR. PAPJANT .OR. PTSTSTT .OR. MORSEDPT) THEN
            DO J2 = 1, NATOMS/2
                WRITE(MYUNIT2,'(3f25.15)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2 = 1, NATOMS/2
                WRITE(MYUNIT2,'(3f25.15)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

        ELSE IF (GBT.OR.GBDT.OR.GBDPT.OR.MSGBT) THEN
            DO J2 = 1, NATOMS/2
                WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2 = 1, NATOMS/2
                WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (GEMT) THEN
            DO J2 = 1, NATOMS
                WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO

        ELSE IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
!
! this writes 'lowest' in xyz (Xmakemol) format
!
            WRITE(MYUNIT,'(A,I6,A)') ' in finalio BLN block MYUNIT2=',MYUNIT2
            DO J2=1,NATOMS
                WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
        ELSE IF(QALCST.OR.MULTIPERMT.OR.MLJT.OR.MGUPTAT.OR.MSCT) THEN
           !
           CALL SET_ATOMLISTS(QMINT(J1,1:NATOMS),1)
           IF(SPECMASST) THEN 
              EDUMMY = SUM(ATMASS)
              WRITE(MYUNIT,'(A,I6,A,E20.10)') &
                   'finalio> Min ',J1,' has total mass ', EDUMMY
           ENDIF
           !
           IF(STRESST) CALL CALC_STRESS(QMINP(J1,1:3*NATOMS),EDUMMY)
           !
           DO J2=1,NSPECIES(0) !ds656> Should not exceed 10
              IF(SPECLABELST) THEN
                 ATOM_TYPE = SPECLABELS(J2)
              ELSE
                 WRITE(ATOM_TYPE,'(I1)') J2
                 ATOM_TYPE='L' // TRIM(ADJUSTL(ATOM_TYPE))
              ENDIF
              DO J3=1,ATOMLISTS(J2,1,0)
                 J4=ATOMLISTS(J2,1,J3) ! Actual atom index
                 DO J5=1,3
                    P(J5) = QMINP(J1,3*(J4-1)+J5)                  
                    IF(PERIODIC) THEN ! wrap back into the box
                       P(J5)=P(J5)-BOX3D(J5)*ANINT(P(J5)/BOX3D(J5))
                    ENDIF
                 ENDDO
                 WRITE(MYUNIT2,'(A,3(1X,F20.10))',ADVANCE='NO') &
                      ATOM_TYPE,(P(J5), J5=1,3)
                 IF(STRESST) THEN ! Tack on local stresses
                    IF(STRESS_MODE==2) THEN
                       ! Just print the uniqe elements of the stress
                       ! tensor
                       WRITE(MYUNIT2,'(6(1X,F10.5))',ADVANCE='YES') &
                            STRESS(J4,1,1),STRESS(J4,2,2),STRESS(J4,3,3),&
                            STRESS(J4,1,2),STRESS(J4,1,3),STRESS(J4,2,3)
                    ELSEIF(STRESS_MODE==1) THEN
                       WRITE(MYUNIT2,'(2(1X,F10.5))',ADVANCE='YES') &
                            -(STRESS(J4,1,1)+STRESS(J4,2,2)+STRESS(J4,3,3))/3.0D0, &
                            DSQRT( (STRESS(J4,1,1)-STRESS(J4,2,2))**2 + &
                            (STRESS(J4,2,2)-STRESS(J4,3,3))**2 + &
                            (STRESS(J4,3,3)-STRESS(J4,1,1))**2 )/3.0D0
                    ELSE
                       WRITE(MYUNIT2,'(1X,F10.5)',ADVANCE='YES') &
                            -(STRESS(0,1,1)+STRESS(0,2,2)+STRESS(0,3,3))/&
                            DBLE(3*NATOMS)
                    ENDIF
                 ELSE ! Just advance to next line
                    WRITE(MYUNIT2,*)
                 ENDIF
              ENDDO
           ENDDO
        ELSE
            IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
                WRITE(MYUNIT3,'(I6)') NATOMS
                WRITE(MYUNIT3,'(A,I6,2(A,G20.10))') 'averaged structure for final solution ',J1, &
                ' CSM=',QMINAV(J1),' CSM for reference structure=',QMIN(J1)
                WRITE(MYUNIT3,30) (QMINPCSMAV(J1,J2),J2=1,3*(NATOMS-NS))
            ENDIF
            WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30          FORMAT('SI ',3F20.10)
        ENDIF

        !|gd351>
        IF (ASAOOS) THEN
            LUNIT=GETUNIT()
            OPEN(LUNIT,file='particles.xyz')
            WRITE(LUNIT,*) NATOMS
            WRITE(LUNIT,*) ' '
            DO J2=1,NATOMS
                WRITE(LUNIT,'(A,3F20.10)') 'H  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            CLOSE(LUNIT)
        END IF
        !<gd351|

        IF ((NS.GT.0).AND.(.NOT.(WELCH.OR.TOSI))) THEN
            IF (MSORIGT.OR.FRAUSIT) THEN
                WRITE(MYUNIT2,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
40              FORMAT('Si',3F20.10)
            ELSE IF (MSTRANST) THEN
                WRITE(MYUNIT2,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
            ELSE
                WRITE(MYUNIT2,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
50              FORMAT('LB',3F20.10)
            ENDIF
        ENDIF
        IF (AMBER) CALL AMBERDUMP(J1,QMINP)
    ENDDO

!
! End of loop over dump to file lowest or equivalent.
!
    CLOSE(MYUNIT2)
    IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) CLOSE(MYUNIT3)
    !
    !     csw34> New loop for dumping interaction energy files if A9INTE is specified
    !     Added the missing IF block to test for A9INTE 9/12/09 DJW
    !
    IF (A9INTET) THEN
        IF (MPIT) THEN
            WRITE (ISTR, '(I10)') MYUNIT-22980+1
            MYUNIT2=(MYUNIT-22980+1)+100
            MYFILENAME2="intelowest."//trim(adjustl(istr))
            OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
        ELSE
            MYUNIT2=25
            OPEN(MYUNIT2,FILE='intelowest',STATUS='UNKNOWN')
        ENDIF
    ENDIF
    !
    !     csw34> loop structure copied from the ELSEIF(AMBERT) block above
    !
    IF (A9INTET.AND.AMBERT) THEN
        DO J1=1,NSAVEINTE
            WRITE(MYUNIT2,*) NATOMS
            !     csw34> write header to intelowest for current minimum
            WRITE(MYUNIT2,10) J1, INTEQMIN(J1), INTEFF(J1)
            !     sf344> write out coordinates
            COORDS1(1:3*NATOMS) = INTEQMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
                CALL INTEFINALIO(j1,MYUNIT2,AMBFINALIO_NODE,'0',0,COORDS1(1:3*NATOMS))
                WRITE(J1CHAR2,'(I3)') J1
                WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
                OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')

                DO J2=1,NATOMS
                    WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
                ENDDO
                CLOSE(226)
                WRITE(J1CHAR2,'(I3)') J1
                WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
                OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
              
                DO J2=1,NATOMS
                    WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
                ENDDO
                CLOSE(226)

            ELSE
                DO I1=1,NATOMS
                    WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
                ENDDO
            ENDIF
        ENDDO
    ENDIF
    !
    !  End of loop over dump to file intelowest
    !

! csw34> Output for the HBONDMATRIX keyword (subroutine in hbondmatrix.f90)
    IF ((AMBERT.OR.AMBER12T.OR.(CUDAT.AND.(CUDAPOT.EQ.'A'))).AND.HBONDMATRIX) CALL HBONDMATRIXFINALIO()

    CLOSE(MYUNIT2)

    !     csw34> Edits to the RMS keyword
    IF (CHRMMT.AND.RMST) THEN
        !        IF (PROGRESS) THEN
        !           DCOORDS(1:3*NATOMS)=RMSCOOR(1,1:3*NATOMS)
        !           IF(RMSBEST(1,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'closestrms')
        !           WRITE(MYUNIT,'(A9,F8.5)') 'RMSDmin= ',RMSBEST(1,1)
        !        ELSE
        ! REMEMBER TO RE-INDENT THE BELOW IF UNCOMMENTING ABOVE!
        OPEN(UNIT=MYUNIT2,FILE='rmsbest.'//TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
        DO J2=1,RMSSAVE
            WRITE(MYUNIT2,'(I6,F6.3,F15.5)')J2,RMSBEST(J2,1),RMSBEST(J2,2)
            WRITE(CRMS,'(I6)') J2
            DCOORDS(1:3*NATOMS)=RMSCOOR(J2,1:3*NATOMS)
            IF(RMSBEST(J2,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'rms.'//TRIM(ADJUSTL(ISTR))//'.'//TRIM(ADJUSTL(CRMS)))
        ENDDO
        CLOSE(MYUNIT2)
    !        ENDIF
    ENDIF

    IF (LJCOULT) THEN
        OPEN(UNIT=AMHUNIT1,FILE='ljcoul.xyz',STATUS='UNKNOWN')
        DO J1=1,NSAVE
            WRITE(AMHUNIT1,'(I6)') NATOMS
            WRITE(AMHUNIT1,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS
                !              Use "O" atom type to highlight charged particles and "N" for the neutral ones.
                IF (J2.LE.COULN) THEN
                    WRITE(AMHUNIT1,'(A4,3F18.10,A12,3F18.10)') 'O ',&
                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
                ELSE
                    WRITE(AMHUNIT1,'(A4,3F18.10,A12,3F18.10)') 'N ',&
                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
                END IF
            ENDDO
        ENDDO
        CLOSE(AMHUNIT1)


    ELSE IF (STOCKT) THEN
        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT,FILE='stock.xyz',STATUS='UNKNOWN')
        DO J1=1,NSAVE
            WRITE(LUNIT,'(I6)') NATOMS/2
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
                WRITE (LUNIT,'(A4,3F18.10,A12,3F18.10)') 'LA ', &
                QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3), &
                'atom_vector', &
                SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)), &
                SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)), &
                COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))
            ENDDO
        ENDDO
        CLOSE(LUNIT)

    ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT.OR.PYT) THEN
        DO J1=1,NSAVE
            WRITE(J1CHAR2,'(I3)') J1
            IF (MPIT) THEN
                WRITE (ISTR, '(I10)') MYNODE+1
                MYUNIT2=GETUNIT()
                MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))//'.'//trim(adjustl(istr))
                OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
            ELSE
                MYUNIT2=GETUNIT()
                MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))
                OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
            END IF

            DO J2=1,NATOMS
                WRITE(MYUNIT2,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
            ENDDO

            CLOSE(MYUNIT2)

        ENDDO
!
!  Write out lowest NSAVE structures to xmakemol ellipsoid format for
!  clusters of ellipsoids of revolution
!

        IF (MPIT) THEN
            WRITE (ISTR, '(I10)') MYNODE+1
            MYUNIT2=GETUNIT()
            MYFILENAME2="ellipsoid."//trim(adjustl(istr))//".xyz"
            OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
        ELSE
            MYUNIT2=GETUNIT()
            OPEN(MYUNIT2,FILE='ellipsoid.xyz',STATUS='UNKNOWN')
        ENDIF

        IF(PYT) THEN
            CALL PY_OUTPUT(MYUNIT2)
        ELSE

        do J1=1,NSAVE
            WRITE(MYUNIT2,*) NATOMS/2
            WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)

            DO J2=1,NATOMS/2
               
                IF (GAYBERNET) THEN
                    CALL EllipsoidsAAtoPolar(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhi,EulerPsi,EulerTheta,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
!  EulerPhiDeg = 90-EulerPhiDeg    ! EulerPhiDeg returned from EllipsoidsAAtoPolar is in fact the alpha angle
!                                  ! defined by me (angle of vector with the xy plane)
                    WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') '0',QMINP(J1,3*(J2-1)+1),&
                    QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                    'ellipse ',1.0D0,1.0D0,GBANISOTROPYR,&
                    EulerPsiDeg,EulerPhiDeg,0.0D0,&! this is in degrees
                    'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
                    CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                    WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),&
                    QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                    'ellipse ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,&
                    'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                ELSE IF (GBT.OR.GBDT) THEN
                    CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                    WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),&
                    QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                    'ellipse ',GBKAPPA,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,&
                    'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                ELSE IF (LJCAPSIDT) THEN
                    CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                    WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),&
                    QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                    'ellipse ',1.0D0,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,&
                    'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),&
                    QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                ENDIF
            ENDDO
        ENDDO

        ENDIF ! PYT

        CLOSE(MYUNIT2)

    ELSE IF (TIP) THEN
        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT,FILE='tip.xyz',STATUS='UNKNOWN')
        DO J1=1,NSAVE
            WRITE(LUNIT,'(I6)') (NATOMS/2)*3
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
                CALL TIPIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),&
                RBCOORDS)
                WRITE(LUNIT,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
                WRITE(LUNIT,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
                WRITE(LUNIT,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
        ENDDO
        CLOSE(LUNIT)
    ELSE IF (CAPSID) THEN
        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT,FILE='capsid.xyz',STATUS='UNKNOWN')
        DO J1=1,NSAVE
            WRITE(LUNIT,'(I6)') (NATOMS/2)*6
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
                CALL CAPSIDIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),&
                RBCOORDS,RAD,HEIGHT)
                DO J3=1,5
                    WRITE(LUNIT,'(A4,3F20.10)') 'C1 ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
                ENDDO
                WRITE(LUNIT,'(A4,3F20.10)') 'C4  ',RBCOORDS(16),RBCOORDS(17),RBCOORDS(18)
            ENDDO
        ENDDO
        CLOSE(LUNIT)

    ELSE IF (DBPT .OR. DMBLPYT) THEN

        LUNIT=GETUNIT()
        IF (DBPT) OPEN(UNIT=LUNIT, FILE='dbp.xyz', STATUS='UNKNOWN')
        IF (DMBLPYT) OPEN(UNIT=LUNIT, FILE='dmblpy.xyz', STATUS='UNKNOWN')

        GTEST = .FALSE.
        DU    = (/0.D0, 1.D0, 0.D0/)
        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES 
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                P(:) = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                DO J4 = 1, NRBSITES

                    IF (J4 == 1) THEN
                        RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                        WRITE(LUNIT,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSEIF (J4 == 2) THEN
                        RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                        WRITE(LUNIT,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSE
                        RBCOORDS(1:3) = MATMUL(RMI(:,:),DU)
                        WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                        'H', QMINP(J1,J3-2), QMINP(J1,J3-1), QMINP(J1,J3),&
                        'atom_vector', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDIF

                ENDDO

            ENDDO

        ENDDO
        CLOSE(LUNIT)

        RETURN

    ELSE IF (DBPTDT) THEN

        CALL VIEWDMBLTD()
        RETURN

    ELSE IF (DMBLMT) THEN

        CALL VIEWDMBL()
        RETURN

    ELSE IF (GBT .OR. GBDT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='gbe.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.
        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') NATOMS/2
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                RBCOORDS(1:3) = QMINP(J1,J3-2:J3)
                P(:)          = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                PHI   = DATAN2(RMI(2,3),RMI(1,3))
                IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI

                THT   = DACOS(RMI(3,3))

                PHI   = PHI*180.D0/PI
                THT   = THT*180.D0/PI

                WRITE(LUNIT,'(a5,2x,3f20.10,2x,a8,6f20.10)')&
                'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),&
                'ellipse', 2.D0*ESA(1), 2.D0*ESA(2), 2.D0*ESA(3), PHI, THT, 0.D0

            ENDDO

        ENDDO

        RETURN

    ELSE IF (LINRODT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='linrod.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3    = 3*J2
                J5    = 3*NATOMS/2 + J3
                P(:)  = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)


                DO J4 = 1, NRBSITES

                    RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

                ENDDO

                DO J4 = 1, NRBSITES

                    J3 = J4 + 1
                    IF (J4 == NRBSITES) J3 = 1
                    P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                    WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                    'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

                ENDDO

            ENDDO

        ENDDO

        CLOSE(UNIT=LUNIT)

        RETURN

    ELSE IF (LWOTPT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='lwotp.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.
         
        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3    = 3*J2
                J5    = 3*NATOMS/2 + J3
                P(:)  = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)


                DO J4 = 1, NRBSITES

                    RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

                ENDDO

                DO J4 = 1, NRBSITES

                    J3 = J4 + 1
                    IF (J4 == NRBSITES) J3 = 1
                    P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                    WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                    'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

                ENDDO

            ENDDO

        ENDDO

        CLOSE(UNIT=LUNIT)

        RETURN

    ELSE IF (MSTBINT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='mstbin.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') NPS*NRBSITES1 + (NATOMS/2 - NPS)*(NRBSITES - NRBSITES1)
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3    = 3*J2
                J5    = 3*NATOMS/2 + J3
                P(:)  = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                IF (J2 <= NPS) THEN
                    NRBS1 = 1
                    NRBS2 = NRBSITES1
                ELSE
                    NRBS1 = NRBSITES1 + 1
                    NRBS2 = NRBSITES
                ENDIF

                DO J4 = NRBS1, NRBS2

                    RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

                ENDDO

                DO J4 = NRBS1, NRBS2

                    J3 = J4 + 1
                 
                    IF (J2 <= NPS) THEN
                        IF (J4 == NRBSITES1) J3 = 1
                    ELSE
                        IF (J4 == NRBSITES) J3 = NRBSITES1 + 1
                    ENDIF
 
                    P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                    IF (J2 <= NPS) THEN
                        WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                        'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)
                    ELSE
                        WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                        'N', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)
                    ENDIF

                ENDDO

            ENDDO

        ENDDO

        CLOSE(UNIT=LUNIT)

        RETURN

    ELSE IF (MSSTOCKT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='msstock.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.
         
        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                P(:) = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                DO J4 = 1, NRBSITES

                    RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                    P(:)          = MATMUL(RMI(:,:),RBUV(J4,:))
                    WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                    'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3), 'atom_vector', P(1), P(2), P(3)

                ENDDO

            ENDDO

        ENDDO

        CLOSE(UNIT=LUNIT)

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='msstktr.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES !(NRBSITES - 1)
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3    = 3*J2
                J5    = 3*NATOMS/2 + J3
                P(:)  = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                DO J4 = 1, NRBSITES

                    RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

                ENDDO

                DO J4 = 1, NRBSITES !- 1

                    J3 = J4 + 1
                    IF (J4 == NRBSITES) J3 = 1
                    !                  IF (J4 == NRBSITES - 1) J3 = 1
                    P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                    WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)')&
                    'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

                ENDDO

            ENDDO

        ENDDO

        CLOSE(UNIT=LUNIT)

        RETURN

    ELSE IF (MULTPAHAT) THEN

        CALL VIEWMULTPAHA()
        RETURN

    ELSE IF (NPAHT .OR. PAHAT .OR. PAHW99T) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='rigid.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.

        IF (PAHW99T) THEN
            NCPHST = NCARBON + (NRBSITES-NCARBON)/2
            DO J1 = 1, (NCPHST-NCARBON) 
                SITE(NCARBON+J1,:) = SITE(NCPHST+J1,:)
            ENDDO
        ELSE
            NCPHST = NRBSITES
        ENDIF

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                P(:) = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                DO J4 = 1, NCPHST

                    RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                    IF (J4 <= NCARBON) THEN
                        WRITE(LUNIT,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSE
                        WRITE(LUNIT,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDIF

                ENDDO

            ENDDO

        ENDDO

        RETURN
    ELSE IF (CAPBINT) THEN
        CALL VIEWCAPBIN()
        RETURN

    ELSE IF (NCAPT) THEN

        CALL VIEWNEWCAPSID()
        RETURN

    ELSE IF (NTIPT) THEN

        CALL VIEWNEWTIP()
        RETURN

    ELSE IF (MWFILMT) THEN

        CALL MWDRAW(NATOMS, PERIODIC, LAT)
        RETURN

    ELSE IF (PAPT) THEN

        CALL VIEWPAP()
        RETURN

    ELSE IF (PAPBINT) THEN

        CALL VIEWPAPBIN()
        RETURN

    ELSE IF (PAPJANT) THEN

        CALL VIEWPAPJANUS()
        RETURN

    ELSE IF (PTSTSTT) THEN

        CALL VIEWPTSTST()
        RETURN

    !|gd351>

    ELSE IF (PATCHY) THEN

        CALL VIEWPATCHY()
        RETURN

    !<gd351|
    ELSE IF (SANDBOXT) THEN

        CALL SANDBOX_OUTPUT()

    ELSE IF (SILANET) THEN

        CALL VIEWSILANE()
        RETURN

    ELSE IF (STOCKAAT .OR. MORSEDPT) THEN

        LUNIT=GETUNIT()
        IF (STOCKAAT) OPEN(UNIT=LUNIT, FILE='stockaa.xyz', STATUS='UNKNOWN')
        IF (MORSEDPT) OPEN(UNIT=LUNIT, FILE='morsedp.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.
        DU    = (/0.D0, 0.D0, 1.D0/)

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') NATOMS/2
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                P(:) = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                RBCOORDS(1:3) = MATMUL(RMI(:,:),DU(:))
                WRITE(LUNIT,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'O', QMINP(J1,J3-2), QMINP(J1,J3-1), QMINP(J1,J3),&
                'atom_vector', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)


            ENDDO

        ENDDO

        CLOSE (UNIT=LUNIT)

        RETURN

    ELSE IF (TDHDT) THEN

        CALL VIEWTDHD()
        RETURN

    ELSE IF (DDMT) THEN

        CALL VIEWDDM()
        RETURN

    ELSE IF (WATERDCT .OR. WATERKZT) THEN

        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT, FILE='rigid.xyz', STATUS='UNKNOWN')
        GTEST = .FALSE.

        DO J1 = 1, NSAVE

            WRITE(LUNIT,'(I6)') (NATOMS/2)*(NRBSITES - 1)
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

                J3   = 3*J2
                J5   = 3*NATOMS/2 + J3
                P(:) = QMINP(J1,J5-2:J5)

                CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

                DO J4 = 1, NRBSITES - 1

                    RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                    IF (J4 == 1) THEN
                        WRITE(LUNIT,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSE
                        WRITE(LUNIT,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDIF

                ENDDO

            ENDDO

        ENDDO

        RETURN

    ELSE IF (RIGID) THEN
        LUNIT=GETUNIT()
        OPEN(UNIT=LUNIT,FILE='rigid.xyz',STATUS='UNKNOWN')
        DO J1=1,NSAVE
            WRITE(LUNIT,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(LUNIT,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
                CALL RBIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),&
                QMINP(J1,3*(NATOMS/2+J2-1)+1),&
                QMINP(J1,3*(NATOMS/2+J2-1)+2),&
                QMINP(J1,3*(NATOMS/2+J2-1)+3),&
                RBCOORDS,NRBSITES,SITE)
                DO J3=1,NRBSITES
                    WRITE(LUNIT,'(A4,3F20.10)') 'LA ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
                ENDDO
            ENDDO
        ENDDO
        CLOSE(LUNIT)

    ENDIF

    IF(ALLOCATED(DBNAME)) DEALLOCATE(DBNAME)

    IF (AMBER12T) THEN
       CALL AMBER12_FINISH()
    END IF

    CALL CPU_TIME(TEND)
    WRITE(MYUNIT,"(A,F18.1,A)") "time elapsed ", TEND - TSTART, " seconds"
    WRITE(MYUNIT,"(A,I18)") "Number of potential calls ", NPCALL

    RETURN
END SUBROUTINE FINALIO

SUBROUTINE AMBERDUMP(J1,QMINP)
    USE COMMONS
    USE MODAMBER
    IMPLICIT NONE


    CHARACTER(LEN=25) COORDFILE
    CHARACTER(LEN=2) FNAME
    INTEGER J1, LUNIT, GETUNIT
    DOUBLE PRECISION QMINP(NSAVE,3*NATOMS)

    IF (J1.LT.10) THEN
        WRITE (FNAME,'(I1)') J1
    ELSE
        WRITE (FNAME,'(I2)') J1
    ENDIF

    DO A=1,ATOMS
        X(A)=QMINP(J1,3*A-2)
        Y(A)=QMINP(J1,3*A-1)
        Z(A)=QMINP(J1,3*A)
    END DO

    COORDFILE='acoords.dump.'//FNAME

 
    LUNIT=GETUNIT()
    OPEN (UNIT=LUNIT,IOSTAT=IOS,FILE=coordfile,STATUS='UNKNOWN')

    DO a=1,atoms
        WRITE (UNIT=LUNIT,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F7.3,3X,F7.3,3X,F7.3)') label(a),typech(a),&
        a,bondedto(a),x(a),y(a),z(a)
    ENDDO

    WRITE (UNIT=LUNIT,FMT='(A3)') 'end'
    WRITE (UNIT=LUNIT,FMT='(A)') ' '
    WRITE (UNIT=LUNIT,FMT='(A4,7X,I2)') 'loop',rings

    DO a=1,rings
        WRITE (UNIT=LUNIT,FMT='(I3,4X,I3)') loopatom(2*a-1),loopatom(2*a)
    END DO

    WRITE (UNIT=LUNIT,FMT='(A)') ' '
    WRITE (UNIT=LUNIT,FMT='(A7)') 'charges'

    DO a=1,atoms
        q(a)=q(a)/18.2223
        WRITE (UNIT=LUNIT,FMT='(I3,2X,F7.4)') a,q(a)
    END DO

    WRITE (UNIT=LUNIT,FMT='(A3)') 'end'
    CLOSE(LUNIT)

    RETURN

END SUBROUTINE AMBERDUMP
!
!  SUBROUTINE to convert capsid CofM and DV coordinates to penatgons.
!
SUBROUTINE CAPSIDIO(X1, Y1, Z1, L1, M1, N1,COORDS,RAD,HEIGHT)
    IMPLICIT NONE
    DOUBLE PRECISION X1, Y1, Z1, COORDS(*), HEIGHT, C2A1,&
    M1, L1, N1, ALPHA1, RAD, CA1, S1, C3A1,&
    NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

    NUM1=-(1.0D0+SQRT(5.0D0))/4.0D0
    NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)/2.0D0
    NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)/2.0D0
    NUM4=(SQRT(5.0D0)-1.0D0)/4.0D0
    NUM5=-(1.0D0+SQRT(5.0D0))/4.0D0

    L12=L1**2
    M12=M1**2
    N12=N1**2
    ALPHA1=SQRT(L12+M12+N12)
    CA1=COS(ALPHA1)
    C2A1=RAD*CA1
    IF (ALPHA1.LT.0.0001D0) THEN
        !        C3A1=RAD*(-ALPHA1/2+ALPHA1**3/24)
        C3A1=RAD*(-0.5D0+ALPHA1**2/24.0D0)
        S1=RAD*(1.0D0-ALPHA1**2/6)
    ELSE
        C3A1=RAD*(CA1-1.0D0)/ALPHA1**2
        S1=RAD*SIN(ALPHA1)/ALPHA1
    ENDIF

    COORDS(1) =     c2a1 - c3a1*l12 + x1
    COORDS(2) =     -(c3a1*l1*m1) - n1*s1 + y1
    COORDS(3) =     -(c3a1*l1*n1) + m1*s1 + z1
    COORDS(4) =     c2a1*num4 - c3a1*l1*(m1*num3 + l1*num4) + n1*num3*s1 + x1
    COORDS(5) =     c2a1*num3 - c3a1*m1*(m1*num3 + l1*num4) - n1*num4*s1 + y1
    COORDS(6) =     -(c3a1*n1*(m1*num3 + l1*num4)) - l1*num3*s1 + m1*num4*s1 + z1
    COORDS(7) =     c2a1*num1 - c3a1*l1*(l1*num1 + m1*num2) + n1*num2*s1 + x1
    COORDS(8) = c2a1*num2 - c3a1*m1*(l1*num1 + m1*num2) - n1*num5*s1 + y1
    COORDS(9) = -(c3a1*n1*(l1*num1 + m1*num2)) + m1*num1*s1 - l1*num2*s1 + z1
    COORDS(10) = c2a1*num1 + c3a1*l1*(-(l1*num1) + m1*num2) - n1*num2*s1 + x1
    COORDS(11) = -(c2a1*num2) + c3a1*m1*(-(l1*num1) + m1*num2) - n1*num5*s1 + y1
    COORDS(12) = -(c3a1*l1*n1*num1) + c3a1*m1*n1*num2 + m1*num1*s1 + l1*num2*s1 + z1
    COORDS(13) = c2a1*num4 + c3a1*l1*(m1*num3 - l1*num4) - n1*num3*s1 + x1
    COORDS(14) = -(c2a1*num3) + c3a1*m1*(m1*num3 - l1*num4) - n1*num4*s1 + y1
    COORDS(15) = c3a1*n1*(m1*num3 - l1*num4) + l1*num3*s1 + m1*num4*s1 + z1
    !     COORDS(16)= (-(c3a1*l1*n1) - m1*s1 + 2*x1)/2.
    !     COORDS(17)= -(c3a1*m1*n1)/2. + (l1*s1)/2. + y1
    !     COORDS(18)= (c2a1 - c3a1*n12 + 2*z1)/2.
    COORDS(16)= -(c3a1*height*l1*n1) - height*m1*s1 + x1
    COORDS(17)= -(c3a1*height*m1*n1) + height*l1*s1 + y1
    COORDS(18)= c2a1*height - c3a1*height*n12 + z1

    RETURN
END SUBROUTINE CAPSIDIO
!
!  Subroutine to convert rigid body CofM and DV coordinates to molecular sites.
!
SUBROUTINE RBIO(X1, Y1, Z1, L1, M1, N1, COORDS, NRBSITES, SITE)
    IMPLICIT NONE
    INTEGER NRBSITES
    DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, SITE(NRBSITES,3),&
    M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12
    INTEGER J1

    L12=L1**2
    M12=M1**2
    N12=N1**2
    ALPHA1=SQRT(L12+M12+N12)
    CA1=COS(ALPHA1)
    C2A1=CA1
    IF (ALPHA1.LT.0.0001D0) THEN
        !        C3A1=(-ALPHA1/2+ALPHA1**3/24)
        C3A1=(-0.5D0+ALPHA1**2/24.0D0)
        S1=(1.0D0-ALPHA1**2/6)
    ELSE
        C3A1=(CA1-1.0D0)/ALPHA1**2
        S1=SIN(ALPHA1)/ALPHA1
    ENDIF
   
    DO J1=1,NRBSITES
        COORDS(3*(J1-1)+1)=c2a1*SITE(J1,1) + s1*(n1*SITE(J1,2) - m1*SITE(J1,3)) - &
        c3a1*l1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + X1
        COORDS(3*(J1-1)+2)=c2a1*SITE(J1,2) + s1*(-(n1*SITE(J1,1)) + l1*SITE(J1,3)) &
        - c3a1*m1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Y1
        COORDS(3*(J1-1)+3)=s1*(m1*SITE(J1,1) - l1*SITE(J1,2)) + c2a1*SITE(J1,3) &
        - c3a1*n1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Z1
    ENDDO

    RETURN
END SUBROUTINE RBIO
!
!  SUBROUTINE to convert TIP oxygen and DV coordinates to Cartesians.
!
SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
    IMPLICIT NONE
    DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

    L12=L1**2
    M12=M1**2
    N12=N1**2
    ALPHA1=SQRT(L12+M12+N12)
    CA1=COS(ALPHA1)
    C2A1=CA1
    IF (ALPHA1.LT.0.0001D0) THEN
        !        C3A1=(-ALPHA1/2+ALPHA1**3/24)
        C3A1=(-0.5D0+ALPHA1**2/24.0D0)
        S1=(1.0D0-ALPHA1**2/6)
    ELSE
        C3A1=(CA1-1.0D0)/ALPHA1**2
        S1=SIN(ALPHA1)/ALPHA1
    ENDIF

    COORDS(1) = X1
    COORDS(2) = Y1
    COORDS(3) = Z1
    COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
    COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
    COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
    COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
    COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
    COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

    RETURN
END SUBROUTINE TIPIO
