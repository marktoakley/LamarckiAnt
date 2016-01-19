!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2012 David J. Wales
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
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!   USA
!

! csw34> This file contains the functions and subroutines used by the HBONDMATRIX keyword
! GETHBONDMATRIX(TOPFILE,RSTFILE,DONORFILE,RESFILE,MATRIX)
! Calls hbond_matrix.sh (found in SCRIPTS/AMBER) to return a hydrogen-bond matrix
! INTEGER FUNCTION SUMSQUAREDIFF(A,B,WIDTH)
! Takes two matricies (A and B) of dimension WIDTH and returns the sum of the squared differences of their elements
! SUBROUTINE DUMPRST(COORDS,NATOMS,FNAME)
! Dumps an AMBER restart format (no velocities) file FNAME.rst for the supplied coordinates
! SUBROUTINE WRITEHBONDMATRIX(MATRIX,NRES,FNAME)
! Dumps the provided hydrogen-bond matrix of size NRES to file FNAME.mat
! SUBROUTINE HBONDMATRIXFINALIO()
! Writes the final output for the hydrogen-bond matrix run - including concaternated quenchX.pdb files for each group if DUMPQU is used
!
! It currently only works with AMBER - but that is easy enough to change.
! You would only need to write an equivilent script to hbond_matrix.sh for
! CHARMM say, and a new subroutine to call it like the one below. It would then
! be a case of just putting some ELSE IF bits in mc.F when the matrix is
! constructed - everything else could remain the same.

SUBROUTINE GETHBONDMATRIX(TOPFILE,RSTFILE,DONORFILE,RESFILE,DCUT,ACUT,MATRIX)
    ! Subroutine to call the hbond_matrix.sh script for use in GMIN
    ! Takes an AMBER topology file (TOPFILE), an AMBER restart file (RSTFILE), a
    ! file containing the ptraj format list of donor and acceptor atoms (DONORFILE)
    ! and a file containing the list of residues you're interested in (RESFILE).

    ! You should give it a matrix of the size required (MATRIX)
    USE PORFUNCS
    USE COMMONS, ONLY : HBONDNRES
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: TOPFILE, RSTFILE, DONORFILE, RESFILE, DCUT, ACUT
    INTEGER, INTENT(INOUT) :: MATRIX(HBONDNRES,HBONDNRES)
    INTEGER :: ISTAT

    ! DEBUG PRINTING
    !PRINT *,"In GETHBONDMATRIX - TOPFILE,RSTFILE,DONORFILE,RESFILE:"
    !PRINT *,TOPFILE,RSTFILE,DONORFILE,RESFILE

    ! Call the script which calculates the matrix (this could be replaced)
    CALL SYSTEM_SUBR('bash hbond_matrix.sh ' // TRIM(ADJUSTL(TOPFILE)) // ' ' &
    // TRIM(ADJUSTL(RSTFILE))  // ' ' // TRIM(ADJUSTL(DONORFILE)) &
    // ' ' // TRIM(ADJUSTL(RESFILE)) // ' ' // TRIM(ADJUSTL(DCUT)) // ' ' &
    // TRIM(ADJUSTL(ACUT)) // ' > temp.mat',ISTAT)

    ! Read matrix into MATRIX
    OPEN(UNIT=20,FILE='temp.mat',STATUS='OLD')
    READ(20,*) MATRIX(:,:)
    CLOSE(20)

    ! Remove temporaty matrix file
    CALL SYSTEM_SUBR('rm temp.mat',ISTAT)

END SUBROUTINE

! This function returns the sum of the squared differences between the elements
! of the input matricies A and B
INTEGER FUNCTION SUMSQUAREDIFF(A,B,WIDTH)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: WIDTH,A(WIDTH,WIDTH),B(WIDTH,WIDTH)
    INTEGER :: I,J,AIJ,BIJ
    SUMSQUAREDIFF=0

    DO I=1,WIDTH
        DO J=1,WIDTH
            IF(I.EQ.J) CYCLE
            SUMSQUAREDIFF = SUMSQUAREDIFF + (A(I,J)-B(I,J))**2
        END DO
    END DO
    RETURN
END 

! This function returns the sum of the squared differences between the elements
! of the 1D input matricies A and B
INTEGER FUNCTION SUMSQUAREDIFF1D(A,B,LENGTH)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LENGTH,A(LENGTH),B(LENGTH)
    INTEGER :: I

    SUMSQUAREDIFF1D=0
    DO I=1,LENGTH
       SUMSQUAREDIFF1D = SUMSQUAREDIFF1D + (A(I)-B(I))**2
    END DO

    RETURN
END 

! This subroutine dumps an AMBER .rst format file (no velocities)  
SUBROUTINE DUMPRST(COORDS,NATOMS,FNAME)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NATOMS
    DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
    CHARACTER(LEN=*), INTENT(IN) :: FNAME

    ! DEBUG PRINTING
    !PRINT *,"In DUMPRST - NATOMS,FNAME:"
    !PRINT *,NATOMS,FNAME

    OPEN(UNIT=20,FILE=TRIM(ADJUSTL(FNAME))//TRIM(ADJUSTL('.rst')),STATUS='UNKNOWN')
    WRITE(20,'(a20)') TRIM(ADJUSTL(FNAME))
    WRITE(20,'(i5)') NATOMS
    WRITE(20,'(6f12.7)') COORDS(:)
    CLOSE(20)

END SUBROUTINE

! This subroutine dumps the hydrogen bond matrix passed in!
SUBROUTINE WRITEHBONDMATRIX(MATRIX,NRES,FNAME)

    IMPLICIT NONE

    INTEGER :: I,NRES
    INTEGER, INTENT(IN) :: MATRIX(NRES,NRES)
    CHARACTER(LEN=*), INTENT(IN) :: FNAME

    ! Write out the matrix
    OPEN(UNIT=20,FILE=TRIM(ADJUSTL(FNAME))//TRIM(ADJUSTL('.mat')),STATUS='UNKNOWN')
    DO I=1,NRES
        ! Use an upper bound to the number of columns to get around being unsure how
        ! wide it should be! Thanks to David for this :)
        WRITE(20,'(1000I3)') MATRIX(I,:)
    END DO
    CLOSE(20)

END SUBROUTINE

! This subroutine does the final output for the hydrogen-bond matrix run
SUBROUTINE HBONDMATRIXFINALIO()

    USE COMMONS, ONLY : NATOMS,HBONDBESTCOORDS,HBONDGROUPPOP,HBONDBEST,NHBONDGROUPS, &
    HBONDGROUPIDS,MCSTEPS,DUMPQUT,HBONDMAXE,HBONDQUE,HBONDQUZEROE,AMBERT,AMBER12T
    USE AMBER12_INTERFACE_MOD, ONLY : AMBER12_WRITE_PDB
    USE PORFUNCS

    IMPLICIT NONE

    INTEGER :: I,J,ISTAT
    CHARACTER(LEN=30) :: GROUPNUMBERCHAR,QUNUMBERCHAR
    CHARACTER(LEN=30) :: GROUPNAME,QUNAME
    DOUBLE PRECISION :: HBONDAVERAGEE(NHBONDGROUPS)

    ! Work out the average energy for each group
    ! Initialise the array and include Qu 0 in group 1
    HBONDAVERAGEE(:)=0.0D0
    HBONDAVERAGEE(1)=HBONDQUZEROE
    ! Loop over groups
    DO I=1,NHBONDGROUPS
       DO J=1,MCSTEPS(1)
       ! Was the structure at quench (step) J part of group I?
          IF (HBONDGROUPIDS(J).EQ.I) HBONDAVERAGEE(I)=HBONDAVERAGEE(I)+HBONDQUE(J)
       ENDDO
       HBONDAVERAGEE(I)=HBONDAVERAGEE(I)/HBONDGROUPPOP(I)
    ENDDO

    ! Write the key information for each group to a file
    OPEN(UNIT=20,FILE='hbondgroupinfo',STATUS='UNKNOWN')
    ! Write header for table
    WRITE(20,'(A)') ' ID  Size      E(min)              E(max)            E(average)'
    DO I=1,NHBONDGROUPS
        WRITE(20,'(I4,I4,3G20.10)') I,HBONDGROUPPOP(I),HBONDBEST(I),HBONDMAXE(I),HBONDAVERAGEE(I)
    ENDDO
    CLOSE(20)

    ! Dump the best (lowest E) structure for each group
    DO I=1,NHBONDGROUPS
        ! Construct the group minimum file name
        WRITE(GROUPNUMBERCHAR,*) I
        GROUPNAME='min.group'//TRIM(ADJUSTL(GROUPNUMBERCHAR))
        ! Dump .rst file
        CALL DUMPRST(HBONDBESTCOORDS(I,:),NATOMS,GROUPNAME)
        ! Dump .pdb file
        IF (AMBERT) THEN
           CALL A9DUMPPDB(HBONDBESTCOORDS(I,:),GROUPNAME)
        ELSEIF (AMBER12T) THEN
           CALL AMBER12_WRITE_PDB(HBONDBESTCOORDS(I,:),TRIM(ADJUSTL(GROUPNAME))//'.pdb', &
     &                            LEN(TRIM(ADJUSTL(GROUPNAME))//'.pdb'))
        ELSE
           PRINT *,"ERROR: You shouldn't be here...are you running HBONDMATRIX for a non AMBER potential?"
           STOP
        ENDIF
    ENDDO

    ! If DUMPQUT - produce cat'd .pdb file for each group and a reference file
    ! connecting each quench to the appropriate group
    IF (DUMPQUT) THEN
        DO I=1,NHBONDGROUPS
            ! Construct group .pdb file name
            WRITE(GROUPNUMBERCHAR,*) I
            GROUPNAME='allgroup'//TRIM(ADJUSTL(GROUPNUMBERCHAR))//'.pdb'
            DO J=1,MCSTEPS(1)
                ! Was the structure at quench (step) J part of group I?
                IF (HBONDGROUPIDS(J).EQ.I) THEN
                    ! Construct quench .pdb file name
                    WRITE(QUNUMBERCHAR,*) J
                    QUNAME='quench'// TRIM(ADJUSTL(QUNUMBERCHAR))//'.pdb'
                    ! Add a remark line so we know which quench the coordinates came from
                    CALL SYSTEM_SUBR('echo "REMARK ' // TRIM(ADJUSTL(QUNAME)) // '" >> ' &
                    // TRIM(ADJUSTL(GROUPNAME)),ISTAT)
                    ! cat the quench .pdb file onto the bottom of the group .pdb file
                    CALL SYSTEM_SUBR('cat ' // TRIM(ADJUSTL(QUNAME)) // ' >> ' &
                    // TRIM(ADJUSTL(GROUPNAME)),ISTAT)
                ENDIF
            ENDDO
        ENDDO
        ! Produce reference file connecting each quench to its group
        OPEN(UNIT=20,FILE='quenchtogroupref',STATUS='UNKNOWN')
        WRITE(20,'(A)') 'Quench  Group ID  Energy '
        WRITE(20,'(A,G20.10)') '     0     1',HBONDQUZEROE
        DO J=1,MCSTEPS(1)
            WRITE(20,'(2I6,G20.10)') J,HBONDGROUPIDS(J),HBONDQUE(J)
        ENDDO
        CLOSE(20)
    ENDIF
END SUBROUTINE

