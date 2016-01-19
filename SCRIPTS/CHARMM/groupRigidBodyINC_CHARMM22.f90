!
! This program constructs the list of groups for the generalised rigid body minimization 
! Halim Kusumaatmaja
! Modified for use with the all-atom CHARMM22 forcefield by Chris Whittleston
!!! MAKE SURE YOU UNCOMMENT THE READ LINE FOR YOUR INPUT FORMAT !!!

PROGRAM MAIN
  IMPLICIT NONE
  INTEGER, PARAMETER :: TNATOMS = 10000
  INTEGER ATOMNO(TNATOMS), RESIDUENO(TNATOMS), J1, J2, J3, J4, J5, J6, J7, IOSTATUS, TRESIDUE, TATOM, TGROUP
  INTEGER GROUPMEMBERS(TNATOMS, TNATOMS), GROUPCOLOUR(TNATOMS)
  INTEGER DUMMYINT
  CHARACTER(LEN=6) RECORDNAME
  CHARACTER(LEN=4) ATOMNAME(TNATOMS)
  CHARACTER(LEN=1) DUMMY1, DUMMY2
  CHARACTER(LEN=4) RESIDUENAME(TNATOMS)
  DOUBLE PRECISION X(TNATOMS), Y(TNATOMS), Z(TNATOMS)
  DOUBLE PRECISION DUMMYDP
  CHARACTER(LEN=80) FNAME, MYLINE
  CHARACTER(LEN=1) SCREENIN, DEFAULTAMP
  CHARACTER(LEN=3) RESIDUEDATABASE(30)
  CHARACTER(LEN=5) RESIDUEDATABASEMEMBER(30)
  LOGICAL LTEST, NONEFROMFILE, file_exists

  INTEGER ATOMCENTRE
  DOUBLE PRECISION :: RADIUSNONE, RADIUSALL, DIST
  INTEGER RESIDUENONE(TNATOMS), RESIDUEALL(TNATOMS)
  INTEGER JNONE, JALL

! Additional variables for outputting GROUPROTATION file
  INTEGER AXIS1, AXIS2
  CHARACTER(LEN=6) PROB(30), AMP(30), PROBA, AMPA
  DOUBLE PRECISION :: PROBD, AMPD
  INTEGER JDB, JUP, JSE
  INTEGER SELECTGRMODE, CURRENTRESIDUE

! Default amplitudes
  CHARACTER(LEN=6) AMP0, PROB0
  CHARACTER(LEN=6) ARGAMPAB, ARGAMPBG, ARGAMPGD
  CHARACTER(LEN=6) HISAMPAB, HISAMPBG
  CHARACTER(LEN=6) LYSAMPAB, LYSAMPBG, LYSAMPGD
  CHARACTER(LEN=6) ASPAMPAB, ASPAMPBG
  CHARACTER(LEN=6) GLUAMPAB, GLUAMPBG, GLUAMPGD
  CHARACTER(LEN=6) SERAMPAB, THRAMPAB
  CHARACTER(LEN=6) ASNAMPAB, ASNAMPBG
  CHARACTER(LEN=6) GLNAMPAB, GLNAMPBG, GLNAMPGD
  CHARACTER(LEN=6) CYSAMPAB, ALAAMPAB, VALAMPAB
  CHARACTER(LEN=6) ILEAMPAB, ILEAMPBG
  CHARACTER(LEN=6) LEUAMPAB, LEUAMPBG
  CHARACTER(LEN=6) METAMPAB, METAMPBG
  CHARACTER(LEN=6) PHEAMPAB, PHEAMPBG
  CHARACTER(LEN=6) TYRAMPAB, TYRAMPBG
  CHARACTER(LEN=6) TRPAMPAB, TRPAMPBG

  AMP0 = '0.1'
  PROB0 = '0.025'
  ARGAMPAB = '0.2'; ARGAMPBG = '0.3'; ARGAMPGD = '0.5'
  HISAMPAB = '0.3'; HISAMPBG = '0.5'
  LYSAMPAB = '0.2'; LYSAMPBG = '0.3'; LYSAMPGD = '0.5'
  ASPAMPAB = '0.5'; ASPAMPBG = '1.0'
  GLUAMPAB = '0.3'; GLUAMPBG = '0.5'; GLUAMPGD = '1.0'
  SERAMPAB = '1.0'; THRAMPAB = '1.0'
  ASNAMPAB = '0.5'; ASNAMPBG = '1.0'
  GLNAMPAB = '0.3'; GLNAMPBG = '0.5'; GLNAMPGD = '1.0'
  CYSAMPAB = '1.0'; ALAAMPAB = '1.0'; VALAMPAB = '1.0'
  ILEAMPAB = '0.5'; ILEAMPBG = '1.0'
  LEUAMPAB = '0.5'; LEUAMPBG = '1.0'
  METAMPAB = '0.5'; METAMPBG = '0.7'
  PHEAMPAB = '0.3'; PHEAMPBG = '0.5'
  TYRAMPAB = '0.3'; TYRAMPBG = '0.5'
  TRPAMPAB = '0.3'; TRPAMPBG = '0.4'
! Finish setting default amplitude  

  NONEFROMFILE=.FALSE.
  GROUPMEMBERS(1:TNATOMS,1:TNATOMS) = 0
  RESIDUEDATABASE(:) = 'OOO'

! -------------------------------------------------------------
! hk286 - input
! -------------------------------------------------------------

! Check that the input.crd file is present
INQUIRE(FILE='input.crd', EXIST=file_exists)
IF (file_exists) THEN
! If it is, proceed to read the comment lines and the number of atoms
   OPEN(UNIT=1,FILE='input.crd',STATUS='UNKNOWN')
   DO
      READ(1,*) MYLINE
      IF (MYLINE(1:1)=='*') THEN ! Check for the CHARMM comment line
         CYCLE
      ELSE
         READ(MYLINE,*) TATOM
         EXIT
      ENDIF
   ENDDO
! Now - read in the atom and residue info. 

!!! MAKE SURE YOU UNCOMMENT ONLY THE READ LINES BELOW WHICH MATCH YOUR INPUT FORMAT !!!

   DO J1=1,TATOM
! Example line from CHARMM input.crd:
!  244   16 ARG  HH11   5.77906  -0.47873  -6.23936 A    16     0.22450
!     READ(1,'(I5,I5,1X,A4,1X,A4,3F10.5)',IOSTAT=iostatus) &
!     ATOMNO(J1), RESIDUENO(J1), RESIDUENAME(J1), ATOMNAME(J1), X(J1), Y(J1), Z(J1)

! Example line from GMIN dbase.X format file (if used for input.crd):
!   12    1 SER  C           4.952733624192000        0.788045157454000       -1.540070962536000  A   1
      READ(1,*,IOSTAT=iostatus) &
      ATOMNO(J1), RESIDUENO(J1), RESIDUENAME(J1), ATOMNAME(J1), X(J1), Y(J1), Z(J1)

! Clean up whitespace in character strings
      ATOMNAME(J1)=ADJUSTL(ATOMNAME(J1))
      RESIDUENAME(J1)=ADJUSTL(RESIDUENAME(J1))
   ENDDO
   CLOSE(1)
! Set the number of residues
   TRESIDUE = RESIDUENO(TATOM)
ELSE
! If not, give a useful error and die
   PRINT *,"ERROR: input.crd not found in current directory"
   STOP
ENDIF  

! -------------------------------------------------------------
! hk286 - write coordsinirigid
! -------------------------------------------------------------

OPEN (UNIT = 28, FILE = 'coordsinirigid')
DO J1 = 1, TATOM
  WRITE (28, *), X(J1), Y(J1), Z(J1)
ENDDO
CLOSE (UNIT = 28)

! -------------------------------------------------------------
! hk286 - MAKE LIST OF RESIDUE BASED ON DISTANCE FROM ATOMCENTRE
! -------------------------------------------------------------

! Initialise counter variables
  JNONE = 0
  JALL = 0
  
  PRINT *, "ENTER THE ATOM INDEX WHERE IT IS THE CENTRE OF THE GROUPINGS"
  READ (5, *) ATOMCENTRE
  
! Read residues to be treated as all-atom from file if present
  INQUIRE(FILE='rigidbody.none', EXIST=file_exists)
  IF (file_exists) THEN
     PRINT *,'Reading in RESIDUENONE from file'
     NONEFROMFILE=.TRUE.
     OPEN(UNIT=223,FILE='rigidbody.none',STATUS='UNKNOWN')
     DO
        READ(223,*,IOSTAT=iostatus) DUMMYINT
        IF (iostatus<0) THEN
           CLOSE(223)
           EXIT
        ELSE
           JNONE = JNONE + 1
           RESIDUENONE(JNONE) = DUMMYINT
        ENDIF
     ENDDO
     CLOSE (UNIT=223)
! Otherwise these residues are identified using a radial cutoff
  ELSE
     PRINT *, "ENTER THE RADIUS BELOW WHICH WE DON'T WANT ANY GROUPINGS"
     READ (5, *) RADIUSNONE
  ENDIF
  
  PRINT *, "ENTER THE RADIUS ABOVE WHICH WE GROUP EVERYTHING INTO ONE"
  READ (5, *) RADIUSALL


! Loop over all atoms
  DO J1 = 1, TATOM
! Calculate the distance to the centre
     DIST = (X(J1)-X(ATOMCENTRE))**2 + (Y(J1)-Y(ATOMCENTRE))**2 + (Z(J1)-Z(ATOMCENTRE))**2
     DIST = SQRT(DIST)
! If the RESIDUENONE have not been read from a file above, read them in here
     IF (.NOT.NONEFROMFILE) THEN
! If the residue is within RADIUSNONE - treat it as all-atom
        IF (DIST < RADIUSNONE) THEN
           IF (J1 > 1) THEN
! If the atom belongs to the same residue as the previous atom - skip
              IF (RESIDUENO(J1) .EQ. RESIDUENO(J1-1)) THEN
              ELSE           
                 JNONE = JNONE + 1
! Add the residue number to the RESIDUENONE array
                 RESIDUENONE(JNONE) = RESIDUENO(J1)
              ENDIF
              ELSE
              JNONE = JNONE + 1
              RESIDUENONE(JNONE) = RESIDUENO(J1)           
           ENDIF
        ELSE
           IF ( J1 > 1 .AND. (RESIDUENO(J1) .EQ. RESIDUENO(J1-1)) .AND. (RESIDUENONE(JNONE) .EQ. RESIDUENO(J1))) THEN
              RESIDUENONE(JNONE) = 0
              JNONE = JNONE -1
           ENDIF
        ENDIF
     ENDIF
! If the residue lies outside RADIUSALL - treat it as part of the 'lump'
     IF (DIST > RADIUSALL) THEN
        IF (J1 > 1) THEN
           IF (RESIDUENO(J1) .EQ. RESIDUENO(J1-1)) THEN
           ELSE           
              JALL = JALL + 1
              RESIDUEALL(JALL)   = RESIDUENO(J1)
           ENDIF
        ELSE
           JALL = JALL + 1
           RESIDUEALL(JALL)   = RESIDUENO(J1)
        ENDIF
     ELSE
        IF ( J1 > 1 .AND. (RESIDUENO(J1) .EQ. RESIDUENO(J1-1)) .AND. (RESIDUENONE(JALL) .EQ. RESIDUENO(J1))) THEN
           RESIDUEALL(JALL) = 0
           JALL = JALL -1
        ENDIF
     ENDIF
  ENDDO
  
! DEBUG PRINTING
! PRINT *,'JNONE=',JNONE
! PRINT *,'RESIDUENONE(1:JNONE)='
! PRINT *,RESIDUENONE(1:JNONE)
! PRINT *,'JALL=',JALL
! PRINT *,'RESIDUEALL(1:JALL)='
! PRINT *,RESIDUEALL(1:JALL)

! Initialise group number counter (J2)
  J2 =1
! -------------------------------------------------------------
! hk286 - GROUPING PEPTIDE BOND
! -------------------------------------------------------------

  PRINT *, "Should I group PEPTIDE BOND C-O-N-H as a rigid body? (Y/N)"
  READ (5,*) SCREENIN

  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM        
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        IF (LTEST .EQ. .TRUE.) THEN
           IF (ATOMNAME(J1) .EQ. 'C') THEN
              GROUPMEMBERS(J2,1) = J1
           ELSEIF (ATOMNAME(J1) .EQ. 'O') THEN
              GROUPMEMBERS(J2,2) = J1
           ELSEIF ((ATOMNAME(J1) .EQ. 'N').AND.(.NOT.(RESIDUENAME(J1) .EQ. 'PRO')).AND.(RESIDUENO(J1)-1.EQ.RESIDUENO(GROUPMEMBERS(J2,2)))) THEN
              GROUPMEMBERS(J2,3) = J1
           ELSEIF ((ATOMNAME(J1) .EQ. 'HN').AND.(RESIDUENO(J1)-1.EQ.RESIDUENO(GROUPMEMBERS(J2,2)))) THEN
              GROUPMEMBERS(J2,4) = J1
              GROUPCOLOUR(J2) = 0
              J2 = J2+1
           ENDIF
        ENDIF
     ENDDO
     IF ((GROUPMEMBERS(J2,3).EQ.0).OR.(GROUPMEMBERS(J2,4).EQ.0)) THEN
        J2 = J2 - 1
     ENDIF
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING PROLINE
! -------------------------------------------------------------

  PRINT *, "Should I group PROLINE N-CD-CG-CB-CA as a rigid body? (Y/N)"
  READ (5,*) SCREENIN

  J3 = 1  
  IF (SCREENIN .EQ. 'Y') THEN
     
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        IF (LTEST .EQ. .TRUE.) THEN
           IF (RESIDUENAME(J1) .EQ. 'PRO') THEN
              IF (ATOMNAME(J1) .EQ. 'N') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CA') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CB') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CG') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 1
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     
  ENDIF

!! -------------------------------------------------------------
!! hk286 - GROUPING ARG
!! -------------------------------------------------------------

  PRINT *, "Should I group ARGININE NE-HE-CZ-NH1-HH11-HH12-NH2-HH21-HH22 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO

        IF ( LTEST .EQ. .TRUE. ) THEN
           IF ( RESIDUENAME(J1) .EQ. 'ARG' ) THEN
              IF (ATOMNAME(J1) .EQ. 'NE') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'NH1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH11') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH12') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'NH2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH21') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH22') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 2
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

!! -------------------------------------------------------------
!! hk286 - GROUPING HISTIDINE
!! -------------------------------------------------------------

  PRINT *, "Should I group HISTIDINE (HIS, HIE, HID) CG-ND1-CE1-NE2-CD2 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM
     
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO

        IF ( LTEST .EQ. .TRUE. ) THEN
           IF ((RESIDUENAME(J1) .EQ. 'HIS') .OR. (RESIDUENAME(J1) .EQ. 'HIE') .OR. (RESIDUENAME(J1) .EQ. 'HID')) THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'ND1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'NE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 3
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING LYSINE
! -------------------------------------------------------------

  PRINT *, "Should I group LYSINE NZ-HZ1-HZ2-HZ3 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM
     
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'LYS') THEN
              IF (ATOMNAME(J1) .EQ. 'NZ') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 4
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF


! -------------------------------------------------------------
! hk286 - GROUPING ASPARAGINE
! -------------------------------------------------------------

  PRINT *, "Should I group ASPARAGINE CB-CG-OD1-OD2 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN
     
     DO J1 = 1, TATOM

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN                
           IF (RESIDUENAME(J1) .EQ. 'ASP') THEN
              IF (ATOMNAME(J1) .EQ. 'CB') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CG') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'OD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'OD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 5
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING GLUTAMIC ACID
! -------------------------------------------------------------

  PRINT *, "Should I group GLUTAMIC ACID CG-CD-OE1-OE2 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN
     
     DO J1 = 1, TATOM

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN                
           IF (RESIDUENAME(J1) .EQ. 'GLU') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'OE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'OE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 6
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING ASPARAGINE
! -------------------------------------------------------------

  PRINT *, "Should I group ASPARAGINE CB-CG-OD1-ND2-HD21-HD22 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'ASN') THEN
              IF (ATOMNAME(J1) .EQ. 'CB') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CG') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'OD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'ND2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HD21') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HD22') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 7
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING GLUTAMINE
! -------------------------------------------------------------

  PRINT *, "Should I group GLUTAMINE CG-CD-OE1-NE2-HE21-HE22 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM     

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'GLN') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'OE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'NE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE21') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE22') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 8
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING PHENYLALANINE
! -------------------------------------------------------------

  PRINT *, "Should I group PHENYLALANINE CG-CD1-HD1-CE1-HE1-CZ-HZ-CE2-HE2-CD2-HD2 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM     

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'PHE') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 9
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING TYROSINE
! -------------------------------------------------------------

  PRINT *, "Should I group TYROSINE CG-CD1-HD1-CE1-HE1-CZ-OH-CE2-HE2-CD2-HD2 as a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     DO J1 = 1, TATOM     

        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'TYR') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'OH') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 10
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING TRYPTOPHAN
! -------------------------------------------------------------

  PRINT *, "Should I group TRYPTOPHAN as a 1 (CG-CD1-HD1-NE1-HE1-CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2) or 2 (CG-CD1-HD1-NE1-HE1)-(CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2) rigid bodies? (1/2/N)"
  READ (5,*) SCREENIN

  IF (SCREENIN .EQ. 'Y') THEN
     PRINT *,"Invalid option - please enter (1/2/N):"
     READ (5,*) SCREENIN
  ENDIF
   
  IF (SCREENIN .EQ. '1') THEN

     DO J1 = 1, TATOM     
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'TRP') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'NE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CE2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CH2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CE3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 11
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

  IF (SCREENIN .EQ. '2') THEN

     DO J1 = 1, TATOM     
        LTEST = .TRUE.
        DO J7 = 1, JNONE
           IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
! First group for TRP
        IF ( LTEST .EQ. .TRUE. ) THEN
           IF (RESIDUENAME(J1) .EQ. 'TRP') THEN
              IF (ATOMNAME(J1) .EQ. 'CG') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HD1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'NE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE1') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 11
                 J3 = 1
                 J2 = J2+1
              ENDIF
! Second group for TRP
              IF (ATOMNAME(J1) .EQ. 'CE2') THEN
                 J3 = 1
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CH2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CD2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1       
              ELSEIF (ATOMNAME(J1) .EQ. 'CZ3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HZ3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'CE3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HE3') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 J3 = J3+1
              ELSEIF (ATOMNAME(J1) .EQ. 'HH2') THEN
                 GROUPMEMBERS(J2,J3) = J1
                 GROUPCOLOUR(J2) = 1 
                 J3 = 1
                 J2 = J2+1
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING USER-DEFINED
! -------------------------------------------------------------
  PRINT *, "ANYTHING ELSE TO GROUP? (Y/N)"
  PRINT *, "WRITE IN USRRIGID FILE IN THE FOLLOWING FORMAT:"
  PRINT *, "TYR 6"
  PRINT *, "CG"
  PRINT *, "CD1"
  PRINT *, "CE1"
  PRINT *, "CZ"
  PRINT *, "CE2"
  PRINT *, "CD2"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     OPEN(UNIT=222,FILE='USRRIGID',status='old')
     DO
        READ(222,*,IOSTAT=iostatus) RESIDUEDATABASE(1), J5
        IF (iostatus<0) THEN
           CLOSE(222)
           EXIT
        ELSE 
           DO J4 = 1, J5
              READ(222,*) RESIDUEDATABASEMEMBER(J4)
           ENDDO
              
           DO J1 = 1, TATOM
              LTEST = .TRUE.
              DO J7 = 1, JNONE
                 IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                    LTEST = .FALSE.
                 ENDIF
              ENDDO
              DO J7 = 1, JALL
                 IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
                    LTEST = .FALSE.
                 ENDIF
              ENDDO
              
              IF ( LTEST .EQ. .TRUE. ) THEN
                 IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(1)) THEN
                    IF (ATOMNAME(J1) .EQ. TRIM(ADJUSTL(RESIDUEDATABASEMEMBER(1)))) THEN
                       J3 = 1
                       GROUPMEMBERS(J2,J3) = J1
                       J3 = J3+1
                    ELSEIF (ATOMNAME(J1) .EQ. TRIM(ADJUSTL(RESIDUEDATABASEMEMBER(J5)))) THEN
                       GROUPMEMBERS(J2,J3) = J1
                       GROUPCOLOUR(J2) = 12
                       J3 = 1
                       J2 = J2+1
                    ELSE
                       DO J6 = 2, J5-1
                          IF (ATOMNAME(J1) .EQ. TRIM(ADJUSTL(RESIDUEDATABASEMEMBER(J6)))) THEN
                             GROUPMEMBERS(J2,J3) = J1
                             J3 = J3+1
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     CLOSE(222)
  ENDIF

! -------------------------------------------------------------
! hk286 - GROUPING ALL OUTSIDE RADIUSALL AS ONE RIGID BODY
! -------------------------------------------------------------

  PRINT *, "Should I group all atoms outside RADIUSALL a rigid body? (Y/N)"
  READ (5,*) SCREENIN
  IF (SCREENIN .EQ. 'Y') THEN

     J3 = 0
     DO J1 = 1, TATOM     
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              J3 = J3 + 1
              GROUPMEMBERS(J2,J3) = J1
           ENDIF
        ENDDO
     ENDDO
     GROUPCOLOUR(J2) = 14                      
     J2 = J2 + 1
  ENDIF

! -------------------------------------------------------------
! hk286 - OUTPUT RIGID BODY GROUPINGS
! -------------------------------------------------------------
    
  TGROUP = J2-1

  OPEN(UNIT = 28, FILE = 'rbodyconfig')            
  J4 = 0
  DO J1 = 1, TGROUP
     J3 = 0
     DO J2 = 1, TATOM
        IF (GROUPMEMBERS(J1,J2) > 0) THEN
           J3 = J2
           J4 = J4 + 1
        ENDIF
     ENDDO
     WRITE (28,('(A6,I5)')), "GROUP ", J3
     DO J2 = 1, TATOM
        IF (GROUPMEMBERS(J1,J2) > 0) THEN
           WRITE (28, '(I5)'), GROUPMEMBERS(J1,J2) !ATOMNAME(GROUPMEMBERS(J1,J2)), RESIDUENAME(GROUPMEMBERS(J1,J2)), RESIDUENO(GROUPMEMBERS(J1,J2))
        ENDIF
     ENDDO
  ENDDO
  !PRINT *, J4
  CLOSE (UNIT=28)
  
  OPEN(UNIT=111,FILE='viewrbody.tcl',STATUS='UNKNOWN')

  IF (TGROUP < 1) THEN
     PRINT *, "NOTHING IS GROUPED! Check you uncommented the correct READ line for your input"
  ENDIF

  DO J1 = 1, TGROUP
     J3 = 0
     DO J2 = 1, TATOM
        IF (GROUPMEMBERS(J1,J2) > 0) THEN
           J3 = J2
        ENDIF
     ENDDO
     ! hk286
     WRITE(111,*) 'mol color ColorID ', GROUPCOLOUR(J1)
     WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
     WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J1,1:J3)-1
     WRITE(111,*) 'mol material Opaque'
     WRITE(111,*) 'mol addrep 0'
     ! hk286
  ENDDO
  CLOSE (UNIT = 111)
! -------------------------------------------------------------
! hk286 - AUTOMATIC GENERATION OF GROUPROTATION FILE
! -------------------------------------------------------------
  
  OPEN (UNIT = 28, FILE = 'atomgroups')
  
  PRINT *, "What is the uniform probability you want to choose? (0.01-0.1 seems sensible)"
  READ (5,*) PROB0

! -------------------------------------------------------------
! hk286 - ROTATION WITH N-CA as axis
! -------------------------------------------------------------

  RESIDUEDATABASE(:) = 'OOO'
  PROB(:) = '0.00D0'
  AMP(:) = '0.00D0'
  
  PRINT *, "Do you want group rotations of the side chain with N-CA as the axis? (Y/N)"
  PRINT *, "FYI I am not doing anything for PRO"
  READ (5,*) SCREENIN
  
  IF (SCREENIN .EQ. 'Y') THEN

     PRINT *, "Do you want to use the default values? (Y/N)"
     READ (5,*) DEFAULTAMP
     
     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
     
     OPEN(UNIT=111,FILE='viewNCA.tcl',STATUS='UNKNOWN')
     
     GROUPMEMBERS(:,:)=0
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        IF ( SELECTGRMODE .EQ. 2 ) THEN 
           DO J7 = 1, JNONE
              IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                 LTEST = .FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ((LTEST .EQ. .TRUE.).AND.(.NOT.(RESIDUENAME(J1).EQ.'PRO'))) THEN
           DO J2 = 2, TRESIDUE-1
              IF (RESIDUENO(J1).EQ.J2) THEN
                 IF (ATOMNAME(J1) .EQ. 'N')THEN
                    AXIS1 = J1
                    J3 = 1
                 ELSEIF (ATOMNAME(J1) .EQ. 'CA') THEN
                    AXIS2 = J1
                 ELSEIF (ATOMNAME(J1) .EQ. 'C') THEN
                    
                    JSE = 0
                    JUP = 0
                    DO JDB = 1, 30
                       IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(JDB)) THEN
                          JSE = JDB
                       ENDIF
                       IF (RESIDUEDATABASE(31-JDB) .EQ. 'OOO') THEN
                          JUP = 31 - JDB
                       ENDIF
                    ENDDO
                    
                    IF ( JSE.EQ.0 ) THEN
                       IF (  DEFAULTAMP .EQ. 'Y' ) THEN
                          PROB(JUP) = PROB0
                          AMP(JUP) = AMP0
                       ELSE
                          PRINT *, "What are the probabilities and amplitude of rotation? ", RESIDUENAME(J1), JUP
                          READ (5,*) PROB(JUP), AMP(JUP)
                       ENDIF
                       RESIDUEDATABASE(JUP) = RESIDUENAME(J1)
                       PROBA = PROB(JUP)
                       AMPA  = AMP(JUP)
                    ELSE
                       PROBA = PROB(JSE)
                       AMPA  = AMP(JSE)
                    ENDIF
                    
                    WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,A6,A1,A6)'), "GROUP ", RESIDUENO(J1), "SIDENC ", AXIS1, AXIS2, J3-1, " ", AMPA, " ", PROBA
                    DO J4 = 1, J3-1
                       WRITE (28, '(I5)'), GROUPMEMBERS(J2,J4)
                    ENDDO
                    
                    WRITE(111,*) 'mol color ColorID 4'
                    WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
                    WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J2,1:J3-1)-1
                    WRITE(111,*) 'mol material Opaque'
                    WRITE(111,*) 'mol addrep 0'
                    
                    J3 = 0
                 ELSEIF (.NOT.((ATOMNAME(J1) .EQ. 'N').OR.(ATOMNAME(J1) .EQ. 'H').OR.(ATOMNAME(J1) .EQ. 'C').OR.(ATOMNAME(J1) .EQ. 'O').OR.(ATOMNAME(J1) .EQ. 'CA'))) THEN ! .OR.(ATOMNAME(J1) .EQ. 'HA'))) THEN
                    IF (.NOT.(J3 .EQ. 0)) THEN
                       GROUPMEMBERS(J2, J3) = J1
                       J3 = J3+1
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        
     ENDDO
     
     CLOSE(UNIT = 111)
  ENDIF

! -------------------------------------------------------------
! hk286 - ROTATION WITH C-CA as axis
! -------------------------------------------------------------

  RESIDUEDATABASE(:) = 'OOO'
  PROB(:) = '0.00D0'
  AMP(:) = '0.00D0'
  
  PRINT *, "Do you want group rotations of the side chain with C-CA as the axis? (Y/N)"
  PRINT *, "FYI I am not doing anything for PRO"
  READ (5,*) SCREENIN
  
  IF (SCREENIN .EQ. 'Y') THEN

     PRINT *, "Do you want to use the default values? (Y/N)"
     READ (5,*) DEFAULTAMP

     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
     
     OPEN(UNIT=111,FILE='viewCCA.tcl',STATUS='UNKNOWN')
     
     GROUPMEMBERS(:,:)=0
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        IF ( SELECTGRMODE .EQ. 2 ) THEN 
           DO J7 = 1, JNONE
              IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                 LTEST = .FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
         
        IF ((LTEST .EQ. .TRUE.).AND.(.NOT.(RESIDUENAME(J1).EQ.'PRO'))) THEN
           DO J2 = 2, TRESIDUE-1
              IF (RESIDUENO(J1).EQ.J2) THEN
                 IF (ATOMNAME(J1) .EQ. 'CA')THEN
                    AXIS2 = J1
                    J3 = 1
                 ELSEIF (ATOMNAME(J1) .EQ. 'C') THEN
                    AXIS1 = J1
                    
                    JSE = 0
                    JUP = 0
                    DO JDB = 1, 30
                       IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(JDB)) THEN
                          JSE = JDB
                       ENDIF
                       IF (RESIDUEDATABASE(31-JDB) .EQ. 'OOO') THEN
                          JUP = 31 - JDB
                       ENDIF
                    ENDDO
                    
                    IF ( JSE.EQ.0 ) THEN
                       IF (  DEFAULTAMP .EQ. 'Y' ) THEN
                          PROB(JUP) = PROB0
                          AMP(JUP) = AMP0
                       ELSE
                          PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                          READ (5,*) PROB(JUP), AMP(JUP)
                       ENDIF
                       RESIDUEDATABASE(JUP) = RESIDUENAME(J1)
                       PROBA = PROB(JUP)
                       AMPA  = AMP(JUP)
                    ELSE
                       PROBA = PROB(JSE)
                       AMPA  = AMP(JSE)
                    ENDIF
                    
                    WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,A6,A1,A6)'), "GROUP ", RESIDUENO(J1), "SIDECC ", AXIS1, AXIS2, J3-1, " ", AMPA, " ", PROBA
                    DO J4 = 1, J3-1
                       WRITE (28, '(I5)'), GROUPMEMBERS(J2,J4)
                    ENDDO
                    
                    WRITE(111,*) 'mol color ColorID 6'
                    WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
                    WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J2,1:J3-1)-1
                    WRITE(111,*) 'mol material Opaque'
                    WRITE(111,*) 'mol addrep 0'
                    
                    J3 = 0
                 ELSEIF (.NOT.((ATOMNAME(J1) .EQ. 'N').OR.(ATOMNAME(J1) .EQ. 'H').OR.(ATOMNAME(J1) .EQ. 'C').OR.(ATOMNAME(J1) .EQ. 'O').OR.(ATOMNAME(J1) .EQ. 'CA'))) THEN !.OR.(ATOMNAME(J1) .EQ. 'HA'))) THEN
                    IF (.NOT.(J3 .EQ. 0)) THEN
                       GROUPMEMBERS(J2, J3) = J1
                       J3 = J3+1
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        
     ENDDO
     
     CLOSE(UNIT = 111)
     
  ENDIF

! -------------------------------------------------------------
! hk286 - ROTATION WITH CA-CB as axis
! -------------------------------------------------------------

  RESIDUEDATABASE(:) = 'OOO'
  PROB(:) = '0.00D0'
  AMP(:) = '0.00D0'
  
  PRINT *, "Do you want group rotations of the side chain with CA-CB as the axis? (Y/N)"
  PRINT *, "FYI I am not doing anything for ALA, GLY, PRO" ! SER, CYS ?
  READ (5,*) SCREENIN

  IF (SCREENIN .EQ. 'Y') THEN

     PRINT *, "Do you want to use the default values? (Y/N)"
     READ (5,*) DEFAULTAMP

     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
     
     OPEN(UNIT=111,FILE='viewCACB.tcl',STATUS='UNKNOWN')
     
     GROUPMEMBERS(:,:)=0
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        IF ( SELECTGRMODE .EQ. 2 ) THEN 
           DO J7 = 1, JNONE
              IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                 LTEST = .FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ((LTEST .EQ. .TRUE.).AND.(.NOT.((RESIDUENAME(J1).EQ.'ALA').OR.(RESIDUENAME(J1).EQ.'GLY').OR.(RESIDUENAME(J1).EQ.'PRO')))) THEN
!.OR.(RESIDUENAME(J1).EQ.'SER').OR.(RESIDUENAME(J1).EQ.'CYS')
           DO J2 = 2, TRESIDUE-1
              IF (RESIDUENO(J1).EQ.J2) THEN
                 IF (ATOMNAME(J1) .EQ. 'CA')THEN
                    AXIS1 = J1
                    J3 = 1
                 ELSEIF (ATOMNAME(J1) .EQ. 'CB') THEN
                    AXIS2 = J1
                 ELSEIF (ATOMNAME(J1) .EQ. 'C') THEN
                    
                    JSE = 0
                    JUP = 0
                    DO JDB = 1, 30
                       IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(JDB)) THEN
                          JSE = JDB
                       ENDIF
                       IF (RESIDUEDATABASE(31-JDB) .EQ. 'OOO') THEN
                          JUP = 31 - JDB
                       ENDIF
                    ENDDO
                    
                    IF ( JSE.EQ.0 ) THEN
                       IF (  DEFAULTAMP .EQ. 'Y' ) THEN
                          PROB(JUP) = PROB0
                          IF (RESIDUENAME(J1) .EQ. 'ARG') THEN                            
                             AMP(JUP) = ARGAMPAB
                          ELSE IF ( (RESIDUENAME(J1) .EQ. 'HIS') .OR. (RESIDUENAME(J1) .EQ. 'HID' ) .OR. (RESIDUENAME(J1) .EQ. 'HIE' ) ) THEN
                             AMP(JUP) = HISAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'LYS' )  THEN
                             AMP(JUP) = LYSAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ASP' )  THEN
                             AMP(JUP) = ASPAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLU' )  THEN
                             AMP(JUP) = GLUAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'SER' )  THEN
                             AMP(JUP) = SERAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'THR' )  THEN
                             AMP(JUP) = THRAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ASN' )  THEN
                             AMP(JUP) = ASNAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLN' )  THEN
                             AMP(JUP) = GLNAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'CYS' )  THEN
                             AMP(JUP) = CYSAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ALA' )  THEN
                             AMP(JUP) = ALAAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'VAL' )  THEN
                             AMP(JUP) = VALAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ILE' )  THEN
                             AMP(JUP) = ILEAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'LEU' )  THEN
                             AMP(JUP) = LEUAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'MET' )  THEN
                             AMP(JUP) = METAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'PHE' )  THEN
                             AMP(JUP) = PHEAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'TYR' )  THEN
                             AMP(JUP) = TYRAMPAB
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'TRP' )  THEN
                             AMP(JUP) = TRPAMPAB
                          ELSE
                             PRINT *, "Can't find any default value!"
                             PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                             READ (5,*) PROB(JUP), AMP(JUP)
                          ENDIF
                       ELSE
                          PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                          READ (5,*) PROB(JUP), AMP(JUP)
                       ENDIF
                       RESIDUEDATABASE(JUP) = RESIDUENAME(J1)
                       PROBA = PROB(JUP)
                       AMPA  = AMP(JUP)
                    ELSE
                       PROBA = PROB(JSE)
                       AMPA  = AMP(JSE)
                    ENDIF
                    
                    WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,A6,A1,A6)'), "GROUP ", RESIDUENO(J1), "SIDEAB ", AXIS1, AXIS2, J3-1, " ", AMPA, " ", PROBA
                    DO J4 = 1, J3-1
                       WRITE (28, '(I5)'), GROUPMEMBERS(J2,J4)
                    ENDDO
                    
                    WRITE(111,*) 'mol color ColorID 7'
                    WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
                    WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J2,1:J3-1)-1
                    WRITE(111,*) 'mol material Opaque'
                    WRITE(111,*) 'mol addrep 0'
                    
                    J3 = 0
                 ELSEIF (.NOT.((ATOMNAME(J1) .EQ. 'N').OR.(ATOMNAME(J1) .EQ. 'H').OR.(ATOMNAME(J1) .EQ. 'C').OR.(ATOMNAME(J1) .EQ. 'O').OR.(ATOMNAME(J1) .EQ. 'CA').OR.(ATOMNAME(J1) .EQ. 'HA').OR.(ATOMNAME(J1) .EQ. 'CB'))) THEN
                    IF (.NOT.(J3 .EQ. 0)) THEN
                       GROUPMEMBERS(J2, J3) = J1
                       J3 = J3+1
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        
     ENDDO
     
     CLOSE (UNIT = 111)
     
  ENDIF

! -------------------------------------------------------------
! hk286 - ROTATION WITH CB-CG as axis
! -------------------------------------------------------------
 
  RESIDUEDATABASE(:) = 'OOO'
  PROB(:) = '0.00D0'
  AMP(:) = '0.00D0'
  
  PRINT *, "Do you want group rotations of the side chain with CB-CG as the axis? (Y/N)"
  PRINT *, "FYI I am not doing anything for ALA, GLY, PRO, SER, CYS, THR, VAL"
  READ (5,*) SCREENIN

  IF (SCREENIN .EQ. 'Y') THEN
     
     PRINT *, "Do you want to use the default values? (Y/N)"
     READ (5,*) DEFAULTAMP

     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
     
     OPEN (UNIT = 111, FILE = 'viewCBCG.tcl', STATUS = 'UNKNOWN')
     
     GROUPMEMBERS(:,:)=0
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        IF ( SELECTGRMODE .EQ. 2 ) THEN 
           DO J7 = 1, JNONE
              IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                 LTEST = .FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ((LTEST .EQ. .TRUE.).AND.(.NOT.((RESIDUENAME(J1).EQ.'ALA').OR.(RESIDUENAME(J1).EQ.'GLY').OR.(RESIDUENAME(J1).EQ.'PRO') &
             .OR.(RESIDUENAME(J1).EQ.'SER').OR.(RESIDUENAME(J1).EQ.'CYS').OR.(RESIDUENAME(J1).EQ.'THR') .OR.(RESIDUENAME(J1).EQ.'VAL')))) THEN
           DO J2 = 2, TRESIDUE-1
              IF (RESIDUENO(J1).EQ.J2) THEN
                 IF (ATOMNAME(J1) .EQ. 'CB')THEN
                    AXIS1 = J1
                    J3 = 1
                 ELSEIF ( (ATOMNAME(J1) .EQ. 'CG') .OR. (ATOMNAME(J1) .EQ. 'CG1') )THEN
                    AXIS2 = J1
                 ELSEIF (ATOMNAME(J1) .EQ. 'C') THEN
                    
                    JSE = 0
                    JUP = 0
                    DO JDB = 1, 30
                       IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(JDB)) THEN
                          JSE = JDB
                       ENDIF
                       IF (RESIDUEDATABASE(31-JDB) .EQ. 'OOO') THEN
                          JUP = 31 - JDB
                       ENDIF
                    ENDDO
                    
                    IF ( JSE.EQ.0 ) THEN
                       IF (  DEFAULTAMP .EQ. 'Y' ) THEN
                          PROB(JUP) = PROB0
                          IF (RESIDUENAME(J1) .EQ. 'ARG') THEN                            
                             AMP(JUP) = ARGAMPBG
                          ELSE IF ( (RESIDUENAME(J1) .EQ. 'HIS') .OR. (RESIDUENAME(J1) .EQ. 'HID' ) .OR. (RESIDUENAME(J1) .EQ. 'HIE' ) ) THEN
                             AMP(JUP) = HISAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'LYS' )  THEN
                             AMP(JUP) = LYSAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ASP' )  THEN
                             AMP(JUP) = ASPAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLU' )  THEN
                             AMP(JUP) = GLUAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ASN' )  THEN
                             AMP(JUP) = ASNAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLN' )  THEN
                             AMP(JUP) = GLNAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'ILE' )  THEN
                             AMP(JUP) = ILEAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'LEU' )  THEN
                             AMP(JUP) = LEUAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'MET' )  THEN
                             AMP(JUP) = METAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'PHE' )  THEN
                             AMP(JUP) = PHEAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'TYR' )  THEN
                             AMP(JUP) = TYRAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'TRP' )  THEN
                             AMP(JUP) = TRPAMPBG
                          ELSE
                             PRINT *, "Can't find any default value!"
                             PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                             READ (5,*) PROB(JUP), AMP(JUP)
                          ENDIF
                       ELSE
                          PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                          READ (5,*) PROB(JUP), AMP(JUP)
                       ENDIF
                       RESIDUEDATABASE(JUP) = RESIDUENAME(J1)
                       PROBA = PROB(JUP)
                       AMPA  = AMP(JUP)
                    ELSE
                       PROBA = PROB(JSE)
                       AMPA  = AMP(JSE)
                    ENDIF
                    
                    WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,A6,A1,A6)'), "GROUP ", RESIDUENO(J1), "SIDEBG ", AXIS1, AXIS2, J3-1, " ", AMPA, " ", PROBA
                    DO J4 = 1, J3-1
                       WRITE (28, '(I5)'), GROUPMEMBERS(J2,J4)
                    ENDDO
                    
                    ! hk286
                    WRITE(111,*) 'mol color ColorID 0'
                    WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
                    WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J2,1:J3-1)-1
                    WRITE(111,*) 'mol material Opaque'
                    WRITE(111,*) 'mol addrep 0'
                    ! hk286
                    
                    J3 = 0
                 ELSEIF (.NOT.((ATOMNAME(J1) .EQ. 'N').OR.(ATOMNAME(J1) .EQ. 'H').OR.(ATOMNAME(J1) .EQ. 'C').OR.(ATOMNAME(J1) .EQ. 'O') &
                      .OR.(ATOMNAME(J1) .EQ. 'CA').OR.(ATOMNAME(J1) .EQ. 'CB').OR.(ATOMNAME(J1) .EQ. 'CG').OR.(ATOMNAME(J1) .EQ. 'CG1') &
                      .OR.(ATOMNAME(J1) .EQ. 'HA').OR.(ATOMNAME(J1) .EQ. 'HB').OR.(ATOMNAME(J1) .EQ. 'HB1').OR.(ATOMNAME(J1) .EQ. 'HB2') & 
                      .OR.(ATOMNAME(J1) .EQ. 'HB3') &
                      .OR.(ATOMNAME(J1) .EQ. 'CG2').OR.(ATOMNAME(J1) .EQ. 'HG21').OR.(ATOMNAME(J1) .EQ. 'HG22').OR.(ATOMNAME(J1) .EQ. 'HG23'))) THEN
                    IF (.NOT.(J3 .EQ. 0)) THEN
                       GROUPMEMBERS(J2, J3) = J1
                       J3 = J3+1
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        
     ENDDO
     
     CLOSE (UNIT = 111)
  ENDIF
  
! -------------------------------------------------------------
! hk286 - ROTATION WITH CG-CD as axis
! -------------------------------------------------------------

  RESIDUEDATABASE(:) = 'OOO'
  PROB(:) = '0.00D0'
  AMP(:) = '0.00D0'
  
  PRINT *, "Do you want group rotations of the side chain with CG-CD as the axis? (Y/N)"
  PRINT *, "FYI I am not doing anything for all residues except for ARG, LYS, GLU, GLN"
  READ (5,*) SCREENIN
  
  IF (SCREENIN .EQ. 'Y') THEN
     
     PRINT *, "Do you want to use the default values? (Y/N)"
     READ (5,*) DEFAULTAMP

     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
 
     OPEN(UNIT=111,FILE='viewCGCD.tcl',STATUS='UNKNOWN')
     
     GROUPMEMBERS(:,:)=0
     DO J1 = 1, TATOM
        
        LTEST = .TRUE.
        IF ( SELECTGRMODE .EQ. 2 ) THEN 
           DO J7 = 1, JNONE
              IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                 LTEST = .FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO J7 = 1, JALL
           IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
              LTEST = .FALSE.
           ENDIF
        ENDDO
        
        IF ((LTEST .EQ. .TRUE.).AND.((RESIDUENAME(J1).EQ.'ARG').OR.(RESIDUENAME(J1).EQ.'LYS').OR.(RESIDUENAME(J1).EQ.'GLU').OR.(RESIDUENAME(J1).EQ.'GLN'))) THEN
! .OR.(RESIDUENAME(J1).EQ.'GLU').OR.(RESIDUENAME(J1).EQ.'GLN')
           DO J2 = 2, TRESIDUE-1
              IF (RESIDUENO(J1).EQ.J2) THEN
                 IF (ATOMNAME(J1) .EQ. 'CG')THEN
                    AXIS1 = J1
                    J3 = 1
                 ELSEIF (ATOMNAME(J1) .EQ. 'CD') THEN
                    AXIS2 = J1
                 ELSEIF (ATOMNAME(J1) .EQ. 'C') THEN
                    
                    JSE = 0
                    JUP = 0
                    DO JDB = 1, 30
                       IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(JDB)) THEN
                          JSE = JDB
                       ENDIF
                       IF (RESIDUEDATABASE(31-JDB) .EQ. 'OOO') THEN
                          JUP = 31 - JDB
                       ENDIF
                    ENDDO
                    
                    IF ( JSE.EQ.0 ) THEN
                       IF (  DEFAULTAMP .EQ. 'Y' ) THEN
                          PROB(JUP) = PROB0
                          IF (RESIDUENAME(J1) .EQ. 'ARG') THEN                            
                             AMP(JUP) = ARGAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'LYS' )  THEN
                             AMP(JUP) = LYSAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLU' )  THEN
                             AMP(JUP) = GLUAMPBG
                          ELSE IF ( RESIDUENAME(J1) .EQ. 'GLN' )  THEN
                             AMP(JUP) = GLNAMPBG
                          ELSE
                             PRINT *, "Can't find any default value!"
                             PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                             READ (5,*) PROB(JUP), AMP(JUP)
                          ENDIF
                       ELSE
                          PRINT *, "What are the probabilities and amplitude of rotation?", RESIDUENAME(J1), JUP
                          READ (5,*) PROB(JUP), AMP(JUP)
                       ENDIF
                       RESIDUEDATABASE(JUP) = RESIDUENAME(J1)
                       PROBA = PROB(JUP)
                       AMPA  = AMP(JUP)
                    ELSE
                       PROBA = PROB(JSE)
                       AMPA  = AMP(JSE)
                    ENDIF
                    
                    WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,A6,A1,A6)'), "GROUP ", RESIDUENO(J1), "SIDEGD ", AXIS1, AXIS2, J3-1, " ", AMPA, " ", PROBA
                    DO J4 = 1, J3-1
                       WRITE (28, '(I5)'), GROUPMEMBERS(J2,J4)
                    ENDDO
                    
                    ! hk286
                    WRITE(111,*) 'mol color ColorID 13'
                    WRITE(111,*) 'mol representation CPK 1.000000 0.300000 10.000000 10.000000'
                    WRITE(111,'(A, 100(I))') 'mol selection index', GROUPMEMBERS(J2,1:J3-1)-1
                    WRITE(111,*) 'mol material Opaque'
                    WRITE(111,*) 'mol addrep 0'                                       
                    ! hk286
                    
                    J3 = 0
                 ELSEIF (.NOT.((ATOMNAME(J1) .EQ. 'N').OR.(ATOMNAME(J1) .EQ. 'H').OR.(ATOMNAME(J1) .EQ. 'C').OR.(ATOMNAME(J1) .EQ. 'O') &
                      .OR.(ATOMNAME(J1) .EQ. 'CA').OR.(ATOMNAME(J1) .EQ. 'CB').OR.(ATOMNAME(J1) .EQ. 'CG').OR.(ATOMNAME(J1) .EQ. 'CD') &
                      .OR.(ATOMNAME(J1) .EQ. 'HA').OR.(ATOMNAME(J1) .EQ. 'HB').OR.(ATOMNAME(J1) .EQ. 'HB1').OR.(ATOMNAME(J1) .EQ. 'HB2') & 
                      .OR.(ATOMNAME(J1) .EQ. 'HB3').OR.(ATOMNAME(J1) .EQ. 'HG1').OR.(ATOMNAME(J1) .EQ. 'HG2').OR.(ATOMNAME(J1) .EQ. 'HG3') &
                      .OR.(ATOMNAME(J1) .EQ. 'HG11').OR.(ATOMNAME(J1) .EQ. 'HG12').OR.(ATOMNAME(J1) .EQ. 'HG13') &
                      .OR.(ATOMNAME(J1) .EQ. 'HG21').OR.(ATOMNAME(J1) .EQ. 'HG22').OR.(ATOMNAME(J1) .EQ. 'HG23'))) THEN
                    IF (.NOT.(J3 .EQ. 0)) THEN
                       GROUPMEMBERS(J2, J3) = J1
                       J3 = J3+1
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO     
     CLOSE(UNIT = 111)     
  ENDIF

! -------------------------------------------------------------
! hk286 - USER-DEFINED
! -------------------------------------------------------------

  PRINT *, "ANYTHING ELSE? (Y/N)"
  PRINT *, "WRITE IN AUTOGR FILE IN THE FOLLOWING FORMAT:"
  PRINT *, "RESIDUENAME NO_MEMBER AMP PROB"
  PRINT *, "AXIS1"
  PRINT *, "AXIS2"
  PRINT *, "MEMBER1"
  PRINT *, "MEMBER2"
  PRINT *, "MEMBER3, ETC"
  READ (5,*) SCREENIN

  IF (SCREENIN .EQ. 'Y') THEN

     PRINT *, "For which regions do you want me to generate the GROUPROTATION FILE?"
     PRINT *, "(1) JUST THE ALL ATOM REGION"
     PRINT *, "(2) JUST THE MINIMAL RIGID BODY REGION"
     PRINT *, "(3) BOTH REGIONS"     
     READ (5,*) SELECTGRMODE
     
     OPEN(UNIT=222,FILE='AUTOGR',status='old')
     DO
        READ(222,*,IOSTAT=iostatus) RESIDUEDATABASE(1), J5, AMPD, PROBD
        IF (iostatus<0) THEN
           CLOSE(222)
           EXIT
        ELSE 
           DO J4 = 1, J5+2
              READ(222,*) RESIDUEDATABASEMEMBER(J4)
           ENDDO
           J2 = 1
           CURRENTRESIDUE = 0
           DO J1 = 1, TATOM
              LTEST = .TRUE.
              IF (SELECTGRMODE .EQ. 2) THEN
                 DO J7 = 1, JNONE
                    IF (RESIDUENO(J1) .EQ. RESIDUENONE(J7)) THEN
                       LTEST = .FALSE.
                    ENDIF
                 ENDDO
              ENDIF
              DO J7 = 1, JALL
                 IF (RESIDUENO(J1) .EQ. RESIDUEALL(J7)) THEN
                    LTEST = .FALSE.
                 ENDIF
              ENDDO
              
              IF ( LTEST .EQ. .TRUE. ) THEN
                 IF (RESIDUENAME(J1) .EQ. RESIDUEDATABASE(1)) THEN
                    DO J4 = 1, J5 + 2                       
                       IF ( (.NOT.(CURRENTRESIDUE .EQ. RESIDUENO(J1))) .AND. (J4 .EQ. 1) ) THEN
                          IF (CURRENTRESIDUE .EQ. 0) THEN
                             CURRENTRESIDUE = RESIDUENO(J1)
                          ELSE
                             WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,F8.4,A1,F8.4)'), "GROUP ", CURRENTRESIDUE, "USERSD ", GROUPMEMBERS(J2,1), GROUPMEMBERS(J2,2), J3, " ", AMPD, " ", PROBD
                             DO J6 = 3, J5 + 2
                                WRITE (28, '(I5)'), GROUPMEMBERS(J2,J6)
                             ENDDO                             
                             CURRENTRESIDUE = RESIDUENO(J1)                             
                          ENDIF
                       ENDIF
                       IF (ATOMNAME(J1) .EQ. TRIM(ADJUSTL(RESIDUEDATABASEMEMBER(J4)))) THEN
                          GROUPMEMBERS(J2,J4) = J1
                       ENDIF
                    ENDDO
                 ENDIF
              ENDIF
           ENDDO
           WRITE (28, '(A6,I5,A7,I5,I5,I5,A1,F8.4,A1,F8.4)'), "GROUP ", RESIDUENO(GROUPMEMBERS(J2,1)), "USERSD ", GROUPMEMBERS(J2,1), GROUPMEMBERS(J2,2), J3, " ", AMPD, " ", PROBD
           DO J6 = 3, J5 + 2
              WRITE (28, '(I5)'), GROUPMEMBERS(J2,J6)
           ENDDO           
        ENDIF
     ENDDO
     CLOSE(222)
  ENDIF

! -------------------------------------------------------------
! hk286 - CLOSING THE FILE
! -------------------------------------------------------------  
  CLOSE (UNIT = 28)
  
END PROGRAM MAIN


