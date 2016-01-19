      SUBROUTINE CHREBUILD(Q, BOND1, BOND2, THET1, THET2, PIC)
      use modamber9
      use commons, only : NATOMS
      IMPLICIT NONE

      DOUBLE PRECISION:: Q(3*NATOMS)
      DOUBLE PRECISION, INTENT(IN)::BOND1(*),BOND2(*),THET1(*),THET2(*),PIC(*)
      INTEGER :: ISEG, i1,I
      INTEGER :: at1, at2, at3, ires, lres
      DOUBLE PRECISION:: TEMP_AR(9)

      ires =1 
      DO ISEG =1, NSEG
         I1 =1
10       IF (ISEG .GT. 1) THEN
            IRES = ires +NICTOT(ISEG-1)  !first residue in current segment
         ENDIF

         CALL GETSEEDATM("N","CA","C",at1,at2,at3,IRES)

         IF (at1.LE.0 .OR. at2.LE.0. OR. at3.LE.0) THEN
             I1=I1+1
             GOTO 10
         ENDIF 

         TEMP_AR(1)= Q(3*(at1-1)+1);TEMP_AR(2)= Q(3*(at1-1)+2)
         TEMP_AR(3)= Q(3*(at1-1)+3);TEMP_AR(4)= Q(3*(at2-1)+1)
         TEMP_AR(5)= Q(3*(at2-1)+2);TEMP_AR(6)= Q(3*(at2-1)+3)
         TEMP_AR(7)= Q(3*(at3-1)+1);TEMP_AR(8)= Q(3*(at3-1)+2)
         TEMP_AR(9)= Q(3*(at3-1)+3)

         lres = IRES + NICTOT(ISEG) -1 !last residue in next segment

         Q(3*(ix(i02+IRES-1)-1)+1: 3*(ix(i02+lres)-1)) = 9999.0000

         DO i= 1,3
            Q(3*(at1-1)+i)= TEMP_AR(i)
            Q(3*(at2-1)+i)= TEMP_AR(3+i)
            Q(3*(at3-1)+i)= TEMP_AR(6+i)
         ENDDO

         !Note: to make this work for segments, need to change NSTART/NSTOP
         !ie make an array containing the amount of lenics per segment
         CALL CHARMMBILDC(1,lenic,Q,BOND1,BOND2,THET1,THET2,PIC)

      ENDDO    

      END SUBROUTINE CHREBUILD


C *************************************************************************




      SUBROUTINE CHARMMBILDC(NSTART,NSTOP,Q,
     4  BOND1,BOND2,THET1,THET2,PIC)
      use modamber9
      use commons, only : NATOMS 
      IMPLICIT NONE 
 
      INTEGER, INTENT(IN):: NSTART,NSTOP
      DOUBLE PRECISION Q(*)
      DOUBLE PRECISION BOND1(*),BOND2(*),THET1(*),THET2(*),PIC(*)
!      INTEGER  IC_COORDS(nres*95+5*nres,4) 
!      CHARACTER(LEN=3) IC_IMPROP(*)
!      INTEGER NATOM
C
      DOUBLE PRECISION R,THETA,PHI
      INTEGER INC,KK,II,I,J,K,L,ITEMP
      LOGICAL RETRY,JUNK,UNKN,CH_OK
      INTEGER CHANGE
      INTEGER OLDLEV
      DOUBLE PRECISION anum
      data anum /9999.00000000/
C
      IF(NSTOP.LT.NSTART) THEN
        WRITE(6,*) 'BILDC called with a null IC table'
        RETURN
      ENDIF
C
      INC=1
      KK=NSTART-1
C
      CHANGE=2
C msb50
C this loop runs forward and backward over all LENIC coordinates and calculates
C the internal coordinates of the last (L) of the 4 atoms per coordinate
C from the cartesian coordinates of the first three atoms (I,J,K)
C the loop is needed for the first terminus, as eg for N terminus the cartesian
C coordinates of the first three atoms are not known
C so the backwards loop tries to find the cartesian coordinates of I from
C the coordinates of I,J,K
C and then this goes on until nothing happens in two successive runs (change=2)
      DO WHILE(CHANGE.GT.0)
         CHANGE=CHANGE-1
         RETRY=.FALSE.
         DO II=NSTART,NSTOP
            KK=KK+INC
C
C Find atom whose coordinates are to be generated (atom l).
C
            I=IC_COORDS(KK,1)
            J=IC_COORDS(KK,2)
            K=IC_COORDS(KK,3)
            L=IC_COORDS(KK,4)
            PHI=(PIC(KK))
            R=(BOND2(KK))
            THETA=(THET2(KK))
C          PRINT '(2i3,a40,4i3,3f7.2)', ii,kk,"msb50 in bildc,i,j,k,l,phi,
C     &      r,theta",i,j,k,l,phi,r,theta
C
            IF(INC.LT.0) THEN
               R=(BOND1(KK))
               THETA=(THET1(KK))
               ITEMP=I
               I=L
               L=ITEMP
               IF(IC_IMPROP(KK) .EQ. "IMP") THEN
                  PHI=-PHI
               ELSE
                  ITEMP=J
                  J=K
                  K=ITEMP
               ENDIF
            ENDIF
C
C SEE IF COORDINATES ARE ALREADY KNOWN
C
            IF(I.LE.0) THEN
               CONTINUE
            ELSE IF(J.LE.0) THEN
               CONTINUE
            ELSE IF(K.LE.0) THEN
               CONTINUE
            ELSE IF(L.LE.0) THEN
               CONTINUE
            ELSE IF(R.LT.0.001) THEN
               CONTINUE
            ELSE IF(THETA.LT.0.001) THEN
               CONTINUE
            ELSE IF(I.EQ.K) THEN
               CONTINUE
            ELSE IF(I.EQ.J) THEN
               CONTINUE
            ELSE IF(K.EQ.J) THEN
               CONTINUE
            ELSE
C
C Check to see if all antecedents are known
C
               UNKN=Q(3*(L-1)+1).EQ.ANUM !X(L)
               JUNK=(Q(3*(I-1)+1).EQ.ANUM.OR.Q(3*(J-1)+1)
     &           .EQ.ANUM.OR.Q(3*(K-1)+1).EQ.ANUM)
               RETRY=RETRY.OR.JUNK.OR.UNKN
               IF (UNKN.AND..NOT.JUNK) THEN
C
C Set geometrical parameters
C
                  CALL CH_CARTCV(Q,I,J,K,L,R,THETA,PHI,CH_OK)
C                 PRINT '(i3,a10, 3f8.4)', ii,"l coords", X(L),Y(L),Z(L)
                  IF(CH_OK) CHANGE=2
               ENDIF
            ENDIF
         ENDDO
         KK=KK+INC
         INC=-INC
C
C Retry will be false if all coordinates are known
      ENDDO
C
      IF(RETRY) THEN
         WRITE(6,*) 'SOME COORDINATES NOT BUILT' 
      ENDIF
      KK=0
      DO I=1,NATOMS
         IF(X(I).EQ.ANUM) KK=KK+1
      ENDDO
      IF (KK.NE.0) WRITE(6,124) KK
  124 FORMAT(' ****  WARNING  ****',I5,
     1          ' COORDINATES ARE STILL UNDEFINED')
C
      RETURN
      !END
      END SUBROUTINE CHARMMBILDC

      SUBROUTINE CH_CARTCV(Q,I,J,K,L,RX,TX,PX,AT_OK)
C taken from CHARMM source/manip/intcor2.src
C
C     THIS ROUTINE FINDS THE POSITION OF L FROM THE POSITIONS
C     OF I,J,K AND R,THETA,PHI. THE COORDINATE ARRAY IS MODIFIED
C     ANGLES ARE IN DEGREES.
C
C     By Bernard R. Brooks    1982
C
      use modamber9
      use commons, only : NATOMS

      IMPLICIT NONE
      DOUBLE PRECISION Q(*), PI,RADDEG,DEGRAD !from source/fcm/consta.fcm
      PARAMETER(PI=3.141592653589793D0)
      PARAMETER (RADDEG=180.0D0/PI)
      PARAMETER (DEGRAD=PI/180.0D0)

      INTEGER I,J,K,L
      DOUBLE PRECISION RX,TX,PX
      LOGICAL AT_OK
C
      DOUBLE PRECISION R,T1,P,CST,SNT,CSP,SNP,XA,YA,ZA,XB,YB,ZB
      DOUBLE PRECISION RA,XC,YC,ZC,RC,WA,WB,WC
C
CRCZ 91/10/24 set trap for unknown atoms
      AT_OK=I.LE.0 .OR. J.LE.0 .OR. K.LE.0 .OR. L.LE.0 
      AT_OK=AT_OK .OR. I.GT.NATOMS .OR. J.GT.NATOMS .OR.
     $           K.GT.NATOMS .OR. L.GT.NATOMS
      IF(AT_OK) THEN 
           WRITE(6,'(A,4(1X,I5))')
     $        ' CARTCV> ERROR: Unknown atoms I,J,K,L=',I,J,K,L
           STOP !msb50
      ENDIF
CRCZ 91/10/24
      R=RX
      T1=DEGRAD*TX
      P=DEGRAD*PX
      CST=COS(T1)
      SNT=SIN(T1)
      CSP=COS(P)
      SNP=SIN(P)
C
      XA=Q(3*(J-1)+1)-Q(3*(K-1)+1)  !X(J) - X(K)
      YA=Q(3*(J-1)+2)-Q(3*(K-1)+2)  !Y(J)-Y(K)
      ZA=Q(3*(J-1)+3)-Q(3*(K-1)+3)  !Z(J)-Z(K)
      XB=Q(3*(I-1)+1)-Q(3*(J-1)+1)  !X(I)-X(J)
      YB=Q(3*(I-1)+2)-Q(3*(J-1)+2)  !Y(I)-Y(J)
      ZB=Q(3*(I-1)+3)-Q(3*(J-1)+3)  !Z(I)-Z(J)
C
      RA=SQRT(XA*XA+YA*YA+ZA*ZA)
      XA=XA/RA
      YA=YA/RA
      ZA=ZA/RA
C
      XC=YB*ZA-ZB*YA
      YC=ZB*XA-XB*ZA
      ZC=XB*YA-YB*XA
C
      RC=SQRT(XC*XC+YC*YC+ZC*ZC)
C
      IF (RC.LT.1.0D-09) THEN
         WRITE(6,22) I,J,K,L
  22     FORMAT(' Note from CARTCV: I-J-K is linear.',4I5,
     1       ', L not built.')
        AT_OK=.FALSE.
        RETURN
      ENDIF
C
      XC=XC/RC
      YC=YC/RC
      ZC=ZC/RC
C
      XB=YA*ZC-ZA*YC
      YB=ZA*XC-XA*ZC
      ZB=XA*YC-YA*XC
C
      WA=R*CST
      WB=R*SNT*CSP
      WC=R*SNT*SNP
C
      Q(3*(L-1)+1)=Q(3*(K-1)+1)+WA*XA+WB*XB+WC*XC
      Q(3*(L-1)+2)=Q(3*(K-1)+2)+WA*YA+WB*YB+WC*YC
      Q(3*(L-1)+3)=Q(3*(K-1)+3)+WA*ZA+WB*ZB+WC*ZC
C
      AT_OK=.TRUE.
      RETURN
      END SUBROUTINE CH_CARTCV



C**************************************************************


      SUBROUTINE CH_SEED(I,J,K,Q,BOND1,BOND2, ANGLE1,
     &   ANGLE2, PHI)
      use modamber9
      use commons, only : NATOMS
      IMPLICIT NONE
 
      DOUBLE PRECISION, INTENT(OUT) :: Q(*)
      DOUBLE PRECISION, INTENT(IN)  :: BOND1(*),BOND2(*),ANGLE1(*),ANGLE2(*),PHI(*)
!      INTEGER, INTENT(IN) :: IC_COORDS(nres*95+5*nres,4)
!      CHARACTER(LEN=3)    :: IC_IMPROP(nres*95+5*nres)
      INTEGER, INTENT(IN) :: I,J,K
      DOUBLE PRECISION              :: RIJ, RJK, THETA
      INTEGER             :: iic
      DOUBLE PRECISION              :: PI
      data PI /3.1415926535897931/

      RIJ = 0.0
      RJK = 0.0
      THETA = 0.0
      DO iic = 1,LENIC
         IF(BOND1(iic).GT.0.0) THEN
            IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,2).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,2).EQ.I.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,3).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,3).EQ.I.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,2).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,2).EQ.K.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,3).EQ.J.AND. 
     &          IC_IMPROP(iic).EQ."IMP")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,3).EQ.K.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RJK=BOND1(iic)
         ENDIF
         IF(BOND2(iic).GT.0.0) THEN
            IF(IC_COORDS(iic,4).EQ.I.AND.IC_COORDS(iic,3).EQ.J) RIJ=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.J.AND.IC_COORDS(iic,3).EQ.I) RIJ=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.K.AND.IC_COORDS(iic,3).EQ.J) RJK=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.J.AND.IC_COORDS(iic,3).EQ.K) RJK=BOND2(iic)
         ENDIF
C
         IF(IC_COORDS(iic,3).EQ.J) THEN
            IF(ANGLE2(iic).GT.0.0) THEN
               IF(IC_COORDS(iic,4).EQ.I.AND.IC_COORDS(iic,2).EQ.K) 
     &                 THETA=ANGLE2(iic)
               IF(IC_COORDS(iic,4).EQ.K.AND.IC_COORDS(iic,2).EQ.I) 
     &                 THETA=ANGLE2(iic)
            ENDIF
            IF(ANGLE1(iic).GT.0.0) THEN
               IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,2).EQ.K.AND.
     &                IC_IMPROP(iic).EQ."IMP")
     &                  THETA=ANGLE1(iic)
               IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,2).EQ.I.AND.
     &                IC_IMPROP(iic).EQ."IMP")
     &                  THETA=ANGLE1(iic)
            ENDIF
         ELSE
            IF(IC_COORDS(iic,2).EQ.J.AND.IC_IMPROP(iic).EQ."NOT"
     &           .AND.ANGLE1(iic).GT.0.0)THEN
               IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,3).EQ.K)
     &                THETA=ANGLE1(iic)
               IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,3).EQ.I) 
     &                THETA=ANGLE1(iic)
            ENDIF
         ENDIF
      ENDDO
C
      IF(RIJ.EQ.0.0.OR.RJK.EQ.0.0.OR.THETA.EQ.0.0) THEN
         PRINT*, "Error in SEED"
         PRINT '(A,3(1X,I5),3(1X,F7.2))',
     &    ' IC SEED> I,J,K,RIJ,RJK,THETA=',I,J,K,RIJ,RJK,THETA
         RETURN
      ENDIF
C
      Q(3*(I-1)+1) = 0.0
      Q(3*(I-1)+2) = 0.0
      Q(3*(I-1)+3) = 0.0 
      Q(3*(J-1)+1) = RIJ
      Q(3*(J-1)+2) = 0.0
      Q(3*(J-1)+3) = 0.0
      THETA=THETA*(PI/180.0)
      Q(3*(K-1)+1) = RIJ-RJK*COS(THETA)
      Q(3*(K-1)+2) = RJK*SIN(THETA)
      Q(3*(K-1)+3) = 0.0

      RETURN 
      END SUBROUTINE CH_SEED


C *******************************************************************************     
C from TAKESTEPCH in twist.src

      SUBROUTINE TAKESTEPAMM(Q)
      use modamber9
      use commons
C      USE KEY, ONLY: CHRIGIDT
      IMPLICIT NONE
      DOUBLE PRECISION    Q(3*NATOM)
      DOUBLE PRECISION    CHPMIN,CHPMAX
      INTEGER             ISEG,CHNMIN,CHNMAX,ISEED

      CHPMIN=0.2D0
      CHPMAX=0.4D0
      CHNMIN=1
      CHNMAX=0
      !msb50 correction for only twisting angles with difference > perthresh
      IF (AMBPERTT) THEN
            CHNMAX = NTW_ANGLES
            !PRINT*, "CHNMAX", CHNMAX
      ELSE
         DO ISEG=1,NSEG
            CHNMAX=CHNMAX+NPHIPSI(ISEG)+NSIDECHAIN(ISEG)
         ENDDO
      ENDIF
      CHNMAX=INT(CHNMAX/2.5D0) 
      ISEED=-1

      CALL PERTDIHAM(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
C      IF (CHRIGIDT .AND. NSEG.GT.1) THEN
C         CALL MKRIGIDTRANS(COORDS)
C         CALL MKRIGIDROT(COORDS)
C      ENDIF

      RETURN
      END SUBROUTINE TAKESTEPAMM

C  ***************************************************************************



C from PERTDIHE in twist.src

      SUBROUTINE PERTDIHAM(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
      use modamber9
      use commons
C      use KEY, ONLY: BHDEBUG, BHSTEPSIZE
      IMPLICIT NONE
      DOUBLE PRECISION::  Q(3*NATOMS), BHSTEPSIZE
      DOUBLE PRECISION::  P,ANGLE,DPRAND,MYRANDOM
      DOUBLE PRECISION::  CHPMIN,CHPMAX
      INTEGER::           ATOT,A,B,C,I1,J1,IICD, III
      INTEGER::           RESNUMseg,J3,ISEG,IRES,TOTPHIPSI,TOTSIDECHAIN
      INTEGER::           CHNMIN,CHNMAX,ISEED
!     DOUBLE PRECISION::    X(NATOMS),Y(NATOMS),Z(NATMS)
      LOGICAL::   TPP(NATOMS),TS(NATOMS),TO(NATOMS),TT,BHDEBUG
      INTEGER::   II,JJ,KK,LL, lenic2
      DOUBLE PRECISION :: BOND1(nbonh+nbona+nbper), BOND2(nbonh+nbona+nbper),
     1       THET1(ntheth+ntheta+ngper), THET2(ntheth+ntheta+ngper),
     2       PHI(nphia+nphih)
      DOUBLE PRECISION BB1,BB2,TT1,TT2,PP
      DOUBLE PRECISION :: ANGLE_SCALE(lenic)
      INTEGER:: count, P_RESPAR, P_RESMAX
      DOUBLE PRECISION:: ANGMAX,ANGMIN
      DOUBLE PRECISION :: PROBABIL
      DOUBLE PRECISION :: slope

      count=0
C      write(*,*)'CHPMIN,CHPMAX,CHNMIN,CHNMAX= ',CHPMIN,CHPMAX,CHNMIN,CHNMAX

      !msb50 - preliminary -  how many residues have a lin de/increasing probability
       P_RESPAR = 6 !number of residues for which twisting prob decreases
       P_RESMAX = 15
      !ANGMAX, ANGMIN for angles rescaling - twist more when at chain end
       ANGMAX= 1.2
       ANGMIN= 0.6

C  fill IC table with actual Cartesians
      CALL CHGETICVAL(Q,BOND1, BOND2, THET1, THET2, PHI, .FALSE.)
      
C initialise random nmber generator with input seed
      IF(ISEED.GT.-1) CALL SDPRND(ISEED)
C
C will be sent back to 192 if too many or too few dihedrals are altered
C as determined by CHNMIN and CHNMAX
192   CONTINUE

      count = count +1 !too avoid endless loop
C
C phi/psi angles, omega and sidechain dihedrals are stored in separate lists
C
C phi/psi angles
      B=0
      ATOT=0
      DO ISEG=1,NSEG !not at the moment - check twist.src for amendments
         DO A=1,NPHIPSI(ISEG)
            TPP(ATOT+A)=.FALSE.
C
C  Calculate P, the probability of twisting
C
            IF (AMBOLDPERTT) THEN
                IF (REAL(A).LE.(0.5*(NPHIPSI(ISEG)+1))) THEN
                   !P=CHPMAX-A*((CHPMAX-CHPMIN)/(NPHIPSI(ISEG)*0.5))
                   slope = (CHPMAX-CHPMIN)/(0.5*(NPHIPSI(ISEG)-1)-1)
                   P=-slope*A+(CHPMAX+slope)
                ELSE
                   !P=CHPMIN+(A-0.5*NPHIPSI(ISEG))*((CHPMAX-CHPMIN)/(NPHIPSI(ISEG)*0.5))
                   slope = (CHPMAX-CHPMIN)/(0.5*(NPHIPSI(ISEG)+1))
                   P = slope*A +(CHPMAX-slope*NPHIPSI(ISEG))
                END IF
            ELSE
            !msb50
            !take 2* P_RESMAX as NPHIPSI(ISEG) = 2*RESNUMseg(ISEG) 
            !  one phi, one psi angle per chain
                P = PROBABIL(A,NPHIPSI(ISEG), CHPMAX,CHPMIN,2*P_RESMAX,2*P_RESPAR)
               !msb50 - calculation of P depends on position of angle in chain only
               ! can't rely on TW_ANGLES being space out equally along the chain -->
               ! would not have P=CHNMIN in the middle of the chain if P= P(NTW_ANGLES)
                IF (AMBPERTT) THEN
                   !PRINT '(a15, 2f7.3)', "AMBERPERT probs",P,  TW_DIFFP(PHIPSI(A))
                   P = (P + TW_DIFFP(PHIPSI(A)))/2!rescaled prob+difference prob-then rescale
                ENDIF

                ! Calculate scaling factor for angle
                ANGLE_SCALE(PHIPSI(ATOT+A)) = PROBABIL(A,NPHIPSI(ISEG),
     &                        ANGMAX,ANGMIN,2*P_RESMAX,2*P_RESPAR)
            ENDIF
C            IF(BHDEBUG) WRITE(*,*)'P phipsi = ',P,ATOT+A
            MYRANDOM=DPRAND()
            !PRINT '(i4,a24,2f7.3)',PHIPSI(A), "nphipsi probabilities", MYRANDOM,P
            MYRANDOM=DPRAND()
            IF (MYRANDOM.LT.P) THEN
               TPP(ATOT+A)=.TRUE.
               B=B+1
            ENDIF
         ENDDO
         ATOT=ATOT+NPHIPSI(ISEG)
      ENDDO
      TOTPHIPSI=ATOT
C      IF(BHDEBUG) WRITE(*,*)'PHIPSI: TOT=',ATOT


C amide bond
      IF (TOMEGAC) THEN
         ATOT=0
         DO ISEG=1,NSEG
            DO A=1,NOMEGAC(ISEG)
               TO(ATOT+A)=.FALSE.
C
C  Calculate P, the probability of twisting
C
               IF (REAL(A).LE.(0.5*NOMEGAC(ISEG))) THEN
                  P=CHPMAX-A*((CHPMAX-CHPMIN)/(NOMEGAC(ISEG)*0.5))
               ELSE
                  P=CHPMIN+(A-0.5*NOMEGAC(ISEG))*((CHPMAX-CHPMIN)/(NOMEGAC(ISEG)*0.5))
               END IF
C               IF(BHDEBUG) WRITE(*,*)'P omega= ',P
               MYRANDOM=DPRAND()
               IF (MYRANDOM.LT.P) THEN
                  TO(ATOT+A)=.TRUE.
                  B=B+1
               ENDIF
            ENDDO
            ATOT=ATOT+NOMEGAC(ISEG)
         ENDDO
C         IF(BHDEBUG) WRITE(*,*)'OMEGA: TOT=',ATOT
      ENDIF

C  sidechains
      ATOT=0
      DO ISEG=1,NSEG
         DO A=1,NSIDECHAIN(ISEG)
            TS(ATOT+A)=.FALSE.
C
C  Calculate P, the probability of twisting
C  Make the probability dependent on the residue number
C  this means that all dihedrals in the same sidechain have the same P
C
            IICD=TW_SIDECHAIN(ATOT+A)
            j3 = IC_COORDS(IICD,2)
            DO III = 1,NRES
               IF (ix(i02+III-1).LE.j3 .AND.(ix(i02+III).GT.j3)) THEN
               ires = III
               EXIT
            ENDIF
            ENDDO
            IF (ISEG .GT. 1) THEN !NICTOT different from CHARMM!!
               RESNUMseg=NICTOT(ISEG)-NICTOT(ISEG-1)     ! number of residues in this segment
               IRES=IRES-NICTOT(ISEG-1) !no of residue in segment
            ELSE 
               RESNUMseg=NICTOT(ISEG)
            ENDIF
C            IF(BHDEBUG) WRITE(*,*)'IICD, JAR, IRES, RESNUM : ',IICD,JAR,IRES,RESNUM
            IF (AMBOLDPERTT) THEN
                IF (REAL(IRES).LE.(0.5*(RESNUMseg+1))) THEN
                   P=CHPMAX-IRES*(CHPMAX-CHPMIN)/(RESNUMseg*0.5)
                   slope = (CHPMAX-CHPMIN)/(0.5*(RESNUMseg-1)-1) !msb50
                   P=-slope*IRES+(CHPMAX+slope)
                ELSE
                   P=CHPMIN+(IRES-0.5*RESNUMseg)*(CHPMAX-CHPMIN)/(RESNUMseg*0.5)
                   slope = (CHPMAX-CHPMIN)/(0.5*(RESNUMseg+1))
                   P = slope*IRES +(CHPMAX-slope*RESNUMseg)
                END IF
            ELSE
                P = PROBABIL(IRES,RESNUMseg,CHPMAX,CHPMIN,P_RESMAX, P_RESPAR)
                IF (AMBPERTT) THEN
                   !PRINT '(a20, 2f7.3)', "AMBERPERT probs", P, TW_DIFFP(IICD)
                   P = (P + TW_DIFFP((IICD)))/2 !rescaled prob + difference prob - then rescale
                ENDIF
                !  Calculate the angle scaling factor
                ANGLE_SCALE(IICD) = PROBABIL(IRES,RESNUMseg,
     &                        ANGMAX,ANGMIN,P_RESMAX,P_RESPAR)
            ENDIF
            !PRINT '(i5,a20,i3,a,2f10.7)',IICD, "side probabs for", ires,":", MYRANDOM,P

            MYRANDOM=DPRAND()

C            IF(BHDEBUG) WRITE(*,*)'P sidechain =',P,ATOT+A
            MYRANDOM=DPRAND()
            IF (MYRANDOM.LT.P) THEN
               TS(ATOT+A)=.TRUE.
               B=B+1
            ENDIF
         ENDDO
         ATOT=ATOT+NSIDECHAIN(ISEG)
      ENDDO
      TOTSIDECHAIN=ATOT
C      IF(BHDEBUG) WRITE(*,*)'SIDECHAIN: TOT=',ATOT
C
C      shifting b dihedrals, should be NTEST1 < B < NTEST2
       IF (B.LT.CHNMIN .OR. B.GT.CHNMAX) THEN
         WRITE (*,'(A)') 'Too many dihedrals shifted - retrying'
         IF (count<11) THEN
            GOTO 192
         ELSEIF (11.LE.count .AND. count.LE.10000) THEN
            CHNMAX = CHNMAX +1 !msb50 - quite often tries to twist too many     
            !if AMBPERTONLY as CHNMAX = NTW_ANGLES << NPHIPSI+NSIDE
            GOTO 192
         ELSE
             PRINT*, "loop in perdiham failed"
             STOP
         ENDIF
       ENDIF

C  twisting phi/psi angles
       ATOT=0
       DO ISEG=1,NSEG
          DO A=1,NPHIPSI(ISEG)
             IF (TPP(ATOT+A)) THEN
                IICD=PHIPSI(ATOT+A)
C                IF(BHDEBUG) PRINT *,'PERTDIHE> changing phipsi ',IICD
                IF (AMBOLDPERTT) THEN
                    ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
                    !PRINT '(a20,i4,2f10.3)', "phipsi change", IICD,PHI(IICD),PHI(IICD)+ANGLE
                ELSE
                    ANGLE=((DPRAND()-0.5)*2.0*BHSTEPSIZE)*ANGLE_SCALE(IICD)
                    !PRINT '(i4,a6,f7.3,a6,f10.3)',IICD,"SCALE",ANGLE_SCALE(IICD),"ANGLE", ANGLE
                ENDIF
                PHI(IICD) = PHI(IICD) + ANGLE
             ENDIF
          ENDDO
          ATOT=ATOT+A
       ENDDO
C
C  twisting amide bond`
       IF (TOMEGAC) THEN
          ATOT=0
          DO ISEG=1,NSEG
             DO A=1,NOMEGAC(ISEG)
               IF (TO(ATOT+A)) THEN
                  IICD=OMEGAC(ATOT+A)
C                  IF(BHDEBUG) PRINT *,'PERTDIHE> changing omega ',IICD
                  ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
                  PHI(IICD) = PHI(IICD) +ANGLE
                  !PRINT '(a15,i4,2f10.3)', "omega change",IICD, PHI(IICD)-ANGLE,PHI(IICD)
               ENDIF
             ENDDO
             ATOT=ATOT+A
          ENDDO
       ENDIF
C
C  twisting sidechains
       ATOT=0
       DO ISEG=1,NSEG
          DO A=1,NSIDECHAIN(ISEG)
             IF (TS(ATOT+A)) THEN
                IICD=TW_SIDECHAIN(ATOT+A)
C                IF(BHDEBUG) PRINT *,'PERTDIHE> changing sidechain ',IICD
                IF (AMBOLDPERTT) THEN
                   ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
                   !PRINT '(a20,i4,2f10.3)', "side change", IICD,PHI(IICD),PHI(IICD)+ANGLE
                ELSE
                   ANGLE=((DPRAND()-0.5)*2.0*BHSTEPSIZE)*ANGLE_SCALE(IICD)
                   !PRINT '(i4,a6,f7.3,a6,f10.3)',IICD,"SCALE",ANGLE_SCALE(IICD),"ANGLE", ANGLE
                ENDIF
                PHI(IICD) = PHI(IICD) +ANGLE
             ENDIF
          ENDDO
          ATOT=ATOT+A
       ENDDO

       IF(PHI(IICD).LT.-180.0) PHI(IICD)=PHI(IICD)+360.0
       IF(PHI(IICD).GT.180.0) PHI(IICD)=PHI(IICD)-360.0
 
       CALL CHREBUILD(Q, BOND1, BOND2,THET1,
     &        THET2, PHI)
C      DO II =1, NATOMS
C         PRINT '(a4,4f11.5)',ih(m04+II-1),Q(3*(II-1)+1),
C     &          Q(3*(II-1)+2),Q(3*(II-1)+3)
C      ENDDO


      END SUBROUTINE PERTDIHAM

! **********************************************************************

      DOUBLE PRECISION FUNCTION PROBABIL(x, xmax, CHPMAX,CHPMIN, 
     & max_x_func, slope_par)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: CHPMAX, CHPMIN
      INTEGER, INTENT(IN) :: x, xmax, max_x_func, slope_par
      DOUBLE PRECISION :: P, slope
      DOUBLE PRECISION :: P_RESPAR

! msb50  - this calculates probability for twisting/scaling of how much angles are twisted
!    if xmax < max_x_func                        else:
!             (a certain no of residues)
!      P |                                         P |
!        |  x           x                            |  x                         x
!        |    x       x                              |    x                     x
!        |      x   x                                |      x                 x       
!        |        x                                  |        x x x x x x x x
!        |__________________                         |_____________________________
!           |     |     |                               |     |             |     |
!           1          xmax                             1                        xmax

      IF (xmax .GT. max_x_func) THEN
         slope = (CHPMAX-CHPMIN)/(slope_par)
         IF (x.LE.slope_par) THEN
            P = -slope*x+(CHPMAX+slope)
         ELSEIF ((slope_par.LT.x).AND.(x.LE.xmax-slope_par)) THEN
            P = CHPMIN
         ELSE 
            P = slope*x +(CHPMAX-slope*xmax)
         ENDIF
      ELSE
            slope = (CHPMAX-CHPMIN)/(0.5*(xmax-1))
         IF (REAL(x).LE.(0.5*(xmax+1))) THEN
            !P=CHPMAX-x*(CHPMAX-CHPMIN)/(xmax*0.5)
            P=-slope*x+(CHPMAX+slope)
         ELSE
            !P=CHPMIN+(x-0.5*xmax)*(CHPMAX-CHPMIN)/(xmax*0.5)
            P = slope*x +(CHPMAX-slope*xmax)
         ENDIF
      ENDIF

      PROBABIL = P
      RETURN
      END FUNCTION PROBABIL

! **********************************************************************


      SUBROUTINE SETDIHEAM()
      use modamber9
      use commons
      IMPLICIT NONE
      INTEGER IICD,ISLCT,OLDUSD,I3, J3,K3, L3,IRES,JRES,ISEG,A, NSEGLAST
      INTEGER AP,AO,AS,AC
      INTEGER i,j
      INTEGER SEG_START(NATOMS)
      LOGICAL LPHIPSI,LOMEGAC,LSIDECHAIN,LCHIRAL,LASTSIDECHAIN
      INTEGER NPHIPSITOT, NOMEGACTOT,NSIDECHAINTOT, NCHIRALTOT
      CHARACTER*4 TYPEI,TYPEJ,TYPEK
      CHARACTER*8 RESLAB

! msb50 - Warning: NICTOT defined differently from CHARMM!
!         CHARMM - residues per segment: NICTOT(ISEG+1) - NICTOT(ISEG) 
!              i.e. NICTOT(ISEG) is at which residue new segment starts
!         AMBER  - residues per segment: NICTOT(ISEG) - NICTOT(ISEG-1) 
!              i.e. NICTOT(ISEG) is at which residue segment finishes

      IF (.NOT. ALLOCATED(NICTOT)) ALLOCATE(NICTOT(NATOMS))

      IF (.NOT. ALLOCATED(IC_COORDS)) THEN
         CALL GETICCOORDS()
         PRINT*, "lenic", lenic
      ENDIF

      IF (.NOT. ALLOCATED(IS_SIDECHAIN)) THEN
          ALLOCATE(IS_SIDECHAIN(lenic))
          DO IICD=1, lenic
             i3 = IC_COORDS(IICD,1); j3 = IC_COORDS(IICD,2)
             k3 = IC_COORDS(IICD,3); l3 = IC_COORDS(IICD,4)
             CALL CHECK_SIDECHAIN(i3,j3,k3,l3,IICD,IS_SIDECHAIN)
          ENDDO
      ENDIF

      nseg = 0
      i = 1
      !SEG_START _ first atom of new segment
      SEG_START(1) =1
      DO WHILE (i .LE. NATOMS)
      IF (ih(m04+i-1).EQ.'OXT ') THEN !CTerm
         nseg = nseg +1
         DO j = 1, nres  !ix(i02+1) = no of atoms in 1st residue +1
            IF (ix(i02+j-1).LT.i .AND. ix(i02+j).GT.i) THEN
            NICTOT(nseg) = j !number of residues prior to segment + in segment
            EXIT
            ENDIF
         ENDDO

         i=i+1
         IF (i .LT. NATOMS) THEN
            IF (ih(m04+i-1).EQ.'N') THEN
               SEG_START(nseg+1) = i
            ELSEIF (ih(m04+i).EQ.'CH3 ') THEN !starts at i+1+1
               SEG_START(nseg+1) = i-1
            ELSE
               PRINT*, "setdiham: segment unknown start", ih(m04+i-1), ih(m04+i-1)
            ENDIF
         ENDIF

      ELSEIF (ih(m04+i-1).EQ.'CH3 ') THEN !C2Term
         nseg = nseg +1
         DO j = 1, nres  !ix(i02+1) 0 no of atoms in 1st residue +1
            IF (ix(i02+j-1).LT.i .AND. ix(i02+j).GT.i) THEN
            NICTOT(nseg) = j !number of residues per segment 
            EXIT
            ENDIF
         ENDDO
         
         i=i+4
         IF (i .LT. NATOMS) THEN
            IF (ih(m04+i-1).EQ.'N') THEN
               SEG_START(nseg+1) = i
            ELSEIF (ih(m04+i).EQ.'CH3 ') THEN !starts at i+1+1
               SEG_START(nseg+1) = i-1
            ELSE
               PRINT*, "setdiham: segment unknown start", ih(m04+i-1), ih(m04+i-1)
            ENDIF
         ENDIF

C     ELSEIF other terminals
      
      ELSE !no terminal
         i=i+1
      ENDIF   
      ENDDO

      IF (nseg .EQ. 0) THEN
         PRINT*, "in setdiham: no terminal found - unknown type"
         STOP
      ENDIF

       write(*,*)'NSEG=',NSEG
       IF (.NOT. ALLOCATED(NPHIPSI)) ALLOCATE(NPHIPSI(NSEG))
       IF (.NOT. ALLOCATED(NOMEGAC)) ALLOCATE(NOMEGAC(NSEG))
       IF (.NOT. ALLOCATED(NSIDECHAIN)) ALLOCATE(NSIDECHAIN(NSEG))
       IF (.NOT. ALLOCATED(NCHIRAL)) ALLOCATE(NCHIRAL(NSEG))
       IF (.NOT. ALLOCATED(PHIPSI)) ALLOCATE(PHIPSI(NATOMS))
       IF (.NOT. ALLOCATED(OMEGAC)) ALLOCATE(OMEGAC(NATOMS))
       IF (.NOT. ALLOCATED(TW_SIDECHAIN)) ALLOCATE(TW_SIDECHAIN(NATOMS))
       IF (.NOT. ALLOCATED(CHIRAL)) ALLOCATE(CHIRAL(NATOMS))

       NPHIPSI(1:NSEG)=0 
       NOMEGAC(1:NSEG)=0
       NSIDECHAIN(1:NSEG)=0
       NCHIRAL(1:NSEG)=0
       PHIPSI(1:NATOMS)=0
       OMEGAC(1:NATOMS)=0
       TW_SIDECHAIN(1:NATOMS)=0
       CHIRAL(1:NATOMS)=0

       AP=0
       AO=0
       AS=0
       AC=0

       NSEGLAST=1
       ISEG = 1
       DO IICD=1,LENIC
          LPHIPSI=.FALSE.
          LOMEGAC=.FALSE.
          LSIDECHAIN=.FALSE.
          LCHIRAL=.FALSE.
          i3 = IC_COORDS(IICD,1)
          DO i = 1, NSEG
             IF (i3 .GE. SEG_START(i)) THEN
             ISEG =i
             ENDIF
          ENDDO
          IF (ISEG .GT.NSEGLAST) THEN
             AP=AP+NPHIPSI(NSEGLAST)
             AO=AO+NOMEGAC(NSEGLAST)
             AS=AS+NSIDECHAIN(NSEGLAST)
             AC=AC+NCHIRAL(NSEGLAST)
             NSEGLAST=ISEG
          ENDIF
          CALL ICTYPECHECKAM(LPHIPSI,LOMEGAC,LSIDECHAIN, LCHIRAL,IICD)
          IF (LPHIPSI) THEN
             NPHIPSI(ISEG)=NPHIPSI(ISEG)+1
             PHIPSI(AP+NPHIPSI(ISEG))=IICD
C            print *,'LP: IICD=', IICD,AP+NPHIPSI(ISEG)
          ENDIF
          IF (LOMEGAC) THEN
             NOMEGAC(ISEG)=NOMEGAC(ISEG)+1
             OMEGAC(AO+NOMEGAC(ISEG))=IICD
C            print *,'LO: IICD=', IICD,AO+NOMEGAC(ISEG)
          ENDIF
          IF (LSIDECHAIN) THEN
             NSIDECHAIN(ISEG)=NSIDECHAIN(ISEG)+1
             TW_SIDECHAIN(AS+NSIDECHAIN(ISEG))=IICD
C            print *,'LS: IICD=', IICD, AS+NSIDECHAIN(ISEG)
          ENDIF
          IF (LCHIRAL) THEN
             NCHIRAL(ISEG)=NCHIRAL(ISEG)+1
             CHIRAL(AC+NCHIRAL(ISEG))=IICD
C            print *,'LC: IICD=', IICD, AC+NCHIRAL(ISEG)
          ENDIF
       ENDDO

       NPHIPSITOT=0
       NOMEGACTOT=0
       NSIDECHAINTOT=0
       NCHIRALTOT=0
       DO ISEG=1,NSEG
          NPHIPSITOT=NPHIPSITOT+NPHIPSI(ISEG)
          NOMEGACTOT=NOMEGACTOT+NOMEGAC(ISEG)
          NSIDECHAINTOT=NSIDECHAINTOT+NSIDECHAIN(ISEG)
          NCHIRALTOT=NCHIRALTOT+NCHIRAL(ISEG)
          WRITE(*,'(A,I4)')'number of internal coordinates for segment ',ISEG
          print *,'setdiheam NPHIPSI',NPHIPSI(ISEG)
          print *,'setdiheam NOMEGAC',NOMEGAC(ISEG)
          print *,'setdiheam NSIDECHAIN',NSIDECHAIN(ISEG)
          print *,'setdiheam NCHIRAL',NCHIRAL(ISEG)
       ENDDO
       WRITE(*,'(A)')'total number of internal coordinates'
       print *,'setdiheam> NPHIPSITOT= ',NPHIPSITOT
       print *,'setdiheam> NOMEGACTOT= ',NOMEGACTOT
       print *,'setdiheam> NSIDECHAINTOT= ',NSIDECHAINTOT
       print *,'setdiheam> NCHIRALTOT= ',NCHIRALTOT
C set NSEGATOMS, ie. number of atoms per segment
C       ALLOCATE(NSEGATOMS(NSEG))
C       DO ISEG=2,NSEG+1
C          NSEGATOMS(ISEG-1)=IBASE(NICTOT(ISEG)+1)-IBASE(NICTOT(ISEG-1)+1)
C       ENDDO
C       DO ISEG=1,NSEG
C          WRITE(*,'(A,I4,I6)')'SEG-NR,NSEGATOMS : ',ISEG,NSEGATOMS(ISEG)
C       ENDDO

      END SUBROUTINE SETDIHEAM



C ****************************************************************************     
C from SUBROUTINE ICCHECKPP in setdihe.src

      SUBROUTINE ICTYPECHECKAM(LPHIPSI,LOMEGAC,LSIDECHAIN, LCHIRAL,IICD)
      use modamber9
      use commons
      IMPLICIT NONE
      INTEGER :: IICD, RESNUM, KRESNUM
      INTEGER :: i
      INTEGER :: i3,j3,k3,l3
      CHARACTER*4:: TYPEJ,TYPEK,TYPEI,TYPEL
      CHARACTER*8::RESLAB,KRESLAB
      LOGICAL :: SIDECH(nphia+nphih)
!      CHARACTER*8:: SEGID(*),RESID(*),ih(m04+*),RES(*)
      LOGICAL LPHIPSI,LOMEGAC,LSIDECHAIN,LCHIRAL
      INTEGER SIGN_CHIRAL(NATOMS)
 
      i3 = IC_COORDS(IICD,1)
      j3 = IC_COORDS(IICD,2)
      k3 = IC_COORDS(IICD,3)
      l3 = IC_COORDS(IICD,4)
      TYPEI=ih(m04+i3-1)
      TYPEJ=ih(m04+j3-1)
      TYPEK=ih(m04+k3-1)
      TYPEL=ih(m04+l3-1)
     
C      PRINT*, IICD,"in ICTYPECHECKAM", i3,j3,k3,l3
      DO i=1, nres
         IF (ix(i02+i-1).LE.j3 .AND.(ix(i02+i).GT.j3)) THEN
            resnum = i
            EXIT
         ENDIF
      ENDDO
      DO i=1, nres
         IF (ix(i02+i-1).LE.k3 .AND.(ix(i02+i).GT.k3)) THEN
            kresnum = i
            EXIT
         ENDIF
      ENDDO
      reslab = ih(m02+resnum-1)
      kreslab = ih(m02+kresnum-1)

C ***** chiral ******
      IF (IC_IMPROP(IICD) .EQ. "IMP") THEN
         IF ((TYPEI.EQ.'N   ').AND.(TYPEJ.EQ.'C   ').AND.
     1        (TYPEK.EQ.'CA  ').AND.(TYPEL.EQ.'CB  ')) LCHIRAL=.TRUE.
         IF ((TYPEI.EQ.'OG1 ').AND.(TYPEJ.EQ.'CA  ').AND.
     1        (TYPEK.EQ.'CB  ').AND.(TYPEL.EQ.'CG2 ')) LCHIRAL=.TRUE.
         IF ((TYPEI.EQ.'CG1 ').AND.(TYPEJ.EQ.'CA  ').AND.
     1        (TYPEK.EQ.'CB  ').AND.(TYPEL.EQ.'CG2 ')) LCHIRAL=.TRUE.
         RETURN
C not interested in any impropers for twisting, not even if they are sidechains


      ELSE
         IF (.NOT. ALLOCATED(IS_SIDECHAIN)) THEN
             ALLOCATE(IS_SIDECHAIN(lenic))
             CALL CHECK_SIDECHAIN(i3,j3,k3,l3,IICD,IS_SIDECHAIN)
         ENDIF
         LSIDECHAIN = IS_SIDECHAIN(IICD)

      IF (.NOT. LSIDECHAIN) THEN

         IF (RESLAB.NE.'PRO')  THEN
         ! don't twist phi in proline
C ******* phi *******     
         ! the standard line works for c19 with ACE capping
            IF ((TYPEI.EQ.'C').AND.(TYPEJ.EQ.'N').AND.(TYPEK.EQ.'CA')
     &          .AND.(TYPEL.EQ.'C')) LPHIPSI=.TRUE.
         ! for first residue of c19 ic tables - no capping group
            IF ((TYPEI .EQ.'H1').AND.(TYPEJ.EQ.'N').AND.(TYPEK.EQ.'CA')
     &          .AND.(TYPEL.EQ.'C')) LPHIPSI=.TRUE.
         ENDIF
     
C ******** psi *********
         ! the standard line works for c19 with CBX capping
         IF ((TYPEI.EQ.'N').AND.(TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'C')
     1          .AND.(TYPEL.EQ.'N')) LPHIPSI=.TRUE.
         ! for last residue
         IF ((TYPEI.EQ.'N').AND.(TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'C')
     1          .AND.(TYPEL.EQ.'OXT')) LPHIPSI=.TRUE.


C ******** omega *********
         IF ((TYPEI.EQ.'CA').AND.(TYPEJ.EQ.'C').AND.(TYPEK.EQ.'N')
     1          .AND.(TYPEL.EQ.'CA')) LOMEGAC=.TRUE.
         ! different cappings - ACE, AMN, NIT

      ELSE 

C ******** backbone **********
         ! don't want to twist dihedrals in proline
         IF (RESLAB.EQ.'PRO'.OR.RESLAB.EQ.'NPRO'.OR.RESLAB.EQ.'CPRO') THEN
            LSIDECHAIN = .FALSE.
            RETURN
         ENDIF

! ********** sidechain **************************************8
         IF (RESLAB.EQ.'ARG'.OR.RESLAB.EQ.'NARG'.OR.RESLAB.EQ.'CARG')THEN
            IF ((TYPEJ.EQ.'NE').AND.(TYPEK.EQ.'CZ')) LSIDECHAIN=.FALSE.
            IF ((TYPEJ.EQ.'CZ').AND.(TYPEK.EQ.'NH1')) LSIDECHAIN=.FALSE.
            IF ((TYPEJ.EQ.'CZ').AND.(TYPEK.EQ.'NH2')) LSIDECHAIN=.FALSE.
         ELSEIF (RESLAB.EQ.'ASN'.OR.RESLAB.EQ.'NASN'.OR.RESLAB.EQ.'CASN') THEN
            IF ((TYPEJ.EQ.'CG').AND.(TYPEK.EQ.'ND2')) LSIDECHAIN=.FALSE.
            ! We do not twist the CG-ND2 bond since it is a peptide bond
         ELSEIF (reslab.EQ."HIE".OR.reslab.EQ.'NHIE'.OR.reslab.EQ.'CHIE') THEN
            LSIDECHAIN = .FALSE.
            IF ((TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'CB')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CB').AND.(TYPEK.EQ.'CG')) LSIDECHAIN=.TRUE.
            IF (.NOT.LSIDECHAIN)
     &       PRINT *, 'WARNING: His ring dihedrals will not be twisted.'
         ELSEIF (RESLAB.EQ.'GLN'.OR.RESLAB.EQ.'NGLN'.OR.RESLAB.EQ.'CGLN') THEN
            LSIDECHAIN = .FALSE.
            IF ((TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'CB')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CB').AND.(TYPEK.EQ.'CG')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CG').AND.(TYPEK.EQ.'CD')) LSIDECHAIN=.TRUE.
            ! We do not twist the CD-NE2 bond since it is a peptide bond
         ELSEIF (RESLAB.EQ.'PHE'.OR.RESLAB.EQ.'NPHE'.OR.RESLAB.EQ.'CPHE') THEN
            LSIDECHAIN = .FALSE.
            IF ((TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'CB')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CB').AND.(TYPEK.EQ.'CG')) LSIDECHAIN=.TRUE.
            ! None of the dihedrals in the phenylalanine ring are twisted
         ELSEIF (RESLAB.EQ.'TRP'.OR.RESLAB.EQ.'NTRP'.OR.RESLAB.EQ.'CTRP') THEN
            LSIDECHAIN = .FALSE.
            IF ((TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'CB')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CB').AND.(TYPEK.EQ.'CG')) LSIDECHAIN=.TRUE.
            ! Currently I am not twisting any of the dihedrals in the 
            ! tryptophan ring.
            ! I think this is correct, but I am not yet certain -- jdb
           PRINT*,'WARNING:None of the dihedrals in the tryptophan ring'
           PRINT*,'   will being twisted.'
         ELSEIF (RESLAB.EQ.'TYR'.OR.RESLAB.EQ.'NTYR'.OR.RESLAB.EQ.'CTYR') THEN
            LSIDECHAIN = .FALSE.
            IF ((TYPEJ.EQ.'CA').AND.(TYPEK.EQ.'CB')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CB').AND.(TYPEK.EQ.'CG')) LSIDECHAIN=.TRUE.
            IF ((TYPEJ.EQ.'CZ').AND.(TYPEK.EQ.'OH')) LSIDECHAIN=.TRUE.
            PRINT*,'WARNING: It is probably not of much consequence,'
            PRINT*,'but the oxygen in the ring,'
            PRINT*,'of tyrosine may not be twisted.'
         ENDIF

      ENDIF
      ENDIF

      IF (LPHIPSI) THEN
         IF (((TYPEI.EQ.'C').AND.(TYPEJ.EQ.'N').AND.
     1       (TYPEK.EQ.'CA').AND.(TYPEL.EQ.'C'))          .OR.
     &       ((TYPEJ.EQ.'N').AND.
     3       (TYPEK.EQ.'CA').AND.(TYPEL.EQ.'C'))) THEN
             !PRINT *,'PHI',IICD
         ENDIF
         IF (((TYPEI.EQ.'N').AND.(TYPEJ.EQ.'CA').AND.
     1      (TYPEK.EQ.'C').AND.(TYPEL.EQ.'N'))          .OR. 
     2       ((TYPEI.EQ.'N').AND.(TYPEJ.EQ.'CA').AND.
     3       (TYPEK.EQ.'C').AND.(TYPEL.EQ.'NT'))         .OR.
     &       ((TYPEI.EQ.'N').AND.(TYPEJ.EQ.'CA')
     5        .AND.(TYPEK.EQ.'C').AND.(TYPEL.EQ.'OXT'))) THEN
             !PRINT *,'PSI',IICD 
         ENDIF
      ENDIF

      END SUBROUTINE ICTYPECHECKAM


! *******************************************************************
      SUBROUTINE GET_TWISTABLE(ST_PHI,FI_PHI)
      use modamber9
      IMPLICIT NONE
C This subroutine is important if perturbation is done in internal coordinates and 
C if AMBPERTT is true
C calculates number of total twistable angles, NTW_ANGLES, sets elements which correspond
C to twistable coordinates in array TW_ANGLES to true and calculates P_DIFF, the 
C probability of twisting according to difference between start and end structure
      DOUBLE PRECISION, INTENT(IN)::ST_PHI(nphia+nphih), FI_PHI(nphia+nphih)
      INTEGER :: A_SIDE, A_PHIPSI, A_OMEGA
      INTEGER :: i
      DOUBLE PRECISION :: DIFF_PHI(nphia+nphih)
      DOUBLE PRECISION :: MAX_DIFF

      IF (AMBPERTT) THEN
         IF (.NOT. ALLOCATED(TW_DIFFP)) ALLOCATE(TW_DIFFP(lenic))
         TW_DIFFP = -0.2D0
         NTW_ANGLES = 0
         A_SIDE = 1
         A_PHIPSI = 1
         A_OMEGA = 1
         DO i=1, lenic
            DIFF_PHI(i) = FI_PHI(i)-ST_PHI(i)
            IF (DIFF_PHI(i) .GT.180.0D0) DIFF_PHI(i) =  360.0D0 -DIFF_PHI(i)
            IF (DIFF_PHI(i) .LT.-180.0D0) DIFF_PHI(I) = 360.0D0 +DIFF_PHI(i)
         ENDDO
         MAX_DIFF = MAXVAL(DABS(DIFF_PHI))
         DO i=1, lenic
            !take into account 2 pi circularity
            IF (PHIPSI(A_PHIPSI)==i) THEN
               IF (DABS(DIFF_PHI(i)).GT.PERTHRESH) THEN
                  TW_DIFFP(i) = (DABS(DIFF_PHI(i))/MAX_DIFF)/5.0D0 +0.20D0
                  !PRINT '(a15,i4,f7.3)',"twist phipsi",i, TW_DIFFP(i)
               ! to count how many twistable angles I have
                  NTW_ANGLES = NTW_ANGLES +1
               ENDIF
               A_PHIPSI = A_PHIPSI +1
            ELSEIF (TW_SIDECHAIN(A_SIDE)==i) THEN
               IF (DABS(DIFF_PHI(i)).GT.PERTHRESH) THEN
                  TW_DIFFP(i) =  (DABS(DIFF_PHI(i))/MAX_DIFF)/5.0D0 +0.20D0
                  !PRINT '(a15,i4,f7.3)',"twist side",i,TW_DIFFP(i)/5.0+0.2
                  NTW_ANGLES = NTW_ANGLES +1
               ENDIF
               A_SIDE = A_SIDE + 1
            ELSEIF (TOMEGAC) THEN
               IF (OMEGAC(A_OMEGA) ==i) THEN
                   IF (DABS(DIFF_PHI(i)).GT.PERTHRESH) THEN
                      TW_DIFFP(i) = DABS(DIFF_PHI(i))/MAX_DIFF
                      NTW_ANGLES = NTW_ANGLES +1
                   ENDIF
                   A_OMEGA = A_OMEGA +1
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      PRINT*, "NTW_ANGLES", NTW_ANGLES
      END SUBROUTINE GET_TWISTABLE
