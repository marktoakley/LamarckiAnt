SUBROUTINE MWFILM(COORDS, NATOMS, G, V, GTEST, BOXLX, BOXLY, BOXLZ, LAT)
! calculates energy and gradient for monoatomic water confined between hydrophobic plates
IMPLICIT NONE
INTEGER                       :: NATOMS, I, K
DOUBLE PRECISION              :: COORDS(3*NATOMS), BOXLX, BOXLY, BOXLZ, BOX(3), LAT(3,3), RI(3), ROUND
DOUBLE PRECISION              :: G(3*NATOMS), V
DOUBLE PRECISION              :: GMW(3*NATOMS), VMW
DOUBLE PRECISION              :: GFLM(3*NATOMS), VFLM
LOGICAL                       :: GTEST, PERIODIC(3)
INTEGER                       :: P, Q
DOUBLE PRECISION              :: AA, BB, A, CSTHTA, LMBDA, GMA, SGMAMW, EPSLNMW, SGMAFLM, EPSLNFLM

! MONOATOMIC WATER PARAMETERS
AA=7.049556277D0
BB=0.6022245584D0
P=4
Q=0
A=1.8D0
CSTHTA=-1.D0/3.D0
LMBDA=23.15D0
GMA=1.2D0
SGMAMW=2.3925D0 ! AA
EPSLNMW=25.89D0 ! kJ mol-1
! HYDROPHOBIC FILM PARAMETERS
SGMAFLM=3.56D0  ! AA
EPSLNFLM=0.569D0! kJ mol-1
! PERIODICITY PARAMETERS
PERIODIC=(/.TRUE., .TRUE., .FALSE./)
BOX=(/BOXLX, BOXLY, BOXLZ/)

VFLM=0.D0
GFLM=0.D0

IF (ANY(PERIODIC)) THEN
   DO I=1,NATOMS
      RI(:)=COORDS(I*3-2:I*3)
      DO K=1,3
         ! keep particles between plates
         IF (.NOT.(PERIODIC(K))) THEN
            IF (RI(K)<0.D0) RI(K)=BOX(K)*(1.D-10)
            IF (RI(K)>BOX(K)) RI(K)=BOX(K)*(1.D0-1.D-10)
         ENDIF
      ENDDO
      COORDS(I*3-2:I*3)=RI(:)
   ENDDO
ENDIF

CALL MW (COORDS/SGMAMW, NATOMS, VMW, GMW, GTEST, AA, BB, P, Q, A, CSTHTA, LMBDA, GMA, PERIODIC, BOX/SGMAMW, LAT/SGMAMW)
CALL AAAFILM (COORDS/SGMAFLM, NATOMS, VFLM, GFLM, GTEST, LAT(3,3)/SGMAFLM)
V=VMW*EPSLNMW+VFLM*EPSLNFLM
G=GMW*EPSLNMW/SGMAMW+GFLM*EPSLNFLM/SGMAFLM

END SUBROUTINE MWFILM

SUBROUTINE MW (COORDS, NATOMS, V, G, GTEST, AA, BB, P, Q, A, CSTHTA, LMBDA, GMA, PERIODIC, BOX, LAT)
! calculates energy and gradient for Stillinger-Weber type potentials in reduced units
! MW = monoatomic water 
IMPLICIT NONE
INTEGER                       :: NATOMS, I, J, K
DOUBLE PRECISION              :: COORDS(NATOMS*3), DMAT(NATOMS,NATOMS), RMAT(NATOMS,NATOMS,3), RI(3), RIJ(3), RET(10)
DOUBLE PRECISION, INTENT(OUT) :: G(NATOMS*3), V
LOGICAL                       :: GTEST, PERIODIC(3)
DOUBLE PRECISION              :: AA, BB, A, CSTHTA, LMBDA, GMA
INTEGER                       :: P, Q
DOUBLE PRECISION              :: BOX(3), LAT(3,3), ROUND

V=0.D0
G=0.D0
DMAT=0.D0
RMAT=0.D0

IF (ANY(PERIODIC)) THEN
! get minimum vectors and distances
   DO I=1,NATOMS
      DO J=I+1,NATOMS
         RIJ=COORDS(J*3-2:J*3)-COORDS(I*3-2:I*3)
         DO K=1,3
            IF (PERIODIC(K)) RIJ=RIJ-LAT(:,K)*ROUND(RIJ(K)/LAT(K,K))
         ENDDO
         DMAT(I,J)=SQRT(DOT_PRODUCT(RIJ,RIJ))
         DMAT(J,I)=DMAT(I,J)
         RMAT(I,J,:)=RIJ
         RMAT(J,I,:)=-1.D0*RIJ
      ENDDO
   ENDDO
ENDIF

! calculate MW energy and gradient
DO I=1,NATOMS
   DO J=1,NATOMS
      DO K=I,NATOMS
         IF ((K==I).AND.(I<J)) THEN !two-body
            CALL MWPHI2(DMAT(I,J), RET(1), RET(2:), GTEST, AA, BB, P, Q, A)
            V=V+RET(1)
            IF (GTEST) THEN
               G(I*3-2:I*3)=G(I*3-2:I*3)+RET(2)*RMAT(J,I,:)/DMAT(J,I)
               G(J*3-2:J*3)=G(J*3-2:J*3)+RET(2)*RMAT(I,J,:)/DMAT(I,J)
            ENDIF
         ELSEIF ((I/=J).AND.(J/=K).AND.(I/=K)) THEN !three-body
            CALL MWPHI3(RMAT(J,I,:), RMAT(J,K,:), DMAT(J,I), DMAT(J,K), &
                      RET(1), RET(2:), GTEST, A, LMBDA, CSTHTA, GMA)
            V=V+RET(1)
            IF (GTEST) THEN
               G(I*3-2:I*3)=G(I*3-2:I*3)+RET(2:4)
               G(J*3-2:J*3)=G(J*3-2:J*3)+RET(5:7)
               G(K*3-2:K*3)=G(K*3-2:K*3)+RET(8:10)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE MW

SUBROUTINE MWPHI2 (DIJ, V, G, GTEST, AA, BB, P, Q, A)
! calculates two-body contributions to the energy and gradient
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: DIJ
LOGICAL, INTENT(IN)           :: GTEST
DOUBLE PRECISION, INTENT(IN)  :: AA, BB, A
INTEGER, INTENT(IN)           :: P, Q
DOUBLE PRECISION, INTENT(OUT) :: V, G(9)
DOUBLE PRECISION              :: DAS, DIJ1

V=0.D0
G=0.D0
IF (DIJ>A) RETURN
DAS=1.D0/(DIJ-A)
DIJ1=1.D0/DIJ

V=AA*(BB*DIJ1**P-DIJ1**Q)*EXP(DAS)

IF (.NOT.GTEST) RETURN

G(1)=AA*DIJ1*(Q*DIJ1**Q-BB*P*DIJ1**P)*EXP(DAS)-V*DAS**2

RETURN
END SUBROUTINE MWPHI2

SUBROUTINE MWPHI3 (RJI, RJK, DJI, DJK, V, G, GTEST, A, LMBDA, CSTHTA, GMA)
! calculates three-body contributions to the energy and gradient
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: RJI(3), RJK(3), DJI, DJK
LOGICAL, INTENT(IN)           :: GTEST
DOUBLE PRECISION, INTENT(IN)  :: A, LMBDA, CSTHTA, GMA
DOUBLE PRECISION, INTENT(OUT) :: V, G(9)
DOUBLE PRECISION              :: CSTHTAIJK, DJIAS, DJKAS

V=0.D0
G=0.D0
IF ((DJI>A).OR.(DJK>A)) RETURN
CSTHTAIJK=DOT_PRODUCT(RJI,RJK)/(DJI*DJK)
DJIAS=1.D0/(DJI-A)
DJKAS=1.D0/(DJK-A)
V=LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA)**2

IF (.NOT.GTEST) RETURN
G(1:3)=-1.D0*V*GMA*RJI/DJI*DJIAS**2 &
       +2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *(RJK/(DJI*DJK)-CSTHTAIJK*RJI/DJI**2)

G(4:6)=+1.D0*V*GMA*(RJI/DJI*DJIAS**2+RJK/DJK*DJKAS**2) &
       -2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *((RJI+RJK)/(DJI*DJK)-CSTHTAIJK*(RJI/DJI**2+RJK/DJK**2))

G(7:9)=-1.D0*V*GMA*RJK/DJK*DJKAS**2 &
       +2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *(RJI/(DJI*DJK)-CSTHTAIJK*RJK/DJK**2)

RETURN
END SUBROUTINE MWPHI3

SUBROUTINE AAAFILM (COORDS, NATOMS, V, G, GTEST, HEIGHT)
! calculates component of energy and gradient corresponding to water-film interaction in reduced units
IMPLICIT NONE
INTEGER                       :: NATOMS
DOUBLE PRECISION              :: HEIGHT, COORDS(NATOMS*3)
LOGICAL                       :: GTEST
DOUBLE PRECISION              :: G(NATOMS*3), V
INTEGER                       :: I
DOUBLE PRECISION              :: RIZ1, RIZ2

V=0.D0
G=0.D0

DO I=1,NATOMS
   RIZ1=1.D0/MAX(COORDS(I*3),TINY(V))
   RIZ2=1.D0/MAX((HEIGHT-COORDS(I*3)),TINY(V))
   V=V+(((2.D0/15.D0)*RIZ1**9-RIZ1**3)+((2.D0/15.D0)*RIZ2**9-RIZ2**3))
   IF (GTEST) THEN
      G(I*3)=((3.D0*RIZ1**4-9.D0*(2.D0/15.D0)*RIZ1**10)-(3.D0*RIZ2**4-9.D0*(2.D0/15.D0)*RIZ2**10))
   ENDIF
ENDDO

RETURN
END SUBROUTINE AAAFILM

SUBROUTINE MWSTEP(COORDS,JP,NPAR,NATOMS,STEP,BOXLX,BOXLY,BOXLZ,LAT)
! takes a random step along each of the lattice vectors
IMPLICIT NONE
INTEGER :: NATOMS, I, JP, NPAR
DOUBLE PRECISION :: COORDS(NATOMS,NPAR), STEP, BOXLX, BOXLY, BOXLZ, DPRAND, LAT(3,3)

DO I=1,NATOMS
   COORDS(I*3-2:I*3,JP)=COORDS(I*3-2:I*3,JP)+LAT(:,1)/BOXLX*STEP*(DPRAND()-0.5D0)*2.D0
   COORDS(I*3-2:I*3,JP)=COORDS(I*3-2:I*3,JP)+LAT(:,2)/BOXLY*STEP*(DPRAND()-0.5D0)*2.D0
   COORDS(I*3-2:I*3,JP)=COORDS(I*3-2:I*3,JP)+LAT(:,3)/BOXLZ*STEP*(DPRAND()-0.5D0)*2.D0
ENDDO

END SUBROUTINE MWSTEP

SUBROUTINE MWDRAW(NATOMS, PERIODIC, LAT)
! writes an xyz file for viewing multiple cells in VMD
! writes cif files for each of NSAVE lowest energy structures
USE VEC3
USE COMMONS, ONLY: BOXLX, BOXLY, BOXLZ, PALPHA, PBETA, PGAMMA, NSAVE
USE QMODULE
IMPLICIT NONE
INTEGER :: NATOMS, J1, J2, K
DOUBLE PRECISION :: COORDS(NATOMS*3), LAT(3,3), LATINV(3,3), RI(3), ROUND, PI
LOGICAL :: PERIODIC
CHARACTER(LEN=99) :: FNAME

PI=4.D0*DATAN(1.D0)
PALPHA=180.*PALPHA/PI
PBETA=180.*PBETA/PI
PGAMMA=180.*PGAMMA/PI

CALL INVERT3x3(LAT(:,:),LATINV(:,:))
OPEN(UNIT=26, FILE='mW.xyz', STATUS='UNKNOWN')

DO J1=1,NSAVE
   WRITE(26,'(I6)')NATOMS
   IF (PERIODIC) THEN
      WRITE(26,'(A,I5,A,F15.8,A,I8,A,6F7.2,A)')'E',J1,'=',QMIN(J1),' ff',FF(J1),' pbc set {',&
&                                              BOXLX,BOXLY,BOXLZ,PALPHA,PBETA,PGAMMA,'}'
   ELSE
      WRITE(26,'(A,I5,A,F15.8,A,I8)')'E',J1,'=',QMIN(J1),' ff',FF(J1)
   ENDIF
   DO J2=1,QMINNATOMS(J1)
      RI(:)=QMINP(J1,J2*3-2:J2*3)
      IF (PERIODIC) THEN
          RI=MATMUL(LATINV(:,:),RI)
          DO K=1,3
             RI(K)=RI(K)-ROUND(RI(K))
             IF (RI(K)<0.D0) RI(K)=RI(K)+1.D0
          ENDDO
          RI=MATMUL(LAT(:,:),RI)
      ENDIF
      ! vmd uses upper triangular lattice matrix: switch a&c, alpha&gamma
      WRITE(26,'(A5,3F20.10)')'C',RI(1),RI(2),RI(3)
   ENDDO
ENDDO

CLOSE (UNIT=26)

IF (PERIODIC) THEN
   DO J1=1,NSAVE
      WRITE(UNIT=FNAME,FMT='("mW",I0,".cif")')J1
      OPEN(UNIT=26, FILE=FNAME, STATUS='UNKNOWN')
      WRITE(26,*)"_cell_length_a ",BOXLX
      WRITE(26,*)"_cell_length_b ",BOXLY
      WRITE(26,*)"_cell_length_c ",BOXLZ
      WRITE(26,*)"_cell_angle_alpha ",PALPHA
      WRITE(26,*)"_cell_angle_beta ",PBETA
      WRITE(26,*)"_cell_angle_gamma ",PGAMMA
      WRITE(26,*)"loop_"
      WRITE(26,*)"_atom_site_label"
      WRITE(26,*)"_atom_site_fract_x"
      WRITE(26,*)"_atom_site_fract_y"
      WRITE(26,*)"_atom_site_fract_z"
      DO J2=1,QMINNATOMS(J1)
         RI(:)=QMINP(J1,J2*3-2:J2*3)
         RI=MATMUL(LATINV(:,:),RI)
         DO K=1,3
            RI(K)=RI(K)-ROUND(RI(K))
            IF (RI(K)<0.D0) RI(K)=RI(K)+1.D0
         ENDDO
         WRITE(26,'(A5,3F20.10)')"C",RI(:)
      ENDDO
      CLOSE(26)
   ENDDO
ENDIF

END SUBROUTINE MWDRAW

FUNCTION ROUND(X)
! rounds to nearest integer (no universal intrinsic for this)
! no facility for edge cases (probably very unlikely)
IMPLICIT NONE
DOUBLE PRECISION :: X,ROUND
ROUND=SIGN(1.D0*INT(INT(2.D0*SIGN(X,1.D0)+1)/2),X)
END FUNCTION ROUND
