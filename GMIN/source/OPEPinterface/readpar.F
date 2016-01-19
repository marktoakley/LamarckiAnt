      SUBROUTINE RDPDB(nf,X,text2,Id_atom,text3,text4,numres)

      implicit none

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA

      character*100 line
      character*7 text2(*)
      character*5 text3(*)
      character*7 text4(*)
      integer     numres(*),Id_atom(*)
      double precision X(*)
      integer nf
      integer i

      i = 1
      do while (i < natom + 1)
        read(nf,'(a)') line
        if(line(1:4) .eq. 'ATOM' .or. line(1:6) .eq. 'HETATM') then
          read(line,5000) text2(i),Id_atom(i),text3(i),text4(i),
     $      numres(i), x(i*3-2),x(i*3-1),x(i*3)
          i = i + 1
        endif
      enddo
 5000 format(a,I4,a,a,I3,f12.3,2f8.3) 
      return
      end

C-----------------------------------------------------------------------
      SUBROUTINE RDTOP(nf,CHRG)

      parameter (MAXPRE = 1500)                !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)           !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)            !! maximum number of cart coord 
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)    !! max number of SC-SC interactions 
      parameter (MAXBO  = MAXNAT)             !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)            !! maximum number of bond angles 
      parameter (MAXPHI = MAXNAT*4)           !! maximum number of torsional angles
      parameter (MAXTTY = 50000)              !! 
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs 

      implicit double precision (a-h,o-z)

      CHARACTER*4 IGRAPH,LABRES

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA

      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno

      COMMON/MISC3/IPRES(MAXPRE),IGRAPH(MAXNAT),LABRES(MAXNAT)

      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC
      COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),
     $             PK(MAXPHI),PN(MAXPHI),
     $             PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),
     $             GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)

      COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),
     $             ICBH(MAXBO)
      COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),
     $             JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
      COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),
     1             ICP(MAXPHI)
      COMMON/ENER4/
     1 IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)

      COMMON/RUNHED/ITITL(20),ITITL1(20)

      COMMON/NONBON/NUMEX(MAXNAT),NATEX(MAXTTY)
      COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
      COMMON/NPAIR/NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)

      REAL*8 IDIEL
      DIMENSION NPAIR(MAXNAT,MAXNAT)
      double precision CHRG
      DIMENSION CHRG(MAXNAT)
      character*20 ctitl

c      OPEN(UNIT=12,FILE="ost_rdparm.top",status="unknown")


c      read(nf,9108) (ITITL(I),I=1,20)
      read(nf,*) ctitl
      read(nf,9118) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     $              NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     $          NUMBND,NUMANG,NPTRA,NATYP,NPHB,IDUM,IDUM,IDUM,IDUM,IDUM,
     $          IDUM,IDUM,IFBOX,NMXRS,IFCAP

      NTYPE = NTYPES*NTYPES
      NTTYP = (NTYPES*(NTYPES+1))/2

      natom3 = natom*3
c      print *, "natom3", natom3
c      print *, numbnd, numbnd

c      print*,' MAX NUMBER of RESIDUES ',MAXPRE 
c      print*,' MAX NUMBER of ATOMS ',MAXNAT 
c      print*,' MAX NUMBER of Cartesian Coordinates ',MAXXC 
c      print*,' MAX NUMBER of SC-SC interactions ',MAXPNB 
c      print*,' MAX NUMBER of bonds ',MAXBO 
c      print*,' MAX NUMBER of bond-angles ',MAXTH
c      print*,' MAX NUMBER of torsions ',MAXPHI
c      print*,' MAX NUMBER of nonbonded pairs ',MAXPAI

c
      IF(NATOM .GE. MAXNAT .or. (NBONH+MBONA) .GE. MAXBO .or.
     1 (NTHETH+MTHETA) .GE. MAXTH .or. (NPHIH+MPHIA) .GE. MAXPHI)THEN
      PRINT*,' ERRORS in the dimensions of the variables '
      STOP
      ENDIF 

C
C     ----- read THE SYMBOLS AND THE CHARGES AND THE MASSES -----
C
      ns = 6
      read(nf,9108) (igraph(i),I = 1,NATOM)
      read(nf,9128) (CHRG(I),I = 1,NATOM)
      read(nf,9128) (amass(i),I = 1,NATOM)
!      do i = 1, natom
!#ifdef H2D
!          amass(i) = 2.0d0*amass(i)
!#endif
!        write(*, *) i, amass(i)
!      end do
!       write(*,*) 'amass', amass(1)
      read(nf,9118) (iac(I),I = 1,NATOM)
      read(nf,9118) (numex(i),I = 1,NATOM)
      read(nf,9118) (nno(i),I = 1,NTYPE)
      read(nf,9108) (labres(i),I=1,NRES)
      read(nf,9118) (ipres(i),I=1,NRES)


C     ----- read THE PARAMETERS -----

      read(nf,9128) (RK(I),    I = 1,NUMBND)
      read(nf,9128) (REQ(I),   I = 1,NUMBND)
      read(nf,9128) (TK(I),    I = 1,NUMANG)
      read(nf,9128) (TEQ(I),   I = 1,NUMANG)
      read(nf,9128) (PK(I),    I = 1,NPTRA)
      read(nf,9128) (PN(I),    I = 1,NPTRA)
      read(nf,9128) (PHASE(I), I = 1,NPTRA)
      read(nf,9128) (SOLTY(I), I = 1,NATYP)
      read(nf,9128) (CN1(I),   I = 1,NTTYP)
      read(nf,9128) (CN2(I),   I = 1,NTTYP)

C
C     ----- read THE BONDING InfORMATIONS -----
C
      
      read(nf,9118) (ibh(I),jbh(I),icbh(I),I = 1,NBONH)
      read(nf,9118) (ib(I),jb(I),icb(i),I = 1,NBONA)

      read(nf,9118) (ith(I),jth(I),kth(i),icth(i), I=1, NTHETH)
      read(nf,9118) (it(I),jt(I),kt(i),ict(i), I=1, NTHETA)

      read(nf,9118) (iph(i),jph(i),kph(i),lph(i),icph(i), 
     1               I = 1,NPHIH)
      read(nf,9118) (ip(i),jp(i),kp(i),lp(i),icp(i),I=1,NPHIA) 

      read(nf,9118) (natex(i),I=1,NNB)


c      print *, ' END OF TOPOLOGY FILE ', NNB


C     ----- INVERSE OF THE MASS -----

       amass(1:natom) = 1.0d0/amass(1:natom)


C     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
!       write(*,*) 'dielc', dielc   ! OK
      IF (DIELC .GT. 1.0d0) THEN
           DUMD = SQRT(DIELC)
           DO 140 I = 1,NATOM
                CHRG(I) = CHRG(I)/DUMD
  140      CONTINUE
      ENDIF


C     --- Some parameters for vector ephi ---

      CALL DIHPAR(NPTRA,PK,PN,PHASE,GAMC,GAMS,IPN,FMN)


C     ---- DETERMINE THE PAIRLIST of NONBONDED INTERACTIONS
c      print*,' DETERMINE THE PAIRLIST of NONBONDED INTERACTIONS'
      npair(1:natom,1:natom) = 0
      DO I1=1,NATOM
      DO J1=I1+1,NATOM
      NPAIR(I1,J1) = 1 
      ENDDO
      ENDDO


      DO I1=1,NATOM
      DO J1=I1+1,NATOM
      DO K1=1,NBONH
      IF (I1 .eq. (IBH(K1)/3)+1 .AND. J1 .eq. (JBH(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
      ENDDO !! NBONH
      DO K1=1,NBONA
      IF (I1 .eq. (IB(K1)/3)+1 .AND. J1 .eq. (JB(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
      ENDDO !! NBONA

      DO K1=1,NTHETH
      IF (I1 .eq. (ITH(K1)/3)+1 .AND. J1 .eq. (KTH(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
      ENDDO !! NTHETH

      DO K1=1,NTHETA
      IF (TK(ICT(K1)) .ge. 6.1 .and. TK(ICT(K1)) .le. 6.20) THEN
      GO TO 100  !! THE CA-CA-CA ANGLE NOT SKIPPED
      ENDIF 
      IF (I1 .eq. (IT(K1)/3)+1 .AND. J1 .eq. (KT(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
 100  CONTINUE
      ENDDO !! NTHETA 

      DO K1=1,NPHIH
      IF (I1 .eq. (IPH(K1)/3)+1 .AND. J1 .eq. ABS(LPH(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
      ENDDO !! NPHIH 

      DO K1=1,NPHIA
      IF (PK(ICP(K1)) .ge. 3.05 .and. PK(ICP(K1)) .le. 3.15) THEN
      GO TO 200  !! THE CA-CA-CA ANGLE NOT SKIPPED
      ENDIF 
      IF (I1 .eq. (IP(K1)/3)+1 .AND. J1 .eq. ABS(LP(K1)/3)+1)THEN
      NPAIR(I1,J1) = 0
      ENDIF
 200  CONTINUE
      ENDDO !! NPHIA 
      ENDDO
      ENDDO

      npair2 = 0
      do i1=1,natom
      do j1=i1+1,natom
      if (npair(i1,j1) .eq. 1) then
      npair2 = npair2 + 1
      ipair(npair2) = i1
      jpair(npair2) = j1
      endif
      enddo
      enddo
c      print*,' TOTAL NUMBER of N-BONDS',npair2 
      
 9108 FORMAT(20A4)
 9118 FORMAT(12I6)
 9128 FORMAT(5E16.8)
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE DIHPAR(NUMPHI,PK,PN,PHASE,GAMC,GAMS,IPN,FMN)

      implicit double precision (a-h,o-z)
C
C     GET ADDITIONAL PARAMETERS FOR THE VECTOR EPHI 
C
      DIMENSION PK(*),PN(*),PHASE(*),GAMC(*),GAMS(*),IPN(*)
      DIMENSION FMN(*)

      DATA PI/3.141592653589793d+00/

      DO 100 I = 1,NUMPHI
        DUM = PHASE(I)
        IF(DABS(DUM-PI) .LE. 1.0d-03) DUM = SIGN(PI,DUM)

        DUMC = DCOS(DUM)
        DUMS = DSIN(DUM)
        IF(DABS(DUMC) .LE. 1.0d-06) DUMC = 0.0d0
        IF(DABS(DUMS) .LE. 1.0d-06) DUMS = 0.0d0

        GAMC(I) = DUMC*PK(I)
        GAMS(I) = DUMS*PK(I)

        FMN(I) = 1.0d0
        IF(PN(I) .LE. 0.0d0) FMN(I) = 0.0d0

        PN(I) = DABS(PN(I))
        IPN(I) = INT(PN(I)+1.0d-03)
  100 CONTINUE
      RETURN
      END
