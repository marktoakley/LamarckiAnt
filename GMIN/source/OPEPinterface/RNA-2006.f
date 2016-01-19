      module restraint_params
        integer :: Nrests,  Nposres
! Distance restraints
        integer, allocatable, dimension(:) :: resti, restj, trest
        double precision, allocatable, dimension(:) :: restk, restl,
     &   drestl

! Position restraints
        integer, allocatable, dimension(:) :: pri
        double precision, allocatable, dimension(:) :: prk
        double precision, allocatable, dimension(:,:) :: prx, prdx
        save
      endmodule restraint_params

      module RNAnb
        double precision, allocatable, dimension(:,:) :: nbcoef,nbct2
        double precision, allocatable, dimension(:) :: chrg
        integer, allocatable, dimension(:,:) :: nbscore
        integer bcoef(4,4)
        character basenum(4)
        save
      end module RNAnb

      module rnabase
        integer, allocatable, dimension(:) :: blist, btype
        save
      endmodule rnabase

      module RNAHBparams  
        double precision, dimension(6,4,4) :: dREF, alpam, alpbm, s
        integer :: Nparam(4,4)
        double precision :: wc = 1.3
        save
      endmodule RNAHBparams

      SUBROUTINE INITIALISE_RNA(P_B_C,BoxL,CofM,NATOM_CHECK,
     &                              XFORT,AMASSES,
     &                              ATOMIC_TYPE,
     &                              a_scaling_factor,
     &                              use_qbug_)
C
C   1.READ THE TOPOLOGY FILE (INTERNAL COORDINATES,
C     PARAMETERS OF THE OPEP FUNCTION) AND THE STARTING
C     PDB FILE: (UNIT 14).
C
C   2.READ THE PARAMETERS FOR THE SIDE-CHAIN SIDE-CHAIN INTERACTIONS
C     (UNIT 28).
C
C   3.COMPUTE ENERGY and FORCES (NO CUTOFF FOR THE NONBONDED-INTERACTIONS. 
C
C     THIS PROGRAM WAS WRITTEN BY P. DERREUMAUX while at IBPC, UPR 9080
C     CNRS, PARIS, FRANCE.  DECEMBER 1998
C
      use restraint_params
      use RNAnb
      use RNAHBparams
      use fileio
      use rnabase

      parameter (MAXPRE = 1500)    !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord 
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)!! max number of SC-SC interactions 
      parameter (MAXBO  = MAXNAT)  !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles 
      parameter (MAXPHI = MAXNAT*4)  !! maximum number of torsional angles
      parameter (MAXTTY = 50000)       !! maximum number of residue name types 
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs 

      implicit double precision (a-h,o-z)
      double precision :: xfort(3*NATOM_CHECK), amasses(NATOM_CHECK)
      character*5  ATOMIC_TYPE(NATOM_CHECK)

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC

      COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),
     $             ICBH(MAXBO)
      COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),
     $             JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
      COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),
     1             ICP(MAXPHI)
      COMMON/ENER4/
     1 IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)
      
      COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),
     $             PK(MAXPHI),PN(MAXPHI),
     $             PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),
     $             GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)
      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno

      common/scalingfactor/scaling_factor
      common/charge/cg,k1
!-------------------------------------------------------------------------------
      common/PBC_R/periodicBC,CM
!-------------------------------------------------------------------------------

C      common/texte/header,text
      common/texte/header


      COMMON/REWARD/ EHHB1
      COMMON/PROPE/ RESP(MAXPRE)
      COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
      common/pos/x_natv(maxxc)
c      common/forcek/force_const,rmsd
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom

      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)

      integer nfrag,lenfrag,ichain 

      common/scor/score(272),score_RNA(17)
      common/JANV06/NHB,ref_NHB,NHBH,NHBB

      logical qbug
      common/debug/qbug
      logical use_qbug_
C---------------------------------------------------------------------	  
      logical P_B_C,CofM
      real(8) BoxL
C---------------------------------------------------------------------	  

      real(8) IDIEL
      DIMENSION X(MAXXC),CG(MAXNAT)

      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)
      integer     numres(MAXNAT),Id_atom(MAXNAT)
      integer  ref_NHB

      CHARACTER*3 RESP 

      real(8) score_temp

      common/cacascsc/ct0lj(maxpnb),ct2lj(maxpnb)

      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $               rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $               rcut2_4b_out, rcut2_4b_in,
     $               rcut2_lj_out, rcut2_lj_in

      real(8) box_length, inv_box_length

      common/pbcBL/box_length, inv_box_length

      character(50) dummy

      integer i, j, idum

      logical fileexists

C________________________________________________________________________
      character*30 path_for_inits
      logical      single_conf_init, periodicBC,CM
      periodicBC = P_B_C
      CM = CofM
      write(*,*) 'PBC and CM are : ',P_B_C, CofM
C________________________________________________________________________

      qbug = use_qbug_

C --- SET THE PREFACTOR FOR SCALING THE POTENTIAL - ESSENTIAL FOR MD
      scaling_factor = a_scaling_factor !!3.2 for MD with alpha and
c      beta OPEP propensities 

C --- SET SOME PARAMETERS
      SCNB = 80.0      !! divide the 1-4 VDW interactions by 8.0
      SCEE = 80.0      !! divide the 1-4 ELEC interactions by 2.0
      DIELC = 1.0     !! 2.0 dielectric constant epsilon = 2r
      IDIEL = 0.0    !! ...............................
      CUT = 100.0      !! NO CUTOFF
      NHB = 13

C--   save data for optimization
C      OPEN(UNIT=56,FILE="proteinA-opt.dat",status="unknown")
c      OPEN(UNIT=56,FILE="Ab42-opt.dat",status="unknown")


C---  read TOPOLOGY file
      nf=21
C     This was the original configuration
      OPEN(UNIT=nf,FILE="parametres.top",status="unknown")
      CALL RDTOP(nf,CG)
      CLOSE(UNIT=nf)
 

C---  read INPUT PDB file
      CALL init_remd(path_for_inits,single_conf_init)
      if(.not. single_conf_init)  then 
         nf=41
         OPEN(UNIT=nf,FILE=path_for_inits,status="unknown")
         CALL RDPDB(nf,X,text2,Id_atom,text3,text4,numres) 
         close(nf)
      else  
         nf=14
         OPEN(UNIT=nf,FILE="conf_initiale_RNA.pdb",status="old")
         CALL RDPDB(nf,X,text2,Id_atom,text3,text4,numres) 
         close(nf)
      endif



c---  NEW JANV05
C---- Set-up the vector associating each atom with a fragment of the protein
C----    NFRAG : number of fragments
C----    length_fragment: table containing the length of each fragment
C----    ichain : the fragment to which atom idatom belongs to

C---  read ichain file
      nf=55
      OPEN(UNIT=nf,FILE="ichain.dat",status="old")
      read(55,*)  NFRAG 
      do i=1, NFRAG
         read(55,*) idum, lenfrag(i)
      enddo
      idatom = 1
      do i=1, NFRAG
         do j=1, lenfrag(i)
            ichain(idatom) = i
            idatom = idatom +1
         end do
      end do
      close(55)


C---  read the restraints
      inquire(file="restraints.dat", exist=fileexists)
      if(fileexists .eqv. .true.) then
c      if(testfile("restraints.dat") .eqv. .true.) then
        open(unit=55, file="restraints.dat", status="unknown")
        read(55,*) Nrests, Nposres
        allocate(resti(Nrests), restj(Nrests), trest(Nrests))
        allocate(restk(Nrests), restl(Nrests), drestl(Nrests))
        allocate(pri(Nposres), prk(Nposres))
        allocate(prx(3,Nposres), prdx(3,Nposres))
!       read a blank line for comments at the top of the list
        read(55, *) 
        read(55, *) (resti(i), restj(i), restk(i), restl(i),
     $   drestl(i), trest(i), i = 1,Nrests)
!       read a blank line for comments at the top of the posres list
        if(Nposres .gt. 0) then
          read(55, *) 
          read(55, *) (pri(i), prk(i), prx(1,i), prx(2,i), prx(3,i),
     $       prdx(1,i), prdx(2,i), prdx(3,i), i = 1,Nposres)
        endif
      endif

c        read(55,'(2i4, f15.10)') (resti(I), restj(i), restlens(i),
c     $      i = 1,Nrests)
      if( Nrests .gt. 0) then
        write (*,*) "Distance restraints"
        write (*,"(2a5,3a10,a5)") "ri","rj","k","length","dl","type"
        do i = 1, Nrests
          write(*,"(2i5, 3f10.6, i2)") resti(i), restj(i), restk(i), 
     &      restl(i), drestl(i), trest(i)
        enddo
      endif
      if(Nposres .gt. 0) then
        write (*,*) "Position restraints"
        write (*,"(2a5, 6a10)") "pi", "pk", "px", "py", "pz",
     &     "pdx", "pdy", "pdz"
        do i = 1, Nposres
          write(*,"(i5, 7f10.6)") pri(i), prk(i),
     &      prx(1,i), prx(2,i), prx(3,i),
     &      prdx(1,i), prdx(2,i), prdx(3,i)
        enddo
      endif


C---  read list of bases, for RNA
      if (.not. allocated(blist)) then
        allocate( blist(NRES))
        allocate( btype(NRES))
      endif
      open(unit=38,file="baselist.dat",status="old")
      read(38,*) (blist(i), btype(i), i =1, NRES)
      close(38)

C---  read weight file for RNA
      OPEN(UNIT=55,FILE="scale_RNA.dat",status="old")
      do i=1,17
         read(55,'(i4,f15.10,A)') idum, score_temp, dummy
         score_RNA(i)=dble(score_temp)
      enddo
      close(55)

      ! NTYPES is the number of particle types, found in the .top file
      ! Currently, it's taken to be :
      ! 1: C5*  2: O5*  3: P  4: CA  5: CY
      ! 6,7: G1,2  8,9: A1,2  10: U1  11: C1
      ! 12: D 13: MG
      if (.not. allocated(nbcoef)) then
        allocate( nbcoef(NTYPES, NTYPES))
        allocate( nbct2(NTYPES, NTYPES))
        allocate( nbscore(NTYPES, NTYPES))
        allocate( chrg(NTYPES))
      endif
      nbcoef = 1.0
      nbct2 = 3.6 ! 3.2
      nbscore = 4
      chrg = 0
!      ! bases beads size
      nbct2(6:11,6:11) = 3.2
      ! P-P interactions
      nbcoef(3,3) = 1.0
      nbct2(3,3) = 4.0
      nbscore(3,3) = 5
!     dummy residues are larger
      if (NTYPES .GT. 11) then
        nbct2(12,:) = 8
        nbct2(:,12) = 8
      endif
      ! Mg++ interactions
      nbscore(3,13) = 5
      nbscore(13,3) = 5
      nbscore(13,13) = 5
      chrg(3) = -1
      chrg(13) = 2

      do i = 1, NTYPES
        do j = 1, NTYPES
          nbcoef(i,j) = nbcoef(i,j)*score_RNA(nbscore(i,j))
        enddo
      enddo

      alpam = 2
      alpbm = 2
      bcoef = 35

      basenum = (/ 'G', 'A', 'C', 'U' /) 

      !                G   A   C  U
      bcoef(1,:) = (/ 34, 20, 19, 21 /)
      bcoef(2,:) = (/ 20, 32, 22, 18 /)
      bcoef(3,:) = (/ 19, 22, 33, 27 /)
      bcoef(4,:) = (/ 21, 18, 27, 35 /)


      ! Fill RNAHBparams

      !! dREF: distance reference, alpa/alpb: angles reference,
      !! s: strength of interaction
      
!     A-A
      dREF(1:6,2,2) =  (/ 5.43, 6.44, 6.86, 6.44, 5.86, 5.86/)
      alpam(1:6,2,2) = (/ 2.55, 1.01, 0.98, -1.73, 3.05, 2.33/)
      alpbm(1:6,2,2) = (/ 2.55, -1.73, 0.98, 1.01, 2.33, 3.05/)
      s(1:6,2,2) = (/ 2*wc, 1.0d0, 2.0d0, 1.0d0, 1*wc, 1*wc /)
      Nparam(2,2) = 6

!     A-C
      dREF(1:3,2,3) = (/ 4.8, 5.60, 6.94 /)
      alpam(1:3,2,3) = (/ 2.60, 2.40, 2.07 /)
      alpbm(1:3,2,3) = (/ 2.05, 1.82, 1.48 /)
      s(1:3,2,3) = (/ 1.0d0*wc, 1.0d0, 1.0d0 /)
      Nparam(2,3) = 3
      dREF(:,3,2) = dREF(:,2,3)
      alpam(:,3,2) = alpbm(:,2,3)
      alpbm(:,3,2) = alpam(:,2,3)
      s(:,3,2) = s(:,2,3)
      Nparam(3,2) = Nparam(2,3)

!     A-G
      dREF(1:4,2,1) = (/ 5.07, 6.26, 7.02, 7.54 /)
      alpam(1:4,2,1) = (/ 2.95, 1.20, 2.06, 0.80 /)
      alpbm(1:4,2,1) = (/ 2.67, -1.48, 1.93, -2.14 /)
      s(1:4,2,1) = (/ 2*wc, 2.0d0, 1.0d0, 1.0d0 /)
      Nparam(2,1) = 4
      dREF(:,1,2) = dREF(:,2,1)
      alpam(:,1,2) = alpbm(:,2,1)
      alpbm(:,1,2) = alpam(:,2,1)
      s(:,1,2) = s(:,2,1)
      Nparam(1,2) = Nparam(2,1)

!     A-U
      dREF(1:2,2,4) = (/ 4.96, 6.43 /)
      alpam(1:2,2,4) = (/ 2.92, 0.84 /)
      alpbm(1:2,2,4) = (/ 2.23, 1.95 /)
      s(1:2,2,4) = (/ 2*wc, 2.0d0 /)
      Nparam(2,4) = 2
      dREF(:,4,2) = dREF(:,2,4)
      alpam(:,4,2) = alpbm(:,2,4)
      alpbm(:,4,2) = alpam(:,2,4)
      s(:,4,2) = s(:,2,4)
      Nparam(4,2) = Nparam(2,4)

!     C-C
      dREF(1:1,3,3) =  (/ 4.91 /)
      alpam(1:1,3,3) = (/ 2.22 /)
      alpbm(1:1,3,3) = (/ 2.24 /)
      s(1:1,3,3) = (/ 2.0d0/)
      Nparam(3,3) = 1

!     C-G
      dREF(1:2,3,1) = (/ 4.80, 7.40 /)
      alpam(1:2,3,1) = (/ 2.17, 2.89 /)
!      alpbm(1:2,3,1) = (/ 2.76, 1.28 /)
      alpbm(1:2,3,1) = (/ 2.53, 1.28 /)
      s(1:2,3,1) = (/ 3*wc, 1.0d0 /)
!      s(1:2,3,1) = (/ 3, 1 /)
      Nparam(3,1) = 2
      dREF(:,1,3) = dREF(:,3,1)
      alpam(:,1,3) = alpbm(:,3,1)
      alpbm(:,1,3) = alpam(:,3,1)
      s(:,1,3) = s(:,3,1)
      Nparam(1,3) = Nparam(3,1)

!     C-U
      dREF(1:3,3,4) = (/ 5.2, 5.62, 7.48 /)
      alpam(1:3,3,4) = (/ 2.18, 1.75, 2.88 /)
      alpbm(1:3,3,4) = (/ 2.18, 2.64, -2.71 /)
      s(1:3,3,4) = (/ 2*wc, 1.0d0, 1.0d0 /)
      Nparam(3,4) = 3
      dREF(:,4,3) = dREF(:,3,4)
      alpam(:,4,3) = alpbm(:,3,4)
      alpbm(:,4,3) = alpam(:,3,4)
      s(:,4,3) = s(:,3,4)
      Nparam(4,3) = Nparam(3,4)

!     G-G
      dREF(1:4,1,1) =  (/ 6.24, 7.33, 6.24, 7.33/)
      alpam(1:4,1,1) = (/ 2.77, -2.93, 1.20, 1.27/)
      alpbm(1:4,1,1) = (/ 1.20, 1.27, 2.77, -2.93/)
      s(1:4,1,1) =   (/  2.0d0,  1.0d0,  2.0d0,  1.0d0 /)
      Nparam(1,1) = 4

!     G-U
      dREF(1:2,1,4) = (/ 5.20, 7.05 /)
      alpam(1:2,1,4) = (/ 2.53, -2.85 /)
      alpbm(1:2,1,4) = (/ 1.91, 1.43 /)
      s(1:2,1,4) = (/ 2*wc, 1.0d0 /)
      Nparam(1,4) = 2
      dREF(:,4,1) = dREF(:,1,4)
      alpam(:,4,1) = alpbm(:,1,4)
      alpbm(:,4,1) = alpam(:,1,4)
      s(:,4,1) = s(:,1,4)
      Nparam(4,1) = Nparam(1,4)

!     U-U
      dREF(1:4,4,4) =  (/ 5.39, 6.18, 5.39, 6.18 /)
      alpam(1:4,4,4) = (/ 1.71, 2.56, 2.61, -2.39 /)
      alpbm(1:4,4,4) = (/ 2.61, -2.39, 1.71, 2.56 /)
!      s(1:2,4,4) = (/ 2.0d0, 1 /)
      s(1:4,4,4) =   (/  1.0d0,  1.0d0,  1.0d0,  1.0d0 /)
      Nparam(4,4) = 3 !!! 3, deactivating the 4th (angluar params too close




C---------------------------------------------------------------

      if (natom .ne. natom_check) then
        write(*,*) 'ERROR : Number of atom does not match input file'
        stop
      endif

      do i=1, 3*natom_check
        xfort(i) = x(i)
      enddo

      do i=1, natom_check
        amasses(i) = amass(i)
        atomic_type(i) = text3(i)
      end do 

      open(57, file="cutoffs.dat")
      read(57, *)  rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $             rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $             rcut2_4b_out, rcut2_4b_in,
     $             rcut2_lj_out, rcut2_lj_in
      close(57)

c     We copy the box length and define the inverse box-length
      box_length = BoxL
      inv_box_length = 1.0d0 / box_length

      return
      end
      
