      SUBROUTINE INITIALISE_PROTEIN(P_B_C,BoxL,CofM,NATOM_CHECK,
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

      use ion_pair                 !! Added by YC for OPEPv5        

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
      dimension xfort(3*NATOM_CHECK), amasses(NATOM_CHECK)
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

cdx   COMMON/NPAIRCB/NPAIR2,IPAIR(MAXPAI), JPAIR(MAXPAI),
cdx  1 NPAIR(MAXNAT,MAXNAT)

      COMMON/HYDRO/ rncoe(maxpnb),vamax(maxpnb),
     1 ni(maxpnb),nj(maxpnb),nstep,nb,
     2 ivi(maxpnb),ivj(maxpnb),
     $ epshb_mcmc(maxpnb)
      double precision rncoe, vamax, epshb_mcmc
      integer*8 ni, nj, ivi, ivj, nstep, nb

      COMMON/REWARD/ EHHB1
      COMMON/PROPE/ RESP(MAXPRE)
      COMMON/VP22/VPNE(MAXPRE)
      COMMON/INDHB/indh(MAXPRE)
      COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
      common/pos/x_natv(maxxc)
c      common/forcek/force_const,rmsd
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom

      common/propens/ialpha(MAXNAT),ibeta(MAXNAT),icoeff(MAXPAI),
     1 foal(20),walpha(20),fobe(20),wbeta(20)
      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
      common/propens2/walpha_foal(20), wbeta_fobe(20)

      common/scor/score(272),score_RNA(17)
      common/JANV06/NHB,ref_NHB,NHBH,NHBB

      logical qbug
      common/debug/qbug
      logical use_qbug_
C---------------------------------------------------------------------	  
      logical P_B_C,CofM
      real*8 BoxL
C---------------------------------------------------------------------	  

      REAL*8 IDIEL
      DIMENSION X(MAXXC),CG(MAXNAT)
      DIMENSION FOAL1(20),FOBE1(20)

      character*10 noms
      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)
      integer     numres(MAXNAT),Id_atom(MAXNAT)
      integer  ref_NHB

      CHARACTER*3 RESP 
      DATA ZERO/0.0d0/

      real*8 score_temp

      common/cacascsc/ct0lj(maxpnb),ct2lj(maxpnb)

      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $               rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $               rcut2_4b_out, rcut2_4b_in,
     $               rcut2_lj_out, rcut2_lj_in

      real*8 box_length, inv_box_length

      common/pbcBL/box_length, inv_box_length

      character(50) dummy
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
!       write(*,*) 'CG', CG(1)


C---  read INPUT PDB file
      CALL init_remd(path_for_inits,single_conf_init)
      if(.not. single_conf_init)  then 
         nf=41
         ns=43 
         OPEN(UNIT=nf,FILE=path_for_inits,status="unknown")
         CALL RDPDB(nf,X,text2,Id_atom,text3,text4,numres) 
         close(nf)
      else  
         nf=14
         ns=34
         OPEN(UNIT=nf,FILE="conf_initiale.pdb",status="unknown")
         CALL RDPDB(nf,X,text2,Id_atom,text3,text4,numres) 
         close(nf)
      endif

c --- Ion-Pair-Control
c --- Selection Ion Pair
c --- Read the potential and derivative on a grid       

      IF(ion_pair_control) THEN
         call ip_select(natom,text3,numres)
         call ip_potential_init
      ENDIF



C --- READ THE PARAMETERS FOR THE SIDE-CHAIN SIDE-CHAIN INTERACTIONS
      open(unit=28,file='parametres.list',status='unknown')
      do ikl = 1,nres 
      read(28,7205) resp(ikl),ialpha(ikl),ibeta(ikl)
      enddo 
      nb = 0
 7201 continue
      nb = nb + 1
      if (nb .eq. 1) then
      do i=1,20
      walpha(i) = 1.0d0
      wbeta(i) = 1.0d0
      enddo
      endif

      read(28,7200) ni(nb),ivi(nb),nj(nb),ivj(nb),rncoe(nb),vamax(nb),
     1 icoeff(nb)
 7200 format(i4,4x,2i4,4x,i4,f8.3,f13.3,i9)
 7205 format(1x,a,2i7)

      if(ni(nb) .ne. -999) goto 7201
      nb = nb - 1
      read(28,'(A)') noms
      read(28,*) NVALR  !! nbre residus with phi et psi
      close(28)

c      NVALR = NVALR + 1
      IF (NB .GE. MAXPNB) THEN
      PRINT*,' ERRORS in the dimension of the variables associated'
      PRINT*,' with the side-chain side-chain interactions '
      STOP
      ENDIF 

      call prop(foal1,fobe1)
      do i=1,20
      foal(i) = foal1(i)
      fobe(i) = fobe1(i)
      enddo
c---  NEW JANV05
C---- Set-up the vector associating each atom with a fragment of the protein
C----    NFRAG : number of fragments
C----    length_fragment: table containing the length of each fragment
C----    ichain : the fragment to which atom idatom belongs to

C---  read ichain file
      nf=55
      OPEN(UNIT=nf,FILE="ichain.dat",status="unknown")
      read(55,*)  NFRAG 
      do i=1, NFRAG
         read(55,*) id, lenfrag(i)
      enddo
      idatom = 1
      do i=1, NFRAG
         do j=1, lenfrag(i)
            ichain(idatom) = i
            idatom = idatom +1
         end do
      end do
      close(55)
c      do ikl=1,NATOM
c      enddo

C---  read weight file for proper OPEP
      OPEN(UNIT=55,FILE="scale.dat",status="unknown")
      do ih=1,272
         read(55,4900) ikl, score_temp, dummy
         score(ih)=dble(score_temp)
      enddo
 4900 format(i4,f9.10,A)
      close(55)
     
c -- scale for MD, 30/01/06  - rescaled for MD 21/02/2006
c    with the parameters of bestVect.final-24July06
      do ih=1,266
        score(ih) = dble(scaling_factor)*score(ih) 
      enddo
      score(270) = dble(scaling_factor)*score(270)

c --- NEW SEPT 06
      score(267) = score(267)*3.0d0/1.6d0
      score(268) = score(268)*4.0d0 ! 3.0 before
c -- end scale MD

      
c -- end scale MD  

c -- the following lines must be used for ART and MD
      do ih=225,244 
         walpha(ih-224) = score(ih)
      enddo
      do ih=245,264 
         wbeta(ih-244) = score(ih) 
      enddo
c --  end modif for ART and MD

      
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
      if(periodicBC)then
        box_length = BoxL
!       write(*,*) 'box_length', box_length
       inv_box_length = 1.0d0 / box_length
      endif
      
      WEIHB14 = score(222)      !! ONE  !! weight for H-bonds helix
      WEIHB15 = score(223)      !! ONE  !! weight for others intra and inter
      
      WEICA = score(224)        !! 1.0d0 !! weight for CA-CA VdW
      
      do i = 1, nb
        ct0lj(i) = rncoe(i)
        ct2lj(i) = vamax(i)*vamax(i)
        if (icoeff(i).eq.-1) then
           ct0lj(i) = ct0lj(i) * WEICA !! ponderation CA-CA
        else if (icoeff(i) /= -2) then
           if (rncoe(i) .gt. 0.) then
              ct0lj(i) = 1.5D0*score(ICOEFF(I)) !! WEISCSC(ICOEFF(I))= 1.5D0*score(ICOEFF(I))
           else
              ct0lj(i) = score(ICOEFF(I)) !! WEISCSC(ICOEFF(I))= score(ICOEFF(I))
           endif
           ct0lj(i) = rncoe(i) *  ct0lj(i) !! WEISCSC(ICOEFF(I))   !! ponderation Sc-Sc
        end if
        
        
        
c     ---     first Set-up the right parameters for each H-bond
        if (ichain(ni(i)) == ichain(nj(i))) then !! INTRA VS. INTER-CHAINS
           IF (ivi(i) /= (ivj(i)-4) .and. ivi(i) /= (ivj(i)+4) ) then
              EPSHB = 2.0d0 * WEIHB15 !! this corresponds to j > i+4
           else
              EPSHB = 1.25d0 * WEIHB14 !! this corresponds to j = i+4
           endif
           IF (ivi(i) == (ivj(i)-5)) then
              EPSHB = 2.25d0 *  WEIHB15 !! 1.75 23-dec to PI-helix
           endif
c     ---       bug resolved 8 dec 05 Ph.D.
           idif = abs(ivi(i)-ivj(i))
           IF (idif >= 15 .and. idif <= 30)then
              EPSHB = 0.75d0 *  WEIHB15 !! pour 1FSD
           endif
        else
           EPSHB = 2.0d0 * WEIHB15 !! INTER-CHAIN INTERACTIONS
        endif                   !! INTRA VS INTER-CHAINS
        epshb_mcmc(i) = epshb
      end do
      
      do i = 1, 20
        walpha_foal(i) = walpha(i)*foal(i)
        wbeta_fobe(i) = wbeta(i)*fobe(i)
      end do
      return
      end
      
C -------------------------------------------------------

      subroutine prop(foal,fobe)



c ------- energies of residues for alpha and beta
c ------- positive values are required for bad propensities

      implicit double precision (a-h,o-z)

      DIMENSION fobe(20),foal(20)

      i=1  !! CYS
      foal(i) = 0.6 
      fobe(i) = 0.2  

      i=2  !! LEU 
      foal(i) = 0.19
      fobe(i) = 0.2
c      fobe(i) = 0.0
     
      i=3  !! VAL
      foal(i) = 0.46
      fobe(i) = -0.3
     
      i=4  !! ILE 
      foal(i) = 0.35
      fobe(i) = -0.3
     
      i=5  !! MET 
      foal(i) = 0.21
      fobe(i) = 0.3
     
      i=6  !! PHE 
      foal(i) = 0.47
      fobe(i) = -0.3
     
      i=7  !! TYR 
      foal(i) = 0.47
      fobe(i) = -0.3
     
      i=8  !! LYS 
      foal(i) = 0.15
      fobe(i) = 0.3

      i=9  !! ARG 
      foal(i) = 0.06
      fobe(i) = 0.3
     
      i=10  !! PRO 
      foal(i) = 1.3 
      fobe(i) = 0.3
     
      i=11  !! GLY 
      foal(i) = 1.10 
      fobe(i) = 0.3
c      fobe(i) = 0.0
     
      i=12  !! ALA 
      foal(i) = 0.1 
      fobe(i) = 0.3  
c      fobe(i) = 0.0 !! 
     
      i=13  !! GLN 
      foal(i) = 0.32
      fobe(i) = 0.3

      i=14  !! HIS 
      foal(i) = 0.62
      fobe(i) = 0.3

      i=15  !! ASN 
      foal(i) = 0.60
      fobe(i) = 0.3
c      fobe(i) = 0.0  !!

      i=16  !! ASP 
      foal(i) = 0.59
      fobe(i) = 0.3

      i=17  !! GLU 
      foal(i) = 0.34
      fobe(i) = 0.3

      i=18  !! SER 
      foal(i) = 0.52
      fobe(i) = 0.3

      i=19  !! THR 
      foal(i) = 0.57
      fobe(i) = -0.3

      i=20  !! TRP 
      foal(i) = 0.47
      fobe(i) = 0.3

      return 
      end

