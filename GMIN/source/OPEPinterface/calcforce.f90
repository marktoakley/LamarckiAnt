      module calcforces
      contains
      subroutine calcforce_protein(scale,x,f,etot)
      use ion_pair
 
      parameter (MAXPRE = 1500)    !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord 
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)!! max number of SC-SC interactions 
      parameter (MAXBO  = MAXNAT)  !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles 
      parameter (MAXPHI = MAXNAT*4)  !! maximum number of torsional angles
      parameter (MAXTTY = 50000)       !! maximum number of residue name types 

      implicit double precision (a-h,o-z)

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3, &
                  NPHIA,NNB,NTYPES, MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3, &
                  NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC

      COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO), ICBH(MAXBO)
      COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH), JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
      COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI), ICP(MAXPHI)
      COMMON/ENER4/IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)
      
      COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH), &
                  PK(MAXPHI),PN(MAXPHI), &
                  PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60), &
                  GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)
      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno

      common/charge/cg,k1

      common/REWARD/ EHHB1
     
      COMMON/VP22/VPNE(MAXPRE)
      COMMON/INDHB/INDH(MAXPRE)
      COMMON/EHB/EHBT
      COMMON/TORS/QTOR

!      common/pos/x_natv(maxxc)
      common/nnumres/numres,Id_atom
      common/ncall/LTER
      real(8)  x,f 
      real(8)  scale     !! lambda for hamiltonian replica exchange it scales H-bond attraction
      real(8) IDIEL 
      DIMENSION X(MAXXC),F(MAXXC),CG(MAXNAT)
      DIMENSION FOHIG(MAXXC),FOLOW(MAXXC),fphi(maxxc)

      integer     numres(MAXNAT),Id_atom(MAXNAT)
      LOGICAL QTOR,QDAT,QBUG
      DATA ZERO/0.0d0/

      real(8) box_length, inv_box_length
      common/pbcBL/box_length, inv_box_length

      COMMON/DEBUG/QBUG
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
  

      dimension FIP(MAXXC)   !! Added by YC for OPEPv5


!---  COMPUTE THE FORCES F(*) and the ENERGY
!     FOR THE BOND LENGHTS, BOND ANGLES, DIHEDRALS, HYDROPHOBE-HYDROPHILIC
!     FOR THE NONBONDED AND PENALTIES (PROPENSITY, ETC...) SUCCESSFULLY


!--- Interface data Ion Pair   !! Added by YC for OPEPv5

      integer index_ip(mxip,3)

!---  reset the forces to zero


      f(1:natom3) = zero
      fip(1:natom3) = zero   !! Added by YC for OPEPv5

      EBONH = ZERO 
      EBONA = ZERO 
      ETHH = ZERO 
      ETHA = ZERO 
      EPH = ZERO 
      ENBPH = ZERO 
      EELPH = ZERO 
      ECNH = ZERO 
      EPA = ZERO 
      ENBPA = ZERO 
      EELPA = ZERO 
      ECNA = ZERO 
      EVDW = ZERO 
      EHBT = ZERO       !not in etot
      ELEC = ZERO 
      EHYDRO = ZERO 
      ETOT = ZERO
      EHHB1 = ZERO
      ener_phi = zero    
      ener_psi = zero !! NEW SHANGHAI    
      eextra = zero 
      e_ip = zero     ! salt bridge contribution
       ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA &
             +EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA &
             +ECNH+EHHB1+ECNA + ener_phi + eextra + ener_psi+e_ip
      QDAT = .FALSE.
!      LTER = LTER + 1
!      LTER = 0

      if (QDAT) THEN
      LTER = LTER + 1
!--   save data for optimization
      OPEN(UNIT=56,FILE="proteinA-opt.dat",status="unknown")
      REWIND(56)     
      ENDIF

! Transform the coordinates by translation and rotation
!      call Rotation(x_natv,x,delr,rmsd_natv,npart)


! ---- BOND LENTHS
      IF (NBONH .gt. 0) then
      CALL EBOND(NBONH,IBH,JBH,ICBH,X,F,EBONH,RK,REQ)
      ENDIF
      IF (NBONA .gt. 0) then
      CALL EBOND(NBONA,IB,JB,ICB,X,F,EBONA,RK,REQ)
      ENDIF
!     write(37,*) ' EBONH + EBONA ',EBONH+EBONA
!     stop !! SHANGHAI

! ---- BOND ANGLES 
      IF (NTHETH .gt. 0) then
      CALL ETHETA(MAXTH,NTHETH,ITH,JTH,KTH,ICTH,X,F,ETHH,TK,TEQ)
      ENDIF
      IF (NTHETA .gt. 0) then
      CALL ETHETA(MAXTH,NTHETA,IT,JT,KT,ICT,X,F,ETHA,TK,TEQ)
      ENDIF

!      open(unit=64,file="beta32.edih",status="unknown")
 
! ---- TORSIONS
      IF (NPHIH .gt. 0) then
      QTOR = .FALSE. 
      CALL ETORS(scale,MAXPHI,NPHIH,IPH,JPH,KPH,LPH,ICPH,CG,IAC,X, &
        F,EPH,ENBPH,EELPH,NPHIH,ECNH,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

      IF (NPHIA .gt. 0) then
      QTOR = .TRUE. 
      CALL ETORS(scale,MAXPHI,NPHIA,IP,JP,KP,LP,ICP,CG,IAC,X,F,EPA, &
        ENBPA,EELPA,NPHIA,ECNA,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

! --- SEPARATION SLOW AND HIGH
      fohig(1:natom3) = f(1:natom3)
      f(1:natom3) = zero
      
!      IF (NPHIA .gt. 0) then
!      QTOR = .TRUE. 
!      CALL ETORS(MAXPHI,NPHIA,IP,JP,KP,LP,ICP,CG,IAC,X,F,EPA,ENBPA,
!     +   EELPA,NPHIA,ECNA,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
!      ENDIF

!       close(64)
!       stop


! ---- Interface with module data


         call ip_interface(index_ip,n_ip)   !! Added by YC for OPEPv5


! ---- NONBONDED INTERACTIONS
       CALL ENBOND(scale,NATOM,CG,IAC,X,F,CN1,CN2,EVDW,ELEC,index_ip,&
          n_ip,ion_pair_control)     !! Modified by YC for OPEPv5

!----- HYDROPHOBIC-HYDROPHILIC INTERACTIONS
       CALL HYDROP(scale,X,F,EHYDRO,index_ip,n_ip,ion_pair_control)       !! Modified by YC for OPEPv5

! --- SEPARATION SLOW AND HIGH
      folow(1:natom3) = f(1:natom3)
      f(1:natom3) = zero

      call ephipv(x,f,ener_phi,NATOM) 
      call epsipv(x,f,ener_psi,NATOM)
!     call ephipvha(x,f,ener_phi,NATOM) 
!     call epsipvha(x,f,ener_psi,NATOM)


      fphi(1:natom3) = f(1:natom3)
      f(1:natom3)=zero


!----- SALT-BRIDGES INTERACTIONS    !! Added by YC for OPEPv5


      IF(ion_pair_control) THEN
         call ion_pair_force(x,f,e_ip)
         fip(1:natom3)=f(1:natom3)
      ENDIF

      
      f(1:natom3) = fphi(1:natom3) + folow(1:natom3) + fohig(1:natom3)&
          +fip(1:natom3)   !! Modified by YC for OPEPv5

!     Calculate the energy and force for the constrained energy
!     THE POTENTIAL IS EXPRESSED BY: RK_con * rmsd**2


!---  PRINT OUT THE ENERGY AND THE FORCE VECTOR
       ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
             EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
             ECNH+EHHB1+ECNA + ener_phi + eextra + ener_psi  + e_ip   !! Modified by YC for OPEPv5

!--- data to write out for test IP contribution ----  !! Added by YC for OPEPv5
       IF(ion_pair_control) THEN
          energy_ip(1)=ETOT
          energy_ip(2)=e_ip
       ENDIF

      if (QBUG) then
         write (*,*) '========================'
         write (*,*) 'Ehydro   ', EHYDRO
         write (*,*) 'Evdw     ', EVDW
         write (*,*) 'Eelec    ', ELEC
         write (*,*) 'Eph      ', EPH
         write (*,*) 'Epa      ', EPA
         write (*,*) 'Ethh     ', ETHH
         write (*,*) 'Etha     ', ETHA
         write (*,*) 'Ebonh    ', EBONH
         write (*,*) 'Ebona    ', EBONA
         write (*,*) 'Enbph    ', ENBPH
         write (*,*) 'Enbpa    ', ENBPA
         write (*,*) 'Eelph    ', EELPH
         write (*,*) 'Eelpa    ', EELPA
         write (*,*) 'Ecnh     ', ECNH
         write (*,*) 'Ehhb1    ', EHHB1
         write (*,*) 'Ecna     ', ECNA
         write (*,*) 'ener_phi ', ener_phi
         write (*,*) 'eextra   ', eextra
         write (*,*) 'ener_psi ', ener_psi
         write (*,*) 'e_ion    ', e_ip    !! added by YC for OPEPv5
         write (*,*) 'Etot     ', ETOT
         STOP
      endif

      return
      end subroutine calcforce_protein


!---------------------------------------------------------------------------
      subroutine calcforce_RNA(scale,x,f,etot,logenervar,simuPercent)
      
      use RNAnb
      use energies

      parameter (MAXPRE = 1500)    !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord 
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)!! max number of SC-SC interactions 
      parameter (MAXBO  = MAXNAT)  !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles 
      parameter (MAXPHI = MAXNAT*4)  !! maximum number of torsional angles
      parameter (MAXTTY = 50000)       !! maximum number of residue name types 

      implicit double precision (a-h,o-z)

      double precision simuPercent

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3, &
                  NPHIA,NNB,NTYPES, MBONA,MTHETA,MPHIA
      integer NATOM
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC

      COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO), ICBH(MAXBO)
      COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH), JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
      COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI), ICP(MAXPHI)
      COMMON/ENER4/IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)
      
      COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH), &
                  PK(MAXPHI),PN(MAXPHI), &
                  PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60), &
                  GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)
      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno

      common/charge/cg,k1

      COMMON/VP22/VPNE(MAXPRE)
      COMMON/INDHB/INDH(MAXPRE)
      COMMON/EHB/EHBT
      COMMON/TORS/QTOR

      common/scalingfactor/scaling_factor

!      common/pos/x_natv(maxxc)
      common/nnumres/numres,Id_atom
      common/ncall/LTER
      real(8)  x,f 
      real(8)  scale     !! lambda for hamiltonian replica exchange it scales H-bond attraction
      real(8) IDIEL 
      DIMENSION X(MAXXC),F(MAXXC),CG(MAXNAT)

      integer     numres(MAXNAT),Id_atom(MAXNAT)
      LOGICAL QTOR,QDAT,QBUG

      real(8) box_length, inv_box_length
      common/pbcBL/box_length, inv_box_length

      COMMON/DEBUG/QBUG
      integer i
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      logical logenervar

      double precision Frest(natom3)


!---  COMPUTE THE FORCES F(*) and the ENERGY
!     FOR THE BOND LENGHTS, BOND ANGLES, DIHEDRALS, HYDROPHOBE-HYDROPHILIC
!     FOR THE NONBONDED AND PENALTIES (PROPENSITY, ETC...) SUCCESSFULLY

!---  reset the forces to zero


      f(1:natom3) = 0.0d0

      EBONH = 0.0d0 
      EBONA = 0.0d0 
      ETHH = 0.0d0 
      ETHA = 0.0d0 
      EPH = 0.0d0 
      ENBPH = 0.0d0 
      EELPH = 0.0d0 
      ECNH = 0.0d0 
      EPA = 0.0d0 
      ENBPA = 0.0d0 
      EELPA = 0.0d0 
      ECNA = 0.0d0 
      EVDW = 0.0d0 
      EHBT = 0.0d0       !not in etot
      ELEC = 0.0d0 
      EHYDRO = 0.0d0 
      ETOT = 0.0d0
      EHHB = 0.0d0
      ener_phi = 0.0d0    
      ener_psi = 0.0d0 !! NEW SHANGHAI    
      eextra = 0.0d0 
       ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
             EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
             ECNH+EHHB+ECNA + ener_phi + eextra + ener_psi 
      QDAT = .FALSE. .or. qbug
!      LTER = LTER + 1
!      LTER = 0

      if (QDAT) THEN
      LTER = LTER + 1
!--   save data for optimization
      OPEN(UNIT=56,FILE="proteinA-opt.dat",status="unknown")
      REWIND(56)     
      ENDIF

! Transform the coordinates by translation and rotation
!      call Rotation(x_natv,x,delr,rmsd_natv,npart)


! ---- BOND LENTHS
      IF (NBONH .gt. 0) then
      call RNA_EBOND(NBONH,IBH,JBH,ICBH,X,F,EBONH,RK,REQ)
      ENDIF
      IF (NBONA .gt. 0) then
      call RNA_EBOND(NBONA,IB,JB,ICB,X,F,EBONA,RK,REQ)
      ENDIF

! ---- BOND ANGLES 
      IF (NTHETH .gt. 0) then
      call RNA_ETHETA(MAXTH,NTHETH,ITH,JTH,KTH,ICTH,X,F,ETHH,TK,TEQ)
      ENDIF
      IF (NTHETA .gt. 0) then
      call RNA_ETHETA(MAXTH,NTHETA,IT,JT,KT,ICT,X,F,ETHA,TK,TEQ)
      ENDIF

 
! ---- TORSIONS
      IF (NPHIH .gt. 0) then
      QTOR = .FALSE. 
      call RNA_ETORS(scale,NPHIH,IPH,JPH,KPH,LPH,ICPH,CG,IAC,X, &
        F,EPH,ENBPH,EELPH,ECNH,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

      IF (NPHIA .gt. 0) then
      QTOR = .TRUE. 
      call RNA_ETORS(scale,NPHIA,IP,JP,KP,LP,ICP,CG,IAC,X,F,EPA, &
        ENBPA,EELPA,ECNA,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

      

!----- HYDROPHOBIC-HYDROPHILIC INTERACTIONS
       call RNA_HYDROP(scale,X,F)

!---  PRINT OUT THE ENERGY AND THE FORCE VECTOR
       ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
             EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
             ECNH+ECNA + ener_phi + eextra + ener_psi

       Frest = 0
       if (simuPercent .ge. 1e-7) then
         call RNA_RESTRAINTS(X, Frest, Erest, nFconst, simuPercent)
         call RNA_POSRES(X, Frest, Eposres, simuPercent)
         Econst = Erest + Eposres
       endif

       ETOT = ETOT + Econst
       f(1:natom3) = f(1:natom3) + Frest(1:natom3)



       if (QDAT) then
         evec(1) = ebona
         evec(2) = etha
         evec(3) = epa
         evec(10) = ehhb
         do i=1, 17
           write(56,'(i4,4x,f8.3,4x,f15.10,A)') i, 1.0, evec(i)
         enddo
         write(56,'(A,f15.10)') 'Etot     ', ETOT
       endif
       if (QBUG) then
          write (*,*) '========================'
          write (*,*) 'bond     ', Ebona
          write (*,*) 'angle    ', Etha
          write (*,*) 'torsion  ', Epa
          write (*,*) 'Evdw     ', Evdw
          write (*,*) 'Ehhb     ', Ehhb
          write (*,*) 'Ebarrier ', Ehbr
          write (*,*) 'Ecoop    ', Ecoop
          write (*,*) 'Estak    ', Estak
          write (*,*) 'Ehydrop  ', Ehydro
          write (*,*) 'Econst   ', Econst
          write (*,*) 'Etot     ', ETOT
          do i=1, 4
            do j=1,4
              write (*,*) basenum(i), basenum(j), evec(bcoef(i,j))
            enddo
          enddo
          STOP
       endif

       if(logenervar) then 
         write(62,"(13(f10.3, 1x))") Ebona, etha, epa, &
          evec(5), evec(4), ehhb, ehbr, Ecoop, Estak, Econst, nFconst, &
          Etot, Etot*scaling_factor
       endif
 
       f(1:natom3) = f(1:natom3) * scaling_factor
       ETOT = ETOT * scaling_factor

      return
      end subroutine calcforce_RNA
!---------------------------------------------------------------------------
      subroutine ephipv(x,f,ener_phi,NATOM) 
      
!     Subroutine used to introduce Phi penalty energy outside regions.
!     term E=k*(phi-phi0)**2, and the corresponding force. 
!     if -160 < phi < -60, phi0 = phi, E=0; else, phi0=0, E=k*phi**2 
!     a smaller force constant k is used for Gly and Asp
!     no term for Pro (L) or Pro (D)
      
      implicit none
      
      real(8) pbc_mic
      
      integer MAXPRE,MAXNAT,MAXXC,NATOM
      parameter (MAXPRE = 1500)  !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT) !! maximum number of cart coord
      
      real(8) force_k, forceg_k, radian, rad2
      real(8) scaling_factor
!     parameter (force_k = (5.0d0)) !! *0.76d0) !scale 10 from OPEP3.0 - for MD
!     parameter (forceg_k = (1.5d0)) !! 3.0for Gly and ASP from OPEP3.0 - for MD
      parameter (force_k = (1.1d0)) !! best-vecteur 24 July06
      parameter (forceg_k = (0.5d0)) !! 	   
!     parameter (force_k = 10.0d0) !! *0.76d0) !scale 10 from OPEP3.0
!     parameter (forceg_k = 3.0d0) !! 3.0 for GlY and ASP from OPEP3.0
      parameter (radian = 57.29577951308232088d0) 
      parameter (rad2 = radian*radian)
      
      common/scalingfactor/scaling_factor
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom
      
      
      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
      common/propens2/walpha_foal(20),wbeta_fobe(20)
      
      integer nfrag, lenfrag,ichain 
      real(8)  walpha_foal, wbeta_fobe
      
      integer i,i3 
      integer Id_atom(MAXNAT), numres(MAXNAT) !Id of atom and residue
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(8)  x,f,deg           ! deg: first derivative of phi>0 energy term 
      dimension x(MAXXC), f(MAXXC), deg(MAXXC) 
      real(8) xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(8) xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(8) xca,yca,zca,xdb,ydb,zdb
      real(8) xt,yt,zt,dot_prq,rcb,rprq
      real(8) xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(8) phi_lower,phi_upper,phi,sin_phi,cos_phi,phi0  
      real(8) dt,dt2,e_phi,ener_phi,dedphi
      real(8) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(8) dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(8) dedxc,dedyc,dedzc,dedxd,dedyd,dedzd 
      real(8) t1,t2
      real(8) dr,drtr,drur,fphi,fdphi
      
      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
      
      integer natom3
      natom3 = natom*3
!     open(unit=38,file="beta32.ephi",status="unknown")
      
!     zero out the restraint energy term and first derivatives
      ener_phi = 0.0d0
      f(1:natom3) = 0.0d0
      deg(1:natom3) = 0.0d0
      
      phi_lower = -160.0d0
      phi_upper = -60.0d0 
      
      
!     Calculating the phi>0 penalty energy and the first derivative of this term 
      do i = 5, NATOM-4
        
        if (text4(i).eq. ' PRO   ') then
           go to 5200
        endif 
        
        if ( text3(i).eq.'  CA ' .and. (text4(i).ne. ' DPR   ')) then
           
           if (ichain(i-4) .eq. ichain(i+2)) then 
              
              i3 = i*3
              
              ixa = i3 - 14
              iya = i3 - 13
              iza = i3 - 12
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
               
!     The x, y, and z coordinates for N2 
              ixb = i3 - 8
              iyb = i3 - 7 
              izb = i3 - 6 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              ixc = i3 - 2
              iyc = i3 - 1 
              izc = i3 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for C2
              if (text4(i).ne. ' GLY   ') then
                 ixd = i3 + 4
                 iyd = i3 + 5 
                 izd = i3 + 6
              else
                 ixd = i3 + 1
                 iyd = i3 + 2 
                 izd = i3 + 3
              endif
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
              
!      If required we apply periodic boundary conditions
              if (periodicBC) then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )
                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif
               
!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_phi = dot_prq / rtru
               cos_phi = min(1.0d0,max(-1.0d0,cos_phi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_phi = rprq / (rcb*rtru)
               
               
!     calculating phi 
               phi = dacos(cos_phi) * radian 
               if (sin_phi .lt. 0.0d0) phi = -phi 
               if ((phi .gt. phi_lower) .and. (phi .lt. phi_upper)) then
                  phi0 = phi
               else if(phi.gt.phi_lower.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else if(phi.lt.phi_upper.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else
                  t1 = phi - phi_lower
                  t2 = phi - phi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     phi0 = phi_lower
                  else
                     phi0 = phi_upper
                  end if
               endif 
               dt = phi - phi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if 
!     Using E(phi) = k * (phi-phi0)**2
               dt2 = dt * dt
!     e_phi = force_k * dt2 !in degrees
               if (text4(i).ne.' GLY   '.and.text4(i).ne. ' ASP   ')then
                  fphi = force_k * scaling_factor / (rad2)
                  fdphi = fphi * radian
                  e_phi = fphi * dt2 !in radian
                  dedphi = 2.0d0 * fdphi * dt
               else
                  fphi = forceg_k * scaling_factor / (rad2)
                  fdphi = fphi * radian
                  e_phi = fphi * dt2 !in radian
                  dedphi = 2.0d0 * fdphi * dt
               endif
               
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha  
               xca = xc - xa
               yca = yc - ya
               zca = zc - za
!     x,y,z component of vector from C2 to N2 
               xdb = xd - xb
               ydb = yd - yb
               zdb = zd - zb
               
!     If required we apply periodic boundary conditions                 
               if(periodicBC)then
                  
                  xca = pbc_mic( xca ) 
                  yca = pbc_mic( yca ) 
                  zca = pbc_mic( zca ) 
                  
                  xdb = pbc_mic( xdb ) 
                  ydb = pbc_mic( ydb ) 
                  zdb = pbc_mic( zdb ) 
               endif
               
               dr = dedphi/rcb
               drtr = dr/rt2
               drur = dr/ru2
               dedxt =  drtr * (yt*zcb - ycb*zt)
               dedyt =  drtr * (zt*xcb - zcb*xt)
               dedzt =  drtr * (xt*ycb - xcb*yt)
               dedxu = -drur * (yu*zcb - ycb*zu)
               dedyu = -drur * (zu*xcb - zcb*xu)
               dedzu = -drur * (xu*ycb - xcb*yu)      
               
!     compute derivative components for this interaction
               dedxa = zcb*dedyt - ycb*dedzt
               dedya = xcb*dedzt - zcb*dedxt
               dedza = ycb*dedxt - xcb*dedyt
               dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxd = zcb*dedyu - ycb*dedzu
               dedyd = xcb*dedzu - zcb*dedxu
               dedzd = ycb*dedxu - xcb*dedyu
               
!     increment the overall energy term and derivatives
               
               ener_phi = ener_phi + e_phi !total energy of phi>0 
               
               deg(ixa) = deg(ixa) + dedxa
               deg(iya) = deg(iya) + dedya
               deg(iza) = deg(iza) + dedza
               deg(ixb) = deg(ixb) + dedxb
               deg(iyb) = deg(iyb) + dedyb
               deg(izb) = deg(izb) + dedzb
               deg(ixc) = deg(ixc) + dedxc
               deg(iyc) = deg(iyc) + dedyc
               deg(izc) = deg(izc) + dedzc
               deg(ixd) = deg(ixd) + dedxd
               deg(iyd) = deg(iyd) + dedyd
               deg(izd) = deg(izd) + dedzd
               
             endif
           endif
!     enddo  !! loop j
        endif 
 5200   continue
      enddo
      
      f(1:natom3) = -deg(1:natom3)
      
      
      return
      end subroutine ephipv 
      
      
      subroutine epsipv(x,f,ener_psi,NATOM) 

!     This subroutine is used to restrict the Psi region 
!     term E=k*(psi-psi0)**2, and the corresponding force. 
!     if -60 < psi < 160, phi0 = phi, E=0; else, phi0=0, E=k*phi**2 
!     a smaller force constant k is used for Gly and Asp
!     no term for Pro (L) or Pro (D)
      
      implicit none
      
      real(8) pbc_mic
       
      integer MAXPRE,MAXNAT,MAXXC,NATOM
      parameter (MAXPRE = 1500)  !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT) !! maximum number of cart coord
      
      real(8) force_k, forceg_k, radian, rad2
      real(8) scaling_factor
!     parameter (force_k = (5.0d0)) ! *0.76d0) ! scale 10 from OPEP3.0   for MD
!     parameter (forceg_k = (1.5d0)) ! 3.0 for Gly from OPEP3.0           for MD
      parameter (force_k = (1.1d0)) !! best-vecteur 24 July06
      parameter (forceg_k = (0.5d0)) !! 
!     parameter (force_k = 10.0d0) ! *0.76d0) ! scale 10 from OPEP3.0 
!     parameter (forceg_k = 3.0d0) ! 3.0 for Gly and Asp from OPEP3.0 
      parameter (radian = 57.29577951308232088d0) 
      parameter (rad2 = radian*radian)
      
      common/scalingfactor/scaling_factor
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom
      
      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
      integer nfrag, lenfrag,ichain 
      
      
      integer i,i3 
      integer Id_atom(MAXNAT), numres(MAXNAT) !Id of atom and residue
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(8)  x,f,deg           ! deg: first derivative of psi>0 energy term 
      dimension x(MAXXC), f(MAXXC), deg(MAXXC) 
      real(8) xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(8) xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(8) xca,yca,zca,xdb,ydb,zdb
      real(8) xt,yt,zt,dot_prq,rcb,rprq
      real(8) xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(8) psi_lower,psi_upper,psi,sin_psi,cos_psi,psi0  
      real(8) dt,dt2,e_psi, ener_psi,dedpsi
      real(8) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(8) dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(8) dedxc,dedyc,dedzc,dedxd,dedyd,dedzd 
      real(8) t1,t2
      real(8) dr,drtr,drur,fpsi,fdpsi
      
      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
      
      integer natom3
      natom3 = natom*3
!     open(unit=38,file="beta32.epsi",status="unknown")
      
! zero out the restraint energy term and first derivatives
      ener_psi = 0.0d0
      
      deg(1:natom3) = 0.0d0
      
      psi_lower = -60.0D0
      psi_upper = 160.0d0 
      
      
!     Calculating the psi>0 penalty energy and the first derivative of this term 
      do i = 3, NATOM-4
        if (text4(i).eq. ' PRO   ') then
           go to 5200
        endif 
        
        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' DPR   ' ) then !!           
           
           if (ichain(i-2) .eq. ichain(i+4)) then 
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following: 
!     The x, y, and z coordinates for C1
              ixa = i3-8
              iya = i3-7
              iza = i3-6
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for N2 
              ixb = i3-2
              iyb = i3-1 
              izb = i3 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              if (text4(i).ne. ' GLY   ') then
                 ixc = i3+4
                 iyc = i3+5 
                 izc = i3+6
              else
                 ixc = i3+1
                 iyc = i3+2 
                 izc = i3+3 
              endif
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
              
!     The x, y, and z coordinates for C2
              if (text4(i).ne. ' GLY   ') then
                 ixd = i3+10
                 iyd = i3+11 
                 izd = i3+12
              else
                 ixd = i3+7
                 iyd = i3+8 
                 izd = i3+9 
              endif
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
              
!     If required, we apply periodic boundary conditions
              if(periodicBC)then
                 xba = pbc_mic( xba ) 
                 yba = pbc_mic( yba ) 
                 zba = pbc_mic( zba ) 
                 xcb = pbc_mic( xcb ) 
                 ycb = pbc_mic( ycb ) 
                 zcb = pbc_mic( zcb ) 
                 xdc = pbc_mic( xdc ) 
                 ydc = pbc_mic( ydc ) 
                 zdc = pbc_mic( zdc ) 
              endif

!        cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_psi = dot_prq / rtru
               cos_psi = min(1.0d0,max(-1.0d0,cos_psi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_psi = rprq / (rcb*rtru)
!     calculating psi 
               psi = dacos(cos_psi) * radian 
               if (sin_psi .lt. 0.0d0) psi = -psi 
               if ((psi .gt. psi_lower) .and. (psi .lt. psi_upper)) then
                  psi0 = psi
               else if(psi.gt.psi_lower.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else if(psi.lt.psi_upper.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else
                  t1 = psi - psi_lower
                  t2 = psi - psi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     psi0 = psi_lower
                  else
                     psi0 = psi_upper
                  end if
               endif
               dt = psi - psi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
                 
!     Using E(psi) = k * (psi-psi0)**2
               dt2 = dt * dt
!     e_psi = force_k * dt2 !in degrees
               if (text4(i).ne.' GLY   '.and.text4(i).ne. ' ASP   ')then
                    fpsi = force_k * scaling_factor / (rad2)
                    fdpsi = fpsi * radian
                    e_psi = fpsi * dt2 !in radian
                    dedpsi = 2.0d0 * fdpsi * dt
                 else
                    fpsi = forceg_k * scaling_factor / (rad2)
                    fdpsi = fpsi * radian
                    e_psi = fpsi * dt2 !in radian
                    dedpsi = 2.0d0 * fdpsi * dt
                 endif
                
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha  
                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2 
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb
                 
!     if required, we apply periodic boundary conditions
                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )
                    
                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif
                 
                 dr = dedpsi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)           
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu
                 
!     increment the overall energy term and derivatives
                 
                 ener_psi = ener_psi + e_psi !total energy of psi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd
                 
              endif
           endif
!     enddo  !! loop j 
        endif 
 5200   continue
      enddo
      f(1:natom3) = f(1:natom3) - deg(1:natom3)
      return
      end subroutine epsipv
      
!---------------------------------------------------------------------------
      subroutine ephipvha(x,f,ener_phi,NATOM) 
      
!     This subroutine is used to calculate Phi positive penalty energy 
!     term E=k*(phi-phi0)**2, and the cooresponding force. 
!     if -180 < phi < 0, phi0 = phi, E=0; else, phi0=0, E=k*phi**2 
      
      implicit none
      
      real(8) pbc_mic
      
      integer MAXPRE,MAXNAT,MAXXC,NATOM
      parameter (MAXPRE = 1500)  !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT) !! maximum number of cart coord
      
      real(8) force_k, radian, rad2
      parameter (force_k = 7.6d0) ! 10.0d0*0.76d0) !scale 10 from OPEP3.0
      parameter (radian = 57.29577951308232088d0) 
      parameter (rad2 = radian*radian)
      
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom
      
      
      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)

      integer nfrag, lenfrag,ichain 

      integer i,i3 
      integer Id_atom(MAXNAT), numres(MAXNAT) !Id of atom and residue
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(8)  x,f,deg           ! deg: first derivative of phi>0 energy term 
      dimension x(MAXXC), f(MAXXC), deg(MAXXC) 
      real(8) xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(8) xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(8) xca,yca,zca,xdb,ydb,zdb
      real(8) xt,yt,zt,dot_prq,rcb,rprq
      real(8) xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(8) phi_lower,phi_upper,phi,sin_phi,cos_phi,phi0  
      real(8) dt,dt2,e_phi, ener_phi,dedphi
      real(8) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(8) dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(8) dedxc,dedyc,dedzc,dedxd,dedyd,dedzd 
      real(8) t1,t2
      real(8) dr,drtr,drur,fphi,fdphi

      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)

      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      integer  natom3
      natom3 = natom*3
      

!     zero out the restraint energy term and first derivatives
      ener_phi = 0.0d0
      f(1:natom3) = 0.0d0
      deg(1:natom3) = 0.0d0
      
      phi_lower = -160.0d0
      phi_upper = -60.0d0 
      

      do i = 7, NATOM-4
        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' GLY ') then
           if (ichain(i-4) .eq. ichain(i+2)) then 
              
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following: 
!     The x, y, and z coordinates for C1
              ixa = i3 - 14
              iya = i3 - 13
              iza = i3 - 12
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for N2 
              ixb = i3 - 8
              iyb = i3 - 7 
              izb = i3 - 6 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              ixc = i3 - 2
              iyc = i3 - 1 
              izc = i3 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for C2
              ixd = i3 + 7
              iyd = i3 + 8 
              izd = i3 + 9
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
!-----------------------------------------------------------RL    ----------------------
              
              if(periodicBC)then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )

                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif

!-----------------------------------------------------------RL    ----------------------
!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_phi = dot_prq / rtru
               cos_phi = min(1.0d0,max(-1.0d0,cos_phi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_phi = rprq / (rcb*rtru)
!     calculating phi 
               phi = dacos(cos_phi) * radian 
               if (sin_phi .lt. 0.0d0) phi = -phi 
               if ((phi .gt. phi_lower) .and. (phi .lt. phi_upper)) then
                  phi0 = phi
               else if(phi.gt.phi_lower.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else if(phi.lt.phi_upper.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else
                  t1 = phi - phi_lower
                  t2 = phi - phi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     phi0 = phi_lower
                  else
                     phi0 = phi_upper
                  end if
               endif
               dt = phi - phi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               
!     Using E(phi) = k * (phi-phi0)**2
               dt2 = dt * dt
!     e_phi = force_k * dt2 !in degrees
               fphi = force_k / (rad2)
               fdphi = fphi * radian
               e_phi = fphi * dt2 !in radian
               dedphi = 2.0d0 * fdphi * dt
               
                
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha  


                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2 
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb
!-----------------------------------------------------------RL    ----------------------

                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )

                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif

!-----------------------------------------------------------RL    ----------------------

                 dr = dedphi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)    
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu

!     increment the overall energy term and derivatives

                 ener_phi = ener_phi + e_phi !total energy of phi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd

              endif
           endif
!     enddo  !! loop j
        endif 
      enddo
      
      f(1:natom3) = -deg(1:natom3)

      return
      end subroutine ephipvha

!     --------------------------------------------------------------------
      subroutine epsipvha(x,f,ener_psi,NATOM) 

!     This subroutine is used to calculate Phi positive penalty energy 
!     term E=k*(phi-phi0)**2, and the cooresponding force. 
!     if -180 < phi < 0, phi0 = phi, E=0; else, phi0=0, E=k*phi**2 

      implicit none
      
      real(8) pbc_mic

      integer MAXPRE,MAXNAT,MAXXC,NATOM
      parameter (MAXPRE = 1500)  !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT) !! maximum number of cart coord
      
      real(8) force_k, radian, rad2
      parameter (force_k = 7.6d0) !10*0.76d0) ! scale 10 from OPEP3.0 
      parameter (radian = 57.29577951308232088d0) 
      parameter (rad2 = radian*radian)

      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom

      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)

      integer nfrag, lenfrag,ichain 


      integer i,i3 
      integer Id_atom(MAXNAT), numres(MAXNAT) !Id of atom and residue
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(8)  x,f,deg           ! deg: first derivative of psi>0 energy term 
      dimension x(MAXXC), f(MAXXC), deg(MAXXC) 
      real(8) xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(8) xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(8) xca,yca,zca,xdb,ydb,zdb
      real(8) xt,yt,zt,dot_prq,rcb,rprq
      real(8) xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(8) psi_lower,psi_upper,psi,sin_psi,cos_psi,psi0  
      real(8) dt,dt2,e_psi, ener_psi,dedpsi
      real(8) dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(8) dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(8) dedxc,dedyc,dedzc,dedxd,dedyd,dedzd 
      real(8) t1,t2
      real(8) dr,drtr,drur,fpsi,fdpsi

      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT)

      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      integer natom3
      natom3 = natom*3
      

!     zero out the restraint energy term and first derivatives
      ener_psi = 0.0d0
      deg(1:natom3) = 0.0d0
      
      psi_lower = -60.0D0
      psi_upper = 160.0d0 
      

      do i = 1, NATOM-5

        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' GLY ') then
           if (ichain(i-2) .eq. ichain(i+4)) then 
              
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following: 
!     The x, y, and z coordinates for N1
              ixa = i3 - 8
              iya = i3 - 7
              iza = i3 - 6
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for C_alpha 
              ixb = i3 - 2
              iyb = i3 - 1 
              izb = i3 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C(O)
              ixc = i3 + 7
              iyc = i3 + 8 
              izc = i3 + 9 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for N2 
              ixd = i3 + 13
              iyd = i3 + 14 
              izd = i3 + 15
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc

!-----------------------------------------------------------RL    ----------------------
              
              if(periodicBC)then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )

                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif

!-----------------------------------------------------------RL    ----------------------

!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_psi = dot_prq / rtru
               cos_psi = min(1.0d0,max(-1.0d0,cos_psi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_psi = rprq / (rcb*rtru)
!     calculating psi 
               psi = dacos(cos_psi) * radian 
               if (sin_psi .lt. 0.0d0) psi = -psi
               if ((psi .gt. psi_lower) .and. (psi .lt. psi_upper)) then
                  psi0 = psi
               else if(psi.gt.psi_lower.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else if(psi.lt.psi_upper.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else
                    t1 = psi - psi_lower
                    t2 = psi - psi_upper
                    if (t1 .gt. 180.0d0) then
                       t1 = t1 - 360.0d0
                    else if (t1 .lt. -180.0d0) then
                       t1 = t1 + 360.0d0
                    end if
                    if (t2 .gt. 180.0d0) then
                       t2 = t2 - 360.0d0
                    else if (t2 .lt. -180.0d0) then
                       t2 = t2 + 360.0d0
                    end if
                    if (abs(t1) .lt. abs(t2)) then
                       psi0 = psi_lower
                    else
                       psi0 = psi_upper
                    end if
                 endif
                 dt = psi - psi0
                 if (dt .gt. 180.0d0) then
                    dt = dt - 360.0d0
                 else if (dt .lt. -180.0d0) then
                    dt = dt + 360.0d0
                 end if
                 
!     Using E(psi) = k * (psi-psi0)**2
                 dt2 = dt * dt
!     e_psi = force_k * dt2 !in degrees
                 fpsi = force_k / (rad2)
                 fdpsi = fpsi * radian
                 e_psi = fpsi * dt2 !in radian
                 dedpsi = 2.0d0 * fdpsi * dt
                 
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha  
                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2 
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb

!-----------------------------------------------------------RL    ----------------------

                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )

                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif

!-----------------------------------------------------------RL    ----------------------


                 dr = dedpsi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)    
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu

!     increment the overall energy term and derivatives

                 ener_psi = ener_psi + e_psi !total energy of psi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd

              endif
           endif
!     enddo  !! loop j 
        endif 
      enddo
      f(1:natom3) = f(1:natom3) - deg(1:natom3)

      return
      end subroutine epsipvha
!----------------------------------------------------



!=========================================================================
!=======================      writing the pdb file       =================
!=======================       with center of mass       =================       
!=======================         inside the box          =================
!=========================================================================


      subroutine writeCM(NATOM,X)

      implicit none

      integer MAXPRE,MAXNAT,MAXXC,MAXTTY
      parameter (MAXPRE = 1500)    !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord
      parameter (MAXTTY = 50000)       !! maximum number of residue name types 

      integer NATOM, ii, jj, chain_length
      integer a,b, ll
      
      real(8) AMASS_2(MAXNAT)
      real(8) x(MAXXC)
      real(8) xx1(MAXNAT)
      real(8) yy1(MAXNAT)
      real(8) zz1(MAXNAT)

      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno
  
      real(8) box_length, inv_box_length
      common/pbcBL/box_length, inv_box_length


      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)

      integer nfrag, lenfrag,ichain


      real(8)  cm1, cm2, cm3, total_mass
      dimension cm1(MAXPRE),cm2(MAXPRE),cm3(MAXPRE)


      ll = 0

      xx1 = x(1:3*NATOM:3)
      yy1 = x(2:3*NATOM:3)
      zz1 = x(3:3*NATOM:3)
  
      AMASS_2(1:NATOM) = 1.0d0/AMASS(1:NATOM)

      do ii = 1, nfrag
        chain_length = lenfrag(ii)
        total_mass = 0.0d0
        cm1(ii) = 0.0d0
        cm2(ii) = 0.0d0
        cm3(ii) = 0.0d0

        do jj = 1, chain_length
          total_mass = total_mass + AMASS_2(jj+ll)
          cm1(ii)=cm1(ii)+(AMASS_2(jj+ll)*xx1(jj+ll))
          cm2(ii)=cm2(ii)+(AMASS_2(jj+ll)*yy1(jj+ll))
          cm3(ii)=cm3(ii)+(AMASS_2(jj+ll)*zz1(jj+ll))
        end do 

        a = ll + 1
        b = ll + lenfrag(ii)
        cm1(ii) = cm1(ii) / total_mass
        cm2(ii) = cm2(ii) / total_mass
        cm3(ii) = cm3(ii) / total_mass
        xx1(a:b)=xx1(a:b)-box_length*dnint(cm1(ii)*inv_box_length)
        yy1(a:b)=yy1(a:b)-box_length*dnint(cm2(ii)*inv_box_length)
        zz1(a:b)=zz1(a:b)-box_length*dnint(cm3(ii)*inv_box_length)

        ll = ll + chain_length
      end do

      x(1:3*NATOM:3) = xx1
      x(2:3*NATOM:3) = yy1
      x(3:3*NATOM:3) = zz1
    
      return
 
      end subroutine writeCM

      endmodule calcforces
