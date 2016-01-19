C-----------------------------------------------------------------------
c  copy/rewrite from rdtop(nf,chrg) in other-2005-last.f
c  as a helper subroutine
c
c  Xiao Dong, 10/2006
c
c
      SUBROUTINE readtop1(nbondh,nbonda,nbonds)


      implicit double precision (a-h,o-z)


      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA


      COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP

      nbondh = nbonh
      nbonda = nbona
      nbonds = numbnd
      
      RETURN
      END

      subroutine readtop2(nbondt,bia,bib,beq,beq2)
!    $ redu1_mass, redu2_mass)
      
      parameter (MAXPRE = 1500)    !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord 
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)!! max number of SC-SC interactions 
      parameter (MAXBO  = MAXNAT)  !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles 
      parameter (MAXPHI = MAXNAT*4)  !! maximum number of torsional angles
      parameter (MAXTTY = 50000)  !! 
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs 

      implicit double precision (a-h,o-z)


      COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),
     $             ICBH(MAXBO)

      COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),
     $             PK(MAXPHI),PN(MAXPHI),
     $             PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),
     $             GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA


      integer :: it
      integer :: nbondt
      integer, dimension(nbondt) :: bia, bib
      double precision, dimension(nbondt) :: beq,beq2
      
      it = 0
      do i = 1, nbonh
        it = it + 1
        i3 = ibh(i);  i3 = i3+1;  i3 = i3+2
        j3 = jbh(i);  j3 = j3+1;  j3 = j3+2
        if (mod(i3,3)/=0) stop "i3"
        if (mod(j3,3)/=0) stop "j3"
        i3=i3/3; bia(it) = i3
        j3=j3/3; bib(it) = j3
        
        ic = icbh(i)
        beq(it) = req(ic)
        beq2(it) = beq(it)*beq(it)
      end do


      do i = 1, nbona
        it = it + 1
        i3 = ib(i);  i3 = i3+1;  i3 = i3+2
        j3 = jb(i);  j3 = j3+1;  j3 = j3+2
        if (mod(i3,3)/=0) stop "i3"
        if (mod(j3,3)/=0) stop "j3"   
        i3=i3/3; bia(it) = i3
        j3=j3/3; bib(it) = j3

        ic = icb(i)
        beq(it) = req(ic)
        beq2(it) = beq(it)*beq(it)
      end do
      

      RETURN
      END

