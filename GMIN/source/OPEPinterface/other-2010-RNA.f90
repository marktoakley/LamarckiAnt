module energies
      double precision evdw,elec,eph,epa,ethh,etha,ebonh,ebona, &
        enbph,enbpa,eelph,eelpa,ehydro, &
        Ehhb, Ehbr, Ecoop, Estak, Eplane, & 
        Econst, Erest, Eposres, nFconst
!     energy vector, for the optimization
      double precision Evec(35)
end module energies


!> @brief
!
subroutine RNA_RESTRAINTS(X, F, Econ, nFcon, simtime) 

  use geometric_corrections
  use restraint_params
  implicit none

  double precision, intent(in) :: X(*), simtime
  double precision, intent(out) :: F(*), Econ, nFcon
  double precision v, escale

  integer idx
  double precision diff(3), dlen, df(3)

  escale = exp(-80*cos(simtime*3.14159/5)**2)

  Econ = 0
  nFcon = 0
  do idx = 1, Nrests
    diff = X(resti(idx)*3-2:resti(idx)*3)-X(restj(idx)*3-2:restj(idx)*3)
    dlen = euc_norm(diff)
!     Econ = Econ + restk(idx)*(dlen - restl(idx))**2
!     df = 2*restk(idx)*(dlen - restl(idx))*diff/dlen

!    linear for d > 2
    v = dlen - restl(idx)
    if((trest(idx) .ge. 2) .and. (abs(v) .gt. 2)) then
      Econ = Econ + restk(idx)*(4*abs(v)-4)
      df = restk(idx)*4*v/abs(v)*diff/dlen
    else
      Econ = Econ + restk(idx)*v**2
      df = 2*restk(idx)*v*diff/dlen
    endif

    if(trest(idx) .eq. 2) then
      Econ = Econ*escale
      df = df*escale
    endif

    nFcon = nFcon + euc_norm(df)

    F(resti(idx)*3-2 : resti(idx)*3) = F(resti(idx)*3-2 : resti(idx)*3) - df
    F(restj(idx)*3-2 : restj(idx)*3) = F(restj(idx)*3-2 : restj(idx)*3) + df

    restl(idx) = restl(idx) + drestl(idx)
  enddo

  end subroutine RNA_RESTRAINTS


!> @brief
!
  subroutine RNA_POSRES(X, F, Econ, simtime)

  use geometric_corrections
  use restraint_params

  implicit none

  double precision, intent(in) :: X(*), simtime
  double precision, intent(out) :: F(*), Econ
  double precision v, escale

  integer idx
  double precision diff(3), df(3)

  Econ = 0
  do idx = 1, Nposres
    diff = X(pri(idx)*3-2:pri(idx)*3) - prx(:,idx)

    v = dot_product(diff, diff)

    Econ = Econ + prk(idx)*v
    df = 2*prk(idx)*diff

    F(pri(idx)*3-2 : pri(idx)*3) = F(pri(idx)*3-2 : pri(idx)*3) - df
  enddo
  if (Nposres .gt. 0) then
    prx = prx + prdx
  endif


  end subroutine RNA_POSRES


!> @brief
!> @todo Update the indentation
      subroutine RNA_EBOND(NBON,IB,JB,ICB,X,F,EBON,RK,REQ)
!
!     THE POTENTIAL IS EXPRESSED BY: RK * (RIJ-REQ)**2
!
      implicit none

      integer NBON,IB(*),JB(*),ICB(*)
      double precision X(*),F(*),EBON,RK(*),REQ(*)

      real(8) pbc_mic

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      LOGICAL qbug
      COMMON/debug/qbug

      integer*8 jn, I3, J3, ic
      double precision xa(3), rij, da, df, enerb

      logical QDET

      QDET = .FALSE. .or. qbug
      if(QDET) then
         open(unit=7,file="beta32.ebond",status="unknown",position="append")
         write(7,*) '  i ','  j ','   rij  ','   req  ','  enerb ','  force '
      endif

      EBON = 0.0d0

      DO JN = 1,NBON
      I3 = IB(JN)
      J3 = JB(JN)
      xa = X((I3+1):(I3+3))-X((J3+1):(J3+3))

        if(periodicBC)then
           xa = pbc_mic(xa)   !; print *, "xij", xij
        endif

      RIJ = dsqrt(dot_product(xa,xa))
      IC = ICB(JN)
      DA = RIJ-REQ(IC)
      DF = RK(IC)*DA*score_RNA(1)
      ENERB = DF*DA
      if(QDET) then
         if (enerb .ge. 0.5) then
            write(7,'(i4,i4,f8.3,f8.3,f8.3,f8.3)') I3/3+1,J3/3+1,RIJ,REQ(IC),ENERB,DF
         endif
      endif
      DF = (DF+DF)/RIJ
      xa = DF*xa
      F((I3+1):(I3+3)) = F((I3+1):(I3+3)) - xa
      F((J3+1):(J3+3)) = F((J3+1):(J3+3)) + xa
      EBON = EBON + ENERB
      ENDDO

      if(QDET) then
         close(7)
      endif
      RETURN

end subroutine RNA_EBOND


!*********************************************************************
!>
!>
!*********************************************************************
      subroutine RNA_ETHETA(MAXTT,NTHETH,IT,JT,KT,ICT,X,F,ETHH,TK,TEQ)
!
!     THE POTENTIAL IS EXPRESSED BY: TK * (ANG(I,J,K)-TEQ)**2
!
      implicit none

      real(8) pbc_mic

      integer, intent(in) :: IT(*),JT(*),KT(*),ICT(*), NTHETH, MAXTT
      double precision, intent(in) :: X(*), TK(*),TEQ(*)
      double precision, intent(out) :: ETHH, F(*)


      LOGICAL qbug
      COMMON/debug/qbug

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
     
      integer :: i3, j3, k3, ic, jn, p1, p2, p3

      double precision :: ct0, ct1, ct2
      
      double precision rIJ(3), rKJ(3)
      double precision rDI(3), rDJ(3), rDK(3)
      double precision ANT, RKJ0, RIJ0, RIK0, EAW, DFW, DA, DF

      double precision pt999
      DATA pt999 /0.9990d0/
      logical QDET

      QDET = .FALSE. .or. qbug

      if(QDET) then
         open(unit=7,file="beta32.etheta",status="unknown")
         write(7,*) "   P1  ","   P2  ","   P3  ", "     t    ", &
                   "    teq   ","   diff   ","  energy  "
      endif
      ETHH = 0.0d0

      DO JN = 1,NTHETH
        IC = ICT(JN)
        if (TK(IC) .ge. 2.0d0) then
           I3 = IT(JN)
           J3 = JT(JN)
           K3 = KT(JN)
           rIJ = X(I3+1:I3+3)-X(J3+1:J3+3)
           rKJ = X(K3+1:K3+3)-X(J3+1:J3+3)

           if(periodicBC)then
              rIJ = pbc_mic(rIJ)
              rKJ = pbc_mic(rKJ)
           endif

           RIJ0 = dot_product(rIJ, rIJ)
           RKJ0 = dot_product(rKJ, rKJ)
           RIK0 = dsqrt(RIJ0*RKJ0)
           CT0 = dot_product(rIJ, rKJ)/RIK0
           CT1 = MAX(-pt999,CT0)
           CT2 = MIN(pt999,CT1)
           ANT = DACOS(CT2)


!     ENERGY
           DA = ANT-TEQ(IC)
           DF = TK(IC)*DA*score_RNA(2)
           EAW = DF*DA
           DFW = -(2*DF)/DSIN(ANT)
           if(QDET) then
               if (EAW .ge. 5.0d0) then
                  P1 = IT(JN)/3 +1
                  P2 = JT(JN)/3 +1
                  P3 = KT(JN)/3 +1
                  write(7,'(i7,i7,i7, f10.3,f10.3,f10.3,f10.3)') P1,P2,P3, ANT*180/3.14, &
                      TEQ(IC)*180/3.14,DA*180/3.14, EAW
               endif
            endif
            ETHH = ETHH + EAW
!     FORCE

            rDI = DFW*(rKJ/RIK0-CT2*rIJ/RIJ0)
            rDK = DFW*(rIJ/RIK0-CT2*rKJ/RKJ0)
            rDJ = -rDI-rDK
            F(I3+1:I3+3) = F(I3+1:I3+3)-rDI
            F(J3+1:J3+3) = F(J3+1:J3+3)-rDJ
            F(K3+1:K3+3) = F(K3+1:K3+3)-rDK
         endif
      ENDDO

      if(QDET) then
         close(7)
      endif
!      write(71,*)ETHH

      RETURN
      end subroutine RNA_ETHETA



!-----------------------------------------------------------------------
      subroutine RNA_ETORS(lambda,NPHI,IP,JP,KP,LP,ICP,CG,IAC, &
          X,F,EP,ENBP,EELP,ECN,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)

      use geometric_corrections
      implicit none

      real(8) pbc_mic

      logical QDET

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      double precision CUT,SCNB,SCEE,IDIEL,DIELC
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC

      LOGICAL qbug
      COMMON/debug/qbug

      integer :: i3, j3, k3, l3, k3t, l3t, ic, ii, jj, jn
      integer :: p1, p2, p3, p4
      integer :: ibig, isml, idumi, iduml, inc, kdiv

      double precision rIJ(3), rKJ(3), rKL(3), rD(3), rG(3), rIL(3)
      double precision lenD, lenG, dotDG
      double precision vfmul, FMULN, vCPHI, vSPHI
      double precision vEPW, vDF, ENW, EEW, F1, F2
      double precision COSNP, SINNP, CT0, CT1, DF0, DF1, DFLIM, DFN
      double precision SCNB0, SCEE0
      double precision AP1, dums, r2, r6, rfac, z10, z20, z12
      double precision rFI(3), rFJ(3), rFK(3), rFL(3)
      double precision rDC(3), rDC2(3), rDR1(3), rDR2(3), rDR(3)
      double precision rA(3)


      integer, intent(in) :: IP(*),JP(*),KP(*),LP(*),ICP(*),IPN(*),IAC(*), NPHI
      double precision, intent(in) :: CG(*), X(*), GAMC(*), GAMS(*), FMN(*), CN1(*), CN2(*), PK(*), PN(*)
      double precision, intent(out) :: F(*), EP, ENBP, EELP, ECN

      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      double precision eqangle, curangle

      double precision :: GMUL(10), tm24, tm06, tenm3, pi
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00, &
          0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      DATA PI/3.141592653589793d+00/

      real(8) lambda

      QDET = .FALSE. .or. qbug
      if(QDET) then
         open(unit=7,file="beta32.etors",status="unknown")
         write(7,*) "   P1  ", "   P2  ", "   P3  ", "   P4  ", &
          "    PK   ","    PN   ","   AP1   ","   GAMC  ","   Etors "
      endif
      EP = 0
      ECN = 0
      EELP = 0
      SCNB0 = 1.0d0/SCNB
      SCEE0 = 1.0d0/SCEE

      DO JN = 1,NPHI
        I3 = IP(JN)
        J3 = JP(JN)
        K3T = KP(JN)
        L3T = LP(JN)
        K3 = IABS(K3T)
        L3 = IABS(L3T)
        rIJ = X(I3+1:I3+3) - X(J3+1:J3+3)
        rKJ = X(K3+1:K3+3) - X(J3+1:J3+3)
        rKL = X(K3+1:K3+3) - X(L3+1:L3+3)

        if(periodicBC)then
           rIJ = pbc_mic( rIJ )
           rKJ = pbc_mic( rKJ )
           rKL = pbc_mic( rKL )
        endif

        rD = crossproduct(rIJ, rKJ)
        rG = crossproduct(rKL, rKJ)

        lenD = dsqrt(dot_product(rD,rD)+TM24)
        lenG = dsqrt(dot_product(rG,rG)+TM24)
        dotDG = dot_product(rD,rG)

        z10 = 1.0d0/lenD
        z20 = 1.0d0/lenG
        if (tenm3 .gt. lenD) z10 = 0
        if (tenm3 .gt. lenG) z20 = 0
        Z12 = Z10*Z20
        vFMUL = 0
        if (z12 .ne. 0.0d0) vFMUL = 1.0d0

        CT0 = MIN(1.0d0,dotDG*Z12)
        CT1 = MAX(-1.0d0,CT0)

        AP1 = PI-DSIGN(DACOS(CT1),dot_product(rKJ,crossproduct(rG,rD)))
!        vCPHI = DCOS(AP1)
        vCPHI = -CT1
        vSPHI = DSIN(AP1)

!     ----- ENERGY AND THE DERIVATIVES WITH RESPECT TO
!           COSPHI -----

        IC = ICP(JN)
        INC = IPN(IC)
        CT0 = PN(IC)*AP1
        COSNP = DCOS(CT0)
        SINNP = DSIN(CT0)
        vEPW= (PK(IC)+COSNP*GAMC(IC)+SINNP*GAMS(IC))*vFMUL !! might be revised
        DF0 = PN(IC)*(GAMC(IC)*SINNP-GAMS(IC)*COSNP)
        DUMS = vSPHI+SIGN(TM24,vSPHI)
        DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
        df1 = df0/dums
        if(tm06.gt.abs(dums)) df1 = dflim
        vDF = DF1*vFMUL

        vEPW = vEPW*score_RNA(3)
        vDF = vDF*score_RNA(3)

!       if(PK(IC)> 0.0d0) print *,'..',DFLIM,GAMC(IC)*PN(IC)

        if(QDET) then
           if(vEPW .ge. 5) then
             P1 = I3/3+1
             P2 = J3/3+1
             P3 = K3/3+1
             P4 = L3/3+1
             eqangle = atan2(gams(IC)/pk(IC), gamc(IC)/pk(IC))*180/PI
             curangle = atan2(sinnp, cosnp)*180/PI
             write(7,'(f7.0, f7.0, f7.0, f7.0, f9.3,f9.3,f9.3,f9.3,f9.3)') &
             P1, P2, P3, P4, PK(IC),PN(IC), curangle, eqangle, vEPW
           endif
        endif
!     END ENERGY WITH RESPECT TO COSPHI


!     ----- DC = FIRST DER. OF COSPHI W/RESPECT
!           TO THE CARTESIAN DIFFERENCES T -----
        rDC = -rG*Z12-vCPHI*rD*Z10**2
        rDC2 = rD*Z12+vCPHI*rG*Z20**2
!     ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
        rDR1 = vDF*(crossproduct(rKJ,rDC))
        rDR2 = vDF*(crossproduct(rKJ,rDC2))
        rDR = vDF*(crossproduct(rIJ,rDC) + crossproduct(rDC2, rKL))
        rFI = - rDR1
        rFJ = - rDR + rDR1
        rFK = + rDR + rDR2
        rFL = - rDR2
!     ----- CALCULATE 1-4 NONBONDED CONTRIBUTIONS

        rIL = X(I3+1:I3+3)-X(L3+1:L3+3)

        if(periodicBC)then
           rIL = pbc_mic(rIL)
        endif

        IDUMI = SIGN(1,K3T)
        IDUML = SIGN(1,L3T)
        KDIV = (2+IDUMI+IDUML)/4
        FMULN = dble(kdiv)*FMN(ICP(JN))
        II = (I3+3)/3
        JJ = (L3+3)/3
        IBIG = MAX0(IAC(II),IAC(JJ))
        ISML = MIN0(IAC(II),IAC(JJ))
        IC = IBIG*(IBIG-1)/2+ISML
        R2 = FMULN/dot_product(rIL,rIL)
        R6 = R2**3
        rfac = R6*score_RNA(3)
        F1 = CN1(IC)*R6*rfac
        F2 = CN2(IC)*rfac
        ENW = F1-F2
        if (IDIEL.gt.0) then
          EEW = CG(II)*CG(JJ)*dsqrt(R2)*SCEE0
          DFN =((-12.0d0*F1+6.0d0*F2)*SCNB0-EEW)*R2
        else
          EEW = CG(II)*CG(JJ)*R2*SCEE0
          DFN =((-12.0d0*F1+6.0d0*F2)*SCNB0-(2*EEW))*R2
        endif
        rA = rIL*DFN
        rFI = rFI - rA
        rFL = rFL + rA

        enbp = enbp + enw  !! 1-4 nb
        eelp = eelp + eew  !! 1-4 elec
!      ----- THE TOTAL FORCE VECTOR -----
!
        F(I3+1:I3+3) = F(I3+1:I3+3) + (rFI)
        F(J3+1:J3+3) = F(J3+1:J3+3) + (rFJ)
        F(K3+1:K3+3) = F(K3+1:K3+3) + (rFK)
        F(L3+1:L3+3) = F(L3+1:L3+3) + (rFL)

        ep   = ep + vEPW  !! torsions
      enddo
      ENBP = ENBP*SCNB0
      if(QDET) then
         close(7)
      endif
!      write(72,*)ep
      RETURN
      end subroutine RNA_ETORS



!-----------------------------------------------------------------------
!
!>     @brief Routine to calculate the hydrophobic/hydrophilic forces 
!>     and the H-bond forces.
!
!>     @details
!>     The analytic form includes now the propensity of residues to
!>     prefer alpha or beta states and the weights for all contributions
!>     has been rescaled (to be published in 2006).
!
!-----------------------------------------------------------------------
      subroutine RNA_HYDROP(lambda,X,F)

      use RNAnb
      use rnabase
      use energies

      implicit none

      real(8)  lambda     !! lambda, for hamiltonian replica exchange it scales H-bond attraction
      double precision X(*),F(*)

      real(8) pbc_mic

      integer MAXPRE, MAXNAT, MAXTTY
      parameter (MAXPRE = 1500)    !! maximum number of residues
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXTTY = 50000)       !! maximum number of residue name types 

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3, &
          NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA
      integer NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3, &
                  NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA


      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
      integer nfrag, lenfrag,ichain 

      double precision ehhb1, estak_t, evdw_t
!     double precision ehbrp, e4b

      logical QDET
      logical qbug
      common/debug/qbug

      COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
      double precision amass
      integer iac, nno

      double precision rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                    rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                    rcut2_4b_out, rcut2_4b_in, &
                    rcut2_lj_out, rcut2_lj_in
      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                    rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                    rcut2_4b_out, rcut2_4b_in, &
                    rcut2_lj_out, rcut2_lj_in
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
      
      integer ti, tj, tk, tl, i, j, k, l
      integer prev_res
!#      integer bpairs(15,MAXPRE), nbpairs(MAXPRE)

      double precision df, da2, a(3),dx(3) 
 
      logical hbexist

!#      nbpairs = 0
!#      bpairs = 0
      evec = 0

      QDET = .false. .or. qbug

      if (QDET) then
      open(unit=37,file="beta32.hydrop",status="unknown")
      open(unit=77,file="beta32.lj",status="unknown")
      write(77,*) ' ni ',' nj ', '    ct0lj', '   vamax', &
       '   distance','        evdw  ','        df    ', 'icoeff'
      open(unit=78,file="beta32.hb",status="unknown")
      write(78,*) ' ni ',' nj ','    ti  ','    tj  ','  dho2  ', &
       '  shb2  ','ebarrier','  ehhb  ','  Vangl ', '  ca0a  ','  ca0b  ',' alpa   ',' alpb   '
      open(unit=79,file="ehyd",status="unknown")
      write(79,*) '    Ecumul   ','   ', '     Ehyd    '
      open(unit=94,file="beta32.hbarrier",status="unknown")
      write(94,"(3a4, 6a8, a15)") 'I','J','K','dhoa','R1','dhob','R2', &
          'gaussa', 'gaussb','EHBR'
      open(unit=90,file="beta32.coop",status="unknown")
      write(90,"(4a4, 11a12)")  'I','J','K','L','ECOOP','DH1','DH2', &
        'R1','R2','DD1','DD2','VHB1','VHB2','CP','CM'
      endif


!--   set-up some parameters
      evdw = 0
      ESTAK = 0

      prev_res = 1
      do i = 1, NRES

!        ! Intra-base nb interactions
!        do j = blist(i), blist(i)-l, -1
!          tj = iac(j)
!          do k = j-3, prev_res, -1
!            tk = iac(k)
!            a = x(j*3-2:j*3) - x(k*3-2:k*3)
!            da2 = dot_product(a,a)
!            call RNA_lj(da2,evdw_t,df,nbct2(tj,tk),
!     $           rcut2_caca_scsc_in,rcut2_caca_scsc_out)
!            evdw = evdw + evdw_t*nbcoef(tj,tk)
!            dx = df*nbcoef(tj,tk)*a
!            F((j*3-2):j*3) = F((j*3-2):j*3) - dx
!            F((k*3-2):k*3) = F((k*3-2):k*3) + dx
!            evec(nbscore(tj,tk)) = evec(nbscore(tj,tk)) + evdw_t
!            
!            if(qdet .and. evdw_t >= 5) then
!                  write(77,1200) tj, tk, nbcoef(tj,tk), nbct2(tj,tk),
!     $              dsqrt(da2), evdw_t, df, nbscore(tj,tk)
!BER1200             format(i4,i4,f10.3,f10.3,f10.3,f15.3,f15.3,i4)
!            endif
!          enddo
!        enddo

        ! Inter-base nb interactions
        j = i
        do 
          j = j + 1
          if (j .gt. NRES) then
            exit
          endif


          k = prev_res
          l = blist(j-1)+1
          a = x(k*3-2:k*3) - x(l*3-2:l*3)
          if (periodicBC) then
            a = pbc_mic(a)
          endif
          da2 = dot_product(a,a)
          tk = iac(k)
          tl = iac(l)
          call RNA_debye(da2, evdw_t, df, chrg(tk), chrg(tl))

          evdw_t = evdw_t * nbcoef(tk,tl)
          df = df * nbcoef(tk,tl)
          evdw = evdw + evdw_t
          dx = df*a
          F((k*3-2):k*3) = F((k*3-2):k*3) - dx
          F((l*3-2):l*3) = F((l*3-2):l*3) + dx
          evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
          ! check the P-P distance, if it's big enough,
          ! just skips the whole residue altogether

!          if ( da2 > 6400) then
!            j = j + 10
!            cycle
!          else if ( da2 > 4900) then
!            j = j + 8
!            cycle
!          else if ( da2 > 3600) then
!            j = j + 6
!            cycle
!          else if ( da2 > 2500) then
!            j = j + 4
!            cycle
!          else if ( da2 > 1600) then
!            j = j + 2
!            cycle
          if ( da2 > 400) then
            cycle
          endif

!--------------------- STACKING -------------------------------------------
          call RNA_Stackv(blist(i),blist(j),F,X,estak_t)
          ESTAK = ESTAK + estak_t

          do k = prev_res, blist(i)
            tk = iac(k)
            
            l = blist(j-1)
            do while (l .lt. blist(j))
              l = l+1
!            do l = blist(j-1)+1, blist(j)
              tl = iac(l)
              a = x(k*3-2:k*3) - x(l*3-2:l*3)
              if (periodicBC) then
                a = pbc_mic(a)
              endif
              da2 = dot_product(a,a)
 
              if (da2 > rcut2_caca_scsc_out) then
                cycle
              endif

              call RNA_lj( da2, evdw_t, df, nbct2(tk,tl), rcut2_caca_scsc_in, rcut2_caca_scsc_out)

              evdw_t = evdw_t * nbcoef(tk,tl)
              df = df * nbcoef(tk,tl)
              ! if it's the next base, scale down the interaction
!              if( j .eq. i+1) then
!                evdw_t = evdw_t / 2
!                df = df / 2
!              endif

              evdw = evdw + evdw_t
              dx = df*a
              F((k*3-2):k*3) = F((k*3-2):k*3) - dx
              F((l*3-2):l*3) = F((l*3-2):l*3) + dx
              evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
              if(qdet .and. evdw_t >= 5) then
                !write(77,1200) k, l, nbcoef(tk,tl), nbct2(tk,tl), dsqrt(da2), evdw_t, df, nbscore(tk,tl)
                write(77,'(i4,i4,f10.3,f10.3,f10.3,f15.3,f15.3,i4)') &
                      k, l, nbcoef(tk,tl), nbct2(tk,tl), dsqrt(da2), evdw_t, df, nbscore(tk,tl)
              endif
            enddo
          enddo
        enddo
        prev_res = blist(i)+1
      enddo


!      evdw = 0
!
!      prev_res = 1
!      do i = 1, NRES
!        if (btype(i) <= 2) then
!          k = 1
!        else
!          k = 0
!        endif
!        do j = prev_res, blist(i)
!          tj = iac(j)
!          do l = j+3, NATOM
!            ! we consider i-i backbone-base interactions
!            if (l < blist(i)-k) then
!              cycle
!            endif
!            tl = iac(l)
!            a = x(j*3-2:j*3) - x(l*3-2:l*3)
!            da2 = dot_product(a,a)
!
!            call RNA_lj(da2,evdw_t,df,nbct2(tj,tl),
!     $           rcut2_caca_scsc_in,rcut2_caca_scsc_out)
!            evdw = evdw + evdw_t*nbcoef(tj,tl)
!            dx = df*nbcoef(tj,tl)*a
!            F((j*3-2):j*3) = F((j*3-2):j*3) - dx
!            F((l*3-2):l*3) = F((l*3-2):l*3) + dx
!            evec(nbscore(tj,tl)) = evec(nbscore(tj,tl)) + evdw_t
!          enddo
!        enddo
!        prev_res = blist(i)+1
!      enddo


!-----------------------------------------
! -- Interaction between nucleobases
!-----------------------------------------
      EHHB = 0
      do i = 1, NRES-1
        ti = btype(i)
        if (ti .gt. 4) then
          cycle
        endif
        do j = i+1, NRES
          tj = btype(j)

          if (tj .gt. 4) then
            cycle
          endif
          call RNA_BB(blist(i), ti, blist(j)-1, tj, score_RNA(10), X, F, Ehhb1, hbexist )
          if (.not. hbexist) then
            cycle
          endif
 
          if (qdet) then
             evec(int(bcoef(ti,tj))) = evec(int(bcoef(ti,tj))) + ehhb1

!             if (abs(EHHB1) .ge. 2) then
!               write(37,*) '  '
!               write(37,*) 'energy is counted'
!               write(37,*) i,j+1,dsqrt(DHO2),EHHB1,
!     $              DACOS(CA0a)*180/3.14,DACOS(CA0b)*180/3.14,
!     $              rncoe(i)
!               write(37,*) '  '
!             endif
          endif

!--------------------- TOTAL ENERGY HB -------------------------------------------
          EHHB = EHHB + EHHB1

!#          nbpairs(i) = nbpairs(i) + 1
!#          nbpairs(j) = nbpairs(j) + 1
!#          bpairs(nbpairs(i), i) = j
!#          bpairs(nbpairs(j), j) = i
       enddo
      enddo

!c--------------------- STACKING -------------------------------------------
!      ESTAK = 0
!      do i = 1, NRES - 1
!        do j = i+1, NRES
!          call RNA_Stackv(blist(i),blist(j),F,X,estak_t)
!          ESTAK = ESTAK + estak_t
!        enddo
!      enddo


!#!---------------------- BIFURCATION BARRIER --------------------------------------
!#      EHBR = 0
!#      do i = 1, NRES
!#        ! If we have more than two base pairs for this residue...
!#
!#        if (nbpairs(i) .lt. 2) then
!#          cycle
!#        endif
!#        li = blist(i)
!#        ti = btype(i)
!#        ! call the barrier function for /I\
!#        !                              J   K
!#        do j = 1, nbpairs(i)-1
!#          lj = blist(bpairs(j,i))
!#          tj = btype(bpairs(j,i))
!#          do k = j+1, nbpairs(i)
!#            lk = blist(bpairs(k,i))
!#            tk = btype(bpairs(k,i))
!#c            if (abs(lj-lk) .ge. 24) then
!#c            print *,'*',lj,li,lk
!#              call RNA_HBarrier(lj,li,lk,X,F,EHBRp,tj,ti,tk)
!#              EHBR = EHBR + EHBRp
!#c            endif
!#          enddo
!#        enddo
!#      enddo
!#
!#!---------------------------------   4 BODY HB ---------------------------------------
!#      Ecoop = 0
!#      do i = 1, NRES-1
!#        li = blist(i)
!#        ti = btype(i)
!#        do j= 1, nbpairs(i)
!#          lj = blist(bpairs(j,i))
!#          tj = btype(bpairs(j,i))
!#           if (lj .le. li) then
!#            cycle
!#          endif
!#!          do k = i+1, NRES-3
!#          do k = i+1, NRES
!#            lk = blist(k)
!#            tk = btype(k)
!#            ! this is to avoid trying 1 - 3 3 - 5
!#            if (lj .eq. lk) then
!#              cycle
!#            endif
!#            do l = 1, nbpairs(k)
!#              ll = blist(bpairs(l,k))
!#              tl = btype(bpairs(l,k))
!#              ! this is to avoid trying 1 - 5 3 - 5
!#              if (ll.le.lk .or. lj.eq.ll) then
!#                cycle
!#              endif
!#              ! i -> j or i -> j  ?
!#              ! k -> l    l <- k
!#              call RNA_Coop(li,lj,lk,ll,X,F,E4B,ti,tj,tk,tl)
!#              Ecoop = Ecoop + E4B
!#              call RNA_Coop(li,lj,ll,lk,X,F,E4B,ti,tj,tl,tk)
!#              Ecoop = Ecoop + E4B
!#            enddo
!#          enddo
!#        enddo
!#      enddo

      ! Mg - all VdW interactions
      do i = blist(NRES)+1, NATOM
        do j = 1, NATOM
          if(j .gt. blist(NRES) .and. j .le. i) then
            cycle
          endif
          a = x(i*3-2:i*3) - x(j*3-2:j*3)
          if (periodicBC) then
            a = pbc_mic(a)
          endif
          da2 = dot_product(a,a)
          tk = iac(i)
          tl = iac(j)
          call RNA_lj(da2, evdw_t, df, nbct2(tk,tl), rcut2_caca_scsc_in, rcut2_caca_scsc_out)

          evdw_t = evdw_t * nbcoef(tk,tl)
          df = df * nbcoef(tk,tl)
          evdw = evdw + evdw_t
          dx = df*a
          F((i*3-2):i*3) = F((i*3-2):i*3) - dx
          F((j*3-2):j*3) = F((j*3-2):j*3) + dx
          evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
        enddo

        ! Mg++ - all debye interaction
        do j = 1, NATOM
          if(j .gt. blist(NRES) .and. j .le. i) then
            cycle
          endif
          tl = iac(j)
          ! skip atoms besides P
          if(tl .ne. 3 .and. tl .ne. 13) then
            cycle
          endif
          tk = iac(i)
          a = x(i*3-2:i*3) - x(j*3-2:j*3)
          if (periodicBC) then
            a = pbc_mic(a)
          endif
          da2 = dot_product(a,a)
          call RNA_debye(da2, evdw_t, df, chrg(tk), chrg(tl))

          evdw_t = evdw_t * nbcoef(tk,tl)
!          print *, da2, evdw_t, nbcoef(tk,tl)
          df = df * nbcoef(tk,tl)
 !         print *, df
          evdw = evdw + evdw_t
          dx = df*a
          F((i*3-2):i*3) = F((i*3-2):i*3) - dx
          F((j*3-2):j*3) = F((j*3-2):j*3) + dx
          evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
        enddo
      enddo



      if(qbug) then
         evec(17) = Ehbr
         evec(7) = Estak
         evec(11) = Ecoop
      endif

      Ehydro = Ehhb + Ehbr + Estak + Ecoop

      RETURN
      end subroutine RNA_HYDROP


!-----------------------------------------------------------------------
      subroutine RNA_switch_cutoff(r2,ene_switched,for_switched,ri2,ro2)

        implicit none

        double precision r2,ene_switched,for_switched,ri2,ro2
        double precision rd6
        double precision sw_func,d_sw_func

        rd6 = 1.0d0/(ro2-ri2)**3

        sw_func = (ro2-r2)**2*(ro2+2.0d0*r2-3.0d0*ri2)*rd6
        d_sw_func = 12.0d0*(ro2-r2)*(ri2-r2)*rd6 !*r1

!dx     for_switched = for_switched*r1
        for_switched = for_switched*sw_func - ene_switched*d_sw_func !/r1
!dx     for_switched = for_switched/r1

        ene_switched = ene_switched*sw_func

      return
      end subroutine RNA_switch_cutoff


!-----------------------------------------------------------------------
      subroutine RNA_lj(da2,eahyd,df,ct2,r2in,r2out) !!NEW - NO LONGER LJ
        implicit none

        double precision, intent(in) ::  da2,ct2,r2in,r2out
        double precision, intent(out) :: eahyd, df
        double precision r

        r = dsqrt(da2)
!        Eahyd = exp(2.0*(ct2-r)) * (ct2/r)**6
!        DF = -2*Eahyd*(3+r)/da2
!        Eahyd = exp(4.0*(ct2-r))
!        DF = -4.0*Eahyd/r
!        Eahyd = exp(4.0*(ct2-r)) - exp(-r)
!        DF = -(4.0*exp(4.0*(ct2-r)) - exp(-r))/r 
!        Eahyd = exp(4.0*(ct2+0.2-r)) - exp(3.0*(ct2+0.1-r))
!        DF = -(4.0*exp(4.0*(ct2+0.2-r)) - 3.0*exp(3.0*(ct2+0.1-r)))/r 
        Eahyd = exp(4.0*(ct2-r)) - exp(3.0*(ct2-0.2-r))
        DF = -(4.0*exp(4.0*(ct2-r)) - 3.0*exp(3.0*(ct2-0.2-r)))/r 

! ------- store energy in ehydro and forces
        if (da2>=r2in) then
          call RNA_switch_cutoff(DA2,eahyd,DF,r2in,r2out)
        endif

      return
      end subroutine RNA_lj

      subroutine RNA_debye(da2,eahyd,df,qi,qj)
        implicit none

        double precision score, score_RNA
        common/scor/score(272),score_RNA(17)
        double precision, intent(in) ::  da2, qi, qj
        double precision, intent(out) :: eahyd, df
        double precision r
        double precision pi, Dlength
        !data pi/3.141592653589793d+00/

        Dlength = score_RNA(6)

        r = dsqrt(da2)
        ! Debye-Huckel prefactor : 1/(4*pi*e_0*e_r)
        !   1/(4*pi*e_0) ~= 332.056 kcal/mol * A
        !   1/(4*pi*e_0*e_r) ~= 332.056/80 = 4.1507
        ! Debye length = 1/sqrt(8*pi*l_b*I)
        !   l_b ~= 7
        !   I : Ionic strength
        !   DL ~= 3.04 / sqrt(I)
        Eahyd = (4.1507*qi*qj/r)*exp(-r/Dlength)
        DF = -Eahyd*(1/Dlength + 1/r)/r

      return
      end subroutine RNA_debye


!-----------------------------------------------------------------------
!>     @details
!>     This function takes care of the Base-Base interaction,
!>     including hydrogen bonding, stacking and cooperativity
!>     the 3 last atoms of each bases are used.
!>
!>     I is the last atom's index for the first base
!>     (so B1 for A and G, CY for C and U)
!>     J is the central atom's index for base 2
!>     TI and TJ are the bases' types
!>
!>     X is the system's coordinates vector
!>     F is the system's force vector
!>
!>     EHHB is the hydrogen bonding energy
!>
!>
!>     I-2 == I-1 == I  - - -  J+1 == J == J-1
!>
!> F : FI3    Fi2   Fi1        Fj1   FJ2   Fj3
!>
!-----------------------------------------------------------------------
      subroutine RNA_BB(I, TI, J, TJ, epshb, X, F, EHHB, hbexist)
      implicit none

      integer, intent(in) :: I, J, TI, TJ
      double precision, intent(in) :: epshb, X(*)
      double precision, intent(inout) :: F(*)
      logical, intent(out) :: hbexist

      integer idx, id
      double precision Ehhb, Enp1, Enp2, Etemp
      double precision, dimension(3,3) :: Ftemp
      double precision, dimension(3,3) :: Fi, Fj,  Fhb_i, Fhb_j, Fnp1_i, Fnp1_j, Fnp2_i, Fnp2_j
      double precision Ftemp_o(3)

      Fi = 0
      Fj = 0
      Fhb_i = 0
      Fhb_j = 0
 
      call RNA_hbnew( I, TI, J, TJ, epshb, X, EHHB, hbexist, Fhb_i, Fhb_j )
!     $            Fhb_i(:,2),Fhb_i(:,1),Fhb_j(:,1),Fhb_j(:,2))

      if (.not. hbexist) then
        return
      endif

      Enp1 = 0
      Fnp1_i = 0
      Fnp1_j = 0
      Enp2 = 0
      Fnp2_i = 0
      Fnp2_j = 0

!!      do idx = 0,5
      do idx = 0,2
        Etemp = 0
        Ftemp = 0
        Ftemp_o = 0
        call RNA_NewPlanev( I-2, I-1, I, J-idx+1 , X, Etemp, Ftemp(:,3), Ftemp(:,2), Ftemp(:,1), Ftemp_o )
        Fnp1_i = Fnp1_i + Ftemp
!        Fnp1_i(:,4:6) = Fnp1_i(:,4:6) + Ftemp
        Fnp1_j(:,idx+1) = Fnp1_j(:,idx+1) + Ftemp_o
        Enp1 = Enp1 + Etemp

        Etemp = 0
        Ftemp = 0
        Ftemp_o = 0
        call RNA_NewPlanev(J-1,J,J+1,I-idx , X, Etemp, Ftemp(:,3),Ftemp(:,2),Ftemp(:,1),Ftemp_o)
        Fnp2_j = Fnp2_j + Ftemp
!        Fnp2_i(:,4:6) = Fnp2_i(:,4:6) + Ftemp
        Fnp2_i(:,idx+1) = Fnp2_i(:,idx+1) + Ftemp_o
        Enp2 = Enp2 + Etemp
      enddo


!      if (ehhb .le. -2) then
!!        print *, Ehhb, Enp1, Enp2
!      endif
!
!
      Fi = Fhb_i*Enp1*Enp2 + Ehhb*Fnp1_i*Enp2 + Ehhb*Enp1*Fnp2_i
      Fj = Fhb_j*Enp1*Enp2 + Ehhb*Fnp1_j*Enp2 + Ehhb*Enp1*Fnp2_j
      Ehhb = Ehhb * Enp1 * Enp2

!      Fi = Fhb_i*(Enp1+Enp2)+ Ehhb*(Fnp1_i+Fnp2_i)
!      Fj = Fhb_j*(Enp1+Enp2)+ Ehhb*(Fnp1_j+Fnp2_j)
!      Ehhb = Ehhb * (Enp1 + Enp2)

!      Fi = Fhb_i
!      Fj = Fhb_j

!      F(I*3-2:I*3  )=F(I*3-2:I*3  ) + Fi(:,1)
!      F(I*3-5:I*3-3)=F(I*3-5:I*3-3) + Fi(:,2)
!      F(I*3-8:I*3-6)=F(I*3-8:I*3-6) + Fi(:,3)
!
!      F(J*3+1:J*3+3)=F(J*3+1:J*3+3) + Fj(:,1)
!      F(J*3-2:J*3  )=F(J*3-2:J*3  ) + Fj(:,2)
!      F(J*3-5:J*3-3)=F(J*3-5:J*3-3) + Fj(:,3)

      do idx = 1,3
        id = I - idx + 1
        F(id*3-2:id*3  )=F(id*3-2:id*3  ) + Fi(:,idx)

        id = J - idx + 2
        F(id*3-2:id*3  )=F(id*3-2:id*3  ) + Fj(:,idx)
      enddo

      end subroutine RNA_BB

!-----------------------------------------------------------------------
!>  @brief
!>  This routine calculates the h-bond energies and forces between two bases.
!
!>  @details
!>  hbexist is a boolean whose value will depend on the presence of an h-bond
!>
!>    Diagram:
!>
!>      va    ua           ub    vb
!>   a1 -- a2 -- a3 - - b3 -- b2 -- b1
!>              anga   angb
!>
!-----------------------------------------------------------------------
      subroutine RNA_hbnew( idxa, tya, idxb, tyb, epshb, X, EHHB, hbexist, fa, fb )
!     &                   FI,FJ,FK,FL)

      use geometric_corrections
      use RNAHBparams
      implicit none
      real(8) pbc_mic

      ! index of last particle in bases a,b and types of bases a,b
      integer, intent(in) :: idxa, tya, idxb, tyb
      double precision, intent(in) :: epshb, X(*)
      logical, intent(out) :: hbexist
      double precision, intent(out) :: EHHB, fa(3,3), fb(3,3)
!      double precision, intent(out) :: EHHB, fi(3), fj(3), fk(3), fl(3)

      logical qbug
      common/debug/qbug
      logical periodicBC,CM
      common/PBC_R/periodicBC,CM
      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                  rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                  rcut2_4b_out, rcut2_4b_in, &
                  rcut2_lj_out, rcut2_lj_in
      double precision rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                  rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                  rcut2_4b_out, rcut2_4b_in, &
                  rcut2_lj_out, rcut2_lj_in
 
      double precision, dimension(3) :: a1, a2, a3, b1, b2, b3
      double precision sighb, d2, Ehha, Ehb, dEhb(3), Vangl, y, &
      p, cosa, alpa, cosb, alpb, sina, sinb
      double precision anga, angb, &
      rba(3), dba, rba0(3), ralpa(3), ralpb(3), &
      ua(3), dua, ua0(3), va(3), dva, va0(3), &
      ub(3), dub, ub0(3), vb(3), dvb, vb0(3), &
      na(3), dna, na0(3), ma(3), dma, ma0(3), ra(3), dra, ra0(3), &
      nb(3), dnb, nb0(3), mb(3), dmb, mb0(3), rb(3), drb, rb0(3)

!      double precision, dimension(:), pointer :: dREF, alpam, alpbm, s
!     integer, allocatable, dimension(:) :: s
      integer par
      double precision str
      double precision pi
      data pi/3.141592653589793d+00/

      p = score_RNA(12)
      y = score_RNA(13)

      EHHB = 0
      fa = 0
      fb = 0
      hbexist = .false.


      a1 = x(idxa*3-8:idxa*3-6)
      a2 = x(idxa*3-5:idxa*3-3)
      a3 = x(idxa*3-2:idxa*3)

      b1 = x(idxb*3-5:idxb*3-3)
      b2 = x(idxb*3-2:idxb*3)
      b3 = x(idxb*3+1:idxb*3+3)


      ua = a3 - a2
      va = a1 - a2
      ub = b3 - b2
      vb = b1 - b2
      rba = a3 - b3
      if (periodicBC) then
        rba = pbc_mic( rba )
      endif
      dba = euc_norm(rba)
      rba0 = rba / dba

      if(dba**2 >= rcut2_hb_mcmc_out) then
        return
      endif
!      hbexist = .true.

      dua = euc_norm(ua)
      dva = euc_norm(va)
      dub = euc_norm(ub)
      dvb = euc_norm(vb)
      ua0 = ua / dua
      va0 = va / dva
      ub0 = ub / dub
      vb0 = vb / dvb

      na = crossproduct(ua, va)
      nb = crossproduct(ub, vb)
      dna = euc_norm(na)
      dnb = euc_norm(nb)
      na0 = na / dna
      nb0 = nb / dnb

      ma = crossproduct(na, ua)
      mb = crossproduct(nb, ub)
      dma = euc_norm(ma)
      dmb = euc_norm(mb)
      ma0 = ma / dma
      mb0 = mb / dmb

      ra = -rba - na0*dot_product(-rba, na0)
      dra = euc_norm(ra)
      ra0 = ra / dra
      rb = rba - nb0*dot_product(rba, nb0)
      drb = euc_norm(rb)
      rb0 = rb / drb

      cosa = dot_product(ra0, ua0)
      sina = dot_product(ra0, ma0)
      cosb = dot_product(rb0, ub0)
      sinb = dot_product(rb0, mb0)

!      print *, na, ma
!      print *, cosa, sina
!      print *, nb, mb
!      print *, cosb, sinb

      !	 selection of parameter for each base pair
      do par = 1, Nparam(tya,tyb)
        SIGHB = dREF(par,tya,tyb)
        alpa = alpam(par,tya,tyb)+pi
        alpb = alpbm(par,tya,tyb)+pi
        str = s(par,tya,tyb)

        d2 = (dba - SIGHB)/y
        Ehha = -epshb * str * dexp(-d2**2)

        if(ehha .ge. -1e-8) then
          cycle
        endif
 
 
        anga = cosa*DCOS(alpa) + sina*DSIN(alpa)
        ralpa = DCOS(alpa)*ua0 + DSIN(alpa)*ma0
        angb = cosb*DCOS(alpb) + sinb*DSIN(alpb)
        ralpb = DCOS(alpb)*ub0 + DSIN(alpb)*mb0

        if(anga < 0.0d0 .or. angb < 0.0d0) then
          cycle
        endif
!      hbexist = .true.

        Vangl = anga**p*angb**p

        Ehb = Ehha*Vangl
        dEhb = -Ehb * 2*d2/y * rba0

        if (Ehb .ge. -1e-7) then
          cycle
        endif

        Ehhb = Ehhb + Ehb

        fa(:,3) = fa(:,3)        - Ehb*(p/anga)* (&
!                   d anga / d ma0
          DSIN(alpa)*(crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)))/dma &
!                   d anga / d ra0
          -(crossproduct(ralpa-ra0*anga, ua)*dot_product(-rba, na0)/dna)/dra)
        fa(:,2) = fa(:,2)        - Ehb*(p/anga)* (&
!                   d anga / d ua0
          -DCOS(alpa)*(ra0- ua0*cosa)/dua &
!                   d anga / d ma0
          -DSIN(alpa)*(crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)))/dma &
-DSIN(alpa)*(crossproduct(va, crossproduct(ua, ra0-sina*ma0))+crossproduct(ra0-sina*ma0, na))/dma &
!                   d anga / d ra0
          -(crossproduct(ua-va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)
        fa(:,1) = fa(:,1) - dEhb - Ehb*(p/anga)* (&
!                   d anga / d ua0
          DCOS(alpa)*(ra0- ua0*cosa)/dua +&
!                   d anga / d ma0
DSIN(alpa)*(crossproduct(va, crossproduct(ua, ra0-sina*ma0))+crossproduct(ra0-sina*ma0, na))/dma &
!                   d anga / d ra0
          -(ralpa-ra0*anga+crossproduct(va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)&
!                   d angb / d rb0
          - Ehb*(p/angb)*(ralpb-rb0*angb)/drb
        fb(:,1) = fb(:,1) + dEhb - Ehb*(p/angb)* (&
!                   d angb / d ub0
          DCOS(alpb)*(rb0- ub0*cosb)/dub +&
!                   d angb / d mb0
DSIN(alpb)*(crossproduct(vb, crossproduct(ub, rb0-sinb*mb0))+crossproduct(rb0-sinb*mb0, nb))/dmb &
!                   d angb / d rb0
          -(ralpb-rb0*angb+crossproduct(vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)&
!                   d anga / d ra0
          - Ehb*(p/anga)*(ralpa-ra0*anga)/dra
        fb(:,2) = fb(:,2)        - Ehb*(p/angb)* (&
!                   d angb / d ub0
          -DCOS(alpb)*(rb0- ub0*cosb)/dub &
!                   d angb / d mb0
          -DSIN(alpb)*(crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)))/dmb &
-DSIN(alpb)*(crossproduct(vb, crossproduct(ub, rb0-sinb*mb0))+crossproduct(rb0-sinb*mb0, nb))/dmb &
!                   d angb / d rb0
          -(crossproduct(ub-vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)
        fb(:,3) = fb(:,3)        - Ehb*(p/angb)* (&
!                   d angb / d mb0
          DSIN(alpb)*(crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)))/dmb &
!                   d angb / d rb0
          -(crossproduct(ralpb-rb0*angb, ub)*dot_product(rba, nb0)/dnb)/drb)

        !c  ---  DEBUG
!       if (qbug) then
!            if (abs(EHB) .ge. 0.2 .or. (i .eq. 32 .and. j .eq. 50)) then
!            write(78,"(4i4, 10f8.3)") i,j,ti,tj,dkj,Ehb,Ehhb,Vangl,&
!              cosa, sina, cosb, sinb
!                 DACOS(cosa)*180/3.14,DACOS(cosb)*180/3.14,&
!                alpa*180/3.14,alpb*180/3.14
!          endif
!       endif

      enddo

      if(ehhb .ge. -1e-7) then      !! TEST HB EXISTANCE 06-04-2012
        return
      endif
      hbexist = .true.

      

      end subroutine RNA_hbnew


!-----------------------------------------------------------------------
      subroutine RNA_HBarrier(I,J,K,X,F,EHBR,ti,tj,tk)

      implicit none

      integer I, J, K, ti, tj, tk
      double precision X(*), F(*), EHBR, DRA,DRB

      double precision ra(3),rb(3),DHOA,DHOB,R1,R2
      double precision dfNI(3),dfNK(3)

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      double precision eta, A
      logical qbug
      common/debug/qbug

      eta = score_RNA(16)
      A   = score_RNA(17)

      ra = x((i*3-2):i*3) - x((j*3-2):(j*3))
      DHOA = dsqrt(dot_product(ra,ra))
 
      rb = x((k*3-2):k*3) - x((j*3-2):(j*3))
      DHOB = dsqrt(dot_product(rb,rb))

      R1 = 4.8
      if(ti .eq. 4 .and. tj .eq. 1) then
         R1 = 5.2
      endif
      if(tj .eq. 4 .and. ti .eq. 1) then
         R1 = 5.2
      endif
      if(ti .eq. 4 .and. tj .eq. 2) then
         R1 = 5.0
      endif
      if(tj .eq. 4 .and. ti .eq. 2) then
         R1 = 5.0
      endif
      R2 = 4.8
      if(tk .eq. 4 .and. tj .eq. 1) then
         R2 = 5.2
      endif
      if(tj .eq. 4 .and. tk .eq. 1) then
         R2 = 5.2
      endif
      if(tk .eq. 4 .and. tj .eq. 2) then
         R2 = 5.0
      endif
      if(tj .eq. 4 .and. tk .eq. 2) then
         R2 = 5.0
      endif

!      EHBR = A*exp(-eta*((DHOA-R1)**2+(DHOB-R2)**2+(DHOA-DHOB)**2))
      EHBR = A*exp(-eta*((DHOA-DHOB - (R1-R2))**2))
!      DRA = -2*eta*EHBR*(2*DHOA-DHOB-R1)
      DRA = -2*eta*EHBR*(DHOA-DHOB - (R1-R2))
!      DRB = -2*eta*EHBR*(DHOB-DHOA-R2)
      DRB = -2*eta*EHBR*(DHOA-DHOB - (R1-R2))
      
      if(qbug .and. (ehbr .ge. A*0.1 .or. ehbr .lt. 0.0)) then
         write(94,"(3i4, 6f8.3, f15.3)") I,J,K,dhoa,R1,dhob,R2, &
          exp(-eta*((DHOA-R1)**2)), exp(-eta*((DHOB-R2)**2)),EHBR
      endif


      dfNI = DRA*ra/DHOA

!      dfNK = DRB*rb/DHOB
      dfNK = -DRB*rb/DHOB

      F((I*3-2):I*3) = F((I*3-2):I*3) - dfNI
      F((J*3-2):J*3) = F((J*3-2):J*3) + dfNK + dfNI
      F((K*3-2):K*3) = F((K*3-2):K*3) - dfNK

      RETURN
      end subroutine RNA_HBarrier


!-------------------------------------------------------------------------
!>    @brief
!>     Computes the distance between one point and the plane defined by 3 other points.
!>     distance(l, plane(i,j,k))
!-------------------------------------------------------------------------
      subroutine RNA_NewPlanev(I,J,K,L,X,Enewpl,FI,FJ,FK,FL)

      use geometric_corrections

      implicit none

      integer I,J,K,L
      double precision X(*),Enewpl
      double precision, dimension(3) :: FI,FJ,FK,FL

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)

      logical qbug
      common/debug/qbug

      double precision, dimension(3) :: ri, rj, rk, rl, v1, v2, normal, t, dndq
      double precision dist, nnorm, dedd

      ri = X(i*3-2:i*3)
      rj = X(j*3-2:j*3)
      rk = X(k*3-2:k*3)
      rl = X(l*3-2:l*3)

      v1 = ri - rj
      v2 = rk - rj
      t = rl - rj

      normal = crossproduct(v1, v2)
      nnorm = dsqrt(dot_product(normal,normal))

      dist = dot_product(normal/nnorm,t)

      Enewpl = score_RNA(15)*exp(-(dist/score_RNA(14))**2)/6.0d0
      dedd = +2*(dist/score_RNA(14)**2)*Enewpl                           ! FORCE !!!

      FL = dedd*normal/nnorm

!     dn / d ri
      FI(1:3) = dedd*(crossproduct(v2, t) - dot_product(normal, t)*crossproduct(v2,normal)/nnorm**2)/nnorm

!     dn / d rk
      FK(1:3) = dedd*(crossproduct(t, v1) - dot_product(normal, t)*crossproduct(normal,v1)/nnorm**2)/nnorm

!     dn / d rj
      dndq = v2 - v1
      FJ(1:3) = dedd*(crossproduct(t, dndq)-normal - dot_product(normal, t)*crossproduct(normal,dndq)/nnorm**2)/nnorm

      end subroutine RNA_NewPlanev

!-----------------------------------------------------------------------------------
      subroutine RNA_Stackv(I,J,F,X,Estk)
      use geometric_corrections

!       i-2    i   j-2    j
!         \   /      \   /
!          \ /        \ /
!          i-1        j-1

      implicit none

      integer, intent(in) :: I, J
      double precision, intent(in) :: X(*)
      double precision, intent(inout) :: F(*)
      double precision, intent(out) :: Estk

      double precision a(3), b(3), c(3), d(3), r(3)
      double precision axb(3), cxd(3), Da(3)
      double precision axb0(3), cxd0(3), r0(3)
      double precision r1,r2,DotP,Dvr(3), Dvrij(3), SK, dot1, dot2
      double precision VA,VC, dot1w, dot2w, dot1wd, dot2wd, dotw, dotwd

      double precision eq, wid
      ! used in energy function cos(x)**cosp
      integer cosp

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      
      eq = score_RNA(8)
      wid = score_RNA(9)
      cosp = 2

      SK = score_RNA(7)

      a = X(i*3-8:i*3-6) - X(i*3-5:i*3-3)   !! vector a : I-1 -> I-2
      b = X(i*3-2:i*3  ) - X(i*3-5:i*3-3)    !! vector b : I-1 -> I

      c = X(j*3-8:j*3-6) - X(j*3-5:j*3-3)    !! vector c : J-1 -> J-2
      d = X(j*3-2:j*3  ) - X(j*3-5:j*3-3)    !! vector d : J-1 -> J

      axb = crossproduct(a, b)
      cxd = crossproduct(c, d)
      VA = euc_norm(axb)
      VC = euc_norm(cxd)
      axb0 = axb/VA
      cxd0 = cxd/VC

!      r = (X(i*3-2:i*3)+X(i*3-5:i*3-3)
!     $    -X(j*3-2:j*3)-X(j*3-5:j*3-3))/2
      r = (X(i*3-2:i*3)+X(i*3-5:i*3-3)+X(i*3-8:i*3-6) - X(j*3-2:j*3)-X(j*3-5:j*3-3)-X(j*3-8:j*3-6))/3
      r1 = euc_norm(r)
      r0 = r/r1
      
      DotP = dot_product(axb0, cxd0)
      dotw = 1 - (1 -DotP**2)**2
      dotwd = 2*(1-DotP**2) / (2-dotP**2)
      dot1 = dot_product(r0, axb0)
      dot1w = 1 - (1 -dot1**2)**2
      dot1wd = 2*(1-dot1**2) / (2-dot1**2)
      dot2 = dot_product(r0, cxd0)
      dot2w = 1 - (1 -dot2**2)**2
      dot2wd = 2*(1-dot2**2) / (2-dot2**2)

      r2 = (r1-eq)/wid
!      Estk = -SK * DotP**cosp * dexp(-r2**2)
!      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1**cosp * dot2**cosp
      ! Bypass derivatives calculation if E is very small
      ! This also prevents unstabilities arising from cos(x) ~= 0
      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1w * dot2w
      Estk = -SK * dotw * dexp(-r2**2) * dot1w * dot2w
      if (Estk > -1e-6 .or. abs(dot1) < 1e-6 .or. abs(dot2) < 1e-6 .or. abs(dotP) < 1e-6) then
        Estk = 0
        return
      endif
!      Estk = -SK * DotP**4 * (eq/r1)**6
!      Estk = -SK * DotP**4 * dexp(-3*(r1-3))
      Dvr = -Estk * 1/3 * 2*r2**1/wid * r/r1
!      Dvr = -Estk * 1/2 * 6/r1  *r/r1
!      Dvr = -Estk * 1/3 * (-3*(r1-3))  * -3*r/r1
!      Dvrij = Estk*cosp * 1/(3*r1) * (axb0/dot1 + cxd0/dot2 - 2*r0)
      Dvrij = Estk*cosp * 1/(3*r1) * (dot1wd*axb0/dot1 + dot2wd*cxd0/dot2 - r0*(dot1wd+dot2wd))
      Dvr = Dvr + Dvrij

!------- Derivatives on the 6 particles   -- 

!      Da = (cxd0/DotP - axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + r0/dot1 - 2*axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + dot1wd*r0/dot1 - axb0*(1+dot1wd))*cosp*Estk/VA
      Da = (dotwd*cxd0/DotP + dot1wd*r0/dot1 - axb0*(dotwd+dot1wd))*2*Estk/VA
!      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - crossproduct(b, Da)
      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - Dvr - crossproduct(b, Da)
      F(i*3-5:i*3-3) = F(i*3-5:i*3-3) - Dvr - crossproduct(a-b, Da)
      F(i*3-2:i*3  ) = F(i*3-2:i*3  ) - Dvr - crossproduct(Da, a)


!      Da = (axb0/DotP - cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + r0/dot2 - 2*cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + dot2wd*r0/dot2 - cxd0*(1+dot2wd))*cosp*Estk/VC
      Da = (dotwd*axb0/DotP + dot2wd*r0/dot2 - cxd0*(dotwd+dot2wd))*2*Estk/VC
!      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) - crossproduct(d, Da)
      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) + Dvr - crossproduct(d, Da)
      F(j*3-5:j*3-3) = F(j*3-5:j*3-3) + Dvr - crossproduct(c-d, Da)
      F(j*3-2:j*3  ) = F(j*3-2:j*3  ) + Dvr - crossproduct(Da, c)


      RETURN
      end subroutine RNA_Stackv
!-----------------------------------------------------------------------
      subroutine RNA_Coop(I,J,K,L,X,F,ECOOP,ti,tj,tk,tl)

!             I ------ J
!             |        |
!             |        |
!             K ------ L

      implicit none

      double precision X(*), F(*), ecoop
      integer i, j, k, l
      integer ti, tj, tk, tl

      double precision R1, R2, R0
      double precision, dimension(3) :: H1,H2,D1,D2
      double precision DE, DH1,DH2,DD1,DD2
      double precision VHB1,VHB2,CP,CM
      double precision DHB1,DHB2,DCP,DCM
      double precision, dimension(3) :: drh1,drh2, drd1,drd2

      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                  rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                  rcut2_4b_out, rcut2_4b_in, &
                  rcut2_lj_out, rcut2_lj_in
      double precision rcut2_caca_scsc_out, rcut2_caca_scsc_in, &
                  rcut2_hb_mcmc_out, rcut2_hb_mcmc_in, &
                  rcut2_4b_out, rcut2_4b_in, &
                  rcut2_lj_out, rcut2_lj_in

      logical qbug
      common/debug/qbug

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      double precision gam, delta, lambda, A

      gam = 4.0d0
      delta = 0.5d0
      lambda = 5.0d0  !! Delta/4
      A   = score_RNA(11)

      ECOOP = 0

      ! distance HB I---K
      D1 = x((i*3-2):i*3) - x((k*3-2):k*3)
      DD1 = dsqrt(dot_product(D1,D1))
      if (DD1**2 .gt. rcut2_4b_out) then
        return
      endif



      R0 = 5.0
      R1 = 4.8
      if((ti .eq. 4 .and. tj .eq. 1) .or. &
          (tj .eq. 4 .and. ti .eq. 1) .or. &
          (ti .eq. 3 .and. tj .eq. 2) .or. &
          (tj .eq. 3 .and. ti .eq. 2)) then
         R1 = 5.2
      endif

      R2 = 4.8
      if((tk .eq. 4 .and. tl .eq. 1) .or. &
          (tl .eq. 4 .and. tk .eq. 1) .or. &
          (tk .eq. 3 .and. tl .eq. 2) .or. &
          (tl .eq. 3 .and. tk .eq. 2)) then
         R2 = 5.2
      endif

      ! distance HB I---J
      H1 = x((i*3-2):i*3) - x((j*3-2):j*3)
      DH1 = dsqrt(dot_product(H1,H1))

      ! distance HB K---L
      H2 = x((k*3-2):k*3) - x((l*3-2):l*3)
      DH2 = dsqrt(dot_product(H2,H2))

      ! distance HB J---L
      D2 = x((j*3-2):j*3) - x((l*3-2):l*3)
      DD2 = dsqrt(dot_product(D2,D2))

!      CP = 1.0d0
!      DCP = 0.0d0

      VHB1 = dexp(-gam*(DH1-R1)**2)
      VHB2 = dexp(-gam*(DH2-R2)**2)
      CP = dexp(-delta*(DD1-DD2)**2)
!      CM = exp(-lambda*(DD1+DD2-2*R0)*(DD1+DD2-2*R0))
 
      if(DD1+DD2 .le. 2*R0) then
         CM = 1.0d0
         DCM = 0.0d0
      else
         CM = dexp(-lambda*(DD1+DD2-2*R0)**2)
         DCM = lambda*(DD1+DD2-2*R0)
      endif

      DE = 2*A*VHB1*VHB2*CP*CM

      DCM = DE*DCM
      DCP = DE*delta*(DD1-DD2)
      DHB1 = DE*gam*(DH1-R1)
      DHB2 = DE*gam*(DH2-R2)

      ECOOP = -DE/2

      if(qbug) then
      if(DH1 .le. R1+1 .and. DH2 .le. R2+1)then
      if(ECOOP .le. -0.001) then
          write(90,*) I,J,K,L
         write(90,"(4i4, 11f12.3)") I, J, K, L, ECOOP, DH1, DH2, R1, R2, DD1, DD2, VHB1, VHB2, CP, CM
      endif
      endif
      endif

!      DCP = 2*delta*CP*(DD1-DD2)*A*VHB1*VHB2*CM
!      DCM = 2*lambda*CM*(DD1+DD2-2*R0)*A*VHB1*VHB2*CP
    
      drh1 = DHB1*H1/DH1
  
      drh2 = DHB2*H2/DH2

      drd1 = (DCM+DCP)*D1/DD1
 
      drd2 = (DCM-DCP)*D2/DD2    

      F(I*3-2:I*3) = F(I*3-2:I*3) - drh1 - drd1

      F(J*3-2:J*3) = F(J*3-2:J*3) + drh1 - drd2

      F(K*3-2:K*3) = F(K*3-2:K*3) - drh2 + drd1

      F(L*3-2:L*3) = F(L*3-2:L*3) + drh2 + drd2

      RETURN
      end subroutine RNA_Coop



