      SUBROUTINE EBOND(NBON,IB,JB,ICB,X,F,EBON,RK,REQ)
C
C     THE POTENTIAL IS EXPRESSED BY: RK * (RIJ-REQ)**2 
C
      implicit double precision (a-h,o-z)

      real*8 pbc_mic

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      DIMENSION IB(*),JB(*),ICB(*),X(*),F(*),RK(*),REQ(*)
      DATA ZERO/0.0d0/

c     open(unit=7,file="beta32.ebond",status="unknown")
      EBON = ZERO 

      DO JN = 1,NBON
        I3 = IB(JN)
        J3 = JB(JN)
        XIJ = X(I3+1)-X(J3+1)
        YIJ = X(I3+2)-X(J3+2)
        ZIJ = X(I3+3)-X(J3+3)
        
        if(periodicBC)then
           XIJ = pbc_mic(xij)   !; print *, "xij", xij
           YIJ = pbc_mic(yij)   !; print *, "yij", yij
           ZIJ = pbc_mic(zij)   !; print *, "zij", zij
        endif
        
        RIJ0 = XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
        RIJ = SQRT(RIJ0)
        IC = ICB(JN)
        RIJ0 = RIJ
        DA = RIJ-REQ(IC)
        DF = RK(IC)*DA*score(267)
        ENERB = DF*DA
c     if (enerb .ge. 0.05) then
c     write(7,1200) I3,J3,RIJ,REQ(IC),ENERB,RK(IC)
c     1200 format(' I3 ',i4,' J3 ',i4, ' RIJ ',f8.3,' REQ ',f8.3,
c     1 ' ENERB ',f8.3, ' force ',f8.3)
c     endif
        DF = (DF+DF)/RIJ
        XA = DF*XIJ 
        YA = DF*YIJ 
        ZA = DF*ZIJ 
        F(I3+1) = F(I3+1)-XA
        F(I3+2) = F(I3+2)-YA
        F(I3+3) = F(I3+3)-ZA
        F(J3+1) = F(J3+1)+XA
        F(J3+2) = F(J3+2)+YA
        F(J3+3) = F(J3+3)+ZA
        EBON = EBON + ENERB
      ENDDO 
c     close(7)
c     stop
      RETURN
      END
C-------------------------------------------------------
      SUBROUTINE ETHETA(MAXTT,NTHETH,IT,JT,KT,ICT,X,F,ETHH,TK,TEQ)
C
C     THE POTENTIAL IS EXPRESSED BY: TK * (ANG(I,J,K)-TEQ)**2
C
      implicit double precision (a-h,o-z)

      real*8 pbc_mic

      parameter (MAXPRE = 1500)    !! maximum number of residues
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles

      LOGICAL qbug
      COMMON/debug/qbug

      common/scor/score(272),score_RNA(17)
      
      common/PBC_R/periodicBC,CM
      logical periodicBC,CM
      
      DIMENSION IT(*),JT(*),KT(*),ICT(*),X(*),F(*),TK(*),TEQ(*)
      
      DIMENSION XIJ(MAXTH),YIJ(MAXTH),ZIJ(MAXTH),XKJ(MAXTH),YKJ(MAXTH)
      DIMENSION ZKJ(MAXTH),CST(MAXTH),ANT(MAXTH),RKJ(MAXTH),RIJ(MAXTH)
      DIMENSION RIK(MAXTH),EAW(MAXTH),DFW(MAXTH)

      DATA pt999 /0.9990d0/
      DATA ZERO/0.0d0/

c     open(unit=7,file="beta32.ebond",status="unknown")
      ETHH = ZERO 

      DO JN = 1,NTHETH
        I3 = IT(JN)
        J3 = JT(JN)
        K3 = KT(JN)
        XIJ(JN) = X(I3+1)-X(J3+1)
        YIJ(JN) = X(I3+2)-X(J3+2)
        ZIJ(JN) = X(I3+3)-X(J3+3)
        XKJ(JN) = X(K3+1)-X(J3+1)
        YKJ(JN) = X(K3+2)-X(J3+2)
        ZKJ(JN) = X(K3+3)-X(J3+3)

        if(periodicBC)then
           XIJ(JN) = pbc_mic(XIJ(JN))
           YIJ(JN) = pbc_mic(YIJ(JN))
           ZIJ(JN) = pbc_mic(ZIJ(JN))
           XKJ(JN) = pbc_mic(XKJ(JN))
           YKJ(JN) = pbc_mic(YKJ(JN))
           ZKJ(JN) = pbc_mic(ZKJ(JN))
        endif
      ENDDO
      
      DO JN = 1,NTHETH
        RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
        RKJ0 = XKJ(JN)*XKJ(JN)+YKJ(JN)*YKJ(JN)+ZKJ(JN)*ZKJ(JN)
        RIK0 = SQRT(RIJ0*RKJ0)
        CT0 = (XIJ(JN)*XKJ(JN)+YIJ(JN)*YKJ(JN)+ZIJ(JN)*ZKJ(JN))/RIK0
        CT1 = MAX(-pt999,CT0)
        CT2 = MIN(pt999,CT1)
        CST(JN) = CT2
        ANT(JN) = DACOS(CT2)

        RIJ(JN) = RIJ0
        RKJ(JN) = RKJ0
        RIK(JN) = RIK0
      ENDDO

C     ENERGY
      DO JN = 1,NTHETH
        IC = ICT(JN)
        ANT0 = ANT(JN)
        DA = ANT0-TEQ(IC)
        DF = TK(IC)*DA*score(268)
        EAW(JN) = DF*DA
        DFW(JN) = -(DF+DF)/DSIN(ANT0)
        
        ETHH = ETHH + EAW(JN)
      ENDDO 

C     FORCE 
      DO JN = 1,NTHETH
        I3 = IT(JN)
        J3 = JT(JN)
        K3 = KT(JN)
        ST = DFW(JN)
        STH = ST*CST(JN)
        CIK = ST/RIK(JN)
        CII = STH/RIJ(JN)
        CKK = STH/RKJ(JN)
        DT1 = CIK*XKJ(JN)-CII*XIJ(JN)
        DT2 = CIK*YKJ(JN)-CII*YIJ(JN)
        DT3 = CIK*ZKJ(JN)-CII*ZIJ(JN)
        DT7 = CIK*XIJ(JN)-CKK*XKJ(JN)
        DT8 = CIK*YIJ(JN)-CKK*YKJ(JN)
        DT9 = CIK*ZIJ(JN)-CKK*ZKJ(JN)
        DT4 = -DT1-DT7
        DT5 = -DT2-DT8
        DT6 = -DT3-DT9
        F(I3+1) = F(I3+1)-DT1
        F(I3+2) = F(I3+2)-DT2
        F(I3+3) = F(I3+3)-DT3
        F(J3+1) = F(J3+1)-DT4
        F(J3+2) = F(J3+2)-DT5
        F(J3+3) = F(J3+3)-DT6
        F(K3+1) = F(K3+1)-DT7
        F(K3+2) = F(K3+2)-DT8
        F(K3+3) = F(K3+3)-DT9
      ENDDO 
c     write(7,*) jn,ethh 

      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ETORS(lambda,MAXPPI,NPHI,IP,JP,KP,LP,ICP,CG,IAC,X,F,
     +     EP,ENBP,EELP,MPHI,ECN,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)

      implicit double precision (a-h,o-z)

      real*8 pbc_mic

      parameter (MAXPRE = 1500)  !! maximum number of residues
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXTH = MAXNAT*3) !! maximum number of bond angles
      parameter (MAXPHI = MAXNAT*4) !! maximum number of torsional angles
      parameter (MAXTTY = 50000) !!
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs
      parameter (scal = 1.0)  !! scale
      parameter (scal14 = (1.476)) ! - weigth 261 from bestvect 24JULY06
      parameter (scalt = (0.232)) !  - weigth 258 from bestvect 24JULY06 
c     parameter (scal14 = 1.253) ! scale 261 
c     parameter (scalt = 0.983) ! scale torsion 258 

      LOGICAL DIELD
      LOGICAL QTOR
      logical QDET


      common/scalingfactor/scaling_factor
      common/scor/score(272),score_RNA(17)
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC
      COMMON/REWARD/ EHHB1

      COMMON/INDHB/indh(MAXPRE)
      COMMON/VP22/VPNE(MAXPRE)
      COMMON/ANGLES/ CTPHI(MAXPHI),CT1PHI(MAXPHI),NUPHI
      COMMON/TORS/QTOR

      LOGICAL qbug
      COMMON/debug/qbug

      real*8 scaling_factor
      DIMENSION XIJ(MAXPHI),YIJ(MAXPHI),ZIJ(MAXPHI),XKJ(MAXPHI)
      DIMENSION ZKJ(MAXPHI),XKL(MAXPHI),YKL(MAXPHI),ZKL(MAXPHI)
      DIMENSION GX(MAXPHI),GY(MAXPHI),GZ(MAXPHI),CT(MAXPHI),CPHI(MAXPHI)

      DIMENSION Z2(MAXPHI),FXI(MAXPHI),FYI(MAXPHI),FZI(MAXPHI),
     + FXJ(MAXPHI),DX(MAXPHI),
     + FXK(MAXPHI),FYK(MAXPHI),FZK(MAXPHI),FXL(MAXPHI),FYL(MAXPHI),
     + YKJ(MAXPHI),DY(MAXPHI),DZ(MAXPHI),SPHI(MAXPHI),Z1(MAXPHI),
     + FYJ(MAXPHI),FZJ(MAXPHI),FZL(MAXPHI),DF(MAXPHI)

      DIMENSION XIL(MAXPHI),YIL(MAXPHI),ZIL(MAXPHI),RRW(MAXPHI)
      DIMENSION EPW(MAXPHI),ENW(MAXPHI),FMUL(MAXPHI),BI(MAXPHI)
      DIMENSION GMUL(10),EEW(MAXPHI),BK(MAXPHI)

      DIMENSION IP(*),JP(*),KP(*),LP(*),ICP(*),IPN(*)
      DIMENSION CG(*),IAC(*),X(*),F(*),CN1(*),CN2(*),PK(*),PN(*)
      DIMENSION GAMC(*),GAMS(*),FMN(*)
      DIMENSION CTPHIB(MAXPHI)

      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

      REAL*8 IDIEL
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     +     0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      DATA PI/3.141592653589793d+00/
      data zero,one,four,six,twelve/0.0d0,1.0d0,4.0d0,6.0d0,12.0d0/

      real*8 lambda
     

c     open(unit=64,file="beta32.edih",status="unknown")
      EPL = zero
      EP = zero
      ECNL = zero
      ENBPL = zero
      EELPL = zero
      EV1 = ZERO
      EV2 = ZERO
      EV3 = ZERO
      IRES = 0
      SCNB0 = one/SCNB
      SCEE0 = one/SCEE
      DIELD = IDIEL.LE.0

      DO JN = 1,MPHI
        I3 = IP(JN)
        J3 = JP(JN)
        K3T = KP(JN)
        L3T = LP(JN)
        K3 = IABS(K3T)
        L3 = IABS(L3T)
        XIJ(JN) = X(I3+1)-X(J3+1)
        YIJ(JN) = X(I3+2)-X(J3+2)
        ZIJ(JN) = X(I3+3)-X(J3+3)
        XKJ(JN) = X(K3+1)-X(J3+1)
        YKJ(JN) = X(K3+2)-X(J3+2)
        ZKJ(JN) = X(K3+3)-X(J3+3)
        XKL(JN) = X(K3+1)-X(L3+1)
        YKL(JN) = X(K3+2)-X(L3+2)
        ZKL(JN) = X(K3+3)-X(L3+3)

        if(periodicBC)then
           XIJ(JN) = pbc_mic( XIJ(JN) )
           YIJ(JN) = pbc_mic( YIJ(JN) )
           ZIJ(JN) = pbc_mic( ZIJ(JN) )
           XKJ(JN) = pbc_mic( XKJ(JN) )
           YKJ(JN) = pbc_mic( YKJ(JN) )
           ZKJ(JN) = pbc_mic( ZKJ(JN) )
           XKL(JN) = pbc_mic( XKL(JN) )
           YKL(JN) = pbc_mic( YKL(JN) )
           ZKL(JN) = pbc_mic( ZKL(JN) )
        endif
c     write(64,*) jn,i3,j3,k3,l3 
c     dx   ENDDO 
c     dx   DO JN = 1,MPHI
        DX(JN) = YIJ(JN)*ZKJ(JN)-ZIJ(JN)*YKJ(JN)
        DY(JN) = ZIJ(JN)*XKJ(JN)-XIJ(JN)*ZKJ(JN)
        DZ(JN) = XIJ(JN)*YKJ(JN)-YIJ(JN)*XKJ(JN)
        GX(JN) = ZKJ(JN)*YKL(JN)-YKJ(JN)*ZKL(JN)
        GY(JN) = XKJ(JN)*ZKL(JN)-ZKJ(JN)*XKL(JN)
        GZ(JN) = YKJ(JN)*XKL(JN)-XKJ(JN)*YKL(JN)
c     dx   ENDDO 
c     dx   DO JN = 1,MPHI
        BI(JN) = SQRT(DX(JN)*DX(JN)+DY(JN)*DY(JN)+DZ(JN)*DZ(JN)+TM24)
        BK(JN) = SQRT(GX(JN)*GX(JN)+GY(JN)*GY(JN)+GZ(JN)*GZ(JN)+TM24)
        CT(JN) = DX(JN)*GX(JN)+DY(JN)*GY(JN)+DZ(JN)*GZ(JN)
        
        ctphib(jn) = ct(jn)/(bi(jn)*bk(jn))
        if(ctphib(jn) .gt. 0.9999999d0) then
           ctphib(jn) = 1.0d0
        else if(ctphib(jn) .lt. -0.9999999d0) then
           ctphib(jn) = -1.0d0
        else  
           ctphib(jn) = -ctphib(jn)
        end if  
        
        ctphib(jn) = acos(ctphib(jn))*57.29577951d0
        sigphi = xij(jn)*gx(jn)+yij(jn)*gy(jn)+zij(jn)*gz(jn)
        if(sigphi .gt. 0.0d0) ctphib(jn) = 360.0d0-ctphib(jn)
        if(ctphib(jn) .gt. 180.) ctphib(jn) = ctphib(jn) - 360.0d0
c     dx   ENDDO 
c     dx   DO JN = 1,MPHI
        z10 = one/bi(jn)
        z20 = one/bk(jn)
        if (tenm3 .gt. bi(jn)) z10 = zero
        if (tenm3 .gt. bk(jn)) z20 = zero
        Z12 = Z10*Z20
        Z1(JN) = Z10
        Z2(JN) = Z20
        ftem = zero
        if (z12 .ne. zero) ftem = one
        FMUL(JN) = FTEM
        
        CT0 = MIN(one,CT(JN)*Z12)
        CT1 = MAX(-one,CT0)
        
        S = XKJ(JN)*(DZ(JN)*GY(JN)-DY(JN)*GZ(JN))+
     +       YKJ(JN)*(DX(JN)*GZ(JN)-DZ(JN)*GX(JN))+
     +       ZKJ(JN)*(DY(JN)*GX(JN)-DX(JN)*GY(JN))
        AP0 = DACOS(CT1)
        AP1 = PI-DSIGN(AP0,S)
        CT(JN) = AP1
        CPHI(JN) = DCOS(AP1)
        SPHI(JN) = DSIN(AP1)
c     dx   ENDDO 

C     ----- ENERGY AND THE DERIVATIVES WITH RESPECT TO
C     COSPHI -----

c     dx   DO JN = 1,MPHI
        IC = ICP(JN)
        INC = IPN(IC)
        CT0 = PN(IC)*CT(JN)
        COSNP = DCOS(CT0)
        SINNP = DSIN(CT0)
        EPW(JN) = (PK(IC)+COSNP*GAMC(IC)+SINNP*GAMS(IC))*FMUL(JN)
        DF0 = PN(IC)*(GAMC(IC)*SINNP-GAMS(IC)*COSNP)
        DUMS = SPHI(JN)+SIGN(TM24,SPHI(JN))
        DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*CPHI(JN))
        df1 = df0/dums
        if(tm06.gt.abs(dums)) df1 = dflim
        DF(JN) = DF1*FMUL(JN)
        
        EPW(JN) = EPW(JN)*scalt*scaling_factor  
        DF(JN) = DF(JN)*scalt*scaling_factor   
        
        I3 = IP(JN)
        J3 = JP(JN)
        K3T = KP(JN)
        L3T = LP(JN)
        K3 = IABS(K3T)
        L3 = IABS(L3T)
c     dx   ENDDO 
C     END ENERGY WITH RESPECT TO COSPHI
        
        
c     dx   DO JN = 1,MPHI
C     ----- DC = FIRST DER. OF COSPHI W/RESPECT
C     TO THE CARTESIAN DIFFERENCES T -----
        Z11 = Z1(JN)*Z1(JN)
        Z12 = Z1(JN)*Z2(JN)
        Z22 = Z2(JN)*Z2(JN)
        DC1 = -GX(JN)*Z12-CPHI(JN)*DX(JN)*Z11
        DC2 = -GY(JN)*Z12-CPHI(JN)*DY(JN)*Z11
        DC3 = -GZ(JN)*Z12-CPHI(JN)*DZ(JN)*Z11
        DC4 =  DX(JN)*Z12+CPHI(JN)*GX(JN)*Z22
        DC5 =  DY(JN)*Z12+CPHI(JN)*GY(JN)*Z22
        DC6 =  DZ(JN)*Z12+CPHI(JN)*GZ(JN)*Z22
C     ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
        DR1 = DF(JN)*(DC3*YKJ(JN) - DC2*ZKJ(JN))
        DR2 = DF(JN)*(DC1*ZKJ(JN) - DC3*XKJ(JN))
        DR3 = DF(JN)*(DC2*XKJ(JN) - DC1*YKJ(JN))
        DR4 = DF(JN)*(DC6*YKJ(JN) - DC5*ZKJ(JN))
        DR5 = DF(JN)*(DC4*ZKJ(JN) - DC6*XKJ(JN))
        DR6 = DF(JN)*(DC5*XKJ(JN) - DC4*YKJ(JN))
        DRX = DF(JN)*(-DC2*ZIJ(JN) + DC3*YIJ(JN) +
     +       DC5*ZKL(JN) - DC6*YKL(JN))
        DRY = DF(JN)*( DC1*ZIJ(JN) - DC3*XIJ(JN) -
     +       DC4*ZKL(JN) + DC6*XKL(JN))
        DRZ = DF(JN)*(-DC1*YIJ(JN) + DC2*XIJ(JN) +
     +       DC4*YKL(JN) - DC5*XKL(JN))
        FXI(JN) = - DR1
        FYI(JN) = - DR2
        FZI(JN) = - DR3
        FXJ(JN) = - DRX + DR1
        FYJ(JN) = - DRY + DR2
        FZJ(JN) = - DRZ + DR3
        FXK(JN) = + DRX + DR4
        FYK(JN) = + DRY + DR5
        FZK(JN) = + DRZ + DR6
        FXL(JN) = - DR4
        FYL(JN) = - DR5
        FZL(JN) = - DR6
c     dx   ENDDO 
C     ----- CALCULATE 1-4 NONBONDED CONTRIBUTIONS
c     dx   DO JN = 1,MPHI
        I3 = IP(JN)
        L3T = LP(JN)
        L3 = IABS(L3T)
        XIL(JN) = X(I3+1)-X(L3+1)
        YIL(JN) = X(I3+2)-X(L3+2)
        ZIL(JN) = X(I3+3)-X(L3+3)

        if(periodicBC)then
           XIL(JN) = pbc_mic(xil(jn))
           YIL(JN) = pbc_mic(yil(jn))
           ZIL(JN) = pbc_mic(zil(jn))
        endif
c     ENDDO
c     DO JN = 1,MPHI
        RRW(JN) = XIL(JN)*XIL(JN)+YIL(JN)*YIL(JN)+ZIL(JN)*ZIL(JN)
      ENDDO

c     dx   IF(DIELD) GO TO 720
      if (.not. DIELD) then
         DO 700 JN = 1,MPHI
           I3 = IP(JN)
           K3T = KP(JN)
           L3T = LP(JN)
           IC0 = ICP(JN)
           IDUMI = ISIGN(1,K3T)
           IDUML = ISIGN(1,L3T)
           KDIV = (2+IDUMI+IDUML)/4
           L3 = IABS(L3T)
c     dx       FMULN = FLOAT(KDIV)*FMN(IC0)
           FMULN = dble(kdiv)*FMN(IC0)
           II = (I3+3)/3
           JJ = (L3+3)/3
           IA1 = IAC(II)
           IA2 = IAC(JJ)
           IBIG = MAX0(IA1,IA2)
           ISML = MIN0(IA1,IA2)   
           IC = IBIG*(IBIG-1)/2+ISML
           R2 = FMULN/RRW(JN)
           R1 = SQRT(R2)
           G = CG(II)*CG(JJ)*R1
           EEW(JN) = G
           R6 = R2*R2*R2
           R12 = R6*R6
           F1 = CN1(IC)*R12*scal14*scaling_factor 
           F2 = CN2(IC)*R6*scal14*scaling_factor*lambda
           ENW(JN) = (F1-F2)
           DFN =((-twelve*F1+six*F2)*SCNB0-G*SCEE0)*R2
           XA = XIL(JN)*DFN
           YA = YIL(JN)*DFN 
           ZA = ZIL(JN)*DFN
           FXI(JN) = FXI(JN)-XA
           FYI(JN) = FYI(JN)-YA
           FZI(JN) = FZI(JN)-ZA
           FXL(JN) = FXL(JN)+XA
           FYL(JN) = FYL(JN)+YA
           FZL(JN) = FZL(JN)+ZA
 700     CONTINUE
c     dx   GO TO 760
      else
C     
C     ----- DISTANCE DEPENDENT DIELECTRIC -----
C     
 720     CONTINUE
         DO 740 JN = 1,MPHI
           I3 = IP(JN)
           K3T = KP(JN)
           L3T = LP(JN)
           IC0 = ICP(JN)
           IDUMI = ISIGN(1,K3T)
           IDUML = ISIGN(1,L3T)
           KDIV = (2+IDUMI+IDUML)/4
           L3 = IABS(L3T)
c     dx       FMULN = FLOAT(KDIV)*FMN(IC0)
           FMULN = dble(KDIV)*FMN(IC0)
           II = (I3+3)/3  
           JJ = (L3+3)/3
           IA1 = IAC(II)
           IA2 = IAC(JJ)
           IBIG = MAX0(IA1,IA2)
           ISML = MIN0(IA1,IA2)
           IC = IBIG*(IBIG-1)/2+ISML
           R2 = FMULN/RRW(JN)
c     print*,' R2 ',R2
           G = CG(II)*CG(JJ)*R2
!            write(*,*) 'CG', II
           EEW(JN) = G
           R6 = R2*R2*R2
           R12 = R6*R6
           F1 = CN1(IC)*R12*scal14*scaling_factor  
           F2 = CN2(IC)*R6*scal14*scaling_factor*lambda  
           ENW(JN) = F1-F2
           DFN =((-twelve*F1+six*F2)*SCNB0-(G+G)*SCEE0)*R2
           XA = XIL(JN)*DFN
           YA = YIL(JN)*DFN
           ZA = ZIL(JN)*DFN
           FXI(JN) = FXI(JN)-XA
           FYI(JN) = FYI(JN)-YA
           FZI(JN) = FZI(JN)-ZA
           FXL(JN) = FXL(JN)+XA
           FYL(JN) = FYL(JN)+YA
           FZL(JN) = FZL(JN)+ZA
 740     CONTINUE
 760     CONTINUE
      end if
C     ----- THE TOTAL FORCE VECTOR -----

      DO JN = 1,MPHI
        I3 = IP(JN)
        J3 = JP(JN)
        K3 = IABS(KP(JN))
        L3 = IABS(LP(JN))
C     
        F(I3+1) = F(I3+1) + (FXI(JN)) ! *scal)
        F(I3+2) = F(I3+2) + (FYI(JN)) ! *scal)
        F(I3+3) = F(I3+3) + (FZI(JN)) ! *scal)
        F(J3+1) = F(J3+1) + (FXJ(JN)) ! *scal) 
        F(J3+2) = F(J3+2) + (FYJ(JN)) ! *scal)
        F(J3+3) = F(J3+3) + (FZJ(JN)) ! *scal)
        F(K3+1) = F(K3+1) + (FXK(JN)) ! *scal)
        F(K3+2) = F(K3+2) + (FYK(JN)) ! *scal)
        F(K3+3) = F(K3+3) + (FZK(JN)) ! *scal)
        F(L3+1) = F(L3+1) + (FXL(JN)) ! *scal)
        F(L3+2) = F(L3+2) + (FYL(JN)) ! *scal)
        F(L3+3) = F(L3+3) + (FZL(JN)) ! *scal)
c     dx   ENDDO 
c     dx   do ksum = 1,MPHI
        enbpl = enbpl+enw(JN)   !! 1-4 nb
        eelpl = eelpl+eew(JN)   !! 1-4 elec
!         write(*,*) 'eelpl', epl 
        epl   = epl  +epw(JN)   !! torsions
      enddo 
      ENBP = ENBPL*SCNB0
      EELP = EELPL*SCEE0
      EP   = EPL
      ECN = ECNL
c     close(64)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ENBOND(lambda,NATOM,CG,IAC,X,F,CN1,CN2,EVDW,ELEC
     $     ,index_ip,n_ip,ip_control)     ! YC

      implicit double precision (a-h,o-z)

      real*8 pbc_mic

      logical dield
      parameter (MAXPRE = 1500)  !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6) !! maximum number of atoms
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2) !! max number of nonbonded-pairs 

      common/scor/score(272),score_RNA(17)
      COMMON/NPAIR/NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)
      COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC

      COMMON/EHB/EHBT

      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $     rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $     rcut2_4b_out, rcut2_4b_in,
     $     rcut2_lj_out, rcut2_lj_in

      logical periodicBC,CM
      common/PBC_R/periodicBC,CM
            
      REAL*8 IDIEL
      DIMENSION CG(*),IAC(*),X(*),F(*),CN1(*),CN2(*)
      data zero,one,four,six,ten,twelve/0.0d0,1.0d0,4.0d0,6.0d0,
     1     10.0d0,12.0d0/
      
      real*8 lambda
!      OPEN(UNIT=72,FILE="enbond.list",status="unknown")

c----- Test for salt bridge : removing specific interaction

      integer index_ip(500,3),n_ip
      logical ip_control

c----------------------------------------------------------


      Scnb0 = ONE 
      SCEE0 = ONE 
      EVDW = ZERO 
      ELEC = ZERO 
      EHBT = ZERO 
      
      DIELD = IDIEL.LE.0
      IF (.NOT. DIELD) THEN
       PRINT*,' NONBONDED POTENTIEL FOR INDEPENDENT-DISTANCE DIELECTRIC'
       PRINT*,' NOT DEFINED'
       STOP
      ENDIF

      rcut_2_lj_in=144
      rcut_2_lj_out=225

!       write(*,*) 'npair2', npair2
      DO JN = 1,NPAIR2



         G=ZERO
         ENW=ZERO
         F1=ZERO
         F2=ZERO
         DFN=ZERO
         
        II = IPAIR(JN)          !! (I3+3)/3  
        JJ = JPAIR(JN)          !! (L3+3)/3
        I3 = II*3;   J3 = JJ*3;   ZA = X(I3)-X(J3)
        I2 = I3 - 1; J2 = J3 - 1; YA = X(I2)-X(J2)
        I1 = I3 - 2; J1 = J3 - 2; XA = X(I1)-X(J1)

        if (periodicBC) then
           xa = pbc_mic( xa )
           ya = pbc_mic( ya )
           za = pbc_mic( za )
        endif
        DA = XA*XA+YA*YA+ZA*ZA
!         write(*,*) 'DA', DA

        if (DA < rcut2_lj_out) then
           IA1 = IAC(II)
           IA2 = IAC(JJ)
           IBIG = MAX0(IA1,IA2)
           ISML = MIN0(IA1,IA2)
           IC = IBIG*(IBIG-1)/2+ISML

           R2 = 1.0d0/DA; R6 = R2*R2*R2; R12 = R6*R6
           F1 = CN1(IC)*score(270)*R12
           F2 = CN2(IC)*score(270)*R6*lambda      !!!lambda
           ENW = F1-F2
!            write(*,*) 'EVDW-F1', enw, F1,F2  !! OK
!           write(72,1205) II,JJ,IC,CN1(IC),CN2(IC) 
! 1205 format('II',i4,'JJ',i4, 'IC', i4, ' CN1 ',f12.4,' CN2 ',f12.4)
!      ; EEW = G; ELEC = ELEC + EEW
c     dx BUG!          
c     dx       G = CG(II)*CG(JJ)*R2
c     dx
           G = CG(II)*CG(JJ)*sqrt(R2)
!            write(*,*) 'II', II
           DFN =((-twelve*F1+six*F2)-(G+G))*R2
           
           EEW = G
           ELEC = ELEC + EEW
!            write(*,*) 'ELEC-G', G

           if (DA > rcut2_lj_in) 
     $          call switch_cutoff(DA,enw,dfn,rcut2_lj_in,rcut2_lj_out)
           
           EVDW = EVDW + ENW
!            write(*,*) 'EVDW', evdw
           DXA = DFN*XA; F(I1) = F(I1) - DXA; F(J1) = F(J1) + DXA
           DYA = DFN*YA; F(I2) = F(I2) - DYA; F(J2) = F(J2) + DYA
           DZA = DFN*ZA; F(I3) = F(I3) - DZA; F(J3) = F(J3) + DZA
        end if
      ENDDO
!     $omp end parallel do
!      STOP
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE HYDROP(lambda,X,F,EHYDRO,index_ip,n_ip,ip_control)
c      SUBROUTINE HYDROP(lambda,X,F,EHYDRO)
C     
      implicit double precision (a-h,o-z)

      real*8 pbc_mic
C     
C     ----- ROUTINE TO CALCULATE THE HYDROPHOBIC/HYDROPHILIC FORCES -----
C     --  and THE H-BOND FORCES
C     -- the analytic form includes now the propensity of residues to
C     prefer alpha or beta states and the weights for all contributions
c     has been rescaled (to be published in 2006)

      parameter (MAXPRE = 1500)    !! maximum number of residues
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord
      parameter (MAXPNB = 3*MAXPRE*MAXPRE)!! max number of SC-SC interactions
      parameter (MAXBO  = MAXNAT)  !! maximum number of bonds
      parameter (MAXTH = MAXNAT*3)  !! maximum number of bond angles
      parameter (MAXPHI = MAXNAT*4)  !! maximum number of torsional angles
      parameter (MAXTTY = 50000)  !!
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! NEW 21 JANV05

      character*3 resp
      LOGICAL QDAT,QPRINT 
      integer ref_NHB

      double precision score, score_RNA
      common/scor/score(272),score_RNA(17)
      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $     NPHIA,NNB,NTYPES,
     $     MBONA,MTHETA,MPHIA

      COMMON/HYDRO/ rncoe(maxpnb),vamax(maxpnb),
     1     ni(maxpnb),nj(maxpnb),nstep,nb,
     2     ivi(maxpnb),ivj(maxpnb),
     $     epshb_mcmc(maxpnb)
      double precision rncoe, vamax, epshb_mcmc
      integer*8 ni, nj, ivi, ivj, nstep, nb

      common/propens/ialpha(MAXNAT),ibeta(MAXNAT),icoeff(MAXPAI),
     1     foal(20),walpha(20),fobe(20),wbeta(20),nfrag,
     2     lenfrag(MAXPRE),ichain(MAXNAT)
      common/propens2/walpha_foal(20), wbeta_fobe(20)

      common/ncall/LTER


      COMMON/PROPE/ resp(MAXPRE)
      COMMON/REWARD/ EHHB1

      COMMON/ENEHB/ EHHB3
      common/JANV06/NHB,ref_NHB,NHBH,NHBB

      DIMENSION X(*),F(*), EH11(MAXNAT),EPOU(MAXPHI) 
c     dx   dimension x(*),f(*)
c     dx   dimension eh11(maxpnb), epou(maxpnb)

c     dx   DIMENSION DFOHX(MAXNAT),DFOHY(MAXNAT),DFOHZ(MAXNAT),dfohr(maxnat)
      dimension dho2s(maxnat)
      dimension dfo(maxnat), hx(maxnat), hy(maxnat), hz(maxnat)
c     dx   DIMENSION WEISCSC(300),EL(MAXPNB),ELT(MAXPNB)
      common/cacascsc/ct0lj(maxpnb),ct2lj(maxpnb)
      
      DIMENSION IHCO(MAXPNB),JHCO(MAXPNB),RNC(MAXPNB),
     1     ini(maxpnb),inj(maxpnb) 

      DATA zero/0.0d0/
      DATA one/1.0d0/
      DATA two/2.0d0/
      DATA three/3.0d0/
      DATA four/4.0d0/
      DATA five/5.0d0/
      DATA six/6.0d0/
      DATA ten/10.0d0/
      DATA twelve/12.0d0/
      DATA ninety/90.0d0/
      data pt999 /0.9990d0/

      real*8 parah, parah1, parah2, parah3, parah4
      logical logic1, logic2, logic3, logic4
      logical logic5, logic6, logic7, logic8,logic9
      logical logic21, logic22, logic23, logic24, logic25, logic26
      logical logic31, logic32, logic33, logic34
      logical logic41, logic42, logic43
      logical logic51, logic52, logic53, logic54
      logical logic55, logic56, logic57, logic58, logic59
      logical logic1333, logic6039, logic6300
      logical logic1i, logic21i, logic51i, logic58i

      parameter (rad2deg=180.0d0/3.141592653d0)
      
      common/cutoffs/rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $               rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $               rcut2_4b_out, rcut2_4b_in,
     $               rcut2_lj_out, rcut2_lj_in
      logical periodicBC,CM
      common/PBC_R/periodicBC,CM
      
      logical lcont
      real*8  lambda     !! lambda, for hamiltonian replica exchange it scales H-bond attraction

c----- Test for salt bridge : removing specific interaction


      integer   index_ip(500,3),n_ip
      logical   ip_control

C--   SET-up SOME PARAMETERS
      QDAT = .FALSE.
      QPRINT = .FALSE.          !! .FALSE.

      NHBT = NHB
      NHB = 0
      NHBH = 0
      NHBB = 0
      ivnu = 0
      inbre = 0
      inbreb = 0
      ehhb = zero
      ehhb1 = zero
      ehhb3 = zero
      ehydro = zero
      ecaca = zero
      etest = zero
      EHB14 = ZERO
      EHB15P = ZERO
      EHBCOA = -1.25D0          !! -0.5D0 !! estimate for alpha-helices
      EHBCOAB = -2.0D0*score(266) !! estimate for anti-parallel beta-sheets 
      EHBCOBP = -2.0D0*score(266) !! estimate for parallel beta-sheets
      DHOREF = 1.8D0            !! 1.80 optimal distance hydrogen bond
      SIGHB = 1.8D0             !! 1.80 optimal distance hydrogen bond 
      SHB2 = SIGHB*SIGHB
      SHB5 = SHB2*SHB2*SIGHB
      SHB10 = SHB5*SHB5
      SHB12 = SHB10*SHB2

      WEICA = score(224)        !! 1.0d0 !! weight for CA-CA VdW
      WEIHB14 = score(222)      !! ONE  !! weight for H-bonds helix
      WEIHB15 = score(223)      !! ONE  !! weight for others intra and inter

      WCOA = score(265)         !! weight for H-bond coop in helix without propensities
      wcoa_EHBCOA = wcoa*EHBCOA


      epou = 0.0d0

      DO 100 I = 1, nb          !! all two-body interactions
        NI3 = ni(i)*3
        NI2 = NI3 - 1
        NI1 = NI3 - 2
        NI4 = NI3 - 5
        NI5 = NI3 - 4
        NI6 = NI3 - 3

        NJ3 = nj(i)*3
        NJ2 = NJ3 - 1
        NJ1 = NJ3 - 2
        NJ4 = NJ3 + 1
        NJ5 = NJ3 + 2
        NJ6 = NJ3 + 3

        XA = X(NI1)-X(NJ1)
        YA = X(NI2)-X(NJ2) 
        ZA = X(NI3)-X(NJ3)

        if(periodicBC)then
           xa = pbc_mic( xa )
           ya = pbc_mic( ya )
           za = pbc_mic( za )
        endif
        DA2 = XA*XA+YA*YA+ZA*ZA


C-----------------------------------
C---  THE CA-CA and Sc-Sc terms are calculated 
C-----------------------------------
        IF (icoeff(i) /= -2 .and. 
     $       rncoe(i)<=9.0d0 .and. DA2<=rcut2_caca_scsc_out) then 
           call ljcasc(lambda,da2,eahyd,df,rncoe(i),ct0lj(i),ct2lj(i),
     $          rcut2_caca_scsc_in,rcut2_caca_scsc_out,vamax(i))
     
c--- I remove contribution from LJ for Ion-Pair. Later more elegant exlusion
!            write(*,*) 'ehydro', eahyd
           IF(ip_control) THEN
              DO icont=1,n_ip
                 ip1=index_ip(icont,1)
                 ip2=index_ip(icont,2)
                 IF(ip1.eq.ni(i).AND.ip2.eq.nj(i)) THEN
                    DF=0
                    eahyd=0
                 ENDIF
                 
                 IF(ip1.eq.nj(i).AND.ip2.eq.ni(i)) THEN
                    DF=0
                    eahyd=0
                 ENDIF
                 
              ENDDO
           ENDIF

           ehydro = ehydro + eahyd

           DXA = DF*XA
           DYA = DF*YA
           DZA = DF*ZA
           F(NI1) = F(NI1) - DXA
           F(NI2) = F(NI2) - DYA
           F(NI3) = F(NI3) - DZA
           F(NJ1) = F(NJ1) + DXA
           F(NJ2) = F(NJ2) + DYA
           F(NJ3) = F(NJ3) + DZA
           
        ELSE
c-----------------------------------------
c     --  H-bonds from MAIN CHAIN to MAIN CHAIN
c-----------------------------------------
           XA = X(NI1)-X(NJ4)
           YA = X(NI2)-X(NJ5)
           ZA = X(NI3)-X(NJ6)

           if(periodicBC)then
              xa = pbc_mic( xa )
              ya = pbc_mic( ya )
              za = pbc_mic( za )
           endif
           
           DHO2 = XA*XA+YA*YA+ZA*ZA

           if (icoeff(i) /= -1 .and. 
     $          rncoe(i)>=25.0d0 .and. rncoe(i)<=29.0d0 
     $          .and. DHO2<=rcut2_hb_mcmc_out) then

c     ---     calculate angle N-H..O
              XIJ = XA
              YIJ = YA 
              ZIJ = ZA 
              XKJ = X(NJ1) - X(NJ4)
              YKJ = X(NJ2) - X(NJ5)
              ZKJ = X(NJ3) - X(NJ6)

              if(periodicBC)then
                 xkj = pbc_mic( xkj )
                 ykj = pbc_mic( ykj )
                 zkj = pbc_mic( zkj )
              endif

              RIJ0 = DHO2
c     dx       RIJ0 = XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
              RKJ0 = XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ
              RIK0 = SQRT(RIJ0*RKJ0)
              dp = XIJ*XKJ+YIJ*YKJ+ZIJ*ZKJ
              CA0 = dp/RIK0

              EPSHB = epshb_mcmc(i)

              lcont = .true.
              do iw2 = 1, ivnu  ! while (iw2<=ivnu)
                abei = epou(iw2)
                lcont = lcont .and. (.not. 
     $               ( (abei>0.01d0 .or. abei<-0.01d0) .and.
     $               ( ni(i)==ihco(iw2) .or. nj(i)==jhco(iw2) ) ))
              enddo

c     ----   hb energy function only considering cos(NHO), JULY 2002      
c              if (ca0 < 0.0d0 .and. lcont) then !!1ST DISCONTINUITY 
              if (ca0 < 0.0d0) then !! NEW 2009

                 
                 call hbmcmc(lambda,dho2,ehhb,
     $                epshb,shb2,ca0,dp,rij0,rkj0,
     $                rcut2_hb_mcmc_in,rcut2_hb_mcmc_out,
     $                fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3,
     $                xij,yij,zij,xkj,ykj,zkj)
                 

                 EHHB1 = EHHB1 + EHHB

c     ---       force for Oxygen    ! -[dcosa2/dx]*mu*switch -[cosa2]*[d(mu*switch)/dr]*x/r
                 F(NI1)=F(NI1) -fx1     +fx2
                 F(NI2)=F(NI2) -fy1     +fy2 
                 F(NI3)=F(NI3) -fz1     +fz2 
c     ---       force for Hydrogen
                 F(NJ4)=F(NJ4) +fx1+fx3 -fx2 
                 F(NJ5)=F(NJ5) +fy1+fy3 -fy2 
                 F(NJ6)=F(NJ6) +fz1+fz3 -fz2 
c     ---       force for Nitrogen
                 F(NJ1)=F(NJ1)     -fx3
                 F(NJ2)=F(NJ2)     -fy3
                 F(NJ3)=F(NJ3)     -fz3

c               if (dho2 <= rcut2_4b_out) then !!2ND DISCONTINUITY
              if (dho2 <= rcut2_4b_out .and. lcont)then !!NEW2009
                    dho = sqrt(dho2)
c     ---         precalculate the gaussians for cooperativity
                    ivnu = ivnu+1
                    ihco(ivnu) = ni(i)
                    jhco(ivnu) = nj(i)
                    rnc(ivnu) = rncoe(i)
                    ini(ivnu) = ivi(i)
                    inj(ivnu) = ivj(i)
                    epou(ivnu) = ehhb
                    DFNOT = DHOREF-DHO
                    EH11(ivnu) = dexp(-0.5d0*DFNOT*DFNOT) 

                    dho2s(ivnu) = dho2
                    dfo(ivnu) = -dfnot/dho
                    hx(ivnu) = xa
                    hy(ivnu) = ya
                    hz(ivnu) = za
                 endif          !! DHO < 8.0d0
              end if 
           end if
        ENDIF                   !! RNCOE for hb

 100  CONTINUE                  !! end loop on two-body interactions

      do i=1,ivnu-1
        i1b = ini(i)
        i2b = inj(i)
        i1 = i1b+1
        i2 = i1b+2
        i3 = i1b+3
        i4 = i1b+4
        j1 = i2b+1
        j2 = i2b+2
        iabi = abs(i1b-i2b)
        icb1 = ichain(ihco(i))
        icb3 = ichain(jhco(i)+2)
        logic1i = iabi==5
        logic21i = iabi==4
        logic22 = (i1>=1) .and. (i4<=nres)
        logic42 = (i1b>=1 .and. i1b<=nres)
        logic43 = (i2b>=1 .and. i2b<=nres)

        logic51i = (i1b<=i2b)
        logic58i = (i1b<i2b)

        do k=i+1,ivnu
          EHHB2 = zero
c     dx       EUH = zero
          
c     dx       i1b = ini(i)
c     dx       i2b = inj(i)
          i3b = ini(k)
          i4b = inj(k)
c     dx       i1 = i1b+1
c     dx       i2 = i1b+2
c     dx       i3 = i1b+3
c     dx       i4 = i1b+4
c     dx       j1 = i2b+1
c     dx       j2 = i2b+2
c     dx       iabi = abs(i1b-i2b)
          iabk = abs(i3b-i4b)
          ki31 = i3b - i1b
          ki42 = i4b - i2b
c     dx       kiabi = abs(i1b-i3b)
c     dx       kiabj = abs(i2b-i4b)
          kiabi = abs(ki31)
          kiabj = abs(ki42)
          
c     dx       icb1 = ichain(ihco(i))
          icb2 = ichain(ihco(k))
c     dx       icb3 = ichain(jhco(i)+2)
          icb4 = ichain(jhco(k))

c--   SET PI-helix to zero
c     dx       logic1 = (iabi==5) .and. (iabk==5)
          logic1 = logic1i .and. (iabk==5)
          logic2 = (ki31==1) .and. (ki42==1) !(i3b == i1) .and. (i4b == j1)
          logic3 = (ki31==2) .and. (ki42==2) !(i3b == i2) .and. (i4b == j2)
          logic6039 = logic1 .and. (logic2 .or. logic3)
c     dx       iabi == 5 && iabk == 5 && (i3b == 1 + i1b && i4b == 1 + i2b).or.(i3b == 2 + i1b && i4b == 2 + i2b)
          
c     dx       logic21 = (iabi==4) .and. (iabk==4)
          logic21 = logic21i .and. (iabk==4)
c     dx       logic22 = (i1>=1) .and. (i4<=nres)
          logic23 = icb2==icb1
          logic4 = logic23 .and. logic2 .and. logic21 .and. logic22
c     dxt      logic4 = logic21 .and. logic22 .and. logic23

          logic24 = logic23     !! icb2 == icb1
          logic25 = icb4==icb3
          logic26 = logic24 .and. logic25
          ILIM = 1
          if (logic26) ILIM = 5

          logic31 = (i3b==i2b) .and. (i4b==i1b)
          logic32 = iabi >= ILIM
          logic33 = iabk >= ILIM
          logic34 = (kiabi==2) .and. (kiabj==2)
          logic5 = (logic31 .and. logic32) .or. 
     $         (logic34 .and. logic32 .and. logic33)
c     dxt      logic5 = logic32 .and. (logic31 .or. (logic33 .and. logic34))

          logic41 = logic31     !! (i1b == i4b) .and. (i2b == i3b)
c     dx       logic42 = (i1b>=1 .and. i1b<=nres)
c     dx       logic43 = (i2b>=1 .and. i2b<=nres)
          logic6 = logic41 .and. logic42 .and. logic43

          logic51 = logic51i .and. (i3b<=i4b)
          logic52 = (i1b==(i3b-2)) .and. (i2b==(i4b-2))
          logic7 = logic51 .and. logic52

          logic8 = logic42 .and. logic43
          logic53 = logic23     !! icb2 == icb1
          logic54 = logic25     !! icb4 == icb3
          logic55 = logic53 .and. logic54 .and. (iabk<30)
c     dxt      logic55 = (icb1/=icb2).or.(icb3/=icb4).or.(iabk>=30)
          logic56 = (i3b==i2b) .and. (i4b==i1b+2)
          logic57 = (.not. logic55) .and. logic56
c     dx       logic58 = (i1b<i2b) .and. (i3b>i4b)
          logic58 = logic58i .and. (i3b>i4b)
c     dxt      logic58 = (i1b>=i2b) .or. (i3b<=i4b)
          logic59 = logic57 .and. (.not. logic58)
          logic1333 = logic59 .and. logic8
          logic9 = logic1333
c     dxt      logic9 = logic55 .and. logic56 .and. logic58 .and. logic8
          


          if (.not. logic6039) then
c     ----------------------------------------------  
c-----COOPERATIVE CONTRIBUTION for a-HELIX
             if (logic4) then
c     --          for a-helix, inclusion propensities
                parah1 = walpha_foal(ialpha(i1))
                parah2 = walpha_foal(ialpha(i2))
                parah3 = walpha_foal(ialpha(i3))
                parah4 = walpha_foal(ialpha(i4))
                parah = WCOA_EHBCOA + parah1 + parah2 + parah3 + parah4
                EHHB2 = parah*EH11(i)*EH11(k) !! EUH
c-----COOPERATIVE CONTRIBUTION for antiparallel intra and inter beta-sheet
c     ---- tested on tetramer on 28 JANVO5 by Ph.D
             else if (logic5) then
                parab = EHBCOAB
                if (logic6) then
                   parab = parab + wbeta_fobe(ibeta(i1b)) 
                   parab = parab + wbeta_fobe(ibeta(i2b))
                   EHHB2 = parab*EH11(i)*EH11(k)
                else if (.not. logic7) then
c---  beta sheet propensities are not included
c     --  new for correct treatment of constraints on the
c     sequence separation between residues in parallel sheets
                   EHHB2 = parab*EH11(i)*EH11(k)
                endif
c-----COOPERATIVE CONTRIBUTION for parallel beta-sheet both for 
c     --- intra-chain and inter-chain interactions, 28JANV05 Ph.D.
             else if (logic9) then
                parab = EHBCOBP
                parab = parab + wbeta_fobe(ibeta(i1b))
                parab = parab + wbeta_fobe(ibeta(i2b))
                EHHB2 = parab*EH11(i)*EH11(k)
             end if    
 1333        continue
 7044        continue

             ehhb2i = ehhb2
             fhhb2i = dfo(i)*ehhb2i
             ehhb2k = ehhb2
             fhhb2k = dfo(k)*ehhb2k
c     dx     the EHHB2 was calculated for cut-off dhref = 8.0
             if (dho2s(i)>=rcut2_4b_in)
     $            call switch_cutoff(dho2s(i),ehhb2i,fhhb2i,
     $            rcut2_4b_in,rcut2_4b_out) 

             if (dho2s(k)>=rcut2_4b_in)
     $            call switch_cutoff(dho2s(k),ehhb2k,fhhb2k,
     $            rcut2_4b_in,rcut2_4b_out)

             EHHB1 = EHHB1 + EHHB2i 
             EHHB3 = EHHB3 + EHHB2k

             IHI3 = 3*ihco(i)
             IHI2 = IHI3 - 1
             IHI1 = IHI3 - 2

c             JHI3 = 3*jhco(i) !!3RD DISCONTINUITY
c             JHI2 = JHI3 - 1 !!3RD
c             JHI1 = JHI3 - 2 !!3RD

              JHI1 = 3*jhco(i)+1 ! NEW 2009 
              JHI2 = 3*jhco(i)+2 ! NEW 2009 
              JHI3 = 3*jhco(i)+3 ! NEW2009 

             IHK3 = 3*ihco(k)
             IHK2 = IHK3 - 1
             IHK1 = IHK3 - 2

c             JHK3 = 3*jhco(k) !!3RD DISCONTINUITY
c             JHK2 = JHK3 - 1 !! 3RD
c             JHK1 = JHK3 - 2 !! 3RD

              JHK1 = 3*jhco(k)+1 ! NEW2009 
              JHK2 = 3*jhco(k)+2 ! NEW2009 
              JHK3 = 3*jhco(k)+3 ! NEW2009 

c     dx         DXEI = DFOHX(i)*EHHB2
c     dx         DYEI = DFOHY(i)*EHHB2
c     dx         DZEI = DFOHZ(i)*EHHB2
             DXEI =    HX(i)*fhhb2i
             DYEI =    HY(i)*fhhb2i
             DZEI =    HZ(i)*fhhb2i
c     dx         DXEK = DFOHX(k)*EHHB2
c     dx         DYEK = DFOHY(k)*EHHB2
c     dx         DZEK = DFOHZ(k)*EHHB2
             DXEK =    HX(k)*fhhb2k
             DYEK =    HY(k)*fhhb2k
             DZEK =    HZ(k)*fhhb2k

             F(IHI1) = F(IHI1) + DXEI 
             F(IHI2) = F(IHI2) + DYEI
             F(IHI3) = F(IHI3) + DZEI

             F(JHI1) = F(JHI1) - DXEI
             F(JHI2) = F(JHI2) - DYEI
             F(JHI3) = F(JHI3) - DZEI

             F(IHK1) = F(IHK1) + DXEK 
             F(IHK2) = F(IHK2) + DYEK
             F(IHK3) = F(IHK3) + DZEK

             F(JHK1) = F(JHK1) - DXEK
             F(JHK2) = F(JHK2) - DYEK
             F(JHK3) = F(JHK3) - DZEK
 6039        continue
          end if
        enddo
 6027   continue
      enddo

 1450 continue       

      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE REWA4(ires,ncalls,jn,ic,ctphib,pk)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      parameter (MAXPRE = 1500) 

      logical qmetro,qsol
      character*3 resp

      COMMON/PROPE/ resp(MAXPRE)
      COMMON/MCMIN/ ENPOT,ERE2,NVALR,QMETRO
      COMMON/REWARD/ EHHB1
      COMMON/INDHB/indh(MAXPRE)
      COMMON/VP22/VPNE(MAXPRE)
      COMMON/HYDBO/ NRESI(300),NNRES,QSOL,IMEM

      dimension ctphib(*),pk(*)
      DATA ZERO/0.0d0/

c     open(unit=7,file="beta32.wa4",status="unknown")

      bo1 = abs(ctphib(jn))
      bo2 = abs(ctphib(jn+1))
      iu = 0
      indh(ires) = 0 
      VPNE(ires) = 0.0D0

      if (pk(ic) .ge. 45.19d0 .and. pk(ic).le. 45.21d0) then
         if (ctphib(jn) .ge. 0.0d0 .and. bo1 .le. 150.0d0) then 
            VPNE(ires) = 2.0 
            ev3 = ev3 + VPNE(ires)
            indh(ires) = -1
c     if(ncalls .eq. 1 .and. QMETRO) then
 1457       format(' substract 2kc for X-Pro',a,i4, ' phi>0 ')
c     endif
         endif
      endif
ccc   end ----  substract +2 if X in X-Pro >0 


ccc   substract +6 if phi_Pro >0 
      if (pk(ic) .ge. 44.19d0 .and. pk(ic).le. 44.21d0) then

         if (ctphib(jn) .ge. 0.0d0) then 
            VPNE(ires) = 6.0d0 
            ev3 = ev3 + VPNE(ires) 
            indh(ires) = -1
            GO TO 10000
         endif
c     if(ncalls .eq. 1 .and. QMETRO) then

 3457    format(' substract 6kc for Pro',a,i4, ' phi>0   ')
c     endif
      endif
ccc   end ----substract +6 if phi_Pro >0 

      
      
ccc   substract +1 if phi >0 for certain residues 
      if (pk(ic) .ge. 41.19d0 .and. pk(ic).le. 41.21d0 .or.
     1     pk(ic) .ge. 41.59d0 .and. pk(ic) .le. 41.65d0 ) then 
         if (ctphib(jn) .ge. 0.0d0 .and. resp(ires) .ne. 'XXX' .and.  
     1        ( bo1 .le. 110.0d0 .and. bo2 .le. 110.0d0) ) then 
            VPNE(ires) = 1.0d0 
            ev3 = ev3 + VPNE(ires) 
            indh(ires) = -1
c     if(ncalls .eq. 1 .and. QMETRO) then
 1467       format(' sub v3 for',a,i4, ' if phi>0  -1  ',f7.2,2x,f12.6)
c     endif
c     GO TO 10000
         endif
      endif
ccc   end ---- substract +1 if phi >0 for certain residues 

      
ccc   substract +0.5 if phi >0 for other residues 
      if (pk(ic) .ge. 42.19d0 .and. pk(ic).le. 42.21d0 ) then 
         if (ctphib(jn) .ge. 0.0d0 .and. resp(ires+1).ne.'GLY' 
     1        .and. bo1 .le. 110.0d0 ) then 
            VPNE(ires) = 0.5d0 
            ev3 = ev3 + VPNE(ires)
            indh(ires) = -1
c     if(ncalls .eq. 1 .and. QMETRO) then
c     write(7,5467) resp(ires), ires,VP0,ctphib(jn)
 5467       format(' substract 0.5kc for',a,i4, ' phi>0 ',f7.2,2x,f12.6)
c     endif
            GO TO 10000
         endif 
      endif 
      




ccc   reward epsilon for certain residues in beta 
c     if(ctphib(jn) .le. 0. .and. ( (bo1 .ge. 60. .and.  
c     1 ctphib(jn+1).ge. 30.) .or. (bo1 .ge. 60. .and. 
c     2 ctphib(jn+1).lt. 0. .and. bo2 .ge. 140.))) then  

      if( (ctphib(jn) .le. 0.0d0 .and. bo1 .ge. 60.0d0 .and.
     1     (ctphib(jn+1).ge. 40.0d0 .or. bo2 .ge. 140.0d0))
     2     .or. (ctphib(jn) .gt. 0.0d0 .and.
     3     bo1 .ge. 130.0d0 .and. bo2 .ge. 130.0d0)  ) then



         VPNE(ires) = 0.0d0
         if (pk(ic) .ge. 45.19d0 .and. pk(ic).le. 45.21d0) then
            VPNE(ires) = 0.0d0 
            go to 5200
         endif 
         if( resp(ires) .eq. 'TRP' .or. resp(ires) .eq.'SER') then
            VPNE(ires) = -0.15d0 
            go to 5200
         endif
         if( resp(ires) .eq. 'GLU') then 
            VPNE(ires) = -0.15d0
            go to 5200
         endif
         if( resp(ires) .eq. 'THR' .or. resp(ires) .eq.'ILE') then
            VPNE(ires) = -0.5d0 
            go to 5200
         endif
         if( resp(ires) .eq. 'VAL' .or. resp(ires) .eq.'TYR') then
            VPNE(ires) = -0.5d0
            go to 5200
         endif
         if( resp(ires) .eq. 'PHE') then
            VPNE(ires) = -0.5d0 
            go to 5200
         endif
         if( resp(ires) .eq. 'PRO') then
            VPNE(ires) = 2.5d0 
            go to 5200
         endif
         if( resp(ires) .eq. 'XXX') then
            VPNE(ires) = 0.0d0 
            go to 5200
         endif
         if(resp(ires).eq.'GLY') then 
            VPNE(ires) = 0.75d0  
         endif

 5200    continue
         indh(ires) = 1
         ev1 = ev1 + VPNE(ires)
c     if(ncalls .eq. 1 .and. QMETRO ) then
c     write(7,1451)ires,resp(ires),ctphib(jn),
c     1 ctphib(jn+1),VPNE(ires)
 1451    format('RES in BET',i4, 'RESP',a3,'phi ',f7.2,
     1        ' psi ',f7.2,' VPNE ',f7.2)
c     endif
         GO TO 10000
      endif
ccc   end ----- reward epsilon for certain residues in beta 



c     substract epsilon for Pro within an helix

      if(resp(ires) .eq. 'PRO') then
         if (ctphib(jn) .le. 0.0d0 .and. ctphib(jn+1).le. 0.0d0 .and.
     1        bo1 .le. 90.0d0 .and. bo2 .le. 70.0d0)then
            VPNE(ires) = 2.8d0 
            ev2 = ev2 + VPNE(ires)
            indh(ires) = 2 
c     if(ncalls .eq. 1 .and. QMETRO) then
c     write(6,1452) resp(ires), ires,ctphib(jn),ctphib(jn+1)
 1452     format(' substract v3 for ',a,i4, ' in alpha 2.0  ',f9.3,f9.3)
c     endif
         endif
         GO TO 10000
      endif 

c     end ---- substract epsilon for Pro if within an helix



c     substract epsilon for Gly within an helix

      if(resp(ires) .eq. 'GLY') then
         if (ctphib(jn) .le. 0.0d0 .and. ctphib(jn+1) .le. 0.0d0 .and. 
     1        bo1 .le. 90.0d0 .and. bo2 .le. 80.0d0 .and. 
     2        bo1 .ge. 30.0d0 .and. bo2 .ge. 10.0d0) then
            VPNE(ires) = 1.40d0  
            ev2 = ev2 + VPNE(ires)
            indh(ires) = 2 
c     if(ncalls .eq. 1 .and. QMETRO) then
c     write(6,1452) resp(ires), ires,ctphib(jn),ctphib(jn+1)
c     endif
         endif
         GO TO 10000
      endif 

c     end --- substract epsilon for Gly within an helix


c     substract epsilon for other residues within an helix

      IF (ctphib(jn) .le. 0.0d0 .and. ctphib(jn+1).le. 0.0d0 .and.
     1     bo1 .le. 80.0d0 .and. bo2 .le. 80.0d0 .and.
     2     bo1 .ge. 20.0d0 .and. bo2 .ge. 20.0d0) then 

         if (resp(ires-1) .eq. 'GLY' .and. 
     1        pk(ic) .ge. 42.19d0 .and. pk(ic).le. 42.21d0 ) then 
            VPNE(ires) = 0.0d0 
            go to 1205
         endif

         if (resp(ires) .eq. 'ALA') then
            VPNE(ires) = -0.4d0
            go to 1205
         endif 
         if(resp(ires) .eq. 'ARG' .or. resp(ires) .eq. 'LEU') then
            VPNE(ires) = -0.2d0  
            go to 1205
         endif 
         if(resp(ires) .eq. 'GLU')  then
            VPNE(ires) = -0.2d0 
            go to 1205
         endif
         if(resp(ires) .eq. 'LYS')  then
            VPNE(ires) = -0.2d0 
            go to 1205
         endif
         if (resp(ires) .eq. 'MET') then
            VPNE(ires) = -0.1d0 
            go to 1205
         endif 
         if (resp(ires) .eq. 'GLN') then  
            VPNE(ires) = -0.1d0 
            go to 1205
         endif 
         if (resp(ires) .eq. 'XXX') then
            VPNE(ires) = 0.0d0 
            go to 1205
         endif 
         if (resp(ires) .eq. 'CYS') then
            VPNE(ires) = 0.3d0 
            go to 1205
         endif 
         if(resp(ires) .eq. 'ILE' .or. resp(ires) .eq. 'SER') then
            VPNE(ires) = 0.1d0 
            go to 1205
         endif
         if(resp(ires) .eq. 'TRP') then
            VPNE(ires) = 0.5d0 
            go to 1205
         endif
         if(resp(ires) .eq. 'TYR' .or. resp(ires) .eq. 'PHE') then
            VPNE(ires) = 0.5d0 
            go to 1205
         endif
         if(resp(ires) .eq. 'VAL' .or. resp(ires) .eq. 'HID') then
            VPNE(ires) = 0.5d0  
            go to 1205
         endif
         if(resp(ires) .eq. 'ASN' .or. resp(ires) .eq. 'THR') then
            VPNE(ires) = 0.5d0   
            go to 1205
         endif
         if(resp(ires) .eq. 'ASP') then
            VPNE(ires) = 0.5d0 
            go to 1205
         endif
         if(resp(ires) .eq. 'GLY' .or. resp(ires) .eq. 'PRO') then
c     VP0 = 0.0   !! already considered 
            VPNE(ires) = 0.0d0 
            go to 1205
         endif
         VPNE(ires) = 0.0d0 
 1205    continue
         ev2 = ev2 + VPNE(ires)
         indh(ires) = 2
c     if(ncalls .eq. 1 .and. QMETRO) then
c     write(7,2451)ires,resp(ires),ctphib(jn),ctphib(jn+1),VPNE(ires)
 2451    format(' for residue in alpha',i4, '   ',a,' phi ',f7.2,
     1        ' psi ',2f7.2)
c     endif
      ENDIF

c     end ---  substract epsilon for other residues within an helix


cccc  TO BE FURTHER EXAMINED

      if (pk(ic) .ge. 45.19d0 .and. pk(ic).le. 45.21d0 ) then 
         IF (ctphib(jn) .le. 0.0d0 .and. ctphib(jn+1).le. 0.0d0 .and.
     1        bo1 .le. 110.0d0 .and. bo2 .le. 110.0d0) then
            VPNE(ires) = 1.0    !! GUESS 
            ev2 = ev2 + VPNE(IRES)
            indh(ires) = 2
c     if(ncalls .eq. 1 .and. QMETRO) then
c     write(6,4451)ires,resp(ires),ctphib(jn),ctphib(jn+1),VPNE(ires)
c     4451  format(' for X-Pro in alpha',i4, '   ',a,' phi ',f7.2,
c     1 ' psi ',2f7.2)
c     endif
         ENDIF
      endif

10000 CONTINUE 
c     write(7,*) ev1,ev2,ev3
c     stop

c     ev3 = ZERO

      RETURN
      END

      subroutine switch_cutoff(r2,ene_switched,for_switched,ri2,ro2)

        implicit none

        real*8 r2,ene_switched,for_switched
        real*8 ri2,ro2
cdx     real*8 r1
        real*8 rd6
        real*8 sw_func,d_sw_func
        
cdx     r1 = sqrt(r2)
        rd6 = ro2-ri2
        rd6 = 1.0d0/(rd6*rd6*rd6)
        
        sw_func = (ro2-r2)*(ro2-r2)*(ro2+2.0d0*r2-3.0d0*ri2)*rd6
        d_sw_func = 12.0d0*(ro2-r2)*(ri2-r2)*rd6 !*r1

cdx     for_switched = for_switched*r1
        for_switched = for_switched*sw_func - ene_switched*d_sw_func !/r1
cdx     for_switched = for_switched/r1

        ene_switched = ene_switched*sw_func
        
      return
      end

      subroutine ljcasc(lambda,da2,eahyd,df,rncoe,ct0,ct2,
     $          r2in,r2out,vamax)
        implicit none
        
        real*8 da2,eahyd,df,rncoe,ct0,ct2,r2in,r2out
        real*8 cd8,cd2,vamax,da6
        real*8 eight,G,Rmin,V,r,DV1,DV2,G6,KG
        parameter(KG=0.7,eight=8.0d0)
        real*8 lambda

C---- 12-6 potential if attractive, otherwise 6-repulsive potential
          CD2 = ct2/da2
          CD8 = CD2*CD2*CD2*CD2
          Rmin = vamax

          if (rncoe .gt. 0.0d0) then
            Rmin = vamax
            r = sqrt(da2)
            da6=da2*da2*da2
            G = -KG*exp(2.*(Rmin-0.5)/5.0)*(Rmin-0.5)
            G6 = G**6
       V=G6*exp(-2*r)/da6+lambda*0.6563701*(TanH(2*(r-Rmin-0.5))-1)
c       modif lambda Oct 2011
            Eahyd = ct0*V
            DV1 = -2*exp(-2*r)*G6*(3+r)/(da6*r)
            DV2 = lambda/CosH(1+2*Rmin-2*r)
c       modif lambda Oct 2011
            DF = ct0*(DV1+1.312740*DV2*DV2)/r
           else if (rncoe .lt. 0.0d0) then
               Eahyd = -ct0*CD8   !! orig sign -1
               DF =  eight*ct0*CD8/DA2   !! orig sign +1
          endif

c ------- store energy in ehydro and forces
          if (da2>=r2in)
     $    call switch_cutoff(DA2,eahyd,DF,r2in,r2out)

      return
      end


      subroutine hbmcmc(lambda,dho2,ehhb,
     $      epshb,shb2,ca0,dp,rij0,rkj0,
     $      r2in,r2out,
     $      fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3,
     $      xij,yij,zij,xkj,ykj,zkj)

        implicit none
        
        real*8 dho2,ehhb,dfact1,dfact2,dfact3
        real*8 shb2,ca0,dp,rij0,rkj0,r2in,r2out
        
        real*8 dint2, dint10, urij
        real*8 ehha, epshb
        real*8 dmu_dr

        real*8 asq
        real*8 dpehha, rijk0, rijk1

        real*8 fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3
        real*8 xij,yij,zij,xkj,ykj,zkj

        real*8 dfact1x1,dfact1y1,dfact1z1,dfact1x2,dfact1y2,dfact1z2
        real*8 dfact2x,dfact2y,dfact2z,dfact3x,dfact3y,dfact3z

        real*8 lambda


        dint2 = shb2/dho2
        dint10 = dint2*dint2*dint2*dint2*dint2
        urij = dint10 * (5.0d0*dint2 - lambda*6.0d0)

        ehha = EPSHB * URIJ
        dmu_dr = epshb * 60.0d0 * dint10 * (dint2 - lambda*1.0d0) / dho2

        if (dho2>=r2in) call switch_cutoff(dho2,ehha,dmu_dr,r2in,r2out)

        ASQ  = CA0 * CA0
        EHHB = ehha * asq           ![cosa2]*[mu*switch]
        dmu_dr = dmu_dr * asq       ![cosa2]*[d(mu*switch)/dr]/r

        rijk1 = 2.0d0*dp*ehha/(rij0*rkj0)
cdx     dfact1 = 2.0d0*dp   *ehha/(rij0*rkj0)
cdx     dfact2 = 2.0d0*dp*dp*ehha/(rij0*rij0*rkj0)
cdx     dfact3 = 2.0d0*dp*dp*ehha/(     rij0*rkj0*rkj0)
        dfact1 = rijk1
        dfact2 = rijk1*dp/rij0
        dfact3 = rijk1*dp/rkj0

        dfact1x1 = dfact1*xkj; dfact1x2 = dfact1*xij
        dfact1y1 = dfact1*ykj; dfact1y2 = dfact1*yij
        dfact1z1 = dfact1*zkj; dfact1z2 = dfact1*zij
                               dfact2x  = dfact2*xij
                               dfact2y  = dfact2*yij
                               dfact2z  = dfact2*zij
        dfact3x  = dfact3*xkj
        dfact3y  = dfact3*ykj
        dfact3z  = dfact3*zkj

        fx1 = (dfact1x1-dfact2x)
        fy1 = (dfact1y1-dfact2y)
        fz1 = (dfact1z1-dfact2z)
        fx2 = dmu_dr*xij
        fy2 = dmu_dr*yij
        fz2 = dmu_dr*zij
        fx3 = (dfact1x2-dfact3x)
        fy3 = (dfact1y2-dfact3y)
        fz3 = (dfact1z2-dfact3z)

        return
      end

C     Application the nearest-image convention.
C     As usual in Fortran, the atoms are place in a box going from 
C     [-0.5*box_length, 0.5*box_length]

C==============================================================================
C=======================    periodic boundary condition  ======================
C=======================        for each atom            ======================
C==============================================================================


      function pbc_mic( x )
      
      implicit none

      integer MAXPRE,MAXNAT,MAXPAI,MAXXC,MAXTTY
      parameter (MAXPRE = 1500)    !! maximum number of residus 
      parameter (MAXNAT = MAXPRE*6)  !! maximum number of atoms
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord 
      
      real*8 pbc_mic, x
      
      real*8 box_length, inv_box_length

      common/pbcBL/box_length, inv_box_length
      
       pbc_mic = x - box_length * dnint( x * inv_box_length)

      
      end





      
cdxIn[1]:= f[r]=(ro^2-r^2)^2*(ro^2+2*r^2-3*ri^2)/(ro^2-ri^2)^3
cdx
cdx           2     2 2     2       2     2
cdx        (-r  + ro )  (2 r  - 3 ri  + ro )
cdxOut[1]= ---------------------------------
cdx                      2     2 3
cdx                  (-ri  + ro )
cdx
cdxIn[2]:= D[f[r],r]
cdx
cdx               2     2 2          2     2      2       2     2
cdx        4 r (-r  + ro )    4 r (-r  + ro ) (2 r  - 3 ri  + ro )
cdxOut[2]= ---------------- - ------------------------------------
cdx             2     2 3                    2     2 3
cdx         (-ri  + ro )                 (-ri  + ro )
cdx
cdxIn[3]:= Simplify[%]
cdx
cdx                2     2    2     2
cdx        -12 r (r  - ri ) (r  - ro )
cdxOut[3]= ---------------------------
cdx                  2     2 3
cdx               (ri  - ro )
cdx

