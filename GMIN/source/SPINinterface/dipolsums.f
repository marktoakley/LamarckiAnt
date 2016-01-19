      subroutine dipsum(a0x0,a0y0,adx,ady,adz,asum)
c
c   CALCULATION OF THE DIPOLE INTERACTION OF A PARTICULAR SPIN (1) 
c   IN THE UNIT CELL WITH ANOTHER SPIN (2) IN THE SAME UNIT CELL, 
c   PLUS ALL SPINS OF THE UNIT CELL COPIES FORMING AN INFINITE 2D
c   RECTANGULAR ARRAY OF UNIT CELLS. THE INTERACTION OF SPIN (1)
c   WITH ITS PICTURES IN THE NEIGHBORING UNIT CELLS IS ALSO INCLUDED 
c
c   2D ARRAY  periodic array of unit cell with infinite lateral 
c                 extension, spans the XY- plane 
c             finite unit cell thickness along Z-axis, which is the 
c                 array normal, hence no periodic extension along Z
c             spins (1) and (2) need not to be in same plane
c  
c   A0X0    unit cell extension in X- direction 
c   A0Y0    unit cell extension in Y- direction 
c           (A0X0 = A0Y0:  square unit cell) 
c
c           Absolute distances of spins (1) and (2) in unit cell: 
c   ADX     separation in X- direction 
c   ADY     separation in Y- direction 
c   ADZ     separation in Z- direction (thickness of unit cell) 
c           (the subroutine distinguishes between the cases 
c            ADZ = 0  and  |ADZ| > 0) 
c
c   ADEX0   relative distances in X- direction in unit of A0X0
c   ADEY0   relative distances in Y- direction in unit of A0Y0
c           (NECESSARY: |ADEX0| < 1 and |ADEY0| < 1) 
c
c   The following lattice sums are determined: 
c
c     S00 = ADZ * ADZ * sum_(i,j) 1 / R_ij**5 
c     Sxx = sum_(i,j) X_i * X_i / R_ij**5 
c     Syy = sum_(i,j) Y_j * Y_j / R_ij**5 
c     Sxy = sum_(i,j) X_i * Y_j / R_ij**5 
c     Sx0 = ADZ * sum_(i,j) X_i / R_ij**5 
c     S0y = ADZ * sum_(i,j) Y_j / R_ij**5 
c
c     with:    X_i = i * A0X0 + ADX
c              Y_j = j * A0Y0 + ADY
c
c     absolute distance   R_ij = SQRT (X_i**2 + Y_j**2 + ADZ**2) 
c
c     the indices i and j are running from -infinity to +infinity, 
c     the terms with R_ij = 0 are omitted from the sums if ADZ = 0
c
c     vector ASUM(n), six components: 
c              n = 1: z*z*S00,  n = 2: Sxx,  n = 3: Syy,  n = 4: Sxy, 
c              n = 5: z*Sx0,    n = 6: z*S0y 
c
c
      implicit double precision (a,b)
      dimension asum(6),asurr(6),asuxx(6),asurxx(6)

      api=datan(1.d0)*4.        ! pi
      nja=10
      njb=1000
      njc=100000

      do ix=1,6
         asum(ix)=0.d0
      enddo                     ! ix

      adex0=adx/a0x0
      adey0=ady/a0y0
      az=dabs(adz)
      if(dabs(adex0).ge.1.d0.or.dabs(adey0).ge.1.d0) goto 901 

      adex00=adex0
      adey00=adey0
      if(dabs(adex0).gt.0.5d0) adex00=adex0-dsign(1.d0,adex0)
      if(dabs(adey0).gt.0.5d0) adey00=adey0-dsign(1.d0,adey0)

      a0x=a0x0
      a0y=a0y0
      adex=adex00
      adey=adey00
      ipro=0
      if(dabs(adex00).lt.dabs(adey00)) then 
         a0x=a0y0
         a0y=a0x0
         adex=adey00
         adey=adex00
         ipro=1
      endif

      if(az.lt.1.d-4) then 
! ================= adz=0 =============================
         asur03=0.d0
         do iks=1,njb
            do ik=-nja*iks,nja*iks
               if(iks.gt.1.and.ik.ge.-nja*(iks-1).and.ik.le.
     #              nja*(iks-1)) goto 102
               agk=a0x*(ik-adex)
               if(dabs(agk).lt.1.d-12) goto 102
               asux=a0y/api/dabs(agk)
               asurr(2)=asux*asux/2.
               asurr(3)=0.25/agk/agk
               asurr(4)=0.d0
               do il=1,njc
                  avv=api*il/a0y
                  auu=2.*dabs(agk)*avv
                  asc=dcos(2.*api*adey*il)
                  ass=dsin(2.*api*adey*il)
                  abk0=besk0(auu,abk1,abk2)
                  asurr(2)=asurr(2)+2.*il*il*asc*abk2
                  asurr(3)=asurr(3)+asc*(abk1/dabs(agk)
     #                 -2.*avv*abk0)*avv
                  asurr(4)=asurr(4)-2.*il*il*ass*abk1
                  if(il.gt.5) then
                     asum0=dabs(asurr(4)-asuxx(4))+
     #               dabs(asurr(2)-asuxx(2))+dabs(asurr(3)-asuxx(3))
                     if(asum0.lt.1.d-9) goto 101                     
                  endif
                  do ix=2,4
                     asuxx(ix)=asurr(ix)
                  enddo         ! ix
               enddo            ! il
               write(6,*) '# accuracy problem in L loop !'
 101           asum(4)=asum(4)+asurr(4)*dsign(1.d0,agk)
               asum(2)=asum(2)+asurr(2)
               asum(3)=asum(3)+asurr(3)
               if(dabs(adex).lt.1.d-7) then 
                  arr=a0y*dabs(ik-adey)
                  if(arr.gt.1.d-7) asur03=asur03+1.d0/arr/arr/arr
               endif
 102           continue 
            enddo               ! ik
            aba=(nja*iks+1)*1.d0
            att=aba*a0y*a0y/api/api/a0x/a0x/(aba*aba-adex*adex)
            if(iks.gt.1) then
               asum0=dabs(asum(2)+att-asurxx(2)-attx)
     #              +dabs(asur03-asur0x)-dabs(asum(2)-asurxx(2))
               do ix=2,4
                  asum0=asum0+dabs(asum(ix)-asurxx(ix))
               enddo            ! ix
               if(asum0.lt.1.d-7) goto 103
            endif
            do ix=2,4
               asurxx(ix)=asum(ix)
            enddo               ! ix
            asur0x=asur03
            attx=att
         enddo                  ! iks
         write(6,*) '# accuracy problem in K loop'
 103     asum(2)=(asum(2)+att)*8.*api*api/a0y/a0y/a0y/3.   ! Sxx
         asum(3)=asum(3)*8./a0y/3.+asur03                  ! Syy
         asum(4)=asum(4)*8.*api*api/a0y/a0y/a0y/3.         ! Sxy
      else

! ================= |adz|>0 ===========================
         njd=10
         if(az.lt.0.5d0) njd=min0(5000,nint(5./az))
         asum(1)=0.25d0/az
         asum(2)=0.25d0/az
         asum(3)=0.25d0/az
         do ik=1,njd
            auxk=api*ik/a0x
            acxk=dcos(2.*api*ik*adex)
            asxk=dsin(2.*api*ik*adex)
            aex=dexp(-2.*az*auxk)
            auyk=api*ik/a0y
            acyk=dcos(2.*api*ik*adey)
            asyk=dsin(2.*api*ik*adey)
            aey=dexp(-2.*az*auyk)
            do ii=1,6
               asurr(ii)=0.d0
            enddo               ! ii
            do il=1,njd
               auyl=api*il/a0y
               acxl=dcos(2.*api*il*adex)
               aukl=dsqrt(auxk*auxk+auyl*auyl)
               aekl=dexp(-2.*az*aukl)
               auxl=api*il/a0x
               acyl=dcos(2.*api*il*adey)
               asyl=dsin(2.*api*il*adey)
               aulk=dsqrt(auyk*auyk+auxl*auxl)
               aelk=dexp(-2.*az*aulk)                                    
               asurr(1)=asurr(1)+acyl*aekl*(aukl+0.5/az)           ! S00
               asurr(2)=asurr(2)+acyl*aekl*(0.5/az-auxk*auxk/aukl) ! Sxx
               asurr(3)=asurr(3)+acxl*aelk*(0.5/az-auyk*auyk/aulk) ! Syy
               asurr(4)=asurr(4)+auyl*asyl*aekl/aukl               ! Sxy
               asurr(5)=asurr(5)+acyl*aekl                         ! Sx0
               asurr(6)=asurr(6)+acxl*aelk                         ! S0y
            enddo               ! il
            asum(1)=asum(1)+acyk*aey*(auyk+0.5/az)
     #           +acxk*(2.*asurr(1)+(0.5/az+auxk)*aex)
            asum(2)=asum(2)+acyk*aey/az/2.
     #           +acxk*(2.*asurr(2)+(0.5/az-auxk)*aex)
            asum(3)=asum(3)+acxk*aex/az/2.
     #           +acyk*(2.*asurr(3)+(0.5/az-auyk)*aey)
            asum(4)=asum(4)+2.*auxk*asxk*asurr(4)
            asum(5)=asum(5)-asxk*auxk*(aex+2.*asurr(5))
            asum(6)=asum(6)-asyk*auyk*(aey+2.*asurr(6))
         enddo                  ! ik
         axy=8.*api/a0x/a0y/3.
         do ix=1,6
            asum(ix)=asum(ix)*axy
         enddo                  ! ix
      endif

! ============================================================
      if(ipro.eq.1) then 
         asurx=asum(3)
         asum(3)=asum(2)
         asum(2)=asurx
         asurx=asum(6)
         asum(6)=asum(5)
         asum(5)=asurx
      endif
      return
 901  write(6,*) '# ADEX0 or ADEY0 equal to or larger than 1'
      stop
      end
!
      double precision function besk0(x,besk1,besk2) 
!  =================================================
!   calculates modified Bessel functions K_\nu(x) 
!   (Abramovitz & Stegun, 9.8.5 - 9.8.8, 9.6.26)
!  =================================================
      implicit double precision (b,x,y) 
      if(x.le.0.d0) goto 600 
      if(dabs(x)-2.d0.le.0.d0) then
         y=x*x/4.
         yy=x*x/14.0625
         bi1=x*(0.5d0+yy*(0.87890594+yy*(0.51498869+yy*(0.15084934+
     #        yy*(0.02658733+yy*(0.00301532+yy*0.00032411))))))
         besk1=dlog(x/2.)*bi1+(1.d0+y*(0.15443144+            ! K_1
     #        y*(-0.67278579+y*(-0.18156897+y*(-0.01919402+
     #        y*(-0.00110404-y*0.00004686))))))/x
         bi0=1.d0+yy*(3.5156229+yy*(3.0899424+yy*(1.2067492+
     #        yy*(0.2659732+yy*(0.0360768+yy*0.0045813)))))
         besk0=-dlog(x/2.)*bi0+(-0.57721566+y*(0.4227842+     ! K_0
     #        y*(0.23069756+y*(0.0348859+y*(0.00262698+
     #        y*(0.0001075+y*0.0000074))))))
      else
         y=2./x
         besk1=(dexp(-x)/dsqrt(x))*(1.25331414+y*(0.23498619+
     #        y*(-0.0365562+y*(0.01504268+y*(-0.00780353+
     #        y*(0.00325614-y*0.00068245))))))                ! K_1
         besk0=(dexp(-x)/dsqrt(x))*(1.25331414+y*(-0.07832358+
     #        y*(0.02189568+y*(-0.01062446+y*(0.00587872+
     #        y*(-0.0025154+y*0.00053208))))))                ! K_0
      endif
      besk2=besk0+2.*besk1/x                                  ! K_2
      return
 600  write(6,*) ' big problem with bessel: ',x
      stop
      end