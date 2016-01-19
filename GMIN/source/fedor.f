C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C

           SUBROUTINE ASAR1(R,Ex,Ex1)
c          REAL*8  FUNCTION UTN(R,I,IP)
C  TERM X0G+ AR2 (AZIZ - JCP 99 (1993) 4518)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      save
      DATA E/143.235d0/,Rm/3.757d0/, A/9.03228328d0/,B/-2.37132823d0/,
     * C6,C8,C10/1.09309955d0,.51568309d0,.32521242d0/,AA/8.73933927D4/
     *,C12,C14/.27818156d0,.31111959d0/,ro/1.107d0/
     *,NCALL/1/,Ekau/3.1668d-6/,Rau/0.52918d0/
      IF(NCALL.EQ.1) THEN
      Rm=Rm/Rau
      E=E*Ekau
      g61=2.1d0/6.d0
      g62=0.109d0/sqrt(6.d0)
      g81=2.1d0/8.d0
      g82=0.109d0/sqrt(8.d0)
      g101=2.1d-1
      g102=0.109d0/sqrt(10.d0)
      g121=2.1d0/12.d0
      g122=0.109d0/sqrt(12.d0)
      g141=0.15d0
      g142=0.109d0/sqrt(14.d0)
      NCALL=2
                     END IF
      X=R/Rm
      X2=X*X
      X4=X2*X2
      X6=X4*X2
      X8=X4*X4
      X10=X6*X4
      X12=X6*X6
      X14=X8*X6
      C6X6=C6/X6
      C8X8=C8/X8
      C10X10=C10/X10
      C12X12=C12/X12
      C14X14=C14/X14
      rr=ro*R
      F=1.d0-rr**1.68*exp(-0.78d0*rr)
      g6s=1.d0-exp(-g61*rr-g62*rr*rr)
      g6d=g6s*g6s
      g6t=g6d*g6d
      g6=g6t*g6d
      g8s=1.d0-exp(-g81*rr-g82*rr*rr)
      g8d=g8s*g8s
      g8t=g8d*g8d
      g8=g8t*g8t
      g10s=1.d0-exp(-g101*rr-g102*rr*rr)
      g10d=g10s*g10s
      g10t=g10d*g10d
      g10=g10t*g10t*g10d
      g12s=1.d0-exp(-g121*rr-g122*rr*rr)
      g12d=g12s*g12s
      g12t=g12d*g12d
      g12=g12t*g12t*g12t
      g14s=1.d0-exp(-g141*rr-g142*rr*rr)
      g14d=g14s*g14s
      g14t=g14d*g14d
      g14o=g14t*g14t
      g14=g14o*g14t*g14d
      Vscf=AA*EXP((-A+B*X)*X)
      Vdisp=-F*(C6X6*g6+C8X8*g8+C10X10*g10+C12X12*g12+C14X14*g14)

        Ex=(Vscf+Vdisp)*E

      Vscf1=Vscf*(-A+2.d0*B*X)
      F1=(F-1.d0)*(1.68d0/rr-0.78d0)
      g6p=6.d0*g6/g6s*(1.d0-g6s)*(g61+2.d0*g62*rr)
      g8p=8.d0*g8/g8s*(1.d0-g8s)*(g81+2.d0*g82*rr)
      g10p=10.d0*g10/g10s*(1.d0-g10s)*(g101+2.d0*g102*rr)
      g12p=12.d0*g12/g12s*(1.d0-g12s)*(g121+2.d0*g122*rr)
      g14p=14.d0*g14/g14s*(1.d0-g14s)*(g141+2.d0*g142*rr)
      Vdisp1=ro*Rm*(F1*Vdisp/F
     &-F*(C6X6*g6p+C8X8*g8p+C10X10*g10p+C12X12*g12p+C14X14*g14p))
     &+F*(6.d0*C6X6*g6+8.d0*C8X8*g8+10.d0*C10X10*g10
     &+12.d0*C12X12*g12+14.d0*C14X14*g14)/X

        Ex1=(Vscf1+Vdisp1)*E/Rm

      RETURN
      END



      subroutine energyar(r,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
      implicit NONE
      DOUBLE PRECISION r,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1
      DOUBLE PRECISION r0,ut1,ar2ut0
      save

      r0=r
      Esu=ar2ut0(r0,1,ut1)
      Esu1=ut1
        Esg=ar2ut0(r0,2,ut1)  
      Esg1=ut1
        Epu=ar2ut0(r0,3,ut1)  
      Epu1=ut1
        Epg=ar2ut0(r0,4,ut1)  
      Epg1=ut1
      call asar1(r,Ex,Ex1)
      return
      end

          DOUBLE PRECISION  FUNCTION ar2UT0(RR,I,UT1)
      implicit NONE
      INTEGER NB(8)
      INTEGER, DIMENSION(8) :: NB0=(/ 17,17,17,17,17,17,17,17 /)
      INTEGER :: NF=8
      INTEGER :: NCALL=1
      INTEGER nend(8)
C    * ,L(8)
      DOUBLE PRECISION COEF(4,33,8),RR,UT1,SOAT,EAS,RMIN,RFNBI,AAE,AE,PPVALU
      INTEGER I,NF1,J,J1,JN,JJ,N1,NN,IBCBEG,K,NENDI,NBI,NB0I
C    * ,BREAK(33,8),S(17,6)
CCC  AR2* (Spiegelmann,1984)
      DOUBLE PRECISION, DIMENSION(136) ::  RF= (/
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15. /)
      DOUBLE PRECISION, DIMENSION(136) :: ENF=(/
     *-.55165,-.57103,-.57143,-.57120,-.57051,-.56819,-.56444,-.56055,    
     *-.55401,-.55017,-.54784,-.54667,-.54615,-.54586,-.54582,-.54573,
     *-.54551,
     *-.51143,-.53534,-.53663,-.53730,-.53754,-.53725,-.53642,-.53591,    
     *-.53664,-.53912,-.54116,-.54277,-.54392,-.54518,-.54563,-.54572,
     *-.54551,
     *-.47586,-.49780,-.49872,-.50108,-.50994,-.51996,-.52839,-.53382,    
     *-.53961,-.54275,-.54402,-.54474,-.54518,-.54565,-.54578,-.54571,
     *-.54551,
     *-.50651,-.52661,-.52702,-.52676,-.52599,-.53231,-.53737,-.54033,    
     *-.54303,-.54456,-.54498,-.54524,-.54545,-.54572,-.54580,-.54571,
     *-.54551,
     *-.54822,-.56754,-.56792,-.56765,-.56691,-.56448,-.56055,-.55647,    
     *-.54951,-.54525,-.54261,-.54120,-.54048,-.53992,-.53971,-.53941,
     *-.53889,
     *-.50609,-.52941,-.53058,-.53113,-.53123,-.53060,-.52935,-.52849,    
     *-.52893,-.53143,-.53341,-.53499,-.53617,-.53760,-.53826,-.53866,
     *-.53889,
     *-.47526,-.49699,-.49769,-.49774,-.50045,-.51213,-.52052,-.52595,    
     *-.53176,-.53493,-.53627,-.53710,-.53768,-.53839,-.53873,-.53888,
     *-.53889,
     *-.50591,-.52605,-.52647,-.52620,-.52543,-.52757,-.53254,-.53525,    
     *-.53755,-.53869,-.53888,-.53900,-.53912,-.53931,-.53937,-.53924,
     *-.53889 /)
        save

      ar2ut0=0.0
      ut1=0.0
      IF(I.GT.NF) return
      IF(NCALL.GT.1) GO TO 5

      NCALL=2
      soat=soat/8065/27.212
      NB(1)=1
      NF1=NF-1
      DO J=1,NF1
        NB(J+1)=NB(J)+NB0(J)
      END DO
c     EAS=0.
      EAS=ENF(NB0(1))
      RMIN=0.
      DO 2 J=1,NF
      J1=NB(J)
      JN=J1+NB0(J)-1
      nend(j)=jn
      RMIN=MAX(RMIN,RF(J1))
C     EAS=ENF(JN)
      DO 1 JJ=J1,JN
      ENF(JJ)=ENF(JJ)-EAS 
    1 continue
    2 CONTINUE

      DO 4 J=1,NF
      N1=NB(J)
      NN=NB0(J)
      IBCBEG=0
      IF(J.EQ.3) THEN
       IBCBEG=2
       COEF(2,1,J)=0.
               END IF
      COEF(2,NN,J)=0.
      DO K=1,NN
        COEF(1,K,J)=ENF(N1+K-1)
      END DO
C     L(J)=2*NB0(J)-1
      CALL CUBSPL(RF(N1),COEF(1,1,J),NN,IBCBEG,1)
C     CALL TAUTSP(RF(NB(J)),ENF(NB(J)),NB0(J),GAM,S,BREAK(1,J),
C    *                                  COEF(1,1,J),L(J),K,IFLAG)
    4 CONTINUE

    5 continue
        NENDI=NEND(I)
      IF(RR.LT.RF(NENDI)) then
        nbi=nb(i)
        nb0i=nb0(i)
        rfnbi=rf(nbi)
      if(rr.lt.rfnbi) then
      AAe=enf(nbi) +0.1
      ae=-PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,rfnbi,1)/AAe
      ar2ut0=AAe*exp(-ae*(rr-rfnbi)) -0.1
      ut1=-ae*(ar2ut0 +0.1)
                  else
C     UT0=PPVALU(BREAK(1,I),COEF(1,1,I),L(I),4,RR,0)
      ar2UT0=PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,RR,0)
      UT1=PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,RR,1) 
                  end if
                      else
      ar2UT0=ENF(NENDI)
      UT1=0.0
                      end if
      RETURN
        END


      subroutine rgnii(n,x,grad,ereal,gradt) !  ,h0,h1,ee,ev,w,natoms))
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      save
        DOUBLE PRECISION x(3,n),h0(9*n*(n-1)/2,9*n*(n-1)/2)
     &        ,h1(9*n*(n-1)/2,9*n*(n-1)/2)
     &        ,grad(n*3),ee(9*n*(n-1)/2),ev(9*n*(n-1)/2,9*n*(n-1)/2)
     &        ,w(9*n*(n-1)/2,1)
      logical gradt
        INTEGER INFO, NCALL
        DOUBLE PRECISION TEMPA(9*N)
      data iev/1/, ncall/1/, shift/10.0d0/
      n3=9*n*(n-1)/2

        id=1
        it=0
      call hmatd(x,n,n3,h0,h1,h2,0)
      if(ncall.eq.1.or.it.eq.0) then
         if(id.eq.0) then
C            call jacob2(h0,ev,w(1,1),ee,w(1,2),n3,iev)
                     else
C       call tred2(h0,n3,n3,ee,w(1,1),iev)
C       call tqli(ee,w(1,1),n3,n3,h0,iev)

        CALL DSYEV('V','U',N3,H0,N3,EE,TEMPA,9*N,INFO)
        IF (INFO.NE.0) THEN
           PRINT*,'WARNING - INFO=',INFO,' in DSYEV'
        ENDIF

      do i=1,n3
      do j=1,n3
      ev(i,j)=h0(i,j)
      enddo
      enddo
                     end if
      call sortv(ee,ev,n3,iev)
      ereal=ee(n3)
      ncall=2
                                else
      do i=1,n3
      h0(i,i)=h0(i,i)-shift
      enddo
C       call iteig(h0,ev(1,n3),ereal,w(1,1),w(1,2),n3)
      ereal=ereal+shift
                                end if

      if(gradt) then

      do i=1,n*3
      call hmatd(x,n,n3,h0,h1,h2,i)
      vhvi=0.0d0
      do j=1,n3
      vh=0.0d0
      do k=1,n3
      vh=vh+ev(k,n3)*h1(k,j)
      enddo
      vhvi=vhvi+vh*ev(j,n3)
      enddo
      grad(i)=vhvi
      enddo

              end if
      return
      end
c ------------------------------------------------------------
c      real*8 function energy(r0,i,ip)
c      real*8 r0
c      r=r0
c      energy=ut(r,i,ip)
c      return
c      end
      subroutine hmatd(x,n,n3,h0,h1,h2,igrad)
      USE commons
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER N, N3
      save
C        common/test/itest
      DOUBLE PRECISION x(3,n),h0(9*n*(n-1)/2,9*n*(n-1)/2)
     1               ,h1(9*n*(n-1)/2,9*n*(n-1)/2)
      INTEGER J0, K, KI, KI0, IJKD, NN9, KJ0, IR, IJK0, IGRAD, I, IJ, IJ0, IJK, J

      data RauA/0.52918d0/ ,Runit/1.d0/
c           isu/1/,isg/2/,ipu/3/,ipg/4/,ix/5/


C       IF (9*n*(n-1)/2.GT.3*MXATMS) THEN
C          PRINT*,'Increase MXATMS for hmatd'
C          STOP
C       ENDIF
        IF (NEON) Runit=RauA
      
      if(igrad.ne.0) go to 10

      nn9=9*n*(n-1)/2
      do i=1,nn9
      do j=1,nn9
      h0(i,j)=0.0d0
      enddo
      enddo

      ExS=0.0d0
      do k=1,n-1
      do l=k+1,n
      rkl=
     1      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
cc      Ex=energy(rkl,ix,0)
        IF (NEON) THEN
           Ex=utn(rkl/RauA,5,0)
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        ENDIF
      ExS=ExS+Ex
      enddo
      enddo

      do i=1,n-1
      do j=i+1,n

      ij = (2*n-i)*(i-1)/2 +j-i
      rij=
     1      sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2+(x(3,i)-x(3,j))**2)
cc      Exij=energy(rij,ix,0)
        IF (NEON) THEN
           Exij=utn(rij/RauA,5,0)
        ELSE
           CALL ASAR1(rij,Exij,Exij1)
        ENDIF

      do ki=1,3
      do kj=1,3

      ijk = 9*(ij-1) +3*(ki-1) +kj
      h0(ijk,ijk) = h0(ijk,ijk) +ExS -Exij +Runit/rij

      if(j.eq.n) go to 2
      do j0=j+1,n

      dx=x(1,j0)-x(1,j)
      dy=x(2,j0)-x(2,j)
      dz=x(3,j0)-x(3,j)
      dxy=sqrt(dx*dx+dy*dy)
      rjj0=sqrt(dxy*dxy+dz*dz)
      ct=dz/rjj0
      st=dxy/rjj0
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if
        IF (NEON) THEN
           call fenergy(rjj0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rjj0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rjj0,isu,0)
cc      Esg=energy(rjj0,isg,0)
cc      Epu=energy(rjj0,ipu,0)
cc      Epg=energy(rjj0,ipg,0)
      Ess=0.5*(Esu+Esg)
      Esd=0.5*(Esu-Esg)
      Eps=0.5*(Epu+Epg)
      Epd=0.5*(Epu-Epg)
cc      Ex=energy(rjj0,ix,0)
      
      ij0 = (2*n-i)*(i-1)/2 +j0-i

      do kj0=1,3

      ijk0 = 9*(ij0-1) +3*(ki-1) +kj0

      if(kj.eq.1) then
      if(kj0.eq.1) then
      hs=Ess*(st*cf)**2+Eps*((ct*cf)**2+sf**2)
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(st*cf)**2+Epd*((ct*cf)**2+sf**2)      
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.2) then
        hs=(Ess-Eps)*st**2*sf*cf
        h0(ijk,ijk+1)=h0(ijk,ijk+1)+hs
        h0(ijk+1,ijk)=h0(ijk+1,ijk)+hs
      hd=(Esd-Epd)*st**2*sf*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.3) then
        hs=(Ess-Eps)*st*ct*cf
        h0(ijk,ijk+2)=h0(ijk,ijk+2)+hs
        h0(ijk+2,ijk)=h0(ijk+2,ijk)+hs
      hd=(Esd-Epd)*st*ct*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      if(kj.eq.2) then
      if(kj0.eq.2) then
      hs=Ess*(st*sf)**2+Eps*((ct*sf)**2+cf**2)
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(st*sf)**2+Epd*((ct*sf)**2+cf**2)
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.1) then
        hs=(Ess-Eps)*st**2*sf*cf
        h0(ijk0,ijk0+1)=h0(ijk0,ijk0+1)+hs
        h0(ijk0+1,ijk0)=h0(ijk0+1,ijk0)+hs
      hd=(Esd-Epd)*st**2*sf*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.3) then
        hs=(Ess-Eps)*st*ct*sf
        h0(ijk,ijk+1)=h0(ijk,ijk+1)+hs
        h0(ijk+1,ijk)=h0(ijk+1,ijk)+hs
      hd=(Esd-Epd)*st*ct*sf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      if(kj.eq.3) then
      if(kj0.eq.3) then
      hs=Ess*(ct)**2   +Eps*(st)**2
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(ct)**2   +Epd*(st)**2
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.1) then
        hs=(Ess-Eps)*st*ct*cf
        h0(ijk0,ijk0+2)=h0(ijk0,ijk0+2)+hs
        h0(ijk0+2,ijk0)=h0(ijk0+2,ijk0)+hs
      hd=(Esd-Epd)*st*ct*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(kj0.eq.2) then
        hs=(Ess-Eps)*st*ct*sf
        h0(ijk0,ijk0+1)=h0(ijk0,ijk0+1)+hs
        h0(ijk0+1,ijk0)=h0(ijk0+1,ijk0)+hs
      hd=(Esd-Epd)*st*ct*sf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      enddo
      enddo

   2      continue
      do i0=i+1,n
      if(i0.eq.j) go to 1

      dx=x(1,i0)-x(1,i)
      dy=x(2,i0)-x(2,i)
      dz=x(3,i0)-x(3,i)
      dxy=sqrt(dx*dx+dy*dy)
      rii0=sqrt(dxy*dxy+dz*dz)
      ct=dz/rii0
      st=dxy/rii0
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if
        IF (NEON) THEN
           call fenergy(rii0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rii0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rii0,isu,0)
cc      Esg=energy(rii0,isg,0)
cc      Epu=energy(rii0,ipu,0)
cc      Epg=energy(rii0,ipg,0)
      Ess=0.5*(Esu+Esg)
      Esd=0.5*(Esu-Esg)
      Eps=0.5*(Epu+Epg)
      Epd=0.5*(Epu-Epg)
cc      Ex=energy(rii0,ix,0)
      
      if(i0.lt.j) ij0 = (2*n-i0)*(i0-1)/2 +j-i0
      if(i0.gt.j) ij0 = (2*n-j)*(j-1)/2 +i0-j

      do ki0=1,3

      if(i0.lt.j) ijk0 = 9*(ij0-1) +3*(ki0-1) +kj
      if(i0.gt.j) ijk0 = 9*(ij0-1) +3*(kj-1) +ki0

      if(ki.eq.1) then
      if(ki0.eq.1) then
      hs=Ess*(st*cf)**2+Eps*((ct*cf)**2+sf**2)
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(st*cf)**2+Epd*((ct*cf)**2+sf**2)      
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.2) then
        hs=(Ess-Eps)*st**2*sf*cf
        h0(ijk,ijk+3)=h0(ijk,ijk+3)+hs
        h0(ijk+3,ijk)=h0(ijk+3,ijk)+hs
      hd=(Esd-Epd)*st**2*sf*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.3) then
        hs=(Ess-Eps)*st*ct*cf
        h0(ijk,ijk+6)=h0(ijk,ijk+6)+hs
        h0(ijk+6,ijk)=h0(ijk+6,ijk)+hs
      hd=(Esd-Epd)*st*ct*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      if(ki.eq.2) then
      if(ki0.eq.2) then
      hs=Ess*(st*sf)**2+Eps*((ct*sf)**2+cf**2)
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(st*sf)**2+Epd*((ct*sf)**2+cf**2)
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.1) then
        hs=(Ess-Eps)*st**2*sf*cf
      ijkd=1
      if(i0.lt.j) ijkd=3
        h0(ijk0,ijk0+ijkd)=h0(ijk0,ijk0+ijkd)+hs
        h0(ijk0+ijkd,ijk0)=h0(ijk0+ijkd,ijk0)+hs
      hd=(Esd-Epd)*st**2*sf*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.3) then
        hs=(Ess-Eps)*st*ct*sf
        h0(ijk,ijk+3)=h0(ijk,ijk+3)+hs
        h0(ijk+3,ijk)=h0(ijk+3,ijk)+hs
      hd=(Esd-Epd)*st*ct*sf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      if(ki.eq.3) then
      if(ki0.eq.3) then
      hs=Ess*(ct)**2   +Eps*(st)**2
      h0(ijk,ijk)=h0(ijk,ijk)+hs-Ex
      h0(ijk0,ijk0)=h0(ijk0,ijk0)+hs-Ex
      hd=Esd*(ct)**2   +Epd*(st)**2
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.1) then
        hs=(Ess-Eps)*st*ct*cf
      ijkd=2
      if(i0.lt.j) ijkd=6
        h0(ijk0,ijk0+ijkd)=h0(ijk0,ijk0+ijkd)+hs
        h0(ijk0+ijkd,ijk0)=h0(ijk0+ijkd,ijk0)+hs
      hd=(Esd-Epd)*st*ct*cf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
      if(ki0.eq.2) then
        hs=(Ess-Eps)*st*ct*sf
      ijkd=1
      if(i0.lt.j) ijkd=3
        h0(ijk0,ijk0+ijkd)=h0(ijk0,ijk0+ijkd)+hs
        h0(ijk0+ijkd,ijk0)=h0(ijk0+ijkd,ijk0)+hs
      hd=(Esd-Epd)*st*ct*sf
      h0(ijk,ijk0)=hd
      h0(ijk0,ijk)=hd
                    end if
                    end if
      enddo
   1      continue
      enddo

      enddo
      enddo

      enddo
      enddo

c ---------------------------------------------------------------------

  10      continue
      if(igrad.le.0) go to 20

      nn9=9*n*(n-1)/2
      do i=1,nn9
      do j=1,nn9
      h1(i,j)=0.0d0
      enddo
      enddo

      m=int(igrad/3)+1
      ir=mod(igrad,3)
      if(ir.eq.0) then
      ir=3
      m=m-1
                  end if

      ExS1=0.0d0
      do k=1,n-1
      do l=k+1,n
      if(l.eq.m.or.k.eq.m) then
      rkl=
     1      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
cc      Ex1= energy(rkl,ix,1)*(x(ir,l)-x(ir,k))/rkl
        IF (NEON) THEN
        Ex1= utn(rkl/RauA,5,1)/RauA*(x(ir,l)-x(ir,k))/rkl
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        Ex1= Ex1*(x(ir,l)-x(ir,k))/rkl
        ENDIF
      if(k.eq.m) Ex1=-Ex1
      ExS1=ExS1+Ex1
                       end if
      enddo
      enddo

      do i=1,n-1
      do j=i+1,n

      ij = (2*n-i)*(i-1)/2 +j-i

      if(i.eq.m.or.j.eq.m) then

      rij=
     1      sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2+(x(3,i)-x(3,j))**2)
      drijdx=(x(ir,i)-x(ir,j))/rij
      if(j.eq.m) drijdx=-drijdx
cc      Exij1= energy(rij,ix,1)*drijdx
        IF (NEON) THEN
        Exij1= utn(rij/RauA,5,1)/RauA*drijdx
        ELSE
           CALL ASAR1(rij,Exij,Exij1)
        Exij1= Exij1*drijdx
        ENDIF

                           end if

      do ki=1,3
      do kj=1,3

      ijk = 9*(ij-1) +3*(ki-1) +kj
      h1(ijk,ijk) = h1(ijk,ijk) +ExS1

      if(i.eq.m.or.j.eq.m) then

      h1(ijk,ijk) = h1(ijk,ijk) -Exij1 -Runit/(rij*rij)*drijdx

                           end if

      if(j.eq.n) go to 22
      do j0=j+1,n

      if(j.eq.m.or.j0.eq.m) then

      dx=x(1,j0)-x(1,j)
      dy=x(2,j0)-x(2,j)
      dz=x(3,j0)-x(3,j)
      dxy=sqrt(dx*dx+dy*dy)
      rjj0=sqrt(dxy*dxy+dz*dz)
      ct=dz/rjj0
      st=dxy/rjj0
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if

      drdx= (x(ir,j0)-x(ir,j))/rjj0
      if(j.eq.m) drdx=-drdx

      ct1=-ct/rjj0*drdx
      if(ir.eq.3) then
      if(m.eq.j0) ct1=ct1+1.0d0/rjj0
      if(m.eq.j) ct1=ct1-1.0d0/rjj0
                end if
      st1=-st/rjj0*drdx
      dxydx=0.0d0
      if(ir.eq.1.or.ir.eq.2) then
c      dxydx=1.0d0
      if(dxy.gt.0.0d0) then
      dxydx= (x(ir,j0)-x(ir,j))/dxy
      if(j.eq.m) dxydx=-dxydx
                   end if
      st1=st1+dxydx/rjj0
                             end if
      cf1=0.0d0
      sf1=0.0d0
      if(dxy.gt.0.0d0) then
      if(ir.eq.1.or.ir.eq.2) then
      cf1=-cf/dxy*dxydx
      if(ir.eq.1) then
      if(m.eq.j0) cf1=cf1+1.0d0/dxy
      if(m.eq.j) cf1=cf1-1.0d0/dxy
                    end if
      sf1=-sf/dxy*dxydx
      if(ir.eq.2) then
      if(m.eq.j0) sf1=sf1+1.0d0/dxy
      if(m.eq.j) sf1=sf1-1.0d0/dxy
                    end if
                             end if
                   end if

        IF (NEON) THEN
           call fenergy(rjj0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rjj0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rjj0,isu,0)
cc      Esg=energy(rjj0,isg,0)
cc      Epu=energy(rjj0,ipu,0)
cc      Epg=energy(rjj0,ipg,0)
      Ess=0.5*(Esu+Esg)
      Esd=0.5*(Esu-Esg)
      Eps=0.5*(Epu+Epg)
      Epd=0.5*(Epu-Epg)
cc      Ex=energy(rjj0,ix,0)

cc      Esu1=energy(rjj0,isu,1)*drdx
cc      Esg1=energy(rjj0,isg,1)*drdx
cc      Epu1=energy(rjj0,ipu,1)*drdx
cc      Epg1=energy(rjj0,ipg,1)*drdx
        Esu1=Esu1*drdx
        Esg1=Esg1*drdx
        Epu1=Epu1*drdx
        Epg1=Epg1*drdx
      Ess1=0.5*(Esu1+Esg1)
      Esd1=0.5*(Esu1-Esg1)
      Eps1=0.5*(Epu1+Epg1)
      Epd1=0.5*(Epu1-Epg1)
cc      Ex1=energy(rjj0,ix,1)*drdx
      Ex1=Ex1*drdx
      
      ij0 = (2*n-i)*(i-1)/2 +j0-i

      do kj0=1,3

      ijk0 = 9*(ij0-1) +3*(ki-1) +kj0

      if(kj.eq.1) then
      if(kj0.eq.1) then
           trig1=(st*cf)**2
           trig2=(ct*cf)**2+sf**2
           trig3=2.0d0*st*cf*(cf*st1+st*cf1)
           trig4=(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
cc      hs1=Ess1*(st*cf)**2+Eps1*((ct*cf)**2+sf**2) +
cc     1          Ess*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     1      Eps*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
        hs1= Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(st*cf)**2+Epd1*((ct*cf)**2+sf**2) +
cc     1          Esd*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     1      Epd*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
        hd1= Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.2) then
           trig1=st**2*sf*cf
           trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     1      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)) 
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+1)=h1(ijk,ijk+1)+hs1
        h1(ijk+1,ijk)=h1(ijk+1,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     1      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.3) then
           trig1=st*ct*cf
           trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     1      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+2)=h1(ijk,ijk+2)+hs1
        h1(ijk+2,ijk)=h1(ijk+2,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     1      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      if(kj.eq.2) then
      if(kj0.eq.2) then
           trig1=(st*sf)**2
           trig2=(ct*sf)**2+cf**2
           trig3=2.0d0*st*sf*(sf*st1+st*sf1)
           trig4=2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1
cc      hs1=Ess1*(st*sf)**2+Eps1*((ct*sf)**2+cf**2) +
cc     1          Ess*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     1      Eps*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
        hs1=Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(st*sf)**2+Epd1*((ct*sf)**2+cf**2) +
cc     1          Esd*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     1      Epd*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
        hd1=Esd1*trig1+Epd1*trig2+Esd*trig3+Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.1) then
           trig1=st**2*sf*cf
           trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     1      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk0,ijk0+1)=h1(ijk0,ijk0+1)+hs1
        h1(ijk0+1,ijk0)=h1(ijk0+1,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     1      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.3) then
           trig1=st*ct*sf
           trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
cc        hs1=(Ess1-Eps1)*st*ct*sf +
cc     1      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+1)=h1(ijk,ijk+1)+hs1
        h1(ijk+1,ijk)=h1(ijk+1,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     1      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      if(kj.eq.3) then
      if(kj0.eq.3) then
           trig1=(ct)**2
           trig2=(st)**2
           trig3=2.0d0*ct*ct1
           trig4=2.0d0*st*st1
cc      hs1=Ess1*(ct)**2   +Eps1*(st)**2 +
cc     1      Ess*2.0d0*ct*ct1   +Eps*2.0d0*st*st1
        hs1=Ess1*trig1 + Eps1*trig2 + Ess*trig3 + Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(ct)**2   +Epd1*(st)**2 +
cc     1      Esd*2.0d0*ct*ct1   +Epd*2.0d0*st*st1
        hd1=Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.1) then
           trig1=st*ct*cf
           trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     1      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk0,ijk0+2)=h1(ijk0,ijk0+2)+hs1
        h1(ijk0+2,ijk0)=h1(ijk0+2,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     1      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(kj0.eq.2) then
           trig1=st*ct*sf
           trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
cc        hs1=(Ess1-Eps1)*st*ct*sf + 
cc     1      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk0,ijk0+1)=h1(ijk0,ijk0+1)+hs1
        h1(ijk0+1,ijk0)=h1(ijk0+1,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     1      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      enddo
                       end if
      enddo

  22      continue
      do i0=i+1,n
      if(i0.eq.j) go to 11

      if(i.eq.m.or.i0.eq.m) then

      dx=x(1,i0)-x(1,i)
      dy=x(2,i0)-x(2,i)
      dz=x(3,i0)-x(3,i)
      dxy=sqrt(dx*dx+dy*dy)
      rii0=sqrt(dxy*dxy+dz*dz)
      ct=dz/rii0
      st=dxy/rii0
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if

      drdx= (x(ir,i0)-x(ir,i))/rii0
      if(i.eq.m) drdx=-drdx

      ct1=-ct/rii0*drdx
      if(ir.eq.3) then
      if(m.eq.i0) ct1=ct1+1.0d0/rii0
      if(m.eq.i) ct1=ct1-1.0d0/rii0
                end if
      st1=-st/rii0*drdx
      dxydx=0.0d0
      if(ir.eq.1.or.ir.eq.2) then
c      dxydx=1.0d0
      if(dxy.gt.0.0d0) then
      dxydx= (x(ir,i0)-x(ir,i))/dxy
      if(i.eq.m) dxydx=-dxydx
                   end if
      st1=st1+dxydx/rii0
                             end if
      cf1=0.0d0
      sf1=0.0d0
      if(dxy.gt.0.0d0) then
      if(ir.eq.1.or.ir.eq.2) then
      cf1=-cf/dxy*dxydx
      if(ir.eq.1) then
      if(m.eq.i0) cf1=cf1+1.0d0/dxy
      if(m.eq.i) cf1=cf1-1.0d0/dxy
                    end if
      sf1=-sf/dxy*dxydx
      if(ir.eq.2) then
      if(m.eq.i0) sf1=sf1+1.0d0/dxy
      if(m.eq.i) sf1=sf1-1.0d0/dxy
                    end if
                             end if
                   end if

        IF (NEON) THEN
           call fenergy(rii0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rii0,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rii0,isu,0)
cc      Esg=energy(rii0,isg,0)
cc      Epu=energy(rii0,ipu,0)
cc      Epg=energy(rii0,ipg,0)
      Ess=0.5*(Esu+Esg)
      Esd=0.5*(Esu-Esg)
      Eps=0.5*(Epu+Epg)
      Epd=0.5*(Epu-Epg)
cc      Ex=energy(rii0,ix,0)

cc      Esu1=energy(rii0,isu,1)*drdx
cc      Esg1=energy(rii0,isg,1)*drdx
cc      Epu1=energy(rii0,ipu,1)*drdx
cc      Epg1=energy(rii0,ipg,1)*drdx
        Esu1=Esu1*drdx
        Esg1=Esg1*drdx
        Epu1=Epu1*drdx
        Epg1=Epg1*drdx
      Ess1=0.5*(Esu1+Esg1)
      Esd1=0.5*(Esu1-Esg1)
      Eps1=0.5*(Epu1+Epg1)
      Epd1=0.5*(Epu1-Epg1)
cc      Ex1=energy(rii0,ix,1)*drdx
      Ex1=Ex1*drdx      

      if(i0.lt.j) ij0 = (2*n-i0)*(i0-1)/2 +j-i0
      if(i0.gt.j) ij0 = (2*n-j)*(j-1)/2 +i0-j

      do ki0=1,3

      if(i0.lt.j) ijk0 = 9*(ij0-1) +3*(ki0-1) +kj
      if(i0.gt.j) ijk0 = 9*(ij0-1) +3*(kj-1) +ki0

      if(ki.eq.1) then
      if(ki0.eq.1) then
           trig1=(st*cf)**2
           trig2=(ct*cf)**2+sf**2
           trig3=2.0d0*st*cf*(cf*st1+st*cf1)
           trig4=(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
cc      hs1=Ess1*(st*cf)**2+Eps1*((ct*cf)**2+sf**2) +
cc     1          Ess*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     1      Eps*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
        hs1= Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(st*cf)**2+Epd1*((ct*cf)**2+sf**2) +
cc     1          Esd*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     1      Epd*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
        hd1= Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.2) then
           trig1=st**2*sf*cf
           trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     1      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)) 
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+3)=h1(ijk,ijk+3)+hs1
        h1(ijk+3,ijk)=h1(ijk+3,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     1      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.3) then
           trig1=st*ct*cf
           trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     1      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+6)=h1(ijk,ijk+6)+hs1
        h1(ijk+6,ijk)=h1(ijk+6,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     1      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      if(ki.eq.2) then
      if(ki0.eq.2) then
           trig1=(st*sf)**2
           trig2=(ct*sf)**2+cf**2
           trig3=2.0d0*st*sf*(sf*st1+st*sf1)
           trig4=2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1
cc      hs1=Ess1*(st*sf)**2+Eps1*((ct*sf)**2+cf**2) +
cc     1          Ess*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     1      Eps*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
        hs1=Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(st*sf)**2+Epd1*((ct*sf)**2+cf**2) +
cc     1          Esd*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     1      Epd*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
        hd1=Esd1*trig1+Epd1*trig2+Esd*trig3+Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.1) then
           trig1=st**2*sf*cf
           trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     1      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
      ijkd=1
      if(i0.lt.j) ijkd=3
        h1(ijk0,ijk0+ijkd)=h1(ijk0,ijk0+ijkd)+hs1
        h1(ijk0+ijkd,ijk0)=h1(ijk0+ijkd,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     1      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.3) then
           trig1=st*ct*sf
           trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
cc        hs1=(Ess1-Eps1)*st*ct*sf +
cc     1      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(ijk,ijk+3)=h1(ijk,ijk+3)+hs1
        h1(ijk+3,ijk)=h1(ijk+3,ijk)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     1      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      if(ki.eq.3) then
      if(ki0.eq.3) then
           trig1=(ct)**2
           trig2=(st)**2
           trig3=2.0d0*ct*ct1
           trig4=2.0d0*st*st1
cc      hs1=Ess1*(ct)**2   +Eps1*(st)**2 +
cc     1      Ess*2.0d0*ct*ct1   +Eps*2.0d0*st*st1
        hs1=Ess1*trig1 + Eps1*trig2 + Ess*trig3 + Eps*trig4
      h1(ijk,ijk)=h1(ijk,ijk)+hs1-Ex1
      h1(ijk0,ijk0)=h1(ijk0,ijk0)+hs1-Ex1
cc      hd1=Esd1*(ct)**2   +Epd1*(st)**2 +
cc     1      Esd*2.0d0*ct*ct1   +Epd*2.0d0*st*st1
        hd1=Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.1) then
           trig1=st*ct*cf
           trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     1      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
      ijkd=2
      if(i0.lt.j) ijkd=6
        h1(ijk0,ijk0+ijkd)=h1(ijk0,ijk0+ijkd)+hs1
        h1(ijk0+ijkd,ijk0)=h1(ijk0+ijkd,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     1      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
      if(ki0.eq.2) then
           trig1=st*ct*sf
           trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
cc        hs1=(Ess1-Eps1)*st*ct*sf + 
cc     1      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
      ijkd=1
      if(i0.lt.j) ijkd=3
        h1(ijk0,ijk0+ijkd)=h1(ijk0,ijk0+ijkd)+hs1
        h1(ijk0+ijkd,ijk0)=h1(ijk0+ijkd,ijk0)+hs1
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     1      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(ijk,ijk0)=hd1
      h1(ijk0,ijk)=hd1
                    end if
                    end if

      enddo
                       end if
  11      continue
      enddo

      enddo
      enddo

      enddo
      enddo

  20      continue

      return
      end
      subroutine grnd(n,x,grad,ereal,gradt)
      USE commons
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION x(3,n),grad(n*3)
      logical gradt

      save
      data RauA/0.52918d0/

      ExS=0.0d0
      do k=1,n-1
      do l=k+1,n
      rkl=
     1      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
        IF (NEON) THEN
           Ex=utn(rkl/RauA,5,0)
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        ENDIF
      ExS=ExS+Ex
      enddo
      enddo
      ereal=ExS


      if(gradt) then

      do m=1,n
      do ir=1,3
      ic=(m-1)*3+ir

      ExS1=0.0d0
      do k=1,n-1
      do l=k+1,n
      if(l.eq.m.or.k.eq.m) then
      rkl=
     1      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
        IF (NEON) THEN
      Ex1= utn(rkl/RauA,5,1)/RauA*(x(ir,l)-x(ir,k))/rkl
      if(k.eq.m) Ex1=-Ex1
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        Ex1= Ex1*(x(ir,l)-x(ir,k))/rkl
        if(k.eq.m) Ex1=-Ex1
        ENDIF
      ExS1=ExS1+Ex1
                       end if
      enddo
      enddo

      grad(ic)=ExS1

      enddo
      enddo

                  end if

      return
      end
             SUBROUTINE CUBSPL(TAU,C,N,IBCBEG,IBCEND)
      implicit double precision (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER IBCBEG,IBCEND,N,I,J,L,M,J1
      DOUBLE PRECISION C(4,N),TAU(N),DIVDF1,DIVDF3,DTAU,G
      L=N-1
      DO M=2,N
        C(3,M)=TAU(M)-TAU(M-1)
        C(4,M)=(C(1,M)-C(1,M-1))/C(3,M)
      END DO
      IF(IBCBEG-1 .LT. 0) THEN
        GOTO 11
      ELSE IF (IBCBEG-1 .EQ. 0) THEN
        GOTO 15
      ELSE 
        GOTO 16
      END IF
   11 IF(N.GT.2) GO TO 12
      C(4,1)=1.d0
      C(3,1)=1.d0
      C(2,1)=2.d0*C(4,2)
      GO TO 25
   12 C(4,1)=C(3,3)
      C(3,1)=C(3,2)+C(3,3)
      C(2,1)=((C(3,2)+2.d0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))
     1      /C(3,1)
      GO TO 19
   15 C(4,1)=1.d0
      C(3,1)=0.d0
      GO TO 18
   16 C(4,1)=2.d0
      C(3,1)=1.d0
      C(2,1)=3.d0*C(4,2)-C(3,2)/2.d0*C(2,1)
   18 IF(N.EQ.2) GO TO 25
   19 DO M=2,L
        G=-C(3,M+1)/C(4,M-1)
        C(2,M)=G*C(2,M-1)+3.d0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
        C(4,M)=G*C(3,M-1)+2.d0*(C(3,M)+C(3,M+1))
      END DO
      IF(IBCEND-1 .LT. 0) THEN
        GOTO 21
      ELSE IF(IBCEND-1 .EQ. 0) THEN
        GOTO 30
      ELSE
        GOTO 24
      END IF
   21 IF(N.EQ.3.AND.IBCBEG.EQ.0) GO TO 22
      G=C(3,N-1)+C(3,N)
      C(2,N)=( (C(3,N)+2.d0*G)*C(4,N)*C(3,N-1)+C(3,N)**2*
     1         (C(1,N-1)-C(1,N-2))/C(3,N-1) )/G
      G=-G/C(4,N-1)
      C(4,N)=C(3,N-1)
      GO TO 29
   22 C(2,N)=2.d0*C(4,N)
      C(4,N)=1.d0
      GO TO 28
   24 C(2,N)=3.d0*C(4,N)+C(3,N)/2.d0*C(2,N)
      C(4,N)=2.d0
      GO TO 28
   25 IF(IBCEND-1 .LT. 0) THEN
        GOTO 26
      ELSE IF(IBCEND-1 .EQ. 0) THEN
        GOTO 30
      ELSE
        GOTO 24
      END IF
   26 IF(IBCBEG.GT.0) GO TO 22
      C(2,N)=C(4,N)
      GO TO 30
   28 G=-1.d0/C(4,N-1)
   29 C(4,N)=G*C(3,N-1)+C(4,N)
      C(2,N)=(G*C(2,N-1)+C(2,N))/C(4,N)
   30 DO J1=1,L
        J=L+1-J1
        C(2,J)=(C(2,J)-C(3,J)*C(2,J+1))/C(4,J)
      END DO
      DO I=2,N
        DTAU=C(3,I)
        DIVDF1=(C(1,I)-C(1,I-1))/DTAU
        DIVDF3=C(2,I-1)+C(2,I)-2.d0*DIVDF1
        C(3,I-1)=2.d0*(DIVDF1-C(2,I-1)-DIVDF3)/DTAU
        C(4,I-1)=(DIVDF3/DTAU)*(6.d0/DTAU)
      END DO
      RETURN
      END
      SUBROUTINE INTERV( XT, LXT, X, LEFT, MFLAG)
COMPUTES LEFT = MAX( I, 1 .LE. I .LE. LXT .AND. XT(I) .LE. X)
C
C******  I N P U T  ******
C  XT  ... A REAL SEQUENCE, OF LENGTH  LXT, ASSUMED TO BE NONDECREASING
C  LXT ... NUMBER OF TERMS IN THE SEQUENCE  XT.
C  X   ... THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE  XT
C          IS TO BE DETERMINED.
C
C******  O U T P U T  ******
C  LEFT, MFLAG ... BOTH INTEGER, WHOSE VALUE IS
C
C   1    -1      IF  XT(LXT) .LE. X .LT.  XT(1)
C   1     0      IF   XT(I)  .LE. X .LT. XT(I+1)
C  LXT    1      IF  XT(LXT) .LE. X
C
C     IN PARTICULAR, MFLAG = 0 IS THE 'USUAL' CASE. MFLAG .NE. 0
C     INDICATES THAT  X  LIES OUTSIDE THE HALFOPEN INTERVAL
C     XT(1) .LE. Y .LT. XT(LXT) . THE ASSYMETRIC TREATMENT OF THE
C     INTERVAL IS DUE TO THE DECISION TO MAKE ALL  PP  FUNCTIONS CONT-
C     INUOUS FROM THE RIGHT.
C
C******  M E T H O D  ******
C  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT
C  IT IS CALLED REPEATEDLY, WITH  X  TAKEN FROM AN INCREASING OR DECREA-
C  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE
C  GRAPHED. THE FIRST GUESS FOR  LEFT  IS THEREFORE TAKEN TO BE THE VAL-
C  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE  L O C A L  VARIA
C  BLE  ILO . A FIRST CHECK ASCERTAIN THAT  ILO .LT. LXT (THIS IS NEC-
C  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI-
C  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT =
C  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.
C     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO
C  WHILE ALSO MOVING  ILO  IND  IHI  IN THE DIRECTION OF  X , UNTIL
C             XT(ILO) .LE. X .LT. XT(IHI) ,
C  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .
C  LEFT = ILO IS THEN RETURNED.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER LEFT,LXT,MFLAG,
     *        IHI,ILO,ISTEP,MIDDLE
      DOUBLE PRECISION  X,XT(LXT)
      save
      DATA ILO /1/
C     SAVE ILO  (A VALID FORTRAN STATEMENT IN THE NEW 1977 STANDARD)
      IHI = ILO + 1
      IF (IHI .LT. LXT)                GO TO 20
        IF (X .GE. XT(LXT))            GO TO 110
          IF (LXT .LE. 1)              GO TO 90
            ILO = LXT - 1
            IHI = LXT
C
   20 IF (X .GE. XT(IHI))              GO TO 40
        IF (X .GE. XT(ILO))            GO TO 100
C      **** NOW X .LT. XT(ILO) . DECREASE  ILO  TO CAPTURE  X .
            ISTEP = 1
   31    IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .LE. 1)               GO TO 35
           IF (X .GE. XT(ILO))         GO TO 50
             ISTEP = ISTEP*2
                                       GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1))                GO TO 90
                                       GO TO 50
C     **** NOW X .GE. XT(IHI) . INCREASE  IHI TO CAPTURE  X .
   40 ISTEP = 1
   41   ILO = IHI
        IHI = ILO + ISTEP
        IF (IHI .GE. LXT)              GO TO 45
          IF (X .LT. XT(IHI))          GO TO 50
            ISTEP =ISTEP*2
                                       GO TO 41
   45 IF (X .GE. XT(LXT))              GO TO 110
        IHI = LXT
C
C       **** NOW XT(ILO) .LE. x .LT. XT(IHI) . NPRROW THE INTERVAL.
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO)             GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO + 1 .
         IF (X .LT. XT(MIDDLE))         GO TO 53
           ILO = MIDDLE
                                        GO TO 50
   53      IHI = MIDDLE
                                        GO TO 50
C**** SET OUTPUT AND RETURN.
   90 MFLAG = -1
      LEFT  =  1
                                        GO TO 99999
  100 MFLAG =  0
      LEFT  = ILO
                                        GO TO 99999
  110 MFLAG =  1
      LEFT  = LXT
99999                                   RETURN
      END

      DOUBLE PRECISION FUNCTION PPVALU( BREAK, COEF, L, K, X, JDERIV)
C   CALLS INTERV
CALCULATES VALUE AT  X  OF JDERIV-TH DERIVATIVE OF  PP  FROM PP-REPR
C
C******  I N P U T  ******
C  BREAK, COEF, L, K ... FORMS THE PP-REPRESENTATION OF THE FUNCTION  F
C         TO BE EVALUATED.SPECIFICALLY, THE J-TH DERIVATIVE OF  F  IS
C         GIVEN BY
C
C   D(**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I)+
C                            + H* COEF(  K,I) / (K-J-1))/(K-J-2) . )/2)
C
C   WITH  H = X - BREAK(I) , AND
C
C   I = MAX ( 1, MAX( J, BREAK(J) .LE. X, 1 .LE. J .LE. L ) ).
C
C  X      ... THE POINT AT WHICH TO EVALUATE.
C  JDERIV ... INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUATE-
C             ED. A S S U M E D  TO BE ZERO OR POSITIVE.
C
C******  O U T P U T  ******
C  PPVALU ... THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF  F  AT  X.
C
C******  M E T H O D  ******
C     THE INTERVAL INDEX  I, APPROPRIATE FOR X, IS FOUND THROUGH A
C  CALL TO INTERV. THE FORMULA ABOVE FOR THE  JDERIV-TH DERIVATIVE
C  OF  F  IS THEN EVALUATED (BY NESTED MULTIPLICATION).
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER JDERIV,K,L,
     *        I,M,NDUMMY,JDRP1
      DOUBLE PRECISION    BREAK(L),COEF(K,L)
c      DOUBLE PRECISION  FMMJDR,H,PPSUM,FDBLE,YD,X,FSNGL,YS
      PPSUM = 0.0
      FMMJDR = (FLOAT(K - JDERIV))
C     DERIVATIVES OF ORDER  K  OR HIGHER ARE IDENTICALLY ZERO.
      IF (FMMJDR .LE. 0.0)             GO TO 99
C     FIND INDEX  I  OF LARGEST BREAKPOINT TO THE LEFT OF  X.
      CALL INTERV( BREAK, L, X, I, NDUMMY)
C     EVALUATE JDERIV-TH DERIVATIVE OF  I-TH  POLYNOMIAL PIECE AT  X.
      H = X - BREAK(I)
      JDRP1 = JDERIV + 1
      M = K
      DO NDUMMY=JDRP1,K
         PPSUM = (PPSUM/FMMJDR)*H + COEF(M,I)
         FMMJDR = FMMJDR - 1.0
         M = M - 1
      END DO
   99 PPVALU = PPSUM
                                       RETURN
      END

        subroutine dipoles(x,n,h0,h1,igrad)
        implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
cc      save
      COMMON/DIPOL/ adip,Runit,ReX

        DOUBLE PRECISION x(3,n),h0(n*3,n*3)
     &               ,h1(n*3,n*3)

      if(igrad.gt.0) go to 10

      do i=1,n

        x1i=x(1,i)/Runit
        x2i=x(2,i)/Runit
        x3i=x(3,i)/Runit
      Eddi=0.d0

      do j=1,n-1
      if(j.eq.i) go to 1

      x1j=x(1,j)/Runit
        x2j=x(2,j)/Runit
        x3j=x(3,j)/Runit
      x1ij=x1j-x1i
      x2ij=x2j-x2i
      x3ij=x3j-x3i
      rij2= x1ij*x1ij + x2ij*x2ij + x3ij*x3ij
      rij=sqrt(rij2)
      Dij=adip/rij2

      do k=j+1,n
      if(k.eq.i.or.k.eq.j) go to 2 

        x1k=x(1,k)/Runit
        x2k=x(2,k)/Runit
        x3k=x(3,k)/Runit
        x1ik=x1k-x1i
        x2ik=x2k-x2i
        x3ik=x3k-x3i
        rik2= x1ik*x1ik + x2ik*x2ik + x3ik*x3ik        
        rik=sqrt(rik2) 
        Dik=adip/rik2

        x1jk=x1k-x1j                                                   
        x2jk=x2k-x2j                                                   
        x3jk=x3k-x3j                                                   
        rjk2= x1jk*x1jk + x2jk*x2jk + x3jk*x3jk                                 
      rjk=sqrt(rjk2)

      if(rjk.lt.0.5*ReX) go to 2

      DijDik=Dij*Dik*(x1ij*x1ik + x2ij*x2ik + x3ij*x3ik)/(rij*rik)
        Dijrjk=Dij*(x1ij*x1jk + x2ij*x2jk + x3ij*x3jk)/rij
        Dikrjk=Dik*(x1ik*x1jk + x2ik*x2jk + x3ik*x3jk)/rik          

      Eddjki=(rjk2*DijDik -3.d0*Dijrjk*Dikrjk)/(rjk2*rjk2*rjk)
      Eddi=Eddi+Eddjki

  2     continue
        enddo

  1      continue
        enddo

        i3=3*(i-1)
        do l=1,3
        i3l=i3+l
        h0(i3l,i3l)=h0(i3l,i3l)+Eddi
        enddo

      enddo

  10      continue
c----------------------------------------------------------------------
      if(igrad.le.0) go to 20

        m=int(igrad/3)+1
        ir=mod(igrad,3)
        if(ir.eq.0) then
        ir=3
        m=m-1
                    end if

        do i=1,n

        x1i=x(1,i)/Runit
        x2i=x(2,i)/Runit
        x3i=x(3,i)/Runit
        Eddi1=0.d0

        do j=1,n-1
        if(j.eq.i) go to 11

        x1j=x(1,j)/Runit
        x2j=x(2,j)/Runit
        x3j=x(3,j)/Runit
        x1ij=x1j-x1i
        x2ij=x2j-x2i
        x3ij=x3j-x3i
        rij2= x1ij*x1ij + x2ij*x2ij + x3ij*x3ij
        rij=sqrt(rij2)
        Dij=adip/rij2

      rij1=0.d0
      Dij1=0.d0
      if(i.eq.m.or.j.eq.m) then
      rij1=(x(ir,j)-x(ir,i))/rij/Runit
      if(i.eq.m) rij1=-rij1
      Dij1=-2.d0*Dij/rij*rij1
                       end if

        do k=j+1,n
        if(k.eq.i.or.k.eq.j) go to 22

        x1k=x(1,k)/Runit
        x2k=x(2,k)/Runit
        x3k=x(3,k)/Runit
        x1ik=x1k-x1i
        x2ik=x2k-x2i
        x3ik=x3k-x3i
        rik2= x1ik*x1ik + x2ik*x2ik + x3ik*x3ik
        rik=sqrt(rik2)
        Dik=adip/rik2

      rik1=0.d0
      Dik1=0.d0
        if(i.eq.m.or.k.eq.m) then
        rik1=(x(ir,k)-x(ir,i))/rik/Runit
        if(i.eq.m) rik1=-rik1      
        Dik1=-2.d0*Dik/rik*rik1
                             end if

        x1jk=x1k-x1j
        x2jk=x2k-x2j
        x3jk=x3k-x3j
        rjk2= x1jk*x1jk + x2jk*x2jk + x3jk*x3jk
        rjk=sqrt(rjk2)

      if (rjk.lt.0.5*ReX) go to 22

      xxjki=x1ij*x1ik + x2ij*x2ik + x3ij*x3ik
        DijDik=Dij*Dik*xxjki/(rij*rik)
      xxikj=x1ij*x1jk + x2ij*x2jk + x3ij*x3jk
        Dijrjk=Dij*xxikj/rij
      xxijk=x1ik*x1jk + x2ik*x2jk + x3ik*x3jk
        Dikrjk=Dik*xxijk/rik

      rjk1=0.d0
      if(j.eq.m.or.k.eq.m) then
        rjk1=(x(ir,k)-x(ir,j))/rjk/Runit
        if(j.eq.m) rjk1=-rjk1
                       end if

      xxjki1=0.d0
      DijDik1=0.d0
        if(i.eq.m.or.j.eq.m.or.k.eq.m) then
      if(j.eq.m) xxjki1=(x(ir,k)-x(ir,i))/Runit
        if(k.eq.m) xxjki1=(x(ir,j)-x(ir,i))/Runit
      if(i.eq.m) xxjki1=-(x(ir,k)-x(ir,i)+x(ir,j)-x(ir,i))/Runit
      DijDik1=(Dij1*Dik*xxjki+Dij*Dik1*xxjki+Dij*Dik*xxjki1)/(rij*rik)
     &             -DijDik*(rij1/rij+rik1/rik)
                                          end if
        xxikj1=0.d0
        Dijrjk1=0.d0
        if(i.eq.m.or.j.eq.m.or.k.eq.m) then
        if(i.eq.m) xxikj1=-(x(ir,k)-x(ir,j))/Runit
        if(k.eq.m) xxikj1=(x(ir,j)-x(ir,i))/Runit
        if(j.eq.m) xxikj1=(x(ir,k)-x(ir,j)-(x(ir,j)-x(ir,i)))/Runit
        Dijrjk1=(Dij1*xxikj+Dij*xxikj1)/rij
     &         -Dijrjk*rij1/rij
                                       end if
        xxijk1=0.d0  
        Dikrjk1=0.d0
        if(i.eq.m.or.j.eq.m.or.k.eq.m) then
        if(i.eq.m) xxijk1=-(x(ir,k)-x(ir,j))/Runit
        if(j.eq.m) xxijk1=-(x(ir,k)-x(ir,i))/Runit
        if(k.eq.m) xxijk1=(x(ir,k)-x(ir,j)+x(ir,k)-x(ir,i))/Runit 
        Dikrjk1=(Dik1*xxijk+Dik*xxijk1)/rik                             
     &         -Dikrjk*rik1/rik                            
                                       end if


        Eddjki=(rjk2*DijDik -3.d0*Dijrjk*Dikrjk)/(rjk2*rjk2*rjk)
        Eddjki1=(rjk2*DijDik1 -3.d0*(Dijrjk1*Dikrjk+Dijrjk*Dikrjk1)
     &             +2.d0*rjk*rjk1*DijDik)/(rjk2*rjk2*rjk)
     &             - 5.d0*Eddjki/rjk*rjk1
        Eddi1=Eddi1+Eddjki1

  22    continue
        enddo

  11    continue
        enddo

        i3=3*(i-1)
        do l=1,3
        i3l=i3+l
        h1(i3l,i3l)=h1(i3l,i3l)+Eddi1/Runit
        enddo

        enddo

  20    continue

      return
      end


c   2  2  2  2  2
c ne2+su.dat_qne2+sg.dat_qne2+pu.dat_qne2+pg.dat_qne2x.dat_q
c ne2+su.dat  ne2+sg.dat  ne2+pu.dat  ne2+pg.dat  ne2x.dat   

      subroutine fenergy(r,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)

      COMMON/DIPOL/ adip,Runit,rex

cc        real r0,ut,ut1
      save
      data RauA/0.52918d0/

      Runit=RauA
      adip=2.68d0

        r0=r
        Esu=ut(r0,1,ut1)
        Esu1=ut1
        Esg=ut(r0,2,ut1)
        Esg1=ut1
        Epu=ut(r0,3,ut1)
        Epu1=ut1
        Epg=ut(r0,4,ut1)
        Epg1=ut1
      Ex =utn(r/RauA,5,0)
      Ex1=utn(r/RauA,5,1)/RauA
        return
        end

      DOUBLE PRECISION function ut(rr,i0,ut1)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      save
      common/points/ru,Eu, rg,Eg, rpu,Epu, rpg,Epg, rx,Ex,
     1              coefu,coefg,coefpu,coefpg,coefx, ir1,ir2,ir3,ir4,ir5
      character(LEN=12) name1,name2,name3,name4,name5
      DOUBLE PRECISION ru(50),Eu(50), rg(50),Eg(50), 
     1     rpu(50),Epu(50), rpg(50),Epg(50), rx(50),Ex(50),
     1     coefu(4,50), coefg(4,50), coefpu(4,50), coefpg(4,50), 
     1     coefx(4,50)
      INTEGER IND1, IND2, IND3, IND4, IND5, NC, NC3, NC4, NC5, NC2, NC1, IR1, IR2, IR3, IR4, IR5
      data nc/1/,nc1/1/,nc2/1/,nc3/1/,nc4/1/,nc5/1/, rauA/0.52918d0/
      if(nc.eq.1) then
      nc=2
C     open(110,file='rg2+_3+n1.f')
!     open(110,file='energy.f')
      ind1=2  
      ind2=2  
      ind3=2  
      ind4=2  
      ind5=2
!     read(110,'(2x,5i3)') ind1,ind2,ind3,ind4,ind5
cc      print '(2x,5i3)', ind1,ind2,ind3,ind4,ind5
!     read(110,'(2x,5a12)') name1,name2,name3,name4,name5
cc      print '(2x,5a12)', name1,name2,name3,name4,name5
!      name1='ne2+su.dat'  
!      name2='ne2+sg.dat'
!      name3='ne2+pu.dat'
!      name4='ne2+pg.dat'
       name1='ne2+su.dat_q'  
       name2='ne2+sg.dat_q'
       name3='ne2+pu.dat_q'
       name4='ne2+pg.dat_q'
       name5='ne2x.dat'
                  end if
      ut=0.0
c
      if(i0.eq.1) then
      if(ind1.eq.0.or.ind1.eq.2) then
      if(nc1.eq.1) then
      nc1=2
      open(11,file=name1,STATUS='OLD')
      ir1=0
   10 continue
      ir1=ir1+1
      read(11,*,end=11) ru(ir1), Eu(ir1)
cc      print *, ru(ir1), Eu(ir1)
      go to 10
   11 continue
      ir1=ir1-1
      if(ind1.eq.2) then
      Easu=Eu(ir1)
      do i=1,ir1
      coefu(1,i)=Eu(i) -Easu
      enddo
      coefu(2,ir1)=0.0
      call cubspl(ru,coefu,ir1,0,1)
                    end if
                   end if
      if(rr.lt.ru(ir1)) then
      if(ind1.eq.0) then
      do 1 i=1,ir1-1
      if(.not.( ru(i).le.rr.and.rr.le.ru(i+1) )) go to 1
      ut = Eu(i) + (Eu(i+1)-Eu(i))/(ru(i+1)-ru(i))*(rr-ru(i)) -Easu
    1 continue
                    end if
      if(ind1.eq.2) then
      ut=ppvalu(ru,coefu,ir1,4,rr,0)
        ut1=ppvalu(ru,coefu,ir1,4,rr,1) 
                end if
                        else
      ut=Eu(ir1) -Easu
      ut1=0.0
                        end if
                    end if
      if(ind1.eq.1) ut=utn(rr/rauA,i0,0)
      return
                 end if
c
      if(i0.eq.2) then
      if(ind2.eq.0.or.ind2.eq.2) then
      if(nc2.eq.1) then
      nc2=2
      open(22,file=name2,STATUS='OLD')
      ir2=0
   20 continue
      ir2=ir2+1
      read(22,*,end=22) rg(ir2), Eg(ir2)
cc      print *, rg(ir2), Eg(ir2)
      go to 20
   22 continue
      ir2=ir2-1
      if(ind2.eq.2) then
      Easg=Eg(ir2)
      do i=1,ir2
      coefg(1,i)=Eg(i) -Easg
      enddo
      coefg(2,ir2)=0.0
      call cubspl(rg,coefg,ir2,0,1)
                    end if
                   end if
      if(rr.lt.rg(ir2)) then
      if(ind2.eq.0) then
      do 2 i=1,ir2-1
      if(.not.( rg(i).le.rr.and.rr.le.rg(i+1) )) go to 2
      ut = Eg(i) + (Eg(i+1)-Eg(i))/(rg(i+1)-rg(i))*(rr-rg(i)) -Easg
    2 continue
                    end if
      if(ind2.eq.2) then
      ut=ppvalu(rg,coefg,ir2,4,rr,0)
        ut1=ppvalu(rg,coefg,ir2,4,rr,1) 
                end if
                        else
      ut=Eg(ir2) -Easg
      ut1=0.0
                        end if
                    end if 
      if(ind2.eq.1) ut=utn(rr/rauA,i0,0)
      return
                 end if
c
      if(i0.eq.3) then
      if(ind3.eq.0.or.ind3.eq.2) then
      if(nc3.eq.1) then
      nc3=2
      open(33,file=name3,STATUS='OLD')
      ir3=0
   30 continue
      ir3=ir3+1
      read(33,*,end=33) rpu(ir3), Epu(ir3)
cc      print *, rpu(ir3), Epu(ir3)
      go to 30
   33 continue
      ir3=ir3-1
      if(ind3.eq.2) then
      Easpu=Epu(ir3)
      do i=1,ir3
      coefpu(1,i)=Epu(i) -Easpu
      enddo
      coefpu(2,ir3)=0.0
      call cubspl(rpu,coefpu,ir3,0,1)
                    end if
                   end if
      if(rr.lt.rpu(ir3)) then
      if(ind3.eq.0) then  
      do 3 i=1,ir3-1
      if(.not.( rpu(i).le.rr.and.rr.le.rpu(i+1) )) go to 3
      ut = Epu(i) + (Epu(i+1)-Epu(i))/(rpu(i+1)-rpu(i))*(rr-rpu(i))
     1     -Easpu
    3 continue
                    end if
      if(ind3.eq.2) then
      ut=ppvalu(rpu,coefpu,ir3,4,rr,0) 
        ut1=ppvalu(rpu,coefpu,ir3,4,rr,1) 
                end if
                        else
      ut=Epu(ir3) -Easpu
      ut1=0.0
                        end if
                    end if
      if(ind3.eq.1) ut=utn(rr/rauA,i0,0)
      return
               end if
c
      if(i0.eq.4) then
      if(ind4.eq.0.or.ind4.eq.2) then
      if(nc4.eq.1) then
      nc4=2
      open(44,file=name4,STATUS='OLD')
      ir4=0
   40 continue
      ir4=ir4+1
      read(44,*,end=44) rpg(ir4), Epg(ir4)
cc      print *, rpg(ir4), Epg(ir4)
      go to 40
   44 continue
      ir4=ir4-1
      if(ind4.eq.2) then
      Easpg=Epg(ir4)
      do i=1,ir4
      coefpg(1,i)=Epg(i) -Easpg
      enddo
      coefpg(2,ir4)=0.0
      call cubspl(rpg,coefpg,ir4,0,1)
                    end if
                   end if
      if(rr.lt.rpg(ir4)) then
      if(ind4.eq.0) then
      do 4 i=1,ir4-1
      if(.not.( rpg(i).le.rr.and.rr.le.rpg(i+1) )) go to 4
      ut = Epg(i) + (Epg(i+1)-Epg(i))/(rpg(i+1)-rpg(i))*(rr-rpg(i))
     1     -Easpg
    4 continue
                    end if
      if(ind4.eq.2) then
      ut=ppvalu(rpg,coefpg,ir4,4,rr,0)
        ut1=ppvalu(rpg,coefpg,ir4,4,rr,1)
                end if
                        else
      ut=Epg(ir4) -Easpg
      ut1=0.0
                        end if
                    end if 
      if(ind4.eq.1) ut=utn(rr/rauA,i0,0)
      return
                 end if
c
      if(i0.eq.5) then
      if(ind5.eq.0.or.ind5.eq.2) then
      if(nc5.eq.1) then
      nc5=2
      open(55,file=name5,STATUS='OLD')
      ir5=0
   50 continue
      ir5=ir5+1
      read(55,*,end=55) rx(ir5), Ex(ir5)
cc      print *, rx(ir5), Ex(ir5)
      go to 50
   55 continue
      ir5=ir5-1
      if(ind5.eq.2) then
      Easx=Ex(ir5)
      do i=1,ir5
      coefx(1,i)=Ex(i) -Easx
      enddo
      coefx(2,ir5)=0.0
      call cubspl(rx,coefx,ir5,0,1)
                    end if
                   end if
      if(rr.lt.rx(ir5)) then
      if(ind5.eq.0) then
      do 5 i=1,ir5-1
      if(.not.( rx(i).le.rr.and.rr.le.rx(i+1) )) go to 5
      ut = Ex(i) + (Ex(i+1)-Ex(i))/(rx(i+1)-rx(i))*(rr-rx(i)) -Easx
    5 continue
                    end if
      if(ind5.eq.2) then
      ut=ppvalu(rx,coefx,ir5,4,rr,0)
      ut1=ppvalu(rx,coefx,ir5,4,rr,1)
                end if
                        else 
      ut=Ex(ir5) -Easx
      ut1=0.0
                        end if
                    end if
      if(ind5.eq.1) ut=utn(rr/rauA,i0,0)
      return
                 end if
      return 
      end


      subroutine rgni(n,x,grad,ereal,gradt) !  ,h0,h1,ee,ev,w,natoms)
        implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
        INTEGER INFO, NEMAX
        DOUBLE PRECISION TEMPA(9*N)

      parameter (nemax=100)
        integer iw(nemax,3),  IWORK(3*5*N), IFAIL(3*N), NCALL, N, I, ID, IEV, IT, J, K, N3

        DOUBLE PRECISION x(3,n),h0(n*3,n*3)
     &        ,h1(n*3,n*3)
     &        ,grad(n*3),ee(n*3),ev(n*3,n*3),w(n*3,n*3),WORK(3*8*N)
      logical gradt
      data iev/1/, ncall/1/, shift/10.0d0/
     &  , Elow/-1.d3/,eps/1.d-6/,nstep/1000/,message/6/,nwork/777/
     &  , iflag/-2/ 
      save

      n3=n*3
cc      if(gradt) igrad=1

        id=1
        it=0
c       it=1
C       it=2
      call hmat(x,n,n3,h0,h1,h2,0)
      if(ncall.eq.1.or.it.eq.0) then
         if(id.eq.0) then
C            call jacob2(h0,ev,w(1,1),ee,w(1,2),n3,iev)
                     else
C       call tred2(h0,n3,n3,ee,w(1,1),iev)
C       call tqli(ee,w(1,1),n3,n3,h0,iev)

        CALL DSYEV('V','U',N3,H0,N3,EE,TEMPA,9*N,INFO)
        IF (INFO.NE.0) THEN
           PRINT*,'WARNING - INFO=',INFO,' in DSYEV'
        ENDIF

C       ABSTOL=0.0D0
C
C  This is loop is only needed for DSYEVX when only one eigenvalue is found.
C  The ev=h0 loop should be commented is DSYEVX is used.
C
C       DO J1=1,N3
C          EE(J1)=1.0D6
C       ENDDO
C       CALL DSYEVX('V','I','U',N3,H0,N3,DUMMY1,DUMMY2,1,1,ABSTOL,M,EE,EV,N3,WORK,3*8*N,IWORK,IFAIL,INFO)

        do i=1,n3
        do j=1,n3
        ev(i,j)=h0(i,j)
        enddo
        enddo
                     end if
        call sortv(ee,ev,n3,iev)
        ereal=ee(n3)

c      if(it.eq.2) open(nwork,file='itlan.tmp')

      ncall=2
                                else
      if(it.eq.1) then
      do i=1,n3
      h0(i,i)=h0(i,i)-shift
      enddo
C        call iteig(h0,ev(1,n3),ereal,w(1,1),w(1,2),n3)
      ereal=ereal+shift
                end if

      if(it.eq.2) then

        Ehigh=dmax1(0.d0,2.d0*ereal)
      do istep=1,nstep
c      call itlane(n3,Elow,Ehigh,eps,n,n3,nstep,message,nwork
c     &      ,iflag,ev,ev(1,n3),ee,iw,ne,w,w(1,2),iw(1,3),h1,h1(1,2))
      if(iflag.eq.0) go to 1
      if(iflag.eq.1) then
      do i=1,n3
      hvi=0.d0
      do j=1,n3
      hvi=hvi+h0(i,j)*ev(j,n3)
      enddo
      ev(i,1)=ev(i,1) +hvi
      enddo
                          end if
        if(iflag.eq.2) then
        print *,' itlane: iflag=', iflag
        stop
                       end if
      enddo
  1     ereal=ee(1)
      if(iev.ne.0) then
      nw=0
      PRINT '(A)','error, ne has not been set!'
      STOP
      do i=1,ne
         nw=max0(nw,iw(1,i))
      enddo
c        call itlanv(n3,nstep,message,nwork
c     &      ,ee,iw,1,w,w(1,2),n3,nw,jflag,ev(1,n3),h1,h1(1,2))
                 end if
C     if(jflag.gt.0) then
C        print *,' itlanv: jflag= ',jflag
C        stop
C     end if

                end if

                                end if

      if(gradt) then

      do i=1,n3
      call hmat(x,n,n3,h0,h1,h2,i)
      vhvi=0.0d0
      do j=1,n3
      vh=0.0d0
      do k=1,n3
cc      vh=vh+ev(k,n3)*h1(k,j,i)
      vh=vh+ev(k,n3)*h1(k,j)
      enddo
      vhvi=vhvi+vh*ev(j,n3)
      enddo
      grad(i)=vhvi
      enddo

              end if
      return
      end
c ------------------------------------------------------------
c      DOUBLE PRECISION function energy(r0,i,ip)
c      DOUBLE PRECISION r0
c      r=r0
c      energy=ut(r,i,ip)
c      return
c      end

              subroutine sortv(ee,ev,n,iev)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION ee(n),ev(n,n),em,evnn
      DO 17 M=1,N-1
      NN=N+1-M
      EM=EE(1)
      IND=1
      DO 15 K=2,NN
      IF(EE(K).GE.EM) GO TO 15
      EM=EE(K)
      IND=K
   15 CONTINUE
      IF(IND.EQ.NN) GO TO 17
      EE(IND)=EE(NN)
      EE(NN)=EM
      if(iev.gt.0) then
      do 16 k=1,n
      evnn=ev(k,nn)
      ev(k,nn)=ev(k,ind)
      ev(k,ind)=evnn
   16 continue
                   end if
   17 CONTINUE
      RETURN
      END

      subroutine hmat(x,n,n3,h0,h1,h2,igrad)
      USE commons
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
C       common/test/itest

      DOUBLE PRECISION x(3,n),h0(n*3,n*3)
     &               ,h1(n*3,n*3)

      save
      data RauA/0.52918d0/
cc     &      ,isu/1/,isg/2/,ipu/3/,ipg/4/,ix/5/

      if(igrad.ne.0) go to 10

      do i=1,n3
      do j=1,n3
      h0(i,j)=0.0d0
      enddo
      enddo

      ExS=0.0d0
      do k=1,n-1
      do l=k+1,n
      rkl=
     &      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
cc      Ex=energy(rkl,ix,0)
        IF (NEON) THEN
           Ex=utn(rkl/RauA,5,0)
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        ENDIF
      ExS=ExS+Ex
      enddo
      enddo

      do i=1,n3-1
      h0(i,i)=h0(i,i)+ExS
      do j=i+1,n3

      kk=int(i/3)
      if(mod(i,3).gt.0) kk=kk+1
      ll=int(j/3)
      if(mod(j,3).gt.0) ll=ll+1
      if(ll.eq.kk) go to 1

      dx=x(1,ll)-x(1,kk)
      dy=x(2,ll)-x(2,kk)
      dz=x(3,ll)-x(3,kk)
      dxy=sqrt(dx*dx+dy*dy)
      rkl=sqrt(dxy*dxy+dz*dz)
      ct=dz/rkl
      st=dxy/rkl
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if
        IF (NEON) THEN
           call fenergy(rkl,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rkl,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rkl,isu,0)
cc      Esg=energy(rkl,isg,0)
cc      Epu=energy(rkl,ipu,0)
cc      Epg=energy(rkl,ipg,0)
      Ess=0.5d0*(Esu+Esg)
      Esd=0.5d0*(Esu-Esg)
      Eps=0.5d0*(Epu+Epg)
      Epd=0.5d0*(Epu-Epg)
cc      Ex=energy(rkl,ix,0)
      
      if(mod(i,3).eq.1) then
      if(mod(j,3).eq.1) then
         trig1=(st*cf)**2
         trig2=(ct*cf)**2+sf**2
cc      hs=Ess*(st*cf)**2+Eps*((ct*cf)**2+sf**2)
      hs=Ess*trig1+Eps*trig2
      h0(i,i)=h0(i,i)+hs-Ex
      h0(j,j)=h0(j,j)+hs-Ex
cc      hd=Esd*(st*cf)**2+Epd*((ct*cf)**2+sf**2)      
      hd=Esd*trig1+Epd*trig2
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.2) then
         trig=st**2*sf*cf
cc        hs=(Ess-Eps)*st**2*sf*cf
        hs=(Ess-Eps)*trig
        h0(i,i+1)=h0(i,i+1)+hs
        h0(i+1,i)=h0(i+1,i)+hs
cc      hd=(Esd-Epd)*st**2*sf*cf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.0) then
         trig=st*ct*cf
cc        hs=(Ess-Eps)*st*ct*cf
        hs=(Ess-Eps)*trig
        h0(i,i+2)=h0(i,i+2)+hs
        h0(i+2,i)=h0(i+2,i)+hs
cc      hd=(Esd-Epd)*st*ct*cf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
                    end if

      if(mod(i,3).eq.2) then
      if(mod(j,3).eq.2) then
         trig1=(st*sf)**2
         trig2=(ct*sf)**2+cf**2
cc      hs=Ess*(st*sf)**2+Eps*((ct*sf)**2+cf**2)
      hs=Ess*trig1+Eps*trig2
      h0(i,i)=h0(i,i)+hs-Ex
      h0(j,j)=h0(j,j)+hs-Ex
cc      hd=Esd*(st*sf)**2+Epd*((ct*sf)**2+cf**2)
      hd=Esd*trig1+Epd*trig2
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.1) then
         trig=st**2*sf*cf
cc        hs=(Ess-Eps)*st**2*sf*cf
        hs=(Ess-Eps)*trig
        h0(j,j+1)=h0(j,j+1)+hs
        h0(j+1,j)=h0(j+1,j)+hs
cc      hd=(Esd-Epd)*st**2*sf*cf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.0) then
         trig=st*ct*sf
cc        hs=(Ess-Eps)*st*ct*sf
        hs=(Ess-Eps)*trig
        h0(i,i+1)=h0(i,i+1)+hs
        h0(i+1,i)=h0(i+1,i)+hs
cc      hd=(Esd-Epd)*st*ct*sf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
                    end if

      if(mod(i,3).eq.0) then
      if(mod(j,3).eq.0) then
         trig1=(ct)**2
         trig2=(st)**2
cc      hs=Ess*(ct)**2   +Eps*(st)**2
      hs=Ess*trig1   +Eps*trig2
      h0(i,i)=h0(i,i)+hs-Ex
      h0(j,j)=h0(j,j)+hs-Ex
cc      hd=Esd*(ct)**2   +Epd*(st)**2
      hd=Esd*trig1   +Epd*trig2
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.1) then
         trig=st*ct*cf
cc        hs=(Ess-Eps)*st*ct*cf
        hs=(Ess-Eps)*trig
        h0(j,j+2)=h0(j,j+2)+hs
        h0(j+2,j)=h0(j+2,j)+hs
cc      hd=(Esd-Epd)*st*ct*cf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
      if(mod(j,3).eq.2) then
         trig=st*ct*sf
cc        hs=(Ess-Eps)*st*ct*sf
        hs=(Ess-Eps)*trig
        h0(j,j+1)=h0(j,j+1)+hs
        h0(j+1,j)=h0(j+1,j)+hs
cc      hd=(Esd-Epd)*st*ct*sf
      hd=(Esd-Epd)*trig
      h0(i,j)=hd
      h0(j,i)=hd
                    end if
                    end if

   1      continue
      enddo
      enddo
      h0(n3,n3)=h0(n3,n3)+ExS

      IF(DIPOLE) CALL dipoles(x,n,h0,h1,igrad)

c ---------------------------------------------------------------------

  10      continue
      if(igrad.le.0) go to 2

      do i=1,n3
      do j=1,n3
cc      do k=1,n3
cc      h1(i,j,k)=0.0d0
cc      enddo
      h1(i,j)=0.0d0
      enddo
      enddo

cc      do m=1,n
cc      do ir=1,3
cc      ic=(m-1)*3+ir

      m=int(igrad/3)+1
      ir=mod(igrad,3)
      if(ir.eq.0) then
      ir=3
      m=m-1
                  end if

      ExS1=0.0d0
      do k=1,n-1
      do l=k+1,n
      if(l.eq.m.or.k.eq.m) then
      rkl=
     &      sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
cc      if(l.eq.m) Ex1= energy(rkl,ix,1)*(x(ir,l)-x(ir,k))/rkl
cc      if(k.eq.m) Ex1=-energy(rkl,ix,1)*(x(ir,l)-x(ir,k))/rkl
        IF (NEON) THEN
      if(l.eq.m) Ex1= utn(rkl/RauA,5,1)/RauA*(x(ir,l)-x(ir,k))/rkl
      if(k.eq.m) Ex1=-utn(rkl/RauA,5,1)/RauA*(x(ir,l)-x(ir,k))/rkl
        ELSE
           CALL ASAR1(rkl,Ex,Ex1)
        if(l.eq.m) Ex1= Ex1*(x(ir,l)-x(ir,k))/rkl
        if(k.eq.m) Ex1=-Ex1*(x(ir,l)-x(ir,k))/rkl
        ENDIF
      ExS1=ExS1+Ex1
                       end if
      enddo
      enddo

      do i=1,n3-1
      h1(i,i)=h1(i,i)+ExS1
      do j=i+1,n3

      kk=int(i/3)
      if(mod(i,3).gt.0) kk=kk+1
      ll=int(j/3)
      if(mod(j,3).gt.0) ll=ll+1
      if(ll.eq.kk) go to 11

      if(ll.eq.m.or.kk.eq.m) then

      dx=x(1,ll)-x(1,kk)
      dy=x(2,ll)-x(2,kk)
      dz=x(3,ll)-x(3,kk)
      dxy=sqrt(dx*dx+dy*dy)
      rkl=sqrt(dxy*dxy+dz*dz)
      ct=dz/rkl
      st=dxy/rkl
      cf=1.0d0
      sf=0.0d0
      if(dxy.gt.0.0d0) then
      cf=dx/dxy
      sf=dy/dxy
                   end if

      if(ll.eq.m) drdx= (x(ir,ll)-x(ir,kk))/rkl
      if(kk.eq.m) drdx=-(x(ir,ll)-x(ir,kk))/rkl

      ct1=-ct/rkl*drdx
      if(ir.eq.3) then
      if(m.eq.ll) ct1=ct1+1.0d0/rkl
      if(m.eq.kk) ct1=ct1-1.0d0/rkl
                end if
      st1=-st/rkl*drdx
      dxydx=0.0d0
      if(ir.eq.1.or.ir.eq.2) then
c      dxydx=1.0d0
      if(dxy.gt.0.0d0) then
      if(ll.eq.m) dxydx= (x(ir,ll)-x(ir,kk))/dxy
      if(kk.eq.m) dxydx=-(x(ir,ll)-x(ir,kk))/dxy
                   end if
      st1=st1+dxydx/rkl
                             end if
      cf1=0.0d0
      sf1=0.0d0
      if(dxy.gt.0.0d0) then
      if(ir.eq.1.or.ir.eq.2) then
      cf1=-cf/dxy*dxydx
      if(ir.eq.1) then
      if(m.eq.ll) cf1=cf1+1.0d0/dxy
      if(m.eq.kk) cf1=cf1-1.0d0/dxy
                    end if
      sf1=-sf/dxy*dxydx
      if(ir.eq.2) then
      if(m.eq.ll) sf1=sf1+1.0d0/dxy
      if(m.eq.kk) sf1=sf1-1.0d0/dxy
                    end if
                             end if
                   end if

        IF (NEON) THEN
           call fenergy(rkl,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ELSE
         call energyar(rkl,Esu,Esg,Epu,Epg,Ex,Esu1,Esg1,Epu1,Epg1,Ex1)
        ENDIF
cc      Esu=energy(rkl,isu,0)
cc      Esg=energy(rkl,isg,0)
cc      Epu=energy(rkl,ipu,0)
cc      Epg=energy(rkl,ipg,0)
      Ess=0.5d0*(Esu+Esg)
      Esd=0.5d0*(Esu-Esg)
      Eps=0.5d0*(Epu+Epg)
      Epd=0.5d0*(Epu-Epg)
cc      Ex=energy(rkl,ix,0)

cc      Esu1=energy(rkl,isu,1)*drdx
cc      Esg1=energy(rkl,isg,1)*drdx
cc      Epu1=energy(rkl,ipu,1)*drdx
cc      Epg1=energy(rkl,ipg,1)*drdx
      Esu1=Esu1*drdx
      Esg1=Esg1*drdx
      Epu1=Epu1*drdx
      Epg1=Epg1*drdx
      Ess1=0.5d0*(Esu1+Esg1)
      Esd1=0.5d0*(Esu1-Esg1)
      Eps1=0.5d0*(Epu1+Epg1)
      Epd1=0.5d0*(Epu1-Epg1)
cc      Ex1=energy(rkl,ix,1)*drdx
      Ex1=Ex1*drdx
      
      if(mod(i,3).eq.1) then
      if(mod(j,3).eq.1) then
         trig1=(st*cf)**2
         trig2=(ct*cf)**2+sf**2
         trig3=2.0d0*st*cf*(cf*st1+st*cf1)
         trig4=(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
c      hs=Ess*(st*cf)**2+Eps*((ct*cf)**2+sf**2)
cc      hs1=Ess1*(st*cf)**2+Eps1*((ct*cf)**2+sf**2) +
cc     &          Ess*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     &      Eps*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
      hs1= Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(i,i)=h1(i,i)+hs1-Ex1
      h1(j,j)=h1(j,j)+hs1-Ex1
c      hd=Esd*(st*cf)**2+Epd*((ct*cf)**2+sf**2)      
cc      hd1=Esd1*(st*cf)**2+Epd1*((ct*cf)**2+sf**2) +
cc     &          Esd*2.0d0*st*cf*(cf*st1+st*cf1)+
cc     &      Epd*(2.0d0*ct*cf*(cf*ct1+ct*cf1)+2.0d0*sf*sf1)
      hd1= Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.2) then
         trig1=st**2*sf*cf
         trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
c        hs=(Ess-Eps)*st**2*sf*cf
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     &      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)) 
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(i,i+1)=h1(i,i+1)+hs1
        h1(i+1,i)=h1(i+1,i)+hs1
c      hd=(Esd-Epd)*st**2*sf*cf
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     &      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.0) then
         trig1=st*ct*cf
         trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
c        hs=(Ess-Eps)*st*ct*cf
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     &      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(i,i+2)=h1(i,i+2)+hs1
        h1(i+2,i)=h1(i+2,i)+hs1
c      hd=(Esd-Epd)*st*ct*cf
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     &      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
                    end if

      if(mod(i,3).eq.2) then
      if(mod(j,3).eq.2) then
         trig1=(st*sf)**2
         trig2=(ct*sf)**2+cf**2
         trig3=2.0d0*st*sf*(sf*st1+st*sf1)
         trig4=2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1
c      hs=Ess*(st*sf)**2+Eps*((ct*sf)**2+cf**2)
cc      hs1=Ess1*(st*sf)**2+Eps1*((ct*sf)**2+cf**2) +
cc     &          Ess*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     &      Eps*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
      hs1=Ess1*trig1+Eps1*trig2+Ess*trig3+Eps*trig4
      h1(i,i)=h1(i,i)+hs1-Ex1
      h1(j,j)=h1(j,j)+hs1-Ex1
c      hd=Esd*(st*sf)**2+Epd*((ct*sf)**2+cf**2)
cc      hd1=Esd1*(st*sf)**2+Epd1*((ct*sf)**2+cf**2) +
cc     &          Esd*2.0d0*st*sf*(sf*st1+st*sf1)+
cc     &      Epd*(2.0d0*ct*sf*(sf*ct1+ct*sf1)+2.0d0*cf*cf1)
      hd1=Esd1*trig1+Epd1*trig2+Esd*trig3+Epd*trig4
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.1) then
         trig1=st**2*sf*cf
         trig2=2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1)
c        hs=(Ess-Eps)*st**2*sf*cf
cc        hs1=(Ess1-Eps1)*st**2*sf*cf +
cc     &      (Ess-Eps)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(j,j+1)=h1(j,j+1)+hs1
        h1(j+1,j)=h1(j+1,j)+hs1
c      hd=(Esd-Epd)*st**2*sf*cf
cc      hd1=(Esd1-Epd1)*st**2*sf*cf + 
cc     &      (Esd-Epd)*(2.0d0*st*st1*sf*cf+st**2*(sf1*cf+sf*cf1))
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.0) then
         trig1=st*ct*sf
         trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
c        hs=(Ess-Eps)*st*ct*sf
cc        hs1=(Ess1-Eps1)*st*ct*sf +
cc     &      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(i,i+1)=h1(i,i+1)+hs1
        h1(i+1,i)=h1(i+1,i)+hs1
c      hd=(Esd-Epd)*st*ct*sf
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     &      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
                    end if

      if(mod(i,3).eq.0) then
      if(mod(j,3).eq.0) then
         trig1=(ct)**2
         trig2=(st)**2
         trig3=2.0d0*ct*ct1
         trig4=2.0d0*st*st1
c      hs=Ess*(ct)**2   +Eps*(st)**2
cc      hs1=Ess1*(ct)**2   +Eps1*(st)**2 +
cc     &      Ess*2.0d0*ct*ct1   +Eps*2.0d0*st*st1
      hs1=Ess1*trig1 + Eps1*trig2 + Ess*trig3 + Eps*trig4
      h1(i,i)=h1(i,i)+hs1-Ex1
      h1(j,j)=h1(j,j)+hs1-Ex1
c      hd=Esd*(ct)**2   +Epd*(st)**2
cc      hd1=Esd1*(ct)**2   +Epd1*(st)**2 +
cc     &      Esd*2.0d0*ct*ct1   +Epd*2.0d0*st*st1
      hd1=Esd1*trig1 + Epd1*trig2 + Esd*trig3 + Epd*trig4
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.1) then
         trig1=st*ct*cf
         trig2=st1*ct*cf+st*ct1*cf+st*ct*cf1
c        hs=(Ess-Eps)*st*ct*cf
cc        hs1=(Ess1-Eps1)*st*ct*cf + 
cc     &      (Ess-Eps)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(j,j+2)=h1(j,j+2)+hs1
        h1(j+2,j)=h1(j+2,j)+hs1
c      hd=(Esd-Epd)*st*ct*cf
cc      hd1=(Esd1-Epd1)*st*ct*cf + 
cc     &      (Esd-Epd)*(st1*ct*cf+st*ct1*cf+st*ct*cf1)
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
      if(mod(j,3).eq.2) then
         trig1=st*ct*sf
         trig2=st1*ct*sf+st*ct1*sf+st*ct*sf1
c        hs=(Ess-Eps)*st*ct*sf
cc        hs1=(Ess1-Eps1)*st*ct*sf + 
cc     &      (Ess-Eps)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
        hs1=(Ess1-Eps1)*trig1 + (Ess-Eps)*trig2
        h1(j,j+1)=h1(j,j+1)+hs1
        h1(j+1,j)=h1(j+1,j)+hs1
c      hd=(Esd-Epd)*st*ct*sf
cc      hd1=(Esd1-Epd1)*st*ct*sf + 
cc     &      (Esd-Epd)*(st1*ct*sf+st*ct1*sf+st*ct*sf1)
      hd1=(Esd1-Epd1)*trig1 + (Esd-Epd)*trig2
      h1(i,j)=hd1
      h1(j,i)=hd1
                    end if
                    end if

                       end if

  11      continue
      enddo
      enddo
      h1(n3,n3)=h1(n3,n3)+ExS1

cc      enddo
cc      enddo

        IF(DIPOLE) CALL dipoles(x,n,h0,h1,igrad)

   2      continue

      return
      end

cc           real*4 function utn(r,i,ip)
cc           utn=0.0d0
cc           if(i.ne.5) return
cc           utn=asne1(r,i,ip)
cc           return
cc           end
cc                 FUNCTION ASNE1(R,I,IP)
                 DOUBLE PRECISION FUNCTION UTN(R,I,IP)
           implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
           save
             COMMON/DIPOL/ adip,Runit,ReX                                            
           DATA E/1.34d-4/,RM/5.841d0/, A/13.86434671d0/,B/-.12993822d0/, 
     & D/1.36d0/,
     * C6,C8,C10/1.21317545d0,.53222749d0,.24570703d0/, AA/8.9571795d5/
      ReX=RM
           X=R/RM
           F=1.d0
           IF(X.LT.D) F=EXP(-(D/X-1.d0)**2)

           X2=X*X
           X4=X2*X2
           X6=X4*X2
           X8=X4*X4
           X10=X6*X4
           if(ip.eq.0) then
           V=AA*EXP((-A+B*X)*X)-F*(C6/X6+C8/X8+C10/X10)
cc           ASNE1=V*E
           UTN=V*E
           RETURN
                       end if
           if(ip.eq.1) then
           F1=0.d0
           IF(X.LT.D) F1=F*2.d0*(D/X-1.d0)*D/X2
           V=AA*EXP((-A+B*X)*X) *(-A+2.d0*B*X)
     * +F*(6.d0*C6/X6+8.d0*C8/X8+10.d0*C10/X10)/X
           IF(X.LT.D) V=V-F1*(C6/X6+C8/X8+C10/X10)
cc           ASNE1=V*E/RM
           UTN=V*E/RM
                       end if
           RETURN
           END

        subroutine rgnx2(n,x,grad,ereal,gradt)
        implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
        DOUBLE PRECISION x(3,n),grad(n*3)
        logical gradt
        save

      data hrx2/0.995D0/, EaueV/27.212D0/, RauA/0.52918D0/

      Ereal=0.d0
        if(n.le.0) return
      do i=1,n*3
      grad(i)=0.d0
      enddo

        if(n.gt.1) then
        do k=1,n-1
        do l=k+1,n
        rkl=
     1  sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
        call ashe1(rkl/RauA,Err,Err1)
        Ereal=Ereal+Err*EaueV

      if(gradt) then
      do m=1,3
      dEdxm=Err1*EaueV/RauA*(x(m,l)-x(m,k))/rkl
      lgrad=3*(l-1)+m
      grad(lgrad)=grad(lgrad)+dEdxm
        kgrad=3*(k-1)+m
        grad(kgrad)=grad(kgrad)-dEdxm
      enddo
              end if
        enddo
        enddo
                   end if

      do i=1,n
      d2xy=x(1,i)**2+x(2,i)**2
      Rxa=sqrt(d2xy + (x(3,i)-hrx2)**2)
        Rxb=sqrt(d2xy + (x(3,i)+hrx2)**2)                
      cxa=(x(3,i)-hrx2)/Rxa
      cxb=(x(3,i)+hrx2)/Rxb
      call ehecl2(Rxa,cxa,Erxa,Erxa1,Erxa1c)
        call ehecl2(Rxb,-cxb,Erxb,Erxb1,Erxb1c)
      Ereal=Ereal+Erxa+Erxb

      if(gradt) then
      do m=1,3
      if(m.le.2) then
      dEadxm=(Erxa1 - Erxa1c*cxa/Rxa)*x(m,i)/Rxa
        dEbdxm=(Erxb1 + Erxb1c*cxb/Rxb)*x(m,i)/Rxb
               else
        dEadxm=(Erxa1 - Erxa1c*cxa/Rxa)*(x(m,i)-hrx2)/Rxa + Erxa1c/Rxa
        dEbdxm=(Erxb1 + Erxb1c*cxb/Rxb)*(x(m,i)+hrx2)/Rxb - Erxb1c/Rxb
               end if
      igrad=3*(i-1)+m
        grad(igrad)=grad(igrad)+dEadxm+dEbdxm
      enddo
              end if
      enddo

      return
      end


      subroutine ehecl2(r,c,e,e1,e1c)
        implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      common/init/Es,Ep,Epb
      data hre/0.995D0/, nc/1/
      save

      Es = hecl2ut0(r,1,Es1)
      Ep = hecl2ut0(r,2,Ep1)
        Epb= hecl2ut0(r,3,Epb1)

        ratio=hre/r
        A= 1.d0 + ratio*ratio
        B= 2.*ratio
        A1= -2.d0*(A-1.d0)/r
        B1=-B/r

        c2= c*c
      s2= 1.d0 - c2
      ABc= A + B*c
      s02= s2/ABc
      s021= -s02/ABc*(A1 + B1*c)
      s021c= -(2.d0*c + s02*B)/ABc
      Ept=Epb
      Ept1=Epb1
      if(s02.gt.0.5d0) then
      factor= 4.d0*s02*(1.d0 - s02)
        factor1= 4.d0*s021*(1.d0 - 2.d0*s02)
      factor1c= 4.d0*s021c*(1.d0 - 2.d0*s02)
      Ept= Ep + (Epb - Ep)*factor
      Ept1= Ep1 + (Epb1 - Ep1)*factor + (Epb - Ep)*factor1
        Ept1c= (Epb - Ep)*factor1c 
                   end if

      e = Es*c2 + Ept*s2 
      e1= Es1*c2 + Ept1*s2
      e1c= 2.d0*c*(Es - Ept) + Ept1c*s2

      return
      end

      DOUBLE PRECISION function hecl2ut0(rr,i0,ut1)
        implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER IR
      parameter(ir=41)
      DOUBLE PRECISION, DIMENSION(41) :: ru=(/
     1      2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,     
     1      3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,
     1      4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90,
     1      5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90,
     1      6.00 /)
      DOUBLE PRECISION, DIMENSION(41) :: Eu=(/
     1       0.3609288,  0.2394266,  0.1526561,  0.0937677,  0.0559396,
     1        0.0322421,  0.0168441,  0.0070019,  0.0010727, -0.0022567,
     1      -0.0039128, -0.0046519, -0.0048347, -0.0046815, -0.0043490,
     1      -0.0039317, -0.0035008, -0.0030874, -0.0027041, -0.0023572,
     1      -0.0020523, -0.0017784, -0.0015446, -0.0013409, -0.0011662,
     1       -0.0010178, -0.0008888, -0.0007792, -0.0006841, -0.0006028,
     1      -0.0005304, -0.0004683, -0.0004153, -0.0003680, -0.0003261,
     1      -0.0002896, -0.0002588, -0.0002309, -0.0002080, -0.0001890,
     1      -0.0001690/)
      DOUBLE PRECISION, DIMENSION(41) :: Eg=(/
     1       0.6226088,  0.4446164,  0.3092623,  0.2090408,  0.1372434,
     1       0.0879939,  0.0557593,  0.0348669,  0.0208967,  0.0116730,
     1       0.0057593,  0.0021143, -0.0000627, -0.0013357, -0.0020061,
     1      -0.0023021, -0.0023595, -0.0022811, -0.0021323, -0.0019467,
     1      -0.0017533, -0.0015611, -0.0013833, -0.0012201, -0.0010725,
     1      -0.0009412, -0.0008269, -0.0007237, -0.0006362, -0.0005610,
     1      -0.0004941, -0.0004352, -0.0003848, -0.0003403, -0.0003010,
     1      -0.0002675, -0.0002369, -0.0002114, -0.0001892, -0.0001696,
     1      -0.0001528/)
      DOUBLE PRECISION, DIMENSION(41) :: Epu=(/
     1       0.8823669,  0.6168400,  0.4294591,  0.2970549,  0.2044612,
     1       0.1390856,  0.0935656,  0.0623482,  0.0408527,  0.0261192,
     1       0.0161569,  0.0095486,  0.0051923,  0.0023638,  0.0005976,
     1      -0.0004888, -0.0011015, -0.0014176, -0.0015430, -0.0015577,
     1      -0.0015120, -0.0014252, -0.0013083, -0.0011809, -0.0010553,
     1      -0.0009371, -0.0008288, -0.0007346, -0.0006490, -0.0005784,
     1      -0.0005115, -0.0004544, -0.0004049, -0.0003628, -0.0003226,
     1      -0.0002862, -0.0002551, -0.0002274, -0.0001998, -0.0001736,
     1      -0.0001570/)
      DOUBLE PRECISION coefu(4,ir), coefg(4,ir), coefpu(4,ir)
      data nc1/1/,nc2/1/,nc3/1/
      save

      hecl2ut0=0.d0
      ut1=0.d0

      if(i0.eq.1) then
      if(nc1.eq.1) then
      nc1=2
      do i=1,ir
      coefu(1,i)=Eu(i) 
      enddo
      coefu(2,ir)=-6.d0*Eu(ir)/ru(ir)
      call cubspl(ru,coefu,ir,0,1)
        r0u=ru(1)
        AAu=Eu(1)
        alphau=-coefu(2,1)/AAu
                   end if
        if(rr.lt.ru(ir)) then
        if(rr.lt.r0u) then
        hecl2ut0=AAu*exp(-alphau*(rr-r0u))
      ut1=-alphau*hecl2ut0
                        else
      ind=int(10.d0*(rr-r0u)) +1
      dr=rr-ru(ind)
      hecl2ut0=coefu(1,ind)+dr*(coefu(2,ind)
     1      +0.5d0*dr*(coefu(3,ind)+dr*coefu(4,ind)/3.d0))
        ut1=coefu(2,ind) + dr*(coefu(3,ind) + 0.5d0*dr*coefu(4,ind))
                        end if
                         else
            hecl2ut0=Eu(ir)*(ru(ir)/rr)**6
      ut1=-6.d0*hecl2ut0/rr
                         end if
      return
                  end if

      if(i0.eq.2) then
      if(nc2.eq.1) then
      nc2=2
      do i=1,ir
      coefg(1,i)=Eg(i)
      enddo
      coefg(2,ir)=-6.d0*Eg(ir)/ru(ir)
      call cubspl(ru,coefg,ir,0,1)
        r0g=ru(1)
        AAg=Eg(1)
        alphag=-coefg(2,1)/AAg
                   end if
        if(rr.lt.ru(ir)) then
        if(rr.lt.r0g) then
        hecl2ut0=AAg*exp(-alphag*(rr-r0g))
        ut1=-alphag*hecl2ut0
                        else
        ind=int(10.d0*(rr-r0g)) +1
        dr=rr-ru(ind)
        hecl2ut0=coefg(1,ind)+dr*(coefg(2,ind)
     1      +0.5d0*dr*(coefg(3,ind)+dr*coefg(4,ind)/3.d0))        
        ut1=coefg(2,ind) + dr*(coefg(3,ind) + 0.5d0*dr*coefg(4,ind))
                        end if
                         else
        hecl2ut0=Eg(ir)*(ru(ir)/rr)**6
        ut1=-6.d0*hecl2ut0/rr
                         end if
      return
                  end if

      if(i0.eq.3) then
      if(nc3.eq.1) then
      nc3=2
      do i=1,ir
      coefpu(1,i)=Epu(i)
      enddo
      coefpu(2,ir)=-6.d0*Epu(ir)/ru(ir)
      call cubspl(ru,coefpu,ir,0,1)
        r0pu=ru(1)
        AApu=Epu(1)
        alphapu=-coefpu(2,1)/AApu
                   end if
       if(rr.lt.ru(ir)) then
       if(rr.lt.r0pu) then
        hecl2ut0=AApu*exp(-alphapu*(rr-r0pu))
        ut1=-alphapu*hecl2ut0
                        else
        ind=int(10.d0*(rr-r0pu)) +1
        dr=rr-ru(ind)
        hecl2ut0=coefpu(1,ind)+dr*(coefpu(2,ind)
     1      +0.5d0*dr*(coefpu(3,ind)+dr*coefpu(4,ind)/3.d0))
        ut1=coefpu(2,ind) + dr*(coefpu(3,ind) + 0.5d0*dr*coefpu(4,ind))
                        end if
                         else
        hecl2ut0=Epu(ir)*(ru(ir)/rr)**6
        ut1=-6.d0*hecl2ut0/rr                               
                         end if
      return
                  end if
      return 
      end

           SUBROUTINE ASHE1(R,Ex,Ex1)
C  Potential of He2 (Aziz et al, 1987)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER NCALL
      save
      DATA E/10.948D0/,RM/2.963D0/, A/10.43329537D0/,B/-2.27965105D0/, D/1.4826D0/
     *    ,C6,C8,C10/1.36745214D0,0.42123807D0,0.17473318D0/, AA/1.8443101D5/
     *    ,NCALL/1/,Ekau/3.1668D-6/,Rau/0.52918D0/
      IF(NCALL.EQ.1) THEN
      RM=RM/Rau
      E=E*Ekau
      NCALL=2
                     END IF
      X=R/RM
        X2=X*X
      F=1.d0
      F1=0.d0
      IF(X.LT.D) then
      dx1=D/X-1.d0
      F=EXP(-dx1*dx1)
      F1=F*2.d0*dx1*D/X2
             end if
      X6=X2*X2*X2
      VdW=(C6+(C8+C10/X2)/X2)/X6
      Aexabx=AA*EXP((-A+B*X)*X)
      V=Aexabx - F*VdW
      V1=(-A+2.d0*B*X)*Aexabx - F1*VdW + 
     1            F*(6.d0*C6+(8.d0*C8+10.d0*C10/X2)/X2)/X6/X
      Ex=V*E
      Ex1=V1*E/RM
      RETURN
      END
        subroutine rgnxy(n,x,grad,ereal,gradt,socouple)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
        save
        DOUBLE PRECISION x(3,n),grad(n*3)
        logical gradt,socouple

      data EaueV/27.212D0/, RauA/0.52918D0/
     1      hsplit/0.0074267D0/

      Ereal=0.d0
        if(n.le.0) return
      do i=1,n*3
      grad(i)=0.d0
      enddo

        if(n.gt.1) then
        do k=1,n-1
        do l=k+1,n
        rkl=
     1  sqrt((x(1,l)-x(1,k))**2+(x(2,l)-x(2,k))**2+(x(3,l)-x(3,k))**2)
        call asar1(rkl/RauA,Err,Err1)
        Ereal=Ereal+Err*EaueV

      if(gradt) then
      do m=1,3
      dEdxm=Err1*EaueV/RauA*(x(m,l)-x(m,k))/rkl
      lgrad=3*(l-1)+m
      grad(lgrad)=grad(lgrad)+dEdxm
        kgrad=3*(k-1)+m
        grad(kgrad)=grad(kgrad)-dEdxm
      enddo
              end if
        enddo
        enddo
                   end if

      E1r=0.d0
      E2r=0.d0
      E12r=0.d0
      do i=1,n

      x1i=x(1,i)
      x2i=x(2,i)
      x3i=x(3,i)
      dxy2=x1i*x1i+x2i*x2i
      ri=sqrt(dxy2+x3i*x3i)
      ci=x3i/ri
      call arnoenergy(ri,ci,E1ri,E1ri1,E1ri1c,E2ri,E2ri1,E2ri1c)
      dxy=sqrt(dxy2)
      cai=1.d0
      sai=0.d0
      if(dxy.gt.0.d0) then
      cai=x1i/dxy
      sai=x2i/dxy
                  end if
      dEri=E1ri-E2ri
      cai2=cai*cai
      E1r=E1r +dEri*cai2 +E2ri
      E2r=E2r -dEri*cai2 +E1ri
      E12r=E12r +cai*sai*dEri

        enddo

        E12r2=E12r*E12r
      if(socouple) E12r2=E12r2 + hsplit*hsplit

      D = sqrt((E1r-E2r)**2 + 4.d0*E12r2)
      Ereal= Ereal + 0.5d0*(E1r+E2r + D)
      if(socouple) Ereal= Ereal - hsplit

      if(gradt) then

        do i=1,n

      x1i=x(1,i)
      x2i=x(2,i)
      x3i=x(3,i)
      dxy2=x1i*x1i+x2i*x2i
      ri=sqrt(dxy2+x3i*x3i)
      ci=x3i/ri
        call arnoenergy(ri,ci,E1ri,E1ri1,E1ri1c,E2ri,E2ri1,E2ri1c)
      dxy=sqrt(dxy2)
      cai=1.d0
      sai=0.d0
        if(dxy.gt.0.d0) then
      cai=x1i/dxy
      sai=x2i/dxy
                  end if

      dEri1=E1ri1-E2ri1
      cai2=cai*cai
      E1r1= dEri1*cai2 + E2ri1
      E2r1=-dEri1*cai2 + E1ri1
      E12r1= cai*sai*dEri1
      dEdri=0.5d0*( E1r1 + E2r1 + 
     1            ((E1r-E2r)*(E1r1-E2r1) + 4.d0*E12r*E12r1)/D )

      dEri1c=E1ri1c-E2ri1c
      E1r1c= dEri1c*cai2 + E2ri1c
      E2r1c=-dEri1c*cai2 + E1ri1c
      E12r1c= cai*sai*dEri1c
        dEdci=0.5d0*( E1r1c + E2r1c +
     1          ((E1r-E2r)*(E1r1c-E2r1c) + 4.d0*E12r*E12r1c)/D )

        dEri=E1ri-E2ri
      E1r1ca= 2.d0*cai*dEri
cc      E2r1ca= -E1r1ca
      E12r1ca= 0.d0
      if(sai.ne.0.d0) E12r1ca= (sai - cai2/sai)*dEri
        dEdcai= ((E1r-E2r)*E1r1ca + 2.d0*E12r*E12r1ca)/D 

      do m=1,3
      dridxm= x(m,i)/ri
      dcidxm= -ci/ri*dridxm
      if(m.eq.3) dcidxm= dcidxm +1.d0/ri
      dcaidxm= 0.d0
      if(m.le.2.and.dxy.gt.0.d0) then
      dcaidxm= -cai*x(m,i)/dxy2
      if(m.eq.1) dcaidxm= dcaidxm +1.d0/dxy
               end if
        dEdxm= dEdri*dridxm + dEdci*dcidxm + dEdcai*dcaidxm
      igrad=3*(i-1)+m
        grad(igrad)=grad(igrad)+dEdxm
      enddo

        enddo

              end if


      return
      end

      subroutine arnoenergy(r,c1,V,Vd,Vdc,V1,V1d,V1dc)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER N
      parameter (n=9)
      save
      DOUBLE PRECISION P(n),Vi(n),Vid(n),P1(n)

      if(abs(c1).gt.1.d0) c1=c1/abs(c1)

      V0 = ut0(r,1,V0d)
      n1=n-1
      do i=1,n1
        Vi(i) = ut0(r,i+1,Vdi)
      Vid(i)=Vdi
      enddo

        call polegn(n1,c1,P,P1)

        V = V0 
      Vd = V0d
      Vdc = 0.d0
      do i=1,n1
      V = V + Vi(i)*P(i)
      Vd = Vd + Vid(i)*P(i)
      Vdc = Vdc + Vi(i)*P1(i)
      enddo

        V0 = ut01(r,1,V0d)
      do i=1,n1
      Vi(i) = ut01(r,i+1,Vdi)
      Vid(i) = Vdi
      enddo

      V1 = V0
      V1d = V0d
      V1dc = 0.d0
      do i=1,n1
      V1 = V1 + Vi(i)*P(i)
      V1d = V1d + Vid(i)*P(i)
      V1dc = V1dc + Vi(i)*P1(i)
      enddo

      return
      end


      DOUBLE PRECISION function ut0(rr,i0,ut0d)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER N0, N1, N2, N3, N4, N5, N6, N7, N8
      save
      parameter(n0=24,n1=24,n2=24,n3=20,n4=20,n5=20,n6=20,n7=17,n8=18)

      DOUBLE PRECISION, DIMENSION(n0) :: ru=(/ 
     1 0.00, 0.25, 0.50, 0.75, 1.00,
     1 1.25, 1.50, 1.75, 2.00, 2.25,
     1 2.50, 2.75,
     1 3.00, 3.25, 3.50, 3.75, 4.00, 
     1 4.25, 4.50, 5.00, 5.50, 6.00, 
     1 7.00, 8.00/)

      DOUBLE PRECISION, DIMENSION(n0) :: Eu=(/
     123.5757810,19.0536413,15.1520331,11.8249968, 9.0265728,
     1 6.7107999, 4.8317169, 3.3433641, 2.1997803, 1.3550050,
     1  .7630777,  .3780378,
     1  .1539246,  .0447773,  .0045563, -.0075913, -.0095564,
     1 -.0083419, -.0064951, -.0035558, -.0019332, -.0010942,
     1 -.0004037, -.0001733/)
      DOUBLE PRECISION, DIMENSION(n1) :: Eg=(/
     1 2.9072499, 2.4029129, 1.9608503, 1.5769362, 1.2470459,
     1  .9670573,  .7328453,  .5402856,  .3852534,  .2636248,
     1  .1712753,  .1040805,
     1  .0579161,  .0286581,  .0122599,  .0046522,  .0013999, 
     1  .0001305, -.0002864, -.0003208, -.0001932, -.0001146,
     1 -.0000367, -.0000163/)
      DOUBLE PRECISION, DIMENSION(n2) :: Epu=(/
     136.1544630,29.2137453,23.2280250,18.1262749,13.8374705,
     110.2905828, 7.4145843, 5.1384474, 3.3911439, 2.1016460,
     1 1.1989267,  .6119579,
     1  .2697118,  .1011612,  .0353422,  .0104119,  .0016621,
     1 -.0009844, -.0014811, -.0010339, -.0005530, -.0002932,
     1 -.0000931, -.0000396/)
      DOUBLE PRECISION, DIMENSION(n3) :: E3=(/  
     1 5.7704509, 4.6198981, 3.6335576, 2.7987752, 2.1029035,
     1 1.5332943, 1.0772957,  .7222591,  .4555346,  .2644726,
     1  .1364226,  .0587354,
     1  .0187612,  .0038494,  .0012187,  .0002678, -.0000507,
     1 -.0001010, -.0000828, -.0000240/)
      DOUBLE PRECISION, DIMENSION(n4) :: E4=(/  
     1 6.4428066, 5.1725998, 4.0815667, 3.1560647, 2.3824570,
     1 1.7471023, 1.2363623,  .8365934,  .5341561,  .3154102,
     1  .1667152,  .0744307,
     1  .0249163,  .0045315, -.0003323, -.0010465, -.0007592,
     1 -.0004024, -.0001781, -.0000020/)
      DOUBLE PRECISION, DIMENSION(n5) :: E5=(/  
     1-4.6869956,-3.7386937,-2.9276822,-2.2432417,-1.6746429,
     1-1.2111622, -.8420752, -.5566573, -.3441831, -.1939286,
     1 -.0951689, -.0371794,
     1 -.0092352, -.0006112, -.0004883, -.0002925, -.0001497,
     1 -.0000855, -.0000426, -.0000077/)
      DOUBLE PRECISION, DIMENSION(n6) :: E6=(/ 
     1-4.1197525,-3.2840736,-2.5697842,-1.9673653,-1.4672967,
     1-1.0600520, -.7361093, -.4859449, -.3000351, -.1688565,
     1 -.0828862, -.0326009,
     1 -.0084771, -.0009921, -.0007325, -.0004486, -.0002404,
     1 -.0001257, -.0000563, -.0000099/)
      DOUBLE PRECISION, DIMENSION(n7) :: E7=(/  
     1 2.0150222, 1.6081385, 1.2600261,  .9661046,  .7217950,
     1  .5225175,  .3636947,  .2407501,  .1491035,  .0841764,
     1  .0413916,  .0161697,
     1  .0039326,  .0001014,  .0000209, -.0000055, -.0000049/)
      DOUBLE PRECISION, DIMENSION(n8) :: E8=(/
     1 2.8090209, 2.2393604, 1.7523957, 1.3416469, 1.0006275,
     1  .7228596,  .5018594,  .3311493,  .2042477,  .1146736,
     1  .0559460,  .0215842,
     1  .0051075,  .0000354, -.0000178, -.0000116, -.0000069,
     1 -.0000047/)

      DOUBLE PRECISION coefu(4,n0), coefg(4,n1), coefpu(4,n2), coef3(4,n3),
     1coef4(4,n4), coef5(4,n5), coef6(4,n6), coef7(4,n7), coef8(4,n8)
      data 
     1 nc0/1/,nc1/1/,nc2/1/,nc3/1/,nc4/1/,nc5/1/,nc6/1/,nc7/1/,nc8/1/

      ut0=0.d0

      if(i0.eq.1) then
      if(nc0.eq.1) then
      nc0=2
      do i=1,n0
      coefu(1,i)=Eu(i) 
      enddo
      coefu(2,n0)=-6.d0*Eu(n0)/ru(n0)
      call cubspl(ru,coefu,n0,0,1)
        r0u=ru(1)
        AAu=Eu(1)
        alphau=-coefu(2,1)/AAu
                   end if

cc      rr=dmax1(rr0,r0u)
        if(rr.lt.ru(n0)) then
        if(rr.lt.r0u) then
        ut0=AAu*exp(-alphau*(rr-r0u))
      ut0d=-alphau*ut0
                        else
      do i=1,n0-1
      ind=i
      if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 10
      enddo
   10      continue
      dr=rr-ru(ind)
      ut0=coefu(1,ind)+dr*(coefu(2,ind)
     1      +0.5d0*dr*(coefu(3,ind)+dr*coefu(4,ind)/3.d0))
      ut0d=coefu(2,ind) + dr*(coefu(3,ind) + 0.5d0*dr*coefu(4,ind))
                        end if
                         else
            ut0=Eu(n0)*(ru(n0)/rr)**6
        ut0d=-6.d0*ut0/rr
                   end if
      return
                  end if

      if(i0.eq.2) then
      if(nc1.eq.1) then
      nc1=2
      do i=1,n1
      coefg(1,i)=Eg(i)
      enddo
      coefg(2,n1)=-6.d0*Eg(n1)/ru(n1)
      call cubspl(ru,coefg,n1,0,1)
        r0g=ru(1)
        AAg=Eg(1)
        alphag=-coefg(2,1)/AAg
                   end if

cc        rr=dmax1(rr0,r0g)
        if(rr.lt.ru(n1)) then
        if(rr.lt.r0g) then
        ut0=AAg*exp(-alphag*(rr-r0g))
      ut0d=-alphag*ut0
                        else
        do i=1,n1-1
        ind=i
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 11
        enddo
   11   continue
        dr=rr-ru(ind)
        ut0=coefg(1,ind)+dr*(coefg(2,ind)
     1      +0.5d0*dr*(coefg(3,ind)+dr*coefg(4,ind)/3.d0))        
      ut0d=coefg(2,ind) + dr*(coefg(3,ind) + 0.5d0*dr*coefg(4,ind))
                        end if
                         else
        ut0=Eg(n1)*(ru(n1)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.3) then
      if(nc2.eq.1) then
      nc2=2
      do i=1,n2
      coefpu(1,i)=Epu(i)
      enddo
      coefpu(2,n2)=-6.d0*Epu(n2)/ru(n2)
      call cubspl(ru,coefpu,n2,0,1)
        r0pu=ru(1)
        AApu=Epu(1)
        alphapu=-coefpu(2,1)/AApu
                   end if
cc        rr=dmax1(rr0,r0pu)
      if(rr.lt.ru(n2)) then
       if(rr.lt.r0pu) then
        ut0=AApu*exp(-alphapu*(rr-r0pu))
      ut0d=-alphapu*ut0
                        else
        do i=1,n2-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 12
        enddo                     
   12    continue                     
        dr=rr-ru(ind)
        ut0=coefpu(1,ind)+dr*(coefpu(2,ind)
     1      +0.5d0*dr*(coefpu(3,ind)+dr*coefpu(4,ind)/3.d0))
      ut0d=coefpu(2,ind) + dr*(coefpu(3,ind) + 0.5d0*dr*coefpu(4,ind))
                        end if
                         else
        ut0=Epu(n2)*(ru(n2)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.4) then
      if(nc3.eq.1) then
      nc3=2
      do i=1,n3
      coef3(1,i)=E3(i)
      enddo
      coef3(2,n3)=-6.d0*E3(n3)/ru(n3)
      call cubspl(ru,coef3,n3,0,1)
        r03=ru(1)
        AA3=E3(1)
        alpha3=-coef3(2,1)/AA3
                   end if

cc        rr=dmax1(rr0,r03)
        if(rr.lt.ru(n3)) then
        if(rr.lt.r03) then
        ut0=AA3*exp(-alpha3*(rr-r03))
      ut0d=-alpha3*ut0
                        else
        do i=1,n3-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 13
        enddo                     
   13   continue                     
        dr=rr-ru(ind)
        ut0=coef3(1,ind)+dr*(coef3(2,ind)
     1  +0.5d0*dr*(coef3(3,ind)+dr*coef3(4,ind)/3.d0))
      ut0d=coef3(2,ind) + dr*(coef3(3,ind) + 0.5d0*dr*coef3(4,ind))
                        end if
                         else
        ut0=E3(n3)*(ru(n3)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.5) then
      if(nc4.eq.1) then
      nc4=2
      do i=1,n4
      coef4(1,i)=E4(i)
      enddo
      coef4(2,n4)=-6.d0*E4(n4)/ru(n4)
      call cubspl(ru,coef4,n4,0,1)
        r04=ru(1)
        AA4=E4(1)
        alpha4=-coef4(2,1)/AA4
                   end if

cc        rr=dmax1(rr0,r04)
        if(rr.lt.ru(n4)) then
        if(rr.lt.r04) then
        ut0=AA4*exp(-alpha4*(rr-r04))
      ut0d=-alpha4*ut0
                        else
        do i=1,n4-1
        ind=i
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 14
        enddo
   14   continue 
        dr=rr-ru(ind)
        ut0=coef4(1,ind)+dr*(coef4(2,ind)
     1  +0.5d0*dr*(coef4(3,ind)+dr*coef4(4,ind)/3.d0))
      ut0d=coef4(2,ind) + dr*(coef4(3,ind) + 0.5d0*dr*coef4(4,ind))
                        end if
                         else
        ut0=E4(n4)*(ru(n4)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.6) then
      if(nc5.eq.1) then
      nc5=2
      do i=1,n5
      coef5(1,i)=E5(i)
      enddo
      coef5(2,n5)=-6.d0*E5(n5)/ru(n5)
      call cubspl(ru,coef5,n5,0,1)
        r05=ru(1)
        AA5=E5(1)
        alpha5=-coef5(2,1)/AA5
                   end if

cc        rr=dmax1(rr0,r05)
        if(rr.lt.ru(n5)) then
        if(rr.lt.r05) then
        ut0=AA5*exp(-alpha5*(rr-r05))
      ut0d=-alpha5*ut0
                        else
        do i=1,n5-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 15                      
        enddo                     
   15   continue                      
        dr=rr-ru(ind)
        ut0=coef5(1,ind)+dr*(coef5(2,ind)
     1  +0.5d0*dr*(coef5(3,ind)+dr*coef5(4,ind)/3.d0))
      ut0d=coef5(2,ind) + dr*(coef5(3,ind) + 0.5d0*dr*coef5(4,ind))
                        end if
                         else
        ut0=E5(n5)*(ru(n5)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.7) then
      if(nc6.eq.1) then
      nc6=2
      do i=1,n6
      coef6(1,i)=E6(i)
      enddo
      coef6(2,n6)=-6.d0*E6(n6)/ru(n6)
      call cubspl(ru,coef6,n6,0,1)
        r06=ru(1)
        AA6=E6(1)
        alpha6=-coef6(2,1)/AA6
                   end if

cc        rr=dmax1(rr0,r06)
      if(rr.lt.ru(n6)) then
        if(rr.lt.r06) then
        ut0=AA6*exp(-alpha6*(rr-r06))
      ut0d=-alpha6*ut0
                        else
        do i=1,n6-1                                          
        ind=i                                          
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 16
        enddo                                          
   16   continue                                           
        dr=rr-ru(ind)
        ut0=coef6(1,ind)+dr*(coef6(2,ind)
     1  +0.5d0*dr*(coef6(3,ind)+dr*coef6(4,ind)/3.d0))
      ut0d=coef6(2,ind) + dr*(coef6(3,ind) + 0.5d0*dr*coef6(4,ind))
                        end if
                         else
        ut0=E6(n6)*(ru(n6)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.8) then
      if(nc7.eq.1) then
      nc7=2
      do i=1,n7
      coef7(1,i)=E7(i)
      enddo
      coef7(2,n7)=-6.d0*E7(n7)/ru(n7)
      call cubspl(ru,coef7,n7,0,1)
        r07=ru(1)
        AA7=E7(1)
        alpha7=-coef7(2,1)/AA7
                   end if

cc      rr=dmax1(rr0,r07)
      if(rr.lt.ru(n7)) then
        if(rr.lt.r07) then
        ut0=AA7*exp(-alpha7*(rr-r07))
      ut0d=-alpha7*ut0
                        else
        do i=1,n7-1                                                             
        ind=i                                                               
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 17                     
        enddo                                                               
   17   continue                                                                
        dr=rr-ru(ind)
        ut0=coef7(1,ind)+dr*(coef7(2,ind)
     1  +0.5d0*dr*(coef7(3,ind)+dr*coef7(4,ind)/3.d0))
      ut0d=coef7(2,ind) + dr*(coef7(3,ind) + 0.5d0*dr*coef7(4,ind))
                        end if
                         else
        ut0=E7(n7)*(ru(n7)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      if(i0.eq.9) then
      if(nc8.eq.1) then
      nc8=2
      do i=1,n8
      coef8(1,i)=E8(i)
      enddo
      coef8(2,n8)=-6.d0*E8(n8)/ru(n8)
      call cubspl(ru,coef8,n8,0,1)
        r08=ru(1)
        AA8=E8(1)
        alpha8=-coef8(2,1)/AA8
                   end if

cc      rr=dmax1(rr0,r08)
      if(rr.lt.ru(n8)) then
        if(rr.lt.r08) then
        ut0=AA8*exp(-alpha8*(rr-r08))
      ut0d=-alpha8*ut0
                        else
        do i=1,n8-1  
        ind=i        
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 18                              
        enddo        
   18   continue     
        dr=rr-ru(ind)
        ut0=coef8(1,ind)+dr*(coef8(2,ind)
     1  +0.5d0*dr*(coef8(3,ind)+dr*coef8(4,ind)/3.d0))
      ut0d=coef8(2,ind) + dr*(coef8(3,ind) + 0.5d0*dr*coef8(4,ind))
                        end if
                         else
        ut0=E8(n8)*(ru(n8)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      return
                  end if

      return 
      end

      subroutine polegn(n,x,P,P1)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION P(n),P1(n)
      P(1)=x
      P1(1)=1.d0
      if(n.gt.1) then
      Pn0=1.d0
      Pn=x
      Pn0d=0.d0
      Pnd=1.d0
      n1=n-1
      do i=1,n1
      ri=1.0D0*i
      Pn1=((2.d0*ri+1.d0)*x*Pn - ri*Pn0)/(ri+1.d0)
      Pn1d=((2.d0*ri+1.d0)*(Pn+x*Pnd) - ri*Pn0d)/(ri+1.d0)
      Pn0=Pn
      Pn=Pn1
      P(i+1)=Pn1
      Pn0d=Pnd
      Pnd=Pn1d
      P1(i+1)=Pn1d
      enddo
               end if
      return
      end

c------------------------------------------------------------------------

      DOUBLE PRECISION function ut01(rr,i0,ut01d)
      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT INTEGER(I-N)
      INTEGER N0, N1, N2, N3, N4, N5, N6, N7, N8
      save
      parameter(n0=13,n1=12,n2=12,n3=10,n4=10,n5=8,n6=8,n7=6,n8=6)

      DOUBLE PRECISION, DIMENSION(n0) :: ru=(/ 
     1 3.00, 3.25, 3.50, 3.75, 4.00, 
     1 4.25, 4.50, 5.00, 5.50, 6.00, 
     1 7.00, 8.00, 10.0/)

      DOUBLE PRECISION, DIMENSION(n0) :: Eu=(/
     1  .1220345,  .0289952, -.0023568, -.0104175, -.0105767,
     1 -.0086238, -.0065077, -.0034779, -.0018770, -.0010615,
     1 -.0003881, -.0001650, -.0000238/)
      DOUBLE PRECISION, DIMENSION(n1) :: Eg=(/
     1  .0275582,  .0126120,  .0047627,  .0012348, -.0001393,
     1 -.0005516, -.0005758, -.0003482, -.0001767, -.0000959,
     1 -.0000267, -.0000039/)
      DOUBLE PRECISION, DIMENSION(n2) :: Epu=(/
     1  .2295942,  .0839463,  .0271378,  .0064960, -.0002063,
     1 -.0018637, -.0018941, -.0011290, -.0005799, -.0003031,
     1 -.0000955, -.0000272/)
      DOUBLE PRECISION, DIMENSION(n3) :: E3=(/  
     1  .0347188,  .0157193,  .0067817,  .0027936,  .0010881,
     1  .0004036,  .0001308, -.0000088, -.0000128, -.0000039/)
      DOUBLE PRECISION, DIMENSION(n4) :: E4=(/  
     1  .0799795,  .0316737,  .0122230,  .0045691,  .0016355,
     1  .0005520,  .0001638,  .0000002, -.0000114, -.0000065/)
      DOUBLE PRECISION, DIMENSION(n5) :: E5=(/  
     1  .0078080,  .0032611,  .0013291,  .0005546,  .0002360,
     1  .0000854,  .0000302,  .0000003/)
      DOUBLE PRECISION, DIMENSION(n6) :: E6=(/  
     1  .0121196,  .0044330,  .0016626,  .0006155,  .0002380,
     1  .0000723,  .0000263,  .0000035/)
      DOUBLE PRECISION, DIMENSION(n7) :: E7=(/  
     1  .0012897,  .0004052,  .0001378,  .0000390,  .0000099,
     1  .0000014/)
      DOUBLE PRECISION, DIMENSION(n8) :: E8=(/
     1  .0014554,  .0004653,  .0001494,  .0000507,  .0000087,
     1  .0000041/)

      DOUBLE PRECISION coefu(4,n0), coefg(4,n1), coefpu(4,n2), coef3(4,n3),
     1coef4(4,n4), coef5(4,n5), coef6(4,n6), coef7(4,n7), coef8(4,n8)
      data 
     1 nc0/1/,nc1/1/,nc2/1/,nc3/1/,nc4/1/,nc5/1/,nc6/1/,nc7/1/,nc8/1/

      ut0=0.d0
      ut0d=0.d0

      if(i0.eq.1) then
      if(nc0.eq.1) then
      nc0=2
      do i=1,n0
      coefu(1,i)=Eu(i) 
      enddo
      coefu(2,n0)=-6.d0*Eu(n0)/ru(n0)
      call cubspl(ru,coefu,n0,0,1)
        r0u=ru(1)
        AAu=Eu(1)
        alphau=-coefu(2,1)/AAu
                   end if

cc      rr=dmax1(rr0,r0u)
        if(rr.lt.ru(n0)) then
        if(rr.lt.r0u) then
        ut0=AAu*exp(-alphau*(rr-r0u))
      ut0d=-alphau*ut0
                        else
      do i=1,n0-1
      ind=i
      if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 10
      enddo
   10      continue
      dr=rr-ru(ind)
      ut0=coefu(1,ind)+dr*(coefu(2,ind)
     1      +0.5d0*dr*(coefu(3,ind)+dr*coefu(4,ind)/3.d0))
      ut0d=coefu(2,ind) + dr*(coefu(3,ind) + 0.5d0*dr*coefu(4,ind))
                        end if
                         else
            ut0=Eu(n0)*(ru(n0)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
      ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.2) then
      if(nc1.eq.1) then
      nc1=2
      do i=1,n1
      coefg(1,i)=Eg(i)
      enddo
      coefg(2,n1)=-6.d0*Eg(n1)/ru(n1)
      call cubspl(ru,coefg,n1,0,1)
        r0g=ru(1)
        AAg=Eg(1)
        alphag=-coefg(2,1)/AAg
                   end if

cc        rr=dmax1(rr0,r0g)
        if(rr.lt.ru(n1)) then
        if(rr.lt.r0g) then
        ut0=AAg*exp(-alphag*(rr-r0g))
      ut0d=-alphag*ut0
                        else
        do i=1,n1-1
        ind=i
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 11
        enddo
   11   continue
        dr=rr-ru(ind)
        ut0=coefg(1,ind)+dr*(coefg(2,ind)
     1      +0.5d0*dr*(coefg(3,ind)+dr*coefg(4,ind)/3.d0))        
      ut0d=coefg(2,ind) + dr*(coefg(3,ind) + 0.5d0*dr*coefg(4,ind))
                        end if
                         else
        ut0=Eg(n1)*(ru(n1)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.3) then
      if(nc2.eq.1) then
      nc2=2
      do i=1,n2
      coefpu(1,i)=Epu(i)
      enddo
      coefpu(2,n2)=-6.d0*Epu(n2)/ru(n2)
      call cubspl(ru,coefpu,n2,0,1)
        r0pu=ru(1)
        AApu=Epu(1)
        alphapu=-coefpu(2,1)/AApu
                   end if
cc        rr=dmax1(rr0,r0pu)
      if(rr.lt.ru(n2)) then
       if(rr.lt.r0pu) then
        ut0=AApu*exp(-alphapu*(rr-r0pu))
      ut0d=-alphapu*ut0
                        else
        do i=1,n2-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 12
        enddo                     
   12    continue                     
        dr=rr-ru(ind)
        ut0=coefpu(1,ind)+dr*(coefpu(2,ind)
     1      +0.5d0*dr*(coefpu(3,ind)+dr*coefpu(4,ind)/3.d0))
      ut0d=coefpu(2,ind) + dr*(coefpu(3,ind) + 0.5d0*dr*coefpu(4,ind))
                        end if
                         else
        ut0=Epu(n2)*(ru(n2)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.4) then
      if(nc3.eq.1) then
      nc3=2
      do i=1,n3
      coef3(1,i)=E3(i)
      enddo
      coef3(2,n3)=-6.d0*E3(n3)/ru(n3)
      call cubspl(ru,coef3,n3,0,1)
        r03=ru(1)
        AA3=E3(1)
        alpha3=-coef3(2,1)/AA3
                   end if

cc        rr=dmax1(rr0,r03)
        if(rr.lt.ru(n3)) then
        if(rr.lt.r03) then
        ut0=AA3*exp(-alpha3*(rr-r03))
      ut0d=-alpha3*ut0
                        else
        do i=1,n3-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 13
        enddo                     
   13   continue                     
        dr=rr-ru(ind)
        ut0=coef3(1,ind)+dr*(coef3(2,ind)
     1  +0.5d0*dr*(coef3(3,ind)+dr*coef3(4,ind)/3.d0))
      ut0d=coef3(2,ind) + dr*(coef3(3,ind) + 0.5d0*dr*coef3(4,ind))
                        end if
                         else
        ut0=E3(n3)*(ru(n3)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.5) then
      if(nc4.eq.1) then
      nc4=2
      do i=1,n4
      coef4(1,i)=E4(i)
      enddo
      coef4(2,n4)=-6.d0*E4(n4)/ru(n4)
      call cubspl(ru,coef4,n4,0,1)
        r04=ru(1)
        AA4=E4(1)
        alpha4=-coef4(2,1)/AA4
                   end if

cc        rr=dmax1(rr0,r04)
        if(rr.lt.ru(n4)) then
        if(rr.lt.r04) then
        ut0=AA4*exp(-alpha4*(rr-r04))
      ut0d=-alpha4*ut0
                        else
        do i=1,n4-1
        ind=i
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 14
        enddo
   14   continue 
        dr=rr-ru(ind)
        ut0=coef4(1,ind)+dr*(coef4(2,ind)
     1  +0.5d0*dr*(coef4(3,ind)+dr*coef4(4,ind)/3.d0))
      ut0d=coef4(2,ind) + dr*(coef4(3,ind) + 0.5d0*dr*coef4(4,ind))
                        end if
                         else
        ut0=E4(n4)*(ru(n4)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.6) then
      if(nc5.eq.1) then
      nc5=2
      do i=1,n5
      coef5(1,i)=E5(i)
      enddo
      coef5(2,n5)=-6.d0*E5(n5)/ru(n5)
      call cubspl(ru,coef5,n5,0,1)
        r05=ru(1)
        AA5=E5(1)
        alpha5=-coef5(2,1)/AA5
                   end if

cc        rr=dmax1(rr0,r05)
        if(rr.lt.ru(n5)) then
        if(rr.lt.r05) then
        ut0=AA5*exp(-alpha5*(rr-r05))
      ut0d=-alpha5*ut0
                        else
        do i=1,n5-1                     
        ind=i                     
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 15                      
        enddo                     
   15   continue                      
        dr=rr-ru(ind)
        ut0=coef5(1,ind)+dr*(coef5(2,ind)
     1  +0.5d0*dr*(coef5(3,ind)+dr*coef5(4,ind)/3.d0))
      ut0d=coef5(2,ind) + dr*(coef5(3,ind) + 0.5d0*dr*coef5(4,ind))
                        end if
                         else
        ut0=E5(n5)*(ru(n5)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.7) then
      if(nc6.eq.1) then
      nc6=2
      do i=1,n6
      coef6(1,i)=E6(i)
      enddo
      coef6(2,n6)=-6.d0*E6(n6)/ru(n6)
      call cubspl(ru,coef6,n6,0,1)
        r06=ru(1)
        AA6=E6(1)
        alpha6=-coef6(2,1)/AA6
                   end if

cc        rr=dmax1(rr0,r06)
      if(rr.lt.ru(n6)) then
        if(rr.lt.r06) then
        ut0=AA6*exp(-alpha6*(rr-r06))
      ut0d=-alpha6*ut0
                        else
        do i=1,n6-1                                          
        ind=i                                          
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 16
        enddo                                          
   16   continue                                           
        dr=rr-ru(ind)
        ut0=coef6(1,ind)+dr*(coef6(2,ind)
     1  +0.5d0*dr*(coef6(3,ind)+dr*coef6(4,ind)/3.d0))
      ut0d=coef6(2,ind) + dr*(coef6(3,ind) + 0.5d0*dr*coef6(4,ind))
                        end if
                         else
        ut0=E6(n6)*(ru(n6)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.8) then
      if(nc7.eq.1) then
      nc7=2
      do i=1,n7
      coef7(1,i)=E7(i)
      enddo
      coef7(2,n7)=-6.d0*E7(n7)/ru(n7)
      call cubspl(ru,coef7,n7,0,1)
        r07=ru(1)
        AA7=E7(1)
        alpha7=-coef7(2,1)/AA7
                   end if

cc      rr=dmax1(rr0,r07)
      if(rr.lt.ru(n7)) then
        if(rr.lt.r07) then
        ut0=AA7*exp(-alpha7*(rr-r07))
      ut0d=-alpha7*ut0
                        else
        do i=1,n7-1                                                             
        ind=i                                                               
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 17                     
        enddo                                                               
   17   continue                                                                
        dr=rr-ru(ind)
        ut0=coef7(1,ind)+dr*(coef7(2,ind)
     1  +0.5d0*dr*(coef7(3,ind)+dr*coef7(4,ind)/3.d0))
      ut0d=coef7(2,ind) + dr*(coef7(3,ind) + 0.5d0*dr*coef7(4,ind))
                        end if
                         else
        ut0=E7(n7)*(ru(n7)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

      if(i0.eq.9) then
      if(nc8.eq.1) then
      nc8=2
      do i=1,n8
      coef8(1,i)=E8(i)
      enddo
      coef8(2,n8)=-6.d0*E8(n8)/ru(n8)
      call cubspl(ru,coef8,n8,0,1)
        r08=ru(1)
        AA8=E8(1)
        alpha8=-coef8(2,1)/AA8
                   end if

cc      rr=dmax1(rr0,r08)
      if(rr.lt.ru(n8)) then
        if(rr.lt.r08) then
        ut0=AA8*exp(-alpha8*(rr-r08))
      ut0d=-alpha8*ut0
                        else
        do i=1,n8-1  
        ind=i        
        if(ru(i).le.rr.and.rr.lt.ru(i+1)) go to 18                              
        enddo        
   18   continue     
        dr=rr-ru(ind)
        ut0=coef8(1,ind)+dr*(coef8(2,ind)
     1  +0.5d0*dr*(coef8(3,ind)+dr*coef8(4,ind)/3.d0))
      ut0d=coef8(2,ind) + dr*(coef8(3,ind) + 0.5d0*dr*coef8(4,ind))
                        end if
                         else
        ut0=E8(n8)*(ru(n8)/rr)**6
      ut0d=-6.d0*ut0/rr
                         end if
        ut01=ut0
      ut01d=ut0d
      return
                  end if

        ut01=ut0
      ut01d=ut0d
      return 
      end


