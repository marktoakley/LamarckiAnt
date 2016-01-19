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
C*******************************************************************************
C
C This subroutine calculates the Al/Ni Farkas energy and first derivative
C
C*******************************************************************************

      SUBROUTINE FARKAS(X,gradFark,eFark,GRADT,NATOMS)
C      USE commons
      implicit none
C     integer amax

      integer j1, j2, lfden, lfembed, lfpair, natoms
      double precision fdenx(500), fdeny(500), fdeny2(500),
     1                 fembedx(500), fembedy(500), fembedy2(500),
     2                 fpairx(500), fpairy(500), fpairy2(500),
     3                 vpair, f, rho(NATOMS), x(3*NATOMS), eFark, 
     4                 rhotemp, dist, rcut, gradFark(3*NATOMS), df(NATOMS), 
     5                 dvtemp, dvtemp1, dvtemp2, dvtemp3, drho,
     6                 dvpair, vptemp, rcutsq, dvplo, rmin
      logical GRADT
      common /cfarkas/ lfden, lfembed, lfpair, 
     1                fdenx, fdeny, fdeny2,
     2                fembedx, fembedy, fembedy2,
     3                fpairx, fpairy, fpairy2

C Calculate energy

      rcut=fpairx(lfpair)
      rmin=fpairx(1)
      rcutsq=rcut**2
      call dsplint(fpairx,fpairy,fpairy2,lfpair,fpairx(1),dvplo)
      eFark=0.0d0
      DO 22 J1=1,NATOMS
         RHO(J1)=0.0d0
         vpair=0.0d0
         DO 23 J2=1,NATOMS
            if (j1.ne.j2) then
               DIST=(X(3*(J2-1)+1)-X(3*(J1-1)+1))**2 +
     1              (X(3*(J2-1)+2)-X(3*(J1-1)+2))**2 +
     2              (X(3*(J2-1)+3)-X(3*(J1-1)+3))**2
               if (dist.lt.rcutsq) then
                 dist=dsqrt(dist)
                 if (dist.gt.rmin) then
                   call splint(fdenx,fdeny,fdeny2,lfden,dist,
     1                         rhotemp)
                 else
                   rhotemp=fdeny(1)
                 endif
                 RHO(J1)=RHO(J1)+rhotemp
                 if (j1.lt.j2) then
                   if (dist.gt.rmin) then
                     call splint(fpairx,fpairy,fpairy2,lfpair,
     1                         dist,vptemp)
                   else
                     vptemp=fpairy(1)+(dist-rmin)*dvplo
C                     print*,'ED: r outside range', dist
                   endif
                   vpair=vpair+vptemp
                 endif
               endif
            endif
23       CONTINUE
         if (rho(j1).lt.fembedx(lfembed)) then
           call splint(fembedx,fembedy,fembedy2,lfembed,rho(j1),f)
         else
C           print*, 'ED: rho outside of range', rho(j1)
           f=fembedy(lfembed)
         endif
         eFark=eFark+vpair+f
22    CONTINUE

C      print*,eFark

C Calculate gradient

      IF (.not.GRADT) RETURN

      DO 30 J1=1,NATOMS
        if (rho(j1).lt.fembedx(lfembed)) then
          call dsplint(fembedx,fembedy,fembedy2,lfembed,rho(j1),
     1               df(j1))
        else
C          print*, 'ED: rho outside of range', rho(j1)
          df(j1)=0.0d0
        endif
30    CONTINUE    
      
      DO 50 J1=1,NATOMS
         DVTEMP1=0.0d0
         DVTEMP2=0.0d0
         DVTEMP3=0.0d0
         DO 40 J2=1,NATOMS
            if (j1.ne.j2) then
               DIST=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1              ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2              ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
               if (dist.lt.rcutsq) then
                 dist=dsqrt(dist)
                 if (dist.gt.rmin) then
                   call dsplint(fpairx,fpairy,fpairy2,lfpair,dist,
     1                        dvpair)
                   call dsplint(fdenx,fdeny,fdeny2,lfden,dist,drho)
                 else
C                   print*, 'ED: r outside of range', dist
                   dvpair=dvplo 
                   drho=0.0d0
                 endif
                 dvtemp=(dvpair+(df(j1)+df(j2))*drho)/dist
                 DVTEMP1=DVTEMP1+DVTEMP*(X(3*(J1-1)+1)-X(3*(J2-1)+1))
                 DVTEMP2=DVTEMP2+DVTEMP*(X(3*(J1-1)+2)-X(3*(J2-1)+2))
                 DVTEMP3=DVTEMP3+DVTEMP*(X(3*(J1-1)+3)-X(3*(J2-1)+3))
               endif
            endif
40       CONTINUE
         gradFark(3*(J1-1)+1)=DVTEMP1
         gradFark(3*(J1-1)+2)=DVTEMP2
         gradFark(3*(J1-1)+3)=DVTEMP3
50    CONTINUE

      RETURN
      END

C*******************************************************************************
C
C This subroutine initializes the arrays of 2nd derivatives for the 
C cubic spline interpolations for Ni.
C
C*******************************************************************************

      SUBROUTINE NIINIT

      implicit none

      integer i, lfden, lfembed, lfpair
      double precision fdenx(500), fdeny(500), fdeny2(500),
     1                 fembedx(500), fembedy(500), fembedy2(500),
     2                 fpairx(500), fpairy(500), fpairy2(500),
     3                 dylo, dyhi
      common /cfarkas/ lfden, lfembed, lfpair, 
     1                fdenx, fdeny, fdeny2,
     2                fembedx, fembedy, fembedy2,
     3                fpairx, fpairy, fpairy2

      open(unit=8,file='ni.den',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 20 i=1,500
        read(8,*,end=30) fdenx(i), fdeny(i)
20    continue
      print*,'WARNING: increase the dimensions in NIINIT'
      stop

30    lfden=i-1
      close(unit=8)
      print*,'ni.den has ', lfden, ' entries'
      print*,'First derivatives at end points ', dylo, dyhi

      call splinegmin(fdenx,fdeny,lfden,dylo,dyhi,fdeny2)

      open(unit=8,file='ni.pair',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 70 i=1,500
        read(8,*,end=80) fpairx(i), fpairy(i)
70    continue
      print*,'WARNING: increase the dimensions in NIINIT'
      stop

80    lfpair=i-1
      close(unit=8)
      print*,'ni.pair has ', lfpair, ' entries'
      print*,'second derivatives at end points ', dylo, dyhi

      call splinegmin(fpairx,fpairy,lfpair,dylo,dyhi,fpairy2)

      open(unit=8,file='ni.embed',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 120 i=1,500
        read(8,*,end=130) fembedx(i), fembedy(i)
120   continue
      print*,'WARNING: increase 500'
      stop

130   lfembed=i-1
      close(unit=8)
      print*,'ni.embed has ', lfembed, ' entries'
      print*,'First derivatives at end points ',dylo,dyhi

      call splinegmin(fembedx,fembedy,lfembed,dylo,dyhi,
     1            fembedy2)

      RETURN
      END

C*******************************************************************************
C
C This subroutine initializes the arrays of 2nd derivatives for the 
C cubic spline interpolations for Al.
C
C*******************************************************************************

      SUBROUTINE ALINIT

      implicit none

      integer i, lfden, lfembed, lfpair
      double precision fdenx(500), fdeny(500), fdeny2(500),
     1                 fembedx(500), fembedy(500), fembedy2(500),
     2                 fpairx(500), fpairy(500), fpairy2(500), 
     3                 dylo, dyhi
      common /cfarkas/ lfden, lfembed, lfpair, 
     1                fdenx, fdeny, fdeny2,
     2                fembedx, fembedy, fembedy2,
     3                fpairx, fpairy, fpairy2

      open(unit=8,file='al.den',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 20 i=1,500
        read(8,*,end=30) fdenx(i), fdeny(i)
20    continue
      print*,'WARNING: increase the dimensions in ALINIT'
      stop

30    lfden=i-1
      close(unit=8)
      print*,'al.den has ', lfden, ' entries'
      print*,'First derivatives at end points ',dylo,dyhi

      call splinegmin(fdenx,fdeny,lfden,dylo,dyhi,fdeny2)

      open(unit=8,file='al.pair',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 70 i=1,500
        read(8,*,end=80) fpairx(i), fpairy(i)
70    continue
      print*,'WARNING: increase the dimensions in ALINIT'
      stop

80    lfpair=i-1
      close(unit=8)
      print*,'al.pair has ', lfpair, ' entries'
      print*,'First derivatives at end points ',dylo,dyhi

      call splinegmin(fpairx,fpairy,lfpair,dylo,dyhi,fpairy2)

      open(unit=8,file='al.embed',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*) dylo, dyhi
      read(8,*)
      do 120 i=1,500
        read(8,*,end=130) fembedx(i), fembedy(i)
120   continue
      print*,'WARNING: increase the dimensions in ALINIT'
      stop

130   lfembed=i-1
      close(unit=8)
      print*,'al.embed has ', lfembed, ' entries'
      print*,'First derivatives at end points ',dylo,dyhi

      call splinegmin(fembedx,fembedy,lfembed,dylo,dyhi,fembedy2)

      RETURN
      END

      SUBROUTINE splinegmin(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=5000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER khi,klo
      DOUBLE PRECISION a,b,h,pos

      if ((x.lt.xa(1)).or.(x.gt.xa(n))) then
        print*, 'WARNING: x out of range in SPLINT', x, xa(1),xa(n)
        stop
      endif

C works if intervals equally spaced 

      pos=1+(n-1)*(x-xa(1))/(xa(n)-xa(1))
      klo=int(pos)
      khi=int(pos)+1

C Check if on end point: unlikely but possible

      if (klo.eq.n) then
        klo=n-1
        khi=n
      endif

C interval bisection if not
C      klo=1
C      khi=n
C1     if (khi-klo.gt.1) then
C        k=(khi+klo)/2
C        if(xa(k).gt.x)then
C          khi=k
C        else
C          klo=k
C        endif
C      goto 1
C      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) print*, 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
      END

      SUBROUTINE dsplint(xa,ya,y2a,n,x,dydx)
      INTEGER n
      DOUBLE PRECISION x,dydx,xa(n),y2a(n),ya(n)
      INTEGER khi,klo
      DOUBLE PRECISION a,b,h,pos

      if ((x.lt.xa(1)).or.(x.gt.xa(n))) then
        print*, 'WARNING: x out of range in DSPLINT', x, xa(1),xa(n)
        stop
      endif

C works if intervals equally spaced 

      pos=1+(n-1)*(x-xa(1))/(xa(n)-xa(1))
      klo=int(pos)
      khi=int(pos)+1

C Check if on end point: unlikely but possible

      if (klo.eq.n) then
        klo=n-1
        khi=n
      endif

C interval bisection if not
C      klo=1
C      khi=n
C1     if (khi-klo.gt.1) then
C        k=(khi+klo)/2
C        if(xa(k).gt.x)then
C          khi=k
C        else
C          klo=k
C        endif
C      goto 1
C      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) print*, 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      dydx=(ya(khi)-ya(klo))/(xa(khi)-xa(klo))-
     1     (3*a**2-1)*(xa(khi)-xa(klo))*y2a(klo)/6.0d0+
     2     (3*b**2-1)*(xa(khi)-xa(klo))*y2a(khi)/6.0d0
      return
      END
