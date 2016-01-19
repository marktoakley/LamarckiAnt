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
      SUBROUTINE ODESD(ITMAX,VARS,MFLAG,NSTP,ENERGY,NP)
      USE commons
      use odesdsavear
      IMPLICIT NONE
      INTEGER I, NOK, NBAD, NSTP, ITMAX, NOPT, NP
      LOGICAL MFLAG
      DOUBLE PRECISION VARS(3*NATOMS), YSCAL(3*NATOMS), ENERGY, MXSTP, STRY,
     1                 DYDX(3*NATOMS), STRYDID, EPS, STRYNEXT, DUMMY
      DOUBLE PRECISION TINY, STRYMIN, SLENGTH
      PARAMETER (TINY=1.D-30)
      COMMON /BSNEW/ SLENGTH,NOK,NBAD,EPS

      STRYMIN=0.0D0
      NOPT=3*NATOMS
      MXSTP=0.00001D0
      STRY=MXSTP
      NSTP=1
10    CALL POTENTIAL(VARS,DYDX,ENERGY,.TRUE.,.FALSE.)

      !FIXIMAGE=.TRUE.

      IF (RMS.LT.GMAX) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            RETURN
         ENDIF
      ENDIF

      IF (NSTP.EQ.ITMAX) THEN
         MFLAG=.FALSE.
         FIXIMAGE=.FALSE.
         RETURN
      ENDIF

      DUMMY=RMS*SQRT(1.0D0*NOPT)
      DO I=1,NOPT
         YSCAL(I)=ABS(VARS(I))+ABS(STRY*DYDX(I))+TINY
         DYDX(I)=-DYDX(I)/MAX(DUMMY,1.0D0)
C        DYDX(I)=-DYDX(I)
      ENDDO

      IF (BSMIN) CALL BSSTEP(VARS,DYDX,SLENGTH,STRY,EPS,YSCAL,STRYDID,STRYNEXT,ENERGY,NOPT)
      IF (RKMIN) CALL RKQS(VARS,DYDX,SLENGTH,STRY,EPS,YSCAL,STRYDID,STRYNEXT,ENERGY,NOPT)

      WRITE(*,'(A,G20.10,A,G20.10)') 'Step length=',STRYDID,' next estimated step size=',STRYNEXT

      IF (STRYDID.EQ.STRY) THEN
         NOK=NOK+1
      ELSE
         NBAD=NBAD+1
      ENDIF
      PRINT '(A,2G20.10,2I8)','STRYDID,STRY,NOK,NBAD=',STRYDID,STRY,NOK,NBAD
      IF (ABS(STRYNEXT).LT.STRYMIN) THEN
         PRINT*, ' WARNING, stepsize < 0 in odesd'
         OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
         PRINT*,' intractable discontinuity in main - quit '
         WRITE(96,'(A)') 'intractable discontinuity'
         CLOSE(96)
         WRITE(*,'(A,2F20.10,A,I6,A,F15.10)') ' Energy and RMS force=',ENERGY,RMS
         STOP
      ENDIF
      IF (DEBUG) WRITE(*,'(A,2F20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',NSTP,' SD steps'

C     STRY=MIN(STRYNEXT,MXSTP)
      STRY=STRYNEXT
      MFLAG=.FALSE.
      FIXIMAGE=.FALSE.
      NSTP=NSTP+1
      GOTO 10

      RETURN
      END

      SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,energy,nv)
      USE commons
      use odesdsavear
      IMPLICIT NONE
      INTEGER nv,KMAXX,IMAX
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(3*NATOMS),y(3*NATOMS),yscal(3*NATOMS)
     *,SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     *REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.9d0)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      DOUBLE PRECISION eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
     *,xest,xnew,energy,eold,VNEW(3*NATOMS),
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(3*NATOMS),ysav(3*NATOMS),
     *yseq(3*NATOMS)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      allocate(qcol(3*NATOMS,IMAX),d(3*NATOMS,IMAX))

C     PRINT*,'in bsstep, x,eps,epsold=',x,eps,epsold
      eold=energy
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if ((xnew.eq.x).OR.(EPS.LT.1.0D-20)) THEN
           print*,'stepsize underflow in bsstep'
           OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
           PRINT*,' intractable discontinuity in bssstep - quit '
           WRITE(96,'(A)') 'intractable discontinuity'
           CLOSE(96)
!          STOP
        ENDIF
        call mmid(ysav,dydx,x,h,nseq(k),yseq,nv)
        xest=(h/nseq(k))**2
        call rzextr(k,xest,yseq,y,yerr,nv,natoms)
C       call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          IF (ERRMAX.LT.1.0D0) THEN
             CALL POTENTIAL(Y,VNEW,ENERGY,.FALSE.,.FALSE.)
             IF (ENERGY.GT.EOLD+1.0D-12) THEN
                ERRMAX=2.0D0
                EPS=EPS/2.0D0
                WRITE(*,'(A,4G20.10)') 'Y,ENERGY,EOLD,EPS=',Y(1),ENERGY,EOLD,EPS
             ELSE
C               EPS=EPS*1.1D0
             ENDIF
          ENDIF
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

      SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,NOPT)
      USE commons
      IMPLICIT NONE
      INTEGER nstep, NOPT
      DOUBLE PRECISION htot,xs,dydx(3*NATOMS),y(3*NATOMS),yout(3*NATOMS),ENERGY
      INTEGER i,n
      DOUBLE PRECISION h,h2,swap,x,ym(3*NATOMS),yn(3*NATOMS)

C     PRINT*,'in mmid'
      h=htot/nstep
      do 11 i=1,NOPT
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      CALL POTENTIAL(YN,YOUT,ENERGY,.TRUE.,.FALSE.)
      RMS=RMS*SQRT(1.0D0*NOPT)
      DO I=1,NOPT
         YOUT(I)=-YOUT(I)/MAX(RMS,1.0D0)
C        YOUT(I)=-YOUT(I)
      ENDDO
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,NOPT
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        CALL POTENTIAL(YN,YOUT,ENERGY,.TRUE.,.FALSE.)
        RMS=RMS*SQRT(1.0D0*NOPT)
        DO I=1,NOPT
           YOUT(I)=-YOUT(I)/MAX(RMS,1.0D0)
C          YOUT(I)=-YOUT(I)
        ENDDO
13    continue
      do 14 i=1,NOPT
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END

      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv,natoms)
      use odesdsavear
      IMPLICIT NONE
      INTEGER iest,nv,IMAX, natoms
      DOUBLE PRECISION xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13)
      INTEGER j,k1
      DOUBLE PRECISION delta,f1,f2,q,x(IMAX)
      SAVE x

C     PRINT*,'in pzextr'
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j,1)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j,1)-q
            dy(j)=f1*delta
            d(j,1)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
      SUBROUTINE rzextr(iest,xest,yest,yz,dy,nv,natoms)
      use odesdsavear
      IMPLICIT NONE
      INTEGER iest,nv,IMAX, natoms
      DOUBLE PRECISION xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13)
      INTEGER j,k
      DOUBLE PRECISION b,b1,c,ddy,v,yy,fx(IMAX),x(IMAX)
      SAVE x

C     PRINT*,'in rzextr'
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.d0) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      END

      SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,energy,n)
      USE commons
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(3*NATOMS),y(3*NATOMS),yscal(3*NATOMS)
CU    USES rkck
      INTEGER i
      DOUBLE PRECISION errmax,h,xnew,yerr(3*NATOMS),ytemp(3*NATOMS),SAFETY,PGROW
     *,PSHRNK,
     *ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
      DOUBLE PRECISION EOLD, ENERGY, VNEW(3*NATOMS)
      
      eold=energy
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps

      PRINT*,'ERRMAX,h in rkqs=',ERRMAX,H
      IF (ERRMAX.LT.1.0D0) THEN
         CALL POTENTIAL(YTEMP,VNEW,ENERGY,.FALSE.,.FALSE.)
         PRINT*,'ENERGY,EOLD=',ENERGY,EOLD

         IF (ENERGY.GT.EOLD+1.0D-12) THEN
            ERRMAX=2.0D0
C           EPS=EPS/2.0D0
            WRITE(*,'(A,4G20.10)') 'YTEMP(1),ENERGY,EOLD,EPS=',YTEMP(1),ENERGY,EOLD,EPS
         ELSE
C           EPS=EPS*1.1D0
         ENDIF
      ENDIF

      if(errmax.gt.1.d0)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1d0*h)then
          h=.1d0*h
        endif
        xnew=x+h
        if ((xnew.eq.x).OR.(EPS.LT.1.0D-20)) THEN
           PRINT*,' stepsize underflow in rkqs'
           OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
           PRINT '(A,3G20.10)',' intractable discontinuity in rkqs x,xnew,eps=',x,xnew,eps
           PRINT '(A,3G20.10)',' x,xnew=',x,xnew
           PRINT '(A,3G20.10)',' eps=',eps
           WRITE(96,'(A)') 'intractable discontinuity'
           CLOSE(96)
!          STOP
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr)
      USE commons
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION h,x,dydx(n),y(n),yerr(n),yout(n)
      INTEGER i
      DOUBLE PRECISION ak2(3*NATOMS),ak3(3*NATOMS),ak4(3*NATOMS),ak5(3*NATOMS),ak6(3*NATOMS)
     *,
     *ytemp(3*NATOMS),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31
     *=3.d0/40.d0,
     *B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52
     *=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0
     */512.d0,
     *B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1
     *=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0
     */27648.d0,
     *DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0
     */14336.d0,
     *DC6=C6-.25d0)
      DOUBLE PRECISION E1, E2, E3,E4,E5

      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
C     call derivs(x+A2*h,ytemp,ak2)

      CALL POTENTIAL(YTEMP,AK2,E1,.TRUE.,.FALSE.)

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
         AK2(I)=-AK2(I)/MAX(RMS,1.0D0)
      ENDDO

      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
C     call derivs(x+A3*h,ytemp,ak3)

      CALL POTENTIAL(YTEMP,AK3,E2,.TRUE.,.FALSE.)

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
         AK3(I)=-AK3(I)/MAX(RMS,1.0D0)
      ENDDO

      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
C     call derivs(x+A4*h,ytemp,ak4)

      CALL POTENTIAL(YTEMP,AK4,E3,.TRUE.,.FALSE.)

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
         AK4(I)=-AK4(I)/MAX(RMS,1.0D0)
      ENDDO

      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
C     call derivs(x+A5*h,ytemp,ak5)

      CALL POTENTIAL(YTEMP,AK5,E4,.TRUE.,.FALSE.)

      RMS=RMS*SQRT(1.0D0*N)

       DO I=1,N
          AK5(I)=-AK5(I)/MAX(RMS,1.0D0)
       ENDDO

      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
15    continue
C     call derivs(x+A6*h,ytemp,ak6)

      CALL POTENTIAL(YTEMP,AK6,E5,.TRUE.,.FALSE.)

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
         AK6(I)=-AK6(I)/MAX(RMS,1.0D0)
      ENDDO

      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
17    continue
      PRINT '(A,5F20.10)','energies=',E1,E2,E3,E4,E5
      return
      END

