!Computes Carlson's elliptic integral of the rst kind, RF (x; y; z). x, y, and z must be
!nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
!underflow limit, BIG at most one fth the machine overflow limit.

SUBROUTINE Carlsonrf(rf,x,y,z)

      DOUBLE PRECISION :: rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      DOUBLE PRECISION :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

      ERRTOL=.0025D0
      TINY=1.5e-38
      BIG=3.E37
      THIRD=1.D0/3.
      C1=1.D0/24.
      C2=.1D0
      C3=3.D0/44.
      C4=1.D0/14.

      IF(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or. max(x,y,z).gt.BIG) THEN
         PRINT *, x, y, z
         STOP 'invalid arguments in rf'
      ENDIF

      xt=x
      yt=y
      zt=z
 1    continue

      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      ave=THIRD*(xt+yt+zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave

      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return

    END SUBROUTINE Carlsonrf


!Computes Carlson's elliptic integral of the second kind, RD(x; y; z). x and y must be
!nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
!the negative 2/3 power of the machine overflow limit. BIG must be at most 0:1ERRTOL
!times the negative 2/3 power of the machine underflow limit.

SUBROUTINE Carlsonrd(rd,x,y,z)

      DOUBLE PRECISION :: rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      DOUBLE PRECISION :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt

      ERRTOL=.0015D0
      TINY=1.e-25
      BIG=4.5E21
      C1=3.D0/14.
      C2=1.D0/6.
      C3=9.D0/22.
      C4=3.D0/26.
      C5=.25D0*C3
      C6=1.5D0*C4

      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG) STOP 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
 1    continue

      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      sum=sum+fac/(sqrtz*(zt+alamb))
      fac=.25*fac
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      ave=.2*(xt+yt+3.*zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
     
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
      return 
    
END SUBROUTINE Carlsonrd

SUBROUTINE EllipIntegral (Felint, Selint, phi, k)

      DOUBLE PRECISION :: Felint, Selint, phi, k, rf, rd, x, y, z
      DOUBLE PRECISION :: pi, CFelint, CSelint
      INTEGER :: N2PI, PosNeg

      CALL CompleteEllipIntegral (CFelint, CSelint, k)

      pi = 4.0D0 * ATAN(1.0D0)
      IF (phi < 0.0D0) THEN
         PosNeg = -1
         phi = -1.0D0*phi
      ENDIF
      N2PI = FLOOR(phi/2.0D0/pi)
      phi = phi - N2PI*2.0D0*pi
      
      IF ( (phi < 0.000001D0) .AND. (phi > -0.000001D0) ) THEN 
         Felint = 0.0D0
         Selint = 0.0D0
      ELSEIF ( (phi-pi < 0.000001D0) .AND. (phi-pi > -0.000001D0) ) THEN 
         Felint = 2.0D0 * CFelint
         Selint = 2.0D0 * CSelint
      ELSEIF ( (phi-2.0D0*pi < 0.000001D0) .AND. (phi-2.0D0*pi > -0.000001D0) ) THEN 
         Felint = 4.0D0 * CFelint
         Selint = 4.0D0 * CSelint
      ELSE
         x = 1.0D0/sin(phi)/sin(phi) - 1.0D0
         y = 1.0D0/sin(phi)/sin(phi) - k*k
         z = 1.0D0/sin(phi)/sin(phi)
         CALL Carlsonrf (rf, x, y, z)
         CALL Carlsonrd (rd, x, y, z)
         Felint = rf
         Selint = rf - 1.0D0/3.0D0*k*k*rd        
         
         IF ( phi .LE. pi/2.0D0 ) THEN
         ELSEIF ( phi > pi/2.0D0 .AND. phi .LE. pi) THEN
            Felint = 2.0D0*CFelint - Felint
            Selint = 2.0D0*CSelint - Selint
         ELSEIF ( phi > pi .AND. phi .LE. 3.0D0*pi/2.0D0 ) THEN
            Felint = 2.0D0*CFelint + Felint
            Selint = 2.0D0*CSelint + Selint
         ELSEIF ( phi > 3.0D0*pi/2.0D0 .AND. phi .LE. 2.0D0*pi) THEN
            Felint = 4.0D0*CFelint - Felint
            Selint = 4.0D0*CSelint - Selint
         ENDIF
         Felint = Felint + N2PI * 4.0D0*CFelint
         Selint = Selint + N2PI * 4.0D0*CSelint
         
         IF (PosNeg .EQ. -1) THEN
            Felint = -Felint
            Selint = -Selint
         ENDIF

      ENDIF

END SUBROUTINE EllipIntegral



SUBROUTINE CompleteEllipIntegral (Felint, Selint, k)

      DOUBLE PRECISION :: Felint, Selint, phi, k, rf, rd, x, y, z

      phi = ATAN(1.0D0)*2.0D0
      x = 1.0D0/sin(phi)/sin(phi) - 1.0D0
      y = 1.0D0/sin(phi)/sin(phi) - k*k
      z = 1.0D0/sin(phi)/sin(phi)

      CALL Carlsonrf (rf, x, y, z)
      CALL Carlsonrd (rd, x, y, z)

      Felint = rf
      Selint = rf - 1.0D0/3.0D0*k*k*rd

END SUBROUTINE CompleteEllipIntegral


SUBROUTINE TESTTHOMSON()

      DOUBLE PRECISION :: c, a, pi, mu, k, m, n, u1, v1, u2, v2, uu1, vv1, x1, z1, xx1, zz1, x2, z2, Felint, Selint, u
      DOUBLE PRECISION :: r1(3), r2(3), rr1(3), dr1, dr2
      DOUBLE PRECISION :: DE1, DE2, drdu(3)

!      c = 0.7D0
!      a = 0.1915D0
      pi = 4.0D0*ATAN(1.0D0)

!      mu = 2.0D0/(a+c)
!      k = SQRT(1-(a/c)**2)
!      m = (c**2-a**2)/2.0D0
!      n = (c**2+a**2)/2.0D0

      u = 0.50D0
      k = 0.15D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

      u = -0.40D0
      k = 0.20D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

      u = pi/2+0.30D0
      k = 0.25D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

      u = pi+0.20D0
      k = 0.30D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

      u = 2*pi-0.10D0
      k = 0.35D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

      u = -8*pi-0.90D0
      k = 0.40D0
      PRINT *, u, k
      CALL EllipIntegral(Felint, Selint, u, k)
      PRINT *, Felint, Selint

!      u1 = 0.5D0
!      v1 = 1.9D0
!      u2 = 1.8D0
!      v2 = 1.0D0
!      uu1 = 0.7001D0
!      vv1 = 1.900D0

!      CALL EllipIntegral(Felint, Selint, (mu*u1/2.0D0-pi/4.0D0), k)
!      x1 = a*Felint + c*Selint
!      z1 = SQRT(m*SIN(mu*u1)+n)

!      CALL EllipIntegral(Felint, Selint, (mu*uu1/2.0D0-pi/4.0D0), k)
!      xx1 = a*Felint + c*Selint
!      zz1 = SQRT(m*SIN(mu*uu1)+n)

!      CALL EllipIntegral(Felint, Selint, (mu*u2/2.0D0-pi/4.0D0), k)
!      x2 = a*Felint + c*Selint
!      z2 = SQRT(m*SIN(mu*u2)+n)

!      r1(1) = x1
!      r1(2) = z1*COS(v1)
!      r1(3) = z1*SIN(v1)

!      r2(1) = x2
!      r2(2) = z2*COS(v2)
!      r2(3) = z2*SIN(v2)

!      rr1(1) = xx1
!      rr1(2) = zz1*COS(vv1)
!      rr1(3) = zz1*SIN(vv1)
      
!      dr1 = SQRT(DOT_PRODUCT( r1-r2,  r1-r2))
!      dr2 = SQRT(DOT_PRODUCT(rr1-r2, rr1-r2))

!      DE1 = (1/dr2 - 1/dr1)/(uu1-u1)

!      drdu(1) = 1.0D0*(mu/2.0D0)*(a/SQRT(1-k**2*(SIN(mu*u1/2.0D0-pi/4.0D0))**2) + c*SQRT(1-k**2*(SIN(mu*u1/2.0D0-pi/4.0D0))**2))
!      drdu(2) = 0.5D0 / SQRT(m*SIN(mu*u1)+n) * m * mu * COS(mu*u1) * COS(v1)
!      drdu(3) = 0.5D0 / SQRT(m*SIN(mu*u1)+n) * m * mu * COS(mu*u1) * SIN(v1)
      
!      DE2 = -1.0D0 / (dr1**3) * DOT_PRODUCT(r1-r2,drdu)

!      PRINT *, r1
!      PRINT *, r2
!      PRINT *, rr1

!      PRINT *, dr1, dr2, DE1
!      PRINT *, dr1, dr2, DE2

!      PRINT *, (rr1-r1)/(uu1-u1)
!      PRINT *, drdu
!      PRINT *, xx1, x1
!      PRINT *, mu*u1/2-pi/4
!      CALL EllipIntegral(Felint, Selint, (mu*u1/2.0D0-pi/4.0D0), k)
!      PRINT *, Felint, Selint
!      PRINT *, k

END SUBROUTINE TESTTHOMSON
