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
C       POTENTIAL ENERGY OF A (C60)N SYSTEM
C       PACHECO-PRATES-RAMALHO POTENTIAL (PRL 1997)
C
C      X,Y,Z: vectors

      SUBROUTINE PRC60(LCOORDS,V,EPPR,GTEST)
      USE commons
      IMPLICIT NONE
C       PARAMETER(CAT=4752000.D0)
        INTEGER J1, J2, I, J
      DOUBLE PRECISION x(NATOMS),y(NATOMS),z(NATOMS),LCOORDS(3*NATOMS),V2B,EPPR,V(3*NATOMS)
      DOUBLE PRECISION mij,fij,wij,xij,yij,zij,rij,rij2,vij,pmorse,wtot,rij23,
     1                   fx(NATOMS),fy(NATOMS),fZ(NATOMS),ermi,dwij,dfij,dmij
        double precision dmu, delta
      parameter(dmu=10.05D0)
      parameter(delta=1.04D0)
        double precision dM0, tau, D0
      parameter(dM0=0.3D0)
      parameter(tau=9.75D0)
      parameter(D0=10.3D0)
        double precision C6, C8, C10, C12
      parameter(C6=75600.D0)
      parameter(C8=9122400.D0)
      parameter(C10=2.09D8)
      parameter(C12=7.78D10)
        LOGICAL GTEST

        DO J1=1,NATOMS
           VT(J1)=0.0D0
           J2=3*(J1-1)
           X(J1)=LCOORDS(J2+1)
           Y(J1)=LCOORDS(J2+2)
           Z(J1)=LCOORDS(J2+3)
        ENDDO

      V2B=0.D0

        IF (.NOT.GTEST) THEN
      do i=1,NATOMS
         do j=i+1,NATOMS
C
C       2-body interaction
C
            xij=x(j)-x(i)
            yij=y(j)-y(i)
            zij=z(j)-z(i)
            rij2=xij*xij+yij*yij+zij*zij
            rij=sqrt(rij2)
            
            fij=1.0D0/(1.D0+exp((rij-dmu)/delta))
              pmorse=exp(tau*(1.D0-rij/d0))
            mij=dM0*pmorse*(pmorse-2.D0)
            wij=-(C6+(C8+(C10+C12/rij2)/rij2)/rij2)/rij2**3
            
            vij=fij*mij+(1.D0-fij)*wij
            
            v2b=v2b+vij
              VT(i)=VT(i)+vij
              VT(j)=VT(j)+vij
         enddo
      enddo

        ELSE

      do i=1,NATOMS
           fx(i)=0.D0
           fy(i)=0.D0
           fz(i)=0.D0
        enddo
      do i=1,NATOMS
         do j=i+1,NATOMS
C
C       2-body interaction
C
            xij=x(j)-x(i)
            yij=y(j)-y(i)
            zij=z(j)-z(i)
            rij2=xij*xij+yij*yij+zij*zij
            rij=sqrt(rij2)
            
              ermi=exp((rij-dmu)/delta)
              fij=1.0D0/(1.D0+ermi)
              dfij=-ermi/(delta*(1.D0+ermi)**2)

              pmorse=exp(tau*(1.D0-rij/d0))
            mij=dM0*pmorse*(pmorse-2.D0)
              dmij=(2.D0*tau*dM0*pmorse*(1.D0-pmorse))/d0

              rij23=rij2**3
            wij=-(C6+(C8+(C10+C12/rij2)/rij2)/rij2)/rij23
              dwij=(6*C6+(8*C8+(10*C10+12*C12/rij2)/rij2)/rij2)/(rij*rij23)
            
            vij=fij*mij+(1.D0-fij)*wij
            
            v2b=v2b+vij
              VT(i)=VT(i)+vij
              VT(j)=VT(j)+vij

              wtot=mij*dfij+fij*dmij+(1.D0-fij)*dwij-dfij*wij
              wtot=-wtot

              fx(i)=fx(i)+(x(i)-x(j))*wtot/rij
              fy(i)=fy(i)+(y(i)-y(j))*wtot/rij
              fz(i)=fz(i)+(z(i)-z(j))*wtot/rij

              fx(j)=fx(j)+(x(j)-x(i))*wtot/rij
              fy(j)=fy(j)+(y(j)-y(i))*wtot/rij
              fz(j)=fz(j)+(z(j)-z(i))*wtot/rij
       
         enddo
      enddo
        ENDIF
      
      EPPR=v2b

        IF (GTEST) THEN
           DO J1=1,NATOMS
              J2=3*(J1-1)
              V(J2+1)=-FX(J1)
              V(J2+2)=-FY(J1)
              V(J2+3)=-FZ(J1)
C             WRITE(*,'(A,I3,3F20.10)') 'J1,FX,FY,FZ=',J1,FX(J1),FY(J1),FZ(J1)
           ENDDO
        ENDIF
      
      return
      end
