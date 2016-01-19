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
C*************************************************************************
C
C  Here we calculate the Aluminium glue potential and gradient
C                                        
C*************************************************************************

      SUBROUTINE MGGLUE (X,V,EMG,GRADT)
      USE commons
      IMPLICIT NONE 
      INTEGER J1, J2
      DOUBLE PRECISION X(3*NATOMS), EMG, DIST, ETEMP, ETEMP2, 
     1                 V(3*NATOMS), rho1, RTEMP, dutemp(NATOMS),
     2                 drtemp, dvtemp, vtemp, vtemp1, vtemp2, vtemp3, 
     3                 rrc, rrcsq
      logical GRADT

      EMG=0.0D0
      RRC=6.67993511355129d0
      RRCSQ=RRC**2

      DO J1=1,NATOMS
        VT(J1)=0.0d0
      ENDDO

      ETEMP=0.0d0
      DO 22 J1=1,NATOMS
         RTEMP=0.0d0
         DO 23 J2=1,NATOMS
           if (j1.ne.j2) then
             DIST=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1             ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2             ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
             if (dist.lt.rrcsq) then
               dist=dsqrt(dist)
               call mgrh(dist,rho1,0)
               RTEMP=RTEMP+rho1
               if (j1.lt.j2) then
                 call mgv2(dist,Etemp2,0)
                 ETEMP=ETEMP+Etemp2
                 VT(J1)=VT(J1)+Etemp2/2.0d0
                 VT(J2)=VT(J2)+Etemp2/2.0d0
               endif
             endif
           endif
23       CONTINUE
         call mguu(RTEMP,Etemp2,0)
         EMG=EMG+ETEMP2
         VT(J1)=VT(J1)+Etemp2
         call mguu(RTEMP,dutemp(J1),1)
22    CONTINUE
      EMG=EMG+ETEMP

C
C Now calculate the gradient analytically.
C
      if (gradt) then 

      DO J1=1,NATOMS
         VTEMP1=0.0D0
         VTEMP2=0.0D0
         VTEMP3=0.0D0
         DO J2=1,NATOMS
           if (j1.ne.j2) then
             DIST=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1              ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2              ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
             if (dist.lt.rrcsq) then
               dist=dsqrt(dist)
               call mgv2(dist,dvtemp,1)
               call mgrh (dist,drtemp,1)
               VTEMP=(dvtemp+(dutemp(j1)+dutemp(j2))*drtemp)/dist
               VTEMP1=VTEMP1+VTEMP*(X(3*(J1-1)+1)-X(3*(J2-1)+1))
               VTEMP2=VTEMP2+VTEMP*(X(3*(J1-1)+2)-X(3*(J2-1)+2))
               VTEMP3=VTEMP3+VTEMP*(X(3*(J1-1)+3)-X(3*(J2-1)+3))
             endif
           endif
         ENDDO
         V(3*(J1-1)+1)=VTEMP1
         V(3*(J1-1)+2)=VTEMP2
         V(3*(J1-1)+3)=VTEMP3
      ENDDO

      endif

      RETURN
      END


*     Liu-Adams-Ercolessi-Moriarty glue potential for Mg.
*     Ref.: X.-Y. Liu, J. B. Adams, F. Ercolessi  and J. A. Moriarty,
*     Modelling Simul. Mater. Sci. Eng. 4, 293 (1996).
*     Potential home page: http://www.sissa.it/furio/potentials/Mg/

      SUBROUTINE mgv2(arg,func,dflag)
*        Magnesium : pair potential   and its first two derivatives.
*        Generated automatically by PoCo, version 26-nov-93
*        Hamiltonian type #  1, run on 95/03/23 at 14.05.29
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NV2
      parameter (nv2= 17)
      parameter (argmax=   .667993511355129D+01)
      integer dflag
      double precision xv2(nv2),yv2(nv2),bv2(nv2),cv2(nv2),dv2(nv2)
      save xv2,yv2,bv2,cv2,dv2
      data xv2(  1) /   .209371100573996D+01 /
      data xv2(  2) /   .249251310207138D+01 /
      data xv2(  3) /   .283149488395308D+01 /
      data xv2(  4) /   .317047666583479D+01 /
      data xv2(  5) /   .350945844771650D+01 /
      data xv2(  6) /   .384844022959820D+01 /
      data xv2(  7) /   .418742201147991D+01 /
      data xv2(  8) /   .452640379336162D+01 /
      data xv2(  9) /   .486538557524332D+01 /
      data xv2( 10) /   .520436735712503D+01 /
      data xv2( 11) /   .554334913900674D+01 /
      data xv2( 12) /   .588233092088845D+01 /
      data xv2( 13) /   .622131270277015D+01 /
      data xv2( 14) /   .656029448465186D+01 /
      data xv2( 15) /   .667993511355129D+01 /
      data xv2( 16) /   .667996901172947D+01 /
      data xv2( 17) /   .668000290990766D+01 /
      data yv2(  1) /   .153928517802236D+01 /
      data yv2(  2) /   .302534885675981D+00 /
      data yv2(  3) /  -.475831443182775D-01 /
      data yv2(  4) /  -.126245209044009D+00 /
      data yv2(  5) /  -.113018124801363D+00 /
      data yv2(  6) /  -.793483963982679D-01 /
      data yv2(  7) /  -.547284563621586D-01 /
      data yv2(  8) /  -.318233872544895D-01 /
      data yv2(  9) /  -.971224623229114D-02 /
      data yv2( 10) /   .151183119534004D-02 /
      data yv2( 11) /   .146823073635825D-02 /
      data yv2( 12) /  -.188223305893613D-02 /
      data yv2( 13) /  -.354784297292924D-02 /
      data yv2( 14) /  -.110924151560009D-02 /
      data yv2( 15) /   .000000000000000D+00 /
      data yv2( 16) /   .000000000000000D+00 /
      data yv2( 17) /   .000000000000000D+00 /
      data bv2(  1) /  -.467533787577523D+01 /
      data bv2(  2) /  -.176659850116244D+01 /
      data bv2(  3) /  -.496070082140902D+00 /
      data bv2(  4) /  -.438392138791983D-01 /
      data bv2(  5) /   .923251925445662D-01 /
      data bv2(  6) /   .895767892906107D-01 /
      data bv2(  7) /   .652331470961059D-01 /
      data bv2(  8) /   .700885277341982D-01 /
      data bv2(  9) /   .528076929875946D-01 /
      data bv2( 10) /   .136984879804464D-01 /
      data bv2( 11) /  -.865405333157377D-02 /
      data bv2( 12) /  -.911985574864718D-02 /
      data bv2( 13) /   .741057818215175D-03 /
      data bv2( 14) /   .129966252432477D-01 /
      data bv2( 15) /   .251830318859266D-05 /
      data bv2( 16) /  -.839682282456011D-06 /
      data bv2( 17) /   .840425941231385D-06 /
      data cv2(  1) /   .454808379415846D+01 /
      data cv2(  2) /   .274560749124421D+01 /
      data cv2(  3) /   .100246537558882D+01 /
      data cv2(  4) /   .331620679866272D+00 /
      data cv2(  5) /   .700658227767553D-01 /
      data cv2(  6) /  -.781736427271344D-01 /
      data cv2(  7) /   .635962942125154D-02 /
      data cv2(  8) /   .796379708013682D-02 /
      data cv2(  9) /  -.589424504186043D-01 /
      data cv2( 10) /  -.564301362438733D-01 /
      data cv2( 11) /  -.951010747052316D-02 /
      data cv2( 12) /   .813598519619999D-02 /
      data cv2( 13) /   .209538187222403D-01 /
      data cv2( 14) /   .152002405224814D-01 /
      data cv2( 15) /  -.123809724270535D+00 /
      data cv2( 16) /   .247487820445084D-01 /
      data cv2( 17) /   .248145960925009D-01 /
      data dv2(  1) /  -.150657541630735D+01 /
      data dv2(  2) /  -.171409616369264D+01 /
      data dv2(  3) /  -.659666420221412D+00 /
      data dv2(  4) /  -.257196572273601D+00 /
      data dv2(  5) /  -.145769353425626D+00 /
      data dv2(  6) /   .831246954523014D-01 /
      data dv2(  7) /   .157743743629373D-02 /
      data dv2(  8) /  -.657913896604715D-01 /
      data dv2(  9) /   .247045151982798D-02 /
      data dv2( 10) /   .461382010097950D-01 /
      data dv2( 11) /   .173520560788141D-01 /
      data dv2( 12) /   .126042501506401D-01 /
      data dv2( 13) /  -.565770640909810D-02 /
      data dv2( 14) /  -.387298657297743D+00 /
      data dv2( 15) /   .146083078065090D+04 /
      data dv2( 16) /   .647173894591949D+00 /
      data dv2( 17) /   .647173894591949D+00 /
      if (arg.ge.argmax) then
         func =   .000000000000000D+00
      else
         call seval(nv2,arg,xv2,yv2,bv2,cv2,dv2,func,dflag)
      endif
      end

      SUBROUTINE mgrh(arg,func,dflag)
*        Magnesium : atomic density   and its first two derivatives.
*        Generated automatically by PoCo, version 26-nov-93
*        Hamiltonian type #  1, run on 95/03/23 at 14.05.29
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NRH
      parameter (nrh= 10)
      parameter (argmax=   .667993511355129D+01)
      integer dflag
      double precision xrh(nrh),yrh(nrh),brh(nrh),crh(nrh),drh(nrh)
      save xrh,yrh,brh,crh,drh
      data xrh(  1) /   .209371100573996D+01 /
      data xrh(  2) /   .249251310207138D+01 /
      data xrh(  3) /   .330606937858747D+01 /
      data xrh(  4) /   .411962565510357D+01 /
      data xrh(  5) /   .493318193161967D+01 /
      data xrh(  6) /   .574673820813576D+01 /
      data xrh(  7) /   .656029448465186D+01 /
      data xrh(  8) /   .667993511355129D+01 /
      data xrh(  9) /   .668001646917894D+01 /
      data xrh( 10) /   .668009782480659D+01 /
      data yrh(  1) /   .143187372690200D+00 /
      data yrh(  2) /   .891051954501756D-01 /
      data yrh(  3) /   .497245267820974D-01 /
      data yrh(  4) /   .346372124961649D-01 /
      data yrh(  5) /   .235194366080153D-01 /
      data yrh(  6) /   .302069379049224D-02 /
      data yrh(  7) /   .802148619196398D-03 /
      data yrh(  8) /   .000000000000000D+00 /
      data yrh(  9) /   .000000000000000D+00 /
      data yrh( 10) /   .000000000000000D+00 /
      data brh(  1) /  -.178789353139449D+00 /
      data brh(  2) /  -.966399582905674D-01 /
      data brh(  3) /  -.228583230701083D-01 /
      data brh(  4) /  -.127781884490251D-01 /
      data brh(  5) /  -.226605548863010D-01 /
      data brh(  6) /  -.131659411679867D-01 /
      data brh(  7) /  -.844599349041348D-02 /
      data brh(  8) /  -.475726389590590D-05 /
      data brh(  9) /   .158678593730407D-05 /
      data brh( 10) /  -.158987985331039D-05 /
      data crh(  1) /   .118815732796611D+00 /
      data crh(  2) /   .871746459014311D-01 /
      data crh(  3) /   .351562021809925D-02 /
      data crh(  4) /   .887459163577463D-02 /
      data crh(  5) /  -.210217125203868D-01 /
      data crh(  6) /   .326922189037599D-01 /
      data crh(  7) /  -.268905948319926D-01 /
      data crh(  8) /   .974455250774157D-01 /
      data crh(  9) /  -.194662872995196D-01 /
      data crh( 10) /  -.195803758793374D-01 /
      data drh(  1) /  -.264467741320369D-01 /
      data drh(  2) /  -.342770902264586D-01 /
      data drh(  3) /   .219569789755835D-02 /
      data drh(  4) /  -.122492260294464D-01 /
      data drh(  5) /   .220078737327456D-01 /
      data drh(  6) /  -.244124942369915D-01 /
      data drh(  7) /   .346415541980947D+00 /
      data drh(  8) /  -.479015468881384D+03 /
      data drh(  9) /  -.467448014399572D+00 /
      data drh( 10) /  -.467448014399572D+00 /
      if (arg.ge.argmax) then
         func =   .000000000000000D+00
      else
         call seval(nrh,arg,xrh,yrh,brh,crh,drh,func,dflag)
      endif
      end

      SUBROUTINE mguu(arg,func,dflag)
*        Magnesium : glue function    and its first two derivatives.
*        Generated automatically by PoCo, version 26-nov-93
*        Hamiltonian type #  1, run on 95/03/23 at 14.05.29
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NUU
      parameter (nuu= 11)
      parameter (argmin=   .000000000000000D+00)
      integer dflag
      double precision xuu(nuu),yuu(nuu),buu(nuu),cuu(nuu),duu(nuu)
      save xuu,yuu,buu,cuu,duu
      data xuu(  1) /   .000000000000000D+00 /
      data xuu(  2) /   .100000000000000D+00 /
      data xuu(  3) /   .237500000000000D+00 /
      data xuu(  4) /   .375000000000000D+00 /
      data xuu(  5) /   .512500000000000D+00 /
      data xuu(  6) /   .650000000000000D+00 /
      data xuu(  7) /   .787500000000000D+00 /
      data xuu(  8) /   .925000000000000D+00 /
      data xuu(  9) /   .106250000000000D+01 /
      data xuu( 10) /   .120000000000000D+01 /
      data xuu( 11) /   .140000000000000D+01 /
      data yuu(  1) /   .000000000000000D+00 /
      data yuu(  2) /  -.206805227959338D+00 /
      data yuu(  3) /  -.290481511485127D+00 /
      data yuu(  4) /  -.403009316082408D+00 /
      data yuu(  5) /  -.556469868697284D+00 /
      data yuu(  6) /  -.609101336623323D+00 /
      data yuu(  7) /  -.626377489625343D+00 /
      data yuu(  8) /  -.635377221699801D+00 /
      data yuu(  9) /  -.635910883069410D+00 /
      data yuu( 10) /  -.627564569808675D+00 /
      data yuu( 11) /  -.565414923439952D+00 /
      data buu(  1) /  -.315967367115533D+01 /
      data buu(  2) /  -.116065136867398D+01 /
      data buu(  3) /  -.498734722118105D+00 /
      data buu(  4) /  -.112522621099331D+01 /
      data buu(  5) /  -.803742773083007D+00 /
      data buu(  6) /  -.156355872112811D+00 /
      data buu(  7) /  -.960909223506751D-01 /
      data buu(  8) /  -.325724765167350D-01 /
      data buu(  9) /   .183795169106987D-01 /
      data buu( 10) /   .129512268316687D+00 /
      data buu( 11) /   .534598088628589D+00 /
      data cuu(  1) /   .127584187220449D+02 /
      data cuu(  2) /   .723180430276862D+01 /
      data cuu(  3) /  -.241786505508954D+01 /
      data cuu(  4) /  -.213843668218466D+01 /
      data cuu(  5) /   .447649804880503D+01 /
      data cuu(  6) /   .231770321887313D+00 /
      data cuu(  7) /   .206520221837309D+00 /
      data cuu(  8) /   .255432111500437D+00 /
      data cuu(  9) /   .115127840699081D+00 /
      data cuu( 10) /   .693110351344470D+00 /
      data cuu( 11) /   .133231875021504D+01 /
      data duu(  1) /  -.184220480642544D+02 /
      data duu(  2) /  -.233931378372319D+02 /
      data duu(  3) /   .677402116133043D+00 /
      data duu(  4) /   .160362054084599D+02 /
      data duu(  5) /  -.102902490349520D+02 /
      data duu(  6) /  -.612123637575849D-01 /
      data duu(  7) /   .118574277971218D+00 /
      data duu(  8) /  -.340131565579043D+00 /
      data duu(  9) /   .140116972277670D+01 /
      data duu( 10) /   .106534733145095D+01 /
      data duu( 11) /   .106534733145095D+01 /
      if (arg.le.argmin) then
         func =   .000000000000000D+00
      else
         call seval(nuu,arg,xuu,yuu,buu,cuu,duu,func,dflag)
      endif
      end

