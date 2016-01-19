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

      SUBROUTINE ALGLUE (X,V,EAL,GRADT)
      USE commons
      IMPLICIT NONE 
      INTEGER J1, J2
      DOUBLE PRECISION X(3*NATOMS), EAL, DIST, ETEMP, ETEMP2, 
     1                 V(3*NATOMS), rho1, RTEMP, dutemp(NATOMS),
     2                 drtemp, dvtemp, vtemp, vtemp1, vtemp2, vtemp3, 
     3                 rrc, rrcsq
      logical GRADT

      EAL=0.0D0
      RRC=5.55805441821810D0
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
               call alrh(dist,rho1,0)
               RTEMP=RTEMP+rho1
               if (j1.lt.j2) then
                 call alv2(dist,Etemp2,0)
                 ETEMP=ETEMP+Etemp2
                 VT(J1)=VT(J1)+Etemp2/2.0d0
                 VT(J2)=VT(J2)+Etemp2/2.0d0
               endif
             endif
           endif
23       CONTINUE
         call aluu(RTEMP,Etemp2,0)
         EAL=EAL+ETEMP2
         VT(J1)=VT(J1)+Etemp2
         call aluu(RTEMP,dutemp(J1),1)
22    CONTINUE
      EAL=EAL+ETEMP

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
               call alv2(dist,dvtemp,1)
               call alrh (dist,drtemp,1)
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

*     Ercolessi-Adams glue potential for Al.
*     Ref.: F. Ercolessi and J. B. Adams, Europhys. Lett. 26, 583 (1994).
*     Potential home page: http://www.sissa.it/furio/potentials/Al/

      SUBROUTINE alv2(arg,func,dflag)
*        Aluminum  : pair potential   and its first two derivatives.
*        Generated automatically by PoCo, version 04-may-93           
*        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NV2
      parameter (nv2= 17)
      parameter (argmax=   .555805441821810D+01)
      integer dflag
      double precision xv2(nv2),yv2(nv2),bv2(nv2),cv2(nv2),dv2(nv2)
      save xv2,yv2,bv2,cv2,dv2
      data xv2(  1) /   .202111069753385D+01 /
      data xv2(  2) /   .227374953472558D+01 /
      data xv2(  3) /   .252638837191732D+01 /
      data xv2(  4) /   .277902720910905D+01 /
      data xv2(  5) /   .303166604630078D+01 /
      data xv2(  6) /   .328430488349251D+01 /
      data xv2(  7) /   .353694372068424D+01 /
      data xv2(  8) /   .378958255787597D+01 /
      data xv2(  9) /   .404222139506771D+01 /
      data xv2( 10) /   .429486023225944D+01 /
      data xv2( 11) /   .454749906945117D+01 /
      data xv2( 12) /   .480013790664290D+01 /
      data xv2( 13) /   .505277674383463D+01 /
      data xv2( 14) /   .530541558102636D+01 /
      data xv2( 15) /   .555805441821810D+01 /
      data xv2( 16) /   .555807968210182D+01 /
      data xv2( 17) /   .555810494598553D+01 /
      data yv2(  1) /   .196016472197158D+01 /
      data yv2(  2) /   .682724240745344D+00 /
      data yv2(  3) /   .147370824539188D+00 /
      data yv2(  4) /  -.188188235860390D-01 /
      data yv2(  5) /  -.576011902692490D-01 /
      data yv2(  6) /  -.519846499644276D-01 /
      data yv2(  7) /  -.376352484845919D-01 /
      data yv2(  8) /  -.373737879689433D-01 /
      data yv2(  9) /  -.531351030124350D-01 /
      data yv2( 10) /  -.632864983555742D-01 /
      data yv2( 11) /  -.548103623840369D-01 /
      data yv2( 12) /  -.372889232343935D-01 /
      data yv2( 13) /  -.188876517630154D-01 /
      data yv2( 14) /  -.585239362533525D-02 /
      data yv2( 15) /   .000000000000000D+00 /
      data yv2( 16) /   .000000000000000D+00 /
      data yv2( 17) /   .000000000000000D+00 /
      data bv2(  1) /  -.702739315585347D+01 /
      data bv2(  2) /  -.333140549270729D+01 /
      data bv2(  3) /  -.117329394261502D+01 /
      data bv2(  4) /  -.306003283486901D+00 /
      data bv2(  5) /  -.366656699104026D-01 /
      data bv2(  6) /   .588330899204400D-01 /
      data bv2(  7) /   .384220572312032D-01 /
      data bv2(  8) /  -.390223173707191D-01 /
      data bv2(  9) /  -.663882722510521D-01 /
      data bv2( 10) /  -.312918894386669D-02 /
      data bv2( 11) /   .590118945294245D-01 /
      data bv2( 12) /   .757939459148246D-01 /
      data bv2( 13) /   .643822548468606D-01 /
      data bv2( 14) /   .399750987463792D-01 /
      data bv2( 15) /   .177103852679117D-05 /
      data bv2( 16) /  -.590423369301474D-06 /
      data bv2( 17) /   .590654950414731D-06 /
      data cv2(  1) /   .877545959718548D+01 /
      data cv2(  2) /   .585407125495837D+01 /
      data cv2(  3) /   .268820820643116D+01 /
      data cv2(  4) /   .744718689404422D+00 /
      data cv2(  5) /   .321378734769888D+00 /
      data cv2(  6) /   .566263292669091D-01 /
      data cv2(  7) /  -.137417679148505D+00 /
      data cv2(  8) /  -.169124163201523D+00 /
      data cv2(  9) /   .608037039066423D-01 /
      data cv2( 10) /   .189589640245655D+00 /
      data cv2( 11) /   .563784150384640D-01 /
      data cv2( 12) /   .100486298765028D-01 /
      data cv2( 13) /  -.552186092621482D-01 /
      data cv2( 14) /  -.413902746758285D-01 /
      data cv2( 15) /  -.116832934994489D+00 /
      data cv2( 16) /   .233610871054729D-01 /
      data cv2( 17) /   .233885865725971D-01 /
      data dv2(  1) /  -.385449887634130D+01 /
      data dv2(  2) /  -.417706040200591D+01 /
      data dv2(  3) /  -.256425277368288D+01 /
      data dv2(  4) /  -.558557503589276D+00 /
      data dv2(  5) /  -.349316054551627D+00 /
      data dv2(  6) /  -.256022933201611D+00 /
      data dv2(  7) /  -.418337423301704D-01 /
      data dv2(  8) /   .303368330939646D+00 /
      data dv2(  9) /   .169921006301015D+00 /
      data dv2( 10) /  -.175759761362548D+00 /
      data dv2( 11) /  -.611278214082881D-01 /
      data dv2( 12) /  -.861140219824535D-01 /
      data dv2( 13) /   .182451950513387D-01 /
      data dv2( 14) /  -.995395392057973D-01 /
      data dv2( 15) /   .184972909229936D+04 /
      data dv2( 16) /   .362829766922787D+00 /
      data dv2( 17) /   .362829766922787D+00 /

      if (arg.ge.argmax) then
         func =   .000000000000000D+00
      else 
         call seval(nv2,arg,xv2,yv2,bv2,cv2,dv2,func,dflag)
      endif

      end

      SUBROUTINE alrh(arg,func,dflag)
*        Aluminum  : atomic density   and its first two derivatives.
*        Generated automatically by PoCo, version 04-may-93           
*        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NRH
      parameter (nrh= 17)
      parameter (argmax=   .555805441821810D+01)
      integer dflag
      double precision xrh(nrh),yrh(nrh),brh(nrh),crh(nrh),drh(nrh)
      save xrh,yrh,brh,crh,drh
      data xrh(  1) /   .202111069753385D+01 /
      data xrh(  2) /   .227374953472558D+01 /
      data xrh(  3) /   .252638837191732D+01 /
      data xrh(  4) /   .277902720910905D+01 /
      data xrh(  5) /   .303166604630078D+01 /
      data xrh(  6) /   .328430488349251D+01 /
      data xrh(  7) /   .353694372068424D+01 /
      data xrh(  8) /   .378958255787597D+01 /
      data xrh(  9) /   .404222139506771D+01 /
      data xrh( 10) /   .429486023225944D+01 /
      data xrh( 11) /   .454749906945117D+01 /
      data xrh( 12) /   .480013790664290D+01 /
      data xrh( 13) /   .505277674383463D+01 /
      data xrh( 14) /   .530541558102636D+01 /
      data xrh( 15) /   .555805441821810D+01 /
      data xrh( 16) /   .555807968210182D+01 /
      data xrh( 17) /   .555810494598553D+01 /
      data yrh(  1) /   .865674623712589D-01 /
      data yrh(  2) /   .925214702944478D-01 /
      data yrh(  3) /   .862003123832002D-01 /
      data yrh(  4) /   .762736292751052D-01 /
      data yrh(  5) /   .606481841271735D-01 /
      data yrh(  6) /   .466030959588197D-01 /
      data yrh(  7) /   .338740138848363D-01 /
      data yrh(  8) /   .232572661705343D-01 /
      data yrh(  9) /   .109046405489829D-01 /
      data yrh( 10) /   .524910605677597D-02 /
      data yrh( 11) /   .391702419142291D-02 /
      data yrh( 12) /   .308277776293383D-02 /
      data yrh( 13) /   .250214745349505D-02 /
      data yrh( 14) /   .147220513798186D-02 /
      data yrh( 15) /   .000000000000000D+00 /
      data yrh( 16) /   .000000000000000D+00 /
      data yrh( 17) /   .000000000000000D+00 /
      data brh(  1) /   .608555214104682D-01 /
      data brh(  2) /  -.800158928716306D-02 /
      data brh(  3) /  -.332089451111092D-01 /
      data brh(  4) /  -.521001991705069D-01 /
      data brh(  5) /  -.618130637429111D-01 /
      data brh(  6) /  -.529750064268036D-01 /
      data brh(  7) /  -.442210477548108D-01 /
      data brh(  8) /  -.473645664984640D-01 /
      data brh(  9) /  -.390741582571631D-01 /
      data brh( 10) /  -.101795580610560D-01 /
      data brh( 11) /  -.318316981110289D-02 /
      data brh( 12) /  -.281217210746153D-02 /
      data brh( 13) /  -.236932031483360D-02 /
      data brh( 14) /  -.683554708271547D-02 /
      data brh( 15) /  -.638718204858808D-06 /
      data brh( 16) /   .212925486831149D-06 /
      data brh( 17) /  -.212983742465787D-06 /
      data crh(  1) /  -.170233687052940D+00 /
      data crh(  2) /  -.102317878901959D+00 /
      data crh(  3) /   .254162872544396D-02 /
      data crh(  4) /  -.773173610292656D-01 /
      data crh(  5) /   .388717099948882D-01 /
      data crh(  6) /  -.388873819867093D-02 /
      data crh(  7) /   .385388290924526D-01 /
      data crh(  8) /  -.509815666327127D-01 /
      data crh(  9) /   .837968231208082D-01 /
      data crh( 10) /   .305743500420042D-01 /
      data crh( 11) /  -.288110886134041D-02 /
      data crh( 12) /   .434959924771674D-02 /
      data crh( 13) /  -.259669459714693D-02 /
      data crh( 14) /  -.150816117849093D-01 /
      data crh( 15) /   .421356801161513D-01 /
      data crh( 16) /  -.842575249165724D-02 /
      data crh( 17) /  -.843267014952237D-02 /
      data drh(  1) /   .896085612514625D-01 /
      data drh(  2) /   .138352319847830D+00 /
      data drh(  3) /  -.105366473134009D+00 /
      data drh(  4) /   .153300619856764D+00 /
      data drh(  5) /  -.564184148788224D-01 /
      data drh(  6) /   .559792096400504D-01 /
      data drh(  7) /  -.118113795329664D+00 /
      data drh(  8) /   .177827488509794D+00 /
      data drh(  9) /  -.702220789044304D-01 /
      data drh( 10) /  -.441413511810337D-01 /
      data drh( 11) /   .954024354744484D-02 /
      data drh( 12) /  -.916498550800407D-02 /
      data drh( 13) /  -.164726813535368D-01 /
      data drh( 14) /   .754928689733184D-01 /
      data drh( 15) /  -.667110847110954D+03 /
      data drh( 16) /  -.912720300911022D-01 /
      data drh( 17) /  -.912720300911022D-01 /
      if (arg.ge.argmax) then
         func =   .000000000000000D+00
      else
         call seval(nrh,arg,xrh,yrh,brh,crh,drh,func,dflag)
      endif
      end

      SUBROUTINE aluu(arg,func,dflag)
*        Aluminum  : glue function    and its first two derivatives.
*        Generated automatically by PoCo, version 04-may-93           
*        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
*        Uses subroutine seval from netlib@ornl.gov [to get it,
*        use 'send seval from sfmm'], trivially modified to
*        compute also dfunc and d2func and use double precision.
      implicit double precision (a-h,o-z)
      INTEGER NUU
      parameter (nuu= 13)
      parameter (argmin=   .000000000000000D+00)
      integer dflag
      double precision xuu(nuu),yuu(nuu),buu(nuu),cuu(nuu),duu(nuu)
      save xuu,yuu,buu,cuu,duu
      data xuu(  1) /   .000000000000000D+00 /
      data xuu(  2) /   .100000000000000D+00 /
      data xuu(  3) /   .200000000000000D+00 /
      data xuu(  4) /   .300000000000000D+00 /
      data xuu(  5) /   .400000000000000D+00 /
      data xuu(  6) /   .500000000000000D+00 /
      data xuu(  7) /   .600000000000000D+00 /
      data xuu(  8) /   .700000000000000D+00 /
      data xuu(  9) /   .800000000000000D+00 /
      data xuu( 10) /   .900000000000000D+00 /
      data xuu( 11) /   .100000000000000D+01 /
      data xuu( 12) /   .110000000000000D+01 /
      data xuu( 13) /   .120000000000000D+01 /
      data yuu(  1) /   .000000000000000D+00 /
      data yuu(  2) /  -.113953324143752D+01 /
      data yuu(  3) /  -.145709859805864D+01 /
      data yuu(  4) /  -.174913308002738D+01 /
      data yuu(  5) /  -.202960322136630D+01 /
      data yuu(  6) /  -.225202324967546D+01 /
      data yuu(  7) /  -.242723053979436D+01 /
      data yuu(  8) /  -.255171976467357D+01 /
      data yuu(  9) /  -.260521638832322D+01 /
      data yuu( 10) /  -.264397894381693D+01 /
      data yuu( 11) /  -.265707884842034D+01 /
      data yuu( 12) /  -.264564149400021D+01 /
      data yuu( 13) /  -.260870604452106D+01 /
      data buu(  1) /  -.183757286015853D+02 /
      data buu(  2) /  -.574233124410516D+01 /
      data buu(  3) /  -.236790436375322D+01 /
      data buu(  4) /  -.307404645857774D+01 /
      data buu(  5) /  -.251104850116555D+01 /
      data buu(  6) /  -.196846462620234D+01 /
      data buu(  7) /  -.154391254686695D+01 /
      data buu(  8) /  -.846780636273251D+00 /
      data buu(  9) /  -.408540363905760D+00 /
      data buu( 10) /  -.286833282404628D+00 /
      data buu( 11) /  -.309389414590161D-06 /
      data buu( 12) /   .236958014464143D+00 /
      data buu( 13) /   .503352368511243D+00 /
      data cuu(  1) /   .830779120415016D+02 /
      data cuu(  2) /   .432560615333001D+02 /
      data cuu(  3) /  -.951179272978074D+01 /
      data cuu(  4) /   .245037178153561D+01 /
      data cuu(  5) /   .317960779258630D+01 /
      data cuu(  6) /   .224623095704576D+01 /
      data cuu(  7) /   .199928983630817D+01 /
      data cuu(  8) /   .497202926962879D+01 /
      data cuu(  9) /  -.589626545953876D+00 /
      data cuu( 10) /   .180669736096520D+01 /
      data cuu( 11) /   .106163236918694D+01 /
      data cuu( 12) /   .130795086934864D+01 /
      data cuu( 13) /   .135599267112235D+01 /
      data duu(  1) /  -.132739501694005D+03 /
      data duu(  2) /  -.175892847543603D+03 /
      data duu(  3) /   .398738817043878D+02 /
      data duu(  4) /   .243078670350231D+01 /
      data duu(  5) /  -.311125611846847D+01 /
      data duu(  6) /  -.823137069125319D+00 /
      data duu(  7) /   .990913144440207D+01 /
      data duu(  8) /  -.185388527186089D+02 /
      data duu(  9) /   .798774635639692D+01 /
      data duu( 10) /  -.248354997259420D+01 /
      data duu( 11) /   .821061667205675D+00 /
      data duu( 12) /   .160139339245701D+00 /
      data duu( 13) /   .160139339245701D+00 /
      if (arg.le.argmin) then
         func =   .000000000000000D+00
         dfunc  = 0.d0
         d2func = 0.d0
      else
         call seval(nuu,arg,xuu,yuu,buu,cuu,duu,func,dflag)
      endif
      end

      subroutine seval(n, u, x, y, b, c, d, f, dflag)
      implicit none
      integer n, dflag
      double precision  u, x(n), y(n), b(n), c(n), d(n)
      double precision f
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      if (dflag.eq.0) then
        f = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      else if (dflag.eq.1) then
        f = b(i) + dx*(2.d0*c(i) + 3.d0*dx*d(i))
      else
        print*,'WARNING: inappropriate value of dflag'
        stop
      endif

      return
      end
