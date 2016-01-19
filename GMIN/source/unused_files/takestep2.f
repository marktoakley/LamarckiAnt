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
      SUBROUTINE TAKESTEP2(NP)
      USE COMMONS
      IMPLICIT NONE

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS), RDOTN, XL, YL, ZL,
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS), XC, YC, ZC, ANGLE, COST, SINT,
     2                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ, RVEC(3), TX, TY, TZ
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, NP, J3, JMAX2, RANATOM, INDEXD(NATOMS)
C
C  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
C  outside the permitted radius. Then it may be impossible to take a step that keeps all the
C  atoms inside.
C
      IF ((.NOT.NORESET).AND.(.NOT.PERMUTE).AND.(.NOT.DIFFRACTT).AND.(.NOT.BLNT).AND.(.NOT.PERIODIC)
     &     .AND.(.NOT.GAUSST)) THEN
         DO J1=1,NATOMS
            IF ((.NOT.RIGID).OR.(J1.LE.NATOMS/2)) THEN
               J2=3*J1
               DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
               IF (DUMMY2.GT.RADIUS) THEN
                  WRITE(MYUNIT,'(A,I5,5F20.10)') 'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                  WRITE(MYUNIT,'(A)') 'initial coordinate outside container - increase container radius'
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      DO J1=1,3*(NATOMS-NSEED)
         COORDSO(J1,NP)=COORDS(J1,NP)
      ENDDO
      DO J1=1,NATOMS
         VATO(J1,NP)=VAT(J1,NP)
      ENDDO
      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS
C
C  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
C  and the pair energy of the most tightly bound atom, VMIN. An angular step is
C  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
C  DMAX (or CMMAX from CM of the cluster).
C
      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,NP).GT.VMAX) THEN
            VMAX=VAT(J1,NP)
            JMAX=J1
         ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
            VMAX2=VAT(J1,NP)
            JMAX2=J1
         ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO
C
C  SHELLMOVES is true: try an angular move for all the free atoms 
C  about the centre of coordinates of the frozen set.
C
      NSURFMOVES=NSURFMOVES+1
      IF (DPRAND().GT.(1.0D0-SHELLPROB)) THEN 
         IF (NCORE.EQ.0) THEN
            WRITE(MYUNIT,'(A)') 'takestep> ERROR - NCORE=0 - turning off surface moves'
            SHELLMOVES=.FALSE.
         ENDIF
         WRITE(MYUNIT,'(A,I8)') 'takestep> shell move number ',NSURFMOVES
         XC=0.0D0; YC=0.0D0; ZC=0.0D0
         DO J1=1,NCORE
            XC=XC+COORDS(3*(J1-1)+1,NP)
            YC=YC+COORDS(3*(J1-1)+2,NP)
            ZC=ZC+COORDS(3*(J1-1)+3,NP)
         ENDDO
         XC=XC/NCORE; YC=YC/NCORE; ZC=ZC/NCORE
         IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'takestep> centre of coordinates for frozen atoms: ',XC, YC, ZC

         RVEC(1)=(DPRAND()-0.5D0)*2.0D0
         RVEC(2)=(DPRAND()-0.5D0)*2.0D0
         RVEC(3)=(DPRAND()-0.5D0)*2.0D0
         DUMMY=SQRT(RVEC(1)**2+RVEC(2)**2+RVEC(3)**2)
         RVEC(1)=RVEC(1)/DUMMY; RVEC(2)=RVEC(2)/DUMMY; RVEC(3)=RVEC(3)/DUMMY
         ANGLE=DPRAND()*PI*2.0D0
         COST=COS(ANGLE)
         SINT=SIN(ANGLE)

         WRITE(MYUNIT,'(A,F10.2,A,3F12.4)') 'takestep> angle=',ANGLE,' axis: ',RVEC(1:3)
C
C  Rotate all the non-core atoms through ANGLE about RVEC. Use rotation formula
C  from Goldstein p. 165.
C  
         DO J1=NCORE+1,NATOMS
            XL=COORDS(3*(J1-1)+1,NP); YL=COORDS(3*(J1-1)+2,NP); ZL=COORDS(3*(J1-1)+3,NP)
            DUMMY=SQRT((XL-XC)**2+(YL-YC)**2+(ZL-ZC)**2)

            RDOTN=(XL-XC)*RVEC(1)+(YL-YC)*RVEC(2)+(ZL-ZC)*RVEC(3)
            TX=(XL-XC)*COST + RVEC(1)*RDOTN*(1.0D0-COST)-((YL-YC)*RVEC(3)-(ZL-ZC)*RVEC(2))*SINT
            TY=(YL-YC)*COST + RVEC(2)*RDOTN*(1.0D0-COST)-((ZL-ZC)*RVEC(1)-(XL-XC)*RVEC(3))*SINT
            TZ=(ZL-ZC)*COST + RVEC(3)*RDOTN*(1.0D0-COST)-((XL-XC)*RVEC(2)-(YL-YC)*RVEC(1))*SINT
            IF (DUMMY.GT.0.1D0) THEN
               COORDS(3*(J1-1)+1,NP)=(XC+TX)*(DUMMY+0.5D0)/DUMMY
               COORDS(3*(J1-1)+2,NP)=(YC+TY)*(DUMMY+0.5D0)/DUMMY
               COORDS(3*(J1-1)+3,NP)=(ZC+TZ)*(DUMMY+0.5D0)/DUMMY
            ELSE 
               COORDS(3*(J1-1)+1,NP)=(XC+TX)
               COORDS(3*(J1-1)+2,NP)=(YC+TY)
               COORDS(3*(J1-1)+3,NP)=(ZC+TZ)
            ENDIF
         ENDDO
         RETURN
      ENDIF
      IF (NSURFMOVES.GE.SHELLMOVEMAX) THEN
         SHELLMOVES=.FALSE.
         NCORE=0
      ENDIF

      DO J1=1,NATOMS
         IF (SHELLMOVES.AND.(J1.LE.NCORE)) CYCLE
10       J2=3*J1
         LOCALSTEP=STEP(NP)

         IF (RIGID.AND.(J1.GT.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
         ELSE IF (RIGID.AND.(J1.LE.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
         ENDIF
C
C  Angular move block.
C  If NORESET is .TRUE. then VAT won;t be set, so we should skip this block.
C
         IF ((((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX))).AND.(.NOT.RIGID).AND.(.NOT.BLNT).AND. 
     &         (.NOT.DIFFRACTT).AND.(.NOT.GAUSST) 
     &        .AND.(.NOT.NORESET).AND.(.NOT.PERIODIC).AND.(.NOT.THOMSONT)) THEN
 
            IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 'angular move for atom ',J1, 
     &           ' V=',VMAX,' Vmin=',VMIN,' next most weakly bound atom is ',JMAX2,' V=',VMAX2

           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
C
C  Evaporation is judged from the origin, not the centre of mass. We don't want the
C  angular move to cause evaporation. Obviously this will cause problems if we have a cluster that drifts
C  away from the origin up to the container radius.  
C
              COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
              COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
              COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
              IF (DUMMY.GT.RADIUS) THEN
                 DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
                 COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
                 COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
                 COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
              ENDIF
         ELSE IF (.NOT.FIXD) THEN
C
C  Tried changing to scale steps according to distance from CM. Maximum
C  allowed shift is linear with this distance. Worse.
C  Now try moving atoms according to how strongly bound they are.
C
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
C
C Stop atoms leaving the container in this step
C
           IF ((.NOT.PERIODIC).AND.(.NOT.AMBER).AND.(.NOT.(RIGID.AND.((J1.GT.NATOMS/2)))).AND.(.NOT.BLNT)
     1                        .AND.(.NOT.DIFFRACTT).AND.(.NOT.THOMSONT).AND.(.NOT.GAUSST)) THEN
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
C
C  Simply rescaling the radius of an atom that leaves the container will bias the sampling
C  of configuration space.
C
              IF (DUMMY.GT.RADIUS) THEN
                 COORDS(J2-2,NP)=COORDSO(J2-2,NP)
                 COORDS(J2-1,NP)=COORDSO(J2-1,NP)
                 COORDS(J2,NP)=COORDSO(J2,NP)
                 DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
                 IF (DEBUG) WRITE(MYUNIT,'(A,I5,3F20.10)') 'J1,DUMMY,RADIUS=',J1,DUMMY,RADIUS,DUMMY2
                 GOTO 10
              ENDIF
           ENDIF
         ENDIF
         IF (TOSI.OR.WELCH) THEN
            DO J3=1,J1-1
               DUMMY=(COORDS(J2-2,NP)-COORDS(3*(J3-1)+1,NP))**2
     1              +(COORDS(J2-1,NP)-COORDS(3*(J3-1)+2,NP))**2
     2              +(COORDS(J2  ,NP)-COORDS(3*(J3-1)+3,NP))**2
               IF (DUMMY.LT.1.0D0) GOTO 10 
            ENDDO
         ENDIF
      ENDDO
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))

      RETURN
      END
