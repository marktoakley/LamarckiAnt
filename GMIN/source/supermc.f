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
      SUBROUTINE SUPERMC(SPOTEL,SCOORDS,NSUPERCOUNT,POTEL)
      USE commons
      IMPLICIT NONE
      

      INTEGER J1, NSSUCCESS, NSFAIL, NSFAILT, NSSUCCESST, J2, INDEX(NSUPER), NSUPERCOUNT
      DOUBLE PRECISION RANDOM, SPOTEL(NSUPER), SCOORDS(3*NATOMS,NSUPER), SEPREV,
     1                 SCOORDSO(3*NATOMS), POTEL
      LOGICAL ASTEST
      SAVE
C
C  Save the NSUPER previous energies and coordinates.
C
      SPOTEL(NSUPERCOUNT)=POTEL
      DO J2=1,3*NATOMS
         SCOORDS(J2,NSUPERCOUNT)=COORDS(J2,1)
      ENDDO
      NSUPERCOUNT=NSUPERCOUNT-1
      IF (NSUPERCOUNT.NE.0) RETURN

      CALL SORT5(NSUPER,NSUPER,SPOTEL,INDEX)
      IF (NSUPERSTEP.EQ.0) THEN
         NSSUCCESS=0
         NSFAIL=0
         NSSUCCESST=0
         NSFAILT=0
C
C  Calculate the initial energy and save in SEPREV.
C
         SEPREV=SPOTEL(1)
         DO J1=1,3*NATOMS
            SCOORDSO(J1)=SCOORDS(J1,INDEX(1))
         ENDDO
      ELSE
         CALL STRANSITION(SPOTEL(1),SEPREV,ASTEST,RANDOM)
         IF (ASTEST) THEN
C           IF (DEBUG) THEN
               WRITE(*,34) RANDOM,SPOTEL(1),SEPREV,NSSUCCESS,NSFAIL
34             FORMAT('Super RAN,SPOTEL,SEPREV,NSSUC,NSFAIL=',3F15.7,2I6,' ACC')
C           ENDIF
            NSSUCCESS=NSSUCCESS+1
            SEPREV=SPOTEL(1)
            DO J2=1,3*(NATOMS-NSEED)
               SCOORDSO(J2)=SCOORDS(J2,INDEX(1))
            ENDDO
C           EPREV(1)=1.0D100
            EPREV(1)=SPOTEL(NSUPER)
C           CALL STAKESTEP(SCOORDS,INDEX)
            DO J2=1,3*(NATOMS-NSEED)
C              COORDS(J2,1)=SCOORDS(J2,INDEX(1))
C              COORDSO(J2,1)=SCOORDS(J2,INDEX(1))
               COORDS(J2,1)=SCOORDS(J2,INDEX(NSUPER))
               COORDSO(J2,1)=SCOORDS(J2,INDEX(NSUPER))
            ENDDO
            PRINT*,'Resetting coordinates to ',EPREV(1)
         ELSE
            NSFAIL=NSFAIL+1
            EPREV(1)=SEPREV
            DO J2=1,3*(NATOMS-NSEED)
               COORDS(J2,1)=SCOORDSO(J2)
               COORDSO(J2,1)=SCOORDS(J2,INDEX(1))
            ENDDO
C           DO J2=1,NATOMS    
C              VAT(J2,INDEX(1)?)=VATO(J2,INDEX(1)?)   ?????
C           ENDDO
C           IF (DEBUG) THEN
               WRITE(*,36) RANDOM,SPOTEL(1),SEPREV,NSSUCCESS,NSFAIL
36             FORMAT('Super RAN,SPOTEL,SEPREV,NSSUC,NSFAIL=',3F15.7,2I6,' REJ')
C           ENDIF
         ENDIF
C
C  Check the acceptance ratio.
C
         IF ((MOD(NSUPERSTEP,NSACCEPT).EQ.0).AND.(NSEED.EQ.0)) THEN
            IF (1.0D0*(NSSUCCESS)/(1.0D0*(NSSUCCESS+NSFAIL)).GT.SACCRAT) THEN
C              SUPSTEP=SUPSTEP*1.1D0
               TEMPS=TEMPS/1.1D0
            ELSE
C              SUPSTEP=SUPSTEP/1.1D0
               TEMPS=TEMPS*1.1D0
            ENDIF
C           WRITE(*,'(A,F15.7)') 'Superstep scaling factor is now ',SUPSTEP
            WRITE(*,'(A,F15.7)') 'Superstep temperature is now ',TEMPS
            WRITE(*,'(A,I4,A,F15.7)') 
     1        'Acceptance ratio for previous ',NSACCEPT,' supersteps=',1.0D0*(NSSUCCESS)/(1.0D0*(NSSUCCESS+NSFAIL))
            NSSUCCESST=NSSUCCESST+NSSUCCESS
            NSFAILT=NSFAILT+NSFAIL
            NSSUCCESS=0
            NSFAIL=0 
         ENDIF
      ENDIF
      WRITE(*,'(A,I6,A,F20.10)') 'Energy in superstep chain at step ',NSUPERSTEP,' is ',SPOTEL(1)

      NSUPERSTEP=NSUPERSTEP+1
      NSUPERCOUNT=NSUPER
     
      RETURN
      END
C
C     This subprogram performs a sort on the input data and
C     arranges it from smallest to biggest. The exchange-sort
C     algorithm is used.
C
      SUBROUTINE SORT5(N,J3,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NA(J3), NTEMP
      DOUBLE PRECISION TEMP, A(J3)
 
      DO J1=1,N
         NA(J1)=J1
      ENDDO
      DO J1=1,N-1
         L=J1
         DO J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
         ENDDO
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
      ENDDO

      RETURN
      END
C
C************************************************************************************************
C
      SUBROUTINE STRANSITION(ENEW,EOLD,ASTEST,RANDOM)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION ENEW, EOLD, DPRAND, RANDOM
      LOGICAL ASTEST

      IF (ENEW.LT.EOLD) THEN
         RANDOM=0.0D0
         ASTEST=.TRUE.
      ELSE
         RANDOM=DPRAND()
         IF (DEXP(-(ENEW-EOLD)/TEMPS).GT.RANDOM) THEN
            ASTEST=.TRUE.
         ELSE
            ASTEST=.FALSE.
         ENDIF
      ENDIF

      RETURN
      END
C
C************************************************************************************************
C
      SUBROUTINE STAKESTEP(SCOORDS,INDEX)
      USE commons
      IMPLICIT NONE
      

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, 
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN,
     2                 THETA, PHI, PI, DUMMY, SCOORDS(3*NATOMS,NSUPER)
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, INDEX(NSUPER)

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMIN=1.0D6
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)=DSQRT(SCOORDS(J2-2,INDEX(1))**2+SCOORDS(J2-1,INDEX(1))**2+SCOORDS(J2,INDEX(1))**2)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,INDEX(1)).GT.VMAX) THEN
            VMAX=VAT(J1,INDEX(1))
            JMAX=J1
         ENDIF
         IF (VAT(J1,INDEX(1)).LT.VMIN) VMIN=VAT(J1,INDEX(1))
      ENDDO

      DO J1=1,NATOMS-NSEED
         J2=3*J1
         IF ((((VAT(J1,INDEX(1)).GT.SUPSTEP*ASTEP(1)*VMIN).AND.
     1       (J1.EQ.JMAX)).OR.(NATOMS-NSEED.EQ.1))) THEN
           IF (DEBUG) PRINT*,'angular move for point ',J1
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           SCOORDS(J2-2,INDEX(1))=DMAX*DSIN(THETA)*DCOS(PHI)
           SCOORDS(J2-1,INDEX(1))=DMAX*DSIN(THETA)*DSIN(PHI)
           SCOORDS(J2,INDEX(1))=  DMAX*DCOS(THETA)
         ELSE
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           SCOORDS(J2-2,INDEX(1))=SCOORDS(J2-2,INDEX(1))+SUPSTEP*STEP(1)*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           SCOORDS(J2-1,INDEX(1))=SCOORDS(J2-1,INDEX(1))+SUPSTEP*STEP(1)*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           SCOORDS(J2,INDEX(1))=SCOORDS(J2,INDEX(1))+SUPSTEP*STEP(1)*RANDOM
C
C Stop atoms leaving the container in this step
C
           IF (.NOT.PERIODIC) THEN
              DUMMY=SCOORDS(J2-2,INDEX(1))**2+SCOORDS(J2-1,INDEX(1))**2+SCOORDS(J2,INDEX(1))**2
              IF (DUMMY.GT.RADIUS) THEN
                 SCOORDS(J2-2,INDEX(1))=SCOORDS(J2-2,INDEX(1))*DSQRT(RADIUS/DUMMY)
                 SCOORDS(J2-1,INDEX(1))=SCOORDS(J2-1,INDEX(1))*DSQRT(RADIUS/DUMMY)
                 SCOORDS(J2,INDEX(1))=SCOORDS(J2,INDEX(1))*DSQRT(RADIUS/DUMMY)
              ENDIF
           ENDIF
         ENDIF
      ENDDO
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO J1=1,NATOMS
            J2=3*J1
            XMASS=XMASS+SCOORDS(J2-2,INDEX(1))
            YMASS=YMASS+SCOORDS(J2-1,INDEX(1))
            ZMASS=ZMASS+SCOORDS(J2,INDEX(1))
         ENDDO
         XMASS=XMASS/(NATOMS)
         YMASS=YMASS/(NATOMS)
         ZMASS=ZMASS/(NATOMS)
         DO J1=1,NATOMS
            J2=3*J1
            SCOORDS(J2-2,INDEX(1))=SCOORDS(J2-2,INDEX(1))-XMASS
            SCOORDS(J2-1,INDEX(1))=SCOORDS(J2-1,INDEX(1))-YMASS
            SCOORDS(J2,INDEX(1))=  SCOORDS(J2,INDEX(1))-ZMASS
         ENDDO
      ENDIF

      RETURN
      END
