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
      SUBROUTINE TABOO(EREAL,POTEL,P,NP,RES)
      USE commons
      IMPLICIT NONE
      INTEGER NP, J1, J2, J3
      DOUBLE PRECISION P(3*NATOMS), EREAL, POTEL, XIP
      LOGICAL RES

      CALL NEWINERTIA(P,NATOMS,NATOMS,XIP)
C
C  First check if the latest energy has been found before in any of the
C  other parallel runs. If so, combine the lists for the run which found
C  it first and reseed the other one.
C
      DO J1=1,NPAR
         IF (J1.NE.NP) THEN
            DO J2=1,NT(J1)
               IF (DABS(EREAL-ESAVE(J2,J1)).LT.ECONV) THEN
                  IF (2.0D0*DABS(XIP-XINSAVE(J2,J1))/(XIP+XINSAVE(J2,J1)).LT.1.0D-2) THEN
C                    IF (J2.EQ.NP) RETURN 
                     DO J3=1,MIN(NT(NP),NTAB-NT(J1))
                        ESAVE(NT(J1)+J3,J1)=ESAVE(J3,NP)
                        XINSAVE(NT(J1)+J3,J1)=XINSAVE(J3,NP)
                     ENDDO
                     NT(J1)=MIN(NT(J1)+NT(NP),NTAB)
                     CALL GSORT(NT(J1),ESAVE,XINSAVE,J1,NPAR,NTAB)  
                     NT(J1)=MIN(NT(J1),NTAB)
                     NT(NP)=0
C                    CALL RESEED(NATOMS,P,RADIUS)
                     CALL HSMOVE(P,1,15)
                     DO J3=1,3*NATOMS
                        COORDS(J3,NP)=P(J3)
                        COORDSO(J3,NP)=P(J3)
                     ENDDO
                     EPREV(NP)=1.0D100
                     EREAL=1.0D100
                     POTEL=1.0D100
                     WRITE(*,'(A,I2,A,I2)') 'Parallel run ',NP,' reseeded after overlap with run ',J1
                     RES=.TRUE.
                     RETURN
                  ELSE
                     PRINT*,'Energies nearly degenerate:',EREAL,ESAVE(J2,J1)
                     PRINT*,'But  different  structures:',XIP,XINSAVE(J2,J1)
                  ENDIF
               ENDIF
               IF (EREAL.LT.ESAVE(J2,J1)) GOTO 10
            ENDDO
10          CONTINUE
         ENDIF
      ENDDO

      IF (NT(NP).EQ.0) THEN
         NT(NP)=1
         ESAVE(1,NP)=EREAL
         XINSAVE(1,NP)=XIP
         RETURN
      ENDIF

      SELFT=.FALSE.
      DO J2=1,NT(NP)
         IF (DABS(EREAL-ESAVE(J2,NP)).LT.ECONV) THEN
C
C Self-avoiding walk if NPAR=1.
C
            IF (NPAR.EQ.1) THEN
               IF (DABS(XIP-XINSAVE(J2,NP))/(XIP+XINSAVE(J2,NP)).LT.1.0D-2) THEN
                  SELFT=.TRUE.
                  RETURN
               ELSE
                  PRINT*,'Energies nearly degenerate:',EREAL,ESAVE(J2,NP)
                  PRINT*,'But  different  structures:',XIP,XINSAVE(J2,NP)
               ENDIF
            ELSE
               RETURN
            ENDIF
         ENDIF
         IF (EREAL.LT.ESAVE(J2,NP)) THEN
            NT(NP)=MIN(NT(NP)+1,NTAB)
            DO J3=NT(NP),J2+1,-1
               ESAVE(J3,NP)=ESAVE(J3-1,NP)
               XINSAVE(J3,NP)=XINSAVE(J3-1,NP)
            ENDDO
            ESAVE(J2,NP)=EREAL
            XINSAVE(J2,NP)=XIP
            RETURN
         ENDIF
      ENDDO

      NT(NP)=NT(NP)+1
      ESAVE(NT(NP),NP)=EREAL
      XINSAVE(NT(NP),NP)=XIP

      RETURN
      END
