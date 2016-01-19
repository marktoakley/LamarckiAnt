
C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C***********************************************************************
C
C Subroutine for Ackland metal potentials. Periodic boundary conditions!
C
C***********************************************************************
C
      SUBROUTINE ACK(X,VNEW,PSC,GTEST)
      USE COMMONS, ONLY : ACKLANDID, NATOMS, BOXLX, BOXLY, BOXLZ, CUTOFF, VT
      IMPLICIT NONE
      INTEGER N, J1, J2, J, I, MZ, MY, MX 
      DOUBLE PRECISION X(3*NATOMS), POTA, POTB, DIST, PSC, VNEW(3*NATOMS),Rc
      COMMON /param_cut_off/Rc	

      DOUBLE PRECISION RHO(3*NATOMS)
      DOUBLE PRECISION VEC(NATOMS,NATOMS,3),Rneigh(NATOMS,NATOMS,3),force(3,NATOMS)
      double precision rbuf,r,zero,norm,rho_temp,vpot_temp,Fembed_d_i
      double precision Vpot,Vpot_d,rho_pot,rho_pot_d,Fembed,Fembed_d
      integer ipot,icount,ja
      integer ic(NATOMS),neigh_type(NATOMS,NATOMS),ndir(NATOMS)
      LOGICAL GTEST

      N=NATOMS 
      Rc=CUTOFF
!
! IF ACKLANDID is negative we use cluster boundary conditions!
!
      ipot=ABS(ACKLANDID)
      zero=1.0D-12
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      icount=0
      ic(:)=0
      
      IF (ACKLANDID.GT.0) THEN
         DO 25 J1=1,N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
	    
	    icount=ic(J1)
            
	    DO 15 J2=J1+1,N
               VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
               VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
               VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
               MX=NINT(VEC(J2,J1,1)/BOXLX)
               MY=NINT(VEC(J2,J1,2)/BOXLY)
               MZ=NINT(VEC(J2,J1,3)/BOXLZ)
               VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(MX)
               VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(MY)
               VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(MZ)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
	       norm=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 + VEC(J1,J2,3)**2
	       if (norm < Rc*Rc.and.norm>zero) then
	       icount=icount+1
	       ic(J2)=ic(J2)+1
	       Rneigh(J2,ic(J2),1)=VEC(J1,J2,1)
	       Rneigh(J2,ic(J2),2)=VEC(J1,J2,2)
	       Rneigh(J2,ic(J2),3)=VEC(J1,J2,3)
	       
	       Rneigh(J1,icount,1)=VEC(J2,J1,1)
	       Rneigh(J1,icount,2)=VEC(J2,J1,2)
	       Rneigh(J1,icount,3)=VEC(J2,J1,3)
	       
	       neigh_type(J1,icount)=J2
	       neigh_type(J2,ic(J2))=J1
	       
	       end if
15          CONTINUE
         ndir(J1)=icount
25       CONTINUE
      ELSE ! cluster case !   
         DO J1=1,N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
	    
	    icount=ic(J1)
            
	    DO J2=J1+1,N
               VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
               VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
               VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
	       norm=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 + VEC(J1,J2,3)**2
	       IF (NORM < Rc*Rc.AND.NORM>ZERO) THEN
	          ICOUNT=ICOUNT+1
	          IC(J2)=IC(J2)+1
	          Rneigh(J2,ic(J2),1)=VEC(J1,J2,1)
	          Rneigh(J2,ic(J2),2)=VEC(J1,J2,2)
	          Rneigh(J2,ic(J2),3)=VEC(J1,J2,3)
	       
	          Rneigh(J1,icount,1)=VEC(J2,J1,1)
	          Rneigh(J1,icount,2)=VEC(J2,J1,2)
	          Rneigh(J1,icount,3)=VEC(J2,J1,3)
	       
	          neigh_type(J1,icount)=J2
	          neigh_type(J2,ic(J2))=J1
	       
	       end if
            ENDDO
            ndir(J1)=icount
         ENDDO

      ENDIF
C
C Store density matrix: In the case of the perfect fcc lattice,
C the infinitely extended crystal implies that every RHO(J) is
C equal to RHO(1).
C
      DO 11 I=1,N
	 
	 rho_temp=0.d0
         vpot_temp=0.d0
	 
	 DO 122 J=1,ndir(I)
	     
	     ja=neigh_type(I,J)
	     
	     rbuf=dsqrt(Rneigh(I,J,1)**2+Rneigh(I,J,2)**2+Rneigh(I,J,3)**2)
	     
	     rho_temp  = rho_temp  +  rho_pot(ipot,rbuf)
	     vpot_temp = vpot_temp +  Vpot(ipot,rbuf)
	     
122      CONTINUE
       RHO(I)=rho_temp
       VT(I)=vpot_temp
!C        write(*,*) I, RHO(I)
11    CONTINUE
C
C calculate the potential energy:
C
      POTA=0.0D0
      POTB=0.0D0
      DO 13 I=1,N
        POTA=POTA+VT(I)
        POTB=POTB - Fembed(ipot,RHO(I))
13    CONTINUE
      PSC=POTA - POTB

Cdebug      write(*,*) POTA, -POTB, PSC
C      stop
      
      force(:,:)=0.d0
      
      if (GTEST) then
         do I=1,N
	    Fembed_d_i=Fembed_d(ipot,RHO(I))
	    do J=1,ndir(I)
	       ja=neigh_type(I,J)
	       rbuf=dsqrt(Rneigh(I,J,1)**2+Rneigh(I,J,2)**2+Rneigh(I,J,3)**2)
   
               force(1,I)= force(1,I)+Rneigh(I,J,1)*(2*Vpot_d(ipot,rbuf)+      !&
     1              (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                 !&
     1		     rho_pot_d(ipot,rbuf)  ) /rbuf
               force(2,I)= force(2,I)+Rneigh(I,J,2)*(2*Vpot_d(ipot,rbuf)+      !&
     1               (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                !&
     1		     rho_pot_d(ipot,rbuf)  ) /rbuf
               force(3,I)= force(3,I)+Rneigh(I,J,3)*(2*Vpot_d(ipot,rbuf)+      !&
     1               (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                !&
     1		     rho_pot_d(ipot,rbuf)  ) /rbuf
            end do  
	    Vnew(3*(I-1)+1)=-force(1,I)   
	    Vnew(3*(I-1)+2)=-force(2,I)   
	    Vnew(3*(I-1)+3)=-force(3,I) 
	 end do    
      ENDIF
	  

      RETURN
      END
