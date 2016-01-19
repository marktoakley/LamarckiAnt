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
      SUBROUTINE Gsort(n,arr,arr2,np,NPAR,NTAB)
      INTEGER n,M,NSTACK,np,NPAR,NTAB,a2
      DOUBLE PRECISION arr(NTAB,NPAR), arr2(NTAB,NPAR)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a,temp,temp2
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j,np)
C******************************
          a2=arr2(j,np)
C******************************
          do 11 i=j-1,1,-1
            if(arr(i,np).le.a)goto 2
            arr(i+1,np)=arr(i,np)
C******************************
            arr2(i+1,np)=arr2(i,np)
C******************************
11        continue
          i=0
2         arr(i+1,np)=a
C******************************
          arr2(i+1,np)=a2
C******************************
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k,np)
        arr(k,np)=arr(l+1,np)
        arr(l+1,np)=temp
C******************************
        temp2=arr2(k,np)
        arr2(k,np)=arr2(l+1,np)
        arr2(l+1,np)=temp2
C******************************
        if(arr(l+1,np).gt.arr(ir,np))then
          temp=arr(l+1,np)
          arr(l+1,np)=arr(ir,np)
          arr(ir,np)=temp
C******************************
          temp2=arr2(l+1,np)
          arr2(l+1,np)=arr2(ir,np)
          arr2(ir,np)=temp2
C******************************
        endif
        if(arr(l,np).gt.arr(ir,np))then
          temp=arr(l,np)
          arr(l,np)=arr(ir,np)
          arr(ir,np)=temp
C******************************
          temp2=arr2(l,np)
          arr2(l,np)=arr2(ir,np)
          arr2(ir,np)=temp2
C******************************
        endif
        if(arr(l+1,np).gt.arr(l,np))then
          temp=arr(l+1,np)
          arr(l+1,np)=arr(l,np)
          arr(l,np)=temp
C******************************
          temp2=arr2(l+1,np)
          arr2(l+1,np)=arr2(l,np)
          arr2(l,np)=temp2
C******************************
        endif
        i=l+1
        j=ir
        a=arr(l,np)
C******************************
        a2=arr2(l,np)
C******************************
3       continue
          i=i+1
        if(arr(i,np).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j,np).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i,np)
        arr(i,np)=arr(j,np)
        arr(j,np)=temp
C******************************
        temp2=arr2(i,np)
        arr2(i,np)=arr2(j,np)
        arr2(j,np)=temp2
C******************************
        goto 3
5       arr(l,np)=arr(j,np)
        arr(j,np)=a
C******************************
        arr2(l,np)=arr2(j,np)
        arr2(j,np)=a2
C******************************
        jstack=jstack+2
        if(jstack.gt.NSTACK) print*,'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
