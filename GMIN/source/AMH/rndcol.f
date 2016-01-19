
c     --------------------- rndcol ----------------------

      subroutine rndcol(nmres,maxpro,numpro,
     *                  ires,procnt,numcrd,prcord,iseed_amh)

c     ---------------------------------------------------

c     RNDCOL generate random coil configurations

c     arguments:

c        maxsiz- maximum target protein size (i)
c        maxpro- maximum number of trial structures (i)
c        numpro- actual number of trial structures (i)
c        procnt- number of random coil configurations to
c                be constructed (i)
c        numcrd- number of coordinate types (i)
c        prcord- random coil strucures (o)
c        bondln- bond lengths for coordinates types (i)

c     ---------------------------------------------------

      use amhglobals,  only:maxsiz,pdb,quench_crd,x_mcp,aminoa

      implicit none

c     passed argument declarations:

      integer nmres,numpro,procnt,numcrd,maxpro,ires(maxsiz),gly_c,iseed_amh(4)

      double precision prcord(maxsiz,3,maxpro,numcrd)

      integer temp_mcp_pass
      logical touch

ccccccccccccccccccccccccccccccc
c     internal variables
       double precision pi,pi2,dist1,dist2,
     *     nitcord(maxsiz,3),cprcord(maxsiz,3),phi,cphi,sphi,cpsi,spsi,
     *     xprime(3),yprime(3),zprime(3),xlen,ylen,zlen

         integer try(maxsiz,2), iii, jjj

         character*3 res_type

c        --- do loop indices ---

         integer i501,i503,i506,i507,i1

c My changes below   - Skip
c       double precision ran_num_array(1),ran_num
       double precision ran_num
       real ran_num_array(1)

       external SLARNV

c     --------------------- begin -----------------------
        pi=3.141592653589793238462643383279502884197
        pi2=2.0D0*pi
        do 600 i1=1,maxsiz
          try(i1,1)=0
          try(i1,2)=0
600        continue


        procnt=procnt + 1
        do 503 i503=procnt,numpro

c        generate initial configurations for residue 1

        prcord(1,1,i503,1)=0.0D0
        prcord(1,2,i503,1)=0.0D0
        prcord(1,3,i503,1)=0.0D0
        cprcord(1,1)=0.0D0
        cprcord(1,2)=0.0D0
        cprcord(1,3)=1.53D0
        nitcord(1,1)=0.0D0
        nitcord(1,2)=-1.3859293D0
        nitcord(1,3)=-0.4899999D0
        prcord(1,1,i503,2)=-1.2574047D0
        prcord(1,2,i503,2)=0.7259629D0
        prcord(1,3,i503,2)=-0.5133333D0
        
        if (ires(1).eq.8) then
          prcord(1,1,i503,2)=prcord(1,1,i503,1)
          prcord(1,2,i503,2)=prcord(1,2,i503,1)
          prcord(1,3,i503,2)=prcord(1,3,i503,1)
        endif

        i501=1

         do while (i501.lt.nmres)

c           generate random numbers for x- and
c           y-coordinates

105      zprime(1)=(cprcord(i501,1)-prcord(i501,1,i503,1))
         zprime(2)=(cprcord(i501,2)-prcord(i501,2,i503,1))
         zprime(3)=(cprcord(i501,3)-prcord(i501,3,i503,1))

         zlen=dsqrt(zprime(1)*zprime(1)+zprime(2)*zprime(2)+zprime(3)*zprime(3))

         zprime(1)=zprime(1)/zlen
         zprime(2)=zprime(2)/zlen
         zprime(3)=zprime(3)/zlen

         yprime(1)=(prcord(i501,1,i503,1)-nitcord(i501,1)-zprime(1)*0.4899999D0)
         yprime(2)=(prcord(i501,2,i503,1)-nitcord(i501,2)-zprime(2)*0.4899999D0)
         yprime(3)=(prcord(i501,3,i503,1)-nitcord(i501,3)-zprime(3)*0.4899999D0)

         ylen=dsqrt(yprime(1)*yprime(1)+yprime(2)*yprime(2)+yprime(3)*yprime(3))

         yprime(1)=yprime(1)/ylen
         yprime(2)=yprime(2)/ylen
         yprime(3)=yprime(3)/ylen

         xprime(1)=yprime(2)*zprime(3)-yprime(3)*zprime(2)
         xprime(2)=yprime(3)*zprime(1)-yprime(1)*zprime(3)
         xprime(3)=yprime(1)*zprime(2)-yprime(2)*zprime(1)

         xlen=dsqrt(xprime(1)*xprime(1)+xprime(2)*xprime(2)+xprime(3)*xprime(3))

         xprime(1)=xprime(1)/xlen
         xprime(2)=xprime(2)/xlen
         xprime(3)=xprime(3)/xlen

c 110               phi=rand()*1.92+1.57
              
         temp_mcp_pass=1

         call SLARNV(temp_mcp_pass, iseed_amh,temp_mcp_pass,ran_num_array)
               ran_num=ran_num_array(1)

110            phi=ran_num*1.92D0+1.57D0

               if (phi.gt.pi) phi=phi+2.09D0
               cphi=dcos(phi)
               sphi=dsin(phi)

               prcord(i501,1,i503,3)
     *                  =-1.0515796D0*sphi*xprime(1)
     *                   +1.0515796D0*cphi*yprime(1)
     *                   +0.6570998D0*zprime(1)
     *                   +cprcord(i501,1)
               prcord(i501,2,i503,3)
     *                  =-1.0515796D0*sphi*xprime(2)
     *                   +1.0515796D0*cphi*yprime(2)
     *                   +0.6570998D0*zprime(2)
     *                   +cprcord(i501,2)
               prcord(i501,3,i503,3)
     *                  =-1.0515796D0*sphi*xprime(3)
     *                   +1.0515796D0*cphi*yprime(3)
     *                   +0.6570998D0*zprime(3)
     *                   +cprcord(i501,3)

               nitcord(i501+1,1)
     *                   =1.20588D0*sphi*xprime(1)
     *                   -1.20588D0*cphi*yprime(1)
     *                   +0.5368923D0*zprime(1)
     *                   +cprcord(i501,1)
               nitcord(i501+1,2)
     *                   =1.20588D0*sphi*xprime(2)
     *                   -1.20588D0*cphi*yprime(2)
     *                   +0.5368923D0*zprime(2)
     *                   +cprcord(i501,2)
               nitcord(i501+1,3)
     *                   =1.20588D0*sphi*xprime(3)
     *                   -1.20588D0*cphi*yprime(3)
     *                   +0.5368923D0*zprime(3)
     *                   +cprcord(i501,3)

               prcord(i501+1,1,i503,1)
     *                   =1.4358387D0*sphi*xprime(1)
     *                   -1.4358387D0*cphi*yprime(1)
     *                   +1.9887942D0*zprime(1)
     *                   +cprcord(i501,1)
               prcord(i501+1,2,i503,1)
     *                   =1.4358387D0*sphi*xprime(2)
     *                   -1.4358387D0*cphi*yprime(2)
     *                   +1.9887942D0*zprime(2)
     *                   +cprcord(i501,2)
               prcord(i501+1,3,i503,1)
     *                   =1.4358387D0*sphi*xprime(3)
     *                   -1.4358387D0*cphi*yprime(3)
     *                   +1.9887942D0*zprime(3)
     *                   +cprcord(i501,3)

c      make chain self-avoiding

            touch=.false.

            do 506 i506=1,i501-1

             dist1=dsqrt( (prcord(i501+1,1,i503,1)-prcord(i506,1,i503,1))**2
     *                   +(prcord(i501+1,2,i503,1)-prcord(i506,2,i503,1))**2
     *                   +(prcord(i501+1,3,i503,1)-prcord(i506,3,i503,1))**2)

             dist2=dsqrt( (prcord(i501+1,1,i503,1)-prcord(i506,1,i503,2))**2
     *                   +(prcord(i501+1,2,i503,1)-prcord(i506,2,i503,2))**2
     *                   +(prcord(i501+1,3,i503,1)-prcord(i506,3,i503,2))**2)

                if ( (dist1.lt.5.0D0).or.(dist2.lt.4.0D0) ) 
     *             touch=.true.

506              continue

                if (touch) then

                   try(i501,1)=try(i501,1)+1
                   if (try(i501,1).lt.11) goto 110
                   
c           failed 10 times-backtrack to previous good bond

c                   write (6,*) i501,1,try(i501,1)

                   do 301 i501=i501-1,1,-1

                   try(i501+1,1)=0
                   try(i501,2)=try(i501,2)+1
                   if (try(i501,2).lt.11) goto 115 

                   try(i501,2)=0
                   try(i501,1)=try(i501,1)+1
                   if (try(i501,1).lt.11) goto 105
c                   write (6,*) i501,1,try(i501,1)

301                continue        

                   stop ' failed in rndcol '

                endif

115    zprime(1)=(prcord(i501+1,1,i503,1)-nitcord(i501+1,1))/1.47D0
       zprime(2)=(prcord(i501+1,2,i503,1)-nitcord(i501+1,2))/1.47D0
       zprime(3)=(prcord(i501+1,3,i503,1)-nitcord(i501+1,3))/1.47D0

       zlen=dsqrt(zprime(1)*zprime(1)+zprime(2)*zprime(2)+zprime(3)*zprime(3))

               zprime(1)=zprime(1)/zlen
               zprime(2)=zprime(2)/zlen
               zprime(3)=zprime(3)/zlen

               yprime(1)=(nitcord(i501+1,1)-cprcord(i501,1)
     *                     -zprime(1)*0.7189235D0)
     *                     /1.1070451D0
               yprime(2)=(nitcord(i501+1,2)-cprcord(i501,2)
     *                     -zprime(2)*0.7189235D0)
     *                     /1.1070451D0
               yprime(3)=(nitcord(i501+1,3)-cprcord(i501,3)
     *                     -zprime(3)*0.7189235D0)
     *                     /1.1070451D0

        ylen=sqrt(yprime(1)*yprime(1)+yprime(2)*yprime(2)+yprime(3)*yprime(3))

        yprime(1)=yprime(1)/ylen
        yprime(2)=yprime(2)/ylen
        yprime(3)=yprime(3)/ylen

        xprime(1)=yprime(2)*zprime(3)-yprime(3)*zprime(2)
        xprime(2)=yprime(3)*zprime(1)-yprime(1)*zprime(3)
        xprime(3)=yprime(1)*zprime(2)-yprime(2)*zprime(1)

        xlen=sqrt(xprime(1)*xprime(1)+xprime(2)*xprime(2)+xprime(3)*xprime(3))

        xprime(1)=xprime(1)/xlen
        xprime(2)=xprime(2)/xlen
        xprime(3)=xprime(3)/xlen

c  120               phi=rand()*1.75+3.49

               call SLARNV(1, iseed_amh,1,ran_num_array(1))
               ran_num=ran_num_array(1)

120            phi=ran_num*1.75D0+3.49D0
               cphi=dcos(phi)
               sphi=dsin(phi)
               cpsi=dcos(phi-pi2/3.0D0)
               spsi=dsin(phi-pi2/3.0D0)

               prcord(i501+1,1,i503,2)
     *                  =1.4519259D0*spsi*xprime(1)-1.4519259D0*cpsi*yprime(1)
     *                  +0.5133333D0*zprime(1)+prcord(i501+1,1,i503,1)
               prcord(i501+1,2,i503,2)
     *                  =1.4519259D0*spsi*xprime(2)-1.4519259D0*cpsi*yprime(2)
     *                  +0.5133333D0*zprime(2)+prcord(i501+1,2,i503,1)
               prcord(i501+1,3,i503,2)
     *                  =1.4519259D0*spsi*xprime(3)-1.4519259D0*cpsi*yprime(3)
     *                  +0.5133333D0*zprime(3)+prcord(i501+1,3,i503,1)

        if (ires(i501+1).eq.8) then
          prcord(i501+1,1,i503,2)=prcord(i501+1,1,i503,1)
          prcord(i501+1,2,i503,2)=prcord(i501+1,2,i503,1)
          prcord(i501+1,3,i503,2)=prcord(i501+1,3,i503,1)
        endif

c      make chain self-avoiding

            touch=.false.

            do 507 i507=1,i501-1

                dist1=dsqrt( (prcord(i501+1,1,i503,2)
     *                      -prcord(i507,1,i503,1))**2
     *                     +(prcord(i501+1,2,i503,2)
     *                      -prcord(i507,2,i503,1))**2
     *                     +(prcord(i501+1,3,i503,2)
     *                      -prcord(i507,3,i503,1))**2)

                dist2=dsqrt( (prcord(i501+1,1,i503,2)
     *                      -prcord(i507,1,i503,2))**2
     *                     +(prcord(i501+1,2,i503,2)
     *                      -prcord(i507,2,i503,2))**2
     *                     +(prcord(i501+1,3,i503,2)
     *                      -prcord(i507,3,i503,2))**2)

                if ( (dist1.lt.4.0D0).or.(dist2.lt.3.0D0) ) 
     *             touch=.true.

507        continue

                if (touch) then

                   try(i501,2)=try(i501,2)+1
                   if (try(i501,2).lt.11) goto 120

c           failed 10 times-backtrack to previous good bond

c                   write (6,*) i501,2,try(i501,2)
                   try(i501,2)=0
                   try(i501,1)=try(i501,1)+1
                   if (try(i501,1).lt.11) goto 105
c                   write (6,*) i501,1,try(i501,1)

                   do 302 i501=i501-1,1,-1

                   try(i501+1,1)=0
                   try(i501,2)=try(i501,2)+1
                   if (try(i501,2).lt.11) goto 115
c                   write (6,*) i501,2,try(i501,2)

                   try(i501,2)=0
                   try(i501,1)=try(i501,1)+1
                   if (try(i501,1).lt.11) goto 105
                   write (6,*) i501,1,try(i501,1)

302             continue

                   stop ' failed in rndcol '

                endif

       cprcord(i501+1,1)=1.4424978D0*sphi*xprime(1)-1.4424978D0*cphi*yprime(1)
     *                  +0.51D0*zprime(1)+prcord(i501+1,1,i503,1)
       cprcord(i501+1,2)=1.4424978D0*sphi*xprime(2)-1.4424978D0*cphi*yprime(2)
     *                  +0.51D0*zprime(2)+prcord(i501+1,2,i503,1)
       cprcord(i501+1,3)=1.4424978D0*sphi*xprime(3)-1.4424978D0*cphi*yprime(3)
     *                  +0.51D0*zprime(3)+prcord(i501+1,3,i503,1)

        i501=i501+1
       end do

 503   continue

         open(unit=pdb,file='0.pdb',form='formatted',status = 'unknown')

       do 664 iii = 1,nmres
            res_type = 'ALA'
             if (ires(iii) .eq. 8) then
              res_type = 'GLY'
             end if

CACACACACACACACACACACAC
       write(pdb,665) iii, res_type,iii,(prcord(iii, jjj, 1, 1), jjj =1,3), iii
665    format('ATOM    ',i3,'  CA  ',a3,'   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)
C'C'C'C'C'C'C'C'C'C'C'C'
       write(pdb,671) iii, res_type,iii, (cprcord(iii, jjj),jjj =1,3), iii
671   format('ATOM    ',i3,'  C   ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
COOOOOOOOOOOOOOOOOOOOOOOO
       if (iii .lt. nmres) then
        write(pdb,670) iii, res_type,iii, (prcord(iii, jjj, 1,3),jjj =1,3), iii
670     format('ATOM    ',i3,'  O   ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
       end if
CNNNNNNNNNNNNNNNNNNNNNNN
        write(pdb,669) iii, res_type, iii, (nitcord(iii, jjj),jjj =1,3), iii
669     format('ATOM    ',i3,'  N   ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)

CBCBCBCBCBCBCBCBCBCBCBCB
       if (ires(iii) .ne. 8) then
        write(pdb,668) iii, iii, (prcord(iii, jjj, 1, 2),jjj =1,3), iii
668     format('ATOM    ',i3,'  CB  ALA   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
       end if
ccccccccccccccccccccccccc

664     continue
        close(pdb)

c  MCP structure for GMIN

       gly_c = 0

       do 764 iii = 1,nmres

       if (ires(iii).eq.8) gly_c = gly_c +1

c  prcord(iii, jjj, 1, 1)  = prcord(res id , coord id, num pro id , atom id)
c  convert nmres to num of atoms 
c  assume x(CA(x,y,z),CB(x,y,z),O(x,y,z))  

       res_type = aminoa(ires(iii))

       if (ires(iii).eq.8) then
        x_mcp(9*(iii-1)+1-(gly_c-1)*3) = (prcord(iii, 1, 1, 1))   !  CA X
        x_mcp(9*(iii-1)+2-(gly_c-1)*3) = (prcord(iii, 2, 1, 1))   !  CA Y
        x_mcp(9*(iii-1)+3-(gly_c-1)*3) = (prcord(iii, 3, 1, 1))   !  CA Z
!       x_mcp(9*(iii-1)+4) = (prcord(iii, 1, 1, 2))   !  CB X
!       x_mcp(9*(iii-1)+5) = (prcord(iii, 2, 1, 2))   !  CB Y
!       x_mcp(9*(iii-1)+6) = (prcord(iii, 3, 1, 2))   !  CB Z
        x_mcp(9*(iii-1)+4-(gly_c-1)*3) = (prcord(iii, 1, 1, 3))   !   O X
        x_mcp(9*(iii-1)+5-(gly_c-1)*3) = (prcord(iii, 2, 1, 3))   !   O Y
        x_mcp(9*(iii-1)+6-(gly_c-1)*3) = (prcord(iii, 3, 1, 3))   !   O Z
      else
        x_mcp(9*(iii-1)+1-gly_c*3) = (prcord(iii, 1, 1, 1))   !  CA X
        x_mcp(9*(iii-1)+2-gly_c*3) = (prcord(iii, 2, 1, 1))   !  CA Y
        x_mcp(9*(iii-1)+3-gly_c*3) = (prcord(iii, 3, 1, 1))   !  CA Z
        x_mcp(9*(iii-1)+4-gly_c*3) = (prcord(iii, 1, 1, 2))   !  CB X
        x_mcp(9*(iii-1)+5-gly_c*3) = (prcord(iii, 2, 1, 2))   !  CB Y
        x_mcp(9*(iii-1)+6-gly_c*3) = (prcord(iii, 3, 1, 2))   !  CB Z
        x_mcp(9*(iii-1)+7-gly_c*3) = (prcord(iii, 1, 1, 3))   !   O X
        x_mcp(9*(iii-1)+8-gly_c*3) = (prcord(iii, 2, 1, 3))   !   O Y
        x_mcp(9*(iii-1)+9-gly_c*3) = (prcord(iii, 3, 1, 3))   !   O Z
      endif

CACACACACACACACACACACAC
      if (ires(iii) .ne. 8) then
c        write(6,765) iii, res_type,iii,x_mcp(9*(iii-1)+1-gly_c*3),x_mcp(9*(iii-1)+2-gly_c*3),x_mcp(9*(iii-1)+3-gly_c*3), iii
765     format('ATOM    ',i3,'  CA  ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
      endif

      if (ires(iii) .eq. 8) then
c      write(6,865)iii,res_type,iii,x_mcp(9*(iii-1)+1-(gly_c-1)*3),x_mcp(9*(iii-1)+2-(gly_c-1)*3),x_mcp(9*(iii-1)+3-(gly_c-1)*3),iii
865     format('ATOM    ',i3,'  XA  ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00     TPDB ',i3)
      endif

CBCBCBCBCBCBCBCBCBCBCBCB
       if (ires(iii) .ne. 8) then
c        write(6,768) iii, res_type,iii, x_mcp(9*(iii-1)+4- gly_c*3),x_mcp(9*(iii-1)+5-gly_c*3),x_mcp(9*(iii-1)+6-gly_c*3) , iii
768     format('ATOM    ',i3,'  CB  ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
       end if

COOOOOOOOOOOOOOOOOOOOOOOO
c       if (iii .lt. nmres) then
      if (ires(iii) .ne. 8) then
c        write(6,770) iii, res_type,iii, x_mcp(9*(iii-1)+7-gly_c*3),x_mcp(9*(iii-1)+8-gly_c*3),x_mcp(9*(iii-1)+9-gly_c*3), iii
770     format('ATOM    ',i3,'  O   ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
       endif

      if (ires(iii) .eq. 8) then
c       write(6,965)iii,res_type,iii,x_mcp(9*(iii-1)+4-(gly_c-1)*3),x_mcp(9*(iii-1)+5-(gly_c-1)*3),x_mcp(9*(iii-1)+6-(gly_c-1)*3),iii
965     format('ATOM    ',i3,'  X0  ',a3,'   ',i3,'    ',3(f8.3),'  1.00  0.00      TPDB ',i3)
      endif

ccccccccccccccccccccccccc

764     continue

c      write(6,*)'out rndcol'

c       do 364 iii = 1,nmres*3*3
c           write(6,*)'coord   res',  x_mcp(iii), iii
c364    continue 

   
      quench_crd(:,:,:,:,1)=prcord

      return
      end
