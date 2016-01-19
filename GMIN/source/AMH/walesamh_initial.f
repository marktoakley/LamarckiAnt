      subroutine walesamh_initial()

      use amhglobals, only: SO, nmtemp,itgrd,temgrd,temtur,ictemp,ctemp,
     *  iscltab,nmres,oarchv,nmstep,numpro,idigns,maxpro,maxcrd,prcord,ires,oconv,
     *  omovi,omoviseg,quench,nquench,quench_crd

      use commons, only: coords

      implicit none

c     subroutines required by main program

       external gentab,initil,intstr,scltab,zero_amh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     internal variables:

         integer jstrt,jfins,i_quench,len,ishkit,nmdifv,gly_count
         integer iii

         character*10 save_name


      call zero_amh

c     --------------------- begin -----------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     open required files, read input parameter file, and generate header file
c      call read_input_alt() ! called BEFORE initil : Johan
      call initil
c      call read_altgamma()  ! called AFTER  initil : Johan
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up temperature-annealing schedule
c      call annsch(nmtemp,itgrd,temgrd,temtur,ictemp,ctemp)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate requisite force/potential tables
      call gentab
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(SO,*) 'in scltab',iscltab
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scale tables
      if (iscltab) call scltab
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(SO,*) 'out scltab'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate initial structures
      quench_crd=0.0D0
      call intstr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        write(SO,*) 'out intstr'

        if (.not. quench) nquench=1
        do i_quench = 1,nquench
        prcord=quench_crd(:,:,:,:,i_quench)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set indicies for the first and last residues
c     which are not fixed in crystal conformation
      jstrt=1
      jfins=nmres

c     set subsegment length
      len=jfins - jstrt + 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c++++++++++++++++++++++++++++++johan
!       call read_input_alt()
!       call read_altgamma() ! must be called after initil

c------------------------------johan
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- diagnostics ---
c      write(oarchv,121)jstrt,jfins,len
c  121 format(/' start ',i3,' end ',i3,' length ',i3)
c      write(oarchv,122)jstrt,jfins,nmstep,numpro
c  122 format('jstrt ',i3,' jfins ',i3,' mutations/T ',
c     *       i3,' numpro ',i3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     generate initial ensemble of proteins
c     find configuration which satisfies the constraints

      idigns=.false.

        enddo ! i_quench

        gly_count = 0
       do 764 iii = 1,nmres
       if (ires(iii).eq.8) then
         coords(9*(iii-1)+1 - gly_count*3,1) = (prcord(iii, 1, 1, 1)) !  CA X
         coords(9*(iii-1)+2 - gly_count*3,1) = (prcord(iii, 2, 1, 1)) !  CA Y
         coords(9*(iii-1)+3 - gly_count*3,1) = (prcord(iii, 3, 1, 1)) !  CA Z
         coords(9*(iii-1)+4 - gly_count*3,1) = (prcord(iii, 1, 1, 3)) !  O X
         coords(9*(iii-1)+5 - gly_count*3,1) = (prcord(iii, 2, 1, 3)) !  O Y
         coords(9*(iii-1)+6 - gly_count*3,1) = (prcord(iii, 3, 1, 3)) !  O Z
         gly_count = gly_count + 1
       else
         coords(9*(iii-1)+1 - gly_count*3,1) = (prcord(iii, 1, 1, 1)) !  CA X
         coords(9*(iii-1)+2 - gly_count*3,1) = (prcord(iii, 2, 1, 1)) !  CA Y
         coords(9*(iii-1)+3 - gly_count*3,1) = (prcord(iii, 3, 1, 1)) !  CA Z
         coords(9*(iii-1)+4 - gly_count*3,1) = (prcord(iii, 1, 1, 2)) !  CB X
         coords(9*(iii-1)+5 - gly_count*3,1) = (prcord(iii, 2, 1, 2)) !  CB Y
         coords(9*(iii-1)+6 - gly_count*3,1) = (prcord(iii, 3, 1, 2)) !  CB Z
         coords(9*(iii-1)+7 - gly_count*3,1) = (prcord(iii, 1, 1, 3)) !  O X
         coords(9*(iii-1)+8 - gly_count*3,1) = (prcord(iii, 2, 1, 3)) !  O Y
         coords(9*(iii-1)+9 - gly_count*3,1) = (prcord(iii, 3, 1, 3)) !  O Z
       endif
764    continue



      end
