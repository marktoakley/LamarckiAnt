         program qscore
c---------------------------------------------------------------
c	This program will find the qscores related to contacts 
c       qloc refers to 8 angstrom shells of local contacts
c       qscore is a global score  
c       pe is the total potential energy
c       rg is the radius of gyration
c       ifort qscore.f calc_rama.f -o qscore_lin.out
c--------------------------------------------------------------
        
        implicit none
	
        integer maxmem,nQ,maxres
        parameter (maxres=2100)
	parameter (maxmem=5100)
        parameter (nQ=4)

	integer imem,nmem,tgnres,count,contact,
     *		tgires(maxres),icoord,ires,i_axis,i_atom,
     *		i_res,i_targ,i,j,k,l,m,mires,nmovies,init1,
     *          init2,i1,i2,iskip,iQ,atom_id,id,
     *          best_q_id,best_q_id_cut,dummyi,idummy, 
     *          gap(maxres),tempgap,open_status,
     *          num_foldon,foldstrt_min(10),foldstrt_max(10),
     *          seglist(maxres),numconst,verytemp,i_test,
     *          qpos_short_cut,qpos_medium_cut,qpos_long_cut,
     *          best_rg_id

	real del,delseg, memcord(maxres,3,5,maxmem),
     *       tgdist(maxres,maxres),
     *	     memdist(maxmem,maxmem),width0,widthexp,qpos,
     *	     qtot,qseg,qsegtot,q(maxmem),width(maxmem),
     *       q_short(maxmem),q_medium(maxmem),q_long(maxmem),
     *       qpos_short,qpos_medium,qpos_long,
     *       qtot_short,qtot_medium,qtot_long,
     *       qposlocal_short,qposlocal_medium,qposlocal_long,
     *       qtotlocal_short,qtotlocal_medium,qtotlocal_long,
     *	     ntemp(maxmem),pe(maxmem),dummy,power(1:nQ),
     *       nres,rcut,rg(1:maxmem),rgtarg,dens(1:maxmem),
     *       denstarg,qglobal_short(maxmem,nQ),
     *       qglobal_long(maxmem,nQ),qlocal_short(maxmem),
     *       qlocal_medium(maxmem), 
     *       qlocal_long(maxmem),ncontact(maxmem),
     *       norm_short(nQ),phipsiwidth, targcords(maxres,3,3),
     *       finalcords(maxres,3,3),
     *       qtot_phi, delphi,finalphi(maxres),targphi(maxres),
     *       qtot_psi, delpsi,finalpsi(maxres),targpsi(maxres),
     *       q_phi(maxmem), q_psi(maxmem), best_q,best_temp,
     *       best_q_cut,best_temp_cut,
     *       memsegcord(maxres,3,3,1),memsegdist(maxres,maxres),	
     *       dummyf,cm(maxres),ave_q,ave_q_temp, 
     *       qtot_cut,qpos_cut,
     *       qtot_short_cut,qtot_medium_cut,qtot_long_cut,
     *       q_cut(maxmem),
     *       q_short_cut(maxmem),q_medium_cut(maxmem),
     *       q_long_cut(maxmem),
     *       best_rg, best_temp_rg


       character tgname*5,pheader*25,relqpheader*33,
     *           aaname(20)*3,atom_type*2,
     *           res_type*3,ccount*8

        data aaname /'ALA','ARG','ASN','ASP','CYS',   
     *   'GLN','GLU','GLY','HIS','ILE','LEU','LYS',    
     *    'MET','PHE','PRO','SER','THR','TRP','TYR',   
     *    'VAL'/

        external calc_rama

        data phipsiwidth  / 0.52359 /  

	logical CA,CB,OX,rama,relativeq,localq,movieseg,qinterface

        data width0,widthexp / 1.0, 0.15 /  ! was 1.0, 0.5

	pheader='/home/mp466/amh/proteins/'
        relqpheader='/usr17/mprentis/data005/proteins/'

        rcut=8.0

c        write(6,*)'opening qcut.dat'
        open(90,file='qcut.dat',status='unknown')
        write(90,*) 'cutoff for short range Qs is',rcut
        close(90)

c---------------------------------------------------------------
c	Read in the target protein and memories
c--------------------------------------------------------------

c        write(6,*)'opening targlistqscore'

	open(12,file=
     *'targlistqscore',
     * status='old')
	    read(12,5)tgname
c           write(6,*)tgname
5	    format(a5)
	  read(12,*)CA
	  read(12,*)CB
          read(12,*)OX
          read(12,*)rama
          read(12,*)relativeq
          read(12,*)localq
          read(12,*)movieseg
          read(12,*)num_foldon

            if (num_foldon .gt. 10) STOP 'Too many Foldons: Max = 10'
               do i = 1,num_foldon
                 read(12,*) foldstrt_min(i),foldstrt_max(i)
c                 write(6,*) foldstrt_min(i),foldstrt_max(i)
               enddo
	close(12)

          do i = 1, maxres 
             gap(i) = 0
          enddo

          if ( num_foldon .ne. 1) then
             do i = 1,num_foldon - 1
                  gap(i) = foldstrt_min(i+1) - foldstrt_max(i) - 1
             enddo
          endif

          do i = 1,num_foldon
          if ( i .eq. 1) then
            do j = foldstrt_min(i), foldstrt_max(i)
                seglist(j - (foldstrt_min(i)-1)) = j
                numconst = j - (foldstrt_min(i)-1)
            enddo
          endif

          tempgap = 0

          if ( i .ne. 1) then
              do k = 1 , i - 1
                  tempgap = tempgap + gap(k)
              enddo
              do j = foldstrt_min(i), foldstrt_max(i)
                  seglist(j - tempgap -  (foldstrt_min(1)-1)) = j
                  numconst = j - tempgap  - (foldstrt_min(1)-1)
              enddo
          endif   ! if ( i .ne. 1) then
          enddo

c--------------------------------------------------------------
c	Open protein files
c--------------------------------------------------------------

      open(13,file='qscore_w.plot',status='unknown')
      open(33,file='qscore_o.plot',status='unknown')
      open(16,file='best_q.plot',status='unknown')
 
       if (.not.relativeq)then
c        write(6,*)'opening target structure'

c        open(12,file=tgname,status='old',iostat=open_status)

	open(12,file=pheader//tgname,status='old',iostat=open_status)
         if (open_status.ne.0) then
            write(6,*) 'problem opening', tgname
            write(6,*) 'open_status',open_status
            stop
          endif

         read (12,*)
         read (12,*) tgnres

          nres=real(tgnres)
          if (tgnres .gt. maxres) then
            write(6,*)'tgnres .gt. maxres '
             stop 
          endif
          read (12,1022) (tgires(i1),i1=1,tgnres)
1022	  format(25(i2,1x))
          do 130 icoord=1,3
            read (12,1023)(finalcords(i1,icoord,1),i1=1,tgnres)
c            write(6,*)(finalcords(i1,icoord,1),i1=1,tgnres)
1023        format(8(f8.3,1x))
130       continue
          do 140 icoord=1,3
            read (12,1023)(finalcords(i1,icoord,2),i1=1,tgnres)
140       continue
          do 150 icoord=1,3
            read (12,1023)(finalcords(i1,icoord,3),i1=1,tgnres)
150       continue
        close (12)
        endif

c   read in lpsf file format

c        write(6,*)'rel ',relativeq
c        write(6,*)'trap name ',tgname//"_trap"

        if (relativeq) then
	   open(12,file=tgname//"_trap",status='old',iostat=open_status)
          if (open_status.ne.0) then
            write(6,*) 'problem opening',pheader//tgname
            write(6,*) 'open_status',open_status
            stop
          endif

           read (12,*)
           read (12,*) tgnres
           nres=real(tgnres)
           read (12,1022) (tgires(i1),i1=1,tgnres)
           do 330 icoord=1,3
             read (12,1023)(finalcords(i1,icoord,1),i1=1,tgnres)
330        continue
           do 340 icoord=1,3
             read (12,1023)(finalcords(i1,icoord,2),i1=1,tgnres)
340        continue
           do 350 icoord=1,3
             read (12,1023)(finalcords(i1,icoord,3),i1=1,tgnres)
350        continue
           close (12)
        endif

c    read in coordinates of protein during folding run from movie file
c    memcord = atom number, coordinate, atom type, memory number

        nmem=0

c        write(6,*)'opening movie'

	open(11,file='movie',status='old',iostat=open_status)
          if (open_status.ne.0) then
            write(6,*) 'problem opening movie'
            write(6,*) 'open_status',open_status
            stop
          endif
          read(11,40) tgnres,idummy,idummy,nmovies
40        format(4(i8,1x)) 
	  do 160 imem=1,nmovies
	    nmem = nmem + 1
	    read(11,45)ntemp(imem)
c	    write(6,45)ntemp(imem)
45	    format(22x,f8.4)
             do 175 i=1,tgnres
                 read(11,20)  
     * memcord(i,1,1,imem),memcord(i,2,1,imem),memcord(i,3,1,imem),
     * memcord(i,1,2,imem),memcord(i,2,2,imem),memcord(i,3,2,imem),
     * memcord(i,1,3,imem),memcord(i,2,3,imem),memcord(i,3,3,imem)

20            format(3(3x,3(1x,f8.4),1x))
175         continue !tgnres
160	  continue !imem
	close(11)

      if (movieseg) then
        open(93,file='movieseg',status='old')
          read(93,40)dummyi
          do 360 imem=1,1
            read(93,45)dummyf
             do 375 i=1,tgnres
                 read(93,20)
     * memsegcord(i,1,1,imem),memsegcord(i,2,1,imem),
     * memsegcord(i,3,1,imem),
     * memsegcord(i,1,2,imem),memsegcord(i,2,2,imem),
     * memsegcord(i,3,2,imem),
     * memsegcord(i,1,3,imem),memsegcord(i,2,3,imem),
     * memsegcord(i,3,3,imem)
375         continue !tgnres
360       continue !imem
        close(93)
      endif

        if (CA) i2=1
        if (CB) i2=2
        if (OX) i2=3
      
c       determine distances

	do 205 init1 = 1,maxres
           do 206 init2 = 1,maxres
        	tgdist(init1,init2) = 0.0
206	   continue
205 	continue

c        write(6,*)'opening rijtarg.dat'

        open(55,file='rijtarg2.dat',status='unknown')
        rgtarg=0.0
        do 200 i=1,tgnres-1
          do 210 j=i+1,tgnres
            tgdist(i,j)=sqrt( (finalcords(i,1,i2)
     *                      -finalcords(j,1,i2))**2
     *                    + (finalcords(i,2,i2)
     *                      -finalcords(j,2,i2))**2
     *                    + (finalcords(i,3,i2)
     *                      -finalcords(j,3,i2))**2  )
            tgdist(j,i)=tgdist(i,j)
 
            if (tgdist(i,j).lt.58.0) contact=contact+1
            count=count+1
            rgtarg=rgtarg+tgdist(i,j)*tgdist(i,j)
210       continue
200     continue
c        write(6,*) contact,count
        close(55)

        rgtarg=sqrt(rgtarg)/nres
        denstarg=1.0/rgtarg**3

	do 220 imem=1, nmem !loop over movie file

        do 235 init1 = 1,maxres
             do 236 init2 = 1,maxres
               memdist(init1,init2) = 0.0
236         continue
235     continue

c         write(6,*)'memdist loop'
          rg(imem)=0.0 
          do 230 i=1,tgnres-1
            do 240 j=i+1,tgnres
              memdist(i,j) = sqrt( (memcord(i,1,i2,imem)
     *                      -memcord(j,1,i2,imem))**2
     *                    + (memcord(i,2,i2,imem)
     *                      -memcord(j,2,i2,imem))**2
     *                    + (memcord(i,3,i2,imem)
     *                      -memcord(j,3,i2,imem))**2)
              memdist(j,i)=memdist(i,j)

            if (imem .eq. 1)then
             memsegdist(i,j) = sqrt( (memsegcord(i,1,i2,imem)
     *                      -memsegcord(j,1,i2,imem))**2
     *                    + (memsegcord(i,2,i2,imem)
     *                      -memsegcord(j,2,i2,imem))**2
     *                    + (memsegcord(i,3,i2,imem)
     *                      -memsegcord(j,3,i2,imem))**2  )
             memsegdist(j,i)=memsegdist(i,j)
            endif

              rg(imem)=rg(imem)+memdist(i,j)*memdist(i,j)

240         continue
230       continue
          rg(imem)=sqrt(rg(imem))/nres
          dens(imem)=1.0/rg(imem)**3

c         write(6,*)'call rama'
        if (rama) then
           delphi=0
           delpsi=0
           qtot_phi=0
           qtot_psi=0
       call calc_rama(tgnres,finalcords,finalphi,finalpsi)
        do 940 i_res = 1,  tgnres
          do 100 icoord=1,3
            targcords(i_res,icoord,1) = memcord(i_res,icoord,1,imem)
100       continue
          do 101 icoord=1,3
             targcords(i_res,icoord,2) = memcord(i_res,icoord,2,imem)
101	  continue
          do 102 icoord=1,3
             targcords(i_res,icoord,3) = memcord(i_res,icoord,3,imem)
102	  continue  
940     continue
       call calc_rama(tgnres,targcords,targphi,targpsi)
         do 112  i_targ=2, tgnres-1
           delphi=(finalphi(i_targ)-targphi(i_targ))*phipsiwidth
           delpsi=(finalpsi(i_targ)-targpsi(i_targ))*phipsiwidth
           qtot_phi=qtot_phi+exp(-delphi*delphi*0.5)
           qtot_psi=qtot_psi+exp(-delpsi*delpsi*0.5)
112      continue  
           q_phi(imem)=qtot_phi/(tgnres-2)
           q_psi(imem)=qtot_psi/(tgnres-2)
        endif ! rama

c---------------------------------------------------------
c       calculate width of gaussian for q-scores
c---------------------------------------------------------
          do 260 iskip=1,maxres
            del=real(iskip)
            width(iskip)=1.0/(width0*(del)**widthexp)
260       continue

c---------------------------------------------------------
c      determine qpos 
c---------------------------------------------------------

        qpos=tgnres*(tgnres-1.0)/2.0-tgnres+1.0
        nres = real(numconst)

c----------------------------------------------------------
c	determine qscore
c----------------------------------------------------------
c   first calculate the dels
c   note that q(imem) calculated here is from the
c   original code, and should be the same as the
c   new qglobal_long with power=1

         qtot=0.0
         qsegtot=0.0
         qtot_short=0.0
         qtot_medium=0.0
         qtot_long=0.0
         qpos_short=0.0
         qpos_medium=0.0
         qpos_long=0.0
         qtotlocal_short=0.0
         qtotlocal_medium=0.0
         qtotlocal_long=0.0
         qposlocal_short=0.0
         qposlocal_medium=0.0
         qposlocal_long=0.0
         qtot_cut=0.0
         qtot_short_cut=0.0
         qtot_medium_cut=0.0
         qtot_long_cut=0.0
         qpos_cut=0
         qpos_short_cut=0
         qpos_medium_cut=0
         qpos_long_cut=0

       do 270 k = 1,  numconst - 2
               i = seglist(k)
        do 280 l = k+2, numconst
               j = seglist(l)
             del=(tgdist(i,j)-memdist(i,j))*width(j-i)
             qtot=qtot+exp(-del*del*0.5)

             delseg = (tgdist(i,j)-memsegdist(i,j))
     *             *width(j-i)
             qsegtot=qsegtot+exp(-delseg*delseg*0.5)

c$$$        if (localq) then
c$$$             m = j - 2
c$$$
c$$$            if (m-i.lt.1) then
c$$$               qtotlocal_short = qtotlocal_short+exp(-del*del*0.5)
c$$$               qposlocal_short = qposlocal_short+1.0
c$$$             elseif (m-i.lt.2) then
c$$$               qtotlocal_medium = qtotlocal_medium+exp(-del*del*0.5)
c$$$               qposlocal_medium = qposlocal_medium+1.0
c$$$             elseif (m-i.lt.3) then
c$$$               qtotlocal_long = qtotlocal_long + exp(-del*del*0.5)
c$$$               qposlocal_long = qposlocal_long + 1.0
c$$$             endif
c$$$         endif

             if (tgdist(i,j).le.rcut) then
                 qtot_cut=qtot_cut+exp(-del*del*0.5)
                 qpos_cut=qpos_cut+1.0
             if (j-i.lt.5) then
               qtot_short_cut=qtot_short_cut+exp(-del*del*0.5)
               qpos_short_cut=qpos_short_cut+1.0
             elseif (j-i.lt.13) then
               qtot_medium_cut=qtot_medium_cut+exp(-del*del*0.5)
               qpos_medium_cut=qpos_medium_cut+1.0
             else
               qtot_long_cut=qtot_long_cut+exp(-del*del*0.5)
               qpos_long_cut=qpos_long_cut+1.0
             endif

             endif

             if (j-i.lt.5) then
               qtot_short=qtot_short+exp(-del*del*0.5)
               qpos_short=qpos_short+1.0
             elseif (j-i.lt.13) then
               qtot_medium=qtot_medium+exp(-del*del*0.5)
               qpos_medium=qpos_medium+1.0
             else
               qtot_long=qtot_long+exp(-del*del*0.5)
               qpos_long=qpos_long+1.0
             endif
280        continue
270      continue

c        write(6,*)'normalize by proximity class'
c        write(6,*)'qtot ',qtot
c        write(6,*)'qpos ',qpos
c        write(6,*)'qtot_short ',qtot_short
c        write(6,*)'qpos_short ',qpos_short
c        write(6,*)'qtott_medium ',qtot_medium
c        write(6,*)'qpos_medium ',qpos_medium
c        write(6,*)'qtot_long ',qtot_long
c        write(6,*)'qpos_long ',qpos_long
c        write(6,*)'imem ',imem

         qseg=qsegtot/qpos
         q(imem)=qtot/qpos
         q_short(imem)=qtot_short/qpos_short
         q_medium(imem)=qtot_medium/qpos_medium
         q_long(imem)=qtot_long/qpos_long
         q_cut(imem)=qtot_cut/qpos_cut
         q_short_cut(imem)=qtot_short_cut/qpos_short_cut
         q_medium_cut(imem)=qtot_medium_cut/qpos_medium_cut
         q_long_cut(imem)=qtot_long_cut/qpos_long_cut

220	continue   ! movie loop

c           write(6,*)'opening qphipsi '
        if (rama) then
           open(14,file='qphipsi.plot',status='unknown')
               do imem=1,nmem 
	          write(14,*)ntemp(imem),q_phi(imem),q_psi(imem)
               enddo
	   close(14)
        endif ! rama

         do imem=1,nmem 
	  write(13,4023)ntemp(imem),q(imem),q_short(imem),
     *           q_medium(imem),q_long(imem)! ,pe(imem),rg(imem)
          write(33,4023)ntemp(imem),q_cut(imem),q_short_cut(imem),
     *     q_medium_cut(imem),q_long_cut(imem)! ,pe(imem),rg(imem)
4023        format(5(f8.5,1x),(f8.1,1x),f8.5)

         enddo

	close(13)
        close(33)

          best_q = 0.0 
          best_rg = 100.0 
          ave_q_temp = 0.0
          best_temp = 0.0
          best_q_cut = 0.0
          best_temp_cut = 0.0

         do imem=1,nmem
            ave_q_temp = ave_q_temp + q(imem)
         if (q(imem) .gt. best_q) then
            best_q = q(imem)
            best_q_id = imem
            best_temp = ntemp(imem)
         endif ! if (q(imem) .gt. best_q) then

c         if (rg(imem) .lt. best_rg) then
c            write(6,*)best_rg, rg(imem), nmem, imem
c            best_rg = rg(imem)
c            best_rg_id = imem
c            best_temp_rg = ntemp(imem)
c         endif ! if (rg(imem) .lt. best_rg) then

         if (q_cut(imem) .gt. best_q_cut) then
            best_q_cut = q_cut(imem)
            best_q_id_cut = imem
            best_temp_cut = ntemp(imem)
         endif 
            ave_q = ave_q_temp / real(nmem)
         enddo ! do imem=1,nmem

         do imem=1,nmem
             if (rg(imem) .lt. best_rg) then
               best_rg = rg(imem)
               best_rg_id = imem
               best_temp_rg = ntemp(imem)
             endif 
c           write(6,*)best_rg, rg(imem), nmem, imem
         enddo ! do imem=1,nmem


         write(6,7878)best_q_id, best_q, best_temp
         write(6,7879)best_q_id_cut, best_q_cut, best_temp_cut
         write(6,7880)best_rg_id, best_rg, best_temp_rg

7878     format('best wolynes i.d q temp',i5,3x,f8.4,3x,f8.4)
7879     format('best onuchic i.d q temp',i5,3x,f8.4,3x,f8.4)
7880     format('best rg      i.d q temp',i5,3x,f8.4,3x,f8.4)

         write(16,7878)best_q_id, best_q, best_temp
         write(16,7879)best_q_id_cut, best_q_cut, best_temp_cut
         write(16,7880)best_rg_id, best_rg, best_temp_rg

         open(13,file='qscore_w.plot',status='unknown')
         open(33,file='qscore_o.plot',status='unknown')

        open(14,file='pdbq_wolynes_best.ent',status='unknown')
        open(15,file='pdbq_onuchic_best.ent',status='unknown')

         do i_test=1,2

              if (i_test .eq. 1) imem =  best_q_id
              if (i_test .eq. 2) imem =  best_q_id_cut

!         do i_axis=1,3
!               cm(i_axis)=0.0
!            do ires=1,tgnres
!                cm(i_axis)=cm(i_axis)+memcord(ires,i_axis,1,imem)
!            enddo
!               cm(i_axis)=cm(i_axis)/real(tgnres)
!          enddo
      
!         do i_axis=1,3
!            do i_atom=1,3
!                do ires=1,tgnres
!              memcord(ires,i_axis,i_atom,imem)=
!     *                     memcord(ires,i_axis,i_atom,imem)-cm(i_axis)
!                enddo
!             enddo
!         enddo

                do ires = 1+1,tgnres      ! Nitrogens    
                memcord(ires,1,4,imem) =               
     *         0.4831806*memcord(ires-1,1,1,imem) +   
     *         0.7032820*memcord(ires,1,1,imem) -     
     *         0.1864626*memcord(ires-1,1,3,imem)
                
                memcord(ires,2,4,imem) =               
     *         0.4831806*memcord(ires-1,2,1,imem) +   
     *         0.7032820*memcord(ires,2,1,imem) -     
     *         0.1864626*memcord(ires-1,2,3,imem)                                      
                memcord(ires,3,4,imem) =               
     *         0.4831806*memcord(ires-1,3,1,imem) +   
     *         0.7032820*memcord(ires,3,1,imem) -     
     *         0.1864626*memcord(ires-1,3,3,imem)
                enddo
                                                                       
                do ires = 1,tgnres-1
                
                memcord(ires,1,5,imem) =   
     *         0.4436538*memcord(ires,1,1,imem) +     
     *         0.2352006*memcord(ires+1,1,1,imem) +   
     *         0.3211455*memcord(ires,1,3,imem)

                memcord(ires,2,5,imem) =               
     *         0.4436538*memcord(ires,2,1,imem) +     
     *         0.2352006*memcord(ires+1,2,1,imem) +   
     *         0.3211455*memcord(ires,2,3,imem)
                
                memcord(ires,3,5,imem) =               
     *         0.4436538*memcord(ires,3,1,imem) +    
     *         0.2352006*memcord(ires+1,3,1,imem) +   
     *         0.3211455*memcord(ires,3,3,imem)
                enddo

         enddo ! do i_test

        atom_id = 0 

        do ires = 1,tgnres

        res_type = aaname(tgires(ires))
        atom_id   = atom_id + 1
c        write(6,*)'atom_id ',atom_id
        if (ires .ne. 1) then
        atom_type='N'
        id = 4

        write(14,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id),
     *  memcord(ires,2,id,best_q_id),
     *  memcord(ires,3,id,best_q_id),atom_id

        write(15,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id_cut),
     *  memcord(ires,2,id,best_q_id_cut),
     *  memcord(ires,3,id,best_q_id_cut),atom_id

60      format('ATOM',4x,I3,2x,a2,2x,a3,3x,I3,4x,f8.3,f8.3,f8.3, 
     *          2x,'1.00',2x,'0.00',6x,'TPDB',1x,I3)
        endif

c        if (ires .ne. 1) then
        atom_type='CA'
        id = 1

        write(14,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id),
     *  memcord(ires,2,id,best_q_id),
     *  memcord(ires,3,id,best_q_id),atom_id

        write(15,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id_cut),
     *  memcord(ires,2,id,best_q_id_cut),
     *  memcord(ires,3,id,best_q_id_cut),atom_id

c        indif

        if (ires .ne. tgnres) then
           atom_type='C'
           id = 5

        write(14,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id),
     *  memcord(ires,2,id,best_q_id),
     *  memcord(ires,3,id,best_q_id),atom_id

        write(15,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id_cut),
     *  memcord(ires,2,id,best_q_id_cut),
     *  memcord(ires,3,id,best_q_id_cut),atom_id

        endif

c        if (res_type .ne. 'GLY') then
            atom_type='CB'
            id = 2
        write(14,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id),
     *  memcord(ires,2,id,best_q_id),
     *  memcord(ires,3,id,best_q_id),atom_id

        write(15,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id_cut),
     *  memcord(ires,2,id,best_q_id_cut),
     *  memcord(ires,3,id,best_q_id_cut),atom_id

c        if (ires .ne. tgnres) then
           atom_type='O'
           id = 3

        write(14,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id),
     *  memcord(ires,2,id,best_q_id),
     *  memcord(ires,3,id,best_q_id),atom_id

        write(15,60)atom_id,atom_type,res_type,ires,
     *  memcord(ires,1,id,best_q_id_cut),
     *  memcord(ires,2,id,best_q_id_cut),
     *  memcord(ires,3,id,best_q_id_cut),atom_id
c        endif

        enddo ! do ires = 1,tgnres

        close(14)
        close(15)


c       endif ! .not. relativeq 

        if ( movieseg ) then
         open(13,file='qscorevsqmoverg.plot',status='unknown')
          do imem=1,nmem
           write(13,*)q(imem),qseg
          enddo
         close(13)
        endif ! movieseg

c!        if ( relativeq ) then
c!         open(13,file='relativeqscore.plot',status='unknown')

c!           do imem=1,nmem 
c!	     write(13,*)ntemp(imem),q(imem),q_short(imem),
c     *           q_medium(imem),q_long(imem)!,pe(imem),rg(imem
c           enddo
c	 close(13)
c        endif !  relativeq 

        end
