	program Rg

c---------------------------------------------------------------
c	This program will find the  radius of gyration 
c        f77 Rg.f -o Rg.out
c        ifort -132 Rg.f -o Rg.out
c--------------------------------------------------------------

	integer maxsiz,maxmem,maxtg

	parameter (maxsiz=5000)
	parameter (maxmem=8000)
	parameter (maxtg=5)

	integer itg,ntg,imem,nmem(maxtg),tgnres(maxtg),
     *		tgires(maxsiz,maxtg),icoord,memnres(maxmem),
     *		memires(maxsiz,maxtg),match(maxsiz,maxmem),
     *		nmode,ialigin,jalign,i,j,ires,nmovies,i_axis,
     *          i_atom,imov,zscor(maxsiz,maxmem),nat_zscor

	real coord(maxsiz),tgcord(maxsiz,3,2,maxtg),del,
     *	     memcord(maxsiz,3,3,maxmem),tgdist(maxsiz,maxsiz),
     *	     memdist(maxsiz,maxsiz),width0,widthexp,qpos,
     *	     qtot,q(maxmem,maxtg),width(maxsiz),qs(maxsiz),x,y,z,
     *	     ntemp(maxmem),rnorm,radg(maxmem,maxtg),centm(maxsiz),
     *       cm(maxmem,maxtg)

	real predicted_norm, radg_native, nat_radg, rij

	character tgname(maxtg)*5,memname(maxmem,maxtg)*10,
     *		  pheader*16,header*43,signal*4,atom*2,res*3

	logical CA,CB,OX

        data width0,widthexp / 1.0, 0.15 /

	external sset

	pheader='/home/mprentis/amh/proteins/'
	header='/usr1/kristin/programs/muself/caseAB/align/'

c---------------------------------------------------------------
c	Read in the target protein and memories
c--------------------------------------------------------------

	open(12,file='targlist_rg',status='old')
	  read(12,*)ntg
c          write(6,*)ntg
	  do 100 itg=1,ntg
	    read(12,5)tgname(itg)
c            write(6,*)tgname(itg)
5	    format(a5)
100	  continue
	close(12)

       itg=1

       open(12,file='rg_input',status='old')
          read(12,*)tgnres(1)
c          write(6,*)tgnres(1) 
          read(12,*)CA
c          write(6,*)CA
          read(12,*)CB
c          write(6,*)CB
          read(12,*)OX
c          write(6,*)OX
       close(12)

c	CA=.true.
c	CB=.false.

c--------------------------------------------------------------
c	Open protein files
c--------------------------------------------------------------

        rnorm=1.0/(tgnres(1))**2

	open(12,file='/home/mp466/amh/proteins/'//tgname(itg),status='old')
          read (12,*)
          read (12,*) tgnres(itg)
	  predicted_norm = 2.2*(real(tgnres(itg)))**0.38 
         read (12,1022) (tgires(i1,itg),
     *            i1=1,tgnres(itg))
1022	  format(25(i2,1x))
          do 130 icoord=1,3
            read (12,1023)
     *           (tgcord(i1,icoord,1,itg),i1=1,tgnres(itg))
1023        format(8(f8.3,1x))
130       continue
          do 140 icoord=1,3
            read (12,1023)
     *           (tgcord(i1,icoord,2,itg),i1=1,tgnres(itg))
140       continue
          do 150 icoord=1,3
            read (12,1023)
     *           (coord(i1),i1=1,tgnres(itg))
150       continue
        close (12)

        nat_radg=0.0
        nat_zscor=0
        do i = 1,tgnres(itg)-1
        do j = i+1,tgnres(itg)

        rij = 0.0
        rij = (tgcord(j,1,1,itg)-tgcord(i,1,1,itg))**2
     *  + (tgcord(j,2,1,itg)-tgcord(i,2,1,itg))**2
     *  + (tgcord(j,3,1,itg)-tgcord(i,3,1,itg))**2
        nat_radg=nat_radg+rij
        rij = sqrt(rij)

        if  (rij .lt. 8.0)nat_zscor=nat_zscor+2.0

        enddo
        enddo

         nat_radg=nat_radg*rnorm     
c        write(6,*)'nat_radg, rnorm, tgnres', radg(imem,itg),tgnres(1)
         nat_radg = sqrt(nat_radg)
         write(6,*)' asdf ', nat_zscor, nat_radg

c---------------------------------------------------------
c  	read in coordinates of protein during folding run
c 	from movie file
c---------------------------------------------------------

	open(11,file='movie',status='old')

	  read(11,40)nmovies
c          write(6,*)'nmovies', nmovies
40	  format(30x,i5)
	  do 160 imem=1,nmovies
	    nmem(itg)=nmem(itg)+1
	    read(11,45)ntemp(imem)
c            write(6,*)'ntemp', ntemp(imem)
45	    format(21x,f7.4)
c            do 165 i=1,tgnres(itg)
c	      read(11,*)
c165	    continue
c	    read(11,45)ntemp(imem)
	    do 175 i=1,tgnres(itg)
              if (CA) read(11,20)memcord(i,1,1,imem),
     *		memcord(i,2,1,imem),memcord(i,3,1,imem)
20     format('CA: ',3(f8.3,1x),'CB: ',3(f8.3,1x),'Ox: ', 3(f8.3,1x))

 
c              write(6,20)memcord(i,1,1,imem),
c    *          memcord(i,2,1,imem),memcord(i,3,1,imem)

	      if (CB) read(11,22)memcord(i,1,2,imem),
     *		memcord(i,2,2,imem),memcord(i,3,2,imem)
22	      format(34x,3(1x,f8.4))
175         continue
160	  continue

	close(11)

c----------------------------------------------------------------
c   CALCULATE RG OF EACH snapshot 
c-----------------------------------------------------------------

        do 400 imem = 1,nmovies

       do i_axis=1,3
          cm(i_axis,imov)=0.0
         do ires=1,tgnres(itg)
              cm(i_axis,imov)=
     *           cm(i_axis,imov)+memcord(ires,i_axis,1,imov)
         enddo
              cm(i_axis,imov)=cm(i_axis,imov)/real(tgnres(itg))
       enddo

         do i_axis=1,3
           do i_atom=1,3
             do ires=1,tgnres(itg)
           memcord(ires,i_axis,i_atom,imov)=
     *              memcord(ires,i_axis,i_atom,imov)-cm(i_axis,imov)
             enddo
           enddo
         enddo
400      continue

	do 300 imem = 1,nmovies

	radg(imem,itg)=0.0
        zscor(imem,itg)=0 
	do i = 1,tgnres(itg)-1 
	do j = i+1,tgnres(itg)

	rij = 0.0
	rij =  (memcord(j,1,1,imem)-memcord(i,1,1,imem))**2
     *  + (memcord(j,2,1,imem)-memcord(i,2,1,imem))**2 
     *  + (memcord(j,3,1,imem)-memcord(i,3,1,imem))**2 
c        write(6,*)rij
	radg(imem,itg)=radg(imem,itg)+rij
        rij = sqrt(rij)

        if  (rij .lt. 8.0)zscor(imem,itg)=zscor(imem,itg)+2.0

	enddo
	enddo

	radg(imem,itg)=radg(imem,itg)*rnorm	
c        write(6,*)'radg, rnorm, tgnres', radg(imem,itg),tgnres(1)
	radg(imem,itg)=sqrt(radg(imem,itg))
300      continue

        open(13,file='rg.plot',status='unknown')
	     
	do 290 imem=1, nmovies
          write(13,25)ntemp(imem),radg(imem,itg),radg(imem,itg)/nat_radg
25	  format(f10.5,10x,f10.5,10x,f10.5)
290	continue
	close(13)

        open(13,file='num_contacts.plot',status='unknown')
        do 390 imem=1, nmovies
          write(13,26)ntemp(imem),zscor(imem,itg),zscor(imem,itg)/real(nat_zscor)
26        format(f10.5,10x,i5,10x,f10.5)
390     continue
        close(13)

c	print *, radg_native

	stop
	end

