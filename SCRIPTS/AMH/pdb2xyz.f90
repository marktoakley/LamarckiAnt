	program pdb2xyz

!   This program will convert a pdb file into a protein file that can be
!   read by the Wolynes group's programs.  It will convert the AA sequence
!   from three letter representation to numerical representation as well
!   as record the CA, CB and O coordinates

!   INPUT  file targlist : where the first line is the number of proteins
!                          converted, the number of residues, 
!                          followed by the name of the protein
!               example  : 856
!                          56
!                          1
!                          3GB1 
!   OUTPUT file          : 3GB1ZZZ

!   ifort pdb2wolynes.f90 -o pdb2wolynes.out

	integer maxsiz,maxtg

	parameter (maxsiz=1600)
	parameter (maxtg=500)

	integer i,j,isig,ires,nres,ntgres(maxsiz),first,last,i1,i2,i3, &
                numofres,numoflines, numofgly,tgsurf(maxsiz), &
                tgks(maxsiz),tempval,numcrd,numpro,nmdifv

	real x,y,z,tgcoord(maxsiz,3,3)

	character filname(maxtg)*20,tgfile(maxtg)*5,signal*4,res*3, &
     		  tgres(maxsiz)*3,atom*2,name*60,sub*1,set(maxtg)*1, &
                   tempfile*20, pdbfilename(maxtg)*100,temptempfile*20

!    The name of the pdb file and the name of the new file will be read
!    in interactively
           numcrd=3
           numpro=1
           nmdifv=1
           numofgly = 0
         open(11,file='targlist',status='old')
            read(11,*)nlines
            write(6,*)nlines
            read(11,*)nres
            write(6,*)nres
            read(11,*)ntg
            write(6,*)ntg

	  do  itg=1,ntg
	    read(11,211)pdbfilename(itg)
             tempfile=trim(pdbfilename(itg))
!              tempval=len_trim(tempfile)
!              temptempfile=tempfile(1:5)
!            write(6,33)trim(tempfile)
33          format('tempfile',a10)
!            write(6,*)'temptempfile ',temptempfile
211         format(a100)
          enddo

!---------------------------------------------------------------------------
!   The pdb file will be opened and the sequence and coordinates will be
!   read in
!---------------------------------------------------------------------------
        open(10,file=tempfile,status='old')
!	  read(10,*)
!          read(10,6) name
!6         format(10x,a60)
          do 100 isig=1,10000
            read(10,15)signal
15          format(a4)
            if (signal.eq.'ATOM') goto 110
100       continue
110       backspace(10)

!---------------------------------------------------------------------------
!   Now the CA and CB coordinates are read in for whichever subunit was
!   picked
!--------------------------------------------------------------------------

!          read(10,17)sub
!17        format(21x,1a)
!          if (sub.eq.set(itg)) goto 111
!          do 135 i1=1,100000
!	    read(10,17)sub
!	    if (sub.eq.set(itg)) goto 111
!135	  continue

!111	  backspace(10)	
	  last=0
	  do 200 i1=1,maxsiz
	    tgres(i1)='   '
	    do 210 i2=1,3
	      do 220 i3=1,3
	 	tgcoord(i1,i2,i3)=0.0
220	      continue
210	    continue
                tgks(i1)=0
                tgsurf(i1)=0
200	  continue
	  
          do 120 i=1, nlines
            read(10,20)signal,atom,res,ires,x,y,z
20          format(a4,9x,a2,2x,a3,3x,i3,4x,3(f8.3))
!            write(6,20)signal,atom,res,ires,x,y,z
            if (signal.eq.'TER '.or.signal.eq.'HETA'.or. &
     		  signal.eq.'END ') goto 130

               if (res.eq.'GLY') then
                   numofgly = numofgly +1
                endif

            if (atom.eq.'CA') then
!	      if (ires.ne.last+1)
!     *	   	write(6,*)'nonconsecutive residues: ',last,ires,
!     *          filname(itg)
!                last=ires
                tgres(ires)=res
                tgcoord(ires,1,1)=x
                tgcoord(ires,2,1)=y
                tgcoord(ires,3,1)=z

                if (tgres(ires).eq.'GLY') then
                  tgcoord(ires,1,2)=tgcoord(ires,1,1)
                  tgcoord(ires,2,2)=tgcoord(ires,2,1)
                  tgcoord(ires,3,2)=tgcoord(ires,3,1)
                endif
              endif

              if (atom.eq.'CB') then
                 tgcoord(ires,1,2)=x
                 tgcoord(ires,2,2)=y
                 tgcoord(ires,3,2)=z
              endif

              if (i.eq. numoflines) then
                  tgcoord(ires,1,3)=0.0
                  tgcoord(ires,2,3)=0.0
                  tgcoord(ires,3,3)=0.0
                  goto 130 
              endif

	      if (atom.eq.'O ') then
		  tgcoord(ires,1,3)=x
		  tgcoord(ires,2,3)=y
		  tgcoord(ires,3,3)=z
	      endif
120        continue
130        continue
!	   nres=last-first+1

        close(10)

!----------------------------------------------------------------------------
!   Now the new protein file will be written out
!---------------------------------------------------------------------------

	open(12,file=trim(tempfile)//'XYZ',status='unknown')   
            DO ires=1,nres
                WRITE(12,'(3F12.3)')tgcoord(ires,1,1), tgcoord(ires,2,1), tgcoord(ires,3,1)
                WRITE(12,'(3F12.3)')tgcoord(ires,1,2), tgcoord(ires,2,2), tgcoord(ires,3,2)
                WRITE(12,'(3F12.3)')tgcoord(ires,1,3), tgcoord(ires,2,3), tgcoord(ires,3,3)
            ENDDO
	close(12)

        open(66,file='movieseg'//trim(tempfile),status='new')

        write(66,334)nres,numcrd,numpro,nmdifv
334     format(4(i8,1x),' nmres nmcrd numpro nmsnap')
        write(66,683)1,1,1,1.0,1
683     format(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')

        do ires = 1,nres

        write(66,50)tgcoord(ires,1,1),tgcoord(ires,2,1),tgcoord(ires,3,1),&
                    tgcoord(ires,1,2),tgcoord(ires,2,2),tgcoord(ires,3,2),&
                    tgcoord(ires,1,3),tgcoord(ires,2,3),tgcoord(ires,3,3)
50       format('CA: ',3(f8.3,1x),'CB: ',3(f8.3,1x),'Ox: ', 3(f8.3,1x))

        enddo

        close(66)




95	continue
        close(11)

	stop
	end
