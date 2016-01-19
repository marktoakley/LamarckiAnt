
c     --------------------- getmem ----------------------

      subroutine getmem(protnm,nmrss,maxres,ires,imemri,
     *                  numcrd,ywork,secseq,
     *                  oarchv,sa,
     *                  iwork,passi, 
     *                  i_mem)

c     ---------------------------------------------------

c     GETMEM gets the coordinates for the next memory
c            protein and sets the hydrophobicity file

c     arguments:

c 	 protnm - name of protein in pdb code (i)
c        nmrss - length of input protein (o)
c        maxres- maximum database protein length (i)
c        ires  - coded amino acid sequence (o)
c        imemri- tape number for coordinate file (i)
c        numcrd- number of coordinate types (i)
c        ywork - coordinates (o)
c        secseq - secondary structure from protein file (o)
c        oarchv- output/diagnostic file number (i)

c     ---------------------------------------------------
cError 355 : In program unit GETMEM variable TARFL has not been given a type

      use amhglobals,  only: known,imem_cons,mem_cons

      implicit none


c     set required parameters

c     if echo, then echo input

      logical echo
      parameter(echo=.false.)

c     rewrite database excluding 'bad' proteins

      logical sortpb
      parameter(sortpb=.false.)

c     argument declarations:
c1234567890123456789012345678901234567890123456789012345678901234567890
	 character protnm*5,confile*41

         logical passi

         integer nmrss,maxres, ires(maxres), imemri,
     *           oarchv,numcrd, secseq(maxres), sa(maxres),
     *           iwork(maxres), i_mem, nmrsstemp,
     *           irestemp(maxres)

         integer num_coord_to_read,open_status

         double precision ywork(maxres,3,numcrd)

c     internal variables:


         integer itrack

c        --- do loop indices ---

         integer i_coord,i_res,i_axis

c        --- implied do loop indices ---

         integer i1,i517

c------------------- begin -----------------------
c      write(6,*)'mem_cons imem_cons',mem_cons, imem_cons 

      if(( .not. mem_cons) .or. ( i_mem.eq.0 ))then
c      write(6,*)'mem_cons false read native seq '
      read(imemri,*)
      read(imemri,101)nmrss
  101 format(i5)

      if( echo )write(oarchv,101)nmrss

c     check enough storage has been reserved
      if( nmrss.gt.maxres )then
         write(oarchv,110)protnm,nmrss,maxres
  110    format('Getmem: protein ',a5,' too large ',i4,
     *          ' reserved space ',i4)
         stop
      endif

c     read in coded primary sequence
      read(imemri,103)(ires(i1),i1=1,nmrss)
  103 format(25(i2,1x))
c      write(6,103)(ires(i1),i1=1,nmrss)
        if( echo )write(oarchv,103)(ires(i1),i1=1,nmrss)

       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read in consensus sequences
      if((mem_cons) .and.  (i_mem.ne.0))then

        confile =
     *       '/home/mprentis/amh/md_input/memcons/'//protnm
c           123456789012345678901234567890123456789012345678

c         write(6,*)'mem_cons directory = ' , confile
         open(imem_cons,file=confile,status='old',iostat=open_status)
               if (open_status.ne.0) then
                 write(6,*) 'failure to open file ',confile
                 write(6,*) 'error number ',open_status
                 stop
               endif

        read(imem_cons,*)
        read(imem_cons,101)nmrss
      
       if( echo )write(oarchv,101)nmrss
c     check enough storage has been reserved
                                                                                
      if( nmrss.gt.maxres )then
         write(oarchv,110)protnm,nmrss,maxres
         stop
      endif
                                                                                
      read(imem_cons,103)(ires(i1),i1=1,nmrss)
c      write(6,103)(ires(i1),i1=1,nmrss)
      if( echo )write(oarchv,103)(ires(i1),i1=1,nmrss)
      close(imem_cons)

      endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     check that residues fall within correct range

      if( passi )then
         itrack=0
c         do 517 i517=1,nmrss
c            if( (ires(i517).lt.1).or.(ires(i517).gt.20) )then
c               itrack=itrack + 1
c               iwork(itrack)=i517
c            endif
c  517    continue
 
c        print message if there is 'bad' residue

         if( itrack.gt.0  )then
            write(oarchv,322)protnm
  322       format(/'Getmem: Warning ',a5,' has invalid',
     *             ' residue type')
            write(oarchv,355)(iwork(i1),i1=1,itrack)
  355       format(20(i3,1x))
         elseif( sortpb )then

            write(42,150)protnm
  150       format(a5,' zzz')
            write(42,151)nmrss
  151       format(i5,' xxx')
            write(42,103)(ires(i1),i1=1,nmrss)
         endif
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read in protein coordinates 

        if (i_mem .ne. 0 .or.  known ) then
              num_coord_to_read = 3

c  read but ignore

      if((mem_cons) .and. (i_mem .ne. 0))then

      read(imemri,*)
      read(imemri,901)nmrsstemp
901     format(i5)
      read(imemri,903)(irestemp(i1),i1=1,nmrss)
903     format(25(i2,1x))

      endif

      do 500 i_coord=1, num_coord_to_read 
        do 501 i_axis=1,3

          read(imemri,104)(ywork(i_res,i_axis,i_coord),i_res=1,nmrss)

               if( sortpb.and.(itrack.eq.0) )
     *    write(42,104)(ywork(i_res,i_axis,i_coord),i_res=1,nmrss)
  104          format(8(f8.3,1x))

  501    continue
  500  continue

       read(imemri,103)(sa(i1),i1=1,nmrss)
       do 600 i1=1,nmrss
	 sa(i1)=2-sa(i1)
600    continue
       read(imemri,103)(secseq(i1),i1=1,nmrss)

        endif

c     ---------------------- done -----------------------

      return
      end
