
c     --------------------- pttarg ----------------------
 
      subroutine pttarg(maxsiz,nmres,numcrd,target,
     *                  maxres,ywork,ires,jres,shydro,
     *                  hydseq,tarpre,mempre,bondln,
     *                  oarchv,passi)

c     ---------------------------------------------------

c     PTTARG copies the target coordinates, coded 
c            sequence, hydrophobicity profile
c            secondary structure predictions into
c            arrays reserved for target; the nearest-
c            neighbor distances between alpha-carbons
c            and the distance between the alpha- and
c            beta-carbon of the same residue are 
c            also computed

c     arguments:

c        maxsiz- maximum number of protein residues (i)
c        nmres - actual number of residues (i)
c        numcrd- number of coordinate types (i)
c        target- target protein coordinates (o)
c        maxres- maximum protein database length (i)
c        ywork - target coordinates (i)
c        ires  - target coded sequence (o)
c        jres  - target coded sequence (i)
c        shydro- target hydrophobicity profile (o)
c        hydseq- target hydrophobicity profile (i)
c        tarpre- target secondary structure profile (o)
c        mempre- target secondary structure profile (i)
c        bondln- bond lenghts for target coordinates;
c                bond i is the bond between residues i-1
c                and i (o)
c        oarchv- archive/diagnostic output file (i)
c        passi - true on first call to pttarg (i)

c     ---------------------------------------------------

c     set required parameters

      implicit none

c     argument declarations:

         logical passi

         integer maxsiz,nmres,ires(maxsiz),jres(maxsiz),
     *           tarpre(maxsiz),mempre(maxsiz),oarchv,numcrd,maxres

         double precision target(maxsiz,3,numcrd),hydseq(maxres,2),
     *      ywork(maxres,3,numcrd),shydro(maxsiz,2),bondln(maxsiz,numcrd)

c     internal variables:

         integer indx1,indx2,indx3,indx4,isum,iisum

c        --- do loop indices ---

         integer i1,i500,i_coord,i_axis,i_res,i505,i510,i511

c     --------------------- begin -----------------------

c     --- diagnostics ---

c     echo scalar input arguments

c      write(oarchv,100)maxsiz,nmres
c  100 format('Pttarg:maxsiz ',i3,' nmres ',i3)

c     --- end diagnostics ---

c     copy coordinates,

      do 502 i_coord=1,numcrd
         do 503 i_axis=1,3
            do 504 i_res=1,nmres
               target(i_res,i_axis,i_coord)=ywork(i_res,i_axis,i_coord)
  504       continue
  503    continue
  502 continue

c     copy coded sequence, h profile and
c     secondary structure profile

      do 500 i500=1,nmres

         ires(i500)=jres(i500)
         shydro(i500,1)=hydseq(i500,1)
         shydro(i500,2)=hydseq(i500,2)
         tarpre(i500)=mempre(i500)

  500 continue

c     print target h profile

      if( passi )then
         indx2=10
         indx1=int( float(nmres)/float(indx2) ) + 1
c         write(oarchv,712)
c  712    format(/'Hydrophobicity profile and row sum ')
         iisum=0
         do 510 i510=1,indx1
            indx4=(i510-1)*indx2 + 1
            if( indx1.eq.i510 )then
               indx3=nmres
            else
               indx3=indx4 + indx2 - 1
            endif
            isum=0
            do 511 i511=indx4,indx3
               isum=isum + int(shydro(i511,1))
  511       continue
            iisum=iisum + isum
c            write(oarchv,711)(int(shydro(i1,1)),i1=indx4,indx3),
c     *                        isum
  711       format(11(i3,1x))
  510    continue
c         write(oarchv,713)iisum
c  713    format('net h ',i3)
      endif

c     --- diagnostics ---

c     echo h profile and coded primary sequence

c      write(oarchv,107)prot,nmres
c  107 format(/'target ',a4,' residues ',i3)
c      do 504 i_res=1,nmrs
c         write(oarchv,108)i_res,hydscl(ires(i_res)),
c     *                    (shydro(i_res,i1),i1=1,2)
c  108    format(i3,3(1x,f7.2))
c  504 continue

c     --- end diagnostics ---

c      set bond lengths for alpha(i+1)-alpha(i) and
c      alpha(i)-beta(i)

      do 505 i505=1,nmres
             bondln(i505,1)=3.8004D0
             bondln(i505,2)=1.54D0
             bondln(i505,3)=2.42677D0
  505 continue
      bondln(1,1)=0.0D0

c     set glycine alpha-beta bond

      do 512 i_res=1,nmres
         if( jres(i_res).eq.8 )
     *      bondln(i_res,2)=0.0D0
  512 continue

c     --- diagnostics ---
c
c      if( passi )then
c         write(oarchv,110)
c  110    format(/'Pttarg:site      bondln',22x,
c     *           'coordinates')
c         do 509 i509=1,min(nmres,109)
c
c            write(oarchv,109)i509,bondln(i509,1),
c     *                       (target(i509,i1,1),i1=1,3)
c  109       format(/8x,i3,1x,4(1x,1pe10.3))
c
c            write(oarchv,111)bondln(i509,2),
c     *                       (target(i509,i1,2),i1=1,3)
c  111       format(12x,4(1x,1pe10.3))
c  509    continue
c      endif
c
c     --- end diagnostics ---

c     ---------------------- done -----------------------

      return
      end
