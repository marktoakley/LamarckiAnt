
c     --------------------- hprofl ----------------------

      subroutine hprofl(maxsiz,nmres,ires,hydseq,hydscl)

c     ---------------------------------------------------

c     HPROFL set hydrophobicity profile for coded amino
c            acid sequence in the array ires

c     arguments:

c        maxsiz- maximum protein length (i)
c        nmres - length of current protein (i)
c        ires  - coded amino acid sequence (i)
c        hydseq- hydrophobicity profile (o)
c        hydscl- hydrophobicity scale (i)

c     ---------------------------------------------------

c     set required parameters

      implicit none

c     argument declarations:

         integer maxsiz,nmres,ires(maxsiz)

         double precision hydseq(maxsiz),hydscl(0:21)

c     internal variables:

         integer i500,i501
     
c     --------------------- begin -----------------------

c     set hydrophobicity profile

c     put h value from array hydscl in 'utility'

      do 500 i500=1,nmres
         hydseq(i500)=hydscl(ires(i500))
  500 continue

c     translate h values to binary scale

      do 501 i501=1,nmres
         if( hydseq(i501).lt.0.0D0 )then
            hydseq(i501)=-1.0D0
         elseif( hydseq(i501).gt.0.0D0 )then
            hydseq(i501)=1.0D0
         endif
  501 continue

c     ---------------------- done -----------------------

      return
      end
