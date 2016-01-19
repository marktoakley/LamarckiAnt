
c     --------------------- potind ----------------------

      subroutine potind(maxsiz,maxcnt,numlng,ilong,
     *                  jbegn,jend,nmres,i_tab)

c     ---------------------------------------------------

c     POTIND sets up constraint list for specified
c            parameters (lstrt,lfins)

c     arguments:

c        nmres - actual number of residues

c     ---------------------------------------------------

c     set required parameters

      use amhglobals,  only: ixn_from_site,numconst_a,seglist_a,maxres,
     *                i_ixn_Qbias_b , i_ixn_Qbias_a,
     *                numconst_b, seglist_b
      implicit none

c     argument declarations:

         integer maxsiz,maxcnt,numlng(0:maxsiz),
     *           ilong(maxcnt,2),nmres,
     *           jbegn,jend,i_tab,
     *           isit1,isit2,i_temp   
  
c     internal variables:

         integer idd

c        --- do loop indices ---

         integer i_res,i,j
 
c        --- implied do loop indices ---

         integer res_a, res_b,arseflap

c     --------------------- begin -----------------------


c     set up indices for pairs of constraints utilized in
c     potential

      numlng(0)=0
      idd=0
      arseflap=0
      do 500 res_a=1,nmres-jbegn
         do 901 res_b=jbegn,jend-res_a 
            arseflap=arseflap+1

c  assignment of  ilong

            ilong(idd+res_b-jbegn+1,1)=res_a
            ilong(idd+res_b-jbegn+1,2)=res_a + res_b

            ixn_from_site(res_a,res_a+res_b,i_tab)
     *                             =idd+res_b-jbegn+1
           

  901    continue 
         idd=idd + jend - res_a - jbegn  + 1
         numlng(res_a)=idd 

  500 continue

c          write(6,*)'arseflap ',arseflap

      do 502 i_res=nmres-jbegn+1,nmres
         numlng(i_res)=numlng(nmres-jbegn)
  502 continue

c finds proper interaction from list for biasing of segments
c  creat table once

       if (i_tab .eq. 1)then

       do 30  i = 1, maxres
         do 40 j = 1, maxres
                i_ixn_Qbias_a(i,j)=0
40     continue
30     continue   

       do 50 i = 1,  numconst_a - 2
             isit1=seglist_a(i)
        do 100 j = i+2, numconst_a
                isit2 = seglist_a(j)
             do i_temp=1,numlng(nmres)
                if (isit1 .eq. ilong(i_temp,1) .and.
     *                (isit2 ).eq. ilong(i_temp,2)) then
                        i_ixn_Qbias_a(isit1,isit2) = i_temp
                endif
             enddo ! i_temp=1,numlng(nmres)
100         continue
50          continue

       endif ! if (i_tab .eq. 1)then

c finds proper interaction from list for biasing of segments
c  creat table once

      if (i_tab .eq. 1)then
       do 301 i = 1, maxres
         do 401 j = 1, maxres
             i_ixn_Qbias_b(i,j)=0
401     continue
301     continue

       do 501 i = 1,  numconst_b - 2
             isit1=seglist_b(i)
        do 1001 j = i+2, numconst_b
                isit2 = seglist_b(j)
             do i_temp=1,numlng(nmres)
                if (isit1 .eq. ilong(i_temp,1) .and.
     *                (isit2 ).eq. ilong(i_temp,2)) then
                    i_ixn_Qbias_b(isit1,isit2) = i_temp
                endif
             enddo ! i_temp=1,numlng(nmres)
1001     continue
501    continue
      endif ! if (i_tab .eq. 1)then

c     ---------------------- done -----------------------

      return
      end
