 
c     --------------------- scltab ----------------------

      subroutine scltab

c     ---------------------------------------------------

c     SCLTAB scales the potential energy according to
c            the specified energy/residue for the
c            target conformation

c     ---------------------------------------------------

      use amhglobals,  only:SO, maxtab,maxr,oarchv,i521,numcrd,i523,
     *  prcord,zrcord,avep,vpotnt,forse,maxsiz,nmres,nummem,
     *  iwork,pexcld,numlng,rinc,rincinv,iprolstscl,
     *  rincsq,ilong,crdixn,idigns,ires,iprolst,
     *  eqdist,hbond,oxexcldv,i1,numtab,srplr,eres,
     *  i502,i507,maxs,i505,i511,i512,i2,i508,i509,i516,i515,
     *  i525,target,igomb,ngomb,iclassweight,test_site,
     *  ixn_from_site,i_V_test,minmr,maxmr,ovscltab

      implicit none


c     internal variables:

         logical tempav,scl_call

         integer i,i_test_tab,i_test_ixn,i_range

         integer i_ixn, i_coord, i_mem, i_tab,count, k,j

         double precision rnorm,rmax(maxtab),rmin(maxtab),trgeng(maxtab,3),trgscl(3),srplr_tot,
     *        vpotint_test(0:maxr+1),trgengtot(3)

         character*10 ccount1,ccount2,ccount3
         character*5 tarfl,memfl

        integer Nres,memires(maxsiz),memsurf(maxsiz),memSS(maxsiz),open_status

        double precision memcord(maxsiz,3,3)
 
c     required subroutines:

         external force,num_to_char

c     --------------------- begin -----------------------
c        write(SO,*) 'scltab'c     get target's raw table potential  energy
c     find potential energy for each interaction type

         trgengtot(i_range)=0.0D0
      tempav=.true.
c ***************

c   read target (which is considered 0th memory)
          rewind iprolst

          read (iprolst,1000)tarfl
       do 609 i_mem=1,nummem
1000     format(a5)
          read (iprolst,1000)memfl

          open(17,file='proteins/'//memfl,
     *        status='old',iostat=open_status)

         if (open_status.ne.0) then
           write(SO,*) 'failure to open protein file ', memfl
           write(SO,*) 'in scale tab'
           write(SO,*) 'error number ',open_status
           stop
         endif

           read(17,*)
           read(17,*)Nres
           read(17,201)(memires(j),j=1,Nres)
201        format (25(i2,1x))
           do j = 1,3
              read(17,202)(memcord(k,j,1),k=1,Nres)
202           format(8(f8.3,1x))
           enddo
           do j = 1,3
              read(17,202)(memcord(k,j,2),k=1,Nres)
           enddo
           do j = 1,3
              read(17,202)(memcord(k,j,3),k=1,Nres)
           enddo
           read(17,201)(memsurf(k),k=1,Nres)
           read(17,201)(memSS(k),k=1,Nres)
        close(17)
 
      write(oarchv,*)

        do 521 i521=1,3                ! set protein coords to native str
          do 522 i_coord=1,numcrd
            do 523 i523=1,nmres
c               prcord(i523,i521,1,i_coord)=target(i523,i521,i_coord)
              prcord(i523,i521,1,i_coord)=memcord(i523,i521,i_coord)
 523            continue
 522      continue
 521    continue

      scl_call=.true.

         call force(1,prcord,zrcord,avep,tempav,maxr,vpotnt,forse,trgeng,pexcld,
     *     numlng,nmres,rincinv,rincsq,ilong,crdixn,ires,eqdist,hbond,oxexcldv,scl_call)

         idigns=.false.

       do i_range=1,3
         do i_tab=1,numtab
          trgengtot(i_range)=trgengtot(i_range)+trgeng(i_tab,i_range)
           write(SO,*) i_tab,i_range,trgeng(i_tab,i_range),trgengtot(i_range)
         enddo
       enddo 

609     continue

       do i_range=1,3
           write(SO,*)'before normalization ',trgengtot(i_range) 
           trgengtot(i_range)=trgengtot(i_range)/dble(nummem)
           write(SO,*)'after normalization ',trgengtot(i_range)
       enddo

       srplr_tot=srplr(1)+srplr(2)+srplr(3)
 
       if (iclassweight) then
         do i_range=1,3
           trgscl(i_range)= dabs( (srplr(i_range)/srplr_tot)*4.0D0*
     *                        eres*dble(nmres)/trgengtot(i_range) )
         enddo
       else
         do i_range=1,3
           if (igomb) then
            trgscl(i_range)=(dabs( 4.0D0*eres*dble(nmres)/avep(1,1,5)))**(1.0D0/ngomb)
           else
            trgscl(i_range)=(dabs( 4.0D0*eres*dble(nmres)/avep(1,1,5)))  
           endif
         enddo
       endif

         write(SO,*) 'short, med and long range conts to potential are'
         write(SO,*) (trgeng(i1,1),i1=1,maxtab)
         write(SO,*) (trgeng(i1,2),i1=1,maxtab)
         write(SO,*) (trgeng(i1,3),i1=1,maxtab)
         write(SO,*) 'totals per class:'
         write(SO,*) (trgengtot(i_range),i_range=1,3)
         write(SO,*) 'scaling factors' 
         write(SO,*) (trgscl(i_range),i_range=1,3)

c     weight table energies for each coordinate type

      do 504 i_tab=1,numtab
         do 502 i502=1,numlng(nmres,i_tab)  ! set iwork=j-i
            iwork(i502)=ilong(i502,2,i_tab) - ilong(i502,1,i_tab)
502      continue


c        loop over interaction pairs
         do 506 i_ixn=1,numlng(nmres,i_tab)

             if( iwork(i_ixn).lt.minmr )then

c              short-range ixns

             do 507 i507=0,maxs+1
                 vpotnt(i507,i_ixn,i_tab)=trgscl(1)*vpotnt(i507,i_ixn,i_tab)
                 forse(i507,i_ixn,i_tab)=trgscl(1)*forse(i507,i_ixn,i_tab)
507          continue

            elseif( iwork(i_ixn).gt.maxmr ) then

c              long-range ixns

               do 525 i525=0,maxs+1
                  vpotnt(i525,i_ixn,i_tab)=trgscl(3)*vpotnt(i525,i_ixn,i_tab)
                  forse(i525,i_ixn,i_tab)=trgscl(3)*forse(i525,i_ixn,i_tab)
525            continue
            else

c              med-range ixns

               do 505 i505=0,maxs+1
                  vpotnt(i505,i_ixn,i_tab)=trgscl(2)*vpotnt(i505,i_ixn,i_tab) 
                  forse(i505,i_ixn,i_tab)=trgscl(2)*forse(i505,i_ixn,i_tab)
  505          continue
            endif
  506    continue


       if (i_V_test) then

         i_test_tab=i_tab
         i_test_ixn=ixn_from_site(test_site(1),
     *                         test_site(2),i_test_tab)
         if ( (test_site(1).ne.ilong(i_test_ixn,1,i_test_tab))
     *    .or.(test_site(2).ne.ilong(i_test_ixn,2,i_test_tab)) )
     *      then
            write(SO,*) 'error in scltab'
            write(SO,*) 'test sites for table',i_test_tab,'wrong'
            write(SO,*) 'should be',test_site(1),test_site(2)
            write(SO,*) 'find from looking up interaction they are'
            write(SO,*) ilong(i_test_ixn,1,i_test_tab),
     *                    ilong(i_test_ixn,1,i_test_tab)
            goto 1111
          endif

         vpotint_test=0.0D0

         count=test_site(1)
         call num_to_char(count,ccount1)
         count=test_site(2)
         call num_to_char(count,ccount2)
         count=i_test_tab
         call num_to_char(count,ccount3)
         open(ovscltab,file='V_scltab_'//trim(ccount1)//'_'//
     *    trim(ccount2)//'.'//trim(ccount3),status='unknown',action='write')

       write(ovscltab,*) 'r,vpotnt,forse and int force for ixn,table'
     *                                 ,i_test_ixn,i_test_tab
       write(ovscltab,*) 'sites',ilong(i_test_ixn,1,i_test_tab)
     *                               ,ilong(i_test_ixn,2,i_test_tab)
       do i=maxs,0,-1 
         vpotint_test(i)=vpotint_test(i+1)+
     *                  (1.0D0*rinc(i_test_ixn,i_test_tab))*0.5D0*
     *                  (forse(i+1,i_test_ixn,i_test_tab)*dble(i+1)+
     *                   forse(i,i_test_ixn,i_test_tab)*
     *                   dble(i))*rinc(i_test_ixn,i_test_tab)
       enddo
       do i=0,maxs+1
         write(ovscltab,*) dble(i)*rinc(i_test_ixn,i_test_tab),
     *                   vpotnt(i,i_test_ixn,i_test_tab),
     *   (forse(i,i_test_ixn,i_test_tab)*rinc(i_test_ixn,i_test_tab))
     *                          *dble(i),vpotint_test(i) 
       enddo

1111  endif !end of writing out test potential
       
       
  504 continue

c     --- diagnostics --- 

      rnorm=0.0D0
      do 511 i511=1,numtab
         rnorm=rnorm + trgeng(i511,1) + trgeng(i511,2)
  511 continue
      write(oarchv,345)rnorm
  345 format(/'Target Energies and Scaling Factors ',1pe10.3)
      do 512 i512=1,numtab
         write(oarchv,130)i512,(trgeng(i512,i2),i2=1,2),
     *                    (trgscl(i2),i2=1,2)
  130    format('Table ',i2,1x,4(1pe10.3,1x))
  512 continue

c     print sample of ixn pairs

      write(oarchv,132)
  132 format(/4(3x,'sites',9x,'V',2x))
      do 508 i508=1,nmres,nmres-10
         do 509 i509=numlng(i508-1,1)+1,
     *               min( numlng(i508-1,1)+5,numlng(i508,1))

            do 516 i516=1,numtab
               rmax(i516)=-100000.0D0
               rmin(i516)=1000000.0D0
               do 515 i515=1,maxs
                  rmax(i516)=
     *               max( rmax(i516),vpotnt(i515,i509,i516) )
                  rmin(i516)=
     *               min( rmin(i516),vpotnt(i515,i509,i516) )
  515           continue
  516       continue
            write(oarchv,103)(ilong(i509,1,i1),ilong(i509,2,i1),
     *                       vpotnt(1,i509,i1),i1=1,min(4,numtab))
  103       format(/4(2(1x,i3),1x,1pe9.2,2x))
            write(oarchv,106)(vpotnt(maxs,i509,i1),i1=1,4)
            write(oarchv,106)(rmin(i1),i1=1,numtab)
            write(oarchv,106)(rmax(i1),i1=1,numtab)
  106       format(4(2(1x,3x),1x,1pe9.2,2x))

  509    continue
  508 continue
 
c     ---------------------- done -----------------------

      return
      end
