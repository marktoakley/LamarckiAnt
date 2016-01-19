      subroutine ev_set_up

      use amhglobals,  only:SO, exvmin_gamma,nmres,ires,numlng,exvminS,exvmin,&
            ccev_dist,ilong,crdixn,iexcld_beta,iexcld_gamma,c_of_m_dist,&
            exvminS_beta,o_exvminS,o_exvmin,ooev_dist,ievgamma,ievbeta

      implicit none

      double precision hrdrad_gamma(1:20),hrdrad_beta(1:20)
      double precision, parameter, dimension(20):: c_of_m_in_angstroms=(/1.52D0,3.56D0,2.09D0,2.09D0,1.88D0,&
                                                             2.63D0,2.69D0,0.00D0,2.39D0,2.19D0,&
                                                             2.31D0,2.99D0,2.55D0,2.52D0,1.42D0,&
                                                             1.78D0,1.89D0,2.73D0,2.67D0,1.91D0/)

      integer i,j,i_tab,i_ixn,i_atom,j_atom,open_status,i_dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! centre of mass distances for `gamma' excluded volume

      do i=1,nmres   !in units of CA-CB bond length
        c_of_m_dist(i)=c_of_m_in_angstroms(ires(i))/1.52D0
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  gamma excluded volume

      exvmin_gamma=0.0D0
      
      if (iexcld_gamma) then
         open(unit=ievgamma,file='params/r_ev_gamma.dat',&
           status='old',action='read',iostat=open_status)
         if (open_status.ne.0) then
             write(SO,*) 'error opening gamma ev file'
             stop
         endif
         do i=1,20
           read(ievgamma,*) i_dummy,hrdrad_gamma(i)
         enddo
      endif


      do i=1,nmres-1
      do j=i+1,nmres
        exvmin_gamma(i,j)=hrdrad_gamma(ires(i))+hrdrad_gamma(ires(j))
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (iexcld_beta) then !set *beta* values. Only used if iexcld_beta=.true.
         open(unit=ievbeta,file='params/r_ev_beta.dat',&
           status='old',action='read',iostat=open_status)
         if (open_status.ne.0) then
             write(SO,*) 'error opening beta ev file'
             stop
         endif
         do i=1,20
           read(ievbeta,*) i_dummy,hrdrad_beta(i)
         enddo
      endif

         do i_tab=1,3             !leave beta-beta as special case
         do i_ixn=1,numlng(nmres,i_tab)
            i=ilong(i_ixn,1,i_tab)
            j=ilong(i_ixn,2,i_tab)

            i_atom=crdixn(i_tab,1)
            j_atom=crdixn(i_tab,2)

            if( (i_atom.eq.2).and.(ires(i).eq.8) ) cycle
            if( (j_atom.eq.2).and.(ires(j).eq.8) ) cycle

!              set equilibrium distance

            if (iabs(j-i).lt.5) then
              ccev_dist(i_ixn,i_tab)=exvminS 
            else
              ccev_dist(i_ixn,i_tab)=exvmin
            endif

          enddo
          enddo

          i_tab=4
          do i_ixn=1,numlng(nmres,i_tab)
            i=ilong(i_ixn,1,i_tab)
            j=ilong(i_ixn,2,i_tab)

            i_atom=crdixn(i_tab,1)
            j_atom=crdixn(i_tab,2)

            if( (i_atom.eq.2).and.(ires(i).eq.8) ) cycle
            if( (j_atom.eq.2).and.(ires(j).eq.8) ) cycle


!              set equilibrium distance

            if (iabs(j-i).lt.5) then
              if (iabs(j-i).lt.1) then
                write(SO,*) 'problem setting up ev',iabs(j-i),i_ixn
                stop
              endif
              ccev_dist(i_ixn,i_tab)=exvminS_beta(iabs(j-i)) 
            elseif (iabs(j-i).lt.13) then
              ccev_dist(i_ixn,i_tab)=exvmin
            elseif(.not.iexcld_beta) then
              ccev_dist(i_ixn,i_tab)=exvmin
            else 
              ccev_dist(i_ixn,i_tab)=hrdrad_beta(ires(i))+hrdrad_beta(ires(j))
            endif
          enddo


! do O-O excluded volume

          i_tab=2
          do i_ixn=1,numlng(nmres,i_tab)

            i=ilong(i_ixn,1,i_tab)
            j=ilong(i_ixn,2,i_tab)

            if (iabs(j-i).lt.5) then
              ooev_dist(i_ixn)=o_exvminS(iabs(j-i))
            else
              ooev_dist(i_ixn)=o_exvmin
            endif
          enddo

          write(SO,*) '*****************************'
          write(SO,*) 'o_exvmin=',o_exvmin
          write(SO,*) 'o_exvminS=',o_exvminS
          write(SO,*) '*****************************'
          write(SO,*) 'exvmin=',exvmin
          write(SO,*) 'exvminS=',exvminS
          write(SO,*) 'exvminS_beta=',exvminS_beta
          write(SO,*) '*****************************'

          return
          end
      
