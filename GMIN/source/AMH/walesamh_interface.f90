
      subroutine walesamh_interface(coord_mcp,grad_for_wales,e_for_wales)

!     ---------------------------------------------------
!     calculates energies and forces for gmin via a scltab lookup
!     to avoid dynamics routines
!     ---------------------------------------------------

     use amhglobals, only: maxtab,maxr,prcord,zrcord,avep,vpotnt, &
       forse,maxsiz,nmres,pexcld,numlng,rincinv, rincsq,ilong,crdixn,ires, &
       eqdist,hbond,oxexcldv,avep,maxpro,maxcrd,x_mcp, weight_P_AP
! aminoa, weight_P_AP

     use commons

      implicit none

      DOUBLE PRECISION grad_for_wales(3*3*nmres),e_for_wales, &
                       coord_mcp(3*3*nmres), f_cord(maxsiz,3,maxcrd), & 
                       E_temp_mcp,trgeng(maxtab,3), E_PE_no_bias,E_P_AP   

      logical tempav,scl_call

      integer jstrt,jfins, gly_count, iii

      character*3 res_type

      LOGICAL :: CALLED=.FALSE.

!     required subroutines:
      external force, harm_spring

!     --------------------- begin -----------------------

      tempav=.true.
        res_type='ALA'

!       write(6,*) 'MCP mc_count mylbgfs_count ', mc_count, mylbgfs_count 

!     IF (mc_count.EQ.0 .and. mylbgfs_count.EQ.0 ) THEN
      IF (.NOT.CALLED) THEN

        gly_count = 0

      do iii = 1,nmres
       if (ires(iii).eq.8) then
        prcord(iii,1,1,1) = x_mcp(9*(iii-1) + 1 - gly_count*3) !  CA X
        prcord(iii,2,1,1) = x_mcp(9*(iii-1) + 2 - gly_count*3) !  CA Y
        prcord(iii,3,1,1) = x_mcp(9*(iii-1) + 3 - gly_count*3) !  CA Z
!    SWAP  CA for CB
        prcord(iii,1,1,2) = x_mcp(9*(iii-1) + 1 - gly_count*3) !  CB X
        prcord(iii,2,1,2) = x_mcp(9*(iii-1) + 2 - gly_count*3) !  CB Y
        prcord(iii,3,1,2) = x_mcp(9*(iii-1) + 3 - gly_count*3) !  CB Z
        prcord(iii,1,1,3) = x_mcp(9*(iii-1) + 4 - gly_count*3) !  O X
        prcord(iii,2,1,3) = x_mcp(9*(iii-1) + 5 - gly_count*3) !  O Y
        prcord(iii,3,1,3) = x_mcp(9*(iii-1) + 6 - gly_count*3) !  O Z

        gly_count = gly_count +1

!    write(6,*)'prcord ',prcord(iii,1,1,1),prcord(iii,2,1,1),prcord(iii,3,1,1)
!    write(6,*)'prcord ',prcord(iii,1,1,2),prcord(iii,2,1,2),prcord(iii,3,1,2)
!    write(6,*)'prcord ',prcord(iii,1,1,3),prcord(iii,2,1,3),prcord(iii,3,1,3)
      else
        prcord(iii,1,1,1) = x_mcp(9*(iii-1) + 1 - gly_count*3) !  CA X
        prcord(iii,2,1,1) = x_mcp(9*(iii-1) + 2 - gly_count*3) !  CA Y
        prcord(iii,3,1,1) = x_mcp(9*(iii-1) + 3 - gly_count*3) !  CA Z
        prcord(iii,1,1,2) = x_mcp(9*(iii-1) + 4 - gly_count*3) !  CB X
        prcord(iii,2,1,2) = x_mcp(9*(iii-1) + 5 - gly_count*3) !  CB Y
        prcord(iii,3,1,2) = x_mcp(9*(iii-1) + 6 - gly_count*3) !  CB Z
        prcord(iii,1,1,3) = x_mcp(9*(iii-1) + 7 - gly_count*3) !  O X
        prcord(iii,2,1,3) = x_mcp(9*(iii-1) + 8 - gly_count*3) !  O Y
        prcord(iii,3,1,3) = x_mcp(9*(iii-1) + 9 - gly_count*3) !  O Z

!        write(6,*)'prcord ',prcord(iii,1,1,1),prcord(iii,2,1,1),prcord(iii,3,1,1)
!        write(6,*)'prcord ',prcord(iii,1,1,2),prcord(iii,2,1,2),prcord(iii,3,1,2)
!        write(6,*)'prcord ',prcord(iii,1,1,3),prcord(iii,2,1,3),prcord(iii,3,1,3)
       endif
      enddo

         CALLED=.TRUE.
       
         ELSE
      
        gly_count = 0

       do 1964 iii = 1,nmres
        if (ires(iii).eq.8) then
         prcord(iii, 1, 1, 1) = (coord_mcp(9*(iii-1)+1- gly_count*3)) !  CA X
         prcord(iii, 2, 1, 1) = (coord_mcp(9*(iii-1)+2- gly_count*3)) !  CA Y
         prcord(iii, 3, 1, 1) = (coord_mcp(9*(iii-1)+3- gly_count*3)) !  CA Z
!    SWAP  CA for CB
         prcord(iii, 1, 1, 2) = (coord_mcp(9*(iii-1)+1- gly_count*3)) !  CB X
         prcord(iii, 2, 1, 2) = (coord_mcp(9*(iii-1)+2- gly_count*3)) !  CB Y
         prcord(iii, 3, 1, 2) = (coord_mcp(9*(iii-1)+3- gly_count*3)) !  CB Z
         prcord(iii, 1, 1, 3) = (coord_mcp(9*(iii-1)+4- gly_count*3)) !  O X
         prcord(iii, 2, 1, 3) = (coord_mcp(9*(iii-1)+5- gly_count*3)) !  O Y
         prcord(iii, 3, 1, 3) = (coord_mcp(9*(iii-1)+6- gly_count*3)) !  O Z
         gly_count = gly_count +1
        else
         prcord(iii, 1, 1, 1) = (coord_mcp(9*(iii-1)+1- gly_count*3)) !  CA X
         prcord(iii, 2, 1, 1) = (coord_mcp(9*(iii-1)+2- gly_count*3)) !  CA Y
         prcord(iii, 3, 1, 1) = (coord_mcp(9*(iii-1)+3- gly_count*3)) !  CA Z
         prcord(iii, 1, 1, 2) = (coord_mcp(9*(iii-1)+4- gly_count*3)) !  CB X
         prcord(iii, 2, 1, 2) = (coord_mcp(9*(iii-1)+5- gly_count*3)) !  CB Y
         prcord(iii, 3, 1, 2) = (coord_mcp(9*(iii-1)+6- gly_count*3)) !  CB Z
         prcord(iii, 1, 1, 3) = (coord_mcp(9*(iii-1)+7- gly_count*3)) !  O X
         prcord(iii, 2, 1, 3) = (coord_mcp(9*(iii-1)+8- gly_count*3)) !  O Y
         prcord(iii, 3, 1, 3) = (coord_mcp(9*(iii-1)+9- gly_count*3)) !  O Z
        endif

1964     continue
        end if 

      scl_call=.true.

      jstrt=1
      jfins=nmres

        E_temp_mcp=0.0D0
        f_cord=0.0D0
        E_P_AP=0.0D0           
        e_for_wales=0.0D0 
        grad_for_wales(:)=0.0D0
        E_PE_no_bias=0.0D0
        zrcord(:,:,1,:)=0.D0

       call harm_spring(prcord,jstrt,jfins,maxpro,maxcrd,ires,f_cord,E_temp_mcp)

       call force(1,prcord,zrcord,avep,tempav,maxr,vpotnt,forse,trgeng, &
                  pexcld,numlng,nmres,rincinv,rincsq,ilong,crdixn,ires, & 
                  eqdist,hbond,oxexcldv,scl_call)

        zrcord(:,:,1,:)=zrcord(:,:,1,:)+f_cord

!        E_P_AP=-weight_P_AP(1)*avep(1,1,40)-weight_P_AP(2)*avep(1,1,41)-weight_P_AP(3)*avep(1,1,42)

         e_for_wales =                  & !
!                       avep(1,1,1)  +   & ! hydrn potential
                        avep(1,1,2)  +   & ! rama potential
                        avep(1,1,3)  +   & ! peptide bond
                        avep(1,1,4)  +   & ! chirality
                        avep(1,1,5)  +   & ! total amh
!                       avep(1,1,6)  +   & ! ?????  empty
!                       avep(1,1,7)  +   & ! amhsr
!                       avep(1,1,8)  +   & ! amhlr
                        avep(1,1,9)  +   & ! carbon ex vol
                        avep(1,1,10) +   & ! q_bias_a
                        avep(1,1,11) +   & ! oxy ex vol
                        avep(1,1,12) +   & ! obias_Rg
!                       avep(1,1,13) +   & ! amhmr
!                       avep(1,1,14) +   & ! ohdrgn_s   avep(i_pro,1,21:24)
!                       avep(1,1,15) +   & ! ohdrgn_m   avep(i_pro,1,25:28) 
!                       avep(1,1,16) +   & ! ohdrgn_l   avep(i_pro,1,29:32)
!                       avep(1,1,17) +   & ! non_add  
!                       avep(1,1,18) +   & ! q_bias_b
!                       avep(1,1,19) +   & ! ?????  empty
!                       avep(1,1,20) +   & ! ?????  empty
!                       avep(1,1,40) +   & ! E_P_AP_1
!                       avep(1,1,41) +   & ! E_P_AP_2
!                       avep(1,1,42) +   & ! E_P_AP_3
!                       E_P_AP       +   & ! E_P_AP
!                       avep(1,1,43) +   & ! replica energy term
                        E_temp_mcp

!     i_pro                   index over proteins
!     i_step                  index over time steps per temperature
!     i_temp                  index over temperatures
!     itcnt                   increment T (temperature) counter

        gly_count = 0
       do 764 iii = 1,nmres
       if (ires(iii).eq.8) then
         coord_mcp(9*(iii-1)+1 - gly_count*3) = (prcord(iii, 1, 1, 1)) !  CA X
         coord_mcp(9*(iii-1)+2 - gly_count*3) = (prcord(iii, 2, 1, 1)) !  CA Y
         coord_mcp(9*(iii-1)+3 - gly_count*3) = (prcord(iii, 3, 1, 1)) !  CA Z
         coord_mcp(9*(iii-1)+4 - gly_count*3) = (prcord(iii, 1, 1, 3)) !  O X
         coord_mcp(9*(iii-1)+5 - gly_count*3) = (prcord(iii, 2, 1, 3)) !  O Y
         coord_mcp(9*(iii-1)+6 - gly_count*3) = (prcord(iii, 3, 1, 3)) !  O Z
         gly_count = gly_count + 1 
       else
         coord_mcp(9*(iii-1)+1 - gly_count*3) = (prcord(iii, 1, 1, 1)) !  CA X
         coord_mcp(9*(iii-1)+2 - gly_count*3) = (prcord(iii, 2, 1, 1)) !  CA Y
         coord_mcp(9*(iii-1)+3 - gly_count*3) = (prcord(iii, 3, 1, 1)) !  CA Z
         coord_mcp(9*(iii-1)+4 - gly_count*3) = (prcord(iii, 1, 1, 2)) !  CB X
         coord_mcp(9*(iii-1)+5 - gly_count*3) = (prcord(iii, 2, 1, 2)) !  CB Y
         coord_mcp(9*(iii-1)+6 - gly_count*3) = (prcord(iii, 3, 1, 2)) !  CB Z
         coord_mcp(9*(iii-1)+7 - gly_count*3) = (prcord(iii, 1, 1, 3)) !  O X
         coord_mcp(9*(iii-1)+8 - gly_count*3) = (prcord(iii, 2, 1, 3)) !  O Y
         coord_mcp(9*(iii-1)+9 - gly_count*3) = (prcord(iii, 3, 1, 3)) !  O Z
       endif
764    continue

!  avep(1,1,1,itcnt)  = hydrogen bond potential
!  avep(1,1,2,itcnt)  = rama potential
!  avep(1,1,3,itcnt)  = peptide bond plainarity potential
!  avep(1,1,4,itcnt)  = chirality potential
!  avep(1,1,5,itcnt)  = total amh short + medium + long
!  avep(1,1,6,itcnt)  = ------
!  avep(1,1,7,itcnt)  = amh short range
!  avep(1,1,8,itcnt)  = amh long range
!  avep(1,1,9,itcnt)  = potential due to carbon (a/b) excluded volume
!  avep(1,1,10,itcnt) =  Q biasing potental
!  avep(1,1,11,itcnt) = oxy potential
!  avep(1,1,12,itcnt) =  radius of Gyration bias
!  avep(1,1,13,itcnt) = amh medium range
!  avep(1,1,14,itcnt) = hydrogen bonds short range
!  avep(1,1,15,itcnt) = hydrogen bonds medium range
!  avep(1,1,16,itcnt) = hydrogen bonds long range
!  avep(1,1,17,itcnt) = nonaddative
!  avep(1,1,18,itcnt) = Q biasing potental seg

!        write(6,*)'E shake springs       ', E_temp_mcp
!        write(6,*)'e 1  ohdrgn           ',avep(1,1,1)
!        write(6,*)'e 2  orama            ',avep(1,1,2)
!        write(6,*)'e 3  ooxy             ',avep(1,1,3)
!        write(6,*)'e 4  ochiral          ',avep(1,1,4)
!        write(6,*)'e 5  oamh             ',avep(1,1,5)
!        write(6,*)'e 6  ----             ',avep(1,1,6)
!        write(6,*)'e 7  oamhsr           ',avep(1,1,7)
!        write(6,*)'e 8  oamhlr           ',avep(1,1,8)
!        write(6,*)'e 9  occev            ',avep(1,1,9)
!        write(6,*)'e 10 Qbias_a          ',avep(1,1,10)
!        write(6,*)'e 11 ooev             ',avep(1,1,11)
!        write(6,*)'e 12 rg_bias          ',avep(1,1,12)
!        write(6,*)'e 13 amhmr            ',avep(1,1,13)
!        write(6,*)'e 14 hdrgns           ',avep(1,1,14)
!        write(6,*)'e 15 hdrgnm           ',avep(1,1,15)
!        write(6,*)'e 16 hdrgnl           ',avep(1,1,16)
!        write(6,*)'e 17 non add contact  ',avep(1,1,17)
!        write(6,*)'e 18 Qbias_b          ',avep(1,1,18)
!        write(6,*)'e E_P_AP              ',E_P_AP 

        gly_count = 0

      do 513 iii=1,nmres
         if (ires(iii).eq.8) then
          grad_for_wales(9*(iii-1)+1 - gly_count*3) = -(zrcord(iii,1,1,1)) ! CA x
          grad_for_wales(9*(iii-1)+2 - gly_count*3) = -(zrcord(iii,2,1,1)) ! CA y 
          grad_for_wales(9*(iii-1)+3 - gly_count*3) = -(zrcord(iii,3,1,1)) ! CA z 
          grad_for_wales(9*(iii-1)+4 - gly_count*3) = -(zrcord(iii,1,1,3)) ! O  x
          grad_for_wales(9*(iii-1)+5 - gly_count*3) = -(zrcord(iii,2,1,3)) ! O  y 
          grad_for_wales(9*(iii-1)+6 - gly_count*3) = -(zrcord(iii,3,1,3)) ! O  z 
          gly_count = gly_count +1
         else  
          grad_for_wales(9*(iii-1)+1 - gly_count*3) = -(zrcord(iii,1,1,1)) ! CA x
          grad_for_wales(9*(iii-1)+2 - gly_count*3) = -(zrcord(iii,2,1,1)) ! CA y 
          grad_for_wales(9*(iii-1)+3 - gly_count*3) = -(zrcord(iii,3,1,1)) ! CA z 
          grad_for_wales(9*(iii-1)+4 - gly_count*3) = -(zrcord(iii,1,1,2)) ! CB x 
          grad_for_wales(9*(iii-1)+5 - gly_count*3) = -(zrcord(iii,2,1,2)) ! CB y 
          grad_for_wales(9*(iii-1)+6 - gly_count*3) = -(zrcord(iii,3,1,2)) ! CB z 
          grad_for_wales(9*(iii-1)+7 - gly_count*3) = -(zrcord(iii,1,1,3)) ! O  x
          grad_for_wales(9*(iii-1)+8 - gly_count*3) = -(zrcord(iii,2,1,3)) ! O  y 
          grad_for_wales(9*(iii-1)+9 - gly_count*3) = -(zrcord(iii,3,1,3)) ! O  z 
         endif
513    continue

      return
      end
