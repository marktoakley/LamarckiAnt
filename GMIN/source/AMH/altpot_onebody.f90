      subroutine init_onebody()
  
        use globals_alt, only : onebody_gamma,hp_scale,kappa_OB,Alimits_OB,onebodyflag &
             ,max_letters,OB_density,OB_dns_count
        implicit none
        integer i

        onebodyflag=0

        Alimits_OB(1,1)=0.0D0
        Alimits_OB(2,1)=3.0D0
        Alimits_OB(1,2)=3.0D0
        Alimits_OB(2,2)=6.0D0
        Alimits_OB(1,2)=6.0D0
        Alimits_OB(2,2)=9.0D0

        kappa_OB=3.0D0
        do i=1,max_letters,1
           onebody_gamma(i,1)=1.- hp_scale(i)
           onebody_gamma(i,2)=hp_scale(i)*.5D0 +.25D0
           onebody_gamma(i,3)=hp_scale(i)
        enddo

        OB_density=0.0D0
        OB_dns_count=0.0D0

        return
      end subroutine init_onebody

      !----------------------------

      subroutine calc_OB_density(A,nmres)

        use amhglobals,  only : maxsiz,tgsequences_amw,numseq_amw
        use globals_alt, only :OB_dns_count,OB_density,Alimits_OB
        implicit none
        double precision A(maxsiz),A_i
        integer well,itype,nmres,i,iseq

        OB_dns_count=OB_dns_count+1.0D0
        OB_density=0.0D0

       do iseq=1,numseq_amw
        do i=1,nmres
           itype=tgsequences_amw(i,iseq)
           A_i=A(i)
           well=1
           if(A_i .gt. Alimits_OB(1,2) .and. A_i .le. Alimits_OB(2,2) )well=2
           if(A_i .gt. Alimits_OB(1,3) .and. A_i .le. Alimits_OB(2,3) )well=3
           OB_density(itype,well)=OB_density(itype,well)+1.0
        enddo
       enddo

          OB_density=OB_density/real(numseq_amw)

        return
      end subroutine calc_OB_density

      !-------------------------------------------
      !
      ! Energy_OB = -gamma*f(A)
      ! Force_OB  = +gamma*df(A)/dr = +gamma*[df(A)/dA] * [DA/dr] = 
      !  gamma*sum(i,k) {[dtheta/ddrij] }

      subroutine calc_onebody_force(f_cord,theta_dot,xyz_unit_vect,nmres,A)

        use amhglobals,  only : maxsiz,maxcrd,tgsequences_amw,numseq_amw
        use globals_alt, only : max_well_alt,accumulated_time
! debugging   use globals_alt, only :  OB_density,OB_dns_count,  
        use altpot_interfaces, only: calc_OB_density, calc_ddA_onebody


        implicit none
       double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
!        double precision f_cord(maxsiz,3,maxcrd)
        integer, intent(in):: nmres
        double precision A(maxsiz),theta_dot(maxsiz,maxsiz,max_well_alt)
        double precision f_OB(maxsiz,3),xyz_unit_vect(maxsiz,maxsiz,3)
        double precision ddAmat(maxsiz)
! debugging   double precision OB_force(maxsiz)
        integer itype,i,j,iatom,iseq

        double precision start_time, end_time

!        external cpu_time 

        call CPU_TIME(start_time)

        f_OB=0.0D0

       do iseq=1,numseq_amw,1
        do i=1,nmres,1
           itype=tgsequences_amw(i,iseq)
           ddAmat(i)=calc_ddA_onebody(A(i),itype)
        enddo
       enddo


        do j=1,nmres,1
           do i=1,nmres,1
              if(i.ne.j) then
                 f_OB(j,:)=f_OB(j,:)-ddAmat(i)*theta_dot(j,i,1)*xyz_unit_vect(j,i,:)
              endif
              if(abs(i-j).gt.1) then
                 f_OB(j,:)=f_OB(j,:)-ddAmat(j)*theta_dot(j,i,1)*xyz_unit_vect(j,i,:)
              endif
           enddo
        enddo

        ! For debugging, should be switched off
        !      open(unit=1001,file='kolla_OB_force',status='unknown')  
        !       if(OB_dns_count.lt.3.)rewind(1001)
        !       write(1001,100)OB_dns_count,(f_OB(i,1),i=1,30)
        !      close(1001)xyz_unit_vect(i,j,:)
        !100   format(200(1x,e12.6))

        do i=1,nmres,1                  
           iatom=2
!    hack for consistency
           if (tgsequences_amw(i,1) .eq. 8) iatom=1
           f_cord(i,:,iatom)=f_cord(i,:,iatom)+f_OB(i,:)
        enddo

        call calc_OB_density(A,nmres)


        call CPU_TIME(end_time)
        accumulated_time(3)=accumulated_time(3)+end_time-start_time
        return
      end subroutine calc_onebody_force



      !----------------------------------


      double precision function calc_onebody_pot(A,nmres,E_OB)

        use amhglobals,  only : maxsiz,tgsequences_amw,numseq_amw
        use globals_alt, only : onebody_gamma,kappa_OB,Alimits_OB,onebody_type

        implicit none
        double precision A(maxsiz),obpot1,obpot2,obpot3,ENERGY,E_OB(3)
        integer itype,i,nmres,iseq

        ENERGY=0.0D0
        E_OB=0.0D0

        if(onebody_type.eq.1)then
            do iseq=1,numseq_amw
             do i=1,nmres
                 itype=tgsequences_amw(i,iseq)
                 ENERGY=ENERGY-onebody_gamma(itype,1)*A(i)
             enddo
            enddo
           E_OB(1)=ENERGY/real(numseq_amw)
        elseif(onebody_type.eq.2)then
          do iseq=1,numseq_amw
           do i=1,nmres
               itype=tgsequences_amw(i,iseq)
!              itype=ires(i)
              obpot1=(tanh(kappa_OB*(A(i)-Alimits_OB(1,1)))+ &
                   tanh(kappa_OB*(Alimits_OB(2,1)-A(i))))*onebody_gamma(itype,1)
              obpot2=(tanh(kappa_OB*(A(i)-Alimits_OB(1,2)))+ &
                   tanh(kappa_OB*(Alimits_OB(2,2)-A(i))))*onebody_gamma(itype,2)
              obpot3=(tanh(kappa_OB*(A(i)-Alimits_OB(1,3)))+ &
                   tanh(kappa_OB*(Alimits_OB(2,3)-A(i))))*onebody_gamma(itype,3)
              E_OB(1)=E_OB(1)-obpot1*0.5D0
              E_OB(2)=E_OB(2)-obpot2*0.5D0
              E_OB(3)=E_OB(3)-obpot3*0.5D0
              !       if(i.eq.30) write(*,*) "OB Details: ", A(i)
           enddo
          enddo
            ENERGY=E_OB(1)/real(numseq_amw)+ & 
                        E_OB(2)/real(numseq_amw)+ & 
                        E_OB(3)/real(numseq_amw)
           !     write(*,*) "OB 2: ", E_OB(1),E_OB(2),E_OB(3),ENERGY
        else
           ENERGY=0.0D0
        endif

        calc_onebody_pot=ENERGY
        return
      end function calc_onebody_pot


      !---------------------------------------
      ! 
      double precision function calc_ddA_onebody(A_i,itype)

        use globals_alt, only : onebody_gamma,kappa_OB,Alimits_OB,onebody_type
        implicit none
        integer itype
        double precision A_i,t11,t12,t21,t22,t13,t23,ddA1,ddA2,ddA3

        if(onebody_type.eq.1)then
           calc_ddA_onebody=onebody_gamma(itype,1)
        elseif(onebody_type.eq.2)then
           t11=tanh(kappa_OB*(A_i-Alimits_OB(1,1)))
           t21=tanh(kappa_OB*(Alimits_OB(2,1)-A_i))
           t12=tanh(kappa_OB*(A_i-Alimits_OB(1,2)))
           t22=tanh(kappa_OB*(Alimits_OB(2,2)-A_i))
           t13=tanh(kappa_OB*(A_i-Alimits_OB(1,3)))
           t23=tanh(kappa_OB*(Alimits_OB(2,3)-A_i))
           ddA1=(t21**2 - t11**2 )*kappa_OB*onebody_gamma(itype,1)*0.5D0
           ddA2=(t22**2 - t12**2 )*kappa_OB*onebody_gamma(itype,2)*0.5D0
           ddA3=(t23**2 - t13**2 )*kappa_OB*onebody_gamma(itype,3)*0.5D0
           calc_ddA_onebody=-1*(ddA1+ddA2+ddA3)
        else
           calc_ddA_onebody=0.0D0
        endif

        return
      end function calc_ddA_onebody
