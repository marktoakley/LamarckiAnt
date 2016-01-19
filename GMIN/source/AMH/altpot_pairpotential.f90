      double precision function calc_energy_alt(theta,sigma,E_alt,nmres,i_well)

        use amhglobals,  only : maxsiz,CUTOFF_CONT_LOW,tgsequences_amw,numseq_amw
        use globals_alt, only : max_well_alt,altgamma

        implicit none
        integer i,j,i_well,nmres,itype,jtype,iseq
        double precision theta(maxsiz,maxsiz,max_well_alt),E_alt(2,max_well_alt)
        double precision sigma(maxsiz,maxsiz),tij

        E_alt(1,i_well)=0.0D0
        E_alt(2,i_well)=0.0D0

        do iseq=1,numseq_amw
          do i=1,nmres-CUTOFF_CONT_LOW 
!               itype=ires(i)
                itype=tgsequences_amw(i,iseq)
             do j=i+CUTOFF_CONT_LOW,nmres
!               jtype=ires(j)
                jtype=tgsequences_amw(j,iseq)

                tij=theta(i,j,i_well)
                E_alt(1,i_well)= E_alt(1,i_well) - &
                  tij*(1.0-sigma(i,j)) * altgamma(itype,jtype,1,i_well)
                E_alt(2,i_well)= E_alt(2,i_well) - & 
                  tij*sigma(i,j) * altgamma(itype,jtype,2,i_well)
              if(sigma(i,j).gt.1.01D0)then
                 write(*,*)'sigma ',i,' ',j,' bigger than one',sigma(i,j)
                 stop
              elseif(sigma(i,j).lt. -0.01D0)then
                 write(*,*)'sigma ',i,' ',j,' less than zero',sigma(i,j)
                 stop
              endif
           enddo
        enddo
       enddo

         E_alt(1,i_well)=E_alt(1,i_well)/real(numseq_amw)
         E_alt(2,i_well)=E_alt(2,i_well)/real(numseq_amw)

        calc_energy_alt=E_alt(1,i_well)+E_alt(2,i_well)

        return

      end function calc_energy_alt

      !------------------------------



      subroutine calc_force_alt(f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                , i_well, A, kappa, treshold)
        use amhglobals,  only : maxsiz,maxcrd, tgsequences_amw,numseq_amw
        use globals_alt, only : max_well_alt,debug_numer_forces, do_send_output,accumulated_time
        use altpot_interfaces, only: calc_heaviside_alt, calc_force_alt_1, calc_force_alt_2, calc_force_alt_3

        implicit none
        integer, intent(in) :: i_well, nmres
         double precision, intent(in) :: kappa, treshold
         double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz), xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
         double precision, intent(out):: f_cord(maxsiz,3,maxcrd)

        !Local
        integer i,j,iatom,iseq,itype
        double precision forces1(maxsiz,3,maxcrd),forces2(maxsiz,3,maxcrd),forces3(maxsiz,3,maxcrd)
        double precision heaviside(maxsiz), heaviside_dot(maxsiz)

        double precision start_time, end_time

!        external cpu_time   

        call CPU_TIME(start_time)

        forces1=0.D0
        forces2=0.D0
        forces3=0.D0

        call calc_heaviside_alt(heaviside,heaviside_dot,A,treshold,kappa)

        call calc_force_alt_1 (forces1,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                              , i_well, A, kappa, treshold, heaviside, heaviside_dot)
        call calc_force_alt_2 (forces2,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                              , i_well, A, kappa, treshold, heaviside, heaviside_dot)
        call calc_force_alt_3 (forces3,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                              , i_well, A, kappa, treshold, heaviside, heaviside_dot)

        
        f_cord=f_cord+forces1/real(numseq_amw)+forces2/real(numseq_amw)+forces3/real(numseq_amw)


        if(do_send_output .and. debug_numer_forces) then
           open(unit=1561,file='alt_anal_force_detailed_1_pp',status='unknown')  
           open(unit=1562,file='alt_anal_force_detailed_2_pp',status='unknown')  
           open(unit=1563,file='alt_anal_force_detailed_3_pp',status='unknown')  

          do iseq=1,numseq_amw
           do i=1,nmres
              itype=tgsequences_amw(i,iseq)
              iatom=2
              if (itype .eq. 8) iatom=1
              write(1561,1100) (forces1(i,j,iatom),j=1,3), A(i)
              write(1562,1100) (forces2(i,j,iatom),j=1,3), A(i)
              write(1563,1100) (forces3(i,j,iatom),j=1,3), A(i)
           enddo
          enddo
           close(1561)
           close(1562)
           close(1563)
1100       format(400(1x,e12.6))
        endif

        call CPU_TIME(end_time)
        accumulated_time(2)=accumulated_time(2)+(end_time-start_time)

        return
      end subroutine calc_force_alt

      !-----------------------------------------------

      subroutine calc_force_alt_1 (f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
           , i_well, A, kappa, treshold, heaviside, heaviside_dot)
        use amhglobals,  only : maxsiz,maxcrd,ires, & 
                         CUTOFF_CONT_LOW,tgsequences_amw,numseq_amw
        use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time

        implicit none
        integer, intent(in) :: i_well, nmres
         double precision, intent(in) :: kappa, treshold
         double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: theta_dot(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz)
         double precision, intent(in) :: xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
         double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
         double precision, intent(out):: f_cord(maxsiz,3,maxcrd)


        ! Local
        integer k,i,j,katom,itype,jtype,iseq
        double precision gamma1,gamma2,deltagamma_theta,start_time, end_time

         double precision, parameter :: epsilon=0.0001D0

!        external cpu_time


        call CPU_TIME(start_time)

       do iseq=1,numseq_amw
         do i=1,nmres-CUTOFF_CONT_LOW,1 
          itype=tgsequences_amw(i,iseq)
!           itype=ires(i)
           do j=i+CUTOFF_CONT_LOW,nmres,1
! This optimizes the subroutine 15-fold
              if(abs(theta(i,j,i_well))<epsilon) cycle 
              do k=1,nmres,1 
                 if(i.eq.k .or. j.eq.k) cycle
                 katom=2
                 if (ires(k) .eq. 8) katom=1
                 jtype=tgsequences_amw(j,iseq)
!                 jtype=ires(j)
                 gamma1=altgamma(itype,jtype,1,i_well)
                 gamma2=altgamma(itype,jtype,2,i_well)
                 deltagamma_theta=(gamma2-gamma1)*theta(i,j,i_well)
                 
                 if(abs(i-k)>1 .and. abs(theta_dot(i,k,1))>epsilon) then
                    f_cord(k,:,katom)=f_cord(k,:,katom) + &
                       deltagamma_theta*heaviside_dot(i)*heaviside(j)* &
                        theta_dot(i,k,1)*xyz_unit_vect(k,i,:)
                 endif
                 if(abs(j-k)>1 .and. abs(theta_dot(j,k,1))>epsilon) then
                    f_cord(k,:,katom)=f_cord(k,:,katom) + & 
                       deltagamma_theta*heaviside(i)*heaviside_dot(j)* & 
                        theta_dot(j,k,1)*xyz_unit_vect(k,j,:) 
                 endif
              enddo
           enddo
        enddo
       enddo
        
        call CPU_TIME(end_time)
        accumulated_time(11)= accumulated_time(11) + end_time-start_time

        return
      end subroutine calc_force_alt_1

      !-----------------------------------------------


      subroutine calc_force_alt_2 (f_cord,pro_cord,theta,theta_dot,sigma, & 
                xyz_unit_vect,nmres &
           , i_well, A, kappa, treshold, heaviside, heaviside_dot)
        use amhglobals,  only : maxsiz,maxcrd,ires, CUTOFF_CONT_LOW,tgsequences_amw,numseq_amw
        use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time

        implicit none
        integer, intent(in) :: i_well, nmres
         double precision, intent(in) :: kappa, treshold
         double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: theta_dot(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz)
         double precision, intent(in) ::  xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
         double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
         double precision, intent(out):: f_cord(maxsiz,3,maxcrd)


        ! Local
        integer k,i,katom,iatom,itype,ktype,iseq
        double precision gamma1,gamma2,sum_gamma_sigma_theta_dot
        double precision force(3)

        double precision start_time, end_time

!        external cpu_time

        call CPU_TIME(start_time)

       do iseq=1,numseq_amw
 
        do k=1,nmres-CUTOFF_CONT_LOW,1
           katom=2

!           if (ires(k) .eq. 8) katom=1
!           ktype=ires(k)
!   hack for consistency
            if (tgsequences_amw(k,1) .eq. 8) katom=1

            ktype=tgsequences_amw(k,iseq)

           do i=k+CUTOFF_CONT_LOW,nmres,1
              iatom=2

              if (ires(i) .eq. 8) iatom=1
!              itype=ires(i)
!   hack for consistency
            if (tgsequences_amw(i,1) .eq. 8) iatom=1

              itype=tgsequences_amw(i,iseq)
              gamma1=altgamma(ktype,itype,1,i_well)
              gamma2=altgamma(ktype,itype,2,i_well)
              sum_gamma_sigma_theta_dot=(gamma1*(1.0D0-sigma(k,i)) + & 
                               gamma2*sigma(k,i))*theta_dot(k,i,i_well)
              force(:)=xyz_unit_vect(k,i,:)*sum_gamma_sigma_theta_dot
              f_cord(k,:,katom)=f_cord(k,:,katom)+force(:)
              f_cord(i,:,iatom)=f_cord(i,:,iatom)-force(:)
           enddo
        enddo
       enddo

        call CPU_TIME(end_time)
        accumulated_time(12)= accumulated_time(12)+(end_time-start_time)
        return

      end subroutine calc_force_alt_2

      !-----------------------------------------------



      subroutine calc_force_alt_3 (f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
           , i_well, A, kappa, treshold, heaviside, heaviside_dot)
        use amhglobals,  only : maxsiz,maxcrd, CUTOFF_CONT_LOW,tgsequences_amw,numseq_amw
        use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time
        use altpot_interfaces, only: calc_sum_theta_dot_alt 

        implicit none
        integer, intent(in) :: i_well, nmres
         double precision, intent(in) :: kappa, treshold
         double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: theta_dot(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz)
         double precision, intent(in) :: xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
         double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
         double precision, intent(out):: f_cord(maxsiz,3,maxcrd)


        ! Local
        integer k,i,katom,iatom,itype,ktype,iseq
        double precision gamma1,gamma2,deltagamma_theta
        double precision force_k(3), force_i(3)
        double precision sum_theta_dot(maxsiz,3)


        double precision start_time, end_time

!        external cpu_time

        call CPU_TIME(start_time)

        call calc_sum_theta_dot_alt(sum_theta_dot,theta_dot,xyz_unit_vect,nmres)

      do iseq=1,numseq_amw

        do k=1,nmres-CUTOFF_CONT_LOW,1
           katom=2

           if (tgsequences_amw(k,1) .eq. 8) katom=1
           ktype=tgsequences_amw(k,iseq)
                                                                                                     
!           if (ires(k) .eq. 8) katom=1
!           ktype=ires(k)
           do i=k+CUTOFF_CONT_LOW,nmres,1
              iatom=2

           if (tgsequences_amw(i,1) .eq. 8) iatom=1
           itype=tgsequences_amw(i,iseq)

!              if (ires(i) .eq. 8) iatom=1
!              itype=ires(i)

              gamma1=altgamma(ktype,itype,1,i_well)
              gamma2=altgamma(ktype,itype,2,i_well)
              deltagamma_theta=(gamma2-gamma1)*theta(k,i,i_well) 
              force_k=0.D0
              force_i=0.D0
              force_k(:)=force_k(:)+xyz_unit_vect(k,i,:)*deltagamma_theta*heaviside_dot(i)*heaviside(k)*theta_dot(k,i,1) 
              force_i(:)=force_i(:)+xyz_unit_vect(i,k,:)*deltagamma_theta*heaviside_dot(k)*heaviside(i)*theta_dot(i,k,1) 
              force_k(:)=force_k(:)+deltagamma_theta*heaviside(i)*heaviside_dot(k)*sum_theta_dot(k,:)
              force_i(:)=force_i(:)+deltagamma_theta*heaviside(k)*heaviside_dot(i)*sum_theta_dot(i,:)
              f_cord(k,:,katom)=f_cord(k,:,katom)+force_k(:)
              f_cord(i,:,iatom)=f_cord(i,:,iatom)+force_i(:)
           enddo
        enddo
      enddo

        call CPU_TIME(end_time)
        accumulated_time(13)= accumulated_time(13) + (end_time-start_time) 

        return
      end subroutine calc_force_alt_3

      !-----------------------------------------------

