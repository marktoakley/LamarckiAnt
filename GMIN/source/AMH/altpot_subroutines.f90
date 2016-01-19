!-----------------------------------------
!        A L T P O T    R E A D    M E ! 
!
!   Requires small alterations in "qchrgmk.f" and "nnet.f" and "initil.f" from standard AMH
!     
!
!    Most PARAMETERS are stored in MODULE "globals_alt"
!
!    ALTPOT_MASTER is the routine that generally should be called by a main program
!    READ_INPUT_ALT initializes the parameters and should be called
!    only  ONCE but PRIOR to any use of ALTPOT_MASTER.
!    
!    the others are "slaves" (compare with private and public in C++)
!
!    PUBLIC:
!     subroutine ALTPOT_MASTER
!     subroutine READ_INPUT_ALT
!    PRIVATE:
!     subroutine READ_ALTGAMMA    (slave)
!     subroutine default_alt()    (slave)
!     subroutine CALC_FORCE_ALT   (slave)
!     function   CALC_ENERGY_ALT  (slave)
!     subroutine CALC_XYZ         (slave)
!     subroutine CALC_THETA_ALT   (slave)
!     subroutine CALC_A_ALT       (slave)
!     subroutine CALC_SIGMA_ALT   (slave)
!
! DESCRIPTION: ALTPOT refers to potential functions used for describing the use of
!  interaction-"gamma"-parameters whose weights alters in response to different environments.
!  This is supposed to mimick a many-body potential and in particular to lett gammas vary
! depending on if residue pairs are located in the core or at the surface.
! 
! Variables: A(i) = density around site i ; Theta(i,j,well))=number of guys ij in well 
!           sigma=measure of how buried -> weighting function
!
! PARAMETERS:
!  debugflag-  set statically during compilation in module.
!  altpotflag- is set in file 'input'. Format=i1,i2,keyword "altpotflag" i1 and i1 refers to 
!              well 1 and 2 and non-zero values -> TRUE for that particular well. e.g.
!              "1,1,altpot" 
!  kappa_alt : steepness of step-function for sigma, set in file 'input' e.g. "7.0,kappa_alt"
!  kappa_well : steepness of step-functions for wells, set in file 'input' e.g. "7.0,kappa_well"
!  treshold_alt : Treshold value for contacts determines wether in or out.e.g. "3,treshold_alt"
!  kappa_OB : steepness of the wells for the Onebodypotential. set in 'input', keyword "kappa_OB"
!  onebodyflag,onebody_type : ff flag is not 0 -> onebodypot included.onebody_type makes for different
!                            types of OB_pot, see "onebody.f90"  keyword "onebodyflag"
!
! DATA INPUT: 
!       gamma.dat: 1                  2              3             4       
!       gamma.dat: gamma(direct) gamma(exposed)  i_well            keyword="altpot" somewhere in line
!       OENBODY            1                        2
!       gamma.dat: Onebodygamma(exposed)   Onebodygamma(buried)   keyword="onebody" somewhere in line
!-----------------------------------------

      subroutine altpot_master(pro_cord,f_cord,E,tempav,nmres)
  
        use amhglobals,  only : maxsiz,maxcrd,r_min,r_max, ires
        use globals_alt, only : debug_numer_forces, debugflag,kappa_alt,treshold_alt &
          ,altpotflag,max_well_alt,kappa_well, do_send_output , onebodyflag, accumulated_time

        use altpot_interfaces, only: calc_xyz, calc_theta_alt, calc_A_alt, calc_sigma_alt, &
             calc_force_alt, calc_energy_alt, calc_onebody_force, calc_onebody_pot, &
             send_output_alt, calc_numerical_force_alt

        implicit none
        integer, intent(in):: nmres
         double precision, intent(in) ::  pro_cord(maxsiz,3,maxcrd)
        logical, intent(in):: tempav
        !       double precision, intent(out) :: f_cord(maxsiz,3,maxcrd),E
        double precision f_cord(maxsiz,3,maxcrd),E

        ! Local Variables:
        integer i,j,i_well,maxloop
        double precision start_time, end_time1, end_time2
        double precision xyz_dist(maxsiz,maxsiz),xyz_unit_vect(maxsiz,maxsiz,3)
        double precision theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
        double precision A(maxsiz),sigma(maxsiz,maxsiz),E_alt(2,max_well_alt),E_OB(3) 

        integer iatom
        double precision pp_1_f_cord(maxsiz,3,maxcrd), pp_2_f_cord(maxsiz,3,maxcrd), ob_f_cord(maxsiz,3,maxcrd)
        double precision f_alt_numer_pair_pot(maxsiz,3,maxcrd), f_alt_numer_ob_pot(maxsiz,3,maxcrd)

!        external CPU_TIME

        call CPU_TIME(start_time)

        f_cord=0.0D0
        E=0.0D0
        E_alt=0.0D0
        E_OB=0.0D0
        pp_1_f_cord=0.D0;
        pp_2_f_cord=0.D0;

        !                                START SEGMENT_ALTPOTCALC
        !        IMPORTANT!!! Note the ordering is not arbitrary and should be:      IMPORTANT!!! 
        !        calc_xyz , calc_theta , calc_A , calc_sigma and calc_force
        ! 
        if((altpotflag(1).gt.0).or.(altpotflag(2).gt.0))then
           call calc_xyz (xyz_dist,xyz_unit_vect,pro_cord,nmres)

           call calc_theta_alt(theta,theta_dot,xyz_dist,r_min(1),r_max(1),kappa_well,nmres,1)
           if(altpotflag(2).gt.0)then
              call calc_theta_alt(theta,theta_dot,xyz_dist,r_min(2),r_max(2),kappa_well,nmres,2)
           endif
           call calc_A_alt(A,theta,nmres)
           call calc_sigma_alt(sigma,A,treshold_alt,kappa_alt,nmres) ! kappa OK
        endif

        call CPU_TIME(end_time1)
        accumulated_time(10)=accumulated_time(10)+(end_time1-start_time)

        maxloop=1                        ! NOTE HARDWIRE!!  max_well is  
        if(altpotflag(2).gt.0)maxloop=2  !  currently 2 , need be changed for more wells
        do i_well=1,maxloop,1
           if(altpotflag(i_well).gt.0)then
              call calc_force_alt &
                   (f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres,i_well &
                   ,A, kappa_alt, treshold_alt)
              if(tempav.or.do_send_output)E=E+calc_energy_alt(theta,sigma,E_alt,nmres,i_well)
           endif
           if(do_send_output .and. debug_numer_forces) then 
              if(i_well .eq. 1) pp_1_f_cord=f_cord
              if(i_well .eq. 2) pp_2_f_cord=f_cord-pp_1_f_cord
           endif
        enddo

        !                             ONE BODY POTENTIAL

        if(onebodyflag.eq.1)then
           call calc_onebody_force(f_cord,theta_dot,xyz_unit_vect,nmres,A)
           if(tempav.or.do_send_output)E=E+calc_onebody_pot(A,nmres,E_OB)
        endif

        !                             END SEGMENT_ALTPOTCALC


        !DEBUG FORCES:
        if(do_send_output .and. debug_numer_forces) then
           call calc_numerical_force_alt(pro_cord,f_alt_numer_pair_pot,f_alt_numer_ob_pot,nmres)
           open(unit=1282,file='alt_numer_force_pp',status='unknown')  
           open(unit=1283,file='alt_analyt_force_pp',status='unknown')  
           if(onebodyflag.eq.1)then
              ob_f_cord=f_cord-pp_1_f_cord-pp_2_f_cord
              open(unit=1284,file='alt_numer_force_ob',status='unknown')  
              open(unit=1285,file='alt_analyt_force_ob',status='unknown')  
           endif
           do i=1,nmres
              iatom=2
              if (ires(i) .eq. 8) iatom=1
              write(1282,100) (f_alt_numer_pair_pot(i,j,iatom),j=1,3)
              write(1283,100) (pp_2_f_cord(i,j,iatom),j=1,3)
              if(onebodyflag.eq.1)then
                 write(1284,100) (f_alt_numer_ob_pot(i,j,iatom),j=1,3)
                 write(1285,100) (ob_f_cord(i,j,iatom),j=1,3)
              endif
           enddo
           close(1282)
           close(1283)
           if(onebodyflag.eq.1)then
              close(1284)
              close(1285)
           endif
100        format(400(1x,e12.6))
        endif

        call CPU_TIME(end_time2)
        accumulated_time(1)=accumulated_time(1)+(end_time2-start_time)

        call send_output_alt(A,E_alt,nmres,E_OB)

      end subroutine altpot_master

      ! This is for debugging analytical forces
      subroutine calc_numerical_force_alt (pro_cord,f_alt_numer_pair_pot,f_alt_numer_ob_pot,nmres)

        use amhglobals,  only : maxsiz,maxcrd,r_min,r_max, ires
        use globals_alt, only : kappa_alt,treshold_alt &
             ,altpotflag,max_well_alt,kappa_well, onebodyflag,accumulated_time

        use altpot_interfaces, only: calc_xyz, calc_theta_alt, calc_A_alt, calc_sigma_alt, &
                           calc_force_alt, calc_energy_alt, calc_onebody_force, calc_onebody_pot

        implicit none
        integer, intent(in) :: nmres
        double precision pro_cord(maxsiz,3,maxcrd)
         double precision, intent(out) :: f_alt_numer_pair_pot(maxsiz,3,maxcrd)
         double precision, intent(out) ::  f_alt_numer_ob_pot(maxsiz,3,maxcrd)


        ! Local Variables (JU):
        double precision xyz_dist(maxsiz,maxsiz),xyz_unit_vect(maxsiz,maxsiz,3)
        double precision theta(maxsiz,maxsiz,max_well_alt), theta_dot(maxsiz,maxsiz,max_well_alt)
        double precision A(maxsiz),sigma(maxsiz,maxsiz),E_alt(2,max_well_alt),E_OB(3)

        ! Local Variables (GAP):

        integer, parameter :: r16 = SELECTED_REAL_KIND(16, 50) ! at least 16 decimal digits of precision and the exponent range 200
        !      double precision (kind = r16) :: E_pair_pot(2), E_ob_pot(2)
        double precision :: E_pair_pot(2), E_ob_pot(2)
         double precision, parameter :: step_size=0.005D0

        integer k,katom,xyz,step

        integer, parameter :: WELL=2 !NOTE: 2nd well is hardwired for numerical forces, see also printing in altpot_master

        double precision save_pro_cord(maxsiz,3,maxcrd), saved_cord, start_time, end_time

!        external CPU_TIME

        call CPU_TIME(start_time)


        save_pro_cord=pro_cord

        f_alt_numer_pair_pot=0.D0
        f_alt_numer_ob_pot=0.D0

        do k=1,nmres
           katom=2
           if (ires(k) .eq. 8) katom=1
           do xyz=1,3
              pro_cord=save_pro_cord
              saved_cord=pro_cord(k,xyz,katom)
              do step=-1,1,2
                 pro_cord(k,xyz,katom)=saved_cord+step_size*step

                 E_alt=0.0D0
                 E_OB=0.0D0

                 call calc_xyz (xyz_dist,xyz_unit_vect,pro_cord,nmres)

                 call calc_theta_alt(theta,theta_dot,xyz_dist,r_min(1),r_max(1),kappa_well,nmres,1)
                 call calc_theta_alt(theta,theta_dot,xyz_dist,r_min(2),r_max(2),kappa_well,nmres,2)
                 call calc_A_alt(A,theta,nmres)
                 call calc_sigma_alt(sigma,A,treshold_alt,kappa_alt,nmres) 

                 if(altpotflag(WELL).gt.0)then
                    E_pair_pot((step+3)/2)=calc_energy_alt(theta,sigma,E_alt,nmres,WELL)
                 endif

                 if(onebodyflag.eq.1)then
                    E_ob_pot((step+3)/2)=calc_onebody_pot(A,nmres,E_OB)
                    !                  write(*,*) "NumDiff 2: ",k,xyz,step,E_ob_pot((step+3)/2)
                 endif

              enddo

              if(altpotflag(WELL).gt.0)then
                 f_alt_numer_pair_pot(k,xyz,katom)=-1*(E_pair_pot(2)-E_pair_pot(1))/2/step_size
                 !           write(*,664) WELL,k,xyz,E_pair_pot(1),E_pair_pot(2),f_alt_numer_pair_pot(k,xyz,katom)
                 !664        format('NumDiff PP: WELL=',i3,' k=', i3,' xyz=',i3,' E1=',E16.8,' E2=',E16.8,' F=',E16.8)
              endif
              if(onebodyflag.eq.1)then
                 f_alt_numer_ob_pot(k,xyz,katom)=-1*(E_ob_pot(2)-E_ob_pot(1))/2/step_size
                 !           write(*,665) k,xyz,E_ob_pot(1),E_ob_pot(2),f_alt_numer_ob_pot(k,xyz,katom)
                 !665        format('NumDiff OB: k=', i3,' xyz=',i3,' E1=',E16.8,' E2=',E16.8,' F=',E16.8)
              endif
           enddo
        enddo

        pro_cord=save_pro_cord

        call CPU_TIME(end_time)
        accumulated_time(4) = accumulated_time(4) + (end_time-start_time)

        return
      end subroutine calc_numerical_force_alt

