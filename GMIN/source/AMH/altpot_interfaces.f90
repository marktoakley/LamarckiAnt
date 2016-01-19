      module altpot_interfaces
        
        
!!!!!!!!!!!!!!! PUBLIC !!!!!!!!!!!!!!!!!!!!!!
        interface
           subroutine read_input_alt()
             use globals_alt, only : altpotflag,kappa_alt,treshold_alt,kappa_well,onebodyflag &
                  , onebody_type,kappa_OB,Alimits_OB, debug_numer_forces
             implicit none
           end subroutine read_input_alt
        end interface

        interface
           subroutine altpot_master(pro_cord,f_cord,E,tempav,nmres)
             use amhglobals,  only : maxsiz,maxcrd,r_min,r_max, ires
             use globals_alt, only : debug_numer_forces, debugflag,kappa_alt,treshold_alt &
                  ,altpotflag,max_well_alt,kappa_well,output_stepsize_alt &
                  , output_step_alt ,do_send_output , onebodyflag, onebody_type 
             implicit none
             integer, intent(in):: nmres
              double precision, intent(in) ::  pro_cord(maxsiz,3,maxcrd)
             logical, intent(in):: tempav
!              double precision, intent(out) :: f_cord(maxsiz,3,maxcrd),E(:,:)
              double precision f_cord(maxsiz,3,maxcrd),E
           end subroutine altpot_master
        end interface

!!!!!!!!!!!!!!! PRIVATE !!!!!!!!!!!!!!!!!!!!!!

        interface
           subroutine calc_theta_alt(theta, theta_dot, xyz_dist, rmin, rmax,kappa,nmres,  i_well)
             use amhglobals,  only : maxsiz,ires
             use globals_alt, only : max_well_alt
             implicit none
             integer, intent(in) ::  i_well,nmres
              double precision, intent(in) ::  kappa,rmin,rmax
              double precision, intent(in) ::  xyz_dist(maxsiz,maxsiz)
              double precision, intent(out) ::  theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
           end subroutine calc_theta_alt
        end interface



        interface
           subroutine calc_sum_theta_dot_alt(sum_theta_dot,theta_dot,xyz_unit_vect,nmres)
             use globals_alt, only : max_well_alt
             use amhglobals,  only : maxsiz
             implicit none
             integer, intent(in) ::  nmres
              double precision, intent(in) ::  theta_dot(maxsiz,maxsiz,max_well_alt)
              double precision, intent(in) :: xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(out) :: sum_theta_dot(maxsiz,3)
           end subroutine calc_sum_theta_dot_alt
        end interface


        interface
           subroutine calc_A_alt(A,theta,nmres)
             use globals_alt, only : max_well_alt
             use amhglobals,  only : maxsiz
             implicit none
             integer, intent(in) ::  nmres
              double precision, intent(in) ::   theta(maxsiz,maxsiz,max_well_alt)
              double precision, intent(out) ::  A(maxsiz)
           end subroutine calc_A_alt
        end interface

        interface
           double precision function calc_ddA_onebody(A_i,itype)
             use globals_alt, only : onebody_gamma,kappa_OB,Alimits_OB,onebody_type
             use amhglobals,  only : maxsiz
             implicit none
!             integer, intent(in) :: itype
!             double precision, intent(in) :: A_i(maxsiz)
!             integer itype
             double precision A_i
             integer itype
           end function calc_ddA_onebody
        end interface

        interface
           double precision function calc_energy_alt(theta,sigma,E_alt,nmres,i_well)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : max_well_alt,altgamma
             implicit none
             integer nmres,i_well
             double precision theta(maxsiz,maxsiz,max_well_alt),E_alt(2,max_well_alt)
             double precision sigma(maxsiz,maxsiz),tij
           end function calc_energy_alt
        end interface


        interface
           subroutine calc_force_alt(f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                , i_well, A, kappa, treshold)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time
             implicit none
             integer, intent(in) :: i_well, nmres
              double precision, intent(in) :: kappa, treshold
              double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
              double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz), xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
           end subroutine calc_force_alt
        end interface


        interface
           subroutine calc_force_alt_1(f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                , i_well, A, kappa, treshold, heaviside, heaviside_dot)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time
             implicit none
             integer, intent(in) :: i_well, nmres
              double precision, intent(in) :: kappa, treshold
              double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
              double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz), xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
           end subroutine calc_force_alt_1
         end interface


        interface
           subroutine calc_force_alt_2(f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                , i_well, A, kappa, treshold, heaviside, heaviside_dot)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time
             implicit none
             integer, intent(in) :: i_well, nmres
              double precision, intent(in) :: kappa, treshold
              double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
              double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz), xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
           end subroutine calc_force_alt_2
        end interface



        interface
           subroutine calc_force_alt_3(f_cord,pro_cord,theta,theta_dot,sigma,xyz_unit_vect,nmres &
                , i_well, A, kappa, treshold, heaviside, heaviside_dot)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : altgamma,max_well_alt,debugflag,accumulated_time
             implicit none
             integer, intent(in) :: i_well, nmres
              double precision, intent(in) :: kappa, treshold
              double precision, intent(in) :: theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
              double precision, intent(in) :: sigma(maxsiz,maxsiz), A(maxsiz), xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(in) :: heaviside(maxsiz), heaviside_dot(maxsiz)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
           end subroutine calc_force_alt_3
        end interface



        interface
           subroutine calc_heaviside_alt(heaviside,heaviside_dot,A,treshold,kappa)
             use amhglobals,  only : maxsiz
             implicit none
              double precision, intent(out) :: heaviside(maxsiz), heaviside_dot(maxsiz)
              double precision, intent(in) :: kappa,treshold
              double precision, intent(in) ::  A(maxsiz)
           end subroutine calc_heaviside_alt
        end interface



        interface
           subroutine calc_numerical_force_alt(pro_cord,f_alt_numer_pair_pot,f_alt_numer_ob_pot,nmres)
             use amhglobals,  only : maxsiz,maxcrd,r_min,r_max, ires
             use globals_alt, only : kappa_alt,treshold_alt &
                  ,altpotflag,max_well_alt,kappa_well, onebodyflag,accumulated_time
             implicit none
             integer, intent(in) :: nmres
             double precision pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out) :: f_alt_numer_pair_pot(maxsiz,3,maxcrd), f_alt_numer_ob_pot(maxsiz,3,maxcrd)
           end subroutine calc_numerical_force_alt
        end interface



        interface
           subroutine calc_OB_density(A,nmres)
             use amhglobals,  only : ires,maxsiz
             use globals_alt, only :max_letters,OB_dns_count,OB_density,Alimits_OB
             implicit none
!            integer, intent(in) ::   nmres
!               double precision, intent(in) ::  A(maxsiz)
              integer  nmres
              double precision  A(maxsiz)
           end subroutine calc_OB_density
        end interface



        interface
           subroutine calc_onebody_force(f_cord,theta_dot,xyz_unit_vect,nmres,A)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : max_well_alt,OB_density,OB_dns_count,max_letters,accumulated_time
             implicit none
             integer, intent(in) :: nmres

              double precision  A(maxsiz),theta_dot(maxsiz,maxsiz,max_well_alt)
!              double precision, intent(in) ::  A(maxsiz),theta_dot(maxsiz,maxsiz,max_well_alt)
             double precision xyz_unit_vect(maxsiz,maxsiz,3)
              double precision, intent(out) :: f_cord(maxsiz,3,maxcrd)
           end subroutine calc_onebody_force
        end interface



        interface
           double precision function calc_onebody_pot(A,nmres,E_OB)
             use amhglobals,  only : maxsiz,ires
             use globals_alt, only : onebody_gamma,kappa_OB,Alimits_OB,onebody_type
             implicit none
!             integer, intent(in) ::  nmres
!              double precision, intent(in) ::  A(maxsiz)
!              double precision, intent(out) ::  E_OB(3)
              integer  nmres
              double precision A(maxsiz)
              double precision E_OB(3)
           end function calc_onebody_pot
        end interface



        interface
           subroutine calc_sigma_alt(sigma,A,treshold,kappa,nmres)
             use amhglobals,  only : maxsiz
             use globals_alt, only : max_well_alt
             implicit none
             integer, intent(in) ::   nmres
              double precision, intent(in) ::   kappa,treshold
              double precision, intent(out) ::   sigma(maxsiz,maxsiz),A(maxsiz)
           end subroutine calc_sigma_alt
        end interface




        interface
           subroutine calc_xyz(xyz_dist,xyz_unit_vect,pro_cord,nmres)
             use amhglobals,  only : maxsiz,maxcrd,ires
             use globals_alt, only : OB_dns_count,outfile1_alt
             implicit none
             integer, intent(in) :: nmres
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out) :: xyz_dist(maxsiz,maxsiz),xyz_unit_vect(maxsiz,maxsiz,3)
           end subroutine calc_xyz
        end interface



        interface
           subroutine default_alt()
             use globals_alt, only : altpotflag,kappa_alt,treshold_alt,kappa_well,kappa_OB, & 
                  debug_numer_forces,outfile1_alt, output_step_alt, output_stepsize_alt, & 
                  do_send_output,onebodyflag, onebody_type , timings_file_alt, accumulated_time
             implicit none
           end subroutine default_alt
        end interface



        interface
           subroutine finalize_alt()
             use globals_alt, only : outfile1_alt
             implicit none
           end subroutine finalize_alt
        end interface



        interface
           subroutine init_onebody()
             use globals_alt, only : onebody_gamma,hp_scale,kappa_OB,Alimits_OB,onebodyflag &
                  ,max_letters
           end subroutine init_onebody
        end interface



        interface
           subroutine longscale()
             use amhglobals,  only:SO,  ab_c_of_n_old,ab_c_of_n_new,max_well,alpha_c_of_n &
                  ,num_well,nmres,n_letters_con
             use globals_alt, only : altgamma
             implicit none
           end subroutine longscale
        end interface



        interface
           subroutine read_altgamma()
             use amhglobals,  only : n_letters_con
             use globals_alt, only : altgamma,altpotflag,max_well_alt,onebody_gamma,onebodyflag
             implicit none
           end subroutine read_altgamma
        end interface



        interface
           subroutine send_output_alt(A,E_alt,nmres,E_OB)
             use amhglobals,  only:SO, maxsiz
             use globals_alt, only : output_step_alt, output_stepsize_alt,outfile1_alt,do_send_output &
                  ,count_alt,T_alt,max_well_alt,treshold_alt,OB_density,OB_dns_count &
                  , max_letters,aminoacids,accumulated_time
             implicit none
             integer, intent(in) ::  nmres
              double precision, intent(in) :: A(maxsiz),E_alt(2,max_well_alt),E_OB(3)
           end subroutine send_output_alt
        end interface



      end module altpot_interfaces

