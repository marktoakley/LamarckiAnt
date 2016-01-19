      subroutine  contact_P_AP(pro_cord,nmres,f_cord,lambda_P_AP)

! does a (non-additive) contact potential that weights C_alpha:C_alpha contacts
! differently if there is also a contact between the C_alpha atoms i_diff_P_AP along chain.
! (eg if i,j and i-i_diff_P_AP,j-i_diff_P_AP or i+i_diff_P_AP,j+i_diff_P_AP
!  or i-i_diff_P_AP,j+i_diff_P_AP or i+i_diff_P_AP,j-i_diff_P_AP is also in contact)
! this is intended to help beta strands line up better even before H_bonds
! come in fully.
! there are different weights according to the seq separation of the contact pair
! and whether the adjacent contact indicates a parallel or antiparallel `double contact'

      use amhglobals,  only : maxsiz,maxcrd,r_cut_P_AP,weight_P_AP,&
        maxmr,minmr,i_diff_P_AP,i_atom_P_AP

      implicit none

       double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
      integer, intent(in):: nmres

       double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
       double precision, intent(out):: lambda_P_AP(3)



      double precision :: dist,dist_factor(maxsiz,maxsiz,3),theta(maxsiz,maxsiz),&
       theta_dot(maxsiz,maxsiz),t_max
      integer :: i,j,i_med_max,i_med_min,k

!     zero force and energy

      f_cord=0.0D0
      lambda_P_AP=0.0D0


      k=i_atom_P_AP    !the atom (alpha/beta/oxygen = 1/2/3) potential is on

      i_med_max=maxmr
      i_med_min=minmr    !make edges of classes commensurate with AMH definition
                         ! however note that two pair contacts are required to
                         ! define a parallel or anti-parallel contact. we thus make
                         ! the convention that the class of the AP contact is determined
                         ! by the *shorter* seq separated pair (for P case both pairs
                         ! are at same distance so there is no complication)

      do i=1,nmres-i_med_min   !work out theta values between all pairs i<j in long/med classes
      do j=i+i_med_min,nmres   
        dist=dsqrt (  &
              (pro_cord(j,1,k)-pro_cord(i,1,k))**2 + &
              (pro_cord(j,2,k)-pro_cord(i,2,k))**2 + &
              (pro_cord(j,3,k)-pro_cord(i,3,k))**2 )
        dist_factor(i,j,:)=(pro_cord(i,:,k)-pro_cord(j,:,k))/dist
        t_max=DTANH(7.0D0*(r_cut_P_AP-dist))
        theta(i,j)=-0.5D0*(1.0D0+t_max)   !note always negative (having a function that changes sign doesn't make sense)
        theta_dot(i,j)=3.5D0*(1.0D0-t_max**2)
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! work out the three lambda values corresponding to the weight terms -weight*lambda=energy
! and also get forces 
! note that have written f_cord=f_cord+(-weight)*... to emphasize
! that the potential is of form -weight*sum(theta*theta)
!
! lambdas must be positive
! expect (but don't require) weight positive so 
! expect energy contribution to be negative
! (this is in line with dssp_hdrgn convention)
! the three classes for the lambdas/weights are
! (1) AP med range
! (2) AP long range
! (3) P  long range

      
      do i=1,nmres-(i_med_min+2*i_diff_P_AP)   !AP med range. Note loop over i,j (with j>i) and examine *shorter* seq sperarated
      do j=i+(i_med_min+2*i_diff_P_AP),min(i+i_med_max+2*i_diff_P_AP,nmres)   
                                                   !neighbour i+i_diff,j-i_diff (if did i-i_diff,j+i_diff too would overcount)
                                                   ! so pairs are [i,j] and [i+i_diff,j-i_diff] with the restriction that
                                                   ! the shorter separation (ie j-i-2*i_diff) is in the medium range class
                                                   ! so j-i-2*i_diff>=i_med_min  and  j-i-2*i_diff<=i_med_max
        lambda_P_AP(1)=lambda_P_AP(1)+theta(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(i,:,k)=f_cord(i,:,k)-(-weight_P_AP(1))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(j,:,k)=f_cord(j,:,k)+(-weight_P_AP(1))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(i+i_diff_P_AP,:,k)=f_cord(i+i_diff_P_AP,:,k)-(-weight_P_AP(1))*dist_factor(i+i_diff_P_AP,j-i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j-i_diff_P_AP)*theta(i,j)
        f_cord(j-i_diff_P_AP,:,k)=f_cord(j-i_diff_P_AP,:,k)+(-weight_P_AP(1))*dist_factor(i+i_diff_P_AP,j-i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j-i_diff_P_AP)*theta(i,j)
      enddo
      enddo

   
      do i=1,nmres-(i_med_max+2*i_diff_P_AP+1)  !AP long range  ( [i,j] + [i+i_diff,j-i_diff], j-i-2*i_diff>i_med_max )
      do j=i+(i_med_max+2*i_diff_P_AP+1),nmres           
        lambda_P_AP(2)=lambda_P_AP(2)+theta(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(i,:,k)=f_cord(i,:,k)-(-weight_P_AP(2))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(j,:,k)=f_cord(j,:,k)+(-weight_P_AP(2))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j-i_diff_P_AP)
        f_cord(i+i_diff_P_AP,:,k)=f_cord(i+i_diff_P_AP,:,k)-(-weight_P_AP(2))*dist_factor(i+i_diff_P_AP,j-i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j-i_diff_P_AP)*theta(i,j)
        f_cord(j-i_diff_P_AP,:,k)=f_cord(j-i_diff_P_AP,:,k)+(-weight_P_AP(2))*dist_factor(i+i_diff_P_AP,j-i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j-i_diff_P_AP)*theta(i,j)
      enddo
      enddo


      do i=1,nmres-(i_med_max+1+i_diff_P_AP)  !P long range ( [i,j] + [i+i_diff,j+i_diff], j-i>i_med_max )
      do j=i+(i_med_max+1),nmres-i_diff_P_AP
        lambda_P_AP(3)=lambda_P_AP(3)+theta(i,j)*theta(i+i_diff_P_AP,j+i_diff_P_AP)
        f_cord(i,:,k)=f_cord(i,:,k)-(-weight_P_AP(3))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j+i_diff_P_AP)
        f_cord(j,:,k)=f_cord(j,:,k)+(-weight_P_AP(3))*dist_factor(i,j,:)*theta_dot(i,j)*theta(i+i_diff_P_AP,j+i_diff_P_AP)
        f_cord(i+i_diff_P_AP,:,k)=f_cord(i+i_diff_P_AP,:,k)-(-weight_P_AP(3))*dist_factor(i+i_diff_P_AP,j+i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j+i_diff_P_AP)*theta(i,j)
        f_cord(j+i_diff_P_AP,:,k)=f_cord(j+i_diff_P_AP,:,k)+(-weight_P_AP(3))*dist_factor(i+i_diff_P_AP,j+i_diff_P_AP,:) &
                                                                            *theta_dot(i+i_diff_P_AP,j+i_diff_P_AP)*theta(i,j)
      enddo
      enddo


      end
