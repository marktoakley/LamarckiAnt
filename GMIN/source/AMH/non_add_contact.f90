      subroutine  non_add_contact(pro_cord,nmres,tempav,f_cord,E)

! does non-additive contact potential

      use amhglobals,  only : maxsiz,maxcrd,ires,r_min,r_max,sort_non_add,&
            gamma_non_add,class_of_res_2

      implicit none

       double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
      integer, intent(in):: nmres
      logical, intent(in):: tempav

       double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
       double precision, intent(out):: E

      double precision :: dist,dist_factor(maxsiz,maxsiz,3),theta(maxsiz,maxsiz,3),A(maxsiz,2,3),&
       weight(maxsiz,maxsiz,3),weight2(maxsiz,2,3),f_non_add(maxsiz,maxsiz),&
       theta_dot(maxsiz,maxsiz,3),t_min,t_max
      integer :: i,j,iatom,jatom,i_well,i_well2,T_i,T_j,igam1,igam2,i_type

!     zero force and energy

      f_cord=0.0D0
      E=0.0D0

      do i=1,nmres-1   !work out theta values between all pairs
      do j=i+1,nmres
        iatom=2
        jatom=2
        if (ires(i) .eq. 8) iatom=1
        if (ires(j) .eq. 8) jatom=1
        dist=sqrt (  &
              (pro_cord(j,1,jatom)-pro_cord(i,1,iatom))**2 + &
              (pro_cord(j,2,jatom)-pro_cord(i,2,iatom))**2 + &
              (pro_cord(j,3,jatom)-pro_cord(i,3,iatom))**2 )
        dist_factor(i,j,:)=(pro_cord(i,:,iatom)-pro_cord(j,:,jatom))/dist
        do i_well=1,2
          t_min=tanh(7.0D0*(dist-r_min(i_well)))
          t_max=tanh(7.0D0*(r_max(i_well)-dist))
          theta(i,j,i_well) = 0.25D0*( 1.0D0+t_min )*( 1.0D0+t_max ) 
          theta_dot(i,j,i_well)=7.0D0*theta(i,j,i_well)*(t_max-t_min) 
        enddo
      enddo
      enddo


      A=0.0D0
      do i=1,nmres-1     !calculate local densities, A
      do j=i+1,nmres
        do i_well=1,2
          A(i,class_of_res_2(j),i_well)=A(i,class_of_res_2(j),i_well)+theta(i,j,i_well)
          A(j,class_of_res_2(i),i_well)=A(j,class_of_res_2(i),i_well)+theta(i,j,i_well)
        enddo
      enddo 
      enddo


      weight=0.0D0
      do i=1,nmres-1     !calculate weights given to bare interactions
      T_i=class_of_res_2(i)
      do j=i+1,nmres
      T_j=class_of_res_2(j)
        igam1=sort_non_add(T_i,T_j,1,1,1)
        igam2=sort_non_add(T_j,T_i,1,1,1)
        do i_well=1,2
          do i_type=1,2
          do i_well2=1,2
            weight(i,j,i_well)=weight(i,j,i_well)+ &
                 gamma_non_add(igam1)*A(i,i_type,i_well2)+ &
                 gamma_non_add(igam2)*A(j,i_type,i_well2)
            igam1=igam1+1
            igam2=igam2+1
          enddo
          enddo
        enddo
      enddo
      enddo



      weight2=0.0D0
      do i=1,nmres-13
      T_i=class_of_res_2(i)
      do j=i+13,nmres
      T_j=class_of_res_2(j)
      igam1=sort_non_add(T_i,T_j,1,1,1)
      igam2=sort_non_add(T_j,T_i,1,1,1)
      do i_well=1,2
      do i_type=1,2
      do i_well2=1,2
        weight2(i,i_type,i_well2)=weight2(i,i_type,i_well2)+  &
                    theta(i,j,i_well)*gamma_non_add(igam1)
        weight2(j,i_type,i_well2)=weight2(j,i_type,i_well2)+  &
                    theta(i,j,i_well)*gamma_non_add(igam2)
        igam1=igam1+1
        igam2=igam2+1
      enddo
      enddo
      enddo
      enddo
      enddo
  

      if (tempav) then        !calculate E
        E=0.0D0
        do i=1,nmres-13
        do j=i+13,nmres
        do i_well=1,2
          E=E-theta(i,j,i_well)*weight(i,j,i_well)
        enddo
        enddo
        enddo
      endif

      f_non_add=0.0D0    !calculate 1st contribution to F_ij
      do i=1,nmres-13
      do j=i+13,nmres
      do i_well=1,2
        f_non_add(i,j)=f_non_add(i,j)+theta_dot(i,j,i_well)*weight(i,j,i_well)
      enddo
      enddo
      enddo

      do i=1,nmres-1   !add 2nd contribution to F_ij
      T_i=class_of_res_2(i)
      do j=i+1,nmres
      T_j=class_of_res_2(j)
      do i_well2=1,2
        f_non_add(i,j)=f_non_add(i,j)+theta_dot(i,j,i_well2)* &
                ( weight2(i,T_j,i_well2) + weight2(j,T_i,i_well2) )
      enddo
      enddo
      enddo


      do i=1,nmres-1                  !calculate f_cord
      iatom=2
      if (ires(i) .eq. 8) iatom=1
      do j=i+1,nmres
      jatom=2
      if (ires(j) .eq. 8) jatom=1
        f_cord(i,:,iatom)=f_cord(i,:,iatom)+dist_factor(i,j,:)*f_non_add(i,j)
        f_cord(j,:,jatom)=f_cord(j,:,jatom)-dist_factor(i,j,:)*f_non_add(i,j)
      enddo
      enddo


      end
