
!                             DISTANCE             
      subroutine calc_xyz(xyz_dist,xyz_unit_vect,pro_cord,nmres)
        use amhglobals,  only : maxsiz,maxcrd,tgsequences_amw
!        use globals_alt, only : OB_dns_count,outfile1_alt

        implicit none
        integer i,j,iatom,jatom
        integer, intent(in) :: nmres
        double precision xyz_dist(maxsiz,maxsiz),xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
        double precision dist_temp
        !      double precision dens(maxsiz)  !  DEBUGGING


        !      dens=0.0  !  DEBUGGING
        xyz_dist=0.0D0
        xyz_unit_vect=0.0D0
        do i=1,nmres-1   
           do j=i+1,nmres
              iatom=2
              jatom=2
!              if (ires(i) .eq. 8) iatom=1
!              if (ires(j) .eq. 8) jatom=1
!  bit of a hack only look at iseq 1              
!              if (tgsequences_amw(i,iseq) .eq. 8) iatom=1
!              if (tgsequences_amw(j,iseq) .eq. 8) jatom=1
              if (tgsequences_amw(i,1) .eq. 8) iatom=1
              if (tgsequences_amw(j,1) .eq. 8) jatom=1

              xyz_dist(i,j)=sqrt (  &
                   (pro_cord(j,1,jatom)-pro_cord(i,1,iatom))**2 + &
                   (pro_cord(j,2,jatom)-pro_cord(i,2,iatom))**2 + &
                   (pro_cord(j,3,jatom)-pro_cord(i,3,iatom))**2 )
              dist_temp=xyz_dist(i,j)
              xyz_unit_vect(i,j,:)=(pro_cord(i,:,iatom)-pro_cord(j,:,jatom))/dist_temp
              xyz_unit_vect(j,i,:)=-xyz_unit_vect(i,j,:)
              !          if((dist_temp.lt.6.5) .and. (dist_temp.ge.4.5))dens(i)=dens(i)+1.0 !  DEBUGGING
              !          if((dist_temp.lt.6.5) .and. (dist_temp.ge.4.5))dens(j)=dens(j)+1.0 !  DEBUGGING
           enddo
        enddo
!              write(6,*)'dist_temp calc xyz', dist_temp
!              write(6,*)'pro_cord 1 2  ',pro_cord(1,1,1), pro_cord(2,1,1)

        !      if(mod(OB_dns_count,1000.).lt.0.5)write(outfile1_alt(7),1000)'RESIDUE DENSITY',dens  !  DEBUGGING
        return
        ! 1000   format(a,200(1x,f6.3)) !  DEBUGGING

      end subroutine calc_xyz

      !---------------------------

      subroutine calc_theta_alt(theta, theta_dot, xyz_dist, rmin, rmax,kappa,nmres,  i_well)
        ! Calculates theta for alternative potential
        ! E = theta * [sigma(w)*gamma(w) + (1-sigma(w))*gamma(d)] 

        use amhglobals,  only : maxsiz
        use globals_alt, only : max_well_alt

        !     kappa = steepness of "stepfunction" tanh 
        implicit none
        integer i,j,i_well,nmres
        double precision kappa,t_min,t_max,rmin,rmax
        double precision theta(maxsiz,maxsiz,max_well_alt),theta_dot(maxsiz,maxsiz,max_well_alt)
        double precision xyz_dist(maxsiz,maxsiz)

        theta(:,:,i_well)=0.D0
        theta_dot(:,:,i_well)=0.D0

!C              write(6,*)'RMIN AND MAX I_WELL'rmin,rmax,i_well
               
        do i=1,nmres-2   ! work out theta values between all pairs
           do j=i+2,nmres   ! theta(ij)=k*(1+tmin)*(1+tmax))
              t_min=tanh(kappa*(xyz_dist(i,j)-rmin))
               
              t_max=tanh(kappa*(rmax-xyz_dist(i,j) ))
              theta(i,j,i_well) = 0.25D0*( 1.0D0+t_min )*( 1.0D0+t_max ) 
              theta(j,i,i_well) =  theta(i,j,i_well) 
              theta_dot(i,j,i_well)=kappa*theta(i,j,i_well)*(t_max-t_min)
              theta_dot(j,i,i_well)=theta_dot(i,j,i_well)
           enddo
!           write(6,*)'xyzdist theta calc ',xyz_dist(i,j)
        enddo

        return
      end subroutine calc_theta_alt

      !----------------------------

      subroutine calc_sum_theta_dot_alt(sum_theta_dot,theta_dot,xyz_unit_vect,nmres)
        use globals_alt, only : max_well_alt
        use amhglobals,  only : maxsiz

        implicit none
        integer, intent(in) ::  nmres
         double precision, intent(in) ::  theta_dot(maxsiz,maxsiz,max_well_alt)
         double precision, intent(in) :: xyz_unit_vect(maxsiz,maxsiz,3)
         double precision, intent(out) :: sum_theta_dot(maxsiz,3)

        integer p,k

        sum_theta_dot=0.D0

        do k=1,nmres-2,1     !calculate densities, A
           do p=k+2,nmres,1
              sum_theta_dot(k,:)=sum_theta_dot(k,:)+theta_dot(p,k,1)*xyz_unit_vect(k,p,:)
              sum_theta_dot(p,:)=sum_theta_dot(p,:)+theta_dot(k,p,1)*xyz_unit_vect(p,k,:)
           enddo
        enddo

        return
      end subroutine calc_sum_theta_dot_alt

      subroutine calc_A_alt(A,theta,nmres)

        use globals_alt, only : max_well_alt
        use amhglobals,  only : maxsiz

        implicit none
        integer i,j,nmres
        double precision theta(maxsiz,maxsiz,max_well_alt),A(maxsiz)

        A=0.D0

        do i=1,nmres-2,1     !calculate densities, A
           do j=i+2,nmres,1
              
              A(i)=A(i)+theta(i,j,1)
              A(j)=A(j)+theta(i,j,1)
!              write(6,*)'Ai Aj theta(i,j,1)',A(i),A(j),theta(i,j,1)
           enddo
        enddo

        return
      end subroutine calc_A_alt


      !--------------------------------




      subroutine calc_sigma_alt(sigma,A,treshold,kappa,nmres)
        ! Calculates sigma = weighting functions for alternative potential
        ! E = theta * [sigma(w)*gamma(w) + sigma(d)*gamma(d)] 


        use amhglobals,  only : maxsiz, CUTOFF_CONT_LOW
!        use globals_alt, only : max_well_alt

        implicit none
        integer i,j,nmres
        double precision g(maxsiz),kappa,treshold
        double precision sigma(maxsiz,maxsiz),A(maxsiz)
        double precision heaviside(maxsiz)

        sigma=0.D0

        g=kappa*(A(:) - treshold) ! calculate g

        heaviside=0.5D0*(1.0D0-tanh(g) )
!      open(unit=1383,file='alt_sigma_data',status='unknown',access='append')
!           write(1383,*)'treshold kappa nmres ' , treshold, kappa 
!           write(1383,*)'CUTOFF_CONT_LOW ', CUTOFF_CONT_LOW
!           write(1383,291) (sigma(32,j),j=1,nmres)
!           write(1383,291) (sigma(32,j),j=1,nmres)

        do i=1,nmres-CUTOFF_CONT_LOW,1
           do j=i+CUTOFF_CONT_LOW,nmres,1
              sigma(i,j)=heaviside(i)*heaviside(j)
              sigma(j,i)=sigma(i,j)
           enddo
        enddo

!      open(unit=1382,file='alt_A_data',status='unknown',access='append')  
!           write(1382,291) (A(i),i=1,nmres)

!           write(1383,291) (sigma(32,j),j=1,nmres)
!           write(1383,291) (heaviside(j),j=1,nmres)
!           291 format(110(1x,e12.6))

!              close (1382)
!             close (1383)


        return
      end subroutine calc_sigma_alt


      subroutine calc_heaviside_alt(heaviside,heaviside_dot,A,treshold,kappa)

        use amhglobals,  only : maxsiz

        implicit none
         double precision, intent(out) :: heaviside(maxsiz), heaviside_dot(maxsiz)
         double precision, intent(in) :: kappa,treshold,A(maxsiz)

        !Local
        double precision g(maxsiz)

        heaviside=0.D0
        heaviside_dot=0.D0

        g=kappa*(A(:) - treshold) ! calculate g
        heaviside=0.5D0*(1.0D0-tanh(g) )
        heaviside_dot=-0.5D0*kappa*( 1.0D0 - tanh(g)*tanh(g) ) ! Note: gdot is explicitly calculated below
        return
      end subroutine calc_heaviside_alt

