      subroutine Rg_bias(pro_cord,f_cord,E,tempav)


!     calculates the contribution to forces (ie f_cord) 
!     and energies due to
!     a potential that is polynomial in radius of gyration
!     (note Rg calc from C-alpha carbons only)


      use amhglobals,  only:SO, maxsiz,maxcrd,n_Rg_bias,Rg_biaspoly,nmres,rg_bounds,&
           i_rg_corey, rg_shift, rg_scl,&
           i_rg_garyk, i_rg_first, D_rg, T_rg, delR_rg, M_rg, kappa_rg, &
           oarchv

     
!     argument declarations

      implicit none
      
      logical, intent(in):: tempav
       double precision, intent(in), dimension(maxsiz,3,maxcrd)::  pro_cord
       double precision, intent(out), dimension(maxsiz,3,maxcrd):: f_cord
       double precision, intent(out):: E(:,:)

!     internal variables
      
      integer :: i_axis,i_res,i
      double precision, dimension(3):: r_cm=0.0D0
      double precision :: Rg,V,dV_dRg,factor,rg_predicted
      double precision :: alpha, beta, D, D_critical, M_rg_critical
      double precision ::  Rg_max, delR_rg_critical,  V_max

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     zero force and energy

      f_cord=0.0D0
      E=0.0D0
      D=0.0D0
     
!           Calculate centre of mass and radius of gyration


      
      do i_axis=1,3
        r_cm(i_axis)=0.0D0
        do i_res=1,nmres
          r_cm(i_axis)=r_cm(i_axis)+pro_cord(i_res,i_axis,1)
!         write(SO,*) i_res,pro_cord(i_res,i_axis,1),r_cm(i_axis)
        enddo
        r_cm(i_axis)=r_cm(i_axis)/dble(nmres)
!        write(SO,*) 'final c of m',i_axis,r_cm(i_axis)
      enddo
!       write(SO,*) 'final c of m',(r_cm(i_axis),i_axis=1,3)

      Rg=0.0D0
      do i_axis=1,3
      do i_res=1,nmres
        Rg=Rg+(pro_cord(i_res,i_axis,1)-r_cm(i_axis))**2
      enddo
      enddo
      Rg=dsqrt(Rg/dble(nmres))

      rg_predicted = rg_shift*2.2D0*(dble(nmres)**0.38D0)
      
      if(i_rg_garyk) then
         ! to request a Maple worksheet for this potential, send an email to gpapoian@unc.edu 
         ! V=(D+alpha*(R-R_0)^2)/(1+beta*(R-R_0)^4), minimum at R_0, goes to 0 at large R
         !i_rg_garyk, D_rg, T_rg, delR_rg, M_rg, kappa_rg
         D=D_rg*dble(nmres)
         ! Two strategies to deal for the pole which develop for shallow D values: 
         ! 1. to increase delR_rg, i.e. allowing larger Rg fluctuations, or  
         ! 2. to decrease M_rg, i.e. decreasing the capture radius of the potential
         ! An unphysically large value of M_rg (>100 instead of ~3) in the input file 
         ! indicates that it needs to be calculated here, otherwise delR_rg is adjusted
         if(M_rg.lt.100) then ! Physical M_rg, the user wants to adjust delR_rg, if necessary
            delR_rg_critical=dsqrt((T_rg*(1+M_rg**2-2*M_rg))/(2*D*(kappa_rg-1)))
            if(delR_rg.lt.delR_rg_critical) delR_rg=delR_rg_critical+0.01D0
         endif
         M_rg_critical=1+dsqrt((2.0D0*kappa_rg*delR_rg**2*D-2.0D0*D*delR_rg**2)/(T_rg))
         if ((M_rg+0.1D0).gt.M_rg_critical) M_rg=M_rg_critical-0.1D0
         alpha=T_rg/(2.0D0*(delR_rg**2)*(rg_predicted**2))
         beta=-1.0D0*(2.0D0*kappa_rg*delR_rg**2*D-2.0D0*D*delR_rg**2-T_rg*M_rg**2+2.0D0*T_rg*M_rg-T_rg)&
              /(2.0D0*kappa_rg*delR_rg**2*D*rg_predicted**4 &
                *(M_rg**4-4.0D0*M_rg**3+6.0D0*M_rg**2-4.0D0*M_rg+1))
         Rg_max=(2.0D0*alpha*beta*rg_predicted+2* &
                 dsqrt(-1.0D0*alpha*beta**2*D+alpha*beta*dsqrt(beta**2*D**2+alpha**2*beta)))&
                 /(2.0D0*alpha*beta)
         D_critical=-1.0D0*(T_rg*(1+M_rg**2-2.0D0*M_rg))/(2.0D0*(1-kappa_rg)*delR_rg**2)
         V_max=(D+alpha*(Rg_max-rg_predicted)**2)/(1+beta*(Rg_max-rg_predicted)**4)
         if(i_rg_first) then
            write(oarchv,99) 'Rg_bias: D_rg, delR_rg, alpha, beta, M_rg, Rg_min, Rg_max', & 
                              D_rg, delR_rg, alpha, beta, M_rg, rg_predicted, Rg_max
            write(oarchv,98) 'Rg_bias: V_min-V_max=', D-V_max
98          format(a,2x,f16.8)
            if(Rg_max<3.0D0*rg_predicted) then 
               write (oarchv,99)'Rg_bias Warning: Capture radius may be too short, Rg_max/Rg_min=',Rg_max/rg_predicted
               write (oarchv,99) 'Rg_bias Warning: Consider setting M_rg above 100'
            endif
99          format(a,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)
            i_rg_first=.false.
         endif
         if(D.gt.D_critical) then
            write(SO,100) 'Serious Bug: by design D_rg can not be less than D_rg_critical', D/real(nmres), D_critical/real(nmres)
100         format(a,2x,f12.8,f12.8)
            stop
         endif
      endif
         
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!           Calculate V(Rg) (if needed) 

      if (tempav) then
         V=0.0D0
         if (i_rg_corey) then
            if (Rg/rg_predicted .lt. rg_bounds(1)) then
               V=rg_scl*10.0D0*(rg_bounds(1)*rg_predicted - rg_predicted)**2
               dV_dRg=0.0D0      
            elseif (Rg/rg_predicted .gt. rg_bounds(2)) then
               V=rg_scl*10.0D0*(rg_bounds(2)*rg_predicted - rg_predicted)**2
               dV_dRg=0.0D0      
            else
               V=rg_scl*10.0D0*(Rg - rg_predicted)**2
               dV_dRg= rg_scl*20.0D0*(Rg - rg_predicted)
            endif
         elseif(i_rg_garyk) then
            if(Rg.lt.Rg_max) then
               V=(D+alpha*(Rg-rg_predicted)**2)/(1+beta*(Rg-rg_predicted)**4)
               dV_dRg=(2*alpha*(Rg-rg_predicted))/(1+beta*(Rg-rg_predicted)**4)-&
                    (4*beta*(D+alpha*(Rg-rg_predicted)**2)*(Rg-rg_predicted)**3)/ &
                                      (1+beta*(Rg-rg_predicted)**4)**2
            else
               V=(D+alpha*(Rg_max-rg_predicted)**2)/ &
                         (1+beta*(Rg_max-rg_predicted)**4)
               dV_dRg=0.0D0
            endif
            write(SO,1000) 'Rg, rg_predicted, Rg_max, V, dV_dRg', &
                          Rg, rg_predicted, Rg_max, V, dV_dRg
1000   format(a,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3)
         else
            do i=1,n_Rg_bias
               V=V+Rg_biaspoly(i)*Rg**i
            enddo
            dV_dRg=0.0D0      
            do i=1,n_Rg_bias
               dV_dRg=dV_dRg+dble(i)*Rg_biaspoly(i)*Rg**(i-1)
            enddo
         endif
         E(1,12)=E(1,12)+V
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate dV(Rg)/dRg  


      dV_dRg=0.0D0      
      if (i_rg_corey) then
         if ( (Rg/rg_predicted .gt. rg_bounds(1)) .and.&
              (Rg/rg_predicted .lt. rg_bounds(2)) ) then
            dV_dRg= rg_scl*20.0D0*(Rg - rg_predicted)
         endif
      elseif(i_rg_garyk) then
         if(Rg.lt.Rg_max) then
            dV_dRg=(2*alpha*(Rg-rg_predicted))/(1+beta*(Rg-rg_predicted)**4)-&
                 (4*beta*(D+alpha*(Rg-rg_predicted)**2)*&
                    (Rg-rg_predicted)**3)/(1+beta*(Rg-rg_predicted)**4)**2
         endif
      else
         do i=1,n_Rg_bias
            dV_dRg=dV_dRg+dble(i)*Rg_biaspoly(i)*Rg**(i-1)
         enddo
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        calculate contribution to force, and add to f_cord

      factor=-dV_dRg/(dble(nmres)*Rg)
!     write(SO,*) (f_cord(50,i_axis,1),i_axis=1,3)
      do i_axis=1,3
      do i_res=1,nmres
        f_cord(i_res,i_axis,1)=f_cord(i_res,i_axis,1)+          &
                 (pro_cord(i_res,i_axis,1)-r_cm(i_axis))*factor
      enddo
      enddo


      end
