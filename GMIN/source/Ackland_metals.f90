!***********************************************************************
!SUBROUTINES:                                                          *  
!            BUILD_RHO_SITE                                            *
!            BUILD_V_SITE	                                       *
!***********************************************************************
! FUNCTIONS:  rho_pot,rho_pot_d,rho_pot_dd                             *
!             Vpot,Vpot_d,Vpot_dd				       *             
!             Fembed,Fembed_d,Fembed_dd                                *
!             fcut, fcut_2,fcut_dd                                     *
!             Mfunc,Mfunc_d, Mfunc_dd,                                 *
!             gfunc,gfunc_d,gfunc_dd                                   *
!             Hfunc                                                    *
!             delta_dirac                                              *
!                                                                      *             
!****|******************************************************************|

!****|******************************************************************|
        double precision function rho_pot(ipot,R)
	implicit double precision (a-h,o-z)
	double precision Rc
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h' 
        COMMON /param_cut_off/Rc	
	integer ipot		
	double precision  R
	
!debug       	write(*,*) '================================'
!debug	write(*,*) Rc
!debig	stop 
		
        if(ipot.eq.1) then   	
           rho_pot=exp(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)	   
        elseif(ipot.eq.2) then	
           rho_pot=(Ro/R)**(2*q)*fcut(ipot,R,Rc,delta)		   
        elseif(ipot.eq.3) then	
            rho_pot=exp(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)		   
	elseif(ipot.eq.4) then			
	rho_pot=gfunc(R,a,beta1,beta2,R03,R04)*fcut(ipot,R,Rc,h)	
	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	rho_pot = fpsi(R)
	elseif ((ipot.eq.12).or.(ipot.eq.13)) then
	rho_pot =  (   0.77718711248373d0 * (5.6d0-R)**4 &
	            -  0.48102928454986d0 * (5.6d0-R)**5 &  
		    +  0.14501312593993d0 * (5.6d0-R)**6 &
		    -  0.021292226813959d0* (5.6d0-R)**7 &
		    +  0.001220921762567d0* (5.6d0-R)**8) * Hfunc(5.6d0-R)
	else
	write(*,*) 'erreur de ipot'	
	
	endif 
	   
        Return
        End    

!****|******************************************************************|	
!****|******************************************************************|	
       double precision function rho_pot_d(ipot,R)
	implicit double precision (a-h,o-z)  
        double precision Rc     
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h' 
	COMMON /param_cut_off/Rc
	double precision R		
	integer ipot	      


        if(ipot.eq.1.) then     	   	     	        
      rho_pot_d=-(2*q)/Ro*EXP(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)+     &
               EXP(-2*q*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)
     
     
        elseif(ipot.eq.2) then
    
        rho_pot_d=-2*q/R*(Ro/R)**(2*q)*fcut(ipot,R,Rc,delta)+               &
               (Ro/R)**(2*q)*fcut_d(ipot,R,Rc,delta)	
     
        elseif(ipot.eq.3) then 
	   	     	        
      rho_pot_d=-(2*q)/Ro*EXP(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)+     &
               EXP(-2*q*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)	
	
        elseif(ipot.eq.4) then 	
   			    
       rho_pot_d=gfunc_d(R,a,beta1,beta2,R03,R04)*fcut(ipot,R,Rc,h)    &
     	        +gfunc(R,a,beta1,beta2,R03,R04)*fcut_d(ipot,R,Rc,h)   
	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	
	rho_pot_d=fpsi_d(R)
	elseif ((ipot.eq.12).or.(ipot.eq.13)) then
	rho_pot_d =  (- 0.77718711248373d0 * 4.d0* (5.6d0-R)**3 & 
	              + 0.48102928454986d0 * 5.d0* (5.6d0-R)**4 &
	              - 0.14501312593993d0 * 6.d0* (5.6d0-R)**5 &
		      + 0.021292226813959d0* 7.d0* (5.6d0-R)**6 &
		      - 0.001220921762567d0* 8.d0* (5.6d0-R)**7) * Hfunc(5.6d0-R)
	else
	write(*,*) 'erreur de ipot'
	endif
		
		        
        Return
        End  
!****|******************************************************************|	
        double precision  function rho_pot_dd(ipot,R)
	implicit  double precision(a-h,o-z)  
	double precision Rc      
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h' 
	COMMON /param_cut_off/Rc	
	double precision R		
	integer ipot	      


       if(ipot.eq.1) then   
       rho_pot_dd=(2*q/Ro)**2*EXP(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)        &
              -2.0d0*(2*q/Ro)*EXP(-2*q*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)      &
                             +EXP(-2*q*(R/Ro-1.0d0))*fcut_dd(ipot,R,Rc,delta)      
	     
        elseif(ipot.eq.2) then
        rho_pot_dd=2*q*(2*q+1)/R**2*(Ro/R)**(2*q)*fcut(ipot,R,Rc,delta)          &
                  -   2*(2*q)/R*(Ro/R)**(2*q)*fcut_d(ipot,R,Rc,delta)            &
                   + (Ro/R)**(2*q)*fcut_dd(ipot,R,Rc,delta)
       
	
        elseif(ipot.eq.3) then
       rho_pot_dd=(2*q/Ro)**2*EXP(-2*q*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)         &
              -2.0d0*(2*q/Ro)*EXP(-2*q*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)       &
                             +EXP(-2*q*(R/Ro-1.0d0))*fcut_dd(ipot,R,Rc,delta)           
	          	   	
        elseif(ipot.eq.4) then 	
!****|******************************************************************|	   			    
       rho_pot_dd=gfunc_dd(R,a,beta1,beta2,R03,R04)*fcut(ipot,R,Rc,h)     &    
             +2*gfunc_d(R,a,beta1,beta2,R03,R04)*fcut_d(ipot,R,Rc,h)      &
                + gfunc(R,a,beta1,beta2,R03,R04)*fcut_dd(ipot,R,Rc,h)  
!****|******************************************************************|	   
	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	rho_pot_dd = fpsi_dd(R)
	
	elseif ((ipot.eq.12).or.(ipot.eq.13)) then
	rho_pot_dd =  (  0.77718711248373d0 *12.d0* (5.6d0-R)**2 & 
	               - 0.48102928454986d0 *20.d0* (5.6d0-R)**3 &
	               + 0.14501312593993d0 *30.d0* (5.6d0-R)**4 &
		       - 0.021292226813959d0*42.d0* (5.6d0-R)**5 &
		       + 0.001220921762567d0*56.d0* (5.6d0-R)**6) * Hfunc(5.6d0-R)
	else
	write(*,*) 'erreur de ipot'
	endif	        
        Return
        End  	
	
!****|******************************************************************|
        double precision function Vpot(ipot,R)
	implicit double precision (a-h,o-z)
	double precision Rc
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h' 
	COMMON /param_cut_off/Rc	
	integer ipot
	double precision Mfunc,R	
	
        if(ipot.eq.1) then   	
           Vpot=exp(-p*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)
           Vpot=a_r*Vpot	   	    
        elseif(ipot.eq.2) then	
           Vpot=(Ro/R)**p*fcut(ipot,R,Rc,delta)
           Vpot=a_r*Vpot	   		   
        elseif(ipot.eq.3) then		
           Vpot=(Ro/R)**p*fcut(ipot,R,Rc,delta)	
           Vpot=a_r*Vpot	   			   	   
	elseif(ipot.eq.4) then	
		
	Vpot=(E1*Mfunc(R,R01,alpha1)+E2*Mfunc(R,R02,alpha2)+dd)*      &
             fcut(ipot,R,Rc,h)                                        &
            -S1*Hfunc(Rs1-R)*(Rs1-R)**4                               &
            -S2*Hfunc(Rs2-R)*(Rs2-R)**4                               &   
            -S3*Hfunc(Rs3-R)*(Rs3-R)**4	
	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	
	Vpot = 0.5d0*fvarphi(R)
	
	elseif (ipot.eq.12) then
	
	Vpot = exp(  12.33339230761400d0	         &
	           - 10.84732196908600d0*R	         &
	           +  4.57335244245080d0*R**2            &
		   -  0.85266291445935d0*R**3)*            Hfunc(2.3d0-R) * Hfunc0(R-1.d0)  & 
		    
	         +(- 14.261501929757d0*   (3.5d0-R)**4   &
		   + 15.850036758176d0*   (3.5d0-R)**5   &
	           - 11.325102264291d0*   (3.5d0-R)**6   &
		   - 4.0971114831366d0*   (3.5d0-R)**7   &
	           + 3.6739378016909d0*   (3.5d0-R)**8  ) * Hfunc(3.5d0-R)*Hfunc(R-2.3d0)  &
		 +(  1.3066813393823d0*   (6.0d0-R)**4   &
		   - 0.60542710718094d0*  (6.0d0-R)**5   &
		   + 1.0055527194350d0 *  (6.0d0-R)**6   &
		   - 0.14918186777562d0*  (6.0d0-R)**7   &
		   + 0.032773112059590d0* (6.0d0-R)**8  ) * Hfunc(6.d0-R)*Hfunc(R-2.3d0)   &
		 +(  0.011433120304691d0* (7.6d0-R)**4   &
		   - 0.021982172508973d0* (7.6d0-R)**5   &
		   - 0.012542439692607d0* (7.6d0-R)**6   &
		   + 0.025062673874258d0* (7.6d0-R)**7   &
		   - 0.0075442887837418d0*(7.6d0-R)**8  ) * Hfunc(7.6d0-R)*Hfunc(R-2.3d0) 

	Vpot=0.5d0*Vpot
	elseif (ipot.eq.13) then
        
	Vpot = exp(  12.8822300381920d0    		     &
	           - 12.1838501578140d0*R  		     &
	           +  5.5998956281737d0*R**2		     &
		   -  1.0915156420318d0*R**3)               *  Hfunc(2.3d0-R)*Hfunc0(R-1.d0)  & 
		   
		 +(   8.4670497139946d0*     (3.5d0-R)**4    &
		   - 46.183472786003d0*      (3.5d0-R)**5    &
		   + 79.633499844770d0*      (3.5d0-R)**6    &
		   - 64.847634731465d0*      (3.5d0-R)**7    &
		   + 19.454623850774d0*      (3.5d0-R)**8  ) * Hfunc(3.5d0-R)*Hfunc(R-2.3d0)  & 
		   
		 +(-  0.097845860135187d0*   (6.0d0-R)**4    &
		   -  0.47537134413743d0*    (6.0d0-R)**5    &
		   -  0.00096806164225329d0* (6.0d0-R)**6    &
		   -  0.16355187497617d0*    (6.0d0-R)**7    &
		   -  0.00090914903435333d0* (6.0d0-R)**8  ) * Hfunc(6.d0-R)*Hfunc(R-2.3d0)   &   
		     
		 +(-  0.022038480751134d0*   (7.6d0-R)**4    &
		   -  0.060955465943384d0*   (7.6d0-R)**5    &
		   +  0.11573689045653d0*    (7.6d0-R)**6    &
		   -  0.062697675088029d0*   (7.6d0-R)**7    &
		   +  0.011273545085049d0*   (7.6d0-R)**8  ) * Hfunc(7.6d0-R)*Hfunc(R-2.3d0) 

        Vpot=0.5d0*Vpot
	else
	write(*,*) 'erreur de ipot'	
	
	endif    
	   
        Return
        End    
!****|******************************************************************|
        double precision function Vpot_d(ipot,R)
	implicit double precision (a-h,o-z)	
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h'	
	double precision Rc
	COMMON /param_cut_off/Rc	 
	double precision Mfunc,Mfunc_d,R		
	integer ipot	      


        if(ipot.eq.1) then
      	   	     	        
        Vpot_d=-p/Ro*DEXP(-p*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)+    &
               DEXP(-p*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)
        Vpot_d=a_r*Vpot_d     
     
        elseif(ipot.eq.2) then
    
	Vpot_d=-p/R*(Ro/R)**p*fcut(ipot,R,Rc,delta)+                 &
               (Ro/R)**p*fcut_d(ipot,R,Rc,delta)
        Vpot_d=a_r*Vpot_d     
        elseif(ipot.eq.3) then
    
	Vpot_d=-p/R*(Ro/R)**p*fcut(ipot,R,Rc,delta)+                 &
               (Ro/R)**p*fcut_d(ipot,R,Rc,delta)        	
        Vpot_d=a_r*Vpot_d    
        elseif(ipot.eq.4) then 

      Vpot_d=(E1*Mfunc_d(R,R01,alpha1)+E2*Mfunc_d(R,R02,alpha2))*    &
            fcut(ipot,R,Rc,h)                                        &
            +4*S1*Hfunc(Rs1-R)*(Rs1-R)**3                            &
            +4*S2*Hfunc(Rs2-R)*(Rs2-R)**3                            &    
            +4*S3*Hfunc(Rs3-R)*(Rs3-R)**3                            &
      +    (E1*Mfunc(R,R01,alpha1)+E2*Mfunc(R,R02,alpha2)+dd)*       &
            fcut_d(ipot,R,Rc,h)     		 	   
	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	
	Vpot_d = 0.5d0*fvarphi_d(R)
     
	elseif (ipot.eq.12) then
	
	Vpot_d= exp( 12.333392307614d0	   -  10.847321969086d0*R	    &
	            +4.5733524424508d0*R**2 -  0.85266291445935d0*R**3)     &
		*( -  10.847321969086d0	   +   4.5733524424508d0*2.d0*R     &
		   -  0.85266291445935d0*3.d0*R**2)                         &   
		    * Hfunc(2.3d0-R) * Hfunc(R-1.d0)                        &  
	         +(  14.261501929757d0*   4.d0*  (3.5d0-R)**3   &
		   - 15.850036758176d0*   5.d0*  (3.5d0-R)**4   &
	           + 11.325102264291d0*   6.d0*  (3.5d0-R)**5   &
		   + 4.0971114831366d0*   7.d0*  (3.5d0-R)**6   &
	           - 3.6739378016909d0*   8.d0*  (3.5d0-R)**7  ) * Hfunc(3.5d0-R)*Hfunc(R-2.3d0)  &
		 +(- 1.3066813393823d0*   4.d0*  (6.0d0-R)**3   &
		   + 0.60542710718094d0*  5.d0*  (6.0d0-R)**4   &
		   - 1.0055527194350d0*   6.d0*  (6.0d0-R)**5   &
		   + 0.14918186777562d0*  7.d0*  (6.0d0-R)**6   &
		   - 0.032773112059590*   8.d0*  (6.0d0-R)**7  ) * Hfunc(6.d0-R)*Hfunc(R-2.3d0)	&
		 +(- 0.011433120304691d0* 4.d0*  (7.6d0-R)**3   &
		   + 0.021982172508973d0* 5.d0*  (7.6d0-R)**4   &
		   + 0.012542439692607d0* 6.d0*  (7.6d0-R)**5   &
		   - 0.025062673874258d0* 7.d0*  (7.6d0-R)**6   &
		   + 0.0075442887837418d0*8.d0*  (7.6d0-R)**7  ) * Hfunc(7.6d0-R)*Hfunc(R-2.3d0) 

	Vpot_d=0.5d0*Vpot_d
	elseif (ipot.eq.13) then
        
	Vpot_d = exp( 12.882230038192d0      -   12.183850157814*R               &
	           +5.5998956281737d0*R**2 -    1.0915156420318d0*R**3)        &
		 *(-12.183850157814        +    5.5998956281737d0*2.d0*R       &
		   -1.0915156420318d0*3.d0*R**2)                               &
		   * Hfunc(2.3d0-R) * Hfunc0(R-1.d0)                            &
		 +(- 8.4670497139946d0*     4.d0*  (3.5d0-R)**3    &
		   + 46.183472786003d0*     5.d0*  (3.5d0-R)**4    &
		   - 79.633499844770d0*     6.d0*  (3.5d0-R)**5    &
		   + 64.847634731465d0*     7.d0*  (3.5d0-R)**6    &
		   - 19.454623850774d0*     8.d0*  (3.5d0-R)**7  ) * Hfunc(3.5d0-R)*Hfunc(R-2.3d0)  & 
		 +(+ 0.097845860135187d0*   4.d0*  (6.0d0-R)**3    &
		   + 0.47537134413743d0*    5.d0*  (6.0d0-R)**4    &
		   + 0.00096806164225329d0* 6.d0*  (6.0d0-R)**5    &
		   + 0.16355187497617d0*    7.d0*  (6.0d0-R)**6    &
		   + 0.00090914903435333d0* 8.d0*  (6.0d0-R)**7  ) * Hfunc(6.d0-R)*Hfunc(R-2.3d0)   &	  
		 +(+ 0.022038480751134d0*   4.d0*  (7.6d0-R)**3    &
		   + 0.060955465943384d0*   5.d0*  (7.6d0-R)**4    &
		   - 0.11573689045653d0*    6.d0*  (7.6d0-R)**5    &
		   + 0.062697675088029d0*   7.d0*  (7.6d0-R)**6    &
		   - 0.011273545085049d0*   8.d0*  (7.6d0-R)**7  ) * Hfunc(7.6d0-R)*Hfunc(R-2.3d0) 
	Vpot_d=0.5d0*Vpot_d

	else
	write(*,*) 'erreur de ipot'
	endif	                         
    
        Return
        End    
	
!****|******************************************************************|
        double precision function Vpot_dd(ipot,R)
	implicit double precision (a-h,o-z)	
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h' 
	double precision Rc
	COMMON /param_cut_off/Rc	
	double precision Mfunc,Mfunc_d,Mfunc_dd,R		
	integer ipot	      


        if(ipot.eq.1) then   
       Vpot_dd=(p/Ro)**2*EXP(-p*(R/Ro-1.0d0))*fcut(ipot,R,Rc,delta)      &
          -2.0d0*(p/Ro)*EXP(-p*(R/Ro-1.0d0))*fcut_d(ipot,R,Rc,delta)     &
              +EXP(-p*(R/Ro-1.0d0))*fcut_dd(ipot,R,Rc,delta)
         Vpot_dd=a_r*Vpot_dd   
	     
        elseif(ipot.eq.2) then
        Vpot_dd=p*(p+1)/R**2*(Ro/R)**p*fcut(ipot,R,Rc,delta)          &
                      -2*p/R*(Ro/R)**p*fcut_d(ipot,R,Rc,delta)        &
                            +(Ro/R)**p*fcut_dd(ipot,R,Rc,delta)       
        Vpot_dd=a_r*Vpot_dd  
	
        elseif(ipot.eq.3) then
       Vpot_dd=p*(p+1)/R**2*(Ro/R)**p*fcut(ipot,R,Rc,delta)           &
                      -2*p/R*(Ro/R)**p*fcut_d(ipot,R,Rc,delta)        &
                            +(Ro/R)**p*fcut_dd(ipot,R,Rc,delta)       
         Vpot_dd=a_r*Vpot_dd  
	          	  
        elseif(ipot.eq.4) then 

!****|******************************************************************|	   			    
      Vpot_dd=(E1*Mfunc_dd(R,R01,alpha1)+E2*Mfunc_dd(R,R02,alpha2))*fcut(ipot,R,Rc,h)       &
       +    2*(E1*Mfunc_d(R,R01,alpha1)+E2*Mfunc_d(R,R02,alpha2))*fcut_d(ipot,R,Rc,h)       &
       +     (E1*Mfunc(R,R01,alpha1)+E2*Mfunc(R,R02,alpha2)+dd)*fcut_dd(ipot,R,Rc,h)        &   
       -      12*S1*Hfunc(Rs1-R)*(Rs1-R)**2            					    &
       -      12*S2*Hfunc(Rs2-R)*(Rs2-R)**2       					    &
       -      12*S3*Hfunc(Rs3-R)*(Rs3-R)**2      					          
     		 
!****|******************************************************************|	   
 	elseif ((ipot.eq.5).or.(ipot.eq.6)) then
	
	Vpot_dd = 0.5d0*fvarphi_dd(R)
	
	else
	write(*,*) 'erreur de ipot'
	endif	                         
    
        Return
        End    		
	
	  
!****|******************************************************************|	
	double precision  function fcut(ipot,R,Rc,width)
	implicit double precision (a-h,o-z)	
	double precision R,Rc,width,test,x
	
	if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then
           test=(R-Rc)/width
	   if(test.lt.-100.d0)  then
           fcut=1.0d0	
	   elseif(test.gt.100.d0) then
           fcut=0.0d0	
	   else	
	   x=test	
           fcut=1.0d0/(1.0d0+EXP(x))
	   endif
	elseif(ipot.eq.4) then	   
	    test=(R-Rc)
	   if(test.gt.0.d0)  then
           fcut=0.0d0
	   elseif(test.le.0.d0) then
	   x=test/width	   
           fcut= x**4/(1+x**4)	   	   
	   endif
	elseif ((ipot.eq.5).or.(ipot.eq.6).or.(ipot.eq.12).or.(ipot.eq.13)) then
	fcut=1.d0   		   
	endif
	   			
        Return
        End   	
	
!****|******************************************************************|	
	double precision  function fcut_d(ipot,R,Rc,width)
	implicit double precision (a-h,o-z)	
	double precision R,Rc,width,test,x
	
	if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then	
          test=(R-Rc)/width
	   if(test.lt.-100.d0)  then
             fcut_d=0.0d0	
	   elseif(test.gt.100.d0) then
             fcut_d=0.0d0		
	    else
	     x=test		
             fcut_d=-EXP(x)/(1+EXP(x))**2
             fcut_d=fcut_d/width          
	   endif		
	elseif(ipot.eq.4) then	
	    test=(R-Rc)
	   if(test.gt.0.d0)  then
	   fcut_d=0.0d0	
	   elseif(test.le.0.d0) then
	   x=test/width
           fcut_d= 4*x**3/(1+x**4)**2  	 
           fcut_d=fcut_d/width  	     
	   endif		
	elseif ((ipot.eq.5).or.(ipot.eq.6).or.(ipot.eq.12).or.(ipot.eq.13)) then
	fcut_d=1.d0   		   
	endif	
		
        Return
        End   	
	   		
!****|******************************************************************|	
	double precision function fcut_dd(ipot,R,Rc,width)
	implicit  double precision(a-h,o-z)	
	double precision R,Rc,width,test,x
	
	if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then	
          test=(R-Rc)/width
	   if(test.lt.-100.d0)  then
             fcut_dd=0.0d0	
	   elseif(test.gt.100.d0) then
             fcut_dd=0.0d0		
	    else
	     x=test		
             fcut_dd= -EXP(x)/( (1.0d0+EXP(x) ) )**2 +2.0d0*EXP(2*x)/(1.0d0+EXP(x) )**3    
           fcut_dd=fcut_dd/width**2     
	   endif		
	elseif(ipot.eq.4) then	
	    test=(R-Rc)
	   if(test.gt.0.d0)  then
	   fcut_dd=0.0d0	
	   elseif(test.le.0.d0) then
	   x=test/width	   
           fcut_dd=12*x**2/(1+x**4)**2-32*x**6/(1+x**4)**3
           fcut_dd=fcut_dd/width**2 	      	   
	   endif		
	elseif ((ipot.eq.5).or.(ipot.eq.6).or.(ipot.eq.12).or.(ipot.eq.13)) then
	fcut_dd=1.d0   		   
	
	endif	
		
        Return
        End  	
		   	 	 	
!****|******************************************************************|
        double precision function Fembed(ipot,x)
        USE COMMONS, ONLY : ACK1, ACK2
	implicit double precision (a-h,o-z)	
        double precision x
	integer ipot
	double precision Rc
	COMMON /param_cut_off/Rc	
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h'
	include 'ackland_mendelev_fe.h'	

		
	  if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then	
             Fembed=x**(alpha)
             Fembed=-b_a*Fembed	  
	        
	  elseif(ipot.eq.4) then	     
	  
	  if(x.le.1) then	  	  
	   Fembed=F0+0.5*F2*(x-1)**2+q1*(x-1)**3+q2*(x-1)**4 +q3*(x-1)**5 +q4*(x-1)**6 
        
	  elseif(x.gt.1) then	  	  	
	   Fembed=(F0+0.5*F2*(x-1)**2+q1*(x-1)**3+qq1*(x-1)**4)/(1+qq2*(x-1)**3)    	   
	  endif	   	  
	  	  
	  elseif(ipot.eq.5) then
	  
	  Fembed = -dsqrt(x) + aphi*x*x + aphi2*x*x*x*x  
	  
	  elseif(ipot.eq.6) then
	   if (x.lt.1.d0) then
	     Fembed = -aphi*dsqrt(x) - aphi2*(1-dsqrt(x))*dlog(2.d0-x)/dlog(2.d0) 
	     else
	     Fembed = -aphi*dsqrt(x) 
	   endif

	  elseif(ipot.eq.12) then
	   Fembed = -     dsqrt(x)                                         &
	            - 1.9162462126235d0*1.d-7*(x-60.d0)**4*Hfunc(x-60.d0) &
	            + 4.6418727035037d0*1.d-7*(x-70.d0)**4*Hfunc(x-70.d0) &
	            + 6.6448294272955d0*1.d-7*(x-80.d0)**4*Hfunc(x-80.d0) &
		    - 2.0680252960229d0*1.d-6*(x-85.d0)**4*Hfunc(x-85.d0) &
		    + 1.1387131464983d0*1.d-6*(x-90.d0)**4*Hfunc(x-90.d0)
	  elseif(ipot.eq.13) then
	   Fembed = - dsqrt(x)                                         &
	            + 3.2283012597866d0*1.d-7*(x-60.d0)**4*Hfunc(x-60.d0) &		
		    - 1.1552813894483d0*1.d-6*(x-70.d0)**4*Hfunc(x-70.d0) &
		    + 2.3747280268355d0*1.d-6*(x-80.d0)**4*Hfunc(x-80.d0) &
		    - 2.0379550826523d0*1.d-6*(x-85.d0)**4*Hfunc(x-85.d0) & 
		    + 4.9758343293936d0*1.d-7*(x-90.d0)**4*Hfunc(x-90.d0)


	 end if 
	return
        end
!****|******************************************************************|
        double precision function Fembed_d(ipot,x)
        USE COMMONS, ONLY : ACK1, ACK2
	implicit double precision (a-h,o-z)	
	double precision x
	double precision Rc
	COMMON /param_cut_off/Rc	
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h'
	include 'ackland_mendelev_fe.h'
	
          
	  if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then	
             Fembed_d=alpha*x**(alpha-1.0d0)
             Fembed_d=-b_a*Fembed_d	     	  
	  elseif(ipot.eq.4) then	     	  
	  if(x.le.1) then	  	  
	   Fembed_d=F2*(x-1)+3*q1*(x-1)**2+4*q2*(x-1)**3+5*q3*(x-1)**4+6*q4*(x-1)**5 
        
	  elseif(x.gt.1) then	  	  	
	 Fembed_d=(F2*(x-1)+3*q1*(x-1)**2+4*qq1*(x-1)**3)/(1+qq2*(x-1)**3)       &
                  -3*qq2*(x-1)**2*(F0+0.5*F2*(x-1)**2+q1*(x-1)**3+qq1*(x-1)**4)  &
                  /(1+qq2*(x-1)**3)**2      
          
	      	   
	   endif	   	  	  

	  elseif (ipot.eq.5) then
	  
	  Fembed_d = -0.5d0/dsqrt(x) + 2.d0*aphi*x  &
	              + 4.d0*aphi2*x*x*x


	  elseif(ipot.eq.6) then
	  if (x.lt.1.d0) then
	    Fembed_d = -aphi*0.5d0/dsqrt(x)                   &
	               + aphi2*(0.5d0*dlog(2.d0-x)/dsqrt(x)  &
	               + (1-dsqrt(x))/(2.d0-x))/dlog(2.d0) 
	   else
	    Fembed_d = -0.5d0*aphi/dsqrt(x)     
          end if 	    
	    
	  elseif(ipot.eq.12) then
	   Fembed_d = - 0.5d0 / dsqrt(x)                                 &
	            + 4.d0*(                                           &
		    - 1.9162462126235d0*1.d-7*(x-60.d0)**3*Hfunc(x-60.d0) &
	            + 4.6418727035037d0*1.d-7*(x-70.d0)**3*Hfunc(x-70.d0) &
	            + 6.6448294272955d0*1.d-7*(x-80.d0)**3*Hfunc(x-80.d0) &
		    - 2.0680252960229d0*1.d-6*(x-85.d0)**3*Hfunc(x-85.d0) &
		    + 1.1387131464983d0*1.d-6*(x-90.d0)**3*Hfunc(x-90.d0) &
		           )
	  elseif(ipot.eq.13) then
	   Fembed_d = - 0.5d0 / dsqrt(x)                                 &
	            + 4.d0*(                                           &
	            + 3.2283012597866d0*1.d-7*(x-60.d0)**3*Hfunc(x-60.d0) &		
		    - 1.1552813894483d0*1.d-6*(x-70.d0)**3*Hfunc(x-70.d0) &
		    + 2.3747280268355d0*1.d-6*(x-80.d0)**3*Hfunc(x-80.d0) &
		    - 2.0379550826523d0*1.d-6*(x-85.d0)**3*Hfunc(x-85.d0) & 
		    + 4.9758343293936d0*1.d-7*(x-90.d0)**3*Hfunc(x-90.d0) &
		           )
        
	
	endif

        return
        end
	
!****|******************************************************************|	
        double precision function Fembed_dd(ipot,x)
        USE COMMONS, ONLY : ACK1, ACK2
	implicit double precision (a-h,o-z)	
        include 'ackland_sma.h'  
        include 'ackland_mishin_cu.h'
	double precision Rc
	COMMON /param_cut_off/Rc	
	double precision x,denom
	include 'ackland_mendelev_fe.h'
          
	  if(ipot.eq.1.or.ipot.eq.2.or.ipot.eq.3) then	
             Fembed_dd=alpha*(alpha-1)*x**(alpha-2.0d0)
             Fembed_dd=-b_a*Fembed_dd	     	  
	  elseif(ipot.eq.4) then	     	  
	  if(x<1) then
	   Fembed_dd=F2+6*q1*(x-1)+12*q2*(x-1)**2+20*q3*(x-1)**3+30*q4*(x-1)**4 
        
	  elseif(x>1) then
	 denom=1+qq2*(x-1)**3
	 	  	  	  	
	 Fembed_dd=(F2+6*q1*(x-1)+12*qq1*(x-1)**2)/denom                              &
         -6*qq2*(x-1)**2*(F2*(x-1)+3*q1*(x-1)**2+4*qq1*(x-1)**3)/denom**2             &
         -6*qq2*(x-1)*(F0+0.5*F2*(x-1)**2+q1*(x-1)**3+qq1*(x-1)**4)/denom**2          &
          +18*qq2**2*(x-1)**4*(F0+0.5*F2*(x-1)**2+q1*(x-1)**3+qq1*(x-1)**4)/denom**3  
                        	   
	   endif	   	  	  
        	    
	  elseif (ipot.eq.5) then
	    Fembed_dd = 0.25d0/dsqrt(x**3) + 2.d0*aphi  &
	              + 12.d0*aphi2*x*x
		      
          elseif (ipot.eq.6) then
	   if (x.lt.1.d0) then
	     Fembed_dd = 0.25d0*aphi/dsqrt(x**3)  + aphi2/dlog(2.d0)*   &
	             (                                              &
		       -0.25d0*dlog(2.d0-x)/dsqrt(x**3)             &
		       -1.0/((2.d0-x)*dsqrt(x))                     &
		       +(1.0-dsqrt(x))/(2.d0-x)**2                  &
		      ) 
	     else
	      Fembed_dd = 0.25d0*aphi/dsqrt(x**3)     
            end if 
	 
        endif

        return
        end
!****|******************************************************************|
        double precision function delta_dirac(i,j)
        integer i,j

         delta_dirac=0.0d0
        if(i.eq.j) then
         delta_dirac=1.0d0
        endif
        return
        end
!****|******************************************************************|	
	double precision function Hfunc(x)
	double precision x
	
	if(x.lt.0.0) then
	  Hfunc=0.0d0
	elseif(x.ge.0.0) then	  
	  Hfunc=1.0d0	 
	 endif
        return
        end	 
!****|******************************************************************|	
	double precision function Hfunc0(x)
	double precision x
	
	if(x.le.0.0) then
	  Hfunc0=0.0d0
	elseif(x.gt.0.0) then	  
	  Hfunc0=1.0d0	 
	 endif
        return
        end	 
!****|******************************************************************|
	double precision function Mfunc(R,R0,alpha)
	double precision R,R0,alpha
	
	Mfunc=exp(-2*alpha*(R-R0))-2*exp(-alpha*(R-R0))
		
        return
        end	 
!****|******************************************************************|
	double precision function Mfunc_d(R,R0,alpha)
	double precision R,R0,alpha
	
       Mfunc_d=-2*alpha*exp(-2*alpha*(R-R0))+2*alpha*exp(-alpha*(R-R0))
		
        return
        end
!****|******************************************************************|	
        double precision function Mfunc_dd(R,R0,alpha)
	double precision R,R0,alpha
	
         Mfunc_dd=+4*alpha**2*exp(-2*alpha*(R-R0))-2*alpha**2*exp(-alpha*(R-R0))
		
        return
        end		
!****|******************************************************************|
	double precision function gfunc(x,a,b1,b2,x1,x2)
	double precision x,a,b1,b2,x1,x2
	
	gfunc=a*exp(-b1*(x-x1)**2)+exp(-b2*(x-x2))
		
        return
        end	 		
!****|******************************************************************|
	double precision function gfunc_d(x,a,b1,b2,x1,x2)
	double precision x,a,b1,b2,x1,x2
	
	gfunc_d=-2*b1*(x-x1)*a*exp(-b1*(x-x1)**2)-b2*exp(-b2*(x-x2))
		
        return
        end	
!****|******************************************************************|
	double precision function gfunc_dd(x,a,b1,b2,x1,x2)
	double precision x,a,b1,b2,x1,x2
	
	gfunc_dd=-2*a*b1*exp(-b1*(x-x1)**2)+4*a*b1**2*(x-x1)**2*exp(-b1*(x-x1)**2)   &
                        +b2**2*exp(-b2*(x-x2)) 		
        return
        end	
	 		
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|


!****|******************************************************************|
!****|******************************************************************|
       double precision function fpsi(x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z) 
       include 'ackland_mendelev_fe.h'
       
       temp = 0.d0
       do i=1,npsi
         temp = temp + ap(i)*Hfunc(rp(i)-x)*(rp(i)-x)**3
       end do
       fpsi = temp
       return
       end  
!****|******************************************************************|
       double precision function fpsi_d(x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z) 
       include 'ackland_mendelev_fe.h'
       
       temp = 0.d0
       do i=1,npsi
         temp = temp - 3.d0*ap(i)*Hfunc(rp(i)-x)*(rp(i)-x)**2
       end do
       fpsi_d = temp
       return
       end  
!****|******************************************************************|
       double precision function fpsi_dd(x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z) 
       include 'ackland_mendelev_fe.h'
       
       temp = 0.d0
       do i=1,npsi
         temp = temp + 6.d0*ap(i)*Hfunc(rp(i)-x)*(rp(i)-x)
       end do
       fpsi_dd = temp
       return
       end  
!****|******************************************************************|
!****|******************************************************************|
       double precision function fvarphi (x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z)
       include 'ackland_mendelev_fe.h'
        
	 threei   =1.d0/3.d0      
        ZnFE2   = ZnFE*ZnFE
       AU_TO_EV = hart
        AU_TO_A = abohr
         temp   = 0.d0

       rs = 0.88534d0*abohr*ZnFE**(-threei)/dsqrt(2.d0)
       rx = x/rs
       
        if (x.lt.r1) then 
          fvarphi = ZnFE2 * fphi(rx) * AU_TO_A* AU_TO_EV/x 
        else if ((x.ge.r1).and.(x.lt.r2)) then 
          fvarphi = dexp (  bFE0 + bFE1*x + bFE2*x*x + bFE3*x*x*x )
        else if (x.ge.r2) then
          do i=1,nvarphi
             temp = temp + af(i)*Hfunc(rf(i)-x)*(rf(i)-x)**3
          end do
          fvarphi = temp
        end if              
        
       return
       end 
!****|******************************************************************|
       double precision function fvarphi_d (x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z)
       include 'ackland_mendelev_fe.h'
         threei   = 1.d0/3.d0
        ZnFE2   = ZnFE*ZnFE
       AU_TO_EV = hart
        AU_TO_A = abohr
         temp   = 0.d0

       rs = 0.88534d0*abohr*ZnFE**(-threei)/dsqrt(2.d0)
       rx = x/rs

        if (x.lt.r1) then 
          fvarphi_d =  - ZnFE2*fphi(rx)*AU_TO_EV*AU_TO_A/(x*x) &
                 + ZnFE2 * fphi_d(rx) * AU_TO_EV*AU_TO_A / (x*rs)     
        else if ((x.ge.r1).and.(x.lt.r2)) then 
          fvarphi_d = (bFE1 + 2.0d0*bFE2*x + 3.0d0*bFE3*x*x) &
                     *dexp (bFE0 + bFE1*x + bFE2*x*x + bFE3*x*x*x)
        else if (x.ge.r2) then
          do i=1,nvarphi
             temp = temp - 3.d0*af(i)*Hfunc(rf(i)-x)*(rf(i)-x)**2
          end do
          fvarphi_d = temp
        end if              
        
       return
       end
!****|******************************************************************|
       double precision function fvarphi_dd (x)
       USE COMMONS, ONLY : ACK1, ACK2
       implicit double precision (a-h,o-z)
       include 'ackland_mendelev_fe.h'
         threei   = 1.d0/3.d0
        ZnFE2   = ZnFE*ZnFE
       AU_TO_EV = hart
        AU_TO_A = abohr
         temp   = 0.d0

       rs = 0.88534d0*abohr*ZnFE**(-threei)/dsqrt(2.d0)
       rx = x/rs

        if (x.lt.r1) then 
	fvarphi_dd = 0.5d0*ZnFE2*fphi(rx)*AU_TO_EV*AU_TO_A/(x*x*x) &
	             - ZnFE2*fphi_d(rx)*AU_TO_EV*AU_TO_A/(x*x*rs)  &
		     - ZnFE2 * fphi_d(rx) * AU_TO_EV*AU_TO_A / (x*x*rs) &
		     + ZnFE2 * fphi_dd(rx) * AU_TO_EV*AU_TO_A / (x*rs*rs)	 		 
        else if ((x.ge.r1).and.(x.lt.r2)) then 
          fvarphi_dd = (2.0d0*bFE2 + 6.0d0*bFE3*x) &
                     *dexp (bFE0 + bFE1*x + bFE2*x*x + bFE3*x*x*x)&
		     + &
                    (bFE1 + 2.0d0*bFE2*x + 3.0d0*bFE3*x*x)**2 &
                     *dexp (bFE0 + bFE1*x + bFE2*x*x + bFE3*x*x*x)		     
        else if (x.ge.r2) then
          do i=1,nvarphi
             temp = temp + 6.d0*af(i)*Hfunc(rf(i)-x)*(rf(i)-x)
          end do
          fvarphi_dd = temp
        end if              
        
       return
       end
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|
       double precision function fphi(x)
       double precision x
  
        fphi = 0.1818d0*dexp(-3.2d0*x) & 
            +  0.5099d0*dexp(-0.9423d0*x) &
            +  0.2802d0*dexp(-0.4029d0*x) &
            +  0.02817*dexp(-0.2016*x)

       return
       end
!****|******************************************************************|
       double precision function fphi_d(x)
       double precision x
  
       fphi_d = -  0.58176d0*dexp(-3.2d0*x) &
                -  0.480479d0*dexp(-0.9423d0*x) & 
                -  0.112893d0*dexp(-0.4029d0*x) &
                -  0.00567907d0*dexp(-0.2016*x)

       return
       end

!****|******************************************************************|
       double precision function fphi_dd(x)
       double precision x
  
       fphi_dd =    1.86163d0*dexp(-3.2d0*x) &
                +  0.452755d0*dexp(-0.9423d0*x) & 
                +  0.0454846d0*dexp(-0.4029d0*x) &
                +  0.0011449d0*dexp(-0.2016*x)

       return
       end

!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|

      SUBROUTINE BUILD_RHO_SITE(ipot,rho_site,Rn,ndir,nat_up,ndir_max)
        implicit double precision (a-h,o-z)
        double precision Rn(nat_up,ndir_max,3), rho_site(nat_up),    &
                          rho_temp,normR,R,Rtemp(3)
        integer ndir(nat_up),nat_up,ndir_max,ipot    

  
       do i=1,nat_up  
          rho_temp=0.0D0
           do ni=1,ndir(i)
	    Rtemp(:)=Rn(i,ni,:)	   
            R=sqrt(DOT_PRODUCT(Rtemp,Rtemp)) 
            rho_temp=rho_temp+rho_pot(ipot,R)	    
           end do

            rho_site(i)=rho_temp	
        end do
 
       RETURN
       END
!****|******************************************************************|       
       SUBROUTINE BUILD_V_SITE(ipot,V_site,Rn,ndir,nat_up,ndir_max)
        implicit double precision (a-h,o-z)
        double precision Rn(nat_up,ndir_max,3),V_site(nat_up),   &
                         V_temp,normR,R,Rtemp(3)
        integer ndir(nat_up),nat_up,nat_max,ndir_max,ipot

  
       do i=1,nat_up  
          V_temp=0.0D0
           do ni=1,ndir(i)	   
	    Rtemp(:)=Rn(i,ni,:)	   
            R=sqrt(DOT_PRODUCT(Rtemp,Rtemp))           	   
            V_temp=V_temp+Vpot(ipot,R)	    
           end do
              V_site(i)=V_temp	
        end do
 
       RETURN
       END      

      

