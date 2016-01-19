      subroutine default_alt()

        use globals_alt, only : altpotflag,kappa_alt,treshold_alt,kappa_well,kappa_OB, & 
             debug_numer_forces,outfile1_alt, output_step_alt, output_stepsize_alt, & 
             do_send_output,onebodyflag, onebody_type , timings_file_alt, accumulated_time

        use altpot_interfaces, only: init_onebody

        implicit none
        integer i

        treshold_alt=0.70D0
        kappa_alt=7.00D0
        kappa_well=7.00D0
        kappa_OB=7.00D0
        debug_numer_forces=.false.
        altpotflag=0
        output_step_alt=0
        output_stepsize_alt=10
        do_send_output=.false.
        onebodyflag=0
        onebody_type=1
        accumulated_time=0.00D0


        call init_onebody()

        do i=1,10,1
           outfile1_alt(i)=102+i
        enddo

!        open(unit=outfile1_alt(1),file='check_alt_file_1',status='unknown')
!        open(unit=outfile1_alt(2),file='check_alt_file_2',status='unknown')
!        open(unit=outfile1_alt(3),file='check_alt_file_3',status='unknown')
!        open(unit=outfile1_alt(4),file='check_alt_file_4',status='unknown')
!        open(unit=outfile1_alt(5),file='OB_density',status='unknown')
!        open(unit=outfile1_alt(6),file='OB_Energy',status='unknown')
!        open(unit=outfile1_alt(7),file='residue_density',status='unknown')
!        open(unit=timings_file_alt,file='alt_timings',status='unknown')

        return
      end subroutine default_alt

      ! read_input calls default_alt and read_altgamma and reads input, 
      ! Should be called after initil
      subroutine read_input_alt()

        use globals_alt, only : altpotflag,kappa_alt,treshold_alt,kappa_well,onebodyflag &
             , onebody_type,kappa_OB,Alimits_OB, debug_numer_forces

        use altpot_interfaces, only: default_alt
        use amhglobals,  only : SO

        implicit none
        integer ninput,ioer,n_output,dumi
        parameter(ninput=130)
        parameter(n_output=131)
        character str200*200,inputfile*200,paramcheckfile*200
        inputfile='input_amh'
        paramcheckfile='AMW_output'
        open(unit=ninput,file=inputfile,status='old')
!        open(unit=n_output,file=paramcheckfile,status='unknown')

        ! Initialize altpot values
        call default_alt()


        !++++++ START INPUTLOOP
!        write(6,*) 'start default alt'
100     read(ninput,'(a)',end=999)str200

        if(index(str200,'altpotflag').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) altpotflag(1),altpotflag(2)
        elseif(index(str200,'kappa_alt').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) kappa_alt
        elseif(index(str200,'treshold_alt').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) treshold_alt
        elseif(index(str200,'kappa_well').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) kappa_well
        elseif(index(str200,'onebodyflag').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) onebodyflag,onebody_type
        elseif(index(str200,'kappa_OB').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900) kappa_OB
        elseif(index(str200,'densitylimits_OB').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900)Alimits_OB
        elseif(index(str200,'debug_numer_forces').gt.0)then
           backspace(ninput)
           read(ninput,*,iostat=ioer,err=900)dumi
           if(dumi.gt.0)debug_numer_forces=.true.
        endif
        goto 100
        !------ END INPUTLOOP
!        return

        !  REPORT VALUES to file PARAMCHECK_alt:
999      write(6,*)
!999     write(n_output,*)'altpotflags are : ',altpotflag(1),altpotflag(2)
!        write(n_output,*)'kappa_alt is : ',kappa_alt
!        write(n_output,*)'treshold_alt is : ',treshold_alt
!        write(n_output,*)'kappa_well is : ',kappa_well
!1        write(n_output,*)'onebodyflag,onebody_type is : ',onebodyflag,onebody_type
!        write(n_output,*)'kappa_OB is : ',kappa_OB
!        write(n_output,*)'Alimits are : ',Alimits_OB

900     if(ioer.ne.0)then
           write(*,'(a)')'problem with reading input subroutine: read_input_alt'
        endif

        close(ninput)
!        close(n_output)

        return


      end subroutine read_input_alt

      !------------------------------------------------

      !------------------------------------

      subroutine read_altgamma()

        ! reads in altgamma and onebody_gamma from gamma.dat
        ! Format in gamma.dat is expected to be 
        ! gamma1   gamma2  well  keyword=altpot
        ! altpot is to be used for anisotropic potential
        !      use amhglobals,  only : maxsiz,maxcrd,ires,r_min,r_max,sort_non_add,
        !            gamma_non_add,class_of_res_2

        use amhglobals,  only : SO,n_letters_con
        use globals_alt, only : altgamma,altpotflag,max_well_alt,onebody_gamma,&
               onebodyflag

        use altpot_interfaces, only: longscale

        implicit none 

        integer i,j,iwell,altpot_index,did_read_gamma(max_well_alt)
        integer read_status,did_read_onebody_gamma,onebody_index
            double precision x1,x2,x3
        character  str200*200,gammafil*200


        did_read_gamma=0
        did_read_onebody_gamma=0

        gammafil='gamma.dat'
        open(unit=130,file=gammafil,iostat=read_status,status='old')
!        open(unit=131,file='alt_gammacheck',iostat=read_status,status='unknown')
!        open(unit=132,file='onebody_gammacheck',iostat=read_status,status='unknown')

        !  Loop and search for keyword altpot in gamma.dat 
        !                                              !!!!!!!!!!!!!!!!!!! ! START READING LOOP
100     read (130,'(a)',err=99,end=999)str200                
        altpot_index=index(str200,'altpot')
        onebody_index=index(str200,'onebody')
        if(altpot_index.gt.0)then                   ! ALTPOT
           backspace(130)
           do i=1,n_letters_con,1
              do j=i,n_letters_con,1
                read (130,*)x1,x2,iwell
                altgamma(i,j,1,iwell)=x1
                altgamma(i,j,2,iwell)=x2
                altgamma(j,i,1,iwell)=altgamma(i,j,1,iwell)
                altgamma(j,i,2,iwell)=altgamma(i,j,2,iwell)
!                write(131,*)altgamma(i,j,1,iwell),altgamma(i,j,2,iwell), &
!                iwell
              enddo
           enddo
           did_read_gamma(iwell)=1
        endif

        if(onebody_index.gt.0)then                   ! ONEBODY
           backspace(130)
           do i=1,n_letters_con,1
              read (130,*)x1,x2,x3
              onebody_gamma(i,1)=x1
              onebody_gamma(i,2)=x2
              onebody_gamma(i,3)=x3
!              write(132,*)onebody_gamma(i,1),onebody_gamma(i,2), & 
!                          onebody_gamma(i,3),i
           enddo
           did_read_onebody_gamma=1
        endif

        goto 100                                            
        !                                            !!!!!!!!!!!!!!!!!! ! END READING LOOP

999     close (130)
!        close (131)
!        close (132)

        do i=1,max_well_alt    ! ERRORCHECK INPUT
           if((altpotflag(i).gt.0).and.(did_read_gamma(i).eq.0))then
              write(*,*)'ALTPOTFLAG  ',i,' SET BUT NOT PRESENT IN gamma.dat'
              stop
           endif
        enddo

        if((onebodyflag.gt.0) .and. (did_read_onebody_gamma.eq.0))then
         write(SO,*)'ONEBODYFLAG SET BUT NOT PRESENT IN gamma.dat will use HP_scale as ONEBODYGAMMAS'
           do i=1,n_letters_con,1
              write(132,*)onebody_gamma(i,1),onebody_gamma(i,2),onebody_gamma(i,3),i
           enddo
        endif

        call longscale()          ! The large N fudgefactor

        !                               ! ERROR 

99      if(read_status.ne.0)then   
           write(SO,*)'something wrong with input alternative potential'
           stop
        endif

        return

      end subroutine read_altgamma

      !--------------------------------------------- 
      
        subroutine finalize_alt()

        use globals_alt, only : outfile1_alt,timings_file_alt

        implicit none

!        close (outfile1_alt(1))
!        close (outfile1_alt(2))
!        close (outfile1_alt(3))
!        close (outfile1_alt(4))
!        close (outfile1_alt(5))
!        close (outfile1_alt(6))
!        close (outfile1_alt(7))
!        close (timings_file_alt)

        return
      end subroutine finalize_alt


      !....


      subroutine send_output_alt(A,E_alt,nmres,E_OB)

        use amhglobals,  only: maxsiz
        use globals_alt, only : output_step_alt, output_stepsize_alt,outfile1_alt,do_send_output &
             ,count_alt,T_alt,max_well_alt,treshold_alt,OB_density,OB_dns_count &
             , max_letters,aminoacids,timings_file_alt,accumulated_time
        implicit none
        double precision, intent(in) :: A(maxsiz),E_alt(2,max_well_alt),E_OB(3)
        integer, intent(in) ::  nmres
        integer n_in,i
   
       double precision rnmres,E_alt_sum


        output_step_alt=output_step_alt+1
        n_in=0
        rnmres=dble(nmres)

        !      str200="# step T  n_in/nmres  E_alt[ 11   12    21   22] end"
        !      if(output_step_alt.eq.1)write(outfile1_alt(1),'(a)')str200
        !      if(output_step_alt.eq.1)write(outfile1_alt(1),'(a)')str200(1:min(200,index('end',str200)))
!        if(output_step_alt.eq.1)then
!           open(unit=120,file='OUTPUT_GUIDE',status='unknown')
!           write(120,*)'Guide for check_alt_files'
!1           write(120,*)'check_alt_file_1: step   T  n_in/nmres  E_alt[ 11   21    12   22]  E_alt_tot'
!           !           write(120,*)'check_alt_file_2: step   T  n_in/nmres  A[i]  '
!           close(120)
!        endif

        if(do_send_output)then
           do i=1,nmres,1
              if(A(i).gt.treshold_alt)n_in=n_in+1
           enddo
           E_alt_sum=E_alt(1,1)+E_alt(1,2)+E_alt(2,1)+E_alt(2,2)
!           write(outfile1_alt(1),1000)count_alt,T_alt,dble(n_in)/rnmres,E_alt(1,1),E_alt(2,1) &
!                ,E_alt(1,2),E_alt(2,2),E_alt_sum
           !         rewind(outfile1_alt(2))
!           write(outfile1_alt(2),1000)count_alt,T_alt,(A(i),i=1,nmres)
!           write(timings_file_alt,800) count_alt,T_alt,(accumulated_time(i),i=1,4),(accumulated_time(i),i=10,13)
!800        format("Count T  Timings:[altpot][pp][ob][num] [altpot_util][pp1][pp2][pp3]", & 
!                                       i8,1x,f8.3,4x,4(f8.3,1x),2x,4(f8.3,1x))

        endif

        !      if((OB_dns_count.gt.1.0) .and. do_send_output)write(outfile1_alt(5),1000)count_alt,T_alt,OB_density/OB_dns_count
!        if((OB_dns_count.gt.1.0) .and. do_send_output)then
!           rewind(outfile1_alt(5))
!           do i=1,max_letters,1
!              write(outfile1_alt(5),*)1.0,OB_density(i,1),aminoacids(i)
!              write(outfile1_alt(5),*)2.0,OB_density(i,2),aminoacids(i)
!              write(outfile1_alt(5),*)3.0,OB_density(i,3),aminoacids(i)
!              write(outfile1_alt(5),2000)
!           enddo
!        endif
!        if(do_send_output)write(outfile1_alt(6),1000)count_alt,T_alt,E_OB


        if(mod(output_step_alt,output_stepsize_alt).ne.0)return

        !      write(outfile1_alt(1),60)

        !60    format(20(1x,e12.5))
1000    format(i8,2x,200(f8.3,2x))
2000    format()

        return

      end subroutine send_output_alt

      subroutine longscale()

        use amhglobals,  only:SO,  ab_c_of_n_old,ab_c_of_n_new,max_well, & 
              alpha_c_of_n,num_well,nmres,n_letters_con
        use globals_alt, only : altgamma

        implicit none
        integer i,j
            double precision long_nfactor(max_well),rnmres


        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !c   set long_nfactor which is some fudge factor for the
        !c   contact potentials (it partially allows for the effect of 
        !c   size -- which the fraction of contacts at a certain distance
        !c   depends on). It should be the same as included in the
        !c   optimisation. In corey's code it is in qchrgmk.f instead.
        !c   Johan: added long_nfactor for 2-wells 

        long_nfactor=1.0D0
        rnmres=dble(nmres)

        !c   start if alpha_c_of_n 
        if ( alpha_c_of_n ) then

           if (num_well .eq. 2) then
              long_nfactor(2)=1.0D0/(rnmres*0.012D0 + 0.87D0)
           elseif (num_well .eq. 3) then
              long_nfactor(2)=1.0D0/(rnmres*0.0065D0 + 0.87D0)
              long_nfactor(3)=1.0D0/(rnmres*0.04187261D0 + 0.1256658D0)
           elseif (num_well .eq. 10) then
              long_nfactor(1)=1.0D0/(rnmres*0.0008D0 + 0.09D0)
              long_nfactor(2)=1.0D0/(rnmres*0.0009D0 + 0.16D0)
              long_nfactor(3)=1.0D0/(rnmres*0.001D0 + 0.20D0)
              long_nfactor(4)=1.0D0/(rnmres*0.003D0 + 0.29D0)
              long_nfactor(5)=1.0D0/(rnmres*0.004D0 + 0.53D0)
              long_nfactor(6)=1.0D0/(rnmres*0.004D0 + 0.76D0)
              long_nfactor(7)=1.0D0/(rnmres*0.005D0 + 0.77D0)
              long_nfactor(8)=1.0D0/(rnmres*0.005D0 + 0.94D0)
              long_nfactor(9)=1.0D0/(rnmres*0.006D0 + 1.18D0)
              long_nfactor(10)=1.0D0/(rnmres*0.013D0 + 1.8D0)
           else
              write(SO,*) 'unsupported no of wells (alpha_c_of_n)',num_well
              stop
           endif

        endif
        !c   end if alpha_c_of_n 

        if (ab_c_of_n_new) then
           rnmres=dble(nmres)
           if (num_well .eq. 2) then
            long_nfactor(1)= ( 0.035D0*rnmres)/(rnmres*0.043D0 + 1.0D0)
            long_nfactor(2)=( 0.07D0*rnmres)/(rnmres*0.023D0 + 1.0D0)
            long_nfactor(1)=1.0D0/(long_nfactor(1))
            long_nfactor(2)=1.0D0/(long_nfactor(2))
           elseif (num_well .eq. 3) then
             long_nfactor(1)= ( 0.0843467D0*rnmres)/(rnmres*0.0453928D0 + 1.0D0)
             long_nfactor(2)=( 0.0669808D0*rnmres)/(rnmres*0.025112D0 + 1.0D0)
             long_nfactor(3)=( 0.18665D0*rnmres) /(rnmres*0.0107983D0 + 1.0D0)
             long_nfactor(1)=1.0D0/(long_nfactor(1))
             long_nfactor(2)=1.0D0/(long_nfactor(2))
             long_nfactor(3)=1.0D0/(long_nfactor(3))
           elseif (num_well .eq. 5) then

            long_nfactor(1)=( 0.0297375D0*rnmres) /(rnmres*0.02977935D0 + 1.0D0)
            long_nfactor(2)=( 0.0389704D0*rnmres) /(rnmres*0.021101D0 + 1.0D0)
            long_nfactor(3)=( 0.0596751D0*rnmres) /(rnmres*0.0133269D0  + 1.0D0)
            long_nfactor(4)=( 0.0681322D0*rnmres) /(rnmres*0.0100256D0 + 1.0D0)
            long_nfactor(5)=( 0.0729201D0*rnmres) /(rnmres*0.00347563D0 + 1.0D0)

              long_nfactor(1)=1.0D0/(long_nfactor(1))
              long_nfactor(2)=1.0D0/(long_nfactor(2))
              long_nfactor(3)=1.0D0/(long_nfactor(3))
              long_nfactor(4)=1.0D0/(long_nfactor(4))
              long_nfactor(5)=1.0D0/(long_nfactor(5))

           elseif (num_well .eq. 10) then

            long_nfactor(1)=( 0.0785047D0*rnmres) /(rnmres*0.245032D0 + 1.0D0)
            long_nfactor(2)=( 0.0152761D0*rnmres) /(rnmres*0.0283803D0 + 1.0D0)
            long_nfactor(3)=( 0.205481D0*rnmres) /(rnmres*0.385813D0  + 1.0D0)
            long_nfactor(4)=( 0.0174765D0*rnmres) /(rnmres*0.0174638D0 + 1.0D0)
            long_nfactor(5)=( 0.0352685D0*rnmres) /(rnmres*0.0269838D0 + 1.0D0)
            long_nfactor(6)=( 0.0474026D0*rnmres) /(rnmres*0.0249249D0 + 1.0D0)
            long_nfactor(7)=( 0.18665D0*rnmres) /(rnmres*0.0107983D0 + 1.0D0)
            long_nfactor(8)=( 0.0390303D0*rnmres) /(rnmres*0.0140943D0 + 1.0D0)
            long_nfactor(9)=( 0.0327411D0*rnmres) /(rnmres*0.00812347D0 + 1.0D0)
            long_nfactor(10)=( 0.0561461D0*rnmres) /(rnmres*0.00743991D0 + 1.0D0)

              long_nfactor(1)=1.0D0/(long_nfactor(1))
              long_nfactor(2)=1.0D0/(long_nfactor(2))
              long_nfactor(3)=1.0D0/(long_nfactor(3))
              long_nfactor(4)=1.0D0/(long_nfactor(4))
              long_nfactor(5)=1.0D0/(long_nfactor(5))
              long_nfactor(6)=1.0D0/(long_nfactor(6))
              long_nfactor(7)=1.0D0/(long_nfactor(7))
              long_nfactor(8)=1.0D0/(long_nfactor(8))
              long_nfactor(9)=1.0D0/(long_nfactor(9))
              long_nfactor(10)=1.0D0/(long_nfactor(10))


           else
              write(SO,*) 'num_well problem for ab_c_of_n_new',num_well
              stop 
           endif
        endif


        if (ab_c_of_n_old) then
           rnmres=dble(nmres)
           long_nfactor(1)=1.0D0/(rnmres*0.0015D0 + 1.94D0)
           long_nfactor(2)=1.0D0/(rnmres*0.0032D0 + 1.83D0)
           long_nfactor(3)=1.0D0/(rnmres*0.022D0 + 7.77D0)

           if (num_well.ne.3) then 
              write(SO,*) 'numwell must be 3 for ab_c_of_n_*old*',num_well
              stop
           endif

        endif

        do i=1,n_letters_con,1
           do j=i,n_letters_con,1
              altgamma(i,j,1,1)=altgamma(i,j,1,1)*long_nfactor(1)
              altgamma(i,j,2,1)=altgamma(i,j,2,1)*long_nfactor(1)
              altgamma(i,j,1,2)=altgamma(i,j,1,2)*long_nfactor(2)
              altgamma(i,j,2,2)=altgamma(i,j,2,2)*long_nfactor(2)
              altgamma(j,i,1,1)=altgamma(i,j,1,1)
              altgamma(j,i,2,1)=altgamma(i,j,2,1)
              altgamma(j,i,1,2)=altgamma(i,j,1,2)
              altgamma(j,i,2,2)=altgamma(i,j,2,2)
           enddo
        enddo

        return
      end subroutine longscale
