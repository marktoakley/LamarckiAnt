
!     --------------------- hdrgn.f ----------------------

      subroutine hdrgn(pro_cord,f_cord,tempav,E)
 
!     --------------------------------------------------

!     hdrgn finds the  potential due to hydrogen bonds 
!     between N and O     

!     ---------------------------------------------------

      use amhglobals,  only: maxsiz,hbscl,nmres,ho_zero,NO_zero,ires,&
                       sigma_NO,sigma_h,hbond,maxcrd

      implicit none

!     argument declarations:

       double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
       double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
       logical, intent(in):: tempav
          
!     internal variables:

       integer isit1,isit2,add,i_class,i_res
       double precision distances(maxsiz,maxsiz,2),rNO(9),hpot(9), &
                        hpot_tot,nitcord(maxsiz,3),rHO(9),h_cord(maxsiz,3), &
                        hval1,hval2,hval3,nval1,nval2,nval3

!         parameter(sigma_h=1.732)  ! Thornton = 0.19  
!         parameter(sigma_NO=0.71) ! Thornton = 0.17
!         parameter(ho_zero=2.2)  ! Thornton = 2.06
!         parameter(NO_zero=2.960) ! Thornton = 2.98

         external hforce
!     --------------------- begin -----------------------

!     zero force and energy

      f_cord=0.0D0
      E=0.0D0


      nval1=0.483D0
      nval2=0.703D0
      nval3=0.186D0

      hval1=0.8409657D0
      hval2=0.8929599D0 
      hval3=0.7339256D0


!     calculate coordinates and distances

          do  i_res=2,nmres

        nitcord(i_res,1)=nval1*pro_cord(i_res-1,1,1) &
                        +nval2*pro_cord(i_res,1,1)  &
                        -nval3*pro_cord(i_res-1,1,3)

        nitcord(i_res,2)=nval1*pro_cord(i_res-1,2,1) &
                        +nval2*pro_cord(i_res,2,1)  &
                        -nval3*pro_cord(i_res-1,2,3)

        nitcord(i_res,3)=nval1*pro_cord(i_res-1,3,1) &
                        +nval2*pro_cord(i_res,3,1)  &
                        -nval3*pro_cord(i_res-1,3,3)


        h_cord(i_res,1) = hval1*pro_cord(i_res-1,1,1) &
                         +hval2*pro_cord(i_res,1,1)  &
                         -hval3*pro_cord(i_res-1,1,3)

        h_cord(i_res,2) = hval1*pro_cord(i_res-1,2,1) &
                         +hval2*pro_cord(i_res,2,1)  &
                         -hval3*pro_cord(i_res-1,2,3)

        h_cord(i_res,3) = hval1*pro_cord(i_res-1,3,1) &
                         +hval2*pro_cord(i_res,3,1)  &
                         -hval3*pro_cord(i_res-1,3,3)

   enddo


        do isit1 = 1,nmres
        do isit2 = 1,nmres

        distances(isit1,isit2,1) = dsqrt (           &
        (pro_cord(isit1,1,3)-nitcord(isit2,1))**2 + &
        (pro_cord(isit1,2,3)-nitcord(isit2,2))**2 + &
        (pro_cord(isit1,3,3)-nitcord(isit2,3))**2 )

        distances(isit1,isit2,2) = dsqrt(            &
        (pro_cord(isit1,1,3)-h_cord(isit2,1))**2 +  &
        (pro_cord(isit1,2,3)-h_cord(isit2,2))**2 +  &
        (pro_cord(isit1,3,3)-h_cord(isit2,3))**2 )

        enddo
        enddo


        do  isit1=3,nmres-2

        do  isit2 = 3,nmres-2

        if ((ires(isit1) .ne. 15) .and. (ires(isit2) .ne. 15) ) then
        if (abs(isit2-isit1) .gt. 2)  then
        if (distances(isit1,isit2,1) .lt. 7.0D0) then

        rNO(1)=distances(isit1,isit2,1)
        rHO(1)=distances(isit1,isit2,2)

        rNO(2)=distances(isit1-1,isit2,1)
        rHO(2)=distances(isit1-1,isit2,2)

        rNO(3)=distances(isit1+1,isit2,1)
        rHO(3)=distances(isit1+1,isit2,2)

        rNO(4)=distances(isit1-2,isit2,1)
        rHO(4)=distances(isit1-2,isit2,2)

        rNO(5)=distances(isit1+2,isit2,1)
        rHO(5)=distances(isit1+2,isit2,2)

        rNO(6)=distances(isit1,isit2-1,1)
        rHO(6)=distances(isit1,isit2-1,2)

        rNO(7)=distances(isit1,isit2+1,1)
        rHO(7)=distances(isit1,isit2+1,2)

        rNO(8)=distances(isit1,isit2-2,1)
        rHO(8)=distances(isit1,isit2-2,2)

        rNO(9)=distances(isit1,isit2+2,1)
        rHO(9)=distances(isit1,isit2+2,2)


          i_class = 3
          if (abs(isit2-isit1) .lt. 13) i_class = 2
          if (abs(isit2-isit1) .lt. 5)  i_class = 1


          do add = 1,9

          hpot(add)= -1.0D0*( exp(-(rNO(add)-NO_zero)**2/(2.0D0*(sigma_NO**2)) - &
          (rHO(add) - ho_zero)**2/(2.0D0*(sigma_h**2)) ) )

          enddo 

        hpot_tot = hpot(1) 

        do add = 2,9

        hpot_tot = hpot_tot + hpot(1)*hpot(add)

        enddo

        hpot_tot = hbscl(i_class)*hpot_tot

        if (tempav)then
          E(1,1)=E(1,1)+hpot_tot
          E(1,13+i_class)=E(1,13+i_class)+hpot_tot
        endif


!     find force due to hbonds

         if (hbond) then


        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),1.0D0,i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(2),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1-1,isit2,rNO(2),rHO(2),hpot(2),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(3),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1+1,isit2,rNO(3),rHO(3),hpot(3),hpot(1),i_class,f_cord) 

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(4),i_class,f_cord) 
        call hforce(h_cord,nitcord,isit1-2,isit2,rNO(4),rHO(4),hpot(4),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1+2,isit2,rNO(1),rHO(1),hpot(1),hpot(5),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1,isit2-1,rNO(5),rHO(5),hpot(5),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(6),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1,isit2-1,rNO(6),rHO(6),hpot(6),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(7),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1,isit2+1,rNO(7),rHO(7),hpot(7),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(8),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1,isit2-2,rNO(8),rHO(8),hpot(8),hpot(1),i_class,f_cord)

        call hforce(h_cord,nitcord,isit1,isit2,rNO(1),rHO(1),hpot(1),hpot(9),i_class,f_cord)
        call hforce(h_cord,nitcord,isit1,isit2+2,rNO(9),rHO(9),hpot(9),hpot(1),i_class,f_cord)

        
        endif ! endif hbond

        endif ! proline (ires = 15)
        endif ! distances cut-off
        endif ! isit2-isit1 .gt. 2

        enddo   ! isit2        


        enddo  ! isit1


!     ----------------------- done -----------------------

      return
      end
