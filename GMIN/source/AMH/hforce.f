
c     --------------------- hforce.f ----------------------

      subroutine hforce(h_cord,nitcord,idx1,idx2,r1,r2,pot,factor,
     *  i_class,f_cord)
 
c     --------------------------------------------------

c     hdrgn finds the  potential due to hydrogen bonds between N and O     

c     ---------------------------------------------------

      use amhglobals,  only:maxsiz,maxcrd, prcord,ho_zero,NO_zero,sigma_NO,sigma_h,hbscl

      implicit none


c     argument declarations:

          
              double precision  h_cord(maxsiz,3),f_cord(maxsiz,3,maxcrd)

c     internal variables:

         integer idx1,idx2 
c        --- do loop indices ---

         integer i_axis,i_class 


         double precision  r1,factor,
     *         pot,
     *         nitcord(maxsiz,3),
     *         r2,
     *         dV_drNO,dV_drHO,drNO_dO(3),
     *         drNO_dN(3),drHO_dO(3),drHO_dH(3)

c     --------------------- begin -----------------------


c     find force due to hbonds

c do *NOT* zero f_cord here because the caller (hdrgn) is working out running total

        dV_drNO = -hbscl(i_class)*pot*((r1 - NO_zero)/
     *                                    (sigma_NO**2))*factor
        dV_drHO = -hbscl(i_class)*pot*((r2 - ho_zero)/
     *                                    (sigma_h**2))*factor

        do i_axis = 1,3
        
        drNO_dO(i_axis) = 
     *  (prcord(idx1,i_axis,1,3)-nitcord(idx2,i_axis)) /r1 

        drNO_dN(i_axis) = 
     *  -(prcord(idx1,i_axis,1,3)-nitcord(idx2,i_axis)) /r1 

        drHO_dO(i_axis) =
     *  (prcord(idx1,i_axis,1,3)-h_cord(idx2,i_axis))/r2

        drHO_dH(i_axis) = 
     *  -(prcord(idx1,i_axis,1,3)-h_cord(idx2,i_axis))/r2

c        vec1 is the force vector acting on N

c        vec1(i_axis) = -1.0*drNO_dN(i_axis)*dV_drNO
c        vec2(i_axis) = -1.0*dV_drHO*drHO_dH(i_axis)
c        vec3(i_axis)=
c    *  prcord(idx1,i_axis,1,3)-nitcord(idx2,i_axis)

        f_cord(idx2,i_axis,1) = f_cord(idx2,i_axis,1) -
     *  dV_drNO*drNO_dN(i_axis)*0.7032820        
        f_cord(idx2,i_axis,1) = f_cord(idx2,i_axis,1) -
     *  dV_drHO*drHO_dH(i_axis)*0.8929599

        f_cord(idx2-1,i_axis,1) = f_cord(idx2-1,i_axis,1) -     
     *  dV_drNO*drNO_dN(i_axis)*0.4831806   
        f_cord(idx2-1,i_axis,1) = f_cord(idx2-1,i_axis,1) -
     *  dV_drHO*drHO_dH(i_axis)*0.8409657

        f_cord(idx1,i_axis,3) = f_cord(idx1,i_axis,3) -     
     *  dV_drNO*drNO_dO(i_axis)   
        f_cord(idx1,i_axis,3) = f_cord(idx1,i_axis,3) -
     *  dV_drHO*drHO_dO(i_axis)

        f_cord(idx2-1,i_axis,3) = f_cord(idx2-1,i_axis,3) + 
     *  dV_drNO*drNO_dN(i_axis)*0.1864626   
        f_cord(idx2-1,i_axis,3) = f_cord(idx2-1,i_axis,3) +
     *  dV_drHO*drHO_dH(i_axis)*0.7338894

        enddo  ! i_axis        


c     ----------------------- done -----------------------

      return
      end
