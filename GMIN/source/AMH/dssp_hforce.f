
c     --------------------- hforce.f ----------------------

       subroutine dssp_hforce(h_cord,nitcord,idx1,idx2,r1,r2,pot,
     * factor,lambda_hb,sigmaNO,sigmaH,pro_cord,f_cord)
 
c     --------------------------------------------------

c     hdrgn finds the  potential due to hydrogen bonds 
c     between N and O     

c     ---------------------------------------------------

      use amhglobals,  only: maxsiz,maxcrd,ho_zero,NO_zero

      implicit none

c     argument declarations:
          
       double precision  h_cord(maxsiz,3),lambda_hb,sigmaNO,
     *      sigmaH,f_cord(maxsiz,3,maxcrd),pro_cord(maxsiz,3,maxcrd)

c     internal variables:

         integer idx1,idx2 
c        --- do loop indices ---

         integer i_axis 


         double precision  r1,factor,pot,nitcord(maxsiz,3),r2,
     *         dV_drNO,dV_drHO,drNO_dO(3),drNO_dN(3),drHO_dO(3),drHO_dH(3)

c     --------------------- begin -----------------------

c WARNING: do not zero f_cord (the force) here as dssp_hdrgn is compiling running total


c     find force due to hbonds

        dV_drNO = -lambda_hb*pot*((r1 - NO_zero)/(sigmaNO**2))*factor
        dV_drHO = -lambda_hb*pot*((r2 - ho_zero)/(sigmaH**2))*factor

        do i_axis = 1,3
        
        drNO_dO(i_axis) = 
     *  (pro_cord(idx1,i_axis,3)-nitcord(idx2,i_axis)) /r1 

        drNO_dN(i_axis) = 
     *  -(pro_cord(idx1,i_axis,3)-nitcord(idx2,i_axis)) /r1 

        drHO_dO(i_axis) =
     *  (pro_cord(idx1,i_axis,3)-h_cord(idx2,i_axis))/r2

        drHO_dH(i_axis) = 
     *  -(pro_cord(idx1,i_axis,3)-h_cord(idx2,i_axis))/r2

c        apply force

        f_cord(idx2,i_axis,1) = f_cord(idx2,i_axis,1) -
     *  dV_drNO*drNO_dN(i_axis)*0.7032820D0
        f_cord(idx2,i_axis,1) = f_cord(idx2,i_axis,1) -
     *  dV_drHO*drHO_dH(i_axis)*0.8929599D0

        f_cord(idx2-1,i_axis,1) = f_cord(idx2-1,i_axis,1) -     
     *  dV_drNO*drNO_dN(i_axis)*0.4831806D0   
        f_cord(idx2-1,i_axis,1) = f_cord(idx2-1,i_axis,1) -
     *  dV_drHO*drHO_dH(i_axis)*0.8409657D0

        f_cord(idx1,i_axis,3) = f_cord(idx1,i_axis,3) -     
     *  dV_drNO*drNO_dO(i_axis)   
        f_cord(idx1,i_axis,3) = f_cord(idx1,i_axis,3) -
     *  dV_drHO*drHO_dO(i_axis)

        f_cord(idx2-1,i_axis,3) = f_cord(idx2-1,i_axis,3) + 
     *  dV_drNO*drNO_dN(i_axis)*0.1864626D0   
        f_cord(idx2-1,i_axis,3) = f_cord(idx2-1,i_axis,3) +
     *  dV_drHO*drHO_dH(i_axis)*0.7338894D0

        enddo  ! i_axis        

c     ----------------------- done -----------------------

      return
      end
