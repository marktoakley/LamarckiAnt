      subroutine rep_bias(prcord,frcord,E)

!     this subroutine is called by force.f and will
!     calculate the force on  replicas introduced 
!     by the interreplica potential
!     this force is matched with the experimental 
!     phivalues read in earlier

!     declaring variables
      use amhglobals,  only: maxsiz,maxcrd,nmres,numpro,maxpro,rep_lambda,rep_phi_exp,&
        rep_cut_off,rep_tol,n_rep_con,rep_con_2_res

       double precision, intent(in):: prcord(maxsiz,3,maxpro,maxcrd)
       double precision, intent(out):: frcord(maxsiz,3,maxpro,maxcrd),E

      integer :: i_res,n_rep,j_res,i_con,i_rep
      double precision :: r(maxpro,maxsiz,maxsiz),r0,phi_sim(nmres),kernel(maxpro,maxsiz)

!     i_res,j_res,n_rep are integers do loop over residues or number of replicas
!     maxsiz= max size for array for holding residues
!     maxcrd = 3 x,y,z
!     nmres  number of res of protein
!     numpro number of replicas
!     maxpro  parameter maximum number of replicas possible
!     rep_lambda is the lambda coefficient determining the uncertainties in the force
!     rep_phi_exp are the experimental phi values read in in rep_contact
!     rep_cut_off and rep_tol determine r0 the radius used for determining a contact 
      r0=rep_cut_off+rep_tol
!     n_rep_con is the amount total number of contacts made for residue i
!     rep_con_2_res specifies the residue label of the made contact 
!     prcord are the coordinates of the protein
!     frcord is the resulting force
!     tempav is a flag that stores energies or not
!     E is the total replica energy
!     r is an radius array of residues in specified replicas
!     phi_sim are the calculated phivalues for a certain residue
!     kernel is explained later, it ia phi*lambda

!     we want Cb distances, so we set mxcrd array entry to 2

!     we need to calculate all the
!     distances between residues during each folding run
!     we get the distances from prcord


      do n_rep=1,numpro
      do i_res=1,nmres
      do j_res=1,nmres
              r(n_rep,i_res,j_res)=sqrt(sum((prcord(i_res,:,n_rep,2)-prcord(j_res,:,n_rep,2))**2))
      end do
      end do
      end do





!     calculating the force, we have a constant piece for 
!     the force, called kernel ,and a multiplicative term
!     let us calculate the kernel first 
!     we need phi values to calculate the kernel
 
!     setting kernel and phisim equl to 0

      kernel=0.0D0
      phi_sim=0.0D0

!     loops to calculate kernel and phisim

      do i_rep=1,numpro
      do i_res=1,nmres
      do i_con=1,n_rep_con(i_res)
                j_res=rep_con_2_res(i_res,i_con)
                phi_sim(i_res)=phi_sim(i_res)+((0.50D0/n_rep_con(i_res))*(1.0D0+tanh(5.0D0*(r0-r(i_rep,i_res,j_res)))))/numpro
      end do
      end do
      end do


      do i_rep=1,numpro
      do i_res=1,nmres
          kernel(i_rep,i_res)=rep_lambda(i_res)*(phi_sim(i_res)-rep_phi_exp(i_res))
      end do
      end do

    


!     setting the force equal to zero

      frcord=0.0D0

!     force term in specified direction
!     we simply need to set up arrays that hold all the forces in x,y,z 
!     direction for residue i_res

      do n_rep=1,numpro
      do i_res=1,nmres
      if (n_rep_con(i_res) .eq.0) then
          frcord(i_res,:,n_rep,2)=frcord(i_res,:,n_rep,2)+0.0D0
      else
      do i_con=1,n_rep_con(i_res)
          j_res=rep_con_2_res(i_res,i_con)
          frcord(i_res,:,n_rep,2)=frcord(i_res,:,n_rep,2)+(prcord(i_res,:,n_rep,2)-prcord(j_res,:,n_rep,2))&
          *((kernel(n_rep,i_res)+kernel(n_rep,j_res))/r(n_rep,i_res,j_res))&
          *5.0D0*(1.0D0-(tanh(5.0D0*(r0-r(n_rep,i_res,j_res))))**2)
      end do
      end if
    
      end do
      end do
 





!     we also need to calculate the potential energy E
      E=0.0D0

      do i_res=1,nmres
      do i_rep=1,numpro
        E=E+(kernel(i_rep,i_res)**2)/rep_lambda(i_res)
      end do
      end do


      end






