      subroutine rep_contact(xcord)

!     this subroutine will read in the experimental phivalues 
!     and will also provide the setup for the array that 
!     contains the information on contacts
!     it will also read in lambdas

!     declaring variables
      use amhglobals,  only: maxsiz,maxcrd,nmres,rep_phi_exp, &
             n_rep_con,rep_con_2_res,rep_lambda,ires, &
         rep_cut_off,i_rep_lambda_uniform,iphi,  &
         icontacts
       double precision, intent(in):: xcord(maxsiz,3,maxcrd)
      double precision :: r
      integer :: i,i_res,j_res,i_atom,j_atom

!     read in phi values and store them in rep_phi_exp(i)
      open(iphi,file='phivalues',status='old')
      do 10 i=1,nmres
      read(iphi,*) rep_phi_exp(i)
10    continue
      close(iphi)

!     read in lambda values and store them in lambda(i)
      if (i_rep_lambda_uniform) then
        rep_lambda(:)=rep_lambda(1)
      else
        open(2,file='lambdas',status='old')
        do 15 i=1,nmres
        read(2,*) rep_lambda(i)
15      continue
        close(2)
      endif


!     now we need to calculate the contacts made in the crystal structure 
!     and store them in the array rep_con_2_res
!     we use as cutoff distance rep_cut_off angstrom
!     we have the radius of the crystal structure


      n_rep_con=0

      do 20 i_res=1,nmres
      if (ires(i_res).eq.8) then   !note 'ires' is identity, and 'i_res' is residue counter
        i_atom=1     ! base distance on C_alpha for glycine
      else
        i_atom=2     ! base distance on C_beta otherwise
      endif
      do 30 j_res=1,i_res-2
      if (ires(j_res).eq.8) then
        j_atom=1     
      else
        j_atom=2     
      endif
      r=sqrt(sum((xcord(i_res,:,i_atom)-xcord(j_res,:,j_atom))**2))
      if ( r <= rep_cut_off) then
      n_rep_con(i_res)=n_rep_con(i_res)+1
      rep_con_2_res(i_res,n_rep_con(i_res))=j_res
      end if
30    continue
      do 40 j_res=i_res+2,nmres
      if (ires(j_res).eq.8) then
        j_atom=1     
      else
        j_atom=2     
      endif
      r=sqrt(sum((xcord(i_res,:,i_atom)-xcord(j_res,:,j_atom))**2))
      if ( r <= rep_cut_off) then
      n_rep_con(i_res)=n_rep_con(i_res)+1
      rep_con_2_res(i_res,n_rep_con(i_res))=j_res
      end if
40    continue
20    continue

      open(icontacts,file='contacts', status='new')
      do i=1,nmres
      write(icontacts,*) n_rep_con(i),rep_con_2_res(i,:)
      enddo
      close(icontacts)
     

      end
