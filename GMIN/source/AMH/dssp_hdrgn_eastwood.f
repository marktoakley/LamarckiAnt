
c     --------------------- dssp_hdrgn.f ----------------------
       subroutine dssp_hdrgn_eastwood(pro_cord,f_cord,tempav,E)
c     --------------------------------------------------

c     hdrgn finds the  potential due to hydrogen bonds 
c     between N and O     

c     ---------------------------------------------------

      use amhglobals,  only: maxsiz,hbscl,nmres,prcord,ho_zero,NO_zero,
     *                  sigma_NO,sigma_h,para_HB,anti_HB,
     *                  anti_NHB,para_one,anti_one,maxcrd,
     *                  numseq_hb,tgsequences_hb,ave_seq_hb,ires
      implicit none

      double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
      double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
      logical, intent(in)::  tempav
                                                 
c     argument declarations:

      double precision theta(5), h_cord(maxsiz,3)

c     internal variables:

      integer isit1,isit2,i,hb_class,iseq,i_res,iter_tgseq
      logical i_repulsive,i_P,i_AP

c        --- do loop indices ---

      double precision  rNO(5), distances(maxsiz,maxsiz,2),
     *         hpot_tot,beta_scaling,nitcord(maxsiz,3),
     *         rHO(5),lambda(4),sigma(3,2),theta_seq_anti_HB(2),
     *         theta_seq_anti_NHB(2),theta_seq_para_HB(2),
     *         hval1,hval2,hval3,nval1,nval2,nval3

       external dssp_hforce

c     --------------------- begin -----------------------

c     zero force and energy
                                                                  
      f_cord=0.0D0
      E=0.0D0
      theta_seq_anti_NHB=0.0D0 
      theta_seq_para_HB=0.0D0
      theta_seq_anti_HB=0.0D0

      beta_scaling=1.0D0

      nval1=0.483D0
      nval2=0.703D0
      nval3=0.186D0

      hval1=0.8409657D0 
      hval2=0.8929599D0 
      hval3=0.7339256D0

c      write(6,*)'beta_scaling for hbonds =  ',beta_scaling 

c      do iseq = 1,numseq_hb
c         write(6,*)'seq in hbonds'
c         write(6,999)(tgsequences_hb(id3,iseq),id3=1,nmres)
c999      format(25(i2,1x))
c      enddo

c     calculate nitrogen and hydrogen coordinates
c     for residues 2 to N (we ignore terminal N)

        if (ave_seq_hb)then
             iter_tgseq = numseq_hb
        else
             iter_tgseq = 1
        endif

c          write(6,*)'dssp numseq_hb ', numseq_hb
c          write(6,*)'dssp ave_seq_hb ',ave_seq_hb
c          write(6,*)'iter_tgseq ', iter_tgseq

	  do  i_res=2,nmres

        nitcord(i_res,1)=nval1*prcord(i_res-1,1,1,1)
     *                   +nval2*prcord(i_res,1,1,1)
     *                   -nval3*prcord(i_res-1,1,1,3)

        nitcord(i_res,2)=nval1*prcord(i_res-1,2,1,1)
     *                   +nval2*prcord(i_res,2,1,1)
     *                   -nval3*prcord(i_res-1,2,1,3)

        nitcord(i_res,3)=nval1*prcord(i_res-1,3,1,1)
     *                   +nval2*prcord(i_res,3,1,1)
     *                   -nval3*prcord(i_res-1,3,1,3)


        h_cord(i_res,1) = hval1*prcord(i_res-1,1,1,1) 
     *                    +hval2*prcord(i_res,1,1,1) 
     *                    -hval3*prcord(i_res-1,1,1,3)
        h_cord(i_res,2) = hval1*prcord(i_res-1,2,1,1) 
     *                    +hval2*prcord(i_res,2,1,1) 
     *                    -hval3*prcord(i_res-1,2,1,3)
        h_cord(i_res,3) = hval1*prcord(i_res-1,3,1,1) 
     *                    +hval2*prcord(i_res,3,1,1) 
     *                    -hval3*prcord(i_res-1,3,1,3)

	enddo

! calculate r_ij (from O on i [isit1] to N and H on j [isit2])
! we ignore terminal oxygen (i position N) and terminal nitrogen (j position 1)

	do isit1 = 1,nmres-1
        do isit2 = 2,nmres

        distances(isit1,isit2,1) = dsqrt (
     *  (prcord(isit1,1,1,3)-nitcord(isit2,1))**2 +
     *  (prcord(isit1,2,1,3)-nitcord(isit2,2))**2 +
     *  (prcord(isit1,3,1,3)-nitcord(isit2,3))**2 )

        distances(isit1,isit2,2) = dsqrt(
     *  (prcord(isit1,1,1,3)-h_cord(isit2,1))**2 +
     *  (prcord(isit1,2,1,3)-h_cord(isit2,2))**2 +
     *  (prcord(isit1,3,1,3)-h_cord(isit2,3))**2 )

        enddo
        enddo

! reminder about the H-bond weights hbscl
! hbscl(1) short range additive
! hbscl(2) short range repulsive
! hbscl(3) med range additive
! hbscl(4) med range repulsive
! hbscl(5) med range anti-parallel non-add
! hbscl(6) med range parallel non-add
! hbscl(7) long range additive
! hbscl(8) long range repulsive
! hbscl(9) long range anti-parallel non-add
! hbscl(10) long range parallel non-add
!   seq-dep terms (that depend non-additively on identities of residues in H-bond pair)
! hbscl(11) med range seq-dep anti-P for H bonded pair
! hbscl(12) med range seq-dep anti-P for non-H bonded pair
! hbscl(13) long range seq-dep anti-P for H bonded pair
! hbscl(14) long range seq-dep anti-P for non-H bonded pair
! hbscl(15) long range seq-dep parallel (all cross strand pairs same for parallel)
!   seq-dep terms (that depend *additively* on identities of residues in H-bond pair)
! hbscl(16) seq-dep anti-P weight  
! hbscl(17) seq-dep parallel weight
 
        do  isit1=1,nmres-1
        do  isit2=2,nmres
c          reasons to increment loop
c          if (ires(isit2).eq.15) cycle !proline doen't have N-H

      lambda(1) = 0.0D0
	  lambda(2) = 0.0D0
	  lambda(3) = 0.0D0
	  lambda(4) = 0.0D0

           if (ires(isit2).eq.15) cycle !proline doen't have N-H
	   if (abs(isit2-isit1) .le. 2)  cycle
	   if (distances(isit1,isit2,1) .gt. 7.0D0) cycle

c discontunity attempt

	   if ((distances(isit1,isit2,1) .gt. 6.8D0) .and. 
     *	         (distances(isit1,isit2,1) .lt. 7.0D0)) beta_scaling=0.1D0

	   if  ((distances(isit1,isit2,1) .gt. 6.6D0) .and. 
     *	         (distances(isit1,isit2,1) .lt. 6.8D0)) beta_scaling=0.3D0

	   if ((distances(isit1,isit2,1) .gt. 6.4D0) .and. 
     *	         (distances(isit1,isit2,1) .lt. 6.6D0)) beta_scaling=0.5D0

	   if ((distances(isit1,isit2,1) .gt. 6.2D0) .and. 
     *	         (distances(isit1,isit2,1) .lt. 6.4D0)) beta_scaling=0.7D0

	   if ((distances(isit1,isit2,1) .gt. 6.0D0) .and. 
     *	         (distances(isit1,isit2,1) .lt. 6.2D0)) beta_scaling=0.9D0

!          terms involving a second H-bond (in addition to i->j)
!          that should be cut-off because they involve out-of-range atoms
!          or proline N-H

           i_repulsive=.true.
           i_AP=.true.
           i_P=.true.

      if (isit2.eq.nmres.or.ires(isit2+1).eq.15) i_repulsive=.false.

      if (isit1.eq.1 .or. isit2.eq.nmres .or. ires(isit1).eq.15) i_AP=.false.
      if (isit1.eq.nmres-1 .or. isit2.eq.nmres .or.ires(isit1+2).eq.15) i_P=.false.

!   although not currently being used, we have the ability to set
!   seperate gaussion parameters for alpha,anti-parallel beta and parallel beta
!      hydrogen bonds
          
            sigma(1,1) = sigma_NO
            sigma(1,2) = sigma_h
            sigma(2,1) = sigma_NO
            sigma(2,2) = sigma_h
            sigma(3,1) = sigma_NO
            sigma(3,2) = sigma_h

! first work out sequence-identity weights
! note an anti-P hbond (ie from i->j and j->i simultaneously)
! brings make i and j an anti parallel H-bonded pair
! and probably makes both j-1,i+1 and j+1,i-1 anti-par NON H-bonded pairs
!
! while a parallel H-bond from (ie from i->j and j->i+2 simultaneously)
! makes j and i+1 a parallel pair
 
       if (i_AP) then

      do iseq=1,iter_tgseq

       theta_seq_anti_HB(:)=  theta_seq_anti_HB(:) +  0.5D0*
     * anti_HB(tgsequences_hb(isit1,iseq),tgsequences_hb(isit2,iseq),:)

      theta_seq_anti_NHB(:)= theta_seq_anti_NHB(:) + 0.25D0*(
     * anti_NHB(tgsequences_hb(isit1+1,iseq),tgsequences_hb(isit2-1,iseq),:) +
     * anti_NHB(tgsequences_hb(isit1-1,iseq),tgsequences_hb(isit2+1,iseq),:))

       enddo

      if (ave_seq_hb) theta_seq_anti_HB(:) = theta_seq_anti_HB(:)/dble(iter_tgseq) 
      if (ave_seq_hb) theta_seq_anti_NHB(:) = theta_seq_anti_NHB(:)/dble(iter_tgseq) 

       endif !  if (i_AP) then

       if (i_P) then
         
       do iseq=1,iter_tgseq
        theta_seq_para_HB(:)= theta_seq_para_HB(:) + 
     *     para_HB(tgsequences_hb(isit1+1,iseq),tgsequences_hb(isit2,iseq),:)
       enddo

      if(ave_seq_hb)theta_seq_para_HB(:)=theta_seq_para_HB(:)/dble(iter_tgseq) 

          endif  !    if (i_P) then

! note that the 0.5 and 0.25 are conventions, but make sense, since
! the loop isit1 and isit2 (ie i and j) both goes over all sites, so
! will get the H_bonded pair twice and the non-H-bonded pair [which
! is a neighbour of *two* H-bonded pairs] four times. The parallel
! pairs are only picked up once.

      lambda(1) = -hbscl(7)
	  lambda(2) = -hbscl(8)

       do iseq=1,iter_tgseq
          lambda(3) =  lambda(3) + (-hbscl(9)
     *       -hbscl(13)*theta_seq_anti_HB(2)
     *       -hbscl(14)*theta_seq_anti_NHB(2)
     *       -hbscl(16)*(anti_one(tgsequences_hb(isit1,iseq))+ 
     *                      anti_one(tgsequences_hb(isit2,iseq))))

	  lambda(4) =  lambda(4) + (-hbscl(10) 
     *       -hbscl(15)*theta_seq_para_HB(2)
     *       -hbscl(17)*(para_one(tgsequences_hb(isit1+1,iseq))+
     *                      para_one(tgsequences_hb(isit2,iseq))))
       enddo
 
         if (ave_seq_hb) lambda(3) = lambda(3)/dble(iter_tgseq)
         if (ave_seq_hb) lambda(4) = lambda(4)/dble(iter_tgseq)

          hb_class=3
 
         if (abs(isit2-isit1) .lt. 20) then
           lambda(1) = -hbscl(3)
	       lambda(2) = -hbscl(4)

           do iseq=1, iter_tgseq 
             lambda(3) =  lambda(3) + (-hbscl(5)
     *         -hbscl(11)*theta_seq_anti_HB(2)
     *         -hbscl(12)*theta_seq_anti_NHB(2)
     *         -hbscl(16)*(anti_one(tgsequences_hb(isit1,iseq))+ 
     *                     anti_one(tgsequences_hb(isit2,iseq))))
           enddo

      if (ave_seq_hb) lambda(3) = lambda(3)/dble(iter_tgseq)

	      lambda(4) = -hbscl(6)
              hb_class=2

	  endif

          if (abs(isit2-isit1) .lt. 5) then
          lambda(1) = -hbscl(1)
	  lambda(2) = -hbscl(2)
	  lambda(3) = 0.0D0
	  lambda(4) = 0.0D0
          hb_class=1
	  endif

	rNO(1)=distances(isit1,isit2,1)
	rHO(1)=distances(isit1,isit2,2)

	rNO(2)=distances(isit1,isit2+1,1)
        rHO(2)=distances(isit1,isit2+1,2)

	rNO(4)=distances(isit2,isit1,1)
        rHO(4)=distances(isit2,isit1,2)

	rNO(5)=distances(isit2,isit1+2,1)
        rHO(5)=distances(isit2,isit1+2,2)

	  hpot_tot = 0.0D0

	   do i = 1,5   
             if (i.eq.3) cycle !note i=3 redundant (originally part of repulsive term)
      theta(i)=exp(-(rNO(i)-NO_zero)**2/(2.0D0*(sigma(hb_class,1))**2) -
     *             (rHO(i)-ho_zero)**2 /(2.0D0*(sigma(hb_class,2))**2) )

	   enddo

          lambda(1) = lambda(1) * beta_scaling 
          lambda(2) = lambda(2) * beta_scaling
          lambda(3) = lambda(3) * beta_scaling
          lambda(4) = lambda(4) * beta_scaling

! sometimes non-additive terms switched off 
!      (because involve out of bounds residues) as determined above

          if (.not.i_repulsive) theta(2)=0.0D0
          if (.not.i_AP) theta(4)=0.0D0
          if (.not.i_P) theta(5)=0.0D0
	   
	  hpot_tot= lambda(1)*theta(1) + 
     *    lambda(2)*theta(1)*theta(2) +
     *    lambda(3)*theta(1)*theta(4) +
     *    lambda(4)*theta(1)*theta(5)

	if (tempav) then
           E(1,1)=E(1,1)+hpot_tot  !total h-bond E
           E(1,13+hb_class)=E(1,13+hb_class)+hpot_tot !h-bond E by class

c$$$  now calculate frequencies for optimisation prog (these are -E/gamma)
c$$$  note that if gamma values in input file (ie hbscl) are +ve (which
c$$$  expect except for repulsive hbscl(2,6,10) then since all theta things 
c$$$  calculate here are positive the H-bond energies are mostly negative
c$$$  but the thing stored in E (-E/gamma) are positive

       E(1,21+(hb_class-1)*4)=E(1,21+(hb_class-1)*4)+theta(1) 
       E(1,22+(hb_class-1)*4)=E(1,22+(hb_class-1)*4)+theta(1)*theta(2)
       E(1,23+(hb_class-1)*4)=E(1,23+(hb_class-1)*4)+theta(1)*theta(4) 
       E(1,24+(hb_class-1)*4)=E(1,24+(hb_class-1)*4)+theta(1)*theta(5)

        if (hb_class.eq.2) then
           E(1,33)=E(1,33)+theta(1)*theta(4)*theta_seq_anti_HB(1)
           E(1,34)=E(1,34)+theta(1)*theta(4)*theta_seq_anti_NHB(1)
           elseif (hb_class.eq.3) then
           E(1,35)=E(1,35)+theta(1)*theta(5)*theta_seq_anti_HB(2)
           E(1,36)=E(1,36)+theta(1)*theta(5)*theta_seq_anti_NHB(2)
           E(1,37)=E(1,37)+theta(1)*theta(5)*theta_seq_para_HB(2)
           endif
 
        do iseq=1,iter_tgseq 
           E(1,38)=E(1,38)+ theta(1)*theta(4)*
     *          ( anti_one(tgsequences_hb(isit1,iseq)) + 
     *            anti_one(tgsequences_hb(isit2,iseq)))
       enddo

           if (hb_class.eq.3) then
             do iseq=1, iter_tgseq  
               E(1,39)=E(1,39)+ theta(1)*theta(5)*
     *              ( para_one(tgsequences_hb(isit1,iseq)) + 
     *                para_one(tgsequences_hb(isit2,iseq)))
             enddo
           endif

         endif !tempav

c     find force due to hbonds

        call dssp_hforce(h_cord,nitcord,isit1,isit2,
     *  rNO(1),rHO(1),theta(1),1.0D0,lambda(1),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)

        if (i_repulsive) then
        call dssp_hforce(h_cord,nitcord,isit1,isit2,
     *  rNO(1),rHO(1),theta(1),theta(2),lambda(2),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)


c    subroutine dssp_hforce(h_cord,nitcord,idx1,idx2,r1,r2,pot,
c     * factor,lambda_hb,sigmaNO,sigmaH,pro_cord,f_cord)


        call dssp_hforce(h_cord,nitcord,isit1,
     *  isit2+1,rNO(2),rHO(2),theta(2),theta(1),lambda(2),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)
        endif

        if (i_AP) then
        call dssp_hforce(h_cord,nitcord,isit1,isit2,
     *  rNO(1),rHO(1),theta(1),theta(4),lambda(3),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)

        call dssp_hforce(h_cord,nitcord,isit2,
     *  isit1,rNO(4),rHO(4),theta(4),theta(1),lambda(3),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)
        endif

        if (i_P) then
        call dssp_hforce(h_cord,nitcord,isit1,isit2,
     *  rNO(1),rHO(1),theta(1),theta(5),lambda(4),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)

        call dssp_hforce(h_cord,nitcord,isit2,
     *  isit1+2,rNO(5),rHO(5),theta(5),theta(1),lambda(4),
     *  sigma(hb_class,1),sigma(hb_class,2),pro_cord,f_cord)
        endif

        enddo   ! isit2	
	enddo  ! isit1

         if (ave_seq_hb) E(1,38)=E(1,38)/dble(iter_tgseq)
         if (ave_seq_hb) E(1,39)=E(1,39)/dble(iter_tgseq)

c     ----------------------- done -----------------------

      return
      end
