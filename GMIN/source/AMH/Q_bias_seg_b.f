      subroutine Q_bias_seg_b(distne,f_cord,nmres,E,
     *                  xdiff,ydiff,zdiff)

c     calculates the contribution to forces (ie zrcord) due to
c     a potential that is polynomial in Q (Q depends on Ca positions
c     only, so it only gives forces on these atoms, but a dimension
c     for table number has been included in many arrays to hopefully
c     allow generalization

c format 1  -  start with q structure and constain to q with poly_nomial
c format 2  -  start with q structure and expand quartic polynomial
c format 3  -  start with random structure and expand quartic polynomial
 
      use amhglobals,  only: maxsiz,maxtab,maxcnt,maxcrd,i_bias_native_b,
     *   del_r_b,Q_ij_b,dq_dr_ij_b,qbiaspoly_b,targ_dist,n_Qbias_b,
     *   i_Q_format_b,Q0_b,Q_weight_b,Q_clip_b,i_ixn_Qbias_b,
     *   numconst_b, seglist_b,foldstrt_min_b,ss_b,Qvalue_b,
     *   foldstrt_max_b,n_divs_max,ss_dist,ss_pattern_b
c     argument declarations

      implicit none

       double precision, intent(in)::distne(maxcnt,maxtab),
     * xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab),
     * zdiff(maxcnt,maxtab)

       double precision, intent(out)::f_cord(maxsiz,3,maxcrd),E(:,:)
      integer, intent(in):: nmres

c     internal variables

      integer i_ixn,i_diff,i_dist,i_cord,i_res,i_tab,
     *        isit1,isit2,k,i,j,f
      double precision Q,r_dist,factor,Q_force(maxsiz,3,maxtab),
     *     V,dV_dQ  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        f_cord=0.0D0
        E=0.0D0

c      initialize Q_force

       do i_res=1,nmres
         do i_cord=1,3
           do i_tab=1,4
             Q_force(i_res,i_cord,i_tab)=0.0D0
            enddo
         enddo
       enddo
 
c      calculate Q and 'unscaled' contribution to zrcord

        Q=0.0D0
       
       do 50 i = 1,  numconst_b 
c   this cycle command prevent going outside array boundries
            if (i > nmres )cycle 

             isit1=seglist_b(i)

c identify the secondary structure unit
   
      if (ss_b)then
       f = 0
      if (seglist_b(i) <= foldstrt_max_b(1) .and. 
     *                    seglist_b(i) > foldstrt_min_b(1)) f = 1
      if (seglist_b(i) <= foldstrt_max_b(2) .and.
     *                    seglist_b(i) > foldstrt_min_b(2)) f = 2
      if (seglist_b(i) <= foldstrt_max_b(3) .and.
     *                    seglist_b(i) > foldstrt_min_b(3)) f = 3
      if (seglist_b(i) <= foldstrt_max_b(4) .and.
     *                    seglist_b(i) > foldstrt_min_b(4)) f = 4
      if (seglist_b(i) <= foldstrt_max_b(5) .and.
     *                    seglist_b(i) > foldstrt_min_b(5)) f = 5
      if (seglist_b(i) <= foldstrt_max_b(6) .and.
     *                    seglist_b(i) > foldstrt_min_b(6)) f = 6
      if (seglist_b(i) <= foldstrt_max_b(7) .and.
     *                    seglist_b(i) > foldstrt_min_b(7)) f = 7
      if (seglist_b(i) <= foldstrt_max_b(8) .and.
     *                    seglist_b(i) > foldstrt_min_b(8)) f = 8
      if (seglist_b(i) <= foldstrt_max_b(9) .and.
     *                    seglist_b(i) > foldstrt_min_b(9)) f = 9
      if (seglist_b(i) <= foldstrt_max_b(10) .and.
     *                    seglist_b(i) > foldstrt_min_b(10)) f = 10
      if (seglist_b(i) <= foldstrt_max_b(11) .and.
     *                    seglist_b(i) > foldstrt_min_b(10)) f = 11
      if (seglist_b(i) <= foldstrt_max_b(12) .and.
     *                    seglist_b(i) > foldstrt_min_b(12)) f = 12
      if (seglist_b(i) <= foldstrt_max_b(13) .and.
     *                    seglist_b(i) > foldstrt_min_b(13)) f = 13
      if (seglist_b(i) <= foldstrt_max_b(14) .and.
     *                    seglist_b(i) > foldstrt_min_b(14)) f = 14
      if (seglist_b(i) <= foldstrt_max_b(15) .and.
     *                    seglist_b(i) > foldstrt_min_b(15)) f = 15
      if (seglist_b(i) <= foldstrt_max_b(16) .and.
     *                    seglist_b(i) > foldstrt_min_b(16)) f = 16
      if (seglist_b(i) <= foldstrt_max_b(17) .and.
     *                    seglist_b(i) > foldstrt_min_b(17)) f = 17
      if (seglist_b(i) <= foldstrt_max_b(18) .and.
     *                    seglist_b(i) > foldstrt_min_b(18)) f = 18
      if (seglist_b(i) <= foldstrt_max_b(19) .and.
     *                    seglist_b(i) > foldstrt_min_b(19)) f = 19
      if (seglist_b(i) <= foldstrt_max_b(20) .and.
     *                    seglist_b(i) > foldstrt_min_b(20)) f = 20
      if (seglist_b(i) <= foldstrt_max_b(21) .and.
     *                    seglist_b(i) > foldstrt_min_b(21)) f = 21
      if (seglist_b(i) <= foldstrt_max_b(22) .and.
     *                    seglist_b(i) > foldstrt_min_b(22)) f = 22
      if (seglist_b(i) <= foldstrt_max_b(23) .and.
     *                    seglist_b(i) > foldstrt_min_b(23)) f = 23
        endif  ! ss_b

      do 100 j = i+2, numconst_b

       if ((ss_b) .and. (seglist_b(j)  >=  foldstrt_max_b(f) .or. 
     *              seglist_b(j)  < foldstrt_min_b(f)) )  cycle

            isit2 = seglist_b(j)
            i_diff = isit2-isit1
            i_ixn = i_ixn_Qbias_b(isit1,isit2)

       if (ss_b) then        
        r_dist=distne(i_ixn,1)-
     *   ss_dist(seglist_b(i),seglist_b(j),ss_pattern_b(f))
      endif

       if (.not. ss_b .and.  i_bias_native_b) then
        r_dist=distne(i_ixn,1)-targ_dist(i_ixn,1)
c         write(6,*)'rdist distne targ', 
c     *          r_dist, distne(i_ixn,1),targ_dist(i_ixn,1)
       endif 

       if (.not. ss_b  .and. .not. i_bias_native_b ) then
        r_dist=distne(i_ixn,1)-
     *   ss_dist(seglist_b(i),seglist_b(j),3)
       endif

        i_dist=int(abs(r_dist)/del_r_b(i_diff))+1
        if (i_dist.ge.n_divs_max) goto 100

        Q=Q+Q_ij_b(i_dist-1,i_diff)
     *       +dq_dr_ij_b(i_dist,i_diff)*(abs(r_dist)!linearly interp to get Q
     *       -del_r_b(i_diff)*float(i_dist-1) )         

       factor=dq_dr_ij_b(i_dist,i_diff)/sign(distne(i_ixn,1),r_dist)

        Q_force(isit1,1,1)=Q_force(isit1,1,1)
     *                           +factor*xdiff(i_ixn,1)
        Q_force(isit2,1,1)=Q_force(isit2,1,1)
     *                           -factor*xdiff(i_ixn,1)
        Q_force(isit1,2,1)=Q_force(isit1,2,1)
     *                           +factor*ydiff(i_ixn,1)
        Q_force(isit2,2,1)=Q_force(isit2,2,1)
     *                           -factor*ydiff(i_ixn,1)
        Q_force(isit1,3,1)=Q_force(isit1,3,1)
     *                           +factor*zdiff(i_ixn,1)
        Q_force(isit2,3,1)=Q_force(isit2,3,1)
     *                           -factor*zdiff(i_ixn,1)

100     enddo    ! j
50    enddo    ! i

c     calc Q, and hence V(Q) and dV(Q)/dQ
c        write(6,*)'Q_weight_b ',Q_weight_b
c        write(6,*)'n_Qbias_b ',n_Qbias_b
c        write(6,*)'i_Q_format_b ',i_Q_format_b
c        write(6,*)'Q_clip_b ',Q_clip_b
c        write(6,*)'Q0_b ',Q0_b
c        write(6,*)'Q ',Q

      Qvalue_b=Q

      V=0.0D0
      dV_dQ=0.0D0
      if (i_Q_format_b.eq.1) then
        do i=1,n_Qbias_b
          V=V+qbiaspoly_b(i)*Q**i
          dV_dQ=dV_dQ+float(i)*qbiaspoly_b(i)*Q**(i-1)
        enddo
      elseif (i_Q_format_b.eq.2 .or.
     *     (i_Q_format_b.eq.3.and.abs(Q-Q0_b).le.Q_clip_b) ) then

        V= Q_weight_b*(Q)**4  + (-4)*Q_weight_b*Q0_b*(Q)**3
     * + 6*Q_weight_b*(Q0_b)**2*(Q)**2 + (-4)*Q_weight_b*Q0_b**3*(Q)**1

        dV_dQ=dble(n_Qbias_b)*Q_weight_b*(Q-Q0_b)**(n_Qbias_b-1)
      elseif (i_Q_format_b.eq.3.and.(Q-Q0_b.gt.Q_clip_b)) then
        dV_dQ=dble(n_Qbias_b)*Q_weight_b*(Q_clip_b)**(n_Qbias_b-1)
        V=Q_weight_b*(Q_clip_b)**n_Qbias_b + dV_dQ*(Q-Q0_b-Q_clip_b)
      elseif (i_Q_format_b.eq.3.and.(Q-Q0_b.lt.-Q_clip_b)) then
        dV_dQ=dble(n_Qbias_b)*Q_weight_b*(-Q_clip_b)**(n_Qbias_b-1)
        V=Q_weight_b*(-Q_clip_b)**n_Qbias_b - dV_dQ*(Q0_b-Q-Q_clip_b)
      else
       write(6,*) 'i_Q_format_b wrong Q_bias',i_Q_format_b
       stop
      endif

      E(1,18)=E(1,18)+V

c     calc 'properly scaled' contribution to zrcord

       do k =1, numconst_b
          i_res  = seglist_b(k)
        do i_cord=1,3
          f_cord(i_res,i_cord,1)=f_cord(i_res,i_cord,1)-
     +                               dV_dQ*Q_force(i_res,i_cord,1)
        enddo
      enddo

      return
      end
