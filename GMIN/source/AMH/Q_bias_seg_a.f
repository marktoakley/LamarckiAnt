      subroutine Q_bias_seg_a(distne,f_cord,nmres,E,
     *                  xdiff,ydiff,zdiff)

c     calculates the contribution to forces (ie zrcord) due to
c     a potential that is polynomial in Q (Q depends on Ca positions
c     only, so it only gives forces on these atoms, but a dimension
c     for table number has been included in many arrays to hopefully
c     allow generalization

c format 1  -  start with q structure and constain to q with poly_nomial
c format 2  -  start with q structure and expand quartic polynomial
c format 3  -  start with random structure and expand quartic polynomial
 
      use amhglobals,  only: maxsiz,maxtab,maxcnt,maxcrd,
     *   del_r_a,Q_ij_a,dq_dr_ij_a,qbiaspoly_a,targ_dist,n_Qbias_a,
     *   i_Q_format_a,Q0_a,Q_weight_a,Q_clip_a,i_ixn_Qbias_a,
     *   numconst_a, seglist_a,foldstrt_min_a,ss_a,Qvalue_a,
     *   foldstrt_max_a,n_divs_max,ss_dist,ss_pattern_a
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
c        write(6,*)'Q_bias_seg_a  numconst_a',numconst_a
 
c       do 50 i = 1,  numconst_a,5
c   this cycle command prevent going outside array boundries

       do 50 i = 1,  numconst_a
            if (i > nmres )cycle 

             isit1=seglist_a(i)

c identify the secondary structure unit
   
       if (ss_a)then
        f = 0
       if (seglist_a(i) <= foldstrt_max_a(1) .and. 
     *                    seglist_a(i) > foldstrt_min_a(1)) f = 1
       if (seglist_a(i) <= foldstrt_max_a(2) .and.
     *                    seglist_a(i) > foldstrt_min_a(2)) f = 2
       if (seglist_a(i) <= foldstrt_max_a(3) .and.
     *                    seglist_a(i) > foldstrt_min_a(3)) f = 3
       if (seglist_a(i) <= foldstrt_max_a(4) .and.
     *                    seglist_a(i) > foldstrt_min_a(4)) f = 4
       if (seglist_a(i) <= foldstrt_max_a(5) .and.
     *                    seglist_a(i) > foldstrt_min_a(5)) f = 5
       if (seglist_a(i) <= foldstrt_max_a(6) .and.
     *                    seglist_a(i) > foldstrt_min_a(6)) f = 6
       if (seglist_a(i) <= foldstrt_max_a(7) .and.
     *                    seglist_a(i) > foldstrt_min_a(7)) f = 7
       if (seglist_a(i) <= foldstrt_max_a(8) .and.
     *                    seglist_a(i) > foldstrt_min_a(8)) f = 8
       if (seglist_a(i) <= foldstrt_max_a(9) .and.
     *                    seglist_a(i) > foldstrt_min_a(9)) f = 9
       if (seglist_a(i) <= foldstrt_max_a(10) .and.
     *                    seglist_a(i) > foldstrt_min_a(10)) f = 10
      if (seglist_a(i) <= foldstrt_max_a(11) .and.
     *                    seglist_a(i) > foldstrt_min_a(10)) f = 11
       if (seglist_a(i) <= foldstrt_max_a(12) .and.
     *                    seglist_a(i) > foldstrt_min_a(12)) f = 12
       if (seglist_a(i) <= foldstrt_max_a(13) .and.
     *                    seglist_a(i) > foldstrt_min_a(13)) f = 13
       if (seglist_a(i) <= foldstrt_max_a(14) .and.
     *                    seglist_a(i) > foldstrt_min_a(13)) f = 13
       if (seglist_a(i) <= foldstrt_max_a(15) .and.
     *                    seglist_a(i) > foldstrt_min_a(15)) f = 15
       if (seglist_a(i) <= foldstrt_max_a(16) .and.
     *                    seglist_a(i) > foldstrt_min_a(16)) f = 16
       if (seglist_a(i) <= foldstrt_max_a(17) .and.
     *                    seglist_a(i) > foldstrt_min_a(17)) f = 17
       if (seglist_a(i) <= foldstrt_max_a(18) .and.
     *                    seglist_a(i) > foldstrt_min_a(18)) f = 18
       if (seglist_a(i) <= foldstrt_max_a(19) .and.
     *                    seglist_a(i) > foldstrt_min_a(19)) f = 19
       if (seglist_a(i) <= foldstrt_max_a(20) .and.
     *                    seglist_a(i) > foldstrt_min_a(20)) f = 20
       if (seglist_a(i) <= foldstrt_max_a(21) .and.
     *                    seglist_a(i) > foldstrt_min_a(21)) f = 21
       if (seglist_a(i) <= foldstrt_max_a(22) .and.
     *                    seglist_a(i) > foldstrt_min_a(22)) f = 22
       if (seglist_a(i) <= foldstrt_max_a(23) .and.
     *                    seglist_a(i) > foldstrt_min_a(23)) f = 23
        endif  ! ss_a

c      do 100 j = i+2, numconst_a,4

      do 100 j = i+2, numconst_a
c        write(6,*)'Q_bias_seg_a i, j  ', i, j 

       if ((ss_a) .and. (seglist_a(j)  >=  foldstrt_max_a(f) .or. 
     *              seglist_a(j)  < foldstrt_min_a(f)) )  cycle

            isit2 = seglist_a(j)
            i_diff = isit2-isit1
            i_ixn = i_ixn_Qbias_a(isit1,isit2)

       if (ss_a) then        
        r_dist=distne(i_ixn,1)-
     *   ss_dist(seglist_a(i),seglist_a(j),ss_pattern_a(f))
      endif

       if (.not. ss_a) then
          r_dist=distne(i_ixn,1)-targ_dist(i_ixn,1)
       endif 
 
        i_dist=int(abs(r_dist)/del_r_a(i_diff))+1
        if (i_dist.ge.n_divs_max) goto 100

        Q=Q+Q_ij_a(i_dist-1,i_diff)
     *       +dq_dr_ij_a(i_dist,i_diff)*(abs(r_dist)!linearly interp to get Q
     *       -del_r_a(i_diff)*float(i_dist-1) )         

       factor=dq_dr_ij_a(i_dist,i_diff)/sign(distne(i_ixn,1),r_dist)

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
c        write(6,*)'Q_weight_a ',Q_weight_a
c        write(6,*)'n_Qbias_a ',n_Qbias_a
c        write(6,*)'i_Q_format_a ',i_Q_format_a
c        write(6,*)'Q_clip_a ',Q_clip_a
c        write(6,*)'Q0_a ',Q0_a
c        write(6,*)'Q ',Q

      Qvalue_a=Q
      V=0.0D0
      dV_dQ=0.0D0
      if (i_Q_format_a.eq.1) then
        do i=1,n_Qbias_a
          V=V+qbiaspoly_a(i)*Q**i
          dV_dQ=dV_dQ+float(i)*qbiaspoly_a(i)*Q**(i-1)
        enddo
      elseif (i_Q_format_a.eq.2 .or.
     *     (i_Q_format_a.eq.3.and.abs(Q-Q0_a).le.Q_clip_a) ) then
        V= Q_weight_a*(Q)**4  + (-4)*Q_weight_a*Q0_a*(Q)**3
     * + 6*Q_weight_a*(Q0_a)**2*(Q)**2 + (-4)*Q_weight_a*Q0_a**3*(Q)**1
        dV_dQ=dble(n_Qbias_a)*Q_weight_a*(Q-Q0_a)**(n_Qbias_a-1)
      elseif (i_Q_format_a.eq.3.and.(Q-Q0_a.gt.Q_clip_a)) then
        dV_dQ=dble(n_Qbias_a)*Q_weight_a*(Q_clip_a)**(n_Qbias_a-1)
        V=Q_weight_a*(Q_clip_a)**n_Qbias_a + dV_dQ*(Q-Q0_a-Q_clip_a)
      elseif (i_Q_format_a.eq.3.and.(Q-Q0_a.lt.-Q_clip_a)) then
        dV_dQ=dble(n_Qbias_a)*Q_weight_a*(-Q_clip_a)**(n_Qbias_a-1)
        V=Q_weight_a*(-Q_clip_a)**n_Qbias_a - dV_dQ*(Q0_a-Q-Q_clip_a)
      else
       write(6,*) 'i_Q_format_a wrong Q_bias',i_Q_format_a
       stop
      endif

      E(1,10)=E(1,10)+V

c     calc 'properly scaled' contribution to zrcord

       do k =1, numconst_a
          i_res  = seglist_a(k)

        do i_cord=1,3
          f_cord(i_res,i_cord,1)=f_cord(i_res,i_cord,1)-
     +                               dV_dQ*Q_force(i_res,i_cord,1)

        enddo
      enddo

      return
      end
