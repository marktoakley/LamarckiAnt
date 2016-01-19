      subroutine ev_gamma(pro_cord,f_cord,E,tempav)

! excluded volume between gamma atoms

      use amhglobals,  only:SO, maxsiz,maxcrd,nmres,pexcld_gamma, & 
           exvmin_gamma, i_type_ev_gamma,c_of_m_dist,ires

      implicit none

       double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
       double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
      logical, intent(in):: tempav

      double precision :: gamma_cord(maxsiz,3),diff(3),dist,factor(3),xi
      double precision :: pexcld_mcp, exvmin_mcp,distmcp


      integer :: i,j

!     zero force and energy

      f_cord=0.0D0
      E=0.0D0
      
! calculate gamma positions

      do i=1,nmres
        gamma_cord(i,:)=pro_cord(i,:,1)+c_of_m_dist(i)*(pro_cord(i,:,2)-pro_cord(i,:,1))
      enddo
      
      if (i_type_ev_gamma.eq.1) then !quadratic (hard function)
        do i=1,nmres-13
        do j=i+13,nmres 
          diff(:)=gamma_cord(i,:)-gamma_cord(j,:)
          dist=dsqrt( dot_product( diff,diff ) )
          if (dist.gt.exvmin_gamma(i,j)) cycle
          factor(:)=-2.0D0*pexcld_gamma*( dist-exvmin_gamma(i,j) )*diff(:)/dist
          f_cord(i,:,1)=f_cord(i,:,1)+factor*(1.0D0-c_of_m_dist(i))
          f_cord(j,:,1)=f_cord(j,:,1)-factor*(1.0D0-c_of_m_dist(j))
          f_cord(i,:,2)=f_cord(i,:,2)+factor*c_of_m_dist(i)
          f_cord(j,:,2)=f_cord(j,:,2)-factor*c_of_m_dist(j)
          if (tempav) E(1,9)=E(1,9)+pexcld_gamma*( dist-exvmin_gamma(i,j) )**2
        enddo
        enddo
      elseif (i_type_ev_gamma.eq.2) then !tanh (soft function)
        xi=0.5D0
        do i=1,nmres-1
        do j=i+1,nmres
          diff(:)=gamma_cord(i,:)-gamma_cord(j,:)
          dist=dsqrt( dot_product( diff,diff ) )
          factor(:)=(pexcld_gamma/(2.0D0*xi))*( cosh((exvmin_gamma(i,j)-dist)/xi)**(-2) ) *diff(:)/dist
          f_cord(i,:,1)=f_cord(i,:,1)+factor*(1.0D0-c_of_m_dist(i))
          f_cord(j,:,1)=f_cord(j,:,1)-factor*(1.0D0-c_of_m_dist(j))
          f_cord(i,:,2)=f_cord(i,:,2)+factor*c_of_m_dist(i)
          f_cord(j,:,2)=f_cord(j,:,2)-factor*c_of_m_dist(j)
          if (tempav) E(1,9)=E(1,9)+0.5D0*pexcld_gamma*(1.0D0+tanh((exvmin_gamma(i,j)-dist)/xi))
        enddo
        enddo
      elseif (i_type_ev_gamma.eq.3) then !tanh (soft function MCP)


        xi=0.5D0
        pexcld_mcp=10.0D0
        exvmin_mcp=2.0D0

        do i=1,nmres-10
        do j=i+10,nmres

!          if( (ires(j).eq.8) .or. (ires(i).eq.8) ) cycle

          diff(1)=pro_cord(i,1,2)-pro_cord(j,1,2)
          diff(2)=pro_cord(i,2,2)-pro_cord(j,2,2)
          diff(3)=pro_cord(i,3,2)-pro_cord(j,3,2)
          distmcp=dsqrt( diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3) )


         factor(1)=(pexcld_mcp/(2.0D0*xi))*(cosh((exvmin_mcp-distmcp)/xi)**(-2))*diff(1)/distmcp
         factor(2)=(pexcld_mcp/(2.0D0*xi))*(cosh((exvmin_mcp-distmcp)/xi)**(-2))*diff(2)/distmcp
         factor(3)=(pexcld_mcp/(2.0D0*xi))*(cosh((exvmin_mcp-distmcp)/xi)**(-2))*diff(3)/distmcp

          f_cord(i,1,2)=f_cord(i,1,2)+factor(1)
          f_cord(j,1,2)=f_cord(j,1,2)-factor(1)

          f_cord(i,2,2)=f_cord(i,2,2)+factor(2)
          f_cord(j,2,2)=f_cord(j,2,2)-factor(2)

          f_cord(i,3,2)=f_cord(i,3,2)+factor(3)
          f_cord(j,3,2)=f_cord(j,3,2)-factor(3)

          E(1,9)=E(1,9)+0.5D0*pexcld_mcp*(1.0D0+tanh((exvmin_mcp-distmcp)/xi))
        enddo
        enddo

!        xi=0.5D0

!        do i=1,nmres
!          gamma_cord(i,:)=pro_cord(i,:,2)
!        enddo

!        do i=1,nmres-3
!        do j=i+3,nmres
!          diff(:)=gamma_cord(i,:)-gamma_cord(j,:)
!          dist=dsqrt( dot_product( diff,diff ) )
!          factor(:)=(pexcld_gamma/(2.0D0*xi))*( cosh((exvmin_gamma(i,j)-dist)/xi)**(-2) ) *diff(:)/dist
!          f_cord(i,:,1)=f_cord(i,:,1)+factor(:)
!          f_cord(j,:,1)=f_cord(j,:,1)-factor(:)
!          f_cord(i,:,2)=f_cord(i,:,2)+factor(:)
!          f_cord(j,:,2)=f_cord(j,:,2)-factor(:)
!          if (tempav) E(1,9)=E(1,9)+0.5D0*pexcld_gamma*(1.0D0+tanh((exvmin_gamma(i,j)-dist)/xi))
!        enddo
!        enddo

      else
        write(SO,*) 'i_type_ev_gamma cock-up',i_type_ev_gamma
        stop
      endif

      return
      end
