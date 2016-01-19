module ion_pair

  implicit none

  integer,PARAMETER                        :: mxip=500,mxgrid=200,mxip_type=4 !
  
  TYPE rdf_t
     real(8),dimension(200)                :: rdf_sph,v_sph,mix_1,mix_2
     real(8)                               :: dens
  END TYPE rdf_t

  TYPE ip_potential_t
     real(8)                               :: rmin,rmax,dr,box_ref,rc
     real(8),DIMENSION(mxgrid,mxip_type,2) :: v_grid,dv_grid
     real(8),DIMENSION(mxip_type,2)        :: alfa,beta
  END TYPE ip_potential_t

  TYPE ip_info_t
     integer,DIMENSION(mxip,3)             :: index_ip
     integer                               :: n_ip
     integer,DIMENSION(mxip,2)             :: index_ca
  END TYPE ip_info_t


  type(ip_potential_t)                     :: potential
  type(rdf_t)                              :: rdf
  type(ip_info_t)                          :: ip_info


  real(8)                                  :: temp_ip=0.59227! here set a value for ambient temp.
  real(8),DIMENSION(3)                     :: energy_ip
  real(8)                                  :: ion_pair_scaling
  logical                                  :: ion_pair_control

  
contains

  subroutine ip_select(natom,text3,numres)


    integer                                          :: natom
    integer,DIMENSION(*)                             :: numres 
    character(len=5),DIMENSION(*)                    :: text3

    integer                                          :: ipp,inn,iat,i_cp,i_cn
    integer                                          :: ii,jj,kk,nn,np,npair,ires
    
    integer,dimension(mxip)                          :: index_ip_plus,index_ip_minus
    integer,dimension(mxip)                          :: itype_ip_plus,itype_ip_minus
    integer,dimension(mxip)                          :: resindex_ip_plus,resindex_ip_minus
    character(len=5)                                 :: name
    character(len=3),DIMENSION(2)                    :: crg_plus=(/ 'ARG','LYS' /),crg_minus=(/ 'ASP','GLU' /)

!    data crg_plus /'ARG','LYS'/
!    data crg_minus /'ASP','GLU'/


      ipp=0
      inn=0
      do iat=1,natom
         name=text3(iat)
         ires=numres(iat)
         do i_cp=1,2
            IF(name(3:5).eq.crg_plus(i_cp)) THEN 
               ipp=ipp+1
               index_ip_plus(ipp)=iat
               itype_ip_plus(ipp)=i_cp
               resindex_ip_plus(ipp)=ires
            ENDIF
         enddo
         do i_cn=1,2
            IF(name(3:5).eq.crg_minus(i_cn)) THEN  
               inn=inn+1
               index_ip_minus(inn)=iat
               itype_ip_minus(inn)=i_cn
               resindex_ip_minus(inn)=ires
            ENDIF
         enddo
      enddo
      

      np=ipp
      nn=inn
      npair=0
      
      !--- the 3 value define the type of potential according to the 2*2 matrix ION-PAIR(i_cp,i_cn)
      !--- (Asp/Arg,Asp/Lys)
      !--- (Glu/Arg,Glu/Lys)
      !--- index_ip(3,npair)=ii+jj*(jj-1) pick the type of potential to be read
 
      do i_cp=1,np
         do i_cn=1,nn
            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)+1)) cycle
            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)-1)) cycle
            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)+2)) cycle
            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)-2)) cycle
!            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)+3)) cycle
!            IF(resindex_ip_plus(i_cp).eq.(resindex_ip_minus(i_cn)-3)) cycle


            npair=npair+1
            ip_info%index_ip(npair,1)=index_ip_plus(i_cp)
            ip_info%index_ip(npair,2)=index_ip_minus(i_cn)
            ii=itype_ip_plus(i_cp)
            jj=itype_ip_minus(i_cn)
            ip_info%index_ip(npair,3)=ii+jj*(jj-1)

            ip_info%index_ca(npair,1)=ip_info%index_ip(npair,1)-1
            ip_info%index_ca(npair,2)=ip_info%index_ip(npair,2)-1
         enddo
      enddo
      ip_info%n_ip=npair


      print *,'-------------------------------------------'
      print *,'---> Tot Ion pairs :', ip_info%n_ip
      do kk=1,ip_info%n_ip
         print *,'---> Pairs N.',kk
         print *,'---> Elements:',text3(ip_info%index_ip(kk,1))&
              &,ip_info%index_ip(kk,1),text3(ip_info%index_ip(kk,2))&
              &,ip_info%index_ip(kk,2)
         print *,'---> Ion-Pair type:',ip_info%index_ip(kk,3)
      enddo

  end subroutine ip_select


  subroutine ip_potential_init
  


    integer,PARAMETER                      :: nf=27
    integer                                :: ip_type,ir,i,j
    real(8)                                :: rmin,rmax,dr,box,r
    character*10 type_str


!--- Executive statements

!--- Open I/O

      open(unit=nf,file='ip_potential.dat',status='unknown')

      read(nf,*) dr,rmin,rmax,box
      ir=INT((rmax-rmin)/dr)+1
      print *,'ir',ir
      do ip_type=1,mxip_type
         read(nf,*) type_str
        
         read(nf,*) (potential%alfa(ip_type,j),potential%beta(ip_type,j),j=1,2)
         do i=1,ir
            read(nf,*) r,(potential%v_grid(i,ip_type,j),potential%dv_grid(i,ip_type,j),j=1,2)            
            print *,r,ip_type,i,potential%v_grid(i,ip_type,1),potential%v_grid(i,ip_type,2)
         enddo
         read(nf,*)
      enddo
      close(nf)

    
    potential%rmin=rmin
    potential%rmax=rmax
    potential%dr=dr
    potential%box_ref=box
    potential%rc=rmax-2            ! to define a input fleg later

    print *,'-------------------------------------------'
    print *,'---> Ion-Pair cutoff:',potential%rc
    print *,'---> IP potential DONE'
    print *,'-------------------------------------------'



  end subroutine ip_potential_init

  subroutine rdf_init

!--- local variables

    integer                                :: imax
    real(8)                                :: rmax,dr
                               
!--- executive statement

    rmax=potential%rmax
    dr=potential%dr
    imax=INT(rmax/dr)+1
        
!    allocate(rdf%rdf_sph(imax))
!    allocate(rdf%v_sph(imax))

    rdf%rdf_sph=0.D0
    rdf%v_sph=0.D0
    rdf%mix_1=0.D0
    rdf%mix_2=0.D0

  end subroutine rdf_init

  subroutine rdf_out(nstep)

!--- local variables  
    integer                        :: imax,i,nstep,imin
    integer                        :: nrdf=38
    real(8)                        :: pi,v,r1,r2,dr,box
        
!--- executive statement

    pi=atan(1.0)*4.0
    imax=(potential%rmax/potential%dr)+1
    imin=(potential%rmin/potential%dr)+1
    dr=potential%dr
    box=potential%box_ref

    rdf%dens=1/(box**3)
!    v=pi*potential%rmax**3*(4/3.D0)
!    rdf%dens=1/v

!--- volume shells
    
    do i=imin,imax
       r1=(i-1)*dr
       r2=(i)*dr
       v=pi*(r2**3-r1**3)*(4/3.D0)
       rdf%v_sph(i)=v
    enddo
    
!--- print rdf

    open(unit=nrdf,file='rdf_ip.dat')
    do i=imin,imax
       write(nrdf,*) (i-1)*dr,rdf%rdf_sph(i)/(rdf%v_sph(i)*rdf%dens*DBLE(nstep))
    enddo
    close(nrdf)

!    open(unit=nrdf,file='rdf_scaling.dat')
!    do i=imin,imax
!       write(nrdf,*) (i-1)*dr,rdf%mix_1(i)/rdf%rdf_sph(i),rdf%mix_2(i)/rdf%rdf_sph(i)
!    enddo
!    close(nrdf)

  end subroutine rdf_out


!-------------------------------------------------------------------------------
!- subroutine ion_pair_force. SALT BRIDGE INTERACTIONS.
!- 10/5/2011. 
!- fabio.sterpone@ibpc.fr
!-------------------------------------------------------------------------------


  subroutine ion_pair_force(x,f,e_ip)
      
      implicit none

!----- Common variables

      common/PBC_R/periodicBC,CM
      logical periodicBC,CM

!----- main variables

      integer index_ip_l(mxip,3),n_ip,index_ca_l(mxip,2)
      integer nstate,istate
      real(8) e_ip,rc_ip
      real(8) x(*),f(*)
      

!----- local variables 

      integer i,I1,I2,IT,id
      integer ica1,ica2
      real(8) x1,x2,y1,y2,z1,z2,d,d2,x12,y12,z12
      real(8) fx,fy,fz,e,rc,dr
      real(8) pbc_mic
      real(8) alfa,beta,cteta,inv_mix,mix
      real(8) xca1,yca1,zca1,xca2,yca2,zca2
      real(8) nx1,ny1,nz1,nx2,ny2,nz2,nn1,nn2


      index_ip_l=ip_info%index_ip
      index_ca_l=ip_info%index_ca 
      n_ip=ip_info%n_ip
      rc=potential%rc
      dr=potential%dr

      e_ip=0.D0


      alfa=0
      beta=0
      cteta=0
      inv_mix=0
      mix=0

      do i=1,n_ip
         fx=0
         fy=0
         fz=0

         I1=index_ip_l(i,1)
         I2=index_ip_l(i,2)
         IT=index_ip_l(i,3)         
         x1=x(I1*3-2)
         x2=x(I2*3-2)
         y1=x(I1*3-1)
         y2=x(I2*3-1)
         z1=x(I1*3)
         z2=x(I2*3)
         x12=x2-x1
         y12=y2-y1
         z12=z2-z1
         if(periodicBC)then
            x12 = pbc_mic( x12 )
            y12 = pbc_mic( y12 )
            z12 = pbc_mic( z12 )
         endif

         d2=x12*x12+y12*y12+z12*z12
         d=dsqrt(d2)

!         id=INT(d/dr)
!         rdf%rdf_sph(id)=rdf%rdf_sph(id)+1


         IF(d.gt.rc) cycle   ! -- check this is correctly understood in .f
!         IF(d.gt.6.5) cycle   ! -- check this is correctly understood in .f

! here check if is type for two state (arg-xx) or is type for one state
! if arg-xx nstate=2 else nstate=1
! create a do loop for istate=1,nstate 
! call ip_potential (x12,y12,z12,d,it,fx,fy,fz,e,istate)
! ** change ip_potential v dv dimension
! *** add coefficient for mixing alfa1(istate) alfa2 (istate) to read in 
! *** check how case nstate=1 do not mix


         nstate=1
!         IF(d.lt.7) THEN
            IF(IT.eq.1) nstate=2    !--- asp-ARG
            IF(IT.eq.3) nstate=2    !--- glu-ARG 2
            IF(IT.eq.2) nstate=1
            IF(IT.eq.4) nstate=1
!         ENDIF


!--- cos(teta)  ca--->sc  sc<---ca


         IF(nstate.eq.2) THEN
            ICA1=index_ca_l(i,1)
            ICA2=index_ca_l(i,2)

            xca1=x(ICA1*3-2)
            yca1=x(ICA1*3-1)
            zca1=x(ICA1*3)
         
            xca2=x(ICA2*3-2)
            yca2=x(ICA2*3-1)
            zca2=x(ICA2*3)


            nx1=x1-xca1
            ny1=y1-yca1
            nz1=z1-zca1
            nn1=nx1*nx1+ny1*ny1+nz1*nz1

            nx2=x2-xca1
            ny2=y2-yca2
            nz2=z2-zca2
            nn2=nx2*nx2+ny2*ny2+nz2*nz2
         
            cteta=nx1*nx2+ny1*ny2+nz1*nz2
            cteta=cteta/dsqrt(nn1*nn2)
         ENDIF

!---


         do istate=1,nstate

            
            IF(nstate.eq.2) THEN
               IF(d.gt.10) THEN 
                  mix=0.5
               ELSE
                  alfa=potential%alfa(it,istate)
                  beta=potential%beta(it,istate)
                  
                  inv_mix=(exp(alfa*(cteta-beta))+1)
                  mix=1/inv_mix
               ENDIF
            ELSE
               mix=1
            ENDIF
!            IF(istate.eq.1) rdf%mix_1(id)=rdf%mix_1(id)+mix
!            IF(istate.eq.2) rdf%mix_2(id)=rdf%mix_2(id)+mix


            call ip_potential_grid(x12,y12,z12,d,it,fx,fy,fz,e,istate)
         
         

            F(I1*3-2) = F(I1*3-2)-mix*fx
            F(I1*3-1) = F(I1*3-1)-mix*fy
            F(I1*3)   = F(I1*3)-mix*fz
            F(I2*3-2) = F(I2*3-2)+mix*fx
            F(I2*3-1) = F(I2*3-1)+mix*fy
            F(I2*3)   = F(I2*3)+mix*fz
            e_ip=e_ip+mix*e
         

         enddo

      enddo


      return
    end subroutine ion_pair_force

!-------------------------------------------------------------------------------
!- subroutine ip_potential. SALT BRIDGE potential.
!- 10/5/2011. 
!- fabio.sterpone@ibpc.fr
!-------------------------------------------------------------------------------
      
    subroutine ip_potential_grid(x12,y12,z12,r,it,fx,fy,fz,e,i_state)
             
      integer it
      real(8) x12,y12,z12
      real(8) fx,fy,fz,e,r
      real(8) v_grid(200,4,2),dv_grid(200,4,2),dr_potential
      real(8) rmin,rmax
      real(8) temp_ip
      
!---- Function 
  
      real(8) interpolation

!---- Local Variables ------------------

      integer ir,i_state
      real(8) g,f,f1,f2
      real(8) factor,dr,ra,dva,dvb,va,vb

!---- test variable
      integer ite
      real(8) ftest,rtest


      dr=potential%dr
      rmin=potential%rmin

!---- factor scaling energy

      factor=ion_pair_scaling

!---- grad v(r) 

      ir=INT((r-rmin)/dr)+1
      ra=(ir-1)*dr+rmin

      dva=potential%dv_grid(ir,it,i_state)
      dvb=potential%dv_grid(ir+1,it,i_state)



      f1=interpolation(dr,ra,r,dva,dvb)


!---- f_x2=-grad v(r)*dr/dx2

      fx=-factor*f1*(x12/r)
      fy=-factor*f1*(y12/r)
      fz=-factor*f1*(z12/r)

!---- v(r)

      va=potential%v_grid(ir,it,i_state)
      vb=potential%v_grid(ir+1,it,i_state)
      e=interpolation(dr,ra,r,va,vb)

!      print *,'ir interpolation',ir,r,dr
!      print *,'derivative',dva,dvb
!      print *,'potential',f1
!      print *,'pot',va,vb
!      print *,'e',e,'f1',f1,x12/r,'ra',ra
!      stop
      e=e*factor

      
      return
    end subroutine ip_potential_grid


!-------------------------------------------------------------------------------
!- Subroutine interface. Define variables for excluding interctions in 
!- HYDROP and ENBOND calculations
!- 10/10/2011. 
!- fabio.sterpone@ibpc.fr
!-------------------------------------------------------------------------------


    subroutine ip_interface(index_ip,n_ip)

      integer,dimension(mxip,3)                 :: index_ip
      integer                                   :: n_ip


      index_ip=ip_info%index_ip
      n_ip=ip_info%n_ip


    end subroutine ip_interface
    
!-------------------------------------------------------------------------------
!- Subroutine write_ene_ip. SALT BRIDGE energy output. 
!- 10/10/2011. 
!- fabio.sterpone@ibpc.fr
!-------------------------------------------------------------------------------

    subroutine write_ene_ip(nstep)

      integer                               :: nstep
      real(8)                               :: e
      logical                               :: status

!--      energy_ip(1)   E_POT
!--      energy_ip(2)   E_IP
!--      energy_ip(3)   E_K

      e=energy_ip(1)+energy_ip(3)  !K+U
      INQUIRE(unit=80,opened=status)
      IF(STATUS) THEN
         write(80,*) nstep,energy_ip(1),energy_ip(2),energy_ip(3)
      ELSE
         open(unit=80,file='ene_ip.dat')
         write(80,*) nstep,energy_ip(1),energy_ip(2),energy_ip(3)
      ENDIF
      
      return
    end subroutine write_ene_ip




end module ion_pair

!-------------------------------------------------------------------------------
!- function interpolation. SALT BRIDGE potential. Linear Interpolation 
!- 10/10/2011. 
!- fabio.sterpone@ibpc.fr
!-------------------------------------------------------------------------------


  function interpolation(dx,xa,x,fxa,fxb)


    implicit none
  
    real(8),INTENT(IN)                     :: fxa,fxb,dx,xa,x
    real(8)                                :: a,b
    real(8)                                :: interpolation
      
    a=(fxb-fxa)/dx
    b=fxa-a*xa
    
    interpolation=a*x+b
      
  end function interpolation



