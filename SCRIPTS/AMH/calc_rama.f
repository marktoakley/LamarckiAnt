c     ------------------- rama ----------------------
      subroutine calc_rama(Nres,cords,phi,psi)
c     -------------------------------------------------


c       number of residues, coordinates , return phi psi


      implicit none

c	a	vector from CA (prcord) to C'
c	b	vector from N to CA
c	c	vector from C' to N
c
c	adb	a dot b,  
c       bdc	b dot c
c	adc	a dot c
c
c	bxa	b cross a
c
c      sign : Returns the magnitude of its first argument with the sign of its
c            second argument.

c     argument declarations:
          
          integer maxres
          parameter (maxres=8000)
 
         integer i_res,i1,Nres,ires

         real nitcord(maxres,3),cprcord(maxres,3),
     *        a(3),b(3),c(3),pa(3),pb(3),pc(3),
     *      pa2,pb2,pc2,padb,padc,pbdc,pbxa(3),
     *      pcdbxa,pchi,pquad,a2,b2,c2,adb,adc,bdc,bxa(3),cdbxa,
     *        chi,quad,phi(maxres),psi(maxres),pkappa1,pkappa2,kappa1,
     *        kappa2,kappa,pkappa,cords(maxres,3,3)

c     --------------------- begin -----------------------

        do 600 i1=1,maxres
          phi(i1)=0.0
          psi(i1)=0.0
600     continue
c------------------------------------------------------------------
c     cords= atom number,coordinate(x,y,z),atom(1 = CA,2 = CB,3 = O)) 
c     calccoords= atom number,coordinate(x,y,z),atom(2 = N, 1 = C'))  
c------------------------------------------------------------------
cc     calculate nitrogen, cprime position

       do 500 i_res = ires, Nres 
          cprcord(i_res,1)=0.4436538*cords(i_res,1,1)
     *                   +0.2352006*cords(i_res+1,1,1)
     *                   +0.3211455*cords(i_res,1,3)
          nitcord(i_res+1,1)=0.4831806*cords(i_res,1,1)
     *                   +0.7032820*cords(i_res+1,1,1)
     *                   -0.1864626*cords(i_res,1,3)
          cprcord(i_res,2)=0.4436538*cords(i_res,2,1)
     *                   +0.2352006*cords(i_res+1,2,1)
     *                   +0.3211455*cords(i_res,2,3)
          nitcord(i_res+1,2)=0.4831806*cords(i_res,2,1)
     *                   +0.7032820*cords(i_res+1,2,1)
     *                   -0.1864626*cords(i_res,2,3)
          cprcord(i_res,3)=0.4436538*cords(i_res,3,1)
     *                   +0.2352006*cords(i_res+1,3,1)
     *                   +0.3211455*cords(i_res,3,3)
          nitcord(i_res+1,3)=0.4831806*cords(i_res,3,1)
     *                   +0.7032820*cords(i_res+1,3,1)
     *                   -0.1864626*cords(i_res,3,3)
         
500       continue
c if (startstruct-1.eq.0)then
c elseif (endstruct.eq.lastres)then
  
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHI
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 501 i_res=1+ires, Nres-1
c       calculate vectors	

            pa(1)=  cprcord(i_res,1)-cords(i_res,1,1)
            pa(2)=  cprcord(i_res,2)-cords(i_res,2,1)
            pa(3)=  cprcord(i_res,3)-cords(i_res,3,1)

             pb(1)=  cords(i_res,1,1)-nitcord(i_res,1)
            pb(2)=  cords(i_res,2,1)-nitcord(i_res,2)
            pb(3)=  cords(i_res,3,1)-nitcord(i_res,3)

            pc(1)=  nitcord(i_res,1)-cprcord(i_res-1,1)
            pc(2)=  nitcord(i_res,2)-cprcord(i_res-1,2)
            pc(3)=  nitcord(i_res,3)-cprcord(i_res-1,3)

c     calculate dot products 

          padb=pa(1)*pb(1)+pa(2)*pb(2)+pa(3)*pb(3)
          pbdc=pb(1)*pc(1)+pb(2)*pc(2)+pb(3)*pc(3)
           padc=pa(1)*pc(1)+pa(2)*pc(2)+pa(3)*pc(3)

         pbxa(1)=pb(2)*pa(3)-pb(3)*pa(2)
         pbxa(2)=pb(3)*pa(1)-pb(1)*pa(3)
        pbxa(3)=pb(1)*pa(2)-pb(2)*pa(1)

        pcdbxa=pc(1)*pbxa(1)+pc(2)*pbxa(2)+pc(3)*pbxa(3)

        pa2=pa(1)*pa(1)+pa(2)*pa(2)+pa(3)*pa(3)
        pb2=pb(1)*pb(1)+pb(2)*pb(2)+pb(3)*pb(3)
        pc2=pc(1)*pc(1)+pc(2)*pc(2)+pc(3)*pc(3)

        pkappa1=sqrt( pa2*pb2-padb*padb )
        pkappa2=sqrt( pb2*pc2-pbdc*pbdc )
         pkappa=pkappa1*pkappa2

          pchi=(padb*pbdc-padc*pb2)/pkappa
          pchi=max(min(pchi,0.9999999),-0.9999999)

          pquad=sign(1.0,pcdbxa)
           phi(i_res)=pquad*acos(pchi)


ccccccccccccccccccccccccccccccccccccccccccccccccc
c  PSI
ccccccccccccccccccccccccccccccccccccccccccccccccc

            a(1)=  nitcord(i_res+1,1)-cprcord(i_res,1)
            a(2)=  nitcord(i_res+1,2)-cprcord(i_res,2)
            a(3)=  nitcord(i_res+1,3)-cprcord(i_res,3)

            b(1)=  cprcord(i_res,1)-cords(i_res,1,1)
            b(2)=  cprcord(i_res,2)-cords(i_res,2,1)
            b(3)=  cprcord(i_res,3)-cords(i_res,3,1)

            c(1)=  cords(i_res,1,1)-nitcord(i_res,1)
            c(2)=  cords(i_res,2,1)-nitcord(i_res,2)
            c(3)=  cords(i_res,3,1)-nitcord(i_res,3)

c     calculate dot products

        adb=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
        bdc=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
        adc=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)

        bxa(1)=b(2)*a(3)-b(3)*a(2)
        bxa(2)=b(3)*a(1)-b(1)*a(3)
        bxa(3)=b(1)*a(2)-b(2)*a(1)

        cdbxa=c(1)*bxa(1)+c(2)*bxa(2)+c(3)*bxa(3)

        a2=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
        b2=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
        c2=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)

        kappa1=sqrt( a2*b2-adb*adb )
        kappa2=sqrt( b2*c2-bdc*bdc )
        kappa=kappa1*kappa2

        chi=(adb*bdc-adc*b2)/kappa
        chi=max(min(chi,0.9999999),-0.9999999)
        quad=sign(1.0,cdbxa)
        psi(i_res)=quad*acos(chi)

        if (psi(i_res).ge.3.14) psi(i_res)=-psi(i_res)
        if (phi(i_res).ge.3.14) phi(i_res)=-phi(i_res)

        phi(i_res) =  phi(i_res)*(180/3.14)            
        psi(i_res) =  psi(i_res)*(180/3.14)

501     continue
         
         return
      end
  
