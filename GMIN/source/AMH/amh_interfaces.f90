      module amh_interfaces


        interface
           subroutine additive_ev(distne,f_cord,nmres,E,numlng,ilong,tempav,crdixn, &
                xdiff,ydiff,zdiff,ccev_dist,pexcld)
             use amhglobals,  only:SO, maxsiz,maxtab,maxcnt,maxcrd
              double precision, intent(in):: distne(maxcnt,maxtab),xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab),&
                  zdiff(maxcnt,maxtab),ccev_dist(maxcnt,maxtab),pexcld
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
             integer, intent(in):: nmres,numlng(0:maxsiz,maxtab),ilong(maxcnt,2,maxtab),crdixn(maxtab,2)
             logical, intent(in):: tempav
           end subroutine additive_ev
        end interface


        interface
           subroutine  contact_P_AP(pro_cord,nmres,f_cord,lambda_P_AP)
             use amhglobals,  only : maxsiz,maxcrd
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
             integer, intent(in):: nmres
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: lambda_P_AP(3)
           end subroutine contact_P_AP
        end interface

        interface
           subroutine dssp_hdrgn_eastwood(pro_cord,f_cord,tempav,E)
             use amhglobals,  only:SO, maxsiz,maxcrd
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
             logical, intent(in)::  tempav
           end subroutine dssp_hdrgn_eastwood
        end interface

        interface
           subroutine E_write(avep,T,numpro,i_temp)
             implicit none
              double precision, intent(in):: avep(:,:,:),T
             integer, intent(in):: numpro,i_temp
           end subroutine E_write
        end interface

        interface
           subroutine ev_gamma(pro_cord,f_cord,E,tempav)
             use amhglobals,  only:SO, maxsiz,maxcrd
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
             logical, intent(in):: tempav
           end subroutine ev_gamma
        end interface

        interface
           subroutine gomb(  eta,index,tempav,ngomb,A_to_nmo,E,ibiasgauss,bias_av, &
                bias_var,bias_prefactor,bias_F,ibiaspoly,nbiaspoly,biaspoly)
             use amhglobals,  only:SO, maxcnt,maxsiz
              double precision, intent(in):: eta(1:maxcnt,1:4),ngomb,bias_av,bias_var,bias_prefactor,biaspoly(1:100)
              double precision, intent(out):: A_to_nmo(1:maxsiz,1:2),E(:,:),bias_F
             integer, intent(in):: index(1:maxcnt,1:4),nbiaspoly
             logical, intent(in):: tempav,ibiasgauss,ibiaspoly
           end subroutine gomb
        end interface


        interface
           subroutine hdrgn(pro_cord,f_cord,tempav,E)
             use amhglobals,  only:SO, maxsiz,maxcrd
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:)
             logical, intent(in)::  tempav
           end subroutine hdrgn
        end interface


        interface
           subroutine lookup(maxr,crdixn,vpotnt,forse,pro_cord,f_cord, &
                trgeng,numlng,nmres,rincinv,rincsq,       &
                igomb,ngomb,tempav,maxtab,ilong,          &
                E,ibiasgauss,bias_av,bias_var,            & 
                bias_prefactor,ibiaspoly,nbiaspoly,       &
                biaspoly,i_Qbias_a,i_Qbias_b,ccev_dist,pexcld)
             use amhglobals,  only:SO,maxcnt,maxsiz,maxcrd
             integer, intent(in):: maxr,maxtab 
             integer, intent(in):: crdixn(maxtab,2),numlng(0:maxsiz,maxtab),nmres, &
                  ilong(maxcnt,2,maxtab),nbiaspoly
              double precision, intent(in):: vpotnt(0:maxr+1,maxcnt,maxtab),  &
                  forse(0:maxr+1,maxcnt,maxtab),pro_cord(maxsiz,3,maxcrd),  &
                  rincinv(maxcnt,maxtab),rincsq(maxcnt,maxtab),ngomb,bias_av,bias_var,bias_prefactor, &
                  biaspoly(1:100),ccev_dist(maxcnt,maxtab),pexcld
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd),E(:,:),trgeng(maxtab,3)
             logical, intent(in):: igomb,tempav,ibiasgauss,ibiaspoly,i_Qbias_a,i_Qbias_b
           end subroutine lookup
        end interface

        interface
           subroutine  non_add_contact(pro_cord,nmres,tempav,f_cord,E)
             use amhglobals,  only : maxsiz,maxcrd
              double precision, intent(in):: pro_cord(maxsiz,3,maxcrd)
             integer, intent(in):: nmres
             logical, intent(in):: tempav
              double precision, intent(out):: f_cord(maxsiz,3,maxcrd)
              double precision, intent(out):: E
           end subroutine non_add_contact
        end interface

        interface
           subroutine oxy(maxtab,ires, jstrt,jfins,pro_cord,f_cord, &
                eqdist,oxscl,chrlscl,ramascl,             &  
                oxexcldv,numlng,ilong,nmres,E)
             use amhglobals,  only:SO, maxsiz,maxcrd,maxcnt
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd),eqdist(:),oxscl, &
                  chrlscl,ramascl
              double precision, intent(out) :: f_cord(maxsiz,3,maxcrd),E(:,:)
             integer, intent(in):: maxtab,ires(maxsiz),jstrt,jfins, &
                  numlng(0:maxsiz,maxtab),ilong(maxcnt,2,4),nmres
             logical, intent(in):: oxexcldv
           end subroutine oxy
        end interface


        interface
           subroutine Q_bias_seg_a(distne,f_cord,nmres,E,xdiff,ydiff,zdiff)
             use amhglobals,  only:SO, maxsiz,maxtab,maxcnt,maxcrd
              double precision, intent(in)::distne(maxcnt,maxtab), &
                       xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab), &
                       zdiff(maxcnt,maxtab)
              double precision, intent(out)::f_cord(maxsiz,3,maxcrd),E(:,:)
             integer, intent(in):: nmres
           end subroutine Q_bias_seg_a
        end interface

        interface
           subroutine Q_bias_seg_b(distne,f_cord,nmres,E,xdiff,ydiff,zdiff)
             use amhglobals,  only:SO, maxsiz,maxtab,maxcnt,maxcrd
              double precision, intent(in)::distne(maxcnt,maxtab), &
                       xdiff(maxcnt,maxtab),ydiff(maxcnt,maxtab), &
                       zdiff(maxcnt,maxtab)
              double precision, intent(out)::f_cord(maxsiz,3,maxcrd),E(:,:)
             integer, intent(in):: nmres
           end subroutine Q_bias_seg_b
        end interface

        interface
           subroutine rama(ires,jstrt,jfins,pro_cord,f_cord, &
                ramascl,ramapot,nitcord,cprcord)
             use amhglobals,  only:SO, maxsiz,maxcrd
              double precision, intent(in) :: pro_cord(maxsiz,3,maxcrd),ramascl, &
                  nitcord(maxsiz,3),cprcord(maxsiz,3)
              double precision, intent(out) :: f_cord(maxsiz,3,maxcrd),ramapot(maxsiz)
             integer, intent(in):: ires(maxsiz),jstrt,jfins
           end subroutine rama
        end interface

        interface
           subroutine rep_bias(prcord,frcord,E)
             use amhglobals,  only:SO, maxsiz,maxcrd,maxpro
              double precision, intent(in):: prcord(maxsiz,3,maxpro,maxcrd)
              double precision, intent(out):: frcord(maxsiz,3,maxpro,maxcrd),E
           end subroutine rep_bias
        end interface

        interface
           subroutine rep_contact(xcord)
             use amhglobals,  only:SO, maxsiz,maxcrd
              double precision, intent(in):: xcord(maxsiz,3,maxcrd)
           end subroutine rep_contact
        end interface


        interface
           subroutine Rg_bias(pro_cord,f_cord,E,tempav)
             use amhglobals,  only:SO, maxsiz,maxcrd
             logical, intent(in):: tempav
              double precision, intent(in), dimension(maxsiz,3,maxcrd)::  pro_cord
              double precision, intent(out), dimension(maxsiz,3,maxcrd):: f_cord
              double precision, intent(out):: E(:,:)
           end subroutine Rg_bias
        end interface


      end module amh_interfaces

