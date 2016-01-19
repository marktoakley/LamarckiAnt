
c     --------------------- opnfil ----------------------

      subroutine opnfil(iprolst,iprolstscl,input_amh,oarchv,
     *          omovi,ohdrgn,ohdrgn_s,ohdrgn_m,ohdrgn_l,ohdrgn_seq,
     *          orama,ooxy,ochiral,oamh,oamhsr,orep,
     *          oPE_no_bias,oPE_with_bias,oPE_backbone,
     *          oPE_backbone_norama,oamhlr,oamhmr,ononadd,
     *          ocon_P_AP,occev,ooev,obias_Rg,omoviseg,
     *          oobiassega,oobiassegb,oharmspring)

      use amhglobals,  only:SO

c     ---------------------------------------------------

c     OPNFIL opens the required i/o files

c     input_amh:

c        iprolst  - protein list (input_amh) (i)
c        input_amh - input_amh parameters (input_amh) (i)
c        oarchv- output/diagnostic file (output) (i)
c        omovi - intermediate structures 
c                throughout run (output) (i)

c     ---------------------------------------------------

      implicit none

c     argument declarations:

       integer iprolst,iprolstscl,input_amh,oarchv,omovi,
     *           ohdrgn,orama,ooxy,ochiral,oamh,oamhsr,orep,
     *           oKE,oamhlr,occev,ooev,oamhmr,obias_Rg,omoviseg,ohdrgn_s,
     *           ohdrgn_m,ohdrgn_l,ohdrgn_seq,ononadd,ocon_P_AP,
     *           oPE_no_bias,oPE_with_bias,oPE_plus_KE,oPE_backbone,
     *           oPE_backbone_norama,oobiassega,oobiassegb,oharmspring
c     --------------------- begin -----------------------

c     open specified files

c     input_amh files:

c        target proteins and memories

       open(unit=iprolst,file='pro.list',status='unknown',form='formatted')
c      open(unit=iprolstscl,file='pro.list.scl',status='unknown',form='formatted')

c        input_amh parameters

       open(unit=input_amh,file='input_amh',status='unknown',form='formatted')

c     output files:

c        archive file

       open(unit=oarchv,file='archive_amh',status='unknown',form='formatted')

c        hydrogen bond potential plot file

       open(unit=ohdrgn,file='hdrgn.plot',status='unknown')
       open(unit=ohdrgn_s,file='hdrgn_s.plot',status='unknown')
       open(unit=ohdrgn_m,file='hdrgn_m.plot',status='unknown')
       open(unit=ohdrgn_l,file='hdrgn_l.plot',status='unknown')
       open(unit=ohdrgn_seq,file='hdrgn_seq.plot',status='unknown')

c        rama potential plot file

       open(unit=orama,file='rama.plot',status='unknown')
        
c        oxy potential plot file

       open(unit=ooxy,file='oxy.plot',status='unknown')

c        chiral potential plot file

       open(unit=ochiral,file='chiral.plot',status='unknown')

c        amh potential plot file

       open(unit=oamh,file='amh.plot',status='unknown')

       open(unit=oamhsr,file='amhsr.plot',status='unknown')

c   replica interaction energy 
c        open(unit=orep,file='rep_bias.plot',status='unknown')
 
c        total potential at every time step

        open(unit=oPE_no_bias,file='PE_no_bias.plot',status='unknown')
        open(unit=oPE_with_bias,file='PE_with_bias.plot',status='unknown')
c        open(unit=oPE_plus_KE,file='PE_and_KE.plot',status='unknown')
        open(unit=oPE_backbone,file='backbone.plot',status='unknown')
         open(unit=oPE_backbone_norama,file='backbone_norama.plot',status='unknown')


        open(unit=oamhlr,file='amhlr.plot',status='unknown')
        open(unit=oamhmr,file='amhmr.plot',status='unknown')

c       non-additive contact potential 
        open(unit=ononadd,file='nonadd.plot',status='unknown')

c       P-AP contact potential 
        open(unit=ocon_P_AP,file='con_P_AP.plot',status='unknown')

c       potential due to carbon (a/b) excluded volume

        open(unit=occev,file='ccev.plot',status='unknown')

c       potential due to oxy excluded volume

        open(unit=ooev,file='ooev.plot',status='unknown')

c       potential due to biasing potential

c               open(unit=obias,file='bias.plot',status='unknown')

c       potential due to Rg-dependent biasing potential

        open(unit=obias_Rg,file='Rg_bias.plot',status='unknown')

c       partial biasing potential
        open(unit=oobiassega,file='biasseg_a.plot',status='unknown')
        open(unit=oobiassegb,file='biasseg_b.plot',status='unknown')

        open(unit=oharmspring,file='harm_spring.plot',status='unknown')


c        movie configurations
 
           open(unit=omovi,file='movie',status='unknown',
     *        form='formatted')
 
c        final movie configuration
 
c           open(unit=omoviseg,file='movieseg_final',status='unknown',
c     *        form='formatted')

c    
c          write(6,*)'SO', SO 
           open(unit=SO,file='output_amh',status='unknown')

c     end i/o tape openings

c     ---------------------- done -----------------------

      return
      end
