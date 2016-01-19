      subroutine zero_amh

      use amhglobals

      implicit none

c     set hydrophobicity scale for each amino acid

      hydscl(0:21,1)=(/0.0D0,0.25D0,-1.8D0,-0.64D0,-0.72D0,
     *         0.04D0,-0.69D0,-0.62D0,0.16D0,-0.4D0,0.73D0,0.53D0,-1.1D0,
     *         0.26D0,0.61D0,-0.07D0,-0.26D0,-0.18D0,0.37D0,0.02D0,0.54D0,-1.1D0/)

c     random assignment of h scale

      hydscl(0:21,2)=(/0.0D0,1.0D0,1.0D0,-1.0D0,1.0D0,
     *             1.0D0,-1.0D0,-1.0D0,0.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0,
     *            -1.0D0,1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,1.0D0,-1.0D0,-1.0D0/)

       eqdist=(/2.45798D0,2.50665D0,2.44973D0,2.42677D0,2.82146D0/)

c     input_amh files

      icon=31
      imat=32
      ihomol=33
      ihbond=34
      iphi=35
      icontacts=37

      input_amh=41
      imemri=43
      iran=44
      iprolst=45
      iprolstscl=40 
      ievgamma=46
      ievbeta=47
      igamma=48
      iss_bias=49
      iss_struct=50
      itarg_seq=51
      imem_cons=52

c     output files

      oarchv=30
      oconv=76
      omovi=55
      ohdrgn=60
      ohdrgn_s=81
      ohdrgn_m=82
      ohdrgn_l=83
      ohdrgn_seq=85
      ononadd=84
      ocon_P_AP=86
      oobiassega=57
      oobiassegb=58

      orama=61
      ooxy=62
      ochiral=63
      oamh=64
      oamhsr=66
      oPE_no_bias=67
      oPE_with_bias=87
      oPE_plus_KE=88
      oPE_backbone=90
      oPE_backbone_norama=91

      oamhlr=68
      oamhmr=78
      oran=70
      occev=71
c     obias=72
      obias_Rg=77
      oKE=73
      ooev=74
      omoviseg=75
      orep=89
      ibiasfile=79

      ovgaspot=95
      ovscltab=96
c      pnames=97
      pdb=98
      SO=97
      oharmspring=120


      itgrd=0
      nmstep=0
      incmov=0

      t1=0.0D0
      t2=0.0D0
      t3=0.0D0
      t4=0.0D0
      t5=0.0D0

      nmtemp=0
      ictemp=.false.
      ctemp=0.0D0
      iresrt=.false.
      movanal=.false.

      temtur=0.0D0
      temgrd=0.0D0
      prcord=0.0D0
      bondln=0.0D0
      target=0.0D0
      zrcord=0.0D0
      avep=0.0D0
      chrlscl=0.0D0
      oxscl=0.0D0
      ramascl=0.0D0
      ccev_dist=0.0D0

      end
      
