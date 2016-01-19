C-----------------------------------------------------------------------
      SUBROUTINE writepdb(name,natoms,X,header)


      implicit double precision (a-h,o-z)

      parameter (MAXPRE = 1500)                  !! maximum number of residus
      parameter (MAXNAT = MAXPRE*6)             !! maximum number of atoms
      parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)  !! max number of nonbonded-pairs 
      parameter (MAXXC = 3*MAXNAT)  !! maximum number of cart coord

      COMMON/MISC1/NATOM,NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,
     $             NPHIA,NNB,NTYPES,
     $             MBONA,MTHETA,MPHIA
      common/textt/text2,text3,text4
      common/nnumres/numres,Id_atom

      common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
      integer nfrag, lenfrag, ichain

      common/PBC_R/periodicBC,CM
      logical periodicBC, CM

      CHARACTER*100 HEADER
      character*7 text2(MAXNAT)
      character*5 text3(MAXNAT)
      character*7 text4(MAXNAT) !! 
      integer  numres(MAXNAT),Id_atom(MAXNAT)
      CHARACTER*100 name
      DIMENSION X(*)

      integer  atom_i3
      real(8) X_a(MAXXC)
     
      integer ns
!      logical UNITOK, UNITOP

cATOM      3  ASP ASP     1       1.718   5.627 -12.238
cATOM      4  C   ASP     1       1.364   3.436  -9.969
cATOM      5  O   ASP     1       1.311   2.290  -9.529

      integer idatom,i,j
      character*7 string7
      
      X_a(1:natom3) = X(1:natom3)

      if(periodicBC .and. .not. CM)then
        do atom_i3 = 1, natom3
            X_a(atom_i3) = pbc_mic( X_a(atom_i3) )
        end do
      endif

      ns=11
!      inquire (unit=ns,exist=UNITOK,opened=UNITOP)
!      if (UNITOK .and. .not. UNITOP) then
        open(unit=ns,file=name,status='unknown',position='append')
!      endif
      write(ns,'(a)') "MODEL"
      write(ns,'(a)') header


      if (periodicBC .and. CM) then
        call writeCM(natoms,X_a)
      end if

      idatom = 0
      do i=1,nfrag
       do j = 1, lenfrag(i)
        idatom = idatom + 1
        string7 = text4(idatom)
        string7(6:6) = char( ichain(idatom)+64 )

        write(ns,'(a,I4,a,a,I3,f12.3,2f8.3)') text2(idatom), id_atom(idatom), text3(idatom),
     $  string7, numres(idatom), 
     $  X_a(idatom*3-2), X_a(idatom*3-1), X_a(idatom*3)

       end do
      enddo
      
      write(ns,'(a)') "ENDMDL"
      close(ns)
      RETURN
      END

