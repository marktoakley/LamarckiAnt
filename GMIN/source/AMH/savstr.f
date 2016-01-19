
c     --------------------- savstr ----------------------

      subroutine savstr(nmres,numpro,maxpro,
     *            maxcrd,prcord,ires,save_name,oconv)

c     --------------------------------------------------

c     SAVSTR save protein structures, eg. final.pdb

c     arguments:
c        nmres  - number of residues (i)
c        maxsiz - maximum number of residues (i)
c        numpro - number of trial protein structures (i)
c        maxpro - maximum number of trial protein 
c                 structures (i)
c        maxcrd - maximum number of atoms/ residue (i)
c        prcord - trial structures (i)
c        oconv  - unit id for file to which structures
c                 are to be written (i)

c     ---------------------------------------------------

      use amhglobals,  only: maxsiz


      implicit none

c     argument declarations:

       integer nmres,numpro,maxpro,
     *           maxcrd,oconv


        double precision prcord(maxsiz,3,maxpro,maxcrd)

        integer ires(maxsiz)

c     internal variables:

        double precision cprcord(maxsiz,3), nitcord(maxsiz,3)

        integer i_pro, i_res, i_axis, atom_no

        character*10 save_name
        integer nl

        character*3 res_type(maxsiz)

        external get_res_name

c     --------------------- begin -----------------------


        do 11 nl = 10, 1, -1

          if (save_name(nl:nl) .ne. ' ') then

             go to 12

          end if

11      continue
12      continue


         open(unit=oconv,file=save_name(1:nl)//'.pdb',
     *        status='new',form='formatted')


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calcultate N and C' positions

      do 490 i_pro=1,numpro
          do 491 i_res=1,nmres
              do 492 i_axis = 1,3
              
        cprcord(i_res,i_axis)=0.4436538*prcord(i_res,i_axis,i_pro,1)
     *                   +0.2352006*prcord(i_res+1,i_axis,i_pro,1)
     *                   +0.3211455*prcord(i_res,i_axis,i_pro,3)
        nitcord(i_res+1,i_axis)=0.4831806*prcord(i_res,i_axis,i_pro,1)
     *                   +0.7032820*prcord(i_res+1,i_axis,i_pro,1)
     *                   -0.1864626*prcord(i_res,i_axis,i_pro,3)


492           continue
491        continue
490     continue


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Write out final structure

        atom_no = 1

      do 502 i_pro=1,numpro
            do 500 i_res=1,nmres

            call get_res_name(ires(i_res), res_type(i_res))

cNNNNNNNNNNNNNNN
            if( i_res .gt. 1) then            ! write N position

         write(oconv,665) atom_no, res_type(i_res), i_res,
     *        (nitcord(i_res, i_axis), i_axis =1,3), atom_no

665        format('ATOM    ',i3,'  N   ', a3, '   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)

            atom_no = atom_no + 1


            end if

cCACACACACCACACACACA
         write(oconv,666) atom_no, res_type(i_res), i_res, 
     *        (prcord(i_res, i_axis, 1, 1), i_axis =1,3), atom_no
666        format('ATOM    ',i3,'  CA  ', a3, '   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)
            atom_no = atom_no + 1


cC'C'C'C'C'C'C'C'C'C'

            if( i_res .lt. nmres) then            ! write C' position

         write(oconv,667) atom_no, res_type(i_res), i_res,
     *        (cprcord(i_res, i_axis), i_axis =1,3), atom_no

667        format('ATOM    ',i3,'  C   ', a3, '   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)

            atom_no = atom_no + 1

cOOOOOOOOOOOOO                                  ! write O position

         write(oconv,668) atom_no, res_type(i_res), i_res,
     *        (prcord(i_res, i_axis, 1,3), i_axis =1,3), atom_no

668        format('ATOM    ',i3,'  O   ', a3, '   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)

            atom_no = atom_no + 1

            end if


             if (ires(i_res) .ne. 8) then
         write(oconv,669) atom_no, res_type(i_res), i_res,
     *        (prcord(i_res, i_axis, 1, 2), i_axis =1,3), atom_no
669        format('ATOM    ',i3,'  CB  ', a3, '   ',i3,'    ',3(f8.3),
     *             '  1.00  0.00      TPDB ',i3)
            atom_no = atom_no + 1

             end if 



  500       continue
  502 continue


      close(oconv)

c     ---------------------- done -----------------------

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_res_name(res_number, res_name)



        implicit none


        integer res_number

        character*3 res_name

        if (res_number .eq. 1) then 
                res_name =  "ALA" 
        endif
        if (res_number .eq. 2) then 
                res_name =  "ARG" 
        endif
        if (res_number .eq. 3) then 
                res_name =  "ASN" 
        endif
        if (res_number .eq. 4) then 
                res_name =  "ASP" 
        endif
        if (res_number .eq. 5) then 
                res_name =  "CYS" 
        endif
        if (res_number .eq. 6) then 
                res_name =  "GLN" 
        endif
        if (res_number .eq. 7) then 
                res_name =  "GLU" 
        endif
        if (res_number .eq. 8) then 
                res_name =  "GLY" 
        endif
        if (res_number .eq. 9) then 
                res_name =  "HIS" 
        endif
        if (res_number .eq. 10) then 
                res_name =  "ILE" 
        endif
        if (res_number .eq. 11) then 
                res_name =  "LEU" 
        endif
        if (res_number .eq. 12) then 
                res_name =  "LYS" 
        endif
        if (res_number .eq. 13) then 
                res_name =  "MET" 
        endif
        if (res_number .eq. 14) then 
                res_name =  "PHE" 
        endif
        if (res_number .eq. 15) then 
                res_name =  "PRO" 
        endif
        if (res_number .eq. 16) then 
                res_name =  "SER" 
        endif
        if (res_number .eq. 17) then 
                res_name =  "THR" 
        endif
        if (res_number .eq. 18) then 
                res_name =  "TRP" 
        endif
        if (res_number .eq. 19) then 
                res_name =  "TYR" 
        endif
        if (res_number .eq. 20) then 
                res_name =  "VAL" 
        endif


        if ((res_number .gt. 20) .or. (res_number .lt. 1)) then
           write(*,*) 'Residue out of Range'
        end if
             

        return

        end   

