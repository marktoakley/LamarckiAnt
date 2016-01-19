module writetofile
  use geometric_corrections
  implicit none

contains

subroutine  write_to_file(etype,id,natoms, pos, posref, ndigits, aname, counter, value,  &
                          energy,save_unrotated,singlefile)
  use defs, only : usextc, SIMULATION_TYPE
  implicit none

  integer, intent(in)                      :: natoms, counter, ndigits, id
  logical, intent(in)                      :: save_unrotated, singlefile
  character(len=*), intent(in)             :: aname
  character(len=*), intent(in)             :: etype
  double precision, intent(in)                      :: value, energy
  double precision, dimension(3*natoms), intent(in) ::pos, posref

  integer                      :: npart
  double precision                      :: delr, rmsd
  double precision, dimension(3*natoms) :: posa, posb
  character(len=100)           :: fname, header
  character(len=10)  :: textstring 
  character(len=16)  :: textstring2
  character(len=30)  :: scounter, chaine
  character(len=4)  :: PDB_EXT   = '.pdb'

  logical, save :: openf = .true.
  character(5) :: mname
  integer, save :: xd, ret
  real :: box(9), prec
  !     from xdrfile doc:
  !     prec ~= 1/precision, i.e prec = 1000 => precision 1e-3
  prec = 1000
  box = 0

  ! In all cases we re-orient the conformation
  posa = pos
  chaine = ''

  if (.not. save_unrotated) then 
    posb = posref
    call minimize_rmsd(natoms,posb,posa,delr,rmsd,npart)
  endif

  ! We now save conformation
  if (etype .eq. 'relaxe') then
     fname = trim(aname)
     call convert_to_chain(id,chaine,"0")
     call real_to_chain(energy,16,6,textstring2)
     header = "HEADER - id: " // chaine(29:30) // " - configuration energy : "  &
     &     // textstring2
     call writepdb(fname,natoms,posa,header)
!     call mywritextc(fname,natoms,posa,counter,value)

  else if (etype .eq. 'thermalize') then
     call convert_to_chain(counter,scounter,"0")
     fname =   trim(aname) // scounter(31-ndigits:30) // PDB_EXT
     call convert_to_chain(id,chaine,"0")
     call real_to_chain(value,10,6,textstring)
     call real_to_chain(energy,16,6,textstring2)
     header = "HEADER - id: " // chaine(29:30) // " - thermalize - Temperature : " &
     &   // textstring // " - Total energy: " // textstring2
     call writepdb(fname,natoms,posa,header)
  else if (etype .eq. 'production') then

     if (SINGLEFILE) then
       fname = trim(aname)
     else
       ! We need to write the configuration in a min.... file
       call convert_to_chain(counter,scounter,"0")
       fname =   trim(aname) // scounter(31-ndigits:30) // PDB_EXT
     endif
     if(usextc) then
       if(openf) then
         if(SIMULATION_TYPE .eq. "tempering") then
           fname = "min.xtc" //char(0)
         else
           fname = trim(aname) // ".xtc" // char(0)
         endif
         mname = "a" // char(0)
         call xdropen(xd, fname,mname)
         openf = .false.
       endif
     call writextc(xd,natoms,counter,real(value),box, posa, prec,ret)
!     call xdrclose(xd,ret)
     else
       call convert_to_chain(id,chaine,"0")
       call real_to_chain(value,10,4,textstring)
       call real_to_chain(energy,16,6,textstring2)
       header = "HEADER - id: " // chaine(29:30) // " - simulation time (ns) : " // textstring // " - Total energy: " // textstring2
       call writepdb(fname,natoms,posa,header)
     endif
  else
     write(*,*) 'etype : ', etype
     stop 'Wrong type in write_to_file'
  endif

end subroutine write_to_file

end module writetofile
