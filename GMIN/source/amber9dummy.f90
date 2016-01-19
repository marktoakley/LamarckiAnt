SUBROUTINE SETUPAMB(filename)
implicit none
        CHARACTER(*) :: filename
        WRITE(*,*) 'sf344> shouldn`t be in this routine!  Compiled correctly?'
END SUBROUTINE SETUPAMB

SUBROUTINE AMBERINTERFACE(dummy1,dummy2,dummy3,dummy4)
implicit none
        INTEGER         :: dummy1, dummy2, dummy4
character(len=*) :: dummy3

WRITE(*,*)'dummy AMBERINTERFACE> AMBER9 keyword or coords.inpcrd file present but the executable is not AMBGMIN'
WRITE(*,*)'dummy AMBERINTERFACE> Stopping'
WRITE(*,*)'dummy AMBERINTERFACE> Please consult the Makefile, recompile and rerun'
STOP

END SUBROUTINE AMBERINTERFACE

subroutine amber_readcoords(ostring)
implicit none
character(len=20) ostring

end subroutine amber_readcoords

SUBROUTINE A9RESTOATOM(frozenres,frozen,nfreeze,unfreeze)
implicit none
logical frozenres(*), frozen(*), unfreeze
integer nfreeze
END SUBROUTINE A9RESTOATOM

SUBROUTINE AMBERFINALIO(nsave,nunit,node,ostring,filth,coords2)
implicit none

integer nsave,nunit,filth,node
character(len=*) ostring
double precision coords2(*)

END SUBROUTINE AMBERFINALIO

SUBROUTINE AMBERSECDER(OLDX,STEST)
implicit none
DOUBLE PRECISION  :: OLDX(*)
LOGICAL    :: STEST
END SUBROUTINE AMBERSECDER

SUBROUTINE AMBERENERGIES(y,grad,ereal,gradt,stest)
implicit none
double precision :: y(*),grad(*),ereal
logical          :: gradt,stest

END SUBROUTINE AMBERENERGIES

subroutine check_cistrans_rna(coords,natoms,atomlabels,goodstructure)
implicit none

integer,intent(in) :: natoms
character(len=5),intent(in) :: atomlabels(natoms)
double precision,intent(in) :: coords(3*natoms)
logical,intent(out) :: goodstructure

goodstructure=.false.

end subroutine check_cistrans_rna

subroutine check_cistrans_dna(coords,natoms,atomlabels,goodstructure)
implicit none

integer,intent(in) :: natoms
character(len=5),intent(in) :: atomlabels(natoms)
double precision,intent(in) :: coords(3*natoms)
logical,intent(out) :: goodstructure

goodstructure=.false.

end subroutine check_cistrans_dna

SUBROUTINE MME2WRAPPER(COORDS,ENERGY,VNEW,TEMPHESS,ATMASS,GRAD1)
implicit none

double precision coords(*),energy,vnew(*),temphess(*),atmass(*),grad1(*)

END SUBROUTINE MME2WRAPPER

SUBROUTINE MMEINITWRAPPER(prmtop,igb,saltcon,rgbmax,cut)
implicit none

character       :: prmtop
integer         :: igb
double precision :: saltcon,rgbmax,cut

END SUBROUTINE MMEINITWRAPPER

SUBROUTINE CHECK_CISTRANS_PROTEIN(coords,natoms,goodstructure,minomega,cisarray)
implicit none

! stuff for detecting cis-trans isomerisation of omega angles in peptide bonds
logical,intent(out) :: goodstructure
integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms),minomega
integer :: cisarray(natoms)

goodstructure=.false.

END SUBROUTINE CHECK_CISTRANS_PROTEIN

SUBROUTINE SET_CHECK_CHIRAL(coords,natoms,goodstructure,chiarray)
implicit none

! stuff for detecting cis-trans isomerisation of omega angles in peptide bonds
logical,intent(out) :: goodstructure
integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms)
integer :: chiarray(natoms)

goodstructure=.false.

END SUBROUTINE SET_CHECK_CHIRAL

subroutine check_chirality(coords,natoms,goodstructure)
implicit none

integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms)
logical :: goodstructure
end subroutine check_chirality

subroutine amber_deallocate_stacks()

end subroutine amber_deallocate_stacks

SUBROUTINE TAKESTEPAMBER(NP,COORDS,MOVABLEATOMLIST,NMOVABLEATOMS,LIGMOVET,MDSTEPT,RANDOMSEEDT, &
     &                                  BLOCKMOVET,NBLOCKS,ATOMSINBLOCK)
IMPLICIT NONE
DOUBLE PRECISION COORDS(*)
INTEGER NP, MOVABLEATOMLIST(*),NMOVABLEATOMS, NBLOCKS,ATOMSINBLOCK(*)
LOGICAL LIGMOVET, MDSTEPT, RANDOMSEEDT, BLOCKMOVET
END SUBROUTINE TAKESTEPAMBER

subroutine takestepamm(q, debug, bhstepsize)
implicit none
double precision :: q(3*5000)
logical :: debug
double precision :: bhstepsize
end subroutine takestepamm

subroutine  GETSEEDATM(dummy1,dummy2,dummy3,at1,at2,at3,IRES)

character :: dummy1, dummy2, dummy3
integer :: at1, at2, at3, IRES

end subroutine GETSEEDATM

subroutine GETICCOORDS()

end subroutine GETICCOORDS

subroutine CHECK_SIDECHAIN(i3,j3,k3,l3,IICD,IS_SIDECHAIN)
integer :: i3,j3,k3,l3,IICD
logical :: IS_SIDECHAIN(100)

end subroutine CHECK_SIDECHAIN

subroutine CHGETICVAL(Q,BOND1, BOND2, THET1, THET2, PHI)
double precision :: Q(3*1000), BOND1(1000), BOND2(1000), THET1(1000)
double precision :: THET2(1000), PHI(1000)
end subroutine CHGETICVAL 

SUBROUTINE A9DUMPPDB(DUMPCOORDS,FNAMEF)
CHARACTER :: FNAMEF
DOUBLE PRECISION :: DUMPCOORDS(*)
END SUBROUTINE A9DUMPPDB

SUBROUTINE INTEFINALIO(nsave,nunit,node,ostring,filth,coords2)
implicit none
integer nsave,nunit,filth,node
character(len=*) ostring
double precision coords2(*)
END SUBROUTINE INTEFINALIO

SUBROUTINE check_chirality_generic(P,NATOMS,GOODSTRUCTURE,init)
implicit none
integer natoms
double precision P(3*natoms)
logical goodstructure, init
END SUBROUTINE check_chirality_generic


