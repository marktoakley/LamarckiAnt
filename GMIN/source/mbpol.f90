module mbpolmod
implicit none

! does not compile with pgi due to sqrt - DJW
! double precision, parameter :: maxdOO=3.5D0**2, maxdOH=2.45D0**2, maxcosOOH=sqrt(3.D0)/2.D0
double precision, parameter :: maxdOO=3.5D0**2, maxdOH=2.45D0**2, maxcosOOH=0.8660254037844386D0
integer :: nnode

public :: mbpolinit, mbpol, mbpolstep
private
contains

subroutine mbpolinit()
use commons, only : natoms
implicit none
nnode=natoms/3
end subroutine mbpolinit

subroutine mbpol(coordinates, energy, gradient, gradt)
implicit none
         logical, intent(in)  :: gradt
double precision, intent(in)  :: coordinates(nnode*9)
double precision, intent(out) :: gradient(nnode*9), energy

if (gradt) then
  call mbpolenergygradient(nnode,energy,coordinates,gradient)
else
  call mbpolenergy(nnode,energy,coordinates)
endif

end subroutine mbpol

subroutine mbpolstep(coordinates,displacement,angle,translate,rotate)
implicit none
         logical, intent(in)    :: translate, rotate
double precision, intent(in)    :: displacement, angle
double precision, intent(inout) :: coordinates(nnode*9)

call mbpolhbondgeometric(coordinates)

end subroutine mbpolstep

subroutine mbpolhbondgeometric(coordinates)
use graph_mod, only: digraph_adj_cycle
implicit none
double precision :: coordinates(nnode*9)
integer :: adjacency(nnode,nnode)
double precision :: hydrogenbonds(nnode*6,2), dOO, dOH, cosOOH, dprand
integer :: mol1, mol2, nmol, hydrogen, oxygen
! arguments for digraphadjcycle
integer :: adj2(nnode,nnode), dad(nnode), order(nnode)

adjacency=0

do mol1=1,nnode
  do mol2=mol1+1,nnode
    dOO=sum((coordinates(mol1*9-8:mol1*9-6)-coordinates(mol2*9-8:mol2*9-6))**2)
    if (dOO<maxdOO) then
      dOO=sqrt(dOO)
      !1H1O2
      dOH=sum((coordinates(mol1*9-5:mol1*9-3)-coordinates(mol2*9-8:mol2*9-6))**2)
      if (dOH<maxdOH) then
        dOH=sqrt(dOH)
        cosOOH=dot_product((coordinates(mol1*9-8:mol1*9-6)-coordinates(mol2*9-8:mol2*9-6)),&
                (coordinates(mol1*9-8:mol1*9-6)-coordinates(mol1*9-5:mol1*9-3)))/(dOH*dOO)
        if (cosOOH<maxcosOOH) then
          write(*,*)(mol1-1)*3,(mol1-1)*3+1,(mol2-1)*3
          adjacency(mol1,mol2)=(mol1-1)*3+1
        endif
      endif
      !2H1O2
      dOH=sum((coordinates(mol1*9-2:mol1*9-0)-coordinates(mol2*9-8:mol2*9-6))**2)
      if (dOH<maxdOH) then
        dOH=sqrt(dOH)
        cosOOH=dot_product((coordinates(mol1*9-8:mol1*9-6)-coordinates(mol2*9-8:mol2*9-6)),&
                (coordinates(mol1*9-8:mol1*9-6)-coordinates(mol1*9-2:mol1*9-0)))/(dOH*dOO)
        if (cosOOH<maxcosOOH) then
          write(*,*)(mol1-1)*3,(mol1-1)*3+2,(mol2-1)*3
          adjacency(mol1,mol2)=(mol1-1)*3+2
        endif
      endif
      !1H2O1
      dOH=sum((coordinates(mol2*9-5:mol2*9-3)-coordinates(mol1*9-8:mol1*9-6))**2)
      if (dOH<maxdOH) then
        dOH=sqrt(dOH)
        cosOOH=dot_product((coordinates(mol2*9-8:mol2*9-6)-coordinates(mol1*9-8:mol1*9-6)),&
                (coordinates(mol2*9-8:mol2*9-6)-coordinates(mol2*9-5:mol2*9-3)))/(dOH*dOO)
        if (cosOOH<maxcosOOH) then
          write(*,*)(mol2-1)*3,(mol2-1)*3+1,(mol1-1)*3
          adjacency(mol2,mol1)=(mol2-1)*3+1
        endif
      endif  
      !2H1O1
      dOH=sum((coordinates(mol2*9-2:mol2*9-0)-coordinates(mol1*9-8:mol1*9-6))**2)
      if (dOH<maxdOH) then
        dOH=sqrt(dOH)
        cosOOH=dot_product((coordinates(mol2*9-8:mol2*9-6)-coordinates(mol1*9-8:mol1*9-6)),&
                (coordinates(mol2*9-8:mol2*9-6)-coordinates(mol2*9-2:mol2*9-0)))/(dOH*dOO)
        if (cosOOH<maxcosOOH) then
          write(*,*)(mol2-1)*3,(mol2-1)*3+2,(mol1-1)*3
          adjacency(mol2,mol1)=(mol2-1)*3+2
        endif
      endif
    endif
  enddo
enddo

write(*,"(10(I2,1X))")adjacency
call digraph_adj_cycle(adjacency,nnode,nnode,adj2,dad,order)
write(*,*)
write(*,"(10(I2,1X))")adj2
write(*,*)
write(*,"(10(I2,1X))")dad
write(*,*)
write(*,"(10(I2,1X))")order

end subroutine mbpolhbondgeometric

end module mbpolmod
