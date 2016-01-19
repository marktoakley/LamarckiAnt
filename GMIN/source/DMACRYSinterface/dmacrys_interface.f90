! TODO: ID6 calculation added +6 for nag
module dmacrys_interface
      integer :: num_atoms, num_rigid, num_lattice
      double precision, allocatable :: local_axes(:,:,:)
      double precision, allocatable :: cell_shift(:)

      Real (kind=8), allocatable, private :: w1(:,:)
      Real (kind=8), allocatable, private :: g1(:,:)
      Real (kind=8), allocatable, private :: wrk(:)
      Real (kind=8), allocatable, private :: wrk1(:)
      Real (kind=8), allocatable, private :: del(:)
      Real (kind=8), allocatable, private :: wstiff(:,:)
      Real (kind=8), allocatable, private :: gstiff(:,:)  ! gradient for molecule
      Real (kind=8), allocatable, private :: wrkpol(:)
      Real (kind=8), allocatable, private :: gfield(:)
      Real (kind=8), allocatable, private :: wrkelc(:)
      Real (kind=8), allocatable, private :: s(:,:)
      Real (kind=8), allocatable, private :: t(:,:)
      Real (kind=8), allocatable, private :: u(:,:)
      Real (kind=8), allocatable, private :: v(:,:)
      Real (kind=8), allocatable, private :: w(:)
      Real (kind=8), allocatable, private :: wa1(:,:)
      Real (kind=8), allocatable, private :: ga1(:,:)
      Real (kind=8), allocatable, private :: delsymm(:)
      Real (kind=8), allocatable, private :: gg(:,:)
      logical, private :: LOPEN,LFIX,LPRINT,LFRCD
      logical, private :: ITWSW
      integer, private ::  ID6, icomm,id7,id5i
      integer, private :: iupd ! not used
      logical :: use_remove_overlap_potential
      double precision ener

      double precision coord_scaling
      double precision initial_volume

      public dmacrys_get_lattice
      logical dmadidstep

      integer interaction_flags

      integer flag_coulomb, flag_shortrange, flag_3body, flag_torsions, flag_ellipsoid, flag_dma
      parameter(flag_coulomb=0)
      parameter(flag_shortrange=1)
      parameter(flag_3body=2)
      parameter(flag_torsions=3)
      parameter(flag_ellispoid=4)
      parameter(flag_dma=5)

contains

! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! Initialize all the necessary variables in dmacrys for energy and
! derivative calculations. This routine assumes that from now on
! everything in the data/odata file is DMACRYS input
! ----------------------------------------------------------------------
subroutine dmacrys_initialize
    include 'plutoincludes'
    integer imol, i
!      INCLUDE 'include/cntrl'
      LOGWR=.TRUE.
      logwr_=.true.
!      IF(LOGWR)WRITE(6,100)
!   call CASCAD here. DMACRYS reads input file and basically
!   runs a whole setup and optimization if specified in input. After DMACRYS returns
!   all structures should be initialized and can be tweaked by us
    print *,"DMACRYS> running DMACRYS parser"
    call CASCAD
    print *,"DMACRYS> DMACRYS returned from CASCAD, starting GMIN setup"

! some checks about stuff I don't know
    if(NUMION /= 0) Then
        print *,"NUMION is not zero, what does it mean?"
        stop
    endif
    if(NSHL /= 0) Then
        print *,"NSHL is not zero, what does it mean?"
        stop
    endif
    if(NEFREE /= 0) Then
        print *,"NEFREE is not zero, what does it mean?"
        stop
    endif

    if(NELIP /= 0) Then
        print *,"NELIP is not zero, what does it mean?"
        stop
    endif
    if(NUMION /= 0) Then
        print *,"NUMION is not zero, what does it mean?"
        stop
    endif

!   Run some stuff fromdmacrys_interface% a pluto simulation. This sets up
!   everything for calculating energies and gradients
    call copy_paste_from_pluto_init
!    print *,"number atoms",natoms,MOLS,
    num_atoms = NFREE !natoms
    num_rigid = MOLS
    num_lattice = 6

    ! allocate memory to store local axes
    allocate(local_axes(3,3,num_rigid))
    do i = 1, num_rigid
        local_axes(1:3,1,i) = EDM(1:3,1,i)
        local_axes(1:3,2,i) = EDM(1:3,2,i)
        local_axes(1:3,3,i) = EDM(1:3,3,i)
    end do

    print *,"DMACRYS> interface ready to calculate energies"
    call dmacrys_energy
    print *,"DMACRYS> the current energy is", ener
    dmadidstep=.true.

    use_remove_overlap_potential = .false.
    interaction_flags = 63
    IF(NOION) interaction_flags = IBCLR(interaction_flags, flag_coulomb)
    !flag_coulomb, flag_shortrange, flag_3body, flag_torsions, flag_ellipsoid, flag_dma
    !interaction_flags = IBCLR(interaction_flags, flag_coulomb)
    !interaction_flags = IBCLR(interaction_flags, flag_shortrange)
    !interaction_flags = IBCLR(interaction_flags, flag_3body)
    !interaction_flags = IBCLR(interaction_flags, flag_torsions)
    !interaction_flags = IBCLR(interaction_flags, flag_dma)
end subroutine

! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! calculate the potential
! this updates the atom positions based on rigid body coordinates
! and adjusts (hopefully) all the dmacrys structures
! ----------------------------------------------------------------------
subroutine dmacrys_calc_potential(rigid_coords,grad,ereal,gradt)
    use system, only: atoms ! from DMACRYS
    use commons ! from GMIN
    use genrigid
    use cell_reduction
    use polarizabilities
    use dma_overlap
    implicit none
    include 'include/maxbas'
    include 'include/maxmol'
    include 'include/ivec'
    include 'include/dma'
    include 'include/molcul'
    include 'include/rtable'
    include 'include/lattis'

    double precision, intent(in)  :: rigid_coords(1:3*NATOMS) ! that is the natoms hack
    double precision :: rigid_coords_n(DEGFREEDOMS) ! that is the natoms hack
    double precision, intent(out) :: ereal
    double precision, intent(out) :: grad(1:3*NATOMS)
    logical, intent(in) :: GRADT

    double precision :: lattice_coords(6)
    double precision :: mlattice(3,3)
    double precision xtmp(1:3*NATOMS)
    character(len=10) :: enerstr

    double precision R_GMIN(3,3), G_dma(3), R_tmp(3,3), E(6)
    double precision a(3), b(3), c(3)

    integer i, offset, imol, k
    double precision :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      ! the 6 elements of the symmetric strain matrix
      !  ( E1    E6/2  E5/2 )
      ! (  E6/2  E2    E4/2  )
      !  ( E5/2  E4/2  E3   )
      double precision strain(6)


    lattice_coords(1:6) = rigid_coords(6*num_rigid+1:6*num_rigid+6)
    !print *,lattice_coords(1:6)
!    print *,"strain"
!    print *,strain
    ! construct the strain transformation

    call get_lattice_matrix(lattice_coords, mlattice)

    XLAT=mlattice * coord_scaling
!    print *,"XLAT in calc pot"
!    print *,XLAT
!    print *, "angles", dot_product(XLAT(:,1),XLAT(:,3))
!    print *,"lattice"
!    print *,xlat
! call DMACRYS routines to update new coordinates
! lattice
    CALL RCPLAT(0)
    CALL PRMRST

    rigid_coords_n(1:DEGFREEDOMS) = rigid_coords(1:DEGFREEDOMS)
    do i=1,3*num_rigid
      rigid_coords_n(i)=rigid_coords_n(i) - entier(rigid_coords_n(i))
    end do
!   build cartesian coordinates
    rigid_coords_n(1:3*num_rigid) = rigid_coords_n(1:3*num_rigid) + cell_shift(1:3*num_rigid)
    call TRANSFORMRIGIDTOC (1, num_rigid, xtmp, rigid_coords_n)
    rigid_coords_n(1:3*num_rigid) = rigid_coords_n(1:3*num_rigid) - cell_shift(1:3*num_rigid)

    if(use_remove_overlap_potential) then
        ereal = overlap_potential_wrapper(xtmp, grad, rigid_coords_n, mlattice, gradt)
        return
    endif
!    open(unit=333,file="latest.raw")
!    write (333,*) xtmp
!    close(333)

! now copy atom coordinates
!    print *,num_atoms
 !   print *,"title"
    do i=1,num_atoms
        atoms(i)%coord(1:3) = xtmp(3*(i-1)+1:3*i) * coord_scaling
!        WRITE (*,'(A,3F8.4)') "A ",xtmp(3*(i-1)+1:3*i)
    end do
!    print *,"end"
    ! adjust center of mass
    do imol = 1, num_rigid
        COFMAS(1:3,imol) = matmul(XLAT,rigid_coords_n(3*imol-2:3*imol))
    end do

! now update local axis systems of molecules
    do i = 1, num_rigid
        offset   = 3*NRIGIDBODY + 3*i
        P(:) = rigid_coords(offset-2:offset)

        ! get the rotation matrix in angle axis (RMI)
        ! we don't need derivatives (DRMI?) here
        CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .false.)

        ! TODO: check, is this correct or do we have to use transform RMI?
        EDM(1:3,1,i) = MATMUL(RMI(:,:), local_axes(1:3,1,i))
        EDM(1:3,2,i) = MATMUL(RMI(:,:), local_axes(1:3,2,i))
        EDM(1:3,3,i) = MATMUL(RMI(:,:), local_axes(1:3,3,i))
    end do

    ! adjust axis system for dmacrys rigid bodies
    do imol = 1,num_rigid
      do I=MEMBRS(imol)+1,MEMBRS(imol+1)
        RSLJ(1,I)= 0.0D+0
        RSLJ(2,I)= 0.0D+0
        RSLJ(3,I)= 0.0D+0
        do K=1,3
            RSLJ(1,I)= RSLJ(1,I)+SLJ(K,I)*EDM(1,K,imol)
            RSLJ(2,I)= RSLJ(2,I)+SLJ(K,I)*EDM(2,K,imol)
            RSLJ(3,I)= RSLJ(3,I)+SLJ(K,I)*EDM(3,K,imol)
        end do
      end do
    end do


! update coordinates in DMACRYS
    CALL CATOCR

    if(INDUCE) print*,"do induce"
!   whatever this is
    IF ((LHPOLS.AND.ICOMM.LE.2).AND..NOT.INDUCE) CALL TBLCNT

! now calculate potential
    call dmacrys_energy
    ereal = ener
!    WRITE (enerstr,'(F8.6)') ereal

    !call dmacrys_dump_cif("conf_"//trim(enerstr)//".cif", rigid_coords)


    if(.not.GRADT) return
! now back-transform gradient
! rigid body translation
    !print *,gstiff(1:30,1)
    grad(1:6*num_rigid+6)=0

    ! copy center of mass forces and transform to initial lattice
    do i=1,num_rigid
        grad(3*i-2:3*i) = matmul(transpose(mlattice), gstiff(3*i-2:3*i,1)) * coord_scaling
        !grad(3*i-2:3*i) = gstiff(3*i-2:3*i,1) * coord_scaling
    end do
    !grad(1:3*num_rigid) = gstiff(1:3*num_rigid,1) * coord_scaling

    ! orientational gradients
    do i=1,num_rigid
        offset = 3*num_rigid + 3*i-2

        ! gradient from DMACRYS, given for rotations
        ! along unrotated state (P_gmin = (0,0,0))
        G_dma(1:3) = gstiff(offset:offset+2,1)
        P(:) = rigid_coords(offset:offset+2)

        ! gradients have to be converted to gradients in angle axis
        ! coordinates P
        ! dU / dPk = \sum_i dU/dTheta_i dTheta(dPk)
        ! R_k R^(-1) = dTheta
        !
        ! with (no guarantee for the sign here)
        !
        !            (      0    -dTheta3  dTheta2   )
        ! dTheta =  (    dTheta3     0    -dTheta1    )
        !            (  -dTheta2  dTheta1     0      )
        !
        !
        !
        CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .true.)
        CALL RMDRVT(-P, RMI, DRMI1, DRMI2, DRMI3, .false.)

        DRMI1=matmul(DRMI1,RMI)
        DRMI2=matmul(DRMI2,RMI)
        DRMI3=matmul(DRMI3,RMI)
        grad(offset+0) = 0.5*G_dma(1)*(DRMI1(3,2) - DRMI1(2,3)) &
                       + 0.5*G_dma(2)*(DRMI1(1,3) - DRMI1(3,1)) &
                       + 0.5*G_dma(3)*(DRMI1(2,1) - DRMI1(1,2))
        grad(offset+1) = 0.5*G_dma(1)*(DRMI2(3,2) - DRMI2(2,3)) &
                       + 0.5*G_dma(2)*(DRMI2(1,3) - DRMI2(3,1)) &
                       + 0.5*G_dma(3)*(DRMI2(2,1) - DRMI2(1,2))
        grad(offset+2) = 0.5*G_dma(1)*(DRMI3(3,2) - DRMI3(2,3)) &
                       + 0.5*G_dma(2)*(DRMI3(1,3) - DRMI3(3,1)) &
                       + 0.5*G_dma(3)*(DRMI3(2,1) - DRMI3(1,2))

    end do
    E=gstiff(6*num_rigid+1:6*num_rigid+6,1)! / coord_scaling

    a = mlattice(:,1)
    b = mlattice(:,2)
    c = mlattice(:,3)

    ! a1
    grad(6*num_rigid+1) = E(1)/a(1)
    ! a23
    grad(6*num_rigid+2) = E(6) / a(1)
    ! a3
    grad(6*num_rigid+3) = E(5) / a(1)
    ! b2
    grad(6*num_rigid+4) = E(2) / b(2)  - E(6) * a(2) / (a(1) * b(2))
    ! b3
    grad(6*num_rigid+5) = - E(5) * a(2) / (a(1) * b(2)) + E(4) / b(2)
    ! c3
    grad(6*num_rigid+6) = E(5) * (a(2) * b(3) - a(3) * b(2))/(a(1)*b(2)*c(3)) &
     - b(3)/(b(2)*c(3))*E(4) + E(3) / c(3)

    do i=1,num_rigid
        offset = 3*num_rigid + 3*i-2

        grad(6*num_rigid+2) = grad(6*num_rigid+2) &
                    - 0.5d0/a(1) * gstiff(offset+2,1)
        grad(6*num_rigid+3) = grad(6*num_rigid+3) &
                    + 0.5d0/a(1) * gstiff(offset+1,1)

        grad(6*num_rigid+4) = grad(6*num_rigid+4) &
                    + 0.5d0*a(2) / (a(1)*b(2)) * gstiff(offset+2,1)
        grad(6*num_rigid+5) = grad(6*num_rigid+5) &
                    - 0.5d0 * a(2) / (a(1)*b(2)) * gstiff(offset+1,1) &
                    - 0.5d0 /b(2) * gstiff(offset+0,1)

        grad(6*num_rigid+6) = grad(6*num_rigid+6) &
                    + 0.5d0*(a(2)*b(3) - a(3)*b(2))/(a(1)*b(2)*c(3)) * gstiff(offset+1,1) &
                    + 0.5d0*b(3)/(b(2)*c(3))* gstiff(offset+0,1)
    end do
    grad = 2*grad
    ! we abuse RMI here for lattice calculation
!    call strain_to_matrix(E, RMI)
    !print *,E
    !print *,RMI
    !print *,mlattice
!    call invert3x3(mlattice, mtmp)
    !print *,matmul(mtmp, RMI)
    ! to lazy here to explicitly write it for strain values E1-E6
    !call matrix_to_strain(matmul(RMI, mlattice), grad(6*num_rigid+1:6*num_rigid+6))
    !print *,matmul(mlattice, RMI)
!    RMI = matmul(RMI,mtmp)-mtmp
!    RMI(1,1)=RMI(1,1)+1
!    RMI(2,2)=RMI(2,2)+1
!    RMI(3,3)=RMI(3,3)+1
!    call matrix_to_strain(RMI, grad(6*num_rigid+1:6*num_rigid+6))
    !grad(6*num_rigid+1:6*num_rigid+6) = E
    !!grad(6*num_rigid+1:6*num_rigid+6) = E
    !grad(6*num_rigid+1:6*num_rigid+6)=E
!    print *,"lattice gradient"
!    print *,grad(6*num_rigid+1:6*num_rigid+6)
!    print *,"irset"
!    print *,irset


end subroutine

subroutine strain_to_matrix(E, M)
    implicit none
    double precision, intent(in) :: E(6)
    double precision, intent(out) :: M(3,3)

    M(1,1) = 1d0 + E(1)
    M(2,2) = 1d0 + E(2)
    M(3,3) = 1d0 + E(3)

    M(1,2)=0.5d0*E(6)
    M(1,3)=0.5d0*E(5)
    M(2,3)=0.5d0*E(4)
    M(2,1)=M(1,2)
    M(3,1)=M(1,3)
    M(3,2)=M(2,3)
end subroutine

subroutine matrix_to_strain(M, E)
    implicit none
    double precision, intent(out) :: E(6)
    double precision, intent(in) :: M(3,3)

    E(1) = M(1,1) - 1d0
    E(2) = M(2,2) - 1d0
    E(3) = M(3,3) - 1d0

    E(6) = 2d0*M(1,2)
    E(5) = 2d0*M(1,3)
    E(4) = 2d0*M(2,3)
end subroutine

! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! initializes gmin structures and setup rigid body framework
! ----------------------------------------------------------------------
subroutine dmacrys_init_gmin
    use system !, only: atoms ! from DMACRYS
    use atom_library
    use commons ! from GMIN
    use genrigid
    use vec3
    use dma_overlap

    implicit none
      INCLUDE 'include/maxel'
      INCLUDE 'include/maxmol'
      INCLUDE 'include/maxbas'
      INCLUDE 'include/maxpot'
      INCLUDE 'include/maxtor'
      INCLUDE 'include/maxvec'
      INCLUDE 'include/maxdef'
      INCLUDE 'include/maxine'
      include 'include/ivec'
    include 'include/bascas'
    include 'include/dma'
    include 'include/molcul'
    include 'include/rtable'
    include 'include/lattis'
    include 'include/potent'

    double precision :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
    double precision :: mlattice(3,3)
    double precision :: rcoords(100)

    integer iatom,imol,i,counter,maxcounter

    print *, "DMACRYS> Filling GMIN structures"
    print *, "DMACRYS> Number of atoms: ", num_atoms
    print *, "DMACRYS> Number of rigid bodies: ", num_rigid
    print *, "DMACRYS> Number of lattice vectors: ", num_lattice

    print *, "DMACRYS> Initializing GMIN coordinates"

    coord_scaling = 1d0  / RLSCAL

    do i=1,num_atoms
        coords(3*(i-1)+1:3*i, 1) = atoms(i)%coord(1:3) / coord_scaling
        print *, "elem", atoms(i)%name
    end do
    print *, "DMACRYS> Initializing GMIN rigid bodies"
    RIGIDINIT=.true.
    has_lattice_coords=.true.
! vr274> count maximum number of atoms per rigid body
    maxcounter=0
    do imol=1,num_rigid
       counter=0
       do i=MEMBRS(imol)+1,MEMBRS(imol+1)
           iatom=MAPX(i)
           counter=counter+1
       end do
       maxcounter=max(maxcounter,counter)
    end do

    call GENRIGID_ALLOCATE(num_rigid, maxcounter)

    do i=1,num_atoms
        GR_WEIGHTS(i) = AMASS(iabs(lbas(i)))
    end do
   
    do imol=1,num_rigid
       NSITEPERBODY(imol) = MEMBRS(imol+1) - MEMBRS(imol)
       do i=1,NSITEPERBODY(imol)
           iatom=MAPX(i+MEMBRS(imol))
           RIGIDGROUPS(i,imol)=iatom
       end do
    end do

    ! this is for debugging purpose
    do imol=1,num_rigid
       !P(1) = 0 !3.141592654/4.
       !P(2) = 0
       !P(3) = 3.141592654 
       !CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .false.)
       call invert3x3(local_axes(1:3,1:3,imol), rmi)
       do i=1,NSITEPERBODY(imol)
           iatom=MAPX(i+MEMBRS(imol))
           P=coords(3*iatom-2:3*iatom, 1) - COFMAS(1:3,imol)
           coords(3*iatom-2:3*iatom, 1) = matmul(RMI, P) + COFMAS(1:3,imol)
       end do
       local_axes(1:3,1,imol) = MATMUL(RMI(:,:), local_axes(1:3,1,imol))
       local_axes(1:3,2,imol) = MATMUL(RMI(:,:), local_axes(1:3,2,imol))
       local_axes(1:3,3,imol) = MATMUL(RMI(:,:), local_axes(1:3,3,imol))
       !print *,local_axes
    end do

!           print *,i,imol
    CALL GENRIGID_INITIALISE(coords(:, 1))
    do i=1,num_atoms
        coords(3*(i-1)+1:3*i, 1) = atoms(i)%coord(1:3) / coord_scaling
    end do

    ! initialize the lattice coordinates
    coords(3*num_atoms+1, 1) = XLAT(1,1) / coord_scaling
    coords(3*num_atoms+2, 1) = XLAT(2,1) / coord_scaling
    coords(3*num_atoms+3, 1) = XLAT(3,1) / coord_scaling
    coords(3*num_atoms+4, 1) = XLAT(2,2) / coord_scaling
    coords(3*num_atoms+5, 1) = XLAT(3,2) / coord_scaling
    coords(3*num_atoms+6, 1) = XLAT(3,3) / coord_scaling

!    open(unit=333,file="input.raw")
!    read (333,*) coords(:,1)
!    close(333)
!    print *,coords(:,1)
    call get_lattice_matrix(coords(3*num_atoms+1:3*num_atoms+6, 1), mlattice)
    print *,"Initial lattice vectors"
    print *,"a:", mlattice(:,1)
    print *,"b:", mlattice(:,2)
    print *,"c:", mlattice(:,3)

    initial_volume = dot_product(mlattice(:,1), vec_cross(mlattice(:,2),mlattice(:,3)))

    call dmacrys_determine_cell_shift(coords(:,1))
    if(DMACRYS_RANDOMSTART) then
        call dmacrys_genrandom(coords(:, 1))
    endif

!    open(unit=665, file="dmp")
!    read(665,*) coords(:,1)
!    close(665)
!    call transformctorigid(coords,rcoords(1:degfreedoms))
!    open(unit=665, file="dmp_rigid")
!    write(665,*) rcoords(1:degfreedoms)
!    close(665)

end subroutine

subroutine dmacrys_determine_cell_shift(coords)
    use genrigid
    use commons, only: natoms
    use vec3
    implicit none
    include 'include/maxbas'
    include 'include/maxmol'
    include 'include/ivec'
    include 'include/dma'
    include 'include/molcul'
    include 'include/rtable'
    include 'include/lattis'

    double precision :: coords(3*NATOMS)
    double precision :: rcoords(degfreedoms)
    double precision :: minv(3,3)
    double precision :: diff(3)
    integer i,j

    call transformctorigid(coords,rcoords)
    call invert3x3(XLAT,minv)

    allocate(cell_shift(3*NRIGIDBODY))

    do i=1,NRIGIDBODY
        diff = rcoords(3*i-2:3*i)-matmul(minv,COFMAS(1:3,i))
        do j=1,3
            if(abs(diff(j) - dint(diff(j)+0.1)) > 1e-8) then
                print *,"DMAGMIN ERROR: the center of mass calculation seems wrong"
                print *,"Molecule ", i
            endif
        end do
        cell_shift(3*i-2:3*i) = dint(diff)
    end do

end subroutine
! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! Call DMACRYS to calculate energies.
! This is called from dmacrys_calc_potential. All coordinates and
! DMACRYS structures have to be updated first!
! ----------------------------------------------------------------------
subroutine dmacrys_energy
    use cascad_export
    include 'plutoincludes'
INCLUDE 'include/fitidc'
INCLUDE 'include/lattis'
INCLUDE 'include/esq'
    integer ICTRL2
    double precision e2, to_eV
!    implicit none
    EXTERNAL RCBK,PERSCT,WTWOF,SYMTWO, &
      THBVBK,THBVBS,TORVBS,TORVBK,FRCSYM

!    if(allocated(JBARR)) print*,"JBARR is allocated"
!    if(allocated(IBARR)) print*,"IBARR is allocated"
 !allocate (ibarr(9,1),stat=ialloc(2))
!allocate (itarr(14,1),stat=ialloc(3))

    !print *, "dmacrys_energt"

  IGWSW = 1 ! Cell energy and 1st derivatives only.
  ICTRL2 = 2 ! Co-ordinate and strain derivatives.

    gstiff=0
  LOGWR=.FALSE. ! suppress output
  !print *,"fl",interaction_flags
  CALL FRCNST(JBARR,IBARR,ITARR, &
      W1,WSTIFF,G1(:,ING), &
      GSTIFF(:,ING),WRK,WRK1, &
      WRKELC,WRKPOL,GFIELD, &
      E2,ID1,ID5,ID6,LFRCD,ICTRL2,interaction_flags, &
      RCBK,PERSCT,WTWOF,SYMTWO,THBVBK,THBVBS,TORVBS,TORVBK,FRCSYM)

   to_eV = (ESQ/2.0D+0/RLSCAL)

   GSTIFF = GSTIFF*to_eV
   ener = e2*to_eV
   !print *,"energy",ener

   !print *,EQQS+EPOLSS,ESHRT,ETHB,ETOR,UES,EANIS!,PRESSURE_intunits*VOLUM
   !if(ener < -1000) THEN
!      ener=100
!      GSTIFF=0
!   endif
end subroutine


! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! This routine is the whole pluto.f90 initialization copy&pasted
! till the do loop for minimization begins with
!     DO WHILE (iter.le.2)
!
! The variable definitions were extracted and moved to the module.
! One if had to be removed
!    IF(.NOT.PLREAD)THEN
! ----------------------------------------------------------------------
subroutine copy_paste_from_pluto_init
    include 'plutoincludes'
INCLUDE 'include/fitidc'

    integer i, j,jsw

    print *,"----- Running copy_paste_from_pluto_init"
10  FORMAT(1X,'ERROR - THE DYNAMIC MEMORY REQUIREMENT OF PLUTO', &
    ' SECTION HAS EXCEEDED THE ',I8,' WORDS OF MEMORY AVAILABLE.', &
    ' PROGRAM TERMINATING.')


    !
    !     LOOP OVER RANGES AND NUMBER OF POTENTIALS TO FILL CUTOFFS^2
    !
    DO 641 I=1,MAXRNG
        DO 6411 J=1,MAXPOT
            RMAX2(I,J)=RMAX(I,J)*RMAX(I,J)
6411    CONTINUE
641 CONTINUE
    !
    CUTPO2=CUTPOT*CUTPOT
    CALL INPT(JSW)
    IF(JSW.NE.1)THEN
        CALL OSETPL
    ELSE
        CALL SETPL
        !CCC
        !CCC  REMEMBER ORIGINAL STRUCTURE FOR COMPARISON LATER
        !CCC
        CALL REMEM
        IUPD=0
    ENDIF
    !
    IF (.NOT.LCLUST) CALL INTRA2
    LFIX=ISTR.GE.3
    CALL CHSET(Q,QSORT,LDM,LFIX)
    IPLWR(1)=MPRINT(5)-MPRINT(5)/10*10
    IPLWR(2)=MPRINT(5)/10-MPRINT(5)/100*10
    IPLWR(3)=MPRINT(5)/100-MPRINT(5)/1000*10
    IPLWR(4)=MPRINT(5)/1000-MPRINT(5)/10000*10
    IPLWR(5)=MPRINT(5)/10000
    !
    !     COMMENTED OUT BIT 1 MOVED TO BOTTOM OF FILE
    !
    !     SET MAXIMUM RANGE OF ALL POTENTIALS TO CUTPOT + A SMALL AMOUNT
    !     DO NOT DO FOR QUINTIC SPLINE, DONT BOTHER IF NRANGE ZERO!
    !
    DO 101 I=1,NBP
        IF(ISWPOT(1,I).NE.15.AND.NRANGE(I).NE.0)THEN
            RMAX(NRANGE(I),I)=CUTPOT+1.0D+0
        ENDIF
101 CONTINUE
    IF(LSYMM)THEN
        ID7=INDREP(6,IREPA1)
    ELSE
        ID7=1
    ENDIF
    !
! ------------------------------------------------
!    vr274: commented this out from original pluto
!    IF(.NOT.PLREAD)THEN
! ------------------------------------------------
        !
        !     Set dynamic core control arrays for properties calculation.
        !
        IF(IGWSW.GT.0.or.induce)THEN
            !
            !     Need derivatives case
            !
            IF(LFIX)THEN
                ID1=3*NFREE+6+6*NEFREE
                ID2=3*NSHL
                ID4=3*NBAS
                ID5=3*MOLION+3*(MOLION-NUMION)+6
                ID5i=3*MOLION+3*(MOLION-NUMION)+6
                ID6=6*NFREE+6+6*NEFREE+6
            ELSE
                !ccc
                !ccc  DJ Willock's changes to include GSTIFF and WSTIFF arrays for rigid molecules
                !ccc  adapted for ions 7/93
                !ccc
                ID1=3*NBAS+6+6*NELIP
                ID2=3*NSHL
                ID4=3*NBAS
                ID5=3*MOLION+3*(MOLION-NUMION)+6
                ID5i=3*MOLION+3*(MOLION-NUMION)+6
                ID6=6*NBAS+6+6*NELIP+6
            ENDIF
        ELSE
            !
            !     energy only calc case.
            !
            ID1=1
            ID2=1
            ID4=1
            ID5=1
            ID5i=1
            ID6=1
        ENDIF
!        print *,"id5 1", id5
!        print *,"id6 1", id6
        allocate (w1(id1,id1),stat=ialloc(1))
        allocate (g1(id1,2),stat=ialloc(2))
        allocate (wrk(6*id1),stat=ialloc(4))
        allocate (wrk1(id1),stat=ialloc(5))
        allocate (del(id1),stat=ialloc(6))
        allocate (wstiff(id5,id5),stat=ialloc(7))
        allocate (gstiff(id5i,2),stat=ialloc(8))
        allocate (wrkpol(id6),stat=ialloc(10))
        allocate (wrkelc(id1),stat=ialloc(11))
        allocate (s(id5,id5),stat=ialloc(12))
        allocate (gfield(id1),stat=ialloc(13))
        if(.not.lsymm)then
            lg=1
        endif
        allocate (t(id5,lg),stat=ialloc(14))
        allocate (u(lg,lg),stat=ialloc(15))
        allocate (v(lg,id5),stat=ialloc(16))
        allocate (w(id5),stat=ialloc(17))
        allocate (wa1(id7,id7),stat=ialloc(18))
        allocate (ga1(id7,2),stat=ialloc(19))
        allocate (delsymm(id7),stat=ialloc(20))
        allocate (gg(id5,8),stat=ialloc(21))
        call check_alloc(10,21)

        !
        INQUIRE(UNIT=8,OPENED=LOPEN)
        IF(.NOT.LOPEN)THEN
            OPEN(UNIT=8,FORM='FORMATTED')
        ENDIF
        !
        !     Calculate the first and second derivatives of the lattice energy
        !     with respect to the unit cell coordinates and bulk strains.
        !
        DO 990 I=1,NBAS
            DO 9901 J=1,3
                BASTRN(I,J)=0.0D+0
9901        CONTINUE
990     CONTINUE
        DO 993 I=1,6
            BKSTRN(I)=0.0D+0
993     CONTINUE
        IEXIT=0
        !
        ICOMM=ICOM
end subroutine

! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! Calculate derivatives using finite differences
! ----------------------------------------------------------------------
subroutine dmacrys_numerical_derivative(n,x,grad)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: x(n)
    double precision, intent(out)   :: grad(n)

    integer :: i
    double precision :: grad_tmp(n)
    double precision e1,e2
    double precision tmp
    double precision :: eps = 1d-7
    print *,n
    do i=1,6*num_rigid+6
        tmp = x(i)
        x(i)=x(i)-0.5*eps
        call dmacrys_calc_potential(x,grad_tmp,e1,.false.)
        x(i)=tmp+0.5*eps
        call dmacrys_calc_potential(x,grad_tmp,e2,.false.)
        x(i)=tmp !- 0.5*eps

        grad(i) = (e2-e1) / eps / 2d0
    end do
end subroutine

! ----------------------------------------------------------------------
! --------------------------- vr274 ------------------------------------
! Routine to check analytical derivatives against numerically calculated
! ones. This should be run from time to time for new systems to ensure
! that all used DMACRYS functionality is properly interfaced
! ----------------------------------------------------------------------
subroutine dmacrys_check_derivatives(n,x,grad)
    implicit none

    integer, intent(in) :: n
    double precision, intent(inout) :: x(n)
    double precision, intent(in)    :: grad(n)

    double precision grad_num(n) ! temporary array for numerical gradient
    double precision e1
    !call dmacrys_calc_potential(x,grad_num,e1,.false.)

    print *,"coords",num_rigid
    print *,grad(1:6*num_rigid+6)
    print *,"BEGIN NUMDIFF"
    call dmacrys_numerical_derivative(n,x,grad_num)
    print *,grad_num(1:6*num_rigid+6) !- grad_num(1:6*num_rigid+6)
    stop
!    print *,"analytical"
!    print *,grad(1:3*num_rigid)
!    print *,"-"
!    print *,grad(3*num_rigid+1:6*num_rigid)
!    print *,"-"
!    print *,grad(6*num_rigid+1:6*num_rigid+6)
!    print *,"numerical"
!    print *,grad_num(1:3*num_rigid)
!    print *,"-"
!    print *,grad_num(3*num_rigid+1:6*num_rigid)
!    print *,"-"
!    print *,grad_num(6*num_rigid+1:6*num_rigid+6)
    print *,"END NUMDIFF"

end subroutine


! Calculate the goldstein rotation matrix
!subroutine dmacrys_goldstein(P, A)
!    implicit none
!    double precision theta
!    double precision , intent(in) :: P(3)
!    double precision , intent(out) :: A(3,3)
!    double precision sint, cost
!
!    double precision ratio
!    double precision Ps(3)
!    double precision tmp
!    integer i
!    integer ii, jj, sgn
!
!    theta=0.0D+0
!
!    theta = sqrt(dot_product(P, P))
!    sint = sin(theta/2.0D+0)
!    cost = cos(theta/2.0D+0)
!
!    if (theta > 1.0D-06) then
!        ratio = sint/theta
!    else
!        ratio=0.5D+0-theta*theta/48.0D+0
!    endif
!
!    Ps(1:3)=ratio*P(1:3)
!
!    tmp = cost*cost - dot_product(Ps, Ps)
!
!    ! calculate the diagonals
!    do i=1,3
!        A(i,i)= tmp + 2.*Ps(i)*Ps(i)
!    end do
!
!    sgn=1
!    do I=1,3
!        ii=(i-1)/2+1
!        jj=(i+2)/2+1
!        sgn=-1*sgn
!        A(ii,jj)= 2.0D+0*(Ps(ii)*Ps(jj) + dble(sgn)*cost*Ps(4-i))
!        A(jj,ii)= 2.0D+0*(Ps(ii)*Ps(jj) - dble(sgn)*cost*Ps(4-i))
!    end do
!end subroutine

subroutine dmacrys_dump_cif(filename, rigid_coords, cifname)
    use commons ! from GMIN
    use genrigid
    use atom_library, only : elem_libr
    use system, only : atoms
    implicit none
    double precision mlattice(3,3)
    character*(*) :: filename, cifname
    double precision, intent(in)  :: rigid_coords(1:3*NATOMS) ! that is the natoms hack
    double precision xtmp(1:3*NATOMS)

    double precision PI
    parameter (PI=3.141592654D0)


! calculate the current lattice vectors
    call get_lattice_matrix(rigid_coords(6*num_rigid+1:6*num_rigid+6), mlattice)

!   build cartesian coordinates
    call TRANSFORMRIGIDTOC (1, num_rigid, xtmp, rigid_coords)

    call dmacrys_dump_cif_no_transform(filename, NATOMS-2,xtmp(1:3*NATOMS-6), mlattice, cifname)
end subroutine

subroutine dmacrys_dump_cif_no_transform(filename, natoms, atomcoords, lattice, cifname)
    use cell_reduction
    use atom_library, only : elem_libr
    use porfuncs, only: flush
    use system, only : atoms
    use vec3
    implicit none

    character*(*) :: filename
    character*(*) :: cifname
    integer natoms
    double precision atomcoords(1:3*NATOMS)
    double precision lattice(3,3), lattice_(3,3)
    double precision a(3), b(3), c(3)
    double precision invlattice(3,3)

    double precision PI
    parameter (PI=3.141592654D0)
    integer out, i

    out=822
    open(unit=out,file=filename)

    !call cell_set_lattice(lattice_)
    !call cell_get_lattice_for_cif(lattice)

! the lattice vectors
    a = lattice(:,1)
    b = lattice(:,2)
    c = lattice(:,3)

! title section
    write (out,'(A,A)') "data_", cifname
! symmetry, we always use triclinic so far
    write(out,'(A)') "_symmetry_cell_setting           triclinic"
    write(out,'(A)') "_symmetry_space_group_name_H-M   'P 1'"
    write(out,'(A)') "_symmetry_Int_Tables_number      1"
    write(out,'(A)') "loop_"
    write(out,'(A)') "_symmetry_equiv_pos_site_id"
    write(out,'(A)') "_symmetry_equiv_pos_as_xyz"
    write(out,'(A)') "1 x,y,z"

! write cell parameters
    write(out, '(A,F8.4)') "_cell_length_a                   ", vec_len(a)
    write(out, '(A,F8.4)') "_cell_length_b                   ", vec_len(b)
    write(out, '(A,F8.4)') "_cell_length_c                   ", vec_len(c)

    write(out, '(A,F8.4)') "_cell_angle_alpha                ", vec_angle(b, c)/PI*180D0
    write(out, '(A,F8.4)') "_cell_angle_beta                 ", vec_angle(a, c)/PI*180D0
    write(out, '(A,F8.4)') "_cell_angle_gamma                ", vec_angle(a, b)/PI*180D0
    write(out, '(A,F8.6)') "_cell_volume                     ", 0d0

! coordinates
    write(out,'(A)') "loop_"
    write(out,'(A)') "_atom_site_label"
    write(out,'(A)') "_atom_site_type_symbol"
    write(out,'(A)') "_atom_site_fract_x"
    write(out,'(A)') "_atom_site_fract_y"
    write(out,'(A)') "_atom_site_fract_z"

!    lattice(1,2:3)=0.
!    lattice(2,1)=0.
!    lattice(2,3)=0.
!    lattice(3,1:2)=0.
!    print *,"lattice",lattice
    call invert3x3(lattice, invlattice)
    do i=1,NATOMS
        write(out,*) atoms(i)%name, " ", atoms(i)%name(1:1), matmul(invlattice, atomcoords(3*i-2:3*i))
    end do
    call flush(out)
    close(out)
end subroutine

subroutine remove_overlap(xcoords)
    use commons !, only: natoms
    use dma_overlap
    implicit none
    double precision, intent(inout):: xcoords(3*natoms)

    integer iter
    double precision ereal
    logical cflag

    !coords(:,1) = xcoords
    !print *, "mindist", minimum_distance(coords)
    use_remove_overlap_potential=.true.
    CALL MYLBFGS(3*natoms,MUPDATE,xcoords,.FALSE.,1d-3,CFLAG,EREAL,100,ITER,.TRUE.,1)
    !print *, "mindist", minimum_distance(coords)
    !print *, "done", iter
    use_remove_overlap_potential=.false.
end subroutine

end module dmacrys_interface
