!  This subroutine works when we want to have different initial files
!  for replica exchange, the initial files must be in "conf_init" folder
!  with the names conf01 to conf02 ... correspond to the replica's id 
!  (see also lines 153 to 160 in protein-2006.f  and 131-136 in read_parameters.f90 )
!  Rozita Laghaei 1st Aug. 2007
!
subroutine init_remd(path_inits,single_init)

  use defs
  implicit none
  
  character(len=30) chain_Tid,path_inits
  logical      single_init

if ((.not. init_single_file)  .and. (SIMULATION_TYPE .eq. 'replica') ) then
   single_init = .false.
   call convert_to_chain(T_id,chain_Tid,"0")
   path_inits = trim("conf_init")// '/' //"conf"// chain_Tid(29:30)// PDB_EXT
else 
   single_init = .true.
endif

end subroutine init_remd
