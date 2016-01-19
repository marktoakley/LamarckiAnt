#include "../SLOOP_LIB/lib_header.h"
#include "../SLOOP_LIB/lib_atoms.h"
/* 
*      David Burke March 2010
* June 2011 Added an option to iterate through all combinations rather than random search. 
* June  9 2011 Stopped selecting Glycine and Alanine for rotamer search!!! 
* June 10 2011 Added an option to select filename 
*/

#define MaxGroups 20 

extern FRAGMENTS * AddDunbrackSidechains (FRAGMENTS * ,int ,int,int,int *) ;
void TransformLIG (RESIDUE * ThisRes, int StartAtom, float Angle) ;
int AdjustLIGANDAngles (FRAGMENTS *) ; 
extern float MY_RANDOM_FLOAT ( float );
extern int MY_RANDOM_INT ( int );
void OutputVariables (void);

void ReadLIGDEF ( char * filename ) ; 

char NewRotLabel[MaxLine]=""; 

int ExhaustiveSearch=-999 ; // Random or iterative exhaustive search
int Offset=0; // Offset for new structure number 
char filename[MaxLine];
float MinBFactor = 9999; 
int MaxNstructures; 
int MaxRepeats; 
int UsingLovell = 1 ; 
int LigRotamerStates = 4; 
int limit = 99999;
int MaxNligchange=9999;
int MaxNrotchange=9999;
int MaxNrotchangeOrig=9999;
int MaxNrotcombine=0;
int ReduceMaxRotChange=0;
int RandGroup=0;
FRAGMENTS * PDB,* LIGPDB, * ROTPDB;
int startres=-1;
int endres=-1;
int RotamerTries=0;
int MaxRotamerTries=1;
int TryHarder=0;
float RmsdThreshold=9999;
float VDWThreshold;
float Tweakchi=-1;     /** tweak chi angles from original **/
float Tweakxyz=-1;     /** tweak xyz coords from original **/
float Tweaklig=-1;     /** tweak LIG angles from original **/
int CentreResidue=-1;
int LigandResidue=-1;
float DistanceThreshold = 9999;
float AcceptDistanceThreshold = 9999;
int MaxNligchangeOrig=9999;
float RotThreshold=0.95;
float ChiThreshold=0.95;
float XYZThreshold=0.95;
float LIGThreshold=0.95;
int Outpdb=0; // output Pdb format file 
int Outxyz=0; // output XYZ format file 
int DEFgroups=-1;
float SelectedResidues[MaxRes];
int LIGDEFtable[MaxGroups][MaxHet] ; // List of atoms within a "GROUP" for Sugar rotations 
float LIGDEFscale[MaxGroups]; // List of scale factors for each group 
float LIGDEFtries[MaxGroups]; // List of attempts for exhaustive search 
float LIGDEFdoubles[MaxGroups]; // List of groups which are allowed double moves for exhaustive search 
char LIGDEFfile[MaxLine]="";
int MadeLigChange=0;
int MadeLigSubsetChange=0;
float DeltaXYZ[MaxRes];
int i;
char remark[MaxLine]="";
int LIGdouble1=-1;
int LIGdouble2=-1;

/********/
/* MAIN */
/********/

int  main (const int argc, char *argv[]) { 
void ParseCommandLine (const int , char **);

void help();
extern void OutputPDB ( FRAGMENTS * ,char * ,char * ) ;
extern FRAGMENTS * CopyFRAGMENT(FRAGMENTS *);
extern void FreeFRAGMENT(FRAGMENTS *);
int OutputXYZ ( FRAGMENTS * Pdb , char * filename) ;
extern void SetSLoopEnvironment(void);
extern void ReadRotamers(void);
extern void ReadPenultimateRotamers(void);
extern int UpdateDerivedPDB  ( FRAGMENTS *  ) ;
extern void BuildSidechains(FRAGMENTS *,int,int);
extern void build_hydrogens ( FRAGMENTS * , int  ,int ) ;
 extern void SetDummyResidue (void) ;

FRAGMENTS * AdjustSidechains ( FRAGMENTS * Structure, int startres, int endres)  ;
void srandom(unsigned int );

int res;
int repeats; 
int MaxNrotchange2;
int Nstructures; 

   globals.amberpdb = 1;
   globals.ignore_water = JIVE_TRUE ;
   globals.use_het = 3 ;
   options.nrotamers = MaxRotamers;
//   options.verbose = 3;

/* Sets up some environment variables */
   SetSLoopEnvironment();
// Set default for tweakxyz
   for(i=0;i<MaxRes;i++) { DeltaXYZ[i]=-999;}
   for(i=0;i<MaxGroups;i++) { LIGDEFdoubles[i]=-999;}

// Set default for Lovell library 
if ( UsingLovell == 0 ) { /* Reads the Dunbrack rotamer library */
   sprintf(globals.BbdDir,"/software/bbdep/bbdep02.May.sortlib");
}
if ( UsingLovell == 1 ) { /* Reads the Lovell rotamer library */
   sprintf(globals.BbdDir,"/home1/dave/work/PdbRotamerSearch/penultimate.lib");
}

//printf("startres= %d end %d\n",startres,endres); 
  SetDummyResidue () ;
  MaxNstructures=1;
  MaxRepeats=1;

   ParseCommandLine ( argc, argv);
    
if ( UsingLovell == 0 ) { /* Reads the Dunbrack rotamer library */
   ReadRotamers();
}
else if ( UsingLovell == 1 ) { /* Reads the Dunbrack rotamer library */
   ReadPenultimateRotamers();
}

   if ( PDB == NULL ) { printf("NO PDBFILE !!!!\n"); help();return;}
   UpdateDerivedPDB (PDB) ;   /** update stuff like phi/psi/chi angles **/
   
// Check some variables 
   if ( startres < 0 ) { startres=0;}
   if ( endres < 0 ) { endres=MaxRes;}
   if ( endres > PDB->numres ) { endres=PDB->numres;}
   if ( strlen(LIGDEFfile)>0 ) ReadLIGDEF(LIGDEFfile);
   if ( CentreResidue > PDB->numres ) { printf("ERROR residue selected for distance check %d is larger than numres %d\n",CentreResidue, PDB->numres); CentreResidue=-1;}
   if ( LigandResidue > PDB->numres ) { printf("ERROR residue selected for Ligand %d is larger than numres %d\n",LigandResidue, PDB->numres); LigandResidue=-1;}
   // If no ligand centre supplied, use centre residue 
   //if ( LigandResidue == -1 && CentreResidue != -1 ) { LigandResidue=CentreResidue;}

   //printf("srandom(%s)\n",time(NULL) );
   srandom( (unsigned int)time(NULL) );
   strcpy(PDB->Label,strtok(PDB->Label,"."));
   if ( strlen(NewRotLabel)==0) { strcpy(NewRotLabel,PDB->Label);}

   MaxNligchangeOrig=MaxNligchange;
   MaxNrotchangeOrig=MaxNrotchange;
   if ( MaxNrotcombine==1) MaxNligchangeOrig=MaxNrotchange;// combine both 

   LIGPDB = CopyFRAGMENT(PDB); 

   OutputVariables ();
  // Make MaxRepeats adjustments to the Ligand 
  while ( repeats < MaxRepeats ) { 
   MaxNligchange=MaxNligchangeOrig; // Reset to original max
   // For each ligand repeat ( ie double moves) reset the exhaustive ligand moves 
   for(i=0;i<MaxGroups;i++) { LIGDEFtries[i]=LigRotamerStates;  }
   // Reducing MaxNligchange and MaxNrotchange
   if ( ReduceMaxRotChange == 1) { 
      MaxNligchange=MY_RANDOM_INT(MaxNligchange);
      printf("Reduced MaxNligchange to %d\n",MaxNligchange);
   }
   // Adjust Ligand angles 
   if ( MaxNligchange > 0 ) {
     if ( ExhaustiveSearch < 0 ) { MadeLigChange=AdjustLIGANDAngles  ( LIGPDB ) ; }
     else { 
         //MadeLigSubsetChange=AdjustLIGANDAnglesSubset ( LIGPDB ) ;    // Allow double moves for a subset of the groups for exhaustive search
         //if ( MadeLigSubsetChange == 0 ) { break;}
     }
   }
   // if combined then reduce Nrotchange to use up lig moves 
   if ( MaxNrotcombine==1) MaxNrotchangeOrig=MaxNligchange;
   // Make a copy of PDB and MaxNrotchange for multiple structures to be produced
   ROTPDB=CopyFRAGMENT(LIGPDB); 
   // For every LIGAND adjustment, produce MaxNstructures rotamer searches 
   Nstructures=0;
   while ( Nstructures < MaxNstructures ) { 
     printf("Produced nstructures %d out of maxnstructures %d\n",Nstructures,MaxNstructures);
     MaxNrotchange=MaxNrotchangeOrig; // Reset to original max
     if ( ExhaustiveSearch>=0 ) MaxNligchange=MaxNligchangeOrig; // Reset to original max
     if ( ReduceMaxRotChange == 1) { 
        MaxNrotchange=MY_RANDOM_INT(MaxNrotchange);
        printf("Reduced MaxNrotchange to %d\n",MaxNrotchange);
     }
     // Adjust sidechain of selected residues
     printf("###############################################################\n");
     printf("###############################################################\n");
     printf("#########  About to Do AA rotamer search ######################\n");
     printf("###############################################################\n");
     printf("###############################################################\n");
     ROTPDB=AdjustSidechains ( ROTPDB,  startres,  endres)  ;
     if ( ROTPDB == NULL ) { 
       printf("ROTPDB is null MadeLigChange %d MaxNligchange %d \n ",MadeLigChange, MaxNligchange);
       //if ( MaxNligchangeOrig == 0 && ExhaustiveSearch > -1 ) { 
       // If no AA changes then stop unless ligand changes were made
       if ( MadeLigChange == 0 && ExhaustiveSearch > -1 ) { 
         printf("No new changes made \n"); 
         break;
       }
       printf("ROTPDB is null but copying ligpdb MadeLigSubsetChange==%d MadeLigChange %d \n ", MadeLigSubsetChange,MadeLigChange);
       //assume ligand move is a unique structure but only once!! Then allow new rotamer moves 
       ROTPDB=CopyFRAGMENT(LIGPDB); 
       MaxNrotchange=MaxNrotchangeOrig; // Reset to original max
     }
     // Rebuild xyz coords 
     sprintf(filename,"%s_newrot%d",NewRotLabel,repeats*MaxNstructures+Nstructures+Offset);
     printf("### LIST OF CHANGED SIDECHAINS %s\n",filename);
     (void) BuildSidechains(ROTPDB,startres,endres);
     printf("### END OF LIST OF CHANGED SIDECHAINS\n");
     /** Build the hydrogen atoms **/
     if ( options.hydrogens == JIVE_TRUE) (void) build_hydrogens(ROTPDB,startres,endres);
     Nstructures++;
     if ( Outpdb==1 ) { 
      sprintf(filename,"%s_newrot%d.pdb",NewRotLabel,repeats*MaxNstructures+Nstructures+Offset);
      printf("Writing Outputpdb %s %s\n",filename,PDB->Label);
      OutputPDB(ROTPDB,filename,remark);
      strcpy(remark,"");
     }
     if ( Outxyz==1 ) { 
       sprintf(filename,"%s_newrot%d.xyz",NewRotLabel,repeats*MaxNstructures+Nstructures+Offset);
       printf("Writing Outputxyz %s %s\n",filename,PDB->Label);
       OutputXYZ(ROTPDB,filename);
     }
     // RESET to copy of PDB after ligand move 
     FreeFRAGMENT(ROTPDB);
     ROTPDB=CopyFRAGMENT(LIGPDB); 
   }
   // RESET to original copy of PDB before ligand move 
   FreeFRAGMENT(LIGPDB); 
   LIGPDB=CopyFRAGMENT(PDB); 
   MaxNrotchange=MaxNrotchangeOrig;
   MaxNligchange=MaxNligchangeOrig;
   PDB->res[LigandResidue]->isfixed=0;
   LIGPDB->res[LigandResidue]->isfixed=0;
   repeats++;
  } 
}

void  ParseCommandLine ( const int argc, char *argv[]) {
extern FRAGMENTS * ReadPDB ( char *  ,int) ;
void ReadLIGDEF ( char * filename ) ; 

void help();
int narg;
int res; 
char filename[MaxLine];
    /** Parse Arguments */
    narg=1;
    while ( narg < argc ) { printf("%s ",argv[narg]); narg++;}
    printf("\n");

    narg=1;
    while ( narg < argc ) {
printf("Analysing command line option %s\n",argv[narg]); 
      if ( strncmp(argv[narg],"-startres",5) == 0 ) {
        if ( narg++ && narg < argc) { startres=atoi(argv[narg]); }
        else {  break ;  }
      }
      else if (strcmp ( argv[narg], "-hydrogens") == 0 ) {
         options.hydrogens=JIVE_TRUE;
      }
      else if (strncmp ( argv[narg], "-reducemax",10) == 0 ) {
          ReduceMaxRotChange=1;
      }
      else if (strncmp ( argv[narg], "-use_water",10) == 0 ) {
          globals.ignore_water = JIVE_FALSE ;
printf("reading in waters\n");
      }
      else if (strncmp ( argv[narg], "-noreducemax",12) == 0 ) {
          ReduceMaxRotChange=0;
      }
      else if (strncmp ( argv[narg], "-ligstates",10) == 0 ) {
        if ( narg++ && narg < argc) { LigRotamerStates=atoi(argv[narg]); }
        else {  break ;  }
        printf("LigRotamerStates=%d\n",LigRotamerStates);
      }
      else if (strncmp ( argv[narg], "-ligdouble",10) == 0 ) {
        if ( narg++ && narg < argc) {
         res=atoi(argv[narg])-1;
         LIGDEFdoubles[res]=1;
         if ( LIGdouble1==-1 ) LIGdouble1=res; 
         else if ( LIGdouble2==-1 ) LIGdouble2=res; 
         printf("group %d enabled for double ligand moves %d %d \n",res,LIGdouble1,LIGdouble2);
        }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-randgroup",10) == 0 ) {
          RandGroup=1;
      }
      else if (strncmp ( argv[narg], "-exhaustive",11) == 0 ) {
        printf("Using exhaustive search\n");
        ExhaustiveSearch=0;
      }
      else if (strncmp ( argv[narg], "-norandgroup",12) == 0 ) {
          RandGroup=0;
      }
      else if ( strncmp(argv[narg],"-maxrmsdtries",13) == 0 ) {
        if ( narg++ && narg < argc) { MaxRotamerTries=atoi(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-nrepeat",8) == 0 ) {
        if ( narg++ && narg < argc) { MaxRepeats=atoi(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-nstructures",12) == 0 ) {
        if ( narg++ && narg < argc) { MaxNstructures=atoi(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-endres",5) == 0 ) {
        if ( narg++ && narg < argc) { endres=atoi(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-all",4) == 0 ) {
        MinBFactor=-10000;
      }
      else if ( strncmp(argv[narg],"-outpdb",7) == 0 ) { 
        Outpdb = 1;
      }
      else if ( strncmp(argv[narg],"-outxyz",7) == 0 ) { 
        Outxyz = 1;
      }
      else if ( strncmp(argv[narg],"-pdb",4) == 0 ) { 
        if ( narg++ && narg < argc) { 
          PDB =ReadPDB(argv[narg],MaxRes); 
          if ( PDB == NULL ) { printf("############# ERROR NO PDB %s ###################\n",argv[narg]);exit;}
        }
      }
      else if ( strncmp(argv[narg],"-bbdscrwl",9) == 0 ) { 
        UsingLovell = 0;
      }
      else if ( strncmp(argv[narg],"-bbdlovell",10) == 0 ) { 
        UsingLovell = 1;
      }
      else if ( strncmp(argv[narg],"-bbdlib",7) == 0 ) { 
        if ( narg++ && narg < argc) {
         strcpy(globals.BbdDir,argv[narg]);
         printf("Setting location for BBD to %s\n",globals.BbdDir);
        }
      }
      else if ( strncmp(argv[narg],"-bbdminprob",11) == 0 ) { 
        if ( narg++ && narg < argc) {
          globals.MinRotamerProb=atof(argv[narg]);
        }
      }
      else if ( strncmp(argv[narg],"-file",5) == 0 ) { 
        if ( narg++ && narg < argc) {
         strcpy(globals.BbdDir,argv[narg]);
         printf("Setting location for BBD to %s\n",globals.BbdDir);
        }
      }
      else if ( strncmp(argv[narg],"-freeze",7) == 0 ) { 
        if ( narg++ && narg < argc) {
         res=atoi(argv[narg])-1;
         if ( PDB!=NULL && PDB->numres>res) PDB->res[res]->isfixed=2; 
         printf("Freezing %s %d\n",argv[narg],PDB->res[res]->isfixed);
        }
      }
      else if ( strncmp(argv[narg],"-markres",8) == 0 ) { 
        if ( narg++ && narg < argc) {
         res=atoi(argv[narg])-1;
         if ( PDB!=NULL && PDB->numres>res) PDB->res[res]->atom[CAATOM]->b_factor = MinBFactor+999 ; 
         printf("Selecting %s for rotamer search using B-factor \n",argv[narg]);
        }
      }
      else if ( strncmp(argv[narg],"-moveres",8) == 0 ) { 
        if ( narg++ && narg < argc) {
         res=atoi(argv[narg])-1;
         if ( res<MaxRes) {
           DeltaXYZ[res]=1; 
           printf("Selecting %d for moving \n",res);
         }
        }
      }
      else if ( strncmp(argv[narg],"-file",5) == 0 ) { 
        if ( narg++ && narg < argc) {
         strcpy(filename,argv[narg]);
        }
      }
      else if (strncmp ( argv[narg], "-tryharder",10) == 0 ) {
          TryHarder=1;
      }
      else if (strncmp ( argv[narg], "-notryharder",12) == 0 ) {
          TryHarder=0;
      }
      else if (strncmp ( argv[narg], "-tweaklig",9) == 0 ) {
        if ( narg++ && narg < argc) { Tweaklig =DegreeToRadian*atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-tweakchi",9) == 0 ) {
        if ( narg++ && narg < argc) { Tweakchi =DegreeToRadian*atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-tweakxyz",9) == 0 ) {
        if ( narg++ && narg < argc) { Tweakxyz =atof(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-limit",6) == 0 ) { 
        if ( narg++ && narg < argc) { limit = atoi(argv[narg]); }
      }
      else if ( strncmp(argv[narg],"-label",6) == 0 ) { 
        if ( narg++ && narg < argc) { strcpy(NewRotLabel,argv[narg]); }
      }
      else if ( strncmp(argv[narg],"-bfactor",8) == 0 ) { 
        if ( narg++ && narg < argc) { MinBFactor = atof(argv[narg]); }
      }
      else if ( strncmp(argv[narg],"-residuecentre",14) == 0 ) { 
        if ( narg++ && narg < argc) { CentreResidue = atoi(argv[narg])-1; }
      }
      else if ( strncmp(argv[narg],"-residueligand",14) == 0 ) { 
        if ( narg++ && narg < argc) { LigandResidue = atoi(argv[narg])-1; }
      }
      else if ( strncmp(argv[narg],"-rotreject",10) == 0 ) { 
        if ( narg++ && narg < argc) { AcceptDistanceThreshold  = (atof(argv[narg])); }
      }
      else if (strncmp ( argv[narg], "-rotamer",8) == 0 ) {
        if ( narg++ && narg < argc) { options.nrotamers= atoi(argv[narg]); }
        else {  break ;  }
        if ( options.nrotamers <0 || options.nrotamers>MaxRotamers ) options.nrotamers = MaxRotamers;
      }
      else if (strncmp ( argv[narg], "-maxligchange",13) == 0 ) {
        if ( narg++ && narg < argc) { MaxNligchange= atoi(argv[narg]); }
        else {  break ;  }
        if ( MaxNligchange <0)  MaxNligchange = 9999;
      }
      else if (strncmp ( argv[narg], "-maxrotchange",13) == 0 ) {
        if ( narg++ && narg < argc) { MaxNrotchange= atoi(argv[narg]); }
        else {  break ;  }
        if ( MaxNrotchange <0)  MaxNrotchange = 9999;
      }
      else if ( strncmp(argv[narg],"-maxrotcombine",14) == 0 ) {
          MaxNrotcombine=1;
      }
      else if ( strncmp(argv[narg],"-nomaxrotcombine",16) == 0 ) {
          MaxNrotcombine=0;
      }
      else if ( strncmp(argv[narg],"-vdwthresh",10) == 0 ) {
        if ( narg++ && narg < argc) { VDWThreshold=atof(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-rmsdthresh",11) == 0 ) {
        if ( narg++ && narg < argc) { RmsdThreshold=atof(argv[narg]); }
        else {  break ;  }
      }
      else if ( strncmp(argv[narg],"-distthresh",11) == 0 ) { 
        if ( narg++ && narg < argc) { DistanceThreshold  = (atof(argv[narg])); }
      }
      else if (strncmp ( argv[narg], "-rotthresh",10) == 0 ) {
        if ( narg++ && narg < argc) { RotThreshold= atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-chithresh",10) == 0 ) {
        if ( narg++ && narg < argc) { ChiThreshold= atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-xyzthresh",10) == 0 ) {
        if ( narg++ && narg < argc) { XYZThreshold= atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strncmp ( argv[narg], "-ligthresh",10) == 0 ) {
        if ( narg++ && narg < argc) { LIGThreshold= atof(argv[narg]); }
        else {  break ;  }
      }
      else if (strcmp ( argv[narg], "-ligdef") == 0 ) {
        if ( narg++ && narg < argc) { strcpy(LIGDEFfile,argv[narg]); }
        else {  break ;  }
      }

      else if ( strncmp(argv[narg],"-h",2) == 0 ) { help(); return ; }
      else if ( strncmp(argv[narg],"-v",2) == 0 ) { options.verbose++;}
      else if ( strncmp(argv[narg],"-V",2) == 0 ) { options.verbose++;}
      else if ( strncmp(argv[narg],"-offset",7) == 0 ) { 
        if ( narg++ && narg < argc) { Offset= atof(argv[narg]); }
        else {  break ;  }
      }
      narg++;
    }


}


/** Adjust Residue Sidechain to random rotamer **/
FRAGMENTS * AdjustSidechains ( FRAGMENTS * Structure, int startres, int endres)  {
extern void CloneFRAGMENT (FRAGMENTS * NEWFRAG,FRAGMENTS * OrigFrag) ;
extern float FindClash (FRAGMENTS * Frag, int res);
extern float MY_RANDOM_FLOAT ( float );
extern int MY_RANDOM_INT ( int );
extern  float CalcMinDistance (FRAGMENTS * Frag,int res1,int res2,int startatom, int endatom,int startatom2,int endatom2) ;
extern void BuildSidechains(FRAGMENTS *,int,int);
extern FRAGMENTS * CopyFRAGMENT(FRAGMENTS *);
extern void FreeFRAGMENT(FRAGMENTS *);
int AdjustLIGANDAngles (FRAGMENTS *) ; 
float CalcRotamerRMSD( FRAGMENTS * Pdb1,FRAGMENTS * Pdb2, int res ) ;
int  tweak_xyz_full (FRAGMENTS *  NewStructure,int res,float amount ) ;
int  tweak_xyz (FRAGMENTS *  NewStructure,int res ) ;
int  tweak_chi (FRAGMENTS *  NewStructure,int res ) ;

float rand;
float orig_distance;
float lig_distance;
float new_distance;
float vdw_distance;
float chi1,chi2,chi3,chi4;
float rmsd;
int res; 
float mindist; 
int minres; 
int group;
int number_selected=0;
int Altered=0;
int MadeLigandChange=0;
FRAGMENTS * OldStructure; 

       if ( startres > Structure->numres) startres = Structure->numres; 
       if ( endres   > Structure->numres) endres   = Structure->numres; 
       if ( startres < 0 )  startres=0;
       if ( endres   < 0 )  endres=Structure->numres;

       OldStructure=CopyFRAGMENT(Structure);
       printf("Creating list for sidechain rotamer moves\n"); 
       //Create a list of selected residues 
       for ( res=startres;res<endres;res++) {
          SelectedResidues[res]=-1; 
          if ( Structure->res[res]->isfixed==2 ) { printf("Found Freeze Residue %d \n",res); continue; }
          // Check if selected by bfactor 
          if ( DistanceThreshold < 9999 ) { 
           // Check how close to centre
           orig_distance=DistanceThreshold+1;
           lig_distance=DistanceThreshold+1;
           if ( CentreResidue >= 0 && CentreResidue<Structure->numres) orig_distance = (CalcMinDistance(Structure,CentreResidue,res,ATOM_CA,999,ATOM_CA,999));
           //// THIS DOESN'T WORK  IF LigRes is set ( even if no moves for lig ) then this screws up - instead of increasing the number of residues to be moved it shrinks them. 
           ////if ( LigandResidue >= 0 && LigandResidue<Structure->numres) lig_distance  = (CalcMinDistance(Structure,LigandResidue,res,ATOM_CA,999,ATOM_CA,999));
           //printf("Distance of %d to %d is %f threshold %f\n",res,CentreResidue, (orig_distance),(DistanceThreshold ));
           if ( orig_distance< DistanceThreshold && Structure->res[res]->Sequence != 'G' && Structure->res[res]->Sequence != 'A') { 
              SelectedResidues[res]=orig_distance; 
              number_selected++;
              printf("Marking res %d %s %s :%c: as possible residue to choose - distance to centre %f total=%d\n",res,Structure->res[res]->Resname,Structure->res[res]->Seqno,Structure->res[res]->Sequence,orig_distance,number_selected); 
           } 
           // Check how close to Ligand  - this may screw up scaling by distance later 
           else if ( lig_distance< DistanceThreshold && Structure->res[res]->Sequence!= 'G' && Structure->res[res]->Sequence != 'A') { 
              SelectedResidues[res]=lig_distance; 
              number_selected++;
              printf("Marking LIG res %d %s %s :%c: as possible residue to choose - distance to ligand %f total=%d\n",res,Structure->res[res]->Resname,Structure->res[res]->Seqno,Structure->res[res]->Sequence,orig_distance,number_selected); 
           } 
          } 
          if ( Structure->res[res]->atom[CAATOM]->b_factor >= MinBFactor ) { 
             SelectedResidues[res]=0.0001;
             number_selected++;
             printf("FOUND Marked residue  %d\n",res); 
          }
          if ( DeltaXYZ[res] > -2 ) { 
             SelectedResidues[res]=0.0001;
             number_selected++;
             printf("FOUND Moving residue  %d\n",res); 
          }
       } 

       Altered=0; 
       // Sort residues by the distance/b_factor  -- This is very clunky 
       // Loop through these residues based on distance 
       minres=9999;
       while ( MaxNrotchange>=0 && minres>-1 && number_selected > 0 ) { 
         mindist = 9999; minres  = -1; 
         // Rather than sequentially choosing groups, choose randomnly as long as distance is within threshold ( SelectedResidues[res]>0)
         if ( RandGroup == 1 ) { 
             res=MY_RANDOM_INT(endres-1);
             printf("Selecting res %d as possible residue to choose - distance %f \n",res,SelectedResidues[res]); 
             while ( SelectedResidues[res]<0 ) {
                  res=MY_RANDOM_INT(endres-1);
                  printf("Selecting res %d as possible residue to choose - distance %f \n",res,SelectedResidues[res]); 
             }
             mindist=SelectedResidues[res]; minres=res; 
         }
         else { 
          // Find next closest residue 
          res=startres; 
          while ( res<endres ) { 
           //printf("looking for next min dist %d %f current min %d %f\n",res,SelectedResidues[res],minres,mindist);
           if (SelectedResidues[res]>=0 && SelectedResidues[res]< mindist ) { 
             mindist=SelectedResidues[res]; minres=res; 
           }      
           res++;
          }
         }      
         printf("#####################################################################\n");
         printf("### Selected residue %d distance %f MaxNrotchange %d number left within threshold  %d\n",minres,SelectedResidues[minres],MaxNrotchange,number_selected);
         printf("#####################################################################\n");
         if ( mindist == 9999 ) break;
         orig_distance = SelectedResidues[minres];
         if ( orig_distance == 0 ) orig_distance=1;
         res = minres; 
         RotamerTries=0;
// Alter sidechain rotamer 
         if ( MaxNrotchange > 0 ) {
           //printf("Residue %d is within the distance or b-factor threshold residue left to change %d\n",res, MaxNrotchange);
// For residues which are close, rand decreases ie more likely to be picked (by  sqrt(DistanceThreshold))
// For residues which are far away,  rand increases ie less likely to be picked 
           rand = MY_RANDOM_FLOAT(1);
           if ( DistanceThreshold < 9999 ) rand=rand/sqrt(DistanceThreshold/orig_distance);
           printf("Considering res %d for dRotamer.threshold=%f rand=%f == ran(%f)/ratio(%f)\n",res,RotThreshold,rand,rand*sqrt(DistanceThreshold/orig_distance),sqrt(DistanceThreshold/orig_distance));
           if ( rand<RotThreshold || ExhaustiveSearch>=0 ) { 
            while ( RotamerTries<MaxRotamerTries) { 
             // Save old values 
             //chi1=Structure->res[res]->chi1; chi2=Structure->res[res]->chi2; 
             //chi3=Structure->res[res]->chi3; chi4=Structure->res[res]->chi4; 
             // Reset to undefined so function will change it
             //FreeFRAGMENT(OldStructure); 
             // keep a copy of original struct so we can reject move
             CloneFRAGMENT (OldStructure,Structure) ;
/////////////////////////////////////////
// Make the moves (XYZ/Ligand/Rotamer) 
/////////////////////////////////////////
printf("MaxNligchange %d ExhaustiveSearch %d res %d LigandResidue %d DeltaXYZ[res] %f\n",MaxNligchangeOrig,ExhaustiveSearch,res,LigandResidue,DeltaXYZ[res]);
             MadeLigChange=0;
             // Exhaustive XYZ moves are 1*Tweakxyz, then -1*Tweakxyz
             if ( ExhaustiveSearch>=0 && Tweakxyz>0 && DeltaXYZ[res] > -2 ) { 
                 printf("Considering res %d for discrete dXYZ. rand=%f \n",res,Tweakxyz);
                 tweak_xyz_full (Structure,res,(Tweakxyz*DeltaXYZ[res]));
                 sprintf(remark,"%sREMARK DiscreteXYZ Adjusting res %d %f \n",remark,res,(Tweakxyz*DeltaXYZ[res]));
                 // Also tweak neighbouring residues by half the amount to stop crazy stuff
                 if ( res>0) { 
                      tweak_xyz_full (Structure,res-1,(Tweakxyz*DeltaXYZ[res])/2);
                      sprintf(remark,"%sREMARK DiscreteXYZ Adjusting res %d %f \n",remark,res-1,(Tweakxyz*DeltaXYZ[res])/2);
                 }
                 if ( res+1 < Structure->numres ) {
                      tweak_xyz_full (Structure,res+1,(Tweakxyz*DeltaXYZ[res])/2);
                      sprintf(remark,"%sREMARK DiscreteXYZ Adjusting res %d %f \n",remark,res+1,(Tweakxyz*DeltaXYZ[res])/2);
                 }

                 DeltaXYZ[res]=DeltaXYZ[res]-1;
                 if ( DeltaXYZ[res]==0 ) DeltaXYZ[res]=DeltaXYZ[res]-1;
                 MadeLigChange=1;
             }
             else if ( ExhaustiveSearch>=0 && MaxNligchangeOrig > 0 && res==LigandResidue ) { 
                  printf("Considering res %d for discrete LIGAND moves. \n",res);
                  MadeLigChange=AdjustLIGANDAngles  ( Structure ) ; 
                  // If all single ligand moves have been done then try doubles 
                  if ( MadeLigChange==0 && LIGdouble1>-1 && LIGdouble2>-1 ) { 
                    MadeLigChange=AdjustLIGANDAnglesSubset ( LIGPDB ) ;    // Allow double moves for a subset of the groups for exhaustive search
                  }
             }
             // If no ligand or BB movement ( or current residue has no rotamers left ) then tweak residue  
             if ( MadeLigChange==0 && res!= LigandResidue ) { 
               printf("Considering res %d for Rotamer moves. \n",res);
               printf("rand within threshold -Tweaking rotamer for residue %d max_rot %d current chi1 %f chi2 %f chi3 %f chi4 %f \n",res,options.nrotamers, Structure->res[res]->chi1*RadianToDegree,Structure->res[res]->chi2*RadianToDegree,Structure->res[res]->chi3*RadianToDegree,Structure->res[res]->chi4*RadianToDegree); 
               Structure->res[res]->chi1=UNDEFINED_ANGLE; Structure->res[res]->chi2=UNDEFINED_ANGLE; 
               Structure->res[res]->chi3=UNDEFINED_ANGLE; Structure->res[res]->chi4=UNDEFINED_ANGLE; 
               printf("res res=%d SelectedResidues[minres=%d]=%f\n",res,minres,SelectedResidues[minres]);
               printf("RES %d max_rot %d before AddDunbrack %d\n",res,options.nrotamers,ExhaustiveSearch);
               Structure=AddDunbrackSidechains(Structure,res,res+1,options.nrotamers,&ExhaustiveSearch);
               MadeLigChange=1; 
               if ( Structure->res[res]->Resindex == RES_GLY || Structure->res[res]->Resindex == RES_ALA ) { ExhaustiveSearch=-1;} 
               ExhaustiveSearch++;
               if ( ExhaustiveSearch==0 ) { MadeLigChange=0;} 
               printf("RES %d max_rot %d after AddDunbrack %d\n",res,options.nrotamers,ExhaustiveSearch);
               if ( MadeLigChange==1) sprintf(remark,"%sREMARK Tweaking rotamer res %d chi1 %f chi2 %f chi3 %f chi4 %f \n",remark,res,options.nrotamers, Structure->res[res]->chi1*RadianToDegree,Structure->res[res]->chi2*RadianToDegree,Structure->res[res]->chi3*RadianToDegree,Structure->res[res]->chi4*RadianToDegree); 
             }
/////////////////////////////////////
// Moves exhausted? 
/////////////////////////////////////
             printf("Finished moves - ExhaustiveSearch=%d DeltaXYZ[res]=%f MadeLigChange=%d\n",ExhaustiveSearch,DeltaXYZ[res],MadeLigChange );
             if ( MadeLigChange==0  ) { 
                // Once residue has exhausted its moves - Have to freeze residue in global variable pdb structure !!
                LIGPDB->res[res]->isfixed=2 ;
                ROTPDB->res[res]->isfixed=2 ;
                PDB->res[res]->isfixed=2 ;
                printf("Exhausted res res=%d - Now freezing %s %s \n",res,ROTPDB->res[res]->Resname,PDB->res[res]->Seqno);
                break; // Number of rotamers exceeded. 
             }
/////////////////////
// Build new rotamer 
/////////////////////
             vdw_distance=VDWThreshold; 
             rmsd=RmsdThreshold;
             new_distance=-1;  
             if ( res != LigandResidue ) { 
              Structure->res[res]->isfixed=0; 
              printf("Assigned new rotamer for res %d Rot %d nrot %d chi1 %f chi2 %f chi3 %f chi4 %f \n",res,ExhaustiveSearch-1,options.nrotamers, Structure->res[res]->chi1*RadianToDegree,Structure->res[res]->chi2*RadianToDegree,Structure->res[res]->chi3*RadianToDegree,Structure->res[res]->chi4*RadianToDegree); 
              (void) BuildSidechains(Structure,res,res+1);
              rmsd=CalcRotamerRMSD(OldStructure,Structure,  res ) ; 
              new_distance = (CalcMinDistance(Structure,CentreResidue,res,ATOM_CA,999,ATOM_CA,999));
              vdw_distance=FindClash (Structure, res) ;
              printf("Res %d Minimum VDW distance for %f\n", res,vdw_distance);
              printf("Res %d Minimum distance to centre BEFORE rot %f AFTER rot %f DIFF %f \n", res,orig_distance,new_distance, orig_distance-new_distance ) ; 
             }
////////////////////////////
// Check Acceptance criteria
////////////////////////////
             // Need to add Chi threshold check for non-ligand residue 
             // RMSD Threshold
             if ( rmsd < RmsdThreshold) { 
               // Reset to old values - should just use CopyFragmentInfo ?
               printf("Rejecting new rotamer: rmsd too small %f thresh %f tries %d\n",rmsd,RmsdThreshold,RotamerTries);
               CloneFRAGMENT (Structure,OldStructure) ;
               RotamerTries++; 
               //if ( ExhaustiveSearch>0 ) ExhaustiveSearch++;
               continue;
             }  
             // We only accept this move if it doesn't screw up a close contact
             if ( VDWThreshold>0 && vdw_distance < VDWThreshold ) { 
               printf("Rejecting new rotamer: vdw overlap too small vdwdist %f\n",vdw_distance);
                CloneFRAGMENT (Structure,OldStructure) ;
                RotamerTries++; 
                //if ( ExhaustiveSearch>0 ) ExhaustiveSearch++;
                continue;
             }
             else if ( orig_distance <= AcceptDistanceThreshold && new_distance > AcceptDistanceThreshold) { 
                // Reset to old values 
                printf("Rejecting new rotamer: screws up good contact dist %f dist2 %f\n",orig_distance,new_distance);
                CloneFRAGMENT (Structure,OldStructure) ;
                printf("Reset distance is %f\n", (CalcMinDistance(Structure,CentreResidue,res,ATOM_CA,999,ATOM_CA,999)));
                RotamerTries++; 
                //if ( ExhaustiveSearch>0 ) ExhaustiveSearch++;
                continue;
             }
             else if ( TryHarder == 1 && new_distance > AcceptDistanceThreshold && RotamerTries < MaxRotamerTries-1 ) { 
                printf("Rejecting new rotamer: dist to centre too large -trying harder - dist %f dist2 %f\n",orig_distance,new_distance);
                CloneFRAGMENT (Structure,OldStructure) ;
                printf("Reseted distance is %f\n", (CalcMinDistance(Structure,CentreResidue,res,ATOM_CA,999,ATOM_CA,999)));
                RotamerTries++; 
                //if ( ExhaustiveSearch>0 ) ExhaustiveSearch++;
                continue;
             }
             else {
               printf("Accepting new rotamer %f %f \n",orig_distance,new_distance);
               Structure->res[res]->isfixed=0; 
               CloneFRAGMENT (OldStructure,Structure) ;
               MaxNrotchange--;
               number_selected--;
               Altered=1; 
               // If selecting by random order, consider another residue 
               SelectedResidues[minres]=-1; 
               break;
             }
           } // reached end of number of re-tries
          }// end of change rotamer 
         } // reached max number of rotamer changes
         //  Tweak Sidechain ANGLE 
         if ( Tweakchi > 0) { 
             rand = MY_RANDOM_FLOAT(1);
             if ( DistanceThreshold < 9999 ) rand=rand/sqrt(DistanceThreshold/orig_distance);
             printf("Considering res %d for dChi. threshold=%f rand=%f == ran(%f)/ratio(%f)\n",res,ChiThreshold,rand,rand*sqrt(DistanceThreshold/orig_distance),sqrt(DistanceThreshold/orig_distance));
             if ( rand< ChiThreshold ) { 
               printf("rand within threshold for tweaking chi res %d\n",res);
               tweak_chi (Structure,res);
               Structure->res[res]->isfixed=0; 
               (void) BuildSidechains(PDB,res,res+1);
               Altered=1; 
             }
         }
         //  Tweak XYZ coords 
         if ( ExhaustiveSearch < 0 && Tweakxyz > 0 ) { 
            rand = MY_RANDOM_FLOAT(1);
            if ( DistanceThreshold < 9999 ) rand=rand/sqrt(DistanceThreshold/orig_distance);
            printf("Considering res %d for dXYZ. threshold=%f rand=%f == ran(%f)/ratio(%f)\n",res,XYZThreshold,rand,rand*sqrt(DistanceThreshold/orig_distance),sqrt(DistanceThreshold/orig_distance));
            if ( rand<XYZThreshold ) { 
              printf("rand with threshold for tweaking xyz res %d \n",res);
              tweak_xyz (Structure,res);
              Altered=1; 
            }
         }
         // If selecting by distance order, consider another residue 
         if (RandGroup == 0 ) SelectedResidues[minres]=-1; 
         // If exhausted all rotamers in res consider another residue 
         if (ExhaustiveSearch>=0 ) SelectedResidues[minres]=-1; 
       }
       if ( Altered == 0 ) return NULL; 
       else { return Structure; }
}

// Tweak the xyz position of residue res by random amount
int  tweak_xyz (FRAGMENTS *  NewStructure,int res ) {
extern float MY_RANDOM_FLOAT ( float );
RESIDUE * ThisRes ;
int natom;
float rand;

       ThisRes =  NewStructure->res[res];
       for ( natom = 0 ; natom < ThisRes->numatom;natom++ ) {
          rand = MY_RANDOM_FLOAT(2)-1.0;
          ThisRes->atom[natom]->x += Tweakxyz * rand ;
          printf("Altering res %d atom %d ",res,natom);
          printf("x by %f ",Tweakxyz * rand);
          rand = MY_RANDOM_FLOAT(2)-1.0;
          ThisRes->atom[natom]->y += Tweakxyz * rand ;
          printf("y by %f ",Tweakxyz * rand);
          rand = MY_RANDOM_FLOAT(2)-1.0;
          ThisRes->atom[natom]->z += Tweakxyz * rand ;
          printf("z by %f\n",Tweakxyz * rand);
       }
}
// Tweak the xyz position of residue res by discrete amount
int  tweak_xyz_full (FRAGMENTS *  NewStructure,int res,float amount ) {
RESIDUE * ThisRes ;
int natom;

       ThisRes =  NewStructure->res[res];
       for ( natom = 0 ; natom < ThisRes->numatom;natom++ ) {
          ThisRes->atom[natom]->x += amount  ;
          printf("Altering res %d atom %d ",res,natom);
          printf("x by %f ",amount);
          ThisRes->atom[natom]->y += amount ;
          printf("y by %f ",amount );
          ThisRes->atom[natom]->z += amount ;
          printf("z by %f\n",amount );
       }
}

// Tweak the chi angles of residue res
int  tweak_chi (FRAGMENTS *  NewStructure,int res ) {
extern float MY_RANDOM_FLOAT ( float );
RESIDUE * ThisRes ;
float rand;
//int natom;
//float chi_orig;

       ThisRes =  NewStructure->res[res];
       //for ( natom = 0 ; natom < ThisRes->numatom;natom++ ) {
        // ThisRes->atom[natom]-> occupancy = 1.00 ;
       //}
       if ( ThisRes->chi1 != UNDEFINED_ANGLE )  {
          rand = MY_RANDOM_FLOAT(2)-1.0;
          ThisRes->chi1 += Tweakchi * rand ;
          if (options.verbose>0) printf("tweaking chi1 rand = %f ",rand* Tweakchi*RadianToDegree);
       }
       if ( ThisRes->chi2 != UNDEFINED_ANGLE )  {
         rand = MY_RANDOM_FLOAT(2)-1.0;
         ThisRes->chi2 += Tweakchi * rand ;
         if (options.verbose>0) printf("tweaking chi2 rand = %f ",rand* Tweakchi*RadianToDegree);
       }
       if ( ThisRes->chi3 != UNDEFINED_ANGLE )  {
         rand = MY_RANDOM_FLOAT(2)-1.0;
         ThisRes->chi3 += Tweakchi * rand ;
         if (options.verbose>0) printf("tweaking chi3 rand = %f ",rand* Tweakchi*RadianToDegree);
       }
       if ( ThisRes->chi4 != UNDEFINED_ANGLE )  {
         rand = MY_RANDOM_FLOAT(2)-1.0;
         ThisRes->chi4 += Tweakchi * rand ;
         if (options.verbose>0) printf("tweaking chi4 rand = %f ",rand* Tweakchi*RadianToDegree);
       }
       printf("\n");
   return res;
}


void help ( ) {
       printf("A program which changes the sidechain conformation using a rotamer library\n");
       printf("Usage: PdbRotamerSearch [options]\n");
       printf("VERSION 2.0 July 2011 \n "); 
       printf("Now added Exhaustive searching\n");
       printf("[options]\n");
       printf ("-use_water          : Treat water molecules as crystallographic waters ie read them NEEDS TO BEFORE -pdb flag  \n");
       printf ("-pdb pdbfile        : pdbfile \n");
       printf ("-label string       : Filenames has the name ${string}_newrot{n}.pdb\n");
       printf ("                    : Default ${pdbfile:r}_newrot{n}.pdb \n");
       printf ("=======================\n");
       printf ("Selection of residues to change \n");
       printf ("=======================\n");
       printf ("-all                : Change sidechains for all residues ( same as -bfactor -1000)\n");
       printf ("-startres int       : Change sidechain starting at residue {int} \n");
       printf ("-endres int         : Change sidechain ending   at residue {int} \n");
       printf ("-bfactor       {num} : Change residues with b_Factor > num  default 9999 \n");
       printf ("-freeze       {res}  : DO NOT move this residue\n");
       printf ("-markres      {res}  : DO rotamer search this residue (uses bfactor - put after any -bfactor changes)\n");
       printf ("-moveres      {res}  : DO Move XYZ for this residue ( needs -tweakxyz {n})\n");
       printf ("-residuecentre {res}*** : Select residue res as centroid - default undefined \n");
       printf ("-residueligand {res}*** : Select residue res as ligand - default same as residuecentre  \n");
       printf ("NOTE:\n");
       printf ("The probability of selecting a residue for a move is increased by sqrt(dist_thresh/distance)\n");
       printf ("ie For residues close to centre eg (distance~2) increase is sqrt(10/2) ie sqrt(5)= ~2\n"); 
       printf ("   For residues at the distance limit, the probability is unchanged ie ratio==1\n");
       printf ("=======================\n");
       printf ("==== LIGAND MOVES =====\n");
       printf ("=======================\n");
       printf ("-ligdef file        : Define LIG groups to rotate \n");
       printf ("-ligdoubles   {group} : DO Allow double ligand moves for this group \n");
       printf ("=======================\n");
       printf ("Thresholds \n");
       printf ("=======================\n");
       printf ("-distthreshold  dist  : Change sidechains only for residues with distance of < num from residue {res} \n");
       printf ("-vdwthreshold  thresh : Only accept rotamer if there are no vdw contacts < thresh \n");
       printf ("-rmsdthreshold thresh : Only accept rotamer if rmsd is greater than this value\n");
       printf ("-rotthreshold thresh  : Only change rotamer if rand() < threshold (<1)\n");
       printf ("-chithreshold thresh  : Only change chi angle if rand() < threshold (<1)\n");
       printf ("-xyzthreshold thresh  : Only change xyz angle if rand() < threshold (<1) \n");
       printf ("-ligthreshold thresh  : Only change LIG angle if rand() < threshold (<1)\n");
       printf ("-bbdminprob   thresh  : Only use rotamer is observed freq in sclib is higher than thresh ( 0.01 ie 1%)\n");
       printf ("=======================\n");
       printf ("Limits\n");
       printf ("=======================\n");
       printf ("-maxligchange N     : Only change at most, N groups within ligand/residuecentre \n");
       printf ("-maxrotchange N     : Only change at most, N amino acids \n");
       printf ("-reducemax          : For each structure, scale maxligchange and maxrotchange by rand(1)\n");
       printf ("-randgroup          : Randomnly choose lig group/residue rather than cycling through them\n");
       printf ("-maxrotcombine      : Only change at most, -maxrotchange N groups OR amino acids combined \n");
       printf ("-rotamer N          : Allow top N likely rotamers (max 20) for SCRWL library\n");
       printf ("-nrepeat  N         : Produce N searches \n ");
       printf ("-nstructures n      : For every ligand search, produce n rotamer searches \n ");
       printf ("                    : The total number of structures produced is Nxn      \n ");
       printf ("-maxrmsdtries    n  : Try n attempts to change rotamer and then give up\n");
       printf ("-rotreject  distance: Reject rotamer if orig is within {distance} from lig and new rotamer is further than {distance} ie 4.2\n");
       printf ("-tryharder          : Only accept rotamer if within distance threshold or number of tries reached \n");
       printf ("-offset      n      : Start from n in the output filename ie -offset 100 starts from file_newrot101.pdb \n");
       printf ("=======================\n");
       printf ("Output \n");
       printf ("=======================\n");
       printf ("-outxyz             : Output XYZ format file \n");
       printf ("-outpdb             : Output PDB format file \n");
       printf ("-hydrogens    {num} : Add hydrogens NOT WORKING \n");
       printf ("=======================\n");
       printf ("Moves \n");
       printf ("=======================\n");
       printf ("-tweakchi     {num} : Adjust sidechain dihedral angles of residue by up to +-{num} degrees \n");
       printf ("-tweakxyz     {num} : Adjust xyz coords of residue by up to +-{num} Angstroms \n");
       printf ("-tweaklig     {num} : Adjust dihedral angles of LIGAND by up to +-{num} degrees \n");
       printf ("\n");
       printf ("=======================\n");
       printf ("Sidechain Libraries: \n");
       printf ("========================\n");
       printf ("-bbdscwrl           : Use Dunbrack SCRWL library format \n");
       printf ("-bbdlovell          : Use Lovell penultimate library format\n");
       printf ("-bbdlib       {num} : Location of Lovell or Dunbrack BackBone Dependant library\n");
       printf ("SLOOP_BBD       # BBDEP library file for sidechains %s\n",globals.BbdDir);


       printf ("\n");
} 









// THIS CAN BE SPEEDED UP BY USING MATRIX MULTIPLICATION
// Rotate LIG ring about a point 
void TransformLIG (RESIDUE * ThisRes, int group, float Angle) {

float Xcentre, Ycentre, Zcentre;
float Xcentre2, Ycentre2, Zcentre2;
float Xcentre3, Ycentre3, Zcentre3;
float rotation00, rotation01, rotation02, rotation10, rotation11, rotation12, rotation20, rotation21, rotation22;
float zrotaxis00, zrotaxis01, zrotaxis02, zrotaxis10, zrotaxis11, zrotaxis12, zrotaxis20, zrotaxis21, zrotaxis22;
float yrotaxis00, yrotaxis01, yrotaxis02, yrotaxis10, yrotaxis11, yrotaxis12, yrotaxis20, yrotaxis21, yrotaxis22;
float cosA,sinA;
float hypsqr;
int Atom; 
int index;
int firstatom ;
int secondatom ;
ATOM * ThisAtom;
 float x,y,z;
 float x2,y2,z2;

//printf("Transforming LIG group %d by %f\n",group,Angle*RadianToDegree);

// Define centres
   firstatom=LIGDEFtable[group][0];
   secondatom=LIGDEFtable[group][1];
   if ( firstatom>=ThisRes->numatom || secondatom>=ThisRes->numatom) return; 
// Blank Rotation about Z axis
   zrotaxis00=1; zrotaxis01=0; zrotaxis02=0; 
   zrotaxis10=0; zrotaxis11=1; zrotaxis12=0; 
   zrotaxis20=0; zrotaxis21=0; zrotaxis22=1; 

//////////////////////////////////////////////////
// rotate in Z so that second atom is aligned to X axis 
//////////////////////////////////////////////////
   Xcentre=ThisRes->atom[firstatom]->x; Ycentre=ThisRes->atom[firstatom]->y; Zcentre=ThisRes->atom[firstatom]->z; 
   Xcentre2=ThisRes->atom[secondatom]->x; Ycentre2=ThisRes->atom[secondatom]->y; Zcentre2=ThisRes->atom[secondatom]->z; 
   hypsqr = Sqr(ThisRes->atom[secondatom]->x-Xcentre)+Sqr(ThisRes->atom[secondatom]->y-Ycentre);
   if (hypsqr > 0 ) { 
     cosA = (ThisRes->atom[secondatom]->x-Xcentre)/sqrt(hypsqr); 
     sinA = (ThisRes->atom[secondatom]->y-Ycentre)/sqrt(hypsqr); 
//printf("cos angle is %f %f \n",cosA,acos(cosA)*RadianToDegree);
//printf("sin angle is %f %f \n",sinA,asin(sinA)*RadianToDegree);
     // Rotate about second atom in z plane to give zero in xaxis
     zrotaxis00=cosA;    zrotaxis01=sinA; zrotaxis02=0; 
     zrotaxis10=-sinA; zrotaxis11=cosA; zrotaxis12=0; 
     zrotaxis20=0;       zrotaxis21=0;    zrotaxis22=1; 
   }
   //printf("zrot matrix is %f %f %f\n",zrotaxis00, zrotaxis01,zrotaxis02);
   //printf("zrot matrix is %f %f %f\n",zrotaxis10, zrotaxis11,zrotaxis12);
   //printf("zrot matrix is %f %f %f\n",zrotaxis20, zrotaxis21,zrotaxis22);
   secondatom=LIGDEFtable[group][1];
//printf("    coords for first  atom %f %f %f\n",Xcentre,Ycentre,Zcentre);
//printf("ORG coords for second atom %f %f %f\n",Xcentre2,Ycentre2,Zcentre2);
//printf("CNTDcoords for second atom %f %f %f\n",Xcentre2-Xcentre,Ycentre2-Ycentre,Zcentre2-Zcentre);

   index=1; Atom=LIGDEFtable[group][index];
   while ( Atom >= 0 ) { 
    ThisAtom=ThisRes->atom[Atom];
    /* centre molecule about second atom  */
    if (ThisAtom != NULL ) {
      // centre on first atom 
      x= ThisAtom->x-Xcentre ; y= ThisAtom->y -Ycentre ; z= ThisAtom->z -Zcentre ;
      // Rotate about x axis to centre second atom on axis
      ThisRes->atom[Atom]->x  = x*zrotaxis00 + y*zrotaxis01 + z*zrotaxis02 ;
      ThisRes->atom[Atom]->y  = x*zrotaxis10 + y*zrotaxis11 + z*zrotaxis12 ;
      ThisRes->atom[Atom]->z  = x*zrotaxis20 + y*zrotaxis21 + z*zrotaxis22 ;
    }
    index++; Atom=LIGDEFtable[group][index];
   }
   secondatom=LIGDEFtable[group][1];
   Xcentre2=ThisRes->atom[secondatom]->x; Ycentre2=ThisRes->atom[secondatom]->y; Zcentre2=ThisRes->atom[secondatom]->z; 
//printf("new coords for second atom %f %f %f\n",Xcentre2,Ycentre2,Zcentre2);

//////////////////////////////////////////////////
// rotate in Y so that second atom is aligned to X axis 
//////////////////////////////////////////////////
// Blank Rotation about Y axis
   yrotaxis00=1; yrotaxis01=0; yrotaxis02=0; 
   yrotaxis10=0; yrotaxis11=1; yrotaxis12=0; 
   yrotaxis20=0; yrotaxis21=0; yrotaxis22=1; 

   hypsqr = Sqr(ThisRes->atom[secondatom]->z)+Sqr(ThisRes->atom[secondatom]->x);
   if (hypsqr > 0 ) { 
     cosA = (ThisRes->atom[secondatom]->x)/sqrt(hypsqr); 
     sinA = (ThisRes->atom[secondatom]->z)/sqrt(hypsqr); 
//printf("cos angle is %f %f \n",cosA,acos(cosA)*RadianToDegree);
//printf("sin angle is %f %f \n",sinA,asin(sinA)*RadianToDegree);
     // Rotate about second atom in Y plane to give zero in xaxis
     yrotaxis00=cosA;    yrotaxis01=0; yrotaxis02=sinA; 
     yrotaxis10=0;       yrotaxis11=1; yrotaxis12=0; 
     yrotaxis20=-sinA; yrotaxis21=0; yrotaxis22=cosA; 
   }
   //printf("yrot matrix is %f %f %f\n",yrotaxis00, yrotaxis01,yrotaxis02);
   //printf("yrot matrix is %f %f %f\n",yrotaxis10, yrotaxis11,yrotaxis12);
   //printf("yrot matrix is %f %f %f\n",yrotaxis20, yrotaxis21,yrotaxis22);
   secondatom=LIGDEFtable[group][1];

   index=1; Atom=LIGDEFtable[group][index];
   while ( Atom >= 0 ) { 
    ThisAtom=ThisRes->atom[Atom];
    /* centre molecule about second atom  */
    if (ThisAtom != NULL ) {
      // centre on first atom 
      x= ThisAtom->x ; y= ThisAtom->y  ; z= ThisAtom->z ;
      // Rotate about Y axis to centre second atom on axis
      ThisRes->atom[Atom]->x  = x*yrotaxis00 + y*yrotaxis01 + z*yrotaxis02 ;
      ThisRes->atom[Atom]->y  = x*yrotaxis10 + y*yrotaxis11 + z*yrotaxis12 ;
      ThisRes->atom[Atom]->z  = x*yrotaxis20 + y*yrotaxis21 + z*yrotaxis22 ;
    }
    index++; Atom=LIGDEFtable[group][index];
   }
   Xcentre3=ThisRes->atom[secondatom]->x; Ycentre3=ThisRes->atom[secondatom]->y; Zcentre3=ThisRes->atom[secondatom]->z; 
//printf("ORG coords for second atom %f %f %f\n",Xcentre2,Ycentre2,Zcentre2);
//printf("CNTDcoords for second atom %f %f %f\n",Xcentre3,Ycentre3,Zcentre3);

//////////////////////////////////////////////////////////////
// Rotation about xaxis, centred on second atom by angle Angle
//////////////////////////////////////////////////
   rotation00=1; rotation01=0;          rotation02=0; 
   rotation10=0; rotation11=cos(Angle); rotation12=-sin(Angle); 
   rotation20=0; rotation21=sin(Angle); rotation22=cos(Angle); 
//   printf("xrot matrix is %f %f %f\n",rotation00, rotation01,rotation02);
//   printf("xrot matrix is %f %f %f\n",rotation10, rotation11,rotation12);
//   printf("xrot matrix is %f %f %f\n",rotation20, rotation21,rotation22);

//printf("XCentre1                   %f %f %f\n",Xcentre,Ycentre,Zcentre);
//printf("XCentre2                   %f %f %f\n",Xcentre2,Ycentre2,Zcentre2);
//printf("XCentre3                   %f %f %f\n",Xcentre3,Ycentre3,Zcentre3);
// Rotate about x axis and uncentre and unrotate
   index=1; Atom=LIGDEFtable[group][index];
   while ( Atom >= 0 ) { 
    ThisAtom=ThisRes->atom[Atom];
    /* centre molecule about second atom  */
    if (ThisAtom != NULL ) {
      x= ThisAtom->x-Xcentre3 ; y= ThisAtom->y -Ycentre3 ; z= ThisAtom->z -Zcentre3 ;
      // Rotate Angle about x axis and uncentre on second atom
      x2 = x*rotation00 + y*rotation01 + z*rotation02 + Xcentre3;
      y2 = x*rotation10 + y*rotation11 + z*rotation12 + Ycentre3;
      z2 = x*rotation20 + y*rotation21 + z*rotation22 + Zcentre3;

      // UnRotate about y axis and re-centre on second atom 
      x= x2*yrotaxis00 + y2*yrotaxis01 - z2*yrotaxis02 ;
      y= x2*yrotaxis10 + y2*yrotaxis11 + z2*yrotaxis12 ;
      z= -x2*yrotaxis20 + y2*yrotaxis21 + z2*yrotaxis22;

      // UnRotate about z axis and re-centre on first  atom 
      ThisRes->atom[Atom]->x = x*zrotaxis00 - y*zrotaxis01 + z*zrotaxis02 + Xcentre;
      ThisRes->atom[Atom]->y = -x*zrotaxis10 + y*zrotaxis11 + z*zrotaxis12 + Ycentre;
      ThisRes->atom[Atom]->z = x*zrotaxis20 + y*zrotaxis21 + z*zrotaxis22 + Zcentre;
    }
    index++; Atom=LIGDEFtable[group][index];
   }

}

/*
#keyword group_name scale_factor id id id .. space 
group c2n2 1 16 18 19 20 21 22 23 24 25 
group c6o6 1 9 11 12 13 14 15 
group o7 1 78 80 81  
group o8 1 82 84 85 
group o9 1 86 89 90 
group c9o9 1 82 86 87 88 89 90 
*/

void ReadLIGDEF ( char * filename ) { 

FILE * File; 
char line[MaxLine];
char copy[MaxLine];
char  id[MaxLine];
int index=0;
 if ( ( File = fopen(filename,"r")) == NULL ) {
   printf("Error - can't open LIGDEF File %s\n",filename );
 }
 while ( fgets(line,LengthOfLine,File) != NULL ) {
    if (line[0]  != '#') {
        if ( strncmp(line,"group",5)== 0 ){ DEFgroups++;index=0;} 
        strcpy(copy,line);
// Get first token
        strtok(copy, " ");
// Get second token ( group_name )
        printf("Reading group %s ",strtok(NULL, " "));
	LIGDEFscale[DEFgroups]=atof(strtok(NULL, " "));
	LIGDEFtries[DEFgroups]=LigRotamerStates;  // For determinstic moves, rand=(LIGDEFtries[DEFgroups]/2)-1 ie (1,0.5,0,-0.5,-1)*angle*scale
	if ( LIGDEFdoubles[DEFgroups] == 1 ) LIGDEFdoubles[DEFgroups]=LigRotamerStates; 
        printf(" with scale factor %f\n  ",LIGDEFscale[DEFgroups]); 
        strcpy(id, strtok(NULL, " "));
        while ( isspace(id[0]) == 0  ) {
           printf("Reading %s into group %d\n",(id),DEFgroups);
           LIGDEFtable[DEFgroups][index] = atoi(id)-1 ; // List of atoms within a "GROUP" for Sugar rotations 
           printf("LIGDEF[%d][%d]=%d %s\n",DEFgroups,index,LIGDEFtable[DEFgroups][index],id); 
           strcpy(id, strtok(NULL, " "));
           index++;
        }
        LIGDEFtable[DEFgroups][index] = -1 ; // Terminate the end of atom list with -1
    }
 }
 if ( File != NULL ) fclose(File);


}

//  Tweak LIGAND ANGLE 
int AdjustLIGANDAngles ( FRAGMENTS * Structure)  {
extern float MY_RANDOM_FLOAT ( float );
extern int MY_RANDOM_INT ( int );
int group,i;
float rand;
int MadeChange; 
     MadeChange=0;
     if ( LigandResidue<0 ) { return 0; }
     if ( DEFgroups >= 0 && Tweaklig > 0 ) { 
         printf("############# ADJUST LIGAND ANGLE #######################\n");
         group=0;
         // Rather than sequentially choosing groups, choose randomnly
         if ( RandGroup == 1 ) { group=MY_RANDOM_INT(DEFgroups);}
         while ( MaxNligchange > 0 && group<=DEFgroups ) { 
             rand = MY_RANDOM_FLOAT(1);
             printf("Considering group %d Ligand res %d for rotations. threshold=%f rand=%f MaxNligchange=%d Tries=%f\n",group,LigandResidue,LIGThreshold,rand,MaxNligchange,LIGDEFtries[group]);
//rand=(LIGDEFtries[DEFgroups]/2)-1 ie (1,0.5,0,-0.5,-1)*angle*scale
             if ( (ExhaustiveSearch<0 && rand<LIGThreshold) || ( ExhaustiveSearch>=0 && LIGDEFtries[group]>= 0 ))  { 
              //group = MY_RANDOM_INT(DEFgroups);
              if ( ExhaustiveSearch>=0 ) {
                 if ( LIGDEFtries[group]==(LigRotamerStates/2) ) LIGDEFtries[group]--;  // don't bother with zero move
                 rand=(LIGDEFtries[group]/(LigRotamerStates/2))-1; 
                 printf("Angle for exhaustive ligand rotations %f\n",rand);
                 LIGDEFtries[group]--;
              }
              else { 
                rand=MY_RANDOM_FLOAT(2)-1; 
                printf("rand within threshold for ligand rotations\n");
              }
              printf("Adjusting LIGAND Angles rand=%f maxtweaklig =%f group %d groupscale=%f total angle=%f\n",rand,RadianToDegree*Tweaklig,group,LIGDEFscale[group], RadianToDegree*rand*Tweaklig*LIGDEFscale[group]);
              sprintf(remark,"%sREMARK AdjustLIGAND Angles group %d groupscale=%f total angle=%f\n",remark,group,LIGDEFscale[group], RadianToDegree*rand*Tweaklig*LIGDEFscale[group]);
              TransformLIG (Structure->res[LigandResidue],group, rand*Tweaklig*LIGDEFscale[group]) ;
              MadeChange=1;
              MaxNligchange--;
              if ( MaxNligchange == 0 ) break;
             }
             if ( RandGroup == 1 ) { group=MY_RANDOM_INT(DEFgroups);}
             else { group++;}
         }
     }
     printf("############# END ADJUST LIGAND ANGLE #######################\n");
     return MadeChange;
}

//  Tweak LIGAND ANGLE for a subset
int AdjustLIGANDAnglesSubset ( FRAGMENTS * Structure)  {
extern float MY_RANDOM_FLOAT ( float );
extern int MY_RANDOM_INT ( int );
int group,i;
float rand;
int MadeChange; 
     MadeChange=0;
     if ( LigandResidue<0 ) { return 0; }
     if ( DEFgroups >= 0 && Tweaklig > 0 ) { 
         printf("############# ADJUST LIGAND ANGLE SUBSET #######################\n");
         group=0;
         // Rather than sequentially choosing groups, choose randomnly
         if ( RandGroup == 1 ) { group=MY_RANDOM_INT(DEFgroups);}
         while ( MaxNligchange > 0 && group<=DEFgroups ) { 
             rand = MY_RANDOM_FLOAT(1);
             printf("Considering group %d Ligand res %d for rotations. threshold=%f rand=%f MaxNligchange=%d Doubles=%f\n",group,LigandResidue,LIGThreshold,rand,MaxNligchange,LIGDEFdoubles[group]);
             if ( ExhaustiveSearch>=0 && LIGDEFdoubles[group]>=0)  { 
              if ( ExhaustiveSearch>=0 ) {
                 // THIS VERSION ALLOWS ZERO ie don't move group and do single moves
                 rand=(LIGDEFdoubles[group]/(LigRotamerStates/2))-1; 
                 printf("Angle for exhaustive ligand rotations %f\n",RadianToDegree*rand*Tweaklig*LIGDEFscale[group]);
                 if ( group==LIGdouble1) LIGDEFdoubles[group]--;
              }
              else { 
                rand=MY_RANDOM_FLOAT(2)-1; 
                printf("rand within threshold for ligand rotations\n");
              }
              printf("Adjusting LIGAND Angles rand=%f maxtweaklig =%f group %d groupscale=%f total angle=%f\n",rand,RadianToDegree*Tweaklig,group,LIGDEFscale[group], RadianToDegree*rand*Tweaklig*LIGDEFscale[group]);
              sprintf(remark,"%sREMARK AdjustLIGAND Angles SUBSET group %d groupscale=%f total angle=%f\n",remark,group,LIGDEFscale[group], RadianToDegree*rand*Tweaklig*LIGDEFscale[group]);
              TransformLIG (Structure->res[LigandResidue],group, rand*Tweaklig*LIGDEFscale[group]) ;
              //MaxNligchange--;
              //if ( MaxNligchange == 0 ) break;
              if ( MadeChange==2) break;
              MadeChange++;
             }
             if ( RandGroup == 1 ) { group=MY_RANDOM_INT(DEFgroups);}
             else { group++;}
         }
printf("group1 %d doubles %f group2 %d %f \n",LIGdouble1,LIGDEFdoubles[LIGdouble1],LIGdouble2,LIGDEFdoubles[LIGdouble2]);
         // Need to reset first group and reduce second
         if ( LIGDEFdoubles[LIGdouble1]<=0) { 
printf("RESETING group1 %f group2 %f \n",LIGDEFdoubles[LIGdouble1],LIGDEFdoubles[LIGdouble2]);
            LIGDEFdoubles[LIGdouble1]=LigRotamerStates;
            LIGDEFdoubles[LIGdouble2]--;
         }
     }
     printf("############# END ADJUST LIGAND ANGLE  SUBSET #######################\n");
     return MadeChange-1;
}

/* Writes pdb format for each FRAGMENTS * */
int OutputXYZ ( FRAGMENTS * Pdb , char * filename) {
void OutputRESIDUEXYZ ( RESIDUE *  ,FILE * ,int * ) ;
FILE * AtmFile;
int current_atom,nres;
int modelno;
char chain=' ';
    //sprintf(filename,"%s_sup.pdb",Pdb->Label);
    printf("REMARK About to write XYZ %s numres=%d \n",filename,Pdb->numres);
    if ( strlen(filename)== 0 ) {
     AtmFile = stdout;
    }
    else if (  (AtmFile = fopen ( filename, "w" ))  == NULL ) {
        printf("Error opening %s for writing",filename);
        return JIVE_FALSE;
    }

    current_atom=1;
    modelno=UNDEFINED;
    chain=' ';
    for ( nres=0; nres < Pdb->numres; nres++) { 
        if ( Pdb->res[nres]->chain != chain) { 
             if ( chain!=' ') fprintf(AtmFile,"TER\n");
             chain=Pdb->res[nres]->chain ;
        }
        if ( Pdb->res[nres]->modelno != modelno) { 
          if ( modelno != UNDEFINED) { 
            fprintf(AtmFile,"TER\n"); fprintf(AtmFile,"ENDMDL\n");
          }
          fprintf(AtmFile,"MODEL %d\n",Pdb->res[nres]->modelno);
          /** Only select the first model Should this be in Baton?*/
          modelno=Pdb->res[nres]->modelno;
        }
        OutputRESIDUEXYZ ( Pdb->res[nres], AtmFile,&current_atom); 
    }
    current_atom=1;
    chain=' ';
    modelno=UNDEFINED;

    //DFB  Output Hetatm 
    for ( nres=0; nres < Pdb->numhetatm; nres++) { 
        if ( Pdb->hetatm[nres]->chain != chain) { 
             chain=Pdb->hetatm[nres]->chain ;
        }
        if ( Pdb->hetatm[nres]->modelno != modelno) { 
          if ( modelno != UNDEFINED) { 
            fprintf(AtmFile,"ENDMDL\n");
          }
          fprintf(AtmFile,"MODEL %d\n",Pdb->hetatm[nres]->modelno);
          /** Only select the first model Should this be in Baton?*/
          modelno=Pdb->hetatm[nres]->modelno;
        }
        OutputRESIDUEXYZ ( Pdb->hetatm[nres], AtmFile,&current_atom); 
    }
    if ( modelno != UNDEFINED ) { fprintf(AtmFile,"ENDMDL\n"); }
    fclose(AtmFile);
    return JIVE_TRUE;
}

void OutputRESIDUEXYZ ( RESIDUE * ThisRes ,FILE * AtmFile,int * current_atom ) {
int natom;
ATOM * ThisAtom;

    if ( ThisRes == NULL ) { printf("OutputRESIDUEL:NULL ThisRes\n");return;}
//printf("outputresidue natom %d\n",ThisRes->numatom);
    for ( natom=0; natom < ThisRes->numatom; natom++) {
      ThisAtom = ThisRes->atom[natom];
      if ( ThisAtom == NULL ) { printf("OutputRESIDUEL:NULL ThisAtom %d\n",natom);continue;}
      if ( ThisAtom->occupancy == UNDEFINED ) { printf("OutputRESIDUEL:NULL ThisAtom-occupancy %d\n",natom);continue;}
      if ( ThisAtom->x == UNDEFINED_COORD ) { printf("OutputRESIDUEL:NULL ThisAtom-xyz %d\n",natom);continue;}
      if ( isnan(ThisAtom->x) || isnan(ThisAtom->y) || isnan(ThisAtom->z) ) {
        printf("Residue %s %s atom[%d] %f %f %f has nan for xyz coordinates\n",ThisRes->Resname,ThisRes->Seqno,natom,ThisAtom->x,ThisAtom->y,ThisAtom->z);continue;
      }
      fprintf(AtmFile,"%8f ",ThisAtom->x);
      fprintf(AtmFile,"%8f ",ThisAtom->y);
      fprintf(AtmFile,"%8f ",ThisAtom->z);
      fprintf(AtmFile,"\n");
      (*current_atom)++;
    }
}


//Calcs RMSD of sidechain of Pdb1 before rotamer and Pdb2 after rotamer
float CalcRotamerRMSD( FRAGMENTS * Pdb1,FRAGMENTS * Pdb2, int res ) { 
int j;
int count=0;
float xsqr,ysqr,zsqr;
      j=ATOM_CB; 
      xsqr=ysqr=zsqr=0;
      while  ( j< Pdb1->res[res]->numatom && j < Pdb2->res[res]->numatom ) {
         count++;
         xsqr += Sqr(Pdb1->res[res]->atom[j]->x - Pdb2->res[res]->atom[j]->x);
         ysqr += Sqr(Pdb1->res[res]->atom[j]->y - Pdb2->res[res]->atom[j]->y);
         zsqr += Sqr(Pdb1->res[res]->atom[j]->z - Pdb2->res[res]->atom[j]->z);
//printf("atom %d %s %s count %d x %f y %f z %f\n",j,Pdb1->res[res]->atom[j]->atomname,Pdb2->res[res]->atom[j]->atomname,count,xsqr,ysqr,zsqr);
         j++;
       }
       if ( options.verbose > 3 ) printf("sum=%f sqrtsum = %f count=%d\n",xsqr+ysqr+zsqr,sqrt (xsqr+ysqr+zsqr),count);
       if ( count == 0 ) count = 1;
       return ( sqrt ( (float)(xsqr+ysqr+zsqr)/count) );
}


void OutputVariables (void){
   printf("start res %d endres %d \n",startres,endres); 
   printf("Checking Rand sequence \n"); 
   printf("LABEL LENGTH %d\n",strlen(NewRotLabel)); 
   printf("Ligand Residue is %d\n",LigandResidue);
   printf("Centre Residue is %d\n",CentreResidue);
   printf("PDBLABEL %s\n",PDB->Label);
   printf("NEWPDBLABEL %s\n",NewRotLabel);
   printf("DEFgroups %d Tweaklig %f\n",DEFgroups, Tweaklig); 
}
