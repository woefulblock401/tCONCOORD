#include <tconcoord.h>



void get_time(real sec,char *timestr)
{
   
  int d,h,m,s;
  real rest = 1;

  d = (int) sec / 86400.;
  
  rest = (int) sec % 86400;
  
  h = (int) rest / 3600.;
  
  rest = (int)rest % 3600;
  
  m = (int) rest / 60.;
  
  s = (int) rest % 60;
  
  sprintf(timestr,"Total computation time: %d days %d hours %d minutes %d seconds\n",d,h,m,s);
}

  
    


/*===================================================================*/


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "tdisco generates structure ensembles based on",
    "geometrical constraints."
  };


  char logstr[STRLEN];
  
  t_atomlist *al = NULL;
  t_resl *rl;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb = NULL;
  t_types *tp;
  t_bounds *b;
  t_dihed *impr;
  t_topology *top;
  t_atoms *atoms = NULL; 
  real       rlong = 0.8;
  rvec *x;
  int natom;

  gmx_rng_t rng;
  int seed=-1;
  int succ;
  int count;
  int impr_func=1;

  int i,k;

  bool bVerbose = FALSE;
  int maxit = 5000;
  int nstructs = 100;
  int nunsucc = 50;
  bool bPDB = FALSE;
  bool bPDBTRX = FALSE;
  bool bRandom = TRUE;
  bool bPert = FALSE;
  bool bRot = FALSE;
  int nrot = 5;
  bool bCheck = FALSE;
  real nbfac = .5;
  
  char *outf_base=NULL,*outf_ext=NULL;
  char out_file[STRLEN];
  char fname[STRLEN];

  int funcnr = 1;
  real bump_tol = 0.3;
  int dfunc = 0;
  bool pRand = FALSE;
  atom_id *index;
  char    *grpnames;
  int isize;
  bool bGroup = FALSE;
  bool bDoRef = TRUE;
  bool bTarget = FALSE;
  atom_id *targ_index;
  int targ_n;
  int nl = 10;
  bool bFit = TRUE;
  int cfreq = 1000;
  int start = 1;
  bool bTPX;
  real damp = 0;
  real maxf = 0.;
  int nsteps = 10;
  real gyr  = 0.;
  real gyr_tol = 0.;
  t_idxgroups *rot = NULL;
  char ref_title[STRLEN];
  t_topology ref_top;
  rvec *ref_x;
  matrix ref_box;

  t_topology targ_top;
  rvec *targ_x = NULL;
  char targ_title[STRLEN];
  matrix targ_box;

  atom_id *posr = NULL;
  int npos = 0;
  t_atomlist *tal = NULL;
  
  t_idxgroups *don = NULL;
  t_idxgroups *acc = NULL;
  t_idxgroups *phob = NULL;
  t_idxgroups *pl = NULL;
  t_idxgroups *imp = NULL;
  rvec bx;  
  t_boundtrack *bt = NULL;
  
  FILE *log = NULL;
  char *host = NULL;
  char *user = NULL;
  time_t t;
  struct tm *tt = NULL;
  char *ltime = NULL;

  FILE *ptrx = NULL;

  real *mass;
  rvec *xold, *xref;
  rvec co;
  rvec *newx;

  int unsucc = 0;

  int xtc;
  time_t tstart;
  real rmsd;
  time_t tend;  
  char timestr[STRLEN];

  t_gridmap *gp = NULL;
  
  t_pargs pa[] = {
     { "-v",   FALSE, etBOOL, {&bVerbose}, 
       "Make noise" },
     { "-n",   FALSE, etINT, {&nstructs}, 
       "Number of structures to be generated" },
     { "-i",   FALSE, etINT, {&maxit}, 
       "Maximum number of iterations" },
       { "-nsteps",   FALSE, etINT, {&nsteps},   
         "Number of EM steps" },  
     { "-start",   FALSE, etINT, {&start}, 
       "Start numbering with #" },
     {
       "-rand",FALSE,etBOOL,{&bRandom}, 
       "Start from random coordinates" },
     {
       "-pert",FALSE,etBOOL,{&bPert}, 
       "Only perturb starting configuration" },
     {
       "-nrot",FALSE,etINT,{&nrot}, 
       "Do # random sidechain rotations per round" },
      { 
       "-seed",FALSE,etINT,{&seed}, 
       "Initial random seed (-1 means make a seed)" },
     {
       "-nu",FALSE,etINT,{&nunsucc}, 
       "Maximum number of unsuccessful trials" },
     {
       "-check",FALSE,etBOOL,{&bCheck}, 
       "Check bounds" },
     {
       "-cfreq",FALSE,etINT,{&cfreq}, 
       "Check for improvement every # steps" },
     { 
       "-refine",FALSE,etBOOL,{&bDoRef},
       "Do refinement steps" },
     {
       "-fit",FALSE,etBOOL,{&bFit}, 
       "Fit structures to reference" },
     {
       "-gyr",FALSE,etREAL,{&gyr}, 
       "Target radius of gyration" },
     {
       "-gyr_tol",FALSE,etREAL,{&gyr_tol}, 
       "Radius of gyration tolerance" },
     
     {
       "-btol",FALSE,etREAL,{&bump_tol}, 
       "Bump tolerance" },
     {
       "-damp",FALSE,etREAL,{&damp}, 
       "damp factor for slow converging cases" },
       
  };


  t_filenm fnm[] = {
    { efCTP,"-top","tdist", ffREAD},
    { efTPX,"-s","topol", ffOPTRD}, 
    { efSTO,"-c","coords.pdb", ffOPTRD}, 
    { efXTC,"-x","tdisco", ffWRITE},
    { efDAT,"-d","tdist", ffREAD},
    { efNDX,"-grp","index", ffOPTRD},
    { efPDB,"-op","tdisco",ffOPTWR},
    { efPDB,"-on","tdisco_trn",ffOPTWR},
    { efPDB,"-ref","ref",ffWRITE},
    { efDAT,"-rot","rot",ffOPTRD},
    { efLOG,"-log","tdisco",ffWRITE},
    { efPDB,"-target","target",ffOPTRD},
    { efNDX,"-tidx","t_index",ffOPTRD}
  };
  
#define NFILE asize(fnm)

  cnc_copyright(argv[0]);
 
  parse_common_args(&argc,argv,0,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  bPDB  = opt2bSet("-op",NFILE,fnm);
  bPDBTRX  = opt2bSet("-on",NFILE,fnm);
  bRot = opt2bSet("-rot",NFILE,fnm);


  al=atomlist_init();


  log = ffopen(opt2fn("-log",NFILE,fnm),"w");

  host = getenv("HOSTNAME");
  user = getenv("USER");
  t = time(0);
  tt = localtime(&t);
  ltime = asctime(tt); 
  fprintf(log,"#=====================================\n");
  fprintf(log,"# tDISCO LOG FILE\n");
  fprintf(log,"#=====================================\n");
  fprintf(log,"# USER.......: %s\n",user);
  fprintf(log,"# TIME.......: %s",ltime);
  fprintf(log,"# HOST.......: %s\n",host);
  fprintf(log,"# CMDLINE....: %s\n",command_line());
#if (defined BUILD_MACHINE && defined BUILD_TIME && defined BUILD_USER)
  fprintf(log,
          "# tdisco was built on %s by\n"
          "# %s (%s)\n",BUILD_TIME,BUILD_USER,BUILD_MACHINE);
#endif
  fprintf(log,"#=====================================\n");
  fflush(log);
  
  don = idx_init();
  acc = idx_init();
  phob = idx_init();
  pl = idx_init();
  imp = idx_init();
  gp = gridmap_init();
  

  sprintf(logstr,"tCNC__log_> Reading Topology: %s\n",opt2fn("-top",NFILE,fnm));
  CNClog(log,logstr);
  read_cnctop(opt2fn("-top",NFILE,fnm),al,acc,don,phob,pl,imp);
  occ_to_one(al);
  sprintf(logstr,"tCNC__log_> Number of atoms: %d\n",al->natoms);
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of planar groups: %d\n",pl->n);
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of impropers: %d\n",imp->n);
  CNClog(log,logstr);

  
  if(opt2bSet("-c",NFILE,fnm)){

	read_tps_conf(opt2fn("-c",NFILE,fnm),ref_title,&ref_top,&ref_x,NULL,
                  ref_box,TRUE);
    rvec2al(al, ref_x);
    sprintf(logstr,"tCNC__log_> Reading reference coordinates from:  %s\n",opt2fn("-c",NFILE,fnm));
    CNClog(log,logstr);
  }
  
  snew(posr,al->natoms);

  for(i=0;i<al->natoms;i++){
    if(al->isposres[i]) {
      posr[npos] = i;
      npos++;
    }
  }
  

  if(opt2bSet("-grp",NFILE,fnm)){
    pRand = TRUE;
    atoms = al2atoms(atoms,al,&x);
    fprintf(stderr,"\nSelect a group to randomise:\n");
    get_index(atoms,opt2fn("-grp",NFILE,fnm),1,&isize,&index,&grpnames);
  }




  if(opt2bSet("-target",NFILE,fnm)){
	bTarget = TRUE;
	read_tps_conf(opt2fn("-target",NFILE,fnm),targ_title,&targ_top,&targ_x,NULL,
                targ_box,TRUE);
  
  tal = al_from_atoms(&(targ_top.atoms),targ_x);
  occ_to_one(tal);
  
    fprintf(stderr,"\nSelect the target group:\n");
    get_index(atoms,opt2fn("-tidx",NFILE,fnm),1,&targ_n,&targ_index,&grpnames);
  }


  rl = resl_init();
  rl->nres = count_res(al);
  sprintf(logstr,"tCNC__log_> Found %d atoms in %d residues\n",al->natoms,rl->nres);
  
  rl = resl_realloc(rl,rl->nres); 
  fill_resl(rl,al); 

  if(!bVerbose) {
    sprintf(logstr,"tCNC__log_> Reading bounds \n");
    CNClog(log,logstr);
  }
  
  b = read_bounds(opt2fn("-d",NFILE,fnm),bVerbose);
  sprintf(logstr,"tCNC__log_> Read %d constraints from: %s \n",b->n,opt2fn("-d",NFILE,fnm));
  CNClog(log,logstr);
    

  if (damp!=0) {
    sprintf(logstr,"tCNC__log_> Damping bounds: %4.2f\n",damp);
    CNClog(log,logstr);
    damp_bounds(b,damp);
  }

  bt = boundtrack_init();
  bt = boundtrack_realloc(bt,al->natoms);
  fill_boundtrack(b,bt);
  if(bCheck) {
    sprintf(logstr,"tCNC__log_> Checking bounds\n");
    check_bounds(log,b,al,pl,imp); 
    gp = nb(al, gp, rlong);
    fill_neighborlist(log, al, gp, rlong, UPDATE_NONBONDED);
    gp = reset_grid(gp);
    bumps_check(log,al,bt,bump_tol);
  }
  
  
  if(seed == -1)
    seed = gmx_rng_make_seed();
  rng = gmx_rng_init(seed);
  sprintf(logstr,"tCNC__log_> Initial random seed: %d\n",seed);
  CNClog(log,logstr);


  if(bPDBTRX)
  {
    ptrx = ffopen(opt2fn("-on",NFILE,fnm),"w");
    sprintf(logstr,"tCNC__log_> Opening pdb trajectory file: %s\n",opt2fn("-on",NFILE,fnm));
    CNClog(log,logstr);
  }

  if(bPDB){
    strcpy(out_file,opt2fn("-op",NFILE,fnm));
    outf_ext = strrchr(out_file,'.');
    if (outf_ext == NULL)
      gmx_fatal(FARGS,"Output file name '%s' does not contain a '.'",out_file);
    outf_base = strdup(out_file);
    outf_base[outf_ext - out_file] = '\0';
    sprintf(logstr,"tCNC__log_> Output each structure as pdb: %s\n",opt2fn("-op",NFILE,fnm));
    CNClog(log,logstr);
  }
  
  count = 0;

  snew(mass,al->natoms);
  snew(xold,al->natoms);
  snew(xref,al->natoms);
  snew(newx,al->natoms);
  
   if(npos){
     sprintf(logstr,"tCNC__log_> Position constraints on %d atoms\n",npos);
     CNClog(log,logstr);
     sprintf(logstr,"tCNC__log_> Ignoring these for center of mass calculation\n");
     CNClog(log,logstr);
     sprintf(logstr,"tCNC__log_> Getting center of mass\n");
     CNClog(log,logstr);
     com_frag(al,npos,posr,newx,co); 
   } 

   else{
     sprintf(logstr,"tCNC__log_> Getting center of mass\n");
     calc_cent(al,co);
   }
  
 
/*  if(bTarget){ */
/*    sprintf(logstr,"tCNC__log_> Shifting target coordinates\n"); */
/*    CNClog(log,logstr); */
/*    for(i=0;i<tal->natoms;i++){ */
/*      for(k=0;k<DIM;k++){ */
/*        tal->x[i][k]-=co[k]; */
/*      } */
/*    } */
/*  } */

/*  sprintf(logstr,"tCNC__log_> Shifting coordinates\n"); */
/*  CNClog(log,logstr); */

/*   for(i=0;i<al->natoms;i++){ */
/*     for(k=0;k<DIM;k++){ */
/*       al->x[i][k]-=co[k]; */
/*     } */
/*   } */
      
  for(i=0;i<al->natoms;i++){
    copy_rvec(al->x[i],xold[i]);
    copy_rvec(al->x[i],xref[i]);
    mass[i] = al->m[i];
  }


  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      xref[i][k]-=co[k];
    }
  }
      



  if(bTarget){
    calc_cent(tal,co);
  }
  sprintf(logstr,"tCNC__log_> Writing structure file: %s\n",opt2fn("-ref",NFILE,fnm));
  CNClog(log,logstr);
  write_pdb(al,opt2fn("-ref",NFILE,fnm)); 
  if(bTarget){
    sprintf(logstr,"tCNC__log_> Writing target structure file: target_out.pdb\n");
    CNClog(log,logstr);
    write_pdb(tal,"target_out.pdb");
    for(k=0;k<targ_n;k++){
      copy_rvec(tal->x[k],al->x[targ_index[k]]);
    }
  }
  
  if(bRot)
  {
    sprintf(logstr,"tCNC__log_> Reading rotation file: %s\n", opt2fn("-rot",NFILE,fnm));
    CNClog(log,logstr);
    rot = read_rotations(opt2fn("-rot",NFILE,fnm),bVerbose);
  }
  
  sprintf(logstr,"tCNC__log_> Opening xtc file: %s\n", opt2fn("-x",NFILE,fnm));
  CNClog(log,logstr);
  xtc = open_xtc(opt2fn("-x",NFILE,fnm),"w");
  

  sprintf(logstr,"tCNC__log_> Structure perturbation settings\n");
  CNClog(log,logstr);
  if(bRandom && !pRand && !bPert && !bRot)
  {
    sprintf(logstr,"tCNC__log_> Applying complete randomization prior to each run\n");
    CNClog(log,logstr);
  }
  else if(pRand)
  {
    sprintf(logstr,"tCNC__log_> Applying partly randomization prior to each run\n");
    CNClog(log,logstr);
  }
  else if(bPert){
    sprintf(logstr,"tCNC__log_> Applying random perturbations prior to each run\n");
    CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Coordinates will not be reset\n");
    CNClog(log,logstr);
  }
  else if(bRot){
    sprintf(logstr,"tCNC__log_> Applying %d random rotations prior to each run\n", nrot);
    CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Coordinates will not be reset\n");
    CNClog(log,logstr);
  }
  else if(!bRandom && !pRand && !bPert && !bRot){
    sprintf(logstr,"tCNC__log_> Coordinates remain unchanged\n");
    CNClog(log,logstr);
  }
  sprintf(logstr,"tCNC__log_> Generating %d structures\n",nstructs);
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> -----------------------------------------\n");  
  CNClog(log,logstr);
/*****************************************************************************/

  /* MAIN tCONCOORD LOOP */


  /* get the time */

  time(&tstart);
  

  int ntry = 1;
  
  while(count <  nstructs){
    sprintf(logstr,"tCNC__log_> Starting run %d\n",ntry);
    CNClog(log,logstr);
    
    if(bRandom && !pRand && !bPert && !bRot)
      random_al(al,rng,bTarget,co);  
    else if(pRand){
      for(i=0;i<al->natoms;i++){
        copy_rvec(xold[i],al->x[i]);
      }
      random_al_idx(al,isize,index,rng);
    }
    else if(bPert){
      perturb_al(al,rng);
    }
    else if(bRot){
      do_random_rotations(al, rot, nrot, rng);
/*          write_pdb(al,"start.pdb");     */
/*          exit(0);    */
      
    }
    
    else{
      for(i=0;i<al->natoms;i++){
        copy_rvec(xold[i],al->x[i]);
      }
    }
    
    
    succ = 0;    
    succ = do_disco(log,al,b,bt,impr,imp,rl,bump_tol,funcnr,dfunc,pl,rng,gp,
/*                     top,fr,md,  */
/*                     ir,grps,&nrnb,cr,cgs,nsb, */
                    rlong,nl,maxit,impr_func,bVerbose,
                    bTarget,tal,targ_n,targ_index,xold,cfreq,npos,posr,maxf,nbfac,bDoRef,nsteps,
                    gyr,gyr_tol);
    ntry+=1;
    char conv[STRLEN];
    if(succ) strcpy(conv,"successful");
    else strcpy(conv,"not successful");
/*     fprintf(log,"Run %5d Status..................: %s\n",ntry,conv); */
    
    if(succ){
      unsucc = 0;
      sprintf(logstr,"tCNC__log_> Generated Structure...............: %d (%d)\n",count+start,nstructs);
      CNClog(log,logstr);
      if(bFit){
        calc_cent(al, co);
        com(al);
        do_fit(al->natoms,mass,xref,al->x);
        rmsd = calc_rmsd(xref,al->x,al->natoms);
        sprintf(logstr,"tCNC__log_> RMSD..............................: %4.2f A\n",rmsd);
        CNClog(log,logstr);
        for(i=0;i<al->natoms;i++){
          for(k=0;k<DIM;k++){
            al->x[i][k]+=co[k];
          }
        }
      }
/*       else { */
/*           for(i=0;i<al->natoms;i++){ */
/*             for(k=0;k<DIM;k++){ */
/*               al->x[i][k]+=co[k]; */
/*             } */
/*           } */
/*       } */
      
      if(bPert || bRot){
        for(i=0;i<al->natoms;i++){
          copy_rvec(al->x[i],xold[i]);
        }
      }
      
      sprintf(logstr,"tCNC__log_> -----------------------------------------\n");  
      CNClog(log,logstr);
 
      al2xtc(al,xtc,count);
      if(ptrx) 
        write_pdb_frame(ptrx,al,count+start); 
      
      
      if(bPDB){ 
        sprintf(fname,"%s%d.pdb",outf_base,count+start); 
        write_pdb(al,fname); 
      } 
      count++;
    }
    else{
      sprintf(logstr,"tCNC__log_> -----------------------------------------\n");  
      CNClog(log,logstr);

      unsucc++;
      if(bPert || bRot){
        for(i=0;i<al->natoms;i++){
          copy_rvec(xold[i],al->x[i]);
        }
        unsucc--;
      }

      if(unsucc==nunsucc){
        
        sprintf(logstr,"tCNC__log_> Reached maximum number of unsuccessful trials\n");
        CNClog(log,logstr);
        sprintf(logstr,"tCNC__log_> Exiting....\n");
        CNClog(log,logstr);
        sprintf(logstr,"tCNC__log_> -----------------------------------------\n");
        CNClog(log,logstr);
        exit(1);
      }
    }

  }

  time(&tend);

  get_time(difftime(tend,tstart),timestr);
  
  sprintf(logstr,"tCNC__log_> %s\n",timestr);
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> -----------------------------------------\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Job terminated sucessfully\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Exiting\n");
  CNClog(log,logstr);
  return 0;
}

