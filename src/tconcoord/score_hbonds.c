#include<tconcoord.h>

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


int main(int argc, char **argv)
{
  static char *desc[] = {
   "Calculate Lazaridis-Karplus solvation energy",
     "........"
  };

  t_atomlist *al;
  t_contab *ct = NULL;
  t_resl *rl;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb;
  t_vdwcomb *vdw14comb;
  t_types *tp;
  t_forcerec *fr = NULL;
  t_topology *top = NULL; 
  t_atoms *atoms = NULL; 
  t_nsborder *nsb;
  t_mdatoms  *md;
  t_block    *cgs;
  t_inputrec *ir;
  t_nrnb     nrnb;
  t_commrec  *cr;
  t_groups   *grps;
  t_nblist   *nlist;
  bool bVerbose = FALSE;
  char title[STRLEN];
  bool bTop = FALSE;
  int params = 1;
  
  int        *cg_index;
  int        i,m,natoms;
  int step;
  real timestep;
  real prec;
  bool bOk;
  
  ivec       *nFreeze;
  rvec       box_size;
  real       lambda=0,dvdlambda=0;
  real       rlong = 0.8;
 
  matrix box;
  rvec box_space;
  rvec *x;
  real bInd = FALSE;
  bool bIndex = FALSE;
  real cutoff = 4.;
  

  t_pargs pa[] = {
     { "-v",   FALSE, etBOOL, {&bVerbose}, 
       "Make noise" }
/*      { "-ind",   FALSE, etBOOL, {&bInd},  */
/*        "Calculate packing score for group of atoms" }, */
/*      { "-cutoff",   FALSE, etREAL, {&cutoff},  */
/*        "Calculate contacts for group of atoms" } */
     
  };
  
  
  t_filenm fnm[] = {
    { efTPS,"-s","topol", ffREAD},
    { efXTC,"-f",NULL, ffOPTRD},
/*     { efNDX,"-n",NULL, ffOPTRD}, */
    { efDAT,"-o","hbonds.dat", ffWRITE},
    
  };
  
#define NFILE asize(fnm)
  
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  
  snew(top,1);
  
  read_tps_conf(opt2fn("-s",NFILE,fnm),title,top,&x,NULL,
                box,TRUE);
  bTop=fn2bTPX(opt2fn("-s",NFILE,fnm));
  
  al = al_from_atoms(&(top->atoms),x);
  get_symbol(stderr,al);
  if(bTop)
    bonds_from_idef(al,&(top->idef));

  t_gridmap *gp = gridmap_init(); 
  gp = nb(al, gp,9);   
  if(bTop)
    fill_neighborlist(stderr, al, gp, 9.,UPDATE_NONBONDED);   
  else
    fill_neighborlist(stderr, al, gp, 9.,MAKE_FULL_NEIGHBORLIST);   
  gp = reset_grid(gp);   

  int size;
  atom_id *atomlist;
  char *names;
  
    
  rename_at(al); 
  get_order(al);
  vdwcomb = read_vdwcomb(params,FALSE);
  vdw14comb = read_vdwcomb(params,TRUE);
  tp = read_atom_types(params);
  get_types(stderr,al,tp);
  vdw = read_vdw_radii(params);  
  get_vdw_radii(stderr,al,vdw,vdwcomb,vdw14comb);
  
  get_hybrid2(al,tp);
  renumber_atoms(al,1);
  renumber_residues(al,1);
  get_rosetta_types(al);
  get_cnc_solv(al);
  FILE *of = ffopen(opt2fn("-o",NFILE,fnm),"w");
  t_idxgroups *don = idx_init(); 
  t_idxgroups *acc = idx_init();
  don = idx_realloc(don,1);
  acc = idx_realloc(acc,1); 
  hbonds(al,don,acc); 
  
  if(opt2bSet("-f",NFILE,fnm)){
  
    int fp = open_xtc(opt2fn("-f",NFILE,fnm),"r");
  
/*   fprintf(of,"#============================\n"); */
/*   fprintf(of,"# Solvation energy\n"); */
/*   fprintf(of,"#============================\n"); */
/*   fprintf(of,"# time      energy score  \n"); */
  
    if (!read_first_xtc(fp,&natoms,&step,&timestep,box,&x,
                        &prec,&bOk)) {
      fprintf(stderr,"XTC sucks.....\n");
      exit(1);
    }
    do{
      rvec2al(al,x);
      gp = nb(al, gp,9);   
      fill_neighborlist(stderr,al, gp, 9., UPDATE_NONBONDED);   
      gp = reset_grid(gp);
      get_hbonds(al, 2.6, 125., 6.);  
      
      fprintf(stderr,"Reading frame %10.3fps\r",timestep);
      fprintf(of,"# TIME %g\n",timestep);
      print_hbonds_to_file(of,al);  
      free_hbonds(al);   
/*       snew(al->hbonds,1); */
      
/*     fprintf(of,"%10.3f %10.3f %10.3f %10.3f \n",timestep,solv_en, lj_en, solv_en+lj_en);      */
      
    }while(read_next_xtc(fp,natoms,&step,&timestep,box,x,&prec,&bOk));
  }
  
  else {
    get_hbonds(al, 2.6, 125., 6.); 
    print_hbonds_to_file(of,al); 
  }
  
  return 0;
}

