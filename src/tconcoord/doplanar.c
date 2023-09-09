#include <tconcoord.h>

void print_group_to_file(FILE *fp,t_atomlist *al,int *idx, int n)

/* self explaining */

{
  int i;
  fprintf(fp,"group: %d atoms\n",n);
  for(i=0;i<n;i++){
    fprintf(fp,"%d (%s/%d/%s)\n",al->id[idx[i]],al->name[idx[i]],
           al->resid[idx[i]],al->resname[idx[i]]);
  }
}

int main(int argc, char **argv)
  
{
  
  t_atomlist *al;
  t_contab *ct;
  t_resl *rl;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb;
  t_vdwcomb *vdw14comb;
  t_types *tp;
  t_bounds *b;
  t_dihed *impr;
  t_topology top;
  t_atoms *atoms = NULL; 
  char title[STRLEN];
  rvec       box_size;
  rvec *x;
  matrix box;
  bool bVerbose = FALSE;
  
  int i,k,j;
  real plan;
  gmx_rng_t rng;
  real def = 0.5;
  bool bGrp=FALSE;
  atom_id *index;
  char    *grpnames;
  char *posnames;
  int isize;
  int psize = 0;
  atom_id *pos_id;
  int seed = -1;
  int merge = 3;
  bool bNOH = FALSE;
  bool bpH = FALSE;
  bool bIgn = FALSE;
  bool bPosres = FALSE;
  atom_id *flex;
  int flex_size = 0;
  atom_id *target;
  int ntarg = 0;
  char warn[STRLEN];
  bool bTop;

  
  static char *desc[] = {  
  };
  
  
  t_pargs pa[] = {
  };
  
  
  
  t_filenm fnm[] = {
/*     { efPDB,"-f","protein", ffOPTRD },  */
    { efTPS,"-s","topol", ffREAD }, 
    { efDAT,"-o","pldata", ffWRITE } 
  };
  
  
#define NFILE asize(fnm)
  
/*   CopyRight(stderr,argv[0]); */
  cnc_copyright(argv[0]);
  
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  read_tps_conf(opt2fn("-s",NFILE,fnm),title,&top,&x,NULL,
                box,TRUE);
  bTop=fn2bTPX(opt2fn("-s",NFILE,fnm));
  
  al = al_from_atoms(&(top.atoms),x);
  
  if(bTop)
    bonds_from_idef(al,&(top.idef));
  
  get_symbol(stderr,al);
  
  ct = contab_init();
  ct=contab_realloc(ct,al->natoms); 
  
  nb_search(al,ct,0.8,!bTop);
  
  rename_at(al);
  get_order(al);
  vdwcomb = read_vdwcomb(1,FALSE);
  vdw14comb = read_vdwcomb(1,TRUE);
  tp = read_atom_types(1);
  get_types(stderr,al,tp);
  vdw = read_vdw_radii(1);  
  get_vdw_radii(stderr,al,vdw,vdwcomb,vdw14comb);
  
  get_hybrid2(al,tp);
  renumber_atoms(al,1);
  renumber_residues(al,1);
  get_cnc_solv(al);
  
  
  t_bondlist *bl = bl_init();
  rotatable_bonds(al,bl);
  fprintf(stderr,"number of bonds..................: %d\n",bl->n);
  
  
  
  
/* get planar groups */
  rl = resl_init();
  rl->nres = count_res(al);
  rl = resl_realloc(rl,rl->nres);
  fill_resl(rl,al);
  t_idxgroups *pln = idx_init();
  get_planar_groups(stderr,al,rl,pln,0);
  
  /* check non-protein residues */
  
  for(i=0;i<bl->n;i++){
    if(!IsProtein(al,bl->at1[i])){
      if(bl->type[i] == DOUBLE ||
         bl->type[i] == CNDOUBLE ||
         bl->type[i] == CNOMEGA ||
         bl->type[i] == NNDOUBLE ||
         bl->type[i] == TRIPLE ||
         bl->type[i] == CNTRIPLE){
        pln->n++;
        pln = idx_realloc(pln,pln->n);
        add_to_group(pln,pln->n-1,bl->at1[i]);
        add_to_group(pln,pln->n-1,bl->at2[i]);
        pln->flag[pln->n-1] = TRUE;
        for(k=0;k<al->nbonds[bl->at1[i]];k++){
          add_to_group(pln,pln->n-1,al->bonds[bl->at1[i]][k]);
        }
        for(k=0;k<al->nbonds[bl->at2[i]];k++){
          add_to_group(pln,pln->n-1,al->bonds[bl->at2[i]][k]);
        }
        
        if(!is_planar_group(al,pln->natoms[pln->n-1],pln->atoms[pln->n-1],
                            0.03,&plan)){
          pln->n--;
        }
        else {
          if (plan < 0.001) {
            plan = 0.001;
          }
          
          pln->val[pln->n-1] = plan;
        }
        
      }
    }
  }
  
  
  
  for(i=0;i<al->natoms;i++){
    if(!IsProtein(al,i) && strcmp(al->hyb[i],"sp2") == 0 
       && al->nbonds[i] >= 2){
      pln->n++;
      pln = idx_realloc(pln,pln->n);
      add_to_group(pln,pln->n-1,i);
      pln->val[pln->n-1] = 0.03;
      for(k=0;k<al->nbonds[i];k++){
        add_to_group(pln,pln->n-1,al->bonds[i][k]);
      }
    }
  }
  
  pln = remove_redundant_groups(pln);
  
  
  
  /* get rings */
  
  t_idxgroups *rings = idx_init();
  for(i=0;i<al->natoms;i++){
    if(!IsProtein(al,i) && strcmp(al->symbol[i],"H") != 0)
      get_ring(al,i,rings,bl);
  }
  
  rings = remove_redundant_groups(rings);
  
  
  for(i=0;i<rings->n;i++){
    if(is_planar_group(al,rings->natoms[i],rings->atoms[i],0.03,&plan)){
      pln->n++;
      pln = idx_realloc(pln,pln->n);
      for(k=0;k<rings->natoms[i];k++){
        add_to_group(pln,pln->n-1,rings->atoms[i][k]);
      }
    }
  }
  
  pln = remove_redundant_groups(pln);
  pln = merge_groups(pln,3,0);  
  pln = remove_redundant_groups(pln);  
  

  FILE *out =   ffopen(opt2fn("-o",NFILE,fnm),"w");
  
  for(i=0;i<pln->n;i++){
    is_planar_group(al,pln->natoms[i],pln->atoms[i],0.03,&plan);
    fprintf(out,"//\n");
    print_group_to_file(out,al,pln->atoms[i],pln->natoms[i]);
    fprintf(out,"pl = %g\n",plan);
    
  }
  
  
  return 0;
}
