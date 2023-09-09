#include <tconcoord.h>

int main(int argc, char **argv)
{
  static char *desc[] = {
    "macht tolles packing.....",
    ".....vielleicht..."
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
  real time;
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
       "Make noise" },
     { "-ind",   FALSE, etBOOL, {&bInd}, 
       "Calculate packing score for group of atoms" },
     { "-cutoff",   FALSE, etREAL, {&cutoff}, 
       "Calculate contacts for group of atoms" }

  };
  

  t_filenm fnm[] = {
    { efTPS,"-s","topol", ffREAD},
    { efXTC,"-f",NULL, ffOPTRD},
    { efNDX,"-n",NULL, ffOPTRD},
    { efDAT,"-o","packing.dat", ffWRITE},

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
  ct = contab_init();
  ct=contab_realloc(ct,al->natoms); 

  int size;
  atom_id *atomlist;
  char *names;
  

  if(opt2bSet("-n",NFILE,fnm)){
    fprintf(stderr,"Select a group for packing calculation\n");
    get_index(atoms,opt2fn("-n",NFILE,fnm),1,&size,&atomlist,&names);
    bIndex = TRUE;
  }
  printf("size = %d\n",size);
  
  
  if(bTop)
    bonds_from_idef(al,&(top->idef));
  
  
  atoms = al2atoms(atoms,al,&x); 
  clear_rvec(box_space);
  clear_mat(box);
  gen_box(0,atoms->nr,x,box,box_space,FALSE);
  
  /* avoid periodic stuff */
  for(i=0;i<DIM;i++){
    box[i][i] = box[i][i]+rlong+1;
  }
  
  
  natoms = atoms->nr; 
  
  snew(cg_index,natoms);
  for(i=0; (i<natoms); i++)
    cg_index[i]=i;
  
/*   snew(top,1); */
/*   init_top(top); */
  stupid_fill(&(top->blocks[ebCGS]),natoms,FALSE);
  memcpy(&(top->atoms),atoms,sizeof(*atoms));
  stupid_fill(&(top->atoms.excl),natoms,FALSE);
  top->atoms.grps[egcENER].nr = 1;
  
  /* Some nasty shortcuts */
  cgs  = &(top->blocks[ebCGS]);
  
  top->idef.ntypes = 1;
  top->idef.nodeid = 0;
  top->idef.atnr   = 1;
  snew(top->idef.functype,1);
  snew(top->idef.iparams,1);
  top->idef.iparams[0].lj.c6  = 1;
  top->idef.iparams[0].lj.c12 = 1;
  
  /* mdatoms structure */
  snew(nFreeze,2);
  md = atoms2md(debug,atoms,nFreeze,eiMD,0,0,NULL,FALSE,FALSE);
  sfree(nFreeze);
  
  /* nsborder struct */
  snew(nsb,1);
  nsb->nodeid  = 0;
  nsb->nnodes  = 1;
  calc_nsb(debug,&(top->blocks[ebCGS]),1,nsb,0,&(top->idef));
  if (debug)
    print_nsb(debug,"nsborder",nsb);
  
  /* inputrec structure */
  snew(ir,1);
  ir->coulombtype = eelCUT;
  ir->vdwtype     = evdwCUT;
  ir->ndelta      = 2;
  ir->ns_type     = ensGRID;
  snew(ir->opts.egp_flags,1);
  
  /* forcerec structure */
  if (fr == NULL)
    fr = mk_forcerec();
  snew(cr,1);
  cr->nnodes   = 1;
  cr->nthreads = 1;
  
  ir->rlist       = ir->rcoulomb = ir->rvdw = rlong;
/*   printf("Neighborsearching with a cut-off of %g\n",rlong); */
  init_forcerec(stdout,fr,ir,top,cr,md,nsb,box,FALSE,NULL,NULL,TRUE);
  fr->cg0 = 0;
  fr->hcg = top->blocks[ebCGS].nr;
  fr->nWatMol = 0;
  if (debug)
    pr_forcerec(debug,fr,cr);
  
  /* Prepare for neighboursearching */
  init_nrnb(&nrnb);
  
  /* Group stuff */
  snew(grps,1);
  
/* Init things dependent on parameters */  
  ir->rlist       = ir->rcoulomb = ir->rvdw = rlong;
/*   printf("Neighborsearching with a cut-off of %g\n",rlong); */
  init_forcerec(debug,fr,ir,top,cr,md,nsb,box,FALSE,NULL,NULL,TRUE);
  
  /* Calculate new stuff dependent on coords and box */
  for(m=0; (m<DIM); m++)
    box_size[m] = box[m][m];
  calc_shifts(box,fr->shift_vec);
  put_charge_groups_in_box(NULL,0,cgs->nr,box,cgs,x,fr->cg_cm);
  
  /* Do the actual neighboursearching */
  init_neighbor_list(NULL,fr,HOMENR(nsb));
  
  
  
  search_neighbours(stderr,fr,x,box,top,grps,cr,nsb,&nrnb,md,lambda,&dvdlambda,
                    TRUE,FALSE);
  nlist = &(fr->nblists[0].nlist_sr[eNL_VDW]);
  fill_nl(al,nlist,ct);  
  
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
  
  FILE *of = ffopen(opt2fn("-o",NFILE,fnm),"a");

  if(opt2bSet("-f",NFILE,fnm)){
    
    int fp = open_xtc(opt2fn("-f",NFILE,fnm),"r");
    
    
    fprintf(of,"#============================\n");
    fprintf(of,"# packing score calculation\n");
    fprintf(of,"#============================\n");
    fprintf(of,"# time      packing score  \n");
    
    if (!read_first_xtc(fp,&natoms,&step,&time,box,&x,
                        &prec,&bOk)) {
      fprintf(stderr,"XTC sucks.....\n");
      exit(1);
    }
    do{
      rvec2al(al,x);
      
      update_neighborlist(stderr,al, x, top,md, fr,ir,   
                          grps, &nrnb, cr, cgs,nsb, 
                          rlong,0,0,NULL);   
      fprintf(stderr,"Reading frame %10.3fps\r",time);
      
      if(bIndex){
        fprintf(of,"%10.3f %10.3f %10d\n",time,packing_score_subset(al,atomlist,size),count_contacts(al,atomlist,size,cutoff));     
      }
      else{
        fprintf(of,"%10.3f %10.3f\n",time,packing_score(al));     
      }
      
    }while(read_next_xtc(fp,natoms,&step,&time,box,x,&prec,&bOk));
    fprintf(stderr,"\n\nHabe fertig\n");
  }
  else{
    if(bIndex){
      fprintf(of,"%-15s  %10.3f %10d\n",opt2fn("-s",NFILE,fnm),packing_score_subset(al,atomlist,size),
              count_contacts(al,atomlist,size,cutoff));
    }
    else {
      fprintf(of,"%-15s  %10.3f\n",opt2fn("-s",NFILE,fnm),packing_score(al));     
    }
  }
  
    
    return 0;
}

  
