#include <tconcoord.h>
#include <rmpbc.h>


int main(int argc, char **argv)
{
  static char *desc[] = {
    "Analyses hbonds",
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
  bool bDrop = FALSE;
  real rad = 3.;
  real ang = 120.;
  
  
  int        *cg_index;
  int        i,m,natoms;
  int step;
  real time;
  real prec;
  bool bOk;
  real d;
  real a;
  int at;
  
  ivec       *nFreeze;
  rvec       box_size;
  real       lambda=0,dvdlambda=0;
  real       rlong = 0.8;
 
  matrix box;
  rvec box_space;
  rvec *x;
  int k,j,don,acc;
  int atid = 0;
  

  
  t_pargs pa[] = {
    { "-v",   FALSE, etBOOL, {&bVerbose},"Make noise"} ,
    {"-drop",   FALSE, etBOOL, {&bDrop},"drop water nl"},
    {"-dist",   FALSE, etREAL, {&rad},"max dist"},
    {"-ang",   FALSE, etREAL, {&ang},"min angle"}
    
    };

  

  t_filenm fnm[] = {
    { efTPS,"-s","topol", ffOPTRD},
    { efXTC,"-f",NULL, ffREAD},
    { efNDX,"-n","index", ffOPTRD}, 
    { efDAT,"-o","hb.dat", ffWRITE}

  };
  
  #define NFILE asize(fnm)

  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  atid-=1;
  snew(top,1);
  
  read_tps_conf(opt2fn("-s",NFILE,fnm),title,top,&x,NULL,
                box,TRUE);
  bTop=fn2bTPX(opt2fn("-s",NFILE,fnm));

  al = al_from_atoms(&(top->atoms),x);
  get_symbol(stderr,al);
  ct = contab_init();
  ct=contab_realloc(ct,al->natoms); 

  
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
  
   int size; 
   atom_id *id; 
   char *names; 
   get_index(atoms,opt2fn("-n",NFILE,fnm),1,&size,&id,&names); 
  
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

  int xtc = open_xtc(opt2fn("-f",NFILE,fnm),"r");
/*   FILE *of = ffopen(opt2fn("-o",NFILE,fnm),"w"); */
  
/*   fprintf(of,"#============================\n"); */
/*   fprintf(of,"# packing score calculation\n"); */
/*   fprintf(of,"#============================\n"); */
/*   fprintf(of,"# time      packing score  \n"); */



  for(i=0;i<size;i++){
    atid = id[i];
    fprintf(stderr,"Calculating hbonds for atom %d %s (%d %s)\n",
            al->id[atid],al->name[atid],al->resid[atid],al->resname[atid]);
  }
  


  if (!read_first_xtc(xtc,&natoms,&step,&time,box,&x,
                      &prec,&bOk)) {
    fprintf(stderr,"XTC sucks.....\n");
    exit(1);
  }
/*   real d; */
/*   real a; */
/*   int at; */
  FILE *fp[size];
  char outn[STRLEN];
  char dum[STRLEN];
  
  for(i=0;i<size;i++){
    strcpy(outn,"hbonds_");
    sprintf(dum,"%d.dat",al->id[id[i]]);
    strcat(outn,dum);
    fp[i] = fopen(outn,"w");
  }
  
  
  do{


    real lambda, dvdlambda;
    lambda = dvdlambda = 0;
    t_nblist *nlist;
    calc_shifts(box,fr->shift_vec); 
    put_charge_groups_in_box(NULL,0,cgs->nr,box,cgs,x,fr->cg_cm); 

    search_neighbours(stderr,fr,x,box,top,grps,cr,
                      nsb,&nrnb,md,lambda,&dvdlambda,
                      TRUE,FALSE);
    nlist = &(fr->nblists[0].nlist_sr[eNL_VDW]);
    rvec2al(al,x);    
    update_nl(al,nlist);

    fprintf(stderr,"Reading frame %10.3fps\r",time);
    
      

    for(k=0;k<size;k++){
      fprintf(fp[k],"time = %g\n",time);
/*       fprintf(fp[k]," %g %g %g\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]); */
      
      atid = id[k];
      don = al->bonds[atid][0];
      for(i=0;i<al->nnb[atid];i++){
        at = al->nb[atid][i];
        /* printf("name = '%s'\n",al->name[i]); */
        
        if(strcmp(al->name[at]," OW ") == 0){
          
          d = DIST(al,atid,at);
          /*  a = RAD2DEG*ANGLE(al,atid,at,don); */
          a = RAD2DEG*angle_ij_ik(al->x[atid],al->x[at],al->x[don]);
          
/*         printf("d = %g a = %g  (%s-%s-%s)\n",d,a,al->name[atid],al->name[at],al->name[don]); */
          if (!bDrop){
            if (d< rad && a > ang){
              fprintf(fp[k],"%8d %8.3f %8.3f\n",al->id[at],d,a);
            }
          }
          else {
            fprintf(fp[k],"%8d %8.3f %8.3f\n",al->id[at],d,a);
          }
          
        }
      }
      fprintf(fp[k],"//\n");
      fflush(fp[k]);
    }
    



    
/*     fprintf(of,"%10.3f %10.3f\n",time,packing_score(al));      */
    
  }while(read_next_xtc(xtc,natoms,&step,&time,box,x,&prec,&bOk));
/*   fprintf(stderr,"\n\nHabe fertig\n"); */
  


  return 0;
}

  
