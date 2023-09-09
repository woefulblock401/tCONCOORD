#include <tconcoord.h>

void print_cell(int *cell, int nat)
{
  int i;
  t_idxgroups *dd = idx_init();
  dd->n = 1;
  for(i=0;i<nat;i++){
    add_to_group(dd,0,cell[i]+1);
  }
  write_group_script("ids.pml","id",dd);
  
/*   FILE *fp = fopen("ids.pml","w"); */
/*   fprintf(fp,"select xx, id "); */
/*   for(i=0;i<nat;i++){ */
/*     fprintf(fp,"%d+",cell[i]+1); */
/*   } */
}


/*============================================================*/
void fill_nl(t_atomlist *al,t_nblist *nlist,t_contab *ct)
{
  int i,j,k;
  int j0,j1;
  int inr,jnr;
  real d;
  
  /* do bonds */

  for(i=0;i<nlist->nri;i++){ 
    inr = nlist->iinr[i];
    j0  = nlist->jindex[i];
    j1  = nlist->jindex[i+1];
    for(j=j0; (j<j1); j++) {
      jnr = nlist->jjnr[j]; 
      d = DIST2(al,inr,jnr);
       if(d < sqr(al->bcontr[inr]+al->bcontr[jnr]) &&
          (strcmp(al->symbol[inr],"H")!=0 || strcmp(al->symbol[jnr],"H") !=0)){
        add_bond(al,inr,jnr,ct,BOND,BOND);
      }
    }
  }
  /* fill 1-3 array */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      inr = al->bonds[i][k];
      for(j=0;j<al->nbonds[inr];j++){
        if(al->bonds[inr][j] < i){
          add_bond(al,i,al->bonds[inr][j],ct,B13,B13);
        }
      }
    }
  }
  
  /* fill 1-4 array */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nb13[i];k++){
      inr = al->b13[i][k];
      for(j=0;j<al->nbonds[inr];j++){
        jnr = al->bonds[inr][j];
        if(i<jnr && !connected(ct,i,jnr)){
          add_bond(al,i,jnr,ct,B14,B14);   
        }
      }
    }
  }
  /* put the rest to the neighborlist */

  for(i=0;i<nlist->nri;i++){ 
    inr = nlist->iinr[i];
    j0  = nlist->jindex[i];
    j1  = nlist->jindex[i+1];
    for(j=j0; (j<j1); j++) {
      jnr = nlist->jjnr[j]; 
      if(!connected(ct,inr,jnr)){   
        add_bond(al,inr,jnr,ct,NL,NL);   
      } 
     } 
  }
}
/*============================================================*/

void do_simple_search(t_atomlist *al, t_contab *ct,real max, bool bUpdate)
/* do a quick neighbor search for small molecules */
{
  int i,k,j,inr,jnr;
  real d;

  if(!bUpdate) {
    if(ct!=NULL)
      ct=contab_realloc(ct,al->natoms); 
    
    
    for(i=0;i<al->natoms;i++){
      for(k=i+1;k<al->natoms;k++){
        d = DIST2(al,i,k);
        if(d < sqr(al->bcontr[i]+al->bcontr[k]) &&
           (strcmp(al->symbol[i],"H")!=0 || strcmp(al->symbol[k],"H") !=0)){
          add_bond(al,i,k,ct,BOND,BOND);
        }
      }
    }
    
    
    /* fill 1-3 array */
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nbonds[i];k++){
        inr = al->bonds[i][k];
        for(j=0;j<al->nbonds[inr];j++){
          if(al->bonds[inr][j] < i){
            add_bond(al,i,al->bonds[inr][j],ct,B13,B13);
          }
        }
      }
    }
    
    /* fill 1-4 array */
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nb13[i];k++){
        inr = al->b13[i][k];
        for(j=0;j<al->nbonds[inr];j++){
          jnr = al->bonds[inr][j];
          if(i<jnr && !al_connected(al,i,jnr)){
            add_bond(al,i,jnr,ct,B14,B14);   
          }
        }
      }
    }
  }
  
  /* put the rest to the neighborlist */
  max=sqr(max*10.);
  for(i=0;i<al->natoms;i++){
    for(k=i+1;k<al->natoms;k++){
      d = DIST2(al,i,k);
      if(d < max && !al_connected(al,i,k)){
        add_bond(al,i,k,ct,NL,NL);
      }
    }
  }
  
}



/*============================================================*/
  

void nb_search(t_atomlist *al,
               t_contab *ct,real rlong, bool bAll)

/* Neighbor search. For small systems a simple 
   search algorithm is used. Larger systems use
   the gromacs grid based neighbor search.
*/

{
  if(al->natoms < 500)
    do_simple_search(al,ct,rlong, FALSE);

  else{
    t_forcerec *fr = NULL;
    t_topology *top;
    t_atoms *atoms = NULL; 
    t_nsborder *nsb;
    t_mdatoms  *md;
    t_block    *cgs;
    t_inputrec *ir;
    t_nrnb     nrnb;
    t_commrec  *cr;
    t_groups   *grps;
    t_nblist   *nlist;
    
    int        *cg_index;
    int        i,m,natoms;
    ivec       *nFreeze;
    rvec       box_size;
    real       lambda=0,dvdlambda=0;
    
    matrix box;
    rvec box_space;
    rvec *x;
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
    
    snew(top,1);
    init_top(top);
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

    if(bAll)
      fill_nl(al,nlist,ct);  
    else
      update_nl(al,nlist);
    
    sfree(grps);
    sfree(ir->opts.egp_flags);
    sfree(top->idef.functype);
    sfree(top->idef.iparams);
    sfree(nsb);
    sfree(cr);
    sfree(cg_index);
    sfree(x);
    sfree(atoms->atom);
    sfree(atoms);
    
  }
}
/*============================================================*/


void update_nl(t_atomlist *al, t_nblist *nlist)
{
  int i,j,inr,jnr,k;
  int j0,j1,nb1,nb2;
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      al->nb[i][k] = 0;
    }
    al->nnb[i] = 0;
/*     srenew(al->nb[i],MAX_NL); */
    
  }
  
  for(i=0;i<nlist->nri;i++){ 
    inr = nlist->iinr[i];
    j0  = nlist->jindex[i];
    j1  = nlist->jindex[i+1];
    for(j=j0; (j<j1); j++) {
      jnr = nlist->jjnr[j]; 
      if(!al_connected(al,inr,jnr)){   
        al->nnb[inr]+=1;
        al->nnb[jnr]+=1;
        nb1 = al->nnb[inr];
        nb2 = al->nnb[jnr];
        if(nb1 > MAX_NL*al->memcount[inr]){
          al->memcount[inr]+=1;
          srenew(al->nb[inr],MAX_NL*al->memcount[inr]);
        }
        if(nb2 > MAX_NL*al->memcount[jnr]) {
          al->memcount[jnr]+=1;
          srenew(al->nb[jnr],MAX_NL*al->memcount[jnr]);
        }
        al->nb[inr][nb1-1] = jnr;
        al->nb[jnr][nb2-1] = inr;
      }
    }
  }
}
/*============================================================*/  
void update_nl2(t_atomlist *al, t_nblist *nlist, bool trunc)
{
  int i,j,inr,jnr,k;
  int j0,j1,nb1,nb2;
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      al->nb[i][k] = 0;
    }
    al->nnb[i] = 0;
/*     srenew(al->nb[i],MAX_NL);  */
    
  }
  
  for(i=0;i<nlist->nri;i++){ 
    inr = nlist->iinr[i];
    j0  = nlist->jindex[i];
    j1  = nlist->jindex[i+1];
    for(j=j0; (j<j1); j++) {
      jnr = nlist->jjnr[j]; 
      if(!al_connected(al,inr,jnr)){   
        bool nb1Ok = FALSE;
        bool nb2Ok = FALSE;
        


        nb1 = al->nnb[inr];
        nb2 = al->nnb[jnr];
        if(nb1 > MAX_NL){
          if(!trunc){
            al->nnb[inr]+=1;
            srenew(al->nb[inr],nb1);
            nb1Ok = TRUE;
          }
          else{
            nb1Ok = FALSE;
/*             printf("Truncated neighborlist\n"); */
          }
        }
        else {
          al->nnb[inr]+=1;
          nb1 = TRUE;
        }
        if(nb2 > MAX_NL) {
          if(!trunc){
            al->nnb[jnr]+=1;
            srenew(al->nb[jnr],nb2);
            nb2Ok = TRUE;
          }
          else{
            nb2Ok = FALSE;
/*             printf("Truncated neighborlist\n"); */
          }
        }
        else{
          al->nnb[jnr]+=1;
          nb2Ok = TRUE;
        }
        if(nb1Ok)
          al->nb[inr][nb1-1] = jnr;
        if(nb2Ok)
          al->nb[jnr][nb2-1] = inr;
      }
    }
  }
}


/*============================================================*/
void dump_nl(FILE *fp,t_atomlist *al,real max)
{
  int i,k;
  real d;
  for(i=0;i<al->natoms;i++){
    fprintf(fp,"id: %d name: %s res: %d\n",al->id[i],al->name[i],al->resid[i]);
    for(k=0;k<al->nnb[i];k++){
      d = DIST(al,i,al->nb[i][k]);
      if(d < max)
        fprintf(fp,"\tid: %d name: %s res: %d dist: %g\n",al->id[al->nb[i][k]],
                al->name[al->nb[i][k]],al->resid[al->nb[i][k]],d);
    }
  }
}

/*============================================================*/

void update_neighborlist(FILE *log,t_atomlist *al, rvec *x, t_topology *top,
                         t_mdatoms *md, t_forcerec *fr,
                         t_inputrec *ir,
                         t_groups *grps, t_nrnb *nrnb,
                         t_commrec *cr, t_block *cgs,
                         t_nsborder *nsb, real rlong, bool bSimple,
                         int nat, atom_id *ids)
{
  if(bSimple){
/*     printf("Doing simple search\n"); */
    int i,k;
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nnb[i];k++){
        al->nb[i][k] = 0;
      }
      al->nnb[i] = 0;
      srenew(al->nb[i],300);
    }

    real d;
    real rlsq = sqr(rlong);
    for(i=0;i<nat;i++){
      for(k=0;k<al->natoms;k++){
        d = dist2(al->x[ids[i]],al->x[k]);
        if(d < rlsq && !al_connected(al,ids[i],k)) {
             al->nnb[ids[i]]+=1;
             al->nnb[k]+=1;
             al->nb[ids[i]][al->nnb[ids[i]]-1] = k;
             al->nb[k][al->nnb[k]-1] = ids[i];
        }
      }
    }
  }
  
  else 
  {
    
    int m;
    real lambda, dvdlambda;
    lambda = dvdlambda = 0;
    
    t_nblist *nlist;
    matrix box;
    rvec box_space;
    rvec box_size;
    clear_rvec(box_size);
    clear_rvec(box_space);
    al2rvec(x,al);
    matrix mat2;
    max_coords(al,mat2);
    
    clear_mat(box);
    gen_box(0,al->natoms,x,box,box_space,FALSE);
    /*  write_pdb(al,"xx.pdb"); */
    for(m=0; (m<DIM); m++){
      box[m][m] += rlong+1;
      if(box[m][m] < 2.){
        box[m][m] = 2.;
      }
      box_size[m] = box[m][m];  
    }
    
    calc_shifts(box,fr->shift_vec); 
    put_charge_groups_in_box(NULL,0,cgs->nr,box,cgs,x,fr->cg_cm); 

/*     printf("nodeid = %d\n",cr->nodeid); */
    
    search_neighbours(log,fr,x,box,top,grps,cr,
                      nsb,nrnb,md,lambda,&dvdlambda,
                      TRUE,FALSE);


    nlist = &(fr->nblists[0].nlist_sr[eNL_VDW]);
    
    update_nl(al,nlist); 
  }
}

/*===================================================================*/




void max_crd(rvec *x, int n, matrix max)
{
  /* matrix X stores minimum
     matrix Z stores maximu
     this max[YY][XX] is minimum y value
  */

  int i,k;
  for(i=0;i<DIM;i++){
    max[i][XX] = 99999;
    max[i][YY] = 0.;
    max[i][ZZ] = -99999.;
  }
  for(i=0;i<n;i++){
    if(x[i][XX] < max[XX][XX]) max[XX][XX] = x[i][XX];
    if(x[i][XX] > max[XX][ZZ]) max[XX][ZZ] = x[i][XX];

    if(x[i][YY] < max[YY][XX]) max[YY][XX] = x[i][YY];
    if(x[i][YY] > max[YY][ZZ]) max[YY][ZZ] = x[i][YY];

    if(x[i][ZZ] < max[ZZ][XX]) max[ZZ][XX] = x[i][ZZ];
    if(x[i][ZZ] > max[ZZ][ZZ]) max[ZZ][ZZ] = x[i][ZZ];
  }
}

int gridp(real x, real origin, real inv_spacing, int max)
{
  real n;
  int point;
  n = (x-origin)*inv_spacing;
  point = (int) n;
  return point;
}

t_gridmap *nb(t_atomlist *al, t_gridmap *gp, real cutoff)
{
  int i,k;
  matrix max;
  int xdim, ydim, zdim;
  int gpx, gpy, gpz;
  int nx, ny, nz;
  int idx;
  if(al->natoms < 300) {
/*     printf("\tUsing non-grid-based simple neighborsearch for small systems\n"); */
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nnb[i];k++){
        al->nb[i][k] = 0;
      }
      al->nnb[i] = 0;
      al->memcount[i] = 1;
      srenew(al->nb[i],MAX_NL);
    }
    return NULL;
  }
  max_crd(al->x,al->natoms,max);
  gp->spacing = cutoff;
  gp->inv_spacing = 1./cutoff;
  
  for(i=0;i<DIM;i++){
    gp->origin[i] = max[i][XX];
  }
  xdim = (int) (max[XX][ZZ] - max[XX][XX])*gp->inv_spacing +1;
  ydim = (int) (max[YY][ZZ] - max[YY][XX])*gp->inv_spacing +1;
  zdim = (int) (max[ZZ][ZZ] - max[ZZ][XX])*gp->inv_spacing +1;
  
  gp->n[XX] = xdim ;
  gp->n[YY] = ydim ;
  gp->n[ZZ] = zdim ;
  nx = gp->n[XX];
  ny = gp->n[YY];
  nz = gp->n[ZZ];
  
  gp->nelem = gp->n[XX]*gp->n[YY]*gp->n[ZZ];
  srenew(gp->cell,gp->nelem);
  srenew(gp->natom,gp->nelem);
  for(i=0;i<gp->nelem;i++){
    snew(gp->cell[i],1);
    gp->natom[i] = 0;
  }
  for(i=0;i<al->natoms;i++){
    gpx = gridp(al->x[i][XX],gp->origin[XX],gp->inv_spacing, xdim);
    gpy = gridp(al->x[i][YY],gp->origin[YY],gp->inv_spacing, ydim);
    gpz = gridp(al->x[i][ZZ],gp->origin[ZZ],gp->inv_spacing, zdim);
    al->cell[i][XX] = gpx;
    al->cell[i][YY] = gpy;
    al->cell[i][ZZ] = gpz;
    idx = gpz*nx*ny + gpy*nx + gpx;
    gp->natom[idx]++;
    srenew(gp->cell[idx], gp->natom[idx]);
    gp->cell[idx][gp->natom[idx]-1] = i;
  }
/*   printf("ncells = %d\n",gp->nelem); */
/*   printf("x = %d y = %d z = %d\n",gp->n[XX],gp->n[YY],gp->n[ZZ]); */
  
/*   printf("Returning gridmap\n"); */
/*   for(i=0;i<al->natoms;i++){ */
/*     printf("cell = %d %d %d\n",al->cell[i][XX],al->cell[i][XX],al->cell[i][XX]); */
/*   } */
  
  return gp;
}

void free_gridmap(t_gridmap *gp)
{
  int i;
  for(i=0;i<gp->nelem;i++){
    sfree(gp->cell[i]);
  }
  sfree(gp->cell);
  sfree(gp->values);
  sfree(gp->natom);
  sfree(gp);
}



bool is_bonded(FILE *log, t_atomlist *al, int i, int k, real d2)
{
  char wstr[STRLEN];
  
  if(d2 < sqr(al->bcontr[i]+al->bcontr[k]) &&
     (strcmp(al->symbol[i],"H")!=0 || strcmp(al->symbol[k],"H") !=0)) {
    if(IsProtein(al,i) && IsProtein(al,k)) 
       /* check whether we have a suspicious bond */
      {
        if(al->resid[i] != al->resid[k] && !IS_OMEGA(al,i,k) && 
           !IS_SSBOND(al,i,k)) 
        {
          sprintf(wstr,"tCNC__wrn_> Encountered short distance between atom %d(%s)-%d%s and atom %d(%s)-%d%s\n",
                  al->id[i],al->name[i],al->resid[i],al->resname[i],
                  al->id[k],al->name[k],al->resid[k],al->resname[k]);
          CNClog(log,wstr);
          sprintf(wstr,"tCNC__wrn_> dist: %4.3g\n",sqrt(d2));
          CNClog(log,wstr);
          sprintf(wstr,"tCNC__wrn_> This will not be considered as bond!!\n");
          CNClog(log,wstr);
          return FALSE;
        }
        else 
        {
          return TRUE;
        }
      }
    else 
    {
      return TRUE;
    }
  }

  return FALSE;
}

void make_bonded(FILE *log, t_atomlist *al, int **nlist, int *natom, real **dlist, bool bUpdate)
{
  int i,k,j;
  int atom_id, atom2_id;
  real d2;
  if(bUpdate == MAKE_FULL_NEIGHBORLIST){
    
    for(i=0;i<al->natoms;i++){
      for(k=0;k<natom[i];k++){
        d2 = dlist[i][k];
        atom_id = nlist[i][k];
        if(is_bonded(log,al,i,atom_id,dlist[i][k])){
          add_bond(al, i, atom_id, NULL, BOND, BOND);
        }
      }
    }
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nbonds[i];k++){
        atom_id = al->bonds[i][k];
        for(j=0;j<al->nbonds[atom_id];j++){
          if(al->bonds[atom_id][j] < i){
            add_bond(al,i,al->bonds[atom_id][j],NULL,B13,B13);
          }
        }
      }
    }
    for(i=0;i<al->natoms;i++){
      for(k=0;k<al->nb13[i];k++){
        atom_id = al->b13[i][k];
        for(j=0;j<al->nbonds[atom_id];j++){
          atom2_id = al->bonds[atom_id][j];
          if(i<atom2_id && !al_connected(al,i,atom2_id)){
            add_bond(al,i,atom2_id,NULL,B14,B14);
          }
        }
      }
    }
  }
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      al->nb[i][k] = 0;
    }
    al->nnb[i] = 0;
  }
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<natom[i];k++){
      d2 = dlist[i][k];
      atom_id = nlist[i][k];
      if(i<atom_id && !al_connected(al,i,atom_id)){
        add_bond(al,i,atom_id,NULL,NL,NL);
      }
    }
  }
}

  

void fill_neighborlist(FILE *log, t_atomlist *al, t_gridmap *gp,  real cutoff, bool bUpdate)
{
  
  int i,k;
  int x, y, z;
  real d;
  rvec diff;
  real cut2 = sqr(cutoff);
  int idx;
  int nx, ny, nz;
  int natom;
  int atom_id;
  int gpx, gpy, gpz;
  int *nlist[al->natoms];
  real *dlist[al->natoms];
  int nnatoms[al->natoms];
  if(gp==NULL) {
    do_simple_search(al,NULL,cutoff,bUpdate);
  }
  else {
    
    nx = gp->n[XX];
    ny = gp->n[YY];
    nz = gp->n[ZZ];
    
    x = -1;
    y = -1;
    z = -1;

/*     gpx = al->cell[0][XX]; */
/*     gpy = al->cell[0][YY]; */
/*     gpz = al->cell[0][ZZ]; */
/*     idx = gpz * nx*ny + gpy*nx + gpx; */
/*     print_cell(gp->cell[idx],gp->natom[idx]); */
/*     exit(0); */
    


    
    for(i=0;i<al->natoms;i++){
      snew(nlist[i],1);
      snew(dlist[i],1);
      nnatoms[i] = 0;
    }
    for(i=0;i<al->natoms;i++){
      x = -1;
      while(x<2){
        y = -1;
        while(y<2){
          z = -1;
          while(z<2){
            gpx = x+al->cell[i][XX];
            gpy = y+al->cell[i][YY];
            gpz = z+al->cell[i][ZZ];
            if( (gpx < gp->n[XX]) && (gpy < gp->n[YY]) && (gpz < gp->n[ZZ])
                && (gpx >=0) && (gpy >= 0) && (gpz >= 0))
            {
              idx = gpz*nx*ny + gpy*nx + gpx;
              natom = gp->natom[idx];
              for(k=0;k<natom;k++){
                atom_id = gp->cell[idx][k];
                if(atom_id > i){
                  rvec_sub(al->x[i],al->x[atom_id],diff);
                  d = norm2(diff);
                  if(d < cut2){
                    nnatoms[i]+=1;
                    srenew(nlist[i],nnatoms[i]+1); 
                    srenew(dlist[i],nnatoms[i]+1); 
                    nlist[i][nnatoms[i]-1] = atom_id; 
                    dlist[i][nnatoms[i]-1] = d; 
                  }
                }
              }
            }
            z++;
          }
          y++;
        }
        x++;
      }
    }
    
    make_bonded(log, al, nlist, nnatoms, dlist, bUpdate); 
    for(i=0;i<al->natoms;i++){
      sfree(nlist[i]);
      sfree(dlist[i]);
    }
  }
  
}

                
