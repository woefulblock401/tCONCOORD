
#include <tconcoord.h>

/*============================================================*/
t_atomlist *atomlist_init(void)
/* initialise atomlist structure */
{
  t_atomlist *al = NULL;
  snew(al,1);
  al->natoms = 0;
  al->ngroups = 0;
  al->id = NULL;
  al->name = NULL;
  al->resname = NULL;
  al->resid = NULL;
  al->chain = NULL;
  al->x = NULL;
  al->f = NULL;
  al->occ = NULL;
  al->bfac = NULL;
  al->m = NULL;
  al->q = NULL;
  al->altloc = NULL;
  al->hyb = NULL;
  al->symbol = NULL;
  al->type = NULL;
  al->order = NULL;
  al->bonds = NULL;
  al->nbonds = NULL;
  al->b13 = NULL;
  al->nb13 = NULL;
  al->b14 = NULL;
  al->nb14 = NULL;
  al->nb = NULL;
  al->nnb = NULL;
  al->memcount = NULL;
  al->restype = NULL;
  al->ngrps = NULL;
  al->grpnr = NULL;
  al->isdon = NULL;
  al->isacc = NULL;
  al->ishphob = NULL;
  al->isposres = NULL;
  al->isflex = NULL;
  al->bcontr = NULL;
  al->vdw = NULL;
  al->vdw14 = NULL;
  al->vdwtab = NULL;
  al->vdw14tab = NULL;
  al->ptype = NULL;
  al->rs_type = NULL;
  al->rs_rad = NULL;
  al->rs_eps = NULL;
  al->rs_G = NULL;
  al->rs_Gref = NULL;
  al->rs_V = NULL;
  al->rs_lambda = NULL;
  al->rs_comb = NULL;
  al->cnc_solv = NULL;
  al->nrs_comb = 0;
  al->hbonds  = NULL;
  al->nhbonds = 0;
  al->cell = NULL;
  al->sasa_type = NULL;
  al->sasa_Ri = NULL;
  al->sasa_pi = NULL;
  al->sasa_sigi = NULL;
  return al;
}
/*============================================================*/
t_atomlist *al_realloc(t_atomlist *al,int natoms)
/* reallocate memory for atomlist */
{
  al->natoms = natoms;
  snew(al->id,natoms);
  snew(al->name,natoms);
  snew(al->altloc,natoms);
  snew(al->resname,natoms);
  snew(al->chain,natoms);
  snew(al->resid,natoms);
  snew(al->x,natoms);
  snew(al->f,natoms);
  snew(al->m,natoms);
  snew(al->occ,natoms);
  snew(al->hyb,natoms);
  snew(al->bfac,natoms);
  snew(al->symbol,natoms);
  snew(al->bonds,natoms);
  snew(al->nbonds,natoms);
  snew(al->b13,natoms);
  snew(al->nb13,natoms);
  snew(al->b14,natoms);
  snew(al->nb14,natoms);
  snew(al->nb,natoms);
  snew(al->nnb,natoms);
  snew(al->memcount,natoms);
  snew(al->restype,natoms);
  snew(al->ngrps,natoms);
  snew(al->grpnr,natoms);
  snew(al->bcontr,natoms);
  snew(al->vdw,natoms);
  snew(al->vdw14,natoms);
  snew(al->type,natoms);
  snew(al->order,natoms);
  snew(al->isdon,natoms);
  snew(al->isacc,natoms);
  snew(al->ishphob,natoms);
  snew(al->isposres,natoms);
  snew(al->isflex,natoms);
  snew(al->q,natoms);
  snew(al->ptype,natoms);
  snew(al->cell,natoms);
  int i;
  for(i=0;i<al->natoms;i++){
    al->nbonds[i] = 0;
    al->nb13[i] = 0;
    al->nb14[i] = 0;
    al->nnb[i] = 0;
    
    srenew(al->bonds[i],MAX_BOND);
    srenew(al->b13[i],MAX_B13);
    srenew(al->b14[i],MAX_B14);
    srenew(al->nb[i],MAX_NL);
    al->memcount[i] = 1;
  }
  snew(al->rs_type,natoms);
  snew(al->rs_rad,natoms);
  snew(al->rs_G,natoms);
  snew(al->rs_Gref,natoms);
  snew(al->rs_V,natoms);
  snew(al->rs_lambda,natoms);
  snew(al->rs_eps,natoms);
  snew(al->cnc_solv,natoms);
  snew(al->sasa_type,natoms);
  snew(al->sasa_Ri, natoms);
  snew(al->sasa_pi, natoms);
  snew(al->sasa_sigi, natoms);
  return al;
}
/*============================================================*/
void copy_atom(t_atomlist *al1, int i, t_atomlist *al2, int k)
{
  int n;
  al2->id[k] = al1->id[i];
  strcpy(al2->name[k],al1->name[i]);
  strcpy(al2->resname[k],al1->resname[i]);
  strcpy(al2->chain[k],al1->chain[i]);
  strcpy(al2->altloc[k],al1->altloc[i]);
  al2->resid[k] = al1->resid[i];
  copy_rvec(al1->x[i],al2->x[k]);
  al2->occ[k] = al1->occ[i];
  al2->bfac[k] = al1->bfac[i];
  strcpy(al2->type[k],al1->type[i]);
  copy_rvec(al1->f[i],al2->f[k]);
  strcpy(al2->symbol[k],al1->symbol[i]);
  strcpy(al2->hyb[k],al1->hyb[i]);
  al2->order[k] = al1->order[i];
  al2->nbonds[k] = al1->nbonds[i];
  al2->nb13[k] = al1->nb13[i];
  al2->nb14[k] = al1->nb14[i];
  al2->nnb[k] = al1->nnb[i];
  al2->restype[k] = al1->restype[i];
  al2->isdon[k] = al1->isdon[i];
  al2->isacc[k] = al1->isacc[i];
  al2->ishphob[k] = al1->ishphob[i];
  al2->isflex[k] = al1->isflex[i];
  al2->isposres[k] = al1->isposres[i];
  al2->bcontr[k] = al1->bcontr[i];
  al2->vdw[k] = al1->vdw[i];
  al2->vdw14[k] = al1->vdw14[i];
  al2->ptype[k] = al1->ptype[i];
  al2->m[k] = al1->m[i];
  al2->q[k] = al1->q[i];
  al2->rs_type[k] = al1->rs_type[i];
  al2->rs_rad[k] = al1->rs_type[i];
  al2->rs_G[k] = al1->rs_G[i];
  al2->rs_Gref[k] = al1->rs_Gref[i];
  al2->rs_V[k] = al1->rs_V[i];
  al2->rs_eps[k] = al1->rs_eps[i];
  al2->rs_lambda[k] = al1->rs_lambda[i];
  
  /* bonds and stuff */
  al2->nbonds[k] = al1->nbonds[i];
  al2->nb13[k] = al1->nb13[i];
  al2->nb14[k] = al1->nb14[i];
  al2->nnb[k] = al1->nnb[i];
 
  for(n=0;n<al1->nbonds[i];n++){
    al2->bonds[k][n] = al1->bonds[i][n];
  }
  for(n=0;n<al1->nb13[i];n++){
    al2->b13[k][n] = al1->b13[i][n];
  }
  for(n=0;n<al1->nb14[i];n++){
    al2->b14[k][n] = al1->b14[i][n];
  }
  for(n=0;n<al1->nnb[i];n++){
    al2->nb[k][n] = al1->nb[i][n];
  }

  

 
}
/*============================================================*/
void copy_al(t_atomlist *src, t_atomlist *dest)
{
  int i,k;
/*   dest->natoms = src->natoms; */
  dest->nvdw = src->nvdw;
  /* malloc vdwtab */
  snew(dest->vdwtab,dest->nvdw);
  for(i=0;i<dest->nvdw;i++){
    snew(dest->vdwtab[i],dest->nvdw);
  }
  /* copy vdwtab */
  for(i=0;i<dest->nvdw;i++){
    for(k=0;k<dest->nvdw;k++){
      dest->vdwtab[i][k] = src->vdwtab[i][k];
    }
  }
  for(i=0;i<src->natoms;i++){
    copy_atom(src,i,dest,i);
  }
}

/*============================================================*/

void free_al(t_atomlist *al)

{
  int i;
  sfree(al->id);
  sfree(al->resid);
  sfree(al->name);
  sfree(al->resname);
  sfree(al->chain);
  sfree(al->type);
  sfree(al->altloc);
  sfree(al->x);
  sfree(al->f);
  sfree(al->occ);
  sfree(al->bfac);
  sfree(al->symbol);
  sfree(al->hyb);
  sfree(al->order);

  if(al!=NULL){
    for(i=0;i<al->natoms;i++){
      sfree(al->bonds[i]);
      sfree(al->b13[i]);
      sfree(al->b14[i]);
      sfree(al->nb[i]);
    }
  }
  sfree(al->nbonds);
  sfree(al->bonds);
  sfree(al->nb13);
  sfree(al->b13);
  sfree(al->nb14);
  sfree(al->b14);
  sfree(al->nnb);
  sfree(al->nb);
  
  sfree(al->restype);
  sfree(al->ngrps);
  for(i=0;i<al->natoms;i++){
    sfree(al->grpnr[i]);
  }
  sfree(al->grpnr);

  sfree(al->isdon);
  sfree(al->isacc);
  sfree(al->ishphob);
  sfree(al->isposres);
  sfree(al->isflex);
  sfree(al->bcontr);
  sfree(al->vdw);
  sfree(al->vdw14);
  
  sfree(al->q);
  sfree(al->m);
  sfree(al->ptype);
  for(i=0;i<al->nvdw;i++){
    sfree(al->vdwtab[i]);
    sfree(al->vdw14tab[i]);
  }
  sfree(al->vdwtab);
  sfree(al->vdw14tab);

  sfree(al->rs_type);
  sfree(al->rs_rad);
  sfree(al->rs_G);
  sfree(al->rs_Gref);
  sfree(al->rs_V);
  sfree(al->rs_lambda);
  sfree(al->rs_eps);
  for(i=0;i<al->nrs_comb;i++){
    sfree(al->rs_comb[i]);
  }
  sfree(al->rs_comb);
  sfree(al->cnc_solv);

  sfree(al);
}

/*============================================================*/

void get_restype(t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    if(IsProtein(al,i)) al->restype[i] = rtPROTEIN;
    else if(IsNucac(al,i)) al->restype[i] = rtNUC; 
     else if(IsIon(al,i)) al->restype[i] = rtION;  
    else if(IsSol(al,i)) al->restype[i] = rtSOL; 
    else al->restype[i] = rtOTH;
  }
}


/*============================================================*/
int restype(char *resname)
{
  int i;
  int ret = rtOTH;

  for(i=0;i<asize(prot_res);i++){
    if(strcmp(resname,prot_res[i])==0){
      ret = rtPROTEIN;
      return ret;
    }
  }
  for(i=0;i<asize(nucac_res);i++){
    if(strcmp(resname,nucac_res[i])==0){
      ret = rtNUC;
      return ret;
    }
  }
  for(i=0;i<asize(ion_res);i++){
    if(strcmp(resname,ion_res[i])==0){
      ret = rtION;
      return ret;
    }
  }
  for(i=0;i<asize(sol_res);i++){
    if(strcmp(resname,sol_res[i])==0){
      ret = rtSOL;
      return ret;
    }
  }
  return ret;
}  
/*============================================================*/
t_atoms *al2atoms(t_atoms *atoms,t_atomlist *al, rvec **x)
/* build t_atoms structure from atomlist */
{
  int i;
  snew(atoms,1);
  init_t_atoms(atoms,al->natoms,TRUE);
  atoms->nr = al->natoms;
  snew(*x,al->natoms);
  for(i=0;i<al->natoms;i++){ 
    al2rvec((*x),al); 
    atoms->atom[i].chain = al->chain[i][0]; 
    snew(atoms->atomname[i],1); 
    *atoms->atomname[i]=al->name[i];   
  } 
  return atoms;
}

/*============================================================*/
void rvec2al(t_atomlist *al,rvec *x)
/* copy rvec to atomlist */
{
  int i,k;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      al->x[i][k] = x[i][k]*10.;
    }
  }
}
/*============================================================*/
void al2rvec(rvec *x,t_atomlist *al)
/* copy al->x to rvec */
{
  int i,k;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      x[i][k] = al->x[i][k]*0.1;
    }
  }
}
/*============================================================*/
void get_symbol(FILE *log,t_atomlist *al)
{
  char warn[STRLEN];
  char dum[STRLEN];
  
  int i;
  for(i=0;i<al->natoms;i++){
    strcpy(al->symbol[i],"");
    strcpy(dum,al->name[i]);
    trim(dum);
    switch(dum[0])
    {
      case('H'):
        strcpy(al->symbol[i],"H");
        break;
      case('\''):
        if(al->name[i][1]=='H') {
          strcpy(al->symbol[i],"H");
        }
        break;
      case('C'):
        if (IsProtein(al,i) || IsNucac(al,i))
          strcpy(al->symbol[i],"C");
        else if(IsIon(al,i) && al->name[i][2]=='A')
          strcpy(al->symbol[i],"CA");
        else {
          if(al->name[i][2]=='L') {
            strcpy(al->symbol[i],"CL");
          }
          else {
            strcpy(al->symbol[i],"C");
          }
        }
        break;
      case('N'):
        if(IsIon(al,i) && al->name[i][2]=='I')
          strcpy(al->symbol[i],"NI");
        else strcpy(al->symbol[i],"N");
        break;
      case('O'):
        strcpy(al->symbol[i],"O");
        break;
      case('S'):
        strcpy(al->symbol[i],"S");
        break;
      case('F'):
        if(al->name[i][2]=='E')
          strcpy(al->symbol[i],"FE");
        else
          strcpy(al->symbol[i],"F");
        break;
      case('B'):
        if(al->name[i][2]=='R')
          strcpy(al->symbol[i],"BR");
        break;
      case('P'):
        strcpy(al->symbol[i],"P");
        break;
      case('M'):
        if(al->name[i][2]=='N')
          strcpy(al->symbol[i],"MN");
        else if(al->name[i][2]=='G')
          strcpy(al->symbol[i],"MG");
        else if(al->name[i][2]=='W')
          strcpy(al->symbol[i],"DM");
        break;
      case('Z'):
        if(al->name[i][2]=='N')
          strcpy(al->symbol[i],"ZN");
        break;
      case('I'):
        strcpy(al->symbol[i],"I");
        break;
      case('D'):
        strcpy(al->symbol[i],"D");
        break;

    }
    if(strcmp(al->symbol[i],"") == 0)
    {
      sprintf(warn,"Could not determine element type for %s(%s)\n",
              al->name[i],al->resname[i]);
      CNCerr(log,warn);
    }
  }
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->symbol[i],"H")==0 ||
       strcmp(al->symbol[i],"DM") == 0){
        al->bcontr[i]=HCONTR;
        al->m[i] = 1.;
    }
    else if(strcmp(al->symbol[i],"C")==0) {
      al->bcontr[i]=CCONTR;
    }
    else if(strcmp(al->symbol[i],"O")==0) {
      al->bcontr[i]=OCONTR;
    }
    else if(strcmp(al->symbol[i],"N")==0) {
      al->bcontr[i]=NCONTR;
    }
    else if(strcmp(al->symbol[i],"S")==0) {
      al->bcontr[i]=SCONTR;
    }
    else if(strcmp(al->symbol[i],"P")==0) {
      al->bcontr[i]=PCONTR;
    }
    else if(strcmp(al->symbol[i],"I")==0) {
      al->bcontr[i]=ICONTR;
    }
    else if(strcmp(al->symbol[i],"F")==0) {
      al->bcontr[i]=FCONTR;
    }
    else if(strcmp(al->symbol[i],"BR")==0) {
      al->bcontr[i]=BRCONTR;
    }
    else if(strcmp(al->symbol[i],"CL")==0) {
      al->bcontr[i]=CLCONTR;
    }
    else if(IsIon(al,i)) {
      al->bcontr[i] = MCONTR;
    }
    else if(strcmp(al->symbol[i],"D")==0) {
      switch(al->name[i][1]){
        case 'C':
          al->bcontr[i]=CCONTR;
          break;
        case 'N':
          al->bcontr[i]=NCONTR;
          break;
        case 'H':
          al->bcontr[i]=HCONTR;
          break;
        case 'O':
          al->bcontr[i]=OCONTR;
          break;
        case 'S':
          al->bcontr[i]=SCONTR;
          break;
        default:
          al->bcontr[i] = 1.3;
      }
    }
    

  }
  get_restype(al);   
  for(i=0;i<al->natoms;i++){
    if(al->m[i] < 1.0){
      sprintf(warn,"Unknown element %s %g\n",al->name[i],al->m[i]); 
      CNCwarn(log,warn);
    }
  }
  
}


/*============================================================*/
void get_symbol2(FILE *log,t_atomlist *al)
/* get element symbol, mass and bond contribution */
{
  char warn[STRLEN];
  char dum[STRLEN];
  
  int i;
  for(i=0;i<al->natoms;i++){
    strcpy(al->symbol[i],"");
    printf("\'%s\'\n",al->name[i]);
    strcpy(dum,al->name[i]);
    trim(dum);
    switch(dum[0]) {
      case('1'):
      case('2'):
      case('3'):
      case('H'):
        strcpy(al->symbol[i],"H");
        break;
      case('C'):
        if (IsProtein(al,i) || IsNucac(al,i))
          strcpy(al->symbol[i],"C");
        else if(IsIon(al,i)){
          if (dum[1]=='A')
            strcpy(al->symbol[i],"CA");
          else if (dum[1]=='L' || dum[1]=='l')
            strcpy(al->symbol[i],"CL");
        }
        else strcpy(al->symbol[i],"C");
        break;
      case('F'):
        if(dum[1]=='E' || dum[1]=='e')
          strcpy(al->symbol[i],"FE");
        else
          strcpy(al->symbol[i],"F");
        break;
      case('B'):
        if(dum[1]=='R' || dum[1]=='r')
          strcpy(al->symbol[i],"BR");
        else
          strcpy(al->symbol[i],"B");
        break;
      case('M'):
        if(dum[1]=='N' || dum[1]=='n')
          strcpy(al->symbol[i],"MN");
        else if(dum[1]=='G' || dum[1]=='g')
          strcpy(al->symbol[i],"MG");
        break;
      case('Z'):
        if(dum[1]=='N' || dum[1]=='n')
          strcpy(al->symbol[i],"ZN");
        break;
      case('I'):
        strcpy(al->symbol[i],"I");
        break;
    }
    if(strcmp(al->symbol[i],"") == 0)
    {
      al->symbol[i][0]=al->name[i][1];
      al->symbol[i][1]=0;
    }
    if(IsIon(al,i)){
      shstr dumm;
      sscanf(al->name[i],"%s",dumm);
      strcpy(al->symbol[i],dumm);
    }
    
  }
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->symbol[i],"")!=0){
      if(strcmp(al->symbol[i],"H")==0) {
        al->bcontr[i]=HCONTR;
        al->m[i] = 1.; 
      }
      else if(strcmp(al->symbol[i],"C")==0) {
        al->bcontr[i]=CCONTR;
        al->m[i] = 12.;
      }
      else if(strcmp(al->symbol[i],"O")==0) {
        al->bcontr[i]=OCONTR;
        al->m[i] = 16.;
      }
      else if(strcmp(al->symbol[i],"N")==0) {
        al->bcontr[i]=NCONTR;
        al->m[i] = 14.;
      }
      else if(strcmp(al->symbol[i],"S")==0) {
        al->bcontr[i]=SCONTR;
        al->m[i] = 32.;
      }
      else if(strcmp(al->symbol[i],"P")==0) {
        al->bcontr[i]=PCONTR;
        al->m[i] = 31.;
      }
      else if(strcmp(al->symbol[i],"I")==0) {
        al->bcontr[i]=ICONTR;
        al->m[i] = 127.;
      }
      else if(strcmp(al->symbol[i],"F")==0) {
        al->bcontr[i]=FCONTR;
        al->m[i] = 19.;
      }
      else if(strcmp(al->symbol[i],"BR")==0) {
        al->bcontr[i]=BRCONTR;
        al->m[i] = 80.;
      }
      else if(strcmp(al->symbol[i],"CL")==0) {
        al->bcontr[i]=CLCONTR;
        al->m[i] = 35.5;
      }
      else if(IsIon(al,i)) {
        al->bcontr[i] = MCONTR;
        al->m[i] = 25.;
      }
    }
  }
  get_restype(al);   
  for(i=0;i<al->natoms;i++){
    if(al->m[i] < 1.0){
      al->m[i] = 1.;
      sprintf(warn,"Unknown element %s\n",al->name[i]); 
      CNCwarn(log,warn);
      sprintf(warn,"Setting mass of %s(%s) to 1.0\n",al->name[i],al->resname[i]);
      CNCerr(log,warn);
    }
  }
}
/*============================================================*/
void get_bconstr(t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->symbol[i],"")!=0){
      if(strcmp(al->symbol[i],"H")==0) {
        al->bcontr[i]=HCONTR;
      }
      else if(strcmp(al->symbol[i],"C")==0) {
        al->bcontr[i]=CCONTR;
      }
      else if(strcmp(al->symbol[i],"O")==0) {
        al->bcontr[i]=OCONTR;
      }
      else if(strcmp(al->symbol[i],"N")==0) {
        al->bcontr[i]=NCONTR;
      }
      else if(strcmp(al->symbol[i],"S")==0) {
        al->bcontr[i]=SCONTR;
      }
      else if(strcmp(al->symbol[i],"P")==0) {
        al->bcontr[i]=PCONTR;
      }
      else if(strcmp(al->symbol[i],"I")==0) {
        al->bcontr[i]=ICONTR;
      }
      else if(strcmp(al->symbol[i],"F")==0) {
        al->bcontr[i]=FCONTR;
      }
      
      else if(IsIon(al,i)) {
        al->bcontr[i] = MCONTR;
      }
    }
  }
}

/*============================================================*/

void com(t_atomlist *al)
/* translate atomlist coordinates to the center of mass */
{
  int i,k;
  rvec com;
  clear_rvec(com);
  real mtot = 0;
  for(i=0;i<al->natoms;i++) {
    com[XX]+= al->x[i][XX]*al->m[i];
    com[YY]+= al->x[i][YY]*al->m[i];
    com[ZZ]+= al->x[i][ZZ]*al->m[i];
    mtot += al->m[i];
  }
  com[XX] /=mtot;
  com[YY] /=mtot;
  com[ZZ] /=mtot;
  for (i=0; i<al->natoms; i++){
    for(k=0;k<DIM;k++){
      al->x[i][k]-=com[k];
    }
  }
}
/*============================================================*/
void calc_cent(t_atomlist *al, rvec com)
/* translate atomlist coordinates to the center of mass */
{
  int i,k;
  clear_rvec(com);
  real mtot = 0;
  for(i=0;i<al->natoms;i++) {
    com[XX]+= al->x[i][XX]*al->m[i];
    com[YY]+= al->x[i][YY]*al->m[i];
    com[ZZ]+= al->x[i][ZZ]*al->m[i];
    mtot += al->m[i];
  }
  com[XX] /=mtot;
  com[YY] /=mtot;
  com[ZZ] /=mtot;
}
/*============================================================*/
void com_frag(t_atomlist *al, int n,int *idx,rvec *new,rvec cm)
/* calculate the center of a part of the atomlist */
{
  int i,k;
  real mtot = 0;
  clear_rvec(cm);
  for(i=0;i<n;i++) {
    cm[XX]+= al->x[idx[i]][XX]*al->m[idx[i]];
    cm[YY]+= al->x[idx[i]][YY]*al->m[idx[i]];
    cm[ZZ]+= al->x[idx[i]][ZZ]*al->m[idx[i]];
    mtot += al->m[idx[i]];
  }
  cm[XX] /=mtot;
  cm[YY] /=mtot;
  cm[ZZ] /=mtot;
  for (i=0; i<n; i++){
    for(k=0;k<DIM;k++){
      new[i][k] = al->x[idx[i]][k]-cm[k];
    }
  }
}
  
/*============================================================*/

bool al_connected(t_atomlist *al, int i, int k)
{
  /* check whether atom k is 
     1-2, 1-3 or 1-4 connected 
     to atom i */

  int j;
  if (i == k) return TRUE;
  
  for(j=0;j<al->nbonds[i];j++){
    if(al->bonds[i][j] == k)
      return TRUE;
  }
  for(j=0;j<al->nb13[i];j++){
    if(al->b13[i][j] == k)
      return TRUE;
  }
  for(j=0;j<al->nb14[i];j++){
    if(al->b14[i][j] == k)
      return TRUE;
  }
  return FALSE;
}

/*============================================================*/  

int count_res(t_atomlist *al)
{
  int i;
  int count = 1;
  for(i=1;i<al->natoms;i++){
    if(al->resid[i-1] != al->resid[i]) count++;
  }
  return count;
}


/*============================================================*/
real get_vdw(t_atomlist *al, int i, int k, int flag)
{
  /* look up vdw table or vdw14 table if flag is TRUE */

  if(!flag)
    return al->vdwtab[al->ptype[i]][al->ptype[k]];
  else
    return al->vdw14tab[al->ptype[i]][al->ptype[k]];
}
/*============================================================*/
void print_atom(FILE *fp,t_atomlist *al,int i)
{
  /* write atom as pdb string */

  char chain[STRLEN];
  char res[STRLEN];
  if(strlen(al->resname[i]) > 3){
    strncpy(res,al->resname[i],3);
    res[4] = '\0';
  }
  else strcpy(res,al->resname[i]);
  
  if(strcmp(al->chain[i],"0") == 0) strcpy(chain," ");
  else strcpy(chain,al->chain[i]);
  fprintf(fp,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f%8.3f%8.3f %5.2f %5.2f\n",
          al->id[i],al->name[i],al->altloc[i],res,chain,
          al->resid[i],al->x[i][0],al->x[i][1],al->x[i][2],
          al->occ[i],al->bfac[i]);
}

/*============================================================*/
void write_pdb(t_atomlist *al,char *filename)
{
  /* write pdb file */

  int i;
  FILE *fp = ffopen(filename,"w");
  fprintf(fp,"HEADER    FILE GENERATED BY tCONCOORD\n");
  fprintf(fp,"MODEL 1\n");
  for(i=0;i<al->natoms;i++){
    print_atom(fp,al,i);
  }
  fprintf(fp,"ENDMDL\n");
  fclose(fp);
}
/*============================================================*/
void write_pdb_frame(FILE *fp,t_atomlist *al, int n)
{
  /* write pdb trajectory frame */

  int i;
  fprintf(fp,"MODEL%5d\n",n);
  for(i=0;i<al->natoms;i++){
    print_atom(fp,al,i);
  }
  fprintf(fp,"ENDMDL\n");
  fflush(fp);
}
/*============================================================*/
void al2xtc(t_atomlist *al, int xtc, int frame)
{
  /* write coordinates as to xtc file */

  int natoms = al->natoms;
  real time = frame;
  matrix box;
  clear_mat(box);
  real prec = 1000.;
  rvec *x;
  snew(x,al->natoms);
  al2rvec(x,al);
  if(!write_xtc(xtc,natoms,frame,time,box,x,prec))
    fatal_error("writing xtc sucks....\n");
  sfree(x);
}
/*============================================================*/
void renumber_atoms(t_atomlist *al, int start)
{
  /* does pretty much what is says */

  int i;
  for(i=0;i<al->natoms;i++){
    al->id[i] = start;
    start++;
  }
}
/*============================================================*/
void renumber_residues(t_atomlist *al, int start)
{
  /* does pretty much what is says */

  int i,k;
  int beg = 0;
  for(i=0;i<al->natoms;i++)
  {
    if(al->resid[i] != al->resid[i+1] ||
       i+1 == al->natoms){
      for(k=beg;k<i+1;k++){
        al->resid[k] = start;
      }
      start++;
      beg = i+1;
    }
  }
}


/*============================================================*/
void rename_at(t_atomlist *al)
{
  
  /* here we rename atoms to assign
     the correct parameters in
     the lib files.
     The element symbol must be in
     pos 1 of the name array.
  */
  
  char dum;
  int i,k;
  
  for(i=0;i<al->natoms;i++){
    if(IsProtein(al,i)){
      /* do c-terminus here */
      if(strcmp(al->name[i]," O1 ")==0) strcpy(al->name[i]," O  ");
      else if(strcmp(al->name[i]," O2 ") ==0 ) strcpy(al->name[i]," OXT");
      else if(strcmp(al->name[i]," O\' ") ==0 ) strcpy(al->name[i]," O  ");
      else if(strcmp(al->name[i]," O\'\'") ==0 ) strcpy(al->name[i]," OXT");
      /* gromacs ILE */
      if(strcmp(al->name[i]," CD ") ==0 && strcmp(al->resname[i],"ILE")==0)
        strcpy(al->name[i]," CD1");
      if(strcmp(al->name[i]," HD1") ==0 && strcmp(al->resname[i],"ILE")==0)
        strcpy(al->name[i],"1HD1");
      if(strcmp(al->name[i]," HD2") ==0 && strcmp(al->resname[i],"ILE")==0)
        strcpy(al->name[i],"2HD1");
      if(strcmp(al->name[i]," HD3") ==0 && strcmp(al->resname[i],"ILE")==0)
        strcpy(al->name[i],"3HD1");
    
    
    /* check for hydrogens */
    trim(al->name[i]);
    if(strlen(al->name[i]) == 4 &&
       isalpha(al->name[i][0])){
      dum = al->name[i][3];
      al->name[i][3] = al->name[i][2];
      al->name[i][2] = al->name[i][1];
      al->name[i][1] = al->name[i][0];
      al->name[i][0] = dum;
    }
    else if(strlen(al->name[i]) == 3){
      if(al->name[i][0] == 'H'){
        for(k=0;k<al->nb13[i];k++){
          if(strcmp(al->symbol[al->b13[i][k]],"H") == 0){
            extend_name(al->name[i]);
            al->name[i][0] = al->name[i][3];
            al->name[i][3] = ' ';
            break;
          }
        }
      }
    }
    else if(strlen(al->name[i]) == 2){
      if(al->name[i][0] == 'H' && isdigit(al->name[i][1])){
        for(k=0;k<al->nb13[i];k++){
          if(strcmp(al->symbol[al->b13[i][k]],"H") == 0){
/*             printf("name %s resn %s\n",al->name[i], al->resname[i]); */
            extend_name(al->name[i]);
            al->name[i][0] = al->name[i][2];
            al->name[i][2] = ' ';
            break;
          }
        }
      }
    }
    extend_name(al->name[i]);
    }
  }
  
}



/*============================================================*/
void get_order(t_atomlist *al)
{
  /* here we determine the number of bonds
     the particular atom is away from
     the backbone. Correct assignment is needed
     for swapping sidechains in tdisco.
  */

  int i;
  for(i=0;i<al->natoms;i++){
    if(IsProtein(al,i)){
      switch(al->name[i][2]){
        case(' '):
        case('A'):
          if(strcmp(al->symbol[i],"H")!= 0 &&
             strcmp(al->symbol[i],"O")!=0 ) al->order[i] = 0;
          else al->order[i]=1;
          break;
        case('B'):
          if(strcmp(al->symbol[i],"H")!= 0) al->order[i] = 1;
          else al->order[i]=2;
          break;
        case('G'):
          if(strcmp(al->symbol[i],"H")!= 0) al->order[i] = 2;
          else al->order[i]=3;
          break;
        case('D'):
          if(strcmp(al->symbol[i],"H")!= 0) al->order[i] = 3;
          else al->order[i]=4;
          break;
        case('E'):
          if(strcmp(al->symbol[i],"H")!= 0) al->order[i] = 4;
          else al->order[i]=5;
          break;
        case('Z'):
          if(strcmp(al->symbol[i],"H")!= 0) al->order[i] = 5;
          else al->order[i]=6;
          break;
        case('H'):
          if(strcmp(al->symbol[i],"H")!=0) al->order[i] = 6;
          else al->order[i]=7;
          break;
        case('X'): /* c-terminus */
          al->order[i]=1;
          break;
      }
    }
    
  }
}
/*============================================================*/
void occ_to_one(t_atomlist *al)
{
  /* set the occupancy value to 1. */

  int i;
  for(i=0;i<al->natoms;i++){
    al->occ[i] = 1.;
  }
}

/*============================================================*/

t_atomlist *al_from_atoms(t_atoms *atoms,
                          rvec *x)
{
  /* build t_atomlist from t_atoms */

  int i,nat;
  t_atomlist *al = atomlist_init();
  nat = atoms->nr;
  
  al->natoms = nat;
  al = al_realloc(al,nat); 
  rvec2al(al,x); 

  for(i=0;i<nat;i++){
    strcpy(al->name[i],*(atoms->atomname[i]));
    al->id[i] = i+1;
    al->m[i] = atoms->atom[i].m;
    al->q[i] = atoms->atom[i].q;
    al->resid[i] = atoms->atom[i].resnr+1;
    strcpy(al->resname[i],*(atoms->resname[atoms->atom[i].resnr]));
    extend_name(al->name[i]);
    al->chain[i][0] = atoms->atom[i].chain;
    al->chain[i][1] = '\0';
    al->occ[i] = 1.;
    al->bfac[i] = 0.;
  }
  return al;
}

/*============================================================*/

void bonds_from_idef(t_atomlist *al, t_idef *idef)
{
  /* take bonds, angles and dihedrals from
     idef structure
  */

  int i,k,na;
  int at1,at2,at3,at4;
  t_iatom *ia = NULL; 
  

  t_contab *ct = contab_init();
  ct = contab_realloc(ct,al->natoms);
  

  if(idef->ntypes!=-1){
    /* read bonds */
    for(i=0;i<F_NRE;i++){
      if(IS_CHEMBOND(i)  && idef->il[i].nr ){
        int nratoms=interaction_function[i].nratoms;
        k=0;
        while(k<idef->il[i].nr){
          at1 = idef->il[i].iatoms[k+1];
          at2 = idef->il[i].iatoms[k+2];
          add_bond(al,at1,at2,ct,BOND,BOND); 
          k+=nratoms+1;
        }
      }
    }
    /* read angles */
    na = idef->il[F_ANGLES].nr; 
    ia = idef->il[F_ANGLES].iatoms;
    for(i=0;i<na;i+=4){
      at1 = ia[i+1];
      at2 = ia[i+2];
      at3 = ia[i+3];
      add_bond(al,at1,at3,ct,B13,B13);
    }

    /* Ryckaert-Bell. Dihedrals : 
       works for OPLSAA only
    */

    na = idef->il[F_RBDIHS].nr; 
    ia = idef->il[F_RBDIHS].iatoms;
    for(i=0;i<na;i+=5){
       at1 = ia[i+1]; 
       at2 = ia[i+2]; 
       at3 = ia[i+3]; 
       at4 = ia[i+4];
       add_bond(al,at1,at4,ct,B14,B14); 
    }
  }
}
/*============================================================*/
void copy_al2(t_atomlist *src, t_atomlist *dest)
{
  int i,k;
/*   dest->natoms = src->natoms; */
  dest->nvdw = src->nvdw;
  /* malloc vdwtab */
  snew(dest->vdwtab,dest->nvdw);
  for(i=0;i<dest->nvdw;i++){
    snew(dest->vdwtab[i],dest->nvdw);
  }
  /* copy vdwtab */
  for(i=0;i<dest->nvdw;i++){
    for(k=0;k<dest->nvdw;k++){
      dest->vdwtab[i][k] = src->vdwtab[i][k];
    }
  }
  for(i=0;i<src->natoms;i++){
    copy_atom(src,i,dest,i);
  }
}

t_atomlist *fuse_al(t_atomlist *al1, t_atomlist *al2)
{
  int i,k;
  int natoms;
  natoms = al1->natoms+al2->natoms;
  t_atomlist *al = atomlist_init();
  al = al_realloc(al,natoms);
  copy_al2(al1,al);
  for(i=0;i<al2->natoms;i++){
    copy_atom(al2,i,al,i+al1->natoms);
  }
  al->natoms = natoms;

  /* assign new id's to bonds/b13 ... */

  for(i=al1->natoms;i<natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      al->bonds[i][k]+=al1->natoms;
    }
    for(k=0;k<al->nb13[i];k++){
      al->b13[i][k]+=al1->natoms;
    }
    for(k=0;k<al->nb14[i];k++){
      al->b13[i][k]+=al1->natoms;
    }
  }
  return al;
}
