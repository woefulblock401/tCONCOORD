#include <tconcoord.h>
/*================================================================*/
void input_error(int n, char *line)
{
  char error[STRLEN];
  sprintf(error,"Error while reading bounds: line %d -> %s\n",n,line);
  CNCerr(stderr,error);
}

/*================================================================*/

t_bounds *bounds_init(int n)
{
  t_bounds *b = NULL;
  snew(b,1);
  b->n = 0;
  b->at1 = b->at2 = b->at3 = b->at4 = NULL;
  b->av = b->lb = b->ub = b->ang = b->dih = b->sig = NULL;
  b->isang = b->isdih = b->isbond = NULL;
  if(n){
    snew(b->at1,n);
    snew(b->at2,n);
    snew(b->at3,n);
    snew(b->at4,n);
    snew(b->av,n);
    snew(b->lb,n);
    snew(b->ub,n);
    snew(b->ang,n);
    snew(b->sig,n);
    snew(b->dih,n);
    snew(b->isbond,n);
    snew(b->isang,n);
    snew(b->isdih,n);
    snew(b->bCheck,n);
    snew(b->pald,n);
    snew(b->pdla,n);
    snew(b->phplhp,n);
    snew(b->ac_ac,n);
    snew(b->don_ac,n);
    snew(b->don_don,n);
    
    int i;
    for(i=0;i<n;i++){
      b->at1[i] = 0;
      b->at2[i] = 0;
      b->at3[i] = 0;
      b->at4[i] = 0;
      b->av[i] = 0;
      b->lb[i] = 0;
      b->ub[i] = 0;
      b->ang[i] = 0;
      b->sig[i] = 0;
      b->dih[i] = 0;
      b->isbond[i] = FALSE;
      b->isang[i] = FALSE;
      b->isdih[i] = FALSE;
      b->bCheck[i] = TRUE;
      b->pald[i] = FALSE;
      b->pdla[i] = FALSE;
      b->phplhp[i] = FALSE;
      b->ac_ac[i] = FALSE;
      b->don_ac[i] = FALSE;
      b->don_don[i] = FALSE;
    }
  }
  return b;
}
/*================================================================*/
t_bounds *bounds_realloc(t_bounds *b, int n)
{
  srenew(b->at1,n);
  b->at1[n-1] = 0;
  srenew(b->at2,n);
  b->at2[n-1] = 0;
  srenew(b->at3,n);
  b->at3[n-1] = 0;
  srenew(b->at4,n);
  b->at4[n-1] = 0;
  srenew(b->av,n);
  b->av[n-1] = 0;
  srenew(b->lb,n);
  b->lb[n-1] = 0;
  srenew(b->ub,n);
  b->ub[n-1] = 0;
  srenew(b->ang,n);
  b->ang[n-1] = 0;
  srenew(b->dih,n);
  b->dih[n-1] = 0;
  srenew(b->sig,n);
  b->sig[n-1] = 0;
  srenew(b->isang,n);
  b->isang[n-1] = FALSE;
  srenew(b->isdih,n);
  b->isdih[n-1] = FALSE;
  srenew(b->isbond,n);
  b->isbond[n-1] = FALSE;
  srenew(b->bCheck,n);
  b->bCheck[n-1] = TRUE;
  srenew(b->pdla,n);
  b->pdla[n-1] = FALSE;
  srenew(b->pald,n);
  b->pald[n-1] = FALSE;
  srenew(b->phplhp,n);
  b->phplhp[n-1] = FALSE;
  srenew(b->ac_ac,n);
  b->ac_ac[n-1] = FALSE;
  srenew(b->don_ac,n);
  b->don_ac[n-1] = FALSE;
  srenew(b->don_don,n);
  b->don_don[n-1] = FALSE;
  return b;
}
/*================================================================*/
int count_bounds(FILE *fp)
{
  int count = 0;
  char line[STRLEN];
  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')==NULL) count++;
  }
  return count;
}
  

/*================================================================*/
t_bounds *read_bounds(char *filename, bool bVerbose)
{
  FILE *fp = ffopen(filename,"r");

  char line[STRLEN];
  int i = 0;
  bool isbond = FALSE;
  bool isang = FALSE;
  bool isdih = FALSE;
  bool pald = FALSE;
  bool pdla = FALSE;
  bool phplhp = FALSE;
  
  int dum;
  int nbounds;
  nbounds = count_bounds(fp);
  t_bounds *b = bounds_init(nbounds);
  rewind(fp);
/*   if(!bVerbose) */
/*     fprintf(stderr,"Reading bounds.......\n"); */
  
  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      if(strcmp(line,"[ BONDS ]") == 0){
        isbond = TRUE;
        isang = FALSE;
        isdih = FALSE;
        pald = FALSE;
        pdla = FALSE;
        phplhp = FALSE;
        
      }
      else if(strcmp(line,"[ ANGLES ]") == 0){
        isbond = FALSE;
        isang = TRUE;
        isdih = FALSE;
        pald = FALSE;
        pdla = FALSE;
        phplhp = FALSE;
        }
      else if(strcmp(line,"[ 1-4 PAIRS ]") == 0){
        isbond = FALSE;
        isang = FALSE;
        isdih = TRUE;
        pald = FALSE;
        pdla = FALSE;
        phplhp = FALSE;
      }

/*       else if(strcmp(line,"# pa/ld") == 0){ */
/*         isbond = FALSE; */
/*         isang = FALSE; */
/*         isdih = FALSE; */
/*         pald = TRUE; */
/*         pdla = FALSE; */
/*         phplhp = FALSE; */
/*       } */
/*       else if(strcmp(line,"# pd/la") == 0){ */
/*         isbond = FALSE; */
/*         isang = FALSE; */
/*         isdih = FALSE; */
/*         pald = FALSE; */
/*         pdla = TRUE; */
/*         phplhp = FALSE; */
/*       } */
/*       else if(strcmp(line,"# php/lhp") == 0){ */
/*         isbond = FALSE; */
/*         isang = FALSE; */
/*         isdih = FALSE; */
/*         pald = FALSE; */
/*         pdla = FALSE; */
/*         phplhp = TRUE; */
/*       } */

      else{
        isbond = FALSE;
        isang = FALSE;
        isdih = FALSE;
        pald = FALSE;
        pdla = FALSE;
        phplhp = FALSE;
      }
    }
    else{
      i++;
      if(bVerbose)
        progress(stderr,"tCNC__log_> Reading bounds ",i,nbounds);
      
      if(isbond || isdih){
#ifdef GMX_DOUBLE
        if(sscanf(line,"%d %d %lf %lf %lf",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#else
        if(sscanf(line,"%d %d %f %f %f",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#endif
        input_error(i,line);
/*         if(b->at1[i-1] > b->at2[i-1]){ */
/*           dum = b->at1[i-1]; */
/*           b->at1[i-1] = b->at2[i-1]; */
/*           b->at2[i-1] = dum; */
/*         } */
        b->at1[i-1]-=1;
        b->at2[i-1]-=1;
        b->at3[i-1]-=1;
        b->at4[i-1]-=1;
/*         b->av[i-1] = sqr(b->av[i-1]); */
/*         b->lb[i-1] = sqr(b->lb[i-1]); */
/*         b->ub[i-1] = sqr(b->ub[i-1]); */
        b->isbond[i-1] = isbond;
        b->isang[i-1] = isang;
        b->isdih[i-1] = isdih;
        b->n++;
      }
      else if(isang){
#ifdef GMX_DOUBLE
        if(sscanf(line,"%d %d %lf %lf %lf %d %lf %lf",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1],&b->at3[i-1],
                  &b->ang[i-1],&b->sig[i-1])!= 8)
#else
        if(sscanf(line,"%d %d %f %f %f %d %f %f",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1],&b->at3[i-1],
                  &b->ang[i-1],&b->sig[i-1])!= 8)
#endif
        input_error(i,line);
/*         if(b->at1[i-1] > b->at2[i-1]){ */
/*           dum = b->at1[i-1]; */
/*           b->at1[i-1] = b->at2[i-1]; */
/*           b->at2[i-1] = dum; */
/*         } */
        b->at1[i-1]-=1;
        b->at2[i-1]-=1;
        b->at3[i-1]-=1;
        b->at4[i-1]-=1;
/*         b->av[i-1] = sqr(b->av[i-1]); */
/*         b->lb[i-1] = sqr(b->lb[i-1]); */
/*         b->ub[i-1] = sqr(b->ub[i-1]); */
        b->isbond[i-1] = FALSE;
        b->isang[i-1] = TRUE;
        b->isdih[i-1] = FALSE;
        b->n++;
      }
/*       else if(isdih){ */
/*         if(sscanf(line,"%d %d %d %d %f %f",&b->at1[i-1],&b->at2[i-1], */
/*                   &b->at3[i-1],&b->at4[i-1],&b->ang[i-1],&b->sig[i-1]) != 6) */
/*         input_error(i); */
/*         b->at1[i-1]-=1; */
/*         b->at2[i-1]-=1; */
/*         b->at3[i-1]-=1; */
/*         b->at4[i-1]-=1; */
/*         b->isbond[i-1] = FALSE; */
/*         b->isang[i-1] = FALSE; */
/*         b->isdih[i-1] = TRUE; */
/*         b->n++; */
        
/*       } */
      
      else {
#ifdef GMX_DOUBLE
        if(sscanf(line,"%d %d %lf %lf %lf",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#else
        if(sscanf(line,"%d %d %f %f %f",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#endif
        input_error(i,line);
/*         if(b->at1[i-1] > b->at2[i-1]){ */
/*           dum = b->at1[i-1]; */
/*           b->at1[i-1] = b->at2[i-1]; */
/*           b->at2[i-1] = dum; */
/*         } */
        b->at1[i-1]-=1;
        b->at2[i-1]-=1;
        b->at3[i-1]-=1;
        b->at4[i-1]-=1;
/*         b->av[i-1] = sqr(b->av[i-1]); */
/*         b->lb[i-1] = sqr(b->lb[i-1]); */
/*         b->ub[i-1] = sqr(b->ub[i-1]); */
        b->isbond[i-1] = FALSE;
        b->isang[i-1] = FALSE;
        b->isdih[i-1] = FALSE;
        if(pald) b->pald[i-1] = TRUE;
        if(pdla) b->pdla[i-1] = TRUE;
        if(phplhp) b->phplhp[i-1] = TRUE;
        b->n++;
      }
      if(b->at1[i-1] == b->at2[i-1]) {
        input_error(i,line);
      }
    }
  }
  if (bVerbose) {
    fprintf(stderr,"\n");
  }
  
  return b;
}

/*================================================================*/
void dump_bounds(FILE *fp,t_bounds *b)
{
  int i;
  fprintf(fp,"=========================================\n");
  fprintf(fp,"#\t constraints\n");
  fprintf(fp,"# at1  at2  at3  at4  av  lb  ub  ang  dih   sig\n");
  fprintf(fp,"=========================================\n");
  for(i=0;i<b->n;i++){
    fprintf(fp,"%8d %8d %8d %8d %g %g %g %g %g %g\n",
            b->at1[i],b->at2[i],b->at3[i],b->at4[i],
            b->av[i],b->lb[i],b->ub[i],b->ang[i],b->dih[i],
            b->sig[i]);
  }
}
/*================================================================*/
void copy_bound(t_bounds *b1, int i, t_bounds *b2, int k)
{
  b1->at1[i] = b2->at1[k];
  b1->at2[i] = b2->at2[k];
  b1->at3[i] = b2->at3[k];
  b1->at4[i] = b2->at4[k];
  b1->av[i] = b2->av[k];
  b1->lb[i] = b2->lb[k];
  b1->ub[i] = b2->ub[k];
  b1->ang[i] = b2->ang[k];
  b1->sig[i] = b2->sig[k];
  b1->isbond[i] = b2->isbond[k];
  b1->isang[i] = b2->isang[k];
  b1->isdih[i] = b2->isdih[k];
  b1->bCheck[i] = b2->bCheck[k];
  b1->pdla[i] = b2->pdla[k];
  b1->pald[i] = b2->pald[k];
  b1->phplhp[i] = b2->phplhp[k];
  b1->ac_ac[i] = b2->ac_ac[k];
  b1->don_ac[i] = b2->don_ac[k];
  b1->don_don[i] = b2->don_don[k];
  
}

/*================================================================*/
t_bounds *sort_bounds(FILE *fp,t_bounds *b, int natoms)
{
  t_bounds *new = bounds_init(b->n);
  
  int count = 0;
  int i,k;
  i=0;
  while(i<natoms){
    if(fp)
      progress(fp,"Sorting bounds ",i,natoms);
    for(k=0;k<b->n;k++){
      if(b->at1[k]==i){
        copy_bound(new,count,b,k);
        count++;
      }
    }
    i++;
  }
  new->n = b->n;
  printf("\n");
  return new;
}

/*================================================================*/

void check_bounds(FILE *log, t_bounds *b,t_atomlist *al, t_idxgroups *pln,
                  t_idxgroups *imp)
{
  int i,k;
  real di,ang;
  char err[STRLEN];
  bool bOk = TRUE;
  bool aOk = TRUE;
  char wstr[STRLEN];
  char bflag[STRLEN];
  real viol = 0;
  
/*   fprintf(stderr,"\nChecking bounds......\n"); */
  
  for(i=0;i<b->n;i++){
    if(b->at1[i] == b->at2[i])
    {
      sprintf(err,"error while reading bounds: %d %d\n",b->at1[i],b->at2[i]);
      CNCerr(log,err);
    }
    else if(b->at1[i] ==-1 ||
            b->at2[i] ==-1){
      sprintf(err,"error while reading bounds: %d %d\n",b->at1[i],b->at2[i]);
      CNCerr(log,err);

    }
    di = DIST(al,b->at1[i],b->at2[i]);
    if( di < b->lb[i] || di > b->ub[i]){
      if (di > b->ub[i]) viol = di-b->ub[i];
      else if(di < b->lb[i]) viol = b->lb[i]-di;
      
      if(b->isbond[i]) strcpy(bflag,"bond");
      else if(b->isang[i]) strcpy(bflag,"angle");
      else if(b->isdih[i]) strcpy(bflag,"1-4 pair");
      else strcpy(bflag,"other");
      bOk = FALSE;
      sprintf(wstr,"Dist. out of bounds -> d: %6.3f av: %6.3f lb: %6.3f ub: %6.3f (%d-%4s-%d(%3s)<->%d-%4s-%d(%3s)) v=%g !! %s\n",
              di,b->av[i],b->lb[i],b->ub[i],al->id[b->at1[i]],al->name[b->at1[i]],al->resid[b->at1[i]],al->resname[b->at1[i]],
              al->id[b->at2[i]],al->name[b->at2[i]],al->resid[b->at2[i]],al->resname[b->at2[i]],viol,bflag);
      CNCwarn(log,wstr);
      
    }
    if(b->isang[i]){
      ang = RAD2DEG*angle_ij_ik(al->x[b->at3[i]],al->x[b->at2[i]],al->x[b->at1[i]]);
      if(ang > b->ang[i]+b->sig[i] || ang < b->ang[i]-b->sig[i]) {
        sprintf(wstr,"Angle out of bounds -> a: %6.3f  av: %6.3f sig: %6.3f (%d%s-%d%s-%d%s)\n",
                ang, b->ang[i],b->sig[i],al->id[b->at1[i]],al->name[b->at1[i]],
                al->id[b->at3[i]],al->name[b->at3[i]],
                al->id[b->at2[i]],al->name[b->at2[i]]);
        
        CNCwarn(log,wstr);
        aOk = FALSE;
      }
    }
    

    if(b->lb[i] > b->ub[i]){
      sprintf(wstr," lb > ub  %d--%d\n",b->at1[i]+1,b->at2[i]+1);
      CNCerr(log,wstr);
    }
    if ( (b->ub[i] - b->lb[i]) < 0.03) {
      sprintf(wstr,"tCNC__wrn_> Overtight bound   %d--%d\n",b->at1[i]+1,b->at2[i]+1);
      CNClog(log,wstr);
    }
    
  }

  if(bOk) {
    sprintf(wstr,"tCNC__log_> Distance Constraints........: Ok\n");
  }
  
  else{
    sprintf(wstr,"tCNC__wrn_> Distance Constraints........: Not Ok !!\n");
  }
  CNClog(log,wstr);
  
  if(aOk)
  {
    sprintf(wstr,"tCNC__log_> Angle Constraints...........: Ok\n");
  }
  else{
    sprintf(wstr,"tCNC__wrn_> Angle Constraints...........: Not Ok !!\n");
  }
  CNClog(log,wstr);
  
  int at1,at2;
  bool pOk = TRUE;
  
/*   fprintf(stderr,"Checking planar groups.......\n"); */
  for(i=0;i<pln->n;i++){
    real plan;
    bool plnOk = is_planar_group(al,pln->natoms[i],pln->atoms[i],pln->val[i],&plan);
    
    if(!plnOk){
      pOk = FALSE;
      char wstr[STRLEN];
      sprintf(wstr,"tCNC__wrn_> Planarity of group %d not ok\n",i);
      CNClog(log,wstr);
      sprintf(wstr,"tCNC__wrn_> tol = %g plan = %g\n",pln->val[i],plan);
/*       print_group(al,pln->atoms[i],pln->natoms[i]); */
      CNClog(log,wstr);
    }
    fflush(log);
    
  }

  if(pOk)
  {
    sprintf(wstr,"tCNC__log_> Planar Groups...............: Ok\n");
    CNClog(log,wstr);
  }
  
  else
  {
    sprintf(wstr,"tCNC__wrn_> Planar Groups...............: Not Ok !!\n");
    CNClog(log,wstr);
  }
  

/*   fprintf(stderr,"Checking chiral groups.......\n"); */
  int imps[imp->n];
  bool cOk = TRUE;
  
  check_impr(al,imp,imps);
  for(i=0;i<imp->n;i++){
    if(imps[i]){
      cOk = FALSE;
      char wstr[STRLEN];
      sprintf(wstr,"tCNC__wrn_> Chirality of group %d not ok\n",i);
      CNClog(log,wstr);
/*       print_group(al,imp->atoms[i],imp->natoms[i]); */
    }
  }
  if(cOk) {
    sprintf(wstr,"tCNC__log_> Chiral Groups...............: Ok\n");
    CNClog(log,wstr);
  }
  
  else{
    sprintf(wstr,"tCNC__wrn_> Chiral Groups...............: Not Ok !!\n");
    CNClog(log,wstr);
  }
  
  sprintf(wstr,"tCNC__log_> Checking Bounds....... done\n\n");
  
}

/*=========================================================*/

void damp_bounds(t_bounds *b, real damp)
{
  int i;
  real range,tol;
  
  for(i=0;i<b->n;i++){
    
    if(!b->isbond[i] && !b->isang[i] && b->isdih[i])
    {
      range = b->ub[i]-b->lb[i];
      if(range > 4.)
      {
        tol = range/(damp*2);
        if(tol < 2.) tol = 2.;
        b->ub[i] = b->av[i]+tol;
        b->lb[i] = b->av[i]-tol;
      }
    }
  }
}
