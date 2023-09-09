#include <tconcoord.h>

/*============================================================*/

t_contab *contab_init(void)
/* initialise connection table */
{
  t_contab *ct = NULL;
  snew(ct,1);
  return ct;
}
/*============================================================*/
t_contab *contab_realloc(t_contab *ct, int n)
/* reallocate connection table */
{
  int i; 
  ct->n = n; 
  snew(ct->ncon,n);
  snew(ct->con,n);
  snew(ct->type,n);
  for(i=0;i<ct->n;i++){
    ct->ncon[i] = 0;
    snew(ct->con[i],MAX_CON);
    snew(ct->type[i],MAX_CON);
  }
  return ct;
}
/*============================================================*/
void free_contab(t_contab *ct)
/* free connection table */
{
  int i;
  for(i=0;i<ct->n;i++){
    sfree(ct->con[i]);
    sfree(ct->type[i]);
  }
  sfree(ct->ncon);
  sfree(ct->con);
  sfree(ct->type);
  sfree(ct);
}

/*============================================================*/
bool ct_connected(t_contab *ct, int i, int k, int flag1, int flag2)
{
  int j;
  for(j=0;j<ct->ncon[i];j++){
    if(ct->con[i][j] == k &&
       ct->type[i][j] == flag1) return TRUE;
  }
  for(j=0;j<ct->ncon[k];j++){
    if(ct->con[k][j] == i &&
       ct->type[k][j] == flag2) return TRUE;
  }
  return FALSE;
}

/*============================================================*/
void add2contab(t_contab *ct, int i, int k, int flag1, int flag2,bool check)
{
  /* first check if connection is already listed */
  int j;
  bool doit=TRUE;
  if(check)
    doit = !connected(ct,i,k) && !connected(ct,k,i);  
  else{
    doit = !ct_connected(ct,i,k,flag1,flag2);
  }
  if(doit){ 

    ct->ncon[i]+=1; 
    ct->ncon[k]+=1; 
    if(ct->ncon[i] > MAX_CON){
      srenew(ct->con[i],ct->ncon[i]); 
      srenew(ct->type[i],ct->ncon[i]);
    }
    if(ct->ncon[k] > MAX_CON){
      srenew(ct->con[k],ct->ncon[k]);
      srenew(ct->type[k],ct->ncon[k]);
    }
    ct->con[i][ct->ncon[i]-1] = k;
    ct->type[i][ct->ncon[i]-1] = flag1;
    ct->con[k][ct->ncon[k]-1] = i;
    ct->type[k][ct->ncon[k]-1] = flag2;
  }  
}
                
/*============================================================*/

void add_bond(t_atomlist *al,int i, int k, t_contab *ct,int flag1, int flag2)
{
  /* here we add an atom to the appropriate
     array of another one and also add
     this info to the connection table
  */

  int nb1,nb2;
  switch(flag1){
    case BOND:
/*       printf("here now\n"); */
      
      al->nbonds[i]++;
      al->nbonds[k]++;
      nb1 = al->nbonds[i];
      nb2 = al->nbonds[k];
      if(nb1 > MAX_BOND){
        srenew(al->bonds[i],nb1);
      }
      if(nb2 > MAX_BOND) {
        srenew(al->bonds[k],nb2);
      }
      al->bonds[i][nb1-1] = k;
      al->bonds[k][nb2-1] = i;
      break;
    case B13:
      al->nb13[i]++;
      al->nb13[k]++;
      nb1 = al->nb13[i];
      nb2 = al->nb13[k];
      if(nb1 > MAX_B13){
        srenew(al->b13[i],nb1);
      }
      if(nb2 > MAX_B13) {
        srenew(al->b13[k],nb2);
      }
      al->b13[i][nb1-1] = k;
      al->b13[k][nb2-1] = i;
      break;
    case B14:
      al->nb14[i]++;
      al->nb14[k]++;
      nb1 = al->nb14[i];
      nb2 = al->nb14[k];
      if(nb1 > MAX_B14){
        srenew(al->b14[i],nb1);
      }
      if(nb2 > MAX_B14) {
        srenew(al->b14[k],nb2);
      }
      al->b14[i][nb1-1] = k;
      al->b14[k][nb2-1] = i;
      break;
    case NL:
      al->nnb[i]++;
      al->nnb[k]++;
      nb1 = al->nnb[i];
      nb2 = al->nnb[k];
      if(nb1 > MAX_NL*al->memcount[i]){
        al->memcount[i]+=1;
        printf("reallocate nlist %d\n",al->memcount[i]);
        
        srenew(al->nb[i],MAX_NL*al->memcount[i]);
      }
      if(nb2 > MAX_NL*al->memcount[k]) {
        al->memcount[k]+=1;
        printf("reallocate nlist %d\n",al->memcount[k]);
        srenew(al->nb[k],MAX_NL*al->memcount[k]);
      }
      al->nb[i][nb1-1] = k;
      al->nb[k][nb2-1] = i;
      break;
  }
  
/*   printf("adding\n"); */
  if(ct != NULL)
    add2contab(ct,i,k,flag1,flag2,1);  
/* printf("adding2\n"); */
}

/*============================================================*/
bool connected(t_contab *ct, int i, int k)
{
  int j;
  for(j=0;j<ct->ncon[i];j++){
    if(k == ct->con[i][j]) return TRUE;
  }
  return FALSE;
}
/*============================================================*/

void print_contab(t_contab *rct, char *filename)
{
  int i,k;
  shstr dum;
  FILE *fp = ffopen(filename,"w");

  for(i=0;i<rct->n;i++){
    for(k=0;k<rct->ncon[i];k++){
      switch(rct->type[i][k]){
        case eCOV:
          strcpy(dum,"COV");
          break;
        case eSUL:
          strcpy(dum,"SUL");
          break;
        case eBBHD:
          strcpy(dum,"BBHD");
          break;
        case eBBHA:
          strcpy(dum,"BBHA");
          break;
        case eSBH_BD:
          strcpy(dum,"SBH_BD");
          break;
        case eSBH_BA:
          strcpy(dum,"SBH_BA");
          break;
        case eSBH_SD:
          strcpy(dum,"SBH_SD");
          break;
        case eSBH_SA:
          strcpy(dum,"SBH_SA");
          break;
        case eSSHD:
          strcpy(dum,"SSHD");
          break;
        case eSSHA:
          strcpy(dum,"SSHA");
          break;
        case eNPR:
          strcpy(dum,"NPR");
          break;
        case ePHO:
          strcpy(dum,"PHO");
          break;
        case ePACK:
          strcpy(dum,"PACK");
          break;
        case eMET:
          strcpy(dum,"MET");
          break;
        case eEX:
          strcpy(dum,"EX");
          break;
        case eFO:
          strcpy(dum,"FO");
          break;
      }
      fprintf(fp,"%5d %5d %6s\n",i+1,rct->con[i][k]+1,dum);
    }
  }
}

      
  

