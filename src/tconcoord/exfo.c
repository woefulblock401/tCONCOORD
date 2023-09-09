#include <tconcoord.h>
/*=============================================================*/
t_excl *excl_init(void)

/* initialises exclusion structure */ 
{
  t_excl *ex = NULL;
  snew(ex,1);
  ex->n = 0;
  ex->id1 = NULL;
  ex->id2 = NULL;
  return ex;
}
/*=============================================================*/
t_force *force_init(void)
  /* initialises forced constraints structure */

{
  t_force *fo = NULL;
  snew(fo,1);
  fo->n = 0;
  fo->id1 = NULL;
  fo->id2 = NULL;
  fo->lb = NULL;
  fo->ub = NULL;
  return fo;
}
/*=============================================================*/


void read_exfo(char *filename, t_excl *ex, t_force *fo)
{
  
#define MAXPTR 254

  int ninp;
  int i,k;
  int nel;
  char *ptr[MAXPTR];
  char error[STRLEN];
  t_inpfile *inf = read_inpfile(filename,&ninp);
  char dum[STRLEN][3];
  
  for(i=0;i<ninp;i++){
    if(strcmp(inf[i].name,"excl") == 0){
      nel = str_nelem(inf[i].value,MAXPTR,ptr);
      if(nel != 2){
        sprintf(error,"Invalid input: %s = %s\n",inf[i].name,inf[i].value);
        fatal_error(error);
      }
      ex->n+=1;
      srenew(ex->id1,ex->n);
      srenew(ex->id2,ex->n);
      sscanf(inf[i].value,"%d %d",&ex->id1[ex->n-1],&ex->id2[ex->n-1]);
    }
    else if(strcmp(inf[i].name,"forc") == 0){
      nel = str_nelem(inf[i].value,MAXPTR,ptr);
      if(nel != 4){
        sprintf(error,"Invalid input: %s = %s\n",inf[i].name,inf[i].value);
        fatal_error(error);
      }
      fo->n+=1;
      srenew(fo->id1,fo->n);
      srenew(fo->id2,fo->n);
      srenew(fo->lb,fo->n);
      srenew(fo->ub,fo->n);
#ifdef GMX_DOUBLE
      sscanf(inf[i].value,"%d %d %lf %lf",&fo->id1[fo->n-1],&fo->id2[fo->n-1],
             &fo->lb[fo->n-1],&fo->ub[fo->n-1]);
#else
      sscanf(inf[i].value,"%d %d %f %f",&fo->id1[fo->n-1],&fo->id2[fo->n-1],
             &fo->lb[fo->n-1],&fo->ub[fo->n-1]);
#endif
    }
  }
}
/*==========================================================*/
bool is_excluded(t_excl *ex, t_atomlist *al,int id1, int id2)

/* returns TRUE if atom id1 and 
   atom id2 are excluded from constraint 
   definition*/

{
  int i;
  for(i=0;i<ex->n;i++){
    if((ex->id1[i] == al->resid[id1] &&
       ex->id2[i] == al->resid[id2]) ||
       (ex->id1[i] == al->resid[id2] &&
       ex->id2[i] == al->resid[id1])) return TRUE;
  }
  return FALSE;
}
/*=============================================================*/
bool is_forced(t_force *fo, t_atomlist *al,int id1, int id2)

/* returns TRUE if atom id1 and 
   atom id2 are forced constraints */
{
  int i;
  for(i=0;i<fo->n;i++){
    if((fo->id1[i] == al->resid[id1] &&
       fo->id2[i] == al->resid[id2]) ||
       (fo->id1[i] == al->resid[id2] &&
       fo->id2[i] == al->resid[id1])) return TRUE;
  }
  return FALSE;
}

/*=============================================================*/

