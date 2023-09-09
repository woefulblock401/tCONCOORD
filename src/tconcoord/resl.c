#include <tconcoord.h>
/*============================================*/
t_resl *resl_init(void)
{
  t_resl *rl = NULL;
  snew(rl,1);
  rl->nres = 0;
  rl->id = NULL;
  rl->resname = NULL;
  rl->j0 = NULL;
  rl->j1 = NULL;
  rl->phi_restr = NULL;
  rl->psi_restr = NULL;
  rl->sc_restr = NULL;
  rl->ngrps = NULL;
  rl->grps = NULL;
  return rl;
}
/*============================================*/

t_resl *resl_realloc(t_resl *rl, int n)
{
  snew(rl->id,n);
  snew(rl->resname,n);
  snew(rl->j0,n);
  snew(rl->j1,n);
  snew(rl->phi_restr,n);
  snew(rl->psi_restr,n);
  snew(rl->sc_restr,n);
  snew(rl->ngrps,n);
  snew(rl->grps,n);
  return rl;
}

/*============================================*/

void fill_resl(t_resl *rl,t_atomlist *al)
{
  int i,k;
  k = 0;
  rl->j0[0] = 0;
  for(i=1;i<al->natoms;i++){
    if(al->resid[i] != al->resid[i-1]){
      rl->id[k] = al->resid[i-1];
      strcpy(rl->resname[k],al->resname[i-1]);
      rl->j1[k] = i-1;
      rl->j0[k+1] = i;
      k++;
    }
  }
  strcpy(rl->resname[rl->nres-1],al->resname[al->natoms-1]);
  rl->j1[rl->nres-1] = al->natoms - 1;
  rl->id[rl->nres-1] = al->resid[al->natoms-1];
}
/*============================================*/

int get_atom_idx(t_resl *rl, int residx, t_atomlist *al, char *name)
{
  int i;
  if(residx >= rl->nres) return -1;
  
  for(i=rl->j0[residx];i<=rl->j1[residx];i++){
    if(strcmp(al->name[i],name) == 0) return i;
  }
  return -1;
}
