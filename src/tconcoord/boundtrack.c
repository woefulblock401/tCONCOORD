#include <tconcoord.h>


/*========================================================*/
t_boundtrack *boundtrack_init(void)
{
  t_boundtrack *bt = NULL;
  snew(bt,1);
  bt->natoms = 0;
  bt->id = NULL;
  bt->n = NULL;
  bt->con = NULL;
  return bt;
}
/*========================================================*/

t_boundtrack *boundtrack_realloc(t_boundtrack *bt,int n)
{
  int i;
  bt->natoms = n;

  srenew(bt->id,n);
  srenew(bt->n,n);
  srenew(bt->con,n);
  for(i=0;i<bt->natoms;i++){
    bt->id[i] = 0;
    bt->n[i] = 0;
    bt->con[i] = NULL;
  }

  return bt;
}
/*========================================================*/

void add_bound(t_boundtrack *bt, int i,int k)
{
  bt->n[i]+=1;
  bt->n[k]+=1;
  srenew(bt->con[i],bt->n[i]);
  srenew(bt->con[k],bt->n[k]);
  bt->con[i][bt->n[i]-1] = k;
  bt->con[k][bt->n[k]-1] = i;
}
/*========================================================*/

bool is_bound(t_boundtrack *bt, int i, int k)
{
  int check = FALSE;
  int j;
  for(j=0;j<bt->n[i];j++){
    if(k == bt->con[i][j]){
      check = TRUE;
      break;
    }
  }
  return check;
}
/*========================================================*/

void fill_boundtrack(t_bounds *b, t_boundtrack *bt)
{
  int i;
  for(i=0;i<b->n;i++){
    add_bound(bt,b->at1[i],b->at2[i]);
  }
}
