#include <tconcoord.h>
/*================================================*/
t_dihed *dihed_init(void)
{
  t_dihed *dh = NULL;
  snew(dh,1);
  dh->n = 0;
  dh->center = NULL;
  dh->at1=dh->at2=dh->at3=dh->at4 = NULL;
  dh->tp1=dh->tp2=dh->tp3=dh->tp4 = NULL;
  dh->phi = NULL;
  dh->flag = NULL;
  dh->restr = NULL;
  return dh;
}
/*================================================*/

t_dihed *dihed_realloc(t_dihed *dh, int n)
{
  dh->n = n;
  srenew(dh->center,n);
  srenew(dh->at1,n);
  srenew(dh->at2,n);
  srenew(dh->at3,n);
  srenew(dh->at4,n);
  srenew(dh->tp1,n);
  srenew(dh->tp2,n);
  srenew(dh->tp3,n);
  srenew(dh->tp4,n);
  srenew(dh->phi,n);
  srenew(dh->flag,n);
  srenew(dh->restr,n);
  return dh;
}
/*================================================*/

void fill_dihed(t_dihed *dh, t_atomlist *al)
{
  int i;
  for(i=0;i<dh->n;i++){
    dh->phi[i] = DIHED(al,dh->at1[i],
                       dh->at2[i],
                       dh->at3[i],
                       dh->at4[i]);
    if(IS_OMEGA(al,dh->at2[i],dh->at3[i]))
      dh->flag[i] = OMEGA;
    else if(IS_PSI(al,dh->at2[i],dh->at3[i]))
      dh->flag[i] = PSI;
    else if(IS_PHI(al,dh->at2[i],dh->at3[i]))
      dh->flag[i] = PHI;
    else
      dh->flag[i] = FREE;
  }
}
/*================================================*/
void get_dihedrals(t_atomlist *al, t_dihed *dihed)

/* fill dihedral structure from atomlist */

{
  int at1,at2,at3,at4;
  int i,k,j,l,m;
  bool check;
  int count = 0;
  real d;
  
  for(i=0;i<al->natoms;i++){
    at1 = i;
    for(k=0;k<al->nbonds[i];k++){
      at2 = al->bonds[i][k];
      for(l=0;l<al->nbonds[at2];l++){
        at3 = al->bonds[at2][l];
        for(j=0;j<al->nbonds[at3];j++){
          at4 = al->bonds[at3][j];
          check = FALSE;
          for(m=0;m<al->nb14[i];m++){
            if(at4 == al->b14[i][m] && at4 > at1) check = TRUE;
          }
          if(check){
            count++;
            dihed->n = count;
            dihed = dihed_realloc(dihed,count);
            dihed->at1[count-1] = at1;
            dihed->at2[count-1] = at2;
            dihed->at3[count-1] = at3;
            dihed->at4[count-1] = at4;
            d = RAD2DEG*dihedral(al->x[at1],al->x[at2],al->x[at3],al->x[at4]);
            dihed->phi[count-1] = d;
            if(IsProtein(al,at2)){
              if(IS_OMEGA(al,at2,at3)) dihed->flag[count-1] = OMEGA;
              else if(IS_PHI(al,at2,at3)) dihed->flag[count-1] = PHI;
              else if(IS_PSI(al,at2,at3)) dihed->flag[count-1] = PSI;
              else dihed->flag[count-1] = FREE;
            }
            else dihed->flag[count-1] = FREE;
          }
        }
      }
    }
  }
}
/*=============================================================*/

