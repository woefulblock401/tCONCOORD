#include <tconcoord.h>

/*=============================================================*/

t_bondlist *bl_init(void)

/* initialise bondlist */

{
  t_bondlist *bl = NULL;
  snew(bl,1);
  bl->n = 0;
  bl->at1 = NULL;
  bl->at2 = NULL;
  bl->blen = NULL;
  bl->type = NULL;
  bl->restricted = NULL;
  return bl;
}
/*=============================================================*/
t_bondlist *bl_realloc(t_bondlist *bl, int n)

/* reallocate bondlist */

{
  srenew(bl->at1,n);
  srenew(bl->at2,n);
  srenew(bl->restricted,n);
  srenew(bl->blen,n);
  srenew(bl->type,n);
  return bl;
}
/*=============================================================*/

void restrict_bond(t_bondlist *bl, int at1, int at2)

/* self explaining */

{
  int i;
  for(i=0;i<bl->n;i++){
    if((bl->at1[i] == at1 && bl->at2[i] == at2) ||
       (bl->at1[i] == at2 && bl->at2[i] == at1)){
      bl->restricted[i] = TRUE;
      break;
    }
  }
}
/*=============================================================*/
bool is_bond(t_atomlist *al, int at1, int at2)

/* check whether at1 is bound to at2 */

{
  int i;
  for(i=0;i<al->nbonds[at1];i++){
    if(al->bonds[at1][i]==at2) return TRUE;
  }
  return FALSE;
}

/*=============================================================*/


void rotatable_bonds(t_atomlist *al,t_bondlist *bl)

/* check which bonds are rotatable */

{
  int i,k,j,l;
  int at,at2,at3;
  real d;
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      if(i > al->bonds[i][k]){
        at = al->bonds[i][k];
        bl->n++;
        bl = bl_realloc(bl,bl->n);
        bl->at1[bl->n-1] = i;
        bl->at2[bl->n-1] = at;
        bl->restricted[bl->n-1] = FALSE;
        d = DIST(al,i,at);
        bl->blen[bl->n-1] = d;
  
        if(al->nbonds[i] == 1 ||
           al->nbonds[at] == 1){
          bl->restricted[bl->n-1] = TRUE;
          bl->type[bl->n-1] = SINGLE;
        }
        
        if((strcmp(al->hyb[i],"sp3") == 0 ||
            strcmp(al->hyb[at],"sp3") == 0) &&
           (al->nbonds[i] > 1 &&
            al->nbonds[at] > 1)){
          bl->restricted[bl->n-1] = FALSE;
          bl->type[bl->n-1] = SINGLE;
          
        }
        else if(strcmp(al->hyb[i],"sp2") == 0 &&
                strcmp(al->hyb[at],"sp2") == 0){
          if(strcmp(al->symbol[i],"C") == 0 &&
              strcmp(al->symbol[at],"C") == 0)
          {
            if(d < 1.46 && d >= 1.30){
              bl->restricted[bl->n-1] = TRUE;
              bl->type[bl->n-1] = DOUBLE;
            }
            else if(d < 1.33){
              bl->restricted[bl->n-1] = TRUE;
              bl->type[bl->n-1] = TRIPLE;
            }
            else if(d >= 1.46){
              bl->restricted[bl->n-1] = FALSE;
              bl->type[bl->n-1] = SINGLE;
            }
          }
          else if((strcmp(al->symbol[i],"C") == 0 &&
                   strcmp(al->symbol[at],"N") == 0) ||
                  (strcmp(al->symbol[i],"N") == 0 &&
                   strcmp(al->symbol[at],"C") == 0))
          {

            if(d > 1.42 && d < 2.){
              bl->restricted[bl->n-1] = FALSE;
              bl->type[bl->n-1] = CNSINGLE;
            }
            else if(d <=1.42 && d>1.25){
              bl->restricted[bl->n-1] = TRUE;
              bl->type[bl->n-1] = CNOMEGA;
            }
            else if(d<=1.25){
              bl->restricted[bl->n-1] = TRUE;
              bl->type[bl->n-1] = CNDOUBLE;
            }
            else
            {
              bl->restricted[bl->n-1] = TRUE;
              bl->type[bl->n-1] = CNOMEGA;
              printf("Defining bond %d-%d\n",i,at);
            }
              
          }
             
        }
      }
    }
  }
}

/*=============================================================*/
void print_rotatable_bonds(FILE *log,t_atomlist *al, t_bondlist *bl)

/* print a list of all rotatable bonds */

{
  int i;
  for(i=0;i<bl->n;i++){
    if(!bl->restricted[i]){
      fprintf(log,"%s(%d%s) -- %s(%d%s)\n",al->name[bl->at1[i]],al->resid[bl->at1[i]],
              al->resname[bl->at1[i]],al->name[bl->at2[i]],al->resid[bl->at2[i]],
              al->resname[bl->at2[i]]);
    }
  }
}

