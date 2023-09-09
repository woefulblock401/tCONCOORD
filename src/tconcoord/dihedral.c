#include <tconcoord.h>

/* THIS FILE IS OUT OF DATE AND SHOULD NOT BE USED !!!! */

/*=============================================================*/
int check_dihedral(rvec at1, rvec at2, rvec at3, rvec at4,
                   real ang, real sig)
{
/*   real phi; */
  real dihed = RAD2DEG*dihedral(at1,at2,at3,at4);
/*   printf("dihed = %g ang = %g\n",dihed,ang);      */
  if(ang < 10.0 && dihed > 180){
    dihed = dihed - 360;
  }
  
  if(dihed > ang + sig ||
     dihed < ang - sig)
  {

    if(dihed > ang + sig) {
      rotate_atom(at2,at3,at4,0.5*DEG2RAD);
      rotate_atom(at3,at2,at1,0.5 *DEG2RAD);
    }
    else{
      rotate_atom(at2,at3,at4,0.5*DEG2RAD);
      rotate_atom(at3,at2,at1,0.5*DEG2RAD);
    }
     real dihed = RAD2DEG*dihedral(at1,at2,at3,at4); 
/*     printf("\tdihed = %g\n",dihed);       */

    return 1;
  }
  else 
    return 0;
}
/*=============================================================*/

int check_omega(t_atomlist *al)
{
  int i,k,j,l;
  int at1, at2, at3, at4;

  int ncheck = 0;
  
  for(i=0;i<al->natoms;i++){
    at1 = at2 = at3 = at4 = -1;
    
    if(strcmp(al->name[i]," CA ")==0){
      at1 = i;
      for(k=0;k<al->nbonds[at1];k++){
        if(strcmp(al->name[al->bonds[at1][k]]," C  ") == 0){
          at2 = al->bonds[at1][k];
          for(j=0;j<al->nbonds[at2];j++){
            if(strcmp(al->name[al->bonds[at2][j]]," N  ") == 0){
              at3 = al->bonds[at2][j];
              for(l=0;l<al->nbonds[at3];l++){
                if(strcmp(al->name[al->bonds[at3][l]]," CA ") ==0){
                  at4 = al->bonds[at3][l];
                  ncheck+=check_dihedral(al->x[at1],al->x[at2],al->x[at3],
                                         al->x[at4],180.0,12.);
                }
                else if(strcmp(al->name[al->bonds[at3][l]]," H  ") ==0){
                  at4 = al->bonds[at3][l];
                  ncheck+=check_dihedral(al->x[at1],al->x[at2],al->x[at3],
                                 al->x[at4],0.,12.);

                }
              }
            }
          }
        }
      }
    }
  }

  for(i=0;i<al->natoms;i++){
    at1 = at2 = at3 = at4 = -1;
    
    if(strcmp(al->name[i]," H  ")==0){
      at1 = i;
      for(k=0;k<al->nbonds[at1];k++){
        if(strcmp(al->name[al->bonds[at1][k]]," N  ") == 0){
          at2 = al->bonds[at1][k];
          for(j=0;j<al->nbonds[at2];j++){
            if(strcmp(al->name[al->bonds[at2][j]]," C  ") == 0){
              at3 = al->bonds[at2][j];
              for(l=0;l<al->nbonds[at3];l++){
                if(strcmp(al->name[al->bonds[at3][l]]," O  ") ==0){
                  at4 = al->bonds[at3][l];
                  ncheck+=check_dihedral(al->x[at1],al->x[at2],al->x[at3],
                                 al->x[at4],180.0,12.);

                }
              }
            }
          }
        }
      }
    }
  }
  for(i=0;i<al->natoms;i++){
    at1 = at2 = at3 = at4 = -1;
    
    if(strcmp(al->name[i]," O  ")==0){
      at1 = i;
      for(k=0;k<al->nbonds[at1];k++){
        if(strcmp(al->name[al->bonds[at1][k]]," C  ") == 0){
          at2 = al->bonds[at1][k];
          for(j=0;j<al->nbonds[at2];j++){
            if(strcmp(al->name[al->bonds[at2][j]]," N  ") == 0){
              at3 = al->bonds[at2][j];
              for(l=0;l<al->nbonds[at3];l++){
                if(strcmp(al->name[al->bonds[at3][l]]," CA  ") ==0){
                  at4 = al->bonds[at3][l];
                  ncheck+=check_dihedral(al->x[at1],al->x[at2],al->x[at3],
                                         al->x[at4],0.,12.);
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  return ncheck;
}

/*=============================================================*/


