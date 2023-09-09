#include <tconcoord.h>

/*=========================================================*/
int random_int(gmx_rng_t rng, int max)
{
  int iran;
  real ran = gmx_rng_uniform_real(rng);
  if(ran == 1) ran = 0.99999;
  ran = ran*max;
  iran = (int) ran;
  return iran;
}
/*=========================================================*/

int random_int_in_range(gmx_rng_t rng, int min, int max)
{
  int iran;
  real ran = gmx_rng_uniform_real(rng);
  ran = ran*(max-min)+min;
  iran = (int)ran;
  return iran;
}

/*=========================================================*/

void rand_array_int(gmx_rng_t rng, int n, int *array)
{
  int i;
  int idx1, idx2;
  for(i=0;i<n;i++){
    array[i]=i;
  }
  
  i = 0;
  while(i<n)
  {
    idx1 = random_int(rng,n);
    idx2 = random_int(rng,n);
    switch_int(&array[idx1],&array[idx2]);
    i++;
  }
}
/*=========================================================*/

void simple_array_int(int *array, int n)
{
  int i;
  for(i=0;i<n;i++){
    array[i] = i;
  }
}

/*=========================================================*/
void random_al_simple(t_atomlist *al, gmx_rng_t rng)
{
  int i,k;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      al->x[i][k] = (gmx_rng_uniform_real(rng)-0.5)   *100;
    }
  }
}

/*=========================================================*/

void random_al(t_atomlist *al,gmx_rng_t rng, bool bTarget, rvec tco)
{
  int i,k;
  
/* if we have posres, we use the center of
   the position restraint atoms 
*/
  int size = 0;
  rvec co;
  clear_rvec(co);
  atom_id *index = NULL;
  
  for(i=0;i<al->natoms;i++){
    if(al->isposres[i]){
      size++;
      srenew(index,size);
      index[size-1] = i;
    }
  }
  /**/
  if(size && !bTarget){
    rvec dum[size];
    com_frag(al,size,index,dum,co);
  }
  else if(bTarget){
/*     printf("Using Target center of mass\n"); */
    
    copy_rvec(tco,co);
  }
  
/*   printf("co = %g %g %g\n",co[0],co[1],co[2]); */
  
  for(i=0;i<al->natoms;i++){
    if(!al->isposres[i]){
      for(k=0;k<DIM;k++){
        al->x[i][k] = (gmx_rng_uniform_real(rng)-0.5)   *100   + co[k];  
/*         al->x[i][k] = (gmx_rng_uniform_real(rng)-0.5); */
        

      }
    }
  }
}

/*=========================================================*/

void perturb_al(t_atomlist *al,gmx_rng_t rng)
{
  int i,k;
  for(i=0;i<al->natoms;i++){
    if(!al->isposres[i]){
      for(k=0;k<DIM;k++){
        al->x[i][k] += (gmx_rng_uniform_real(rng)-0.5)*4.;
      }
    }
  }
}

/*=========================================================*/
void random_al_idx(t_atomlist *al,int size, atom_id *index, gmx_rng_t rng)
{
  int i,k;
  rvec cm;
  rvec new[size];
  com_frag(al,size,index,new,cm);

  for(i=0;i<size;i++){
    for(k=0;k<DIM;k++){
      al->x[index[i]][k] = cm[k]+(gmx_rng_uniform_real(rng)-0.5)*100;

    }
  }
}
/*=========================================================*/

void random_al_idx2(t_atomlist *al,int size, atom_id *index, 
                    gmx_rng_t rng, rvec cent, rvec bsize)
{
  int i,k;
/*   rvec co; */
/*   calc_com(al,co); */
  
  
  for(i=0;i<size;i++){
    for(k=0;k<DIM;k++){
      al->x[index[i]][k] = cent[k]+(gmx_rng_uniform_real(rng)-0.5)*bsize[k]; 
/*       al->x[index[i]][k] = co[k]+gmx_rng_uniform_real(rng)*100.-50.; */
    }
/*      al->x[index[i]][ZZ] = -50.;  */
    
  }
}
/*=========================================================*/

void randomise_group(t_atomlist *al, gmx_rng_t rng, int grpnr)
{
  int i,k,j;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->ngrps[i];k++){
      if(al->grpnr[i][k] == grpnr){
        for(j=0;j<DIM;j++){
        al->x[i][j] = gmx_rng_uniform_real(rng)*100.;
        }
      }
    }
  }
} 

