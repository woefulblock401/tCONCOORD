#include <tconcoord.h>


static int NVEC = 30;

static rvec vectors[] = {
  /* Vectors for octahedral triangulation */
  {  0.00000000,   0.00000000,   1.00000000}, 
  {  0.00000000,   0.00000000,  -1.00000000},
  {  0.70710678,   0.70710678,   0.00000000},
  {  0.70710678,  -0.70710678,   0.00000000},
  { -0.70710678,   0.70710678,   0.00000000},
  { -0.70710678,  -0.70710678,   0.00000000},
  {  0.86602540,   0.28867513,   0.40824829},
  {  0.86602540,  -0.28867513,   0.40824829},
  {  0.57735027,   0.00000000,   0.81649658},
  {  0.86602540,   0.28867513,  -0.40824829},
  {  0.86602540,  -0.28867513,  -0.40824829},
  {  0.57735027,   0.00000000,  -0.81649658},
  {  0.28867513,  -0.86602540,   0.40824829},
  { -0.28867513,  -0.86602540,   0.40824829},
  {  0.00000000,  -0.57735027,   0.81649658},
  {  0.28867513,  -0.86602540,  -0.40824829},
  { -0.28867513,  -0.86602540,  -0.40824829},
  {  0.00000000,  -0.57735027,  -0.81649658},
  { -0.86602540,   0.28867513,   0.40824829},
  { -0.86602540,  -0.28867513,   0.40824829},
  { -0.57735027,   0.00000000,   0.81649658},
  { -0.86602540,   0.28867513,  -0.40824829},
  { -0.86602540,  -0.28867513,  -0.40824829},
  { -0.57735027,   0.00000000,  -0.81649658},
  { -0.28867513,   0.86602540,   0.40824829},
  {  0.28867513,   0.86602540,   0.40824829},
  {  0.00000000,   0.57735027,   0.81649658},
  { -0.28867513,   0.86602540,  -0.40824829},
  {  0.28867513,   0.86602540,  -0.40824829},
  {  0.00000000,   0.57735027,  -0.81649658}
};

bool excluded(int i, int *exclude, int nex)
{
  int k;
  if(exclude == NULL) return FALSE;
  
  for(k=0;k<nex;k++){
    if (i==exclude[k]) return TRUE;
  }
  return FALSE;
}


bool find_atom(t_atomlist *al, rvec x, real cut2, int *exclude, int nex)
{
  /* cut2 is the squared cutoff distance */

  rvec diff;
  int i;
  
  for(i=0;i<al->natoms;i++){
    if (strcmp(al->symbol[i],"H") != 0 || excluded(i,exclude,nex)){
      
      if(sqr(al->x[i][XX]-x[XX]) < cut2 &&
         sqr(al->x[i][YY]-x[YY]) < cut2 &&
         sqr(al->x[i][ZZ]-x[ZZ]) < cut2){
        rvec_sub(al->x[i],x,diff);
        if(norm2(diff) < cut2) return TRUE;
      }
    }
  }
  return FALSE;
}

real ray_search(t_atomlist *al,rvec x, int *exclude, int nex)
{
  /* search for atoms along the ray vectors */
  
  real cutoff = .3;
  real cut2; 
  int i,k;
  rvec probe;
  real count = 0;
  
  for(i=0;i<NVEC;i++){
    copy_rvec(x,probe);
    cutoff = .3;
    cut2 = sqr(cutoff);

    if(exclude == NULL && find_atom(al,probe,cut2,exclude,nex)) count+=1.;
    else {
      for(k=1;k<7;k++){
        cutoff+=.3; 
        cut2 =  sqr(cutoff); 
        probe[XX]+=vectors[i][XX];
        probe[YY]+=vectors[i][YY];
        probe[ZZ]+=vectors[i][ZZ];

        if(find_atom(al,probe,cut2,exclude,nex))
        {
          count+=1./sqr((real) k); 
/*           count+=1.; */
          break; /* take next vector */
        }
      }
    }
  }

  
  return count;
}
