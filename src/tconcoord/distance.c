#include <tconcoord.h>


void correct_dist(rvec x1, rvec x2, rvec diff, real d,
                  t_bounds *b, int i, int dfunc, gmx_rng_t rng,
                  bool notouch1, bool notouch2)
{

  /* correct distance */

  rvec mid,u;
  real lb, ub, tol, newlen;
  real ran;
  real sig;
  real mean;

  lb = b->lb[i];
  ub = b->ub[i];
  
  mid[XX] = x1[XX]+diff[XX]*.5;
  mid[YY] = x1[YY]+diff[YY]*.5;
  mid[ZZ] = x1[ZZ]+diff[ZZ]*.5;
  
  u[XX]=diff[XX]/d;
  u[YY]=diff[YY]/d;
  u[ZZ]=diff[ZZ]/d;
  
  newlen = 0;

  if(dfunc == 0){
    if(b->isbond[i]){
      tol = (ub-lb)*.5; 
 /*      tol = (ub-lb); */
      ran = gmx_rng_uniform_real(rng);
      tol = tol*ran;
      newlen = lb + tol*1.5;
    }
    else{
      tol = ub-lb;
      ran = gmx_rng_uniform_real(rng);
      tol = tol*ran;
      newlen = lb + tol;
    }  
  }
/*   else if(dfunc == 1){ */
/*     sig = (ub-lb)/4.; */
/*     mean = b->av[i]; */
/*     ran = gmx_rng_gaussian_real(rng); */
/*     tol = sig*tol; */
/*     newlen = mean + tol; */
/*   } */

  if(!notouch1 && !notouch2){
    x1[XX] = mid[XX] - newlen*0.5*u[XX];
    x1[YY] = mid[YY] - newlen*0.5*u[YY];
    x1[ZZ] = mid[ZZ] - newlen*0.5*u[ZZ];
    
    x2[XX] = mid[XX] + newlen*0.5*u[XX];
    x2[YY] = mid[YY] + newlen*0.5*u[YY];
    x2[ZZ] = mid[ZZ] + newlen*0.5*u[ZZ];
    
    
   }
  
  else if(notouch1 && !notouch2){
    x2[XX] = x1[XX]+newlen*u[XX];
    x2[YY] = x1[YY]+newlen*u[YY];
    x2[ZZ] = x1[ZZ]+newlen*u[ZZ];
  }
  else if(notouch2 && !notouch1){
    x1[XX] = x2[XX] - newlen*u[XX];
    x1[YY] = x2[YY] - newlen*u[YY];
    x1[ZZ] = x2[ZZ] - newlen*u[ZZ];
  }
}


