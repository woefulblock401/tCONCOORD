#include <tconcoord.h>

/*===============================================================*/

real do_angle_force(t_atomlist *al, t_bounds *b, real weight)
{
  int i,k;
  real viol = 0;
  bool bDec = FALSE;
  bool bInc = FALSE;
  real phi;
  real di;
  rvec diff,shift1,shift2;
  real lb;
  int at;
  real f = 0.;
  real tot = 0;
  real scale = 0;
  bool bDDec = FALSE;
  bool bDInc = FALSE;
  int nviol = 0;
  
  
  for(i=0;i<b->n;i++){
    if(b->isang[i]){
      rvec_sub(al->x[b->at2[i]],al->x[b->at1[i]],diff);
      bDec = bInc = FALSE;
      di = norm2(diff);
      if(di < sqr(b->ub[i]) && di > sqr(b->lb[i])) {
        phi = RAD2DEG*angle_ij_ik(al->x[b->at3[i]],al->x[b->at1[i]],al->x[b->at2[i]]);
        if(phi > b->ang[i]+b->sig[i]){
          viol = phi-(b->ang[i]+b->sig[i]);
          scale = viol;
/*         fprintf(stderr,"scale = %g\n",scale); */
          
          bDec = TRUE;
/*         f = 0.01; */
        }
        else if(phi < b->ang[i]-b->sig[i]){
          viol = (b->ang[i]-b->sig[i]) - phi;
          scale = viol;
          bInc = TRUE;
        }
        
        if(bDec || bInc) {
          nviol++;
          
          get_angle_vec(al->x[b->at3[i]],al->x[b->at2[i]],al->x[b->at1[i]],shift1);
          get_angle_vec(al->x[b->at3[i]],al->x[b->at1[i]],al->x[b->at2[i]],shift2);        
          tot+=6*weight*scale;
/*           fprintf(stderr,"phi = %g viol = %g di = %g corr = %g shift = %g\n",phi,viol,di, weight*scale);     */
/*           printf("ang force = %g\n",weight*norm(shift1)*scale); */
          
          if(!al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
            if(bDec){
/*           fprintf(stderr,"af = %g\n",weight*f*.5); */
              al->f[b->at1[i]][XX] += weight*shift1[XX]*scale*.5;
              al->f[b->at1[i]][YY] += weight*shift1[YY]*scale*.5;
              al->f[b->at1[i]][ZZ] += weight*shift1[ZZ]*scale*.5;
              al->f[b->at2[i]][XX] += weight*shift2[XX]*scale*.5;
              al->f[b->at2[i]][YY] += weight*shift2[YY]*scale*.5;
              al->f[b->at2[i]][ZZ] += weight*shift2[ZZ]*scale*.5;
            }
            else if(bInc) {
              al->f[b->at1[i]][XX] -= weight*shift1[XX]*scale*.5;
              al->f[b->at1[i]][YY] -= weight*shift1[YY]*scale*.5;
              al->f[b->at1[i]][ZZ] -= weight*shift1[ZZ]*scale*.5;
              al->f[b->at2[i]][XX] -= weight*shift2[XX]*scale*.5;
              al->f[b->at2[i]][YY] -= weight*shift2[YY]*scale*.5;
              al->f[b->at2[i]][ZZ] -= weight*shift2[ZZ]*scale*.5;
            }
          }
          else if(al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
            if(bDec){
              al->f[b->at2[i]][XX] += weight*shift2[XX]*scale;
              al->f[b->at2[i]][YY] += weight*shift2[YY]*scale;
              al->f[b->at2[i]][ZZ] += weight*shift2[ZZ]*scale;
            }
            else if(bInc) {
              al->f[b->at2[i]][XX] -= weight*shift2[XX]*scale;
              al->f[b->at2[i]][YY] -= weight*shift2[YY]*scale;
              al->f[b->at2[i]][ZZ] -= weight*shift2[ZZ]*scale;
            }
          }
          else if(!al->isposres[b->at1[i]] && al->isposres[b->at2[i]]){
            if(bDec){
              al->f[b->at1[i]][XX] += weight*shift1[XX]*scale;
              al->f[b->at1[i]][YY] += weight*shift1[YY]*scale;
              al->f[b->at1[i]][ZZ] += weight*shift1[ZZ]*scale;
            }
            else if(bInc) {
              al->f[b->at1[i]][XX] -= weight*shift1[XX]*scale;
              al->f[b->at1[i]][YY] -= weight*shift1[YY]*scale;
              al->f[b->at1[i]][ZZ] -= weight*shift1[ZZ]*scale;
            }
          }
        }
      }
    }
  }
/*   fprintf(stderr,"f (angles) per constraint: %g\n",tot/(real)nviol);   */
  return tot;
  
}




real do_bound_force(t_atomlist *al, t_bounds *b, 
                    real bweight, real aweight,real dweight,
                    gmx_rng_t rng) 
{ 
  /* calculate fake forces for correction */ 

  int i,k;
  real di;
  rvec diff;
  real lb;
  int at;
  real f = 0.;
  rvec u;
  real viol = 0;
  real weight = 0.2;
  bool bDec = FALSE;
  bool bInc = FALSE;
  real phi;
  int nviol = 0;
  real r = 0;
  
/*   af = do_angle_force(al,b,aweight);  */

  for(i=0;i<b->n;i++){
    if(b->isbond[i]  || b->isang[i]  || b->isdih[i]) {
      f = 0.;
      rvec_sub(al->x[b->at2[i]],al->x[b->at1[i]],diff);
      bDec = bInc = FALSE;
      di = norm(diff);

      if(b->isbond[i]) weight = bweight;
      else if(b->isang[i]) weight = aweight;
      else if(b->isdih[i]) weight = dweight;
      

      if(di > b->ub[i]) {
        bDec = TRUE;
        f = di - b->ub[i];
/*         fprintf(stderr,"f_cons = %g\n",f); */
        
/*         r = gmx_rng_gaussian_real(rng); */
/*         r*=f*.05; */
/*         f+=r; */
/*         fprintf(stderr,"f_app = %g\n",f); */
      } 
      else if(di < b->lb[i]) {
        bInc = TRUE;
        f = b->lb[i] - di;
/*         r = gmx_rng_gaussian_real(rng); */
/*         r*=f*.05; */
/*         f+=r; */
      }
      



      if(bDec || bInc) {
        nviol++;
        
        u[XX]=diff[XX]/di; 
        u[YY]=diff[YY]/di; 
        u[ZZ]=diff[ZZ]/di; 
        viol+=6*weight*f;
        
        if(!al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
          if(bDec){
/*             if(b->isang[i]) */
/*               fprintf(stderr,"f = %g\n",weight*f*.5); */
              
             al->f[b->at1[i]][XX] += weight*u[XX]*f*0.5; 
             al->f[b->at1[i]][YY] += weight*u[YY]*f*0.5; 
             al->f[b->at1[i]][ZZ] += weight*u[ZZ]*f*0.5; 
            
             al->f[b->at2[i]][XX] -= weight*u[XX]*f*0.5; 
             al->f[b->at2[i]][YY] -= weight*u[YY]*f*0.5; 
             al->f[b->at2[i]][ZZ] -= weight*u[ZZ]*f*0.5; 
          }
          else if(bInc) {
             al->f[b->at1[i]][XX] -= weight*u[XX]*f*0.5; 
             al->f[b->at1[i]][YY] -= weight*u[YY]*f*0.5; 
             al->f[b->at1[i]][ZZ] -= weight*u[ZZ]*f*0.5; 
            
             al->f[b->at2[i]][XX] += weight*u[XX]*f*0.5; 
             al->f[b->at2[i]][YY] += weight*u[YY]*f*0.5; 
             al->f[b->at2[i]][ZZ] += weight*u[ZZ]*f*0.5; 
          }
        }
        else if(al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
          if(bDec){

             al->f[b->at2[i]][XX] -= weight*u[XX]*f; 
             al->f[b->at2[i]][YY] -= weight*u[YY]*f; 
             al->f[b->at2[i]][ZZ] -= weight*u[ZZ]*f; 
          }
          else if(bInc) {
             al->f[b->at2[i]][XX] += weight*u[XX]*f; 
             al->f[b->at2[i]][YY] += weight*u[YY]*f; 
             al->f[b->at2[i]][ZZ] += weight*u[ZZ]*f; 
          }
        }
        else if(!al->isposres[b->at1[i]] && al->isposres[b->at2[i]]){
          if(bDec){

             al->f[b->at1[i]][XX] += weight*u[XX]*f; 
             al->f[b->at1[i]][YY] += weight*u[YY]*f; 
             al->f[b->at1[i]][ZZ] += weight*u[ZZ]*f; 
          }
          else if(bInc) {
             al->f[b->at1[i]][XX] -= weight*u[XX]*f; 
             al->f[b->at1[i]][YY] -= weight*u[YY]*f; 
             al->f[b->at1[i]][ZZ] -= weight*u[ZZ]*f; 
          }
        }
        
      }
    }
  }
/*   fprintf(stderr,"f per constraint: %g\n",viol/(real)nviol); */
  
  return viol;
}
/*===============================================================*/

real do_bound_force2(t_atomlist *al, t_bounds *b, 
                     real weight)
{ 
  /* calculate fake forces for correction */ 

  int i,k;
  real di;
  rvec diff;
  real lb;
  int at;
  real f = 0.;
  rvec u;
  real viol = 0;
  bool bDec = FALSE;
  bool bInc = FALSE;
  real phi;
  

  for(i=0;i<b->n;i++){
    
    if(!b->isbond[i] && !b->isang[i] && !b->isdih[i]) {
      f = 0.;
      
      rvec_sub(al->x[b->at2[i]],al->x[b->at1[i]],diff);
      bDec = bInc = FALSE;
      di = norm(diff);

      if(di > b->ub[i]) {
        bDec = TRUE;
        f = di - b->ub[i];
/*         viol += f; */
      } 
      else if(di < b->lb[i]) {
        bInc = TRUE;
        f = b->lb[i] - di;
/*         viol += f; */
      }
      
      if(bDec || bInc) {
        
        u[XX]=diff[XX]/di;
        u[YY]=diff[YY]/di;
        u[ZZ]=diff[ZZ]/di;
        viol+=6*weight*f;
        
        if(!al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
          if(bDec){
            al->f[b->at1[i]][XX] += weight*u[XX]*f*0.5;
            al->f[b->at1[i]][YY] += weight*u[YY]*f*0.5;
            al->f[b->at1[i]][ZZ] += weight*u[ZZ]*f*0.5;
            
            al->f[b->at2[i]][XX] -= weight*u[XX]*f*0.5;
            al->f[b->at2[i]][YY] -= weight*u[YY]*f*0.5;
            al->f[b->at2[i]][ZZ] -= weight*u[ZZ]*f*0.5;
          }
          else if(bInc) {
            al->f[b->at1[i]][XX] -= weight*u[XX]*f*0.5;
            al->f[b->at1[i]][YY] -= weight*u[YY]*f*0.5;
            al->f[b->at1[i]][ZZ] -= weight*u[ZZ]*f*0.5;
            
            al->f[b->at2[i]][XX] += weight*u[XX]*f*0.5;
            al->f[b->at2[i]][YY] += weight*u[YY]*f*0.5;
            al->f[b->at2[i]][ZZ] += weight*u[ZZ]*f*0.5;
          }
        }
        else if(al->isposres[b->at1[i]] && !al->isposres[b->at2[i]]){
          if(bDec){
            al->f[b->at2[i]][XX] -= weight*u[XX]*f;
            al->f[b->at2[i]][YY] -= weight*u[YY]*f;
            al->f[b->at2[i]][ZZ] -= weight*u[ZZ]*f;
          }
          else if(bInc) {
            al->f[b->at2[i]][XX] += weight*u[XX]*f;
            al->f[b->at2[i]][YY] += weight*u[YY]*f;
            al->f[b->at2[i]][ZZ] += weight*u[ZZ]*f;
          }
        }
        else if(!al->isposres[b->at1[i]] && al->isposres[b->at2[i]]){
          if(bDec){
            al->f[b->at1[i]][XX] += weight*u[XX]*f;
            al->f[b->at1[i]][YY] += weight*u[YY]*f;
            al->f[b->at1[i]][ZZ] += weight*u[ZZ]*f;
          }
          else if(bInc) {
            al->f[b->at1[i]][XX] -= weight*u[XX]*f;
            al->f[b->at1[i]][YY] -= weight*u[YY]*f;
            al->f[b->at1[i]][ZZ] -= weight*u[ZZ]*f;
          }
        }
        
      }
    }
  }
/*   exit(0); */
  
  return viol;
}
/*===============================================================*/

