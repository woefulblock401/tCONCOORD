#include <tconcoord.h>

void get_angle_vec(rvec x1,rvec x2, rvec x3, rvec shift)
 {

   rvec v1,v2,dumvec;
   real dum1,dum2;
  
   clear_rvec(v1);
   clear_rvec(v2);
   clear_rvec(shift);
  
   rvec_sub(x2,x1,v1);
   rvec_sub(x3,x1,v2);
  
   dum1 = scalarProduct(v2,v1);
   dum2 = scalarProduct(v2,v2);
  
   v2[XX]*=dum1/dum2;
   v2[YY]*=dum1/dum2;
   v2[ZZ]*=dum1/dum2;

   rvec_sub(v1,v2,dumvec);
   scaleVector(dumvec,shift);
 }
/*========================================================*/

/* void get_angle_vec(rvec x1,rvec x2, rvec x3, rvec shift)  */
/* {  */
 
/*   rvec v1,v2,dumvec;  */
/*   real dum1,dum2;  */
/*   rvec x;  */
/*   rvec cross_x_v2;  */
  
/*   clear_rvec(v1);  */
/*   clear_rvec(v2);  */
/*   clear_rvec(shift);  */
  
/*   rvec_sub(x2,x1,v1);  */
/*   rvec_sub(x3,x1,v2);  */

/*   vectorProduct(v1,v2,x);  */
/*   vectorProduct(v2,x,cross_x_v2);  */
/*   scaleVector(cross_x_v2,shift); */

/* } */

/*========================================================*/
 
int check_angle(rvec at3, rvec at1, rvec at2, t_bounds *b, 
                int i, bool notouch1, bool notouch2) 
{
  /* check angle and apply correction */

  real phi;
  rvec shift1,shift2;
  real bef,aft;
  real scale;
  real viol;
  

  
  phi = RAD2DEG*angle_ij_ik(at3,at1,at2);
/*   bef = phi; */
  
/*   fprintf(stderr,"before = %g\n",phi);  */
  
  if(phi > b->ang[i]+b->sig[i] ||
     phi < b->ang[i]-b->sig[i]){
    get_angle_vec(at3,at2,at1,shift1);
    get_angle_vec(at3,at1,at2,shift2);
    
    if(phi < b->ang[i]-b->sig[i])
    {
      viol = (b->ang[i]-b->sig[i]) - phi;
      scale = viol*.001;
/*       rvec a_b,u; */
/*       rvec_sub(at1,at2,a_b); */
/*       scaleVector(a_b,u); */

      if(!notouch1){
        at1[XX]-=shift1[XX]*scale;
        at1[YY]-=shift1[YY]*scale;
        at1[ZZ]-=shift1[ZZ]*scale;        

/*         at1[0]+=u[0]*0.05; */
/*         at1[1]+=u[1]*0.05; */
/*         at1[2]+=u[2]*0.05; */
      }
      if(!notouch2){
        at2[XX]-=shift2[XX]*scale;
        at2[YY]-=shift2[YY]*scale;
        at2[ZZ]-=shift2[ZZ]*scale;        

/*         at2[0]-=u[0]*0.05; */
/*         at2[1]-=u[1]*0.05; */
/*         at2[2]-=u[2]*0.05; */
      }
      
    }
    
    else if(phi > b->ang[i]+b->sig[i])
    {
      viol = phi-(b->ang[i]+b->sig[i]);
      scale = viol*.001;
      
/*       rvec a_b,u; */
/*       rvec_sub(at1,at2,a_b); */
/*       scaleVector(a_b,u); */
      if(!notouch1){
        at1[XX]+=shift1[XX]*scale;
        at1[YY]+=shift1[YY]*scale;
        at1[ZZ]+=shift1[ZZ]*scale;        
/*         at1[0]-=u[0]*0.05; */
/*         at1[1]-=u[1]*0.05; */
/*         at1[2]-=u[2]*0.05; */
      }
      if(!notouch2){
        at2[XX]+=shift2[XX]*scale;
        at2[YY]+=shift2[YY]*scale;
        at2[ZZ]+=shift2[ZZ]*scale;        
/*         at2[0]+=u[0]*0.05; */
/*         at2[1]+=u[1]*0.05; */
/*         at2[2]+=u[2]*0.05; */
      }
/*       phi = RAD2DEG*angle_ij_ik(at3,at1,at2); */
/*       if(bef < phi) fprintf(stderr,"bef \n"); */
      
    }

/*     phi = RAD2DEG*angle_ij_ik(at3,at1,at2);  */
/*     fprintf(stderr,"after = %g\n",phi);  */
  
    return 0;

  }
  else
    return 1;
}

/*========================================================*/

int check_angle2(rvec at3, rvec at1, rvec at2, t_bounds *b, 
                int i, bool notouch1, bool notouch2) 
{
  /* check angle and apply correction */

  real phi;
  matrix tm1,tm2;
  phi = RAD2DEG*angle_ij_ik(at3,at1,at2);
  
  if(phi > b->ang[i]+b->sig[i] ||
     phi < b->ang[i]-b->sig[i]){
    if(phi < b->ang[i]-b->sig[i])
    {
      rvec a_b,u;
      rvec_sub(at1,at2,a_b);
      scaleVector(a_b,u);

      if(!notouch1){
        at1[0]+=u[0]*0.05;
        at1[1]+=u[1]*0.05;
        at1[2]+=u[2]*0.05;
      }
      if(!notouch2){
        at2[0]-=u[0]*0.05;
        at2[1]-=u[1]*0.05;
        at2[2]-=u[2]*0.05;
      }
      
    }
    
    else if(phi > b->ang[i]+b->sig[i])
    {

      rvec a_b,u;
      rvec_sub(at1,at2,a_b);
      scaleVector(a_b,u);
      if(!notouch1){
        at1[0]-=u[0]*0.05;
        at1[1]-=u[1]*0.05;
        at1[2]-=u[2]*0.05;
      }
      if(!notouch2){
        at2[0]+=u[0]*0.05;
        at2[1]+=u[1]*0.05;
        at2[2]+=u[2]*0.05;
      }

    }
    return 0;
  }
  else
    return 1;
}
