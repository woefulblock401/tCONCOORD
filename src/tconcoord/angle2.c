#include <tconcoord.h>

/*========================================================*/
 
int check_angle(rvec at3, rvec at1, rvec at2, t_bounds *b, 
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

/*========================================================*/

