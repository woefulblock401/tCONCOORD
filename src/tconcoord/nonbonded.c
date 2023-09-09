#include <tconcoord.h>
/*======================================================*/
real do_nb_force(t_atomlist *al, t_boundtrack *bt, 
                        real tol, real weight)
{
  /* calculate fake forces for debumping */
  int i,k;
  real vdw,di;
  rvec diff;
  real lb;
  int at;
  real f;
  rvec u;
  real viol = 0;

  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      if(i > al->nb[i][k] && !is_bound(bt,i,al->nb[i][k]))
      {
        at = al->nb[i][k];
        vdw = get_vdw(al,i,at,FALSE);
        lb = vdw - tol;
        rvec_sub(al->x[at],al->x[i],diff);
        
        if(diff[XX] < lb && diff[YY] < lb && diff[ZZ] < lb)
        {
          di = norm(diff);

          if(di < lb && (!al->isposres[i] || !al->isposres[at])){


            
            f = lb-di;
            viol+=weight*f*6;
            u[XX]=diff[XX]/di;
            u[YY]=diff[YY]/di;
            u[ZZ]=diff[ZZ]/di;

          
            if(!al->isposres[i] && !al->isposres[at]){
              al->f[i][XX]-=u[XX]*weight*f;
              al->f[i][YY]-=u[YY]*weight*f;
              al->f[i][ZZ]-=u[ZZ]*weight*f;
              al->f[at][XX]+=u[XX]*weight*f;
              al->f[at][YY]+=u[YY]*weight*f;
              al->f[at][ZZ]+=u[ZZ]*weight*f;
            }
            
            else if(!al->isposres[i] && al->isposres[at]){
              al->f[i][XX]-=u[XX]*weight*f*2;
              al->f[i][YY]-=u[YY]*weight*f*2;
              al->f[i][ZZ]-=u[ZZ]*weight*f*2;
            }
            else if(!al->isposres[at] && al->isposres[i]){
              al->f[at][XX]+=u[XX]*weight*f*2;
              al->f[at][YY]+=u[YY]*weight*f*2;
              al->f[at][ZZ]+=u[ZZ]*weight*f*2;
            }
          }
        }
      }
    }
  }
  return viol;
}
/*======================================================*/



bool count_bumps(t_atomlist *al, t_boundtrack *bt, 
                 real tol, real *bsum, real max, 
                 bool bVerbose)
{

  /* do bump check */

  int i,k;
  real vdw, d;
  rvec diff;
  real lb, ub;
  rvec mid;
  rvec u;
  real worst = 0;
  real dworst = 0.;
  real bump;
  *bsum = 0;
  int at1,at2;
  at1 = at2 = -1;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      if((!al->isposres[i] || !al->isposres[al->nb[i][k]]) 
         && !is_bound(bt,i,al->nb[i][k]))
      {
        vdw = get_vdw(al,i,al->nb[i][k],FALSE);
        rvec_sub(al->x[al->nb[i][k]],al->x[i],diff);
        if(diff[XX] < vdw && diff[YY] < vdw && diff[ZZ] < vdw){
          d = norm(diff);
          if(d < vdw - tol){
            bump = (vdw - tol) - d;
            if(bump > worst) {
              worst = bump;
              dworst = d;
              at1 = i;
              at2 = al->nb[i][k];
            }
            *bsum+=bump;
            
          }
        }
      }
    }
  }
  if(bVerbose && at1!=-1 && at2 !=-1 && bVerbose!=22)
    fprintf(stderr,"Worst bump: %g(%g vs. %g) (%g-%g)  | (%d%s-%s%d) >< (%d%s-%s%d)\n",worst,dworst,get_vdw(al,at1,at2,FALSE),al->vdw[at1],al->vdw[at2],
            al->id[at1],al->name[at1],al->resname[at1],al->resid[at1],
            al->id[at2],al->name[at2],al->resname[at2],al->resid[at2]); 
  
/*   *bsum = *bsum/(real) al->natoms; */
  if(worst > max) {
    return FALSE;
  }
  
  else return TRUE;
}

void bumps_check(FILE *of, t_atomlist *al, t_boundtrack *bt, 
                 real tol)
{

  /* do bump check */


  int i,k;
  bool check;
  real d,vdw;
  rvec diff;
  int at;
  char logstr[STRLEN];
  
  sprintf(logstr,"tCNC__log_> Checking bumps\n");
  CNClog(of,logstr);
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      at = al->nb[i][k];
      
      if((!al->isposres[i] || !al->isposres[at]))
      {
        check = TRUE;
        if(bt==NULL) check = TRUE;
        else 
        {
          check = !is_bound(bt,i,at);
        }
        if(check){
          vdw = get_vdw(al,i,at,FALSE);
          rvec_sub(al->x[at],al->x[i],diff);
          if(diff[XX] < vdw && diff[YY] < vdw && diff[ZZ] < vdw){
            d = norm(diff);
            if(d < vdw - tol){
              sprintf(logstr,"tCNC__log_> BUMP: %d-%s-%d-%s >><< %d-%s-%d-%s : %8.3f (%8.3f) %g %g\n",
                      al->id[i],al->name[i],al->resid[i],al->resname[i],
                      al->id[at],al->name[at],al->resid[at],al->resname[at],
                      d,vdw,al->vdw[i],al->vdw[at]);
              CNClog(of,logstr);
            }
          }
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Done checking bumps\n");
  CNClog(of,logstr);
}

