#include <tconcoord.h>


/*================================================================*/

void dump_forces(t_atomlist *al, char *filename)
{
  int i;
  rvec f;
  real fo;
  
  real max = 0;
  
  FILE *fp = fopen(filename,"w");
  for(i=0;i<al->natoms;i++){
    fo = norm(al->f[i]);
    if (fo>max) max = fo;
    fprintf(fp,"%d%s%d%s  %g  [%g,%g,%g]\n",al->id[i],al->name[i],al->resid[i],al->resname[i],
            fo,al->f[i][XX],al->f[i][YY],al->f[i][ZZ]);
  }
  fclose(fp);
  
}
/*================================================================*/

void apply_f(t_atomlist *al, real weight, int *tags)
{
  /* apply the forces to atoms */

  int i;
  real di;
  int n = 0;
  real sum = 0;
  
  for(i=0;i<al->natoms;i++){
    if(al->isposres[i]) clear_rvec(al->f[i]);
    else{
      di = norm(al->f[i]); 
/*       printf("Moving %g\n",di); */
      
/*       if(di > .1){ */
/*         al->f[i][XX]/=di; */
/*         al->f[i][YY]/=di; */
/*         al->f[i][ZZ]/=di; */
/*       } */
      
      if(di > 0.){
        n++;
        al->x[i][XX]+=al->f[i][XX]*weight;
        al->x[i][YY]+=al->f[i][YY]*weight;
        al->x[i][ZZ]+=al->f[i][ZZ]*weight;
/*         sum+=norm(al->f[i])*weight; */
        
        tags[i] = TRUE;
      }
/*     else tags[i] = FALSE; */
    }
  }
/*   fprintf(stderr,"Moved %d of %d atoms by %g A\n",n,al->natoms,sum); */
  
}
/*================================================================*/
void apply_f2(t_atomlist *al, real weight, int *tags, real max)
{
  /* apply the forces to atoms */

  int i;
  real di;

  for(i=0;i<al->natoms;i++){
    if(al->isposres[i]) clear_rvec(al->f[i]);
    else{
      di = norm(al->f[i]); 
      if (di > max){
        
        al->f[i][XX]/=di * max; 
        al->f[i][YY]/=di * max; 
        al->f[i][ZZ]/=di* max; 
      }
      
      if(di > 0.){
        al->x[i][XX]+=al->f[i][XX]*weight;
        al->x[i][YY]+=al->f[i][YY]*weight;
        al->x[i][ZZ]+=al->f[i][ZZ]*weight;
        tags[i] = TRUE;
      }
/*     else tags[i] = FALSE; */
    }
  }
}
/*================================================================*/
void clear_forces(t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    clear_rvec(al->f[i]);
  }
}
/*================================================================*/
real sum_forces(t_atomlist *al)
{
  int i;
  real sum = 0.;
  for(i=0;i<al->natoms;i++){
    sum+=norm(al->f[i]);
  }
  return sum;
}

/*================================================================*/
void backbone_opt(t_atomlist *al, int *nt, int *tags)
{
  int i,k,j,l;
  real dij,eij,rij,r;
  rvec u,diff;
  bool bDec,bInc;
  real atr,rep,ener;
  real max  = -1000;
  int imax,jmax;
  imax = 0;
  
  real en;
  
  ener = 0;
  
  for(i=0;i<al->natoms;i++){
    if(al->order[i] < 2){
      atr = rep = 0;
      for(k=0;k<al->nb14[i];k++){
        bInc = bDec = FALSE;
        j = al->b14[i][k];
        if(al->order[j] < 2){
          rvec_sub(al->x[i],al->x[j],diff);
          dij = norm(diff);
          rij = al->rs_comb[al->rs_type[i]][al->rs_type[j]]; 
          eij = sqrt(al->rs_eps[i]*al->rs_eps[j]);
          
          r = rij / dij *.9;
          if(r < 1.12) atr += (pow(r,12) - 2*pow(r,6)) * eij;
          else{
            en = (pow(r,12) - 2*pow(r,6)) * eij;
            rep+=en;
     /*        if(en > max) { */
/*               max = en; */
/*               imax = i; */
/*               jmax = j; */
/*             } */
            
          }
          u[XX]=diff[XX]/dij;
          u[YY]=diff[YY]/dij;
          u[ZZ]=diff[ZZ]/dij;
          
        }
      }
      en = rep+atr;
      if(en > max){
        imax = i;
        max = en;
      }
      ener+=en;
    }
    
  }
  printf("imax = %d\n",imax);
  
/*   printf("lj-energy: %s(%d%s) = %g\n",al->name[i],al->resid[i],al->resname[i],rep+atr);  */
  printf("lj-energy  = %g max(%s/%d%s) = %g \n",ener,al->name[imax],al->resid[imax],al->resname[imax],max); 
}

