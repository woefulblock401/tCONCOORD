#include <tconcoord.h>


/*============================================================*/
void switch_rvec(rvec x1, rvec x2)
{

  /* this causes somhow trouble when doing target runs */
  rvec dumvec;
  clear_rvec(dumvec);
  
  dumvec[0] = x1[0];
  dumvec[1] = x1[1];
  dumvec[2] = x1[2];
  
  x1[0] = x2[0];
  x1[1] = x2[1];
  x1[2] = x2[2];
  
  x2[0] = dumvec[0];
  x2[1] = dumvec[1];
  x2[2] = dumvec[2];
  

/*   copy_rvec(x1,dum); */
/*   copy_rvec(x2,x1); */
/*   copy_rvec(dum,x2); */
}


/*============================================================*/
void get_impropers(t_atomlist *al, t_dihed *impr)
{
  int i,k;
  real di;

  for(i=0;i<al->natoms;i++){
    if(al->nbonds[i]==4){
      bool doit = FALSE;
      if(al->isposres[i]){
        for(k=0;k<al->nbonds[i];k++){
          if(!al->isposres[al->bonds[i][k]]) doit = TRUE;
        }
      }
      else doit = TRUE;
      if(doit){
        impr->n++;
        impr = dihed_realloc(impr,impr->n);
        impr->center[impr->n-1] = i;
        impr->at1[impr->n-1] = al->bonds[i][0];
        impr->at2[impr->n-1] = al->bonds[i][1];
        impr->at3[impr->n-1] = al->bonds[i][2];
        impr->at4[impr->n-1] = al->bonds[i][3];
        impr->phi[impr->n-1] = DIHED(al,al->bonds[i][0],
                                     al->bonds[i][1],
                                     al->bonds[i][2],
                                     al->bonds[i][3]);
      }
    }
  }
}
/*============================================================*/
void print_improps(t_atomlist *al, t_dihed *impr)
{
  int i;
  for(i=0;i<impr->n;i++){
    printf("center = %i%s-%s%d\n",al->id[impr->center[i]],al->name[impr->center[i]],al->resname[impr->center[i]],
           al->resid[impr->center[i]]);
    printf("at1 = %i%s-%s%d\n",al->id[impr->at1[i]],al->name[impr->at1[i]],al->resname[impr->at1[i]],
           al->resid[impr->at1[i]]);
    printf("at2 = %i%s-%s%d\n",al->id[impr->at2[i]],al->name[impr->at2[i]],al->resname[impr->at2[i]],
           al->resid[impr->at2[i]]);
    printf("at3 = %i%s-%s%d\n",al->id[impr->at3[i]],al->name[impr->at3[i]],al->resname[impr->at3[i]],
           al->resid[impr->at3[i]]);
    printf("at4 = %i%s-%s%d\n",al->id[impr->at4[i]],al->name[impr->at4[i]],al->resname[impr->at4[i]],
           al->resid[impr->at4[i]]);
  }
}
/*============================================================*/
void swap_sidechain(t_atomlist *al, t_resl *rl, int center, 
                    int *idx, int *tags)
{
  /* swap the coordinates of 
     idx[2] and idx[3] and 
     all their bonds
  */
  int i,k;
  int res;
  rvec diff;
  
  if(al->resid[idx[2]] != al->resid[idx[3]]){
    fatal_error("Function swap_sidechain\n");
  }

  res = al->resid[idx[2]]-1;
  
  if(strcmp(al->symbol[idx[2]],"H") == 0 &&
     strcmp(al->symbol[idx[3]],"H") == 0){
    switch_rvec(al->x[idx[2]],al->x[idx[3]]);
    tags[idx[2]] = TRUE;
    tags[idx[3]] = TRUE;
  }
  else{
    rvec_sub(al->x[idx[3]],al->x[idx[2]],diff);
/*     write_pdb(al,"before.pdb"); */
        
 /*    real phi = ANGLE(al,center,idx[2],idx[3]); */
    real phi = angle_ij_ik(al->x[center],al->x[idx[2]],al->x[idx[3]]);
    
/*     printf("ang = %g\n",phi); */
    rvec vec1, vec2, rotvec, s_rotvec;
    rvec_sub(al->x[center],al->x[idx[2]],vec1);
    rvec_sub(al->x[center],al->x[idx[3]],vec2);
    vectorProduct(vec1, vec2, rotvec);
    scaleVector(rotvec,s_rotvec);
    matrix tm1,tm2;
    createRotationMatrix1(s_rotvec,tm1);
    createRotationMatrix2(s_rotvec,tm2);
    
    switch_rvec(al->x[idx[2]],al->x[idx[3]]);   


     for(k=rl->j0[res];k<=rl->j1[res];k++){ 
       if(al->order[k] > al->order[idx[2]] && 
          al->name[k][3] == al->name[idx[2]][3]) 
       { 
         rotateVectorAroundVector(al->x[k],al->x[center],tm1, 
                                  tm2,phi); 
         tags[k] = TRUE; 
       } 
       else if(al->order[k] > al->order[idx[3]] && 
               al->name[k][3] == al->name[idx[3]][3]) 
       { 
         rotateVectorAroundVector(al->x[k],al->x[center],tm1, 
                                  tm2,-phi); 
         tags[k] = TRUE; 
       } 
     } 

    
  }
}
/*============================================================*/

int check_impr(t_atomlist *al, t_idxgroups *imp, int *imps)
{
  int k,l;
  int count = 0;
  real dihed;
  bool doit;
    
  for(k=0;k<imp->n;k++){
    imps[k] = FALSE;
    dihed = DIHED(al,imp->atoms[k][0],imp->atoms[k][1],
                  imp->atoms[k][2],imp->atoms[k][3]);
/*     printf("k = %d dihed = %g\n",k,dihed); */
    if(dihed * imp->val[k] < 0.){
      count++;
      imps[k] = TRUE;
    }
  }
  return count;
}

/*============================================================*/
int do_improp2(t_atomlist *al, t_dihed *impr, t_resl *rl,int *tags)
{
  int i,k;
  int count = 0;
  int idx[4];
  real order[4];
  
  real dihed;
  for(i=0;i<impr->n;i++){
    dihed = DIHED(al,impr->at1[i],
                  impr->at2[i],
                  impr->at3[i],
                  impr->at4[i]);
    if(dihed*impr->phi[i] < 0){
      idx[0] = impr->at1[i];
      idx[1] = impr->at2[i];
      idx[2] = impr->at3[i];
      idx[3] = impr->at4[i];
      order[0] = al->order[impr->at1[i]];
      order[1] = al->order[impr->at2[i]];
      order[2] = al->order[impr->at3[i]];
      order[3] = al->order[impr->at4[i]];
      array_sort(idx,order,4);
  /*     printf("%s- %s -%s -%s\n",al->name[idx[0]],al->name[idx[1]],al->name[idx[2]],al->name[idx[3]]); */
      
  
      swap_sidechain(al,rl, impr->center[i],idx, tags); 
      count++;
      
    }
  }
  return count;
  
}
/*============================================================*/

/*============================================================*/

int do_improp(t_atomlist *al, t_idxgroups *imp, int *tags, 
               gmx_rng_t rng)
{
  
  int at1,at2,at3,at4;
  real dihed;
  int i;
  rvec o,p,v1,v2,v3,v4,pp;
  real l1,l2,i12;
  real rk,rl,a1,a2;
  real p3,p4;
  int nimp = 0;
  
  for(i=0;i<imp->n;i++){
    at1 = imp->atoms[i][0];
    at2 = imp->atoms[i][1];
    at3 = imp->atoms[i][2];
    at4 = imp->atoms[i][3];
    dihed = DIHED(al,at1,at2,at3,at4);
    if(dihed * imp->val[i] < 0){
      nimp++;
      copy_rvec(al->x[at1],o);
      copy_rvec(al->x[at4],p);
      
      rvec_sub(al->x[at2],al->x[at1],v1);
      rvec_sub(al->x[at3],al->x[at1],v2);
      
      l1 = DIST(al,at1,at2);
      l2 = DIST(al,at1,at3);
      i12 = iprod(v1,v2);
      
      rk = -i12/sqr(l1);
      rl = sqr(rk)*sqr(l1)+sqr(l2)+2.*rk*i12;
      
      a2 = sqrt(rl)/rl;
      a1 = rk*a2;
      
      v3[XX] = v1[XX]/l1;
      v3[YY] = v1[YY]/l1;
      v3[ZZ] = v1[ZZ]/l1;
      
      v4[XX] = a1*v1[XX]+a2*v2[XX];
      v4[YY] = a1*v1[YY]+a2*v2[YY];
      v4[ZZ] = a1*v1[ZZ]+a2*v2[ZZ];
      
      p3 = (p[XX]-o[XX])*v3[XX]+
        (p[YY]-o[YY])*v3[YY]+
        (p[ZZ]-o[ZZ])*v3[ZZ];
      
      p4 = (p[XX]-o[XX])*v4[XX]+
        (p[YY]-o[YY])*v4[YY]+
        (p[ZZ]-o[ZZ])*v4[ZZ];
      
      pp[XX] = o[XX]+p3*v3[XX]+p4*v4[XX];
      pp[YY] = o[YY]+p3*v3[YY]+p4*v4[YY];
      pp[ZZ] = o[ZZ]+p3*v3[ZZ]+p4*v4[ZZ];
      
      al->x[at4][XX] = 2*pp[XX]-p[XX];
      al->x[at4][YY] = 2*pp[YY]-p[YY];
      al->x[at4][ZZ] = 2*pp[ZZ]-p[ZZ];
      tags[at4] = TRUE;
      
    }
  }
  return nimp;
}


/*============================================================*/

int do_chiral(t_atomlist *al, t_idxgroups *imp, int *tags)
{
  int i,k;
  real dihed;
  int at1, at2, at3, at4;
  int atom_idx[3];
  rvec cm;
  rvec new[3];
  rvec diff;
  int nviol = 0;
  
  for(i=0;i<imp->n;i++){
    at1 = imp->atoms[i][0];
    at2 = imp->atoms[i][1];
    at3 = imp->atoms[i][2];
    at4 = imp->atoms[i][3];
    dihed = DIHED(al,at1,at2,at3,at4);
    if(dihed*imp->val[i]<0){
/*       printf("dihed before = %g (%g)\n", dihed, imp->val[i]); */
      
      nviol++;
      /* correct chirality */
      atom_idx[0] = at2;
      atom_idx[1] = at3;
      atom_idx[2] = at4;
      com_frag(al,3,atom_idx,new,cm);
      rvec_sub(cm,al->x[at1],diff);
      al->x[at1][XX]+=diff[XX];
      al->x[at1][YY]+=diff[YY];
      al->x[at1][ZZ]+=diff[ZZ];
      for(k=1;k<4;k++) {
        al->x[imp->atoms[i][k]][XX]-=diff[XX];
        al->x[imp->atoms[i][k]][YY]-=diff[YY];
        al->x[imp->atoms[i][k]][ZZ]-=diff[ZZ];
      }
      tags[at1] = TRUE;
      tags[at2] = TRUE;
      tags[at3] = TRUE;
      tags[at4] = TRUE;
      dihed = DIHED(al,at1,at2,at3,at4);
/*       printf("dihed after = %g\n", dihed); */

    }

  }
  return nviol;
}



void mirror(t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    al->x[i][ZZ] = -al->x[i][ZZ];
  }
}

/*============================================================*/

real do_chiral_force(t_atomlist *al, t_idxgroups *imp, real weight)
{
  int i,k;
  real dihed;
  int at1, at2, at3, at4;
  int atom_idx[3];
  rvec cm;
  rvec new[3];
  rvec diff;
  real nviol = 0;
  
  for(i=0;i<imp->n;i++){
    at1 = imp->atoms[i][0];
    at2 = imp->atoms[i][1];
    at3 = imp->atoms[i][2];
    at4 = imp->atoms[i][3];
    dihed = DIHED(al,at1,at2,at3,at4);
    if(dihed*imp->val[i]<0){

      
      /* correct chirality */
      atom_idx[0] = at2;
      atom_idx[1] = at3;
      atom_idx[2] = at4;
      com_frag(al,3,atom_idx,new,cm);
      rvec_sub(cm,al->x[at1],diff);
/*       printf("correcting chirality with %g\n", norm(diff)*weight); */
      nviol+=norm(diff)*weight;
      al->f[at1][XX]+=diff[XX]*weight;
      al->f[at1][YY]+=diff[YY]*weight;
      al->f[at1][ZZ]+=diff[ZZ]*weight;
      for(k=1;k<4;k++) {
        al->f[imp->atoms[i][k]][XX]-=diff[XX]*weight;
        al->f[imp->atoms[i][k]][YY]-=diff[YY]*weight;
        al->f[imp->atoms[i][k]][ZZ]-=diff[ZZ]*weight;
      }
    }

  }
  return nviol;

}

