#include <tconcoord.h>
/*===========================================================*/
bool check_planar(t_atomlist *al, int n, int *idx, real tol)
{
  int i,m,k;
  rvec x[n];
  matrix mat,trans;
  rvec dd;
  rvec co;
  bool ret = FALSE;

  
  clear_mat(mat);
  clear_mat(trans);

  for(i=0;i<n;i++){
    copy_rvec(al->x[idx[i]],x[i]);
  }
  com_frag(al,n,idx,x,co);
  clear_rvec(dd);
  princ_comp(n,x,mat,dd);
  if (det(mat) < 0) {
    for(m=0; (m<DIM); m++)
      mat[ZZ][m] = -mat[ZZ][m];
  }

/*   for(i=0;i<DIM;i++){ */
/*     printf("%8.3f %8.3f %8.3f\n",mat[i][XX],mat[i][YY],mat[i][ZZ]); */
/*   } */
/*   printf("\n"); */
  

  
  rotate_atoms(n,x,mat);
  real adev = 0.;
  for(k=0;k<n;k++){
    adev+=sqrt(sqr(x[k][0]));
  }
  adev/=(real)n;

  if(adev > tol){
    ret  = TRUE;
    /* make it flat */
    for(k=0;k<n;k++){
      x[k][0] = 0.;
    }
    matrix trans;
    transpose(mat,trans);

    rotate_atoms(n,x,trans);

    for(i=0;i<n;i++){
      copy_rvec(x[i],al->x[idx[i]]);
      al->x[idx[i]][0]+=co[0];
      al->x[idx[i]][1]+=co[1];
      al->x[idx[i]][2]+=co[2];
    }
  }
  return ret;
}

/*===========================================================*/


bool count_planar(t_atomlist *al, t_idxgroups *pl, real *sum,
                  int *count, real max, real *pworst)
{
  int i,m,k,n,l;
  
  matrix mat,trans;
  rvec dd;
  rvec co;
  *sum = 0;
  *count = 0;
  real dev[pl->n];
  *pworst = 0;
  real ptol;
  real adev;
    
  for(k=0;k<pl->n;k++){ 

    ptol = pl->val[k];
    is_planar_group(al,pl->natoms[k],pl->atoms[k],ptol,&adev);

    if(adev > ptol){
      *sum+=adev - ptol;
      dev[k] = (adev/ptol-1)*100;
      *count+=1;
      
      if(dev[k] > *pworst ) *pworst = dev[k];
      
    }
    if(adev > ptol) *count+=1;
  }
/*   if(*pworst < max && *count < pl->n/100.) */
  if(*pworst < max)
    return TRUE;
  else return FALSE;

}
/*===========================================================*/

real do_planar_force(t_atomlist *al, int n, int *idx, real tol, real weight)
{
  int i,m,k;
  rvec x[n];
  rvec f[n];
  matrix mat,trans;
  rvec dd;
  rvec co;
  bool ret = FALSE;
  real force = 0;
  
  
  clear_mat(mat);
  clear_mat(trans);
/*   printf("nat = %d\n",n); */
  for(i=0;i<n;i++){
    copy_rvec(al->x[idx[i]],x[i]);
    copy_rvec(al->f[idx[i]],f[i]);
  }
  com_frag(al,n,idx,x,co);
  clear_rvec(dd);
  princ_comp(n,x,mat,dd);
/*   printf("after princ\n"); */
  if (det(mat) < 0) {
    for(m=0; (m<DIM); m++)
      mat[ZZ][m] = -mat[ZZ][m];
  }
  
  rotate_atoms(n,x,mat);
  rotate_atoms(n,f,mat);

  real adev = 0.;
  for(k=0;k<n;k++){
    adev+=sqrt(sqr(x[k][0]));
  }
  adev/=(real)n;
  if(adev > tol){
    ret  = TRUE;
    /* make it flat */
    for(k=0;k<n;k++){
/*       x[k][0] = 0.; */
      f[k][0] = -x[k][0];
      force+=sqrt(sqr(f[k][0]));
    }
    matrix trans;
    transpose(mat,trans);

    rotate_atoms(n,f,trans);

    for(i=0;i<n;i++){

      al->f[idx[i]][XX]+=f[i][XX]*weight;
      al->f[idx[i]][YY]+=f[i][YY]*weight;
      al->f[idx[i]][ZZ]+=f[i][ZZ]*weight;
    }
  }
  return force;
  
}

/*===========================================================*/

/* bool is_omega_group(t_atomlist *al, int *idx, int n) */
/* { */
/*   int i; */
/*   bool have_ca = FALSE; */
/*   bool have_o = FALSE; */
/*   bool have_n = FALSE; */
/*   bool have_c = FALSE; */
  
/*   for(i=0;i<n;i++){ */
/*     if(strcmp(al->name[idx[i]]," CA ") == 0) have_ca = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," N  ") == 0) have_n = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," O  ") == 0) have_o = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," C  ") == 0) have_c = TRUE; */
/*   } */
/*   if(have_ca && have_o && have_n && have_c) return TRUE; */
/*   else return FALSE; */
/* } */
 /*=============================================================*/ 
/* bool is_arg_head_group(t_atomlist *al, int *idx, int n) */
/* { */
/*   int i; */
/*   bool have_cd = FALSE; */
/*   bool have_ne = FALSE; */
/*   bool have_cz = FALSE; */
/*   bool have_nh1 = FALSE; */
/*   bool have_nh2 = FALSE; */
  
/*   for(i=0;i<n;i++){ */
/*     if(strcmp(al->name[idx[i]]," CD ") == 0) have_cd = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," NE ") == 0) have_ne = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," CZ ") == 0) have_cz = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," NH1") == 0) have_nh1 = TRUE; */
/*     else if(strcmp(al->name[idx[i]]," NH2") == 0) have_nh2 = TRUE; */
/*   } */
/*   if(have_cd && have_ne && have_cz && have_nh1 && have_nh2) return TRUE; */
/*   else return FALSE; */
/* } */

/*=============================================================*/

bool is_planar_group(t_atomlist *al,int n,int *atoms, 
                     const real tol, real *plan)

/* check whether a group of atoms is planar */

{
  int k,m;
  rvec x[n];
  for(k=0;k<n;k++){
    x[k][XX] = al->x[atoms[k]][XX];
    x[k][YY] = al->x[atoms[k]][YY];
    x[k][ZZ] = al->x[atoms[k]][ZZ];
  }
  matrix mat;
  rvec dd;
  rvec co;
  com_frag(al,n,atoms,x,co);
  clear_mat(mat);
  clear_rvec(dd);
  princ_comp(n,x,mat,dd);
  
  if (det(mat) < 0) {
    for(m=0; (m<DIM); m++)
      mat[ZZ][m] = -mat[ZZ][m];
  }
  rotate_atoms(n,x,mat);

  real adev = 0.;
  for(k=0;k<n;k++){
    adev+=sqrt(sqr(x[k][XX]));
  }
  *plan = adev/(real)n;
  if(*plan < tol) return TRUE;
  return FALSE;
}
/*=============================================================*/

int get_planar_groups(FILE *log,t_atomlist *al, t_resl *rl, 
                       t_idxgroups *pln, bool bIgn)
{
  int i,k;
  int j,l,m;
  int n;
  real omtol = 0;
  real plan;
  char wstr[STRLEN];
  int plviol = 0;
  
  FILE *fp = cnclib("Planar.dat");
  
  t_namegroups *ng = read_namegroups(fp);
  
  for(i=0;i<rl->nres;i++){
    for(k=0;k<ng->n;k++){
      if(strcmp(rl->resname[i],ng->resname[k]) == 0){
        pln->n++;
        pln = idx_realloc(pln,pln->n);
        n = pln->n-1;
        pln->val[n] = ng->val[k];
        
        for(j=rl->j0[i];j<=rl->j1[i];j++){
          for(l=0;l<ng->natoms[k];l++){
            if(strcmp(al->name[j],ng->atomnames[k][l]) == 0)
            {
              add_to_group(pln,n,j);
              for(m=0;m<al->nbonds[j];m++){
/*               printf("adding %d%s\n",al->resid[al->bonds[j][m]],al->name[al->bonds[j][m]]); */

                add_to_group(pln,n,al->bonds[j][m]);
              }
            }
          }
        }
        if(pln->natoms[n] < 4){
          sprintf(wstr,"Group has less than 4 atoms -> cannot create planar group\n");
          plviol++;
          CNCwarn(log,wstr);
          pln->n--;
        }
        
        else{
          
          if(!is_planar_group(al,pln->natoms[n],pln->atoms[n],
                              pln->val[n], &plan))
          {
            sprintf(wstr,"Planarity violation: ");

/*             sprintf(wstr,"Planarity violation %d%s \n", */
/*                     al->resid[pln->atoms[n][0]],al->resname[pln->atoms[n][0]]); */

            plviol++;
            char ss[STRLEN];
            for(l=0;l<pln->natoms[n];l++){
              char ndum[STRLEN];
              char rdum[STRLEN];
              strcpy(ndum,al->name[pln->atoms[n][l]]);
              strcpy(rdum,al->resname[pln->atoms[n][l]]);
              trim(ndum);
              trim(rdum);
              if (l!=0) {
                sprintf(ss,"  \t\t\t\t%s-%d-%s\n",ndum,al->resid[pln->atoms[n][l]],rdum);
              }
              else {
                sprintf(ss,"%s-%d-%s\n",ndum,al->resid[pln->atoms[n][l]],rdum);
              }
              strcat(wstr,ss);
            }
           
/*             sprintf(wstr,"(%d%s %d%s ..)... %g -> %g\n",  */
/*                     al->id[pln->atoms[n][0]],al->name[pln->atoms[n][0]],  */
/*                     al->id[pln->atoms[n][1]],al->name[pln->atoms[n][1]],  */
/*                     pln->val[n],plan);  */
            
            CNCwarn(log,wstr); 
            if(!bIgn) {
              sprintf(wstr,"Increasing tolerance from %g to %g (lib = %g)\n",plan,plan*1.01,pln->val[n]);
              CNCwarn(log,wstr);
              pln->val[n] = plan*1.01;
            }
          }
        }
      }
    }
  }
  
  /* omega guys */
  
  for(k=0;k<ng->n;k++){
    if(strcmp(ng->resname[k],"OMEGA") == 0) omtol = ng->val[k];
  }



  int omega_n,omega_c;
  for(i=0;i<rl->nres-1;i++){
    if(IsProtein(al,rl->j0[i])){
      omega_c = -1;omega_n=-1;
      for(k=rl->j0[i];k<=rl->j1[i];k++){
        if(strcmp(al->name[k]," C  ") == 0){
          omega_c =  k;
        }
      }
      for(k=0;k<al->nbonds[omega_c];k++){
        if(strcmp(al->name[al->bonds[omega_c][k]]," N  ") == 0){
          omega_n =  al->bonds[omega_c][k];
        }
      }
      if(omega_n!=-1 && omega_c!=-1){
        pln->n++;
        pln = idx_realloc(pln,pln->n);
        n = pln->n-1;
        pln->val[n] = omtol;
        add_to_group(pln,n,omega_c);
        for(m=0;m<al->nbonds[omega_c];m++){
          add_to_group(pln,n,al->bonds[omega_c][m]);
        }
        for(m=0;m<al->nbonds[omega_n];m++){
          add_to_group(pln,n,al->bonds[omega_n][m]);
        }
        if(pln->natoms[n] < 4){
          sprintf(wstr,"Group has less than 4 atoms -> cannot create planar group\n");
          CNCwarn(log,wstr);
          plviol++;
          pln->n--;
        }
        else{
          
          if(!is_planar_group(al,pln->natoms[n],pln->atoms[n],omtol,&plan))
          {
            sprintf(wstr,"Planarity violation: ");
            
            char ss[STRLEN];
            for(l=0;l<pln->natoms[n];l++){
              char ndum[STRLEN];
              char rdum[STRLEN];
              strcpy(ndum,al->name[pln->atoms[n][l]]);
              strcpy(rdum,al->resname[pln->atoms[n][l]]);
              trim(ndum);
              trim(rdum);
              if (l!=0) {
                sprintf(ss,"  \t\t\t\t%s-%d-%s\n",ndum,al->resid[pln->atoms[n][l]],rdum);
              }
              else {
                sprintf(ss,"%s-%d-%s\n",ndum,al->resid[pln->atoms[n][l]],rdum);
              }
              
              strcat(wstr,ss);
            }

           CNCwarn(log,wstr);
           plviol++; 
            if(!bIgn) {
              sprintf(wstr,"Increasing tolerance from %g to %g (lib = %g)\n",plan,plan*1.01,pln->val[n]);
              CNCwarn(log,wstr);
              pln->val[n] = plan*1.01;
            }
          }
        }
        

      }
    }
  }
  return plviol;
  
}

