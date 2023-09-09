#include <tconcoord.h>

/*=============================================================*/
t_idxgroups *idx_init(void)
/* initialise group structure */
{
  t_idxgroups *grp = NULL;
  snew(grp,1);
  grp->n = 0;
/*   grp->natoms = NULL; */
/*   grp->atoms = NULL; */
/*   grp->val = NULL; */
/*   grp->flag = NULL; */
  snew(grp->natoms,1);
  snew(grp->atoms,1);
  snew(grp->val,1);
  snew(grp->flag,1);
  return grp;
}
/*=============================================================*/
t_idxgroups *idx_realloc(t_idxgroups *grp, int n)
/* reallocate memory for n groups */
{
  srenew(grp->natoms,n);
  srenew(grp->atoms,n);
  srenew(grp->val,n);
  srenew(grp->flag,n);
  grp->natoms[n-1] = 0;
  grp->val[n-1] = 0.;
  grp->flag[n-1] = 0;
/*   grp->atoms[n-1] = NULL; */
  snew(grp->atoms[n-1],1);
  
  return grp;
}
/*=============================================================*/
void add_to_group(t_idxgroups *grp, int ir, int id)
/* add atom id to group ir
   if id is not present in group yet.
   Then reorder group*/
{
  int i;
  bool check = FALSE;
  for(i=0;i<grp->natoms[ir];i++){
    if(grp->atoms[ir][i] == id) check = TRUE;
  }
  if(!check){
    grp->natoms[ir]+=1;
     srenew(grp->atoms[ir],grp->natoms[ir]);
    grp->atoms[ir][grp->natoms[ir]-1] = id;
    i_sort(grp->atoms[ir],grp->natoms[ir]); 
  }
}

/*=============================================================*/
void add_to_group_nosort(t_idxgroups *grp, int ir, int id)
/* add atom id to group but don't reorder */

{
  int i;
  bool check = FALSE;
  for(i=0;i<grp->natoms[ir];i++){
    if(grp->atoms[ir][i] == id) check = TRUE;
  }
  if(!check){
    grp->natoms[ir]+=1;
    srenew(grp->atoms[ir],grp->natoms[ir]);
    grp->atoms[ir][grp->natoms[ir]-1] = id;
/*     i_sort(grp->atoms[ir],grp->natoms[ir]);  */
  }
}



/*=============================================================*/
void free_groups(t_idxgroups *grp)

/* free memory for idxgroups structure */

{
  int i;
  for(i=0;i<grp->n;i++){
    sfree(grp->atoms[i]);
  }
  sfree(grp->atoms);
  sfree(grp->natoms);
  sfree(grp->val);
  sfree(grp->flag);
  sfree(grp);
  grp = NULL;
  
}
  
/*=============================================================*/
t_idxgroups *del_group(t_idxgroups *grp, int i)

/* delete a group (needs speed up ! )
 */
{
  int k,j;
  for(k=i;k<grp->n-1;k++){
    if(grp->natoms[k+1] > grp->natoms[k]){
      srenew(grp->atoms[k],grp->natoms[k+1]);
    }
    grp->natoms[k] = grp->natoms[k+1];
    grp->val[k] = grp->val[k+1];
    grp->flag[k] = grp->flag[k+1];
    for(j=0;j<grp->natoms[k];j++){
      grp->atoms[k][j] = grp->atoms[k+1][j];
    }
  }
  sfree(grp->atoms[grp->n-1]);
  grp->n-=1;
  return grp;
}
/*=============================================================*/
void sort_group_by_order(t_atomlist *al, t_idxgroups *idx)

/* sort group entries by atom order */

{  
  int i,k;
  int dum;
  
  for(k=0;k<idx->n;k++){
    bool changed = TRUE;
    while(changed){
      changed = FALSE;
      for(i=0;i<idx->natoms[k]-1;i++){
        if(al->order[idx->atoms[k][i]] > al->order[idx->atoms[k][i+1]]){
          dum = idx->atoms[k][i];
          idx->atoms[k][i] = idx->atoms[k][i+1];
          idx->atoms[k][i+1] = dum;
          changed = TRUE;
          
          
        }
      }
    }
    idx->val[k] = DIHED(al,idx->atoms[k][0],idx->atoms[k][1],
                        idx->atoms[k][2],idx->atoms[k][3]);
  }
}
/*=============================================================*/
void sort_group_by_posres(t_atomlist *al, t_idxgroups *idx)

/* check if we have posres in this group */

{  
  int i,k;
  int dum;
  
  for(k=0;k<idx->n;k++){
    bool changed = TRUE;
    while(changed){
      changed = FALSE;
      for(i=0;i<idx->natoms[k]-1;i++){
        if(!al->isposres[idx->atoms[k][i]] && al->isposres[idx->atoms[k][i+1]]){
          dum = idx->atoms[k][i];
          idx->atoms[k][i] = idx->atoms[k][i+1];
          idx->atoms[k][i+1] = dum;
          changed = TRUE;
          
          
        }
      }
    }
    idx->val[k] = DIHED(al,idx->atoms[k][0],idx->atoms[k][1],
                        idx->atoms[k][2],idx->atoms[k][3]);
  }
}

/*=============================================================*/
bool is_same_group(t_idxgroups *grp, int i, int k)

/* check whether two group are identical */

{
  if(grp->natoms[i] != grp->natoms[k]) return FALSE;
  else{
    int j;
    for(j=0;j<grp->natoms[i];j++){
      if(grp->atoms[i][j]!=grp->atoms[k][j]) return FALSE;
    }
  }
  return TRUE;
}
/*=============================================================*/
bool is_redundant_group(t_idxgroups *grp, int i, int k)

/* check whether group k is a subgroup of group i. 
   That means, all atoms of group k have to be listed 
   in group i.
*/
{

  if(grp->natoms[k] > grp->natoms[i]) return FALSE;
  int j,l;
  
  for(j=0;j<grp->natoms[k];j++){
    bool check = FALSE;
    for(l=0;l<grp->natoms[i];l++){
      if(grp->atoms[k][j]==grp->atoms[i][l]) check = TRUE;
    }
    if(!check) {
      return FALSE;
    }
    
  }
  return TRUE;
}


/*=============================================================*/
void print_groups(t_atomlist *al,t_idxgroups *grp)

/* self explaining */

{
  int i,k;
  for(i=0;i<grp->n;i++){
    printf("grp-nr: %d\n[",i);
    for(k=0;k<grp->natoms[i];k++){
      printf("%d (%s/%d%s)[%g %g %g]",al->id[grp->atoms[i][k]],al->name[grp->atoms[i][k]],
             al->resid[grp->atoms[i][k]],al->resname[grp->atoms[i][k]],
             al->x[grp->atoms[i][k]][0],al->x[grp->atoms[i][k]][1],al->x[grp->atoms[i][k]][2]);
    }
    printf("]\n");
  }
}

/*=============================================================*/
void print_group(t_atomlist *al,int *idx, int n)

/* self explaining */

{
  int i;
  printf("group: %d atoms\n",n);
  for(i=0;i<n;i++){
    printf("%d (%s/%d%s)\n",al->id[idx[i]],al->name[idx[i]],
           al->resid[idx[i]],al->resname[idx[i]]);
  }
}

/*=============================================================*/
/* void print_simple_group(int *idx, int n) */


/* { */
/*   int i; */
/*   printf("group: %d id's\n",n); */
/*   for(i=0;i<n;i++){ */
/*     printf("%d ",idx[i]); */
/*     if (i % 10 == 0) printf("\n"); */
/*   } */
/*   printf("\n"); */
/* } */

/*=============================================================*/
void copy_group(t_idxgroups *grp1,int idx1,t_idxgroups *grp2,int idx2)

/* copy grp1 to grp2  */

{
  int i;
/*   sfree(grp2->atoms[idx2]);   */
  snew(grp2->atoms[idx2],grp1->natoms[idx1]);
  for(i=0;i<grp1->natoms[idx1];i++){
    grp2->atoms[idx2][i] = grp1->atoms[idx1][i];
  }
  grp2->natoms[idx2] = grp1->natoms[idx1];
  grp2->flag[idx2] = grp1->flag[idx1];
  grp2->val[idx2] = grp1->val[idx1];

}
/*=============================================================*/
/*=============================================================*/
t_idxgroups *remove_redundant_groups(t_idxgroups *idx)



{
  int count = 1;
  int nrounds = 0;
/*   time_t tstart; */
/*   time(&tstart); */
    t_idxgroups *idx2 = idx_init();
/*   do */
/*   { */
    bool remove[idx->n];
    nrounds+=1;
/*     fprintf(stderr,"\rremoving groups..... round %d\n",nrounds);  */
    
    int i,k;
    int ng0 = idx->n;
    for(i=0;i<idx->n;i++) remove[i] = FALSE;
    for(i=0;i<idx->n;i++){
      for(k=i+1;k<idx->n;k++){
        if(idx->natoms[i] > idx->natoms[k]){
          if(!remove[k] && !remove[i])
            remove[k] = is_redundant_group(idx,i,k);
        }
        else{
          if(!remove[i] && !remove[k])
            remove[i] = is_redundant_group(idx,k,i);
        }
      }
    }
    int r = 0;
    for(i=0;i<idx->n;i++){
      if(remove[i]) {
        idx->flag[i] = TRUE;
      }
      else{
        r++;
      }
    }
    

    idx2->n = r;
    idx2 = idx_realloc(idx2,r+1);
    count = 0;
    
    for(i=0;i<ng0;i++){
      if(!remove[i]){
        copy_group(idx,i,idx2,count);
        count+=1;
      }
    }
    

/*     count = 0; */
/*     for(i=0;i<ng0;i++){ */
/*       if(remove[i]) */
/*       { */
/*         del_group(idx,i-count); */
/*         count++; */
/*       } */
/*     } */
/*   } */
/*   while(count); */
/*   fprintf(stderr,"\n");  */
/*   time_t tend; */
/*   time(&tend); */
/*   fprintf(stderr,"Total computation time: %gs\n",difftime(tend,tstart)); */
/*   free_groups(idx);  */
  
  return idx2;
}

/*=============================================================*/

t_idxgroups *merge_groups(t_idxgroups *grp, int min, bool bVerbose)

/* fuse groups with a least 'min' members in common */

{
  int i,k,j,l;
  int count;
  bool merged = TRUE;
  int level = 0;
  while(merged){
    level+=1;
     grp=remove_redundant_groups(grp); 
    int n0 = grp->n;
    merged = FALSE;
    for(i=0;i<grp->n-1;i++){
      if(bVerbose){
        fprintf(stderr,"\rtCNC__log_> Reducing groups level ... %4d (%8d groups) %d",level,grp->n,i+2); 
        fflush(stderr);
      }
      for(k=i+1;k<grp->n;k++){
/*         if(grp->natoms[k] < 7) { */
          
          count = 0;
          for(j=0;j<grp->natoms[i];j++){
            for(l=0;l<grp->natoms[k];l++){
              if(grp->atoms[i][j] == grp->atoms[k][l]){
                count++;
              }
            }
          }
          if(count >=min){
            merged = TRUE;
            for(j=0;j<grp->natoms[k];j++){
              add_to_group(grp,i,grp->atoms[k][j]);
              grp->flag[k] = TRUE;
            }
          }
/*         } */
        
        }
      }
    if(bVerbose)
      fprintf(stderr,"\n");
  }
  return grp;
  
}
/*=============================================================*/
void add_to_group_blind(t_idxgroups *grp, int ir, int id)

/* add atom to group without check */

{
  int i;
  grp->natoms[ir]+=1;
  srenew(grp->atoms[ir],grp->natoms[ir]);
  grp->atoms[ir][grp->natoms[ir]-1] = id;
  i_sort(grp->atoms[ir],grp->natoms[ir]); 
}
/*=============================================================*/
int is_in_group(t_idxgroups *grp, int idx, int item)

/* check whether atom item is member of group idx */

{
  int i;
  int count = 0;
  for(i=0;i<grp->natoms[idx];i++){
    if(grp->atoms[idx][i] == item) count++;
  }
  return count;
}


/*=============================================================*/
void sort_groups(t_idxgroups *grp)

/* sort groups by size */

{
  int i;
  int idum;
  int *ptr;
  bool changed = TRUE;
  while(changed){
    changed = FALSE;
    for(i=0;i<grp->n-1;i++){
      if(grp->natoms[i]<grp->natoms[i+1]){
        idum = grp->natoms[i];
        snew(ptr,grp->natoms[i]);
        copy_iarr(ptr,grp->atoms[i],grp->natoms[i]);
        sfree(grp->atoms[i]);
        snew(grp->atoms[i],grp->natoms[i+1]);
        copy_iarr(grp->atoms[i],grp->atoms[i+1],grp->natoms[i+1]);
        grp->natoms[i] = grp->natoms[i+1];
        grp->natoms[i+1] = idum;
        sfree(grp->atoms[i+1]);
        snew(grp->atoms[i+1],idum);
        copy_iarr(grp->atoms[i+1],ptr,idum);
        changed = TRUE;
      }
    }
  }
}

/*=============================================================*/
void write_group_script(char *filename, char *sel, t_idxgroups *grp)

/* write PyMOL script for visualising group definition */

{
  FILE *fp = ffopen(filename,"w");
  int i,k,j;
  int beg,end;

  for(i=0;i<grp->n;i++){
    fprintf(fp,"select grp_%d, %s ",i,sel);
    if(grp->natoms[i] < 5){
      for(k=0;k<grp->natoms[i]-1;k++){
        fprintf(fp,"%d+",grp->atoms[i][k]+1);
      }
      fprintf(fp,"%d\n",grp->atoms[i][grp->natoms[i]-1]+1);
    }
    else{
      beg = grp->atoms[i][0];
      end = 0;
      for(k=0;k<grp->natoms[i]-1;k++){
        if(grp->atoms[i][k] != grp->atoms[i][k+1] -1)
        {
          end = grp->atoms[i][k];
          if(end - beg == 0){
            fprintf(fp,"%d+",grp->atoms[i][k]+1);
          }
          else{
            fprintf(fp,"%d-%d+",beg+1,end+1);
          }
          beg = grp->atoms[i][k+1];
        }
        else end = grp->atoms[i][k+1];
      }
      if(end - beg == 0){
        fprintf(fp,"%d+",grp->atoms[i][k]+1);
      }
      else{
        fprintf(fp,"%d-%d+",beg+1,end+1);
      }
      
      fprintf(fp,"%d\n",grp->atoms[i][grp->natoms[i]-1]+1);
 /*      for(k=0;k<grp->natoms[i];k++){ */
/*         fprintf(fp,"%d ",grp->atoms[i][k]+1); */
/*       } */
/*       fprintf(fp,"\n"); */

    }
  }
  fprintf(fp,"deselect\n");
  fprintf(fp,"hide\n");
  fprintf(fp,"set sphere_scale, 0.3\n");
  fprintf(fp,"color grey\n");
  fprintf(fp,"show cartoon\n"); 
  fclose(fp);
}


