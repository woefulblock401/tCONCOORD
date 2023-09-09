#include <tconcoord.h>

t_rotation *rotation_init(void)
{
  t_rotation *rot = NULL;
  snew(rot,1);
  rot->n = 0;
  strcpy(rot->resname,"");
  snew(rot->names,1);
  return rot;
}

t_rotation *rotation_realloc(t_rotation *rot, int n)
{
  srenew(rot->names,n);
  return rot;
}

t_rotdata *rotdata_init(void)
{
  t_rotdata *rt = NULL;
  snew(rt,1);
  snew(rt->resname,1);
  rt->n = 0;
  snew(rt->rots,1);
  return rt;
}
  
t_rotdata *rotdata_realloc(t_rotdata *rt, int n)
{
  srenew(rt->rots,n);
  srenew(rt->resname,n);
  return rt;
}



t_rotdata *read_rotation_lib(void)
{
  FILE *fp = cnclib("rotations.dat");
  char line[STRLEN];
  int i,k, nrot;
  nrot = 0;
  
  char resn[STRLEN];
  t_rotdata *rt = rotdata_init();
#define MAXPTR 254
  char *ptr[MAXPTR];
  int nel;
  
  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      sscanf(line+1,"%s",resn);
    }
    else{
      rt->n++;
      nrot = rt->n-1;
      rt = rotdata_realloc(rt,nrot+1);
      rt->rots[nrot] = rotation_init();
      nel = str_nelem(line,MAXPTR,ptr);
      rt->rots[nrot]->n = nel;
      rt->rots[nrot] = rotation_realloc(rt->rots[nrot],nel);
      strcpy(rt->resname[nrot],resn);
      for(i=0;i<nel;i++){
        sscanf(ptr[i],"%s",rt->rots[nrot]->names[i]);
      }
    }
  }
   for(i=0;i<nrot+1;i++){ 
     for(k=0;k<rt->rots[i]->n;k++){ 
       extend_name(rt->rots[i]->names[k]); 
     } 
   } 
   return rt;
}

bool is_cys_with_ss_bond(t_atomlist *al, t_resl *rl, int j)
{
  int i,k;
  for(i=rl->j0[j];i<=rl->j1[j];i++){
    if(strcmp(al->name[i]," SG ") == 0){
      for(k=0;k<al->nbonds[i];k++){
        if(strcmp(al->name[i]," SG ") == 0) return TRUE;
      }
    }
  }
  return FALSE;
}

bool posres_on_sidechain(t_atomlist *al, t_resl *rl, int j)
{
  int i,k;
  for(i=rl->j0[j];i<=rl->j1[j];i++){
    if(al->order[i] > 1 && al->isposres[i]) return TRUE;
  }
  return FALSE;
}
  

void write_rotations(FILE *fp,t_atomlist *al, t_resl *rl, t_rotdata *rt)
{
  int i,k,j, l;
  char resn[STRLEN];
  char name1[STRLEN];
  char name2[STRLEN];
  char name3[STRLEN];
  int id1, id2, id3;
  int *arr;
  longstr *namearr;
  int size = 0;
  
  for(i=0;i<rl->nres;i++){
    if(IsProtein(al,rl->j0[i]) &&
       strcmp(rl->resname[i],"GLY") != 0 &&
       strcmp(rl->resname[i],"PRO") != 0 && 
       !is_cys_with_ss_bond(al,rl,i) && !posres_on_sidechain(al,rl,i)
       && !rl->sc_restr[i])
    {
      /* go through rotations */
      strcpy(resn,rl->resname[i]);
      
      for(k=0;k<rt->n;k++){
        if(strcmp(rt->resname[k],resn) == 0){
          strcpy(name1, rt->rots[k]->names[0]); /* first atom in bond */
          strcpy(name2, rt->rots[k]->names[1]); /* second atom in bond */
          id1 = id2 = -1;
          for(j=rl->j0[i];j<=rl->j1[i];j++){
            if(strcmp(al->name[j],name1) == 0) id1 = j+1;
            else if(strcmp(al->name[j],name2) == 0) id2 = j+1;
          }
          if(id1 == -1){
            printf("could not find atom \"%s\" %s\n",name1, resn);
          }
          if(id2 == -1){
            printf("could not find atom \"%s\" %s\n",name2, resn);
          }
          snew(arr,2);
          snew(namearr,2);
          arr[0] = id1;
          strcpy(namearr[0],al->name[id1-1]);
          arr[1] = id2;
          strcpy(namearr[1],al->name[id2-1]);
          size=2;
/*           printf("Searching for res %s/%d bond %s-%s\n",resn,i+1,name1,name2); */
          
          for(l=2;l<rt->rots[k]->n;l++){
            bool bFound = FALSE;
            for(j=rl->j0[i];j<=rl->j1[i];j++){
              if(strcmp(rt->rots[k]->names[l],al->name[j]) == 0){
                bFound = TRUE;
              
                size++;
                srenew(arr,size);
                srenew(namearr,size);
                arr[size-1] = j+1;
                strcpy(namearr[size-1],al->name[j]);
              }
            }
/*             if(!bFound){ */
/*               printf("not found \"%s\"\n",rt->rots[k]->names[l] ); */
/*               printf("resn = %s\n",rl->resname[i]); */
              
/*               for(j=rl->j0[i];j<=rl->j1[i];j++){ */
/*                 printf("\t %s\n",al->name[j]); */
/*               } */
              
/*             } */
            
          }
          if(size > 2)
          {
            for(l=0;l<size;l++){
              fprintf(fp," %d ",arr[l]);
            }
/*             fprintf(fp,"; "); */
/*             for(l=0;l<size;l++){ */
/*               fprintf(fp," %s ",namearr[l]); */
/*             } */
            fprintf(fp,"\n");
          }
          sfree(arr);
          sfree(namearr);
        }
      }
    }
  }
}

t_idxgroups *read_rotations(char *filename, bool bVerbose)
{
  int i,k;
  int nel, id;
  char *ptr[MAXPTR];
  if(bVerbose){
    fprintf(stderr,"Reading rotation file %s\n",filename);
  }
  
  FILE *fp = ffopen(filename,"r");
  char line[STRLEN];
  t_idxgroups *rot = idx_init();
#define MAXPTR 254        
  while(get_a_line(fp,line,STRLEN)){
    nel = str_nelem(line,MAXPTR,ptr);
    rot->n++;
    rot = idx_realloc(rot,rot->n);
    for(i=0;i<nel;i++){
      sscanf(ptr[i],"%d",&id);
      id-=1;
      add_to_group_nosort(rot,rot->n-1,id);
    }
  }
  if(bVerbose){
    fprintf(stderr,"Initialized %d rotation index groups\n",rot->n);
  }
  
  return rot;
}

void make_bb_rotations(FILE *fp, t_atomlist *al, t_resl *rl, int nrot, bool bVerbose)
{
  
  int i,k,m;
  int this_res = 0;
  bool bOk = FALSE;
  int *storage = NULL;
  int n = 0;
  int idx;
  t_idxgroups *rot = idx_init();
  
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->name[i]," CA ") == 0 && !al->isposres[i])
    {
      n = 0;
      if(storage != NULL) sfree(storage);
      this_res = al->resid[i]-1;
      if(rl->psi_restr[this_res] || rl->phi_restr[this_res+nrot+1]) {
        bOk = FALSE;
      }
      else {
        bOk = TRUE;
      }
      snew(storage,1);
      storage[n] = i;
      n++;
      idx = get_atom_idx(rl,this_res+nrot+1,al," CA ");
      if (idx != -1 && strcmp(al->resname[idx],"PRO") != 0) {
        srenew(storage,n+1);
        storage[n] = idx;
        n++;
      }
      else {
        bOk = FALSE;
      }
      for(k=1;k<=nrot;k++) {
        
        if(is_cys_with_ss_bond(al,rl,this_res+k) ||
           rl->phi_restr[this_res+k] || rl->psi_restr[this_res+k])
        {
          bOk = FALSE;
        }
      }
      
      
      if(bOk) 
      {
        /* all atom from i except N and H  */

        for(k=rl->j0[this_res];k<=rl->j1[this_res];k++){
          if(strcmp(al->name[k]," N  ") != 0 &&
             strcmp(al->name[k]," CA ") != 0 &&
             strcmp(al->name[k]," H  ") != 0)
          {
            if(al->isposres[k]) {
              bOk = FALSE;
            }
            else {
              srenew(storage,n+1);
              storage[n] = k;
              n++;
            }
            
          }
        }

        /* all atoms from res i+1 */
        if(bOk) {

          m = 1;
          do {
            
            for(k=rl->j0[this_res+m];k<=rl->j1[this_res+m];k++){
              if(al->isposres[k]) {
                bOk = FALSE;
              }
              else {
                srenew(storage,n+1);
                storage[n] = k;
                n++;
              }
            }
            m++;
          }while (m<=nrot);
          
        }
        if(bOk) {
        
          for(k=rl->j0[this_res+nrot+1];k<=rl->j1[this_res+nrot+1];k++){
            if(strcmp(al->name[k]," C  ") != 0 &&
               strcmp(al->name[k]," CA ") != 0 &&
               strcmp(al->name[k]," O  ") != 0)
            {
              if(al->isposres[k]) {
                bOk = FALSE;
              }
              else {
                srenew(storage,n+1);
                storage[n] = k;
                n++;
              }
            }
            
          }
        }
        if(bOk) {
          rot->n++;
          rot = idx_realloc(rot,rot->n);
          for(k=0;k<n;k++){
            add_to_group_nosort(rot,rot->n-1,storage[k]);
          }
        }
        
      }
    }
  }

  for(i=0;i<rot->n;i++){
    for(k=0;k<rot->natoms[i];k++){
      fprintf(fp," %d", rot->atoms[i][k]+1);
    }
    fprintf(fp,"\n");
  }
  free_groups(rot);
  
}

      

void rotate_group(t_atomlist *al, t_idxgroups *rot, int idx, real phi)
{
  /* the bond to rotate is idx->atoms[1] - idx->atoms[0] 
     the rest of the group is rotated around this vector 
  */

  int i;
  rvec diff;
  rvec scaled;
  matrix tm1, tm2;
  rvec_sub(al->x[rot->atoms[idx][1]],al->x[rot->atoms[idx][0]], diff);
/*   printf("Moving %d %d\n",rot->atoms[idx][1]+1, rot->atoms[idx][0]+1); */
  
  scaleVector(diff,scaled);
  createRotationMatrix1(scaled,tm1);
  createRotationMatrix2(scaled,tm2);
  for(i=2;i<rot->natoms[idx];i++){
    rotateVectorAroundVector(al->x[rot->atoms[idx][i]],al->x[rot->atoms[idx][1]]
                             ,tm1,tm2,phi);
  }
}

void do_random_rotations(t_atomlist *al, t_idxgroups *rot, int n, gmx_rng_t rng)
{
  int i;
  real phi;
  int idx;
  int nrot;
  
/*   int nrot = (int) rot->n * (real) n/100.; */
  if(n > rot->n) {
    nrot=rot->n;
  }
  else
  {
    nrot = n;
  }
  
/*   printf("Will perturb %d of %d rotations\n",nrot,rot->n); */
  
  for(i = 0;i<nrot;i++){

    idx = random_int_in_range(rng, 0, rot->n-1);
/*     printf("rotating group %d \n",idx); */
    
    phi = gmx_rng_uniform_real(rng)*2*PI;
    rotate_group(al, rot, idx, phi);
  }
}


