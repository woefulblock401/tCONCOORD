#include <tconcoord.h>

t_idxgroups *analyze_core(t_idxgroups *phb, t_idxgroups *resgroups, t_resl *rl)
{
  int i,k,j,l;
  int res;
  int counter[resgroups->n];
  int *merge;
  snew(merge,1);
  int nmerge = 0;
  

  for(i=0;i<resgroups->n;i++){
    counter[i] = 0;
  }

  
  
  for(i=0;i<phb->n;i++){
    if(phb->natoms[i] > 3){
/*       printf("Cluster : %d\n",i); */
      nmerge = 0;
      for(j=0;j<phb->natoms[i];j++){
        res = phb->atoms[i][j];
        for(k=0;k<rl->ngrps[res];k++){
          counter[rl->grps[res][k]]+=1;
/*           printf("resi = %d grpnr = %d\n",res+1,rl->grps[res][k]); */
        }
      }
      for(j=0;j<resgroups->n;j++){
/*         printf("group[%d] = %d\n",j,counter[j]); */
        if(counter[j]  > 2) { 
          nmerge++;
          srenew(merge,nmerge);
          merge[nmerge-1] = j;
        }
        counter[j] = 0;
      }
      if(nmerge!=0){
        for(j=1;j<nmerge;j++){
          for(k=0;k<resgroups->natoms[merge[j]];k++){
            add_to_group(resgroups,merge[0],resgroups->atoms[merge[j]][k]);
          }
        }
        for(j=0;j<phb->natoms[i];j++){
          add_to_group(resgroups,merge[0],phb->atoms[i][j]);
        }
        
      }
      
    }
  }
  sfree(merge);
  return resgroups;
}



t_idxgroups *get_hydrophobic_patches(t_idxgroups *resgroups, 
                                     t_idxgroups *phb, t_contab *rct)
{
  int i,k,j,l,m,n,o,x;
  int membersA, membersB;
  int n_sc_hbonds;
  int res1,res2;
  
  bool merged = TRUE;
  
/* go through the list of resgroups and find connections */
  while (merged) {
    merged = FALSE;
    resgroups = remove_redundant_groups(resgroups);  
    for(i=0;i<resgroups->n-1;i++){
      for(k=0;k<phb->n;k++){
        membersA = 0;     

        
        for(j=0;j<resgroups->natoms[i];j++){
          for(l=0;l<phb->natoms[k];l++){
            if(phb->atoms[k][l] == resgroups->atoms[i][j]) {
              membersA++;
            }
          }
        }
        /* now we have the number af members of group i in patch j */
        /* now we search all other clusters */
        for(m=i+1;m<resgroups->n;m++){
          n_sc_hbonds = 0;
          membersB=0;
          for(j=0;j<resgroups->natoms[m];j++){
            res2 = resgroups->atoms[m][j];
            for(l=0;l<phb->natoms[k];l++){
              if(phb->atoms[k][l] == resgroups->atoms[m][j]) {
                membersB++;
              }
            }
            for(x=0;x<resgroups->natoms[i];x++){
              res1 = resgroups->atoms[i][x];
              for(o=0;o<rct->ncon[res1];o++){
                if(rct->con[res1][o] == res2){
                  if(rct->type[res1][o] == eSBH_BD ||
                     rct->type[res1][o] == eSBH_BA ||
                     rct->type[res1][o] == eSBH_SD ||
                     rct->type[res1][o] == eSBH_SA ||
                     rct->type[res1][o] == eSSHD ||
                     rct->type[res1][o] == eSSHA ||
                     rct->type[res1][o] == eMET) {
                    n_sc_hbonds++;
                  }
                }
                
              }
            }
            
          }
          /* now we have the number or shared residues */
/*           printf("group %d and group %d contribute with %d and %d residues\n",i,m,membersA,membersB); */
          if((membersA > 2 &&  membersB > 2) ||  n_sc_hbonds > 2 ){
            merged = TRUE;
            printf("Merging group %d and %d (%d)\n",i,m, n_sc_hbonds);
            for(j=0;j<resgroups->natoms[m];j++){
              add_to_group(resgroups,i,resgroups->atoms[m][j]);
            }
          }
          
        }
      }
    }
  }
  return resgroups;
}

t_idxgroups *sidechain_hbonds(FILE *log, t_idxgroups *resgroups, t_contab *rct)
{
  int i,j,k,l,m;
  int res1, res2;
  int n_sc_hbonds;
  int *merge;
  snew(merge,1);
  int nmerge = 0;
  bool merged = TRUE;
  char logstr[STRLEN];
  
  while (merged){
    merged = FALSE;
    resgroups = remove_redundant_groups(resgroups);  
    for(i=0;i<resgroups->n-1;i++){
      nmerge = 0;
      srenew(merge,1);
      for(k=i+1;k<resgroups->n;k++){
        n_sc_hbonds = 0;
        for(j=0;j<resgroups->natoms[i];j++){
          res1 = resgroups->atoms[i][j];
          for(l=0;l<resgroups->natoms[k];l++){
            res2 = resgroups->atoms[k][l];
            for(m=0;m<rct->ncon[res1];m++){
              if(rct->con[res1][m] == res2){
                if(rct->type[res1][m] == eSBH_BD ||
                   rct->type[res1][m] == eSBH_BA ||
                   rct->type[res1][m] == eSBH_SD ||
                   rct->type[res1][m] == eSBH_SA ||
                   rct->type[res1][m] == eSSHD ||
                   rct->type[res1][m] == eSSHA ||
                   rct->type[res1][m] == eMET) {
                  n_sc_hbonds++;
                }
              }
            }
          }
        }
        if(n_sc_hbonds > 2) {
          nmerge++;
          srenew(merge,nmerge);
          merge[nmerge-1] = k;
        }
      }
      if(nmerge!=0){
        merged = TRUE;
        for(j=0;j<nmerge;j++){
          sprintf(logstr,"tCNC__log_> Merging residue groups %d and %d\n",i,merge[j]);
          CNClog(log,logstr);
          sprintf(logstr,"tCNC__log_> due to more than 2 stable sidechain hbonds\n");
          CNClog(log,logstr);
          for(k=0;k<resgroups->natoms[merge[j]];k++){
            add_to_group(resgroups,i,resgroups->atoms[merge[j]][k]);
          }
        }
        
      }
    }
  }
  
  return resgroups;
  
}


    



void calc_packing(t_atomlist *al)
{
  int i,k;
  real d;
  rvec diff;
  rvec scaled;
  rvec sumvec;
  real sum;
  
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->symbol[i],"H") != 0){ 
      clear_rvec(sumvec);
      for(k=0;k<al->nnb[i];k++){
        if(strcmp(al->symbol[al->nb[i][k]],"H") != 0){ 
          
          rvec_sub(al->x[i],al->x[al->nb[i][k]],diff);
/*       scaleVector(diff, scaled); */
          rvec_add(sumvec,diff,sumvec);
         } 
      }
      d = norm(sumvec)/(real) al->nnb[i];
      al->bfac[i] = d;
      } 
      else al->bfac[i] = 0.; 
    
  }
}

void calc_hbond_packing(t_atomlist *al)
{
  int i,k;
  real d;
  rvec diff;
  rvec scaled;
  rvec sumvec;
  real sum;
  int don;
  
  for(i=0;i<al->nhbonds;i++){
    don = al->hbonds[i]->don;
    clear_rvec(sumvec);
    for(k=0;k<al->nnb[don];k++){
      rvec_sub(al->x[don],al->x[al->nb[don][k]],diff);
      rvec_add(sumvec,diff,sumvec);
    }
    d = norm(sumvec)/(real) al->nnb[don];
    al->hbonds[i]->packing = d;
  } 
}

void calc_buried(t_atomlist *al)
{
  int i, k;
  
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->symbol[i],"H") != 0){
      
      int nex = 1+al->nbonds[i]+al->nb13[i]+al->nb14[i];
      int excluded[nex];
      excluded[0] = i;
      for(k=0;k<al->nbonds[i];k++){
        excluded[k+1] = al->bonds[i][k];
      }
      for(k=0;k<al->nb13[i];k++){
        excluded[al->nbonds[i]+k+1] = al->b13[i][k];
      }
      for(k=0;k<al->nb14[i];k++){
        excluded[al->nbonds[i]+al->nb13[i]+k+1] = al->b14[i][k];
      }
      al->bfac[i] = ray_search(al, al->x[i],excluded,nex); 
    }
  }
}

void packing_per_residue(t_atomlist *al, t_resl *rl, real *respack)
{
  int i,k;
  real sum;
  int nat;
  
  FILE *fp = ffopen("buried.dat","w");
  
  for(i=0;i<rl->nres;i++){
    sum = 0;
    nat = 0;
    
    for(k=rl->j0[i];k<=rl->j1[i];k++){
      if(strcmp(al->symbol[k],"H") != 0) {
        sum+=al->bfac[k];
        nat++;
      }
    }
    respack[i] = sum/(real) nat;
    fprintf(fp,"%d %s %g\n",i+1,rl->resname[i],respack[i]);
  }
}



static void copy2block(int n,atom_id *index,t_block *block)
{
  int i,n0;

  block->nr++;
  n0=block->nra;
  block->nra=n0+n;
  srenew(block->index,block->nr+1);
  block->index[block->nr]=n0+n;
  srenew(block->a,n0+n);
  for(i=0; (i<n); i++)
    block->a[n0+i]=index[i];
}


void print_simple_group(int *idx, int n)


{
  int i;
  printf("XXgroup: %d id's\n",n); 
  for(i=0;i<n;i++){
    printf("XX %6d \n",idx[i]+1);
/*     if (i && i % 10 == 0) printf("\n"); */
  }
  printf("//\n");
}
/*=============================================================*/
void read_additional_connections(char *filename, t_atomlist *al, 
                                 t_contab *ct)
{
  char line[STRLEN];
  char dum1[STRLEN];
  char btype[STRLEN];
  int at1,at2;
  
  FILE *fp = ffopen(filename,"r");
  while(get_a_line(fp,line,STRLEN)){
/*     printf("line  = %s\n",line); */
    
    sscanf(line,"%s %d %d",dum1,&at1,&at2);
    if(strcmp(dum1,"BOND") == 0){
      strcpy(btype,dum1);
      
      printf("ADDING BOND = %d-%d (%s)\n",at1,at2,btype);
      at1-=1;
      at2-=1;
      add_bond(al,at1,at2,ct,BOND,BOND); 
    }
  }
}
/*=============================================================*/

bool id_in_index(int id, int size, int *array)
{
  int i;
  for(i=0;i<size;i++){
    if(array[i] == id) return TRUE;
  }
  return FALSE;
}


/*=============================================================*/
bool is_omega_group(t_atomlist *al, int *idx, int n)
{
  int i;
  bool have_ca = FALSE;
  bool have_o = FALSE;
  bool have_n = FALSE;
  bool have_c = FALSE;
  
  for(i=0;i<n;i++){
    if(strcmp(al->name[idx[i]]," CA ") == 0) have_ca = TRUE;
    else if(strcmp(al->name[idx[i]]," N  ") == 0) have_n = TRUE;
    else if(strcmp(al->name[idx[i]]," O  ") == 0) have_o = TRUE;
    else if(strcmp(al->name[idx[i]]," C  ") == 0) have_c = TRUE;
  }
  if(have_ca && have_o && have_n && have_c) return TRUE;
  else return FALSE;
}
/*=============================================================*/
bool is_arg_head_group(t_atomlist *al, int *idx, int n)
{
  int i;
  bool have_cd = FALSE;
  bool have_ne = FALSE;
  bool have_cz = FALSE;
  bool have_nh1 = FALSE;
  bool have_nh2 = FALSE;
  
  for(i=0;i<n;i++){
    if(strcmp(al->name[idx[i]]," CD ") == 0) have_cd = TRUE;
    else if(strcmp(al->name[idx[i]]," NE ") == 0) have_ne = TRUE;
    else if(strcmp(al->name[idx[i]]," CZ ") == 0) have_cz = TRUE;
    else if(strcmp(al->name[idx[i]]," NH1") == 0) have_nh1 = TRUE;
    else if(strcmp(al->name[idx[i]]," NH2") == 0) have_nh2 = TRUE;
  }
  if(have_cd && have_ne && have_cz && have_nh1 && have_nh2) return TRUE;
  else return FALSE;
}

/*=============================================================*/

void posr2al(t_atomlist *al, atom_id *pos_id, int npos)

/* tell atoms in atomlist whether they are position restraint */

{
  int i;
  for(i=0;i<al->natoms;i++){
    al->isposres[i] = FALSE;
  }
  for(i=0;i<npos;i++){
    al->isposres[pos_id[i]] = TRUE;
  }
}
/*=============================================================*/
void flex2al(t_atomlist *al, atom_id *pos_id, int npos)

/* tell atoms in atomlist whether they are position restraint */

{
  int i;
  for(i=0;i<al->natoms;i++){
    al->isflex[i] = FALSE;
  }
  for(i=0;i<npos;i++){
    al->isflex[pos_id[i]] = TRUE;
  }
}

/*=============================================================*/

void get_restricted(t_atomlist *al, t_contab *rct, 
                    t_resl *rl, t_dihed *dihed)

/* get backbone dihedral restrictions */

{
  int i,j,k;
  bool phi_restr;
  bool psi_restr;
  bool sc_restr;
  bool bbdon;
  bool bbacc;
  bool sc;

  for(i=0;i<rct->n;i++){
    bbdon = FALSE;
    bbacc = FALSE;
    sc = FALSE;
    for(k=0;k<rct->ncon[i];k++){
      switch(rct->type[i][k]){
        case(eBBHD):
        case(eSBH_BD):
          bbdon = TRUE;
          break;
        case(eBBHA):
        case(eSBH_BA):
          bbacc = TRUE;
          break;
        case(eSSHD):
        case(eSSHA):
        case(ePHO):
          sc = TRUE;
          break;
      }
    }
    if(bbdon && bbacc){
      rl->phi_restr[i] = TRUE;
      rl->psi_restr[i] = TRUE;
    }
    if(sc){
      rl->sc_restr[i] = TRUE;
    }
    
  }
  for(i=0;i<dihed->n;i++){
    if(IS_PSI(al,dihed->at2[i],dihed->at3[i])){
      if(rl->psi_restr[al->resid[dihed->at2[i]]-1]){
        dihed->flag[i] = TRUE;
      }
    }
    else if(IS_PHI(al,dihed->at2[i],dihed->at3[i])){
      if(rl->phi_restr[al->resid[dihed->at2[i]]-1]){
        dihed->flag[i] = TRUE;
      }
    }
    else {
      if(rl->sc_restr[al->resid[dihed->at2[i]]-1]){
        dihed->flag[i] = TRUE;
      }
    }
    if( (al->isflex[dihed->at2[i]] || al->isflex[dihed->at3[i]]) && 
        (!IS_OMEGA(al,dihed->at2[i],dihed->at3[i]))) {
      dihed->flag[i] = FALSE;
    }
    if(IS_OMEGA(al,dihed->at2[i],dihed->at3[i])) dihed->flag[i] = TRUE;
    
    
  }
  

}
/*=============================================================*/
void print_connections(t_atomlist *al, t_contab *rct, t_resl *rl)
{
  int i,k;
  char dum[STRLEN];
  int resid;
  
  for(i=0;i<rct->n;i++){
     printf("Connections of %d phi = %d psi = %d\n",i+1, 
            rl->phi_restr[i],rl->psi_restr[i]); 
    for(k=0;k<rct->ncon[i];k++){
      resid = rct->con[i][k];
      switch(rct->type[i][k]){
        case eNPR:
          strcpy(dum,"NPR");
          break;
        case eCOV:
          strcpy(dum,"COV");
          break;
        case eSUL:
          strcpy(dum,"SUL");
          break;
        case eBBHD:
          strcpy(dum,"BBHD");
          break;
        case eBBHA:
          strcpy(dum,"BBHA");
          break;
        case eSBH_BD:
          strcpy(dum,"SBH_BD");
          break;
        case eSBH_BA:
          strcpy(dum,"SBH_BA");
          break;
        case eSBH_SD:
          strcpy(dum,"SBH_SD");
          break;
        case eSBH_SA:
          strcpy(dum,"SBH_SA");
          break;
        case eSSHA:
          strcpy(dum,"SSHD");
          break;
        case eSSHD:
          strcpy(dum,"SSHA");
          break;
        case ePHO:
          strcpy(dum,"PHO");
          break;
      }
      
      printf("res = %d type = %s\n",resid+1,dum); 
    }
  }
}


  


/*=============================================================*/
void do_resnet(t_contab *rct, t_idxgroups *resgroup)

/* build residue network */

{
  int i,j,k,l,m,n,o;
  int res1,res2,res3,res4, res5, res6;
  int count = 0;
  
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50);
    fflush(stderr);
    /* 4 - Rings */

   /* take a look into the connections of i */

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] != eNPR && 
         rct->type[i][k] != eEX &&
         rct->type[i][k] != eFO){
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] != eNPR && 
             rct->type[i][j] != eEX &&
             rct->type[i][j] != eFO){
            res2 = rct->con[i][j];

            /* now we have two connections of i. if they
               have a connection which is not i, the residues
               i, res1, res2 and the fourth one build a 4-ring
            */

            if(res1 > res2){
              for(l=0;l<rct->ncon[res1];l++){
                for(m=0;m<rct->ncon[res2];m++){
                  if(rct->con[res1][l] == rct->con[res2][m] &&
                     rct->con[res1][l] != i && rct->type[res1][l] != eNPR
                     && rct->type[res2][m] != eEX){

                    /* build the group */

                    resgroup->n+=1;
                    resgroup = idx_realloc(resgroup,resgroup->n);
                    add_to_group(resgroup,resgroup->n-1,i);
                    add_to_group(resgroup,resgroup->n-1,res1);
                    add_to_group(resgroup,resgroup->n-1,res2);
                    add_to_group(resgroup,resgroup->n-1,rct->con[res1][l]);
                    count++;
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  count= 0;
  
  /* 5 - Rings */
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50+50);
    fflush(stderr);
    
    /* browse the connections of i */ 

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] != eNPR && 
         rct->type[i][k] != eEX &&
         rct->type[i][k] != eFO){
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] != eNPR && 
             rct->type[i][j] != eEX &&
             rct->type[i][j] != eFO){
            res2 = rct->con[i][j];
            if(res1 != res2){

              /* now we have residue i with two neighbors res1 and res2. 
               let's find a neighbor of res1 who is not i and not an neighbor
              of res2.*/

              for(l=0;l<rct->ncon[res1];l++){
                if(rct->con[res1][l] != i && rct->type[res1][l] != eNPR && 
                   rct->type[res1][l] != eEX && rct->type[res1][l] != eFO &&
                   res2 != rct->con[res1][l]){
                  res3 = rct->con[res1][l];

                  /* got it. now find a neighbor of res2 which is not i and not a neighbor of res1 */

                  for(m=0;m<rct->ncon[res2];m++){
                    if(rct->con[res2][m] != i && rct->type[res2][m] != eNPR && rct->type[res2][m] != eEX &&
                       res1 != rct->con[res2][m]){
                      res4 = rct->con[res2][m];

                      /* ok. if now res3 and res4 are connected, we have found a 5-ring */

                      for(n=0;n<rct->ncon[res3];n++){
                        if(rct->con[res3][n] == res4 && rct->type[res3][n] !=eNPR && rct->type[res3][n] !=eEX){

                          /* store this group */

                          resgroup->n+=1;
                          resgroup = idx_realloc(resgroup,resgroup->n);
                          add_to_group(resgroup,resgroup->n-1,i);
                          add_to_group(resgroup,resgroup->n-1,res1);
                          add_to_group(resgroup,resgroup->n-1,res2);
                          add_to_group(resgroup,resgroup->n-1,res3);
                          add_to_group(resgroup,resgroup->n-1,res4);
                          count++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  count= 0;
  
  /* 6 - Rings */
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50+50);
    fflush(stderr);
    
    /* browse the connections of i */ 

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] != eNPR && 
         rct->type[i][k] != eEX &&
         rct->type[i][k] != eFO){
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] != eNPR && 
             rct->type[i][j] != eEX &&
             rct->type[i][j] != eFO){
            res2 = rct->con[i][j];
            if(res1 != res2){

              /* now we have residue i with two neighbors res1 and res2. 
               let's find a neighbor of res1 who is not i and not an neighbor
              of res2.*/

              for(l=0;l<rct->ncon[res1];l++){
                if(rct->con[res1][l] != i && rct->type[res1][l] != eNPR && 
                   rct->type[res1][l] != eEX && rct->type[res1][l] != eFO &&
                   res2 != rct->con[res1][l]){
                  res3 = rct->con[res1][l];

                  /* got it. now find a neighbor of res2 which is not i and not a neighbor of res1 */

                  for(m=0;m<rct->ncon[res2];m++){
                    if(rct->con[res2][m] != i && rct->type[res2][m] != eNPR && rct->type[res2][m] != eEX &&
                       res1 != rct->con[res2][m]){
                      res4 = rct->con[res2][m];

                      /* ok. now we check all neihbors of res3 and look if they are also neighbors of res 4 */
                      for(n=0;n<rct->ncon[res3];n++){
                        if(rct->type[res3][n] !=eNPR && rct->type[res3][n] !=eEX ){
                          res5 = rct->con[res3][n];
                          for(o=0;o<rct->ncon[res4];o++){
                            if(rct->type[res4][o] !=eNPR && rct->type[res4][o] !=eEX ){
                              res6 = rct->con[res4][o];
                              if(res5 == res6) {
                                /* store this group */
                                resgroup->n+=1;
                                resgroup = idx_realloc(resgroup,resgroup->n);
                                add_to_group(resgroup,resgroup->n-1,i);
                                add_to_group(resgroup,resgroup->n-1,res1);
                                add_to_group(resgroup,resgroup->n-1,res2);
                                add_to_group(resgroup,resgroup->n-1,res3);
                                add_to_group(resgroup,resgroup->n-1,res4);
                                add_to_group(resgroup,resgroup->n-1,res5);
/*                                 printf("found group %d %d %d %d %d %d\n",i+1,res1+1,res2+1,res3+1,res4+1,res5+1); */
                                
                                count++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
/*   fprintf(stderr,"Found %d 5-Rings\n",count);  */
  fprintf(stderr,"\n");
}

/************************************************************/
void do_hbond_resnet(t_contab *rct, t_idxgroups *resgroup)

/* build residue network */

{
  int i,j,k,l,m,n,o;
  int res1,res2,res3,res4, res5, res6;
  int count = 0;
  
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50);
    fflush(stderr);
    /* 4 - Rings */

   /* take a look into the connections of i */

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] == eCOV || rct->type[i][k] == eSUL ||
         rct->type[i][k] == eBBHD || rct->type[i][k] == eBBHA) {
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] == eCOV || rct->type[i][j] == eSUL ||
             rct->type[i][j] == eBBHD || rct->type[i][j] == eBBHA) {
            res2 = rct->con[i][j];

            /* now we have two connections of i. if they
               have a connection which is not i, the residues
               i, res1, res2 and the fourth one build a 4-ring
            */

            if(res1 > res2){
              for(l=0;l<rct->ncon[res1];l++){
                for(m=0;m<rct->ncon[res2];m++){
                  if(rct->con[res1][l] == rct->con[res2][m] &&
                     rct->con[res1][l] != i) {
                    if( (rct->type[res1][l] == eCOV ||
                         rct->type[res1][l] == eSUL ||
                         rct->type[res1][l] == eBBHA ||
                         rct->type[res1][l] == eBBHD) &&
                        (rct->type[res2][m] == eCOV ||
                         rct->type[res2][m] == eSUL ||
                         rct->type[res2][m] == eBBHA ||
                         rct->type[res2][m] == eBBHD)) {

                    /* build the group */

                      resgroup->n+=1;
                      resgroup = idx_realloc(resgroup,resgroup->n);
                      add_to_group(resgroup,resgroup->n-1,i);
                      add_to_group(resgroup,resgroup->n-1,res1);
                      add_to_group(resgroup,resgroup->n-1,res2);
                      add_to_group(resgroup,resgroup->n-1,rct->con[res1][l]);
                      count++;
                    }
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  count= 0;
  
  /* 5 - Rings */
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50+50);
    fflush(stderr);
    
    /* browse the connections of i */ 

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] == eCOV || rct->type[i][k] == eSUL ||
         rct->type[i][k] == eBBHD || rct->type[i][k] == eBBHA) {
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] == eCOV || rct->type[i][j] == eSUL ||
             rct->type[i][j] == eBBHD || rct->type[i][j] == eBBHA) {
            res2 = rct->con[i][j];
            if(res1 != res2){

              /* now we have residue i with two neighbors res1 and res2. 
               let's find a neighbor of res1 who is not i and not an neighbor
              of res2.*/

              for(l=0;l<rct->ncon[res1];l++){
                if(rct->type[res1][l] == eCOV || rct->type[res1][l] == eSUL ||
                   rct->type[res1][l] == eBBHD || rct->type[res1][l] == eBBHA) {
                  res3 = rct->con[res1][l];

                  /* got it. now find a neighbor of res2 which is not i and not a neighbor of res1 */
                  for(m=0;m<rct->ncon[res2];m++){
                    if(rct->type[res2][m] == eCOV || rct->type[res2][m] == eSUL ||
                       rct->type[res2][m] == eBBHD || rct->type[res2][m] == eBBHA) {
                      res4 = rct->con[res2][m];

                      /* ok. if now res3 and res4 are connected, we have found a 5-ring */

                      for(n=0;n<rct->ncon[res3];n++){
                        if(rct->con[res3][n] == res4 ){
                          if(rct->type[res3][n] == eCOV ||
                             rct->type[res3][n] == eSUL ||
                             rct->type[res3][n] == eBBHD ||
                             rct->type[res3][n] == eBBHA) 
                          {

                          /* store this group */
                            
                            resgroup->n+=1;
                            resgroup = idx_realloc(resgroup,resgroup->n);
                            add_to_group(resgroup,resgroup->n-1,i);
                            add_to_group(resgroup,resgroup->n-1,res1);
                            add_to_group(resgroup,resgroup->n-1,res2);
                            add_to_group(resgroup,resgroup->n-1,res3);
                            add_to_group(resgroup,resgroup->n-1,res4);
                            count++;
                          }
                          
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  count= 0;
  
  /* 6 - Rings */
  for(i=0;i<rct->n;i++){
    fprintf(stderr,"\rtCNC__log_> Building Residue Network.....%3.0f%%",(real)(i/(real)rct->n)*50+50);
    fflush(stderr);
    
    /* browse the connections of i */ 

    for(k=0;k<rct->ncon[i]-1;k++){
      if(rct->type[i][k] == eCOV || rct->type[i][k] == eSUL ||
         rct->type[i][k] == eBBHD || rct->type[i][k] == eBBHA) {
        res1 = rct->con[i][k];
        for(j=k+1;j<rct->ncon[i];j++){
          if(rct->type[i][j] == eCOV || rct->type[i][j] == eSUL ||
             rct->type[i][j] == eBBHD || rct->type[i][j] == eBBHA) {
            res2 = rct->con[i][j];
            if(res1 != res2){
              /* now we have residue i with two neighbors res1 and res2. 
               let's find a neighbor of res1 who is not i and not an neighbor
              of res2.*/

              for(l=0;l<rct->ncon[res1];l++){
                if( (rct->type[res1][l] == eCOV || rct->type[res1][l] == eSUL ||
                     rct->type[res1][l] == eBBHD || rct->type[res1][l] == eBBHA ) &&
                    rct->con[res1][l] != i && rct->con[res1][l] != res2) {
                  res3 = rct->con[res1][l];

                  /* got it. now find a neighbor of res2 which is not i and not a neighbor of res1 */

                  for(m=0;m<rct->ncon[res2];m++){
                    if((rct->type[res2][m] == eCOV || rct->type[res2][m] == eSUL ||
                        rct->type[res2][m] == eBBHD || rct->type[res2][m] == eBBHA) &&
                       rct->con[res2][m] != i && rct->con[res2][m] != res1 && rct->con[res2][m] != res3) {
                      res4 = rct->con[res2][m];

                      /* ok. now we check all neihbors of res3 and look if they are also neighbors of res 4 */
                      for(n=0;n<rct->ncon[res3];n++){
                        if((rct->type[res3][n] == eCOV || rct->type[res3][n] == eSUL ||
                            rct->type[res3][n] == eBBHD || rct->type[res3][n] == eBBHA) && 
                           rct->con[res3][n] != i &&  
                           rct->con[res3][n] != res1 &&  
                           rct->con[res3][n] != res2 &&
                           rct->con[res3][n] != res4 ) {
                          res5 = rct->con[res3][n];
/*                           printf(" %d %d %d %d - %d\n",i,res1,res2,res4,res5); */
                          
                          for(o=0;o<rct->ncon[res4];o++){
                            if(rct->type[res4][o] == eCOV || rct->type[res4][o] == eSUL ||
                               rct->type[res4][o] == eBBHD || rct->type[res4][o] == eBBHA) {
                              res6 = rct->con[res4][o];
                              if(res5 == res6) {
                                /* store this group */
                                resgroup->n+=1;
                                resgroup = idx_realloc(resgroup,resgroup->n);
                                add_to_group(resgroup,resgroup->n-1,i);
                                add_to_group(resgroup,resgroup->n-1,res1);
                                add_to_group(resgroup,resgroup->n-1,res2);
                                add_to_group(resgroup,resgroup->n-1,res3);
                                add_to_group(resgroup,resgroup->n-1,res4);
                                add_to_group(resgroup,resgroup->n-1,res5);
/*                                  printf("found group %d %d %d %d %d %d\n",i+1,res1+1,res2+1,res3+1,res4+1,res5+1);  */
                                
                                count++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
/*   fprintf(stderr,"Found %d 5-Rings\n",count);  */
  fprintf(stderr,"\n");
}



/*=============================================================*/
/* void grp_speed_up(t_idxgroups *grp) */
/* { */
/*   int i,k,j; */
  
/*=============================================================*/

void find_lonely_res(t_atomlist *al, t_idxgroups *grp, t_resl *rl,
                     t_contab *rct)
{
  int i,k,j,l;
  for(i=0;i<rl->nres;i++){
    rl->ngrps[i] = 0;
    for(k=0;k<grp->n;k++){
      for(j=0;j<grp->natoms[k];j++){
        if(i == grp->atoms[k][j]){
          rl->ngrps[i]+=1;
          srenew(rl->grps[i],rl->ngrps[i]);
          rl->grps[i][rl->ngrps[i]-1] = k;
        }
      }
    }
  }


  int nrounds = 0;
  do{
    int res;
    int mem[rl->nres];
    int ingrp[rl->nres];
    int count = 0;
    nrounds++;
    for(i=0;i<rl->nres;i++){
      if(rl->ngrps[i] == 0){
        t_idxgroups *tmp = idx_init();
        tmp = idx_realloc(tmp,1);
        for(k=0;k<rct->ncon[i];k++){
          if(rct->type[i][k] != eNPR){
            res= rct->con[i][k];
            for(j=0;j<rl->ngrps[res];j++){
              add_to_group_blind(tmp,0,rl->grps[res][j]);
            }
          }
        }
        for(k=0;k<tmp->natoms[0];k++){
          if(is_in_group(tmp,0,tmp->atoms[0][k] > 1)){
            mem[count] = i;
            ingrp[count] = tmp->atoms[0][k];
            count++;
/*           add_to_group(grp,grp->atoms[0][k],i); */
            /*         rl->ngrps[i]+=1; */
/*           srenew(rl->grps[i],rl->ngrps[i]); */
/*           rl->grps[i][rl->ngrps[i]-1] = grp->atoms[0][k]; */
          }
        }
        free_groups(tmp);
      }
    }
    for(i=0;i<count;i++){
      add_to_group(grp,ingrp[i],mem[i]);
      bool check = FALSE;
      for(k=0;k<rl->ngrps[mem[i]];k++){
        if(ingrp[i] == rl->grps[i][k]) check = TRUE;
      }
      if(!check){
        rl->ngrps[mem[i]]+=1;
        srenew(rl->grps[mem[i]],rl->ngrps[mem[i]]);
        rl->grps[mem[i]][rl->ngrps[mem[i]]-1] = ingrp[i];
      }
    } 
  }while(nrounds < 2);
}
/*=============================================================*/

int get_vdw_params(void)

/* ask for vdw parameter set */

{
  fprintf(stderr,"\nselect a set of Van der Waals parameters:\n");
  fprintf(stderr,"1: CONCOORD all atom\n");
  fprintf(stderr,"2: CONCOORD polar hydrogens\n");
  fprintf(stderr,"3: CONCOORD no hydrogens\n");
  fprintf(stderr,"4: OPLS - all atom\n");
  fprintf(stderr,"5: Yamber2 all atom\n");
  fprintf(stderr,"6: Li et al.\n");
  fprintf(stderr," >> ");
  char xx[STRLEN];
  fgets(xx,STRLEN,stdin);
  int param = atoi(xx);
  return param;
  
}
/*=============================================================*/
int get_bond_lib(void)

/* ask for bond parameters */

{
  fprintf(stderr,"\nselect a set of bond/angle parameters:\n");
  fprintf(stderr,"1: CONCOORD parameters\n");
  fprintf(stderr,"2: Engh-Huber parameters\n");
  fprintf(stderr," >> ");
  char xx[STRLEN];
  fgets(xx,STRLEN,stdin);
  int param = atoi(xx);
  return param;
}
/*=============================================================*/
t_atomlist *remove_H(t_atomlist *al, int polFlag)
{
  int i,k;
  int nat = 0;
  int count = 0;
  int new[al->natoms];
  int idx[al->natoms];

  /* new contains the new atom indices 
     idx contains the atoms indices
     that are kept
  */
  
  for(i=0;i<al->natoms;i++){
    new[i] = -1;
    idx[i] = -1;
  }

  if(!polFlag){
    for(i=0;i<al->natoms;i++){
      if(strcmp(al->symbol[i],"H") != 0){
        nat++;
        idx[nat-1] = i;
        new[i] = count;
        count++;
      }
    }
    
  }
  else{
    for(i=0;i<al->natoms;i++){
      if(strcmp(al->symbol[i],"H") != 0 ||
         is_polar_H(al,i)){
        nat++;
        idx[nat-1] = i;
        new[i] = count;
        count++;
      }
    }
  }

  
  t_atomlist *al2 = atomlist_init();
  al2 = al_realloc(al2,nat);
  al2->natoms = nat;


  for(i=0;i<nat;i++){
    copy_atom(al,idx[i],al2,i);            
  }
  
  for(i=0;i<nat;i++){
    int n = 0;
    al2->nbonds[i] = 0;
    for(k=0;k<al->nbonds[idx[i]];k++){
      if(!polFlag){
        if(strcmp(al->symbol[al->bonds[idx[i]][k]],"H")!= 0){
          al2->bonds[i][n] = new[al->bonds[idx[i]][k]];
          al2->nbonds[i]+=1;
          n++;
        }
      }
      else {
        if(strcmp(al->symbol[al->bonds[idx[i]][k]],"H")!= 0 ||
           is_polar_H(al,al->bonds[idx[i]][k])){
          al2->bonds[i][n] = new[al->bonds[idx[i]][k]];
          al2->nbonds[i]+=1;
          n++;
        }
      }
    }
  }
  
  for(i=0;i<nat;i++){
    int n = 0;
    al2->nb13[i] = 0;
    for(k=0;k<al->nb13[idx[i]];k++){
      if(!polFlag){
        if(strcmp(al->symbol[al->b13[idx[i]][k]],"H")!= 0){
          al2->b13[i][n] = new[al->b13[idx[i]][k]];
          al2->nb13[i]+=1;
          n++;
        }
      }
      else {
        if(strcmp(al->symbol[al->b13[idx[i]][k]],"H")!= 0 ||
           is_polar_H(al,al->b13[idx[i]][k])){
          al2->b13[i][n] = new[al->b13[idx[i]][k]];
          al2->nb13[i]+=1;
          n++;
        }
      }
    }
  }

  for(i=0;i<nat;i++){
    int n = 0;
    al2->nb14[i] = 0;
    for(k=0;k<al->nb14[idx[i]];k++){
      if(!polFlag){
        if(strcmp(al->symbol[al->b14[idx[i]][k]],"H")!= 0){
          al2->b14[i][n] = new[al->b14[idx[i]][k]];
          al2->nb14[i]+=1;
          n++;
        }
      }
      else {
        if(strcmp(al->symbol[al->b14[idx[i]][k]],"H")!= 0 ||
           is_polar_H(al,al->b14[idx[i]][k])){
          al2->b14[i][n] = new[al->b14[idx[i]][k]];
          al2->nb14[i]+=1;
          n++;
        }
      }
    }
  }

  for(i=0;i<nat;i++){
    int n = 0;
    al2->nnb[i] = 0;
    sfree(al2->nb[i]);
    snew(al2->nb[i],al->nnb[idx[i]]);
    
    for(k=0;k<al->nnb[idx[i]];k++){
      if(!polFlag){
        if(strcmp(al->symbol[al->nb[idx[i]][k]],"H")!= 0){
          al2->nb[i][n] = new[al->nb[idx[i]][k]];
          al2->nnb[i]+=1;
          n++;
        }
      }
      else {
        if(strcmp(al->symbol[al->nb[idx[i]][k]],"H")!= 0 ||
           is_polar_H(al,al->nb[idx[i]][k])){
          al2->nb[i][n] = new[al->nb[idx[i]][k]];
          al2->nnb[i]+=1;
          n++;
        }
      }
    }
  }

  al2->nvdw = al->nvdw;

  snew(al2->vdwtab,al->nvdw);
  snew(al2->vdw14tab,al->nvdw);
  for(i=0;i<al->nvdw;i++){
    snew(al2->vdwtab[i],al->nvdw);
    snew(al2->vdw14tab[i],al->nvdw);
  }
  for(i=0;i<al->nvdw;i++){
    for(k=0;k<al->nvdw;k++){
      al2->vdwtab[i][k] = al->vdwtab[i][k];
      al2->vdw14tab[i][k] = al->vdw14tab[i][k];
    }
  }



  /* hbonds */

  al2->nhbonds = al->nhbonds;
  snew(al2->hbonds,al2->nhbonds);
  
  for(i=0;i<al2->nhbonds;i++){
    al2->hbonds[i] = copy_hbond(al->hbonds[i]);
    if(polFlag)
      al2->hbonds[i]->don = new[al->hbonds[i]->don];
    else
      al2->hbonds[i]->don = -1;
    
    al2->hbonds[i]->acc = new[al->hbonds[i]->acc];
    al2->hbonds[i]->don_b = new[al->hbonds[i]->don_b];
    al2->hbonds[i]->acc_b = new[al->hbonds[i]->acc_b];
  }
      

  return al2;


}

void write_dihed(char *fname, t_atomlist *al, t_dihed *dihed)
{
  FILE *fp = ffopen(fname,"w");
  int i;
  char rot[STRLEN];
  char xx[STRLEN];
  
  int at1, at2, at3, at4;
  for(i=0;i<dihed->n;i++){
    at1 = dihed->at1[i];
    at2 = dihed->at2[i];
    at3 = dihed->at3[i];
    at4 = dihed->at4[i];
    if(dihed->flag[i]) strcpy(rot,"no");
    else strcpy(rot,"yes");
    if(IS_PSI(al, at2, at3)) strcpy(xx,"PSI");
    else if(IS_PHI(al, at2, at3)) strcpy(xx, "PHI");
    else if(IS_OMEGA(al, at2, at3)) strcpy(xx, "OMEGA");
    else strcpy(xx, "SC");
    
    fprintf(fp, "%d %s %s %d %s %s %d %s %s %d %s %s -> %s %s\n",
            al->resid[at1], al->name[at1], al->resname[at1],
            al->resid[at2], al->name[at2], al->resname[at2],
            al->resid[at3], al->name[at3], al->resname[at3],
            al->resid[at4], al->name[at4], al->resname[at4],rot,xx);
  }
}

/*=============================================================*/
/*=============================================================*/
/*=============================================================*/
/*=============================================================*/

int main(int argc, char **argv)

/*=============================================================*/
/*                                                             */
/*                   PROGRAM tDIST                             */
/*    Constraint definiton for macromolecular                  */
/*                   structures                                */
/*                                                             */
/*=============================================================*/


{
  static char *desc[] = {  

    "tdist analysis a protein structure and",
    "defines geometrical constraints like bonds,"
    "angles, hydrogen bonds.......",
    "The output file tdist.dat contains the defined",
    "constraints. The coordinates, atom types, and",
    "lots of other stuff is stored in a topology file",
    "tdist.ctp. tdist requires an input pdb file and",
    "a parameter file containing information about",
    "how constraints should be defined. This file is called",
    "input.cpf by default.\n",
    "If a tpr file (use OPLS-AA) is used as input, tdist takes bond, angle\n",
    "and dihedral definitions from the tpr file.\n",
    "Otherwise distance criteria are used.\n"
/*     "Command line options:\n", */
/*     "-def:    expects a value between 0 and 1. This option controles\n", */
/*     "         the number of constraints defined by tdist. A higher\n", */
/*     "         number means more constraints.\n", */
/*     "-grp:    tdist defines groups that are restricted by interactions.\n", */
/*     "         With this option the group definition can by done\n", */
/*     "         manually by reading an index file.\n", */
  };
  
  
  t_atomlist *al;
  t_contab *ct;
  t_resl *rl;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb;
  t_vdwcomb *vdw14comb;
  t_types *tp;
  t_bounds *b;
  t_dihed *impr;
  t_topology top;
  t_atoms *atoms = NULL; 
  char title[STRLEN];
  rvec       box_size;
  rvec *x;
  matrix box;
  bool bVerbose = FALSE;
  
  int i,k,j;
  real plan;
  gmx_rng_t rng;
  real def = 0.5;
  bool bGrp=FALSE;
  atom_id *index;
  char    *grpnames;
  char *posnames;
  int isize;
  int psize = 0;
  atom_id *pos_id;
  int seed = -1;
  int merge = 3;
  bool bNOH = FALSE;
  bool bpH = FALSE;
  bool bIgn = FALSE;
  bool bPosres = FALSE;
  atom_id *flex;
  int flex_size = 0;
  atom_id *target;
  int ntarg = 0;
  char warn[STRLEN];
  bool bTop;
  bool bDihed = FALSE;
  real packed = 1.5;
  char logstr[STRLEN];
  
  
  t_idxgroups *resgroups = idx_init();
  
  
  t_pargs pa[] = {
/*     { "-v",   FALSE, etBOOL, {&bVerbose},"Make Noise"},  */
    { "-def",   FALSE, etREAL, {&def},"Define how many constraints should be defined (0-1)"},
    { "-seed",   FALSE, etINT, {&seed},"Initial random seed (-1 means make a seed)"},
    { "-merge",   FALSE, etINT, {&merge},"Merge groups that have at least # members in common"},
    { "-resample_badpacked",   FALSE, etREAL, {&packed},"Resample residues with poorer packing than #"},
    { "-H",   FALSE, etBOOL, {&bNOH},"Remove Hydrogens"},
    { "-pH",   FALSE, etBOOL, {&bpH},"Remove non-polar hydrogens"},
    { "-ign",   FALSE, etBOOL, {&bIgn},"Use lib entries for bonds and angles"} 
  };
  
  
  t_filenm fnm[] = {
    { efPDBQT,"-p","protein", ffOPTRD },  
    { efTPS,"-s","topol", ffOPTRD }, 
    { efCPF,"-inp","input", ffREAD},
    { efPDB,"-op","tdist", ffWRITE},
    { efDAT,"-od","tdist", ffWRITE},
    { efCTP,"-top","tdist", ffWRITE},
    { efLOG,"-log","tdist", ffWRITE},
    { efPDB,"-clean","clean", ffOPTWR},
    { efPDB,"-packing","packing", ffOPTWR},
    { efNDX,"-resample","resample", ffOPTWR}, 
    { efNDX,"-flexible","flexible", ffOPTWR}, 
    { efDAT,"-rot","rot.dat", ffOPTWR}, 
    { efDAT,"-mvbb","mvbb.dat", ffOPTWR}, 
    { efNDX,"-grp","groups", ffOPTRD},
    { efNDX,"-grpdef","grpdef", ffWRITE},
    { efNDX,"-pos","posres", ffOPTRD},
    { efCPF,"-out","tdist_out", ffOPTWR},
    { efDAT,"-idd","idd", ffOPTWR}, 
    { efDAT,"-hbonds","hbonds", ffWRITE},
    { efNDX,"-flex","flex", ffOPTRD},
    { efNDX,"-tidx","t_index.ndx", ffOPTRD},
    { efDAT,"-ef","exfo.dat",ffOPTRD},
    { efDAT,"-add","conect.dat",ffOPTRD},
    { efDAT,"-dihed","dihed.dat",ffOPTWR}
  };
  
#define NFILE asize(fnm)
  
/*   CopyRight(stderr,argv[0]); */
  cnc_copyright(argv[0]);
  
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  
  int params = 1;/* get_vdw_params(); */
/*   fprintf(stderr,"Selected parameter set %d\n",params); */
  int bonds = 1;/* get_bond_lib(); */
/*   fprintf(stderr,"Selected parameter set %d\n",bonds); */
  
  
  if(opt2bSet("-pos",NFILE,fnm)){
    fprintf(stderr,"Select a group for fixed positions\n");
    get_index(atoms,opt2fn("-pos",NFILE,fnm),1,&psize,&pos_id,&posnames);
    bPosres = TRUE;
  }
  
  if(opt2bSet("-flex",NFILE,fnm)){
    fprintf(stderr,"Select a group for unconstrained atoms\n");
    get_index(atoms,opt2fn("-flex",NFILE,fnm),1,&flex_size,&flex,&posnames);
  }
  
  if(opt2bSet("-tidx",NFILE,fnm)){
    fprintf(stderr,"Select a group of target atoms\n");
    get_index(atoms,opt2fn("-tidx",NFILE,fnm),1,&ntarg,&target,&posnames);
  }
  if(opt2bSet("-dihed",NFILE,fnm)) {
    bDihed = TRUE;
  }
  
/*   t_gridmap *rama = read_rama(); */
/*   printf("%12.8f\n", val_from_2d_table(rama, -54, -17)); */

/*   exit(0); */
  
  
  FILE *log = ffopen(opt2fn("-log",NFILE,fnm),"w");
  char *host = getenv("HOSTNAME");
  char *user = getenv("USER");
  time_t t = time(0);
  struct tm *tt = localtime(&t);
  char *ltime = asctime(tt); 
  fprintf(log,"#=====================================\n");
  fprintf(log,"# tDIST LOG FILE\n");
  fprintf(log,"#=====================================\n");
  fprintf(log,"# USER.......: %s\n",user);
  fprintf(log,"# TIME.......: %s",ltime);
  fprintf(log,"# HOST.......: %s\n",host);
  fprintf(log,"# CMDLINE....: %s\n",command_line());
#if (defined BUILD_MACHINE && defined BUILD_TIME && defined BUILD_USER)
  fprintf(log,
          "# tdist was built on %s by\n"
          "# %s (%s)\n",BUILD_TIME,BUILD_USER,BUILD_MACHINE);
#endif
  fprintf(log,"#=====================================\n");
  if(seed == -1)
    seed = gmx_rng_make_seed();
  rng = gmx_rng_init(seed);
  papers_log(log);
  sprintf(logstr,"tCNC__log_> Initial random seed = %d\n",seed);
  
  fflush(log);
  t_input *inp = read_cnc_input(log,opt2fn("-inp",NFILE,fnm));
  sprintf(logstr,"tCNC__log_> Input parameters:\n");
  CNClog(log,logstr);
  print_input_log(log,inp); 
  print_input_log(stderr,inp); 


  
  if(opt2bSet("-out",NFILE,fnm)){
    sprintf(logstr,"tCNC__log_> Writing parameter file: %s\n",opt2fn("-out",NFILE,fnm));
    CNClog(log,logstr);
    print_input_file(ffopen(opt2fn("-out",NFILE,fnm),"w"),inp); 
  }
  
  if(opt2bSet("-p",NFILE,fnm)){
    sprintf(logstr,"tCNC__log_> Reading AutoDock pdbqt file: %s\n",opt2fn("-p",NFILE,fnm));
    CNClog(log,logstr);
    al = read_pdbqt(opt2fn("-p",NFILE,fnm),log,bVerbose);
    bTop = FALSE;
    get_symbol2(log,al);
  }
  else{
    sprintf(logstr,"tCNC__log_> Reading structure file: %s\n",opt2fn("-s",NFILE,fnm));
    CNClog(log,logstr);
    read_tps_conf(opt2fn("-s",NFILE,fnm),title,&top,&x,NULL,
                  box,TRUE);
    bTop=fn2bTPX(opt2fn("-s",NFILE,fnm));
    al = al_from_atoms(&(top.atoms),x);
    get_symbol(log,al);
  }
  sprintf(logstr,"tCNC__log_> Created atomlist with %d atoms\n",al->natoms);
  CNClog(log,logstr);
  if(bTop) 
  {
    sprintf(logstr,"tCNC__log_> Bonds are taken from topology->idef\n");
    CNClog(log,logstr);
    bonds_from_idef(al,&(top.idef));
  }
  if(psize != 0)
  {
    sprintf(logstr,"tCNC__log_> Applying position constraints on %d atoms\n",psize);
    CNClog(log,logstr);
  }
  posr2al(al,pos_id,psize);   
  if (flex_size != 0)
  {
    sprintf(logstr,"tCNC__log_> Defining %d flexible atoms\n",flex_size);
    CNClog(log,logstr);
  }
  flex2al(al,flex,flex_size); 
  
  /* get the groups from index file */
  
  if(opt2bSet("-grp",NFILE,fnm)){
    atoms = al2atoms(atoms,al,&x);
    fprintf(stderr,"How many groups do you want to select? ");
    int ngrp;
    fscanf(stdin,"%d",&ngrp);
    int count = 0;
    while(count < ngrp){
      fprintf(stderr,"\nSelect a group:\n");
      get_index(atoms,opt2fn("-grp",NFILE,fnm),1,&isize,&index,&grpnames);
      resgroups->n+=1;
      resgroups = idx_realloc(resgroups,resgroups->n);
      for(i=0;i<isize;i++){
        add_to_group(resgroups,resgroups->n-1,al->resid[index[i]]-1);
      }
      count++;
    }
    sfree(index);
    sfree(grpnames);
  }


  ct = contab_init();
  ct=contab_realloc(ct,al->natoms); 
  if(opt2bSet("-add",NFILE,fnm))
    read_additional_connections(opt2fn("-add",NFILE,fnm),al,ct);

  t_gridmap *gp = gridmap_init(); 
  gp = nb(al, gp,9);   
  if(bTop)
  {
    fill_neighborlist(log, al, gp, 9.,UPDATE_NONBONDED);   
  }
  else 
  {
    fill_neighborlist(log, al, gp, 9.,MAKE_FULL_NEIGHBORLIST);   
  }
  
  gp = reset_grid(gp);    
/*   dump_nl(stdout,al,4); */
/*   exit(0); */
  
/*   nb_search(al,ct,0.8,!bTop); */

  sprintf(logstr,"tCNC__log_> Assigning parameters\n");
  CNClog(log,logstr);
  rename_at(al);
  get_order(al);
  vdwcomb = read_vdwcomb(params,FALSE);
  vdw14comb = read_vdwcomb(params,TRUE);
  tp = read_atom_types(params);
  get_types(log,al,tp);
  vdw = read_vdw_radii(params);  
  get_vdw_radii(log,al,vdw,vdwcomb,vdw14comb);
  
  get_hybrid2(al,tp);
  renumber_atoms(al,1);
  renumber_residues(al,1);
/*   get_rosetta_types(al); */
  get_cnc_solv(al);
  
  rl = resl_init();
  rl->nres = count_res(al);
  rl = resl_realloc(rl,rl->nres); 
  fill_resl(rl,al);
  char seq[rl->nres];
  sequence(rl,seq);

  
  
  if(opt2bSet("-idd",NFILE,fnm)){
    t_histogram *h = new_histogram(0,150,1);
    FILE *xx = ffopen(opt2fn("-idd",NFILE,fnm),"w");
    fprintf(xx,"# Interatomic distances\n");
    for(i=0;i<al->natoms-1;i++){
      for(k=i+1;k<al->natoms;k++){
        real d = DIST(al,i,k);
        add_value(h,d,FALSE,0);
        
/*         fprintf(xx, "%g\n",d); */
      }
    }
    norm_histogram(h);
    print_histogram(xx,h);
    
  }

  if(opt2bSet("-clean",NFILE,fnm)){
    sprintf(logstr,"tCNC__log_> Writing structure file: %s\n",opt2fn("-clean",NFILE,fnm));
    CNClog(log,logstr);
    write_pdb(al,opt2fn("-clean",NFILE,fnm));
    exit(0);
  }
  
  
  
  
  t_idxgroups *don = idx_init();
  t_idxgroups *acc = idx_init();
  don = idx_realloc(don,1);
  acc = idx_realloc(acc,1);
  don->n = acc->n = 1;
  t_idxgroups *phob = idx_init();
  phob = idx_realloc(phob,1);
  
  
  /* read lib entries */
  
  t_bondlib *blib = read_bonds_and_angles(bonds);
  
  
/*   if(inp->use_hbonds) */
    hbonds(al,don,acc);
  
/*   if(inp->use_hydrophobics) { */
    FILE *phobics = cnclib("Hydrophobic.dat");
    hydrophobics(al,phob,phobics);
/*   } */
  
  
  sprintf(logstr,"tCNC__log_> Processing structure.............: %s\n",opt2fn("-s",NFILE,fnm));
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of atoms..................: %d\n",al->natoms);
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of residues...............: %d\n",rl->nres);
  CNClog(log,logstr);
  print_sequence_log(stderr,seq, rl->nres);  
  print_sequence_log(log,seq, rl->nres);  
  
  t_contab *rct = contab_init();  
  rct = contab_realloc(rct,rl->nres);  
  
  t_excl *ex = excl_init();
  t_force *fo = force_init();
  if(opt2bSet("-ef",NFILE,fnm)){
    sprintf(logstr,"tCNC__log_> Reading enforced/exlcuded constraints: %s\n",opt2fn("-ef",NFILE,fnm));
    CNClog(log,logstr);
    read_exfo(opt2fn("-ef",NFILE,fnm),ex,fo);
  }
  
  
  int caid[rl->nres];  
  int nca;  
  
  FILE *hb = NULL;
  
  if(inp->use_hbonds){
    sprintf(logstr,"tCNC__log_> Analyzing hydrogen bonds\n");
    CNClog(log,logstr);
    int bbhb, bbhb_rej,scbbhb, scbbhb_rej, scschb, scschb_rej;
    do_hbonds(al,rct,NULL,inp->use_hbond_protection,
              inp->hbond_max_dist, inp->hbond_min_angle,
              inp->solvation_max,inp->solvation_rad,
              &bbhb,&bbhb_rej,&scbbhb,&scbbhb_rej, 
              &scschb,&scschb_rej);  
    sprintf(logstr,"tCNC__log_> Number of hbond donors...........: %d\n",don->natoms[0]);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Number of hbond acceptors........: %d\n",don->natoms[0]);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Number of BB-BB hydrogen bonds...: %d\n",bbhb);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Rejected BB-BB hydrogen bonds....: %d\n",bbhb_rej);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Number of SC-BB hydrogen bonds...: %d\n",scbbhb);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Rejected SC-BB hydrogen bonds....: %d\n",scbbhb_rej);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Number of SC-SC hydrogen bonds...: %d\n",scschb);  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Rejected SC-SC hydrogen bonds....: %d\n",scschb_rej);  CNClog(log,logstr);
/*     hb = ffopen(opt2fn("-hbonds",NFILE,fnm),"w"); */
/*     print_hbonds_to_file(hb,al); */
  }
  
  t_idxgroups *phb = idx_init();
  
  if(inp->use_hydrophobics) {
    sprintf(logstr,"tCNC__log_> Analyzing hydrophobic clusters\n");  CNClog(log,logstr);
    sprintf(logstr,"tCNC__log_> Number of hydrophobic atoms......: %d\n",phob->natoms[0]);  CNClog(log,logstr);
    do_hydrophobics(al,rct,inp->hphob_dist,phb);    

    phb = merge_groups(phb,1,FALSE); 
    phb = remove_redundant_groups(phb); 
    sort_groups(phb); 

    /* dirty! we keep only those groups with more than 3 members */

    for(i=0;i<phb->n;i++){
      if(phb->natoms[i] < 4)
      {
        phb->n = i;
        break;
      }
    }


/*     for(i=0;i<phb->n;i++){ */
/*       for(k=0;k<phb->natoms[i]-1;k++){ */
/*         for(j=k+1;j<phb->natoms[i];j++){ */
/*           add2contab(rct,phb->atoms[i][k],phb->atoms[i][j],ePHO,ePHO,0); */
/*         } */
/*       } */
/*     } */
    /********************************/
    sprintf(logstr,"tCNC__log_> Defined %d hydrophobic clusters\n", phb->n);  CNClog(log,logstr);

    write_group_script("phobcs.pml","resi",phb); 
/*     exit(0); */
    
    
/*     if(opt2bSet("-hbonds",NFILE,fnm)){ */
/*       print_hydrophobics_to_file(hb,al,phb); */
/*     } */
  }
  
  t_atomlist *al2 = NULL;
  
  if(bNOH || bpH) {
    bool polFlag = FALSE;
    if(bpH) polFlag = TRUE;
    sprintf(logstr,"tCNC__log_> Removing Hydrogen Atoms........\n");  CNClog(log,logstr);
    al = remove_H(al,polFlag);
    sprintf(logstr,"tCNC__log_> Now there are %d atoms........\n",al->natoms);  CNClog(log,logstr);
    renumber_atoms(al,1);
    fill_resl(rl,al);
  }
  
/*   if(inp->use_hbonds){ */
/*     hb = ffopen(opt2fn("-hbonds",NFILE,fnm),"w"); */
/*     print_hbonds_to_file(hb,al); */
/*   } */
  
/*   if(inp->use_hydrophobics)  */
/*       print_hydrophobics_to_file(hb,al,phb); */


  

/*   char seq[rl->nres]; */
/*   sequence(rl,seq); */
/*   fprintf(stderr,"tCNC__log_> %s\n",seq);    */
  
  t_dihed *dihed = dihed_init();
  get_dihedrals(al,dihed);
  
  
  
  t_bondlist *bl = bl_init();
  rotatable_bonds(al,bl);
/*   fprintf(stderr,"number of bonds..................: %d\n",bl->n); */
  do_cov(al,bl,rct,ex);
  
  
  
  
  /* get planar groups */
  
  sprintf(logstr,"tCNC__log_> Analyzing planar groups\n");  CNClog(log,logstr);

  t_idxgroups *pln = idx_init();

  int plviol = get_planar_groups(log,al,rl,pln,bIgn);
  
  /* check non-protein residues */
  
  for(i=0;i<bl->n;i++){
    if(!IsProtein(al,bl->at1[i])){
      if(bl->type[i] == DOUBLE ||
         bl->type[i] == CNDOUBLE ||
         bl->type[i] == CNOMEGA ||
         bl->type[i] == NNDOUBLE ||
         bl->type[i] == TRIPLE ||
         bl->type[i] == CNTRIPLE){
        pln->n++;
        pln = idx_realloc(pln,pln->n);
        add_to_group(pln,pln->n-1,bl->at1[i]);
        add_to_group(pln,pln->n-1,bl->at2[i]);
        pln->flag[pln->n-1] = TRUE;
        for(k=0;k<al->nbonds[bl->at1[i]];k++){
          add_to_group(pln,pln->n-1,al->bonds[bl->at1[i]][k]);
        }
        for(k=0;k<al->nbonds[bl->at2[i]];k++){
          add_to_group(pln,pln->n-1,al->bonds[bl->at2[i]][k]);
        }
        
/*         if(!is_planar_group(al,pln->natoms[pln->n-1],pln->atoms[pln->n-1], */
/*                             inp->plan_tol,&plan)){ */
/*           pln->n--; */
/*         } */
/*         else { */
/*           if (plan < 0.001) { */
/*             plan = 0.001; */
/*           } */
          
/*           pln->val[pln->n-1] = plan; */
/*         } */
        
      }
    }
  }
  
  
  
  for(i=0;i<al->natoms;i++){
    if(!IsProtein(al,i) && strcmp(al->hyb[i],"sp2") == 0 
       && al->nbonds[i] >= 2){
      pln->n++;
      pln = idx_realloc(pln,pln->n);
      add_to_group(pln,pln->n-1,i);
      pln->val[pln->n-1] = inp->plan_tol;
      for(k=0;k<al->nbonds[i];k++){
        add_to_group(pln,pln->n-1,al->bonds[i][k]);
      }
    }
  }
  
  pln = remove_redundant_groups(pln);
  
  
  
  /* get rings */
  
  t_idxgroups *rings = idx_init();
  for(i=0;i<al->natoms;i++){
    if(!IsProtein(al,i) && strcmp(al->symbol[i],"H") != 0)
      get_ring(al,i,rings,bl);
  }
  
  rings = remove_redundant_groups(rings);
  
  
  for(i=0;i<rings->n;i++){
    if(is_planar_group(al,rings->natoms[i],rings->atoms[i],inp->plan_tol,&plan)){
      pln->n++;
      pln = idx_realloc(pln,pln->n);
      for(k=0;k<rings->natoms[i];k++){
        add_to_group(pln,pln->n-1,rings->atoms[i][k]);
      }
    }
  }
  
  pln = remove_redundant_groups(pln);
  pln = merge_groups(pln,3,0);  
  pln = remove_redundant_groups(pln);  
  
  char wstr[STRLEN];
  
  for(i=0;i<pln->n;i++){
    is_planar_group(al,pln->natoms[i],pln->atoms[i],inp->plan_tol,&plan);
    if(pln->val[i] < plan) {
      sprintf(wstr,"Planarity violation %d%s \n(%d%s %d%s ..)... %g -> %g\n", 
              al->resid[pln->atoms[i][0]],al->resname[pln->atoms[i][0]], 
              al->id[pln->atoms[i][0]],al->name[pln->atoms[i][0]], 
              al->id[pln->atoms[i][1]],al->name[pln->atoms[i][1]], 
              pln->val[i],plan); 
      CNCwarn(log,wstr); 
      if(plan > inp->plan_tol) plan*=1.01;
      else plan = inp->plan_tol;
      pln->val[i] = plan;
      sprintf(wstr,"Setting tolerance to %g\n",plan);
      CNCwarn(log,wstr);
    }
  }
  
      
/*     if(pln->val[i] <  0.001){ */
/*       pln->val[i] = 0.001; */
/*     } */
/*     is_planar_group(al,pln->natoms[i],pln->atoms[i],inp->plan_tol,&plan); */
/*     if(plan > pln->val[i]){ */
/*       pln->val[i] = plan*1.05; */
/*     } */
/*   } */
  
  
  
  
  for(i=0;i<pln->n;i++){
    int plcheck = 0;
    for(k=0;k<pln->natoms[i];k++){
      if(al->isposres[pln->atoms[i][k]]) plcheck++;
      
      if(ntarg){  
        if(id_in_index(pln->atoms[i][k],ntarg,target)) plcheck++; 
      } 
      
      
    }
    if(pln->natoms[i] == plcheck) pln->flag[i] = FALSE;
    else pln->flag[i] = TRUE;
  }
  
  
  
  
  int ndih = 0;
  for(i=0;i<dihed->n;i++){
    for(k=0;k<bl->n;k++){
      if((dihed->at2[i] == bl->at1[k] &&
          dihed->at3[i] == bl->at2[k]) ||
         (dihed->at3[i] == bl->at1[k] &&
          dihed->at2[i] == bl->at2[k])){
        dihed->flag[i] = bl->restricted[k];
        
      }
    }
    if(dihed->flag[i]) ndih++;
    }


  sprintf(logstr,"tCNC__log_> Analyzing dihedral angles\n");  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of dihedrals...........................: %d\n",dihed->n);  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Number of chemically restricted  dihedrals....: %d\n",ndih);  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Analyzing hydrogen bond restrictions\n");  CNClog(log,logstr);
  
  
  get_restricted(al,rct,rl,dihed);    
  ndih = 0;
  for(i=0;i<dihed->n;i++){
    if(dihed->flag[i]) ndih++;
  }
  sprintf(logstr,"tCNC__log_> Total number of restricted dihedrals...........: %d\n",ndih);  CNClog(log,logstr);
  if(bDihed){
    write_dihed(opt2fn("-dihed",NFILE,fnm), al, dihed);
  }
  
  
  t_idxgroups *imp = idx_init();
  
/*   for(i=0;i<al->natoms;i++){ */
/*     t_idxgroups *tmp = idx_init(); */
/*     tmp = idx_realloc(tmp,1); */
/*     if(al->nbonds[i] > 2){ */
/*       srenew(tmp->atoms[0],4); */
/*       tmp->atoms[0][0] = i; */
/*       tmp->natoms[0] = 4; */
/*       for(k=0;k<3;k++){ */
/*         tmp->atoms[0][k+1] = al->bonds[i][k]; */
/*       } */
/*       if(!is_planar_group(al,tmp->natoms[0],tmp->atoms[0],0.003,&plan) && */
/*          !is_omega_group(al,tmp->atoms[0],tmp->natoms[0])) */
/*       { */
/*         imp->n+=1; */
/*         imp = idx_realloc(imp,imp->n); */
/*         tmp->val[0] = DIHED(al,tmp->atoms[0][0],tmp->atoms[0][1],tmp->atoms[0][2],tmp->atoms[0][3]); */
/*         copy_group(tmp,0,imp,imp->n-1); */
/*       } */
/*     } */
/*     free_groups(tmp); */
/*   } */


  for(i=0;i<al->natoms;i++){
    t_idxgroups *tmp = idx_init();
    tmp = idx_realloc(tmp,1);
    if(strcmp(al->hyb[i],"sp3") == 0 && al->nbonds[i] > 2){
      srenew(tmp->atoms[0],4);
      tmp->atoms[0][0] = i;
      tmp->natoms[0] = 4;
      for(k=0;k<3;k++){
        tmp->atoms[0][k+1] = al->bonds[i][k];
      }
/*       if(!is_planar_group(al,tmp->natoms[0],tmp->atoms[0],0.003,&plan) && */
/*          !is_omega_group(al,tmp->atoms[0],tmp->natoms[0])) */
/*       { */
        imp->n+=1;
        imp = idx_realloc(imp,imp->n);
        tmp->val[0] = DIHED(al,tmp->atoms[0][0],tmp->atoms[0][1],tmp->atoms[0][2],tmp->atoms[0][3]);
        copy_group(tmp,0,imp,imp->n-1);
/*       } */
    }
    free_groups(tmp);
  }




  sort_group_by_order(al,imp);  



  
  for(i=0;i<imp->n;i++){
    int imcheck = 0;
    for(k=0;k<imp->natoms[i];k++){
      if(al->isposres[imp->atoms[i][k]]) imcheck++;
    }
    if(imcheck == imp->natoms[i]) imp->flag[i] = FALSE; 
/*     if(imcheck) imp->flag[i] = FALSE; */
    else imp->flag[i] = TRUE;
  }
  
   sort_group_by_posres(al,imp);   
/*   sort_group_by_order(al,imp);   */

   sprintf(logstr,"tCNC__log_> Analyzing chirality\n");  CNClog(log,logstr);
   sprintf(logstr,"tCNC__log_> Number of impropers..............: %d\n",imp->n);  CNClog(log,logstr);
  
  
  
  
  if(fo->n > 0){
    for(i=0;i<fo->n;i++){
      sprintf(logstr,"tCNC__log_> Applying forced constraints\n");  CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Will force constraints...........: %d%s -- %d%s  %g  %g\n",fo->id1[i],
              rl->resname[fo->id1[i]-1],fo->id2[i],rl->resname[fo->id2[i]-1],
              fo->lb[i],fo->ub[i]);  CNClog(log,logstr);
      add2contab(rct,fo->id1[i]-1,fo->id2[i]-1,eFO,eFO,0);
    }
  }
  if(ex->n > 0){
    for(i=0;i<ex->n;i++){
      sprintf(logstr,"tCNC__log_> Applying exclusions\n");  CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Will exclude constraints.........: %d%s -- %d%s\n",ex->id1[i],
              rl->resname[ex->id1[i]-1],ex->id2[i],rl->resname[ex->id2[i]-1]);  CNClog(log,logstr);
      for(k=0;k<al->nhbonds;k++){
        if(is_excluded(ex,al,al->hbonds[k]->don,al->hbonds[k]->acc)){
          sprintf(logstr,"tCNC__log_> Excluding H-Bond %d-%s -> %d-%s\n",
                  al->id[al->hbonds[k]->don],al->name[al->hbonds[k]->don],
                  al->id[al->hbonds[k]->acc],al->name[al->hbonds[k]->acc]);  CNClog(log,logstr);
          al->hbonds[k]->constr = FALSE;
        }
      }
    }
  }
  for(i=0;i<rct->n;i++){
    for(k=0;k<rct->ncon[i];k++){
      if(is_excluded(ex,al,rl->j0[i],rl->j0[rct->con[i][k]])){
        sprintf(logstr,"tCNC__log_> Changing type %d to type %d\n",rct->type[i][k],eEX);  CNClog(log,logstr);
        rct->type[i][k] = eEX;
        
      }
    }
  }
  
  
  
/*   if(inp->use_hbonds){ */
/*     hb = ffopen(opt2fn("-hbonds",NFILE,fnm),"w"); */
/*     print_hbonds_to_file(hb,al); */
/*   } */
  
/*   if(inp->use_hydrophobics)  */
/*       print_hydrophobics_to_file(hb,al,phb); */

  
  sprintf(logstr,"tCNC__log_> Writing structure file: %s\n",opt2fn("-op",NFILE,fnm));  CNClog(log,logstr);
  write_pdb(al,opt2fn("-op",NFILE,fnm));
  
  
  if(inp->use_network){
    if(!bGrp){
      sprintf(logstr,"tCNC__log_> Analyzing hydrogen bond network\n");  CNClog(log,logstr);
      do_hbond_resnet(rct,resgroups);  
/*        do_resnet(rct,resgroups);   */
/*       print_connections(al,rct,rl);   */
      
      resgroups = remove_redundant_groups(resgroups);  
      resgroups = merge_groups(resgroups,merge,TRUE);      
      for(i=0;i<rl->nres;i++){
        rl->ngrps[i] = 0;
        for(k=0;k<resgroups->n;k++){
          for(j=0;j<resgroups->natoms[k];j++){
            if(i == resgroups->atoms[k][j]){
              rl->ngrps[i]+=1;
              srenew(rl->grps[i],rl->ngrps[i]);
              rl->grps[i][rl->ngrps[i]-1] = k;
            }
          }
        }
      }
      sprintf(logstr,"tCNC__log_> Reduced hydrogen bond network to %d groups\n",resgroups->n);  CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Analyzing core\n");  CNClog(log,logstr);
      resgroups = analyze_core(phb,resgroups,rl);
/*       exit(0); */
      
/*       sprintf(logstr,"Merging with merge = %d\n",merge); */
      /* add groups from hydrophobic clusters */
/*       for(i=0;i<phb->n;i++){  */
/*         if(phb->natoms[i] > 3 ){  */
/*           resgroups->n++;  */
/*           resgroups = idx_realloc(resgroups,resgroups->n);  */
/*           copy_group(phb,i,resgroups,resgroups->n-1);  */
/*         }  */
/*       }   */

      sprintf(logstr,"tCNC__log_> Merging residue groups\n");  CNClog(log,logstr);
      resgroups = merge_groups(resgroups,3,TRUE);       
      sprintf(logstr,"tCNC__log_> Reduced residue network to %d groups\n",resgroups->n);


      
/*       resgroups = get_hydrophobic_patches(resgroups, phb, rct);  */
      
      
      
      sprintf(logstr,"tCNC__log_> Analyzing sidechain hydrogen bonds\n");  CNClog(log,logstr);
      resgroups =  sidechain_hbonds(log, resgroups, rct); 
      resgroups = remove_redundant_groups(resgroups);         
      sprintf(logstr,"tCNC__log_> Reduced residue network to %d groups\n",resgroups->n);  CNClog(log,logstr);
      sort_groups(resgroups);


/*       for(i=0;i<resgroups->n;i++){ */
/*         print_simple_group(resgroups->atoms[i],resgroups->natoms[i]); */
/*       } */
      
    }   
    
    for(i=0;i<rl->nres;i++){
      rl->ngrps[i] = 0;
      for(k=0;k<resgroups->n;k++){
        for(j=0;j<resgroups->natoms[k];j++){
          if(i == resgroups->atoms[k][j]){
            rl->ngrps[i]+=1;
            srenew(rl->grps[i],rl->ngrps[i]);
            rl->grps[i][rl->ngrps[i]-1] = k;
          }
        }
      }
    }
    for(i=0;i<rl->nres;i++){
      for(k=rl->j0[i];k<=rl->j1[i];k++){
        al->ngrps[k] = rl->ngrps[i];
        al->ngroups++;
        
        snew(al->grpnr[k],al->ngrps[k]);
        for(j=0;j<rl->ngrps[i];j++){
          al->grpnr[k][j] = rl->grps[i][j];
        }
      }
    }
    

    t_block *block;
    char **gnames;
    atom_id *newgrp;
    snew(newgrp,al->natoms);
    int nn= 0;
    char grpname[STRLEN];    
    block = new_block();
    gnames = NULL;
    snew(gnames,1);
    analyse(&(top.atoms),block,&gnames,FALSE,FALSE);

    for(i=0;i<resgroups->n;i++){
      nn = 0;
      sprintf(grpname,"Group_%d",i);
      for(j=0;j<resgroups->natoms[i];j++){
        /* go over residue ids */
        int rid = resgroups->atoms[i][j];
        for(k=rl->j0[rid];k<=rl->j1[rid];k++){
          newgrp[nn] = al->id[k]; 
          nn++; 
        }
      }
/*       printf("adding group %s with %d atoms\n",grpname,nn); */
      add_grp(block,&gnames,nn,newgrp,grpname); 
    }
    write_index(opt2fn("-grpdef",NFILE,fnm),block,gnames); 
    write_group_script("grp.pml","resi",resgroups); 
  }
  
/*   FILE *gp = ffopen("groups.dat","w"); */
/*   for(i=0;i<resgroups->n;i++){ */
/*     fprintf(gp,"[ GROUP %d ]\n",i); */
/*     for(j=0;j<resgroups->natoms[i];j++){ */
/*       fprintf(gp, "%d\n",resgroups->atoms[i][j]); */
/*     } */
/*   } */
  


  occ_to_one(al);
  
/*   fake_charge(al); */

  sprintf(logstr,"tCNC__log_> Writing topology file: %s\n", opt2fn("-top",NFILE,fnm));  CNClog(log,logstr);
  write_cnctop(opt2fn("-top",NFILE,fnm),al,don,acc,phob,pln,imp,inp,
               vdw,vdwcomb,psize,pos_id); 
  
  
  
  t_boundtrack *bt = boundtrack_init();
  bt = boundtrack_realloc(bt,al->natoms);  

  sprintf(logstr,"tCNC__log_> Writing constraint file: %s\n", opt2fn("-od",NFILE,fnm));  CNClog(log,logstr);
  FILE *fp  = ffopen(opt2fn("-od",NFILE,fnm),"w");
  
  sprintf(logstr,"tCNC__log_>--------------------------------------------------------\n");  
  CNClog(log,logstr);
  int ccount = 0;
  
  if(bTop){
    ccount+=write_bonds_from_tpx(fp,log,al,bt,inp,ex,&top,bIgn);
    ccount+=write_angles_from_tpx(fp,log,al,bt,inp,ex,&top,bIgn);
  }
  else{
    ccount+=write_bonds(fp,log,al,bt,inp,ex,blib,bIgn); 
    ccount+=write_angles(fp,log,al,bt,inp,ex,blib,bIgn);  
  }
  
  if(inp->use_hbond_angles && !bNOH)
    ccount+=write_hbond_angles(fp,log,al,bt,ex);   
  ccount+=write_dihedrals(fp,log,al,bt,dihed,inp,ex);    

  ccount+=write_planar(fp,log,al,bt,pln,inp,ex);  

  if(fo->n > 0)
    ccount+=write_forced(fp,log,al,rl,fo,inp,bt);
  

/*   if(inp->use_sidechain){ */
/*     ccount+=write_sidechain(fp,al,rl,bt);      */
/*   } */

  


  if(inp->use_hbonds)
    ccount+=write_hbonds(fp,log,al,bt,ex);   
  if(inp->use_packing) {
    calc_packing(al);
    ccount+=write_packing_constraints(fp,log,al,bt,ex,inp->pack_limit,inp);   
  }

  
  
  if(inp->use_hydrophobics)
    ccount+=write_hydrophobics(fp,log,al,bt,rct,rl,phb,inp,rng,def,ex,bpH || bNOH);  
  
/*   if(inp->use_packing) { */
/*     calc_packing(al); */
/*     ccount+=write_packing_constraints(fp,al,bt,ex,inp->pack_limit,inp);    */
/*   } */
  if(inp->use_network){
    ccount+=write_network_restrictions(fp,log,al,bt,resgroups,rl,
                                       inp->max,rng,
                                       inp,def,40,ex,bpH || bNOH);        
  }
/*   if(inp->use_sidechain){ */
/*     ccount+=write_sidechain(fp,al,rl,bt);      */
/*   } */
  if(inp->use_sidechain){ 
    ccount+=write_sidechain(fp,log,al,rl,bt, 50);      
  } 
  
  if(inp->use_close_pairs) 
    ccount+=write_close_pairs(fp,log,al,bt,inp,inp->min, 
                              inp->close_pairs_fix_dist,ex);   
  if(inp->use_long_range)
    ccount+=long_range(fp,log,al,bt,rl,inp->max,rng,ex,inp);      
  

  if(inp->use_neighbor_ca)
    ccount+=write_neighbor_res(fp, log,al,bt,ex,bIgn);    
  

  
/*   for(i=0;i<rl->nres;i++){ */
/*     fprintf(xx, "resid %d (%s) restr = %d\n",rl->id[i], rl->resname[i], rl->sc_restr[i]); */
/*   } */

  real respack[rl->nres];
  if(opt2bSet("-packing",NFILE,fnm)){ 
    calc_packing(al);
    packing_per_residue(al,rl,respack); 
    write_pdb(al,opt2fn("-packing",NFILE,fnm));
  }


  if(opt2bSet("-flexible",NFILE,fnm)){ 
    FILE *xx = ffopen(opt2fn("-flexible",NFILE,fnm),"w");
    fprintf(xx,"[ resample_posre ]\n");
    int cc = 0;
    for(i=0;i<rl->nres;i++){
      if(rl->ngrps[i] != 0) {
        for(k=rl->j0[i];k<=rl->j1[i];k++){
          fprintf(xx, "%d ",al->id[k]);
          cc++; 
         if(cc % 15 == 0) fprintf(xx,"\n"); 
        }
      }
    }
    fprintf(xx,"\n[ resample_flexible ]\n");
    cc = 0;
    
    for(i=0;i<rl->nres;i++){
      if(rl->ngrps[i] == 0) {
        for(k=rl->j0[i];k<=rl->j1[i];k++){
          fprintf(xx, "%d ",al->id[k]);
          cc++; 
         if(cc % 15 == 0) fprintf(xx,"\n"); 
        }
      }
    }
    
/*     for(i=0;i<al->natoms;i++){ */
/*       if(al->ngrps[i] != 0) { */
/*         fprintf(xx, "%d ",al->id[i]); */
/*         cc++; */
/*         if(cc % 15 == 0) fprintf(xx,"\n"); */
/*       } */
/*     } */
/*     fprintf(xx,"\n[ resample_flexible ]\n"); */
/*     cc = 0; */
/*     for(i=0;i<al->natoms;i++){ */
/*       if(al->ngrps[i] == 0) { */
/*         fprintf(xx, "%d ",al->id[i]); */
/*         cc++; */
/*         if(cc % 15 == 0) fprintf(xx,"\n"); */
/*       } */
/*     } */
    fprintf(xx,"\n");
    fclose(xx);
  }

/*     t_idxgroups *flx = idx_init(); */
/*     flx->n = 1; */
    
/*     for(i=0;i<rl->nres;i++){ */
/*       if(rl->ngrps[i] == 0) { */
/*         add_to_group(flx,0,i+1); */
/*       } */
/*     } */
/*     printf("natoms = %d\n",flx->natoms[0]); */
    
/*     write_group_script(opt2fn("-flexible",NFILE,fnm),"resi",flx); */
    

  

      
  for(i=0;i<rct->n;i++){
    rl->sc_restr[i] = FALSE;
    if(respack[i] < packed ) rl->sc_restr[i] = TRUE;

    for(k=0;k<rct->ncon[i];k++){
      switch(rct->type[i][k]){
        case eSUL:
        case eSSHD:
        case eSSHA:
        case eSBH_SD:
        case eSBH_SA:          
/*         case ePHO: */
          rl->sc_restr[i] = TRUE;
          break;
      }
    }
  }
/*   if(opt2bSet("-packing",NFILE,fnm)){  */
/*     calc_packing(al); */
/*     real respack[rl->nres]; */
/*     packing_per_residue(al,rl,respack);  */
/*     write_pdb(al,opt2fn("-packing",NFILE,fnm)); */
/*   } */

  if(opt2bSet("-resample",NFILE,fnm)){
    FILE *xx = ffopen(opt2fn("-resample",NFILE,fnm),"w");
    fprintf(xx,"[ resample_posre ]\n");
    int cc = 0;
    for(i=0;i<rl->nres;i++){
      for(k=rl->j0[i];k<=rl->j1[i];k++){
        if(al->order[k] < 2 || rl->sc_restr[i]) {
          fprintf(xx, "%d ",al->id[k]);
          cc++;
          if(cc % 15 == 0) fprintf(xx,"\n");
        }
      }
    }
    fprintf(xx,"\n");
    fclose(xx);
  }
  
/*   if(opt2bSet("-packing",NFILE,fnm)){  */
/*     calc_packing(al); */
/*     real respack[rl->nres]; */
/*     packing_per_residue(al,rl,respack);  */
    
/*     for(i=0;i<al->natoms;i++){ */
/*       int nex = 1+al->nbonds[i]+al->nb13[i]+al->nb14[i]; */
/*       int excluded[nex]; */
/*       excluded[0] = i; */
/*       for(k=0;k<al->nbonds[i];k++){ */
/*         excluded[k+1] = al->bonds[i][k]; */
/*       } */
/*       for(k=0;k<al->nb13[i];k++){ */
/*         excluded[al->nbonds[i]+k+1] = al->b13[i][k]; */
/*       } */
/*       for(k=0;k<al->nb14[i];k++){ */
/*         excluded[al->nbonds[i]+al->nb13[i]+k+1] = al->b14[i][k]; */
/*       } */
      
/*       if(strcmp(al->symbol[i],"H")!=0){ */
/*         real ps = ray_search(al, al->x[i],excluded,nex); */
/*         al->bfac[i] = ps; */
/*       } */
/*       else al->bfac[i] = 0; */
      
/*     } */
/*     write_pdb(al,opt2fn("-packing",NFILE,fnm)); */
/*   } */
  
/*     t_histogram *h = new_histogram(0,150,1); */
/*     FILE *xx = ffopen(opt2fn("-idd",NFILE,fnm),"w"); */
/*     fprintf(xx,"# Interatomic distances\n"); */
/*     for(i=0;i<al->natoms-1;i++){ */
/*       for(k=i+1;k<al->natoms;k++){ */
/*         real d = DIST(al,i,k); */
/*         add_value(h,d,FALSE,0); */
        
/*       } */
/*     } */
/*     norm_histogram(h); */
/*     print_histogram(xx,h); */
    
/*   } */
  if(opt2bSet("-rot",NFILE,fnm)){
  
    t_rotdata *rt = read_rotation_lib(); 
    FILE *xx = ffopen(opt2fn("-rot",NFILE,fnm),"w"); 
    write_rotations(xx,al,rl,rt); 
  }

  if(opt2bSet("-mvbb",NFILE,fnm)){
    FILE *xx = ffopen(opt2fn("-mvbb",NFILE,fnm),"w"); 
    for(i=1;i<5;i++){
      make_bb_rotations(xx,al,rl,i,bVerbose);
    }
    
  }
  
  
/*   t_idxgroups *rot = read_rotations("rot.dat"); */
/*   for(i=0;i<rot->n;i++){ */
/*     for(k=0;k<rot->natoms[i];k++){ */
/*       printf("%d %s ",al->id[rot->atoms[i][k]],al->name[rot->atoms[i][k]]); */
/*     } */
/*     printf("\n"); */
/*   } */
  

/*   calc_hbond_packing(al); */
  
/*   if(inp->use_hbonds){ */
    hb = ffopen(opt2fn("-hbonds",NFILE,fnm),"w");
    print_hbonds_to_file(hb,al);
/*   } */
  
/*   if(inp->use_hydrophobics)  */
      print_hydrophobics_to_file(hb,al,phb);
  sprintf(logstr,"tCNC__log_>--------------------------------------------------------\n");

  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Total sum of constraints %d\n",ccount);  CNClog(log,logstr);

/*   fprintf(stderr,"Pl-check____> %d planar groups violated library entries\n", plviol); */
  

  
  FILE *of = ffopen("constr.dat","w"); 
  print_constr_per_atom(of,bt,al); 
  print_contab(rct,"contab.dat"); 


  sprintf(logstr,"tCNC__log_>--------------------------------------------------------\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Job terminated sucessfully\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> Exiting\n");
  CNClog(log,logstr);
  
  return 0;
}
  

