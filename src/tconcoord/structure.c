#include <tconcoord.h>

/*=============================================================*/

bool is_polar_H(t_atomlist *al, int i)
{
  if(strcmp(al->symbol[i],"H") == 0 &&
     strcmp(al->symbol[al->bonds[i][0]],"C") != 0)
    return TRUE;
  else return FALSE;
}

/*============================================================*/
bool IsProtein(t_atomlist *al, int k)
{
  int i;
  for(i=0;i<asize(prot_res);i++){
    if(!strcmp(al->resname[k],prot_res[i])) 
      return TRUE;
  }
  return FALSE;
}
/*============================================================*/
bool IsNucac(t_atomlist *al, int k)
{
  int i;
  for(i=0;i<asize(nucac_res);i++){
    if(!strcmp(al->resname[k],nucac_res[i])) 
      return TRUE;
  }
  return FALSE;
}
/*============================================================*/
bool IsIon(t_atomlist *al, int k)
{
  int i;
  for(i=0;i<asize(ion_res);i++){
    char dum1[STRLEN];
    char dum2[STRLEN];
    strcpy(dum1,al->resname[k]);
    strcpy(dum2,ion_res[i]);
    trim(dum1);
    trim(dum2);
/*     printf("dum1 = %s dum2 = %s\n",dum1,dum2); */
    
    if(!strcmp(al->resname[k],ion_res[i]) ||
       !strcmp(dum1,dum2)) 
      return TRUE;
  }
  return FALSE;
}
/*============================================================*/
bool IsSol(t_atomlist *al, int k)
{
  int i;
  for(i=0;i<asize(sol_res);i++){
    if(!strcmp(al->resname[k],sol_res[i])) 
      return TRUE;
  }
  return FALSE;
}


/*=============================================================*/
void get_ring(t_atomlist *al, int i, t_idxgroups *tr, t_bondlist *bl)

/* check whether atom i is part of a ring */

{
  int k,j,l,m;
  int at1,at2;
  int count = 0;
  int mem[1000];
  for(j=0;j<1000;j++){
    mem[j] = -1;
  }
  /* check for six rings */
  int nat = 0;
  bool check = FALSE;
  for(k=0;k<al->nb14[i];k++){
    at1 = al->b14[i][k];
    count = 0;
    nat = 0;
    for(j=0;j<al->nbonds[at1];j++){
      for(l=0;l<al->nb13[i];l++){
        if(al->bonds[at1][j] == al->b13[i][l] &&
           strcmp(al->symbol[al->bonds[at1][j]],"H") != 0) {
          mem[nat] = al->bonds[at1][j];
          nat++;
          count++;
        }
      }
    }
    if(count == 2) {
      restrict_bond(bl,at1,mem[0]);
      restrict_bond(bl,at1,mem[1]);
      mem[nat] = at1;nat++;
      mem[nat] = i;nat++;
      /* now do reverse to get the last 2 members */
      for(l=0;l<al->nbonds[i];l++){
        for(j=0;j<al->nb13[at1];j++){
          if(al->bonds[i][l] == al->b13[at1][j] &&
             strcmp(al->symbol[al->bonds[i][l]],"H") != 0){
            mem[nat] = al->bonds[i][l];
            restrict_bond(bl,i,al->bonds[i][l]);
            nat++;
          }
        }
      }
      restrict_bond(bl,mem[0],mem[4]);
      restrict_bond(bl,mem[0],mem[5]);
      restrict_bond(bl,mem[1],mem[4]);
      restrict_bond(bl,mem[2],mem[5]);
      
      tr->n++;
      tr = idx_realloc(tr,tr->n);
      for(l=0;l<nat;l++){
        add_to_group(tr,tr->n-1,mem[l]);
      }
      /* add connected atoms to ring */
      int size = tr->natoms[tr->n-1];
      int dum[size];
      for(l=0;l<size;l++){
        dum[l] = tr->atoms[tr->n-1][l];
      }
      for(l=0;l<size;l++){
        for(m=0;m<al->nbonds[dum[l]];m++){
          add_to_group(tr,tr->n-1,al->bonds[dum[l]][m]);
        }
      }
    }
  }
  
  /* check for five rings */

  for(j=0;j<1000;j++){
    mem[j] = -1;
  }
  nat = 0;
  check = 0;
  for(k=0;k<al->nb13[i];k++){
    at1 = al->b13[i][k];
    for(j=0;j<al->nbonds[at1];j++){
      at2 = al->bonds[at1][j];
      for(l=0;l<al->nb13[i];l++){
        if(al->b13[i][l] == at2 && strcmp(al->symbol[at2],"H") !=0 &&
           strcmp(al->symbol[at1],"H") !=0 && strcmp(al->symbol[i],"H") !=0) {
          check++;
          if(check==1){
            mem[nat] = i;nat++;
          }
          mem[nat] = at1;nat++;
          mem[nat] = at2;nat++;
        }

      }
    }
  }
  
  if(check > 1){
    for(k=0;k<al->nbonds[i];k++){
      for(l=0;l<al->nbonds[mem[1]];l++){
        if(al->bonds[i][k] == al->bonds[mem[1]][l]){
          mem[3] = al->bonds[i][k];
        }
      }
      for(l=0;l<al->nbonds[mem[2]];l++){
        if(al->bonds[i][k] == al->bonds[mem[2]][l]){
          mem[4] = al->bonds[i][k];
        }
      }
    }
    tr->n++;
    tr = idx_realloc(tr,tr->n);
    for(l=0;l<nat;l++){
      add_to_group(tr,tr->n-1,mem[l]);
    }
    for(k=0;k<nat-1;k++){
      for(l=k+1;l<nat;l++){
        if(is_bond(al,mem[k],mem[l])){
          restrict_bond(bl,mem[k],mem[l]);
        }
      }
    }
    
    int size = tr->natoms[tr->n-1];
    int dum[size];
    for(l=0;l<size;l++){
      dum[l] = tr->atoms[tr->n-1][l];
    }
    for(l=0;l<size;l++){
      for(m=0;m<al->nbonds[dum[l]];m++){
        add_to_group(tr,tr->n-1,al->bonds[dum[l]][m]);
      }
    }
    
  }  

}
/*=============================================================*/
bool is_in_five_ring(t_atomlist *al, int i)

/* check whether atom i is part of a 5-ring */

{
  int k,j,l;
  int at1,at2;
  
  bool check = FALSE;
  for(k=0;k<al->nb13[i];k++){
    at1 = al->b13[i][k];
    for(j=0;j<al->nbonds[at1];j++){
      at2 = al->bonds[at1][j];
      for(l=0;l<al->nb13[i];l++){
        if(al->b13[i][l] == at2) check = TRUE;
      }
    }
  }
  return check;
}

/*=============================================================*/
bool is_in_six_ring(t_atomlist *al, int i)

/* check whether atom i is part of a 6-ring */

{
  int k,j,l;
  int at1,at2;
  int count = 0;
  bool check = FALSE;
  for(k=0;k<al->nb14[i];k++){
    at1 = al->b14[i][k];
    count = 0;
    for(j=0;j<al->nbonds[at1];j++){
      for(l=0;l<al->nb13[i];l++){
        if(al->bonds[at1][j] == al->b13[i][l]) {
          count++;
        }
      }
    }
    if(count == 2) {
      check = TRUE;
    }
  }
  return check;
}




/*============================================================*/
void hbonds(t_atomlist *al, t_idxgroups *don, t_idxgroups *acc)

/* check who is a hydrogen bond donor or acceptor */
{
  int i;
  for(i=0;i<al->natoms;i++){
    al->isdon[i] = al->isacc[i] = FALSE;
    if(strcmp(al->symbol[i],"H") ==0){
      if(strcmp(al->symbol[al->bonds[i][0]],"N") ==  0 ||
         strcmp(al->symbol[al->bonds[i][0]],"O") ==  0){
        add_to_group(don,0,i); 
        al->isdon[i] = TRUE;
      }
    }
    else if(strcmp(al->symbol[i],"O") == 0){
      add_to_group(acc,0,i);  
      al->isacc[i] = TRUE;
    }
    else if(strcmp(al->symbol[i],"N") == 0){
      if(strcmp(al->hyb[i],"sp2") == 0 &&
         al->nbonds[i] == 2){
        add_to_group(acc,0,i);
        al->isacc[i] = TRUE; 
      }
    }
    
  }
}
/*============================================================*/

void fake_charge(t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->resname[i],"ARG") == 0 &&
       (strcmp(al->name[i]," NH1") == 0 ||
        strcmp(al->name[i]," NH2") == 0)) al->q[i] = 1.;
    else if(strcmp(al->resname[i],"LYS") == 0 &&
            strcmp(al->name[i]," NZ ") == 0) al->q[i] = 1.;
    else if(strcmp(al->resname[i],"GLU") == 0 &&
       (strcmp(al->name[i]," OE1") == 0 ||
        strcmp(al->name[i]," OE2") == 0)) al->q[i] = -1.;
    else if(strcmp(al->resname[i],"ASP") == 0 &&
       (strcmp(al->name[i]," OD1") == 0 ||
        strcmp(al->name[i]," OD2") == 0)) al->q[i] = -1.;
  }
}
/*============================================================*/

void hydrophobics(t_atomlist *al, t_idxgroups *phob, FILE *fp)
/* assign hydrophobic flags to atoms */
{
  int i,k,j;

/* if no file is given we try to estimate */
  if (fp==NULL){
    for(i=0;i<al->natoms;i++){
      bool check = FALSE;
      al->ishphob[i] = FALSE;
      if(strcmp(al->symbol[i],"C") == 0){
        for(k=0;k<al->nbonds[i];k++){
          if((strcmp(al->symbol[al->bonds[i][k]],"N") == 0) ||
             (strcmp(al->symbol[al->bonds[i][k]],"O") == 0)) {
            check=TRUE;
            break;
          }
        }
      for(k=0;k<al->nb13[i];k++){
        if((strcmp(al->symbol[al->b13[i][k]],"N") == 0) ||
           (strcmp(al->symbol[al->b13[i][k]],"O") == 0)) {
          check=TRUE;
          break;
        }
      }
        if(!check){
          add_to_group(phob,0,i);
          al->ishphob[i] = TRUE;
        }
      }
    }
  }

  /* read definitions from file */

  else{
    t_namegroups *ng = read_namegroups(fp);
    for(i=0;i<al->natoms;i++){
      al->ishphob[i] = FALSE;
      for(k=0;k<ng->n;k++){
        if(strcmp(al->resname[i],ng->resname[k]) == 0){
          for(j=0;j<ng->natoms[k];j++){
            if(strcmp(al->name[i],ng->atomnames[k][j]) == 0){
              add_to_group(phob,0,i);
              al->ishphob[i] = TRUE;
            }
          }
        }
      }
    }
  }
}


/*============================================================*/

bool protected(t_atomlist *al, int i, int minHP, real d)

/* check the number of hydrophobic atoms around 
   a hydrogen bond */

{
  int k;
  int count = 0;
  for(k=0;k<al->nnb[i];k++){
    if(al->ishphob[al->nb[i][k]]){
      if(DIST(al,i,al->nb[i][k]) < d) count++;
    }
  }
  if(count >= minHP) return TRUE;
  else return FALSE;
}

/*============================================================*/
bool protected2(t_atomlist *al, int i, int at1,int minHP, real rad)

/* check the number of hydrophobic atoms around 
   a hydrogen bond */

{
  int k,j,at,at2,l;
  real d,d2,d3;
  int count;
  int n;
  int c;
  n = al->bonds[i][0];
  c = al->bonds[at1][0];
  
  count = 0;
  for(k=0;k<al->nnb[n];k++){
    j = al->nb[n][k];
    if(al->rs_type[j] == 3 ||
       al->rs_type[j] == 4 ||
       al->rs_type[j] == 5 ||
       al->rs_type[j] == 6 ||
       al->rs_type[j] == 7){
      d2 = dist_ij(al->x[n],al->x[j]);
      if(d2 < rad){
        d3 = dist_ij(al->x[c],al->x[j]);
        if(d3 < rad) count++;
      }
    }
  }
/*   printf("count = %d\n",count); */
  
  if(count > minHP) return TRUE;
  else return FALSE;
}
/*============================================================*/
bool protected3(t_atomlist *al, int i, real minHP, real rad)

/* check the number of hydrophobic atoms around 
   a hydrogen bond */

{
  int k,j;
  real d;
  real xx = 0.;
  int count = 0;
  real weight;
  
  for(k=0;k<al->nnb[i];k++){
    j = al->nb[i][k];
    if(al->rs_type[j] == 3 ||
       al->rs_type[j] == 4 ||
       al->rs_type[j] == 5 ||
       al->rs_type[j] == 6 ||
       al->rs_type[j] == 7){
      d = dist_ij(al->x[i],al->x[j]);
      if(d < rad){
        if(al->ishphob[j]) weight = 1.;
        else if(IS_CHARGED_RESI(al,j)){
          weight = 1./(2*al->order[j]);
        }
        
        else {
          weight = 1./(1*al->order[j]);
        }
/*         printf("%s / %s weight = %g\n",al->name[j],al->resname[j],weight); */
        
        count++;
        xx+=1./sqr(d)*weight;
      }
    }
  }
/*   printf("name: %s (%d/%s) xx = %g (%d)\n",al->name[i],al->resid[i],al->resname[i],xx,count); */
  
  if(xx > minHP) return TRUE;
  else return FALSE;
}


/*============================================================*/
bool protected4(t_atomlist *al, int i, real minHP, real rad)
{
  int k;
  real d;
  real solv = 0;
  int nat = 0;
  for(k=0;k<al->nnb[i];k++){
    d = DIST2(al,i,al->nb[i][k]);
    if(d < rad*rad){
      solv+=al->cnc_solv[al->nb[i][k]];
      nat++;
    }
  }
  for(k=0;k<al->nb14[i];k++){
    d = DIST2(al,i,al->b14[i][k]);
    if(d < rad*rad){
      solv+=al->cnc_solv[al->b14[i][k]];
      nat++;
    }
  }
  solv/=nat;
/*   solv*=10.; */
  if(solv < minHP) return TRUE;
  else return FALSE;
}

/*============================================================*/
void do_cov(t_atomlist *al, t_bondlist *bl, t_contab *rct, t_excl *ex)
{
  /* make covalent connections in contab */

  int i,k;
  int at1,at2;
  int res1,res2;
  int flag1,flag2;
  
  /* find interresidual bonds */
  
  for(i=0;i<bl->n;i++){
    at1 = bl->at1[i];
    at2 = bl->at2[i];
    if(al->resid[at1] != al->resid[at2] &&
       !is_excluded(ex,al,at1,at2)){
      res1 = al->resid[at1]-1;
      res2 = al->resid[at2]-1;
      if(strcmp(al->name[at1]," SG ") == 0 &&
         strcmp(al->name[at2]," SG ") == 0) {
        flag1 = flag2 = eSUL;
      }
      else flag1 = flag2 = eCOV;
      add2contab(rct,res1,res2,flag1,flag2,0);
    }
  }
}

  
    
/*============================================================*/
void do_hbonds(t_atomlist *al, t_contab *rct,
               FILE *fp,bool bCheck[], 
               real hbdist,real hbangle, 
               rvec minHP, real *mindist,
               int *bbhb, int *bbhb_rej,
               int *scbbhb, int *scbbhb_rej,
               int *scschb, int *scschb_rej)

/* search hydrogen bonds and estimate stability if recommanded */

{
  int i,k;
  int at1,at2;
  int flag1 = eNON;
  int flag2 = eNON;
  real d, a;
  *bbhb = 0;
  *bbhb_rej  = 0;
  *scbbhb = 0;
  *scbbhb_rej = 0;
  *scschb = 0;
  *scschb_rej = 0;
  
  for(i=0;i<al->natoms;i++){
    if(al->isdon[i]){
      for(k=0;k<al->nnb[i];k++){
        at1 = al->nb[i][k];
        if(al->isacc[at1] && 
           al->resid[i]!=al->resid[at1]){
          at2 = al->bonds[i][0];
          d = DIST(al,i,at1);
          a = RAD2DEG*angle_ij_ik(al->x[i],al->x[at1],al->x[at2]);
          if(d < hbdist &&
             a > hbangle){
            /* is it a backbone-backbone bond ? */
            if(al->order[i] <= 1 &&
               al->order[at1] <= 1) {
              add_hbond(al,i,at1,ehBBBB,mindist[0]);
              flag1 = eBBHD;
              flag2 = eBBHA;
              if(bCheck[0] && (al->restype[i] == rtPROTEIN &&
                               al->restype[at1] == rtPROTEIN)) 
              {
                if(al->hbonds[al->nhbonds-1]->prot > minHP[0]) {
                  al->hbonds[al->nhbonds-1]->constr=FALSE;
                  flag1=eNPR;
                  flag2=eNPR;
                  *bbhb_rej+=1;
                }
                else *bbhb+=1;
              }
              else *bbhb+=1;
            }
            /* or backbone-sidechain */
            else if((al->order[i] > 1 &&
                     al->order[at1] <= 1) ||
                    (al->order[i] <= 1 &&
                     al->order[at1] > 1)) {
              if(al->order[i] <= 1){
                flag1 = eSBH_BD;
                flag2 = eSBH_SA;
                add_hbond(al,i,at1,ehBBSC,mindist[1]);
              }
              else{
                flag1 = eSBH_SD;
                flag2 = eSBH_BA;
                add_hbond(al,i,at1,ehSCBB,mindist[1]);
              }             
              
              if(bCheck[1] && ( al->restype[i] == rtPROTEIN &&
                                al->restype[at1] == rtPROTEIN))
              {
                
                if(al->hbonds[al->nhbonds-1]->prot > minHP[1]) {
                  al->hbonds[al->nhbonds-1]->constr=FALSE;
                  flag1=eNPR;
                  flag2=eNPR;
                  *scbbhb_rej+=1;
                }
                else *scbbhb+=1;
              }
              else *scbbhb+=1;
            }
            /* or sidechain-sidechain ? */
            else if(al->order[i] > 1 &&
                    al->order[at1] > 1) {
              add_hbond(al,i,at1,ehSCSC,mindist[2]);
              flag1 = eSSHD;
              flag2 = eSSHA;
              if(bCheck[2] && (al->restype[i] == rtPROTEIN &&
                               al->restype[at1] == rtPROTEIN)) 
              {
                if(al->hbonds[al->nhbonds-1]->prot > minHP[2]) {
                  al->hbonds[al->nhbonds-1]->constr=FALSE;
                  flag1=eNPR;
                  flag2=eNPR;
                  *scschb_rej+=1;
                }
                else *scschb+=1;
              }
              else *scschb+=1;
            }
            add2contab(rct,al->resid[i]-1,al->resid[at1]-1,flag1,flag2,0);
            if(fp){
              fprintf(fp,"%5d%4s|%4s ==>> %5d%4s|%4s  %4d  %4d \n",al->resid[i],
                      al->resname[i],al->type[i],al->resid[at1],al->resname[at1],
                      al->type[at1],flag1,flag2);
            }
          }
        }
      }
    }
  }
}
/*=============================================================*/
void do_hydrophobics(t_atomlist *al,t_contab *rct,real min,
                     t_idxgroups *phb)
{
  /* a hydrophobic cluster is defined if
     three hydrophobic atoms from three
     different residues are found within
     min A
  */
  int i,k,j;
  int at1,at2;
  int count=0;
  real d;
  phb = idx_realloc(phb,1);
  
  for(i=0;i<al->natoms;i++){
    if(al->ishphob[i]){
      for(k=0;k<al->nnb[i];k++){
        at1 = al->nb[i][k];
        if(al->resid[i] != al->resid[at1] && al->ishphob[at1]){
          if(DIST(al,i,at1) < min){
            for(j=0;j<al->nnb[at1];j++){
              at2=al->nb[at1][j];
              if(al->resid[at2] != al->resid[i] && al->resid[at2]!=al->resid[at1] &&
                 al->ishphob[at2]){
                if(DIST(al,i,at2) < min &&
                   DIST(al,at1,at2) < min){
/*                   add2contab(rct,al->resid[i]-1,al->resid[at1]-1,ePHO,ePHO,0); */
/*                   add2contab(rct,al->resid[i]-1,al->resid[at2]-1,ePHO,ePHO,0); */
/*                   add2contab(rct,al->resid[at1]-1,al->resid[at2]-1,ePHO,ePHO,0); */
                  phb = idx_realloc(phb,count+1);
                  phb->n++;
/*                   add_to_group(phb,count,i); */
/*                   add_to_group(phb,count,at1); */
/*                   add_to_group(phb,count,at2); */
                  add_to_group(phb,count,al->resid[i]-1);
                  add_to_group(phb,count,al->resid[at1]-1);
                  add_to_group(phb,count,al->resid[at2]-1);
                  count++;
                  
/*                   printf("found cluster %d%s--%d%s--%d%s\n",al->resid[i], */
/*                          al->resname[i],al->resid[at1],al->resname[at1], */
/*                          al->resid[at2],al->resname[at2]); */
                }
              }
            }
          }
        }
      }
    }
  }
  phb = remove_redundant_groups(phb);

}
/*=============================================================*/
void print_hydrophobics_to_file(FILE *fp, t_atomlist *al, t_idxgroups *phb)
{
  int i;
  fprintf(fp,"[ HYDROPHOBIC CLUSTERS ]\n");
  for(i=0;i<phb->n;i++){
    fprintf(fp,"%8d %8d %8d\n",phb->atoms[i][0]+1,
            phb->atoms[i][1]+1,
            phb->atoms[i][2]+1);
  }
}



/*=============================================================*/


void get_hybrid(t_atomlist *al)

/* get hybridisation (only sp2 or sp3 is of interest */

{
  int i,k,j,l;
  int at,at2,at3;
  real d;
  int count = 0;
  
  for(i=0;i<al->natoms;i++){
    if(al->nbonds[i]>2 && (strcmp(al->symbol[i],"C") == 0 ||
                            strcmp(al->symbol[i],"N") == 0)){
      t_idxgroups *tmp = idx_init();
      tmp = idx_realloc(tmp,1);
      add_to_group(tmp,0,i);
      for(k=0;k<al->nbonds[i];k++){
        add_to_group(tmp,0,al->bonds[i][k]);
      }
      real tol = 0.005;
      real plan;
      if(is_planar_group(al,tmp->natoms[0],tmp->atoms[0],tol,&plan)){
        strcpy(al->hyb[i],"sp2");           
      }
      else{
        strcpy(al->hyb[i],"sp3");
      }
    }
  }
}
/*============================================================*/
void get_hybrid2(t_atomlist *al, t_types *tp)
{
  int i,k,j,l;
  int at,at2,at3;
  real d;
  int count = 0;
  
  for(i=0;i<al->natoms;i++){
    /* check carbon atoms */
    if(strcmp(al->hyb[i],"") == 0){
      if(strcmp(al->symbol[i],"C") == 0)
      {
        bool has_double = FALSE;
        for(k=0;k<al->nbonds[i];k++){
          at = al->bonds[i][k];
          d = DIST(al,i,at);
          
          if(strcmp(al->symbol[at],"C") == 0)
          {
            /* check for double/aromatic bonds */
            if(d < 1.46){
              has_double = TRUE;
            }
          }
          else if(strcmp(al->symbol[at],"N") == 0)
          {
            /* check for omega bond */
            if(d < 1.42) has_double = TRUE;
          }
          else if(strcmp(al->symbol[at],"O") == 0)
          {
            /* check for C=O */
            if(d < 1.3) has_double = TRUE;
          }
        }
        if(has_double)
          strcpy(al->hyb[i],"sp2");
        else strcpy(al->hyb[i],"sp3");
        if (al->nbonds[i] == 4) strcpy(al->hyb[i],"sp3");
        
      }
      
      /* now check nitrogen atoms */
      
      else if(strcmp(al->symbol[i],"N") == 0)
      {
/*         printf("checking...... %d %s %d %s\n",al->id[i],al->name[i],al->resid[i],al->resname[i]); */
        
        bool has_double = FALSE;
        for(k=0;k<al->nbonds[i];k++){
          at = al->bonds[i][k];
          d = DIST(al,i,at);
          if(strcmp(al->symbol[at],"C") == 0)
          {
            if(d < 1.42) has_double = TRUE;
          }
        }
        if(has_double)
          strcpy(al->hyb[i],"sp2");
        else strcpy(al->hyb[i],"sp3");
      }
    }
  }
}

/*============================================================*/

real radius_of_gyration(t_atomlist *al, rvec gvec)
{
  int i,k,m;
  matrix mat;
  rvec dd;
  rvec comp;
  real dx2;
  real gyro;
  real tm = 0.;
  rvec x[al->natoms];
  copy_coords(al->x,x,al->natoms);
  
  clear_rvec(dd);
  clear_rvec(comp);
  clear_mat(mat);
  com(al);
  princ_comp(al->natoms,al->x,mat,dd);
  
  if (det(mat) < 0) {
    for(m=0; (m<DIM); m++)
      mat[ZZ][m] = -mat[ZZ][m];
  }
  rotate_atoms(al->natoms,al->x,mat);
  
  for(i=0;i<al->natoms;i++){
    tm+=al->m[i];
    for(m=0;m<DIM;m++){
      dx2 = sqr(al->x[i][m]);
      comp[m]+=dx2*al->m[i];
    }
  }
  gyro = comp[XX]+comp[YY]+comp[ZZ];
  for(m=0; (m<DIM); m++)
    gvec[m]=sqrt((gyro-comp[m])/tm);

  copy_coords(x,al->x,al->natoms);
  return sqrt(gyro/tm);
}

/*============================================================*/

void sequence(t_resl *rl, char *seq)
{
  FILE *fp = cnclib("aa.dat");
  char line[STRLEN];
  int n = 0;
  int i,k;
  char rdum[STRLEN];
  shstr three[1000];
  shstr one[1000];
  bool bFound = FALSE;
/*   char xx[STRLEN]; */
  
  while(get_a_line(fp,line,STRLEN)){
    sscanf(line,"%s %s",three[n],one[n]);
    n++;
  }
    
  for(i=0;i<rl->nres;i++){
    strcpy(rdum,rl->resname[i]);
    trim(rdum);
    bFound = FALSE;
    for(k=0;k<n;k++){
      if(strcmp(rdum,three[k]) == 0){
        seq[i] = one[k][0];
        bFound = TRUE;
      }
    }
    if(!bFound)
    {
      seq[i] = 'X';
    }
  }
  seq[rl->nres] = 0;
}

void print_sequence_log(FILE *fp, char *seq, int n)
{
  int i;
  fprintf(fp,"tCNC__seq_> ");
  for(i=0;i<n;i++){
    fprintf(fp,"%c",seq[i]);
    if( (i+1) % 70 == 0) fprintf(fp,"\ntCNC__seq_> ");
  }
  fprintf(fp,"\n");
  fflush(fp);
  
}

