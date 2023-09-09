#include <tconcoord.h>
/*=================================================================*/
t_hbond *hbond_init(void)
{
  t_hbond *hb = NULL;
  snew(hb,1);
  hb->don = 0;
  hb->acc = 0;
  hb->acc_b = 0;
  hb->don_b = 0;
  hb->type = eNON;
  hb->d_a_dist = 0.;
  hb->d_ab_dist = 0.;
  hb->db_a_dist = 0;
  hb->db_ab_dist = 0;
  hb->angle = 0.;
  hb->angle = 0.;
  hb->dihed = 0.;
  hb->energy = 0.;
  hb->prot = 0;
  hb->constr = TRUE;
  hb->packing = 0.;
  return hb;
}

/*=================================================================*/

t_hbond *copy_hbond(t_hbond *hb)
{
  t_hbond *new = hbond_init();
  new->don = hb->don;
  new->acc = hb->acc;
  new->don_b = hb->don_b;
  new->acc_b = hb->acc_b;
  new->type = hb->type;
  new->d_a_dist = hb->d_a_dist;
  new->d_ab_dist = hb->d_ab_dist;
  new->db_a_dist = hb->db_a_dist;
  new->db_ab_dist = hb->db_ab_dist;
  new->angle = hb->angle;
  new->angle2 = hb->angle2;
  new->dihed = hb->dihed;
  new->energy = hb->energy;
  new->prot = hb->prot;
  new->constr = hb->constr;
  new->packing = hb->packing;
  return new;
}

void free_hbonds(t_atomlist *al)
{
  int i;
  for(i=0;i<al->nhbonds;i++){
    sfree(al->hbonds[i]);
  }
  al->nhbonds = 0;
  sfree(al->hbonds);
  snew(al->hbonds,1);
}



/*=================================================================*/

void add_hbond(t_atomlist *al, int don, int acc, int type,
               real pr_rad)
{
  int i,n;
  real d,a;
  int acc_b, don_b;
  char err[STRLEN];
  
  acc_b = don_b = -1;
  
  n = al->nhbonds;
  al->nhbonds+=1;
  srenew(al->hbonds,al->nhbonds);
  al->hbonds[n] = hbond_init();
  
  /* find don bond */
  don_b = al->bonds[don][0];
  
  /* find acc bond */

  if(al->nbonds[acc] == 1) acc_b = al->bonds[acc][0];
  else {
    for(i=0;i<al->nbonds[acc];i++){
      if(strcmp(al->symbol[al->bonds[acc][i]],"H") != 0){
        acc_b = al->bonds[acc][i];
        break;
      }
    }
  }
  if(don_b == -1) {
    sprintf(err,"Hbond donor (%d%s/%d%s) has no bond\n",
            al->id[don],al->name[don],al->resid[don],al->resname[don]);
    CNCerr(stderr,err);
  }

  if(acc_b == -1) {
    sprintf(err,"Hbond acceptor (%d%s/%d%s) has no bond\n",
            al->id[acc],al->name[acc],al->resid[acc],al->resname[acc]);
    CNCerr(stderr,err);
  }
  

  al->hbonds[n]->don = don;
  al->hbonds[n]->acc = acc;
  al->hbonds[n]->don_b = don_b;
  al->hbonds[n]->acc_b = acc_b;
  
  al->hbonds[n]->type = type;
  

  /* calc geometry */

  al->hbonds[n]->d_a_dist = DIST(al,don,acc);
  al->hbonds[n]->d_ab_dist = DIST(al,don,acc_b);
  al->hbonds[n]->db_a_dist = DIST(al,don_b,acc);
  al->hbonds[n]->db_ab_dist = DIST(al,don_b,acc_b);
  
  al->hbonds[n]->angle = RAD2DEG*angle_ij_ik(al->x[don],al->x[don_b],al->x[acc]);
  al->hbonds[n]->angle2 = RAD2DEG*angle_ij_ik(al->x[acc],al->x[don],al->x[acc_b]);
  al->hbonds[n]->dihed = RAD2DEG*dihedral(al->x[don_b],al->x[don],
                              al->x[acc],al->x[acc_b]);

  al->hbonds[n]->energy = hbond_energy(al->hbonds[n]->d_a_dist);
  al->hbonds[n]->prot = hbond_protection(al,don,acc,pr_rad);
  
}
/*=================================================================*/
void get_hbonds(t_atomlist *al, real dmax, real amin, real pr_rad)
{
  int i,k;
  real d,a;
  int at,don,acc,don_b;
  int type;
  t_idxgroups *didx = idx_init();
  t_idxgroups *aidx = idx_init();
  didx = idx_realloc(didx,1);
  aidx = idx_realloc(aidx,1);
  
  hbonds(al,didx,aidx);
  free_groups(didx); 
  free_groups(aidx); 
    
  for(i=0;i<al->natoms;i++){
    if(al->isdon[i]){
      don = i;
      for(k=0;k<al->nnb[don];k++){
        at = al->nb[don][k];
        if(al->isacc[at]){
          acc = at;
          don_b = al->bonds[don][0];
          d = DIST(al,don,acc);
          if (d < dmax) {
            a = RAD2DEG*angle_ij_ik(al->x[don],al->x[acc],al->x[don_b]);
            if (a > amin) {
              
              if(al->order[don] == 1 &&
                 al->order[acc] == 1) {
                type = ehBBBB;
                add_hbond(al,don,acc,type,pr_rad);
              }
              else if(al->order[don] != 1 &&
                      al->order[acc] == 1){
                type = ehSCBB;
                add_hbond(al,don,acc,type,pr_rad);
              }
              else if(al->order[don] == 1 &&
                      al->order[acc] != 1){
                type = ehBBSC;
                add_hbond(al,don,acc,type,pr_rad);
              }
              else if(al->order[don] != 1 &&
                      al->order[acc] != 1){
                type = ehSCSC;
                add_hbond(al,don,acc,type,pr_rad);
              }
            }
          }
        }
      }
    }
  }
}
/*=================================================================*/
real hbond_protection(t_atomlist *al, int don, int acc, real rad)
{
  int k;
  real d;
  real solv = 0;
  int nat = 0;
  real rad2 = sqr(rad);
  t_idxgroups *xx = idx_init();
  xx = idx_realloc(xx,1);
  
      

  for(k=0;k<al->nnb[don];k++){
    d = DIST2(al,don,al->nb[don][k]);
    if(d < rad2){
/*       printf("here\n"); */
      
      add_to_group(xx,0,al->nb[don][k]);
      
/*       solv+=al->cnc_solv[al->nb[don][k]]; */
/*       nat++; */
    }
  }

  for(k=0;k<al->nb14[don];k++){
    d = DIST2(al,don,al->b14[don][k]);
    if(d < rad2){
      add_to_group(xx,0,al->b14[don][k]);
/*       solv+=al->cnc_solv[al->b14[don][k]]; */
/*       nat++; */
    }
  }

  for(k=0;k<al->nnb[acc];k++){
    d = DIST2(al,acc,al->nb[acc][k]);
    if(d < rad2){
      add_to_group(xx,0,al->nb[acc][k]);
      
/*       solv+=al->cnc_solv[al->nb[acc][k]]; */
/*       nat++; */
    }
  }

  for(k=0;k<al->nb14[acc];k++){
    d = DIST2(al,acc,al->b14[acc][k]);
    if(d < rad2){
      add_to_group(xx,0,al->b14[acc][k]);
/*       solv+=al->cnc_solv[al->b14[acc][k]]; */
/*       nat++; */
    }
  }
/*   printf("natoms = %d\n",xx->natoms[0]); */
  
  for(k=0;k<xx->natoms[0];k++){
    solv+=al->cnc_solv[xx->atoms[0][k]];
  }
  solv/=xx->natoms[0];
/*   printf("solv = %g\n",solv);  */
  
/*   print_group(al,xx->atoms[0],xx->natoms[0]); */
  
  free_groups(xx);
  





/*   solv/=nat; */
  return solv;
}
          
/*=================================================================*/
real hbond_energy(real d)
{
  real energy;
  
  energy = -25.36*1000*exp(-3.6*d);
  return energy;
  
}
/*=================================================================*/
void print_hbonds_to_file(FILE *fp, t_atomlist *al)
{
  int i;
  int don,acc,don_b,acc_b;
  real dda,ddba,ddbab,ddab,a1,a2,dih,en,prot,packing;
  bool constr;
  shstr con;
  
  fprintf(fp,";==================================================\n");
  fprintf(fp,"; tCONCOORD HBOND FILE\n");
  fprintf(fp,";==================================================\n");
  fprintf(fp,"[ HBONDS ]\n");
  fprintf(fp,";%16s %16s %16s         %8s %8s %8s %8s  %8s %8s %8s %8s %8s %8s %8s\n",
          "DONOR","ACCEPTOR","DONOR-BOND","ACC-BOND","d(D-A)","d(DB-A)","d(D-AB)","d(DB-AB)","a(DB-D-A)","a(D-A-AB)","dih(DB-D-A-DB)","E(kJ/mol)","PS","CONSTR");
  
  
    
  for(i=0;i<al->nhbonds;i++){
    don = al->hbonds[i]->don;
    acc = al->hbonds[i]->acc;
    don_b = al->hbonds[i]->don_b;
    acc_b = al->hbonds[i]->acc_b;
    dda = al->hbonds[i]->d_a_dist;
    ddba = al->hbonds[i]->db_a_dist;
    ddbab = al->hbonds[i]->db_ab_dist;
    ddab = al->hbonds[i]->d_ab_dist;
    a1 = al->hbonds[i]->angle;
    a2 = al->hbonds[i]->angle2;
    dih = al->hbonds[i]->dihed;
    en = al->hbonds[i]->energy;
    prot = al->hbonds[i]->prot;
    constr = al->hbonds[i]->constr;
    packing = al->hbonds[i]->packing;
    if(constr) strcpy(con,"yes");
    else strcpy(con,"no");
    
/*     fprintf(fp,"%6d %4s %5d %4s -> %6d %4s %5d %4s ||  %6d %4s %5d %4s ->  %6d %4s %5d %4s  || %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %8.5f %3s\n", al->id[don],al->name[don],al->resid[don],al->resname[don], */
/*             al->id[acc],al->name[acc],al->resid[acc],al->resname[acc], */
/*             al->id[don_b],al->name[don_b],al->resid[don_b],al->resname[don_b], */
/*             al->id[acc_b],al->name[acc_b],al->resid[acc_b],al->resname[acc_b], */
/*             dda,ddba,ddab,ddbab,a1,a2,dih,en,prot,packing,con); */
    fprintf(fp,"%6d %4s %5d %4s -> %6d %4s %5d %4s ||  %6d %4s %5d %4s ->  %6d %4s %5d %4s  || %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %3s\n", al->id[don],al->name[don],al->resid[don],al->resname[don],
            al->id[acc],al->name[acc],al->resid[acc],al->resname[acc],
            al->id[don_b],al->name[don_b],al->resid[don_b],al->resname[don_b],
            al->id[acc_b],al->name[acc_b],al->resid[acc_b],al->resname[acc_b],
            dda,ddba,ddab,ddbab,a1,a2,dih,en,prot,con);

/*     fprintf(fp,"%20s %20s %20s %20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f\n", */
/*             donstr,accstr,don_bstr,acc_bstr, */
/*             dda,ddba,ddab,ddbab,a1,a2,dih,en,prot); */
  }
}

/*=================================================================*/

real energy_from_histogram(t_histogram *h, real x)
{

  /* do linear spline interpolation */
  real diff;
  int lower, upper;
  real y_lower, y_upper, x_lower, x_upper;
  real val;

  diff = x - h->x[0];
/*   printf("x = %g diff = %g\n",x,diff); */
  lower = (int) (diff*h->inv_delta);
  upper = lower+1;
  y_lower = h->y[lower];
  x_lower = h->x[lower];
  y_upper = h->y[upper];
  x_upper = h->x[upper];
/*   printf("x_lower = %g x_upper = %g y_lower = %g y_upper = %g\n", */
/*          x_lower, x_upper, y_lower, y_upper); */
  
  val = y_lower + (y_upper - y_lower)/(x_upper-x_lower)*(x-x_lower);
  return val;
}

