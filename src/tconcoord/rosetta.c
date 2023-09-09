#include <tconcoord.h>
/*==========================================================*/

void get_rosetta_types(t_atomlist *al)
{
  FILE *fp = cnclib("rosetta.dat");
  
  int i = 0;
  int k;
  t_types *tp = types_init();
  char line[STRLEN];
  char dummy[STRLEN];

  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      substring(dummy,line,'[',']');
      sscanf(dummy+1,"%s",dummy);
    }
    else{
      cut_string(line,'#');
      sscanf(line,"%s %d",tp->name[i],&tp->rs_type[i]);
      strcpy(tp->resname[i],dummy);
      extend_name(tp->name[i]);
      i++;
    }
  }
  tp->n = i;
  fclose(fp);

  /* assign rs_types to atomlist */

  /* now do the defaults */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(tp->resname[k],"DEFAULT")==0){
        if(strcmp(al->name[i],tp->name[k])==0){
          al->rs_type[i] = tp->rs_type[k];
        }
      }
    }
  }

  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(al->name[i],tp->name[k]) == 0 &&
         strcmp(al->resname[i],tp->resname[k]) == 0){
        al->rs_type[i] = tp->rs_type[k];
      }
    }
  }

  free_types(tp);
  
  fp = cnclib("rosetta_params.dat");
  
  t_vdw *vdw = vdw_init();
  
  i = 0;
  if(!find_key_word(fp,"[ PARAMETERS ]")){
    fatal_error("Shitty file!\n");
  }
  
  while(get_a_line(fp,line,STRLEN) &&
        strchr(line,'[')==NULL)
  {
#ifdef GMX_DOUBLE
    sscanf(line,"%d %lf %lf %lf %lf %lf %lf",&vdw->rs_type[i],
           &vdw->rs_rad[i],&vdw->rs_eps[i],&vdw->rs_G[i],
           &vdw->rs_V[i],&vdw->rs_lambda[i], &vdw->rs_Gref[i]);
#else
    sscanf(line,"%d %f %f %f %f %f %f",&vdw->rs_type[i],
           &vdw->rs_rad[i],&vdw->rs_eps[i],&vdw->rs_G[i],
           &vdw->rs_V[i],&vdw->rs_lambda[i],&vdw->rs_Gref[i]);
#endif
    i++;
  }
  vdw->n = i;

/* assign values to atomlist */

  for(i=0;i<al->natoms;i++){
    for(k=0;k<vdw->n;k++){
      if(al->rs_type[i] == vdw->rs_type[k]){
        al->rs_rad[i] = vdw->rs_rad[k];
        al->rs_eps[i] = vdw->rs_eps[k];
        al->rs_G[i] = vdw->rs_G[k];
        al->rs_Gref[i] = vdw->rs_Gref[k];
        al->rs_V[i] = vdw->rs_V[k];
        al->rs_lambda[i] = vdw->rs_lambda[k];
      }
    }
  }



  t_vdwcomb *comb = read_vdwcomb2(fp, 2);
  al->nrs_comb = vdw->n+1;
  snew(al->rs_comb,al->nrs_comb);
  for(i=0;i<al->nrs_comb;i++){
    snew(al->rs_comb[i],al->nrs_comb);
  }
  for(i=1;i<al->nrs_comb;i++){
    for(k=1;k<al->nrs_comb;k++){
      al->rs_comb[i][k] = vdw->rs_rad[i-1] + vdw->rs_rad[k-1]; 
    }
  }
  for(i=0;i<comb->n;i++){
    al->rs_comb[comb->rs_type1[i]][comb->rs_type2[i]] = comb->rs_comb[i];
    al->rs_comb[comb->rs_type2[i]][comb->rs_type1[i]] = comb->rs_comb[i];
  }

  free_vdw(vdw);  
  fclose(fp);  
}

/*====================================================================*/
/*     ROSETTA ENERGY FUNCTIONS                                       */
/*====================================================================*/

real rs_lj_energy(t_atomlist *al)
{
  
  /* ROSETTA lennard-jones energy */

  int i,k,j;
  int nb;
  real en,dij,rij,eij,r;
  rvec diff;
  real slope_const = -400.59961665224989; /* -12(1.33^13 - 1.33^7) */
  real slope;
  real y_intercept;
  real y_const = 19.56532575351482;  /* 1.33^12 - 2(1.33)^6 */

  
  en = 0;
  real atr = 0.;
  real rep = 0.;

  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      j = al->nb[i][k];
      if(i > j)
      {

        rvec_sub(al->x[i],al->x[j],diff);
        dij = norm(diff);
        rij = al->rs_comb[al->rs_type[i]][al->rs_type[j]]; 
        eij = sqrt(al->rs_eps[i]*al->rs_eps[j]);

        r = rij / dij;
        if(r < 1.12) atr += (pow(r,12) - 2*pow(r,6)) * eij;
        else{
/*           if(r <= 1.33) { */
            
            rep += (pow(r,12) - 2*pow(r,6)) * eij;
/*           } */
/*           else { */
/*             slope = slope_const*eij/rij;  */
/*             y_intercept = -slope*rij/1.33+eij*y_const;  */
/*             rep += y_intercept - dij * slope;  */
/*           }  */
        }
      }
    }
  }
  
/*   printf("atr-lj = %g\n",atr);  */
/*   printf("rep-lj = %g\n",rep);  */
/*   printf("atr+rep lj = %g\n",atr+rep);  */
  return (atr+rep);
}
/*==================================================================*/

real lj_energy_ij(real dij,real rij, real eij)
{
  real r;
  r = rij / dij;
  return (pow(r,12) - 2*pow(r,6)) * eij;
}

  
/*==================================================================*/
real lj_lattice_energy(t_atomlist *al,rvec x,real cutoff, 
                       real eps, real r, int maxorder)
{
  /* calculate lj potential at point x */

  int i,k;
  real dij,rij,eij,cut2;
  rvec diff;
  real en = 0;
  cut2 = sqr(cutoff);
    
  for(i=0;i<al->natoms;i++){
    if(al->rs_rad[i] > 0){
      if(sqr(al->x[i][XX]-x[XX]) < cut2 &&
         sqr(al->x[i][YY]-x[YY]) < cut2 &&
         sqr(al->x[i][ZZ]-x[ZZ]) < cut2){
        rvec_sub(al->x[i],x,diff);
        dij = norm(diff);
        if(dij < cutoff){
          eij = sqrt(al->rs_eps[i]*eps);
          rij = al->rs_rad[i]+r;
          if(maxorder == 0 || al->order[i] <=maxorder)
            en+=lj_energy_ij(dij,rij,eij);
        }
      }
    }
  }
  return en;
}

  


/*==================================================================*/

real rs_solv_energy(t_atomlist *al)
{
    /* ROSETTA solvation energy (Lazaridis-Karplus solvation model) */

  int i,k,j;
  int nb;
  real en,dij,rij,eij,r;
  rvec diff;
  real con = 4*PI*sqrt(PI);

  en = 0;

  for(i=0;i<al->natoms;i++){
    en+=al->rs_Gref[i];
    
    int nn = al->nbonds[i]+al->nb13[i]+al->nb14[i]+al->nnb[i];
    int ndx[nn];
    int count = 0;
    

    for(k=0;k<al->nb14[i];k++){
      j = al->b14[i][k];
      if(i < j) {
        ndx[count] = j;
        count++;
      }
    }
    
    for(k=0;k<al->nnb[i];k++){
      j = al->nb[i][k];

      if(i < j) {
        ndx[count] = j;
        count++;
      }
    }
    real en_per_at = 0;
    
    for(k=0;k<count;k++){
      j = ndx[k];
      rvec_sub(al->x[i],al->x[j],diff);
      dij = norm(diff);
      rij = al->rs_comb[al->rs_type[i]][al->rs_type[j]]; 
      real xi = (dij-al->rs_rad[i])/al->rs_lambda[i];
      real xj = (dij-al->rs_rad[j])/al->rs_lambda[j];
      
      real solv = -2*al->rs_G[i]/(con*al->rs_lambda[i]*rij*rij) * 
        exp(-xi)*al->rs_V[j]
        -2*al->rs_G[j]/(con*al->rs_lambda[j]*rij*rij) * 
        exp(-xj)*al->rs_V[i];
     
      en_per_at+=solv;
      
      en += solv;/* -2*al->rs_G[i]/(con*al->rs_lambda[i]*rij*rij) * exp(-dij*dij)*al->rs_V[j] */
/*           -2*al->rs_G[j]/(con*al->rs_lambda[j]*rij*rij) * exp(-dij*dij)*al->rs_V[i]; */
    }
/*     printf("solv = %g (%s/%s)\n",en_per_at,al->name[i],al->resname[i]);  */
  }
  return en;
}

/*==================================================================*/

real rs_hbond_energy(t_atomlist *al)
{
  /* ROSETTA empirical hydrogen bond energy */

  int i,k,j,l;
  int nb;
  real en,dij,rij,eij,r;
  rvec diff;
  
  real A,B,C,D; /* fitted parameters */
  real lcut = 1.725;
  real ucut = 2.625;
  A = -43.4069;
  B = 288.415;
  C = -631.17;
  D = 454.646;

  /* backbone hbond distance energy */

  en = 0;
  real bben = 0.; 
  for(i=0;i<al->natoms;i++){
    if(al->rs_type[i] == 24){ /* HN */
      for(k=0;k<al->nnb[i];k++){
        j = al->nb[i][k];
        if(al->rs_type[j] == 20){ /* O */
          rvec_sub(al->x[i],al->x[j],diff);
          dij = norm(diff);
          if(dij > lcut && dij < ucut){
/*             printf("type1 = %s type2 = %s\n",al->type[i],al->type[j]); */
            en += A*pow(dij,3) + B*dij*dij + C*dij + D;
            bben += A*pow(dij,3) + B*dij*dij + C*dij + D;
/*             printf("dist = %g energy = %g\n",dij,A*pow(dij,3) + B*dij*dij + C*dij + D); */
          }
        }
      }
    }
  }
/*   printf("bb-bb hb energy = %g\n",bben); */
  /* sc-sc hbond distance energy */
  real scen = 0;
  real scsc_pol[12];
  int pol_size = 12;

  scsc_pol[0]  = 46925.71568;
  scsc_pol[1]  = -161091.865038;
  scsc_pol[2]  = 227015.492901;
  scsc_pol[3]  = -162355.553553;
  scsc_pol[4]  = 53871.6235556;
  scsc_pol[5]  = 527.755117727;
  scsc_pol[6]  = -5788.5283664;
  scsc_pol[7]  = 862.027651106;
  scsc_pol[8]  = 528.676659945;
  scsc_pol[9]  = -239.974865312;
  scsc_pol[10] = 38.637524896;
  scsc_pol[11] = -2.29854866567;
  

  lcut = 1.475;
  ucut = 2.9;


  for(i=0;i<al->natoms;i++){
    if(al->rs_type[i] == 21){ /* polar H */
      for(k=0;k<al->nnb[i];k++){
        j = al->nb[i][k];
        if(al->rs_type[j] == 14 ||
           al->rs_type[j] == 15 ){ /* side-chain oxygen sp2 */
          rvec_sub(al->x[i],al->x[j],diff);
          dij = norm(diff);
          if(dij > lcut && dij < ucut){
            for(l=0;l<pol_size;l++){
              en += scsc_pol[l]*pow(dij,l);
              scen+=scsc_pol[l]*pow(dij,l);
            }
          }
        }
      }
    }
  }
  printf("sc-sc sp2 hb energy = %g\n",scen);


  scen = 0;

  scsc_pol[0]  = -4831.05379014;
  scsc_pol[1]  = 11510.9163482;
  scsc_pol[2]  = -10128.9341084;
  scsc_pol[3]  = 3551.04339365;
  scsc_pol[4]  = -4.83133965718;
  scsc_pol[5]  = -267.356919403;
  scsc_pol[6]  = 39.7418286802;
  scsc_pol[7]  = -14.9890322597;
  scsc_pol[8]  = 14.5382926298;
  scsc_pol[9]  = -4.44513524623;
  scsc_pol[10] = 0.489095106584;
  scsc_pol[11] = -0.0114212858522;

  lcut = 1.625;
  ucut = 2.9;


  for(i=0;i<al->natoms;i++){
    if(al->rs_type[i] == 21){ /* polar H */
      for(k=0;k<al->nnb[i];k++){
        j = al->nb[i][k];
        if(al->rs_type[j] == 13){ /* sp3 oxygen */

          rvec_sub(al->x[i],al->x[j],diff);
          dij = norm(diff);
          if(dij > lcut && dij < ucut){
            for(l=0;l<pol_size;l++){
              en += scsc_pol[l]*pow(dij,l);
              scen+=scsc_pol[l]*pow(dij,l);
            }
          }
        }
      }
    }
  }
/*   printf("sc-sc sp3 hb energy = %g\n",scen); */

  return en;
}

/*****************************************************/
real bb_hbond_energy(t_atomlist *al)
{
/*   FILE *fp = cnclib("bb_dist_energy.dat"); */
/*   t_histogram *hb_dist = read_histogram(fp); */
/*   fp = cnclib("bb_angle_energy.dat"); */
/*   t_histogram *hb_angle = read_histogram(fp); */

  int i,k;
  real tot_energy = 0;
  real energy = 0.;
  real d, a;
  
  /* calculate bb hbond energy */

  for(i=0;i<al->natoms;i++){
    if(strcmp(al->name[i]," H  ") == 0){
      for (k=0;k<al->nnb[i];k++){
        if(strcmp(al->name[al->nb[i][k]]," O  ") == 0){
          d = DIST(al,i,al->nb[i][k]);
          real rij = 1.9;
          real eij = 5.;
          energy = lj_energy_ij(d,rij,eij);
/*           printf("energy = %g\n",energy); */
          
/*           a = RAD2DEG*angle_ij_ik(al->x[i],al->x[al->bonds[i][0]],al->x[al->nb[i][k]]); */
          
/*           if(d < hb_dist->x[0]) { */
/*             energy = hb_dist->y[0]; */
/*           } */
/*           else if(d > hb_dist->x[hb_dist->n-1]){ */
/*             energy = hb_dist->y[hb_dist->n-1]; */
/*           } */
/*           else { */
/*             energy = energy_from_histogram(hb_dist,d); */
/*           } */
/*           tot_energy+=energy; */
/*           if(a < hb_angle->x[0]) { */
/*             energy = hb_angle->y[0]; */
/*           } */
/*           else { */
/*             energy = energy_from_histogram(hb_angle,a); */
/*           } */
          tot_energy+=energy;
          
        }
      }
    }
  }
  return tot_energy;
}

