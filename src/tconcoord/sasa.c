#include <tconcoord.h>

void get_sasa_types(t_atomlist *al)
{
  FILE *fp = cnclib("sasa_types.dat");
  
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

  /* assign sasa_types to atomlist */

  /* now do the defaults */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(tp->resname[k],"DEFAULT")==0){
        if(strcmp(al->name[i],tp->name[k])==0){
          al->sasa_type[i] = tp->rs_type[k];
        }
      }
    }
  }

  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(al->name[i],tp->name[k]) == 0 &&
         strcmp(al->resname[i],tp->resname[k]) == 0){
        al->sasa_type[i] = tp->rs_type[k];
      }
    }
  }

  free_types(tp);
  
  fp = cnclib("sasa_params.dat");
  
  t_vdw *vdw = vdw_init();
  
  i = 0;
  if(!find_key_word(fp,"[ PARAMETERS ]")){
    fatal_error("Shitty file!\n");
  }
  
  while(get_a_line(fp,line,STRLEN) &&
        strchr(line,'[')==NULL)
  {

    /* we read that in the vdw types
       rs_rad -> charmm_rad
       rs_eps -> R_i
       rs_G -> p_i
       rs_V -> sig_i
    */
#ifdef GMX_DOUBLE
    sscanf(line,"%d %lf %lf %lf %lf",&vdw->rs_type[i],
           &vdw->rs_rad[i],&vdw->rs_eps[i],&vdw->rs_G[i],
           &vdw->rs_V[i]); 
#else
    sscanf(line,"%d %f %f %f %f",&vdw->rs_type[i],
           &vdw->rs_rad[i],&vdw->rs_eps[i],&vdw->rs_G[i],
           &vdw->rs_V[i]);
#endif
    i++;
  }
  vdw->n = i;

/* assign values to atomlist */

  for(i=0;i<al->natoms;i++){
    for(k=0;k<vdw->n;k++){
      if(al->sasa_type[i] == vdw->rs_type[k]){
        al->sasa_Ri[i] = vdw->rs_eps[k];
        al->sasa_pi[i] = vdw->rs_G[k];
        al->sasa_sigi[i] = vdw->rs_V[k];
      }
    }
  }

  free_vdw(vdw);  
  fclose(fp);  
}


  
void sasa_per_atom(t_atomlist *al, real *sasa)
{
  /* S_i -> SASA of isolated atom i */

  int i, k, j;
  real Rprobe = 1.4;
  real S_i, R_i, R_j, A_i, rij,pi;
  rvec diff;
  
  for(i=0;i<al->natoms;i++)
  {
    if(al->sasa_type[i] != 22 && 
       al->sasa_type[i] != 23)
    {
      
      R_i = al->sasa_Ri[i];
      S_i = 4*PI*sqr(R_i+Rprobe);
      
      
/* get all bonded atoms plus neighbors */
      
      int nn = al->nbonds[i]+al->nb13[i]+al->nb14[i]+al->nnb[i];
      int ndx[nn];
      int count = 0;
      real pij[nn];
      real bij[nn];
      
      for(k=0;k<al->nbonds[i];k++){
        j = al->bonds[i][k];
        if(al->sasa_type[j] != 22 && 
           al->sasa_type[j] != 23)
        {
          ndx[count] = j;
          pij[count] = 0.8875;
          count++;
        }
        
      }
      for(k=0;k<al->nb13[i];k++){
        j = al->b13[i][k];
        if(al->sasa_type[j] != 22 && 
           al->sasa_type[j] != 23)
        {
          ndx[count] = j;
          pij[count] = 0.3516;
          count++;
        }
        
      }
      for(k=0;k<al->nb14[i];k++){
        j = al->b14[i][k];
        if(al->sasa_type[j] != 22 && 
           al->sasa_type[j] != 23)
        {
          ndx[count] = j;
          pij[count] = 0.3516;
          count++;
        }
        
      }
      for(k=0;k<al->nnb[i];k++){
        j = al->nb[i][k];
        if(al->sasa_type[j] != 22 && 
           al->sasa_type[j] != 23)
        {
          ndx[count] = j;
          pij[count] = 0.3516;
          count++;
        }
        
      }
      
      for(k=0;k<count;k++){
        j = ndx[k];
        R_j = al->sasa_Ri[j];
        rvec_sub(al->x[i],al->x[j],diff);
        rij = norm(diff);
        pi = al->sasa_pi[i];
        
        if(rij > R_i + R_j + 2*Rprobe) 
        {
          bij[k] = 0;
        }
        else 
        {
          bij[k] = PI*(R_i + Rprobe)*(R_i + R_j +2*Rprobe-rij)*(1+(R_j-R_i)/rij);
        }
      }
      A_i = 0;
      real tmp = 1;
      
      for(k=0;k<count;k++) 
      {
        tmp=tmp*(1-pi*pij[k]*bij[k]/S_i);
      }
      A_i = S_i*tmp;
      sasa[i] = A_i;
    }
    else 
    {
      sasa[i] = 0;
    }
  }
}

real sasa(t_atomlist *al)
{
  int i;
  real sasum = 0;
  real sasa_per_at[al->natoms];
  sasa_per_atom(al,sasa_per_at);

  for(i=0;i<al->natoms;i++){
    sasum+=sasa_per_at[i];
  }
  return sasum;
}

real sasa_energy(t_atomlist *al)
{
  int i;
  real en = 0;
  real sasa_per_at[al->natoms];
  sasa_per_atom(al,sasa_per_at);

  for(i=0;i<al->natoms;i++){
    en+=al->sasa_sigi[i]*sasa_per_at[i];
  }
  return en;
}

