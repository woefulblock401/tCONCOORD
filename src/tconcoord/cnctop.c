#include <tconcoord.h>
/*===========================================================*/

void read_planar_groups(FILE *fp, t_idxgroups *pl)
{
  char line[STRLEN];
  int count;
  int i;
  count = 0;
  rewind(fp);
  while(get_a_line(fp,line,STRLEN)){
    if(strncmp(line,"PL",2) == 0){
      count++;
      pl = idx_realloc(pl,count);
#ifdef GMX_DOUBLE
      sscanf(line+2,"%d %lf",&pl->natoms[count-1],&pl->val[count-1]);
#else
      sscanf(line+2,"%d %f",&pl->natoms[count-1],&pl->val[count-1]);
#endif
      snew(pl->atoms[count-1],pl->natoms[count-1]);
      for(i=0;i<pl->natoms[count-1];i++){
        get_a_line(fp,line,STRLEN);
        sscanf(line,"%d",&pl->atoms[count-1][i]);
        pl->atoms[count-1][i]-=1;
      }
    }
  }
  pl->n=count;
/*   return pl; */
}
/*===========================================================*/

void read_impropers(FILE *fp, t_idxgroups *imp)
{
  char line[STRLEN];
  int count;
  int i;
  count = 0;
  rewind(fp);
  while(get_a_line(fp,line,STRLEN)){
     if(strncmp(line,"IM",2) == 0){
      count++;
      imp = idx_realloc(imp,count);
#ifdef GMX_DOUBLE
      sscanf(line+2,"%d %lf",&imp->natoms[count-1],&imp->val[count-1]);
#else
      sscanf(line+2,"%d %f",&imp->natoms[count-1],&imp->val[count-1]);
#endif
      snew(imp->atoms[count-1],imp->natoms[count-1]);
      for(i=0;i<imp->natoms[count-1];i++){
        get_a_line(fp,line,STRLEN);
        sscanf(line,"%d",&imp->atoms[count-1][i]);
        imp->atoms[count-1][i]-=1;
      }
    }
  }
  imp->n=count;
}
/*===========================================================*/

void read_cnctop(char *filename, t_atomlist *al,t_idxgroups *acc, t_idxgroups *don,
                 t_idxgroups *phob, t_idxgroups *pl, t_idxgroups *imp)
{
  
  char line[STRLEN];
  int natoms;
  int i,k,l;
  int at;
  char error[STRLEN];
  int nel;
#define MAXPTR 254
  char *ptr[MAXPTR];
  FILE *fp = ffopen(filename,"r");
  acc = idx_realloc(acc,1);
  don = idx_realloc(don,1);
  phob = idx_realloc(phob,1);
  t_vdwcomb *vdwcomb = vdwcomb_init();    

  /* read atoms */
  if(!find_key_word(fp,"[ ATOMS ]")){
    fatal_error("Could not find key word \'ATOMS\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&natoms);
/*       al = atomlist_init(); */
      al = al_realloc(al,natoms);
      al->natoms = natoms;
      for(i=0;i<natoms;i++){
        get_a_line(fp,line,STRLEN);
#ifdef GMX_DOUBLE
        sscanf(line,"%d %s %s %s %d %lf %lf %lf %s %lf %lf %lf %lf %s %d %d %d",
               &al->id[i],al->name[i],al->resname[i],al->chain[i],&al->resid[i],
               &al->x[i][0],&al->x[i][1],&al->x[i][2],al->type[i],&al->vdw[i],
               &al->vdw14[i],&al->m[i],&al->q[i],al->symbol[i],&al->order[i],&al->isposres[i],&al->ptype[i]);
#else
        sscanf(line,"%d %s %s %s %d %f %f %f %s %f %f %f %f %s %d %d %d",
               &al->id[i],al->name[i],al->resname[i],al->chain[i],&al->resid[i],
               &al->x[i][0],&al->x[i][1],&al->x[i][2],al->type[i],&al->vdw[i],
               &al->vdw14[i],&al->m[i],&al->q[i],al->symbol[i],&al->order[i],&al->isposres[i],&al->ptype[i]);
#endif
        extend_name(al->name[i]);
      }
    }
  }
  /* read bonds */
  rewind(fp);

  
  if(!find_key_word(fp,"[ BONDS ]")){
    fatal_error("Could not find key word \'BONDS\' in topology\n");
  }
  while(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
    int at1,at2;
    sscanf(line,"%d %d",&at1,&at2);
    at1-=1;
    at2-=1;
    al->nbonds[at1]+=1;
    al->nbonds[at2]+=1;
    al->bonds[at1][al->nbonds[at1]-1] = at2;
    al->bonds[at2][al->nbonds[at2]-1] = at1;
  }
  rewind(fp);

  if(!find_key_word(fp,"[ B13 ]")){
    fatal_error("Could not find key word \'B13\' in topology\n");
  }
  while(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
    int at1,at2;
    sscanf(line,"%d %d",&at1,&at2);
    at1-=1;
    at2-=1;
    al->nb13[at1]+=1;
    al->nb13[at2]+=1;
    al->b13[at1][al->nb13[at1]-1] = at2;
    al->b13[at2][al->nb13[at2]-1] = at1;
  }
  rewind(fp);
  
  if(!find_key_word(fp,"[ B14 ]")){
    fatal_error("Could not find key word \'B14\' in topology\n");
  }
  while(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
    int at1,at2;
    sscanf(line,"%d %d",&at1,&at2);
    at1-=1;
    at2-=1;
    al->nb14[at1]+=1;
    al->nb14[at2]+=1;
    al->b14[at1][al->nb14[at1]-1] = at2;
    al->b14[at2][al->nb14[at2]-1] = at1;
  }
  rewind(fp);
  
    
  shstr *sdum = NULL;
  real *rdum = NULL;
  if(!find_key_word(fp,"[ VDWTABLE ]")){
    fatal_error("Could not find key word \'VDWTABLE\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&al->nvdw);
      snew(sdum,al->nvdw);
      snew(rdum,al->nvdw);
      snew(al->vdwtab,al->nvdw);
      for(i=0;i<al->nvdw;i++){
        snew(al->vdwtab[i],al->nvdw);
      }
      for(i=0;i<al->nvdw;i++){
        get_a_line(fp,line,STRLEN);
#ifdef GMX_DOUBLE
        sscanf(line,"%s %lf",sdum[i],&rdum[i]);
#else
        sscanf(line,"%s %f",sdum[i],&rdum[i]);
#endif
      }
    }
  }
  /* build vdwtab */
  for(i=0;i<al->nvdw;i++){
    for(k=0;k<al->nvdw;k++){
      al->vdwtab[i][k] = rdum[i]+rdum[k];
    }
  }
      
  if(!find_key_word(fp,"[ VDWCOMB ]")){
    fatal_error("Could not find key word \'VDWCOMB\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&natoms);
      vdwcomb->n = natoms;
      vdwcomb = vdwcomb_realloc(vdwcomb,natoms);
      for(i=0;i<natoms;i++){
        get_a_line(fp,line,STRLEN);
#ifdef GMX_DOUBLE
        sscanf(line,"%s %s %lf",vdwcomb->type1[i],vdwcomb->type2[i],&vdwcomb->vdw[i]);
#else
        sscanf(line,"%s %s %f",vdwcomb->type1[i],vdwcomb->type2[i],&vdwcomb->vdw[i]);
#endif
      }
    }
  }

  /* correct vdwtab */
  for(i=0;i<al->nvdw;i++){
    for(l=i;l<al->nvdw;l++){
      for(k=0;k<vdwcomb->n;k++){
        if((!strcmp(sdum[i],vdwcomb->type1[k]) &&
            !strcmp(sdum[l],vdwcomb->type2[k])) ||
           (!strcmp(sdum[l],vdwcomb->type1[k]) &&
            !strcmp(sdum[i],vdwcomb->type2[k]))){
          al->vdwtab[i][l] = vdwcomb->vdw[k];
          al->vdwtab[l][i] = vdwcomb->vdw[k];
        }
      }
    }
  }
  sfree(rdum);
  sfree(sdum);


  /* read donors */
  if(!find_key_word(fp,"[ HB DONORS ]")){
    fatal_error("Could not find key word \'HB DONORS\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&natoms);
      for(i=0;i<natoms;i++){
        get_a_line(fp,line,STRLEN);
        if(sscanf(line,"%d",&at)==1)
          add_to_group(don,0,at-1);
        else{
          sprintf(error,"while reading %s\n",filename);
          fatal_error(error);
        }
      }
    }
  }
  /* read acceptors */
  if(!find_key_word(fp,"[ HB ACCEPTORS ]")){
    fatal_error("Could not find key word \'HB ACCEPTORS\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&natoms);
      for(i=0;i<natoms;i++){
        get_a_line(fp,line,STRLEN);
        if(sscanf(line,"%d",&at)==1)
          add_to_group(acc,0,at-1);
        else{
          sprintf(error,"while reading %s\n",filename);
          fatal_error(error);
         }
      }
    }
  }

  /* read hydrophobics */
  if(!find_key_word(fp,"[ HYDROPHOBIC ATOMS ]")){
    fatal_error("Could not find key word \'HYDROPHOBIC ATOMS\' in topology\n");
  }
  else{
    if(get_a_line(fp,line,STRLEN) && strchr(line,'[') == 0){
      sscanf(line,"%d",&natoms);
      for(i=0;i<natoms;i++){
        get_a_line(fp,line,STRLEN);
        if(sscanf(line,"%d",&at)==1)
          add_to_group(phob,0,at-1);
        else{
          sprintf(error,"while reading %s\n",filename);
          fatal_error(error);
         }
      }
    }
  }

  /* read planar groups */
  read_planar_groups(fp,pl);

/* read impropers */
  read_impropers(fp,imp);
  free_vdwcomb(vdwcomb);

  
}
/*===========================================================*/

void write_cnctop(char *filename,t_atomlist *al,t_idxgroups *don,
                  t_idxgroups *acc, t_idxgroups *phob, t_idxgroups *pln,
                  t_idxgroups *imp,t_input *inp, t_vdw *vdw, t_vdwcomb *vdwcomb,
                  int psize, atom_id *pos_id)

/* write tCONCOORD topology file */ 


{
  int i,k;
  FILE *topf = ffopen(filename,"w");
  char *host = getenv("HOSTNAME");
  char *user = getenv("USER");
  time_t t = time(0);
  struct tm *tt = localtime(&t);
  char *ltime = asctime(tt); 

  fprintf(topf,";============================================\n");
  
  fprintf(topf,";  >> tCONCOORD TOPOLOGY FILE <<\n;\n");
  fprintf(topf,";   BUILD ON HOST...........: %s\n",host);
  fprintf(topf,";   USER....................: %s\n",user);
  fprintf(topf,";   TIME....................: %s;\n",ltime);
  
  fprintf(topf,";============================================\n");
  
  fprintf(topf,"\n[ ATOMS ]\n");
  fprintf(topf,"%10d ; number of atoms\n",al->natoms);
  fprintf(topf,";      id  name resname  chain  resid     x       y       z      type     vdw     vdw14      m       q   symbol  order  posres  ptype hyb\n");
  for(i=0;i<al->natoms;i++){
    char chain[STRLEN];
    if(strcmp(al->chain[i]," ") == 0)  strcpy(chain,"0");
    else strcpy(chain,al->chain[i]);
    fprintf(topf,"%8d %6s %6s   %3s %6d %8.3f %8.3f %8.3f %6s %8.3f %8.3f %8.3f %8.3f %4s %6d %6d %6d %5s\n",
            al->id[i],al->name[i],al->resname[i],chain,al->resid[i],al->x[i][0],
            al->x[i][1],al->x[i][2],al->type[i],al->vdw[i],al->vdw14[i],al->m[i],al->q[i],al->symbol[i],al->order[i],al->isposres[i],al->ptype[i],al->hyb[i]);
  }
  fprintf(topf,"\n[ BONDS ]\n");
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      if(i < al->bonds[i][k]){
        fprintf(topf,"%8d %8d\n",al->id[i],al->id[al->bonds[i][k]]);
      }
    }
  }
  fprintf(topf,"\n[ B13 ]\n");
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nb13[i];k++){
      if(i < al->b13[i][k]){
        fprintf(topf,"%8d %8d\n",al->id[i],al->id[al->b13[i][k]]);
      }
    }
  }
  fprintf(topf,"\n[ B14 ]\n");
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nb14[i];k++){
      if(i < al->b14[i][k]){
        fprintf(topf,"%8d %8d\n",al->id[i],al->id[al->b14[i][k]]);
      }
    }
  }
  
  fprintf(topf,"\n[ VDWTABLE ]\n");
  fprintf(topf,"%10d ; size of vdw table\n",vdw->n);
  for(i=0;i<vdw->n;i++){
    fprintf(topf,"%6s %8.3f\n",vdw->type[i],vdw->vdw[i]);
  }
  fprintf(topf,"\n[ VDWCOMB ]\n");
  fprintf(topf,"%10d ; number of vdw combinations\n",vdwcomb->n);
  fprintf(topf,";      type1     type2     vdwcomb\n");
  for(i=0;i<vdwcomb->n;i++){
    fprintf(topf,"%8s %8s %8.3f\n",vdwcomb->type1[i],vdwcomb->type2[i],vdwcomb->vdw[i]);
  }

/*   fprintf(topf,"\n[ GROUPS ]\n"); */
/*   fprintf(topf,"%10d ; number of atoms\n",al->natoms); */
/*   for(i=0;i<al->natoms;i++){ */
/*     fprintf(topf,"%8d %5d ",al->id[i],al->ngrps[i]); */
/*     for(k=0;k<al->ngrps[i];k++){ */
/*       fprintf(topf," %6d ", al->grpnr[i][k]); */
/*     } */
/*     fprintf(topf,"\n"); */
/*   } */

  int ndon = 0;
  for(i=0;i<al->natoms;i++){
    if(al->isdon[i]) ndon++;
  }
  
  fprintf(topf,"\n[ HB DONORS ]\n");
  fprintf(topf,"%10d ; number of donors\n",ndon);
  for(i=0;i<al->natoms;i++){
    if(al->isdon[i])
      fprintf(topf,"%10d ;%s(%d%s)\n",al->id[i],al->name[i],
              al->resid[i],al->resname[i]);
  }

  int nacc = 0;
  for(i=0;i<al->natoms;i++){
    if(al->isacc[i]) nacc++;
  }
  fprintf(topf,"\n[ HB ACCEPTORS ]\n");
  fprintf(topf,"%10d ; number of acceptors\n",nacc);

  for(i=0;i<al->natoms;i++){
    if(al->isacc[i])
      fprintf(topf,"%10d ;%s(%d%s)\n",al->id[i],al->name[i],
              al->resid[i],al->resname[i]);
  }

/*   fprintf(topf,"\n[ HB DONORS ]\n"); */
/*   fprintf(topf,"%10d ; number of donors\n",don->natoms[0]); */
/*   for(i=0;i<don->natoms[0];i++){ */
/*     fprintf(topf,"%10d ;%s(%d%s)\n",al->id[don->atoms[0][i]],al->name[don->atoms[0][i]], */
/*             al->resid[don->atoms[0][i]],al->resname[don->atoms[0][i]]); */
/*   } */

/*   fprintf(topf,"\n[ HB ACCEPTORS ]\n"); */
/*   fprintf(topf,"%10d ; number of accpetors\n",acc->natoms[0]); */
/*   for(i=0;i<acc->natoms[0];i++){ */
/*     fprintf(topf,"%10d ;%s(%d%s)\n",al->id[acc->atoms[0][i]],al->name[acc->atoms[0][i]], */
/*             al->resid[acc->atoms[0][i]],al->resname[acc->atoms[0][i]]); */
/*   } */

  int nphob = 0;
  for(i=0;i<al->natoms;i++){
    if(al->ishphob[i]) nphob++;
  }
  fprintf(topf,"\n[ HYDROPHOBIC ATOMS ]\n");
  fprintf(topf,"%10d ; number of hydrophobic atoms\n",nphob);

  for(i=0;i<al->natoms;i++){
    if(al->ishphob[i])
      fprintf(topf,"%10d ;%s(%d%s)\n",al->id[i],al->name[i],
              al->resid[i],al->resname[i]);
  }
  
  
/*   fprintf(topf,"\n[ HYDROPHOBIC ATOMS ]\n"); */
/*   fprintf(topf,"%10d ; number of hydrophobic atoms\n",phob->natoms[0]); */
/*   for(i=0;i<phob->natoms[0];i++){ */
/*     fprintf(topf,"%10d ;%s(%d%s)\n",al->id[phob->atoms[0][i]],al->name[phob->atoms[0][i]], */
/*             al->resid[phob->atoms[0][i]],al->resname[phob->atoms[0][i]]); */
/*   } */

  fprintf(topf,"\n[ PLANAR GROUPS ]\n");
  for(i=0;i<pln->n;i++){
     if(pln->flag[i]){  
      fprintf(topf,"PL %-8d %g  ;%d(%s)\n",pln->natoms[i],pln->val[i],al->resid[pln->atoms[i][0]],al->resname[pln->atoms[i][0]]);
      for(k=0;k<pln->natoms[i];k++){
        fprintf(topf,"%8d  ; %s(%d%s)\n",al->id[pln->atoms[i][k]],
                al->name[pln->atoms[i][k]],al->resid[pln->atoms[i][k]],al->resname[pln->atoms[i][k]]);
       }  
    }
  }
  fprintf(topf,"\n[ IMPROPERS ]\n");
  for(i=0;i<imp->n;i++){
    if(imp->flag[i]){
      fprintf(topf,"IM %-8d %g  ;%d(%s)\n",imp->natoms[i],imp->val[i],al->resid[imp->atoms[i][0]],al->resname[imp->atoms[i][0]]);
      for(k=0;k<imp->natoms[i];k++){
        fprintf(topf,"%8d  ; %s(%d%s) %d\n",al->id[imp->atoms[i][k]],
                al->name[imp->atoms[i][k]],al->resid[imp->atoms[i][k]],
                al->resname[imp->atoms[i][k]],al->isposres[imp->atoms[i][k]]);
      }
    }
  }

}

