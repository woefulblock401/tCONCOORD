#include <tconcoord.h>
  
/*================================================================*/
FILE *cnclib(char *filename)
{
  FILE *fp = NULL;
  fp=fopen(filename,"r");
  if(fp==NULL){
    longstr dummy;
    char *path = getenv("CNCLIB");
   if(path==NULL){
      fatal_error("environmental variable \'CNCLIB\' not declared");
    }
 
    strcpy(dummy,path);
    strcat(dummy,"/");
    strcat(dummy,filename);
/*     printf("path = %s\n",dummy); */
    fp = fopen(dummy,"r");
    if(fp==NULL){
      char error[STRLEN];
      sprintf(error,"File %s not found\n",dummy);
      fatal_error(error);
    }
    else
    {
      fprintf(stderr,"tCNC__lib_> Reading tCONCOORD lib \'%s\'\n",dummy);
    }
    
  }
  else
  {
    fprintf(stderr,"tCNC__lib_> Reading tCONCOORD lib \'%s\'\n",filename);
  }
  fflush(stderr);
  return fp;
}
/*================================================================*/
  

void get_cnc_solv(t_atomlist *al)
{
  FILE *fp = cnclib("Solvation.dat");
  char line[STRLEN];
  shstr resn[200];
  shstr name[200];
  real solv[200];
  int n = 0;
  while(get_a_line(fp,line,STRLEN)){
#ifdef GMX_DOUBLE
    sscanf(line,"%s %s %lf",resn[n],name[n],&solv[n]);
#else
    sscanf(line,"%s %s %f",resn[n],name[n],&solv[n]);
#endif
    extend_name(name[n]);
/*     printf("reading '%s' '%s' %g\n",resn[n],name[n],solv[n]); */
    n++;
  }
  int i,k;
  for(i=0;i<al->natoms;i++){
    for(k=0;k<n;k++){
      if(strcmp(al->resname[i],resn[k]) == 0 &&
         strcmp(al->name[i],name[k]) == 0) al->cnc_solv[i] = solv[k];
    }
  }
/*   for(i=0;i<al->natoms;i++){ */
/*     printf("name = %s / %s  solv = %g\n",al->name[i],al->resname[i],al->cnc_solv[i]); */
/*   } */
}

/*================================================================*/

void get_types(FILE *log,t_atomlist *al, t_types *tp)
{
  int i,k;
  char warn[STRLEN];
  
  /* first do generic types */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(tp->resname[k],"GENERIC")==0){
        if(strcmp(al->symbol[i],tp->name[k])==0){
          strcpy(al->type[i],tp->type[k]);
        }
      }
    }
  }
  /* now do the defaults */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(tp->resname[k],"DEFAULT")==0){
        if(IsProtein(al,i) && strcmp(al->name[i],tp->name[k])==0){
          strcpy(al->type[i],tp->type[k]);
          strcpy(al->hyb[i],tp->hyb[k]);
        }
      }
    }
  }

  /* overwrite with specific entries */

  for(i=0;i<al->natoms;i++){
    for(k=0;k<tp->n;k++){
      if(strcmp(tp->resname[k],al->resname[i])==0 &&
         strcmp(tp->name[k],al->name[i])==0){
        strcpy(al->type[i],tp->type[k]);
        strcpy(al->hyb[i],tp->hyb[k]);

       }
    }
  }

  /* special cases */
  for(i=0;i<al->natoms;i++){
    if(strcmp(al->resname[i],"HIS")==0 &&
       strcmp(al->name[i]," ND1") ==0){
      if(al->nbonds[i]==3) strcpy(al->type[i],"NH1");
      else strcpy(al->type[i],"NHS");
    }
    if(strcmp(al->resname[i],"HIS")==0 &&
       strcmp(al->name[i]," NE2") ==0){
      if(al->nbonds[i]==3) strcpy(al->type[i],"NH1");
      else strcpy(al->type[i],"NHS");
    }
    if(strcmp(al->name[i]," OXT")==0){
      strcpy(al->type[i],"OC");
      strcpy(al->hyb[i],"sp2");
    }
    if(strcmp(al->name[i],"1H  ")==0 ||
       strcmp(al->name[i],"2H  ")==0 ||
       strcmp(al->name[i],"3H  ")==0)
      strcpy(al->type[i],"HC");
    if(strcmp(al->name[i]," N  ") == 0 &&
       (al->resid[i] == 1 || strcmp(al->chain[i-1],al->chain[i])!=0)){
      strcpy(al->hyb[i],"sp3");
/*       strcpy(al->type[i],"NH3"); */
    }
    
  }
  
  for(i=0;i<al->natoms;i++){
    if(strstr(al->type[i],"GEN")!=0){
      sprintf(warn,"Using generic atom type %s for atom %s(%d)%s\n",
              al->type[i],al->name[i],al->resid[i],al->resname[i]);
      CNCwarn(log,warn);
      
    }
  }
}

/*============================================================*/

void get_vdw_radii(FILE *log,t_atomlist *al, t_vdw *vdw, 
                   t_vdwcomb *vdwcomb, t_vdwcomb *vdw14comb)
{
  int i,k,l;
  bool check;
  char warn[STRLEN];
  
/*   fprintf(stderr,"Dimensions of vdw table %dx%d\n",vdw->n,vdw->n); */

  for(i=0;i<al->natoms;i++){
    al->ptype[i] = -1;
    check = FALSE;
    for(k=0;k<vdw->n;k++){
      if(strcmp(al->type[i],vdw->type[k])==0){
        al->vdw[i] = vdw->vdw[k];
        al->ptype[i] = k;
        al->vdw14[i] = vdw->vdw14[k];
        check = TRUE;
      }
    }
    if(!check){
      sprintf(warn,"Could not assign vdw radius to atom %s %s (%s%d%s)\n",
              al->type[i],al->symbol[i],al->name[i],al->resid[i],al->resname[i]);
      CNCwarn(log,warn);
    }
  }

  /* build lookup table */

  al->nvdw = vdw->n;
  snew(al->vdwtab,vdw->n);
  snew(al->vdw14tab,vdw->n);
  for(i=0;i<vdw->n;i++){
    snew(al->vdwtab[i],vdw->n);
    snew(al->vdw14tab[i],vdw->n);
  }
  for(i=0;i<vdw->n;i++){
    for(k=0;k<vdw->n;k++){
      al->vdwtab[i][k] = vdw->vdw[i]+vdw->vdw[k];
      al->vdw14tab[i][k] = vdw->vdw14[i]+vdw->vdw14[k];
    }
  }

  for(i=0;i<vdw->n-1;i++){
    for(l=i+1;l<vdw->n;l++){
      for(k=0;k<vdwcomb->n;k++){
        if((!strcmp(vdw->type[i],vdwcomb->type1[k]) &&
            !strcmp(vdw->type[l],vdwcomb->type2[k])) ||
           (!strcmp(vdw->type[l],vdwcomb->type1[k]) &&
            !strcmp(vdw->type[i],vdwcomb->type2[k]))){
          al->vdwtab[i][l] = vdwcomb->vdw[k];
          al->vdwtab[l][i] = vdwcomb->vdw[k];
        }
      }
    }
  }

  for(i=0;i<vdw->n-1;i++){
    for(l=i+1;l<vdw->n;l++){
      for(k=0;k<vdw14comb->n;k++){
        if((!strcmp(vdw->type[i],vdw14comb->type1[k]) &&
            !strcmp(vdw->type[l],vdw14comb->type2[k])) ||
           (!strcmp(vdw->type[l],vdw14comb->type1[k]) &&
            !strcmp(vdw->type[i],vdw14comb->type2[k]))){
          al->vdw14tab[i][l] = vdw14comb->vdw[k];
          al->vdw14tab[l][i] = vdw14comb->vdw[k];
          check = TRUE;
        }
      }
    }
  }
}


/*============================================================*/

t_area *read_area(char *filename)
{
  int ninp;
  t_area *ta = NULL;
  snew(ta,1);
  int i;
  t_inpfile *inf = read_inpfile(filename,&ninp);
  for(i=0;i<ninp;i++){
    if(strcmp(inf[i].name,"xmax") == 0)
      ta->xmax = atof(inf[i].value);
    else if(strcmp(inf[i].name,"xmin") == 0)
      ta->xmin = atof(inf[i].value);
    else if(strcmp(inf[i].name,"ymax") == 0)
      ta->ymax = atof(inf[i].value);
    else if(strcmp(inf[i].name,"ymin") == 0)
      ta->ymin = atof(inf[i].value);
    else if(strcmp(inf[i].name,"zmax") == 0)
      ta->zmax = atof(inf[i].value);
    else if(strcmp(inf[i].name,"zmin") == 0)
      ta->zmin = atof(inf[i].value);
  }
  return ta;
}

/*============================================================*/

