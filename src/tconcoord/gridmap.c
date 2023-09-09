#include <tconcoord.h>

/*============================================================*/

t_gridmap *gridmap_init(void)
{
  int i;
  t_gridmap *map = NULL;
  snew(map,1);
  strcpy(map->name,"map");
  for(i=0;i<DIM;i++){
    map->npts[i] = 0;
    map->n[i] = 0;
    map->origin[i] = 0.;
    map->center[i] = 0.;
  }
  map->nelem = 0;
  map->spacing = 0.;
  map->inv_spacing = 0.;
  map->precision = 0.0001;
  strcpy(map->datafile,"");
  strcpy(map->molecule,"");
  strcpy(map->paramfile,"");
  snew(map->cell,1);
  snew(map->natom,1);
  snew(map->values,1);

  return map;
}

t_gridmap *reset_grid(t_gridmap *map)
{
  int i;
  if(map!=NULL) {
    
    for(i=0;i<map->nelem;i++){
      sfree(map->cell[i]);
    }
    sfree(map->cell);
    sfree(map->values);
    sfree(map->natom);
    
    for(i=0;i<DIM;i++){
      map->npts[i] = 0;
      map->n[i] = 0;
      map->origin[i] = 0.;
      map->center[i] = 0.;
    }
    map->nelem = 0;
    map->spacing = 0.;
    map->inv_spacing = 0.;
    map->precision = 0.0001;
    strcpy(map->datafile,"");
    strcpy(map->molecule,"");
    strcpy(map->paramfile,"");
    snew(map->cell,1);
    snew(map->natom,1);
    snew(map->values,1);
  }
  else {
    map = gridmap_init();
  }
  
  return map;
}


t_gridmap *read_map(char *filename, FILE *log,bool bVerbose)
{
  char line[STRLEN];
  int count;
  int i;
  t_gridmap *map = gridmap_init();
  char dum[STRLEN];
  FILE *fp = ffopen(filename,"r");
  
  
  i = 0;count = 0;
  
  if(log==NULL) log = stderr;  
  while(get_a_line(fp,line,STRLEN)){
    if (i==0) {
      sscanf(line,"%s %s",dum,map->paramfile);
    }
    else if (i==1) {
      sscanf(line,"%s %s",dum,map->datafile);
    }
    else if (i==2) {
      sscanf(line,"%s %s",dum,map->molecule);
    }
    else if(i==3) {
#ifdef GMX_DOUBLE
      sscanf(line,"%s %lf",dum,&map->spacing);
#else
      sscanf(line,"%s %f",dum,&map->spacing);
#endif
    }
    else if(i==4){
      sscanf(line,"%s %d %d %d",dum,&map->npts[XX],&map->npts[YY],&map->npts[ZZ]);
    }
    else if(i==5){
#ifdef GMX_DOUBLE
      sscanf(line,"%s %lf %lf %lf",dum,&map->center[XX],&map->center[YY],&map->center[ZZ]);
#else         
      sscanf(line,"%s %f %f %f",dum,&map->center[XX],&map->center[YY],&map->center[ZZ]);
#endif
    }
    else {
      count++;
    }
    if(bVerbose && i<6) {
      fprintf(stderr,"Reading %s\n",dum);
    }
    i++;
  }
  if(bVerbose){
    fprintf(stderr,"Will read %d data points\n",count);
  }
  
  snew(map->values,count);
  rewind(fp);
  i = 0;count = 0;
  while(get_a_line(fp,line,STRLEN)){
    if(i>5){
#ifdef GMX_DOUBLE
      sscanf(line,"%lf",&map->values[count]);
#else
      sscanf(line,"%f",&map->values[count]);
#endif
      count++;
    }
    i++;
  }
  
  for(i=0;i<DIM;i++){
    map->n[i] = map->npts[i]+1;
  }
  map->nelem = map->n[XX]*map->n[YY]*map->n[ZZ];
  map->inv_spacing = 1./map->spacing;
  for(i=0;i<DIM;i++){
    map->origin[i] = map->center[i]-map->npts[i]/2*map->spacing;
  }
  return map;
}


void write_dx_map(char *filename, t_gridmap *map)
{
  int i,j,k;
  int nx,ny,nz;
  int col = 0;
  FILE *fp = ffopen(filename,"w");
  
  nx = map->n[XX];
  ny = map->n[YY];
  nz = map->n[ZZ];
  
  fprintf(fp,"#==================================\n");
  fprintf(fp,"# AutoGrid Map File: %s\n",map->name);
  fprintf(fp,"# Receptor File Name: %s\n",map->molecule);
  fprintf(fp,"#==================================\n");
  fprintf(fp,"object 1 class gridpositions counts %d %d %d\n",nx,ny,nz);
  fprintf(fp,"origin %12.5E %12.5E %12.5E\n",map->origin[XX],map->origin[YY],map->origin[ZZ]);
  fprintf(fp,"delta %12.5E %12.5E %12.5E\n",map->spacing,0.,0.);
  fprintf(fp,"delta %12.5E %12.5E %12.5E\n",0.,map->spacing,0.);
  fprintf(fp,"delta %12.5E %12.5E %12.5E\n",0.,0.,map->spacing);
  fprintf(fp,"object 2 class gridconnections counts %d %d %d\n",nx,ny,nz);
  fprintf(fp,"object 3 class array type double rank 0 items %d data follows\n",map->nelem);
  
  for(k=0;k<nz;k++){
    col = 0;
    for(j=0;j<ny;j++){
      for(i=0;i<nx;i++){
        fprintf(fp," %12.5E",map->values[i*ny*nz + j*nz + k]);
        col++;
        if(col==3) 
        {
          fprintf(fp,"\n");
          col = 0;
        }
      }
    }
  }
  if(col) fprintf(fp,"\n");
  
  fprintf(fp,"attribute \"dep\" string \"positions\"\n");
  fprintf(fp,"object \"regular positions regular connections\" class field\n");
  fprintf(fp,"component \"positions\" value 1\n");
  fprintf(fp,"component \"connections\" value 2\n");
  fprintf(fp,"component \"data\" value 3\n");
  fclose(fp);
}



       
         
