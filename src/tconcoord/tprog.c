#include <tconcoord.h>

#define MAXNB 500


typedef struct t_maps
{
  int n;
  char **types;
  t_gridmap **map;
}t_maps;


void get_com(t_atomlist *al, rvec co)
{
  int i;
  clear_rvec(co);
  real mtot = 0;
  for(i=0;i<al->natoms;i++) {
    co[XX]+= al->x[i][XX]*al->m[i];
    co[YY]+= al->x[i][YY]*al->m[i];
    co[ZZ]+= al->x[i][ZZ]*al->m[i];
    mtot += al->m[i];
  }
  co[XX] /=mtot;
  co[YY] /=mtot;
  co[ZZ] /=mtot;
}

void shift_coords(rvec *x, int n,rvec co)
{
  int i,k;
  for(i=0;i<n;i++){
    for(k=0;k<DIM;k++){
      x[i][k]-=co[k];
    }
  }
}

void shift_coord_set(rvec **coord_set, int n, int natoms,rvec co)
{
  int i;
  for(i=0;i<n;i++){
    shift_coords(coord_set[i],natoms,co);
  }
}


void get_time(real sec,char *timestr)
{
   
  int d,h,m,s;
  real rest = 1;

  d = (int) sec / 86400.;
  
  rest = (int) sec % 86400;
  
  h = (int) rest / 3600.;
  
  rest = (int)rest % 3600;
  
  m = (int) rest / 60.;
  
  s = (int) rest % 60;
  
  sprintf(timestr,"Total computation time: %d days %d hours %d minutes %d seconds\n",d,h,m,s);
}

  
void perturb_atomlist(t_atomlist *al,gmx_rng_t rng, real scale)
{
  int i,k;
  for(i=0;i<al->natoms;i++){
    if(!al->isposres[i]){
      for(k=0;k<DIM;k++){
        al->x[i][k] += (gmx_rng_uniform_real(rng)-0.5)*scale;
      }
    }
  }
}

void random_rot(rvec *x, int n, gmx_rng_t rng)
{
  /* performs a random rotation */
  real phix,phiy,phiz;
  rvec ex,ey,ez,ori;
  int i;
  clear_rvec(ex);
  clear_rvec(ey);
  clear_rvec(ez);
  clear_rvec(ori);
  ex[XX] = 1.;
  ey[YY] = 1.;
  ez[ZZ] = 1.;

  for(i=0;i<DIM;i++){
    phix = gmx_rng_uniform_real(rng)*2*PI;
    phiy = gmx_rng_uniform_real(rng)*2*PI;
    phiz = gmx_rng_uniform_real(rng)*2*PI;
  }
  for(i=0;i<n;i++){
    rotate_atom(ex,ori,x[i],phix);
  }
  for(i=0;i<n;i++){
    rotate_atom(ey,ori,x[i],phiy);
  }
  for(i=0;i<n;i++){
    rotate_atom(ez,ori,x[i],phiz);
  }
  
}



int iround(real x)
{

  if ( x-(int) x >=.5 ) return (int) (x+1);
  else return (int) x;
}

int gridpoint(real x, real origin, real inv_spacing,int max)
{
  real n;
  int point;
  
  n = (x-origin)*inv_spacing;
  if(n - (int) n > .5) point = (int)n +1;
  else point = (int) n;
  if(point >= max) return point-1;
  else return point;
}




void write_conformers(char *filename,t_atomlist *al, rvec **x, bool *write)
{
  int i,k;
  FILE *fp = ffopen(filename,"w");
  for(i=0;i<1000;i++){
    if(write[i]) {
      for(k=0;k<al->natoms;k++){
        copy_rvec(x[i][k],al->x[k]);
      }
      write_pdb_frame(fp,al,i);
    }
  }
}


void max_coord(rvec *x, int n, matrix max)
{
  /* matrix X stores minimum
     matrix Z stores maximu
     this max[YY][XX] is minimum y value
  */

  int i,k;
  for(i=0;i<DIM;i++){
    max[i][XX] = 1000;
    max[i][YY] = 0.;
    max[i][ZZ] = -1000.;
  }
  for(i=0;i<n;i++){
    if(x[i][XX] < max[XX][XX]) max[XX][XX] = x[i][XX];
    if(x[i][XX] > max[XX][ZZ]) max[XX][ZZ] = x[i][XX];

    if(x[i][YY] < max[YY][XX]) max[YY][XX] = x[i][YY];
    if(x[i][YY] > max[YY][ZZ]) max[YY][ZZ] = x[i][YY];

    if(x[i][ZZ] < max[ZZ][XX]) max[ZZ][XX] = x[i][ZZ];
    if(x[i][ZZ] > max[ZZ][ZZ]) max[ZZ][ZZ] = x[i][ZZ];
  }
}

void coord_set2al(rvec *x, t_atomlist *al)
{
  int i;
  for(i=0;i<al->natoms;i++){
    al->x[i][XX] = x[i][XX];
    al->x[i][YY] = x[i][YY];
    al->x[i][ZZ] = x[i][ZZ];
  }
}
void put_molecule_in_box(t_atomlist *al, rvec *x, 
                         real xmin, real xmax,
                         real ymin, real ymax,
                         real zmin, real zmax)
{
  
  matrix max;
  rvec sh;
  int k;
  
  clear_rvec(sh);
  clear_mat(max);
  
  max_coord(x,al->natoms,max);
  if(max[XX][XX] < xmin) sh[XX] = xmin-max[XX][XX];
  else if(max[XX][ZZ] > xmax) sh[XX] = xmax-(max[XX][ZZ]);
  
  if(max[YY][XX] < ymin) sh[YY] = ymin-max[YY][XX];
  else if(max[YY][ZZ] > ymax) sh[YY] = ymax-(max[YY][ZZ]);
  
  if(max[ZZ][XX] < zmin) sh[ZZ] = zmin-max[ZZ][XX];
  else if(max[ZZ][ZZ] > zmax) sh[ZZ] = zmax-(max[ZZ][ZZ]);
  
  
  for(k=0;k<al->natoms;k++){
    x[k][XX]+=sh[XX];
    x[k][YY]+=sh[YY];
    x[k][ZZ]+=sh[ZZ];
  }
}

void flood_binding_site(t_atomlist *al, rvec **coord_set,
                        int nset, gmx_rng_t rng,
                        real xmin,real xmax,
                        real ymin, real ymax,
                        real zmin, real zmax)
{

  /* the molecules must be centered at com */

  int i,k;
  real rx,ry,rz;
  real xdim,ydim,zdim;
  
  xdim = xmax-xmin;
  ydim = ymax-ymin;
  zdim = zmax-zmin;

  /* put a random sample of
     conformations into box
  */

  for(i=0;i<nset;i++){

    /* put molecules to box min coords */

    for(k=0;k<al->natoms;k++){
      coord_set[i][k][XX]+=xmin;
      coord_set[i][k][YY]+=ymin;
      coord_set[i][k][ZZ]+=zmin;
    }
    /* apply random shift into the box */
    rx = gmx_rng_uniform_real(rng)*xdim;
    ry = gmx_rng_uniform_real(rng)*ydim;
    rz = gmx_rng_uniform_real(rng)*zdim;
    for(k=0;k<al->natoms;k++){
      coord_set[i][k][XX]+=rx;
      coord_set[i][k][YY]+=ry;
      coord_set[i][k][ZZ]+=rz;
    }

    /* check whether all atoms are in box */

    put_molecule_in_box(al,coord_set[i],
                        xmin,xmax,
                        ymin,ymax,
                        zmin,zmax);
    
  }
}

void copy_coord_set(rvec **source, rvec **target, int nset, int natoms)
{
  int i,k,j;
/*   snew(target,nset); */
  for(i=0;i<nset;i++){
/*     snew(target[i],natoms); */
    for(k=0;k<natoms;k++){
      copy_rvec(source[i][k],target[i][k]);
    }
  }
}

void get_lowest_energies(real *energies, int n,
                         real *ener, int nb, int *index)
{
  int i,k,idx;
  real dum;
  
  for(i=0;i<nb;i++) ener[i] = 999999.;

  for(i=0;i<n;i++){
    idx = -1;
    if(energies[i] < ener[0]) idx = 0;
    if(idx==-1){
      for(k=0;k<nb-1;k++){
        if(energies[i] > ener[k] &&
           energies[i] <= ener[k+1]) idx = k+1;
      }
    }
    if(idx!=-1) {
      for(k=nb-1;k>idx;k--){
        ener[k] = ener[k-1];
        index[k] = index[k-1];
      }
      ener[idx] = energies[i];
      index[idx] = i;
    }
  }
}

real calc_energy(t_atomlist *al, t_maps *maps)
{
  int i,j,k;
  int mapid;
  real en;
  int xpts, ypts, zpts, idx;
  int nx,ny,nz;
  real inv_spacing;
  int ELEC = -1;
  int DESOLV = -1;
  inv_spacing = maps->map[0]->inv_spacing;
  nx = maps->map[0]->n[XX]; 
  ny = maps->map[0]->n[YY]; 
  nz = maps->map[0]->n[ZZ]; 
  en = 0;
  for(i=0;i<maps->n;i++){
    if(strcmp(maps->types[i],"e") == 0) ELEC = i;
    else if(strcmp(maps->types[i],"d") == 0) DESOLV = i;
  }
  
  for(k=0;k<al->natoms;k++){  
    /* energy per atom */
    mapid = -1;
    for(j=0;j<maps->n;j++){
      if(strcmp(al->type[k],maps->types[j])==0){
        mapid = j;
      }
    }
    if(mapid == -1) {
      fprintf(stderr,"Error: No map for atomtype %s\n",al->type[k]);
      exit(1);
    }
    xpts = gridpoint(al->x[k][XX],maps->map[mapid]->origin[XX],inv_spacing,nx);
    ypts = gridpoint(al->x[k][YY],maps->map[mapid]->origin[YY],inv_spacing,ny);
    zpts = gridpoint(al->x[k][ZZ],maps->map[mapid]->origin[ZZ],inv_spacing,nz);
    idx = zpts*nx*ny + ypts*nx + xpts;
    en += maps->map[mapid]->values[idx]; 
    en += maps->map[ELEC]->values[idx]*al->q[k]; 
    en += maps->map[DESOLV]->values[idx]*fabs(al->q[k]); 
  }
  return en;
}

real eval_energies(t_atomlist *al, t_maps *maps, 
                   rvec **coord_set, int nset,
                   real *ener, int nbest, int *index)
{
  
  real best = 10000.;
  int best_id = -1;
  int i,j,k;
  int mapid;
  real en;
  int xpts, ypts, zpts, idx;
  int nx,ny,nz;
  real inv_spacing;
  real energies[nset];
  int ELEC = -1;
  int DESOLV = -1;
  
  for(i=0;i<nset;i++){
    energies[i] = 0.;
  }
  inv_spacing = maps->map[0]->inv_spacing;
  nx = maps->map[0]->n[XX]; 
  ny = maps->map[0]->n[YY]; 
  nz = maps->map[0]->n[ZZ]; 


  for(i=0;i<maps->n;i++){
    if(strcmp(maps->types[i],"e") == 0) ELEC = i;
    else if(strcmp(maps->types[i],"d") == 0) DESOLV = i;
  }
  

  
  for(i=0;i<nset;i++){
    en = 0;
    for(k=0;k<al->natoms;k++){  
      /* energy per atom */
      mapid = -1;
      for(j=0;j<maps->n;j++){
        if(strcmp(al->type[k],maps->types[j])==0){
          mapid = j;
        }
      }
      if(mapid == -1) {
        fprintf(stderr,"Error: No map for atomtype %s\n",al->type[k]);
        exit(1);
      }
      xpts = gridpoint(coord_set[i][k][XX],maps->map[mapid]->origin[XX],inv_spacing,nx);
      ypts = gridpoint(coord_set[i][k][YY],maps->map[mapid]->origin[YY],inv_spacing,ny);
      zpts = gridpoint(coord_set[i][k][ZZ],maps->map[mapid]->origin[ZZ],inv_spacing,nz);
      idx = zpts*nx*ny + ypts*nx + xpts;
      en += maps->map[mapid]->values[idx]; 
      en += maps->map[ELEC]->values[idx]*al->q[k]; 
      en += maps->map[DESOLV]->values[idx]*fabs(al->q[k]); 
      }
    energies[i] = en;
    if(energies[i] < best ) {
      best_id = i;
      best = energies[i];
    }
  }
  get_lowest_energies(energies,nset,ener,nbest,index);
  
  /* write the best conformation to atomlist */
  coord_set2al(coord_set[best_id],al);   
  return best;
}




int main(int argc, char **argv)
{
  
   static char *desc[] = {
    "Testprogram",
    "........"
  };

   t_atomlist *al = atomlist_init();
  int i,j,k;
  t_idxgroups *don = idx_init();
  t_idxgroups *acc = idx_init();
  t_idxgroups *phob = idx_init();
  t_idxgroups *pl = idx_init();
  t_idxgroups *imp = idx_init();
  t_bounds *b = NULL;
  bool bVerbose = FALSE;
  t_boundtrack *bt = NULL;
  gmx_rng_t rng;
  rvec co;
  int seed=-1;
  int xpts,ypts,zpts; 
  int nx,ny,nz; 
  real en, scalar; 
  char **fnms,*in_file;
  int nfile_in;
  t_maps *maps;
  char *map_ext1=NULL,*map_ext2=NULL;
  int nelem;
  rvec **coord_set = NULL;
  rvec **new_set = NULL;
  int nset = 25;
  int nconf = 1000;
  int ncopies = 0;
  
  int nbest = 10;
  real best_ene[10];
  int index[10];
  

  
  time_t tstart;


  t_filenm fnm[] = {
    { efPDBQT,"-p","receptor", ffREAD},
    { efGRID,"-map",NULL, ffRDMULT}, 
    { efCTP,"-top","ligand", ffREAD}, 
    { efDAT,"-d","tdist", ffREAD}, 
    { efPDB,"-o","pot.pdb", ffWRITE}
  };
  
  
  t_pargs pa[] = {
    {"-v",   FALSE, etBOOL, {&bVerbose},"make noise"} ,
    {"-nset",   FALSE, etINT, {&nset},"Number of coordinate sets"}, 
    {"-nconf",   FALSE, etINT, {&nconf},"Number of conformations"} 
/*     {"-energy",   FALSE, etREAL, {&energy},"energy"}, */
/*     {"-minat",   FALSE, etINT, {&minat},"min # of atoms"}, */
/*     {"-cutoff",   FALSE, etREAL, {&cutoff},"cutoff"} */
    
  };
  
#define NFILE asize(fnm)
  
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if(seed == -1)
    seed = gmx_rng_make_seed();
  rng = gmx_rng_init(seed);
  fprintf(stderr,"Initial random seed = %d\n",seed);


  
  time(&tstart);
  snew(maps,1);

  nfile_in = opt2fns(&fnms,"-map",NFILE,fnm);
  snew(maps->map,nfile_in);
  snew(maps->types,nfile_in);
  
  ncopies = nset*nconf;
  

  char xx[STRLEN];
  char *tmp;
  
  for(i=0;i<nfile_in;i++){
    printf("Reading Map %s\n",fnms[i]);
    snew(maps->types[i],10);
    map_ext1 = strchr(fnms[i],'.');
    strcpy(xx,map_ext1+1);
    map_ext2 = strchr(xx,'.');
    tmp = strdup(xx);
    tmp[map_ext2-xx] = '\0';
    strcpy(maps->types[i],tmp);
    maps->map[i] = read_map(fnms[i],NULL,TRUE); 
    maps->n++;
  }
  
  t_atomlist *qal = read_pdbqt(opt2fn("-p",NFILE,fnm),NULL,FALSE); 
  
  nelem = maps->map[0]->nelem;
  for(i=0;i<maps->n;i++){
    if(maps->map[i]->nelem!=nelem){
      fprintf(stderr,"Error, Maps are not consistent\n");
      exit(1);
    }
  }

  for(i=0;i<maps->n;i++){
    char fn[STRLEN];
    sprintf(fn,"map_%s.dx",maps->types[i]);
    fprintf(stderr,"Writing map %s\n",fn);
    write_dx_map(fn,maps->map[i]);
  }
  

  
  nx = maps->map[0]->n[XX]; 
  ny = maps->map[0]->n[YY]; 
  nz = maps->map[0]->n[ZZ]; 


  real ener;
   
  read_cnctop(opt2fn("-top",NFILE,fnm),al,acc,don,phob,pl,imp);
  get_bconstr(al);  
  occ_to_one(al);
  
  b = read_bounds(opt2fn("-d",NFILE,fnm),bVerbose);
  bt = boundtrack_init();
  bt = boundtrack_realloc(bt,al->natoms);
  fill_boundtrack(b,bt);


  /* make neighborlist */
  
  for(i=0;i<al->natoms;i++){
    strcpy(al->type[i],qal->type[i]);
    snew(al->nb[i],MAXNB);
    al->nnb[i]=0;
  }
  for(i=0;i<al->natoms-1;i++){
    for(k=i+1;k<al->natoms;k++){
      if(!al_connected(al,i,k)){
        al->nb[i][al->nnb[i]]=k;
        al->nb[k][al->nnb[k]]=i;
        al->nnb[i]+=1;
        al->nnb[k]+=1;
      }
    }
  }


  /* generate ensemble of conformations */

  snew(coord_set,ncopies);
  snew(new_set,ncopies);
  i = 0;
  int ntot=0;
/*   printf("Generating Initial Conformers....\n"); */
  while (i<nconf){
    random_al_simple(al,rng);
    int succ = ligand_disco(al,b,imp,pl,500,rng,stderr,bt,TRUE);
    if(succ) {
      com(al);
      for(j=0;j<nset;j++){
        snew(coord_set[ntot],al->natoms);        
        snew(new_set[ntot],al->natoms);        
        random_rot(al->x,al->natoms,rng); 
        for(k=0;k<al->natoms;k++){
          copy_rvec(al->x[k],coord_set[ntot][k]);
        }
        ntot++;
      }
      i++;
       printf("Generating Conformers %d\r",i); 
    }
  }
  printf("\n");


  
  real xmin,xmax, xdim;
  real ymin,ymax, ydim;
  real zmin,zmax, zdim;
  real rx,ry,rz;
  
  xmin = maps->map[0]->origin[XX];
  ymin = maps->map[0]->origin[YY];
  zmin = maps->map[0]->origin[ZZ];
  
  /* calculate box maxima */
  xmax = xmin+maps->map[0]->n[XX]*maps->map[0]->spacing;
  ymax = ymin+maps->map[0]->n[YY]*maps->map[0]->spacing;
  zmax = zmin+maps->map[0]->n[ZZ]*maps->map[0]->spacing;
  


  int runs = 0;
  real bbest = 10000.;
  
  FILE *of = ffopen("traj.pdb","w");
  
  do
  {
    /* outer loop. flood binding site and proceed
       with the best scored conformation
    */

    copy_coord_set(coord_set,new_set,ncopies,al->natoms);
    flood_binding_site(al,new_set,ncopies,rng,
                       xmin,xmax,
                       ymin,ymax,
                       zmin,zmax);
    
    
    
    ener = eval_energies(al,maps,new_set,ncopies,best_ene,nbest,index);
    printf("Lowest energy %g\n",ener);
    
    /* now we have the coords of the best scored conformation
       in the atomlist. Let's disco again around this 
       conformation
    */
    
    rvec *conf;
    snew(conf,al->natoms);
    al2rvec(conf,al);
    
    int nround = 1;
    i = 0;

    real best = 10000.;
    real pert = 16.5;
    
    do{
      
      /* inner loop. start with the best scored conformation    
         and try to optimize
      */
      int ntry = 0;
      
      while(i<1000)
      {
        rvec2al(al,conf);
        perturb_atomlist(al,rng,2.);
        int succ = ligand_disco(al,b,imp,pl,500,rng,stderr,bt,TRUE);
        if(succ){
          ntry++;
          put_molecule_in_box(al,al->x,
                              xmin,xmax,
                              ymin,ymax,
                              zmin,zmax);
          real e = calc_energy(al,maps); 
          if(e < ener){
            printf("energy = %g\n",e);
            al2rvec(conf,al);
            write_pdb_frame(of,al,i);
            i++;ntry = 0;
            ener = e;
            if(e < bbest){
              for(k=0;k<al->natoms;k++){
                copy_rvec(al->x[k],qal->x[k]);
              }
            }
            
          }
        }
        if(ntry%100 == 0) printf("ntry = %d\n",ntry);
        
        if(ntry > 1000) break;
        
      }
      break;
      
      



/*       pert = 15./(real) nround;  */
/*       while (i<500){ */
/*         rvec2al(al,conf); */
/*         perturb_atomlist(al,rng,pert); */
/*         int succ = ligand_disco(al,b,imp,pl,500,rng,stderr,bt,TRUE); */
/*         if(succ) { */
/*           for(k=0;k<al->natoms;k++){ */
/*             copy_rvec(al->x[k],new_set[i][k]); */
/*           } */
/*           put_molecule_in_box(al,new_set[i], */
/*                               xmin,xmax, */
/*                               ymin,ymax, */
/*                               zmin,zmax); */
/*           i++; */
/*         } */
/*       } */
      
/*       ener = eval_energies(al,maps,new_set,500,best_ene,nbest,index); */
/*       if(ener < best) { */
/*         best = ener; */
/*         al2rvec(conf,al); */
/*       } */
/*       if(ener < bbest){ */
/*         bbest = ener; */
/*         for(j=0;j<al->natoms;j++){ */
/*           copy_rvec(al->x[j],qal->x[j]); */
/*         } */
/*       } */
      
/*       printf("Lowest energy %g after round %d (run = %d (%g))\n",ener,nround,runs,pert); */
/*       write_pdb_frame(of,al,nround); */
/*       nround++;i = 0; */
      

      

    }while(nround < 16);
    runs++;
    sfree(conf);
    
    
  }while (runs < 20);



  printf("Lowest energy of all runs = %g\n",bbest);
  

  write_pdb(qal,"best.pdb");  


  time_t tend;
  time(&tend);
  char timestr[STRLEN];
  get_time(difftime(tend,tstart),timestr);
  
  fprintf(stderr,"%s\n",timestr);

  
   return 0;
   
}
