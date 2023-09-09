#include <tconcoord.h>

/* static int NVEC = 30; */

/* static rvec vectors[] = { */
/*   /\* Vectors for octahedral triangulation *\/ */
/*   {  0.00000000,   0.00000000,   1.00000000},  */
/*   {  0.00000000,   0.00000000,  -1.00000000}, */
/*   {  0.70710678,   0.70710678,   0.00000000}, */
/*   {  0.70710678,  -0.70710678,   0.00000000}, */
/*   { -0.70710678,   0.70710678,   0.00000000}, */
/*   { -0.70710678,  -0.70710678,   0.00000000}, */
/*   {  0.86602540,   0.28867513,   0.40824829}, */
/*   {  0.86602540,  -0.28867513,   0.40824829}, */
/*   {  0.57735027,   0.00000000,   0.81649658}, */
/*   {  0.86602540,   0.28867513,  -0.40824829}, */
/*   {  0.86602540,  -0.28867513,  -0.40824829}, */
/*   {  0.57735027,   0.00000000,  -0.81649658}, */
/*   {  0.28867513,  -0.86602540,   0.40824829}, */
/*   { -0.28867513,  -0.86602540,   0.40824829}, */
/*   {  0.00000000,  -0.57735027,   0.81649658}, */
/*   {  0.28867513,  -0.86602540,  -0.40824829}, */
/*   { -0.28867513,  -0.86602540,  -0.40824829}, */
/*   {  0.00000000,  -0.57735027,  -0.81649658}, */
/*   { -0.86602540,   0.28867513,   0.40824829}, */
/*   { -0.86602540,  -0.28867513,   0.40824829}, */
/*   { -0.57735027,   0.00000000,   0.81649658}, */
/*   { -0.86602540,   0.28867513,  -0.40824829}, */
/*   { -0.86602540,  -0.28867513,  -0.40824829}, */
/*   { -0.57735027,   0.00000000,  -0.81649658}, */
/*   { -0.28867513,   0.86602540,   0.40824829}, */
/*   {  0.28867513,   0.86602540,   0.40824829}, */
/*   {  0.00000000,   0.57735027,   0.81649658}, */
/*   { -0.28867513,   0.86602540,  -0.40824829}, */
/*   {  0.28867513,   0.86602540,  -0.40824829}, */
/*   {  0.00000000,   0.57735027,  -0.81649658} */
/* }; */


void sort_clusters(t_idxgroups *grp)

/* sort clusters by score */

{
  int i;
  int idum;
  real rdum;
  int *ptr;
  bool changed = TRUE;
  while(changed){
    changed = FALSE;
    for(i=0;i<grp->n-1;i++){
      if(grp->val[i]<grp->val[i+1]){
        idum = grp->natoms[i];
        rdum = grp->val[i];
        grp->val[i] = grp->val[i+1];
        snew(ptr,grp->natoms[i]);
        copy_iarr(ptr,grp->atoms[i],grp->natoms[i]);
        sfree(grp->atoms[i]);
        snew(grp->atoms[i],grp->natoms[i+1]);
        copy_iarr(grp->atoms[i],grp->atoms[i+1],grp->natoms[i+1]);
        grp->natoms[i] = grp->natoms[i+1];
        grp->natoms[i+1] = idum;
        grp->val[i+1] = rdum;
        sfree(grp->atoms[i+1]);
        snew(grp->atoms[i+1],idum);
        copy_iarr(grp->atoms[i+1],ptr,idum);
        changed = TRUE;
      }
    }
  }
}

void max_coord(t_atomlist *al, matrix max)
{
  int i,k;
  for(i=0;i<DIM;i++){
    max[i][XX] = 1000;
    max[i][YY] = -1000.;
  }
  for(i=0;i<al->natoms;i++){
    if(al->x[i][XX] < max[XX][XX]) max[XX][XX] = al->x[i][XX];
    if(al->x[i][XX] > max[XX][YY]) max[XX][YY] = al->x[i][XX];
    if(al->x[i][YY] < max[YY][XX]) max[YY][XX] = al->x[i][YY];
    if(al->x[i][YY] > max[YY][YY]) max[YY][YY] = al->x[i][YY];
    if(al->x[i][ZZ] < max[ZZ][XX]) max[ZZ][XX] = al->x[i][ZZ];
    if(al->x[i][ZZ] > max[ZZ][YY]) max[ZZ][YY] = al->x[i][ZZ];
  }
  max[XX][ZZ] = max[XX][YY]-max[XX][XX];
  max[YY][ZZ] = max[YY][YY]-max[YY][XX];
  max[ZZ][ZZ] = max[ZZ][YY]-max[ZZ][XX];
}

/* bool find_atom(t_atomlist *al, rvec x, real cut2) */
/* { */

/*   rvec diff; */
/*   int i; */
  
/*   for(i=0;i<al->natoms;i++){ */
/*     if (strcmp(al->symbol[i],"H") != 0){ */
      
/*       if(sqr(al->x[i][XX]-x[XX]) < cut2 && */
/*          sqr(al->x[i][YY]-x[YY]) < cut2 && */
/*          sqr(al->x[i][ZZ]-x[ZZ]) < cut2){ */
/*         rvec_sub(al->x[i],x,diff); */
/*         if(norm2(diff) < cut2) return TRUE; */
/*       } */
/*     } */
/*   } */
/*   return FALSE; */
/* } */

/* real ray_search(t_atomlist *al,rvec x) */
/* { */
/*   /\* search for atoms along the ray vectors *\/ */
  
/*   real cutoff = .3; */
/*   real cut2;  */
/*   int i,k; */
/*   rvec probe; */
/*   real count = 0; */
/* /\*   FILE *fp = ffopen("tmp.pdb","w"); *\/ */
  
/*   for(i=0;i<NVEC;i++){ */
/*     copy_rvec(x,probe); */
/*     cutoff = .3; */
/*     cut2 = sqr(cutoff); */
    
/* /\*     fprintf(fp,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n",  *\/ */
/* /\*             count," DUM"," ","DUM"," ",1,probe[XX],probe[YY],probe[ZZ],1.00,1.00);  *\/ */
/*     if(find_atom(al,probe,cut2)) count+=1.; */
/*     else { */
/*       /\* move along the ray *\/ */
/*       for(k=1;k<11;k++){ */
/*         cutoff+=.3; */
/*         cut2 =  sqr(cutoff); */
/*         probe[XX]+=vectors[i][XX]; */
/*         probe[YY]+=vectors[i][YY]; */
/*         probe[ZZ]+=vectors[i][ZZ]; */
/* /\*         fprintf(fp,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n",  *\/ */
/* /\*                 count," DUM"," ","DUM"," ",1,probe[XX],probe[YY],probe[ZZ],1.00,1.00);  *\/ */

/*         if(find_atom(al,probe,cut2)) */
/*         { */
/*           count+=1./(real) k; */
/*           break; /\* take next vector *\/ */
/*         } */
/*       } */
/*     } */
/*   } */

  
/*   return count; */
/* } */

t_idxgroups *cluster_grid(rvec *stored, int *idx, int size, real min, real *contact)

{
  /* the coordinates are stored in 'stored' */
  /* idx contains the remaining grid points */

  printf("Beginning cluster analysis........\n");
  

  t_idxgroups *cl = idx_init();
  
  int i,k;
  rvec diff;
  real d, cut;
  bool clustered[size];
  for(i=0;i<size;i++) clustered[i] = FALSE;
  cut = sqr(min);
  
  for(i=0;i<size-1;i++){
    cl->n++;
    cl = idx_realloc(cl,cl->n);
    for(k=i+1;k<size;k++){
      rvec_sub(stored[idx[i]],stored[idx[k]],diff);
      d = norm2(diff);
      if(d < cut){
        add_to_group(cl,i,idx[k]); 
      }
    }
  }
  cl = remove_redundant_groups(cl);
  cl = merge_groups(cl,3,TRUE);
  cl = remove_redundant_groups(cl);
  real av = 0.;
  int id;
  for(i=0;i<cl->n;i++){
    for(k=0;k<cl->natoms[i];k++){
      id = cl->atoms[i][k];
      av+=contact[id];
    }
    av/=(real) cl->natoms[i];
    cl->val[i] = av;
  }
  
  sort_clusters(cl);
  return cl;
}

void print_clusters(char *filename, t_idxgroups *cl, rvec *stored, 
                    real *contact, int clsize)
{
  FILE *fp=ffopen(filename,"w");
  int i,k;
  int idx;
  real av;
  int count;
  
  for(i=0;i<cl->n;i++){
    if(cl->natoms[i] >= clsize) 
    {
      
 /*      av = 0.; */
/*       for(k=0;k<cl->natoms[i];k++){ */
/*         idx = cl->atoms[i][k]; */
/*         av+=contact[idx]; */
/*       } */
/*       av/=(real) cl->natoms[i]; */
/*       cl->val[i] = av; */
      
      fprintf(fp,"MODEL%10d\n",i+1);
      fprintf(fp,"REMARK   CLUSTERSCORE = %8.3f\n",cl->val[i]);
      fprintf(fp,"REMARK   CLUSTERSIZE  = %8d\n",cl->natoms[i]);
      count = 1;
      for(k=0;k<cl->natoms[i];k++){
        idx = cl->atoms[i][k];
        fprintf(fp,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n", 
                count," DUM"," ","DUM"," ",1,stored[idx][XX],stored[idx][YY],stored[idx][ZZ],1.00,(real) contact[idx]); 
        count++;
      }
      fprintf(fp,"ENDMDL\n");
    }
  }
}





int main(int argc, char **argv)
{
  
   static char *desc[] = {
    "Find binding site",
    ".....hopefully..."
  };

  t_atomlist *al;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb;
  t_vdwcomb *vdw14comb;
  t_types *tp;
  t_forcerec *fr = NULL;
  t_topology *top = NULL; 
  t_atoms *atoms = NULL; 
  
  
  char title[STRLEN];
  matrix box;
  rvec *x;
  int i,k,j;
  real rad = 1.55;
  real eps = 0.1591;
  real grid = 1.;
  real xx,yy,zz;
  real cutoff = 6.;
  real energy = -1.;
  int maxorder = 0;
  t_area *ta = NULL;
  int minat = 20;
  int clsize = 20.;
  


  t_filenm fnm[] = {
    { efTPS,"-s","topol", ffREAD},
/*     { efPDB,"-pot","pot.pdb", ffWRITE},   */
/*     { efPDB,"-o","filtered.pdb", ffWRITE},  */
    { efPDB,"-cl","clusters.pdb", ffWRITE}
  };
  
  
  t_pargs pa[] = {
    {"-r",   FALSE, etREAL, {&rad},"probe radius"},
    {"-eps",   FALSE, etREAL, {&eps},"probe epsilon"},
    {"-grid",   FALSE, etREAL, {&grid},"grid size"},
    {"-energy",   FALSE, etREAL, {&energy},"energy"},
    {"-minat",   FALSE, etINT, {&minat},"Minimum number of atoms"},
    {"-clsize",   FALSE, etINT, {&clsize},"Minimum cluster size"},
    {"-maxorder",   FALSE, etINT, {&maxorder},"Cut sidechains"},
    {"-cutoff",   FALSE, etREAL, {&cutoff},"cutoff radius"}
    
  };
  
#define NFILE asize(fnm)
  
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  
  snew(top,1);
  read_tps_conf(opt2fn("-s",NFILE,fnm),title,top,&x,NULL,
                box,TRUE);
/*   FILE *fp = ffopen(opt2fn("-o",NFILE,fnm),"w");  */
/*   FILE *fpot = ffopen(opt2fn("-pot",NFILE,fnm),"w");  */

  al = al_from_atoms(&(top->atoms),x);
  get_symbol(stderr,al);
  get_order(al);

  t_contab *ct = contab_init();
  ct=contab_realloc(ct,al->natoms);
  nb_search(al,ct,0.8,TRUE);
  
  rename_at(al);
  get_rosetta_types(al);
  
  

  matrix max;
  max_coord(al,max);
  snew(ta,1);
  
  ta->xmin = max[XX][XX];
  ta->xmax = max[XX][YY];
  ta->ymin = max[YY][XX];
  ta->ymax = max[YY][YY];
  ta->zmin = max[ZZ][XX];
  ta->zmax = max[ZZ][YY];
  

  real en;
  rvec coord;
  int count = 0;

  rvec *stored;
  snew(stored,1);
    
  for(zz = ta->zmin;zz<ta->zmax;zz+=grid){
    for(yy = ta->ymin;yy<ta->ymax;yy+=grid){
      for(xx = ta->xmin;xx<ta->xmax;xx+=grid){
        coord[XX]=xx;
        coord[YY]=yy;
        coord[ZZ]=zz;
        en = lj_lattice_energy(al,coord,cutoff,eps,rad,maxorder);
        if(en < energy) {
          copy_rvec(coord,stored[count]);
/*           fprintf(fpot,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n",    */
/*                   count," DUM"," ","DUM"," ",1,xx,yy,zz,1.00,en);    */
          count++;
          srenew(stored,count+1);
        }
      }
    }
  }
  
  printf("Evaluating %d grid points for buriedness \n",count);
    
  real contact[count];
  for(i=0;i<count;i++) contact[i] = 0.;
  int *idx;
  snew(idx,1);
  int idxcount = 0;
  
  for(i=0;i<count;i++){
    if(i%1000 == 0) printf("Number of Grid Points Completed...: %d\n",i);
    contact[i] = ray_search(al,stored[i],NULL,0);
    if(contact[i] > minat) {
/*        fprintf(fp,"ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n",   */
/*                count," DUM"," ","DUM"," ",1,stored[i][XX],stored[i][YY],stored[i][ZZ],1.00,(real) contact[i]);   */
      idxcount++;
      srenew(idx,idxcount);
      idx[idxcount-1] = i;
    }
  }
  
  printf("Using %d grid points for cluster calculation\n",idxcount);
  
  t_idxgroups *cluster = cluster_grid(stored,idx,idxcount,2.5, contact);
  print_clusters(opt2fn("-cl",NFILE,fnm),cluster,stored,contact,clsize);
  return 0;
  
}
