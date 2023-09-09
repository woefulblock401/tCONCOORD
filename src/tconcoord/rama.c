#include <tconcoord.h>

real val_from_2d_table(t_gridmap *gp, real v1, real v2)
{
  int gpx, gpy, idx;
  gpx = gridp(v1,gp->origin[XX],gp->inv_spacing,gp->n[XX]);
  gpy = gridp(v2,gp->origin[YY],gp->inv_spacing,gp->n[YY]);
  idx = gpx*gp->n[YY] + gpy;
  return gp->values[idx];
}


t_gridmap *read_rama(void)
{
  
  FILE *fp = cnclib("rama.dat");
  int i;
  char line[STRLEN];
  int xdim, ydim, zdim;
  int nx, ny, nz;
  real phi, psi, val;
  t_gridmap *gp = gridmap_init();
  int idx;
  int gpx, gpy, gpz;
  
  gp->spacing = 7.2;
  gp->inv_spacing = 1./gp->spacing;
  for(i=0;i<DIM;i++){
    gp->origin[i] = -180.;
  }
  
  xdim = (int) (360)*gp->inv_spacing ;
  ydim = (int) (360)*gp->inv_spacing ;

  gp->n[XX] = xdim ;
  gp->n[YY] = ydim ;

  nx = gp->n[XX];
  ny = gp->n[YY];

  gp->nelem = gp->n[XX]*gp->n[YY];
  srenew(gp->values,gp->nelem);
  printf("nx = %d, ny = %d\n",nx,ny);
  idx = 0;
  
  while(get_a_line(fp,line,STRLEN))
  {
    sscanf(line,"%lf %lf %lf", &phi, &psi, &val);
/*     gpx = gridp(phi,gp->origin[XX],gp->inv_spacing,xdim); */
/*     gpy = gridp(psi,gp->origin[YY],gp->inv_spacing,ydim); */
/*     idx = gpx*ny + gpy; */
    printf("phi = %g psi = %g gp = %d\n",phi, psi, idx); 
    gp->values[idx] = val;
    idx++;
  }
  return gp;
}

