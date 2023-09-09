#include <tconcoord.h>

/*============================================================*/
t_histogram *new_histogram(real lb, real ub, real delta)
{
  real size;
  int isize;
  int i;
  
  t_histogram *h = NULL;
  snew(h,1);
  
  /* calculate array size */
  
  size = (ub-lb)/delta;
  isize = (int) size + 1;
  
  snew(h->x,isize);
  snew(h->y,isize);
  h->n = isize;
  h->delta = delta;
  h->inv_delta = 1./delta;
  
  /* init arrays */

  for(i=0;i<h->n;i++){
    h->x[i] = lb+delta*i;
    h->y[i] = 0.;
  }
  
  return h;
}
/*============================================================*/
t_histogram *histogram_init(void)
{
  t_histogram *h = NULL;
  snew(h,1);
  snew(h->x,1);
  snew(h->y,1);
  h->delta = 0.;
  h->inv_delta = 0;
  h->n=0;
  return h;
}

/*============================================================*/
t_histogram *histogram_realloc(t_histogram *h, int n)
{
  srenew(h->x,n);
  srenew(h->y,n);
  return h;
}

/*============================================================*/
void free_histogram(t_histogram *h)
{
  sfree(h->x);
  sfree(h->y);
  sfree(h);
}

/*============================================================*/

void add_value(t_histogram *h, real val, bool flag, real weight)
{
  int i;
  
  for(i=0;i<h->n-1;i++){
    if(val >= h->x[i] && val < h->x[i+1]){
      if(flag)
        h->y[i]+=weight;
      else
        h->y[i]+=1;
      break;
    }
  }
}

/*============================================================*/
void print_histogram(FILE *fp,t_histogram *h)
{
  int i;
  
  for(i=0;i<h->n;i++){
    fprintf(fp,"%8.4f %8.4f\n",h->x[i],h->y[i]);
  }
}

/*============================================================*/
void norm_histogram(t_histogram *h)
{
  int i;

  real sum = integrate(h);
    
  for(i=0;i<h->n;i++){
    h->y[i]/=sum;
  }

}
/*============================================================*/
real integrate(t_histogram *h)
{
  int i;
  real sum = 0.;
  
  for(i=0;i<h->n;i++){
    sum += h->delta*h->y[i];
  }
  
  return sum;
}
/*============================================================*/

t_histogram *read_histogram(FILE *fp)
{
  char line[STRLEN];
  t_histogram *h = histogram_init();
  
    
  while(get_a_line(fp,line,STRLEN)){
    h->n++;
    h = histogram_realloc(h,h->n);
#ifdef GMX_DOUBLE
    sscanf(line,"%lf %lf",&h->x[h->n-1],&h->y[h->n-1]);
#else
    sscanf(line,"%f %f",&h->x[h->n-1],&h->y[h->n-1]);
#endif
  }
  h->delta = h->x[1]-h->x[0];
  h->inv_delta = 1./h->delta;
  
  return h;
}


