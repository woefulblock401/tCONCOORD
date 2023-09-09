#include<tconcoord.h>


t_scdata *scdata_init(void)
{
  t_scdata *data = NULL;
  snew(data,1);
  data->n = 0;
  data->resname = NULL;
  data->name1 = NULL;
  data->name2 = NULL;
  data->lb = NULL;
  data->ub = NULL;
  return data;
}

t_scdata *scdata_realloc(t_scdata *data, int n)
{
  srenew(data->resname,n);
  srenew(data->name1,n);
  srenew(data->name2,n);
  srenew(data->lb,n);
  srenew(data->ub,n);
  return data;
}

t_scdata *read_scdata(void)
{
  
  FILE *fp = cnclib("Sidechain.dat");
  char line[STRLEN];
  shstr resname,name1,name2;
  real lb,ub;
  t_scdata *data = scdata_init();
  int i;
  
  while(get_a_line(fp,line,STRLEN)){
#ifdef GMX_DOUBLE
    sscanf(line,"%s %s %s %lf %lf",resname,name1,name2,&lb,&ub);
#else
    sscanf(line,"%s %s %s %f %f",resname,name1,name2,&lb,&ub);
#endif
    extend_name(name1);
    extend_name(name2);
    data->n++;
    data = scdata_realloc(data,data->n);
    i = data->n-1;
    strcpy(data->resname[i],resname);
    strcpy(data->name1[i],name1);
    strcpy(data->name2[i],name2);
    data->lb[i] = lb;
    data->ub[i] = ub;
  }
  return data;
}

