#include <tconcoord.h>
/*============================================*/
t_namegroups *namegroups_init(void)
{
  t_namegroups *ng = NULL;
  snew(ng,1);
  ng->n = 0;
  ng->natoms = NULL;
  ng->atomnames = NULL;
  ng->resname = NULL;
  ng->val = NULL;
  return ng;
}
/*============================================*/
t_namegroups *namegroups_realloc(t_namegroups *ng, int n)
{
  snew(ng->natoms,n);
  snew(ng->atomnames,n);
  snew(ng->resname,n);
  snew(ng->val,n);
  return ng;
}
/*============================================*/
t_namegroups *read_namegroups(FILE *fp)
{
  char line[STRLEN];
  t_namegroups *ng = namegroups_init();
  int max = 2000;
  ng = namegroups_realloc(ng,max); /* dirty */
  int n,i,k;
  char dummy[STRLEN];
  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      sscanf(line+1,"%s",dummy);
      ng->n++;
      n = ng->n-1;
      strcpy(ng->resname[n],dummy);
      get_a_line(fp,line,STRLEN);
#ifdef GMX_DOUBLE
      sscanf(line,"%d %lf", &ng->natoms[n],
             &ng->val[n]);
#else
      sscanf(line,"%d %f", &ng->natoms[n],
             &ng->val[n]);
#endif
      int nat = ng->natoms[n];
      srenew(ng->atomnames[n],nat);
      for(i=0;i<nat;i++){
        get_a_line(fp,line,STRLEN);
        sscanf(line,"%s",ng->atomnames[n][i]); 
      }
    }
  }
  
  for(i=0;i<ng->n;i++){
/*     printf("resn: %s %g\n",ng->resname[i], ng->val[i]);  */
    for(k=0;k<ng->natoms[i];k++){
      extend_name(ng->atomnames[i][k]);
/*       printf("'%s' \n",ng->atomnames[i][k]);  */
    }
  }
  
  return ng;
}
/*============================================*/

