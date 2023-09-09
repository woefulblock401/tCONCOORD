#include <tconcoord.h>
/*============================================*/

t_types *types_init(void)
{
  t_types *tp = NULL;
  snew(tp,1);
  /* do it lazy 
     who needs more ? 
  */
  snew(tp->resname,1000);
  snew(tp->name,1000);
  snew(tp->type,1000);
  snew(tp->hyb,1000);
  snew(tp->rs_type,1000);
  tp->n = 0;
  return tp;

}
/*============================================*/

void free_types(t_types *tp)
{
  sfree(tp->resname);
  sfree(tp->name);
  sfree(tp->type);
  sfree(tp->hyb);
  sfree(tp->rs_type);
  sfree(tp);
}
/*============================================*/

t_types *read_atom_types(int num)
{
  int i = 0;
  t_types *tp = types_init();
  char line[STRLEN];
  char dummy[STRLEN];
  FILE *fp = NULL;

/*   if(num==1) */
    fp = cnclib("Atomtypes.dat");
/*   else */
/*     fp = cnclib("Atomtypes_old.dat"); */

  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      substring(dummy,line,'[',']');
      sscanf(dummy+1,"%s",dummy);
    }
    else{
      cut_string(line,'#');
      sscanf(line,"%s %s %s",tp->name[i],tp->type[i],tp->hyb[i]);
      strcpy(tp->resname[i],dummy);
      if(strcmp(tp->resname[i],"GENERIC")!=0)
        extend_name(tp->name[i]);
      i++;
    }
  }
  tp->n = i;
  fclose(fp);
  
  return tp;
}
/*============================================*/

