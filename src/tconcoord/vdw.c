#include <tconcoord.h>
/*============================================*/
t_vdw *vdw_init(void)
{
  t_vdw *vdw = NULL;
  snew(vdw,1);
  /* yes.... not dynamically again */
  snew(vdw->type,1000);
  snew(vdw->vdw,1000);
  snew(vdw->vdw14,1000);
  snew(vdw->rs_type,1000);
  snew(vdw->rs_rad,1000);
  snew(vdw->rs_G,1000);
  snew(vdw->rs_Gref,1000);
  snew(vdw->rs_V,1000);
  snew(vdw->rs_lambda,1000);
  snew(vdw->rs_eps,1000);
  vdw->n = 0;
  return vdw;
}
/*============================================*/
void free_vdw(t_vdw *vdw)
{
  sfree(vdw->type);
  sfree(vdw->vdw);
  sfree(vdw->vdw14);
  sfree(vdw->rs_type);
  sfree(vdw->rs_rad);
  sfree(vdw->rs_G);
  sfree(vdw->rs_Gref);
  sfree(vdw->rs_V);
  sfree(vdw->rs_lambda);
  sfree(vdw);
}
/*============================================*/

t_vdw *read_vdw_radii(int num)
{
  t_vdw *vdw = vdw_init();
  int i=0;
  char line[STRLEN];
  FILE *fp = NULL;

  switch(num){
    case(1):
      fp = cnclib("Atomradii.dat");
      break;
    case(2):
      fp = cnclib("Atomradii_poh.dat");
      break;
    case(3):
      fp = cnclib("Atomradii_noh.dat");
      break;
    case(4):
      fp = cnclib("Atomradii_oplsaa.dat");
      break;
    case(5):
      fp = cnclib("Atomradii_yamber2.dat");
      break;
    case(6):
      fp = cnclib("Atomradii_li.dat");
      break;
  }
  
  if(!find_key_word(fp,"[ TYPES ]")){
    fatal_error("Shitty file!\n");
  }
  
  while(get_a_line(fp,line,STRLEN) &&
        strchr(line,'[')==NULL)
  {
    cut_string(line,'#');
#ifdef GMX_DOUBLE
    if(sscanf(line,"%s %lf %lf",vdw->type[i],
              &vdw->vdw[i],&vdw->vdw14[i])==3){
#else
    if(sscanf(line,"%s %f %f",vdw->type[i],
              &vdw->vdw[i],&vdw->vdw14[i])==3){
#endif
      i++;
    }
  }
  vdw->n = i;
  fclose(fp);
  
  return vdw;
}

