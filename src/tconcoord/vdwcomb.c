#include <tconcoord.h>

/*============================================*/

t_vdwcomb *vdwcomb_init(void)
{
  t_vdwcomb *vdw = NULL;
  snew(vdw,1);
  vdw->n = 0;
  vdw->type1 = NULL;
  vdw->type2 = NULL;
  vdw->vdw = NULL;
  vdw->rs_type1 = NULL;
  vdw->rs_type1 = NULL;
  vdw->rs_comb = NULL;
  return vdw;
}
/*============================================*/

t_vdwcomb *vdwcomb_realloc(t_vdwcomb *vdw, int n)
{
  srenew(vdw->type1,n);
  srenew(vdw->type2,n);
  srenew(vdw->rs_type1,n);
  srenew(vdw->rs_type2,n);
  srenew(vdw->rs_comb,n);
  srenew(vdw->vdw,n);
  return vdw;
}
/*============================================*/

void free_vdwcomb(t_vdwcomb *vdw)
{
  sfree(vdw->type1);
  sfree(vdw->type2);
  sfree(vdw->vdw);
  sfree(vdw);
}

/*============================================*/

t_vdwcomb *read_vdwcomb(int num, bool flag)
{
  FILE *fp = NULL;
  char sstring[STRLEN];
  if(!flag)
    strcpy(sstring,"[ COMBINATIONS ]"); 
  else strcpy(sstring,"[ 14_COMBINATIONS ]");
    
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



  t_vdwcomb *vdw = vdwcomb_init();
  vdw = vdwcomb_realloc(vdw,200); /* lazy */
  int count = 0;
  rewind(fp);
  char line[STRLEN];
  while(get_a_line(fp,line,STRLEN)) {
    /* remove the # */
    cut_string(line,'#');
    if(strlen(line) > 0){
      if(!strcmp(line,sstring)){
        while(get_a_line(fp,line,STRLEN) && strchr(line,'[')==NULL){
          cut_string(line,'#');
          if(strlen(line) > 0){
#ifdef GMX_DOUBLE
            if(sscanf(line,"%s %s %lf",vdw->type1[count],vdw->type2[count],
                      &vdw->vdw[count]) != 3)
#else
            if(sscanf(line,"%s %s %f",vdw->type1[count],vdw->type2[count],
                      &vdw->vdw[count]) != 3)
#endif
              fatal_error("Lib file corrupted\n");
            count++;
            if(count >= 200) {
              fprintf(stderr,"Eeek: that's a big one...\n");
              fprintf(stderr,"Increase array size in 'read_vdwcomb'");
              fatal_error("ran out of memory\n");
            }
          }
        }
      }
    }
  }
  fclose(fp);
  
  vdw->n = count;
  return vdw;
}
/*============================================*/


extern t_vdwcomb *read_vdwcomb2(FILE *fp, bool flag)
{
  rewind(fp);

  char sstring[STRLEN];
  if(flag==0 || flag == 2)
    strcpy(sstring,"[ COMBINATIONS ]");
  else strcpy(sstring,"[ 14_COMBINATIONS ]");


  t_vdwcomb *vdw = vdwcomb_init();
  vdw = vdwcomb_realloc(vdw,200); /* lazy */
  int count = 0;
/*   rewind(fp); */
  char line[STRLEN];
  while(get_a_line(fp,line,STRLEN)) {
    /* remove the # */
    cut_string(line,'#');
    if(strlen(line) > 0){
      if(!strcmp(line,sstring)){
        while(get_a_line(fp,line,STRLEN) && strchr(line,'[')==NULL){
/*           cut_string(line,'#'); */
/*           if(strlen(line) > 0){ */
/*           printf("line = %s\n",line); */

            if(flag!=2){
#ifdef GMX_DOUBLE
              if(sscanf(line,"%s %s %lf",vdw->type1[count],vdw->type2[count],
                        &vdw->vdw[count]) != 3)
#else
              if(sscanf(line,"%s %s %f",vdw->type1[count],vdw->type2[count],
                        &vdw->vdw[count]) != 3)
#endif
                fatal_error("Lib file corrupted\n");
            }
            else
            {
/*               printf("servus\n"); */
#ifdef GMX_DOUBLE
              if(sscanf(line,"%d %d %lf",&vdw->rs_type1[count],&vdw->rs_type2[count],
                        &vdw->rs_comb[count]) != 3)
#else
              if(sscanf(line,"%d %d %f",&vdw->rs_type1[count],&vdw->rs_type2[count],
                        &vdw->rs_comb[count]) != 3)
#endif
                fatal_error("Lib file corrupted\n");
/*               printf("here\n"); */
            }

            count++;
            if(count >= 200) {
              fprintf(stderr,"Eeek: that's a big one...\n");
              fprintf(stderr,"Increase array size in 'read_vdwcomb'");
              fatal_error("ran out of memory\n");
            }
/*           } */
        }
      }
    }
  }

  vdw->n = count;
  return vdw;
}
