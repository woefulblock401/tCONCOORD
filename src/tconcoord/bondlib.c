#include <tconcoord.h>


/*======================================================*/

t_bondlib *bondlib_init(void)
{
  t_bondlib *bl = NULL;
  snew(bl,1);
  bl->n = 0;
  bl->type1 = bl->type2 = bl->type3 = NULL;
  bl->av = bl->lb = bl->ub = NULL;
  bl->ang = bl->sig = NULL;
  bl->is_bond = bl->is_ang = NULL;
  return bl;
}
/*======================================================*/

t_bondlib *bondlib_realloc(t_bondlib *bl, int n)
{
  srenew(bl->type1,n);
  srenew(bl->type2,n);
  srenew(bl->type3,n);
  srenew(bl->av,n);
  srenew(bl->lb,n);
  srenew(bl->ub,n);
  srenew(bl->ang,n);
  srenew(bl->sig,n);
  srenew(bl->is_bond,n);
  srenew(bl->is_ang,n);
  return bl;
}
/*======================================================*/

t_bondlib *read_bonds_and_angles(int bonds)
{
  char filename[STRLEN];
  FILE *fp = NULL;
  if(bonds == 1)
    fp=cnclib("Constraints.dat");
  else if(bonds == 2)
    fp=cnclib("Constraints_eh.dat");
  else{
    fprintf(stderr,"Fatal Error: Unknown parameter set\n");
    exit(0);
  }
  
  t_bondlib *bl = bondlib_init();
  char line[STRLEN];
  
  while(get_a_line(fp,line,STRLEN) && strchr(line,'[') != NULL)
  {
    if(strcmp(line,"[ BONDS ]") == 0){
      while(get_a_line(fp,line,STRLEN) && strchr(line,'[') == NULL){
        bl->n+=1;
        bl = bondlib_realloc(bl,bl->n);
#ifdef GMX_DOUBLE
        sscanf(line,"%s %s %lf %lf %lf",bl->type1[bl->n-1],bl->type2[bl->n-1],
               &bl->av[bl->n-1],&bl->lb[bl->n-1],&bl->ub[bl->n-1]);
#else
        sscanf(line,"%s %s %f %f %f",bl->type1[bl->n-1],bl->type2[bl->n-1],
               &bl->av[bl->n-1],&bl->lb[bl->n-1],&bl->ub[bl->n-1]);
#endif
        bl->is_bond[bl->n-1] = TRUE;
        bl->is_ang[bl->n-1] = FALSE;
      }
    }
    if(strcmp(line,"[ ANGLES ]") == 0){
      while(get_a_line(fp,line,STRLEN) && strchr(line,'[') == NULL){
        bl->n+=1;
        bl = bondlib_realloc(bl,bl->n);
#ifdef GMX_DOUBLE
        sscanf(line,"%s %s %s %lf %lf %lf %lf %lf",bl->type1[bl->n-1],bl->type2[bl->n-1],
               bl->type3[bl->n-1],&bl->av[bl->n-1],&bl->lb[bl->n-1],&bl->ub[bl->n-1],
               &bl->ang[bl->n-1],&bl->sig[bl->n-1]);
#else
        sscanf(line,"%s %s %s %f %f %f %f %f",bl->type1[bl->n-1],bl->type2[bl->n-1],
               bl->type3[bl->n-1],&bl->av[bl->n-1],&bl->lb[bl->n-1],&bl->ub[bl->n-1],
               &bl->ang[bl->n-1],&bl->sig[bl->n-1]);
#endif
        bl->is_ang[bl->n-1] = TRUE;
        bl->is_bond[bl->n-1] = FALSE;
      }
    }
  }
  return bl;
}


/*=============================================================*/
bool get_bond(char *type1, char *type2, t_bondlib *bl, real *av, real *lb, real *ub)

/* get lib entries for a bond */

{
  int i;
  for(i=0;i<bl->n;i++){
    if(bl->is_bond[i]){
      if((strcmp(type1,bl->type1[i])==0 && strcmp(type2,bl->type2[i]) == 0) ||
         (strcmp(type1,bl->type2[i])==0 && strcmp(type2,bl->type1[i]) == 0)){
/*         fprintf(stderr,"Reading from lib.... %s-%s\n",type1,type2); */
        
        *av = bl->av[i];
        *lb = bl->lb[i];
        *ub = bl->ub[i];
        return TRUE;
      }
    }
  }

  return FALSE;
}
/*=============================================================*/
bool get_angle(char *type1, char *type2, char *type3, t_bondlib *bl,
               real *av, real *lb, real *ub, real *avang, real *sig)

/* get lib entries for an angle */ 

{
  int i;
  for(i=0;i<bl->n;i++){
    if(bl->is_ang[i]){
      if((strcmp(type1,bl->type1[i])==0 && strcmp(type2,bl->type2[i]) == 0 && strcmp(type3,bl->type3[i]) == 0) ||
         (strcmp(type1,bl->type3[i])==0 && strcmp(type2,bl->type2[i]) == 0 && strcmp(type3,bl->type1[i]) == 0) ){
        *av = bl->av[i];
        *lb = bl->lb[i];
        *ub = bl->ub[i];
        *avang = bl->ang[i];
        *sig = bl->sig[i];
        return TRUE;
      }
    }
  }

  return FALSE;
}
/*=============================================================*/

