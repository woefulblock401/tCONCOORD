#include <tconcoord.h>

/*==============================================================*/

void CNCwarn(FILE *log,char *warn)
  
{
  char string[STRLEN];
  strcpy(string,"tCNC__wrn_> ");
  strcat(string,warn);
  fprintf(log,string);
/*   fprintf(stderr,string); */
}

/*==============================================================*/

void CNClog(FILE *log,char *string)
  
{
  fprintf(log,string);
  fprintf(stderr,string); 
  fflush(log);
  fflush(stderr);
  
}
/*==============================================================*/
void CNCerr(FILE *log,char *err)
  
{
  char string[STRLEN];
  strcpy(string,"tCNC__err_> ");
  strcat(string,err);
  fprintf(log,string);
  fprintf(stderr,string);
  exit(1);
}
  
/*==============================================================*/

void fatal_error(char *error)
{
  CNCerr(stderr,error);
  
}

