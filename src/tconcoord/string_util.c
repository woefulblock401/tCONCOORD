
#include <tconcoord.h>
/*============================================================*/

void cut_string(char *in, char c)
{
  int i;
  for(i=0;i<strlen(in);i++){
    if(in[i]==c){
      in[i]=0;
      break;
    }
  }
}
/*============================================================*/

void slice_string(char *dest, char *src, int beg, int end){
  int i;
  for(i=beg;i<end;i++){
    dest[i-beg]=src[i];
  }
  dest[end-beg]=0;
}

/*============================================================*/

void switch_string(char *s1,char *s2)
{
  longstr dummy;
  strcpy(dummy,s1);
  strcpy(s1,s2);
  strcpy(s2,dummy);
}
/*============================================================*/

bool find_key_word(FILE *fp,char *key)
{
  char line[STRLEN];
  int size = strlen(key);
  while(get_a_line(fp,line,STRLEN))
  {
    if(strncmp(key,line,size)==0)
      return TRUE;
  }
  return FALSE;
}

/*============================================================*/

void strip(char *out,char *in)
{
    sscanf(in,"%s",out);
}

/*============================================================*/

void substring(char *out, char *in, char begin, char end)
{
  int i,k;
  for(i=0;i<strlen(in);i++){
    if(in[i]==begin){
      if(end=='|'){
        for(k=i;k<strlen(in);k++){
          out[k-i]=in[k];
        }
        out[k-i]=0;
      }
      else{
        for(k=i;k<strlen(in);k++){
          out[k-i]=in[k];
          if(in[k]==end) {
            out[k-i+1]=0;
            break;
          }
        }
      }
    }
  }
}
/*============================================================*/

void insert_string(char *string, char *insert, int pos, char *flag)
{
  int i,k,len,newlen;
  len=strlen(string)-1;
  newlen=strlen(string)+strlen(insert)-1;
  if(!strcmp(flag,">")){
    i=len;
    while(i>=pos){
      string[i+strlen(insert)]=string[i]; 
      i--;
    }
    for(k=0;k<strlen(insert);k++){
      string[pos+k]=insert[k];
    }
    string[newlen+1]=0;
  }
  else if(!strcmp(flag,"|")){
    for(k=0;k<strlen(insert);k++){
      string[pos+k]=insert[k];
    }
    string[len+1]=0;
  }
}

/*============================================================*/


void extend_name(char *name)
{
  if(strlen(name) == 1){
    insert_string(name," ",0,">");
    insert_string(name,"  ",2,">");
  }
  else if(strlen(name) == 2){
    if(!isdigit(name[0])){
      insert_string(name," ",0,">");
      insert_string(name," ",3,">");
    }
    else{
      insert_string(name,"  ",2,">");
    }
  }
  else if(strlen(name) == 3){
    if(!isdigit(name[0])){
      insert_string(name," ",0,">");
    }
    else{
      insert_string(name," ",3,">");
    }
  }
  
}
/*============================================================*/

void progress(FILE *fp,char *string, int now, int full)
{
  real stat = now/(real) full;
  fprintf(fp,"\r%-30s ===> %3.0f %%",string,stat*100);
  fflush(fp);
}
/*============================================================*/


void cnc_copyright(char *prgname)
{
  fprintf(stderr,"\n=======================================================\n");
  fprintf(stderr,"\n\t tCONCOORD Version 1.0\n");
  fprintf(stderr,"\t Program : %s\n",prgname);
  fprintf(stderr,"\t (c) Daniel Seeliger ( dseelig@gwdg.de )\n");
  fprintf(stderr,"\t Please cite:\n");
  fprintf(stderr,"\t Daniel Seeliger, Juergen Haas and Bert L. de Groot\n");
  fprintf(stderr,"\t \"Geometry-based Sampling of Conformational Transitions\n");
  fprintf(stderr,"\t in Proteins\", Structure 15, 1482-1492 (2007)\n");
  fprintf(stderr,"\n=======================================================\n");
}

/*============================================================*/

void papers_log(FILE *log)
{
  char logstr[STRLEN];
  sprintf(logstr,"tCNC__log_>--------------------------------------------------------\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> ++++ Original tCONCOORD publication: ++++\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> [1] Daniel Seeliger, Juergen Haas and Bert L. de Groot\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   \"Geometry-based Sampling of Conformational Transitions\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   in Proteins\", Structure 15, 1482-1492 (2007)\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> ++++ Parameter set used for this simulation: ++++\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> [2] Daniel Seeliger and Bert L. de Groot\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   \"Atomic contacts in protein structure. A detailed analysis\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   of atomic radii, packing and overlaps\", Proteins 68, 595-601 (2007)\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> ++++ Graphical frontend for PyMOL: ++++\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_> [3] Daniel Seeliger and Bert L. de Groot\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   \"tCONCOORD-GUI. Visually supported conformational sampling\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>   of bioactive molecules\", J. Comp. Chem. in press (2008)\n");
  CNClog(log,logstr);
  sprintf(logstr,"tCNC__log_>--------------------------------------------------------\n");
  CNClog(log,logstr);
};



int str_nelem(char *str,int maxptr,char *ptr[])
{
  int  np=0;
  char *copy0,*copy;
  
  copy0=strdup(str); 
  copy=copy0;
  ltrim(copy);
  while (*copy != '\0') {
    if (np >= maxptr)
      gmx_fatal(FARGS,"Too many groups on line: '%s' (max is %d)",
                  str,maxptr);
    if (ptr) 
      ptr[np]=copy;
    np++;
    while ((*copy != '\0') && !isspace(*copy))
      copy++;
    if (*copy != '\0') {
      *copy='\0';
      copy++;
    }
    ltrim(copy);
  }
  if (ptr == NULL)
    sfree(copy0);

  return np;
}

