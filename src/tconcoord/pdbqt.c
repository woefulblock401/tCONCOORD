#include <tconcoord.h>

/*============================================================*/

t_atomlist *read_pdbqt(char *filename, FILE *log, bool bVerbose)
{
  FILE *fp = ffopen(filename,"r");
  char line[STRLEN];
  char error[STRLEN];
  int count = 0;
  t_atomlist *al = atomlist_init();
  int  i,j,k;
  char nc='\0';
  char anr[12],anm[12],altloc,resnm[12],chain[12],resnr[12];
  char xc[12],yc[12],zc[12],occup[12],bfac[12],pdbresnr[12];
  char charge[12], atype[12];
  
  if(log==NULL) log = stderr;
  
  /* count atoms */
  while(get_a_line(fp,line,STRLEN)){
    if(strncmp(line,"ATOM",4) ==0 ||
       strncmp(line,"HETATM",6) == 0){
      count++;
    }
  }
  if(bVerbose)
    fprintf(stderr,"Read %d atoms from %s\n",count,filename);
    
  al = al_realloc(al,count);
  rewind(fp);
  i = 0;
  while(get_a_line(fp,line,STRLEN)){
    if(strncmp(line,"ATOM",4) ==0 ||
       strncmp(line,"HETATM",6) == 0){
      
      /* Skip over type */  
      j=6;
      for(k=0; (k<5); k++,j++) anr[k]=line[j];
      anr[k]=nc;
      trim(anr);
      j++;
      for(k=0; (k<4); k++,j++) anm[k]=line[j];
      anm[k]=nc;
      trim(anm);
      altloc=line[j];
      j++;
      for(k=0; (k<4); k++,j++) 
        resnm[k]=line[j];
      resnm[k]=nc;
      trim(resnm);
      
      for(k=0; (k<1); k++,j++)
        chain[k]=line[j];
      chain[k]=nc;
      
      for(k=0; (k<4); k++,j++) {
        resnr[k]=line[j];
        pdbresnr[k]=line[j];
      }
      resnr[k]=nc;
      trim(resnr);
      pdbresnr[k]=line[j];
      pdbresnr[k+1]=nc;
      trim(pdbresnr);
      j+=4;
      
      /* X,Y,Z Coordinate */
      for(k=0; (k<8); k++,j++) xc[k]=line[j];
      xc[k]=nc;
      for(k=0; (k<8); k++,j++) yc[k]=line[j];
      yc[k]=nc;
      for(k=0; (k<8); k++,j++) zc[k]=line[j];
      zc[k]=nc;
      
      /* Weight */
      for(k=0; (k<6); k++,j++) occup[k]=line[j];
      occup[k]=nc;
      /* B-Factor */
      for(k=0; (k<7); k++,j++) bfac[k]=line[j];
      bfac[k]=nc;
      
      for(k=0; (k<9); k++,j++) charge[k]=line[j];
      charge[k]=nc;

      for(k=0; (k<4); k++,j++) atype[k]=line[j];
      atype[k]=nc;
      trim(atype);
      
      al->id[i] = atoi(anr);
      strcpy(al->name[i],anm);
      extend_name(al->name[i]);
      strcpy(al->resname[i],resnm);
      strcpy(al->chain[i],chain);
      al->resid[i] = atoi(pdbresnr);
      al->x[i][XX] = atof(xc);
      al->x[i][YY] = atof(yc);
      al->x[i][ZZ] = atof(zc);
      al->occ[i] = atof(occup);
      al->bfac[i] = atof(bfac);
      al->q[i] = atof(charge);
      strcpy(al->type[i],atype);
      i++;
    }
  }
  return al;
}


      
