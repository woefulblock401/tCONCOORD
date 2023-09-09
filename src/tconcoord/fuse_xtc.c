#include <tconcoord.h>

int main(int argc, char **argv)
{
  static char *desc[] = {
    "Concatenate xtc files",
    "...without losing frames..."
  };

 t_filenm fnm[] = {
      { efTRX, "-f",     NULL,      ffRDMULT },
      { efTRX, "-o",     "comb", ffWRITE }
  };

 t_pargs pa[] = {
 };
 
 #define NFILE asize(fnm)

  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);


  char        **fnms,*in_file;
  int nfile_in;
  nfile_in = opt2fns(&fnms,"-f",NFILE,fnm);
  int xtc;
  int i;
  int natoms;
  int step;
  real time;
  matrix box;
  rvec *x;
  real prec;
  bool bOk;
  int newstep = 0;
  real newtime = 0.;
  

  int fp = open_xtc(opt2fn("-o",NFILE,fnm),"w");
  for(i=0;i<nfile_in;i++){
    printf("%s\n",fnms[i]);
    xtc = open_xtc(fnms[i],"r");
    if (!read_first_xtc(xtc,&natoms,&step,&time,box,&x,
                        &prec,&bOk)) {
      fprintf(stderr,"Failed to read %s.....\n", fnms[i]);
/*       exit(1); */
    }
    do{
      newtime+=1.0;
      newstep++;
      printf("writing step = %d\n",step);
      write_xtc(fp,natoms,newstep,newtime,box,x,prec);
    }while(read_next_xtc(xtc,natoms,&step,&time,box,x,&prec,&bOk));
  }
  printf("Wrote %d xtc frames to %s\n",newstep,opt2fn("-o",NFILE,fnm));
  return 0;
}


/*   int xtc = open_xtc(opt2fn("-f",NFILE,fnm),"r"); */
