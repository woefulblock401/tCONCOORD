#include <tconcoord.h>


t_bounds *read_pharmacophore(char *filename)
{
  FILE *fp = ffopen(filename,"r");
  char line[STRLEN];
  int i = 0;
  int dum;
  int nbounds;
  nbounds = count_bounds(fp);
  t_bounds *b = bounds_init(nbounds);
  rewind(fp);
  bool acc_acc = FALSE;
  bool don_acc = FALSE;
  bool don_don = FALSE;
  

  while(get_a_line(fp,line,STRLEN)){
    if(strchr(line,'[')!=NULL){
      if(strcmp(line,"[ ACC_ACC ]") == 0){
        acc_acc = TRUE;
        don_acc = FALSE;
        don_don = FALSE;
      }
      
      else if(strcmp(line,"[ DON_ACC ]") == 0){
        acc_acc = FALSE;
        don_acc = TRUE;
        don_don = FALSE;
      }
        
      else if(strcmp(line,"[ DON_DON ]") == 0){
        acc_acc = FALSE;
        don_acc = FALSE;
        don_don = TRUE;
      }
    }
    else {
      i++;
#ifdef GMX_DOUBLE
      if(sscanf(line,"%d %d %lf %lf %lf",&b->at1[i-1],&b->at2[i-1],
                &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#else
        if(sscanf(line,"%d %d %f %f %f",&b->at1[i-1],&b->at2[i-1],
                  &b->av[i-1],&b->lb[i-1],&b->ub[i-1]) != 5)
#endif
          input_error(i,line);
        b->ac_ac[i-1] = acc_acc;
        b->don_ac[i-1] = don_acc;
        b->don_don[i-1] = don_don;
        b->n++;
    }
  }
  return b;
}

bool check_acc_acc(t_atomlist *al, t_idxgroups *acc,t_bounds *b)
{
  
  int i,k,l;
  int count = 0;
  int at1, at2;
  real d;

  for(l=0;l<b->n;l++){
    if(b->ac_ac[l]) {
      b->bCheck[l] = TRUE;
    }
  }
  
  for(i=0;i<acc->natoms[0]-1;i++){
    for(k=i+1;k<acc->natoms[0];k++){
      at1 = acc->atoms[0][i];
      at2 = acc->atoms[0][k];
      d = DIST(al,at1,at2);
      for(l=0;l<b->n;l++){
        if(b->ac_ac[l]) {
          if(d >= b->lb[l] && d <= b->ub[l]) {
            count++;
            b->at3[l] = at1;
            b->at4[l] = at2;
            b->bCheck[l] = FALSE;
          }
        }
      }
    }
  }

  for(l=0;l<b->n;l++){
    if(b->ac_ac[l] && b->bCheck[l]) {
      return FALSE;
    }
  }
  return TRUE;
}

bool check_don_don(t_atomlist *al, t_idxgroups *don,t_bounds *b)
{
  
  int i,k,l;
  int count = 0;
  int at1, at2;
  real d;
  for(l=0;l<b->n;l++){
    if(b->don_don[l]) {
      b->bCheck[l] = TRUE;
    }
  }
  
  for(i=0;i<don->natoms[0]-1;i++){
    for(k=i+1;k<don->natoms[0];k++){
      at1 = don->atoms[0][i];
      at2 = don->atoms[0][k];
      d = DIST(al,at1,at2);
      for(l=0;l<b->n;l++){
        if(b->don_don[l]) {
          if(d >= b->lb[l] && d <= b->ub[l]) {
            count++;
            b->at3[l] = at1;
            b->at4[l] = at2;
            b->bCheck[l] = FALSE;
          }
        }
      }
    }
  }
  for(l=0;l<b->n;l++){
    if(b->don_don[l] && b->bCheck[l]) {
      return FALSE;
    }
  }
  return TRUE;
}

bool check_don_acc(t_atomlist *al, t_idxgroups *acc, t_idxgroups *don, t_bounds *b)
{
  
  int i,k,l;
  int count = 0;
  int at1, at2;
  real d;
  for(l=0;l<b->n;l++){
    if(b->don_ac[l]) {
      b->bCheck[l] = TRUE;
    }
  }
  
  for(i=0;i<don->natoms[0];i++){
    for(k=0;k<acc->natoms[0];k++){
      at1 = don->atoms[0][i];
      at2 = acc->atoms[0][k];
      d = DIST(al,at1,at2);
      for(l=0;l<b->n;l++){
        if(b->don_ac[l]) {
          if(d >= b->lb[l] && d <= b->ub[l]) {
            b->at3[l] = at1;
            b->at4[l] = at2;
            b->bCheck[l] = FALSE;
            count++;
          }
        }
      }
    }
  }
  for(l=0;l<b->n;l++){
    if(b->don_ac[l] && b->bCheck[l]) {
      return FALSE;
    }
  }
  return TRUE;
}

int count_acc_acc(t_bounds *b)
{
  int i;
  int count = 0;
  
  for(i=0;i<b->n;i++){
    if(b->ac_ac[i]) count++;
  }
  return count;
}

int count_don_acc(t_bounds *b)
{
  int i;
  int count = 0;
  
  for(i=0;i<b->n;i++){
    if(b->don_ac[i]) count++;
  }
  return count;
}

int count_don_don(t_bounds *b)
{
  int i;
  int count = 0;
  
  for(i=0;i<b->n;i++){
    if(b->don_don[i]) count++;
  }
  return count;
}


int get_atnum(t_bounds *b)
{
  int i;
  int tmp = 0;
  for(i=0;i<b->n;i++){
    if(b->at2[i] > tmp) tmp = b->at2[i];
  }
  return tmp;
}

int pick_atom(t_idxgroups *idx,int start,int not)
{
  int i;
  for(i=start;i<idx->natoms[0];i++){
    if(idx->atoms[0][i]!=not) return idx->atoms[0][i];
  }
  return -1;
}




int main(int argc,char *argv[])
{

  FILE *log = NULL;

  t_atomlist *al = atomlist_init();
  t_idxgroups *don = idx_init();
  t_idxgroups *acc = idx_init();
  t_idxgroups *phob = idx_init();
  t_idxgroups *pl = idx_init();
  t_idxgroups *imp = idx_init();
  t_trxframe fr;
  int        status;
  int        flags = TRX_READ_X;
  int n_ac_ac = 0;
  int n_don_ac = 0;
  int n_don_don = 0;
  real d;
  


  static char *desc[] = {
    ""
  };

  t_pargs pa[] = {
/*      { "-v",   FALSE, etBOOL, {&bVerbose},  */
/*        "Make noise" }, */
/*      { "-n",   FALSE, etINT, {&nstructs},  */
/*        "Number of structures to be generated" }, */
/*      { "-i",   FALSE, etINT, {&maxit} */
     };
     
     
  t_filenm fnm[] = {
    { efCTP,"-top","topol", ffREAD},
    { efTRX,"-f","traj", ffREAD},
    { efDAT,"-p","pharm", ffREAD},
    { efPDB,"-o","hits", ffWRITE},
    { efLOG,"-log","pharm", ffWRITE}
  };


#define NFILE asize(fnm)
     
     
  parse_common_args(&argc,argv,0,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

     
       

  log = ffopen(opt2fn("-log",NFILE,fnm),"w"); 
  t_bounds *b = read_pharmacophore(opt2fn("-p",NFILE,fnm));
  fprintf(stderr,"Reading Topology.....\n");
  read_cnctop(opt2fn("-top",NFILE,fnm),al,acc,don,phob,pl,imp);
  get_bconstr(al);  
  occ_to_one(al);
  int i,k,l;
  
  int nn = get_atnum(b);
  printf("using %d point pharmacophore model\n",nn);
  int points[nn];
  int atoms[nn];
  
  for(i=0;i<nn;i++){
    points[i] = -1;
    atoms[i] = -1;
  }
  
  /* check whether points are donors or acceptors */

  for(i=0;i<b->n;i++){
    if(b->ac_ac[i]) {
      if(points[b->at1[i]-1]==-1) points[b->at1[i]-1] = 0;
      if(points[b->at2[i]-1]==-1) points[b->at2[i]-1] = 0;
    }
    else if(b->don_ac[i]) {
      if(points[b->at1[i]-1]==-1) points[b->at1[i]-1] = 0;
      if(points[b->at2[i]-1]==-1) points[b->at2[i]-1] = 1;
    }
    else if(b->don_don[i]) {
      if(points[b->at1[i]-1]==-1) points[b->at1[i]-1] = 1;
      if(points[b->at2[i]-1]==-1) points[b->at2[i]-1] = 1;
    }
  }
  
  for(i=0;i<nn;i++){
    printf("point %d = %d\n",i+1,points[i]);
  }
  
  /* fill atoms with points */
  


  /* find all permutations and store in group */

  t_idxgroups *perms = idx_init();
  

  int p[nn];
  for(i=0;i<nn;i++) p[i] = 0;
  
  if(points[0] == 0){
    /* pick first atoms */
    for(i=0;i<acc->natoms[0];i++){
      atoms[0] = acc->atoms[0][i];
      p[0]++;
      if(points[1] == 0){
        atoms[1] = 0;
        p[1]=0;
        
        do
        {
          atoms[1] = pick_atom(acc,p[1],atoms[0]);
          p[1]++;
          if(atoms[1]!=-1){
            if(points[2] == 0){
              atoms[2] = 0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(acc,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
                
              }while(atoms[2]!=-1);
            }
            else if(points[2]==1){
              atoms[2]=0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(don,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
                
              }while(atoms[2]!=-1);
            }
          }
        }while(atoms[1]!=-1);
      }
      
      else if(points[1] == 1){
        atoms[1] = 0;
        p[1]=0;
        
        do
        {
          atoms[1] = pick_atom(don,p[1],atoms[0]);
          p[1]++;
          if(atoms[1]!=-1){
            if(points[2] == 0){
              atoms[2] = 0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(acc,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
              }while(atoms[2]!=-1);
            }
            else if(points[2]==1){
              atoms[2]=0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(don,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
                
              }while(atoms[2]!=-1);
            }
          }
        }while(atoms[1]!=-1);
      }
    }
  }

  if(points[0] == 1){
    /* pick first atoms */
    for(i=0;i<don->natoms[0];i++){
      atoms[0] = don->atoms[0][i];
      p[0]++;
      if(points[1] == 0){
        atoms[1] = 0;
        p[1]=0;
        
        do
        {
          atoms[1] = pick_atom(acc,p[1],atoms[0]);
          p[1]++;
          if(atoms[1]!=-1){
            if(points[2] == 0){
              atoms[2] = 0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(acc,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
              }while(atoms[2]!=-1);
            }
            else if(points[2]==1){
              atoms[2]=0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(don,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
                
              }while(atoms[2]!=-1);
            }
          }
        }while(atoms[1]!=-1);
      }
      
      else if(points[1] == 1){
        atoms[1] = 0;
        p[1]=0;
        
        do
        {
          atoms[1] = pick_atom(don,p[1],atoms[0]);
          p[1]++;
          if(atoms[1]!=-1){
            if(points[2] == 0){
              atoms[2] = 0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(acc,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
              }while(atoms[2]!=-1);
            }
            else if(points[2]==1){
              atoms[2]=0;
              p[2] = 0;
              do
              {
                atoms[2] = pick_atom(don,p[2],atoms[1]);
                p[2]++;
                if(atoms[2]!=-1 && atoms[2]!=atoms[0]){
                  perms->n+=1;
                  perms = idx_realloc(perms,perms->n);
                  add_to_group_nosort(perms,perms->n-1,atoms[0]);
                  add_to_group_nosort(perms,perms->n-1,atoms[1]);
                  add_to_group_nosort(perms,perms->n-1,atoms[2]);
                }
                
              }while(atoms[2]!=-1);
            }
          }
        }while(atoms[1]!=-1);
      }
    }
  }


  int nsucc = 0;
  FILE *fp = ffopen(opt2fn("-o",NFILE,fnm),"w");   
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags); 
  
  do { 
    rvec2al(al,fr.x); 
    /* now rebuild bounds */
    
    for(k=0;k<perms->n;k++){
      for(i=0;i<b->n;i++){
        if(b->at1[i]-1 == 0) b->at3[i] = perms->atoms[k][0];
        else if(b->at1[i]-1 == 1) b->at3[i] = perms->atoms[k][1];
        else if(b->at1[i]-1 == 2) b->at3[i] = perms->atoms[k][2];
        if(b->at2[i]-1 == 0) b->at4[i] = perms->atoms[k][0];
        else if(b->at2[i]-1 == 1) b->at4[i] = perms->atoms[k][1];
        else if(b->at2[i]-1 == 2) b->at4[i] = perms->atoms[k][2];
        b->bCheck[i] = TRUE;
      }
      for(i=0;i<b->n;i++){
        d = DIST(al,b->at3[i],b->at4[i]);
/*         printf("Checking constraint %s-%s %g %g %g ",al->name[b->at3[i]],  */
/*                al->name[b->at4[i]],d,b->lb[i],b->ub[i]); */
        
        
        if(d >=b->lb[i] && d <=b->ub[i]) {
/*           printf("Ok\n"); */
          b->bCheck[i] = FALSE;
        }
/*         else printf("\n"); */
        
      }
/*       printf("\n"); */
      
      bool bAll = TRUE;
      for(i=0;i<b->n;i++){
        if(b->bCheck[i] == TRUE) bAll = FALSE;
      }
      if(bAll) {
        nsucc++;
        fprintf(log,"======================================\n"); 
        fprintf(log,"MODEL %d\n",nsucc);
        fprintf(log,"ATOM %d %s\n",al->id[perms->atoms[k][0]],al->name[perms->atoms[k][0]]);
        fprintf(log,"ATOM %d %s\n",al->id[perms->atoms[k][1]],al->name[perms->atoms[k][1]]);
        fprintf(log,"ATOM %d %s\n",al->id[perms->atoms[k][2]],al->name[perms->atoms[k][2]]);
        write_pdb_frame(fp,al,nsucc); 
      }
    }
  } while(read_next_frame(status,&fr));   
  
  


/*   FILE *fp = ffopen(opt2fn("-o",NFILE,fnm),"w"); */
/*   int nsucc = 0; */
/*   real d; */
/*   int i; */

  

  
  
/*   do { */
/*     int nad = 0; */
/*     int naa = 0; */
/*     int ndd = 0; */
    
/*     rvec2al(al,fr.x); */
/*     naa = check_acc_acc(al,acc,b); */
/*     ndd = check_don_don(al,don,b);  */
/*     nad = check_don_acc(al,acc,don,b); */
    
  /*   if(naa >= n_ac_ac && */
/*        nad >= n_don_ac && */
/*        ndd >= n_don_don) */
/*     { */
/*     if(naa && nad && ndd) */
/*     { */
      
/*       nsucc++; */
/*       printf("Fulfilled constraints\n"); */
/*       fprintf(log,"======================================\n"); */
/*       fprintf(log,"MODEL %d\n",nsucc); */
/*       for(i=0;i<b->n;i++){ */
/*         d = DIST(al,b->at1[i],b->at2[i]); */
/*         fprintf(log,"%6d (%s) %6d (%s) %8.3f %8.3f %8.3f -> %8.3f\n", */
/*                 al->id[b->at1[i]],al->name[b->at1[i]], */
/*                 al->id[b->at2[i]],al->name[b->at2[i]], */
/*                 b->lb[i],b->av[i],b->ub[i],d); */
/*       } */
/*       fprintf(log,"======================================\n"); */
/*       write_pdb_frame(fp,al,nsucc); */
/*     } */
    

/*     printf("Coordinates at t=%8.3f : %8.5f %8.5f %8.5f\n",fr.time,fr.x[n][XX],fr.x[n][YY],fr.x[n][ZZ]); */
/*   } while(read_next_frame(status,&fr)); */


  
  return 0;
}



