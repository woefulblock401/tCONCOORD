#include <tconcoord.h>

enum {
  FTOT,FB,FANG,FBO,FNB,FP,FCHI,FNR
  };


/*=======================================================*/
volatile int GotTermSignal = 0;

static RETSIGTYPE signal_handler(int n)
{
  fprintf(stderr,"Press Ctrl+C again to quit\n");
  GotTermSignal+=1;
/*   fprintf(stderr,"GotTermSig = %d\n",GotTermSignal); */
}

/*=======================================================*/

void clear_tags(int *tags, int n)
{
  int i;
  for(i=0;i<n;i++){
    tags[i] = FALSE;
  }
}

/*=======================================================*/

void update_tags(int *src, int *dest, int n)
{
  int i;
  for(i=0;i<n;i++){
    dest[i] = src[i];
  }
}

/*=======================================================*/
  
void shift_system(t_atomlist *al, int *notouch)
{
  int i,k;
  int size = 0;
  rvec co;
  clear_rvec(co);
  atom_id *index = NULL;

  for(i=0;i<al->natoms;i++){
    if(notouch[i]){
      size++;
      srenew(index,size);
      index[size-1] = i;
    }
  }
  rvec dum[size];
  com_frag(al,size,index,dum,co);
  
  for(i=0;i<al->natoms;i++){
    if(!notouch[i]){
      for(k=0;k<DIM;k++){
        al->x[i][k]-=co[k];
      }
    }
  }
  sfree(index);
}


/*=======================================================*/
void do_refine(t_atomlist *al, t_bounds *b,real bump_tol, real nbfac,real bfac, 
               t_idxgroups *pl, t_idxgroups *imp,t_boundtrack *bt, int *tags,
               int *maxid, real *maxf, real *F, bool bDo, bool bVerbose,
               gmx_rng_t rng, real angfac)
{
  
  clear_forces(al); 
  real bforce = 0;
  real bforce2 = 0;
  real nbforce = 0;
  real plforce = 0;
  
  *maxf = 0;
  *maxid = -1;
  int i,k;
  
  for(i=0;i<FNR;i++) F[i] = 0.;
  
  F[FB] = do_bound_force(al, b, 1., 1.,1.,rng);  
  F[FANG] = do_angle_force(al,b,angfac);  
  F[FBO] = do_bound_force2(al, b, .5);  
  F[FNB] = do_nb_force(al,bt,bump_tol,1.);      
/*   F[FCHI] = do_chiral_force(al,imp,1.);       */
  
  for(k=0;k<pl->n;k++){
    F[FP] += do_planar_force(al,pl->natoms[k],pl->atoms[k],pl->val[k],1.);
  }
  
  F[FTOT] = sum_forces(al);  
  
  *maxf = 0.;
  for(i=0;i<al->natoms;i++){
    real ff = norm(al->f[i])*bfac;
    
    if (ff>*maxf) {
      *maxf = ff;
      *maxid = i;
    }
  }

  
  if(bDo) {
    apply_f(al,bfac,tags); 
    plforce = 0; 
  
    clear_forces(al); 
    for(k=0;k<pl->n;k++){ 
      plforce += do_planar_force(al,pl->natoms[k],pl->atoms[k],pl->val[k],1); 
    } 
  
    real bviol = do_bound_force(al, b, 1., 1.,1.,rng); 
    bviol = do_angle_force(al,b,angfac);
/*     bviol = do_chiral_force(al,imp,1.);        */

    apply_f(al,bfac,tags); 
    clear_forces(al); 
  }
}



/*=======================================================*/
/* void stretch_coords(t_atomlist *al, int axis, real scale) */
/* { */
/*   int i; */
/*   for(i=0;i<al->natoms;i++){ */
/*     if(!al->isposres[i]) */
/*       al->x[i][axis]*=scale; */
/*   } */
/* } */
void stretch_coords(t_atomlist *al, real scale)
{
  int i;
  rvec cm;
  rvec vec;
  real d;
  
  calc_cent(al, cm);
  
  for(i=0;i<al->natoms;i++){
    if(!al->isposres[i]) {
      rvec_sub(cm,al->x[i],vec);
/*       d = absVec(vec); */
/*       scaleVector(vec, vec_scale); */
/*       printf("Moving %g A\n",scale*d); */
      
      al->x[i][XX] = cm[XX] - scale*vec[XX];
      al->x[i][YY] = cm[YY] - scale*vec[YY];
      al->x[i][ZZ] = cm[ZZ] - scale*vec[ZZ];
    }
  }
}

/*=======================================================*/
bool check_gyr(t_atomlist *al, real gyr, real gyr_tol, bool bVerbose)
{
  int m;
  real Rg;
  rvec gvec;
  bool bCorr = TRUE;
  
  Rg = radius_of_gyration(al,gvec);
  
/*   printf("before =  %g -> %g (%g)\n",Rg, gyr,gyr_tol);  */

  if (Rg > gyr + gyr_tol){
    bCorr = FALSE;
    if(bVerbose)
      printf("Correcting Radius of Gyration....\n"); 
    stretch_coords(al, 0.99);
/*     for(m=0;m<DIM;m++){ */
/*       stretch_coords(al,m,0.99); */
/*     } */
  }
  else if (Rg < gyr - gyr_tol){
    bCorr = FALSE;
    if(bVerbose)
      printf("Correcting Radius of Gyration....\n");
    stretch_coords(al, 1.01);
/*     for(m=0;m<DIM;m++){ */
/*       stretch_coords(al,m,1.01); */
/*     } */
  }
  
/*   Rg = radius_of_gyration(al,gvec);   */
/*   printf("after =  %g -> %g (%g)\n",Rg, gyr,gyr_tol);   */

  return bCorr;
}


/*=======================================================*/

bool cross_check_dist(FILE *log,t_atomlist *al, t_bounds *b, real *dviol, 
                      real *aviol,gmx_rng_t rng, t_idxgroups *pl, real btol, 
                      t_boundtrack *bt, bool bVerbose,
                      int nstep, int nrefst, real gyr, real gyr_tol)
{

  int i,k;
  real d,ang,lb,ub;
  real viol;
  int nang = 0;
  int ang_viol = 0;
  int ndist = 0;
  int dist_viol = 0;
  real bond_viol = 0;
  char logstr[STRLEN];
  
  int ndih = 0;
  int dih_viol = 0;
  int noth = 0;
  
  int oth_viol = 0;
  int ntot = 0;
  *dviol = *aviol = 0;
  real worst = 0.; 
  real dworst = 0.;
  real aworst = 0.;
  real range;
  int mem = 0;
  int amem = 0;
  real dmem = 0;
  real bsum;
  real plsum;
  real wang = 0;
  int pln;
  real maxpl = 5.;
  real max_bond = 0.1;
  real max_ang_d = 0.2;
  real max_ang = 20.;
  real max_dih = 0.1;
  bool bondOk = TRUE;
  bool angOk = TRUE; 
  bool dihOk = TRUE;
  bool gyrOk = TRUE;
  rvec gvec;
  real Rg;
  bool this_angle_ok = TRUE;
  
  real tol = 0.2;
  
  real pworst;
  

  /* check bounds */

  for(i=0;i<b->n;i++){
    this_angle_ok = TRUE;
    if(b->isbond[i]) ndist++;
    else if(b->isang[i]) nang++;
    else if(b->isdih[i]) ndih++;
    else noth++;
    d = DIST(al,b->at1[i],b->at2[i]);
    range = b->ub[i] - b->lb[i];
    if(d < b->lb[i]){
      viol = b->lb[i]-d;
      *dviol+=viol;
      ntot++;
      if(b->isbond[i]) {
        dist_viol++;
        bond_viol+=viol;
        if(viol > max_bond) bondOk = FALSE;
      }
      else if(b->isang[i]) {
        if(viol > max_ang_d) {
          angOk = FALSE;
          ang_viol++;
          this_angle_ok = FALSE;
        }
      }
      
      else if(b->isdih[i]) {
        dih_viol++;
        if(viol > max_dih) dihOk = FALSE;
      }
      else oth_viol++;
      

      if(viol > dworst) dworst = viol;
      if(viol/range > worst) {
        worst = viol/range;
        mem = i;
        dmem = d;
      }
        
    }
    else if(d > b->ub[i]){
      viol = (d-b->ub[i]);
      if(b->isbond[i]) {
        dist_viol++;
        bond_viol+=viol;
        if(viol > max_bond) bondOk = FALSE;
      }
      else if(b->isang[i]) {
        if(viol > max_ang_d) {
          angOk = FALSE;
          this_angle_ok = FALSE;
          ang_viol++;
        }
      }
      else if(b->isdih[i]) {
        dih_viol++;
        if(viol > max_dih) dihOk = FALSE;
      }
      else oth_viol++;
      
      if(viol > dworst) dworst = viol;
      *dviol+=viol; 
      ntot++;
      if(viol/range > worst) {
        worst = viol/range;
        mem = i;
        dmem = d;
      }
    }

    /* check angles */

    if(b->isang[i] && this_angle_ok){

      ang = RAD2DEG*angle_ij_ik(al->x[b->at3[i]],
                                al->x[b->at1[i]],
                                al->x[b->at2[i]]);
      range = 2*b->sig[i];

      lb = b->ang[i] - b->sig[i];
      ub = b->ang[i] + b->sig[i];
      
      if(ang  < lb)
      {
        viol = lb - ang;
        if (viol > max_ang) angOk = FALSE;
        
        *aviol+= viol;
        ang_viol++;
        ntot++;
        if(viol > aworst) {
          aworst = viol;
          amem = i;
          wang = ang;
          
        }
      }
      else if(ang > ub)
      {
        viol = ang - ub;
        if (viol > max_ang) angOk = FALSE;
        *aviol+=viol;
        ang_viol++;
        ntot++;
        if(viol> aworst) {
          aworst = viol;
          amem = i;
          wang = ang;
          
        }
      }
    }
  }

  
  plsum = 0;
  bool plOk;
  
  plOk = count_planar(al,pl,&plsum,&pln,maxpl,&pworst);
/*   if(bVerbose && bVerbose!=22){  */
/*     fprintf(stderr,"Checking bounds at step = %d\n",nstep);  */
/*     fprintf(stderr,"Number of planarity violations = %d\n",pln);  */
/*     fprintf(stderr,"Maximum violation = %10.3f%%\n",pworst);  */
/*     fprintf(stderr,"Sum of violations = %10.3f A\n",plsum);  */
/*     fprintf(stderr,"Planarity constraints fulfilled = (%d of %d)\n",  */
/*             pl->n - pln, pl->n);  */
/*   }  */
  
  bsum = 0;
  int bump = count_bumps(al,bt,btol,&bsum,0.5,bVerbose);
  
  char bflag[STRLEN];
  strcpy(bflag,"nobond");
  if(b->isbond[mem]) strcpy(bflag,"bond");
  if(b->isang[mem]) strcpy(bflag,"angle");
  if(b->isdih[mem]) strcpy(bflag,"dihedral");
 

  bool ret = TRUE;
  char dum[STRLEN];
  
  if(bondOk) {
    strcpy(dum,"Ok");
  }
  else {
    strcpy(dum,"Not Ok");
    ret = FALSE;
  }
  if(bVerbose && bVerbose!=22)
    fprintf(stderr,"Bonds.................: %s (%5.2f%%) (%6.5f A per bond) Worst violation: %6.5f\n",
            dum,100.-dist_viol/(real) ndist*100, bond_viol/(real) ndist, dworst);
  if(angOk) {
    strcpy(dum,"Ok");
  }
  else {
    strcpy(dum,"Not Ok");
    ret = FALSE;
  }
  if(bVerbose && bVerbose!=22)
    fprintf(stderr,"Angles................: %s (%5.2f%%) (%7.3f deg per angle) Worst violation: %7.3f\n"
            ,dum,100.-ang_viol/(real) nang*100,*aviol/(real) nang, aworst);
  if(dihOk) {
    strcpy(dum,"Ok");
  }
  else {
    strcpy(dum,"Not Ok");
    ret = FALSE;
  }

  if(bVerbose && bVerbose!=22)
    fprintf(stderr,"1-4 pairs.............: %s (%5.2f%%)\n",dum,100.-dih_viol/(real) ndih*100);
  if(plOk) {
    strcpy(dum,"Ok");
  }
  else {
    strcpy(dum,"Not Ok");
    ret = FALSE;
  }
  if(bVerbose && bVerbose!=22)
    fprintf(stderr,"Planarity.............: %s (%5.2f%%)\n",dum,100.-pln/(real) pl->n*100);
  if(bump) {
    strcpy(dum,"Ok");
  }
  else {
    strcpy(dum,"Not Ok");
    ret = FALSE;
  }
  if(bVerbose && bVerbose!=22) {
    
    fprintf(stderr,"Bumps.................: %s\n",dum);
    fprintf(stderr,"Chirality.............: Ok\n");

  }
  if (gyr != 0.)
  {
    gyrOk = check_gyr(al, gyr, gyr_tol, FALSE);
      if (gyrOk){
        strcpy(dum,"Ok");
      }
      else {
        strcpy(dum,"Not Ok");
        ret = FALSE;
      }
      if(bVerbose && bVerbose!=22) {
        fprintf(stderr,"Radius of Gyration....: %s\n", dum);
      }
  }

  

  if(!bondOk && bVerbose && bVerbose!=22)
    fprintf(stderr,"\nWorst: %d -- %d  lb:  %g  ub:  %g  dist:  %g (%s)\n", 
            b->at1[mem],b->at2[mem],b->lb[mem],b->ub[mem],dmem,bflag);  
  if(!angOk && bVerbose && bVerbose!=22)
    fprintf(stderr,"Angle: %d -- %d  angle:  (%g)%g  sig:  %g  dev.:  %g dist: %g lb: %g ub: %g\n\n", 
            b->at1[amem],b->at2[amem],wang,b->ang[amem],b->sig[amem],aworst,dist_ij(al->x[b->at1[amem]],al->x[b->at2[amem]]),
            b->lb[amem],b->ub[amem]);  
  
  if(!ret) return FALSE;  


  else {
    if(bVerbose!=22) {
      Rg = radius_of_gyration(al,gvec);
      sprintf(logstr,"tCNC__log_> Reached convergence after %d steps (%d ref. steps)\n", nstep, nrefst);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Sum of Violations.................: %7.4f A\n",*dviol);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Number of bond violations.........: %d\n",dist_viol);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Number of angle violations........: %d\n",ang_viol);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Number of dihedral violations.....: %d\n",dih_viol);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Bump violations...................: %7.4f A\n",bsum);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Planarity violations..............: %7.4f A\n",plsum);
      CNClog(log,logstr);
      sprintf(logstr,"tCNC__log_> Number of improper violations.....: 0\n");
      CNClog(log,logstr);      
      sprintf(logstr,"tCNC__log_> Radius of Gyration................: %6.3f A\n",Rg);
      CNClog(log,logstr);      
      sprintf(logstr,"tCNC__log_> Total number of violations........: %d\n",ntot);
      CNClog(log,logstr);
    }
    
  }
  return TRUE;
}
/*=======================================================*/

void tell_reject(FILE *out, int nit, int nviol, int plan, int nang, int nimp)
{
  char logstr[STRLEN];
  sprintf(logstr,"tCNC__log_> Structure rejected after %d steps..\n",nit);
  CNClog(out,logstr);
  sprintf(logstr,"tCNC__log_> Number of violations................: %d \n",nviol);
  CNClog(out,logstr);
  sprintf(logstr,"tCNC__log_> Number of planarity violations......: %d \n",plan);
  CNClog(out,logstr);
  sprintf(logstr,"tCNC__log_> Number of angle violations..........: %d \n",nang);
  CNClog(out,logstr);
  sprintf(logstr,"tCNC__log_> Number of improper violations.......: %d \n",nimp);
  CNClog(out,logstr);

}
/*=======================================================*/
bool atoms_moved(rvec *x, rvec *xold, int size, real rlong, real max)
{
  /* here we check how far atoms have moved and 
     whether we have to update the neighborlist
  */
  int i;
  rvec diff;
  real d;
  real cut;
  cut = sqr(rlong*10.-max);
/*   printf("cutoff is %g\n",sqrt(cut)); */
  
  
  for(i=0;i<size;i++){
    rvec_sub(x[i],xold[i],diff);
    d = norm2(diff);
/*     printf("d = %d\n",sqrt(d)); */
    
    if(d>cut) return TRUE;
  }
  return FALSE;
}
/*=======================================================*/
void copy_coords(rvec *source, rvec *target, int natoms)
{
  int i;
  for(i=0;i<natoms;i++){
    copy_rvec(source[i],target[i]);
  }
}

/*=======================================================*/

int do_disco(FILE *log,t_atomlist *al, t_bounds *b,
             t_boundtrack *bt,t_dihed *impr, t_idxgroups *imp,t_resl *rl, 
             real bump_tol, int funcnr, int dfunc,
             t_idxgroups *pl,gmx_rng_t rng,t_gridmap *gp,
/*              t_topology *top, t_forcerec *fr, */
/*              t_mdatoms *md, t_inputrec *ir, */
/*              t_groups *grps, t_nrnb *nrnb, */
/*              t_commrec *cr, t_block *cgs, t_nsborder *nsb, */
             real rlong, int nl, int maxit, int impr_func,bool bVerbose,
             bool bTarget,
             t_atomlist *tal, int targ_n, atom_id *targ_index,
             rvec *oldx, int cfreq,int npos, atom_id *posr,
             real maxff, real nbfac, bool bDoRef, int nsteps,
             real gyr, real gyr_tol)

{
#ifdef ENABLE_SIGNAL
  signal(SIGINT,signal_handler);
  signal(SIGTERM,signal_handler);
  signal(SIGUSR1,signal_handler);
#endif
  real maxf=0;
  int maxid=0;
  
  rvec xold[al->natoms];
  
  int i,k,l;
  real d;
  real mean;
  real sig;
  real ran;
  rvec diff;
  real scale;
  int nviol = 0;
  bool again = TRUE;
  bool tags[al->natoms];
  bool oldtags[al->natoms];
  
  char filename[STRLEN];
  int nit = 0;
  int r_array[b->n];
  int low;
  int high;
  int tagcount = 0;
  int nreject = 0;
  rvec mid;
  real newlen;
  real tol;
  real lb, ub;
  rvec u;
  rvec x[al->natoms];
  real angfac = .008;
  
  matrix box;
  clear_mat(box);
  int idx;
  int dih = 0;
  bool final = FALSE;
  real sumviol;
  int local;
  int nba = 0;
  int nimp_tot = 0;
  int nbond = 0;
  int nang = 0;
  int noth = 0;
  int ndih = 0;
  int ncheck = 0;
  int plan = 0;
  int nimpr=0,nbumps = 0;
  real bsum = 0;
  int ntry;
  int mem = 10000000;
  bool notouch[al->natoms];
  real plsum;
  real bviol = 0.;
  real fforce = 0.;
  real limit = (real)al->natoms/1000.;
  if (limit < 1.5) limit = 1.5;
  if (limit > 10.) limit = 10.;
  
  int nrefst = 0;  
  
  real F[FNR];
  bool bSimple = FALSE;
  int nem = 0;

  gp = nb(al, gp, 9);  
  fill_neighborlist(log, al, gp, 9.,UPDATE_NONBONDED);  
  gp = reset_grid(gp);
  
 

/*    int xtc = open_xtc("progress.xtc","w");  */
/*    int xtc_count = 1;  */
  
/*   FILE *of = ffopen("conv.log","w"); */
  



  for(k=0;k<al->natoms;k++){
    notouch[k] = al->isposres[k];
  }

  if(bTarget){
    for(k=0;k<targ_n;k++){
      notouch[targ_index[k]] = TRUE;
      copy_rvec(tal->x[k],al->x[targ_index[k]]);
    }
    /* shift the system */
     shift_system(al,notouch);      
 
  }


  for(k=0;k<b->n;k++){
    if(b->isbond[k] ||
       b->isang[k] || 
       b->isdih[k]) nba++;
  }

  
  for(i=0;i<al->natoms;i++){
    tags[i] = TRUE;
  }

  /* write_pdb(al,"start.pdb"); */
  
  /*  CHECK LOOP */

/*     write_pdb(al,"start.pdb");         */
  copy_coords(al->x,xold,al->natoms);

  while(again){

    
    nit++;
    if(nit > maxit) {
          
        tell_reject(log,nit,nviol,plan,nang,nimpr);
        return 0;
    }


    if(nit % 10 ==0 && gyr != 0.) {
      bool b_gyrOk = check_gyr(al,gyr,gyr_tol, bVerbose);
      if(!b_gyrOk) {
        for(i=0;i<al->natoms;i++){
          tags[i] = TRUE;
        }
      }
    }
    


    sumviol = 0;
    rand_array_int(rng,b->n,r_array);      
    for(k=0;k<al->natoms;k++){
      oldtags[k] = tags[k];
    }
    
    for(k=0;k<al->natoms;k++){
      tags[k] = FALSE;
    }

    if(bTarget){
      for(k=0;k<targ_n;k++){
        copy_rvec(tal->x[k],al->x[targ_index[k]]);
      }
    }
    
    for(k=0;k<al->natoms;k++){
      if(al->isposres[k]) copy_rvec(oldx[k],al->x[k]);
    }
/*     write_pdb(al,"start2.pdb");     */
/*      exit(0);  */
    




    
    nviol = 0;
    bool end = FALSE;
    i=0;

    nbond = nang = noth = ndih = 0;
    ncheck = 0;
    int barray[nba];
    int bcount = 0;

    while(!end){
      idx = r_array[i];
      if(oldtags[b->at1[idx]] || oldtags[b->at2[idx]]){
        ncheck++;
        
        if(b->isang[idx]){
          
          if(!check_angle(al->x[b->at3[idx]],al->x[b->at1[idx]],al->x[b->at2[idx]],b,idx,
                          notouch[b->at1[idx]],notouch[b->at2[idx]])){
            nang++;nviol++;
            tags[b->at1[idx]] = tags[b->at2[idx]]  = tags[b->at3[idx]] = TRUE;
/*              al2xtc(al,xtc,xtc_count);   */
/*                xtc_count++;   */
            
          }
        }
        rvec_sub(al->x[b->at2[idx]],al->x[b->at1[idx]],diff);

        d = norm2(diff);
        if(d > sqr(b->ub[idx]) || d < sqr(b->lb[idx])){
          d = norm(diff);
          if(b->isbond[idx]) nbond++;
          else if(b->isdih[idx]) ndih++;
          else if(b->isang[idx]) nang++;
          else noth++;
          
          if(d < b->ub[idx]){
            sumviol+=b->ub[idx] - d;
          }else {
            sumviol+=d - b->lb[idx];
          }
          correct_dist(al->x[b->at1[idx]],al->x[b->at2[idx]],
                       diff,d,b,idx,dfunc,rng,
                       notouch[b->at1[idx]],notouch[b->at2[idx]]);
          tags[b->at1[idx]] = TRUE;
          tags[b->at2[idx]] = TRUE;
/*            al2xtc(al,xtc,xtc_count);  */
/*            xtc_count++;  */
          
          nviol++;
        }
      }
      i++;
      if(i == b->n){
        end = TRUE;
      }
    }
  

    if(nit % cfreq == 0) {
      if(bVerbose)
        fprintf(stderr,"nviol = %d mem = %d\n",nviol,mem);

      if(nviol > mem) {
        tell_reject(log,nit,nviol,plan,nang,nimpr);
        return 0;
      }
      else 
        mem = nviol;
    }
/*             al2xtc(al,xtc,xtc_count);   */
/*             xtc_count++;   */

/*      if (GotTermSignal) {  */
/*        write_pdb(al,"sig.pdb"); */
/*        exit(0); */
/*      } */
     
    if(GotTermSignal == 1){
      GotTermSignal = 0;
      write_pdb(al,"sig.pdb");
      
      return 0;
    }
    else if(GotTermSignal > 1){
      exit(0);
    }
  
    



    /* END CHECK LOOP */

    nimpr=0;
    nbumps = 0;

   /* check improps */
    
    
    int imps[imp->n];
    nimp_tot = impr->n;
    int tmp = 0;
    
    if(impr_func == 1){
      tmp = check_impr(al,imp,imps);
      if(tmp > imp->n/2. && !bTarget && !npos) mirror(al);
       nimpr = do_improp(al,imp,tags,rng); 
/*       nimpr = do_chiral(al,imp,tags); */
    }
    
    else if(impr_func == 2){
      tmp = check_impr(al,imp,imps);
      nimpr = do_improp2(al,impr,rl,tags);         
      
      if(nimpr > impr->n/2. && !bTarget && !npos) {
        mirror(al);
      }
    }
/*     al2xtc(al,xtc,xtc_count);   */
/*     xtc_count++;   */
    
    
    nviol+=nimpr;
    

    /* PLANAR GROUPS */
    

    plan = 0;

    int plarr[pl->n];
    rand_array_int(rng,pl->n,plarr);
    int xxx;
    real planx;
    real plmem = 0;
    rvec mem[al->natoms];
    rvec mem2[al->natoms];
    bool plOk;
    real plsum,plmax,plworst;
    int plcount;

/*     plOk = count_planar(al, pl, &plsum, */
/*                         &plcount, 1.1, &plworst); */
/*     printf("nplanviol = %d\n",plcount); */
    
    

    for(k=0;k<pl->n;k++){
      plan += check_planar(al,pl->natoms[plarr[k]],pl->atoms[plarr[k]],
                            pl->val[plarr[k]]);    
      
    }

/*     al2xtc(al,xtc,xtc_count);   */
/*     xtc_count++;   */

    plOk = count_planar(al, pl, &plsum,
                        &plcount, 5., &plworst);
/*     printf("nplanviol = %d\n",plcount); */
/*     exit(0); */
    plan = plcount;
    
    nviol+=plan;







/*     fprintf(of,"%10d %10d %g %10d %10d %10d %10d %10d\n", nit,nviol,sumviol,nimpr,plan,nbond,nang,ndih); */
    

    if(((nviol <= nba/80.) && (plan <= pl->n/80.))  
       || (nl && (nit % nl == 0)))     
    
    
    { 

    /* FINAL CHECKS */
    /*=================================================*/
      
      real dviol,aviol,worst;
      real lastf;
      
      if(atoms_moved(al->x,xold,al->natoms,rlong,4.))
      {
        if(bVerbose)
          printf("Update nl at step %d\n",nit);
        gp = nb(al, gp,9);
        fill_neighborlist(log, al, gp, 9.,UPDATE_NONBONDED);  
        gp = reset_grid(gp);


/*         update_neighborlist(log,al, x, top,md, fr,ir,   */
/*                             grps, nrnb, cr, cgs,nsb, */
/*                             rlong,bSimple,npos,posr);   */
        copy_coords(al->x,xold,al->natoms);
      }
      

      fforce = 0;
      maxf = 0.;
      maxid = -1;
      
      
      nimpr = check_impr(al,imp,imps);
      
      
      if(nimpr == 0 &&  cross_check_dist(log,al,b,&dviol,&aviol,
                                         rng,pl,bump_tol,bt,
                                         bVerbose,nit,nrefst, gyr, gyr_tol)) { 
        final = TRUE; 
      }
      else { 
        do_refine(al,b,bump_tol,nbfac,1.,pl,imp,bt,tags, 
                   &maxid,&maxf,F,TRUE,bVerbose,rng,angfac);  
        nimpr = check_impr(al,imp,imps); 
        nimpr = do_improp(al,imp,tags,rng); 
/*         nimpr = do_chiral(al,imp,tags); */
        
        nrefst++; 
        bool bMin = TRUE;
        if ((maxf < limit*10) && bDoRef  /* && nimpr == 0 */ ) {
/*           printf("maxf %g (limit = %g)  numpr = %d\n",maxf,limit,nimpr); */
          int xx = 1;
          do {
            if(atoms_moved(al->x,xold,al->natoms,rlong,4.))
            {
              if(bVerbose)
                printf("Update nl at refstep %d\n",nrefst);
              gp = nb(al, gp, 9);
              fill_neighborlist(log, al, gp, 9.,UPDATE_NONBONDED);
              gp = reset_grid(gp);

/*               update_neighborlist(log,al, x, top,md, fr,ir,   */
/*                                   grps, nrnb, cr, cgs,nsb, */
/*                                   rlong,bSimple,npos,posr); */
              copy_coords(al->x,xold,al->natoms);
            }
            
            do_refine(al,b,bump_tol,nbfac*.1,.1,pl,imp,bt,tags,&maxid,
                      &maxf,F,TRUE,bVerbose,rng,angfac); 
            nrefst++;
            if(impr_func == 1){
              tmp = check_impr(al,imp,imps);
              if(tmp > imp->n/2. && !bTarget && !npos) mirror(al);
               nimpr = do_improp(al,imp,tags,rng); 
/*               nimpr = do_chiral(al,imp,tags); */
            }
            else if(impr_func == 2){
              tmp = check_impr(al,imp,imps);
              nimpr = do_improp2(al,impr,rl,tags);         
              if(nimpr > impr->n/2. && !bTarget && !npos) {
                mirror(al);
              }
            }
            if(bVerbose){
              fprintf(stderr,"FTOT = %8.3f FB = %8.3f FANG = %8.3f FPL = %8.3f FBO = %8.3f FNB = %8.3f FCHI = %8.3f CHI = %8d LEVEL = %2d\n",
                      F[FTOT],F[FB],F[FANG],F[FP],F[FBO],F[FNB],F[FCHI],tmp,xx);
            }
            if (maxf > 4.) { 
              bMin = FALSE; 
              break; 
            } 
            xx++;
          }
          while (xx < 11);
          xx=1;
          if(bMin) {
            lastf = maxf;
            
            do {
              if(atoms_moved(al->x,xold,al->natoms,rlong,4.))
              {
                if(bVerbose)
                  printf("Update nl at refstep %d\n",nrefst);
                gp = nb(al, gp, 9);

                fill_neighborlist(log, al, gp, 9.,UPDATE_NONBONDED);  
                gp= reset_grid(gp);
                
/*                 update_neighborlist(log,al, x, top,md, fr,ir,   */
/*                                     grps, nrnb, cr, cgs,nsb, */
/*                                     rlong,bSimple,npos,posr);  */
                copy_coords(al->x,xold,al->natoms); 
              }
              do_refine(al,b,bump_tol,nbfac*.1,.1,pl,imp,bt,tags,
                        &maxid,&maxf,F,TRUE,bVerbose,rng,angfac); 
              nrefst++;
              if(impr_func == 1){
                tmp = check_impr(al,imp,imps);
                if(tmp > imp->n/2. && !bTarget && !npos) mirror(al);
/*                 nimpr = do_chiral(al,imp,tags); */
                 nimpr = do_improp(al,imp,tags,rng); 
              }
              else if(impr_func == 2){
                tmp = check_impr(al,imp,imps);
                nimpr = do_improp2(al,impr,rl,tags);         
                if(nimpr > impr->n/2. && !bTarget && !npos) {
                  mirror(al);
                }
              }
              if(bVerbose){
                fprintf(stderr,"FTOT = %8.3f FB = %8.3f FANG = %8.3f FPL = %8.3f FBO = %8.3f FNB = %8.3f FCHI = %8.3f CHI = %8d LEVEL = %2d\n",
                        F[FTOT],F[FB],F[FANG],F[FP],F[FBO],F[FNB],F[FCHI],tmp,xx);

              }
               if (maxf > 4.) { 
                 bMin = FALSE; 
                 break; 
               }
               if(tmp == 0 && lastf > maxf && F[FTOT] > .1) {
                 xx--;
                 lastf = maxf;
               }
               
               
               xx++;
            }
            while (xx < nsteps);
          }
          do_refine(al,b,bump_tol,nbfac,1.,pl,imp,bt,tags,&maxid,&maxf,F,FALSE,bVerbose,rng,angfac); 
/*           if(nit == 1000) { */
/*             write_pdb(al,"step1000.pdb"); exit(0); */
/*           } */
          
          
          if(!check_impr(al,imp,imps) && cross_check_dist(log,al,b,&dviol,&aviol,
                                                          rng,pl,bump_tol,bt,
                                                          bVerbose,nit,nrefst, gyr, gyr_tol) )
          {
            return TRUE;
          }
        }
/*          else  */
/*          {  */
/*            apply_f(al,1,tags);  */
/*            clear_forces(al);  */
/*          }  */
        
/*         al2xtc(al,xtc,xtc_count);   */
/*         xtc_count++;   */

        final = FALSE;
      }
      if(final){
        again = FALSE;
        break;
      }
      else{
        for(k=0;k<al->natoms;k++){
          tags[k] = TRUE;
        }
        final = TRUE;
      }
    }
    
    tagcount = 0;
    for(k=0;k<al->natoms;k++){
      tagcount+=tags[k];
    }
    
    if(bVerbose && nit % 10 == 0){
      real Rg;
      rvec gvec;
      Rg = radius_of_gyration(al,gvec);

      fprintf(stderr,"st=%d nv = %d sum: %g im = %d (%d) b= %d a= %d dih= %d ot= %d pln= %d chk= %d(%d) f = %g maxf=%g (%s %d %s) Rg = %g\n",  
              nit,nviol,sumviol,nimpr,imp->n,nbond,nang,ndih,noth,plan,ncheck,b->n,F[FTOT],maxf,al->name[maxid],al->resid[maxid],al->resname[maxid], Rg);  
      fflush(stderr);
    }
    
  } 
/*   printf("\n"); */
/*     al2xtc(al,xtc,xtc_count);  */
/*     xtc_count++;  */
/*    close_xtc(xtc);  */
/*     exit(0);  */
  
  return TRUE;
}

/*=======================================================*/

bool ligand_disco(t_atomlist *al, t_bounds *b, t_idxgroups *imp, t_idxgroups *pl,
                  int maxit, gmx_rng_t rng, FILE *log, t_boundtrack *bt, bool bVerbose)
{
  /* optimized disco routine for small molecules.
     We don't do neighborsearching since 
     all atoms have the full atom set in their
     neighborlist.
  */

  bool again,end;
  int nit, ncheck;
  real sumviol,d,ang;
  int i,k,idx;
  rvec diff;
  bool tags[al->natoms],oldtags[al->natoms];
  int r_array[b->n];
  bool notouch[al->natoms];
  int imps[imp->n];
  int imp_count;
  int plarr[pl->n];
  int plan;
  int maxid;
  real maxf;
  real F[FNR];
  int nrefst = 0;
  real limit = (real)al->natoms/1000.;
  bool bMin;
  real dviol, aviol;
  int tagcount = 0;
  bool final = FALSE;
  int nviol;
  int nimpr;
  
    
  if (limit < 1.5) limit = 1.5;
  if (limit > 10.) limit = 10.;
  

  nit = 0;
  again = TRUE;
    
  
  for(i=0;i<al->natoms;i++) notouch[i]=FALSE;
  
  while(again){
    nit++;
    if(nit > maxit) {
      return FALSE;
    }
    sumviol = 0;
    rand_array_int(rng,b->n,r_array);      
    for(k=0;k<al->natoms;k++){
      oldtags[k] = tags[k];
    }
    
    for(k=0;k<al->natoms;k++){
      tags[k] = FALSE;
    }
    
    nviol = 0;
    end = FALSE;
    i=0;
    ncheck = 0;
    
    /* do distance and angle constraints */
    
    while(!end) {
      idx = r_array[i];
      if(oldtags[b->at1[idx]] || oldtags[b->at2[idx]]){
        ncheck++;
        if(b->isang[idx]){
          if(!check_angle(al->x[b->at3[idx]],al->x[b->at1[idx]],al->x[b->at2[idx]],b,idx,
                          notouch[b->at1[idx]],notouch[b->at2[idx]])){
            nviol++;
            tags[b->at1[idx]] = tags[b->at2[idx]]  = tags[b->at3[idx]] = TRUE;
          }
        }
        rvec_sub(al->x[b->at2[idx]],al->x[b->at1[idx]],diff);
        
        d = norm2(diff);
        if(d > sqr(b->ub[idx]) || d < sqr(b->lb[idx])){
          d = norm(diff);
          if(d < b->ub[idx]){
            sumviol+=b->ub[idx] - d;
          }else {
            sumviol+=d - b->lb[idx];
          }
          correct_dist(al->x[b->at1[idx]],al->x[b->at2[idx]],
                       diff,d,b,idx,0,rng,
                       notouch[b->at1[idx]],notouch[b->at2[idx]]);
          tags[b->at1[idx]] = TRUE;
          tags[b->at2[idx]] = TRUE;
          nviol++;
        }
      }
      i++;
      if(i == b->n){
        end = TRUE;
      }
    }
    
    /* check chirality */
    
    imp_count = check_impr(al,imp,imps);
/*     if(imp_count > imp->n/2.) mirror(al); */
    nimpr = do_improp(al,imp,tags,rng);
    
    /* check planar groups */
    
    rand_array_int(rng,pl->n,plarr);
    plan = 0;
    for(k=0;k<pl->n;k++){
      plan += check_planar(al,pl->natoms[plarr[k]],pl->atoms[plarr[k]],
                           pl->val[plarr[k]]);    
      
    }
    /* do nonbonded stuff */
    if(nit % 10 ==0){
      nimpr = check_impr(al,imp,imps);
      if(nimpr == 0 &&  cross_check_dist(log,al,b,&dviol,&aviol,
                                         rng,pl,.3,bt,
                                         22,nit,nrefst, 0, 0)) { 
        final = TRUE; 
      }
      else { 
        do_refine(al,b,.3,.5,1.,pl,imp,bt,tags,
                  &maxid,&maxf,F,TRUE,bVerbose,rng,0.0012); 
        
        nrefst++;
        bMin = TRUE;
        
        if (maxf < limit){
          int xx = 1;
          do {
            do_refine(al,b,.3,.5*.1*xx,.1*xx,pl,imp,bt,tags,&maxid,
                      &maxf,F,TRUE,bVerbose,rng,0.0012); 
            nrefst++;
            imp_count = check_impr(al,imp,imps);
/*             if(imp_count > imp->n/2.) mirror(al); */
            nimpr = do_improp(al,imp,tags,rng);
            if (maxf > 3.) {
              bMin = FALSE;
              break;
            }
            xx++;
          }
          while (xx < 11);
          xx=1;
          if(bMin) {
            do {
              do_refine(al,b,.3,.5,1.,pl,imp,bt,tags,
                        &maxid,&maxf,F,TRUE,bVerbose,rng,0.0012); 
              
              nrefst++;
              imp_count = check_impr(al,imp,imps);
/*               if(imp_count > imp->n/2.) mirror(al); */
              nimpr = do_improp(al,imp,tags,rng);
              if (maxf > 3.) {
                bMin = FALSE;
                break;
              }
              xx++;
            }
            while (xx < 20);
          }
          do_refine(al,b,.3,.5,1.,pl,imp,bt,tags,&maxid,&maxf,F,FALSE,bVerbose,rng,0.0012); 
          
          if(!check_impr(al,imp,imps) && cross_check_dist(log,al,b,&dviol,&aviol,
                                                          rng,pl,.3,bt,
                                                          22,nit,nrefst,0, 0))
          {
            return TRUE;
          }
        }
        final = FALSE;
      }
      if(final){
        again = FALSE;
        break;
      }
      else{
        for(k=0;k<al->natoms;k++){
          tags[k] = TRUE;
        }
        final = TRUE;
      }
    }
    tagcount = 0;
    for(k=0;k<al->natoms;k++){
      tagcount+=tags[k];
    }
  }
  return TRUE;
}
  

  
     
