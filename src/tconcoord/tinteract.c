#include <tconcoord.h>


int main(int argc, char **argv)
{
  
   static char *desc[] = {
     "Visualise interactions between receptor and ligand"
   };
   
   
  t_atomlist *al;
  t_nblist   *nlist;
  t_contab *ct;
  t_resl *rl;
  t_vdw *vdw;
  t_vdwcomb *vdwcomb;
  t_types *tp;
  t_bounds *b;
  t_dihed *impr;
  t_forcerec *fr = NULL;
  t_topology top;
  t_atoms *atoms = NULL; 
  t_nsborder *nsb;
  t_mdatoms  *md;
  t_block    *cgs;
  t_inputrec *ir;
  t_nrnb     nrnb;
  t_commrec  *cr;
  t_groups   *grps;
  char title[STRLEN];

  int        *cg_index;
  int        m,natoms;
  ivec       *nFreeze;
  rvec       box_size;
  real       lambda=0,dvdlambda=0;
  real       rlong = 0.8;
  rvec box_space;
  rvec *x;
  matrix box;
  int natom;
  int nalt;
  bool bVerbose = FALSE;
  bool bLig = FALSE;
  int i,k,j;
  real plan;

  atom_id *ndx1;
  atom_id *ndx2;
  int nat1;
  int nat2;
  char *names;
  

  t_pargs pa[] = {
    { "-v",   FALSE, etBOOL, {&bVerbose},"Make Noise"}
  };


  t_filenm fnm[] = {
    { efPDB,"-s","protein", ffREAD }, 
    { efNDX,"-n","index", ffREAD},
    { efPDB,"-o","interact", ffWRITE}, 
    { efPDB,"-op","out", ffWRITE}, 
    { efDAT,"-od","pr_lg", ffWRITE} 

  };


  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,0,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

/*   pdb = pdb2pdbrecord(opt2fn("-f",NFILE,fnm)); */
/*   if(pdb->n > 1) */
/*   { */
/*     fprintf(stderr,"Warning! Input file contains more than 1 model\n"); */
/*     fprintf(stderr,"---> Will take first model.\n"); */
/*   } */

  read_tps_conf(opt2fn("-s",NFILE,fnm),title,&top,&x,NULL,
                box,TRUE);
/*   bool bTop=fn2bTPX(opt2fn("-s",NFILE,fnm));  */
  
  al = al_from_atoms(&(top.atoms),x);
  
/*   al = pdb->al[0]; */

  fprintf(stderr,"Select first group:\n");
  
  get_index(atoms,opt2fn("-n",NFILE,fnm),1,&nat1,&ndx1,&names);
  
  fprintf(stderr,"Select second group:\n");
  
  get_index(atoms,opt2fn("-n",NFILE,fnm),1,&nat2,&ndx2,&names);

  get_symbol(stdout,al);

  ct = contab_init();
  ct = contab_realloc(ct,al->natoms);
  
  nb_search(al,ct,0.8,1);
  rename_at(al);
  
  get_order(al);

    
  vdwcomb = read_vdwcomb(1,FALSE); 
  t_vdwcomb *vdw14comb = read_vdwcomb(1,TRUE);
  tp = read_atom_types(1); 
  get_types(stdout,al,tp); 
  vdw = read_vdw_radii(1);   


  get_vdw_radii(stdout,al,vdw,vdwcomb,vdw14comb); 
  get_hybrid2(al,tp);
  renumber_atoms(al,1);
  renumber_residues(al,1);


  t_idxgroups *don = idx_init();
  t_idxgroups *acc = idx_init();
  don = idx_realloc(don,1);
  acc = idx_realloc(acc,1);
  don->n = acc->n = 1;
  hbonds(al,don,acc);
  
  t_idxgroups *phob = idx_init();
  phob = idx_realloc(phob,1);
  FILE *phobics = cnclib("Hydrophobic.dat");
  hydrophobics(al,phob,phobics);


  t_idxgroups *npho = idx_init();
  npho = idx_realloc(npho,1);
  t_idxgroups *ndon = idx_init();
  ndon = idx_realloc(ndon,1);
  t_idxgroups *nacc = idx_init();
  nacc = idx_realloc(nacc,1);
  
  /* fake charges */
  
  fake_charge(al);

  t_idxgroups *pos = idx_init();
  pos = idx_realloc(pos,1);

  t_idxgroups *neg = idx_init();
  neg = idx_realloc(neg,1);

  t_idxgroups *npos = idx_init();
  npos = idx_realloc(npos,1);

  t_idxgroups *nneg = idx_init();
  nneg = idx_realloc(nneg,1);
  
  for(i=0;i<al->natoms;i++){
    if(al->q[i]==-1.) add_to_group(neg,0,i);
    else if(al->q[i] == 1.) add_to_group(pos,0,i);
  }

  for(i=0;i<al->natoms;i++){
    al->ngrps[i] = 1;
    snew(al->grpnr[i],1);
    for(k=0;k<nat1;k++){
      if(i == ndx1[k]) al->grpnr[i][0] = 0;
    }
    for(k=0;k<nat2;k++){
      if(i == ndx2[k]) al->grpnr[i][0] = 1;
    }
  }
  
  /* get hbonds */

  rvec diff;
  real dist,angle;
  real lb,ub;
  int idx,at1,at2;
  char name[STRLEN];
  FILE *od = ffopen(opt2fn("-od",NFILE,fnm),"w");
  FILE *of = ffopen(opt2fn("-o",NFILE,fnm),"w");
  fprintf(od,";Possible protein-ligand constraints\n;===================================\n; protein acceptor / ligand donor combinations\n[ pa/ld ]\n");
  
  int count = 1;
  
  for(i=0;i<nat1;i++){
    idx = ndx1[i];
    
    if(al->isacc[idx]){
      for(k=0;k<al->nnb[idx];k++){
        at1 = al->nb[idx][k];
        if(al->grpnr[at1][0] == 1 && al->isdon[at1]){
          rvec_sub(al->x[idx],al->x[at1],diff);
          dist = norm(diff);
          if(dist < 2.5){
            at2 = al->bonds[at1][0];
            angle = RAD2DEG*angle_ij_ik(al->x[at1],al->x[idx],al->x[at2]);
            if(angle > 130.){
              strcpy(name," ACC");
              fprintf(of,pdbstr,count,name," ",al->resname[idx],
                      al->chain[idx],al->resid[idx],al->x[idx][0],
                      al->x[idx][1],al->x[idx][2],1.0,0.0);
              fprintf(of,pdbstr,count+1,name," ",al->resname[at1],
                      al->chain[at1],al->resid[at1],al->x[at1][0],
                      al->x[at1][1],al->x[at1][2],1.0,0.0);
              fprintf(of,"CONECT%5d%5d\n",count,count+1);
              fprintf(od,"%8d %8d %9.4f %9.4f %9.4f\n",al->id[idx],al->id[at1]-nat1,1.9,1.7,2.5);
              count+=2;
            }
          }
        }
      }
    }
  }
  fprintf(od,"; protein donor / ligand acceptor combinations\n[ pd/la ]\n");
  
  for(i=0;i<nat1;i++){
    idx = ndx1[i];
    if(al->isdon[idx]){
      for(k=0;k<al->nnb[idx];k++){
        at1 = al->nb[idx][k];
        if(al->grpnr[at1][0] == 1 && al->isacc[at1]){
          rvec_sub(al->x[idx],al->x[at1],diff);
          dist = norm(diff);
          if(dist < 2.5){
            at2 = al->bonds[idx][0];
            angle = RAD2DEG*angle_ij_ik(al->x[idx],al->x[at1],al->x[at2]);
            if(angle > 130.){
              strcpy(name," DON");
              fprintf(of,pdbstr,count,name," ",al->resname[idx],
                      al->chain[idx],al->resid[idx],al->x[idx][0],
                      al->x[idx][1],al->x[idx][2],1.0,0.0);
              fprintf(of,pdbstr,count+1,name," ",al->resname[at1],
                      al->chain[at1],al->resid[at1],al->x[at1][0],
                      al->x[at1][1],al->x[at1][2],1.0,0.0);
              fprintf(of,"CONECT%5d%5d\n",count,count+1);
              fprintf(od,"%8d %8d %9.4f %9.4f %9.4f\n",al->id[idx],al->id[at1]-nat1,1.9,1.7,2.5);
              count+=2;
            }
          }
        }
      }
    }
  }
  fprintf(od,"; Hydrophobic combinations\n[ php/lhp ]\n");
  for(i=0;i<nat1;i++){
    idx = ndx1[i];
    if(al->ishphob[idx]){
      for(k=0;k<al->nnb[idx];k++){
        at1 = al->nb[idx][k];
        if(al->grpnr[at1][0] == 1 && al->ishphob[at1]){
          rvec_sub(al->x[idx],al->x[at1],diff);
          dist = norm(diff);
          if(dist < 4.){
/*             at2 = al->bonds[at1][0]; */
/*             angle = RAD2DEG*angle_ij_ik(al->x[at1],al->x[idx],al->x[at2]); */
/*             if(angle > 130.){ */
            strcpy(name," PHO");
            fprintf(of,pdbstr,count,name," ",al->resname[idx],
                    al->chain[idx],al->resid[idx],al->x[idx][0],
                    al->x[idx][1],al->x[idx][2],1.0,0.0);
            fprintf(of,pdbstr,count+1,name," ",al->resname[at1],
                    al->chain[at1],al->resid[at1],al->x[at1][0],
                    al->x[at1][1],al->x[at1][2],1.0,0.0);
            fprintf(of,"CONECT%5d%5d\n",count,count+1);
            lb = 3.2;
            if (dist < lb) lb = dist*.9;
            ub = 5.5;
            fprintf(od,"%8d %8d %9.4f %9.4f %9.4f\n",al->id[idx],al->id[at1]-nat1,dist,lb,ub);
            count+=2;
/*             } */
          }
        }
      }
    }
  }


  write_pdb(al,opt2fn("-op",NFILE,fnm));
  
  return 0;
}
