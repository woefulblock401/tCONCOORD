#include <tconcoord.h>
/*==============================================================*/
int write_bonds(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                t_input *inp,t_excl *ex,t_bondlib *bl, bool bIgn)

/* write the bonds to dist file */

{
  int i,k;
  int counter = 0;
  
  real di;
  real av,lb,ub;
  fprintf(fp,"[ BONDS ]\n");
  char warn[STRLEN];
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      if(i > al->bonds[i][k] && !is_excluded(ex,al,i,al->bonds[i][k]) 
         && (!al->isposres[i] || !al->isposres[al->bonds[i][k]]))
      {
        add_bound(bt,i,al->bonds[i][k]); 
        di = DIST(al,i,al->bonds[i][k]);
        if(!get_bond(al->type[i],al->type[al->bonds[i][k]],bl,&av,&lb,&ub)){
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f ; from input\n",al->id[i],al->id[al->bonds[i][k]],di,di-inp->bond_tol
                  ,   di+inp->bond_tol);
          sprintf(warn,"Using input geometry for bond %s - %s\n",al->type[i],al->type[al->bonds[i][k]]);
          CNCwarn(log,warn);
          counter++;
        }
        else{
          if(!bIgn){
            lb = av - inp->bond_tol;
            ub = av + inp->bond_tol;
            if(di < lb){
              lb = di*0.99;
            }
            if(di > ub){
              ub = di*1.01;
            }
          }
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[i],al->id[al->bonds[i][k]],av,lb,ub);
          counter++;
        }
        
      }
    }
  }
  sprintf(warn,"tCNC__log_> Wrote %10d bond constraints\n",counter);
  CNClog(log,warn);
  return counter;
}
/*=============================================================*/
int write_bonds_from_tpx(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                         t_input *inp,t_excl *ex,t_topology *top, bool bIgn)
{
  int i,k;
  int at1,at2,p;
  real d,av,lb,ub,tol;
  int counter = 0;
  char logstr[STRLEN];
  
  t_iparams *ip = top->idef.iparams; 
  tol = inp->bond_tol;
  
  fprintf(fp,"[ BONDS ]\n");
  
  for(i=0;i<F_NRE;i++){ 
    if(IS_CHEMBOND(i)  && top->idef.il[i].nr ){ 
      int nratoms=interaction_function[i].nratoms; 
      k=0; 
      while(k<top->idef.il[i].nr){ 
        p =  top->idef.il[i].iatoms[k];
        at1 = top->idef.il[i].iatoms[k+1]; 
        at2 = top->idef.il[i].iatoms[k+2]; 
        if ((!al->isposres[at1] || !al->isposres[at2]) 
            && !is_excluded(ex,al,at1,at2)) {
          av = ip[p].bham.a*10.;
          lb = av-tol;
          ub = av+tol;
          
          if(!bIgn){
            d = DIST(al,at1,at2);
            if(d > ub) ub = d*1.01;
            else if(d < lb) lb = d*0.99;
          }
          add_bound(bt,at1,at2); 
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],al->id[at2],av,lb,ub);
          counter++;
        }
        k+=nratoms+1; 
      } 
    } 
  } 
  sprintf(logstr,"tCNC__log_> Wrote %10d bond constraints\n",counter);
  CNClog(log,logstr);

  return counter;
}




/*=============================================================*/
int write_angles(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                 t_input *inp,t_excl *ex, t_bondlib *bl, bool bIgn)

/* write angles to dist file */

{
  int counter = 0;
  fprintf(fp,"[ ANGLES ]\n");
  int i,k,j;
  real av,lb,ub,avang,sig;
  real di;
  real angle;
  int at1,at2;
  real tol = inp->angle_tol;
  real sigma = inp->angle_sig;
  char warn[STRLEN];
  

  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nbonds[i];k++){
      at1 = al->bonds[i][k];
      for(j=0;j<al->nbonds[at1];j++){
        at2 = al->bonds[at1][j];
        if(i > at2 && !is_excluded(ex,al,i,at2) && (!al->isposres[i] || !al->isposres[at1] || !al->isposres[at2]))
        {
          add_bound(bt,i,at2); 
          di = DIST(al,i,at2);
          angle = RAD2DEG*angle_ij_ik(al->x[at1],al->x[i],al->x[at2]);
          if(!get_angle(al->type[i],al->type[at1],al->type[at2],bl,&av,&lb,&ub,
                        &avang,&sig)){
            fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f %8d %9.4f %9.4f\n",
                    al->id[i],al->id[at2],di,di-tol,di+tol,al->id[at1],angle,sigma);
            sprintf(warn,"Using input geometry for angle %s - %s - %s\n",al->type[i],al->type[at1],al->type[at2]);
            CNCwarn(log,warn);
            counter++;
          }


          else{
/*             sig = sig+sigma; */
            if(!bIgn){
              sig = sigma;
              if(di < lb){
                lb = di*0.99;
              }
              if(di > ub){
                ub = di*1.01;
              }
              if(angle < avang-sig){
                sig = avang - angle+0.1;
                sprintf(warn,"Unusual angle %s-%s-%s (%d-%d-%d) %8.3f (%8.3f)\n",
                        al->name[i],al->name[at1],al->name[at2],i+1,at1+1,at2+1,
                        angle,avang);
                CNCwarn(log,warn);
                
              }
              if(angle > avang+sig){
                sig = angle - avang+0.1;
                sprintf(warn,"Unusual angle %s-%s-%s (%d-%d-%d) %8.3f (%8.3f)\n",
                        al->name[i],al->name[at1],al->name[at2],i+1,at1+1,at2+1,
                        angle,avang);
                CNCwarn(log,warn);

              }
              
            }
            fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f %8d %9.4f %9.4f\n",
                    al->id[i],al->id[at2],av,lb-tol,ub+tol,al->id[at1],avang,sig);
            counter++;
          }
        }
      }
    }
  }
  sprintf(warn,"tCNC__log_> Wrote %10d angle constraints\n",counter);
  CNClog(log,warn);

  return counter;
}
/*=============================================================*/
int write_angles_from_tpx(FILE *fp, FILE *log,t_atomlist *al, 
                          t_boundtrack *bt, t_input *inp,t_excl *ex, 
                          t_topology *top, bool bIgn)
{

  int i,k,na;
  int at1,at2,at3,p;
  real d,avang,lb,ub,tol,sig,sigma,angle;
  int counter = 0;
  char logstr[STRLEN];
  
  t_iparams *ip = top->idef.iparams; 
  t_iatom *ia = top->idef.il[F_ANGLES].iatoms; 
  
  na = top->idef.il[F_ANGLES].nr;  

  tol = inp->angle_tol;
  sigma = inp->angle_sig;
  
  fprintf(fp,"[ ANGLES ]\n");
  for(i=0;i<na;i+=4){ 
    avang = ip[ia[i]].bham.a;
    at1 = ia[i+1]; 
    at2 = ia[i+2]; 
    at3 = ia[i+3]; 
    if(!is_excluded(ex,al,at1,at3) && (!al->isposres[at1] 
                                       || !al->isposres[at2] 
                                       || !al->isposres[at3]))
    {
      
      d = DIST(al,at1,at3);
      lb = d - tol;
      ub = d + tol;
      angle = RAD2DEG*angle_ij_ik(al->x[at2],al->x[at1],al->x[at3]);
      sig = sigma;      
      if(!bIgn){
        if(d < lb){
          lb = d*0.99;
        }
        if(d > ub){
          ub = d*1.01;
        }
        if(angle < avang-sig){
          sig = avang - angle+0.1;
        }
        if(angle > avang+sig){
          sig = angle - avang+0.1;
        }
      }
      add_bound(bt,at1,at3); 
      fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f %8d %9.4f %9.4f\n",
              al->id[at1],al->id[at3],d,lb,ub,al->id[at2],avang,sig);
      counter++;
      
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d angle constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}






/*=============================================================*/
int write_hbond_angles(FILE *fp, FILE *log, t_atomlist *al, 
                              t_boundtrack *bt, t_excl *ex)
{
  int i,k;
  int don,acc,don_b, acc_b;
  real di,a;
  real ang = 160;
  real sig = 30;
  real lb, ub;
  int counter = 0;
  real angle;
  char logstr[STRLEN];
  
  for(i=0;i<al->nhbonds;i++){
    if(al->hbonds[i]->type == ehBBBB && al->hbonds[i]->constr == TRUE){
      sig = 30;
      
      don = al->hbonds[i]->don;
      acc = al->hbonds[i]->acc;
      don_b = al->hbonds[i]->don_b;
      acc_b = al->hbonds[i]->acc_b;
      di = al->hbonds[i]->db_a_dist;
      angle = al->hbonds[i]->angle;
      if(angle < 130) sig = 160-angle;
      
      lb = di - 0.3;
      ub = di + 0.3;
      if(!is_excluded(ex,al,don_b,acc) && (!al->isposres[don_b] 
                                           || !al->isposres[acc] 
                                           || !al->isposres[don]) &&
         !al->isflex[don_b] && !al->isflex[acc] &&
         !is_bound(bt,don_b,acc)){
        
        add_bound(bt,don_b,acc); 
        fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f %8d %9.4f %9.4f\n",
                al->id[don_b],al->id[acc],di,lb,ub,al->id[don],ang,sig);
        counter++;
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d bb-hbond angle constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}


      
/*=============================================================*/
int write_dihedrals(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                    t_dihed *dihed,
                     t_input *inp,t_excl *ex)

/* write 1-4 pairs to dist file */

{
  real di,lb,ub;
  int counter = 0;
  char logstr[STRLEN];
  int i,k;
  real vdwr;
  fprintf(fp,"[ 1-4 PAIRS ]\n");
  for(i=0;i<dihed->n;i++){
    if(!is_bound(bt,dihed->at1[i],dihed->at4[i]) &&
       !is_excluded(ex,al,dihed->at1[i],dihed->at4[i])
       && (!al->isposres[dihed->at1[i]] || !al->isposres[dihed->at4[i]])){
      
      di = DIST(al,dihed->at1[i],dihed->at4[i]);
      vdwr = get_vdw(al,dihed->at1[i],dihed->at4[i],TRUE);
      if(!dihed->flag[i]){ 
        lb = vdwr - inp->bump_14tol;
        ub = vdwr*(1+inp->non_restr_dihed_tol);
        if(di > ub && di < 10.) ub = di*1.2;
        
      }
      else{
        
        lb = di*(1-inp->restr_dihed_tol);
        ub = di*(1+inp->restr_dihed_tol);

      }

      if(lb < vdwr - inp->bump_14tol){
        lb = vdwr - inp->bump_14tol;
      }
      if(lb > di) lb = di*0.99;
      

      if(di > 10){
        lb = vdwr - inp->bump_14tol;
        ub = vdwr*(1+inp->non_restr_dihed_tol);
        di = lb + (ub-lb)/2.;
      }

      add_bound(bt,dihed->at1[i],dihed->at4[i]); 
      fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[dihed->at1[i]],al->id[dihed->at4[i]],di,lb,ub);  
      counter++;

      
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d 1-4 constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}
/*=============================================================*/
int write_planar(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                 t_idxgroups *pln,t_input *inp,t_excl *ex)

/* write planar groups to dist file */

{
  int i,j,k;
  int counter = 0;
  real di,lb,ub;
  char logstr[STRLEN];
  
  fprintf(fp,"[ PLANAR GROUPS ]\n");
  for(i=0;i<pln->n;i++){
/*     printf("i = %d natoms = %d\n",i,pln->natoms[i]); */
/*     print_group(al,pln->atoms[i],pln->natoms[i]); */
    
    for(k=0;k<pln->natoms[i];k++){
      for(j=k+1;j<pln->natoms[i];j++){
        if(!is_bound(bt,pln->atoms[i][j],pln->atoms[i][k]) &&
           !is_excluded(ex,al,pln->atoms[i][j],pln->atoms[i][k]) &&
           !al->isposres[pln->atoms[i][j]] && 
           !al->isposres[pln->atoms[i][k]]){

          add_bound(bt,pln->atoms[i][j],pln->atoms[i][k]);
          di = DIST(al,pln->atoms[i][j],pln->atoms[i][k]);
          lb = di-0.3;
          ub = di+0.3;
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[pln->atoms[i][j]],al->id[pln->atoms[i][k]],di,lb,ub);  
          counter++;
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d planar group constraints\n",counter);
  CNClog(log,logstr);
  
  return counter;
}

/*=============================================================*/
int write_hbonds(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt,
                 t_excl *ex) 

{
  
  int counter = 0;
  int i,k,j,l,id;
  int at1,at2;
  real di,lb,ub;
  real vdwr;
  int res1, res2;
  real tol;
  char logstr[STRLEN];
  
  fprintf(fp,"[ HBONDS ]\n");
  for(i=0;i<al->nhbonds;i++){
    
    if((al->hbonds[i]->type == ehBBBB ||
        al->hbonds[i]->type == ehBBSC ||
        al->hbonds[i]->type == ehSCBB ||
        al->hbonds[i]->type == ehSCSC) &&
       al->hbonds[i]->constr == TRUE)
    {
      /* donor - acceptor pair */
      
      at1 = al->hbonds[i]->don;
      at2 = al->hbonds[i]->acc;
      if(at1!=-1 && !is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) &&
         (!al->isposres[at1] || !al->isposres[at2]) &&
         (!al->isflex[at1] && !al->isflex[at2]))
      {
        di = al->hbonds[i]->d_a_dist;
        
        if(di <= 2.2) {
          lb = 1.7;
          ub = 2.4;
          if(lb > di) lb = di*0.99;
        }
        else {
          lb = 1.7;
          ub = 2.6;
          if(ub < di) ub = di*1.01;
        }
        
        fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],
                al->id[at2],di,lb,ub);
        add_bound(bt,at1,at2);
        counter++;
      }
      
      /* donor - acceptor_bond pair */

      at1 = al->hbonds[i]->don;
      at2 = al->hbonds[i]->acc_b;
      if(at1!=-1 && !is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) &&
         (!al->isposres[at1] || !al->isposres[at2]) &&
         (!al->isflex[at1] && !al->isflex[at2]))
      {
        di = al->hbonds[i]->d_ab_dist;
        lb = di*0.85;
        ub = di*1.15;
        fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],
                al->id[at2],di,lb,ub);
        add_bound(bt,at1,at2);
        counter++;
      }
      
      /* donor_bond - acceptor pair */
      at1 = al->hbonds[i]->don_b;
      at2 = al->hbonds[i]->acc;
      if(!is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) &&
         (!al->isposres[at1] || !al->isposres[at2]) &&
         (!al->isflex[at1] && !al->isflex[at2]))
      {
        di = al->hbonds[i]->db_a_dist;
        lb = di*0.85;
        ub = di*1.15;
        fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],
                al->id[at2],di,lb,ub);
        add_bound(bt,at1,at2);
        counter++;
      }


      
      /* donor_bond - acceptor_bond pair */
      at1 = al->hbonds[i]->don_b;
      at2 = al->hbonds[i]->acc_b;
      if(!is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) &&
         (!al->isposres[at1] || !al->isposres[at2]) &&
         (!al->isflex[at1] && !al->isflex[at2]))
      {
        di = al->hbonds[i]->db_ab_dist;
        lb = di*0.85;
        ub = di*1.15;
        fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],
                al->id[at2],di,lb,ub);
        add_bound(bt,at1,at2);
        counter++;
      }
    }
    
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d hydrogen bond constraints\n",counter);
  CNClog(log,logstr);
 
  return counter;
}


/*=============================================================*/

int write_hydrophobics(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                       t_contab *rct, t_resl *rl, t_idxgroups *phb,
                       t_input *inp, gmx_rng_t rng, real def, t_excl *ex,
                       bool bnoH)

/* write hydrophobic restrictions to dist file */

{
  int counter = 0;
  int i,k,j,l,m,n;
  int at;
  
  real di,lb,ub;
  real vdwr;
  int res1, res2;
  real tol;
  real fixdist;
  char logstr[STRLEN];
  
  if(bnoH) fixdist = 8.;
  else fixdist = 5.;

 
  fprintf(fp,"[ HYDROPHOBICS ]\n");
  
  for(i=0;i<phb->n;i++){
    if(phb->natoms[i] > 3) {
      for(k=0;k<phb->natoms[i]-1;k++){
        res1 = phb->atoms[i][k];
        for(l=k+1;l<phb->natoms[i];l++){
          res2 = phb->atoms[i][l];
          for(m=rl->j0[res1];m<=rl->j1[res1];m++){
            for(n=rl->j0[res2];n<=rl->j1[res2];n++){
              
              if(!is_bound(bt,m,n) && !is_excluded(ex,al,m,n)
                 && (!al->isposres[m] || !al->isposres[n]) 
                 && (!al->isflex[m] && !al->isflex[n])){


/*                 if(al->order[m] >= 2  && */
/*                    al->order[n] >= 2 ){ */
/*                   tol = 0.2; */
/*                 } */
/*                 else{ */
/*                   tol = 0.5-0.05*(al->order[m]+al->order[n]); */
/*                 } */
                di = DIST(al,m,n);
                
                lb = di*.8;
                ub = di*1.2;
                
/*                 lb = di - di*tol;  */
/*                 ub = di + di*tol;  */
                
                vdwr = get_vdw(al,m,n,FALSE);
                if(lb < vdwr - inp->bump_tol) lb = vdwr - inp->bump_tol;
                if (lb > di) lb = di*.99;
                
                if(gmx_rng_uniform_real(rng)*2. < def || di < fixdist){
                  add_bound(bt,m,n);
                  fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[m],al->id[n],di,lb,ub);
                  counter++;
                }
              }
            }
          }
        }
      }
    }
  }
  
/*     return counter; */
    





/*   for(i=0;i<rct->n;i++){ */
/*     for(k=0;k<rct->ncon[i];k++){ */
/*       if(rct->type[i][k] == ePHO){ */
/*         res1 = i; */
/*         res2 = rct->con[i][k]; */
/*         for(j=rl->j0[res1];j<=rl->j1[res1];j++){ */
/*           for(l=rl->j0[res2];l<=rl->j1[res2];l++){ */
/*             if(!is_bound(bt,j,l) && !is_excluded(ex,al,j,l) */
/*                && (!al->isposres[j] || !al->isposres[l])  */
/*                && (!al->isflex[j] && !al->isflex[l])){ */
/*               if(al->order[l] >= 2  && */
/*                  al->order[j] >= 2 ){ */
/*                 tol = 0.2; */
/*               } */
/*               else{ */
/*                 tol = 0.5-0.05*(al->order[j]+al->order[l]); */
/*               } */
/*               di = DIST(al,j,l); */
/*               lb = di - di*tol;  */
/*               ub = di + di*tol;  */

/*               vdwr = get_vdw(al,j,l,FALSE); */
/*               if(lb < vdwr - inp->bump_tol) lb = vdwr - inp->bump_tol; */
/*               if (lb > di) lb = di*.99; */

/*               if(gmx_rng_uniform_real(rng)*2. < def || di < 4.){ */
/*                add_bound(bt,j,l); */
/*                fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[j],al->id[l],di,lb,ub); */
/*                counter++; */
/*               } */
/*             } */
/*           } */
/*         } */
/*       } */
/*     } */
/*   } */
  sprintf(logstr,"tCNC__log_> Wrote %10d hydrophobic cluster constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}

/*=============================================================*/
int write_sidechain(FILE *fp, FILE *log, t_atomlist *al, t_resl *rl, 
                    t_boundtrack *bt, int min)
{
  int i,j,k,l,beg,end;
  real di,lb,ub;
  int counter = 0;
  char logstr[STRLEN];
  
  t_scdata *data = read_scdata();
  fprintf(fp,"[ SIDECHAIN ]\n");
  for (i=0;i<rl->nres;i++){
    beg = rl->j0[i];
    end = rl->j1[i];
    for(j=beg;j<end-1;j++){
      for(k=j+1;k<end;k++){
        for(l=0;l<data->n;l++){
          if(strcmp(data->resname[l],al->resname[j]) == 0 &&
             ((strcmp(al->name[j],data->name1[l])==0 &&
              strcmp(al->name[k],data->name2[l])==0) || 
             (strcmp(al->name[j],data->name2[l])==0 &&
              strcmp(al->name[k],data->name1[l])==0)) && 
             !is_bound(bt,j,k) && (!al->isposres[j] || !al->isposres[k]) 
             && (bt->n[j] < min || bt->n[k] < min)){
            di = dist_ij(al->x[j],al->x[k]);
            lb = data->lb[l];
            ub = data->ub[l];
            if(di < lb) lb = di*0.99;
            if(di > ub) ub = di*1.01;
            if (ub - lb < 0.6) {
/*               fprintf(stderr,"%s (%s-%s) %g %g\n",al->resname[j],al->name[j],al->name[k],lb,ub); */
              
              real tmp = 0.6 - (ub - lb);
              ub+=tmp;
              lb-=tmp;
/*               fprintf(stderr,"Strechting bound\n"); */
/*               fprintf(stderr,"new %s (%s-%s) %g %g\n",al->resname[j],al->name[j],al->name[k],lb,ub); */

            }
            

            add_bound(bt,j,k);
            fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[j],al->id[k],di,lb,ub);
            counter++;
          }
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d sidechain constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}


/*=============================================================*/
int write_close_pairs(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                      t_input *inp,int min, real fixdist, t_excl *ex)

/* write close pairs to dist file */

{
  int i,k;
  real di,lb,ub;
  real vdwr;
  char logstr[STRLEN];
  
  
  int counter = 0;
  fprintf(fp,"[ CLOSE PAIRS ]\n");
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      if(i > al->nb[i][k]){
        if(!is_bound(bt,i,al->nb[i][k]) && 
           !is_excluded(ex,al,i,al->nb[i][k]) && (!al->isposres[i] || !al->isposres[al->nb[i][k]])){
          di = DIST(al,i,al->nb[i][k]);
          if(di < fixdist && (bt->n[i] < min ||
                              bt->n[al->nb[i][k]] < min) ){
            
            
            vdwr = get_vdw(al,i,al->nb[i][k],FALSE);
            
 /*            if( fabs(al->resid[i]-al->resid[al->nb[i][k]]) < 5) {  */
              lb = di*.8;
 /*            } */
            
/*             else{ */
/*                 lb = di*.5; */
/*               } */
              
            if(lb < vdwr - inp->bump_tol)
              lb = vdwr - inp->bump_tol;
            
/*             if (di < lb) */
/*               lb = di*0.99; */

/*             lb = di * 0.5 - 5.;   */
            if(lb > di) lb=di*0.99;
/*             ub = di*1.2;   */
            
/*             if (fabs(al->resid[i]-al->resid[al->nb[i][k]]) < 5) */
              ub = di*1.2;
 /*            else { */
/*               ub = di*2; */
/*             } */
            
            add_bound(bt,i,al->nb[i][k]);
            fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f ; %d %d (%d, %d)\n",
                    al->id[i],al->id[al->nb[i][k]],di,lb,ub,
                    al->resid[i],al->resid[al->nb[i][k]],
                    al->ngrps[i],al->ngrps[al->nb[i][k]]);
            counter++;
          }
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d close pair constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}
/*=============================================================*/
int write_packing_constraints(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                              t_excl *ex, real pack_limit, t_input *inp)


/* write well packed pairs to dist file */

{
  int i,k;
  real di,lb,ub;
  real vdwr;
  char logstr[STRLEN];
  int counter = 0;
  fprintf(fp,"[ PACKING CONSTRAINTS ]\n");
  for(i=0;i<al->natoms;i++){
    if(al->bfac[i] < pack_limit) {
      for(k=0;k<al->nnb[i];k++){
        if(i > al->nb[i][k] && al->bfac[al->nb[i][k]] < pack_limit){
          if(!is_bound(bt,i,al->nb[i][k]) && 
             !is_excluded(ex,al,i,al->nb[i][k]) && (!al->isposres[i] || !al->isposres[al->nb[i][k]])){
            di = DIST(al,i,al->nb[i][k]);
            if(di < 6.) {
              vdwr = get_vdw(al,i,al->nb[i][k],FALSE);
              lb = di*.8;
              if(lb < vdwr - inp->bump_tol)
                lb = vdwr - inp->bump_tol;
              if(lb > di) lb=di*0.99;
              ub = di*1.2;  
              add_bound(bt,i,al->nb[i][k]);
              fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[i],al->id[al->nb[i][k]],di,lb,ub);
              counter++;
            }
            
          }
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d packing constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}
/*=============================================================*/
int write_forced(FILE *fp, FILE *log, t_atomlist *al, t_resl *rl, t_force *fo, 
                 t_input *inp, t_boundtrack *bt)

/* write forced constraints to dist file */

{
  int counter = 0;
  char logstr[STRLEN];
  int i,k,j;
  int res1,res2;
  real di;
  real lb,ub;
  real vdw;
  fprintf(fp,"[ FORCED ]\n");
  for(i=0;i<fo->n;i++){
    res1 = fo->id1[i] -1;
    res2 = fo->id2[i] -1;
    for(k=rl->j0[res1];k<=rl->j1[res1];k++){
      for(j=rl->j0[res2];j<=rl->j1[res2];j++){
        if(!is_bound(bt,j,k)){
          di = DIST(al,k,j);
          lb = di+di*fo->lb[i];
          ub = di+di*fo->ub[i];
          vdw = get_vdw(al,k,j,FALSE);
          if(lb < vdw-inp->bump_tol) lb = vdw-inp->bump_tol;
          add_bound(bt,j,k);
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[j],al->id[k],di,lb,ub);
          counter++;
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d forced constraints\n",counter);
  CNClog(log,logstr);
  
  return counter;
}
/*=============================================================*/
int write_network_restrictions(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                               t_idxgroups *resgroups,
                               t_resl *rl,int max, gmx_rng_t rng,
                               t_input *inp, real def, int maxn,
                               t_excl *ex, bool bnoH)

/* write network restrictions to dist file */

{
  int i,k,j,l,id;
  char logstr[STRLEN];
  
  real di,lb,ub;
  real vdwr;
  int at1, at2;
  int res1,res2;
  int counter = 0;
  int nconstr[al->natoms];
  real fixdist;
  if(bnoH) fixdist = 8.;
  else fixdist = 5.;
  bool in_same_group;
  
  
  for(i=0;i<al->natoms;i++){
    nconstr[i] = 0;
  }
  int r_array[resgroups->n];
  rand_array_int(rng,resgroups->n,r_array);

  fprintf(fp,"[ NETWORK ]\n");

  /* do close pairs in network first */
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      in_same_group = FALSE;
      for(j=0;j<al->ngrps[i];j++){
        for(l=0;l<al->ngrps[al->nb[i][k]];l++){
          if(al->grpnr[i][j] == al->grpnr[al->nb[i][k]][l]) 
            in_same_group = TRUE;
        }
      }
      if(in_same_group){
        at1 = i;
        at2 = al->nb[i][k];
        if(!is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) 
           && (!al->isposres[at1] || !al->isposres[at2]) 
           && (!al->isflex[at1] && !al->isflex[at2])){
          di = DIST(al,at1,at2);
          lb = di*inp->network_tol[0];
          vdwr = get_vdw(al,at1,at2,FALSE);
          if(lb < vdwr - inp->bump_tol){
            lb = vdwr - inp->bump_tol;
          }
          if(lb > di)
            lb = di*0.99; 
          ub = di *inp->network_tol[1];
          if(di < fixdist){ 
            nconstr[at1]++;
            nconstr[at2]++;
            add_bound(bt,at1,at2);
            fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],al->id[at2],di,lb,ub);
            counter++;
          }
        }
      }
    }
  }
  

  /* now do longer guys */

  for(i=0;i<resgroups->n;i++){
    id = r_array[i];
    for(k=0;k<resgroups->natoms[id];k++){
      
      int array[resgroups->natoms[id]];
      rand_array_int(rng,resgroups->natoms[id],array);
      for(j=k+1;j<resgroups->natoms[id];j++){  
        int fixed = 0;
        do{
          res1 = resgroups->atoms[id][array[j]];
          res2 = resgroups->atoms[id][k];
          if(res1 != res2){
            at1 = random_int_in_range(rng,rl->j0[res1],rl->j1[res1]);
            at2 = random_int_in_range(rng,rl->j0[res2],rl->j1[res2]);
/*             if((bt->n[at1] < max || bt->n[at2] < max) && !is_bound(bt,at1,at2)){ */
            if(nconstr[at1] < maxn && nconstr[at2] < maxn &&
               !is_bound(bt,at1,at2) && !is_excluded(ex,al,at1,at2) 
               && (!al->isposres[at1] || !al->isposres[at2]) 
               && (!al->isflex[at1] && !al->isflex[at2])){
              di = DIST(al,at1,at2);
              lb = di*inp->network_tol[0];
              vdwr = get_vdw(al,at1,at2,FALSE);
              if(lb < vdwr - inp->bump_tol){
                lb = vdwr - inp->bump_tol;
              }
              if(lb > di)
                lb = di*0.99; 
              ub = di *inp->network_tol[1];
              if(di < fixdist){
                nconstr[at1]++;
                nconstr[at2]++;
                add_bound(bt,at1,at2);
                fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],al->id[at2],di,lb,ub);
                counter++;
              }
              else{
                if(gmx_rng_uniform_real(rng) < def && bt->n[at1] < max && bt->n[at2] < max){
                  nconstr[at1]++;
                  nconstr[at2]++;
                  add_bound(bt,at1,at2);
                  fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[at1],al->id[at2],di,lb,ub);
                  counter++;
                }
              }
              
            }
         
          }fixed+=1;
          
        }while(fixed<6);
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d network constraints\n",counter);
  CNClog(log,logstr);
  return counter;
}

/*=============================================================*/
int long_range(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,t_resl *rl, 
               int max,gmx_rng_t rng,
                t_excl *ex, t_input *inp)

/* add long-range constraints to dist file */

{
  int counter = 0;
  char logstr[STRLEN];
  
  int i,j,k;
  real di,lb,ub;
  real vdw;
  int at1,at2;
  fprintf(fp,"[ LONG RANGE ]\n");

  int r_array[al->natoms];
  rand_array_int(rng,al->natoms,r_array);
  for(i=0;i<al->natoms;i++){
    int n = 0;
    int ntry = 0;
    do
    {
      int rdn = random_int_in_range(rng, 0,al->natoms);
      ntry++;
      if(r_array[i] != rdn && !is_bound(bt,r_array[i],rdn) && !is_excluded(ex,al,rdn,r_array[i])
         && (!al->isposres[r_array[i]] || !al->isposres[rdn]) && (bt->n[r_array[i]] < max && bt->n[rdn] < max) 
         && (!al->isflex[r_array[i]] && !al->isflex[rdn]) ){
        di = DIST(al,r_array[i],rdn);

/*         if(di < 10. || gmx_rng_uniform_real(rng) < 10./di) {     */
/*         if(di < 10. || gmx_rng_uniform_real(rng) < 100./sqr(di)) {   */
          
          lb = di * inp->lr_tol[0];
          vdw = get_vdw(al,r_array[i],rdn,FALSE);
          if (lb < vdw) lb = vdw;
          ub = di * inp->lr_tol[1];
/*           if(strcmp(al->chain[r_array[i]],al->chain[rdn]) == 0){ */
/*             int diff = abs(al->resid[r_array[i]]-al->resid[rdn]); */
/*             real max = 3.8*diff+10.; */
/*             if(ub > max) { */
/*               ub = max; */
/*             } */
/*           } */
          
          if(lb > di) lb = di*0.99;
          add_bound(bt,r_array[i],rdn);
          fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[r_array[i]],al->id[rdn],di,lb,ub);
          counter++;
          n++;
        }
/*       }     */
      
    }while(ntry < 70);
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d long range constraints\n",counter);
  CNClog(log,logstr);
  
  return counter;
}
/*=============================================================*/
int write_neighbor_res(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                       t_excl *ex, bool bIgn)
{
  int i,k;
  real di,lb,ub;
  real vdw;
  int counter = 0;
  char logstr[STRLEN];
  
  fprintf(fp,"[ CA-TRACE ]\n");
  
  for(i=0;i<al->natoms;i++){
    if(IsProtein(al,i)){
      if(strcmp(al->name[i]," CA ") == 0){
        for(k=i+1;k<al->natoms;k++){
          if(al->resid[i] == al->resid[k]-2 &&
             strcmp(al->chain[i],al->chain[k]) == 0 &&
             strcmp(al->name[k]," CA ") == 0 && IsProtein(al,k)){
            if(!is_bound(bt,i,k) && !is_excluded(ex,al,i,k) &&
               (!al->isposres[i] || !al->isposres[k])){
              di = DIST(al,i,k);
              lb = 5.;
              ub = 7.2;
              if(!bIgn) {
                if(lb > di) lb = di*0.99;
                if(ub < di) ub = di*1.01;
              }
              
              add_bound(bt,i,k);
              fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[i],al->id[k],di,lb,ub);
              counter++;
            }
          }
        }
      }
      else if(strcmp(al->name[i]," N  ") == 0){
        for(k=i+1;k<al->natoms;k++){
          if(al->resid[i] == al->resid[k]-2 &&
             strcmp(al->chain[i],al->chain[k]) == 0 &&
             strcmp(al->name[k]," N  ") == 0 && IsProtein(al,k)){
            if(!is_bound(bt,i,k) && !is_excluded(ex,al,i,k) &&
               (!al->isposres[i] || !al->isposres[k])){
              di = DIST(al,i,k);
              lb = 4.;
              ub = 7.5;
              if(!bIgn){
                if(lb > di) lb = di*0.99;
                if(ub < di) ub = di*1.01;
              }
              
              add_bound(bt,i,k);
              fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[i],al->id[k],di,lb,ub);
              counter++;
            }
          }
        }
      }
      
      else if(strcmp(al->name[i]," C  ") == 0){
        for(k=i+1;k<al->natoms;k++){
          if(al->resid[i] == al->resid[k]-2 &&
             strcmp(al->chain[i],al->chain[k]) == 0 &&
             strcmp(al->name[k]," C  ") == 0 && IsProtein(al,k)){
            if(!is_bound(bt,i,k) && !is_excluded(ex,al,i,k) &&
               (!al->isposres[i] || !al->isposres[k])){
              di = DIST(al,i,k);
              lb = 4.;
              ub = 7.5;
              if(!bIgn){
                if(lb > di) lb = di*0.99;
                if(ub < di) ub = di*1.01;
              }
              add_bound(bt,i,k);
              fprintf(fp,"%8d %7d %9.4f %9.4f %9.4f\n",al->id[i],al->id[k],di,lb,ub);
              counter++;
            }
          }
        }
      }
    }
  }
  sprintf(logstr,"tCNC__log_> Wrote %10d short range constraints\n",counter);
  CNClog(log,logstr);
  
  return counter;
}

/*=============================================================*/

void print_constr_per_atom(FILE *fp, t_boundtrack *bt, t_atomlist *al)
{
  int i;
  
  for(i=0;i<bt->natoms;i++){
    if(!al->isposres[i])
    fprintf(fp," %8d %10d ; %s%d%s\n",i+1,bt->n[i],al->name[i],al->resid[i],al->resname[i]);
  }
}

