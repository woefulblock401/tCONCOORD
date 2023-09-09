/*
 * $Id: calcpot.c,v 1.54.2.4 2008/02/29 07:02:54 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec.h"
#include "calcpot.h"
#include "nrnb.h"
#include "mdebin.h"
#include "mshift.h"
#include "smalloc.h"
#include "force.h"
#include "main.h"
#include "filenm.h"
#include "gmx_fatal.h"
#include "mdrun.h"
#include "ns.h"
#include "txtdump.h"

static void c_tabpot(real tabscale,   real VFtab[],
		     int  nri,        int  iinr[],
		     int  shift[],
		     int  jindex[],   int  jjnr[],
		     real pos[],      
		     real facel,      real charge[],
		     real pot[],      real shiftvec[])
{
  /* Local variables */
  const real nul = 0.000000;

  /* Table stuff */
  const real one = 1.000000;
  const real two = 2.000000;
  real r1,r1t,fijC,eps,eps2,Y,F,Fp,Geps,Heps2,VV,FF;
  int  n0,n1,nnn;

  /* General and coulomb stuff */
  int  ii,k,n,jnr,ii3,nj0,nj1,is3,j3,ggid;
  real fxJ,fyJ,fzJ,fxO,fyO,fzO;
  real ixO,iyO,izO,dxO,dyO,dzO;
  real txO,tyO,tzO,vcO,fsO,qO,rsqO,rinv1O,rinv2O;
  real qqO,qj;
  real jx,jy,jz,shX,shY,shZ,poti;

  /* Outer loop (over i particles) starts here */
  for(n=0; (n<nri); n++) {

    /* Unpack shift vector */
    is3               = 3*shift[n];
    shX               = shiftvec[is3];
    shY               = shiftvec[is3+1];
    shZ               = shiftvec[is3+2];

    /* Unpack I particle */
    ii                = iinr[n];
    ii3               = 3*ii;

    /* Charge of i particle(s) divided by 4 pi eps0 */
    qO                = facel*charge[ii];

    /* Bounds for the innerloop */
    nj0               = jindex[n];
    nj1               = jindex[n+1];

    /* Compute shifted i position */
    ixO               = shX + pos[ii3];
    iyO               = shY + pos[ii3+1];
    izO               = shZ + pos[ii3+2];
    poti              = nul;
    
    /* Inner loop (over j-particles) starts right here */
    for(k=nj0; (k<nj1); k++) {

      /* Unpack neighbourlist */
      jnr               = jjnr[k];
      j3                = 3*jnr;
      qj                = facel*charge[jnr];
      jx                = pos[j3];
      jy                = pos[j3+1];
      jz                = pos[j3+2];

      /* First one is for oxygen, with LJ */
      dxO               = ixO - jx;
      dyO               = iyO - jy;
      dzO               = izO - jz;
      rsqO              = dxO*dxO + dyO*dyO + dzO*dzO;

      /* Doing fast invsqrt */
      rinv1O            = invsqrt(rsqO);

      /* O block */
      r1                = one/rinv1O;
      r1t               = r1*tabscale;
      n0                = r1t;
      n1                = 12*n0;
      eps               = r1t-n0;
      eps2              = eps*eps;
      nnn               = n1;
      Y                 = VFtab[nnn];
      F                 = VFtab[nnn+1];
      Geps              = eps*VFtab[nnn+2];
      Heps2             = eps2*VFtab[nnn+3];
      Fp                = F+Geps+Heps2;
      VV                = Y+eps*Fp;

      pot[jnr]         += VV*qO;
      poti             += VV*qj;
      
    }
    pot[ii] += poti;
  }
}

static void low_calc_pot(FILE *log,int nl_type,t_forcerec *fr,
			 rvec x[],t_mdatoms *mdatoms,real pot[])
{
  t_nblist *nlist;
  
  nlist = &fr->nblists[0].nlist_sr[nl_type];
  
  c_tabpot(fr->nblists[0].tab.scale,fr->nblists[0].tab.tab,
	   nlist->nri,nlist->iinr,
	   nlist->shift,nlist->jindex,nlist->jjnr,
	   x[0],fr->epsfac,mdatoms->chargeA,pot,fr->shift_vec[0]);

  fprintf(log,"There were %d interactions\n",nlist->nrj);
}

void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
	      t_inputrec *inputrec,t_topology *top,rvec x[],t_forcerec *fr,
	      t_mdatoms *mdatoms,real pot[],matrix box,t_graph *graph)
{
  static bool        bFirst=TRUE;
  static t_nrnb      nrnb;
  static rvec        *f;
  real        lam=0,dum=0;
  rvec        box_size;
  int         i,m;

  /* Calc the force */
  fprintf(stderr,"Doing single force calculation...\n");

  if (bFirst) {
    snew(f,   nsb->natoms);
    
    bFirst = FALSE;
  }
  /* Reset long range forces if necessary */
  if (fr->bTwinRange) {
    clear_rvecs(nsb->natoms,fr->f_twin);
    clear_rvecs(SHIFTS,fr->fshift_twin);
  }
  if (inputrec->ePBC != epbcNONE)
    calc_shifts(box,fr->shift_vec);
  put_charge_groups_in_box(stdlog,0,top->blocks[ebCGS].nr,
			   box,&(top->blocks[ebCGS]),x,
			   fr->cg_cm);
  if (graph)
    mk_mshift(stdlog,graph,box,x);
  /* Do the actual neighbour searching and if twin range electrostatics
   * also do the calculation of long range forces and energies.
   */
  
  ns(logf,fr,x,f,box,grps,&(inputrec->opts),top,mdatoms,cr,
     &nrnb,nsb,0,lam,&dum,TRUE,FALSE);
  for(m=0; (m<DIM); m++)
    box_size[m] = box[m][m];
  for(i=0; (i<mdatoms->nr); i++)
    pot[i] = 0;
  if (debug) {
    pr_rvecs(debug,0,"x",x,mdatoms->nr);
    pr_rvecs(debug,0,"cgcm",fr->cg_cm,top->blocks[ebCGS].nr);
  }
  /* electrostatics from any atom to atoms without LJ */
  low_calc_pot(logf,eNL_QQ,fr,x,mdatoms,pot);
    /* electrostatics from any atom to atoms without charge */
  low_calc_pot(logf,eNL_VDW,fr,x,mdatoms,pot);
  /* electrostatics from any atom to atoms with LJ */
  low_calc_pot(logf,eNL_VDWQQ,fr,x,mdatoms,pot);
}

void init_calcpot(char *log,char *tpx,char *table,t_topology *top,
		  t_inputrec *inputrec,t_commrec *cr,
		  t_graph **graph,t_mdatoms **mdatoms,
		  t_nsborder *nsb,t_groups *grps,
		  t_forcerec **fr,real **pot,
		  matrix box,rvec **x)
{
  real     t,t0,lam,lam0;
  bool     bNEMD,bSA;
  char     *traj,*xtc_traj;
  t_state  *state;
  rvec     mutot;
  t_nrnb   nrnb;
  t_mdebin *mdebin;
  t_vcm    *vcm=NULL;
  int      fp_ene,m;
  rvec     box_size;
  tensor   force_vir,shake_vir;
  
  /* Initiate */
  cr->nnodes = 1; cr->nodeid    = 0; cr->left   = 0; cr->right  = 1;
  cr->nthreads = 1 ; cr->threadid = 0;
  open_log(log,cr);

  if (inputrec->efep) {
    fprintf(stderr,"WARNING: turning of free energy, will use lambda=0\n");
    inputrec->efep = 0;
  }

  init_nrnb(&nrnb);
  snew(state,1);
  init_single(stdlog,inputrec,tpx,top,state,mdatoms,nsb);
  init_md(cr,inputrec,&t,&t0,&lam,&lam0,
	  &nrnb,top,-1,NULL,&traj,&xtc_traj,&fp_ene,NULL,NULL,
	  &mdebin,grps,force_vir,
	  shake_vir,*mdatoms,mutot,&bNEMD,&bSA,&vcm,nsb);
  init_groups(stdlog,*mdatoms,&(inputrec->opts),grps);  

  if (inputrec->ePBC == epbcXYZ) {
    /* Calculate intramolecular shift vectors to make molecules whole again */
    *graph = mk_graph(&(top->idef),top->atoms.nr,FALSE,FALSE);
    mk_mshift(stdlog,*graph,state->box,state->x);
  } else {
    *graph = NULL;
  }

  /* Turn off twin range if appropriate */
  inputrec->rvdw  = inputrec->rcoulomb;
  inputrec->rlist = inputrec->rcoulomb;
  fprintf(stderr,"Using a coulomb cut-off of %g nm\n",inputrec->rcoulomb); 
  
  /* Turn off free energy computation */
  inputrec->efep = 0;

  /* Set vanderwaals to shift, to force tables */
  inputrec->vdwtype     = evdwSHIFT;
  inputrec->rvdw_switch = 0.0;
    
  /* Initiate forcerecord */
  *fr = mk_forcerec();
  init_forcerec(stdlog,*fr,inputrec,top,cr,*mdatoms,
		nsb,state->box,FALSE,table,table,TRUE);

  /* Remove periodicity */  
  for(m=0; (m<DIM); m++)
    box_size[m] = state->box[m][m];
  if (inputrec->ePBC != epbcNONE)
    do_pbc_first(stdlog,state->box,*fr,*graph,state->x);

  copy_mat(state->box,box);
  *x = state->x;
  state->x = NULL;
  done_state(state);
  sfree(state);

  snew(*pot,nsb->natoms);
}
