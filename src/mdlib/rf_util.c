/*
 * $Id: rf_util.c,v 1.5.2.4 2008/02/29 07:02:52 spoel Exp $
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

#include "sysstuff.h"
#include "typedefs.h"
#include "force.h"
#include "names.h"
#include "vec.h"
#include "physics.h"
#include "copyrite.h"
#include "pbc.h"
#include "network.h"

real RF_excl_correction(FILE *log,const t_nsborder *nsb,
			const t_commrec *cr,
			const t_forcerec *fr,t_graph *g,
			const t_mdatoms *mdatoms,const t_block *excl,
			rvec x[],rvec f[],rvec *fshift,const t_pbc *pbc,
			real lambda,real *dvdlambda)
{
  /* Calculate the reaction-field energy correction for this node:
   * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
   * and force correction for all excluded pairs, including self pairs.
   */
  static bool bFirst=TRUE;
  static bool bSumForces=FALSE;
  int    top,i,j,j1,j2,k,ki;
  double q2sumA,q2sumB,ener;
  const real *chargeA,*chargeB;
  real   ek,ec,L1,qiA,qiB,qqA,qqB,qqL,v;
  rvec   dx,df;
  atom_id *AA;
  ivec   dt;
  int    start = START(nsb);
  int    end   = start+HOMENR(nsb);
  bool   bFullPBC;

  if (fr->bTPI)
    /* For test particle insertion we only correct for the test particle */
    start = mdatoms->nr - 1;

  ek = fr->epsfac*fr->k_rf;
  ec = fr->epsfac*fr->c_rf;
  chargeA = mdatoms->chargeA;
  chargeB = mdatoms->chargeB;
  AA = excl->a;
  bFullPBC = (fr->ePBC == epbcFULL);
  ki = CENTRAL;

  q2sumA = 0;
  q2sumB = 0;
  ener = 0;
  if (fr->efep == efepNO) {
    for(i=start; i<end; i++) {
      qiA = chargeA[i];
      q2sumA += qiA*qiA;
      /* Do the exclusions */
      j1  = excl->index[i];
      j2  = excl->index[i+1];
      for(j=j1; j<j2; j++) {
	k = AA[j];
	if (k > i) {
	  bSumForces = bSumForces || (k >= end);
	  qqA = qiA*chargeA[k];
	  if (qqA != 0) {
	    if (g) {
	      rvec_sub(x[i],x[k],dx);
	      ivec_sub(SHIFT_IVEC(g,i),SHIFT_IVEC(g,k),dt);
	      ki=IVEC2IS(dt);
	    } else if (bFullPBC) {
	      ki = pbc_dx(pbc,x[i],x[k],dx);
	    } else
	      rvec_sub(x[i],x[k],dx);
	    ener += qqA*(ek*norm2(dx) - ec);
	    svmul(-2*qqA*ek,dx,df);
	    rvec_inc(f[i],df);
	    rvec_dec(f[k],df);
	    rvec_inc(fshift[ki],df);
	    rvec_dec(fshift[CENTRAL],df);
	  }
	}
      }
    }
    ener += -0.5*ec*q2sumA;
  } else {
    L1 = 1.0 - lambda;
    for(i=start; i<end; i++) {
      qiA = chargeA[i];
      qiB = chargeB[i];
      q2sumA += qiA*qiA;
      q2sumB += qiB*qiB;
      /* Do the exclusions */
      j1  = excl->index[i];
      j2  = excl->index[i+1];
      for(j=j1; j<j2; j++) {
	k = AA[j];
	if (k > i) {
	  qqA = qiA*chargeA[k];
	  qqB = qiB*chargeB[k];
	  if (qqA != 0 || qqB != 0) {
	    qqL = L1*qqA + lambda*qqB;
	    if (g) {
	      rvec_sub(x[i],x[k],dx);
	      ivec_sub(SHIFT_IVEC(g,i),SHIFT_IVEC(g,k),dt);
	      ki=IVEC2IS(dt);
	    } else if (bFullPBC) {
	      ki = pbc_dx(pbc,x[i],x[k],dx);
	    } else
	      rvec_sub(x[i],x[k],dx);
	    v = ek*norm2(dx) - ec;
	    ener += qqL*v;
	    svmul(-2*qqL*ek,dx,df);
	    rvec_inc(f[i],df);
	    rvec_dec(f[k],df);
	    rvec_inc(fshift[ki],df);
	    rvec_dec(fshift[CENTRAL],df);
	    *dvdlambda += (qqB - qqA)*v;
	  }
	}
      }
    }
    ener += -0.5*ec*(L1*q2sumA + lambda*q2sumB);
    *dvdlambda += -0.5*ec*(q2sumB - q2sumA);
  }
  if (bFirst) {
    if (PAR(cr))
      gmx_sumi(1,&bSumForces,cr);
    else 
      bSumForces = FALSE;
    bFirst = FALSE;
  }
  if (bSumForces) {
    /* This is necessary if molecules are split over processors. Should
       be optimized! */
    gmx_sum(nsb->natoms*DIM,fr->f_el_recip[0],cr);
  }
    

  if(debug)
    fprintf(debug,"RF exclusion energy: %g\n",ener);
  
  return ener;
}

void calc_rffac(FILE *log,int eel,real eps_r,real eps_rf,real Rc,real Temp,
		real zsq,matrix box,
		real *kappa,real *krf,real *crf)
{
  /* Compute constants for Generalized reaction field */
  static bool bFirst=TRUE;
  real   k1,k2,I,vol,rmin;
  
  if (EEL_RF(eel)) {
    vol     = det(box);
    I       = zsq/vol;
    if (eel == eelGRF) {
      /* Consistency check */
      if (Temp <= 0.0)
	gmx_fatal(FARGS,"Temperature is %f while using"
		    " Generalized Reaction Field\n",Temp);
      /* Ionic strength (only needed for eelGRF */
      *kappa  = sqrt(2*I/(EPSILON0*eps_rf*BOLTZ*Temp));
    }
    else
      *kappa = 0;

    /* eps == 0 signals infinite dielectric */
    if (eps_rf == 0) {
      *krf = 1/(2*Rc*Rc*Rc);
    } else {
      k1   = 1 + *kappa*Rc;
      k2   = eps_rf*sqr((real)(*kappa*Rc));
      
      *krf = ((eps_rf - eps_r)*k1 + 0.5*k2)/((2*eps_rf + eps_r)*k1 + k2)/(Rc*Rc*Rc);
    }
    *crf   = 1/Rc + *krf*Rc*Rc;
    rmin   = pow(*krf*2.0,-1.0/3.0);
    
    if (bFirst) {
      if (eel == eelGRF)
	please_cite(log,"Tironi95a");
      fprintf(log,"%s:\n"
	      "epsRF = %10g, I   = %10g, volume = %10g, kappa  = %10g\n"
	      "rc    = %10g, krf = %10g, crf    = %10g, epsfac = %10g\n",
	      eel_names[eel],eps_rf,I,vol,*kappa,Rc,*krf,*crf,
	      ONE_4PI_EPS0/eps_r);
      fprintf(log,
	      "The electrostatics potential has its minimum at rc = %g\n",
	      rmin);
      
      bFirst=FALSE;
    }
  }
}
