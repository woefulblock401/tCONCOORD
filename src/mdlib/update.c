/*
 * $Id: update.c,v 1.77.2.2 2008/02/29 07:02:52 spoel Exp $
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


#include <stdio.h>
#include <math.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "nrnb.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "main.h"
#include "confio.h"
#include "update.h"
#include "gmx_random.h"
#include "futil.h"
#include "mshift.h"
#include "tgroup.h"
#include "force.h"
#include "names.h"
#include "txtdump.h"
#include "mdrun.h"
#include "copyrite.h"
#include "constr.h"
#include "edsam.h"
#include "pull.h"

typedef struct {
  real gdt;
  real eph;
  real emh;
  real em;
  real b;
  real c;
  real d;
  real vvcorr;
} t_sdconst;





static t_sdconst *sdc;

static void do_update_md(int start,int homenr,double dt,
                         t_grp_tcstat *tcstat,t_grp_acc *gstat,real nh_xi[],
                         rvec accel[],ivec nFreeze[],real invmass[],
                         unsigned short ptype[],unsigned short cFREEZE[],
                         unsigned short cACC[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],rvec vold[],
                         rvec f[],matrix M,bool bExtended)
{
  double imass,w_dt;
  int    gf,ga,gt;
  rvec   vrel;
  real   vn,vv,va,vb,vnrel;
  real   lg,xi,u;
  int    n,d;

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parrinello-Rahman coupling
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */
    for(n=start; n<start+homenr; n++) {
      imass = invmass[n];
      gf   = cFREEZE[n];
      ga   = cACC[n];
      gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      xi   = nh_xi[gt];

      rvec_sub(v[n],gstat[ga].u,vrel);

      for(d=0; d<DIM; d++) {
        vold[n][d]     = v[n][d];

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
          vnrel = lg*(vrel[d] + dt*(imass*f[n][d] - 0.5*xi*vrel[d]
				    - iprod(M[d],vrel)))/(1 + 0.5*xi*dt);  
          /* do not scale the mean velocities u */
          vn             = gstat[ga].u[d] + accel[ga][d]*dt + vnrel; 
          v[n][d]        = vn;
          xprime[n][d]   = x[n][d]+vn*dt;
        } else
          xprime[n][d]   = x[n][d];
      }
    }

  } else {
    /* Classic version of update, used with berendsen coupling */
    for(n=start; n<start+homenr; n++) {
      w_dt = invmass[n]*dt;
      gf   = cFREEZE[n];
      ga   = cACC[n];
      gt   = cTC[n];
      lg   = tcstat[gt].lambda;

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];
        vold[n][d]     = vn;

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
          vv             = lg*(vn + f[n][d]*w_dt);

          /* do not scale the mean velocities u */
          u              = gstat[ga].u[d];
          va             = vv + accel[ga][d]*dt;
          vb             = va + (1.0-lg)*u;
          v[n][d]        = vb;
          xprime[n][d]   = x[n][d]+vb*dt;
        } else {
          v[n][d]        = 0.0;
          xprime[n][d]   = x[n][d];
        }
      }
    }
  }
}

static void do_update_visc(int start,int homenr,double dt,
			   t_grp_tcstat *tcstat,real invmass[],real nh_xi[],
                           unsigned short ptype[],unsigned short cTC[],
                           rvec x[],rvec xprime[],rvec v[],rvec vold[],
                           rvec f[],matrix M,matrix box,real
                           cos_accel,real vcos,bool bExtended)
{
  double imass,w_dt;
  int    gt;
  real   vn,vc;
  real   lg,xi,vv;
  real   fac,cosz;
  rvec   vrel;
  int    n,d;

  fac = 2*M_PI/(box[ZZ][ZZ]);

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parrinello-Rahman coupling
     */
    for(n=start; n<start+homenr; n++) {
      imass = invmass[n];
      gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      cosz = cos(fac*x[n][ZZ]);

      copy_rvec(v[n],vold[n]);
      copy_rvec(v[n],vrel);

      vc            = cosz*vcos;
      vrel[XX]     -= vc;
      xi           = nh_xi[gt];

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
          vn    = lg*(vrel[d] + dt*(imass*f[n][d] - 0.5*xi*vrel[d]
				    - iprod(M[d],vrel)))/(1 + 0.5*xi*dt);
          if(d == XX)
            vn += vc + dt*cosz*cos_accel;

          v[n][d]        = vn;
          xprime[n][d]   = x[n][d]+vn*dt;
        } else
          xprime[n][d]   = x[n][d];
      }
    }

  } else {
    /* Classic version of update, used with berendsen coupling */
    for(n=start; n<start+homenr; n++) {
      w_dt = invmass[n]*dt;
      gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      cosz = cos(fac*x[n][ZZ]);

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];
        vold[n][d]     = vn;

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
          if(d == XX) {
            vc           = cosz*vcos;
            /* Do not scale the cosine velocity profile */
            vv           = vc + lg*(vn - vc + f[n][d]*w_dt);
            /* Add the cosine accelaration profile */
            vv          += dt*cosz*cos_accel;
          } else
            vv           = lg*(vn + f[n][d]*w_dt);

          v[n][d]        = vv;
          xprime[n][d]   = x[n][d]+vv*dt;
        } else {
          v[n][d]        = 0.0;
          xprime[n][d]   = x[n][d];
        }
      }
    }
  }
}


void init_sd_consts(int ngtc,real tau_t[],real dt)
{
  int  n;
  real y;

  snew(sdc,ngtc);

  for(n=0; n<ngtc; n++) {
    sdc[n].gdt = dt/tau_t[n];
    sdc[n].eph = exp(sdc[n].gdt/2);
    sdc[n].emh = exp(-sdc[n].gdt/2);
    sdc[n].em  = exp(-sdc[n].gdt);
    if(sdc[n].gdt >= 0.25) {
      sdc[n].b = sdc[n].gdt*(sqr(sdc[n].eph)-1) - 4*sqr(sdc[n].eph-1);
      sdc[n].c = sdc[n].gdt - 3 + 4*sdc[n].emh - sdc[n].em;
      sdc[n].d = 2 - sdc[n].eph - sdc[n].emh;
    } else {
      y = sdc[n].gdt/2;
      /* Seventh order expansions for small y */
      sdc[n].b = y*y*y*y*(1/3.0+y*(1/3.0+y*(17/90.0+y*7/9.0)));
      sdc[n].c = y*y*y*(2/3.0+y*(-1/2.0+y*(7/30.0+y*(-1/12.0+y*31/1260.0))));
      sdc[n].d = y*y*(-1+y*y*(-1/12.0-y*y/360.0));
    }
    /* The missing velocity correlation over one MD step */
    sdc[n].vvcorr = 0.5*(1 - sdc[n].em);
    if(debug)
      fprintf(debug,"SD const tc-grp %d: b %g  c %g  d %g  vvcorr %g\n",
              n,sdc[n].b,sdc[n].c,sdc[n].d,sdc[n].vvcorr);
  }
}

static void do_update_sd(int start,int homenr,
                         rvec accel[],ivec nFreeze[],
                         real invmass[],unsigned short ptype[],
                         unsigned short cFREEZE[],unsigned short cACC[],
                         unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],rvec vold[],rvec f[],
			 rvec sd_X[],
                         int ngtc,real tau_t[],real ref_t[],
                         gmx_rng_t gaussrand, bool bFirstHalf)
{
  typedef struct {
    real V;
    real X;
    real Yv;
    real Yx;
  } t_sd_sigmas;

  static bool bFirst = TRUE;
  static t_sd_sigmas *sig=NULL;
  /* The random part of the velocity update, generated in the first
   * half of the update, needs to be remembered for the second half.
   */
  static rvec *sd_V;
  real   kT;
  int    gf,ga,gt;
  real   vn=0,Vmh,Xmh;
  real   ism;
  int    n,d;
  unsigned long  jran;

  if(sig == NULL) {
    snew(sig,ngtc);
    snew(sd_V,homenr);
  }

  if(bFirstHalf) {
    for(n=0; n<ngtc; n++) {
      kT = BOLTZ*ref_t[n];
      /* The mass is encounted for later, since this differs per atom */
      sig[n].V  = sqrt(kT*(1-sdc[n].em));
      sig[n].X  = sqrt(kT*sqr(tau_t[n])*sdc[n].c);
      sig[n].Yv = sqrt(kT*sdc[n].b/sdc[n].c);
      sig[n].Yx = sqrt(kT*sqr(tau_t[n])*sdc[n].b/(1-sdc[n].em));
    }
  }

  for(n=start; n<start+homenr; n++) {
    ism = sqrt(invmass[n]);
    gf  = cFREEZE[n];
    ga  = cACC[n];
    gt  = cTC[n];

    for(d=0; d<DIM; d++) {
      if(bFirstHalf) {
        vn             = v[n][d];
        vold[n][d]     = vn;
      }
      if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
        if(bFirstHalf) {

          if(bFirst)
            sd_X[n][d] = ism*sig[gt].X*gmx_rng_gaussian_table(gaussrand);

          Vmh = sd_X[n][d]*sdc[gt].d/(tau_t[gt]*sdc[gt].c) 
                + ism*sig[gt].Yv*gmx_rng_gaussian_table(gaussrand);
          sd_V[n-start][d] = ism*sig[gt].V*gmx_rng_gaussian_table(gaussrand);

          v[n][d] = vn*sdc[gt].em 
                    + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
                    + sd_V[n-start][d] - sdc[gt].em*Vmh;

          xprime[n][d] = x[n][d] + v[n][d]*tau_t[gt]*(sdc[gt].eph - sdc[gt].emh); 

        } else {

          /* Correct the velocities for the constraints */
          v[n][d] = 
          (xprime[n][d] - x[n][d])/(tau_t[gt]*(sdc[gt].eph - sdc[gt].emh));  

          Xmh = sd_V[n-start][d]*tau_t[gt]*sdc[gt].d/(sdc[gt].em-1) 
                + ism*sig[gt].Yx*gmx_rng_gaussian_table(gaussrand);
          sd_X[n][d] = ism*sig[gt].X*gmx_rng_gaussian_table(gaussrand);

          xprime[n][d] += sd_X[n][d] - Xmh;

        }
      } else {
        if(bFirstHalf) {
          v[n][d]        = 0.0;
          xprime[n][d]   = x[n][d];
        }
      }
    }
  }

  bFirst = FALSE;
}

static void do_update_bd(int start,int homenr,double dt,
                         ivec nFreeze[],
                         real invmass[],unsigned short ptype[],
                         unsigned short cFREEZE[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],rvec vold[],
                         rvec f[],real friction_coefficient,
                         int ngtc,real tau_t[],real ref_t[],
                         gmx_rng_t gaussrand)
{
  int    gf,gt;
  real   vn;
  static real *rf=NULL;
  real   invfr=0;
  int    n,d;
  unsigned long  jran;

  if (rf == NULL)
    snew(rf,ngtc);

  if (friction_coefficient != 0) {
    invfr = 1.0/friction_coefficient;
    for(n=0; n<ngtc; n++)
      rf[n] = sqrt(2.0*BOLTZ*ref_t[n]/(friction_coefficient*dt));
  } else
    for(n=0; n<ngtc; n++)
      rf[n] = sqrt(2.0*BOLTZ*ref_t[n]);

  for(n=start; (n<start+homenr); n++) {
    gf = cFREEZE[n];
    gt = cTC[n];
    for(d=0; (d<DIM); d++) {
      vold[n][d]     = v[n][d];
      if((ptype[n]!=eptVSite) && (ptype[n]!=eptShell) && !nFreeze[gf][d]) {
        if (friction_coefficient != 0)
          vn = invfr*f[n][d] + rf[gt]*gmx_rng_gaussian_table(gaussrand);
        else
          /* NOTE: invmass = 1/(mass*friction_constant*dt) */
          vn = invmass[n]*f[n][d]*dt 
	    + sqrt(invmass[n])*rf[gt]*gmx_rng_gaussian_table(gaussrand);

        v[n][d]      = vn;
        xprime[n][d] = x[n][d]+vn*dt;
      } else {
        v[n][d]      = 0.0;
        xprime[n][d] = x[n][d];
      }
    }
  }
}

static void dump_it_all(FILE *fp,char *title,
                        int natoms,rvec x[],rvec xp[],rvec v[],
                        rvec vold[],rvec f[])
{
#ifdef DEBUG
  fprintf(fp,"%s\n",title);
  pr_rvecs(fp,0,"x",x,natoms);
  pr_rvecs(fp,0,"xp",xp,natoms);
  pr_rvecs(fp,0,"v",v,natoms);
  pr_rvecs(fp,0,"vold",vold,natoms);
  pr_rvecs(fp,0,"f",f,natoms);
#endif
}

void calc_ke_part(int start,int homenr,rvec v[],
                  t_grpopts *opts,t_mdatoms *md,t_groups *grps,
                  t_nrnb *nrnb,real lambda)
{
  int          g,d,n,ga,gt;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_grp_acc    *grpstat=grps->grpstat;
  real         dekindl;

  /* group velocities are calculated in update_grps and
   * accumulated in acumulate_groups.
   * Now the partial global and groups ekin.
   */
  for(g=0; (g<opts->ngtc); g++) {
    copy_mat(grps->tcstat[g].ekinh,grps->tcstat[g].ekinh_old);
    clear_mat(grps->tcstat[g].ekinh);
  }
  grps->dekindl_old = grps->dekindl;

  dekindl = 0;
  for(n=start; (n<start+homenr); n++) {
    ga   = md->cACC[n];
    gt   = md->cTC[n];
    hm   = 0.5*md->massT[n];

    for(d=0; (d<DIM); d++) {
      v_corrt[d] = v[n][d] - grpstat[ga].u[d];
    }
    for(d=0; (d<DIM); d++) {
      tcstat[gt].ekinh[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekinh[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekinh[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if (md->nPerturbed && md->bPerturbed[n])
      dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
  }
  grps->dekindl = dekindl;

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

void calc_ke_part_visc(int start,int homenr,
                       matrix box,rvec x[],rvec v[],
                       t_grpopts *opts,t_mdatoms *md,t_groups *grps,
                       t_nrnb *nrnb,real lambda)
{
  int          g,d,n,gt;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_cos_acc    *cosacc=&(grps->cosacc);
  real         dekindl;
  real         fac,cosz;
  double       mvcos;

  for(g=0; g<opts->ngtc; g++) {
    copy_mat(grps->tcstat[g].ekinh,grps->tcstat[g].ekinh_old);
    clear_mat(grps->tcstat[g].ekinh);
  }
  grps->dekindl_old = grps->dekindl;

  fac = 2*M_PI/box[ZZ][ZZ];
  mvcos = 0;
  dekindl = 0;
  for(n=start; n<start+homenr; n++) {
    gt   = md->cTC[n];
    hm   = 0.5*md->massT[n];

    /* Note that the times of x and v differ by half a step */
    cosz         = cos(fac*x[n][ZZ]);
    /* Calculate the amplitude of the new velocity profile */
    mvcos       += 2*cosz*md->massT[n]*v[n][XX];

    copy_rvec(v[n],v_corrt);
    /* Subtract the profile for the kinetic energy */
    v_corrt[XX] -= cosz*cosacc->vcos;
    for(d=0; (d<DIM); d++) {
      tcstat[gt].ekinh[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekinh[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekinh[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if(md->nPerturbed && md->bPerturbed[n])
      dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
  }
  grps->dekindl = dekindl;
  cosacc->mvcos = mvcos;

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

/* Static variables for the deform algorithm */
static int step_store;
static matrix box_store;

static void deform_store(matrix box,const t_inputrec *ir,
			 int step,bool bFirstStep)
{
  if (bFirstStep ||
      ir->init_step+step == ir->nsteps ||
      (do_per_step(step,ir->nstxout)
       && (do_per_step(step,ir->nstvout) || ir->eI == eiBD))) {
    /* Store the structure, so we can avoid rounding problems
     * during slow deformation.
     * We need to store the structure again at steps from which
     * exact restarts can be performed.
     */
    step_store = step;
    copy_mat(box,box_store);
  }
}

static void deform(int start,int homenr,rvec x[],matrix box,
		   const t_inputrec *ir,int step)
{
  matrix new,invbox,mu;
  real   elapsed_time;
  int    i,j;  

  elapsed_time = (step + 1 - step_store)*ir->delta_t;
  copy_mat(box,new);
  for(i=0; i<DIM; i++)
    for(j=0; j<DIM; j++)
      if (ir->deform[i][j] != 0) {
	new[i][j] = box_store[i][j] + elapsed_time*ir->deform[i][j];
      }
  /* We correct the off-diagonal elements,
   * which can grow indefinitely during shearing,
   * so the shifts do not get messed up.
   */
  for(i=1; i<DIM; i++) {
    for(j=i-1; j>=0; j--) {
      while (new[i][j] - box[i][j] > 0.5*new[j][j])
	rvec_dec(new[i],new[j]);
      while (new[i][j] - box[i][j] < -0.5*new[j][j])
	rvec_inc(new[i],new[j]);
    }
  }
  m_inv_lowerleft0(box,invbox);
  copy_mat(new,box);
  mmul_lowerleft0(box,invbox,mu);
  
  for(i=start; i<start+homenr; i++) {
    x[i][XX] = mu[XX][XX]*x[i][XX]+mu[YY][XX]*x[i][YY]+mu[ZZ][XX]*x[i][ZZ];
    x[i][YY] = mu[YY][YY]*x[i][YY]+mu[ZZ][YY]*x[i][ZZ];
    x[i][ZZ] = mu[ZZ][ZZ]*x[i][ZZ];
  }
}


void update(int          natoms,  /* number of atoms in simulation */
            int          start,
            int          homenr,  /* number of home particles 	*/
            int          step,
            real         *dvdlambda,    /* FEP stuff */
            t_inputrec   *inputrec,      /* input record and box stuff	*/
            t_mdatoms    *md,
	    t_state      *state,
            t_graph      *graph,  
            rvec         force[],   /* forces on home particles 	*/
            rvec         vold[],  /* Old velocities		   */
            t_topology   *top,
            t_groups     *grps,
            tensor       vir_part,
            t_commrec    *cr,
            t_nrnb       *nrnb,
            t_edsamyn    *edyn,
            t_pull       *pulldata,
            bool         bNEMD,
	    bool         bDoUpdate,
	    bool         bFirstStep,
	    rvec         *shakefirst_x,
	    tensor       pres)
{
  static bool      bFirst=TRUE;
  static rvec      *xprime=NULL;
  static int       ngtc,ngacc;
  static bool      bHaveConstr,bExtended;
  double           dt;
  real             dt_1;
  int              i,n,m,g;
  matrix           M;
  tensor           vir_con;
  static gmx_rng_t sd_gaussrand=NULL;
  
  if(bFirst) {
    bHaveConstr = init_constraints(stdlog,top,inputrec,md,start,homenr,
                                   inputrec->eI!=eiSteep,cr);
    bHaveConstr = bHaveConstr || pulldata->bPull;
    bExtended   = (inputrec->etc==etcNOSEHOOVER) || (inputrec->epc==epcPARRINELLORAHMAN);

    /* if edyn->bEdsam == FALSE --- this has no effect */
    bHaveConstr = bHaveConstr || ed_constraints(edyn);  
    do_first_edsam(stdlog,top,md,start,homenr,cr,state->x,state->box,edyn,bHaveConstr);

    /* Initiate random number generator for stochastic and brownian dynamic integrators */
    if(inputrec->eI==eiSD || inputrec->eI==eiBD) 
      sd_gaussrand = gmx_rng_init(inputrec->ld_seed);
    
    /* Allocate memory for xold and for xprime. */
    snew(xprime,natoms);
    /* Copy the pointer to the external acceleration in the opts */
    ngacc=inputrec->opts.ngacc;    
    ngtc=inputrec->opts.ngtc;    
    /* Set Berendsen tcoupl lambda's to 1, 
     * so runs without Berendsen coupling are not affected.
     */
    for(i=0; i<ngtc; i++)
      grps->tcstat[i].lambda = 1.0;
  
    /* done with initializing */
    bFirst=FALSE;
  }

  dt   = inputrec->delta_t;
  dt_1 = 1.0/dt;
  if (bDoUpdate) {
    clear_mat(M);

    if (inputrec->etc==etcBERENDSEN)
      berendsen_tcoupl(&(inputrec->opts),grps,inputrec->delta_t);
    if (!bFirstStep) {
      if (inputrec->etc==etcNOSEHOOVER)
	nosehoover_tcoupl(&(inputrec->opts),grps,inputrec->delta_t,
			  state->nosehoover_xi);
      if (inputrec->epc == epcBERENDSEN)
	berendsen_pcoupl(inputrec,step,pres,state->box,state->pcoupl_mu);
    }
    if (inputrec->epc == epcPARRINELLORAHMAN)
      parrinellorahman_pcoupl(inputrec,step,pres,
			      state->box,state->boxv,M,bFirstStep);

    /* Now do the actual update of velocities and positions */
    where();
    dump_it_all(stdlog,"Before update",
		natoms,state->x,xprime,state->v,vold,force);
    if (inputrec->eI == eiMD) {
      if (grps->cosacc.cos_accel == 0)
        /* use normal version of update */
        do_update_md(start,homenr,dt,
		     grps->tcstat,grps->grpstat,state->nosehoover_xi,
                     inputrec->opts.acc,inputrec->opts.nFreeze,md->invmass,md->ptype,
                     md->cFREEZE,md->cACC,md->cTC,
		     state->x,xprime,state->v,vold,force,M,
                     bExtended);
      else
        do_update_visc(start,homenr,dt,
		       grps->tcstat,md->invmass,state->nosehoover_xi,
                       md->ptype,md->cTC,state->x,xprime,state->v,vold,force,M,
                       state->box,grps->cosacc.cos_accel,grps->cosacc.vcos,bExtended);
    } else if (inputrec->eI == eiSD) {
      /* The SD update is done in 2 parts, because an extra constraint step
       * is needed 
       */
      do_update_sd(start,homenr,
                   inputrec->opts.acc,inputrec->opts.nFreeze,
                   md->invmass,md->ptype,
                   md->cFREEZE,md->cACC,md->cTC,
                   state->x,xprime,state->v,vold,force,state->sd_X,
                   inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
                   sd_gaussrand,TRUE);
    } else if (inputrec->eI == eiBD) {
      do_update_bd(start,homenr,dt,
                   inputrec->opts.nFreeze,md->invmass,md->ptype,
                   md->cFREEZE,md->cTC,
                   state->x,xprime,state->v,vold,force,
                   inputrec->bd_fric,
                   inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
                   sd_gaussrand);
    } else {
      gmx_fatal(FARGS,"Don't know how to update coordinates");
    }
    where();
    inc_nrnb(nrnb, bExtended ? eNR_EXTUPDATE : eNR_UPDATE, homenr);
    dump_it_all(stdlog,"After update",
		natoms,state->x,xprime,state->v,vold,force);
  } else {
    /* If we're not updating we're doing shakefirst!
     */
    for(n=start; (n<start+homenr); n++) {
      copy_rvec(shakefirst_x[n],xprime[n]);
    }
  }

  /* 
   *  Steps (7C, 8C)
   *  APPLY CONSTRAINTS:
   *  BLOCK SHAKE 
   */

  /* When doing PR pressure coupling we have to constrain the
   * bonds in each iteration. If we are only using Nose-Hoover tcoupling
   * it is enough to do this once though, since the relative velocities 
   * after this will be normal to the bond vector
   */
  prepare_edsam(step,start,homenr,cr,xprime,edyn);
  
  if (bHaveConstr) {
    /* Constrain the coordinates xprime */
    constrain(stdlog,top,inputrec,step,md,start,homenr,state->x,xprime,NULL,
              state->box,state->lambda,dvdlambda,&vir_con,nrnb,TRUE);
    if (inputrec->eI == eiSD) {
      /* A correction factor eph is needed for the SD constraint force */
      /* Here we can, unfortunately, not have proper corrections
       * for different friction constants, so we use the first one.
       */
      for(i=0; i<DIM; i++)
	for(m=0; m<DIM; m++)
	  vir_part[i][m] += sdc[0].eph*vir_con[i][m];
    } else {
      m_add(vir_part,vir_con,vir_part);
    }
    where();

    dump_it_all(stdlog,"After Shake",
		natoms,state->x,xprime,state->v,vold,force);

    /* Apply Essential Dynamics constraints when required.
     * Note that the reported ED and pull forces can be off
     * with SD and BD with high friction, as the forces
     * are no longer proportional to positional differences.
     */

    /* Apply Essential Dynamics constraints when required. */
    if (edyn->bEdsam)
      do_edsam(stdlog,top,inputrec,step,md,start,homenr,cr,xprime,state->x,
               force,state->box,edyn,bDoUpdate);

    /* apply pull constraints when required. Act on xprime, the SHAKED
       coordinates. Don't do anything to f */
    if (pulldata->bPull && pulldata->runtype == eConstraint)
      pull(pulldata,xprime,force,vir_part,state->box,
	   top,dt,step,inputrec->init_t+step*dt,
	   md,start,homenr,cr);

    where();      

    if (bDoUpdate) {
      /* Note that with SD with ED or pull constraints the velocity
       * corrections for these constraints is missing (this only
       * involves a few degress of freedom though).
       */
      if (inputrec->eI != eiSD) {
        /* The constraint virial and the velocities are incorrect for BD */
        for(n=start; n<start+homenr; n++) {
          for(i=0; i<DIM; i++) {
            state->v[n][i] = (xprime[n][i] - state->x[n][i])*dt_1;
          }
        }
	inc_nrnb(nrnb,eNR_CONSTR_V,homenr);
        where();
      }
      dump_it_all(stdlog,"After Shake-V",
		  natoms,state->x,xprime,state->v,vold,force);
      where();

      if (debug)
	pr_rvecs(debug,0,"constraint virial",vir_part,DIM);
    }
  }  else if (edyn->bEdsam) {     
      /* no constraints but still edsam - yes, that can happen */
    do_edsam(stdlog,top,inputrec,step,md,start,homenr,cr,xprime,state->x,
 	     force,state->box,edyn,bDoUpdate);
  };  
  
  where();
  if (bDoUpdate) {
    if (inputrec->eI == eiSD) {
      /* The second part of the SD integration */
      do_update_sd(start,homenr,
		   inputrec->opts.acc,inputrec->opts.nFreeze,
		   md->invmass,md->ptype,
		   md->cFREEZE,md->cACC,md->cTC,
		   state->x,xprime,state->v,vold,force,state->sd_X,
		   inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
		   sd_gaussrand,FALSE);
      inc_nrnb(nrnb, eNR_UPDATE, homenr);
      
      if (bHaveConstr) {
	/* Constrain the coordinates xprime */
	constrain(stdlog,top,inputrec,step,md,start,homenr,state->x,xprime,NULL,
		  state->box,state->lambda,dvdlambda,NULL,nrnb,TRUE);
      }
    }
    
    /* We must always unshift here, also if we did not shake
     * x was shifted in do_force */
    if ((inputrec->ePBC == epbcXYZ) && (graph->nnodes > 0)) {
      unshift_x(graph,state->box,state->x,xprime);
      if (TRICLINIC(state->box))
	inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
      else
	inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);    
      for(n=start; (n<graph->start); n++)
	copy_rvec(xprime[n],state->x[n]);
      for(n=graph->start+graph->nnodes; (n<start+homenr); n++)
	copy_rvec(xprime[n],state->x[n]);
    } else {
      for(n=start; (n<start+homenr); n++)
	copy_rvec(xprime[n],state->x[n]);
    }
  } else {
    for(n=start; (n<start+homenr); n++)
      copy_rvec(xprime[n],shakefirst_x[n]);
  }
  dump_it_all(stdlog,"After unshift",
	      natoms,state->x,xprime,state->v,vold,force);
  where();

  if (bDoUpdate) {
    update_grps(start,homenr,grps,&(inputrec->opts),state->v,md,bNEMD);
    if (DEFORM(*inputrec))
      deform_store(state->box,inputrec,step,bFirstStep);
    if (inputrec->epc == epcBERENDSEN) {
      berendsen_pscale(state->pcoupl_mu,state->box,start,homenr,state->x,
		       md->cFREEZE,nrnb,inputrec->opts.nFreeze);
    } else if (inputrec->epc == epcPARRINELLORAHMAN) {
      /* The box velocities were updated in do_pr_pcoupl in the update
       * iteration, but we dont change the box vectors until we get here
       * since we need to be able to shift/unshift above.
       */
      for(i=0;i<DIM;i++)
        for(m=0;m<=i;m++)
          state->box[i][m] += dt*state->boxv[i][m];
    }
    if (DEFORM(*inputrec))
      deform(start,homenr,state->x,state->box,inputrec,step);
    where();
  }
}


void correct_ekin(FILE *log,int start,int end,rvec v[],rvec vcm,real mass[],
                  real tmass,tensor ekin)
{
  /* 
   * This is a debugging routine. It should not be called for production code
   *
   * The kinetic energy should calculated according to:
   *   Ekin = 1/2 m (v-vcm)^2
   * However the correction is not always applied, since vcm may not be
   * known in time and we compute
   *   Ekin' = 1/2 m v^2 instead
   * This can be corrected afterwards by computing
   *   Ekin = Ekin' + 1/2 m ( -2 v vcm + vcm^2)
   * or in hsorthand:
   *   Ekin = Ekin' - m v vcm + 1/2 m vcm^2
   */
  int    i,j,k;
  real   m,tm;
  rvec   hvcm,mv;
  tensor dekin;

  /* Local particles */  
  clear_rvec(mv);

  /* Processor dependent part. */
  tm = 0;
  for(i=start; (i<end); i++) {
    m      = mass[i];
    tm    += m;
    for(j=0; (j<DIM); j++)
      mv[j] += m*v[i][j];
  }
  /* Shortcut */ 
  svmul(1/tmass,vcm,vcm); 
  svmul(0.5,vcm,hvcm);
  clear_mat(dekin);
  for(j=0; (j<DIM); j++)
    for(k=0; (k<DIM); k++)
      dekin[j][k] += vcm[k]*(tm*hvcm[j]-mv[j]);

  pr_rvecs(log,0,"dekin",dekin,DIM);
  pr_rvecs(log,0," ekin", ekin,DIM);
  fprintf(log,"dekin = %g, ekin = %g  vcm = (%8.4f %8.4f %8.4f)\n",
          trace(dekin),trace(ekin),vcm[XX],vcm[YY],vcm[ZZ]);
  fprintf(log,"mv = (%8.4f %8.4f %8.4f)\n",
          mv[XX],mv[YY],mv[ZZ]);
}
