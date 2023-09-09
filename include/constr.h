/*
 * $Id: constr.h,v 1.23.2.1 2008/02/29 07:02:41 spoel Exp $
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
#include<config.h>
#endif

extern bool bshakef(FILE *log,		/* Log file			*/
		    int natoms,		/* Total number of atoms	*/
		    real invmass[],	/* Atomic masses		*/
		    int nblocks,	/* The number of shake blocks	*/
		    int sblock[],       /* The shake blocks             */
		    t_idef *idef,	/* The interaction def		*/
		    t_inputrec *ir,	/* Input record		        */
		    matrix box,		/* The box			*/
		    rvec x_s[],		/* Coords before update		*/
		    rvec xp[],		/* Output coords		*/
		    t_nrnb *nrnb,       /* Performance measure          */
		    real lambda,        /* FEP lambda                   */
		    real *dvdlambda,    /* FEP force                    */
		    bool bCalcVir,      /* Calculate r x m delta_r      */
		    tensor rmdr,        /* sum r x m delta_r            */
		    bool bDumpOnError); /* Dump debugging stuff on error*/
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0 
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from 
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return TRUE when OK, FALSE when shake-error
 */
extern void csettle(FILE *log,
		    int nshake,		/* Number of water molecules 	*/
		    int owptr[],	/* pointer to Oxygen in b4 & after */
		    real b4[],		/* Old coordinates		*/
		    real after[],	/* New coords, to be settled	*/
		    real dOH,		/* Constraint length Ox-Hyd	*/
		    real dHH, 		/* Constraint length Hyd-Hyd	*/
		    real mO,  		/* Mass of Oxygen		*/
		    real mH, 		/* Mass of Hydrogen		*/
		    bool bCalcVir,      /* Calculate r x m delta_r      */
		    tensor rmdr,        /* sum r x m delta_r            */
		    int *xerror);

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],real omega,
		   real invmass[],real tt[],real lagr[],int *nerror);
/* Regular iterative shake */

extern bool constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
		      t_mdatoms *md,int start,int homenr,
		      rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		      real lambda,real *dvdlambda,tensor *vir,
		      t_nrnb *nrnb,bool bCoordinates);
/*
 * When bCoordinates=TRUE constrains coordinates xprime using th
 * directions in x, min_proj is not used.
 *
 * When bCoordinates=FALSE, calculates the components xprime in
 * the constraint directions and subtracts these components from min_proj.
 * So when min_proj=xprime, the constraint components are projected out.
 *
 * If vir!=NULL calculate the constraint virial.
 *
 * Init_constraints must have be called once, before calling constrain.
 *
 * Return TRUE if OK, FALSE in case of shake error
 *
 */

extern int count_constraints(t_topology *top,t_commrec *cr);
/* Returns the total number of constraints in the system */

extern int init_constraints(FILE *log,t_topology *top,t_inputrec *ir,
			    t_mdatoms *md,int start,int homenr,
			    bool bOnlyCoords,t_commrec *cr);
/* Initialize constraints stuff */


/* LINCS stuff */
typedef struct {
  int  nc;       /* the number of constraints */
  int  nzerolen; /* the number of constraints with zero length */
  int  ncc;      /* the number of constraint connections */
  real matlam;   /* the FE lambda value used for filling blc and blcc */
  real *bllen0;  /* the reference distance in topology A */
  real *ddist;   /* the reference distance in top B - the r.d. in top A */
  int  *bla;     /* the atom pairs involved in the constraints */
  real *blc;     /* 1/sqrt(invmass1 + invmass2) */
  int  *blnr;    /* index into blbnb and blcc */
  int  *blbnb;   /* list of bond connections */
  real *blcc;    /* bond coupling coefficient matrix */
  real *bllen;   /* the reference bond length */
  /* arrays for temporary storage in the LINCS algorithm */
  rvec *tmpv;
  real *tmpncc;
  real *tmp1;
  real *tmp2;
  real *tmp3;
  real *lambda;  /* the Lagrange multipliers */
} t_lincsdata;

extern t_lincsdata *init_lincs(FILE *log,t_idef *idef,int start,int homenr,
			       bool bDynamics);
/* Initialize lincs stuff */

extern void set_lincs_matrix(t_lincsdata *li,real *invmass);
/* Sets the elements of the LINCS constraint coupling matrix */

extern bool constrain_lincs(FILE *log,t_inputrec *ir,
			    int step,t_lincsdata *lincsd,t_mdatoms *md,
			    rvec *x,rvec *xprime,rvec *min_proj,matrix box,
			    real lambda,real *dvdlambda,
			    bool bCalcVir,tensor rmdr,
			    bool bCoordinates,
			    t_nrnb *nrnb,bool bDumpOnError);
/* Returns if the constraining succeeded */
