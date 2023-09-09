/*
 * $Id: group.h,v 1.12.2.1 2008/02/29 07:02:42 spoel Exp $
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

enum { egCOULSR, egLJSR, egBHAMSR, egCOULLR, egLJLR, egBHAMLR,
       egCOUL14, egLJ14, egNR };
	
typedef struct {
  real    Th;		/* Temperature at half step        */
  real    T;		/* Temperature at full step        */
  tensor  ekinh;	/* Kinetic energy at half step     */
  tensor  ekinh_old;	/* Kinetic energy at old half step */
  tensor  ekin; 	/* Kinetic energy at full step     */
  real    lambda;       /* Berendsen coupling lambda       */
} t_grp_tcstat;

typedef struct {
  int           nn;             /* Number of terms 			*/
  real 		*ee[egNR];	/* Arrays of energy terms for THIS 	*/
  				/* group with ALL other groups		*/
} t_grp_ener;

typedef struct {
  int     nat;		/* Number of atoms in this group		*/
  atom_id *aid;		/* Atom ids of the atoms in this group		*/
  real    M;		/* Total mass of group				*/
  rvec	  u;           	/* Mean velocities of home particles    	*/
} t_grp_acc;

typedef struct {
  real    cos_accel;    /* The acceleration for the cosine profile      */
  real    mvcos;        /* The cos momenta of home particles            */
  real    vcos;         /* The velocity of the cosine profile           */
} t_cos_acc;

typedef struct {
  t_grp_ener   estat;		/* Energy logging stuff			*/
  t_grp_tcstat *tcstat;         /* T-coupling data 			*/
  t_grp_acc    *grpstat;	/* Acceleration data			*/
  real         dekindl;         /* dEkin/dlambda at half step           */
  real         dekindl_old;     /* dEkin/dlambda at old half step       */
  t_cos_acc    cosacc;          /* Cosine acceleration data             */
} t_groups;

#define GID(igid,jgid,gnr) ((igid < jgid) ? (igid*gnr+jgid) : (jgid*gnr+igid))

