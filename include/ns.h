/*
 * $Id: ns.h,v 1.17.2.1 2008/02/29 07:02:42 spoel Exp $
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

#ifndef _ns_h
#define _ns_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "pbc.h"
#include "tgroup.h"
#include "nsb.h"
#include "network.h"

/****************************************************
 *
 *    U T I L I T I E S May be found in ns.c
 *
 ****************************************************/
extern void correct_box(tensor box,t_forcerec *fr,t_graph *g);
/* Corrects the box by subtracting box vector when the box is too skewed.
 * The integer shift vectors in the graph and 
 * the shift indices of the short-range neighborlists are changed accordingly.
 */

extern void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr);
/* 
 * nn is the number of energy terms in the energy matrix
 * (ngener*(ngener-1))/2
 * start is the first atom on this processor
 * homenr is the number of atoms on this processor
 */
 
extern int calc_naaj(int icg,int cgtot);
/* Calculate the number of charge groups to interact with for icg */

/****************************************************
 *
 *    N E I G H B O R  S E A R C H I N G
 *
 *    Calls either ns5_core (when grid selected in .mdp file)
 *    or ns_simple_core (when simple selected in .mdp file)
 *
 *    Return total number of pairs searched 
 *
 ****************************************************/
extern int search_neighbours(FILE *log,t_forcerec *fr,
			     rvec x[],matrix box,
			     t_topology *top,t_groups *grps,
			     t_commrec *cr,t_nsborder *nsb,t_nrnb *nrnb,
			     t_mdatoms *md,real lambda,real *dvdlambda,
			     bool bFillGrid,bool bDoForces);
 

#endif	/* _ns_h */
