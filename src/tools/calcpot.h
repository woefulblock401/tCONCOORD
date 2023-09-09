/*
 * $Id: calcpot.h,v 1.12.2.2 2008/02/29 07:02:54 spoel Exp $
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
	
extern void init_calcpot(char *log,char *tpx,char *table,
			 t_topology *top,t_inputrec *inputrec,t_commrec *cr,
			 t_graph **graph,t_mdatoms **mdatoms,
			 t_nsborder *nsb,t_groups *grps,
			 t_forcerec **fr,real **coulomb,
			 matrix box,rvec **x);

extern void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
		     t_inputrec *inputrec,t_topology *top,rvec x[],t_forcerec *fr,
		     t_mdatoms *mdatoms,real coulomb[],matrix box,t_graph *graph);

extern void write_pdb_coul();

extern void delete_atom(t_topology *top,int inr);
/* Delete an atom from a topology */

extern void replace_atom(t_topology *top,int inr,char *anm,char *resnm,
			 real q,real m,int type);
/* Replace an atom in a topology by someting else */

