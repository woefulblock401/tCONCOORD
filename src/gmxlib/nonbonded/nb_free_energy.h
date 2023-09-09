/*
 * $Id: nb_free_energy.h,v 1.3.2.2 2008/02/29 07:02:44 spoel Exp $
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

#ifndef _nb_free_energy_h
#define _nb_free_energy_h

#include <typedefs.h>

void
gmx_nb_free_energy_kernel(int                  icoul,
                          int                  ivdw,
                          int                  nri,
                          int *                iinr,
                          int *                jindex,
                          int *                jjnr,
                          int *                shift,
                          real *               shiftvec,
                          real *               fshift,
                          int *                gid,
                          real *               x,
                          real *               f,
                          real *               chargeA,
                          real *               chargeB,
                          real                 facel,
                          real                 krf,
                          real                 crf,
			  real                 ewc,
                          real *               Vc,
                          int *                typeA,
                          int *                typeB,
                          int                  ntype,
                          real *               nbfp,
                          real *               Vvdw,
                          real                 tabscale,
                          real *               VFtab,
                          real                 lambda,
                          real *               dvdlambda,
                          real                 alpha,
			  int                  lam_power,
                          real                 def_sigma6,
                          int *                outeriter,
                          int *                inneriter);

#endif

