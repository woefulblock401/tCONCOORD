/*
 * $Id: gmx_parallel_3dfft.h,v 1.1.2.3 2008/02/29 07:02:41 spoel Exp $
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

#ifndef _gmx_parallel_3dfft_h
#define _gmx_parallel_3dfft_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_MPI

#include "types/simple.h"
#include "gmxcomplex.h"
#include "gmx_fft.h"

/* We NEED MPI here. */
#include <mpi.h>

typedef struct gmx_parallel_3dfft *
gmx_parallel_3dfft_t;



/*! \brief Initialize parallel MPI-based 3D-FFT.
 *
 *  This routine performs real-to-complex and complex-to-real parallel 3D FFTs,
 *  but not complex-to-complex.
 *
 *  The routine is optimized for small-to-medium size FFTs used for PME and
 *  PPPM algorithms, and do allocate extra workspace whenever it might improve
 *  performance. 
 *
 *  \param pfft_setup  Pointer to parallel 3dfft setup structure, previously
 *                     allocated or with automatic storage.
 *  \param ngridx      Global number of grid cells in the x direction. Must be
 *                     divisible by the number of nodes.
 *  \param ngridy      Global number of grid cells in the y direction. Must be
 *                     divisible by the number of nodes.
 *  \param ngridz      Global number of grid cells in the z direction.
 *  \param comm        MPI communicator, must have been initialized. 
 *    
 *  \return 0 or a standard error code.
 */
int
gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                           int                       ngridx,
                           int                       ngridy,
                           int                       ngridz,
                           MPI_Comm                  comm);
                           




/*! \brief Get direct space grid index limits
 *
 *  The z dimension is never distributed. In the direct space, the x dimension
 *  is distributed over nodes, and after the real-to-complex FFT we work with
 *  a transposed grid where the y dimension is partitioned over nodes.
 *
 *  The order is determined from the rank in the communicator used at
 *  initialization.
 */
int
gmx_parallel_3dfft_limits(gmx_parallel_3dfft_t      pfft_setup,
                          int *                     local_x_start,
                          int *                     local_nx,
                          int *                     local_y_start,
                          int *                     local_ny);


int
gmx_parallel_transpose(t_complex *   data,
                       t_complex *   work,
                       int           nx,
                       int           ny,
                       int           local_x_start,
                       int           local_nx,
                       int           local_y_start,
                       int           local_ny,
                       int           nelem,
                       MPI_Comm      comm);


/*! \brief Perform forward parallel MPI FFT.
 *
 *  Direction is either GMX_FFT_REAL_TO_COMPLEX or GMX_FFT_COMPLEX_TO_REAL.
 *
 *  If input and output arrays are separate there is no packing to consider.
 *  Input is simply nx*ny*nz in real, and output ny*nx*nzc in complex.
 *
 *  In they are identical we need to make sure there is room for the complex
 *  (length nzc=nz/2+1) in the array, so the _real_ space dimensions is
 *  always padded to nzc*2.
 *  In this case, the real dimensions are nx*ny*(nzc*2) while the complex
 *  dimensions is ny*nx*nzc (of type complex).
 *
 *  Note that the X and Y dimensions are transposed in the reciprocal space
 *  to avoid extra communication!
 */
int
gmx_parallel_3dfft(gmx_parallel_3dfft_t    pfft_setup,
                   enum gmx_fft_direction  dir,
                   void *                  in_data,
                   void *                  out_data);



/*! \brief Release all data in parallel fft setup
 *
 *  All temporary storage and FFT plans are released. The structure itself
 *  is not released, but the contents is invalid after this call.
 *
 *  \param pfft_setup Parallel 3dfft setup.
 *
 *  \return 0 or a standard error code.
 */
int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup);

#endif /* GMX_MPI */

#endif /* _gmx_parallel_3dfft_h_ */

