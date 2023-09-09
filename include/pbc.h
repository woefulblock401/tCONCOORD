/*
 * $Id: pbc.h,v 1.34.2.3 2008/02/29 07:02:42 spoel Exp $
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

#ifndef _pbc_h
#define _pbc_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"

#ifdef CPLUSPLUS
extern "C" { 
#endif

#define BOX_MARGIN 0.5001
  /* margin factor for checking if the box is too skewed */

#define TRICLINIC(box) (box[YY][XX]!=0 || box[ZZ][XX]!=0 || box[ZZ][YY]!=0)

#define NTRICIMG 14
#define NCUCVERT 24
#define NCUCEDGE 36

  enum {
    ecenterTRIC, /* 0.5*(a+b+c)                  */
    ecenterRECT, /* (0.5*a[x],0.5*b[y],0.5*c[z]) */
    ecenterZERO, /* (0,0,0)                      */
    ecenterDEF = ecenterTRIC
  };

  extern void dump_pbc(FILE *fp,t_pbc *pbc);
  /* Dump the contents of the pbc structure to the file */
  
  extern char *check_box(matrix box);
  /* Returns NULL if the box is supported by Gromacs.
   * Otherwise is returns a string with the problem.
   */

  extern real max_cutoff2(matrix box);
  /* Returns the square of the maximum cut-off allowed for the box,
   * taking into account that the grid neighborsearch code and pbc_dx
   * only check combinations of single box-vector shifts.
   */
  
  extern void set_pbc(t_pbc *pbc,matrix box);
  /* Initiate the periodic boundary conditions.
   * pbc_dx will not use pbc and return the normal difference vector
   * when one or more of the diagonal elements of box is zero.
   */
  extern void set_pbc_ss(t_pbc *pbc,matrix box);
  /* As pbc_dx, but additionally sets that correct distances can be
   * obtained using (combinations of) single box-vector shifts.
   * In this case pbc_dx is slightly more efficient.
   */

  extern int pbc_dx(const t_pbc *pbc,const rvec x1, const rvec x2, rvec dx);
  /* Calculate the correct distance vector from x2 to x1 and put it in dx.
   * Returns the ishift required to shift x1 at closest distance to x2;
   * i.e. if 0<=ishift<SHIFTS then x1 - x2 + shift_vec[ishift] = dx
   * (see calc_shifts below on how to obtain shift_vec)
   * set_pbc must be called before ever calling this routine.
   *
   * For triclinic boxes pbc_dx does not necessarily return the shortest
   * distance vector. If pbc->bLimitDistance=TRUE an atom pair with
   * distance vector dx with norm2(dx) > pbc->limit_distance2 could
   * have a shorter distance, but not shorter than sqrt(pbc->limit_distance2).
   * pbc->limit_distance2 is always larger than max_cutoff2(box).
   * For the standard rhombic dodecahedron and truncated octahedron
   * pbc->bLimitDistance=FALSE and thus all distances are correct.
   */

  extern bool image_rect(ivec xi,ivec xj,ivec box_size,
			 real rlong2,int *shift,real *r2);
  /* Calculate the distance between xi and xj for a rectangular box.
   * When the distance is SMALLER than rlong2 return TRUE, return
   * the shift code in shift and the distance in r2. When the distance is
   * >= rlong2 return FALSE;
   * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
   */

  extern bool image_tri(ivec xi,ivec xj,imatrix box,
			real rlong2,int *shift,real *r2);
  /* Calculate the distance between xi and xj for a triclinic box.
   * When the distance is SMALLER than rlong2 return TRUE, return
   * the shift code in shift and the distance in r2. When the distance is
   * >= rlong2 return FALSE;
   * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
   */
  
  extern bool image_cylindric(ivec xi,ivec xj,ivec box_size,real rlong2,
			      int *shift,real *r2);
  /* Calculate the distance between xi and xj for a rectangular box
   * using a cylindric cutoff for long-range only.
   * When the distance is SMALLER than rlong2 (in X and Y dir.)
   * return TRUE, return
   * the shift code in shift and the distance in r2. When the distance is
   * >= rlong2 return FALSE;
   * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
   */

  extern void calc_shifts(matrix box,rvec shift_vec[]);
  /* This routine calculates ths shift vectors necessary to use the
   * ns routine.
   */

  extern void calc_cgcm(FILE *log,int cg0,int cg1,t_block *cgs,
			rvec pos[],rvec cg_cm[]);
  /* Routine to compute centers of geometry of charge groups. No periodicity
   * is used.
   */
  
  extern void put_charge_groups_in_box (FILE *log,int cg0,int cg1,
					matrix box,t_block *cgs,
					rvec pos[],
					rvec cg_cm[]);
			    
  /* This routine puts charge groups in the periodic box, keeping them
   * together.
   */

  extern void calc_box_center(int ecenter,matrix box,rvec box_center);
  /* Calculates the center of the box.
   * See the description for the enum ecenter above.
   */

  extern void calc_triclinic_images(matrix box,rvec img[]);
  /* Calculates the NTRICIMG box images */

  extern void calc_compact_unitcell_vertices(int ecenter,matrix box,
					     rvec vert[]);
  /* Calculates the NCUCVERT vertices of a compact unitcell */
  
  extern int *compact_unitcell_edges(void);
  /* Return an array of unitcell edges of length NCUCEDGE*2,
   * this is an index in vert[], which is calculated by calc_unitcell_vertices.
   * The index consists of NCUCEDGE pairs of vertex indices.
   * The index does not change, so it needs to be retrieved only once.
   */
  extern void put_atom_in_box(matrix box,rvec x);

  extern void put_atoms_in_box(matrix box,int natoms,rvec x[]);
  /* These routines puts ONE or ALL atoms in the box, not caring 
   * about charge groups!
   * Also works for triclinic cells.
   */
  
  extern void put_atoms_in_triclinic_unitcell(int ecenter,matrix box,
					      int natoms,rvec x[]);
  /* This puts ALL atoms in the triclinic unit cell, centered around the
   * box center as calculated by calc_box_center.
   */

  extern char *put_atoms_in_compact_unitcell(int ecenter,matrix box,
					     int natoms,rvec x[]);
  /* This puts ALL atoms at the closest distance for the center of the box
   * as calculated by calc_box_center.
   * Will return NULL is everything went ok and a warning string if not
   * all atoms could be placed in the unitcell. This can happen for some
   * triclinic unitcells, see the comment at pbc_dx above.
   */
  
#ifdef CPLUSPLUS
}
#endif

#endif	/* _pbc_h */
