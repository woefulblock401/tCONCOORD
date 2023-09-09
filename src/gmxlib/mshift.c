/*
 * $Id: mshift.c,v 1.45.2.8 2008/02/29 07:02:44 spoel Exp $
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

#include <string.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "mshift.h"
#include "main.h"
#include "pbc.h"

/************************************************************
 *
 *      S H I F T   U T I L I T I E S
 *
 ************************************************************/
 

/************************************************************
 *
 *      G R A P H   G E N E R A T I O N    C O D E
 *
 ************************************************************/

static void add_gbond(t_graph *g,atom_id a0,atom_id a1)
{
  int     i,j,k;
  atom_id inda0,inda1,indb,indc;
  bool    bFound;

  inda0 = a0 - g->start;
  inda1 = a1 - g->start;
  bFound = FALSE;
  /* Search for a direct edge between a0 and a1.
   * All egdes are bidirectional, so we only need to search one way.
   */
  for(i=0; (i<g->nedge[inda0] && !bFound); i++)
    bFound = (g->edge[inda0][i] == a1);
  /* Search for connections via two or three edges. */
  for(i=0; (i<g->nedge[inda0] && !bFound); i++) {
    indb = g->edge[inda0][i] - g->start;
    for(j=0; (j<g->nedge[indb] && !bFound); j++) {
      bFound = (g->edge[indb][j] == a1);
      indc =  g->edge[indb][j] - g->start;
      for(k=0; (k<g->nedge[indc] && !bFound); k++)
	bFound = (g->edge[indc][k] == a1);
    }
  }

  if (!bFound) {
    g->edge[inda0][g->nedge[inda0]++] = a1;
    g->edge[inda1][g->nedge[inda1]++] = a0;
  }
}

static void mk_igraph(t_graph *g,t_functype ftype[],t_ilist *il,int natoms)
{
  t_iatom *ia;
  t_iatom tp;
  int     i,j,np;
  int     end;

  end=il->nr;
  ia=il->iatoms;

  i = 0;
  while (i < end) {
    tp = ftype[ia[0]];
    np = interaction_function[tp].nratoms;
    
    if (ia[1] < natoms) {
      if (ia[np] >= natoms)
	gmx_fatal(FARGS,
		  "Molecule in topology has atom numbers below and "
		  "above natoms (%d).\n"
		  "You are probably trying to use a trajectory which does "
		  "not match the first %d atoms of the run input file.\n"
		  "You can make a matching run input file with tpbconv.",
		  natoms,natoms);
      if (tp == F_SETTLE) {
	/* Bond all the atoms in the settle */
	add_gbond(g,ia[1],ia[1]+1);
	add_gbond(g,ia[1],ia[1]+2);
      } else {
	for(j=1; j<np; j++)
	  add_gbond(g,ia[j],ia[j+1]);
      }
    }
    ia+=np+1;
    i+=np+1;
  }
}

static void g_error(int line,char *file)
{
  gmx_fatal(FARGS,"Tring to print non existant graph (file %s, line %d)",
	      file,line);
}
#define GCHECK(g) if (g == NULL) g_error(__LINE__,__FILE__)

void p_graph(FILE *log,char *title,t_graph *g)
{
  int  i,j;
  char *cc[egcolNR] = { "W", "G", "B" };
  
  GCHECK(g);
  fprintf(log,"graph:  %s\n",title);
  fprintf(log,"nnodes: %d\n",g->nnodes);
  fprintf(log,"nbound: %d\n",g->nbound);
  fprintf(log,"start:  %d\n",g->start);
  fprintf(log,"end:    %d\n",g->end);
  fprintf(log," atom shiftx shifty shiftz C nedg    e1    e2 etc.\n");
  for(i=0; (i<g->nnodes); i++)
    if (g->nedge[i] > 0) {
      fprintf(log,"%5d%7d%7d%7d %1s%5d",g->start+i+1,
	      g->ishift[i][XX],g->ishift[i][YY],
	      g->ishift[i][ZZ],
	      (g->negc > 0) ? cc[g->egc[i]] : " ",
	      g->nedge[i]);
      for(j=0; (j<g->nedge[i]); j++)
	fprintf(log," %5u",g->edge[i][j]+1);
      fprintf(log,"\n");
    }
  fflush(log);
}

static void calc_1se(t_graph *g,t_ilist *il,t_functype ftype[],
		     int nbond[],int natoms)
{
  int     k,nratoms,tp,end,j;
  t_iatom *ia,iaa;

  end=il->nr;

  ia=il->iatoms;
  for(j=0; (j<end); j+=nratoms+1,ia+=nratoms+1) {
    tp      = ftype[ia[0]];
    nratoms = interaction_function[tp].nratoms;
    
    if (tp == F_SETTLE) {
      iaa          = ia[1];
      if (iaa<natoms) {
	nbond[iaa]   += 2;
	nbond[iaa+1] += 1;
	nbond[iaa+2] += 1;
	g->start      = min(g->start,iaa);
	g->end        = max(g->end,iaa+2);
      }
    } else {
      for(k=1; (k<=nratoms); k++) {
	iaa=ia[k];
	if (iaa<natoms) {
	  g->start=min(g->start,iaa);
	  g->end  =max(g->end,  iaa);
	  /*if (interaction_function[tp].flags & IF_CHEMBOND)*/
	  nbond[iaa]++;
	}
      }
    }
  }
}

static int calc_start_end(t_graph *g,t_idef *idef,int natoms,
			  int nbond[])
{
  int   i,nnb,nbtot;
  
  g->start=natoms;
  g->end=0;

  /* First add all the real bonds: they should determine the molecular 
   * graph.
   */
  for(i=0; (i<F_NRE); i++)
    if (interaction_function[i].flags & IF_CHEMBOND)
      calc_1se(g,&idef->il[i],idef->functype,nbond,natoms);
  /* Then add all the other interactions in fixed lists, but first
   * check to see what's there already.
   */
  for(i=0; (i<F_NRE); i++) {
    if (!(interaction_function[i].flags & IF_CHEMBOND)) {
      calc_1se(g,&idef->il[i],idef->functype,nbond,natoms);
    }
  }
  
  nnb   = 0;
  nbtot = 0;
  for(i=g->start; (i<=g->end); i++) {
    nbtot += nbond[i];
    nnb    = max(nnb,nbond[i]);
  }
  if (stdlog) {
    fprintf(stdlog,"Max number of connections per atom is %d\n",nnb);
    fprintf(stdlog,"Total number of connections is %d\n",nbtot);
  }
  return nbtot;
}



static void compact_graph(t_graph *g)
{
  int i,j,n,max_nedge;
  atom_id *e;

  n = 0;
  e = g->edge[0];
  max_nedge = 0;
  for(i=0; i<g->nnodes; i++) {
    for(j=0; j<g->nedge[i]; j++)
      e[n++] = g->edge[i][j];
    max_nedge = max(max_nedge,g->nedge[i]);
  }
  srenew(g->edge[0],n);
  for(i=1; i<g->nnodes; i++) {
    g->edge[i] = g->edge[i-1] + g->nedge[i-1];
  }
  if (stdlog) {
    fprintf(stdlog,"Max number of graph edges per atom is %d\n",
	    max_nedge);
    fprintf(stdlog,"Total number of graph edges is %d\n",n);
  }
}

t_graph *mk_graph(t_idef *idef,int natoms,bool bShakeOnly,bool bSettle)
{
  t_graph *g;
  int     *nbond;
  int     i,nbtot;
  
  snew(g,1);

  snew(nbond,natoms);
  nbtot = calc_start_end(g,idef,natoms,nbond);
  
  if (g->start >= g->end) {
    g->nnodes=0;
  }
  else {
    g->nnodes=g->end-g->start+1;
    snew(g->ishift,g->nnodes);
    snew(g->nedge,g->nnodes);
  
    /* To prevent malloc problems with many small arrays, using realloc,
     * we allocate some more memory, and divide it ourselves.
     * We calculate pointers... (Yuck Yuck)
     */
    snew(g->edge,g->nnodes);
    snew(g->edge[0],nbtot);

    for(i=1; (i<g->nnodes); i++)
      g->edge[i]=g->edge[i-1]+nbond[g->start+i-1];

    if (!bShakeOnly) {
      /* First add all the real bonds: they should determine the molecular 
       * graph.
       */
      for(i=0; (i<F_NRE); i++)
	if (interaction_function[i].flags & IF_CHEMBOND)
	  mk_igraph(g,idef->functype,&(idef->il[i]),natoms);
      /* Then add all the other interactions in fixed lists, but first
       * check to see what's there already.
       */
      for(i=0; (i<F_NRE); i++) {
	if (!(interaction_function[i].flags & IF_CHEMBOND)) {
	  mk_igraph(g,idef->functype,&(idef->il[i]),natoms);
	}
      }
      
      /* Removed all the unused space from the edge array */
      compact_graph(g);
    }
    else {
      /* This is a special thing used in grompp to generate shake-blocks */
      mk_igraph(g,idef->functype,&(idef->il[F_SHAKE]),natoms);
      if (bSettle)
	mk_igraph(g,idef->functype,&(idef->il[F_SETTLE]),natoms);
    }
    g->nbound=0;
    for(i=0; (i<g->nnodes); i++)
      if (g->nedge[i] > 0)
        g->nbound++;
  }
  if (debug)
    p_graph(debug,"graph",g);

  g->negc = 0;
  g->egc = NULL;
  
  sfree(nbond);
  
  return g;
}

void done_graph(t_graph *g)
{
  GCHECK(g);
  if (g->nnodes > 0) {
    sfree(g->ishift);
    sfree(g->nedge);
    /* This is malloced in a NASTY way, see above */
    sfree(g->edge[0]);
    sfree(g->edge);
    sfree(g->egc);
  }
}

/************************************************************
 *
 *      S H I F T   C A L C U L A T I O N   C O D E
 *
 ************************************************************/

static void mk_1shift_tric(matrix box,rvec hbox,rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for triclinic box... */
  int  m,d;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  for(m=DIM-1; (m>=0); m--) {
    /* If dx < hbox, then xj will be reduced by box, so that
     * xi - xj will be bigger
     */
    if (dx[m] < -hbox[m]) {
      mj[m]=mi[m]-1;
      for(d=m-1; d>=0; d--)
	dx[d]+=box[m][d];
    } else if (dx[m] >= hbox[m]) {
      mj[m]=mi[m]+1;
      for(d=m-1; d>=0; d--)
	dx[d]-=box[m][d];
    } else
      mj[m]=mi[m];
  }
}

static void mk_1shift(rvec hbox,rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for rectangular box... */
  int  m;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  for(m=0; (m<DIM); m++) {
    /* If dx < hbox, then xj will be reduced by box, so that
     * xi - xj will be bigger
     */
    if (dx[m] < -hbox[m])
      mj[m]=mi[m]-1;
    else if (dx[m] >= hbox[m])
      mj[m]=mi[m]+1;
    else
      mj[m]=mi[m];
  }
}

static int mk_grey(FILE *log,int nnodes,egCol egc[],t_graph *g,int *AtomI,
		   matrix box,rvec x[],int *nerror)
{
  int      m,j,ng,ai,aj,g0;
  rvec     hbox;
  bool     bTriclinic;
  ivec     is_aj;
  
  for(m=0; (m<DIM); m++)
    hbox[m]=box[m][m]*0.5;
  bTriclinic = TRICLINIC(box);
  
  ng=0;
  ai=*AtomI;
  
  g0=g->start;
  /* Loop over all the bonds */
  for(j=0; (j<g->nedge[ai]); j++) {
    aj=g->edge[ai][j]-g0;
    /* If there is a white one, make it gray and set pbc */
    if (bTriclinic)
      mk_1shift_tric(box,hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    else
      mk_1shift(hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    
    if (egc[aj] == egcolWhite) {
      if (aj < *AtomI)
	*AtomI = aj;
      egc[aj] = egcolGrey;
      
      copy_ivec(is_aj,g->ishift[aj]);

      ng++;
    }
    else if ((is_aj[XX] != g->ishift[aj][XX]) ||
	     (is_aj[YY] != g->ishift[aj][YY]) ||
	     (is_aj[ZZ] != g->ishift[aj][ZZ])) {
      (*nerror)++;
    }
  }
  return ng;
}

static int first_colour(int fC,egCol Col,t_graph *g,egCol egc[])
/* Return the first node with colour Col starting at fC.
 * return -1 if none found.
 */
{
  int i;
  
  for(i=fC; (i<g->nnodes); i++)
    if ((g->nedge[i] > 0) && (egc[i]==Col))
      return i;
  
  return -1;
}

void mk_mshift(FILE *log,t_graph *g,matrix box,rvec x[])
{
  static int nerror_tot = 0;
  int    ng,nnodes,i;
  int    nW,nG,nB;		/* Number of Grey, Black, White	*/
  int    fW,fG;			/* First of each category	*/
  int    nerror=0;

  GCHECK(g);
  /* This puts everything in the central box, that is does not move it 
   * at all. If we return without doing this for a system without bonds
   * (i.e. only settles) all water molecules are moved to the opposite octant
   */
  for(i=0; (i<g->nnodes); i++) {
      g->ishift[i][XX]=g->ishift[i][YY]=g->ishift[i][ZZ]=0;
  }
    
  if (!g->nbound)
    return;

  nnodes=g->nnodes;
  if (nnodes > g->negc) {
    g->negc = nnodes;
    srenew(g->egc,g->negc);
  }
  memset(g->egc,0,(size_t)(nnodes*sizeof(g->egc[0])));

  nW=g->nbound;
  nG=0;
  nB=0;

  fW=0;

  /* We even have a loop invariant:
   * nW+nG+nB == g->nbound
   */
#ifdef DEBUG2
  fprintf(log,"Starting W loop\n");
#endif
  while (nW > 0) {
    /* Find the first white, this will allways be a larger
     * number than before, because no nodes are made white
     * in the loop
     */
    if ((fW=first_colour(fW,egcolWhite,g,g->egc)) == -1) 
      gmx_fatal(FARGS,"No WHITE nodes found while nW=%d\n",nW);
    
    /* Make the first white node grey */
    g->egc[fW]=egcolGrey;
    nG++;
    nW--;

    /* Initial value for the first grey */
    fG=fW;
#ifdef DEBUG2
    fprintf(log,"Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n",
	    nW,nG,nB,nW+nG+nB);
#endif
    while (nG > 0) {
      if ((fG=first_colour(fG,egcolGrey,g,g->egc)) == -1)
	gmx_fatal(FARGS,"No GREY nodes found while nG=%d\n",nG);
      
      /* Make the first grey node black */
      g->egc[fG]=egcolBlack;
      nB++;
      nG--;

      /* Make all the neighbours of this black node grey
       * and set their periodicity 
       */
      ng=mk_grey(log,nnodes,g->egc,g,&fG,box,x,&nerror);
      /* ng is the number of white nodes made grey */
      nG+=ng;
      nW-=ng;
    }
  }
  if (nerror > 0) {
    nerror_tot++;
    if (nerror_tot <= 100)
      fprintf(log,"There were %d inconsistent shifts. Check your topology\n",
	      nerror);
    if (nerror_tot == 100)
      fprintf(log,"Will stop reporting inconsistent shifts\n");
  }
}

/************************************************************
 *
 *      A C T U A L   S H I F T   C O D E
 *
 ************************************************************/
 
static void shift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai)
{
  int tx,ty,tz;
  
  tx=(g->ishift[ai-g->start])[XX];
  ty=(g->ishift[ai-g->start])[YY];
  tz=(g->ishift[ai-g->start])[ZZ];

  x_s[ai][XX]=x[ai][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
  x_s[ai][YY]=x[ai][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
  x_s[ai][ZZ]=x[ai][ZZ]+tz*box[ZZ][ZZ];
}
 
static void unshift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai)
{
  int tx,ty,tz;
  
  tx=(g->ishift[ai-g->start])[XX];
  ty=(g->ishift[ai-g->start])[YY];
  tz=(g->ishift[ai-g->start])[ZZ];

  x_s[ai][XX]=x[ai][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
  x_s[ai][YY]=x[ai][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
  x_s[ai][ZZ]=x[ai][ZZ]-tz*box[ZZ][ZZ];
}

void shift_x(t_graph *g,matrix box,rvec x[],rvec x_s[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  GCHECK(g);
  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  if(TRICLINIC(box)) {
     for(i=0,j=g0; (i<gn); i++,j++) { 
	 tx=is[i][XX];
	 ty=is[i][YY];
	 tz=is[i][ZZ];
	 
	 x_s[j][XX]=x[j][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
	 x_s[j][YY]=x[j][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
	 x_s[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
     }
  } else {
     for(i=0,j=g0; (i<gn); i++,j++) { 
	 tx=is[i][XX];
	 ty=is[i][YY];
	 tz=is[i][ZZ];
	 
	 x_s[j][XX]=x[j][XX]+tx*box[XX][XX];
	 x_s[j][YY]=x[j][YY]+ty*box[YY][YY];
	 x_s[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
     }
  }       
     
}

void shift_self(t_graph *g,matrix box,rvec x[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;

#ifdef DEBUG
  fprintf(stderr,"Shifting atoms %d to %d\n",g0,g0+gn);
#endif
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) { 
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
	  x[j][YY]=x[j][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
	  x[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) { 
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]+tx*box[XX][XX];
	  x[j][YY]=x[j][YY]+ty*box[YY][YY];
	  x[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
      }
  }       
  
}

void unshift_x(t_graph *g,matrix box,rvec x[],rvec x_s[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x_s[j][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
	  x[j][YY]=x_s[j][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
	  x[j][ZZ]=x_s[j][ZZ]-tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x_s[j][XX]-tx*box[XX][XX];
	  x[j][YY]=x_s[j][YY]-ty*box[YY][YY];
	  x[j][ZZ]=x_s[j][ZZ]-tz*box[ZZ][ZZ];
      }
  }
}

void unshift_self(t_graph *g,matrix box,rvec x[])
{
  ivec *is;
  int g0,gn;
  int i,j,tx,ty,tz;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
	  x[j][YY]=x[j][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
	  x[j][ZZ]=x[j][ZZ]-tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]-tx*box[XX][XX];
	  x[j][YY]=x[j][YY]-ty*box[YY][YY];
	  x[j][ZZ]=x[j][ZZ]-tz*box[ZZ][ZZ];
      }
  }
}
#undef GCHECK

#ifdef DEBUGMSHIFT
void main(int argc,char *argv[])
{
  FILE         *out;
  t_args       targ;
  t_topology   top;
  t_statheader sh;
  rvec         *x;
  ivec         *mshift;
  matrix       box;

  t_graph      *g;
  int          i,idum,pid;
  real         rdum;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,&targ,PCA_NEED_INOUT,NULL);
  if (argc > 1)
    pid=atoi(argv[1]);
  else
    pid=0;
  
  read_status_header(targ.infile,&sh);
  snew(x,sh.natoms);
  snew(mshift,sh.natoms);

  fprintf(stderr,"Reading Status %s\n",targ.infile);
  read_status(targ.infile,&idum,&rdum,&rdum,NULL,
	      box,NULL,NULL,&idum,x,NULL,NULL,&idum,NULL,&top);

  fprintf(stderr,"Making Graph Structure...\n");
  g=mk_graph(&(top.idef),top.atoms.nr,FALSE,FALSE);

  out=ffopen(targ.outfile,"w");

  fprintf(stderr,"Making Shift...\n");
  mk_mshift(out,g,box,x,mshift);

  p_graph(out,"In Den Haag daar woont een graaf...",g);
  fclose(out);
}
#endif

