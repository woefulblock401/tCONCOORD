/*
 * $Id: toppush.c,v 1.44.2.18 2008/02/29 07:02:51 spoel Exp $
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

#include <ctype.h>
#include <math.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "string2.h"
#include "names.h"
#include "toputil.h"
#include "toppush.h"
#include "topdirs.h"
#include "symtab.h"
#include "gmx_fatal.h"


void generate_nbparams(int comb,int ftype,t_params *plist,t_atomtype *atype)
{
  int   i,j,k=-1,nf;
  int   nr,nrfp;
  real  c,bi,bj;
  char  errbuf[256];

  /* Lean mean shortcuts */
  nr   = atype->nr;
  nrfp = NRFP(ftype);
  snew(plist->param,nr*nr);
  plist->nr=nr*nr;
  
  /* Fill the matrix with force parameters */
  switch(ftype) {
  case F_LJ:
    switch (comb) {
    case eCOMB_GEOMETRIC:
      /* Gromos rules */
      for(i=k=0; (i<nr); i++)
	for(j=0; (j<nr); j++,k++)
	  for(nf=0; (nf<nrfp); nf++) {
	    c = sqrt(atype->nb[i].c[nf] * atype->nb[j].c[nf]);
	    plist->param[k].c[nf]      = c;
	  }
      break;
    
    case eCOMB_ARITHMETIC:
      /* c0 and c1 are epsilon and sigma */
      for (i=k=0; (i < nr); i++)
	for (j=0; (j < nr); j++,k++) {
	  plist->param[k].c[0] = (atype->nb[i].C0+atype->nb[j].C0)*0.5;
	  plist->param[k].c[1] = sqrt(atype->nb[i].C1*atype->nb[j].C1);
	}
      
      break;
    case eCOMB_GEOM_SIG_EPS:
      /* c0 and c1 are epsilon and sigma */
      for (i=k=0; (i < nr); i++)
	for (j=0; (j < nr); j++,k++) {
	  plist->param[k].c[0] = sqrt(atype->nb[i].C0*atype->nb[j].C0);
	  plist->param[k].c[1] = sqrt(atype->nb[i].C1*atype->nb[j].C1);
	}
      
      break;
    default:
      gmx_fatal(FARGS,"No such combination rule %d",comb);
    }
    if (plist->nr != k)
      gmx_incons("Topology processing, generate nb parameters");
    break;
    
  case F_BHAM:
    /* Buckingham rules */
    for(i=k=0; (i<nr); i++)
      for(j=0; (j<nr); j++,k++) {
	plist->param[k].c[0] = sqrt(atype->nb[i].c[0] * atype->nb[j].c[0]);
	bi = atype->nb[i].c[1];
	bj = atype->nb[j].c[1];
	if ((bi == 0) || (bj == 0))
	  plist->param[k].c[1] = 0;
	else
	  plist->param[k].c[1] = 2.0/(1/bi+1/bj);
	plist->param[k].c[2] = sqrt(atype->nb[i].c[2] * atype->nb[j].c[2]);
      }
    
    break;
  default:
    sprintf(errbuf,"Invalid nonbonded type %s",
	    interaction_function[ftype].longname);
    warning_error(errbuf);
  }
}

void push_at (t_symtab *symtab, t_atomtype *at, t_bond_atomtype *bat,
	      char *line,int nb_funct,
	      t_nbparam ***nbparam,t_nbparam ***pair)
{
  typedef struct {
    char *entry;
    int  ptype;
  } t_xlate;
  t_xlate xl[eptNR] = {
    { "A",   eptAtom }, 
    { "N",   eptNucleus }, 
    { "S",   eptShell }, 
    { "B",   eptBond },
    { "V",   eptVSite }, 
  };
  
  int    nr,i,nfields,j,pt,nfp0=-1;
  int    batype_nr,nread;
  char   type[STRLEN],btype[STRLEN],ptype[STRLEN];
  double m,q;
  double c[MAXFORCEPARAM];
  double radius,vol,surftens;
  char   tmpfield[12][100]; /* Max 12 fields of width 100 */
  char   errbuf[256];
  int    atomnr;
  bool   have_atomic_number;
  bool   have_bonded_type;
  
  /* First assign input line to temporary array */
  nfields=sscanf(line,"%s%s%s%s%s%s%s%s%s%s%s%s",
 		 tmpfield[0],tmpfield[1],tmpfield[2],tmpfield[3],tmpfield[4],tmpfield[5],
		 tmpfield[6],tmpfield[7],tmpfield[8],tmpfield[9],tmpfield[10],tmpfield[11]);

  /* Comments on optional fields in the atomtypes section:
   *
   * The force field format is getting a bit old. For OPLS-AA we needed
   * to add a special bonded atomtype, and for Gerrit Groenhofs QM/MM stuff
   * we also needed the atomic numbers.
   * To avoid making all old or user-generated force fields unusable we
   * have introduced both these quantities as optional columns, and do some
   * acrobatics to check whether they are present or not.
   * This will all look much nicer when we switch to XML... sigh.
   *
   * Field 0 (mandatory) is the nonbonded type name. (string)
   * Field 1 (optional)  is the bonded type (string)
   * Field 2 (optional)  is the atomic number (int)
   * Field 3 (mandatory) is the mass (numerical)
   * Field 4 (mandatory) is the charge (numerical)
   * Field 5 (mandatory) is the particle type (single character)
   * This is followed by a number of nonbonded parameters.
   * 
   * The safest way to identify the format is the particle type field.
   *
   * So, here is what we do:
   *
   * A. Read in the first six fields as strings
   * B. If field 3 (starting from 0) is a single char, we have neither
   *    bonded_type or atomic numbers.
   * C. If field 5 is a single char we have both.
   * D. If field 4 is a single char we check field 1. If this begins with
   *    an alphabetical character we have bonded types, otherwise atomic numbers.
   */  
  
  if(nfields<6)
  {
      too_few();
      return;
  }
  
  if( (strlen(tmpfield[5])==1) && isalpha(tmpfield[5][0]) )
  {
      have_bonded_type   = TRUE;
      have_atomic_number = TRUE;
  }
  else if( (strlen(tmpfield[3])==1) && isalpha(tmpfield[3][0]) )
  {
      have_bonded_type   = FALSE;
      have_atomic_number = FALSE;
  }
  else
  {
      have_bonded_type   = ( isalpha(tmpfield[1][0]) != 0 );
      have_atomic_number = !have_bonded_type; 
  }
  
  /* optional fields */
  surftens = -1;
  vol      =  0;
  radius   =  0;
  atomnr   = -1;
  
  switch (nb_funct) {
      
  case F_LJ:
    nfp0 = 2;

    if( have_atomic_number )
    {
        if ( have_bonded_type )
        {
            nread = sscanf(line,"%s%s%d%lf%lf%s%lf%lf%lf%lf%lf",
                           type,btype,&atomnr,&m,&q,ptype,&c[0],&c[1],
                           &radius,&vol,&surftens);
            if(nread < 8)
            {
                too_few();
                return;
            }            
        }
        else
        {
            /* have_atomic_number && !have_bonded_type */
            nread = sscanf(line,"%s%d%lf%lf%s%lf%lf%lf%lf%lf",
                           type,&atomnr,&m,&q,ptype,&c[0],&c[1],
                           &radius,&vol,&surftens);
            if(nread < 7)
            {
                too_few();
                return;
            }            
        }
    }
    else
    {
        if ( have_bonded_type )
        {
            /* !have_atomic_number && have_bonded_type */
            nread = sscanf(line,"%s%s%lf%lf%s%lf%lf%lf%lf%lf",
                           type,btype,&m,&q,ptype,&c[0],&c[1],
                           &radius,&vol,&surftens);
            if(nread < 7)
            {
                too_few();
                return;
            }            
        }
        else
        {
            /* !have_atomic_number && !have_bonded_type */
            nread = sscanf(line,"%s%lf%lf%s%lf%lf%lf%lf%lf",
                           type,&m,&q,ptype,&c[0],&c[1],
                           &radius,&vol,&surftens);
            if(nread < 6)
            {
                too_few();
                return;
            }            
        }    
    }
    
    if( !have_bonded_type )
    {
        strcpy(btype,type);
    }

    if( !have_atomic_number )
    {
        atomnr = -1;
    }
        
    break;
      
  case F_BHAM:
    nfp0 = 3;

      if( have_atomic_number )
      {
          if ( have_bonded_type )
          {
              nread = sscanf(line,"%s%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf",
                             type,btype,&atomnr,&m,&q,ptype,&c[0],&c[1],&c[2],
                             &radius,&vol,&surftens);
              if(nread < 9)
              {
                  too_few();
                  return;
              }            
          }
          else
          {
              /* have_atomic_number && !have_bonded_type */
              nread = sscanf(line,"%s%d%lf%lf%s%lf%lf%lf%lf%lf%lf",
                             type,&atomnr,&m,&q,ptype,&c[0],&c[1],&c[2],
                             &radius,&vol,&surftens);
              if(nread < 8)
              {
                  too_few();
                  return;
              }            
          }
      }
      else
      {
          if ( have_bonded_type )
          {
              /* !have_atomic_number && have_bonded_type */
              nread = sscanf(line,"%s%s%lf%lf%s%lf%lf%lf%lf%lf%lf",
                             type,btype,&m,&q,ptype,&c[0],&c[1],&c[2],
                             &radius,&vol,&surftens);
              if(nread < 8)
              {
                  too_few();
                  return;
              }            
          }
          else
          {
              /* !have_atomic_number && !have_bonded_type */
              nread = sscanf(line,"%s%lf%lf%s%lf%lf%lf%lf%lf%lf",
                             type,&m,&q,ptype,&c[0],&c[1],&c[2],
                             &radius,&vol,&surftens);
              if(nread < 7)
              {
                  too_few();
                  return;
              }            
          }        
      }
          
      if( !have_bonded_type )
      {
          strcpy(btype,type);
      }
          
      if( !have_atomic_number )
      {
          atomnr = -1;
      }
          
      break;

      default:
          gmx_fatal(FARGS,"Invalid function type %d in push_at %s %d",nb_funct,
                    __FILE__,__LINE__);
  }
  for(j=nfp0; (j<MAXFORCEPARAM); j++)
    c[j]=0.0;
  
  if(strlen(type)==1 && isdigit(type[0])) 
     gmx_fatal(FARGS,"Atom type names can't be single digits.");

  if(strlen(btype)==1 && isdigit(btype[0])) 
     gmx_fatal(FARGS,"Bond atom type names can't be single digits.");

  /* Hack to read old topologies */
  if (strcasecmp(ptype,"D") == 0)
    sprintf(ptype,"V");
  for(j=0; (j<eptNR); j++)
    if (strcasecmp(ptype,xl[j].entry) == 0)
      break;
  if (j == eptNR)
    gmx_fatal(FARGS,"Invalid particle type %s on line %s",
		ptype,line);
  pt=xl[j].ptype;
  if (debug)
    fprintf(debug,"ptype: %s\n",ptype_str[pt]);
  
  nr = bat->nr;
  /* First test if we have a new bond_atomtype */
  for (i=0; ((i<nr) && (strcasecmp(*(bat->atomname[i]),btype) != 0)); i++)
    ;
  if ((i==nr) || (nr==0)) {
    /* New bond_atomtype */
    srenew(bat->atomname,nr+1);
    bat->nr++; 
    bat->atomname[nr] = put_symtab(symtab,btype);
  } 
  batype_nr=i;

  nr = at->nr;
  /* Test if this atomtype overwrites another */
  for (i=0; ((i<nr) && (strcasecmp(*(at->atomname[i]),type) != 0)); i++)
    ;

  if ((i==nr) || (nr==0)) {
    /* New atomtype, get new space for arrays */
    srenew(at->atom,nr+1);
    srenew(at->atomname,nr+1);
    srenew(at->bondatomtype,nr+1);
    srenew(at->nb,nr+1);
    srenew(at->radius,nr+1);
    srenew(at->vol,nr+1);
    srenew(at->surftens,nr+1);
    srenew(at->atomnumber,nr+1);
    at->nr++;
    /* Add space in the non-bonded parameters matrix */
    srenew(*nbparam,at->nr);
    snew((*nbparam)[at->nr-1],at->nr);
    if (pair) {
      srenew(*pair,at->nr);
      snew((*pair)[at->nr-1],at->nr);
    }
  }
  else {
    sprintf(errbuf,"Overriding atomtype %s",type);
    warning(errbuf);
    nr = i;
  } 

  /* fill the arrays */
  at->atomname[nr] = put_symtab(symtab,type);
  at->bondatomtype[nr] = batype_nr;
  at->atom[nr].ptype = pt;
  at->atom[nr].m = m;
  at->atom[nr].q = q;
  at->radius[nr] = radius;
  at->vol[nr] = vol;
  at->surftens[nr] = surftens;
  at->atomnumber[nr] = atomnr;

  for (i=0; (i<MAXFORCEPARAM); i++)
    at->nb[nr].c[i] = c[i];
}

static void push_bondtype(t_params *       bt,
                          t_param *        b,
                          int              nral,
                          int              ftype,
                          bool             bAllowRepeat,
                          char *           line)
{
    int  i,j;
    bool bTest,bFound,bCont,bId;
    int  nr   = bt->nr;
    int  nrfp = NRFP(ftype);
    char  errbuf[256];
    
    /* If bAllowRepeat is TRUE, we allow multiple entries as long as they
        are on directly _adjacent_ lines.
     */
    
    /* First check if our atomtypes are _identical_ (not reversed) to the previous
       entry. If they are not identical we search for earlier duplicates. If they are
       we can skip it, since we already searched for the first line
       in this group.
     */
    
    bFound=FALSE;
    bCont =FALSE;
    
    if(bAllowRepeat && nr>1)
    {
        for (j=0,bCont=TRUE; (j < nral); j++)
        {
            bCont=bCont && (b->a[j] == bt->param[nr-2].a[j]);
        }
    }
    
    /* Search for earlier duplicates if this entry was not a continuation
        from the previous line.
        */
    if(!bCont)
    {
        bFound=FALSE;
        for (i=0; (i < nr); i++) {
            bTest = TRUE;
            for (j=0; (j < nral); j++)
                bTest=(bTest && (b->a[j] == bt->param[i].a[j]));
            if (!bTest) {
                bTest=TRUE;
                for(j=0; (j<nral); j++)
                    bTest=(bTest && (b->a[nral-1-j] == bt->param[i].a[j]));
            }
            if (bTest) {
                if (!bFound) {
                    bId = TRUE;
                    for (j=0; (j < nrfp); j++)
                        bId = bId && (bt->param[i].c[j] == b->c[j]);
                    if (!bId) {
                        sprintf(errbuf,"Overriding %s parameters.%s",
                                interaction_function[ftype].longname,
                                (ftype==F_PDIHS) ? "\nUse dihedraltype 4 to allow several multiplicity terms." : "");
                        warning(errbuf);
                        fprintf(stderr,"  old:");
                        for (j=0; (j < nrfp); j++)
                            fprintf(stderr," %g",bt->param[i].c[j]);
                        fprintf(stderr," \n  new: %s\n\n",line);
                    }
                }
                /* Overwrite it! */
                for (j=0; (j < nrfp); j++)
                    bt->param[i].c[j] = b->c[j];
                bFound=TRUE;
            }
        }
        if (!bFound) {
            /* alloc */
            pr_alloc (2,bt);
            
            /* fill the arrays up and down */
            memcpy(bt->param[bt->nr].c,  b->c,sizeof(b->c));
            memcpy(bt->param[bt->nr].a,  b->a,sizeof(b->a));
            memcpy(bt->param[bt->nr+1].c,b->c,sizeof(b->c));
            
            for (j=0; (j < nral); j++) 
                bt->param[bt->nr+1].a[j] = b->a[nral-1-j];
            
            bt->nr += 2;
        }
    }
}
 

void push_bt(directive d,t_params bt[],int nral,char ***typenames, int ntypes,char *line)
{
  const char *formal[MAXATOMLIST+1] = {
    "%s",
    "%s%s",
    "%s%s%s",
    "%s%s%s%s",
    "%s%s%s%s%s",
    "%s%s%s%s%s%s"
  };
  const char *formnl[MAXATOMLIST+1] = {
    "%*s",
    "%*s%*s",
    "%*s%*s%*s",
    "%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s%*s"
  };
  const char *formlf = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
  int      i,ft,ftype,nn,nrfp,nrfpA,nrfpB;
  char     f1[STRLEN];
  char     alc[MAXATOMLIST+1][20];
  /* One force parameter more, so we can check if we read too many */
  double   c[MAXFORCEPARAM+1];
  t_param  p;
  char  errbuf[256];

  /* Make format string (nral ints+functype) */
  if ((nn=sscanf(line,formal[nral],
		 alc[0],alc[1],alc[2],alc[3],alc[4],alc[5])) != nral+1) {
    sprintf(errbuf,"Not enough atomtypes (%d instead of %d)",nn-1,nral);
    warning_error(errbuf);
    return;
  }
  
  ft    = atoi(alc[nral]);
  ftype = ifunc_index(d,ft);
  nrfp  = NRFP(ftype);
  nrfpA = interaction_function[ftype].nrfpA;
  nrfpB = interaction_function[ftype].nrfpB;
  strcpy(f1,formnl[nral]);
  strcat(f1,formlf);
  if ((nn=sscanf(line,f1,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9],&c[10],&c[11],&c[12])) 
      != nrfp) {
    if (nn == nrfpA) {
      /* Copy the B-state from the A-state */
      for(i=0; (i<nrfpB); i++)
	c[nrfpA+i] = c[i];
    } else {
      if (nn < nrfpA) {
      warning_error("Not enough parameters");
      } else if (nn > nrfpA && nn < nrfp) {
	warning_error("Too many parameters or not enough parameters for topology B");
      } else if (nn > nrfp) {
	warning_error("Too many parameters");
      }
      for(i=nn; (i<nrfp); i++)
	c[i] = 0.0;
    }
  }
  for(i=0; (i<nral); i++)
    p.a[i]=name2index(alc[i],typenames,ntypes); 
  for(i=0; (i<nrfp); i++)
    p.c[i]=c[i];
  push_bondtype (&(bt[ftype]),&p,nral,ftype,FALSE,line);
}


void push_dihedraltype(directive d,t_params bt[],char ***typenames,int ntypes,char *line)
{
  const char *formal[MAXATOMLIST+1] = {
    "%s",
    "%s%s",
    "%s%s%s",
    "%s%s%s%s",
    "%s%s%s%s%s",
    "%s%s%s%s%s%s"
  };
  const char *formnl[MAXATOMLIST+1] = {
    "%*s",
    "%*s%*s",
    "%*s%*s%*s",
    "%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s%*s"
  };
  const char *formlf[MAXFORCEPARAM] = {
    "%lf",
    "%lf%lf",
    "%lf%lf%lf",
    "%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
  };
  int      i,ft,ftype,nn,nrfp,nrfpA,nrfpB,nral;
  char     f1[STRLEN];
  char     alc[MAXATOMLIST+1][20];
  double   c[MAXFORCEPARAM];
  t_param  p;
  bool     bAllowRepeat;
  char  errbuf[256];

  /* This routine accepts dihedraltypes defined from either 2 or 4 atoms.
   *
   * We first check for 2 atoms with the 3th column being an integer 
   * defining the type. If this isn't the case, we try it with 4 atoms
   * and the 5th column defining the dihedral type.
   */
  nn=sscanf(line,formal[4],alc[0],alc[1],alc[2],alc[3],alc[4]);
  if(nn>=3 && strlen(alc[2])==1 && isdigit(alc[2][0])) {
    nral=2;
    ft    = atoi(alc[nral]);
    /* Move atom types around a bit and use 'X' for wildcard atoms
     * to create a 4-atom dihedral definition with arbitrary atoms in
     * position 1 and 4.
     */
    if(alc[2][0]=='2') {
      /* improper - the two atomtypes are 1,4. Use wildcards for 2,3 */
      strcpy(alc[3],alc[1]);
      sprintf(alc[2],"X");
      sprintf(alc[1],"X");
      /* alc[0] stays put */
    } else {
      /* proper - the two atomtypes are 2,3. Use wildcards for 1,4 */
      sprintf(alc[3],"X");
      strcpy(alc[2],alc[1]);
      strcpy(alc[1],alc[0]);
      sprintf(alc[0],"X");
    }
  } else if(nn==5 && strlen(alc[4])==1 && isdigit(alc[4][0])) {
    nral=4;
    ft    = atoi(alc[nral]);
  } else {
    sprintf(errbuf,"Incorrect number of atomtypes for dihedral (%d instead of 2 or 4)",nn-1);
    warning_error(errbuf);
    return;
  }
  
  
  if(ft == 9)
  {
      /* Previously, we have always overwritten parameters if e.g. a torsion
      with the same atomtypes occurs on multiple lines. However, CHARMM and
      some other force fields specify multiple dihedrals over some bonds, 
      including cosines with multiplicity 6 and somethimes even higher.
      Thus, they cannot be represented with Ryckaert-Bellemans terms.
      To add support for these force fields, Dihedral type 4 is identical to
      normal proper dihedrals, but repeated entries are allowed. 
      */
      bAllowRepeat = TRUE;
      ft = 1;
  }
  else
  {
      bAllowRepeat = FALSE;
  }  
  
  ftype = ifunc_index(d,ft);
  nrfp  = NRFP(ftype);
  nrfpA = interaction_function[ftype].nrfpA;
  nrfpB = interaction_function[ftype].nrfpB;
 
  strcpy(f1,formnl[nral]);
  strcat(f1,formlf[nrfp-1]);
  
  /* Check number of parameters given */
  if ((nn=sscanf(line,f1,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9],&c[10],&c[11]))
      != nrfp) {
    if (nn == nrfpA) {
      /* Copy the B-state from the A-state */
      for(i=0; (i<nrfpB); i++)
	c[nrfpA+i] = c[i];
    } else {
      if (nn < nrfpA) {
      warning_error("Not enough parameters");
      } else if (nn > nrfpA && nn < nrfp) {
	warning_error("Too many parameters or not enough parameters for topology B");
      } else if (nn > nrfp) {
	warning_error("Too many parameters");
      }
      for(i=nn; (i<nrfp); i++)
	c[i] = 0.0;
    }
  }
  
  for(i=0; (i<4); i++) {
    if(!strcmp(alc[i],"X"))
      p.a[i]=-1;
    else
    p.a[i]=name2index(alc[i],typenames,ntypes); 
  }
  for(i=0; (i<nrfp); i++)
    p.c[i]=c[i];
  /* Always use 4 atoms here, since we created two wildcard atoms
   * if there wasn't of them 4 already.
   */
  push_bondtype (&(bt[ftype]),&p,4,ftype,bAllowRepeat,line);
}


void push_nbt(directive d,t_nbparam **nbt,t_atomtype *atype,
	      char *pline,int nb_funct)
{
  /* swap the atoms */
  const char *form2="%*s%*s%*s%lf%lf";
  const char *form3="%*s%*s%*s%lf%lf%lf";
  const char *form4="%*s%*s%*s%lf%lf%lf%lf";
  char    a0[80],a1[80];
  int     i,f,n,ftype,atnr,nrfp;
  double  c[4];
  real    cr[4],sig6;
  atom_id ai,aj;
  t_nbparam *nbp;
  bool    bId;
  char  errbuf[256];

  if (sscanf (pline,"%s%s%d",a0,a1,&f) != 3) {
    too_few();
    return;
  }
  
  ftype=ifunc_index(d,f);
  
  if (ftype != nb_funct) {
    sprintf(errbuf,"Trying to add %s while the default nonbond type is %s",
	    interaction_function[ftype].longname,
	    interaction_function[nb_funct].longname);
    warning_error(errbuf);
    return;
  }
  
  /* Get the force parameters */
  nrfp = NRFP(ftype);
  if (ftype == F_LJ14) {
    n = sscanf(pline,form4,&c[0],&c[1],&c[2],&c[3]);
    if (n < 2) {
      too_few();
      return;
    }
    /* When the B topology parameters are not set,
     * copy them from topology A
     */ 
    for(i=n; i<nrfp; i++)
      c[i] = c[i-2];
  }
  else if (nrfp == 2) {
    if (sscanf(pline,form2,&c[0],&c[1]) != 2) {
      too_few();
      return;
    }
  }
  else if (nrfp == 3) {
    if (sscanf(pline,form3,&c[0],&c[1],&c[2]) != 3) {
      too_few();
      return;
    }
  }
  else {
    gmx_fatal(FARGS,"Number of force parameters for nonbonded interactions is %d"
		" in file %s, line %d",nrfp,__FILE__,__LINE__);
  }
  for(i=0; (i<nrfp); i++)
    cr[i] = c[i];

  /* Put the parameters in the matrix */
  ai = at2type (a0,atype);
  aj = at2type (a1,atype);
  nbp = &(nbt[max(ai,aj)][min(ai,aj)]);
  
  if (nbp->bSet) {
    bId = TRUE;
    for (i=0; i<nrfp; i++)
      bId = bId && (nbp->c[i] == cr[i]);
    if (!bId) {
      sprintf(errbuf,"Overriding non-bonded parameters,");
      warning(errbuf);
      fprintf(stderr,"  old:");
      for (i=0; i<nrfp; i++)
	fprintf(stderr," %g",nbp->c[i]);
      fprintf(stderr," new\n%s\n",pline);
    }
  }
  nbp->bSet = TRUE;
  for (i=0; i<nrfp; i++)
    nbp->c[i] = cr[i];
}

static void push_atom_now(t_symtab *symtab,t_atoms *at,int atomnr,
			  int type,char *ctype,int ptype,
			  int resnumber,int cgnumber,
			  char *resname,char *name,real m0,real q0,
			  int typeB,char *ctypeB,real mB,real qB)
{
  int j,resnr_diff=-1;
  int nr = at->nr;

  if (((nr==0) && (atomnr != 1)) || (nr && (atomnr != at->nr+1)))
    gmx_fatal(FARGS,"Atoms in the .top are not numbered consecutively from 1 (rather, atomnr = %d, while at->nr = %d)",atomnr,at->nr);
  if (nr)
    resnr_diff = resnumber - (at->atom[at->nr-1].resnr+1);
  if (((nr==0) && (resnumber != 1)) || 
      (nr && (resnr_diff != 0) && (resnr_diff != 1))) 
    gmx_fatal(FARGS,"Residue numbers in the .top are not numbered consecutively from 1 (rather, resnumber = %d and resnr_diff = %d)",resnumber,resnr_diff);
  
  /* New atom instance
   * get new space for arrays 
   */
  srenew(at->atom,nr+1);
  srenew(at->atomname,nr+1);
  srenew(at->atomtype,nr+1);
  srenew(at->atomtypeB,nr+1);

  if (resnumber > at->nres) {
    at->nres=resnumber;
    srenew(at->resname,resnumber);
    at->resname[resnumber-1] = put_symtab(symtab,resname);
  }
  /* fill the list */
  at->atom[nr].type  = type;
  at->atom[nr].ptype = ptype;
  at->atom[nr].q     = q0;
  at->atom[nr].m     = m0;
  at->atom[nr].typeB = typeB;
  at->atom[nr].qB    = qB;
  at->atom[nr].mB    = mB;
  
  for(j=0; (j<egcNR); j++)
    at->atom[nr].grpnr[j] = -1;
  at->atom[nr].resnr = resnumber-1;
  at->atomname[nr] = put_symtab(symtab,name);
  at->atomtype[nr] = put_symtab(symtab,ctype);
  at->atomtypeB[nr] = put_symtab(symtab,ctypeB);
  at->nr++;
}

void push_cg(t_block *block, int *lastindex, int index, int a)
{
  if (debug)
    fprintf (debug,"Index %d, Atom %d\n",index,a);

  if (((block->nr) && (*lastindex != index)) || (!block->nr)) {
    /* add a new block */
    block->nr++;
    srenew(block->index,block->nr+1);
  }
  srenew(block->a,block->nra+1);
  block->a[block->nra++]=a;
  block->index[block->nr]=block->nra;
  *lastindex = index;
}

void push_atom(t_symtab *symtab,t_block *cgs,
	       t_atoms *at,t_atomtype *atype,char *line,int *lastcg)
{
  int 		nr,ptype;
  int 		resnumber,cgnumber,atomnr,type,typeB,nscan;
  char 		id[STRLEN],ctype[STRLEN],ctypeB[STRLEN],
       		resname[STRLEN],name[STRLEN],check[STRLEN];
  double        m,q,mb,qb;
  real          m0,q0,mB,qB;

  /* Make a shortcut for writing in this molecule  */
  nr = at->nr;

  /* Fixed parameters */
  if (sscanf(line,"%s%s%d%s%s%d",
	     id,ctype,&resnumber,resname,name,&cgnumber) != 6) {
    too_few();
    return;
  }
  sscanf(id,"%d",&atomnr);
  type  = at2type(ctype,atype);
  ptype = atype->atom[type].ptype;

  /* Set default from type */  
  q0    = atype->atom[type].q;
  m0    = atype->atom[type].m;
  typeB = type;
  qB    = q0;
  mB    = m0;
  
  /* Optional parameters */
  nscan=sscanf(line,"%*s%*s%*s%*s%*s%*s%lf%lf%s%lf%lf%s",
	       &q,&m,ctypeB,&qb,&mb,check);
  
  /* Nasty switch that falls thru all the way down! */
  if (nscan > 0) {
    q0 = qB = q;
    if (nscan > 1) {
      m0 = mB = m;
      if (nscan > 2) {
	typeB=at2type(ctypeB,atype);
	qB = atype->atom[typeB].q;
	mB = atype->atom[typeB].m;
	if (nscan > 3) {
	  qB = qb;
	  if (nscan > 4) {
	    mB = mb;
	    if (nscan > 5)
	      warning_error("Too many parameters");
	  }
	}
      }
    }
  }
  if (debug) 
    fprintf(debug,"mB=%g, qB=%g, typeB=%d\n",mB,qB,typeB);
  
  push_cg(cgs,lastcg,cgnumber,nr);

  push_atom_now(symtab,at,atomnr,type,ctype,ptype,resnumber,cgnumber,
		resname,name,m0,q0,typeB,
		typeB==type ? ctype : ctypeB,mB,qB);
}

void push_molt(t_symtab *symtab,int *nmol,t_molinfo **mol,char *line)
{
  char type[STRLEN];
  int nrexcl,i;
  t_molinfo *newmol;

  if ((sscanf(line,"%s%d",type,&nrexcl)) != 2) {
    warning_error("Expected a molecule type name and nrexcl");
  }
  
  /* Test if this atomtype overwrites another */
  i = 0;
  while (i < *nmol) {
    if (strcasecmp(*((*mol)[i].name),type) == 0)
      gmx_fatal(FARGS,"moleculetype %s is redefined",type);
    i++;
  }
  
  (*nmol)++;
  srenew(*mol,*nmol);
  newmol = &((*mol)[*nmol-1]);
  init_molinfo(newmol);
  
  /* Fill in the values */
  newmol->name     = put_symtab(symtab,type);
  newmol->nrexcl   = nrexcl;
  newmol->excl_set = FALSE;
}

static bool default_nb_params(int ftype,t_params bt[],t_atoms *at,
			      t_param *p,bool bB,bool bGenPairs)
{
  int      i,j,ti,tj,ntype;
  bool     bFound;
  t_param  *pi=NULL;
  int      nr   = bt[ftype].nr;
  int      nral = NRAL(ftype);
  int      nrfp = interaction_function[ftype].nrfpA;
  int      nrfpB= interaction_function[ftype].nrfpB;

  if ((!bB && nrfp == 0) || (bB && nrfpB == 0))
    return TRUE;

  bFound=FALSE;
  if(bGenPairs) {
    /* First test the generated-pair position to save
     * time when we have 1000*1000 entries for e.g. OPLS...
     */
    ntype=sqrt(nr);
    if(bB) {
      ti=at->atom[p->a[0]].typeB;
      tj=at->atom[p->a[1]].typeB;
    } else {
      ti=at->atom[p->a[0]].type;
      tj=at->atom[p->a[1]].type;
    } 
    pi=&(bt[ftype].param[ntype*ti+tj]);
    bFound=((ti == pi->a[0]) && (tj == pi->a[1]));
  }

  /* Search explicitly if we didnt find it */
  if(!bFound) {
    for (i=0; ((i < nr) && !bFound); i++) {
      pi=&(bt[ftype].param[i]);
      if (bB)
	for (j=0; ((j < nral) && 
		   (at->atom[p->a[j]].typeB == pi->a[j])); j++);
      else
	for (j=0; ((j < nral) && 
		   (at->atom[p->a[j]].type == pi->a[j])); j++);
      bFound=(j==nral);
    }
  }
  
  if (bFound) {
    if (bB) {
      if (nrfp+nrfpB > MAXFORCEPARAM) {
	gmx_incons("Too many force parameters");
      }
      for (j=0; (j < nrfpB); j++)
	p->c[nrfp+j] = pi->c[j];
    }
    else
      for (j=0; (j < nrfp); j++)
	p->c[j] = pi->c[j];
  }
  else {
    for (j=0; (j < nrfp); j++)
      p->c[j] = 0.0;
  }
  return bFound;
}


static bool default_params(int ftype,t_params bt[],
			   t_atoms *at,t_atomtype *atype,
			   t_param *p,bool bB,
			   t_param **param_def)
{
  int      i,j;
  bool     bFound;
  t_param  *pi=NULL;
  int      nr   = bt[ftype].nr;
  int      nral = NRAL(ftype);
  int      nrfpA= interaction_function[ftype].nrfpA;
  int      nrfpB= interaction_function[ftype].nrfpB;
  int     *batype=atype->bondatomtype;

  if ((!bB && nrfpA == 0) || (bB && nrfpB == 0))
    return TRUE;


  /* We allow wildcards now. The first type (with or without wildcards) that
   * fits is used, so you should probably put the wildcarded bondtypes
   * at the end of each section.
   */
  bFound=FALSE;
  /* OPLS uses 1000s of dihedraltypes, so in order to speed up the scanning we have a 
   * special case for this. Check for B state outside loop to speed it up.
   */
  if( ftype==F_PDIHS || ftype == F_RBDIHS || ftype == F_IDIHS) {
    if (bB) {
      for (i=0; ((i < nr) && !bFound); i++) {
	pi=&(bt[ftype].param[i]);
	bFound=( ((pi->AI==-1) || (batype[at->atom[p->AI].typeB]==pi->AI)) &&
		 ((pi->AJ==-1) || (batype[at->atom[p->AJ].typeB]==pi->AJ)) &&
		 ((pi->AK==-1) || (batype[at->atom[p->AK].typeB]==pi->AK)) &&
		 ((pi->AL==-1) || (batype[at->atom[p->AL].typeB]==pi->AL)) );
      }
    } else {
      for (i=0; ((i < nr) && !bFound); i++) {
	pi=&(bt[ftype].param[i]);
	bFound=( ((pi->AI==-1) || (batype[at->atom[p->AI].type]==pi->AI)) &&
		 ((pi->AJ==-1) || (batype[at->atom[p->AJ].type]==pi->AJ)) &&
		 ((pi->AK==-1) || (batype[at->atom[p->AK].type]==pi->AK)) &&
		 ((pi->AL==-1) || (batype[at->atom[p->AL].type]==pi->AL)) );
      }
    }
  } else { /* Not a dihedral */
    for (i=0; ((i < nr) && !bFound); i++) {
      pi=&(bt[ftype].param[i]);
      if (bB)
	for (j=0; ((j < nral) && 
		   (atype->bondatomtype[at->atom[p->a[j]].typeB] == pi->a[j])); j++);
      else
	for (j=0; ((j < nral) && 
		   (atype->bondatomtype[at->atom[p->a[j]].type] == pi->a[j])); j++);
      bFound=(j==nral);
    }
  }

  *param_def = pi;

  return bFound;
}



void push_bondnow(t_params *bond, t_param *b)
{
  int j;
  
  if (debug) {
    fprintf(debug,"push_bondnow: nr = %d\n",bond->nr);
    fflush(debug);
  }
  /* allocate one position extra */
  pr_alloc (1,bond);

  /* fill the arrays */
  for (j=0; (j < MAXFORCEPARAM); j++) 
    bond->param[bond->nr].c[j]   = b->c[j];
  for (j=0; (j < MAXATOMLIST); j++) 
    bond->param[bond->nr].a[j]   = b->a[j];
  
  bond->nr++;
}

void push_bond(directive d,t_params bondtype[],t_params bond[],
	       t_atoms *at,t_atomtype *atype,char *line,
	       bool bBonded,bool bGenPairs,bool bZero,bool *bWarn_copy_A_B)
{
  const char *aaformat[MAXATOMLIST]= {
    "%d%d",
    "%d%d%d",
    "%d%d%d%d",
    "%d%d%d%d%d",
    "%d%d%d%d%d%d"
  };
  const char *asformat[MAXATOMLIST]= {
    "%*s%*s",
    "%*s%*s%*s",
    "%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s%*s"
  };
  const char *ccformat="%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
  int      nr,i,j,nral,nread,ftype;
  char     format[STRLEN];
  /* One force parameter more, so we can check if we read too many */
  double   cc[MAXFORCEPARAM+1];
  int      aa[MAXATOMLIST+1];
  t_param  param,paramB,*param_def;
  bool     bFoundA,bFoundB,bDef,bPert,bSwapParity=FALSE;
  char  errbuf[256];

  ftype = ifunc_index(d,1);
  nral  = NRAL(ftype);
  for(j=0; j<MAXATOMLIST; j++)
    aa[j]=NOTSET;
  bDef = (NRFP(ftype) > 0);
  
  nread = sscanf(line,aaformat[nral-1],
		 &aa[0],&aa[1],&aa[2],&aa[3],&aa[4],&aa[5]);
  if (nread < nral) {
    too_few();
    return;
  } else if (nread == nral) 
    ftype = ifunc_index(d,1);
  else {
    /* this is a hack to allow for virtual sites with swapped parity */
    bSwapParity = (aa[nral]<0);
    if (bSwapParity)
      aa[nral] = -aa[nral];
    ftype = ifunc_index(d,aa[nral]);
    if (bSwapParity)
      switch(ftype) {
      case F_VSITE3FAD:
      case F_VSITE3OUT:
	break;
      default:
	gmx_fatal(FARGS,"Negative function types only allowed for %s and %s",
		    interaction_function[F_VSITE3FAD].longname,
		    interaction_function[F_VSITE3OUT].longname);
      }
  }
  
  /* Check for double atoms and atoms out of bounds */
  for(i=0; (i<nral); i++) {
    if ( aa[i] < 1 || aa[i] > at->nr )
      gmx_fatal(FARGS,"[ file %s, line %d ]:\n"
		"Atom index (%d) in %s out of bounds (1-%d).\n"
		"This probably means that you have inserted topology section \"%s\"\n"
		"in a part belonging to a different molecule than you intended to.\n" 
		"In that case move the \"%s\" section to the right molecule.",
		get_warning_file(),get_warning_line(),
		aa[i],dir2str(d),at->nr,dir2str(d),dir2str(d));
    for(j=i+1; (j<nral); j++)
      if (aa[i] == aa[j]) {
	sprintf(errbuf,"Duplicate atom index (%d) in %s",aa[i],dir2str(d));
	warning(errbuf);
      }
  }
  if (ftype == F_SETTLE)
    if (aa[0]+2 > at->nr)
      gmx_fatal(FARGS,"[ file %s, line %d ]:\n"
		  "             Atom index (%d) in %s out of bounds (1-%d)\n"
		  "             Settle works on atoms %d, %d and %d",
		  get_warning_file(),get_warning_line(),
		  aa[0],dir2str(d),at->nr,
		  aa[0],aa[0]+1,aa[0]+2);
  
  /* default force parameters  */
  for(j=0; (j<MAXATOMLIST); j++)
    param.a[j] = aa[j]-1;
  for(j=0; (j<MAXFORCEPARAM); j++)
    param.c[j] = 0.0;
  
  /* Get force params for normal and free energy perturbation
   * studies, as determined by types!
   */
  if(bBonded) {
    bFoundA = default_params(ftype,bondtype,at,atype,&param,FALSE,&param_def);
    if (bFoundA) {
      /* Copy the A-state and B-state default parameters */
      for(j=0; (j<NRFPA(ftype)+NRFPB(ftype)); j++)
	param.c[j] = param_def->c[j];
    }
    bFoundB = default_params(ftype,bondtype,at,atype,&param,TRUE,&param_def);
    if (bFoundB) {
      /* Copy only the B-state default parameters */
      for(j=NRFPA(ftype); (j<NRFP(ftype)); j++)
	param.c[j] = param_def->c[j];
    }
  } else {
    bFoundA = default_nb_params(ftype,bondtype,at,&param,FALSE,bGenPairs);
    bFoundB = default_nb_params(ftype,bondtype,at,&param,TRUE,bGenPairs);
  }
  
  if (nread > nral) {  
    strcpy(format,asformat[nral-1]);
    strcat(format,ccformat);
    
    nread = sscanf(line,format,&cc[0],&cc[1],&cc[2],&cc[3],&cc[4],&cc[5],
		   &cc[6],&cc[7],&cc[8],&cc[9],&cc[10],&cc[11],&cc[12]);

    if ((nread == NRFPA(ftype)) && (NRFPB(ftype) != 0))
      {
	/* We only have to issue a warning if these atoms are perturbed! */
	bPert = FALSE;
	for(j=0; (j<nral); j++)
	  bPert = bPert || PERTURBED(at->atom[param.a[j]]);

	if (bPert && *bWarn_copy_A_B)
	  {
	    sprintf(errbuf,
		    "Some parameters for bonded interaction involving perturbed atoms are specified explicitly in state A, but not B - copying A to B");
	    warning(errbuf);
	    *bWarn_copy_A_B = FALSE;
	  }
	
	/* If only the A parameters were specified, copy them to the B state */
	/* The B-state parameters correspond to the first nrfpB
	 * A-state parameters.
	 */
	for(j=0; (j<NRFPB(ftype)); j++)
	  cc[nread++] = cc[j];
      }
    
    /* If nread was 0 or EOF, no parameters were read => use defaults.
     * If nread was nrfpA we copied above so nread=nrfp.
     * If nread was nrfp we are cool.
     * Anything else is an error!
     */	
    if ((nread != 0) && (nread != EOF) && (nread != NRFP(ftype)))
      {
	gmx_fatal(FARGS,"Incorrect number of parameters - found %d, expected %d or %d for %s.",
		  nread,NRFPA(ftype),NRFP(ftype),
		  interaction_function[ftype].longname);	
      }
    
    for(j=0; (j<nread); j++)
      param.c[j]=cc[j];
      
    /* Check whether we have to use the defaults */
    if (nread == NRFP(ftype))
      bDef = FALSE;
  } 
  else
  {
    nread = 0;
  }
   	
  /* nread now holds the number of force parameters read! */
  
  if (bDef) {
    /* Use defaults */

    if (nread > 0 && nread < NRFPA(ftype)) {
      /* Issue an error, do not use defaults */
      sprintf(errbuf,"Not enough parameters, there should be at least %d (or 0 for defaults)",NRFPA(ftype));
      warning_error(errbuf);
    } else {
      if (nread == 0 || nread == EOF) {
	if (!bFoundA) {
	  if (interaction_function[ftype].flags & IF_VSITE) {
	    /* set them to NOTSET, will be calculated later */
	    for(j=0; (j<MAXFORCEPARAM); j++)
	      param.c[j] = NOTSET;
	    if (bSwapParity)
	      param.C1 = -1; /* flag to swap parity of vsite construction */
	  } else {
	    if (bZero) {
	      fprintf(stderr,"NOTE: No default %s types, using zeroes\n",
		      interaction_function[ftype].longname);
	    } else {
	      sprintf(errbuf,"No default %s types",
		      interaction_function[ftype].longname);
	      warning_error(errbuf);
	    }
	  }
	} else
	  if (bSwapParity)
	    switch(ftype) {
	    case F_VSITE3FAD:
	      param.C0 = 360 - param.C0;
	      break;
	    case F_VSITE3OUT:
	      param.C2 = -param.C2;
	      break;
	    }
      }
      if (!bFoundB) {
	/* We only have to issue a warning if these atoms are perturbed! */
	bPert = FALSE;
	for(j=0; (j<nral); j++)
	  bPert = bPert || PERTURBED(at->atom[param.a[j]]);
	
	if (bPert) {
	  sprintf(errbuf,"No default %s types for perturbed atoms, "
		  "using normal values",interaction_function[ftype].longname);
	  warning(errbuf);
	}
      }
    }
  }
  
  if ((ftype==F_PDIHS || ftype==F_ANGRES || ftype==F_ANGRESZ)
      && param.c[5]!=param.c[2])
    gmx_fatal(FARGS,"[ file %s, line %d ]:\n"
		"             %s multiplicity can not be perturbed %f!=%f",
		get_warning_file(),get_warning_line(),
		interaction_function[ftype].longname,
		param.c[2],param.c[5]);

  /* Dont add R-B dihedrals where all parameters are zero (no interaction) */
  if (ftype==F_RBDIHS) {
    nr=0;
    for(i=0;i<NRFP(ftype);i++) {
      if(param.c[i]!=0)
	nr++;
    }
    if(nr==0)
      return;
  }

  /* Put the values in the appropriate arrays */
  push_bondnow (&bond[ftype],&param);
}




void push_mol(int nrmols,t_molinfo mols[],char *pline,int *whichmol,
		  int *nrcopies)
{
  int  i,copies;
  char type[STRLEN];

  *nrcopies=0;
  if (sscanf(pline,"%s%d",type,&copies) != 2) {
    too_few();
    return;
  }
  
  /* search moleculename */
  for (i=0; ((i<nrmols) && strcasecmp(type,*(mols[i].name))); i++)
    ;

  if (i<nrmols) {
    *nrcopies        = copies;
    *whichmol        = i;
  } else
    gmx_fatal(FARGS,"No such moleculetype %s",type);
}

void init_block2(t_block2 *b2, int natom)
{
  int i;

  b2->nr=natom;
  snew(b2->nra,b2->nr);
  snew(b2->a,b2->nr);
  for(i=0; (i<b2->nr); i++)
    b2->a[i]=NULL;
}

void done_block2(t_block2 *b2)
{
  int i;
  
  if (b2->nr) {
    for(i=0; (i<b2->nr); i++)
      sfree(b2->a[i]);
    sfree(b2->a);
    sfree(b2->nra);
    b2->nr=0;
  }
}

void push_excl(char *line, t_block2 *b2)
{
  int  i,j;
  int  n;
  char base[STRLEN],format[STRLEN];

  if (sscanf(line,"%d",&i)==0)
    return;
    
  if ((1 <= i) && (i <= b2->nr))
    i--;
  else {
    if (debug)
      fprintf(debug,"Unbound atom %d\n",i-1);
    return;
  }
  strcpy(base,"%*d");
  do {
    strcpy(format,base);
    strcat(format,"%d");
    n=sscanf(line,format,&j);
    if (n == 1) {
      if ((1 <= j) && (j <= b2->nr)) {
	j--;
	srenew(b2->a[i],++(b2->nra[i]));
	b2->a[i][b2->nra[i]-1]=j;
	/* also add the reverse exclusion! */
	srenew(b2->a[j],++(b2->nra[j]));
	b2->a[j][b2->nra[j]-1]=i;
	strcat(base,"%*d");
      }
      else 
	gmx_fatal(FARGS,"Invalid Atomnr j: %d, b2->nr: %d\n",j,b2->nr);
    }
  } while (n == 1);
}

void b_to_b2(t_block *b, t_block2 *b2)
{
  int     i;
  atom_id j,a;

  for(i=0; (i<b->nr); i++)
    for(j=b->index[i]; (j<b->index[i+1]); j++) {
      a=b->a[j];
      srenew(b2->a[i],++b2->nra[i]);
      b2->a[i][b2->nra[i]-1]=a;
    }
}

void b2_to_b(t_block2 *b2, t_block *b)
{
  int     i,nra;
  atom_id j;

  nra=0;
  for(i=0; (i<b2->nr); i++) {
    b->index[i]=nra;
    for(j=0; (j<b2->nra[i]); j++)
      b->a[nra+j]=b2->a[i][j];
    nra+=b2->nra[i];
  }
  /* terminate list */
  b->index[i]=nra;
}

static int icomp(const void *v1, const void *v2)
{
  return (*((atom_id *) v1))-(*((atom_id *) v2));
}

void merge_excl(t_block *excl, t_block2 *b2)
{
  int     i,k;
  atom_id j;
  int     nra;

  if (!b2->nr)
    return;
  else if (b2->nr != excl->nr) {
    gmx_fatal(FARGS,"DEATH HORROR: b2->nr = %d, while excl->nr = %d",
		b2->nr,excl->nr);
  }
  else if (debug)
    fprintf(debug,"Entering merge_excl\n");

  /* First copy all entries from excl to b2 */
  b_to_b2(excl,b2);

  /* Count and sort the exclusions */
  nra=0;
  for(i=0; (i<b2->nr); i++) {
    if (b2->nra[i] > 0) {
      /* remove double entries */
      qsort(b2->a[i],(size_t)b2->nra[i],(size_t)sizeof(b2->a[i][0]),icomp);
      k=1;
      for(j=1; (j<b2->nra[i]); j++)
	if (b2->a[i][j]!=b2->a[i][k-1]) {
	  b2->a[i][k]=b2->a[i][j];
	  k++;
	}
      b2->nra[i]=k;
      nra+=b2->nra[i];
    }
  }
  excl->nra=nra;
  srenew(excl->a,excl->nra);

  b2_to_b(b2,excl);
}

