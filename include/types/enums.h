/*
 * $Id: enums.h,v 1.49.2.3 2008/02/29 07:02:42 spoel Exp $
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

/* note: these enums should correspond to the names in gmxlib/names.c */

enum {
  ebCGS,ebMOLS,ebSBLOCKS,ebNR
};

enum {
  epbcXYZ, epbcNONE, epbcFULL, epbcNR
};

enum {
  etcNO, etcBERENDSEN, etcNOSEHOOVER, etcYES, etcANDERSEN, etcANDERSENINTERVAL, etcNR
}; /* yes is an alias for berendsen */

enum {
  epcNO, epcBERENDSEN, epcPARRINELLORAHMAN, epcISOTROPIC, epcNR
}; /* isotropic is an alias for berendsen */

enum {
  epctISOTROPIC, epctSEMIISOTROPIC, epctANISOTROPIC,
  epctSURFACETENSION, epctNR
};

enum {
  eelCUT,     eelRF,     eelGRF,   eelPME,  eelEWALD,  eelPPPM, 
  eelPOISSON, eelSWITCH, eelSHIFT, eelUSER, eelGB, eelRF_NEC, eelENCADSHIFT, 
  eelPMEUSER, eelNR
};

/* Ewald geometry */
enum { 
  eewg3D, eewg3DC, eewgNR
};

#define EEL_RF(e) ((e == eelRF) || (e == eelGRF) || (e == eelRF_NEC))

#define EEL_FULL(e) ((e == eelPPPM) || (e == eelPOISSON) || (e ==  eelPME) || (e ==  eelPMEUSER) || (e == eelEWALD))

enum {
  evdwCUT, evdwSWITCH, evdwSHIFT, evdwUSER, evdwENCADSHIFT, evdwNR
};

enum { 
  ensGRID, ensSIMPLE, ensNR
};

enum {
  eiMD, eiSteep, eiCG, eiBD, eiSD, eiNM, eiLBFGS, eiTPI, eiNR
};

#define EI_DYNAMICS(e) ((e == eiMD) || (e == eiSD) || (e == eiBD))
#define EI_ENERGY_MINIMIZATION(e) ((e == eiSteep) || (e == eiCG) || (e == eiLBFGS))

enum {
  estLINCS, estSHAKE, estNR
};

enum {
  edrNone, edrSimple, edrEnsemble, edrNR
};

enum {
  edrwConservative, edrwEqual, edrwNR
};

/* Combination rule things */
enum { 
  eCOMB_NONE, eCOMB_GEOMETRIC, eCOMB_ARITHMETIC, eCOMB_GEOM_SIG_EPS, eCOMB_NR 
};

/* NBF selection */
enum { 
  eNBF_NONE, eNBF_LJ, eNBF_BHAM, eNBF_NR 
};

/* FEP selection */
enum {
  efepNO, efepYES, efepNR
};

/* Solvent model */
enum {
  esolNO, esolSPC, esolTIP4P, esolNR
};

/* Neighborlist type */
enum {
  enlistATOM, enlistWATER, enlistWATERWATER, enlistNR
};

/* Dispersion correction */
enum {
  edispcNO, edispcEnerPres, edispcEner, edispcAllEnerPres, edispcAllEner, edispcNR
}; 

/* Shell types, for completion stuff */
enum {
  eshellCSH, eshellBASH, eshellZSH, eshellNR
}; 

/* Center of mass motion selection */
enum { 
  ecmLINEAR, ecmANGULAR, ecmNO, ecmNR 
};

/* New version of simulated annealing */
enum { 
  eannNO, eannSINGLE, eannPERIODIC, eannNR 
};

/* Algorithms for calculating GB radii */
enum { 
  egbSTILL, egbKARPLUS, egbNR 
};

/* Implicit solvent algorithms */
enum { 
  eisNO, eisLCPO, eisNR 
};

enum {
  eQMmethodAM1, eQMmethodPM3, eQMmethodRHF, 
  eQMmethodUHF, eQMmethodDFT, eQMmethodB3LYP, eQMmethodMP2, eQMmethodCASSCF, eQMmethodB3LYPLAN,
  eQMmethodDIRECT, eQMmethodNR
};

enum {
  eQMbasisSTO3G, eQMbasisSTO3G2, eQMbasis321G, 
  eQMbasis321Gp, eQMbasis321dGp, eQMbasis621G,
  eQMbasis631G, eQMbasis631Gp, eQMbasis631dGp, 
  eQMbasis6311G, eQMbasisNR
};

enum {
  eQMMMschemenormal,eQMMMschemeoniom,eQMMMschemeNR
};

