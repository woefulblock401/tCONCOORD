/*
 * $Id: gmx_fft_fftw3.c,v 1.1.2.6 2008/02/29 07:02:51 spoel Exp $
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

#include <errno.h>
#include <stdlib.h>

#include <fftw3.h>

#include "gmx_fft.h"
#include "gmx_fatal.h"


#ifdef GMX_DOUBLE
#define FFTWPREFIX(name) fftw_ ## name
#else
#define FFTWPREFIX(name) fftwf_ ## name
#endif

struct gmx_fft
{
    /* Three alternatives (unaligned/aligned, out-of-place/in-place, forward/backward)
     * results in 8 different FFTW plans. Keep track of them with 3 array indices:
     * first index:   0=unaligned, 1=aligned
     * second index:  0=out-of-place, 1=in-place
     * third index:   0=backward, 1=forward
     */
    FFTWPREFIX(plan)         plan[2][2][2];
    /* Catch user mistakes */
    int                      real_transform;
    int                      ndim;
};



int
gmx_fft_init_1d(gmx_fft_t *   pfft,
                int           nx) 
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1,*p2,*up1,*up2;
    size_t                 pc;
    int                    i,j,k;
    
    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx+2));
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx+2));
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    /* make unaligned pointers. 
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     *
     * Do NOT use size_t, since that is 32 bytes on some NEC UX machines.
     */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (FFTWPREFIX(complex) *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (FFTWPREFIX(complex) *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_1d)(nx,up1,up2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_1d)(nx,up1,up2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_1d)(nx,up1,up1,FFTW_BACKWARD,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_1d)(nx,up1,up1,FFTW_FORWARD,FFTW_MEASURE);  
    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_1d)(nx,p1,p2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_1d)(nx,p1,p2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_1d)(nx,p1,p1,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_1d)(nx,p1,p1,FFTW_FORWARD,FFTW_MEASURE); 


    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    } 
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 0;
    fft->ndim           = 1;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_1d_real(gmx_fft_t *   pfft,
                     int           nx) 
{
    gmx_fft_t              fft;
    real                 *p1,*p2,*up1,*up2;
    size_t                pc;
    int                   i,j,k;

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(real)*(nx+2));
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(real)*(nx+2));
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    /* make unaligned pointers. 
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (real *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (real *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_c2r_1d)(nx,(FFTWPREFIX(complex) *)up1,up2,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_r2c_1d)(nx,up1,(FFTWPREFIX(complex) *)up2,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_c2r_1d)(nx,(FFTWPREFIX(complex) *)up1,up1,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_r2c_1d)(nx,up1,(FFTWPREFIX(complex) *)up1,FFTW_MEASURE);  

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_c2r_1d)(nx,(FFTWPREFIX(complex) *)p1,p2,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_r2c_1d)(nx,p1,(FFTWPREFIX(complex) *)p2,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_c2r_1d)(nx,(FFTWPREFIX(complex) *)p1,p1,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_r2c_1d)(nx,p1,(FFTWPREFIX(complex) *)p1,FFTW_MEASURE); 


    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    }
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 1;
    fft->ndim           = 1;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_2d(gmx_fft_t *   pfft,
                int           nx, 
                int           ny) 
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1,*p2,*up1,*up2;
    size_t                 pc;
    int                   i,j,k;

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx*ny+2));
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx*ny+2));
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    /* make unaligned pointers. 
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (FFTWPREFIX(complex) *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (FFTWPREFIX(complex) *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_2d)(nx,ny,up1,up2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_2d)(nx,ny,up1,up2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_2d)(nx,ny,up1,up1,FFTW_BACKWARD,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_2d)(nx,ny,up1,up1,FFTW_FORWARD,FFTW_MEASURE);  

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_2d)(nx,ny,p1,p2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_2d)(nx,ny,p1,p2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_2d)(nx,ny,p1,p1,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_2d)(nx,ny,p1,p1,FFTW_FORWARD,FFTW_MEASURE); 
    

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    }
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 0;
    fft->ndim           = 2;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_2d_real(gmx_fft_t *   pfft,
                     int           nx, 
                     int           ny) 
{
    gmx_fft_t              fft;
    real            *p1,*p2,*up1,*up2;
    size_t                pc;
    int                   i,j,k;

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(real)*( nx*(ny/2+1)*2 + 2) );
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(real)*( nx*(ny/2+1)*2 + 2) );
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }

    /* make unaligned pointers. 
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (real *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (real *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx,ny,(FFTWPREFIX(complex) *)up1,up2,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx,ny,up1,(FFTWPREFIX(complex) *)up2,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx,ny,(FFTWPREFIX(complex) *)up1,up1,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx,ny,up1,(FFTWPREFIX(complex) *)up1,FFTW_MEASURE);  
    
    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx,ny,(FFTWPREFIX(complex) *)p1,p2,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx,ny,p1,(FFTWPREFIX(complex) *)p2,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx,ny,(FFTWPREFIX(complex) *)p1,p1,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx,ny,p1,(FFTWPREFIX(complex) *)p1,FFTW_MEASURE); 
    

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    }
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 1;
    fft->ndim           = 2;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_3d(gmx_fft_t *   pfft,
                int           nx, 
                int           ny,
                int           nz) 
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1,*p2,*up1,*up2;
    size_t                 pc;
    int                   i,j,k;

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
    
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx*ny*nz+2));
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx*ny*nz+2));
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    /* make unaligned pointers. 
        * In double precision the actual complex datatype will be 16 bytes,
        * so go to a char pointer and force an offset of 8 bytes instead.
        */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (FFTWPREFIX(complex) *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (FFTWPREFIX(complex) *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,up1,up2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,up1,up2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,up1,up1,FFTW_BACKWARD,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,up1,up1,FFTW_FORWARD,FFTW_MEASURE);  

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,p1,p2,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,p1,p2,FFTW_FORWARD,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,p1,p1,FFTW_BACKWARD,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_3d)(nx,ny,nz,p1,p1,FFTW_FORWARD,FFTW_MEASURE); 
    

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    }
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 0;
    fft->ndim           = 3;
    
    *pfft = fft;
    return 0;
}



int
gmx_fft_init_3d_real(gmx_fft_t *   pfft,
                     int           nx, 
                     int           ny,
                     int           nz) 
{
    gmx_fft_t             fft;
    real            *p1,*p2,*up1,*up2;
    size_t                pc;
    int                   i,j,k;

    if(pfft==NULL)
    {
        gmx_fatal(FARGS,"Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;
        
    if( (fft = FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }    
    
    /* allocate aligned, and extra memory to make it unaligned */
    p1  = FFTWPREFIX(malloc)(sizeof(real)*( nx*ny*(nz/2+1)*2 + 2) );
    if(p1==NULL)
    {
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    p2  = FFTWPREFIX(malloc)(sizeof(real)*( nx*ny*(nz/2+1)*2 + 2) );
    if(p2==NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        return ENOMEM;
    }
    
    /* make unaligned pointers. 
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a void pointer and force an offset of 8 bytes instead.
     */
    pc = (size_t)p1;
    pc += 8; 
    up1 = (real *)pc;
    
    pc = (size_t)p2;
    pc += 8; 
    up2 = (real *)pc;
    
    
    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx,ny,nz,(FFTWPREFIX(complex) *)up1,up2,FFTW_MEASURE); 
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx,ny,nz,up1,(FFTWPREFIX(complex) *)up2,FFTW_MEASURE); 
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx,ny,nz,(FFTWPREFIX(complex) *)up1,up1,FFTW_MEASURE);  
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx,ny,nz,up1,(FFTWPREFIX(complex) *)up1,FFTW_MEASURE);  
    
    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx,ny,nz,(FFTWPREFIX(complex) *)p1,p2,FFTW_MEASURE); 
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx,ny,nz,p1,(FFTWPREFIX(complex) *)p2,FFTW_MEASURE); 
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx,ny,nz,(FFTWPREFIX(complex) *)p1,p1,FFTW_MEASURE); 
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx,ny,nz,p1,(FFTWPREFIX(complex) *)p1,FFTW_MEASURE); 
    

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            for(k=0;k<2;k++)
            {
                if(fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS,"Error initializing FFTW3 plan.");
                    gmx_fft_destroy(fft);
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    return -1;
                }
            }
        }
    }
    
    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);
    
    fft->real_transform = 1;
    fft->ndim           = 3;
    
    *pfft = fft;
    return 0;
}


int 
gmx_fft_1d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);
    
    /* Some checks */
    if( (fft->real_transform == 1) || (fft->ndim != 1) ||
        ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    

    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            in_data,
                            out_data);
    
    return 0;
}


int 
gmx_fft_1d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);
    
    /* Some checks */    
    if( (fft->real_transform != 1) || (fft->ndim != 1) ||
        ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    if(isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    in_data,out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    in_data,out_data);
    }
    
    return 0;
}


int 
gmx_fft_2d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);
    
    /* Some checks */
    if( (fft->real_transform == 1) || (fft->ndim != 2) ||
        ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    

    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            in_data,
                            out_data);
    
    return 0;
}


int 
gmx_fft_2d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);
    
    /* Some checks */
    if( (fft->real_transform != 1) || (fft->ndim != 2) ||
        ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    if(isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    in_data,
                                    out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    in_data,
                                    out_data);
    }
    
    
    return 0;
}


int 
gmx_fft_3d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);
    
    /* Some checks */
    if( (fft->real_transform == 1) || (fft->ndim != 3) ||
        ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }    
    
    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            in_data,
                            out_data);
    
    return 0;
}


int 
gmx_fft_3d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = (((size_t)in_data & (size_t)out_data & 0xf)==0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);
    
    /* Some checks */
    if( (fft->real_transform != 1) || (fft->ndim != 3) ||
        ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }
    
    if(isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    in_data,
                                    out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    in_data,
                                    out_data);
    }
    
    
    return 0;
}


void
gmx_fft_destroy(gmx_fft_t      fft)
{
    int                   i,j,k;

    if(fft != NULL)
    {
        for(i=0;i<2;i++)
        {
            for(j=0;j<2;j++)
            {
                for(k=0;k<2;k++)
                {
                    if(fft->plan[i][j][k] != NULL)
                    {
                        FFTWPREFIX(destroy_plan)(fft->plan[i][j][k]);
                        fft->plan[i][j][k] = NULL;
                    }
                }
            }
        }
        FFTWPREFIX(free)(fft);
    }

}

