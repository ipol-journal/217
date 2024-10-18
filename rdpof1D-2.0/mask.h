// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López  <nmonzon@ctim.es>
// Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>
// All rights reserved.

#ifndef MASK_H
#define MASK_H

#include <omp.h>


/**
 *
 * Function to apply a 3x3 mask to an image
 *
 */
void
mask3x3 (const float *input,	//input image
         float *output,		//output image
         const int nx,		//image width
         const int ny,		//image height
         const int nz,          // number of color channels in the image
         const float *mask	//mask to be applied
        );

/**
 *
 * Compute the second order X derivative
 *
 */
void
Dxx (const float *I,		//input image
     float *Ixx,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image
    );


/**
 *
 * Compute the second order Y derivative
 *
 */
void
Dyy (const float *I,		//input image
     float *Iyy,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image
    );


/**
 *
 * Compute the second order XY derivative
 *
 */
void
Dxy (const float *I,		//input image
     float *Ixy,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image
    );


/**
 *
 * Compute the gradient with central differences
 *
 */
void
gradient (const float *input,	//input image
          float *dx,		//computed x derivative
          float *dy,		//computed y derivative
          const int nx,		//image width
          const int ny,		//image height
          const int nz          //number of color channels in the image
         );


#endif
