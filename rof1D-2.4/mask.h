// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
//               2016, Tristan Dagobert, CMLA
// All rights reserved.


#ifndef MASK_H
#define MASK_H

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
     const int ny		//image height
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
     const int ny		//image height
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
     const int ny		//image height
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
          const int ny		//image height
         );


/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void
psi_divergence (const float *psi,	//robust functional
                float *psi1,		//coefficients of divergence
                float *psi2,		//coefficients of divergence
                float *psi3,		//coefficients of divergence
                float *psi4,		//coefficients of divergence
                const int nx,		//image width
                const int ny		//image height
               );



/**
 *
 * Compute the divergence of the optical flow
 *
 */
void
divergence_u (const float *u,		//x component of optical flow
              const float *v,		//y component of optical flow
              const float *psi1,	//coefficients of divergence
              const float *psi2,	//coefficients of divergence
              const float *psi3,	//coefficients of divergence
              const float *psi4,	//coefficients of divergence
              float *div_u,		//computed divergence for u
              float *div_v,		//computed divergence for v
              const int nx,		//image width
              const int ny		//image height
             );

#endif
