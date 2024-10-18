// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>
// All rights reserved.

#ifndef SMOOTHNESS_H
#define SMOOTHNESS_H

/**
  *
  * Compute the coefficients of the robust functional (smoothness term)
  *
**/
void psi_smooth(
    const float *ux, //gradient of x component of the optical flow
    const float *uy, //gradient of x component of the optical flow
    const float *expo, //exponential smoothing factor
    const int   size_flow,
    float       *psi       //output coefficients
);

void max_gradients(
    const float *Ix, // Computed Image 1 derivative in x
    const float *Iy, // Computed Image 1 derivative in y
    const int   size, // Total image size (height * weight * nchannels)
    const int   nz, // nº channels
    float *maximum_gradients_per_pixel // vector with the maximum gradients per pixel

);
/**
  * Calculate the lambda optimum using the maximum gradient from all the multi-channel image
  * It also return the maximum gradient per pixel
**/
float lambda_optimum_using_maximum_gradient_per_pixel(
    const float *Ix, // Computed Image 1 derivative in x
    const float *Iy, // Computed Image 1 derivative in y
    const int   size, // Total image size (height * weight * nchannels)
    const int   size_flow, // Total flow size (height * weight)
    const int   nz, // nº channels
    const float alpha,  // smoothness weight
    float *lambda_per_pixel, // local lambda per pixel
    float *maximum_gradients_per_pixel // vector with the maximum gradients per pixel
);

/**
 **  Calculate the exponential values.
 **
**/
void exponential_calculation(
    const float *Ix,    	 // Computed Image 1 derivative in x
    const float *Iy,    	 // Computed Image 1 derivative in y
    const int   size_flow,   // size of the flow field
    const int   size,    	 // size of the multi-channel image
    const int   nz,     	 // nº of image channels
    const float alpha,  	 // smoothness weight
    const float lambda,      // Coeffient for decreasing function
    const int   method_type, // (1 = DF, 2 = DF_BETA, 3 = DF_AUTO)
    float       *expo        // e^(lambda * DI)
);

#endif
