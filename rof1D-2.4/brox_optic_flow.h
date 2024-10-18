// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
//               2016, Tristan Dagobert, CMLA
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

/**
  *
  * Compute the coefficients of the robust functional (data term)
  *
**/
void
psi_data (const float *I1,	//first image
          const float *I2,	//second image
          const float *I2x,	//gradient of the second image
          const float *I2y,	//gradient of the second image
          const float *du,	//motion increment
          const float *dv,	//motion increment
          float *psip,		//output coefficients
          const int nx,		//image width
          const int ny		//image height
         );

/**
  *
  * Compute the coefficients of the robust functional (gradient term)
  *
**/
void
psi_gradient (const float *I1x,	//gradient of the first image
              const float *I1y,		//gradient of the first image
              const float *I2x,		//gradient of the second image
              const float *I2y,		//gradient of the second image
              const float *I2xx,	//second derivatives of the second image
              const float *I2xy,	//second derivatives of the second image
              const float *I2yy,	//second derivatives of the second image
              const float *du,		//motion increment
              const float *dv,		//motion increment
              float *psip,		//output coefficients
              const int nx,		//image width
              const int ny		//image height
             );


/**
  *
  * Compute the coefficients of the robust functional (smoothness term)
  *
**/
void
psi_smooth (const float *ux,	//gradient of x component of the optical flow
            const float *uy,	//gradient of x component of the optical flow
            const float *vx,	//gradient of y component of the optical flow
            const float *vy,	//gradient of y component of the optical flow
            float *psi,		//output coefficients
            const int nx,	//image width
            const int ny	//image height
           );


/**
 *
 *  SOR iteration in one position
 *
 */
inline float
sor_iteration (const float *Au,	//constant part of the numerator of u
               const float *Av,		//constant part of the numerator of v
               const float *Du,		//denominator of u
               const float *Dv,		//denominator of v
               const float *D,		//constant part of the numerator
               float *du,		//x component of the motion increment
               float *dv,		//y component of the motion increment
               const float alpha,	//alpha smoothness parameter
               const float *psi1,	//coefficients of the divergence
               const float *psi2,
               const float *psi3,
               const float *psi4,
               const int i,		//current row
               const int i0,		//previous row
               const int i1,		//following row
               const int j,		//current column
               const int nx,		//number of columns
               const int j0,		//previous column
               const int j1		//following column
              );


/**
  *
  * Compute the optic flow with the Brox spatial method
  *
  *
**/
void
brox_optic_flow (const float *I1,	//first image
                 const float *I2,	//second image
                 float *u,		//x component of the optical flow
                 float *v,		//y component of the optical flow
                 const int nx,		//image width
                 const int ny,		//image height
                 const float alpha,	//smoothness parameter
                 const float gamma,	//gradient term parameter
                 const float TOL,	//stopping criterion threshold
                 const int inner_iter,	//number of inner iterations
                 const int outer_iter,	//number of outer iterations
                 const int    number_of_threads, // number of threads for the parallel code
                 const bool verbose	//switch on messages
                );

/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
void
image_normalization (const float *I1,	//input image 1
                     const float *I2,	//input image 2
                     float *I1n,	//normalized output image 1
                     float *I2n,	//normalized output image 2
                     int size		//size of the image
                    );


/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void
brox_optic_flow (const float *I1,	//first image
                 const float *I2,	//second image
                 float *u,		//x component of the optical flow
                 float *v,		//y component of the optical flow
                 const int nxx,		//image width
                 const int nyy,		//image height
                 const float alpha,	//smoothness parameter
                 const float gamma,	//gradient term parameter
                 const int nscales,	//number of scales
                 const float nu,	//downsampling factor
                 const float TOL,	//stopping criterion threshold
                 const int inner_iter,	//number of inner iterations
                 const int outer_iter,	//number of outer iterations
                 const bool verbose	//switch on messages
                );

#endif
