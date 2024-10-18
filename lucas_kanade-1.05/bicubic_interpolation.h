/* 
 * bicubic_interpolation.h
 *
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the simplified BSD License. You should have received a
 * copy of this license along this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 *
 * Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
 * Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
 * Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>
 * All rights reserved.
 */

#ifndef BICUBIC_INTERPOLATION_H
#define BICUBIC_INTERPOLATION_H

// 0 Neumann

/**
  *
  * Neumann boundary condition test
  *
**/
int
neumann_bc (int x, int nx, int * out);

/**
  *
  * Bicubic interpolation in one dimension
  *
**/
double
cubic_interpolation (double v[4],	//interpolation points
                     double x		//point to be interpolated
                    );


/**
  *
  * Bicubic interpolation in two dimension
  *
**/
double
bicubic_interpolation (double p[4][4],	//array containing the interpolation points
                       double x,	//x position to be interpolated
                       double y		//y position to be interpolated
                      );

/**
  *
  * Compute the bicubic interpolation of a point in an image.
  * Detects if the point goes outside the image domain
  *
**/
double
bicubic_interp_point (const double *input,	//image to be interpolated
                       const double uu,		//x component of the vector field
                       const double vv,		//y component of the vector field
                       const int nx,		//width of the image
                       const int ny,		//height of the image
                       const int nz,            //number of channels of the image
                       const int k,  		//actual channel
                       const int border_out	//if true, put zeros outside the region
                      );

#endif
