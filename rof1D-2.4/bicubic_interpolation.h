// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
//               2016, Tristan Dagobert, CMLA
// All rights reserved.


#ifndef BICUBIC_INTERPOLATION_H
#define BICUBIC_INTERPOLATION_H


//0 Neumann
//1 Periodic
//2 Symmetric

/**
  *
  * Neumann boundary condition test
  *
**/
int
neumann_bc (int x, int nx, bool & out);

/**
  *
  * Periodic boundary condition test
  *
**/
int
periodic_bc (int x, int nx, bool & out);

/**
  *
  * Symmetric boundary condition test
  *
**/
int
symmetric_bc (int x, int nx, bool & out);

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
float
bicubic_interpolation (const float *input,	//image to be interpolated
                       const float uu,		//x component of the vector field
                       const float vv,		//y component of the vector field
                       const int nx,		//width of the image
                       const int ny,		//height of the image
                       const bool border_out = false	//if true, put zeros outside the region
                      );


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void
bicubic_interpolation (const float *input,	//image to be warped
                       const float *u,		//x component of the vector field
                       const float *v,		//y component of the vector field
                       float *output,		//warped output image with bicubic interpolation
                       const int nx,		//width of the image
                       const int ny,		//height of the image
                       bool border_out = false	//if true, put zeros outside the region
                      );


#endif
