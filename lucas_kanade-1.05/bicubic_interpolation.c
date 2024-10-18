/*
 * bicubic_interpolation.c
 *
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the simplified BSD License. You should have received a
 * copy of this license along this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 *
 * Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
 * Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
 * All rights reserved.
 */

#include "bicubic_interpolation.h"

#define BOUNDARY_CONDITION 0

//0 Neumann

/**
  *
  * Neumann boundary condition test
  *
**/
int
neumann_bc (int x, int nx, int * out)
{
    if (x < 0)
    {
        x = 0;
        *out = 1;
    }
    else if (x >= nx)
    {
        x = nx - 1;
        *out = 1;
    }

    return x;
}

/**
  *
  * Bicubic interpolation in one dimension
  *
**/
inline double
cubic_interpolation (double v[4],	//interpolation points
                     double x		//point to be interpolated
                    )
{
	double  val;
	val = v[1] + 0.5 * x * (v[2] - v[0]
							+ x * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3]
								   + x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
	return val;
}


/**
  *
  * Bicubic interpolation in two dimension
  *
**/
inline double
bicubic_interpolation (double p[4][4],	//array containing the interpolation points
                       double x,	//x position to be interpolated
                       double y		//y position to be interpolated
                      )
{
    double v[4];
    v[0] = cubic_interpolation (p[0], y);
    v[1] = cubic_interpolation (p[1], y);
    v[2] = cubic_interpolation (p[2], y);
    v[3] = cubic_interpolation (p[3], y);
    return cubic_interpolation (v, x);
}

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
					  const int border_out //if true, put zeros outside the region
	)
{
    const int sx = (uu < 0) ? -1 : 1;
    const int sy = (vv < 0) ? -1 : 1;

    int x, y, mx, my, dx, dy, ddx, ddy;
    int out = 0;
	
	x = neumann_bc ((int) uu, nx, &out);
	y = neumann_bc ((int) vv, ny, &out);
	mx = neumann_bc ((int) uu - sx, nx, &out);
	my = neumann_bc ((int) vv - sx, ny, &out);
	dx = neumann_bc ((int) uu + sx, nx, &out);
	dy = neumann_bc ((int) vv + sy, ny, &out);
	ddx = neumann_bc ((int) uu + 2 * sx, nx, &out);
	ddy = neumann_bc ((int) vv + 2 * sy, ny, &out);
	

    if (out && border_out)
        return 0.0;
    else
    {
        //obtain the interpolation points of the image
        const double p11 = input[(mx  + nx * my) * nz + k];
        const double p12 = input[(x   + nx * my) * nz + k];
        const double p13 = input[(dx  + nx * my) * nz + k];
        const double p14 = input[(ddx + nx * my) * nz + k];

        const double p21 = input[(mx  + nx * y) * nz + k];
        const double p22 = input[(x   + nx * y) * nz + k];
        const double p23 = input[(dx  + nx * y) * nz + k];
        const double p24 = input[(ddx + nx * y) * nz + k];

        const double p31 = input[(mx  + nx * dy) * nz + k];
        const double p32 = input[(x   + nx * dy) * nz + k];
        const double p33 = input[(dx  + nx * dy) * nz + k];
        const double p34 = input[(ddx + nx * dy) * nz + k];

        const double p41 = input[(mx  + nx * ddy) * nz + k];
        const double p42 = input[(x   + nx * ddy) * nz + k];
        const double p43 = input[(dx  + nx * ddy) * nz + k];
        const double p44 = input[(ddx + nx * ddy) * nz + k];

        //create array
        double pol[4][4] =
        {
            {p11, p21, p31, p41}, {p12, p22, p32, p42},
            {p13, p23, p33, p43}, {p14, p24, p34, p44}
        };

        //return interpolation
        return bicubic_interpolation (pol, (double) uu - x, (double) vv - y);
    }
}


