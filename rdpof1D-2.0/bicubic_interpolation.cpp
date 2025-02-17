// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#include "bicubic_interpolation.h"

#define BOUNDARY_CONDITION 0

//0 Neumann
//1 Periodic
//2 Symmetric

/**
  *
  * Neumann boundary condition test
  *
**/
int
neumann_bc (int x, int nx, bool & out)
{
    if (x < 0)
    {
        x = 0;
        out = true;
    }
    else if (x >= nx)
    {
        x = nx - 1;
        out = true;
    }

    return x;
}

/**
  *
  * Periodic boundary condition test
  *
**/
int
periodic_bc (int x, int nx, bool & out)
{
    if (x < 0)
    {
        const int n = 1 - (int) (x / (nx + 1));
        const int ixx = x + n * nx;

        x = ixx % nx;
        out = true;
    }
    else if (x >= nx)
    {
        x = x % nx;
        out = true;
    }

    return x;
}


/**
  *
  * Symmetric boundary condition test
  *
**/
int
symmetric_bc (int x, int nx, bool & out)
{
    if (x < 0)
    {
        const int border = nx - 1;
        const int xx = -x;
        const int n = (int) (xx / border) % 2;

        if (n)
            x = border - (xx % border);
        else
            x = xx % border;
        out = true;
    }
    else if (x >= nx)
    {
        const int border = nx - 1;
        const int n = (int) (x / border) % 2;

        if (n)
            x = border - (x % border);
        else
            x = x % border;
        out = true;
    }

    return x;
}


/**
  *
  * Bicubic interpolation in one dimension
  *
**/
double
cubic_interpolation (double v[4],	//interpolation points
                     double x		//point to be interpolated
                    )
{
    return v[1] + 0.5 * x * (v[2] - v[0]
                             + x * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3]
                                    + x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
  *
  * Bicubic interpolation in two dimension
  *
**/
double
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
float
bicubic_interpolation (const float *input,	//image to be interpolated
                       const float uu,		//x component of the vector field
                       const float vv,		//y component of the vector field
                       const int nx,		//width of the image
                       const int ny,		//height of the image
                       const int nz,            //number of channels of the image
                       const int k,  		//actual channel
                       const bool border_out //if true, put zeros outside the region
                      )
{
    const int sx = (uu < 0) ? -1 : 1;
    const int sy = (vv < 0) ? -1 : 1;

    int x, y, mx, my, dx, dy, ddx, ddy;
    bool out = false;

    switch (BOUNDARY_CONDITION)
    {

    case 0:
        x = neumann_bc ((int) uu, nx, out);
        y = neumann_bc ((int) vv, ny, out);
        mx = neumann_bc ((int) uu - sx, nx, out);
        my = neumann_bc ((int) vv - sx, ny, out);
        dx = neumann_bc ((int) uu + sx, nx, out);
        dy = neumann_bc ((int) vv + sy, ny, out);
        ddx = neumann_bc ((int) uu + 2 * sx, nx, out);
        ddy = neumann_bc ((int) vv + 2 * sy, ny, out);
        break;

    case 1:
        x = periodic_bc ((int) uu, nx, out);
        y = periodic_bc ((int) vv, ny, out);
        mx = periodic_bc ((int) uu - sx, nx, out);
        my = periodic_bc ((int) vv - sx, ny, out);
        dx = periodic_bc ((int) uu + sx, nx, out);
        dy = periodic_bc ((int) vv + sy, ny, out);
        ddx = periodic_bc ((int) uu + 2 * sx, nx, out);
        ddy = periodic_bc ((int) vv + 2 * sy, ny, out);
        break;

    case 2:
        x = symmetric_bc ((int) uu, nx, out);
        y = symmetric_bc ((int) vv, ny, out);
        mx = symmetric_bc ((int) uu - sx, nx, out);
        my = symmetric_bc ((int) vv - sx, ny, out);
        dx = symmetric_bc ((int) uu + sx, nx, out);
        dy = symmetric_bc ((int) vv + sy, ny, out);
        ddx = symmetric_bc ((int) uu + 2 * sx, nx, out);
        ddy = symmetric_bc ((int) vv + 2 * sy, ny, out);
        break;

    default:
        x = neumann_bc ((int) uu, nx, out);
        y = neumann_bc ((int) vv, ny, out);
        mx = neumann_bc ((int) uu - sx, nx, out);
        my = neumann_bc ((int) vv - sx, ny, out);
        dx = neumann_bc ((int) uu + sx, nx, out);
        dy = neumann_bc ((int) vv + sy, ny, out);
        ddx = neumann_bc ((int) uu + 2 * sx, nx, out);
        ddy = neumann_bc ((int) vv + 2 * sy, ny, out);
        break;
    }

    if (out && border_out)

        return 0.0;

    else
    {
        //obtain the interpolation points of the image
        const float p11 = input[(mx  + nx * my) * nz + k];
        const float p12 = input[(x   + nx * my) * nz + k];
        const float p13 = input[(dx  + nx * my) * nz + k];
        const float p14 = input[(ddx + nx * my) * nz + k];

        const float p21 = input[(mx  + nx * y) * nz + k];
        const float p22 = input[(x   + nx * y) * nz + k];
        const float p23 = input[(dx  + nx * y) * nz + k];
        const float p24 = input[(ddx + nx * y) * nz + k];

        const float p31 = input[(mx  + nx * dy) * nz + k];
        const float p32 = input[(x   + nx * dy) * nz + k];
        const float p33 = input[(dx  + nx * dy) * nz + k];
        const float p34 = input[(ddx + nx * dy) * nz + k];

        const float p41 = input[(mx  + nx * ddy) * nz + k];
        const float p42 = input[(x   + nx * ddy) * nz + k];
        const float p43 = input[(dx  + nx * ddy) * nz + k];
        const float p44 = input[(ddx + nx * ddy) * nz + k];

        //create array
        double pol[4][4] =
        {
            {p11, p21, p31, p41}, {p12, p22, p32, p42},
            {p13, p23, p33, p43}, {p14, p24, p34, p44}
        };

        //return interpolation
        return bicubic_interpolation (pol, (float) uu - x, (float) vv - y);
    }
}


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void
bicubic_interpolation (const float *input,	//image to be warped
                       const float *u,		//x component of the vector field
                       float *output,		//warped output image with bicubic interpolation
                       const int nx,		//width of the image
                       const int ny,		//height of the image
                       const int nz,		//number of channels of the image
                       bool border_out	//if true, put zeros outside the region
                      )
{
    for(int k = 0; k < nz; k++)
    {
        #pragma omp parallel for
        for (int i = 0; i < ny; i++)
            for (int j = 0; j < nx; j++)
            {

                const int p = i * nx + j;
                const float uu = (float) (j + u[p]);
                const float vv = (float) (i); //(float) (i + v[p]);

                //obtain the bicubic interpolation at position (uu, vv)
                output[p * nz + k] = bicubic_interpolation (input, uu, vv, nx, ny, nz, k, border_out);

            }
    } // end multi-channel loop
}
