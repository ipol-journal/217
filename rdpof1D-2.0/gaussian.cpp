// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#include <cmath>
#include <iostream>

#include "gaussian.h"

/**
 *
 * Convolution with a Gaussian
 *
 */
void
gaussian (float *I,		//input/output image
          const int xdim,	//image width
          const int ydim,	//image height
          const int zdim,       // number of color channels in the image
          const double sigma,	//Gaussian sigma
          const int bc,	//boundary condition
          const int precision //defines the size of the window
         )
{

    int i, j, k;

    const double den = 2 * sigma * sigma;
    const int size = (int) (precision * sigma) + 1;
    const int bdx = xdim + size;
    const int bdy = ydim + size;

    if (bc && size > xdim)
    {
        std::cerr << "GaussianSmooth: sigma too large for this bc\n" << std::endl;
        throw 1;
    }

    // compute the coefficients of the 1D convolution kernel
    double *B = new double[size];
    for (int i = 0; i < size; i++)
        B[i] = 1 / (sigma * sqrt (2.0 * 3.1415926)) * exp (-i * i / den);

    double norm = 0;

    // normalize the 1D convolution kernel
    for (int i = 0; i < size; i++)
        norm += B[i];

    norm *= 2;

    norm -= B[0];

    for (int i = 0; i < size; i++)
        B[i] /= norm;



    double *R = new double[size + xdim + size];
    double *T = new double[size + ydim + size];

    //Loop for every channel
    for(int index_multichannel = 0; index_multichannel < zdim; index_multichannel++)
    {

        // convolution of each line of the input image
        for (k = 0; k < ydim; k++)
        {
            for (i = size; i < bdx; i++)
                R[i] = I[(k * xdim + i - size) * zdim + index_multichannel];
            switch (bc)
            {
            case 0:		// Dirichlet boundary conditions

                for (i = 0, j = bdx; i < size; i++, j++)
                    R[i] = R[j] = 0;
                break;

            case 1:		// Reflecting boundary conditions

                for (i = 0, j = bdx; i < size; i++, j++)
                {
                    R[i] = I[(k * xdim + size - i ) * zdim + index_multichannel];
                    R[j] = I[(k * xdim + xdim - i - 1) * zdim + index_multichannel ];
                }
                break;

            case 2:		// Periodic boundary conditions

                for (i = 0, j = bdx; i < size; i++, j++)
                {
                    R[i] = I[(k * xdim + xdim - size + i) * zdim + index_multichannel];
                    R[j] = I[(k * xdim + i) * zdim + index_multichannel];
                }
                break;
            }

            for (i = size; i < bdx; i++)
            {

                double sum = B[0] * R[i];

                for (int j = 1; j < size; j++)
                    sum += B[j] * (R[i - j] + R[i + j]);

                I[(k * xdim + i - size) * zdim + index_multichannel] = sum;

            }
        }

        // convolution of each column of the input image

        for (k = 0; k < xdim; k++)
        {
            for (i = size; i < bdy; i++)
                T[i] = I[((i - size) * xdim + k) * zdim + index_multichannel];

            switch (bc)
            {
            case 0:		// Dirichlet boundary conditions

                for (i = 0, j = bdy; i < size; i++, j++)
                    T[i] = T[j] = 0;
                break;

            case 1:		// Reflecting boundary conditions

                for (i = 0, j = bdy; i < size; i++, j++)
                {
                    T[i] = I[((size - i) * xdim + k) * zdim + index_multichannel];
                    T[j] = I[((ydim - i - 1) * xdim + k) * zdim + index_multichannel];
                }
                break;

            case 2:		// Periodic boundary conditions

                for (i = 0, j = bdx; i < size; i++, j++)
                {
                    T[i] = I[((ydim - size + i) * xdim + k) * zdim + index_multichannel];
                    T[j] = I[(i * xdim + k) * zdim + index_multichannel];
                }
                break;
            }

            for (i = size; i < bdy; i++)
            {
                double sum = B[0] * T[i];

                for (j = 1; j < size; j++)
                    sum += B[j] * (T[i - j] + T[i + j]);

                I[((i - size) * xdim + k) * zdim + index_multichannel] = sum;

            }
        }

    }

    delete[]B;
    delete[]R;
    delete[]T;
}

