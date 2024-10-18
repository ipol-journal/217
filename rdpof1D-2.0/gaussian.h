// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>
// All rights reserved.


#ifndef GAUSSIAN_H
#define GAUSSIAN_H

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
          const int bc = 1,	//boundary condition
          const int precision = 5	//defines the size of the window
         );

#endif
