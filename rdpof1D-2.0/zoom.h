// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>

// All rights reserved.


#ifndef ZOOM_H
#define ZOOM_H

/**
  *
  * Compute the size of a zoomed image from the zoom factor
  *
**/
void
zoom_size (int nx,		//width of the orignal image
           int ny,		//height of the orignal image
           int &nxx,		//width of the zoomed image
           int &nyy,		//height of the zoomed image
           float factor = 0.5	//zoom factor between 0 and 1
          );

/**
  *
  * Function to downsample the image
  *
**/
void
zoom_out (const float *I,	//input image
          float *Iout,		//output image
          const int nx,		//image width
          const int ny,		//image height
          const int nz,		// number of color channels in image
          const float factor = 0.5	//zoom factor between 0 and 1
         );


/**
  *
  * Function to upsample the flow
  *
**/

void
zoom_in_flow (const float *I,	//input image
              float *Iout,		//output image
              int nx,		//width of the original image
              int ny,		//height of the original image
              int nxx,		//width of the zoomed image
              int nyy		//height of the zoomed image
             );


#endif
