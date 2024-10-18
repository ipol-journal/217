// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
//               2016, Tristan Dagobert, CMLA
// All rights reserved.

/*===========================================================================*/
int
comparaison (const void *p1, const void *p2);
/*===========================================================================*/
void leftright_test(float *left_dx, float *right_dx, int ncol, int nlig,
                    float threshold);
/*===========================================================================*/
/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
bool
read_image (const char *fname, float **f, int *w, int *h);
/*===========================================================================*/
