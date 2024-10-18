// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2013,2015 Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// All rights reserved.

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>


extern "C"
{
#include "iio.h"
}

#include "occultation.h"
#include "robust_expo_methods.h"
#include "rdpof1D.h"

#define DEFAULT_NPROC 1
#define DEFAULT_METHOD 1
#define DEFAULT_ALPHA 50
#define DEFAULT_GAMMA 10
#define DEFAULT_LAMBDA 0.2
#define DEFAULT_NSCALES 10
#define DEFAULT_ZFACTOR 0.5
#define DEFAULT_TOL 0.0001
#define DEFAULT_INNER_ITER 1
#define DEFAULT_OUTER_ITER 15
#define DEFAULT_VERBOSE 0
#define DEFAULT_TEST_LR        1
#define DEFAULT_NAN_THRESHOLD  0.01
#define DEFAULT_NAN_MASK "mask_nan.png"
using namespace std;

/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
//bool read_image(const char *fname, float **f, int &nx, int &ny, int &nz){
//
//	*f = iio_read_image_float_vec(fname, &nx, &ny, &nz);
//
//	return *f ? true : false;
//}
//
bool
read_image (const char *fname, float **f, int *w, int *h, int *c)
{
    *f = iio_read_image_float_vec (fname, w, h, c);

    int nlig, ncol, i, j;
    double somme = 0.0;

    nlig = *h;
    ncol = *w;
    for(i=0; i<nlig; i++)
    {
        for(j=0; j<ncol; j++)
        {
            if (isnormal((*f)[i*ncol + j]) )
            {
                somme += (*f)[i*ncol + j];
            }
        }
    }

    // moyenne
    somme /= (nlig * ncol);

    // qui va remplacer les valeurs à la noix
    for(i=0; i<nlig; i++)
    {
        for(j=0; j<ncol; j++)
        {
            if (! isnormal((*f)[i*ncol + j]) )
            {
                printf("remplacement de %f (*f)[i*ncol + j]) par %8.4e\n",
                       (*f)[i*ncol + j], somme);

                (*f)[i*ncol + j] = somme;
            }
        }
    }
    return *f ? true : false;
}
/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -I1          first image
 *   -I2          second image
 *   -out_file    name of the output flow field
 *   -processors  number of threads used with the OpenMP library
 *   -method_type Choose Diffusion Tensor to use with the method (DF, DF-Beta or DF-Auto)
 *   -alpha       smoothing parameter
 *   -gamma       gradient constancy parameter
 *   -lambda	  Coeficient for exponential smoothing factor
 *   -nscales     number of scales for the pyramidal approach
 *   -zoom_factor reduction factor for creating the scales
 *   -TOL         stopping criterion threshold for the iterative process
 *   -inner_iter  number of inner iterations
 *   -outer_iter  number of outer iterations
 *   -verbose     Switch on/off messages
 *   -test_lr
 *   -nan_threshold
 *   -nan_mask
 */
int main (int argc, char *argv[])
{

    if (argc < 3)
    {

//    cout << "Usage: " << argv[0]
//	      << " I1 I2 [out_file processors"
//	      << " method_type alpha gamma lambda"
//	      << " nscales zoom_factor TOL"
//	      << " inner_iter outer_iter verbose]"
//	      << endl;
//
        cout << "Usage: " << argv[0]
             << " I1 I2 out_file "
             << " processors " << DEFAULT_NPROC
             << " method_type " << DEFAULT_METHOD
             << " alpha " << DEFAULT_ALPHA
             << " gamma " << DEFAULT_GAMMA
             << " lambda " << DEFAULT_LAMBDA
             << " nscales " << DEFAULT_NSCALES
             << " zoom_factor " << DEFAULT_ZFACTOR
             << " TOL " << DEFAULT_TOL
             << " inner_iter " << DEFAULT_INNER_ITER
             << " outer_iter " << DEFAULT_OUTER_ITER
             << " verbose " << DEFAULT_VERBOSE
             << " test_lr " << DEFAULT_TEST_LR
             << " nan_threshold " << DEFAULT_NAN_THRESHOLD
             << " nan_mask "<< DEFAULT_NAN_MASK
             << endl;
    }
    else
    {

        int i = 1, nx, ny, nz, nx1, ny1, nz1;

        //read parameters from the console
        const char *image1  = argv[i++];
        const char *image2  = argv[i++];
        const char *outfile = (argc >= 4) ? argv[i++] : "flow.flo";
        int   nproc         = (argc > i)  ? atoi (argv[i++]) : DEFAULT_NPROC;

        int   method_type = (argc > i) ? atoi (argv[i++]) : DEFAULT_METHOD;
        float alpha       = (argc > i) ? atof (argv[i++]) : DEFAULT_ALPHA;
        float gamma       = (argc > i) ? atof (argv[i++]) : DEFAULT_GAMMA;
        float lambda     	= (argc > i) ? atof (argv[i++])  : DEFAULT_LAMBDA;

        int   nscales     = (argc > i) ? atoi (argv[i++]) : DEFAULT_NSCALES;
        float zfactor     = (argc > i) ? atof (argv[i++]) : DEFAULT_ZFACTOR;
        float TOL         = (argc > i) ? atof (argv[i++]) : DEFAULT_TOL;
        int   initer      = (argc > i) ? atoi (argv[i++]) : DEFAULT_INNER_ITER;
        int   outiter     = (argc > i) ? atoi (argv[i++]) : DEFAULT_OUTER_ITER;
        int   verbose	= (argc > i) ? atoi (argv[i++])   : DEFAULT_VERBOSE;
        int test_lr = (argc > i) ? atoi (argv[i++]) : DEFAULT_TEST_LR;
        double seuil = (argc > i) ? atof (argv[i++]) : DEFAULT_NAN_THRESHOLD;
        const char *masque = (argc > i) ? argv[i++] : DEFAULT_NAN_MASK;

        //check parameters
        if (nproc > 0)                          omp_set_num_threads (nproc);
        if (method_type < 1 || method_type > 3) method_type = DEFAULT_METHOD;
        if (alpha <= 0) 	                      alpha       = DEFAULT_ALPHA;
        if (gamma < 0)                          gamma       = DEFAULT_GAMMA;
        if (lambda < 0)                         lambda      = DEFAULT_LAMBDA;
        if (nscales <= 0)                       nscales     = DEFAULT_NSCALES;
        if (zfactor <= 0 || zfactor >= 1)       zfactor     = DEFAULT_ZFACTOR;
        if (TOL <= 0)                           TOL         = DEFAULT_TOL;
        if (initer <= 0)                        initer      = DEFAULT_INNER_ITER;
        if (outiter <= 0)                       outiter     = DEFAULT_OUTER_ITER;

        float *I1, *I2;

        //read the input images
        bool correct1 = read_image (image1, &I1, &nx, &ny, &nz);
        bool correct2 = read_image (image2, &I2, &nx1, &ny1, &nz1);

        // if the images are correct, compute the optical flow
        if (correct1 && correct2 && nx == nx1 && ny == ny1 && nz == nz1)
        {

            //set the number of scales according to the size of the
            //images.  The value N is computed to assure that the smaller
            //images of the pyramid don't have a size smaller than 16x16
            const float N =  1 + log (std::min (nx, ny) / 16.) / log (1. / zfactor);
            if ((int) N < nscales) nscales = (int) N;

            cout  << endl
                  << " ncores:" << nproc   << " method_type:" << method_type
                  << " alpha:"  << alpha   << " gamma:"       << gamma << " lambda:" << lambda
                  << " scales:" << nscales << " nu:"          << zfactor << " TOL:"  << TOL
                  << " inner:"  << initer  << " outer:"       << outiter << " LR:" << test_lr
                  << " seuil:"  << seuil   << " masque:"      << masque
                  << endl;

            //allocate memory for the flow
            float *u1 = new float[nx * ny];

            //compute the optic flow
            robust_expo_methods(
                I1, I2, u1, nx, ny, nz,
                method_type, alpha, gamma, lambda,
                nscales, zfactor, TOL, initer, outiter, verbose
            );
            if (test_lr)
            {
                cout << "méthode RDP OF avec test gauche/droite\n" << endl;
                // allocate memory for the flow
                float *u2 = new float[nx * ny];

                // compute the optic flow between I_right and I_left
                //compute the optic flow
                robust_expo_methods(I2, I1, u2, nx, ny, nz,
                                    method_type, alpha, gamma, lambda,
                                    nscales, zfactor, TOL, initer, outiter, verbose
                                   );

                // determination of the occulted pixels of the left disparity map
                set_nan_disparities(u1, u2, seuil, nx, ny);
                float *mask = new float[nx * ny];
                int c, l;
                for (l=0; l<ny; l++)
                {
                    for (c=0; c<nx; c++)
                    {
                        mask[l*nx + c] = (isnan(u1[l*nx + c])) ? 0 : 1;
                    }
                }
                iio_save_image_float ((char *) masque, mask, nx, ny);

                delete[]mask;
                delete[]u2;
            }
            else
            {
                cout << "méthode RDP OF sans test gauche/droite\n" << endl;
            }

            iio_save_image_float ((char *) outfile, u1, nx, ny);


            //free dynamic memory
            free (I1);
            free (I2);
            delete[]u1;

//          //save the flow
//	  float *f = new float[nx * ny * 2];
//	  for (int i = 0; i < nx * ny; i++){
//	      f[2 * i] = u[i];
//	      f[2 * i + 1] = v[i];
//          }
//	  iio_save_image_float_vec ((char *) outfile, f, nx, ny, 2);
//
//	  //free dynamic memory
//	  free (I1);
//	  free (I2);
//	  delete[]u;
//	  delete[]v;
//	  delete[]f;
//
        }
        else cerr << "Cannot read the images or the size of the images are not equal" << endl;

    }

    return 0;
}
