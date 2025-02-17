// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

extern "C"
{
#include "iio.h"
}

#include "occultation.h"
#include "brox_optic_flow.h"
#include "rof1D.h"

#define DEFAULT_NPROC          0
#define DEFAULT_ALPHA         18
#define DEFAULT_GAMMA          7
#define DEFAULT_NSCALES      100
#define DEFAULT_ZFACTOR        0.75
#define DEFAULT_TOL            0.0001
#define DEFAULT_INNER_ITER     1
#define DEFAULT_OUTER_ITER    15
#define DEFAULT_VERBOSE        0
#define DEFAULT_TEST_LR        1
#define DEFAULT_NAN_THRESHOLD  0.01

using namespace std;
/*===========================================================================*/
int
comparaison (const void *p1, const void *p2)
{
    float *a = (float *) p1;
    float *b = (float *) p2;

    if (*a < *b)
        return -1;
    if (*a == *b)
        return 0;
    if (*a > *b)
        return +1;

    return 0;
}
/*===========================================================================*/
void leftright_test(float *left_dx, float *right_dx, int ncol, int nlig,
                    float threshold)
{
    int p, left_p;
    int left_x;
    float right_p;

    for(int y=0; y<nlig; y++)
    {
        for(int x=0; x<ncol; x++)
        {
            p = x+y*ncol;

            left_x = round(x + left_dx[p]);
            if(0<=left_x && left_x<ncol)
            {
                left_p = left_x + y*ncol;
                right_p = left_x + right_dx[left_p];
                if (fabs(right_p-x) > threshold)
                {
                    left_dx[p]  = nanf("");
                }
            }
            else
            {
                left_dx[p]  = nanf("");
            }
        }
    }
    return;
}

/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
bool
read_image (const char *fname, float **f, int *w, int *h)
{
    *f = iio_read_image_float (fname, w, h);

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
 *   -alpha       smoothing parameter
 *   -gamma       gradient constancy parameter
 *   -nscales     number of scales for the pyramidal approach
 *   -zoom_factor reduction factor for creating the scales
 *   -TOL         stopping criterion threshold for the iterative process
 *   -inner_iter  number of inner iterations
 *   -outer_iter  number of outer iterations
 *   -verbose     switch on/off messages
 *
 */
int
main (int argc, char *argv[])
{
    if (argc < 3)
    {
        cout << "Usage: " << argv[0]
             << " I1 I2 out_file "
             << " processors " << DEFAULT_NPROC
             << " alpha " << DEFAULT_ALPHA
             << " gamma " << DEFAULT_GAMMA
             << " nscales " << DEFAULT_NSCALES
             << " zoom_factor " << DEFAULT_ZFACTOR
             << " TOL " << DEFAULT_TOL
             << " inner_iter " << DEFAULT_INNER_ITER
             << " outer_iter " << DEFAULT_OUTER_ITER
             << " verbose " << DEFAULT_VERBOSE
             << " test_lr " << DEFAULT_TEST_LR
             << " nan_threshold " << DEFAULT_NAN_THRESHOLD
             << endl;
    }
    else
    {
        int i = 1;

        //read parameters from the console
        const char *image1 = argv[i++];
        const char *image2 = argv[i++];
        const char *outfile = (argc >= 4) ? argv[i++] : "flow.flo";
        int nproc = (argc > i) ? atoi (argv[i++]) : DEFAULT_NPROC;
        float alpha = (argc > i) ? atof (argv[i++]) : DEFAULT_ALPHA;
        float gamma = (argc > i) ? atof (argv[i++]) : DEFAULT_GAMMA;
        int nscales = (argc > i) ? atoi (argv[i++]) : DEFAULT_NSCALES;
        float zfactor = (argc > i) ? atof (argv[i++]) : DEFAULT_ZFACTOR;
        float TOL = (argc > i) ? atof (argv[i++]) : DEFAULT_TOL;
        int initer = (argc > i) ? atoi (argv[i++]) : DEFAULT_INNER_ITER;
        int outiter = (argc > i) ? atoi (argv[i++]) : DEFAULT_OUTER_ITER;
        int verbose = (argc > i) ? atoi (argv[i++]) : DEFAULT_VERBOSE;
        int test_lr = (argc > i) ? atoi (argv[i++]) : DEFAULT_TEST_LR;
        double seuil = (argc > i) ? atof (argv[i++]) : DEFAULT_NAN_THRESHOLD;


        //check parameters
        if (nproc > 0)
            omp_set_num_threads (nproc);
        if (alpha <= 0)
            alpha = DEFAULT_ALPHA;
        if (gamma < 0)
            gamma = DEFAULT_GAMMA;
        if (nscales <= 0)
            nscales = DEFAULT_NSCALES;
        if (zfactor <= 0 || zfactor >= 1)
            zfactor = DEFAULT_ZFACTOR;
        if (TOL <= 0)
            TOL = DEFAULT_TOL;
        if (initer <= 0)
            initer = DEFAULT_INNER_ITER;
        if (outiter <= 0)
            outiter = DEFAULT_OUTER_ITER;

        int nx, ny, nx1, ny1;

        float *I1, *I2;

        //read the input images
        bool correct1 = read_image (image1, &I1, &nx, &ny);
        bool correct2 = read_image (image2, &I2, &nx1, &ny1);

        // if the images are correct, compute the optical flow
        if (correct1 && correct2 && nx == nx1 && ny == ny1)
        {
            //set the number of scales according to the size of the
            //images.  The value N is computed to assure that the smaller
            //images of the pyramid don't have a size smaller than 16x16
            const float N =
                1 + log (std::min (nx, ny) / 16.) / log (1. / zfactor);
            if ((int) N < nscales)
                nscales = (int) N;

            cout << endl << "ncores:" << nproc << " alpha:" << alpha
                 << " gamma:" << gamma << " scales:" << nscales << " nu:" <<
                 zfactor << " TOL:" << TOL << " inner:" << initer << " outer:" <<
                 outiter << endl;

            //allocate memory for the flow
            float *u1 = new float[nx * ny];
            float *v1 = new float[nx * ny];

            //compute the optic flow
            brox_optic_flow (I1, I2, u1, v1, nx, ny, alpha, gamma,
                             nscales, zfactor, TOL, initer, outiter, verbose);
            if (test_lr)
            {
                cout << "Brox method with left / right test\n" << endl;
                // allocate memory for the flow
                float *u2 = new float[nx * ny];
                float *v2 = new float[nx * ny];

                // compute the optic flow between I_right and I_left
                brox_optic_flow (I2, I1, u2, v2, nx, ny, alpha, gamma,
                                 nscales, zfactor, TOL, initer, outiter, verbose);

                // determination of the occulted pixels of the left disparity map
                set_nan_disparities(u1, u2, seuil, nx, ny);

                delete[]u2;
                delete[]v2;
            }
            else
            {
                cout << "Brox method without left / right test\n" << endl;
            }

            iio_write_image_float ((char *) outfile, u1, nx, ny);


            //free dynamic memory
            free (I1);
            free (I2);
            delete[]u1;
            delete[]v1;

        }
        else
            cerr <<
                 "Cannot read the images or the size of the images are not equal" <<
                 endl;
    }

    return 0;
}
