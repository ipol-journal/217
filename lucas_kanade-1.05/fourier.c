
/*
* fourier.c
*
* Copyright (C) 2016, Tristan Dagobert, CMLA, École Normale Supérieure de Cachan.
*
* This software is a computer program.[describe
* functionalities and technical features of your software].
*
* This software is governed by the CeCILL-C license under French law and
* abiding by the rules of distribution of free software.  You can  use,
* modify and/ or redistribute the software under the terms of the CeCILL-C
* license as circulated by CEA, CNRS and INRIA at the following URL
* "http://www.cecill.info".
*
* As a counterpart to the access to the source code and  rights to copy,
* modify and redistribute granted by the license, users are provided only
* with a limited warranty  and the software's author,  the holder of the
* economic rights,  and the successive licensors  have only  limited
* liability.
*
* In this respect, the user's attention is drawn to the risks associated
* with loading,  using,  modifying and/or developing or reproducing the
* software by the user in light of its specific status of free software,
* that may mean  that it is complicated to manipulate,  and  that  also
* therefore means  that it is reserved for developers  and  experienced
* professionals having in-depth computer knowledge. Users are therefore
* encouraged to load and test the software's suitability as regards their
* requirements in conditions enabling the security of their systems and/or
* data to be ensured and,  more generally, to use and operate it in the
* same conditions as regards security.
*
* The fact that you are presently reading this means that you have had
* knowledge of the CeCILL-C license and that you accept its terms.
*
*/


#include "fourier.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/*===========================================================================*/
void
fft2 (fftw_complex ** output, void *input, int nlig, int ncol, int type)
{
    fftw_complex *spatial;
    fftw_complex *frequency;
    fftw_plan plan;

    float *f;
    double *d;
    long int i;
    long int npix;

    npix = nlig * ncol;

    /* allocation of the output image */
    frequency = *output;
    if (frequency == NULL)
    {
        frequency = (fftw_complex *) fftw_malloc (npix * sizeof (fftw_complex));
        if (frequency == NULL)
        {
            perror ("fft2()\n");
            exit (EXIT_FAILURE);
        }
        *output = frequency;
    }

    /* depending on the input image */
    switch (type)
    {
    case FLOAT_INPUT:
    {
        f = (float *) input;
        spatial = fftw_malloc (npix * sizeof (fftw_complex));

        if (spatial == NULL)
        {
            perror ("fft2()\n");
            exit (EXIT_FAILURE);
        }

        /* filling the structure which will be used by fftw */
        for (i = 0; i < npix; i++)
        {

            spatial[i][0] = (double) f[i];
            spatial[i][1] = 0;
        }
        break;
    }
    case DOUBLE_INPUT:
    {
        d = (double *) input;
        spatial = fftw_malloc (npix * sizeof (fftw_complex));

        if (spatial == NULL)
        {
            perror ("fft2()\n");
            exit (EXIT_FAILURE);
        }

        /* filling the structure which will be used by fftw */
        for (i = 0; i < npix; i++)
        {

            spatial[i][0] = d[i];
            spatial[i][1] = 0;
        }
        break;
    }

    case FFTW_COMPLEX_INPUT:
    {
        spatial = (fftw_complex *) input;
        break;
    }
    default:
    {
        perror ("fft2()\n");
        exit (EXIT_FAILURE);
    }
    }
    /* computation of the execution plan */
    plan =
        fftw_plan_dft_2d (nlig, ncol, spatial, frequency, FFTW_FORWARD,
                          FFTW_ESTIMATE);

    /* Fourier transform */
    fftw_execute (plan);

    /* structures destruction */
    fftw_destroy_plan (plan);

    if (type == FLOAT_INPUT || type == DOUBLE_INPUT)
    {
        fftw_free (spatial);
    }

    return;
}

/*===========================================================================*/
void
ifft2 (fftw_complex ** output, void *input, int nlig, int ncol, int type)
{
    fftw_complex *spatial;
    fftw_complex *frequency;
    fftw_plan plan;

    long int npix;
    long int i;
    float *f;

    npix = nlig * ncol;

    /* depending on the input image */
    switch (type)
    {
    case FLOAT_INPUT:
    {
        f = (float *) input;
        frequency = fftw_malloc (npix * sizeof (fftw_complex));

        if (frequency == NULL)
        {
            perror ("ifft2()\n");
            exit (EXIT_FAILURE);
        }

        /* filling the structure which will be used by fftw */
        for (i = 0; i < npix; i++)
        {

            frequency[i][0] = (double) f[i];
            frequency[i][1] = 0;
        }
        break;
    }

    case FFTW_COMPLEX_INPUT:
    {
        frequency = (fftw_complex *) input;
        break;
    }
    default:
    {
        perror ("ifft2()\n");
        exit (EXIT_FAILURE);
    }
    }


    /* allocation of output image */
    spatial = *output;
    if (spatial == NULL)
    {
        spatial = fftw_malloc (npix * sizeof (fftw_complex));
        if (spatial == NULL)
        {
            perror ("fft2()\n");
            exit (EXIT_FAILURE);
        }
        *output = spatial;
    }

    /* computation of the execution plan */
    plan =
        fftw_plan_dft_2d (nlig, ncol, frequency, spatial, FFTW_BACKWARD,
                          FFTW_ESTIMATE);

    /* inverse Fourier transform */
    fftw_execute (plan);

    /* suitable normalization */
    for (i = 0; i < npix; i++)
    {
        spatial[i][0] /= npix;
        spatial[i][1] /= npix;
    }

    /* structures destruction */
    fftw_destroy_plan (plan);

    if (type == FLOAT_INPUT)
    {
        fftw_free (frequency);
    }

    return;
}

/*===========================================================================*/
void
fftshift (void *input, int nlig, int ncol, int type)
{

    switch (type)
    {
    case FFTW_COMPLEX_INPUT:
    {
        fftw_complex *temporaire;
        fftw_complex *t2 = NULL;
        fftw_complex *A, *neoA, *B;
        fftw_complex *fw_input = (fftw_complex *) input;

        int i, ncolp1, nligp1;
        int nboctets, nboctetsplus1;

        temporaire =
            (fftw_complex *) fftw_malloc (ncol * sizeof (fftw_complex));
        if (temporaire == NULL)
        {
            perror ("fftshift()\n");
            exit (EXIT_FAILURE);
        }

        /* Be an image I such that
           .   A  B                      C  D
           I = D  C   then fftshift(I) = B  A
         */
        nboctets = (ncol / 2) * sizeof (fftw_complex);
        ncolp1 = round (ncol / 2.0);
        nboctetsplus1 = ncolp1 * sizeof (fftw_complex);

        for (i = 0; i < nlig; i++)
        {
            // A <--> B
            // D <--> C
            A = &(fw_input[i * ncol]);
            neoA = &(fw_input[i * ncol + ncol - ncolp1]);

            B = &(fw_input[i * ncol + ncolp1]);
            memcpy (temporaire, A, nboctetsplus1);
            memcpy (A, B, nboctets);
            memcpy (neoA, temporaire, nboctetsplus1);
        }

        nboctets = ncol * sizeof (fftw_complex);
        nligp1 = round (nlig / 2.0);
        if (nlig % 2 == 1)
        {
            t2 = (fftw_complex *) fftw_malloc (ncol * sizeof (fftw_complex));
            if (t2 == NULL)
            {
                perror ("fftshift()\n");
                exit (EXIT_FAILURE);
            }
            memcpy (t2, &(fw_input[(nlig / 2) * ncol]), nboctets);
        }
        for (i = 0; i < nlig / 2; i++)
        {
            /* B A */
            /* ^ ^ */
            /* | | */
            /* V V */
            /* C D */
            A = &(fw_input[i * ncol]);
            neoA = &(fw_input[(i + (nlig - nligp1)) * ncol]);
            //C = &(input[(i+nlig/2)*ncol + ncol/2]);
            B = &(fw_input[(i + nligp1) * ncol]);
            memcpy (temporaire, A, nboctets);
            memcpy (A, B, nboctets);
            memcpy (neoA, temporaire, nboctets);
        }
        if (nlig % 2 == 1)
        {
            B = &(fw_input[(nlig - 1) * ncol]);
            memcpy (B, t2, nboctets);
            fftw_free (t2);
        }

        fftw_free (temporaire);
        break;
    }
    case FLOAT_INPUT:
    {
        float *t1;
        float *t3 = NULL;
        float *a, *neoa, *b;	// *C *B, *D;
        float *f_input = (float *) input;
        int i, ncolp1, nligp1;
        int nboctets, nboctetsplus1;

        t1 = (float *) malloc (ncol * sizeof (float));
        if (t1 == NULL)
        {
            perror ("fftshift()\n");
            exit (EXIT_FAILURE);
        }

        /* Be an image I such that
           .   A  B                      C  D
           I = D  C   then fftshift(I) = B  A
         */
        nboctets = (ncol / 2) * sizeof (float);
        ncolp1 = round (ncol / 2.0);
        nboctetsplus1 = ncolp1 * sizeof (float);

        for (i = 0; i < nlig; i++)
        {
            // A <--> B
            // D <--> C
            a = &(f_input[i * ncol]);
            neoa = &(f_input[i * ncol + ncol - ncolp1]);

            b = &(f_input[i * ncol + ncolp1]);
            memcpy (t1, a, nboctetsplus1);
            memcpy (a, b, nboctets);
            memcpy (neoa, t1, nboctetsplus1);
        }

        nboctets = ncol * sizeof (float);
        nligp1 = round (nlig / 2.0);
        if (nlig % 2 == 1)
        {
            t3 = (float *) malloc (ncol * sizeof (float));
            if (t3 == NULL)
            {
                perror ("fftshift()\n");
                exit (EXIT_FAILURE);
            }
            memcpy (t3, &(f_input[(nlig / 2) * ncol]), nboctets);
        }
        for (i = 0; i < nlig / 2; i++)
        {
            /* B A */
            /* ^ ^ */
            /* | | */
            /* V V */
            /* C D */
            a = &(f_input[i * ncol]);
            neoa = &(f_input[(i + (nlig - nligp1)) * ncol]);
            b = &(f_input[(i + nligp1) * ncol]);
            memcpy (t1, a, nboctets);
            memcpy (a, b, nboctets);
            memcpy (neoa, t1, nboctets);
        }
        if (nlig % 2 == 1)
        {
            b = &(f_input[(nlig - 1) * ncol]);
            memcpy (b, t3, nboctets);
            free (t3);
        }

        free (t1);

        break;
    }
    default:
    {
        perror ("fftshift()\n");
        exit (EXIT_FAILURE);
    }
    }

    return;
}

/*===========================================================================*/
void
modulus (float **output, fftw_complex * input, int npix)
{
    int i;
    float *im_modulus;
    im_modulus = *output;

    /* allocation of the modulus image */
    if (im_modulus == NULL)
    {
        im_modulus = (float *) malloc (npix * sizeof (float));
        if (im_modulus == NULL)
        {
            perror ("modulus()\n");
            exit (EXIT_FAILURE);
        }
        *output = im_modulus;
    }
    for (i = 0; i < npix; i++)
    {
        im_modulus[i] = (float)
                        sqrt (input[i][0] * input[i][0] + input[i][1] * input[i][1]);
    }
    return;
}

/*===========================================================================*/
void
reel (double **neo, fftw_complex * ante, int npix)
{
    int i;
    double *im = *neo;

    /* allocation of the modulus image */
    if (im == NULL)
    {
        im = (double *) malloc (npix * sizeof (double));
        if (im == NULL)
        {
            perror ("reel()\n");
            exit (EXIT_FAILURE);
        }
        *neo = im;
    }
    for (i = 0; i < npix; i++)
    {
        im[i] = ante[i][0];
    }
    return;
}

/*===========================================================================*/
