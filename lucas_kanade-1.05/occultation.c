
/*
* occultation.c
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


#include "occultation.h"
#include <math.h>
#include <stdio.h>
/*===========================================================================*/
#define ABSOLUTE_TEST 1
/*===========================================================================*/
void
setting_nan_disparities (double *disp_s_left, double *disp_s_right,
                         double seuil, int test_type, int ncol, int nlig)
{
    int l, c;

    /* computation of bilinear interpolation coefficients */
    cerce_t *cerces;
    //cerces = computing_bilinear_splines(disp_s_right, nlig, ncol, 1);
    cerces =
        compute_splines (disp_s_right, nlig, ncol, 1, gsl_interp2d_bilinear);

    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;

    spline = cerces[0].spline;
    xacc = cerces[0].xacc;
    yacc = cerces[0].yacc;

    /* computation of the right interpolated disparity map */
    double *disp_interp_right;
    disp_interp_right = (double *) malloc (nlig * ncol * sizeof (double));

    double x, y;
    int xmin, xmax;

    xmin = 0;
    xmax = ncol - 1;

    for (l = 0; l < nlig; l++)
    {
        for (c = 0; c < ncol; c++)
        {
            x = c + disp_s_left[l * ncol + c];
            y = l;
            if (xmin <= x && x <= xmax)
            {
                disp_interp_right[l * ncol + c] =
                    gsl_spline2d_eval (spline, x, y, xacc, yacc);
            }
            else
            {
                disp_interp_right[l * ncol + c] = nan ("");
            }
        }
    }

    /* computation of the relative difference */
    double ecart, moy;

    if (test_type == ABSOLUTE_TEST)
    {
        printf ("ICI %f !!\n", seuil);
        for (l = 0; l < nlig; l++)
        {
            for (c = 0; c < ncol; c++)
            {
                ecart =
                    fabs (disp_s_left[l * ncol + c] +
                          disp_interp_right[l * ncol + c]);

                if (ecart > seuil || isnan (ecart))
                {
                    disp_s_left[l * ncol + c] = nanf ("");
                }
            }
        }
    }
    else
    {
        /* computation of the relative difference */
        for (l = 0; l < nlig; l++)
        {
            for (c = 0; c < ncol; c++)
            {
                ecart =
                    fabs (disp_s_left[l * ncol + c] +
                          disp_interp_right[l * ncol + c]);
                moy =
                    fabs (disp_s_left[l * ncol + c] -
                          disp_interp_right[l * ncol + c]) / 2.0;
                if (ecart / moy > seuil || isnan (ecart))
                {
                    disp_s_left[l * ncol + c] = nanf ("");
                }
            }
        }
    }

    free (disp_interp_right);
    free_splines (cerces, 1);
}

/*===========================================================================*/
