
/*
* interpolations.c
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


#include "interpolations.h"
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
/*===========================================================================*/
cerce_t *
compute_splines (double *im, int nlig, int ncol, int ncan,
                   const gsl_interp2d_type * T)
{
    cerce_t *cerces = NULL;
    cerces = (cerce_t *) malloc (ncan * sizeof (cerce_t));

    int n, l, c;
    double *xa = (double *) malloc (ncol * sizeof (double));
    double *ya = (double *) malloc (nlig * sizeof (double));
    double *za = NULL;
    for (c = 0; c < ncol; c++)
    {
        xa[c] = 1.0 * c;
    }
    for (l = 0; l < nlig; l++)
    {
        ya[l] = 1.0 * l;
    }

    for (n = 0; n < ncan; n++)
    {
        /* ad hoc allocation */
		cerces[n].im = malloc (ncol * nlig * sizeof (double));
		memcpy(cerces[n].im, im + (ncol * nlig * n), ncol * nlig * sizeof(double));
		cerces[n].nlig = nlig;
		cerces[n].ncol = ncol;
        cerces[n].spline = gsl_spline2d_alloc (T, ncol, nlig);
        cerces[n].xacc = gsl_interp_accel_alloc ();
        cerces[n].yacc = gsl_interp_accel_alloc ();
        cerces[n].xmin = 0;
        cerces[n].xmax = ncol - 1;

        /* ad hoc initialization */
        za = im + ncol * nlig * n;
        gsl_spline2d_init (cerces[n].spline, xa, ya, za, ncol, nlig);
    }

    free (xa);
    free (ya);

    return cerces;
}

/*===========================================================================*/
void
free_splines (cerce_t * cerces, int ncan)
{
    int n;

    for (n = 0; n < ncan; n++)
    {
		free(cerces[n].im);
        gsl_spline2d_free (cerces[n].spline);
        gsl_interp_accel_free (cerces[n].xacc);
        gsl_interp_accel_free (cerces[n].yacc);
    }
    free (cerces);
    cerces = NULL;
    return;
}

/*===========================================================================*/
