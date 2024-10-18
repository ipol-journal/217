
/*
* interpolations.h
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


#if !(defined INTERPOLATIONS_H)
#define INTERPOLATIONS_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/*===========================================================================*/
struct cerce_s
{
	double * im;
	int nlig;
	int ncol;
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    int xmin;
    int xmax;
};
typedef struct cerce_s cerce_t;
/*===========================================================================*/
cerce_t *compute_splines (double *im, int nlig, int ncol, int ncan,
                            const gsl_interp2d_type * T);
void free_splines (cerce_t * cerces, int ncan);
/*===========================================================================*/
#endif
