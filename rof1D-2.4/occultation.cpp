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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <fftw3.h>

#include "iio.h"
#include "occultation.h"
#include "bicubic_interpolation.h"

#define ABSOLUTE_TEST 1

/*===========================================================================*/
void set_nan_disparities(float * map_disp_left, float * map_disp_right,
                         double seuil, int ncol, int nlig)
{
    int l, c;
    float disp_left, disp_interp_right, neo_x, neo_y, err, moy;
    for(l=0; l< nlig; l++)
    {
        for(c=0; c< ncol; c++)
        {
            /* disparity value at pixel (l, c) of the left disparity map */
            disp_left = map_disp_left[l * ncol + c];

            if (!isfinite(disp_left))
            {
                map_disp_left[l*ncol + c] = nanf("");
                continue;
            }

            /* subpixel location Q in the right image */
            neo_x = c + disp_left;
            neo_y = l;

            /* computation of the disparity at location Q of the right map */
            disp_interp_right = bicubic_interpolation (map_disp_right, neo_x,
                                neo_y, ncol, nlig);

            /* bias (there is a « + » because terms have opposite sens…) */
            err = abs(disp_left + disp_interp_right);
            moy = abs((disp_left - disp_interp_right) / 2.0);

            /* computation of the absolute difference */
            if (ABSOLUTE_TEST == 1)
            {
                if (err > seuil)
                {
                    map_disp_left[l*ncol + c] = nanf("");
                }
            }
            /* computation of the relative difference */
            else
            {
                if (err > seuil * moy)
                {
                    map_disp_left[l*ncol + c] = nanf("");
                }
            }
        }
    }
    return;
}
