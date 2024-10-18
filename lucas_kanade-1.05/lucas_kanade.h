
/*
* lucas_kanade.h
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


#if !(defined LUCAS_KANADE_H)
#define LUCAS_KANADE_H
/*===========================================================================*/
#include "interpolations.h"
#include "fourier.h"
#include "iio.h"
#include "occultation.h"

#if defined(WITH_OPENMP)
#include <pthread.h>
#endif
/*===========================================================================*/
struct pyramide_s
{
    double **images;
    int *ncol;
    int *nlig;
    int ncan;
};
typedef struct pyramide_s pyramide_t;

/*===========================================================================*/
struct param_pyr_s
{
#if defined(WITH_OPENMP)
    pthread_mutex_t *mutex;
#endif
    int tid;
    double *im;
    int ncol;
    int nlig;
    int ncan;
    int nb_scales;
    pyramide_t *pyr;
};
typedef struct param_pyr_s param_pyr_t;
/*===========================================================================*/
struct parametres_s
{
#if defined(WITH_OPENMP)
    pthread_mutex_t *mutex;
#endif
    int tid;
    int *indice;

    cerce_t *cerces;
    double *im_l;
    double *im_r;
    double *u_gradx;
    double *dh_init_p;
    double *dispx;

    double *noyau;
    int rayon_noy;
    int diam_noy;

    int nlig;
    int ncol;
    int ncan;

    int half_vcol;
    int half_vlig;
    int nb_iter;
    double sigma_dv;
};
typedef struct parametres_s parametres_t;
/*===========================================================================*/
double *lucas_kanade (pyramide_t * pyramide_g, pyramide_t * pyramide_d,
                      int demi_vx, int demi_vy, int nb_procs,
                      int nb_iter, float sigma_dv, int nb_scales);

void *build_pyramid (pyramide_t ** pyramide_g, pyramide_t ** pyramide_d,
                        float *fim_l, float *fim_r, int nlig, int ncol,
                        int ncan, int nb_scales);
void *build_pyramid_thread (void *args);
void free_pyramid (pyramide_t * pyr, int nb_scales);

double *redimension (double *im, int nlig, int ncol, int neo_nlig,
                         int neo_ncol, int ncan,
                         const gsl_interp2d_type * type);

//void *estimating_shift_thread (void *args);
void *
estimate_shift_thread (int indice, cerce_t *v_cerces, double * im_l, double * im_r,
						 double *u_gradx, double *dispx_init, double *dispx,
						 int nlig, int ncol, int ncan, int half_vcol, int half_vlig,
						 int nb_iter);
double *compute_gradient_x (double *im_l, int nlig, int ncol, int ncan);
double *convolve_with_gaussian (double *im, int nlig, int ncol, int ncan,
                                  double sigma);
double estimate_translation (double *bloc_u, double *bloc_Du,
                               cerce_t * v_cerces, int ligdeb, int coldeb,
                               int nb_iter, double *bloc_dispx_init, int ncol,
                               int nlig, int ncan, double *bloc_v_k,
                               double *bloc_dt);
void shift_block (double *bloc_v_k, double *bloc_dh_init, double dh,
                     cerce_t * v_cerces, int ligdeb, int coldeb, int ncol,
                     int nlig, int ncan);

void write_image (char *nom, double *d, int nlig, int ncol, int ncan);

int algorithm (char *ficim_l, char *ficim_r, char *fic_disp, int half_vcol,
               int half_vlig, int nb_iter, int nb_scales, double sigma_dv,
               long int nb_procs, double seuil, int type_test);
/*===========================================================================*/
#endif
