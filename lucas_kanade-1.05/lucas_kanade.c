
/*
* lucas_kanade.c
*
* Copyright (C) 2016, Tristan Dagobert, CMLA, École Normale Supérieure de Cachan.
*
* This software is a computer program.
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

#include "lucas_kanade.h"
#include "bicubic_interpolation.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(WITH_OPENMP)
#include <pthread.h>
#endif
/*===========================================================================*/
void *
estimate_shift_thread (int indice, cerce_t *v_cerces, double * im_l, double * im_r,
						 double *u_gradx, double *dispx_init, double *dispx,
						 int nlig, int ncol, int ncan, int half_vcol, int half_vlig,
						 int nb_iter)
{
	(void) im_r;

    int i, j;
    int n, l, c, ii;

    int ligdeb;
    int coldeb;
    int ligfin;
    int colfin;

    double dh;

    /* blocks allocation */
    double *bloc_u, *bloc_u_dx, *bloc_dispx_init;
    double *bloc_v_k, *bloc_dt;
    int blig, bcol;

    /* maximal dimension of neighborhood */
    blig = 2 * half_vlig + 1;
    bcol = 2 * half_vcol + 1;

    bloc_u = (double *) malloc (ncan * blig * bcol * sizeof (double));
    bloc_u_dx = (double *) malloc (ncan * blig * bcol * sizeof (double));
    bloc_dispx_init = (double *) malloc (blig * bcol * sizeof (double));

    bloc_v_k = (double *) malloc (bcol * blig * ncan * sizeof (double));
    bloc_dt = (double *) malloc (bcol * blig * ncan * sizeof (double));

    if (bloc_u == NULL || bloc_u_dx == NULL || bloc_dispx_init == NULL
		|| bloc_v_k == NULL || bloc_dt == NULL)
    {
        perror ("estimate_shift_thread()\n");
        exit (EXIT_FAILURE);
    }

	i = (indice) / ncol;
	j = (indice) % ncol;

	/* beginning of the treatment... */

	/* neigborhood dimension */
	ligdeb = i - half_vlig;
	coldeb = j - half_vcol;
	ligfin = i + half_vlig;
	colfin = j + half_vcol;

	/* rectification of the neigborhood size when close to the bound */
	ligdeb = (ligdeb > 0) ? ligdeb : 0;
	coldeb = (coldeb > 0) ? coldeb : 0;

	ligfin = (ligfin < nlig - 1) ? ligfin : nlig - 1;
	colfin = (colfin < ncol - 1) ? colfin : ncol - 1;

	blig = ligfin - ligdeb + 1;
	bcol = colfin - coldeb + 1;

	/* blocks selection */
	for (l = ligdeb, ii = 0; l <= ligfin; l++)
	{
		for (c = coldeb; c <= colfin; c++, ii++)
		{
			bloc_dispx_init[ii] = dispx_init[l * ncol + c];
		}
	}

	for (n = 0, ii = 0; n < ncan; n++)
	{
		for (l = ligdeb; l <= ligfin; l++)
		{
			for (c = coldeb; c <= colfin; c++, ii++)
			{
				bloc_u[ii] = im_l[ncol * nlig * n + l * ncol + c];

				bloc_u_dx[ii] = u_gradx[ncol * nlig * n + l * ncol + c];
			}
		}
	}

	/* estimation of displacement dh */
	dh = estimate_translation (bloc_u, bloc_u_dx, v_cerces,
								 ligdeb, coldeb, nb_iter, bloc_dispx_init,
								 bcol, blig, ncan, bloc_v_k, bloc_dt);

	dispx[i * ncol + j] = dh;

	free (bloc_u);
	free (bloc_u_dx);
	free (bloc_dispx_init);
	free (bloc_v_k);
	free (bloc_dt);

    return (NULL);
}

/*===========================================================================*/
double *
compute_gradient_x (double *im, int nlig, int ncol, int ncan)
{
    double *gradx = NULL;
    int n, l, c;
    int pos;
    gradx = (double *) calloc (nlig * ncol * ncan, sizeof (double));
    if (gradx == NULL)
    {
        perror ("compute_gradient_x()\n");
        exit (EXIT_FAILURE);
    }

    for (n = 0; n < ncan; n++)
    {
        for (l = 0; l < nlig; l++)
        {
            pos = ncol * nlig * n + l * ncol;

            for (c = 1; c < ncol - 1; c++)
            {
                gradx[pos + c] = (im[pos + c + 1] - im[pos + c - 1]) / 2.0;
            }
        }
    }

    return gradx;
}

/*===========================================================================*/
double *
convolve_with_gaussian (double *im, int nlig, int ncol, int ncan,
                          double sigma)
{
    fftw_complex *G = NULL;
    fftw_complex *tfim = NULL;
    fftw_complex *F = NULL;
    fftw_complex *f = NULL;
    double *g = NULL;
    double *neo = NULL;
    double *neo_x4 = NULL;
    double *im_x4 = NULL;

    int x, y;
    int n, c, l;
    long int i, npix;
    double deux_sigma2 = 2 * sigma * sigma;
    double val;

    /* mirror symetrization to avoid rebound effects */
    int xcol = 2 * ncol;
    int xlig = 2 * nlig;

    npix = xlig * xcol;

    g = (double *) malloc (npix * sizeof (double));
    neo = (double *) malloc (nlig * ncol * ncan * sizeof (double));
    neo_x4 = (double *) malloc (npix * sizeof (double));
    im_x4 = (double *) malloc (xlig * xcol * sizeof (double));
    F = (fftw_complex *) malloc (xlig * xcol * sizeof (fftw_complex));

    if (im_x4 == NULL || neo == NULL || neo_x4 == NULL || F == NULL
            || g == NULL)
    {
        perror ("convolve_with_gaussian()\n");
        exit (EXIT_FAILURE);
    }


    /* kernel definition */
    double somme = 0.0;

    for (l = 0, x = -xlig / 2; l < xlig; l++, x++)
    {
        for (c = 0, y = -xcol / 2; c < xcol; c++, y++)
        {
            g[l * xcol + c] =
                exp (-(x * x + y * y) / deux_sigma2) / (deux_sigma2 * M_PI);
            somme += g[l * xcol + c];
        }
    }
    /* do not forget normalization ! */
    for (i = 0; i < npix; i++)
    {
        g[i] /= somme;
    }

    /* then Fourier transform... */
    fft2 (&G, (void *) g, xlig, xcol, DOUBLE_INPUT);

    /* convolution... */
    for (n = 0; n < ncan; n++)
    {

        for (l = 0; l < nlig; l++)
        {
            for (c = 0; c < ncol; c++)
            {
                val = im[ncol * nlig * n + l * ncol + c];
                im_x4[l * xcol + c] = val;
                im_x4[l * xcol + (xcol - 1 - c)] = val;
                im_x4[(xlig - 1 - l) * xcol + c] = val;
                im_x4[(xlig - 1 - l) * xcol + (xcol - 1 - c)] = val;
            }
        }

        fft2 (&tfim, (void *) im_x4, xlig, xcol, DOUBLE_INPUT);
        for (i = 0; i < npix; i++)
        {
            /* (a + ib) x (c + id) = (a*c - b*d) + i(a*d + b*c) */
            F[i][0] = G[i][0] * tfim[i][0] - G[i][1] * tfim[i][1];
            F[i][1] = G[i][0] * tfim[i][1] - G[i][1] * tfim[i][0];
        }

        ifft2 (&f, (void *) F, xlig, xcol, FFTW_COMPLEX_INPUT);
        fftshift ((void *) f, xlig, xcol, FFTW_COMPLEX_INPUT);
        reel (&g, f, npix);

        /* copy in the suitable Panama's channel */
        memcpy (neo_x4, g, npix * sizeof (double));

        /* pick up the useful data */
        for (l = 0; l < nlig; l++)
        {
            for (c = 0; c < ncol; c++)
            {
                val = neo_x4[l * xcol + c];
                neo[ncol * nlig * n + l * ncol + c] = val;
            }
        }
    }


    free (neo_x4);
    free (im_x4);
    free (g);

    fftw_free (G);
    fftw_free (tfim);
    fftw_free (F);
    fftw_free (f);

    return neo;
}


/*===========================================================================*/
double
estimate_translation (double *bloc_u, double *bloc_Du, cerce_t * v_cerces,
                        int ligdeb, int coldeb, int nb_iter,
                        double *bloc_dh_init, int ncol, int nlig, int ncan,
                        double *bloc_v_k, double *bloc_dt)
{
    double dh = 0.0;

    int ii, n;
    double energie_k, energie_k_minus_1;
    double d, nume, denom;

    for (ii = 0; ii < nb_iter; ii++)
    {
        /* image shifting */
        shift_block (bloc_v_k, bloc_dh_init, dh, v_cerces,
                        ligdeb, coldeb, ncol, nlig, ncan);

        /* temporal gradient */
        for (n = 0; n < ncan * ncol * nlig; n++)
        {
            bloc_dt[n] = bloc_v_k[n] - bloc_u[n];
        }

        /* energy computation */
        energie_k = 0.0;
        for (n = 0; n < ncan * ncol * nlig; n++)
        {
            energie_k += bloc_dt[n] * bloc_dt[n];
        }
        if (ii == 0)
        {
            energie_k_minus_1 = energie_k;
        }
        else
        {
            /* verify energy's decreasing */
            if (energie_k <= energie_k_minus_1)
            {
                energie_k_minus_1 = energie_k;
            }
            else
            {
                /* if one diverges, the we go out, cancelling the
                 * previous shift */
                dh = dh - d;
                break;
            }
        }

        /* displacement estimation */
        d = 0.0;
        nume = 0.0;
        denom = 0.0;

        for (n = 0; n < ncan * ncol * nlig; n++)
        {
            nume += bloc_Du[n] * bloc_dt[n];
            denom += bloc_Du[n] * bloc_Du[n];
        }
        if (nume == 0.0 && denom == 0.0)
        {
            d = 0;
        }
        else
        {
            d = nume / denom;
        }

        dh = dh + d;
    }

    return dh;
}

/*===========================================================================*/
void
shift_block (double *bloc_v_k, double *bloc_dh_init, double dh,
                cerce_t * v_cerces, int ligdeb, int coldeb, int ncol,
                int nlig, int ncan)
{
    int c, l, n;
    double x, y;
    double dh_init;
    double z;

    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;

    for (l = 0; l < nlig; l++)
    {
        y = ligdeb + l;
        for (c = 0; c < ncol; c++)
        {
            /* compute the interpolate value */
            dh_init = bloc_dh_init[l * ncol + c];
            x = coldeb + c - dh_init - dh;

            for (n = 0; n < ncan; n++)
            {
                /* pick up the coefficients of the spline interpolation */
                spline = v_cerces[n].spline;
                xacc = v_cerces[n].xacc;
                yacc = v_cerces[n].yacc;
				(void) spline; (void) xacc; (void) yacc;
				
#if defined(WITH_BICUBIC)				
				z = bicubic_interp_point (v_cerces[n].im, x, y,
										  v_cerces[n].ncol, v_cerces[n].nlig, 1,
										  0, 1);
#else									  
                if (v_cerces[n].xmin <= x && x <= v_cerces[n].xmax)
                {
                    z = gsl_spline2d_eval (spline, x, y, xacc, yacc);
                }
                else
                {
                    z = 0.0;
                }
#endif				
                bloc_v_k[ncol * nlig * n + l * ncol + c] = z;
            }
        }
    }

    return;
}

/*===========================================================================*/
void
free_pyramid (pyramide_t * pyr, int nb_scales)
{
    int n;

    for (n = 0; n < nb_scales; n++)
    {
        free (pyr->images[n]);
    }
    free (pyr->images);
    free (pyr->nlig);
    free (pyr->ncol);
    free (pyr);
    pyr = NULL;
    return;
}

/*===========================================================================*/
void *
build_pyramid (pyramide_t ** pyramide_g, pyramide_t ** pyramide_d,
                  float *fim_l, float *fim_r, int nlig, int ncol, int ncan,
                  int nb_scales)
{
    long int i;

    param_pyr_t *param_pyr;
    double *im_l, *im_r;

    im_l = (double *) malloc (ncol * nlig * ncan * sizeof (double));
    im_r = (double *) malloc (ncol * nlig * ncan * sizeof (double));
    if (im_l == NULL || im_r == NULL)
    {
        printf ("[FATAL] Allocation failed\n");
		exit(EXIT_FAILURE);
    }

    for (i = 0; i < ncol * nlig * ncan; i++)
    {
        im_l[i] = (double) fim_l[i];
        im_r[i] = (double) fim_r[i];
    }

    /* encapsulating parameters */
    param_pyr = (param_pyr_t *) malloc (2 * sizeof (param_pyr_t));

#if defined(WITH_OPENMP)	
	pthread_mutex_t bille = PTHREAD_MUTEX_INITIALIZER;
#endif

#if defined(WITH_OPENMP)	
	param_pyr[0].mutex = &bille;
#endif	
    param_pyr[0].tid = 0;
    param_pyr[0].im = im_l;
    param_pyr[0].ncol = ncol;
    param_pyr[0].nlig = nlig;
    param_pyr[0].ncan = ncan;
    param_pyr[0].nb_scales = nb_scales;
    param_pyr[0].pyr = NULL;

#if defined(WITH_OPENMP)	
	param_pyr[1].mutex = &bille;
#endif	
    param_pyr[1].tid = 1;
    param_pyr[1].im = im_r;
    param_pyr[1].ncol = ncol;
    param_pyr[1].nlig = nlig;
    param_pyr[1].ncan = ncan;
    param_pyr[1].nb_scales = nb_scales;
    param_pyr[1].pyr = NULL;

#if defined(WITH_OPENMP)
#ifdef _OPENMP	
    #pragma omp parallel for
#endif			
#endif		
	for (int ii=0; ii < 2; ii++)
	{
		build_pyramid_thread ((void *) (param_pyr + ii));
	}
	
    /* pick up the pyramids */
    *pyramide_g = param_pyr[0].pyr;
    *pyramide_d = param_pyr[1].pyr;
    free (param_pyr);

    free (im_l);
    free (im_r);

    return NULL;
}

/*===========================================================================*/
void *
build_pyramid_thread (void *args)
{
    printf ("[INFO] building pyramid...\n");
    param_pyr_t *p = (param_pyr_t *) args;
    int tid = p->tid;
    double *im = p->im;
    int ncol = p->ncol;
    int nlig = p->nlig;
    int ncan = p->ncan;
    int nb_scales = p->nb_scales;

    pyramide_t *kheops;

    double **images;
    int *ncols;
    int *nligs;

    kheops = (pyramide_t *) malloc (sizeof (pyramide_t));
    images = (double **) malloc (nb_scales * sizeof (double *));
    nligs = (int *) malloc (nb_scales * sizeof (int));
    ncols = (int *) malloc (nb_scales * sizeof (int));

    if (kheops == NULL || nligs == NULL || ncols == NULL || images == NULL)
    {
        perror ("build_pyramid_thread()\n");
        exit (EXIT_FAILURE);
    }

    p->pyr = kheops;
    kheops->ncan = ncan;
    kheops->nlig = nligs;
    kheops->ncol = ncols;
    kheops->images = images;

    /* do not forget the original image */
    nligs[0] = nlig;
    ncols[0] = ncol;
    images[0] = malloc (ncol * nlig * ncan * sizeof (double));
    if (images[0] == NULL)
    {
        perror ("build_pyramid_thread()\n");
        exit (EXIT_FAILURE);
    }

    memcpy (images[0], im, ncol * nlig * ncan * sizeof (double));

    int e = 0;
    double *ante;

    printf ("[OK] thread %d building at scale %d image of size %d x %d\n", tid,
            e, ncol, nlig);

    for (e = 1; e < nb_scales; e++)
    {
        nligs[e] = nligs[e - 1] / 2;
        ncols[e] = ncols[e - 1] / 2;

		/* One has to use a mutex here because the FFTW3 plan is a structure
		 * which cannot be call by several threads simultaneously. If not
		 * it induces locking and segmentation fault… 
		 */
		/* image convolution at previous scale with sigma=1.6 */
#if defined(WITH_OPENMP)
		pthread_mutex_lock (p->mutex);
#endif
        ante =
            convolve_with_gaussian (images[e - 1], nligs[e - 1], ncols[e - 1],
                                      ncan, 1.6);
#if defined(WITH_OPENMP)		
		pthread_mutex_unlock (p->mutex);
#endif		
        /* under sampling of the image by a bicubic interpolation */
        images[e] =
            redimension (ante, nligs[e - 1], ncols[e - 1], nligs[e], ncols[e],
                             ncan, gsl_interp2d_bicubic);
        free (ante);
    }
    return NULL;
}

/*===========================================================================*/
void
write_image (char *nom, double *d, int nlig, int ncol, int ncan)
{
    float *f;
    int j;
    int npix = ncol * nlig * ncan;

    f = (float *) malloc (npix * sizeof (float));
    for (j = 0; j < npix; j++)
        f[j] = (float) d[j];
    iio_save_image_float_split (nom, f, ncol, nlig, ncan);

    free (f);
    return;
}

/*===========================================================================*/
double *
redimension (double *im, int nlig, int ncol, int neo_nlig, int neo_ncol,
                 int ncan, const gsl_interp2d_type * type)
{

    cerce_t *cerces = NULL;

    cerces = compute_splines (im, nlig, ncol, ncan, type);

    /* computation of the interpolated values */
    int c, l, n;
    double *neo;
    int xlig = neo_nlig;
    int xcol = neo_ncol;
    double x, y, z;

    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;

    neo = (double *) malloc (xlig * xcol * ncan * sizeof (double));

    for (l = 0; l < xlig; l++)
    {
        y = (nlig - 1) * 1.0 * l / (xlig - 1);

        for (c = 0; c < xcol; c++)
        {
            x = (ncol - 1) * 1.0 * c / (xcol - 1);

            for (n = 0; n < ncan; n++)
            {
                /* pick up the coefficients of spline interpolation */
                spline = cerces[n].spline;
                xacc = cerces[n].xacc;
                yacc = cerces[n].yacc;

                if (0.0 <= x && x <= (double) (ncol - 1) && 0.0 <= y
                        && y <= (double) (nlig - 1))
                {
                    z = gsl_spline2d_eval (spline, x, y, xacc, yacc);
                }
                else
                {
                    z = 0.0;
                }
                neo[xcol * xlig * n + l * xcol + c] = z;
            }
        }
    }
    free_splines (cerces, ncan);
    return neo;
}

/*===========================================================================*/
double *
lucas_kanade (pyramide_t * pyramide_g, pyramide_t * pyramide_d,
              int half_vcol, int half_vlig, int nb_procs,
              int nb_iter, float sigma_dv, int nb_scales)
{
    printf ("[INFO] entering Lucas Kanade\n");
	(void) nb_procs;
    double coef;
    int i;

    /* coefficients of spline's interpolation */
    cerce_t *cerces = NULL;
    double *tmp;
    double *u_gradx;

    /* relative displacement between image n and n+1 */
    double *disp_e = NULL;
    double *disp_cumul = NULL;
    double *disp_init;

    int *ncols, *nligs;
    int ncol_e, nlig_e;
    int ncan;
    int nc, nl;

    double **images_g;
    double **images_d;

    int e;

    /* image column, line, channel */
    int c, l; //, n;

    for (i = 0; i < nb_scales; i++)
    {
        nl = pyramide_g->nlig[i];
        nc = pyramide_g->ncol[i];
        printf ("[INFO] image %d = %d lines x %d col.\n", i, nl, nc);
    }

    nlig_e = pyramide_g->nlig[nb_scales - 1];
    ncol_e = pyramide_g->ncol[nb_scales - 1];

    nligs = pyramide_g->nlig;
    ncols = pyramide_g->ncol;
    ncan = pyramide_g->ncan;

    images_g = pyramide_g->images;
    images_d = pyramide_d->images;

    /* initial disparity */
    disp_init = calloc (nlig_e * ncol_e, sizeof (double));
    if (disp_init == NULL)
    {
        perror ("lucas_kanade()\n");
        exit (EXIT_FAILURE);
    }

    for (e = nb_scales - 1; e >= 0; e--)
    {
        printf ("[INFO] processing at scale %d\n", e);
        free (disp_cumul);
        disp_cumul = calloc (nligs[e] * ncols[e], sizeof (double));
        if (disp_cumul == NULL)
        {
            perror ("lucas_kanade()\n");
            exit (EXIT_FAILURE);
        }

        /* images blurring */
        tmp =
            convolve_with_gaussian (images_g[e], nligs[e], ncols[e], ncan,
                                      sigma_dv);
        free (images_g[e]);
        images_g[e] = tmp;
        tmp =
            convolve_with_gaussian (images_d[e], nligs[e], ncols[e], ncan,
                                      sigma_dv);
        free (images_d[e]);
        images_d[e] = tmp;

        /* computation of gradient x of left image */
        u_gradx = compute_gradient_x (images_g[e], nligs[e], ncols[e], ncan);

        /* computation of the coefficients of bicubic spline interpolation */
        cerces = compute_splines (images_d[e], nligs[e], ncols[e], ncan,
                                    gsl_interp2d_bicubic);

        disp_e = calloc (nligs[e] * ncols[e], sizeof (double));
        if (disp_e == NULL)
        {
            perror ("main()\n");
            exit (EXIT_FAILURE);
        }

#if defined(WITH_OPENMP)
#ifdef _OPENMP		
		#pragma omp parallel for
#endif		
#endif		
		
		for (int indice = 0; indice < nligs[e]*ncols[e]; indice++)
		{
			estimate_shift_thread (indice, cerces, images_g[e], images_d[e], u_gradx,
									 disp_init, disp_e, nligs[e], ncols[e], ncan,
									 half_vcol, half_vlig, nb_iter);
		}
		
        free_splines (cerces, ncan);

        for (l = 0; l < nligs[e]; l++)
        {
            for (c = 0; c < ncols[e]; c++)
            {
                disp_cumul[l * ncols[e] + c] =
                    disp_init[l * ncols[e] + c] + disp_e[l * ncols[e] + c];
            }
        }

        free (disp_e);
        free (u_gradx);
        free (disp_init);

        /* bilinear interpolation */
        if (e > 0)
        {
            disp_init =
                redimension (disp_cumul, nligs[e], ncols[e], nligs[e - 1],
                                 ncols[e - 1], 1, gsl_interp2d_bilinear);
            coef = (double) (ncols[e - 1] - 1) / (ncols[e] - 1);

            for (l = 0; l < nligs[e - 1]; l++)
            {
                for (c = 0; c < ncols[e - 1]; c++)
                {
                    disp_init[l * ncols[e - 1] + c] =
                        coef * disp_init[l * ncols[e - 1] + c];
                }
            }
        }

        printf ("[OK] estimation of optical flow at scale %d\n", e);
    }
    // end of an era, end of the multiscale analysis

    /* do not forget to set the good displacement direction ! */
    e = 0;
    for (l = 0; l < nligs[e]; l++)
    {
        for (c = 0; c < ncols[e]; c++)
        {
            disp_cumul[l * ncols[e] + c] = -disp_cumul[l * ncols[e] + c];
        }
    }

    return (disp_cumul);
}

/*===========================================================================*/
int
algorithm (char *ficim_l, char *ficim_r, char *fic_disp, int half_vcol,
           int half_vlig, int nb_iter, int nb_scales, double sigma_dv,
           long int nb_procs, double seuil, int type_test)
{
    float *fim_l;
    float *fim_r;
    int nlig, ncol, ncan;

    /* reading files */
    fim_l = iio_read_image_float_split (ficim_l, &ncol, &nlig, &ncan);
    fim_r = iio_read_image_float_split (ficim_r, &ncol, &nlig, &ncan);
    if (fim_l == NULL || fim_r == NULL)
    {
        printf ("[ERROR] file unreadable %p %p\n", fim_l, fim_r);
        exit(EXIT_FAILURE);
    }
    printf ("[OK] loading image of size %d x %d x %d\n", ncol, nlig, ncan);

    /* two pyramids building */
    pyramide_t *pyramide_l;
    pyramide_t *pyramide_r;
    build_pyramid (&pyramide_l, &pyramide_r, fim_l, fim_r,
                      nlig, ncol, ncan, nb_scales);

    double *disp_left, *disp_right;
    /* matching left/right */
    printf ("[INFO] matching left/right\n");
    disp_left = lucas_kanade (pyramide_l, pyramide_r,
                              half_vcol, half_vlig, nb_procs, nb_iter, sigma_dv,
                              nb_scales);

    /* matching right/left */
    printf ("[INFO] matching right/left\n");
    disp_right = lucas_kanade (pyramide_r, pyramide_l,
                               half_vcol, half_vlig, nb_procs, nb_iter,
                               sigma_dv, nb_scales);

    free_pyramid (pyramide_l, nb_scales);
    free_pyramid (pyramide_r, nb_scales);

    /* computation of occulted pixels of the left disparity map if needed */
    if (seuil > 0.0)
    {
        printf ("[INFO] setting NaN disparities\n");
        setting_nan_disparities (disp_left, disp_right, seuil, type_test, ncol,
                                 nlig);
    }

    free (fim_l);
    free (fim_r);

    write_image (fic_disp, disp_left, nlig, ncol, 1);
    free (disp_left);
    free (disp_right);

    printf ("[INFO] end of processing !\n");
    return 0;
}
