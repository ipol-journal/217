
/*
* main.c
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


#include "lucas_kanade.h"

#include <string.h>
#include <stdio.h>
/*===========================================================================*/
#define VERSION "1.05"
/*===========================================================================*/
void display_error (char **argv);
void display_help (char **argv);
void display_usage (char **argv);
void display_version (void);
/*===========================================================================*/
void
display_error (char **argv)
{
    printf ("%s : Incorrect parameters…\n", argv[0]);
    display_usage (argv);
    printf ("To know more, type : %s --help\n", argv[0]);
    return;
}

/*===========================================================================*/
void
display_usage (char **argv)
{
    printf ("Usage: %s --version | --help | ", argv[0]);
    printf
    ("FIC_IM_L FIC_IM_R FIC_DISP NB_PROCS VX VY NB_SCALES SIGMA_DV NB_ITER"
     " THRESHOLD TYPE_TEST\n");
    return;
}

/*===========================================================================*/
void
display_version (void)
{

    printf ("IPOL lucas_kanade.exe %s\n", VERSION);
    printf
    ("Copyright (c) 2016 Tristan Dagobert, CMLA, École Normale Supérieure "
     "Paris-Saclay\n\n");
    printf ("This program is distributed in the hope that it will be useful,\n"
            "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
            "GNU General Public License for more details.\n");
    return;
}

/*===========================================================================*/
void
display_help (char **argv)
{
    display_usage (argv);
    printf ("\n");
    printf
    ("Create the disparity map FIC_DISP from left image FIC_IM_L and right "
     "image FIC_IM_R\nconsidered respectively as the reference and the shifted"
     " images. A neighborhood\nof VX by VY pixels is used around the central "
     "pixel to estimate its displacement,\nwhile a presmoothing by a "
     "Gaussian kernel with a SIGMA_DV standard deviation,\nis applied to "
     "both images at each scale. The inner gradient descent stops after\n"
     "NB_ITER iterations maximum. The left right test consistency is "
     "absolute\nTYPE_TEST=1 or relative TYPE_TEST=0 with a threshold "
     "THRESHOLD. In the absolute\ncase, the threshold value corresponds to a "
     "distance in pixel. In the relative case,\nits value corresponds to a "
     "ratio. If THRESHOLD is negative the consistency test is not done : "
     "the disparity map is dense.\n\n");
    printf ("Note : the input and output images are all types supported by the "
            "iio library. \n");

    printf ("Options :\n");
    printf ("--version            display the program version\n");
    printf ("--help               display this help\n\n");
    printf ("\n");
    printf ("Report bugs to: tristan.dagobert@cmla.ens-cachan.fr\n");
    printf ("Home page: <http://www.ipol.im/>\n");
    return;
}

/*===========================================================================*/
int
main (int argc, char **argv)
{
    char *ficim_l, *ficim_r, *fic_disp;
    int half_vcol, half_vlig;
    int nb_iter, nb_scales;
    double sigma_dv;
    long int nb_procs;
    double seuil;
    int type_test;
    int i;
    for (i = 0; i < argc; i++)
        printf ("%s\n", argv[i]);
    /* processing depending on the parameters number */
    switch (argc)
    {
    case 1:
    {
        display_error (argv);
        exit (EXIT_SUCCESS);
        break;
    }
    case 2:
    {
        if (strcmp ("--help", argv[1]) == 0)
        {
            display_help (argv);
            exit (EXIT_SUCCESS);
        }
        else if (strcmp ("--version", argv[1]) == 0)
        {
            display_version ();
            exit (EXIT_SUCCESS);
        }
        else
        {
            display_error (argv);
            exit (EXIT_SUCCESS);
        }
        break;
    }
    case 12:
    {
        break;
    }
    default:
    {
        display_error (argv);
        exit (EXIT_SUCCESS);
    }
    }

    /* pick up arguments */
    i = 1;

    /* file image left */
    ficim_l = argv[i];
    i += 1;

    /* file image right */
    ficim_r = argv[i];
    i += 1;

    /* file estimated disparities */
    fic_disp = argv[i];
    i += 1;

    /* number of processors implied into this fiasco */
    nb_procs = atoi (argv[i]);
    i += 1;

    /* half neigborhood in x */
    half_vcol = atoi (argv[i]) / 2;
    i += 1;

    /* half neigborhood in y */
    half_vlig = atoi (argv[i]) / 2;
    i += 1;

    /* number of scales */
    nb_scales = atoi (argv[i]);
    i += 1;

    /* sigma of the Gaussian kernel of blur and derivation */
    sigma_dv = atof (argv[i]);
    i += 1;

    /* iterations number used during the gradient descent */
    nb_iter = atoi (argv[i]);
    i += 1;

    /* threshold of the left right consistency test */
    seuil = atof (argv[i]);
    i++;

    /* kind of test used in the consistency test */
    type_test = atoi (argv[i]);

    printf ("Program\t\t\t: «%s»\n", argv[0]);
    printf ("Image left\t\t: «%s»\n", ficim_l);
    printf ("Image right\t\t: «%s»\n", ficim_r);
    printf ("Disparity map\t\t: «%s»\n", fic_disp);
    printf ("Number of procs\t\t: «%ld»\n", nb_procs);
    printf ("Half neigborhood x\t: «%d»\n", half_vcol);
    printf ("Half neigborhood y\t: «%d»\n", half_vlig);
    printf ("Number of scales\t: «%d»\n", nb_scales);
    printf ("Sigma deriv.\t\t: «%7.3f»\n", sigma_dv);
    printf ("Number of iter.\t\t: «%d»\n", nb_iter);
    printf ("Consistency threshold\t: «%7.3f»\n", seuil);
    printf ("Kind of test\t: «%s»\n", type_test == 1 ? "ABSOLUTE" : "RELATIVE");
    printf ("\n");

    algorithm (ficim_l, ficim_r, fic_disp, half_vcol, half_vlig, nb_iter,
               nb_scales, sigma_dv, nb_procs, seuil, type_test);


    exit (EXIT_SUCCESS);
}
