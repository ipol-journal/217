// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2016, Tristan Dagobert <tristan.dagobert@cmla.ens-cachan.fr>
// All rights reserved.

#ifndef GENERIC_TENSOR_H
#define GENERIC_TENSOR_H

#include <omp.h>

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void psi_divergence(
    float *psi1,      //coefficients of divergence
    float *psi2,
    float *psi3,
    float *psi4,
    const float *psi, //robust functional
    const int nx,     //image width
    const int ny      //image height
);

/**
 *
 * Compute the divergence of the optical flow
 *
 */
void divergence(
    const float *u,    //component of optical flow
    const float *psi1, //coefficients of divergence
    const float *psi2,
    const float *psi3,
    const float *psi4,
    const int nx,      //image width
    const int ny,       //image height
    float *div      //computed divergence for u
);


#endif
