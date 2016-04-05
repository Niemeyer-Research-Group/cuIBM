/***************************************************************************//**
 * \file N.h
 * \author Chris Minar
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void Nmidx(double *N, double *u, double *dx, double *dy, int nx, int ny);

__global__
void Nbcx(double *N, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny);

__global__
void Nmidy(double *N, double *u, double *dx, double *dy, int nx, int ny);

__global__
void Nbcy(double *N, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny);
} // end of namespace kernels
