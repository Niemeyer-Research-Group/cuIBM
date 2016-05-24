/***************************************************************************//**
 * \file L.h
 * \author Chris Minar
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

//note: for all "Lhat..." uhat and Lhat should be input instead of u and l
namespace kernels
{
__global__
void Lmidx(double *L, double *u, double *dx, double *dy, int nx, int ny, double nu);

__global__
void Lbcx(double *L, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny, double nu);

__global__
void Lmidy(double *L, double *u, double *dx, double *dy, int nx, int ny, double nu);

__global__
void Lbcy(double *L, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny, double nu);

} // end of namespace kernels
