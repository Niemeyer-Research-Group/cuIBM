/***************************************************************************//**
 * \file generateRHS.inl
 * \author Chris Minar
 * \brief Declaration of kernels to generate the right hand side of step 1: solve for uhat
 * \		-G*p -1.5N(u) + 0.5 N(uold) + 0.5 L(u)
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void generateRHS(double *rhs, double *L, double *Nold, double *N, double *u, double *bc1, double dt, int nx, int ny);

__global__
void bc1X(double *u, double *bc1, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, double nu, double dt, int nx, int ny);

__global__
void bc1Y(double *u, double *bc1, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, double nu, double dt, int nx, int ny);
} // end of namespace kernels
