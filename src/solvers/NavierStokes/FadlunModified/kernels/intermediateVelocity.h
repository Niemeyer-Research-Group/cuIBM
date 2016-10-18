/***************************************************************************//**
 * \file
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
void updateBoundaryX(double *u, double *xp, double *dx, double dt, double Uinf, int nx, int ny);

__global__
void updateBoundaryY(double *u, double *xp, double *dx, double dt, double Vinf, int nx, int ny);

__global__
void updateRHS1forIBX(int *hybridTagsUV, int *ghostTagsUV, double *rhs, double *a, double *b, double *uv, int nx, int ny);

__global__
void updateRHS1forIBY(int *hybridTagsUV, int *ghostTagsUV, double *rhs, double *a, double *b, double *uv, int nx, int ny);
} // end of namespace kernels
