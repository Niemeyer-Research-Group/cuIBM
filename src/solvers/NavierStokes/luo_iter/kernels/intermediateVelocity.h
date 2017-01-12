/***************************************************************************//**
 * \file
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
void size_LHS1(int *hybridTagsUV, int *count, int *startI, int *startJ, int width, int height, int nx, int ny);
__global__
void update_rhs1_x(double *rhs1, double *ns_rhs, double *interp_rhs, int nx, int ny);
__global__
void update_rhs1_y(double *rhs1, double *ns_rhs, double *interp_rhs, int nx, int ny);
}
