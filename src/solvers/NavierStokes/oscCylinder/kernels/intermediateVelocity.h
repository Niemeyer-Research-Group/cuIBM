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
void setInsideVelocity(int *ghostTags, double *u, double *uB, double *vB, int nx, int ny);
}
