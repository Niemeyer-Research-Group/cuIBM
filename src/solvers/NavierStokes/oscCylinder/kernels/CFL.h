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
void calculateCFL(double *cfl, double *u, double *dx, double *dy,
		int nx, int ny, double dt);
}
