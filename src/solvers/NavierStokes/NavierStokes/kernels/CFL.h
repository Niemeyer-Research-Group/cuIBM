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
__global__
void testDistance(double *distance,int *ghostTagsUV, int *ghostTagsP, double *xu, double *xv, double *yu, double *yv, double midX, double midY,
					int *i_start, int *j_start, int width, int nx, int ny);
}
