/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side for the initial velocity solve
 */


#include "CFL.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
//size p
__global__
void calculateCFL(double *cfl, double *u, double *dx, double *dy,
		int nx, int ny, double dt)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*ny)
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iu	= (nx-1)*J + I,
		iv	= (nx-1)*ny  +  nx*J +I;

	cfl[ip] = dt*(abs(u[iu])/dx[I] + abs(u[iv])/dy[J]);
}
}
