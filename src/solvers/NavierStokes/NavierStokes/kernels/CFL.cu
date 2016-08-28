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
	if (I==nx-1||J==ny-1)
		return;
	cfl[ip] = dt*(abs(u[iu])/dx[I] + abs(u[iv])/dy[J]);
}

__global__
void testDistance(double *distance,int *ghostTagsUV, int *ghostTagsP, double *xu, double *xv, double *yu, double *yv, double midX, double midY,
					int *i_start, int *j_start, int width, int nx, int ny)
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iu = J*(nx-1) + I,
		iv	= (nx-1)*ny  +  nx*J +I,
		ip = nx*J +I;
	if (idx >= (ny-1)*(nx-1)) //return if we're out of bound
		return;
	if (ghostTagsP[ip]==-1) //return if we're outside body
		return;
	distance[ip] = sqrt(pow(xv[I]-midX,2) + pow(yu[J]-midY,2));
}
}
