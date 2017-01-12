/***************************************************************************//**
 * \file intermediateVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side for the initial velocity solve
 */


#include "intermediateVelocity.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void size_LHS1(int *hybridTagsUV, int *count, int *startI, int *startJ, int width, int height, int nx, int ny)
{
	int iu,iv,counter = 0;
	for (int j=startJ[0]; j<startJ[0]+height; j++)
	{
		for (int i=startI[0]; i<startI[0]+width; i++)
		{
			iu = j*(nx-1)+i;
			if (hybridTagsUV[iu]>0)
			{
				counter+=1;
				count[iu] = counter;
			}
		}
	}
	for (int j=startJ[0]; j<startJ[0]+height; j++)
	{
		for (int i=startI[0]; i<startI[0]+width; i++)
		{
			iv = j*nx+i + (nx-1)*ny;
			if (hybridTagsUV[iv]>0)
			{
				counter+=1;
				count[iv] = counter;
			}
		}
	}
}

__global__
void update_rhs1_x(double *rhs1, double *ns_rhs, double *interp_rhs, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int iu 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= iu % (nx-1),
		J	= iu / (nx-1);
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;
	rhs1[iu] = rhs1[iu] * ns_rhs[iu] + interp_rhs[iu];
}

__global__
void update_rhs1_y(double *rhs1, double *ns_rhs, double *interp_rhs, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iv = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;

	rhs1[iv] = rhs1[iv] * ns_rhs[iv] + interp_rhs[iv];
}

}
