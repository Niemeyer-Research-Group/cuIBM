/***************************************************************************//**
 * \file intermediateVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side for the initial velocity solve
 */


#include "pressure.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void size_LHS2(int *hybridTagsP, int *count, int *startI, int *startJ, int width, int height, int nx, int ny)
{
	int ip,counter = 0;
	for (int j=startJ[0]; j<startJ[0]+height; j++)
	{
		for (int i=startI[0]; i<startI[0]+width; i++)
		{
			ip = j*nx+i;
			if (hybridTagsP[ip]>0)
			{
				counter+=1;
				count[ip] = counter;
			}
		}
	}
}

__global__
void update_rhs2(double *rhs2, double *ns_rhs, double *interp_rhs, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-1)
		return;
	rhs2[ip] = rhs2[ip] * ns_rhs[ip] + interp_rhs[ip];
}

}
