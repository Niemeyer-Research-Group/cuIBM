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
void setInsideVelocity(int *ghostTags, double *u, double *uB, double *vB, int nx, int ny) //flag doesn't need to cover whole domain, could only span over the bounding box
{																  //flag kernel could mess up if the body is too close to the edge because were doing the x values and y values in the same kernel
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x),
			I	= i % (nx-1),
			J	= i / (nx-1),
			iu	= J*(nx-1) + I,
			iv	= J*nx + I + (nx-1)*ny;

	if (iu >= (nx-1)*ny) //flag indexing is janky for doing x and y at the same time
			return;
	//			 not at inside edge             at inside edge
	u[iu] = (ghostTags[iu] != 0) * u[iu] + (ghostTags[iu] == 0) * uB[0]; //flag won't work for rotating bodies because were not getting a local body velocity
	u[iv] = (ghostTags[iv] != 0) * u[iv] + (ghostTags[iv] == 0) * vB[0];
}
}
