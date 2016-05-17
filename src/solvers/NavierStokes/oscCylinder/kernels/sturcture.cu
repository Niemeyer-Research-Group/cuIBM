/***************************************************************************//**
 * \file structure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */


#include "structure.h"

namespace kernels
{

__global__
void update_body_viv(double *y, double *vB, double dy, double vnew, int totalPoints)
{
	int i	= threadIdx.x + (blockDim.x * blockIdx.x);
	if (i > totalPoints)
		return;
	vB[i] = vnew;
	y[i] = y[i] + dy;
}

}
