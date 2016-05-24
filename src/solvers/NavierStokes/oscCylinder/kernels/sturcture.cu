/***************************************************************************//**
 * \file structure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */


#include "structure.h"

namespace kernels
{

/*
 * Updates all the velocities and positions of the body nodes
 * param double y y positions of the nodes
 * param double vB v velocities of body nodes
 * param double dy change in y location
 * param vnew new v velocity
 * param totalPoints total number of body points
 */
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
