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
void zeroInside(int *ghostTags, double *value, int points)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= points)
			return;
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x);

	//		  if not inside
	value[i] = (ghostTags[i] != 0) * value[i];
}

__global__//note dx and dy must be equal and uniform at the point the boundary atm for the second line (forcing for the inside) to work
void updateRHS1_luo_Y(int *hybridTagsUV, int *ghostTagsUV, double *rhs, double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *uv, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= nx*(ny-1))
		return;
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x) +  (nx-1)*ny;

	//		  if not outtag  & if not in tag    rhs				if out tag		outside interpolation
	rhs[i]	= (hybridTagsUV[i]==-1) * (ghostTagsUV[i]<=0) * (rhs[i])   +   (hybridTagsUV[i]!=-1) * distance_between_nodes_at_IB[i]/(distance_from_intersection_to_node[i]+distance_between_nodes_at_IB[i]) * uv[i];

}

__global__
void updateRHS1_luo_X(int *hybridTagsUV, int *ghostTagsUV, double *rhs, double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *uv, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (nx-1)*ny)
			return;
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x);

	//		  if not outtag  & if not in tag    rhs				if out tag		outside interpolation //flag inside interpolation?
	//rhs[i]	= (hybridTagsUV[i]==-1) * (ghostTagsUV[i]<=0) * (rhs[i])   +   (hybridTagsUV[i]!=-1) * distance_between_nodes_at_IB[i]/(distance_from_intersection_to_node[i]+distance_between_nodes_at_IB[i]) * uv[i];
	rhs[i] = (ghostTagsUV[i] == -1) * rhs[i];
}
} // end of namespace kernels
