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
//updates the east boundary for use with a convective boundary condition
__global__
void updateBoundaryX(double *u, double *xp, double *dx, double dt, double Uinf, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= ny)
		return;
	int J = threadIdx.x + (blockDim.x * blockIdx.x),
		I = nx-1,
		i = J*(nx-1) + I;
	double beta = Uinf * dt / dx[I];
	xp[J] = xp[J]*(1-beta) + beta*u[i-1];
}

__global__
void updateBoundaryY(double *u, double *xp, double *dx, double dt, double Vinf, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= ny - 1)
		return;
	int J = threadIdx.x + (blockDim.x * blockIdx.x),
		I = nx,
		i = J*nx + I,
		numU = (nx-1)*ny;
	double beta = Vinf * dt / dx[I-1];
	xp[J+ny] = xp[J+ny]*(1-beta) + beta*u[i + numU-1];
}

__global__
void updateRHS1forIBX(int *tags, int *tagsIn, double *rhs, double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *uv, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (nx-1)*ny)
			return;
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x);

	//		  if not outtag  & if not in tag    rhs				if out tag		outside interpolation //flag inside interpolation?
	rhs[i]	= (tags[i]==-1) * (tagsIn[i]<=0) * (rhs[i])   +   (tags[i]!=-1) * distance_between_nodes_at_IB[i]/(distance_from_intersection_to_node[i]+distance_between_nodes_at_IB[i]) * uv[i];
}

__global__//note dx and dy must be equal and uniform at the point the boundary atm for the second line (forcing for the inside) to work
void updateRHS1forIBY(int *tags, int *tagsIn, double *rhs, double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *uv, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= nx*(ny-1))
		return;
	int 	i 	= threadIdx.x + (blockDim.x * blockIdx.x) +  (nx-1)*ny;

	//		  if not outtag  & if not in tag    rhs				if out tag		outside interpolation
	rhs[i]	= (tags[i]==-1) * (tagsIn[i]<=0) * (rhs[i])   +   (tags[i]!=-1) * distance_between_nodes_at_IB[i]/(distance_from_intersection_to_node[i]+distance_between_nodes_at_IB[i]) * uv[i];
}
} // end of namespace kernels
