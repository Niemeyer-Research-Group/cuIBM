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
void generateRHS(double *rhs, double *L, double *Nold, double *N, double *u, double *bc1, double dt, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (ny-1)*nx + (nx-1)*ny)
			return;
	int i 	= threadIdx.x + (blockDim.x * blockIdx.x);
	rhs[i]  = u[i] + dt*(0.5*Nold[i] - 1.5*N[i] + 0.5*L[i]) + bc1[i];
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

__global__
void bc1X(double *u, double *bc1, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, double nu, double dt, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
			return;

	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);

	double temp = 0;

	//East
	if (I == nx-2)
	{
		temp += xp[J] * 0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
	}

	//West
	if (I == 0)
	{
		temp += xm[J] * 0.5*dt*nu * (1/(dx[I]*(dx[I+1]+dx[I])*0.5));
	}

	//North
	if (J == ny-1)
	{
		temp += (2*yp[I] - u[i]) * nu*dt*0.5 / dy[J] / dy[J];
	}

	//South
	if (J == 0)
	{
		temp += (2*ym[I] - u[i]) * nu*dt*0.5 / dy[J] / dy[J];
	}

	bc1[i] = temp;
}

__global__
void bc1Y(double *u, double *bc1, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, double nu, double dt, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
			return;

	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iv = ip + (nx-1)*ny;

	double temp = 0;

	//East
	if (I == nx-1)
	{
		temp += (2*xp[ny+J] - u[iv]) * nu*dt*0.5 / dx[I] / dx[I];
	}

	//West
	if (I == 0)
	{
		temp += (2*xm[ny + J] - u[iv]) * nu*dt*0.5 / dx[I] / dx[I];
	}

	//North
	if (J == ny-2)
	{
		temp += yp[(nx-1) + I] * 0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
	}

	//South
	if (J == 0)
	{
		temp += ym[(nx-1) + I] * 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
	}
	bc1[iv] = temp;
}
} // end of namespace kernels
