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
/*
 * sums the components of the right hand side of the intermediate velocity equation
 * param rhs right hand side fo the velocity equation
 * param L explicit diffusion terms
 * param N explicit advection terms
 * param u u velocities
 * param bc1 boundary condition terms
 * param dt change in time
 * param nx number of cells in x direction
 * param ny number of cells in y direction
 */
__global__
void generateRHS(double *rhs, double *L, double *Nold, double *N, double *u, double *bc1, double dt, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (ny-1)*nx + (nx-1)*ny)
			return;
	int i 	= threadIdx.x + (blockDim.x * blockIdx.x);
	rhs[i]  = u[i] + dt*(0.5*Nold[i] - 1.5*N[i] + 0.5*L[i]) + bc1[i];
}

/*
 * calculates boundary terms for u intermediate velocity
 * param u u velocities
 * param bc1 boundry terms
 * param ym yminus boundary velocities
 * param yp yplus boundary velocities
 * param xm xminus boundary velocities
 * param xp xplus boundary velocities
 * param dx distance between nodes in the x direction (measured between node sides, where u velocites are stored)
 * param dy distance between nodes in the y direction (measured between node top/bot, where v velocites are stored)
 * param nu viscosity
 * param dt change in time
 * param nx number of cells in x direction
 * param ny number of cells in y direction
 */
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

/*
 * calculates boundary terms for u intermediate velocity
 * param u v velocities
 * param bc1 boundry terms
 * param ym yminus boundary velocities
 * param yp yplus boundary velocities
 * param xm xminus boundary velocities
 * param xp xplus boundary velocities
 * param dx distance between nodes in the x direction (measured between node sides, where u velocites are stored)
 * param dy distance between nodes in the y direction (measured between node top/bot, where v velocites are stored)
 * param nu viscosity
 * param dt change in time
 * param nx number of cells in x direction
 * param ny number of cells in y direction
 */
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
