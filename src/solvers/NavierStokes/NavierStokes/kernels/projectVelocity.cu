/***************************************************************************//**
 * \file projectVelocity.cu
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to update the velocity field
 */


#include "projectVelocity.h"

namespace kernels 
{
/*
 * gets velocity from intermediate velocity and pressure
 * param u velocities
 * param uhat intermediate velocities
 * param uold velocites from previous time step
 * param presure pressure
 * param dx distance between nodes in the x direction (measured between node sides, where u velocites are stored)
 * param dt change in time
 * param nx number of cells in x direction
 * param ny number of cells in y direction
 */
__global__
void project_velocity_X_nobody(double *u, double *uhat, double *uold, double *pressure, double *dx, double dt, int nx, int ny)
{
	int i	= threadIdx.x + (blockDim.x * blockIdx.x),
		I	= i % (nx-1),
		J 	= i / (nx-1),
		ip  = J*nx + I,
		numU	= (nx-1)*ny;

	if (i >= numU)
		return;

	uold[i] = u[i];

	u[i] = uhat[i] - dt*(pressure[ip+1]-pressure[ip]) / (0.5*dx[I+1]+0.5*dx[I]);
}

/*
 * gets velocity from intermediate velocity and pressure
 * param u velocities
 * param uhat intermediate velocities
 * param uold velocites from previous time step
 * param presure pressure
 * param dx distance between nodes in the x direction (measured between node sides, where u velocites are stored)
 * param dt change in time
 * param nx number of cells in x direction
 * param ny number of cells in y direction
 */
__global__
void project_velocity_Y_nobody(double *u, double *uhat, double *uold, double *pressure, double *dy, double dt, int nx, int ny)
{
	int numU= (nx-1)*ny,
		i	= threadIdx.x + (blockDim.x * blockIdx.x),
		I	= i % nx,
		J	= i / nx,
		ip	= J*nx + I,
		numUV	= (ny-1)*nx + numU;

	i += numU;

	if (i >= numUV)
		return;

	uold[i] = u[i];

	u[i] = uhat[i] - dt*(pressure[ip+nx]-pressure[ip]) / (0.5*dy[J+1]+0.5*dy[J]);
}
}//end namespace kernels
