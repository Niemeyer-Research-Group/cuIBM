/***************************************************************************//**
 * \file initialise.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */

#include "initialise.h"

namespace kernels
{
/*
 * sets all the initial u values
 * param u u velocities
 * param xu x locations of where u is stored
 * param yu y locations of where u is stored
 * param uInitial initial u velocity
 * param uPerturb perturbation coefficient
 * param pi... //flag
 * param xmax highest x value in the domain
 * param xmin lowest x value in the domain
 * param ymax highest y value in the domain
 * param ymin lowest yvalue in the domain
 * nx number of nodes in the x direction
 * ny number of nodes in the y direction
 */
__global__
void initialiseU(double *u, double *xu, double *yu, double uInitial, double uPerturb, double pi, double xmax, double xmin,double ymax,double ymin, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (nx-1)*ny)
			return;
	int 	idx 	= threadIdx.x + (blockDim.x * blockIdx.x),
			i = idx%(nx-1),
			j = idx/(nx-1);

	u[idx] = uInitial + uPerturb * cos(0.5*pi*(2*xu[i]-xmax-xmin)/(xmax-xmin)) * sin( pi * (2*yu[j]-ymax-ymin)/(ymax-ymin));
}

/*
 * sets all the initial v values
 * param u v velocities
 * param xv x locations of where v is stored
 * param yv y locations of where v is stored
 * param vInitial initial v velocity
 * param vPerturb perturbation coefficient
 * param pi... //flag
 * param xmax highest x value in the domain
 * param xmin lowest x value in the domain
 * param ymax highest y value in the domain
 * param ymin lowest yvalue in the domain
 * nx number of nodes in the x direction
 * ny number of nodes in the y direction
 */
__global__
void initialiseV(double *u, double *xv, double *yv, double vInitial, double vPerturb, double pi, double xmax, double xmin,double ymax,double ymin, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= nx*(ny-1))
		return;
	int 	idx 	= threadIdx.x + (blockDim.x * blockIdx.x),
			i = idx%nx,
			j = idx/nx;
	idx +=  (nx-1)*ny;

	u[idx] = vInitial + vPerturb * cos(0.5*pi*(2*yv[j]-ymax-ymin)/(ymax-ymin)) * sin( pi * (2*xv[i]-xmax-xmin)/(xmax-xmin));
}
}
