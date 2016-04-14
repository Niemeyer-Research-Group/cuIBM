/***************************************************************************//**
 * \file initialise.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */

#include "initialise.h"

namespace kernels
{
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

__global__
void initialiseV(double *u, double *xv, double *yv, double vInitial, double vPerturb, double pi, double xmax, double xmin,double ymax,double ymin, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= nx*(ny-1))
		return;
	int 	idx 	= threadIdx.x + (blockDim.x * blockIdx.x),
			i = idx%nx,
			j = idx/nx;
	idx +=  (nx-1)*ny;

	u[idx] = vInitial + vPerturb * cos(0.5*pi*(2*yv[i]-ymax-ymin)/(ymax-ymin)) * sin( pi * (2*xv[j]-xmax-xmin)/(xmax-xmin));
}
}
