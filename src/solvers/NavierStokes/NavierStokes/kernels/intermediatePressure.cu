/***************************************************************************//**
 * \file intermediatePressure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side of the poission equation
 */

#include "intermediatePressure.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void intermediatePressureNoBody(double *rhs2, double *uhat, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*ny)
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iu	= (nx-1)*J + I,
		iv	= (nx-1)*ny  +  nx*J +I;

	double temp = 0;

	//EAST
	//if not on the east wall and east is outside the body, add east term
	if (I != nx-1)//not at east boundry
		temp -= uhat[iu]/dx[I];
	else if (I == nx-1)//at east boundry
		temp -= xp[J]/dx[I];

	//WEST
	//if not on west wall and west is outside the body, add west term
	if (I != 0)//not at west boundary
		temp += uhat[iu - 1]/dx[I];
	else if (I == 0)//at the west boundary
		temp += xm[J]/dx[I];

	//NORTH
	//if not on north wall and north is outside the body, add north term
	if (J != ny-1)//not at north boundry
		temp -= uhat[iv]/dy[J];
	else if (J == ny-1)//at north boundry
		temp -= yp[(nx-1)+I]/dy[J];

	//SOUTH
	//if not on south wall and south is outside the body, add south term
	if (J != 0)//not at south boundry
		temp += uhat[iv-nx]/dy[J];
	else if (J == 0)//at south boundry
		temp += ym[(nx-1)+I]/dy[J];

	rhs2[ip] = temp;
}
}
