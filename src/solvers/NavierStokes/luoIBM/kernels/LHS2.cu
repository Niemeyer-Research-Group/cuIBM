/***************************************************************************//**
 * \file LHS2.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the poission solve
 */

#include "LHS2.h"

namespace kernels
{
__global__
void LHS2_mid_luo(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny)
		return;
	int	I	= ip % nx,
		J	= ip / nx;

	if (I == 0 || I == nx-1 || J == 0 || J == ny-1)
		return;

	int numE = nx*4-2 + (J-1)*(nx*5-2) + I*5-1;
	double temp = 0;
	//EAST
	row[numE] = ip;
	col[numE] = ip + 1;
	val[numE] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
	numE++;
	temp 	  += dt/(dx[I]*(dx[I]+dx[I+1])*0.5);

	//WEST
	row[numE] = ip;
	col[numE] = ip - 1;
	val[numE] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
	temp 	  += dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
	numE++;

	//NORTH
	row[numE] = ip;
	col[numE] = ip + nx;
	val[numE] = -dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
	temp += dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
	numE++;

	//SOUTH
	row[numE] = ip;
	col[numE] = ip - nx;
	val[numE] = -dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
	temp 	  += dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
	numE++;
	//MID
	row[numE] = ip;
	col[numE] = ip;
	val[numE] = temp;

	//do some jank so the solver works, although this modifies the matricies it doesn't really change the results //flag
	if(row[numE]==col[numE] && col[numE]==(ny/2)*nx+nx/2)
	{
		//val[numE] += val[numE];
	}
}
}
