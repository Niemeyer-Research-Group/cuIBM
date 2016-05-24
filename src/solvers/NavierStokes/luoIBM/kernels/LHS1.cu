/***************************************************************************//**
 * \file LHS1.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the intermediate velocity solve
 */

#include "LHS1.h"

namespace kernels
{
__global__
void LHS1_mid_luo_X(int *row, int *col, double *val, int *ghostTagsUV, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;

	//int numE = i*5;
	//			top row - corner    mid           sides    current row
	int numE = (nx-1)*4 - 2      + (J-1)*(5*(nx-1)  - 2) + I*5 - 1;

	double temp = 1;
	//EAST
	row[numE] = i;
	col[numE] = i+1;
	val[numE] = -0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
	temp += 0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
	numE++;

	//WEST
	row[numE] = i;
	col[numE] = i-1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5));
	temp += 0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5));
	numE++;

	//NORTH
	row[numE] = i;
	col[numE] = i+(nx-1);
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5));
	temp += 0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5));
	numE++;

	//SOUTH
	row[numE] = i;
	col[numE] = i-(nx-1);
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5));
	temp += 0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5));
	numE++;

	//CENTER
	row[numE] = i;
	col[numE] = i;
	val[numE] = temp;
	numE++;
}

__global__
void LHS1_mid_luo_Y(int *row, int *col, double *val, int *ghostTagsUV, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		i = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;

	int numE = (nx-1)*ny*5 - 2*ny-2*(nx-1)  +  nx*4-2  + (J-1)*(nx*5 - 2) + I*5 - 1;
	double temp = 1;

	//EAST
	row[numE] = i;
	col[numE] = i+1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5));
	temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5));
	numE++;

	//WEST
	row[numE] = i;
	col[numE] = i-1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5));
	temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5));
	numE++;

	//NORTH
	row[numE] = i;
	col[numE] = i + nx;
	val[numE] = -0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
	temp += 0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
	numE++;

	//SOUTH
	row[numE] = i;
	col[numE] = i-nx;
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
	temp += 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
	numE++;

	//CENTER
	row[numE] = i;
	col[numE] = i;
	val[numE] = temp;
	numE++;
}

}//end kernel
