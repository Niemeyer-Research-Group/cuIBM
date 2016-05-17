/***************************************************************************//**
 * \file LHS1.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the intermediate velocity solve
 */

#include "LHS1.h"

namespace kernels
{
__global__
void LHS_BC_X(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);
	if (I != 0 && I != nx-2 && J != 0 && J != ny-1)
		return;

	double temp = 1;
	int numE = 0;
	if (J == 0)
	{
		numE = I*4;
		if (I != 0)
			numE -= 1;
	}
	else if (J == ny-1)
	{
		numE = (nx-1)*4 - 2 + (J-1)*(5*(nx-1)  - 2) + I*4;
		if (I != 0)
			numE-=1;
	}
	else
	{
		if (I == 0)
			numE = (nx-1)*4 - 2 + (J-1)*(5*(nx-1)  - 2) + I*5;
		else
			numE = (nx-1)*4 - 2 + (J-1)*(5*(nx-1)  - 2) + I*5 - 1;
	}

	//EAST
	if(I != nx-2)//check if on east boundary
	{
		row[numE] = i;
		col[numE] = i+1;
		val[numE] = -0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
		temp += 0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5));
	}

	//WEST
	if(I != 0)//check if on west boundary
	{
		row[numE] = i;
		col[numE] = i-1;
		val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5));
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5));
	}

	//NORTH
	if(J != ny-1)//check if on north boundary
	{
		row[numE] = i;
		col[numE] = i+(nx-1);
		val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5));
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J])*0.5));
	}

	//SOUTH
	if(J != 0)//check if on south boundary
	{
		row[numE] = i;
		col[numE] = i-(nx-1);
		val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5));
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J])*0.5));
	}

	//CENTER
	row[numE] = i;
	col[numE] = i;
	val[numE] = temp;
	numE++;
}

__global__
void LHS_BC_Y(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		i = ip + (nx-1)*ny;
	if (I != 0 && I != nx-1 && J != 0 && J != ny-2)
		return;

	int numE = (nx-1)*ny*5 - 2*ny-2*(nx-1);
	if (J == 0)
	{
		numE += I*4;
		if (I != 0)
			numE -= 1;
	}
	else if (J == ny-2)
	{
		numE += nx*4 - 2 + (J-1)*(5*nx - 2) + I*4;
		if (I != 0)
			numE-=1;
	}
	else
	{
		if (I == 0)
			numE += nx*4 - 2 + (J-1)*(5*nx - 2) + I*5;
		else
			numE += nx*4 - 2 + (J-1)*(5*nx - 2) + I*5 - 1;
	}
	double temp = 1;

	//EAST
	if(I != nx-1)//check if on east boundary
	{
		row[numE] = i;
		col[numE] = i+1;
		val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5));
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I])*0.5));
	}

	//WEST
	if(I != 0)//check if  on west boundary
	{
		row[numE] = i;
		col[numE] = i-1;
		val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5));
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I])*0.5));
	}

	//NORTH
	if(J != ny-2)//check if on north boundary
	{
		row[numE] = i;
		col[numE] = i + nx;
		val[numE] = -0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
		temp += 0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5));
	}

	//SOUTH
	if(J != 0)//check if on south boundary
	{
		row[numE] = i;
		col[numE] = i-nx;
		val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
		numE++;
	}
	else
	{
		temp += 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));
	}

	//CENTER
	row[numE] = i;
	col[numE] = i;
	val[numE] = temp;
	numE++;
}

__global__
void LHS_mid_X_nobody(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;

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
void LHS_mid_Y_nobody(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		i = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;

	//         (              numU       )     (row1)    (rows2-before me)  (current row)
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
}
