/***************************************************************************//**
 * \file LHS2.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the poission solve
 */

#include "LHS2.h"

namespace kernels
{
__global__
void LHS2_mid_luo(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt, int *count, double *stencilCoef, double *interpCoef,
					double *detA, int *hybridTagsP, int *ghostTagsP, double *alpha,
					double *xv, double *yu,
					double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
					double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
					/*double *q1, double *q2, double *q3, double *q4,*/ //not used
					bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag, //not currently used
					/*double *x1, double *x2, double *x3, double *x4, //not used
					double *y1, double *y2, double *y3, double *y4,*/ //not used
					int *index1, int *index2, int *index3, int *index4)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny)
		return;
	int	I	= ip % nx,
		J	= ip / nx;

	if (I == 0 || I == nx-1 || J == 0 || J == ny-1)
		return;

	int numE = nx*4-2 + (J-1)*(nx*5-2) + I*5-1;
	double	temp = 0;

	if (hybridTagsP[ip]>0)//if were at hybrid node
	{
	int index[4] = {index1[ip], index2[ip], index3[ip], index4[ip]};
	int cardinal[5] = {ip+nx, ip+1, ip-nx, ip-1, ip};//n e s w p
	double tempInterp[4];
	double tempStencil[5];
	bool interpMatch[4] = {false, false, false, false};
	//bool qflag[4] = {q1flag[ip], q2flag[ip], q3flag[ip], q4flag[ip]};
	double	x = xv[I],
			y = yu[J];

	//calculate the pressure coefficients for the stencil pressure calculation
	tempStencil[0] = dt/(dy[J]*(dy[J]+dy[J+1])*0.5); //N
	tempStencil[1] = dt/(dx[I]*(dx[I]+dx[I+1])*0.5); //E
	tempStencil[2] = dt/(dy[J]*(dy[J]+dy[J-1])*0.5); //s
	tempStencil[3] = dt/(dx[I]*(dx[I]+dx[I-1])*0.5); //w
	tempStencil[4] = tempStencil[0] + tempStencil[1] + tempStencil[2] + tempStencil[3]; //P
	stencilCoef[ip] = tempStencil[4];

	//divide the eq by the pressure coefficient, add weighting
	for (int i=0;i<4;i++)
		tempStencil[i] = (1-alpha[ip])*tempStencil[i]/tempStencil[4];

	//calculate pressure coefficients for the interpolation pressure calculation	//is multiplied by pressure from corner:
	//					a0        a1        a2        a3
	tempInterp[0] = ( b11[ip] + b21[ip]*x + b31[ip]*y + b41[ip]*x*y )/detA[ip];			//1
	tempInterp[1] = ( b12[ip] + b22[ip]*x + b32[ip]*y + b42[ip]*x*y )/detA[ip];			//2
	tempInterp[2] = ( b13[ip] + b23[ip]*x + b33[ip]*y + b43[ip]*x*y )/detA[ip];			//3
	tempInterp[3] = ( b14[ip] + b24[ip]*x + b34[ip]*y + b44[ip]*x*y )/detA[ip];			//4

	//figure out which interpolation term is being multiplied by pressure at ip
	for (int i=0;i<4;i++)
	{
		if (index[i] == ip)
		{
			interpMatch[i] = true;
		}
	}

	//if were at the center node, set the center node
	//else multiply by alpha and divide by center coefficient
	for (int i=0;i<4;i++)
	{
		tempInterp[i] = alpha[ip]*tempInterp[i];// / (1 - interpCoef)
	}

	//combine the stencil coefficients and interpolation coefficients
	//loop through each interpolation corner
	for (int i=0;i<4;i++)
	{
		//loop through each side of the poisson stencil the center (nesw)
		for (int j=0;j<4;j++)
		{
			//if the index of the corner matches the index of the a poisson stencil
			//and if the corner at the index is not being used as a BC
			if (index[i] == cardinal[j])
			{
				tempStencil[j] += tempInterp[i]; //note, this should also check if the corner is multiplied by P. if it is using a neuman condition it shouldn't be added (doesnt matter for HN)
				interpMatch[i] = true;
			}
		}
	}
	//add the 4 sides of the Poisson stencil to the sparse matrix
	for (int i=0;i<4;i++)
	{
		row[numE]=ip;
		col[numE]=cardinal[i];
		val[numE]=-tempStencil[i];
		numE++;
	}
	//add the center of the poisson stencil to the sparse matrix
	row[numE]=ip;
	col[numE]=ip;
	val[numE]=1;
	//look at all of the nodes used for interpolation, if they are not coincident with any part of the poisson stencil, add them the the sparse matrix
	//change numE so that the values are added to the end of the sparse matrix
	numE = nx*ny*5 - nx*2 - ny*2 + count[ip]-1;
	for (int j=0;j<4;j++)
	{
	//check that the value isn't at a stencil node
	if (!interpMatch[j])
	{
		row[numE]=ip;
		col[numE]=index[j];
		val[numE]=-tempInterp[j];
		numE++;
	}
	}
	}
	else //if were not at a hybrid node
	{
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
	}

	/*
	//do some jank so the solver works, although this modifies the matricies it doesn't really change the results //flag
	if(row[numE]==col[numE] && col[numE]==(ny/2)*nx+nx/2)
	{
		//val[numE] += val[numE];
	}*/
}
}
