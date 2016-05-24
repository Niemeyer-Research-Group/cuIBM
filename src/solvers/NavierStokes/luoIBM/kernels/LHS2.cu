/***************************************************************************//**
 * \file LHS2.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the poission solve
 */

#include "LHS2.h"

namespace kernels
{
__global__
void LHS2_mid_luo(int *row, int *col, double *val, double *distance_from_u_to_body, double *distance_from_v_to_body, int  *ghostTagsP, int *hybridTagsP, double *dx, double *dy, int nx, int ny, double dt)
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
	//just outside immersed body
	if (hybridTagsP[ip] != -1)
	{
		//EAST
		//check if east is outside body
		if (ghostTagsP[ip+1] == -1)
		{
			row[numE] = ip;
			col[numE] = ip + 1;
			if (distance_from_u_to_body[ip] > dx[I]/2 && distance_from_u_to_body[ip] < dx[I])
			{
				val[numE] = -dt/dx[I]/distance_from_u_to_body[ip];
				temp 	  += dt/dx[I]/distance_from_u_to_body[ip];
			}
			else
			{
				val[numE] = -dt/dx[I]/dx[I];
				temp 	  += dt/dx[I]/dx[I];
			}
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip + 1;
			val[numE] = 0;
			numE++;
		}

		//WEST
		//check if west pressure node is outside the body
		if(ghostTagsP[ip-1] == -1)
		{
			row[numE] = ip;
			col[numE] = ip - 1;
			if (distance_from_u_to_body[ip] > dx[I]/2 && distance_from_u_to_body[ip] < dx[I])
			{
				val[numE] = -dt/dx[I]/distance_from_u_to_body[ip];
				temp 	  += dt/dx[I]/distance_from_u_to_body[ip];
			}
			else
			{
				val[numE] = -dt/dx[I]/dx[I];
				temp 	  += dt/dx[I]/dx[I];
			}

			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip - 1;
			val[numE] = 0;
			numE++;
		}

		//NORTH
		//check if north pressure node is outside body
		if (ghostTagsP[ip+nx] == -1)
		{
			row[numE] = ip;
			col[numE] = ip + nx;
			if (distance_from_v_to_body[ip] > dy[J]/2 && distance_from_v_to_body[ip] < dy[J])//one of dys in the something/dy^2 term isn't dy anymore, it is smaller
			{
				val[numE] = -dt/dy[J]/distance_from_v_to_body[ip];
				temp 	  += dt/dy[J]/distance_from_v_to_body[ip];
			}
			else
			{
				val[numE] = -dt/dy[J]/dy[J];
				temp 	  += dt/dy[J]/dy[J];
			}
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip + nx;
			val[numE] = 0;
			numE++;
		}

		//SOUTH
		//check if south pressure node is outside body
		if (ghostTagsP[ip-nx] == -1)
		{
			row[numE] = ip;
			col[numE] = ip - nx;
			if (distance_from_v_to_body[ip] > dy[J]/2 && distance_from_v_to_body[ip] < dy[J])
			{
				val[numE] = -dt/dy[J]/distance_from_v_to_body[ip];
				temp 	  += dt/dy[J]/distance_from_v_to_body[ip];
			}
			else
			{
				val[numE] = -dt/dy[J]/dy[J];
				temp 	  += dt/dy[J]/dy[J];
			}
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip - nx;
			val[numE] = 0;
			numE++;
		}

	}
	//end just outside immersed body
	//if just inside body
	else if (ghostTagsP[ip] > 0)
	{
		//EAST
		if (ghostTagsP[ip+1] == 0)
		{
			row[numE] = ip;
			col[numE] = ip + 1;
			val[numE] = -dt/dx[I]/dx[I];
			numE++;
			temp 	  += dt/dx[I]/dx[I];
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip + 1;
			val[numE] = 0;
			numE++;
		}

		//WEST
		if (ghostTagsP[ip-1] == 0)
		{
			row[numE] = ip;
			col[numE]= ip - 1;
			val[numE] = -dt/dx[I]/dx[I];
			temp 	  += dt/dx[I]/dx[I];
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip - 1;
			val[numE] = 0;
			numE++;
		}

		//NORTH
		if (ghostTagsP[ip+nx] == 0)
		{
			row[numE] = ip;
			col[numE] = ip + nx;
			val[numE] = -dt/dy[J]/dy[J];
			temp 	  += dt/dy[J]/dy[J];
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip + nx;
			val[numE] = 0;
			numE++;
		}

		//SOUTH
		if (ghostTagsP[ip-nx] == 0)
		{
			row[numE] = ip;
			col[numE] = ip - nx;
			val[numE] = -dt/dy[J]/dy[J];
			temp 	  += dt/dy[J]/dy[J];
			numE++;
		}
		else
		{
			row[numE] = ip;
			col[numE] = ip - nx;
			val[numE] = 0;
			numE++;
		}
	}
	//end just inside body
	//everywhere else
	else
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
	}//end everywhere else
	//MID
	row[numE] = ip;
	col[numE] = ip;
	val[numE] = temp;

	//do some jank so the solver works, although this modifies the matricies it doesn't really change the results //flag
	if(row[numE]==col[numE] && col[numE]==(ny/2)*nx+nx/2)
	{
		val[numE] += val[numE];
	}
}
}
