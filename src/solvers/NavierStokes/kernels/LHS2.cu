#include "LHS2.h"

namespace kernels
{
__global__
void LHS2_mid(int *row, int *col, double *val, double *distance_from_u_to_body, double *distance_from_v_to_body, int  *tagsP, int *tagsPOut, double *dx, double *dy, int nx, int ny, double dt)
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
	if (tagsPOut[ip] != -1)
	{
		//EAST
		//check if east is outside body
		if (tagsP[ip+1] == -1)
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
		if(tagsP[ip-1] == -1)
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
		if (tagsP[ip+nx] == -1)
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
		if (tagsP[ip-nx] == -1)
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
	else if (tagsP[ip] > 0)
	{
		//EAST
		if (tagsP[ip+1] == 0)
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
		if (tagsP[ip-1] == 0)
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
		if (tagsP[ip+nx] == 0)
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
		if (tagsP[ip-nx] == 0)
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

	//do some jank shit so the solver works, although this modifies the matricies it doesn't really change the results
	if(row[numE]==col[numE] && col[numE]==(ny/2)*nx+nx/2)
	{
		val[numE] += val[numE];
	}
}

__global__
void LHS2_BC(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny)
		return;
	int	I	= ip % nx,
		J	= ip / nx;
	if (I != 0 && I != nx-1 && J != 0 && J != ny-1)
		return;
	int numE = 0;
	if (J == 0)
	{
		numE = I*4;
		if (I!=0)
			numE-=1;
	}
	else if (J == ny-1)
	{
		numE = nx*4-2 + (J-1)*(nx*5-2) + I*4;
		if (I!=0)
			numE-=1;
	}
	else
	{
		numE = nx*4-2   +    (J-1)*(nx*5 - 2) + I*5;
		if (I != 0)
			numE -= 1;
	}

	double temp = 0;

	//EAST
	//if not on the east wall and east is outside the body, add east term
	if (I != nx-1)//not at east boundry
	{
		row[numE] = ip;
		col[numE] = ip + 1;
		val[numE] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
		numE++;
		temp += dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
	}

	//WEST
	//if not on west wall and west is outside the body, add west term
	if (I != 0)//not at west boundary
	{
		row[numE] = ip;
		col[numE] = ip - 1;
		val[numE] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
		temp += dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
		numE++;
	}

	//NORTH
	//if not on north wall and north is outside the body, add north term
	if (J != ny-1)//not at north boundry
	{
		row[numE] = ip;
		col[numE] = ip + nx;
		val[numE] = -dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
		temp += dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
		numE++;
	}

	//SOUTH
	//if not on south wall and south is outside the body, add south term
	if (J != 0)//not at south boundry
	{
		row[numE] = ip;
		col[numE] = ip - nx;
		val[numE] = -dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
		temp += dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
		numE++;
	}

	//MID
	row[numE] = ip;
	col[numE] = ip;
	val[numE] = temp;

	//do some jank shit so the solver works, although this modifies the matricies it doesn't really change the results
	//if (ip == 0)
	//	val[numE] += val[numE];
}

__global__
void LHS2_mid_nobody(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny)
		return;
	int	I	= ip % nx,
		J	= ip / nx;

	if (I == 0 || I == nx-1 || J == 0 || J == ny-1)
		return;

	int numE= nx*4-2   +    (J-1)*(nx*5 - 2) + I*5 - 1;
	double temp = 0;
	//EAST
	row[numE] = ip;
	col[numE] = ip + 1;
	val[numE] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
	numE++;
	temp += dt/(dx[I]*(dx[I]+dx[I+1])*0.5);

	//WEST
	row[numE] = ip;
	col[numE] = ip - 1;
	val[numE] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
	temp += dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
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
	temp += dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
	numE++;

	//MID
	row[numE] = ip;
	col[numE] = ip;
	val[numE] = temp;

	//do some jank shit so the solver works, although this modifies the matricies it doesn't really change the results
	if(row[numE]==col[numE] && col[numE]==(ny/2)*nx+nx/2)
	{
		val[numE] += val[numE];
	}
}























}
