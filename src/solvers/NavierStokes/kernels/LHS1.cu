#include "LHS1.h"

namespace kernels
{
__global__
void LHS_mid_X(int *row, int *col, double *val, int *tags, int *tags2, int *tagsIn, double *a, double *b, double *dx, double *dy, double dt, double nu, int nx, int ny)
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
	if( (tags[i] == -1  && tagsIn[i] == -1) || tagsIn[i] == 0)// if point isn't tagged
	{
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
	//end untagged
	else if (tags[i]!=-1) //if point is tagged
	{
		//ADJACENT POINT
		if (tags[i] == tags2[i] - 1) // right is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i + 1;
		}
		else if (tags[i] == tags2[i] + 1) //left is away from surface
		{
			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i - 1;
		}
		else if (tags[i] == tags2[i] + (nx-1)) // below is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i - (nx-1);
		}
		else if (tags[i] == tags2[i] - (nx-1)) // above is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i + (nx-1);
		}



		row[numE] = i;
		val[numE] = -a[i]/(a[i]+b[i]);
		numE++;

		//CENTER
		row[numE] = i;
		col[numE] = i;
		val[numE] = 1;
		numE++;
	}
	//end tagged
	else if (tagsIn[i] > 0) //inner point is tagged Note:: dx and dy must be uniform in a section with a body...
	{
		//ADJACENT POINT
		if (tags[i+1] != -1)// go right to body
		{
			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i-1;
			val[numE] = -(dx[I]-a[i+1]) / (dx[I] - a[i+1] + b[i+1]);
		}
		else if (tags[i-1] != -1) //go left to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i + 1;
			val[numE] = -(dx[I]-a[i-1]) / (dx[I] - a[i-1] + b[i-1]);
		}
		else if (tags[i+(nx-1)] != -1)//go north to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i - (nx-1);
			val[numE] = -(dx[I]-a[i+(nx-1)]) / (dx[I] - a[i+(nx-1)] + b[i+(nx-1)]);
		}
		else if (tags[i-(nx-1)] != -1)//go south to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE ++;

			row[numE] = i;
			col[numE] = i - (nx-1);
			val[numE] = 0;
			numE ++;

			col[numE] = i + (nx-1);
			val[numE] = -(dx[I]-a[i-(nx-1)]) / (dx[I] - a[i-(nx-1)] + b[i-(nx-1)]);
		}
		row[numE] = i;
		numE++;

		//CENTER
		col[numE] = i;
		row[numE] = i;
		val[numE] = 1;
		numE++;
	}
}

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
void LHS_mid_Y(int *row, int *col, double *val, int *tags, int *tags2, int *tagsIn, double *a, double *b, double *dx, double *dy, double dt, double nu, int nx, int ny)
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

	if((tags[i] == -1 && tagsIn[i] == -1)||tagsIn[i] == 0)	//if not tagged
	{
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
	//end untagged section
	else if (tags[i]>0) //if point is tagged
	{
		if (tags[i] == tags2[i] - 1) //right is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i - nx;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i + 1;
		}
		else if (tags[i] == tags2[i] + 1)//left  is away from surface
		{
			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i - nx;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i - 1;
		}
		else if (tags[i] == tags2[i] + nx)//below is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i - nx;
		}
		else if (tags[i] == tags2[i] - nx)//above is away from surface
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 - nx;
			val[numE] = 0;
			numE++;

			col[numE] = i + nx;
		}
		row[numE] = i;
		val[numE] = -a[i]/(b[i]+a[i]);
		numE++;

		row[numE] = i;
		col[numE] = i;
		val[numE] = 1;
		numE++;
	}//end tagged

	else if (tagsIn[i] != -1) //inner point is tagged Note:: dx and dy must be uniform in a section with a body...
	{
		//setup adjacent point
		if (tags[i+1] != -1)// go right to body
		{
			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i - nx;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i-1;
			val[numE] = -(dx[I]-a[i+1]) / (dx[I] - a[i+1] + b[i+1]);
		}
		else if (tags[i-1] != -1) //go left to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i - nx;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i + 1;
			val[numE] = -(dx[I]-a[i-1]) / (dx[I] - a[i-1] + b[i-1]);
		}
		else if (tags[i+nx] != -1)//go north to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 + nx;
			val[numE] = 0;
			numE++;

			col[numE] = i - nx;
			val[numE] = -(dx[I]-a[i+nx]) / (dx[I] - a[i+nx] + b[i+nx]);
		}
		else if (tags[i-nx] != -1)//go south to body
		{
			row[numE] = i;
			col[numE] = i - 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = i + 1;
			val[numE] = 0;
			numE++;

			row[numE] = i;
			col[numE] = 1 - nx;
			val[numE] = 0;
			numE++;

			col[numE] = i + nx;
			val[numE] = -(dx[I]-a[i-nx]) / (dx[I] - a[i-nx] + b[i-nx]);
		}
		row[numE] = i;
		numE++;

		//setup p
		col[numE] = i;
		row[numE] = i;
		val[numE] = 1;
		numE++;
	}//end inside tag
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
