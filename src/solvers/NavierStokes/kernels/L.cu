#include "L.h"

namespace kernels
{

//input in range from 0 to ~ nx*ny
__global__
void Lmidx(double *L, double *u, double *dx, double *dy, int nx, int ny, double nu)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;


	L[i] = nu*(
				 (u[i+1]     -u[i]) / (dx[I+1]*(dx[I+1]+dx[I])*0.5)//east
				+(u[i-1]     -u[i]) / (dx[I]  *(dx[I+1]+dx[I])*0.5)//west
				+(u[i+(nx-1)]-u[i]) / (dy[J]  *(dy[J+1]+dy[J])*0.5)//north
				+(u[i-(nx-1)]-u[i]) / (dy[J]  *(dy[J-1]+dy[J])*0.5)//south
			   );
}

__global__
void Lbcx(double *L, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny, double nu)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;

	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1);

	if (I != 0 && I != nx-2 && J != 0 && J != ny-1)
			return;

	double temp = 0;

	//East
	if(I != nx-2)
		temp += nu * (u[i+1]     -u[i]) / (dx[I+1]*(dx[I+1]+dx[I])*0.5);
	//East Boundary
	else
		temp += nu * (xp[J]      -u[i]) / (dx[I+1]*(dx[I+1]+dx[I])*0.5);

	//West
	if(I != 0)
		temp += nu * (u[i-1]     -u[i]) / (dx[I]  *(dx[I+1]+dx[I])*0.5);
	//West Boundary
	else
		temp += nu * (xm[J]      -u[i]) / (dx[I]  *(dx[I+1]+dx[I])*0.5);

	//North
	if(J != ny-1)
		temp += nu * (u[i+(nx-1)]-u[i]) / (dy[J]  *(dy[J+1]+dy[J])*0.5);
	//North Boundary
	else
		temp += nu * (2*yp[I]  -2*u[i]) / (dy[J] * dy[J]);

	//South
	if(J != 0)
		temp += nu * (u[i-(nx-1)]-u[i]) / (dy[J]  *(dy[J-1]+dy[J])*0.5);
	//South Boundary
	else
		temp += nu * (2*ym[I]  -2*u[i]) / (dy[J] * dy[J]);

	L[i] = temp;
}

//input in range from 0 to ~ nx*ny
__global__
void Lmidy(double *L, double *u, double *dx, double *dy, int nx, int ny, double nu)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;

	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iv = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;


	L[iv] = nu* (
				  (u[iv+1]  -u[iv]) / (dx[I]  *(dx[I]+dx[I+1])*0.5)//east
				 +(u[iv-1]  -u[iv]) / (dx[I]  *(dx[I]+dx[I-1])*0.5)//west
				 +(u[iv+nx] -u[iv]) / (dy[J+1]*(dy[J]+dy[J+1])*0.5)//north
				 +(u[iv-nx] -u[iv]) / (dy[J]  *(dy[J]+dy[J+1])*0.5)//south
				 );
}

__global__
void Lbcy(double *L, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny, double nu)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;

	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iv = ip + (nx-1)*ny;

	if (I != 0 && I != nx-1 && J != 0 && J != ny-2)
			return;

	double temp = 0;

	//East
	if(I != nx-1)
		temp += nu * (u[iv+1]       -u[iv]) / (dx[I]  *(dx[I]+dx[I+1])*0.5);
	//East Boundary
	else
		temp += nu * (2*xp[ny+J]  -2*u[iv]) / (dx[I] * dx[I]);

	//West
	if(I != 0)
		temp += nu * (u[iv-1]       -u[iv]) / (dx[I]  *(dx[I]+dx[I-1])*0.5);
	//West Boundary
	else
		temp += nu * (2*xm[ny+J]  -2*u[iv]) / (dx[I]  *dx[I]);

	//North
	if(J != ny-2)
		temp += nu * (u[iv+nx]      -u[iv]) / (dy[J+1]*(dy[J]+dy[J+1])*0.5);
	//North Boundary
	else
		temp += nu * (yp[nx-1+I]    -u[iv]) / (dy[J+1]*(dy[J]+dy[J+1])*0.5);

	//South
	if(J != 0)
		temp += nu * (u[iv-nx]      -u[iv]) / (dy[J]  *(dy[J]+dy[J+1])*0.5);
	//South Boundary
	else
		temp += nu * (ym[nx-1+I]    -u[iv]) / (dy[J]  *(dy[J]+dy[J+1])*0.5);

	L[iv] = temp;
}
}
