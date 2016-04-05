/***************************************************************************//**
 * \file N.cu
 * \author Chris Minar
 * \calculates N on the GPU
 */


#include "N.h"


/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{

//input in range from 0 to ~ nx*ny
__global__
void Nmidx(double *N, double *u, double *dx, double *dy, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1),
		iv	= (nx-1)*ny  +  nx*J +I;
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;

	N[i] = u[i]*(
				   u[i+1] * dx[I]/(dx[I+1]*(dx[I]+dx[I+1]))//east
				 + u[i-1] * (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I])//west
				 + u[i]   * (-dx[I]/(dx[I+1]*(dx[I]+dx[I+1])) - (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I]))//center
				) +
		  (((u[iv+1]-u[iv])*dx[I]/(dx[I]+dx[I+1]) + u[iv]) / 2 + ((u[iv-nx+1]-u[iv-nx])*dx[I]/(dx[I]+dx[I+1]) + u[iv-nx]) / 2) //v
		       *(
			       u[i+nx-1] * dy[J]/(dy[J]*(dy[J+1]+dy[J]))//North
			     + u[i-nx+1] * (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J])//South
			     + u[i]      * (-dy[J]/(dy[J]*(dy[J+1]+dy[J])) - (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J]))//more center
		        );
}

//
__global__
void Nbcx(double *N, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;

	int i 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= i % (nx-1),
		J	= i / (nx-1),
		iv	= (nx-1)*ny  +  nx*J +I;
	if (I != 0 && I != nx-2 && J != 0 && J != ny-1)
		return;

	double temp = 0;


	//East
	if (I == nx-2)
		temp += u[i]*(
				   xp[J]  * dx[I]/(dx[I+1]*(dx[I]+dx[I+1]))//east
				 + u[i-1] * (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I])//west
				 + u[i]   * (-dx[I]/(dx[I+1]*(dx[I]+dx[I+1])) - (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I]))//center
				);
	//West
	else if(I == 0)
		temp += u[i]*(
				   u[i+1] * dx[I]/(dx[I+1]*(dx[I]+dx[I+1]))//east
				 + xm[J]  * (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I])//west
				 + u[i]   * (-dx[I]/(dx[I+1]*(dx[I]+dx[I+1])) - (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I]))//center
				);
	//E-W center
	else
		temp += u[i]*(
				   u[i+1] * dx[I]/(dx[I+1]*(dx[I]+dx[I+1]))//east
				 + u[i-1] * (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I])//west
				 + u[i]   * (-dx[I]/(dx[I+1]*(dx[I]+dx[I+1])) - (dx[I]/(dx[I]*(dx[I]+dx[I+1])) - 1/dx[I]))//center
				);

	//North
	if(J == ny-1)
		temp += (((yp[(nx-1)+I+1]-yp[(nx-1)+I])*dx[I]/(dx[I]+dx[I+1]) + yp[(nx-1)+I]) / 2 + ((u[iv-nx+1]-u[iv-nx])*dx[I]/(dx[I]+dx[I+1]) + u[iv-nx]) / 2) //v
			   *(
				   (2*yp[I] - u[i]) * dy[J]   /(dy[J]*(dy[J]+dy[J]))//North
				 + u[i-nx+1]        * (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J])//South
				 + u[i]             * (-dy[J] /(dy[J]*(dy[J]+dy[J])) - (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J]))//more center
				);
	//South
	else if(J == 0)
		temp += (((u[iv+1]-u[iv])*dx[I]/(dx[I]+dx[I+1]) + u[iv]) / 2 + ((ym[(nx-1)+I+1]-ym[(nx-1)+I])*dx[I]/(dx[I]+dx[I+1]) + ym[(nx-1)+I]) / 2) //v
			   *(
				   u[i+nx-1]      * dy[J]   /(dy[J]*(dy[J+1]+dy[J]))//North
				 + (2*ym[I]-u[i]) * (dy[J]  /(dy[J]*(dy[J]+dy[J])) -1/dy[J])//South
				 + u[i]           * (-dy[J] /(dy[J]*(dy[J+1]+dy[J])) - (dy[J]/(dy[J]*(dy[J]+dy[J])) -1/dy[J]))//more center
				);
	//N-S center
	else
		temp += (((u[iv+1]-u[iv])*dx[I]/(dx[I]+dx[I+1]) + u[iv]) / 2 + ((u[iv-nx+1]-u[iv-nx])*dx[I]/(dx[I]+dx[I+1]) + u[i-nx]) / 2) //v
			   *(
				   u[i+nx-1] * dy[J]/(dy[J]*(dy[J+1]+dy[J]))//North
				 + u[i-nx+1] * (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J])//South
				 + u[i]      * (-dy[J]/(dy[J]*(dy[J+1]+dy[J])) - (dy[J-1]/(dy[J]*(dy[J]+dy[J-1])) -1/dy[J]))//more center
				);

	N[i] = temp;
}

//input in range from 0 to ~ nx*ny
__global__
void Nmidy(double *N, double *u, double *dx, double *dy, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;

	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iu	= (nx-1)*J + I,
		iv = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;

	N[iv] = (((u[iu+(nx-1)]-u[iu])*dy[J]/(dy[J]+dy[J+1]) + u[iu]) / 2 + ((u[iu+(nx-1)-1]-u[iu-1])*dy[J]/(dy[J]+dy[J+1]) + u[iu-1]) / 2) //u
			   *(
				    u[iv+1] * dx[I]/(dx[I]*(dx[I]+dx[I+1]))//east
			 	  + u[iv-1] * (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I])//west
				  + u[iv]   * (-dx[I]/(dx[I]*(dx[I]+dx[I+1])) - (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I]))//center
				 )+
			u[iv]
			    *(
				    u[iv+nx] * dy[J]/(dy[J+1]*(dy[J+1]+dy[J]))//North
				  + u[iv-nx] *(dy[J]/(dy[J]  *(dy[J+1]+dy[J])) - 1/dy[J])//South
				  + u[iv]    * (-dy[J]/(dy[J+1]*(dy[J+1]+dy[J])) - (dy[J]/(dy[J]*(dy[J+1]+dy[J])) - 1/dy[J]))//more center
				 );
}

__global__
void Nbcy(double *N, double *u, double *dx, double *dy, double *ym, double *yp, double *xm, double *xp, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;

	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iu	= (nx-1)*J + I,
		iv = ip + (nx-1)*ny;

	if (I != 0 && I != nx-1 && J != 0 && J != ny-2)
			return;

	double temp = 0;

	//East
	if (I == nx-1)
		temp += (((xp[J+1]-xp[J])*dy[J]/(dy[J]+dy[J+1]) + xp[J]) / 2 + ((u[iu+(nx-1)-1]-u[iu-1])*dy[J]/(dy[J]+dy[J+1]) + u[iu-1]) / 2) //u
			   *(
					(2*xp[ny+J] - u[iv]) * dx[I]   /(dx[I]*(dx[I]+dx[I]))//east
				  + u[iv-1]              * (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I])//west
				  + u[iv]                * (-dx[I] /(dx[I]*(dx[I]+dx[I])) - (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I]))//center
				 );
	//West
	else if(I == 0)
		temp += (((u[iu+(nx-1)]-u[iu])*dy[J]/(dy[J]+dy[J+1]) + u[iu]) / 2 + ((xm[J+1]-xm[J])*dy[J]/(dy[J]+dy[J+1]) + xm[J]) / 2) //u
			   *(
					u[iv+1]              * dx[I]  /(dx[I]*(dx[I]+dx[I+1]))//east
				  + (2*xm[ny+J] - u[iv]) * (dx[I] /(dx[I]*(dx[I]+dx[I])) - 1/dx[I])//west
				  + u[iv]                * (-dx[I]/(dx[I]*(dx[I]+dx[I+1])) - (dx[I]/(dx[I]*(dx[I]+dx[I])) - 1/dx[I]))//center
				 );
	//E-W center
	else
		temp += (((u[iu+(nx-1)]-u[iu])*dy[J]/(dy[J]+dy[J+1]) + u[iu]) / 2 + ((u[iu+(nx-1)-1]-u[iu-1])*dy[J]/(dy[J]+dy[J+1]) + u[iu-1]) / 2) //u
			   *(
					u[iv+1] * dx[I]/(dx[I]*(dx[I]+dx[I+1]))//east
				  + u[iv-1] * (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I])//west
				  + u[iv]   * (-dx[I]/(dx[I]*(dx[I]+dx[I+1])) - (dx[I-1]/(dx[I]*(dx[I]+dx[I-1])) - 1/dx[I]))//center
				 );

	//North
	if(J == ny-2)
		temp += u[iv]
				    *(
					    yp[(nx-1)+I] * dy[J]/(dy[J+1]*(dy[J+1]+dy[J]))//North
					  + u[iv-nx] *(dy[J]/(dy[J]  *(dy[J+1]+dy[J])) - 1/dy[J])//South
					  + u[iv]    * (-dy[J]/(dy[J+1]*(dy[J+1]+dy[J])) - (dy[J]/(dy[J]*(dy[J+1]+dy[J])) - 1/dy[J]))//more center
					 );
	//South
	else if(J == 0)
		temp += u[iv]
				    *(
					    u[iv+nx] * dy[J]/(dy[J+1]*(dy[J+1]+dy[J]))//North
					  + ym[(nx-1)+I] *(dy[J]/(dy[J]  *(dy[J+1]+dy[J])) - 1/dy[J])//South
					  + u[iv]    * (-dy[J]/(dy[J+1]*(dy[J+1]+dy[J])) - (dy[J]/(dy[J]*(dy[J+1]+dy[J])) - 1/dy[J]))//more center
					 );
	//N-S center
	else
		temp += u[iv]
				    *(
					    u[iv+nx] * dy[J]/(dy[J+1]*(dy[J+1]+dy[J]))//North
					  + u[iv-nx] *(dy[J]/(dy[J]  *(dy[J+1]+dy[J])) - 1/dy[J])//South
					  + u[iv]    * (-dy[J]/(dy[J+1]*(dy[J+1]+dy[J])) - (dy[J]/(dy[J]*(dy[J+1]+dy[J])) - 1/dy[J]))//more center
					 );

	N[iv] = temp;
}
}
