/***************************************************************************//**
 * \file projectVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \CPU Author, Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

#include "tagPoints.h"

namespace kernels
{
__global__
void interpolateVelocityX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
							double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints)
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
		iu = J*(nx-1) + I,
		//ip = J*nx + I,//flag not used
		ii= I-5,
		jj = J-5;
	if (iu > J*(nx-1) + I) //return if we're out of bound
		return;
	if (ghostTagsUV[iu]<=0) //return if we're not at an interpolation point
		return;

	/*
	 *   (x1,y1)__________(x2,y2)
	 *   |						|
	 *   | 		*(ip_x,ip_y)	|
	 *   |						|
	 *   |						|
	 *   |						|
	 *   (x3,y3)__________(x4,y4)
	 */
	int minl;
	double x1, x2, x3, x4, y1, y2, y3, y4, q1, q2, q3, q4, minS;
	double	a30 = 30,
			a20 = 30,
			a10 = 30,
			a0 = 1,
			a1 = 1,
			a2 = 1,
			a3 = 1,
			count = 0,
			tol = 0.001,
			alpha = 0.9;

	//find the x and y locations that bound the image point
	while (xu[ii] < image_point_x[iu])
		ii++;
	x1 = xu[ii];
	x2 = 1;
	x3 = xu[ii+1];
	x4 = x3;
	while (yu[jj] <image_point_y[iu])
		jj++;
	y1 = yu[jj];
	y2 = y1;
	y3 = yu[jj+1];
	y4 = y3;

	//find q1,q2,q3,q4
	q1 = u[jj*(nx-1)+ii];
	q2 = u[jj*(nx-1)+ii+1];
	q3 = u[(jj+1)*(nx-1)+ii];
	q4 = u[(jj+1)*(nx-1)+ii+1];
	//check if any points are inside of the body
	//point 1
	if (ghostTagsUV[jj*(nx-1)+ii] != -1)
	{
		//find closest body node to the point that needs replacement
		minS = 1;
		for (int l = 0; l < totalPoints; l++)
		{
			if (minS > sqrt(pow(bx[l]-x1,2)+pow(by[l]-y1,2)))
			{
				minS = sqrt(pow(bx[l]-x1,2)+pow(by[l]-y1,2));
				minl = l;
			}
		}
		x1 = bx[minl];
		y1 = by[minl];
		q1 = uB[minl];
	}
	//point 2
	if (ghostTagsUV[jj*(nx-1)+ii+1] != -1)
	{
		//find closest body node to the point that needs replacement
		minS = 1;
		for (int l = 0; l < totalPoints; l++)
		{
			if (minS > sqrt(pow(bx[l]-x2,2)+pow(by[l]-y2,2)))
			{
				minS = sqrt(pow(bx[l]-x2,2)+pow(by[l]-y2,2));
				minl = l;
			}
		}
		x2 = bx[minl];
		y2 = by[minl];
		q2 = uB[minl];
	}
	//point 3
	if (ghostTagsUV[(jj+1)*(nx-1)+ii] != -1)
	{
		//find closest body node to the point that needs replacement
		minS = 1;
		for (int l = 0; l < totalPoints; l++)
		{
			if (minS > sqrt(pow(bx[l]-x3,2)+pow(by[l]-y3,2)))
			{
				minS = sqrt(pow(bx[l]-x3,2)+pow(by[l]-y3,2));
				minl = l;
			}
		}
		x3 = bx[minl];
		y3 = by[minl];
		q3 = uB[minl];
	}
	//point 4
	if (ghostTagsUV[(jj+1)*(nx-1)+ii+1] != -1)
	{
		//find closest body node to the point that needs replacement
		minS = 1;
		for (int l = 0; l < totalPoints; l++)
		{
			if (minS > sqrt(pow(bx[l]-x4,2)+pow(by[l]-y4,2)))
			{
				minS = sqrt(pow(bx[l]-x4,2)+pow(by[l]-y4,2));
				minl = l;
			}
		}
		x4 = bx[minl];
		y4 = by[minl];
		q4 = uB[minl];
	}
	//interpolate
	//http://math.stackexchange.com/questions/828392/spatial-interpolation-for-irregular-grid
	while ( ( ( abs(a3 - a30) < tol ) && abs( (a2 - a20) < tol ) && ( (a1 - a10) < tol ) ) || count > 10000)
	{
		a0 = ((q1     -a1*x1  - a2*y1 - a3*x1*y1)        )* alpha + (1-alpha)*a0;
		a1 = ((q2 -a0         - a2*y2 - a3*x2*y2)/x2     )* alpha + (1-alpha)*a1;
		a2 = ((q3 -a0 - a1*x3         - a3*x3*y3)/y3     )* alpha + (1-alpha)*a2;
		a3 = ((q4 -a0 - a1*x4 - a2*y4           )/(y4*x4))* alpha + (1-alpha)*a3;

		count = count +1;
		a30 = a3;
		a20 = a2;
		a10 = a1;
	}

	//u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
	u[iu] =  2*uB[0] - a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];

}

__global__
void interpolateVelocityY(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
							double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints)
{
}
}
