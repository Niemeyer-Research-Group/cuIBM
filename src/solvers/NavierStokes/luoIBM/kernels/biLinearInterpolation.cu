/***************************************************************************//**
 * \file .cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \CPU Author, Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

#include "tagPoints.h"

namespace kernels
{
__global__
void interpolateVelocityToGhostNodeX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int *i_start, int *j_start, int width, int nx, int ny,
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u)//testing variables
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iu = J*(nx-1) + I,
		ii= I-5,
		jj = J-5;
	if (iu > J*(nx-1) + I) //return if we're out of bound
		return;
	if (ghostTagsUV[iu]<=0) //return if we're not at an interpolation point
		return;

	/*
	 *   (x3,y3)__________(x4,y4)
	 *   |						|
	 *   | 		*(ip_x,ip_y)	|
	 *   |						|
	 *   |						|
	 *   |						|
	 *   (x1,y1)__________(x2,y2)
	 */

	//find x and y of nodes that bound the image point
	while (xu[ii] < image_point_x[iu])
		ii++;
	x1[iu] = xu[ii-1];
	x2[iu] = xu[ii];
	x3[iu] = x1[iu];
	x4[iu] = x2[iu];
	while (yu[jj] <image_point_y[iu])
		jj++;
	y1[iu] = yu[jj-1];
	y2[iu] = y1[iu];
	y3[iu] = yu[jj];
	y4[iu] = y3[iu];

	//find q1,q2,q3,q4
	q1[iu] = u[(jj-1)*(nx-1)+ii-1];
	q2[iu] = u[(jj-1)*(nx-1)+ii];
	q3[iu] = u[jj*(nx-1)+ii-1];
	q4[iu] = u[jj*(nx-1)+ii];
	//check if any points are inside of the body, then move them to the body intercept
	//point 1
	if (ghostTagsUV[(jj-1)*(nx-1)+ii-1] > 0)
	{
		x1[iu] = body_intercept_x[(jj-1)*(nx-1)+ii-1];
		y1[iu] = body_intercept_y[(jj-1)*(nx-1)+ii-1];
		q1[iu] = uB[0];
	}
	if (ghostTagsUV[(jj-1)*(nx-1)+ii] > 0)
	{
		x2[iu] = body_intercept_x[(jj-1)*(nx-1)+ii];
		y2[iu] = body_intercept_y[(jj-1)*(nx-1)+ii];
		q2[iu] = uB[0];
	}
	if (ghostTagsUV[jj*(nx-1)+ii-1] > 0)
	{
		x3[iu] = body_intercept_x[jj*(nx-1)+ii-1];
		y3[iu] = body_intercept_y[jj*(nx-1)+ii-1];
		q3[iu] = uB[0];
	}
	if (ghostTagsUV[jj*(nx-1)+ii] > 0)
	{
		x4[iu] = body_intercept_x[jj*(nx-1)+ii];
		y4[iu] = body_intercept_y[jj*(nx-1)+ii];
		q4[iu] = uB[0];
	}
	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */
	double a12 = x1[iu],  a13 = y1[iu], a14 = x1[iu]*y1[iu];
	double a22 = x2[iu],  a23 = y2[iu], a24 = x2[iu]*y2[iu];
	double a32 = x3[iu],  a33 = y3[iu], a34 = x3[iu]*y3[iu];
	double a42 = x4[iu],  a43 = y4[iu], a44 = x4[iu]*y4[iu];

	double
	detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
	      +a12*1*a34*a43 + a12*a23*1*a44 + a12*a24*a33*1
	      +a13*1*a32*a44 + a13*a22*a34*1 + a13*a24*1*a42
	      +a14*1*a33*a42 + a14*a22*1*a43 + a14*a23*a32*1
	      -1*a22*a34*a43 - 1*a23*a32*a44 - 1*a24*a33*a42
	      -a12*1*a33*a44 - a12*a23*a34*1 - a12*a24*1*a43
	      -a13*1*a34*a42 - a13*a22*1*a44 - a13*a24*a32*1
	      -a14*1*a32*a43 - a14*a22*a33*1 - a14*a23*1*a42;

	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = 1*a34*a43 + a23*1*a44 + a24*a33*1 - 1*a33*a44 - a23*a34*1 - a24*1*a43;
	double b22 = 1*a33*a44 + a13*a34*1 + a14*1*a43 - 1*a34*a43 - a13*1*a44 - a14*a33*1;
	double b23 = 1*a24*a43 + a13*1*a44 + a14*a23*1 - 1*a23*a44 - a13*a24*1 - a14*1*a43;
	double b24 = 1*a23*a34 + a13*a24*1 + a14*1*a33 - 1*a24*a33 - a13*1*a34 - a14*a23*1;
	double b31 = 1*a32*a44 + a22*a34*1 + a24*1*a42 - 1*a34*a42 - a22*1*a44 - a24*a32*1;
	double b32 = 1*a34*a42 + a12*1*a44 + a14*a32*1 - 1*a32*a44 - a12*a34*1 - a14*1*a42;
	double b33 = 1*a22*a44 + a12*a24*1 + a14*1*a42 - 1*a24*a42 - a12*1*a44 - a14*a22*1;
	double b34 = 1*a24*a32 + a12*1*a34 + a14*a22*1 - 1*a22*a34 - a12*a24*1 - a14*1*a32;
	double b41 = 1*a33*a42 + a22*1*a43 + a23*a32*1 - 1*a32*a43 - a22*a33*1 - a23*1*a42;
	double b42 = 1*a32*a43 + a12*a33*1 + a13*1*a42 - 1*a33*a42 - a12*1*a43 - a13*a32*1;
	double b43 = 1*a23*a42 + a12*1*a43 + a13*a22*1 - 1*a22*a43 - a12*a23*1 - a13*1*a42;
	double b44 = 1*a22*a33 + a12*a23*1 + a13*1*a32 - 1*a23*a32 - a12*1*a33 - a13*a22*1;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	 double a0 = b11/detA*q1[iu]  +  b12/detA*q2[iu]  +  b13/detA*q3[iu]  +  b14/detA*q4[iu];
	 double a1 = b21/detA*q1[iu]  +  b22/detA*q2[iu]  +  b23/detA*q3[iu]  +  b24/detA*q4[iu];
	 double a2 = b31/detA*q1[iu]  +  b32/detA*q2[iu]  +  b33/detA*q3[iu]  +  b34/detA*q4[iu];
	 double a3 = b41/detA*q1[iu]  +  b42/detA*q2[iu]  +  b43/detA*q3[iu]  +  b44/detA*q4[iu];
	 image_point_u[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
	 u[iu] =  2*uB[0] - image_point_u[iu]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToGhostNodeY(double *u, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int *i_start, int *j_start, int width, int nx, int ny,
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u)//testing variables
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iv = J*nx + I + (nx-1)*ny,
		ii= I-5,
		jj = J-5;
	if (J*nx + I > nx*(ny-1)) //return if we're out of bound
		return;
	if (ghostTagsUV[iv]<=0) //return if we're not at an interpolation point
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

	//find x and y of nodes that bound the image point
	while (xv[ii] < image_point_x[iv])
		ii++;
	x1[iv] = xv[ii-1];
	x2[iv] = xv[ii];
	x3[iv] = x1[iv];
	x4[iv] = x2[iv];
	while (yv[jj] <image_point_y[iv])
		jj++;
	y1[iv] = yv[jj-1];
	y2[iv] = y1[iv];
	y3[iv] = yv[jj];
	y4[iv] = y3[iv];

	//find q1,q2,q3,q4
	q1[iv] = u[(jj-1)*nx+ii-1+ (nx-1)*ny];
	q2[iv] = u[(jj-1)*nx+ii+ (nx-1)*ny];
	q3[iv] = u[jj*nx+ii-1+ (nx-1)*ny];
	q4[iv] = u[jj*nx+ii+ (nx-1)*ny];
	//check if any points are inside of the body, then move them to the body intercept
	//point 1
	if (ghostTagsUV[(jj-1)*nx+ii-1 + (nx-1)*ny] > 0)
	{
		x1[iv] = body_intercept_x[(jj-1)*nx+ii-1+ (nx-1)*ny];
		y1[iv] = body_intercept_y[(jj-1)*nx+ii-1+ (nx-1)*ny];
		q1[iv] = vB[0];
	}
	if (ghostTagsUV[(jj-1)*nx+ii+ (nx-1)*ny] > 0)
	{
		x2[iv] = body_intercept_x[(jj-1)*nx+ii+ (nx-1)*ny];
		y2[iv] = body_intercept_y[(jj-1)*nx+ii+ (nx-1)*ny];
		q2[iv] = vB[0];
	}
	if (ghostTagsUV[jj*nx+ii-1+ (nx-1)*ny] > 0)
	{
		x3[iv] = body_intercept_x[jj*nx+ii-1+ (nx-1)*ny];
		y3[iv] = body_intercept_y[jj*nx+ii-1+ (nx-1)*ny];
		q3[iv] = vB[0];
	}
	if (ghostTagsUV[jj*nx+ii+ (nx-1)*ny] > 0)
	{
		x4[iv] = body_intercept_x[jj*nx+ii+ (nx-1)*ny];
		y4[iv] = body_intercept_y[jj*nx+ii+ (nx-1)*ny];
		q4[iv] = vB[0];
	}
	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */
	double a12 = x1[iv],  a13 = y1[iv], a14 = x1[iv]*y1[iv];
	double a22 = x2[iv],  a23 = y2[iv], a24 = x2[iv]*y2[iv];
	double a32 = x3[iv],  a33 = y3[iv], a34 = x3[iv]*y3[iv];
	double a42 = x4[iv],  a43 = y4[iv], a44 = x4[iv]*y4[iv];

	double
	detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
	      +a12*1*a34*a43 + a12*a23*1*a44 + a12*a24*a33*1
	      +a13*1*a32*a44 + a13*a22*a34*1 + a13*a24*1*a42
	      +a14*1*a33*a42 + a14*a22*1*a43 + a14*a23*a32*1
	      -1*a22*a34*a43 - 1*a23*a32*a44 - 1*a24*a33*a42
	      -a12*1*a33*a44 - a12*a23*a34*1 - a12*a24*1*a43
	      -a13*1*a34*a42 - a13*a22*1*a44 - a13*a24*a32*1
	      -a14*1*a32*a43 - a14*a22*a33*1 - a14*a23*1*a42;

	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = 1*a34*a43 + a23*1*a44 + a24*a33*1 - 1*a33*a44 - a23*a34*1 - a24*1*a43;
	double b22 = 1*a33*a44 + a13*a34*1 + a14*1*a43 - 1*a34*a43 - a13*1*a44 - a14*a33*1;
	double b23 = 1*a24*a43 + a13*1*a44 + a14*a23*1 - 1*a23*a44 - a13*a24*1 - a14*1*a43;
	double b24 = 1*a23*a34 + a13*a24*1 + a14*1*a33 - 1*a24*a33 - a13*1*a34 - a14*a23*1;
	double b31 = 1*a32*a44 + a22*a34*1 + a24*1*a42 - 1*a34*a42 - a22*1*a44 - a24*a32*1;
	double b32 = 1*a34*a42 + a12*1*a44 + a14*a32*1 - 1*a32*a44 - a12*a34*1 - a14*1*a42;
	double b33 = 1*a22*a44 + a12*a24*1 + a14*1*a42 - 1*a24*a42 - a12*1*a44 - a14*a22*1;
	double b34 = 1*a24*a32 + a12*1*a34 + a14*a22*1 - 1*a22*a34 - a12*a24*1 - a14*1*a32;
	double b41 = 1*a33*a42 + a22*1*a43 + a23*a32*1 - 1*a32*a43 - a22*a33*1 - a23*1*a42;
	double b42 = 1*a32*a43 + a12*a33*1 + a13*1*a42 - 1*a33*a42 - a12*1*a43 - a13*a32*1;
	double b43 = 1*a23*a42 + a12*1*a43 + a13*a22*1 - 1*a22*a43 - a12*a23*1 - a13*1*a42;
	double b44 = 1*a22*a33 + a12*a23*1 + a13*1*a32 - 1*a23*a32 - a12*1*a33 - a13*a22*1;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	 double a0 = b11/detA*q1[iv]  +  b12/detA*q2[iv]  +  b13/detA*q3[iv]  +  b14/detA*q4[iv];
	 double a1 = b21/detA*q1[iv]  +  b22/detA*q2[iv]  +  b23/detA*q3[iv]  +  b24/detA*q4[iv];
	 double a2 = b31/detA*q1[iv]  +  b32/detA*q2[iv]  +  b33/detA*q3[iv]  +  b34/detA*q4[iv];
	 double a3 = b41/detA*q1[iv]  +  b42/detA*q2[iv]  +  b43/detA*q3[iv]  +  b44/detA*q4[iv];
	 image_point_u[iv] = a0 + a1*image_point_x[iv] + a2*image_point_y[iv] + a3*image_point_x[iv]*image_point_y[iv];
	 u[iv] =  2*vB[0] - image_point_u[iv]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u)//test
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iu = J*(nx-1) + I,
		ii= I-5,
		jj = J-5;
	if (iu > J*(nx-1) + I) //return if we're out of bound
		return;
	if (hybridTagsUV[iu]<=0) //return if we're not at an interpolation point
		return;

	/*
	 *   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *   	|					   |
	 *   	(x1,y1)__________(x2,y2)
	 *
	 *   *(BI_x,BI_y)
	 *
	 *
	 */

	//find x and y of nodes that bound the image point
	while (xu[ii] < image_point_x[iu])
		ii++;
	x1[iu] = xu[ii-1];
	x2[iu] = xu[ii];
	x3[iu] = x1[iu];
	x4[iu] = x2[iu];
	while (yu[jj] <image_point_y[iu])
		jj++;
	y1[iu] = yu[jj-1];
	y2[iu] = y1[iu];
	y3[iu] = yu[jj];
	y4[iu] = y3[iu];

	//find q1,q2,q3,q4
	q1[iu] = u[(jj-1)*(nx-1)+ii-1];
	q2[iu] = u[(jj-1)*(nx-1)+ii];
	q3[iu] = u[jj*(nx-1)+ii-1];
	q4[iu] = u[jj*(nx-1)+ii];

	//check if any points are inside of the body, then move them to the body intercept
	//point 1
	if (hybridTagsUV[(jj-1)*(nx-1)+ii-1] == iu)
	{
		x1[iu] = body_intercept_x[iu];
		y1[iu] = body_intercept_y[iu];
		q1[iu] = uB[0];
	}
	if (hybridTagsUV[(jj-1)*(nx-1)+ii] == iu)
	{
		x2[iu] = body_intercept_x[iu];
		y2[iu] = body_intercept_y[iu];
		q2[iu] = uB[0];
	}
	if (hybridTagsUV[jj*(nx-1)+ii-1] == iu)
	{
		x3[iu] = body_intercept_x[iu];
		y3[iu] = body_intercept_y[iu];
		q3[iu] = uB[0];
	}
	if (hybridTagsUV[jj*(nx-1)+ii] == iu)
	{
		x4[iu] = body_intercept_x[iu];
		y4[iu] = body_intercept_y[iu];
		q4[iu] = uB[0];
	}

	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */
	double a12 = x1[iu],  a13 = y1[iu], a14 = x1[iu]*y1[iu];
	double a22 = x2[iu],  a23 = y2[iu], a24 = x2[iu]*y2[iu];
	double a32 = x3[iu],  a33 = y3[iu], a34 = x3[iu]*y3[iu];
	double a42 = x4[iu],  a43 = y4[iu], a44 = x4[iu]*y4[iu];

	double
	detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
		  +a12*1*a34*a43 + a12*a23*1*a44 + a12*a24*a33*1
		  +a13*1*a32*a44 + a13*a22*a34*1 + a13*a24*1*a42
		  +a14*1*a33*a42 + a14*a22*1*a43 + a14*a23*a32*1
		  -1*a22*a34*a43 - 1*a23*a32*a44 - 1*a24*a33*a42
		  -a12*1*a33*a44 - a12*a23*a34*1 - a12*a24*1*a43
		  -a13*1*a34*a42 - a13*a22*1*a44 - a13*a24*a32*1
		  -a14*1*a32*a43 - a14*a22*a33*1 - a14*a23*1*a42;

	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = 1*a34*a43 + a23*1*a44 + a24*a33*1 - 1*a33*a44 - a23*a34*1 - a24*1*a43;
	double b22 = 1*a33*a44 + a13*a34*1 + a14*1*a43 - 1*a34*a43 - a13*1*a44 - a14*a33*1;
	double b23 = 1*a24*a43 + a13*1*a44 + a14*a23*1 - 1*a23*a44 - a13*a24*1 - a14*1*a43;
	double b24 = 1*a23*a34 + a13*a24*1 + a14*1*a33 - 1*a24*a33 - a13*1*a34 - a14*a23*1;
	double b31 = 1*a32*a44 + a22*a34*1 + a24*1*a42 - 1*a34*a42 - a22*1*a44 - a24*a32*1;
	double b32 = 1*a34*a42 + a12*1*a44 + a14*a32*1 - 1*a32*a44 - a12*a34*1 - a14*1*a42;
	double b33 = 1*a22*a44 + a12*a24*1 + a14*1*a42 - 1*a24*a42 - a12*1*a44 - a14*a22*1;
	double b34 = 1*a24*a32 + a12*1*a34 + a14*a22*1 - 1*a22*a34 - a12*a24*1 - a14*1*a32;
	double b41 = 1*a33*a42 + a22*1*a43 + a23*a32*1 - 1*a32*a43 - a22*a33*1 - a23*1*a42;
	double b42 = 1*a32*a43 + a12*a33*1 + a13*1*a42 - 1*a33*a42 - a12*1*a43 - a13*a32*1;
	double b43 = 1*a23*a42 + a12*1*a43 + a13*a22*1 - 1*a22*a43 - a12*a23*1 - a13*1*a42;
	double b44 = 1*a22*a33 + a12*a23*1 + a13*1*a32 - 1*a23*a32 - a12*1*a33 - a13*a22*1;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	 double a0 = b11/detA*q1[iu]  +  b12/detA*q2[iu]  +  b13/detA*q3[iu]  +  b14/detA*q4[iu];
	 double a1 = b21/detA*q1[iu]  +  b22/detA*q2[iu]  +  b23/detA*q3[iu]  +  b24/detA*q4[iu];
	 double a2 = b31/detA*q1[iu]  +  b32/detA*q2[iu]  +  b33/detA*q3[iu]  +  b34/detA*q4[iu];
	 double a3 = b41/detA*q1[iu]  +  b42/detA*q2[iu]  +  b43/detA*q3[iu]  +  b44/detA*q4[iu];
	 ustar[iu] = a0 + a1*xu[I] + a2*yu[J] + a3*yu[J]*xu[I];
	 image_point_u[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
}

__global__
void interpolateVelocityToHybridNodeY(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u)//test
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iv = J*nx + I + (nx-1)*ny,
		ii= I-5,
		jj = J-5;
	if (J*nx + I > nx*(ny-1)) //return if we're out of bound
		return;
	if (hybridTagsUV[iv]<=0) //return if we're not at an interpolation point
		return;

	/*
	 *   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *   	|					   |
	 *   	(x1,y1)__________(x2,y2)
	 *
	 *   *(BI_x,BI_y)
	 *
	 *
	 */

	//find x and y of nodes that bound the image point
	while (xv[ii] < image_point_x[iv])
		ii++;
	x1[iv] = xv[ii-1];
	x2[iv] = xv[ii];
	x3[iv] = x1[iv];
	x4[iv] = x2[iv];
	while (yv[jj] <image_point_y[iv])
		jj++;
	y1[iv] = yv[jj-1];
	y2[iv] = y1[iv];
	y3[iv] = yv[jj];
	y4[iv] = y3[iv];

	//find q1,q2,q3,q4
	q1[iv] = u[(jj-1)*nx+ii-1 + (nx-1)*ny];
	q2[iv] = u[(jj-1)*nx+ii + (nx-1)*ny];
	q3[iv] = u[jj*nx+ii-1 + (nx-1)*ny];
	q4[iv] = u[jj*nx+ii + (nx-1)*ny];

	//check if any points are inside of the body, then move them to the body intercept
	//point 1
	if (hybridTagsUV[(jj-1)*nx+ii-1 + (nx-1)*ny] == iv)
	{
		x1[iv] = body_intercept_x[iv];
		y1[iv] = body_intercept_y[iv];
		q1[iv] = vB[0];
	}
	if (hybridTagsUV[(jj-1)*nx+ii + (nx-1)*ny] == iv)
	{
		x2[iv] = body_intercept_x[iv];
		y2[iv] = body_intercept_y[iv];
		q2[iv] = vB[0];
	}
	if (hybridTagsUV[jj*nx+ii-1 + (nx-1)*ny] == iv)
	{
		x3[iv] = body_intercept_x[iv];
		y3[iv] = body_intercept_y[iv];
		q3[iv] = vB[0];
	}
	if (hybridTagsUV[jj*nx+ii + (nx-1)*ny] == iv)
	{
		x4[iv] = body_intercept_x[iv];
		y4[iv] = body_intercept_y[iv];
		q4[iv] = vB[0];
	}
	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */
	double a12 = x1[iv],  a13 = y1[iv], a14 = x1[iv]*y1[iv];
	double a22 = x2[iv],  a23 = y2[iv], a24 = x2[iv]*y2[iv];
	double a32 = x3[iv],  a33 = y3[iv], a34 = x3[iv]*y3[iv];
	double a42 = x4[iv],  a43 = y4[iv], a44 = x4[iv]*y4[iv];

	double
	detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
		  +a12*1*a34*a43 + a12*a23*1*a44 + a12*a24*a33*1
		  +a13*1*a32*a44 + a13*a22*a34*1 + a13*a24*1*a42
		  +a14*1*a33*a42 + a14*a22*1*a43 + a14*a23*a32*1
		  -1*a22*a34*a43 - 1*a23*a32*a44 - 1*a24*a33*a42
		  -a12*1*a33*a44 - a12*a23*a34*1 - a12*a24*1*a43
		  -a13*1*a34*a42 - a13*a22*1*a44 - a13*a24*a32*1
		  -a14*1*a32*a43 - a14*a22*a33*1 - a14*a23*1*a42;

	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = 1*a34*a43 + a23*1*a44 + a24*a33*1 - 1*a33*a44 - a23*a34*1 - a24*1*a43;
	double b22 = 1*a33*a44 + a13*a34*1 + a14*1*a43 - 1*a34*a43 - a13*1*a44 - a14*a33*1;
	double b23 = 1*a24*a43 + a13*1*a44 + a14*a23*1 - 1*a23*a44 - a13*a24*1 - a14*1*a43;
	double b24 = 1*a23*a34 + a13*a24*1 + a14*1*a33 - 1*a24*a33 - a13*1*a34 - a14*a23*1;
	double b31 = 1*a32*a44 + a22*a34*1 + a24*1*a42 - 1*a34*a42 - a22*1*a44 - a24*a32*1;
	double b32 = 1*a34*a42 + a12*1*a44 + a14*a32*1 - 1*a32*a44 - a12*a34*1 - a14*1*a42;
	double b33 = 1*a22*a44 + a12*a24*1 + a14*1*a42 - 1*a24*a42 - a12*1*a44 - a14*a22*1;
	double b34 = 1*a24*a32 + a12*1*a34 + a14*a22*1 - 1*a22*a34 - a12*a24*1 - a14*1*a32;
	double b41 = 1*a33*a42 + a22*1*a43 + a23*a32*1 - 1*a32*a43 - a22*a33*1 - a23*1*a42;
	double b42 = 1*a32*a43 + a12*a33*1 + a13*1*a42 - 1*a33*a42 - a12*1*a43 - a13*a32*1;
	double b43 = 1*a23*a42 + a12*1*a43 + a13*a22*1 - 1*a22*a43 - a12*a23*1 - a13*1*a42;
	double b44 = 1*a22*a33 + a12*a23*1 + a13*1*a32 - 1*a23*a32 - a12*1*a33 - a13*a22*1;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	 double a0 = b11/detA*q1[iv]  +  b12/detA*q2[iv]  +  b13/detA*q3[iv]  +  b14/detA*q4[iv];
	 double a1 = b21/detA*q1[iv]  +  b22/detA*q2[iv]  +  b23/detA*q3[iv]  +  b24/detA*q4[iv];
	 double a2 = b31/detA*q1[iv]  +  b32/detA*q2[iv]  +  b33/detA*q3[iv]  +  b34/detA*q4[iv];
	 double a3 = b41/detA*q1[iv]  +  b42/detA*q2[iv]  +  b43/detA*q3[iv]  +  b44/detA*q4[iv];
	 ustar[iv] = a0 + a1*xv[I] + a2*yv[J] + a3*yv[J]*xv[I];
	 image_point_u[iv] = a0 + a1*image_point_x[iv] + a2*image_point_y[iv] + a3*image_point_x[iv]*image_point_y[iv];
}

__global__
void interpolatePressureToHybridNode(double *pressure, double *pressureStar, double *u, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny, double dt, double totalPoints,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, int timeStep)//test
{//flag u not used anymore
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		ip = J*nx + I,
		ii= I-5,
		jj = J-5;
	if (ip > J*nx + I) //return if we're out of bound
		return;
	if (hybridTagsP[ip]<=0) //return if we're not at an interpolation point
		return;

	double	n_x,
			n_y,
			nl;
	/*
	 *   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *   	|					   |
	 *   	(x1,y1)__________(x2,y2)
	 *
	 *   *(BI_x,BI_y)
	 *
	 *
	 */

	//find x and y of nodes that bound the image point
	while (xv[ii] < image_point_p_x[ip])
		ii++;
	x1[ip] = xv[ii-1];
	x2[ip] = xv[ii];
	x3[ip] = x1[ip];
	x4[ip] = x2[ip];
	while (yu[jj] <image_point_p_y[ip])
		jj++;
	y1[ip] = yu[jj-1];
	y2[ip] = y1[ip];
	y3[ip] = yu[jj];
	y4[ip] = y3[ip];

	//find q1,q2,q3,q4
	int	index1 = (jj-1)*nx+ii-1,
		index2 = (jj-1)*nx+ii,
		index3 = jj*nx+ii-1,
		index4 = jj*nx+ii;

	q1[ip] = pressure[index1];
	q2[ip] = pressure[index2];
	q3[ip] = pressure[index3];
	q4[ip] = pressure[index4];

	double a11 = 1, a12 = x1[ip],  a13 = y1[ip], a14 = x1[ip]*y1[ip];
	double a21 = 1, a22 = x2[ip],  a23 = y2[ip], a24 = x2[ip]*y2[ip];
	double a31 = 1, a32 = x3[ip],  a33 = y3[ip], a34 = x3[ip]*y3[ip];
	double a41 = 1, a42 = x4[ip],  a43 = y4[ip], a44 = x4[ip]*y4[ip];

	//setup for neuman BC
	double distance, distance2, min, min2, matDi, matDj, matD2i, matD2j, matDBIi, matDBIj;
	int bodyindex, bodyindex2;
	//move the closes node to the body to the surface then calculate the neuman boundary condition for it
	//point 1
	if (hybridTagsP[index1] == ip)
	{
		//setup
		x1[ip] = body_intercept_p_x[ip];
		y1[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x1[ip];
		n_y = image_point_p_y[ip] - y1[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x1[ip],2)+pow(by[k]-y1[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x1[ip],2)+pow(by[k]-y1[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		q1[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[ip]+a12*y1[ip];
	}
	//point 2
	else if (hybridTagsP[index2] == ip)
	{
		x2[ip] = body_intercept_p_x[ip];
		y2[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x2[ip];
		n_y = image_point_p_y[ip] - y2[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x2[ip],2)+pow(by[k]-y2[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x2[ip],2)+pow(by[k]-y2[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		q2[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[ip]+a22*y2[ip];
	}
	//point 3
	else if (hybridTagsP[index3] == ip)
	{
		x3[ip] = body_intercept_p_x[ip];
		y3[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x3[ip];
		n_y = image_point_p_y[ip] - y3[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x3[ip],2)+pow(by[k]-y3[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x3[ip],2)+pow(by[k]-y3[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		q3[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x3[ip]+a32*y3[ip];
	}
	//4
	if (hybridTagsP[index4] == ip)
	{
		x4[ip] = body_intercept_p_x[ip];
		y4[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x4[ip];
		n_y = image_point_p_y[ip] - y4[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x4[ip],2)+pow(by[k]-y4[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x4[ip],2)+pow(by[k]-y4[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		q4[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a41 = 0;
		a42 = n_x/nl;
		a43 = n_y/nl;
		a44 = a43*x4[ip]+a42*y4[ip];
	}
	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html  //for solving a 4x4 matrix exactly
	//https://www.physicsforums.com/threads/is-normal-derivative-a-definition.706458/   //for dealing with the normal at the boundary
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */

	double
	detA = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43
		  +a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41
		  +a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42
		  +a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41
		  -a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42
		  -a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43
		  -a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41
		  -a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;

	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
	double b22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
	double b23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
	double b24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
	double b31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
	double b32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
	double b33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
	double b34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
	double b41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
	double b42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
	double b43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
	double b44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	 a0[ip] = b11/detA*q1[ip]  +  b12/detA*q2[ip]  +  b13/detA*q3[ip]  +  b14/detA*q4[ip];
	 a1[ip] = b21/detA*q1[ip]  +  b22/detA*q2[ip]  +  b23/detA*q3[ip]  +  b24/detA*q4[ip];
	 a2[ip] = b31/detA*q1[ip]  +  b32/detA*q2[ip]  +  b33/detA*q3[ip]  +  b34/detA*q4[ip];
	 a3[ip] = b41/detA*q1[ip]  +  b42/detA*q2[ip]  +  b43/detA*q3[ip]  +  b44/detA*q4[ip];

	//dudt[ip] = du_dt;
	//ududx[ip] = u_du_dx;
	pressureStar[ip] = a0[ip] + a1[ip]*xv[I] + a2[ip]*yu[J] + a3[ip]*xv[I]*yu[J];
}

//flag this function is a mess
__global__
void interpolatePressureToGhostNode(double *pressure, double *u, int *ghostTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y, double *body_intercept_p,
									int *i_start, int *j_start, int width, int nx, int ny, double dt, double totalPoints,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4)//test
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		ip = J*nx + I,
		ii= I-5,
		jj = J-5;
	if (ip > J*nx + I) //return if we're out of bound
		return;
	if (ghostTagsP[ip]<=0) //return if we're not at an interpolation point
		return;

	double	n_x,
			n_y,
			nl,
			matDClose;
	int close_index,
		index1,
		index2,
		index3,
		index4;
	/*
	 *   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *   	|					   |
	 *   	(x1,y1)__________(x2,y2)
	 *
	 *   *(BI_x,BI_y)
	 *
	 *
	 */

	//find x and y of nodes that bound the image point (find points 1, 2, 3 and 4)
	while (xv[ii] < image_point_p_x[ip])
		ii++;
	x1[ip] = xv[ii-1];
	x2[ip] = xv[ii];
	x3[ip] = x1[ip];
	x4[ip] = x2[ip];
	while (yu[jj] <image_point_p_y[ip])
		jj++;
	y1[ip] = yu[jj-1];
	y2[ip] = y1[ip];
	y3[ip] = yu[jj];
	y4[ip] = y3[ip];

	//find q1,q2,q3,q4
	index1 = (jj-1)*nx+ii-1;
	index2 = (jj-1)*nx+ii;
	index3 = jj*nx+ii-1;
	index4 = jj*nx+ii;

	q1[ip] = pressure[index1];
	q2[ip] = pressure[index2];
	q3[ip] = pressure[index3];
	q4[ip] = pressure[index4];

	double a11 = 1, a12 = x1[ip],  a13 = y1[ip], a14 = x1[ip]*y1[ip];
	double a21 = 1, a22 = x2[ip],  a23 = y2[ip], a24 = x2[ip]*y2[ip];
	double a31 = 1, a32 = x3[ip],  a33 = y3[ip], a34 = x3[ip]*y3[ip];
	double a41 = 1, a42 = x4[ip],  a43 = y4[ip], a44 = x4[ip]*y4[ip];

	//find the closest node to the body intercept (e.g. find node 1 in the above image)
	double min = 1.0;
	double s;
	s = sqrt(pow(x1[ip]-body_intercept_p_x[ip],2) + pow(y1[ip]-body_intercept_p_y[ip],2));
	if (s < min)
	{
		min = s;
		close_index = index1;
	}
	s = sqrt(pow(x2[ip]-body_intercept_p_x[ip],2) + pow(y2[ip]-body_intercept_p_y[ip],2));
	if (s < min)
	{
		min = s;
		close_index = index2;
	}
	s = sqrt(pow(x3[ip]-body_intercept_p_x[ip],2) + pow(y3[ip]-body_intercept_p_y[ip],2));
	if (s < min)
	{
		min = s;
		close_index = index3;
	}
	s = sqrt(pow(x4[ip]-body_intercept_p_x[ip],2) + pow(y4[ip]-body_intercept_p_y[ip],2));
	if (s < min)
	{
		min = s;
		close_index = index4;
	}

	//setup for neuman BC
	double distance, distance2, min2, matDi, matDj, matD2i, matD2j, matDBIi, matDBIj;
	int bodyindex, bodyindex2;

	//if the node is inside the body, move it to the surface then set it to be a neuman condition
	//point 1
	if (ghostTagsP[index1] != -1 || index1 == close_index)
	{
		//setup
		x1[ip] = body_intercept_p_x[index1];
		y1[ip] = body_intercept_p_y[index1];
		n_x = image_point_p_x[index1] - x1[ip];
		n_y = image_point_p_y[index1] - y1[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x1[ip],2)+pow(by[k]-y1[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x1[ip],2)+pow(by[k]-y1[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;


		if (index1 == close_index)
		{
			matDClose = ( matDBIi*n_x/nl + matDBIj*n_y/nl );
		}
		q1[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl );

		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[ip]+a12*y1[ip];
	}
	//point 2
	if (ghostTagsP[index2] != -1 || index2 == close_index)
	{
		x2[ip] = body_intercept_p_x[index2];
		y2[ip] = body_intercept_p_y[index2];
		n_x = image_point_p_x[index2] - x2[ip];
		n_y = image_point_p_y[index2] - y2[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x2[ip],2)+pow(by[k]-y2[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x2[ip],2)+pow(by[k]-y2[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		//if this is the close index, set matDclose
		if (index2 == close_index)
		{
			matDClose = ( matDBIi*n_x/nl + matDBIj*n_y/nl );
			dudt[ip] = uB[0];
			ududx[ip] = uB0[0];
			vdudy[ip] = 0;
			dvdt[ip] = 0;
			udvdx[ip] = 0;
		}
		q2[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl );

		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[ip]+a22*y2[ip];
	}
	//point 3
	if (ghostTagsP[index3] != -1 || index3 == close_index)
	{
		x3[ip] = body_intercept_p_x[index3];
		y3[ip] = body_intercept_p_y[index3];
		n_x = image_point_p_x[index3] - x3[ip];
		n_y = image_point_p_y[index3] - y3[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x3[ip],2)+pow(by[k]-y3[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x3[ip],2)+pow(by[k]-y3[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		//if this is the close index, set matDclose
		if (index3 == close_index)
		{
			matDClose = ( matDBIi*n_x/nl + matDBIj*n_y/nl );
		}
		q3[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl );

		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x3[ip]+a32*y3[ip];
	}
	//4
	if (ghostTagsP[index4] != -1 || index4 == close_index)
	{
		x4[ip] = body_intercept_p_x[index4];
		y4[ip] = body_intercept_p_y[index4];
		n_x = image_point_p_x[index4] - x4[ip];
		n_y = image_point_p_y[index4] - y4[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find two closest body nodes
		min = 1;
		min2 = 1;
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x4[ip],2)+pow(by[k]-y4[ip],2));
			if (distance<min)
			{
				min = distance;
				bodyindex = k;
			}
		}
		for (int k=0; k<totalPoints; k++)
		{
			distance = sqrt(pow(bx[k]-x4[ip],2)+pow(by[k]-y4[ip],2));
			distance2 = sqrt(pow(bx[bodyindex]-bx[k],2)+pow(by[bodyindex]-bx[k],2));
			if (distance<min2 && distance2>0)
			{
				min2 = distance;
				bodyindex2 = k;
			}
		}

		//calc Du/Dt at body nodes
		matDi = (uB[bodyindex]-uB0[bodyindex])/dt;
		matDj = (vB[bodyindex]-vB0[bodyindex])/dt;
		matD2i = (uB[bodyindex2]-uB0[bodyindex2])/dt;
		matD2j = (vB[bodyindex2]-vB0[bodyindex2])/dt;

		//interp to BI
		matDBIi = matDi + (matD2i-matDi)/(min+min2)*min;
		matDBIj = matDj + (matD2j-matDj)/(min+min2)*min;

		//if this is the close index, set matDclose
		if (index4 == close_index)
		{
			matDClose = ( matDBIi*n_x/nl + matDBIj*n_y/nl );
		}
		q4[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl );

		a41 = 0;
		a42 = n_x/nl;
		a43 = n_y/nl;
		a44 = a43*x4[ip]+a42*y4[ip];
	}
	//solve equation for bilinear interpolation of values to image point
	//http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html  //for solving a 4x4 matrix exactly
	//https://www.physicsforums.com/threads/is-normal-derivative-a-definition.706458/   //for dealing with the normal at the boundary (how to calculate a normal deriviative)
	//df/dn = grad(f) dot n
	//f = a0 + a1x + a2y + a3xy
	//df/dn = ((a1+a3y)i + (a2+a3x)j) dot ((n_x/nl) i+ (n_y/nl)j)
	//where n_x, n_y and nl are the normal vector lengths in the x, y and magnitude respectively
	//solve for a
	/*  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 *  |0  n_x/nl n_y/nl (n_y*x+n_x*y)/nl|   |  |    =   |q | replace one row with this depending on which node is the closes to the body intercept<-
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 */

	double
	detA = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43
		  +a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41
		  +a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42
		  +a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41
		  -a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42
		  -a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43
		  -a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41
		  -a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;

	/*	       B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 */
	double b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	double b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	double b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	double b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	double b21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
	double b22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
	double b23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
	double b24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
	double b31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
	double b32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
	double b33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
	double b34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
	double b41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
	double b42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
	double b43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
	double b44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;

	/* Solve A*a = q for a
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 * interpolate for a value using the newly formed function
	 * p= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;
	 */
	a0[ip] = b11/detA*q1[ip]  +  b12/detA*q2[ip]  +  b13/detA*q3[ip]  +  b14/detA*q4[ip];
	a1[ip] = b21/detA*q1[ip]  +  b22/detA*q2[ip]  +  b23/detA*q3[ip]  +  b24/detA*q4[ip];
	a2[ip] = b31/detA*q1[ip]  +  b32/detA*q2[ip]  +  b33/detA*q3[ip]  +  b34/detA*q4[ip];
	a3[ip] = b41/detA*q1[ip]  +  b42/detA*q2[ip]  +  b43/detA*q3[ip]  +  b44/detA*q4[ip];
	//pressure at the image point
	double image_point_pressure = a0[ip] + a1[ip]*image_point_p_x[ip]    + a2[ip]*image_point_p_y[ip]    + a3[ip] * image_point_p_y[ip]   *image_point_p_x[ip];
	body_intercept_p[ip]        = a0[ip] + a1[ip]*body_intercept_p_x[ip] + a2[ip]*body_intercept_p_y[ip] + a3[ip] * body_intercept_p_x[ip]*body_intercept_p_y[ip]; //used for force calc

	//extrapolate pressure to the ghost node
	pressure[ip] = image_point_pressure + sqrt(pow(image_point_p_x[ip]-xv[I],2)+pow(image_point_p_y[ip]-yu[J],2))*matDClose;
}
}
