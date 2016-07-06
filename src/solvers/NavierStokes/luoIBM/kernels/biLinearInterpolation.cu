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
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, int timeStep)//test
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
	if (hybridTagsP[ip]<=0) //return if we're not at an interpolation point
		return;

	double	n_x,
			n_y,
			nl,
			du_dt,
			u_du_dx,
			v_du_dy,
			dv_dt,
			u_dv_dx,
			v_dv_dy;
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
	double X1u,X2u,X3u,X4u,Y1u,Y2u,Y3u,Y4u,velTemp,lTemp;
	double X1v,X2v,X3v,X4v,Y1v,Y2v,Y3v,Y4v;
	int i1u, i2u, i3u, i4u, i1v, i2v, i3v, i4v;
	//check if any points are inside of the body, then calculate the neuman boundary condition for them
	//point 1
	if (hybridTagsP[index1] == ip)
	{
		//setup
		x1[ip] = body_intercept_p_x[ip];
		y1[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x1[ip]; //flag this seems questionable
		n_y = image_point_p_y[ip] - y1[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[ip])
			ii++;
		while (yu[jj] < body_intercept_p_y[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[ip])
			ii++;
		while (yv[jj] < body_intercept_p_y[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//calc time derivatives //flag this doesn't work for rotating bodies because it is only using body index 0
		du_dt = (uB[0] - uB0[0])/dt;
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[ip] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(y1[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(y1[ip]-Y2u)/(Y4u-Y2u);					//flag which order do the subtractions go in
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-x1[ip]);

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[ip] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u + (nx-1)] + (u[i4u + (nx-1)] - u[i3u + (nx-1)])*(x1[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(x1[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y1[ip]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[ip] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(y1[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(y1[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x1[ip]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[ip] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(x1[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(x1[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y1[ip]);

		q1[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt + u_dv_dx + v_dv_dy));

		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[ip]+a12*y1[ip];
	}
	//point 2
	if (hybridTagsP[index2] == ip)
	{
		x2[ip] = body_intercept_p_x[ip];
		y2[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x2[ip];
		n_y = image_point_p_y[ip] - y2[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[ip])
			ii++;
		while (yu[jj] < body_intercept_p_y[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[ip])
			ii++;
		while (yv[jj] < body_intercept_p_y[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);
		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[ip] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(y2[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(y2[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-x2[ip]); //flag were using - uB here because the slope is wrong, there should be many deriviatives that are not dudx that have the wrong slope as well, fix them

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[ip] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u+(nx-1)] + (u[i4u+(nx-1)] - u[i3u+(nx-1)])*(x2[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(x2[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y2[ip]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[ip] ) < (X2v-X1v)*0.75 )
		{
			velTemp =u[i1v-1] + (u[i3v-1] - u[i1v-1])*(y2[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(y2[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x2[ip]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[ip] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(x2[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(x2[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y2[ip]);

		q2[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));

		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[ip]+a22*y2[ip];
	}
	//point 3
	if (hybridTagsP[index3] == ip)
	{
		x3[ip] = body_intercept_p_x[ip];
		y3[ip] = body_intercept_p_y[ip];
		n_x = image_point_p_x[ip] - x3[ip];
		n_y = image_point_p_y[ip] - y3[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[ip])
			ii++;
		while (yu[jj] < body_intercept_p_y[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[ip])
			ii++;
		while (yv[jj] < body_intercept_p_y[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[ip] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(y3[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(y3[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-x3[ip]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y1u-body_intercept_p_y[ip] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(x3[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(x3[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y3[ip]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[ip] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(y3[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(y3[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x3[ip]);

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[ip] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(x3[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(x3[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y3[ip]);

		q3[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));

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
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[ip])
			ii++;
		while (yu[jj] < body_intercept_p_y[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[ip])
			ii++;
		while (yv[jj] < body_intercept_p_y[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);
		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[ip] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(y4[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(y4[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-x4[ip]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[ip] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(x4[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(x4[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y4[ip]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[ip] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i1v-1] + (u[i3v-1] - u[i1v-1])*(y4[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(y4[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x4[ip]);
		//if (timeStep == 1)
		//	q4[ip] = u_dv_dx;

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[ip] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(x4[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(x4[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y4[ip]);

		q4[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));

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

	dudt[ip] = du_dt;
	ududx[ip] = u_du_dx;
	vdudy[ip] = v_du_dy;
	dvdt[ip] = dv_dt;
	udvdx[ip] = u_dv_dx;
	vdvdy[ip] = v_dv_dy;
	pressureStar[ip] = a0[ip] + a1[ip]*xv[I] + a2[ip]*yu[J] + a3[ip]*xv[I]*yu[J];
}

//flag this function is a mess
__global__
void interpolatePressureToGhostNode(double *pressure, double *u, int *ghostTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y, double *body_intercept_p,
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
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
			du_dt,
			u_du_dx,
			v_du_dy,
			dv_dt,
			u_dv_dx,
			v_dv_dy;
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
	double X1u,X2u,X3u,X4u,Y1u,Y2u,Y3u,Y4u,velTemp,lTemp;
	double X1v,X2v,X3v,X4v,Y1v,Y2v,Y3v,Y4v;
	int i1u, i2u, i3u, i4u, i1v, i2v, i3v, i4v;

	//if the node is inside the body, set it to be a neuman condition
	//point 1
	if (ghostTagsP[index1] != -1)
	{
		//the material derivatve is calculated as follows
		//-(n_x/nl * (du/dt + u*du/dx + v*du/dy) + n_y/nl * (dv/dt + u*dv/dx + v*dv/dy));
		//it is discritized around the body intercept(as opposed to the image point or something)
		//the derivatvies are found as follows:
		//1 find the four u and v nodes that surround the BI
		//2 interpolate between two of those nodes to get a u or v value at an x or y location that is horizontal or vertical to the BI
		//3 use interpolated value and body velocity to calc derivative

		//setup
		x1[ip] = body_intercept_p_x[index1];
		y1[ip] = body_intercept_p_y[index1];
		n_x = image_point_p_x[index1] - x1[ip];
		n_y = image_point_p_y[index1] - y1[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < x1[ip])
			ii++;
		while (yu[jj] < y1[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < x1[ip])
			ii++;
		while (yv[jj] < y1[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//calc time derivatives //flag this doesn't work for rotating bodies because it is only using body index 0
		du_dt = (uB[0] - uB0[0])/dt;
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[index1] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(y1[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(y1[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-x1[ip]);

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index1] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u + (nx-1)] + (u[i4u + (nx-1)] - u[i3u + (nx-1)])*(x1[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(x1[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y1[ip]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[index1] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(y1[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(y1[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x1[ip]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[index1] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(x1[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(x1[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y1[ip]);

		q1[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt + u_dv_dx + v_dv_dy));
		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[ip]+a12*y1[ip];
	}
	//point 2
	if (ghostTagsP[index2] != -1)
	{
		x2[ip] = body_intercept_p_x[index2];
		y2[ip] = body_intercept_p_y[index2];
		n_x = image_point_p_x[index2] - x2[ip];
		n_y = image_point_p_y[index2] - y2[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < x2[ip])
			ii++;
		while (yu[jj] < y2[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < x2[ip])
			ii++;
		while (yv[jj] < y2[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);
		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[index2] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(y2[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(y2[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-x2[ip]);

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index2] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u+(nx-1)] + (u[i4u+(nx-1)] - u[i3u+(nx-1)])*(x2[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(x2[ip]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y2[ip]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[index2] ) < (X2v-X1v)*0.75 )
		{
			velTemp =u[i1v-1] + (u[i3v-1] - u[i1v-1])*(y2[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(y2[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x2[ip]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[index2] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(x2[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(x2[ip]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y2[ip]);

		q2[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[ip]+a22*y2[ip];
	}
	//point 3
	if (ghostTagsP[index3] != -1)
	{
		x3[ip] = body_intercept_p_x[index3];
		y3[ip] = body_intercept_p_y[index3];
		n_x = image_point_p_x[index3] - x3[ip];
		n_y = image_point_p_y[index3] - y3[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < x3[ip])
			ii++;
		while (yu[jj] < y3[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < x3[ip])
			ii++;
		while (yv[jj] < y3[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[index3] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(y3[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(y3[ip]-Y2u)/(Y4u-Y2u);
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-x3[ip]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y1u-body_intercept_p_y[index3] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(x3[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(x3[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y3[ip]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[index3] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(y3[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(y3[ip]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x3[ip]);

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[index3] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(x3[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(x3[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y3[ip]);

		q3[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x3[ip]+a32*y3[ip];
	}
	//4
	if (ghostTagsP[index4] != -1)
	{
		x4[ip] = body_intercept_p_x[index4];
		y4[ip] = body_intercept_p_y[index4];
		n_x = image_point_p_x[index4] - x4[ip];
		n_y = image_point_p_y[index4] - y4[ip];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < x4[ip])
			ii++;
		while (yu[jj] < y4[ip])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < x4[ip])
			ii++;
		while (yv[jj] < y4[ip])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);
		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[index4] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(y4[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(y4[ip]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-x4[ip]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index4] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(x4[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(x4[ip]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-y4[ip]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[index4] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i1v-1] + (u[i3v-1] - u[i1v-1])*(y4[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(y4[ip]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-x4[ip]);

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[index4] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(x4[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(x4[ip]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-y4[ip]);

		q4[ip] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
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
	//interpolate pressure to the ghost node
	double matD = 0;
	if (close_index == index1)
	{
		n_x = image_point_p_x[index1] - body_intercept_p_x[index1];
		n_y = image_point_p_y[index1] - body_intercept_p_y[index1];
		nl = sqrt(n_x*n_x+n_y*n_y);

		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[index1])
			ii++;
		while (yu[jj] < body_intercept_p_y[index1])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[index1])
			ii++;
		while (yv[jj] < body_intercept_p_y[index1])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//calc time derivatives //flag this doesn't work for rotating bodies because it is only using body index 0
		du_dt = (uB[0] - uB0[0])/dt;
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[index1] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(body_intercept_p_y[index1]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(body_intercept_p_y[index1]-Y2u)/(Y4u-Y2u);
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_x[index1]);

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index1] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u + (nx-1)] + (u[i4u + (nx-1)] - u[i3u + (nx-1)])*(body_intercept_p_x[index1]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(body_intercept_p_x[index1]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_y[index1]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[index1] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(body_intercept_p_y[index1]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(body_intercept_p_y[index1]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-body_intercept_p_x[index1]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[index1] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(body_intercept_p_x[index1]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(body_intercept_p_x[index1]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-body_intercept_p_y[index1]);

		matD = (n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt + u_dv_dx + v_dv_dy));
	}
	else if (close_index == index2)
	{
		n_x = image_point_p_x[index2] - body_intercept_p_x[index2];
		n_y = image_point_p_y[index2] - body_intercept_p_y[index2];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[index2])
			ii++;
		while (yu[jj] < body_intercept_p_y[index2])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[index2])
			ii++;
		while (yv[jj] < body_intercept_p_y[index2])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[index2] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(body_intercept_p_y[index2]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(body_intercept_p_y[index2]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_x[index2]);

		//find du/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index2] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i3u+(nx-1)] + (u[i4u+(nx-1)] - u[i3u+(nx-1)])*(body_intercept_p_x[index2]-X3u)/(X4u-X3u);
			lTemp = Y3u + Y3u - Y1u;
		}
		else
		{
			velTemp = u[i3u] + (u[i4u] - u[i3u])*(body_intercept_p_x[index2]-X3u)/(X4u-X3u);
			lTemp = Y3u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_y[index2]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[index2] ) < (X2v-X1v)*0.75 )
		{
			velTemp =u[i1v-1] + (u[i3v-1] - u[i1v-1])*(body_intercept_p_y[index2]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(body_intercept_p_y[index2]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-body_intercept_p_x[index2]);

		//find dv/dy
		//U_3 + (U_4-U_3)*(XBI-X3)/(X4-X3)
		if ( abs( Y3v-body_intercept_p_y[index2] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i3v+nx] + (u[i4v+nx] - u[i3v+nx])*(body_intercept_p_x[index2]-X3v)/(X4v-X3v);
			lTemp = Y3v+Y3v-Y1v;
		}
		else
		{
			velTemp = u[i3v] + (u[i4v] - u[i3v])*(body_intercept_p_x[index2]-X3v)/(X4v-X3v);
			lTemp = Y3v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-body_intercept_p_y[index2]);

		matD = (n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
	}
	else if (close_index == index3)
	{
		n_x = image_point_p_x[index3] - body_intercept_p_x[index3];
		n_y = image_point_p_y[index3] - body_intercept_p_y[index3];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[index3])
			ii++;
		while (yu[jj] < body_intercept_p_y[index3])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[index3])
			ii++;
		while (yv[jj] < body_intercept_p_y[index3])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_2 + (U_4-U_2)*(YBI-Y2)/(Y4-Y2)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X2u-body_intercept_p_x[index3] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i2u+1] + (u[i4u+1] - u[i2u+1])*(body_intercept_p_y[index3]-Y2u)/(Y4u-Y2u);
			lTemp = X2u + X2u - X1u;
		}
		else
		{
			velTemp = u[i2u] + (u[i4u] - u[i2u])*(body_intercept_p_y[index3]-Y2u)/(Y4u-Y2u);
			lTemp = X2u;
		}
		u_du_dx = uB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_x[index3]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y1u-body_intercept_p_y[index3] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(body_intercept_p_x[index3]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(body_intercept_p_x[index3]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_y[index3]);

		//find dv/dx
		//V_2 + (V_4-V_2)(YBI-Y2)/(Y4-Y2)
		//check if were too close to the v node in the x direction
		if ( abs( X2v-body_intercept_p_x[index3] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i2v+1] + (u[i4v+1] - u[i2v+1])*(body_intercept_p_y[index3]-Y2v)/(Y4v-Y2v);
			lTemp = X2v+X2v-X1v;
		}
		else
		{
			velTemp = u[i2v] + (u[i4v] - u[i2v])*(body_intercept_p_y[index3]-Y2v)/(Y4v-Y2v);
			lTemp = X2v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-body_intercept_p_x[index3]);

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[index3] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(body_intercept_p_x[index3]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(body_intercept_p_x[index3]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-body_intercept_p_y[index3]);

		matD = (n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
	}
	else if (close_index == index4)
	{
		n_x = image_point_p_x[index4] - body_intercept_p_x[index4];
		n_y = image_point_p_y[index4] - body_intercept_p_y[index4];
		nl = sqrt(n_x*n_x+n_y*n_y);
		//find the four u velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xu[ii] < body_intercept_p_x[index4])
			ii++;
		while (yu[jj] < body_intercept_p_y[index4])
			jj++;
		X3u = xu[ii-1]; X4u = xu[ii];
		X1u = xu[ii-1]; X2u = xu[ii];

		Y3u = yu[jj];   Y4u = yu[jj];
		Y1u = yu[jj-1]; Y2u = yu[jj-1];

		i3u = jj*(nx-1) + ii -1;	i4u = jj*(nx-1) + ii;
		i1u = (jj-1)*(nx-1)+ii-1;	i2u = (jj-1)*(nx-1)+ii;

		//find the four v velocity nodes around the body intercept
		ii = I-5; jj = J-5;
		while (xv[ii] < body_intercept_p_x[index4])
			ii++;
		while (yv[jj] < body_intercept_p_y[index4])
			jj++;
		X3v = xv[ii-1];	X4v = xv[ii];
		X1v = xv[ii-1];	X2v = xv[ii];

		Y3v = yv[jj];	Y4v = Y3v;
		Y1v = yv[jj-1];	Y2v = Y1v;

		i3v = jj*nx+ii-1 + ny*(nx-1);		i4v = jj*nx+ii + ny*(nx-1);
		i1v = (jj-1)*(nx)+ii-1 + ny*(nx-1); i2v = (jj-1)*nx+ii + ny*(nx-1);

		//time derivatives
		du_dt = (uB[0] - uB0[0])/dt;//flag this doesn't work for rotating bodies because it is only using body index 0
		dv_dt = (vB[0] - vB0[0])/dt;

		//find du/dx
		//U_1 + (U_3-U_1)*(YBI-Y1)/(Y3-Y1)
		//check if were too close to the u node in x direction to get good values, if we are: interpolate from u nodes farther away
		if ( abs( X1u-body_intercept_p_x[index4] ) < (X2u-X1u)*0.75 )
		{
			velTemp = u[i1u-1] + (u[i3u-1] - u[i1u-1])*(body_intercept_p_y[index4]-Y1u)/(Y3u-Y1u);
			lTemp = X1u + X1u - X2u;
		}
		else
		{
			velTemp = u[i1u] + (u[i3u] - u[i1u])*(body_intercept_p_y[index4]-Y1u)/(Y3u-Y1u);
			lTemp = X1u;
		}
		u_du_dx = -uB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_x[index4]);

		//find du/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		//check if were too close to u node in y direction
		if ( abs( Y3u-body_intercept_p_y[index4] ) < (Y3u-Y1u)*0.75 )
		{
			velTemp = u[i1u-(nx-1)] + (u[i2u-(nx-1)] - u[i2u-(nx-1)])*(body_intercept_p_x[index4]-X1u)/(X2u-X1u);
			lTemp = Y1u + Y1u - Y3u;
		}
		else
		{
			velTemp = u[i1u] + (u[i2u] - u[i2u])*(body_intercept_p_x[index4]-X1u)/(X2u-X1u);
			lTemp = Y1u;
		}
		v_du_dy = vB[0]  *  (velTemp - uB[0])/(lTemp-body_intercept_p_y[index4]);

		//find dv/dx
		//V_1 + (V_3-V_1)(YBI-Y1)/(Y3-Y1)
		//check if were too close to the v node in the x direction
		if ( abs( X1v-body_intercept_p_x[index4] ) < (X2v-X1v)*0.75 )
		{
			velTemp = u[i1v-1] + (u[i3v-1] - u[i1v-1])*(body_intercept_p_y[index4]-Y1v)/(Y3v-Y1v);
			lTemp = X1v+X1v-X2v;
		}
		else
		{
			velTemp = u[i1v] + (u[i3v] - u[i1v])*(body_intercept_p_y[index4]-Y1v)/(Y3v-Y1v);
			lTemp = X1v;
		}
		u_dv_dx = uB[0]  *  (velTemp-vB[0])/(lTemp-body_intercept_p_x[index4]);

		//find dv/dy
		//U_1 + (U_2-U_1)*(XBI-X1)/(X2-X1)
		if ( abs( Y1v-body_intercept_p_y[index4] ) < (Y3v-Y1v)*0.75 )
		{
			velTemp = u[i1v-nx] + (u[i2v-nx] - u[i1v-nx])*(body_intercept_p_x[index4]-X1v)/(X2v-X1v);
			lTemp = Y1v+Y1v-Y3v;
		}
		else
		{
			velTemp = u[i1v] + (u[i2v] - u[i1v])*(body_intercept_p_x[index4]-X1v)/(X2v-X1v);
			lTemp = Y1v;
		}
		v_dv_dy = vB[0]  *  (velTemp - vB[0])/(lTemp-body_intercept_p_y[index4]);

		matD = (n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
	}
	dudt[ip] = du_dt;
	ududx[ip] = u_du_dx;
	vdudy[ip] = v_du_dy;
	dvdt[ip] = dv_dt;
	udvdx[ip] = u_dv_dx;
	vdvdy[ip] = v_dv_dy;
	pressure[ip] = image_point_pressure + sqrt(pow(image_point_p_x[ip]-xv[I],2)+pow(image_point_p_y[ip]-yu[J],2))*matD;
}

}
