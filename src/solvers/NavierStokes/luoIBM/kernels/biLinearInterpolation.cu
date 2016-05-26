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
void interpolateVelocityToGhostNodeX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,//flag vB needed?
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u)//testing variables
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
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
	 ip_u[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
	 u[iu] =  2*uB[0] - ip_u[iu]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToGhostNodeY(double *u, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u)//testing variables
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
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
	 ip_u[iv] = a0 + a1*image_point_x[iv] + a2*image_point_y[iv] + a3*image_point_x[iv]*image_point_y[iv];
	 u[iv] =  2*vB[0] - ip_u[iv]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4)//test
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
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
	x1[iu] = ii;//xu[ii-1];
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
	/*if (hybridTagsUV[(jj-1)*(nx-1)+ii-1] == iu)
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
	}*/

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
	 ustar[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
}







}
