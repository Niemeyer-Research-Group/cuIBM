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
void interpolateVelocityToGhostNodeX(double *u, bool set, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int *i_start, int *j_start, int width, int nx, int ny,
							int *index1, int *index2, int *index3, int *index4,
							double *q1coef, double *q2coef, double *q3coef, double *q4coef,
							double *x1, double *x2, double *x3, double *x4,
							double *y1, double *y2, double *y3, double *y4,
							double *q1, double *q2, double *q3, double *q4, double *image_point_u)//testing variables
{//In the luo et al method they only move corners coincident to the GN to the boundary. We are moving all corners inside to the boundary
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
	while (yu[jj] <image_point_y[iu])
		jj++;
	double x[4] = {xu[ii-1], xu[ii], xu[ii-1], xu[ii]};
	double y[4] = {yu[jj-1], yu[jj-1], yu[jj], yu[jj]};
	//find index at corners and the u value at the corners
	int index[4] = {(jj-1)*(nx-1)+ii-1,   (jj-1)*(nx-1)+ii,   jj*(nx-1)+ii-1,   jj*(nx-1)+ii};
	double q[4] = {u[index[0]], u[index[1]], u[index[2]], u[index[3]]};

	//find the closest corner to the body intercept
	double min = 1.0;
	double s;
	int close_index;
	bool inflag = false; //a boolean that is true if there is a node inside the body
	for (int l=0;l<4;l++)
	{
		//find the closest node to the BI
		s = sqrt(pow(x[l]-body_intercept_x[iu],2) + pow(y[l]-body_intercept_y[iu],2));
		if (s<min)
		{
			min = s;
			close_index = index[l];
		}
		//check if any of the points are inside the body
		if (ghostTagsUV[index[l]]>0)
				inflag = true;
	}

	//if point is inside of the body
	//or if no points are inside the body and the node is the closest to the BI
	//	then move them to the body intercept
	//point 1
	for (int l=0;l<4;l++)
	{
	if ( ghostTagsUV[index[l]] > 0)
	{
		x[l] = body_intercept_x[index[l]];
		y[l] = body_intercept_y[index[l]];
		q[l] = uB[0];
	}
	/*else if ( index[l]==close_index && !inflag ) //uncomment this if you want to move the closest node outside of the body to the body
	{
		x[l] = body_intercept_x[iu];
		y[l] = body_intercept_y[iu];
		q[l] = uB[0];
	}*/
	}

	x1[iu] = x[0];
	x2[iu] = x[1];
	x3[iu] = x[2];
	x4[iu] = x[3];
	y1[iu] = y[0];
	y2[iu] = y[1];
	y3[iu] = y[2];
	y4[iu] = y[3];
	q1[iu] = q[0];
	q2[iu] = q[1];
	q3[iu] = q[2];
	q4[iu] = q[3];
	index1[iu] = index[0];
	index2[iu] = index[1];
	index3[iu] = index[2];
	index4[iu] = index[3];

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
	 q1coef[iu] = (b11+b21*image_point_x[iu]+b31*image_point_y[iu]+b41*image_point_x[iu]*image_point_y[iu])/detA;
	 q2coef[iu] = (b12+b22*image_point_x[iu]+b32*image_point_y[iu]+b42*image_point_x[iu]*image_point_y[iu])/detA;
	 q3coef[iu] = (b13+b23*image_point_x[iu]+b33*image_point_y[iu]+b43*image_point_x[iu]*image_point_y[iu])/detA;
	 q4coef[iu] = (b14+b24*image_point_x[iu]+b34*image_point_y[iu]+b44*image_point_x[iu]*image_point_y[iu])/detA;
	 image_point_u[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
	 if (set)
		 u[iu] =  2*uB[0] - image_point_u[iu]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToGhostNodeY(double *u, bool set, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int *i_start, int *j_start, int width, int nx, int ny,
							int *index1, int *index2, int *index3, int *index4,
							double *q1coef, double *q2coef, double *q3coef, double *q4coef,
							double *x1, double *x2, double *x3, double *x4,
							double *y1, double *y2, double *y3, double *y4,
							double *q1, double *q2, double *q3, double *q4, double *image_point_u)//testing variables
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
	while (yv[jj] <image_point_y[iv])
		jj++;
	double x[4] = {xv[ii-1], xv[ii], xv[ii-1], xv[ii]};
	double y[4] = {yv[jj-1], yv[jj-1], yv[jj], yv[jj]};
	//find index at corners and the u value at the corners
	int index[4] = {(jj-1)*nx+ii-1 + (nx-1)*ny,   (jj-1)*nx+ii + (nx-1)*ny,   jj*nx+ii-1 + (nx-1)*ny,   jj*nx+ii + (nx-1)*ny};
	double q[4] = {u[index[0]], u[index[1]], u[index[2]], u[index[3]]};

	//find the closest corner to the body intercept
	double min = 1.0;
	double s;
	int close_index;
	bool inflag = false; //a boolean that is true if there is a node inside the body
	for (int l=0;l<4;l++)
	{
		//find the closest node to the BI
		s = sqrt(pow(x[l]-body_intercept_x[iv],2) + pow(y[l]-body_intercept_y[iv],2));
		if (s<min)
		{
			min = s;
			close_index = index[l];
		}
		//check if any of the points are inside the body
		if (ghostTagsUV[index[l]]>0)
				inflag = true;
	}

	//if point is inside of the body
	//or if no points are inside the body and the node is the closest to the BI
	//	then move them to the body intercept
	//point 1
	for (int l=0;l<4;l++)
	{
	if ( ghostTagsUV[index[l]] > 0)
	{
		x[l] = body_intercept_x[index[l]];
		y[l] = body_intercept_y[index[l]];
		q[l] = vB[0];
	}
	/*else if ( index[l]==close_index && !inflag ) //uncomment this if you want to move the closest node outside of the body to the body
	{
		x[l] = body_intercept_x[iv];
		y[l] = body_intercept_y[iv];
		q[l] = vB[0];
	}*/
	}

	x1[iv] = x[0];
	x2[iv] = x[1];
	x3[iv] = x[2];
	x4[iv] = x[3];
	y1[iv] = y[0];
	y2[iv] = y[1];
	y3[iv] = y[2];
	y4[iv] = y[3];
	q1[iv] = q[0];
	q2[iv] = q[1];
	q3[iv] = q[2];
	q4[iv] = q[3];
	index1[iv] = index[0];
	index2[iv] = index[1];
	index3[iv] = index[2];
	index4[iv] = index[3];

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
	 q1coef[iv] = (b11+b21*image_point_x[iv]+b31*image_point_y[iv]+b41*image_point_x[iv]*image_point_y[iv])/detA;
	 q2coef[iv] = (b12+b22*image_point_x[iv]+b32*image_point_y[iv]+b42*image_point_x[iv]*image_point_y[iv])/detA;
	 q3coef[iv] = (b13+b23*image_point_x[iv]+b33*image_point_y[iv]+b43*image_point_x[iv]*image_point_y[iv])/detA;
	 q4coef[iv] = (b14+b24*image_point_x[iv]+b34*image_point_y[iv]+b44*image_point_x[iv]*image_point_y[iv])/detA;
	 image_point_u[iv] = a0 + a1*image_point_x[iv] + a2*image_point_y[iv] + a3*image_point_x[iv]*image_point_y[iv];
	 if (set)
		 u[iv] =  2*vB[0] - image_point_u[iv]; //u_gn = 2*u_BI  - u_IP //flag doesn't currently work with a rotating body because of uB[0], need to use the actual u at the body intercept
}

__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV,
										double *bx, double *by, double *uB, double *yu, double *xu,
										double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
										int *i_start, int *j_start, int width, int nx, int ny,
										int *index1, int *index2, int *index3, int *index4,
										double *q1coef, double *q2coef, double *q3coef, double *q4coef,
										double *x1, double *x2, double *x3, double *x4,
										double *y1, double *y2, double *y3, double *y4,
										double *q1, double *q2, double *q3, double *q4, double *image_point_u)
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

	/*   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *   	|					   |
	 *   	(x1,y1)__________(x2,y2)
	 *
	 *   *(BI_x,BI_y)
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

	index1[iu] = (jj-1)*(nx-1)+ii-1;
	index2[iu] = (jj-1)*(nx-1)+ii;
	index3[iu] = jj*(nx-1)+ii-1;
	index4[iu] = jj*(nx-1)+ii;

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

	double detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
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
	 double a0 = (b11*q1[iu]  +  b12*q2[iu]  +  b13*q3[iu]  +  b14*q4[iu])/detA;
	 double a1 = (b21*q1[iu]  +  b22*q2[iu]  +  b23*q3[iu]  +  b24*q4[iu])/detA;
	 double a2 = (b31*q1[iu]  +  b32*q2[iu]  +  b33*q3[iu]  +  b34*q4[iu])/detA;
	 double a3 = (b41*q1[iu]  +  b42*q2[iu]  +  b43*q3[iu]  +  b44*q4[iu])/detA;
	 q1coef[iu] = (b11+b21*xu[I]+b31*yu[J]+b41*xu[I]*yu[J])/detA;
	 q2coef[iu] = (b12+b22*xu[I]+b32*yu[J]+b42*xu[I]*yu[J])/detA;
	 q3coef[iu] = (b13+b23*xu[I]+b33*yu[J]+b43*xu[I]*yu[J])/detA;
	 q4coef[iu] = (b14+b24*xu[I]+b34*yu[J]+b44*xu[I]*yu[J])/detA;

	 ustar[iu] = a0 + a1*xu[I] + a2*yu[J] + a3*yu[J]*xu[I];
	 image_point_u[iu] = a0 + a1*image_point_x[iu] + a2*image_point_y[iu] + a3*image_point_x[iu]*image_point_y[iu];
}

__global__
void interpolateVelocityToHybridNodeY(double *u, double *ustar, int *hybridTagsUV,
										double *bx, double *by, double *vB, double *yv, double *xv,
										double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
										int *i_start, int *j_start, int width, int nx, int ny,
										int *index1, int *index2, int *index3, int *index4,
										double *q1coef, double *q2coef, double *q3coef, double *q4coef,
										double *x1, double *x2, double *x3, double *x4,
										double *y1, double *y2, double *y3, double *y4,
										double *q1, double *q2, double *q3, double *q4, double *image_point_u)
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
	index1[iv] = (jj-1)*nx+ii-1 + (nx-1)*ny;
	index2[iv] = (jj-1)*nx+ii + (nx-1)*ny;
	index3[iv] = jj*nx+ii-1 + (nx-1)*ny;
	index4[iv] = jj*nx+ii + (nx-1)*ny;

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

	double detA = 1*a22*a33*a44 + 1*a23*a34*a42 + 1*a24*a32*a43
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
	 double a0 = (b11*q1[iv]  +  b12*q2[iv]  +  b13*q3[iv]  +  b14*q4[iv])/detA;
	 double a1 = (b21*q1[iv]  +  b22*q2[iv]  +  b23*q3[iv]  +  b24*q4[iv])/detA;
	 double a2 = (b31*q1[iv]  +  b32*q2[iv]  +  b33*q3[iv]  +  b34*q4[iv])/detA;
	 double a3 = (b41*q1[iv]  +  b42*q2[iv]  +  b43*q3[iv]  +  b44*q4[iv])/detA;
	 q1coef[iv] = (b11+b21*xv[I]+b31*yv[J]+b41*xv[I]*yv[J])/detA;
	 q2coef[iv] = (b12+b22*xv[I]+b32*yv[J]+b42*xv[I]*yv[J])/detA;
	 q3coef[iv] = (b13+b23*xv[I]+b33*yv[J]+b43*xv[I]*yv[J])/detA;
	 q4coef[iv] = (b14+b24*xv[I]+b34*yv[J]+b44*xv[I]*yv[J])/detA;

	 ustar[iv] = a0 + a1*xv[I] + a2*yv[J] + a3*yv[J]*xv[I];
	 image_point_u[iv] = a0 + a1*image_point_x[iv] + a2*image_point_y[iv] + a3*image_point_x[iv]*image_point_y[iv];
}

__global__
void interpolatePressureToHybridNode(double *pressure, double *pressureStar, double *u, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv, //yv xu not used?
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
									int *index1, int *index2, int *index3, int *index4,
									double *q1coef, double *q2coef, double *q3coef, double *q4coef,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4,
									double *y1, double *y2, double *y3, double *y4,
									double *q1, double *q2, double *q3, double *q4)
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
	double x[4];
	double y[4];
	double q[4];
	int index[4];
	while (xv[ii] < image_point_p_x[ip])
		ii++;
	x[0] = xv[ii-1];
	x[1] = xv[ii];
	x[2] = xv[ii-1];
	x[3] = xv[ii];
	while (yu[jj] <image_point_p_y[ip])
		jj++;
	y[0] = yu[jj-1];
	y[1] = yu[jj-1];
	y[2] = yu[jj];
	y[3] = yu[jj];

	//find q1,q2,q3,q4
	index[0] = (jj-1)*nx+ii-1;
	index[1] = (jj-1)*nx+ii;
	index[2] = jj*nx+ii-1;
	index[3] = jj*nx+ii;

	for (int l=0; l<4; l++)
		q[l] = pressure[index[l]];

	double a[16] = {1, x[0], y[0], x[0]*y[0],
					1, x[1], y[1], x[1]*y[1],
					1, x[2], y[2], x[2]*y[2],
					1, x[3], y[3], x[3]*y[3]};

	//setup for neuman BC


	double dudt;
	double dvdt;
	for (int l = 0; l<4; l++)
	{
		//move the closes node to the body to the surface then calculate the neuman boundary condition for it
		if ( hybridTagsP[index[l]] == ip )
		{
			dudt = (uB[0] - uB0[0])/dt;
			dvdt = (vB[0] - vB0[0])/dt;
			x[l] = body_intercept_p_x[ip];
			y[l] = body_intercept_p_y[ip];
			n_x = image_point_p_x[ip] - x[l];
			n_y = image_point_p_y[ip] - y[l];
			nl = sqrt(n_x*n_x + n_y*n_y);

			q[l] = - ( n_x / nl * dudt + n_y/nl * dvdt);
			a[l*4] = 0;
			a[l*4+1] = n_x/nl;
			a[l*4+2] = n_y/nl;
			a[l*4+3] = n_y/nl*x[l] + n_y/nl*y[l];
		}
	}

	x1[ip] = x[0];
	x2[ip] = x[1];
	x3[ip] = x[2];
	x4[ip] = x[3];
	y1[ip] = y[0];
	y2[ip] = y[1];
	y3[ip] = y[2];
	y4[ip] = y[3];
	q1[ip] = q[0];
	q2[ip] = q[1];
	q3[ip] = q[2];
	q4[ip] = q[3];
	index1[ip] = index[0];
	index2[ip] = index[1];
	index3[ip] = index[2];
	index4[ip] = index[3];
	double a11 = a[0], a12 = a[1],  a13 = a[2], a14 = a[3];
	double a21 = a[4], a22 = a[5],  a23 = a[6], a24 = a[7];
	double a31 = a[8], a32 = a[9],  a33 = a[10], a34 = a[11];
	double a41 = a[12], a42 = a[13],  a43 = a[14], a44 = a[15];

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

	 q1coef[ip] = (b11 + b21*xv[I] + b31*yu[J] + b41*xv[I]*yu[J])/detA;
	 q2coef[ip] = (b12 + b22*xv[I] + b32*yu[J] + b42*xv[I]*yu[J])/detA;
	 q3coef[ip] = (b13 + b23*xv[I] + b33*yu[J] + b43*xv[I]*yu[J])/detA;
	 q4coef[ip] = (b14 + b24*xv[I] + b34*yu[J] + b44*xv[I]*yu[J])/detA;

	pressureStar[ip] = a0[ip] + a1[ip]*xv[I] + a2[ip]*yu[J] + a3[ip]*xv[I]*yu[J];
}

__global__
void interpolatePressureToGhostNode(double *pressure, bool set, double *u, int *ghostTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y, double *body_intercept_p,
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
									int *index1, int *index2, int *index3, int *index4,
									double *q1coef, double *q2coef, double *q3coef, double *q4coef,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4,
									double *y1, double *y2, double *y3, double *y4,
									double *q1, double *q2, double *q3, double *q4)//test
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
	int close_index;
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
	while (yu[jj] < image_point_p_y[ip])
		jj++;
	double x[4] = {xv[ii-1], xv[ii], xv[ii-1], xv[ii]};
	double y[4] = {yu[jj-1], yu[jj-1], yu[jj], yu[jj]};
	//find index at corners and the u value at the corners
	int index[4] = {(jj-1)*nx+ii-1,   (jj-1)*nx+ii,   jj*nx+ii-1,   jj*nx+ii};
	double q[4] = {pressure[index[0]], pressure[index[1]], pressure[index[2]], pressure[index[3]]};

	double a[16] = {1, x[0], y[0], x[0]*y[0],
					1, x[1], y[1], x[1]*y[1],
					1, x[2], y[2], x[2]*y[2],
					1, x[3], y[3], x[3]*y[3]};

	//find the closest corner to the body intercept
	double min = 1.0;
	double s;
	for (int l=0;l<4;l++)
	{
		//find the closest node to the BI
		s = sqrt(pow(x[l]-body_intercept_p_x[ip],2) + pow(y[l]-body_intercept_p_y[ip],2));
		if (s<min)
		{
			min = s;
			close_index = index[l];
		}
	}

	//setup for neuman BC
	double dudt;
	double dvdt;
	for (int l=0; l<4; l++)
	{
		//set nodes inside the body to neuman bc
		if ( ghostTagsP[index[l]] > 0 )
		{
			dudt = (uB[0] - uB0[0])/dt;
			dvdt = (vB[0] - vB0[0])/dt;
			x[l] = body_intercept_p_x[index[l]];
			y[l] = body_intercept_p_y[index[l]];
			n_x = image_point_p_x[index[l]] - x[l];
			n_y = image_point_p_y[index[l]] - y[l];
			nl = sqrt(n_x*n_x + n_y*n_y);

			q[l] = - ( n_x / nl * dudt + n_y/nl * dvdt);
			a[l*4] = 0;
			a[l*4+1] = n_x/nl;
			a[l*4+2] = n_y/nl;
			a[l*4+3] = n_y/nl*x[l] + n_y/nl*y[l];
		}
		//if the node is the closest to the body, set the closeMatD
		if (index[l] == close_index)
		{
			dudt = (uB[0] - uB0[0])/dt;
			dvdt = (vB[0] - vB0[0])/dt;
			x[l] = body_intercept_p_x[ip];
			y[l] = body_intercept_p_y[ip];
			n_x = image_point_p_x[ip] - x[l];
			n_y = image_point_p_y[ip] - y[l];
			nl = sqrt(n_x*n_x + n_y*n_y);

			matDClose = ( n_x / nl * dudt + n_y/nl * dvdt);
		}
	}
	x1[ip] = x[0];
	x2[ip] = x[1];
	x3[ip] = x[2];
	x4[ip] = x[3];
	y1[ip] = y[0];
	y2[ip] = y[1];
	y3[ip] = y[2];
	y4[ip] = y[3];
	q1[ip] = q[0];
	q2[ip] = q[1];
	q3[ip] = q[2];
	q4[ip] = q[3];
	index1[ip] = index[0];
	index2[ip] = index[1];
	index3[ip] = index[2];
	index4[ip] = index[3];
	double a11 = a[0],  a12 = a[1],   a13 = a[2],  a14 = a[3];
	double a21 = a[4],  a22 = a[5],   a23 = a[6],  a24 = a[7];
	double a31 = a[8],  a32 = a[9],   a33 = a[10], a34 = a[11];
	double a41 = a[12], a42 = a[13],  a43 = a[14], a44 = a[15];
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
	a0[ip] = b11/detA*q[0]  +  b12/detA*q[1]  +  b13/detA*q[2]  +  b14/detA*q[3];
	a1[ip] = b21/detA*q[0]  +  b22/detA*q[1]  +  b23/detA*q[2]  +  b24/detA*q[3];
	a2[ip] = b31/detA*q[0]  +  b32/detA*q[1]  +  b33/detA*q[2]  +  b34/detA*q[3];
	a3[ip] = b41/detA*q[0]  +  b42/detA*q[1]  +  b43/detA*q[2]  +  b44/detA*q[3];

	q1coef[ip] = (b11 + b21*image_point_p_x[ip] + b31*image_point_p_y[ip] + b41*image_point_p_x[ip]*image_point_p_y[ip])/detA;
	q2coef[ip] = (b12 + b22*image_point_p_x[ip] + b32*image_point_p_y[ip] + b42*image_point_p_x[ip]*image_point_p_y[ip])/detA;
	q3coef[ip] = (b13 + b23*image_point_p_x[ip] + b33*image_point_p_y[ip] + b43*image_point_p_x[ip]*image_point_p_y[ip])/detA;
	q4coef[ip] = (b14 + b24*image_point_p_x[ip] + b34*image_point_p_y[ip] + b44*image_point_p_x[ip]*image_point_p_y[ip])/detA;

	//pressure at the image point
	double image_point_pressure = a0[ip] + a1[ip]*image_point_p_x[ip]    + a2[ip]*image_point_p_y[ip]    + a3[ip] * image_point_p_y[ip]   *image_point_p_x[ip];
	body_intercept_p[ip]        = a0[ip] + a1[ip]*body_intercept_p_x[ip] + a2[ip]*body_intercept_p_y[ip] + a3[ip] * body_intercept_p_x[ip]*body_intercept_p_y[ip]; //used for force calc

	//extrapolate pressure to the ghost node
	if (set)
		pressure[ip] = image_point_pressure ;//+ sqrt(pow(image_point_p_x[ip]-xv[I],2)+pow(image_point_p_y[ip]-yu[J],2))*matDClose;
}
}
