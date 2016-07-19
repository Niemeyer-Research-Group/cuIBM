/***************************************************************************//**
 * \file intermediatePressure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side of the poission equation
 */

#include "intermediatePressure.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void intermediatePressure_luo(double *rhs2, double *detA, int *hybridTagsP, double *alpha, double *stencilCoef,
								double *xv, double *yu,
								double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
								double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
								double *q1, double *q2, double *q3, double *q4,
								bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag,
								int nx, int ny)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny || hybridTagsP[ip]<=0)
		return;
	int	I	= ip % nx,
		J	= ip / nx;

	double temp = 0;
	double	x = xv[I],
			y = yu[J];

	if (q1flag[ip] == true)
	{
		temp = (b11[ip] + b21[ip]*x + b31[ip]*y + b41[ip]*x*y)*q1[ip]/detA[ip];
	}
	else if (q2flag[ip] == true)
	{
		temp = (b12[ip] + b22[ip]*x + b32[ip]*y + b42[ip]*x*y)*q2[ip]/detA[ip];
	}
	else if (q3flag[ip] == true)
	{
		temp = (b13[ip] + b23[ip]*x + b33[ip]*y + b43[ip]*x*y)*q3[ip]/detA[ip];
	}
	else if (q4flag[ip] == true)
	{
		temp = (b14[ip] + b24[ip]*x + b34[ip]*y + b44[ip]*x*y)*q4[ip]/detA[ip];
	}
	rhs2[ip] = (1-alpha[ip])*rhs2[ip]/stencilCoef[ip] + alpha[ip]*temp;
}

//calc b, det(A), q, qflag, index
__global__
void interpolate_P_HN_setup(double *detA, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0,
									double *yu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny, double dt, double totalPoints,
									double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
									double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
									double *q1, double *q2, double *q3, double *q4,
									bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag,
									int *index1, int *index2, int *index3, int *index4,
									double *x1, double *x2, double *x3, double *x4,
									double *y1, double *y2, double *y3, double *y4,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy)//test
{//flag u not used anymore
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		ip = J*nx + I,
		ii= I-5,
		jj = J-5;
	//if (ip > J*nx + I) //return if we're out of bound
	if (ip >= nx*ny)
		return;
	if (hybridTagsP[ip]<=0) //return if we're not at an interpolation point
		return;

	double	n_x,				//distance from ip_x to BI_x
			n_y,				//distance from ip_y to BI_y
			nl,					//distance from ip to BI
			distance,			//placeholder
			distance2,			//placeholder
			min,				//minimum distance from BI to body node (used in finding closest body node to the BI)
			min2,				//second closest distance
			matDi,				//du/dt at BN 1
			matDj,				//dv/dt at BN 1
			matD2i,				//du/dt at BN 2
			matD2j,				//dv/dt at BN 2
			matDBIi,			//du/dt at BI
			matDBIj;			//dv/dt at BI

	int		bodyindex,			//index of body node 1
			bodyindex2;			//index of body node 2

	double	a11, a12, a13, a14,
			a21, a22, a23, a24,
			a31, a32, a33, a34,
			a41, a42, a43, a44;
	/*
	 * an example of a node setup you might find, this would occur on the top right section of the body
	 * In this case, we would use the field data from points 2,3 and 4 alongside the bc at BI to interpolate for a value at point 1
	 * In the general case, field data from 3 nodes will be used alongside the bc at the node closest to the body (node ip in this kernel)
	 *   	(x3,y3)__________(x4,y4)
	 *   	|					   |
	 *   	| 					   |
	 *   	|					   |
	 *   	|	 *ip			   |
	 *\   	|					   |
	 * \  	(x1,y1)__________(x2,y2)
	 *  \
	 *   *(BI_x,BI_y)
	 *    \
	 *
	 *Here are some references for solving the equation for bilinear interpolation of values
	 *http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html  //for solving a 4x4 matrix exactly
	 *https://www.physicsforums.com/threads/is-normal-derivative-a-definition.706458/   //for dealing with the normal at the boundary
	 *We want to solving the following equation for the a coefficients using the data from the for corenrs
	 *p= @(X,Y) a0 + a1*X + a2*Y + a3*X*Y;
	 *This results in a system of equations (before applying the neuman BC) that looks like this:
	 *  	   A             a			 q
	 *  |1	x1	y1	x1y1|	|a0|	=	|q1|
	 *  |1	x2	y2	x2y2|	|a1|	=	|q2|
	 *  |1	x3	y3	x3y3|	|a2|	=	|q3|
	 *  |1	x4	y4	x4y4|	|a3|	=	|q4|
	 *
	 * when there is a neumann BC, the equation will become:
	 * |0  nx  ny  ny*x+nx*y|  |a|  |Du/Dt . n|
	 * nx and ny represent the unit vector compoents at the BI (in the code they are n_x/nl and n_y/nl

	 *
	 *         A
	 *  |a11	a12		a13		a14|
	 *  |a21	a22		a23		a24|
	 *  |a31	a13		a33		a34|
	 *  |a41	a14		a43		a44|
	 *
	 *
	 *  	   B
	 * |b11 b12 b13 b14|
	 * |b21 b22 b23 b24|
	 * |b31 b32 b33 b34|
	 * |b41 b42 b43 b44|
	 *
	 * Ainv = B/det(A)
	 * a = Ainv*q';
	 */

	//find x and y locations of nodes 1,2,3,4 by finding the nodes that bound the image point
	while (xv[ii] < image_point_p_x[ip])
		ii++;
	while (yu[jj] <image_point_p_y[ip])
		jj++;

	//set x values at corners
	x1[ip] = xv[ii-1];
	x2[ip] = xv[ii];
	x3[ip] = x1[ip];
	x4[ip] = x2[ip];

	//set y values at corners
	y1[ip] = yu[jj-1];
	y2[ip] = y1[ip];
	y3[ip] = yu[jj];
	y4[ip] = y3[ip];

	//set index values
	index1[ip] = (jj-1)*nx+ii-1,
	index2[ip] = (jj-1)*nx+ii,
	index3[ip] = jj*nx+ii-1,
	index4[ip] = jj*nx+ii;

	//
	q1flag[ip] = false;
	q2flag[ip] = false;
	q3flag[ip] = false;
	q4flag[ip] = false;

	a11 = 1,  a12 = x1[ip],  a13 = y1[ip],  a14 = x1[ip]*y1[ip];
	a21 = 1,  a22 = x2[ip],  a23 = y2[ip],  a24 = x2[ip]*y2[ip];
	a31 = 1,  a32 = x3[ip],  a33 = y3[ip],  a34 = x3[ip]*y3[ip];
	a41 = 1,  a42 = x4[ip],  a43 = y4[ip],  a44 = x4[ip]*y4[ip];

	//setup for neuman BC

	//move the closes node to the body to the surface then calculate the neuman boundary condition for it
	//point 1
	if (hybridTagsP[index1[ip]] == ip)
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

		q1flag[ip] = true;
		q1[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[ip]+a12*y1[ip];
	}
	//point 2
	else if (hybridTagsP[index2[ip]] == ip)
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

		q2flag[ip] = true;
		q2[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[ip]+a22*y2[ip];
	}
	//point 3
	else if (hybridTagsP[index3[ip]] == ip)
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

		q3flag[ip] = true;
		q3[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x3[ip]+a32*y3[ip];
	}
	//4
	if (hybridTagsP[index4[ip]] == ip)
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

		q4flag[ip] = true;
		q4[ip] = - ( matDBIi*n_x/nl + matDBIj*n_y/nl ) ;

		a41 = 0;
		a42 = n_x/nl;
		a43 = n_y/nl;
		a44 = a43*x4[ip]+a42*y4[ip];
	}

	detA[ip] = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43
		  +a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41
		  +a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42
		  +a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41
		  -a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42
		  -a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43
		  -a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41
		  -a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;

	b11[ip] = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	b12[ip] = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	b13[ip] = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	b14[ip] = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	b21[ip] = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
	b22[ip] = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
	b23[ip] = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
	b24[ip] = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
	b31[ip] = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
	b32[ip] = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
	b33[ip] = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
	b34[ip] = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
	b41[ip] = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
	b42[ip] = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
	b43[ip] = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
	b44[ip] = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
}

//this kernel is bad
__global__
void hybridPressureNodeCount(int *countD, int *index1, int *index2, int *index3, int *index4, int *hybridTagsP,
								int *i_start, int *j_start, int width, int height, int nx, int ny)
{
	int ip = 0, count = 0;
	for (int j=j_start[0]; j<j_start[0]+height; j++)
	{
	for (int i=i_start[0]; i<i_start[0]+width; i++)
	{
		ip = j*nx+i;
		if (hybridTagsP[ip]>0)
		{
			if (index1[ip] != ip+nx && index1[ip] != ip-1 && index1[ip] != ip && index1[ip] != ip+1 && index1[ip] != ip-nx)
			{
				count+=1;
			}
			if (index2[ip] != ip+nx && index2[ip] != ip-1 && index2[ip] != ip && index2[ip] != ip+1 && index2[ip] != ip-nx)
			{
				count+=1;
			}
			if (index3[ip] != ip+nx && index3[ip] != ip-1 && index3[ip] != ip && index3[ip] != ip+1 && index3[ip] != ip-nx)
			{
				count+=1;
			}
			if (index4[ip] != ip+nx && index4[ip] != ip-1 && index4[ip] != ip && index4[ip] != ip+1 && index4[ip] != ip-nx)
			{
				count+=1;
			}
			countD[ip] = count;
		}
	}
	}
}


}
