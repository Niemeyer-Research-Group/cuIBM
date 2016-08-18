/***************************************************************************//**
 * \file calculateForce.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief
 */

#include "calculateForce.h"

namespace kernels
{
__global__//kernel should be of size totalPoints
void force_pressure(double *force_pressure, double *body_intercept_p,
					double *body_intercept_p_x, double *body_intercept_p_y,
					double *bx, double *by, double *xv, double *yu, int *ghostTagsP,
					int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY)
{
	//initialise
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >= totalPoints)
		return;
	int ii = i_start[0],
		jj = j_start[0],
		ip;
	int I0,
		If,
		J0,
		Jf,
		count = 0,
		thetaID,
		minID,
		maxID;
	//find 16 surrounding nodes
	double	theta,
			thetaNode,
			min=-10,
			max=10;

	while (xv[ii] < bx[idx])
		ii++;
	I0 = ii-2;
	If = ii+1;
	while (yu[jj] <by[idx])
		jj++;
	J0=jj-2;
	Jf=jj+1;
	thetaNode = asin((by[idx]-midY)/sqrt(pow(bx[idx]-midX,2)+pow(by[idx]-midY,2)));
	if (bx[idx] < midX)//this janky if statement forces theta to be continuous
	{
		thetaNode = M_PI-thetaNode;
	}
	if (thetaNode > M_PI*5/4 || thetaNode < -M_PI/4)
	{
		if(bx[idx]>midX)
			thetaNode += 2*M_PI;
	}
	//sweep over nodes calculating theta
	//find theta above and below node
	for (int i = I0; i<=If; i++)
	{
		for(int j=J0;j<=Jf;j++)
		{
			ip = j*nx+i;

			if (ghostTagsP[ip]>0)
			{
				theta = asin((body_intercept_p_y[ip]-midY)/sqrt(pow(body_intercept_p_x[ip]-midX,2)+pow(body_intercept_p_y[ip]-midY,2)));
				if (body_intercept_p_x[ip]<midX)
				{
					theta = M_PI-theta;
				}
				if (thetaNode > M_PI*5/4 || thetaNode < -M_PI/4)
					{
						if(body_intercept_p_x[ip]>midX)
							theta += 2*M_PI;
					}
				thetaID = ip;
				if (theta > thetaNode && theta < max)
				{
					max = theta;
					maxID = thetaID;
				}
				if (theta<thetaNode && theta > min)
				{
					min = theta;
					minID = thetaID;
				}
			}
			count ++;
		}
	}
	//interp for node
	force_pressure[idx] = body_intercept_p[minID] + (body_intercept_p[maxID] - body_intercept_p[minID]) * (thetaNode-min) / (max-min);
}

__global__
void force_velocity_x(double *force_dudx, double *uB, double *u,
						double *bx, double *by, double *xu, double *yu,
						int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY, double dx)
{
	//initialise
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >= totalPoints)
		return;

	int ii = i_start[0],
		jj = j_start[0];
	double            y3,        q3, q4,
		   x1, x2,    y1,        q1, q2;

	//find extended image point
	double	rise = by[idx]-midY,
			run = bx[idx]-midX,
			radius = sqrt(rise*rise+run*run);

	double	dn = dx*sqrt(2.0),//distance from body to calc normal at, needs to be at least sqrt(2)*dx to place it a full node away from the body
			ipx = bx[idx] + dn/radius*run,
			ipy = by[idx] + dn/radius*rise;

	//find points bounding extended image point
	while (xu[ii] < ipx)
		ii++;
	x1 = xu[ii-1];	x2 = xu[ii];
	while (yu[jj] < ipy)
		jj++;
	y3 = yu[jj];
	y1 = yu[jj-1];

	q3 = u[(jj)*(nx-1)+(ii-1)];     q4 = u[(jj)*(nx-1)+(ii)];
	q1 = u[(jj-1)*(nx-1)+(ii-1)];   q2 = u[(jj-1)*(nx-1)+(ii)];

	//interp for u at extended image point //flag grid must be uniform
	//http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
	double	topleft = (x2-ipx)*(y1-ipy)/(x2-x1)/(y1-y3)*q3,
			topright = (ipx-x1)*(y1-ipy)/(x2-x1)/(y1-y3)*q4,
			botleft = (x2-ipx)*(ipy-y3)/(x2-x1)/(y1-y3)*q1,
			botright = (ipx-x1)*(ipy-y3)/(x2-x1)/(y1-y3)*q2,
			ipu = botleft + botright + topleft + topright;

	//calc normal derivative
	force_dudx[idx] = (ipu-uB[0])/dn;
}

__global__
void force_velocity_y(double *force_dvdx, double *vB, double *u,
						double *bx, double *by, double *xv, double *yv,
						int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY, double dx)
{
	//initialise
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >= totalPoints)
		return;

	int ii = i_start[0],
		jj = j_start[0];
	double            y3,        q3, q4,
		   x1, x2,    y1,        q1, q2;

	//find extended image point
	double	rise = by[idx]-midY,
			run = bx[idx]-midX,
			radius = sqrt(rise*rise+run*run);

	double	dn = dx*sqrt(2.0),//distance from body to calc normal at, needs to be at least sqrt(2)*dx to place it a full node away from the body
			ipx = bx[idx] + dn/radius*run,
			ipy = by[idx] + dn/radius*rise;

	//find points bounding extended image point
	while (xv[ii] < ipx)
		ii++;
	x1 = xv[ii-1];	x2 = xv[ii];
	while (yv[jj] < ipy)
		jj++;
	y3 = yv[jj];
	y1 = yv[jj-1];

	q3 = u[(jj)*(nx-1)+(ii-1) + ny*(nx-1)];     q4 = u[(jj)*(nx-1)+(ii) + ny*(nx-1)];
	q1 = u[(jj-1)*(nx-1)+(ii-1) + ny*(nx-1)];   q2 = u[(jj-1)*(nx-1)+(ii) + ny*(nx-1)];

	//interp for u at extended image point //flag grid must be uniform
	//http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
	double	topleft = (x2-ipx)*(y1-ipy)/(x2-x1)/(y1-y3)*q3,
			topright = (ipx-x1)*(y1-ipy)/(x2-x1)/(y1-y3)*q4,
			botleft = (x2-ipx)*(ipy-y3)/(x2-x1)/(y1-y3)*q1,
			botright = (ipx-x1)*(ipy-y3)/(x2-x1)/(y1-y3)*q2,
			ipu = botleft + botright + topleft + topright;

	//calc normal derivative
	force_dvdx[idx] = (ipu-vB[0])/dn;
}

__global__
void force(double *force_x, double *force_y, double *pressure, double *dudn, double *dvdn,
			double *bx, double *by,
			int totalPoints, double midX, double midY, double nu)
{
	//initialise
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >= totalPoints)
		return;
	//get area
	double area = sqrt(pow(bx[0]-bx[1],2) + pow(by[0]-by[1],2));
	//get normal vector
	double	h =  sqrt(pow(by[idx]-midY,2) + pow(bx[idx]-midX,2)),
			n1 = (bx[idx]-midX)/h,
			n2 = (by[idx]-midY)/h;
	//calc tang stress
	double	mu = nu,
			tau_x = mu*((1-n1*n1)*dudn[idx]+(-n1*n2)*dvdn[idx]),
			tau_y = mu*((-n1*n2)*dudn[idx]+(1-n2*n2)*dvdn[idx]);
	//integrate
	force_x[idx] = area * tau_x - area * n1 * pressure[idx];
	force_y[idx] = area * tau_y - area * n2 * pressure[idx];
}
}
