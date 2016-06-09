/***************************************************************************//**
 * \file calculateForce.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief
 */

#include "calculateForce.h"

namespace kernels
{
__global__
void pressure_at_BI(double *force_pressure, double *pressure, double *u, int *ghostTagsP, int *hybridTagsP, double *bx, double *by,
						double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *xv,
						double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4,
						double *point_x, double *point_y, double *point2_x, double *point2_y, double *point3_x, double *point3_y,
						int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY)
{
	//flag this kernel is going to have terrible divergance issues. could be helped by using some sort of padding function. i.e. set it up so we only get one real node per warp
	int idx	= threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >= totalPoints) //return if we're out of bound
		return;
	//find the closest ghost node UV
	int I = 0,
		J = 0,
		ip[4], //index of p nodes
		Ip[4], //I index of p
		Jp[4]; //J index of p
	/*
	 *   	__________(iv3)______________(iv4)_______
	 *   	|					   |				|
	 *   	| 					   |				|
	 *   	(iu3)	  (ip3)      (iu4)   (ip4)   	|
	 *   	|	 				   |				|
	 *   	|					*  |				|
	 *   	|_________(iv2)________|_____(iv2)______|
	 *   	|					   |
	 * 		|					   |
	 * 		(iu1)	  (ip1)		 (iu2)   (ip2)
	 */

	double BI_x, BI_y, BI_u, BI_v, BI_u0, BI_v0;
	if (idx == totalPoints-1)
	{
		BI_x = (bx[idx]+bx[0])/2;
		BI_y = (by[idx]+by[0])/2;
		BI_u = (uB[idx]+uB[0])/2;
		BI_v = (vB[idx]+vB[0])/2;
		BI_u0 = (uB0[idx]+uB0[0])/2;
		BI_v0 = (vB0[idx]+vB0[0])/2;
	}
	else
	{
		BI_x = (bx[idx]+bx[idx+1])/2;
		BI_y = (by[idx]+by[idx+1])/2;
		BI_u = (uB[idx]+uB[idx+1])/2;
		BI_v = (vB[idx]+vB[idx+1])/2;
		BI_u0 = (uB0[idx]+uB0[idx+1])/2;
		BI_v0 = (vB0[idx]+vB0[idx+1])/2;
	}

	//find the four pressure indices surrounding the body intercept as well as the I and J indices for those points
	I = 0;
	J = 0;
	while (xv[I] < BI_x)
		I+=1;
	while (yu[J] < BI_y)
		J+=1;
	Ip[0] = I-1; Ip[1]=I;    Ip[2]=I-1;   Ip[3]=I;
	Jp[0] = J-1; Jp[1]=J-1;  Jp[2]=J;     Jp[3]=J;
	for (int i=0;i<4;i++)
		ip[i]=Jp[i]*nx+Ip[i];

	point_x[idx] = BI_x;
	point_y[idx] = BI_y;
	//find closest hybrid pressure node to the body
	double	s = 1,
			minIDX,
			distance;
	for (int i=0; i<4;i++)
	{
		distance = sqrt( pow(BI_x-xv[Ip[i]],2) + pow(BI_y-yu[Jp[i]],2) );
		if ( distance < s && hybridTagsP[ip[i]] > 0 )
		{
			minIDX = ip[i];
			s = distance;
			point2_x[idx] = xv[Ip[i]];
			point2_y[idx] = yu[Jp[i]];
		}
	}

	//using that distance, s, find the image point that is normal to the point
	double theta;
	if (idx == totalPoints-1)
	{
		distance = sqrt( pow(bx[idx]-bx[0],2) + pow(by[idx]-by[0],2) );
		theta = asin((by[0]-by[idx])/distance);
	}
	else
	{
		distance = sqrt( pow(bx[idx]-bx[idx+1],2) + pow(by[idx]-by[idx+1],2) );
		theta = asin((by[idx+1]-by[idx])/distance);
	}
	s = sqrt( pow(xv[Ip[0]]-xv[Ip[1]],2) + pow(yu[Jp[0]]-yu[Jp[2]],2) );
	double ip_x, ip_y;
	if (BI_y>midY)
	{
		ip_x = cos(M_PI/2-theta)*s+BI_x;
		ip_y = sin(M_PI/2-theta)*s+BI_y;
	}
	else
	{
		ip_x = BI_x-cos(M_PI/2+theta)*s;
		ip_y = BI_y-sin(M_PI/2+theta)*s;
	}
	point3_x[idx] = ip_x;
	point3_y[idx] = ip_y;

	//find the four nodes that surround the image point
	int		ii = Ip[0]-2,
			jj = Jp[0]-2;
	while (xv[ii] < ip_x)
		ii++;
	x1[idx] = xv[ii-1];
	x2[idx] = xv[ii];
	x3[idx] = x1[idx];
	x4[idx] = x2[idx];
	while (yu[jj] < ip_y)
		jj++;
	y1[idx] = yu[jj-1];
	y2[idx] = y1[idx];
	y3[idx] = yu[jj];
	y4[idx] = y3[idx];

	int		Ip2[4],
			Jp2[4],
			ip2[4];
	Ip2[0] = ii-1; Ip2[1] = ii; Ip2[2] = ii-1; Ip2[3] = ii;
	Jp2[0] = jj-1; Jp2[1] = jj-1; Jp2[2] = jj; Jp2[3] = jj;
	for (int i=0;i<4;i++)
	{
		ip2[i] = Jp2[i]*nx + Ip2[i];
	}
	q1[idx] = pressure[ip2[0]];
	q2[idx] = pressure[ip2[1]];
	q3[idx] = pressure[ip2[2]];
	q4[idx] = pressure[ip2[3]];

	//bilinear interpolation for pressure
	double a11 = 1,  a12 = xv[Ip2[0]],  a13 = yu[Jp2[0]],  a14 = xv[Ip2[0]]*yu[Jp2[0]];
	double a21 = 1,  a22 = xv[Ip2[1]],  a23 = yu[Jp2[1]],  a24 = xv[Ip2[1]]*yu[Jp2[1]];
	double a41 = 1,  a42 = xv[Ip2[3]],  a43 = yu[Jp2[3]],  a44 = xv[Ip2[3]]*yu[Jp2[3]];
	double a31 = 1,  a32 = xv[Ip2[2]],  a33 = yu[Jp2[2]],  a34 = xv[Ip2[2]]*yu[Jp2[2]];

	//find the cloest node to the BI
	s=1;
	for (int i=0; i<4;i++)
	{
		distance = sqrt( pow(BI_x-xv[Ip2[i]],2) + pow(BI_y-yu[Jp2[i]],2) );
		if ( distance < s)
		{
			minIDX = ip2[i];
			s = distance;
		}
	}
	//if the node is the cloesest to the body, move it to the BI and set it to be a neuman condition
	//point 1
	double n_x, n_y, nl,du_dt, u_du_dx, v_du_dy, dv_dt, u_dv_dx, v_dv_dy;
	int iu, iv;
	if (ip2[0] == minIDX)
	{
		iu = Jp2[0]*(nx-1) + Ip2[0];
		iv = Jp2[0]*nx + Ip2[0] + (nx-1)*ny;
		x1[idx] = BI_x;
		y1[idx] = BI_y;
		n_x = ip_x - x1[idx];
		n_y = ip_y - y1[idx];
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = BI_u - BI_u0;
		u_du_dx = BI_u * ( (u[iu]+u[iu-1]+u[iu-(nx-1)]+u[iu-(nx-1)-1]) /4 - BI_u )  /  (xv[Ip2[0]]-x1[idx]);
		v_du_dy = BI_v * (u[iu-1]-BI_u)/(yu[Jp2[0]]-y1[idx]);
		dv_dt = BI_v - BI_v0;
		u_dv_dx = BI_u*(u[iv-nx]-BI_v)/(xv[Ip2[0]]-x1[idx]);
		v_dv_dy = BI_v*((u[iv]+u[iv-nx]+u[iv-1]+u[iv-nx-1])/4 - BI_v)/(yu[Jp2[0]]-y1[idx]);
		q1[idx] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt + u_dv_dx + v_dv_dy));
		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x1[idx]+a12*y1[idx];
	}
	//point 2
	if (ip2[1] == minIDX)
	{
		iu = Jp2[1]*(nx-1) + Ip2[1],
		iv = Jp2[1]*nx + Ip2[1] + (nx-1)*ny;
		x2[idx] = BI_x;
		y2[idx] = BI_y;
		n_x = ip_x - x2[idx];
		n_y = ip_y - y2[idx];
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = BI_u - BI_u0;
		u_du_dx = BI_u*(BI_u - (u[iu]+u[iu-1]+u[iu-(nx-1)]+u[iu-(nx-1)-1])/4)/(x2[idx]-xv[Ip2[1]]); //flag this approximation of du/dx might be too rough as we are calculating our u values at different heights. One point is the body intercept and the second point is the v node between poitns 1 and 2
		v_du_dy = BI_v*(u[iu] - BI_u)/(yu[Jp2[1]]-y2[idx]);
		dv_dt = BI_v - BI_v0;
		u_dv_dx = BI_u*(BI_v-u[iv-nx])/(x2[idx]-xv[Ip2[1]]);
		v_dv_dy = BI_v*((u[iv]+u[iv-nx]+u[iv+1]+u[iv-nx+1])/4 - BI_v)/(yu[Jp2[1]]-y2[idx]);
		q2[idx] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x2[idx]+a22*y2[idx];
	}
	//point 3
	if (ip2[2] == minIDX)
	{
		iu = Jp2[2]*(nx-1) + Ip2[2];
		iv = Jp2[2]*nx + Ip2[2] + (nx-1)*ny;
		x3[idx] = BI_x;
		y3[idx] = BI_y;
		n_x = ip_x - x3[idx];
		n_y = ip_y - y3[idx];
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = BI_u - BI_u0;
		u_du_dx = BI_u*((u[iu]+u[iu-1]+u[iu+(nx-1)]+u[iu+(nx-1)-1])/4 - BI_u)/(xv[Ip2[2]]-x3[idx]);
		v_du_dy = BI_v*(BI_u - u[iu-1])/(y3[idx] - yu[Jp2[2]]);
		dv_dt = BI_v - BI_v0;
		u_dv_dx = BI_u*(u[iv]-BI_v)/(xv[Ip2[2]]-x3[idx]);
		v_dv_dy = BI_v*(BI_v- (u[iv]+u[iv-nx]+u[iv-1]+u[iv-nx-1])/4)/(y3[idx]-yu[Jp2[2]]);
		q3[idx] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x3[idx]+a32*y3[idx];
	}
	//4
	if (ip2[3] == minIDX)
	{
		iu = Jp2[3]*(nx-1) + Ip2[3];
		iv = Jp2[3]*nx + Ip2[3] + (nx-1)*ny;
		x4[idx] = BI_x;
		y4[idx] = BI_y;
		n_x = ip_x - x4[idx];
		n_y = ip_y - y4[idx];
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = BI_u - BI_u0;
		u_du_dx = BI_u*(BI_u - (u[iu]+u[iu-1]+u[iu+(nx-1)]+u[iu+(nx-1)-1])/4)/(x4[idx]-xv[Ip2[3]]); //flag this approximation of du/dx might be too rough as we are calculating our u values at different heights. One point is the body intercept and the second point is the v node between poitns 1 and 2
		v_du_dy = BI_v*(BI_u - u[iu])/(y4[idx]-yu[Jp2[3]]);
		dv_dt = BI_v - BI_v0;
		u_dv_dx = BI_u*(BI_v-u[iv])/(x4[idx]-xv[Ip2[3]]);
		v_dv_dy = BI_v*(BI_v- (u[iv]+u[iv-nx]+u[iv+1]+u[iv-nx+1])/4)/(y4[idx]-yu[Jp2[3]]);
		q4[idx] = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a41 = 0;
		a42 = n_x/nl;
		a43 = n_y/nl;
		a44 = a43*x4[idx]+a42*y4[idx];
	}

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

	double a0 = b11/detA*q1[idx]  +  b12/detA*q2[idx]  +  b13/detA*q3[idx]  +  b14/detA*q4[idx];
	double a1 = b21/detA*q1[idx]  +  b22/detA*q2[idx]  +  b23/detA*q3[idx]  +  b24/detA*q4[idx];
	double a2 = b31/detA*q1[idx]  +  b32/detA*q2[idx]  +  b33/detA*q3[idx]  +  b34/detA*q4[idx];
	double a3 = b41/detA*q1[idx]  +  b42/detA*q2[idx]  +  b43/detA*q3[idx]  +  b44/detA*q4[idx];
	//interpolate pressure to the body intercept
	double matD = 0;
	if (minIDX == ip2[0])
		matD = -q1[idx];
	if (minIDX == ip2[1])
		matD = -q2[idx];
	if (minIDX == ip2[2])
		matD = -q3[idx];
	if (minIDX == ip2[3])
		matD = -q4[idx];
	force_pressure[idx] = a0 + a1*ip_x + a2*ip_y + a3*ip_y*ip_x
			   	        + sqrt(pow(ip_x-BI_x,2)+pow(ip_y-BI_y,2))*matD; //flag don't need to get ip then interpolate, can just use the interpoalter at the BI

}

/*
__global__
void luoForce()
{
	//flag this kernel is going to have terrible divergance issues. could be helped by using some sort of padding function. i.e. set it up so we only get one real node per warp
	int idx	= threadIdx.x + blockDim.x * blockIdx.x;
	if (idx > totalPoints) //return if we're out of bound
		return;
	//find the closest ghost node UV
	int I = 0,
		J = 0,
		iu[4], //index of u nodes
		iv[4], //index of v nodes
		ip[4], //index of p nodes
		Iu[4], //I index of u nodes
		Iv[4], //I index of v
		Ip[4], //I index of p
		Ju[4], //J index of u nodes
		Jv[4], //J index of v
		Jp[4]; //J index of p
	/*
	 *   	__________(iv3)______________(iv4)_______
	 *   	|					   |				|
	 *   	| 					   |				|
	 *   	(iu3)	  (ip3)      (iu4)   (ip4)   	|
	 *   	|	 				   |				|
	 *   	|					*  |				|
	 *   	|_________(iv2)________|_____(iv2)______|
	 *   	|					   |
	 * 		|					   |
	 * 		(iu1)	  (ip1)		 (iu2)   (ip2)
	 */

/*
	//find ius
	while (xu[I] < bx[idx])
		I+=1;
	while (yu[J] < by[idx])
		J+=1;
	Iu[0] = I; Iu[1]=I+1; Iu[2]=I;   Iu[3]=I+1;
	Ju[0] = J; Ju[1]=J;   Ju[2]=J+1; Ju[3]=J+1;
	for (int i=0;i<4;i++)
		iu[i]=Ju[i]*(nx-1)+Iu[i];

	//find ivs
	I = 0;
	J = 0;
	while (xv[I] < bx[idx])
		I+=1;
	while (yv[J] < by[idx])
		J+=1;
	Iv[0] = I; Iv[1]=I+1; Iv[2]=I;   Iv[3]=I+1;
	Jv[0] = J; Jv[1]=J;   Jv[2]=J+1; Jv[3]=J+1;
	for (int i=0;i<4;i++)
		iv[i]=Jv[i]*nx+Iv[i] + (nx-1)*ny;

	//find ips
	I = 0;
	J = 0;
	while (xv[I] < bx[idx])
		I+=1;
	while (yu[J] < by[idx])
		J+=1;
	Ip[0] = I; Ip[1]=I+1; Ip[2]=I;   Ip[3]=I+1;
	Jp[0] = J; Jp[1]=J;   Jp[2]=J+1; Jp[3]=J+1;
	for (int i=0;i<4;i++)
		ip[i]=Jp[i]*nx+Ip[i];


	//find closest ghost node iu
	double	s[3],
			minIDX[3],
			distance;
	s[0] = 1; s[1] = 1; s[2] = 1;
	for (int i=0; i<4;i++)
	{
		distance = sqrt( pow(bx[idx]-xu[Iu[i]],2) + pow(by[idx]-yu[Ju[i]],2) );
		if (distance < s[0] && ghostTagsUV[iu[i]] > 0)
			minIDX[0] = iu[i];
		distance = sqrt( pow(bx[idx]-xv[Iu[i]],2) + pow(by[idx]-yv[Ju[i]],2) );
		if (distance < s[1] && ghostTagsUV[iv[i]] > 0)
			minIDX[1] = iv[i];
		distance = sqrt( pow(bx[idx]-xv[Iu[i]],2) + pow(by[idx]-yu[Ju[i]],2) );
		if (distance < s[2] && ghostTagsP[ip[i]] > 0)
			minIDX[2] = ip[i];
	}
	//bilinear interpolation for pressure at the nodes
	double	q1 = pressure[ip[0]],
			q2 = pressure[ip[1]],
			q3 = pressure[ip[2]],
			q4 = pressure[ip[3]];

	double a11 = 1, a12 = xv[Ip[0]],  a13 = yu[Jp[0]], a14 = xv[Ip[0]]*yu[Jp[0]];
	double a21 = 1, a22 = xv[Ip[1]],  a23 = yu[Jp[1]], a24 = xv[Ip[1]]*yu[Jp[1]];
	double a31 = 1, a32 = xv[Ip[2]],  a33 = yu[Jp[2]], a34 = xv[Ip[2]]*yu[Jp[2]];
	double a41 = 1, a42 = xv[Ip[3]],  a43 = yu[Jp[3]], a44 = xv[Ip[3]]*yu[Jp[3]];

	//if the node is inside the body, set it to be a neuman condition
	//point 1
	double x, y, n_x,n_y,nl_du_dt,u_du_dx,v_du_dy,dv_dt,u_dv_dx_v_dv_dy;
	if (ghostTagsP[ip[0]] != -1)
	{
		double idxu = Ju[0]*(nx-1) + Iu[0],
			   idxv = Ju[0]*nx + Iu[0] + (nx-1)*ny;
		x = body_intercept_p_x[ip[0]];
		y = body_intercept_p_y[ip[0]];
		n_x = image_point_p_x[ip[0]] - x[ip[0]];
		n_y = image_point_p_y[ip[0]] - y[ip[0]];
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = uB[idx] - uB0[idx];
		u_du_dx = uB[idx]*((u[idxu]+u[idxu-1]+u[idxu-(nx-1)]+u[idxu-(nx-1)-1])/4 - uB[idx])/(xv[Ip[0]]-x);
		v_du_dy = vB[idx]*(u[idxu-1]-uB[idx])/(yu[Jp[0]]-y[ip[0]]);
		dv_dt = vB[idx] - vB0[idx];
		u_dv_dx = uB[idx]*(u[idxv-nx]-vB[idx])/(xv[Ip[0]]-x);
		v_dv_dy = vB[idx]*((u[idxv]+u[idxv-nx]+u[idxv-1]+u[idxv-nx-1])/4 - vB[idx])/(yu[Jp[0]]-y);
		q1 = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt + u_dv_dx + v_dv_dy));
		a11 = 0;
		a12 = n_x/nl;
		a13 = n_y/nl;
		a14 = a13*x+a12*y;
	}
	//point 2
	if (ghostTagsP[ip[1]] != -1)
	{
		double idxu = Ju[1]*(nx-1) + Iu[1],
			   idxv = Ju[1]*nx + Iu[1] + (nx-1)*ny;
		x = body_intercept_p_x[ip[1]];
		y = body_intercept_p_y[ip[1]];
		n_x = image_point_p_x[ip[1]] - x;
		n_y = image_point_p_y[ip[1]] - y;
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = uB[idx] - uB0[idx];
		u_du_dx = uB[idx]*(uB[idx] - (u[idxu]+u[idxu-1]+u[idxu-(nx-1)]+u[idxu-(nx-1)-1])/4)/(x-xv[Ip[1]]); //flag this approximation of du/dx might be too rough as we are calculating our u values at different heights. One point is the body intercept and the second point is the v node between poitns 1 and 2
		v_du_dy = vB[idx]*(u[idxu] - uB[idx])/(yu[Jp[1]]-y);
		dv_dt = vB[idx] - vB0[idx];
		u_dv_dx = uB[idx]*(vB[idx]-u[idxv-nx])/(x-xv[Ip[1]]);
		v_dv_dy = vB[idx]*((u[idxv]+u[idxv-nx]+u[idxv+1]+u[idxv-nx+1])/4 - vB[idx])/(yu[Jp[1]]-y);
		q2 = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a21 = 0;
		a22 = n_x/nl;
		a23 = n_y/nl;
		a24 = a23*x+a22*y;
	}
	//point 3
	if (ghostTagsP[ip[2]] != -1)
	{
		double idxu = Ju[2]*(nx-1) + Iu[2],
			   idxv = Ju[2]*nx + Iu[2] + (nx-1)*ny;
		x = body_intercept_p_x[ip[2]];
		y = body_intercept_p_y[ip[2]];
		n_x = image_point_p_x[ip[2]] - x;
		n_y = image_point_p_y[ip[2]] - y;
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = uB[idx] - uB0[idx];
		u_du_dx = uB[idx]*((u[idxu]+u[idxu-1]+u[idxu+(nx-1)]+u[idxu+(nx-1)-1])/4 - uB[idx])/(xv[Ip[2]]-x); //flag this approximation of du/dx might be too rough as we are calculating our u values at different heights. One point is the body intercept and the second point is the v node between poitns 1 and 2
		v_du_dy = vB[idx]*(uB[idx] - u[idxu-1])/(y - yu[Jp[2]]);
		dv_dt = vB[idx] - vB0[idx];
		u_dv_dx = uB[idx]*(u[idxv]-vB[idx])/(xv[Ip[2]]-x);
		v_dv_dy = vB[idx]*(vB[idx]- (u[idxv]+u[idxv-nx]+u[idxv-1]+u[idxv-nx-1])/4)/(y-yu[Jp[2]]);
		q3 = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a31 = 0;
		a32 = n_x/nl;
		a33 = n_y/nl;
		a34 = a33*x+a32*y;
	}
	//4
	if (ghostTagsP[ip[3]] != -1)
	{
		x = body_intercept_p_x[ip[3]];
		y = body_intercept_p_y[ip[3]];
		n_x = image_point_p_x[ip[3]] - x;
		n_y = image_point_p_y[ip[3]] - y;
		nl = sqrt(n_x*n_x+n_y*n_y);
		du_dt = uB[idx] - uB0[idx];
		u_du_dx = uB[idx]*(uB[idx] - (u[idxu]+u[idxu-1]+u[idxu+(nx-1)]+u[idxu+(nx-1)-1])/4)/(x-xv[Ip[3]]); //flag this approximation of du/dx might be too rough as we are calculating our u values at different heights. One point is the body intercept and the second point is the v node between poitns 1 and 2
		v_du_dy = vB[idx]*(uB[idx] - u[idxu])/(y-yu[Jp[3]]);
		dv_dt = vB[idx] - vB0[idx];
		u_dv_dx = uB[idx]*(vB[idx]-u[idxv])/(x-xv[Ip[3]]);
		v_dv_dy = vB[idx]*(vB[idx]- (u[idxv]+u[idxv-nx]+u[idxv+1]+u[idxv-nx+1])/4)/(y-yu[Jp[3]]);
		q4 = -(n_x/nl*(du_dt+u_du_dx+v_du_dy) + n_y/nl*(dv_dt+u_dv_dx + v_dv_dy));
		a41 = 0;
		a42 = n_x/nl;
		a43 = n_y/nl;
		a44 = a43*x+a42*y;
	}

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

	 double a0 = b11/detA*q1  +  b12/detA*q2  +  b13/detA*q3  +  b14/detA*q4;
	 double a1 = b21/detA*q1  +  b22/detA*q2  +  b23/detA*q3  +  b24/detA*q4;
	 double a2 = b31/detA*q1  +  b32/detA*q2  +  b33/detA*q3  +  b34/detA*q4;
	 double a3 = b41/detA*q1  +  b42/detA*q2  +  b43/detA*q3  +  b44/detA*q4;
	 //interpolate pressure to the body intercept
	 double matD = 0;
	 if (minIDX[2] == ip[0])
		 matD = -q1;
	 if (minIDX[2] == ip[1])
		 matD = -q2;
	 if (minIDX[2] == ip[2])
		 matD = -q3;
	 if (minIDX[2] == ip[3])
		 matD = -q4;
	 double pressure = a0 + a1*image_point_p_x[ip] + a2*image_point_p_y[ip] + a3*image_point_p_y[ip]*image_point_p_x[ip]
	                 + sqrt(pow(image_point_p_x[ip]-bx[idx],2)+pow(image_point_p_y[ip]-by[idx],2))*matD; //flag don't need to get ip then interpolate, can just use the interpoalter at the BI

	//solve for u velocity at the image point
	//find x
	double	x1 = xu[Iu[0]],
			x2 = xu[Iu[1]],
			x3 = xu[Iu[2]],
			x4 = xu[Iu[3]];

	double	y1 = yu[Ju[0]],
			y2 = yu[Ju[1]],
			y3 = yu[Ju[2]],
			y4 = yu[Ju[3]];
	//find q1,q2,q3,q4
	q1 = u[Iu[0]];
	q2 = u[Iu[1]];
	q3 = u[Iu[2]];
	q4 = u[Iu[3]];

	//check if any points are inside of the body, then move them to the body intercept
	//point 1
	if (hybridTagsUV[Iu[0]] == iu[0])
	{
		x1 = bx[idx];
		y1 = by[idx];
		q1 = uB[idx];
	}
	if (hybridTagsUV[Iu[1]] == iu[1])
	{
		x2 = bx[idx];
		y2 = by[idx];
		q2 = uB[idx];
	}
	if (hybridTagsUV[Iu[2]] == iu[2])
	{
		x3 = bx[idx];
		y3 = by[idx];
		q3 = uB[idx];
	}
	if (hybridTagsUV[Iu[3]] == iu[3])
	{
		x4 = bx[idx];
		y4 = by[idx];
		q4 = uB[idx];
	}
//YOU ARE HERE
//we don't have the image point locations for the normal point atm.
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
	 /*double a0 = b11/detA*q1[iv]  +  b12/detA*q2[iv]  +  b13/detA*q3[iv]  +  b14/detA*q4[iv];
	 double a1 = b21/detA*q1[iv]  +  b22/detA*q2[iv]  +  b23/detA*q3[iv]  +  b24/detA*q4[iv];
	 double a2 = b31/detA*q1[iv]  +  b32/detA*q2[iv]  +  b33/detA*q3[iv]  +  b34/detA*q4[iv];
	 double a3 = b41/detA*q1[iv]  +  b42/detA*q2[iv]  +  b43/detA*q3[iv]  +  b44/detA*q4[iv];
	 ustar[iv] = a0 + a1*xv[I] + a2*yv[J] + a3*yv[J]*xv[I];







}*/
}
