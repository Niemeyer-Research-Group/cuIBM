/***************************************************************************//**
 * \file weight.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */

#include "weight.h"

namespace kernels
{
__global__
void weightX(double *uhat, double *ustar, int *ghostTagsUV, int *hybridTagsUV, double *yu, double *xu,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int *i_start, int *j_start, int width, int nx, int ny)
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iu = J*(nx-1) + I;
	if (iu > ny*(nx-1)) //return if we're out of bounds
		return;
	if (hybridTagsUV[iu]<=0) //return if we're not at a hybrid node point
		return;

	double  delta_1 = 0,
			delta_2 = 0,
			alpha = 1,
			dx = xu[I]-xu[I-1],
			dy = yu[J]-yu[J-1];
	//find ghost node in x direction
	//west is inside
	if (ghostTagsUV[iu-1]>0)
		delta_1 = sqrt( pow(body_intercept_x[iu-1]-xu[I-1],2 ) + pow( (body_intercept_y[iu-1]-yu[J]), 2 ) );
	//east is inside
	else if(ghostTagsUV[iu+1]>0)
		delta_1 = sqrt( pow( (body_intercept_x[iu+1]-xu[I+1]),2 ) + pow( (body_intercept_y[iu+1]-yu[J]), 2 ) );
	//find ghost node in y direction
	//south is inside
	if (ghostTagsUV[iu-(nx-1)]>0)
		delta_2 = sqrt( pow( body_intercept_x[iu-(nx-1)]-xu[I],2 ) + pow( body_intercept_y[iu-(nx-1)]-yu[J-1], 2 ) );
	//north is inside
	else if (ghostTagsUV[iu+(nx-1)]>0)
		delta_2 = sqrt( pow( body_intercept_x[iu+(nx-1)]-xu[I],2 ) + pow( body_intercept_y[iu+(nx-1)]-yu[J+1], 2 ) );
	//calculate alpha
	alpha = sqrt( pow( delta_1/dx , 2 ) + pow( delta_2/dy , 2 ) );
	//blend uhat
	uhat[iu] = (1-alpha)*uhat[iu] + alpha*ustar[iu];
}

__global__
void weightY(double *uhat, double *ustar, int *ghostTagsUV, int *hybridTagsUV, double *yv, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int *i_start, int *j_start, int width, int nx, int ny)
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		iv = J*nx + I + (nx-1)*ny;
	if (J*nx + I > nx*(ny-1)) //return if we're out of bound
		return;
	if (hybridTagsUV[iv]<=0) //return if we're not at an interpolation point
		return;

	double  delta_1 = 0,
			delta_2 = 0,
			alpha = 1,
			dx = xv[I]-xv[I-1],
			dy = yv[J]-yv[J-1];
	//find ghost node in x direction
	//west is inside
	if (ghostTagsUV[iv-1]>0)
		delta_1 = sqrt( pow( body_intercept_x[iv-1]-xv[I-1],2 ) + pow( body_intercept_y[iv-1]-yv[J], 2 ) );
	//east is inside
	else if(ghostTagsUV[iv+1]>0)
		delta_1 = sqrt( pow( body_intercept_x[iv+1]-xv[I+1],2 ) + pow( body_intercept_y[iv+1]-yv[J], 2 ) );
	//find ghost node in y direction
	//south is inside
	if (ghostTagsUV[iv-nx]>0)
		delta_2 = sqrt( pow( body_intercept_x[iv-nx]-xv[I],2 ) + pow( body_intercept_y[iv-nx]-yv[J-1], 2 ) );
	//north is inside
	else if (ghostTagsUV[iv+nx]>0)
		delta_2 = sqrt( pow( body_intercept_x[iv+nx]-xv[I],2 ) + pow( body_intercept_y[iv+nx]-yv[J+1], 2 ) );
	//calculate alpha
	alpha = sqrt( pow( delta_1/dx , 2 ) + pow( delta_2/dy , 2 ) );
	//blend uhat
	uhat[iv] = (1-alpha)*uhat[iv] + alpha*ustar[iv];
}

__global__
void alpha_(double *alpha, int *ghostTagsP, int *hybridTagsP, double *yu, double *xv,
				double *body_intercept_x, double *body_intercept_y,
				int *i_start, int *j_start, int width, int nx, int ny)
{
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (width),
		j	= idx / (width),
		I	= i_start[0] + i,
		J	= j_start[0] + j,
		ip = J*nx + I;
	if (J*nx + I > nx*ny) //return if we're out of bound
		return;
	if (hybridTagsP[ip]<=0) //return if we're not at an interpolation point
		return;

	double  delta_1 = 0,
			delta_2 = 0,
			dx = xv[I]-xv[I-1],
			dy = yu[J]-yu[J-1];
	alpha[ip] = 1;
	//find ghost node in x direction
	//west is inside
	if (ghostTagsP[ip-1]>0)
		delta_1 = sqrt( pow( body_intercept_x[ip-1]-xv[I-1],2 ) + pow( body_intercept_y[ip-1]-yu[J], 2 ) );
	//east is inside
	else if(ghostTagsP[ip+1]>0)
		delta_1 = sqrt( pow( body_intercept_x[ip+1]-xv[I+1],2 ) + pow( body_intercept_y[ip+1]-yu[J], 2 ) );
	//find ghost node in y direction
	//south is inside
	if (ghostTagsP[ip-nx]>0)
		delta_2 = sqrt( pow( body_intercept_x[ip-nx]-xv[I],2 ) + pow( body_intercept_y[ip-nx]-yu[J-1], 2 ) );
	//north is inside
	else if (ghostTagsP[ip+nx]>0)
		delta_2 = sqrt( pow( body_intercept_x[ip+nx]-xv[I],2 ) + pow( body_intercept_y[ip+nx]-yu[J+1], 2 ) );
	//calculate alpha
	alpha[ip] = sqrt( pow( delta_1/dx , 2 ) + pow( delta_2/dy , 2 ) );
}

}
