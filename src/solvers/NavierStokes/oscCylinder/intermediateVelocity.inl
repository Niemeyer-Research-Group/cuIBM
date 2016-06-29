/***************************************************************************//**
 * \file intermediateVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/oscCylinder/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/luoIBM/kernels/biLinearInterpolation.h> //interpolate


void oscCylinder::preRHS1Interpolation()
{
	logger.startTimer("RHS1 Interpolation");
	double	*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*ustar_r	= thrust::raw_pointer_cast ( &(ustar[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),//not sure if these are on the host or not
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) );
	
	double	*x1_r = thrust::raw_pointer_cast ( &(x1[0]) ),
			*x2_r = thrust::raw_pointer_cast ( &(x2[0]) ),
			*x3_r = thrust::raw_pointer_cast ( &(x3[0]) ),
			*x4_r = thrust::raw_pointer_cast ( &(x4[0]) ),
			*y1_r = thrust::raw_pointer_cast ( &(y1[0]) ),
			*y2_r = thrust::raw_pointer_cast ( &(y2[0]) ),
			*y3_r = thrust::raw_pointer_cast ( &(y3[0]) ),
			*y4_r = thrust::raw_pointer_cast ( &(y4[0]) ),
			*q1_r = thrust::raw_pointer_cast ( &(q1[0]) ),
			*q2_r = thrust::raw_pointer_cast ( &(q2[0]) ),
			*q3_r = thrust::raw_pointer_cast ( &(q3[0]) ),
			*q4_r = thrust::raw_pointer_cast ( &(q4[0]) ),
			*ip_u_r = thrust::raw_pointer_cast( &(ip_u[0]) );
	
	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
			*hybridTagsUV_r		= thrust::raw_pointer_cast ( &(hybridTagsUV[0]) );
		
	parameterDB  &db = *paramDB;
	double	nu = db["flow"]["nu"].get<double>();
	double	dt = db["simulation"]["dt"].get<double>();
		
	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		//width_i = B.numCellsX[0],
		//height_j = B.numCellsY[0];
		width_i = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height_j=B.numCellsYHost;  //this is done because we need the value on the host to calculate the grid size, but copying it to the host every TS is expensive

	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) );
	
	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	//interpolate velocity to image point and ghost node
	
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
	setVelocityInside();
	//interpolate velocity to hybrid node
	cusp::blas::fill(ustar,0.0);//flag not needed?
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r);
	kernels::interpolateVelocityToHybridNodeY<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r);
	//testInterpX();
	//testInterpY();
	logger.stopTimer("RHS1 Interpolation");
}

void oscCylinder::setVelocityInside()
{
	double	*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) );
	
	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) );
	
	int nx = domInfo ->nx,
		ny = domInfo ->ny;
			
	const int blocksize = 256;
	dim3 grid( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::setInsideVelocity<<<grid,block>>>(ghostTagsUV_r, u_r, uB_r, vB_r, nx, ny);
}