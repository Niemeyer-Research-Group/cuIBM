/***************************************************************************//**
 * \file projectVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that update the velocity
 */
#include <solvers/NavierStokes/luoIBM/kernels/projectVelocity.h>
#include <solvers/NavierStokes/luoIBM/kernels/biLinearInterpolation.h>
#include <solvers/NavierStokes/oscCylinder/kernels/intermediateVelocity.h>//set inside
void luoIBM::velocityProjection()
{
	logger.startTimer("Velocity Projection");
	double 	*u_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::u[0]) ),
			*u_old	= thrust::raw_pointer_cast( &(NavierStokesSolver::uold[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::uhat[0]) ),
			*pressure_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::pressure[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo-> dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo-> dy[0]) );
	
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double	dt = db["simulation"]["dt"].get<double>();

	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;

	const int blocksize = 256;
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	//project velocity
	kernels::project_velocity_luo_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, u_old, pressure_r, dx_r, dt, nx, ny);
	kernels::project_velocity_luo_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, u_old, pressure_r, dy_r, dt, nx, ny);
	
	//force inside velocity to be body velocity...
	/*double	*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
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
			*image_point_u_r = thrust::raw_pointer_cast( &(image_point_u[0]) );
	
	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) );
		
	double	nu = db["flow"]["nu"].get<double>();

	int width_i = B.numCellsXHost,
		height_j=B.numCellsYHost;
	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) );
	
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r, image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start_r, j_start_r, width_i, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r, image_point_u_r);
	
	dim3 grid2( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	kernels::setInsideVelocity<<<grid2,block>>>(ghostTagsUV_r, u_r, uB_r, vB_r, nx, ny);
	//testInterpX();
*/
	logger.stopTimer("Velocity Projection");
}
