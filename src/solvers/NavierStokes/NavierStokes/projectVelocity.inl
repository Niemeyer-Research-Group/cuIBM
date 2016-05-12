/***************************************************************************//**
 * \file projectVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that update the velocity
 */

#include <solvers/NavierStokes/kernels/projectVelocity.h>

void NavierStokesSolver::velocityProjection()
{
	logger.startTimer("Velocity Projection");
	double 	*u_r	= thrust::raw_pointer_cast( &(u[0]) ),
			*u_old	= thrust::raw_pointer_cast( &(uold[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(uhat[0]) ),
			*pressure_r	= thrust::raw_pointer_cast( &(pressure[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo-> dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo-> dy[0]) );
	
	int *tagsP_r = thrust::raw_pointer_cast( &(tagsP[0]) ),
		*tagsin_r = thrust::raw_pointer_cast( &(tagsIn[0]) );

	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	const int blocksize = 256;
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	if (B.numBodies > 0)
	{
		kernels::project_velocity_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, u_old, pressure_r, tagsP_r, tagsin_r, dx_r, dt, nx, ny);
		kernels::project_velocity_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, u_old, pressure_r, tagsP_r, tagsin_r, dy_r, dt, nx, ny);
	}
	else
	{
		kernels::project_velocity_X_nobody<<<dimGridU,dimBlockU>>>(u_r, uhat_r, u_old, pressure_r, dx_r, dt, nx, ny);
		kernels::project_velocity_Y_nobody<<<dimGridV,dimBlockV>>>(u_r, uhat_r, u_old, pressure_r, dy_r, dt, nx, ny);
	}
	logger.stopTimer("Velocity Projection");
}
