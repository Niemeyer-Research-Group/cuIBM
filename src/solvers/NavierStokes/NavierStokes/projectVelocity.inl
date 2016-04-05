/***************************************************************************//**
 * \file projectVelocity.inl
 * \author Chris Minar
 * \brief Implementation of the methods to calulate step 3: solve for u_l+1
 * \		
 */

#include <solvers/NavierStokes/kernels/projectVelocity.h>

void NavierStokesSolver::velocityProjection()
{
	logger.startTimer("Velocity Projection");
	double 	*u_r	= thrust::raw_pointer_cast( &(u[0]) ),
			*u_old	= thrust::raw_pointer_cast( &(uold[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(uhat[0]) ),
			*pressure_r	= thrust::raw_pointer_cast( &(pressure[0]) ),
			*dxD_r	= thrust::raw_pointer_cast( &(domInfo-> dxD[0]) ),
			*dyD_r	= thrust::raw_pointer_cast( &(domInfo-> dyD[0]) );
	
	int *tagsP_r = thrust::raw_pointer_cast( &(tagsPD[0]) ),
		*tagsin_r = thrust::raw_pointer_cast( &(tagsInD[0]) );

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
		kernels::project_velocity_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, u_old, pressure_r, tagsP_r, tagsin_r, dxD_r, dt, nx, ny);
		kernels::project_velocity_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, u_old, pressure_r, tagsP_r, tagsin_r, dyD_r, dt, nx, ny);
	}
	else
	{
		kernels::project_velocity_X_nobody<<<dimGridU,dimBlockU>>>(u_r, uhat_r, u_old, pressure_r, dxD_r, dt, nx, ny);
		kernels::project_velocity_Y_nobody<<<dimGridV,dimBlockV>>>(u_r, uhat_r, u_old, pressure_r, dyD_r, dt, nx, ny);
	}
	logger.stopTimer("Velocity Projection");
}
