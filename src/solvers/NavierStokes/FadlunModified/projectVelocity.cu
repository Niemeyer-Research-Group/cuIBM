/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that update the velocity
 */
#include <solvers/NavierStokes/fadlunModified.h>

#include <solvers/NavierStokes/FadlunModified/kernels/projectVelocity.h>

void fadlunModified::velocityProjection()
{
	logger.startTimer("Velocity Projection");

	const int blocksize = 256;

	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	kernels::project_velocity_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, uold_r, pressure_r, tagsP_r, tagsIn_r, dx_r, dt, nx, ny);
	kernels::project_velocity_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, uold_r, pressure_r, tagsP_r, tagsIn_r, dy_r, dt, nx, ny);

	logger.stopTimer("Velocity Projection");
}
