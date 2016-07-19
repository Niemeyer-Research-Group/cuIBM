/***************************************************************************//**
 * \file projectVelocity
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief Functions that call kernels to solve for the final velocity
 */
#include <solvers/NavierStokes/NavierStokesSolver.h>

#include <solvers/NavierStokes/NavierStokes/kernels/projectVelocity.h>

void NavierStokesSolver::velocityProjection()
{
	logger.startTimer("Velocity Projection");

	const int blocksize = 256;

	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	kernels::project_velocity_X_nobody<<<dimGridU,dimBlockU>>>(u_r, uhat_r, uold_r, pressure_r, dx_r, dt, nx, ny);
	kernels::project_velocity_Y_nobody<<<dimGridV,dimBlockV>>>(u_r, uhat_r, uold_r, pressure_r, dy_r, dt, nx, ny);

	logger.stopTimer("Velocity Projection");
}
