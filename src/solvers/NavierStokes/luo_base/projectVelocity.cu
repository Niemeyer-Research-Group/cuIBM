/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that update the velocity
 */
#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/luo_base/kernels/projectVelocity.h>
#include <solvers/NavierStokes/luo_base/kernels/biLinearInterpolation.h>
#include <solvers/NavierStokes/luo_base/kernels/intermediateVelocity.h>//set inside
void luo_base::_project_velocity()
{
	logger.startTimer("Velocity Projection");

	const int blocksize = 256;
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	//project velocity
	kernels::project_velocity_luo_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, uold_r, pressure_r, dx_r, dt, nx, ny);
	kernels::project_velocity_luo_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, uold_r, pressure_r, dy_r, dt, nx, ny);

	//force inside velocity to be body velocity...

	dim3 block(blocksize, 1);
	dim3 grid2( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);

	kernels::setInsideVelocity<<<grid2,block>>>(ghostTagsUV_r, u_r, B.uB_r, B.vB_r, nx, ny);
	//testInterpX();

	logger.stopTimer("Velocity Projection");
}
