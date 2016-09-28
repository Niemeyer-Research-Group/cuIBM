/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */
#include <solvers/NavierStokes/NavierStokesSolver.h>

#include <solvers/NavierStokes/NavierStokes/kernels/CFL.h>

void NavierStokesSolver::CFL()
{
	logger.startTimer("CFL");

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::calculateCFL<<<grid,block>>>(cfl_r, u_r, dx_r, dy_r, nx, ny, dt);

	thrust::device_vector<double>::iterator iter = thrust::max_element(cfl.begin(),cfl.end());
	unsigned int position = iter - cfl.begin();
	double max_val = *iter;
	if (max_val > cfl_max)
	{
		if (timeStep>100 && max_val > cfl_max*1.2)
			crash();
		cfl_max = max_val;
		cfl_I = position%nx;
		cfl_J = int(position/nx);
		cfl_ts = timeStep;
	}
	logger.stopTimer("CFL");
}
